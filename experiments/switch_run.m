function row = switch_run(N, method, opts)
%SWITCH_RUN  One (size, method) of Example 3, switching diffusion (SWITCH_PROBLEM; Fig. 9).
%   Meant to run in its own MATLAB process (SWITCH_DRIVE) so that the peak memory is the operating system's
%   number (read from the .NET process object, Windows).  The outputs kept by every method are the solution on
%   the slice of cells through the target and its mean over the cube, for every state.
%
%   ROW = SWITCH_RUN(N, METHOD, OPTS)
%     METHOD  'gfadi'    one-sided g-fADI (GFADI_MEGA), Chebyshev zeros on E, a priori number of steps
%             'lanczos'  one-pass block Lanczos Galerkin, same outputs (TB_MEGA_LANCZOS)
%             'arnoldi'  one-pass block Arnoldi (TB_ADV_ARNOLDI), for the swirl
%             'eksm'     extended Krylov with one sparse LU of A
%             'gs'       block Gauss-Seidel over the states, warm-started CG (BiCGStab with swirl) per state;
%                        projected beyond opts.tcap
%             'bs'       Bartels-Stewart, complex Schur form of L: one shifted solve (sparse LU) timed, times n
%             'check'    (small N) g-fADI outputs against a converged full solution, and its residual
%     OPTS    .tol (1e-8) .tcap (1800 s) .threads (8) .swirl (0) .Ebottom ('coarse', or 'cw': rigorous E(1))
%   Appends ROW to results/switch<tag><cwtag>.csv (tag = _swirl<W> if swirl > 0, cwtag = _cw for 'cw').  'gfadi'
%   saves the outputs to results/switch_ref/N<N><tag>.mat (N<N><tag>_cw.mat for 'cw'), which 'lanczos', 'arnoldi',
%   'eksm', 'gs' and 'check' compare against; 'check' also saves results/switch_ref/N<N><tag>_full.mat.
if nargin < 3, opts = struct(); end
o = struct('tol', 1e-8, 'tcap', 1800, 'threads', 8, 'swirl', 0, 'Ebottom', 'coarse');
fn = fieldnames(opts);
for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
setup_paths();
here = app_root();
maxNumCompThreads(o.threads);
if isempty(gcp('nocreate')), parpool('Threads', o.threads); end
tb = tic;
P = switch_problem(N, struct('tol', o.tol, 'swirl', o.swirl, 'Ebottom', o.Ebottom, 'sparse', any(strcmp(method, {'bs', 'check', 'eksm'}))));
tbuild = toc(tb);
tag = '';
if o.swirl > 0, tag = sprintf('_swirl%g', o.swirl); end
cwtag = '';
if strcmp(o.Ebottom, 'cw'), cwtag = '_cw'; end                  % rigorous E(1): separate csv and outputs
if strcmp(o.Ebottom, 'cw')
    fprintf('switch N = %d: lambda_min(H) >= %.4g (torsion %.4g, one inverse step %.4g), %d products, %.1f s\n', ...
            N, P.E(1), P.Einfo.bound1, P.Einfo.bound2, P.Einfo.products, P.tE);
end
m = P.m;  n = P.n;  ks = round(P.o.c(3)*N + 0.5);                          % the slice z through the centre of the target
outf = @(Z) [reshape(Z(:, :, ks, :), N*N, []); sum(reshape(Z, m, []), 1)/m];
fprintf('switch N = %d %s: built in %.1f s, E = [%.3g %.3g], f = %.3g (Re W(L) >= %.3g), p = %d\n', N, method, tbuild, P.E, P.f, P.WL(1), P.p);
row = struct('N', N, 'method', method, 'm', m, 'n', n, 'unknowns', m*n, 'p', P.p, 'f', P.f, 'time', NaN, 'products', NaN, ...
             'k', NaN, 's', NaN, 'err', NaN, 'res', NaN, 'peakGB', NaN, 'fillGB', NaN, 'projected', 0, 'note', '');
reff = fullfile(here, 'results', 'switch_ref', sprintf('N%d%s.mat', N, tag));
try
switch method
    case 'gfadi'
        [Zo, Qy, I] = gfadi_mega(P.op, P.U, P.Bop, P.V, P.al, outf, 'verbose', 600);
        row.time = I.time + P.tE + P.tshift;  row.products = P.p;  row.s = I.s;      % the bound on E and the shift design counted
        Xo = (double(Zo)*Qy').*P.back';                                     % back to u: columns scaled
        row.note = sprintf('products %.0f s, small %.0f s, outputs %.0f s, shift design %.1f s', I.tproduct, I.tsmall, I.toutput, P.tshift);
        if strcmp(o.Ebottom, 'cw')
            row.products = P.p + P.Einfo.products;
            row.note = sprintf('%s; E(1) = %.4g by Collatz-Wielandt, %d CG products, %.1f s', row.note, P.E(1), P.Einfo.products, P.tE);
            if exist(reff, 'file'), R = load(reff);  row.err = norm(Xo - R.Xo, 'fro')/norm(R.Xo, 'fro'); end   % against the earlier run
            reff = strrep(reff, '.mat', '_cw.mat');
        end
        if ~exist(fileparts(reff), 'dir'), mkdir(fileparts(reff)); end
        save(reff, 'Xo', 'ks');
    case 'lanczos'
        R = load(reff);
        [Yo, I] = tb_mega_lanczos(P.op, P.U, P.L, P.V, outf, o.tol, struct('tcap', o.tcap));
        Xo = real(Yo).*P.back';
        row.time = I.time;  row.products = I.products;  row.k = I.k;  row.res = I.res;  row.projected = I.capped;
        row.err = norm(Xo - R.Xo, 'fro')/norm(R.Xo, 'fro');
        row.note = sprintf('k = %d: products %.0f s, vectors %.0f s, projected %.0f s, outputs %.0f s', I.k, I.tproduct, I.tvec, I.tsmall, I.tout);
    case 'arnoldi'
        R = load(reff);
        [Yo, I] = tb_adv_arnoldi(P.op, P.U, P.L, P.V, outf, o.tol, struct('tcap', o.tcap, 'kcap', ceil(0.5*P.p)));
        Xo = real(Yo).*P.back';
        row.time = I.time;  row.products = I.products;  row.k = I.k;  row.res = I.res;  row.projected = I.capped;
        row.err = norm(Xo - R.Xo, 'fro')/norm(R.Xo, 'fro');
        row.note = sprintf('k = %d: products %.0f s, orthogonalisation %.0f s, basis %.2f GB', I.k, I.tproduct, I.torth, I.basisGB);
    case 'eksm'
        R = load(reff);
        [Yo, I] = eksm_outputs(P, outf, o);
        Xo = Yo.*P.back';
        row.time = I.time;  row.products = I.products;  row.k = I.k;  row.res = I.res;  row.projected = I.capped;
        row.fillGB = I.fillGB;  row.err = norm(Xo - R.Xo, 'fro')/norm(R.Xo, 'fro');
        row.note = sprintf('k = %d: one real LU %.1f s (%.2f GB), products %.0f s, solves %.0f s, orthogonalisation %.0f s, projected %.0f s, basis %.2f GB', ...
                           I.k, I.tlu, I.fillGB, I.tprod, I.tsolve, I.torth, I.tproj, I.basisGB);
    case 'gs'
        [Xo, I] = block_gs(P, outf, o);
        row.time = I.time;  row.products = I.products;  row.k = I.sweeps;  row.projected = I.projected;
        if exist(reff, 'file'), R = load(reff);  row.err = norm(Xo - R.Xo, 'fro')/norm(R.Xo, 'fro'); end
        row.note = sprintf('%d sweeps done in %.0f s (of %d needed, contraction %.4f per sweep), %d products so far', I.done, I.tdone, I.sweeps, I.rate, I.pdone);
    case 'bs'
        lam = eig(P.L);  lc = lam(abs(imag(lam)) > 0);  [~, j] = min(abs(real(lc) - median(real(lc))));  a = lc(j);   % a typical complex shift
        t = tic;  [Lf, Uf, Pp, Qp] = lu(P.As + a*speye(m));  tfac = toc(t);  %#ok<ASGLU>
        row.fillGB = (nnz(Lf) + nnz(Uf))*16/2^30;  clear Lf Uf Pp Qp
        row.time = n*tfac;  row.projected = 1;
        row.note = sprintf('one complex LU %.2f s (shift %.3g%+.3gi), times n = %d', tfac, real(a), imag(a), n);
    case 'check'
        R = load(reff);
        Ud = full_solution(P);
        Xd = outf(reshape(Ud, [P.sz n]));
        row.err = norm(R.Xo - Xd, 'fro')/norm(Xd, 'fro');
        row.res = norm(P.As*Ud.*P.Dk' + Ud*P.G.' - P.b(:)*P.Dk', 'fro')/norm(P.b(:)*P.Dk', 'fro');
        save(fullfile(here, 'results', 'switch_ref', sprintf('N%d%s_full.mat', N, tag)), 'Ud');
        row.note = 'g-fADI outputs against a converged block Gauss-Seidel solution with sparse direct inner solves';
end
catch err
    row.note = ['ERROR: ' err.message];  fprintf(2, '%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
end
pk = System.Diagnostics.Process.GetCurrentProcess();             % peak memory of this process (Windows .NET)
row.peakGB = double(pk.PeakWorkingSet64)/2^30;
append_row(fullfile(here, 'results', ['switch' tag cwtag '.csv']), row);
fprintf('switch N = %d %s: %.1f s, products %g, err %.1e, res %.1e, peak %.1f GB, %s\n', N, method, row.time, row.products, row.err, row.res, row.peakGB, row.note);
end

function [Xo, I] = block_gs(P, outf, o)
% block Gauss-Seidel over the states:  (D_k A + G_kk) u_k = D_k b - sum_{j ~= k} G_kj u_j,  CG for each, warm started
m = P.m;  n = P.n;  sz = P.sz;  G = P.G;  Dk = P.Dk;  b = P.b(:);  U = zeros(m, n);  prods = 0;
t0 = tic;  done = 0;  ch = [];  lastt = 0;  nb = norm(b);
for sweep = 1:100000
    ts = tic;  dmax = 0;  umax = 0;
    for k = 1:n
        rhs = Dk(k)*b - U*G(k, :).' + G(k, k)*U(:, k);
        if o.swirl > 0
            [x, it] = bicg_shift(P.op, sz, Dk(k), G(k, k), rhs, U(:, k), 1e-2*o.tol);
        else
            [x, it] = cg_shift(P.op, sz, Dk(k), G(k, k), rhs, U(:, k), 1e-2*o.tol);
        end
        prods = prods + it;
        dmax = max(dmax, norm(x - U(:, k)));  umax = max(umax, norm(x));  U(:, k) = x;
    end
    done = sweep;  lastt = toc(ts);  ch(end + 1) = dmax/max(umax, realmin);   %#ok<AGROW>
    fprintf('   sweep %d: %.0f s, change %.2e, %d products so far\n', sweep, lastt, ch(end), prods);
    if ch(end) < o.tol, break, end
    if toc(t0) > o.tcap && sweep >= 3, break, end
end
tdone = toc(t0);  pdone = prods;
rate = NaN;  if numel(ch) >= 3, rate = ch(end)/ch(end - 1); end
if ch(end) < o.tol
    I = struct('time', tdone, 'products', pdone, 'projected', 0, 'sweeps', done, 'done', done, 'tdone', tdone, 'pdone', pdone, 'rate', rate);
else
    need = done + ceil(log(o.tol/ch(end))/log(rate));
    I = struct('time', tdone + (need - done)*lastt, 'products', pdone + (need - done)*pdone/done, 'projected', 1, 'sweeps', need, ...
               'done', done, 'tdone', tdone, 'pdone', pdone, 'rate', rate);
end
Xo = outf(reshape(U, [sz n]));
end

function [Yo, I] = eksm_outputs(P, outf, o)
% Extended Krylov (Simoncini 2007) for A Y + Y L.' = U V' with the sparse A: one sparse LU of A, the space
% span{U, A^{-1}U, A U, A^{-2}U, ...} with full orthogonalization, Galerkin with the small side kept exact.
A = P.As;  At = A';  m = P.m;  b = double(P.U(:));  c = double(P.V);  L = P.L;  t0 = tic;
t = tic;  [Lf, Uf, Pp, Qp] = lu(A);  tlu = toc(t);  fillGB = (nnz(Lf) + nnz(Uf))*8/2^30;
solveA = @(y) Qp*(Uf\(Lf\(Pp*y)));
kcap = 600;  chunk = 64;  V = zeros(m, chunk);  TT = zeros(2*kcap + 2);  tprod = 0;  tsolve = 0;  torth = 0;  tproj = 0;  nprod = 0;
t = tic;  w = solveA(b);  tsolve = tsolve + toc(t);
[Q0, ~] = qr([b, w], 0);  V(:, 1:2) = Q0;  ncol = 2;  nrmF = norm(b)*norm(c);  res = Inf;  capped = 1;  Yk = [];
for k = 1:kcap
    j = ncol - 1:ncol;                                                       % the newest block
    t = tic;  AVj = A*V(:, j);  tprod = tprod + toc(t);  nprod = nprod + 2;
    t = tic;  wn = solveA(V(:, j(2)));  tsolve = tsolve + toc(t);
    t = tic;  W = [AVj(:, 1), wn];
    for pass = 1:2, W = W - V(:, 1:ncol)*(V(:, 1:ncol)'*W); end            % full, twice
    [Qn, ~] = qr(W, 0);
    if ncol + 2 > size(V, 2), V(:, end + chunk) = 0; end                   % grow the basis in chunks
    V(:, ncol + 1:ncol + 2) = Qn;  torth = torth + toc(t);
    t = tic;  AtQ = At*Qn;  tprod = tprod + toc(t);  nprod = nprod + 2;
    t = tic;
    TT(1:ncol + 2, j) = V(:, 1:ncol + 2)'*AVj;                              % new columns of V' A V
    if ncol > 2, TT(ncol + 1:ncol + 2, 1:ncol - 2) = AtQ'*V(:, 1:ncol - 2); end   % new rows: Qn' A V_old
    T = TT(1:ncol, 1:ncol);  rhs = (V(:, 1:ncol)'*b)*c';
    Yk = sylvester(T, L.', rhs);
    res = norm(TT(ncol + 1:ncol + 2, 1:ncol)*Yk, 'fro')/nrmF;  tproj = tproj + toc(t);
    ncol = ncol + 2;
    if mod(k, 10) == 0, fprintf('   eksm k = %d (%d vectors): residual %.2e, %.0f s\n', k, ncol, res, toc(t0)); end
    if res <= o.tol, capped = 0;  break, end
    if toc(t0) > o.tcap, break, end
end
nk = size(Yk, 1);
Yo = outf(reshape(V(:, 1:nk), [P.sz nk]))*Yk;
I = struct('time', toc(t0), 'products', nprod, 'k', nk, 'res', res, 'capped', capped, 'tlu', tlu, 'fillGB', fillGB, ...
           'tprod', tprod, 'tsolve', tsolve, 'torth', torth, 'tproj', tproj, 'basisGB', nk*m*8/2^30);
end
function [x, it] = cg_shift(op, sz, D, g, b, x, tol)
% CG for (D A + g I) x = b, matrix-free:  op.mv(x, a, s) = s (A - a I) x
mv = @(v) reshape(op.mv(reshape(v, sz), -g/D, D), [], 1);
r = b - mv(x);  pv = r;  rr = r'*r;  nb = max(norm(b), realmin);  it = 1;
while sqrt(rr) > tol*nb && it < 20000
    Ap = mv(pv);  it = it + 1;  al = rr/(pv'*Ap);  x = x + al*pv;  r = r - al*Ap;  rn = r'*r;  pv = r + (rn/rr)*pv;  rr = rn;
end
end

function [x, it] = bicg_shift(op, sz, D, g, b, x, tol)
% BiCGStab for (D A + g I) x = b with a non-symmetric A, matrix-free (two products per iteration)
mv = @(v) reshape(op.mv(reshape(v, sz), -g/D, D), [], 1);
r = b - mv(x);  rt = r;  rho0 = 1;  al = 1;  om = 1;  v = zeros(size(b));  pv = v;  nb = max(norm(b), realmin);  it = 1;
while norm(r) > tol*nb && it < 40000
    rho1 = rt'*r;  be = (rho1/rho0)*(al/om);  pv = r + be*(pv - om*v);
    v = mv(pv);  al = rho1/(rt'*v);  s = r - al*v;
    t = mv(s);  it = it + 2;  om = (t'*s)/(t'*t);
    x = x + al*pv + om*s;  r = s - om*t;  rho0 = rho1;
end
end

function U = full_solution(P)
% the coupled system by block Gauss-Seidel with sparse direct inner solves, to 1e-13 (small N only)
m = P.m;  n = P.n;  G = P.G;  Dk = P.Dk;  b = P.b(:);  U = zeros(m, n);  F = cell(n, 1);
for k = 1:n, F{k} = decomposition(Dk(k)*P.As + G(k, k)*speye(m), 'lu'); end
for sweep = 1:5000
    dmax = 0;
    for k = 1:n
        x = F{k}\(Dk(k)*b - U*G(k, :).' + G(k, k)*U(:, k));  dmax = max(dmax, norm(x - U(:, k))/max(norm(x), realmin));  U(:, k) = x;
    end
    if dmax < 1e-13, break, end
end
end

function append_row(f, row)
newf = ~exist(f, 'file');  fid = fopen(f, 'a');  fn = fieldnames(row);
if newf, fprintf(fid, '%s\n', strjoin(fn', ',')); end
v = cell(1, numel(fn));
for i = 1:numel(fn)
    x = row.(fn{i});
    if ischar(x), v{i} = ['"' strrep(x, '"', '''') '"']; else, v{i} = sprintf('%.6g', x); end
end
fprintf(fid, '%s\n', strjoin(v, ','));  fclose(fid);
end
