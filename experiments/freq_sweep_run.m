function row = freq_sweep_run(N, method, opts)
%FREQ_SWEEP_RUN  Example 4: a frequency sweep (A + i w_j I) x_j = b, j = 1..nw, by one method (Fig. 10(d,e)).
%   A = -div(a grad) + w.grad on N^3 interior points, Dirichlet, centered differences for the drift (A is not
%   symmetric).  a = 1..10 smooth, or (opts.C > 0) 27 spherical inclusions of a = C in a = 1.  The frequencies
%   are log-spaced over [1e-2 lambda_min(H), lambda_max(H)].  Meant for one (N, method) per MATLAB process, so
%   that the peak memory is the process's own.
%
%   ROW = FREQ_SWEEP_RUN(N, METHOD, OPTS)
%     METHOD  all solutions:  'gfadi' (one-sided g-fADI, GFADI_LOWRANK), 'msfom' (multishift FOM, stored basis),
%                             'rsfom' (restarted multishift FOM(mr)), 'gmres' (GMRES(100) per frequency),
%                             'lu' (sparse LU per frequency)
%             outputs only (the response at ell sensors): 'gfadi_out', 'msfom_out', 'rsfom_out', 'bicg_out'
%                             (multishift BiCG), 'gmres_out'
%     OPTS    .nw (1000) .W (drift, 1) .tol (1e-8) .mr (50) .tcap (7200 s) .threads (8) .C (0) .ell (100)
%   Appends ROW to results/freq_sweep.csv, or results/freq_sweep_C<C>.csv when C > 0.  Accuracy: relative residual
%   <= tol for every frequency (all-solutions methods), or the sensor error at 12 frequencies (outputs only).
%   Peak memory is read from the .NET process object (Windows).
%
%   Fig. 10 reads results/freq_sweep_C1000.csv, written for each N by
%     freq_sweep_run(N, 'gfadi_out', struct('C', 1000))
%     freq_sweep_run(N, 'msfom_out', struct('C', 1000))
%     freq_sweep_run(N, 'rsfom_out', struct('C', 1000, 'mr', 200))
if nargin < 3, opts = struct(); end
o = struct('nw', 1000, 'W', 1, 'tol', 1e-8, 'mr', 50, 'tcap', 7200, 'threads', 8, 'C', 0, 'ell', 100);   % C > 0: inclusions of a = C in a = 1
fn = fieldnames(opts);
for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
setup_paths();
here = app_root();
maxNumCompThreads(o.threads);
ts = tic;
[A, H, b, E, c, ymax] = build_op(N, o.W, o.C);
m = N^3;
tsetup = toc(ts);
om = logspace(log10(1e-2*E(1)), log10(E(2)), o.nw).';  nw = o.nw;
row = struct('N', N, 'method', method, 'm', m, 'nw', nw, 'setup', tsetup, 'time', NaN, 'products', NaN, 'k', NaN, 'maxres', NaN, ...
             'peakGB', NaN, 'storeGB', NaN, 'note', '');
fprintf('sweep N = %d %s: E = [%.3g %.3g], c = %.3g, nw = %d\n', N, method, E, c, nw);
try
t0 = tic;
switch method
    case 'gfadi'
        [p, al] = adv_shifts(E, c, ymax, [0 0], o.tol/4);  al = leja_order(al, E');   % margin: the bound is on ||residual|| of the whole block
        B = -1i*spdiags(om, 0, nw, nw);                                  % A X - X B = b 1^T  <=>  (A + i w_j) x_j = b
        ytol = min(1e-13, 0.1*o.tol*E(1)/E(2));                          % X to ytol*||X|| keeps the residual below tol
        [Z, Y] = gfadi_lowrank(gfadi_op(A), gfadi_op(B), b, ones(nw, 1), al, [], 'ytol', ytol);
        row.time = toc(t0);  row.peakGB = peak();  row.products = p;  row.k = size(Z, 2);  row.storeGB = (numel(Z)*8 + numel(Y)*16)/2^30;
        row.maxres = check_res(A, b, om, @(J) Z*Y(J, :)');
        row.note = sprintf('p = %d, rank %d', p, size(Z, 2));
    case 'gmres_out'
        % each frequency on its own, GMRES(100), keeping only the sensor values; stops at the time cap and records how
        % many frequencies were done (a measured outcome, no extrapolation)
        rng(3);  Ci = sort(randperm(m, o.ell));  np = 0;  done = 0;  mx = 0;  Im = speye(m);
        for j = 1:nw
            Aj = A + 1i*om(j)*Im;  [x, ~, ~, it] = gmres(Aj, b, 100, o.tol, ceil(20000/100));
            np = np + (it(1) - 1)*100 + it(2);  done = j;  mx = max(mx, norm(Aj*x - b)/norm(b));
            if toc(t0) > o.tcap, break, end
        end
        row.time = toc(t0);  row.peakGB = peak();  row.products = np;  row.maxres = mx;  row.storeGB = 103*m*16/2^30;
        row.note = sprintf('GMRES(100) per frequency: %d of %d frequencies in %.0f s, %d products', done, nw, row.time, np);
        if done < nw, row.note = ['STOPPED: ' row.note]; end
    case {'gfadi_out', 'msfom_out', 'rsfom_out', 'bicg_out'}
        % only outputs wanted: the response at ell sensors (grid points) at every frequency, Y = C X (ell x nw)
        rng(3);  Ci = sort(randperm(m, o.ell));
        if strcmp(method, 'gfadi_out')
            [p, al] = adv_shifts(E, c, ymax, [0 0], o.tol/4);  al = leja_order(al, E');
            Yout = gfadi_outputs(A, b, -1i*om, al, Ci);           % right steps only, plain loop (see GFADI_SIMPLE)
            row.time = toc(t0);  row.peakGB = peak();  row.products = p;  row.storeGB = 3*m*16/2^30;
            row.note = sprintf('p = %d, %d sensors', p, o.ell);
        elseif strcmp(method, 'bicg_out')
            [Yout, np, k, I] = bicg_outputs(A, b, 1i*om, o.tol, Ci, o.tcap);   % short recurrences, no basis
            row.time = toc(t0);  row.peakGB = peak();  row.products = np;  row.k = k;  row.storeGB = 6*m*8/2^30;
            row.note = sprintf('multi-shift BiCG: %d steps, %d products (A and A^T), %d of %d frequencies converged, %d sensors', ...
                               k, np, I.nconv, nw, o.ell);
            if I.nconv < nw, row.note = ['STOPPED: ' row.note]; end
        elseif strcmp(method, 'rsfom_out')
            [Yout, np] = rsfom(A, b, om, o.tol, o.mr, o.tcap, Ci);   % restarted: keeps only the sensor values
            row.time = toc(t0);  row.peakGB = peak();  row.products = np;  row.storeGB = m*(o.mr + 1)*8/2^30;
            row.note = sprintf('FOM(%d), %d products, %d sensors', o.mr, np, o.ell);
        else
            [V, Yk, k] = msfom(A, b, om, o.tol);  Yout = V(Ci, :)*Yk;
            row.time = toc(t0);  row.peakGB = peak();  row.products = k;  row.k = k;  row.storeGB = numel(V)*8/2^30;
            row.note = sprintf('k = %d, basis %.2f GB, %d sensors', k, row.storeGB, o.ell);
        end
        % accuracy: exact outputs at 12 frequencies (sparse direct), relative error of the sensor responses
        J = unique(round(linspace(1, nw, 12)));  mx = 0;
        for j = J, x = (A + 1i*om(j)*speye(m))\b;  mx = max(mx, norm(Yout(:, j) - x(Ci))/norm(x(Ci))); end
        row.maxres = mx;  row.note = [row.note sprintf(', output error %.1e at %d frequencies', mx, numel(J))];
    case 'msfom'
        [V, Yk, k] = msfom(A, b, om, o.tol);
        row.time = toc(t0);  row.peakGB = peak();  row.products = k;  row.k = k;  row.storeGB = numel(V)*8/2^30;
        row.maxres = check_res(A, b, om, @(J) V*Yk(:, J));
        row.note = sprintf('k = %d, basis %.2f GB', k, row.storeGB);
    case 'rsfom'
        [X, np] = rsfom(A, b, om, o.tol, o.mr, o.tcap);
        row.time = toc(t0);  row.peakGB = peak();  row.products = np;  row.storeGB = (numel(X)*16 + m*(o.mr + 1)*8)/2^30;
        row.maxres = check_res(A, b, om, @(J) X(:, J));
        row.note = sprintf('FOM(%d), %d products', o.mr, np);
    case 'gmres'
        np = 0;  mx = 0;
        for j = 1:nw
            Aj = A + 1i*om(j)*speye(m);
            [x, ~, ~, it] = gmres(Aj, b, 100, o.tol, 200);  np = np + (it(1) - 1)*100 + it(2);
            mx = max(mx, norm(Aj*x - b)/norm(b));
            if toc(t0) > o.tcap, error('time cap after %d of %d frequencies', j, nw); end
        end
        row.time = toc(t0);  row.products = np;  row.maxres = mx;  row.storeGB = 103*m*16/2^30;
        row.note = sprintf('%d products', np);
    case 'lu'
        mx = 0;
        for j = 1:nw
            Aj = A + 1i*om(j)*speye(m);  x = Aj\b;  mx = max(mx, norm(Aj*x - b)/norm(b));
            if toc(t0) > o.tcap, error('time cap after %d of %d frequencies', j, nw); end
        end
        row.time = toc(t0);  row.maxres = mx;  row.note = 'sparse LU per frequency (backslash)';
end
catch err
    row.note = ['ERROR: ' err.message];  fprintf(2, '%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
end
if isnan(row.peakGB), row.peakGB = peak(); end      % peak memory of the solve, measured before the residual check
csvf = 'freq_sweep.csv';
if o.C > 0, csvf = sprintf('freq_sweep_C%g.csv', o.C); end
append_row(fullfile(here, 'results', csvf), row);
fprintf('sweep N = %d %s: %.1f s, products %g, max residual %.1e, peak %.1f GB, %s\n', N, method, row.time, row.products, row.maxres, row.peakGB, row.note);
end

function [A, H, b, E, c, ymax] = build_op(N, W, C)
h = 1/(N + 1);  g = (1:N)'*h;  e = ones(N, 1);  I = speye(N);
D = spdiags([-e e], [-1 0], N + 1, N)/h;                              % forward differences incl. the Dirichlet faces
gf = ((0:N)' + 0.5)*h;                                                % face positions
if C > 0     % composite: 27 spheres (radius 0.11, centres on a jittered 3x3x3 lattice) of conductivity C in a = 1
    ctr = 0.25*[kron([1;1;1], kron([1;1;1], (1:3)')), kron([1;1;1], kron((1:3)', [1;1;1])), kron((1:3)', kron([1;1;1], [1;1;1]))];
    ctr = ctr + 0.04*sin((1:27)'*[1.3 2.1 3.7]);
    af = @(x, y, z) 1 + (C - 1)*incl(x, y, z, ctr, 0.11);
else
    af = @(x, y, z) 10.^((1 + sin(2*pi*x + 1).*sin(3*pi*y).*sin(2*z + 0.5))/2);   % 1 .. 10, smooth, not separable
end
[X, Y, Z] = ndgrid(gf, g, g);  Ax = af(X, Y, Z);  [X, Y, Z] = ndgrid(g, gf, g);  Ay = af(X, Y, Z);  [X, Y, Z] = ndgrid(g, g, gf);  Az = af(X, Y, Z);
Dx = kron(I, kron(I, D));  Dy = kron(I, kron(D, I));  Dz = kron(D, kron(I, I));
H = Dx'*spdiags(Ax(:), 0, numel(Ax), numel(Ax))*Dx + Dy'*spdiags(Ay(:), 0, numel(Ay), numel(Ay))*Dy + Dz'*spdiags(Az(:), 0, numel(Az), numel(Az))*Dz;
wd = W*[1 0.6 0.3]/norm([1 0.6 0.3]);  C1 = spdiags([-e 0*e e], -1:1, N, N)/(2*h);
S = wd(1)*kron(I, kron(I, C1)) + wd(2)*kron(I, kron(C1, I)) + wd(3)*kron(C1, kron(I, I));
A = H + S;
[X, Y, Z] = ndgrid(g, g, g);  b = exp(-((X - 0.3).^2 + (Y - 0.4).^2 + (Z - 0.5).^2)/0.01);  b = b(:)/norm(b(:));
% a >= 1, so H >= -Lap_h and lambda_min(H) >= 3 (4/h^2) sin^2(pi h/2): a free, rigorous bound (no eigs, no factorization)
E = [3*(4/h^2)*sin(pi*h/2)^2, full(max(sum(abs(H), 2)))];  c = W;  ymax = sum(abs(wd))/h;
end

function Y = gfadi_outputs(A, b, sB, al, Ci)
% One-sided g-fADI for A X - X diag(sB) = b 1^T, keeping only Y = X(Ci, :): right steps with zeros al.
% Right step (GFADI_SIMPLE): Z_j = -Zt, Y_j = (B' - conj(a) I)^{-1} Yt, Zt <- -(A - a I) Z_j, Yt <- Y_j.
Zt = b;  yt = ones(numel(sB), 1);  Y = zeros(numel(Ci), numel(sB));  dB = conj(sB(:));
nb = 256;  Zb = zeros(numel(Ci), nb);  Yb = complex(zeros(numel(sB), nb));  ib = 0;   % outputs added in blocks (one GEMM)
for j = 1:numel(al)
    a = al(j);
    yj = yt./(dB - conj(a));                    % the solve with the diagonal B
    ib = ib + 1;  Zb(:, ib) = Zt(Ci);  Yb(:, ib) = yj;     % the outputs of Z_j Y_j^*  (Z_j = -Zt)
    if ib == nb || j == numel(al), Y = Y - Zb(:, 1:ib)*Yb(:, 1:ib)';  ib = 0; end
    Zt = A*Zt - a*Zt;                           % Zt <- (A - a I) Zt  (the two sign flips cancel)
    s = norm(Zt);  Zt = Zt/s;  yt = yj*s;       % rescale, carry the scale on the small side
end
end

function g = peak()
% peak memory of this process in GB (Windows .NET)
pk = System.Diagnostics.Process.GetCurrentProcess();
pk.Refresh();
g = double(pk.PeakWorkingSet64)/2^30;
end

function s = incl(x, y, z, ctr, r)
s = false(size(x));
for i = 1:size(ctr, 1), s = s | ((x - ctr(i, 1)).^2 + (y - ctr(i, 2)).^2 + (z - ctr(i, 3)).^2 < r^2); end
s = double(s);
end

function mx = check_res(A, b, om, Xfun)
m = numel(b);  mx = 0;  nb = norm(b);
for J = 1:10:numel(om)
    JJ = J:min(J + 9, numel(om));  X = Xfun(JJ);
    R = A*X + X.*(1i*om(JJ).') - b;  mx = max(mx, max(sqrt(sum(abs(R).^2, 1)))/nb);
end
end

function [V, Yk, k] = msfom(A, b, om, tol)
% FOM for all shifts in one Arnoldi basis: (H_k + i w I) y = ||b|| e1; residual of shift j = h_{k+1,k} |e_k' y_j|
m = numel(b);  kmax = 30000;  V = zeros(m, 256);  Hh = zeros(257, 256);  nb = norm(b);  V(:, 1) = b/nb;
tp = 0;  to = 0;  tc = 0;  gc = zeros(kmax, 1);  gs = zeros(kmax, 1);  gres = 1;   % gres: GMRES residual / ||b||
for k = 1:kmax
    t = tic;  w = A*V(:, k);  tp = tp + toc(t);  t = tic;
    for pass = 1:2, h = V'*w;  h(k + 1:end) = 0;  w = w - V*h;  Hh(1:k, k) = Hh(1:k, k) + h(1:k); end   % whole V: no copy (its unused columns are zero)
    Hh(k + 1, k) = norm(w);
    if k + 1 > size(V, 2) && (size(V, 2) + 256)*m*8 > 24*2^30, error('basis of %d vectors would exceed 24 GB (32 GB machine)', k); end
    if k + 1 > size(V, 2), V = [V, zeros(m, 256)];  Hh(end + 256, end + 256) = 0; end   %#ok<AGROW> chunks of 256
    V(:, k + 1) = w/Hh(k + 1, k);
    to = to + toc(t);  t = tic;
    % O(k) Givens update (as in GMRES) for the hardest, smallest frequency: FOM residual = GMRES residual / |c_k|
    hc = Hh(1:k + 1, k);  hc(k) = hc(k) + 1i*om(1);
    for i = 1:k - 1, tg = gc(i)*hc(i) + gs(i)*hc(i + 1);  hc(i + 1) = -conj(gs(i))*hc(i) + gc(i)*hc(i + 1);  hc(i) = tg; end
    rr = hypot(abs(hc(k)), abs(hc(k + 1)));  gc(k) = abs(hc(k))/rr;  gs(k) = (hc(k)/max(abs(hc(k)), realmin))*conj(hc(k + 1))/rr;
    if abs(hc(k)) == 0, gc(k) = 0;  gs(k) = 1; end
    gres = gres*abs(-conj(gs(k)));  fomres = gres/max(gc(k), realmin);
    if fomres <= tol && mod(k, 10) == 0          % then all frequencies, exactly
        [Q, T] = schur(Hh(1:k, 1:k), 'complex');  e1 = Q'*eye(k, 1);
        Zs = zeros(k, numel(om));  sg = 1i*om(:).';  dT = diag(T);      % back substitution for all shifts at once
        for i = k:-1:1, Zs(i, :) = (nb*e1(i) - T(i, i + 1:k)*Zs(i + 1:k, :))./(dT(i) + sg); end
        Yk = Q*Zs;  rmax = max(Hh(k + 1, k)*abs(Yk(k, :))/nb);
        if rmax <= tol, tc = tc + toc(t);  break, end
    end
    tc = tc + toc(t);
end
V = V(:, 1:k);
fprintf('msfom: products %.1f s, orthogonalisation %.1f s, small problems %.1f s\n', tp, to, tc);
end

function [X, np] = rsfom(A, b, om, tol, mr, tcap, Ci)
% restarted multishift FOM: every residual is a multiple of the same vector, r_j = beta_j v
% (with CI: only the rows CI of the solutions are kept, X is numel(CI) x nw)
m = numel(b);  nw = numel(om);
if nargin < 7, Ci = 1:m; end
X = complex(zeros(numel(Ci), nw));  v = b;  bet = norm(b)*ones(nw, 1);  v = v/norm(v);  np = 0;  nb = norm(b);  t0 = tic;
active = true(nw, 1);
while any(active)
    V = zeros(m, mr + 1);  Hh = zeros(mr + 1, mr);  V(:, 1) = v;
    for k = 1:mr
        w = A*V(:, k);  np = np + 1;
        for pass = 1:2, h = V'*w;  h(k + 1:end) = 0;  w = w - V*h;  Hh(1:k, k) = Hh(1:k, k) + h(1:k); end   % no copies
        Hh(k + 1, k) = norm(w);  V(:, k + 1) = w/Hh(k + 1, k);
    end
    for j = find(active).'
        y = (Hh(1:mr, 1:mr) + 1i*om(j)*eye(mr))\(bet(j)*eye(mr, 1));
        X(:, j) = X(:, j) + V(Ci, 1:mr)*y;  bet(j) = -Hh(mr + 1, mr)*y(mr);
        if abs(bet(j))/nb <= tol, active(j) = false; end
    end
    v = V(:, mr + 1);
    if toc(t0) > tcap, error('time cap with %d of %d shifts unconverged', nnz(active), nw); end
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

function [Y, np, k, I] = bicg_outputs(A, b, sig, tol, Ci, tcap)
% Multishift BiCG for (A + sig_j I) x_j = b, j = 1..nw, keeping only Y = X(Ci, :).
% Seed: BiCG on the real system A x = b (shadow residual b), all seed quantities real.  The residual polynomials satisfy
%   R_{k+1}(t) = (1 + c_k - al_k t) R_k(t) - c_k R_{k-1}(t),   c_k = al_k be_{k-1}/al_{k-1},
% and the shifted residuals are collinear with the seed's:  r^s_k = r_k/pi_k,  pi_k = R_k(-sig) (so R^s_k(0) = 1).  Then
%   al^s_k = al_k pi_k/pi_{k+1},  be^s_k = (pi_k/pi_{k+1})^2 be_k,  p^s_{k+1} = r_{k+1}/pi_{k+1} + be^s_k p^s_k,
% and since only C x^s is wanted, only C p^s (ell x nw) is kept.  A frequency stops once ||r_k||/|pi_k| <= tol ||b||
% (frozen before |pi| can overflow); the seed itself (sig = 0, the hardest system) need not converge.
t0 = tic;  At = A.';  sig = sig(:).';  nw = numel(sig);  nb = norm(b);
r = b;  rt = b;  p = r;  pt = rt;  rho = rt'*r;
pim = ones(1, nw);  pi0 = ones(1, nw);  alold = 1;  beold = 0;
cx = complex(zeros(numel(Ci), nw));  cp = complex(repmat(r(Ci), 1, nw));
act = true(1, nw);  k = 0;
while any(act)
    k = k + 1;
    Ap = A*p;  Atp = At*pt;  den = pt'*Ap;
    if den == 0 || rho == 0, error('bicg_outputs:breakdown', 'BiCG breakdown at step %d', k); end
    al = rho/den;  c = al*beold/alold;
    pin = pi0;  pin(act) = (1 + c + al*sig(act)).*pi0(act) - c*pim(act);           % pi_{k+1}, active shifts only
    J = find(act);
    cx(:, J) = cx(:, J) + (al*pi0(J)./pin(J)).*cp(:, J);                          % C x^s += al^s C p^s
    r = r - al*Ap;  rt = rt - al*Atp;  rhon = rt'*r;  be = rhon/rho;
    cp(:, J) = r(Ci)./pin(J) + ((pi0(J)./pin(J)).^2*be).*cp(:, J);               % C p^s_{k+1}
    p = r + be*p;  pt = rt + be*pt;
    act(J) = norm(r)./abs(pin(J)) > tol*nb;                                        % shifted residual norms
    pim = pi0;  pi0 = pin;  alold = al;  beold = be;  rho = rhon;
    if mod(k, 500) == 0, fprintf('   bicg k = %d: %d of %d frequencies left, %.0f s\n', k, nnz(act), nw, toc(t0)); end
    if toc(t0) > tcap, break, end
end
Y = cx;  np = 2*k;  I = struct('nconv', nw - nnz(act), 'time', toc(t0));
end
