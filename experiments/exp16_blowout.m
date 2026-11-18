function R = exp16_blowout(Nlist, n, tol, opts)
%EXP16_BLOWOUT  Example 2 (structured population in 3-D): cost of g-fADI against the alternatives (Fig. 8).
%   A X + X L^T = U V^T from BUILD_AGE3D: N^3 x n unknowns, A the SPD 7-point matrix of -div(kappa grad),
%   L tridiagonal and nonnormal, rank-2 right-hand side.  Methods: one-sided g-fADI with compression on the
%   small side (GFADI_LOWRANK), Bartels-Stewart with n shifted solves (SYLV_SCHUR_DIRECT, also the reference),
%   the eigendecomposition route, classical fADI, and Galerkin projection onto the g-fADI spaces.  The
%   expensive methods run only up to the size caps in OPTS; one shifted solve is timed up to maxfact.
%
%   R = EXP16_BLOWOUT(NLIST, n, TOL, OPTS)
%     NLIST  grid sizes (default [24 32 48 64 96 128]),  n  size of L (400),  TOL  target (1e-8)
%     OPTS   .maxdirect (32) .maxgal (48) .maxclass (48) .maxfact (96), .noeig (skip the eig route),
%            .tag (suffix of the output files, e.g. '_gfadi', '_ext', '_ext96', '_ext40')
%   Writes results/exp16_blowout<tag>.mat (R, saved after every N) and results/exp16_blowout<tag>.csv.
%   Fig. 8 (MAKE_PAPER_FIGS('cost')) reads the untagged run and merges four tagged runs, each in its own
%   MATLAB session.  The options of the '_ext' and '_ext96' runs are reconstructed from the methods they
%   measured; the other calls are as run:
%     exp16_blowout()                                                            % untagged base run
%     exp16_blowout([24 32 48 64 96 128], 400, 1e-8, struct('maxdirect', 32, 'maxgal', 0, 'maxclass', 0, ...
%                   'maxfact', 0, 'tag', '_gfadi'))                              % g-fADI times and memory
%     exp16_blowout([48 64], 400, 1e-8, struct('maxdirect', 64, 'maxgal', 64, 'maxclass', 64, 'tag', '_ext'))
%     exp16_blowout(96, 400, 1e-8, struct('maxdirect', 0, 'maxgal', 0, 'maxclass', 96, 'tag', '_ext96'))
%     exp16_blowout([16 40], 400, 1e-8, struct('maxdirect', 40, 'maxgal', 40, 'maxclass', 40, 'tag', '_ext40'))
if nargin < 1, Nlist = [24 32 48 64 96 128]; end
if nargin < 2, n = 400; end
if nargin < 3, tol = 1e-8; end
if nargin < 4, opts = struct(); end
setup_paths();
here = app_root();
maxdirect = 32;  maxgal = 48;  maxclass = 48;  maxfact = 96;
if isfield(opts, 'maxdirect'), maxdirect = opts.maxdirect; end
if isfield(opts, 'maxgal'), maxgal = opts.maxgal; end
if isfield(opts, 'maxclass'), maxclass = opts.maxclass; end
if isfield(opts, 'maxfact'), maxfact = opts.maxfact; end
tag = '';
if isfield(opts, 'tag'), tag = opts.tag; end       % e.g. '_gfadi' for a g-fADI-only rerun
Rc = {};
for N = Nlist
    P = build_age3d(N, n);
    m = P.m;  U = P.U;  V = P.V;  rho = size(U, 2);
    Aop = gfadi_op(P.A);  Bop = gfadi_op(-P.L.');
    Eint = [P.lminA P.lmaxA];  Fint = [-P.reLmax -P.reLmin];
    fprintf('\n=== N = %d: m = %d, n = %d, %.2e unknowns;  spec(A) in [%.3g, %.3g], Re W(L) in [%.3g, %.3g]\n', ...
            N, m, n, m*n, Eint, P.reLmin, P.reLmax);
    r = struct('N', N, 'm', m, 'n', n, 'unknowns', m*n);
    % ---------------------------------------------------------------- g-fADI
    [p, al, ~, xi] = shift_search_real(Eint, Fint, 0, tol);
    al = leja_order(al, Eint');
    [Zc, Yc, ig] = gfadi_lowrank(Aop, Bop, U, V, al, []);
    resg = lowrank_residual(P, Zc, Yc);
    r.p = p;  r.xi = xi;  r.t_gfadi = ig.time;  r.t_comp = ig.tcomp;  r.rank = ig.rank;  r.s = ig.s;
    r.mem_gfadi_GB = (ig.maxcols + 4*rho)*m*8/2^30;  r.res_gfadi = resg;
    fprintf('g-fADI     : p = %d (xi %.1e), %.1f s, basis %d (numerical rank %d), memory %.2f GB, residual %.1e\n', ...
            p, xi, ig.time, ig.s, r.rank, r.mem_gfadi_GB, resg);
    % ---------------------------------------------------------------- unit cost of one shifted solve (Cholesky of A + sI)
    r.t_fact = NaN;
    if N <= maxfact
        sh = sqrt(P.reLmin*P.reLmax);
        tf = tic;  [Lc, ~, Pp] = chol(P.A + sh*speye(m), 'lower', 'vector');  r.t_fact = toc(tf);
        xs = zeros(m, 1);  tb = tic;  xs(Pp) = Lc'\(Lc\U(Pp, 1));  r.t_backsolve = toc(tb);
        r.fill_GB = nnz(Lc)*16/2^30;  clear Lc
        fprintf('one sparse factorisation of A + sI: %.2f s (+%.2f s solve), factor %.2f GB\n', r.t_fact, r.t_backsolve, r.fill_GB);
    else
        r.t_backsolve = NaN;  r.fill_GB = NaN;
    end
    % ---------------------------------------------------------------- direct (Schur) and eig
    r.err_gfadi = NaN;  r.t_direct = NaN;  r.err_eig = NaN;  r.kappaW = NaN;
    if N <= maxdirect
        [X, idr] = sylv_schur_direct(P.A, P.L, U, V);
        r.t_direct = idr.time;
        nX = norm(X, 'fro');
        r.err_gfadi = norm(Zc*Yc' - X, 'fro')/nX;
        fprintf('direct     : %.1f s (n = %d factorisations);  g-fADI error %.2e\n', r.t_direct, n, r.err_gfadi);
        % eig route (exact solves, to isolate the effect of the eigenvector basis); opts.noeig skips it
        if isfield(opts, 'noeig') && opts.noeig
            clear X;
        else
            [W, D] = eig(full(P.L.'));  lam = diag(D);  r.kappaW = cond(W);
            Fw = U*(V.'*W);  Xw = zeros(m, n);
            for j = 1:n, Xw(:, j) = (P.A + lam(j)*speye(m))\Fw(:, j); end
            Xe = real(Xw/W);  clear Xw
            r.err_eig = norm(Xe - X, 'fro')/nX;
            fprintf('eig route  : cond(W) = %.1e, error %.2e\n', r.kappaW, r.err_eig);
            clear Xe
        end
    end
    % ---------------------------------------------------------------- classical fADI
    [kc, alc, bec] = classical_shifts_real(Eint, Fint, tol);
    r.k_classical = kc;
    if N <= maxclass
        [alc, bec, scc] = order_shifts(alc, bec, Eint, Fint);
        tc = tic;  [Zk, Yk, icl] = gfadi(Aop, Bop, U, V, alc, bec, 'sched', scc);  r.t_classical = toc(tc);
        r.res_classical = lowrank_residual(P, Zk, Yk);  r.classical_measured = true;
        if exist('X', 'var') && N <= maxdirect, r.err_classical = norm(Zk*Yk' - X, 'fro')/nX; else, r.err_classical = NaN; end
        clear Zk Yk
    else
        r.t_classical = NaN;  r.res_classical = NaN;  r.err_classical = NaN;  r.classical_measured = false;
    end
    fprintf('classical  : k = %d factorisations, %.1f s%s\n', kc, r.t_classical, ifelse(r.classical_measured, '', ' (not run)'));
    % ---------------------------------------------------------------- Galerkin on the same spaces
    r.k_gal = NaN;  r.t_gal = NaN;  r.mem_gal_GB = NaN;  r.res_gal = NaN;
    if N <= maxgal
        target = max(resg, 1e-12);                       % smallest k reaching the g-fADI residual (doubling, then bisection)
        resk = @(k) gal_res(Aop, Bop, U, V, al, k, P);
        lo = 0;  hi = min(p, 200);
        while resk(hi) > target && hi < p, lo = hi; hi = min(p, 2*hi); end
        while hi - lo > max(1, round(0.05*hi))           % bisect to within 5%
            mid = floor((lo + hi)/2);
            if resk(mid) <= target, hi = mid; else, lo = mid; end
        end
        [QA, QB, Yg, igal] = sylvester_galerkin(Aop, Bop, U, V, inf(1, hi), al(1:hi).', 'residual', false);
        r.k_gal = hi;  r.t_gal = igal.time;  r.mem_gal_GB = 2*size(QA, 1)*size(QA, 2)*8/2^30;
        r.res_gal = lowrank_residual(P, QA*Yg, QB);
        fprintf('Galerkin   : k = %d block iterations (%d vectors), %.1f s, basis memory %.2f GB, residual %.1e\n', ...
                hi, size(QA, 2), r.t_gal, r.mem_gal_GB, r.res_gal);
        clear QA QB
    end
    Rc{end+1} = r;                                                                         %#ok<AGROW>
    R = [Rc{:}];  save(fullfile(here, 'results', ['exp16_blowout' tag '.mat']), 'R');
    clear P X Zc Yc
end
R = [Rc{:}];  T = struct2table(R);  disp(T)
writetable(T, fullfile(here, 'results', ['exp16_blowout' tag '.csv']));
end

function res = lowrank_residual(P, Z, Y)
% || A Z Y' + Z Y' L' - U V' ||_F / || U V' ||_F, exactly, from thin QR factors
Lf = [P.A*Z, Z, -P.U];  Rf = [Y, P.L*Y, P.V];
[~, R1] = qr(Lf, 0);  [~, R2] = qr(Rf, 0);
res = norm(R1*R2', 'fro')/sqrt(abs(trace((P.U'*P.U)*(P.V'*P.V))));
end

function r = gal_res(Aop, Bop, U, V, al, k, P)
[QA, QB, Y] = sylvester_galerkin(Aop, Bop, U, V, inf(1, k), al(1:k).', 'residual', false);
r = lowrank_residual(P, QA*Y, QB);
end

function y = ifelse(c, a, b)
if c, y = a; else, y = b; end
end

