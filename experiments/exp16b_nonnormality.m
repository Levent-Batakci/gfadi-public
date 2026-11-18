function R = exp16b_nonnormality(N, n, epslist, tol)
%EXP16B_NONNORMALITY  Example 2 at fixed size with an increasingly nonnormal B (Fig. 7(b)).
%   The problem of EXP16_BLOWOUT (BUILD_AGE3D).  Smaller eps means a more drift-dominated, more nonnormal L
%   (eigenvector condition number ~ exp(v/eps)).  For each eps: one-sided g-fADI (shifts from Re W(L)), the
%   eigendecomposition route, classical fADI (shifts for the real intervals), and the stable direct method
%   SYLV_SCHUR_DIRECT (the reference).
%
%   R = EXP16B_NONNORMALITY(N, n, EPSLIST, TOL)
%     N (32), n (400), EPSLIST ([0.1 0.05 0.02 0.01 0.005 0.003 0.002]), TOL (1e-8)
%   Writes results/exp16b_N<N>.mat (R) and results/exp16b_N<N>.csv.  Fig. 7(b) is EXP16B_NONNORMALITY() with
%   the defaults, read from results/exp16b_N32.csv by MAKE_PAPER_FIGS('bound').
if nargin < 1, N = 32; end
if nargin < 2, n = 400; end
if nargin < 3, epslist = [0.1 0.05 0.02 0.01 0.005 0.003 0.002]; end
if nargin < 4, tol = 1e-8; end
setup_paths();
here = app_root();
Rc = {};
for ep = epslist
    P = build_age3d(N, n, struct('eps', ep));
    m = P.m;  U = P.U;  V = P.V;
    Aop = gfadi_op(P.A);  Bop = gfadi_op(-P.L.');
    Eint = [P.lminA P.lmaxA];  Fint = [-P.reLmax -P.reLmin];
    % reference: Bartels-Stewart (complex Schur form when rounding gives complex eigenvalue pairs)
    F = U*V.';  rr = @(Xa) norm(P.A*Xa + Xa*P.L.' - F, 'fro')/norm(F, 'fro');
    [X, idr] = sylv_schur_direct(P.A, P.L, U, V);
    t_direct = idr.time;
    nX = norm(X, 'fro');
    res_direct = rr(X);
    % g-fADI
    [p, al] = shift_search_real(Eint, Fint, 0, tol);
    al = leja_order(al, Eint');
    [Zc, Yc, ig] = gfadi_lowrank(Aop, Bop, U, V, al, []);
    eg = norm(Zc*Yc' - X, 'fro')/nX;
    % eig route
    [W, D] = eig(full(P.L.'));  lam = diag(D);  Fw = U*(V.'*W);  Xw = zeros(m, n);
    for j = 1:n, Xw(:, j) = (P.A + lam(j)*speye(m))\Fw(:, j); end
    ee = norm(real(Xw/W) - X, 'fro')/nX;
    kW = cond(W);
    clear Xw
    % classical
    [kc, alc, bec] = classical_shifts_real(Eint, Fint, tol);
    [alc, bec, scc] = order_shifts(alc, bec, Eint, Fint);
    [Zk, Yk] = gfadi(Aop, Bop, U, V, alc, bec, 'sched', scc);
    ec = norm(Zk*Yk' - X, 'fro')/nX;
    fprintf(['eps = %.3f: cond(W) %.1e | g-fADI p = %d, %.2f s, error %.1e | eig route error %.1e | classical (%d,%d) error %.1e', ...
             ' | direct %.0f s (complex Schur %d), residual %.1e\n'], ep, kW, p, ig.time, eg, ee, kc, kc, ec, t_direct, idr.complex, res_direct);
    Rc{end+1} = struct('eps', ep, 'kappaW', kW, 'p', p, 't_gfadi', ig.time, 'err_gfadi', eg, ...
                       'err_eig', ee, 'k_classical', kc, 'err_classical', ec, 't_direct', t_direct, ...
                       'direct_complex', idr.complex, 'res_direct', res_direct);                       %#ok<AGROW>
end
R = [Rc{:}];  disp(struct2table(R));
save(fullfile(here, 'results', sprintf('exp16b_N%d.mat', N)), 'R');
writetable(struct2table(R), fullfile(here, 'results', sprintf('exp16b_N%d.csv', N)));
end
