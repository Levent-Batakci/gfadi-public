function R = exp18_analysis(which)
%EXP18_ANALYSIS  The field-of-values error bound for one-sided g-fADI, checked against experiment (Fig. 7(a)).
%   For A X - X B = U V^T, A SPD with spectrum in E, W(B) in Re z <= -f, a one-sided run satisfies
%   ||X_k - X|| <= max_E|r_k| / |r_k(-f)| * ||X||, r_k(z) = prod_{j<=k} (z - alpha_j)  (the bound of Example 2).
%   Along the run it records the true error, this bound, its sharpness and the normal-theory estimate; for
%   classical fADI the ratio on Re W(B) against W(B); and the singular values of X against the rank bound.
%
%   R = EXP18_ANALYSIS(WHICH)
%     WHICH  'age' (the Example 2 problem of EXP16_BLOWOUT, N = 24, eps = 0.02 and 0.005),
%            'frac' (time-fractional problem of BUILD_FRAC3D, N = 24), or 'both' (default)
%   Writes results/exp18_analysis.mat: R (one entry per case) and F (classical fADI on the eps sweep of
%   EXP16B_NONNORMALITY, N = 32).  Fig. 7(a) (MAKE_PAPER_FIGS('bound')) draws R(1), the case eps = 0.02.
if nargin < 1, which = 'both'; end
setup_paths();
here = app_root();
cases = {};
if any(strcmp(which, {'age', 'both'}))
    for ep = [0.02 0.005]
        P = build_age3d(24, 400, struct('eps', ep));
        % exact eigenvalues of the tridiagonal Toeplitz L (the computed ones are not, cond(W) ~ 1e10-1e31)
        ht = 1/(P.n + 1);  dg = P.opts.mu + 2*ep/ht^2;  sp = (ep/ht^2)^2 - (P.opts.v/(2*ht))^2;
        lamL = dg - 2*sqrt(sp)*cos((1:P.n)'*pi*ht);
        cases{end+1} = struct('name', sprintf('age, eps = %g', ep), 'A', P.A, 'L', P.L, 'U', P.U, 'V', P.V, ...
                              'E', [P.lminA P.lmaxA], 'f', P.reLmin, 'fmax', P.reLmax, 'lamB', -lamL);  %#ok<AGROW>
    end
end
if any(strcmp(which, {'frac', 'both'}))
    P = build_frac3d(24, 400);
    cases{end+1} = struct('name', 'fractional, alpha = 0.5', 'A', P.A, 'L', P.T, 'U', P.U, 'V', P.V, ...
                          'E', [P.lminA P.lmaxA], 'f', P.reTmin, 'fmax', P.reTmax, 'lamB', -diag(P.T));
end
R = struct([]);
for c = 1:numel(cases)
    C = cases{c};  m = size(C.A, 1);  n = size(C.L, 1);  rho = size(C.U, 2);
    B = -full(C.L).';  E = C.E;  f = C.f;
    fprintf('\n=== %s: E = [%.3g, %.3g], Re W(B) <= -%.4g\n', C.name, E, f);
    % reference solution (complex Schur when needed)
    X = sylv_schur_direct(C.A, C.L, C.U, C.V);
    nX2 = norm(X);
    nXF = norm(X, 'fro');
    % ---- one-sided run, all factors kept
    [p, al] = shift_search_real(E, [-C.fmax -f], 0, 1e-8);
    al = leja_order(al, E');
    [Z, Y] = gfadi(gfadi_op(C.A), gfadi_op(B), C.U, C.V, al, []);
    xs = E(1) + (E(2) - E(1))*(1 - cos(pi*(0:20000)/20000))/2;      % dense in E, clustered at the ends
    lamB = C.lamB;                                                  % exact eigenvalues of B
    smin = min(svd(B + f*eye(n)));                                  % is -f an approximate eigenvalue?
    fprintf('eigenvalues of B in [%.4g, %.4g]; W(B) reaches Re z = %.4g; sigma_min(B + f I) = %.1e (||B|| = %.1e)\n', ...
            min(lamB), max(lamB), -f, smin, norm(B));
    % err2/errF: true error;  bnd: max_E|r_k| / |r_k(-f)|;  sharp: ||r_k(B)^{-1}|| |r_k(-f)| (<= 1);
    % eigp: max_E|r_k| / min_{lambda(B)}|r_k| (normal theory)
    D = X;  err2 = zeros(p, 1);  errF = zeros(p, 1);  bnd = zeros(p, 1);  eigp = zeros(p, 1);  sharp = zeros(p, 1);
    lw = zeros(size(xs));  lf = 0;  ll = zeros(size(lamB));  M = eye(n);  lM = 0;
    for k = 1:p
        D = D - Z(:, (k-1)*rho+1:k*rho)*Y(:, (k-1)*rho+1:k*rho)';
        if k <= 60 || mod(k, 5) == 0 || k == p, err2(k) = norm(D)/nX2; else, err2(k) = NaN; end
        errF(k) = norm(D, 'fro')/nXF;
        lw = lw + log(abs(xs - al(k)));  lf = lf + log(al(k) + f);  ll = ll + log(abs(lamB - al(k)));
        bnd(k) = exp(max(lw) - lf);  eigp(k) = exp(max(lw) - min(ll));
        M = (B - al(k)*eye(n))\M;  sc = norm(M);  lM = lM + log(sc);  M = M/sc;   % r_k(B)^{-1}, scale kept in log
        sharp(k) = exp(lM + lf);
    end
    fprintf('one-sided: p = %d; final error %.1e (2-norm), %.1e (F); bound %.1e; normal-theory estimate %.1e\n', ...
            p, err2(p), errF(p), bnd(p), eigp(p));
    fprintf('           max_k err/bnd = %.2f; ||r_p(B)^-1|| |r_p(-f)| = %.3f (min over k %.3f)\n', ...
            max(err2(~isnan(err2))./bnd(~isnan(err2))), sharp(p), min(sharp));
    % ---- classical fADI: designed (real F) vs field of values
    zW = fov_boundary(B, 720);
    [kc, alc, bec] = classical_shifts_real(E, [-C.fmax -f], 1e-8);
    [alc, bec, scc] = order_shifts(alc, bec, E, [-C.fmax -f]);
    [Zk, Yk] = gfadi(gfadi_op(C.A), gfadi_op(B), C.U, C.V, alc, bec, 'sched', scc);
    errc = norm(Zk*Yk' - X)/nX2;
    rE = @(z, a, b) sum(log(abs(z(:) - a(:).')), 2) - sum(log(abs(z(:) - b(:).')), 2);
    yR = -f - (C.fmax - f)*(1 - cos(pi*(0:20000)/20000))/2;
    xiR = exp(max(rE(xs, alc, bec)) - min(rE(yR, alc, bec)));
    xiW = exp(max(rE(xs, alc, bec)) - min(rE(zW, alc, bec)));
    fprintf('classical (%d,%d): error %.1e | ratio on Re W(B): %.1e | on W(B): %.1e (bound with Crouzeix factor 2: %.1e)\n', ...
            kc, kc, errc, xiR, xiW, 2*xiW);
    % ---- singular values of X against the rank bound from W(B)
    sv = svd(X);  ks = 1:16;  rb = zeros(size(ks));
    for i = ks
        [a1, b1] = opt_shifts_real(E, [-C.fmax -f], i, i);
        rb(i) = 2*exp(max(rE(xs, a1, b1)) - min(rE(zW, a1, b1)));
    end
    fprintf('sigma_{2k+1}/sigma_1 vs bound, k = 2,4,8,12,16: ');
    fprintf('%.1e/%.1e  ', [sv(2*ks([2 4 8 12 16]) + 1)'/sv(1); rb([2 4 8 12 16])]);
    fprintf('\n');
    R(c).name = C.name;  R(c).p = p;  R(c).err2 = err2;  R(c).errF = errF;  R(c).bnd = bnd;  R(c).eigp = eigp;
    R(c).sharp = sharp;  R(c).zW = zW;  R(c).lamB = lamB;  R(c).kc = kc;  R(c).errc = errc;  R(c).xiR = xiR;
    R(c).xiW = xiW;  R(c).sv = sv;  R(c).rankbound = rb;  R(c).E = E;  R(c).f = f;  R(c).smin = smin;
end
% ---- classical fADI on the EXP16B sweep (N = 32): designed ratio vs ratio on W(B)
F = struct([]);
if any(strcmp(which, {'age', 'both'}))
    for ep = [0.1 0.05 0.02 0.01 0.005]
        P = build_age3d(32, 400, struct('eps', ep));  E = [P.lminA P.lmaxA];  Fi = [-P.reLmax -P.reLmin];
        B = -full(P.L).';  zW = fov_boundary(B, 1440);
        [kc, alc, bec] = classical_shifts_real(E, Fi, 1e-8);
        xs = E(1) + (E(2) - E(1))*(1 - cos(pi*(0:20000)/20000))/2;
        yR = Fi(2) - (Fi(2) - Fi(1))*(1 - cos(pi*(0:20000)/20000))/2;
        rE = @(z) sum(log(abs(z(:) - alc(:).')), 2) - sum(log(abs(z(:) - bec(:).')), 2);
        F(end+1).eps = ep;  F(end).k = kc;                                                        %#ok<AGROW>
        F(end).xiR = exp(max(rE(xs)) - min(rE(yR)));  F(end).xiW = exp(max(rE(xs)) - min(rE(zW)));
        F(end).imW = max(abs(imag(zW)));
        fprintf('eps = %.3f: classical k = %d, ratio on Re W = %.1e, on W(B) = %.1e (max |Im W| = %.0f)\n', ...
                ep, kc, F(end).xiR, F(end).xiW, F(end).imW);
    end
end
save(fullfile(here, 'results', 'exp18_analysis.mat'), 'R', 'F');
end
