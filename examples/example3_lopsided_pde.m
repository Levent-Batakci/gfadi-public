%EXAMPLE3_LOPSIDED_PDE  A lopsided 3-D PDE: A is large and sparse, B is small and nonnormal.
%
%   The problem of Example 2 of the paper (Section 4.2) at a small size: the steady state of a population that diffuses
%   in x (heterogeneous coefficient, 7-point stencil) and drifts and disperses in a structure
%   variable t. The discretization is A X - X B = U V' with A of size N^3 and B = -L' of size n
%   (BUILD_AGE3D). One-sided g-fADI uses products with A, tridiagonal solves with B, and
%   compresses the solution on the small side (GFADI_LOWRANK), so its memory does not grow
%   with the number of steps. Runs in under a minute.

setup_paths();

N = 16;  n = 100;
P = build_age3d(N, n);                               % fields used: A, L, U, V, lminA, lmaxA, reLmin, reLmax

% ---- shifts: E bounds the spectrum of A, F bounds Re W(B) (B = -L' is far from normal)
E = [P.lminA, P.lmaxA];
F = [-P.reLmax, -P.reLmin];
[p, al] = shift_search_real(E, F, 0, 1e-8);
al = leja_order(al, E');

% ---- one-sided g-fADI with compression on the small side
[Z, Y, info] = gfadi_lowrank(gfadi_op(P.A), gfadi_op(-P.L.'), P.U, P.V, al, []);

% ---- check the residual (small enough here to form X)
X = Z*Y';
res = norm(P.A*X + X*P.L.' - P.U*P.V', 'fro')/norm(P.U*P.V', 'fro');
fprintf('N = %d (m = %d), n = %d: %d steps, residual %.1e\n', N, N^3, n, p, res);
fprintf('the solution is stored as %d columns of length m, not %d\n', info.s, 2*p);
