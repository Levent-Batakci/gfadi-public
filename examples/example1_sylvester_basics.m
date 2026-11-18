%EXAMPLE1_SYLVESTER_BASICS  Solve a small Sylvester equation A X - X B = U V' with g-fADI.
%
%   A is the 5-point Laplacian on a 30 x 30 grid (m = 900), B = -L' with L a 1-D diffusion
%   matrix (n = 50), a small version of the basic test case of the paper (Example 1, in two
%   dimensions). We pick shifts for a target error of 1e-8 and compare three types of
%   iteration against MATLAB's dense SYLVESTER:
%     q = 0   one-sided: right steps only (products with A, shifted solves with the small B)
%     q = 3   lopsided: a few left steps (shifted solves with A)
%     q = p   classical fADI: one left step per right step
%   Runs in a few seconds.

setup_paths();
rng(1);

% ---- the problem
N = 30;  h = 1/(N + 1);  e = ones(N, 1);
T = spdiags([-e 2*e -e], -1:1, N, N)/h^2;          % 1-D second-order difference matrix
A = kron(speye(N), T) + kron(T, speye(N));          % 2-D Laplacian, m = N^2
n = 50;  en = ones(n, 1);
L = speye(n) + (1e3 - 1)/4*spdiags([-en 2*en -en], -1:1, n, n);   % spectrum in [1, 1e3]
B = -L';
U = randn(N^2, 2);  V = randn(n, 2);                % rank-2 right-hand side

% ---- the two spectral intervals the shifts are designed for
% E contains the spectrum of A (known exactly here), F contains the spectrum of B.
E = 2*4/h^2*[sin(pi*h/2)^2, cos(pi*h/2)^2];
eL = eig(full(L));
F = [-max(eL), -min(eL)];

% ---- reference solution
X = sylvester(full(A), -full(B), U*V');

% ---- three types of g-fADI
tol = 1e-8;
fprintf('%-18s %6s %6s %12s\n', 'iteration', 'p', 'q', 'rel. error');
for q = [0 3 NaN]
    if isnan(q)
        [~, al, be] = classical_shifts_real(E, F, tol);  % type (k,k): k zeros and k poles
        name = 'classical fADI';
    else
        [~, al, be] = shift_search_real(E, F, q, tol);   % type (p,q): smallest p for this q
        name = sprintf('g-fADI, q = %d', q);
    end
    al = leja_order(al, E');                          % a stable order of the shifts
    be = leja_order(be, F');
    [Z, Y] = gfadi(A, B, U, V, al, be);               % X ~ Z*Y'
    err = norm(Z*Y' - X, 'fro')/norm(X, 'fro');
    fprintf('%-18s %6d %6d %12.2e\n', name, numel(al), numel(be), err);
end
fprintf('\nEvery right step costs a product with A; every left step a shifted solve with A.\n');
