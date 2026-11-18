%EXAMPLE2_SHIFTED_SYSTEMS  Many shifted systems (A + s_j I) x_j = b with products by A only.
%
%   The k systems are one Sylvester equation, A X - X B = b 1', with B = diag(-s).
%   One-sided g-fADI then needs one product with A per step and divisions by the s_j,
%   and no shifted solves with A at all. This is the reformulation of Example 4 of the paper
%   (Section 4.4), here with a symmetric A and real shifts. We compare with k sparse direct
%   solves. Runs in about ten seconds.

setup_paths();
rng(2);

% ---- the problem: 2-D Laplacian on a 150 x 150 grid, 100 shifts over its spectrum
N = 150;  h = 1/(N + 1);  e = ones(N, 1);
T = spdiags([-e 2*e -e], -1:1, N, N)/h^2;
A = kron(speye(N), T) + kron(T, speye(N));
m = N^2;  k = 100;
E = 2*4/h^2*[sin(pi*h/2)^2, cos(pi*h/2)^2];         % spectrum of A
s = logspace(log10(E(1)), log10(E(2)), k)';         % the shifts
b = randn(m, 1);  b = b/norm(b);

% ---- shifts for g-fADI: zeros on E, chosen so that the error bound is 1e-8
F = [-max(s), -min(s)];                              % spectrum of B = diag(-s)
[p, al] = shift_search_real(E, F, 0, 1e-8);          % q = 0: Chebyshev points of E
al = leja_order(al, E');

% ---- one-sided g-fADI: A appears only through products
Aop = gfadi_op(@(x) A*x, 'free');                    % no solves with A are possible
Bop = gfadi_op(-s, 'diag');
t = tic;
[Z, Y, info] = gfadi_lowrank(Aop, Bop, b, ones(k, 1), al, []);
tg = toc(t);
X = Z*Y';                                            % column j is x_j

% ---- k direct solves, for comparison
t = tic;
Xd = zeros(m, k);
for j = 1:k, Xd(:, j) = (A + s(j)*speye(m))\b; end
td = toc(t);

err = max(sqrt(sum(abs(X - Xd).^2, 1))./sqrt(sum(abs(Xd).^2, 1)));
fprintf('g-fADI: %d steps (products with A), %.2f s, compressed rank %d\n', p, tg, info.s);
fprintf('direct: %d sparse solves, %.2f s\n', k, td);
fprintf('largest relative error over the %d systems: %.1e\n', k, err);
