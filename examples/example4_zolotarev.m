%EXAMPLE4_ZOLOTAREV  Compute Zolotarev rational functions with ZOLO3 and look at the ratio they reach.
%
%   For two sets E and F, the Zolotarev number of type (p,q) is the smallest
%   max_E |r| / min_F |r| over rational r with p zeros and q poles. In g-fADI the zeros of r
%   are the right-step shifts (products with A) and the poles are the left-step shifts
%   (shifted solves with A), so a type (p,q) function with q < p needs fewer solves with A.
%   ZOLO3(E, F, N, 'delta', D) computes a function of type (N-D, N), with D fewer ZEROS.
%   To get D fewer POLES, swap E and F: 1/r is then small on E with N-D poles and N zeros,
%   and the ratio is the same. Runs in under a minute and draws one figure.

setup_paths();

% ---- two real intervals, sampled densely (no repeated points: AAA needs distinct samples)
E = unique([linspace(1, 100, 400)'; logspace(0, 2, 400)']);
F = -E;

% ---- type (k,k) against type (k, k-2) for k = 3..8
K = 3:8;  sig = zeros(numel(K), 2);
for i = 1:numel(K)
    o = zolo3(E, F, K(i));                           % type (k,k)
    sig(i, 1) = o.sigma;
    o = zolo3(F, E, K(i), 'delta', 2);               % swapped sets: type (k, k-2) for E, F
    sig(i, 2) = o.sigma;
end
disp(table(K', sig(:, 1), sig(:, 2), 'VariableNames', {'k', 'type_k_k', 'type_k_km2'}));

% ---- the closed form for q = 0 (one-sided): Chebyshev points of E
p = 20;
[al, xi] = cheb_shifts(1, 100, p, -100, -1);
fprintf('type (%d,0) by Chebyshev points: ratio %.2e\n', p, xi);

% ---- one function, drawn: |r| is small on E and large on F
o = zolo3(E, F, 6);
x = linspace(-120, 120, 2000)';
semilogy(x, abs(o.r3(x)), 'k-', E, abs(o.r3(E)), 'b.', F, abs(o.r3(F)), 'r.');
xlabel('x');  ylabel('|r(x)|');  title(sprintf('type (6,6), ratio %.1e', o.sigma));
legend('|r|', 'on E', 'on F', 'Location', 'best');
