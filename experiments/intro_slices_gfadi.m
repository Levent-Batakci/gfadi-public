function intro_slices_gfadi(N, omegas, check)
%INTRO_SLICES_GFADI  Example 4: the solution on the plane z = 1/2 through the source (Fig. 10(b,c)).
%   Computed with one-sided g-fADI (products with A only), so that the slice can be drawn on a fine grid.
%   Same operator as FREQ_SWEEP_RUN with C = 1e3, W = 1: 27 inclusions of conductivity 1e3 in a = 1, drift |w| = 1,
%   Dirichlet; (A + i w I) x = b for each w in OMEGAS, i.e. A X - X diag(-i w) = b 1^T.  The shifts are those of the
%   sweep: ADV_SHIFTS on the region Omega_A containing W(A), with F = {0}, bound tol/4, Leja order.
%
%   INTRO_SLICES_GFADI(N, OMEGAS, CHECK)
%     N       grid points per direction; must be odd, so that z = 1/2 is a grid plane (h = 1/(N+1))
%     OMEGAS  frequencies (default [3e2 3e3 3e4])
%     CHECK   true: also solve directly and report the relative error on the slice (small N only; default false)
%   Writes results/intro_slices_N<N>.mat (the slices and the settings), which MAKE_SWEEP_FIG reads.
%   Examples:  intro_slices_gfadi(31, [3e2 3e3 3e4], true)    intro_slices_gfadi(127)   (Fig. 10)
if nargin < 2 || isempty(omegas), omegas = [3e2 3e3 3e4]; end
if nargin < 3, check = false; end
assert(mod(N, 2) == 1, 'N must be odd so that z = 1/2 is a grid plane');
setup_paths();  addpath(fullfile(app_root(), 'testbed'));  maxNumCompThreads(8);
t0 = tic;
[A, b, E, c, ymax] = build_op(N);
m = N^3;  k = (N + 1)/2;                                   % z = k h = 1/2
ci = (k - 1)*N^2 + (1:N^2);                                % the plane z = 1/2 in the column-major (x, y, z) ordering
tol = 1e-8;
[p, al] = adv_shifts(E, c, ymax, [0 0], tol/4);  al = leja_order(al, E');
fprintf('N = %d (m = %d): E = [%.3g %.3g], p = %d steps, set-up %.1f s\n', N, m, E, p, toc(t0));
sB = -1i*omegas(:);                                        % A X - X diag(sB) = b 1^T  <=>  (A + i w I) x = b
t1 = tic;
Ys = gfadi_outputs(A, b, sB, al, ci);
tg = toc(t1);
fprintf('g-fADI: %.1f s\n', tg);
U = reshape(Ys, N, N, numel(omegas));                       % U(:, :, j) = u(x, y, 1/2) at omegas(j)
err = NaN(1, numel(omegas));
if check
    for j = 1:numel(omegas)
        x = (A + 1i*omegas(j)*speye(m))\b;  ud = reshape(x(ci), N, N);
        err(j) = norm(ud(:) - reshape(U(:, :, j), [], 1))/norm(ud(:));
        fprintf('omega = %g: relative error on the slice %.1e\n', omegas(j), err(j));
    end
end
f = fullfile(app_root(), 'results', sprintf('intro_slices_N%d.mat', N));
save(f, 'N', 'omegas', 'U', 'p', 'tol', 'tg', 'err');
fprintf('wrote %s\n', f);
end

function [A, b, E, c, ymax] = build_op(N)
% FREQ_SWEEP_RUN's operator with W = 1, C = 1e3 (copied, so this file is self-contained)
W = 1;  C = 1e3;
h = 1/(N + 1);  g = (1:N)'*h;  e = ones(N, 1);  I = speye(N);
D = spdiags([-e e], [-1 0], N + 1, N)/h;
gf = ((0:N)' + 0.5)*h;
ctr = 0.25*[kron([1;1;1], kron([1;1;1], (1:3)')), kron([1;1;1], kron((1:3)', [1;1;1])), kron((1:3)', kron([1;1;1], [1;1;1]))];
ctr = ctr + 0.04*sin((1:27)'*[1.3 2.1 3.7]);
af = @(x, y, z) 1 + (C - 1)*incl(x, y, z, ctr, 0.11);
[X, Y, Z] = ndgrid(gf, g, g);  Ax = af(X, Y, Z);  [X, Y, Z] = ndgrid(g, gf, g);  Ay = af(X, Y, Z);  [X, Y, Z] = ndgrid(g, g, gf);  Az = af(X, Y, Z);
Dx = kron(I, kron(I, D));  Dy = kron(I, kron(D, I));  Dz = kron(D, kron(I, I));
H = Dx'*spdiags(Ax(:), 0, numel(Ax), numel(Ax))*Dx + Dy'*spdiags(Ay(:), 0, numel(Ay), numel(Ay))*Dy + Dz'*spdiags(Az(:), 0, numel(Az), numel(Az))*Dz;
wd = W*[1 0.6 0.3]/norm([1 0.6 0.3]);  C1 = spdiags([-e 0*e e], -1:1, N, N)/(2*h);
S = wd(1)*kron(I, kron(I, C1)) + wd(2)*kron(I, kron(C1, I)) + wd(3)*kron(C1, kron(I, I));
A = H + S;
[X, Y, Z] = ndgrid(g, g, g);  b = exp(-((X - 0.3).^2 + (Y - 0.4).^2 + (Z - 0.5).^2)/0.01);  b = b(:)/norm(b(:));
E = [3*(4/h^2)*sin(pi*h/2)^2, full(max(sum(abs(H), 2)))];  c = W;  ymax = sum(abs(wd))/h;
end

function s = incl(x, y, z, ctr, r)
s = false(size(x));
for i = 1:size(ctr, 1), s = s | ((x - ctr(i, 1)).^2 + (y - ctr(i, 2)).^2 + (z - ctr(i, 3)).^2 < r^2); end
end

function Y = gfadi_outputs(A, b, sB, al, Ci)
% FREQ_SWEEP_RUN's output-only one-sided g-fADI (copied), with a progress line every 10% of the steps
Zt = b;  yt = ones(numel(sB), 1);  Y = zeros(numel(Ci), numel(sB));  dB = conj(sB(:));
nb = 256;  Zb = zeros(numel(Ci), nb);  Yb = complex(zeros(numel(sB), nb));  ib = 0;
p = numel(al);  t = tic;
for j = 1:p
    a = al(j);
    yj = yt./(dB - conj(a));
    ib = ib + 1;  Zb(:, ib) = Zt(Ci);  Yb(:, ib) = yj;
    if ib == nb || j == p, Y = Y - Zb(:, 1:ib)*Yb(:, 1:ib)';  ib = 0; end
    Zt = A*Zt - a*Zt;
    s = norm(Zt);  Zt = Zt/s;  yt = yj*s;
    if mod(j, ceil(p/10)) == 0, fprintf('  step %d / %d, %.0f s\n', j, p, toc(t)); end
end
end
