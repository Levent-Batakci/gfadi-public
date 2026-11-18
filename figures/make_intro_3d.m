function make_intro_3d(N, omega)
%MAKE_INTRO_3D  3-D render of the composite material of Example 4 (Fig. 10(a) of the paper).
%   Draws the 27 inclusions of conductivity 1e3 in a = 1 (the geometry of FREQ_SWEEP_RUN), the heat
%   source and the 100 sensors.  MAKE_SWEEP_FIG reads the PNG and places it as panel (a) of Fig. 10,
%   so run this first.
%
%   Call:   make_intro_3d()  or  make_intro_3d(N, omega)   (defaults N = 48, omega = 30; the picture
%           shows the set-up only, so omega does not affect it)
%   Reads:  nothing (everything is computed here).
%   Writes: figs/paper_intro_3d.{pdf,png}.
if nargin < 1, N = 48; end
if nargin < 2, omega = 30; end
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultAxesFontSize', 8);

% ---- the operator (as in FREQ_SWEEP_RUN, C = 1e3, |w| = 1)
C = 1e3;
W = 1;
h = 1/(N + 1);
g = (1:N)'*h;
e = ones(N, 1);
I = speye(N);
D = spdiags([-e e], [-1 0], N + 1, N)/h;
gf = ((0:N)' + 0.5)*h;
ctr = 0.25*[kron([1;1;1], kron([1;1;1], (1:3)')), kron([1;1;1], kron((1:3)', [1;1;1])), kron((1:3)', kron([1;1;1], [1;1;1]))];
ctr = ctr + 0.04*sin((1:27)'*[1.3 2.1 3.7]);
r0 = 0.11;
af = @(x, y, z) 1 + (C - 1)*incl(x, y, z, ctr, r0);
[X, Y, Z] = ndgrid(gf, g, g);
Ax = af(X, Y, Z);
[X, Y, Z] = ndgrid(g, gf, g);
Ay = af(X, Y, Z);
[X, Y, Z] = ndgrid(g, g, gf);
Az = af(X, Y, Z);
Dx = kron(I, kron(I, D));
Dy = kron(I, kron(D, I));
Dz = kron(D, kron(I, I));
H = Dx'*spdiags(Ax(:), 0, numel(Ax), numel(Ax))*Dx + Dy'*spdiags(Ay(:), 0, numel(Ay), numel(Ay))*Dy + Dz'*spdiags(Az(:), 0, numel(Az), numel(Az))*Dz;
wd = W*[1 0.6 0.3]/norm([1 0.6 0.3]);
C1 = spdiags([-e 0*e e], -1:1, N, N)/(2*h);
A = H + wd(1)*kron(I, kron(I, C1)) + wd(2)*kron(I, kron(C1, I)) + wd(3)*kron(C1, kron(I, I));
[X, Y, Z] = ndgrid(g, g, g);
x0 = [0.3 0.4 0.5];
b = exp(-((X - x0(1)).^2 + (Y - x0(2)).^2 + (Z - x0(3)).^2)/0.01);
b = b(:)/norm(b(:));
m = N^3;
rng(3);
Ci = sort(randperm(m, 100));
[si, sj, sk] = ind2sub([N N N], Ci);
S = [g(si) g(sj) g(sk)];

% ---- the picture (the set-up only, no solution)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.2 3.6], 'Visible', 'off');
ax = axes(fig, 'Position', [0.02 0.04 0.86 0.92]);
hold(ax, 'on');
% the inclusions: smooth spheres, light grey
[sx, sy, sz] = sphere(40);
for i = 1:27
    surf(ax, ctr(i, 1) + r0*sx, ctr(i, 2) + r0*sy, ctr(i, 3) + r0*sz, 'FaceColor', [0.86 0.86 0.88], 'EdgeColor', 'none', ...
         'FaceAlpha', 0.45, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.6, 'DiffuseStrength', 0.6, 'SpecularStrength', 0.4);
end
% the source: a small orange ball; the sensors: black dots
surf(ax, x0(1) + 0.04*sx, x0(2) + 0.04*sy, x0(3) + 0.04*sz, 'FaceColor', [0.91 0.45 0.13], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
plot3(ax, S(:, 1), S(:, 2), S(:, 3), 'o', 'Color', [0.15 0.15 0.15], 'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 1.6);
% the unit cube
cube = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1; 0 0 1];
plot3(ax, cube(:, 1), cube(:, 2), cube(:, 3), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
for c = [1 0; 1 1; 0 1]', plot3(ax, [c(1) c(1)], [c(2) c(2)], [0 1], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5); end
axis(ax, 'equal');
axis(ax, 'off');
view(ax, [-35 30]);
camlight(ax, 'headlight');
lighting(ax, 'gouraud');
xlim(ax, [0 1]);
ylim(ax, [0 1]);
zlim(ax, [0 1]);
set(ax, 'Position', [0.02 0.02 0.96 0.88]);
% the panel letter (MAKE_SWEEP_FIG crops it off and sets its own)
annotation(fig, 'textbox', [0 0.92 1 0.08], 'String', '(a)', 'Interpreter', 'latex', 'FontSize', 20, ...
           'HorizontalAlignment', 'center', 'EdgeColor', 'none');
f = fullfile(app_root(), 'figs', 'paper_intro_3d');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
exportgraphics(fig, [f '.pdf'], 'ContentType', 'image', 'Resolution', 500);
close(fig);
fprintf('intro 3d: set-up picture written\n');
end

function s = incl(x, y, z, ctr, r)
% indicator (as double) of the union of the balls |x - ctr(i,:)| < r
s = false(size(x));
for i = 1:size(ctr, 1), s = s | ((x - ctr(i, 1)).^2 + (y - ctr(i, 2)).^2 + (z - ctr(i, 3)).^2 < r^2); end
s = double(s);
end

function c = infer_map(n)
% a dark-to-light blue-to-orange map (inferno-like, perceptually ordered)
k = [0 0 0.02 0.10; 0.25 0.26 0.07 0.43; 0.5 0.73 0.21 0.33; 0.75 0.98 0.55 0.04; 1 0.99 1.00 0.64];
c = interp1(k(:, 1), k(:, 2:4), linspace(0, 1, n));
end
