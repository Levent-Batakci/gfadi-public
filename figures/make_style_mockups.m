function make_style_mockups()
%MAKE_STYLE_MOCKUPS  Figs. 4 and 5 of the paper: Zolotarev rationals on real intervals.
%   Paper figures (local functions, computed here with ZOLO3, no runs needed):
%     mock_intervals   Fig. 4 (ms_zolo_intervals): E = [-1, -1/8], F = [1/8, 3], types (10,10) and
%                      (15,5), with the sign-approximation error   -> figs/mock/mock_intervals.{pdf,png}
%     mock_ratio_asym  Fig. 5 (ms_ratio_asym): (a) the efficiency ratio F(1)/F(1/2) of type (p,0)
%                      against type (p,p), (b) E = [-1/4, -1/16], F = [1/8, 2], types (4,4), (8,0),
%                      (0,8), (6,2)                                -> figs/mock/mock_ratio_asym.{pdf,png}
%   The other local functions (mock_ratio, mock_asym, mock_gallery) draw variants that are not in the paper.
%
%   Call:   make_style_mockups()   draws Figs. 4 and 5 (mock_intervals, mock_ratio_asym)
%           To draw other local functions, set the environment variable MOCK_ONLY to a comma-separated
%           list of their names, e.g. in MATLAB
%               setenv('MOCK_ONLY', 'mock_ratio,mock_asym'); make_style_mockups(); setenv('MOCK_ONLY', '')
%   Writes: figs/mock/mock_<name>.{pdf,png} (mock_gallery writes gallery_*.{pdf,png} via ZOLO_GALLERY).
%   The colorbar label of Fig. 5(a), F(1)/F(1/2) in script F, is set in the LaTeX source, since MATLAB
%   cannot render \mathscr.
setup_paths();
out = fullfile(app_root(), 'figs', 'mock');
if ~exist(out, 'dir'), mkdir(out); end
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex', 'defaultColorbarTickLabelInterpreter', 'latex', ...
           'defaultAxesFontSize', 8, 'defaultTextFontSize', 8);
only = getenv('MOCK_ONLY');
if ~isempty(only), for s = strsplit(only, ','), feval(s{1}, out); end, return, end
mock_intervals(out);
mock_ratio_asym(out);
end

% ---------------------------------------------------------------- the house style
function c = col(k)
% the palette, in order: blue, red, teal, purple, amber, grey
P = [0.04 0.37 0.65;  0.82 0.29 0.36;  0.16 0.62 0.56;  0.56 0.37 0.64;  0.88 0.58 0.18;  0.42 0.46 0.49];
c = P(k, :);
end

function sty(ax)
% full box, ticks inside, light major grid, LaTeX ticks
set(ax, 'FontSize', 8, 'LineWidth', 0.6, 'Box', 'on', 'TickDir', 'in', 'TickLength', [0.012 0.012], ...
        'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.12, 'Layer', 'top', ...
        'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'TickLabelInterpreter', 'latex');
end

function lg = boxlegend(ax, h, names, loc)
% a small boxed legend inside the axes
lg = legend(ax, h, names, 'Location', loc, 'Interpreter', 'latex', 'FontSize', 7.5, 'Box', 'on', ...
            'EdgeColor', [0.75 0.75 0.75], 'Color', 'w');
lg.ItemTokenSize = [14 8];
end

function cm = divmap(n)
% deep blue - off-white - deep red (as in ZOLO_GALLERY)
a = [8 22 52; 22 63 122; 61 121 184; 156 195 228; 238 241 242; 243 201 166; 223 124 69; 168 42 38; 61 7 17]/255;
cm = interp1(linspace(0, 1, size(a, 1)), a, linspace(0, 1, n));
end

function save_mock(fig, out, name)
f = fullfile(out, ['mock_' name]);
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

% ---------------------------------------------------------------- F(1)/F(1/2) alone (Fig. 5(a) is drawn by mock_ratio_asym)
function mock_ratio(out)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.3 2.0], 'Visible', 'off');
ax = axes(fig, 'Position', [0.09 0.17 0.74 0.78]);
draw_ratio(ax, [0.86 0.17 0.025 0.78]);
f = fullfile(out, 'mock_ratio');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');     % under -batch, 'image' comes out scaled and clipped
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function draw_ratio(ax, cbpos)
% F(1)/F(1/2) for E = [-1, 1], F = [1+d, 1+d+L]: x = l(E)/l(F) = 2/L, y = log10(d/L);
% colorbar at CBPOS (normalized figure units)
x = linspace(0.004, 0.996, 400);
y = linspace(-3, 3, 361);
[Xg, Yg] = meshgrid(x, y);
L = 2./Xg;
d = L.*10.^Yg;
k2 = 2*L./((2 + d).*(d + L));
[K, ~] = ellipke(k2);
[Kp, ~] = ellipke(1 - k2);
R = log(1 + d + sqrt(d.^2 + 2*d))./(pi*Kp./(2*K));
hold(ax, 'on');
% band edges stroked in their own color, so the vector PDF shows no hairline seams between bands
contourf(ax, x, y, R, linspace(0, 2, 41), 'LineColor', 'flat', 'LineWidth', 0.5);
colormap(ax, divmap(256));
caxis(ax, [0 2]);
lo = [0.2 0.4 0.6 0.8];
hi = [1.1 1.2 1.3 1.4];
[c1, h1] = contour(ax, x, y, R, lo, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.5);
[c2, h2] = contour(ax, x, y, R, hi, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.5);
[c0, h0] = contour(ax, x, y, R, [1 1], '-', 'Color', 'k', 'LineWidth', 1.1);
clabel(c1, h1, 'FontSize', 6.5, 'Interpreter', 'latex', 'LabelSpacing', 300, 'Color', [0.15 0.15 0.15]);
clabel(c2, h2, 'FontSize', 6.5, 'Interpreter', 'latex', 'LabelSpacing', 200, 'Color', [0.15 0.15 0.15]);
clabel(c0, h0, 'FontSize', 6.5, 'Interpreter', 'latex', 'LabelSpacing', 400);
sty(ax);
set(ax, 'XLim', [0 1], 'YLim', [-3 3], 'XTick', 0:0.2:1, 'YTick', -3:3);
grid(ax, 'off');
xlabel(ax, '$\ell(E)/\ell(F)$');
ylabel(ax, '$\log_{10}\bigl(d(E,F)/\ell(F)\bigr)$');
cb = colorbar(ax, 'Position', cbpos);
cb.TickLabelInterpreter = 'latex';
cb.Ticks = 0:0.5:2;
% no colorbar label: the paper names the ratio \mathscr{F}(1)/\mathscr{F}(1/2), which MATLAB's LaTeX
% interpreter cannot render, so the paper sets the label beside the image
end

% ---------------------------------------------------------------- Zolotarev rationals on two intervals
function [r, g, sgm] = zolotype(E, F, nz, np)
% the type-(nz, np) rational, scaled so that max_E|r| = sqrt(sigma) = 1/min_F|r|, and its sign approximant g
if nz <= np
    o = zolo3(E, F, np, 'delta', np - nz, 'lawson', 400, 'damping', 0.9);
    sgm = o.sigma;
    r = @(z) o.r3(z)/sqrt(sgm);
    g = @(z) o.r4(z);
else
    o = zolo3(F, E, nz, 'delta', nz - np, 'lawson', 400, 'damping', 0.9);
    sgm = o.sigma;
    r = @(z) sqrt(sgm)./o.r3(z);
    g = @(z) -o.r4(z);
end
end

function h = onsets(ax, xs, y, e, f, c, lw)
% |r| solid on E and F, dashed between and outside them; returns the solid line for the legend
in = (xs >= e(1) & xs <= e(2)) | (xs >= f(1) & xs <= f(2));
ys = y;
ys(~in) = NaN;
yd = y;
yd(in & ~ismember(xs, [e(:); f(:)])) = NaN;
plot(ax, xs, yd, '--', 'Color', c, 'LineWidth', 0.7);
h = plot(ax, xs, ys, '-', 'Color', c, 'LineWidth', lw);
end

function [E, F] = intervals(e, f, M)
% M Chebyshev-clustered points on each of the intervals e and f
c = (1 - cos(pi*(0:M-1)'/(M-1)))/2;
E = e(1) + (e(2) - e(1))*c;
F = f(1) + (f(2) - f(1))*c;
end

function mock_intervals(out)
% Fig. 4: E = [-1, -1/8], F = [1/8, 3], types (10,10) and (15,5)
[E, F] = intervals([-1 -1/8], [1/8 3], 900);
ty = [10 10; 15 5];
xs = unique([linspace(-1, 3, 4000)'; -1; -1/8; 1/8; 3]);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.4 2.2], 'Visible', 'off');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl);
hold(a1, 'on');
a2 = nexttile(tl);
hold(a2, 'on');
h = gobjects(2, 1);
nm = cell(2, 1);
for k = 1:2
    [r, g, sg] = zolotype(E, F, ty(k, 1), ty(k, 2));
    h(k) = onsets(a1, xs, abs(r(xs)), [-1 -1/8], [1/8 3], col(k), 1.0 + 0.3*(k == 1));
    nm{k} = sprintf('type $(%d,%d)$, $\\xi=%s$', ty(k, 1), ty(k, 2), texnum(sg));
    plot(a2, E, real(g(E)) + 1, '-', 'Color', col(k), 'LineWidth', 1.0 + 0.3*(k == 1));
    plot(a2, F, real(g(F)) - 1, '-', 'Color', col(k), 'LineWidth', 1.0 + 0.3*(k == 1));
    fprintf('intervals: type (%d,%d), xi = %.2e\n', ty(k, 1), ty(k, 2), sg);
end
sty(a1);
set(a1, 'YScale', 'log', 'XLim', [-1 3], 'YLim', [1e-8 1e8], 'YTick', 10.^(-8:4:8), 'XTickLabel', []);
ylabel(a1, '$|r(x)|$');
boxlegend(a1, h, nm, 'east');
text(a1, -0.56, 1e-3, '$E$', 'HorizontalAlignment', 'center', 'FontSize', 9);
text(a1, 0.6, 1e2, '$F$', 'HorizontalAlignment', 'center', 'FontSize', 9);
sty(a2);
set(a2, 'XLim', [-1 3], 'YLim', [-1.25e-5 1.25e-5]);
ylabel(a2, '$g(x)-\mathrm{sgn}(x)$');
xlabel(a2, '$x$');
grid(a1, 'off');
grid(a2, 'off');
save_mock(fig, out, 'intervals');
end

function mock_asym(out)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.6 1.65], 'Visible', 'off');
ax = axes(fig, 'Position', [0.1 0.2 0.86 0.74]);
draw_asym(ax);
save_mock(fig, out, 'asym');
end

function draw_asym(ax)
% E = [-1/4, -1/16], F = [1/8, 2], types (4,4), (8,0), (0,8), (6,2)
[E, F] = intervals([-1/4 -1/16], [1/8 2], 900);
ty = [4 4; 8 0; 0 8; 6 2];
xs = unique([linspace(-0.3, 2, 4000)'; -1/4; -1/16; 1/8; 2]);
hold(ax, 'on');
h = gobjects(4, 1);
nm = cell(4, 1);
for k = 1:4
    try
        [r, ~, sg] = zolotype(E, F, ty(k, 1), ty(k, 2));
        h(k) = onsets(ax, xs, abs(r(xs)), [-1/4 -1/16], [1/8 2], col(k), 1.0);
        nm{k} = sprintf('$(%d,%d)$: $\\xi=%s$', ty(k, 1), ty(k, 2), texnum(sg));
        fprintf('asym: type (%d,%d), xi = %.2e\n', ty(k, 1), ty(k, 2), sg);
    catch e
        fprintf(2, 'asym: type (%d,%d) failed: %s\n', ty(k, 1), ty(k, 2), e.message);
        h(k) = plot(ax, NaN, NaN);
        nm{k} = '';
    end
end
ok = ~cellfun(@isempty, nm);
sty(ax);
set(ax, 'YScale', 'log', 'XLim', [-0.3 2], 'YLim', [1e-5 1e5], 'YTick', 10.^(-4:2:4));
grid(ax, 'off');
text(ax, -0.156, 0.6, '$E$', 'HorizontalAlignment', 'center', 'FontSize', 9);
text(ax, 0.42, 1e-2, '$F$', 'HorizontalAlignment', 'center', 'FontSize', 9);       % clear of the legend
xlabel(ax, '$x$');
ylabel(ax, '$|r(x)|$');
boxlegend(ax, h(ok), nm(ok), 'southeast');
end

function mock_gallery(out)
% both gallery sets with the names of the condensers, and every condenser alone (for a three-panel figure)
old = cd(out);
c = onCleanup(@() cd(old));
set(groot, 'defaultTextInterpreter', 'tex', 'defaultLegendInterpreter', 'tex', 'defaultAxesTickLabelInterpreter', 'tex');
for k = 1:2
    zolo_gallery(k, 'save', true, 'markersize', 4, 'linewidth', 1.4, 'title', true, 'grid', 600, 'name', sprintf('gallery_%d_named', k), 'raster', true);
    close all
end
zolo_gallery('all', 'save', true, 'markersize', 3.5, 'linewidth', 1.2, 'title', true, 'grid', 500, 'size', [1.72 1.72], 'fontsize', 7, 'raster', true);
close all
set(groot, 'defaultTextInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex');
end

function s = texnum(v)
% v as c\cdot10^{e} with one decimal
e = floor(log10(v));
s = sprintf('%.1f\\cdot10^{%d}', v/10^e, e);
end

% ---------------------------------------------------------------- Fig. 5: (a) above (b)
function mock_ratio_asym(out)
% the two panels in one figure, plot areas aligned (same left edge and width); the colorbar label
% F(1)/F(1/2) is set by LaTeX beside the image, since MATLAB cannot render \mathscr
W = 4.05;  H = 3.40;  L = 0.6;  w = 2.95;  hb = 1.35;  ya = 0.36 + hb + 0.47;   % inches; (b) is tall enough for its legend
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 W H], 'Visible', 'off');
axb = axes(fig, 'Units', 'inches', 'Position', [L 0.36 w hb]);
axa = axes(fig, 'Units', 'inches', 'Position', [L ya w 1.15]);
draw_ratio(axa, [(L + w + 0.08)/W, ya/H, 0.11/W, 1.15/H]);
draw_asym(axb);
text(axa, -0.19, 1, '(a)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
text(axb, -0.19, 1, '(b)', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
save_mock(fig, out, 'ratio_asym');
end
