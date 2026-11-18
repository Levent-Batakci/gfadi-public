function out = zolo_gallery(which, varargin)
%ZOLO_GALLERY  Lopsided Zolotarev rationals on ten condensers, drawn as heatmaps.
%   OUT = ZOLO_GALLERY(WHICH,...) computes each panel with ZOLO3 and draws
%   u(z) = 1 - 2*log|r(z)|/log(sigma)  (u <= -1 on E, u >= 1 on F, dashed u = 0).
%   WHICH: 1 or 2 (gallery of panels a-e or f-j), a panel name or letter, a cell of
%   them (one figure each), or 'all'.  Panels: eclipse, galaxy, yinyang, gears, cage,
%   sierpinski, mandelbrot, combs, atom, mandala (letters a-j).
%   OUT: struct array per panel (name, label, type, nlaw, E, F, oE, oF, sigma, zer, pol, r, lw).
%   Type (m,n) = (#zeros, #poles); m > n is ZOLO3(F,E,...) with r -> sigma/r.
%   Options (defaults): 'grid' 700, 'save' false (writes .png/.pdf to the current folder),
%     'draw' true, 'lawson' [] (per-panel), 'type' [] (override), 'row' false (panels in one
%     row), 'markersize' 15, 'linewidth' 6, 'halo' 0.8, 'polewidth' 1.6, 'zerowidth' 0.9,
%     'pad' 1.3, 'raster' false; single panels: 'size' [3.4 3.4], 'fontsize' 9,
%     'colorbar' false, 'title' false, 'sigmabox' false, 'name' '' (file stem).
%   Fig. 2 of the paper (gears, yin-yang, combs in one row) is drawn by MAKE_ZOLO_FIG.
%   Example:  zolo_gallery('mandala', 'save', true);

if nargin < 1, which = 1; end
opt = struct('grid', 700, 'save', false, 'draw', true, 'lawson', [], ...
             'markersize', 15, 'linewidth', 6, 'halo', 0.8, ...
             'size', [3.4 3.4], 'fontsize', 9, 'colorbar', false, ...
             'title', false, 'sigmabox', false, 'name', '', 'raster', false, 'row', false, 'type', [], ...
             'polewidth', 1.6, 'zerowidth', 0.9, 'pad', 1.3);
for k = 1:2:numel(varargin), opt.(lower(varargin{k})) = varargin{k+1}; end

%            name          label                    builder          type     lawson  linewidth
PA = {'eclipse',    '(a) eclipse',        @g_eclipse,    [ 8 16],  400, 0.8
      'galaxy',     '(b) galaxy',         @g_galaxy,     [10 20],  400, 1.3
      'yinyang',    '(c) yin-yang',       @g_yinyang,    [24 16],  400, 0.8
      'gears',      '(d) gears',          @g_gears,      [15 10],  400, 0.8
      'cage',       '(e) snowflake cage', @g_cage,       [18  6],  400, 0.8
      'sierpinski', '(f) Sierpinski',     @g_sierpinski, [ 6 18],  400, 0.8
      'mandelbrot', '(g) Mandelbrot',     @g_mandelbrot, [16  8],  400, 0.8
      'combs',      '(h) combs',          @g_combs,      [16 32], 1500, 0.8
      'atom',       '(i) atom',           @g_atom,       [ 8 24],  400, 1.0
      'mandala',    '(j) mandala',        @g_mandala,    [18  6], 1500, 0.8};

[idx, single, stem] = resolve(which, PA);   % which panels, and in which layout
P = PA(idx, :);

fprintf('\n%-12s %10s %12s %8s %8s %8s\n', 'panel', 'type', 'sigma', '#zeros', '#poles', 'time');
out = struct([]);
for k = 1:size(P, 1)
    tic;
    [E, F, oE, oF] = P{k,3}();
    ty = P{k,4};  if ~isempty(opt.type), ty = opt.type; end   % 'type' overrides the table (one panel at a time)
    nlaw = P{k,5};  if ~isempty(opt.lawson), nlaw = opt.lawson; end
    [r, zer, pol, sigma] = zolotype(E, F, ty(1), ty(2), nlaw);
    sc  = max(abs([E; F]));
    zer = zer(abs(zer) < 10*sc);          % drop the numerically infinite eigenvalues
    pol = pol(abs(pol) < 10*sc);
    out(k).name = P{k,1};  out(k).label = P{k,2};  out(k).type = ty;  out(k).nlaw = nlaw;
    out(k).E = E;  out(k).F = F;  out(k).oE = oE;  out(k).oF = oF;
    out(k).sigma = sigma;  out(k).zer = zer;  out(k).pol = pol;  out(k).r = r;
    out(k).lw = P{k,6};
    fprintf('%-12s %4d,%-5d %12.3e %8d %8d %7.1fs\n', P{k,1}, ty(1), ty(2), sigma, ...
            numel(zer), numel(pol), toc);
end

if opt.draw
    if opt.row                                     % the chosen panels side by side, relettered, key below
        nm = opt.name;  if isempty(nm), nm = 'zolo_row'; end
        draw_row(out, nm, opt);
    elseif single
        for k = 1:numel(out)                       % one figure per panel
            nm = opt.name;  if isempty(nm), nm = out(k).name; end
            draw_single(out(k), nm, opt);
        end
    else
        nm = opt.name;  if isempty(nm), nm = stem; end
        draw_gallery(out, nm, opt);
    end
end
end

% ============================================================ panel selection

function [idx, single, stem] = resolve(which, PA)
%RESOLVE  Turn the first argument into row indices into PA, plus the layout.
%   Numeric 1 or 2 keeps the gallery layout; anything naming panels draws them
%   one at a time.
names = PA(:,1);
if isnumeric(which) && isscalar(which) && any(which == [1 2])
    idx = (1:5) + 5*(which - 1);  single = false;
    stem = sprintf('zolo_gallery_%d', which);
    return
end
single = true;  stem = '';
if ischar(which) && strcmpi(which, 'all')
    idx = 1:numel(names);  return
end
if ~iscell(which), which = {which}; end
idx = zeros(1, numel(which));
for k = 1:numel(which)
    w = which{k};
    if isnumeric(w)
        j = w;                                      % panel number 1..10
    elseif numel(w) == 1 && any(lower(w) >= 'a' & lower(w) <= 'j')
        j = double(lower(w)) - double('a') + 1;     % panel letter
    else
        j = find(strcmpi(names, w), 1);             % panel name
    end
    if isempty(j) || j < 1 || j > numel(names)
        error('zolo_gallery:panel', ...
              'Unknown panel "%s".  Use 1, 2, ''all'', a letter a-j, or one of: %s.', ...
              num2str(w), strjoin(names', ', '));
    end
    idx(k) = j;
end
end

% =============================================================== solver wrapper

function [r, zer, pol, sigma] = zolotype(E, F, nz, np, nlaw)
if nz <= np
    o = zolo3(E, F, np, 'delta', np - nz, 'lawson', nlaw, 'damping', 0.9);
    r = o.r3;  zer = o.zer;  pol = o.pol;  sigma = o.sigma;
else
    o = zolo3(F, E, nz, 'delta', nz - np, 'lawson', nlaw, 'damping', 0.9);
    sigma = o.sigma;
    r3 = o.r3;
    r = @(z) sigma./r3(z);                % max_E|r| = sigma, min_F|r| = 1
    zer = o.pol;  pol = o.zer;
end
end

% ==================================================================== drawing

function draw_gallery(out, stem, opt)
UMAX = 1.8;  LEV = linspace(-UMAX, UMAX, 37);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7 4.9]);
set(fig, 'Colormap', zolocmap(256));
W = 0.29;  H = 0.415;  gx = 0.035;  gy = 0.07;  x0 = 0.03;  y0 = 0.035;
for k = 1:numel(out)
    rr = floor((k-1)/3);  cc = mod(k-1, 3);
    ax = axes('Position', [x0 + cc*(W+gx), y0 + (1-rr)*(H+gy), W, H]);
    panel(ax, out(k), opt, LEV, UMAX);
end

% key: colour bar and marker legend in the sixth cell
axk = axes('Position', [x0 + 2*(W+gx), y0, W, H], 'Visible', 'off');
hold(axk, 'on');  caxis(axk, [-UMAX UMAX]);
cb = colorbar(axk);
set(cb, 'Position', [x0 + 2*(W+gx) + 0.03, y0 + 0.14, 0.022, H - 0.17]);
try, set(cb, 'Ticks', -1.5:0.5:1.5); catch, set(cb, 'YTick', -1.5:0.5:1.5); end
set(cb, 'TickLabelInterpreter', 'latex', 'TickLabels', mathlabels(-1.5:0.5:1.5));
ylabel(cb, '$1-2\log|r(z)|/\log\xi$', 'FontSize', 8, 'Interpreter', 'latex');
ms = opt.markersize;
h(1) = plot(axk, NaN, NaN, 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
h(2) = plot(axk, NaN, NaN, 'x', 'MarkerSize', 1.4*ms, 'Color', 'k', 'LineWidth', 1.6);
h(3) = plot(axk, NaN, NaN, '-', 'Color', 'k', 'LineWidth', 0.8*opt.linewidth);
h(4) = plot(axk, NaN, NaN, '--', 'Color', 'k', 'LineWidth', 0.6);
lg = legend(h, {'zeros of $r$', 'poles of $r$', '$\partial E$, $\partial F$', '$|r|=\sqrt{\xi}$'}, 'Interpreter', 'latex');
set(lg, 'Box', 'off', 'FontSize', 8, ...
        'Position', [x0 + 2*(W+gx) + 0.135, y0 + 0.17, 0.15, 0.2]);

if opt.save, writefig(fig, stem, [7 4.9], opt.raster); end
end

% --------------------------------------------------------- one panel, one file

function draw_single(o, stem, opt)
%DRAW_SINGLE  One panel in its own square figure, sized for a paper column.
UMAX = 1.8;  LEV = linspace(-UMAX, UMAX, 37);
fs = opt.fontsize;
wfig = opt.size(1);  hfig = opt.size(2);
cbw = 0;  if opt.colorbar, cbw = 0.75 + 0.022*fs; end   % room for bar and labels
wfig = wfig + cbw;

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 wfig hfig]);
set(fig, 'Colormap', zolocmap(256));
mL = 0.04;  mR = 0.04 + cbw/wfig;  mB = 0.04;
mT = 0.04;  if opt.title, mT = 0.055 + 1.7*fs/72/hfig; end
ax = axes('Position', [mL, mB, 1 - mL - mR, 1 - mB - mT]);
panel(ax, o, opt, LEV, UMAX, fs);

if opt.colorbar
    % peer the bar to a hidden axes: colorbar(ax) would resize the panel itself
    axk = axes('Position', [1 - mR + 0.02, mB, 0.5/wfig, 1 - mB - mT], 'Visible', 'off');
    caxis(axk, [-UMAX UMAX]);
    cb = colorbar(axk);
    set(cb, 'Position', [1 - mR + 0.02, mB + 0.05, 0.4/wfig, 1 - mB - mT - 0.1], ...
            'FontSize', fs - 1);
    try, set(cb, 'Ticks', -1.5:0.5:1.5); catch, set(cb, 'YTick', -1.5:0.5:1.5); end
    set(cb, 'TickLabelInterpreter', 'latex', 'TickLabels', mathlabels(-1.5:0.5:1.5));
    ylabel(cb, '$1-2\log|r(z)|/\log\xi$', 'FontSize', fs, 'Interpreter', 'latex');
end
if opt.save, writefig(fig, ['zolo_gallery_' stem], [wfig hfig], opt.raster); end
end

% ------------------------------------------------- a few panels in one row, key below

function draw_row(out, stem, opt)
%DRAW_ROW  The chosen panels side by side, lettered (a), (b), ..., with the colour bar and the marker key in a
%   strip below; sized for 0.85 of the paper's text width.
UMAX = 1.8;  LEV = linspace(-UMAX, UMAX, 37);  fs = opt.fontsize;  n = numel(out);
wfig = 4.36;  hfig = 1.9;   % printed size, 0.85 textwidth
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 wfig hfig]);
set(fig, 'Colormap', zolocmap(256));
gap = 0.035;  W = (1 - (n + 1)*gap)/n;  H = W*wfig/hfig;  y0 = 1 - H - 0.085;
for k = 1:n
    ax = axes('Position', [gap + (k-1)*(W + gap), y0, W, H]);
    o = out(k);  nm = regexprep(o.label, '^\(\w\)\s*', '');
    o.label = sprintf('(%c) %s', 'a' + k - 1, nm);
    panel(ax, o, setfield(opt, 'title', true), LEV, UMAX, fs);           %#ok<SFLD>
end
% the key: a colour strip under panel (a) with its label above it, the markers under the other panels
yk = 0.075;  hk = 0.045;
axc = axes('Position', [gap, yk, W, hk]);
imagesc(axc, linspace(-UMAX, UMAX, 256), [0 1], linspace(-UMAX, UMAX, 256));  colormap(axc, zolocmap(256));
set(axc, 'YTick', [], 'XTick', -1.5:0.75:1.5, 'XTickLabel', mathlabels(-1.5:0.75:1.5), 'XTickLabelRotation', 0, 'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.6, 'FontSize', fs - 1, ...
         'TickLabelInterpreter', 'latex', 'Layer', 'top');
text(axc, 0.5, 1.12, '$1-2\log|r(z)|/\log\xi$', 'Units', 'normalized', 'Interpreter', 'latex', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs - 1);
axk = axes('Position', [gap + W + gap, 0.02, 1 - 3*gap - W, 0.18], 'Visible', 'off');  hold(axk, 'on');
ms = opt.markersize;
h(1) = plot(axk, NaN, NaN, 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', opt.zerowidth);
h(2) = plot(axk, NaN, NaN, 'x', 'MarkerSize', 1.4*ms, 'Color', 'k', 'LineWidth', opt.polewidth);
h(3) = plot(axk, NaN, NaN, '-', 'Color', 'k', 'LineWidth', 0.8*opt.linewidth);
h(4) = plot(axk, NaN, NaN, '--', 'Color', 'k', 'LineWidth', 0.6);
lg = legend(axk, h, {'zeros of $r$', 'poles of $r$', '$\partial E$, $\partial F$', '$|r|=\sqrt{\xi}$'}, ...
            'Interpreter', 'latex', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', fs - 1);
lg.ItemTokenSize = [14 8];  set(lg, 'Units', 'normalized');
p = lg.Position;  set(lg, 'Position', [gap + W + gap + (1 - 3*gap - W - p(3))/2, yk + hk/2 - p(4)/2 + 0.02, p(3), p(4)]);
if opt.save, writefig(fig, stem, [wfig hfig], opt.raster); end
end

% ------------------------------------------------------------------ file output

function writefig(fig, stem, sz, raster)
if nargin < 4, raster = false; end
set(fig, 'PaperUnits', 'inches', 'PaperSize', sz, 'PaperPosition', [0 0 sz]);
print(fig, [stem '.png'], '-dpng', '-r300');
try
    if exist('exportgraphics', 'file') && raster      % filled contours as an image: no seams in any viewer
        exportgraphics(fig, [stem '.pdf'], 'ContentType', 'image', 'Resolution', 600, 'BackgroundColor', 'white');
    elseif exist('exportgraphics', 'file')    % MATLAB R2020a+
        exportgraphics(fig, [stem '.pdf'], 'ContentType', 'vector');
    else                                  % Octave: needs Ghostscript
        print(fig, [stem '.pdf'], '-dpdf');
    end
catch err
    warning('zolo_gallery:pdf', 'PDF export failed (%s); PNG was written.', err.message);
end
end

function panel(ax, o, opt, LEV, UMAX, fs)
if nargin < 6, fs = 9; end               % gallery sizes; draw_single passes its own
ng = opt.grid;
allz = [o.E; o.F];
cx = (min(real(allz)) + max(real(allz)))/2;
cy = (min(imag(allz)) + max(imag(allz)))/2;
w  = opt.pad*max(max(real(allz)) - min(real(allz)), max(imag(allz)) - min(imag(allz)))/2;
xg = linspace(cx - w, cx + w, ng);
yg = linspace(cy - w, cy + w, ng);
[X, Y] = meshgrid(xg, yg);
Z  = X + 1i*Y;
L  = zeros(size(Z));
ch = 40000;
for i0 = 1:ch:numel(Z)                    % chunked: bounded memory in baryratio
    i1 = min(i0 + ch - 1, numel(Z));
    L(i0:i1) = log10(abs(o.r(Z(i0:i1))));
end
U = 1 - 2*L/log10(o.sigma);
U(isnan(U)) = UMAX;
U = min(max(U, -UMAX), UMAX);

hold(ax, 'on');
contourf(ax, xg, yg, U, LEV, 'LineColor', 'flat', 'LineWidth', 0.5);   % band edges in their own color: no seams in a vector PDF
[~, hc] = contour(ax, xg, yg, U, LEV(1:3:end));
set(hc, 'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.2);
[~, h0] = contour(ax, xg, yg, U, [0 0]);
set(h0, 'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
caxis(ax, [LEV(1) LEV(end)]);

lw = opt.linewidth*o.lw;  halo = lw + opt.halo;
ms = opt.markersize;
cur = [o.oE(:); o.oF(:)];
for j = 1:numel(cur)                      % thin dense curves: renderers dash them
    c = cur{j};  st = max(1, floor(numel(c)/600));
    cur{j} = c([1:st:end-1, end]);
end
if opt.halo > 0
    for j = 1:numel(cur)
        plot(ax, real(cur{j}), imag(cur{j}), '-', 'Color', 'w', 'LineWidth', halo, ...
             'LineJoin', 'round');
    end
end
for j = 1:numel(cur)
    plot(ax, real(cur{j}), imag(cur{j}), '-', 'Color', 'k', 'LineWidth', lw, ...
         'LineJoin', 'round');
end
plot(ax, real(o.zer), imag(o.zer), 'o', 'MarkerSize', ms, 'MarkerFaceColor', 'w', ...
     'MarkerEdgeColor', 'k', 'LineWidth', opt.zerowidth);
if opt.halo > 0
    plot(ax, real(o.pol), imag(o.pol), 'x', 'MarkerSize', 1.4*ms, 'Color', 'w', ...
         'LineWidth', opt.polewidth + opt.halo);
end
plot(ax, real(o.pol), imag(o.pol), 'x', 'MarkerSize', 1.4*ms, 'Color', 'k', ...
     'LineWidth', opt.polewidth);

axis(ax, 'equal');  xlim(ax, xg([1 end]));  ylim(ax, yg([1 end]));
set(ax, 'XTick', [], 'YTick', [], 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.6);
if ~isequal(opt.title, false)
    lb = o.label;  if ischar(opt.title) && strcmpi(opt.title, 'letter'), lb = strtok(lb); end
    text(ax, 0, 1.015, lb, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
         'FontSize', fs, 'Interpreter', 'latex');
    text(ax, 1, 1.015, sprintf('type $(%d,%d)$', o.type(1), o.type(2)), 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs - 1, 'Interpreter', 'latex');
end
if opt.sigmabox                           % off by default: sigma goes in the caption
    e = floor(log10(o.sigma));
    text(ax, 0.042, 0.055, sprintf('$\\xi = %.1f\\cdot10^{%d}$', o.sigma/10^e, e), 'Interpreter', 'latex', ...
         'Units', 'normalized', 'BackgroundColor', 'w', 'Margin', 1.5, 'FontSize', fs - 1.5);
end
end

function c = mathlabels(t)
%MATHLABELS  Tick labels in math mode, so that negative values get a true minus sign under the LaTeX interpreter.
c = arrayfun(@(x) sprintf('$%g$', x), t, 'UniformOutput', false);
end

function cm = zolocmap(n)
%ZOLOCMAP  Diverging map: deep blue (E side) through off-white to deep red (F side).
anchors = [  8  22  52;  22  63 122;  61 121 184; 156 195 228; 238 241 242
           243 201 166; 223 124  69; 168  42  38;  61   7  17]/255;
x  = linspace(0, 1, size(anchors, 1));
cm = interp1(x, anchors, linspace(0, 1, n));
end

% ================================================================= geometries

function s = cheb01(N)
s = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
end

function z = arc(c, R, a0, a1, N)
z = c + R*exp(1i*(a0 + (a1 - a0)*cheb01(N)));
end

function z = circ(c, R, N)
z = c + R*exp(2i*pi*(0:N-1)'/N);
end

function z = closec(z)
z = [z(:); z(1)];
end

function z = polysample(P, density, nmin)
%POLYSAMPLE  Sample a closed polygon, Chebyshev-clustered towards each vertex.
Q = [P(:); P(1)];  z = [];
for j = 1:numel(Q)-1
    a = Q(j);  b = Q(j+1);  L = abs(b - a);
    if L < 1e-12, continue; end
    k = max(nmin, ceil(density*L));
    s = cheb01(k + 1);
    z = [z; a + (b - a)*s(1:end-1)];     %#ok<AGROW>
end
end

% (a) crescent E wrapped around a disk F
function [E, F, oE, oF] = g_eclipse()
c2 = 0.22;  r2 = 0.88;
xi = (1 - r2^2 + c2^2)/(2*c2);  yi = sqrt(1 - xi^2);
a1 = atan2(yi, xi);  b1 = atan2(yi, xi - c2);
outer = arc(0, 1, a1, 2*pi - a1, 700);
inner = arc(c2, r2, 2*pi - b1, b1, 700);
E = [outer; inner(2:end-1)];
F = circ(0.36, 0.56, 700);
oE = {closec([outer; inner])};  oF = {closec(F)};
end

% (b) two interleaved logarithmic spiral arcs; E has the extra outer turn
function [E, F, oE, oF] = g_galaxy()
b  = 0.21;
tE = 0.25*pi + 2.75*pi*cheb01(1400);
tF = 0.25*pi + 1.85*pi*cheb01(1000);
E  = exp((b + 1i)*tE);  F = -exp((b + 1i)*tF);
s  = max(abs(E));  E = E/s;  F = F/s;
oE = {E};  oF = {F};
end

% (c) yin-yang, eroded by ep on either side of every interface, eyes included
function [E, F, oE, oF] = g_yinyang()
ep = 0.11;  rho = 0.2;  N = 500;
R1 = 1 - ep;
yt = (0.5 + ep)^2 - R1^2 - 0.25;  xt = -sqrt(R1^2 - yt^2);   % tail tip
at = atan2(yt, xt) + 2*pi;  bt = atan2(yt + 0.5, xt) + 2*pi;
o1 = arc(0, R1, pi/2, at, N);                 % outer rim, head -> tail
o2 = arc(-0.5i, 0.5 + ep, bt, pi/2, N);       % S-curve, lower half
o3 = arc(0.5i, 0.5 - ep, -pi/2, pi/2, N);     % S-curve, upper half
loop = [o1; o2(2:end); o3(2:end-1)];
hole = circ(0.5i, rho + ep, 300);
eyeE = circ(-0.5i, rho - ep, 200);
E  = [loop; hole; eyeE];  F = -E;
oE = {closec(loop), closec(hole), closec(eyeE)};
oF = cellfun(@(c) -c, oE, 'UniformOutput', false);
end

% (d) small gear E meshing with a large gear F
function z = gear(c, R, h, k, phase, N)
duty = 0.55;  beta = 5;
th = 2*pi*(0:N-1)'/N;
r  = R + h*tanh(beta*(cos(k*(th - phase)) - duty))/tanh(beta*(1 - duty));
z  = c + r.*exp(1i*th);
end

function [E, F, oE, oF] = g_gears()
R1 = 0.5;  k1 = 9;  R2 = 1.0;  k2 = 18;  h = 0.075;
D  = R1 + R2 + 0.075;
E  = gear(0, R1, h, k1, 0, 1400);
F  = gear(D, R2, h, k2, pi - pi/k2, 2800);
E  = E - 0.42*D;  F = F - 0.42*D;
oE = {closec(E)};  oF = {closec(F)};
end

% (e) Koch curve F enclosing a small off-centre Koch region E
function [dense, P] = koch(level, R, rot, c, perseg)
P = R*exp(1i*(pi/2 + rot + 2*pi*(0:3)'/3));
for l = 1:level
    a = P(1:end-1);  b = P(2:end);  d = (b - a)/3;
    p1 = a + d;  p3 = a + 2*d;  p2 = p1 + d*exp(-1i*pi/3);
    nw = [a.'; p1.'; p2.'; p3.'];
    P = [nw(:); P(end)];
end
s = (0:perseg-1)/perseg;
D = P(1:end-1) + (P(2:end) - P(1:end-1))*s;
dense = c + reshape(D.', [], 1);
P = c + P;
end

function [E, F, oE, oF] = g_cage()
[F, Fo] = koch(4, 1.0, 0, 0, 3);
[E, Eo] = koch(3, 0.27, pi/6, 0.13 + 0.08i, 4);
oE = {Eo};  oF = {Fo};
end

% (f) level-4 Sierpinski gasket E (shrunk triangles) around the central hole F
function [E, F, oE, oF] = g_sierpinski()
level = 4;  shrink = 0.80;  fshrink = 0.50;
A = 0;  B = 1;  C = 0.5 + 1i*sqrt(3)/2;
tris = [A B C];  holes = zeros(0, 3);
for gen = 0:level-1
    nw = zeros(3*size(tris, 1), 3);
    for j = 1:size(tris, 1)
        a = tris(j,1);  b = tris(j,2);  c = tris(j,3);
        ab = (a + b)/2;  bc = (b + c)/2;  ca = (c + a)/2;
        if gen < 1, holes = [holes; ab bc ca]; end     %#ok<AGROW>
        nw(3*j-2:3*j, :) = [a ab ca; ab b bc; ca bc c];
    end
    tris = nw;
end
g0 = (A + B + C)/3;
shr = @(t, s) mean(t) + s*(t - mean(t));
E = [];  oE = {};
for j = 1:size(tris, 1)
    t = 2*(shr(tris(j,:), shrink) - g0);
    E = [E; polysample(t, 90, 4)];  oE{end+1} = closec(t.');        %#ok<AGROW>
end
F = [];  oF = {};
for j = 1:size(holes, 1)
    t = 2*(shr(holes(j,:), fshrink) - g0);
    F = [F; polysample(t, 160, 10)];  oF{end+1} = closec(t.');      %#ok<AGROW>
end
end

% (g) F = lemniscate {|p_8(c)| = 2} of the Mandelbrot iteration, E = a disk
function z = mandel_lemniscate(k, R, npts)
xs = linspace(-2.25, 0.75, 1800);  ys = linspace(-1.35, 1.35, 1620);
[X, Y] = meshgrid(xs, ys);  Cg = X + 1i*Y;  Z = zeros(size(Cg));
for j = 1:k
    Z = Z.^2 + Cg;
    Z(abs(Z) > 1e6) = 1e6;
end
Cm = contourc(xs, ys, log(max(abs(Z), 1e-300)), log(R)*[1 1]);
best = [];  i = 1;
while i < size(Cm, 2)                          % keep the longest contour piece
    np = Cm(2, i);
    seg = Cm(1, i+1:i+np) + 1i*Cm(2, i+1:i+np);
    if numel(seg) > numel(best), best = seg; end
    i = i + np + 1;
end
z = best(:);
z = z([true; abs(diff(z)) > 0]);
s = [0; cumsum(abs(diff(z)))];
si = linspace(0, s(end), npts + 1)';  si = si(1:end-1);
z = interp1(s, real(z), si) + 1i*interp1(s, imag(z), si);
end

function [E, F, oE, oF] = g_mandelbrot()
F = mandel_lemniscate(8, 2, 3000) + 0.45;
E = circ(1.35, 0.22, 600);
oE = {closec(E)};  oF = {closec(F)};
end

% (h) interdigitated combs: E has three fingers and the long bar, F two fingers
function P = comb(xl, xr, yb, yt, centres, w, ytip)
P = [xl + 1i*yb; xr + 1i*yb; xr + 1i*yt];
for c = sort(centres, 'descend')
    P = [P; c + w + 1i*yt; c + w + 1i*ytip; c - w + 1i*ytip; c - w + 1i*yt];  %#ok<AGROW>
end
P = [P; xl + 1i*yt];
P = P([abs(diff(P)) > 1e-12; true]);
end

function [E, F, oE, oF] = g_combs()
w  = 0.1;
Pe = comb(-1.3, 1.3, -1.05, -0.85, [-1.2 0 1.2], w, 0.45);
Pf = conj(comb(-0.7, 0.7, -1.05, -0.85, [-0.6 0.6], w, 0.45));   % fingers down
E  = polysample(Pe, 110, 8);  F = polysample(Pf, 110, 8);
% the condenser is symmetric about the imaginary axis; make the samples exactly so, so that the computed r is too
mir = @(z) unique([z; -conj(z)]);
if getenv('ZOLO_COMBS_SYM') ~= "0", E = mir(E);  F = mir(F); end
oE = {closec(Pe)};  oF = {closec(Pf)};
end

% (i) three elliptic orbits plus electrons (E) around a four-nucleon nucleus (F)
function [E, F, oE, oF] = g_atom()
t   = 2*pi*(0:1499)'/1500;
ell = cos(t) + 0.3i*sin(t);
orb = {ell, ell*exp(1i*pi/3), ell*exp(2i*pi/3)};
fr  = [0.10 0.45 0.78];
el  = cell(1, 3);
for k = 1:3
    o = orb{k};  el{k} = circ(o(floor(1500*fr(k)) + 1), 0.055, 120);
end
nuc = cell(1, 4);
for k = 1:4
    nuc{k} = circ(0.11*exp(1i*(pi/4 + pi/2*(k-1))), 0.075, 240);
end
E  = vertcat(orb{:}, el{:});  F = vertcat(nuc{:});
oE = cellfun(@closec, [orb el], 'UniformOutput', false);
oF = cellfun(@closec, nuc, 'UniformOutput', false);
end

% (j) ring of beads E between a central disk and an outer ring of beads (both F)
function [E, F, oE, oF] = g_mandala()
Em = cell(1, 6);  Fo = cell(1, 6);
for k = 1:6
    Em{k} = circ(0.55*exp(2i*pi*(k-1)/6), 0.14, 220);
    Fo{k} = circ(exp(2i*pi*(k-0.5)/6), 0.16, 240);
end
Fc = {circ(0, 0.2, 300)};
E  = vertcat(Em{:});  F = vertcat(Fc{:}, Fo{:});
oE = cellfun(@closec, Em, 'UniformOutput', false);
oF = cellfun(@closec, [Fc Fo], 'UniformOutput', false);
end
