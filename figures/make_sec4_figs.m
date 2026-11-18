function make_sec4_figs(which)
%MAKE_SEC4_FIGS  Fig. 9 of the paper (paper_switch, Example 3), in the style of MAKE_PAPER_FIGS.
%
%   Call:   make_sec4_figs('switch')   Fig. 9: (a) the solution averaged over the stationary distribution
%                                      of the chain on the slice through the target, (b) run time and
%                                      (c) peak memory against 120 N^3
%           make_sec4_figs()           also tries 'galaxy', 'sponge' and 'big', which draw figures that
%                                      are not in the paper and whose inputs are not produced in this
%                                      repository; their errors are printed and skipped
%   Reads ('switch'): results/switch_swirl1.csv (SWITCH_DRIVE(1)), results/switch_swirl1_cw.csv if present
%           (SWITCH_RUN with opts.Ebottom = 'cw'; its g-fADI rows replace those of the first file), and
%           results/switch_ref/N<N>_swirl1.mat (SWITCH_RUN 'gfadi' at the largest N).
%   Writes: figs/paper_<name>.{pdf,png}.  An error in 'switch' stops.
here = app_root();                       % repository root: results/ and figs/ live here
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', 'defaultAxesFontSize', 8, 'defaultTextFontSize', 8);
if nargin < 1, which = 'all'; end
if any(strcmp(which, {'all', 'galaxy'})), try, fgalaxy(here); catch e, fprintf(2, 'galaxy figure: %s\n', e.message); end, end
if any(strcmp(which, {'all', 'sponge'})), try, fsponge(here); catch e, fprintf(2, 'sponge figure: %s\n', e.message); end, end
if any(strcmp(which, {'all', 'big'})), try, fbig(here); catch e, fprintf(2, 'big figure: %s\n', e.message); end, end
if any(strcmp(which, {'all', 'switch'})), fswitch(here); end
end

function c = col(name)
switch name
    case 'gfadi',  c = [0.04 0.37 0.65];
    case 'lanc',   c = [0.16 0.62 0.56];
    case 'direct', c = [0.42 0.46 0.49];
    case 'expl',   c = [0.91 0.63 0.23];
    case 'class',  c = [0.82 0.29 0.36];
    case 'gal',    c = [0.56 0.37 0.64];
    case 'ink',    c = [0.15 0.15 0.15];
end
end

function sty(ax)
% the house style: full box, ticks inside, light major grid
set(ax, 'FontSize', 8, 'LineWidth', 0.6, 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.12, 'Layer', 'top', 'XColor', col('ink'), 'YColor', col('ink'));
end

function save_fig(fig, here, name)
f = fullfile(here, 'figs', ['paper_' name]);
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function c = infer(n)
% the inferno map, from its polynomial fit (no toolbox needed)
t = linspace(0, 1, n)';
C = [0.00021894 0.00165100 -0.01948090; 0.10651342 0.56395644 3.93271239; 11.60249308 -3.97285397 -15.94239411; -41.70399613 17.43639888 44.35414520; 77.16293570 -33.40235894 -81.80730926; -71.31942824 32.62606426 73.20951986; 25.13112622 -12.24266895 -23.07032500];
c = zeros(n, 3);
for k = 7:-1:1, c = c.*t + C(k, :); end
c = min(1, max(0, c));
end

function fgalaxy(here)
% Not in the paper.  (a), (b) the face-on maps at t = 700 Myr, 150 MHz and 15 GHz; (c) the light curves after a starburst
res = fullfile(here, 'results');
f = fullfile(res, 'galaxy_full.mat');
if ~exist(f, 'file'), f = fullfile(res, 'galaxy_mid.mat'); end
G = load(f);
nx = G.o.nx;
nz = G.o.nz;
rv = size(G.Z, 2);
Y = double(reshape(G.Y, G.nE, G.nt, rv));
Z = double(G.Z);
[~, jt] = min(abs(G.t - 700));
En = @(nu) sqrt(nu/96e6);
idx = @(nu) find(abs(G.En - En(nu)) == min(abs(G.En - En(nu))), 1);
face = @(nu) reshape(Z(G.ind.face, :)*squeeze(Y(idx(nu), jt, :)), nx, nx);
edge = @(nu) reshape(Z(G.ind.edge, :)*squeeze(Y(idx(nu), jt, :)), nx, nz);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 1.85], 'Visible', 'off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
cm = infer(256);
for k = 1:2
    nu = [150e6 15e9];
    ax = nexttile(tl);
    v = max(face(nu(k)), 0);
    pc = sort(v(1:7:end));
    ref = pc(round(0.995*numel(pc)));
    imagesc(ax, G.xc, G.xc, log10(max(v'/ref, 10^-2.2)));
    axis(ax, 'image');
    axis(ax, 'xy');
    colormap(ax, cm);
    caxis(ax, [-2.2 0.4]);
    set(ax, 'XLim', [-15 15], 'YLim', [-15 15], 'XTick', [-10 0 10], 'YTick', [-10 0 10], 'FontSize', 7, 'TickDir', 'out', 'Box', 'off');
    xlabel(ax, '$x$ (kpc)');
    if k == 1, ylabel(ax, '$y$ (kpc)'); end
    title(ax, sprintf('(%s) %s', char('a' + k - 1), ternary(k == 1, '150 MHz', '15 GHz')), 'FontWeight', 'normal', 'FontSize', 8);
end
% (c) light curves after a starburst: each frequency labelled next to its own peak
ax = nexttile(tl);
hold(ax, 'on');
tot = reshape(Z(G.ind.total, :)*reshape(permute(Y, [3 1 2]), rv, []), G.nE, G.nt);
j0 = find(G.t >= G.bursts.tb - 60, 1);
cc = {col('class'), col('direct'), col('gfadi')};
nu3 = [150e6 1.4e9 15e9];
lb = {'150 MHz', '1.4 GHz', '15 GHz'};
for k = 1:3
    c = tot(idx(nu3(k)), :)/tot(idx(nu3(k)), j0);
    plot(ax, G.t, c, 'Color', cc{k}, 'LineWidth', 1.1 + 0.3*(k == 3));
end
set(ax, 'XLim', [600 1200], 'XTick', [600 900 1200], 'YLim', [0.8 3]);
sty(ax);
xlabel(ax, '$t$ (Myr)');
ylabel(ax, 'total radio power');
title(ax, '(c)', 'FontWeight', 'normal', 'FontSize', 8);
text(ax, 835, 2.65, lb{3}, 'Color', cc{3}, 'FontSize', 7.5);
text(ax, 850, 2.1, lb{2}, 'Color', cc{2}, 'FontSize', 7.5);
text(ax, 1195, 1.42, lb{1}, 'Color', cc{1}, 'FontSize', 7.5, 'HorizontalAlignment', 'right');
f = fullfile(here, 'figs', 'paper_galaxy');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'image', 'Resolution', 500);
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function fsponge(here)
% Not in the paper.  (a) a picture of the Menger sponge; (b) the local exponent gamma(t) for each level, and the two limits
A = jsondecode(fileread(fullfile(here, 'results', 'menger_analysis.json')));
cs = A.curves;
if isstruct(cs), cs = num2cell(cs); end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 6.4 2.15], 'Visible', 'off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
tl.TileIndexing = 'rowmajor';
% (a) the picture, cropped to a centred square
ax = nexttile(tl, 1);
im = imread(fullfile(here, 'figs', 'thumbs', 'sponge.jpg'));
sz = size(im);
c = round(sz(2)/2);
w = min(round(0.9*sz(1)), floor(sz(2)/2));
image(ax, im(max(1, round(sz(1)/2) - w + 1):min(sz(1), round(sz(1)/2) + w), c - w + 1:c + w, :));
axis(ax, 'image', 'off');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) one curve per level, the finest in the g-fADI colour
ax = nexttile(tl, 2, [1 2]);
hold(ax, 'on');
shades = [0.72 0.72 0.72; 0.55 0.70 0.82; 0.35 0.58 0.78; 0.20 0.47 0.72; 0.10 0.40 0.68; col('gfadi')];
nc = numel(cs);
sh = shades(end - nc + 1:end, :);
for i = 1:nc
    c = cs{i};
    semilogx(ax, c.t, c.s, 'Color', sh(i, :), 'LineWidth', 0.7 + 0.5*(i == nc));
end
gp = (A.df - A.db)/A.dwExit(end);
plot(ax, [1e-7 0.5], [gp gp], '--', 'Color', col('class'), 'LineWidth', 0.8);
plot(ax, [1e-7 0.5], [0.5 0.5], ':', 'Color', col('ink'), 'LineWidth', 0.8);
set(ax, 'XScale', 'log', 'XLim', [1e-7 0.5], 'YLim', [0.2 0.7], 'XTick', 10.^(-7:0));
sty(ax);
xlabel(ax, '$t$');
ylabel(ax, 'local exponent $\gamma(t)$');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
text(ax, 2e-7, 0.515, 'solid block, $\gamma=1/2$', 'FontSize', 7, 'Color', col('ink'));
text(ax, 2e-7, gp - 0.018, sprintf('$(d_f-d_b)/d_w=%.3f$', gp), 'FontSize', 7, 'Color', col('class'));
for i = 1:nc                                                        % level labels where each curve crosses 0.62
    c = cs{i};
    [~, b] = min(abs(c.s - 0.62));
    text(ax, c.t(b)*1.25, 0.62 + 0.02*(mod(i, 2) - 0.5), sprintf('$%d$', c.L), 'FontSize', 7, 'Color', sh(i, :)*0.8);
end
f = fullfile(here, 'figs', 'paper_sponge');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'image', 'Resolution', 500);
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function fbig(here)
% Not in the paper.  (a) time and (b) peak memory against N for the graded family; Cholesky marching and
% explicit Euler are estimates (dashed)
f = fullfile(here, 'results', 'testbed_mega.csv');
if ~exist(f, 'file'), return; end
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'family', 'method', 'cls', 'note'}, 'string');
T = readtable(f, o);
T = T(T.family == "graded" & T.f == 0, :);
N = unique(T.N);
n = numel(N);
g = nan(n, 1);
gm = g;
l = g;
lm = g;
fill = g;
tch = g;
ex = g;
fr = NaN;
flr = NaN;
tf = NaN;
Nr = NaN;
for k = 1:n
    S = T(T.N == N(k), :);
    a = S(S.method == "gfadi", :);
    b = S(S.method == "lanczos", :);
    h = S(S.method == "chol", :);
    if ~isempty(a)
        g(k) = a.time(end);
        gm(k) = a.peakGB(end);
        ex(k) = sum(logspace(-4, -1, a.n(end)))*a.lmax(end)/2*a.time(end)/a.products(end)*a.rho(end);
    end
    if ~isempty(b) && b.projected(end) == 0, l(k) = b.time(end);  lm(k) = b.peakGB(end); end
    if ~isempty(h) && isfinite(h.fillGB(end))
        fr = h.fillGB(end);
        flr = h.flops(end);
        Nr = N(k);
        if isfinite(h.factor_s(end)), tf = h.factor_s(end)/h.flops(end); end
    end
    if isfinite(fr)
        fill(k) = fr*(N(k)/Nr)^4;
        if isfinite(tf), tch(k) = 200*tf*flr*(N(k)/Nr)^6; end
    end
end
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 6.4 2.3], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
% (a) time
ax = nexttile(tl);
hold(ax, 'on');
pl = @(y, c, s) loglog(ax, N(isfinite(y)), y(isfinite(y)), s, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3.5, 'LineWidth', 1.1);
pl(g, col('gfadi'), '-o');
pl(l, col('lanc'), '-s');
pl(tch, col('direct'), '--^');
pl(ex, col('expl'), '--d');
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [N(1)*0.9 N(end)*1.1], 'XTick', N');
sty(ax);
xlabel(ax, '$N$');
ylabel(ax, 'time (s)');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
legend(ax, {'g-fADI', 'Lanczos, one pass', 'Cholesky marching', 'explicit Euler'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
% (b) peak memory
ax = nexttile(tl);
hold(ax, 'on');
pl = @(y, c, s) loglog(ax, N(isfinite(y)), y(isfinite(y)), s, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3.5, 'LineWidth', 1.1);
pl(gm, col('gfadi'), '-o');
pl(lm, col('lanc'), '-s');
pl(fill, col('direct'), '--^');
plot(ax, [N(1)*0.9 N(end)*1.1], [32 32], ':', 'Color', col('ink'), 'LineWidth', 0.8);
text(ax, N(1), 42, 'memory of the machine', 'FontSize', 7);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [N(1)*0.9 N(end)*1.1], 'XTick', N');
sty(ax);
xlabel(ax, '$N$');
ylabel(ax, 'peak memory (GB)');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
save_fig(fig, here, 'big');
end

function fswitch(here, W)
% Example 3 (switching diffusion): (a) the averaged probability on the slice through the target, (b) time and (c) peak
% memory against the number of unknowns 120 N^3, from results/switch<tag>.csv (SWITCH_RUN) and results/switch_ref.
% W: the swirl (default 1, tag '_swirl1'); W = 0 is the symmetric problem (results/switch.csv, Lanczos for Arnoldi).
if nargin < 2, W = 1; end
tag = '';
if W > 0, tag = sprintf('_swirl%g', W); end
f = fullfile(here, 'results', ['switch' tag '.csv']);
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'method', 'note'}, 'string');
T = readtable(f, o);
fb = fullfile(here, 'results', ['switch' tag '_cw.csv']);                   % g-fADI with the rigorous bound (LMIN_CW), all counted
if exist(fb, 'file')
    B = readtable(fb, setvartype(detectImportOptions(fb, 'Delimiter', ','), {'method', 'note'}, 'string'));
    B = B(B.method == "gfadi" & isfinite(B.time), :);
    [~, iu] = unique(B.N, 'last');
    B = B(iu, :);
    if ~isempty(B), T = [T(~(T.method == "gfadi" & ismember(T.N, B.N)), :); B(:, T.Properties.VariableNames)]; end
    fprintf('switch figure: g-fADI rows with the rigorous bound at N = %s\n', mat2str(B.N'));
end
get = @(me) sortrows(T(T.method == me & isfinite(T.time), :), 'N');
g = get("gfadi");
l = get("arnoldi");
if W == 0, l = get("lanczos"); end
s = get("gs");
b = get("bs");
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.1 1.95], 'Visible', 'off');
tl = tiledlayout(fig, 1, 7, 'TileSpacing', 'compact', 'Padding', 'tight');
% (a) the map: the probability averaged over the states (stationary proportions), on the slice through the target
Nm = max(g.N);
R = load(fullfile(here, 'results', 'switch_ref', sprintf('N%d%s.mat', Nm, tag)));
P0 = switch_problem(8);
u = reshape(R.Xo(1:Nm^2, :)*P0.pi, Nm, Nm);
xc = ((1:Nm) - 0.5)/Nm;
[X1, X2] = ndgrid(xc, xc);
zc = (R.ks - 0.5)/Nm;
inT = (X1 - P0.o.c(1)).^2 + (X2 - P0.o.c(2)).^2 + (zc - P0.o.c(3))^2 < P0.o.R^2;
u(inT) = 1;
ax = nexttile(tl, [1 2]);
lo = floor(20*min(u(:)))/20;
imagesc(ax, xc, xc, u');
axis(ax, 'image');
axis(ax, 'xy');
colormap(ax, infer(256));
caxis(ax, [lo 1]);
hold(ax, 'on');
lv = lo + 0.05:0.05:0.95;
[cm, hc] = contour(ax, xc, xc, u', lv, 'Color', [1 1 1]*0.95, 'LineWidth', 0.4);
clabel(cm, hc, lv(abs(lv - 0.65) < 1e-9 | abs(lv - 0.75) < 1e-9), 'Color', [1 1 1]*0.95, 'FontSize', 7, 'LabelSpacing', 400, 'Interpreter', 'latex');
set(ax, 'XTick', [0 0.5 1], 'YTick', [0 0.5 1], 'FontSize', 8, 'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.6);
xlabel(ax, '$x$');
ylabel(ax, '$y$');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) time.  Measured points only: g-fADI, Arnoldi (Lanczos if W = 0) and extended Krylov, runs completed to the end
ax = nexttile(tl, [1 3]);
hold(ax, 'on');
sty(ax);
e = sortrows(T(T.method == "eksm" & isfinite(T.time), :), 'N');
ef = e(e.projected > 0, :);
e = e(e.projected == 0, :);
ser = {g, 'g-fADI', col('gfadi'), 'o';  l, ternary(W > 0, 'Arnoldi', 'Lanczos'), col('lanc'), 's';  e, 'extended Krylov', col('direct'), '^'};
xg = g.unknowns;
ytop = 3e7;
hl = gobjects(size(ser, 1), 1);                % one legend for (b) and (c), which share colours and markers
for i = 1:size(ser, 1), hl(i) = plot(ax, NaN, NaN, ['-' ser{i, 4}], 'Color', ser{i, 3}, 'MarkerFaceColor', ser{i, 3}, 'MarkerSize', 3.2, 'LineWidth', 1.0 + 0.5*(i == 1)); end
for i = 1:size(ser, 1)
    D = ser{i, 1};
    if isempty(D), continue, end
    D = D(D.projected == 0, :);
    x = D.unknowns;
    y = D.time;
    c = ser{i, 3};
    lw = 1.0 + 0.5*(i == 1);
    plot(ax, x, y, ['-' ser{i, 4}], 'Color', c, 'LineWidth', lw, 'MarkerSize', 3.2, 'MarkerFaceColor', c);
end
xl = [min(T.unknowns(T.method == "gfadi"))*0.7, max(g.unknowns)*1.6];
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', [0.1 2e4], 'XTick', 10.^(6:8), ...
        'YTick', [1 60 3600], 'YTickLabel', {'1 s', '1 min', '1 h'});
xlabel(ax, 'unknowns $120N^3$');
ylabel(ax, 'time');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
lg = legend(hl, ser(:, 2), 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8, 'Interpreter', 'latex');
lg.Layout.Tile = 'north';
lg.ItemTokenSize = [14 8];
% (c) memory: the measured peak of the process, same runs as (b)
ax = nexttile(tl, [1 2]);
hold(ax, 'on');
sty(ax);
for i = 1:size(ser, 1)
    D = ser{i, 1};
    if isempty(D), continue, end
    D = D(D.projected == 0, :);
    x = D.unknowns;
    c = ser{i, 3};
    y = D.peakGB;
    plot(ax, x, y, ['-' ser{i, 4}], 'Color', c, 'LineWidth', 1.0 + 0.5*(i == 1), 'MarkerSize', 3.2, 'MarkerFaceColor', c);
end
plot(ax, xl, [32 32], ':', 'Color', col('ink'), 'LineWidth', 0.8);
text(ax, xl(1)*1.3, 40, '32 GB', 'FontSize', 8);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', [0.8 60], 'XTick', 10.^(6:8), 'YTick', [1 10]);
xlabel(ax, 'unknowns $120N^3$');
ylabel(ax, 'peak memory (GB)');
title(ax, '(c)', 'FontWeight', 'normal', 'FontSize', 8);
f = fullfile(here, 'figs', 'paper_switch');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function place_lab(ax, ends)
% end labels, spaced so that 7.5 pt text does not overlap on a log axis
yl = ylim(ax);
u = ax.Units;
ax.Units = 'inches';
hin = ax.Position(4);
ax.Units = u;
r = 10^(1.4*7.5/72*log10(yl(2)/yl(1))/hin);
[~, i] = sort(cell2mat(ends(:, 2)));
ends = ends(i, :);
y = cell2mat(ends(:, 2));
for k = 2:numel(y), y(k) = max(y(k), y(k-1)*r); end
for k = 1:numel(y), text(ax, ends{k, 1}, y(k), ends{k, 3}, 'Color', ends{k, 4}, 'FontSize', 8); end
end

function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end
