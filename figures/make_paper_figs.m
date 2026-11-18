function make_paper_figs(which)
%MAKE_PAPER_FIGS  Figs. 7 and 8 of the paper (Example 2), in the paper's style: LaTeX text, one palette.
%
%   Call:   make_paper_figs('bound')   Fig. 7 (paper_bound): (a) error and bound along a one-sided run,
%                                      (b) increasing nonnormality of B, (c) Example 1 with an
%                                      advection-diffusion L, against the Peclet number
%           make_paper_figs('cost')    Fig. 8 (paper_cost): (a) run time and (b) memory against N^3 n
%           make_paper_figs()          both
%   Reads (results/), by figure:
%     bound    exp18_analysis.mat (EXP18_ANALYSIS), exp16b_N32.csv (EXP16B_NONNORMALITY), testbed.csv
%              (TB_SWEEP, sweep 'Lnonnormal')
%     cost     exp16_blowout.mat, exp16_blowout_gfadi.csv, and if present exp16_blowout_ext.mat,
%              exp16_blowout_ext96.mat, exp16_blowout_ext40.mat (EXP16_BLOWOUT with those tags; the
%              exact calls are in the header of EXP16_BLOWOUT)
%   Writes: figs/paper_<name>.{pdf,png} (with the suffix in getappdata(groot, 'figsuffix') if set).
setup_paths();
here = app_root();                       % repository root: results/ and figs/ live here
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex', 'defaultAxesFontSize', 8, 'defaultTextFontSize', 8);
F = {'cost', @fcost; 'bound', @fbound};
if nargin > 0, F = F(strcmp(F(:, 1), which), :); end
for i = 1:size(F, 1)
    F{i, 2}(here);
    fprintf('paper_%s done\n', F{i, 1});
end
end

% ============================================================================ style
function c = col(name)
switch name
    case 'gfadi',  c = [0.04 0.37 0.65];      % deep blue
    case 'class',  c = [0.82 0.29 0.36];      % red
    case 'gal',    c = [0.56 0.37 0.64];      % purple
    case 'direct', c = [0.42 0.46 0.49];      % grey
    case 'lanc',   c = [0.16 0.62 0.56];      % teal
    case 'eig',    c = [0.91 0.63 0.23];      % amber
    case 'ink',    c = [0.15 0.15 0.15];
end
end

function fig = newfig(w, h)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 w h], 'Visible', 'off');
end

function sty(ax)
% the house style of the whole paper: a full box, ticks inside, a light major grid
set(ax, 'FontSize', 8, 'LineWidth', 0.6, 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Box', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.12, 'Layer', 'top', ...
        'XColor', col('ink'), 'YColor', col('ink'));
end

function lab(x, y, s, c, varargin)
% a 7.5 pt LaTeX label in colour c
text(x, y, s, 'Color', c, 'FontSize', 7.5, 'Interpreter', 'latex', varargin{:});
end

function save_fig(fig, here, name, raster)
% figs/paper_<name>.{pdf,png}; raster = true writes the pdf as a 600 dpi image
if nargin < 4, raster = false; end
f = fullfile(here, 'figs', ['paper_' name]);
if isappdata(groot, 'figsuffix'), f = [f getappdata(groot, 'figsuffix')]; end
if raster
    exportgraphics(fig, [f '.pdf'], 'ContentType', 'image', 'Resolution', 600);
else
    exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
end
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function S = loadr(here, name)
S = load(fullfile(here, 'results', name));
end

% ============================================================================ Fig. 8: the cost of Example 2
function fcost(here)
S = loadr(here, 'exp16_blowout.mat');
R = S.R;
% g-fADI rerun (tag '_gfadi'): its time and memory, and the direct time where it was measured
G = readtable(fullfile(here, 'results', 'exp16_blowout_gfadi.csv'));
for i = 1:numel(R)
    k = find(G.N == R(i).N, 1);
    if ~isempty(k)
        R(i).t_gfadi = G.t_gfadi(k);
        R(i).mem_gfadi_GB = G.mem_gfadi_GB(k);
        if ~isnan(G.t_direct(k)), R(i).t_direct = G.t_direct(k); end
    end
end
% runs that extend the competitors to larger N (tags '_ext', '_ext96'): their measured times and memory replace the
% NaNs of the original runs
fx = {'t_direct', 't_classical', 't_gal', 'k_gal', 'mem_gal_GB', 'fill_GB', 't_fact'};
for tg = {'exp16_blowout_ext.mat', 'exp16_blowout_ext96.mat'}
    fe = fullfile(here, 'results', tg{1});
    if ~exist(fe, 'file'), continue, end
    Ex = load(fe);
    Ex = Ex.R;
    for i = 1:numel(Ex)
        k = find([R.N] == Ex(i).N, 1);
        if isempty(k), continue, end
        for f = fx, if isfield(Ex, f{1}) && isfinite(Ex(i).(f{1})), R(k).(f{1}) = Ex(i).(f{1}); end, end
    end
end
% runs at new sizes (tag '_ext40', N = 16 and 40, every method in one MATLAB process): added as new rows
fe = fullfile(here, 'results', 'exp16_blowout_ext40.mat');
if exist(fe, 'file')
    Ex = load(fe);
    Ex = Ex.R;
    for i = 1:numel(Ex)
        if any([R.N] == Ex(i).N), continue, end
        r = R(1);
        for f = fieldnames(r)', if isfield(Ex, f{1}), r.(f{1}) = Ex(i).(f{1}); else, r.(f{1}) = NaN; end, end
        R(end + 1) = r;  %#ok<AGROW>
    end
    [~, o] = sort([R.N]);
    R = R(o);
end
N = [R.N];
m = [R.m];
n = R(1).n;
unk = [R.unknowns];
rho = 2;
% Measured points only: each method up to the largest size it was run.  Memory: the direct method holds its factor
% and the n solution columns, classical fADI its factor, Galerkin its basis (2 rho k_gal vectors), g-fADI as
% measured.  The direct method is computed here but not drawn.
fl = [R.fill_GB];
kg = [R.k_gal];
tdir = [R.t_direct];
tcl = [R.t_classical];
tgx = [R.t_gal];
memd = fl + m*n*8/2^30;
memd(isnan(tdir)) = NaN;
memc = fl;
memc(isnan(tcl)) = NaN;
memg = 2*m.*rho.*kg*8/2^30;
memg(isnan(tgx)) = NaN;
memf = [R.mem_gfadi_GB];
fig = newfig(5.1, 1.85);   % printed size, full text width
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
series = {tcl, memc, 'class', 'classical';
          tgx, memg, 'gal', 'Galerkin';  [R.t_gfadi], memf, 'gfadi', 'g-fADI'};
lastof = @(v) v(find(isfinite(v), 1, 'last'));
xlast = @(v) unk(find(isfinite(v), 1, 'last'));
% (a) time
ax = nexttile(tl);
hold(ax, 'on');
sty(ax);
for i = 1:3
    t = series{i, 1};
    me = isfinite(t);
    c = col(series{i, 3});
    plot(ax, unk(me), t(me), '-o', 'Color', c, 'LineWidth', 1.1 + 0.6*(i == 3), 'MarkerSize', 3, 'MarkerFaceColor', c);
end
xr = [3e6 2e9];
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xr, 'YLim', [0.1 1e4], 'YTick', [1 60 3600], ...
        'YTickLabel', {'1 s', '1 min', '1 h'}, 'XTick', [1e7 1e8 1e9]);
hl = gobjects(3, 1);
for i = 1:3
    c = col(series{i, 3});
    hl(i) = plot(ax, NaN, NaN, '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3, 'LineWidth', 1.1 + 0.6*(i == 3));
end
lgd = legend(hl, {'classical fADI', 'Galerkin', 'g-fADI'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8, 'Interpreter', 'latex');
lgd.Layout.Tile = 'north';
lgd.ItemTokenSize = [14 8];
xlabel(ax, 'unknowns $N^3n$');
ylabel(ax, 'time');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) memory
ax = nexttile(tl);
hold(ax, 'on');
sty(ax);
plot(ax, xr, [32 32], ':', 'Color', col('ink'), 'LineWidth', 0.8);
lab(4e6, 60, '32 GB', col('ink'));
for i = 1:3
    c = col(series{i, 3});
    y = series{i, 2};
    me = isfinite(y);
    plot(ax, unk(me), y(me), '-o', 'Color', c, 'LineWidth', 1.1 + 0.6*(i == 3), 'MarkerSize', 3, 'MarkerFaceColor', c);
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xr, 'YLim', [0.01 1e3], 'XTick', [1e7 1e8 1e9], 'YTick', [0.01 0.1 1 10 100]);
xlabel(ax, 'unknowns $N^3n$');
ylabel(ax, 'memory (GB)');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
save_fig(fig, here, 'cost');
end

% ============================================================================ Fig. 7: the bound and nonnormality
function fbound(here)
S = loadr(here, 'exp18_analysis.mat');
R = S.R(1);                                                            % the age problem, eps = 0.02, N = 24
fig = newfig(4.36, 2.1);   % printed size, 0.85 text width
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% (a) the measured error and the bound against the step
ax = nexttile(tl);
hold(ax, 'on');
sty(ax);
k = find(isfinite(R.err2));                                            % the steps where the error was measured
plot(ax, [0 R.p], [1e-8 1e-8], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
% the bound oscillates from step to step; its running maximum over the later steps is still a bound, and monotone
bk = R.bnd(k);
env = flipud(cummax(flipud(bk(:))));
plot(ax, k, env, '-', 'Color', col('ink'), 'LineWidth', 1.0);
plot(ax, k, R.err2(k), '-', 'Color', col('gfadi'), 'LineWidth', 1.8);
lab(R.p*0.42, env(find(k >= R.p*0.42, 1))*6, 'bound', col('ink'), 'VerticalAlignment', 'bottom');
lab(R.p*0.03, 1e-10, 'measured error', col('gfadi'), 'VerticalAlignment', 'top');
lab(8, 1e-8/2.2, 'target accuracy', [0.45 0.45 0.45], 'VerticalAlignment', 'top');
set(ax, 'YScale', 'log', 'XLim', [0 R.p], 'YLim', [1e-14 1e3], 'YTick', 10.^(-12:4:0));
xlabel(ax, 'step $k$');
ylabel(ax, '$\|X_k-X\|_2/\|X\|_2$');
grid(ax, 'off');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) increasing nonnormality, on a narrow window about the target: does each method meet the target it was
% designed for?
T = sortrows(readtable(fullfile(here, 'results', 'exp16b_N32.csv')), 'eps', 'descend');
ep = T.eps;
tol = 1e-8;
yl = [2e-9 1e-6];
xl = [0.0018 0.115];
ax = nexttile(tl);
hold(ax, 'on');
sty(ax);
plot(ax, xl, [tol tol], ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
lab(0.0021, tol*1.15, 'target accuracy', [0.45 0.45 0.45], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
plot(ax, ep, T.err_classical, '-d', 'Color', col('class'), 'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', col('class'));
plot(ax, ep, T.err_gfadi, '-o', 'Color', col('gfadi'), 'LineWidth', 1.8, 'MarkerSize', 3.5, 'MarkerFaceColor', col('gfadi'));
lab(0.1, 4e-7, 'classical fADI', col('class'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
lab(0.03, 5e-9, 'g-fADI', col('gfadi'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
tk = [0.1 0.01 0.002];
set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse', 'XLim', xl, 'YLim', yl, 'XTick', fliplr(tk), ...
        'XTickLabel', arrayfun(@(e) sprintf('%g', e), fliplr(tk), 'UniformOutput', false), 'XTickLabelRotation', 0, 'YTick', 10.^(-9:-6));
xlabel(ax, {'$\epsilon$', '(more nonnormal to the right)'});
ylabel(ax, 'relative error');
grid(ax, 'off');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
% (c) the same comparison on Example 1: L an advection-diffusion matrix (centered diffusion, upwind
% advection), against its Peclet number (results/testbed.csv, sweep Lnonnormal: 2-D, N = 256, n = 200)
f = fullfile(here, 'results', 'testbed.csv');
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'sweep', 'param', 'method', 'note'}, 'string');
o = setvartype(o, 'value', 'double');
Tb = readtable(f, o);
Tb = Tb(Tb.sweep == "Lnonnormal" & Tb.value > 0, :);
ax = nexttile(tl);
hold(ax, 'on');
sty(ax);
plot(ax, [7 1500], [tol tol], ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
lab(1300, tol/1.6, 'target accuracy', [0.45 0.45 0.45], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
Rc = sortrows(Tb(Tb.method == "fadi", :), 'value');
Rg = sortrows(Tb(Tb.method == "gfadi", :), 'value');
plot(ax, Rc.value, Rc.res, '-d', 'Color', col('class'), 'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', col('class'));
plot(ax, Rg.value, Rg.res, '-o', 'Color', col('gfadi'), 'LineWidth', 1.8, 'MarkerSize', 3.5, 'MarkerFaceColor', col('gfadi'));
lab(9, 8e-5, 'classical fADI', col('class'), 'VerticalAlignment', 'top');
lab(12, 1.3e-12, 'g-fADI', col('gfadi'), 'VerticalAlignment', 'bottom');
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [7 1500], 'YLim', [1e-12 1e-4], 'XTick', [10 100 1000], 'YTick', 10.^(-12:2:-4));
xlabel(ax, {'P\''eclet number of $B$', '(Example 1)'});
ylabel(ax, 'relative residual');
grid(ax, 'off');
title(ax, '(c)', 'FontWeight', 'normal', 'FontSize', 8);
save_fig(fig, here, 'bound');
end
