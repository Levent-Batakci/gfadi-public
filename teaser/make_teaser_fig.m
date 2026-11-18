function make_teaser_fig(outdir)
%MAKE_TEASER_FIG  Fig. 1: (a) run time and (b) peak memory for the model problem of Section 1.1, against m = N^3.
%
%   MAKE_TEASER_FIG()  or  MAKE_TEASER_FIG(OUTDIR)
%   Reads teaser/teaser.csv (written by TEASER_RUN through teaser/drive.sh, e.g. drive.sh 16 24 32 40 48).
%   Writes OUTDIR/teaser_draft.{pdf,png} (default OUTDIR: the repository's figs/); the paper includes this
%   file as paper_teaser.  Measured points only (rows with projected = 0).
here = fileparts(mfilename('fullpath'));
if nargin < 1, outdir = fullfile(here, '..', 'figs'); end
f = fullfile(here, 'teaser.csv');
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'method', 'note'}, 'string');
T = readtable(f, o);
M = {'bs', 'Bartels--Stewart', col('direct'), '^'; 'fadi', 'classical fADI', col('class'), 'd'; ...
     'eksm', 'extended Krylov', col('gal'), 'v'; 'arnoldi', 'block Arnoldi', col('lanc'), 's'; ...
     'gfadi', 'one-sided g-fADI', col('gfadi'), 'o'};
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.0 1.8], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
fld = {'time', 'peakGB'};
yl = {'time (s)', 'peak memory (GB)'};
ttl = {'(a)', '(b)'};
Ns = unique(T.N);
Nl = Ns(Ns ~= 40);                                              % 40^3 and 48^3 collide at this width
for q = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(size(M, 1), 1);
    for k = 1:size(M, 1)
        R = sortrows(T(T.method == M{k, 1} & T.projected == 0, :), 'N');
        g = k == size(M, 1);
        h(k) = plot(ax, R.m, R.(fld{q}), ['-' M{k, 4}], 'Color', M{k, 3}, 'MarkerFaceColor', M{k, 3}, ...
                    'MarkerSize', 3.5, 'LineWidth', 1.1 + 0.6*g);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [Ns(1)^3/1.4, Ns(end)^3*1.4], 'XTick', Nl.^3, ...
        'XTickLabel', compose('%d^3', Nl));
    sty(ax);
    if q == 1, set(ax, 'YLim', [0.08 2000], 'YTick', 10.^(-1:3), 'YTickLabel', {'0.1', '1', '10', '100', '1000'});
    else, set(ax, 'YLim', [0.01 6], 'YTick', 10.^(-2:0), 'YTickLabel', {'0.01', '0.1', '1'}); end
    xlabel(ax, 'grid points $m$', 'Interpreter', 'latex');
    ylabel(ax, yl{q});
    title(ax, ttl{q}, 'FontWeight', 'normal', 'FontSize', 8);
end
lg = legend(flipud(h), flipud(M(:, 2)), 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8, 'Interpreter', 'latex', 'NumColumns', 3);
lg.Layout.Tile = 'north';
lg.ItemTokenSize = [16 8];
f = fullfile(outdir, 'teaser_draft');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
end
function c = col(name)
switch name
    case 'gfadi',  c = [0.04 0.37 0.65];
    case 'lanc',   c = [0.16 0.62 0.56];
    case 'direct', c = [0.42 0.46 0.49];
    case 'class',  c = [0.82 0.29 0.36];
    case 'gal',    c = [0.56 0.37 0.64];
    case 'ink',    c = [0.15 0.15 0.15];
end
end
function sty(ax)
set(ax, 'FontSize', 8, 'LineWidth', 0.6, 'TickDir', 'in', 'TickLength', [0.012 0.012], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.12, 'Layer', 'top', 'XColor', col('ink'), 'YColor', col('ink'));
end
