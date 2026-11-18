function make_sec5_figs(only)
%MAKE_SEC5_FIGS  Fig. 6 of the paper (paper_testbed, Example 1), in the style of MAKE_SEC4_FIGS.
%   Fig. 6: solve times at the same accuracy, (a) against m = N^2 in 2-D, (b) against m = N^3 in 3-D,
%   (c) in 3-D with N = 48 against the size n of B.  Measured points only.
%
%   Call:   make_sec5_figs()            draws every figure of this file in turn: scaling, memory, testbed
%                                       (Fig. 6), testbed_forest, drift.  An error in one is printed and
%                                       the next is drawn.  Without their results files, scaling and
%                                       drift print errors; memory draws numbers typed into the code.
%           make_sec5_figs('scaling')   that figure alone (any other argument draws all of them)
%   Reads (results/):
%     testbed         testbed.csv (TB_SWEEP), and if present testbed_nsmall_meas.csv and
%                     testbed_size_meas.csv (TB_SWEEP reruns without a time budget, opts.file set to those
%                     names), whose rows replace the matching rows of testbed.csv
%     testbed_forest  testbed.csv (TB_SWEEP); not in the paper
%     scaling, drift  not in the paper; the runs that write their files are not in this repository
%     memory          not in the paper; reads no file
%   Writes: figs/paper_<name>.{pdf,png}.
here = app_root();                       % repository root: results/ and figs/ live here
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', 'defaultAxesFontSize', 8, 'defaultTextFontSize', 8);
if nargin >= 1 && strcmp(only, 'scaling'), fscaling(here); return; end
try, fscaling(here); catch e, fprintf(2, 'scaling figure: %s\n', e.message); end
try, fmemory(here); catch e, fprintf(2, 'memory figure: %s\n', e.message); end
try, ftestbed(here); catch e, fprintf(2, 'testbed figure: %s\n', e.message); end
try, ftestbed_forest(here); catch e, fprintf(2, 'testbed forest figure: %s\n', getReport(e)); end
try, fdrift(here); catch e, fprintf(2, 'drift figure: %s\n', e.message); end
end

function ftestbed(here)
% Fig. 6, figs/paper_testbed: Example 1 (results/testbed.csv), time against (a) m in 2-D, (b) m in 3-D, (c) the size n
% of the small side in 3-D.  Measured points only.
f = fullfile(here, 'results', 'testbed.csv');
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'sweep', 'param', 'method', 'note'}, 'string');
o = setvartype(o, 'value', 'double');
T = readtable(f, o);
% the n-sweep rerun with g-fADI, classical fADI and Bartels-Stewart to completion (results/testbed_nsmall_meas.csv)
% replaces those methods' rows of that sweep, so (c) compares times from one rerun, measured to the end
fm = fullfile(here, 'results', 'testbed_nsmall_meas.csv');
if exist(fm, 'file')
    om = detectImportOptions(fm, 'Delimiter', ',');
    om = setvartype(om, {'sweep', 'param', 'method', 'note'}, 'string');
    om = setvartype(om, 'value', 'double');
    Tm = readtable(fm, om);
    Tm = Tm(:, T.Properties.VariableNames);
    Tm = Tm(Tm.projected == 0, :);                                    % measured to the end only (BS at n = 3200 was projected)
    drop = T.sweep == "nsmall3d" & ismember(T.method, unique(Tm.method));
    T = [T(~drop, :); Tm];
end
% the size sweeps rerun with no time budget at the sizes where the 30 s budget had cut them off
% (results/testbed_size_meas.csv); each new row replaces the row of the same sweep, size and method, so the
% comparison is between runs from one rerun
fs = fullfile(here, 'results', 'testbed_size_meas.csv');
if exist(fs, 'file')
    os = detectImportOptions(fs, 'Delimiter', ',');
    os = setvartype(os, {'sweep', 'param', 'method', 'note'}, 'string');
    os = setvartype(os, 'value', 'double');
    Ts = readtable(fs, os);
    Ts = Ts(:, T.Properties.VariableNames);
    Ts = Ts(Ts.projected == 0, :);
    key = @(S) S.sweep + "|" + string(S.value) + "|" + S.method;
    T = [T(~ismember(key(T), key(Ts)), :); Ts];
    fprintf('testbed: %d rows from the rerun without a time budget\n', height(Ts));
end
% method key, legend name, colour, marker
M = {'gfadi', 'g-fADI', col('gfadi'), 'o'; 'lanczos2', 'Lanczos', col('lanc'), 's'; 'eksm', 'extended Krylov', col('gal'), 'v'; ...
     'fadi', 'classical fADI', col('class'), 'd'; 'bs', 'Bartels--Stewart', col('direct'), '^'};
% one row of three panels; the accuracy comparison (Peclet sweep) is Fig. 7(c), MAKE_PAPER_FIGS('bound')
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.1 1.95], 'Visible', 'off');   % printed size, full text width
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
ttl = {'(a)', '(b)', '(c)'};
sw = {'size2d', 'size3d', 'nsmall3d'};
xl = {'$m=N^2$ (2D)', '$m=N^3$ (3D)', 'size $n$ of $B$'};   % (c): 3D, N = 48, in the caption
yl = repmat({[0.1 3000]}, 1, 3);
yt = repmat({10.^(-1:3)}, 1, 3);
ytl = repmat({{'0.1', '1', '10', '100', '1000'}}, 1, 3);
for q = 1:3
    ax = nexttile(tl);
    hold(ax, 'on');
    S = T(T.sweep == sw{q}, :);
    xa = [];
    if q == 1                                             % one legend for all panels (same colours throughout)
        hl = gobjects(size(M, 1), 1);
        for k = 1:size(M, 1), hl(k) = plot(ax, NaN, NaN, ['-' M{k, 4}], 'Color', M{k, 3}, 'MarkerFaceColor', M{k, 3}, 'MarkerSize', 3.5, 'LineWidth', 1.1 + 0.6*(k == 1)); end
    end
    for k = size(M, 1):-1:1                                                       % absolute run times, g-fADI drawn last
        R = S(S.method == M{k, 1}, :);
        if isempty(R), continue; end
        x = R.m;
        if q == 3, x = R.value; end
        [x, i] = sort(x);
        t = R.time(i);
        ok = R.projected(i) == 0;                                                 % measured only
        x = x(ok);
        t = t(ok);
        xa = [xa; x(:)];   %#ok<AGROW>
        if isempty(x), continue, end
        plot(ax, x, t, ['-' M{k, 4}], 'Color', M{k, 3}, 'LineWidth', 1.1 + 0.6*(k == 1), 'MarkerFaceColor', M{k, 3}, 'MarkerSize', 3.5);
    end
    xr = [min(xa)*0.8, max(xa)*1.25];
    xt = {[1e4 1e5], [1e4 1e5], [1e2 1e3]};
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xr, 'YLim', [0.02 2000], 'XTick', xt{q});
    sty(ax);
    set(ax, 'YTick', 10.^(-1:3), 'YTickLabel', {'0.1', '1', '10', '100', '1000'});
    if q == 3, set(ax, 'YLim', [0.6 3000], 'YTick', 10.^(0:3), 'YTickLabel', {'1', '10', '100', '1000'}); end   % (c): no times below 1 s
    if q == 1, ylabel(ax, 'time (s)'); end
    xlabel(ax, xl{q});
    title(ax, ttl{q}, 'FontWeight', 'normal', 'FontSize', 8);
end
lg = legend(hl, M(:, 2), 'Orientation', 'horizontal', 'NumColumns', 5, 'Box', 'off', 'FontSize', 8, 'Interpreter', 'latex');
lg.Layout.Tile = 'north';
lg.ItemTokenSize = [16 8];
save_fig(fig, here, 'testbed');
end

function ftestbed_forest(here)
% figs/paper_testbed_forest (not in the paper): every sweep of TB_SWEEP in one picture.  For each problem, the speedup of one-sided
% g-fADI over the FASTEST other method that reaches the same accuracy (a competitor counts only if its residual is at
% most 10 times the target), on a log axis with 1 marked; one row per property varied, the dots in the order of
% the parameter, coloured by the method that came closest.  (a) three dimensions, (b) two dimensions.
f = fullfile(here, 'results', 'testbed.csv');
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'sweep', 'param', 'method', 'note', 'value'}, 'string');
T = readtable(f, o);
% competitors: method key, name, colour
C = {'lanczos2', 'Lanczos', col('lanc');  'eksm', 'extended Krylov', col('gal');  'fadi', 'classical fADI', col('class');  'bs', 'Bartels--Stewart', col('direct')};
% rows: sweep, label, which column holds the varied value ('m' or 'v' = value)
R3 = {'size3d', 'grid size $m$', 'm';  'nsmall3d', 'size $n$ of $B$', 'v';  'gap3d', 'gap between spectra', 'v';  'rank3d', 'rank $\rho$', 'v';
      'tol3d', 'target accuracy', 'v';  'Ltype3d', 'kind of $B$', 'v';  'Lnonnormal3d', 'non-normality of $B$', 'v';  'contrast3d', 'contrast in $a(x)$', 'v'};
R2 = {'size2d', 'grid size $m$', 'm';  'nsmall', 'size $n$ of $B$', 'v';  'gap', 'gap between spectra', 'v';  'rank', 'rank $\rho$', 'v';
      'tol', 'target accuracy', 'v';  'Ltype', 'kind of $B$', 'v';  'Lnonnormal', 'non-normality of $B$', 'v';  'contrast', 'contrast in $a(x)$', 'v'};
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.3 3.0], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
hl = gobjects(size(C, 1), 1);
for pnl = 1:2
    if pnl == 1, Rw = R3; else, Rw = R2; end
    ax = nexttile(tl);
    hold(ax, 'on');
    nr = size(Rw, 1);
    plot(ax, [1 1], [0.4 nr + 0.6], '-', 'Color', col('gfadi'), 'LineWidth', 1.6);
    for r = 1:nr
        S = T(T.sweep == Rw{r, 1}, :);
        y = nr + 1 - r;
        if strcmp(Rw{r, 3}, 'm'), key = string(S.m); else, key = S.value; end
        vals = unique(key, 'stable');
        num = str2double(vals);
        if all(isfinite(num)), [~, ord] = sort(num); vals = vals(ord); end
        sp = NaN(numel(vals), 1);
        who = zeros(numel(vals), 1);
        est = false(numel(vals), 1);
        for i = 1:numel(vals)
            Si = S(key == vals(i), :);
            g = Si(Si.method == "gfadi", :);
            if isempty(g), continue; end
            best = Inf;
            for k = 1:size(C, 1)
                c = Si(Si.method == C{k, 1}, :);
                if isempty(c) || ~isfinite(c.time(1)) || c.projected(1) > 0, continue; end   % measured only
                if isfinite(c.res(1)) && c.res(1) > 10*Si.tol(1), continue; end             % missed the target by more than 10x
                if c.time(1) < best, best = c.time(1);  who(i) = k;  est(i) = c.projected(1) > 0; end
            end
            sp(i) = best/g.time(1);
        end
        ok = isfinite(sp) & who > 0;
        idx = find(ok);
        s1 = sp(idx(1));
        s2 = sp(idx(end));
        % the range over the sweep
        plot(ax, [min(sp(idx)) max(sp(idx))], [y y], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
        for j = 1:numel(idx)
            i = idx(j);
            cc = C{who(i), 3};
            fc = cc;
            if est(i), fc = 'w'; end
            plot(ax, sp(i), y, 'o', 'Color', cc, 'MarkerFaceColor', fc, 'MarkerSize', 3.4, 'LineWidth', 0.8);
        end
        % the first and last values of the property, under their dots (one label when the dots are close)
        l1 = label_of(Rw{r, 1}, vals(idx(1)));
        l2 = label_of(Rw{r, 1}, vals(idx(end)));
        if abs(log10(s2/s1)) < 0.35
            text(ax, min(s1, s2)/1.15, y - 0.36, sprintf('%s $\\to$ %s', l1, l2), 'FontSize', 6.5, 'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
        else
            text(ax, s1, y - 0.36, l1, 'FontSize', 6.5, 'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
            text(ax, s2, y - 0.36, l2, 'FontSize', 6.5, 'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
        end
    end
    set(ax, 'XScale', 'log', 'XLim', [2e-3 50], 'YLim', [0.35 nr + 0.6], 'YTick', 1:nr, 'YTickLabel', flipud(Rw(:, 2)), ...
            'XTick', [1e-2 1e-1 1 10], 'XTickLabel', {'0.01', '0.1', '1', '10'});
    sty(ax);
    set(ax, 'YGrid', 'off', 'TickLength', [0.008 0.008]);
    text(ax, 1.25, nr + 0.45, 'g-fADI faster $\to$', 'FontSize', 6.5, 'Color', col('gfadi'));
    text(ax, 0.8, nr + 0.45, '$\leftarrow$ slower', 'FontSize', 6.5, 'Color', col('gfadi'), 'HorizontalAlignment', 'right');
    title(ax, sprintf('(%c) %s', 'a' + pnl - 1, ternary(pnl == 1, 'three dimensions', 'two dimensions')), 'FontWeight', 'normal', 'FontSize', 8);
    if pnl == 2, set(ax, 'YTickLabel', []); end
    if pnl == 1
        for k = 1:size(C, 1), hl(k) = plot(ax, NaN, NaN, 'o', 'Color', C{k, 3}, 'MarkerFaceColor', C{k, 3}, 'MarkerSize', 3.6); end
    end
end
lg = legend(hl(1:3), cellfun(@(s) ['fastest other: ' s], C(1:3, 2), 'UniformOutput', false), ...
            'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7.5, 'Interpreter', 'latex');   % Bartels--Stewart is never the fastest
lg.Layout.Tile = 'north';
lg.ItemTokenSize = [10 8];
xlabel(tl, 'time of the fastest other method that reaches the target, divided by the time of g-fADI', 'Interpreter', 'latex', 'FontSize', 8);
save_fig(fig, here, 'testbed_forest');
end

function s = label_of(sweep, v)
% a parameter value as a reader would say it
x = str2double(v);
sweep = char(sweep);
if ~isfinite(x), s = char(v); return, end
if startsWith(sweep, 'size'), d = 2 + endsWith(sweep, '3d');  s = sprintf('$%d^%d$', round(x^(1/d)), d); return, end
if startsWith(sweep, 'contrast'), if x == 0, s = '1'; else, s = sprintf('$10^{%d}$', x); end, return, end
if x == 0
    s = '0';
elseif abs(x) >= 1e4 || abs(x) < 1e-2
    s = sprintf('$10^{%d}$', round(log10(abs(x))));
else
    s = sprintf('%g', x);
end
end

function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end

function fdrift(here)
% figs/paper_drift (not in the paper): A = -div(a grad) + w.grad, |w| = 5 (results/testbed_adv.csv): (a) time,
% (b) memory added, against N.  Measured points only.  Memory added = peak of the process minus MATLAB with the
% problem built.
f = fullfile(here, 'results', 'testbed_adv.csv');
o = detectImportOptions(f, 'Delimiter', ',');
o = setvartype(o, {'method', 'note'}, 'string');
T = readtable(f, o);
T = T(T.W == 5, :);
% g-fADI with the rigorous Collatz-Wielandt bound on lambda_min (Ebottom 'cw', counted): replaces the earlier rows
fb = fullfile(here, 'results', 'testbed_adv_cw.csv');
if exist(fb, 'file')
    ob = setvartype(detectImportOptions(fb, 'Delimiter', ','), {'method', 'note'}, 'string');
    B = readtable(fb, ob);
    B = B(B.W == 5 & B.method == "gfadi" & isfinite(B.time), :);
    if ~isempty(B), T = [T(~(T.method == "gfadi" & ismember(T.N, B.N)), :); B(:, T.Properties.VariableNames)]; end
    fprintf('drift figure: g-fADI rows with the rigorous bound at N = %s\n', mat2str(B.N'));
end
N = [64 96 128];
get = @(me, Nk, fld) T.(fld)(find(T.method == me & T.N == Nk, 1));
tg = arrayfun(@(k) get("gfadi", k, 'time'), N);
ta = [get("arnoldi", 64, 'time') get("arnoldi", 96, 'time') NaN];
tm = arrayfun(@(k) get("gmresmg", k, 'time'), N);
pm = arrayfun(@(k) get("gmresmg", k, 'projected'), N) > 0;
ka = [get("arnoldi", 64, 'k') get("arnoldi", 96, 'k')];
ex = log(ka(2)/ka(1))/log(96/64);
% memory added, GB (typed in from the measured runs); ma(3) is an estimate, printed only
mg = [0.7 1.0 1.3];
ma = [2.0 10 10*(128/96)^(3 + ex)];
mm = [0.7 1.0 1.3];
fprintf('drift figure: Arnoldi at N = 128, about %.1f GB (k ~ N^%.2f)\n', ma(3), ex);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 1.85], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
% (a) time
ax = nexttile(tl);
hold(ax, 'on');
plot(ax, N, tg, '-o', 'Color', col('gfadi'), 'MarkerFaceColor', col('gfadi'), 'MarkerSize', 3.5, 'LineWidth', 1.3);
plot(ax, N(1:2), ta(1:2), '-s', 'Color', col('lanc'), 'MarkerFaceColor', col('lanc'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
plot(ax, N(~pm), tm(~pm), '-d', 'Color', col('expl'), 'MarkerFaceColor', col('expl'), 'MarkerSize', 3.5, 'LineWidth', 1.1);   % measured only
text(ax, 66, tm(1)*4.5, 'GMRES marching', 'Color', col('expl'), 'FontSize', 7.5);
text(ax, 99, ta(2), 'Arnoldi', 'Color', col('lanc'), 'FontSize', 7.5);
text(ax, 100, tg(2)/2.2, 'g-fADI', 'Color', col('gfadi'), 'FontSize', 7.5);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [60 136], 'YLim', [40 6000], 'XTick', N);
sty(ax);
xlabel(ax, '$N$');
ylabel(ax, 'time (s)');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) memory added
ax = nexttile(tl);
hold(ax, 'on');
plot(ax, [60 136], [32 32], ':', 'Color', col('ink'), 'LineWidth', 0.8);
text(ax, 62, 45, '32 GB', 'FontSize', 7.5);
plot(ax, N(1:2), ma(1:2), '-s', 'Color', col('lanc'), 'MarkerFaceColor', col('lanc'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
plot(ax, N, mg, '-o', 'Color', col('gfadi'), 'MarkerFaceColor', col('gfadi'), 'MarkerSize', 3.5, 'LineWidth', 1.3);
text(ax, 64.5, 4.6, 'Arnoldi', 'Color', col('lanc'), 'FontSize', 7.5);
text(ax, 100, mg(2)*1.7, 'g-fADI', 'Color', col('gfadi'), 'FontSize', 7.5);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [60 136], 'YLim', [0.3 100], 'XTick', N);
sty(ax);
xlabel(ax, '$N$');
ylabel(ax, 'memory added (GB)');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
save_fig(fig, here, 'drift');
end

function fscaling(here)
% figs/paper_scaling (not in the paper): (a) time and (b) products with A against the number of unknowns, graded family.  Measured
% points only: results/testbed_meas.csv (every method run in full); testbed_mega.csv supplies the symbolic analyses
% of Cholesky, from which the estimates printed at the end (not drawn) are made.
rd = @(f) readtable(f, setvartype(setvartype(detectImportOptions(f, 'Delimiter', ','), {'family', 'method', 'cls', 'note'}, 'string'), 'time', 'double'));
T = rd(fullfile(here, 'results', 'testbed_mega.csv'));
T = T(T.family == "graded" & T.rho == 2, :);                           % symbolic analyses only
M = rd(fullfile(here, 'results', 'testbed_meas.csv'));
M = M(M.family == "graded" & isfinite(M.time), :);
fb = fullfile(here, 'results', 'testbed_meas_cw.csv');                 % g-fADI with the rigorous bound on lambda_min
if exist(fb, 'file')
    B = rd(fb);
    B = B(B.family == "graded" & B.method == "gfadi" & isfinite(B.time), :);
    if ~isempty(B), M = [M(~(M.method == "gfadi" & ismember(M.N, B.N)), :); B(:, M.Properties.VariableNames)]; end
    fprintf('scaling figure: g-fADI rows with the rigorous bound at N = %s\n', mat2str(B.N'));
end
Mp = M(M.projected > 0, :);
M = M(M.projected == 0, :);
last = @(S, N, me) S(find(S.N == N & S.method == me, 1, 'last'), :);
unk = @(N) N.^3*200;
% g-fADI: every size run in full, in the same quiet conditions as the competitors
Ng = sort(unique(M.N(M.method == "gfadi")))';
tg = arrayfun(@(N) last(M, N, "gfadi").time, Ng);
mg = arrayfun(@(N) last(M, N, "gfadi").peakGB, Ng);
pp = arrayfun(@(N) last(M, N, "gfadi").p, Ng);
% CG marching: the sizes run to the end
Ncm = sort(M.N(M.method == "cg"))';
tcm = arrayfun(@(N) last(M, N, "cg").time, Ncm);
mcm = arrayfun(@(N) last(M, N, "cg").peakGB, Ncm);
% beyond the full runs (printed, not drawn): time from the trend of the last three of them; peak memory of the
% partial runs at those sizes, measured
Ncx = sort(Mp.N(Mp.method == "cg"))';
Ncx = Ncx(Ncx > max(Ncm));
mcx = arrayfun(@(N) last(Mp, N, "cg").peakGB, Ncx);
k3 = numel(Ncm) - 2:numel(Ncm);
cf = polyfit(log(Ncm(k3)), log(tcm(k3)), 1);
tcx = exp(polyval(cf, log(Ncx)));
fprintf('scaling: CG time ~ N^%.2f from N = %s\n', cf(1), mat2str(Ncm(k3)));
% Cholesky marching: run to the end on the small grids; beyond (printed, not drawn), scaled by the flops and fill of
% the symbolic analysis
Nhm = sort(M.N(M.method == "cholmarch"))';
thm = arrayfun(@(N) last(M, N, "cholmarch").time, Nhm);
mhm = arrayfun(@(N) last(M, N, "cholmarch").peakGB, Nhm);
N0 = Nhm(end);
r0 = last(T, N0, "chol");
Nhx = [128 256];
thx = thm(end)*arrayfun(@(N) last(T, N, "chol").flops, Nhx)/r0.flops;
mhx = mhm(end) - r0.fillGB + arrayfun(@(N) last(T, N, "chol").fillGB, Nhx);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 2.1], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
tc = [tcm tcx];
mc = [mcm mcx];
xc = unk([Ncm Ncx]);
xr = [unk(Ng(1))*0.6, 1e10];
% (a) time
ax = nexttile(tl);
hold(ax, 'on');
series(ax, unk(Nhm), thm, unk(Nhx), thx, col('direct'), '^');
series(ax, unk(Ncm), tcm, unk(Ncx), tcx, col('expl'), 's');
series(ax, unk(Ng), tg, [], [], col('gfadi'), 'o');
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xr, 'YLim', [2 2e4], 'XTick', 10.^(7:10), ...
        'YTick', [1 60 3600], 'YTickLabel', {'1 s', '1 min', '1 h'});
sty(ax);
xlabel(ax, 'unknowns $N^3 n$');
ylabel(ax, 'time');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
ends = {unk(Ng(end))*1.3, tg(end), 'g-fADI', col('gfadi');  unk(Ncm(end))*1.3, tcm(end), 'CG', col('expl')};
drawnow;
place_labels(ax, ends);
text(ax, unk(Nhm(end))*1.3, thm(end)/2, 'Cholesky', 'Color', col('direct'), 'FontSize', 7.5);     % below the CG line
set(ax, 'YTick', [1 60 3600], 'YTickLabel', {'1 s', '1 min', '1 h'});
% (b) products with A: g-fADI (its Collatz-Wielandt bound included) and CG marching; Cholesky marching takes none
ax = nexttile(tl);
hold(ax, 'on');
pg = arrayfun(@(N) last(M, N, "gfadi").products, Ng);
pc = arrayfun(@(N) last(M, N, "cg").products, Ncm);
series(ax, unk(Ncm), pc, [], [], col('expl'), 's');
series(ax, unk(Ng), pg, [], [], col('gfadi'), 'o');
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xr, 'YLim', [400 4e5], 'XTick', 10.^(7:10), 'YTick', 10.^(3:5));
sty(ax);
ends = {unk(Ncm(end))*1.3, pc(end), 'CG', col('expl')};
drawnow;
place_labels(ax, ends);
set(ax, 'YTick', 10.^(3:5));
text(ax, unk(Ng(end))/1.5, pg(end)/2.6, 'g-fADI', 'Color', col('gfadi'), 'FontSize', 7.5, 'HorizontalAlignment', 'right');   % under the line, inside the axes
xlabel(ax, 'unknowns $N^3 n$');
ylabel(ax, 'products with $A$');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
fprintf('scaling: products, g-fADI %s, CG %s (ratio %s)\n', mat2str(pg), mat2str(pc), mat2str(pc./pg(1:numel(pc)), 3));
save_fig(fig, here, 'scaling');
fprintf('scaling: g-fADI N %s: %s s, %s GB, p %s (p/(N+1) = %s); CG N %s: %s s (then %s projected); Cholesky N %s: %s s (then %s estimated, %s GB)\n', ...
        mat2str(Ng), mat2str(round(tg)), mat2str(mg, 3), mat2str(pp), mat2str(pp./(Ng + 1), 3), mat2str(Ncm), mat2str(round(tcm)), mat2str(round(tcx)), ...
        mat2str(Nhm), mat2str(round(thm)), mat2str(round(thx)), mat2str(mhx, 3));
end

function r = pick(a, b)
% the first non-empty of two table rows, as a struct (not called)
if ~isempty(a), r = table2struct(a); else, r = table2struct(b); end
end

function series(ax, x, y, xe, ye, c, mk)
% measured points: solid line, filled markers (the estimates xe, ye are accepted but not drawn)
plot(ax, x, y, ['-' mk], 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3.3, 'LineWidth', 1.0 + 0.4*strcmp(mk, 'o'));
end

function fmemory(here)
% figs/paper_memory (not in the paper): memory added (process peak minus MATLAB with the problem built), GB, by
% g-fADI and one-pass Lanczos, (a) against the outputs kept at every step and (b) against the number of time steps.
% Numbers typed in from the measured runs.
out = [2000 4096 32768 262144 2097152];                 % N = 128: 2000 cells, block means over 8^3, 4^3, 2^3, the whole field
gm  = [0.20 0.21 0.24 0.52 3.5];
lm = [0.30 0.37 1.3 9.3 28];
lest = [0 0 0 0 1];                                     % 1: an estimate (not drawn)
nt  = [1e3 1e4 1e5 1e6];                                % N = 64, 2000 cells kept, n time steps
gn  = [0.1 0.2 0.2 1.4];
ln = [0.2 0.4 0.9 7.3];
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.2 1.95], 'Visible', 'off');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
% (a) against the outputs kept
ax = nexttile(tl);
hold(ax, 'on');
loglog(ax, out, gm, '-o', 'Color', col('gfadi'), 'MarkerFaceColor', col('gfadi'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
m = ~lest;
loglog(ax, out(m), lm(m), '-s', 'Color', col('lanc'), 'MarkerFaceColor', col('lanc'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [1.4e3 3e6], 'YLim', [0.1 15], 'XTick', 10.^(3:6), 'YTick', [0.1 1 10]);
sty(ax);
text(ax, out(4), gm(4)*0.4, 'g-fADI', 'Color', col('gfadi'), 'FontSize', 7.5, 'HorizontalAlignment', 'center');
text(ax, 1.7e3, 2.6, 'one-pass Lanczos', 'Color', col('lanc'), 'FontSize', 7.5);
xlabel(ax, 'cells kept at every step ($N=128$)');
ylabel(ax, 'memory added (GB)');
title(ax, '(a)', 'FontWeight', 'normal', 'FontSize', 8);
% (b) against the number of time steps
ax = nexttile(tl);
hold(ax, 'on');
loglog(ax, nt, gn, '-o', 'Color', col('gfadi'), 'MarkerFaceColor', col('gfadi'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
loglog(ax, nt, ln, '-s', 'Color', col('lanc'), 'MarkerFaceColor', col('lanc'), 'MarkerSize', 3.5, 'LineWidth', 1.1);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [7e2 1.5e6], 'YLim', [0.07 12], 'XTick', 10.^(3:6), 'YTick', [0.1 1 10]);
sty(ax);
text(ax, nt(3)/3, 0.135, 'g-fADI', 'Color', col('gfadi'), 'FontSize', 7.5, 'HorizontalAlignment', 'center');
text(ax, 1.1e3, 0.75, 'one-pass Lanczos', 'Color', col('lanc'), 'FontSize', 7.5);
xlabel(ax, 'time steps $n$ ($N=64$, 2000 cells kept)');
ylabel(ax, 'memory added (GB)');
title(ax, '(b)', 'FontWeight', 'normal', 'FontSize', 8);
save_fig(fig, here, 'memory');
end

function place_labels(ax, ends, r)
% labels at the line ends {x, y, text, colour}, pushed apart so that neighbouring labels are at least a factor r apart
% (log axis); r omitted or []: the factor that keeps 7.5 pt labels clear, from the height of the axes in inches.
% Also widens the y-limits to fit the top label and resets the y-ticks to at most about five decades.
yl = ylim(ax);
u = ax.Units;
ax.Units = 'inches';
hin = ax.Position(4);
ax.Units = u;
rmin = 10^(1.4*7.5/72*log10(yl(2)/yl(1))/hin);
if nargin < 3 || isempty(r), r = rmin; else, r = max(r, rmin); end
[~, i] = sort(cell2mat(ends(:, 2)));
ends = ends(i, :);
y = cell2mat(ends(:, 2));
for k = 2:numel(y), y(k) = max(y(k), y(k-1)*r); end
if y(end)*sqrt(r) > yl(2), ylim(ax, [yl(1) y(end)*sqrt(r)*1.05]); end
yl = ylim(ax);
e = ceil(log10(yl(1))):floor(log10(yl(2)));
e = e(1:ceil(numel(e)/5):end);
set(ax, 'YTick', 10.^e);
for k = 1:numel(y), text(ax, ends{k, 1}, y(k), ends{k, 3}, 'Color', ends{k, 4}, 'FontSize', 7.5); end
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
