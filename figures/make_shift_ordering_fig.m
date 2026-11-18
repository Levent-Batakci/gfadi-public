function make_shift_ordering_fig()
%MAKE_SHIFT_ORDERING_FIG  Fig. 3 of the paper (shift_ordering): g-fADI error under three step orderings.
%   For d = 1/4, 1/16, 1/32, 1/128: A and B are symmetric 300 x 300 with spectra spread over
%   [-1, -d] and [d, 1], and C = u*v' is random of rank 1.  For each q = 1..26, the zeros and poles
%   are those of the type-(q,q) Zolotarev rational for E = [-1, -d] and F = [d, 1], in closed form
%   (Wachspress), and g-fADI takes q right and q left steps: alternating, all left steps first, or
%   all right steps first.  The error ||X_q - X||_2/||X||_2 is plotted against sigma_{q+1}(X)/||X||_2.
%
%   Call:   make_shift_ordering_fig()
%   Reads:  nothing (everything is computed here).
%   Writes: figs/shift_ordering.{pdf,png} and the data in figs/shift_ordering.mat.
setup_paths();
here = app_root();
rng(1);
n = 300;
K = 26;
ds = [1/4 1/16 1/32 1/128];
orders = {'alt', 'leftfirst', 'rightfirst'};
err = NaN(numel(ds), K, numel(orders));
sv = NaN(numel(ds), K);
for i = 1:numel(ds)
    d = ds(i);
    QA = orth(randn(n));
    QB = orth(randn(n));
    A = QA*diag(-linspace(d, 1, n))*QA';
    B = QB*diag(linspace(d, 1, n))*QB';
    A = (A + A')/2;
    B = (B + B')/2;
    U = randn(n, 1);
    V = randn(n, 1);
    X = sylvester(A, -B, U*V');
    s = svd(X);
    sv(i, :) = s(2:K+1)/s(1);
    for k = 1:K
        [alpha, beta] = zolo_sym(d, k);
        for o = 1:numel(orders)
            [Z, Y] = gfadi(A, B, U, V, alpha, beta, 'order', orders{o});
            err(i, k, o) = norm(Z*Y' - X)/s(1);
        end
    end
    fprintf('d = 1/%d: alternating %.1e, left first %.1e, right first %.1e at k = %d\n', ...
            round(1/d), err(i, K, 1), err(i, K, 2), err(i, K, 3), K);
end
save(fullfile(here, 'figs', 'shift_ordering.mat'), 'ds', 'err', 'sv', 'K', 'n');

% ---------------------------------------------------------------- figure, drawn at its printed size
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');
C = [0.04 0.37 0.65; 0.20 0.60 0.24; 0.88 0.58 0.18; 0.56 0.37 0.64];   % blue, green, orange, purple (as in the caption)
mk = {'o', 's', '^'};
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4.4 2.1], 'Visible', 'off');
ax = axes(fig);
hold(ax, 'on');
k = 1:K;
for i = 1:numel(ds)
    for o = numel(orders):-1:1
        plot(ax, k, err(i, :, o), ['-' mk{o}], 'Color', C(i, :), 'MarkerFaceColor', C(i, :), ...
             'MarkerSize', 2.5, 'LineWidth', 0.8);
    end
    plot(ax, k, sv(i, :), '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7);
end
h = gobjects(4, 1);
g = [0.35 0.35 0.35];
for o = 1:numel(orders)
    h(o) = plot(ax, NaN, NaN, ['-' mk{o}], 'Color', g, 'MarkerFaceColor', g, 'MarkerSize', 3, 'LineWidth', 0.8);
end
h(4) = plot(ax, NaN, NaN, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7);
set(ax, 'YScale', 'log', 'FontSize', 8, 'LineWidth', 0.6, 'Box', 'on', 'TickDir', 'in', ...
        'TickLength', [0.012 0.012], 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'Layer', 'top', ...
        'XLim', [0.5 K + 0.5], 'XTick', 0:5:25, 'YLim', [1e-16 1e8], 'YTick', 10.^(-16:4:8), ...
        'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);
xlabel(ax, '$q$');
ylabel(ax, 'relative error');
lg = legend(ax, h, {'alternating', 'left steps first', 'right steps first', '$\sigma_{q+1}(X)/\|X\|_2$'}, ...
            'Location', 'north', 'NumColumns', 2, 'FontSize', 7.5, 'Box', 'on', 'EdgeColor', [0.75 0.75 0.75]);   % top center is empty
lg.ItemTokenSize = [14 8];
f = fullfile(here, 'figs', 'shift_ordering');
exportgraphics(fig, [f '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [f '.png'], 'Resolution', 300);
close(fig);
end

function [alpha, beta] = zolo_sym(d, k)
% type-(k,k) Zolotarev shifts for E = [-1, -d], F = [d, 1]: zeros -q_j in E, poles q_j in F, with
% q_j = dn((2j-1)K'/(2k), k'), k' = sqrt(1 - d^2) (Wachspress's optimal ADI parameters)
m = 1 - d^2;
Kp = ellipke(m);
[~, ~, q] = ellipj((2*(1:k)' - 1)*Kp/(2*k), m);
alpha = -q;
beta = q;
end
