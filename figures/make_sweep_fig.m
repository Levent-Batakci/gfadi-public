function make_sweep_fig(N, omegas)
%MAKE_SWEEP_FIG  Fig. 10 of the paper (paper_sweep): the frequency sweep of Example 4, at printed size
%   (text width 5.1 in, 8 pt text), in two rows:
%     top     (a) the composite material, (b)-(c) |u|/max|u| on the plane z = 1/2 at two frequencies,
%             with isolines
%     bottom  (d) run time and (e) peak memory against N: one-sided g-fADI, multishift FOM, FOM(200).
%
%   Call:   make_sweep_fig()  or  make_sweep_fig(N, omegas)   (defaults N = 127, omegas = [3e2 3e3])
%   Reads:  figs/paper_intro_3d.png (MAKE_INTRO_3D, run it first), results/intro_slices_N<N>.mat
%           (INTRO_SLICES_GFADI) and results/freq_sweep_C1000.csv (FREQ_SWEEP_RUN with opts.C = 1e3).
%   Writes: figs/paper_sweep.pdf, composed with pdflatex (which must be on the system path) from a vector
%           layer and full-resolution images (see the end of the main function), and a preview
%           figs/paper_sweep.png.
if nargin < 1, N = 127; end
if nargin < 2, omegas = [3e2 3e3]; end
here = app_root();
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', ...
    'defaultAxesFontSize', 8, 'defaultTextFontSize', 8, 'defaultColorbarTickLabelInterpreter', 'latex');
S = load(fullfile(here, 'results', sprintf('intro_slices_N%d.mat', N)));
render = imread(fullfile(here, 'figs', 'paper_intro_3d.png'));
render = render(round(0.09*size(render, 1)):end, :, :);            % drop the "(a)" printed in the image
T = readtable(fullfile(here, 'results', 'freq_sweep_C1000.csv'), 'TextType', 'string');
ok = ~startsWith(T.note, "ERROR");
g = sortrows(T(ok & T.method == "gfadi_out", :), 'N');  a = sortrows(T(ok & T.method == "msfom_out", :), 'N');
r = sortrows(T(T.method == "rsfom_out" & contains(T.note, "FOM(200)") & ok, :), 'N');
% layout in inches [left bottom width height]
W = 5.1;  H = 3.5;  s = 1.38;  yt = 1.92;                          % top row: square panels of side s at height yt
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 W H], 'Visible', 'off');
put = @(bx) axes(fig, 'Units', 'inches', 'Position', bx);
ax = put([0.0 yt s s]);  image(ax, render, 'Tag', 'raster');  axis(ax, 'image', 'off');
wr = size(render, 2)/size(render, 1);                              % where axis image puts the render, in inches
if wr <= 1, ra = [(s - s*wr)/2, yt, s*wr, s]; else, ra = [0, yt + (s - s/wr)/2, s, s/wr]; end
text(ax, 0, 0, '(a)', 'Units', 'normalized', 'Position', [0.5 1.02], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
L = 'bc';  xs = [1.62 3.12];  M = 1024;  cm = parula(256);  RGB = cell(1, 2);
for k = 1:2
    [~, j] = min(abs(S.omegas - omegas(k)));
    [Fs, gs] = slice_field(S.U(:, :, j), S.N, M);
    RGB{k} = ind2rgb(flipud(min(256, floor((Fs.' + 4)/4*256) + 1)), cm);   % top row of the PNG is y = 1
    ax = put([xs(k) yt s s]);  slice_panel(ax, Fs, gs, 1/(S.N + 1));
    title(ax, ['(' L(k) ') $\omega=' texnum(omegas(k)) '$'], 'FontWeight', 'normal');
end
cb = colorbar(ax);  cb.Units = 'inches';  cb.Position = [xs(2) + s + 0.06, yt, 0.07, s];  ax.Position = [xs(2) yt s s];
cb.Ticks = -4:2:0;  cb.TickLabels = {'$10^{-4}$', '$10^{-2}$', '$1$'};  cb.FontSize = 8;
% bottom row: wide, short panels
yb = 0.36;  hb = 1.02;
ax = put([0.45 yb 2.05 hb]);  h = timeplot(ax, g, a, r, 'time');  title(ax, '(d)', 'FontWeight', 'normal');
ax = put([3.0 yb 2.05 hb]);   timeplot(ax, g, a, r, 'peakGB');  title(ax, '(e)', 'FontWeight', 'normal');
lg = legend(h, {'g-fADI$\quad$', 'Arnoldi (multishift FOM)$\quad$', 'restarted FOM(200)'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8);
lg.ItemTokenSize = [14 8];  drawnow;  lg.Units = 'inches';  lp = lg.Position;
lg.Position = [(W - lp(3))/2, yb + hb + 0.2, lp(3), lp(4)];       % between the rows, centred
f = fullfile(here, 'figs', 'paper_sweep');
exportgraphics(fig, [f '.png'], 'Resolution', 250);                % preview, everything drawn
% PDF in two layers, because MATLAB downsamples images in vector output to screen resolution (about 110 dpi):
%  1. everything except the three images, as vectors on a transparent background (axes, isolines, circles, text, colorbar)
%  2. the render and the two slices as full-resolution PNGs
% and a small LaTeX file that puts the PNGs under the vector layer at their exact positions in inches.
set(findobj(fig, 'Tag', 'raster'), 'Visible', 'off');
set(fig, 'Color', 'none');  set(findobj(fig, 'Type', 'axes'), 'Color', 'none');
exportgraphics(fig, [f '_vec.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none', 'Padding', 'figure');
close(fig);
imwrite(render, [f '_a.png']);
gr = @(x, y, w, hh, file) ['\put(' sprintf('%.4f,%.4f', x, y) '){\includegraphics[width=' sprintf('%.4f', w) 'in,height=' sprintf('%.4f', hh) 'in]{' file '}}'];
pl = {gr(ra(1), ra(2), ra(3), ra(4), 'paper_sweep_a.png')};
for k = 1:2
    imwrite(RGB{k}, sprintf('%s_%s.png', f, L(k)));
    pl{end + 1} = gr(xs(k), yt, s, s, ['paper_sweep_' L(k) '.png']); %#ok<AGROW>
end
tex = [{'\documentclass{article}', ['\usepackage[paperwidth=' sprintf('%g', W) 'in,paperheight=' sprintf('%g', H) 'in,margin=0pt]{geometry}'], ...
        '\usepackage{graphicx}', '\pagestyle{empty}', '\setlength{\unitlength}{1in}', '\begin{document}', '\noindent', ...
        ['\begin{picture}(' sprintf('%g,%g', W, H) ')']}, pl, ...
       {gr(0, 0, W, H, 'paper_sweep_vec.pdf'), '\end{picture}', '\end{document}'}];
fid = fopen([f '_compose.tex'], 'w');  fprintf(fid, '%s\n', tex{:});  fclose(fid);
old = cd(fullfile(here, 'figs'));  [st, out] = system('pdflatex -interaction=nonstopmode paper_sweep_compose.tex');  cd(old);
if st ~= 0, error('make_sweep_fig: pdflatex failed:\n%s', out); end
movefile([f '_compose.pdf'], [f '.pdf'], 'f');
fprintf('make_sweep_fig: wrote %s.pdf (vector layer + %d-px slices)\n', f, size(RGB{1}, 1));
end

function [Fs, gs] = slice_field(u, N, M)
% log10|u|/max|u| on M x M pixel centres covering [0,1]^2 (bicubic, display only), clamped at -4; u = 0 on the
% boundary (Dirichlet), so the grid is padded with the floor value there and the image fills the whole box
h = 1/(N + 1);  g = [0, (1:N)*h, 1];
F = -4*ones(N + 2);  F(2:end-1, 2:end-1) = max(log10(abs(u)/max(abs(u(:)))), -4);
gs = ((1:M) - 0.5)/M;  [Gx, Gy] = ndgrid(gs, gs);
Fs = max(min(interp2(g, g, F.', Gy, Gx, 'cubic').', 0), -4);      % Fs(i, j) at x = gs(i), y = gs(j)
end

function slice_panel(ax, Fs, gs, h)
hold(ax, 'on');
imagesc(ax, gs, gs, Fs.', 'Tag', 'raster');
% isolines at 1e-1, 1e-2, 1e-3, drawn from a copy smoothed over 0.75 grid cells (display only: the bicubic surface
% follows the kinks of the grid solution at the inclusion edges), and not within 1.5 cells of the walls, where the
% Dirichlet boundary layer would put a line on top of the frame
sg = 0.75*h*numel(gs);  t = -ceil(3*sg):ceil(3*sg);  kk = exp(-t.^2/(2*sg^2));  kk = kk/sum(kk);  p = numel(t);
ix = [ones(1, p), 1:numel(gs), numel(gs)*ones(1, p)];  Fp = Fs(ix, ix);   % replicate the edges
Fc = conv2(kk, kk, Fp, 'same');  Fc = Fc(p+1:end-p, p+1:end-p);
for lev = -3:-1                                     % draw a level unless it traces the Dirichlet boundary layer
    C = contourc(gs, gs, Fc.', [lev lev]);  k = 1;  pts = zeros(2, 0);
    while k < size(C, 2), n = C(2, k);  pts = [pts, C(:, k+1:k+n)];  k = k + n + 1; end %#ok<AGROW>
    dw = min(min(pts, [], 1), min(1 - pts, [], 1));                % distance of each point to the nearest wall
    if isempty(pts) || mean(dw < 2*h) > 0.05, continue, end         % a level that runs along the walls is skipped
    contour(ax, gs, gs, Fc.', [lev lev], 'Color', [0 0 0], 'LineWidth', 0.4);
end
axis(ax, 'xy', 'image');  colormap(ax, parula(256));  caxis(ax, [-4 0]);
ctr = 0.25*[kron([1;1;1], kron([1;1;1], (1:3)')), kron([1;1;1], kron((1:3)', [1;1;1])), kron((1:3)', kron([1;1;1], [1;1;1]))];
ctr = ctr + 0.04*sin((1:27)'*[1.3 2.1 3.7]);  r0 = 0.11;  th = linspace(0, 2*pi, 100);   % the inclusions of FREQ_SWEEP_RUN
for i = find(abs(ctr(:, 3) - 0.5) < r0)'
    rr = sqrt(r0^2 - (ctr(i, 3) - 0.5)^2);  plot(ax, ctr(i, 1) + rr*cos(th), ctr(i, 2) + rr*sin(th), '-', 'Color', 'w', 'LineWidth', 0.5);
end
set(ax, 'XTick', [], 'YTick', [], 'XLim', [0 1], 'YLim', [0 1], 'Box', 'on');
end

function h = timeplot(ax, g, a, r, fld)
cg = [0.04 0.37 0.65];  ca = [0.56 0.37 0.64];  cr = [0.42 0.46 0.49];  ink = [0.15 0.15 0.15];
hold(ax, 'on');  box(ax, 'on');
set(ax, 'LineWidth', 0.6, 'TickDir', 'in', 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridAlpha', 0.12, 'Layer', 'top', 'XColor', ink, 'YColor', ink);
h(1) = plot(ax, g.N, g.(fld), '-o', 'Color', cg, 'MarkerFaceColor', cg, 'MarkerSize', 3, 'LineWidth', 1.4);
h(2) = plot(ax, a.N, a.(fld), '-s', 'Color', ca, 'MarkerFaceColor', ca, 'MarkerSize', 3, 'LineWidth', 1.0);
h(3) = plot(ax, r.N, r.(fld), '-^', 'Color', cr, 'MarkerFaceColor', cr, 'MarkerSize', 3, 'LineWidth', 1.0);
Nt = unique([g.N; a.N]).';
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [min(Nt)/1.12, max(Nt)*1.12], 'XTick', Nt, 'XMinorTick', 'off', 'YMinorTick', 'off');
if strcmp(fld, 'time')
    set(ax, 'YLim', [5 3600], 'YTick', [10 60 600 3600], 'YTickLabel', {'10 s', '1 min', '10 min', '1 h'});   % ticks carry the unit
else
    plot(ax, get(ax, 'XLim'), [32 32], ':', 'Color', ink, 'LineWidth', 0.8);
    set(ax, 'YLim', [0.3 60], 'YTick', [1 10 32]);  ylabel(ax, 'memory (GB)');
end
xlabel(ax, '$N$');
end

function s = texnum(v)
p = floor(log10(v));  c = v/10^p;
if abs(c - 1) < 1e-12, s = sprintf('10^{%d}', p); else, s = sprintf('%g\\cdot10^{%d}', c, p); end
end
