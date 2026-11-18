function make_zolo_fig()
%MAKE_ZOLO_FIG  Fig. 2 of the paper (ms_zolo_gallery): type-(p,q) Zolotarev rationals from ZOLO_GALLERY
%   for the gears, yin-yang and combs condensers in one row, lettered (a)-(c), with the colour bar and
%   the marker key below.
%
%   Call:   make_zolo_fig()
%   Reads:  nothing (everything is computed here).
%   Writes: figs/zolo/zolo_row.{pdf,png} (vector PDF, 'raster' false).
setup_paths();
out = fullfile(app_root(), 'figs', 'zolo');
if ~exist(out, 'dir'), mkdir(out); end
old = cd(out);
c = onCleanup(@() cd(old));
% one thread: the Lawson iteration then lands on the same (symmetric) solution every time
nt = maxNumCompThreads(1);
c2 = onCleanup(@() maxNumCompThreads(nt));
zolo_gallery({'gears', 'yinyang', 'combs'}, 'row', true, 'save', true, 'raster', false, 'name', 'zolo_row', ...
             'markersize', 3.3, 'zerowidth', 0.7, 'polewidth', 1.0, 'halo', 0.55, 'linewidth', 1.2, ...
             'pad', 1.1, 'grid', 600, 'fontsize', 8);
close all
end
