function T = tb_sweep(names, opts)
%TB_SWEEP  Example 1 (basic test case): every method on every configuration of each sweep, results to CSV.
%   T = TB_SWEEP(NAMES, OPTS) runs the sweeps in the cell array NAMES (default: all, see SWEEPS
%   below).  Each sweep varies one parameter around the base problem: 2-D Laplacian, N = 256;
%   L = 1-D diffusion, n = 200, spectrum [1, 1e4]; rank 2, random U; g-fADI target 1e-8.
%   OPTS: .tmax [60 s per method]  .threads [4]  .file ['testbed.csv']
%         .methods [{'gfadi','mixed','fadi','bs','eksm','lanczos2','arnoldi'}]  .values [{}] (run only these)
%   Appends one row per (configuration, method) to app_root()/results/<file> as it goes, so a stop
%   loses nothing.  Columns: sweep,param,value,method,m,n,rho,tol,time,tdesign,its,res,err,memGB,projected,note.
%   T is the table of the rows written.
%   g-fADI runs first and sets the residual target of the Krylov methods (equal accuracy); a completed
%   Bartels-Stewart run is the reference for every method's error.
%   Figures (results/testbed.csv):
%     Fig. 6     sweeps 'size2d', 'size3d', 'nsmall3d'  (MAKE_SEC5_FIGS('testbed'), which also merges reruns
%                written with OPTS.file = 'testbed_size_meas.csv' or 'testbed_nsmall_meas.csv' if present)
%     Fig. 7(c)  sweep 'Lnonnormal'  (MAKE_PAPER_FIGS('bound'))
%   Example:  tb_sweep({'size2d', 'size3d', 'nsmall3d', 'Lnonnormal'})
if nargin < 1 || isempty(names), names = {}; end
if nargin < 2, opts = struct(); end
o = struct('tmax', 60, 'threads', 4, 'file', 'testbed.csv', 'values', {{}}, ...
           'methods', {{'gfadi', 'mixed', 'fadi', 'bs', 'eksm', 'lanczos2', 'arnoldi'}});
fn = fieldnames(opts);
for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
app = app_root();
setup_paths();
maxNumCompThreads(o.threads);
S = sweeps();
if isempty(names), names = fieldnames(S)'; end
out = fullfile(app, 'results', o.file);
rows = {};
for sname = names
    sw = S.(sname{1});
    if ~isempty(o.values), sw.values = o.values; end                 % run only these values
    for v = 1:numel(sw.values)
        spec = sw.base;
        spec.(sw.param) = sw.values{v};
        if isfield(sw, 'tie'), spec = sw.tie(spec); end
        tol = 1e-8;
        if isfield(spec, 'tol'), tol = spec.tol; spec = rmfield(spec, 'tol'); end
        P = tb_problem(spec);
        fprintf('\n=== [%s] %s = %s :: %s, tol %g\n', sname{1}, sw.param, val2str(sw.values{v}), P.name, tol);
        meths = o.methods;
        if isfield(sw, 'methods'), meths = intersect(meths, sw.methods, 'stable'); end
        R = struct([]);
        Xref = [];
        fill = NaN;
        restarget = tol;
        % order: g-fADI (sets the residual target), Bartels-Stewart (reference), the rest
        ord = [intersect({'gfadi'}, meths), intersect({'bs'}, meths), setdiff(meths, {'gfadi', 'bs'}, 'stable')];
        for k = 1:numel(ord)
            mname = ord{k};
            op = struct('restarget', restarget, 'tmax', o.tmax, 'Xref', Xref, 'fillGB', fill);
            r = tb_method(mname, P, tol, op);
            if strcmp(mname, 'gfadi') && isfinite(r.res) && r.res <= 10*tol, restarget = max(r.res, 1e-3*tol); end   % else: the a-priori target
            if strcmp(mname, 'bs') && isfield(r, 'Xfull'), Xref = r.Xfull; end
            if strcmp(mname, 'eksm') && isfield(r, 'fillGB'), fill = r.fillGB; end
            fprintf('  %-9s time %8.3f s (design %.2f)  its %6g  res %8.1e  err %8.1e  mem %6.3f GB  %s%s\n', ...
                    mname, r.time, r.tdesign, r.its, r.res, r.err, r.memGB, r.note, ternary(r.projected, ' [projected]', ''));
            R = [R, rmfield(r, intersect(fieldnames(r), {'Z', 'Y', 'Xfull', 'fillGB'}))];   %#ok<AGROW>
            % keep the g-fADI factors to compute its error once a reference exists
            if strcmp(mname, 'gfadi'), Zg = r.Z; Yg = r.Y; end
            if strcmp(mname, 'bs') && ~isempty(Xref) && exist('Zg', 'var') && ~isempty(Zg)
                R(1).err = norm(Zg*Yg' - Xref, 'fro')/norm(Xref, 'fro');
                fprintf('  (gfadi error vs reference: %.1e)\n', R(1).err);
            end
        end
        clear Zg Yg
        for k = 1:numel(R)
            row = struct('sweep', sname{1}, 'param', sw.param, 'value', val2str(sw.values{v}), 'method', R(k).method, ...
                         'm', P.m, 'n', P.n, 'rho', size(P.U, 2), 'tol', tol, 'time', R(k).time, 'tdesign', R(k).tdesign, ...
                         'its', R(k).its, 'res', R(k).res, 'err', R(k).err, 'memGB', R(k).memGB, ...
                         'projected', R(k).projected, 'note', R(k).note);
            rows{end+1} = row;                                                        %#ok<AGROW>
            append_csv(out, row);
        end
    end
end
T = struct2table([rows{:}]);
end

function S = sweeps()
b = struct('dim', 2, 'N', 256, 'coef', 'lap', 'Ltype', 'tri', 'n', 200, 'f', 1, 'g', 1e4, 'rho', 2, 'rhs', 'rand');
S.size2d   = struct('base', b, 'param', 'N', 'values', {{64, 128, 256, 512}});
% the mixed type (q > 0 left steps) needs a type-(p,q) shift design that takes minutes for
% large p, so it has its own sweep and is left out of the others
S.mixed2d  = struct('base', b, 'param', 'N', 'values', {{64, 128, 256}}, 'methods', {{'gfadi', 'mixed', 'fadi', 'eksm'}});
b3 = b;  b3.dim = 3;
S.size3d   = struct('base', b3, 'param', 'N', 'values', {{16, 24, 32, 48, 64}});
S.gap      = struct('base', b, 'param', 'f', 'values', {{0, 1, 100, 1e4, 1e6}}, 'tie', @(s) setfield(s, 'g', s.f + 1e4)); %#ok<SFLD>
S.Ltype    = struct('base', b, 'param', 'Ltype', 'values', {{'diag', 'tri', 'time', 'dense'}});
ba = b;  ba.Ltype = 'adv';
S.Lnonnormal = struct('base', ba, 'param', 'pe', 'values', {{0, 10, 100, 1000}});
S.nsmall   = struct('base', b, 'param', 'n', 'values', {{50, 200, 800, 3200}});
bd = b;  bd.Ltype = 'dense';
S.ndense   = struct('base', bd, 'param', 'n', 'values', {{100, 400, 1600}});
S.rank     = struct('base', b, 'param', 'rho', 'values', {{1, 4, 16, 32}});
S.tol      = struct('base', b, 'param', 'tol', 'values', {{1e-4, 1e-6, 1e-8, 1e-10, 1e-12}});
S.rhs      = struct('base', b, 'param', 'rhs', 'values', {{'rand', 'smooth', 'point'}});
bh = b;  bh.coef = 'het';
S.contrast = struct('base', bh, 'param', 'contrast', 'values', {{0, 2, 4, 6}});
bw = b;  bw.coef = 'adv';
S.advA     = struct('base', bw, 'param', 'w', 'values', {{0, 30, 300, 3000}});
b3t = b3;  b3t.N = 48;
S.Ltype3d  = struct('base', b3t, 'param', 'Ltype', 'values', {{'diag', 'tri', 'time', 'dense'}});
b3r = b3;  b3r.N = 48;
S.rank3d   = struct('base', b3r, 'param', 'rho', 'values', {{1, 4, 16}});
% the same questions on the 3-D base (N = 48), where shifted solves with A are expensive
S.gap3d    = struct('base', b3r, 'param', 'f', 'values', {{0, 100, 1e4}}, 'tie', @(s) setfield(s, 'g', s.f + 1e4)); %#ok<SFLD>
b3a = b3r;  b3a.Ltype = 'adv';
S.Lnonnormal3d = struct('base', b3a, 'param', 'pe', 'values', {{0, 100, 1000}});
S.nsmall3d = struct('base', b3r, 'param', 'n', 'values', {{50, 800, 3200}});
S.tol3d    = struct('base', b3r, 'param', 'tol', 'values', {{1e-4, 1e-8, 1e-12}});
b3h = b3r;  b3h.coef = 'het';
S.contrast3d = struct('base', b3h, 'param', 'contrast', 'values', {{0, 2}});
end

function s = val2str(v)
if ischar(v), s = v; else, s = num2str(v); end
end

function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end

function append_csv(f, row)
newfile = ~exist(f, 'file');
fid = fopen(f, 'a');
fns = fieldnames(row);
if newfile, fprintf(fid, '%s\n', strjoin(fns', ',')); end
vals = cell(1, numel(fns));
for i = 1:numel(fns)
    x = row.(fns{i});
    if ischar(x), vals{i} = ['"' strrep(x, '"', '''') '"'];
    elseif islogical(x), vals{i} = num2str(double(x));
    else, vals{i} = sprintf('%.6g', x); end
end
fprintf(fid, '%s\n', strjoin(vals, ','));
fclose(fid);
end
