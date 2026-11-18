function switch_drive(W)
%SWITCH_DRIVE  All runs of Example 3 (switching diffusion), each SWITCH_RUN in its own MATLAB process.
%   Each child process starts in the repository root (app_root()), runs setup_paths() and then
%   SWITCH_RUN(N, method, opts); its log goes to logs/switch[_swirl<W>]_N<N>_<method>.log.  Rows already in
%   the results file (same N and method) are skipped, so the driver can be restarted.
%
%   SWITCH_DRIVE or SWITCH_DRIVE(0)  the symmetric problem, results/switch.csv (Krylov competitor: block Lanczos)
%   SWITCH_DRIVE(W), W > 0           the swirl (nonsymmetric A), results/switch_swirl<W>.csv (competitor: Arnoldi)
%   Fig. 9 uses SWITCH_DRIVE(1) (MAKE_SEC4_FIGS('switch')): Arnoldi, extended Krylov, and g-fADI with the
%   rigorous spectral bound (results/switch_swirl1_cw.csv).
if nargin < 1, W = 0; end
setup_paths();                                     % also creates logs/ and results/
app = app_root();
S = {};
if W > 0
    sw = sprintf('''swirl'', %g', W);
    csv = fullfile(app, 'results', sprintf('switch_swirl%g.csv', W));
    for N = [24 32 48 64 96 128 160], S{end + 1} = {N, 'gfadi', sprintf('struct(%s)', sw)}; end              %#ok<AGROW>
    S{end + 1} = {16, 'gfadi', sprintf('struct(%s)', sw)};
    S{end + 1} = {16, 'check', sprintf('struct(%s)', sw)};
    for N = [24 32 48 64], S{end + 1} = {N, 'bs', sprintf('struct(%s)', sw)}; end                        %#ok<AGROW>
    for N = [24 32 48 64 96 128], S{end + 1} = {N, 'arnoldi', sprintf('struct(%s, ''tcap'', 3600)', sw)}; end   %#ok<AGROW>
    for N = [16 24 32], S{end + 1} = {N, 'gs', sprintf('struct(%s, ''tcap'', 120)', sw)}; end         %#ok<AGROW>
    for N = [16 24 32 48 64 96], S{end + 1} = {N, 'eksm', sprintf('struct(%s, ''tcap'', 7200)', sw)}; end   %#ok<AGROW>
    % g-fADI with the rigorous lower bound on the spectrum (LMIN_CW): the g-fADI points of Fig. 9, written to
    % results/switch_swirl<W>_cw.csv, so the skip check uses that file
    cwcsv = fullfile(app, 'results', sprintf('switch_swirl%g_cw.csv', W));
    for N = [16 32 48 64 96 128 160], S{end + 1} = {N, 'gfadi', sprintf('struct(%s, ''Ebottom'', ''cw'')', sw), cwcsv}; end   %#ok<AGROW>
else
    csv = fullfile(app, 'results', 'switch.csv');
    for N = [24 32 48 64 96 128 192], S{end + 1} = {N, 'gfadi', 'struct()'}; end            %#ok<AGROW>
    for N = [24 32], S{end + 1} = {N, 'check', 'struct()'}; end                                  %#ok<AGROW>
    for N = [24 32 48 64], S{end + 1} = {N, 'bs', 'struct()'}; end                               %#ok<AGROW>
    for N = [24 32 48 64], S{end + 1} = {N, 'gs', 'struct(''tcap'', 300)'}; end                  %#ok<AGROW>
    for N = [24 32 48 64 96 128 192], S{end + 1} = {N, 'lanczos', 'struct(''tcap'', 3600)'}; end   %#ok<AGROW>
end
for i = 1:numel(S)
    N = S{i}{1};
    me = S{i}{2};
    f = csv;  sfx = '';
    if numel(S{i}) > 3, f = S{i}{4};  sfx = '_cw'; end              % runs that write their own results file
    key = sprintf('(?m)^%d,"%s",', N, me);
    if exist(f, 'file') && ~isempty(regexp(fileread(f), key, 'once')), fprintf('skip %d %s%s\n', N, me, sfx); continue; end
    lg = fullfile(app, 'logs', sprintf('switch%s_N%d_%s%s.log', ternary(W > 0, sprintf('_swirl%g', W), ''), N, me, sfx));
    t = tic;
    fprintf('%s  %s N=%d ...\n', char(datetime('now', 'Format', 'HH:mm:ss')), me, N);
    cmd = sprintf('matlab -batch "cd(''%s''); setup_paths(); switch_run(%d, ''%s'', %s);" -logfile "%s"', app, N, me, S{i}{3}, lg);
    [st, ~] = system(cmd);
    fprintf('   exit %d after %.0f s\n', st, toc(t));
end
end

function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end
