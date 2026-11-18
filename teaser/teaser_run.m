function teaser_run(N, method, outfile)
%TEASER_RUN  One (N, method) run of the model problem of Section 1.1 (Fig. 1).
%   A X + X L^T = U V^T  (A X - X B = U V*, B = -L^T), A the centered 7-point matrix of -Laplacian + w.grad
%   on N^3 interior points of the unit cube, constant drift w = 5 (1, 0.6, 0.3)/|(1, 0.6, 0.3)|;
%   L = I + (1e4 - 1)/4 tridiag(-1, 2, -1), n = 200; random rank-2 U, V; target 1e-8.
%
%   TEASER_RUN(N, METHOD, OUTFILE) appends one CSV row (N, m, method, time, peakGB, its, res, projected,
%   note) to OUTFILE.
%     METHOD  'gfadi'   one-sided g-fADI, real zeros from ADV_SHIFTS on a region containing W(A), Leja order
%             'fadi'    classical fADI, closed-form type-(k,k) Zolotarev shifts (CLASSICAL_SHIFTS_CLOSED)
%             'bs', 'eksm', 'arnoldi'   Bartels-Stewart, extended Krylov, block Arnoldi (TB_METHOD)
%   Peak memory is the operating system's high-water mark above what is held after the problem is built.
%   Linux: VmHWM after a reset.  Windows: PeakWorkingSet64, which cannot be reset, so the note records the
%   slack (start-up peak above the post-build working set): a method whose extra memory is below that slack
%   reads as the slack.
%   Equal accuracy: 'eksm' and 'arnoldi' stop at the residual g-fADI reached at the same N (read from
%   OUTFILE, so g-fADI must run first, as drive.sh does).
%   Meant for its own MATLAB process (teaser/drive.sh).  Example:  teaser_run(32, 'gfadi', 'teaser.csv')
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_paths'));
maxNumCompThreads(8);
tol = 1e-8;  n = 200;  rho = 2;
h = 1/(N + 1);
e = ones(N, 1);  I = speye(N);
T = spdiags([-e 2*e -e], -1:1, N, N)/h^2;
Dc = spdiags([-e 0*e e], -1:1, N, N)/(2*h);
w = 5*[1 0.6 0.3]/norm([1 0.6 0.3]);
H = kron(I, kron(I, T)) + kron(I, kron(T, I)) + kron(T, kron(I, I));
K = w(1)*kron(I, kron(I, Dc)) + w(2)*kron(I, kron(Dc, I)) + w(3)*kron(Dc, kron(I, I));
A = H + K;
m = N^3;
en = ones(n, 1);
L = speye(n) + (1e4 - 1)/4*spdiags([-en 2*en -en], -1:1, n, n);
rng(1);
U = randn(m, rho);  V = randn(n, rho);
U = U/norm(U, 'fro');  V = V/norm(V, 'fro');
E = 12/h^2*[sin(pi*h/2)^2, cos(pi*h/2)^2];             % spectrum of H, closed form
WL = [1, 1e4];                                         % Re W(L) (L symmetric, Gershgorin)
P = struct('A', A, 'L', L, 'U', U, 'V', V, 'E', E, 'WL', WL, 'm', m, 'n', n);
clear H K T Dc
note = '';
opt = struct('tmax', 1e7, 'memcap', 24, 'nores', true);
if any(strcmp(method, {'eksm', 'arnoldi', 'lanczos2'}))
    opt.restarget = gfadi_res(outfile, N);                      % equal accuracy (see header)
    opt.trunc = min(1e-12, opt.restarget/100);                  % keep the final truncation below the target
    note = sprintf('target %.2e (g-fADI res) ', opt.restarget);
end
reset_peak();
base = rss_now();
slack = peak_rss() - base;
t0 = tic;
switch method
    case 'gfadi'
        td = tic;
        [p, al] = adv_shifts(E, norm(w), sum(abs(w))/h, WL, tol);
        al = leja_order(al, E');
        tdes = toc(td);
        [Z, Y, info] = gfadi_lowrank(gfadi_op(A), gfadi_op(-L.'), U, V, al, []);
        time = toc(t0);
        its = p;
        note = sprintf('p=%d s=%d design=%.2fs', p, info.s, tdes);
        projected = 0;
    case 'fadi'
        td = tic;
        [k, al, be, sg] = classical_shifts_closed(E, [-WL(2), -WL(1)], tol);
        [al, be, sc] = order_shifts(al, be, E, [-WL(2), -WL(1)]);
        tdes = toc(td);
        [Z, Y] = gfadi_lowrank(gfadi_op(A), gfadi_op(-L.'), U, V, al, be, 'sched', sc);
        time = toc(t0);
        its = 2*k;  projected = 0;
        note = sprintf('k=%d sigma=%.1e design=%.2fs (closed-form Zolotarev shifts)', k, sg, tdes);
    otherwise
        r = tb_method(method, P, tol, opt);
        time = toc(t0);
        its = r.its;  projected = r.projected;  note = [note r.note];
        Z = r.Z;  Y = r.Y;
        if isfield(r, 'Xfull'), Z = r.Xfull;  Y = eye(n); end
end
pk = peak_rss() - base;
if ispc, note = sprintf('%s; mem slack %.3f GB', note, slack); end
% the exact residual is a check on the result, not part of any method: it is computed after
% the time and the peak memory are taken
res = NaN;
if ~isempty(Z), res = tbres(A, L, U, V, Z, Y); end
fprintf('N=%d %s: time %.2f s, peak +%.3f GB, res %.2e, %s\n', N, method, time, pk, res, note);
fid = fopen(outfile, 'a');
fprintf(fid, '%d,%d,%s,%.4f,%.4f,%d,%.3e,%d,"%s"\n', N, m, method, time, pk, its, res, projected, note);
fclose(fid);
end

function reset_peak()
if ispc, return; end                                            % Windows: the peak cannot be reset
fid = fopen('/proc/self/clear_refs', 'w');  fprintf(fid, '5');  fclose(fid);
end
function g = peak_rss()
if ispc, p = System.Diagnostics.Process.GetCurrentProcess();  g = double(p.PeakWorkingSet64)/2^30;  return; end
g = read_status('VmHWM');
end
function g = rss_now()
if ispc, p = System.Diagnostics.Process.GetCurrentProcess();  g = double(p.WorkingSet64)/2^30;  return; end
g = read_status('VmRSS');
end
function t = gfadi_res(outfile, N)
% the residual g-fADI reached at this N (its last row in OUTFILE)
txt = splitlines(fileread(outfile));
t = NaN;
for i = 1:numel(txt)
    c = strsplit(txt{i}, ',');
    if numel(c) >= 7 && strcmp(c{3}, 'gfadi') && str2double(c{1}) == N, t = str2double(c{7}); end
end
if ~(t > 0), error('teaser_run:target', 'no g-fADI residual for N=%d in %s: run gfadi first', N, outfile); end
end
function g = read_status(key)
txt = fileread('/proc/self/status');
tk = regexp(txt, [key ':\s+(\d+) kB'], 'tokens', 'once');
g = str2double(tk{1})/2^20;
end
function res = tbres(A, L, U, V, Z, Y)
[~, R1] = qr([A*Z, Z, -U], 0);
[~, R2] = qr([Y, L*Y, V], 0);
res = norm(R1*R2', 'fro')/sqrt(abs(trace((U'*U)*(V'*V))));
end
