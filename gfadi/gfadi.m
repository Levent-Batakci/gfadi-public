function [Z, Y, info] = gfadi(A, B, U, V, alpha, beta, varargin)
%GFADI  Generalized factored ADI (Algorithm 2.1) for  A X - X B = U V'.
%   [Z,Y,INFO] = GFADI(A,B,U,V,ALPHA,BETA,...) returns X ~ Z*Y' after numel(ALPHA)
%   right steps (product with A, shifted solve with B) and numel(BETA) left steps
%   (shifted solve with A, product with B), starting from X0 = 0.  This is type (p,q)
%   with p = numel(ALPHA) zeros and q = numel(BETA) poles.
%   A (m x m), B (n x n): matrices, or operator structs with fields .mv .mvh .sv .svh
%   (see GFADI_OP).  U (m x rho), V (n x rho).
%   Options:  'order'  'alt' (default, compressed pairs) | 'leftfirst' | 'rightfirst'
%             'sched'  explicit step string, e.g. 'RRLRRL' (overrides 'order')
%             'tol'    stop once ||R||_F/||F||_F < tol (default 0: use all shifts)
%             'errfun' @(Z,Y) scalar, evaluated after every step
%   INFO: .Zt,.Yt residual factors (R = Zt*Yt'), .res residual history, .err, .steps,
%         .nsolA .nmvA .nsolB .nmvB (columns pushed through each operation), .time
%   Example:  [Z, Y, info] = gfadi(A, B, U, V, alpha, beta, 'tol', 1e-10);

p = inputParser;
p.addParameter('order', 'alt');
p.addParameter('sched', '');
p.addParameter('tol', 0);
p.addParameter('errfun', []);
p.parse(varargin{:});
opt = p.Results;

if ~isstruct(A), A = gfadi_op(A); end
if ~isstruct(B), B = gfadi_op(B); end
alpha = alpha(:);  beta = beta(:);
np = numel(alpha);  nq = numel(beta);
rho = size(U, 2);
m = size(U, 1);  n = size(V, 1);

% ------------------------------------------------------------- step schedule
% each entry is {kind, index into beta, index into alpha}; kind 'A' = one left
% and one right step merged into a compressed alternating step
if ~isempty(opt.sched)
    opt.order = 'sched';
end
switch lower(opt.order)
    case 'sched'
        % explicit order; shifts are consumed in the order given, adjacent L,R (or R,L)
        % pairs are merged into one compressed step
        sc = opt.sched;
        ja = 0;  jb = 0;  i = 1;
        sched = cell(0, 1);
        assert(sum(sc == 'R') == np && sum(sc == 'L') == nq, 'gfadi:sched', ...
               'schedule does not match the numbers of shifts');
        while i <= numel(sc)
            if i < numel(sc) && sc(i) ~= sc(i+1)
                jb = jb + 1;
                ja = ja + 1;
                sched{end+1} = {'A', jb, ja};   %#ok<AGROW>
                i = i + 2;
            elseif sc(i) == 'L'
                jb = jb + 1;
                sched{end+1} = {'L', jb, 0};    %#ok<AGROW>
                i = i + 1;
            else
                ja = ja + 1;
                sched{end+1} = {'R', 0, ja};    %#ok<AGROW>
                i = i + 1;
            end
        end
    case 'alt'
        k = min(np, nq);
        sched = cell(0, 1);
        for j = 1:k, sched{end+1} = {'A', j, j}; end                         %#ok<AGROW>
        for j = k+1:nq, sched{end+1} = {'L', j, 0}; end                      %#ok<AGROW>
        for j = k+1:np, sched{end+1} = {'R', 0, j}; end                      %#ok<AGROW>
    case 'leftfirst'
        sched = [arrayfun(@(j) {'L', j, 0}, 1:nq, 'UniformOutput', false), ...
                 arrayfun(@(j) {'R', 0, j}, 1:np, 'UniformOutput', false)];
    case 'rightfirst'
        sched = [arrayfun(@(j) {'R', 0, j}, 1:np, 'UniformOutput', false), ...
                 arrayfun(@(j) {'L', j, 0}, 1:nq, 'UniformOutput', false)];
    otherwise
        error('gfadi:order', 'unknown order');
end
ns = numel(sched);

if isreal(U) && isreal(V) && isreal(alpha) && isreal(beta)
    Z = zeros(m, rho*ns);  Y = zeros(n, rho*ns);      % (complex B still makes Y complex on assignment)
else
    Z = complex(zeros(m, rho*ns));  Y = complex(zeros(n, rho*ns));
end
Zt = U;  Yt = V;                        % residual R_0 = F = U V'
nF = sqrt(abs(trace((U'*U)*(V'*V))));

info.res = zeros(1, ns);
info.steps = char(zeros(1, ns));
info.err = [];
cnt = struct('sA', 0, 'mA', 0, 'sB', 0, 'mB', 0);
tstart = tic;
used = 0;  ncol = 0;
for s = 1:ns
    kind = sched{s}{1};  jb = sched{s}{2};  ja = sched{s}{3};
    switch kind
        case 'A'      % alternating step: one left + one right, compressed
            b = beta(jb);  a = alpha(ja);  g = b - a;
            Zn = g*A.sv(b, Zt);                      cnt.sA = cnt.sA + rho;
            Yn = B.svh(a, Yt);                       cnt.sB = cnt.sB + rho;
            Zt = (A.mv(Zn) - a*Zn)/g;                cnt.mA = cnt.mA + rho;
            Yt = B.mvh(Yn) - conj(b)*Yn;             cnt.mB = cnt.mB + rho;
        case 'L'      % left step
            b = beta(jb);
            Zn = A.sv(b, Zt);                        cnt.sA = cnt.sA + rho;
            Yn = Yt;
            Yt = B.mvh(Yn) - conj(b)*Yn;             cnt.mB = cnt.mB + rho;
            Zt = Zn;
        case 'R'      % right step
            a = alpha(ja);
            Zn = -Zt;
            Yn = B.svh(a, Yt);                       cnt.sB = cnt.sB + rho;
            Zt = -(A.mv(Zn) - a*Zn);                 cnt.mA = cnt.mA + rho;
            Yt = Yn;
    end
    % Rescale the residual factors (R = Zt*Yt' is unchanged) so that long one-sided
    % runs do not overflow in Zt and underflow in Yt.  Zt is kept at unit norm.
    sc = norm(Zt, 'fro');
    if sc > 0 && isfinite(sc), Zt = Zt/sc;  Yt = Yt*sc; end
    Z(:, ncol+1:ncol+rho) = Zn;
    Y(:, ncol+1:ncol+rho) = Yn;
    ncol = ncol + rho;
    used = s;
    info.steps(s) = kind;
    info.res(s) = sqrt(abs(trace((Zt'*Zt)*(Yt'*Yt))))/nF;
    if ~isempty(opt.errfun)
        info.err(s) = opt.errfun(Z(:, 1:ncol), Y(:, 1:ncol));  %#ok<AGROW>
    end
    if opt.tol > 0 && info.res(s) < opt.tol, break, end
end
Z = Z(:, 1:ncol);  Y = Y(:, 1:ncol);
info.res = info.res(1:used);
info.steps = info.steps(1:used);
info.Zt = Zt;  info.Yt = Yt;
info.nsolA = cnt.sA;  info.nmvA = cnt.mA;  info.nsolB = cnt.sB;  info.nmvB = cnt.mB;
info.time = toc(tstart);
end
