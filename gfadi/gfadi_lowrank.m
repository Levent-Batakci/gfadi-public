function [Zc, Yc, info] = gfadi_lowrank(A, B, U, V, alpha, beta, varargin)
%GFADI_LOWRANK  g-fADI with on-the-fly compression: X ~ Zc*Yc' in O(m s) memory.
%   [ZC,YC,INFO] = GFADI_LOWRANK(A,B,U,V,ALPHA,BETA,...) runs the same iteration as
%   GFADI, but accumulates X = Zc*Yc' with Yc (n x s) an orthonormal basis of the
%   small side, grown only by directions that matter for X.  Memory is m*s instead of
%   m*rho*(p+q).  With right steps only, this is the small-side compression of
%   Algorithm 4.1 of the paper with all of X kept (Gamma = I).
%   A, B, U, V, ALPHA, BETA: as in GFADI.
%   Options:  'sched'    step string as in GFADI (default: 'RL' pairs, then the rest)
%             'ytol'     basis tolerance (default 1e-13)
%             'ytolmode' 'solution' (default): relative to ||X||;
%                        'increment': relative to the current increment
%             'trunc', 'every', 'verbose' are accepted and ignored.
%   INFO: .res residual history, .s basis size, .rank numerical rank (1e-8),
%         .dropped bound on ||X_computed - X_iterate||_F/||X||_F, .time, .tcomp, .maxcols
%   Example:  [Zc, Yc, info] = gfadi_lowrank(A, B, U, V, alpha, [], 'sched', repmat('R', 1, numel(alpha)));
p = inputParser;
p.addParameter('sched', '');
p.addParameter('ytol', 1e-13);
p.addParameter('trunc', 1e-13);                 % unused, accepted for compatibility
p.addParameter('ytolmode', 'solution');
p.addParameter('every', []);                    % unused, accepted for compatibility
p.addParameter('verbose', false);
p.parse(varargin{:});
opt = p.Results;
if ~isstruct(A), A = gfadi_op(A); end
if ~isstruct(B), B = gfadi_op(B); end
alpha = alpha(:);  beta = beta(:);
np = numel(alpha);  nq = numel(beta);
rho = size(U, 2);  m = size(U, 1);  n = size(V, 1);

% ---- step schedule, as in GFADI ('A' = merged left+right step)
sc = opt.sched;
if isempty(sc), k = min(np, nq); sc = [repmat('RL', 1, k), repmat('L', 1, nq - k), repmat('R', 1, np - k)]; end
sched = cell(0, 1);
ja = 0;  jb = 0;  i = 1;
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
ns = numel(sched);
cplx = ~(isreal(U) && isreal(V) && isreal(alpha) && isreal(beta));
Zt = U;  Yt = V;
nF = sqrt(abs(trace((U'*U)*(V'*V))));
cap = min(n, 64);  s = 0;                       % storage capacity, basis size
Qy = zeros(n, cap);  Zacc = zeros(m, cap);
if cplx, Qy = complex(Qy); Zacc = complex(Zacc); end
info.res = zeros(1, ns);
nbuf = 50;                                      % increments buffered between accumulator updates
Zb = zeros(m, nbuf*rho);
if cplx, Zb = complex(Zb); end
Cb = cell(1, nbuf);  nb = 0;
solmode = strcmp(opt.ytolmode, 'solution');
nXest = 0;  dropped = 0;
tstart = tic;
for st = 1:ns
    kind = sched{st}{1};  jb = sched{st}{2};  ja = sched{st}{3};
    switch kind
        case 'A'
            b = beta(jb);  a = alpha(ja);  g = b - a;
            Zn = g*A.sv(b, Zt);
            Yn = B.svh(a, Yt);
            Zt = (A.mv(Zn) - a*Zn)/g;
            Yt = B.mvh(Yn) - conj(b)*Yn;
        case 'L'
            b = beta(jb);
            Zn = A.sv(b, Zt);
            Yn = Yt;
            Yt = B.mvh(Yn) - conj(b)*Yn;
            Zt = Zn;
        case 'R'
            a = alpha(ja);
            Zn = -Zt;
            Yn = B.svh(a, Yt);
            Zt = -(A.mv(Zn) - a*Zn);
            Yt = Yn;
    end
    % rescale the residual factors (Zt*Yt' unchanged): avoids overflow/underflow
    c = norm(Zt, 'fro');
    if c > 0 && isfinite(c), Zt = Zt/c;  Yt = Yt*c; end

    % ---- accumulate Zn*Yn' in the basis Qy of the small side
    % project Yn onto Qy (twice, classical Gram-Schmidt with reorthogonalisation)
    Qs = Qy(:, 1:s);
    C = Qs'*Yn;  Yp = Yn - Qs*C;
    C2 = Qs'*Yp;  Yp = Yp - Qs*C2;
    C = C + C2;
    nyn = norm(Yn, 'fro');  nzn = norm(Zn, 'fro');  nyp = norm(Yp, 'fro');
    if solmode, nXest = max(nXest, nzn*nyn); thr = opt.ytol*nXest/max(nzn, realmin); else, thr = 0; end
    if nyp > max(opt.ytol*nyn, thr) && s < n
        % rank-revealing SVD: an unpivoted QR would drop components in later columns
        [Qn, Sn] = svd(Yp, 'econ');
        keep = diag(Sn) > max(opt.ytol*nyn, thr);
        Qn = Qn(:, keep);
        % re-orthogonalise against the basis (otherwise it slowly loses orthogonality
        % and grows without bound)
        for pass = 1:2, Qn = Qn - Qs*(Qs'*Qn);  [Qn, ~] = qr(Qn, 0); end
        t = min(size(Qn, 2), n - s);
        Qn = Qn(:, 1:t);
        if t > 0
            if s + t > cap                                  % grow storage
                nc = min(n, max(2*cap, s + t));
                Qy(:, nc) = 0;  Zacc(:, nc) = 0;  cap = nc;
            end
            Cn = Qn'*Yp;  Yp = Yp - Qn*Cn;
            Qy(:, s+1:s+t) = Qn;
            C = [C; Cn];                                    %#ok<AGROW>
            s = s + t;
        end
    end
    dropped = dropped + nzn*norm(Yp, 'fro');                % part of the increment not kept
    % X += Zn*Yn' = Zn*C'*Qy', buffered: the m x s accumulator is touched once every NBUF steps
    nb = nb + 1;
    Zb(:, (nb-1)*rho+1:nb*rho) = Zn;
    Cb{nb} = C;
    if nb == nbuf || st == ns
        Cm = zeros(s, nb*rho);
        if cplx, Cm = complex(Cm); end
        for j = 1:nb, Cm(1:size(Cb{j}, 1), (j-1)*rho+1:j*rho) = Cb{j}; end
        Zacc(:, 1:s) = Zacc(:, 1:s) + Zb(:, 1:nb*rho)*Cm';
        nb = 0;
        if solmode, nXest = norm(Zacc(:, 1:s), 'fro'); end
    end
    info.res(st) = sqrt(abs(trace((Zt'*Zt)*(Yt'*Yt))))/nF;
end
% ---- no truncation (s is already close to the rank); report the numerical rank
t0 = tic;
Zc = Zacc(:, 1:s);  Yc = Qy(:, 1:s);                        % X ~ Zc*Yc'
G = Zc'*Zc;
sv = sqrt(max(sort(eig((G + G')/2), 'descend'), 0));
info.rank = find(sv > 1e-8*sv(1), 1, 'last');               % numerical rank at 1e-8
info.time = toc(tstart);
info.tcomp = toc(t0);
info.s = s;
info.maxcols = s + nbuf*rho;
info.dropped = dropped/max(norm(Zc, 'fro'), realmin);       % relative to ||X||_F
end
