function [Zo, Qy, info] = gfadi_mega(op, U, Bop, V, alpha, outfun, varargin)
%GFADI_MEGA  One-sided g-fADI on very large grids, keeping only chosen outputs of X.
%   [ZO,QY,INFO] = GFADI_MEGA(OP,U,BOP,V,ALPHA,OUTFUN,...) runs right steps only for
%   A X - X B = U V' and never forms X.  The outputs are  X_out = Zo*Qy'  (no x n),
%   Qy an orthonormal basis of the small side (as in GFADI_LOWRANK).  This is
%   Algorithm 4.1 of the paper with Gamma given by OUTFUN (sampled output); used by
%   Example 3 (SWITCH_RUN).
%   OP:     large side, op.mv(x, a, g) = g*(A - a I) x on grid arrays [op.sz rho] (STENCIL3)
%   U:      [op.sz rho] array.   BOP: GFADI_OP of the small side B (n x n).   V: n x rho
%   ALPHA:  the shifts, already ordered.
%   OUTFUN: maps a grid array Z to an (no x rho) matrix of outputs (e.g. values at chosen
%           cells, a coarse volume, projections).
%   Options:  'ytol' (1e-13) basis tolerance relative to the outputs so far
%             'nbuf' (32) output updates buffered over nbuf steps
%             'full' n-by-w W: also accumulate XW = X*W on the full grid ([] = none)
%             'verbose' print progress every so many seconds (300, 0 = silent)
%   INFO: .time .tproduct .tsmall .toutput .s .p .res (estimate every 50 steps) .dropped .XW
q = inputParser;
q.addParameter('ytol', 1e-13);
q.addParameter('nbuf', 32);
q.addParameter('full', []);
q.addParameter('verbose', 300);
q.parse(varargin{:});
o = q.Results;
sz = size(U);
ndg = numel(op.sz);
if numel(sz) == ndg, rho = 1; else, rho = sz(ndg + 1); end
n = size(V, 1);  p = numel(alpha);  alpha = alpha(:);
nF2 = norm(U(:))^2*norm(V, 'fro')^2/rho;                       % a scale for the residual estimate (equal for rho = 1)
Zt = U;  sc = 1;  Yt = V;                                      % the residual factor is sc*Zt, Yt
s = 0;  cap = min(n, 256);
Qy = zeros(n, cap);
c0 = outfun(Zt);
no = size(c0, 1);  cls = class(c0);
Zo = zeros(no, cap, cls);
Zb = zeros(no, o.nbuf*rho, cls);
Cb = cell(1, o.nbuf);  nb = 0;
nXest = 0;  dropped = 0;
W = o.full;
if ~isempty(W), XW = zeros([op.sz size(W, 2)], class(U)); else, XW = []; end
tprod = 0;  tsmall = 0;  tout = 0;
t0 = tic;  tlast = t0;
res = nan(1, ceil(p/50));
for k = 1:p
    a = alpha(k);
    % small side
    ts = tic;
    Yn = Bop.svh(a, Yt);
    tsmall = tsmall + toc(ts);
    % outputs of this increment: Zn = -sc*Zt
    to = tic;
    Co = -sc*outfun(Zt);
    % compress Yn in the basis of the small side (as in GFADI_LOWRANK, 'solution' mode)
    Qs = Qy(:, 1:s);
    C = Qs'*Yn;  Yp = Yn - Qs*C;
    C2 = Qs'*Yp;  Yp = Yp - Qs*C2;
    C = C + C2;
    nyn = norm(Yn, 'fro');  nzn = norm(Co, 'fro');  nyp = norm(Yp, 'fro');
    nXest = max(nXest, nzn*nyn);
    thr = o.ytol*nXest/max(nzn, realmin);
    if nyp > max(o.ytol*nyn, thr) && s < n
        [Qn, Sn] = svd(Yp, 'econ');
        Qn = Qn(:, diag(Sn) > max(o.ytol*nyn, thr));
        for pass = 1:2, Qn = Qn - Qs*(Qs'*Qn);  [Qn, ~] = qr(Qn, 0); end
        t = min(size(Qn, 2), n - s);
        Qn = Qn(:, 1:t);
        if t > 0
            if s + t > cap, nc = min(n, max(2*cap, s + t));  Qy(:, nc) = 0;  Zo(:, nc) = 0;  cap = nc; end
            Cn = Qn'*Yp;  Yp = Yp - Qn*Cn;
            Qy(:, s+1:s+t) = Qn;
            C = [C; Cn];   %#ok<AGROW>
            s = s + t;
        end
    end
    dropped = dropped + nzn*norm(Yp, 'fro');
    nb = nb + 1;
    Zb(:, (nb - 1)*rho + 1:nb*rho) = Co;
    Cb{nb} = C;
    if nb == o.nbuf || k == p
        Cm = zeros(s, nb*rho, cls);
        if ~isreal(Qy), Cm = complex(Cm); end
        for j = 1:nb, Cm(1:size(Cb{j}, 1), (j - 1)*rho + 1:j*rho) = Cb{j}; end
        % update by chunks of columns, so the whole factor is never copied
        nz2 = 0;
        for j0 = 1:64:s
            J = j0:min(j0 + 63, s);
            if nb == o.nbuf, T = Zo(:, J) + Zb*Cm(J, :)'; else, T = Zo(:, J) + Zb(:, 1:nb*rho)*Cm(J, :)'; end
            Zo(:, J) = T;
            nz2 = nz2 + sum(abs(T).^2, 'all');
        end
        clear T;
        nb = 0;
        nXest = max(nXest, sqrt(nz2));
    end
    if ~isempty(W)                                              % X*W += Zn*(Yn'*W)
        G = cast(-sc*(Yn'*W), class(U));
        nw = size(W, 2);
        for w = 1:nw
            if rho == 1
                XW = idxadd(XW, ndg, w, G(1, w)*Zt);
            else
                tmp = reshape(Zt, prod(op.sz), rho)*G(:, w);
                XW = idxadd(XW, ndg, w, reshape(tmp, op.sz));
            end
        end
    end
    tout = tout + toc(to);
    % large side: the new residual factor sc*(A - a I) Zt (sign: -(A-a)Zn)
    tp = tic;
    Zn1 = op.mv(Zt, a, sc);
    tprod = tprod + toc(tp);
    cn = norm(Zn1(:));
    Zt = Zn1;
    clear Zn1;
    sc = 1/cn;
    Yt = Yn*cn;
    if mod(k, 50) == 0 || k == p
        if rho == 1
            r2 = norm(Yt)^2/cn^2*cn^2;                          % |Zt/cn| = 1 (rho = 1)
            r = sqrt(r2);
        else
            Zv = reshape(Zt, prod(op.sz), rho);
            r = sqrt(abs(trace((sc^2*(Zv'*Zv))*(Yt'*Yt))));
            clear Zv;
        end
        res(ceil(k/50)) = r/sqrt(nF2);
        if o.verbose > 0 && toc(tlast) > o.verbose
            fprintf('  gfadi_mega: step %d / %d, %.0f s, residual est. %.1e, basis %d\n', k, p, toc(t0), res(ceil(k/50)), s);
            tlast = tic;
        end
    end
end
Zo = Zo(:, 1:s);  Qy = Qy(:, 1:s);
info = struct('time', toc(t0), 'tproduct', tprod, 'tsmall', tsmall, 'toutput', tout, 's', s, 'p', p, 'res', res, ...
              'dropped', dropped/max(norm(Zo, 'fro'), realmin), 'XW', XW);
end

function XW = idxadd(XW, ndg, w, v)
% XW(..., w) = XW(..., w) + v
switch ndg
    case 3, XW(:, :, :, w) = XW(:, :, :, w) + v;
    case 1, XW(:, w) = XW(:, w) + v;
    case 2, XW(:, :, w) = XW(:, :, w) + v;
    otherwise, error('gfadi_mega:ndg', 'grids of 1, 2 or 3 dimensions');
end
end
