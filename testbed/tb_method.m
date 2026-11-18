function r = tb_method(name, P, tol, opts)
%TB_METHOD  Run one solver on an Example 1 problem  A X + X L^T = U V^T  (see TB_PROBLEM).
%   R = TB_METHOD(NAME, P, TOL, OPTS).  TOL is g-fADI's a priori target.  NAME:
%     'gfadi'    one-sided g-fADI (right steps only), Chebyshev shifts from P.E and Re W(L)
%     'mixed'    g-fADI of type (p,q), q left steps chosen from opts.qlist by measured unit costs
%     'fadi'     classical factored ADI, type (k,k), optimal shifts
%     'bs'       Bartels-Stewart with a sparse large side (n shifted solves with A)
%     'eksm'     extended Krylov, one factorization of A, Galerkin
%     'lanczos2' two-pass block Lanczos Galerkin (SPD A only)
%     'arnoldi'  one-pass block Arnoldi/Lanczos Galerkin, full orthogonalization
%   OPTS: .restarget [TOL] residual at which Krylov methods stop, .tmax [120 s] (beyond it 'fadi'
%   and 'bs' are projected), .Xref [] reference for .err, .qlist [2 4 8], .fillGB [NaN], .memcap [6 GB],
%   .nores [false] skip the exact residual (and .err), so the caller can time and measure without it,
%   .trunc [1e-12] relative SVD truncation of the Krylov results
%   R: .method .time .tdesign .its .res (exact relative residual) .err .memGB (estimated peak)
%      .projected .note, and .Z .Y (or .Xfull for 'bs'), which the caller removes after use.
%   Called by TB_SWEEP (Figs. 6 and 7(c)) and TEASER_RUN (Fig. 1).
if nargin < 4, opts = struct(); end
o = struct('restarget', tol, 'tmax', 120, 'Xref', [], 'qlist', [2 4 8], 'fillGB', NaN, 'memcap', 6, 'nores', false, 'trunc', 1e-12);
fn = fieldnames(opts);
for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
A = P.A;
L = P.L;
U = P.U;
V = P.V;
m = P.m;
n = P.n;
rho = size(U, 2);
r = struct('method', name, 'time', NaN, 'tdesign', 0, 'its', NaN, 'res', NaN, 'err', NaN, ...
           'memGB', NaN, 'projected', false, 'note', '', 'Z', [], 'Y', []);
Eint = P.E;
F = [-P.WL(2), -P.WL(1)];
GB = @(cols) cols*m*8/2^30;
symA = issymmetric(A);
try
switch name
    case 'gfadi'
        t = tic;  [p, al] = shift_search_real(Eint, F, 0, tol);
        if p > 2e5, r.note = sprintf('p = %d too large', p); r.its = p; return; end
        al = leja_order(al, Eint');  r.tdesign = toc(t);
        t = tic;  Bop = small_op(L);
        [Z, Y, info] = gfadi_lowrank(gfadi_op(A), Bop, U, V, al, []);
        r.time = toc(t);
        r.its = p;
        r.memGB = GB(info.maxcols + 4*rho);
        r.Z = Z;
        r.Y = Y;
        r.note = sprintf('p=%d s=%d', p, info.s);
    case 'mixed'
        % measure the unit costs: one right step (product with A + small solve), one left step
        Bop = small_op(L);
        Aop = gfadi_op(A);
        Zt = U;
        Yt = V;
        sa = mean(Eint);
        t = tic;  for k = 1:5, Zt = A*Zt - sa*Zt;  Yt = Bop.svh(sa, Yt); end;  CR = toc(t)/5;
        t = tic;  Ztmp = (A - F(2)*speye(m))\U;  CL = toc(t);  clear Ztmp       %#ok<NASGU>
        t = tic;  best = Inf;
        [p0, al0] = shift_search_real(Eint, F, 0, tol);
        cand = {0, p0, al0, []};
        best = p0*CR;
        if p0 > 3000, r.note = sprintf('skipped: p0 = %d (type-(p,q) design too costly)', p0); return; end
        for q = o.qlist
            if q*CL > best || toc(t) > 90, break; end
            [pq, alq, beq] = shift_search_real(Eint, F, q, tol, round(0.6*p0));
            cost = pq*CR + q*CL;
            if cost < best, best = cost; cand = {q, pq, alq, beq}; end
        end
        r.tdesign = toc(t);
        [q, p, al, be] = cand{:};
        if q == 0, al = leja_order(al, Eint'); sc = '';
        else, [al, be, sc] = order_shifts(al, be, Eint, F); end
        t = tic;
        [Z, Y, info] = gfadi_lowrank(Aop, Bop, U, V, al, be, 'sched', sc);
        r.time = toc(t);
        r.its = p + q;
        r.memGB = GB(info.maxcols + 4*rho) + o.fillGB*(q > 0);
        r.Z = Z;
        r.Y = Y;
        r.note = sprintf('(p,q)=(%d,%d) c=%.0f', p, q, CL/CR);
    case 'fadi'
        t = tic;  [k, al, be] = classical_shifts_real(Eint, F, tol);  [al, be, sc] = order_shifts(al, be, Eint, F);  r.tdesign = toc(t);
        % project if k shifted solves would exceed tmax
        t = tic;  Ztmp = (A - be(1)*speye(m))\U;  t1 = toc(t);  clear Ztmp     %#ok<NASGU>
        if k*t1 > o.tmax
            r.time = k*t1;
            r.projected = true;
            r.its = 2*k;
            r.note = sprintf('k=%d projected', k);
            r.memGB = GB(2*k*rho + 4*rho) + o.fillGB;
            return
        end
        t = tic;
        [Z, Y, info] = gfadi_lowrank(gfadi_op(A), small_op(L), U, V, al, be, 'sched', sc);
        r.time = toc(t);
        r.its = 2*k;
        r.memGB = GB(info.maxcols + 4*rho) + o.fillGB;
        r.Z = Z;
        r.Y = Y;
        r.note = sprintf('k=%d', k);
    case 'bs'
        t = tic;
        Lt = full(L.');
        [Q, T] = schur(Lt, 'real');
        cplx = any(abs(diag(T, -1)) > 0);
        if cplx, [Q, T] = schur(Lt, 'complex'); end
        Ft = U*(V.'*Q);
        if GB(n) > o.memcap, r.note = 'X does not fit'; r.memGB = GB(n); return; end
        Xt = zeros(m, n);
        if cplx, Xt = complex(Xt); end
        I = speye(m);
        done = n;
        for j = 1:n
            Xt(:, j) = (A + T(j, j)*I)\(Ft(:, j) - Xt(:, 1:j-1)*T(1:j-1, j));
            if j >= 3 && toc(t)*n/j > o.tmax, done = j; break; end      % projected from the first solves
        end
        if done < n
            r.time = toc(t)*n/done;
            r.projected = true;
            r.note = sprintf('projected from %d of %d solves', done, n);
        else
            X = Xt*Q';
            if ~cplx, X = real(X); end
            r.time = toc(t);
            r.Xfull = real(X);
        end
        r.its = n;
        r.memGB = GB(n) + o.fillGB;
    case 'eksm'
        [Z, Y, info] = tb_eksm(A, L, U, V, o.restarget, symA, o);
        r.time = info.time;
        r.its = info.k;
        r.memGB = info.mem_GB;
        r.Z = Z;
        r.Y = Y;
        r.projected = info.capped;
        r.note = sprintf('k=%d basis=%d fact=%.2fs', info.k, info.basis, info.tfact);
        r.fillGB = info.fill_GB;
    case 'lanczos2'
        if ~symA, r.note = 'needs SPD A'; return; end
        [Z, Y, info] = lanczos2_sylv(A, L, U, V, o.restarget);
        r.time = info.time;
        r.its = info.k;
        r.memGB = GB(info.rank + 3*rho);
        r.Z = Z;
        r.Y = Y;
        r.note = sprintf('k=%d (%d products)', info.k, info.nmv);
    case 'arnoldi'
        [Z, Y, info] = tb_arnoldi(A, L, U, V, o.restarget, symA, o);
        r.time = info.time;
        r.its = info.k;
        r.memGB = GB(info.basis + rho);
        r.Z = Z;
        r.Y = Y;
        r.projected = info.capped;
        r.note = sprintf('k=%d', info.k);
end
catch err
    r.note = ['ERROR: ' err.message];
    return
end
if o.nores, return; end
if ~isempty(r.Z)
    r.res = tb_res(A, L, U, V, r.Z, r.Y);
    if ~isempty(o.Xref), r.err = norm(r.Z*r.Y' - o.Xref, 'fro')/norm(o.Xref, 'fro'); end
end
if isfield(r, 'Xfull'), r.res = tb_res(A, L, U, V, r.Xfull, eye(n)); end
end

function op = small_op(L)
% the small side of g-fADI: B = -L^T, shifted solves with B (right steps)
n = size(L, 1);
if isdiag(L), op = gfadi_op(-full(diag(L)), 'diag'); return; end
if issparse(L), op = gfadi_op(-L.'); return; end
% dense: Hessenberg form once, then O(n^2) per shifted solve (backslash detects Hessenberg)
[Qh, H] = hess(-L.');
I = eye(n);
Ht = H';
op.mv  = @(X) -(L.'*X);
op.mvh = @(X) -(L*X);
op.sv  = @(s, X) Qh*((H - s*I)\(Qh'*X));
op.svh = @(s, X) Qh*((Ht - conj(s)*I)\(Qh'*X));
op.n = n;
end

function res = tb_res(A, L, U, V, Z, Y)
% exact || A Z Y' + Z Y' L' - U V' ||_F / || U V' ||_F from thin QR factors
[~, R1] = qr([A*Z, Z, -U], 0);
[~, R2] = qr([Y, L*Y, V], 0);
res = norm(R1*R2', 'fro')/sqrt(abs(trace((U'*U)*(V'*V))));
end

function [Z, Vy, info] = tb_eksm(A, L, U, V, tol, symA, o)
% extended Krylov (EKSM) with any A: CHOLMOD if SPD, else sparse LU (decomposition)
t0 = tic;
m = size(A, 1);
rho = size(U, 2);
nF = sqrt(abs(trace((U'*U)*(V'*V))));
kmax = 400;
every = 2;
capped = false;
tf = tic;
if symA, [Lc, ~, pp] = chol(A, 'lower', 'vector');  fill = nnz(Lc)*16/2^30;  Ainv = @(X) cholsolve(Lc, pp, X);
else, dA = decomposition(A, 'lu');  Ainv = @(X) dA\X;  fill = NaN; end
tfact = toc(tf);
if symA && issparse(L), Qs = [];  Ts = []; else, [Qs, Ts] = schur(full(L.'), 'complex'); end
W0 = [U, Ainv(U)];
cap = 2*rho*64;
Q = zeros(m, cap);
AQ = zeros(m, cap);
s = 0;
[Qb, ~] = qr(W0, 0);
[Q, AQ, s] = app(Q, AQ, s, Qb, A);
res = NaN;
Y = [];
for k = 1:kmax
    if mod(k, every) == 0 || k == kmax
        H = Q(:, 1:s)'*AQ(:, 1:s);
        if symA, H = (H + H')/2; end
        Y = proj_solve(H, Qs, Ts, (Q(:, 1:s)'*U)*V.', symA, L);
        [~, Rd] = qr(AQ(:, 1:s) - Q(:, 1:s)*H, 0);
        res = norm(Rd*Y, 'fro')/nF;
        if res <= tol, break; end
        if (s + 2*rho)*m*8*2/2^30 > o.memcap || toc(t0) > 4*o.tmax, capped = true; break; end
    end
    last = Q(:, s-2*rho+1:s);
    Wn = [A*last(:, 1:rho), Ainv(last(:, rho+1:end))];
    for pass = 1:2, Wn = Wn - Q(:, 1:s)*(Q(:, 1:s)'*Wn); end
    [Qb, ~] = qr(Wn, 0);
    [Q, AQ, s] = app(Q, AQ, s, Qb, A);
end
[Uy, Sy, Vy] = svd(Y, 'econ');
sv = diag(Sy);
r = find(sv > o.trunc*sv(1), 1, 'last');   % relative SVD truncation of the result (opts.trunc, default 1e-12)
Z = Q(:, 1:size(Y, 1))*(Uy(:, 1:r)*Sy(1:r, 1:r));
Vy = Vy(:, 1:r);
info = struct('k', k, 'basis', s, 'time', toc(t0), 'tfact', tfact, 'res', res, 'fill_GB', fill, ...
              'mem_GB', 2*s*m*8/2^30 + max(fill, 0), 'capped', capped);
end

function X = cholsolve(Lc, pp, B)
X = zeros(size(B));
X(pp, :) = Lc'\(Lc\B(pp, :));
end

function [Q, AQ, s] = app(Q, AQ, s, Qb, A)
% append the block Qb to the basis Q (and A*Qb to AQ), doubling the storage when full
t = size(Qb, 2);
if s + t > size(Q, 2), Q(:, 2*(s + t)) = 0;  AQ(:, 2*(s + t)) = 0; end
Q(:, s+1:s+t) = Qb;
AQ(:, s+1:s+t) = A*Qb;
s = s + t;
end

function Y = proj_solve(H, Qs, Ts, C, symH, L)
% H Y + Y L^T = C.  Symmetric H and sparse L: H = W diag(lam) W', one sparse solve with L per
% eigenvalue, O(s n) for banded L.  Otherwise L^T = Qs Ts Qs' (complex Schur, precomputed).
n = size(L, 1);
if symH && issparse(L)
    [W, lam] = eig((H + H')/2, 'vector');
    Ch = W'*C;
    Id = speye(n);
    Yh = zeros(size(Ch));
    for i = 1:numel(lam)
        Yh(i, :) = ((L + lam(i)*Id)\Ch(i, :).').';
    end
    Y = W*Yh;
elseif symH
    [W, lam] = eig((H + H')/2, 'vector');
    Ct = W'*C*Qs;
    Yt = complex(zeros(size(Ct)));
    for j = 1:n
        Yt(:, j) = (Ct(:, j) - Yt(:, 1:j-1)*Ts(1:j-1, j))./(lam + Ts(j, j));
    end
    Y = real(W*Yt*Qs');
else
    Y = real(sylvester(H, Qs*Ts*Qs', C));
end
end

function [Z, Vy, info] = tb_arnoldi(A, L, U, V, tol, symA, o)
% one-pass block Arnoldi (Lanczos with full reorthogonalisation if A is SPD) Galerkin:
% X ~ Q Y,  H Y + Y L^T = E1 R0 V^T,  residual || H(k+1,k) Y_k ||
t0 = tic;
m = size(A, 1);
rho = size(U, 2);
nF = sqrt(abs(trace((U'*U)*(V'*V))));
if symA && issparse(L), Qs = [];  Ts = []; else, [Qs, Ts] = schur(full(L.'), 'complex'); end
kmax = floor(o.memcap*2^30/(8*m*rho)) - 1;
capped = false;
cap = rho*64;
Q = zeros(m, cap);
H = zeros(cap + rho, cap);
[Q1, R0] = qr(U, 0);
Q(:, 1:rho) = Q1;
s = rho;
C0 = R0*V.';
next = 10;
Y = [];
res = NaN;
for k = 1:kmax
    blk = s-rho+1:s;
    W = A*Q(:, blk);
    for pass = 1:2, hc = Q(:, 1:s)'*W;  W = W - Q(:, 1:s)*hc;  H(1:s, blk) = H(1:s, blk) + hc; end
    [Qn, Hn] = qr(W, 0);
    if s + rho > size(Q, 2), Q(:, 2*(s + rho)) = 0;  H(2*(s + rho) + rho, 2*(s + rho)) = 0; end
    H(s+1:s+rho, blk) = Hn;
    Q(:, s+1:s+rho) = Qn;
    if k == next || k == kmax
        Hk = H(1:s, 1:s);
        C = zeros(s, size(V, 1));
        C(1:rho, :) = C0;
        Y = proj_solve(Hk, Qs, Ts, C, symA, L);
        res = norm(Hn*Y(blk, :), 'fro')/nF;
        if res <= tol, break; end
        if toc(t0) > 2*o.tmax, capped = true; break; end
        next = max(k + 2, ceil(1.15*k));
    end
    s = s + rho;
end
if k == kmax && res > tol, capped = true; end
sy = size(Y, 1);
[Uy, Sy, Vy] = svd(Y, 'econ');
sv = diag(Sy);
r = find(sv > o.trunc*sv(1), 1, 'last');   % relative SVD truncation of the result (opts.trunc, default 1e-12)
Z = Q(:, 1:sy)*(Uy(:, 1:r)*Sy(1:r, 1:r));
Vy = Vy(:, 1:r);
info = struct('k', k, 'basis', s, 'time', toc(t0), 'res', res, 'capped', capped);
end
