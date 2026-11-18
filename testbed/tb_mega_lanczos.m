function [Xo, info] = tb_mega_lanczos(op, U, L, V, outfun, tol, opts)
%TB_MEGA_LANCZOS  One-pass block Lanczos Galerkin for  A X + X L' = U V'  with restricted outputs.
%   [XO, INFO] = TB_MEGA_LANCZOS(OP, U, L, V, OUTFUN, TOL, OPTS).  OP.mv applies A to [OP.sz rho]
%   arrays, OUTFUN maps a basis block to its outputs.  Keeps the three-term blocks and the outputs of
%   the basis, never the basis itself: no second pass and no basis in memory.  Stops at the exact
%   Galerkin residual TOL (geometric schedule, then bisection).  A must be symmetric; for a
%   nonsymmetric A use TB_ADV_ARNOLDI.
%   OPTS: .kmax [200000] .tcap [Inf s] .kcap [512] (expected iterations, preallocates the outputs)
%   XO: outputs x n.  INFO: .k .products .time .tproduct .tvec .tsmall .tout .res .capped
%   Method 'lanczos' of SWITCH_RUN (Example 3, symmetric case SWITCH_DRIVE(0)).
if nargin < 7, opts = struct(); end
kmax = 200000;
if isfield(opts, 'kmax'), kmax = opts.kmax; end
tcap = Inf;
if isfield(opts, 'tcap'), tcap = opts.tcap; end
sz = size(U);
ndg = numel(op.sz);
if numel(sz) == ndg, rho = 1; else, rho = sz(ndg + 1); end
m = prod(op.sz);
n = size(V, 1);
[Qc, R0] = cholqr2(reshape(U, m, rho));
Lt = L.';
isdiagL = isdiag(Lt);
isup = istriu(Lt);
Qs = [];
if ~(isdiagL || isup), [Qs, Lt] = schur(full(Lt), 'complex'); end               % only if L' is not triangular
C0 = R0*V.';
if ~isempty(Qs), C0 = C0*Qs; end
nF = norm(R0*(V'*V)^(1/2), 'fro');
c1 = outfun(reshape(Qc, [op.sz rho]));
no = size(c1, 1);
cap = 512;
if isfield(opts, 'kcap'), cap = opts.kcap; end      % so that a large output store is not copied as it grows
Cq = zeros(no, cap*rho, class(c1));
Cq(:, 1:rho) = c1;
Ad = zeros(rho, rho, 256);                                           % diagonal blocks of T_k
Bd = zeros(rho, rho, 256);                                           % subdiagonal blocks
Qp = zeros(m, rho, class(Qc));
Bprev = zeros(rho);
t0 = tic;  tprod = 0;  tvec = 0;  tsm = 0;  tout = 0;  res = NaN;  k = 0;  next = 40;  klo = 0;  done = false;
while k < kmax && toc(t0) < tcap
    tp = tic;
    W = op.mv(reshape(Qc, [op.sz rho]));
    W = reshape(W, m, rho);
    tprod = tprod + toc(tp);
    tv = tic;
    W = W - Qp*cast(Bprev', class(Qc));
    Aj = Qc'*W;
    Aj = (Aj + Aj')/2;
    W = W - Qc*Aj;
    [Qn, Bj] = cholqr2(W);
    clear W;
    Aj = double(Aj);
    Bj = double(Bj);
    tvec = tvec + toc(tv);
    k = k + 1;
    if k > size(Ad, 3), Ad(:, :, 2*k) = 0;  Bd(:, :, 2*k) = 0; end
    Ad(:, :, k) = Aj;
    Bd(:, :, k) = Bj;
    Qp = Qc;
    Qc = Qn;
    Bprev = Bj;
    clear Qn
    to = tic;
    if (k + 1)*rho > size(Cq, 2), Cq(:, ceil(1.5*(k + 1))*rho) = 0; end
    Cq(:, k*rho + (1:rho)) = outfun(reshape(Qc, [op.sz rho]));
    tout = tout + toc(to);
    if k == next                                                     % residual check, geometric schedule
        ts = tic;  [~, res] = small_solve(Ad, Bd, k, C0, Lt, rho, n);  res = res/nF;  tsm = tsm + toc(ts);
        if res <= tol, done = true; break; end
        klo = k;
        next = ceil(1.25*k);
    end
end
if done
    % bisection for the smallest k (within 2%) that meets tol
    khi = k;
    while khi - klo > max(2, round(0.02*khi))
        mid = floor((klo + khi)/2);
        ts = tic;  [~, rm] = small_solve(Ad, Bd, mid, C0, Lt, rho, n);  rm = rm/nF;  tsm = tsm + toc(ts);
        if rm <= tol, khi = mid; res = rm; else, klo = mid; end
    end
    k = khi;
end
ts = tic;  [Y, res2] = small_solve(Ad, Bd, k, C0, Lt, rho, n);  if ~done, res = res2/nF; end;  tsm = tsm + toc(ts);
if ~isempty(Qs), Y = Y*Qs'; end
Xo = Cq(:, 1:k*rho)*Y;
if isreal(U) && isreal(L) && isreal(V), Xo = real(Xo); end
info = struct('k', k, 'products', k*rho, 'time', toc(t0), 'tproduct', tprod, 'tvec', tvec, 'tsmall', tsm, 'tout', tout, 'res', res, 'capped', ~done);
end

function [Y, r] = small_solve(Ad, Bd, k, C, Lt, rho, n)
% T_k Y + Y Lt = [C; 0], Lt upper triangular (or diagonal), by columns; r = ||B_k Y_k||_F
K = rho*k;
[I, J, S] = deal(cell(2*k, 1));
for j = 1:k
    [ii, jj] = ndgrid((j - 1)*rho + (1:rho));
    I{j} = ii(:);
    J{j} = jj(:);
    S{j} = reshape(Ad(:, :, j), [], 1);
end
c = k;
for j = 1:k - 1
    [ii, jj] = ndgrid(j*rho + (1:rho), (j - 1)*rho + (1:rho));
    B = Bd(:, :, j);
    c = c + 1;
    I{c} = [ii(:); jj(:)];
    J{c} = [jj(:); ii(:)];
    S{c} = [B(:); B(:)];
end
T = sparse(vertcat(I{:}), vertcat(J{:}), vertcat(S{:}), K, K);      % block tridiagonal T_k
Id = speye(K);
cplx = ~(isreal(C) && isreal(Lt));
Y = zeros(K, n);
if cplx, Y = complex(Y); end
dLt = full(diag(Lt));
sparseL = issparse(Lt);
for j = 1:n
    rhs = zeros(K, 1);
    rhs(1:rho) = C(:, j);
    if sparseL
        % (find on the column, not Lt(1:j-1, j): that indexing is O(n) per column)
        [ci, ~, cv] = find(Lt(:, j));
        kp = ci < j;
        ci = ci(kp);
        cv = cv(kp);
        if ~isempty(ci), rhs = rhs - Y(:, ci)*cv; end
    elseif j > 1
        rhs = rhs - Y(:, 1:j-1)*Lt(1:j-1, j);
    end
    Y(:, j) = (T + dLt(j)*Id)\rhs;
end
r = norm(Bd(:, :, k)*Y(end-rho+1:end, :), 'fro');
end

function [Q, R] = cholqr2(W)
% CholeskyQR2; rho = 1 is a plain normalisation
if size(W, 2) == 1, R = norm(W);  Q = W/R;  return; end
G = W'*W;
R1 = chol((G + G')/2);
Q = W/R1;
G = Q'*Q;
R2 = chol((G + G')/2);
Q = Q/R2;
R = R2*R1;
end
