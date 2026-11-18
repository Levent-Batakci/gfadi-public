function [Xo, info] = tb_adv_arnoldi(op, U, L, V, outfun, tol, opts)
%TB_ADV_ARNOLDI  One-pass block Arnoldi Galerkin for  A X + X L' = U V'  with a nonsymmetric A and restricted outputs.
%   [XO, INFO] = TB_ADV_ARNOLDI(OP, U, L, V, OUTFUN, TOL, OPTS): the Krylov competitor when Lanczos's
%   short recurrence is not available.  OP.mv applies A to [OP.sz rho] arrays, OUTFUN maps a basis block
%   to its outputs.  The basis Q (m x k rho) is kept (every new block is orthogonalized against all
%   earlier ones), and so are the outputs of the basis: no second pass.  Stops at the exact Galerkin
%   residual TOL (geometric schedule, then bisection, as TB_MEGA_LANCZOS).
%   OPTS: .tcap [Inf s] .kcap [256] (expected iterations, preallocates) .kmax [100000]
%   XO: outputs x n.  INFO: .k .products .time .tproduct .torth .tsmall .res .capped .basisGB
%   Method 'arnoldi' of SWITCH_RUN (Example 3 with the swirl, Fig. 9).
if nargin < 7, opts = struct(); end
tcap = Inf;
if isfield(opts, 'tcap'), tcap = opts.tcap; end
kmax = 100000;
if isfield(opts, 'kmax'), kmax = opts.kmax; end
cap = 256;
if isfield(opts, 'kcap'), cap = opts.kcap; end
sz = size(U);
ndg = numel(op.sz);
if numel(sz) == ndg, rho = 1; else, rho = sz(ndg + 1); end
m = prod(op.sz);
n = size(V, 1);
[Q1, R0] = qr(reshape(double(U), m, rho), 0);
Lt = L.';
Qs = [];
if ~istriu(Lt), [Qs, Lt] = schur(full(Lt), 'complex'); end                     % a general small side: its complex Schur form
C0 = R0*V.';
if ~isempty(Qs), C0 = C0*Qs; end
nF = norm(R0*(V'*V)^(1/2), 'fro');
% the basis in chunks of cb blocks: no copy of the basis when it grows or is used (a slice Q(:, 1:K) would copy it every step)
cb = 16;
Qc = {zeros(m, cb*rho)};
Qc{1}(:, 1:rho) = Q1;
c1 = outfun(reshape(Q1, [op.sz rho]));
no = size(c1, 1);
Cq = zeros(no, (cap + 1)*rho);
Cq(:, 1:rho) = c1;
H = zeros((cap + 1)*rho, cap*rho);
t0 = tic;  tp = 0;  to = 0;  ts = 0;  k = 0;  next = 20;  klo = 0;  done = false;  res = NaN;
while k < kmax && toc(t0) < tcap
    K = k*rho;
    J = K + (1:rho);
    [ci, cj] = loc(K, rho, cb);
    a = tic;  W = reshape(op.mv(reshape(Qc{ci}(:, cj), [op.sz rho])), m, rho);  tp = tp + toc(a);
    a = tic;
    Kb = K + rho;
    h = zeros(Kb, rho);
    for pass = 1:2                                                % block classical Gram-Schmidt, twice, chunk by chunk
        for c = 1:numel(Qc)
            hc = Qc{c}'*W;
            W = W - Qc{c}*hc;
            r0 = (c - 1)*cb*rho;
            nr = min(cb*rho, Kb - r0);
            if nr > 0, h(r0 + (1:nr), :) = h(r0 + (1:nr), :) + hc(1:nr, :); end
        end
    end
    [Qn, Rn] = qr(W, 0);
    clear W
    k = k + 1;
    if (k + 1)*rho > size(Cq, 2), nc = ceil(1.5*(k + 1))*rho;  Cq(:, nc) = 0;  H(nc, nc - rho) = 0; end
    [ci, cj] = loc(Kb, rho, cb);
    if ci > numel(Qc), Qc{ci} = zeros(m, cb*rho); end
    H(1:Kb, J) = h;
    H(Kb + (1:rho), J) = Rn;
    Qc{ci}(:, cj) = Qn;
    to = to + toc(a);
    a = tic;  Cq(:, Kb + (1:rho)) = outfun(reshape(Qn, [op.sz rho]));  to = to + toc(a);
    if k == next                                                  % residual check, geometric schedule
        a = tic;  [~, res] = small_solve(H, k, rho, C0, Lt, n);  res = res/nF;  ts = ts + toc(a);
        if res <= tol, done = true;  break; end
        klo = k;
        next = ceil(1.25*k);
    end
end
if done
    % bisection for the smallest k (within 2%) that meets tol
    khi = k;
    while khi - klo > max(2, round(0.02*khi))
        mid = floor((klo + khi)/2);
        a = tic;  [~, rm] = small_solve(H, mid, rho, C0, Lt, n);  rm = rm/nF;  ts = ts + toc(a);
        if rm <= tol, khi = mid;  res = rm; else, klo = mid; end
    end
    k = khi;
end
a = tic;  [Y, r2] = small_solve(H, k, rho, C0, Lt, n);  if ~done, res = r2/nF; end;  ts = ts + toc(a);
if ~isempty(Qs), Y = Y*Qs'; end
Xo = real(Cq(:, 1:k*rho)*Y);
info = struct('k', k, 'products', k*rho, 'time', toc(t0), 'tproduct', tp, 'torth', to, 'tsmall', ts, 'res', res, 'capped', ~done, 'basisGB', m*numel(Qc)*cb*rho*8/2^30);
end

function [ci, cj] = loc(K, rho, cb)
% chunk and columns of the block that starts after column K
ci = floor(K/(cb*rho)) + 1;
c0 = K - (ci - 1)*cb*rho;
cj = c0 + (1:rho);
end

function [Y, r] = small_solve(H, k, rho, C, Lt, n)
% H_k Y + Y Lt = [C; 0]: complex Schur of H_k, then column substitution in the upper triangular Lt
K = k*rho;
Hk = H(1:K, 1:K);
[Z, T] = schur(Hk, 'complex');
B = zeros(K, n);
B(1:rho, :) = C;
B = Z'*B;
Yt = complex(zeros(K, n));
dL = full(diag(Lt));
I = eye(K);
for j = 1:n
    rhs = B(:, j);
    [ci, ~, cv] = find(Lt(1:j-1, j));
    if ~isempty(ci), rhs = rhs - Yt(:, ci)*cv; end
    Yt(:, j) = (T + dL(j)*I)\rhs;
end
Y = Z*Yt;
r = norm(H(K + (1:rho), K - rho + (1:rho))*Y(K - rho + 1:K, :), 'fro');      % ||H_{k+1,k} E_k' Y||_F
end
