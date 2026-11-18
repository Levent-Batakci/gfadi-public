function [Z, Vy, info] = lanczos2_sylv(A, L, U, V, tol, opts)
%LANCZOS2_SYLV  Two-pass block-Lanczos Galerkin for A X + X L^T = U V^T, returning X = Z*Vy'.
%   [Z, Vy, info] = lanczos2_sylv(A, L, U, V, tol, opts)
%   A large SPD (only products A*Q), L small (kept exact).  Pass 1 runs block Lanczos and stops
%   when the exact Galerkin residual is <= tol*||U V^T||_F; pass 2 regenerates the basis to form
%   Z, so memory is m*(rank + 3*rho).
%   Projected equation T_k Y + Y L^T = [C; 0]: for sparse L, T_k (symmetric) is diagonalized and
%   each eigenvalue costs one sparse solve with L, O(k rho n) for banded L; for dense L, the complex
%   Schur form of L^T and column substitution, O(k rho n^2).
%   Method 'lanczos2' ("Lanczos") of the Example 1 test bed (TB_METHOD, Fig. 6).
%   opts: .kmax (20000)
%   info: .k (iterations), .nmv (block products, both passes), .time, .tsmall, .res, .rank
%   Example: [Z, Vy] = lanczos2_sylv(P.A, P.L, P.U, P.V, 1e-8);  X = Z*Vy';
if nargin < 6, opts = struct(); end
kmax = 20000;
if isfield(opts, 'kmax'), kmax = opts.kmax; end
t0 = tic;
m = size(A, 1);
rho = size(U, 2);
n = size(L, 1);
[Q1, R0] = qr(U, 0);
nF = sqrt(abs(trace((U'*U)*(V'*V))));
if issparse(L)
    Sm = struct('L', L, 'Qs', [], 'Ts', []);            % eigen-decomposition of T_k, sparse solves with L
    C = R0*V.';                                         % first block row of the rhs
else
    [Qs, Ts] = schur(full(L.'), 'complex');             % L^T = Qs Ts Qs'
    Sm = struct('L', [], 'Qs', Qs, 'Ts', Ts);
    C = R0*(V.'*Qs);                                    % first block row of the rhs, rotated
end
Ad = zeros(rho, rho, 0);                                % diagonal blocks of T_k
Bd = zeros(rho, rho, 0);                                % subdiagonal blocks of T_k
Qp = zeros(m, rho);
Qc = Q1;
Bprev = zeros(rho);
tsmall = 0;
res = NaN;
k = 0;
next = 40;                                              % next residual check (geometric schedule)
klo = 0;
done = false;
% pass 1: block Lanczos, keep only the small blocks
while k < kmax
    W = A*Qc - Qp*Bprev';
    Aj = Qc'*W;
    Aj = (Aj + Aj')/2;
    W = W - Qc*Aj;
    [Qn, Bj] = qr(W, 0);
    k = k + 1;
    Ad(:, :, k) = Aj;
    Bd(:, :, k) = Bj;
    Qp = Qc;
    Qc = Qn;
    Bprev = Bj;
    if k == next
        ts = tic;
        res = small_res(Ad, Bd, k, C, Sm, rho, n)/nF;
        tsmall = tsmall + toc(ts);
        if res <= tol, done = true; break; end
        klo = k;
        next = ceil(1.25*k);
    end
end
if done                                                 % bisect between klo and k
    khi = k;
    while khi - klo > max(2, round(0.02*khi))
        mid = floor((klo + khi)/2);
        ts = tic;
        rm = small_res(Ad, Bd, mid, C, Sm, rho, n)/nF;
        tsmall = tsmall + toc(ts);
        if rm <= tol, khi = mid; res = rm; else, klo = mid; end
    end
    k = khi;
end
ts = tic;
Yt = small_solve(Ad, Bd, k, C, Sm, rho, n);
if isempty(Sm.L), Y = real(Yt*Sm.Qs'); else, Y = real(Yt); end
tsmall = tsmall + toc(ts);
[Uy, Sy, Vy] = svd(Y, 'econ');
sv = diag(Sy);
r = find(sv > 1e-12*sv(1), 1, 'last');
G = Uy(:, 1:r)*Sy(1:r, 1:r);
Vy = Vy(:, 1:r);
% pass 2: regenerate the same basis, accumulate Z = Q*G
Z = Q1*G(1:rho, :);
Qp = zeros(m, rho);
Qc = Q1;
Bprev = zeros(rho);
for j = 1:k-1
    W = A*Qc - Qp*Bprev';
    W = W - Qc*Ad(:, :, j);
    [Qn, ~] = qr(W, 0);                                % the same operations as pass 1: bitwise the same basis
    Z = Z + Qn*G(j*rho+1:(j+1)*rho, :);
    Qp = Qc;
    Qc = Qn;
    Bprev = Bd(:, :, j);
end
info.k = k;
info.nmv = 2*k - 1;
info.time = toc(t0);
info.tsmall = tsmall;
info.res = res;
info.rank = r;
end

function T = small_T(Ad, Bd, k, rho)
% the block tridiagonal T_k as a sparse matrix
K = rho*k;
[I, J, S] = deal([]);
for j = 1:k
    [ii, jj] = ndgrid((j-1)*rho + (1:rho));
    I = [I; ii(:)];                                    %#ok<AGROW>
    J = [J; jj(:)];                                    %#ok<AGROW>
    S = [S; reshape(Ad(:, :, j), [], 1)];              %#ok<AGROW>
    if j < k
        [ii, jj] = ndgrid(j*rho + (1:rho), (j-1)*rho + (1:rho));
        B = Bd(:, :, j);
        I = [I; ii(:); jj(:)];                         %#ok<AGROW>
        J = [J; jj(:); ii(:)];                         %#ok<AGROW>
        S = [S; B(:); B(:)];                           %#ok<AGROW>
    end
end
T = sparse(I, J, S, K, K);
end

function [Yh, W] = small_solve_eig(Ad, Bd, k, C, L, rho)
% T_k = W diag(th) W': row i of Yh solves (L + th_i I) yh_i^T = (W(1:rho, i)' C)^T, so Y = W*Yh
T = full(small_T(Ad, Bd, k, rho));
[W, th] = eig((T + T')/2, 'vector');
Ch = W(1:rho, :)'*C;
n = size(L, 1);
Id = speye(n);
Yh = zeros(size(Ch), 'like', Ch);
for i = 1:numel(th)
    Yh(i, :) = ((L + th(i)*Id)\Ch(i, :).').';
end
end

function Y = small_solve(Ad, Bd, k, C, Sm, rho, n)
% T_k Y + Y L^T = [C; 0]; Y is in the Schur basis of L^T when L is dense
if ~isempty(Sm.L)
    [Yh, W] = small_solve_eig(Ad, Bd, k, C, Sm.L, rho);
    Y = W*Yh;
    return
end
Ts = Sm.Ts;
K = rho*k;
T = small_T(Ad, Bd, k, rho);
Id = speye(K);
Y = complex(zeros(K, n));
R = complex(zeros(K, n));
R(1:rho, :) = C;
for j = 1:n
    rhs = R(:, j) - Y(:, 1:j-1)*Ts(1:j-1, j);
    Y(:, j) = (T + Ts(j, j)*Id)\rhs;
end
end

function r = small_res(Ad, Bd, k, C, Sm, rho, n)
% exact Galerkin residual ||B_k Y_k||_F after k steps (last block row of Y only)
if ~isempty(Sm.L)
    [Yh, W] = small_solve_eig(Ad, Bd, k, C, Sm.L, rho);
    r = norm(Bd(:, :, k)*(W(end-rho+1:end, :)*Yh), 'fro');
    return
end
Y = small_solve(Ad, Bd, k, C, Sm, rho, n);
r = norm(Bd(:, :, k)*Y(end-rho+1:end, :), 'fro');      % unitary Qs does not change the norm
end
