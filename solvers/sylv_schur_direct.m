function [X, info] = sylv_schur_direct(A, L, U, V)
%SYLV_SCHUR_DIRECT  A X + X L^T = U V^T by Bartels-Stewart with a sparse large side.
%   [X, info] = sylv_schur_direct(A, L, U, V)
%   Schur form L^T = Q T Q^*, then n shifted sparse solves (A + T_jj I) xt_j = ft_j - sum_{i<j} T_ij xt_i,
%   X = Xt Q^*.  The complex Schur form is used when the real one is not triangular.
%   The reference solution of Example 2 (EXP16_BLOWOUT, EXP16B_NONNORMALITY, EXP18_ANALYSIS).
%   info: .complex (true if the complex form was needed), .time, .nfact (= n, the number
%         of shifted solves with A)
%   Example: X = sylv_schur_direct(P.A, P.L, P.U, P.V);
t0 = tic;
n = size(L, 1);
m = size(A, 1);
[Q, T] = schur(full(L.'), 'real');
info.complex = any(abs(diag(T, -1)) > 0);
if info.complex, [Q, T] = schur(full(L.'), 'complex'); end   % ignoring 2-by-2 blocks would be wrong
Ft = U*(V.'*Q);
if info.complex, Xt = complex(zeros(m, n)); else, Xt = zeros(m, n); end
I = speye(m);
for j = 1:n
    Xt(:, j) = (A + T(j, j)*I)\(Ft(:, j) - Xt(:, 1:j-1)*T(1:j-1, j));
end
X = Xt*Q';
if isreal(A) && isreal(L) && isreal(U) && isreal(V), X = real(X); end
info.time = toc(t0);
info.nfact = n;
end
