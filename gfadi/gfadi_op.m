function op = gfadi_op(M, kind)
%GFADI_OP  Operator struct for GFADI from a matrix.
%   OP = GFADI_OP(M)            general (dense or sparse) matrix; every shifted solve
%                               is a call to backslash
%   OP = GFADI_OP(M, 'tri')     dense triangular M: shifted solves by triangular
%                               substitution, O(n^2) each
%   OP = GFADI_OP(d, 'diag')    diagonal matrix with diagonal d
%   OP = GFADI_OP(Mfun, 'free') matrix-free: Mfun(X) = M*X only; adjoint products and
%                               shifted solves are unavailable (calling them is an error)
%   Fields:  .mv(X) = M*X   .mvh(X) = M'*X   .sv(s,X) = (M - s*I)\X
%            .svh(s,X) = (M - s*I)'\X   .n = dimension
%   Example:  B = gfadi_op(diag(Bmat), 'diag');
if nargin < 2, kind = 'matrix'; end
switch kind
    case 'matrix'
        n = size(M, 1);
        if issparse(M), I = speye(n); else, I = eye(n); end
        op.mv  = @(X) M*X;
        op.mvh = @(X) M'*X;
        op.sv  = @(s, X) (M - s*I)\X;
        op.svh = @(s, X) (M - s*I)'\X;
        op.n   = n;
    case 'tri'                                   % dense triangular: substitution only, the
        n = size(M, 1);                          % adjoint stored once, LINSOLVE without checks
        I = eye(n);
        Mt = M';  up = istriu(M);
        op.mv  = @(X) M*X;
        op.mvh = @(X) Mt*X;
        op.sv  = @(s, X) linsolve(M - s*I, X, struct('UT', up, 'LT', ~up));
        op.svh = @(s, X) linsolve(Mt - conj(s)*I, X, struct('UT', ~up, 'LT', up));
        op.n   = n;
    case 'diag'
        d = M(:);
        op.mv  = @(X) d.*X;
        op.mvh = @(X) conj(d).*X;
        op.sv  = @(s, X) X./(d - s);
        op.svh = @(s, X) X./conj(d - s);
        op.n   = numel(d);
    case 'free'
        op.mv  = M;
        op.mvh = @(X) error('gfadi_op:free', 'adjoint product not available');
        op.sv  = @(s, X) error('gfadi_op:free', 'shifted solves not available');
        op.svh = op.sv;
        op.n   = NaN;
    otherwise
        error('gfadi_op:kind', 'unknown kind');
end
end
