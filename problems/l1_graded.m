function [T, t, info] = l1_graded(n, alpha, r, Tend)
%L1_GRADED  The L1 discretisation of the Caputo derivative on a graded mesh, as a matrix.
%   [T, t, info] = l1_graded(n, alpha, r, Tend)
%   Mesh t_j = Tend*(j/n)^r, j = 0..n; T (n x n, lower triangular) with (D^alpha w)(t_k) ~ (T w)_k,
%   w_0 = 0, T_kk = tau_k^(-alpha)/Gamma(2-alpha).  r = 1 is uniform; default r = (2-alpha)/alpha
%   (optimal order 2-alpha for the t^alpha singularity).  Tend defaults to 1.
%   t: the column t_1..t_n
%   info: .reW = [min, max] eigenvalue of (T + T')/2 = Re W(T), .diag = [min, max] of diag(T), .tau
%   Example: [T, t] = l1_graded(100, 0.5);
if nargin < 3 || isempty(r), r = (2 - alpha)/alpha; end
if nargin < 4, Tend = 1; end
tm = Tend*((0:n)'/n).^r;
tau = diff(tm);
tt = tm(2:end);
b = 1 - alpha;
Y = tt - tt.';                                         % Y(k,j) = t_k - t_j
D = zeros(n);
low = Y > 0;
Tj = repmat(tau.', n, 1);
% (Y + tau)^b - Y^b without cancellation: tau_1 can be far below the rounding level of t_k
D(low) = Y(low).^b.*expm1(b*log1p(Tj(low)./Y(low)));
D(1:n+1:end) = tau.^b;                                  % j = k: (t_k - t_{k-1})^b
a = D./(gamma(2 - alpha)*tau.');
T = tril(a - [a(:, 2:end), zeros(n, 1)]);
t = tt;
if nargout > 2
    H = (T + T.')/2;
    e = eig(H);
    info.reW = [min(e), max(e)];
    info.diag = [min(diag(T)), max(diag(T))];
    info.tau = tau;
end
end
