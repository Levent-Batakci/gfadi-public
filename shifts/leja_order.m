function [s, idx] = leja_order(s, S)
%LEJA_ORDER  Leja ordering of the points s with respect to a reference set S.
%   [S,IDX] = LEJA_ORDER(S,SREF): the first point is the one farthest from the centroid of
%   SREF (default: the points themselves); every next one maximises sum_j log|s - s_j|
%   over the points already chosen.  This keeps the partial products prod_{j<=k}(z - s_j)
%   moderate on SREF, so intermediate ADI residuals do not grow.
%   Example:  alpha = leja_order(cheb_shifts(1, 100, 10));
s = s(:);
if nargin < 2, S = s; end
n = numel(s);
if n == 0, idx = zeros(0, 1); return, end
c = mean(S(:));
idx = zeros(n, 1);
[~, i1] = max(abs(s - c));
idx(1) = i1;
lg = zeros(n, 1);
used = false(n, 1);  used(i1) = true;
for k = 2:n
    lg = lg + log(max(abs(s - s(idx(k-1))), realmin));
    lg(used) = -Inf;
    [~, i] = max(lg);
    idx(k) = i;  used(i) = true;
end
s = s(idx);
end
