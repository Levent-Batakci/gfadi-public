function [alpha, beta, sched] = order_shifts(alpha, beta, E, F, place)
%ORDER_SHIFTS  Recommended step order for any type (p,q).
%   [ALPHA,BETA,SCHED] = ORDER_SHIFTS(ALPHA,BETA,E,F,PLACE) reorders the shifts and returns
%   the step string SCHED for GFADI / GFADI_LOWRANK ('sched' option).
%   p >= q: PAIR_SCHED(ALPHA,BETA,E,PLACE), poles paired with their best zeros, the
%           remaining zeros in Leja order with respect to E.
%   p <  q: the mirror image, with respect to F.
%   q = 0 or p = 0: Leja order of the only list.
%   E, F: sample vectors of the spectral sets (or [a b] for intervals).
%   PLACE: where the pairs go, 'first' (default), 'last' or 'spread'.
if nargin < 5, place = 'first'; end
alpha = alpha(:);  beta = beta(:);
if isempty(beta)
    alpha = leja_order(alpha, E(:));
    sched = repmat('R', 1, numel(alpha));
elseif isempty(alpha)
    beta = leja_order(beta, F(:));
    sched = repmat('L', 1, numel(beta));
elseif numel(alpha) >= numel(beta)
    [alpha, beta, sched] = pair_sched(alpha, beta, E, place);
else
    [beta, alpha, s] = pair_sched(beta, alpha, F, place);
    sched = s;
    sched(s == 'R') = 'L';
    sched(s == 'L') = 'R';
end
end
