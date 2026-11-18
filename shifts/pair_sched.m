function [alpha, beta, sched] = pair_sched(alpha, beta, E, place)
%PAIR_SCHED  Stable step order for a lopsided type-(p,q) iteration, p >= q.
%   [ALPHA,BETA,SCHED] = PAIR_SCHED(ALPHA,BETA,E,PLACE) pairs every pole beta_j with the
%   zero alpha minimising max_{z in E} |z - alpha|/|z - beta_j|, so each pair is one
%   compressed alternating step with a factor of modest size on E.  The remaining p - q
%   zeros are right steps in Leja order.  PLACE: 'spread' (default) spreads the pairs
%   evenly among the right steps, 'first' puts them first, 'last' last.
%   E: sample vector of the spectral set of A (or [a b] for an interval).
%   For p < q use ORDER_SHIFTS (it calls PAIR_SCHED with the roles swapped).
if nargin < 4, place = 'spread'; end
if numel(E) == 2 && isreal(E)
    % interval: Chebyshev samples plus samples clustered geometrically at the left end
    a = min(E);  b = max(E);  n = 2000;
    E = [(a + b)/2 + (b - a)/2*cos(pi*(0:n-1)'/(n-1)); a + (b - a)*10.^(-(0:0.02:8)')];
end
alpha = alpha(:);  beta = beta(:);
p = numel(alpha);  q = numel(beta);
if p < q
    error('pair_sched:type', 'for p < q call [b,a,s] = pair_sched(beta, alpha, F, place) and swap R/L in s');
end
% poles in order of their distance to E (closest first, the hardest)
[~, ob] = sort(min(abs(E(:) - beta.'), [], 1));
beta = beta(ob);
free = true(p, 1);  partner = zeros(q, 1);
LE = log(abs(E(:) - alpha.'));
for j = 1:q
    cost = max(LE - log(abs(E(:) - beta(j))), [], 1).';
    cost(~free) = Inf;
    [~, i] = min(cost);
    partner(j) = i;
    free(i) = false;
end
rest = leja_order(alpha(free), E);
switch place
    case 'first', pos = 1:q;
    case 'last',  pos = (p - q) + (1:q);
    otherwise,    pos = round(linspace(1, p, q + 2));  pos = pos(2:end-1);
end
% build the sequence of "units": a unit is either a pair (RL) or a single right step
nunit = p;
isPair = false(nunit, 1);
isPair(pos) = true;
alpha_out = zeros(p, 1);
sched = blanks(0);
ip = 0;  ir = 0;  ka = 0;
for u = 1:nunit
    ka = ka + 1;
    if isPair(u)
        ip = ip + 1;
        alpha_out(ka) = alpha(partner(ip));
        sched = [sched 'RL'];     %#ok<AGROW>
    else
        ir = ir + 1;
        alpha_out(ka) = rest(ir);
        sched = [sched 'R'];      %#ok<AGROW>
    end
end
alpha = alpha_out;
end
