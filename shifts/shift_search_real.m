function [p, alpha, beta, sigma] = shift_search_real(Eint, Fint, q, tol, pguess)
%SHIFT_SEARCH_REAL  Smallest p such that a type-(p,q) rational reaches ratio <= TOL.
%   [P,ALPHA,BETA,SIGMA] = SHIFT_SEARCH_REAL([a b],[c d],Q,TOL,PGUESS) for real intervals
%   E = [a b] and F = [c d]: Q poles (left steps) are fixed, P zeros (right steps) are
%   searched.  ALPHA (P zeros), BETA (Q poles), SIGMA the ratio max_E|r|/min_F|r| reached.
%   q = 0 uses the Chebyshev closed form (CHEB_SHIFTS).  For q > 0, a growth-then-bisection
%   search on p (to ~1%) with OPT_SHIFTS_REAL; after the first (cold) call each call is
%   warm-started from the zeros of the nearest p already solved.
%   PGUESS: starting p (default max(q+1, 20), capped at 60).
%   Example:  [p, alpha, beta, sigma] = shift_search_real([1 1e4], [-1e4 -1], 4, 1e-10);
if q == 0
    p = 1;
    while true
        [alpha, s] = cheb_shifts(Eint(1), Eint(2), p, Fint(1), Fint(2));
        if s <= tol, break, end
        p = p + max(1, floor(p/50));       % coarse then refine
    end
    while p > 1                            % refine downwards
        [a2, s2] = cheb_shifts(Eint(1), Eint(2), p - 1, Fint(1), Fint(2));
        if s2 > tol, break, end
        p = p - 1;  alpha = a2;  s = s2;
    end
    beta = [];  sigma = s;  return
end
if nargin < 5 || isempty(pguess), pguess = max(q + 1, 20); end
cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
ev = @(pp) evalp(pp, Eint, Fint, q, cache);
lo = 0;  hi = [];
pp = max(q + 1, min(round(pguess), 60));          % cheap cold start; warm starts do the rest
for it = 1:24
    s = ev(pp);
    if s <= tol
        hi = pp;
        if lo > 0, break, end
        pn = max(q + 1, floor(0.85*pp));
        if pn == pp, break, end
        s2 = ev(pn);
        if s2 > tol, lo = pn; break, else, hi = pn; pp = pn; end
    else
        lo = pp;  pp = ceil(1.4*pp) + 1;
    end
end
if isempty(hi), error('shift_search_real:fail', 'no p found'); end
while hi - lo > max(1, round(0.01*hi))      % bisection to ~1%
    mid = floor((lo + hi)/2);
    s = ev(mid);
    if s <= tol, hi = mid; else, lo = mid; end
end
p = hi;
r = cache(p);
sigma = r{1};  alpha = r{2};  beta = r{3};
end

function s = evalp(pp, Eint, Fint, q, cache)
% ratio for p = pp, memoised in CACHE (a handle object: {sigma, alpha, beta} per p)
if isKey(cache, pp), r = cache(pp); s = r{1}; return, end
ks = cell2mat(keys(cache));
if isempty(ks)
    [al, be, s] = opt_shifts_real(Eint, Fint, pp, q);
else
    % warm start: resample the Chebyshev angles of the nearest solved p
    [~, i] = min(abs(log(ks/pp)));
    r0 = cache(ks(i));  p0 = ks(i);
    mE = mean(Eint);  hE = diff(Eint)/2;
    th0 = acos(max(min((sort(r0{2}) - mE)/hE, 1), -1));       % angles of the old zeros, ascending in alpha
    u0 = ((1:p0)' - 0.5)/p0;  u = ((1:pp)' - 0.5)/pp;
    th = interp1([0; u0; 1], [0; th0(end:-1:1); pi], u, 'linear');   % theta increases as alpha decreases
    a0 = mE + hE*cos(th);
    [al, be, s] = opt_shifts_real(Eint, Fint, pp, q, 'alpha0', a0, 'beta0', r0{3}, ...
                                  'temps', [32 256 1024], 'maxit', 250);
end
cache(pp) = {s, al, be};
end

