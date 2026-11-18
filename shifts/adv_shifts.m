function [p, al, info] = adv_shifts(E, c, ymax, WL, tol)
%ADV_SHIFTS  Real shifts of one-sided g-fADI for a nonsymmetric large side A = H + S.
%   [P,AL,INFO] = ADV_SHIFTS(E,C,YMAX,WL,TOL) returns P real shifts AL (Chebyshev points
%   of an interval [a', b'] near E) for which the bound
%       sigma = 2 max_Omega |q| / min_F |q|,   q(z) = prod(z - al_j),
%   is <= TOL, with P as small as the search finds.
%   E = [E1 E2]: spectrum of the symmetric part H.  W(A) lies in
%       Omega = { x + i y :  E1 <= x <= E2,  |y| <= min(C sqrt(x), YMAX) }.
%   WL = [min max] of Re W(L) for the small side; F = [-WL(2), -WL(1)].
%   H symmetric, S skew-symmetric.  C = 0 gives the Chebyshev shifts of E (normal A,
%   factor 1).  P = Inf and AL = [] if no interval in the search reaches TOL.
%   INFO: .sigma .aprime .bprime .cr (Crouzeix factor) .psym (P for a symmetric A)
%         .c .ymax .pmax .limit
%   Used for the model problem (Fig. 1), Example 3 with a swirl (SWITCH_PROBLEM) and
%   Example 4 (FREQ_SWEEP_RUN), where C and YMAX are derived; see NOTES.md.
%   Example:  [p, al] = adv_shifts([e1 e2], c, ymax, [l1 l2], 1e-8);
cr = 2;                          % Crouzeix's constant 2 (Jin 2026, Lorist--Schwenninger 2026); Crouzeix--Palencia proved 1 + sqrt 2
if c == 0, cr = 1; end
F1 = -WL(1);                                                      % min_F |q| is attained here for zeros in [a', E2] > 0 > F
psym = search(E(1), E(2), E, 0, ymax, F1, tol, 1, 8, Inf);
best = struct('p', Inf, 'a', NaN);
pmax = 30*psym;
% candidate left ends a' = fa*E1: only a' > c^2/4 can converge (see NOTES.md)
if c > 0, fas = [2 1.6 1.3 1 0.8 0.6 0.45 0.32 0.22 0.15 0.1];  fas = fas(fas*E(1) > 1.02*c^2/4); else, fas = 1; end
% candidate right ends b' = E2 + g*y(E2)^2/(4 a') (see NOTES.md)
yE2 = min(c*sqrt(E(2)), ymax);
gs = [0 0.5 1 2 4];
if c == 0, gs = 0; end
for fa = fas
    for gg = gs
        a1 = fa*E(1);
        b1 = E(2) + gg*yE2^2/(4*a1);
        pp = search(a1, b1, E, c, ymax, F1, tol, cr, max(psym, 1), min(pmax, best.p));
        if pp < best.p, best = struct('p', pp, 'a', a1, 'b', b1); end
    end
end
p = best.p;  al = [];  sig = NaN;
if isfinite(p)
    al = cheb_nodes(best.a, best.b, p);
    sig = ratio(best.a, best.b, p, E, c, ymax, F1, cr, 0.1);
end
if ~isfield(best, 'b'), best.b = NaN; end
info = struct('sigma', sig, 'aprime', best.a, 'bprime', best.b, 'cr', cr, 'psym', psym, 'c', c, 'ymax', ymax, 'pmax', pmax, 'limit', 2*sqrt(E(1)));
end

function p = search(a1, b1, E, c, ymax, F1, tol, cr, p0, pmax)
% smallest p <= pmax with sigma <= tol for Chebyshev points of [a1, b1]: growth by half, then bisection; Inf if none
lo = 0;  hi = max(8, round(p0));
if c == 0 && a1 == E(1) && b1 == E(2), f = @(q) cheb1(E, F1, q); else, f = @(q) ratio(a1, b1, q, E, c, ymax, F1, cr, 0.25); end
while f(hi) > tol
    lo = hi;
    if hi >= pmax, p = Inf; return; end
    hi = min(ceil(1.5*hi), ceil(pmax));
end
while hi - lo > max(1, round(0.01*hi))
    mid = floor((lo + hi)/2);
    if f(mid) <= tol, hi = mid; else, lo = mid; end
end
p = hi;
end

function s = cheb1(E, F1, p)
% the interval itself (a normal A): the closed form of CHEB_SHIFTS
[~, s] = cheb_shifts(E(1), E(2), p, F1, F1);
end

function a = cheb_nodes(a1, b, p)
a = (a1 + b)/2 + (b - a1)/2*cos((2*(1:p)' - 1)*pi/(2*p));
end

function s = ratio(a1, b1, p, E, c, ymax, F1, cr, step)
% cr * max over the upper boundary of Omega of |q| / |q(F1)|, q with the p Chebyshev points
% of [a1, b1] as zeros.  Samples STEP*y apart along the top curve plus the two vertical
% edges; |q| in closed form via T_p, in logs (O(1) per sample for any p).
ytop = @(x) min(c*sqrt(x), ymax);
x = E(1);  xs = zeros(1, 0);
while x < E(2) && numel(xs) < 2e5, xs(end + 1) = x;  x = x + step*ytop(x); end   %#ok<AGROW>
xs = [xs E(2)]';
z = xs + 1i*ytop(xs);
z = [z; E(1) + 1i*ytop(E(1))*linspace(0, 1, 60)'; E(2) + 1i*ytop(E(2))*linspace(0, 1, 60)'];
tz = @(z) (2*z - a1 - b1)/(b1 - a1);
s = cr*exp(max(logT(tz(z), p)) - logT(tz(F1), p));
end

function l = logT(t, p)
% log |T_p(t)| = p log|w| + log|1 + w^(-2p)| - log 2,  w = t + sqrt(t - 1) sqrt(t + 1),  |w| >= 1 (stable for large p)
w = t + sqrt(t - 1).*sqrt(t + 1);
lw = log(complex(w));
l = p*real(lw) + log(abs(1 + exp(-2*p*lw))) - log(2);
end
