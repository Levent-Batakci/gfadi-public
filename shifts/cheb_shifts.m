function [alpha, xi] = cheb_shifts(a, b, p, c, d)
%CHEB_SHIFTS  Optimal type-(p,0) shifts for a real interval E = [a,b].
%   ALPHA = CHEB_SHIFTS(A,B,P) returns the zeros of the degree-P Chebyshev polynomial
%   mapped to E (column vector).
%   [ALPHA,XI] = CHEB_SHIFTS(A,B,P,C,D) also returns the Zolotarev number xi_{p,0} for
%   F = [c,d] (disjoint from E, on either side): 1/|T_p| at the end of F nearest E.
%   Example:  [alpha, xi] = cheb_shifts(1, 100, 10, -100, -1);
j = (1:p)';
xi = NaN;
alpha = (a + b)/2 + (b - a)/2*cos((2*j - 1)*pi/(2*p));
if nargin > 3
    x = @(t) (2*t - a - b)/(b - a);
    xn = min(abs(x([c d])));                 % nearest point of F to E in the map
    xi = 1/real(cosh(p*acosh(xn)));          % |T_p(x)| for |x| >= 1
end
end
