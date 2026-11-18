function [k, alpha, beta, sigma] = classical_shifts_real(Eint, Fint, tol, kmax)
%CLASSICAL_SHIFTS_REAL  Smallest k such that a type-(k,k) rational reaches ratio <= TOL.
%   [K,ALPHA,BETA,SIGMA] = CLASSICAL_SHIFTS_REAL([a b],[c d],TOL,KMAX) gives the
%   balanced (classical fADI) shifts for two real intervals, computed by OPT_SHIFTS_REAL
%   (robust for intervals spanning many decades, where ZOLO3's pole/zero extraction is
%   not).  ALPHA (zeros) are the right shifts, BETA (poles) the left shifts.
%   SIGMA is the ratio max_E|r|/min_F|r| reached.  k runs from 2 to KMAX (default 80);
%   a warning is issued if TOL is not reached.
if nargin < 4, kmax = 80; end
for k = 2:kmax
    [alpha, beta, sigma] = opt_shifts_real(Eint, Fint, k, k);
    if sigma <= tol, return, end
end
warning('classical_shifts_real:tol', 'tolerance not reached with k = %d (sigma = %.2e)', kmax, sigma);
end
