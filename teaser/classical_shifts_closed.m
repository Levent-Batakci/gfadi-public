function [k, alpha, beta, sigma] = classical_shifts_closed(E, F, tol, kmax)
%CLASSICAL_SHIFTS_CLOSED  Zolotarev type-(k,k) shifts for real intervals F < 0 < E, in closed form.
%   [K, ALPHA, BETA, SIGMA] = CLASSICAL_SHIFTS_CLOSED(E, F, TOL, KMAX) returns the smallest k <= KMAX
%   (default 200) with SIGMA = max_E|r|/min_F|r| <= TOL, the zeros ALPHA and poles BETA of r (k each).
%   E = [a b] with 0 < a < b, F = [c d] with c < d < 0.  Needs no Optimization Toolbox.
%   Method: a Moebius map T takes F to [-g,-1] and E to [1,g]; there the zeros are g*dn((2j-1)K/(2k), k')
%   and the poles their negatives; both are mapped back.  SIGMA is measured on fine grids of E and F.
%   Used for classical fADI in TEASER_RUN (Fig. 1).
if nargin < 4, kmax = 200; end
c = F(1); d = F(2); a = E(1); b = E(2);
M = abs((c - a)*(d - b)/((c - b)*(d - a)));            % cross ratio
g = -1 + 2*M + 2*sqrt(M^2 - M);
% Moebius through d -> -1, a -> 1, b -> g (c -> -g follows)
S  = @(z, z1, z2, z3) ((z - z1).*(z2 - z3))./((z - z3).*(z2 - z1));
Si = @(s, z1, z2, z3) (z3.*s.*(z2 - z1) - z1.*(z2 - z3))./(s.*(z2 - z1) - (z2 - z3));
Tinv = @(w) Si(S(w, -1, 1, g), d, a, b);
kp = sqrt(1 - 1/g^2);
K = ellipke(kp^2);
gE = a + (b - a)*(1 - cos(pi*(0:4000)'/4000))/2;  gE = unique([gE; logspace(log10(a), log10(b), 4000)']);
gF = c + (d - c)*(1 - cos(pi*(0:4000)'/4000))/2;  gF = unique([gF; -logspace(log10(-d), log10(-c), 4000)']);
for k = 1:kmax
    [~, ~, dn] = ellipj((2*(1:k)' - 1)*K/(2*k), kp^2);
    wz = g*dn;
    alpha = Tinv(wz);  beta = Tinv(-wz);
    lr = @(z) sum(log(abs(z - alpha.')), 2) - sum(log(abs(z - beta.')), 2);
    sigma = exp(max(lr(gE)) - min(lr(gF)));
    if sigma <= tol, return, end
end
warning('tolerance not reached');
end
