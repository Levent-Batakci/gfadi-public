function out = zolo3(E, F, n, varargin)
%ZOLO3  Third and fourth Zolotarev problems on discrete sets, via constrained AAA-Lawson.
%   OUT = ZOLO3(E,F,N,...) computes a near-optimal degree-N rational r minimising
%       max_{z in E}|r(z)| / min_{z in F}|r(z)|                         (Z3)
%   on sample sets E, F (complex column vectors), and the associated approximant to
%   sign_{E/F} (-1 on E, +1 on F)                                       (Z4).
%   Options:  'delta'   (0)   degree drop D: r has type (N-D, N), i.e. D fewer zeros
%                             (D fewer solves with that ADI matrix); 0 <= D < N+1
%             'lawson'  (300) number of Lawson reweighting steps
%             'damping' (0.9) Lawson damping in (0,1]; 1 often diverges
%             'sweeps'  (3)   inner alternations per Lawson step
%   OUT: .sigma Z3 ratio  .tau Z4 error  .zer (N-D zeros: shifts for the E side)
%        .pol (N poles: shifts for the F side)  .type [N-D, N]
%        .r3 handle, Z3 rational with min_F|r3| = 1   .r4 handle, Z4 approximant
%        .t,.S,.T support points and barycentric coefficients
%        .errE,.errF Z4 error at non-support samples, selected by .keepE,.keepF
%   Example:  x = cos(pi*(0:799)'/799);  out = zolo3(-1.25 + 0.75*x, 2 + x, 10, 'delta', 2);

%   References:
%     E. I. Zolotarev, Zap. Imp. Akad. Nauk St. Petersburg 30 (1877) 1-59.
%     M.-P. Istace & J.-P. Thiran, SIAM J. Numer. Anal. 32 (1995) 249-259.
%     G. Starke, J. Approx. Theory 70 (1992) 115-130.
%     Y. Nakatsukasa & L. N. Trefethen, SIAM J. Sci. Comput. 42 (2020) A3157-A3179.
%     L. N. Trefethen & H. D. Wilber, SIAM J. Sci. Comput. 47 (2025) A2205-A2220.

% ---------------------------------------------------------------- options
p = inputParser;
p.addParameter('delta',   0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('lawson',  300, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('damping', 0.9, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('sweeps',  3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
opt = p.Results;

E = E(:);  F = F(:);
m = n + 1;                            % number of support points
d = round(opt.delta);
if d >= m
    error('zolo3:delta', 'delta must be < n+1.');
end
if numel(E) < m || numel(F) < m
    error('zolo3:samples', 'Need at least n+1 samples in each of E and F.');
end

Z    = [E; F];
f    = [-ones(numel(E),1); ones(numel(F),1)];
nE   = numel(E);

% ------------------------------------------- phase 1: greedy AAA, support points
% Interpolatory AAA only PLACES the support points.  The Trefethen-Wilber singular-vector
% blend (BLENDVEC) is essential: without it the greedy phase makes no progress on sign data.
J   = true(numel(Z),1);
sup = zeros(m,1);
R   = mean(f)*ones(numel(Z),1);
for k = 1:m
    [~, j] = max(abs(f - R).*J);
    sup(k) = j;  J(j) = false;
    ts = Z(sup(1:k));  fs = f(sup(1:k));
    C  = 1./(Z(J) - ts.');
    A  = (f(J) - fs.').*C;
    w  = blendvec(A);
    R       = complex(f);
    R(J)    = (C*(w.*fs))./(C*w);
end
t = Z(sup);

% ------------------------------------------- phase 2: Lawson, non-interpolatory
% r3 = S(z)/T(z) with S, T independent barycentric sums over all support points:
% S small on E, T small on F, rows weighted by 1/|T| (on E) and 1/|S| (on F)
keepE = J(1:nE);  keepF = J(nE+1:end);
ZE = E(keepE);    ZF = F(keepF);
CE = 1./(ZE - t.');
CF = 1./(ZF - t.');

Q = typebasis(t, d);                  % null-space basis enforcing the type constraint

lamE = ones(numel(ZE),1)/numel(ZE);
lamF = ones(numel(ZF),1)/numel(ZF);

S = Q*minvec(CE*Q);
T = minvec(CF);

bestSig = Inf;  bestS = S;  bestT = T;
for k = 1:opt.lawson
    for s = 1:opt.sweeps
        S = Q*minvec((sqrt(lamE)./abs(CE*T)).*CE*Q);
        T =   minvec((sqrt(lamF)./abs(CF*S)).*CF);
    end
    epsE = abs((CE*S)./(CE*T));       % |r3| on E, up to the free scaling
    epsF = abs((CF*T)./(CF*S));       % |1/r3| on F
    sig  = max(epsE)*max(epsF);       % scaling-invariant, = Z3 ratio
    if sig < bestSig
        bestSig = sig;  bestS = S;  bestT = T;
    end
    % damped Lawson reweighting
    lamE = lamE.*((1-opt.damping) + opt.damping*epsE/max(epsE));
    lamE = lamE/sum(lamE);
    lamF = lamF.*((1-opt.damping) + opt.damping*epsF/max(epsF));
    lamF = lamF/sum(lamF);
end
S = bestS;  T = bestT;  sigma = bestSig;

% ------------------------------------------- balancing, normalisation, output
mE  = max(abs((CE*S)./(CE*T)));
mF  = max(abs((CF*T)./(CF*S)));
rho = sqrt(mE/mF);                    % equalises the Z4 error maxima on E and F
T   = rho*T;                          % absorb it, so r3 = S/T and min_F|r3| = 1 below

scal  = 1/min(abs(baryratio(F, t, S, T)));   % normalise so min_F |r3| = 1
S     = scal*S;

out.t     = t;
out.S     = S;
out.T     = T;
out.r3    = @(z) reshape(baryratio(z, t, S, T), size(z));
out.r4    = @(z) reshape(baryratio(z, t, S - scal*T, S + scal*T), size(z));
out.sigma = sigma;
out.tau   = 2*sqrt(sigma)/(1 + sigma);
out.type  = [n-d, n];
out.keepE = keepE;  out.keepF = keepF;
out.errE  = out.r4(ZE) - (-1);
out.errF  = out.r4(ZF) - ( 1);

[out.zer, out.pol] = polezeros(t, S, T);

end
% =====================================================================

function v = baryratio(z, t, cnum, cden)
%BARYRATIO  Evaluate (sum cnum_k/(z-t_k)) / (sum cden_k/(z-t_k)).
%   Both sums have a simple pole at every support point and the ratio does not,
%   so exact hits z = t_k are handled by the limit cnum_k/cden_k rather than
%   returning Inf/Inf = NaN.  This matters because the sample sets E and F
%   contain the support points.
z = z(:);
[hit, loc] = ismember(z, t);
D = z - t.';
D(hit, :) = 1;                        % placeholder; those rows are overwritten
v = sum(cnum.'./D, 2)./sum(cden.'./D, 2);
v(hit) = cnum(loc(hit))./cden(loc(hit));
end

function w = blendvec(A)
%BLENDVEC  Inverse-square blend of right singular vectors (Trefethen-Wilber 'sign').
[~, Sv, V] = svd(A, 0);
s = diag(Sv);
s = max(s, eps*s(1));
w = V*(1./s.^2);
w = w/norm(w);
end

function v = minvec(A)
%MINVEC  Right singular vector of the smallest singular value.
[~, ~, V] = svd(A, 0);
v = V(:,end);
end

function Q = typebasis(t, d)
%TYPEBASIS  Orthonormal basis of {S : sum_k S_k t_k^j = 0, j = 0..d-1}.
m = numel(t);
if d <= 0
    Q = eye(m);
    return
end
c = mean(t);
s = max(abs(t - c));
if s == 0, s = 1; end
tt = (t - c)/s;                        % shift+scale: the raw Vandermonde is hopeless
B  = zeros(d, m);
for j = 0:d-1
    B(j+1,:) = tt.'.^j;
end
Q = null(B);
end

function [zer, pol] = polezeros(t, S, T)
%POLEZEROS  Roots of the two barycentric numerators, by generalised eigenvalue problem.
%   r3 = S(z)/T(z) and the common factor prod(z-t_k) cancels, so the zeros of r3 are
%   the roots of the polynomial numerator of S and the poles those of T.  A degree
%   drop of d shows up as d extra infinite eigenvalues, which are discarded here.
%   This is a direct check that the type constraint took effect.
m = numel(t);
B = eye(m+1);  B(1,1) = 0;
out = cell(1,2);
cs  = {S, T};
for i = 1:2
    Emat = zeros(m+1);
    Emat(1, 2:end)     = cs{i}.';
    Emat(2:end, 1)     = 1;
    Emat(2:end, 2:end) = diag(t);
    ev = eig(Emat, B);
    out{i} = ev(isfinite(ev));
end
zer = out{1};
pol = out{2};
end
