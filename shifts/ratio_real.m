function [sigma, maxE, minF] = ratio_real(Eint, Fint, alpha, beta)
%RATIO_REAL  Certified Zolotarev ratio  max_E|r| / min_F|r|  for real intervals.
%   [SIGMA,MAXE,MINF] = RATIO_REAL([a b],[c d],ALPHA,BETA), r(z) = prod(z - alpha)/prod(z - beta).
%   Each interval is split at the zeros (on E) or poles (on F) inside it, and every piece
%   is sampled with 32 Chebyshev points plus geometric clusters at both ends when it spans
%   more than a decade, so an optimiser cannot hide peaks between grid points.
%   MAXE, MINF are log max_E|r| and log min_F|r| (computed in logs: large p is fine).
alpha = real(alpha(:));  beta = real(beta(:));
lr = @(X) sum(log(abs(X - alpha.')), 2) - sum(log(abs(X - beta.')), 2);
XE = pieces(Eint, alpha);  XF = pieces(Fint, beta);
maxE = -Inf;  minF = Inf;
ch = 4000;
for i = 1:ch:numel(XE), maxE = max(maxE, max(lr(XE(i:min(i+ch-1, end))))); end
for i = 1:ch:numel(XF)
    v = lr(XF(i:min(i+ch-1, end)));  v = v(isfinite(v));
    if ~isempty(v), minF = min(minF, min(v)); end
end
sigma = exp(maxE - minF);
end

function X = pieces(I, br)
a = min(I);  b = max(I);
br = sort([a; b; br(br > a & br < b)]);
t = (1 - cos(pi*(1:31)'/32))/2;                       % interior Chebyshev points of [0,1]
X = cell(numel(br) - 1, 1);
for k = 1:numel(br) - 1
    x0 = br(k);  x1 = br(k+1);  L = x1 - x0;
    if L <= 0, X{k} = zeros(0, 1); continue, end
    xs = x0 + L*t;
    s0 = max(abs(x0), 1e-3*L);  s1 = max(abs(x1), 1e-3*L);
    if L > 10*min(s0, s1)                              % wide piece: add geometric clusters
        g = L*10.^(-(0:1/40:log10(L/(1e-3*min(s0, s1))))');
        xs = [xs; x0 + g; x1 - g];                     %#ok<AGROW>
        xs = xs(xs > x0 & xs < x1);
    end
    X{k} = xs;
end
X = [a; b; vertcat(X{:})];
end
