function [alpha, beta, sigma, info] = opt_shifts_real(Eint, Fint, p, q, varargin)
%OPT_SHIFTS_REAL  Near-optimal type-(p,q) shifts for two disjoint REAL intervals.
%   [ALPHA,BETA,SIGMA,INFO] = OPT_SHIFTS_REAL([a b],[c d],P,Q,...) minimises the Zolotarev
%   ratio max_E|r|/min_F|r|, r = prod(z - alpha_i)/prod(z - beta_j), with real zeros in E
%   and real poles in F, by quasi-Newton (fminunc) on a smoothed (log-sum-exp) objective,
%   sharpened in stages.  SIGMA is the certified ratio at the returned shifts (RATIO_REAL),
%   an upper bound for xi_{p,q}(E,F) and what g-fADI attains.
%   Complements ZOLO3, which loses the type constraint when p and q are far apart.
%   Needs the Optimization Toolbox (fminunc).
%   Options:  'alpha0','beta0' initial shifts          'temps' smoothing schedule ([8 64 512])
%             'maxit' iterations per stage (400)       'nsamp' samples per set (3*(p+q)+200)
%             'nknots' 0 = optimise every zero; default 0 for p < 150, else 24
%   INFO.v: the final optimisation variables.
%   Example:  [alpha, beta, sigma] = opt_shifts_real([1 100], [-100 -1], 8, 4);
pr = inputParser;
pr.addParameter('nsamp', 3*(p + q) + 200);
pr.addParameter('temps', [8 64 512]);
pr.addParameter('beta0', []);
pr.addParameter('alpha0', []);
pr.addParameter('maxit', 400);
pr.addParameter('nknots', -1);
pr.parse(varargin{:});
o = pr.Results;
if o.nknots < 0, o.nknots = 24*(p >= 150); end

a = Eint(1); b = Eint(2); c = Fint(1); d = Fint(2);
mE = (a + b)/2; hE = (b - a)/2; mF = (c + d)/2; hF = (d - c)/2;
% Chebyshev-clustered samples, offset by an irrational fraction so that no sample can
% coincide with an initial Chebyshev zero; the two endpoints are added explicitly
x  = cos(pi*((0:o.nsamp-1)' + 0.3711)/o.nsamp);
x  = [-1; x; 1];
XE = [mE + hE*x; geomsamples(a, b)];
XF = [mF + hF*x; geomsamples(c, d)];

% initial guess: geometric spacing near the close end when an interval is long compared
% with the gap (linear Chebyshev starts end in a poor local minimum there)
dEF = max(c, a) - min(d, b);                         % distance between the intervals
if isempty(o.alpha0)
    th = (2*(1:p)' - 1)*pi/(2*p);
    th = geoinit(th, a, b, c, d, dEF, mE, hE);
else
    th = acos(max(min((o.alpha0(:) - mE)/hE, 1), -1));
end
th = sort(th);
if isempty(o.beta0)
    if q > 0
        ph = (2*(1:q)' - 1)*pi/(2*q);
        ph = geoinit(ph, c, d, a, b, dEF, mF, hF);
    else, ph = zeros(0, 1); end
else
    ph = acos(max(min((o.beta0(:) - mF)/hF, 1), -1));
end

nk = o.nknots;
if nk > 0
    % theta = pi*u + W*delta, W = piecewise-linear hat functions on clustered knots
    u  = ((1:p)' - 0.5)/p;
    kn = (1 - cos(pi*(0:nk)'/nk))/2;
    W  = hatmatrix(u, kn);
    base = pi*u;
    % fit delta to the initial angles (least squares)
    delta0 = W\(th - base);
    v = [delta0; ph];
else
    W = [];  base = [];  v = [th; ph];
end

for t = o.temps
    fun = @(v) objgrad(v, p, q, mE, hE, mF, hF, XE, XF, t, W, base);
    optn = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, ...
                        'Display', 'off', 'MaxIterations', o.maxit, 'MaxFunctionEvaluations', 20*o.maxit, ...
                        'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-12);
    v = fminunc(fun, v, optn);
end
[th, ph] = angles(v, p, q, W, base);
alpha = mE + hE*cos(th);  beta = mF + hF*cos(ph);
% certified ratio: every gap between consecutive zeros (poles) is sampled
sigma = ratio_real(Eint, Fint, alpha, beta);
info.v = v;
end

function s = geomsamples(a, b)
%GEOMSAMPLES  Points clustering geometrically at both ends of [a,b], 40 per decade of the
%   distance to the end (needed for intervals spanning many decades).
L = abs(b - a);
tmin = 1e-3*max(min(abs([a b])), 1e-6*L);       % smallest distance to an end that is sampled
nd = ceil(log10(L/tmin));
t = L*10.^(-(0:1/40:nd)');
s = [a + sign(b - a)*t; b - sign(b - a)*t];
s = s(s > min(a, b) & s < max(a, b));
s = [a; b; s];
end

function th = geoinit(th, a, b, c, d, dEF, m, h)
% map Chebyshev angles of [a,b] to points whose distance to the end of [a,b] nearest
% [c,d] is spread geometrically between dEF and the length of [a,b]
L = b - a;
if L/dEF <= 50, return, end
if b <= c, e0 = b; sg = -1; else, e0 = a; sg = 1; end   % near end, direction into [a,b]
u = (1 - cos(th))/2;                                     % in [0,1]
dist = dEF*(L/dEF).^u;                                   % from dEF to L
dist = min(dist, L);
x = e0 + sg*(dist - dEF*0.5);                           % keep strictly inside
x = min(max(x, a), b);
th = acos(max(min((x - m)/h, 1), -1));
end

function [th, ph] = angles(v, p, q, W, base)
if isempty(W), th = v(1:p); else, th = base + W*v(1:size(W, 2)); end
ph = v(end-q+1:end);  if q == 0, ph = zeros(0, 1); end
end

function W = hatmatrix(u, kn)
% piecewise-linear interpolation matrix (rows: points u, columns: knots kn)
n = numel(u);  nk = numel(kn);
idx = min(max(discretize(u, kn), 1), nk - 1);
t = (u - kn(idx))./(kn(idx+1) - kn(idx));
W = sparse([(1:n)'; (1:n)'], [idx; idx+1], [1 - t; t], n, nk);
end

function [f, g] = objgrad(v, p, q, mE, hE, mF, hF, XE, XF, t, W, base)
[th, ph] = angles(v, p, q, W, base);
al = mE + hE*cos(th);  be = mF + hF*cos(ph);
DE = XE - al.';  DF = XF - al.';                    % samples x p
LE = sum(log(abs(DE)), 2);  LF = sum(log(abs(DF)), 2);
if q > 0
    GE = XE - be.';  GF = XF - be.';
    LE = LE - sum(log(abs(GE)), 2);  LF = LF - sum(log(abs(GF)), 2);
end
% smooth max over E of LE and smooth max over F of -LF
sE = t*LE;  mxE = max(sE);  wE = exp(sE - mxE);  SE = sum(wE);  fE = (mxE + log(SE))/t;  wE = wE/SE;
sF = -t*LF; mxF = max(sF);  wF = exp(sF - mxF);  SF = sum(wF);  fF = (mxF + log(SF))/t;  wF = wF/SF;
f = fE + fF;
if nargout > 1
    % d LE/d alpha_i = -1/(X - alpha_i),  d LE/d beta_j = +1/(X - beta_j)
    gA = -( (wE.'*(1./DE)).' ) + ( (wF.'*(1./DF)).' );     % fF has weight on -LF
    gth = gA.*(-hE*sin(th));
    if isempty(W), g = gth; else, g = W'*gth; end
    if q > 0
        gB = ( (wE.'*(1./GE)).' ) - ( (wF.'*(1./GF)).' );
        g = [g; gB.*(-hF*sin(ph))];
    end
end
end

