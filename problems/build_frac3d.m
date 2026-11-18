function P = build_frac3d(N, n, opts)
%BUILD_FRAC3D  Time-fractional (Caputo) subdiffusion in a heterogeneous 3-D medium.
%   P = build_frac3d(N, n, opts)
%   PDE: D_t^alpha u - div(a(x) grad u) = f(x) g(t), x in (0,1)^3, t in (0,Tend], u(.,0) = u0,
%   Dirichlet BCs in x.  Space: the 7-point operator of BUILD_AGE3D on N^3 nodes.  Time: the L1
%   scheme on the graded mesh t_j = Tend (j/n)^r (L1_GRADED).
%   With w = u - u0, X(:,k) = w(., t_k):  A X + X T^T = U V^T,  U = [f, -A u0],  V = [g(t_k), 1],
%   T dense lower triangular (n x n).  Used by the 'frac' case of EXP18_ANALYSIS; not in a
%   figure of the paper.
%   opts: alpha (0.5), r ((2-alpha)/alpha), Tend (1), contrast (10), afun/amin, nodal (false),
%         u0c/u0w (centre/width of the initial bump), fc/fw (centre/width of the source)
%   P: .A .T .t .U .V .u0 .m .n .N .opts, enclosures .lminA .lmaxA .reTmin .reTmax, .diagT ([min max])
%   Example: P = build_frac3d(12, 40);
if nargin < 3, opts = struct(); end
d = struct('alpha', 0.5, 'r', [], 'contrast', 10, 'afun', [], 'amin', 1, ...
           'u0c', [0.5 0.5 0.5], 'u0w', 0.15, 'fc', [0.3 0.4 0.5], 'fw', 0.1, 'nodal', false, 'Tend', 1);
fn = fieldnames(opts);
for i = 1:numel(fn), d.(fn{i}) = opts.(fn{i}); end
if isempty(d.r), d.r = (2 - d.alpha)/d.alpha; end
S = build_age3d(N, 2, struct('contrast', d.contrast, 'afun', d.afun, 'amin', d.amin, 'nodal', d.nodal));   % spatial operator only
[T, t, in] = l1_graded(n, d.alpha, d.r, d.Tend);
h = 1/(N + 1);
g = (1:N)*h;
[X3, Y3, Z3] = ndgrid(g, g, g);
bump = @(c, w) reshape(exp(-((X3 - c(1)).^2 + (Y3 - c(2)).^2 + (Z3 - c(3)).^2)/w^2), [], 1);
u0 = bump(d.u0c, d.u0w);                                   % initial concentration
f = 50*(0.1/d.fw)^3*bump(d.fc, d.fw);                      % source (fixed mass) ...
gt = exp(-((t - 0.5*d.Tend)/(0.1*d.Tend)).^2);             % ... switched on around t = Tend/2
if d.nodal                                                 % nothing starts inside a sealed wall
    open = S.anode > 0;
    u0 = u0.*open;
    f = f.*open;
end
P.A = S.A;
P.T = T;
P.t = t;
P.m = S.m;
P.n = n;
P.N = N;
P.opts = d;
P.U = [f, -S.A*u0];
P.V = [gt, ones(n, 1)];
P.u0 = u0;
P.lminA = S.lminA;
P.lmaxA = S.lmaxA;
P.reTmin = in.reW(1);
P.reTmax = in.reW(2);
P.diagT = in.diag;
end
