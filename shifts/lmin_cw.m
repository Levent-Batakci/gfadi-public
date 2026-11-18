function [lo, info] = lmin_cw(mv, dg, opts)
%LMIN_CW  Rigorous lower bound on lambda_min of a symmetric Z-matrix H, from products with H only.
%   [LO,INFO] = LMIN_CW(MV,DG,OPTS).  H = D - B with B >= 0 entrywise, D >= 0 diagonal.
%   Collatz-Wielandt: for any z > 0,  lambda_min(H) >= min_i (H z)_i / z_i.
%   z1 ~ H^(-1) 1 and z2 ~ H^(-1) z1 are computed by CG with a pointwise stopping test;
%   the larger of the two bounds is returned, lowered by a rounding margin.
%   MV:   @(x) H*x for x of the shape of DG.   DG: the diagonal of H (any array shape).
%   OPTS: .rtol (0.05) pointwise CG tolerance, .maxit (20000)
%   INFO: .products .bound1 .bound2 .time .res1
if nargin < 3, opts = struct(); end
o = struct('rtol', 0.05, 'maxit', 20000);
fn = fieldnames(opts);
for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
t0 = tic;
sz = size(dg);
d = double(dg(:));
H = @(v) reshape(double(mv(reshape(v, sz))), [], 1);
[z1, it1] = cg_pointwise(H, ones(numel(d), 1), o);
[b1, Hz1] = cw_bound(H, d, z1);
[z2, it2] = cg_pointwise(H, max(z1, 0), o);
b2 = cw_bound(H, d, z2);
lo = max(b1, b2);
info = struct('products', it1 + it2 + 2, 'bound1', b1, 'bound2', b2, 'time', toc(t0), 'res1', norm(Hz1 - 1, Inf));
end

function [b, Hz] = cw_bound(H, d, z)
% min_i (H z)_i / z_i, each ratio lowered by 64 eps (|H| z)_i / z_i,  |H| z = 2 D z - H z
Hz = H(z);
if any(z <= 0), b = -Inf;  return, end
b = min((Hz - 64*eps*(2*d.*z - Hz))./z);
end

function [x, it] = cg_pointwise(H, r, o)
% CG for H x = b (x0 = 0, so r = b), stopped once |b - H x| <= rtol |b| in every entry (checked on the recursive residual)
b = abs(r);
x = zeros(size(r));
p = r;
rr = r'*r;
it = 0;
while it < o.maxit
    Ap = H(p);
    it = it + 1;
    al = rr/(p'*Ap);
    x = x + al*p;
    r = r - al*Ap;
    if all(abs(r) <= o.rtol*b), break, end
    rn = r'*r;
    p = r + (rn/rr)*p;
    rr = rn;
end
end
