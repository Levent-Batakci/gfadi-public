function P = build_age3d(N, n, opts)
%BUILD_AGE3D  Structured-population (age/energy) diffusion in a heterogeneous 3-D medium (Example 2).
%   P = build_age3d(N, n, opts)
%   PDE: -div(a(x) grad u) + v u_t - eps u_tt + mu u = f(x) g(t) on (0,1)^3 x (0,1), Dirichlet BCs.
%   Grid: 7-point FD on N^3 interior nodes in x (h = 1/(N+1)), centred FD on n nodes in t.
%   With X(:,j) = u(., t_j):  A X + X L^T = U V^T,  A sparse SPD (N^3 x N^3), L nonnormal (n x n),
%   i.e. A X - X B = U V^T with B = -L^T.  The paper's Example 2 (Figs. 7a-b and 8) uses n = 400
%   and the default opts.
%   opts: contrast (10, a in [1, contrast]), eps (0.02), v (1), mu (1), afun/amin (custom a >= amin),
%         nodal (false: a at the nodes, each face takes the smaller value; amin may then be 0)
%   P: .A .L .U .V .m .n .N .opts, and spectral enclosures
%      lminA <= lambda(A) <= lmaxA,  Re W(L) in [reLmin, reLmax];  .anode when nodal
%   Example: P = build_age3d(12, 40);
if nargin < 3, opts = struct(); end
d = struct('contrast', 10, 'eps', 0.02, 'v', 1, 'mu', 1, 'afun', [], 'amin', 1, 'nodal', false);
fn = fieldnames(opts);
for i = 1:numel(fn), d.(fn{i}) = opts.(fn{i}); end
h = 1/(N + 1);
m = N^3;
afun = @(x, y, z) d.contrast.^((1 + sin(2*pi*x + 1).*sin(3*pi*y).*sin(2*pi*z + 0.5))/2);
if ~isempty(d.afun), afun = d.afun; end                  % custom coefficient, a >= d.amin
g = (1:N)*h;                                            % nodes
gf = ((0:N) + 0.5)*h;                                   % faces
diagv = zeros(N, N, N);
I = [];
J = [];
S = [];
idx = reshape(1:m, N, N, N);
if d.nodal
    [X, Y, Z] = ndgrid(g, g, g);
    an = afun(X, Y, Z);
    P.anode = an(:);
end
for dirn = 1:3
    switch dirn                                          % face coefficients, (N+1) faces along dirn
        case 1, [X, Y, Z] = ndgrid(gf, g, g);
        case 2, [X, Y, Z] = ndgrid(g, gf, g);
        case 3, [X, Y, Z] = ndgrid(g, g, gf);
    end
    if d.nodal
        p0 = permute(an, [dirn, setdiff(1:3, dirn)]);   % faces along the first dimension
        af = min([p0(1, :, :); p0], [p0; p0(end, :, :)]);
        af = ipermute(af, [dirn, setdiff(1:3, dirn)])/h^2;
    else
        af = afun(X, Y, Z)/h^2;
    end
    lo = cell(1, 3);
    hi = lo;
    for k = 1:3, lo{k} = ':'; hi{k} = ':'; end
    lo{dirn} = 1:N;
    hi{dirn} = 2:N+1;
    diagv = diagv + af(lo{:}) + af(hi{:});
    in = lo;
    in{dirn} = 2:N;                                      % interior faces between nodes i-1 and i
    a1 = idx;
    a2 = idx;
    s1 = {':', ':', ':'};
    s2 = s1;
    s1{dirn} = 1:N-1;
    s2{dirn} = 2:N;
    a1 = a1(s1{:});
    a2 = a2(s2{:});
    w = af(in{:});
    I = [I; a1(:); a2(:)];                               %#ok<AGROW>
    J = [J; a2(:); a1(:)];                               %#ok<AGROW>
    S = [S; -w(:); -w(:)];                               %#ok<AGROW>
end
A = sparse([I; (1:m)'], [J; (1:m)'], [S; diagv(:)], m, m);
ht = 1/(n + 1);
e = ones(n, 1);
T2 = spdiags([-e 2*e -e], -1:1, n, n)/ht^2;
D1 = spdiags([-e 0*e e], -1:1, n, n)/(2*ht);
L = d.eps*T2 + d.v*D1 + d.mu*speye(n);
P.A = A;
P.L = L;
P.m = m;
P.n = n;
P.N = N;
P.opts = d;
P.lminA = d.amin*3*(2 - 2*cos(pi*h))/h^2;              % A >= min(a) * (-Laplacian)
P.lmaxA = full(max(sum(abs(A), 2)));                   % Gershgorin
P.reLmin = d.mu + d.eps*(2 - 2*cos(pi*ht))/ht^2;       % the advection part is skew-symmetric
P.reLmax = d.mu + d.eps*(2 - 2*cos(pi*n*ht))/ht^2;
[X3, Y3, Z3] = ndgrid(g, g, g);
src = @(c, w) exp(-((X3 - c(1)).^2 + (Y3 - c(2)).^2 + (Z3 - c(3)).^2)/w^2);
P.U = [reshape(src([0.3 0.4 0.5], 0.1), [], 1), reshape(src([0.7 0.6 0.4], 0.15), [], 1)];
tt = (1:n)'*ht;
P.V = [exp(-((tt - 0.1)/0.05).^2), exp(-((tt - 0.5)/0.2).^2)];
end
