function P = tb_problem(s)
%TB_PROBLEM  One problem of the basic test case of Example 1:  A X + X L^T = U V^T.
%   P = TB_PROBLEM(S).  Missing fields of the struct S take the defaults in brackets.
%     .dim [2]  .N [256]   large side on the unit square/cube, N^dim nodes, Dirichlet, h = 1/(N+1)
%     .coef ['lap']        'lap' -Laplacian (5/7-point); 'het' -div(a grad), a = 10^(.contrast*field);
%                          'adv' -Laplacian + .w*(d/dx + d/dy [+ d/dz]), central differences
%     .contrast [0]  .w [0]
%     .Ltype ['tri']       small side L: 'diag' 'tri' 'time' 'adv' 'dense' (see below)
%     .n [200]  .f [1]  .g [1e4]   size of L and its spectral range [f, g];  .pe [0] Peclet of 'adv'
%     .rho [2]  .rhs ['rand']      rank of U V^T and kind of U: 'rand' 'smooth' 'point';  .seed [1]
%   L: 'diag' log-spaced in [f, g]; 'tri' f I + (g-f)/4 tridiag(-1,2,-1); 'time' backward Euler,
%   tau = 2/(g-f), plus f I; 'adv' 1-D upwind advection-diffusion shifted by f; 'dense' Q diag Q'.
%   P: .A (sparse) .L .U .V .E = [lower bound lambda_min(A), upper bound lambda_max(A)]
%      .WL = [min max] of Re W(L)  .lamL = eig(L)  .m .n .spec .name
%   The base problem of Example 1 is the default (2-D, N = 256, 'tri', n = 200, [f, g] = [1, 1e4]);
%   Fig. 7(c) uses Ltype 'adv'.  Used by TB_SWEEP (results/testbed.csv, Figs. 6 and 7(c)).
d = struct('dim', 2, 'N', 256, 'coef', 'lap', 'contrast', 0, 'w', 0, 'Ltype', 'tri', ...
           'n', 200, 'f', 1, 'g', 1e4, 'pe', 0, 'rho', 2, 'rhs', 'rand', 'seed', 1);
fn = fieldnames(s);
for i = 1:numel(fn), d.(fn{i}) = s.(fn{i}); end
s = d;
rng(s.seed);
N = s.N;
h = 1/(N + 1);
dim = s.dim;
m = N^dim;
e = ones(N, 1);
D = spdiags([-e e], [-1 0], N + 1, N)/h;          % faces x nodes (Dirichlet: ghost nodes are 0)
I1 = speye(N);
lap1 = D'*D;
lmin1 = 4/h^2*sin(pi*h/2)^2;
lmax1 = 4/h^2*cos(pi*h/2)^2;
switch s.coef
    case {'lap', 'adv'}
        if dim == 2, A = kron(I1, lap1) + kron(lap1, I1);
        else,        A = kron(I1, kron(I1, lap1)) + kron(I1, kron(lap1, I1)) + kron(lap1, kron(I1, I1)); end
        E = dim*[lmin1, lmax1];
        if strcmp(s.coef, 'adv') && s.w ~= 0
            C1 = spdiags([-e e], [-1 1], N, N)/(2*h);            % central first derivative
            if dim == 2, Cv = kron(I1, C1) + kron(C1, I1);
            else,        Cv = kron(I1, kron(I1, C1)) + kron(I1, kron(C1, I1)) + kron(C1, kron(I1, I1)); end
            A = A + s.w*Cv;                                       % sym(A) is still the Laplacian
        end
    case 'het'
        xf = (0:N)'*h + h/2;                                      % face midpoints in the derivative direction
        xn = (1:N)'*h;                                            % node coordinates across
        fld = smooth_field(dim, s.seed);
        if dim == 2
            [Xf, Yn] = ndgrid(xf, xn);
            ax = 10.^(s.contrast*fld(Xf, Yn));                    % x-faces: (N+1) x N
            [Xn, Yf] = ndgrid(xn, xf);
            ay = 10.^(s.contrast*fld(Xn, Yf));                    % y-faces: N x (N+1)
            Dx = kron(I1, D);
            Dy = kron(D, I1);
            A = Dx'*spdiags(ax(:), 0, numel(ax), numel(ax))*Dx + Dy'*spdiags(ay(:), 0, numel(ay), numel(ay))*Dy;
            amin = min([ax(:); ay(:)]);
        else
            [Xf, Yn, Zn] = ndgrid(xf, xn, xn);
            ax = 10.^(s.contrast*fld(Xf, Yn, Zn));
            [Xn, Yf, Zn] = ndgrid(xn, xf, xn);
            ay = 10.^(s.contrast*fld(Xn, Yf, Zn));
            [Xn, Yn, Zf] = ndgrid(xn, xn, xf);
            az = 10.^(s.contrast*fld(Xn, Yn, Zf));
            Dx = kron(I1, kron(I1, D));
            Dy = kron(I1, kron(D, I1));
            Dz = kron(D, kron(I1, I1));
            A = Dx'*spdiags(ax(:), 0, numel(ax), numel(ax))*Dx + Dy'*spdiags(ay(:), 0, numel(ay), numel(ay))*Dy ...
              + Dz'*spdiags(az(:), 0, numel(az), numel(az))*Dz;
            amin = min([ax(:); ay(:); az(:)]);
        end
        E = [amin*dim*lmin1, full(max(sum(abs(A), 2)))];         % min(a)*lambda_min(lap), Gershgorin
end
% ---- small side
n = s.n;
f = s.f;
g = s.g;
switch s.Ltype
    case 'diag'
        lam = logspace(log10(max(f, 1e-12)), log10(g), n)';
        if f == 0, lam(1) = 0; end
        L = spdiags(lam, 0, n, n);
    case 'tri'
        en = ones(n, 1);
        L = f*speye(n) + (g - f)/4*spdiags([-en 2*en -en], -1:1, n, n);
    case 'time'
        tau = 2/(g - f);
        en = ones(n, 1);
        L = spdiags([-en en], [-1 0], n, n)/tau + f*speye(n);
    case 'adv'
        en = ones(n, 1);
        hn = 1/(n + 1);
        nu = 1;
        v = s.pe*2*nu/1;                                          % Peclet over the domain
        K = nu*spdiags([-en 2*en -en], -1:1, n, n)/hn^2 + v*spdiags([-en en], [-1 0], n, n)/hn;   % upwind
        sc = (g - f)/(4*nu/hn^2 + 2*v/hn);
        L = f*speye(n) + sc*K;
    case 'dense'
        lam = logspace(log10(max(f, 1e-12)), log10(g), n)';
        [Q, ~] = qr(randn(n));
        L = Q*diag(lam)*Q';
        L = (L + L')/2;
end
Lf = full(L);
WL = eig((Lf + Lf')/2);
WL = [min(WL), max(WL)];
lamL = eig(Lf);
% ---- right-hand side
rho = s.rho;
switch s.rhs
    case 'rand',   U = randn(m, rho);
    case 'smooth'
        U = zeros(m, rho);
        xn = (1:N)'*h;
        if dim == 2, [X1, X2] = ndgrid(xn, xn); X3 = 0*X1; else, [X1, X2, X3] = ndgrid(xn, xn, xn); end
        for j = 1:rho
            c = 0.2 + 0.6*rand(1, 3);
            r = 0.15;
            u = exp(-((X1 - c(1)).^2 + (X2 - c(2)).^2 + (dim == 3)*(X3 - c(3)).^2)/(2*r^2));
            U(:, j) = u(:);
        end
    case 'point'
        U = zeros(m, rho);
        for j = 1:rho, U(randi(m, 3, 1), j) = 1; end
end
V = randn(n, rho);
U = U/norm(U, 'fro');
V = V/norm(V, 'fro');
P = struct('A', A, 'L', L, 'U', U, 'V', V, 'E', E, 'WL', WL, 'lamL', lamL, 'm', m, 'n', n, 'spec', s);
P.name = sprintf('%dD N=%d %s | L=%s n=%d [%g,%g] | rho=%d %s', dim, N, s.coef, s.Ltype, n, f, g, rho, s.rhs);
end

function fld = smooth_field(dim, seed)
% a smooth random function on the unit square/cube with values in [0, 1]
rs = RandStream('mt19937ar', 'Seed', 1000 + seed);
k = 6;
c = rand(rs, k, 3);
w = 0.12 + 0.2*rand(rs, k, 1);
a = randn(rs, k, 1);
raw = @(x, y, z) sum_g(x, y, z, c, w, a);
% normalise with samples
[xs, ys, zs] = ndgrid(linspace(0, 1, 25));
if dim == 2, zs = 0.5*ones(size(xs)); end
v = raw(xs, ys, zs);
lo = min(v(:));
hi = max(v(:));
if dim == 2, fld = @(x, y) min(1, max(0, (raw(x, y, 0.5*ones(size(x))) - lo)/(hi - lo)));
else,        fld = @(x, y, z) min(1, max(0, (raw(x, y, z) - lo)/(hi - lo))); end
end

function v = sum_g(x, y, z, c, w, a)
v = zeros(size(x));
for i = 1:numel(w)
    v = v + a(i)*exp(-((x - c(i, 1)).^2 + (y - c(i, 2)).^2 + (z - c(i, 3)).^2)/(2*w(i)^2));
end
end
