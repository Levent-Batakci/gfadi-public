function P = switch_problem(N, opts)
%SWITCH_PROBLEM  Markov-switching diffusion (intermittent search) in the unit cube (Example 3).
%   P = switch_problem(N, opts)
%   PDE: D_k Lap u_k + sum_j Q_kj (u_j - u_k) - delta u_k = 0, u = 1 on a spherical target,
%   no flux on the walls; u_k(x) = probability of reaching the target from x in state k.
%   Grid: cell-centred FD on N^3 cells, target cells removed (Dirichlet faces); n = 3*nc states
%   (a driven cycle of nc phases times free/sliding/bound).  After the scaling of Proposition A.1
%   of the paper, A Y + Y L.' = U V.',  A = -Lap_h (matrix-free op), U = b, and u_k = Y(:,k)*back(k).
%   The paper's Example 3 (Fig. 9) uses the defaults with .swirl = 1.
%   opts: .nc (40) .drive (20) .delta (0.2) .Dmode ([1 0.03 1e-3]) .R (0.1) .c ([0.7 0.5 0.5])
%         .tol (1e-8) .class ('double') .sparse (false; true also returns P.As) .swirl (0: drift strength)
%         .Ebottom (unset: eigs at N <= 32, else from the N = 32 problem; 'cw': rigorous lower
%         spectral bound by LMIN_CW); rates .kcyc .kslow .kon .koff .kbind .kunbind
%   P: .op (.mv, .lmax) .U .V .L .Bop .p .al (shifts) .E (spectral interval of A) .f .WL (Re W(L))
%      .S .Dk .G .Q .pi .T (target mask) .b .back .N .n .nc .m .sz .xc .h .o .As (if sparse)
%   Example: P = switch_problem(16, struct('nc', 8));
d = struct('nc', 40, 'drive', 20, 'delta', 0.2, 'Dmode', [1 0.03 1e-3], 'R', 0.1, 'c', [0.7 0.5 0.5], 'tol', 1e-8, ...
           'class', 'double', 'sparse', false, 'kcyc', 60, 'kslow', 4, 'kon', 20, 'koff', 5, 'kbind', 8, 'kunbind', 2, 'swirl', 0);
if nargin > 1
    fn = fieldnames(opts);
    for i = 1:numel(fn), d.(fn{i}) = opts.(fn{i}); end
end
o = d;
cls = o.class;
h = 1/N;
xc = ((1:N) - 0.5)*h;
% ---- the large side: -Lap on the cube, reflecting walls, the target removed with Dirichlet faces
[X1, X2, X3] = ndgrid(xc, xc, xc);
T = (X1 - o.c(1)).^2 + (X2 - o.c(2)).^2 + (X3 - o.c(3)).^2 < o.R^2;           % target cells
F = ~T;
ih2 = 1/h^2;
Ax = cast(ih2*(F(1:end-1, :, :) & F(2:end, :, :)), cls);
Ay = cast(ih2*(F(:, 1:end-1, :) & F(:, 2:end, :)), cls);
Az = cast(ih2*(F(:, :, 1:end-1) & F(:, :, 2:end)), cls);
% faces from a free cell to a target cell: Dirichlet u = 1 there (half-cell offset ignored, first order)
Tx = ih2*(F(1:end-1, :, :) & T(2:end, :, :));
Tx2 = ih2*(T(1:end-1, :, :) & F(2:end, :, :));
Ty = ih2*(F(:, 1:end-1, :) & T(:, 2:end, :));
Ty2 = ih2*(T(:, 1:end-1, :) & F(:, 2:end, :));
Tz = ih2*(F(:, :, 1:end-1) & T(:, :, 2:end));
Tz2 = ih2*(T(:, :, 1:end-1) & F(:, :, 2:end));
bnd = zeros(N, N, N);
bnd(1:end-1, :, :) = bnd(1:end-1, :, :) + Tx;
bnd(2:end, :, :) = bnd(2:end, :, :) + Tx2;
bnd(:, 1:end-1, :) = bnd(:, 1:end-1, :) + Ty;
bnd(:, 2:end, :) = bnd(:, 2:end, :) + Ty2;
bnd(:, :, 1:end-1) = bnd(:, :, 1:end-1) + Tz;
bnd(:, :, 2:end) = bnd(:, :, 2:end) + Tz2;
dg = zeros(N, N, N);
ax = double(Ax);
ay = double(Ay);
az = double(Az);
dg(1:end-1, :, :) = dg(1:end-1, :, :) + ax;
dg(2:end, :, :) = dg(2:end, :, :) + ax;
dg(:, 1:end-1, :) = dg(:, 1:end-1, :) + ay;
dg(:, 2:end, :) = dg(:, 2:end, :) + ay;
dg(:, :, 1:end-1) = dg(:, :, 1:end-1) + az;
dg(:, :, 2:end) = dg(:, :, 2:end) + az;
dg = dg + bnd;
dg(T) = 6*ih2;                                       % target cells: decoupled rows in the middle of the spectrum, X = 0 there
op = stencil3(N, N, N, struct('Ax', Ax, 'Ay', Ay, 'Az', Az, 'd', cast(dg, cls)));
b = bnd;
b(T) = 0;                                            % u = 1 on the target enters through the Dirichlet faces
% ---- the small side: phase i = 1..nc on a driven cycle, mode r = 1 free, 2 sliding, 3 bound;  state k = i + (r-1) nc
nc = o.nc;
n = 3*nc;
ph = (0:nc-1)'/nc;
kf = o.kcyc*ones(nc, 1);
kf(round(0.6*nc)) = o.kslow;                         % one rate-limiting step
kb = kf/o.drive;                                     % backward rates: the cycle runs forward
I = [];
J = [];
W = [];
for r = 1:3
    for i = 1:nc
        k = i + (r - 1)*nc;
        inx = mod(i, nc) + 1 + (r - 1)*nc;
        inp = mod(i - 2, nc) + 1 + (r - 1)*nc;
        I = [I; k; k];                               %#ok<AGROW>
        J = [J; inx; inp];                           %#ok<AGROW>
        W = [W; kf(i); kb(i)];                       %#ok<AGROW>
    end
end
% mode switching gated by the phase: sliding starts in the first half of the cycle, binding from sliding in the second
on = o.kon*(ph < 0.5);
off = o.koff*ones(nc, 1);
bind = o.kbind*(ph >= 0.5);
unbind = o.kunbind*ones(nc, 1);
for i = 1:nc
    f1 = i;
    s1 = i + nc;
    t1 = i + 2*nc;
    I = [I; f1; s1; s1; t1];                         %#ok<AGROW>
    J = [J; s1; f1; t1; s1];                         %#ok<AGROW>
    W = [W; on(i); off(i); bind(i); unbind(i)];      %#ok<AGROW>
end
keep = W > 0;
Q = sparse(I(keep), J(keep), W(keep), n, n);
Q = Q - spdiags(full(sum(Q, 2)), 0, n, n);
Dk = repmat(o.Dmode(:)', nc, 1);
Dk = Dk(:);
G = o.delta*speye(n) - Q;                                                    % nonsingular M-matrix
% ---- the scaling of Proposition A.1 on G' = Dd^(-1/2) G Dd^(-1/2)
Gp = spdiags(1./sqrt(Dk), 0, n, n)*G*spdiags(1./sqrt(Dk), 0, n, n);
xs = Gp\ones(n, 1);
ys = Gp.'\ones(n, 1);
S = sqrt(ys./xs);
f = 0.5*min(1./xs + 1./ys);
L = full(spdiags(S, 0, n, n)*Gp*spdiags(1./S, 0, n, n));
ev = eig((L + L')/2);
WL = [min(ev) max(ev)];                                                      % exact Re W(L), to compare with f
% the stationary distribution of the chain (ignoring degradation, which is the same in every state)
[pv, ~] = eigs(Q.', 1, 'smallestabs');
piq = abs(pv)/sum(abs(pv));
% ---- the right-hand side of the scaled equation:  A Y + Y L.' = U V',  U = b,  V = S Dd^(-1/2) d with d = Dk
U = cast(b, cls);
V = S./sqrt(Dk).*Dk;
% ---- E: Gershgorin above; below, eigs at small N, a coarse grid otherwise, or LMIN_CW (opts.Ebottom = 'cw')
tE = 0;
Einfo = [];
if isfield(o, 'Ebottom') && strcmp(o.Ebottom, 'cw')
    [Elo, Einfo] = lmin_cw(@(x) op.mv(x), dg);
    tE = Einfo.time;
    lmin = NaN;
elseif N <= 32
    As = sparse_of(Ax, Ay, Az, dg);
    lmin = eigs(As, 1, 'smallestabs');
    Elo = 0.98*lmin;
else
    Pc = switch_problem(32, setfield(o, 'sparse', false));   %#ok<SFLD>
    Elo = 0.9*Pc.lmin;
    lmin = NaN;
end
E = [Elo, op.lmax];
opH = op;
Sx = [];
shift = [];
tS = tic;
if o.swirl > 0
    % F = W (-(y - 1/2), x - 1/2, 0.3 sin(2 pi x)) / max|.|, at the faces; divergence free.  Faces into walls/target: none.
    xf = (1:N-1)*h;
    [Xf, Yc, ~] = ndgrid(xf, xc, xc);
    [Xc2, Yf, ~] = ndgrid(xc, xf, xc);
    [Xc3, ~, ~] = ndgrid(xc, xc, xf(1:N-1));
    nrm = max(hypot(hypot(0.5, 0.5), 0.3));
    sc = o.swirl/nrm;
    Fx = -sc*(Yc - 0.5).*(F(1:end-1, :, :) & F(2:end, :, :));
    Fy = sc*(Xc2 - 0.5).*(F(:, 1:end-1, :) & F(:, 2:end, :));
    Fz = sc*0.3*sin(2*pi*Xc3).*(F(:, :, 1:end-1) & F(:, :, 2:end));
    Fx = cast(Fx/(2*h), cls);
    Fy = cast(Fy/(2*h), cls);
    Fz = cast(Fz/(2*h), cls);
    % fused 7-point product, z-slabs on the thread pool (as STENCIL3): face coefficients A -/+ F in the two directions
    Axc = cast(Ax, cls);
    Ayc = cast(Ay, cls);
    Azc = cast(Az, cls);
    dgc = cast(dg, cls);
    z1 = zeros(N, N, 1, cls);
    Azp = cat(3, Azc + Fz, z1);                      % slab k uses Azp(:,:,k), Azm(:,:,k)
    Azm = cat(3, z1, Azc - Fz);
    Axp = Axc + Fx;
    Axm = Axc - Fx;
    Ayp = Ayc + Fy;
    Aym = Ayc - Fy;
    op.mv = @(x, varargin) mv_fused(x, Axp, Axm, Ayp, Aym, Azp, Azm, dgc, varargin{:});
    mF = [max(abs(Fx(:))) max(abs(Fy(:))) max(abs(Fz(:)))]*2*h;
    c = norm(mF);
    ymax = sum(mF)/h;
    [p, al, shift] = adv_shifts(E, c, ymax, [f, WL(2)*1.01], o.tol);
    Sx = struct('Fx', Fx, 'Fy', Fy, 'Fz', Fz, 'c', c, 'ymax', ymax);
else
    [p, al] = shift_search_real(E, [-WL(2)*1.01, -f], 0, o.tol);
end
al = leja_order(al, E');
tshift = toc(tS);                                    % includes building the swirl's face arrays
Bop = struct('svh', @(a, Y) -((L + a*eye(n))\Y), 'n', n);
P = struct('op', op, 'U', U, 'V', V, 'L', L, 'Bop', Bop, 'al', al, 'p', p, 'E', E, 'lmin', lmin, 'f', f, 'WL', WL, 'S', S, ...
           'Dk', Dk, 'G', G, 'Q', Q, 'pi', piq, 'T', T, 'b', b, 'N', N, 'n', n, 'nc', nc, 'm', N^3, 'sz', [N N N], 'xc', xc, 'h', h, 'o', o, ...
           'back', S.^-1./sqrt(Dk), 'opH', opH, 'swirl', Sx, 'shift', shift, 'tE', tE, 'Einfo', Einfo, 'tshift', tshift);    % U(:, k) = Y(:, k)*back(k)
if o.sparse
    P.As = sparse_of(Ax, Ay, Az, dg);
    if o.swirl > 0, P.As = P.As + skew_of(Sx.Fx, Sx.Fy, Sx.Fz); end
end
end

function y = mv_swirl(opH, Fx, Fy, Fz, x, a, g) %#ok<DEFNU>
% g*((H + S) - a I) x with (S x)_i = sum over faces of -/+ F_f/(2h) x_neighbour (skew-symmetric); unfused version of MV_FUSED
if nargin < 6, a = 0; end
if nargin < 7, g = 1; end
if ~isreal(a) && isreal(x), x = complex(x); end
y = opH.mv(x, a, g);
y(1:end-1, :, :, :) = y(1:end-1, :, :, :) - g*Fx.*x(2:end, :, :, :);
y(2:end, :, :, :) = y(2:end, :, :, :) + g*Fx.*x(1:end-1, :, :, :);
y(:, 1:end-1, :, :) = y(:, 1:end-1, :, :) - g*Fy.*x(:, 2:end, :, :);
y(:, 2:end, :, :) = y(:, 2:end, :, :) + g*Fy.*x(:, 1:end-1, :, :);
y(:, :, 1:end-1, :) = y(:, :, 1:end-1, :) - g*Fz.*x(:, :, 2:end, :);
y(:, :, 2:end, :) = y(:, :, 2:end, :) + g*Fz.*x(:, :, 1:end-1, :);
end

function y = mv_fused(x, Axp, Axm, Ayp, Aym, Azp, Azm, dg, a, g)
% g*((H + S) - a I) x:  y_i = d_i x_i - (A + F)_f x_{i+1} - (A - F)_f x_{i-1} in each direction, slab by slab in z
if nargin < 9, a = 0; end
if nargin < 10, g = 1; end
cplx = ~isreal(x) || ~isreal(a);
[n1, n2, n3, nr] = size(x);
y = zeros(n1, n2, n3, nr, 'like', x);
if cplx, y = complex(y); end
for r = 1:nr
    xr = x(:, :, :, r);
    if cplx, yr = complex(zeros(n1, n2, n3, 'like', real(x))); else, yr = zeros(n1, n2, n3, 'like', x); end
    parfor k = 1:n3
        xk = xr(:, :, k);
        s = (dg(:, :, k) - a).*xk;
        s(1:end-1, :) = s(1:end-1, :) - Axp(:, :, k).*xk(2:end, :);
        s(2:end, :) = s(2:end, :) - Axm(:, :, k).*xk(1:end-1, :);
        s(:, 1:end-1) = s(:, 1:end-1) - Ayp(:, :, k).*xk(:, 2:end);
        s(:, 2:end) = s(:, 2:end) - Aym(:, :, k).*xk(:, 1:end-1);
        if k > 1, s = s - Azm(:, :, k).*xr(:, :, k - 1); end                                           %#ok<PFBNS>
        if k < n3, s = s - Azp(:, :, k).*xr(:, :, k + 1); end
        yr(:, :, k) = g*s;
    end
    y(:, :, :, r) = yr;
end
end

function S = skew_of(Fx, Fy, Fz)
% the sparse skew-symmetric drift part (as in MV_FUSED)
[n1, n2, n3] = size(Fx);
n1 = n1 + 1;
M = n1*n2*n3;
id = reshape(1:M, n1, n2, n3);
i1 = id(1:n1-1, :, :);
i2 = id(2:n1, :, :);
j1 = id(:, 1:n2-1, :);
j2 = id(:, 2:n2, :);
k1 = id(:, :, 1:n3-1);
k2 = id(:, :, 2:n3);
fx = double(Fx(:));
fy = double(Fy(:));
fz = double(Fz(:));
S = sparse([i1(:); i2(:); j1(:); j2(:); k1(:); k2(:)], [i2(:); i1(:); j2(:); j1(:); k2(:); k1(:)], [-fx; fx; -fy; fy; -fz; fz], M, M);
end

function A = sparse_of(Ax, Ay, Az, dg)
% the sparse matrix of the 7-point operator with face coefficients Ax, Ay, Az and diagonal dg
[n1, n2, n3] = size(dg);
M = n1*n2*n3;
id = reshape(1:M, n1, n2, n3);
i1 = id(1:n1-1, :, :);
i2 = id(2:n1, :, :);
j1 = id(:, 1:n2-1, :);
j2 = id(:, 2:n2, :);
k1 = id(:, :, 1:n3-1);
k2 = id(:, :, 2:n3);
ax = double(Ax(:));
ay = double(Ay(:));
az = double(Az(:));
A = sparse([i1(:); i2(:); j1(:); j2(:); k1(:); k2(:); id(:)], [i2(:); i1(:); j2(:); j1(:); k2(:); k1(:); id(:)], [-ax; -ax; -ay; -ay; -az; -az; double(dg(:))], M, M);
end
