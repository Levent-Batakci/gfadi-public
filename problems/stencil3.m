function op = stencil3(nx, ny, nz, spec)
%STENCIL3  Matrix-free 7-point operator on an nx x ny x nz grid, applied over z-slabs on a thread pool.
%   op = stencil3(nx, ny, nz, spec)
%   y(i,j,k) = d(i,j,k) x(i,j,k) - sum over the six neighbours of a_face x_neighbour.
%   spec, one of
%     constant:  .cx .cy .cz (coefficient/h^2 per direction), .d0 (added to the diagonal, 0), Dirichlet
%     variable:  .Ax (nx-1 x ny x nz) .Ay .Az (face coefficients / h^2), .d (the full diagonal)
%     cells:     .c (uint8 level per cell), .tab (coefficient/h^2 per level), .f (reaction, 0)
%   op: .mv(x, a, g) = g*(A - a I) x for x of size [nx ny nz] or [nx ny nz rho] (a = 0, g = 1),
%       .lmax (Gershgorin bound), .m, .sz, .diag0 (constant case)
%   Arrays keep the class of the coefficients (single or double); x must match.
%   Example: op = stencil3(8, 8, 8, struct('cx', 81, 'cy', 81, 'cz', 81));  y = op.mv(ones(8, 8, 8));
if isempty(gcp('nocreate')), parpool('Threads', maxNumCompThreads); end
op.sz = [nx ny nz];
op.m = nx*ny*nz;
if isfield(spec, 'cx')
    d = struct('dirichlet', true, 'd0', 0);
    fn = fieldnames(spec);
    for i = 1:numel(fn), d.(fn{i}) = spec.(fn{i}); end
    cx = d.cx;
    cy = d.cy;
    cz = d.cz;
    d0 = 2*(cx + cy + cz) + d.d0;
    if ~d.dirichlet, error('stencil3:bc', 'constant coefficients: only Dirichlet here'); end
    op.mv = @(x, varargin) mv_const(x, cx, cy, cz, d0, varargin{:});
    op.lmax = 2*d0 - d.d0;
    op.diag0 = d0;
elseif isfield(spec, 'c')
    c = spec.c;
    tab = spec.tab;
    f0 = 0;
    if isfield(spec, 'f'), f0 = spec.f; end
    op.mv = @(x, varargin) mv_cell(x, c, tab, f0, varargin{:});
    cmax = 0;
    parfor k = 1:nz, cmax = max(cmax, max(reshape(c(:, :, k), [], 1))); end
    op.lmax = 12*double(tab(double(cmax) + 1)) + f0;
else
    Ax = spec.Ax;
    Ay = spec.Ay;
    Az = spec.Az;
    dg = spec.d;
    op.mv = @(x, varargin) mv_var(x, Ax, Ay, Az, dg, varargin{:});
    op.lmax = gersh(Ax, Ay, Az, dg);       % the diagonal plus the absolute off-diagonals of a row
end
end

function y = mv_const(x, cx, cy, cz, d0, a, g)
if nargin < 6, a = 0; end
if nargin < 7, g = 1; end
[n1, n2, n3, nr] = size(x);
if nr > 1, y = zeros(n1, n2, n3, nr, 'like', x); end
dd = cast(g*(d0 - a), 'like', x);
gx = cast(g*cx, 'like', x);
gy = cast(g*cy, 'like', x);
gz = cast(g*cz, 'like', x);
for r = 1:nr
    if nr == 1, xr = x; else, xr = x(:, :, :, r); end
    yr = zeros(n1, n2, n3, 'like', x);
    parfor k = 1:n3
        xk = xr(:, :, k);
        s = dd*xk;
        s(2:end, :) = s(2:end, :) - gx*xk(1:end-1, :);
        s(1:end-1, :) = s(1:end-1, :) - gx*xk(2:end, :);
        s(:, 2:end) = s(:, 2:end) - gy*xk(:, 1:end-1);
        s(:, 1:end-1) = s(:, 1:end-1) - gy*xk(:, 2:end);
        if k > 1, s = s - gz*xr(:, :, k - 1); end   %#ok<PFBNS>
        if k < n3, s = s - gz*xr(:, :, k + 1); end
        yr(:, :, k) = s;
    end
    if nr == 1, y = yr; else, y(:, :, :, r) = yr; end
end
end

function y = mv_var(x, Ax, Ay, Az, dg, a, g)
if nargin < 6, a = 0; end
if nargin < 7, g = 1; end
[n1, n2, n3, nr] = size(x);
if nr > 1, y = zeros(n1, n2, n3, nr, 'like', x); end
ac = cast(a, 'like', x);
gc = cast(g, 'like', x);
for r = 1:nr
    if nr == 1, xr = x; else, xr = x(:, :, :, r); end
    yr = zeros(n1, n2, n3, 'like', x);
    parfor k = 1:n3
        xk = xr(:, :, k);
        s = (dg(:, :, k) - ac).*xk;
        Axk = Ax(:, :, k);
        Ayk = Ay(:, :, k);
        s(1:end-1, :) = s(1:end-1, :) - Axk.*xk(2:end, :);
        s(2:end, :) = s(2:end, :) - Axk.*xk(1:end-1, :);
        s(:, 1:end-1) = s(:, 1:end-1) - Ayk.*xk(:, 2:end);
        s(:, 2:end) = s(:, 2:end) - Ayk.*xk(:, 1:end-1);
        if k > 1, s = s - Az(:, :, k - 1).*xr(:, :, k - 1); end         %#ok<PFBNS>
        if k < n3, s = s - Az(:, :, k).*xr(:, :, k + 1); end   %#ok<PFBNS>
        yr(:, :, k) = gc*s;
    end
    if nr == 1, y = yr; else, y(:, :, :, r) = yr; end
end
end

function y = mv_cell(x, c, tab, f0, a, g)
% flux form: faces carry the mean of their two cells; Dirichlet (half a cell away) on the box
if nargin < 5, a = 0; end
if nargin < 6, g = 1; end
[n1, n2, n3] = size(x);
y = zeros(n1, n2, n3, 'like', x);
dr = cast(f0 - a, 'like', x);
gc = cast(g, 'like', x);
tb = cast(tab(:), 'like', x);
parfor k = 1:n3
    xk = x(:, :, k);
    ak = tb(uint16(c(:, :, k)) + 1);                                     %#ok<PFBNS>
    ak = reshape(ak, n1, n2);
    q = 0.5*(ak(1:end-1, :) + ak(2:end, :)).*(xk(2:end, :) - xk(1:end-1, :));
    s = dr*xk;
    s(1:end-1, :) = s(1:end-1, :) - q;
    s(2:end, :) = s(2:end, :) + q;
    q = 0.5*(ak(:, 1:end-1) + ak(:, 2:end)).*(xk(:, 2:end) - xk(:, 1:end-1));
    s(:, 1:end-1) = s(:, 1:end-1) - q;
    s(:, 2:end) = s(:, 2:end) + q;
    s(1, :) = s(1, :) + 2*ak(1, :).*xk(1, :);
    s(n1, :) = s(n1, :) + 2*ak(n1, :).*xk(n1, :);
    s(:, 1) = s(:, 1) + 2*ak(:, 1).*xk(:, 1);
    s(:, n2) = s(:, n2) + 2*ak(:, n2).*xk(:, n2);
    if k > 1
        an = reshape(tb(uint16(c(:, :, k - 1)) + 1), n1, n2);   %#ok<PFBNS>
        s = s + 0.5*(ak + an).*(xk - x(:, :, k - 1));   %#ok<PFBNS>
    else
        s = s + 2*ak.*xk;
    end
    if k < n3
        an = reshape(tb(uint16(c(:, :, k + 1)) + 1), n1, n2);
        s = s + 0.5*(ak + an).*(xk - x(:, :, k + 1));
    else
        s = s + 2*ak.*xk;
    end
    y(:, :, k) = gc*s;
end
end

function l = gersh(Ax, Ay, Az, dg)
[n1, n2, n3] = size(dg);
l = 0;
parfor k = 1:n3
    s = abs(dg(:, :, k));
    Axk = abs(Ax(:, :, k));
    Ayk = abs(Ay(:, :, k));
    s(1:end-1, :) = s(1:end-1, :) + Axk;
    s(2:end, :) = s(2:end, :) + Axk;
    s(:, 1:end-1) = s(:, 1:end-1) + Ayk;
    s(:, 2:end) = s(:, 2:end) + Ayk;
    if k > 1, s = s + abs(Az(:, :, k - 1)); end
    if k < n3, s = s + abs(Az(:, :, k)); end   %#ok<PFBNS>
    l = max(l, double(max(s(:))));
end
end
