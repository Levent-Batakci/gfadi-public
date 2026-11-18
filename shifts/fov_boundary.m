function [z, info] = fov_boundary(B, K)
%FOV_BOUNDARY  Points on the boundary of the field of values W(B) (Johnson's method).
%   [Z,INFO] = FOV_BOUNDARY(B,K) returns K (default 720) boundary points: for each angle
%   theta, the top eigenvector x of the Hermitian part of exp(i*theta)*B gives z = x'*B*x.
%   Their polygon is an inner approximation of W(B).  INFO.theta, INFO.lmax: the support
%   function, W(B) lies in the half-planes Re(exp(i*theta_k)*w) <= lmax(k) (outer bound).
if nargin < 2, K = 720; end
B = full(B);
th = 2*pi*(0:K-1)/K;
z = zeros(K, 1);  lmax = zeros(K, 1);
for k = 1:K
    M = exp(1i*th(k))*B;
    H = (M + M')/2;
    [Q, D] = eig(H);
    [lmax(k), j] = max(real(diag(D)));
    x = Q(:, j);
    z(k) = x'*B*x;
end
info.theta = th(:);  info.lmax = lmax;
end
