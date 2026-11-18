function [Z, Y] = gfadi_simple(A, B, U, V, alpha, beta)
%GFADI_SIMPLE  Minimal reference version of g-fADI (alternating order, no options).
%   [Z,Y] = GFADI_SIMPLE(A,B,U,V,ALPHA,BETA) solves A*X - X*B = U*V' approximately,
%   X = Z*Y'.  ALPHA: zeros (right steps: product with A, shifted solve with B),
%   BETA: poles (left steps: shifted solve with A, product with B).  A, B: matrices.
%   See GFADI for the full version (operators, schedules, counters).
Zt = U;  Yt = V;  Z = [];  Y = [];
IA = speye(size(A));  IB = eye(size(B));
p = numel(alpha);  q = numel(beta);
for j = 1:max(p, q)
    if j <= p && j <= q                     % alternating step
        a = alpha(j);  b = beta(j);  g = b - a;
        Zj = g*((A - b*IA)\Zt);   Yj = (B' - conj(a)*IB)\Yt;
        Zt = (A - a*IA)*Zj/g;     Yt = (B' - conj(b)*IB)*Yj;
    elseif j <= q                           % left step
        b = beta(j);
        Zj = (A - b*IA)\Zt;       Yj = Yt;
        Zt = Zj;                  Yt = (B' - conj(b)*IB)*Yj;
    else                                    % right step
        a = alpha(j);
        Zj = -Zt;                 Yj = (B' - conj(a)*IB)\Yt;
        Zt = -(A - a*IA)*Zj;      Yt = Yj;
    end
    s = norm(Zt, 'fro');  Zt = Zt/s;  Yt = Yt*s;   % rescale: the residual Zt*Yt' is unchanged
    Z = [Z, Zj];  Y = [Y, Yj];              %#ok<AGROW>
end
end
