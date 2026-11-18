function [QA, QB, Y, info] = sylvester_galerkin(A, B, U, V, polesA, polesB, varargin)
%SYLVESTER_GALERKIN  (Rational) Krylov Galerkin projection for A X - X B = U V', X ~ QA*Y*QB'.
%   [QA, QB, Y, info] = sylvester_galerkin(A, B, U, V, polesA, polesB, 'tol', tol, 'residual', true)
%   Bases range(QA) = K(A, U; polesA), range(QB) = K(B', V; polesB), one block per iteration
%   (polesA, polesB of equal length k; Inf = a plain product, finite s = a shifted solve).
%   All poles Inf gives PKSM; Inf on A and the g-fADI right shifts on B spans the spaces of
%   one-sided g-fADI; Zolotarev/adaptive poles on both sides give RKSM.
%   A, B: matrices or operator structs (see GFADI_OP).  The Galerkin method of Example 2
%   (EXP16_BLOWOUT, Fig. 8).
%   'tol': stop when the exact ||R||_F/||U V'||_F < tol (default 0: run all k iterations)
%   'residual': compute the residual even when tol = 0 (default true)
%   info: .res residual history, .k iterations, .time, operation counts .nsolA .nmvA .nsolB .nmvB
%   Example: [QA, QB, Y] = sylvester_galerkin(A, B, U, V, inf(1, 30), shifts, 'tol', 1e-8);
p = inputParser;
p.addParameter('tol', 0);
p.addParameter('residual', true);
p.parse(varargin{:});
tol = p.Results.tol;
wantres = p.Results.residual;
if ~isstruct(A), A = gfadi_op(A); end
if ~isstruct(B), B = gfadi_op(B); end
k = numel(polesA);
rho = size(U, 2);
m = size(U, 1);
n = size(V, 1);
nF = sqrt(abs(trace((U'*U)*(V'*V))));
cA = min(m, rho*k);
cB = min(n, rho*k);
QA = zeros(m, cA);
AQ = zeros(m, cA);
QB = zeros(n, cB);
BQ = zeros(n, cB);
cplx = ~isreal(U) || ~isreal(V) || any(imag(polesA(isfinite(polesA)))) || any(imag(polesB(isfinite(polesB))));
if cplx
    QA = complex(QA);
    QB = complex(QB);
    AQ = complex(AQ);
    BQ = complex(BQ);
end
cnt = struct('sA', 0, 'mA', 0, 'sB', 0, 'mB', 0);
tstart = tic;
WA = U;
WB = V;
info.res = [];
rA = 0;
rB = 0;
doneA = false;          % a side stops growing when its Krylov space is exhausted
doneB = false;
for j = 1:k
    if j > 1
        if ~doneA
            if isfinite(polesA(j-1)), WA = A.sv(polesA(j-1), WA);  cnt.sA = cnt.sA + rho;
            else,                      WA = A.mv(WA);               cnt.mA = cnt.mA + rho; end
        end
        if ~doneB
            if isfinite(polesB(j-1)), WB = B.svh(polesB(j-1), WB); cnt.sB = cnt.sB + rho;
            else,                      WB = B.mvh(WB);              cnt.mB = cnt.mB + rho; end
        end
    end
    if ~doneA
        [qa, na] = orthblock(WA, QA);           % zero-padded: no block copy
        if na < 1e-10 || rA + rho > cA
            doneA = true;
        else
            QA(:, rA+1:rA+rho) = qa;
            AQ(:, rA+1:rA+rho) = A.mv(qa);
            cnt.mA = cnt.mA + rho;
            rA = rA + rho;
            WA = qa;                % continuation block: keeps the recurrence well scaled
        end
    end
    if ~doneB
        [qb, nb] = orthblock(WB, QB);
        if nb < 1e-10 || rB + rho > cB
            doneB = true;
        else
            QB(:, rB+1:rB+rho) = qb;
            BQ(:, rB+1:rB+rho) = B.mvh(qb);
            cnt.mB = cnt.mB + rho;
            rB = rB + rho;
            WB = qb;
        end
    end
    % ---- projected equation and exact residual (only when needed)
    if tol > 0 || j == k || (doneA && doneB)
        Ar = QA(:, 1:rA)'*AQ(:, 1:rA);
        Bproj = BQ(:, 1:rB)'*QB(:, 1:rB);              % BQ = B'*QB, so this is QB'*B*QB
        Cr = (QA(:, 1:rA)'*U)*(V'*QB(:, 1:rB));
        Yj = sylvester(Ar, -Bproj, Cr);
        res = NaN;
        if wantres || tol > 0
            % R = AQ*Y*QB' - QA*Y*(BQ)' - U*V'  =  [AQ QA U] * [Y*QB'; -Y*BQ'; -V']
            L  = [AQ(:, 1:rA), QA(:, 1:rA), U];
            Rr = [QB(:, 1:rB)*Yj', -BQ(:, 1:rB)*Yj', -V];
            [~, RL] = qr(L, 0);
            [~, RR] = qr(Rr, 0);
            res = norm(RL*RR', 'fro')/nF;
        end
        info.res(end+1) = res;                                            %#ok<AGROW>
        if (tol > 0 && res < tol) || (doneA && doneB), break, end
    end
end
QA = QA(:, 1:rA);
QB = QB(:, 1:rB);
Y = Yj;
info.k = j;
info.time = toc(tstart);
info.nsolA = cnt.sA;
info.nmvA = cnt.mA;
info.nsolB = cnt.sB;
info.nmvB = cnt.mB;
end

function [q, nrm] = orthblock(W, Q)
%ORTHBLOCK  Orthonormalise the block W against Q (classical Gram-Schmidt twice).
n0 = norm(W, 'fro');
for rep = 1:2
    if ~isempty(Q), W = W - Q*(Q'*W); end   % Q may carry zero columns
end
nrm = norm(W, 'fro')/max(n0, realmin);
[q, ~] = qr(W, 0);
end
