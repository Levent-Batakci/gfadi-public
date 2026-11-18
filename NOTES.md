# NOTES: file-by-file detail

This is the technical companion to the short headers in the code. The headers say what each file does, how to call it and what it returns or writes. This file collects what did not fit there: derivations, design choices, pitfalls, and how the numbers behind each figure are produced. See `README.md` for how to run things.

Repository-wide conventions:
- Paths: `setup_paths` adds every folder; `app_root()` is the repository root. Runs write to `results/`, figures to `figs/`, drivers' logs to `logs/`. None of these is committed.
- Steps: a right step uses a zero `alpha` (product with A, solve with B); a left step uses a pole `beta` (shifted solve with A, product with B). Type (p,q) means p zeros (right steps) and q poles (left steps). Most problems are written `A X + X L' = U V'`, that is `B = -L'`.
- Figures of the paper: Fig. 1 model problem (teaser), Fig. 2 Zolotarev gallery, Fig. 3 shift ordering, Fig. 4 Zolotarev rationals on two intervals, Fig. 5 efficiency ratio and asymmetric intervals, Fig. 6 Example 1 test bed, Fig. 7 bound and nonnormality, Fig. 8 Example 2 cost, Fig. 9 Example 3 switching diffusion, Fig. 10 Example 4 frequency sweep.

## Contents
1. Core library: `gfadi/`, `shifts/`, `zolo/`
2. Competitors and problem builders: `solvers/`, `problems/`
3. Runs behind the figures: `experiments/`, `teaser/`
4. The test bed: `testbed/`
5. Figure scripts: `figures/`

---

## 1. Core library

### gfadi/gfadi.m
- Step types. Left step (pole beta_j): shifted solve with A, product with B. Right step (zero alpha_l): product with A, solve with B. After p right and q left steps from X0 = 0,
  X_{p+q} - X = r(A) (X0 - X) r(B)^{-1},  r(z) = prod(z - alpha)/prod(z - beta),
  so p = numel(alpha) is the number of zeros and q = numel(beta) the number of poles.
- Operator structs: `.mv(X) = M*X`, `.mvh(X) = M'*X`, `.sv(s,X) = (M - s*I)\X`, `.svh(s,X) = (M - s*I)'\X` (build them with `gfadi_op`).
- Orders. `'alt'` (default): alternate left and right steps while both lists last, then finish with the longer list. It uses the compressed (alternating) recurrences for Z and Y of Section 2.4 of the paper, so rank(X) <= rho*max(p,q). `'leftfirst'` / `'rightfirst'`: all left (right) steps first, no compression, so rank(X) <= rho*(p+q). In exact arithmetic all orders give the same result; in double precision they do not (Fig. 3). `'sched'` (e.g. 'RRLRRL'): shifts are used in the order given, and adjacent L,R or R,L pairs are merged into one compressed step.
- Rescaling of the residual factors. R = Zt*Yt' does not change under Zt -> Zt/c, Yt -> c*Yt. Without the rescaling, the product of (A - alpha_j) over a long one-sided run overflows (|A - alpha| ~ 1e4 per step, hundreds of steps), and the compensating factors (B - alpha)^{-1} in Yt underflow. Zt is kept at unit Frobenius norm.
- The counters `.nsolA .nmvA .nsolB .nmvB` count right-hand-side columns (rho per step), not calls.

### gfadi/gfadi_lowrank.m
- Same iteration as GFADI, with compression on the small side. The increments Zn*Yn' are accumulated in a basis of the small side: X = Zacc*Qy', with Qy (n x s) orthonormal. Each new block Yn is projected onto Qy. Only the part Yp that is not yet represented enlarges the basis, and only if it matters for X.
- `'ytolmode' = 'solution'` (default): a direction is added when its contribution to the increment, ||Zn|| * (its size in Yp), exceeds ytol*||X|| (||X|| as of the last flush of the accumulator). Smaller remainders are dropped and their sizes are summed in INFO.dropped, which bounds ||X_computed - X_iterate||_F / ||X||_F. This keeps s at the numerical rank also when late increments are tiny but point in all directions (strongly nonnormal or badly scaled B).
- `'ytolmode' = 'increment'`: a direction is added when its size exceeds ytol*||Yn||, relative to the increment itself, whatever its size.
- Cost: the update is one m x rho x s product per step, and nothing else touches the long dimension. Memory is m*s instead of m*rho*(p+q). A one-sided run on a 2e6-by-400 Sylvester equation would otherwise need about 80 GB.
- Why an SVD for the new directions: an unpivoted QR would drop R(1,2)-type components that are parallel to column 1.
- Why the new directions are reorthogonalized against Qs: when Yp is the small remainder of a strongly cancelled block, its orthogonal factor is orthogonal to Qs only to eps*||Yn||/|R_ii|. Without the second pass the basis loses orthogonality and then grows without bound.
- Buffering (NBUF = 50): the m x s accumulator is memory-bound, and touching it at every step would cost more than the products with A.
- The options `'trunc'`, `'every'` and `'verbose'` are accepted but not used.

### gfadi/gfadi_mega.m
- For grids where even GFADI_LOWRANK's m x s factor does not fit: only the outputs selected by OUTFUN are accumulated (the sampled output of the paper, Gamma X). Memory is one or two grid arrays plus no*s. Time is p products plus the output updates.
- The basis rule is GFADI_LOWRANK's 'solution' mode.
- The output updates go by chunks of 64 columns, so the whole factor is never copied. At 10^6 outputs x 10^3 columns, copies would cost more than the products.
- `'full'`, W (n x w): also accumulates XW = X*W on the whole grid, for example two consecutive time steps, to check the true residual.
- `info.res` is an estimate of the residual every 50 steps, scaled by nF2 = ||U||^2 ||V||^2 / rho (exact for rho = 1).

### shifts/adv_shifts.m
- Setting. A = H + S, H symmetric with spectrum in E = [E1 E2], S skew. W(A) lies in Omega = { x + i y : E1 <= x <= E2, |y| <= min(c sqrt(x), ymax) }. For A = -div(a grad) + w.grad with a >= 1 and central differences (Dirichlet): c = |w|_2 and ymax = sum|w_d|/h, because |x'Sx| <= sum_d |w_d| |D_d x| |x| <= |w|_2 sqrt(x'Hx) |x|. The central difference is the mean of two one-sided ones, so |D_d x|^2 is at most the d-part of x'Hx.
- Bound. Crouzeix: ||q(A)|| <= 2 max_Omega |q| (the constant 2; Crouzeix and Palencia proved 1 + sqrt 2). So the one-sided residual after p steps is at most sigma = 2 max_Omega |q| / min_F |q|, with q(z) = prod(z - alpha_j) and F = -Re W(L) = [-WL(2), -WL(1)]. This is the one-sided bound of Example 2 in the paper with the interval E replaced by Omega and Crouzeix's factor added. For c = 0 the factor is 1.
- Search. The zeros are Chebyshev points of [a', b']. a' and b' are searched on a small grid, and for each pair p is found by growth then bisection, keeping the smallest p with sigma <= tol. sigma is evaluated on dense samples of the upper boundary of Omega (q has real zeros, so |q| is symmetric about the real axis).
- Why a' > c^2/4. Chebyshev points of [a', E2] can only converge on Omega if c < 2 sqrt(a'): far along the top curve the rate is (2 sqrt(a') - c)/sqrt(E2). So only a' > c^2/4 is tried. Above E1, the interval leaves part of Omega outside.
- Right end b' = E2 + g y(E2)^2/(4 a'). Where Omega is ymax tall at E2, Chebyshev points that end at E2 grow like exp(p sqrt(2y/L)) there. Moving the end past E2 by y^2/(4a') brings that growth below the decay at F.
- Samples: the zeros are real and at least the height y of the top curve away from it, so |q| varies there on the scale y, and the samples are STEP*y apart. |q| is a multiple of T_p((2z - a1 - b1)/(b1 - a1)) and is evaluated in closed form in logs: log|T_p(t)| = p log|w| + log|1 + w^(-2p)| - log 2, w = t + sqrt(t-1) sqrt(t+1), |w| >= 1. This is stable for large p.
- `info.limit = 2*sqrt(E1)`: the largest c for which the Chebyshev points of [E1, E2] converge at all.
- `info.psym`: the number of steps for E alone, as if A were symmetric.

### shifts/leja_order.m
- Leja ordering keeps the intermediate ADI residuals from growing (cf. Penzl's ordering of ADI shifts).

### shifts/lmin_cw.m
- Proof of the bound. Let s = max(D) and K = sI - H >= 0. If H z >= c z for some z > 0, then K z <= (s - c) z, so rho(K) <= s - c. For symmetric H, lambda_min(H) = s - lambda_max(K) = s - rho(K) >= c (Perron-Frobenius). The bound holds for every positive z; a good z only makes it sharp.
- Choice of z: z1 ~ H^(-1) 1 (the torsion function) and z2 ~ H^(-1) z1 (one step of inverse iteration), each by plain CG stopped by a pointwise residual test.
- Rounding: (H z)_i is computed in double precision, so each ratio is lowered by the margin 64 eps (|H| z)_i / z_i, with |H| z = 2 D z - H z for a Z-matrix.
- Used by SWITCH_PROBLEM (`Ebottom = 'cw'`), where the bound costs products with H only.

### shifts/opt_shifts_real.m
- For real intervals the optimal zeros lie in E and the optimal poles lie in F (Istace and Thiran). The objective is log max_E|r| - log min_F|r|, smoothed by log-sum-exp at temperatures `temps`, which are increased stage by stage.
- Parametrization: alpha_i = mid(E) + half(E)*cos(theta_i) and beta_j = mid(F) + half(F)*cos(phi_j). This keeps the problem well scaled even when the shifts cluster at the ends.
- `nknots`: for large P the angles theta_i = theta(u_i), u_i = (i-1/2)/P, vary smoothly with the index. Then only a piecewise-linear correction to the Chebyshev angles pi*u_i is optimized, with NKNOTS knots clustered at both ends. This reduces the unknowns from P to NKNOTS and the cost by an order of magnitude, with no visible loss for P >~ 100.
- It complements ZOLO3, which is accurate for moderate degree drops but loses the type constraint when p and q are far apart (e.g. p = 48, q = 2).
- Sampling pitfall (geomsamples). Intervals that span many decades (e.g. [10, 6e7]) need geometric samples, 40 per decade, at both ends. Chebyshev points alone leave the region near the small end essentially unsampled, and a coarser geometric grid lets the optimizer hide narrow peaks of r between samples. Observed: a claimed ratio of 6e-9 that is really 3e-5 on a fine grid. That is why the returned SIGMA is always recomputed by RATIO_REAL.
- Initial guess (geoinit). When an interval is longer than 50 times the gap to the other one, the optimal shifts are spread geometrically in the distance to the near end. Starting from linear Chebyshev points then ends in a poor local minimum.
- The Chebyshev samples are offset by 0.3711 so that no sample coincides with an initial Chebyshev zero.

### shifts/ratio_real.m
- Between consecutive zeros, |r| has a single hump of moderate width relative to the piece. So 32 Chebyshev points per piece, plus 40 per decade clustered towards both ends of wide pieces, resolve the maximum to plotting accuracy.

### shifts/shift_search_real.m
- For q > 0 the Zolotarev ratio decays geometrically in p, so the search needs only a handful of OPT_SHIFTS_REAL calls. Warm starts resample the Chebyshev angles of the zeros of the nearest p already solved (they vary smoothly with the index) and are about an order of magnitude cheaper than a cold start.
- Each p keeps its own optimized sigma; sigma(p) is not made monotone.

### shifts/pair_sched.m
- Why pair. A pole close to the end of E amplifies the low end of the spectrum of the residual factor by about |E|/dist(beta,E) relative to the rest. Applied on its own (a plain left step), this dynamic range is not undone until many right steps later. By then rounding has destroyed some components, and those are amplified. With the spectra of A and B nearly touching (dist/|E| ~ 1e-5), this destroyed the accuracy of every unpaired schedule for q >= 6. Paired, every step has a factor of modest size on E: (A - alpha)(A - beta)^{-1} is a Cayley-type transform. For real intervals the partner is essentially the mirror image of the pole.
- For p < q, PAIR_SCHED raises an error. ORDER_SHIFTS handles that case by calling PAIR_SCHED(beta, alpha, F, place) and swapping R and L.

### zolo/zolo3.m
- Degree drop. 'delta' D constrains the Z3 rational to type (N-D, N): D of its zeros are pushed to infinity. In ADI this means D fewer solves with that matrix per sweep. D > 0 costs accuracy. Under exact circularity the attainable value degrades from exp(-N/cap(E,F)) to roughly exp(-(N - D*w(inf))/cap(E,F)), where w is the harmonic measure of dF in the complement of E u F.
- Method. Write the Z4 approximant as rhat = (S - rho*T)/(S + rho*T), where S(z) = sum_k S_k/(z - t_k) and T(z) = sum_k T_k/(z - t_k) over m = N+1 support points t_k. Then rhat - sign = -2S/(S + rho*T) on E and +2 rho T/(S + rho*T) on F. The Z3 rational is exactly r(z) = (1 + rhat)/(1 - rhat) = S(z)/(rho*T(z)), a quotient of two degree-(m-1) polynomials over the common denominator prod_k (z - t_k), which cancels. So S must be small on E and T small on F, and these are two separate linear Chebyshev problems. This is forced: for sign data the Loewner matrix is exactly block anti-diagonal, so the linearized least-squares problem never sees the relative scaling rho between the two blocks. rho is instead fixed in closed form by balancing the two error maxima, and sigma does not depend on it.
- Two points matter for getting the true minimax:
  (i) The Lawson phase must not use the interpolatory AAA parametrization (alpha_k = w_k f_k). That pins the error to zero at all m support points, which caps the number of alternations at about m+1 instead of the 2N+2 required by the Chebyshev characterization. S and T are therefore independent vectors over all support points.
  (ii) The rows of each block are scaled by 1/|T| (on E) and 1/|S| (on F), so each subproblem is a genuine linear minimax problem for |S/T| and |T/S| rather than for |S| and |T|.
- Type constraint. It is linear in S: deg(numerator) <= N-D iff sum_k S_k t_k^j = 0 for j = 0..D-1. It is imposed at every iteration by projecting onto a null-space basis (TYPEBASIS, on shifted and scaled t, since the raw Vandermonde matrix is hopeless), never by post-processing a converged unconstrained solution. (In the interpolatory AAA parametrization S_k vanishes at the support points in E, and the condition collapses to a sum over the support points in F only.)
- Phase 1. The greedy AAA phase uses the Trefethen-Wilber singular-vector blend ('sign' option, BLENDVEC). Without it, the block anti-diagonal Loewner matrix makes the minimal singular vector supported on one block only, and the greedy phase makes no progress at all (the error stays at 2).
- Accuracy. On two disjoint disks, where the exact answer is known in closed form, it reproduces sigma to full plotting accuracy. On two real intervals the Z4 error attains exactly 2N+2-D alternations, which certifies optimality by the Chebyshev characterization for rational approximation. In double precision the practical floor is sigma ~ 1e-28 (tau ~ 1e-14). Beyond that the error curve is at roundoff level and the alternation structure degrades.
- `out.tau` = 2*sqrt(sigma)/(1 + sigma); 2*sqrt(sigma) is its small-sigma limit.
- Poles and zeros are found by a generalized eigenvalue problem on the barycentric numerators. A degree drop of d shows up as d extra infinite eigenvalues, which are discarded. This is a direct check that the type constraint took effect.

### zolo/zolo_gallery.m
- Output files (when 'save' is true) are written to the current folder, not to app_root()/figs. File stems: zolo_gallery_<1|2>, zolo_gallery_<name>, or zolo_row / 'name'.
- The color scale u = 1 - 2 log|r| / log(sigma) is the same for all panels, so single-panel figures stay comparable with the gallery. It runs in MATLAB and in GNU Octave (PDF export falls back to print -dpdf).
- g_combs: the samples are mirrored so that they are exactly symmetric about the imaginary axis, which makes the computed r symmetric too. Setting the environment variable ZOLO_COMBS_SYM=0 turns this off.

---

## 2. Competitors and problem builders

### solvers/lanczos2_sylv.m
One-sided block Lanczos Galerkin for A X + X L^T = U V^T, two passes (cf. Frommer and Simoncini). A SPD (large, products only), L small (kept exact). This is "Lanczos" in Fig. 6.
- Pass 1 runs block Lanczos, X ~ Q_k Y with T_k Y + Y L^T = E_1 R_0 V^T (T_k block tridiagonal), stores only the small blocks, and stops when ||B_k Y_k||_F <= tol*||U V^T||_F (the exact Galerkin residual), checked on a geometric schedule (first at k = 40, then x1.25) and then bisected down to 2% (at least 2 iterations).
- Pass 2 regenerates the basis with exactly the same operations (bitwise the same basis) and accumulates X = Z*Vy' (Y = Uy*Sy*Vy' truncated at 1e-12), so the memory is m*(rank + 3*rho).
- The small equation is solved in the complex Schur basis of L^T by column substitution.

### solvers/sylv_schur_direct.m
Schur form of the small matrix, L^T = Q T Q^*, then n shifted solves with A:
(A + T_jj I) xt_j = ft_j - sum_{i<j} T_ij xt_i, X = Xt Q^*.
The real Schur form is used when it is triangular. Otherwise (complex eigenvalue pairs, which rounding produces for strongly nonnormal L even when the exact eigenvalues are real) the complex Schur form is used, so that T is triangular and the substitution is exact. Using the real form and ignoring its 2-by-2 blocks is wrong. This is the reference solution of Fig. 7(b).

### solvers/sylvester_galerkin.m
- range(QA) = K(A, U; polesA), range(QB) = K(B', V; polesB), one block per iteration. An entry Inf means a plain product, a finite entry s a solve with (A - sI) resp. (B' - conj(s)I).
- Y solves the projected equation (QA'*A*QA) Y - Y (QB'*B*QB) = (QA'*U)(V'*QB).
- All poles Inf: polynomial Krylov (PKSM). Inf poles on the A side and the g-fADI zeros on the B side: spans the same spaces as one-sided g-fADI (Druskin, Knizhnerman and Simoncini). Zolotarev or adaptive poles on both sides: RKSM.
- If a side's Krylov space is exhausted (e.g. B is small) that side stops growing.
- The residual R = A*X - X*B - U*V' is evaluated exactly (not estimated) from thin QR factors of [AQ QA U] and [QB*Y', -BQ*Y', -V], so the comparison with g-fADI uses the same stopping criterion.
- The next block is built from the last orthonormal block (continuation block), which keeps the recurrence well scaled. Orthogonalization: classical Gram-Schmidt twice.

### problems/build_age3d.m
-div(a(x) grad u) + v du/dt - eps d2u/dt2 + mu u = f(x) g(t), (x,t) in (0,1)^3 x (0,1), Dirichlet conditions; t is a structure variable (age, size, ...) with drift v, dispersion eps and loss rate mu. This is Example 2.
- Space: 7-point finite differences, h = 1/(N+1), N^3 interior nodes, face coefficients from a smooth heterogeneous field a = contrast^((1 + sin(2 pi x + 1) sin(3 pi y) sin(2 pi z + 0.5))/2) in [1, contrast] (or opts.afun, with a >= opts.amin > 0).
- Structure variable: centered differences, ht = 1/(n+1), n unknowns.
- With X(:,j) = u(., t_j): A X + X L^T = U V^T, A (N^3 x N^3, SPD, sparse), L (n x n, nonnormal), L = eps*T2 + v*D1 + mu*I.
- Enclosures: lminA = amin * 3(2 - 2cos(pi h))/h^2 (A >= min(a) times the Laplacian), lmaxA = Gershgorin; Re W(L) in [reLmin, reLmax] because the advection part is skew-symmetric.
- opts.nodal: take a at the nodes and give every face the smaller of its two nodes' values (the boundary faces their node's value), so a = 0 at a node cuts it off completely; then amin may be 0 (impermeable walls; the wall nodes decouple). P.anode holds the nodal values.
- Right-hand side: two Gaussian sources in x (centers [0.3 0.4 0.5], [0.7 0.6 0.4], widths 0.1, 0.15) times two Gaussian pulses in t (centers 0.1, 0.5, widths 0.05, 0.2).

### problems/build_frac3d.m
D_t^alpha u - div(a(x) grad u) = f(x) g(t), x in (0,1)^3, t in (0,Tend], u(.,0) = u0, Dirichlet conditions in x. Used only by the 'frac' case of EXP18_ANALYSIS (not drawn in the paper).
- Space: the 7-point operator of build_age3d (coefficient a in [1, contrast]).
- Time: the L1 scheme on the graded mesh t_j = Tend (j/n)^r (l1_graded), r = (2-alpha)/alpha by default, which resolves the t^alpha singularity caused by the initial data.
- With w = u - u0 and X(:,k) = w(., t_k): A X + X T^T = U V^T, U = [f, -A u0], V = [g(t_k), 1]. T (n x n) is dense lower triangular and strongly nonnormal; its diagonal T_kk = tau_k^-alpha/Gamma(2-alpha) varies over orders of magnitude.
- Enclosures: lminA, lmaxA rigorous as in build_age3d; Re W(T) = [reTmin, reTmax] from the eigenvalues of the symmetric part.
- Source: amplitude 50*(0.1/fw)^3 keeps its mass fixed; switched on around t = Tend/2 (Gaussian of width 0.1 Tend). With opts.nodal, u0 and f are zeroed inside sealed walls.

### problems/l1_graded.m
Mesh t_j = Tend (j/n)^r, j = 0..n, and the lower triangular T with
(D^alpha w)(t_k) ~ sum_{j=1}^k a_{k,j} (w_j - w_{j-1}) = (T w)_k (w_0 = 0),
a_{k,j} = [(t_k - t_{j-1})^(1-alpha) - (t_k - t_j)^(1-alpha)] / (Gamma(2-alpha) tau_j),
so T_kk = a_{k,k} = tau_k^(-alpha)/Gamma(2-alpha) and T_kj = a_{k,j} - a_{k,j+1} (j < k).
r = 1 is the uniform mesh; r = (2-alpha)/alpha restores the optimal order 2-alpha for solutions with the typical t^alpha singularity.

Pitfall: the differences of powers are evaluated as y^b*expm1(b*log1p(tau/y)). On strongly graded meshes tau_1 ~ n^-r is far below the rounding level of t_k, and the naive formula loses all digits.

### problems/stencil3.m
y(i,j,k) = d(i,j,k) x(i,j,k) - sum over the six neighbors of a_face x_neighbor, applied matrix-free over z-slabs by a thread pool: 3 to 5 times faster than MATLAB's sparse product or array slicing in our tests, and it needs no stored matrix, so 10^8 to 10^9 cells fit in memory.

SPEC forms:
- constant: .cx .cy .cz (scalars: coefficient/h^2 per direction), .dirichlet (default true; the no-flux variant is not implemented and raises an error), .d0 (scalar added to the diagonal).
- variable: .Ax (nx-1 x ny x nz), .Ay (nx x ny-1 x nz), .Az (nx x ny x nz-1) face coefficients already divided by h^2, .d (nx x ny x nz) the diagonal with boundary faces and reaction included.
- cells: .c (uint8, nx x ny x nz, values 0..254), .tab (up to 256 values: coefficient/h^2 of each level), .f (reaction, default 0). Faces carry the mean of their two cells; Dirichlet (half a cell away) on the box. One byte per cell instead of four arrays, about twice the time per product. lmax = 12*max(tab) + f.

All arrays keep the class of the coefficients (single or double) and x must match. A thread pool (parpool('Threads')) is started if none is running.

### problems/switch_problem.m
Example 3. A searcher diffuses in the unit cube while its internal state follows a driven, irreversible Markov chain, and its diffusivity depends on the state. It is absorbed by a spherical target and degraded at rate delta. u_k(x) = probability of reaching the target before degradation, starting at x in state k:
    D_k Lap u_k + sum_j Q_kj (u_j - u_k) - delta u_k = 0, u = 1 on the target, no flux on the walls.

- Cell-centered finite differences on N^3 cells (h = 1/N), target cells removed with Dirichlet faces (the half-cell offset is ignored: first order at the target). Target cells keep a decoupled diagonal 6/h^2 (in the middle of the spectrum) and X = 0 there.
- States in the columns: A U Dd + U G.' = b d.', A = -Lap_h (SPD), Dd = diag(D_k), G = delta I - Q (a nonsingular M-matrix), b = the Dirichlet data (rank one).
- With G' = Dd^(-1/2) G Dd^(-1/2), the diagonal scaling of the paper's appendix on the switching diffusion (S = sqrt(y./x), x = G'^-1 1, y = G'^-T 1, f = min(1./x + 1./y)/2) and U = Y Dd^(-1/2) S^(-1): A Y + Y L.' = b (S Dd^(-1/2) d).', L = S G' S^(-1), Re W(L) >= f > 0. P.back = S.^-1./sqrt(Dk) maps back: U(:,k) = Y(:,k)*back(k).
- Internal state: a driven cycle of nc phases (forward rate kcyc with one rate-limiting step kslow at phase round(0.6 nc), backward rates 1/drive of them) times three modes of motion, free, sliding and bound (state k = i + (r-1) nc). Sliding starts (kon) only in the first half of the cycle, binding from sliding (kbind) only in the second; koff, kunbind everywhere.
- E = [Elo, Gershgorin]. Elo = 0.98*lmin by eigs at N <= 32, else 0.9*lmin of the N = 32 problem. opts.Ebottom = 'cw': instead a rigorous lower bound from products with H only (lmin_cw, Collatz-Wielandt with z ~ H^-1 1 and H^-2 1 by CG); its time is P.tE and its products P.Einfo.products, to be counted by the caller.
- Shifts: shift_search_real on E and [-1.01 max Re W(L), -f], Leja ordered.
- opts.swirl = W > 0: a force field F with no potential (a vortex about the vertical axis through the center, with a wobble in z: F = W (-(y - 1/2), x - 1/2, 0.3 sin(2 pi x)) / max|.|, at the faces, divergence free), drift D_k F in state k, so the large side becomes A = H + S = -Lap - F.grad with S skew-symmetric face by face (no faces into the walls or the target). W(A) lies in Re z in [lmin(H), lmax(H)], |Im z| <= min(c sqrt(Re z), sum_d max|F_d|/h), c = sqrt(sum_d max|F_d|^2); the shifts are those of adv_shifts. The product is a fused 7-point product over z-slabs (face coefficients A -/+ F in the two directions). P.tshift includes building these arrays. Fig. 9 uses swirl W = 1.
- mv_swirl (unused) is the unfused version of the same product.

---

## 3. Runs behind the figures

### teaser/teaser_run.m, teaser/drive.sh (Fig. 1)
- The model problem of Section 1: A = 7-point -Laplacian plus a centered drift w = 5 (1, 0.6, 0.3)/|(1, 0.6, 0.3)| on N^3 points, L = I + (1e4 - 1)/4 tridiag(-1, 2, -1), n = 200, rank-2 random right-hand side, target 1e-8.
- g-fADI: one-sided, zeros from ADV_SHIFTS on the field-of-values region, Leja order, compression on the small side. Classical fADI uses closed-form type-(k,k) Zolotarev shifts (CLASSICAL_SHIFTS_CLOSED, no Optimization Toolbox needed).
- Equal accuracy: extended Krylov and block Arnoldi stop at the residual g-fADI reached at the same N, read from the CSV, so g-fADI must run first. drive.sh runs the methods in that order, one MATLAB process per (N, method).
- Peak memory is the operating system's high-water mark above what is held after the problem is built. On Linux: VmHWM after a reset. On Windows: PeakWorkingSet64, which cannot be reset, so the note column records the slack (start-up peak above the post-build working set); a method whose extra memory is below that slack reads as the slack.
- The exact residual is computed after the time and the peak memory are taken.

### experiments/exp16_blowout.m (Fig. 8)
- Problem: Example 2 (BUILD_AGE3D), N^3 x n unknowns (up to 2.1e6 x 400 = 8.4e8). A sparse SPD 7-point, coefficient contrast 10: products are cheap, shifted solves expensive (3-D fill-in). L tridiagonal and nonnormal (drift-dominated in the structure variable). Rank-2 right-hand side.
- g-fADI: one-sided (products with A, tridiagonal solves with L), Chebyshev zeros for the rigorous enclosures, Leja order, compression on the small side (GFADI_LOWRANK): memory O(m s).
- direct: Bartels-Stewart with the Schur form of L (SYLV_SCHUR_DIRECT), n shifted solves with A.
- eig: diagonalize L instead (what a multishift Krylov approach needs): cost as direct, but unstable for nonnormal L. Uses exact shifted solves, to isolate the effect of the eigenvector basis.
- classical: p = q fADI with optimal shifts: k shifted solves with A.
- Galerkin: projection on the same spaces (polynomial in A, rational in L^T): O(m k^2) orthogonalization and O(m k) memory. k is the smallest that reaches the g-fADI residual (doubling then bisection to 5%).
- Beyond the size caps in OPTS the expensive baselines are not run (NaN); one shifted solve is timed up to maxfact. MAKE_PAPER_FIGS plots measured values only and merges the reruns tagged '_gfadi', '_ext', '_ext96', '_ext40'.

### experiments/exp16b_nonnormality.m (Fig. 7(b))
- Smaller eps means more drift-dominated and more nonnormal L (eigenvector condition number ~ exp(v/eps)). For each eps: one-sided g-fADI with zeros from Re W(L), classical fADI with shifts for the real intervals, and SYLV_SCHUR_DIRECT as the reference.

### experiments/exp18_analysis.m (Fig. 7(a))
- Theory: for A X - X B = U V^T with A SPD (spectrum in E = [a, b]) and W(B) in Re z <= -f, a one-sided run with zeros alpha_j in E satisfies X_k - X = -r_k(A) X r_k(B)^{-1}, r_k(z) = prod_{j<=k}(z - alpha_j), and since ||(B - alpha)^{-1}|| <= 1/dist(alpha, W(B)) <= 1/(alpha + f): ||r_k(B)^{-1}|| <= 1/|r_k(-f)| and ||X_k - X|| <= max_E|r_k| / |r_k(-f)| * ||X||.
- Measured along the run: err(k), the true relative error (2-norm and Frobenius; 2-norm only for k <= 60, every 5th step and the last); bnd(k) = max_E|r_k| / |r_k(-f)| (the bound, constant 1); sharp(k) = ||r_k(B)^{-1}|| |r_k(-f)| (<= 1, how sharp the B factor is); eigp(k) = max_E|r_k| / min_{lambda(B)}|r_k| (what normal theory would say).
- Classical fADI with the Zolotarev shifts of E and F = Re W(B): the ratio on the real interval (designed) and on the computed boundary of W(B) (relevant); singular values of X against the rank bound sigma_{rho k + 1}(X) <= (1 + sqrt 2) max_E|r| max_{W(B)}|1/r| ||X||_2.
- Age cases use the exact eigenvalues of the tridiagonal Toeplitz L (the computed ones are inaccurate, cond(W) ~ 1e10 to 1e31).
- Fig. 7(a) uses the first case (age, eps = 0.02, N = 24). The 'frac' case (BUILD_FRAC3D) is computed by default but not drawn.

### experiments/switch_drive.m, experiments/switch_run.m (Fig. 9)
- SWITCH_DRIVE(1) runs everything Fig. 9 needs: g-fADI, Arnoldi and extended Krylov ('eksm') into results/switch_swirl1.csv, and g-fADI with the rigorous lower bound on the spectrum (opts.Ebottom = 'cw', LMIN_CW) into results/switch_swirl1_cw.csv. Fig. 9 draws the 'cw' g-fADI rows.
- Each child is `matlab -batch "cd(<app_root>); setup_paths(); switch_run(...)" -logfile logs/...`, so that each peak memory is the operating system's number and a crash costs one row. Peak memory is read from the .NET process object (Windows only).
- 'gfadi': Chebyshev zeros on E, f from the scaling in SWITCH_PROBLEM, a priori stop. The time counts the bound on E (P.tE) and the shift design (P.tshift). Outputs: the slice through the center of the target (all states) and the mean over the cube, saved as Xo and scaled back to u (columns times P.back).
- 'arnoldi': one-pass block Arnoldi (TB_ADV_ARNOLDI), the competitor for the swirl. 'lanczos': one-pass block Lanczos (TB_MEGA_LANCZOS), for the symmetric problem.
- 'gs': block Gauss-Seidel over the states, (D_k A + G_kk) u_k = D_k b - sum_{j~=k} G_kj u_j, a warm-started CG solve per state and sweep (products only), full field in memory. After opts.tcap seconds (and at least 3 sweeps) the remaining sweeps are projected from the measured contraction per sweep and the last sweep's time (a lower bound when CG iteration counts grow). Not drawn.
- 'bs': Bartels-Stewart with the complex Schur form of L: n complex shifted solves with A. One shifted solve (a typical complex shift, median real part) is timed and its fill measured, then multiplied by n. Not drawn.
- 'check' (small N): true residual with the sparse A, and the g-fADI outputs against a converged block Gauss-Seidel solution with sparse direct inner solves (to 1e-13).
- 'eksm' (EKSM_OUTPUTS): extended Krylov (Simoncini 2007) for the scaled equation A Y + Y L.' = U V' with the sparse A: one sparse LU of A (P A Q = L U), then span{U, A^{-1}U, A U, A^{-2}U, ...}, two vectors at a time with full orthogonalization (twice, A not symmetric), Galerkin with the small side kept exact: T Yk + Yk L.' = (V_k' U) V'. T = V' A V is updated as the space grows (products with A for the new columns, with A' for the new rows), so the basis V is the only m x k array. Residual from A V_k = V_{k+1} Tbar (exact for this space).

### experiments/freq_sweep_run.m (Fig. 10(d,e))
- Example 4: (A + i omega_j I) x_j = b, j = 1..1000. a = 1..10 smooth, or (C > 0) 27 spheres of radius 0.11 on a jittered 3x3x3 lattice with conductivity C in a = 1; Fig. 10 uses C = 1e3. A is nonsymmetric, so no Lanczos. lambda_min(H) >= 3 (4/h^2) sin^2(pi h/2) because a >= 1: a free rigorous bound (no eigs, no shifted solve).
- 'gfadi': real zeros on Omega (ADV_SHIFTS) with margin tol/4 (the bound is on the residual of the whole block); X = Z*Y' compressed on the small side to ytol = min(1e-13, 0.1 tol E1/E2), which keeps the residual below tol.
- 'msfom': multishift FOM, one Arnoldi basis for all shifts (stored), X = V*Yk (the basis is its low-rank factor). O(k) Givens update for the hardest (smallest) frequency: FOM residual = GMRES residual / |c_k|; when that is <= tol (checked every 10 steps), all frequencies are checked exactly via a complex Schur form of H_k. Basis capped at 24 GB (on a 32 GB machine).
- 'rsfom': restarted multishift FOM(mr) (Simoncini 2003): residuals stay collinear, r_j = beta_j v; keeps every x_j (m x nw complex), or only the sensor rows for 'rsfom_out'. Fig. 10 uses mr = 200.
- 'bicg_out': multishift BiCG (Jegerlehner 1996; Frommer 2003), outputs only: one BiCG run on A x = b serves every shift through collinear residuals; per step one product with A and one with A^T, no basis, O(ell*nw) work. Derivation in the BICG_OUTPUTS comment.
- 'gmres_out': each frequency on its own, GMRES(100), sensor values only; stops at the time cap and records how many frequencies were done (measured, no extrapolation).
- Output methods (`*_out`, used in Fig. 10): accuracy is the relative sensor error against sparse direct solutions at 12 frequencies.
- peakGB is measured before the residual check.

### experiments/intro_slices_gfadi.m (Fig. 10(b,c))
- The slices are computed with one-sided g-fADI (products with A only), with the operator and shifts of FREQ_SWEEP_RUN (C = 1e3), so that the plane z = 1/2 can be drawn on a fine grid (N = 127). N must be odd so that z = 1/2 is a grid plane. CHECK = true compares with a direct solve on small grids.

---

## 4. The test bed (Example 1)

### testbed/tb_problem.m
- `coef = 'het'`: -div(a grad), a = 10^(contrast*field), field a smooth random function in [0,1] (a sum of six Gaussians, normalized on a 25^3 sample grid); a is evaluated at the faces.
- `coef = 'adv'`: -Laplacian + w (d/dx + d/dy [+ d/dz]), central differences (nonnormal large side). The cell Peclet number is w*h/2. sym(A) is still the Laplacian.
- `Ltype`:
  - `'diag'`  shifted systems, diag of n log-spaced values in [f, g]
  - `'tri'`   1-D diffusion: f I + (g - f)/4 tridiag(-1, 2, -1) (spectrum in [f, g]); the base case of Example 1
  - `'time'`  backward Euler in time, n steps with tau = 2/(g - f), plus f I (lower bidiagonal: one eigenvalue, a Jordan block, Re W(L) in [~f, g])
  - `'adv'`   1-D advection-diffusion, Peclet 'pe', shifted by f (nonnormality grows with pe); Fig. 7(c)
  - `'dense'` dense normal: Q diag(log-spaced in [f, g]) Q', Q random orthogonal
- `rhs`: 'rand' (iid normal), 'smooth' (sum of a few wide Gaussians), 'point' (spikes).
- `P.E` are the two spectral numbers g-fADI uses: analytic for 'lap' and 'adv', and min(a)*lambda_min(lap) (bottom) and Gershgorin (top) for 'het'.
- `P.WL` is [min, max] of Re W(L), the eigenvalues of the Hermitian part of L.

### testbed/tb_method.m
- `'gfadi'`: one-sided g-fADI, right steps only (products with A, shifted solves with the small L). Chebyshev zeros from E = [lambda_min, lambda_max](A) and Re W(L) (SHIFT_SEARCH_REAL with q = 0), Leja-ordered, compression on the small side (GFADI_LOWRANK). Stops a priori. Refuses p > 2e5.
- `'mixed'`: g-fADI of type (p,q), q > 0 left steps (shifted solves with A - beta I). q is chosen from opts.qlist by the cost model p(q)*C_R + q*C_L, where C_R (one right step: product with A and a small solve, mean of 5) and C_L (one shifted solve with A - F(2) I) are measured. The type-(p,q) design takes minutes for large p, so it is skipped when p0 > 3000 and the q loop stops after 90 s. Pairs first (ORDER_SHIFTS).
- `'fadi'`: classical factored ADI, type (k,k), optimal shifts, k shifted solves with A. One solve is timed; if k times that exceeds opts.tmax, the time is projected as k*t1 and nothing more is run.
- `'bs'`: Bartels-Stewart with a sparse large side (the modified Bartels-Stewart of Example 1): real Schur form of L (complex if it has 2x2 blocks), then n shifted solves with A. After 3 solves, if the projected total toc*n/j exceeds opts.tmax, it stops and reports the projection. Refuses if X (m x n) would exceed opts.memcap GB. When complete, its X is the reference for every other method's error (TB_SWEEP).
- `'eksm'`: extended Krylov (Simoncini 2007), one factorization of A (CHOLMOD if SPD, else sparse LU via `decomposition`), Galerkin, residual checked every 2 steps; capped at memcap GB or 4*tmax.
- `'lanczos2'`: two-pass block Lanczos Galerkin (LANCZOS2_SYLV), matrix-free, O(m rank) memory; SPD A only.
- `'arnoldi'`: one-pass block Arnoldi/Lanczos Galerkin with full orthogonalization (stores the basis; works for nonnormal A). Residual on a geometric schedule (factor 1.15); capped at 2*tmax.
- Krylov methods stop at the exact residual opts.restarget (normally the residual g-fADI reached), so all methods are compared at equal accuracy.
- memGB is an estimate of the peak working set: stored columns of length m times 8 bytes, plus the factorization fill (opts.fillGB, measured by 'eksm') for methods that factorize.

### testbed/tb_sweep.m
- Fig. 6 uses the sweeps `size2d` (a), `size3d` (b) and `nsmall3d` (c, N = 48); Fig. 7(c) uses `Lnonnormal` (2-D, N = 256). The other sweeps are not drawn in the paper.
- The mixed type (q > 0 left steps) needs a type-(p,q) shift design that takes minutes for large p, so it has its own sweep (`mixed2d`) and is left out of the others.
- The `*3d` sweeps ask the same questions on the 3-D base (N = 48), where shifted solves are expensive.
- g-fADI's achieved residual becomes the Krylov target only if it is within 10x of tol (else the a priori target is used), and it is floored at 1e-3*tol.

### testbed/tb_mega_lanczos.m
One-pass block Lanczos Galerkin with restricted outputs; used by SWITCH_RUN ('lanczos', symmetric switching problem).
- Keeps the three-term blocks and C'*Q (the outputs of the basis), never Q itself, so there is no second pass and no basis in memory.
- The projected equation T_k Y + Y L' = E_1 R_0 V' is solved by column substitution in the triangular L'. No Schur form when L is triangular or diagonal (a generic dense Schur form of L would penalize large n unfairly); only a non-triangular L gets a complex Schur form. T_k is block tridiagonal, assembled sparse.
- Exact Galerkin residual ||B_k Y_k||_F, checked first at k = 40 and then at k <- ceil(1.25 k); when it meets tol, bisection between the last two checks finds k to within 2% (as in LANCZOS2_SYLV).
- Vector operations are MATLAB's (multithreaded); block QR by CholeskyQR2.
- opts.kcap (expected iterations) preallocates the output store so that a large store is not copied as it grows.

### testbed/tb_adv_arnoldi.m
One-pass block Arnoldi Galerkin with restricted outputs for a nonsymmetric A; the competitor in Fig. 9 (SWITCH_RUN 'arnoldi').
- The basis is kept in chunks of 16 blocks (cell array) so that it is never copied when it grows or is used (a slice Q(:, 1:K) would copy it every step). Block classical Gram-Schmidt, twice, chunk by chunk.
- The projected equation H_k Y + Y L' = E_1 R_0 V' is solved through a complex Schur form of the small Hessenberg H_k and column substitution in the triangular L' (a non-triangular L gets a complex Schur form).
- Exact Galerkin residual ||H_{k+1,k} E_k' Y||_F, first at k = 20, then k <- ceil(1.25 k); bisection to 2% after.

---

## 5. Figure scripts

All figure scripts read from and write to folders under `app_root()`. Every plotted point is measured; nothing is extrapolated.

### figures/make_zolo_fig.m (Fig. 2)
- One thread (`maxNumCompThreads(1)`) on purpose: the Lawson iteration then lands on the same (symmetric) solution every time.
- ZOLO_GALLERY saves into the current folder, hence the `cd(out)` with an onCleanup back.

### figures/make_shift_ordering_fig.m (Fig. 3)
- The type-(q,q) Zolotarev shifts for E = [-1, -d] and F = [d, 1] are in closed form (Wachspress), so the figure needs no optimization. The data are saved next to the figure in figs/shift_ordering.mat.

### figures/make_style_mockups.m (Figs. 4 and 5)
- Fig. 4 is `mock_intervals`, Fig. 5 is `mock_ratio_asym`. The default run draws both. Run any local function alone with the environment variable MOCK_ONLY (e.g. `setenv('MOCK_ONLY', 'mock_ratio_asym')`).
- Fig. 5(a), the ratio F(1)/F(1/2) for E = [-1, 1], F = [1+d, 1+d+L]: x = l(E)/l(F) = 2/L, y = log10(d/L), with F(1) = log(1 + d + sqrt(d^2 + 2d)) and F(1/2) = pi K(k')/(2K(k)), k^2 = 2L/((2 + d)(d + L)). Band edges are stroked in their own color so the vector PDF shows no seams. The colorbar label is set in LaTeX in the manuscript, since MATLAB cannot render \mathscr.
- `zolotype` returns the type-(nz, np) rational scaled so that max_E|r| = sqrt(sigma) = 1/min_F|r|. For nz > np it computes the type-(np, nz) rational with E and F swapped and inverts it.

### figures/make_paper_figs.m (Figs. 7 and 8)
- Two modes: 'bound' (Fig. 7) and 'cost' (Fig. 8). The calls that produce the 'cost' data are listed in the header of EXP16_BLOWOUT.
- 'cost' data merging: exp16_blowout.mat is the base; exp16_blowout_gfadi.csv (tag '_gfadi') supplies g-fADI time and memory and the direct time where measured; exp16_blowout_ext.mat and _ext96.mat (runs that extend the competitors to larger N) replace NaNs of the original runs with their measured values (fields t_direct, t_classical, t_gal, k_gal, mem_gal_GB, fill_GB, t_fact); exp16_blowout_ext40.mat (N = 16 and 40, every method in one session) adds new rows.
- 'cost' memory: classical fADI holds its factor (fill_GB), Galerkin its basis (2 rho k_gal vectors of length m, rho = 2), g-fADI as measured. The direct method is computed but not drawn.
- 'bound' (a): the bound of the partial products of the Leja-ordered zeros oscillates from step to step; its running maximum over the later steps, max_{j >= k} bound_j, is still a bound at every step k and is monotone. That envelope is what is drawn.
- 'bound' (b): a narrow window about the target, which shows whether each method meets the target it was designed for.
- 'bound' (c): the test bed's Lnonnormal sweep (2-D, N = 256, n = 200): L an advection-diffusion matrix (centered diffusion, upwind advection) against its Peclet number.

### figures/make_sec4_figs.m (Fig. 9)
- Only 'switch' is in the paper; 'galaxy', 'sponge' and 'big' read results that are not produced in this repository.
- 'switch' draws measured points only. Gauss-Seidel and Bartels-Stewart were never run to completion (their rows are projections), so they are not plotted. Extended Krylov is plotted where it finished.
- 'switch' replaces g-fADI rows with those from switch<tag>_cw.csv when present: g-fADI with the rigorous bound on lambda_min (LMIN_CW), all of its cost counted.
- Variables s, b, xg, ytop and the local function place_lab are unused.

### figures/make_sec5_figs.m (Fig. 6)
- Only 'testbed' is in the paper. 'scaling', 'memory', 'testbed_forest' and 'drift' draw figures of an earlier draft whose runs are not in this repository; with no argument the script tries them all and prints an error for each that cannot be drawn.
- 'testbed' merging: testbed_nsmall_meas.csv (the n-sweep rerun with g-fADI, classical fADI and Bartels-Stewart to completion) replaces those methods' rows of the nsmall3d sweep, so (c) compares times from one session, measured to the end (Bartels-Stewart at n = 3200 was projected and is dropped). testbed_size_meas.csv (size sweeps rerun with no time budget at the sizes where the time budget had cut them off) replaces each row with the same sweep, value and method.
- yl, yt, ytl in 'testbed' and the local function pick are unused.

### figures/make_intro_3d.m and figures/make_sweep_fig.m (Fig. 10)
- MAKE_INTRO_3D renders panel (a): the inclusions, the source and the 100 sensors in the unit cube. The operator A and right-hand side b are still built (and `omega` still accepted) but are not used by the picture; the sensor positions depend only on `rng(3)` and `m = N^3`. The local function `infer_map` is unused.
- MAKE_SWEEP_FIG composes the figure: the render from figs/paper_intro_3d.png (its "(a)" is cropped off), the slices of INTRO_SLICES_GFADI for (b) and (c), and FREQ_SWEEP_RUN's rows for (d) and (e): methods gfadi_out and msfom_out (100 sensors) and rsfom_out with note containing "FOM(200)". A multishift FOM run that failed (note starting with "ERROR") is excluded. The PDF is composed with pdflatex from a vector layer and full-resolution images.

### teaser/make_teaser_fig.m (Fig. 1)
- Reads teaser/teaser.csv, plots measured rows only (projected == 0). The tick at N = 40 is omitted because 40^3 and 48^3 collide at the printed width.
