# g-fADI: one-sided and lopsided factored ADI for Sylvester equations

MATLAB code for the paper on generalized factored ADI (g-fADI) for

    A X - X B = U V*,

where A is large and sparse, B is small (or diagonal), and U V* has low rank. g-fADI takes `p` right steps (a product with A and a solve with B) and `q` left steps (a shifted solve with A and a product with B), so a type-(p,q) iteration uses the p zeros and q poles of a rational function. Classical fADI has `p = q`. When shifted solves with A are expensive, a lopsided iteration (`q` much smaller than `p`) or a one-sided one (`q = 0`, products with A only) is much cheaper. The shifts come from Zolotarev-type rational functions on two sets E (containing the spectrum of A) and F (containing that of B).

Requires MATLAB R2023b or later (tested on R2025b), with the Parallel Computing Toolbox (thread pools for the large 3-D products) and the Optimization Toolbox (some shift optimizations). `make_sweep_fig` also needs `pdflatex` on the system path. Peak memory in the large runs is read through .NET, so those runs report memory on Windows only.

## Quick start

```matlab
run('path/to/repository/setup_paths')   % once per session: adds the folders, creates results/, figs/, logs/
example1_sylvester_basics                % start here
```

| Example script | What it shows | Time |
|---|---|---|
| `examples/example1_sylvester_basics.m` | A small Sylvester equation solved with one-sided, lopsided and classical fADI, checked against `sylvester` | seconds |
| `examples/example2_shifted_systems.m` | Shifted systems `(A + s_j I) x_j = b` with products by A only, against one shifted solve per shift | ~10 s |
| `examples/example3_lopsided_pde.m` | A 3-D PDE with a small nonnormal B, with compression on the small side | < 1 min |
| `examples/example4_zolotarev.m` | Zolotarev rational functions of type `(p,q)` with `zolo3` | < 1 min |

These scripts are small demonstrations, not the paper's Examples 1-4.

The pattern in every example:

```matlab
E = [lambda_min(A), lambda_max(A)];                    % bounds on the spectrum of A
F = [-max(spec L), -min(spec L)];                      % bounds on the spectrum (or Re W) of B = -L'
[p, alpha, beta] = shift_search_real(E, F, q, 1e-8);   % shifts for a target of 1e-8 with q left steps
alpha = leja_order(alpha, E');                         % a stable order
[Z, Y] = gfadi(A, B, U, V, alpha, beta);               % X ~ Z*Y'
% large problems: [Z, Y, info] = gfadi_lowrank(gfadi_op(A), gfadi_op(B), U, V, alpha, beta)
```

## Layout

| Folder | Contents |
|---|---|
| `gfadi/` | The iteration: `gfadi`, `gfadi_lowrank` (compression on the small side), `gfadi_mega` (matrix-free, sampled outputs only), `gfadi_op` (operator structs), `gfadi_simple` (a short reference implementation) |
| `shifts/` | Shift design: `cheb_shifts`, `shift_search_real`, `classical_shifts_real`, `opt_shifts_real`, `ratio_real`, `adv_shifts` (nonsymmetric A), `leja_order`, `order_shifts`, `pair_sched`, `fov_boundary`, `lmin_cw` (a rigorous lower bound on the smallest eigenvalue from products only) |
| `zolo/` | `zolo3` (Zolotarev rational functions by constrained AAA-Lawson), `zolo_gallery` |
| `solvers/` | Competitors: two-pass block Lanczos for Sylvester, Galerkin projection, a Schur-based direct solver |
| `problems/` | Problem builders: `build_age3d` (Example 2), `switch_problem` (Example 3), `stencil3`, `build_frac3d`, `l1_graded` |
| `teaser/` | The model problem of Section 1 (Fig. 1) |
| `experiments/` | The runs behind Figs. 7-10. They write to `results/` |
| `testbed/` | The Example 1 test bed (Figs. 6 and 7c), and two Krylov competitors used by Example 3 |
| `figures/` | The scripts that draw the paper's figures from `results/` into `figs/` |
| `examples/` | Small, self-contained examples |

`results/`, `figs/` and `logs/` are outputs and are not in the repository.

## Reproducing the figures

Run the data scripts first, then the figure script. The headers of the run scripts give the exact calls. Large runs take minutes to hours, and the drivers start one MATLAB process per run so that runs cannot affect each other's timings or memory.

| Figure | Content | Figure script | Data |
|---|---|---|---|
| 1 | Model problem: time and memory | `teaser/make_teaser_fig` | `teaser/drive.sh` (calls `teaser_run`) |
| 2 | Zolotarev gallery | `make_zolo_fig` (uses `zolo_gallery`) | computed on the fly |
| 3 | Shift ordering | `make_shift_ordering_fig` | computed on the fly |
| 4 | Zolotarev rationals on two intervals | `make_style_mockups` (local `mock_intervals`) | computed on the fly |
| 5 | Efficiency ratio, asymmetric intervals | `make_style_mockups` (local `mock_ratio_asym`) | computed on the fly |
| 6 | Example 1: test bed | `make_sec5_figs('testbed')` | `tb_sweep` |
| 7 | Bound and nonnormality | `make_paper_figs('bound')` | `exp18_analysis`, `exp16b_nonnormality()`, `tb_sweep` (`Lnonnormal`) |
| 8 | Example 2: cost | `make_paper_figs('cost')` | `exp16_blowout`, five runs listed in its header |
| 9 | Example 3: switching diffusion | `make_sec4_figs('switch')` | `switch_drive(1)` (calls `switch_run`, including the extended Krylov and rigorous-bound runs) |
| 10 | Example 4: frequency sweep | `make_intro_3d`, then `make_sweep_fig(127, [3e2 3e3])` | `freq_sweep_run`, `intro_slices_gfadi(127)` |

`make_sec4_figs` and `make_sec5_figs` have further modes that draw figures from earlier drafts; their headers say which modes belong to the paper. See `NOTES.md` for file-by-file technical detail.

## Conventions

- **Steps and shifts.** A right step uses a zero `alpha_j`: a product with A and a solve with `B* - conj(alpha_j) I`. A left step uses a pole `beta_j`: a shifted solve with `A - beta_j I` and a product with B. The error after the iteration is `r(A) (X0 - X) r(B)^{-1}` with `r(z) = prod(z - alpha_j) / prod(z - beta_j)`.
- **Sign of B.** Most problems are written `A X + X L' = U V'`, which is `A X - X B = U V'` with `B = -L'`. F then bounds the spectrum (or the real part of the field of values) of `-L`.
- **Timings.** Reported times are wall-clock on one laptop (8-core Intel Core Ultra 7 258V, 32 GB), one MATLAB process at a time, with no plotting during a timed run. Every plotted point is measured.
