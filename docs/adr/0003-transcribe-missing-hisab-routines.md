# 0003 — Transcribe the hisab 1.4.0 routines that hisab 3.4.1 lacks

**Status**: Accepted
**Date**: 2026-11-08

## Context

The oracle (kimiya 1.2.0) calls Rust hisab **1.4.0**. The port depends on Cyrius hisab
**3.4.1**, which covers most of that surface — `num_newton`, `num_bisection`, `num_fft`,
`HComplex`, `csr_from_dense`, `pcg32_*`, `svd_compute_checked`, the `ad_tape_*` autodiff — but
not all of it (port plan §6.2, surveyed 2026-11-08):

| Rust hisab 1.4.0 (kimiya caller) | Size | hisab 3.4.1 |
|---|---|---|
| `num::rk4_trajectory` over a **system** `f(t, &[f64], &mut [f64])` (`reaction_dynamics` first/second-order, consecutive, reversible, Michaelis–Menten) | 34 lines | `num_rk4` is one **scalar** step |
| `num::dopri45` over a system (`simulate_adaptive`) | 180 lines | `ode_dopri45(_trajectory)` are **scalar** |
| `num::euler_maruyama` over a system (`simulate_stochastic`) | 36 lines | `ode_euler_maruyama` is one scalar step |
| `num::gaussian_elimination`, partial pivoting (`simulate_stiff`'s Newton step) | 61 lines | no public dense direct solver (only GMRES / BiCGSTAB / PGS) |
| `num::least_squares_poly`, Householder/MGS QR with a 1e-12 rank threshold (`fitting`) | 43 lines + QR helper | none |
| `num::matrix_rank` = SVD + threshold (`network`) | 13 lines over `svd` | `svd_compute_checked` exists |

The parity policy needs the port to reproduce the oracle's numbers. Replacing these with
different algorithms would turn class-B/D rows into unbounded differences (an adaptive
integrator with a different controller takes different steps).

## Decision

Transcribe the five missing routines from
`~/.cargo/registry/src/*/hisab-1.4.0/src/num/{ode,roots,linalg}.rs` into **`src/numerics.cyr`**
as private `_km_*` functions, in the oracle's operation order, at M5:
`_km_rk4_trajectory`, `_km_dopri45`, `_km_euler_maruyama_system`, `_km_gaussian_elimination`,
`_km_least_squares_poly` (with the QR it uses). Their parity target is **class B**, falling back
to class D only where a transcendental is on the path.

Everything else uses hisab 3.4.1:
- `matrix_rank` → count singular values of `svd_compute_checked` above the oracle's threshold;
  the rank is an integer, compared class A. Transcribe 1.4.0's SVD only if a golden row disagrees.
- `newton_raphson` / `bisection` (adiabatic flame temperature) → `num_newton` / `num_bisection`,
  class D. Their 3.x semantics changed (subnormal derivatives, NaN refusals, exact endpoint
  roots); transcribe the 1.4.0 loops (~15 lines each) only if a class-D row fails.
- `fft` / `Complex` → `num_fft` / `HComplex`, class C for n ≥ 256 (3.4.x moved those bits).
- `CsrMatrix::from_dense` → `csr_from_dense`; `Pcg32` → `pcg32_new` (stream 1, as Rust).
- autodiff → `ad_tape_*`, only for `compute_gradient` (ADR 0002 #2 makes the other gradients
  closed-form).

A hisab feature request for system ODE integrators, a dense direct solver, polynomial least
squares and `matrix_rank` is a follow-up for the maintainer; when hisab ships them, a later
kimiya release can switch and re-measure.

## Consequences

- **Positive** — M5 can reach bit parity with the oracle; the transcribed code is small
  (≈ 350 Rust lines) and already audited (it is what 1.2.0 ran).
- **Negative** — kimiya owns numerical code that belongs in hisab until hisab provides it, and
  must not let it drift from the oracle.
- **Neutral** — the transcriptions are private, so moving to hisab later is not an API change.

## Alternatives considered

- **Write new general solvers in kimiya**: rejected — unbounded parity differences and new,
  unaudited numerical code.
- **Implement the gaps in hisab first**: preferred long-term, but it blocks kimiya's M5 on
  another project's release; kept as the follow-up above.
- **Drive the scalar hisab routines per component**: wrong for coupled systems (A→B→C, A⇌B).

## Amendment (M2b, 2026-11-09)

**The flame solver is transcribed too.** The decision above sent `newton_raphson` / `bisection`
(adiabatic flame temperature) to hisab 3.4.1's `num_newton` / `num_bisection`, class D, transcribing the
1.4.0 loops only if a class-D row failed. Porting `thermochem` showed the hisab route can't be made to
match the oracle's semantics, so M2b transcribes them before any row was compared:

- **Callbacks.** `num_newton` / `num_bisection` take `fncall1(f, x)` callbacks with no context argument.
  The flame balance closes over the products, ΔH and T₀, so the callbacks would need module-level mutable
  state — not re-entrant, not thread-safe (Cyrius threads are real), and against kimiya's no-globals style.
- **Semantics changed after 1.4.0** (read in the vendored `lib/hisab.cyr`, 3.4.1). `num_newton` refuses only a
  derivative below `DBL_MIN` whose step is not finite (1.4.0: any |f′| < `EPSILON_F64` = 10⁻¹²).
  `num_bisection` refuses NaN ends and NaN f values, returns an exact endpoint root, tests the bracket by
  comparing signs (1.4.0: `f(a)·f(b) > 0`), takes the midpoint as `lo + (hi − lo)·0.5` (1.4.0: `(lo + hi)·0.5`,
  different bits), and reports `max_iter` exhaustion as `Err` (1.4.0: `Ok` with the midpoint). On the
  bisection fallback the midpoint formula alone moves the result's bits, and the other changes can turn an
  oracle `Ok` into an `Err`, so even class D could diverge there.

So `src/thermochem.cyr` has two private helpers, `_km_flame_newton` and `_km_flame_bisection`, hisab
1.4.0's `num/roots.rs` lines 28–99 operation for operation, specialised to the flame balance: they call
`_km_flame_f` / `_km_flame_df` directly (no function pointers) with the products as arguments. The
derivative threshold is hisab 3.4.1's `EPSILON_F64` (10⁻¹², the value of hisab 1.4.0's), the tolerances and
iteration limits are Rust's (`FLAME_NEWTON_TOL_J` = 10⁻³ J/mol, 50; `FLAME_BISECTION_TOL` = 10⁻⁶, 200),
inline as typed formula constants (port plan §6.4). Newton's `Err` falls through to bisection, as Rust's
`if let Ok`. The helpers stay in thermochem: M5's `src/numerics.cyr` transcribes the general routines
(system ODEs, the dense solver, least squares) for the other modules.

**Consequence.** The whole flame path is + − × ÷ (Shomate is rational), so the golden rows are **class B,
0 ulp**, not class D: 1,836 `Ok` rows on x86_64 ELF, x86_64 PE and aarch64, 1,682 answered by Newton and 154
by the bisection fallback (`docs/development/port-parity.md` §6). `validate.py` recomputes every row, the
solver path included, and shows that each typed constant (10⁻³, 10⁻⁶, the 40 J/(mol·K) guess) changes some
golden row's bits when perturbed; `tests/thermochem.tcyr` pins 10⁻³ and 10⁻⁶ on the private helpers
exactly, and the 10⁻¹² derivative threshold to within ±1 % (no golden row reaches it). The iteration limits
(50, 200) are transcribed but not observable: no row converges near them, so no output depends on them. The status stays **Accepted**; port plan §6.2's row is updated.

## Amendment (M4e, 2026-11-09)

**hisab 1.4.0's FFT is transcribed too.** The decision above sent `num::fft` / `Complex` (spectroscopy's
`fft_spectrum`) to hisab 3.4.1's `num_fft` / `HComplex`, class C for n ≥ 256. Porting `spectroscopy` showed that route
can't reach the corrected oracle (ADR 0006), so M4e transcribes the routine before any row was compared:

- **The twiddles.** hisab 1.4.0 forms wn = (cos θ, sin θ) once per stage, θ = (−2·π)/len (−2·π and the ÷ len are exact),
  and advances w by `w *= wn` across the stage. The release oracle calls glibc's `sincos` there (LLVM joins
  `angle.cos()` and `angle.sin()`; its bits are glibc's `sin` and `cos`, which the golden rows prove). hisab 3.4.1's
  `num_fft` takes them from the stdlib's `f64_cos` / `f64_sin`, which are not correctly rounded (the `math` rows' shadow
  lines: 82 cos and 20 sin rows one ulp off the corrected oracle), and **reseeds** w from `f64_cos(angle·k)` /
  `f64_sin(angle·k)` every 64 steps (`_NUM_FFT_RESEED`), a different recurrence: its bits move from n = 256 on (port
  plan §6.2). Neither is the oracle's operation sequence evaluated with correctly rounded primitives.
- **So** `src/spectroscopy.cyr` has a private `_km_fft(data, n)`: hisab 1.4.0's `src/num/fft.rs` lines 12–58 operation for
  operation — the bit-reversal loop, the butterflies with `Complex`'s `mul` (re = w.re·d.re − w.im·d.im, im = w.re·d.im
  + w.im·d.re), `add`, `sub` and `mul_assign` from `src/num/complex.rs`, no FMA — with the twiddle from `src/mathx.cyr`'s
  correctly rounded `_km_cr_sincos`; `fft_spectrum` then takes `Complex::abs` = √(re·re + im·im) of the first n/2 + 1
  bins, as Rust (not hypot). hisab's own `n <= 1` exit is kept; its power-of-two check is `fft_spectrum`'s, made first
  as in Rust.

**Consequence.** `fft_spectrum` is class C at **bound 0** against the corrected oracle, as every other class-C function:
213 golden rows (every length 0 … 16, every power of two to 2¹³, 2^k ± 1; impulses, constants, ±0, alternating and
Nyquist signals, sinusoids, random, huge, subnormal, NaN and ±∞ samples), 27,512 floats, all 0 ulp on x86_64 ELF, x86_64
PE and aarch64 (`docs/development/port-parity.md` §15). glibc's sin and cos are correctly rounded at the 13 FFT angles
−2π/2^k (MPFR decides), so no row needs a correction; `validate.py` recomputes every row with its own transcription and
glibc's (then MPFR's) twiddles, and perturbing the −2.0 by an ulp changes 41 rows. The status stays **Accepted**; port
plan §6.2's row is updated. **Follow-up for hisab** (`docs/development/state.md`): hisab's FFT twiddles are not correctly
rounded (`f64_cos` / `f64_sin`), so its transforms differ in the last bits by platform and from a correctly rounded
reference; a hisab FFT on correctly rounded twiddles, or with 1.4.0's recurrence as an option, would let kimiya call it.

## Amendment (the overflow audit's B6, 2026-11-11)

**The FFT gains an unbounded-exponent slow path (ADR 0002 #36), and the flame's bisection an exact sign test (ADR 0002
#38).** `_km_fft` is still hisab 1.4.0's `num::fft` operation for operation, now also returning whether its products
provably stayed in binary64's normal range (every twiddle component 0 or at least 2⁻¹⁰⁰, every butterfly's d 0 or at
least 2⁻⁹²¹); where that test, or the retained outputs' finiteness, fails for a finite signal, `fft_spectrum` re-runs the
same transcription with an unbounded exponent (`_km_fft_ux`, on `src/mathx.cyr`'s `_km_ux_*`), each magnitude's √ rounded
once — the user's decision "FFT: both", hisab 1.4.0's transform without the range limit. In-range transforms keep their
bits (the fast path is the transcription itself); the golden rows that change are ADR 0002 #36's (41 of the 213 M4e
rows). `_km_flame_bisection`'s opposite-signs check, hisab 1.4.0's `f(a)·f(b) > 0`, is the exact same-sign comparison,
which differs only where that product underflows to +0 and cannot change a flame result (the flame reaches it only after
its bracket proved opposite signs at the same two points); the bracket test itself is ADR 0002 #38. The status stays
**Accepted**.

## Amendment (M5a, 2026-11-12)

**The ODE integrators and the dense solve are transcribed as decided; PCG32 is hisab 3.4.1's.** `src/numerics.cyr` has
`_km_rk4_trajectory` (with its `_km_rk4_step`), `_km_dopri45`, `_km_gaussian_elimination` and
`_km_euler_maruyama_system`: hisab 1.4.0's `src/num/ode.rs` (`rk4_step`, `rk4_trajectory`, `dopri45`, `euler_maruyama`)
and `src/num/roots.rs` (`gaussian_elimination`) operation for operation — the stage sums left to right as Rust parses
them (`y + hd·a21·k1` is `y + (hd·a21)·k1`), the Neumaier-compensated RK4 increment with both of its branches, the
Butcher tableau as Rust's const quotients of integers (each evaluated as the same correctly rounded division), the step
clamp and the `t += h` time accumulation, the scratch vectors zeroed once and reused, never re-zeroed (a callback that
leaves an entry alone sees its previous value, as in Rust), the strictly-larger pivot rule and its 10⁻¹² threshold (hisab
3.4.1's `EPSILON_F64`, the same bits), whole-row swaps, and back substitution from b down the row. The closures are a
function pointer and a context pointer (`f(ctx, t, y, dy, n)`, called through `callptr`), and the trajectory a flat
buffer of (t, y[n]) points (port plan §4, M5). `backward_euler` is not transcribed: `simulate_stiff` never called it
(it is kimiya's own loop, in `src/reaction_dynamics.cyr`).

- **PCG32.** `pcg32_new(seed)` / `pcg32_next` of hisab 3.4.1 (stream 1: inc = 3, the advance, `state += seed`, the
  advance) give 1.4.0's `Pcg32::new(seed, 1).next_u32()` bit for bit: `tests/golden/numerics.tsv`'s 22 `num::Pcg32::next_u32`
  rows (seeds 0, 1, 2, 7, 42, 2024, 2³², 2⁶³, u64::MAX − 1, u64::MAX and ten seeded ones; 64 or 2,000 draws) are all
  equal, so nothing was transcribed there. 1.4.0's `next_f64` (u / (u32::MAX + 1), exact) and `next_normal` (Box–Muller,
  √(−2·ln max(u₁, 10⁻³⁰⁰))·cos(τ·u₂)) are transcribed on its output, as `_km_pcg32_next_f64` and `_km_pcg32_next_normal`.
- **Classes.** RK4, the elimination and the backward-Euler loop are + − × ÷ and √ only: class B, and every golden row is
  0 ulp. dopri45's step factor 0.9·err^(−0.2) and the Box–Muller ln and cos are transcendental, which this ADR sent to
  class D; under ADR 0006 they are `src/mathx.cyr`'s correctly rounded `_km_cr_pow`, `_km_cr_ln` and `_km_cr_cos`, so
  `num::dopri45`, `num::euler_maruyama`, `num::Pcg32::next_normal`, `simulate_adaptive` and `simulate_stochastic` are
  class C at **bound 0** against the corrected oracle, bit-identical too. glibc misrounds a few of those calls; where
  that changes a golden row's outputs (two pows in a 1,766-point Brusselator run at tol 10⁻¹² move every later point;
  four cos and one ln over three `next_normal` rows) the row gets an ADR 0006 correction
  (`tests/golden/numerics.cr.tsv`, `reaction_dynamics.cr.tsv`), against which the port is exact.
- **The golden rows** call hisab 1.4.0 directly (`scripts/oracle/src/m_numerics.rs`, the generator's one direct
  dependency besides the oracle, the same locked hisab) on the edges kimiya's guards never pass on: t₀ ≠ 0, backward
  integration, t_end = t₀, the empty state, h_init 0, negative, NaN or larger than the span, the 100 000-step limit, NaN
  and ±∞ pivots, pivots at 10⁻¹² ± 1 ulp, ties in the pivot search, Hilbert matrices. `tests/numerics.tcyr` replays them
  (`tests/parity.tcyr` skips the module); `validate.py` §26 recomputes every row with its own transcription.
- `num::least_squares_poly` (fitting) and `matrix_rank` (network) join `src/numerics.cyr` in M5b. The hisab feature
  request (system integrators, a dense solver, …) stands. The status stays **Accepted**.

## Amendment (M5b, 2026-11-12)

**The least squares is transcribed as decided; the rank and the CSR compression are transcribed too, because hisab 3.4.1
disagrees with 1.4.0 on the golden rows; `compute_gradient` runs on hisab 3.4.1's tape, as decided.** `src/numerics.cyr`
gains three private routines, hisab 1.4.0's code operation for operation:

- **`_km_least_squares_poly`** — `num::least_squares_poly` (`src/num/linalg.rs`): the Vandermonde columns by `powi`
  (`_km_powi`, compiler-builtins' loop), `qr_decompose`'s modified Gram–Schmidt (each dot product and norm a `.sum()`
  folded from −0.0, the 10⁻¹² test on the norm, a NaN norm passing), Qᵀy and back substitution. fit_polynomial and
  fit_beer_lambert are class B, fit_arrhenius (a ln per point, an exp of the intercept, `_km_cr_ln` / `_km_cr_exp`)
  class C at bound 0; every golden row is 0 ulp, and glibc 2.44 rounds every one of those calls correctly, so no row
  needs an ADR 0006 correction.
- **`_km_matrix_rank`** — `num::matrix_rank(a, Some(tol))` on `num::svd` (`src/num/svd.rs`): Aᵀ for a wide matrix,
  `svd_tall`'s one-sided Jacobi sweeps on the vectors, each singular value its vector's norm, counted above 10⁻¹⁰ (V, U,
  the basis extension and the sort feed nothing the count reads and are left out). The decision above sent the rank to
  hisab 3.4.1's `svd_compute_checked`, transcribing only if a golden row disagreed. Measured on the 563 golden rank rows
  (a scratch probe through `svd_compute_checked`, the wide ones transposed, since 3.4.1 refuses m < n): it disagrees on
  **191** — on 159 hisab 1.4.0's Jacobi ends in NoConvergence and 3.4.1 answers, on 16 1.4.0 answers and 3.4.1 refuses
  (a non-finite entry, which 1.4.0 counts in a single column's or row's norm, or a singular value past DBL_MAX), and on
  16 both answer with different ranks (presets scaled by 10¹⁰⁰ … 10²⁰⁰, whose unscaled sums of squares overflow in
  1.4.0, which then counts a full rank).
  The NoConvergence rows are 1.4.0's own behaviour: its test |a_pq| ≤ EPSILON_F64²·√(a_pp·a_qq), 10⁻²⁴ relative, lies
  below the rounding of a rotated pair's dot product, so a dense matrix rarely passes it in 100·max(m, n) sweeps —
  Chapman's ozone cycle, Robertson's kinetics, the H₂ + Br₂ chain and about six in ten seeded matrices are kimiya
  1.2.0's ComputationError. The transcription reproduces every row (class A); the defect is the oracle's, recorded in
  the roadmap's backlog for a decision after the port.
- **`_km_csr_from_dense`** — `CsrMatrix::from_dense` (`src/num/sparse.rs`): the entries with |v| > EPSILON_F64 (10⁻¹²),
  row by row, into a hisab 3.4.1 CSR handle (`csr_new`), so `sparse_stoichiometric_matrix` still returns hisab's type.
  The decision above sent it to `csr_from_dense`, which keeps every nonzero entry (a structural test since hisab
  2.15.0): it would store a 10⁻¹³ coefficient and a NaN that 1.4.0 drops (golden rows at the threshold to the ulp).

**autodiff.** `compute_gradient` builds on hisab 3.4.1's tape (`ad_tape_new` / `ad_var` / `ad_add` … / `ad_grad`), the
callback f(ctx, tape, vars, n) recording its function from the parameter nodes. Its 124 golden rows (20 test functions
at physical and seeded points: the fixed formulas' tapes, van 't Hoff, Eyring, Lennard-Jones, Morse, Rosenbrock, …)
agree with 1.4.0's gradients at class D (relative ≤ 10⁻¹², at most 40 ulp), so nothing was transcribed. The two tapes do
differ by design — 3.4.1 folds constants, compensates the sweep's sums (TwoSum), gives a domain miss (ln of x ≤ 0, a
divisor below DBL_MIN) a NaN partial and a 0 value where 1.4.0 has ±∞ or NaN, and takes its exp, ln, sin, cos and pow
from the stdlib — so the rows sit inside every function's domain; the four fixed-formula gradients are closed forms
(ADR 0002 #2's M5b addition), not tapes.

The hisab feature request stands, extended: a Jacobi SVD (or `matrix_rank`) with a reachable convergence test, polynomial
least squares, and a `from_dense` with a threshold. The status stays **Accepted**.
