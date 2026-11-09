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
