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
