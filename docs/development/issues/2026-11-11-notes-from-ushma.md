# Notes for kimiya from ushma's port (2026-11-11)

Two notes from ushma's M0, for kimiya's port to weigh (`ushma/docs/development/port-parity.md` §3,
`ushma-port-plan.md` §8). Neither is a defect report. This directory follows the issues layout the other AGNOS
repositories use (`docs/development/issues/`, `YYYY-MM-DD-slug.md`).

## 1. ganita's LU reproduces Rust hisab 1.4.0's Gaussian elimination, refusals included

kimiya's ADR 0003 plans to transcribe Rust hisab's `num::gaussian_elimination`. ushma's probe 9 compared hisab 1.4.0's
routine with ganita 1.2.15's two solvers on 2,178 seeded systems: row-scaled, unit diagonal, |aᵢⱼ| ≤ 1.
- **`ganita_mat_lu` + `ganita_mat_lu_solve`:** bit-identical on all 2,000 systems hisab solves (37,424 unknowns),
  and refuses exactly the 178 it refuses. Its threshold is 1e-12 · max|A|, which on such matrices is hisab's
  absolute `|pivot| < 1e-12`.
- **`ganita_mat_gaussian_elim`:** refuses 259 of the 2,000. Its threshold also spans the right-hand side (filed
  upstream: ganita `2026-11-11-gaussian-elim-pivot-scale-includes-rhs.md`).

If kimiya's matrices are scaled so that max|A| = 1, LU gives hisab 1.4.0's bits without a transcription. Otherwise
the thresholds differ by that scale. `ganita_mat_lu_solve` allocates 8·n B per call (ganita
`2026-11-11-mat-lu-solve-allocates-per-call.md`).

## 2. `KM_GAS_CONSTANT` is 8.314462618, ushma's R is the exact N_A·k_B

kimiya's `KM_GAS_CONSTANT` is 8.314462618, as its oracle's `element::GAS_CONSTANT` is ("stored to 10 significant
figures"), so it is right for kimiya's parity. ushma 1.4.0 uses the exact N_A·k_B = 8.31446261815324 J/(mol·K). The
two differ by 1.8e-11 relative, in the 11th significant digit, so a value that crosses the ushma ↔ kimiya bridges
(kimiya `bridge.rs`, ushma `bridge.rs`) is computed with two different R. kiran uses both.

This is only information. Whether kimiya keeps its oracle's constant or takes the exact SI value under its accuracy
policy (its ADR 0006) is kimiya's decision. If it changes, ushma's bridge rows would note the agreement.
