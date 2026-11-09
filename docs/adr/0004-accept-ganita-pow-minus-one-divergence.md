# 0004 — Accept ganita's pow(−1, DBL_MAX) divergence

**Status**: Accepted
**Date**: 2026-11-08

## Context

The M0 `math` golden module records Rust std `f64::powf` (glibc 2.44) over special values.
For `powf(-1.0, 1.7976931348623157e308)` glibc returns **1.0** (DBL_MAX is an even integer);
ganita 1.2.15's `ganita_f64_pow` returns **+∞**. It is ganita's |x| == 1 defect for huge even
exponents, already filed upstream by tanmatra
(`ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`) and
accepted there by tanmatra ADR 0007. It is the only class-C row beyond K = 2 ulp on any target
(x86_64 ELF, x86_64 PE, aarch64; port-parity.md §2).

## Decision

Accept that one row: `tests/parity.tcyr` `_par_accepted` lists `f64::powf` with input
`f:bff0000000000000 f:7fefffffffffffff` under ADR 4. It is reported on every run and does not
fail its module. No other row is accepted, and the K = 2 bound is not widened.

## Consequences

- **Positive** — the bound stays tight for every other `powf` row.
- **Negative** — none for kimiya's results: every kimiya `powf` call has a positive base
  (10^x for pH and Antoine, ratios of positive temperatures/pressures, positive concentrations).
- **Neutral** — remove the entry when ganita fixes the defect and the pin moves.

## Alternatives considered

- **Drop the row from the grid**: rejected — the row documents a real primitive divergence.
- **Transcribe glibc's pow**: unnecessary — no kimiya input reaches the defect.
