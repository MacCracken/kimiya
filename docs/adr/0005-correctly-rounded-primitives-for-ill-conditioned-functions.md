# 0005 — Correctly rounded acos, cos, ln and x^2.5 for the ill-conditioned gas functions

**Status**: Accepted
**Date**: 2026-11-09

## Context

The port-parity rule (`docs/development/port-parity.md`, from tanmatra) sets a class-C function's
bound before its first comparison: S is the largest output change, in ulp, that a 1-ulp change in one
of its class-C primitive results causes over its golden grid, and K_F = 2·S·c + 2 (c the class-C calls
on its path). If K_F > 64 the function is ill-conditioned on that grid, and the rule says: make the
primitive bit-identical to the oracle's, or, if that isn't possible, write an ADR.

Measured in M2a (`scripts/oracle/validate.py . --sensitivity`, every golden row, each glibc result
nudged ±1 ulp) with the primitives as the oracle calls them:

| Function | S | Where |
|---|---|---|
| `gas::peng_robinson_z_roots` | 2²⁶ (and 14 root-count flips) | near a double root of the cubic: two roots 10⁻⁸ apart, conditioned like √ε; the flips at P = DBL_MIN, where acos(−1) moves the double root at 0 across Z = B |
| `gas::fugacity_coefficient_pr`, `fugacity_pr` | 318, 492 | φ = exp(ln φ): a 1-ulp change in ln(Z − B) moves φ by |ln φ| ulp, and |ln φ| reaches hundreds at extreme states (N₂ at 2 K) |
| `gas::redlich_kwong_pressure` | 64 | RT/(V − b) − a/(√T·V(V + b)) cancels near P = 0 |

The oracle's primitives are glibc 2.44's acos, cos, log and pow (Rust std). They cannot be transcribed
bit for bit: glibc is not in the oracle's dependency graph (no source to transcribe from, unlike
tanmatra's pure-Rust libm), and on x86_64 it selects FMA variants at load time. glibc is not correctly
rounded either: over 20,000 random arguments each, it is 1 ulp off the correctly rounded value on
0.07 % (acos), 0.16 % (cos), 0.04 % (ln) and 0.08 % (pow(x, 2.5)) of them.

## Decision

The port evaluates those four **correctly rounded** — `_km_acos`, `_km_cos` (for |x| ≤ 8, the
Peng–Robinson arguments being in [−2π/3, π]), `_km_ln` and `_km_pow_2_5` in `src/mathx.cyr`, double-double
evaluation (≈ 2⁻¹⁰⁰ relative) with one final rounding, + − × ÷ and `f64_sqrt` only, so the same bits on
every target — and the Peng–Robinson and Redlich–Kwong functions call them in place of ganita's acos,
`f64_cos`, `f64_ln` and ganita's pow. exp keeps `f64_exp` and 10^x keeps ganita's pow: their functions are
well-conditioned.

The port's result is then the oracle's, bit for bit, wherever glibc's result is correctly rounded, so the
rule's S counts **only the primitive calls whose result can differ between the port and the oracle**:
- every exp and pow(10, ·) result (a stdlib counterpart, within 1 ulp of glibc);
- the acos, cos, ln and x^2.5 results where glibc is not correctly rounded — `validate.py` decides that
  per golden argument against its 60-digit reference;
- never a result that is 0, ±∞ or NaN, or acos(±1) (exact on both sides).

With that, every M2a class-C function has K_F ≤ 10 (port-parity.md §5): `peng_robinson_z` 2,
`peng_robinson_z_roots` 6, `fugacity_coefficient_pr` 6, `fugacity_pr` 10, `redlich_kwong_pressure` 2. The
golden rows keep class C; the bounds are gated like every other.

## Consequences

- **Positive** — the Peng–Robinson roots and Redlich–Kwong pressures are bit-identical to the oracle on
  all 9,261 golden rows, ill-conditioned ones included, on x86_64 ELF, x86_64 PE and aarch64. The
  fugacities differ only through exp (≤ 2 ulp). The port's four primitives are target-independent and
  correctly rounded on their working ranges (x^2.5 for x ≥ 2⁻³⁸⁷; below, the double-double product's
  error term underflows). An independent check (its own 90-digit reference and 48,608 seeded arguments, not
  the golden grid) found no acos, cos or ln result off correct rounding, and the same bits on all three
  targets.
- **Negative** — on an argument where glibc is not correctly rounded, the port's primitive is 1 ulp
  away from the oracle's, and on an ill-conditioned state that 1 ulp is amplified by the state's own S.
  The golden grid has 5 such cos calls and 5 such ln calls, none on an ill-conditioned row; a consumer's
  input can be both. kimiya owns four extra primitives (about 330 lines) and their cost: the three-root
  branch makes one acos and three cos calls in double-double.
- **Neutral** — the decision depends on glibc 2.44's rounding: `validate.py --sensitivity` recomputes the
  set of differing calls, and `--check` already refuses a corpus from another glibc. Re-measure on any
  glibc or grid change, and revisit if glibc's acos, cos or log become correctly rounded (CORE-MATH).

## Alternatives considered

- **Bounds above 64**: refused by the rule — K = 2²⁸ ulp for the roots says nothing.
- **An ADR splitting out the ill-conditioned rows with per-row bounds** (the rule's second option, with
  ganita's and the stdlib's primitives): about 1,140 rows would have needed their own bound (60 root rows
  up to 2²⁸ ulp, 1,067 fugacity rows up to 10³), and root-count flips can't be bounded at all.
- **Transcribing glibc**: not possible — no source in the oracle's dependency graph, and the FMA
  ifunc variants (port plan §0).
- **Transcribing Rust libm's fdlibm acos/cos/log/pow**: those are not the oracle's either (ganita's acos,
  similar, differs from glibc on 28 % of the corpus rows).
