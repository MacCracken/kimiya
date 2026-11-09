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

## Amendment (M3a, 2026-11-09)

**Two more ill-conditioned functions, and log10 is composed, not correctly rounded.** M3a ported `solution` and
`electrochemistry`. Two of their class-C functions cancel by construction, and the golden grids include
near-cancellation rows on purpose (port-parity.md §7):
- `solution::henderson_hasselbalch`, pH = pKa + log10([A⁻]/[HA]), cancels where pKa ≈ −log10(ratio) (pH ≈ 0);
- `electrochemistry::nernst_potential` / `nernst_potential_25c`, E = E° − (RT/nF)·ln Q, cancel where
  E° ≈ (RT/nF)·ln Q (E ≈ 0).

Measured by the rule (`validate.py --sensitivity`, every result of the stdlib-like primitive nudged ±1 ulp), before
any remedy:

| Function | Primitive | S before | Where |
|---|---|---|---|
| `henderson_hasselbalch` | `_km_log10` (fdlibm, within 1 ulp of glibc) | 4.4 × 10¹⁸ | a pH of a few ulp of the log either side of 0 moves to the next representable value off zero: the ulp line from +0 |
| `nernst_potential` | `f64_ln` | 8.8 × 10¹⁸ | the same, for E ≈ 0 |
| `nernst_potential_25c` | `f64_ln` | 8.7 × 10¹⁸ | the same, at 298.15 K |

The replay agreed: with `f64_ln`, aarch64 (the polyfill) has 6 Nernst rows up to 4.3 × 10¹⁸ ulp from the oracle; x86_64
(x87 `fyl2x`) happened to match glibc on every golden argument.

**Nernst: the decision above, unchanged.** The Nernst functions call the correctly rounded `_km_ln`. No golden `ln`
argument of theirs is one where glibc's log is not correctly rounded, so **S = 0, c = 0, K = 2**, and all 2,536 value
rows are 0 ulp on x86_64 ELF, x86_64 PE and aarch64.

**henderson_hasselbalch: a correctly rounded log10 was tried, and fails the grid.** `_km_log10_cr` (`_km_ln` in
double-double times 1/ln 10 in double-double, one rounding) put 6 near-cancellation rows up to 4.4 × 10¹⁸ ulp from the
oracle: at base/acid = 0.8068…/2.6858… glibc's log10 is not correctly rounded (`…69da`, correctly rounded `…69db`), and
a pH that should be the oracle's +0.0 became −1.1 × 10⁻¹⁶. Unlike acos, cos and ln (0.04–0.16 % of arguments), glibc
2.44's log10 is often not correctly rounded: on 22 of 8,000 binade-uniform arguments (0.3 %), 149 of 8,000 between
2⁻⁸⁰ and 2⁸⁰ (1.9 %), 2,340 of 8,000 in [½, 1) (29 %: for x in [½, 2) it is ivln10·log x, rounded twice) and 69 of the
1,497 positive finite golden arguments (4.6 %). A correctly rounded log10 cannot meet the rule on a grid with
near-cancellation rows, and the rule forbids a wider K or a thinner grid.

glibc's log10 is small enough to reproduce exactly. It is Sun's `sysdeps/ieee754/dbl-64/e_log10.c` composition over
glibc's own log: x = 2ᵏ·x′ (x′ in [1, 2) for k ≥ 0, [½, 1) for k < 0), y = k (+1 when k < 0), and
log10(x) = (y·log10_2lo + ivln10·log x′) + y·log10_2hi, with no FMA contraction. Recomputed in Python with glibc's log,
that formula gives glibc's log10 bit for bit on 500,000 seeded arguments and every golden one. So the port's log10 for
`henderson_hasselbalch` is **`_km_log10_ln`**, that composition over the correctly rounded `_km_ln` (`src/mathx.cyr`;
+ − × ÷ and `_km_ln` only, so the same bits on every target). It returns glibc's bits wherever glibc's log at x′ is
correctly rounded. S counts only the calls where it is not: on the golden grid there are none, so **S = 0, c = 0,
K = 2**, and all 831 value rows are 0 ulp on every target.

Verification of `_km_log10_ln` (a scratch probe, not the golden grid): 1,787 golden arguments (the 1,296 `math` log10
rows and the henderson ratios) plus 23,710 seeded ones (xorshift64 over every binade, [½, 1), and 2⁻⁸⁰ … 2⁸⁰): the same
bits on x86_64 ELF, x86_64 PE and aarch64 (apart from the sign of the NaN for a negative argument, which the parity
rules don't compare); it differs from glibc on 13 seeded arguments, and on each of those glibc's log at x′ is not
correctly rounded (the 70-digit reference decides), which is exactly the set the construction predicts. The `math`
corpus shadow line `f64::log10 ~ _km_log10_ln` is 0 ulp on all 1,296 rows, the 48 where glibc is not correctly rounded
included.

`_km_log10` (fdlibm) stays wherever S is small: `ph_from_h_concentration`, `poh_from_oh_concentration`, `weak_acid_ph`
and `weak_base_ph` have S = 1 (K = 4), and `ph_from_h` is a benched path. One primitive does not serve both well:
`_km_log10_ln` costs 1.85 µs a call (its `_km_ln`), `_km_log10` 59 ns.

**Consequences.**
- **Positive** — the Henderson–Hasselbalch and Nernst values are bit-identical to the oracle on all 3,367 golden value
  rows, the near-cancellation ones included, on all three targets.
- **Negative** — `_km_log10_ln` is not correctly rounded: it carries glibc's log10 error (up to about 1.5 ulp), by
  design, because the oracle has it. This is the first port primitive that reproduces glibc's composition rather than
  the mathematical value, and it depends on glibc 2.44's e_log10.c being the Sun formula (`validate.py` checks every
  golden argument; re-measure on any glibc change). The correctly rounded `_km_ln` is slow: `nernst_potential` is
  1.94 µs against Rust's 6.42 ns (about 45 ns with `f64_ln`); a two-phase (Ziv) correctly rounded log, fast path first,
  would recover most of it and is a possible follow-up.
- **Neutral** — as before, a consumer's input can be both ill-conditioned and an argument where glibc's log is not
  correctly rounded; there the port is many ulp from the oracle (but no further from the true value).

The status stays **Accepted**.
