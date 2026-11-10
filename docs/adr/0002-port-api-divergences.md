# 0002 — Accept the planned API changes in the Cyrius 2.0 API

**Status**: Accepted
**Date**: 2026-11-08

## Context

The 2026-11-08 audit found eleven defects whose correct fix changes a public signature. kimiya
1.2.0 documents them but keeps the 1.x API. On 2026-11-08 the maintainer decided there will be
no Rust 2.0: the changes are made in the Cyrius API only. No published dependent (jivanu,
khanij, kiran, jnana, tara) calls any of these functions.

Under ADR 0001 a divergence from the oracle is a defect unless an ADR accepts it. This ADR
accepts these, and fixes how the parity harness still compares each one, so that no function
leaves parity entirely.

## Decision

The Cyrius 2.0.0 API differs from Rust 1.2.0 in exactly these places. Each row's golden rows
are compared through the stated mapping:

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `thermochem::reaction_enthalpy`, `reaction_gibbs_energy` return kJ/mol | return **J/mol** | oracle × 1000, computed as `(Σ kJ) × 1000.0` on both sides — class B, expected 0 ulp |
| 2 | `sensitivity::{equilibrium_constant_gradient, arrhenius_gradient, nernst_gradient, gibbs_temperature_sensitivity}` unvalidated; `compute_gradient` → `Result` never `Err` | the four return `Result`, computed from closed-form derivatives; `compute_gradient` returns its vector directly | finite oracle rows class D; oracle non-finite ↔ Cyrius `Err` |
| 3 | `nuclear::semi_empirical_binding_energy` → NaN for z > a or a = 0 | `Result` | finite rows class B; oracle NaN ↔ `Err` |
| 4 | `phase::antoine_pressure_mmhg`, `antoine_pressure_pa` → NaN at the pole, silent extrapolation | `Result`; `Err` outside the coefficient set's `t_min_c..t_max_c`. No extrapolating variant | in-range rows class C; oracle NaN or out-of-range input ↔ `Err` |
| 6 | `FunctionalGroup::hydrogen_bonding`; `Geometry::TShapedThreeLonePairs` | `functional_group_is_hydrogen_bond_donor` and a new `functional_group_is_hydrogen_bond_acceptor`; `KM_GEOMETRY_T_SHAPED` | donor and geometry rows class A under the new names; acceptor checked by reference tests only |
| 7 | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational oxidation state (numerator / positive denominator, reduced); balancing takes an integer matrix | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer average oxidation state" `Err` ↔ `Ok(p/q)`, verified by reference tests (Fe₃O₄ +8/3, C₃H₈ −8/3); balancing class A after converting the oracle's integer-valued f64 inputs |
| 12 *(2026-11-09)* | `inorganic::cfse_octahedral` computes `t2g · (−0.4) + eg · 0.6`, which misses the exact (3·e_g − 2·t₂g)/5 on 13 of its 22 (d, spin) cases (high-spin d⁵ −2⁻⁵² and d¹⁰ −2⁻⁵¹ instead of 0, d³ −1.2000000000000002) | the exact value: the integer numerator 3·e_g − 2·t₂g converted to f64, then one correctly rounded division by 5 | class A. Where the oracle's value equals the exact one (9 cases), a plain row; elsewhere the golden row keeps the oracle's value, its note names the exact one, and the harness maps the two (q = round(5 × oracle) must be within 10⁻⁹ of it, and the port must return q / 5), reporting each such row by name as accepted under ADR 2 |
| 13 *(2026-11-09, M4c; completed in M4d)* | `nuclear::atoms_after_half_lives` computes `initial * 2.0.powf(-n)`: where 2⁻ⁿ is not a normal double (n > 1022 or n ≤ −1024) it overflows to ∞, underflows to 0 (n = 1075 is a tie, to 0) or is a subnormal of 52 bits or fewer, so the result is ±∞, 0, NaN or up to 2⁵² ulp off where N₀/2ⁿ is representable (1000 after 1075 half-lives: 0, not 2.47 × 10⁻³²¹; 2¹⁰⁰ after 1,100: 0, not 2⁻¹⁰⁰⁰) | there, N₀/2ⁿ **correctly rounded** (rounded once, subnormal results at the 2⁻¹⁰⁷⁴ grid, ties to even): since M4d `src/mathx.cyr`'s `_km_cr_exp2_scaled(N₀, −n)` for every N₀ and n, 2⁻ᶠ from exp's machinery past 53 bits with a rounding test and a triple-double fallback (M4c's version was faithful for a fractional n with N₀ not a power of two). ±0, ±∞, NaN N₀ give N₀ (scaleB). Rust's form stays elsewhere | class C; equal rows plain; the others keep the oracle's value with the port's in the note (`(f:…)`), which the harness requires bit for bit and reports as accepted under ADR 2, counted (1,810 rows) |
| 14 *(2026-11-09, M4c)* | `nuclear::is_secular_equilibrium` tests `λ_p < 0.01·λ_d`: 0.01 is not 1/100 and the product rounds or underflows, so (0.03, 3) is false and (0, 5 × 10⁻³²⁴) false | the exact λ_p < λ_d/100: 100·λ_p against λ_d by an exact double-double (`_km_two_prod` on the significands, the exponents first), Rust's NaN (false), ±∞ and ±0 semantics kept | class A; equal rows plain; the others noted `(i:…)`, mapped and reported by name (208 rows) |
| 15 *(2026-11-09, M4b's inorganic, decided in M4c)* | `inorganic::delta_tet_from_oct` computes `delta_oct * 4.0 / 9.0`, correctly rounded except where Δ·4 overflows (\|Δ\| ≥ 2¹⁰²², above DBL_MAX/4): ±∞ where 4Δ/9 is finite | there, (Δ/9)·4: one rounded division and an exact ×4 that can't overflow, which is again the correctly rounded 4Δ/9 (4·DBL_MAX/9 ≈ 7.99 × 10³⁰⁷). Rust's form elsewhere | class B; noted `(f:…)`, mapped and reported by name (59 rows) |
| 16 *(2026-11-09, M4d)* | `potential::lj_equilibrium_distance` computes `sigma * 2.0.powf(1.0 / 6.0)`: the exponent 1/6 is rounded (a simple exact constant misrepresented in binary) and the product rounds again, so r_min = σ·2^(1/6) is one ulp off on about a fifth of all σ (4,117 of 20,000 seeded). (2^RN(1/6) rounds to RN(2^(1/6)), MPFR checks; the misses are the second rounding) | the correctly rounded σ·2^(1/6), rounded once: σ = s·2^e, s·(H + M + L) by `_km_td_mul_d` (2^(1/6) as a triple-double from `mathx_tables.py`'s exact integer root), then `_km_cr_finish` at the final scale; no rounding test is needed, since 2^(1/6)'s continued fraction proves every σ·2^(1/6) at least 2^-55.1 ulp from a rounding boundary, far beyond the product's 2^-102.6. NaN, ±0, ±∞ give σ | class C (Rust's `powf`); equal rows plain; the others noted `(f:…)`, mapped and counted (1,041 rows) |
| 17 *(2026-11-09, M4d)* | `potential::lennard_jones` computes `4.0 * epsilon * (sr12 - sr6)`, i.e. (4·ε)·d: for ε ≥ 2¹⁰²² the exact 4·ε overflows, so the result is ±∞ (or ∞·0 = NaN at r = σ) where 4·ε·d is representable | there, (ε·d)·4: one rounded product and an exact ×4, which is RN(4·ε·d), the value Rust's (4·ε)·d has wherever 4·ε is finite (4·ε is exact). Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (180 rows) |
| 18 *(2026-11-09, M4d)* | `potential::morse_width_parameter` computes `omega_e * (reduced_mass / (2.0 * d_e)).sqrt()`: for D_e ≥ 2¹⁰²³ the exact 2·D_e overflows, μ/∞ = 0 and a = 0 where ω·√(μ/(2·D_e)) is representable | there, (μ·½)/D_e: RN(μ/(2·D_e)) (μ·½ is exact, or μ < 2⁻¹⁰²¹ and the quotient rounds to 0 either way), the value Rust's μ/(2·D_e) has wherever 2·D_e is finite; then Rust's √ and ×ω. Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (181 rows) |

Deferred with their modules (port plan §2.2), and decided when they port: #5 `bridge::*` →
`Result` and removal of `ionization_to_bond_energy_estimate`; #10 soorat `from_molecule` →
`Result`.

Not in 2.0 (additive, demand-gated, may land in either implementation): #8 Beer–Lambert
intercept, #9 Rydberg nuclear-mass parameter, #11 full `simulate_stiff` trajectory.

The parity harness implements each mapping in one table and reports these rows by name.

## Addition, 2026-11-09 (M4b): #12 and the policy for inexact formulas

**#12, the exact CFSE.** Rust's `cfse_octahedral` writes the crystal field stabilization energy −(2/5)·n(t₂g) +
(3/5)·n(e_g) with the decimal coefficients −0.4 and 0.6, neither exactly representable, and two roundings in the sum.
The quantity is an exact rational, k/5 for the integer k = 3·n(e_g) − 2·n(t₂g), and on 13 of the 22 Ok cases Rust's
bits are not the double nearest it: the textbook zeros of high-spin d⁵ and of d¹⁰ come out as −2.2 × 10⁻¹⁶ and
−4.4 × 10⁻¹⁶, and −6/5 as −1.2000000000000002. The user's direction for the port is the best accuracy available, with
parity in function (ADR 0006: correctly rounded primitives, the corrected oracle as the reference). The same reasoning
applies to an exact-rational formula: the port computes k exactly in integers and returns CR(k/5), one division,
so it is the correctly rounded exact value on every input; the oracle's errors are at most a few ulps, but they turn
exact zeros into nonzero values. Functional parity is unchanged (the same table, the same `Err` codes and check order).

**Golden comparison.** The corpus keeps the oracle's values (it is regenerated only from Rust 1.2.0). Each of the 13
rows whose exact value differs carries it in its note; the harness recovers q = round(5 × oracle), requires the oracle
within 10⁻⁹ of q/5 and the port's value bit-identical to q/5, and counts the row as accepted under ADR 2, printing it by
name. q comes from the oracle, not from the port's tables, so a wrong occupancy or formula still fails. The rows are
class A: the value is exact, so no tolerance applies. `validate.py` checks both the oracle's value (Rust's formula) and
the exact one (exact rational arithmetic) on every row, and the notes.

**Policy (the user's decision, the same day).** Other formulas in M4 whose Rust rounding misses a simple exact
rational or a chemically exact value get the same treatment, each as its own row here, accepted before it lands. A
formula that is already correctly rounded, or whose result is not such a value (a product of rounded constants, a sum of
arbitrary inputs), keeps Rust's operations.

**#6's acceptor criterion (the user's decision, the same day).** `functional_group_is_hydrogen_bond_acceptor` follows
IUPAC 2011's criterion — an electron-rich region such as a lone pair or a π-bonded pair — weak acceptors included, so
every group but the alkane accepts; its sources, one per weak family, are in `src/organic_data.cyr`.

## Addition, 2026-11-09 (M4c): #13, #14 and #15

**The user's decision.** M4c's review under the policy above found two more formulas whose Rust rounding misses an exact
value, `nuclear::atoms_after_half_lives` and `nuclear::is_secular_equilibrium`, and M4b had left one,
`inorganic::delta_tet_from_oct`'s overflow, as a question. The user decided all three the same day: **fix exactly**.

**#13, N₀/2ⁿ.** The value N₀/2ⁿ is exact for an integer n, and Rust's `N₀ · 2.0.powf(−n)` is its correct rounding
wherever 2⁻ⁿ is a normal double (validate.py proves it for every such golden row). Outside that range — n > 1022 or
n ≤ −1024 — Rust's intermediate 2⁻ⁿ overflows, underflows or loses bits as a subnormal, and the answer can be ±∞, 0, NaN
or wrong by up to 2⁵² ulp while N₀/2ⁿ is a perfectly ordinary double. There the port computes it with one rounding: N₀ =
s·2^j (s in [1, 2)) and n = k + f (k = ⌊n⌋), s·CR(2⁻ᶠ) as an exact double-double, then 2^(j − k) applied with a single
final rounding by `src/mathx.cyr`'s `_km_cr_finish` (exact scaling for normal results, one rounding at the subnormal grid,
ties to even, IEEE overflow). For an integer n, CR(2⁻ᶠ) = 1 and the result is the correctly rounded N₀/2ⁿ; for N₀ = ±2^j
with j − n exact it is ±CR(2^(j − n)) from `_km_cr_pow`, also correctly rounded. For any other N₀ and a fractional n the
result is **faithful, not correctly rounded**: CR(2⁻ᶠ) is the one inexact step (relative error ≤ 2⁻⁵³) before the final
rounding, so the error is below one ulp (on the golden rows correctly rounded 743 times in 1,121, nearer N₀/2ⁿ than Rust's
value on 1,087 and farther on 6, by one ulp, where Rust's own two roundings happen to land right). Correct rounding there
would need 2⁻ᶠ beyond 53 bits, a correctly rounded scaled exp2 that `src/mathx.cyr` does not offer; it is noted as a
follow-up. ±0, ±∞ and NaN N₀ return N₀ there (scaleB's semantics: 0 after −1100 half-lives is 0, where Rust's 0·∞ is
NaN), and an infinite n keeps Rust's limits.

**#14, λ_p < λ_d/100.** Rust's 0.01 is 0.01000000000000000020816…, and 0.01·λ_d rounds (or underflows to 0), so its test
disagrees with λ_p < λ_d/100 near the boundary: 66 rows of M4c's first corpus, among them (0.03, 3) and (0, 5 × 10⁻³²⁴).
The port compares 100·λ_p with λ_d exactly: NaN gives false and an infinite operand compares as itself (Rust's results,
since ∞/100 = ∞), zeros and signs decide first, and for two same-signed finite values the significands are compared
after the exponents (100·s_p as `_km_two_prod`'s exact double-double, on [1, 2) so neither overflow nor an underflowing
error term can occur, then p and the sign of e). The threshold remains a heuristic for "≪"; the comparison is now exact.

**#15, 4Δ/9.** Δ·4 is exact unless it overflows, so Rust's (Δ·4)/9 is correctly rounded except for |Δ| ≥ 2¹⁰²², where it
is ±∞. There the port computes (Δ/9)·4: one correctly rounded division and an exact ×4 (|Δ/9|·4 < DBL_MAX), so CR(4Δ/9)
on the whole range.

**Golden comparison** (as #12): the corpus keeps the oracle's values. Where the port's differs, the row's note names it —
`[ADR 0002 #13: the port returns <how> = <x> (f:<bits>)]`, `#14 … = true (i:1)]`, `#15 … (f:<bits>)]` — with the generator
computing it independently of the port (integer arithmetic on the significands for #13 and #14, Rust's (Δ/9)·4 for #15),
and `validate.py` proving every note with exact rational arithmetic (`Fraction`, MPFR for CR(2⁻ᶠ) and CR(2^(j − n))) and
checking that the notes sit exactly where the two values differ. #13's and #15's exact value can't be recovered from the
oracle's (a 0, ±∞ or NaN carries none), so the harness reads it from the note: the port's value must be that one, bit for
bit, and the row is accepted under ADR 2; any other value fails. 1,812 #13 rows, 208 #14 rows and 59 #15 rows are
mapped; the rest are plain ok. The harness counts the mapped rows and prints one summary line per item after its module
(as #4's Antoine rows), not a line per row: the notes are the record, and `validate.py` proves each one.

## Addition, 2026-11-09 (M4d): #13 completed; #16, #17 and #18

**#13, correctly rounded.** M4c left `atoms_after_half_lives`' N₀/2ⁿ faithful for a fractional n with N₀ not a power of
two (s·CR(2⁻ᶠ), then a scaled rounding: two roundings). M4d makes it correctly rounded with one primitive,
`src/mathx.cyr`'s `_km_cr_exp2_scaled(x, y)` = RN(x·2^y): x = s·2^j, y = k + f (k = trunc y); for f = 0 the exact
s·2^(j + k) rounded once by `_km_cr_finish`; otherwise 2^f from exp's fast core (f·ln 2 as a double-double), s·2^f as an
exact product plus a rounded tail (within 1.8·2⁻⁷⁰), exp's rounding test at the final scale (`_km_exp_fin`, normal or
subnormal), and when it fails the triple-double accurate path (exp's accurate core, `_km_td_mul_d`, `_km_cr_finish`;
below 2⁻¹⁵⁰). 2^f is irrational for every double f in (−1, 1) other than 0, so no value is a tie. The golden notes are now
the generator's correctly rounded N₀/2ⁿ (`scripts/oracle/src/exact.rs`: 2⁻ᶠ to 320 bits in integer arithmetic, Ziv's
two-sided test), `validate.py` proves each against MPFR (`mpfr_ref.scaled_exp2`: 2^f to 400 bits, bounded from both
sides) and a 60-digit 2⁻ᶠ with exact rational rounding, and `cr_verify.py` checks `_km_cr_exp2_scaled` against MPFR on
611,518 seeded pairs (100,228 through the accurate path alone) on x86_64 ELF, x86_64 PE and aarch64. 1,810 golden rows
differ from the oracle (564 integer n, 173 N₀ = ±2^j, 1,067 fractional n, 6 scaleB); before, 1,812.

**#16, #17, #18 (M4d's review of `potential` under the policy above).** `lj_equilibrium_distance` is case (a): Rust's
`2.0.powf(1.0 / 6.0)` raises 2 to the rounded 1/6 and the product rounds again; the port returns the correctly rounded
σ·2^(1/6). `lennard_jones`' (4·ε)·d and `morse_width_parameter`'s μ/(2·D_e) are #15's shape, case (b): an exact constant
multiple of an input (4·ε, 2·D_e) whose spurious overflow at the top of the range turns a representable result into
±∞, NaN or 0; there the port reorders exactly ((ε·d)·4, (μ·½)/D_e), which gives the bits Rust's form gives wherever it
does not overflow. **Reviewed and kept as Rust's** (no single exact value to restore, so Rust's sequence is the
definition; reported to the user): the over- and underflow of (σ/r)⁶ and (σ/r)¹² in `lennard_jones` (σ/r past 2^85.3 or
below 2^-89.5), of k·q₁·q₂ in `coulomb` / `coulomb_force` and of r·r in `coulomb_force` (r past 1.3 × 10¹⁵⁴ or below
1.5 × 10⁻¹⁵⁴), of the squared coordinates in `distance_3d`, `bond_angle` and `dihedral_angle` (beyond about 10¹⁵⁴, 10⁷⁷
for the dihedral's fourth powers, or below their reciprocals: an Err, 0, ±∞ or a few-bit result where the angle or
distance is representable, though Rust documents the angles as scale-invariant), `morse`'s (r − r_e) overflow and its
subnormal D_e·x, and the two typed tolerances (10⁻¹², 10⁻²⁰: numerical resolution thresholds, kept like
`nuclear::EQUAL_RATE_REL_TOL`).

**Golden comparison** (as #13–#15): the notes `[ADR 0002 #16: the port returns the correctly rounded σ·2^(1/6) = x
(f:…)]`, `#17 … (ε·d)·4, the correctly rounded 4·ε·d …`, `#18 … ω·√((μ·½)/D_e), with the correctly rounded μ/(2·D_e) …`
sit exactly where the port's value differs; the generator computes #16 with exact integer sixth powers (`exact.rs`) and
#17 / #18 as the reorders; `validate.py` proves #16 with rational arithmetic (x·2^(1/6) > t ⟺ 2x⁶ > t⁶, each midpoint
compared), #17 as RN(4·ε·d) and #18 as RN(μ/(2·D_e)) with `Fraction`, and Rust's forms as those values wherever 4·ε and
2·D_e are finite. The harness maps 1,041 #16, 180 #17 and 181 #18 rows and prints one summary line each after the
module.

## Consequences

- **Positive** — the Cyrius API is internally consistent (SI energies, `Result` for every
  fallible function, chemically meaningful names) from its first release.
- **Negative** — the Rust and Cyrius APIs intentionally differ on these functions; a consumer
  moving from Rust to Cyrius must apply them (2.0.0's CHANGELOG carries the migration notes).
- **Neutral** — the mapped comparisons are slightly more harness code than 1:1 rows.

## Alternatives considered

- **Rust 2.0 first, then a 1:1 port**: offered 2026-11-08, declined by the maintainer.
- **Port the 1.2 API unchanged**: rejected — it would carry known defects (kJ/J mix, NaN
  sentinels, misleading names) into a new codebase.
