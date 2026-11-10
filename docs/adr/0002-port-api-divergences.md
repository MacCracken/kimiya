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
| 19 *(2026-11-09, M4e)* | `spectroscopy::wavenumber_from_wavelength` and `wavelength_from_wavenumber` compute `0.01 / x` for the exact 1/(100·x) (the cm ↔ m factor 1/100): 0.01 is 1/100 + 2.08·10⁻¹⁹ (relative 2⁻⁵⁵·⁴), and the quotient rounds again, so the result is one ulp from RN(1/(100·x)) on about 14 % of all doubles (2,748 of 20,000 seeded) | RN(1/(100·x)), rounded once: `src/mathx.cyr`'s `_km_cr_hundredth_div` (a double-double residual with 1/100 = H + L and a rounding test, an exact 64-bit long division as its fallback and past 2^±900); the special values Rust's | class B; noted `(f:…)`, mapped and counted (827 rows each) |
| 20 *(2026-11-09, M4e-2)* | `potential::distance_3d`, `bond_angle` and `dihedral_angle` square raw coordinates: above about 10¹⁵⁴ (10⁷⁷ for the dihedral's fourth powers) or below their reciprocals a square, cross product or dot product overflows or underflows, and the result is an Err, 0, ±∞ or a degraded value where the true distance or angle is representable (`distance_3d` of (DBL_MAX, 0, 0) and the origin is +∞, the gauche dihedral at 2⁷⁰⁰·Å an Err, the angle of (1, 2⁻¹⁰⁶⁰) 0), so the angles, which Rust documents as scale-invariant, are not | Rust's own operation sequence (every subtraction, product, sum, comparison, `max` and the Err checks, in Rust's order) evaluated in IEEE binary64 round-to-nearest-even with an **unbounded exponent range**: every + − × ÷ √ rounds to 53 bits and nothing overflows or underflows. The exception is the **final** operation (`sqrt` for `distance_3d`, `atan2` for the angles), which is correctly rounded to binary64 itself: its subnormal grid, or ±∞. Comparisons and Err branches use the unbounded values; non-finite inputs keep Rust's behaviour. So every row where Rust's chain stays in the normal range is bit-identical to Rust; the angles are exactly scale-invariant (for every k with all scaled coordinates exact, the same bits at 2^k times the coordinates: every check is homogeneous, and `atan2`'s two arguments have the same degree), and `distance_3d` scales by exactly 2^k wherever the result is normal. Rust's chain with a cheap range test, and `src/mathx.cyr`'s `_km_ux_*` helpers where it fails | class B (`distance_3d`) and C (the angles); noted `(f:…)` or `(E:InvalidInput)`, mapped and counted (854, 1,205 and 1,546 rows; 91 of the dihedral's are the port's Err where Rust's rounded threshold says Ok) |

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
distance is representable, though Rust documents the angles as scale-invariant; fixed in M4e-2 as #20, the user's
decision), `morse`'s (r − r_e) overflow and its
subnormal D_e·x, and the two typed tolerances (10⁻¹², 10⁻²⁰: numerical resolution thresholds, kept like
`nuclear::EQUAL_RATE_REL_TOL`).

**Golden comparison** (as #13–#15): the notes `[ADR 0002 #16: the port returns the correctly rounded σ·2^(1/6) = x
(f:…)]`, `#17 … (ε·d)·4, the correctly rounded 4·ε·d …`, `#18 … ω·√((μ·½)/D_e), with the correctly rounded μ/(2·D_e) …`
sit exactly where the port's value differs; the generator computes #16 with exact integer sixth powers (`exact.rs`) and
#17 / #18 as the reorders; `validate.py` proves #16 with rational arithmetic (x·2^(1/6) > t ⟺ 2x⁶ > t⁶, each midpoint
compared), #17 as RN(4·ε·d) and #18 as RN(μ/(2·D_e)) with `Fraction`, and Rust's forms as those values wherever 4·ε and
2·D_e are finite. The harness maps 1,041 #16, 180 #17 and 181 #18 rows and prints one summary line each after the
module.

## Addition, 2026-11-09 (M4e): #19

**#19, 1/(100·x) (M4e's review of `spectroscopy` under the policy above).** Rust's wavenumber conversions are
ν̃ = 1/(λ in cm) = (1/100)/λ_m and its inverse, written `0.01 / x`: case (a), a simple exact constant (1/100, the
centimetre) misrepresented in binary, then a second rounding. The port returns the correctly rounded 1/(100·x) for every
double x. `src/mathx.cyr`'s `_km_cr_hundredth_div`: for 2⁻⁹⁰⁰ ≤ |x| ≤ 2⁹⁰⁰, y0 = RN(H/x) with 1/100 = H + L + t
(`mathx_tables.py`: H = Rust's 0.01, |t| < 2⁻¹¹⁶), the residual 1/100 − y0·x from an exact TwoProduct to 2⁻¹¹⁰·², its
correction d = RN(r/x), and a rounding test with ε = 2⁻¹⁰⁰·y0 against the proven 2⁻¹⁰²·⁹ (it fails with probability
about 2⁻⁴⁷); otherwise, and everywhere outside that range, an exact restoring division of 2^(b + 53) by 100·m in 64-bit
integers (|x| = m·2^e, 100·m < 2⁶⁰), rounded once with its remainder as the sticky bit at the binade's or the subnormal
grid, ±∞ past DBL_MAX. 1/(100·x) is never a midpoint (100·x has the odd factor 25), so there are no ties. ±0 → ±∞,
±∞ → ±0 and NaN → NaN are Rust's. Rust's value is within one ulp of it on every normal result.

**Golden comparison** (as #16–#18): the notes `[ADR 0002 #19: the port returns the correctly rounded 1/(100·λ) = x
(f:…)]` (and `1/(100·ν̃)`) sit exactly where the port's value differs from the oracle's; the generator computes it in
128-bit integers (`exact::cr_hundredth_div`), `validate.py` proves every note with `Fraction`, `mpfr_ref.hundredth_div`
gives MPFR's value (100·x exact at 64 bits, one division, binary64's range), and `cr_verify.py` checks the primitive on
1,049,032 arguments (100,018 through the integer path alone) on x86_64 ELF, x86_64 PE and aarch64. The rows differ by one
ulp, inside class B's tolerance, so the harness marks a noted row as mapped whatever its class allows
(`parity_core.cyr`'s `_par_force_map`) and counts it as accepted under ADR 2: 827 `wavenumber_from_wavelength` and 827
`wavelength_from_wavenumber` rows.

**Reviewed and kept as Rust's** (reported to the user, the policy's borderline): `reduced_mass(m, m)` — Rust's
m·m/(m + m) misses the homonuclear m/2 by one ulp on 33 of 315 golden equal-mass rows (8.6 % of seeded masses);
`bohr_energy_level(Z, n)`'s ((−R·Z)·Z)/(n·n) misses −R exactly for Z = n on 3 of 12 golden rows (and −R/k² for Z/n = 1/k
on 30 %); the Rydberg factor 1/n₁² − 1/n₂², which cancels for close large n (computing it exactly would move
`rydberg_wavelength` by up to 3.9·10⁷ ulp at n ≈ 2³¹, 39 ulp at n = 100, 939 at n = 10⁴), and `bohr_transition_energy`'s
E(n₂) − E(n₁) likewise; the over- and underflow of products and quotients of inputs (ε·l·c and ε·l in the Beer–Lambert
functions, m₁·m₂ and m₁ + m₂ in `reduced_mass`, k/μ in `vibrational_frequency`, h·B in `rigid_rotor_energy`, 8π²·I in
`rotational_constant`, μ·r·r, the FFT's butterfly sums and re² + im² in its magnitudes, which Rust's `Complex::abs`
computes without hypot), where the result is representable; `find_peaks`' interpolated threshold min·(1 − t) + max·t
(it can round past DBL_MAX); and the isotope pattern's thresholds 0.001 (`>`: RN(0.001) > 1/1000, so only an intensity
of exactly RN(0.001) would be judged differently) and 10⁻¹² (a numerical pruning bound). Each is a chain of rounded
operations on inputs and physical constants, or a heuristic cutoff, which this ADR's policy leaves as Rust's sequence.

## Addition, 2026-11-10 (M4e-2): #20

**The user's decision (2026-11-09), verbatim: "Yes, prescale".** M4d's review had kept as Rust's the over- and underflow of
the squared coordinates in `distance_3d`, `bond_angle` and `dihedral_angle` (no single exact value to restore under the
inexact-formula policy) and reported it: Rust documents the angles as scale-invariant, and an exact power-of-two
prescaling would keep Rust's bits wherever nothing overflows. The user decided to fix it as this row.

**#20, the target.** The result is Rust's own operation sequence (every subtraction, product, sum, comparison, `max` and the
Err checks, in Rust's order) evaluated in IEEE binary64 round-to-nearest-even with an unbounded exponent range: every
+ − × ÷ √ rounds to 53 bits and nothing overflows or underflows. The exception is the final operation (`sqrt` for
`distance_3d`, `atan2` for the angles), which is correctly rounded to binary64 itself: its subnormal grid, or ±∞.
Comparisons and Err branches use the unbounded values. Non-finite inputs keep Rust's behaviour. Consequently every row
where Rust's chain stays in the normal range is bit-identical to Rust (the corpus shows it: only edge rows change); the
angles are exactly scale-invariant (for every k with all scaled coordinates exact, `bond_angle(2^k·a, 2^k·b, 2^k·c)` and
`dihedral_angle(…)` equal the unscaled result bit for bit, since every check compares quantities of the same degree and
atan2's two arguments have the same degree); and `distance_3d` scales by exactly 2^k wherever the result is normal. The two
typed tolerances keep their meaning, now at every scale: a bond below 10⁻¹² of its coordinates is coincident at 2⁹⁰⁰ m too.

**The implementation** (`src/potential.cyr`). Each function runs Rust's chain as before (the fast path) and a cheap test
that proves the chain stayed in the normal range; only where the test fails, and the coordinates are finite, it
re-evaluates with the unbounded-exponent helpers:
- `distance_3d`: s = (dx² + dy²) + dz² in [2⁻⁹⁰⁰, DBL_MAX], or all three differences 0. A sum of squares that large is the
  unbounded one: no square overflowed (that gives +∞), and a square that underflowed is below 2⁻¹⁰²² and is added to a
  partial sum it can't move (under a quarter ulp of it), or the sum would be below 2⁻⁹⁰⁰.
- `bond_angle`: |ba|², |bc|², max(|a|², |b|²), max(|c|², |b|²) and |ba × bc|² in [2⁻⁹⁰⁰, DBL_MAX] and ba · bc finite (the
  same argument for each sum of squares; a cross component that deviates is below 2⁻⁹⁶⁵, so its square is tiny too; and
  where |ba · bc| < 2⁻⁵¹⁰ against |ba × bc| ≥ 2⁻⁴⁵⁰ both atan2s round to RN(π/2), which lies 2⁻⁵⁴·¹⁵ inside its interval).
  The two Rust `is_positive` calls become these comparisons, so the fast path costs nothing measurable. An exactly
  linear or zero-angle geometry (CO₂ on an axis) has ba × bc exactly 0, so |ba × bc|² = 0 fails that test; there, when
  every component of ba and bc is 0 or at least 2⁻⁵¹¹ in magnitude, each cross and dot product is 0 or normal and every
  sum or difference of them is exact wherever it is subnormal, so y is +0 in both evaluations and the result is
  atan2(+0, ba · bc) (added in review: such a geometry took the slow path, 1.06 µs instead of 135 ns).
- `dihedral_angle`: every component of b₁, b₂, b₃ zero or of magnitude in [2⁻²²⁹, 2²⁵⁴) (from below, each component's
  square against 2⁻⁴⁵⁸, which is below it exactly for 0 < |x| < 2⁻²²⁹; from above, |b|² ≤ 2⁵⁰⁶). Then every nonzero product
  of two components is at least 2⁻⁴⁵⁸ and so a multiple of 2⁻⁵¹⁰, each cross component is 0 or in [2⁻⁵¹⁰, 2⁵⁰⁹], b₁ · n₂ is 0
  or at least 2⁻⁷⁹¹, and every product in the chain lies in [2⁻¹⁰²⁰, 2¹⁰²⁰]. (A terminal-value test does not suffice here:
  a deviating small cross component is multiplied by a large b₁ in b₁ · n₂.)

**The helpers** (`src/mathx.cyr`, general: the machinery of the library-wide audit below). A value is a pair (s, e) =
s·2^e with 1 ≤ |s| < 2 (or a signed zero or IEEE special, e = 0). Each operation is one IEEE operation on the significands
whose exact result lies far inside the normal range (× in [1, 4), ÷ in (½, 2), √ of s or 2s, + of s₁ and s₂·2^-d with
d ≤ 55 exact, a multiple of 2⁻¹⁰⁷), so IEEE's rounding is the 53-bit rounding of the exact value, and scaling by 2^e
commutes with it when no exponent bound intervenes; past d = 55 the smaller addend is below a quarter ulp of the larger,
which is then the sum. The final operations round once: a 53-bit value s·2^e is exact in binary64 when normal, ±∞ past
2¹⁰²³ (IEEE overflow is "the 53-bit rounding exceeds DBL_MAX"), and below 2⁻¹⁰²² its second rounding at the 2⁻¹⁰⁷⁴ grid
is wrong only when s·2^e is itself a midpoint of that grid (every such midpoint is a 53-bit number, so s never crosses
one), where the operation's exact error t decides; `_km_ux_round(s, e, t)` rounds (s + t)·2^e once with `_km_cr_finish`,
t from TwoProduct (×), the exact remainder (÷), TwoSum or a sticky ±2⁻⁶⁰ (+) and the residual w − r² (√). The final atan2
needs no scaling: `_km_cr_atan2`'s core works on (significand, exponent difference) pairs, `_km_atan2_tiny` rounding a
quotient below 2⁻⁵⁹ once at any exponent. Entry points: `_km_ux_from`, `_km_ux_neg`, `_km_ux_add`, `_km_ux_sub`,
`_km_ux_mul`, `_km_ux_div`, `_km_ux_sqrt`, `_km_ux_cmp`, `_km_ux_le`, `_km_ux_max`, `_km_ux_is_positive`, `_km_ux_round`,
the final `_km_ux_add_f64`, `_km_ux_sub_f64`, `_km_ux_mul_f64`, `_km_ux_div_f64`, `_km_ux_sqrt_f64` and `_km_ux_atan2_f64`
(`_fast`, `_acc`).

**Golden comparison** (as #16–#19): the corpus keeps the oracle's values; where the port's result differs, the note names it
— `[ADR 0002 #20: the port returns Rust's sequence evaluated with an unbounded exponent (√ rounded once) = <x> (f:…)]`, `…
(atan2 rounded once) = Ok(<x>) (f:…)]` or `… = Err(InvalidInput) (E:InvalidInput)]` — computed by the generator in exact
integer arithmetic (`scripts/oracle/src/exact.rs`'s `Ux`, 53-bit rounding of u128 mantissas at any exponent, and
`ux_atan2_f64`, a 400-bit fixed-point atan by Euler's series with Ziv's two-sided test), with a replica of Rust's f64 chain
checked against the oracle on every row deciding where the chains differ. `validate.py` evaluates every row with finite
coordinates independently with MPFR (`mpfr_ref.UX`: 53 bits at MPFR's widest exponent range, the final operation through
`mpfr_check_range` and `mpfr_subnormalize`) and proves: each noted row's result is MPFR's and differs from the oracle's;
each unnoted row's expected result (the corrected oracle's, or the oracle's) is MPFR's; and no row whose f64 chain stays in
the normal range is noted (1,000, 3,172 and 2,892 such rows: Rust's bits). It also checks the scale invariance through MPFR.
The harness maps 854 `distance_3d`, 1,205 `bond_angle` and 1,546 `dihedral_angle` rows (91 of them the port's Err, where
Rust's threshold (10⁻²⁰·|b₁|²)·|b₂|² rounded 10⁻²⁰·|b₁|² at the subnormal grid before |b₂|² scaled it back up) and prints one
summary line. `cr_verify.py` checks the final operations against MPFR (uxsqrt, uxmul, uxdiv, uxadd, uxatan2: 1,040,039
arguments) and the three functions, public and slow path forced, on 114,000 edge geometries, on x86_64 ELF, x86_64 PE and
aarch64; `tests/potential.tcyr` checks the scale invariance bit for bit on 600 random geometries × 75 scales and the slow
path against the fast one on 5,000 in-range geometries.

**The library-wide audit (the user's decision, the same day).** After M4, every ported module is checked for intermediate
overflow or underflow of products and quotients of inputs, and each such edge fixed the same way: Rust's sequence with an
unbounded exponent, the final operation rounded once, in-range rows keeping Rust's bits. M4d's and M4e's "kept as Rust's"
lists (`lennard_jones`' (σ/r)⁶ and (σ/r)¹², Coulomb's k·q₁·q₂ and r·r, `morse`'s subnormal D_e·x, Beer–Lambert's ε·l·c,
`reduced_mass`'s m₁·m₂ and m₁ + m₂, k/μ, h·B, 8π²·I, μ·r·r, the FFT's sums and re² + im², and the like) are its first
candidates; the helpers above are its machinery. The audit is not part of #20.

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
