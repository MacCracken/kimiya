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
| 2 *(realised 2026-11-12, M5b)* | `sensitivity::{equilibrium_constant_gradient, arrhenius_gradient, nernst_gradient, gibbs_temperature_sensitivity}` unvalidated; `compute_gradient` → `Result` never `Err` | the four return `Result`, validated as `equilibrium_constant`, `arrhenius_rate`, `nernst_potential` and `gibbs_free_energy` (codes and order), Err(ComputationError) for a non-finite entry, computed from closed-form derivatives (each tape's sweep unrolled); `compute_gradient` keeps a `Result` for the Cyrius-only failures of a caller's tape (the M5b addition) | finite oracle rows class D (0 ulp); oracle non-finite, or an input the port refuses ↔ the port's noted Err (227 rows, mapped, counted) |
| 3 | `nuclear::semi_empirical_binding_energy` → NaN for z > a or a = 0 | `Result` | finite rows class B; oracle NaN ↔ `Err` |
| 4 | `phase::antoine_pressure_mmhg`, `antoine_pressure_pa` → NaN at the pole, silent extrapolation | `Result`; `Err` outside the coefficient set's `t_min_c..t_max_c`. No extrapolating variant | in-range rows class C; oracle NaN or out-of-range input ↔ `Err` |
| 6 | `FunctionalGroup::hydrogen_bonding`; `Geometry::TShapedThreeLonePairs` | `functional_group_is_hydrogen_bond_donor` and a new `functional_group_is_hydrogen_bond_acceptor`; `KM_GEOMETRY_T_SHAPED` | donor and geometry rows class A under the new names; acceptor checked by reference tests only |
| 7 *(realised 2026-11-10, M4f)* | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational oxidation state (numerator / positive denominator, reduced; a value outside i32 is Rust's Err); balancing takes an integer matrix (flat row-major, Rust's \|entry\| ≤ 2⁵³ kept) | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer average oxidation state" `Err` ↔ `Ok(p/q)`, noted with the port's whole result (9,090 rows, mapped, counted) and verified by reference tests (Fe₃O₄ +8/3, C₃H₈ −8/3); balancing class A after converting the oracle's integer-valued f64 inputs as Rust's test does; the 78 rows with no integer counterpart (non-integer, non-finite, beyond i64, ragged) keep the oracle's Err(InvalidInput), noted and counted |
| 12 *(2026-11-09)* | `inorganic::cfse_octahedral` computes `t2g · (−0.4) + eg · 0.6`, which misses the exact (3·e_g − 2·t₂g)/5 on 13 of its 22 (d, spin) cases (high-spin d⁵ −2⁻⁵² and d¹⁰ −2⁻⁵¹ instead of 0, d³ −1.2000000000000002) | the exact value: the integer numerator 3·e_g − 2·t₂g converted to f64, then one correctly rounded division by 5 | class A. Where the oracle's value equals the exact one (9 cases), a plain row; elsewhere the golden row keeps the oracle's value, its note names the exact one, and the harness maps the two (q = round(5 × oracle) must be within 10⁻⁹ of it, and the port must return q / 5), reporting each such row by name as accepted under ADR 2 |
| 13 *(2026-11-09, M4c; completed in M4d)* | `nuclear::atoms_after_half_lives` computes `initial * 2.0.powf(-n)`: where 2⁻ⁿ is not a normal double (n > 1022 or n ≤ −1024) it overflows to ∞, underflows to 0 (n = 1075 is a tie, to 0) or is a subnormal of 52 bits or fewer, so the result is ±∞, 0, NaN or up to 2⁵² ulp off where N₀/2ⁿ is representable (1000 after 1075 half-lives: 0, not 2.47 × 10⁻³²¹; 2¹⁰⁰ after 1,100: 0, not 2⁻¹⁰⁰⁰) | there, N₀/2ⁿ **correctly rounded** (rounded once, subnormal results at the 2⁻¹⁰⁷⁴ grid, ties to even): since M4d `src/mathx.cyr`'s `_km_cr_exp2_scaled(N₀, −n)` for every N₀ and n, 2⁻ᶠ from exp's machinery past 53 bits with a rounding test and a triple-double fallback (M4c's version was faithful for a fractional n with N₀ not a power of two). ±0, ±∞, NaN N₀ give N₀ (scaleB). Rust's form stays elsewhere | class C; equal rows plain; the others keep the oracle's value with the port's in the note (`(f:…)`), which the harness requires bit for bit and reports as accepted under ADR 2, counted (1,810 rows) |
| 14 *(2026-11-09, M4c)* | `nuclear::is_secular_equilibrium` tests `λ_p < 0.01·λ_d`: 0.01 is not 1/100 and the product rounds or underflows, so (0.03, 3) is false and (0, 5 × 10⁻³²⁴) false | the exact λ_p < λ_d/100: 100·λ_p against λ_d by an exact double-double (`_km_two_prod` on the significands, the exponents first), Rust's NaN (false), ±∞ and ±0 semantics kept | class A; equal rows plain; the others noted `(i:…)`, mapped and reported by name (208 rows) |
| 15 *(2026-11-09, M4b's inorganic, decided in M4c)* | `inorganic::delta_tet_from_oct` computes `delta_oct * 4.0 / 9.0`, correctly rounded except where Δ·4 overflows (\|Δ\| ≥ 2¹⁰²², above DBL_MAX/4): ±∞ where 4Δ/9 is finite | there, (Δ/9)·4: one rounded division and an exact ×4 that can't overflow, which is again the correctly rounded 4Δ/9 (4·DBL_MAX/9 ≈ 7.99 × 10³⁰⁷). Rust's form elsewhere | class B; noted `(f:…)`, mapped and reported by name (59 rows) |
| 16 *(2026-11-09, M4d)* | `potential::lj_equilibrium_distance` computes `sigma * 2.0.powf(1.0 / 6.0)`: the exponent 1/6 is rounded (a simple exact constant misrepresented in binary) and the product rounds again, so r_min = σ·2^(1/6) is one ulp off on about a fifth of all σ (4,117 of 20,000 seeded). (2^RN(1/6) rounds to RN(2^(1/6)), MPFR checks; the misses are the second rounding) | the correctly rounded σ·2^(1/6), rounded once: σ = s·2^e, s·(H + M + L) by `_km_td_mul_d` (2^(1/6) as a triple-double from `mathx_tables.py`'s exact integer root), then `_km_cr_finish` at the final scale; no rounding test is needed, since 2^(1/6)'s continued fraction proves every σ·2^(1/6) at least 2^-55.1 ulp from a rounding boundary, far beyond the product's 2^-102.6. NaN, ±0, ±∞ give σ | class C (Rust's `powf`); equal rows plain; the others noted `(f:…)`, mapped and counted (1,041 rows) |
| 17 *(2026-11-09, M4d)* | `potential::lennard_jones` computes `4.0 * epsilon * (sr12 - sr6)`, i.e. (4·ε)·d: for ε ≥ 2¹⁰²² the exact 4·ε overflows, so the result is ±∞ (or ∞·0 = NaN at r = σ) where 4·ε·d is representable | there, (ε·d)·4: one rounded product and an exact ×4, which is RN(4·ε·d), the value Rust's (4·ε)·d has wherever 4·ε is finite (4·ε is exact). Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (180 rows). *Folded into #32 by the overflow audit's B3 (2026-11-11): (ε·d)·4 is #32's unbounded-exponent value wherever Rust's d is normal, so all 180 rows keep their value, now noted under #32* |
| 18 *(2026-11-09, M4d)* | `potential::morse_width_parameter` computes `omega_e * (reduced_mass / (2.0 * d_e)).sqrt()`: for D_e ≥ 2¹⁰²³ the exact 2·D_e overflows, μ/∞ = 0 and a = 0 where ω·√(μ/(2·D_e)) is representable | there, (μ·½)/D_e: RN(μ/(2·D_e)) (μ·½ is exact, or μ < 2⁻¹⁰²¹ and the quotient rounds to 0 either way), the value Rust's μ/(2·D_e) has wherever 2·D_e is finite; then Rust's √ and ×ω. Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (181 rows). *Folded into #32 by the overflow audit's B3 (2026-11-11), which rounds the whole ω·√(μ/(2·D_e)) once, everywhere: 128 of the 181 rows keep their value, 42 move by one ulp and 3 by two (√ and ×ω rounding a normal quotient), 8 by up to 1.3·10¹² ulp (a subnormal μ/(2·D_e), degraded before the √); all noted under #32* |
| 19 *(2026-11-09, M4e)* | `spectroscopy::wavenumber_from_wavelength` and `wavelength_from_wavenumber` compute `0.01 / x` for the exact 1/(100·x) (the cm ↔ m factor 1/100): 0.01 is 1/100 + 2.08·10⁻¹⁹ (relative 2⁻⁵⁵·⁴), and the quotient rounds again, so the result is one ulp from RN(1/(100·x)) on about 14 % of all doubles (2,748 of 20,000 seeded) | RN(1/(100·x)), rounded once: `src/mathx.cyr`'s `_km_cr_hundredth_div` (a double-double residual with 1/100 = H + L and a rounding test, an exact 64-bit long division as its fallback and past 2^±900); the special values Rust's | class B; noted `(f:…)`, mapped and counted (827 rows each) |
| 20 *(2026-11-09, M4e-2)* | `potential::distance_3d`, `bond_angle` and `dihedral_angle` square raw coordinates: above about 10¹⁵⁴ (10⁷⁷ for the dihedral's fourth powers) or below their reciprocals a square, cross product or dot product overflows or underflows, and the result is an Err, 0, ±∞ or a degraded value where the true distance or angle is representable (`distance_3d` of (DBL_MAX, 0, 0) and the origin is +∞, the gauche dihedral at 2⁷⁰⁰·Å an Err, the angle of (1, 2⁻¹⁰⁶⁰) 0), so the angles, which Rust documents as scale-invariant, are not | Rust's own operation sequence (every subtraction, product, sum, comparison, `max` and the Err checks, in Rust's order) evaluated in IEEE binary64 round-to-nearest-even with an **unbounded exponent range**: every + − × ÷ √ rounds to 53 bits and nothing overflows or underflows. The exception is the **final** operation (`sqrt` for `distance_3d`, `atan2` for the angles), which is correctly rounded to binary64 itself: its subnormal grid, or ±∞. Comparisons and Err branches use the unbounded values; non-finite inputs keep Rust's behaviour. So every row where Rust's chain stays in the normal range is bit-identical to Rust; the angles are exactly scale-invariant (for every k with all scaled coordinates exact, the same bits at 2^k times the coordinates: every check is homogeneous, and `atan2`'s two arguments have the same degree), and `distance_3d` scales by exactly 2^k wherever the result is normal. Rust's chain with a cheap range test, and `src/mathx.cyr`'s `_km_ux_*` helpers where it fails | class B (`distance_3d`) and C (the angles); noted `(f:…)` or `(E:InvalidInput)`, mapped and counted (854, 1,205 and 1,546 rows; 91 of the dihedral's are the port's Err where Rust's rounded threshold says Ok) |
| 21 *(2026-11-10, M4e-3)* | `spectroscopy::bohr_energy_level` computes ((−R·Z)·Z)/(n·n), `bohr_transition_energy` the difference of two such levels, `rydberg_wavelength` and `SpectralSeries::wavelength` 1/(((R·Z)·Z)·(1/(n₁·n₁) − 1/(n₂·n₂))): three to eight roundings of exact rationals, so −R at Z = n is one ulp off for 86 of the 255 Z, and the levels' difference and 1/n₁² − 1/n₂² cancel for close large n (39 ulp at n = 100, 939 at 10⁴, up to 9.0·10⁷ for the wavelength and 1.8·10⁹ for the transition energy on the corpus, near n = 2³²) | the exact rational with the constant as its double, **rounded once**: −R·Z²/n², R·Z²·(n₂² − n₁²)/(n₁²·n₂²) and n₁²·n₂²/(R·Z²·(n₂² − n₁²)) (R = `KM_RYDBERG_EV`, `KM_RYDBERG`, or Rust's const-evaluated R_H, Z = 1): a double-double fast path with a proven bound and a rounding test (exact operands for n₂ < 2¹³), and the exact slow path `_km_cr_ratio` on multi-word integers, which rounds the 13 exact ties (Z = 255, 2^k·(58,055, 524,544)) to even. Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (3,823 levels, 1,351 transition energies, 1,376 wavelengths, 67 series lines) |
| 22 *(2026-11-10, M4e-3)* | `spectroscopy::reduced_mass` computes (m₁·m₂)/(m₁ + m₂): three roundings, so a homonuclear m/2 is one ulp off on 8.6 % of masses, and m₁·m₂ or m₁ + m₂ overflows or underflows at the range ends (NaN at (DBL_MAX, DBL_MAX), 0 or ∞ where μ is representable) | the correctly rounded m₁·m₂/(m₁ + m₂) for every pair of positive finite doubles, subnormal results and their ties included: with m₁ ≥ m₂ and the exponent gap g, the smaller mass for g ≥ 55 (proven), otherwise s₁s₂/(s₁ + s₂·2⁻ᵍ) scaled by 2^e₂ (no intermediate leaves the range), a double-double fast path and the exact `_km_cr_ratio` | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (2,837 rows) |
| 23 *(2026-11-10, M4e-4)* | `spectroscopy::isotope_abundances`' table rounds the AME2020 masses to five decimals (¹H 1.00783 for 1.007825031898 u, ²H 2.01410 for 2.014101777844; each within 5·10⁻⁶ u) and omits ³⁶S (0.0001), so sulfur's abundances sum to 0.9999; `isotope_pattern` convolves that table | the cited table: the AME2020 masses at full published precision and the TICE 2009 representative abundances (Rust's, unchanged), ³⁶S added; each value the double nearest its published decimal, from the records in `scripts/oracle/isotopes.tsv`. `isotope_pattern`'s algorithm, Err codes and order are Rust's; only its data change | class A (`isotope_abundances`) and B (`isotope_pattern`); the rows whose result changes noted with the port's whole result `(v:n f:… …)`, mapped and counted (7 and 347 rows) |
| 24 *(2026-11-10, M4f-2)* | `stoichiometry::percent_yield` computes (a/t)·100, `theoretical_yield_moles` (m/c)·p and `percent_composition` ((mass·count)/mw)·100: two to four roundings of exact rationals, so about a quarter of results are 1–2 ulp off, and at the range ends the result is 0 or ∞ where the value is representable (`percent_yield(5e-324, 10)` is 0, not 5e-323; the yield's m/c over- or underflows though m·p/c is finite) | each exact rational **rounded once**: RN(100·a/t), RN(m·p/c) with m and c at #25's index, RN(mass·count·100/mw) with mw `molecular_weight`'s own double — `_km_cr_mul_div`: a double-double fast path on the significands with `_km_dd_div_cr`'s rounding test and an exact resolution of the midpoint where it fails, the exact `_km_cr_ratio` for subnormal operands and results that may leave the normal range. Zero operands keep Rust's ±0; Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)` (the yields) or with the whole result (composition), mapped and counted (1,126 `percent_yield`, 240 `theoretical_yield_moles`, 598 `percent_composition`) |
| 25 *(2026-11-10, M4f-2)* | `stoichiometry::limiting_reagent` compares the rounded ratios m/c, so it can pick the wrong reactant when two ratios are within an ulp, or when both overflow to ∞ or underflow (`[DBL_MAX, DBL_MAX]` over `[5e-324, 1e-300]` gives 0; the exact answer is 1) | the exact comparison: i replaces the current minimum i₀ when mᵢ·c_{i₀} < m_{i₀}·cᵢ exactly (the first index on an exact tie, Rust's strict `<`); the rounded ratios still decide wherever they differ (RN is monotone), and equal ones go to `_km_ratio_lt`, an exact comparison of the two products (exact TwoProducts of the significands, the exponents first). Rust's validation and its order are kept | class A; equal rows plain; the others noted `(i:…)`, mapped and counted (920 rows) |
| 26 *(2026-11-10, M4f-2)* | `stoichiometry::empirical_formula` tests \|r·k − round(r·k)\| ≤ 0.1 on rounded mole ratios: inputs at the ±1/10 edge flip between Ok and Err (14 of M4f's 3,383 rows; C 12.011 % with O 17.5989 % is an Err, but exactly CO), and its "non-positive min" and "ratio not finite" Errs are underflow and overflow artefacts | exact evaluation: the moles pctᵢ/massᵢ and their minimum exact, the ratios rᵢ = pctᵢ·mass_j/(massᵢ·pct_j) exact, the multiplier the smallest k in 1 ..= 6 with \|rᵢ·k − round(rᵢ·k)\| ≤ 1/10 exactly for every i, the counts round(rᵢ·k) exactly (past u32::MAX Rust's Err); a ratio of 2³² − ½ or more is the ComputationError (proven), the others a double-double with a tolerance test that is proven or defers to an exact integer test. Rust's validation and its order are kept | class A; equal rows plain; the others noted with the whole result `(v:n i:z i:c …)` or `(E:ComputationError)`, mapped and counted (270 rows) |
| 27 *(2026-11-10, the overflow audit's B1)* | `thermo::heat_transfer` computes (m·c)·ΔT, `specific_heat_from_calorimetry` q/(m·ΔT), `phase::boiling_point_elevation` and `freezing_point_depression` (K·m)·i, `osmotic_pressure` (((i·M)·R)·1000)·T: two to four roundings of one product or quotient of inputs, so about a quarter of results are 1–2 ulp off, and a partial product that over- or underflows gives 0, ±∞ or NaN although the value is representable (`specific_heat_from_calorimetry(4184, 1000, DBL_MAX)` 0, not 2.33·10⁻³⁰⁸; `osmotic_pressure(5e-324, 0.1, 298.15)` 0, not 1.22·10⁻³¹⁸; `heat_transfer(1e300, 1e300, 0)` ∞·0 = NaN, not 0) | each exact product or quotient **rounded once** — RN(m·c·ΔT), RN(q/(m·ΔT)), RN(K·m·i), RN(i·M·R·1000·T) with R `KM_GAS_CONSTANT`'s double and 1000 exact — by `src/mathx.cyr`'s correctly rounded monomial `_km_cr_mono` (M4f-2's `_km_cr_mul_div` generalised: a double-double fast path with `_km_dd_div_cr`'s rounding test and an exact midpoint resolution, the exact `_km_cr_ratio` past the normal range): the sign the factors' xor, a zero factor that signed zero, ties to even, subnormal results on the 2⁻¹⁰⁷⁴ grid, ±∞ only past DBL_MAX. A non-finite input of the three unvalidated functions keeps Rust's expression; Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (151 `heat_transfer`, 171 `specific_heat_from_calorimetry`, 82 + 82 colligative, 196 `osmotic_pressure`: 682) |
| 28 *(2026-11-10, the overflow audit's B1)* | `kinetics::zero_order_half_life` computes [A]₀/(2k), `steady_state_intermediate` (k₁·[A])/k₂, `pre_equilibrium_rate` (k₂·(k₁/k₋₁))·[A], `second_order_concentration` 1/(1/[A]₀ + k·t) and `michaelis_menten` (Vmax·[S])/(Km + [S]): an intermediate over- or underflows where the result is representable (`zero_order_half_life(DBL_MAX, DBL_MAX)` 0, not ½; `second_order_concentration(1, DBL_MAX, 10)` 0, not 5.56·10⁻³¹⁰; `michaelis_menten(100, DBL_MAX, DBL_MAX)` ∞/∞ = NaN, not 50; `pre_equilibrium_rate(1e300, 1, 1e300, 0)` ∞·0 = NaN, not 0), and the monomials' two or three roundings are 1–2 ulp off on a quarter of inputs | the three monomials **rounded once** (#27's helper: RN([A]₀/(2k)) — Rust's own division wherever 2k is finite, where it already is that value —, RN(k₁·[A]/k₂), RN(k₂·k₁·[A]/k₋₁)); the two with a sum inside under **#20's semantics**: Rust's operation sequence with an unbounded exponent (`_km_ux_*`), the final division rounded once into binary64, behind a range test (`_km_second_order_fast`, `_km_michaelis_menten_fast`) that proves Rust's chain stayed normal, so those inputs keep Rust's bits. Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (89 `second_order_concentration`, 90 `zero_order_half_life`, 118 `steady_state_intermediate`, 179 `pre_equilibrium_rate`, 160 `michaelis_menten`: 636) |
| 29 *(2026-11-10, the overflow audit's B2)* | `solution::dilution` computes (M₁·V₁)/V₂, `solution_conductivity` (Λ·c)/1000, `electrochemistry::faraday_mass_deposited` (Q·M)/(n·F), `faraday_charge_required` ((m·n)·F)/M, `faraday_mass_ratio` (m_A·(M_B/n_B))/(M_A/n_A) and `nernst_potential` E° − ((R·T)/(n·F))·ln Q: the five quotients round two to four times, so about a quarter of results are 1–2 ulp off (`dilution(1.5, 0.2, 0.6)` 0.5000000000000001), and a partial product, an equivalent weight M/n or R·T that over- or underflows gives ±∞, 0, NaN or a degraded value although the value is representable (`faraday_mass_deposited(DBL_MAX, 63.546, 2)` ∞, not 5.92·10³⁰⁴; `solution_conductivity(126.4, DBL_MAX)` ∞, not 2.27·10³⁰⁷; `faraday_mass_ratio(1, 5e-324, 3, 1e-300, 1)` ∞, not 6.07·10²³; `nernst_potential(0.342, 2, DBL_MAX, 0.1)` ∞, not 1.78·10³⁰⁴, and at Q = 1 NaN, not E°) | the five monomials **rounded once** — RN(M₁·V₁/V₂), RN(Λ·c/1000), RN(Q·M/(n·F)), RN(m·n·F/M), RN(m_A·M_B·n_A/(n_B·M_A)), F as `KM_FARADAY`'s double, 1000 and the electron counts exact — by #27's `_km_cr_mono` (two new straight-line shapes, `_km_mul2_div2_fast` and `_km_mul3_div2_fast`); `nernst_potential` under **#20's semantics** with the correctly rounded ln Q (`_km_cr_ln`, ADR 0006) as its primitive: Rust's operation sequence with an unbounded exponent (`_km_ux_*`), the final subtraction rounded once into binary64 (`_km_nernst_ux`), behind a range test (`_km_nernst_fast`) that proves Rust's chain stayed normal, so those inputs keep Rust's bits (every `nernst_potential_25c` input among them; corrected in B4: the test's bound is `_km_ux_normal`'s 2⁻¹⁰²¹, since B2's admitted a product c·ln Q that IEEE rounds up to exactly 2⁻¹⁰²² from below, where #20's 53-bit value is 2⁻¹⁰²² − 2⁻¹⁰⁷⁵). The unvalidated `solution_conductivity` keeps Rust's expression for a non-finite input; Rust's checks, Err codes and order are kept | class B (the monomials) and C (`nernst_potential`); equal rows plain; the others noted `(f:…)`, mapped and counted (94 `dilution`, 148 `solution_conductivity`, 904 `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 136 `nernst_potential` (119 in B2, 17 more with B4's correction): 2,773) |
| 30 *(2026-11-10, the overflow audit's B2)* | `spectroscopy::absorbance` computes (ε·l)·c, `concentration_from_absorbance` A/(ε·l), `rigid_rotor_energy` ((h·B)·J)·(J + 1), `rotational_constant` h/(((8·π)·π)·I), `diatomic_moment_of_inertia` (μ·r)·r and `vibrational_frequency` (k/μ).sqrt(): two to four roundings of one product or quotient of inputs and constants (one quotient under a √), so a fraction of results are 1–2 ulp off (HCl's ω = √(516/1.6266·10⁻²⁷) is 563228298138879.9, not …880.0; `rotational_constant(1.45e-46)` 57875971594.410576, not …41057), and a partial product or quotient that over- or underflows gives ±∞, 0, NaN or a degraded value although the value is representable (`vibrational_frequency(500, 2⁻¹⁰²²)` ∞, not 1.50·10¹⁵⁵, and `(10⁻³⁰⁰, 10³⁰)` 0, not 10⁻¹⁶⁵; `rigid_rotor_energy(2³² − 1, 10⁻³⁰⁰)` 0, not 1.22·10⁻³¹⁴; `rotational_constant(10⁻³²⁰)` 8.392090·10²⁸⁴, not 8.392109·10²⁸⁴; `absorbance(1e300, 1e300, 0)` NaN, not 0; `concentration_from_absorbance(1e300, 1e200, 1e200)` 0, not 10⁻¹⁰⁰) | the five monomials **rounded once** — RN(ε·l·c), RN(A/(ε·l)), RN(h·B·J·(J + 1)), RN(h/(8·π·π·I)) (h `KM_PLANCK`'s double, π its double each time Rust uses it, 8, J and J + 1 exact), RN(μ·r·r) — by #27's helper (two new straight-line shapes, `_km_mul4_fast` and `_km_div3_fast`); `vibrational_frequency` the **correctly rounded √ of the exact k/μ**, `src/mathx.cyr`'s `_km_cr_sqrt_div`: a double-double fast path (proven within 5.7u², u = 2⁻⁵³) with `_km_dd_div_cr`'s rounding test, and an exact integer square root on `_km_bn_*` (√(k/μ) is never a rounding boundary of binary64, proven), ±∞ only past DBL_MAX. A non-finite unvalidated input (absorbance's, the moment of inertia's, the rotor's B, concentration_from_absorbance's A) keeps Rust's expression; Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (195 `absorbance`, 219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant`, 118 `diatomic_moment_of_inertia`: 2,106) |
| 31 *(2026-11-11, the overflow audit's B3)* | `gas::ideal_gas_pressure` and `ideal_gas_volume` compute ((n·R)·T)/V and ((n·R)·T)/P, `gas_density` (P·M)/(R·T), `van_der_waals_pressure` (((n·R)·T)/(V − n·b)) − ((a·(n/V))·(n/V)) after testing V ≤ n·b, `peng_robinson_pressure` RT/(V − b) − a·α/(V² + 2bV − b²) on `pr_coefficients`' a·α and b after testing V ≤ b, `redlich_kwong_pressure` RT/(V − b) − a/((√T·V)·(V + b)) with a = (((0.42748·R)·R)·Tc^2.5)/Pc: the three monomials round three times, so about a quarter of results are 1–2 ulp off, and in all six an intermediate that over- or underflows gives ±∞, 0, NaN, the wrong sign or a degraded value although the value is representable (`ideal_gas_volume(DBL_MAX, 300, 101325)` ∞, not 4.43·10³⁰⁶; `gas_density(101325, 0.029, DBL_MAX)` 0, not 1.97·10⁻³⁰⁶; `peng_robinson_pressure` at Tc = 2⁻¹⁰²² NaN, not 99,773.55 Pa, and at V = 2·10¹⁵⁴ +1.25·10⁻¹⁵¹, not −1.49·10⁻⁷; `redlich_kwong_pressure(10¹³⁰, 10³⁰⁰, 300, 1)` −∞, not −1.71·10²⁵), and the co-volume tests decide on a rounded co-volume (a subnormal n·b or b, or one whose subnormal partial product /Pc scales back up) | the three monomials **rounded once** — RN(n·R·T/V), RN(n·R·T/P), RN(P·M/(R·T)), R as `KM_GAS_CONSTANT`'s double — by #27's `_km_cr_mono` (B1's and B2's shapes `_km_mul3_div_fast`, `_km_mul2_div2_fast`); the three equations of state under **#20's semantics**: Rust's operation sequence (`pr_coefficients`' κ and α = u·u included) with an unbounded exponent (`_km_ux_*`), the co-volume tests on the unbounded co-volume (so the port refuses where Rust's rounded one accepts, and the reverse), the final subtraction rounded once, behind range tests (`_km_vdw_fast`, `_km_pr_pressure_fast`, `_km_rk_pressure_fast`, on the new `_km_ux_normal`) that prove Rust's chain stayed normal, so those inputs keep Rust's bits; Redlich–Kwong's Tc^2.5 with an unbounded exponent is the exact rescale `_km_cr_pow(s, 2.5)·2^(5j)` for Tc = s·2^(2j), s in [1, 4) (`_km_rk_tc25_ux`, proven equal to RN₅₃(Tc^2.5)). The Peng–Robinson Z, roots and fugacity keep Rust's chain (the audit's later batches). Rust's checks, Err codes and order are kept | class B (C for Redlich–Kwong); equal rows plain; the others noted `(f:…)` or, where the port's co-volume test refuses, `(E:InvalidInput)`, mapped and counted (146 `ideal_gas_pressure`, 148 `ideal_gas_volume`, 222 `gas_density`, 206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234 `redlich_kwong_pressure`: 1,157, of which 98 the port's Err and 80 its value where Rust's is the Err) |
| 32 *(2026-11-11, the overflow audit's B3)* | `potential::lennard_jones` computes (4·ε)·((σ/r)¹² − (σ/r)⁶) (#17's (ε·d)·4 past 4·ε's overflow), `coulomb` ((k·q₁)·q₂)/r, `coulomb_force` ((k·q₁)·q₂)/(r·r), `morse_width_parameter` ω·(μ/(2·D_e)).sqrt() (#18's (μ·½)/D_e past 2·D_e's overflow), `inorganic::born_lande_lattice_energy` −k·N_A·M·z⁺·z⁻·e·e/r₀·(1 − 1/n) left to right, `nuclear::binding_energy_per_nucleon` ((Z·m_H + N·m_n − M)·931.494…)/A: the monomials round three or four times (1–2 ulp off on a fraction of inputs; the Morse width a √ of a rounded quotient, times ω), and an intermediate that over- or underflows gives ±∞, 0, NaN, the wrong zero or a degraded value although the value is representable (`lennard_jones(10⁻³⁰⁰, 1, 1.72·10⁻²⁶)` +∞, not 5.83·10⁹; at σ/r > 2^170.7 ∞ − ∞ = NaN where the value is +∞ or finite; `coulomb_force(5·10⁻³²⁴, e, 10⁻¹⁰)` 0, not 7.11·10⁻³¹³; `morse_width_parameter(5.634·10¹⁴, DBL_MAX, 7.399·10⁻¹⁹)` ∞, not 6.21·10¹⁷⁷; `born_lande_lattice_energy(DBL_MAX, 1, 1, DBL_MAX, DBL_MAX)` −∞, not −1.39·10⁻⁴; `binding_energy_per_nucleon(1000, 1000, DBL_MAX)` −∞, not −8.37·10³⁰⁷) | `coulomb` and `coulomb_force` **rounded once** — RN(k·q₁·q₂/r), RN(k·q₁·q₂/(r·r)), k as `KM_COULOMB_CONSTANT`'s double — by #27's helper (`_km_mul3_div_fast`, `_km_mul3_div2_fast`), a non-finite charge keeping Rust's expression; `morse_width_parameter` the **correctly rounded ω·√(μ/(2·D_e))** of the exact value, `src/mathx.cyr`'s new `_km_cr_mul_sqrt_ratio` (RN(w·√((x/y)·2^e)): #30's proven double-double root times ω with a rounding test, and an exact integer root on `_km_bn_*`, which rounds the midpoints this value can be, 3ω·2^k, to even); `lennard_jones`, `born_lande_lattice_energy` (whose (1 − 1/n) stays Rust's, M4b's review) and `binding_energy_per_nucleon` (whose mass defect is a sum of inputs) under **#20's semantics**: Rust's operation sequence with an unbounded exponent, the final product or division rounded once, behind range tests (`_km_lj_fast`, `_km_born_lande_fast`, Δm·931.494… finite) that prove Rust's chain stayed normal, so those inputs keep Rust's bits. #17 and #18 are folded in: #17's (ε·d)·4 is the unbounded chain's value wherever Rust's d is normal (all 180 of its rows keep their value); #18's ω·√(RN(μ/(2·D_e))) gives way to the correctly rounded value (128 of its 181 rows keep theirs, 53 move: the B3 addition). Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (582 `lennard_jones`, 244 of them with ε ≥ 2¹⁰²²; 1,860 `morse_width_parameter`, 340 with D_e ≥ 2¹⁰²³; 730 `coulomb`; 679 `coulomb_force`; 175 `born_lande_lattice_energy`; 156 `binding_energy_per_nucleon`: 4,182) |
| 33 *(2026-11-11, the overflow audit's B4)* | `thermo::hess_law` sums its slice from −0.0, `enthalpy_from_formation` computes Σ products − Σ reactants, thermochem's `reaction_enthalpy`, `reaction_gibbs_energy` and `reaction_entropy` Σ coeff·prop (from +0.0) over the products minus the same over the reactants (#1: × 1000 for the two in J/mol), `nuclear::q_value` (Σ m_r − Σ m_p)·931.494…, `solution::ionic_strength` ½·Σ (c·z)·z, `inorganic::born_haber_lattice_energy` ΔH_f − ΔH_sub − IE − ΔH_diss/2 − EA left to right, `reaction::gibbs_free_energy` ΔH − T·ΔS, `solution::kohlrausch_molar_conductivity` Λ° − K·√c, `limiting_molar_conductivity` ν₊·λ₊ + ν₋·λ₋ and `kinetics::zero_order_concentration` [A]₀ − k·t (0 if negative): sums and differences in which a partial sum or a product term passes DBL_MAX and a later term brings it back, so the result is ±∞ or ∞ − ∞ = NaN although the value is finite (`hess_law([DBL_MAX, DBL_MAX, −DBL_MAX])` +∞, not DBL_MAX; `reaction_enthalpy([Al₂O₃ ×10³⁰⁸] → [Al₂O₃ ×10³⁰⁸])` NaN, not 0; `q_value([M, M], [M, 0.999M])` NaN, not 1.67·10³⁰⁸; `ionic_strength([(0.4·DBL_MAX, 2)])` +∞, not 1.44·10³⁰⁸; `gibbs_free_energy(DBL_MAX, 2, 0.75·DBL_MAX)` −∞, not −8.99·10³⁰⁷), or a subnormal product is rounded before the sum | **#20's semantics** for all twelve: Rust's operation sequence, its summation order included, evaluated with an unbounded exponent (`_km_ux_*`; Rust's `.sum()` from −0.0 and thermochem's `total += coeff * prop` from +0.0, `_km_sum_ux`, `_km_sum_property_ux`), the final operation rounded once into binary64 (the last addition, the difference, ½·Σ, × 931.494… or × 1000), comparisons on the unbounded values (`zero_order_concentration`'s clamp), behind a fast-path test that proves Rust's chain stayed normal, so those inputs keep Rust's bits. For a pure sum the test is "result finite": a sum of two doubles is exact below 2⁻¹⁰²¹ (both are multiples of 2⁻¹⁰⁷⁴) and the 53-bit rounding above, so only overflow can leave the range, and ±∞ or NaN persists to the result (proven in the addition below); a product inside takes `_km_ux_normal` or an exact zero. A non-finite input keeps Rust's expression; Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)` (after #1's note for the two J/mol sums), mapped and counted (28 `hess_law`, 27 `enthalpy_from_formation`, 115 `reaction_enthalpy`, 113 `reaction_gibbs_energy`, 66 `reaction_entropy`, 58 `q_value`, 13 `ionic_strength`, 74 `born_haber_lattice_energy`, 57 `gibbs_free_energy`, 54 `kohlrausch_molar_conductivity`, 31 `limiting_molar_conductivity`, 61 `zero_order_concentration`: 697) |
| 34 *(2026-11-11, the overflow audit's B4)* | `thermochem::Shomate::cp` computes (((A + B·t) + (C·t)·t) + ((D·t)·t)·t) + E/(t·t) and `enthalpy_antiderivative` 1000·(t·(A + t·(B/2 + t·(C/3 + t·(D/4)))) − E/t), t = T/1000, for any coefficients and T: as T → 0, t·t underflows and E/(t·t) is ±∞, or 0/0 = NaN for E = 0, where Cp is finite (`Shomate{20.786, 0, 0, 0, 0}.cp(10⁻³⁰⁰)` NaN, not 20.786), and t itself underflows in the antiderivative (at 5·10⁻³²⁴ K NaN, not 1.04·10⁻³²²); for large T, C·t², D·t³ and t⁴ overflow and their difference is ∞ − ∞ = NaN where the value is ±∞ (H₂O(g)'s Cp at 10³⁰⁰ K NaN, not −∞) | **#20's semantics** on Rust's exact evaluation order, the final addition (Cp) or product by 1000 (H) rounded once, behind a fast test: a box (\|T\| in [2⁻³⁰, 2³⁰] and every nonzero coefficient in [2⁻²⁰⁰, 2²⁰⁰], where every intermediate is provably normal: every tabulated segment over its range, every physical use) or, outside it, each product and quotient normal or an exact zero and each sum finite; T = ±0 (every product an exact zero, E/0 IEEE's special in both) and non-finite inputs keep Rust's chain. The data-backed callers (`enthalpy_change_cp`, `adiabatic_flame_temperature`) share the evaluation and stay in the box (proven for every SHOMATE_DATA segment over its range, so the table-row path tests T alone; replayed for the flame's 324,925 golden evaluations) | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (185 `cp`, 25 `enthalpy_antiderivative`: 210) |
| 35 *(2026-11-11, the overflow audit's B5)* | `reaction::equilibrium_constant` computes exp(−ΔG/(R·T)), `solution::debye_huckel_limiting` 10^((((−A)·z)·z)·√I), `debye_huckel_extended` 10^(((((−A)·z)·z)·√I)/(1 + (B·å)·√I)) and `davies_activity` 10^((((−A)·z)·z)·(√I/(1 + √I) − 0.3·I)): an intermediate of the exp or 10^x argument that over- or underflows gives K or γ = 1, 0 or NaN, or a degraded K, although the value is representable (`equilibrium_constant(−DBL_MAX, DBL_MAX)` 1, not e^0.12 = 1.128: R·T overflows; `(−10⁻³²⁰, 10⁻³²²)` 197,357, not 193,216: a subnormal R·T; `debye_huckel_limiting(100000, 0, 10³⁰⁰)` and `davies_activity` there ∞·0 = NaN, not 1; `debye_huckel_extended(1, 1, 10³⁰⁸, 10¹⁵⁵, 10¹⁵⁵)` 1, not 0.977: B·å = ∞; `(100000, 1, 10³⁰⁰, 10¹⁵⁴, 10¹⁵⁴)` 0, not 10⁻¹⁰⁰: A·z² = ∞), and `equilibrium_constant`'s ComputationError decides on that degraded K (Ok where the value overflows, Err where it does not) | **#20's semantics on the argument**: Rust's sequence with an unbounded exponent (`_km_ux_*`), the argument's final quotient or product rounded once into binary64 (`_km_keq_arg_ux`, `_km_dh_limiting_arg_ux`, `_km_dh_extended_arg_ux`, `_km_davies_arg_ux`), then the unchanged correctly rounded `_km_cr_exp` / `_km_cr_exp10`, which of that binary64 argument is the correctly rounded exp / 10^x of the unbounded 53-bit argument (the lemma in the addition below; no new primitive); behind inline fast tests (R·T in [2⁻¹⁰²¹, DBL_MAX]; A·z² finite; the extended law's numerator and (B·å)·√I finite) that prove Rust's chain gives that result, so those inputs keep Rust's bits. `equilibrium_constant`'s ComputationError is Rust's, judged on the port's K; Rust's checks, Err codes and order are kept | class C (the corrected oracle, K = 0); equal rows plain; the others noted `(f:…)` or, where the port's K overflows and Rust's does not, `(E:ComputationError)`, mapped and counted (326 `equilibrium_constant`, 27 of them the port's Err and 33 its K where Rust's overflows; 55 `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity`: 845) |
| 36 *(2026-11-11, the overflow audit's B6)* | `spectroscopy::fft_spectrum` runs hisab 1.4.0's `num::fft` (bit reversal; per stage wn = (cos θ, sin θ), w *= wn; t = w·d, u ± t) and `Complex::abs` = √(re·re + im·im) of bins 0 … n/2 in binary64: a butterfly sum or a product w·d past DBL_MAX, a product or a difference below the normal range, or a square re² or im² that over- or underflows gives ±∞, 0, NaN or a degraded magnitude where the spectrum is representable (`fft_spectrum([1.75·10³⁰⁴])` [∞], not [1.75·10³⁰⁴]; `([−6.96·10⁻³¹⁶])` [0], not [6.96·10⁻³¹⁶]; `([3.3·10³⁰², 3.0·10³⁰⁴])` [∞, ∞], not [3.03·10³⁰⁴, 2.96·10³⁰⁴]) | **#20's semantics over the whole transform** (the user's decision "FFT: both" — the magnitudes and the butterflies): every butterfly and twiddle operation in hisab's order with an unbounded exponent (`_km_ux_*`), each magnitude's √ the final operation, rounded once into binary64; behind a per-operation fast test on Rust's own transform (`_km_fft`: every twiddle component 0 or at least 2⁻¹⁰⁰, every butterfly's d 0 or at least 2⁻⁹²¹, so each product is an exact zero or at least 2⁻¹⁰²¹; every retained output finite, so no sum overflowed; each bin's square sum in [2⁻⁹⁰⁰, DBL_MAX] or from re = im = 0, else that bin's unbounded magnitude), so in-range rows keep Rust's bits; the slow path `_km_fft_ux` on 32·n heap bytes. A NaN or ±∞ sample keeps Rust's IEEE spectrum (B1's rule); Rust's length check and Err kept; ADR 0003's M4e transcription gains the unbounded slow path (its amendment) | class C (bound 0, the corrected oracle: no fft row has a correction); equal rows plain; the others noted with the port's whole spectrum `(v:n f:… …)`, mapped and counted (418 rows: 41 of the 213 rows before B6, 377 of B6's 515 edge rows) |
| 37 *(2026-11-11, the overflow audit's B6)* | `spectroscopy::find_peaks` tests w > RN(RN(min·RN(1 − t)) + RN(max·t)): the rounded threshold decides otherwise for a candidate within a few ulps of the exact min + t·(max − min), at the subnormal grid (`find_peaks([−0, 5·10⁻³²⁴, 0], 0.75)` [], the exact threshold 3.7·10⁻³²⁴), and as far as 2⁵² ulps off where min = −DBL_MAX and max = DBL_MAX (1 − t rounded by 2⁻⁵⁴, times DBL_MAX); `([10, 10.2, 10, 11, 10], 1 − 2⁻⁵³)` [], the exact threshold 11 − 2⁻⁵³ | the **exact predicate** (settled by the main session under the user's policy: "find_peaks' threshold is evaluated exactly"): w > min·(1 − t) + max·t in rational arithmetic, t the given double (1 − t exact), min and max Rust's fold over the finite values, the neighbour tests IEEE, a +∞ candidate above every (finite) threshold; Rust's threshold decides wherever a proven bound separates w from it (\|RN(w − T)\| > B, B = RN(RN(\|min\|·2⁻⁵¹) + RN(\|max\|·2⁻⁵¹)) + 2⁻¹⁰⁷², against the proven error 2.5u(1 + u)(\|min\| + \|max\|) + (1 + u)·2⁻¹⁰⁷⁴), else `_km_peak_above_exact` (the sign of w − min + t·min − t·max, six terms at any exponent, by `src/mathx.cyr`'s new `_km_sign_scaled`). Rust's checks and empty results kept | class A; equal rows plain; the others noted with the port's indices `(v:n i:… …)`, mapped and counted (222 rows: 3 of the 994 rows before B6, 219 of B6's 2,393 edge rows) |
| 38 *(2026-11-11, the overflow audit's B6; the survey's F16, its flame part — B12 extends this row with the weak acid / base Newton and Peng–Robinson's coefficients)* | `thermochem::adiabatic_flame_temperature` tests the bracket as f(T₀)·f(T_max) < 0, a product that underflows to −0 when \|ΔH\| and the moles are tiny, so a reaction whose balance changes sign is "not bracketed", Err(ComputationError) (`(−10⁻³²⁰, [CO₂ ×10⁻³⁰⁰], 298.15)`; CH₄/air scaled by 10⁻²⁰⁰), and hisab 1.4.0's transcribed bisection tests f(a)·f(b) > 0 likewise | the **exact sign test** (the user's decision "Fix all three iterative edges", the flame's being "the flame bracket's exact sign test"): f(T₀) < 0 < f(T_max) or f(T_max) < 0 < f(T₀), Rust's answer wherever the product does not underflow to a zero (NaN, ±∞ and ±0 operands included, proven), and the bisection's same-sign test likewise (it cannot change a result: the flame calls it only after the bracket proved opposite signs at the same points). Newton, its absolute 10⁻³ J/mol tolerance, every Err and Rust's order kept, so such a system gets Newton's initial guess, Rust's own answer at a scale where the product does not underflow (a known, kept Rust behaviour: at any scale below about 10⁻³ the guess is accepted unconverged) | class B; equal rows plain; the others noted `(f:…)` with the solver path, mapped and counted (273 rows: 16 of the 1,937 rows before B6, Rust's Err there; 257 of B6's 443 edge rows) |
| 39 *(2026-11-12, the overflow audit's B7; the survey's F15)* | `stoichiometry::balance_equation` runs its fraction-free Gauss–Jordan elimination on i64 with checked arithmetic: where a product pivot·a or factor·a, their difference, the pivot row's negation, a reduced row's gcd (`i64::try_from` of 2⁶³), the lcm of the pivots or x[c] = −(row[f]·(lcm/piv)) leaves i64, it returns Err(ComputationError), "integer overflow while balancing equation", although the null space can be one-dimensional with a positive u32 vector (`[[3037000500, 1, −3037000501], [1, 3037000500, −3037000501]]` Err, its first product 3037000500² past i64::MAX; the answer [1, 1, 1]) | **Rust's algorithm on unbounded integers** (settled by the main session under the user's policy, not asked: "Balancing: multi-word integers, no size limit", the survey's remedy (a)): every step as written — the pivot search, the swap, reduce_row, the sign flip, row_r ← pivot·row_r − factor·row_pivot then reduce_row, the lcm of the pivots, x[c] = −(row[f]·(lcm/piv)), reduce_row(x), the zero, sign and u32 checks, the i128 verification — with mathematical integers, so every overflow Err disappears and every other Err (ambiguous, only the trivial solution, a zero coefficient, mixed signs, past u32, the verification, each InvalidInput, the cap check after) stays at its place in Rust's order. Rust's i64 code runs first, unchanged (`_km_balance_i64`), and at its first overflow the same steps start again from the original matrix on signed multi-word integers (`_km_balance_big`, `_km_sb_*`: 31-bit limbs, Knuth's division, Lehmer's gcd), their tables sized from a proven Hadamard bound (`_km_balance_bits`), so only allocation limits the size | class A; equal rows plain; the others (Rust's overflow Err, the port's Ok) noted with the port's whole result `(v:n i:… …)`, mapped and counted (413 rows: 104 of the 3,329 rows before B7, 309 of B7's 521 edge rows) |

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

**Reviewed and kept as Rust's** (reported to the user, the policy's borderline; the first four, and `reduced_mass`'s
m₁·m₂ and m₁ + m₂, decided otherwise in M4e-3 as #21 and #22): `reduced_mass(m, m)` — Rust's
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
k/μ, h·B, 8π²·I, μ·r·r, the FFT's sums and re² + im², and the like) are its first candidates; the helpers above are its
machinery. The audit is not part of #20. (`reduced_mass`'s m₁·m₂ and m₁ + m₂ left the list in M4e-3: #22 covers them.)

## Addition, 2026-11-10 (M4e-3): #21 and #22

**The user's decisions (2026-11-10), verbatim.** The questions were the main session's; the answers are the user's.
1. Q: "M4e found that Rust's Rydberg factor 1/n₁² − 1/n₂² cancels badly when n₁ and n₂ are close and large:
   rydberg_wavelength is off by 39 ulp at n = 100, 939 at n = 10⁴ and ~4·10⁷ at n ≈ 2³¹, and bohr_transition_energy's
   E₂ − E₁ behaves the same. Compute the factor exactly ((n₂² − n₁²)/(n₁²·n₂²) in integers, rounded once) as a new ADR 0002
   row?" A: **"Compute exactly (Recommended)"** — the option's text: "Exact integer factor, one rounding; both functions
   become correctly rounded in n. Results move by up to ~4·10⁷ ulp at extreme n, by 0–1 ulp for ordinary lines."
2. Q: "M4e found that two exact values come out 1 ulp wrong: reduced_mass(m, m) misses m/2 on 8.6% of masses (Rust:
   m₁·m₂/(m₁+m₂), three roundings), and bohr_energy_level(Z, n) misses −R·Z²/n² (Z = n misses −R on 3 of 12 rows). Make
   both correctly rounded (rounded once)?" A: **"Round once (Recommended)"** — "reduced_mass and bohr_energy_level return
   the correctly rounded value for every input, so the homonuclear m/2 and Z = n cases come out exact."
3. The main session then told the user: "For the Rydberg wavelength and Bohr transition energy, I plan to round the whole
   expression once, not just the exact factor, to match your 'round once' choice for the Bohr levels. Say so if you'd
   rather keep Rust's remaining steps." The user did not object and replied "Committed, continue". So the whole
   expression is rounded once.

**#21, the target.** Each function returns its exact rational value, the constant taken as its binary64 value, rounded
once to binary64 (nearest, ties to even): `bohr_energy_level` RN(−R·Z²/n²) with R = `KM_RYDBERG_EV`;
`bohr_transition_energy` RN(R·Z²·(n₂² − n₁²)/(n₁²·n₂²)); `rydberg_wavelength` RN(n₁²·n₂²/(R∞·Z²·(n₂² − n₁²))) with
R∞ = `KM_RYDBERG`; `SpectralSeries::wavelength` the same with Rust's const-evaluated R_H and Z = 1. Rust's validation, Err
codes and check order are unchanged (the transition energy: the level checks of n₁, then of n₂, then n₂ ≤ n₁). Every value
is a normal double (|E| ≥ R·2⁻⁶⁴, λ ≤ 2⁷²). On the corpus Rust's levels are at most 2 ulp from the exact value, its series
lines 2 ulp, its wavelengths up to 9.0·10⁷ ulp and its transition energies up to 1.8·10⁹ ulp (the cancellation of two
levels near n = 2³²; the decision's "~4·10⁷" was the wavelength measured in M4e on a smaller grid).

**#22, the target.** `reduced_mass(m₁, m₂)` returns RN(m₁·m₂/(m₁ + m₂)) for every pair of positive finite doubles, subnormal
results and ties included (μ lies in [min/2, min), so it never overflows); no intermediate overflows or underflows. Rust's
NaN at (DBL_MAX, DBL_MAX) and its 0 or ∞ at the range edges are fixed by this; Err for a mass that is not positive and
finite, as Rust's. The library-wide overflow audit's candidate list (#20's addition) loses `reduced_mass`, which this row
covers.

**The implementation** (ADR 0006's two-phase style).
- *The general helper* (`src/mathx.cyr`, for reuse by M4f's stoichiometry and the audit): multi-word unsigned integers,
  8 limbs of 32 bits in i64 words (below 2²⁵⁶; `_km_bn_set`, `_copy`, `_bitlen`, `_cmp`, `_add`, `_sub`, `_mul_small`,
  `_mul`, `_shl`; every limb product plus carry below 2⁶⁴, formed in wrapping 64-bit arithmetic and split with `&` and the
  logical shift), and `_km_cr_ratio(p, q, e)` = RN(p/q·2^e) for p, q < 2²⁵⁴: N = p and D = q aligned to one bit length (and
  N doubled once more if N < D), so p/q = (N/D)·2^t with 1 ≤ N/D < 2; 55 steps of restoring division give Q = ⌊N·2⁵⁴/D⌋
  in [2⁵⁴, 2⁵⁵) and the remainder; then `_km_cr_from_int(2Q + sticky, e + t − 55)` rounds once, ties to even, at the
  binade's or the 2⁻¹⁰⁷⁴ grid, ±∞ past DBL_MAX. That is exact: the value lies in [2Q, 2Q + 2) units, strictly inside when
  the remainder is not 0, and so does 2Q + sticky, while every rounding boundary there (at least 3 of the 56 bits are
  dropped) is a multiple of 4 units, an end point at most.
- *The fast paths* evaluate the exact numerator and denominator as double-doubles and round their quotient with
  `_km_dd_div_cr`: y0 = RN(n_h/d_h), the residual n − y0·d from an exact TwoProduct (n_h − a exact by Sterbenz) to
  7.04u²·n (u = 2⁻⁵³), d = RN(r/d_h), so y0 + d is n/d to 13.12u²; with operands carrying relative errors up to 32u² in
  all, y0 + d is the exact value to 45.13u², and the test RN(y0 + RN(d − ε)) = RN(y0 + RN(d + ε)) with ε = 2⁻¹⁰⁰·y0 =
  64u²·y0 (≥ 63.9u² of the value, the inner roundings 3.04u²) brackets the value: equal ends prove the rounding; exact
  values pass, midpoints never do; it fails with probability about 2⁻⁴⁶. The operands: R·Z² exact (`_km_rz2`: Veltkamp's
  split of R, both halves times Z² < 2¹⁶ exact, FastTwoSum); n² an exact TwoProduct, or for n < 2¹³ an exact double below
  2²⁶, by which the level's quotient is inline (y0·n² split exactly, to 5.04u²); for the Rydberg pair, exact when n₂ < 2¹³
  (Z²·(n₂² − n₁²) < 2⁴² and (n₁n₂)² < 2⁵² are integers, R times the first a TwoProduct), otherwise (n₂ − n₁)(n₂ + n₁),
  n₁·n₂ exact TwoProducts and one `_km_dd_mul` each (within 8.01u², 2⁻¹⁰²·⁹⁹: the two cross products' and their sum's
  roundings, e + c's, the dropped a_l·b_l).
- *The reduced mass*: with m₁ ≥ m₂, m = s·2^e and g = e₁ − e₂ (from the exponent fields, a lower bound for a subnormal
  m₂): for g ≥ 55 the result is m₂, since x = m₂/m₁ < 2⁻⁵⁴ makes m₂ − μ = m₂x/(1 + x) < 2⁻⁵⁴·m₂, under half the gap below
  m₂ (at least 2⁻⁵³·m₂, or 2⁻¹⁰⁷⁴ for a subnormal m₂), with μ < m₂: strictly inside m₂'s rounding interval, never a tie.
  Otherwise μ = 2^e₂·s₁s₂/(s₁ + s₂·2⁻ᵍ): the numerator an exact TwoProduct in [1, 4), the denominator an exact FastTwoSum
  in [1, 4) (s₂·2⁻ᵍ ≥ 2⁻⁵⁴ exact), the quotient F in [½, 2) by `_km_dd_div_cr`, and RN(μ) = RN(F)·2^e₂ exactly for
  e₂ ≥ −1021; below (μ possibly subnormal) and where the test fails, the exact path: m = M·2^E with 2⁵² ≤ M < 2⁵³ (a
  subnormal's significand shifted up), μ = M₁M₂·2^E₁/(M₁·2^G + M₂) by `_km_cr_ratio` (below 2¹⁰⁶ and 2¹⁰⁸).
- *Where ties exist* (the slow path's job). R's odd significands decide: `KM_RYDBERG_EV`'s is 683·11214237642539, squarefree
  with a prime above 2³², so a level is a double exactly when n/Z is a power of two (then −R/4^k) and otherwise has at least
  56 significant bits — no level is a midpoint; `KM_RYDBERG`'s (3·1963825758347083) and R_H's (887·3319205424319) have a
  prime above 2³², which no n₁, n₂ < 2³² cancel, so no wavelength is even dyadic. The transition energy with n₁ = g·x,
  n₂ = g·y (x, y coprime) is dyadic iff odd(x)·odd(y) | odd(Z) and odd(g)² | 683·P·j²·w (w = odd(y² − x²),
  j = odd(Z)/(odd(x)·odd(y))): an exhaustive search over that structure finds 5,657 exact values (16R at (15; 3, 5), 32R at
  (17; 3, 51), …) and 13 ties, (255; 2^k·58,055, 2^k·524,544) = (255; 2^k·683·85, 2^k·683·768) for k = 0 … 12 (x = 85,
  y = 768, w = 683·853, odd(g) = 683: the significand 853·P has 54 bits). The reduced mass has ties in families:
  subnormal ones (m = m at odd multiples of 2⁻¹⁰⁷⁴, (6c, 2c)·2⁻¹⁰⁷⁴) and normal ones (μ = Δαβ·2^(p+q−v) for coprime odd
  α, β with α·2^p + β·2^q = H·2^v, masses HΔα·2^p and HΔβ·2^q; e.g. (2²⁷ − c, 2²⁷ + c) for odd c ≤ 94,906,265).

**Golden comparison** (as #16–#20): the corpus keeps the oracle's values. Where the port's differs the note names it,
`[ADR 0002 #21: the port returns the correctly rounded −R·Z²/n² = <x> (f:…)]` (likewise `R·Z²·(n₂² − n₁²)/(n₁²·n₂²)`,
`n₁²·n₂²/(R∞·Z²·(n₂² − n₁²))`, `n₁²·n₂²/(R_H·(n₂² − n₁²))`) and `[ADR 0002 #22: the port returns the correctly rounded
m₁·m₂/(m₁ + m₂) = <x> (f:…)]`, computed by the generator in exact integer arithmetic (`exact::cr_ratio`, ⌊p·2^k/q⌋ to 60
bits with a sticky remainder; the reduced mass at any exponent gap, without the port's cutoff). `validate.py` proves
every Ok row with `Fraction`: noted exactly where the exact value differs from the oracle's, the note naming it bit for
bit, so every unnoted row is correctly rounded too; it re-derives the structure (R's odd significands factored, the
transition energy's search with its own factorisation) and checks that every exact level (6,869 rows, n = 2^k·Z) and
every transition tie is in the corpus. The harness maps the noted rows (accepted under ADR 2, counted per function):
3,823 `bohr_energy_level`, 1,351 `bohr_transition_energy`, 1,376 `rydberg_wavelength`, 67 `SpectralSeries::wavelength`
(#21: 6,617) and 2,837 `reduced_mass` (#22). The grids gain every exact level, n₁ = n₂ − 1 at every binary scale and next
to u32::MAX, the 2¹³ edge, the 13 ties and each exact-value family's smallest, middle and largest member, the series'
large lines, and for the reduced mass m = m at every power of two and subnormal multiple, the tie families at scales
across the range, exponent gaps 45 … 64 around the proven 55, the 2⁻¹⁰²¹ edge and both range ends (195 tie rows).
`cr_verify.py` checks `_km_cr_ratio` (`ratio`) and the five functions (`bohr`, `trans`, `ryd`, `series`, `rmass`: the public
function, the fast path and the exact path forced) against `Fraction` on x86_64 ELF, x86_64 PE and aarch64.

**Cost** (`tests/kimiya.bcyr`, before → after, the same host): `spectroscopy/bohr_energy_level` 9 → 49 ns (Rust 1.21),
`bohr_transition_energy` 27 → 73, `rydberg_wavelength` 24 → 73, `series_wavelength_balmer_alpha` 28 → 78,
`reduced_mass_hcl` 14 → 73; the exact path on a tie: `bohr_transition_energy_tie` 2.42 µs, `reduced_mass_tie` 2.20 µs.

## Addition, 2026-11-10 (M4e-4): #23

**The user's decision (2026-11-10), verbatim.** The question was the main session's; the answer is the user's.
Q: "Rust's isotope table rounds the AME2020 masses to 5 decimals (¹H 1.00783 vs 1.007825032, ²H 2.01410 vs 2.014101778)
and omits ³⁶S (abundance 0.0001, so sulfur sums to 0.9999). Upgrade it to full-precision AME2020 masses and the complete
IUPAC compositions (³⁶S included)?" A: **"Upgrade the data (Recommended)"**. The agreed form: the full AME2020 masses plus
³⁶S, from a cited source file that the generator reads, recorded as this row.

**The sources** (retrieved 2026-11-10; the hashes are of the files as downloaded):
- **Masses: AME2020**, the unrounded atomic-mass table `mass_1.mas20.txt`
  (https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt, sha256
  `e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307`): M. Wang, W.J. Huang, F.G. Kondev, G. Audi and S.
  Naimi, "The AME 2020 atomic mass evaluation (II). Tables, graphs and references", *Chinese Physics C* 45, 030003 (2021).
  Each mass is read by the file's own Fortran format (its header: the atomic mass is an i3 in u at columns 107–109 and an
  f13.6 in micro-u at 111–123) and kept with its published digits: "1 007825.031898" µu is 1.007825031898 u.
- **Abundances: TICE 2009**, Table 1, column 9, "Representative isotopic composition"
  (https://www.ciaaw.org/pubs/TICE-2009.pdf, sha256 `2145b6762c77f229c4594f5bc4bd4f6d98e1194c0a554f7a5e11881f1db2c893`):
  M. Berglund and M.E. Wieser, "Isotopic compositions of the elements 2009 (IUPAC Technical Report)", *Pure Appl. Chem.*
  83(2), 397–410 (2011), doi:10.1351/PAC-REP-10-06-02. It is Rust 1.2.0's cited source ("IUPAC 2009 representative
  isotopic compositions"); its 16 values for Rust's isotopes are Rust's literals digit for digit, and it gives ³⁶S
  0.0001(1).

**Why TICE 2009 and not CIAAW's current table.** CIAAW's current "Isotopic compositions of the elements 2024"
(https://www.ciaaw.org/isotopic-abundances.htm, retrieved 2026-11-10, sha256
`1ab770a48bb539d26ad2992a94b4c60d51743058b67f106796b754cd5035423b`) gives all seven of these elements intervals, not
values: ¹H [0.999 72, 0.999 99], ¹²C [0.9884, 0.9904], ¹⁴N [0.995 78, 0.996 63], ¹⁶O [0.997 38, 0.997 76], ³²S [0.9441,
0.9529], ³⁵Cl [0.755, 0.761], ⁷⁹Br [0.505, 0.508] (the page cites the 2013 report, J. Meija et al., *Pure Appl. Chem.*
88, 293–306 (2016), and the revisions since). A mass spectrum needs one composition and the current table offers none
for these elements, so the abundances stay TICE 2009's single values — Rust's, with ³⁶S added — and the decision changes
the masses and completes sulfur. As a consistency check, each element's abundance-weighted mean mass lies in its CIAAW
standard-atomic-weight interval where the CIAAW page was retrieved with the sources (H 1.007941 in [1.007 84, 1.008 11],
https://www.ciaaw.org/hydrogen.htm, sha256 `eeb33c101263ab3047789c45758f2b2d63fd60622fd8e135c0af2d1a4e109d3b`; S 32.064787
in [32.059, 32.076], https://www.ciaaw.org/sulfur.htm, sha256
`5103098a7605eccb73932ec5058fd0b3a9a138d5089f2fd7acb7e52b730b5c12`).

**#23, the table** (the 17 isotopes; the masses in u, the abundances as mole fractions):

| Isotope | AME2020 mass (Rust 1.2.0) | TICE 2009 abundance |
|---|---|---|
| ¹H | 1.007825031898 (1.00783) | 0.999885 |
| ²H | 2.014101777844 (2.01410) | 0.000115 |
| ¹²C | 12.0000000 (12.0) | 0.9893 |
| ¹³C | 13.00335483534 (13.00335) | 0.0107 |
| ¹⁴N | 14.00307400425 (14.00307) | 0.99636 |
| ¹⁵N | 15.00010889827 (15.00011) | 0.00364 |
| ¹⁶O | 15.99491461926 (15.99491) | 0.99757 |
| ¹⁷O | 16.99913175595 (16.99913) | 0.00038 |
| ¹⁸O | 17.99915961214 (17.99916) | 0.00205 |
| ³²S | 31.97207117354 (31.97207) | 0.9499 |
| ³³S | 32.97145890862 (32.97146) | 0.0075 |
| ³⁴S | 33.967867011 (33.96787) | 0.0425 |
| ³⁶S | 35.967080692 (—) | 0.0001 (—) |
| ³⁵Cl | 34.968852694 (34.96885) | 0.7576 |
| ³⁷Cl | 36.965902573 (36.96590) | 0.2424 |
| ⁷⁹Br | 78.918337574 (78.91834) | 0.5069 |
| ⁸¹Br | 80.916288197 (80.91629) | 0.4931 |

Each value in `src/spectroscopy_data.cyr` is the double nearest its decimal, written as its shortest round-trip decimal
(12.0000000 is 12.0). `isotope_pattern` is Rust's function operation for operation (the pruned convolution, the stable
sort, the binning, the thresholds, the Err codes and their order); it convolves this table, so every pattern with an
element that has data moves (the masses by up to the rounding Rust's literals dropped, times the atom count; sulfur's
intensities from M+4 on, where ³⁶S joins the bin: S₈'s M+4 is 0.0573616 against 0.0565194 without it). Sulfur's scratch is
4·n records for n pattern peaks (S₈ at most 68, within the 192-record stack tables).

**The source file and its checks.** `scripts/oracle/isotopes.tsv` holds one record per isotope — Z, A, symbol, the mass
in u, the abundance, the TICE value as printed (with its uncertainty) and the AME2020 line verbatim — under a header with
the URLs, sha256s and references; and the two CIAAW intervals. Before writing anything the generator
(`t_spectroscopy::isotope_table`) re-reads each mass from its verbatim line by the Fortran columns and each abundance
from its printed value, and fails unless (a) its elements and isotopes are the oracle's `isotope_abundances(z)` at every
u8 z (an isotope's A its mass rounded), plus ³⁶S as the only addition; (b) every shared abundance is Rust's literal as a
decimal string; (c) every Rust mass is the AME2020 mass rounded to five decimals (exact decimal arithmetic; a tie would
be refused); (d) each element's abundances sum to exactly 1 in decimal arithmetic. `validate.py` repeats (a)–(d) on its
own parse, checks the generated table against the file and the mean masses against the intervals.
`scripts/oracle/isotope_sources.py . DIR` re-verifies the file against fresh downloads: the two sha256s, every recorded
line a whole line of `mass_1.mas20.txt` (its only line for that Z and A), and every TICE value a cell of the PDF's Table 1
in record order (the content streams inflated, the text runs split at the column gaps).

**Golden comparison** (as #16–#22, with the note extended to a whole result): the corpus keeps the oracle's values.
Where the port's result differs the note names it, `[ADR 0002 #23: with the AME2020 masses and the TICE 2009 abundances,
³⁶S included, the port returns [n], m, i, … (v:n f:… f:… …)]` — the trailing parenthesis holds the result's golden
tokens, where #13–#22's hold one `f:` or `i:` token. The generator computes it with `m_spectroscopy::pattern_replica`,
Rust's `isotope_pattern` transcribed with the table a parameter, which must give the oracle's result (bits or Err) on
every row with Rust's table before it is run with this one. `validate.py` proves every note with its own replica
(`r_pattern` with the table swapped), and that a row is noted exactly where the port's result differs from the oracle's,
so every unnoted row's port result is the oracle's; its pinning checks of the three typed thresholds run against the
port's expected results too. The harness compares a noted row's whole result with the note's tokens (`_par_map_toks`),
bit for bit, and counts it as accepted under ADR 2: 7 `isotope_abundances` rows (every element with data) and 347
`isotope_pattern` rows (every Ok row with an element that has data at a count above 0; the 181 Err rows and the 168 Ok
rows with no such element are unchanged). `probe_tables.py` compiles the generated cells against the notes' values.

**Cost** (`tests/kimiya.bcyr`'s isotope benches, the old table against the new, the same host, interleaved):
`spectroscopy/isotope_pattern_glucose` 15.8 µs either way (Rust 2.80); the Cyrius-only `isotope_pattern_cystine`
(C₆H₁₂N₂O₄S₂) 16.9 → 18.5 µs and `isotope_pattern_s8` 8.4 → 18.8 µs: sulfur's fourth isotope lengthens each step's
scratch by a third and widens the pattern (offsets 0, 1, 2, 4 instead of 0, 1, 2), so S₈'s steps carry up to 17 peaks
instead of 15 and its insertion sort moves each ³⁶S record further.

## Addition, 2026-11-10 (M4f): #7 realised

**#7, as planned.** Row #7 has stood since 2026-11-08: `assign_oxidation_states` returns exact rational states, and
`balance_equation` takes an integer matrix. M4f ports `stoichiometry` with both, and settles the points the row left open
(each the option closest to Rust's observable behaviour; the port plan's §6.4 "Decided in M4f" lists them too).

**The rational states.** Each entry's state is a record (Z, numerator, denominator), the denominator positive and the
fraction reduced (`assign_oxidation_states(atoms, n, charge, out, cap)`: n 16-byte (Z, count) records in, n 24-byte
records out, `Ok(n)`). Rust's algorithm is unchanged — the rules, the single unknown element, the peroxide re-solve and
its −2 < os < 0 test, every Err, its code and its order — with Rust's private `exact_oxidation_state` returning the
rational where Rust returned Err "non-integer average oxidation state": C₃H₈'s C −8/3, Fe₃O₄'s Fe +8/3 (the reference
tests), O₂⁻'s −1/2, and in the re-solve branch KO₂'s O −1/2 and KO₃'s ozonide −1/3, the −2 < p/q < 0 test made on the
rational (as p < 0 < p + 2q, formed without overflow). An oracle `Ok(n)` is the port's n/1. Settled:
- **The range is the value's.** Rust's `i32::try_from` refused an integer state outside i32 ("oxidation state out of
  range", ComputationError); the port applies that bound to the rational's value (p/q outside i32::MIN ..= i32::MAX is
  Err(ComputationError)), its numerator and denominator are i64 words. Rust's observable result there is
  Err(ComputationError) either way — its non-integer Err comes first — so these inputs keep Rust's code; they need counts
  near 2³¹ (Mn₃ O_(2³²) is 2³³/3, past i32::MAX; Mn₂ O_(2³¹ − 1) at charge −1 is 2147483646.5, inside). Below that bound
  every non-integer Err of Rust's is the port's `Ok(p/q)`.
- Rust's u8 Z, u32 count and i32 charge are words, never wrapped (port plan §6.4): a Z outside 1 ..= 118 is
  Err(InvalidElement), a count outside 1 ..= u32::MAX Err(InvalidInput) at Rust's `count == 0` check, a charge outside
  i32 Err(InvalidInput), checked first.

**The integer matrix.** `balance_equation(matrix, m, n, out, cap)`: m element rows by n species columns, flat row-major
(entry (i, j) at word i·n + j), the coefficients to the caller's buffer (`Ok(n)`). Rust's fraction-free Gauss–Jordan
elimination is kept operation for operation: checked i64 products, differences and negations (`_km_i64_*_ck`: Cyrius's
+ − × wrap, so each tests its wrapped result; i64::MIN / −1, which would trap, is refused first), Rust's `reduce_row`
with its unsigned gcd to 2⁶³ (a row of gcd 2⁶³ is Rust's `try_from` overflow, a corpus row), the null vector from the
RREF and the lcm of the pivots, the sign convention and the u32 conversion, and the verification of every original row in
i128 with Rust's `checked_add` semantics. `src/mathx.cyr`'s `_km_bn_*` integers are unsigned, so the i128 sum is a small
signed helper (`_km_i128_mul_u32`, `_km_i128_add_ck`: h·2⁶⁴ + m·2³² + l with h an i64, exactly i128's range; the edge
where h₁ + h₂ overflows negatively by one and the carry brings it back to −2¹²⁷ is kept); Rust's i128 `checked_mul`
cannot fail (|entry| ≤ 2⁵³, k < 2³²), nor the product here (below 2⁹⁵). Every Rust overflow Err is an Err at the same step.
Proven fast paths change nothing: four operands below 2³¹ make an elimination step overflow-free, and a row whose i64
sum of entries below 2³¹ never overflows has Rust's i128 sum. Settled:
- **The entry range is Rust's**: |entry| ≤ 2⁵³ (MAX_COMPOSITION_ENTRY), a larger one Err(InvalidInput). Rust's bound
  exists because of f64, but keeping it accepts every integer matrix Rust accepts and refuses every one it refuses, and
  keeps the i64 elimination's overflow points Rust's.
- **The rows with no counterpart.** Rust's f64 entries within 10⁻⁹ of an integer are that rounded integer (the harness
  converts them as Rust's test does, so those rows compare class A: 41 rows); a matrix with an entry Rust's test refuses
  (not within 10⁻⁹ of an integer, or non-finite), an integer beyond i64, or ragged rows has no integer counterpart. Rust
  returns Err(InvalidInput) for all of them (the generator asserts it). Those rows keep the oracle's output and are noted
  `[ADR 0002 #7: no integer-matrix counterpart (<why>): …]`; the harness requires both the note and the oracle's
  Err(InvalidInput), and counts the row as accepted under ADR 2 (78 rows: 55 non-integer, 11 non-finite entries, 7
  integers beyond i64, 5 ragged). An integer-valued f64 entry past 2⁵³ but within i64 (2⁵³ + 2, 2⁶², 10¹⁸) has a
  counterpart, which the port refuses with Rust's code: those rows are plain.
- A null matrix, m < 1, n < 2 or m·n past any buffer is Err(InvalidInput); the cap is checked after Rust's errors, `out`
  untouched (M4e's rule).

**Golden comparison.** The oxidation rows where the oracle's Err is the non-integer one and the port is Ok carry the
port's whole result, `[ADR 0002 #7: the port returns the exact average oxidation states (Z, numerator, denominator) [n],
… (v:n i:z i:p i:q …)]` (#23's whole-result note), from `m_stoichiometry::port_oxidation`, Rust's algorithm on exact
rationals, which must give the oracle's result on every other row (asserted on every run); the harness compares the
note's tokens (`_par_map_toks`, now taking the note's tag) and counts the row as accepted under ADR 2: 9,090 rows.
`validate.py` recomputes every row from Rust's semantics, recomputes the port's states with `Fraction`, and proves each
note is exactly where Rust's Err is the non-integer one and the port is Ok, its tokens the port's (reduced, q > 0); it
checks every balancing row by exact linear algebra and every no-counterpart note against its own conversion.

**The other five functions are ported operation for operation, Rust's arithmetic unchanged,** and measured under the
inexact-formula policy (the user decides; nothing here is fixed). On the golden corpus and seeded grids (the
measurement script, M4f's report): `percent_yield`'s (a/t)·100 is one ulp from RN(100a/t) on 24 % of the corpus's Ok rows
(27 % of 200,000 seeded pairs), up to 50 ulp and 0 for a subnormal result (6 corpus rows; 1,656 of 20,000 underflow-edge
pairs); `theoretical_yield_moles`' (m/c)·p is one ulp from RN(m·p/c) on 15–26 % of seeded triples, and ±∞ or 0 where
m/c over- or underflows though m·p/c is finite (4,773 and 2,722 of 10,000 edge triples each); `limiting_reagent`'s rounded
ratios pick another index than the exact m_i·c_j < m_j·c_i on 227 of 2,852 corpus rows (the deliberate near-ties and
over- and underflows), 21 % of seeded near-ties, about half of the pairs whose ratios both overflow to +∞ or both
underflow, and none of 100,000 random lists; `percent_composition`'s ((mass·count)/mw)·100 is up to 2 ulp from
RN(mass·count·100/mw) given Rust's mw on 26 % of the corpus's entries; `empirical_formula` with exact ratios and an exact
1/10 tolerance changes Ok ↔ Err on 14 of 3,351 corpus rows and 45 of 50,000 formulas' own percentages (never the
multiplier or a count alone); `balance_equation`'s checked i64 elimination fails where a valid u32 solution exists
([[3037000500, 1, −3037000501], [1, 3037000500, −3037000501]]: x = [1, 1, 1], but the first step forms
3037000500² − 1 > i64::MAX), a case for the library-wide overflow audit after M4. The oxidation states' checked sums can
overflow only with 2³⁰ or more entries.

## Addition, 2026-11-10 (M4f-2): #24, #25 and #26

**The user's decisions (2026-11-10), verbatim.** The questions were the main session's; the answers are the user's.
1. Q: "In M4f's stoichiometry, percent_yield (a/t·100), theoretical_yield_moles (m/c·p) and percent_composition
   (mass·count/mw·100) round two to four times. They are 1–2 ulp off on about a quarter of inputs, and at the range edges
   give 0 or ∞ where the answer is representable: percent_yield(5e-324, 10) gives 0, the true value is 5e-323. Return each
   exact rational rounded once? percent_composition would keep molecular_weight's own mw." A: **"Round once
   (Recommended)"** — the option's text: "All three become correctly rounded for every input, range edges included (like
   reduced_mass in #22). About a quarter of results move by 1–2 ulp."
2. Q: "limiting_reagent compares the rounded ratios moles/coefficient. So it can pick the wrong reactant when two ratios
   are within an ulp, or when both overflow to ∞ or underflow: [DBL_MAX, DBL_MAX] over coefficients [5e-324, 1e-300] gives
   index 0, while the exact answer is 1. Compare the ratios exactly (mᵢ·cⱼ against mⱼ·cᵢ, the first index on an exact
   tie)?" A: **"Compare exactly (Recommended)"**.
3. Q: "empirical_formula tests |r·k − round(r·k)| ≤ 0.1 on rounded mole ratios. Evaluated exactly, the result flips
   between Ok and Err only for inputs sitting exactly at the ±0.1 edge: 14 of 3,351 golden rows. For example C 12.011 %
   with O 17.5989 % is an Err in Rust but exactly CO. Evaluate the ratios and the tolerance test exactly?" A: **"Evaluate
   exactly (Recommended)"** — "Exact rational ratios and an exact 1/10 tolerance, as with the heuristic thresholds before.
   Only knife-edge inputs change; multipliers and counts never change otherwise."

**The targets.** Every Rust validation check, Err code and check order stays as it is; only the arithmetic changes.
- **#24**: `percent_yield` returns RN(100·a/t), `theoretical_yield_moles` RN(m·p/c) with m and c at #25's index, and
  `percent_composition` RN(mass·count·100/mw) per atom, mw being `molecular_weight`'s own double — RN to nearest, ties to
  even, subnormal results on the 2⁻¹⁰⁷⁴ grid, ±∞ only where the exact value rounds past DBL_MAX. A zero operand (an
  actual yield or a limiting reagent's moles of ±0, a count of 0) gives Rust's own ±0.
- **#25**: `limiting_reagent`'s minimum i₀ is replaced by i when mᵢ·c_{i₀} < m_{i₀}·cᵢ exactly; Rust's strict `<` keeps the
  first index on an exact tie. An exact ratio is never ∞, so element 0 always seeds the minimum (Rust's +∞ seed is
  unobservable: it leaves index 0 too).
- **#26**: `empirical_formula`'s moles pctᵢ/massᵢ, their minimum (element j), the ratios rᵢ = pctᵢ·mass_j/(massᵢ·pct_j),
  the multiplier — the smallest k in 1 ..= 6 with |rᵢ·k − round(rᵢ·k)| ≤ 1/10 for every i — and the counts round(rᵢ·k)
  (halves away from zero) are exact; past u32::MAX a count is Rust's Err(ComputationError).

**The implementation** (src/stoichiometry.cyr, on src/mathx.cyr's #21/#22 machinery; each bound proven in the source).
- *RN(x·y/z)*, `_km_cr_mul_div`, for positive finite x, y, z (the three #24 functions with (a, 100, t), (m, p, c) and
  (mass, 100·count, mw); 100·count < 2³⁹ is exact). The fast path takes the significands s in [1, 2) and k = eₓ + e_y − e_z:
  N = sₓs_y is an exact TwoProduct p_h + p_l, D = s_z exact, and `_km_dd_div_cr`'s steps (inline, d_l = 0) give F = N/D
  in (½, 4) with exact operands, so its rounding test proves RN(F) whenever it passes. Where it fails, F lies strictly
  between two points rounding to adjacent z₁ < z₂, and the only boundary between is their midpoint μ: the sign of
  N − D·μ = (p_h − a) + p_l − b − h (D·z₁ = a + b a TwoProduct, h = D·(z₂ − z₁)/2 exact, p_h − a exact by Sterbenz) is a
  sum of four multiples of 2⁻¹⁰⁷ below 2⁻⁴⁸, hence of four 64-bit integers below 2⁵⁹ — `_km_mul_div_tie` decides it
  exactly, ties to even. That matters here: midpoints are common (m·3/2 is one whenever 3m has 54 bits — 6.6 % of seeded
  theoretical yields with integer coefficients). RN(F)·2^k is exact for −1021 ≤ k ≤ 1021 (F·2^k > 2⁻¹⁰²², RN(F)·2^k ≤
  2¹⁰²³). A subnormal operand or a k outside that range goes to the exact path `_km_mul_div_acc`: _km_cr_ratio(MₓM_y, M_z,
  Eₓ + E_y − E_z) with M in [2⁵², 2⁵³).
- *The exact comparison*, `_km_ratio_lt(m₁, c₁, m₂, c₂)`: a zero m is the least (two zeros tie); otherwise the products
  m₁·c₂ = P_L·2^E_L and m₂·c₁ = P_R·2^E_R with P in [2¹⁰⁴, 2¹⁰⁶) differ in exponent by two or more (decided), or are
  exact TwoProducts of integer-valued doubles, the larger-exponent side doubled, compared high part first (RN is
  monotone) and then the low parts. `limiting_reagent` keeps Rust's rounded ratio RN(m/c) as the fast path: RN is
  monotone on the extended reals, over- and underflow included, so unequal rounded ratios order the exact ones; only
  equal ones (ties, near-ties, both +∞, both 0) call `_km_ratio_lt`.
- *The empirical formula.* The minimum: Rust's rounded moles, `_km_ratio_lt` on equal ones. Each ratio rᵢ = F·2^s with
  F = (M_pᵢ·M_mⱼ)/(M_mᵢ·M_pⱼ) in (¼, 4) (exact TwoProducts, then `_km_dd_quot`, `_km_dd_div_cr` without its test: to
  13.12u², u = 2⁻⁵³) and s ≥ −1 (rᵢ ≥ 1). **A ratio of 2³² − ½ or more is Err(ComputationError) whichever way Rust's search
  goes**: no multiplier is that Err, and a multiplier k makes the count round(rᵢ·k) ≥ round(rᵢ) ≥ 2³² > u32::MAX, the same
  code; every check after the validation is ComputationError and the port's buffer check comes after them. The port
  returns it as soon as s ≥ 35 (rᵢ > 2³³) or r̃ₕ ≥ 2³² − ¼ (then rᵢ > 2³² − ½); every kept ratio is below 2³², so rᵢ·k <
  2³⁵. The tolerance test: x = rᵢ·k, n = round(k·r̃ₕ), δ = x − n to under 2⁻⁵³; |δ| ≤ 0.1 − 2⁻⁵⁰ proves the pass (and the
  count n), |δ| ≥ 0.1 + 2⁻⁵⁰ the failure (every other integer is ½ − 2⁻¹⁶ away), and between — the knife edge — the exact
  test `_km_emp_close_exact` decides: x = A/B in integers (A < 2¹⁴³, B < 2¹⁰⁸), q = ⌊A/B⌋ by a 36-bit restoring division,
  R = A − qB, the pass 10R ≤ B or 10R ≥ 9B, the count q + [2R ≥ B].

**Rust's empirical Err paths under exact evaluation** (`validate.py` tabulates every corpus row by Rust's path):
- *"mole calculation produced non-positive value"* (the moles of a subnormal percentage round to 0) is gone. Such an input
  is Ok when its exact ratios pass (9 corpus rows: He 2⁻¹⁰⁷⁴ % with He 2·2⁻¹⁰⁷⁴ % is 1 : 2, Og 2⁻¹⁰⁷⁴ % with
  Og 3·2⁻¹⁰⁷⁴ % 1 : 3) and the same Err(ComputationError) otherwise (3 rows: C 5·10⁻³²⁴ % with O 50 %, a ratio past
  2³²).
- *"mole ratio is not finite"* (a ratio overflowing) is gone, but never changes the result: a rounded ratio past DBL_MAX
  means an exact one far past 2³², the ComputationError (5 corpus rows, all kept: C 10⁻³²⁰ % with O 50 %).
- *"no multiplier in 1..=6"*: 188 rows become Ok at the ±1/10 edge (C 12.011 % with O 17.5989 %: r − 1 = 1/10 −
  1/(10·9006636304787571), CO); the rest stay.
- *Ok*: 56 rows become Err(ComputationError) at the edge (C 12.011 % with O 20.7987 %: 3r = 3.9 − 1.0·10⁻¹⁶, so 1/10 and a
  little from 4, where Rust's rounded 3r was within 0.1: 3 : 4), and 17 change
  their multiplier at the edge (C 50 % with C 59 %: 5r = 5.9 exactly is within, so 5 : 6 where Rust's rounded 5.9 failed
  and k = 6 gave 6 : 7); 4,932 are unchanged.
- *"count exceeds u32"*: unchanged (21 rows). The validation Errs are Rust's, in Rust's order.
Of the 3,383 rows M4f generated, 14 change (8 Err → Ok, 6 Ok → Err), the decision's count; the 270 noted rows include
M4f-2's knife-edge, subnormal and 2³²-edge rows, built to sit on the edge. No input changes but those at the ±1/10 edge,
those whose moles or ratios Rust over- or underflows, and those whose Rust ratios are off by enough to cross it.

**Golden comparison** (as #21–#23). The corpus keeps the oracle's values; where the port's result differs the note names
it: `[ADR 0002 #24: the port returns the correctly rounded 100·a/t = <x> (f:…)]` (`m·p/c at the exact limiting reagent`),
`[ADR 0002 #24: the port returns the correctly rounded mass·count·100/mw per entry … (v:n i:z f:… …)]`, `[ADR 0002 #25: the
port returns the exact limiting reagent (mᵢ·c_j < m_j·cᵢ, the first on a tie) = <k> (i:<k>)]` and `[ADR 0002 #26: with
exact ratios and an exact 1/10 tolerance the port returns … (v:n i:z i:c …)]` or `… Err(ComputationError)
(E:ComputationError)]`, computed by the generator in `Big` integers (`exact::cr_mul_div`, `ratio_lt`, `empirical_exact`:
no 2³² cut, the ratios at any size). The harness maps the noted rows (accepted under ADR 2) and counts them per function:
1,126 `percent_yield`, 240 `theoretical_yield_moles`, 598 `percent_composition` (#24: 1,964), 920 `limiting_reagent` (#25),
270 `empirical_formula` (#26). `validate.py` recomputes the port's result of every row of the five with `Fraction` and
proves each note sits exactly where it differs from the oracle's, naming it (so every unnoted row is the port's result
too), checks the 2³² claim and the kept ratios' range on the corpus, tabulates the Err paths above, counts the hard rows
(65 `percent_yield` ties, 17 subnormal; 57 rows at its overflow end and 252 where Rust underflows to 0; 34 and 26 yields
whose m/c over- or underflows, 366 midpoints; 569 exact limiting ties, 1,849 rounded ties, 407 both overflowing, 404 both
underflowing; 504 (row, k) exactly 1/10 off) and pins the exact 1/10 (1/10 − 2⁻⁸⁰ changes 300 rows, the nearest larger
tolerance a row needs, 1/10 + 2.8·10⁻¹⁸, changes one). The grids gain 6,341 rows (`m_stoichiometry::emit_exact`, appended
with their own seed): the ties at RN (subnormal ones, `percent_yield(3·2⁻¹⁰⁷⁴, 200)` among them), the over- and underflow
edges, `limiting_reagent`'s exact ties and near-ties, `empirical_formula`'s knife edges at ±1/10, its subnormal rows and the
2³² edge. `cr_verify.py` checks `muldiv` (640,000 triples), `pyield`, `tyield`, `pcomp`, `limit` and `empir` (2.23 million
arguments in all) against `Fraction` on x86_64 ELF, x86_64 PE and aarch64, the exact paths forced.

**Cost** (`tests/kimiya.bcyr`, before → after, the same host): `stoichiometry/percent_yield` 14 → 81 ns,
`theoretical_yield_moles` 46 → 115, `percent_composition_glucose` 72 → 276 (three entries), `empirical_formula_glucose`
217 → 586, `limiting_reagent` 32 → 33; the new exact-path benches `theoretical_yield_moles_tie` 133 ns (a midpoint,
resolved in the fast path), `percent_yield_subnormal` 2.15 µs (`_km_cr_ratio`), `limiting_reagent_tie` 94 ns,
`empirical_formula_edge` 7.46 µs (two exact tolerance tests).

## Addition, 2026-11-10 (audit B1): #27 and #28

**The user's decisions, verbatim.** The questions were the main session's; the answers are the user's.
1. **The audit (2026-11-09).** Q: "Intermediate products of inputs can overflow or underflow even when the true result is
   representable (…). The same pattern exists in modules ported earlier. How far should the edge-only fix you approved
   for the geometry functions go?" A: **"Library-wide audit (Recommended)"** — "After M4, audit every ported module for
   intermediate over/underflow and fix each edge (in-range rows keep Rust's bits), one ADR 0002 row per function
   family."
2. **Round once in range (2026-11-10).** Q: "… For the audit's formulas that are a single product/quotient of inputs
   (n·R·T/P, ε·l·c, the Faraday products), which?" A: **"Round once (Recommended)"** — "Correctly rounded for every
   input, like #22 and #24: about a quarter of normal-range results move by 1–2 ulp, and each such call gets roughly 5×
   slower (e.g. ~15 → ~80 ns). Formulas with sums inside get the edge fix only."
3. **Settled under the user's policy:** "where Rust gives NaN but the true value is out of range, the port returns the
   correctly rounded ±∞ or 0. Edges that are physically meaningless (T = 10⁻³⁰⁰ K) are fixed too." No such NaN occurs in
   B1's functions: their NaNs from finite inputs are 0·∞ with an exact zero factor, whose value 0 is representable, and
   Michaelis–Menten's ∞/∞ of a value below Vmax; so in B1 the point shows only as the exact 0 or the finite value where
   Rust gives NaN (of the 785 golden rows where Rust's value is 0, ±∞ or NaN and the port's is representable).

The survey behind the audit (sixteen families over every ported module, with the replicas, grids and probes that found
them) proposed the order of the batches; B1 is its families F1 (#27) and F4 (#28), one row each.

**The treatment, per function** (the decisions' rule: one product or quotient rounds once, a sum inside gets #20's edge
fix):

| Function | Rust's formula | Treatment | Why |
|---|---|---|---|
| `thermo::heat_transfer` | (m·c)·ΔT | rounded once, RN(m·c·ΔT) | one product of inputs; no validation, so a non-finite input keeps Rust's expression (IEEE's ±∞, NaN) |
| `thermo::specific_heat_from_calorimetry` | q/(m·ΔT) | rounded once, RN(q/(m·ΔT)) | one quotient |
| `phase::boiling_point_elevation`, `freezing_point_depression` | (K·m)·i | rounded once, RN(K·m·i) | one product; non-finite inputs as `heat_transfer` |
| `phase::osmotic_pressure` | (((i·M)·R)·1000)·T | rounded once, RN(i·M·R·1000·T) | one product; R enters as its double, 1000 exact |
| `kinetics::zero_order_half_life` | [A]₀/(2k) | rounded once, RN([A]₀/(2k)) | one quotient; Rust's division is already that value wherever 2k is finite (2k exact, one IEEE division), so it stays the fast path and the monomial serves k ≥ 2¹⁰²³ |
| `kinetics::steady_state_intermediate` | (k₁·[A])/k₂ | rounded once, RN(k₁·[A]/k₂) | one quotient of a product |
| `kinetics::pre_equilibrium_rate` | (k₂·(k₁/k₋₁))·[A] | rounded once, RN(k₂·k₁·[A]/k₋₁) | one quotient of a product (Rust's k₁/k₋₁ is a quotient inside a product, not a sum) |
| `kinetics::second_order_concentration` | 1/(1/[A]₀ + k·t) | #20's semantics | a sum inside |
| `kinetics::michaelis_menten` | (Vmax·[S])/(Km + [S]) | #20's semantics | a sum inside |

The survey's F4 lists "steady state": that is `steady_state_intermediate`; `steady_state_rate` (k·[X], one operation on
inputs) is in its "no finding" list and is unchanged.

**#27 and the monomials of #28, the machinery** (`src/mathx.cyr`, "correctly rounded monomials"). M4f-2's RN(x·y/z),
`_km_cr_mul_div`, written in `stoichiometry.cyr`, moved to `src/mathx.cyr` and became the 2/1 shape of
`_km_cr_mono(num, nn, den, nd)` = RN(Π num / Π den), for at most six numerator and four denominator doubles, all finite
and the denominators nonzero (otherwise NaN), the exact path's capacity on each side 253 bits of odd significand parts
(four arbitrary doubles and a 41-bit integer; osmotic_pressure's five factors use 219): the sign is the xor of the
factors' signs and a zero factor gives that signed zero (IEEE's rule for products and quotients, so Rust's ±0 wherever
Rust's chain gives one, and the exact 0 where Rust's 0·∞ is NaN). With each factor s·2^e, the value is F·2^k, F = N/D
over the significands other than 1:
- *the fast path* (`_km_cr_mono_fast`): N and D as double-doubles — exact for one or two factors (a TwoProduct), each
  further factor by `_km_dd_mul_d` (proven within 3.01u², u = 2⁻⁵³), at most 18.1u² in all, inside `_km_dd_div_cr`'s
  32u² budget — and its rounding test (`_km_dd_div_ends`). When the test fails, F lies within 2⁻⁹⁸·F of the one midpoint
  μ of two adjacent doubles z₁ < z₂, and for three numerator and two denominator factors or fewer the sign of N − D·μ, an
  exact sum of at most ten doubles (N's exact expansion of TwoProducts, D·z₁'s TwoProducts, D·(u/2) exact), decides it
  (`_km_mono_tie`, by Shewchuk's Grow-Expansion with zero elimination, `_km_exp_sign`): ties go to even. Midpoints are
  common with small integer factors (3·y with 3y of 54 bits), so they stay in the fast path, as M4f-2's
  `_km_mul_div_tie` kept them. The result y = RN(F) scales exactly where y·2^k ≥ 2⁻¹⁰²¹ (the binary64 grid is the 53-bit
  one there), is ±∞ past 2¹⁰²³ (IEEE overflow is the 53-bit rounding past DBL_MAX) and ±0 below 2⁻¹⁰⁷⁶;
- *the exact path* (`_km_cr_mono_acc`): the odd integer significands multiplied in `_km_bn_*`, `_km_cr_ratio` rounding
  P/Q·2^E once (ties to even, the subnormal grid, +∞ past DBL_MAX), the sign after; it serves possibly subnormal results,
  the rare failed test beyond the tie shapes, and the tests' forced exact path.
Through its factor tables and loops the general fast path is slow in this compiler (`heat_transfer` took 148 ns a call
through it, 65 through its own straight-line path), so each shape the
modules call has a straight-line fast path of its own for normal operands, the same arithmetic without the tables:
`_km_mul_div_fast` (M4f-2's, with its 64-bit integer midpoint resolution `_km_mul_div_tie`, now for both signs),
`_km_mul3_fast` and `_km_mul5_fast` (no denominator: the product's double-double rounded by `_km_dd_round_ends`, ε =
2⁻¹⁰⁰·h against a proven 9.1u² at most), `_km_div2_fast` and `_km_mul3_div_fast` (`_km_dd_div_ends`), all scaled by
`_km_mono_fin`; a zero, subnormal or non-finite operand, a failed test (a midpoint other than the 2/1 shape's, or within
2⁻⁹⁸ of one) or a possibly subnormal result falls back to `_km_cr_mono`. `stoichiometry` calls the 2/1 wrapper
`_km_cr_mul_div` with the same results (the same correctly rounded value on every argument; `cr_verify.py`'s `muldiv`,
`pyield`, `tyield` and `pcomp` families pass unchanged); the other wrappers are `_km_cr_mul3`, `_km_cr_div2`,
`_km_cr_mul3_div` and `_km_cr_mul5`.

**#28's unbounded-exponent functions** (`src/kinetics.cyr`), on #20's helpers. `_km_second_order_fast` returns Rust's
1/(1/[A]₀ + k·t) when 1/[A]₀ and a nonzero k·t lie in [2⁻¹⁰²¹, DBL_MAX] (a result that large is the 53-bit rounding of
an exact value in the normal range; k·t = 0 only from a zero factor) and the sum is finite (a sum of two non-negative
doubles at least 2⁻¹⁰²¹); `_km_michaelis_menten_fast` when Vmax·[S] lies in [2⁻¹⁰²¹, DBL_MAX] or [S] = ±0 and Km + [S] is
finite (a sum of non-negative doubles below 2⁻¹⁰²¹ is exact). The final division is then one IEEE division of the same
operands as the unbounded chain's, rounded once in both, so every such input keeps Rust's bits. Otherwise
`_km_second_order_ux` / `_km_michaelis_menten_ux` evaluate Rust's chain with `_km_ux_from`, `_km_ux_div`, `_km_ux_mul`,
`_km_ux_add` and the final `_km_ux_div_f64`. Under these semantics an intermediate keeps its 53-bit rounding:
`second_order_concentration(DBL_MAX, 0, 0)` is 1/RN₅₃(1/DBL_MAX) = 1/(2⁻¹⁰²⁴(1 + 2⁻⁵²)), DBL_MAX's predecessor (Rust's ∞),
and `michaelis_menten(100, DBL_MAX, DBL_MAX)` is exactly 50 (Rust's NaN).

**Golden comparison** (as #16–#26). The corpus keeps the oracle's values; where the port's value differs, the note names
it: `[ADR 0002 #27: the port returns the correctly rounded m·c·ΔT = <x> (f:…)]` (likewise `q/(m·ΔT)`, `K·m·i`,
`i·M·R·1000·T`), `[ADR 0002 #28: the port returns the correctly rounded [A]₀/(2k) = <x> (f:…)]` (likewise `k₁·[A]/k₂`,
`k₂·k₁·[A]/k₋₁`) and `[ADR 0002 #28: the port returns Rust's sequence evaluated with an unbounded exponent (the final
division rounded once) = <x> (f:…)]`. The generator computes them independently of the port (`scripts/oracle/src/audit.rs`:
`exact::cr_monomial` in `Big` integers, `exact::Ux` with its new `div` and `div_f64`), with Rust's validation replicated
and asserted against the oracle on every row, and appends 2,129 edge rows with their own seeds (thermo 537, phase 528,
kinetics 1,064: partial products past DBL_MAX or below 2⁻¹⁰²² with the other factor bringing them back, the true overflow
and underflow ends, 0 against an overflowing product, 53-bit midpoints at scales across the range and the subnormal grid's
ties, k ≥ 2¹⁰²³, 1/[A]₀ and k·t past either end, ∞/∞). The harness maps a noted row (`_par_got_noted`, and
`_par_got_noted_f` for the three f64-valued functions) and counts it per function: **682 #27 rows** (151 `heat_transfer`,
171 `specific_heat_from_calorimetry`, 82 `boiling_point_elevation`, 82 `freezing_point_depression`, 196
`osmotic_pressure`) and **636 #28 rows** (89 `second_order_concentration`, 90 `zero_order_half_life`, 118
`steady_state_intermediate`, 179 `pre_equilibrium_rate`, 160 `michaelis_menten`). Of the rows that existed before B1,
88 of 359, 90 of 374, 33 + 33 of 186 each, 74 of 207, 4 of 987, 2 of 262, 52 of 264, 69 of 283 and 5 of 370 change: the
monomials' in-range quarter, and the edges. `validate.py` §19 recomputes every row of the ten functions — the monomials
with `Fraction`, the two others through `mpfr_ref.UX` — and proves that a row is noted exactly where the port's value
differs from the oracle's, the note naming it bit for bit (so every unnoted row is the port's value too), that no row whose
f64 chain stays in the normal range is noted for the two #20-style functions (1,003 and 328 rows: Rust's bits), and that
`zero_order_half_life` is noted only where 2k overflows; it counts the hard rows (785 where Rust's value is 0, ±∞ or NaN
and the port's representable, 323 53-bit midpoints, 551 subnormal results with 133 ties). `cr_verify.py` checks `mono`,
`monotie` and the ten functions' public, fast and forced slow paths (`heat`, `shc`, `boil`, `frz`, `osm`, `zhl`, `ssi`,
`pre`, `soc`, `mm`) against `Fraction` and `mpfr_ref.UX` on x86_64 ELF, x86_64 PE and aarch64.

**Cost** (`tests/kimiya.bcyr`, before → after, the same host; `docs/benchmarks-rust-v-cyrius.md`'s B1 section):
`thermo/heat_transfer` 5.4 → 65 ns (Rust 0.95), `specific_heat_from_calorimetry` 18 → 86, the colligative pair 6 → 66–67,
`osmotic_pressure` 20 → 139, `steady_state_intermediate` 21 → 88, `pre_equilibrium_rate` 27 → 125, `zero_order_half_life`
16 → 16, `second_order_concentration` 23 → 33, `michaelis_menten` 17 → 24 (Rust 2.92); a midpoint `heat_transfer_tie`
365 ns, a subnormal `osmotic_pressure_subnormal` 4.04 µs, the unbounded-exponent edges 230 and 126 ns; stoichiometry
as in M4f-2 (`percent_yield` 82, `theoretical_yield_moles` 115, `percent_composition_glucose` 266 ns).

## Addition, 2026-11-10 (audit B2): #29 and #30

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B2 unchanged: a formula that
is one product or quotient of inputs and constants is rounded once, correctly rounded for every input; a formula with a
sum inside gets #20's edge fix (Rust's chain with an unbounded exponent, the final operation rounded once, a fast-path
range test so in-range rows keep Rust's bits); where Rust gives NaN and the true value is out of range, the correctly
rounded ±∞ or 0; meaningless edges are fixed too; every Rust check, Err code and order stays. As in B1, no B2 function
has a NaN from finite inputs whose true value is out of range: their NaNs are ∞·0 with an exact zero factor
(`absorbance(1e300, 1e300, 0)`; `nernst_potential` at Q = 1 with R·T past DBL_MAX, whose value is E°) and 0/0 with an
exact zero numerator (`concentration_from_absorbance` with A = 0 and ε·l underflowing, `faraday_mass_ratio` with m_A = 0
and M_A/n_A underflowing), each a representable value. B2 is the survey's families F3 (#29) and F5 (#30), one row each.

**The treatment, per function** (each checked against `rust-old/src/` and the survey):

| Function | Rust's formula | Treatment | Why |
|---|---|---|---|
| `solution::dilution` | (M₁·V₁)/V₂ | rounded once, RN(M₁·V₁/V₂) | one quotient of a product of inputs |
| `solution::solution_conductivity` | (Λ·c)/1000 | rounded once, RN(Λ·c/1000) | one quotient; 1000 is an integer literal of Rust's formula (exact); no validation, so a non-finite input keeps Rust's expression |
| `electrochemistry::faraday_mass_deposited` | (Q·M)/(n·F) | rounded once, RN(Q·M/(n·F)) | one quotient; F enters as its double, n exact |
| `electrochemistry::faraday_charge_required` | ((m·n)·F)/M | rounded once, RN(m·n·F/M) | one quotient |
| `electrochemistry::faraday_mass_ratio` | (m_A·(M_B/n_B))/(M_A/n_A) | rounded once, RN(m_A·M_B·n_A/(n_B·M_A)) | Rust's two equivalent weights are quotients inside one quotient (as #28's `pre_equilibrium_rate`'s k₁/k₋₁), not a sum |
| `electrochemistry::nernst_potential` | E° − ((R·T)/(n·F))·ln Q | #20's semantics, the correctly rounded ln Q as the primitive | a sum (the subtraction) inside |
| `spectroscopy::absorbance` | (ε·l)·c | rounded once, RN(ε·l·c) | one product; no validation, so a non-finite input keeps Rust's expression |
| `spectroscopy::concentration_from_absorbance` | A/(ε·l) | rounded once, RN(A/(ε·l)) | one quotient; A is not validated, so a non-finite A keeps Rust's expression (B1's rule for unvalidated inputs) |
| `spectroscopy::rigid_rotor_energy` | ((h·B)·J)·(J + 1) | rounded once, RN(h·B·J·(J + 1)) | one product; h as its double, J and J + 1 exact (u32: J + 1 ≤ 2³²); a non-finite B keeps Rust's expression |
| `spectroscopy::rotational_constant` | h/(((8·π)·π)·I) | rounded once, RN(h/(8·π·π·I)) | one quotient; π enters as its double each time Rust uses it, 8 exact (computed as (h/8)/(π·π·I), h/8 exact) |
| `spectroscopy::diatomic_moment_of_inertia` | (μ·r)·r | rounded once, RN(μ·r·r) | one product; non-finite inputs keep Rust's expression |
| `spectroscopy::vibrational_frequency` | (k/μ).sqrt() | the correctly rounded √ of the exact k/μ | a monomial under a square root: √ rounded once, as the decision rounds the monomials (the planned treatment; its fallback, #20's semantics with `_km_ux_sqrt_f64`, was not needed) |

The survey's F3 lists `nernst_potential_25c` as "no finding": it calls `nernst_potential` at T = 298.15 K, whose
(R·T)/(n·F) is a normal constant-like factor, so every input takes the fast path and keeps Rust's bits (`validate.py`
proves each golden row's chain normal); `cell_gibbs_energy` and the charge conversions are unchanged.

**#29's and #30's monomials, the machinery** (`src/mathx.cyr`). #27's `_km_cr_mono` serves them, with one straight-line
fast path per new shape, B1's arithmetic without the factor tables (a zero, subnormal or non-finite operand, a failed
test or a possibly subnormal result falls back to `_km_cr_mono`):
- `_km_mul2_div2_fast` (RN(x·y/(z·w)): N and D exact TwoProducts, `_km_dd_div_ends`) — `faraday_mass_deposited`;
- `_km_mul3_div2_fast` (RN(x·y·z/(v·w)): N by `_km_dd_mul_d` within 3.01u², D exact) — `faraday_mass_ratio`;
- `_km_mul4_fast` (RN(a·b·c·d): within 6.1u², `_km_dd_round_ends`) — `rigid_rotor_energy`;
- `_km_div3_fast` (RN(x/(y·z·w)): D within 3.01u²) — `rotational_constant`;
with their wrappers `_km_cr_mul2_div2`, `_km_cr_mul3_div2`, `_km_cr_mul4` and `_km_cr_div3`; `dilution`,
`solution_conductivity` (z = 1000), `faraday_charge_required`, `absorbance`, `concentration_from_absorbance` and
`diatomic_moment_of_inertia` use B1's `_km_cr_mul_div`, `_km_cr_mul3_div`, `_km_cr_mul3` and `_km_cr_div2`. Midpoints
of the 2/2 and 3/2 shapes (3F·y/F with 3F exact, F's odd part having 51 bits; y·w·3/(n_B·w)·n_B) fall back to
`_km_cr_mono`'s fast path, which resolves them exactly (`_km_mono_tie`); the 4/0 and 1/3 shapes have none in practice
(h's odd part has 53 bits, π's squared 100), and their rare failed tests take the exact path.

**#30's square root, the machinery** (`src/mathx.cyr`, "correctly rounded square roots of quotients"):
`_km_cr_sqrt_div(x, y)` = RN(√(x/y)) for positive finite doubles, rounded once (the subnormal grid, +∞ past DBL_MAX):
- *no midpoint*: with x/y = (a/b)·2^E (a, b odd, below 2⁵³) a rounding boundary m = M·2^f (M odd) has √(x/y) = m only if
  a = b·M² and E = 2f, impossible: a 53-bit-grid midpoint has M of 54 bits (M² > a), a subnormal-grid one f = −1075
  (E ≥ −2097). So the fast path's test fails only near a midpoint, and the exact path never breaks a tie;
- *the fast path* `_km_sqrt_div_fast`, normal x and y: the significands with K = eₓ − e_y made even (sₓ doubled), F =
  sₓ/s_y in (½, 4) as q_h + q_l (an exact remainder), y₀ = RN(√q_h), the residual q_h − y₀² by TwoProduct, d = (t₁ +
  q_l)/(2y₀): |y₀ + d − √F| ≤ 5.7u²·√F (the bound in the code: 6.01u² for t, 4.51u² for d, 1.14u² for the Taylor
  term); then `_km_dd_div_cr`'s test (ε = 2⁻¹⁰⁰·y₀, the two ends rounded alike) and `_km_mono_fin`'s exact scale by
  2^(K/2) (normal x and y give at most 2¹⁰²³·√2); it fails with probability about 2⁻⁴⁶, and declines possibly
  subnormal results;
- *the exact path* `_km_sqrt_div_acc`, every positive finite x and y: x = A·2^eₐ, y = B·2^e_b (A, B in [2⁵², 2⁵³)), the
  exponent made even, R = ⌊√(A·2¹¹⁰/B)⌋ (in [2⁵⁴, 2⁵⁶)) by a guess from the fast path's arithmetic corrected by exact
  multi-word comparisons (R + 1 while (R + 1)²·B ≤ A·2¹¹⁰, R − 1 while R²·B > A·2¹¹⁰), sticky = [R²·B ≠ A·2¹¹⁰], and
  `_km_cr_from_int(2R + sticky, h − 1)` — `_km_cr_ratio`'s argument: 2R + sticky has 56 or 57 bits and every rounding
  boundary is a multiple of 4 of its units.
The generator's reference is `exact::cr_sqrt_ratio` (an integer square root of a 121-bit quotient, `round_scaled`),
`validate.py`'s Python's `math.isqrt` and, independently, MPFR's √ of a 256-bit quotient (`mpfr_ref.ux_sqrt_div_final`:
the no-midpoint argument gives every √(x/y) a distance of at least 2⁻¹⁶² from every boundary, far beyond the quotient's
2⁻²⁵⁶).

**#29's Nernst equation** (`src/electrochemistry.cyr`). `_km_nernst_fast` returns Rust's E° − ((R·T)/(n·F))·l (l =
`_km_cr_ln(Q)`) when R·T lies in [2⁻¹⁰²², DBL_MAX], (R·T)/(n·F) is at least 2⁻¹⁰²² and the product with l is at least
2⁻¹⁰²² in magnitude or l = 0 (n·F is always normal, the product always finite): each operation is then the unbounded
chain's, and the final subtraction of two doubles rounds once in both. Otherwise `_km_nernst_ux` evaluates the chain
with `_km_ux_from`, `_km_ux_mul`, `_km_ux_div` and the final `_km_ux_sub_f64`. So
`nernst_potential(0.342, 2, DBL_MAX, 0.1)` is 1.78·10³⁰⁴ (Rust's ∞), `nernst_potential(0.342, 2, DBL_MAX, 1)` E° (Rust's
∞·0 = NaN) and `nernst_potential(0, 1, 10⁻³¹⁰, 10⁻³⁰⁰)` 5.952642933123·10⁻³¹² (Rust's degraded …4057·10⁻³¹²).

**Golden comparison** (as #27 and #28). The corpus keeps the oracle's values; where the port's value differs the note
names it: `[ADR 0002 #29: the port returns the correctly rounded M₁·V₁/V₂ = <x> (f:…)]` (likewise `Λ·c/1000`,
`Q·M/(n·F)`, `m·n·F/M`, `m_A·M_B·n_A/(n_B·M_A)`), `[ADR 0002 #29: the port returns Rust's sequence evaluated with an
unbounded exponent (the final subtraction rounded once) = <x> (f:…)]`, and `[ADR 0002 #30: the port returns the
correctly rounded ε·l·c = <x> (f:…)]` (likewise `A/(ε·l)`, `√(k/μ)`, `h·B·J·(J + 1)`, `h/(8·π·π·I)`, `μ·r·r`). The
generator computes them independently of the port (`scripts/oracle/src/audit.rs`: `exact::cr_monomial`, the new
`exact::cr_sqrt_ratio`, `exact::Ux` with its new `add_f64` / `sub_f64`; the Nernst ln is glibc's, correctly rounded on
every golden row — `validate.py` checks it against MPFR), with Rust's validation replicated and asserted against the
oracle on every row, and appends 2,459 edge rows with their own seeds (solution 458, electrochemistry 823, spectroscopy
1,178): partial products past DBL_MAX or below 2⁻¹⁰²² with the other factors bringing them back, Λ·c in [2¹⁰²⁴, 2¹⁰³⁴),
equivalent weights M/n subnormal or 0, R·T past DBL_MAX and subnormal (R·T)/(n·F), E° against an overflowing product,
k/μ past DBL_MAX or below 2⁻¹⁰⁷⁴, subnormal k, μ and results, exact roots, twelve roots within 2⁻¹⁰⁷ of a 53-bit
midpoint (continued-fraction convergents of m²), h·B and 8π²·I subnormal, the 53-bit midpoints 3·y/2^j, 3F·y/F, y·3·F/F,
y·w·3/(n_B·w)·n_B, 3000·y/1000 and 3·r·r (3r² of 54 bits), the subnormal grid's ties. The harness maps a noted row
(`_par_got_noted`, `_par_got_noted_f`) and counts it per function: **2,756 #29 rows** (94 `dilution`, 148
`solution_conductivity`, 904 `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 119
`nernst_potential`) and **2,106 #30 rows** (195 `absorbance`, 219 `concentration_from_absorbance`, 706
`vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant`, 118 `diatomic_moment_of_inertia`). Of the
rows that existed before B2, 38 of 208, 50 of 186, 870 of 3,975, 565 of 2,168, 716 of 2,199, 2 of 1,485, 131 of 560, 130
of 564, 440 of 3,047, 273 of 1,015, 459 of 1,587 and 111 of 440 change (the same order): the monomials' in-range quarter
(an eighth for the √: Rust's √ of RN(k/μ) misses one input in eight), and the edges. `validate.py` §20 recomputes every
row of the twelve functions — the monomials with `Fraction`, the √ with `math.isqrt` and MPFR, Nernst through
`mpfr_ref.UX` with MPFR's log and `ux_sub_final` — and proves that a row is noted exactly where the port's value differs
from the oracle's, the note naming it bit for bit, that no `nernst_potential` row whose f64 chain stays in the normal
range is noted (1,455 rows: Rust's bits) and that every `nernst_potential_25c` row's chain is normal, that
no √ row is a rounding boundary, and counts the hard rows: 818 where Rust's value is 0, ±∞ or NaN and the port's
representable (56 `dilution`, 64 `solution_conductivity`, 9 `faraday_mass_deposited`, 28 `faraday_charge_required`, 71
`faraday_mass_ratio`, 70 `nernst_potential`, 76 `absorbance`, 79 `concentration_from_absorbance`, 336
`vibrational_frequency`, 29 `rigid_rotor_energy`), 409 midpoints of the 53-bit grid, 1,309 subnormal results with 60 of
the subnormal grid's ties, and for the √ 86 exact roots, 195 quotients past DBL_MAX, 146 below 2⁻¹⁰⁷⁴, 168 subnormal
operands, 3 subnormal results and 13 roots within 2⁻⁹⁸ of a 53-bit midpoint; no row where Rust's NaN stands for an
out-of-range value (none is possible: the decisions paragraph above). `cr_verify.py` checks `sqrtdiv` and the twelve
functions' public, fast and forced slow paths (`dil`, `cond`, `fmass`, `fcharge`, `fratio`, `nernst`, `abs`, `conc`,
`vib`, `rotor`, `rot`, `moment`) against `Fraction`, exact integer roots and `mpfr_ref.UX` on x86_64 ELF, x86_64 PE and
aarch64, roots within 2⁻¹⁰⁰ of a midpoint and the shapes' midpoints included.

**Cost** (`tests/kimiya.bcyr`, before → after, the same host; `docs/benchmarks-rust-v-cyrius.md`'s B2 section):
`spectroscopy/absorbance` 5 → 68 ns (Rust 0.98), `electrochemistry/nernst_potential` 83 → 92 ns (Rust 6.42; its range
test), `dilution` 17 → 86, `solution_conductivity` 6 → 72, `faraday_mass_deposited` 19 → 109, `faraday_charge_required`
21 → 122, `faraday_mass_ratio` 33 → 148, `concentration_from_absorbance` 16 → 90, `vibrational_frequency` 18 → 87,
`rigid_rotor_energy` 14 → 114, `rotational_constant` 18 → 117, `diatomic_moment_of_inertia` 11 → 74; the edges
`nernst_potential_edge` 202 ns (the unbounded slow path), `vibrational_frequency_edge` 88 ns (k/μ past DBL_MAX: the fast
path), `vibrational_frequency_subnormal` 883 ns (the exact root), `faraday_mass_ratio_tie` 402 ns (`_km_mono_tie`).

## Addition, 2026-11-11 (audit B3): #31 and #32

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B3 unchanged: a formula that is
one product or quotient of inputs and constants, possibly under a final square root, is rounded once (B2's
`vibrational_frequency` = RN(√(k/μ)) the precedent for the root); a formula with a sum inside gets #20's edge fix (Rust's
chain with an unbounded exponent, the final operation rounded once, a fast-path range test so in-range rows keep Rust's
bits); where Rust gives NaN and the true value is out of range, the correctly rounded ±∞ or 0; meaningless edges are
fixed too (Tc = 2⁻¹⁰²² K, ω = 10²⁰⁰, atomic masses of 10³⁰⁸ u); every Rust check, Err code and order stays. Unlike B1's and
B2's, B3's functions have such NaNs: `lennard_jones`' ∞ − ∞ past σ/r = 2^170.7 is +∞ where 4ε·d overflows (and a finite
value for a tiny ε), and the equations of state's 0·∞ or ∞ − ∞ stand for values that are finite or ±∞ (232 golden rows
where Rust's NaN is the value ±∞). The unbounded evaluation also gives the correctly signed zero where the value
underflows: `lennard_jones` with σ/r below 2⁻¹⁷⁹ and a huge ε is −0 (or a negative subnormal), Rust's d = 0 − 0 = +0
(137 rows where Rust's special differs and the value is out of range). B3 is the survey's families F2 (#31) and F6
(#32), one row each; the survey's batches B4 and B5.

**The treatment, per function** (each checked against `rust-old/src/`):

| Function | Rust's formula | Treatment | Why |
|---|---|---|---|
| `gas::ideal_gas_pressure`, `ideal_gas_volume` | ((n·R)·T)/V, ((n·R)·T)/P | rounded once, RN(n·R·T/x) | one quotient of a product of inputs and R's double |
| `gas::gas_density` | (P·M)/(R·T) | rounded once, RN(P·M/(R·T)) | one quotient; R·T is a product in the denominator, not a sum |
| `gas::van_der_waals_pressure` | (((n·R)·T)/(V − n·b)) − ((a·(n/V))·(n/V)), after V ≤ n·b | #20's semantics, the co-volume test on the unbounded n·b | sums inside (V − n·b, the final difference) |
| `gas::peng_robinson_pressure` | RT/(V − b) − a·α/((V·V + (2b)·V) − b·b) on `pr_coefficients` (κ = (0.37464 + 1.54226ω) − (0.26992ω)ω, α = (1 + κ(1 − √(T/Tc)))², `powi(2)` = u·u), after V ≤ b | #20's semantics, the whole chain (κ and α included), the co-volume test on the unbounded b | sums inside |
| `gas::redlich_kwong_pressure` | RT/(V − b) − a/((√T·V)·(V + b)), a = (((0.42748·R)·R)·Tc^2.5)/Pc, after V ≤ b | #20's semantics; Tc^2.5 with an unbounded exponent by an exact rescale | sums inside; `powf(2.5)` is the chain's one transcendental (correctly rounded on both sides, ADR 0006) |
| `potential::lennard_jones` | (4·ε)·((σ/r)¹² − (σ/r)⁶), (σ/r)⁶ five products | #20's semantics, #17 folded in | a difference inside |
| `potential::coulomb`, `coulomb_force` | ((k·q₁)·q₂)/r, ((k·q₁)·q₂)/(r·r) | rounded once, RN(k·q₁·q₂/r), RN(k·q₁·q₂/(r·r)) | one quotient each (r·r a product in the denominator); the charges are not validated, so a non-finite one keeps Rust's expression (B1's rule) |
| `potential::morse_width_parameter` | ω·(μ/(2·D_e)).sqrt() | the correctly rounded ω·√(μ/(2·D_e)) of the exact value, #18 folded in | a monomial with a square root: ω·√(μ/(2·D_e)) = √(ω²·μ/(2·D_e)) for ω > 0, rounded once as #30's √(k/μ) |
| `inorganic::born_lande_lattice_energy` | −k·N_A·M·z⁺·z⁻·e·e/r₀·(1 − 1/n), left to right | #20's semantics | (1 − 1/n) is a sum; M4b left it as Rust's (it misses CR((n − 1)/n) on 332 rows, but no exact value exists to reach), and #20's semantics keep exactly that: in-range rows keep Rust's bits |
| `nuclear::binding_energy_per_nucleon` | ((Z·m_H + N·m_n − M)·931.494…)/A | #20's semantics | the mass defect is a sum of inputs (M is one), so not a monomial |

The survey listed these under F2 and F6 with the same treatments; `binding_energy_per_nucleon` was the open question (its
mass defect, Z·m_H + N·m_n − M, is a sum with an input). The Peng–Robinson Z, roots and fugacity share
`pr_coefficients` and keep Rust's chain (`_km_pr_coefficients`, unchanged; the survey's F16, the audit's later batches).

**#31's and #32's monomials** use #27's `_km_cr_mono` through B1's and B2's straight-line shapes (`_km_cr_mul3_div` for
the ideal gas law and `coulomb`, `_km_cr_mul2_div2` for `gas_density`, `_km_cr_mul3_div2` for `coulomb_force`); no new
shape. R's odd significand part has 51 bits and k_e's 52 (B1's comment on `osmotic_pressure` says 53 for R; its
conclusion, no midpoint there, holds, since 1000's 125 adds 7 bits), so q·R for q = 9 … 15 and 5·k_e, 7·k_e are 53-bit
midpoints, rows the corpus carries (ties to even, the shapes' failed tests resolved by `_km_mono_tie`).

**#32's square root, the machinery** (`src/mathx.cyr`, "correctly rounded w·√((x/y)·2^e)"):
`_km_cr_mul_sqrt_ratio(w, x, y, e)` = RN(w·√((x/y)·2^e)) for positive finite doubles and |e| ≤ 2²⁰, rounded once (the
subnormal grid, +∞ past DBL_MAX); `morse_width_parameter` calls it with e = −1.
- *midpoints exist*, unlike #30's √(x/y): with w = W·2^a, x = A·2^b, y = B·2^c (W, A, B odd) a boundary M·2^f (M odd) is
  hit when W²·A = B·M² and the exponents match, e.g. ω = y with 3y of 54 bits, μ = 9, D_e = ½ (v = 3y); so the fast path
  only ever proves a rounding, and the exact path decides every tie;
- *the fast path* `_km_mul_sqrt_ratio_fast`, normal operands: #30's `_km_sqrt_div_dd` for √(sₓ/s_y) (within 5.7u²,
  u = 2⁻⁵³, |d| ≤ 1.52u·y₀, proven there), s_w·y₀ = p_h + p_l by TwoProduct, t = RN(p_l + RN(d·s_w)) (within 4.07u²·p_h
  of s_w(y₀ + d)), so p_h + t is within 9.8u² of s_w·√F, normalised by FastTwoSum; then `_km_dd_round_ends`' test (ε =
  2⁻¹⁰⁰·h, its 32u² budget) and `_km_mono_fin`'s exact scale (±∞ past 2¹⁰²³; a possibly subnormal result declines);
- *the exact path* `_km_mul_sqrt_ratio_acc`, every positive finite operand: v = √(P/B)·2^h with P = W²·A·2⁶ (W, A, B in
  [2⁵², 2⁵³), A doubled for an odd exponent; P below 2¹⁶⁶), R = ⌊√(P/B)⌋ in [2^54.5, 2⁵⁷) from a guess off the fast
  path's arithmetic (8·W·(y₀ + d)) corrected by multi-word comparisons (`_km_bn_*`), sticky = [R²·B ≠ P], and
  `_km_cr_from_int(2R + sticky, h − 1)`: R has at least 55 bits, so every rounding boundary is a multiple of 4 units of
  2R + sticky, an end point of [2R, 2R + 2), reached only by an exact value (sticky = 0): ties to even.
The generator's reference is `exact::cr_sqrt_monomial` (an integer square root of a 121-bit quotient of the mantissas'
products, `round_scaled`), `validate.py`'s Python `math.isqrt` (`b3_sqrt_frac`) and, independently, MPFR's √ of a
320-bit quotient (`mpfr_ref.ux_mul_sqrt_div_final`: exact where the quotient is dyadic, so a tie is rounded exactly;
otherwise at least 2⁻¹⁶² from every boundary's square).

**#31's and #32's unbounded-exponent functions** (`src/gas.cyr`, `src/potential.cyr`, `src/inorganic.cyr`,
`src/nuclear.cyr`), on #20's helpers. Each runs Rust's chain and a range test; only where the test fails, the
unbounded-exponent evaluation (`_km_vdw_ux`, `_km_pr_pressure_ux`, `_km_rk_pressure_ux`, `_km_lj_ux`,
`_km_born_lande_ux`, `_km_bea_ux`). The tests rest on `src/mathx.cyr`'s new `_km_ux_normal(x)`, |x| in [2⁻¹⁰²¹, DBL_MAX]:
an IEEE product, quotient or root that large is the 53-bit rounding of an exact value above 2⁻¹⁰²² (2⁻¹⁰²² itself is
excluded: a value just below it rounds there at the subnormal grid and below it at 53 bits); a sum of two correct
doubles needs only to be finite (a subnormal sum is exact); a zero product of an exact zero factor is exact; √ of a
correct operand is correct. Then every operation is the unbounded chain's, the co-volume test decides alike, and the
final IEEE operation rounds the same operands once.
- `_km_vdw_fast`: n·b normal, an exact zero or +∞ (V ≤ +∞ refuses in both: the unbounded n·b is past DBL_MAX ≥ V); then
  n·R, n·R·T, n/V and the quotient term normal (n ≠ 0), a·(n/V) and its product with n/V normal (a ≠ 0), V − n·b finite.
- `_km_pr_pressure_fast`: every product and quotient of `pr_coefficients` and the pressure (the a, b and α chains, T/Tc,
  R·T, V·V, 2b·V, b·b, the two quotients) normal or an exact zero (ω = 0, κ·(1 − √T_r) with a zero factor, u = 0), κ and
  the sums finite.
- `_km_rk_pressure_fast`: Tc^2.5 (`_km_cr_pow`, RN₅₃(Tc^2.5) wherever it is normal), a's and b's products and quotients,
  R·T, the quotient term, √T·V and (√T·V)·(V + b) normal; V − b and V + b finite.
- `_km_lj_fast`: σ/r and (σ/r)⁶ normal (the four products between lie between them: RN is monotone), and (σ/r)¹² normal
  or (σ/r)⁶ below 2⁻⁵¹² (then (σ/r)¹² is below 2⁻¹⁰²⁴, far under half an ulp of (σ/r)⁶ in both evaluations, so d = −(σ/r)⁶
  in both); the final (4·ε)·d, and for ε ≥ 2¹⁰²² #17's (ε·d)·4 (ε·d is never subnormal for a nonzero d, ×4 exact), each
  RN(4·ε·d).
- `_km_born_lande_fast`: the five partial products and the quotient by r₀ normal (−k·N_A is a normal constant); the
  factor 1 − 1/n is always the unbounded one (where Rust's 1/n is not normal, n > 2¹⁰²¹, both 1/n are below 2⁻⁵⁴ and both
  differences round to 1).
- `binding_energy_per_nucleon`: Δm·931.494… finite. Δm itself is always the unbounded one (Z·m_H and N·m_n are 0 or in
  [1, 2³³·1.01], Δm = s − M with s ≥ 1 is 0 or at least 2⁻⁵³ in magnitude), so only the product's overflow (masses above
  about 1.93·10³⁰⁵ u) leaves the range.
- **Redlich–Kwong's Tc^2.5** (`_km_rk_tc25_ux`). With Tc = s·2^(2j), s in [1, 4) (Tc's unbounded pair, the significand
  doubled for an odd exponent: exact), Tc^2.5 = s^2.5·2^(5j) exactly (real arithmetic, s > 0). s^2.5 lies in [1, 32),
  inside binary64's normal range, where the correctly rounded `_km_cr_pow(s, 2.5)` is RN₅₃(s^2.5) (ADR 0006: exact cases
  and midpoints decided before its accurate path rounds, the general case within its stated 2⁻¹³⁰). With an unbounded
  exponent RN₅₃ commutes with scaling by a power of two: RN₅₃(x·2^k) = RN₅₃(x)·2^k for every x > 0 and integer k (53-bit
  rounding acts on the significand only). So RN₅₃(Tc^2.5) = RN₅₃(s^2.5)·2^(5j) = `_km_cr_pow(s, 2.5)`·2^(5j), the pair
  (its significand, its exponent + 5j), which is the unbounded chain's value. `validate.py` checks the identity on every
  golden row with MPFR on both sides (MPFR's pow of Tc at 53 bits and the widest exponent range against MPFR's
  binary64 pow of s scaled by 2^(5j)), and `cr_verify.py`'s `rkp` family checks the slow path against MPFR's own unbounded
  pow, its ties (Tc = t²·4^j with t⁵ of 54 bits, past either end) included.
So `van_der_waals_pressure(10³⁰⁶, 300, 10³⁰⁰, 0, 10⁻¹⁰)` is 2.49·10⁹ (Rust's +∞), `peng_robinson_pressure` with
Tc = 2⁻¹⁰²² is 99,773.55 Pa (Rust's NaN) and with ω = 10²⁰⁰ at T = Tc 41,775 Pa (Rust's −∞·0 = NaN),
`lennard_jones(2⁻¹⁰⁷⁴, 2¹⁷¹, 1)` is 2⁹⁸⁰ (Rust's NaN), `born_lande_lattice_energy(2⁻¹⁰⁷⁴, 1, 1, 2.82 Å, 8)` −2.13·10⁻³¹⁸
(Rust's −0), and a co-volume that Rust rounds at the subnormal grid decides the other way where the exact one does
(`van_der_waals_pressure(3·2⁻¹⁰⁷⁴, 300, 2⁻¹⁰⁷⁴, 0, 0.2)` is Ok, 18,707.5 Pa, where Rust's n·b rounds up to V: Err; a
subnormal Tc over a tiny Pc scales the rounding error of (0.07780·R)·Tc back up, so V between Rust's b and the exact one
is the port's Err where Rust accepts, or the reverse).

**#17 and #18, folded into #32.** #20's semantics subsume M4d's two reorders (the survey's note): the exact 4·ε and 2·D_e
cannot overflow with an unbounded exponent. #17's rows keep their bits: on every `lennard_jones` row with ε ≥ 2¹⁰²² whose
d chain is normal, (ε·d)·4 is the unbounded chain's RN(4·ε·d) (`validate.py` proves it on all 366 such rows); its 180 M4d
rows are noted under #32 with the same values, and the fast path keeps (ε·d)·4 as its own branch. #18 is not
bit-identical: the Morse width now rounds the whole ω·√(μ/(2·D_e)) once, everywhere (the round-once rule; M4d's
(μ·½)/D_e only rounded the quotient correctly). Of #18's 181 rows (D_e ≥ 2¹⁰²³), 128 keep their value, 42 move by one ulp
and 3 by two (Rust's √ and ×ω rounding a correctly rounded normal quotient), and 8 by up to 1.3·10¹² ulp (a subnormal
μ/(2·D_e), which #18 still rounded at the subnormal grid before the √); 159 more rows with D_e ≥ 2¹⁰²³, where #18 gave
Rust's 0 (μ/(2·D_e) below 2⁻¹⁰⁷⁵), now carry the representable value (`morse_width_parameter(1, 2⁻¹⁰⁷⁴, 2¹⁰²³)` is
2⁻¹⁰⁴⁹). Rows #17 and #18 above record the fold; the harness counts #32's rows in the two regions on their own (244 and
340, B3's edge rows included).

**Golden comparison** (as #27 … #30). The corpus keeps the oracle's values; where the port's result differs the note names
it: `[ADR 0002 #31: the port returns the correctly rounded n·R·T/V = <x> (f:…)]` (likewise `n·R·T/P`, `P·M/(R·T)`), `[ADR
0002 #31: the port returns Rust's sequence evaluated with an unbounded exponent (the final subtraction rounded once) = <x>
(f:…)]` or `… = Err(InvalidInput) (E:InvalidInput)]` (the port's co-volume test refuses), `[ADR 0002 #32: the port returns
Rust's sequence evaluated with an unbounded exponent (the final product rounded once) = <x> (f:…)]` (`lennard_jones`,
`born_lande_lattice_energy`; `binding_energy_per_nucleon` the final division), `[ADR 0002 #32: the port returns the
correctly rounded k_e·q₁·q₂/r = <x> (f:…)]` (likewise `k_e·q₁·q₂/(r·r)`, `ω·√(μ/(2·D_e))`). The generator computes them
independently of the port (`scripts/oracle/src/audit.rs`: `exact::cr_monomial`, the new `exact::cr_sqrt_monomial`,
`exact::Ux` with the new `mul_f64` and `pow25` — Tc^2.5 as RN₅₃(√(Tc⁵)) by an integer square root, not the port's
rescale), with Rust's validation replicated and asserted against the oracle on every row (only the co-volume tests may
differ), and appends 3,804 edge rows with their own seeds (gas 1,680, potential 1,703, inorganic 254, nuclear 167):
partial products past DBL_MAX or below 2⁻¹⁰²² with the other factors bringing them back, the true ends, the 53-bit
midpoints q·R and q·k_e and the subnormal grid's ties, n/V, V − n·b, V·V, Tc·Tc, ω², R·T past DBL_MAX, subnormal Tc,
κ = −∞ at T = Tc, Tc^2.5 past either end and its ties there, the co-volume on its ulps at the subnormal grid and where a
subnormal Tc's error is scaled back up, σ/r past 2^85.3 and 2^170.7 and below 2^-179, the final product's subnormal ties,
the Morse width's quotient past either end, its midpoints 3ω·2^k and subnormal ties and B2's near-midpoint quotients,
Born–Landé's prefix past DBL_MAX and its subnormal M, n past 2¹⁰²¹, Δm·931.494… past DBL_MAX with A ≥ 932 (and below,
where the value itself overflows). Tc^2.5's midpoints in range are not in the corpus: glibc's pow, the oracle's,
misrounds them, which makes them ADR 0006 rows, not B3's (`validate.py` checks that glibc's pow(Tc, 2.5) is correctly
rounded on every Redlich–Kwong row, so the notes, made against glibc's oracle, are against the corrected one too).
The harness maps a noted row (`_par_got_noted`, and the new `_par_got_noted_r` for the three equations of state: the
note's value, or its Err) and counts it per function: **1,157 #31 rows** (146 `ideal_gas_pressure`, 148
`ideal_gas_volume`, 222 `gas_density`, 206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234
`redlich_kwong_pressure`; 98 of them the port's Err where Rust's rounded co-volume accepts, 80 the port's value where
Rust's refuses) and **4,182 #32 rows** (582 `lennard_jones`, 244 with ε ≥ 2¹⁰²²; 1,860 `morse_width_parameter`, 340 with
D_e ≥ 2¹⁰²³; 730 `coulomb`; 679 `coulomb_force`; 175 `born_lande_lattice_energy`; 156 `binding_energy_per_nucleon`).
Of the rows that existed before B3, 73 of 256, 87 of 255, 57 of 205, 1 of 297, 3 of 660, 0 of 623, 100 of 2,105 (280
noted, 180 of them #17's with the same value), 1,552 of 4,509 (1,680 noted, 128 of them #18's with the same value), 519 of
1,930, 437 of 1,930, 23 of 3,369 and 0 of 618 change (the same order): the monomials' in-range quarter to third (the
Morse width's three roundings miss on about a third), and the edges; the #20-style functions only at their edges.
`validate.py` §21 recomputes every row of the twelve functions — the monomials with `Fraction`, the Morse width with
`math.isqrt` and MPFR, the six others through `mpfr_ref.UX` with MPFR's own pow for Tc^2.5 (the new `UX.pow`), the
co-volume tests on UX values and `ux_sub_final`, `ux_mul_final`, `ux_div_final` — and proves that a row is noted exactly
where the port's result differs from the oracle's, the note naming it bit for bit (or its Err), that no row whose f64
chain stays in the normal range is noted for the six #20-style functions (285, 570, 540, 1,480, 3,174 and 544 rows:
Rust's bits; its tracers count 2⁻¹⁰²² itself as out of range), the #17 fold and the #18 counts above, the Tc^2.5 rescale
on the 878 Redlich–Kwong rows it reaches, and counts the hard rows: 2,085 where Rust's value is 0, ±∞ or NaN and
the port's representable (52, 44, 103, 57, 38, 126, 349, 608, 137, 299, 116, 156 in the order above), 232 where Rust's NaN
is the value ±∞, 137 where Rust's other special differs from the out-of-range value (a −0 for Rust's +0), 178 co-volume
decisions the other way (98 + 80), 215 midpoints of the 53-bit grid, 891 subnormal results with 144 of the subnormal grid's
ties, 272 Tc^2.5 past either end; for the √ 18 exact roots, 122 midpoints (ties), 153 quotients μ/(2·D_e) past DBL_MAX,
446 below 2⁻¹⁰⁷⁴, 169 subnormal operands, 300 subnormal results and 18 roots within 2⁻⁹⁸ of a 53-bit midpoint.
`cr_verify.py` checks `wsqrt` (the helper at e in [−4, 4] and far) and the twelve functions' public, fast and forced slow
paths (`igp`, `igv`, `dens`, `vdw`, `prp`, `rkp`, `lj`, `width`, `coul`, `force`, `lande`, `bea`) against `Fraction`, exact
integer roots and `mpfr_ref.UX` on x86_64 ELF, x86_64 PE and aarch64 (2,815,142 arguments), the midpoints and the
co-volume flips included.

**Cost** (`tests/kimiya.bcyr`, before → after, the same host; `docs/benchmarks-rust-v-cyrius.md`'s B3 section):
`gas/ideal_gas_pressure` 18 → 119 ns (Rust 2.18), `nuclear/binding_energy_per_nucleon_fe56` 28 → 29 ns (Rust 2.02: the
inline test), `potential/lennard_jones` 36 → 39 ns; Cyrius-only `ideal_gas_volume` 20 → 119, `gas_density` 20 → 108,
`van_der_waals_pressure` 35 → 55, `peng_robinson_pressure` 78 → 100, `redlich_kwong_pressure` 182 → 187, `coulomb`
16 → 126, `coulomb_force` 18 → 145, `morse_width_parameter` 33 → 135, `born_lande_lattice_energy` 46 → 50 ns; the slow
paths `van_der_waals_pressure_edge` 260 ns, `peng_robinson_pressure_edge` 618 ns, `redlich_kwong_pressure_edge` 584 ns,
`lennard_jones_edge` 208 ns, `born_lande_lattice_energy_edge` 287 ns, `binding_energy_per_nucleon_edge` 117 ns, and the
Morse width's exact root, `morse_width_parameter_subnormal` 1.27 µs and `morse_width_parameter_tie` 1.27 µs. The #20-style
range tests first ask whether the inputs lie in a region where every intermediate is provably normal (`_km_pr_box` and
its inline siblings: Tc, Pc, T, V within 2^±60, every physical state), which halved their cost (Peng–Robinson 147 → 100,
van der Waals 73 → 55, Redlich–Kwong 217 → 187, Born–Landé 67 → 50 ns); inlining `_km_ux_normal`'s bit test was measured
(147 → 156 ns) and not kept.

## Addition, 2026-11-11 (audit B4): #33 and #34

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B4 unchanged. Every B4
function has a sum inside, so each gets #20's edge fix: Rust's operation sequence — its summation order included — with
an unbounded exponent, the final operation rounded once, behind a fast-path test so in-range rows keep Rust's bits;
where Rust gives NaN and the value is out of range (∞ − ∞ of two overflowing terms whose difference overflows too), the
correctly rounded ±∞; meaningless edges are fixed too (sums near DBL_MAX, Shomate at T = 10⁻³⁰⁰ K); every Rust check,
Err code and order stays. F7's ∞ − ∞ can also stand for a finite value (two sums past DBL_MAX that cancel), which the
unbounded evaluation gives. B4 is the survey's families F7 (#33) and F8 (#34), its proposed B6. It also records two
things outside them: the correction of #29's fast path, and the bench split. The latter is the user's decision,
verbatim: Q: "`cyrius bench` stops any bench file that runs longer than its default 300-second deadline.
tests/kimiya.bcyr now takes about 279 s … How should the benches fit?" A: **"Split the bench file (Recommended)"** —
"Move the audit's edge and slow-path benches into a second file, e.g. tests/audit.bcyr. Each file stays well under 300
s, and the deadline keeps catching a hung bench."

**The treatment, per function** (each checked against `rust-old/src/`):

| Function | Rust's formula | Treatment | Why |
|---|---|---|---|
| `thermo::hess_law` | Σ xᵢ, `.sum()` from −0.0 | #20's semantics, the last addition rounded once; fast test "result finite" | a pure sum |
| `thermo::enthalpy_from_formation` | Σ p − Σ r | #20's semantics, the difference rounded once; fast test both sums finite | a difference of two sums |
| `thermochem::reaction_enthalpy`, `reaction_gibbs_energy` | (Σ c·ΔH_f − Σ c·ΔH_f)·1000, `total += c·prop` from +0.0 (#1's × 1000 the port's last operation) | #20's semantics, the product with 1000 rounded once; fast test: each product normal or an exact zero, both sums and the difference finite | sums of products |
| `thermochem::reaction_entropy` | Σ c·S° − Σ c·S° | as above, the difference rounded once | sums of products |
| `nuclear::q_value` | (Σ m_r − Σ m_p)·931.494… | #20's semantics, the product rounded once; fast test the difference finite | a difference of sums |
| `solution::ionic_strength` | ½·Σ (c·z)·z from −0.0 | #20's semantics, ½·Σ rounded once; fast test the sum finite (a term of a finite c and an integer z is exact below 2⁻¹⁰²¹, a multiple of 2⁻¹⁰⁷⁴, so it leaves the range only by overflow) | a sum of products |
| `inorganic::born_haber_lattice_energy` | ΔH_f − ΔH_sub − IE − ΔH_diss/2 − EA, left to right | #20's semantics, the last difference rounded once; fast test the result finite and ΔH_diss/2 exact | a chain of differences (the survey's "a chain of differences") |
| `reaction::gibbs_free_energy` | ΔH − T·ΔS | #20's semantics, the difference rounded once; fast test T·ΔS normal or ΔS = 0 | a difference |
| `solution::kohlrausch_molar_conductivity` | Λ° − K·√c | as `gibbs_free_energy` (√c of a non-negative double is 0 or normal and correctly rounded in both) | a difference |
| `solution::limiting_molar_conductivity` | ν₊·λ₊ + ν₋·λ₋ | #20's semantics, the sum rounded once; fast test both products normal or exact zeros | a sum of products |
| `kinetics::zero_order_concentration` | [A]₀ − k·t, 0 if negative | #20's semantics, the clamp on the unbounded difference, the difference rounded once | a difference; the comparison is #20's "on the unbounded values" |
| `Shomate::cp` | (((A + B·t) + (C·t)·t) + ((D·t)·t)·t) + E/(t·t), t = T/1000 | #20's semantics on Rust's evaluation order, the last addition rounded once | a polynomial: sums of products (#34) |
| `Shomate::enthalpy_antiderivative` | 1000·(t·(A + t·(B/2 + t·(C/3 + t·(D/4)))) − E/t) | as `cp`, the product by 1000 rounded once | a nested polynomial (#34) |

The unvalidated functions (`hess_law`, `enthalpy_from_formation`, the reaction sums' coefficients, `q_value`,
`ionic_strength`, `born_haber_lattice_energy`, `limiting_molar_conductivity`, `zero_order_concentration`, the Shomate
polynomials) keep Rust's expression for a non-finite input (B1's rule); a slice that is no slice keeps its NaN, a charge
outside i32 its NaN, an unknown formula Rust's Err in Rust's order.

**The pure-sum lemma, and the fast tests.** For doubles x and y, x + y is exact wherever |x + y| < 2⁻¹⁰²¹: both are
integer multiples of 2⁻¹⁰⁷⁴, so x + y is one, and a multiple k·2⁻¹⁰⁷⁴ with |k| < 2⁵³ is a double (Hauser). Elsewhere
IEEE's x + y is the 53-bit rounding of the exact sum, and ±∞ exactly where that rounding passes DBL_MAX. So in a chain
of additions of finite doubles each partial sum is the unbounded chain's until one is ±∞, and from there on the chain
stays ±∞ or NaN (a finite addend leaves ±∞ alone; ∞ − ∞ is NaN). A finite result therefore proves every partial sum the
unbounded one, and the final addition is rounded once in both: the fast test of a pure sum is "result finite"; only
overflow can be spurious, which is the survey's observation. A product inside a sum takes B3's `_km_ux_normal` (|p| in
[2⁻¹⁰²¹, DBL_MAX]: the 53-bit rounding of an exact value in the normal range) or an exact zero (a zero factor), and an
integer-valued factor of magnitude ≥ 1 keeps a product of doubles a multiple of 2⁻¹⁰⁷⁴ (`ionic_strength`'s c·z·z: exact
below 2⁻¹⁰²¹, so the sum's test covers it); `born_haber_lattice_energy`'s halving rounds only for an odd multiple of
2⁻¹⁰⁷⁴ below 2⁻¹⁰²¹, tested as (d/2)·2 = d. `validate.py` checks the lemma partial sum by partial sum on every pure-sum
golden row (786 rows with finite sums, 226 with a sum past DBL_MAX).

**#33's machinery** (`src/thermo.cyr`, `src/thermochem.cyr`, `src/nuclear.cyr`, `src/solution.cyr`, `src/inorganic.cyr`,
`src/reaction.cyr`, `src/kinetics.cyr`) is #20's helpers: Rust's `.sum()` from −0.0 as `_km_sum_ux` (and
`_km_sum_ux_f64`, the last addition rounded once), thermochem's fold from +0.0 as `_km_sum_property_ux`, the slow paths
`_km_formation_ux`, `_km_q_value_ux`, `_km_reaction_j_ux` (the difference to 53 bits, × 1000 rounded once),
`_km_reaction_s_ux`, `_km_ionic_strength_ux`, `_km_born_haber_ux`, `_km_gibbs_ux`, `_km_kohlrausch_ux`,
`_km_limiting_ux`, `_km_zero_order_ux` on `_km_ux_add`, `_km_ux_sub`, `_km_ux_mul`, `_km_ux_div`, `_km_ux_sqrt` and the
final `_km_ux_add_f64`, `_km_ux_sub_f64`, `_km_ux_mul_f64`; the fast tests inline, or `_km_kohlrausch_fast` and
`_km_limiting_fast`, and `_km_sum_property` reports per list whether every product agreed and every coefficient is
finite. So `hess_law([DBL_MAX, DBL_MAX, −DBL_MAX])` is DBL_MAX (Rust's +∞), `reaction_enthalpy([Al₂O₃ ×10³⁰⁸] →
[Al₂O₃ ×10³⁰⁸])` +0 (Rust's NaN), `q_value([M, M], [M, 0.999M])` 1.67·10³⁰⁸ (Rust's NaN), `ionic_strength([(0.4·DBL_MAX,
2)])` 1.44·10³⁰⁸ (Rust's +∞), `gibbs_free_energy(DBL_MAX, 2, 0.75·DBL_MAX)` −8.99·10³⁰⁷ (Rust's −∞),
`zero_order_concentration(−DBL_MAX, −DBL_MAX, 1.5)` 8.99·10³⁰⁷ (Rust's +∞). Under the semantics an intermediate keeps
its 53-bit rounding and a zero its sign: `born_haber_lattice_energy(0, 0, 0, 2⁻¹⁰⁷⁴, 0)` is −0 (−2⁻¹⁰⁷⁵ rounded once;
Rust's half is +0), and `zero_order_concentration(−0, 10⁻²⁰⁰, 10⁻²⁰⁰)` +0 (the unbounded difference −10⁻⁴⁰⁰ is negative,
so the clamp gives 0, where Rust's underflowed k·t leaves −0 − 0 = −0).

**#34's machinery** (`src/thermochem.cyr`). `_km_shomate_cp5` and `_km_shomate_h5`, the one evaluation behind
`KmShomate.cp` / `enthalpy_antiderivative` and the table-row path, compute Rust's chain and test first a box:
|T| in [2⁻³⁰, 2³⁰] and every nonzero coefficient in [2⁻²⁰⁰, 2²⁰⁰] (`_km_shomate_box`). There t = T/1000 is in [2⁻⁴⁰, 2²¹]
and every product and quotient is normal or an exact zero: in Cp, B·t, C·t·t, D·t·t·t, t·t and E/(t·t) lie in
[2⁻³²⁰, 2²⁸¹]; in the antiderivative D/4, C/3, B/2 lie in [2⁻²⁰², 2¹⁹⁹], and each nested sum is 0 or at least an ulp of
its smaller term (a multiple of it), so t·(D/4) ≥ 2⁻²⁴², C/3 + t·D/4 is 0 or ≥ 2⁻²⁹⁴, × t ≥ 2⁻³³⁴, B/2 + … 0 or ≥ 2⁻³⁸⁶,
× t ≥ 2⁻⁴²⁶, A + … 0 or ≥ 2⁻⁴⁷⁸, × t ≥ 2⁻⁵¹⁸, E/t ≥ 2⁻²²¹, every value below 2²⁹⁷, each sum finite. Outside the box
`_km_shomate_cp_fast` / `_km_shomate_h_fast` test each product and quotient (normal or an exact zero) and each sum
(finite); T = ±0 (every product an exact zero, E/0 IEEE's special in both) and a non-finite input keep Rust's chain;
otherwise `_km_shomate_cp_ux` / `_km_shomate_h_ux` evaluate Rust's order with an unbounded exponent. So
`Shomate{20.786, 0, 0, 0, 0}.cp(10⁻³⁰⁰)` is 20.786 (Rust's 0/0 = NaN), its antiderivative at 5·10⁻³²⁴ K 1.04·10⁻³²²
(Rust's NaN), H₂O(g)'s Cp at 10³⁰⁰ K −∞ (Rust's ∞ − ∞ = NaN), and at t = 2⁶⁰⁰ with C = 1, D = −2⁻⁶⁰⁰ the value is +0:
Rust's order with 53 bits absorbs A into C·t² before D·t³ cancels it (Rust's NaN) — the semantics are Rust's sequence,
not the exact polynomial. **The data-backed callers** share the evaluation: every SHOMATE_DATA coefficient is 0 or in
[2⁻²⁰⁰, 2²⁰⁰] and every range end in [2⁻³⁰, 2³⁰] (`validate.py` checks the generated table, `tests/thermochem.tcyr`
every row), so the table-row path (`_km_shomate_row_cp`, `_km_shomate_row_h`) tests T alone, inline, and takes the
general evaluation only outside [2⁻³⁰, 2³⁰]; `enthalpy_change_cp`, which evaluates only at temperatures within the
ranges, always gets Rust's bits; the flame's Newton iterates may leave the ranges, and `validate.py` replays every
golden flame row's 324,925 evaluations: all inside the box. Their golden rows are unchanged.

**#29 corrected.** B2's `_km_nernst_fast` accepted an intermediate of at least 2⁻¹⁰²² (`<` tests against 2⁻¹⁰²²). An
IEEE result of exactly 2⁻¹⁰²² can stand for an exact value just below it, which binary64 rounds up at the subnormal grid
but 53-bit rounding does not: for the product c·l = ((R·T)/(n·F))·ln Q in [2⁻¹⁰²² − 2⁻¹⁰⁷⁵, 2⁻¹⁰²² − 2⁻¹⁰⁷⁶) IEEE gives
2⁻¹⁰²² and #20's chain 2⁻¹⁰²² − 2⁻¹⁰⁷⁵, and E° − c·l then differs where E° is an odd multiple of 2⁻¹⁰⁷⁴ (the unbounded
difference is a subnormal tie, rounded to even): `nernst_potential(5·10⁻³²⁴, 1, 5.1625987955063826·10⁻³⁰³,
1.0512872158898858)` is −0x0.ffffffffffffep-1022, where B2 returned Rust's −0x0.fffffffffffffp-1022. The test now uses
`_km_ux_normal`'s bound, 2⁻¹⁰²¹ (inline: R·T and c are positive), as B3's tests do; only values strictly above 2⁻¹⁰²²
are safe, and the slow path gives the same value on [2⁻¹⁰²², 2⁻¹⁰²¹). The factor tests could not misfire in practice —
R·T at 2⁻¹⁰²² makes the factor subnormal, which the next test refuses, and a quotient of two doubles lies within 2⁻⁵³
(relative) below a power of two only when the divisor's significand is 1, never for n·F — but they use the same bound.
`validate.py`'s shared tracer `_t_chk` (§16's #20, §19's B1 and §20's B2 tracers) took 2⁻¹⁰²² as in range too; it now
flags a product or quotient of exactly ±2⁻¹⁰²² that is not its exact value's 53-bit rounding (no count moves). The 23
rows of the generator's `NERNST_2M1022` (the three the main session found, eleven from a seeded search, the first three
with E° = 3, 2, 0 units of 2⁻¹⁰⁷⁴) are golden edge rows, 17 of them noted; `tests/electrochemistry.tcyr`'s
`cyrius_adr29_2m1022` runs them (it fails on B2's code), and `cr_verify.py`'s `nernst` family gains 2,800 such
arguments. **The other fast paths of B1–B3**, audited for the same boundary: B1's `_km_second_order_fast` and
`_km_michaelis_menten_fast` and B3's `_km_lj_fast` test against 2⁻¹⁰²¹; B3's `_km_vdw_fast`, `_km_pr_pressure_fast`,
`_km_rk_pressure_fast` and `_km_born_lande_fast` use `_km_ux_normal`; the binding energy's inline test is "finite" (Δm
itself is never subnormal); the monomial and root shapes (`_km_mono_fin`) decline results whose exponent is below −1021;
their operand tests (`nmin` = 2⁻¹⁰²² on the inputs) are on exact doubles, where 2⁻¹⁰²² is normal; #20's geometry tests
use 2⁻⁹⁰⁰ and 2⁻⁴⁵⁸. None has the boundary.

**Golden comparison** (as #27 … #32). The corpus keeps the oracle's values; where the port's value differs the note
names it: `[ADR 0002 #33: the port returns Rust's sequence evaluated with an unbounded exponent (the final addition
rounded once) = <x> (f:…)]` (`hess_law`, `limiting_molar_conductivity`; "the final subtraction" for
`enthalpy_from_formation`, `reaction_entropy`, `born_haber_lattice_energy`, `gibbs_free_energy`,
`kohlrausch_molar_conductivity`, `zero_order_concentration`; "the final product" for `reaction_enthalpy`,
`reaction_gibbs_energy`, `q_value`, `ionic_strength`), after #1's `[J/mol: …]` for the two J/mol sums, and #34's
likewise (`cp` the final addition, `enthalpy_antiderivative` the final product). The generator computes them
independently of the port (`scripts/oracle/src/audit.rs`: `exact::Ux` in Rust's summation orders, Rust's expression for
a non-finite input and for the Shomate polynomials at T = ±0), with Rust's validation replicated and asserted, and
appends 2,507 edge rows with their own seeds (thermo 184, reaction 172, nuclear 88, solution 430, inorganic 170,
kinetics 171, thermochem 1,292): partial sums and differences past ±DBL_MAX that later terms bring back, two sums past
DBL_MAX that cancel exactly or nearly, DBL_MAX + 2⁹⁷⁰ (the rounding past DBL_MAX itself), products past DBL_MAX or
subnormal, an odd subnormal ΔH_diss, k·t below 2⁻¹⁰⁷⁵ against ±0, the reaction sums' coefficients near DBL_MAX/|prop|
and subnormal, the Shomate polynomials with user coefficients at T → 0, large T, t = 2⁶⁰⁰'s cancellation and every
binade. The harness maps a noted row (`_par_got_noted`, `_par_got_noted_f`) and counts it per function: **697 #33 rows**
(28 `hess_law`, 27 `enthalpy_from_formation`, 115 `reaction_enthalpy`, 113 `reaction_gibbs_energy`, 66
`reaction_entropy`, 58 `q_value`, 13 `ionic_strength`, 74 `born_haber_lattice_energy`, 57 `gibbs_free_energy`, 54
`kohlrausch_molar_conductivity`, 31 `limiting_molar_conductivity`, 61 `zero_order_concentration`) and **210 #34 rows**
(185 `cp`, 25 `enthalpy_antiderivative`). Of the rows that existed before B4 only the survey's change: 2 of 429 for each
reaction sum (Al₂O₃ ×10³⁰⁸ both sides), 70 of 1,188 `cp` (10 at T → 0, 60 of Rust's NaN at 10³⁰⁰ K and DBL_MAX, now ±∞)
and 6 of 1,188 `enthalpy_antiderivative` (3 at 5·10⁻³²⁴ K, 3 at 2⁻¹⁰²² K). The survey's 93 NaN-for-∞ Shomate rows are
60: its other 33 are `enthalpy_antiderivative` at T = −0, where its replica dropped the zero's sign; Rust's +sign(E)·∞
there is IEEE's E/(−0), the unbounded evaluation's too, and stays. `validate.py` §22 recomputes every row of the
fourteen functions through `mpfr_ref.UX` with its new `ux_sum` / `ux_sum_final` (Rust's summation orders) and the final
operations, and proves that a row is noted exactly where the port's value differs from the oracle's, the note naming it
bit for bit, that no row whose f64 chain stays in the normal range is noted (B3's strict tracers), the lemma, and the
data-backed callers' box; it counts the noted rows: 522 where Rust's value is 0, ±∞ or NaN and the port's representable,
165 where Rust's NaN is the value ±∞, 1 where Rust's other special differs from an out-of-range value, 5 where only a
zero's sign differs, 34 degraded (the value normal, Rust's more than 2 ulp off), 180 within 2 ulp or subnormal.
`cr_verify.py` checks the twelve families `hess`, `eff`, `rsum`, `qv`, `ionic`, `haber`, `gibbs`, `kohl`, `lim`, `zoc`,
`shcp`, `shh` (the public function, the Rust-chain fast path where its test passes, the slow path forced; 1,250,069
arguments) against `mpfr_ref.UX` on x86_64 ELF, x86_64 PE and aarch64.

**The bench split** (the user's decision above). `tests/kimiya.bcyr` keeps the benches that mirror
`rust-old/benches/benchmarks.rs` and the Cyrius-only ones from before the audit; `tests/audit.bcyr` holds the audit's
Cyrius-only benches, B1's to B4's (`bench_audit_b1` … `bench_audit_b4`), their typical-input benches with their edge,
slow-path and tie benches (moving the edges alone left kimiya.bcyr above 200 s). The names are unchanged; `cyrius bench`
runs both files. B3's seventeen benches are back at the standard 1,000,000 × 20 calls (the two microsecond exact roots
at 100,000 × 20, as the other files' microsecond rows). On this host `cyrius bench tests/kimiya.bcyr` takes 177 s (261 s
before the split; B3 measured 279 s) and `cyrius bench tests/audit.bcyr` 189 s, each well inside the 300 s per-file
deadline.

**Cost** (before → after, the same host; `docs/benchmarks-rust-v-cyrius.md`'s B4 section): `reaction/gibbs_free_energy`
17 → 18 ns (Rust 1.75), `thermochem/reaction_enthalpy_ch4` 259 → 293 ns (Rust 38.58), `thermochem/shomate_cp` 15 → 32 ns
(Rust 2.34: the box test), `thermochem/enthalpy_change_cp` 171 → 171 ns (Rust 24.00) and
`thermochem/adiabatic_flame_ch4_air` 3.40 → 3.29 µs (Rust 301.48 ns: the table-row path tests T alone, inline),
`electrochemistry/nernst_potential` 90 → 89 ns (Rust 6.42; #29's bound); Cyrius-only `hess_law` 13 → 15,
`enthalpy_from_formation` 20 → 24, `reaction_gibbs_energy` 259 → 299, `reaction_entropy` 259 → 286, `q_value` 24 → 26,
`ionic_strength` 24 → 31, `born_haber_lattice_energy` 10 → 14, `kohlrausch_molar_conductivity` 22 → 28,
`limiting_molar_conductivity` 10 → 24, `zero_order_concentration` 12 → 17, `shomate_enthalpy_antiderivative` 37 → 54 ns;
the slow paths 117–587 ns (`hess_law_edge` 162, `reaction_enthalpy_edge` 587, `gibbs_free_energy_edge` 117), the Shomate
polynomials' 311 and 474 ns.

## Addition, 2026-11-11 (audit B5): #35

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B5 unchanged. None of B5's
four functions is a product or quotient of inputs and constants (the round-once decision's scope: "n·R·T/P, ε·l·c, the
Faraday products"): each is one exp or 10^x of an argument, and the activity models' arguments have sums inside. So each
gets #20's edge fix on its argument, as B2's `nernst_potential` did around its ln Q (the survey's remedy for its family
F9): Rust's sequence for the argument with an unbounded exponent, the argument's final operation rounded once into
binary64, then the unchanged correctly rounded exp or 10^x, behind a fast-path test so in-range rows keep Rust's bits;
where Rust gives NaN and the value is representable (∞·0 at I = 0, ∞/∞), that value; meaningless edges are fixed too (T
= DBL_MAX K, A = 10³⁰⁰); every Rust check, Err code and order stays, `equilibrium_constant`'s ComputationError judged on
the port's K. No new primitive. B5 is the survey's family F9, its proposed B7.

**The treatment, per function** (each checked against `rust-old/src/reaction.rs` and `solution.rs`):

| Function | Rust's formula | The port's argument (then `_km_cr_exp` / `_km_cr_exp10`) | Fast test (inline) |
|---|---|---|---|
| `reaction::equilibrium_constant` | exp((−ΔG)/(R·T)), ComputationError where K is not finite | R·T at 53 bits, the quotient rounded once (`_km_keq_arg_ux`); ComputationError on that K | R·T in [2⁻¹⁰²¹, DBL_MAX] |
| `solution::debye_huckel_limiting` | 10^((((−A)·z)·z)·√I) | ((−A)·z)·z and √I at 53 bits, the product rounded once (`_km_dh_limiting_arg_ux`) | A·z² finite |
| `solution::debye_huckel_extended` | 10^(((((−A)·z)·z)·√I)/(1 + (B·å)·√I)) | the numerator and 1 + (B·å)·√I at 53 bits, the quotient rounded once (`_km_dh_extended_arg_ux`) | the numerator and (B·å)·√I finite |
| `solution::davies_activity` | 10^((((−A)·z)·z)·(√I/(1 + √I) − 0.3·I)) | ((−A)·z)·z and the bracket at 53 bits, the product rounded once (`_km_davies_arg_ux`) | A·z² finite |

Every input of the four is validated (ΔG, T, I, A, B and å finite; z an i32), so no non-finite input reaches the
arithmetic. Signed zeros follow #20 (an exact zero factor gives IEEE's signed zero in both evaluations; the sign of a
zero argument does not reach K or γ).

**The lemma.** Let v be the exact result of the argument's last operation on its 53-bit operands, X = RN₅₃(v) the
unbounded chain's argument (#20's target is the correctly rounded e^X or 10^X: the transcendental is its final
operation) and x = RN(v) in binary64, the port's (`_km_ux_div_f64` / `_km_ux_mul_f64`, or Rust's IEEE operation where
the fast test passes). Then RN(e^x) = RN(e^X) and RN(10^x) = RN(10^X):
- |X| in [2⁻¹⁰²¹, DBL_MAX]: |v| exceeds 2⁻¹⁰²² (RN₅₃ is monotone), where binary64's rounding is the 53-bit one, and does
  not overflow, so x = X;
- |X| > DBL_MAX, so |X| ≥ 2¹⁰²⁴: IEEE overflow is "the 53-bit rounding passes DBL_MAX", so x = ±∞ with X's sign, and
  e^X, 10^X lie past either end of binary64 (2¹⁰²⁴ ≫ 709.79 and 308.26; −2¹⁰²⁴ ≪ −745.14 and −323.61): +∞ or +0, which
  `_km_cr_exp(±∞)` and `_km_cr_exp10(±∞)` return;
- |X| < 2⁻¹⁰²¹ (zero included): |v| < 2⁻¹⁰²¹ and |x| ≤ 2⁻¹⁰²¹; for |t| ≤ 2⁻¹⁰²¹, e^t and 10^t = e^(t·ln 10) lie within
  2⁻¹⁰¹⁹ of 1, far inside 1's rounding interval (1 − 2⁻⁵⁴, 1 + 2⁻⁵³), so both round to 1 (`_km_cr_exp` returns 1 for
  every |t| < 2⁻⁵⁴, `_km_cr_exp10` for every |t| < 2⁻⁵⁶).

So the port needs no exp of an unbounded argument, and the argument's final rounding has no tie that matters: in the
normal range RN and RN₅₃ break a tie alike (to even), and below 2⁻¹⁰²¹ the result is 1 whatever the argument. (A
quotient of two 53-bit values is never a 53-bit midpoint, since the dividend's odd part would need 54 bits; products can
be: `debye_huckel_limiting(1, 9, 2.666…)` with 3·2.666… = (2⁵³ + 1)·2⁻⁵⁰ has the argument −8 in both, the tie rounded to
even, and γ = 10⁻⁸; golden rows.) `validate.py` checks the lemma on every row, and `cr_verify.py` on every argument, by
computing e^X and 10^X themselves with MPFR.

**The fast tests, proven.**
- `equilibrium_constant`: R·T in [2⁻¹⁰²¹, DBL_MAX] (R·T > 0, since T ≥ 2⁻¹⁰⁷⁴). R·T is then the 53-bit rounding of its
  exact value (B3's `_km_ux_normal` argument, inline), −ΔG is exact, and Rust's IEEE quotient of those operands is x
  itself. 2⁻¹⁰²² is excluded, as B4 corrected #29: T = 541,658,533,363,376·2⁻¹⁰⁷⁴ has R·T = (2⁵² − 0.27…)·2⁻¹⁰⁷⁴, which
  IEEE rounds up to 2⁻¹⁰²² and 53-bit rounding to 2⁻¹⁰²² − 2⁻¹⁰⁷⁵; there K(−2⁻¹⁰²²) is e^(1 + 2⁻⁵²) = 2.718281828459046
  (0x4005bf0a8b14576b), two ulp from Rust's e, which a 2⁻¹⁰²² bound would have kept (golden rows,
  `tests/reaction.tcyr`).
- `debye_huckel_limiting` and `davies_activity`: A·z² = ((−A)·z)·z finite. B4's `ionic_strength` argument: a product of
  a double and an integer-valued double is a multiple of 2⁻¹⁰⁷⁴, so it is exact below 2⁻¹⁰²¹ (k·2⁻¹⁰⁷⁴ with |k| < 2⁵³ is
  a double) and the 53-bit rounding above unless it overflows; an overflowed (−A)·z stays ±∞ (z ≠ 0), so a finite A·z²
  proves both products the unbounded ones. √I of a non-negative double is ±0 or at least 2⁻⁵³⁷, normal, correctly
  rounded in both. Davies' bracket is always the unbounded one: 1 + √I is at least 1 and below 2⁵¹²; √I/(1 + √I) is 0 or
  at least min(√I, 1)/2 ≥ 2⁻⁵³⁸, normal; 0.3·I is its 53-bit rounding where IEEE's is at least 2⁻¹⁰²¹, and otherwise I <
  2⁻¹⁰¹⁹ and 0.3·I < 2⁻¹⁰²¹ is below 2⁻⁴⁸³ of √I/(1 + √I), under a quarter of its ulp, so their difference is √I/(1 +
  √I) in both; and a difference of two doubles is exact below 2⁻¹⁰²¹ (Hauser) and the 53-bit rounding above. So the
  final product's operands are the unbounded chain's, and IEEE's product is x. Where A·z² overflows, Rust's argument is
  ±∞ or NaN, but γ differs only where the other factor is 0: √I ≥ 2⁻⁵³⁷, and a nonzero bracket is at least 0.19·√I (I ≤
  1: √I/(1 + √I) ≥ √I/2 and 0.3·I ≤ 0.3·√I) or 2⁻⁵⁴ (I > 1: two doubles of at least 0.25 differ by a multiple of 2⁻⁵⁴),
  so |x| is beyond 2⁴⁸⁴ and γ is Rust's +∞ or +0. Only I = ±0 and Davies' exact zero bracket (I = 1.940363884733242 and
  1.9403638847332425, where √I/(1 + √I) and RN(0.3·I) are the same double) change: Rust's ∞·0 = NaN, γ = 1.
- `debye_huckel_extended`: the numerator n = A·z²·√I and y = (B·å)·√I finite. A finite n proves A·z² finite (∞·√I is ±∞
  or NaN), so A·z² and √I are the unbounded ones; n is its 53-bit rounding where |n| ≥ 2⁻¹⁰²¹ or it is an exact zero,
  and otherwise below 2⁻¹⁰²¹ in both evaluations. A finite y proves B·å finite; B·å is the 53-bit value where it is at
  least 2⁻¹⁰²¹ or an exact zero, otherwise below 2⁻¹⁰²¹ in both, and then y < 2⁻¹⁰²¹·2⁵¹² = 2⁻⁵⁰⁹ in both; the same for
  y itself; so the denominator 1 + y (y ≥ 0 or −0) is the unbounded one in both (1 where y < 2⁻⁵³; 1 + DBL_MAX rounds to
  DBL_MAX in both). The final quotient's operands are then the unbounded chain's (IEEE's quotient is x), or |n| < 2⁻¹⁰²¹
  with a denominator of at least 1, so both arguments are below 2⁻¹⁰²¹ and γ = 1 in both (the lemma's third case). The
  test proves Rust's γ, not always Rust's argument: (k·2⁻¹⁰⁷⁴)·√5 with k = 2,014,070,982,048,630 is (2⁵² −
  0.457…)·2⁻¹⁰⁷⁴, which IEEE rounds up to 2⁻¹⁰²² and RN₅₃ to 2⁻¹⁰²² − 2⁻¹⁰⁷⁵; γ = 1 either way (`tests/solution.tcyr`).

**The machinery** (`src/reaction.cyr`, `src/solution.cyr`) is #20's helpers — `_km_ux_from`, `_km_ux_neg`, `_km_ux_mul`,
`_km_ux_sqrt`, `_km_ux_add`, `_km_ux_sub`, `_km_ux_div` and the final `_km_ux_div_f64`, `_km_ux_mul_f64` — in four slow
paths that return the binary64 argument (`_km_keq_arg_ux`, `_km_dh_limiting_arg_ux`, `_km_dh_extended_arg_ux`,
`_km_davies_arg_ux`); the fast tests are inline compares on Rust's own intermediates. So
`equilibrium_constant(−DBL_MAX, DBL_MAX)` is 1.1278039728561764 (Rust's 1) and `(DBL_MAX, DBL_MAX)` 0.8866789123534372,
`(−10⁻³²⁰, 10⁻³²²)` 193,215.67 (0x410795fd6637a45d; Rust's 197,357), `debye_huckel_limiting(100000, 0, 10³⁰⁰)` and
`davies_activity` there 1 (Rust's NaN), `debye_huckel_extended(1, 1, 10³⁰⁸, 10¹⁵⁵, 10¹⁵⁵)` 0.9772372209558107 (Rust's
1), `(100000, 1, 10³⁰⁰, 10¹⁵⁴, 10¹⁵⁴)` 9.999999999999673·10⁻¹⁰¹ (Rust's 0) and `(1, 0, 1, DBL_MAX, 2)` 1 (Rust's ∞·0 =
NaN). A subnormal R·T can move a row across the ComputationError either way: `equilibrium_constant(−1.3434160979·10⁻³¹³,
2.2764124·10⁻³¹⁷)` is the Err where Rust's K is 1.797687·10³⁰⁸, and `(−2.9138·10⁻³¹⁹, 5·10⁻³²³)`
1.1297907864949920·10³⁰⁸ where Rust's is the Err.

**Golden comparison** (as #27 … #34). The corpus keeps the oracle's values; where the port's result differs the note
names it: `[ADR 0002 #35: the port returns the argument evaluated with an unbounded exponent (its final quotient rounded
once), then the correctly rounded exp = <x> (f:…)]` (`equilibrium_constant`;
`… = Err(ComputationError) (E:ComputationError)]` where only the port's K overflows), `… (its final product rounded
once), then the correctly rounded 10^x = <x> (f:…)]` (the limiting law, Davies) and `… (its final quotient rounded
once), then the correctly rounded 10^x = …` (the extended law). The generator computes them independently of the port
(`scripts/oracle/src/audit.rs`: the argument on `exact::Ux`, its last operation by `div_f64` / `mul_f64`, then glibc's
`exp` and `10.0_f64.powf` of it, as B2's Nernst ln; Rust's validation replicated and asserted, the ComputationError
judged on the port's K) and appends 1,866 edge rows with their own seeds (reaction 676, solution 1,190): R·T past
DBL_MAX and subnormal, R·T at 2⁻¹⁰²² from below (`KEQ_T_2M1022`) and the boundaries T = DBL_MAX/R and 2⁻¹⁰²¹/R ulp by
ulp, ΔG subnormal and huge, x past ±DBL_MAX (K = ∞: the Err, or 0, in both) and subnormal (K = 1 in both), a short
subnormal T with ΔG between the two ends' thresholds as Rust's and the port's R·T place them (the ComputationError flips
both ways, and 0 against a subnormal K); A·z² past DBL_MAX at I = ±0, at DBL_MAX ulp by ulp, against √I, against a tiny
or negative bracket and against Davies' zero bracket; A·z² subnormal; B·å and (B·å)·√I past DBL_MAX, with I = ±0; B·å
past DBL_MAX with the numerator finite (Rust's γ = 1 where log₁₀ γ is ∓L); A·z² and the numerator past DBL_MAX with the
denominator bringing them back (Rust's γ = 0 or ∞, or ∞/∞ = NaN); the numerator and B·å subnormal; B, å = ±0; the
limiting law's argument a 53-bit midpoint (in range and past DBL_MAX); every binade of every input. None of the 2,627
rows that existed before B5 changes (the survey: 1,997 chain-normal, 630 out of range with the same result; the parity
summary was identical before and after the code change). glibc's exp and pow(10, x), the generator's, are not proven
correctly rounded: `validate.py` checks every noted value against MPFR (every one agrees), and the corrected oracle
gains no row (no new row's oracle value is misrounded either). The harness maps a noted row (`_par_got_noted_k`: the
note's K or its ComputationError; `_par_got_noted` for γ) and counts it per function: **845 #35 rows** (326
`equilibrium_constant`, 27 of them the port's ComputationError where Rust's K is finite and 33 its K where Rust's
overflows; 55 `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity`). `validate.py` §23 recomputes
every row of the four functions — the argument on `mpfr_ref.UX` with every operation at 53 bits, the last one too (X),
and MPFR's exp or exp10 of X rounded once (the new `mpfr_ref.ux_exp_final`, `ux_exp10_final`): #20's target itself — and
proves the lemma on every row past validation (3,526 rows with X normal, 131 past DBL_MAX, 759 at most 2⁻¹⁰²¹), that a
row is noted exactly where the port's result differs from its reference (the corrected oracle's outputs where one
exists: 3 rows, none noted), the note naming it bit for bit or its Err, that every row where the port's fast test passes
is unnoted (1,189, 659, 637 and 729 rows in the order above), and that no row whose f64 chain stays in the normal range
is noted (1,220, 582, 541 and 655 rows: B3's strict tracers); it counts the noted rows: 27 where only the port's K
overflows (its ComputationError) and 33 where only Rust's does; 307 where Rust's γ is NaN and the value representable
(∞·0 at I = ±0 and at Davies' zero bracket, the extended law's ∞·0 and ∞/∞); 114 where Rust's value is 0 or ∞ and the
port's representable (a subnormal R·T's K, the extended law's numerator past DBL_MAX that the denominator brings back);
48 where Rust's NaN stands for a γ that is itself out of range (0 or ∞); 283 degraded (the value normal, Rust's more
than 2 ulp off: a subnormal R·T, the extended law's B·å past DBL_MAX); 33 within 2 ulp or subnormal. `cr_verify.py`'s
families `keq`, `dhl`, `dhe`, `dav` (505,982 arguments: the public function, the Rust-chain fast path where its test
passes, the slow path forced on every argument) check the same target against MPFR on x86_64 ELF, x86_64 PE and aarch64.

**Cost** (`bench_audit_b5` in the new `tests/audit2.bcyr`, before → after, the same host;
`docs/benchmarks-rust-v-cyrius.md`'s B5 section): the fast tests are inline compares on Rust's own intermediates, so
typical inputs keep their cost — `reaction/equilibrium_constant` 64 → 64 ns, `solution/debye_huckel_limiting` 86 → 85,
`debye_huckel_extended` 105 → 100, `davies_activity` 98 → 96; the slow paths `equilibrium_constant_edge` 26 → 151 ns
(R·T past DBL_MAX: Rust's exp(−0) is trivial), `equilibrium_constant_subnormal` 80 → 166, `debye_huckel_limiting_edge`
30 → 93, `debye_huckel_extended_edge` 45 → 277 (B·å past DBL_MAX), `debye_huckel_extended_overflow` 47 → 276 (A·z² past
it), `davies_activity_edge` 39 → 162 ns. **The third bench file**: in `tests/audit.bcyr` these benches would have taken
it to 219 s of `cyrius bench`'s 300 s per-file deadline (189 s before), with the audit's later batches still to come, so
they open `tests/audit2.bcyr` (31 s), as B4 split `tests/kimiya.bcyr` under the user's decision there (split rather than
raise the deadline); `tests/audit.bcyr` is unchanged.

## Addition, 2026-11-11 (audit B6): #36, #37 and #38

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B6 unchanged. On the survey's
questions 5 and 6 the user answered (2026-11-10, each the recommended option, as the project's record of the accuracy
policy has them): **"Fix all three iterative edges:** weak acid/base Newton (unbounded slow path), Peng–Robinson
coefficients only, and the flame bracket's exact sign test" and **"FFT: both** the magnitudes and the butterflies,
without the range limit" (hisab 1.4.0's transform, every operation with an unbounded exponent). Question 9 was settled by
the main session under the policy, not asked: "find_peaks' threshold is evaluated exactly". One row per family: #36
(F13), #37 (F14) and #38 (F16, its flame part; B12 adds the weak acid / base Newton and Peng–Robinson's coefficients to
the same row). Every Rust
check, Err code and order stays. B6 is the roadmap's B8.

**The treatment, per function** (each checked against `rust-old/src/spectroscopy.rs` and `thermochem.rs`, and hisab
1.4.0's `src/num/fft.rs` and `roots.rs`):

| Function | Rust's formula | Treatment | Fast test |
|---|---|---|---|
| `spectroscopy::fft_spectrum` | hisab 1.4.0's `num::fft` (bit reversal; per stage wn = (cos θ, sin θ) of θ = (−2·π)/len, w = (1, 0), t = w·d with re = w.re·d.re − w.im·d.im and im = w.re·d.im + w.im·d.re, u ± t, w *= wn), then √(re·re + im·im) of bins 0 … n/2 | #20's semantics over the whole transform: every butterfly and twiddle operation with an unbounded exponent, each √ rounded once (`_km_fft_ux`, `_km_fft_mag_ux`); a NaN or ±∞ sample keeps Rust's IEEE spectrum | `_km_fft`'s operand bounds (twiddle components 0 or ≥ 2⁻¹⁰⁰, each butterfly's d 0 or ≥ 2⁻⁹²¹), every retained output finite; per bin the square sum in [2⁻⁹⁰⁰, DBL_MAX] or re = im = 0, else that bin's unbounded √ |
| `spectroscopy::find_peaks` | w[1] > w[0], w[1] > w[2], w[1] > RN(RN(min·RN(1 − t)) + RN(max·t)) | the exact predicate w > min·(1 − t) + max·t, t the given double; min and max Rust's fold; the neighbour tests IEEE; +∞ above every finite threshold | \|RN(w − T)\| > B, T Rust's threshold, B its proven error bound |
| `thermochem::adiabatic_flame_temperature` | `f(t_lo) * f(t_hi) < 0.0`; hisab's bisection `f_lo * f_hi > 0.0` | the exact sign comparisons; Newton, its tolerances and every Err Rust's | (none: two comparisons) |

**#38, the sign tests, proven equal to Rust's except where the product underflows.** Let a = f(T₀) and b = f(T_max) and
p = RN(a·b). If either is NaN, p is NaN, and p < 0, p > 0 and every comparison of the NaN are false: both tests fail. If
one is ±0 and the other is not NaN, p is ±0 (the other finite) or NaN (0·∞), so p < 0 and p > 0 fail, and ±0 is neither
below nor above 0: both fail. Otherwise both are nonzero and not NaN, p's sign is the xor of theirs (IEEE's rule, also
for ±∞ and for a product that over- or underflows), and |p| = RN(|a·b|) is 0 (underflow: |a·b| ≤ 2⁻¹⁰⁷⁵, the tie at
2⁻¹⁰⁷⁵ to even), subnormal, normal or ∞. So p < 0 ⟺ the signs differ and |p| ≠ 0, the exact test ⟺ the signs differ;
and p > 0 ⟺ the signs agree and |p| ≠ 0, the exact same-sign test ⟺ the signs agree. They differ exactly where the
product underflows to a zero from nonzero operands. `validate.py` checks both pairs of tests on every pair of a set of
38 special and boundary values (±0, ±2⁻¹⁰⁷⁴ … ±DBL_MAX, ±∞, NaN) and 20,000 seeded pairs: they differ on 1,439 (bracket)
and 1,478 (bisection) pairs, each an underflowing product. **The bisection's test cannot change a flame result**: the
flame calls `_km_flame_bisection` only after its bracket proved a and b of strictly opposite signs, and the bisection
evaluates the same f at the same two points (deterministic), so neither form refuses there; it is exact so that the
transcription is (a tiny endothermic ΔH's ends are both negative with a product of +0, which hisab's test would bisect on;
`tests/thermochem.tcyr`). **The kept Rust behaviour**: Newton's tolerance is an absolute 10⁻³ J/mol, so at a scale where
|f| is below it from the start (any system scaled below about 10⁻³), Newton accepts its initial guess T₀ − ΔH/(40·Σn)
unconverged — CH₄/air at 10⁻²⁰⁰ gives 2204.756463878327 K where the scale-free root is 2327.18 K — and that is Rust's own
answer wherever the product does not underflow (the same system at 10⁻¹⁵⁰, within an ulp); the sign fix only makes the
scale range consistent, and the tolerance is not in scope.

**#37, the fast test, proven.** Write a = min, b = max, θ = a·(1 − t) + b·t (exact) and u = 2⁻⁵³. Rust's T is
RN(RN(a·h) + RN(b·t)) with h = RN(1 − t): h is 1 − t exactly for t ≥ ½ (Sterbenz) and within u/2 of it below (1 − t in
(½, 1], its ulp 2⁻⁵³); |RN(x·y) − x·y| ≤ u·|x·y| + 2⁻¹⁰⁷⁵ (the subnormal grid's half step); a sum of doubles rounds by at
most u relative and is exact below 2⁻¹⁰²¹. With |h|, t ≤ 1: |T − θ| ≤ u·(|a| + |b|)·(1 + u) + u·2⁻¹⁰⁷⁴ (the sum) + u·|a| +
u·|b| + 2·2⁻¹⁰⁷⁵ (the products) + (u/2)·|a| (h) ≤ 2.5u·(1 + u)·(|a| + |b|) + (1 + u)·2⁻¹⁰⁷⁴ wherever T is finite. The port
computes B = RN(RN(RN(|a|·2⁻⁵¹) + RN(|b|·2⁻⁵¹)) + 2⁻¹⁰⁷²) ≥ (1 − u)²·(4u·(|a| + |b|) + 3·2⁻¹⁰⁷⁴), above that bound; RN is
monotone, so RN(w − T) > B proves w − T > B and w > θ, and RN(w − T) < −B proves w < θ (an overflowing w − T is ±∞ on
the right side). Otherwise, and where T is ±∞, `_km_peak_above_exact` decides: w − a + t·a − t·b > 0, w and −a as
(significand, exponent) pairs, t·a and −t·b as exact TwoProducts of their significands (in [1, 4), the low parts 0 or at
least 2⁻¹⁰⁴), by `_km_sign_scaled`. A non-finite candidate can only be +∞ (NaN and −∞ fail w > w[0]), above every finite
threshold, as IEEE would compare it; a finite candidate is never below min. `validate.py` checks the bound on every row
past the checks and the fast test's decision against the exact one on every candidate it decides.

**#37's machinery, `_km_sign_scaled`** (`src/mathx.cyr`: the exact sign of Σ cᵢ·2^kᵢ for up to eight doubles cᵢ, each 0 or
of magnitude in [2⁻¹¹⁰, 4), at any integer exponents kᵢ). Round by round, with K the largest kᵢ of a nonzero term: the
terms with kᵢ − K ≥ −900 scale exactly into normal doubles in [2⁻¹⁰¹⁰, 4) and their exact sum E is formed by
Grow-Expansion with zero elimination (as `_km_exp_sign`; no overflow below 32), a nonoverlapping expansion whose largest
component has the sum's sign; each other term is below 2⁻⁸⁹⁸ in that scale, so at most eight of them sum below 2⁻⁸⁹⁵. With
no such term, or E's top component at least 2⁻⁸⁹⁰, E's sign is the sum's; otherwise E's components (exact doubles,
possibly subnormal in that scale) and the deferred terms are the next round's terms, no more than before and all below
2^(K − 888). find_peaks' terms span 2⁻²¹⁵⁰ … 2¹⁰²⁵, so the rounds end within five; the example that needs them is
w − a + t·a − t·b with a = 2⁻¹⁰³, b = DBL_MAX, t = 2⁻¹⁰⁷⁴ and w = 2⁻⁵⁰: the first three terms cancel to 2⁻¹¹⁷⁷ — w is a peak
by 2⁻¹¹⁷⁷ (`tests/spectroscopy.tcyr` checks a three-round case, 2¹⁰⁰⁰ − 2¹⁰⁰⁰ − 2⁻²¹⁰⁰).

**#36, the fast test, proven.** `_km_fft` runs Rust's transform and returns 1 when every component of wn, of each
twiddle w a butterfly uses (checked in the first block of a stage: the recurrence is input-independent, the same in every
block) is 0 or at least 2⁻¹⁰⁰ and every component of each butterfly's d is 0 or at least 2⁻⁹²¹. Then, by induction over
the operations: each product of the transform — w·d, and the recurrence's w·wn — has operands equal to the unbounded
chain's and is either an exact ±0 from a zero factor (the same signed zero in both evaluations) or of exact magnitude at
least 2⁻¹⁰²¹, where binary64's rounding is the 53-bit one unless it overflows; each sum or difference of two doubles is
exact below 2⁻¹⁰²¹ and the 53-bit rounding above unless it overflows (the pure-sum lemma, #33). So every value is the
unbounded chain's unless an operation overflowed. An overflow reaches a retained output as ±∞ or NaN: a non-finite
component of d gives t a non-finite component (w's components are never both 0; a zero one gives 0·∞ = NaN, a nonzero
one ±∞ or NaN), a non-finite u or t gives both u ± t a non-finite component, so a non-finite value at any stage makes
both outputs of the butterfly it enters non-finite, and in the last stage (one block) every value of the previous stage
enters a butterfly whose u + t is a retained bin k < n/2 (the bin n/2 is u − t of k = 0). So "every retained re and im
finite" proves the whole transform the unbounded one; the twiddle that the last k of a block computes is never used.
(The same argument shows a NaN or ±∞ sample makes every output non-finite after log₂ n stages, so such a signal never
passes the test, and B1's rule gives it Rust's IEEE spectrum.) Per bin, the square sum s = RN(RN(re²) + RN(im²)) in
[2⁻⁹⁰⁰, DBL_MAX] is the unbounded one (#20's distance_3d argument in two dimensions: no square overflowed, and a square
below 2⁻¹⁰²¹ in either evaluation is added to one of at least 2⁻⁹⁰¹, under a quarter of its ulp; both below 2⁻¹⁰²¹ would
put s below 2⁻¹⁰²⁰), a normal double whose IEEE √ is the correctly rounded one; re = im = 0 gives +0 in both; any other bin
of a valid transform takes `_km_fft_mag_ux` of its (exact) re and im — the per-bin slow path, cheap — and an invalid
transform of a finite signal the whole slow path. **The twiddles**: the bound never binds in practice (wn's smallest
component is above 2⁻⁶⁰ for every n the port can hold, and the recurrence's smallest nonzero component measures 2⁻⁵⁶ for
every n ≤ 2²²), so the test is the data's: a d below 2⁻⁹²¹ ≈ 5.6·10⁻²⁷⁸, or a sum past DBL_MAX, takes the slow path.
**The 2⁻¹⁰²² pitfall** (B4's #29): the test is on the operands, never on a product's IEEE value; n = 4's second butterfly
has w = (cos(−π/2), −1) = (6.12·10⁻¹⁷, −1), and d = 3.63·10⁻²⁹² puts w.re·d in [2⁻¹⁰²² − 2⁻¹⁰⁷⁵, 2⁻¹⁰²²), which IEEE
rounds up to 2⁻¹⁰²² and 53 bits do not; that d is below 2⁻⁹²¹, so the test refuses (`tests/spectroscopy.tcyr`).

**#36's machinery** (`src/spectroscopy.cyr`): `_km_fft` (hisab's transcription, now with the operand checks and its 1/0),
`_km_fft_spectrum_fast` (Rust's transform and magnitudes in the work table, the per-bin unbounded √ where a bin's square
sum leaves the range; 2, 1 or 0), `_km_fft_ux` (the same operations on `_km_ux_mul`, `_km_ux_add`, `_km_ux_sub`, from
`_km_ux_from` of the samples, signed zeros as IEEE's), `_km_fft_mag_ux` (the squares and their sum at 53 bits,
`_km_ux_sqrt_f64` once) and `_km_fft_spectrum_ux` (the slow path: 32·n heap bytes, then the magnitudes to `out`). The
signal must survive the fast attempt, so the IEEE transform runs in `out` only when cap ≥ 2n and `signal` does not
overlap out's first 2n values (as before, nothing is allocated then); an aliased signal now has a 16·n-byte heap work
table, which the bump allocator never frees (the slow path's 32·n bytes too). So `fft_spectrum([1.75·10³⁰⁴])` is
[1.75·10³⁰⁴] (Rust's [∞]), `([−6.96·10⁻³¹⁶])` [6.96·10⁻³¹⁶] (Rust's [0]), `([3.3·10³⁰², 3.0·10³⁰⁴])` [3.0292·10³⁰⁴,
2.9627·10³⁰⁴] (Rust's [∞, ∞]), `([DBL_MAX])` [DBL_MAX] (Rust's [∞]), `([10²⁰⁰, 1])` [10²⁰⁰, 10²⁰⁰] (Rust's [∞, ∞]: the
transform in range, the per-bin √) and `([∞, 10³⁰⁰])` Rust's [∞, ∞].

**Golden comparison** (as #23 for the whole-result notes, #33 … #35 otherwise). The corpus keeps the oracle's values;
where the port's result differs, the note names it: `[ADR 0002 #36: the port returns hisab 1.4.0's transform and the
magnitudes evaluated with an unbounded exponent (each √ rounded once) = [n], … (v:n f:… …)]`, `[ADR 0002 #37: the port
returns the threshold test evaluated exactly (w > min·(1 − t) + max·t in rational arithmetic) = [n], … (v:n i:… …)]` and
`[ADR 0002 #38: the port returns Rust's solver with an exact bracket sign test (newton, 0 steps) = <x> (f:…)]` (the
port's solver path inside; the oracle's Err rows have none). The generator computes them independently of the port: the
FFT on `exact::Ux` with glibc's cos and sin (the oracle's sincos, correctly rounded at every FFT angle: no fft row has a
correction, MPFR decides), `find_peaks` on `exact::peak_above` (W·2¹⁰⁷⁴ + T·L > L·2¹⁰⁷⁴ + T·H in `Big` integers),
the flame on `m_thermochem.rs`'s solver replica from the public Shomate API with the exact tests (asserting that it
differs from Rust's only where the bracket product underflows to a zero: Rust's ComputationError, the port's Ok). Rust's
validation is replicated and asserted. The edge rows are appended with their own seeds: spectroscopy 2,393 `find_peaks`
(a candidate at the exact threshold's successor and its neighbours ulp by ulp for fourteen (min, max) pairs — normal,
subnormal, mixed signs, ±DBL_MAX — at fourteen t including 0, 1, 1 − 2⁻⁵³, 2⁻¹⁰⁷⁴ and ½ ± an ulp; subnormal data on a
grid of t; ±0 mixes; a +∞ candidate; seeded near-ties over every binade) and 515 `fft_spectrum` (every power of two to
256: DBL_MAX impulses, constants and alternations, √DBL_MAX and 2⁻⁵¹¹ samples, subnormal ones, the data bound 2⁻⁹²¹ ulp by
ulp, huge, subnormal, mixed and near-cancelling seeded signals, n = 4's product below 2⁻¹⁰²², NaN and ±∞ beside huge
samples; seeded signals over every binade); thermochem 443 flame rows (the survey's example and its endothermic mirror,
six real systems scaled by 10⁻³²⁰ … 10⁻¹⁵⁰ at two T₀, the scale where |f(T₀)·f(T_max)| crosses 2⁻¹⁰⁷⁵ ulp by ulp, a
subnormal ΔH against moles from 10⁻⁹ to 10², endothermic tiny systems, seeded fuel mixtures at tiny scales). Of the rows
before B6 exactly the expected change: **41 of 213 `fft_spectrum`** (40 where a bin of Rust's is 0, ±∞ or NaN and the
port's representable — 2,158 such bins — and 1 where Rust's NaN is the value ∞; a magnitude-only fix would have repaired
28 of the 41, so the butterflies matter for 13), **3 of 994 `find_peaks`** (`[−0, 5·10⁻³²⁴, 0]` at t = 0.75 and
1 − 2⁻⁵³, `[10, 10.2, 10, 11, 10]` at 1 − 2⁻⁵³), **16 of 1,937 `adiabatic_flame_temperature`** (four real systems at
10⁻²⁰⁰ and 10⁻¹⁷⁰, each at T₀ = 298.15 and 500 K). The harness maps a noted row (`_par_map_toks` for the two vector
results, `_par_got_noted` for the flame) and counts it: **418 #36 rows, 222 #37 rows, 273 #38 rows**.

`validate.py` §24 proves every note independently. **#38**: every flame row recomputed with its replica `r_flame` in the
port's form (`exact`), noted exactly where the result differs from the oracle's, the note naming the value bit for bit
and the solver path; a noted row's ends have opposite signs and a product underflowing to a zero (Rust's
ComputationError against the port's Ok), and every row whose product is no such zero is unnoted (strict); 2,345 rows past
validation, 47 with an underflowed zero product of equal signs (Err in both), 111 with a subnormal nonzero product
(bracketed in both); B4's replay of the flame's Shomate evaluations now includes the port's own (664,046, all in #34's
box). **#37**: every row with `Fraction`s, noted exactly where the indices differ, the note's tokens bit for bit; on
every row past the checks Rust's T within the proven bound of θ and the port's B above it, the fast test's decision the
exact one on every candidate it decides, and no row noted whose every candidate the fast test decides; 6,546 candidates,
3,993 decided by the fast test, 233 where Rust's threshold decides otherwise (109 at subnormal or zero thresholds, 16 with
±DBL_MAX ends, 108 other normal near-ties). **#36**: every row with a finite signal recomputed on `mpfr_ref.UX` (the
transform as transcribed with MPFR's correctly rounded cos and sin, every operation at 53 bits at MPFR's widest exponent
range, each √ by `ux_sqrt_final`), noted exactly where it differs from the oracle's (no fft row has a correction), the
note's tokens bit for bit; a signal with a NaN or ±∞ never noted; B3's strict tracers over every product, sum and square
(the twiddle recurrence included): no row whose chain stays in the normal range is noted (226 such rows), nor one where
the port's fast test passes with every bin's square sum in range (194); rows longer than 1,024 whose chain stays in range
are covered by that lemma alone (the UX replica is slow). The noted rows by their worst bin: 408 where a bin of Rust's is
0, ±∞ or NaN and the port's representable, 1 where Rust's NaN is the value ±∞, 8 degraded (normal, Rust's more than 2 ulp
off), 1 within 2 ulp or subnormal; 9,177 noted bins.

`cr_verify.py` adds `flame` (adiabatic_flame_temperature and the transcribed bisection's exact same-sign test, 14,010
arguments: tiny scales, around the scale where the product stops underflowing, endothermic and subnormal ΔH, typical,
every Err; reference Rust's solver in Python on the generated SHOMATE table with the exact tests), `pk3`, `pk5`, `pk8`,
`pk16` (find_peaks of 3, 5, 8 and 16 values: the public function, the fast test's scan, every candidate through the exact
predicate forced; reference `Fraction`; 95,000 arrays each with near-ties, subnormal data, ±DBL_MAX ends, ±0 and
non-finite values) and `fft1` … `fft256` (fft_spectrum of n = 1, 2, 4, …, 64 and 256: the public spectrum, the fast path's
where its transform test passes, the unbounded transform forced; reference the transform on `mpfr_ref.UX`; 110,718
signals: per region 3,000 for n ≤ 8, 1,500 for 16, 800 for 32, 400 for 64 and 60 for 256 — the UX reference costs about
16 operations a butterfly — over binade-uniform, typical, huge, subnormal, mixed, near-cancelling, the bounds 2⁻⁹²¹,
2⁻⁵¹¹, 2⁵¹¹ and √DBL_MAX, non-finite samples and n = 4's product below 2⁻¹⁰²²) and `signsc` (`_km_sign_scaled` itself,
140,110 sums of up to eight terms against `Fraction`: random terms over its whole exponent range, exact cancellations
decided one to three rounds down, a round's remainder of 2⁻¹⁰¹⁰ … 2⁻⁹⁰¹ of its scale against a larger deferred term), on
x86_64 ELF, x86_64 PE and aarch64. That last region found a slip in review: the round's bound 2⁻⁸⁹⁰ had been written as
the double 2⁻¹⁰¹⁸, which let a round's tiny remainder decide against a larger deferred term (12,008 of its sums; no
golden row or `find_peaks` argument reached it); corrected, with a three-round case in `tests/spectroscopy.tcyr`.

**Cost** (`tests/kimiya.bcyr`, before → after on the same host, alternating runs of HEAD and the batch;
`docs/benchmarks-rust-v-cyrius.md`'s B6 section): `spectroscopy/fft_spectrum_256` 27.3 → 30.5 µs (+12 %),
`fft_spectrum_1024` 129.2 → 143.1 (+11 %), `fft_spectrum_4096` 605.6 → 669.4 (+11 %): the per-butterfly operand checks
are about 6 % (measured without them: 29.1, 139.6, 638.5), the rest the restructured magnitude pass, its range test and the
extra state; skipping the checks at k = 0 (w = (1, 0), exact products) gained under 1 % and was not kept, merging the
retained-output and magnitude passes gained 1 % and was. `find_peaks_1000` 19.0 → 20.3 µs (+7 %: Rust's threshold and
the bound once, a subtraction and two compares a candidate), `thermochem/adiabatic_flame_ch4_air` 3.37 → 3.40 µs (two
compares for a multiply: unchanged within noise). The new `bench_audit_b6` (`tests/audit2.bcyr`, before → after:
Rust's chain, then the port's): `fft_spectrum_256_huge` 27.7 → 318.1 µs and `fft_spectrum_256_subnormal` 35.4 → 349.8 µs
(the unbounded transform: about 11× Rust's chain), `fft_spectrum_256_squares` 27.4 → 45.1 µs (the per-bin unbounded √ on
all 129 bins), `find_peaks_1000_subnormal` 19.5 → 39.0 µs (about half the candidates through the exact predicate),
`adiabatic_flame_tiny` 0.96 → 1.48 µs (Rust's ComputationError there; here Newton's guess); `tests/audit2.bcyr` 45 s.

## Addition, 2026-11-12 (audit B7): #39

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply to B7 unchanged. The survey's
question 7 was not asked: the main session settled it under the user's policy, with its remedy (a): **"Balancing:
multi-word integers, no size limit"** — Rust's own fraction-free elimination with the integers unbounded, not a modular
null space (its (b)), and no width past which the port gives up (the project's record of the accuracy policy: "balancing
uses Rust's elimination on multi-word integers with no size limit"). One row for the family: #39 (F15). Every Rust
check, Err code and order stays. B7 is the roadmap's B9.

**The treatment** (checked against `rust-old/src/stoichiometry.rs`: `balance_equation`, `reduce_row`, `gcd_u64`,
`balance_overflow`). The result is what Rust's algorithm returns with its i64 and i128 replaced by mathematical integers:
the steps exactly as written — the pivot search over rows rank … m − 1, the swap, reduce_row of the pivot row, its sign
flip, row_r ← pivot·row_r − factor·row_pivot then reduce_row for every other row with a nonzero factor, the lcm of the
pivots by lcm ← (lcm / gcd(lcm, piv))·piv, x[f] = lcm and x[c] = −(row[f]·(lcm/piv)), reduce_row(x), the zero check, the
positive count and its flip, the u32 conversion, and the verification on the original rows. So the overflow Errs
disappear — the elimination's checked_mul, checked_mul and checked_sub, the pivot row's checked_neg, reduce_row's
`i64::try_from(g)`, the lcm's gcd conversion and checked_mul, x's checked_mul and checked_neg, the final flip's
checked_neg — and every other Err stays at its place: ambiguous (a second free column) and only the trivial solution (no
free column), the zero coefficient, mixed signs, a coefficient above u32::MAX, the verification, every
Err(InvalidInput), and the `cap` / `out` check after the computation. All the Errs involved are ComputationError, so the
visible change is Err → Ok; an overflowing input with too small a cap, or a null `out`, now reaches the buffer check and
is Err(InvalidInput). The verification keeps Rust's i128 sums: with u32 coefficients and |entry| ≤ 2⁵³ each partial
sum is below n·2⁵³·(2³² − 1) < 2¹²⁷ for every n below 2⁴² (a row of 32 TiB), so its checked_add cannot fail, as in Rust.

**The structure** (`src/stoichiometry.cyr`). `balance_equation` validates as before, then runs `_km_balance_i64`, the
former body unchanged except that each overflow check returns "overflow" (2) instead of the Err: 0 (x holds the u32
coefficients, past the zero, sign and u32 checks), 1 (one of Rust's other ComputationErrors), 2. While no check fails,
every i64 value is the exact integer, so 0 and 1 are the exact algorithm's outcome. At 2 the multi-word path
`_km_balance_big` runs the same steps from the original matrix; then the shared i128 verification and the buffer check.
`_km_balance_equation(…, force)` skips the fast path (the tests' and `cr_verify.py`'s check of the slow path). In-range
inputs take exactly the old code: `balance_ch4_combustion` 805 → 778 ns, `balance_kmno4_hcl` 2.24 → 2.18 µs
(alternating runs; the fast path moved into its own function).

**The integers** (`_km_sb_*`, beside the i128 helpers in `src/stoichiometry.cyr` rather than in `src/mathx.cyr`:
balancing is their only user, and `src/mathx.cyr`'s `_km_bn_*` are unsigned, fixed at 8 × 32-bit limbs, with no
division by a multi-word divisor). A slot is its signed length (−len for a negative value, 0 for zero) and len limbs of
31 bits, little-endian, the top one nonzero. 31-bit limbs keep every intermediate below 2⁶³ — a limb product
(2³¹ − 1)² plus a limb and a carry is at most 2⁶² − 1, a two-limb numerator below 2⁶² — so the signed words never wrap
and Cyrius's logical `>>` only shifts non-negative words. The operations: from a word (i64::MIN as −2⁶³), copy,
magnitude compare, add and subtract magnitudes, the signed difference, the schoolbook product, Knuth's Algorithm D for
the quotient and remainder (as Hacker's Delight's divmnu: the divisor normalised to a top limb of at least 2³⁰, q̂ from
the top two limbs, below 2³¹ after the two-limb test, whose product q̂·v[n − 2] is formed only then, below 2⁶², the
multiply-subtract with a non-negative carry, the add-back when q̂ was one too large), the exact quotient, and the gcd of
magnitudes by **Lehmer's algorithm** (Cohen's Algorithm 1.3.7): with the larger operand of three limbs or more and the
other within a limb of it, Euclid runs on their top 62 bits (one shift) for both extreme ratios (â + A)/(b̂ + C) and
(â + B)/(b̂ + D) while the quotients agree — so each is the full operands' quotient — below 2³⁰, with the cofactors
below 2³⁰; then (p, q) ← (A·p + B·q, C·p + D·q) in one pass of signed limbs (|t| < 2⁶², the carry an exact division by
2³¹), a unimodular step that keeps the gcd whatever happens (a negative total is negated); where no quotient is admitted,
or the operands differ by two limbs, one full division step; below two limbs, Euclid on words. reduce_row's gcd starts at
the row's entry with the fewest limbs and stops at 1 (gcd is associative and commutative: Rust's left fold's value).
Lehmer took a 39 × 40 dense matrix from 1.8 s (about 844,000 of Euclid's division steps) to 0.7 s.

**Which overflows can become Ok, proven.** At the end of the elimination, with exactly one free column f, every pivot row
k is primitive (reduce_row follows its every change) and has at most two nonzero entries: its pivot p_k > 0 at c_k and
r_k at f (the other pivot columns are eliminated; a free column left of c_k is zero in every row at or below the rank
then, and the later pivot rows are combinations of those rows only). So x = (x[c_k] = −r_k·L/p_k, x[f] = L) with L the
lcm of the p_k, and its gcd is gcd_k((L/p_k)·gcd(p_k, r_k)) = gcd_k(L/p_k) = 1 when no r_k is 0 (for each prime some p_k
carries L's full power): reduce_row(x) changes nothing that passes the zero check, and x[f] stays L. Hence an overflow of the
lcm (L ≥ 2⁶³), of x's product (|x[c]| ≥ 2⁶³) or of its negation (x[c] = 2⁶³) leaves an entry past u32 (or a zero) on
unbounded integers too: Err either way, unnoted. A reduced row whose gcd is 2⁶³ has every nonzero entry −2⁶³ (the only
i64 of that magnitude), so on unbounded integers it reduces to a row of 0s and −1s, orthogonal to no positive x: Err
again. So only the elimination's products, its differences and the pivot row's negation can turn into Ok; all four do in
the corpus (pivot·a 389 rows, factor·a 16, the difference 7, the negation 1: `[[2³¹, 0, −1, 1 − 2³¹], [1, −2³², 2³¹,
2³¹ − 1], [1, 1, −1, −1]]`, whose eliminated row 1 is [0, −2⁶³, 2⁶² + 1, 2⁶² − 1] with gcd 1, gives [1, 1, 1, 1]).
**Unreachable sites**: the lcm's gcd conversion (gcd(lcm, piv) ≤ lcm < 2⁶³), reduce_row(x)'s gcd (≤ x[f] = lcm) and the
final flip (x[f] = lcm > 0 survives reduce_row, so at least one coefficient is positive and the flip never runs); they
are kept as Rust wrote them. `validate.py` checks the gcd claim on every row with a one-dimensional null space and no zero
coefficient (2,588).

**The size bound, proven** (the tables need no growth). Every step multiplies a row by a nonzero integer and adds a
multiple of another, divides it by a positive gcd, negates it, or swaps two rows, so at every point each row is a
nonzero rational multiple of the row of rational Gauss–Jordan elimination it stands for (or both are zero). After k
pivots, with R the pivot rows and S their columns (P the matrix with the rows in their order then), a non-pivot row i
is a multiple of the Schur complement row P_i − P_{i,S}·P_{R,S}⁻¹·P_R, whose entry j times det P_{R,S} is the (k + 1)-minor
det P_{R ∪ {i}, S ∪ {j}} (Schur's determinant formula); a pivot row j is a multiple of row j of P_{R,S}⁻¹·P_R, whose entry
l times det P_{R,S} is, by Cramer's rule, the k-minor of P_{R,S} with column c_j replaced by column l. So a row just
reduced (divided by the gcd of its entries; the sign flip changes no magnitude) is ± the primitive vector of the line
spanned by an integer vector of minors of order at most q = min(m, n), each of its entries at most such a minor; an
original row is its own 1-minors. By Hadamard's inequality a minor is at most the product of its rows' Euclidean
lengths, each at most its whole row's: a row with z nonzero entries of largest magnitude μ has ‖row‖ ≤ √z·μ < 2^b,
b = bitlen(μ) + ⌈bitlen(z)/2⌉ (0 for a zero row), so every minor of order ≤ q is below 2^B, B the sum of the q largest
row b's — or of the column b's, the same argument on the transpose; `_km_balance_bits` takes the smaller. Therefore
every stored entry is below 2^B, and each product pivot·a or factor·a (two stored entries) below 2^(2B), their difference
below 2^(2B + 1). The lcm: pivot row k is (the minor vector w_k)/g_k with w_k's pivot entry D = det P_{R,S}, so
p_k = |D|/g_k divides |D|, L and every partial lcm divide |D| < 2^B; and |x[c_k]| = |r_k|·L/p_k = |w_k[f]|·L/|D| ≤
|w_k[f]| < 2^B, the product row[f]·(lcm/piv) being x[c_k]'s magnitude itself. So a matrix slot holds B/31 + 2 limbs (one
past ⌈B/31⌉: `_km_sb_divexact` writes len(a) − len(d) + 1 quotient limbs before trimming, at most ⌈B/31⌉ + 1), a work
slot (the scratch row, x, the products, the gcd's and the division's operands) 2·(B/31 + 2) + 2 limbs, at least
⌈(2B + 1)/31⌉ + 1, and the division's normalised dividend one word more. `validate.py` checks every stored value, product
and difference of the exact elimination against these bounds on every golden row (the largest stored value 1,823 bits,
the largest product 3,641 bits).

**Heap.** `_km_balance_big` takes m row pointers, m·n matrix slots of B/31 + 3 words and 2n + 11 work slots of
2·(B/31) + 8 words, in one table: a local one of 4,096 words (the survey's 2 × 3, B = 66: 236 words), the bump allocator
past it (which never frees), Err(ComputationError) when the allocation fails (the allocator refuses one past 2 GiB; a
table past 2⁴⁰ words is refused before it is sized). The bench's 30 × 31 matrix of entries below 2⁴⁵ (B = 1,444) takes
423 KB a call, a 40 × 41 of 2⁵³ entries (B = 2,280) 1.1 MB. Rows are swapped through their pointers.

**Golden comparison** (as #36 and #37's whole-result notes). The corpus keeps the oracle's values; where the port's
result differs the note names it: `[ADR 0002 #39: the port returns Rust's fraction-free elimination on unbounded
integers = [n], … (v:n i:… …)]`. The generator computes it independently of the port: `audit::port_balance`, Rust's
algorithm on `exact::Int` (a sign and a `Big` magnitude, with Knuth's division in base 2⁶⁴ and Euclid's gcd), with a
tracer of every value Rust forms with checked arithmetic; it asserts that the port's result differs from Rust's only
where Rust's run overflows, and then Rust's result is the ComputationError. The edge rows are appended with their own
seed (`emit_balance_b7`, 521 rows): the survey's example; Rust's first overflow at each site it can reach (pivot·a,
factor·a, the difference, the pivot row's negation, a reduced row's gcd of 2⁶³, the lcm, x's product and its negation);
matrices orthogonal to a known x (its last entry 1, the last column solved, large entries in the other columns) of 3 …
12 species; x with an entry 0, −1 or −1000, 2³² − 1 (Ok), 2³² or 2³² + 1 (past u32); n − 2 rows (ambiguous) and n or
n + 1 random rows (only the trivial solution), after the overflow; a duplicate, a negated duplicate, the sum of two rows,
a zero row and a rotation; ±2⁵³ entries; and large sizes, dense and sparse, of 16, 20, 24, 32 and 40 species. Of the rows
before B7 exactly the expected change: **104 of the 3,329** `balance_equation` rows, all Rust's ComputationError and the
port's Ok, by (m, n) 55 (2, 3), 10 (8, 9), 6 each (6, 7), (9, 9) and (10, 9), 5 (9, 8), 3 each (8, 8), (7, 8), (8, 7) and
(7, 7), 2 each (6, 6) and (7, 6); Rust overflows on 20 more whose exact result is still the ComputationError (12
ambiguous, 2 trivial, 1 zero coefficient, 5 past u32); the 78 #7 rows keep their notes. Of B7's rows 309 noted. The
harness maps a noted row (`_par_map_toks`) and counts it: **413 #39 rows**.

`validate.py` §25 proves every note independently: each row with an integer counterpart recomputed with the port's
semantics on Python ints (Rust's steps, no i64 check), the result checked against the exact null space (§18's Fraction
RREF: an Ok is the primitive positive basis vector of a one-dimensional null space, an Err a null space of another
dimension or a vector with a zero, mixed signs or an entry past u32), noted exactly where it differs from the oracle's,
the note's tokens the port's result; the strict tracer (`b7_rust_site`, the site of Rust's first i64 overflow): every
noted row is one where Rust overflows, and every row where it does not is unnoted. 3,772 rows; Rust overflows on 620,
413 noted, 207 still the ComputationError (103 ambiguous or trivial, 38 mixed, 45 past u32, 21 zero); by Rust's first
overflow: pivot·a 587 (389 Ok exactly), factor·a 18 (16), the difference 8 (7), the pivot row's negation 1 (1), a reduced
row's gcd of 2⁶³ 3 (0), the lcm 1 (0), x's product 1 (0), x's negation 1 (0); the size bound and the gcd claim on every
row, as above.

`cr_verify.py` adds `bal23`, `bal34`, `bal56`, `bal89`, `bal1213` and `bal2425` (balance_equation of 2 × 3 … 24 × 25
matrices: the public result, the i64 fast path's own outcome and the multi-word path forced, against Rust's algorithm on
Python ints with a tracer of its checked operations — the fast path must report an overflow exactly where Rust's
arithmetic overflows; 54,306 matrices over seven regions: known x with large and small entries, x with a zero, negative
or u32-edge entry, random small and random large entries, ±2⁵³ entries, dependent rows), on x86_64 ELF, x86_64 PE and
aarch64. `tests/stoichiometry.tcyr` adds `cyrius_adr39_unbounded_integers` (the survey's example; each reachable site;
each Err after an overflow; u32::MAX against 2³²; the cap and a null `out` after a multi-word Ok; Rust's error before
the cap), `cyrius_adr39_size_and_forced` (a 32 × 33 dense matrix; 600 seeded matrices, the forced path against the fast)
and `cyrius_adr39_multiword` (the helpers at their edges: −2⁶³, ±2⁵³, the add-back step, Fibonacci gcds, a negative
Lehmer total). No Rust test changes: `balance_intermediate_overflow_is_error`'s and
`balance_large_entries_do_not_panic`'s matrices overflow Rust's i64, and their exact vectors exceed u32 ([1, 1, 2⁵³] for
the first), so both stay the ComputationError.

**Cost** (`tests/kimiya.bcyr`'s balance benches unchanged, as above; the new `bench_audit_b7`, `tests/audit2.bcyr`,
before → after: Rust's ComputationError at its first overflow, then the port's exact result):
`balance_survey_2x3` 0.12 → 3.6 µs, `balance_overflow_10x9` (a golden row) 6.6 → 57 µs, `balance_dense_30x31` (entries
below 2⁴⁵, products past 3,000 bits) 5.4 µs → 158 ms; `tests/audit2.bcyr` 45 → 49 s.

## Addition, 2026-11-12 (audit B10a): the unbounded exponentials (for #40)

**The decisions.** B1's section above quotes the user's decisions verbatim; they apply unchanged. For the survey's
family F10 (e^x times a prefactor) the user approved new primitives, recorded as "New unbounded-exponent primitives
(exp, expm1, exp10, pow, log10) for the e^x × prefactor family … : Rust's chain with an unbounded exponent, last step
rounded once" — #20's semantics, the survey's remedy (a). B10 is the survey's (and the roadmap's) B10; its first part,
B10a, builds and verifies only the exponentials, so no function's result and no golden row changes. Row #40 is reserved
for F10's functions (Arrhenius, collision theory, Eyring ×2, first-order, `atoms_remaining`, the decay chain, van 't
Hoff, Clausius–Clapeyron, Morse, Antoine and `fugacity_pr`), whose behaviour changes in the three later sub-batches
(kinetics; nuclear; van 't Hoff, Clausius, Antoine, Morse and `fugacity_pr`). The semantics below were settled by the
main session under the user's policy (not asked).

**The primitives** (`src/mathx.cyr`, "e^x, e^x − 1 and 10^x with an unbounded exponent"): `_km_cr_exp_ux`,
`_km_cr_expm1_ux` and `_km_cr_exp10_ux` take an unbounded-exponent value x = xs·2^xe (a normalised `_km_ux_*` pair) and
return RN₅₃(e^x), RN₅₃(e^x − 1) and RN₅₃(10^x) as a pair (s, e), 1 ≤ |s| < 2: the significand correctly rounded to 53
bits, ties to even, at any exponent, with no subnormal grid and no overflow at binary64's limits (e^-800 =
1.795…·2^-1155, e^1000 = 1.619…·2^1442). In each F10 chain the exponential is never the final operation, so the chain's
last product (or expm1's subtraction) rounds once into binary64 with the existing `_km_ux_*_f64`. Each has `_fast` (the
fast path alone: (s, e, 1), or (0, 0, 0) where its rounding test fails) and `_acc` (the accurate path alone) entry
points for the tests and `cr_verify.py`. No binary64-argument wrapper: a consumer's double x is `_km_ux_from(x)` first.

**Range and saturation.** Exact wherever the result fits MPFR's widest exponent range (`mpfr_get_emin_min()` = 1 − 2⁶²,
`mpfr_get_emax_max()` = 2⁶² − 1; MPFR's value lies in [2^(E − 1), 2^E), the port's exponent is E − 1), the range
`mpfr_ref.UX` works in, so the references need no special case. The exponents returned lie in [−2⁶² + 154, 2⁶² − 155]
for e^x and e^x − 1 and [−2⁶² + 807, 2⁶² − 808] for 10^x; |x| reaches about 3.2·10¹⁸ (e^x) and 1.4·10¹⁸ (10^x). Past the
range they saturate as MPFR does: (+∞, 0) where e^x or e^x − 1 overflows (x > OVF·2⁶¹, OVF = 0x1.62e42fefa39efp0) or
10^x does (x > OVF₁₀·2⁶⁰, 0x1.34413509f79fep0), and (+0, 0) where e^x (x < −UNF·2⁶¹) or 10^x (x < −UNF₁₀·2⁶⁰)
underflows, UNF = OVF and UNF₁₀ = OVF₁₀. MPFR's round-to-nearest underflow gives its smallest number rather than 0 for a
result in [2^(emin − 2), 2^(emin − 1)); none lands there: the last nonzero result is a regular number (e = −2⁶² + 154)
and the next is far below 2^(emin − 2), which `scripts/oracle/mathx_tables.py` asserts while it bisects the edges with
MPFR (`_KM_CR_EXPUX_OVF` … `_KM_CR_EXP10UX_UNF`). e^x − 1 tends to −1 and never underflows. The input pair is expected
inside MPFR's range too (|xe| < 2⁶², as every chain value is); e^x and 10^x give 1 and e^x − 1 gives x itself for any
exponent below.

**Specials and the tiny end**, as IEEE's exp, expm1 and exp10: ±0 gives 1 (e^x − 1: ±0, the sign kept), +∞ gives +∞, −∞
gives +0 (e^x − 1: −1), NaN gives NaN, each with e = 0. xe < −54 gives 1 for e^x (|x| < 2⁻⁵⁴) and x for e^x − 1 (|x| <
2⁻⁵³: x²/2 + … is below half an ulp of x, and below a quarter of it where |x| is a power of two); xe < −56 gives 1 for
10^x (|x·ln 10| < 2^-54.8). The boundaries are checked exactly: e^(−2⁻⁵⁴) is 2⁻¹⁰⁹ above the midpoint 1 − 2⁻⁵⁴ and e^x −
1 at −2⁻⁵³ 2^-161.6 below its midpoint, both sent to the accurate path.

**The algorithm**, in three tiers of x (a 53-bit number):
- binary64's own range (e^x: |x| ≤ 708; 10^x: |x| ≤ 307; e^x − 1: x ≤ 708): x is a double whose result is a normal
  double, and the result is `_km_cr_exp(x)`, `_km_cr_exp10(x)` or `_km_cr_expm1(x)` itself (their paths and proofs, ADR
  0006), its pair taken from the bits — so the significand is `_km_cr_exp`'s wherever that is normal, by construction
  (and checked);
- the rest of the binary64 cores' range (e^x: 708 < |x| < 1024; 10^x: 307 < |x| < 512; e^x − 1: 708 < x < 1024): the
  existing cores without their final scaling — exp's fast core (`_km_exp_core`; 10^x on exp10's x·ln 10 = yh + yl) and
  the rounding test on hi + lo with ε = 2⁻⁶⁹, the pair by `_km_ux_norm(u, E)`; when it fails, the accurate core
  (`_km_exp_acc_core`) and `_km_rn3` at the same scale. K = RN(x·512/ln 2) stays below 2²⁰ (756,389; 870,800 for x·ln
  10), which the cores' K·L1 needs; the fast core's dropped K·L3 grows from 2^-79.2 to 2^-78.8, inside its bound's
  margin (0.9·2⁻⁷⁰ against ε = 2⁻⁶⁹);
- beyond (9 ≤ xe ≤ 61, to the edge): a new reduction, then the same cores' tails (`_km_exp_core_r`, `_km_exp_acc_r`: the
  code after K and r, unchanged). x·C = K + f with C = 512/ln 2 (10^x: 512·log₂ 10, as 10^x = 2^(x·log₂ 10)), |f| ≤ ½,
  the result 2^E·2^(j/512)·e^r with K = 512·E + j and r = f·ln2/512. The fast path reduces in double-double while |x·C|
  < 2⁵¹ (`_km_expux_dd`: x·C by two TwoProducts and a rounded third term of C as a triple-double, K by the 1.5·2⁵²
  shift, f = x·C − K moved once by ±1 into |f| ≤ ½, r = f·ln2/512 by TwoProduct: within 2⁻¹¹²), and past it, and always
  on the accurate path, in integers (`_km_expux_reduce`, Payne–Hanek's technique without the dropped high bits): the
  table `_KM_CR_EXPUX_C` (`_KM_CR_EXP10UX_C`) holds C·2⁻⁷⁵ with 64 zero bits first, so its 288-bit window from bit xe −
  8 starts above C's leading bit and x·C = M·V·2⁻²⁵⁶ + δ with the integer part whole (M the 53-bit significand, |δ| <
  2⁻²⁰³); M·V in eleven 32-bit columns gives K's 72 bits and the 256-bit fraction, rounded to nearest; its top 192 bits
  as a triple-double times ln2/512 (`_KM_CR_EXPUX_LH`, `LM`, `LL`) give r within 2⁻¹⁶² (measured 2^-168.3). A negative x
  takes −K and −r.

e^x − 1 past x = 708 (E ≥ 1021) rounds as e^x does: the −1 is below 2⁻¹⁰²¹ of e^x, inside the fast test's margin and far
below the accurate path's bound; x ≤ −64 gives −1.

**Error bounds, the exact and tie cases, and the proof.** The fast path's error is the cores' (|hi + lo − e^r·2^(j/512)|
< 0.9·2⁻⁷⁰ with the reduction's 2⁻¹¹² or 2⁻¹⁶²; measured 2^-71.07 and 2^-71.09 over 24,000 arguments each against
500-bit MPFR), so the same ε = 2⁻⁶⁹ and the same test; the accurate path's is below 2⁻¹⁵⁰ relative (the reduction's
2⁻¹⁶², the accurate core's 2⁻¹⁵⁵, `_km_rn3` exact; measured 2^-155.15), as `_km_cr_exp2_scaled`'s. e^x and e^x − 1 are
exact only at x = 0 (Lindemann–Weierstrass), whatever the exponent, so never a 53-bit midpoint. 10^x is rational only
for integer x and dyadic only for x ≥ 0: 10ⁿ = 5ⁿ·2ⁿ is a 53-bit number for 0 ≤ n ≤ 22, a midpoint for n = 23 (5²³ has
54 bits and is odd: the tie, to even, special-cased in `_km_cr_exp10_acc` as before), and neither for n ≥ 24 (5ⁿ odd
with 56 bits or more) or n < 0 (not dyadic); so the unbounded exponent adds no exact or tie case, and all of them are in
the first tier. **Proof:** in the first tier, the binary64 functions' own (ADR 0006): Lefèvre–Muller's worst cases
for e^x, and for 10^x and e^x − 1, which have no published table, an accurate path 18 bits past their expected hardest
case. Elsewhere no worst-case search exists, and B10a follows the precedent of pow and `_km_cr_exp2_scaled` (#13): the fast
path is proven (a failed test costs only the accurate path), the accurate path is proven correct for results not within
2⁻¹⁵⁰ (relative) of a 53-bit midpoint, a random argument is that close with probability about 2⁻⁹⁶, and the domain past
the first tier has fewer than 2⁵⁹ arguments, so the expected number of misrounded arguments in all of it is below 2⁻³⁷.
A provably terminating last resort (a multi-word Ziv loop: e^x is transcendental for rational x ≠ 0) would need an
arbitrary-precision ln 2, exp series and reduction on growable integers — several hundred lines for that 2⁻³⁷ — and is
not built.

**For the consumer sub-batches.** The largest exponent returned is 2⁶² − 155 and the smallest −2⁶² + 154; `_km_ux_mul`
adds i64 exponents, so a product of two such results stays inside i64 but three, or two with a prefactor's exponent, can
wrap past ±2⁶³: a chain with several exponentials must bound them first (any |e| past about 2⁶¹ already decides the
final result's ±∞ or 0 unless another factor compensates). A saturated value is MPFR's (+∞ or +0), so 0·(+∞) in a chain
is NaN where the true value may be 0; what the functions do there is theirs to decide, under the NaN policy (the
correctly rounded ±∞ or 0).

**Verification.** `cr_verify.py` adds `exp_ux`, `expm1_ux` and `exp10_ux` (448,188 arguments: specials and extreme
exponents to ±2⁶³, every binade from 2⁻¹²⁰⁰ to 2⁶⁵ of both signs, binary64's range, the tiny end, each tier's edges,
doubles near multiples of ln2/512 (log₁₀2/512) with K to 2⁷¹, Lefèvre–Muller's hard cases, large arguments to 3.2·10¹⁸,
the saturation edges, expm1's −1 region and exp10's integers −400 … 400 and 10²³) against MPFR at its widest range
(`mpfr_ref.UX.exp`, `.expm1`, `.exp10`, `UX.parts`): the function, its fast path, the accurate path forced on every
argument, and at every normal double x the binary64 function's pair, on x86_64 ELF, x86_64 PE and aarch64; about a
minute of the gate. `tests/mathx.tcyr` adds `cyrius_b10a_exp_ux`, `cyrius_b10a_expm1_ux` and `cyrius_b10a_exp10_ux` (40,
29 and 27 MPFR rows generated as bit patterns, each by all three entry points; the specials; 4,000 seeded arguments each
against the accurate path and the binary64 function).

**Cost** (`tests/audit2.bcyr`'s new `bench_audit_b10a`, beside the binary64 functions at the same arguments): e^x at
−20.05 46 ns (`_km_cr_exp`) → 54 ns (`_km_cr_exp_ux`), e^x − 1 at ln 2 54 → 65 ns, 10^x at −4.75 60 → 69 ns; the second
tier (e^-720.5) 84 ns; the double-double reduction (e^∓10⁴, e^(10¹²), e^(10⁴) − 1, 10^(10⁴)) 158–165 ns; the integer
reduction's fast path (e^(10¹⁵)) 614 ns; the accurate path alone 1.86 µs (−20.05), 2.14 µs (−10⁴) and 2.17 µs
(10^(−10¹²)). The consumers call these only on their slow paths.

## Addition, 2026-11-12 (M5b): #2 realised

**The four fixed-formula gradients** (`src/sensitivity.cyr`) return a `Result`, validated as the module's own Rust
functions validate — the codes and the order of `reaction::equilibrium_constant` (T, then ΔG), `kinetics::arrhenius_rate`
(T, then A, then Ea), `electrochemistry::nernst_potential` (T, then n, Q and E°; n a u8 word, outside 1 ..= 255
InvalidInput, M3a's rule) and `reaction::gibbs_free_energy` (T, then ΔH and ΔS) — and Err(ComputationError) where an
entry is not finite (as `equilibrium_constant` refuses an overflowing K). The vectors go to the caller's buffer of two
or three f64 (`lineweaver_burk`'s form), written only on Ok; `gibbs_temperature_sensitivity` is Ok(−ΔS).

**The closed forms** are each fixed tape's reverse sweep, unrolled: every adjoint from +0.0, each contribution added in
the sweep's node order (`grads[b] -= g·v_a/(v_b·v_b)` for a quotient's divisor, `grads[a] -= g` for a negation), and a
node whose adjoint is 0 skipped, as hisab 1.4.0's `backward` skips it (`if g == 0.0 { continue; }`): so A = 0 or an
underflowing K contributes nothing and never meets a partial that is ±∞ (0·∞ = NaN where (R·T)² underflows), and
∂E/∂T at Q = 1 is +0 (the textbook −R·ln 1/(nF) is −0). The tape's exp and ln are the correctly rounded `_km_cr_exp` /
`_km_cr_ln` (ADR 0006). Every Ok golden row is then the oracle's vector bit for bit (543 rows, 0 ulp; glibc rounds all
their calls correctly), compared at class D as decided.

**The mapping**, "oracle non-finite ↔ Cyrius Err" made exact: a golden row is noted with the port's Err, `[ADR 0002 #2:
the port's Err(V) (E:V)]`, wherever the validating Rust function refuses the input or an entry is not finite; the
harness requires that Err (mapped, accepted under ADR 2, counted) and compares every other row, whose entries are all
finite, at class D. 227 rows: 170 whose oracle output is not finite (107 ComputationError, 63 refused inputs) and 57
whose oracle output is finite at an input the port refuses (a negative or zero T, an infinite ΔG that makes K = 0,
n = 0 …) — equilibrium_constant_gradient 76, arrhenius_gradient 65,
nernst_gradient 66, gibbs_temperature_sensitivity 20. `validate.py` §27 checks every note against the validation rule.

**compute_gradient** keeps a `Result` (Ok(n), the n partials in the caller's buffer). #2's "returns its vector directly"
assumed Rust's failure-free call; in Cyrius the tape is the caller's (hisab 3.4.1's `ad_tape_new(cap)`, reused since the
bump allocator never frees), so a null pointer or negative count, a tape too small for the parameters, and an output
that is not a node of the tape (every hisab op returns AD_FULL once the tape is full) are Err(KM_ERR_INVALID_INPUT), and
an index table that can't be allocated (past 64 parameters) Err(KM_ERR_COMPUTATION_ERROR). Nothing Rust can return is
an Err: every input Rust accepts gives Ok. Its rows compare hisab 3.4.1's tape with 1.4.0's at class D (ADR 0003's M5b
amendment).

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
