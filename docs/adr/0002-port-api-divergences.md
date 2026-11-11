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
| 7 *(realised 2026-11-10, M4f)* | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational oxidation state (numerator / positive denominator, reduced; a value outside i32 is Rust's Err); balancing takes an integer matrix (flat row-major, Rust's \|entry\| ≤ 2⁵³ kept) | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer average oxidation state" `Err` ↔ `Ok(p/q)`, noted with the port's whole result (9,090 rows, mapped, counted) and verified by reference tests (Fe₃O₄ +8/3, C₃H₈ −8/3); balancing class A after converting the oracle's integer-valued f64 inputs as Rust's test does; the 78 rows with no integer counterpart (non-integer, non-finite, beyond i64, ragged) keep the oracle's Err(InvalidInput), noted and counted |
| 12 *(2026-11-09)* | `inorganic::cfse_octahedral` computes `t2g · (−0.4) + eg · 0.6`, which misses the exact (3·e_g − 2·t₂g)/5 on 13 of its 22 (d, spin) cases (high-spin d⁵ −2⁻⁵² and d¹⁰ −2⁻⁵¹ instead of 0, d³ −1.2000000000000002) | the exact value: the integer numerator 3·e_g − 2·t₂g converted to f64, then one correctly rounded division by 5 | class A. Where the oracle's value equals the exact one (9 cases), a plain row; elsewhere the golden row keeps the oracle's value, its note names the exact one, and the harness maps the two (q = round(5 × oracle) must be within 10⁻⁹ of it, and the port must return q / 5), reporting each such row by name as accepted under ADR 2 |
| 13 *(2026-11-09, M4c; completed in M4d)* | `nuclear::atoms_after_half_lives` computes `initial * 2.0.powf(-n)`: where 2⁻ⁿ is not a normal double (n > 1022 or n ≤ −1024) it overflows to ∞, underflows to 0 (n = 1075 is a tie, to 0) or is a subnormal of 52 bits or fewer, so the result is ±∞, 0, NaN or up to 2⁵² ulp off where N₀/2ⁿ is representable (1000 after 1075 half-lives: 0, not 2.47 × 10⁻³²¹; 2¹⁰⁰ after 1,100: 0, not 2⁻¹⁰⁰⁰) | there, N₀/2ⁿ **correctly rounded** (rounded once, subnormal results at the 2⁻¹⁰⁷⁴ grid, ties to even): since M4d `src/mathx.cyr`'s `_km_cr_exp2_scaled(N₀, −n)` for every N₀ and n, 2⁻ᶠ from exp's machinery past 53 bits with a rounding test and a triple-double fallback (M4c's version was faithful for a fractional n with N₀ not a power of two). ±0, ±∞, NaN N₀ give N₀ (scaleB). Rust's form stays elsewhere | class C; equal rows plain; the others keep the oracle's value with the port's in the note (`(f:…)`), which the harness requires bit for bit and reports as accepted under ADR 2, counted (1,810 rows) |
| 14 *(2026-11-09, M4c)* | `nuclear::is_secular_equilibrium` tests `λ_p < 0.01·λ_d`: 0.01 is not 1/100 and the product rounds or underflows, so (0.03, 3) is false and (0, 5 × 10⁻³²⁴) false | the exact λ_p < λ_d/100: 100·λ_p against λ_d by an exact double-double (`_km_two_prod` on the significands, the exponents first), Rust's NaN (false), ±∞ and ±0 semantics kept | class A; equal rows plain; the others noted `(i:…)`, mapped and reported by name (208 rows) |
| 15 *(2026-11-09, M4b's inorganic, decided in M4c)* | `inorganic::delta_tet_from_oct` computes `delta_oct * 4.0 / 9.0`, correctly rounded except where Δ·4 overflows (\|Δ\| ≥ 2¹⁰²², above DBL_MAX/4): ±∞ where 4Δ/9 is finite | there, (Δ/9)·4: one rounded division and an exact ×4 that can't overflow, which is again the correctly rounded 4Δ/9 (4·DBL_MAX/9 ≈ 7.99 × 10³⁰⁷). Rust's form elsewhere | class B; noted `(f:…)`, mapped and reported by name (59 rows) |
| 16 *(2026-11-09, M4d)* | `potential::lj_equilibrium_distance` computes `sigma * 2.0.powf(1.0 / 6.0)`: the exponent 1/6 is rounded (a simple exact constant misrepresented in binary) and the product rounds again, so r_min = σ·2^(1/6) is one ulp off on about a fifth of all σ (4,117 of 20,000 seeded). (2^RN(1/6) rounds to RN(2^(1/6)), MPFR checks; the misses are the second rounding) | the correctly rounded σ·2^(1/6), rounded once: σ = s·2^e, s·(H + M + L) by `_km_td_mul_d` (2^(1/6) as a triple-double from `mathx_tables.py`'s exact integer root), then `_km_cr_finish` at the final scale; no rounding test is needed, since 2^(1/6)'s continued fraction proves every σ·2^(1/6) at least 2^-55.1 ulp from a rounding boundary, far beyond the product's 2^-102.6. NaN, ±0, ±∞ give σ | class C (Rust's `powf`); equal rows plain; the others noted `(f:…)`, mapped and counted (1,041 rows) |
| 17 *(2026-11-09, M4d)* | `potential::lennard_jones` computes `4.0 * epsilon * (sr12 - sr6)`, i.e. (4·ε)·d: for ε ≥ 2¹⁰²² the exact 4·ε overflows, so the result is ±∞ (or ∞·0 = NaN at r = σ) where 4·ε·d is representable | there, (ε·d)·4: one rounded product and an exact ×4, which is RN(4·ε·d), the value Rust's (4·ε)·d has wherever 4·ε is finite (4·ε is exact). Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (180 rows) |
| 18 *(2026-11-09, M4d)* | `potential::morse_width_parameter` computes `omega_e * (reduced_mass / (2.0 * d_e)).sqrt()`: for D_e ≥ 2¹⁰²³ the exact 2·D_e overflows, μ/∞ = 0 and a = 0 where ω·√(μ/(2·D_e)) is representable | there, (μ·½)/D_e: RN(μ/(2·D_e)) (μ·½ is exact, or μ < 2⁻¹⁰²¹ and the quotient rounds to 0 either way), the value Rust's μ/(2·D_e) has wherever 2·D_e is finite; then Rust's √ and ×ω. Rust's form elsewhere | class B; noted `(f:…)`, mapped and counted (181 rows) |
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
| 29 *(2026-11-10, the overflow audit's B2)* | `solution::dilution` computes (M₁·V₁)/V₂, `solution_conductivity` (Λ·c)/1000, `electrochemistry::faraday_mass_deposited` (Q·M)/(n·F), `faraday_charge_required` ((m·n)·F)/M, `faraday_mass_ratio` (m_A·(M_B/n_B))/(M_A/n_A) and `nernst_potential` E° − ((R·T)/(n·F))·ln Q: the five quotients round two to four times, so about a quarter of results are 1–2 ulp off (`dilution(1.5, 0.2, 0.6)` 0.5000000000000001), and a partial product, an equivalent weight M/n or R·T that over- or underflows gives ±∞, 0, NaN or a degraded value although the value is representable (`faraday_mass_deposited(DBL_MAX, 63.546, 2)` ∞, not 5.92·10³⁰⁴; `solution_conductivity(126.4, DBL_MAX)` ∞, not 2.27·10³⁰⁷; `faraday_mass_ratio(1, 5e-324, 3, 1e-300, 1)` ∞, not 6.07·10²³; `nernst_potential(0.342, 2, DBL_MAX, 0.1)` ∞, not 1.78·10³⁰⁴, and at Q = 1 NaN, not E°) | the five monomials **rounded once** — RN(M₁·V₁/V₂), RN(Λ·c/1000), RN(Q·M/(n·F)), RN(m·n·F/M), RN(m_A·M_B·n_A/(n_B·M_A)), F as `KM_FARADAY`'s double, 1000 and the electron counts exact — by #27's `_km_cr_mono` (two new straight-line shapes, `_km_mul2_div2_fast` and `_km_mul3_div2_fast`); `nernst_potential` under **#20's semantics** with the correctly rounded ln Q (`_km_cr_ln`, ADR 0006) as its primitive: Rust's operation sequence with an unbounded exponent (`_km_ux_*`), the final subtraction rounded once into binary64 (`_km_nernst_ux`), behind a range test (`_km_nernst_fast`) that proves Rust's chain stayed normal, so those inputs keep Rust's bits (every `nernst_potential_25c` input among them). The unvalidated `solution_conductivity` keeps Rust's expression for a non-finite input; Rust's checks, Err codes and order are kept | class B (the monomials) and C (`nernst_potential`); equal rows plain; the others noted `(f:…)`, mapped and counted (94 `dilution`, 148 `solution_conductivity`, 904 `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 119 `nernst_potential`: 2,756) |
| 30 *(2026-11-10, the overflow audit's B2)* | `spectroscopy::absorbance` computes (ε·l)·c, `concentration_from_absorbance` A/(ε·l), `rigid_rotor_energy` ((h·B)·J)·(J + 1), `rotational_constant` h/(((8·π)·π)·I), `diatomic_moment_of_inertia` (μ·r)·r and `vibrational_frequency` (k/μ).sqrt(): two to four roundings of one product or quotient of inputs and constants (one quotient under a √), so a fraction of results are 1–2 ulp off (HCl's ω = √(516/1.6266·10⁻²⁷) is 563228298138879.9, not …880.0; `rotational_constant(1.45e-46)` 57875971594.410576, not …41057), and a partial product or quotient that over- or underflows gives ±∞, 0, NaN or a degraded value although the value is representable (`vibrational_frequency(500, 2⁻¹⁰²²)` ∞, not 1.50·10¹⁵⁵, and `(10⁻³⁰⁰, 10³⁰)` 0, not 10⁻¹⁶⁵; `rigid_rotor_energy(2³² − 1, 10⁻³⁰⁰)` 0, not 1.22·10⁻³¹⁴; `rotational_constant(10⁻³²⁰)` 8.392090·10²⁸⁴, not 8.392109·10²⁸⁴; `absorbance(1e300, 1e300, 0)` NaN, not 0; `concentration_from_absorbance(1e300, 1e200, 1e200)` 0, not 10⁻¹⁰⁰) | the five monomials **rounded once** — RN(ε·l·c), RN(A/(ε·l)), RN(h·B·J·(J + 1)), RN(h/(8·π·π·I)) (h `KM_PLANCK`'s double, π its double each time Rust uses it, 8, J and J + 1 exact), RN(μ·r·r) — by #27's helper (two new straight-line shapes, `_km_mul4_fast` and `_km_div3_fast`); `vibrational_frequency` the **correctly rounded √ of the exact k/μ**, `src/mathx.cyr`'s `_km_cr_sqrt_div`: a double-double fast path (proven within 5.7u², u = 2⁻⁵³) with `_km_dd_div_cr`'s rounding test, and an exact integer square root on `_km_bn_*` (√(k/μ) is never a rounding boundary of binary64, proven), ±∞ only past DBL_MAX. A non-finite unvalidated input (absorbance's, the moment of inertia's, the rotor's B, concentration_from_absorbance's A) keeps Rust's expression; Rust's checks, Err codes and order are kept | class B; equal rows plain; the others noted `(f:…)`, mapped and counted (195 `absorbance`, 219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant`, 118 `diatomic_moment_of_inertia`: 2,106) |

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
