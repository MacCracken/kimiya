# kimiya — Current State

> Refreshed every release and at every port milestone. `CLAUDE.md` is
> preferences/process/procedures (durable); this file is **state** (volatile).

## Version

**1.2.0** (`VERSION`) — the last Rust release, tag `1.2.0` (`077e2e2`), frozen as the oracle in
`rust-old/`. The Cyrius port is in progress and ships as **2.0.0**; `VERSION` moves to 2.0.0 at M7.
`cyrius.cyml` reads the version through `version = "${file:VERSION}"`.

## Port progress

Plan of record: [`kimiya-port-plan.md`](kimiya-port-plan.md).

| Milestone | Status |
|---|---|
| M0 Initialize | ✅ **done 2026-11-09.** Steps 0–4: oracle frozen green at tag `1.2.0` (662 + 24 + 2 Rust tests, rustc 1.99.1, glibc 2.44; `rust-old/README.md` restored from the tag for its doctest); `cyrius port` scaffold; CLAUDE.md merged, roadmap restored, Cyrius CI extended, Rust CI/scripts/overview under `rust-old/`, `Cargo.lock` tracked; `cyrius.cyml` with hisab 3.4.1 (`d6703c5`); math semantics probed. **Step 5**: `scripts/oracle/` writes `tests/golden/` (constants 45 + math 18,647 rows, 2.0 MB), deterministic, platform-pinned, `--check`ed, independently validated (`validate.py`, 60-digit reference). **Step 6**: `src/constants.cyr` (14 `KM_` constants) generated and compiler-probed. **Step 7**: `tests/parity.tcyr` + `parity_core.cyr` + `ulp.cyr`. **Step 8**: [`port-parity.md`](port-parity.md) — K = 2 ulp for the math primitives on x86_64 ELF, x86_64 PE and aarch64; 13,340 ok, 0 failed, 5,351 pending, 1 accepted on all three. **Step 9**: ADR 0003 (transcribe the hisab 1.4.0 routines hisab 3.4.1 lacks); ADR 0004 (ganita `pow(-1, DBL_MAX)`) |
| M1 Foundation | ✅ **done 2026-11-09.** `src/error.cyr` (`KmError`), `src/mathx.cyr` (the `_km_` math helpers; `_km_log10` / `_km_expm1` from Rust libm 0.2.16), `src/validate.cyr`, generated `src/element_data.cyr` + `src/element.cyr` (handles, 0 for none), `src/molecule.cyr` (`KmMolecule`); golden `element` 975 + `molecule` 2,472 rows; **22,107 ok, 0 failed, 31 pending, 1 accepted** on x86_64 ELF, x86_64 PE and aarch64 (math 0 pending; log10 / exp_m1 K = 2); 38 Rust tests ported; `[lib]` + `dist/kimiya.cyr`; surface parity 55/55; benches vs Rust in [`benchmarks-rust-v-cyrius.md`](../benchmarks-rust-v-cyrius.md) |
| M2 Thermochemistry & gases | ✅ **done 2026-11-09.** M2a: `src/reaction.cyr`, `src/thermo.cyr`, generated `src/phase_data.cyr` + `src/phase.cyr` (Antoine → `Result`, ADR 0002 #4), generated `src/gas_data.cyr` + `src/gas.cyr` (`KmPengRobinsonParams`, the analytic PR cubic); correctly rounded acos/cos/ln/x^2.5 in `src/mathx.cyr` (ADR 0005). M2b: generated `src/thermochem_data.cyr` (`KmShomate`, 76 THERMOCHEM_DATA rows, 34 Shomate segments) + `src/thermochem.cyr` (reaction sums in J/mol, ADR 0002 #1; van 't Hoff; Shomate Cp and enthalpy; the flame temperature with hisab 1.4.0's Newton + bisection transcribed, ADR 0003 amendment). Golden `reaction` 1,188 + `thermo` 1,107 + `thermochem` 9,487 + `phase` 3,276 + `gas` 19,717 rows (math +7,500); **63,953 ok, 0 failed, 0 measured, 20 pending, 441 accepted** on x86_64 ELF, x86_64 PE and aarch64 (class B 0 ulp, the flame rows too; class C K = 2–10); 152 Rust tests ported; surface parity 140/140; 9 benches vs Rust |
| M3 Solutions & electrochemistry | ✅ **done 2026-11-09.** M3a: `src/solution.cyr` (concentration units, pH, Henderson–Hasselbalch, the weak acid / base charge-balance Newton transcribed as `_km_charge_balance_root`, ionic strength over ion records, Debye–Hückel, Davies, Kohlrausch) with its ten constants generated into `src/constants.cyr`; generated `src/electrochemistry_data.cyr` (49 STANDARD_POTENTIALS rows, the first UTF-8 string blob) + `src/electrochemistry.cyr` (handles, Nernst, cell potentials, ΔG, Faraday's laws). M3b: `src/kinetics.cyr` (Arrhenius, the zero / first / second / nth-order laws and half-lives, steady state and pre-equilibrium, Michaelis–Menten, Lineweaver–Burk into the caller's two-f64 buffer, collision theory, Eyring). Golden `solution` 28,566 + `electrochemistry` 18,574 + `kinetics` 30,146 rows; at M3b 119,044 ok, 0 failed, 22,205 measured (the six kinetics class-C functions, gated since M3c), 10 pending, 441 accepted on all three targets; 131 Rust tests ported; surface parity 203/203; 7 benches vs Rust. Its primitives moved to the correctly rounded ones in M3c |
| M3c Correct rounding (ADR 0006) | ✅ **done 2026-11-09.** M3c-1: the corrected oracle — `scripts/oracle/mpfr_ref.py` (MPFR 4.2 through ctypes), `validate.py --write-corrections` → `tests/golden/<module>.cr.tsv` (665 rows), the harness compares corrected rows against them. M3c-2: correctly rounded `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`, `_km_cr_pow` in `src/mathx.cyr` (two-phase, Ziv), tables generated into `src/mathx_data.cyr` by `scripts/oracle/mathx_tables.py`, verified by `scripts/oracle/cr_verify.py` (≥ 1,000,000 arguments per primitive against MPFR, 0 mismatches, the same bits on every target) and `tests/mathx.tcyr`. M3c-3: every module on them in Rust's operation order (`_km_acos` / `_km_cos` kept; x^2.5 is `_km_cr_pow`); `_km_log10`, `_km_expm1`, `_km_log10_ln`, the `_km_ln` wrapper, `_km_pow_2_5` and `_km_dd_div` removed; every class-C function at bound 0; glibc-matching is a failure; 2,486 hard rows (subnormal exp × prefactor) added to `reaction`, `phase`, `thermochem`, `gas`. **143,071 ok, 0 failed, 0 measured, 10 pending, 1,105 accepted (ADR 6: 665; ADR 0002 #4: 440), 0 glibc-matching** on x86_64 ELF, x86_64 PE and aarch64, every class-C row 0 ulp (sin, atan2 and the 247 cos rows past \|x\| = 8 on the stdlib at 2 ulp); ADR 0004's row a plain ok; 15 test files, 0 warnings; surface parity 203/203 ([`port-parity.md`](port-parity.md) §9–§10) |
| M4 Structure, nuclear & spectroscopy | ✅ **done 2026-11-10.** **M4b ✅ done 2026-11-09** — generated `src/organic_data.cyr` (`KmBondType`, `KmGeometry`, `KmFunctionalGroup` with their property columns: bond energies, ideal angles, the UTF-8 pattern blob, the polar / donor flags and ADR 0002 #6's new acceptor column: IUPAC 2011's criterion, weak acceptors included, every group but the alkane, one source per family) + `src/organic.cyr` (`bond_type_energy_kj`, `enthalpy_from_bonds` over bond records, `predict_geometry` → `Option`, `geometry_ideal_bond_angle`, the `functional_group_*` properties; ADR 0002 #6: `KM_GEOMETRY_T_SHAPED`, `functional_group_is_hydrogen_bond_donor`, `…_acceptor`); generated `src/inorganic_data.cyr` (`KmLigand`, the private octahedral occupancy and unpaired-electron tables read through the API, the 10 Madelung / Born constants, the 30 Shannon IONIC_RADII rows) + `src/inorganic.cyr` (CFSE exact by ADR 0002 #12, (3·e_g − 2·t₂g)/5 in one correctly rounded division, Δ_tet, the spectrochemical series, spin-only moments, Born–Landé, Born–Haber, `lookup_ionic_radius` handles). Golden `organic` 9,104 + `inorganic` 15,107 rows and the last 10 `constants` rows: **167,279 ok, 0 failed, 0 measured, 0 pending, 1,118 accepted (ADR 6: 665; ADR 2: 440 Antoine + 13 exact-CFSE rows, by name), 0 glibc-matching** on x86_64 ELF, x86_64 PE and aarch64, every M4b row 0 ulp from its reference; 51 Rust tests ported (+ 2 integration); surface parity 323/323. **M4c ✅ done 2026-11-09** — generated `src/nuclear_data.cyr` (the 18 ISOTOPES rows: label blob, `u8` Z, `u16` A, half-life, a UTF-8 decay-mode blob) + `src/nuclear.cyr` (decay constants and laws, mass defect and binding energies, the SEMF → `Result` by ADR 0002 #3, Q-values over slices, the two-step decay chain into the caller's buffer with Rust's equal-rate limit, the secular-equilibrium test, `lookup_isotope` handles), every transcendental correctly rounded; by the user's decision on the review, ADR 0002 #13 (N₀/2ⁿ rounded once where Rust's 2⁻ⁿ is not a normal double: correctly rounded for an integer n or N₀ = ±2^j, faithful otherwise), #14 (the exact λ_p < λ_d/100) and, for M4b's inorganic, #15 (the correctly rounded 4Δ/9 past Δ·4's overflow). Golden `nuclear` 54,057 rows (the SEMF over every Z ≤ A ≤ 3Z with A ≤ 300) + 23 corrected rows, `inorganic` 15,162: **219,266 ok, 0 failed, 0 measured, 0 pending, 3,243 accepted (ADR 6: 688; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,812 N₀/2ⁿ + 208 secular + 59 Δ_tet, all but the Antoine ones by name), 0 glibc-matching** on all three targets, every row 0 ulp from its reference; 36 Rust tests ported; surface parity 347/347. **M4a ✅ done 2026-11-09** — correctly rounded `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`, `_km_cr_atan2`, `_km_cr_acos` (two-phase; Cody–Waite and Payne–Hanek reduction, generated tables, `_km_td_norm` for the accurate division and square root), verified against MPFR by `cr_verify.py` on about 1.0–1.1 million arguments each (100,000 through the accurate path) on all three targets. **M4d ✅ done 2026-11-09** — `src/potential.cyr` (Lennard-Jones, Morse, Coulomb, distances, bond and dihedral angles through three-f64 point pointers, on `_km_cr_exp` and `_km_cr_atan2`); gas.cyr on M4a's `_km_cr_acos` / `_km_cr_cos` (every gas row keeps its bits) and ADR 0005's double-double `_km_cos` / `_km_acos` removed; the `math` trig rows at bound 0 (no golden transcendental row on the stdlib); ADR 0002 #13 correctly rounded (`_km_cr_exp2_scaled`), the new #16 (σ·2^(1/6) correctly rounded), #17 (Lennard-Jones past 4·ε's overflow), #18 (the Morse width past 2·D_e's overflow). Golden `potential` 26,938 rows + 1 corrected: **244,803 ok, 0 failed, 0 measured, 0 pending, 4,644 accepted (ADR 6: 689; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,810 N₀/2ⁿ + 208 secular + 59 Δ_tet + 1,041 r_min + 180 Lennard-Jones + 181 Morse width), 0 glibc-matching** on all three targets, every row 0 ulp from its reference; 33 Rust tests ported; surface parity 357/357. **M4e ✅ done 2026-11-09** — generated `src/spectroscopy_data.cyr` (`KmSpectralSeries` and its n_lower column, the isotope abundance table from the oracle's `isotope_abundances(z)` at every u8 z, Rust's private R_H read through the API) + `src/spectroscopy.cyr` (Beer–Lambert on `_km_cr_exp10` / `_km_cr_log10`, the photon conversions, the Bohr levels and the Rydberg formula, the spectral series, the oscillator and rotor, `isotope_pattern`'s pruned convolution with Rust's stable sort, `fft_spectrum` on hisab 1.4.0's `num::fft` transcribed with correctly rounded twiddles (ADR 0003 amendment), `find_peaks`; Vec results into caller buffers); ADR 0002 #19 (the wavenumbers' correctly rounded 1/(100·x), `_km_cr_hundredth_div`). Golden `spectroscopy` 52,985 rows + 118 corrected: **296,016 ok, 0 failed, 0 measured, 0 pending, 6,416 accepted (ADR 6: 807; ADR 2: … + 1,654 1/(100x)), 0 glibc-matching** on all three targets, every row 0 ulp from its reference; 62 Rust tests ported (+ 2 integration); surface parity 394/394. **M4e-2 ✅ done 2026-11-10** — ADR 0002 #20 (the user's decision, "Yes, prescale"): `distance_3d`, `bond_angle` and `dihedral_angle` evaluate Rust's sequence with an unbounded exponent range (the final √ / atan2 rounded once), Rust's chain kept as a fast path behind a proven range test, `src/mathx.cyr`'s general `_km_ux_*` helpers as the slow path (and the machinery of the library-wide audit after M4). Golden `potential` 31,761 rows (4,823 new edge rows, 3,605 mapped under #20) + 3 corrected: **297,232 ok, 0 failed, 0 measured, 0 pending, 10,023 accepted (ADR 6: 809; ADR 2: … + 3,605 unbounded-exponent geometry rows), 0 glibc-matching** on all three targets; no row whose Rust chain stays in the normal range changes. **M4e-3 ✅ done 2026-11-10** — ADR 0002 #21 and #22 (the user's decisions, "Compute exactly" and "Round once"): `bohr_energy_level`, `bohr_transition_energy`, `rydberg_wavelength` and `spectral_series_wavelength` return their exact rationals rounded once, `reduced_mass` the correctly rounded m₁·m₂/(m₁ + m₂) for every pair of positive finite masses; double-double fast paths with a rounding test (`_km_dd_mul`, `_km_dd_div_cr`) and an exact slow path on `src/mathx.cyr`'s new multi-word integers and `_km_cr_ratio(p, q, e)` (RN(p/q·2^e), for reuse by M4f and the overflow audit). Golden `spectroscopy` 66,691 rows (13,706 new: every exact Bohr level, the 13 transition ties of an exhaustive search and its exact-value families, n₁ = n₂ − 1 at every scale, the reduced mass's tie families, gaps and range ends): **301,484 ok, 0 failed, 0 measured, 0 pending, 19,477 accepted (ADR 6: 809; ADR 2: … + 6,617 exact Rydberg rationals + 2,837 exact reduced masses), 0 glibc-matching** on all three targets. **M4e-4 ✅ done 2026-11-10** — ADR 0002 #23 (the user's decision, "Upgrade the data"): the isotope table is the cited AME2020 masses at full precision and the TICE 2009 abundances with ³⁶S, from `scripts/oracle/isotopes.tsv` (each mass with its verbatim AME2020 line; the generator checks it against the oracle's table: the same isotopes plus ³⁶S, Rust's abundances, Rust's masses these to five decimals, each element summing to 1); `isotope_pattern`'s algorithm unchanged. Golden `spectroscopy` keeps the oracle's values, 354 rows (7 `isotope_abundances`, 347 `isotope_pattern`) noted with the port's whole result and mapped: **301,130 ok, 0 failed, 0 measured, 0 pending, 19,831 accepted (ADR 6: 809; ADR 2: … + 354 #23 rows), 0 glibc-matching** on all three targets. **M4f ✅ done 2026-11-10** — `src/stoichiometry.cyr` with ADR 0002 #7 realised: `balance_equation` on an integer matrix (flat row-major, Rust's \|entry\| ≤ 2⁵³ kept; Rust's fraction-free Gauss–Jordan elimination with checked i64 arithmetic and the i128 verification on a signed three-word sum, every overflow an Err at Rust's step), `assign_oxidation_states` with exact rational states ((Z, numerator, denominator) records; Rust's non-integer-average Err is the average, an average outside i32 stays Rust's Err), `limiting_reagent`, `theoretical_yield_moles`, `percent_yield`, `percent_composition`, `empirical_formula` operation for operation. Golden `stoichiometry` 41,435 rows (9,090 oxidation rows noted with the port's exact states, 78 balancing rows with no integer counterpart, both mapped under ADR 2): **333,397 ok, 0 failed, 0 measured, 0 pending, 28,999 accepted (ADR 6: 809; ADR 2: … + 9,168 #7 rows), 0 glibc-matching** on all three targets, every float row 0 ulp; 63 Rust tests ported (7 adapted for #7) + 6 Cyrius-only groups; surface parity 401/401. The five float functions' measured misses of the exact values are reported, not fixed (ADR 0002's M4f addition). **M4f-2 ✅ done 2026-11-10** — ADR 0002 #24, #25 and #26 (the user's decisions, "Round once", "Compare exactly", "Evaluate exactly"): `percent_yield`, `theoretical_yield_moles` and `percent_composition` return their exact rationals rounded once (`_km_cr_mul_div`: a double-double fast path with a rounding test and an exact midpoint resolution, `_km_cr_ratio` past the normal range), `limiting_reagent` compares its ratios exactly (the rounded ones where they differ, `_km_ratio_lt` where they tie), `empirical_formula` evaluates its moles, ratios, 1/10 tolerance and counts exactly (a proven 2³² − ½ cut, a fast tolerance test that defers to an exact integer one at the knife edge); Rust's checks, codes and order kept. Golden `stoichiometry` 47,776 rows (6,341 hard rows appended), the oracle's values kept, 3,154 noted with the port's result and mapped (#24: 1,126 `percent_yield`, 240 `theoretical_yield_moles`, 598 `percent_composition`; #25: 920; #26: 270): **336,584 ok, 0 failed, 0 measured, 0 pending, 32,153 accepted (ADR 6: 809; ADR 2: … + 3,154 #24–#26 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 2.23 million arguments against `Fraction`, 0 mismatches |
| Overflow audit (after M4) | **B1 ✅ done 2026-11-10** — ADR 0002 #27 and #28 (the user's decisions "Library-wide audit" and "Round once"): the survey's F1 (`heat_transfer`, `specific_heat_from_calorimetry`, `boiling_point_elevation`, `freezing_point_depression`, `osmotic_pressure` rounded once) and F4 (`zero_order_half_life`, `steady_state_intermediate`, `pre_equilibrium_rate` rounded once; `second_order_concentration`, `michaelis_menten` with an unbounded exponent behind a range test, Rust's bits in range); M4f-2's `_km_cr_mul_div` generalised into `src/mathx.cyr`'s correctly rounded monomial `_km_cr_mono` (the fast path takes subnormal operands and resolves midpoints by an exact expansion sign; `_km_cr_ratio` for possibly subnormal results). Golden `thermo` 1,644, `phase` 4,174, `kinetics` 31,210 rows (2,129 edge rows appended), the oracle's values kept, 1,318 noted and mapped (#27: 682; #28: 636): **337,395 ok, 0 failed, 0 measured, 0 pending, 33,471 accepted (ADR 6: 809; ADR 2: … + 1,318 #27–#28 rows), 0 glibc-matching** on all three targets. **B2 ✅ done 2026-11-10** — ADR 0002 #29 and #30 (B1's decisions unchanged): the survey's F3 (`dilution`, `solution_conductivity`, `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` rounded once; `nernst_potential` with an unbounded exponent behind a range test, the correctly rounded ln Q its primitive, Rust's bits in range) and F5 (`absorbance`, `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia` rounded once; `vibrational_frequency` the correctly rounded √ of the exact k/μ); four new monomial shapes and `src/mathx.cyr`'s correctly rounded √ of a quotient `_km_cr_sqrt_div` (a proven double-double fast path, an exact integer root on `_km_bn_*`). Golden `solution` 29,024, `electrochemistry` 19,397, `spectroscopy` 67,869 rows (2,459 edge rows appended), the oracle's values kept, 4,862 noted and mapped (#29: 2,756; #30: 2,106): **334,992 ok, 0 failed, 0 measured, 0 pending, 38,333 accepted (ADR 6: 809; ADR 2: … + 4,862 #29–#30 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 3.09 million arguments for B2, 0 mismatches. **B3 ✅ done 2026-11-11** — ADR 0002 #31 and #32 (B1's decisions unchanged): the survey's F2 (`ideal_gas_pressure`, `ideal_gas_volume`, `gas_density` rounded once; `van_der_waals_pressure`, `peng_robinson_pressure`, `redlich_kwong_pressure` with an unbounded exponent behind range tests, the co-volume tests on the unbounded co-volume, Tc^2.5 exactly rescaled) and F6 (`coulomb`, `coulomb_force` rounded once; `morse_width_parameter` the correctly rounded ω·√(μ/(2·D_e)); `lennard_jones`, `born_lande_lattice_energy`, `binding_energy_per_nucleon` with an unbounded exponent, Rust's bits in range), M4d's #17 and #18 folded into #32 (#17's values kept, #18's correctly rounded); `src/mathx.cyr`'s correctly rounded w·√((x/y)·2^e), `_km_cr_mul_sqrt_ratio` (B2's root times w with a rounding test; an exact integer root that rounds its midpoints to even), and `_km_ux_normal`. Golden `gas` 22,773, `potential` 33,464, `inorganic` 15,416, `nuclear` 54,224 rows (3,804 edge rows appended), the oracle's values kept, 5,339 noted and mapped (#31: 1,157, 98 of them the port's Err; #32: 4,182, #17's 180 and #18's 181 among them): **333,818 ok, 0 failed, 0 measured, 0 pending, 43,311 accepted (ADR 6: 809; ADR 2: … + 5,339 #31–#32 rows − 361 former #17 / #18), 0 glibc-matching** on all three targets; `cr_verify.py` 2.82 million arguments for B3, 0 mismatches. **B4 ✅ done 2026-11-11** — ADR 0002 #33 and #34 (B1's decisions unchanged): the survey's F7 (`hess_law`, `enthalpy_from_formation`, the three thermochem reaction sums, `q_value`, `ionic_strength`, `born_haber_lattice_energy`, `gibbs_free_energy`, `kohlrausch_molar_conductivity`, `limiting_molar_conductivity`, `zero_order_concentration`) and F8 (`Shomate::cp`, `enthalpy_antiderivative`), all under #20's semantics in Rust's summation order behind fast tests (a pure sum's "result finite": a sum of doubles is exact below 2⁻¹⁰²¹; the Shomate box, which every tabulated segment is in); #29's fast path corrected (its 2⁻¹⁰²² bound admitted a product IEEE rounds up to 2⁻¹⁰²²; now 2⁻¹⁰²¹, B1–B3's other fast paths audited clean); the audit's benches moved to `tests/audit.bcyr` (the user's decision). Golden `thermo` 1,828, `reaction` 1,730, `nuclear` 54,312, `solution` 29,454, `inorganic` 15,586, `kinetics` 31,381, `thermochem` 11,149, `electrochemistry` 19,420 rows (2,530 edge rows appended), the oracle's values kept, 924 noted and mapped (#33: 697; #34: 210; #29: 17 more): **335,424 ok, 0 failed, 0 measured, 0 pending, 44,235 accepted (ADR 6: 809; ADR 2: … + 924 #33, #34 and #29 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 1.25 million arguments for B4, 0 mismatches. **B5 ✅ done 2026-11-11** — ADR 0002 #35 (B1's decisions unchanged): the survey's F9, the exp / exp10 *arguments* — `equilibrium_constant`'s −ΔG/(R·T) and the activity models' log₁₀ γ (`debye_huckel_limiting`, `debye_huckel_extended`, `davies_activity`) — under #20's semantics, the argument's final quotient or product rounded once into binary64, then the unchanged correctly rounded exp / 10^x (no new primitive: the ADR's lemma proves exp / 10^x of the binary64 argument the correctly rounded exp / 10^x of the 53-bit unbounded one), behind inline fast tests (R·T in [2⁻¹⁰²¹, DBL_MAX]; A·z² finite; the extended law's numerator and (B·å)·√I finite), `equilibrium_constant`'s ComputationError judged on the port's K. Golden `reaction` 2,406, `solution` 30,644 rows (1,866 edge rows appended; none of the 2,627 earlier rows changes), the oracle's values kept, 845 noted and mapped (326 `equilibrium_constant`, 27 of them the port's Err and 33 its K where Rust's overflows; 55 `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity`): **336,445 ok, 0 failed, 0 measured, 0 pending, 45,080 accepted (ADR 6: 809; ADR 2: … + 845 #35 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 505,982 arguments for B5, 0 mismatches. **B6 ✅ done 2026-11-11** — ADR 0002 #36, #37 and #38 (B1's decisions unchanged; the user's "FFT: both" and "fix all three iterative edges", `find_peaks`' exact threshold settled under the policy): the survey's F13 (`fft_spectrum` under #20's semantics over the whole transform — hisab 1.4.0's butterflies and twiddles with an unbounded exponent, each magnitude's √ rounded once — behind a per-operation operand test on Rust's own transform, a bin's unbounded √ where only its square sum leaves the range, 1.2's IEEE spectrum for a NaN or ±∞ sample), F14 (`find_peaks`' threshold test exact, Rust's rounded threshold deciding behind a proven bound, `src/mathx.cyr`'s new `_km_sign_scaled` for the exact sign) and F16's flame part (`adiabatic_flame_temperature`'s bracket and its transcribed bisection compare signs exactly; Newton's absolute 10⁻³ J/mol tolerance, which accepts the initial guess at tiny scales, kept). Golden `spectroscopy` 70,777, `thermochem` 11,592 rows (3,351 edge rows appended), the oracle's values kept, 913 noted and mapped (418 `fft_spectrum`, 41 of the 213 earlier rows; 222 `find_peaks`, 3 of 994; 273 `adiabatic_flame_temperature`, 16 of 1,937): **338,883 ok, 0 failed, 0 measured, 0 pending, 45,993 accepted (ADR 6: 809; ADR 2: … + 913 #36–#38 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 644,838 arguments for B6, 0 mismatches. **B7 ✅ done 2026-11-12** — ADR 0002 #39 (B1's decisions unchanged; "Balancing: multi-word integers, no size limit", the survey's remedy (a)): the survey's F15, `balance_equation` on unbounded integers — Rust's i64 elimination first, unchanged, and at its first overflow the same steps from the original matrix on signed multi-word integers (`_km_sb_*`: 31-bit limbs, Knuth's division, Lehmer's gcd), the tables sized from a proven Hadamard bound (`_km_balance_bits`); every overflow Err gone, every other Err and Rust's order kept. Golden `stoichiometry` 48,297 rows (521 edge rows appended), the oracle's values kept, 413 noted and mapped (104 of the 3,329 earlier `balance_equation` rows, exactly the survey's; 309 edge rows): **338,991 ok, 0 failed, 0 measured, 0 pending, 46,406 accepted (ADR 6: 809; ADR 2: … + 413 #39 rows), 0 glibc-matching** on all three targets; `cr_verify.py` 54,306 matrices for B7, 0 mismatches. **B10a ✅ done 2026-11-12** (the audit's eighth batch, the first part of the survey's B10; the user's "new unbounded-exponent primitives … Rust's chain with an unbounded exponent, last step rounded once", the semantics settled under the policy, not asked) — `src/mathx.cyr`'s `_km_cr_exp_ux`, `_km_cr_expm1_ux` and `_km_cr_exp10_ux`: e^x, e^x − 1 and 10^x of an unbounded-exponent pair, the significand correctly rounded to 53 bits at any exponent, exact up to MPFR's widest exponent range and saturating past it as MPFR does (ADR 0002's "Addition, 2026-11-12 (audit B10a)", for #40); binary64's range is the binary64 functions, the cores' range without their final scaling, and beyond a new reduction (double-double for the fast path below |x·512/ln 2| = 2⁵¹, integer Payne–Hanek-style on 384 generated bits otherwise); no golden row changes (the parity line is B7's); `cr_verify.py` 448,188 arguments for B10a against MPFR at its widest range, 0 mismatches. **Paused after B10a** by the user's decision of 2026-11-12: the audit is post-port work and resumes after M5 / M7 over the whole library, M5's modules included (the roadmap lists M5a's candidate spots). The rest of the survey's B10 (#40's functions, three sub-batches) and its B11–B12 ☐ (its families F10's functions, F11, F12 and F16's other two; its B3, F5, was done in B2, its B4 and B5, F6 and F2, in B3, its B6, F7 and F8, in B4, its B7, F9, in B5, its B8, F13, F14 and the flame, in B6, its B9, F15, in B7) |
| M5 Numerics | **M5a ✅ done 2026-11-12** — `src/numerics.cyr` (hisab 1.4.0's `rk4_trajectory`, `dopri45`, `gaussian_elimination` and `euler_maruyama` transcribed as private `_km_*` routines, ADR 0003's M5a amendment; PCG32 is hisab 3.4.1's `pcg32_new` / `pcg32_next`, equal to 1.4.0's stream 1 on every golden row) and `src/reaction_dynamics.cyr` (the eight `simulate_*`, Rust's checks, codes and order; closures as function pointer + context, trajectories in the caller's flat buffer, the flat Jacobian); golden modules `numerics` (243, replayed by `tests/numerics.tcyr`; 3 ADR 0006 corrections) and `reaction_dynamics` (954; 1 correction): every row bit-identical, RK4 / backward Euler / the elimination class B, dopri45 / Euler–Maruyama / the normal draws class C at bound 0; the parity line 339,944 ok, 0 failed, 0 measured, 0 pending, 46,407 accepted (ADR 6: 810), 0 glibc-matching on x86_64 ELF, x86_64 PE and aarch64; 26 Rust tests + 2 integration tests ported (one adapted), Cyrius-only groups; `tests/numerics.bcyr`; surface parity 410/410. M5b (network, fitting, sensitivity) ☐ |
| M7 Release 2.0.0 | ☐ |

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml` `[package].cyrius`)
- **ganita**: 1.2.15, from the 6.7.18 stdlib bundle (`lib/ganita.cyr`)
- **hisab**: 3.4.1 (`[deps.hisab]`, commit-pinned in `cyrius.lock`), with its sakshi 2.5.10
- **Rust oracle**: kimiya 1.2.0 at `rust-old/`, green on rustc 1.99.1; its transcendental results come
  from the host **glibc 2.44** (Rust std)

## Source

- **Rust reference:** `rust-old/` (16,581 lines incl. tests), frozen; do not edit. `rust-old/docs/`,
  `rust-old/.github/` and `rust-old/scripts/` keep the 1.2.0 architecture overview, CI and cargo scripts;
  `rust-old/benchmarks.md` is the Rust benchmark baseline.
- **Cyrius port:** `[lib]` (bundle `dist/kimiya.cyr`, `cyrius distlib`) = `src/error.cyr`, `src/constants.cyr`
  (generated, 25 `KM_` constants: the 14 physical ones, `KM_WATER_SPECIFIC_HEAT` and the ten solution constants — K_w,
  Debye–Hückel A / B, the seven Λ°), `src/mathx_data.cyr` (generated by `scripts/oracle/mathx_tables.py` from MPFR: the
  correctly rounded primitives' tables and constants), `src/mathx.cyr` (the Rust-semantics helpers and the correctly
  rounded primitives, ADR 0006: the two-phase `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`,
  `_km_cr_pow` (M3c-2), `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`, `_km_cr_atan2`, `_km_cr_acos` (M4a) and
  `_km_cr_exp2_scaled` (M4d) and `_km_cr_hundredth_div` (M4e, ADR 0002 #19), each with `_acc` / `_fast` test entry
  points, and `_km_cr_sixth_root2_mul` (M4d, one proven rounding); since M4e-2 binary64 with an unbounded exponent
  (`_km_ux_*`: (significand, exponent) pairs, + − × ÷ √, cmp / max, the round-once `_km_ux_round` and the final
  `_km_ux_{add,sub,mul,div,sqrt,atan2}_f64`, ADR 0002 #20); since M4e-3 multi-word unsigned integers (`_km_bn_*`, 8 limbs
  of 32 bits), `_km_cr_ratio(p, q, e)` = RN(p/q·2^e) and the double-double `_km_dd_mul` / `_km_dd_div_cr` (ADR 0002 #21,
  #22); since the overflow audit's B1 the correctly rounded monomial `_km_cr_mono` (RN(Π num / Π den), with
  `_km_dd_mul_d`, `_km_dd_div_ends`, `_km_mono_tie`, `_km_exp_sign` and the shapes `_km_cr_mul_div` (moved from
  stoichiometry), `_km_cr_mul3`, `_km_cr_div2`, `_km_cr_mul3_div`, `_km_cr_mul5`, each with a straight-line fast path;
  ADR 0002 #27, #28), since B2 the shapes `_km_cr_mul2_div2`, `_km_cr_mul3_div2`, `_km_cr_mul4`, `_km_cr_div3` and the
  correctly rounded √ of a quotient `_km_cr_sqrt_div` (`_km_sqrt_div_fast`, `_km_sqrt_div_acc`; ADR 0002 #29, #30), since
  B3 the correctly rounded w·√((x/y)·2^e) `_km_cr_mul_sqrt_ratio` (`_km_mul_sqrt_ratio_fast`, `_km_mul_sqrt_ratio_acc`)
  and the range test `_km_ux_normal` (ADR 0002 #31, #32), since B10a e^x, e^x − 1 and 10^x with an unbounded exponent
  `_km_cr_exp_ux`, `_km_cr_expm1_ux`, `_km_cr_exp10_ux` (with `_fast` / `_acc`; the reductions `_km_expux_dd` and
  `_km_expux_reduce`, the cores' tails `_km_exp_core_r` and `_km_exp_acc_r`; for ADR 0002 #40); ADR 0005's double-double
  `_km_acos` / `_km_cos` were removed in M4d), `src/validate.cyr`, `src/element_data.cyr`
  (generated: 118 elements, 32 property rows, `KmElementCategory`), `src/element.cyr`, `src/molecule.cyr`,
  `src/reaction.cyr` (the overflow audit's B4: #33's `gibbs_free_energy`; B5: ADR 0002 #35's unbounded-exponent exp
  argument of `equilibrium_constant`), `src/thermo.cyr`, `src/thermochem_data.cyr` (generated: `KmShomate`, 76 THERMOCHEM_DATA rows
  with formula / name blobs, 34 SHOMATE_DATA segments), `src/thermochem.cyr`, `src/phase_data.cyr` (generated:
  `KmAntoineCoeffs`, 8 Antoine sets, the
  mmHg → Pa factor, K_b/K_f), `src/phase.cyr`, `src/gas_data.cyr` (generated: `KmPengRobinsonParams`, STP, the 5
  `KM_PR_*` presets), `src/gas.cyr` (the overflow audit's B3: ADR 0002 #31's rounded-once ideal gas law and density and
  unbounded-exponent equations of state), `src/solution.cyr` (B5: ADR 0002 #35's unbounded-exponent 10^x arguments of the
  three activity models), `src/electrochemistry_data.cyr` (generated: the 49
  STANDARD_POTENTIALS rows — UTF-8 equation and couple blobs, electrons, E° — and `KM_STANDARD_POTENTIALS_COUNT`),
  `src/electrochemistry.cyr`, `src/kinetics.cyr`, `src/organic_data.cyr` (generated, M4b: `KmBondType`, `KmGeometry`,
  `KmFunctionalGroup` and their property columns, ADR 0002 #6's acceptor column with its sources), `src/organic.cyr`,
  `src/inorganic_data.cyr` (generated, M4b: `KmLigand` and its field strengths, Rust's private octahedral d-electron
  tables, the 10 `KM_MADELUNG_*` / `KM_BORN_EXP_*`, the 30 IONIC_RADII rows), `src/inorganic.cyr`,
  `src/potential.cyr` (M4d: Lennard-Jones, Morse, Coulomb and the 3D geometry; no data file; M4e-2: the geometry's
  unbounded-exponent slow paths, ADR 0002 #20; B3: #32's unbounded-exponent Lennard-Jones, rounded-once Coulomb pair and
  correctly rounded Morse width), `src/nuclear_data.cyr`
  (generated, M4c: the 18 ISOTOPES rows, a UTF-8 decay-mode blob, `KM_ISOTOPES_COUNT`),
  `src/nuclear.cyr` (its `q_value` sums with thermo's `_km_sum`), `src/spectroscopy_data.cyr` (generated, M4e:
  `KmSpectralSeries` and its n_lower column, the isotope abundance pairs of H, C, N, O, S, Cl, Br with per-z start /
  count columns — since M4e-4 ADR 0002 #23's cited table from `scripts/oracle/isotopes.tsv`, AME2020 masses and TICE 2009
  abundances with ³⁶S —, `_KM_RYDBERG_HYDROGEN` and the mass ratio), `src/spectroscopy.cyr` (M4e: calls element's lookups and
  the earlier modules' `_km_is_enum_value`, `_km_is_u8`, `_km_is_u32`; includes `lib/alloc.cyr`; M4e-3: the Bohr model,
  the Rydberg formula and the reduced mass rounded once, ADR 0002 #21, #22), `src/stoichiometry.cyr` (M4f, between kinetics
  and organic_data: balancing on an integer matrix and rational oxidation states, ADR 0002 #7; calls element's lookups and
  `KmMolecule`; its checked i64 helpers `_km_i64_*_ck` and the signed 128-bit `_km_i128_*`; includes `lib/alloc.cyr`; M4f-2:
  the yields and percent composition rounded once (`_km_cr_mul_div`, since B1 `src/mathx.cyr`'s), the exact limiting
  reagent (`_km_ratio_lt`) and the exact empirical formula, ADR 0002 #24–#26, on `src/mathx.cyr`'s `_km_bn_*` and
  `_km_cr_ratio`; the overflow audit's B7: balancing on unbounded integers where Rust's i64 overflows, ADR 0002 #39, on
  its own signed multi-word integers `_km_sb_*` and the size bound `_km_balance_bits`), `src/numerics.cyr` (M5a, after
  spectroscopy: hisab 1.4.0's `_km_rk4_trajectory`, `_km_dopri45`, `_km_gaussian_elimination` and
  `_km_euler_maruyama_system`, the PCG32 draws `_km_pcg32_next_f64` / `_km_pcg32_next_normal` on hisab 3.4.1's
  `pcg32_next`; includes `lib/alloc.cyr` and `lib/hisab.cyr`), `src/reaction_dynamics.cyr` (M5a: the eight `simulate_*` and
  the built-in rate laws' callbacks, the backward-Euler loop); plus the
  `src/main.cyr` stub. Generated files are written by
  `scripts/oracle/` and never edited. Every transcendental the modules call is correctly rounded; `gas` and `nuclear`
  include `lib/ganita.cyr` (cbrt), `gas` `lib/math.cyr` (`f64_clamp`; `mathx.cyr` no longer needs it since M4d). `dist/kimiya.deps`: string, vec, str, math, ganita, alloc (M4e), tagged; the bundle also
  includes `lib/hisab.cyr` (thermo's and thermochem's `EPSILON_F64`; since M5a numerics' PCG32 and pivot threshold), so a
  consumer needs kimiya's `[deps.hisab]`.

## Tests

- `tests/parity.tcyr` (+ `tests/parity_core.cyr`, `tests/ulp.cyr`; since M5a `tests/rd_systems.cyr` and
  `tests/rd_rows.cyr`, the test systems and the trajectory rows): the golden-corpus harness — 19 modules
  (constants, electrochemistry, element, gas, inorganic, kinetics, math, molecule, nuclear, organic, phase, potential,
  reaction, reaction_dynamics, solution, spectroscopy, stoichiometry, thermo, thermochem; `numerics` is
  `tests/numerics.tcyr`'s), 339,944 ok, 0 failed, 0 measured, 0 pending, 46,407 accepted (810
  under ADR 6, against the corrected oracle, M5a's `simulate_adaptive` row among them; under ADR 2 the 440 ADR 0002 #4 mapped Antoine rows, the 13 #12 exact-CFSE and
  23 #3 SEMF rows printed by name, and, counted, the 1,810 #13, 208 #14, 59 #15, 1,041 #16, 1,654 #19,
  3,605 #20, 6,617 #21 (per function: 3,823 levels, 1,351 transition energies, 1,376 wavelengths, 67 series lines),
  2,837 #22 rows whose note names the port's value, 354 #23 rows (7 `isotope_abundances`, 347 `isotope_pattern`) and 9,090
  #7 `assign_oxidation_states` rows whose note names its whole result, the 78 #7 `balance_equation` rows with no
  integer-matrix counterpart, and the 3,154 stoichiometry rows of #24 (1,126 `percent_yield`, 240
  `theoretical_yield_moles`, 598 `percent_composition`), #25 (920 `limiting_reagent`) and #26 (270 `empirical_formula`)
  whose note names the port's value, index or whole result, and the overflow audit's B1: 682 #27 (151 `heat_transfer`,
  171 `specific_heat_from_calorimetry`, 82 + 82 colligative, 196 `osmotic_pressure`) and 636 #28 rows (89
  `second_order_concentration`, 90 `zero_order_half_life`, 118 `steady_state_intermediate`, 179 `pre_equilibrium_rate`,
  160 `michaelis_menten`), and B2's 2,773 #29 (94 `dilution`, 148 `solution_conductivity`, 904 `faraday_mass_deposited`,
  644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 136 `nernst_potential`, 17 of them B4's correction) and 2,106
  #30 rows (195 `absorbance`,
  219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant`,
  118 `diatomic_moment_of_inertia`), and B3's 1,157 #31 (146 `ideal_gas_pressure`, 148 `ideal_gas_volume`, 222
  `gas_density`, 206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234 `redlich_kwong_pressure`; 98 of them the
  port's Err) and 4,182 #32 rows (582 `lennard_jones`, 244 with ε ≥ 2¹⁰²², M4d's 180 #17 rows among them; 1,860
  `morse_width_parameter`, 340 with D_e ≥ 2¹⁰²³, M4d's 181 #18 rows among them; 730 `coulomb`, 679 `coulomb_force`, 175
  `born_lande_lattice_energy`, 156 `binding_energy_per_nucleon`), and B4's 697 #33 (28 `hess_law`, 27
  `enthalpy_from_formation`, 115 `reaction_enthalpy`, 113 `reaction_gibbs_energy`, 66 `reaction_entropy`, 58 `q_value`, 13
  `ionic_strength`, 74 `born_haber_lattice_energy`, 57 `gibbs_free_energy`, 54 `kohlrausch_molar_conductivity`, 31
  `limiting_molar_conductivity`, 61 `zero_order_concentration`) and 210 #34 rows (185 `Shomate::cp`, 25
  `enthalpy_antiderivative`), and B5's 845 #35 rows (326 `equilibrium_constant`, 27 of them the port's
  ComputationError where Rust's K is finite and 33 its K where Rust's overflows; 55 `debye_huckel_limiting`, 385
  `debye_huckel_extended`, 79 `davies_activity`), and B6's 418 #36 `fft_spectrum`, 222 #37 `find_peaks` and 273 #38
  `adiabatic_flame_temperature` rows, and B7's 413 #39 `balance_equation` rows), 0 glibc-matching, on x86_64 ELF,
  x86_64 PE and aarch64
  (`_PAR_MAX` 8448
  tokens since M4e, for the 2¹³-sample FFT rows); ADR 0002 #6's renamed items compared
  under the Rust names (class A, summarised). Every class-C function and every `math` primitive at bound 0 (ADR 0006;
  sin, cos, atan2 and acos on M4a's primitives since M4d)
- `tests/golden/`: `constants.tsv` (45, all compared since M4b), `electrochemistry.tsv` (19,420), `element.tsv` (975),
  `gas.tsv` (22,773), `inorganic.tsv` (15,586), `kinetics.tsv` (31,381), `math.tsv` (26,147), `molecule.tsv` (2,472),
  `nuclear.tsv` (54,312), `organic.tsv` (9,104), `phase.tsv` (4,174), `potential.tsv` (33,464), `reaction.tsv` (2,406),
  `solution.tsv` (30,644), `spectroscopy.tsv` (70,777), `stoichiometry.tsv` (48,297), `thermo.tsv` (1,828),
  `thermochem.tsv` (11,592), since M5a `numerics.tsv` (243) and `reaction_dynamics.tsv` (954), `modules.txt` — generated by
  `scripts/oracle/`; the corrected oracle (ADR 0006) `gas.cr.tsv` (11), `kinetics.cr.tsv` (323), `math.cr.tsv` (96),
  `nuclear.cr.tsv` (23), `phase.cr.tsv` (4), `potential.cr.tsv` (3), `reaction.cr.tsv` (2), `solution.cr.tsv` (229),
  `spectroscopy.cr.tsv` (118), since M5a `numerics.cr.tsv` (3: Box–Muller draws, four cos and one ln misrounded by glibc)
  and `reaction_dynamics.cr.tsv` (1: two pows of a 1,766-point Dormand–Prince run) —
  written by `validate.py --write-corrections`
- `tests/numerics.tcyr` (M5a): `numerics.tsv` through `tests/parity_core.cyr` (240 ok, 3 accepted under ADR 6: PCG32's
  u32 draws class A, the normal draws, dopri45 and Euler–Maruyama class C at bound 0, RK4 and the elimination class B,
  every row 0 ulp), then 5 Cyrius-only groups (the heap work tables, the trajectory capacity, the scratch reuse and
  Neumaier's branch, the elimination's pivot edge, the PCG32 handle; 21 assertions)
- `tests/reaction_dynamics.tcyr` (M5a): the 26 Rust tests of `reaction_dynamics.rs` (`stiff_resized_jacobian_is_error_not_panic`
  adapted: the flat Jacobian can't be resized, so its group checks a Jacobian that writes nothing) and the 2 integration
  tests, plus 5 Cyrius-only groups (the trajectory buffer, out-of-type inputs, the inline constants' bits and MAX_STEPS'
  edges, the check order, the reused Jacobian); 74 assertions
- `tests/mathx.tcyr`: the correctly rounded primitives (M3c-2, M4a, M4d, M4e) — special values, exact and tie cases (10²³,
  the subnormal ties), thresholds, subnormal results, the known hardest-to-round arguments, Payne–Hanek arguments, sincos,
  atan2's sign and zero matrix, acos, x·2^y, 1/(100·x), the fast path's fallbacks, and seeded fast / accurate agreement,
  every expected value MPFR's; since M4e-2 the unbounded-exponent arithmetic (ADR 0002 #20: 48 MPFR values past binary64's
  range, 20,000 seeded pairs against IEEE's own operations); since M4e-3 the multi-word integers, `_km_cr_ratio` (ties,
  subnormal results, the overflow threshold) and the double-double quotient (10,000 seeded quotients, exact against
  double-double); since the overflow audit's B1 the correctly rounded monomial (36 exact cases, the sign and tie
  helpers, 20,000 seeded arguments against the exact path), since B2 its four new shapes and the correctly rounded √ of
  a quotient (hand cases with MPFR's values, 20,000 seeded arguments against the exact paths), since B3 the correctly
  rounded w·√((x/y)·2^e) and `_km_ux_normal` (ties at both grids, the Morse width's edges, the domain; 20,000 seeded
  arguments against the exact path), since B10a the unbounded e^x, e^x − 1 and 10^x (96 MPFR rows at its widest exponent
  range as bit patterns, each by the function, the fast path and the accurate path; the specials; 4,000 seeded arguments
  each against the accurate path and the binary64 function) (36 groups, 1,960 assertions); `scripts/oracle/cr_verify.py` checks ≥ 1,000,000 arguments per primitive
  (611,518 pairs for x·2^y, 1,049,032 for 1/(100·x); the five unbounded-exponent final operations 1,040,039, the three
  geometry functions 114,000 edge geometries, the slow path forced; M4e-3's `ratio` 320,000, `bohr` 260,429, `trans` and
  `ryd` 555,944 each, `series` 55,345 and `rmass` 502,000 against `Fraction`, the exact path forced; M4f-2's `muldiv`
  640,000, `pyield` 397,518, `tyield` 340,006, `pcomp` 223,604, `limit` 400,008 and `empir` 223,899 against `Fraction`,
  the exact paths forced; B1's `mono` 430,010, `monotie` 240,000 and the ten functions' `heat`, `boil`, `frz` 281,000
  each, `shc` 200,010, `osm` 200,008, `zhl` 185,006, `ssi` 220,007, `pre` 220,008, `soc` 200,010, `mm` 200,008 against
  `Fraction` and `mpfr_ref.UX`, the slow paths forced; B2's `sqrtdiv` 326,845 and the twelve functions' `dil` 220,008,
  `cond` 250,121, `fmass` 210,010, `fcharge` 210,007, `fratio` 210,008, `nernst` 200,009, `abs` 281,000, `conc` 200,008,
  `vib` 266,885, `rotor` 200,009, `rot` 200,006, `moment` 311,000 against `Fraction`, exact integer roots and
  `mpfr_ref.UX`, the slow paths forced; B3's `wsqrt` 316,971 and the twelve functions' `igp`, `igv` 220,010 each, `dens`
  200,008, `vdw` 180,009, `prp` 180,007, `rkp` 180,006, `lj` 200,009, `width` 316,695, `coul`, `force` 210,700 each,
  `lande` 200,008, `bea` 180,009 against `Fraction`, exact integer roots and `mpfr_ref.UX`, the slow paths forced; B4's
  twelve families `hess` 110,008, `eff`, `rsum`, `qv`, `haber` 110,005 each, `ionic`, `gibbs`, `kohl` 100,006 each, `lim`
  100,004, `zoc` 100,005, `shcp`, `shh` 100,007 each against `mpfr_ref.UX` with `ux_sum`, the slow paths forced, and
  `nernst` 202,809 with #29's correction; B5's `keq` 144,934, `dhl` 100,012, `dhe` 130,021, `dav` 130,015 against MPFR's
  exp / exp10 of the 53-bit argument on `mpfr_ref.UX` (`ux_exp_final`, `ux_exp10_final`), the slow paths forced; B6's
  `flame` 14,010 against Rust's solver with the exact sign tests (and its bisection's test), `pk3` … `pk16` 95,000 each
  against `Fraction` (every candidate through the exact predicate forced), `fft1` … `fft256` 110,718 signals against the
  transform on `mpfr_ref.UX` (the unbounded transform forced), `signsc` 140,110 sums of `_km_sign_scaled` against
  `Fraction`; B7's `bal23` … `bal2425` 54,306 matrices against Rust's algorithm on Python ints, the fast path's overflow
  report against a tracer of Rust's checked operations, the multi-word path forced; B10a's `exp_ux` 139,453, `expm1_ux`
  168,035, `exp10_ux` 140,700 against MPFR at its widest exponent range (`mpfr_ref.UX.exp`, `.expm1`, `.exp10`), the
  accurate path forced, the binary64 function's pair checked wherever it is normal) on every target
- Reference tests, the Rust tests ported one for one: `tests/error.tcyr` (2 Rust + 1 Cyrius-only group,
  24 assertions), `tests/validate.tcyr` (3 + 1, 25), `tests/element.tcyr` (22 + 1, 96), `tests/molecule.tcyr`
  (11 + 1, 39; `serde_roundtrip` doesn't port), `tests/reaction.tcyr` (10 + 3, 49; B4's #33, B5's #35), `tests/thermo.tcyr`
  (10 + 3, 44; B1's #27, B4's #33), `tests/thermochem.tcyr` (62 + 9, 421; 4 tests assert J/mol for ADR 0002 #1; B4's #33
  and #34, B6's #38), `tests/phase.tcyr`
  (31 + 3, 126; B1's #27), `tests/gas.tcyr` (39 + 3, 181; B3's #31), `tests/solution.tcyr` (49 + 8, 165; B2's #29, B4's #33, B5's #35),
  `tests/electrochemistry.tcyr` (35 + 5, 1,413; B2's #29 and B4's correction), `tests/kinetics.tcyr` (47 + 6, 159; 2
  tests adapted for the
  Lineweaver–Burk buffer; B1's #28, and one Cyrius-only assertion of Rust's 1/∞ = +0 now the port's [A]₀; B4's #33),
  `tests/organic.tcyr` (21 + the 2 organic integration tests + 4, 136; 5 tests use ADR 0002 #6's names; the acceptor
  property's reference tests), `tests/inorganic.tcyr` (30 + 9, 342; the exact CFSE of ADR 0002 #12, the exact Δ_tet of #15, B3's #32, B4's #33),
  `tests/nuclear.tcyr` (36 + 10, 252; `semi_empirical_invalid_is_nan` asserts ADR 0002 #3's Err; #13's correctly
  rounded and #14's exact values, B3's #32, B4's #33), `tests/potential.tcyr` (33 + 9, 192; `morse_width_parameter_matches_force_constant`
  calls `vibrational_frequency` since M4e; null pointers, the helpers, the tolerances to the ulp and the −π → π mapping,
  #16's r_min, #17 / #18 (folded into #32 since B3); M4e-2's #20 edges against MPFR, the scale invariance bit for bit,
  the forced slow path against the fast one; B3's #32),
  `tests/spectroscopy.tcyr` (62 + the 2 spectroscopy integration tests + 13, 403; out-of-type integers, the caller
  buffers,
  the isotope helpers' stable sort and thresholds to the ulp, the FFT transcription, #19, M4e-3's #21 and #22 (exact
  values, ties, the fast paths against the exact ones), M4e-4's #23 (the table cell by cell, S₈ bit for bit against the
  Python replica; `hydrogen_abundances_iupac_2009` asserts the table's pairs and sums sulfur too), R_H, B2's #30, B6's
  #36 (the survey's examples, the fast test's bounds, n = 4's product below 2⁻¹⁰²², 3,000 seeded signals against the
  forced slow path) and #37 (the ties, ±DBL_MAX thresholds, `_km_sign_scaled`, 20,000 seeded arrays against the exact
  predicate forced)),
  `tests/stoichiometry.tcyr` (63 + 11, 350; 7 tests adapted for ADR 0002 #7 — the non-integer entries' integer
  counterparts, the near-integer entries as integers, KO₂ −1/2, O₂⁻ −1/2, C₃H₈ −8/3, Fe₃O₄ +8/3, (O, u32::MAX)×2 at
  i32::MIN −2³⁰/(2³² − 1); the rational states' i32 range, the 2⁵³ bound and the heap tables, the checked i64 / i128
  helpers at their edges, out-of-type words, the buffers; M4f-2's #24 (ties at RN, the range ends, Rust's ±0, the
  midpoint resolution), #25 (over- and underflowing ratios, ties and near-ties, `_km_ratio_lt`'s branches) and #26 (the
  exact 1/10 edges, the count's 2³² edge, subnormal percentages, the fast and exact tolerance tests); the overflow
  audit's B7's #39 (the survey's example, every reachable overflow site, each Err after an overflow, u32::MAX against
  2³², the buffer checks after a multi-word Ok, 32 × 33 dense, 600 seeded matrices against the forced multi-word path,
  the multi-word helpers at their edges); no Rust test's assertion changes)
- `tests/kimiya.tcyr`: smoke test of `dist/kimiya.cyr` as a consumer includes it (47 assertions)
- `tests/kimiya.bcyr`: `element/lookup_by_symbol`, `element/lookup_by_number`, `molecule/molecular_weight_glucose`,
  `reaction/gibbs_free_energy`, `thermo/heat_transfer`, `gas/ideal_gas_pressure`, `gas/peng_robinson_z_n2_1mpa`,
  `thermochem/reaction_enthalpy_ch4`, `thermochem/vant_hoff_k`, `thermochem/enthalpy_change_cp`,
  `thermochem/shomate_cp`, `thermochem/adiabatic_flame_ch4_air`, `solution/ph_from_h`, `solution/weak_acid_ph`,
  `electrochemistry/lookup_half_reaction`, `electrochemistry/nernst_potential`, `kinetics/arrhenius_rate`,
  `kinetics/eyring_rate`, `kinetics/michaelis_menten`, and (Cyrius-only) `mathx/cr_exp`, `cr_expm1`, `cr_exp10`,
  `cr_ln`, `cr_log10`, `cr_pow` beside `ganita_f64_pow` (two arguments each), `organic/enthalpy_from_bonds`,
  `inorganic/lookup_ionic_radius` (M4b; rust-old/benches has no organic or inorganic bench),
  `nuclear/binding_energy_per_nucleon_fe56`, `nuclear/semi_empirical_fe56` (M4c, against Rust) and the Cyrius-only
  `nuclear/decay_chain_two_step`, the M4a trig benches (`mathx/cr_sin`, `cr_cos`, `cr_sincos`, `cr_atan2`, `cr_acos`
  beside `f64_sin`, `f64_cos`, `ganita_f64_atan2`), `potential/dihedral_angle` (M4d, against Rust) and the Cyrius-only
  `potential/bond_angle`, `lennard_jones`, `lj_equilibrium_distance`, `morse`, `mathx/cr_exp2_scaled`, `distance_3d` and
  the #20 slow path's `distance_3d_2p600`, `bond_angle_2p600`, `dihedral_angle_2p600` (M4e-2),
  `spectroscopy/absorbance`, `bohr_energy_level`, `isotope_pattern_glucose` (M4e, against Rust) and the Cyrius-only
  `spectroscopy/fft_spectrum_256` / `_1024` / `_4096`, `find_peaks_1000`, `wavenumber_from_wavelength`, and
  `bohr_transition_energy`, `rydberg_wavelength`, `series_wavelength_balmer_alpha`, `reduced_mass_hcl` and the exact
  path's `bohr_transition_energy_tie`, `reduced_mass_tie` (M4e-3), and the sulfur-bearing `isotope_pattern_cystine`,
  `isotope_pattern_s8` (M4e-4), `stoichiometry/balance_ch4_combustion`, `oxidation_states_kmno4` (M4f, against Rust) and
  the Cyrius-only `balance_kmno4_hcl`, `oxidation_states_fe3o4`, `limiting_reagent`, `theoretical_yield_moles`,
  `percent_yield`, `percent_composition_glucose`, `empirical_formula_glucose`, and ADR 0002 #24–#26's exact paths
  `theoretical_yield_moles_tie`, `percent_yield_subnormal`, `limiting_reagent_tie`, `empirical_formula_edge` (M4f-2)
- `tests/audit.bcyr` (the overflow audit's Cyrius-only benches, moved from `tests/kimiya.bcyr` in B4 so each file stays
  well inside `cyrius bench`'s 300 s deadline; names unchanged): B1's `thermo/specific_heat_from_calorimetry`, `phase/boiling_point_elevation`,
  `freezing_point_depression`, `osmotic_pressure`, `kinetics/second_order_concentration`, `zero_order_half_life`,
  `steady_state_intermediate`, `pre_equilibrium_rate` and the edges `thermo/heat_transfer_tie`,
  `phase/osmotic_pressure_subnormal`, `kinetics/second_order_concentration_edge`, `michaelis_menten_edge`, and B2's:
  `solution/dilution`, `solution_conductivity`, `electrochemistry/faraday_mass_deposited`, `faraday_charge_required`,
  `faraday_mass_ratio`, `spectroscopy/concentration_from_absorbance`, `vibrational_frequency`, `rigid_rotor_energy`,
  `rotational_constant`, `diatomic_moment_of_inertia` and the edges `electrochemistry/nernst_potential_edge`,
  `spectroscopy/vibrational_frequency_edge`, `vibrational_frequency_subnormal`,
  `electrochemistry/faraday_mass_ratio_tie`, and B3's: `gas/ideal_gas_volume`, `gas_density`, `van_der_waals_pressure`,
  `peng_robinson_pressure`, `redlich_kwong_pressure`, `potential/coulomb`, `coulomb_force`, `morse_width_parameter`,
  `inorganic/born_lande_lattice_energy` and the edges `gas/van_der_waals_pressure_edge`, `peng_robinson_pressure_edge`,
  `redlich_kwong_pressure_edge`, `potential/lennard_jones_edge`, `morse_width_parameter_subnormal`,
  `morse_width_parameter_tie`, `inorganic/born_lande_lattice_energy_edge`, `nuclear/binding_energy_per_nucleon_edge`, and
  B4's: `thermo/hess_law`, `enthalpy_from_formation`, `thermochem/reaction_gibbs_energy`, `reaction_entropy`,
  `nuclear/q_value`, `solution/ionic_strength`, `inorganic/born_haber_lattice_energy`,
  `solution/kohlrausch_molar_conductivity`, `limiting_molar_conductivity`, `kinetics/zero_order_concentration`,
  `thermochem/shomate_enthalpy_antiderivative` and the edges `thermo/hess_law_edge`, `enthalpy_from_formation_edge`,
  `thermochem/reaction_enthalpy_edge`, `nuclear/q_value_edge`, `solution/ionic_strength_edge`,
  `inorganic/born_haber_lattice_energy_edge`, `reaction/gibbs_free_energy_edge`,
  `solution/kohlrausch_molar_conductivity_edge`, `limiting_molar_conductivity_edge`,
  `kinetics/zero_order_concentration_edge`, `thermochem/shomate_cp_edge`, `shomate_enthalpy_antiderivative_edge`
- `tests/audit2.bcyr` (the overflow audit's Cyrius-only benches from B5 on: with them `tests/audit.bcyr` would take 219
  s of the 300 s deadline, so the later batches open this third file, as B4 split; 66 s with B10a's): B5's
  `reaction/equilibrium_constant`, `solution/debye_huckel_limiting`, `debye_huckel_extended`, `davies_activity` and the
  edges `reaction/equilibrium_constant_edge`, `equilibrium_constant_subnormal`, `solution/debye_huckel_limiting_edge`,
  `debye_huckel_extended_edge`, `debye_huckel_extended_overflow`, `davies_activity_edge`, and B6's edges
  `spectroscopy/fft_spectrum_256_huge`, `fft_spectrum_256_subnormal`, `fft_spectrum_256_squares`,
  `find_peaks_1000_subnormal`, `thermochem/adiabatic_flame_tiny`, and B7's `stoichiometry/balance_survey_2x3`,
  `balance_overflow_10x9`, `balance_dense_30x31`, and B10a's `mathx/cr_exp`, `cr_exp_ux`, `cr_exp_ux_acc`,
  `cr_exp_ux_subnormal`, `cr_exp_ux_minus_1e4`, `cr_exp_ux_1e12`, `cr_exp_ux_1e15`, `cr_exp_ux_acc_minus_1e4`,
  `cr_expm1`, `cr_expm1_ux`, `cr_expm1_ux_1e4`, `cr_exp10`, `cr_exp10_ux`, `cr_exp10_ux_1e4`,
  `cr_exp10_ux_acc_minus_1e12`
- `tests/numerics.bcyr` (M5a, a new file: 39 s): `dynamics/first_order_100steps` and `dynamics/consecutive_100steps`
  (Rust's benches) and the Cyrius-only `dynamics/adaptive_first_order`, `dynamics/stiff_robertson_400steps`,
  `dynamics/stochastic_100steps`, `numerics/gaussian_elimination_8x8`, `numerics/pcg32_next_normal`
- `tests/kimiya.fcyr`: scaffold fuzz harness

## Dependencies

- stdlib — syscalls, string, alloc, str, fmt, io, vec, args, math, ganita, overflow, tagged, result, fnptr
- hisab 3.4.1 (git tag) → sakshi 2.5.10
- test stdlib — assert, bench
- tooling only (not the shipped library): the system MPFR 4.2 (`libmpfr.so.6`, through ctypes) for
  `scripts/oracle/mpfr_ref.py`, `validate.py`, `mathx_tables.py` and `cr_verify.py`; wine and qemu-aarch64 for the PE and
  aarch64 test runs

## Consumers

- **Rust (crates.io, as of 2026-11-08), all on the 1.x crate:** jivanu, khanij (Arrhenius, pH, first-order
  kinetics, Born–Landé, ionic radii, thermochem fields, `Molecule::new`), jnana (`ELEMENTS`, three constants),
  kiran (re-exports whole modules), tara (declared, unused).
- **Cyrius:** none yet. tanmatra 2.0.0 deferred its kimiya bridges until kimiya ports.

## Next

**M5b** numerics (port plan §4): network (`matrix_rank` from hisab's SVD or transcribed), fitting (hisab 1.4.0's
`least_squares_poly` transcribed into `src/numerics.cyr`), sensitivity (ADR 0002 #2's closed forms, `compute_gradient` on
hisab's tape). M5a (numerics' ODE integrators and dense solve, reaction_dynamics) is done. The library-wide overflow
audit (below) is **paused after B10a** by the user's decision of 2026-11-12: it resumes after M5 / M7 over the whole
library, M5's modules included (the roadmap lists M5a's candidate spots), so M5 is ported at the plan's parity classes.

Follow-ups:
- **The library-wide overflow audit — in progress (B1–B7 done: the survey's B1–B9; B2 took its B2 and B3, B3 its B4 and
  B5, B4 its B6, B5 its B7, B6 its B8, B7 its B9; and B10a, the first part of its B10).** The survey
  (sixteen families over every ported module, with replicas, synthetic grids and probes; the main session's record)
  batches the fixes B1 … B12; the user's decisions: the audit itself (2026-11-09), "Round once" for a formula that is a
  single product or quotient of inputs, #20's edge fix for one with a sum inside (2026-11-10), and the correctly rounded
  ±∞ or 0 where Rust gives NaN for an out-of-range value. **B1 ✅ (2026-11-10)**: F1 → ADR 0002 #27 (thermo's
  calorimetry, phase's colligative products) and F4 → #28 (the kinetics rationals), on the new correctly rounded
  monomial `_km_cr_mono` and the `_km_ux_*` helpers. **B2 ✅ (2026-11-10)**: F3 → ADR 0002 #29 (solution's `dilution` and
  `solution_conductivity`, Faraday's three laws rounded once; `nernst_potential` with an unbounded exponent behind a
  range test) and F5 → #30 (spectroscopy's Beer–Lambert, rotor and moment-of-inertia monomials rounded once;
  `vibrational_frequency` the correctly rounded √(k/μ), the new `_km_cr_sqrt_div`: a proven double-double fast path and
  an exact integer root). **B3 ✅ (2026-11-11)**: F2 → ADR 0002 #31 (gas: the ideal gas law and density rounded once; van
  der Waals, Peng–Robinson and Redlich–Kwong pressure with an unbounded exponent, the co-volume tests on the unbounded
  co-volume, Tc^2.5 exactly rescaled) and F6 → #32 (potential's Coulomb pair rounded once, the correctly rounded Morse
  width on the new `_km_cr_mul_sqrt_ratio`, Lennard-Jones, Born–Landé and the binding energy per nucleon with an
  unbounded exponent; M4d's #17 and #18 folded in). **B4 ✅ (2026-11-11)**: F7 → ADR 0002 #33 (the sums with
  cancellation, in Rust's summation order, behind the pure-sum test "result finite") and F8 → #34 (the Shomate
  polynomials behind a box test), #29's fast path corrected (2⁻¹⁰²¹, not 2⁻¹⁰²²), the benches split into
  `tests/audit.bcyr`. **B5 ✅ (2026-11-11)**: F9 → ADR 0002 #35 (`equilibrium_constant`'s and the three activity
  models' exp / 10^x arguments with an unbounded exponent, the argument's final operation rounded once, then the
  unchanged correctly rounded exp / 10^x; no new primitive; its benches open a third file, `tests/audit2.bcyr`).
  **B6 ✅ (2026-11-11)**: F13 → ADR 0002 #36 (the FFT with an unbounded exponent over the whole transform, behind a
  per-operation operand test), F14 → #37 (`find_peaks`' exact threshold, `_km_sign_scaled`), F16's flame → #38 (the
  bracket's exact sign test; B12 extends the row). **B7 ✅ (2026-11-12)**: F15 → ADR 0002 #39 (`balance_equation` on
  unbounded integers: Rust's i64 code first, the same steps on signed multi-word integers from its first overflow,
  Lehmer's gcd, the tables from a proven Hadamard bound). **B10a ✅ (2026-11-12)**: F10's primitives, for ADR 0002 #40
  (`_km_cr_exp_ux`, `_km_cr_expm1_ux`, `_km_cr_exp10_ux`: 53-bit results at any exponent up to MPFR's widest range; no
  row changes). Next, in the survey's numbering: the rest of B10 (#40's functions: kinetics; nuclear; van 't Hoff,
  Clausius, Antoine, Morse and `fugacity_pr`), B11 (F11, F12; `_km_cr_pow_ux`, `_km_cr_log10_ux`), B12 (the weak acid /
  base Newton and Peng–Robinson's coefficients).
  The history:
- **The library-wide audit of intermediate over- and underflow** (the user's decision, 2026-11-09, with ADR 0002 #20):
  after M4, check every ported module for products and quotients of inputs that over- or underflow in Rust's chain where
  the result is representable, and fix each edge as #20 does — Rust's sequence evaluated with an unbounded exponent, the
  final operation rounded once, rows whose chain stays normal keeping Rust's bits, a cheap fast-path range test. The
  machinery exists: `src/mathx.cyr`'s `_km_ux_*` (pairs, + − × ÷ √, cmp / max, `_km_ux_round` and the final
  `_km_ux_{add,sub,mul,div,sqrt,atan2}_f64`), `mpfr_ref.UX` / `unbounded` / `ux_*_final` for `validate.py`'s proofs,
  `exact::Ux` for the generator's notes, `cr_verify.py`'s `ux*` families. First candidates, from M4d's and M4e's kept
  lists: `lennard_jones`' (σ/r)⁶ and (σ/r)¹², Coulomb's k·q₁·q₂ and r·r, `morse`'s subnormal D_e·x, Beer–Lambert's ε·l·c
  and ε·l, k/μ, h·B, 8π²·I, μ·r·r, the FFT's butterfly sums and re² + im² (`reduced_mass`' m₁·m₂ and m₁ + m₂ left the list
  in M4e-3: ADR 0002 #22 covers them), and since M4f `balance_equation`'s checked i64 elimination where a valid u32
  solution exists (an exact rational or multi-word elimination; done in the audit's B7, ADR 0002 #39). Each fix is an
  ADR 0002 row, accepted before it lands. M4e-3's `_km_cr_ratio`
  (RN(p/q·2^e) of exact integers) is the tool where a quotient of products is the whole formula.
- **M4f, measured and reported to the user; fixed in M4f-2** by the user's decisions of 2026-11-10 (ADR 0002's M4f
  addition has the measurements, its M4f-2 addition the fixes): `percent_yield`, `theoretical_yield_moles` and
  `percent_composition` rounded once (#24), `limiting_reagent`'s exact comparison (#25), `empirical_formula`'s exact
  ratios and 1/10 tolerance (#26). Still open: `balance_equation`'s checked i64 elimination failing where a valid u32
  solution exists (104 corpus rows; [[3037000500, 1, −3037000501], [1, 3037000500, −3037000501]] with x = [1, 1, 1]) —
  for the overflow audit above.
- **hisab's FFT twiddles are not correctly rounded** (M4e, ADR 0003 amendment): hisab 3.4.1's `num_fft` takes them from
  the stdlib's `f64_cos` / `f64_sin` and reseeds them every 64 steps, so its transforms differ from 1.4.0's and from a
  correctly rounded reference in the last bits from n = 256 on, and by platform. kimiya transcribes 1.4.0's `num::fft`
  with `_km_cr_sincos` twiddles; a hisab FFT on correctly rounded twiddles (or with 1.4.0's recurrence as an option) would
  let it call hisab. For the hisab feature request.
- **M4e, reported to the user (kept as Rust's under the policy):** `reduced_mass(m, m)` one ulp from m/2 on 8.6 % of masses
  (33 of 315 golden rows); `bohr_energy_level(Z, Z)` one ulp from −R on 3 of 12 (and Z/n = 1/k from −R/k² on 30 %); the
  Rydberg factor 1/n₁² − 1/n₂² cancelling for close large n (an exact factor would move `rydberg_wavelength` up to 3.9·10⁷
  ulp at n ≈ 2³¹, 39 ulp at n = 100, 939 at n = 10⁴) and `bohr_transition_energy`'s E₂ − E₁ likewise; the over- and underflow of
  products and quotients of inputs (Beer–Lambert's ε·l·c and ε·l, `reduced_mass`'s m₁·m₂ and m₁ + m₂, k/μ, h·B, 8π²·I,
  μ·r·r, the FFT's sums and re² + im², which hypot would avoid); `find_peaks`' interpolated threshold; the isotope
  thresholds 0.001 (only an intensity of exactly RN(0.001) is judged differently from 1/1000) and 10⁻¹². ADR 0002's M4e
  addition lists them. Fixed: #19; in M4e-3, by the user's decisions of 2026-11-10, #21 (the Bohr levels, the transition
  energy, the Rydberg formula and the series rounded once) and #22 (the reduced mass rounded once, its over- and underflow
  included).
- **M4e data, changed in M4e-4** (ADR 0002 #23, the user's decision of 2026-11-10): Rust's isotope masses are AME2020's to
  five decimals and its sulfur omits ³⁶S (summing to 0.9999); the port's table is the full AME2020 masses and the TICE
  2009 abundances with ³⁶S, from `scripts/oracle/isotopes.tsv` (`scripts/oracle/isotope_sources.py` re-verifies it
  against fresh downloads). CIAAW's current table gives these elements intervals only, so TICE 2009's single values stay
  the representative compositions; a later CIAAW single-value table would be another data change.
- **CI: a truncated parity run.** Concurrent `cyrius test --win` runs under wine can print a truncated parity run and
  still exit 0, so the CI gate must check that the parity summary line (`=== parity summary: … ===`) is present (and
  that `cyrius test` reports every test file). The harness cannot detect this itself: whether wine cuts the program's
  output or ends it early and still reports 0, the exit status and the output are what wine hands back, so nothing the
  program prints or returns can mark a cut-off run. Running the three targets one after another (never two wine runs at
  once) avoids it; M4d's gates did.
- **M4d, reported to the user (kept as Rust's under the policy; no single exact value to restore):** in `potential`, the
  over- and underflow of (σ/r)⁶ and (σ/r)¹² (`lennard_jones` for σ/r past 2^85.3 or below 2^-89.5, though 4ε·d may be
  representable), of k·q₁·q₂ and of r·r (`coulomb`, `coulomb_force` for r past 1.3 × 10¹⁵⁴ or below 1.5 × 10⁻¹⁵⁴), of the
  squared coordinates in `distance_3d`, `bond_angle` and `dihedral_angle` (beyond about 10¹⁵⁴, or 10⁷⁷ for the dihedral's
  fourth powers, and below the reciprocals: an Err, 0, ±∞ or a degraded value where the distance or angle is
  representable, though Rust documents the angles as scale-invariant: fixed in M4e-2 as ADR 0002 #20, the user's
  decision), `morse`'s r − r_e overflow and its subnormal D_e·x. Fixed in M4d: #16, #17, #18. The two numerical tolerances (10⁻¹², 10⁻²⁰) are kept, as
  `EQUAL_RATE_REL_TOL`.
- **ADR 0002 #13 is correctly rounded since M4d** (`_km_cr_exp2_scaled`; M4c's faithful fractional-n path is gone).
  `decay_constant`'s LN_2/t½ (1 ulp from CR(ln 2/t½) on 476 of 2,175 rows) is kept under the policy.
- **Two-argument primitives' hard cases.** `_km_cr_exp2_scaled`, like pow and atan2, is proven correct for results not
  within 2⁻¹⁵⁰ (relative) of a midpoint; no exhaustive worst-case search exists for two arguments (a random pair is that
  close with probability about 2⁻⁹⁶). `_km_cr_sixth_root2_mul` needs no such caveat (one argument: its continued-fraction
  bound is proven).
- **M4b, open questions for the user** (port plan §6.4): under ADR 0002 #12's policy, `delta_tet_from_oct` differs from
  CR(4Δ/9) only where Δ·4 overflows (|Δ| > DBL_MAX/4, 5 golden rows: decided in M4c, ADR 0002 #15), and Born–Landé's
  (1 − 1/n) from CR((n − 1)/n) on 332 rows inside a product with no exact value (left as Rust's). `organic.rs`'s bond energies carry no source
  citation in the oracle (roadmap, *Data to re-check*). The IUPAC 2011 acceptor wording in `src/organic_data.cyr` is
  quoted from a secondary copy; the published PDF was unreachable.
- **pow's accurate path has no third phase.** It is proven correct only for results not within 2⁻¹³⁰ (relative) of a
  midpoint; CORE-MATH's pow has a 256-bit third phase that reports rather than misround (src/mathx.cyr's pow header;
  every published hard case and every verified pair are well clear). A random pair is that close with probability
  about 2⁻⁷⁶.
- **log10 got about 10 ns slower** (to about 77 ns) when pow's absolute ln bound was shared with it in M3c-2b; a separate,
  tighter bound for log10 could win it back.
- **The Cyrius issue** `cyrius/docs/development/issues/2026-11-09-exp-polyfill-double-rounds-subnormal-results.md`
  (filed): the aarch64 exp polyfill rounds subnormal results twice; until a Cyrius library ships correctly rounded
  primitives, kimiya owns them (ADR 0006).
- A hisab feature request (ADR 0003, now also noting that `num_newton` / `num_bisection` take no callback context).
- A note to tanmatra about its `parity_core.cyr` per-function table sizing (768 bytes for a 192-slot limit; overflows past
  96 functions).
- A cyrius note that `distlib`'s sidecar verify misses stdlib leaves a bundle only calls from unreached functions (kimiya
  works around it with explicit `include`s, port plan §6.4).
- From M2a: a note to ganita that its acos differs from glibc by up to 2 ulp and gives different bits on aarch64 than on
  x86_64 (port-parity.md §5); kimiya no longer calls it.
