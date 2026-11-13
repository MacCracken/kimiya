# Benchmarks — Rust 1.2.0 vs the Cyrius port

Each row times one benchmark body: the same work in Rust (`rust-old/`, release build) and in
Cyrius (`tests/kimiya.bcyr`; since the overflow audit's B4 `tests/audit.bcyr` and since B5 `tests/audit2.bcyr` for the
audit's Cyrius-only benches; since M5a `tests/numerics.bcyr`; `cyrius bench` runs all four). The release's full table is
[2.0.0: every benchmark](#200-every-benchmark-2026-11-13); the milestone sections after it appended their rows as they
landed. Numbers are host-specific, so compare within a column set, never across hosts.

## Method

- **Rust, criterion:** the *latest* column of [`rust-old/benchmarks.md`](../rust-old/benchmarks.md)
  (criterion 0.8 medians, v1.2 post-audit, rustc 1.99.1 release, 2026-11-08), from
  `cargo bench --manifest-path rust-old/Cargo.toml --bench benchmarks -- <name>`.
- **Cyrius:** `cyrius bench`, mean of 20 batches. Each batch is a tight loop of `per_batch` calls
  (as criterion's `b.iter` loop), inputs read from globals (Rust's `black_box`) and results written
  to a sink. lib/bench.cyr measures its clock-read floor and subtracts it from every sample
  (1.31 µs per read on this host); the per-op figure is lib/bench.cyr's ps mean.
- **Host:** x86_64 Linux, HPET clocksource. Toolchains: rustc 1.99.1 (release); Cyrius 6.7.18
  (ganita 1.2.15).

## 2.0.0: every benchmark (2026-11-13)

The release's numbers: one `cyrius bench` run of all four bench files on 2026-11-13 (the 2.0.0 tree; nothing
re-optimised for it), on the host and toolchains of the Method above (the clock-read floor measured 1.21 µs).
Rust's column is `rust-old/benchmarks.md`'s *latest* criterion median. The milestone sections below describe
each body, the optimizations and why the gap is what it is; their Cyrius numbers are that milestone's run, so
they differ from these by run-to-run noise (a few per cent) and by the later changes they record.

### The 31 Rust benchmarks

Every bench of `rust-old/benches/benchmarks.rs` has its Cyrius counterpart, the same name and body (`dynamics/*`
and `fitting/*` in `tests/numerics.bcyr`, the rest in `tests/kimiya.bcyr`).

| Bench | Rust 1.2.0 | Cyrius 2.0.0 | Cyrius / Rust |
|---|---|---|---|
| `element/lookup_by_number` | 0.84 ns | 4.63 ns | 5.5× |
| `element/lookup_by_symbol` | 44.73 ns | 67.7 ns | 1.5× |
| `molecule/molecular_weight_glucose` | 3.63 ns | 27.8 ns | 7.6× |
| `stoichiometry/balance_ch4_combustion` | 233.90 ns | 783 ns | 3.3× |
| `stoichiometry/oxidation_states_kmno4` | 22.38 ns | 121 ns | 5.4× |
| `gas/ideal_gas_pressure` | 2.18 ns | 115 ns | 52.8× |
| `gas/peng_robinson_z_n2_1mpa` | 61.39 ns | 282 ns | 4.6× |
| `reaction/gibbs_free_energy` | 1.75 ns | 17.9 ns | 10.2× |
| `thermo/heat_transfer` | 0.95 ns | 65.8 ns | 69.2× |
| `thermochem/adiabatic_flame_ch4_air` | 301.48 ns | 3.34 µs | 11.1× |
| `thermochem/enthalpy_change_cp` | 24.00 ns | 166 ns | 6.9× |
| `thermochem/reaction_enthalpy_ch4` | 38.58 ns | 302 ns | 7.8× |
| `thermochem/shomate_cp` | 2.34 ns | 32.9 ns | 14.1× |
| `thermochem/vant_hoff_k` | 9.51 ns | 74.5 ns | 7.8× |
| `kinetics/arrhenius_rate` | 6.78 ns | 65.6 ns | 9.7× |
| `kinetics/eyring_rate` | 6.73 ns | 60.4 ns | 9.0× |
| `kinetics/michaelis_menten` | 2.92 ns | 24.3 ns | 8.3× |
| `dynamics/consecutive_100steps` | 5.04 µs | 18.7 µs | 3.7× |
| `dynamics/first_order_100steps` | 4.41 µs | 9.70 µs | 2.2× |
| `solution/ph_from_h` | 9.73 ns | 86.5 ns | 8.9× |
| `solution/weak_acid_ph` | 34.85 ns | 221 ns | 6.3× |
| `electrochemistry/lookup_half_reaction` | 57.06 ns | 123 ns | 2.2× |
| `electrochemistry/nernst_potential` | 6.42 ns | 89.8 ns | 14.0× |
| `spectroscopy/absorbance` | 0.98 ns | 65.7 ns | 67.0× |
| `spectroscopy/bohr_energy_level` | 1.21 ns | 46.8 ns | 38.7× |
| `spectroscopy/isotope_pattern_glucose` | 2.80 µs | 15.5 µs | 5.5× |
| `nuclear/binding_energy_per_nucleon_fe56` | 2.02 ns | 27.7 ns | 13.7× |
| `nuclear/semi_empirical_fe56` | 18.50 ns | 131 ns | 7.1× |
| `potential/dihedral_angle` | 20.35 ns | 133 ns | 6.5× |
| `fitting/arrhenius_6pts` | 233.40 ns | 951 ns | 4.1× |
| `fitting/polynomial_deg2_20pts` | 500.05 ns | 2.19 µs | 4.4× |

Over the 31: the median Cyrius / Rust ratio is 7.7× and the geometric mean 8.2× (1.5× to 69×). The four
largest are bodies of a few multiplies that the port rounds once: `thermo/heat_transfer` (ADR 0002 #27; 5.4 ns before
the overflow audit's B1, 65 ns after), `spectroscopy/absorbance` (#30; 5 → 68 ns) and `gas/ideal_gas_pressure` (#31;
18 → 119 ns) through the correctly rounded monomial, and `spectroscopy/bohr_energy_level`, the exact Rydberg rational
rounded once (#21; 9 → 49 ns). Rust inlines each into two or three multiplies. Every body computes what the golden
corpus checks: the oracle's bits, or the value an ADR accepts; parity, not speed, is the bar.

### The Cyrius-only benchmarks

No Rust counterpart: the port's own benches, by file. 163 rows, 160 distinct benches (`mathx/cr_exp`, `cr_expm1`
and `cr_exp10` are in two files). Where a milestone section gives an ad hoc Rust timing of the same body, it is there.

**`tests/kimiya.bcyr`** (54): the port's own benches beside Rust's: the correctly rounded primitives against the stdlib's, more bodies per module, the ADR 0002 ties and slow paths.

| Bench | Cyrius 2.0.0 |
|---|---|
| `mathx/cr_exp` | 39.2 ns |
| `mathx/cr_expm1` | 49.8 ns |
| `mathx/cr_exp10` | 54.1 ns |
| `mathx/cr_ln` | 48.3 ns |
| `mathx/cr_log10` | 76.8 ns |
| `mathx/cr_pow` | 133 ns |
| `mathx/ganita_f64_pow` | 188 ns |
| `mathx/cr_pow_1.409e-132` | 171 ns |
| `mathx/ganita_f64_pow_1.409e-132` | 207 ns |
| `organic/enthalpy_from_bonds` | 48.4 ns |
| `inorganic/lookup_ionic_radius` | 108 ns |
| `nuclear/decay_chain_two_step` | 254 ns |
| `mathx/cr_sin` | 73.3 ns |
| `mathx/f64_sin` | 36.6 ns |
| `mathx/cr_cos` | 80.5 ns |
| `mathx/f64_cos` | 39.2 ns |
| `mathx/cr_sincos` | 114 ns |
| `mathx/cr_atan2` | 111 ns |
| `mathx/ganita_f64_atan2` | 51.8 ns |
| `mathx/cr_acos` | 138 ns |
| `potential/bond_angle` | 215 ns |
| `potential/lennard_jones` | 37.0 ns |
| `potential/lj_equilibrium_distance` | 103 ns |
| `potential/morse` | 62.4 ns |
| `mathx/cr_exp2_scaled` | 141 ns |
| `potential/distance_3d` | 16.0 ns |
| `potential/distance_3d_2p600` | 315 ns |
| `potential/bond_angle_2p600` | 1.13 µs |
| `potential/dihedral_angle_2p600` | 1.31 µs |
| `potential/bond_angle_linear` | 138 ns |
| `spectroscopy/bohr_transition_energy` | 71.1 ns |
| `spectroscopy/rydberg_wavelength` | 73.7 ns |
| `spectroscopy/series_wavelength_balmer_alpha` | 78.6 ns |
| `spectroscopy/reduced_mass_hcl` | 73.0 ns |
| `spectroscopy/bohr_transition_energy_tie` | 2.44 µs |
| `spectroscopy/reduced_mass_tie` | 2.20 µs |
| `spectroscopy/isotope_pattern_cystine` | 18.3 µs |
| `spectroscopy/isotope_pattern_s8` | 18.4 µs |
| `spectroscopy/fft_spectrum_256` | 29.1 µs |
| `spectroscopy/fft_spectrum_1024` | 137 µs |
| `spectroscopy/fft_spectrum_4096` | 639 µs |
| `spectroscopy/find_peaks_1000` | 19.7 µs |
| `spectroscopy/wavenumber_from_wavelength` | 47.5 ns |
| `stoichiometry/balance_kmno4_hcl` | 2.19 µs |
| `stoichiometry/oxidation_states_fe3o4` | 110 ns |
| `stoichiometry/limiting_reagent` | 33.1 ns |
| `stoichiometry/theoretical_yield_moles` | 115 ns |
| `stoichiometry/percent_yield` | 82.1 ns |
| `stoichiometry/percent_composition_glucose` | 269 ns |
| `stoichiometry/empirical_formula_glucose` | 592 ns |
| `stoichiometry/theoretical_yield_moles_tie` | 132 ns |
| `stoichiometry/percent_yield_subnormal` | 1.63 µs |
| `stoichiometry/limiting_reagent_tie` | 93.5 ns |
| `stoichiometry/empirical_formula_edge` | 7.49 µs |

**`tests/audit.bcyr`** (66): the overflow audit's B1–B4: the rounded-once functions, their edges, ties and slow paths.

| Bench | Cyrius 2.0.0 |
|---|---|
| `thermo/specific_heat_from_calorimetry` | 86.1 ns |
| `phase/boiling_point_elevation` | 66.6 ns |
| `phase/freezing_point_depression` | 66.9 ns |
| `phase/osmotic_pressure` | 140 ns |
| `kinetics/second_order_concentration` | 32.6 ns |
| `kinetics/zero_order_half_life` | 16.8 ns |
| `kinetics/steady_state_intermediate` | 89.4 ns |
| `kinetics/pre_equilibrium_rate` | 126 ns |
| `thermo/heat_transfer_tie` | 369 ns |
| `phase/osmotic_pressure_subnormal` | 4.05 µs |
| `kinetics/second_order_concentration_edge` | 229 ns |
| `kinetics/michaelis_menten_edge` | 126 ns |
| `solution/dilution` | 84.4 ns |
| `solution/solution_conductivity` | 70.8 ns |
| `electrochemistry/faraday_mass_deposited` | 107 ns |
| `electrochemistry/faraday_charge_required` | 118 ns |
| `electrochemistry/faraday_mass_ratio` | 147 ns |
| `spectroscopy/concentration_from_absorbance` | 87.7 ns |
| `spectroscopy/vibrational_frequency` | 85.1 ns |
| `spectroscopy/rigid_rotor_energy` | 109 ns |
| `spectroscopy/rotational_constant` | 115 ns |
| `spectroscopy/diatomic_moment_of_inertia` | 72.4 ns |
| `electrochemistry/nernst_potential_edge` | 197 ns |
| `spectroscopy/vibrational_frequency_edge` | 85.5 ns |
| `spectroscopy/vibrational_frequency_subnormal` | 884 ns |
| `electrochemistry/faraday_mass_ratio_tie` | 393 ns |
| `gas/ideal_gas_volume` | 116 ns |
| `gas/gas_density` | 105 ns |
| `gas/van_der_waals_pressure` | 52.2 ns |
| `gas/peng_robinson_pressure` | 95.1 ns |
| `gas/redlich_kwong_pressure` | 182 ns |
| `potential/coulomb` | 120 ns |
| `potential/coulomb_force` | 141 ns |
| `potential/morse_width_parameter` | 119 ns |
| `inorganic/born_lande_lattice_energy` | 48.2 ns |
| `gas/van_der_waals_pressure_edge` | 257 ns |
| `gas/peng_robinson_pressure_edge` | 602 ns |
| `gas/redlich_kwong_pressure_edge` | 566 ns |
| `potential/lennard_jones_edge` | 204 ns |
| `potential/morse_width_parameter_subnormal` | 1.22 µs |
| `inorganic/born_lande_lattice_energy_edge` | 258 ns |
| `nuclear/binding_energy_per_nucleon_edge` | 114 ns |
| `potential/morse_width_parameter_tie` | 1.23 µs |
| `thermo/hess_law` | 14.9 ns |
| `thermo/enthalpy_from_formation` | 22.7 ns |
| `thermochem/reaction_gibbs_energy` | 298 ns |
| `thermochem/reaction_entropy` | 284 ns |
| `nuclear/q_value` | 25.7 ns |
| `solution/ionic_strength` | 32.4 ns |
| `inorganic/born_haber_lattice_energy` | 14.4 ns |
| `solution/kohlrausch_molar_conductivity` | 28.3 ns |
| `solution/limiting_molar_conductivity` | 23.9 ns |
| `kinetics/zero_order_concentration` | 17.3 ns |
| `thermochem/shomate_enthalpy_antiderivative` | 54.7 ns |
| `thermo/hess_law_edge` | 165 ns |
| `thermo/enthalpy_from_formation_edge` | 211 ns |
| `thermochem/reaction_enthalpy_edge` | 594 ns |
| `nuclear/q_value_edge` | 261 ns |
| `solution/ionic_strength_edge` | 137 ns |
| `inorganic/born_haber_lattice_energy_edge` | 175 ns |
| `reaction/gibbs_free_energy_edge` | 117 ns |
| `solution/kohlrausch_molar_conductivity_edge` | 127 ns |
| `solution/limiting_molar_conductivity_edge` | 132 ns |
| `kinetics/zero_order_concentration_edge` | 140 ns |
| `thermochem/shomate_cp_edge` | 320 ns |
| `thermochem/shomate_enthalpy_antiderivative_edge` | 474 ns |

**`tests/audit2.bcyr`** (33): the overflow audit's B5–B10a: the exp / exp10 arguments, the FFT, `find_peaks`, the flame, unbounded balancing, the unbounded-exponent primitives (`mathx/cr_exp`, `cr_expm1` and `cr_exp10` repeat `tests/kimiya.bcyr`'s as the B10a baseline).

| Bench | Cyrius 2.0.0 |
|---|---|
| `reaction/equilibrium_constant` | 63.4 ns |
| `solution/debye_huckel_limiting` | 86.4 ns |
| `solution/debye_huckel_extended` | 100 ns |
| `solution/davies_activity` | 96.0 ns |
| `reaction/equilibrium_constant_edge` | 151 ns |
| `reaction/equilibrium_constant_subnormal` | 165 ns |
| `solution/debye_huckel_limiting_edge` | 92.5 ns |
| `solution/debye_huckel_extended_edge` | 277 ns |
| `solution/debye_huckel_extended_overflow` | 277 ns |
| `solution/davies_activity_edge` | 162 ns |
| `spectroscopy/fft_spectrum_256_huge` | 306 µs |
| `spectroscopy/fft_spectrum_256_subnormal` | 340 µs |
| `spectroscopy/fft_spectrum_256_squares` | 43.2 µs |
| `spectroscopy/find_peaks_1000_subnormal` | 37.6 µs |
| `thermochem/adiabatic_flame_tiny` | 1.42 µs |
| `stoichiometry/balance_survey_2x3` | 3.61 µs |
| `stoichiometry/balance_overflow_10x9` | 56.2 µs |
| `stoichiometry/balance_dense_30x31` | 157.4 ms |
| `mathx/cr_exp` | 39.0 ns |
| `mathx/cr_exp_ux` | 49.0 ns |
| `mathx/cr_exp_ux_acc` | 1.71 µs |
| `mathx/cr_exp_ux_subnormal` | 76.1 ns |
| `mathx/cr_exp_ux_minus_1e4` | 145 ns |
| `mathx/cr_exp_ux_1e12` | 145 ns |
| `mathx/cr_exp_ux_1e15` | 568 ns |
| `mathx/cr_exp_ux_acc_minus_1e4` | 1.98 µs |
| `mathx/cr_expm1` | 50.1 ns |
| `mathx/cr_expm1_ux` | 59.5 ns |
| `mathx/cr_expm1_ux_1e4` | 152 ns |
| `mathx/cr_exp10` | 54.1 ns |
| `mathx/cr_exp10_ux` | 64.2 ns |
| `mathx/cr_exp10_ux_1e4` | 147 ns |
| `mathx/cr_exp10_ux_acc_minus_1e12` | 1.99 µs |

**`tests/numerics.bcyr`** (10): M5's systems: the integrators, the elimination, the normal draw, the network rank, the gradients.

| Bench | Cyrius 2.0.0 |
|---|---|
| `dynamics/adaptive_first_order` | 5.04 µs |
| `dynamics/stiff_robertson_400steps` | 425 µs |
| `dynamics/stochastic_100steps` | 21.7 µs |
| `numerics/gaussian_elimination_8x8` | 2.15 µs |
| `numerics/pcg32_next_normal` | 177 ns |
| `network/stoichiometric_matrix_h2_br2` | 277 ns |
| `network/rank_michaelis_menten` | 919 ns |
| `network/rank_chapman` | 138 µs |
| `sensitivity/arrhenius_gradient` | 100 ns |
| `sensitivity/compute_gradient_arrhenius` | 496 ns |
## M1 Foundation (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `element/lookup_by_symbol` | `lookup_by_symbol("Fe")`: a linear scan to row 26 | 44.73 ns | 67.0 ns | 1.5× |
| `element/lookup_by_number` | `lookup_by_number(26)`: range check, one table read | 0.84 ns | 4.54 ns | 5.4× |
| `molecule/molecular_weight_glucose` | `molecular_weight()` of C₆H₁₂O₆, built outside the loop: 3 lookups, 3 multiply-adds | 3.63 ns | 26.7 ns | 7.4× |

All three bodies are bit-identical to Rust (golden rows: `element` 975, `molecule` 2,472, 0 ulp).

**Two optimizations, measured.** The first version of each body measured 92.4 ns, 7.03 ns and
121 ns:
- `molecular_weight` read each atom through `atom_element_number(k)` / `atom_count(k)`, each a
  bounds check, a `vec_get` with its own two checks, and calls. It now walks the vec's buffer
  directly (lib/vec.cyr documents the `{ data, len, cap }` layout) and reads the mass column at
  the checked row: 121 → 57 ns with one `vec_get` per word, → 33 ns walking the buffer.
- `element_lookup_by_number` called `_km_element_row` for its range check; the check is inline
  now: 7.03 → 4.54 ns, and `molecular_weight`, which makes three of them per call, 33 → 26.7 ns.
- `element_lookup_by_symbol` compares the first byte before calling `memeq`, with the probe's
  first byte read once: 92.4 → 67.0 ns.

**Reading the gap.** Rust inlines `lookup_by_number` into a compare and a load, and unrolls the
glucose sum; Cyrius compiles each call as a call (the lookup inside `molecular_weight`, the
`Result` pair return). The symbol scan is the closest, because both sides do the same linear
string search and Rust's per-row cost is a real comparison too. Parity, not speed, is M1's bar.

## M2a reaction, thermo, phase, gas (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `reaction/gibbs_free_energy` | `gibbs_free_energy(-100000, 298, -200)`: three guards, ΔH − T·ΔS, a `Result` | 1.75 ns | 16.9 ns | 9.7× |
| `thermo/heat_transfer` | `heat_transfer(1000, 4.184, 10)`: two multiplies | 0.95 ns | 5.12 ns | 5.4× |
| `gas/ideal_gas_pressure` | `ideal_gas_pressure(1, 300, 0.025)`: three guards, nRT/V, a `Result` | 2.18 ns | 16.7 ns | 7.7× |
| `gas/peng_robinson_z_n2_1mpa` | `peng_robinson_z(PR_NITROGEN, 300 K, 1 MPa)`: validation, a·α and b, A and B, the one-root (Cardano) branch, two Newton steps | 61.39 ns | 270 ns | 4.4× |

All four bodies are bit-identical to Rust (golden rows: `reaction` 1,188, `thermo` 1,107, `gas` 19,717; the PR
roots 0 ulp on every target).

**One optimization, measured.** The NaN/±∞ guards were three nested calls each (`_km_is_positive` →
`_km_is_finite` → `_km_bits`). `_km_is_finite(x)` is now `x − x == 0` and `_km_is_positive` /
`_km_is_non_negative` compare against +∞ directly: `gibbs_free_energy` 23.3 → 16.9 ns (−27%),
`ideal_gas_pressure` 26.7 → 16.7 ns (−37%), `peng_robinson_z_n2_1mpa` 289 → 270 ns (−7%). The M1 benches don't
use the guards and are unchanged (65.3, 4.48, 26.7 ns).

**Reading the gap.** As in M1, Rust inlines each function into a few instructions; Cyrius makes the calls and
returns the `Result` as a pair. In `peng_robinson_z`, ganita's correctly rounded `cbrt` alone takes about 99 ns
(Rust links CORE-MATH's). The three-root branch, not benched in Rust, adds the port's correctly rounded acos
and three cos calls (ADR 0005, double-double): the price of bit-identical roots at the spinodals.

## M2b thermochem (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `thermochem/reaction_enthalpy_ch4` | `reaction_enthalpy([CO2(g) 1, H2O(l) 2], [CH4(g) 1, O2(g) 2])`: four table scans, two sums (and, ADR 0002 #1, one × 1000) | 38.58 ns | 250 ns | 6.5× |
| `thermochem/vant_hoff_k` | `vant_hoff_k(100, -50000, 298.15, 400)`: four guards, one exp | 9.51 ns | 58.9 ns | 6.2× |
| `thermochem/enthalpy_change_cp` | `enthalpy_change_cp("CO2(g)", 300, 600)`: the segment lookup, guards, one segment's antiderivative twice | 24.00 ns | 166 ns | 6.9× |
| `thermochem/shomate_cp` | `s.cp(500)`, s = CO2's first segment (looked up once, outside the loop) | 2.34 ns | 14.5 ns | 6.2× |
| `thermochem/adiabatic_flame_ch4_air` | `adiabatic_flame_temperature(-802300, [CO2(g) 1, H2O(g) 2, N2(g) 7.52], 298.15)`: three segment lookups, the bracket, Newton (3 updates on this input) | 301.48 ns | 3.24 µs | 10.7× |

All five bodies are bit-identical to Rust (golden rows: `thermochem` 9,487, 0 ulp on every class-B row; `vant_hoff_k`
within K = 6).

**Two optimizations tried, neither kept.** Testing the first byte inline before each Shomate formula compare (as
`element_lookup_by_symbol` does) left `enthalpy_change_cp` at 164–166 ns: CO2 is the third row, and the cost is the
`strlen` (16 ns) and the two `memeq` calls (13 ns each) that the run's length needs. Evaluating the antiderivative
straight from the columns, without the six-argument call to the shared `_km_shomate_h5`, took it from 22 to 20 ns
(the flame makes about 90 per call): not worth a second copy of the formula. Of `enthalpy_change_cp`'s 166 ns, the
segment lookup is 70 and the two-antiderivative integral 69.

**Reading the gap.** As in M1 and M2a: Rust inlines the Shomate polynomial, the `Segments` walk and the closures
into straight-line code; Cyrius makes each a call and returns the `Result` as a pair. The flame temperature is the
widest gap because it multiplies the per-call cost: its 3 µs are about six balance evaluations (3 products × up to
3 segments × 2 antiderivatives each) and three derivative evaluations. The port keeps Rust's resolved-segments table
on the stack (32 products), since the bump allocator can't free and the function runs in loops.

## M3a solution, electrochemistry (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `solution/ph_from_h` | `ph_from_h_concentration(1e-4)`: one guard, one log10 (`_km_log10`, fdlibm) | 9.73 ns | 66.2 ns | 6.8× |
| `solution/weak_acid_ph` | `weak_acid_ph(1.8e-5, 0.1)`: two guards, the charge-balance Newton (4 updates on this input), one log10 | 34.85 ns | 182 ns | 5.2× |
| `electrochemistry/lookup_half_reaction` | `lookup_half_reaction("Cu2+/Cu")`: a linear scan to row 27 | 57.06 ns | 116 ns | 2.0× |
| `electrochemistry/nernst_potential` | `nernst_potential(0.342, 2, 298.15, 0.1)`: four guards, one ln — the correctly rounded `_km_ln` | 6.42 ns | 1.94 µs | 302× |

All four bodies are bit-identical to Rust (golden rows: `solution` 28,566, `electrochemistry` 18,574; the Nernst rows 0 ulp
on every target, the pH rows within K = 4).

**The price of the correctly rounded ln.** `nernst_potential` cancels where E° ≈ (RT/nF)·ln Q (S ≈ 9 × 10¹⁸ with the
stdlib ln, port-parity.md §7), so it calls `_km_ln` (ADR 0005, amendment M3a), the double-double atanh series: measured
alone, `_km_ln(0.1)` takes 1.90 µs against 23 ns for `f64_ln`, so the Nernst body would be about 45–60 ns (7–9× Rust)
with the stdlib ln, and is 1.94 µs with the correctly rounded one. `henderson_hasselbalch` (not benched in Rust) pays the
same through `_km_log10_ln` (1.85 µs against 59 ns for `_km_log10`); the four S = 1 log10 functions keep the fast
`_km_log10`. A two-phase (Ziv) correctly rounded log — a fast evaluation with a rounding test, the double-double series
only when the test fails — would recover most of the cost; it is a possible follow-up, not in M3a.

**Reading the gap.** As before, Rust inlines each body; Cyrius makes the guard, primitive and `Result` calls. The lookup is
the closest, because both sides scan the same couple strings (the port tests the first byte before `memeq`, as the
thermochem and element lookups do). `ph_from_h`'s 66 ns are mostly `_km_log10` (59 ns measured alone).

## M3b kinetics (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `kinetics/arrhenius_rate` | `arrhenius_rate(1e13, 50000, 300)`: three guards, −Ea/(RT), one exp (`f64_exp`), A·e | 6.78 ns | 49.4 ns | 7.3× |
| `kinetics/eyring_rate` | `eyring_rate(80000, 298)`: two guards, k_B·T/h, one exp, a multiply | 6.73 ns | 42.4 ns | 6.3× |
| `kinetics/michaelis_menten` | `michaelis_menten(100, 5, 10)`: three guards, Vmax·[S]/(Km + [S]) | 2.92 ns | 16.5 ns | 5.7× |

`michaelis_menten` is class B and bit-identical to Rust (0 ulp on all 370 golden rows). `arrhenius_rate` and
`eyring_rate` are class C and measured, not yet gated (port-parity.md §8: a subnormal exp result scaled by a prefactor
puts their S past the rule's limit; ADR 0006 will set their comparison); on these benched inputs the exp results are
normal, and both bodies return Rust's bits on x86_64 and aarch64.

**Reading the gap.** As before, Rust inlines each body into a few instructions; Cyrius makes the guard calls and returns
the `Result` as a pair. The exp is `f64_exp` (x87 on x86_64); `michaelis_menten`'s three guards and arithmetic take
16.5 ns, so the exp and its call are most of the rest of the Arrhenius and Eyring bodies.

## M3c correct rounding (2026-11-09)

ADR 0006: every transcendental the modules call is now correctly rounded (`_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`,
`_km_cr_ln`, `_km_cr_log10`, `_km_cr_pow`; `_km_acos` / `_km_cos` unchanged). Every bench, before and after the switch,
with Rust's alongside. *Doc* is the number recorded at the bench's milestone above; *before* and *after* were measured
the same day on the same host, `cyrius bench` just before and just after M3c-3's switch (the host ran about 5–10 % slower
than at the earlier milestones, so compare *before* with *after*). Per-op means; ns unless marked.

| Bench | Rust | Cyrius, doc | Cyrius, before | Cyrius, after | after / Rust | What changed |
|---|---|---|---|---|---|---|
| `element/lookup_by_symbol` | 44.73 | 67.0 | 85.9 | 70.2 | 1.6× | nothing (noise) |
| `element/lookup_by_number` | 0.84 | 4.54 | 4.8 | 4.7 | 5.6× | nothing |
| `molecule/molecular_weight_glucose` | 3.63 | 26.7 | 28.8 | 29.1 | 8.0× | nothing |
| `reaction/gibbs_free_energy` | 1.75 | 16.9 | 18.6 | 18.8 | 10.7× | nothing |
| `thermo/heat_transfer` | 0.95 | 5.12 | 5.6 | 5.6 | 5.9× | nothing |
| `gas/ideal_gas_pressure` | 2.18 | 16.7 | 18.4 | 18.1 | 8.3× | nothing |
| `gas/peng_robinson_z_n2_1mpa` | 61.39 | 270 | 297 | 287 | 4.7× | nothing (the one-root branch: cbrt and sqrt) |
| `thermochem/reaction_enthalpy_ch4` | 38.58 | 250 | 276 | 267 | 6.9× | nothing |
| `thermochem/vant_hoff_k` | 9.51 | 58.9 | 65.5 | **77.1** | 8.1× | exp: `f64_exp` → `_km_cr_exp` |
| `thermochem/enthalpy_change_cp` | 24.00 | 166 | 183 | 171 | 7.1× | nothing |
| `thermochem/shomate_cp` | 2.34 | 14.5 | 15.6 | 15.5 | 6.6× | nothing |
| `thermochem/adiabatic_flame_ch4_air` | 301.48 | 3.24 µs | 3.49 µs | 3.43 µs | 11.4× | nothing |
| `solution/ph_from_h` | 9.73 | 66.2 | 69.9 | **87.3** | 9.0× | log10: `_km_log10` (fdlibm) → `_km_cr_log10` |
| `solution/weak_acid_ph` | 34.85 | 182 | 196 | **219** | 6.3× | the same log10 |
| `electrochemistry/lookup_half_reaction` | 57.06 | 116 | 124 | 120 | 2.1× | nothing |
| `electrochemistry/nernst_potential` | 6.42 | 1.94 µs | 86.7 | 83.9 | 13.1× | ln: already `_km_cr_ln` through the `_km_ln` wrapper since M3c-2a (1.94 µs → 87 ns then); now called directly |
| `kinetics/arrhenius_rate` | 6.78 | 49.4 | 52.0 | **65.5** | 9.7× | exp: `f64_exp` → `_km_cr_exp` |
| `kinetics/eyring_rate` | 6.73 | 42.4 | 44.6 | **59.0** | 8.8× | the same exp |
| `kinetics/michaelis_menten` | 2.92 | 16.5 | 17.4 | 17.3 | 5.9× | nothing |
| `mathx/cr_exp` (−50000/(8.314·300)) | — | — | 39.8 | 38.7 | — | Cyrius-only (M3c-2a) |
| `mathx/cr_expm1` (ln 2) | — | — | 51.3 | 49.6 | — | |
| `mathx/cr_exp10` (−4.75) | — | — | 54.4 | 53.5 | — | |
| `mathx/cr_ln` (0.1) | — | — | 49.0 | 49.1 | — | |
| `mathx/cr_log10` (10⁻⁴) | — | — | 77.4 | 78.7 | — | |
| `mathx/cr_pow` (0.1, 1.5) | — | — | 136 | 138 | — | Cyrius-only (M3c-2b) |
| `mathx/ganita_f64_pow` (0.1, 1.5) | — | — | 192 | 198 | — | what the modules called until M3c-3 |
| `mathx/cr_pow_1.409e-132` (1.409·10⁻¹³², 2.451) | — | — | 175 | 181 | — | |
| `mathx/ganita_f64_pow_1.409e-132` | — | — | 204 | 218 | — | |

**The price of correct rounding.** The switch costs 12–15 ns per exp call in a function body (`vant_hoff_k` +11.6,
`arrhenius_rate` +13.5, `eyring_rate` +14.4) and about 20 ns per log10 (`ph_from_h` +17.4, `weak_acid_ph` +23): measured
alone, `_km_cr_exp` takes 38–40 ns against about 11 ns for `f64_exp` (x87) and `_km_cr_log10` 77–79 ns against 59 ns for
`_km_log10` (M3a). pow got faster: `_km_cr_pow` is 136–181 ns against ganita's 192–218 (`nth_order_half_life` is not
benched). Redlich–Kwong's Tc^2.5 went from `_km_pow_2_5` (85 ns, double-double) to `_km_cr_pow(tc, 2.5)` (130 ns); it is
not benched in Rust. `nernst_potential`'s 1.94 µs → 84 ns (23× faster) is M3c-2a's, ADR 0005's double-double ln
replaced by the two-phase one. The benches whose bodies have no transcendental are unchanged within the host's noise.
ADR 0006 asks the fast paths to stay "within a small factor of the stdlib primitives": exp is 3.5×, log10 1.3×, pow 0.7×.

## M4b organic, inorganic (2026-11-09)

`rust-old/benches/benchmarks.rs` has no organic or inorganic bench, so these two are Cyrius-only: there is no Rust
baseline to compare against (they are recorded for later milestones to compare with).

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `organic/enthalpy_from_bonds` | the integration test's CH₄ + 2O₂: `enthalpy_from_bonds([(CH, 4), (OODouble, 2)], [(CODouble, 2), (OH, 4)])`, two lists of two records (each value and count checked), two sums from −0.0, one subtraction | — | 49.5 ns | — |
| `inorganic/lookup_ionic_radius` | `lookup_ionic_radius("Sn", 2)`: a linear scan to the last of the 30 rows (the i8 charge and the first byte decide each row) | — | 107 ns | — |

The same run's other benches were 5–10 % above their M3c figures (`vant_hoff_k` 83.4 ns, `arrhenius_rate` 70.6 ns,
`lookup_half_reaction` 127 ns), none of whose code changed: the host was shared with another build that day, so read
these two numbers with that margin.

## M4c nuclear (2026-11-09)

`rust-old/benches/benchmarks.rs` has two nuclear benches, ported with the same inputs; the decay chain is Cyrius-only.

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `nuclear/binding_energy_per_nucleon_fe56` | `binding_energy_per_nucleon(26, 30, 55.934_936_33)`: the u32 checks and checked add, A = 0, the mass guard, Δm, × m_u·c², one division, a `Result` | 2.02 ns | 27 ns | 13× |
| `nuclear/semi_empirical_fe56` | `semi_empirical_binding_energy(26, 56)`: the checks, one cbrt (ganita's, correctly rounded), one √, the five terms; a `Result` in the port (ADR 0002 #3) | 18.50 ns | 131 ns | 7.1× |
| `nuclear/decay_chain_two_step` | the Rust test's `decay_chain_two_step(1000, 0.1, 0.05, 10)`: three guards, two `_km_cr_exp` and two `_km_cr_expm1` (the Bateman branch), into the caller's buffer | — | 256 ns | — |

`binding_energy_per_nucleon` first measured 41 ns calling `mass_defect` (which checks its u32 words again) and
`binding_energy_mev`; it now calls a private, unchecked `_km_mass_defect` and multiplies by `KM_AMU_MEV` inline: 29 ns
back to back on the same host, 27 ns in the final run (the corpus pins the bits: 618 rows, 0 ulp). The SEMF's cost is
mostly ganita's cbrt. The same run's other benches were within 5 % of their M4b and M3c figures.


## M4a correctly rounded trigonometry (2026-11-09)

Cyrius-only (`mathx/*`, no Rust counterpart: Rust calls glibc's): each correctly rounded primitive beside the stdlib
function it replaces or could replace. sin and cos at 2π/3 + 0.1 (a Peng–Robinson angle), atan2 at (0.866, −0.5) (a 120°
bond angle), acos at −0.35.

| Bench | Correctly rounded | Stdlib | Ratio |
|---|---|---|---|
| sin | `mathx/cr_sin` 74 ns | `f64_sin` 37 ns (x87 on x86_64; not correctly rounded) | 2.0× |
| cos | `mathx/cr_cos` 82 ns | `f64_cos` 40 ns | 2.1× |
| sin and cos | `mathx/cr_sincos` 115 ns | — | — |
| atan2 | `mathx/cr_atan2` 113 ns | `ganita_f64_atan2` 53 ns | 2.1× |
| acos | `mathx/cr_acos` 140 ns | ADR 0005's double-double `_km_acos` 5.8 µs | 0.024× |

ADR 0005's double-double `_km_cos` was 1.68 µs; M4d moved gas.cyr to `_km_cr_cos` and `_km_cr_acos` and removed both.

## M4d potential (2026-11-09)

`rust-old/benches/benchmarks.rs` has one potential bench, ported with the same input; the others are Cyrius-only.

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `potential/dihedral_angle` | Rust's butane-like gauche geometry in ångströms, A (1, 1, 0), B (0, 0, 0), C (1.5, 0, 0), D (2.5, 0, 1), through four pointers: three sub3, two cross, seven dot3, the two collinearity tests, one √ and one atan2, which here is the x = 0 special case (+90°) | 20.35 ns | 130 ns | 6.4× |
| `potential/bond_angle` | water at metre scale (O–H 0.9572 Å, 104.52°): two sub3, a cross, six dot3, the coincidence tests, √ and a full `_km_cr_atan2` (fast path) | — | 215 ns | — |
| `potential/lennard_jones` | Ar at r ≈ r_min: three guards, σ/r, five products, the square, (4·ε)·d | — | 35 ns | — |
| `potential/lj_equilibrium_distance` | Ar's σ: `_km_cr_sixth_root2_mul`, a triple-double product and one `_km_cr_finish` (ADR 0002 #16) | — | 109 ns | — |
| `potential/morse` | HCl at 1.1·r_e: two guards, one `_km_cr_exp`, (D_e·x)·x | — | 66 ns | — |
| `mathx/cr_exp2_scaled` | 1000 after 1075.3 half-lives, a subnormal result (ADR 0002 #13's N₀·2⁻ⁿ): exp's fast core and rounding test at the subnormal grid | — | 151 ns | — |

The same run's M4a trig benches matched the M4a figures within 2 % (`cr_sin` 74 ns, `cr_cos` 81, `cr_atan2` 111, `cr_acos`
142), and the other benches were within 5 % of their M4c figures but `thermochem/adiabatic_flame_ch4_air` (3.80 µs against
3.24; nothing on its path changed, + − × ÷ only) and `gas/peng_robinson_z_n2_1mpa` (283 ns against 270; the one-root
branch, which calls no trig): host noise, as M4b's note describes.

## M4e spectroscopy (2026-11-09)

`rust-old/benches/benchmarks.rs` has three spectroscopy benches, ported with the same inputs; the others are Cyrius-only
(Rust has no FFT, peak-detection or wavenumber bench).

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `spectroscopy/absorbance` | `absorbance(100, 1, 0.01)`: two multiplies | 0.98 ns | 11 ns | 11× |
| `spectroscopy/bohr_energy_level` | `bohr_energy_level(1, 3)`: two range checks, ((−R·Z)·Z)/(n·n), a `Result` | 1.21 ns | 16 ns | 13× |
| `spectroscopy/isotope_pattern_glucose` | C₆H₁₂O₆: 24 convolution steps (the scratch at most 30 records), each a stable insertion sort, the binning, the renormalisation and the pruning, on the two stack tables (nothing allocated); the result into a 64-record buffer | 2.80 µs | 16.8 µs | 6.0× |
| `spectroscopy/fft_spectrum_256` | a seeded signal of 256 samples, in place (cap = 2n): hisab 1.4.0's radix-2 loop, 8 `_km_cr_sincos` twiddles, 129 magnitudes | — | 28.2 µs | — |
| `spectroscopy/fft_spectrum_1024` | the same, 1,024 samples | — | 134 µs | — |
| `spectroscopy/fft_spectrum_4096` | the same, 4,096 samples | — | 628 µs | — |
| `spectroscopy/find_peaks_1000` | 1,000 values at threshold 0.1: the finite min / max pass, the interpolated threshold, the window pass (cap 1000: written in one pass) | — | 22.5 µs | — |
| `spectroscopy/wavenumber_from_wavelength` | Hα (656.47 nm): the guard and ADR 0002 #19's `_km_cr_hundredth_div` fast path (a TwoProduct residual and a rounding test) | — | 58 ns | — |

The FFT costs about 26 ns per butterfly (n/2·log₂n of them: 1,024 at n = 256, 24,576 at 4,096), each 4 multiplies and 6
adds through memory plus the twiddle recurrence's 4 multiplies and 2 adds; the twiddles themselves (one `_km_cr_sincos`,
116 ns, per stage) are under 0.5 % of it. The isotope pattern's 6× is the same shape as the other loops' (each operation
a load and a store); Rust's 2.80 µs is its post-audit sort and binning.

## M4e-2 potential's geometry with an unbounded exponent (2026-11-10)

ADR 0002 #20 adds a range test to Rust's chain (the fast path) and an unbounded-exponent re-evaluation where it fails.
Before → after, the same host, back to back (`cyrius bench tests/kimiya.bcyr`):

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `potential/dihedral_angle` | Rust's gauche geometry: the chain, and the window test on b₁, b₂, b₃ (each component's square against 2⁻⁴⁵⁸ with its zero test, the three \|b\|² against 2⁵⁰⁶; the squares computed inline instead of three `_km_dot3` calls) | 20.35 ns | 131 ns | 133 ns |
| `potential/bond_angle` | water at metre scale: \|ba\|², \|bc\|², the two maxima and \|ba × bc\|² against [2⁻⁹⁰⁰, DBL_MAX] (in place of Rust's two `is_positive` calls), ba · bc finite | — | 213 ns | 213 ns |
| `potential/distance_3d` | water's H–H distance in metres: s against [2⁻⁹⁰⁰, DBL_MAX] | — | 14 ns | 16 ns |
| `potential/distance_3d_2p600` | the same points × 2⁶⁰⁰: the squares overflow, so the slow path (`_km_distance_3d_ux`: six conversions, three subtractions, a dot product, the final √) | — | 14 ns (Rust's +∞) | 314 ns |
| `potential/bond_angle_2p600` | water × 2⁶⁰⁰: the slow path (`_km_bond_angle_ux`: five dot products, a cross product, the final atan2) | — | 52 ns (Rust's Err) | 1.13 µs |
| `potential/dihedral_angle_2p600` | the gauche geometry × 2⁶⁰⁰: the slow path (`_km_dihedral_angle_ux`) | — | 80 ns (Rust's Err) | 1.31 µs |
| `potential/bond_angle_linear` | CO₂ in metres on an axis: ba × bc is exactly 0, so \|ba × bc\|² = 0 fails the window; the linear shortcut (each component of ba and bc 0 or at least 2⁻⁵¹¹, then atan2(+0, ba · bc)) | — | 1.06 µs (the slow path; a 2M-call loop, before the shortcut) | 135 ns |

The slow path runs only where Rust's chain leaves the normal range (no physical geometry; an exactly linear one, whose
\|ba × bc\|² is 0, takes the shortcut above instead); each unbounded operation is a call on (significand, exponent) pairs
and a renormalisation, about 25 ns. The other benches were within noise of the
before run (nothing else changed). A first version tested every coordinate's bit pattern against a window through calls
(+19–21 ns on the two angles); the terminal sums for `distance_3d` and `bond_angle` and the inline squares for the
dihedral brought it to the figures above.

## M4e-3 the exact Rydberg rationals and the exact reduced mass (2026-11-10)

ADR 0002 #21 and #22 round the Bohr levels, the transition energy, the Rydberg formula, the series and the reduced mass
once: a double-double fast path with a rounding test, and the exact multi-word slow path (`_km_cr_ratio`) where the test
fails (the ties, and about 2⁻⁴⁶ of other arguments). Before → after, the same host, back to back (a trimmed copy of
`tests/kimiya.bcyr` with the HEAD `src/spectroscopy.cyr` for "before"):

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `spectroscopy/bohr_energy_level` | `bohr_energy_level(1, 3)`: R·Z² split exactly (Veltkamp), the quotient by n² < 2²⁶ inline with an exact residual and the rounding test, a `Result` | 1.21 ns | 9 ns | 49 ns |
| `spectroscopy/bohr_transition_energy` | `bohr_transition_energy(1, 2, 3)`: R·Z²·(n₂² − n₁²) a TwoProduct, (n₁n₂)² an exact integer (n₂ < 2¹³), `_km_dd_div_cr` | — | 27 ns | 73 ns |
| `spectroscopy/rydberg_wavelength` | `rydberg_wavelength(1, 2, 3)` (Hα with R∞): the same operands, the reciprocal quotient | — | 24 ns | 73 ns |
| `spectroscopy/series_wavelength_balmer_alpha` | `spectral_series_wavelength(Balmer, 1)` (R_H) | — | 28 ns | 78 ns |
| `spectroscopy/reduced_mass_hcl` | `reduced_mass(1.007825, 34.968853)` (¹H³⁵Cl in u): the exponent gap, s₁s₂ a TwoProduct, s₁ + s₂·2⁻ᵍ a FastTwoSum, `_km_dd_div_cr` | — | 14 ns | 73 ns |
| `spectroscopy/bohr_transition_energy_tie` | `bohr_transition_energy(255, 58055, 524544)`, an exact midpoint: the fast path declines, `_km_cr_ratio` on 134- and 128-bit integers | — | 28 ns | 2.42 µs |
| `spectroscopy/reduced_mass_tie` | `reduced_mass(2²⁷ − 1, 2²⁷ + 1)`, an exact midpoint: `_km_cr_ratio` on 106- and 108-bit integers | — | 14 ns | 2.20 µs |

Each f64 operation costs about a nanosecond in this compiler's code and a call about 4 ns (`_km_two_prod` alone 14 ns), so
the fast paths' 25–45 operations set the floor; the first version (TwoProducts through calls everywhere) took 75, 148, 151,
155 and 91 ns, and the exact splits (R·Z² by Veltkamp, n² and (n₁n₂)² as integers below 2⁵², an inline TwoProduct in
`_km_dd_div_cr`) brought it to the figures above. The slow path is 55 steps of bit-wise long division over four or five
32-bit limbs; it runs on the ties and, for the reduced mass, on every result that may be subnormal (masses below 2⁻¹⁰²¹).
The other benches were within noise of the before run.

## M4e-4 the cited isotope table (2026-11-10)

ADR 0002 #23 changes `isotope_abundances`' data (the full AME2020 masses, ³⁶S added; the algorithm is Rust's). Two
sulfur-bearing benches are new and Cyrius-only (rust-old/benches has only glucose, which has no sulfur). The old table
(the HEAD `src/spectroscopy_data.cyr`) against the new, the same host, interleaved runs of a trimmed copy of
`tests/kimiya.bcyr` (each figure the mean of two runs of 20 batches × 10,000 calls):

| Bench | Body | Rust | Old table | New table |
|---|---|---|---|---|
| `spectroscopy/isotope_pattern_glucose` | C₆H₁₂O₆: 24 convolution steps, the scratch at most 30 records (unchanged: no sulfur) | 2.80 µs | 15.8 µs | 15.8 µs |
| `spectroscopy/isotope_pattern_cystine` | C₆H₁₂N₂O₄S₂: 26 steps, two of them sulfur's; the scratch at most 33 records, now 44 | — | 16.9 µs | 18.5 µs |
| `spectroscopy/isotope_pattern_s8` | S₈: 8 sulfur steps; up to 15 peaks a step and 45 scratch records, now 17 and 68 | — | 8.4 µs | 18.8 µs |

Glucose's operations are the same with either table (its masses change, not its peak counts), and it measured the same.
Sulfur's fourth isotope costs in two ways: each step convolves 4 isotopes instead of 3, and ³⁶S's +4 offset widens the
pattern and puts each ³⁶S record further from its sorted place, so the insertion sort moves it further (S₈ 2.2×, cystine
+10 %). Every scratch stays within the 192-record stack tables, so nothing is allocated; the heap path starts past them
(the corpus's S₃₀₀₀ reaches 1,276 records, 951 with Rust's table).

## M4f stoichiometry (2026-11-10)

ADR 0002 #7's API: `balance_equation` on an integer matrix (Rust's bench matrix as integers), `assign_oxidation_states`
with rational states. Rust's two stoichiometry benches keep their names and bodies; the others are Cyrius-only. Before →
after the proven fast paths (an elimination step with four operands below 2³¹ formed directly; each verification row
summed in i64 first; `reduce_row`'s fold stopping at a gcd of 1; the oxidation sums' checks inline), the same host, back
to back (a trimmed copy of `tests/kimiya.bcyr` with only these benches):

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `stoichiometry/balance_ch4_combustion` | `balance_equation` of CH₄ + 2O₂ → CO₂ + 2H₂O (3 × 4): the entry check, 3 pivots of fraction-free elimination with gcd reductions, the null vector from the lcm of the pivots, the sign and u32 checks, the verification | 233.90 ns | 1.20 µs | 770 ns |
| `stoichiometry/oxidation_states_kmno4` | `assign_oxidation_states([(19, 1), (25, 1), (8, 4)], 0)`: three lookups, the rules, Mn solved as +7/1, three 24-byte records | 22.38 ns | 156 ns | 122 ns |
| `stoichiometry/balance_kmno4_hcl` | 2KMnO₄ + 16HCl → … (5 × 6) | — | 3.84 µs | 2.16 µs |
| `stoichiometry/oxidation_states_fe3o4` | Fe₃O₄: Fe +8/3, the gcd-reduced rational (Rust's Err path) | — | 132 ns | 110 ns |
| `stoichiometry/limiting_reagent` | (3, 2) mol over the coefficients (2, 1): two guards and a division per reactant | — | 31 ns | 31 ns |
| `stoichiometry/theoretical_yield_moles` | the same, × 2 | — | 46 ns | 46 ns |
| `stoichiometry/percent_yield` | 8 of 10 | — | 14 ns | 14 ns |
| `stoichiometry/percent_composition_glucose` | C₆H₁₂O₆ as a `KmMolecule`: `molecular_weight`, three records | — | 70 ns | 71 ns |
| `stoichiometry/empirical_formula_glucose` | 40.00 / 6.71 / 53.29 %: the moles, the minimum, the ratios, the multiplier search (k = 1), three records | — | 217 ns | 217 ns |

All bodies are Rust's results (golden `stoichiometry`: 41,435 rows; the float functions 0 ulp). A phase profile of
`balance_ch4_combustion` (early returns, before the last tweaks): the entry check and copy 74 ns, the elimination about
385 ns, the null vector 165 ns, the sign / u32 checks and the verification 150 ns — about a nanosecond per operation and
4 ns per call in this compiler's code, as M4e-3 measured, against Rust's inlined loop. The local work tables are not the
cost (a 16-word table measured the same as the 512-word one).

## M4f-2 exact stoichiometry (2026-11-10)

ADR 0002 #24–#26: the yields and the percent composition rounded once (`_km_cr_mul_div`), the limiting reagent's exact
comparison, the exact empirical formula. All Cyrius-only (rust-old/benches has no bench of these five). Before → after,
the same host, back to back (the trimmed copy of `tests/kimiya.bcyr` with only the stoichiometry benches):

| Bench | Body | Before | After |
|---|---|---|---|
| `stoichiometry/limiting_reagent` | (3, 2) mol over (2, 1): the rounded ratios differ, so they decide (Rust's own path) | 32 ns | 33 ns |
| `stoichiometry/theoretical_yield_moles` | the same, × 2: RN(m·p/c), the fast path (the test passes) | 46 ns | 115 ns |
| `stoichiometry/percent_yield` | 8 of 10: RN(100·a/t), the fast path | 14 ns | 81 ns |
| `stoichiometry/percent_composition_glucose` | C₆H₁₂O₆: three RN(mass·count·100/mw) | 72 ns | 276 ns |
| `stoichiometry/empirical_formula_glucose` | 40.00 / 6.71 / 53.29 %: the exact minimum, three double-double ratios, three fast tolerance tests at k = 1, the counts | 217 ns | 586 ns |
| `stoichiometry/theoretical_yield_moles_tie` | y·3/2 with 3y = 2⁵³ + 1, a midpoint: the test fails, `_km_mul_div_tie` resolves it exactly | — | 133 ns |
| `stoichiometry/percent_yield_subnormal` | 100·3·2⁻¹⁰⁷⁴/200, a subnormal tie: the exact path, `_km_cr_ratio` | — | 2.15 µs |
| `stoichiometry/limiting_reagent_tie` | (1, 2) over (1, 2): equal rounded ratios, `_km_ratio_lt` | — | 94 ns |
| `stoichiometry/empirical_formula_edge` | C 12.011 %, O 17.5989 %: the fast test undecided at k = 1, the exact one twice (the multiplier, the count) | — | 7.46 µs |

The fast path of RN(x·y/z) is about 45 floating-point operations (two TwoProducts, the double-double quotient, its
rounding test), the same work as M4e-3's reduced mass (14 → 73 ns); inlining its helper calls was measured at 81 → 73 ns
for `percent_yield` and not kept. Midpoints, common with integer coefficients (6.6 % of `cr_verify.py`'s typical
theoretical yields), cost 133 rather than 115 ns: the exact resolution stays in the fast path, so only subnormal operands
and results near the range ends reach `_km_cr_ratio`. The empirical formula's exact test is a 36-bit restoring division
on 8-limb integers (about 3.5 µs a test); it runs only for knife-edge inputs, within 2⁻⁵⁰ of n ± 1/10.

## The overflow audit's B1 (2026-11-10)

ADR 0002 #27 and #28: thermo's calorimetry and phase's colligative products rounded once, kinetics' monomials rounded
once, `second_order_concentration` and `michaelis_menten` with an unbounded exponent behind a range test (Rust's chain in
range), on `src/mathx.cyr`'s correctly rounded monomial (M4f-2's `_km_cr_mul_div` generalised; a straight-line fast path
per shape). Before: the full `cyrius bench` on the code before B1, with B1's new benches added first; after: a trimmed
copy of `tests/kimiya.bcyr` with only these benches (it reproduces the full run's numbers), the same host, back to back.
Two benches have a Rust counterpart; the rest are Cyrius-only.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `thermo/heat_transfer` | `heat_transfer(1000, 4.184, 10)`: RN(m·c·ΔT) by `_km_mul3_fast` (two TwoProducts, the rounding test, the scale) | 0.95 ns | 5.4 ns | 65 ns |
| `thermo/specific_heat_from_calorimetry` | (4184 J, 1000 g, 1 K): three guards, RN(q/(m·ΔT)) by `_km_div2_fast` | — | 18 ns | 86 ns |
| `phase/boiling_point_elevation` | (0.512, 1, 2): RN(K·m·i) | — | 6.1 ns | 66 ns |
| `phase/freezing_point_depression` | (1.86, 0.5, 2): RN(K·m·i) | — | 6.4 ns | 67 ns |
| `phase/osmotic_pressure` | (2, 0.1 M, 298.15 K): three guards, RN(i·M·R·1000·T) by `_km_mul5_fast` (four double-double products) | — | 20 ns | 139 ns |
| `kinetics/zero_order_half_life` | (1, 0.1): two guards, Rust's division (2k finite, so already correctly rounded) | — | 16 ns | 16 ns |
| `kinetics/steady_state_intermediate` | (0.1, 1, 10): three guards, RN(k₁·[A]/k₂) by `_km_mul_div_fast` | — | 21 ns | 88 ns |
| `kinetics/pre_equilibrium_rate` | (10, 2, 0.1, 1): four guards, RN(k₂·k₁·[A]/k₋₁) by `_km_mul3_div_fast` | — | 27 ns | 125 ns |
| `kinetics/second_order_concentration` | (1, 0.1, 10): three guards, Rust's chain and its range test (`_km_second_order_fast`) | — | 23 ns | 33 ns |
| `kinetics/michaelis_menten` | (100, 5, 10): three guards, Rust's chain and its range test (`_km_michaelis_menten_fast`) | 2.92 ns | 17 ns | 24 ns |
| `thermo/heat_transfer_tie` | 3·y·1 with 3y = 2⁵³ + 1, a midpoint: `_km_mul3_fast`'s test fails, `_km_cr_mono` resolves it (`_km_mono_tie`) | — | 9.6 ns | 365 ns |
| `phase/osmotic_pressure_subnormal` | i = 2⁻¹⁰⁷⁴: a subnormal result, the exact path (five factors' integer product, `_km_cr_ratio`) | — | 26 ns | 4.04 µs |
| `kinetics/second_order_concentration_edge` | (1, DBL_MAX, 10): k·t overflows, the unbounded-exponent slow path | — | 24 ns | 230 ns |
| `kinetics/michaelis_menten_edge` | (100, DBL_MAX, DBL_MAX): ∞/∞ in Rust, the slow path | — | 22 ns | 126 ns |
| `stoichiometry/percent_yield` | 8 of 10: `_km_cr_mul_div`, now `src/mathx.cyr`'s (M4f-2's fast path, both signs) | — | 81 ns | 82 ns |
| `stoichiometry/theoretical_yield_moles` | (3, 2) over (2, 1), × 2 | — | 115 ns | 115 ns |
| `stoichiometry/percent_composition_glucose` | three RN(mass·count·100/mw) | — | 278 ns | 266 ns |
| `stoichiometry/theoretical_yield_moles_tie` | y·3/2 with 3y = 2⁵³ + 1: `_km_mul_div_tie` | — | 136 ns | 132 ns |
| `stoichiometry/percent_yield_subnormal` | 100·3·2⁻¹⁰⁷⁴/200: the exact path (now `_km_cr_mono_acc`, odd significands) | — | 2.19 µs | 1.62 µs |
| `stoichiometry/limiting_reagent` | (unchanged code) | — | 34 ns | 34 ns |

The rounded-once functions cost 65–139 ns where Rust's two to four operations took 6–27: the decision's "roughly 5×
slower (e.g. ~15 → ~80 ns)", more for the unvalidated products whose Rust body is two multiplies, and for osmotic
pressure's five factors. **One design measured and not kept**: every shape through the general `_km_cr_mono` (its factor
tables and loops) cost `heat_transfer` 148 ns, `specific_heat_from_calorimetry` 110, the colligative pair 107–108,
`osmotic_pressure` 200, `steady_state_intermediate` 115, `pre_equilibrium_rate` 138, and moved stoichiometry's
`percent_yield` to 108, `theoretical_yield_moles` 141, `percent_composition_glucose` 422 and the tie to 262 ns; the
straight-line fast paths per shape (`_km_mul_div_fast`, `_km_mul3_fast`, `_km_div2_fast`, `_km_mul3_div_fast`,
`_km_mul5_fast`) bring them to the numbers above, with stoichiometry back at M4f-2's. Reading a double's bits with
`load64(&x)` instead of a call to `_km_bits` saved 4–15 ns a call; inlining `_km_dd_mul_d`'s TwoProduct about 1. The
#20-style functions pay 7–10 ns for their range test, a function call (`_km_*_fast`, shared with the tests and
`cr_verify.py`).

## The overflow audit's B2 (2026-11-10)

ADR 0002 #29 and #30: solution's and electrochemistry's quotients and spectroscopy's products rounded once on B1's
`_km_cr_mono` (four new straight-line shapes), `vibrational_frequency` on the new correctly rounded √ of a quotient
(`_km_cr_sqrt_div`), `nernst_potential` with an unbounded exponent behind a range test (Rust's chain in range). Before:
B2's new benches added to `tests/kimiya.bcyr` first and run on the code before B2; after: the same trimmed copy of the
bench file (only these benches; it reproduces the full run's numbers) on the code after, the same host, back to back.
Two benches have a Rust counterpart; the rest are Cyrius-only.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `spectroscopy/absorbance` | (100, 1, 0.01): RN(ε·l·c) by `_km_mul3_fast` | 0.98 ns | 5 ns | 68 ns |
| `electrochemistry/nernst_potential` | (0.342, 2, 298.15 K, 0.1): four guards, `_km_cr_ln`, Rust's chain and its range test (`_km_nernst_fast`) | 6.42 ns | 83 ns | 92 ns |
| `solution/dilution` | (2 M, 0.1 L, 0.5 L): three guards, RN(M₁·V₁/V₂) by `_km_mul_div_fast` | — | 17 ns | 86 ns |
| `solution/solution_conductivity` | (126.4, 0.01 M): RN(Λ·c/1000) by `_km_mul_div_fast` | — | 6 ns | 72 ns |
| `electrochemistry/faraday_mass_deposited` | (9648.5 C, 63.546, 2): RN(Q·M/(n·F)) by `_km_mul2_div2_fast` | — | 19 ns | 109 ns |
| `electrochemistry/faraday_charge_required` | (1 g, 63.546, 2): RN(m·n·F/M) by `_km_mul3_div_fast` | — | 21 ns | 122 ns |
| `electrochemistry/faraday_mass_ratio` | (1 g, 63.546, 2, 107.868, 1): RN(m_A·M_B·n_A/(n_B·M_A)) by `_km_mul3_div2_fast` | — | 33 ns | 148 ns |
| `spectroscopy/concentration_from_absorbance` | (0.5, 100, 1 cm): RN(A/(ε·l)) by `_km_div2_fast` | — | 16 ns | 90 ns |
| `spectroscopy/vibrational_frequency` | HCl (516 N/m, 1.6266·10⁻²⁷ kg): RN(√(k/μ)) by `_km_sqrt_div_fast` | — | 18 ns | 87 ns |
| `spectroscopy/rigid_rotor_energy` | (J = 2, 3.19·10¹¹ Hz): RN(h·B·J·(J + 1)) by `_km_mul4_fast` | — | 14 ns | 114 ns |
| `spectroscopy/rotational_constant` | (2.64·10⁻⁴⁷ kg·m²): RN(h/(8π²·I)) by `_km_div3_fast` | — | 18 ns | 117 ns |
| `spectroscopy/diatomic_moment_of_inertia` | (1.6266·10⁻²⁷ kg, 1.2746·10⁻¹⁰ m): RN(μ·r·r) by `_km_mul3_fast` | — | 11 ns | 74 ns |
| `electrochemistry/nernst_potential_edge` | (0.342, 2, DBL_MAX, 0.1): R·T overflows, the unbounded-exponent slow path | — | 90 ns | 202 ns |
| `spectroscopy/vibrational_frequency_edge` | (500, 2⁻¹⁰²²): k/μ overflows in Rust; the √ fast path | — | 19 ns | 88 ns |
| `spectroscopy/vibrational_frequency_subnormal` | (10⁻³⁰⁰, 2⁻¹⁰⁷⁴): a subnormal μ, the exact root (multi-word comparisons) | — | 26 ns | 883 ns |
| `electrochemistry/faraday_mass_ratio_tie` | y·3 with 3y = 2⁵³ + 1, a midpoint: the shape's test fails, `_km_mono_tie` resolves it | — | 33 ns | 402 ns |

The rounded-once functions cost 68–148 ns where Rust's two to four operations took 5–33 ns, as B1's did (the decision's
"roughly 5× slower"; more for the unvalidated products whose Rust body is two multiplies, and for the five-factor
`faraday_mass_ratio`). The √ costs 87 ns: a division, a square root, two TwoProducts and the rounding test; inlining the
two TwoProducts was measured (87 → 86 ns) and not kept. `nernst_potential` pays 9 ns for its range test; its 25 °C form
takes the same path. The exact root (`vibrational_frequency_subnormal`) is a few multi-word products and comparisons
after a guess from the fast path's arithmetic, under a microsecond; it serves subnormal operands, possibly subnormal
results and the rare failed test.

## The overflow audit's B3 (2026-11-11)

ADR 0002 #31 and #32: gas's ideal gas law and density and potential's Coulomb pair rounded once on B1's `_km_cr_mono`
(B1's and B2's straight-line shapes), the Morse width on the new correctly rounded w·√((x/y)·2^e)
(`_km_cr_mul_sqrt_ratio`), the three equations of state, Lennard-Jones, Born–Landé and the binding energy per nucleon with
an unbounded exponent behind range tests (Rust's chain in range). Before: B3's new benches added to `tests/kimiya.bcyr`
first and run, as a trimmed copy of the bench file (only these benches), on the code before B3 (an export of the
committed tree); after: the same copy on the code after, the same host, back to back. Three benches have a Rust
counterpart (`lennard_jones` is M4d's Cyrius-only one); the rest are Cyrius-only.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `gas/ideal_gas_pressure` | (1, 300 K, 0.025 m³): three guards, RN(n·R·T/V) by `_km_mul3_div_fast` | 2.18 ns | 18 ns | 119 ns |
| `nuclear/binding_energy_per_nucleon_fe56` | Fe-56: the u32 checks, Δm, the inline test (Δm·931.494… finite), one division | 2.02 ns | 28 ns | 29 ns |
| `potential/lennard_jones` | Ar at r ≈ r_min: Rust's chain and its range test (`_km_lj_fast`) | — | 36 ns | 39 ns |
| `gas/ideal_gas_volume` | (1, 273.15 K, 101325 Pa): RN(n·R·T/P) by `_km_mul3_div_fast` | — | 20 ns | 119 ns |
| `gas/gas_density` | air at STP: RN(P·M/(R·T)) by `_km_mul2_div2_fast` | — | 20 ns | 108 ns |
| `gas/van_der_waals_pressure` | CO₂ in 1 L at 300 K: Rust's chain, the region test (inside it no per-operation test) | — | 35 ns | 55 ns |
| `gas/peng_robinson_pressure` | N₂ at 300 K, 0.025 m³/mol: `pr_coefficients` and the pressure inline, `_km_pr_box` | — | 78 ns | 100 ns |
| `gas/redlich_kwong_pressure` | N₂ at 300 K, 0.025 m³/mol: `_km_cr_pow(Tc, 2.5)` and the chain, the region test | — | 182 ns | 187 ns |
| `potential/coulomb` | (e, −e, 1 Å): RN(k·q₁·q₂/r) by `_km_mul3_div_fast` | — | 16 ns | 126 ns |
| `potential/coulomb_force` | (e, −e, 1 Å): RN(k·q₁·q₂/(r·r)) by `_km_mul3_div2_fast` | — | 18 ns | 145 ns |
| `potential/morse_width_parameter` | HCl: RN(ω·√(μ/(2·D_e))) by `_km_mul_sqrt_ratio_fast` | — | 33 ns | 135 ns |
| `inorganic/born_lande_lattice_energy` | NaCl: Rust's chain, the region test | — | 46 ns | 50 ns |
| `gas/van_der_waals_pressure_edge` | n = 10³⁰⁶: n·R·T overflows, the unbounded-exponent slow path | — | 36 ns | 260 ns |
| `gas/peng_robinson_pressure_edge` | Tc = 2⁻¹⁰²²: T/Tc overflows (Rust's NaN), the slow path | — | 115 ns | 618 ns |
| `gas/redlich_kwong_pressure_edge` | Tc = 10¹³⁰: Tc^2.5 overflows, the exact rescale on the slow path | — | 134 ns | 584 ns |
| `potential/lennard_jones_edge` | ε = 10⁻³⁰⁰, σ/r = 5.8·10²⁵: (σ/r)¹² overflows, the slow path | — | 46 ns | 208 ns |
| `potential/morse_width_parameter_subnormal` | (1, 2⁻¹⁰⁷⁴, 2¹⁰²³): a subnormal μ, the exact root | — | 42 ns | 1.27 µs |
| `inorganic/born_lande_lattice_energy_edge` | (DBL_MAX, 1, 1, DBL_MAX, DBL_MAX): the prefix overflows, the slow path | — | 66 ns | 287 ns |
| `nuclear/binding_energy_per_nucleon_edge` | (1000, 1000, DBL_MAX): Δm·931.494… overflows, the slow path | — | 39 ns | 117 ns |
| `potential/morse_width_parameter_tie` | ω = y with 3y = 2⁵³ + 1, μ = 9, D_e = ½: a midpoint, the exact root | — | 35 ns | 1.27 µs |

The rounded-once functions cost 108–145 ns where Rust's three or four operations took 16–20 ns, as B1's and B2's did; the
Morse width's 135 ns is B2's √ (87 ns) plus the product by ω (a TwoProduct, FastTwoSum and the rounding test). The
#20-style functions pay 3–22 ns on physical inputs: their range tests first ask whether the inputs lie in a region where
every intermediate is provably normal (Tc, Pc, T, V within 2^±60; `_km_pr_box` and its inline siblings), and only outside
it test each operation with `_km_ux_normal`, which took Peng–Robinson from 147 to 100 ns, van der Waals from 73 to 55,
Redlich–Kwong from 217 to 187 and Born–Landé from 67 to 50; inlining `_km_ux_normal`'s bit test instead was measured
(147 → 156 ns) and not kept. The slow paths are the unbounded-exponent chains (117–618 ns), and the Morse width's exact
root (1.27 µs: a few multi-word products and comparisons) serves subnormal operands, possibly subnormal results, its
midpoints and the rare failed test.

In `tests/kimiya.bcyr` B3's seventeen Cyrius-only benches run 100,000 × 20 calls each (the file's size for its slow
rows), not 1,000,000 × 20: the file took about 255 s before B3, and B3's rows at the larger size (71 s) took it past
`cyrius bench`'s default deadline of 300 s per bench file, so it was killed. With them `cyrius bench` takes 279 s on this
host, and in that full run the rows read 5–10 % above the back-to-back figures above (`ideal_gas_volume` 123 ns,
`peng_robinson_pressure` 103, `redlich_kwong_pressure` 193, `morse_width_parameter` 140, the edges 122–644 ns, the exact
root 1.29–1.32 µs). The next batch's benches need the deadline raised (`[test] timeout` in `cyrius.cyml`) or a second
bench file. (B4 split the file: below.)

## The overflow audit's B4 (2026-11-11)

ADR 0002 #33 and #34: the sums with cancellation and the Shomate polynomials with an unbounded exponent behind fast
tests (Rust's chain in range), and #29's corrected fast-path bound. **The bench split** (the user's decision): the
audit's Cyrius-only benches, B1's to B4's (`bench_audit_b1` … `bench_audit_b4`), live in `tests/audit.bcyr`;
`tests/kimiya.bcyr` keeps the benches that mirror `rust-old/benches/benchmarks.rs` and the Cyrius-only ones from before
the audit. The names are unchanged, and B3's seventeen are back at the standard 1,000,000 × 20 calls (its two exact
roots at 100,000 × 20, as the other microsecond rows). On this host `cyrius bench tests/kimiya.bcyr` takes 177 s (261 s
before the split; 279 s in B3's own run) and `cyrius bench tests/audit.bcyr` 189 s, each well inside `cyrius bench`'s
300 s per-file deadline; moving only the edge, slow-path and tie benches would have left kimiya.bcyr at about 214 s, or
232 s with B3's typical benches restored, so the audit's typical-input benches moved too. B3's benches at the restored
counts read as before (`ideal_gas_volume` 115 ns, `peng_robinson_pressure` 95, `redlich_kwong_pressure` 190,
`morse_width_parameter` 128, the edges 114–601 ns, the exact roots 1.22–1.26 µs).

Before: `tests/kimiya.bcyr` on the code before B4 (the committed tree), and B4's new benches, as a trimmed copy of
`tests/audit.bcyr` (only these benches), run on an export of the committed tree; after: both files on the code after,
the same host.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `reaction/gibbs_free_energy` | (−100000, 298, −200): three guards, ΔH − T·ΔS and its range test (T·ΔS normal, inline) | 1.75 ns | 17 ns | 18 ns |
| `thermochem/reaction_enthalpy_ch4` | CH₄ + 2O₂ → CO₂ + 2H₂O(l): four table scans, two sums with the per-product test inline, × 1000 | 38.58 ns | 259 ns | 293 ns |
| `thermochem/enthalpy_change_cp` | the Rust bench's CO₂(g) call: the table-row path, T inside the box (one inline test) | 24.00 ns | 171 ns | 171 ns |
| `thermochem/shomate_cp` | `s.cp(500)`, CO₂'s first segment: Rust's chain and the box test (`_km_shomate_box`) | 2.34 ns | 15 ns | 32 ns |
| `thermochem/adiabatic_flame_ch4_air` | the flame at φ = 1 in air: about 90 table-row evaluations, each T inside the box | 301.48 ns | 3.40 µs | 3.29 µs |
| `electrochemistry/nernst_potential` | (0.342, 2, 298.15 K, 0.1): #29's tests at 2⁻¹⁰²¹ (inline) | 6.42 ns | 90 ns | 89 ns |
| `thermo/hess_law` | (−100, 50, −30): Rust's sum, "finite" | — | 13 ns | 15 ns |
| `thermo/enthalpy_from_formation` | ([−400], [−200, −100]): the two sums, both finite | — | 20 ns | 24 ns |
| `thermochem/reaction_gibbs_energy` | as `reaction_enthalpy_ch4`, ΔG_f° | — | 259 ns | 299 ns |
| `thermochem/reaction_entropy` | as `reaction_enthalpy_ch4`, S° | — | 259 ns | 286 ns |
| `nuclear/q_value` | D + T → ⁴He + n: the sums, the difference finite, × 931.494… | — | 24 ns | 26 ns |
| `solution/ionic_strength` | 0.1 M of a 1:1 salt: Rust's ½·Σ, the sum finite | — | 24 ns | 31 ns |
| `inorganic/born_haber_lattice_energy` | NaCl's cycle: the result finite, ΔH_diss/2 exact | — | 10 ns | 14 ns |
| `solution/kohlrausch_molar_conductivity` | (126.4, 89.9, 0.01 M): `_km_kohlrausch_fast` | — | 22 ns | 28 ns |
| `solution/limiting_molar_conductivity` | NaCl: `_km_limiting_fast` (two products normal) | — | 10 ns | 24 ns |
| `kinetics/zero_order_concentration` | (1, 0.01, 50): k·t normal, the clamp | — | 12 ns | 17 ns |
| `thermochem/shomate_enthalpy_antiderivative` | CO₂(g)'s segment at 500 K: Rust's chain and the box test | — | 37 ns | 54 ns |
| `thermo/hess_law_edge` | [DBL_MAX, DBL_MAX, −DBL_MAX]: Rust's +∞, the unbounded sum | — | 22 ns | 162 ns |
| `thermo/enthalpy_from_formation_edge` | [DBL_MAX, DBL_MAX] twice: Rust's NaN, the slow path | — | 31 ns | 212 ns |
| `thermochem/reaction_enthalpy_edge` | [Al₂O₃ ×10³⁰⁸] → [Al₂O₃ ×10³⁰⁸]: Rust's NaN, the slow path (the lookups twice) | — | 225 ns | 587 ns |
| `nuclear/q_value_edge` | [M, M] − [M, 0.999M]: Rust's NaN, the slow path | — | 33 ns | 250 ns |
| `solution/ionic_strength_edge` | [(0.4·DBL_MAX, 2)]: Rust's +∞, the slow path | — | 26 ns | 144 ns |
| `inorganic/born_haber_lattice_energy_edge` | (M, −M, M, 0, 0): Rust's +∞, the slow path | — | 10 ns | 175 ns |
| `reaction/gibbs_free_energy_edge` | (M, 2, 0.75M): Rust's −∞, the slow path | — | 27 ns | 117 ns |
| `solution/kohlrausch_molar_conductivity_edge` | (M, M, 1.21): Rust's −∞, the slow path | — | 23 ns | 129 ns |
| `solution/limiting_molar_conductivity_edge` | (M, 1.5, −M, 1): Rust's +∞, the slow path | — | 11 ns | 132 ns |
| `kinetics/zero_order_concentration_edge` | (−M, −M, 1.5): Rust's +∞, the slow path | — | 13 ns | 138 ns |
| `thermochem/shomate_cp_edge` | H₂O(g)'s segment at 10³⁰⁰ K: Rust's ∞ − ∞ = NaN, the per-operation tests, the slow path | — | 29 ns | 311 ns |
| `thermochem/shomate_enthalpy_antiderivative_edge` | A = 20.786 alone at 5·10⁻³²⁴ K: Rust's NaN, the slow path | — | 47 ns | 474 ns |

On physical inputs the fast tests cost 1–7 ns for the scalar functions (`limiting_molar_conductivity` 14 ns: its fast
path is a function, shared with the tests and `cr_verify.py`, testing two products), 2–7 ns a slice for the pure sums
("finite"), and about 9 ns per record for the reaction sums (the per-product test and the coefficient's finiteness,
inline; the table scans dominate). The Shomate box costs 17 ns per `KmShomate.cp` call (the test of T and five
coefficients, a function call); inlining it was measured (34 and 39 ns, with bit and float compares) and not kept. The
table-row path, whose coefficients are all inside the box, tests T alone, inline: `enthalpy_change_cp` is unchanged and
the flame 3% faster than before (the general evaluation was a call, the row path now inline). The slow paths are the
unbounded-exponent chains, 117–587 ns (the Shomate polynomials' per-operation tests before them: 311 and 474 ns).
`nernst_potential`'s bound moved from 2⁻¹⁰²² to 2⁻¹⁰²¹ with the same inline compares (89 ns; `_km_ux_normal` calls
measured 99 ns).

## The overflow audit's B5 (2026-11-11)

ADR 0002 #35: the exp / exp10 *arguments* of `equilibrium_constant` and the three activity models with an unbounded
exponent behind inline fast tests (Rust's chain in range), the argument's final operation rounded once, then the
unchanged correctly rounded exp / 10^x. None of the four had a bench; `bench_audit_b5` in the new `tests/audit2.bcyr`
adds four typical and six edge benches (Cyrius-only: `rust-old/benches` has neither function family).

Before: B5's benches as a copy of `tests/audit.bcyr`'s harness (only these benches) on the code before B5 (the committed
tree); after: the same copy on the code after, two runs each, and `tests/audit2.bcyr`, the same host.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `reaction/equilibrium_constant` | (−5000 J/mol, 298 K): two guards, −ΔG/(R·T) and its range test (R·T normal, inline), `_km_cr_exp` | — | 64 ns | 64 ns |
| `solution/debye_huckel_limiting` | (2, 0.01 M, 0.509): the checks, (((−A)·z)·z)·√I and its test (A·z² finite), `_km_cr_exp10` | — | 86 ns | 85 ns |
| `solution/debye_huckel_extended` | (1, 0.05 M, 0.509, 0.328, 3 Å): the checks, log₁₀ γ and its test (the numerator and (B·å)·√I finite) | — | 105 ns | 100 ns |
| `solution/davies_activity` | (1, 0.1 M, 0.509): the checks, log₁₀ γ and its test (A·z² finite) | — | 98 ns | 96 ns |
| `reaction/equilibrium_constant_edge` | (−DBL_MAX, DBL_MAX): R·T past DBL_MAX, Rust's e^−0 = 1; the slow path, e^0.12 | — | 26 ns | 151 ns |
| `reaction/equilibrium_constant_subnormal` | (−10⁻³²⁰, 10⁻³²²): R·T subnormal; the slow path | — | 80 ns | 166 ns |
| `solution/debye_huckel_limiting_edge` | (100000, 0, 10³⁰⁰): A·z² past DBL_MAX at I = 0, Rust's NaN; the slow path, γ = 1 | — | 30 ns | 93 ns |
| `solution/debye_huckel_extended_edge` | (1, 1, 10³⁰⁸, 10¹⁵⁵, 10¹⁵⁵): B·å past DBL_MAX, Rust's 1; the slow path, 0.977 | — | 45 ns | 277 ns |
| `solution/debye_huckel_extended_overflow` | (100000, 1, 10³⁰⁰, 10¹⁵⁴, 10¹⁵⁴): A·z² past DBL_MAX, Rust's 0; the slow path, 10⁻¹⁰⁰ | — | 47 ns | 276 ns |
| `solution/davies_activity_edge` | (100000, 0, 10³⁰⁰): A·z² past DBL_MAX at I = 0, Rust's NaN; the slow path, γ = 1 | — | 39 ns | 162 ns |

The fast tests are compares on Rust's own intermediates (R·T against two bounds; an intermediate minus itself against
0), so the typical inputs keep their cost, within a few ns either way (the extended law's 105 → 100 ns was the same in
both runs; not investigated further). The slow paths are the unbounded-exponent arguments, 93–277 ns (the extended law's
eight unbounded operations the dearest), then the same exp or 10^x; Rust's own edge results were cheap (exp of ∓0,
10^NaN).

**The third bench file.** Added to `tests/audit.bcyr`, the ten benches took it to 219 s on this host (189 s before),
73 % of `cyrius bench`'s 300 s per-file deadline with the audit's B8–B12 still to add theirs, and a loaded host
stretches every bench; so B5's go into a new `tests/audit2.bcyr` (31 s), the audit's later batches with them, the same
way B4 split `tests/kimiya.bcyr` (the user's decision then: split rather than raise the deadline). `tests/audit.bcyr` is
unchanged (189 s); `cyrius bench` runs all three files.

## The overflow audit's B6 (2026-11-11)

ADR 0002 #36: `fft_spectrum` under #20's semantics over the whole transform, behind a per-operation operand test on Rust's
own transform (and a bin's unbounded √ where only its square sum leaves the range); #37: `find_peaks`' threshold test
exact, Rust's rounded threshold deciding behind a proven bound; #38: `adiabatic_flame_temperature`'s bracket and its
bisection compare signs exactly. The typical benches are `tests/kimiya.bcyr`'s (Cyrius-only for the FFT and the peaks:
`rust-old/benches` has neither); `bench_audit_b6` in `tests/audit2.bcyr` adds five edge benches.

Before: a copy of `tests/kimiya.bcyr`'s harness with only these benches, and of `tests/audit2.bcyr`, on the committed
tree (`git archive HEAD`); after: the same on the batch, alternating runs, the same host.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `spectroscopy/fft_spectrum_256` | a seeded signal in [0, 1), in place: the transform with its operand checks, the bins' range test | — | 27.3 µs | 30.5 µs |
| `spectroscopy/fft_spectrum_1024` | the same, n = 1024 | — | 129.2 µs | 143.1 µs |
| `spectroscopy/fft_spectrum_4096` | the same, n = 4096 | — | 605.6 µs | 669.4 µs |
| `spectroscopy/find_peaks_1000` | 1,000 values, t = 0.1: Rust's threshold and the bound once, a subtraction and two compares a candidate | — | 19.0 µs | 20.3 µs |
| `thermochem/adiabatic_flame_ch4_air` | CH₄ in air: two sign compares for the bracket's multiply | 301.48 ns | 3.37 µs | 3.40 µs |
| `spectroscopy/fft_spectrum_256_huge` | the signal × 10³⁰⁷: sums past DBL_MAX, Rust's ∞; the unbounded transform | — | 27.7 µs | 318.1 µs |
| `spectroscopy/fft_spectrum_256_subnormal` | the signal × 10⁻³¹⁰: data below 2⁻⁹²¹, Rust's degraded; the unbounded transform | — | 35.4 µs | 349.8 µs |
| `spectroscopy/fft_spectrum_256_squares` | the signal × 10²⁰⁰: the transform in range, every square sum past DBL_MAX; the per-bin unbounded √ | — | 27.4 µs | 45.1 µs |
| `spectroscopy/find_peaks_1000_subnormal` | k·2⁻¹⁰⁷⁴ with k < 16, t = ½: about half the candidates within the bound, the exact predicate | — | 19.5 µs | 39.0 µs |
| `thermochem/adiabatic_flame_tiny` | CH₄/air × 10⁻²⁰⁰: Rust's bracket product −0, its Err; here bracketed, Newton's guess | — | 0.96 µs | 1.48 µs |

**The FFT's fast-path cost**, +11–12 %, inside the ~15 % aimed at: measured without the per-butterfly operand checks the
three sizes take 29.1, 139.6 and 638.5 µs, so the checks (two absolute values and compares a butterfly, two more a
butterfly of a stage's first block for the twiddles) are about 6 %, and the rest is the restructured magnitude pass (the
range test of each bin's square sum, which also proves the retained outputs finite) and the transform's extra state.
Tried: skipping the checks at k = 0 (w = (1, 0), whose products are exact) gained under 1 % and was not kept; merging
the retained outputs' finiteness pass into the magnitude pass gained 1 % and was kept. The slow path is about 11× Rust's
chain: each of the n/2·log₂ n butterflies and twiddle steps is sixteen unbounded operations (`_km_ux_*`). A bin's
unbounded √ alone costs about 0.14 µs (the `_squares` row: 129 bins).

`tests/audit2.bcyr` takes 45 s with B6's benches (31 s before); `tests/kimiya.bcyr` and `tests/audit.bcyr` are unchanged
in content, the former's three FFT benches a few seconds longer in all.

## The overflow audit's B7 (2026-11-12)

ADR 0002 #39: `balance_equation` runs Rust's i64 elimination first, unchanged, and at its first overflow the same steps
from the original matrix on signed multi-word integers (31-bit limbs, Knuth's division, Lehmer's gcd), so where Rust's
checked arithmetic overflows the port returns the exact algorithm's result. The typical benches are
`tests/kimiya.bcyr`'s (`balance_ch4_combustion` is Rust's); `bench_audit_b7` in `tests/audit2.bcyr` adds three edge
benches, where Rust's own run stops at its overflow Err.

Before: a copy of the benches' harness with only these benches on the committed tree (`git archive HEAD`); after: the
same on the batch, alternating runs, the same host.

| Bench | Body | Rust | Before | After |
|---|---|---|---|---|
| `stoichiometry/balance_ch4_combustion` | CH₄ + 2O₂ → CO₂ + 2H₂O (3 × 4): the i64 fast path, no overflow | 233.90 ns | 805 ns | 778 ns |
| `stoichiometry/balance_kmno4_hcl` | 2KMnO₄ + 16HCl → … (5 × 6): the same | — | 2.24 µs | 2.18 µs |
| `stoichiometry/balance_survey_2x3` | the survey's example: Rust's first product 3037000500² overflows, its Err; here the multi-word path, [1, 1, 1] | — | 0.12 µs | 3.6 µs |
| `stoichiometry/balance_overflow_10x9` | a 10 × 9 golden row whose elimination overflows i64 (entries ≤ 2,136): Rust's Err; here the exact x | — | 6.6 µs | 57 µs |
| `stoichiometry/balance_dense_30x31` | 30 × 31 seeded entries below 2⁴⁵, the last column solved for x[j] = 1 + j mod 3: Rust's Err at the first product; here products past 3,000 bits, 423 KB of tables | — | 5.4 µs | 158 ms |

**The fast path** is the old code in its own function (each overflow check returns a restart code instead of the Err),
so in-range matrices keep their speed: the 3 % gain is the call's layout, not new work. **The slow path's cost** is the
elimination's own: on the 30 × 31 matrix, 11.0 million limb products in its 870 row updates (pivot·a and factor·a for
31 entries each), 10.7 million limb steps in the 26,800 exact divisions of reduce_row, and Lehmer's gcd (22,000 of its
unimodular steps, 13,000 full division steps, 1.8 million limb steps of its combinations). Plain Euclid in place of
Lehmer's gcd took a 39 × 40 dense matrix 1.8 s (about 844,000 of its division steps; 0.7 s with Lehmer's), and
stepping pointers in the multiply and divide loops instead of indexing gained a further 13 % (193 → 168 ms on the
30 × 31, before the final runs above).

`tests/audit2.bcyr` takes 49 s with B7's benches (45 s before); `tests/kimiya.bcyr` and `tests/audit.bcyr` are
unchanged.

## The overflow audit's B10a (2026-11-12)

The first part of the survey's B10, for ADR 0002 #40: `src/mathx.cyr`'s `_km_cr_exp_ux`, `_km_cr_expm1_ux` and
`_km_cr_exp10_ux`, e^x, e^x − 1 and 10^x of an unbounded-exponent pair, correctly rounded to 53 bits at any exponent.
No function calls them yet (#40's functions will, on their slow paths only), so no other bench moves. Cyrius-only (Rust
has no such function); `bench_audit_b10a` in `tests/audit2.bcyr` times each tier beside the binary64 function at the same
argument, in the same harness and run (the binary64 rows are kimiya.bcyr's `mathx/cr_*` arguments).

| Bench | Argument and path | Binary64 function | Unbounded |
|---|---|---|---|
| `mathx/cr_exp` → `cr_exp_ux` | −50000/(8.314·300): the first tier, `_km_cr_exp` itself and its pair | 46 ns | 54 ns |
| `mathx/cr_exp_ux_acc` | the same, the accurate path alone (the binary64 accurate core, `_km_rn3`) | — | 1.86 µs |
| `mathx/cr_exp_ux_subnormal` | −720.5 (binary64's result subnormal): the second tier, exp's fast core unscaled | — | 84 ns |
| `mathx/cr_exp_ux_minus_1e4` | −10⁴: the double-double reduction, then the fast core's tail | — | 159 ns |
| `mathx/cr_exp_ux_1e12` | 10¹²: the same | — | 158 ns |
| `mathx/cr_exp_ux_1e15` | 10¹⁵ (past 2⁴¹): the integer reduction on the fast path | — | 614 ns |
| `mathx/cr_exp_ux_acc_minus_1e4` | −10⁴, the accurate path alone: the integer reduction and the accurate core's tail | — | 2.14 µs |
| `mathx/cr_expm1` → `cr_expm1_ux` | ln 2: the first tier, `_km_cr_expm1` itself | 54 ns | 65 ns |
| `mathx/cr_expm1_ux_1e4` | 10⁴: e^x − 1 rounds as e^x, the double-double reduction | — | 165 ns |
| `mathx/cr_exp10` → `cr_exp10_ux` | −4.75: the first tier, `_km_cr_exp10` itself | 60 ns | 69 ns |
| `mathx/cr_exp10_ux_1e4` | 10⁴: the double-double reduction with 512·log₂ 10 | — | 159 ns |
| `mathx/cr_exp10_ux_acc_minus_1e12` | −10¹², the accurate path alone | — | 2.17 µs |

**The first tier** costs 8–11 ns over the binary64 function: the range test, the scaling and the pair, inline in the
entry point (79 ns before the inlining, through the shared dispatcher). **The reductions**: the integer one (eleven
columns of 32-bit products, the fraction's six terms summed in triple-double, a triple-double product) is about 0.55 µs,
so the fast path reduces in double-double wherever 1.5·2⁵²'s rounding to an integer can hold K (|x| < 2⁴¹; 2⁴⁰ for
10^x): two TwoProducts, two TwoSums and one more for f·ln2/512, 158 ns in all, and the integer reduction serves only
larger arguments and the accurate path, which runs about once in 50,000 calls past binary64's range. The realistic slow
path of #40's single-exponential chains (|x| below about 1,500, where a binary64 prefactor can bring the product back)
is the first, second or double-double tier.

`tests/audit2.bcyr` takes 66 s with B10a's benches (49 s before); `tests/kimiya.bcyr` and `tests/audit.bcyr` are
unchanged.

## M5a numerics, reaction_dynamics (2026-11-12)

`tests/numerics.bcyr`, a new bench file (39 s; `tests/kimiya.bcyr` and the audit files are unchanged). The two benches
with a Rust counterpart are `rust-old/benches/benchmarks.rs`'s bodies, Rust's criterion medians from
`rust-old/benchmarks.md`; the five others are Cyrius-only benches of the Rust tests' systems, compared with an ad hoc
Rust timing of the same bodies (`std::time::Instant` around a `black_box` loop, the median of 21 batches, release build
against `rust-old/`, the same host and idle; it reads 4.75 and 5.72 µs for the two criterion benches, 8–13 % above
criterion's medians, so its column is a rough guide). Rust returns a fresh `Vec` of trajectory points (one heap vector
per point); the port fills the caller's buffer.

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `dynamics/first_order_100steps` | `simulate_first_order(1, 0.1, 10, 100)`: 100 RK4 steps (4 callbacks each, Neumaier compensation), 101 points | 4.41 µs | 11.2 µs | 2.5× |
| `dynamics/consecutive_100steps` | `simulate_consecutive(1, 0.5, 0.1, 20, 100)`: the same over 3 species | 5.04 µs | 21.4 µs | 4.3× |
| `dynamics/adaptive_first_order` | `simulate_adaptive` of dy = −0.1·y to t = 10, tol 10⁻⁸ (the Rust test): Dormand–Prince, 7 callbacks and one correctly rounded pow per step | 1.80 µs (ad hoc) | 6.19 µs | 3.4× |
| `dynamics/stiff_robertson_400steps` | `simulate_stiff` of Robertson to t = 40 in 400 steps (the Rust test): Newton per step, a 3 × 4 elimination per iteration | 121 µs (ad hoc) | 498 µs | 4.1× |
| `dynamics/stochastic_100steps` | `simulate_stochastic` of dy = −y, σ = 0.1, 100 steps, seed 42: two callbacks and one Box–Muller draw per step | 4.85 µs (ad hoc) | 24.8 µs | 5.1× |
| `numerics/gaussian_elimination_8x8` | `_km_gaussian_elimination` of a fixed 8 × 9 augmented matrix, copied in each call (Rust: `gaussian_elimination` on a clone) | 353 ns (ad hoc) | 2.48 µs | 7.0× |
| `numerics/pcg32_next_normal` | one Box–Muller draw: two PCG32 outputs, ln, √, cos | 22.6 ns (ad hoc) | 200 ns | 8.8× |

Every body is bit-identical to Rust's on the golden rows (`reaction_dynamics` 954, `numerics` 243, 0 ulp; the four rows
where glibc misrounds a pow, ln or cos are the corrected oracle's).

**One optimization, measured.** The elimination addressed each entry as m + 8·(row·(n + 1) + j) inside its inner loops;
it now keeps each row's base pointer and steps by 8 bytes (the same operations, the same bits): the two versions side by
side in one bench run, `gaussian_elimination_8x8` 3.15 → 2.56 µs (−19 %). `simulate_stiff` solves one such system per
Newton iteration.

**Reading the gap.** The RK4 and Euler–Maruyama loops are the port's typical 2.5–5× (each callback is a `callptr` with
five arguments, each component's update a few loads and stores through the work tables, where Rust inlines the closure);
the normal draw is mostly the correctly rounded ln and cos (`mathx/cr_ln` 49 ns, `mathx/cr_cos` 82 ns);
dopri45's correctly rounded pow is one call per step. Parity, not speed, is the bar.

## M5b network, fitting, sensitivity (2026-11-12)

Seven benches appended to `tests/numerics.bcyr` (now 32 s with them). `fitting/arrhenius_6pts` and
`fitting/polynomial_deg2_20pts` are `rust-old/benches/benchmarks.rs`'s bodies, Rust's criterion medians from
`rust-old/benchmarks.md`; the five others are Cyrius-only, compared with M5a's ad hoc Rust timing of the same bodies
(the median of 20 batches; it reads 256 and 537 ns for the two criterion benches, 7–10 % above criterion's medians).
Rust returns fresh `Vec`s (the matrix rows, the coefficients, the gradient, the tape); the port fills the caller's
buffers.

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `fitting/arrhenius_6pts` | `fit_arrhenius` of 6 points, T = 300 … 400 K, k = 10¹⁰·exp(−50 000/(R·T)): 6 correctly rounded ln, the degree-1 QR, one exp | 233.40 ns | 1.01 µs | 4.3× |
| `fitting/polynomial_deg2_20pts` | `fit_polynomial` of y = 1 + 2x + 0.5x² at 20 points, degree 2: 60 powi, the 3-column modified Gram–Schmidt, back substitution | 500.05 ns | 2.32 µs | 4.6× |
| `network/stoichiometric_matrix_h2_br2` | `stoichiometric_matrix` of the H₂ + Br₂ chain, 5 species × 5 reactions: the record checks, 25 zeroed entries, 16 terms | 88.2 ns (ad hoc) | 282 ns | 3.2× |
| `network/rank_michaelis_menten` | `reaction_network_rank` of E + S ⇌ ES → E + P's 4 × 3 matrix: hisab 1.4.0's Jacobi converges in a few sweeps | 594 ns (ad hoc) | 984 ns | 1.7× |
| `network/rank_chapman` | `reaction_network_rank` of Chapman's 3 × 4 cycle: 400 sweeps without convergence, Rust's ComputationError | 35.6 µs (ad hoc) | 147 µs | 4.1× |
| `sensitivity/arrhenius_gradient` | `arrhenius_gradient(1e13, 80 000, 300)`: Rust builds and sweeps a 9-node tape, the port runs the unrolled sweep (ADR 0002 #2), one correctly rounded exp | 213 ns (ad hoc) | 108 ns | 0.51× |
| `sensitivity/compute_gradient_arrhenius` | `compute_gradient` of the same function, a closure on the tape (Rust: a fresh `Tape`; the port: hisab 3.4.1's reused tape, its compensated sweep) | 202 ns (ad hoc) | 522 ns | 2.6× |

Every body but `compute_gradient` is bit-identical to Rust's on the golden rows (`network` 1,296, `fitting` 578, the four
gradients' 543 Ok rows, 0 ulp); `compute_gradient` agrees at class D.

**Reading the gap.** The least squares and the Jacobi are the port's typical 4× (each inner product a loop of loads and
f64 locals, where Rust's iterators vectorise the bounds checks away; the 3-column QR runs 9 such loops over the 20
points); `fit_arrhenius` adds six correctly rounded ln (`mathx/cr_ln` 49 ns each). The Michaelis–Menten rank is closest
because it is mostly Rust's allocation of the transposed copy, V and U. The closed-form gradient is the one body faster
than Rust: Rust allocates a tape of 9 nodes and a gradient vector per call, the port evaluates the sweep's 12
operations directly. Parity, not speed, is the bar.
