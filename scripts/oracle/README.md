# scripts/oracle: what the port derives from the Rust oracle

`oracle` runs the frozen Rust kimiya 1.2.0 in [`rust-old/`](../../rust-old/) and writes
two kinds of output:
- **The golden corpus, `tests/golden/<module>.tsv`** (port plan §3 step 5). It has one row
  per call, with every float as its IEEE-754 bit pattern. The Cyrius parity harness
  (`tests/parity.tcyr`, step 7) replays every row against the port and checks it within its
  tolerance class.
- **The generated Cyrius sources, `src/constants.cyr` and (from M1) `src/<module>_data.cyr`**
  (step 6). These are the constants and tables the port reads, so no value is ever typed by
  hand.

`validate.py` writes a third kind (M3c-1, [ADR 0006](../../docs/adr/0006-correct-rounding-is-the-accuracy-target.md)):
- **The corrected oracle, `tests/golden/<module>.cr.tsv`**: the golden rows whose outputs change when the oracle's
  own operations are evaluated with every primitive result correctly rounded (MPFR, `mpfr_ref.py`). The harness
  compares those rows against the corrected outputs. See [*The corrected oracle*](#the-corrected-oracle-modulecrtsv-adr-0006).

and `mathx_tables.py` a fourth (M3c-2): **`src/mathx_data.cyr`**, the tables and constants of the port's correctly
rounded primitives, from MPFR. `cr_verify.py` checks those primitives against MPFR on every target. See
[*The correctly rounded primitives*](#the-correctly-rounded-primitives-mathx_tablespy-and-cr_verifypy-m3c-2).

The crate, the format and the rules are tanmatra's
(`/home/macro/Repos/tanmatra/scripts/oracle/`), adapted. The one structural difference is
the math library: see [*glibc*](#glibc-not-the-libm-crate).

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write both
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify; exit 1 on any difference
python3 -I scripts/oracle/validate.py .                                            # independent re-check of the corpus and the corrections
python3 -I scripts/oracle/validate.py . --write-corrections                        # write tests/golden/<module>.cr.tsv
python3 -I scripts/oracle/mpfr_ref.py                                              # the MPFR reference's self-checks
python3 -I scripts/oracle/probe_tables.py .                                        # generated .cyr through the Cyrius compiler
python3 -I scripts/oracle/mathx_tables.py                                          # write src/mathx_data.cyr (--check: verify)
python3 -I scripts/oracle/cr_verify.py                                             # the correctly rounded primitives vs MPFR, every target
python3 -I scripts/oracle/isotope_sources.py . DIR                                 # isotopes.tsv (ADR 0002 #23) against fresh downloads in DIR
```

One input is not the oracle's: **`isotopes.tsv`**, the cited source records of the isotope table (M4e-4, ADR 0002 #23;
[*The isotope table*](#the-isotope-table-a-cited-source-m4e-4-adr-0002-23)). The generator reads it, checks it against the
oracle and writes it into `src/spectroscopy_data.cyr`.

## Rules

- **`rust-old/` is never built in place or edited.** This crate depends on it by path, with
  no features (`soorat-compat` and `logging` are deferred, port plan §2.2, and change no value
  the corpus records). `Cargo.lock` was seeded from `rust-old/Cargo.lock` and then pruned by
  cargo, so every shared dependency is the version the oracle was tested with: **hisab 1.4.0**,
  serde 1.0.229, thiserror 2.0.21, … (18 packages, none new but the oracle itself). Build with
  `--locked`.
- **The output is deterministic.** Running it twice gives byte-identical files. Grid inputs use
  seeded xorshift and only exact or correctly rounded operations (+ − × ÷, `floor`,
  `str::parse::<f64>`, bit construction) — never a transcendental std method, timestamps or
  hash-map order. So the *inputs* are the same on any machine; the *outputs* of the
  transcendental rows are glibc's (next section).
- **The corpus is regenerated only from Rust 1.2.0.** A wrong oracle value is fixed in
  Rust: bump to 1.2.x, then regenerate. It is never patched in a `.tsv`. The corrected oracle doesn't patch it
  either: the `.cr.tsv` files sit beside it, are generated from it, and keep the oracle's outputs in their notes.
- **Modules arrive with their milestone.** Each port milestone starts by adding its golden
  emitter (`src/m_<module>.rs`, with the grid from the port plan §4) and its data-module
  emitter (`src/t_<module>.rs`), and registering them in `src/main.rs`.
- **`--check` flags any `.tsv` in `tests/golden/` without an emitter as stale**, and any
  corpus file recorded on another libm platform (next section). A `<module>.cr.tsv` belongs to its module's
  emitter: it is stale here only when that module has none (`validate.py` checks its content), and its header
  records the platform like the corpus's.
- **`tests/golden/modules.txt`**, also generated, lists every module. `tests/parity.tcyr`
  replays exactly that list, so an emitter added here is replayed at once, as pending until
  its functions get dispatcher entries.

| Module | Milestone | Golden rows | Golden emitter | Generated Cyrius |
|---|---|---|---|---|
| `constants` | M0 | 45 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`): the 14 physical constants, `thermo::WATER_SPECIFIC_HEAT` (M2a) and the ten `solution::` constants (M3a: K_w, Debye–Hückel A / B, Λ°). K_b/K_f and STP are generated with their module's data (M2a), the 10 Madelung and Born constants with inorganic's (M4b): every row is compared, none pending |
| `electrochemistry` | M3a (B2, B4) | 19,420 | `src/m_electrochemistry.rs`: the table and `lookup_half_reaction` (every couple, case / space variants, every prefix, the equations, `""`) (A); the cell potentials over all 49 × 49 couple pairs, ΔG, Faraday's laws over every decade of charge and mass, the electron counts 0 … 6, 127, 128, 254, 255, the edges (B); `nernst_potential` / `_25c` (C) over the edges, T and Q decades, every couple's E° and n, a seeded sweep and **near-cancellation rows** (E° = (RT/nF)·ln Q to a few ulp, E ≈ 0 of both signs); **ADR 0002 #29** (the overflow audit's B2, `audit::emit_electrochemistry`, appended with its own seed: 823 rows): R·T past DBL_MAX and a subnormal (R·T)/(n·F), E° against an overflowing product, Faraday's partial products past either end, an equivalent weight M/n underflowing, midpoints; 904 `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio` and 119 `nernst_potential` rows noted with the port's value; and B4's correction of #29 (23 rows: the product c·ln Q that IEEE rounds up to 2⁻¹⁰²² from below, `NERNST_2M1022`), so 136 `nernst_potential` rows noted | `src/electrochemistry_data.cyr` (`src/t_electrochemistry.rs`): the 49 `STANDARD_POTENTIALS` rows (a **UTF-8** equation blob, the couple blob, electrons `u8`, E° f64), `KM_STANDARD_POTENTIALS_COUNT`, `HalfReaction`'s docs and the literal's two CRC Handbook notes |
| `element` | M1 | 975 | `src/m_element.rs`: `lookup_by_number` and `lookup_properties` over every u8 (0..=255), `lookup_by_symbol` over every symbol, its lower- and upper-case forms, every name, and `""`, `"Xx"`, `"FE"`, `"Fe "`, `" Fe"`, a 300-byte string | `src/element_data.cyr` (`src/t_element.rs`): `enum KmElementCategory`, the 118 `ELEMENTS` rows (atomic number, symbol and name blobs, mass, electronegativity flag + value, category) and the 32 `ELEMENT_PROPERTIES` rows (oxidation states flattened over an `i8` table) |
| `math` | M0 step 8 (M2a: +7,500) | 26,147 | `src/m_math.rs`; M2a added dense cos / acos / ln / powf(·, 2.5) sweeps over the Peng–Robinson and Redlich–Kwong argument ranges | — |
| `molecule` | M1 | 2,472 | `src/m_molecule.rs`: the four presets, the cases the 1.2 tests pin, 600 seeded compositions (1–8 atoms, counts 0..=200, repeats, invalid Z 0/119/255, u32::MAX counts); `new`, `molecular_weight` (B), `formula`, `total_atoms` on each | — |
| `nuclear` | M4c (B3, B4) | 54,312 | `src/m_nuclear.rs`: `ISOTOPES` / `ISOTOPES[]` (output `s:symbol i:z i:a f:half_life_s s:decay_mode`), `lookup_isotope` (every label, case / space / prefix variants, near misses), `is_secular_equilibrium` (the edges against each other, λ_p at 0.01·λ_d ± 1–2 ulp, at RN(λ_d/100) ± 1 ulp over 300 signed magnitudes, the underflow and overflow ends; **ADR 0002 #14**: 208 rows noted with the exact λ_p < λ_d/100) (A); `decay_constant` (every decade, t½ = 1 giving LN_2, +∞ → 0), `activity`, `mass_defect` and `binding_energy_per_nucleon` (u32 probes 0 … u32::MAX on each side, the checked-add overflow), `binding_energy_mev`, `semi_empirical_binding_energy` (every Z ≤ A ≤ 3Z with A ≤ 300, Z = 0, the u32 extremes; **ADR 0002 #3**: the 23 NaN rows, z > a or a = 0, noted with the port's `Err(InvalidInput)`), `q_value` (`v:n` + floats twice; the −0.0 folds, ∞ − ∞, seeded) (B); `atoms_remaining`, `atoms_after_half_lives` (`2.0.powf(−n)`: integers −1100 … 1100, 2⁻ⁿ subnormal × N₀ up to DBL_MAX, the tie at n = 1075, N₀ = ±2^j, ±0, ±∞, NaN and 900 seeded n past both ends; **ADR 0002 #13**: 1,812 rows noted with the port's N₀/2ⁿ), `time_for_fraction`, `decay_chain_two_step` (output `f:N_A f:N_B f:N_C`: wide λ × λ × t grids, λ₁ = λ₂ exactly, λ₁ ≈ λ₂ either side of 10⁻¹², **pairs on the threshold to the ulp of 10⁻¹²**, e^(−λ_min·t) subnormal × N₀·λ₁, every isotope pair) (C); **ADR 0002 #32** (the overflow audit's B3, `audit::emit_nuclear`, appended with its own seed: 167 rows): `binding_energy_per_nucleon`'s Δm·931.494… past DBL_MAX with A bringing it back and not, 156 rows noted with the port's unbounded-exponent value; **ADR 0002 #33** (the overflow audit's B4, `audit::emit_nuclear_b4`, appended with its own seed: 88 rows): `q_value`'s two sums past DBL_MAX that cancel, 58 rows noted with the port's unbounded-exponent value | `src/nuclear_data.cyr` (`src/t_nuclear.rs`): the 18 `ISOTOPES` rows (label blob, `u8` Z, `u16` A, half-life, a **UTF-8** decay-mode blob) and `KM_ISOTOPES_COUNT`, with `Isotope`'s docs and the literal's three NUBASE2020 notes; the four nuclear constants are `src/constants.cyr`'s |
| `reaction` | M2a (B4) | 1,730 | `src/m_reaction.rs`: `gibbs_free_energy` (B) and `equilibrium_constant` (C) over each argument's edges, ΔG from target exponents across exp's underflow/overflow ends, the Rust tests' values, seeded sweeps; `is_spontaneous` (A); **ADR 0002 #33** (the overflow audit's B4, `audit::emit_reaction`, appended with its own seed: 172 rows): `gibbs_free_energy`'s T·ΔS past DBL_MAX with ΔH bringing it back, T·ΔS subnormal; 57 rows noted with the port's unbounded-exponent value | — |
| `solution` | M3a (B2, B4) | 29,454 | `src/m_solution.rs`: molarity, molality, dilution, ionic strength (`v:n` + `f:c i:z` pairs), Kohlrausch, conductivity (B; **ADR 0002 #29**, the overflow audit's B2: `audit::emit_solution` appends 458 edge rows with its own seed — M₁·V₁ and Λ·c past either end, midpoints, subnormal ties — and 94 `dilution` and 148 `solution_conductivity` rows are noted with the port's correctly rounded value); pH / pOH over every decade, h from pH across 10^x's ends, the activity models over every charge −4 … 4, ±50 000, i32::MIN / MAX (C); `henderson_hasselbalch` (C) with **near-cancellation rows** (pKa = −log₁₀(base/acid) to a few ulp, pH ≈ 0 of both signs); `weak_acid_ph` / `weak_base_ph` (C) over the Ka (Kb) × C grid of every decade 10⁻⁴⁰ … 10⁴⁰ × 10⁻⁴⁰ … 10²⁰, seeded pairs, the overflow cases, **each row's note naming the Newton's exit** (`[root: relative step, 4 updates]`); **ADR 0002 #33** (the overflow audit's B4, `audit::emit_solution_b4`, appended with its own seed: 430 rows): `ionic_strength`'s terms and sums past DBL_MAX, `kohlrausch_molar_conductivity`'s K·√c past DBL_MAX or subnormal, `limiting_molar_conductivity`'s products past either end; 13, 54 and 31 rows noted with the port's unbounded-exponent value | `src/constants.cyr` gains the ten solution constants (solution has no data file) |
| `stoichiometry` | M4f (M4f-2) | 47,776 | `src/m_stoichiometry.rs`: `balance_equation` (input Rust's `&[Vec<f64>]`, `v:m` then each row `v:len` and its `f:` entries; output `v:n` + `i:` coefficients: the Rust tests, CH₄ in every row order, negated and column-permuted, CuSO₄ + Fe in every row order, KMnO₄ + HCl in every row order and permuted, textbook and ionic equations, every combustion C_aH_bO_c with a ≤ 10, 30- and 70-species chains, 820 seeded consistent matrices, 700 seeded random ones, every error class, the 10⁻⁹ integer tolerance to the ulp at 0, 1, 2, 2²⁰ and 2²⁹, ±2⁵³ and past it, intermediate i64 overflow with and without a valid u32 solution, the gcd-2⁶³ row; **ADR 0002 #7**: 78 rows noted with no integer-matrix counterpart), `assign_oxidation_states` (input `v:n` + n `i:z i:count` pairs, then `i:charge`; output `v:n` + `i:z i:os` pairs: the 1.2 tests, every u8 Z alone over counts and charges, every element's binary oxides, hydrides and fluorides, oxyanions and (su)peroxides, the u32 / i32 extremes and the i32 range's edges, 3,000 seeded species; **ADR 0002 #7**: 9,090 rows noted with the port's exact rational states), `limiting_reagent` (two `v:n` slices; `i:index`), `empirical_formula` (`v:n` + `i:z f:percent` pairs → `v:n` + `i:z i:count`: every multiplier 1 … 6, real compounds, ratios on both sides of 0.1, every Err path, seeded) (A); `theoretical_yield_moles`, `percent_yield`, `percent_composition` (the molecule's atoms → `v:n` + `i:z f:percent`) over the edges, ties, over- and underflow and seeded values (B); **ADR 0002 #24–#26** (M4f-2, `emit_exact`, appended with its own seed): `percent_yield`'s exact ties at RN (subnormal: q·2⁻¹⁰⁷⁴ over 200·2^j, 3·2⁻¹⁰⁷⁴ over 200 among them; normal: 5c and 25c 54-bit odd numbers at scales across the range), its exact values, the overflow and subnormal-result ends to the ulp; `theoretical_yield_moles` where m/c over- or underflows though m·p/c is finite, and its midpoints m·3/1 and m·3/2; `limiting_reagent`'s exact ties, the near-ties whose rounded ratios are equal (both orders), both ratios overflowing or underflowing, 0 against an underflowed ratio; `empirical_formula`'s knife edges at ±1/10 (the same element twice with r = (10n ± 1)/(10k) exactly and one ulp either side; other elements with r·k within two ulps of n ± 1/10; C 12.011 % with O 17.5989 % ± 3 ulps), the subnormal percentages whose moles Rust rounds to 0, the count's 2³² edge; 1,126 `percent_yield`, 240 `theoretical_yield_moles` and 598 `percent_composition` rows noted with the port's correctly rounded value (#24), 920 `limiting_reagent` rows with its exact index (#25), 270 `empirical_formula` rows with its exact result (#26) | — |
| `thermo` | M2a (B1, B4) | 1,828 | `src/m_thermo.rs`: all class B — `heat_transfer`, `hess_law` and `enthalpy_from_formation` (slices as `v:n` + values; −0.0 sums, ∞ − ∞, cancellation, seeded slices), `specific_heat_from_calorimetry` (the 1e-12 threshold ± 1 ulp); **ADR 0002 #27** (the overflow audit's B1, `audit::emit_thermo`, appended with its own seed: 537 rows): partial products past DBL_MAX or below 2⁻¹⁰²² that the other factor brings back, the true overflow and underflow ends, 0 against an overflowing product, signs, 53-bit midpoints and the subnormal grid's ties; 151 `heat_transfer` and 171 `specific_heat_from_calorimetry` rows noted with the port's correctly rounded value; **ADR 0002 #33** (the overflow audit's B4, `audit::emit_thermo_b4`, appended with its own seed: 184 rows): partial sums past ±DBL_MAX that later steps bring back, DBL_MAX + 2⁹⁷⁰; 28 `hess_law` and 27 `enthalpy_from_formation` rows noted with the port's unbounded-exponent value | `src/constants.cyr` gains `KM_WATER_SPECIFIC_HEAT` (thermo has no data file) |
| `organic` | M4b | 9,104 | `src/m_organic.rs`: every enum value × every property — `BondType::energy_kj`, `Geometry::ideal_bond_angle`, `FunctionalGroup::pattern` / `is_polar` / `hydrogen_bonding` (the port's donor function, ADR 0002 #6), input `e:<Enum>.<Variant>` — and `predict_geometry` over the 0 … 15 square and every u8 of each argument against the other's edges (A); `enthalpy_from_bonds` (B), each list `v:n` + n pairs `e:BondType.<Variant> i:count`: empty lists, the Rust and integration tests' lists, every bond alone × 9 counts to u32::MAX on either side, self-cancellation, seeded lists | `src/organic_data.cyr` (`src/t_organic.rs`): `enum KmBondType`, `KmGeometry` (ADR 0002 #6: `KM_GEOMETRY_T_SHAPED`), `KmFunctionalGroup`, each with `KM_<ENUM>_COUNT`, and their property columns (energies with the match's `// Carbon bonds` headings, angles with the arms' notes, a **UTF-8** pattern blob, the polar / donor flags), plus ADR 0002 #6's **acceptor column, the generator's cited classification, not the oracle's** (below) |
| `potential` | M4d (M4e-2, B3) | 33,464 | `src/m_potential.rs`: `lennard_jones` (each argument over the edges, r/σ 0.5–100 at Å and metre scale, (σ/r)⁶ and (σ/r)¹² over- and underflowing, seeded; ε ≥ 2¹⁰²², ADR 0002 #17's 180 rows), `morse_width_parameter` (the edges, the HCl values, every decade, seeded; D_e ≥ 2¹⁰²³, #18's 181 rows), `coulomb`, `coulomb_force`, `distance_3d` (B); **ADR 0002 #32** (the overflow audit's B3, folding #17 and #18 in; `audit::emit_potential`, appended with its own seed: 1,703 rows — σ/r past 2^85.3 and 2^170.7 and below 2^-179, the final product's ties, the Morse width's quotient past either end, its midpoints and near-midpoints, Coulomb's k·q₁·q₂ and r·r past either end, the midpoints q·k_e and subnormal ties): 582 `lennard_jones` (the unbounded-exponent value; #17's rows with their values), 1,860 `morse_width_parameter` (the correctly rounded ω·√(μ/(2·D_e))), 730 `coulomb` and 679 `coulomb_force` rows noted with the port's value; `lj_equilibrium_distance` (σ over every binade, subnormals, the overflow end; **ADR 0002 #16**: 1,041 rows noted with the correctly rounded σ·2^(1/6)), `morse` (exp: targets across exp's range, near-cancellation rows), `bond_angle`, `dihedral_angle` (atan2: Å, metre and femtometre scale, every sign of the zero coordinates of the trans / cis tests, **atoms on the 10⁻¹² and 10⁻²⁰ tolerances to the ulp**, every Err path) (C); **ADR 0002 #20** (M4e-2): the three geometry functions' edge grids (`edges20`: 4,823 rows at 2^±(230 … 1074), mixed and subnormal coordinates, overflowing differences, the tolerances at extreme scales by bisection, tiny angles with subnormal results, trans and cis at scale) and 3,605 rows noted with the unbounded-exponent evaluation. A Rust `[f64; 3]` is three `f:` tokens | — (`COULOMB_CONSTANT` is `src/constants.cyr`'s) |
| `spectroscopy` | M4e (M4e-3, B2) | 67,869 | `src/m_spectroscopy.rs`: `isotope_abundances` (every u8: `v:n` + n (mass, abundance) float pairs), `SpectralSeries::n_lower` (input `e:SpectralSeries.<Variant>`), `find_peaks` (the data `v:n` + floats, then the threshold; output `v:k` + `i:` indices: thresholds 0, 1, −0.0, the ulps inside [0, 1] and every way out, NaN and ±∞; plateaus, non-finite data, fewer than 3 points, seeded) (A); the Beer–Lambert, photon and eV conversions, `bohr_energy_level` (every u8 Z × 20 u32 n), `bohr_transition_energy` and `rydberg_wavelength` (Z × (n₁, n₂) grids, close large n), `SpectralSeries::wavelength` (every series to u32::MAX's `checked_add` overflow), the oscillator, rotor, reduced mass and moment of inertia, **ADR 0002 #19**'s `wavenumber_from_wavelength` / `wavelength_from_wavenumber` (827 rows each noted with the correctly rounded 1/(100·x)), `isotope_pattern` (input `v:n` + n pairs `i:z i:count`, output `v:m` + m (mass, intensity) float pairs or `E:`: the Rust tests' molecules, the empty one, every u8 Z alone, every element with data × counts to 5,000, elements without data to u32::MAX, the Err position, seeded) (B); `transmittance` (exp10's whole range, the integers, subnormal results), `absorbance_from_transmittance` (log10: every decade, T near 1, the round trips), `fft_spectrum` (`v:n` + floats → `v:n/2+1` + floats: every length 0 … 16, every power of two to 2¹³, 2^k ± 1, impulses, constants, ±0, Nyquist, sinusoids (a Taylor series in + − × ÷, not std's sin), random, huge, subnormal, NaN and ±∞ samples) (C) **ADR 0002 #21** (M4e-3): every exact Bohr level (n = 2^k·Z), n₁ = n₂ − 1 at every binary scale and next to u32::MAX, the 2¹³ edge, the 13 transition ties and the exact-value families of `transition_exact_values`' exhaustive search, the series' large lines, and 6,617 rows noted with the correctly rounded exact rational; **ADR 0002 #22**: `reduced_mass`'s m = m at every power of two and subnormal multiple, its tie families, exponent gaps 45 … 64, the range ends, and 2,837 rows noted with the correctly rounded m₁·m₂/(m₁ + m₂); **ADR 0002 #23** (M4e-4): 7 `isotope_abundances` and 347 `isotope_pattern` rows noted with the port's whole result from the cited table. **ADR 0002 #30** (the overflow audit's B2, `audit::emit_spectroscopy`, appended with its own seed: 1,178 rows): k/μ past DBL_MAX or below 2⁻¹⁰⁷⁴, subnormal operands and results, exact roots and twelve roots within 2⁻¹⁰⁷ of a 53-bit midpoint, h·B and 8π²·I subnormal, partial products past either end, midpoints and subnormal ties; 195 `absorbance`, 219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant` and 118 `diatomic_moment_of_inertia` rows noted with the port's value. | `src/spectroscopy_data.cyr` (`src/t_spectroscopy.rs`): `enum KmSpectralSeries` and `KM_SPECTRAL_SERIES_N_LOWER`; the isotope abundances, since M4e-4 **ADR 0002 #23's cited table** from `isotopes.tsv` (AME2020 masses, TICE 2009 abundances, ³⁶S), checked against the oracle's `isotope_abundances(z)` **at every u8 z** (`_KM_ISOTOPE_ABUNDANCES` pairs with per-z start / count columns, each cell labelled with its published decimal); Rust's private `RYDBERG_HYDROGEN` **through the API** (below) and `ELECTRON_PROTON_MASS_RATIO` from its literal |
| `phase` | M2a (B1) | 4,174 | `src/m_phase.rs`: Raoult, Henry, colligative, osmotic (B; **ADR 0002 #27**, the overflow audit's B1: `audit::emit_phase` appends 528 edge rows with its own seed — partial products past either end, midpoints, subnormal results — and 82 + 82 colligative and 196 `osmotic_pressure` rows are noted with the port's correctly rounded value); Clausius–Clapeyron (C); Antoine mmHg / Pa (C) for every set at and one ulp past each range end, 0.5/1/10 °C outside, a 5 °C ladder, seeded points, the pole ± 1 ulp, specials, custom sets (non-finite fields, inverted range, pole inside the range, overflow), and the unit set's class-A rows that pin `PA_PER_MMHG`; `lookup_antoine`, `phase_rule` (all pairs of 13 u32 probes), `ANTOINE_DATA` / `ANTOINE_DATA[]` (A) | `src/phase_data.cyr` (`src/t_phase.rs`): `struct KmAntoineCoeffs`, the 8 sets (name blob + A, B, C, t_min_c, t_max_c columns), `_KM_PA_PER_MMHG` (through the API, see below), `KM_KB_*` / `KM_KF_*` |
| `thermochem` | M2b (B4) | 11,149 | `src/m_thermochem.rs`: the tables and lookups (A); `Shomate::cp` / `enthalpy_antiderivative` per segment incl. range ends, extrapolation, T = 0 / ±∞ / NaN, `enthalpy_change_cp` within and across segments (B); `reaction_enthalpy` / `_gibbs_energy` (B, **J/mol: the oracle's kJ/mol × 1000.0, computed by the generator**, ADR 0002 #1, said in the row note) and `reaction_entropy` over seeded species lists (`v:n` + `s:formula f:coeff` pairs); `vant_hoff_k` (C); `adiabatic_flame_temperature` (B) over fuels, equivalence ratios, oxidisers, T₀, 10⁻²⁰⁰…10³⁰⁰ scalings and every Err path, **each Ok row's note naming the solver that answered** (`[newton, k steps]`, `[bisection after newton: …, k halvings]`); **ADR 0002 #33, #34** (the overflow audit's B4, `audit::emit_thermochem`, appended with its own seed: 1,292 rows): the reaction sums' products past DBL_MAX that cancel and subnormal products; the Shomate polynomials with user coefficients at T → 0, large T, an exact cancellation at t = 2⁶⁰⁰, every binade; 115 `reaction_enthalpy` and 113 `reaction_gibbs_energy` rows noted after #1's note, 66 `reaction_entropy`, 185 `Shomate::cp` and 25 `enthalpy_antiderivative` rows noted with the port's unbounded-exponent value | `src/thermochem_data.cyr` (`src/t_thermochem.rs`): `struct KmShomate`, the 76 `THERMOCHEM_DATA` rows (formula and name blobs, ΔH_f°, ΔG_f°, S°), the 34 `SHOMATE_DATA` segments (formula blob, a … t_max), `KM_THERMOCHEM_COUNT`, `KM_SHOMATE_COUNT`, and the literals' `//` comments |
| `gas` | M2a (M3c, B3) | 22,773 | `src/m_gas.rs`: ideal gas, van der Waals, Dalton, density, `peng_robinson_pressure` (B); `peng_robinson_z`, `_z_roots` (`v:n` + roots), `fugacity_coefficient_pr`, `fugacity_pr` over each preset × 26 temperatures × 29 pressures (one- and three-root regions, the critical point ± 10⁻⁶, 90–200 MPa), seeded states, validation edges and overflow, and `redlich_kwong_pressure` (C); the `PR_*` presets and `joule_thomson_ideal` (A); **ADR 0002 #31** (the overflow audit's B3, `audit::emit_gas`, appended with its own seed: 1,680 rows — n·R·T, P·M, R·T, n/V, V − n·b, V·V, Tc·Tc, ω², Tc^2.5 past either end, subnormal Tc, κ = −∞ at T = Tc, the co-volume at the subnormal grid and where a subnormal Tc's rounding is scaled back up, the midpoints q·R and subnormal ties, Tc^2.5's ties past its range): 146 `ideal_gas_pressure`, 148 `ideal_gas_volume`, 222 `gas_density` (rounded once), 206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234 `redlich_kwong_pressure` rows (the unbounded-exponent result; 98 of them the port's Err) noted | `src/gas_data.cyr` (`src/t_gas.rs`): `struct KmPengRobinsonParams`, `KM_STP_*`, the 5 `KM_PR_*` presets as top-level struct values |
| `inorganic` | M4b (B3, B4) | 15,586 | `src/m_inorganic.rs`: `Ligand::field_strength` / `is_strong_field` (every variant), `unpaired_electrons_oct` (every u8 × both spins), `cfse_octahedral` (every u8 × both spins; ADR 0002 #12's 13 noted rows), `IONIC_RADII` / `IONIC_RADII[]` and `lookup_ionic_radius` (every table symbol × every i8; case, space, prefix and invalid symbols) (A); `delta_tet_from_oct` (since M4c Δ ≥ 2¹⁰²² too; **ADR 0002 #15**: 59 rows noted with the correctly rounded 4Δ/9 where Rust's Δ·4 overflows), `spin_only_moment` (every u8), `born_lande_lattice_energy` (each argument over the edges, 11 u32 charge probes each side, every Madelung × Born exponent and mean × charges × r₀, every tenth decade of r₀ and M, seeded; every Err path), `born_haber_lattice_energy` (B); **ADR 0002 #32** (the overflow audit's B3, `audit::emit_inorganic`, appended with its own seed: 254 rows): Born–Landé's prefix past DBL_MAX, a subnormal M, n past 2¹⁰²¹; 175 rows noted with the port's unbounded-exponent value. A `bool` is `i:0` / `i:1`; **ADR 0002 #33** (the overflow audit's B4, `audit::emit_inorganic_b4`, appended with its own seed: 170 rows): `born_haber_lattice_energy`'s differences past DBL_MAX, an odd subnormal ΔH_diss; 74 rows noted with the port's unbounded-exponent value | `src/inorganic_data.cyr` (`src/t_inorganic.rs`): `enum KmLigand` and `KM_LIGAND_FIELD_STRENGTH`; Rust's private `OCT_OCCUPANCY_*` / `OCT_UNPAIRED_*` tables **through the API** (below); the ten `KM_MADELUNG_*` / `KM_BORN_EXP_*` with their section comments; the 30 `IONIC_RADII` rows (symbol blob, `i8` charge, `u8` coordination, radius) and `KM_IONIC_RADII_COUNT`, with `IonicRadius`'s docs and the literal's Sn²⁺ note |
| `kinetics` | M3b (B1, B4) | 31,381 | `src/m_kinetics.rs`: `half_life_first_order` (every decade of k; k = 1 is LN_2 exactly), the second- and zero-order laws, `zero_order_half_life`, steady state, pre-equilibrium, Michaelis–Menten and `lineweaver_burk` (output `f:slope f:intercept`, the tuple's fields) (B), each argument over the edges and every Err path; `arrhenius_rate`, `first_order_concentration`, `collision_theory_rate`, `eyring_rate`, `eyring_rate_from_activation` (C) over T from the smallest subnormal to DBL_MAX, wide signed energies, the exponent on targets across exp's range, and **on purpose the rows where a subnormal exp result is scaled into the normal range** (A, [A]₀, Z·p up to DBL_MAX; Eyring at low T; one Eyring factor subnormal, the other huge), exp overflow (A·∞, 0·∞ → NaN); `nth_order_half_life` (C) on a fine order grid around 1 (the 10⁻¹² tolerance and the ulps either side), orders to ±300, [A]₀ and k over many decades, **[A]₀^(n−1) subnormal with a finite half-life**, and exact powers at the subnormal grid's half-way points; **ADR 0002 #28** (the overflow audit's B1, `audit::emit_kinetics`, appended with its own seed: 1,064 rows): 1/[A]₀ and k·t past either end, k ≥ 2¹⁰²³, partial products past either end, midpoints and subnormal ties, 0 against an overflow, ∞/∞; 89 `second_order_concentration`, 90 `zero_order_half_life`, 118 `steady_state_intermediate`, 179 `pre_equilibrium_rate` and 160 `michaelis_menten` rows noted with the port's value; **ADR 0002 #33** (the overflow audit's B4, `audit::emit_kinetics_b4`, appended with its own seed: 171 rows): `zero_order_concentration`'s k·t past −DBL_MAX with a negative [A]₀, k·t subnormal or below 2⁻¹⁰⁷⁵ against ±0; 61 rows noted with the port's unbounded-exponent value | — (kinetics' constants are `src/constants.cyr`'s) |

## glibc, not the `libm` crate

tanmatra's oracle called the pinned pure-Rust `libm` crate, so its corpus was the same on
every machine. kimiya 1.2.0 calls the **std f64 methods**, and on Linux those reach several
different implementations. Read off the release oracle's imports (`nm -D`), rustc 1.99.1,
x86_64-unknown-linux-gnu:

| Implementation | Functions | Bits depend on |
|---|---|---|
| **glibc libm** (dynamic) | exp, exp2, expm1, log, log10, pow, sin, cos, sincos (M4e), acos, atan2 | the glibc version, and on x86_64 the CPU: glibc picks an FMA variant (ifunc) when the CPU has FMA and AVX2 |
| **Rust compiler-builtins** (static: the `libm` crate's code) | cbrt (CORE-MATH, correctly rounded), floor, trunc, round, `__powidf2` (powi) | nothing: exact or correctly rounded |
| **inline code** | sqrt (`sqrtsd`), abs, min, max, signum, clamp, total_cmp, is_finite | nothing |

So the header of every `.tsv` records the platform, e.g.
`# oracle: kimiya 1.2.0 (rust-old/), hisab 1.4.0; libm: glibc 2.44 x86_64 fma.` (version from
`gnu_get_libc_version()`, `src/glibc.rs`). **`--check` refuses a corpus recorded on another
platform** with a `PLATFORM` line naming both, before comparing anything: its transcendental
rows can only be checked where they were made. Regenerating on another platform is a
deliberate act that changes the reference; record why.

Two more consequences, both recorded rather than assumed:
- **LLVM rewrites `2.0_f64.powf(y)` to `exp2(y)`** at opt-level ≥ 1 (the oracle binary imports
  `exp2`). kimiya writes that form for half-lives, so `math.tsv` records it as written, under
  `f64::powf_base2`; `10.0_f64.powf(y)` stays `pow` (`f64::powf_base10`). The corpus comes from
  the release-built oracle; Rust's own debug-built tests call `pow` for both.
- **LLVM joins hisab 1.4.0's FFT twiddle `angle.cos()` and `angle.sin()` into `sincos`** (M4e: the oracle binary
  imports it once `fft_spectrum` is called). `validate.py`'s replica calls Python's `math.cos` and `math.sin` and gives the
  golden FFT rows' bits on every row, so glibc's `sincos` is its `sin` and `cos` at those 13 angles (−2π/2^k, k = 1 … 13),
  where all three are correctly rounded.
- **`max`/`min` of +0 and −0 return the second operand** (`maxsd`/`minsd`). Rust leaves that
  case unspecified; the rows record what rustc 1.99.1 does, and the port's helpers must match.

`validate.py` checks every glibc row against its correctly rounded value from MPFR, so the reference
is known good on every row (the histogram is below), and the rows where glibc misrounds get the
corrected oracle (ADR 0006).

## Format 1

### Header

Each file starts with five `#` lines:
- the format and module;
- the oracle provenance (kimiya and hisab versions) and the libm platform;
- the regenerate command;
- the column list;
- the row count.

Lines starting with `#` are comments.

### Rows

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **`fn`** is the Rust path of what was called:
  - `module::function` or `Type::method` (the receiver is the first input) from M1 on;
  - `const` for a constant, whose input is its path, `s:element::AVOGADRO`;
  - `f64::<method>` for a std math reference value, `self` first (`f64::atan2` is
    `y.atan2(x)`, inputs `y x`); `f64::powf_base10` / `f64::powf_base2` are
    `10.0_f64.powf(y)` / `2.0_f64.powf(y)` with input `y`.
- **`class`** is `A`, `B`, `C` or `D` (see *Comparison*).
- **`inputs` and `outputs`** are space-separated tokens, or `-` when empty.
  - A struct is flattened into its fields in declaration order; a tuple in order.
- **`note`** is the same call in readable form, e.g. `f64::exp(-1.5) -> 0.22313016014842982`.
  It is for review only: the harness ignores it, except an ADR 0002 suffix that names the port's value (*Generated Cyrius
  sources*, below), whose trailing parenthesis holds golden tokens: one value, `(f:<bits>)]`, `(i:<n>)]` or
  `(E:<Variant>)]`, or since M4e-4 (#23) a whole result, `(v:n f:… …)]`.

### Tokens

| Token | Meaning |
|---|---|
| `f:<16 hex>` | f64, by IEEE-754 bit pattern (lowercase hex) |
| `i:<decimal>` | integer (u8 / u32 / i64 / usize); booleans are `0` / `1`; `total_cmp` is -1 / 0 / 1 |
| `e:<Enum>.<Variant>` | enum variant by its Rust name. In Cyrius, enum values follow the Rust declaration order and are append-only |
| `s:<text>` | UTF-8 string. Bytes outside printable ASCII (0x21–0x7e), and `%` itself, are `%XX` (uppercase hex) |
| `n:` | `Option::None`. `Some(x)` is just `x` |
| `E:<Variant>` | `Err(KimiyaError::<Variant>)`. Only the variant is compared, not the payload |
| `v:<len>` | a sequence of `len` elements, whose tokens follow in order. The harness knows each function's element shape |

A function returning a `Vec` too long for one row gets a row with its length (`v:<len>`) plus
one row per element under `fn[]`, whose input is the index `i:k`.

### Comparison

**Non-float tokens** (`i`, `e`, `s`, `n`, `E`, `v` lengths) must be equal in every class.

**Floats** are compared by class:

| Class | A float output passes when |
|---|---|
| **A** exact | the bits are identical. Any NaN matches any NaN, because NaN payloads and signs differ between implementations |
| **B** arithmetic | the ulp distance is ≤ 2, measured on the ordered-integer line, where +0 and −0 are the same point. **Expected 0**: a non-zero distance means the operation order differs from Rust |
| **C** transcendental | the ulp distance from the row's reference (the corrected oracle, ADR 0006) is ≤ K, a bound per function (`_par_c_bound` in `tests/parity.tcyr`), set before the function's first comparison and recorded in `docs/development/port-parity.md`. **K = 0** for every function on the port's correctly rounded primitives (all of them since M3c-3, and every `math` primitive since M4d: sin, cos, atan2 and acos on M4a's, the cos rows past \|x\| = 8 included; until then the stdlib's 2 for `f64::sin`, `f64::atan2` and those cos rows). A function with no K yet is *measured*: its histogram is reported but not gated |
| **D** iterative | \|got − want\| ≤ 1e-12 · \|want\| |

In classes B, C and D:
- NaN matches only NaN, and ±∞ match only the same infinity;
- a zero `want` needs a zero `got`, of either sign (class D; in B and C, +0 and −0 are one point
  on the ulp line).

A tolerance is never widened to make a row pass (port plan §5). A divergence accepted by an
ADR in `docs/adr/` is listed by its exact function and input field in `_par_accepted`. It is
reported as `accepted (ADR n)` on every run.

**A row with an ADR 0006 correction** (next section) is compared against the corrected outputs instead of its
golden outputs, under its class and bound. A pass is counted as *accepted (ADR 6)*, quietly. A row that fails
against the corrected outputs but matches the golden (misrounded) outputs under the same rule is *glibc-matching*
and **fails** (since M3c-3; M3c-1 counted it as measured): the port's primitive is not correctly rounded on that
path. A row that matches neither fails, and so does a correction that no golden row takes. The histograms are the
distance from each row's reference, the corrected outputs where there is one.

**Shadow lines.** The harness can also replay a row through a second Cyrius candidate
(`fn ~ candidate`, e.g. `f64::exp ~ _f64_exp_polyfill`, `f64::max ~ f64_max`). A shadow line
reports match / differ counts and the histogram and is never gated: it is evidence for the
step-8 and M1 choices, not a test. It compares against the row's reference, so since ADR 0006 a correctly
rounded candidate is 0 ulp on every row.

## The corrected oracle, `<module>.cr.tsv` (ADR 0006)

[ADR 0006](../../docs/adr/0006-correct-rounding-is-the-accuracy-target.md) makes correct rounding the accuracy
target. A golden row's reference is the oracle's own sequence of operations, as `validate.py` transcribes it from
Rust 1.2.0, evaluated with every transcendental primitive result replaced by its correctly rounded value:
`mpfr_ref.py` calls the system MPFR 4.2 through ctypes at 53 bits in binary64's exponent range, subnormalized, so
exact ties, subnormal results, overflow and the C99 Annex F special cases round as IEEE round-to-nearest-even does.
The primitives are the oracle's calls: exp, expm1 (`exp_m1`), log (`ln`), log10, pow (`powf`; `10.0.powf(y)` is
pow(10, y), whose correctly rounded value is exp10(y), checked; `x.powf(2.5)` is pow(x, 2.5)), cos, acos, exp2
(`2.0_f64.powf(y)`, LLVM's rewrite: the `math` rows' `powf_base2`, from M4c `atoms_after_half_lives`, from M4d
`lj_equilibrium_distance`'s 2^RN(1/6); the port's `_km_cr_pow(2, y)` rounds the same real 2^y), atan2 (from M4d,
potential's bond and dihedral angles), sin (the `math` rows, and from M4e with cos the FFT twiddles: the oracle's
`sincos`, the port's `_km_cr_sincos`). sqrt is IEEE already, and the oracle's cbrt is
Rust's correctly rounded one (`validate.py` proves it with rational arithmetic and against MPFR on every golden
argument; glibc's own cbrt, which the oracle never calls, misrounds about half of them). Following the transcribed
code, a corrected result may take another branch, change the root count or turn `Ok` into `Err`, and the corrected
oracle follows it.

Where glibc rounded every call of a row correctly, the corrected outputs are the golden ones and the row has no
correction. Where they differ, the row gets one, in `tests/golden/<module>.cr.tsv`:

```text
# kimiya corrected oracle (ADR 0006), format 1: solution
# oracle: kimiya 1.2.0 (rust-old/), hisab 1.4.0; libm: glibc 2.44 x86_64 fma. Correctly rounded: MPFR 4.2.2. Generated by scripts/oracle/validate.py; do not edit.
# regenerate: python3 -I scripts/oracle/validate.py . --write-corrections
# columns: fn<TAB>inputs<TAB>outputs<TAB>note (…)
# rows: 229
solution::ph_from_h_concentration	f:3f009147bb80742f	f:401200520787ae65	oracle f:401200520787ae64; log10(f:3f009147bb80742f) glibc f:c01200520787ae64 correctly rounded f:c01200520787ae65
```

- **Header:** five `#` lines, as the golden files': the format and module; the oracle provenance and libm platform
  (the golden file's) and the MPFR version (`mpfr_ref.VERSION`); the regenerate command; the columns; the row count.
- **`fn` and `inputs`** are the golden row's fields, byte for byte. They are the key: the harness finds a row's
  correction by that exact text.
- **`outputs`** are the corrected outputs, in the golden tokens (`f:`, `E:`, `v:`, `i:`).
- **`note`** is `oracle <the golden outputs>`, then `; ` and each call of the oracle's evaluation whose glibc result
  is not correctly rounded, once each in call order: `<name>(<args>) glibc f:… correctly rounded f:…`, where `name`
  is MPFR's (`exp`, `expm1`, `log`, `log10`, `pow`, `cos`, `acos`, `sin`, `atan2`, `exp2`) and the arguments are
  `f:` tokens separated by `, `.
- **Generated, never edited.** `validate.py . --write-corrections` writes one file per module with at least one
  correction, rows in golden-file order, deterministically (a module without corrections has no file, and a stale
  one is removed). It refuses to write off the corpus's recorded libm platform, since the notes hold glibc's
  results.
- **Checked on every `validate.py .` run**, which fails when a `.cr.tsv` is missing, stale, has an extra row or a row
  out of order, or names no golden module. It recomputes every correction: the outputs through the replica, and each
  listed call against MPFR directly (its correctly rounded value is MPFR's, glibc's differs from it, and on the
  recorded platform glibc's is this machine's libm result). Off that platform the notes' glibc values aren't
  compared; the rows and their outputs, which depend only on MPFR, still are.

## The correctly rounded primitives: `mathx_tables.py` and `cr_verify.py` (M3c-2)

Since M3c-3 every transcendental the port calls is correctly rounded (ADR 0006): `_km_cr_exp`, `_km_cr_expm1`,
`_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10` and `_km_cr_pow` in `src/mathx.cyr` (Ziv's two phases: a fast path with a
proven error bound and a rounding test, a triple-double accurate path), since M4a `_km_cr_sin`, `_km_cr_cos`,
`_km_cr_sincos`, `_km_cr_atan2` and `_km_cr_acos` (which replaced ADR 0005's double-double `_km_acos` and `_km_cos` in
M4d), since M4d `_km_cr_exp2_scaled` (x·2^y rounded once, ADR 0002 #13) and `_km_cr_sixth_root2_mul` (x·2^(1/6),
ADR 0002 #16), and since M4e `_km_cr_hundredth_div` (1/(100·x) rounded once, ADR 0002 #19: a double-double residual with
a rounding test, an exact 64-bit long division as its fallback). Their numbers are never typed:

- **`mathx_tables.py`** writes `src/mathx_data.cyr` — exp's table of 2^(j/512) (each as a short sum of doubles, the
  first with at most 26 significant bits), ln's table over 512 intervals of the significand (c_j with r = m·c_j − 1
  exact, and −ln c_j as a multi-word sum), the polynomial and accurate-path coefficients, the split constants (ln 2,
  ln 10 and 1/ln 10, with their Dekker halves) and exp's overflow / underflow thresholds — each value from MPFR at 400
  bits (through `mpfr_ref.py`'s binding) or exact rational arithmetic, split into doubles by round-to-nearest. It
  asserts the properties the error analysis in `src/mathx.cyr` relies on (the bit counts that make products exact, the
  26 + 27-bit splits, r exactly representable with |r| < 2⁻⁹ on every ln interval, the Fast2Sum preconditions), so a
  table that would break a proof is never written. M4a adds sin/cos(i/256), atan(j/256), the bits of 2/π (checked against
  Machin's formula in integers) and π/2's Cody–Waite parts, with the reduction's continued-fraction bound; M4d adds
  2^(1/6) as a triple-double from an exact integer root, with the continued-fraction bound that lets
  `_km_cr_sixth_root2_mul` round its product once without a test. The words are bit patterns in typed `var`s with the
  decimal value as a comment (never a hex `const`). Deterministic; `--check` regenerates in memory and exits 1 if the
  file is stale.
- **`cr_verify.py`** builds a probe program that calls each primitive (the full function, its `_km_cr_<f>_fast` entry
  and its `_km_cr_<f>_acc` accurate path) on seeded argument sets — over 1,000,000 per primitive: binade-uniform over the
  domain, the dense regions where the code changes path or rounding is delicate (subnormal results, the overflow and
  underflow ends, thresholds ± a few thousand ulps, reduction boundaries, near 1, powers of 10, the exact and midpoint
  families of pow), the known hard-to-round cases — runs it on x86_64 ELF, x86_64 PE (wine) and aarch64 (qemu), and
  compares every result bit for bit with MPFR. It also checks pow(10, x) against exp10(x), x^2.5 through `_km_cr_pow`
  (300,007 arguments), pow's C99 special values against MPFR and glibc, the fast path's fallback rate per region, and
  that the three targets agree on every bit and flag. M4a adds sin, cos, sincos, atan2 and acos (about 1.0–1.1 million
  arguments each); M4d adds `exp2s`, `_km_cr_exp2_scaled(x, y)` on 611,518 seeded pairs (nuclear's edge, subnormal and
  overflowing results, exact scalings and their ties, small and large |y|, every binade of x, the special pairs) against
  `mpfr_ref.scaled_exp2`, and drops M2a's `_km_acos` / `_km_cos` comparison with them; M4e adds `h100`,
  `_km_cr_hundredth_div(x)` on 1,049,032 arguments (every binade of both signs, the fast path's edges 2^±900, the
  overflow and subnormal-result ends, the integers, decimal wavelengths k·10⁻ʲ, the specials; 100,018 through the integer
  path alone) against `mpfr_ref.hundredth_div` (100·x exact at 64 bits, one division, MPFR's widest exponent range then
  `mpfr_check_range` and `mpfr_subnormalize`); M4e-2 (ADR 0002 #20) adds the unbounded-exponent final operations of
  `src/mathx.cyr` on (significand, exponent) pairs far outside binary64's range — `uxsqrt`, `uxmul`, `uxdiv`, `uxadd`,
  `uxatan2` (1,040,039 arguments: subnormal results and their ties, the overflow end, sticky exponent gaps,
  cancellations, every exponent difference; atan2's accurate path forced) against `mpfr_ref.py`'s `UX` (53 bits at
  MPFR's widest exponent range, the final operation through `mpfr_check_range` and `mpfr_subnormalize`) — and potential's
  `dist`, `bond`, `dih` (38,000 edge geometries each, the public function and the forced slow path `_km_*_ux`) against
  `mpfr_ref.potential_*`; M4e-3 (ADR 0002 #21, #22) adds `ratio`, `_km_cr_ratio(p, q, e)` on 320,000 triples (random sizes
  to 250 bits, results at the subnormal and overflow ends, exact quotients and midpoints), and `bohr`, `trans`, `ryd`,
  `series`, `rmass` (260,429, 555,944, 555,944, 55,345 and 502,000 arguments: the public function, the fast path and the
  exact path forced), all against exact `Fraction` arithmetic; M4f-2 (ADR 0002 #24–#26) adds stoichiometry's `muldiv`
  (`_km_cr_mul_div` = RN(x·y/z) on 640,000 triples: binade-uniform, the fast path's 2^±1021 edges and past, exact
  quotients and midpoints, values within 2⁻¹⁰⁰ of a midpoint on both sides, subnormal results and ties), `pyield`,
  `tyield`, `pcomp` (the public functions, their fast path and `_km_mul_div_acc` forced), `limit` (the public function and
  the exact comparisons alone: ties, near-ties, both ratios over- or underflowing) and `empir` (the public function and
  `_km_empirical_formula` with every tolerance test exact: formulas' own percentages, r·k at n ± 1/10 ± ulps and exactly,
  subnormal percentages, ratios near 2³²), against `Fraction`; the overflow audit's B1 (ADR 0002 #27, #28) adds `mono`
  (`_km_cr_mono` = RN(Π num / Π den) on 430,010 arguments: every shape to four factors a side and five or six with short
  ones, both signs, zeros, the result at the range ends, the 53-bit midpoints of the 3/0, 2/1, 3/1 and 2/2 shapes and the
  subnormal grid's ties, the domain's NaNs; its fast path and its exact path forced), `monotie` (`_km_mono_tie`'s exact
  sign decision, F below, beside and exactly at the midpoint) and the ten functions' public, fast and forced slow paths
  (`heat`, `shc`, `boil`, `frz`, `osm`, `zhl`, `ssi`, `pre` against `Fraction`; `soc`, `mm` against `mpfr_ref.UX`; their
  cross-target check takes a NaN word as one value, since the unvalidated products keep Rust's (m·c)·ΔT for a non-finite
  input and the hardware's default NaN of ∞·0 is negative on x86_64, positive on aarch64); B2 (ADR 0002 #29, #30) adds
  `sqrtdiv` (`_km_cr_sqrt_div` = RN(√(x/y)) on 326,845 arguments: every binade of both, k/μ past either end, subnormal
  operands and results, exact roots, roots within ulps and within 2⁻¹⁰⁰ of a 53-bit midpoint, the domain; its fast path
  and its exact path forced, against an exact integer square root) and the twelve functions' public, fast and forced
  slow paths (`dil`, `cond`, `fmass`, `fcharge`, `fratio`, `abs`, `conc`, `rotor`, `rot`, `moment` against
  `Fraction`, `vib` against the integer root, `nernst` against `mpfr_ref.UX` with MPFR's log; 3,085,916 arguments in
  all); B3 (ADR 0002 #31, #32) adds `wsqrt` (`_km_cr_mul_sqrt_ratio` = RN(w·√((x/y)·2^e)) on 316,971 arguments: every
  binade of each, e in [−4, 4] and far, the quotient and the result past either end, subnormals, exact roots, its
  midpoints at both grids, near-midpoints, the domain; its fast path and its exact path forced, against an exact integer
  square root) and the twelve functions' public, fast and forced slow paths (`igp`, `igv`, `dens`, `coul`, `force`
  against `Fraction`, `width` against the integer root, `vdw`, `prp`, `rkp`, `lj`, `lande`, `bea` against `mpfr_ref.UX`
  with MPFR's own pow for Tc^2.5; 2,815,142 arguments in all); B4 (ADR 0002 #33, #34) adds the twelve families `hess`,
  `eff`, `rsum`, `qv`, `ionic`, `haber`, `gibbs`, `kohl`, `lim`, `zoc`, `shcp`, `shh` (the public function, the
  Rust-chain fast path where its test passes, the unbounded-exponent slow path forced; against `mpfr_ref.UX` with
  `ux_sum` in Rust's summation order; partial sums and products past DBL_MAX, subnormal products, the Shomate
  polynomials at T → 0 and large T; 1,250,069 arguments) and to `nernst` the products c·ln Q that IEEE rounds up to
  2⁻¹⁰²² (ADR 0002 #29's correction, 2,800 arguments). Exit 1 on a mismatch or a cross-target difference.
  Options:
  `--targets elf`, `--only pow,exp`, `--big` (10,000,000 more per primitive on ELF), `--work DIR` (default
  `build/cr-verify`, gitignored).

`tests/mathx.tcyr` is the fast in-tree check of the same functions (special values, exact and tie cases, thresholds,
subnormal results, hard cases, every expected value MPFR's).

## Generated Cyrius sources

Values come from the oracle's public API (the one exception, since M4e-4, is ADR 0002 #23's isotope table, from its cited
source file and checked against the API). The Rust source is read only for the `///` doc
comments, which carry the citations (SI 2019, CODATA 2018, AME2020, …) into the Cyrius file,
and `t_constants::complete` proves the list is every `pub const NAME: f64` in
`rust-old/src` — no more, no fewer — before anything is written.

### `PA_PER_MMHG`: a private constant through the API

Rust's `phase::PA_PER_MMHG` is private. `t_phase.rs` reads it as `antoine_pressure_pa` of the unit set
(A = B = 0, C = 1, range [0, 0]) at T = 0: `antoine_pressure_mmhg` of that set is 10⁰ = 1.0 exactly (checked),
so the Pa value is the factor itself. Its golden row is class A, and the source literal is a cross-check.

### The flame solver's path, ADR 0002 #1's × 1000, #12's exact CFSE, #3's Err, #13–#19's and #21–#22's exact values, #20's results, #23's isotope table, #7's rational states and integer matrix, #24–#26's exact stoichiometry and #27–#34's audit fixes

Some golden rows carry a bracketed suffix in their note (`Module::row_noted`; the harness ignores notes):
- `thermochem::reaction_enthalpy` / `reaction_gibbs_energy`: `[J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002
  #1)]`. The Cyrius functions return J/mol; the generator multiplies Rust's kJ/mol value by 1000.0, the port's one
  multiply after the same subtraction, so the rows compare directly (class B, 0 ulp).
- `inorganic::cfse_octahedral` (M4b): `[ADR 0002 #12: the port returns the exact (3·e_g − 2·t2g)/5 = -1.2
  (f:bff3333333333333)]` on the 13 rows where Rust's `t2g * (-0.4) + eg * 0.6` misses that value. The row keeps the
  oracle's value; the port returns the exact one, the harness maps the two (q = round(5 × oracle), port = q / 5) and
  counts the row as accepted under ADR 2, and `validate.py` checks the note against exact rational arithmetic.
- `nuclear::semi_empirical_binding_energy` (M4c): `[ADR 0002 #3: the port returns Err(InvalidInput)]` on the 23 rows
  where Rust 1.2 returns NaN (z > a or a = 0). The row keeps the oracle's NaN; the harness maps it to the port's `Err`,
  counts it as accepted under ADR 2 and prints it by name, and `validate.py` checks that the note sits exactly on the
  NaN rows and that those are exactly the z > a or a = 0 ones.
- `nuclear::atoms_after_half_lives`, `nuclear::is_secular_equilibrium` and `inorganic::delta_tet_from_oct` (M4c, the
  user's decision): `[ADR 0002 #13: the port returns <how> = <x> (f:<bits>)]`, `[ADR 0002 #14: the port returns the exact
  λ_p < λ_d/100 = true (i:1)]`, `[ADR 0002 #15: the port returns the correctly rounded 4Δ/9 = <x> (f:<bits>)]`, on exactly
  the rows where the port's exact value differs from the oracle's (1,810, 208, 59). `<how>` says how #13's is obtained:
  `the correctly rounded N₀/2ⁿ` (an integer n), `… (2^(j − n))` (N₀ = ±2^j), `… (2⁻ᶠ to 320 bits)` (a fractional n; M4c's
  was `N₀·CR(2⁻ᶠ)·2⁻ᵏ rounded once (faithful)`) or `N₀ itself (scaleB)`. The generator computes #13 and #14 in integer
  arithmetic on the significands (`exact::round_scaled`, `m_nuclear::secular_exact`; since M4d a fractional n's 2⁻ᶠ to 320
  bits by `exact::cr_exp2_frac`, a fixed-point e^(−f·ln 2) with Ziv's two-sided test; CR(2^(j − n)) from glibc's exp2, which
  `validate.py` checks against MPFR), #15 as (Δ/9)·4. These values can't be derived from the oracle's (a 0, ±∞ or NaN
  carries none), so the harness reads them from the note (the only notes it reads) and requires the port's value bit for
  bit; `validate.py` proves each with `Fraction` and MPFR (`mpfr_ref.scaled_exp2`).
- `potential::lj_equilibrium_distance`, `potential::lennard_jones` and `potential::morse_width_parameter` (M4d):
  `[ADR 0002 #16: the port returns the correctly rounded σ·2^(1/6) = <x> (f:<bits>)]`, `[ADR 0002 #17: the port returns
  (ε·d)·4, the correctly rounded 4·ε·d for Rust's d = (σ/r)¹² − (σ/r)⁶ = <x> (f:<bits>)]` and `[ADR 0002 #18: the port
  returns ω·√((μ·½)/D_e), with the correctly rounded μ/(2·D_e) = <x> (f:<bits>)]`, on exactly the rows where the port's
  value differs (1,041, 180, 181). The generator computes #16 exactly (`exact::cr_root6_2_mul`: x·2^(1/6) > t ⟺ 2x⁶ > t⁶,
  so each candidate's midpoints are compared as integers) and #17 / #18 as the port's reorders; `validate.py` proves #16
  with the same sixth-power test in `Fraction`, #17 as RN(4·ε·d) and #18 as RN(μ/(2·D_e)). Since the overflow audit's B3
  #17 and #18 are folded into #32 (below): their rows carry #32's notes (#17's values unchanged, #18's correctly rounded).
- `spectroscopy::wavenumber_from_wavelength` and `wavelength_from_wavenumber` (M4e): `[ADR 0002 #19: the port returns the
  correctly rounded 1/(100·λ) = <x> (f:<bits>)]` (and `1/(100·ν̃)`) on exactly the rows where Rust's `0.01 / x` misses it
  (827 each). The generator computes it in 128-bit integers (`exact::cr_hundredth_div`: ⌊2¹²⁰/(100·m)⌋ and a sticky bit,
  rounded by `round_scaled`); `validate.py` proves each with `Fraction`. The rows differ by one ulp, inside class B, so the
  harness marks a noted row as mapped (`_par_force_map`) and counts it as accepted, never as a plain ok.
- `potential::distance_3d`, `bond_angle` and `dihedral_angle` (M4e-2, ADR 0002 #20): `[ADR 0002 #20: the port returns
  Rust's sequence evaluated with an unbounded exponent (√ rounded once) = <x> (f:<bits>)]`, `… (atan2 rounded once) =
  Ok(<x>) (f:<bits>)]` or `… = Err(InvalidInput) (E:InvalidInput)]`, on exactly the rows where Rust's chain is not the
  unbounded-exponent evaluation (an Err decision or the final operation's arguments differ, by the f64 replica, which is
  checked against the oracle on every row) and the result differs from the oracle's (854, 1,205, 1,546; 91 dihedral rows
  Err). The generator computes the result in exact integers (`exact::Ux`: 53-bit rounding of u128 mantissas at any
  exponent; `exact::ux_atan2_f64`: a 400-bit fixed-point atan by Euler's series with Ziv's two-sided test);
  `validate.py` evaluates every row with finite coordinates through MPFR at an unbounded exponent (`mpfr_ref.UX`,
  `mpfr_ref.potential_*`), requires each note to be that value and each unnoted row's expected result (the corrected
  oracle's where there is one) to equal it, and checks that no row whose f64 chain stays in the normal range is noted.
  The harness maps a noted row (`_par_got_noted20`).
- `spectroscopy::bohr_energy_level`, `bohr_transition_energy`, `rydberg_wavelength`, `SpectralSeries::wavelength` and
  `reduced_mass` (M4e-3, ADR 0002 #21 and #22): `[ADR 0002 #21: the port returns the correctly rounded −R·Z²/n² = <x>
  (f:<bits>)]` (likewise `R·Z²·(n₂² − n₁²)/(n₁²·n₂²)`, `n₁²·n₂²/(R∞·Z²·(n₂² − n₁²))` and `n₁²·n₂²/(R_H·(n₂² − n₁²))`) and
  `[ADR 0002 #22: the port returns the correctly rounded m₁·m₂/(m₁ + m₂) = <x> (f:<bits>)]`, on exactly the rows where the
  exact rational rounded once differs from Rust's value (3,823, 1,351, 1,376, 67 and 2,837). The generator computes it in
  `Big` integers (`exact::cr_ratio`: ⌊p·2^k/q⌋ to at least 60 bits and a sticky remainder, then `round_big`; the reduced
  mass at any exponent gap, without the port's cutoff of 55); R_H is read through the API as `t_spectroscopy` does.
  `m_spectroscopy::transition_exact_values` searches `bohr_transition_energy`'s whole domain for exact values and ties
  (n₁ = g·x, n₂ = g·y: odd(x)·odd(y) | odd(Z), odd(g)² | 683·P·j²·w, trial division below 2¹⁷ and Miller–Rabin), and the
  grid takes all 13 ties and each exact family's smallest, middle and largest member. `validate.py` proves every Ok row
  with `Fraction` (a note exactly where the rounded value differs, naming it), factors R's odd significands and repeats
  the search independently. The harness maps a noted row through `_par_got_noted` and counts #21 per function.
- `spectroscopy::isotope_abundances` and `spectroscopy::isotope_pattern` (M4e-4, ADR 0002 #23): `[ADR 0002 #23: with the
  AME2020 masses and the TICE 2009 abundances, ³⁶S included, the port returns [n], m, i, … (v:n f:… f:… …)]` on exactly the
  rows whose result changes with the cited table (7 and 347): the readable result, then its golden tokens
  (`tsv::tokens_note`) — the one note whose parenthesis holds a whole result rather than one value. The generator computes
  it with `m_spectroscopy::pattern_replica`, Rust's `isotope_pattern` transcribed with the table a parameter, which must
  give the oracle's tokens on every row with Rust's table (checked on every run). The harness parses the token list and
  requires the port's whole result to be it, bit for bit (`_par_map_toks`), counting the row as accepted under ADR 2;
  `validate.py` proves each with its replica `r_pattern` on its own parse of `isotopes.tsv`.
- `stoichiometry::assign_oxidation_states` and `stoichiometry::balance_equation` (M4f, ADR 0002 #7):
  `[ADR 0002 #7: the port returns the exact average oxidation states (Z, numerator, denominator) [n], … (v:n i:z i:p i:q
  …)]` on exactly the rows where Rust's Err is its "non-integer average oxidation state" (the generator reads the error's
  message) and the port's exact rational states are Ok (9,090), with #23's whole-result token list: three `i:` tokens per
  entry, the fraction reduced, the denominator positive. The generator computes it with `m_stoichiometry::port_oxidation`,
  Rust's algorithm on exact rationals with the i32 range on the value, which must give the oracle's result (or its Err
  variant) on every other row (asserted on every run). `[ADR 0002 #7: no integer-matrix counterpart (<why>): …]` on the
  balancing rows whose f64 matrix has no integer one (`m_stoichiometry::no_counterpart`, Rust's own acceptance test: an
  entry not finite or not within 10⁻⁹ of round(v), an integer beyond i64, ragged rows; 78), each the oracle's
  Err(InvalidInput) (asserted). The harness converts every other row's entries to Rust's rounded integers and replays
  them on the port's integer matrix; it requires a no-counterpart row's note and Err, and counts both kinds as accepted
  under ADR 2. `validate.py` recomputes both functions from Rust's semantics, the port's states with `Fraction`, every
  balancing row by exact linear algebra, and proves each note sits exactly where it should.
- `stoichiometry::percent_yield`, `theoretical_yield_moles`, `percent_composition`, `limiting_reagent` and
  `empirical_formula` (M4f-2, ADR 0002 #24–#26, the user's decisions): `[ADR 0002 #24: the port returns the correctly
  rounded 100·a/t = <x> (f:<bits>)]` (`m·p/c at the exact limiting reagent` for the yield), `[ADR 0002 #24: the port
  returns the correctly rounded mass·count·100/mw per entry [n], z, x, … (v:n i:z f:… …)]`, `[ADR 0002 #25: the port
  returns the exact limiting reagent (mᵢ·c_j < m_j·cᵢ, the first on a tie) = <k> (i:<k>)]` and `[ADR 0002 #26: with exact
  ratios and an exact 1/10 tolerance the port returns [n], z, c, … (v:n i:z i:c …)]` (or `Err(ComputationError)
  (E:ComputationError)]`), on exactly the rows where the port's result differs from the oracle's (1,126, 240, 598, 920,
  270). The generator computes them in `Big` integers, independently of the port's double-double paths:
  `exact::cr_mul_div` (RN(x·y/z) by `cr_ratio`), `exact::ratio_lt` (m₁·c₂ against m₂·c₁, the exponents aligned by a
  shift) behind `m_stoichiometry::port_limiting`, and `exact::empirical_exact` (the ratios as integer fractions at any
  size, the tolerance as 10·R against B and 9·B, no 2³² cut) behind `m_stoichiometry::port_empirical`, whose validation
  must give the oracle's Err (asserted). The harness maps a noted row (`_par_got_noted` for the two yields,
  `_par_map_exact` with an `i:` value for the index, `_par_map_toks` for the two whole results) and counts it per
  function; `validate.py` proves every row of the five with `Fraction`.
- `thermo::heat_transfer`, `specific_heat_from_calorimetry`, `phase::boiling_point_elevation`,
  `freezing_point_depression`, `osmotic_pressure` and `kinetics::zero_order_half_life`, `steady_state_intermediate`,
  `pre_equilibrium_rate`, `second_order_concentration`, `michaelis_menten` (the overflow audit's B1, ADR 0002 #27 and #28,
  the user's decisions): `[ADR 0002 #27: the port returns the correctly rounded m·c·ΔT = <x> (f:<bits>)]` (`q/(m·ΔT)`,
  `K·m·i`, `i·M·R·1000·T`), `[ADR 0002 #28: the port returns the correctly rounded [A]₀/(2k) = <x> (f:<bits>)]`
  (`k₁·[A]/k₂`, `k₂·k₁·[A]/k₋₁`) and `[ADR 0002 #28: the port returns Rust's sequence evaluated with an unbounded exponent
  (the final division rounded once) = <x> (f:<bits>)]`, on exactly the rows where the port's value differs from the
  oracle's (682 and 636). The generator computes them independently of the port (`src/audit.rs`: `exact::cr_monomial` in
  `Big` integers, `exact::Ux` with `div` and `div_f64`), with Rust's validation replicated and asserted against the oracle
  on every row, and appends each module's B1 edge rows with their own seed. The harness maps a noted row
  (`_par_got_noted`; `_par_got_noted_f` for the f64-valued ones) and counts it per function; `validate.py` §19 proves
  every row of the ten functions with `Fraction` and `mpfr_ref.UX`, and that no row of the two #20-style functions whose
  f64 chain stays normal is noted.
- `solution::dilution`, `solution_conductivity`, `electrochemistry::faraday_mass_deposited`, `faraday_charge_required`,
  `faraday_mass_ratio`, `nernst_potential` and `spectroscopy::absorbance`, `concentration_from_absorbance`,
  `vibrational_frequency`, `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia` (the overflow
  audit's B2, ADR 0002 #29 and #30): `[ADR 0002 #29: the port returns the correctly rounded M₁·V₁/V₂ = <x> (f:<bits>)]`
  (`Λ·c/1000`, `Q·M/(n·F)`, `m·n·F/M`, `m_A·M_B·n_A/(n_B·M_A)`), `[ADR 0002 #29: the port returns Rust's sequence
  evaluated with an unbounded exponent (the final subtraction rounded once) = <x> (f:<bits>)]` and `[ADR 0002 #30: the
  port returns the correctly rounded ε·l·c = <x> (f:<bits>)]` (`A/(ε·l)`, `√(k/μ)`, `h·B·J·(J + 1)`, `h/(8·π·π·I)`,
  `μ·r·r`), on exactly the rows where the port's value differs from the oracle's (2,756 and 2,106). The generator
  computes them independently of the port (`src/audit.rs`: `exact::cr_monomial`, the new `exact::cr_sqrt_ratio` and
  `exact::Ux`'s `add_f64` / `sub_f64`; the Nernst ln glibc's, correctly rounded on every row as `validate.py` checks),
  with Rust's validation replicated and asserted against the oracle on every row, and appends each module's B2 edge rows
  with their own seed. The harness maps and counts them as B1's; `validate.py` §20 proves every row of the twelve
  functions (`Fraction`; `math.isqrt` and MPFR's √ of a 256-bit quotient; `mpfr_ref.UX` with MPFR's log and
  `ux_sub_final`), that no `nernst_potential` row whose f64 chain stays normal is noted and every `nernst_potential_25c`
  row's chain is normal, and that no √ row is a rounding boundary.
- `gas::ideal_gas_pressure`, `ideal_gas_volume`, `gas_density`, `van_der_waals_pressure`, `peng_robinson_pressure`,
  `redlich_kwong_pressure` and `potential::lennard_jones`, `morse_width_parameter`, `coulomb`, `coulomb_force`,
  `inorganic::born_lande_lattice_energy`, `nuclear::binding_energy_per_nucleon` (the overflow audit's B3, ADR 0002 #31 and
  #32): `[ADR 0002 #31: the port returns the correctly rounded n·R·T/V = <x> (f:<bits>)]` (`n·R·T/P`, `P·M/(R·T)`), `[ADR
  0002 #31: the port returns Rust's sequence evaluated with an unbounded exponent (the final subtraction rounded once) =
  <x> (f:<bits>)]` or `… = Err(InvalidInput) (E:InvalidInput)]` (the port's co-volume test, on the unbounded co-volume,
  refuses where Rust's accepts; a row whose oracle output is Rust's Err can carry the port's value too), `[ADR 0002 #32:
  the port returns Rust's sequence evaluated with an unbounded exponent (the final product rounded once) = <x>
  (f:<bits>)]` (`lennard_jones`, `born_lande_lattice_energy`; `binding_energy_per_nucleon` the final division) and `[ADR
  0002 #32: the port returns the correctly rounded k_e·q₁·q₂/r = <x> (f:<bits>)]` (`k_e·q₁·q₂/(r·r)`, `ω·√(μ/(2·D_e))`),
  on exactly the rows where the port's result differs from the oracle's (1,157 and 4,182; #17's 180 and #18's 181 rows
  among #32's). The generator computes them independently of the port (`src/audit.rs`: `exact::cr_monomial`, the new
  `exact::cr_sqrt_monomial`, `exact::Ux`'s new `mul_f64` and `pow25`; the Redlich–Kwong in-range check glibc's pow,
  correctly rounded on every row as `validate.py` checks), with Rust's validation replicated and asserted against the
  oracle on every row (only the co-volume tests may differ), and appends each module's B3 edge rows with their own seed.
  The harness maps them (`_par_got_noted`, and `_par_got_noted_r` for a noted value or Err of the equations of state)
  and counts them per function; `validate.py` §21 proves every row of the twelve functions (`Fraction`; `math.isqrt` and
  MPFR's √ of a 320-bit quotient for the Morse width; `mpfr_ref.UX` with MPFR's own pow for Tc^2.5), that no row of the
  six #20-style functions whose f64 chain stays normal is noted, #17's values kept and #18's rows counted, and the
  Tc^2.5 rescale on every Redlich–Kwong row.
- `thermo::hess_law`, `enthalpy_from_formation`, `thermochem::reaction_enthalpy`, `reaction_gibbs_energy`,
  `reaction_entropy`, `nuclear::q_value`, `solution::ionic_strength`, `kohlrausch_molar_conductivity`,
  `limiting_molar_conductivity`, `inorganic::born_haber_lattice_energy`, `reaction::gibbs_free_energy`,
  `kinetics::zero_order_concentration` and `Shomate::cp`, `Shomate::enthalpy_antiderivative` (the overflow audit's B4,
  ADR 0002 #33 and #34): `[ADR 0002 #33: the port returns Rust's sequence evaluated with an unbounded exponent (the
  final addition rounded once) = <x> (f:<bits>)]` (likewise `(the final subtraction rounded once)` and `(the final
  product rounded once)` by the function's last operation; #34 the same for the Shomate polynomials), on exactly the
  rows where the port's value differs from the oracle's (697 and 210). For `reaction_enthalpy` and
  `reaction_gibbs_energy` the note follows ADR 0002 #1's: `… [J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002 #1)]
  [ADR 0002 #33: …]`. The generator computes them independently of the port (`src/audit.rs`: `exact::Ux` with Rust's
  summation orders, `.sum()` from −0.0 and thermochem's `+=` from +0.0; Rust's own expression for a non-finite input and
  for the Shomate polynomials at T = ±0), with Rust's validation replicated and asserted, and appends each module's B4
  edge rows with their own seed. The harness maps and counts them per function; `validate.py` §22 proves every row of
  the fourteen functions through `mpfr_ref.UX` and `ux_sum`, that no row whose f64 chain stays normal is noted, the
  pure-sum lemma (a finite sum of doubles is the unbounded one) partial sum by partial sum on every pure-sum row, and
  that the data-backed Shomate callers stay in the port's box (every SHOMATE_DATA segment; the flame rows' 324,925
  evaluations replayed). B4 also corrects #29 (`nernst_potential`): the 23 rows of `NERNST_2M1022` (17 noted).
- `thermochem::adiabatic_flame_temperature`: the solver that produced an `Ok` value. `m_thermochem::flame_replica`
  rebuilds the solver from the public API (`SHOMATE_DATA`, `Shomate::cp`, `Shomate::enthalpy_antiderivative`, Rust's
  private `Segments` logic, hisab 1.4.0's `newton_raphson` and `bisection` loops), and generation fails unless it
  gives the oracle's bits, or the oracle's error, on every row. `validate.py` recomputes the path independently and
  checks the note. The generator prints the tally (1,682 Newton, 154 bisection, 101 Err).

### The weak acid / base Newton's path (M3a)

`solution::weak_acid_ph` / `weak_base_ph` rows that reach the iteration carry its exit in their note:
`[root: f(x) <= 0, k updates]`, `[root: relative step, k updates]`, `[root: non-finite step, k updates]`,
`[root: 50 iterations]`, with `, checked_root rejected` if the accepted x was not positive and finite.
`m_solution::root_replica` transcribes Rust's private `charge_balance_root` and `checked_root`, and generation fails
unless it gives the oracle's pH bits, or its error, on every row. `validate.py` recomputes the path independently and
checks the note. The generator prints the tally (11,740 `f(x) <= 0`, 6,356 relative step, 74 non-finite step, no
exhaustion or rejection; at most 7 updates).

### Enums with property methods, and the acceptor column (M4b)

kimiya's `BondType`, `Geometry`, `FunctionalGroup` and `Ligand` are `#[non_exhaustive]` enums whose `const fn` methods are
tables written as `match` arms or `matches!` sets. The generator can't match them exhaustively from outside the crate, so:
- `m_organic` / `m_inorganic` list the variants, and `cyr::enum_decl` proves each list against the enum's declaration in
  the source (names, order, count) before writing `enum Km<Enum> { KM_<ENUM>_<VARIANT> = k; … }`, under the Rust docs and
  the enum body's `//` section headings. `validate.py` parses the declarations again and checks the generated enums.
- Each property is read through the API for every variant and written as a column indexed by the value;
  `cyr::match_arms` checks that the method's arms name every variant once and carries their comments — a section heading
  (`// Carbon bonds`) as a comment line inside the initializer (`cyr::table_sectioned`), a trailing note
  (`// equatorial`) as the cell's label. A formula over a column (`is_strong_field`: `field_strength() >= 8`) stays code
  in the port.
- **ADR 0002 #6** renames `Geometry::TShapedThreeLonePairs` (`KM_GEOMETRY_T_SHAPED`) and `hydrogen_bonding` (the donor
  column). The golden rows keep the Rust names; the harness maps them.
- **The acceptor column is not oracle data.** `t_organic::acceptor` classifies each group by IUPAC 2011's criterion
  (Arunan et al., *Pure Appl. Chem.* 83, 1637–1641: an electron-rich region such as a lone pair or a π-bonded pair), the
  weak acceptors included (the user's decision): every group but the alkane, with a source per family — O / N lone
  pairs (McMurry, *Organic Chemistry*, OpenStax 2023, §2.12), X–H···π (Steiner, *Angew. Chem. Int. Ed.* 41, 48 (2002)),
  divalent sulfur (Allen et al., *Acta Cryst.* B53, 696 (1997)), C–X halogens (Brammer et al., *Cryst. Growth Des.* 1,
  277 (2001); Dunitz & Taylor, *Chem. Eur. J.* 3, 89 (1997)) — written with each row's reason into the generated
  comment. It has no golden rows (ADR 0002 #6: reference tests only); `validate.py` recomputes it from the golden
  patterns by the rule (a lone-pair atom — O, N, S, a halogen — or a π bond — `=`, `≡`, an aryl ring) and
  `probe_tables.py` compiles it against the same.

### R_H through the API (M4e)

Rust's private `spectroscopy::RYDBERG_HYDROGEN` is `RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO)`, const-evaluated.
`t_spectroscopy.rs` reads it as `SpectralSeries::Lyman.wavelength(u32::MAX − 1)`: n₂ = u32::MAX, so 1/1 − 1/n₂² rounds to
1.0 (checked) and the line is RN(1/R_H); exactly one double has that reciprocal (its neighbours' differ, checked), and it
must equal the const expression evaluated in f64 on the source's literal ratio. Both are generated
(`_KM_RYDBERG_HYDROGEN`, `_KM_ELECTRON_PROTON_MASS_RATIO`); the golden Lyman row pins R_H, `probe_tables.py` compiles both
through it, and `validate.py` checks the uniqueness.

### The isotope table: a cited source (M4e-4, ADR 0002 #23)

`spectroscopy_data.cyr`'s isotope table is not the oracle's: by the user's decision it is the AME2020 masses at full
published precision and the TICE 2009 representative abundances, ³⁶S included, from **`scripts/oracle/isotopes.tsv`**.
- **The file.** A header with each source's URL, sha256 and reference (AME2020: `mass_1.mas20.txt`, Wang et al., Chinese
  Physics C 45, 030003 (2021); TICE 2009: Berglund and Wieser, Pure Appl. Chem. 83, 397–410 (2011), Table 1, column 9;
  CIAAW's current table, which gives these elements intervals only, for the record), retrieved 2026-11-10. One `isotope`
  record per isotope, tab-separated: Z, A, symbol, the mass in u with its published digits, the abundance, the TICE value
  as printed (with its uncertainty), and the AME2020 line verbatim; two `interval` records, CIAAW's standard atomic weight
  intervals of H and S. Typed from the sources, never from the oracle; edited only to follow a source.
- **The generator's checks** (`t_spectroscopy::isotope_table`, before anything is written): each mass re-read from its
  line by the file's Fortran format (the atomic mass's i3 in u at columns 107–109 and f13.6 in micro-u at 111–123), each
  abundance from its printed value; then (a) the isotopes are the oracle's `isotope_abundances(z)` at every u8 z (A its
  mass rounded, Rust's literal checked against the API's bits) plus ³⁶S, the only addition; (b) every shared abundance is
  Rust's literal as a decimal string; (c) every Rust mass is the AME2020 mass rounded to five decimals (exact decimal
  arithmetic, a tie refused); (d) each element's abundances sum to exactly 1. The emitted values are the doubles nearest
  the decimals (`str::parse`), written as round-trip decimals; each cell's label carries its published decimal.
- **Re-verifying the sources.** Download the URLs into a new directory under their basenames (`mass_1.mas20.txt`,
  `TICE-2009.pdf`, and optionally CIAAW's `hydrogen.htm` and `sulfur.htm`) and run `python3 -I
  scripts/oracle/isotope_sources.py . DIR` from the repository root. It checks the two publications' sha256s against the
  header, that every recorded AME2020 line is a whole line of the file (its only line for that Z and A) and gives the
  record's mass, that every TICE value is a cell of the PDF's Table 1 in record order (the content streams inflated with
  zlib, each text run split at its column gaps), and that each CIAAW page states the recorded interval (a page's own
  sha256 is printed, since web pages change). It reads the downloads as data and needs only the standard library.
- **validate.py** parses the file on its own and repeats (a)–(d), checks the generated table against it and each
  element's abundance-weighted mean mass against its CIAAW interval, and proves the #23 notes (above).

### Private tables through the API (M4b)

Rust's private `OCT_UNPAIRED_HIGH_SPIN` / `_LOW_SPIN` are `unpaired_electrons_oct(n, spin)` for n = 0 … 10 (Err from 11,
checked). The private `OCT_OCCUPANCY_*` (t₂g, e_g) pairs are read through `cfse_octahedral`: for each n and spin exactly
one pair with t₂g + e_g = n (t₂g ≤ 6, e_g ≤ 4) gives `t2g · (−0.4) + eg · 0.6` with the API's bits (uniqueness checked), as
`PA_PER_MMHG` is read through `antoine_pressure_pa`. All four literals are cross-checked, and generated as `_KM_` u8
columns. `probe_tables.py` probes each occupancy cell through what pins it (t₂g + e_g = n, and the formula against the
golden CFSE row).

### Emission rules

Probed on Cyrius 6.7.18 (by tanmatra; re-proven here by `probe_tables.py`), and enforced by
the generator:
- **f64s are decimal literals written positionally.** They use the shortest digits that
  round-trip, and each is checked to parse back to the oracle's bits; the lexer has no
  exponent notation. Decimal literals have been correctly rounded since 6.6.10.
- **Hex is never put in a `const`.** There it is an integer, and `CONST * x` does integer
  arithmetic.
- **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** A NaN is a generator error: Cyrius refuses a
  compile-time NaN, whose bits depend on the target.
- **Tables (from M1) are baked column arrays**, `var KM_<X>_<FIELD>: f64|u8|u16|i8|bool[N] = { .. };`,
  one element per line with a label comment (`cyr::table`, tanmatra's). `i8` arrays sign-extend on
  read (probed on x86_64 ELF/PE and aarch64), so oxidation states keep Rust's `i8`.
- **String columns are a NUL-separated `u8` blob plus a `u16` offset column** (`cyr::string_column`):
  an array initializer can't hold a string literal. `&_KM_X_BYTES + _KM_X_OFF[i]` is row i's C string;
  one blob line per row, its bytes then the text as the label (a row past 120 columns, M2b's longest names, is
  split into lines of 16 bytes with the label on the first). ASCII and NUL-free, or a generator error — except a
  column written with `cyr::string_column_utf8` (M3a: the half-reaction equations), whose blob holds each string's UTF-8
  bytes exactly (NUL-free), so `strlen` counts bytes and a C-string compare is Rust's `&str ==`.
- **Comment lines inside a table initializer** (M4b): a Rust `match`'s section headings (`// Carbon bonds`) are written
  as `# …` lines before the rows they head (`cyr::table_sectioned`); `cyrius fmt` and `lint` accept them.
- **Comments inside a table literal** (M2a's methanol source, M2b's sections and Shomate sources / caveats, M3a's two
  trailing "CRC Handbook" notes, keyed by row and field) travel into the table's comment block, keyed by row and
  word-wrapped; `t_thermochem` and `t_electrochemistry` cross-check each literal's row keys against the API's rows, in
  order.
- **`Option<f64>` cells are a `bool` presence column plus the value** (0.0 for `None`), since a
  compile-time f64 can't be NaN; the accessor returns the value-form `Option`.
- **Per-cell `//` notes in the Rust literal** (the IUPAC conventional and CIAAW 2024 atomic weights)
  travel as that cell's label comment. Only comments are read from the source.
- **The output is `cyrius fmt`- and `cyrius lint`-clean as generated.**

### Names

A bundle is concatenated into its consumers, so every generated global is prefixed `KM_`:

| Oracle / golden token | Generated Cyrius |
|---|---|
| `const` row `s:element::AVOGADRO` | `const KM_AVOGADRO` (the generator refuses two modules' constants with one name) |
| enum `ElementCategory` (M1) | `enum KmElementCategory`, plus `const KM_ELEMENT_CATEGORY_COUNT`; variant `KM_` + snake(enum) + `_` + snake(variant), as tanmatra (`e:ElementCategory.NobleGas` is `KM_ELEMENT_CATEGORY_NOBLE_GAS`). `t_element` proves the variant list against the Rust declaration |
| `ELEMENTS[i].field` (M1) | `KM_ELEMENT_<FIELD>[i]`, row i = Z − 1; symbol and name `_KM_ELEMENT_{SYMBOL,NAME}_{BYTES,OFF}`; `KM_ELEMENT_COUNT` |
| `ELEMENT_PROPERTIES[k].field` (M1) | `KM_ELEMENT_PROPERTIES_<FIELD>[k]`; oxidation states `..._OXIDATION_STATE_START` / `_COUNT` over `..._OXIDATION_STATES`; `KM_ELEMENT_PROPERTIES_COUNT` |
| `THERMOCHEM_DATA[k].field` (M2b) | `KM_THERMOCHEM_<FIELD>[k]`; formula and name `_KM_THERMOCHEM_{FORMULA,NAME}_{BYTES,OFF}`; `KM_THERMOCHEM_COUNT` |
| `SHOMATE_DATA[k]` (M2b) | `KM_SHOMATE_<FIELD>[k]` (A … T_MAX); formula `_KM_SHOMATE_FORMULA_{BYTES,OFF}`; `KM_SHOMATE_COUNT`; `Shomate` is `struct KmShomate` |
| enum `BondType` / `Geometry` / `FunctionalGroup` / `Ligand` (M4b) | `enum KmBondType` … `KmLigand`, `KM_<ENUM>_COUNT`, members by tanmatra's snake rule (`e:BondType.CCDouble` is `KM_BOND_TYPE_CC_DOUBLE`, `HCl` `KM_BOND_TYPE_H_CL`); ADR 0002 #6: `e:Geometry.TShapedThreeLonePairs` is `KM_GEOMETRY_T_SHAPED`. A method's column is `KM_<ENUM>_<METHOD>[v]` (`KM_BOND_TYPE_ENERGY_KJ`, `KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_DONOR`); the pattern `_KM_FUNCTIONAL_GROUP_PATTERN_{BYTES,OFF}` |
| `IONIC_RADII[k].field` (M4b) | `KM_IONIC_RADII_CHARGE[k]` (i8), `_COORDINATION[k]`, `_RADIUS_PM[k]`; symbol `_KM_IONIC_RADII_SYMBOL_{BYTES,OFF}`; `KM_IONIC_RADII_COUNT`; `IonicRadius` has no struct: its fields are `ionic_radius_<field>(h)`. Rust's private d-electron tables are `_KM_OCT_OCCUPANCY_{HIGH,LOW}_SPIN_{T2G,EG}` and `_KM_OCT_UNPAIRED_{HIGH,LOW}_SPIN` |
| `ISOTOPES[k].field` (M4c) | `KM_ISOTOPES_Z[k]` (u8), `_A[k]` (u16), `_HALF_LIFE_S[k]`; label `_KM_ISOTOPES_SYMBOL_{BYTES,OFF}`, decay mode `_KM_ISOTOPES_DECAY_MODE_{BYTES,OFF}` (UTF-8); `KM_ISOTOPES_COUNT`; `Isotope` has no struct: its fields are `isotope_<field>(h)` |
| `STANDARD_POTENTIALS[k].field` (M3a) | `KM_STANDARD_POTENTIALS_ELECTRONS[k]`, `KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL[k]`; equation and couple `_KM_STANDARD_POTENTIALS_{EQUATION,COUPLE}_{BYTES,OFF}`; `KM_STANDARD_POTENTIALS_COUNT`; `HalfReaction` has no struct: its fields are `half_reaction_<field>(h)` |

Each constant's comment block ends with `# Rust: kimiya::<module>::<NAME>`, which
`validate.py` and `probe_tables.py` read to tie it to its golden row.

### `probe_tables.py`

`probe_tables.py` compiles a program that prints every generated constant and every
`element_data.cyr` cell (1,055 in M1: each `ELEMENTS` field, symbol and name read through the blob
offsets, each property and oxidation state, and each category checked against the enum member the
golden variant names), and from M2a every `phase_data.cyr` / `gas_data.cyr` value (65 cells: each Antoine
name and coefficient, `KM_ANTOINE_COUNT`, `_KM_PA_PER_MMHG`, each preset field read through its struct; and
their 10 constants), and from M2b every `thermochem_data.cyr` cell (654: each THERMOCHEM_DATA formula, name and
f64 field, each SHOMATE_DATA formula and the seven columns, both counts), and from M3a every `electrochemistry_data.cyr`
cell (197: each STANDARD_POTENTIALS equation — its UTF-8 bytes — couple, electron count and E°, the count) and the ten
solution constants, and from M4b every `organic_data.cyr` cell (197: each enum member's value and count, each energy,
angle, pattern — UTF-8 — polar and donor flag against its golden row, and the acceptor flag against its rule on the
golden pattern) and `inorganic_data.cyr` cell (220: KmLigand and its strengths, each occupancy pair through t₂g + e_g = n
and the golden CFSE, each unpaired count, the 30 IONIC_RADII rows and count) and the ten inorganic constants, and from M4c every `nuclear_data.cyr` cell (91: each ISOTOPES label, Z, A, half-life
and decay mode — UTF-8 — and the count), and from M4e every `spectroscopy_data.cyr` cell (83 since M4e-4, 81 before: each KmSpectralSeries
member's value and n_lower and the count; for every u8 z the isotope count and each (mass, abundance) pair through the
start / count columns, none past `_KM_ISOTOPE_ABUNDANCE_ZS`, against the row's ADR 0002 #23 note where it has one (the
port's table, M4e-4); R_H and the mass ratio through the golden Lyman row, 1/R_H). Since M4b it also fails if a golden constant has no generated home. It reads the arrays directly, not through `src/element.cyr`, so it tests the
data while the harness tests the accessors. It builds in `build/oracle-probe/` with the project's
manifest, then compares each value bit for bit with the golden rows. It fails on any mismatch, any
compiler warning, a generated constant with no golden row, and a table it doesn't probe. A mass, an
oxidation state or a symbol byte tampered by one digit fails (checked in M1; in M2a an Antoine B, a name
byte, the mmHg factor and a preset's ω).

## `validate.py`

`validate.py` is an independent re-check of the corpus that doesn't use the Rust code. It
needs only Python 3, with no packages (run it with `-I`), and the system MPFR library (`libmpfr.so.6`, 4.x), which
`mpfr_ref.py` calls through ctypes (development only: the shipped library doesn't use it).
- **Format:** the five-line header (with the true row count), five fields per row, a known
  class, well-formed tokens (UTF-8 inside `s:`), and `modules.txt` listing exactly the `.tsv`
  files (the `.cr.tsv` files are not modules).
- **Constants:** `constants.tsv` matches every `pub const NAME: f64` literal in
  `rust-old/src/*.rs` (Python's `float()` is correctly rounded, like rustc), and every
  `const KM_X = <decimal>;` in `src/constants.cyr` is positional and parses to its row's bits.
- **math, exact functions,** recomputed with Rust's semantics: sqrt (`math.sqrt` is correctly
  rounded), abs, floor, trunc, round (half away from zero, exactly), min/max (a NaN operand is
  ignored; for ±0 only the value is checked and the observed operand rule is printed), signum,
  clamp, total_cmp, is_finite, and powi (compiler-builtins' loop). **cbrt is proven correctly
  rounded on every row with exact rational arithmetic**, which is why it is class A.
- **math, class C** (glibc's bits): every row against its **correctly rounded value from MPFR** (`mpfr_ref.py`;
  ADR 0006 replaced the 60-digit `decimal` references, which misjudge exact ties). Each must be within glibc's
  documented bound (1 ulp; log10 2 ulp), and a row where glibc is off gets a correction. cbrt's rows must also be
  MPFR's. Measured on glibc 2.44 x86_64 fma (ulp from correctly rounded; identical to the decimal references' counts):

  | Function | rows | 0 ulp | 1 ulp |
  |---|---|---|---|
  | exp, atan2, powf_base10 | 2,696 | all | 0 |
  | ln | 3,857 | 3,854 | 3 |
  | cos | 2,910 | 2,906 | 4 |
  | acos | 2,620 | 2,619 | 1 |
  | exp_m1 | 810 | 776 | 34 |
  | log10 | 1,296 | 1,248 | 48 |
  | powf | 1,750 | 1,747 | 3 (2 of them powf(x, 2.5)) |
  | powf_base2 (exp2) | 369 | 368 | 1 |
  | sin | 710 | 708 | 2 |

  (Before M2a's sweeps, ln, cos and acos had no 1-ulp row: glibc 2.44 is not correctly rounded on about
  0.1 % of arguments.)

- **element** (M1): every row against the `ELEMENTS` and `ELEMENT_PROPERTIES` literals parsed from
  `rust-old/src/element.rs` (Python's `float()` is correctly rounded, like rustc), with
  `lookup_by_symbol` recomputed as an exact, case-sensitive match; every u8 covered; and the
  generated `KmElementCategory` in the Rust declaration order.
- **molecule** (M1): every row recomputed from those masses with Rust's semantics — the molecular
  weight summed in atom order (`total += mass * count`, bit-exact: Python floats are IEEE doubles),
  the formula (input order, count 1 without a digit, zero counts skipped after the element check),
  saturating `total_atoms`, `InvalidElement` for Z outside 1–118, and the presets parsed from
  `rust-old/src/molecule.rs`.

- **reaction, thermo, phase, gas** (M2a): every row recomputed from Rust 1.2.0's code transcribed in Python
  in its operation order, with the validation guards in Rust's order and IEEE division (Python raises on a zero
  divisor). Class B rows must match bit for bit. The class-C rows call Python's `math`, which is glibc's
  exp, log, pow, acos and cos (and a correctly rounded cbrt, as Rust's compiler-builtins cbrt is), so they too
  match bit for bit when the machine is the corpus's recorded libm platform; elsewhere only their Ok/Err shape
  is checked. The tables and presets (`ANTOINE_DATA`, `PR_*`, `PA_PER_MMHG`) are parsed from the Rust source
  literals, and the generated `src/phase_data.cyr` / `src/gas_data.cyr` literals are checked against them.
- **thermochem** (M2b): every row recomputed from Rust 1.2.0's code transcribed in Python — `sum_property` from
  +0.0 with `+=`, the J/mol × 1000.0 of ADR 0002 #1, `Segments` (lookup, cp, enthalpy_change with Rust's max/min),
  the Shomate polynomial and antiderivative, `vant_hoff_k` (glibc's exp via Python), and the flame temperature with
  hisab 1.4.0's `newton_raphson` and `bisection` loops; the solver path each flame row's note claims is checked
  against the recomputation, at least 100 rows must be bisection rows, and each of Rust's typed flame constants
  (10⁻³, 10⁻⁶, the 40 J/(mol·K) guess) must change some golden row's bits when perturbed. `THERMOCHEM_DATA` and
  `SHOMATE_DATA` are parsed from the Rust literals, and the generated `src/thermochem_data.cyr` columns, string
  blobs and `struct KmShomate` are checked against them.
- **solution, electrochemistry** (M3a): every row recomputed from Rust 1.2.0's code transcribed in Python — the guards in
  Rust's order, Rust's `.sum()` from −0.0 for `ionic_strength`, the weak acid / base Newton (`charge_balance_root`,
  `checked_root`) with its exit path checked against each row's note, and `CHARGE_BALANCE_REL_TOL` (10⁻¹²) required to
  change some golden row's bits when halved and when doubled (the 50-iteration cap is reported as unobservable: no row
  enters more than 7 iterations); `STANDARD_POTENTIALS` parsed from the `electrochemistry.rs`
  literal, and the generated `src/electrochemistry_data.cyr` (the UTF-8 equation blob, the couple blob, electrons, E°,
  the two CRC notes) checked against it.
- **kinetics** (M3b): every row recomputed from Rust 1.2.0's code transcribed in Python (glibc's exp, expm1 and pow via
  Python's `math`), the guards in Rust's order, `lineweaver_burk`'s tuple as two floats; LN_2 (± 1 ulp) and
  `FIRST_ORDER_TOLERANCE` (× ½, × 2) must each change some golden row (± 1 ulp of the tolerance is reported as
  unobservable: near 1 ± 10⁻¹², |n − 1| moves in steps of 2⁻⁵² or 2⁻⁵³).
- **organic, inorganic** (M4b): every row recomputed from the Rust source — the enums parsed from their declarations
  (proving the generated enums complete and in order), each method's match arms or `matches!` set, `predict_geometry`'s
  match parsed and replayed, `enthalpy_from_bonds` as two `.sum()`s from −0.0; the d-electron tables, `IONIC_RADII` and
  the formulas of inorganic in Rust's order with its guards; every grid's coverage (every u8 × spin, every u8 moment,
  every table symbol × i8). Each typed formula constant of inorganic (−0.4, 0.6, 4, 9, the 2.0s, 1.0, the strong-field
  8, the exact CFSE's divisor 5) must change some golden row when perturbed (8 is pinned from above only: no ligand has
  rank 7). ADR 0002 #12: every Ok CFSE row's oracle value is within 10⁻⁹ of (3·e_g − 2·t₂g)/5 and the note names the
  exact value exactly where the two differ (13 rows). Both data files are checked against the source, and ADR 0002 #6's
  acceptor column against its stated rule. `VALIDATED` names every golden
  module, so a module without a section fails.
- **nuclear** (M4c): every row recomputed from Rust 1.2.0's `nuclear.rs` transcribed in its operation order, in the
  M2–M3b loop (glibc's exp, expm1, log and **exp2** via Python's `math`, the release oracle's `2.0.powf` being LLVM's
  exp2), with Rust's private constants (`EQUAL_RATE_REL_TOL`, the `SEMF_*`) parsed from the source and the 18 ISOTOPES
  rows from its literal; the grids' coverage (every Z ≤ A ≤ 3Z with A ≤ 300, every ISOTOPES row and label, u32 counts);
  ADR 0002 #3's notes (exactly the NaN rows, exactly z > a or a = 0); the chain's branches (λ₁ = λ₂ exactly, on the
  threshold, a subnormal e^(−λ_min·t)); each typed formula constant (LN_2, the tolerance — also × ½, × 2 —, the five
  SEMF coefficients, the asymmetry's 4) must change some golden row when perturbed one ulp; the
  generated `src/nuclear_data.cyr` against the literal (both blobs, the u8 / u16 / f64 columns, the NUBASE2020 notes);
  and ADR 0002 #13 / #14 / #15: every row's port value recomputed exactly (`Fraction`; MPFR's CR(2^(j − n)) and, since
  M4d, MPFR's N₀·2⁻ⁿ rounded once, `mpfr_ref.scaled_exp2`), the notes present exactly where it differs from the oracle's
  and naming it; Rust's N₀·2⁻ⁿ proven correctly rounded wherever 2⁻ⁿ is normal and n an integer; every #13 value proven
  correctly rounded by MPFR and by a 60-digit 2⁻ᶠ (M4c's faithful path is gone); seeded exact sweeps beyond the corpus
  (4,000 #13 pairs, #14's comparison replayed against
  `Fraction` on 99,940 boundary pairs, #15's (Δ/9)·4 and Rust's (Δ·4)/9 on 20,000 Δ each); and `decay_constant`'s LN_2/t½
  against the correctly rounded ln 2/t½ (kept, a diagnostic). In the inorganic section, the typed constants'
  pinning (Rust's 4.0 and 9.0 still pin its (Δ·4)/9 branch).
- **potential** (M4d): every row recomputed from Rust 1.2.0's `potential.rs` transcribed in its operation order (glibc's exp,
  atan2 — a new `_prim` kind — and exp2 via Python's `math`), the private tolerances parsed from the source; ADR 0002 #16
  proven exactly (x·2^(1/6) against each midpoint by 2x⁶ against its sixth power) on every row and 20,000 seeded σ,
  Rust's (4·ε)·d proven RN(4·ε·d) wherever 4·ε is finite (#17 and #18, M4d's reorders, are folded into #32 since the
  overflow audit's B3 and proven with it: §21), the notes exactly where the port differs; the tolerances pinned (± 1 ulp, × ½,
  × 2 each change some row), Rust's 4.0, 2.0 and the π of the −π → π mapping too; the Rust tests' geometries present, and
  the rows where atan2's −π becomes π counted.
- **spectroscopy** (M4e): every row recomputed from Rust 1.2.0's `spectroscopy.rs` transcribed in its operation order
  (glibc's pow(10, ·) and log10, and the FFT twiddles' cos and sin — a new `_prim` kind `"sin"` — via Python's `math`),
  hisab 1.4.0's `num::fft` transcribed for `fft_spectrum` (the rows prove the release oracle's `sincos` gives glibc's sin
  and cos bits), `isotope_pattern` with Python's stable sort on total_cmp's key, the private thresholds and the mass ratio
  parsed from the source, the isotope table from the `match atomic_number` literal; ADR 0002 #19 proven with `Fraction` on
  every wavenumber row (the notes exactly where Rust's 0.01/x misses RN(1/(100·x))) and on 20,000 seeded x;
  `src/spectroscopy_data.cyr` against the source (the enum and n_lower, the abundance pairs, start / count columns, R_H =
  RYDBERG / (1 + ratio) and its uniqueness through the Lyman line); the typed constants pinned where a row reaches them
  (the isotope thresholds × ½ and × 2, ħ's 2.0, the oscillator's 0.5, B's 8.0, R_H ± 1 ulp, the FFT's −2.0), the rest
  reported (the thresholds ± 1 ulp and the stable order of equal masses, which `tests/spectroscopy.tcyr` pins on the
  helpers; the rotor's 1.0, unobservable on every input); the grids' coverage; and, since M4e-3, ADR 0002 #21 and #22
  proven with `Fraction` on every Bohr, transition, Rydberg, series and reduced-mass Ok row (each note exactly where the
  exact rational rounded once differs from Rust's value, so every unnoted row is correctly rounded too), with the
  structure behind their hard rows: R's odd significands factored (no level tie, a level exact iff n/Z = 2^k, no
  wavelength dyadic), the transition energy's exhaustive search repeated (13 ties, every one in the corpus), every exact
  level present, and the reduced mass's tie rows counted; since M4e-4, ADR 0002 #23: `isotopes.tsv` parsed independently
  (each mass from its verbatim AME2020 line) and checked against Rust's literal (the same isotopes plus ³⁶S, the shared
  abundances as decimal strings, Rust's masses the AME2020 ones to five decimals, each element summing to exactly 1), the
  mean masses against CIAAW's intervals, the generated table against the file, every isotope row's port result
  (`r_pattern` with the table swapped) noted exactly where it differs from the oracle's and named by the note bit for
  bit, and the thresholds' pinning repeated against the port's results.
- **stoichiometry** (M4f, §18): every row recomputed from Rust 1.2.0's `stoichiometry.rs` — the integer functions on
  Python integers with i64 checks at every Rust `checked_*` call (i128 in the verification), the float functions on IEEE
  doubles in Rust's order; ADR 0002 #7's rational states with `Fraction` and every balancing row by exact linear algebra;
  EMPIRICAL_TOLERANCE, INTEGER_TOLERANCE and MAX_COMPOSITION_ENTRY pinned. Since M4f-2, ADR 0002 #24–#26: the port's
  result of every `percent_yield`, `theoretical_yield_moles`, `percent_composition`, `limiting_reagent` and
  `empirical_formula` row recomputed with `Fraction` (RN(100·a/t), RN(m·p/c) at the exact least m/c, RN(mass·count·100/mw)
  with Rust's mw, the exact index, the exact empirical formula) and each note proven to sit exactly where it differs from
  the oracle's, naming it; the port's claims checked on the corpus (every row with a ratio of 2³² − ½ or more is the
  ComputationError, the kept ratios lie in [1, 2³²)); Rust's empirical Err paths tabulated against the port's results;
  the hard rows counted (percent_yield's ties, subnormal ones included, and range ends; the yield's over- and underflowing
  m/c and its midpoints; limiting_reagent's exact and rounded ties, both overflowing, both underflowing; the rows exactly
  1/10 off); the exact 1/10 pinned (1/10 − 2⁻⁸⁰ and the nearest larger tolerance a row needs each change some row).
- **The corrected oracle** (M3c-1, ADR 0006): every row whose replica calls a primitive is evaluated again with each
  primitive result correctly rounded (the `math` rows of the primitives are their MPFR values), and the corrections
  are checked against `tests/golden/<module>.cr.tsv` as [above](#the-corrected-oracle-modulecrtsv-adr-0006).
  `henderson_hasselbalch`'s oracle call is glibc's log10 itself, so the corrected oracle replaces it whole (M3a's
  check that glibc's log10 is the e_log10.c composition, the basis of the retired `_km_log10_ln`, went with it in
  M3c-3). It prints the corrected rows per module
  and function, the rows whose misrounded calls leave the outputs unchanged, the per-primitive counts of glibc's
  misroundings over every replica call, and glibc's own cbrt on the golden cbrt arguments
  (`docs/development/port-parity.md` §9).
- **`--sensitivity`**, a diagnostic of conditioning since ADR 0006 (M3c-3): `python3 -I scripts/oracle/validate.py .
  --sensitivity` re-evaluates every class-C row with one primitive result nudged one ulp up and down, one call at a
  time — every result but 0, ±∞, NaN and acos(±1) — and prints each function's S (the largest output change, in ulp),
  c (the most nudged calls in one row), the retired port-parity rule's K = 2·S·c + 2 for comparison, the calls where
  glibc's result is not correctly rounded (MPFR decides) and branch flips. For each function whose K passes 64 it lists
  the rows that put it there, how many nudge a subnormal primitive result, the worst three, and S over the rest. It sets
  no bound: the port evaluates the oracle's operations with correctly rounded primitives, as the corrected oracle does,
  so every class-C row is gated at 0 ulp whatever its S (`docs/development/port-parity.md` §10). Until M3c-3 it nudged
  only the calls whose result could differ between port and oracle (ADR 0005) and set the bounds of §5–§8; the
  `_km_log10_ln` port model and the "S before" mode are gone with the primitives they modelled.

Every function in every `.tsv` must have a rule in `validate.py`: a new one fails until it gets one. Extend it
as modules arrive. A check that recomputes a row independently is worth more than one that re-reads the
oracle.
