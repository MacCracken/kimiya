# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Cyrius port, toward **2.0.0** (plan: `docs/development/kimiya-port-plan.md`). Rust 1.2.0 is frozen in
`rust-old/` as the parity oracle (ADR 0001).

### Added
- **M0 initialize** — `cyrius port` scaffold on Cyrius 6.7.18; `cyrius.cyml` with hisab 3.4.1; golden-value
  generator `scripts/oracle/` (constants 45 + math 18,647 rows, glibc 2.44-pinned, `--check`, independent
  `validate.py`); generated `src/constants.cyr` (14 `KM_` physical constants); parity harness
  `tests/parity.tcyr` — 13,340 rows ok, 0 failed on x86_64 ELF, x86_64 PE and aarch64 (`docs/development/port-parity.md`)
- ADRs 0001 (port with Rust oracle), 0002 (API changes in the Cyrius 2.0 API), 0003 (transcribe the hisab 1.4.0
  routines hisab 3.4.1 lacks), 0004 (accept ganita `pow(-1, DBL_MAX)`)

### Changed
- Rust crate, CI, cargo scripts, benchmarks and architecture overview moved to `rust-old/`; CI is now Cyrius
  (`test`, `test --bounds`, `deps --locked`, `audit`)

## [1.2.0] - 2026-11-08

Refresh + audit release: current toolchain and dependencies, then a full-crate audit
and repair sweep in preparation for the Cyrius port. No public signature changed; new
items are additive. Many functions now return corrected values or reject input they
used to accept — see **Behaviour changes** before upgrading.

### Changed — toolchain and dependencies
- MSRV 1.89 → **1.99** (edition 2024); CI MSRV job updated
- hisab 1.1 → **1.4**, now `default-features = false, features = ["num", "autodiff"]` (the `transforms`, `geo` and `calc` modules were unused)
- criterion 0.5 → **0.8** (benchmarks use `std::hint::black_box`); serde 1.0.229, thiserror 2.0.21, tracing-subscriber 0.3.23, serde_json 1.0.154
- CI: actions/checkout v7, codecov-action v7, action-gh-release v3; cargo-audit and cargo-llvm-cov installed as prebuilt binaries (taiki-e/install-action)
- deny.toml: `unused-allowed-license = "allow"` (clean under cargo-deny 0.20)
- `#![warn(missing_docs)]` enabled; every public item documented with units
- README rewritten (all modules, features, dependency stack); its example now runs as a doctest

### Removed
- Unused direct `tracing` dependency (no tracing calls in kimiya; `logging` still uses tracing-subscriber)

### Added
- **gas** — `peng_robinson_z_roots` (all physical Z roots, ascending)
- **thermochem** — `lookup_shomate_at(formula, T)` and `Shomate::enthalpy_antiderivative`; multi-segment NIST Shomate data (N₂, O₂, CO₂, H₂O, H₂, CO, CH₄ to 6000 K); Zn(s), Ti(s), Ba(s), Ag(s), F₂(g) standard-state rows
- **nuclear** — `HYDROGEN_ATOM_MASS_U` (AME2020)
- **solution** — `KW_WATER_25C`
- `validate` (crate-internal) — NaN/±∞-safe guard predicates used by every module
- Benchmarks: balance_equation, oxidation states, isotope_pattern, binding energy, SEMF, Peng–Robinson Z, adiabatic flame, dihedral angle

### Fixed — physics and data
- **nuclear** — `mass_defect` mixed the bare proton mass with neutral-atom masses, so every binding energy was low by Z·mₑc² (He-4 27.27 → 28.296 MeV; Fe-56 BE/A 8.554 → 8.790 MeV); now uses m(¹H)
- **nuclear** — `AMU_MEV` mistyped (931.494 102 22 → 931.494 102 42, CODATA 2018); `PROTON_MASS_U` was CODATA 2014 → 2018
- **nuclear** — semi-empirical mass formula used Povh coefficients with a Z(Z−1) Coulomb term and δ = 12; now self-consistent (Z², δ = 11.2): Fe-56 498.97 → 494.0 MeV (exp. 492.26)
- **nuclear** — half-lives updated to NUBASE2020: I-131 8.0252 d, Pu-239 24 110 y, C-14 5700 y
- **thermochem** — N₂ and O₂ Shomate coefficients were the low-temperature NIST sets mislabelled to 1400/1200 K (N₂ Cp(2000 K) was 127 J/(mol·K); O₂ went negative); CO₂/H₂O lacked high-T segments. Adiabatic flame temperatures were ~550 K low: H₂/air 1974 → 2525 K, CH₄/air 1865 → 2327 K
- **thermochem** — C₂H₄ and Br₂(g) Shomate sets replaced with NIST's; K₂O S° 102.00 → 94.10 (CRC); n-butane row → NBS (ΔHf −126.15, ΔGf −17.03, S° 310.23)
- **thermochem** — `enthalpy_change_cp` / `adiabatic_flame_temperature` silently extrapolated outside the data range; now integrate piecewise (analytic) and return `Err` outside it. A NaN ΔH used to return `Ok(≈299.5 K)`
- **solution** — `weak_acid_ph` / `weak_base_ph` ignored water autoionisation (1e-8 M acid gave pH 8.0, now 6.978) and Newton's absolute tolerance could return the initial guess (Ka 5e-7, C 2e-6: 6.00 → 6.103); now solves the full charge-balance cubic
- **stoichiometry** — `balance_equation` truncated the system to n−1 element rows, so results depended on row order (CuSO₄ + Fe → FeSO₄ + Cu failed in one order), and hid wrong-sign solutions. Rewritten as exact integer null-space elimination with verification
- **stoichiometry** — `assign_oxidation_states` truncated fractional states (propane C −2, Fe₃O₄ +2; now `Err`) and lacked the documented peroxide, OF₂ and metal-hydride rules (H₂O₂, Na₂O₂, OF₂, NaH, LiAlH₄ now correct)
- **stoichiometry** — `empirical_formula` could pick a larger multiplier than needed or return a rounded guess; now the smallest multiplier within 0.1, else `Err`
- **inorganic** — low-spin d⁹ had 3 unpaired electrons (now 1: strong-field Cu²⁺ 3.87 → 1.73 BM); Sn²⁺ radius → Shannon CN 8 (122 pm)
- **gas** — Peng–Robinson Newton from Z = 1 converged to negative roots at high pressure (N₂ 300 K / 90 MPa: −1.65 vs 1.66) and fugacity returned NaN; now analytic cubic roots with Z > B
- **potential** — `morse_width_parameter` used the ordinary-frequency formula with an angular frequency (2π too large; HCl 11.7 → 1.87 Å⁻¹); dihedral/bond-angle degeneracy checks were absolute and rejected SI-metre coordinates
- **spectroscopy** — `isotope_pattern` panicked on large molecules (underflow → NaN → `partial_cmp().unwrap()`; e.g. C₁₆₈); hydrogen `SpectralSeries` now uses R_H (Hα 656.112 → 656.470 nm); H isotope abundances → IUPAC 2009
- **electrochemistry** — table sorted (it was two concatenated runs); K⁺/K −2.924 → −2.931, V²⁺/V −1.130 → −1.175 (CRC); `FARADAY` → exact N_A·e (96 485.332 123 310 01)
- **element** — Li 6.941 → 6.94, Si 28.086 → 28.085, S 32.065 → 32.06 (IUPAC conventional), Zr 91.224 → 91.222 (CIAAW 2024)
- **phase** — mmHg factor 133.322 → 133.322 387 415 Pa; methanol Antoine coefficients and range made consistent (NIST, 15–84 °C)
- **bridge** — used R = 8.314 and F = 96 485.332 instead of the crate constants; `atomic_number_to_valence` returned 2 for every Z ≥ 19
- **reaction_dynamics** — `simulate_stiff` (via hisab `backward_euler`) accepted unconverged Newton iterates and returned `Ok`; now its own backward-Euler step that errors on non-convergence
- **examples/basic.rs** — Nernst example used Q = [Cu²⁺] instead of 1/[Cu²⁺] (0.401 → 0.283 V)

### Fixed — robustness
- Guards written `x <= 0.0` let NaN and ±∞ through and returned `Ok(NaN)`; every validated input in every module now rejects non-finite values
- Removed every `unwrap`/`expect`/`unreachable!` from library code; integer overflow paths (u32 counts, i32 charges, `phase_rule`, SEMF `a − z`, series line index, degree + 1) use checked/saturating/widened arithmetic
- `logging::init()` panicked if a global subscriber was already installed; now uses `try_init`
- `network::rate_from_stoichiometry` silently truncated ragged rows; `reaction_dynamics` caps `steps` at 10⁷ and validates `tol`
- Many documentation errors corrected (units, return semantics, "sparse" vs dense, DOPRI45 "for stiff systems", steady-state condition, 1 atm vs 1 bar)

### Performance
Criterion medians, v1.2 refresh (criterion 0.8, pre-repair) → v1.2 post-repair; full table in `benchmarks.md`.
- thermochem/enthalpy_change_cp: 473 ns → 24.0 ns (−95%) — analytic Shomate enthalpy replaces Simpson integration
- thermochem/adiabatic_flame_ch4_air: 2.76 µs → 301 ns (−89%) — segments resolved once, analytic H(T)
- spectroscopy/isotope_pattern_glucose: 12.7 µs → 2.80 µs (−78%) — pruning + buffer reuse
- stoichiometry/oxidation_states_kmno4: 47.0 ns → 22.4 ns (−52%)
- nuclear/semi_empirical_fe56: 25.1 ns → 18.5 ns (−26%); potential/dihedral_angle: 27.2 ns → 20.4 ns (−25%)
- Cost of correctness (slower, by design):
  - solution/weak_acid_ph: 19.6 → 34.9 ns — full charge-balance cubic instead of the quadratic approximation
  - gas/peng_robinson_z_n2_1mpa: 33.2 → 61.4 ns — analytic cubic roots instead of Newton from Z = 1 (which could return an unphysical root)
  - electrochemistry/lookup_half_reaction: 29.7 → 57.1 ns — linear scan; sorting the table moved Cu²⁺/Cu later
  - Scalar formulas gain ~0.5–1.8 ns from NaN/±∞ guards (e.g. ideal_gas_pressure 1.25 → 2.18 ns, arrhenius_rate 5.80 → 6.78 ns, michaelis_menten 1.16 → 2.92 ns)
- The criterion 0.5 → 0.8 harness change alone shifted several baselines (e.g. michaelis_menten 3.04 → 1.16 ns with no code change); no claims are made for that step

### Behaviour changes (no signature change)
- Corrected values: everything under **Fixed — physics and data** above
- New `Err` cases for input that used to return `Ok`: NaN/±∞ anywhere; `balance_equation` non-integer, ambiguous or mixed-sign input; `assign_oxidation_states` non-integer average states, unknown Z, zero count; Shomate functions outside the data range; `bohr_transition_energy` with n_upper ≤ n_lower; `isotope_pattern` unknown Z; `lennard_jones` ε ≤ 0; negative vdW a/b, Henry x > 1, `phase_rule(_, 0)`, `percent_yield` negative actual
- `decay_chain_two_step` now returns the analytic limit for equal rates (was `Err`)
- `nth_order_half_life` accepts any order ≠ 1 (was order > 1 only)
- Sentinels: `bridge` functions return `NaN` (was `0.0`) for invalid input; Antoine functions return `NaN` at the pole; `semi_empirical_binding_energy` returns `NaN` for z > a or a = 0 (was a panic/wrap)
- `lookup_shomate("N2(g)")` / `("O2(g)")` return the 100 K segment; `SHOMATE_DATA` repeats formulas, one entry per segment
- `predict_geometry(1, n)` returns `Linear` (was `None`); `Molecule::formula` omits zero-count atoms; `find_peaks` threshold is relative to the data's min..max
- soorat `MolecularStructure` bonds form a star around the least electronegative non-H atom (was a chain in input order); radii derive from covalent radii

## [1.1.1] - 2026-03-29

### Fixed
- **license** — corrected SPDX identifier from `GPL-3.0` to `GPL-3.0-only` across Cargo.toml, deny.toml, README, and CLAUDE.md

## [1.1.0]

### Added
- **bridge** — cross-crate primitive-value bridges for ushma (reaction heat release, Arrhenius rate factor, equilibrium to Gibbs), tanmatra (atomic number to valence, ionization to bond energy), bijli (electrochemical current density, ionic conductivity)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `MolecularStructure` (ball-and-stick from `Molecule`), `ReactionNetworkVisualization` (species/reaction graph), `SpectrumVisualization` (x/y data with peak detection), `PhaseDiagramVisualization` (boundary curves with critical/triple points)

### Updated
- hisab 1.1.0 -> 1.3.0, zerocopy 0.8.47 -> 0.8.48

## [1.0.0] - 2026-03-25

Stable release. API frozen. 21 modules, 118 elements, 444 tests.

### Added

- **stoichiometry** module — chemical equation balancing (via hisab gaussian_elimination), oxidation state assignment, limiting reagent, theoretical yield, percent yield, percent composition, empirical formula
- **phase** module — Raoult's law, Henry's law, colligative properties (bp elevation, fp depression, osmotic pressure), Clausius-Clapeyron, Antoine equation with 8 solvents, Gibbs phase rule
- **nuclear** module — radioactive decay law, activity, mass defect, binding energy per nucleon, semi-empirical mass formula (Weizsacker), Q-value, 18 isotopes (H-3 through Np-237)
- **inorganic** module — crystal field theory (CFSE oct/tet), spectrochemical series (11 ligands), spin-only magnetic moment, Born-Landé lattice energy, 5 Madelung constants, 30 Shannon ionic radii
- **potential** module — Lennard-Jones 6-12 potential, Morse potential, Coulomb potential and force
- **reaction_dynamics** module — ODE-powered kinetics simulation (1st/2nd order, consecutive A→B→C, reversible A⇌B, Michaelis-Menten substrate depletion, adaptive DOPRI45) via hisab rk4_trajectory
- **fitting** module — Arrhenius parameter fitting, polynomial calibration curves, Beer-Lambert extinction coefficient fitting via hisab least_squares_poly
- **electrochemistry** — expanded to 49 half-reactions (added MnO₄⁻, Cr₂O₇²⁻, O₂/H₂O, Fe³⁺/Fe²⁺, Ce⁴⁺, H₂O₂, PbO₂, and 23 more)
- **thermochem** — expanded to 71 substances (added 10 elements, 23 inorganic compounds, 7 organic compounds); Shomate Cp(T) expanded to 24 datasets
- **spectroscopy** — wavenumber conversions, reduced mass, quantum harmonic oscillator energy, vibrational frequency from force constant
- **organic** — bond energies expanded to 48 types (added P-H/O/Cl, Si-H/O/C/Cl, B-H/O/N/F, S-S, S=O)
- **kinetics** — zero-order concentration and half-life, nth-order half-life, Michaelis-Menten, Lineweaver-Burk, collision theory, Eyring equation
- **gas** — Peng-Robinson equation of state (pressure, Z-factor), fugacity coefficient, fugacity, 5 substance parameter sets
- **solution** — weak acid/base pH (Newton-Raphson), ionic strength, Debye-Hückel (limiting, extended), Davies equation, Kohlrausch molar conductivity, 7 limiting ionic conductivities
- **element** — full periodic table: 118 elements (H–Og) with Lanthanide and Actinide categories
- **hisab integration** — newton_raphson (8.3x speedup for pH), rk4_trajectory, dopri45, gaussian_elimination, least_squares_poly, integral_simpson, EPSILON_F64

### Changed

- All public functions that can fail now return `Result<T, KimiyaError>` instead of silent sentinel values
- `KimiyaError` now derives `Clone`, `PartialEq`, `Eq`; added `InvalidInput` variant
- `Element` no longer derives `Deserialize` — `&'static str` fields cannot be deserialized
- `lookup_by_number` is O(1) via direct array indexing (was O(n) linear scan)
- Temperature validation added to `ideal_gas_pressure`, `ideal_gas_volume`, `gibbs_free_energy`
- Nuclear constants upgraded to CODATA 2018 full precision
- `hisab` dependency updated from 0.24 to 1.0

### Removed

- `Bond` enum from molecule module (unused dead code)
- `hess_law_enthalpy` from reaction module (duplicate of `thermo::hess_law`)

## [0.1.0] - 2026-03-24

Initial scaffold with real chemistry implementations.

### Added

- **element** — 36 elements (H–Kr), `ElementCategory`, lookup by number/symbol, Avogadro constant, gas constant
- **molecule** — `Molecule` struct, molecular weight, formula string, presets (water, CO₂, methane, glucose)
- **reaction** — Gibbs free energy, equilibrium constant, spontaneity check
- **kinetics** — Arrhenius rate constant, half-life, first/second-order concentration decay
- **gas** — ideal gas law (P and V), Van der Waals equation, partial pressure, gas density, STP constants
- **solution** — molarity, molality, dilution, pH/pOH, Henderson-Hasselbalch equation
- **thermo** — heat transfer (q = mcΔT), Hess's law, formation enthalpy, specific heat from calorimetry
- **error** — `KimiyaError` with `#[non_exhaustive]` and `thiserror` integration

[1.2.0]: https://github.com/MacCracken/kimiya/releases/tag/1.2.0
[1.0.0]: https://github.com/MacCracken/kimiya/releases/tag/v1.0.0
