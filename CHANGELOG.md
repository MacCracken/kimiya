# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.24.3] - 2026-03-25

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

[0.24.3]: https://github.com/MacCracken/kimiya/compare/v0.1.0...v0.24.3
[0.1.0]: https://github.com/MacCracken/kimiya/releases/tag/v0.1.0
