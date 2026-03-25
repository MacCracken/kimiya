# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.24.3] - 2026-03-25

### Added

- **electrochemistry** module — Nernst equation, Faraday's laws (1st/2nd), 19 standard electrode potentials (Li–F₂), galvanic cell EMF, cell Gibbs energy, charge/mole conversions
- **thermochem** module — 31 substances with ΔH_f°, ΔG_f°, S° (NIST/CRC); reaction enthalpy/Gibbs/entropy from stoichiometry; Van't Hoff equation; 15 Shomate Cp(T) datasets; Cp(T) numerical integration; adiabatic flame temperature
- **spectroscopy** module — Beer-Lambert law, photon energy/wavelength/frequency conversions, Bohr model energy levels and transitions, Rydberg formula, hydrogen spectral series (Lyman through Pfund)
- **organic** module — bond energy table (32 types), reaction enthalpy from bond energies, VSEPR geometry prediction (13 geometries), functional groups (16) with polarity and hydrogen bonding classification
- **kinetics** — Michaelis-Menten enzyme kinetics, Lineweaver-Burk linearization, collision theory rate, Eyring equation (transition state theory)
- **element** — expanded to full periodic table: 118 elements (H–Og) with Lanthanide and Actinide categories
- **solution** — weak acid and weak base pH solvers using hisab bisection root finding
- **hisab 1.0 integration** — bisection root finding, Simpson numerical integration, EPSILON_F64 precision constant

### Changed

- All public functions that can fail now return `Result<T, KimiyaError>` instead of silent sentinel values (0.0, f64::INFINITY)
- `KimiyaError` now derives `Clone`, `PartialEq`, `Eq`
- `KimiyaError` gained `InvalidInput` variant for general domain violations
- `Element` no longer derives `Deserialize` — `&'static str` fields cannot be deserialized from external data
- `lookup_by_number` is O(1) via direct array indexing (was O(n) linear scan)
- Temperature validation added to `ideal_gas_pressure`, `ideal_gas_volume`, `gibbs_free_energy`
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
