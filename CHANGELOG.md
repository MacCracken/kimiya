# Changelog

## [Unreleased]

### Added
- **electrochemistry** — Nernst equation, Faraday's laws (1st/2nd), 19 standard electrode potentials (Li–F₂), galvanic cell EMF, cell Gibbs energy, charge/mole conversions
- **thermochem** — 31 substances with ΔH_f°, ΔG_f°, S° (NIST/CRC), reaction enthalpy/Gibbs/entropy from stoichiometry, Van't Hoff equation

### Changed
- All public functions that can fail now return `Result<T, KimiyaError>` instead of silent sentinel values
- `KimiyaError` now derives `Clone`, `PartialEq`, `Eq`
- Added `InvalidInput` variant to `KimiyaError`
- `Element` no longer derives `Deserialize` (`&'static str` fields cannot be deserialized)
- `lookup_by_number` is now O(1) via direct indexing (was O(n) linear scan)
- Temperature validation added to `ideal_gas_pressure`, `ideal_gas_volume`, `gibbs_free_energy`

### Removed
- Unused `Bond` enum from molecule module
- Duplicate `hess_law_enthalpy` from reaction module (kept `thermo::hess_law`)

## [0.1.0] - 2026-03-24

Initial scaffold with real chemistry implementations.

### Modules
- **element** — 36 elements (H–Kr), ElementCategory, lookup by number/symbol, Avogadro, gas constant
- **molecule** — Molecule, molecular weight, formula string, presets (water, CO2, methane, glucose)
- **reaction** — Gibbs free energy, equilibrium constant, spontaneity check
- **kinetics** — Arrhenius rate, half-life, first/second-order concentration
- **gas** — ideal gas law, Van der Waals, partial pressure, gas density, STP constants
- **solution** — molarity, molality, dilution, pH, pOH, Henderson-Hasselbalch
- **thermo** — heat transfer (q=mcΔT), Hess's law, formation enthalpy, specific heat from calorimetry
