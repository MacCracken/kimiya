# Kimiya Roadmap

## Status

**v0.1.0** — Initial scaffold with real chemistry implementations, electrochemistry, and thermochemical data.

## Engineering Backlog

No open items.

## Completed

### Electrochemistry (v0.1.0)
- [x] Nernst equation (cell potential under non-standard conditions)
- [x] Faraday's laws of electrolysis (1st and 2nd)
- [x] Standard electrode potentials table (19 half-reactions)
- [x] Galvanic cell voltage calculator
- [x] Cell Gibbs energy

### Thermodynamics Integration (v0.1.0)
- [x] Standard formation enthalpies table (31 substances)
- [x] Standard Gibbs energy of formation table
- [x] Entropy values table
- [x] Van't Hoff equation (K vs temperature)
- [x] Reaction enthalpy/Gibbs/entropy from stoichiometry

## Future Features (demand-gated)

### Organic Chemistry
- Functional group recognition
- IUPAC naming (simple molecules)
- Molecular geometry (VSEPR)
- Bond energy table

### Spectroscopy
- Beer-Lambert law (absorbance = ε × l × c)
- Emission/absorption wavelengths for elements
- Energy level transitions (Bohr model for hydrogen)

### Advanced Kinetics
- Michaelis-Menten enzyme kinetics
- Collision theory rate estimation
- Transition state theory

### Thermodynamics Integration (continued)
- ushma coupling for heat flow in reactions
- Heat capacity as f(T) — Shomate polynomials
- Adiabatic flame temperature

## v1.0.0 Criteria

- API frozen
- Zero unwrap/panic in library code
- 90%+ test coverage
- All element data verified against IUPAC
- Benchmark history with golden numbers
- 3+ downstream consumers
