# Kimiya Roadmap

## Status

**v0.1.0** — Initial scaffold with real chemistry implementations.

## Engineering Backlog

No open items.

## Future Features (demand-gated)

### Electrochemistry
- Nernst equation (cell potential under non-standard conditions)
- Faraday's laws of electrolysis
- Standard electrode potentials table
- Galvanic cell voltage calculator

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

### Thermodynamics Integration
- Standard formation enthalpies table (common compounds)
- Entropy values table
- ushma coupling for heat flow in reactions

## v1.0.0 Criteria

- API frozen
- Zero unwrap/panic in library code
- 90%+ test coverage
- All element data verified against IUPAC
- Benchmark history with golden numbers
- 3+ downstream consumers
