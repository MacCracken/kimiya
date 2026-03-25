# Kimiya Roadmap

## Status

**v0.24.3** — Full chemistry engine with all planned features implemented.

## Completed

### Core (v0.1.0)
- [x] Element data — 36 elements (H–Kr)
- [x] Molecules — MW, formula, presets
- [x] Reactions — Gibbs, equilibrium K, spontaneity
- [x] Kinetics — Arrhenius, half-life, 1st/2nd order
- [x] Gas laws — ideal gas, Van der Waals, density
- [x] Solution chemistry — pH, molarity, Henderson-Hasselbalch
- [x] Thermochemistry — heat transfer, Hess's law, calorimetry

### Electrochemistry (v0.24.3)
- [x] Nernst equation
- [x] Faraday's laws (1st and 2nd)
- [x] Standard electrode potentials table (19 half-reactions)
- [x] Galvanic cell EMF calculator
- [x] Cell Gibbs energy

### Thermodynamics Integration (v0.24.3)
- [x] Standard formation enthalpies table (31 substances)
- [x] Standard Gibbs energy / entropy tables
- [x] Van't Hoff equation
- [x] Shomate Cp(T) polynomials (15 gases)
- [x] Cp(T) integration via hisab Simpson
- [x] Adiabatic flame temperature

### Spectroscopy (v0.24.3)
- [x] Beer-Lambert law
- [x] Photon energy/wavelength/frequency conversions
- [x] Bohr model energy levels and transitions
- [x] Rydberg wavelength formula
- [x] Hydrogen spectral series (Lyman, Balmer, Paschen, Brackett, Pfund)

### Advanced Kinetics (v0.24.3)
- [x] Michaelis-Menten enzyme kinetics
- [x] Lineweaver-Burk linearization
- [x] Collision theory rate
- [x] Transition state theory (Eyring equation)

### Organic Chemistry (v0.24.3)
- [x] Bond energy table (32 bond types)
- [x] Reaction enthalpy from bond energies
- [x] VSEPR geometry prediction (13 geometries)
- [x] Functional group enum (16 groups) with polarity/H-bonding

### Full Periodic Table (v0.24.3)
- [x] 118 elements (H–Og) with masses, electronegativity, categories
- [x] Lanthanide and Actinide categories
- [x] O(1) lookup by atomic number

### hisab Integration (v0.24.3)
- [x] Bisection root finding for weak acid/base pH
- [x] Simpson integration for Cp(T) enthalpy
- [x] EPSILON_F64 for precision checks
- [x] Bisection for adiabatic flame temperature

## v1.0.0 Criteria

- API frozen
- Zero unwrap/panic in library code
- 90%+ test coverage
- All element data verified against IUPAC
- Benchmark history with golden numbers
- 3+ downstream consumers
