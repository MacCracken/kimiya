# Kimiya Roadmap

## Status

**v0.24.3** — Chemistry engine with 12 modules, 118 elements, 236 tests, hisab 1.0 integration.

## Engineering Backlog

No open items.

## Domain Completeness

### Stoichiometry (new module) — CRITICAL
- Chemical equation balancing (via hisab gaussian_elimination)
- Limiting reagent identification + theoretical yield
- Oxidation state assignment
- Percent yield / percent composition
- Empirical + molecular formula from composition

### Phase Equilibria (new module) — CRITICAL
- Raoult's law (ideal solution vapor pressure)
- Henry's law (gas solubility)
- Colligative properties (bp elevation, fp depression, osmotic pressure)
- Clausius-Clapeyron equation
- Antoine equation + vapor pressure data
- Gibbs phase rule

### Solution Chemistry Extensions — CRITICAL
- Activity coefficients (Debye-Huckel limiting, extended, Davies)
- Conductivity + molar conductivity (Kohlrausch)

### Gas Law Extensions — HIGH
- Fugacity + fugacity coefficients
- Peng-Robinson equation of state
- Redlich-Kwong / SRK equation of state
- Compressibility factor (Z)
- Joule-Thomson coefficient

### Nuclear Chemistry (new module) — HIGH
- Radioactive decay law (N(t), activity)
- Decay chains (Bateman equations via hisab ODE)
- Binding energy per nucleon (semi-empirical mass formula)
- Mass defect, Q-value
- Isotope half-life table (~50 key radioisotopes)

### Inorganic Chemistry (new module) — HIGH
- Crystal field theory (d-orbital splitting, CFSE, spectrochemical series)
- Lattice energy (Born-Lande equation)
- Born-Haber cycle
- Ionic radii table (Shannon)
- Coordination geometry prediction

### Kinetics Extensions — HIGH
- Zero-order kinetics
- nth-order integrated rate law
- Steady-state approximation
- Rate law determination from data (via hisab least_squares_poly)

### Element Data Extensions — HIGH
- Ionization energies (1st–3rd)
- Electron affinities
- Atomic/ionic radii (covalent, van der Waals, Shannon ionic)
- Common oxidation states per element
- Electron configuration generation

### Data Expansion — HIGH
- Electrode potentials: 19 → ~50 (add MnO4⁻, Cr2O7²⁻, O2/H2O, Ce4+, etc.)
- Thermochem substances: 31 → ~100 (add common acids, bases, organics, metal oxides)
- Shomate Cp(T): 15 → ~40 (match thermochem table)
- Bond energies: 32 → ~60 (add P, Si, B, S bonds)

### Spectroscopy Extensions — MEDIUM
- Wavenumber conversions
- Rotational spectroscopy (rigid rotor)
- Vibrational spectroscopy (harmonic oscillator)
- Mass spectrometry isotope patterns

### Computational Chemistry Primitives — MEDIUM
- Lennard-Jones potential
- Morse potential
- Coulomb potential (point charges)
- Molecular geometry from 3D coordinates (bond lengths, angles, dihedrals)

## hisab Deepening

### Reaction Dynamics (new module) — HIGH
- Kinetics trajectory simulation (hisab rk4_trajectory, dopri45)
- Coupled reaction networks
- Michaelis-Menten time evolution (substrate depletion)
- Temperature-programmed kinetics

### Fitting (new module) — HIGH
- Arrhenius parameter fitting from experimental data (hisab levenberg_marquardt)
- Spectroscopy calibration curves (hisab least_squares_poly)
- Rate law order determination
- Shomate coefficient fitting from Cp data

### Sensitivity Analysis — MEDIUM
- d(K_eq)/dT, d(rate)/dT via hisab autodiff Dual numbers
- Uncertainty propagation through calculations
- pH sensitivity to concentration changes

### Performance Upgrades — MEDIUM
- Newton-Raphson for weak_acid_ph / weak_base_ph (3-4x speedup)
- Newton-Raphson for adiabatic_flame_temperature (4x speedup)

### Network Analysis — LOW
- Sparse stoichiometric matrices (hisab CsrMatrix)
- Reaction pathway analysis (hisab SVD)
- Equilibrium composition solver

### Spectral Processing — LOW
- FFT for Fourier-transform spectroscopy (hisab fft)
- Peak detection and fitting

## v1.0.0 Criteria

- [ ] API frozen
- [ ] Zero unwrap/panic in library code
- [ ] 90%+ test coverage (measured)
- [ ] All element data verified against IUPAC 2021
- [ ] Benchmark history with golden numbers
- [ ] 3+ downstream consumers actively depending on kimiya
