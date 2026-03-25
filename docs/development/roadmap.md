# Kimiya Roadmap

## Status

**v0.24.3** — Chemistry engine with 17 modules, 118 elements, 325 tests, 23 benchmarks.

## Engineering Backlog

No open items.

## Remaining Domain Gaps

### Stoichiometry — CRITICAL
- Chemical equation balancing (via hisab gaussian_elimination)
- Oxidation state assignment

### Gas Extensions — HIGH
- Fugacity + fugacity coefficients
- Redlich-Kwong / SRK equation of state
- Joule-Thomson coefficient

### Inorganic Chemistry (new module) — HIGH
- Crystal field theory (d-orbital splitting, CFSE, spectrochemical series)
- Lattice energy (Born-Lande equation)
- Born-Haber cycle
- Ionic radii table (Shannon)
- Coordination geometry prediction

### Element Data Extensions — HIGH
- Ionization energies (1st–3rd)
- Electron affinities
- Atomic/ionic radii (covalent, van der Waals, Shannon ionic)
- Common oxidation states per element
- Electron configuration generation

### Data Expansion — HIGH
- Electrode potentials: 19 → ~50
- Thermochem substances: 31 → ~100
- Shomate Cp(T): 15 → ~40
- Bond energies: 32 → ~60

### Kinetics Extensions — HIGH
- Zero-order kinetics
- nth-order integrated rate law
- Steady-state approximation

### Nuclear Extensions — MEDIUM
- Decay chains (Bateman equations via hisab ODE)
- Q-value of nuclear reactions

### Solution Extensions — MEDIUM
- Conductivity + molar conductivity (Kohlrausch)

### Spectroscopy Extensions — MEDIUM
- Wavenumber conversions
- Rotational spectroscopy (rigid rotor)
- Vibrational spectroscopy (harmonic oscillator)
- Mass spectrometry isotope patterns

### Computational Chemistry Primitives — MEDIUM
- Lennard-Jones potential
- Morse potential
- Coulomb potential (point charges)
- Molecular geometry from 3D coordinates

## hisab Deepening (blocked on hisab 1.2+)

### Sensitivity Analysis — MEDIUM
- Reverse-mode autodiff for reaction network gradients
- Uncertainty propagation

### Stiff ODE Solvers — HIGH
- backward_euler / bdf2 for stiff kinetics
- SDE solvers for stochastic chemical kinetics

### Network Analysis — LOW
- Sparse stoichiometric matrices (CsrMatrix)
- Reaction pathway analysis (SVD)
- Equilibrium composition solver

### Spectral Processing — LOW
- FFT for Fourier-transform spectroscopy
- Peak detection and fitting

## v1.0.0 Criteria

- [ ] API frozen
- [ ] Zero unwrap/panic in library code
- [ ] 90%+ test coverage (measured)
- [ ] All element data verified against IUPAC 2021
- [ ] Benchmark history with golden numbers
- [ ] 3+ downstream consumers actively depending on kimiya
