# Kimiya Roadmap

## Status

**v0.24.3** — Chemistry engine with 19 modules, 118 elements, 411 tests, 23 benchmarks.

## Engineering Backlog

No open items.

## Remaining Domain Gaps

### Kinetics — MEDIUM
- Steady-state approximation for intermediates

### Spectroscopy — MEDIUM
- Mass spectrometry isotope patterns

### Inorganic — MEDIUM
- Born-Haber cycle calculator
- Coordination geometry prediction

### Nuclear — MEDIUM
- Decay chains (Bateman equations via hisab ODE)

## hisab Deepening (blocked on hisab 1.1 publish)

### Stiff ODE Solvers — HIGH
- backward_euler / bdf2 for stiff kinetics
- SDE solvers for stochastic chemical kinetics

### Sensitivity Analysis — MEDIUM
- Reverse-mode autodiff for reaction network gradients
- Uncertainty propagation

### Network Analysis — LOW
- Sparse stoichiometric matrices (CsrMatrix)
- Reaction pathway analysis (SVD)

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
