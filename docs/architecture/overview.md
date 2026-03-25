# Kimiya Architecture

## Module Map

```
kimiya
├── error.rs      — KimiyaError (5 variants)
├── element.rs    — 36 elements, lookup, Avogadro, gas constant
├── molecule.rs   — Molecule, molecular weight, formula, presets
├── reaction.rs   — Gibbs free energy, equilibrium constant, Hess's law
├── kinetics.rs   — Arrhenius, half-life, concentration decay
├── gas.rs        — Ideal gas, Van der Waals, partial pressure, density
├── solution.rs   — Molarity, pH, Henderson-Hasselbalch
└── thermo.rs     — Heat transfer, formation enthalpy, specific heat
```

## Consumers

- **ushma** — thermochemistry coupling (heat of reaction)
- **bijli** — electrochemistry (redox, Nernst equation)
- **joshua** — chemistry simulation in game worlds
- **dravya** — material chemical composition
