# Kimiya

**Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

## Features

- **Elements** — 36 elements (H–Kr) with accurate atomic masses, electronegativity, categories
- **Molecules** — composition, molecular weight, formula strings (H2O, CO2, etc.)
- **Reactions** — Gibbs free energy, equilibrium constants, Hess's law, spontaneity
- **Kinetics** — Arrhenius rate, half-life, first/second-order concentration decay
- **Gas Laws** — ideal gas (PV=nRT), Van der Waals, partial pressure, gas density
- **Solutions** — molarity, molality, dilution, pH, pOH, Henderson-Hasselbalch
- **Thermochemistry** — heat transfer (q=mcΔT), formation enthalpy, specific heat

## Quick Start

```rust
use kimiya::{element, molecule::Molecule, gas, solution};

let fe = element::lookup_by_symbol("Fe").unwrap();
println!("{}: {:.3} g/mol", fe.name, fe.atomic_mass);

let water = Molecule::water();
assert!((water.molecular_weight() - 18.015).abs() < 0.01);

let ph = solution::ph_from_h_concentration(1e-7);
assert!((ph - 7.0).abs() < 0.01); // pure water
```

## License

GPL-3.0-only
