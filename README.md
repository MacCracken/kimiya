# Kimiya

**Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

> **Porting to Cyrius.** The Rust crate's last release, **1.2.0** (tag `1.2.0`), is frozen in
> [`rust-old/`](rust-old/) as the parity oracle; the Cyrius library will ship as **2.0.0**. The
> Rust documentation below describes 1.2.0. Progress: [`docs/development/state.md`](docs/development/state.md);
> plan: [`docs/development/kimiya-port-plan.md`](docs/development/kimiya-port-plan.md).

Pure computation: no I/O, no network, no `unsafe`. Every fallible function returns
`Result<T, KimiyaError>`; NaN, ±∞ and out-of-domain inputs are rejected rather than
propagated.

## Features

| Area | Module | Highlights |
|------|--------|------------|
| Elements | `element` | Full periodic table (118 elements, H–Og): IUPAC atomic weights, Pauling electronegativity, categories, oxidation states; ionization energy, electron affinity and covalent radius for common elements |
| Molecules | `molecule` | Composition, molecular weight, formula strings, presets (water, CO₂, methane, glucose) |
| Stoichiometry | `stoichiometry` | Exact integer equation balancing, oxidation-state assignment, limiting reagent, theoretical/percent yield, percent composition, empirical formula |
| Reactions | `reaction` | Gibbs free energy, equilibrium constant, spontaneity |
| Thermochemistry | `thermo`, `thermochem` | q = mcΔT, Hess's law; ΔHf°/ΔGf°/S° table; multi-segment NIST Shomate Cp(T); van 't Hoff; adiabatic flame temperature |
| Kinetics | `kinetics` | Arrhenius, zero/first/second/nth-order, Michaelis–Menten, Lineweaver–Burk, collision theory, Eyring |
| Reaction dynamics | `reaction_dynamics` | ODE kinetics (RK4, adaptive DOPRI45, implicit stiff, stochastic) via hisab |
| Reaction networks | `network` | Stoichiometric matrices (dense and CSR), rank, conservation laws |
| Gases | `gas` | Ideal gas, Van der Waals, Redlich–Kwong, Peng–Robinson (analytic cubic roots), fugacity |
| Solutions | `solution` | Molarity/molality, pH/pOH, weak acid/base pH (full charge balance incl. Kw), Henderson–Hasselbalch, ionic strength, Debye–Hückel, Davies, Kohlrausch |
| Phase equilibria | `phase` | Raoult, Henry, colligative properties, Clausius–Clapeyron, Antoine, Gibbs phase rule |
| Electrochemistry | `electrochemistry` | 49 standard reduction potentials, Nernst, cell potential, Faraday's laws |
| Spectroscopy | `spectroscopy` | Beer–Lambert, photon energy, Bohr/Rydberg, harmonic oscillator, rigid rotor, isotope patterns, FFT, peak finding |
| Organic | `organic` | Average bond energies, bond-enthalpy ΔH, VSEPR geometry, functional groups |
| Inorganic | `inorganic` | Crystal field theory (CFSE, spin states), spin-only moment, Born–Landé and Born–Haber lattice energy, Shannon ionic radii |
| Nuclear | `nuclear` | Decay law, activity, mass defect, binding energy, semi-empirical mass formula, Q-value, decay chains, isotope table |
| Potentials | `potential` | Lennard-Jones, Morse, Coulomb; distance, bond angle, dihedral |
| Data fitting | `fitting` | Arrhenius parameters, polynomial calibration, Beer–Lambert ε |
| Sensitivity | `sensitivity` | Reverse-mode autodiff gradients of chemical quantities |
| Bridges | `bridge` | Primitive-value conversions for ushma, tanmatra and bijli |
| Visualization | `integration::soorat` | Render-ready structures for soorat (feature `soorat-compat`) |

## Quick Start

```rust
use kimiya::{element, gas, molecule::Molecule, solution};

fn main() -> kimiya::Result<()> {
    let fe = element::lookup_by_symbol("Fe").ok_or_else(|| {
        kimiya::KimiyaError::InvalidElement("Fe".into())
    })?;
    println!("{}: {:.3} g/mol", fe.name, fe.atomic_mass);

    let water = Molecule::water();
    assert!((water.molecular_weight()? - 18.015).abs() < 0.01);

    // 1 mol of ideal gas at 273.15 K and 1 atm occupies ~22.4 L.
    let v = gas::ideal_gas_volume(1.0, gas::STP_TEMPERATURE, gas::STP_PRESSURE)?;
    assert!((v * 1000.0 - 22.414).abs() < 0.01);

    // Weak acid pH solves the full charge balance, so dilute acids stay below 7.
    let ph = solution::weak_acid_ph(1.8e-5, 0.1)?; // 0.1 M acetic acid
    assert!((ph - 2.87).abs() < 0.01);
    Ok(())
}
```

More in [`examples/basic.rs`](examples/basic.rs): `cargo run --example basic`.

## Cargo Features

| Feature | Default | Adds |
|---------|---------|------|
| `logging` | no | `kimiya::logging::init()` — installs a `tracing-subscriber` filtered by `KIMIYA_LOG` |
| `soorat-compat` | no | `kimiya::integration::soorat` visualization data types |

## Dependency Stack

| Crate | Purpose |
|-------|---------|
| [hisab](https://crates.io/crates/hisab) 1.4 (`calc`, `num`, `autodiff`) | Numerical methods: root finding, ODE solvers, least squares, FFT, autodiff, SVD rank |
| serde 1 | `Serialize` on public data types |
| thiserror 2 | `KimiyaError` |
| tracing-subscriber 0.3 (optional) | `logging` feature |

Minimum supported Rust version: **1.99** (edition 2024).

## Consumers

- **ushma** — thermochemistry coupling (heat of reaction)
- **bijli** — electrochemistry (redox, Nernst equation)
- **joshua** — chemistry simulation
- **dravya** — material chemical composition

## Development

```bash
make check   # fmt + clippy + test + audit
make deny    # supply-chain policy
make bench   # ./scripts/bench-history.sh
```

See [CONTRIBUTING.md](CONTRIBUTING.md) and [CLAUDE.md](CLAUDE.md) for the work loop.

## License

GPL-3.0-only
