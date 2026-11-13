# Kimiya

**Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

Kimiya **2.0.0** is a Cyrius library: a parity port of the Rust crate **1.2.0** (tag `1.2.0`, the crates.io release),
which stays frozen in [`rust-old/`](rust-old/) as the reference implementation. Every function is checked against it
value by value:
- **the golden corpus**: 389,119 rows over 22 modules, identical on x86_64 Linux, x86_64 Windows PE and aarch64.
  342,485 match the oracle within their tolerance class (almost all bit for bit). The other 46,634 differ on purpose:
  ADR 0002's API changes, and ADR 0006's correctly rounded values where glibc's are not
  ([port parity](docs/development/port-parity.md));
- **the Rust tests**: all 688 are accounted for, 652 of them ported
  ([Rust-test ledger](docs/development/rust-parity.md));
- **the public surface**: each of `rust-old/src`'s 519 public items is ported, renamed by ADR 0002, not applicable or
  deferred ([the same page](docs/development/rust-parity.md#surface-parity)).

Pure computation: no I/O, no network. Every fallible function returns a value-form `Result` whose Err payload is a
`KM_ERR_*` code, and NaN, ±∞ and out-of-domain inputs are rejected rather than propagated. Every transcendental
kimiya's functions call is correctly rounded ([ADR 0006](docs/adr/0006-correct-rounding-is-the-accuracy-target.md)).

The [port plan](docs/development/kimiya-port-plan.md) records the decisions; [`state.md`](docs/development/state.md)
records the current state.

## Modules

| Area | Module (`src/`) | Highlights |
|------|--------|------------|
| Elements | `element` | Full periodic table (118 elements, H–Og): IUPAC/CIAAW atomic weights, Pauling electronegativity, categories, oxidation states; ionization energy, electron affinity and covalent radius for common elements |
| Molecules | `molecule` | Composition, molecular weight, formula strings, presets (water, CO₂, methane, glucose) |
| Stoichiometry | `stoichiometry` | Exact integer equation balancing (multi-word integers, no size limit), exact rational oxidation states, limiting reagent, theoretical/percent yield, percent composition, empirical formula |
| Reactions | `reaction` | Gibbs free energy, equilibrium constant, spontaneity |
| Thermochemistry | `thermo`, `thermochem` | q = mcΔT, Hess's law; ΔHf°/ΔGf°/S° table; multi-segment NIST Shomate Cp(T); van 't Hoff; adiabatic flame temperature |
| Kinetics | `kinetics` | Arrhenius, zero/first/second/nth-order, Michaelis–Menten, Lineweaver–Burk, collision theory, Eyring |
| Reaction dynamics | `reaction_dynamics`, `numerics` | ODE kinetics (RK4, adaptive Dormand–Prince, implicit stiff, Euler–Maruyama) |
| Reaction networks | `network` | Stoichiometric matrices (dense and CSR), rank, conservation laws |
| Gases | `gas` | Ideal gas, Van der Waals, Redlich–Kwong, Peng–Robinson (analytic cubic roots), fugacity |
| Solutions | `solution` | Molarity/molality, pH/pOH, weak acid/base pH (full charge balance incl. Kw), Henderson–Hasselbalch, ionic strength, Debye–Hückel, Davies, Kohlrausch |
| Phase equilibria | `phase` | Raoult, Henry, colligative properties, Clausius–Clapeyron, Antoine (within each set's fitted range), Gibbs phase rule |
| Electrochemistry | `electrochemistry` | 49 standard reduction potentials, Nernst, cell potential, Faraday's laws |
| Spectroscopy | `spectroscopy` | Beer–Lambert, photon energy, Bohr/Rydberg, harmonic oscillator, rigid rotor, isotope patterns (AME2020 masses, TICE 2009 abundances), FFT, peak finding |
| Organic | `organic` | Average bond energies, bond-enthalpy ΔH, VSEPR geometry, functional groups (hydrogen-bond donors and acceptors) |
| Inorganic | `inorganic` | Crystal field theory (exact CFSE, spin states), spin-only moment, Born–Landé and Born–Haber lattice energy, Shannon (1976) ionic radii |
| Nuclear | `nuclear` | Decay law, activity, mass defect, binding energy, semi-empirical mass formula, Q-value, decay chains, isotope table |
| Potentials | `potential` | Lennard-Jones, Morse, Coulomb; distance, bond angle, dihedral |
| Data fitting | `fitting` | Arrhenius parameters, polynomial calibration, Beer–Lambert ε |
| Sensitivity | `sensitivity` | Closed-form gradients of chemical quantities; reverse-mode autodiff of your own function on hisab's tape |
| Errors and helpers | `error`, `constants`, `mathx`, `validate` | `KM_ERR_*` codes; the physical constants; the correctly rounded primitives; the NaN/±∞-rejecting guards |

The tables (`*_data.cyr`, `constants.cyr`, `mathx_data.cyr`) are generated from the oracle (and MPFR) by
[`scripts/oracle/`](scripts/oracle/), never typed by hand. The Rust crate's bridges (`bridge.rs`), its soorat
integration (feature `soorat-compat`) and its `logging` feature wait on a Cyrius consumer (port plan §2.2).

## Using it

Declare kimiya in your `cyrius.cyml`, with hisab beside it (the bundle calls it), and vendor both with `cyrius deps`:

```toml
[deps]
stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "args", "math", "ganita",
    "overflow", "tagged", "result", "fnptr"]

[deps.hisab]
git = "https://github.com/MacCracken/hisab.git"
tag = "3.4.1"
modules = ["dist/hisab.cyr"]

[deps.kimiya]
git = "https://github.com/MacCracken/kimiya.git"
tag = "2.0.0"
modules = ["dist/kimiya.cyr"]
```

```cyrius
include "lib/kimiya.cyr"

fn main(): i64 {
    alloc_init();

    # Elements are handles (0 for an unknown symbol), read with the element_* accessors.
    var fe = element_lookup_by_symbol("Fe");
    var name = element_name(fe);                          # "Iron", a C string
    var mass: f64 = element_atomic_mass(fe);              # 55.845 g/mol
    println(name);

    # A fallible function returns a value-form Result; its Err payload is a KM_ERR_* code.
    var tw, vw = KmMolecule_water().molecular_weight();
    var mw: f64 = vw;                                     # 18.015 g/mol

    # 1 mol of ideal gas at 273.15 K and 1 atm occupies ~22.4 L.
    var tv, vv = ideal_gas_volume(1.0, KM_STP_TEMPERATURE, KM_STP_PRESSURE);
    var m3: f64 = vv;                                     # bind a payload to a typed f64 first
    var liters: f64 = m3 * 1000.0;                        # 22.414 L

    # Weak acid pH solves the full charge balance: 0.1 M acetic acid, Ka = 1.8e-5.
    var tp, vp = weak_acid_ph(0.000018, 0.1);
    if (is_err_result(tp) == 1) {
        println(km_error_message(vp));
        return 1;
    }
    var ph: f64 = vp;                                     # 2.875
    return 0;
}

var r = main();
sys_exit_group(r);
```

This is the Rust README's Quick Start; `tests/integration.tcyr` (`doc_readme_quick_start`) asserts its values.

**The API follows Rust's** (port plan §6.3–§6.4 has the full mapping):
- **Names.** Free functions keep their Rust names (element's lookups are `element_lookup_by_*`). An enum's methods are
  `snake(type)_method(v)`, as in `bond_type_energy_kj(b)`. Structs are `Km<Name>` with an `impl`, as in `KmMolecule`
  and `KmShomate`. Globals and enum members carry `KM_`.
- **Tables.** A Rust `Option<&'static T>` lookup returns a handle (0 for none), read with accessors such as
  `element_atomic_mass(h)` or `half_reaction_standard_potential(h)`.
- **Results and options.** A fallible function returns a value-form `Result` (`var t, v = f(..); is_ok(t)`) whose Err
  payload is a `KM_ERR_*` code; `km_error_message(code)` gives Rust's message without its payload. An optional value
  is a value-form `Option`.
- **Collections.** A slice is a pointer and a count; a Rust `Vec` result fills the caller's buffer and returns the
  count; a closure is a function pointer with a context.

**Breaking changes from Rust 1.2.0.** [ADR 0002](docs/adr/0002-port-api-divergences.md) lists every one, with how
the parity harness still compares it; the [CHANGELOG](CHANGELOG.md)'s 2.0.0 entry has a migration note for each.
The ones a Rust caller notices first: `reaction_enthalpy` and `reaction_gibbs_energy` return **J/mol**; the Antoine
pressures, `semi_empirical_binding_energy` and the sensitivity gradients return a `Result`; oxidation states are exact
rationals and `balance_equation` takes an integer matrix; `FunctionalGroup::hydrogen_bonding` is
`functional_group_is_hydrogen_bond_donor` and `Geometry::TShapedThreeLonePairs` is `KM_GEOMETRY_T_SHAPED`.

## Building

Cyrius (toolchain pinned in `cyrius.cyml`):

```bash
cyrius deps                                   # vendor the stdlib and hisab into lib/
cyrius build src/main.cyr build/kimiya        # compile
cyrius test                                   # tests/*.tcyr: the golden corpus, the ported Rust tests, Cyrius-only tests
cyrius test --win                             # the same on x86_64 PE (wine); --aarch64 under qemu-aarch64
cyrius bench                                  # tests/*.bcyr
cyrius audit                                  # fmt + lint + test + bench
cyrius distlib --all                          # dist/kimiya.cyr
```

The oracle and its checks:

```bash
cargo test --manifest-path rust-old/Cargo.toml --all-features                       # the Rust 1.2.0 oracle, read-only
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # golden corpus and tables current
python3 -I -B scripts/oracle/validate.py .                                          # an independent check of the corpus
python3 -I -B scripts/surface_parity.py .                                           # every Rust item has its Cyrius form
python3 -I -B scripts/rust_test_ledger.py . --check                                 # every Rust test is accounted for
```

[`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md) compares every Rust benchmark with the port.
See [CONTRIBUTING.md](CONTRIBUTING.md) and [CLAUDE.md](CLAUDE.md) for the work loop.

## Consumers

The published Rust dependents (jivanu, khanij, kiran, jnana, tara) use the Rust crate 1.x. No Cyrius consumer yet.

## License

GPL-3.0-only
