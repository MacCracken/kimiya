# Kimiya Architecture

Flat library crate. Every module is pure computation over `f64`/integer inputs and
static reference tables; there is no I/O, global state, or `unsafe`.

## Module Map

```
kimiya
├── error.rs              — KimiyaError (6 variants), Result<T>
├── validate.rs           — crate-internal NaN/∞-safe guard predicates
│
├── element.rs            — 118 elements, ELEMENT_PROPERTIES, AVOGADRO, GAS_CONSTANT
├── molecule.rs           — Molecule / Atom, molecular weight, formula
├── stoichiometry.rs      — exact integer balancing, oxidation states, yield, composition
├── organic.rs            — bond energies, VSEPR geometry, functional groups
├── inorganic.rs          — crystal field theory, lattice energy, Shannon radii
│
├── reaction.rs           — ΔG = ΔH − TΔS, K = exp(−ΔG/RT)
├── thermo.rs             — calorimetry, Hess's law
├── thermochem.rs         — ΔHf°/ΔGf°/S° table, multi-segment Shomate Cp(T), flame temperature
├── phase.rs              — Raoult, Henry, colligative, Clausius–Clapeyron, Antoine, phase rule
├── gas.rs                — ideal, Van der Waals, Redlich–Kwong, Peng–Robinson, fugacity
├── solution.rs           — concentration, pH, weak acid/base, activity coefficients, conductivity
├── electrochemistry.rs   — standard potentials, Nernst, Faraday
│
├── kinetics.rs           — rate laws, Arrhenius, Michaelis–Menten, collision theory, Eyring
├── reaction_dynamics.rs  — ODE simulation of mechanisms (hisab rk4 / dopri45 / backward Euler / Euler–Maruyama)
├── network.rs            — stoichiometric matrices, rank, conservation laws
├── fitting.rs            — Arrhenius, polynomial, Beer–Lambert least-squares fits
├── sensitivity.rs        — reverse-mode autodiff gradients
│
├── spectroscopy.rs       — Beer–Lambert, photon energy, Bohr/Rydberg, vibrational/rotational levels, isotope patterns, FFT
├── nuclear.rs            — decay, mass defect, binding energy, SEMF, Q-value, isotopes
├── potential.rs          — Lennard-Jones, Morse, Coulomb, molecular geometry
│
├── bridge.rs             — primitive-value bridges for ushma, tanmatra, bijli
├── integration/
│   └── soorat.rs         — visualization structures (feature `soorat-compat`)
└── logging.rs            — tracing-subscriber init (feature `logging`)
```

## Data Flow

```
static tables (element, thermochem, electrochemistry, phase, nuclear, inorganic, organic)
        │ lookup_* → Option<&'static T>
        ▼
domain functions (validate inputs → closed-form or numerical method) → Result<f64 | Vec | struct>
        │                                      │
        │                                      └── hisab::num / calc / autodiff
        ▼
consumers (ushma, bijli, joshua, dravya) · bridge (primitive f64 in/out) · soorat (render data)
```

## Conventions

- **Units are SI unless a name or doc says otherwise.** Exceptions are explicit in the
  name or doc: `thermochem::reaction_enthalpy` and `reaction_gibbs_energy` return kJ/mol,
  Antoine coefficients are in mmHg/°C, bond energies in kJ/mol, ionic radii in pm.
- **Validation:** guards are written `if !is_positive(x) { return Err(..) }` (see
  `validate.rs`) so NaN and ±∞ are rejected; `x <= 0.0` guards let NaN through.
- **Constants:** CODATA 2018 / 2019 SI exact values (R = 8.314462618, N_A = 6.02214076e23,
  F = 96485.33212, k_B = 1.380649e-23, h = 6.62607015e-34).

## hisab Surface Used

Kimiya depends on hisab with only the `num` and `autodiff` features enabled. This is also
the math surface a port must provide (everything else — equation balancing, cubic EOS
roots, weak-acid pH, Shomate enthalpy, backward-Euler stepping — is implemented in kimiya):

| hisab API | Used by |
|-----------|---------|
| `num::newton_raphson`, `num::bisection` | thermochem (adiabatic flame temperature) |
| `num::rk4_trajectory`, `num::dopri45`, `num::euler_maruyama`, `num::Pcg32` | reaction_dynamics |
| `num::gaussian_elimination` | reaction_dynamics (implicit Newton step in `simulate_stiff`) |
| `num::least_squares_poly` | fitting |
| `num::matrix_rank`, `CsrMatrix` | network |
| `num::fft`, `Complex` | spectroscopy |
| `autodiff::{Tape, Var, reverse_gradient}` | sensitivity |
| `EPSILON_F64` | network, thermo |

## Consumers

- **ushma** — thermochemistry coupling (heat of reaction)
- **bijli** — electrochemistry (redox, Nernst equation)
- **joshua** — chemistry simulation in game worlds
- **dravya** — material chemical composition
