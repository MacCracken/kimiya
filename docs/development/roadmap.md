# Kimiya Roadmap

## Status

**v2.0.0** (2026-11-13) — the Cyrius port of Rust 1.2.0, the plan of record
[`kimiya-port-plan.md`](kimiya-port-plan.md) done (M0–M7; tanmatra's template: Rust 1.2.0 frozen in `rust-old/` as the
parity oracle). 22 in-scope modules; the golden corpus 342,485 rows ok and 46,634 accepted (ADR 0002's API changes, ADR
0006's correct rounding), 0 failed, on x86_64 ELF, x86_64 PE and aarch64; the Rust tests and the public surface
accounted for in [`rust-parity.md`](rust-parity.md). Live status in [`state.md`](state.md).

Next, post-release (the Backlog): the library-wide overflow audit, **paused after B10a** by the user's decision of
2026-11-12, resumes over the whole library, M5's modules included; the oracle's unconverging Jacobi rank waits on a
decision; the deferred modules (bridges, soorat, logging; port plan §2.2) wait on a Cyrius consumer.

The Rust crate stays at **v1.2.0** (refresh + audit release; Rust 1.99 / edition 2024, hisab 1.4, criterion 0.8; 23
modules, 118 elements, 688 tests: 662 unit + 24 integration + 2 doc). There is no Rust 2.0.

---

## Completed

### v1.0 – v1.1
- [x] 21 chemistry modules, full periodic table, hisab numerical integration
- [x] **`bridge.rs`** — ushma, tanmatra, bijli primitive-value bridges
- [x] **`integration/soorat.rs`** (`soorat-compat`) — molecular structure, reaction network, spectrum, phase diagram visualization data

### v1.2 — Refresh + audit
- [x] Toolchain/deps: MSRV 1.89 → 1.99, hisab 1.1 → 1.4 (`num` + `autodiff` only), criterion 0.5 → 0.8, CI actions current
- [x] Full-crate audit (4 parallel reviews) and repair sweep — 100+ findings fixed; see CHANGELOG 1.2.0
- [x] NaN/±∞-safe validation everywhere (`validate.rs`)
- [x] No `unwrap`/`expect`/`panic!`/`unreachable!` in library code
- [x] Every public item documented; README example runs as a doctest
- [x] Benchmarks for every repaired hot path

### v2.0.0 — the Cyrius port (2026-11-13)
- [x] M0 initialize: the oracle frozen at tag `1.2.0`, the golden-value generator `scripts/oracle/`, the parity harness
- [x] M1 foundation: error, constants, the math helpers, validate, element, molecule
- [x] M2 thermochemistry and gases: reaction, thermo, thermochem, phase, gas
- [x] M3 solutions and electrochemistry: solution, electrochemistry, kinetics; M3c: every transcendental correctly
  rounded (ADR 0006)
- [x] M4 structure, nuclear and spectroscopy: stoichiometry, organic, inorganic, potential, nuclear, spectroscopy
- [x] The overflow audit's B1–B7 and B10a (ADR 0002 #27–#39, #40's primitives); the rest is post-release (Backlog)
- [x] M5 numerics: reaction_dynamics, network, fitting, sensitivity, on hisab 1.4.0's transcribed routines (ADR 0003)
- [x] M7 release: the Rust-test ledger (688 tests: 652 ported, 1 n/a, 35 deferred), the whole-surface parity (519 items:
  429 ported, 2 renamed by ADR 0002, 39 n/a, 49 deferred), every Rust benchmark beside its Cyrius number, `VERSION` 2.0.0

---

## Cyrius Port Preparation

Kimiya is listed in the AGNOS Cyrius migration plan (Phase 2B, pure-math science crates).
The milestones (M0 initialize → M1 foundation → M2 thermochemistry & gases → M3 solutions &
electrochemistry → M4 structure, nuclear & spectroscopy → M5 numerics → M7 release 2.0.0) and
the parity policy are in [`kimiya-port-plan.md`](kimiya-port-plan.md).

Port-relevant facts after v1.2:
- **hisab surface** is limited to `num` + `autodiff`; see `docs/architecture/overview.md` for
  the exact functions. Balancing, cubic EOS roots, weak-acid pH, Shomate enthalpy and the
  backward-Euler step are now self-contained kimiya code, not hisab calls.
- **Integer semantics are explicit**: overflow-prone arithmetic uses checked/saturating ops or
  i64/i128 widening, so the port does not inherit Rust debug-panic / release-wrap behaviour.
- **Static tables** (`&'static str` keys, `&'static [T]`) port as const arrays; lookups are
  linear string searches (except `element::lookup_by_number`, O(1)).
- **Closures into solvers** remain in `reaction_dynamics` (rk4/dopri45/Euler–Maruyama),
  `thermochem::adiabatic_flame_temperature` (Newton/bisection) and `sensitivity` (autodiff tape).
- **serde** is derive-only; no behaviour depends on it. soorat structs define a JSON shape.

### API changes for the Cyrius port

**Decided 2026-11-08:** there will be no Rust 2.0. The Rust crate stays on 1.x with the
behaviour documented in 1.2; the changes below are made in the **Cyrius API only**, so the
two APIs intentionally differ on these items. No published dependent (jivanu, khanij, kiran,
jnana, tara) calls any of them. Proposed defaults: energies in J/mol, a small exact rational
type for average oxidation states, no separate extrapolating Antoine function. Items marked
*additive* are demand-gated and can land in either implementation without breaking anything.

- [x] `thermochem::reaction_enthalpy` / `reaction_gibbs_energy` return **kJ/mol** while every other energy API uses J/mol — done in the Cyrius API (M2b, ADR 0002 #1): J/mol
- [x] `sensitivity::*_gradient` → `Result` (today invalid T, n = 0 or Q ≤ 0 yield non-finite values); `compute_gradient`'s `Result` is never `Err` — done in the Cyrius API (M5b, ADR 0002 #2): validated closed-form gradients; `compute_gradient`'s Err only for a caller's tape's failures
- [x] `nuclear::semi_empirical_binding_energy` → `Result` (`NaN` for z > a or a = 0 in Rust 1.2) — done in the Cyrius API (M4c, ADR 0002 #3): `Err(KM_ERR_INVALID_INPUT)`
- [x] `phase::antoine_pressure_mmhg` / `_pa` → `Result` with range checking (today `NaN` at the pole, silent extrapolation outside the fitted range) — done in the Cyrius API (M2a, ADR 0002 #4): Err outside `t_min_c..t_max_c`
- [ ] `bridge::*` → `Result` instead of `NaN` sentinels; remove `ionization_to_bond_energy_estimate` (physically unfounded heuristic) — deferred with `bridge.rs` until a Cyrius consumer exists (port plan §2.2; ADR 0002 #5)
- [x] `organic::FunctionalGroup::hydrogen_bonding` → `hydrogen_bond_donor` (+ `hydrogen_bond_acceptor`); `Geometry::TShapedThreeLonePairs` → `TShaped` — done in the Cyrius API (M4b, ADR 0002 #6): `functional_group_is_hydrogen_bond_donor`, the new cited `functional_group_is_hydrogen_bond_acceptor`, `KM_GEOMETRY_T_SHAPED`
- [x] `stoichiometry::assign_oxidation_states` → rational result so average states (Fe₃O₄ +8/3) are representable (today `Err`); `balance_equation` → integer matrix input — done in the Cyrius API (M4f, ADR 0002 #7): (Z, numerator, denominator) records, an average outside i32 Rust's Err; a flat row-major integer matrix with Rust's |entry| ≤ 2⁵³
- [x] M4f's measured candidates under the inexact-formula policy: `percent_yield`, `theoretical_yield_moles`, `percent_composition` → their exact rationals rounded once; `limiting_reagent` → the exact comparison of the ratios; `empirical_formula` → exact ratios and an exact 1/10 tolerance — done in the Cyrius API (M4f-2, ADR 0002 #24, #25 and #26, the user's decisions of 2026-11-10)
- [x] `balance_equation`'s i64 elimination overflowing where a valid u32 solution exists (M4f's measurement; with the overflow audit) — done in the Cyrius API (the audit's B7, ADR 0002 #39): Rust's elimination on unbounded integers
- [x] `inorganic::cfse_octahedral` → the exact (3·e_g − 2·t₂g)/5 instead of `t2g·(−0.4) + eg·0.6` (which misses it on 13 of
  22 cases, d⁵ high-spin and d¹⁰ ≠ 0) — done in the Cyrius API (M4b, ADR 0002 #12, the user's decision)
- [x] `nuclear::atoms_after_half_lives`, `is_secular_equilibrium`, `inorganic::delta_tet_from_oct` → the exact N₀/2ⁿ (correctly
  rounded since M4d), λ_p < λ_d/100 and 4Δ/9 — done in the Cyrius API (M4c, ADR 0002 #13–#15)
- [x] `potential::lj_equilibrium_distance` → the correctly rounded σ·2^(1/6); `lennard_jones` / `morse_width_parameter` → no
  spurious overflow of 4·ε / 2·D_e — done in the Cyrius API (M4d, ADR 0002 #16–#18)
- [x] `spectroscopy::wavenumber_from_wavelength` / `wavelength_from_wavenumber` → the correctly rounded 1/(100·x) instead of
  `0.01 / x` (one ulp off on about 14 % of inputs) — done in the Cyrius API (M4e, ADR 0002 #19)
- [x] `potential::distance_3d` / `bond_angle` / `dihedral_angle` → no over- or underflow of the squared coordinates: Rust's
  sequence evaluated with an unbounded exponent, the final √ / atan2 rounded once, Rust's bits wherever its chain stays
  normal, the angles exactly scale-invariant — done in the Cyrius API (M4e-2, ADR 0002 #20, the user's decision)
- [x] `spectroscopy::bohr_energy_level` / `bohr_transition_energy` / `rydberg_wavelength` / `SpectralSeries::wavelength` →
  the exact rational rounded once (−R exactly at n = Z; no cancellation of 1/n₁² − 1/n₂² or E₂ − E₁ for close large n);
  `reduced_mass` → the correctly rounded m₁·m₂/(m₁ + m₂) for every pair (m/2 exactly; no over- or underflow) — done in the
  Cyrius API (M4e-3, ADR 0002 #21 and #22, the user's decisions)
- The library-wide overflow audit (after M4; B1–B7 and B10a done, ADR 0002 #27–#39 and #40's primitives) and the
  oracle's unconverging Jacobi rank (M5b) are post-release work: see the Backlog's *Post-release* list
- [ ] *additive* — `fitting::fit_beer_lambert` → also return the intercept
- [ ] *additive* — `spectroscopy::rydberg_wavelength` → nuclear-mass parameter for the reduced-mass correction
- [ ] `integration::soorat::MolecularStructure::from_molecule` → `Result` — deferred with soorat (port plan §2.2; ADR 0002 #10)
- [ ] *additive* — `reaction_dynamics::simulate_stiff` → return the full trajectory (today endpoints only)

### Data to re-check against a primary source

- [ ] Standard potentials: Sn⁴⁺/Sn²⁺ 0.154 (CRC 0.151), ClO⁻/Cl⁻ 0.890 (CRC 0.81), Sn²⁺/Sn −0.136 (CRC −0.1375), Cr₂O₇²⁻/Cr³⁺ 1.33 (editions list 1.36 or 1.232)
- [ ] Shomate HF(g) — known defect: Cp below 7R/2 (needs a sourced NIST/JANAF fit)
- [ ] Shomate C₂H₆ (+3.4 % at 298 K) and C₃H₈ (+8.4 %) — no NIST Shomate set exists; replace with a sourced fit (NASA-7 / Burcat)
- [ ] Shomate high-T segments for NH₃ (> 1400 K), H₂S, SO₂, etc. — today out-of-range use returns `Err`
- [ ] Antoine acetone range (−20–77 °C) unconfirmed
- [ ] `organic::BondType::energy_kj` — the 45 average bond energies (kJ/mol) carry no source citation in 1.2.0; tie them to a
  primary compilation (e.g. the CRC Handbook's or Luo's bond dissociation energies) and record the source (M4b review)
- [x] `organic` hydrogen-bond acceptor classification (Cyrius 2.0, ADR 0002 #6) — decided by the user (2026-11-09): IUPAC
  2011's criterion with the weak acceptors (π systems, thiol S, C–X halogens), every group but the alkane, each family
  cited
- [ ] IUPAC 2011 hydrogen-bond definition — check the acceptor sentence quoted in `src/organic_data.cyr` against the
  published *Pure Appl. Chem.* 83, 1637 PDF (quoted from a secondary copy; the PDF was unreachable in M4b)
- [ ] Na₂CO₃ ΔGf° is inconsistent with ΔHf°/S° by 2.55 kJ/mol in NBS/CRC (pinned as a known exception in tests)
- [x] `spectroscopy::isotope_abundances` — the masses were AME2020's to five decimals (¹H 1.00783 for 1.007825032, ¹³C
  13.00335 for 13.003354835; within 5·10⁻⁶ u), and sulfur omitted ³⁶S (0.0001; its abundances summed to 0.9999). Done in the
  Cyrius port (M4e-4, ADR 0002 #23, the user's decision of 2026-11-10): the full AME2020 masses (Wang et al., Chin. Phys. C
  45, 030003 (2021)) and the TICE 2009 abundances (Berglund & Wieser, Pure Appl. Chem. 83, 397 (2011)) with ³⁶S, from the
  cited records in `scripts/oracle/isotopes.tsv`; Rust 1.2.0 keeps its table (no Rust release)

---

## Backlog

### Post-release (after 2.0.0)

Paused or deferred by the user's decisions; each is a Cyrius API change for ADR 0002 when it lands.

- [ ] after M4: the library-wide audit of intermediate over- and underflow of products and quotients of inputs (the user's
  decisions of 2026-11-09, "Library-wide audit", and 2026-11-10, "Round once": a single product or quotient rounds once,
  a sum inside gets #20's edge fix; one ADR 0002 row per family; `state.md` has the survey's sixteen families and the
  batches B1 … B12)
  - [x] B1 — F1 (thermo's calorimetry, phase's colligative products: `heat_transfer`, `specific_heat_from_calorimetry`,
    `boiling_point_elevation`, `freezing_point_depression`, `osmotic_pressure` rounded once) and F4 (kinetics:
    `zero_order_half_life`, `steady_state_intermediate`, `pre_equilibrium_rate` rounded once; `second_order_concentration`,
    `michaelis_menten` with an unbounded exponent) — done in the Cyrius API (ADR 0002 #27, #28; M4f-2's `_km_cr_mul_div`
    generalised to `src/mathx.cyr`'s correctly rounded monomial `_km_cr_mono`)
  - [x] B2 — F3, solution and electrochemistry (`dilution`, `solution_conductivity`, Faraday's three laws rounded once;
    `nernst_potential` with an unbounded exponent, the correctly rounded ln Q its primitive) — done in the Cyrius API
    (ADR 0002 #29)
  - [x] B3 — F5, spectroscopy products (`absorbance`, `concentration_from_absorbance`, `rigid_rotor_energy`,
    `rotational_constant`, `diatomic_moment_of_inertia` rounded once; `vibrational_frequency` the correctly rounded √ of
    the exact k/μ, `src/mathx.cyr`'s new `_km_cr_sqrt_div`) — done with B2, in one batch (ADR 0002 #30)
  - [x] B4 — F6, potentials, Born–Landé, binding energy per nucleon (`coulomb`, `coulomb_force` rounded once;
    `morse_width_parameter` the correctly rounded ω·√(μ/(2·D_e)), `src/mathx.cyr`'s new `_km_cr_mul_sqrt_ratio`;
    `lennard_jones`, `born_lande_lattice_energy`, `binding_energy_per_nucleon` with an unbounded exponent; #17 and #18
    folded in) — done in the Cyrius API (ADR 0002 #32), with B5 in one batch, the audit's B3
  - [x] B5 — F2, gas laws (`ideal_gas_pressure`, `ideal_gas_volume`, `gas_density` rounded once; van der Waals,
    Peng–Robinson and Redlich–Kwong pressure with an unbounded exponent, the co-volume tests on the unbounded co-volume,
    Tc^2.5 exactly rescaled) — done in the Cyrius API (ADR 0002 #31), in the audit's B3
  - [x] B6 — F7 and F8, sums with cancellation and the Shomate polynomial (`hess_law`, `enthalpy_from_formation`, the
    three thermochem reaction sums, `q_value`, `ionic_strength`, `born_haber_lattice_energy`, `gibbs_free_energy`,
    `kohlrausch_molar_conductivity`, `limiting_molar_conductivity`, `zero_order_concentration`, `Shomate::cp`,
    `enthalpy_antiderivative` with an unbounded exponent in Rust's summation order; the data-backed Shomate callers
    unchanged) — done in the Cyrius API (ADR 0002 #33, #34), the audit's B4, with #29's fast-path correction and the
    bench split (`tests/audit.bcyr`)
  - [x] B7 — F9, exp and exp10 arguments (`equilibrium_constant`'s −ΔG/(R·T), the three activity models' log₁₀ γ with an
    unbounded exponent, the argument's final operation rounded once, then the unchanged correctly rounded exp / 10^x;
    no new primitive) — done in the Cyrius API (ADR 0002 #35), the audit's B5
  - [x] B8 — the FFT (both the magnitudes' re² + im² and the butterflies, hisab 1.4.0's transform without the range
    limit, by the user's answer "Both"), the flame's exact bracket sign test, `find_peaks`' threshold evaluated exactly
    — done in the Cyrius API (ADR 0002 #36, #37 and #38's flame part; `src/mathx.cyr`'s new `_km_sign_scaled`), the
    audit's B6
  - [x] B9 — F15, an exact `balance_equation` (Rust's elimination on multi-word integers, no size limit) — done in the
    Cyrius API (ADR 0002 #39; the signed multi-word `_km_sb_*` with Lehmer's gcd, the tables from a proven Hadamard
    bound), the audit's B7
  - **Paused after B10a** (the user's decision of 2026-11-12): the rest of the audit is post-port work; it resumes after
    M5 / M7 over the whole library, M5's modules included (their candidate spots are listed under M5 below)
  - [ ] B10 — F10's unbounded exponentials and the exponential × prefactor functions
    - [x] B10a — the primitives: e^x, e^x − 1 and 10^x correctly rounded to 53 bits with an unbounded exponent
      (`_km_cr_exp_ux`, `_km_cr_expm1_ux`, `_km_cr_exp10_ux`), exact up to MPFR's widest exponent range and verified
      against it on every target — done (ADR 0002's "Addition, 2026-11-12 (audit B10a)", for #40), the audit's B8
    - [ ] #40's functions in three sub-batches: kinetics (Arrhenius, collision, Eyring ×2, first-order); nuclear
      (`atoms_remaining`, the decay chain); van 't Hoff, Clausius–Clapeyron, Antoine, Morse and `fugacity_pr`
  - [ ] B11 — F11 and F12 (`nth_order_half_life`'s pow, `henderson_hasselbalch`'s log10)
  - [ ] B12 — the iterative functions, by the user's answer "Fix all three": the weak acid/base pH Newton with an
    unbounded exponent on a slow path only; Peng–Robinson's coefficient edges only (Rust's cubic unchanged); the flame's
    exact sign test was done in B8 (ADR 0002 #38, which B12 extends with the other two)
  - [ ] M5 (for the resumed audit; M5a ported at the plan's parity classes without the rework, listed here, not fixed):
    - `simulate_michaelis_menten`'s rate Vmax·[S]/(Km + [S]) over- or underflows in Vmax·[S] or Km + [S] where the rate
      is representable (kinetics' `michaelis_menten` is #28's)
    - `simulate_second_order`'s −k·[A]·[A] and the consecutive / reversible closures' k₁[A] − k₂[B] (a product chain, and
      sums with cancellation, as #28 and #33)
    - RK4's increment h/6·(k₁ + 2k₂ + 2k₃ + k₄) and stage points y + ½h·k: the rate sum overflows where h/6 times it is
      representable (huge rates, tiny steps)
    - dopri45's stage combinations Σ aᵢⱼ·kⱼ and its error terms before the × hd (the error norm's squares and the
      tol·|y| scale need nothing: they over- or underflow only far beyond where the step factor clamps, err below
      1.9·10⁻⁴ or above 1,845, and the accept test err ≤ 1)
    - `simulate_stiff`'s Newton tolerances 10⁻¹⁰·‖y₊‖∞ and 10⁻⁶·‖y₊‖∞ underflow for ‖y₊‖∞ below about 10⁻²⁹⁸ (a state
      that small then needs an exact zero correction to converge), and its h·J and h·f products
    - Euler–Maruyama's a·dt + b·dW, and the Gaussian elimination's factor·a products (hisab 1.4.0's, as transcribed)
    - M5b (ported the same way, listed, not fixed): the least squares' Vandermonde powers x^j and the QR's unscaled dot
      products and norms Σ q² over- or underflow (NaN coefficients, or the 10⁻¹² test failing, where the fit of a
      rescaled x is representable), and fit_arrhenius' 1/T of a subnormal T is ∞
    - the Jacobi rank's sums app, aqq, apq, its app·aqq and 1 + τ² (τ past 1.3·10¹⁵⁴) and the singular values' Σ b²
      over- or underflow (a huge matrix counts +∞ singular values, a tiny one none: the 10⁻¹⁰ threshold is absolute)
    - `rate_from_stoichiometry`'s Σ S·v (a sum with cancellation, #33's kind) and `stoichiometric_matrix`'s repeated
      coefficient sums
    - the gradients' K/(R·T), (K·(−ΔG))/((R·T)·(R·T)), A·e and −RT/(nFQ): (R·T)² or a product over- or underflows where
      the derivative is representable (the port's ComputationError, or a 0)
- [ ] *oracle defect (M5b), for a decision after the port* — `network::reaction_network_rank`: hisab 1.4.0's one-sided
  Jacobi tests |a_pq| ≤ EPSILON_F64²·√(a_pp·a_qq) (10⁻²⁴ relative), below the rounding of a rotated pair's dot product,
  so most dense matrices end in NoConvergence after 100·max(m, n) sweeps and kimiya 1.2.0's rank is Err(ComputationError)
  for Chapman's ozone cycle, Robertson's kinetics, the H₂ + Br₂ chain and 159 of the 563 golden matrices (all finite ones). The port
  reproduces it (ADR 0003's M5b amendment); a fix (a reachable test, Rust's answer wherever it converges) is a Cyrius
  API change for ADR 0002

### Demand-gated

- Collision-frequency helper σ·√(8k_BT/πμ)·N_A for `kinetics::collision_theory_rate`
- Liquid/stable-phase selection for Peng–Robinson fugacity (roots are available via `peng_robinson_z_roots`)
- Non-metal hydride rule (NaBH₄, SiH₄) in oxidation-state assignment
- Build CSR stoichiometric matrices directly from triplets
- Optional CODATA 2022 update for measured constants (mₚ, mₙ, u, R∞)
