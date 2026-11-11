# Kimiya Roadmap

## Status

**v1.2.0** — Refresh + audit release. Rust 1.99 / edition 2024, hisab 1.4, criterion 0.8.
23 modules, 118 elements, 688 tests (662 unit + 24 integration + 2 doc), `#![warn(missing_docs)]` clean.
Public API unchanged from 1.1 apart from additions.

Next: **Cyrius port** — plan of record: [`kimiya-port-plan.md`](kimiya-port-plan.md)
(tanmatra template: Rust 1.2.0 frozen in `rust-old/` as the parity oracle, Cyrius ships 2.0.0).
M0–M3 and M3c are done (M3c, ADR 0006: every transcendental the port calls is correctly rounded, and every class-C row
is gated at 0 ulp against the corrected oracle); M4 (structure, nuclear, spectroscopy) is done: M4b (organic,
inorganic), M4c (nuclear), M4a (correctly rounded sin, cos, atan2, acos), M4d (potential, the trig switch), M4e
(spectroscopy, hisab 1.4.0's FFT transcribed), M4e-2 (potential's geometry with an unbounded exponent, ADR 0002 #20),
M4e-3 (the Bohr model, the Rydberg formula and the reduced mass rounded once, ADR 0002 #21 and #22), M4e-4 (the cited
AME2020 / TICE 2009 isotope table with ³⁶S, ADR 0002 #23) and M4f (stoichiometry: the integer matrix and the rational
oxidation states, ADR 0002 #7); M5 (numerics) next — live status in [`state.md`](state.md).

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

- [ ] `thermochem::reaction_enthalpy` / `reaction_gibbs_energy` return **kJ/mol** while every other energy API uses J/mol
- [ ] `sensitivity::*_gradient` → `Result` (today invalid T, n = 0 or Q ≤ 0 yield non-finite values); `compute_gradient`'s `Result` is never `Err`
- [x] `nuclear::semi_empirical_binding_energy` → `Result` (`NaN` for z > a or a = 0 in Rust 1.2) — done in the Cyrius API (M4c, ADR 0002 #3): `Err(KM_ERR_INVALID_INPUT)`
- [ ] `phase::antoine_pressure_mmhg` / `_pa` → `Result` with range checking (today `NaN` at the pole, silent extrapolation outside the fitted range)
- [ ] `bridge::*` → `Result` instead of `NaN` sentinels; remove `ionization_to_bond_energy_estimate` (physically unfounded heuristic)
- [x] `organic::FunctionalGroup::hydrogen_bonding` → `hydrogen_bond_donor` (+ `hydrogen_bond_acceptor`); `Geometry::TShapedThreeLonePairs` → `TShaped` — done in the Cyrius API (M4b, ADR 0002 #6): `functional_group_is_hydrogen_bond_donor`, the new cited `functional_group_is_hydrogen_bond_acceptor`, `KM_GEOMETRY_T_SHAPED`
- [x] `stoichiometry::assign_oxidation_states` → rational result so average states (Fe₃O₄ +8/3) are representable (today `Err`); `balance_equation` → integer matrix input — done in the Cyrius API (M4f, ADR 0002 #7): (Z, numerator, denominator) records, an average outside i32 Rust's Err; a flat row-major integer matrix with Rust's |entry| ≤ 2⁵³
- [ ] M4f's measured candidates under the inexact-formula policy, for the user's decision (ADR 0002's M4f addition): `percent_yield`, `theoretical_yield_moles`, `limiting_reagent`, `percent_composition`, `empirical_formula` against their exact values; `balance_equation`'s i64 elimination overflowing where a valid u32 solution exists (with the overflow audit)
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
- [ ] after M4: the library-wide audit of intermediate over- and underflow of products and quotients of inputs, each edge
  fixed as #20 (the user's decision; the `_km_ux_*` helpers and M4e-3's `_km_cr_ratio` exist; `state.md` lists the
  first candidates, `reduced_mass` no longer among them)
- [ ] *additive* — `fitting::fit_beer_lambert` → also return the intercept
- [ ] *additive* — `spectroscopy::rydberg_wavelength` → nuclear-mass parameter for the reduced-mass correction
- [ ] `integration::soorat::MolecularStructure::from_molecule` → `Result`
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

## Backlog (demand-gated)

- Collision-frequency helper σ·√(8k_BT/πμ)·N_A for `kinetics::collision_theory_rate`
- Liquid/stable-phase selection for Peng–Robinson fugacity (roots are available via `peng_robinson_z_roots`)
- Non-metal hydride rule (NaBH₄, SiH₄) in oxidation-state assignment
- Build CSR stoichiometric matrices directly from triplets
- Optional CODATA 2022 update for measured constants (mₚ, mₙ, u, R∞)
