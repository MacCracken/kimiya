# Kimiya Roadmap

## Status

**v1.2.0** — Refresh + audit release. Rust 1.99 / edition 2024, hisab 1.4, criterion 0.8.
23 modules, 118 elements, 688 tests (662 unit + 24 integration + 2 doc), `#![warn(missing_docs)]` clean.
Public API unchanged from 1.1 apart from additions.

Next: **Cyrius port** (see below).

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
- [ ] `nuclear::semi_empirical_binding_energy` → `Result` (today `NaN` for z > a or a = 0)
- [ ] `phase::antoine_pressure_mmhg` / `_pa` → `Result` with range checking (today `NaN` at the pole, silent extrapolation outside the fitted range)
- [ ] `bridge::*` → `Result` instead of `NaN` sentinels; remove `ionization_to_bond_energy_estimate` (physically unfounded heuristic)
- [ ] `organic::FunctionalGroup::hydrogen_bonding` → `hydrogen_bond_donor` (+ `hydrogen_bond_acceptor`); `Geometry::TShapedThreeLonePairs` → `TShaped`
- [ ] `stoichiometry::assign_oxidation_states` → rational result so average states (Fe₃O₄ +8/3) are representable (today `Err`); `balance_equation` → integer matrix input
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
- [ ] Na₂CO₃ ΔGf° is inconsistent with ΔHf°/S° by 2.55 kJ/mol in NBS/CRC (pinned as a known exception in tests)

---

## Backlog (demand-gated)

- Collision-frequency helper σ·√(8k_BT/πμ)·N_A for `kinetics::collision_theory_rate`
- Liquid/stable-phase selection for Peng–Robinson fugacity (roots are available via `peng_robinson_z_roots`)
- Non-metal hydride rule (NaBH₄, SiH₄) in oxidation-state assignment
- Build CSR stoichiometric matrices directly from triplets
- Optional CODATA 2022 update for measured constants (mₚ, mₙ, u, R∞)
