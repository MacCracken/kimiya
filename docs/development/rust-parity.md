# Rust-test ledger and surface parity

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan M7, the mela
> pattern). Every test of the frozen Rust oracle (`rust-old/`, kimiya 1.2.0: `#[test]` fns, the integration
> tests and the doc tests) is ported to a Cyrius test, folded into golden rows, or recorded here with its
> reason; every public item of `rust-old/src` is ported, renamed by ADR 0002, n/a or deferred. The golden
> corpus (`tests/golden/`, `tests/parity.tcyr`) checks every ported function on its grid besides; this page is
> about the oracle's own assertions and its API.

## Summary

**688 Rust tests** (662 unit, 24 integration, 2 doc): **652 ported** (25 of them adapted to an ADR 0002 change or the flat API), 0 folded into golden rows, 1 not applicable, 35 deferred with their modules (port plan §2.2).
Of the unit tests, the 627 in scope (port plan §2.1) are 626 ported and 1 n/a. M7 ported the last 18: 16 integration tests and the 2 doc tests, in `tests/integration.tcyr`.

A port's `test_group` is the Rust test's name, and the script checks that the Cyrius file's `main` calls it. A port asserts what Rust asserts, in the port's form (port plan §6.4: `Result::unwrap` is `is_ok` + the payload, a lookup's `Option<&T>` a handle, a Vec the caller's buffer, a slice a pointer and a count, a value Rust writes with an exponent written positionally); the adapted ones are listed below.

| Rust file | Cyrius tests | Rust tests | Ported | Adapted | Folded | N/A | Deferred | Cyrius-only groups |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `src/error.rs` | `tests/error.tcyr` | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `src/validate.rs` | `tests/validate.tcyr` | 3 | 3 | 0 | 0 | 0 | 0 | 1 |
| `src/element.rs` | `tests/element.tcyr` | 22 | 22 | 0 | 0 | 0 | 0 | 1 |
| `src/molecule.rs` | `tests/molecule.tcyr` | 12 | 11 | 0 | 0 | 1 | 0 | 1 |
| `src/reaction.rs` | `tests/reaction.tcyr` | 10 | 10 | 0 | 0 | 0 | 0 | 3 |
| `src/thermo.rs` | `tests/thermo.tcyr` | 10 | 10 | 0 | 0 | 0 | 0 | 3 |
| `src/thermochem.rs` | `tests/thermochem.tcyr` | 62 | 62 | 4 | 0 | 0 | 0 | 9 |
| `src/phase.rs` | `tests/phase.tcyr` | 31 | 31 | 1 | 0 | 0 | 0 | 3 |
| `src/gas.rs` | `tests/gas.tcyr` | 39 | 39 | 0 | 0 | 0 | 0 | 3 |
| `src/solution.rs` | `tests/solution.tcyr` | 49 | 49 | 0 | 0 | 0 | 0 | 8 |
| `src/electrochemistry.rs` | `tests/electrochemistry.tcyr` | 35 | 35 | 0 | 0 | 0 | 0 | 5 |
| `src/kinetics.rs` | `tests/kinetics.tcyr` | 47 | 47 | 0 | 0 | 0 | 0 | 6 |
| `src/stoichiometry.rs` | `tests/stoichiometry.tcyr` | 63 | 63 | 7 | 0 | 0 | 0 | 11 |
| `src/organic.rs` | `tests/organic.tcyr` | 21 | 21 | 5 | 0 | 0 | 0 | 4 |
| `src/inorganic.rs` | `tests/inorganic.tcyr` | 30 | 30 | 0 | 0 | 0 | 0 | 9 |
| `src/potential.rs` | `tests/potential.tcyr` | 33 | 33 | 0 | 0 | 0 | 0 | 10 |
| `src/nuclear.rs` | `tests/nuclear.tcyr` | 36 | 36 | 1 | 0 | 0 | 0 | 10 |
| `src/spectroscopy.rs` | `tests/spectroscopy.tcyr` | 62 | 62 | 1 | 0 | 0 | 0 | 13 |
| `src/reaction_dynamics.rs` | `tests/reaction_dynamics.tcyr` | 26 | 26 | 1 | 0 | 0 | 0 | 5 |
| `src/network.rs` | `tests/network.tcyr` | 10 | 10 | 1 | 0 | 0 | 0 | 4 |
| `src/fitting.rs` | `tests/fitting.tcyr` | 16 | 16 | 0 | 0 | 0 | 0 | 4 |
| `src/sensitivity.rs` | `tests/sensitivity.tcyr` | 8 | 8 | 1 | 0 | 0 | 0 | 3 |
| `src/bridge.rs` | — | 22 | 0 | 0 | 0 | 0 | 22 | — |
| `src/integration/soorat.rs` | — | 12 | 0 | 0 | 0 | 0 | 12 | — |
| `src/logging.rs` | — | 1 | 0 | 0 | 0 | 0 | 1 | — |
| `tests/integration.rs` | `tests/integration.tcyr`, `tests/organic.tcyr`, `tests/spectroscopy.tcyr`, `tests/reaction_dynamics.tcyr`, `tests/fitting.tcyr` | 24 | 24 | 2 | 0 | 0 | 0 | 0 |
| doc tests (`src/lib.rs`, `src/thermochem.rs`) | `tests/integration.tcyr` | 2 | 2 | 1 | 0 | 0 | 0 | — |
| **Total** | | **688** | **652** | **25** | **0** | **1** | **35** | **117** |

The Cyrius-only groups (`cyrius_*`, and element's `out_of_range_handles`, validate's `guard_edges`, error's `codes`) test what no Rust test reaches: the port's conventions (handles, caller buffers, words for Rust's integer types), the ADR 0002 changes and the correctly rounded primitives. `tests/mathx.tcyr` and `tests/numerics.tcyr` hold only such tests (the primitives, and hisab 1.4.0's transcribed routines against their golden rows); `tests/parity.tcyr` replays the corpus and `tests/kimiya.tcyr` smoke-tests the bundle.

## Integration and doc tests

`rust-old/tests/integration.rs` has 24 tests. 16 are in `tests/integration.tcyr` (M7), which includes the bundle `dist/kimiya.cyr` as a consumer does; the other 8 were ported beside their module's tests:

- `beer_lambert_and_transmittance`: `tests/spectroscopy.tcyr`
- `hydrogen_balmer_series_visible`: `tests/spectroscopy.tcyr`
- `bond_energy_methane_combustion`: `tests/organic.tcyr`
- `vsepr_predicts_water_bent`: `tests/organic.tcyr`
- `reaction_dynamics_first_order_matches_analytical`: `tests/reaction_dynamics.tcyr`
- `reaction_dynamics_reversible_equilibrium`: `tests/reaction_dynamics.tcyr`
- `fitting_arrhenius_roundtrip`: `tests/fitting.tcyr`
- `fitting_calibration_curve`: `tests/fitting.tcyr`

The doc tests:

- `rust-old/src/lib.rs`, the README's Quick Start (rust-old/README.md, included as `ReadmeDoctests`): `tests/integration.tcyr` `_t_doc_readme_quick_start`.
- `rust-old/src/thermochem.rs`, `adiabatic_flame_temperature`'s example: `tests/integration.tcyr` `_t_doc_adiabatic_flame_temperature`. Adapted: ADR 0002 #1: `reaction_enthalpy` is J/mol already, so the example's `dh_kj * 1000.0` is the J/mol value as it is.

## Adapted ports

| Rust file | Test | Cyrius file | Adaptation |
|---|---|---|---|
| `src/thermochem.rs` | `combustion_of_methane_enthalpy` | `tests/thermochem.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 |
| `src/thermochem.rs` | `combustion_of_methane_gibbs` | `tests/thermochem.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 |
| `src/thermochem.rs` | `haber_process_enthalpy` | `tests/thermochem.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 |
| `src/thermochem.rs` | `added_standard_state_elements` | `tests/thermochem.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 |
| `src/phase.rs` | `antoine_pole_and_non_finite_give_nan` | `tests/phase.tcyr` | ADR 0002 #4: each case is Err where Rust asserts NaN ("just above the pole" lies outside water's fitted range) |
| `src/stoichiometry.rs` | `balance_non_integer_entry_is_error` | `tests/stoichiometry.tcyr` | ADR 0002 #7: the integer matrix: 1.5, NaN and ∞ have no counterpart (the corpus counts those rows); Err(InvalidInput) for the integers Rust's 2⁵³ bound refuses, i64::MAX and i64::MIN standing in for 1e300 |
| `src/stoichiometry.rs` | `balance_near_integer_entry_accepted` | `tests/stoichiometry.tcyr` | ADR 0002 #7: the integer matrix: 2 + 10⁻¹² and −1 − 10⁻¹² are the integers 2 and −1 Rust rounds them to |
| `src/stoichiometry.rs` | `oxidation_state_superoxide_is_error` | `tests/stoichiometry.tcyr` | ADR 0002 #7: Ok with KO₂'s exact O −1/2 where Rust asserts Err |
| `src/stoichiometry.rs` | `oxidation_state_single_element` | `tests/stoichiometry.tcyr` | ADR 0002 #7: O₂⁻'s state is the exact −1/2 |
| `src/stoichiometry.rs` | `oxidation_state_propane_non_integer_is_error` | `tests/stoichiometry.tcyr` | ADR 0002 #7: Ok with C₃H₈'s exact C −8/3 where Rust asserts Err |
| `src/stoichiometry.rs` | `oxidation_state_magnetite_non_integer_is_error` | `tests/stoichiometry.tcyr` | ADR 0002 #7: Ok with Fe₃O₄'s exact Fe +8/3 where Rust asserts Err |
| `src/stoichiometry.rs` | `oxidation_state_huge_count_does_not_panic` | `tests/stoichiometry.tcyr` | ADR 0002 #7: (O, u32::MAX)×2 at i32::MIN is the exact −2³⁰/(2³² − 1) |
| `src/organic.rs` | `vsepr_clf3_t_shaped` | `tests/organic.tcyr` | ADR 0002 #6: the renamed `KM_GEOMETRY_T_SHAPED` (Rust's `TShapedThreeLonePairs`) |
| `src/organic.rs` | `alcohol_polar_and_hbonds` | `tests/organic.tcyr` | ADR 0002 #6: the renamed `functional_group_is_hydrogen_bond_donor` |
| `src/organic.rs` | `ether_polar_not_hbond_donor` | `tests/organic.tcyr` | ADR 0002 #6: the renamed `functional_group_is_hydrogen_bond_donor` |
| `src/organic.rs` | `acceptor_only_groups_are_not_donors` | `tests/organic.tcyr` | ADR 0002 #6: the renamed `functional_group_is_hydrogen_bond_donor` |
| `src/organic.rs` | `carboxylic_acid_hbonds` | `tests/organic.tcyr` | ADR 0002 #6: the renamed `functional_group_is_hydrogen_bond_donor` |
| `src/nuclear.rs` | `semi_empirical_invalid_is_nan` | `tests/nuclear.tcyr` | ADR 0002 #3: `semi_empirical_binding_energy` is a Result, so Err(InvalidInput) where Rust asserts NaN |
| `src/spectroscopy.rs` | `hydrogen_abundances_iupac_2009` | `tests/spectroscopy.tcyr` | ADR 0002 #23: the exact pairs are the cited table's (¹H 1.007825031898, ²H 2.014101777844), and the Σ = 1 loop gains sulfur (Rust's ³⁶S-less sulfur sums to 0.9999) |
| `src/reaction_dynamics.rs` | `stiff_resized_jacobian_is_error_not_panic` | `tests/reaction_dynamics.tcyr` | the flat API: the flat n × n Jacobian cannot be resized, so the group checks the nearest case, a Jacobian callback that writes nothing (J = 0, no panic either) |
| `src/network.rs` | `rate_ragged_matrix_is_error` | `tests/network.tcyr` | the flat API: a flat matrix has one row length, so the group checks the nearest case, a rate vector of another length than every row |
| `src/sensitivity.rs` | `invalid_input_gives_non_finite_gradients` | `tests/sensitivity.tcyr` | ADR 0002 #2: the port's Err (with its code) where Rust asserts a non-finite entry |
| `tests/integration.rs` | `methane_combustion_cross_module` | `tests/integration.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 (−800 kJ/mol is −800,000 J/mol) |
| `tests/integration.rs` | `thermochem_enthalpy_entropy_gibbs_triangle` | `tests/integration.tcyr` | ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and tolerances are asserted × 1000 (ΔH − T·ΔS in J/mol, Rust's 1 kJ/mol tolerance as 1,000 J/mol) |

## Tests without a port

| Rust file | Test | Status |
|---|---|---|
| `src/molecule.rs` | `serde_roundtrip` | n/a: a serde round-trip; serde is not carried over and no Cyrius consumer exchanges the type (port plan §2.4) |
| `src/bridge.rs` | all 22 | deferred: `bridge.rs` (the ushma, tanmatra and bijli bridges) waits on a Cyrius consumer; tanmatra deferred its kimiya bridges likewise, so the two sets port together (port plan §2.2) |
| `src/integration/soorat.rs` | all 12 | deferred: `integration/soorat.rs` (feature `soorat-compat`) waits on soorat's port; its output is a serde JSON shape (port plan §2.2) |
| `src/logging.rs` | all 1 | deferred: `logging.rs` (feature `logging`) wraps tracing-subscriber and kimiya emits no log events; the Cyrius analogue would be sakshi, if a consumer asks (port plan §2.2) |

No test is folded into golden rows: every in-scope Rust test has its own port, and the corpus checks the functions on their grids besides.

## Surface parity

`python3 -I -B scripts/surface_parity.py .` checks every public item of `rust-old/src` (before `#[cfg(test)]`; every `pub` line is one of the items counted) against `src/`: a ported module's item by the Cyrius name the porting rules give (CLAUDE.md, port plan §6.4), or in another form with its reason, or renamed by ADR 0002 (the new name must exist); the §2.2 modules' items are deferred, and `lib.rs`'s module declarations and prelude re-exports n/a (each re-export must name a public item of its module).

| Rust file | Public items | Ported by name | Ported in another form | Renamed (ADR 0002) | N/A | Deferred (§2.2) |
|---|---:|---:|---:|---:|---:|---:|
| `bridge.rs` | 7 | 0 | 0 | 0 | 0 | 7 |
| `electrochemistry.rs` | 20 | 19 | 1 | 0 | 0 | 0 |
| `element.rs` | 31 | 28 | 3 | 0 | 0 | 0 |
| `error.rs` | 8 | 7 | 1 | 0 | 0 | 0 |
| `fitting.rs` | 10 | 10 | 0 | 0 | 0 | 0 |
| `gas.rs` | 23 | 23 | 0 | 0 | 0 | 0 |
| `inorganic.rs` | 37 | 36 | 1 | 0 | 0 | 0 |
| `integration/mod.rs` | 1 | 0 | 0 | 0 | 0 | 1 |
| `integration/soorat.rs` | 37 | 0 | 0 | 0 | 0 | 37 |
| `kinetics.rs` | 16 | 16 | 0 | 0 | 0 | 0 |
| `lib.rs` | 42 | 0 | 0 | 0 | 39 | 3 |
| `logging.rs` | 1 | 0 | 0 | 0 | 0 | 1 |
| `molecule.rs` | 13 | 13 | 0 | 0 | 0 | 0 |
| `network.rs` | 6 | 5 | 1 | 0 | 0 | 0 |
| `nuclear.rs` | 24 | 23 | 1 | 0 | 0 | 0 |
| `organic.rs` | 83 | 81 | 0 | 2 | 0 | 0 |
| `phase.rs` | 27 | 26 | 1 | 0 | 0 | 0 |
| `potential.rs` | 10 | 10 | 0 | 0 | 0 | 0 |
| `reaction.rs` | 3 | 3 | 0 | 0 | 0 | 0 |
| `reaction_dynamics.rs` | 9 | 8 | 1 | 0 | 0 | 0 |
| `sensitivity.rs` | 5 | 5 | 0 | 0 | 0 | 0 |
| `solution.rs` | 27 | 27 | 0 | 0 | 0 | 0 |
| `spectroscopy.rs` | 37 | 37 | 0 | 0 | 0 | 0 |
| `stoichiometry.rs` | 7 | 7 | 0 | 0 | 0 | 0 |
| `thermo.rs` | 5 | 5 | 0 | 0 | 0 | 0 |
| `thermochem.rs` | 27 | 25 | 2 | 0 | 0 | 0 |
| `validate.rs` | 3 | 3 | 0 | 0 | 0 | 0 |
| **Total** | **519** | **417** | **12** | **2** | **39** | **49** |

Renamed by ADR 0002:

- `organic::Geometry::TShapedThreeLonePairs` → `KM_GEOMETRY_T_SHAPED` (ADR 0002 #6, at its Rust position (7): the Rust name is a misnomer, ClF₃ has two lone pairs).
- `organic::FunctionalGroup::hydrogen_bonding` → `functional_group_is_hydrogen_bond_donor` (ADR 0002 #6, with the new functional_group_is_hydrogen_bond_acceptor beside it).

Ported in another form (the port's conventions, port plan §6.3–§6.4):

| Rust item | Cyrius form |
|---|---|
| `error::type Result` | the stdlib value-form Result, Err payload a KM_ERR_* code (port plan §6.3) |
| `element::static ELEMENTS` | the generated columns KM_ELEMENT_* / _KM_ELEMENT_*_BYTES in src/element_data.cyr (row Z − 1, KM_ELEMENT_COUNT rows), read through element_lookup_by_number and the element_* accessors |
| `element::ElementProperties.oxidation_states` | a slice: element_properties_oxidation_state_count(p) + element_properties_oxidation_state(p, k) (Option) |
| `element::static ELEMENT_PROPERTIES` | the generated columns KM_ELEMENT_PROPERTIES_* in src/element_data.cyr (KM_ELEMENT_PROPERTIES_COUNT rows), read through element_lookup_properties and the element_properties_* accessors |
| `thermochem::static THERMOCHEM_DATA` | the generated columns KM_THERMOCHEM_DELTA_HF_KJ / _DELTA_GF_KJ / _S_STANDARD_J and the formula / name blobs in src/thermochem_data.cyr (KM_THERMOCHEM_COUNT rows), read through handles (row + 1): lookup_thermochem and the thermochem_data_* accessors (port plan §6.4) |
| `thermochem::static SHOMATE_DATA` | the generated columns KM_SHOMATE_A … _T_MAX and the formula blob in src/thermochem_data.cyr (KM_SHOMATE_COUNT rows), read through handles (row + 1): lookup_shomate, lookup_shomate_at, shomate_data_formula, shomate_data_coeffs (port plan §6.4) |
| `phase::static ANTOINE_DATA` | the generated columns KM_ANTOINE_A / _B / _C / _T_MIN_C / _T_MAX_C and the name blob _KM_ANTOINE_NAME_BYTES in src/phase_data.cyr (KM_ANTOINE_COUNT rows), read through handles (row + 1): lookup_antoine, antoine_data_name, antoine_data_coeffs (port plan §6.4) |
| `electrochemistry::static STANDARD_POTENTIALS` | the generated columns KM_STANDARD_POTENTIALS_ELECTRONS / _STANDARD_POTENTIAL and the UTF-8 equation and couple blobs in src/electrochemistry_data.cyr (KM_STANDARD_POTENTIALS_COUNT rows), read through handles (row + 1): lookup_half_reaction and the half_reaction_* accessors (port plan §6.4) |
| `inorganic::static IONIC_RADII` | the generated columns KM_IONIC_RADII_CHARGE / _COORDINATION / _RADIUS_PM and the symbol blob _KM_IONIC_RADII_SYMBOL_BYTES in src/inorganic_data.cyr (KM_IONIC_RADII_COUNT rows), read through handles (row + 1): lookup_ionic_radius and the ionic_radius_* accessors (port plan §6.4) |
| `nuclear::static ISOTOPES` | the generated columns KM_ISOTOPES_Z / _A / _HALF_LIFE_S, the symbol blob and the UTF-8 decay-mode blob in src/nuclear_data.cyr (KM_ISOTOPES_COUNT rows), read through handles (row + 1): lookup_isotope and the isotope_* accessors (port plan §6.4) |
| `reaction_dynamics::type TrajectoryPoint` | a point of the caller's flat trajectory buffer: 1 + n f64 (t, then the n concentrations), point k at out + 8·k·(1 + n) (port plan §4, M5) |
| `network::type ReactionDef` | a 32-byte reaction record (the reactants' pointer and count, the products' pointer and count), each term a 16-byte (species index word, coefficient f64) record (src/network.cyr; port plan §6.4, M5b) |

Derived trait impls (not items; 124 on 28 types):

| Trait | Impls | Status |
|---|---:|---|
| `Debug` | 28 | n/a: Rust formatting for `{:?}`; a Cyrius value is a word, a struct or a handle (8 of them on soorat's types: deferred) |
| `Clone` | 28 | n/a: a Cyrius struct or handle is copied by value (8 of them on soorat's types: deferred) |
| `Copy` | 8 | n/a: a Cyrius enum value or struct is copied by value (2 of them on soorat's types: deferred) |
| `PartialEq` | 15 | n/a: a Cyrius enum value or handle is a word, compared with == (8 of them on soorat's types: deferred) |
| `Eq` | 7 | n/a: as PartialEq |
| `Hash` | 3 | n/a: no Cyrius hashed collection holds the enums |
| `Serialize` | 23 | n/a: serde is not carried over (port plan §2.4) (8 of them on soorat's types: deferred) |
| `Deserialize` | 11 | n/a: serde is not carried over (port plan §2.4) (8 of them on soorat's types: deferred) |
| `Error` | 1 | ported in another form: thiserror's Display text without its payload is km_error_message, the variant km_error_name (src/error.cyr; port plan §6.4) |

