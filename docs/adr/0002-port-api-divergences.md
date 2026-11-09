# 0002 — Accept the planned API changes in the Cyrius 2.0 API

**Status**: Accepted
**Date**: 2026-11-08

## Context

The 2026-11-08 audit found eleven defects whose correct fix changes a public signature. kimiya
1.2.0 documents them but keeps the 1.x API. On 2026-11-08 the maintainer decided there will be
no Rust 2.0: the changes are made in the Cyrius API only. No published dependent (jivanu,
khanij, kiran, jnana, tara) calls any of these functions.

Under ADR 0001 a divergence from the oracle is a defect unless an ADR accepts it. This ADR
accepts these, and fixes how the parity harness still compares each one, so that no function
leaves parity entirely.

## Decision

The Cyrius 2.0.0 API differs from Rust 1.2.0 in exactly these places. Each row's golden rows
are compared through the stated mapping:

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `thermochem::reaction_enthalpy`, `reaction_gibbs_energy` return kJ/mol | return **J/mol** | oracle × 1000, computed as `(Σ kJ) × 1000.0` on both sides — class B, expected 0 ulp |
| 2 | `sensitivity::{equilibrium_constant_gradient, arrhenius_gradient, nernst_gradient, gibbs_temperature_sensitivity}` unvalidated; `compute_gradient` → `Result` never `Err` | the four return `Result`, computed from closed-form derivatives; `compute_gradient` returns its vector directly | finite oracle rows class D; oracle non-finite ↔ Cyrius `Err` |
| 3 | `nuclear::semi_empirical_binding_energy` → NaN for z > a or a = 0 | `Result` | finite rows class B; oracle NaN ↔ `Err` |
| 4 | `phase::antoine_pressure_mmhg`, `antoine_pressure_pa` → NaN at the pole, silent extrapolation | `Result`; `Err` outside the coefficient set's `t_min_c..t_max_c`. No extrapolating variant | in-range rows class C; oracle NaN or out-of-range input ↔ `Err` |
| 6 | `FunctionalGroup::hydrogen_bonding`; `Geometry::TShapedThreeLonePairs` | `functional_group_is_hydrogen_bond_donor` and a new `functional_group_is_hydrogen_bond_acceptor`; `KM_GEOMETRY_T_SHAPED` | donor and geometry rows class A under the new names; acceptor checked by reference tests only |
| 7 | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational oxidation state (numerator / positive denominator, reduced); balancing takes an integer matrix | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer average oxidation state" `Err` ↔ `Ok(p/q)`, verified by reference tests (Fe₃O₄ +8/3, C₃H₈ −8/3); balancing class A after converting the oracle's integer-valued f64 inputs |

Deferred with their modules (port plan §2.2), and decided when they port: #5 `bridge::*` →
`Result` and removal of `ionization_to_bond_energy_estimate`; #10 soorat `from_molecule` →
`Result`.

Not in 2.0 (additive, demand-gated, may land in either implementation): #8 Beer–Lambert
intercept, #9 Rydberg nuclear-mass parameter, #11 full `simulate_stiff` trajectory.

The parity harness implements each mapping in one table and reports these rows by name.

## Consequences

- **Positive** — the Cyrius API is internally consistent (SI energies, `Result` for every
  fallible function, chemically meaningful names) from its first release.
- **Negative** — the Rust and Cyrius APIs intentionally differ on these functions; a consumer
  moving from Rust to Cyrius must apply them (2.0.0's CHANGELOG carries the migration notes).
- **Neutral** — the mapped comparisons are slightly more harness code than 1:1 rows.

## Alternatives considered

- **Rust 2.0 first, then a 1:1 port**: offered 2026-11-08, declined by the maintainer.
- **Port the 1.2 API unchanged**: rejected — it would carry known defects (kJ/J mix, NaN
  sentinels, misleading names) into a new codebase.
