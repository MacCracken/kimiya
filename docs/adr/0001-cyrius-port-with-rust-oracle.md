# 0001 — Port to Cyrius, with Rust 1.2.0 frozen as the parity oracle

**Status**: Accepted
**Date**: 2026-11-08

## Context

AGNOS libraries are moving from Rust to Cyrius (agnosticos `cyrius-lang-migration.md`,
Phase 2B lists kimiya). kimiya's math dependency, hisab, is already Cyrius (3.4.1), and
tanmatra completed the same kind of port on 2026-11-08 by the method in its ADR 0006.

kimiya 1.2.0 is the result of a full-crate audit and repair sweep (CHANGELOG 1.2.0): 100+
findings fixed in Rust and pinned by reference-value tests, 688 tests green on rustc 1.99.1.
A translation needs an oracle it can check value by value, and only a repaired crate can be
one.

## Decision

kimiya 2.0.0 is a Cyrius library. The Rust 1.2.0 crate (tag `1.2.0`, commit `077e2e2`) is
frozen in `rust-old/` as the parity oracle, and nothing in it is edited.

- **Initialization.** The port was started with `cyrius port` on Cyrius 6.7.18.
- **Scope, milestones, parity policy.** As `docs/development/kimiya-port-plan.md`: chemistry
  first; `bridge`, `integration::soorat` and `logging` deferred until a Cyrius consumer exists.
- **Parity classes.** Each golden value has a class: A exact, B arithmetic, C transcendental
  with a bound measured per target, or D iterative.
- **Divergences.** A divergence from the oracle beyond its class is a defect, unless an ADR
  accepts it (ADR 0002 accepts the planned API changes).
- **Oracle mechanics.** Generators live in `scripts/oracle/` and depend on `rust-old/` by
  path; `rust-old/Cargo.lock` is committed so the oracle always builds against hisab 1.4.0.

## Consequences

- **Positive**
  - Correctness is measured, not asserted: every in-scope public function is diffed against
    the repaired Rust crate.
  - The audit's reference-value tests travel with the code.
- **Negative**
  - Two toolchains are needed until the port completes: cargo for the oracle and its
    generators, cyrius for the library.
  - The oracle's transcendental functions come from the host glibc (Rust std), so the golden
    corpus is tied to a glibc version (2.44 at freeze) and class C needs per-target bounds.
  - `rust-old/` stays in the tree until the deferred surfaces are ported or dropped.
- **Neutral**
  - `rust-old/README.md` is the 1.2.0 README, restored there after the scaffold moved the
    crate, because `src/lib.rs` runs it as a doctest and `Cargo.toml` names it. It is the
    tagged file, unchanged.
  - The 1.2.0 Rust CI, scripts and architecture overview are kept under `rust-old/`.

## Alternatives considered

- **Port without an oracle** (translate and re-test by hand): rejected — the 1.2 audit showed
  how easily hand-checked values drift, and tanmatra's corpus caught defects a test suite
  would not have.
- **Keep Rust as the primary implementation and expose it to Cyrius via FFI**: rejected — the
  migration plan's goal is a sovereign Cyrius stack, and hisab is already Cyrius.
