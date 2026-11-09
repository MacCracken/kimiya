# kimiya — Claude Code Instructions

> **Core rule**: this file is **preferences, process, and procedures** —
> durable rules that change rarely. Volatile state (current version,
> module line counts, port progress, test counts, consumers) lives in
> [`docs/development/state.md`](docs/development/state.md).
> Do not inline state here.

## Project Identity

**Kimiya** (Arabic: كيمياء — alchemy, the root of "chemistry") — Chemistry engine for AGNOS

- **Type**: Flat library (Cyrius), being ported from the Rust crate 1.2.0; the port ships as 2.0.0
- **License**: GPL-3.0-only
- **Language**: Cyrius. The toolchain is pinned in `cyrius.cyml` `[package].cyrius`
- **Version**: `VERSION` at the project root is the source of truth. Don't inline the number here
- **Rust oracle**: kimiya 1.2.0, frozen at `rust-old/`. It stays the parity reference until the deferred
  bridges and integrations are ported or dropped (port plan §2.2)
- **Plan of record**: [`docs/development/kimiya-port-plan.md`](docs/development/kimiya-port-plan.md)
- **Standards**: [First-Party Standards](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-standards.md) · [First-Party Documentation](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-documentation.md)

## Goal

kimiya owns chemistry for AGNOS: elements and molecules, stoichiometry, thermochemistry,
gases, solutions, electrochemistry, kinetics and reaction dynamics, spectroscopy, nuclear
chemistry and intermolecular potentials. Every value comes from an evaluated primary source.

## Consumers

The published Rust dependents (jivanu, khanij, kiran, jnana, tara) use the Rust crate; the
current list and what each calls is in `state.md`. Cyrius consumers: none yet.

## Quick Start

```sh
cyrius deps                                   # resolve dependencies (vendors lib/, incl. hisab)
cyrius build src/main.cyr build/kimiya        # compile
cyrius test                                   # run tests/*.tcyr (parity + reference tests)
cyrius bench                                  # run tests/*.bcyr
cargo test --manifest-path rust-old/Cargo.toml --all-features   # the Rust oracle, read-only
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # golden corpus up to date
python3 -I scripts/oracle/validate.py .                         # independent check of the corpus
```

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read the roadmap, the port plan, the CHANGELOG and the open issues
1. Test + benchmark sweep of the existing code
2. Cleanliness check (below)
3. Get baseline benchmarks
4. Internal deep review
5. External research — CODATA, NIST, IUPAC, CRC, JANAF, AME/NUBASE
6. Cleanliness check — must be clean after the review
7. Additional tests/benchmarks from the findings
8. Post-review benchmarks
9. Repeat if heavy

### Cleanliness check

```sh
cyrius fmt --check
cyrius lint
cyrius vet
cyrius build src/main.cyr build/kimiya
cyrius test                 # must be warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius distlib --all && cyrius distlib --check   # once [lib] exists (M1)
cyrius deps --locked
cyrius audit
```

`cyrius test` must be **warning-free**. The f64 type diagnostics come only from the
compiles that `test` and `build` drive, never from `cyrius lint`. Each one is a bug:
- *integer arithmetic with an f64 right operand*
- *f64 arithmetic with a non-f64 right operand*
- *an integer stored into an f64 slot*
- *comparison mixes f64 and integer*

### Work Loop (continuous)

1. Work phase — one module at a time, in the order the port plan gives
2. Cleanliness check
3. Test + benchmark additions: golden parity rows, reference-value tests, benches
4. Run benchmarks, recorded against the Rust numbers (`rust-old/benchmarks.md`)
5. Internal review, including a surface-parity diff against `rust-old/src/`
6. Cleanliness check
7. Deeper tests/benchmarks
8. Benchmarks again
9. If the review was heavy, return to step 5
10. Documentation — CHANGELOG, `state.md`, roadmap, the port plan's tracking table
11. Version check — `VERSION` and `cyrius.cyml` (`${file:VERSION}`) in sync
12. Return to step 1

### Task Sizing

- **Low/Medium effort**: Batch freely — multiple items per work loop cycle
- **Large effort**: Small bites only — break into sub-tasks, verify each before moving to the next. Never batch large items together
- **If unsure**: Treat it as large. Smaller bites are always safer than overcommitting

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks
- Never refactor speculatively. Wait for the third instance before extracting an abstraction
- Refactoring is part of the work loop, not a separate phase
- Every refactor must pass the same cleanliness + benchmark gates as new code

## Key Principles

- **Cross-check against `rust-old/`.** The port's correctness bar is "matches what Rust
  did", within the tolerance class the port plan (§5) assigns. Diverge only with an ADR
  in `docs/adr/` (ADR 0002 lists the API changes accepted for 2.0).
- **Correctness over cleverness.** If the Cyrius behavior diverges silently from Rust, the bugs win.
- Test after every change, not after the feature is "done". ONE change at a time.
- **Never skip benchmarks.** Numbers don't lie.
- **Own the stack.** Use hisab (Cyrius) for numerics; where hisab lacks a routine the oracle
  used, follow the hisab-gap ADR rather than inventing a new algorithm.
- **Enums.** Public enums have explicit, append-only numeric values that are never
  renumbered (the Cyrius form of `#[non_exhaustive]`). Property functions `match` with an `_ =>` arm.
- **`#must_use` on all pure functions**, plus `#pure` where it applies.
- **No aborts in library code.** `src/` never uses `+?`, `-?` or `*?` (they exit the
  process), `assert_fatal` or `result_unwrap`. Errors return a value-form `Result` or an
  `err_out` code. Every guard rejects NaN and ±∞ (the 1.2 `validate` rule).
- **f64 hygiene.**
  - Bind every call result to a `: f64` local before using it in an operator.
  - Never put hex in a `const`: it is an integer.
  - Never write `(0 - x)` on an f64: use `-x`.
  - No exponent notation in literals: the lexer has none. Generated files write decimals positionally.
  - Use the canonical `ganita_f64_*` names, not the deprecated short aliases.
- **Tables and constants are generated, never typed by hand.** `scripts/oracle/` writes
  `src/constants.cyr` and `src/*_data.cyr` from the oracle's API. Regenerate them; never edit them.
- **Element data must be accurate** — real atomic masses, real electronegativity values.
- **Physical constants must be precise** — R = 8.314462618, N_A = 6.02214076e23 (exact SI values).
- **Every public global carries the `KM_` prefix** (constants, enum members, tables); private
  helpers use `_km_`. A bundle shares one namespace with its consumers.
- **Function names follow the Rust API**, except where ADR 0002 renames them.
- Build with `cyrius build`, not raw `cat file | cc5` — the manifest auto-resolves deps.
- Source files only need project includes — stdlib auto-resolves from `cyrius.cyml`.
- `var buf[N]` = N **bytes**, not N entries.

## DO NOT

- **Do not commit or push** — the user handles all git operations (including staging)
- **NEVER use `gh` CLI** — use `curl` to the GitHub API only
- **Do not modify `rust-old/`** — it is the parity oracle. Tools that drive it (golden and table
  generators) live in `scripts/oracle/` and depend on it by path
- Do not modify `lib/` files (vendored stdlib / dep symlinks)
- Do not hardcode toolchain versions in CI YAML — `cyrius = "X.Y.Z"` in `cyrius.cyml` is the source of truth
- Do not widen a parity tolerance class to make a row pass. Class-C bounds are set by the
  rule in `docs/development/port-parity.md` before a function's first comparison
- Do not add unnecessary dependencies
- Do not break backward compatibility without a major version bump
- Do not skip benchmarks before claiming performance improvements
- Do not skip tests before claiming changes work
- Do not use approximate/made-up chemistry values — always cite sources

## Documentation Structure

```
Root files (required):
  README.md, CHANGELOG.md, CLAUDE.md, CONTRIBUTING.md, SECURITY.md, CODE_OF_CONDUCT.md, LICENSE

docs/:
  development/kimiya-port-plan.md — the port plan of record
  development/state.md            — live state
  development/roadmap.md          — completed items, backlog, port API changes
  development/port-parity.md      — measured parity baselines (M0 step 8)
  adr/NNNN-title.md               — architecture decision records (template.md)
  architecture/                   — non-obvious constraints
  guides/                         — task-oriented how-tos
rust-old/docs/                    — the Rust 1.2.0 architecture overview, kept for reference
```

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/): `## [Unreleased]`, then per version
`### Added / Changed / Fixed / Removed / Security / Performance`.

- Every change that alters behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers (Rust vs Cyrius where relevant)
- Breaking changes get a **Breaking** section with a migration guide
- Group by module when multiple changes land in one release
- Link the ADR when a change was driven by an architectural decision
