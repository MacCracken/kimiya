# kimiya Rust → Cyrius port plan

> **Plan of record** for porting kimiya to Cyrius. Rust **1.2.0** is the final Rust release
> and the **reference implementation** (parity oracle); the Cyrius port ships as **2.0.0**.
> The approach is tanmatra's ([ADR 0006 there](../../../tanmatra/docs/adr/0006-cyrius-port-with-rust-oracle.md),
> [its port plan](../../../tanmatra/docs/development/tanmatra-port-plan.md)), which completed
> M0–M7 on 2026-11-08. Where this plan says "as tanmatra", follow that document.
>
> ▶ **Start on Cyrius 6.7.18**, after 1.2.0 is committed and tagged (M0 step 0).

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | kimiya **1.2.0**: 6,824 code lines in `src/` (excluding tests and comments), 202 public functions, 23 modules | `Cargo.toml`, `VERSION` |
| Rust health, 2026-11-08 on rustc 1.99.1 | ✅ 662 unit + 24 integration + 2 doc tests; fmt, clippy `-D warnings` (all and no features), rustdoc `-D warnings`, `missing_docs`, `cargo audit`, `cargo deny`, MSRV 1.99 — all clean | local run |
| Correctness baseline | 2026-11-08 audit (4 reviews, 100+ findings) repaired in 1.2.0 and pinned by reference-value tests | `CHANGELOG.md` 1.2.0 |
| Cyrius toolchain | **6.7.18** (installed) | `~/.cyrius/current` |
| ganita | **1.2.15**, bundled in the 6.7.18 stdlib | `~/.cyrius/versions/6.7.18/lib/ganita.cyr` |
| hisab | **3.4.1**, Cyrius, pinned to 6.7.18, tag `3.4.1`, bundle `dist/hisab.cyr` | `/home/macro/Repos/hisab` |
| Oracle's math library | Rust **std**: exp, exp2, expm1, ln, log10, pow, sin, cos, acos, atan2 → host **glibc 2.44**; cbrt, floor, trunc, round, powi → Rust's own built-ins. Release build (LLVM turns `2.0.powf` into `exp2`) | [`port-parity.md`](port-parity.md) §3 |
| Published dependents (crates.io) | jivanu, khanij, kiran, jnana, tara — all Rust, all `^1`; none is in Cyrius | crates.io reverse deps |
| Target | kimiya **2.0.0**, Cyrius | this document |

## 0. Decisions already made

- **No Rust 2.0** (2026-11-08). Rust stays at 1.2.0 with its documented behaviour. The
  breaking API changes from the audit are made in the Cyrius API only; §2.3 lists them and
  how parity is checked across each one. They are accepted by ADR 0002 (M0 step 2).
- **tanmatra template.** Rust 1.2.0 is frozen in `rust-old/` and never edited; a generator
  crate in `scripts/oracle/` drives it to produce golden values; a Cyrius harness replays
  them; parity classes A–D (§5).
- **The port starts from a repaired crate.** Every audit finding is fixed and tested in Rust
  first, so a divergence from the oracle is a port defect unless an ADR accepts it.

## 1. Porting models and references

- **tanmatra 2.0.0**, `/home/macro/Repos/tanmatra` — the template. Copy, then adapt:
  `scripts/oracle/` (generator crate, TSV format 1, `validate.py`, `probe_tables.py`),
  `tests/parity.tcyr`, `tests/parity_core.cyr`, `tests/ulp.cyr`, `docs/development/port-parity.md`,
  the CLAUDE.md principle translation (its plan §2 step 2) and the Rust-test ledger script.
- **hisab 3.4.1**, `/home/macro/Repos/hisab` — the dependency. Read `CLAUDE.md`,
  `docs/guides/migration-3.0.md` (Result convention) and the doc comment of each function
  kimiya calls (§6.2); several changed semantics after Rust 1.4.0.
- **ranga 2.3.2** — differential-oracle model (`scripts/oracle/compare.py`, `ACCEPTED` table).
- **prakash 2.15.1** — `cyrius.cyml` `[lib]` layout and `docs/benchmarks-rust-v-cyrius.md`.
- **mela 1.5.2** — `docs/development/rust-parity.md`, the per-test ledger model for M7.
- **Language and toolchain** (`/home/macro/Repos/cyrius`): `docs/guides/cyrius-guide.md`,
  `docs/stdlib-reference.md`, `docs/development/ecosystem-migration.md` (6.7.x entries),
  `cyrius help manifest`.

---

## 2. Scope

### 2.1 In scope for 2.0.0

| Milestone | Rust module | Code lines | Public fns | Tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M1 | `error.rs` | 18 | — | 2 | `src/error.cyr` |
| M1 | *(constants spread over 6 modules)* | — | (consts) | — | `src/constants.cyr` (generated) |
| M1 | `validate.rs` | 15 | (3 private) | 3 | `src/validate.cyr` |
| M1 | `element.rs` | 1,226 | 3 | 22 | `src/element_data.cyr` (generated), `src/element.cyr` |
| M1 | `molecule.rs` | 80 | 8 | 12 | `src/molecule.cyr` |
| M2 | `reaction.rs` | 42 | 3 | 10 | `src/reaction.cyr` |
| M2 | `thermo.rs` | 36 | 4 | 10 | `src/thermo.cyr` |
| M2 | `thermochem.rs` | 1,201 | 11 | 62 | `src/thermochem_data.cyr` (generated), `src/thermochem.cyr` |
| M2 | `phase.rs` | 245 | 12 | 31 | `src/phase_data.cyr` (generated), `src/phase.cyr` |
| M2 | `gas.rs` | 371 | 12 | 39 | `src/gas.cyr` |
| M3 | `solution.rs` | 249 | 17 | 49 | `src/solution.cyr` |
| M3 | `electrochemistry.rs` | 462 | 12 | 35 | `src/electrochemistry_data.cyr` (generated), `src/electrochemistry.cyr` |
| M3 | `kinetics.rs` | 234 | 15 | 47 | `src/kinetics.cyr` |
| M4 | `stoichiometry.rs` | 458 | 7 | 63 | `src/stoichiometry.cyr` |
| M4 | `organic.rs` | 238 | 7 | 21 | `src/organic.cyr` |
| M4 | `inorganic.rs` | 360 | 9 | 30 | `src/inorganic_data.cyr` (generated), `src/inorganic.cyr` |
| M4 | `potential.rs` | 138 | 9 | 33 | `src/potential.cyr` |
| M4 | `nuclear.rs` | 291 | 13 | 36 | `src/nuclear_data.cyr` (generated), `src/nuclear.cyr` |
| M4 | `spectroscopy.rs` | 339 | 27 | 62 | `src/spectroscopy.cyr` |
| M5 | *(hisab 1.4.0 routines, §6.2)* | — | — | — | `src/numerics.cyr` (private `_km_*`) |
| M5 | `reaction_dynamics.rs` | 293 | 8 | 26 | `src/reaction_dynamics.cyr` |
| M5 | `network.rs` | 71 | 5 | 10 | `src/network.cyr` |
| M5 | `fitting.rs` | 103 | 5 | 16 | `src/fitting.cyr` |
| M5 | `sensitivity.rs` | 86 | 5 | 8 | `src/sensitivity.cyr` |

Counts measured 2026-11-08 on 1.2.0. In scope: 627 of the 662 unit tests, plus the
integration and doc tests that touch in-scope modules.

### 2.2 Deferred until a Cyrius consumer exists

| Rust surface | Why deferred | Unblocked when |
|---|---|---|
| `bridge.rs` (7 fns, 22 tests): ushma, tanmatra, bijli bridges | No Cyrius consumer calls them. tanmatra 2.0.0 deferred its own *kimiya* bridges until kimiya ports, so the two bridge sets should port together | tanmatra's bridge phase, or ushma/bijli in Cyrius |
| `integration/soorat.rs` (feature `soorat-compat`, 12 tests) | soorat is Rust (1.0.0); its output is a serde JSON shape | soorat ports |
| `logging.rs` (feature `logging`) | Wraps tracing-subscriber, and kimiya emits no log events. The Cyrius analogue would be sakshi | a consumer asks for it (likely dropped) |

Every deferred function is a thin wrapper over in-scope chemistry, so each is a small,
parity-testable file later. `rust-old/` keeps the code and tests as the reference.

### 2.3 API changes in the Cyrius version (accepted divergences, ADR 0002)

Decided 2026-11-08: these differ from the oracle on purpose. Each row says how the
harness still compares them, so no function leaves parity entirely.

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `thermochem::reaction_enthalpy` / `reaction_gibbs_energy` return kJ/mol | return **J/mol** | oracle × 1000, computed as `(Σ kJ) × 1000.0` on both sides → class B, expected 0 ulp |
| 2 | `sensitivity::*_gradient` bare `Vec`/`f64`, unvalidated; `compute_gradient` → `Result` that is never `Err` | four gradients → `Result` with closed-form derivatives; `compute_gradient` returns the vector directly | finite oracle rows → class D; oracle non-finite ↔ Cyrius `Err` |
| 3 | `nuclear::semi_empirical_binding_energy` → NaN for z > a or a = 0 | `Result` | finite rows class B; oracle NaN ↔ `Err` |
| 4 | `phase::antoine_pressure_*` → NaN at the pole, silent extrapolation | `Result`, `Err` outside `t_min_c..t_max_c` | in-range rows class C; oracle NaN or out-of-range ↔ `Err` |
| 6 | `FunctionalGroup::hydrogen_bonding`, `Geometry::TShapedThreeLonePairs` | `is_hydrogen_bond_donor` (+ new `is_hydrogen_bond_acceptor`), `TShaped` | renamed fns class A; acceptor has reference tests only |
| 7 | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational state (`num`/`den`); integer matrix input | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer" `Err` ↔ `Ok(p/q)`, checked by reference tests (Fe₃O₄ 8/3, C₃H₈ −8/3); balancing class A |
| 5, 10 | bridge `NaN` sentinels; soorat `from_molecule` defaults | `Result` | deferred with their modules (§2.2) |
| 8, 9, 11 | Beer–Lambert intercept, Rydberg nuclear mass, stiff trajectory | **not in 2.0** — additive, demand-gated | — |

Proposed defaults (#1 J/mol, #7 rational type, no extrapolating Antoine variant) are
confirmed when ADR 0002 is written.

### 2.4 Not carried over

- **serde.** `Serialize` derives exist on most data types, but no behaviour depends on them
  and no Cyrius consumer exchanges them. The 3 serde round-trip tests don't port.
- **Rust-only mechanics:** `#[inline]`, `impl Fn` closures (→ function pointers, §6.3),
  `&'static` lifetimes, `Option<&'static T>` lookups (→ index or 0).

---

## 3. M0 — Initialize with `cyrius port`

`cyrius port` initializes the project; don't hand-build the structure. tanmatra found that
the 6.7.18 scaffold overwrites project files (filed as
`cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`);
expect the same repairs here unless a newer pin fixed it.

0. **Freeze a green oracle (Rust, before the move).** *Maintainer action.*
   - Commit 1.2.0 and tag it `1.2.0` (the release workflow requires the tag to equal `VERSION`).
     Publishing to crates.io is optional but makes 1.2.0 the version dependents see.
   - The P(-1) check is already clean on rustc 1.99.1 (header table). Re-run it on the tagged
     commit; nothing in `src/` changes after this point.
   - Record the glibc version (`ldd --version`) — the oracle's transcendental bits depend on it.
1. **Scaffold.**
   ```sh
   cd /home/macro/Repos/kimiya
   cyrius port --dry-run .
   cyrius port .
   ```
   Expect `Cargo.toml`, `Cargo.lock`, `src/`, `tests/`, `benches/`, `examples/`, `Makefile`,
   `codecov.yml`, `deny.toml` and `target/` to move into `rust-old/` (`--dry-run` doesn't list
   them by name), and `src/main.cyr`, `cyrius.cyml`, `tests/kimiya.{tcyr,bcyr,fcyr}`,
   `docs/{adr,architecture,guides,examples}/`, `docs/development/state.md` to be written.
2. **Repair what the scaffold overwrites.**
   - **`CLAUDE.md`:** merge `git show HEAD:CLAUDE.md` with the scaffold's Cyrius rules. Keep
     identity, process and DO NOT list; translate the Rust principles as tanmatra did
     (`#[non_exhaustive]` → append-only enum values; `#[must_use]` → `#must_use`/`#pure`; no
     `unwrap`/`panic!` → no `+?`, `-?`, `*?`, `assert_fatal`, `result_unwrap` in `src/`).
     Replace the stale consumer list with the real one (header table).
   - **Roadmap** restored from git; **`state.md`** rewritten for 1.2.0 → 2.0.0-dev.
   - **Workflows:** keep the scaffold's Cyrius CI, add `test --bounds`, `deps --locked`,
     `audit`; move the Rust workflows to `rust-old/.github/workflows/`. Append `/rust-old/target/`
     and `/build/` to `.gitignore`.
   - **Move Rust-only material** to `rust-old/`: `scripts/bench-history.sh`,
     `scripts/version-bump.sh`, `benchmarks.md`, `bench-history.csv`, `docs/architecture/overview.md`.
   - **ADRs** (kimiya has none yet): `docs/adr/0001-cyrius-port-with-rust-oracle.md` (this
     approach) and `docs/adr/0002-port-api-divergences.md` (§2.3).
   - Check `version` in `cyrius.cyml` is `${file:VERSION}`, not `0.1.0`, and that `bench` is in
     `[test] stdlib`; document the scaffold's stub functions so `cyrius audit` passes.
3. **`cyrius.cyml`.** Shape (confirm keys with `cyrius help manifest`):
   ```toml
   [package]
   name = "kimiya"
   version = "${file:VERSION}"
   license = "GPL-3.0-only"
   language = "cyrius"
   cyrius = "6.7.18"

   [deps]
   stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita",
             "overflow", "tagged", "result", "fnptr"]

   [deps.hisab]
   git = "https://github.com/MacCracken/hisab.git"
   tag = "3.4.1"
   modules = ["dist/hisab.cyr"]

   [test]
   stdlib = ["assert", "bench"]

   [lib]
   modules = ["src/error.cyr", ...]   # dependency order, below
   ```
   - `[lib]` order: error → constants → mathx → validate → element_data → element → molecule →
     reaction → thermo → thermochem_data → thermochem → phase_data → phase → gas → solution →
     electrochemistry_data → electrochemistry → kinetics → stoichiometry → organic →
     inorganic_data → inorganic → potential → nuclear_data → nuclear → spectroscopy → numerics →
     reaction_dynamics → network → fitting → sensitivity.
   - `constants.cyr` comes first because 14 Rust modules reach into other modules for
     constants (`ELEMENTARY_CHARGE`, `COULOMB_CONSTANT`, `BOLTZMANN`, `PLANCK`, …).
   - `[lib] modules = []` is invalid; `src/error.cyr` is the first entry in M1.
   - No `[features]` in 2.0: both Rust features are deferred (§2.2).
4. **Toolchain and math check.** Confirm in `lib/ganita.cyr` and the builtins:
   - present: `ganita_f64_acos`, `_atan2`, `_cbrt`, `_pow`, `_hypot`; builtins `f64_sqrt`,
     `f64_exp`, `f64_ln`, `f64_sin`, `f64_cos`, `f64_floor`, `f64_trunc`, `f64_min`, `f64_max`;
   - **missing** (need local helpers): `exp_m1` (3 uses), `log10` (4), Rust's round-half-away
     (4 — `f64_round` ties to even), `powi` (1), `signum`/`clamp`/`total_cmp`;
   - **probe:** `f64_min`/`f64_max` NaN semantics (Rust ignores a NaN operand).
5. **Golden-value generator** `scripts/oracle/` — copy tanmatra's crate and adapt.
   - `kimiya-oracle` depends on `rust-old/` by path; `Cargo.lock` seeded from
     `rust-old/Cargo.lock` so hisab is exactly 1.4.0. Same TSV format 1
     (`fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note`, f64 as bit patterns), deterministic,
     `--check` mode, `validate.py`.
   - **glibc, not `libm`.** tanmatra's oracle used the pinned `libm` crate; kimiya's uses std,
     so its transcendental bits come from the host glibc. Write the glibc version into each
     TSV header, and make `--check` refuse a different one.
   - Coverage grows per milestone (grids in §4); start with `constants` and the `math` sweep
     (every std f64 method kimiya calls, over special values and kimiya's ranges).
6. **Table and constant generator** (same crate). Tables are generated through the oracle's
   public API, never re-typed:
   - `constants.cyr` (all physical constants); `element_data.cyr` (118 rows × 6 fields + 32
     property rows incl. oxidation states); `thermochem_data.cyr` (76 rows, 34 Shomate segments);
     `electrochemistry_data.cyr` (49 potentials); `phase_data.cyr` (8 Antoine sets, K_b/K_f);
     `inorganic_data.cyr` (30 ionic radii, Madelung, Born); `nuclear_data.cyr` (18 isotopes);
     organic bond energies, isotope abundances, PR parameters and Λ° inline in their modules.
   - Baked typed column arrays, decimal literals written positionally (no exponent notation
     in the lexer), never a hex pattern in a `const` (both from tanmatra §2 step 6, probed on
     6.7.18). **Probe** baked string columns for the formula/symbol keys.
   - `probe_tables.py` compiles every cell and compares it bit for bit with the golden rows.
7. **Parity harness** — copy `tests/parity.tcyr`, `tests/parity_core.cyr`, `tests/ulp.cyr`.
   Per-module dispatchers; outcomes ok / FAILED / measured / pending; exit code = failing
   modules; plus a **mapping table** for the §2.3 divergences (oracle NaN ↔ `Err`, ×1000, …).
8. **Baseline measurements** → `docs/development/port-parity.md`: the `math` sweep through every
   Cyrius counterpart on x86_64 ELF, x86_64 PE (wine) and aarch64 (qemu), against glibc.
   Set class-C bounds with tanmatra's rule (K_F = 2·S·calls + 2; if K_F > 64, make the
   primitive exact by transcribing it).
9. **Decide the hisab gaps** (ADR 0003) — see §6.2. Recommended: transcribe the Rust hisab
   1.4.0 routines kimiya's oracle used into `src/numerics.cyr`, and file a hisab feature
   request for the general versions.

**Gate:** `cyrius build` and `cyrius bench` pass; `cyrius test` runs the harness (0 modules)
with no warnings; `cyrius lint`, `cyrius fmt --check` clean; generators, `tests/ulp.cyr`,
ADRs 0001–0003 and any upstream issues committed.

---

## 4. Milestones

Every milestone ends with tanmatra's gates (its plan §3): `cyrius build`; `cyrius test` and
`test --bounds` **with no warnings** (the f64 type diagnostics only appear there); `cyrius
bench` recorded against the Rust numbers in `rust-old/benchmarks.md`; `vet`, `lint`, `fmt
--check` on the module's files; `distlib --all` + `distlib --check`; `deps --locked`;
**surface parity** (every in-scope Rust `pub fn` has a Cyrius counterpart or a written
reason); CHANGELOG and `state.md`.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation** | error, constants, validate, element, molecule | `KimiyaError` → code enum (`KM_ERR_*`, append-only), no detail string (§6.4); generated constants and element columns; `lookup_by_symbol` string match, `lookup_by_number` O(1); `Molecule` as a struct over an atom vec; NaN/∞ guard helpers | all 118 elements × every field; lookups incl. invalid symbol/Z; molecular weight and formula for the presets and seeded random compositions |
| **M2 Thermochemistry & gases** | reaction, thermo, thermochem, phase, gas | Shomate segments + analytic enthalpy; flame temperature (Newton + bisection — §6.2); van 't Hoff; Antoine (#4); kJ→J (#1); Peng–Robinson analytic cubic (acos/cos/cbrt → class C) and `peng_robinson_z_roots`; Van der Waals, Redlich–Kwong, fugacity | every thermochem row; Cp and H over 298–6000 K per segment, both sides of each boundary; the 5 tested flame mixtures + error cases; PR roots on a T × P grid for the 5 gases, incl. the 3-root region and 90–200 MPa |
| **M3 Solutions & electrochemistry** | solution, electrochemistry, kinetics | Weak acid/base charge-balance Newton (kimiya's own loop — transcribe exactly, class B/D); ionic strength, Debye–Hückel, Davies, Kohlrausch; sorted potential table, Nernst, Faraday; Arrhenius/Eyring/Michaelis–Menten/nth-order (exp, `exp_m1` → class C) | Ka × C grid 1e-40..1e40 × 1e-40..1e20; all 49 couples; kinetics grids over T 1–5000 K and the edge cases the 1.2 tests pin |
| **M4 Structure, nuclear & spectroscopy** | stoichiometry, organic, inorganic, potential, nuclear, spectroscopy | Exact integer balancing (checked i64/i128 → explicit overflow checks; class A) with integer input (#7); rational oxidation states (#7); VSEPR/bond tables; CFSE tables; Born–Landé; SEMF (#3) with `cbrt`; decay chains incl. equal-λ limit (`exp_m1`); Morse/LJ/Coulomb, atan2 dihedral; isotope pattern (pruned convolution, custom `total_cmp` sort); Rydberg/R_H; `find_peaks`; FFT via hisab `num_fft` | balancing corpus (CH₄, KMnO₄/HCl, CuSO₄ both orders, seeded consistent matrices, every error class); oxidation corpus from the 1.2 tests; all 18 isotopes; SEMF Z ≤ A ≤ 3Z for A ≤ 300; isotope patterns up to C₁₆₈/H₈₈; angles at Å and metre scale; FFT n = 2..2²⁰ |
| **M5 Numerics** | numerics, reaction_dynamics, network, fitting, sensitivity | Transcribed solvers (§6.2); Rust closures → function pointer + context pointer (§6.3); trajectories as flat `(t, y[n])` buffers; Jacobian as a flat n×n buffer (no resize hazard); stiff backward Euler (kimiya's own loop); seeded SDE with PCG32 stream 1; rank via SVD or transcribed; sensitivity closed forms (#2); `compute_gradient` on hisab `ad_tape_*` | every `simulate_*` on the mechanism presets at 10, 100, 1000 steps; Robertson stiff problem; seeded SDE runs; fits on exact synthetic data and noisy seeded data; rank of standard and rank-deficient matrices |
| **M7 Release 2.0.0** | all | Rust-test ledger (mela pattern): each of the 627 + integration + doc tests ported, folded into golden rows, or marked n/a / deferred; full surface-parity diff; `docs/benchmarks-rust-v-cyrius.md`; `cyrius distlib --all`; `cyrius audit`; VERSION 2.0.0 | full corpus green on x86_64 ELF, x86_64 PE, aarch64 |

There is no M6: both Rust features are deferred (§2.2). Suggested order: M1 → M2 → M3 → M4
→ M5 → M7. M5 is last because it carries the hisab-gap work and the callback design.

---

## 5. Parity policy

tanmatra's classes and rules apply unchanged (its plan §4):

| Class | Tolerance | kimiya examples |
|---|---|---|
| **A — exact** | bit-identical / equal integers and strings | constants, every table cell, element and isotope data, balancing coefficients, oxidation states, CFSE/unpaired tables, VSEPR, enum properties |
| **B — arithmetic** | ≤ 2 ulp, expected 0 | + − × ÷ √ and `cbrt` only: ideal gas, Van der Waals, molarity, ionic strength, Debye–Hückel limiting law, SEMF, mass defect, Shomate Cp and enthalpy, Born–Landé, LJ, Coulomb |
| **C — transcendental** | ulp ≤ K, K set per function by the port-parity rule before its first comparison | anything through exp/ln/pow/log10/acos/atan2/cos: Arrhenius, Eyring, Nernst, pH, Clausius–Clapeyron, Antoine, PR roots, Morse, decay, Rydberg, isotope patterns |
| **D — iterative** | relative ≤ 1e-12 | Newton/bisection results, ODE trajectories, least squares, autodiff gradients, weak-acid pH where the iteration count differs |

kimiya specifics:
- **The class-C reference is glibc**, measured per target in M0 step 8.
- **hisab-derived values** are class D unless the routine is transcribed from Rust hisab
  1.4.0, in which case they are expected to reach class B (§6.2).
- **Never widen a class to make a row pass.** A failure is a port defect, a Cyrius/hisab
  defect (file upstream), or a Rust oracle defect (fix Rust in 1.2.x and regenerate). The only
  exception is an ADR-accepted divergence, which the harness still reports by name.
- **The reference-value tests travel with the code** (NIST/JANAF Cp, flame temperatures,
  AME2020 binding energies, NUBASE half-lives, CRC potentials, IUPAC weights): parity proves
  "same as Rust", the reference tests prove "right".
- **Branch points on both sides:** Shomate segment joins; PR one-root vs three-root regions;
  the weak-acid upper bounds; equal vs near-equal decay rates; the isotope-pattern pruning
  threshold; oxidation-state peroxide/hydride rules.

---

## 6. Rust → Cyrius mapping (kimiya specifics)

### 6.1 Math functions

Every std f64 method in `src/` (non-test), 1.2.0:

| Rust (calls) | Cyrius | Parity |
|---|---|---|
| `sqrt` (18), `abs` (13), `floor`, `trunc` | `f64_sqrt`, `f64_abs`, `f64_floor`, `f64_trunc` (builtins) | exact |
| `exp` (18), `ln` (7) | `f64_exp`, `f64_ln` (need `math` on aarch64) | class C; x87 on x86_64, so bits differ across targets |
| `powf` (10) | `ganita_f64_pow` | class C; note ganita's `pow(-1, huge)` defect (tanmatra ADR 0007) — kimiya's bases are positive |
| `cos` (4), `sin` (1), `acos` (1), `atan2` (2) | `f64_cos`, `f64_sin`, `ganita_f64_acos`, `ganita_f64_atan2` | class C (PR trig roots, dihedral) |
| `cbrt` (2) | `ganita_f64_cbrt` (correctly rounded) | exact vs a correctly rounded glibc cbrt — confirm in step 8 |
| `log10` (4), `exp_m1` (3) | **local** `_km_log10`, `_km_expm1` | class C, or transcribe if K > 64 |
| `round` (4) | **local** `_km_round_half_away` (`f64_round` ties to even) | exact |
| `powi` (1) | repeated multiply in Rust's order | class B |
| `max` (13), `min` (7), `clamp`, `signum`, `total_cmp` | builtins / local helpers with Rust's NaN semantics (probe) | exact |

### 6.2 hisab: what the oracle used vs what hisab 3.4.1 provides

| Rust hisab 1.4.0 (used by) | hisab 3.4.1 Cyrius | Plan |
|---|---|---|
| `num::newton_raphson`, `num::bisection` (thermochem flame) | `num_newton`, `num_bisection` — semantics changed in 3.x (subnormal derivatives, NaN refusals, exact endpoint roots) | use hisab; rows class D. Transcribe 1.4.0's ~15-line loops only if D fails |
| `num::rk4_trajectory` — **systems** (dynamics) | `num_rk4` — **scalar** single step only | **gap** → transcribe into `numerics.cyr` |
| `num::dopri45` — **systems** (`simulate_adaptive`) | `ode_dopri45`, `ode_dopri45_trajectory` — **scalar** | **gap** → transcribe |
| `num::euler_maruyama`, `num::Pcg32` (`simulate_stochastic`) | `ode_euler_maruyama` (scalar step), `pcg32_new` (stream 1, as Rust) | PCG32 from hisab; system SDE step local (trivial) |
| `num::gaussian_elimination` (`simulate_stiff` Newton step) | none public (iterative solvers only) | **gap** → transcribe (partial pivoting) |
| `num::least_squares_poly` (fitting) | none | **gap** → transcribe (MGS QR, its 1e-12 rank threshold) |
| `num::matrix_rank` (network) | `svd_compute_checked` | count σ > 1e-10 from hisab's SVD (class A on rank) — or transcribe if a row disagrees |
| `CsrMatrix::from_dense` (network) | `csr_from_dense` | use hisab |
| `num::fft`, `Complex` (spectroscopy) | `num_fft`, `HComplex` — bits move from n ≥ 256 vs 3.3.x | use hisab; class C for n ≥ 256 |
| `autodiff::{Tape, Var, reverse_gradient}` (sensitivity) | `ad_tape_new`, `ad_grad`, … (tape + slot indices) | only `compute_gradient` needs it after #2; redesign its callback around the tape API |
| `EPSILON_F64` (network, thermo) | hisab/stdlib epsilon constant | exact |

**Recommendation for ADR 0003:** transcribe the five gap routines from
`~/.cargo/registry/src/*/hisab-1.4.0/src/num/` into `src/numerics.cyr` as private `_km_*`
functions, so M5 reaches class B against the oracle, and file a hisab feature request for
system ODE integrators, a dense solver, polynomial least squares and `matrix_rank`. When hisab
ships them, a later kimiya release can switch and re-measure.

### 6.3 Types and idioms

| Rust | Cyrius |
|---|---|
| `Result<T, KimiyaError>` | value-form `Result` (stdlib `result`) + `KM_ERR_*` code enum, as tanmatra |
| `Option<&'static T>` lookups | table index (or −1) / 0 handle |
| `&'static str` table keys, `&'static [i8]` lists | baked string columns; flattened start/count/member tables |
| `Vec<f64>`, `Vec<(f64, f64)>` returns | stdlib `vec` / caller-provided out buffers |
| `Vec<Vec<f64>>` matrices (balancing, network, Jacobian) | flat row-major buffer + dimensions |
| `impl Fn(f64, &[f64], &mut [f64])` closures | function pointer + context pointer (stdlib `fnptr`/`callback`), as hisab's ODE API |
| checked / saturating / i64-widened integer ops | stdlib `overflow` operators (`+|`, …) or explicit checks; `i128` sites → explicit two-word or bounded checks |
| `validate::is_positive` and friends | the same predicates in `validate.cyr` |
| `#[must_use]` / `#[non_exhaustive]` | `#must_use` / append-only enum values |

### 6.4 Decisions made at the start of M1 (2026-11-09, probed on 6.7.18)

- **Errors** — `enum KmError` with append-only values in Rust variant order: `KM_ERR_NONE = 0`,
  `INVALID_ELEMENT = 1`, `INVALID_REACTION = 2`, `INVALID_CONCENTRATION = 3`,
  `INVALID_TEMPERATURE = 4`, `INVALID_INPUT = 5`, `COMPUTATION_ERROR = 6`. A fallible function
  returns a value-form `Result` whose `Err` payload is the code; **no detail string** travels with
  it (as tanmatra: the Rust payloads only restate the caller's input). `km_error_name(code)` gives
  the Rust variant name (the golden token `E:InvalidElement`), `km_error_message(code)` the fixed text.
- **String tables** — baked arrays cannot hold string literals ("array initializer elements must
  be constants"). Generated string columns are a NUL-separated `u8` blob plus a `u16` offset
  column (`&_KM_X_BYTES + _KM_X_OFF[i]` is a C string; `streq`/`strlen` work, `--bounds` clean).
- **Lookups** — `Option<&'static T>` becomes a handle: the atomic number (or table row + 1), and
  0 for `None`. Field accessors take the handle; an out-of-range handle returns the field's
  documented "none" value, never reads out of bounds.
- **`Option<f64>` fields** (electronegativity) — the stdlib value-form `Option` if it carries an
  f64 payload cleanly; otherwise a `has_` flag column plus the value.
- **Math helpers** — `src/mathx.cyr`, private `_km_*`: `_km_min` / `_km_max` (Rust's NaN-ignoring
  semantics, including its (+0, −0) operand choice), `_km_signum`, `_km_round_half_away`,
  `_km_powi` (compiler-rt `__powidf2`'s square-and-multiply order), `_km_total_cmp`,
  `_km_is_finite`, and `_km_log10` / `_km_expm1` **transcribed from Rust `libm` 0.2.16**
  (`src/math/{log10,expm1}.rs`, fdlibm/musl, ≤ 1 ulp) so they sit within K = 2 of glibc.
- **Structs** — `Molecule` becomes `KmMolecule` with an `impl` block; atoms are a vec of
  (Z, count) pairs; constructors take a flat (Z, count) pair array and its length.

Settled while implementing M1 (2026-11-09), all within the decisions above:
- **Names.** Free functions get the module prefix where the bare Rust name would be ambiguous in
  one namespace: `element_lookup_by_number` / `_by_symbol` / `_properties`. Struct fields of a
  handle are `element_<field>(h)` and `element_properties_<field>(p)`; the oxidation-state slice is
  `element_properties_oxidation_state_count(p)` + `element_properties_oxidation_state(p, k)`, which
  returns an `Option` (a state can be 0 or negative, so there is no sentinel). `Atom`'s fields are
  `m.atom_element_number(k)` / `m.atom_count(k)` (−1 out of range), `atoms.len()` is `m.atoms_len()`.
- **Electronegativity** is the value-form `Option` (it carries an f64 cleanly); the data is a `bool`
  presence column plus the value, because a compile-time f64 can't be NaN.
- **Out-of-type inputs.** A count is stored mod 2^32 (Rust's u32); an atomic number is stored as
  given, so one outside 1–118 (or past u8) is `Err(INVALID_ELEMENT)`, never aliased mod 256.
- **No aborts.** `KmMolecule` uses the allocator-taking `vec_*_a` / `str_builder_*_a` calls; an
  allocation failure leaves a molecule with no atoms (`atoms` 0), and `molecular_weight` / `formula`
  return `Err(COMPUTATION_ERROR)`.
- **Stdlib includes in `src/`.** `element.cyr` includes `lib/string.cyr`, `molecule.cyr`
  `lib/vec.cyr` and `lib/str.cyr`. `cyrius distlib`'s compile-verify only sees calls a compile
  reaches, and a bundle reaches none, so without them `dist/kimiya.deps` listed only `tagged` and a
  consumer's build failed on `str_builder_*` (probed). No `src/` file includes another.

Language traps already probed on 6.7.18 by tanmatra (its plan §0, §5): no exponent literals;
`(0 - x)` on an f64 is wrong; a hex literal in a `const` is an integer; `f64_round` ties to
even; `+?` aborts the process; `cyrius test` must be warning-free.

---

## 7. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Initialize | ✅ done 2026-11-09 | constants 45 + math 18,647 (2.0 MB): 13,340 ok, 0 failed, 5,351 pending (M1 helpers, module-data constants), 1 accepted (ADR 0004) — identical on x86_64 ELF, x86_64 PE, aarch64 | — | scaffold only | Oracle frozen at tag `1.2.0`; scaffold + repairs; hisab 3.4.1; generator + `constants.cyr` (14 `KM_`) + harness; K = 2 ulp for the math primitives ([`port-parity.md`](port-parity.md)); ADRs 0001–0004 |
| M1 Foundation | ✅ done 2026-11-09 | element 975 + molecule 2,472 + constants 14 (of 45) + math 18,647 (now 0 pending): 22,107 ok, 0 failed, 31 pending (module-data constants), 1 accepted (ADR 0004) — identical on x86_64 ELF, x86_64 PE, aarch64; `molecular_weight` 0 ulp; log10 K = 2 (max 1), exp_m1 K = 2 (max 0) | 38 Rust tests ported (error 2, validate 3, element 22, molecule 11; molecule's `serde_roundtrip` doesn't port), plus 4 Cyrius-only groups: 179 assertions | 3 (`lookup_by_symbol` 1.5× Rust, `lookup_by_number` 5.4×, `molecular_weight_glucose` 7.4×; [benchmarks](../benchmarks-rust-v-cyrius.md)) | `error`, `mathx`, `validate`, generated `element_data` + `element`, `molecule`; `[lib]` + `dist/kimiya.cyr`; surface parity 55/55; settled points in §6.4 |
| M2 Thermochemistry & gases | ☐ | | | | |
| M3 Solutions & electrochemistry | ☐ | | | | |
| M4 Structure, nuclear & spectroscopy | ☐ | | | | |
| M5 Numerics | ☐ | | | | ADR 0003 decides the hisab gaps |
| M7 Release 2.0.0 | ☐ | | | | |
