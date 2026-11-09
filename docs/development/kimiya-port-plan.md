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
     reaction → thermo → thermochem_data → thermochem → phase_data → phase → gas_data → gas →
     solution → electrochemistry_data → electrochemistry → kinetics → stoichiometry → organic →
     inorganic_data → inorganic → potential → nuclear_data → nuclear → spectroscopy → numerics →
     reaction_dynamics → network → fitting → sensitivity. (`gas_data` was added in M2a: the
     `PengRobinsonParams` struct and the `PR_*` presets are generated, so they precede `gas`.)
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
| **C — transcendental** | ulp ≤ K against the corrected oracle; K = 0 once the function is on correctly rounded primitives (ADR 0006) | anything through exp/ln/pow/log10/acos/atan2/cos: Arrhenius, Eyring, Nernst, pH, Clausius–Clapeyron, Antoine, PR roots, Morse, decay, Rydberg, isotope patterns |
| **D — iterative** | relative ≤ 1e-12 | Newton/bisection results, ODE trajectories, least squares, autodiff gradients, weak-acid pH where the iteration count differs |

kimiya specifics:
- **The class-C reference is the corrected oracle** ([ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md),
  decided 2026-11-09, in force since M3c): the oracle's operations, as `scripts/oracle/validate.py` transcribes them,
  with every transcendental primitive result correctly rounded (MPFR, `scripts/oracle/mpfr_ref.py`). Where glibc rounded
  every call of a row correctly that is the golden value; where it didn't, `tests/golden/<module>.cr.tsv` holds the
  corrected outputs, and the port is compared against those (accepted under ADR 6). Every transcendental primitive the
  port calls is correctly rounded and target-independent (`src/mathx.cyr`), so a function on them is compared **exactly**
  (K = 0): any difference is a port defect, and a port value that matches only glibc's misrounded one fails. Parity is
  functional — the oracle's results, `Err` codes, check order and branches — with the best accuracy available. The
  glibc-relative rule (K = 2·S·c + 2, measured per target from M0 step 8) set the bounds until M3c and is history
  (port-parity.md §2–§8); `validate.py --sensitivity` stays as a diagnostic of conditioning. A primitive the port has no
  correctly rounded version of yet (sin, atan2 until M4) keeps the stdlib's bound of 2 ulp.
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
| `exp` (18), `ln` (7) | **local** `_km_cr_exp`, `_km_cr_ln` (correctly rounded, two-phase; M3c-2) | class C, bound 0 against the corrected oracle (ADR 0006); the same bits on every target |
| `powf` (10) | **local** `_km_cr_pow` (correctly rounded, exact and midpoint cases included); `10.0.powf(x)` is `_km_cr_exp10(x)`, the same value; `x.powf(2.5)` is `_km_cr_pow(x, 2.5)` | class C, bound 0. ganita's pow (with its `pow(-1, huge)` defect, ADR 0004) is no longer called |
| `cos` (4), `acos` (1) | **local** `_km_cos` (correctly rounded for \|x\| ≤ 8; f64_cos beyond), `_km_acos` (double-double, ADR 0005) | class C, bound 0 (PR trig roots) |
| `sin` (1), `atan2` (2) | `f64_sin`, `ganita_f64_atan2` for now; correctly rounded `_km_` versions arrive with their first callers in M4 (ADR 0006 item 1) | class C, the stdlib's bound 2 until then (dihedral, spectroscopy) |
| `cbrt` (2) | `ganita_f64_cbrt` (correctly rounded) | exact: the oracle's (compiler-builtins) is correctly rounded on every golden argument (port-parity.md §9) |
| `log10` (4), `exp_m1` (3) | **local** `_km_cr_log10`, `_km_cr_expm1` (correctly rounded, M3c-2). M1's fdlibm `_km_log10` / `_km_expm1` and M3a's `_km_log10_ln` (glibc's composition) were removed in M3c-3 | class C, bound 0 (`henderson_hasselbalch`'s near-cancellation rows included) |
| `round` (4) | **local** `_km_round_half_away` (`f64_round` ties to even) | exact |
| `powi` (1) | repeated multiply in Rust's order | class B |
| `max` (13), `min` (7), `clamp`, `signum`, `total_cmp` | builtins / local helpers with Rust's NaN semantics (probe) | exact |

### 6.2 hisab: what the oracle used vs what hisab 3.4.1 provides

| Rust hisab 1.4.0 (used by) | hisab 3.4.1 Cyrius | Plan |
|---|---|---|
| `num::newton_raphson`, `num::bisection` (thermochem flame) | `num_newton`, `num_bisection` — semantics changed in 3.x (subnormal derivatives, NaN refusals, exact endpoint roots, a `lo + (hi − lo)·0.5` midpoint, `Err` on exhaustion), and their `fncall1` callbacks carry no context | **transcribed in M2b** (ADR 0003 amendment): `_km_flame_newton` / `_km_flame_bisection` in `src/thermochem.cyr`, 1.4.0's loops specialised to the flame balance; rows class B, 0 ulp |
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
  (`src/math/{log10,expm1}.rs`, fdlibm/musl, ≤ 1 ulp) so they sit within K = 2 of glibc (both removed in M3c-3 for the
  correctly rounded `_km_cr_log10` / `_km_cr_expm1`; see "Decided in M3c").
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

Decided in M2a (2026-11-09; reaction, thermo, phase, gas):
- **Struct values** — `PengRobinsonParams` is `struct KmPengRobinsonParams { tc: f64; pc: f64; omega: f64; }`,
  passed by value (as tanmatra's `TmFourMomentum`); `AntoineCoeffs` is `struct KmAntoineCoeffs` (five f64).
  Both are generated with their Rust field docs (`src/gas_data.cyr`, `src/phase_data.cyr`). A `PR_*` preset
  is a top-level struct value, `var KM_PR_WATER = KmPengRobinsonParams { tc: 647.1, … };` (probed: initialised
  data, passes to a by-value parameter). Rust's `{ field: x, ..p }` is a copy and a field store. Struct arrays
  can't be subscripted (6.7.18), so tests loop over per-case helper calls.
- **Generated module data** — K_b/K_f and the ANTOINE_DATA columns (name blob + offsets, A, B, C, t_min_c,
  t_max_c) in `src/phase_data.cyr`; STP and the presets in `src/gas_data.cyr`; `thermo::WATER_SPECIFIC_HEAT` in
  `src/constants.cyr` (thermo has no data file, and constants.cyr precedes it). Rust's *private* `PA_PER_MMHG`
  is generated too, through the API: `antoine_pressure_pa` of A = B = 0, C = 1 at T = 0 is 10⁰ × the factor
  (checked to be exactly 1.0 first), pinned by a class-A golden row and cross-checked against the source
  literal. Typed formula coefficients (PR's 0.45724, …) stay inline in the code, as in Rust.
- **ANTOINE_DATA handles** — `lookup_antoine(name)` returns row + 1 (0 for none); `antoine_data_name(h)` and
  `antoine_data_coeffs(h)` read a row (`""` and NaN fields for an invalid handle, which the Antoine functions
  reject). `KM_ANTOINE_COUNT` is the table's length.
- **ADR 0002 #4 error codes** — `antoine_pressure_mmhg` / `_pa`: a non-finite A, B, C, t_min_c or t_max_c is
  `Err(KM_ERR_INVALID_INPUT)`; T outside `t_min_c ..= t_max_c` (NaN included) or at/below the pole (C + T ≤ 0)
  is `Err(KM_ERR_INVALID_TEMPERATURE)`, the code kimiya uses for a bad temperature elsewhere; in range the value
  is Rust's (10^x may overflow to `Ok(+∞)` for a custom set, as in Rust). The harness maps the golden rows:
  oracle finite and T in range (with finite range ends) ↔ `Ok`, class C; oracle NaN or T out of range ↔ `Err`,
  counted as accepted under ADR 2 and summarised, not printed row by row.
- **Slices and the root vector** — a Rust `&[f64]` is a pointer and a count (`hess_law(xs, n)`,
  `enthalpy_from_formation(p, np, r, nr)`); Rust's f64 `sum` folds from −0.0, so an empty sum is −0.0; a negative
  count or a null pointer with values is NaN (no error channel). `peng_robinson_z_roots(params, t, p, out)` fills
  the caller's three-f64 buffer (`var roots[24];`) and returns `Ok(n)`.
- **`phase_rule`** takes Rust's u32 counts as words; one outside 0 ..= u32::MAX is `Err(KM_ERR_INVALID_INPUT)`,
  never wrapped (as M1 refuses an out-of-range atomic number).
- **Private Rust helpers** (`pr_coefficients`, `pr_dimensionless`, `pr_physical_roots`, `cubic_real_roots`,
  `validate_pr_params`) are `_km_*` functions in the same operation order; multi-value returns carry
  (code, A, B) and the double-double pairs.
- **Math** — `thermo` uses hisab's `EPSILON_F64` (1e-12, the value of hisab 1.4.0's) and includes
  `lib/hisab.cyr` explicitly: `cyrius distlib` keeps a non-stdlib include in place, so a consumer needs kimiya's
  `[deps.hisab]`. Modules include `lib/math.cyr` (f64_exp / f64_cos polyfills, f64_clamp) and `lib/ganita.cyr`
  (pow, cbrt) where they call them, so `dist/kimiya.deps` lists them. The Peng–Robinson and Redlich–Kwong
  functions use correctly rounded acos, cos, ln and x^2.5 (`src/mathx.cyr`, ADR 0005); see port-parity.md §5.
- **Faster guards** — `_km_is_finite(x)` is `x − x == 0`, `_km_is_positive` / `_km_is_non_negative` compare
  against +∞ directly (no nested calls): `gibbs_free_energy` 23.3 → 17.0 ns, `ideal_gas_pressure` 26.7 →
  16.7 ns; the `is_finite` corpus rows and every M1 row are unchanged.

Decided in M2b (2026-11-09; thermochem):
- **Struct methods by value** — `Shomate` is `struct KmShomate` (seven f64, generated with its Rust field docs in
  `src/thermochem_data.cyr`), and `Shomate::cp` / `Shomate::enthalpy_antiderivative` are methods of `impl KmShomate`
  taking `self` by value (probed on 6.7.18: `s.cp(t)` and `KmShomate_cp(s, t)` both work, as `KmMolecule`'s impl),
  so no free-function names were needed. One evaluation each (`_km_shomate_cp5` / `_km_shomate_h5`) serves the
  method and the table-row path (inlining it in the row path saved 2 of 22 ns; not worth a second copy).
- **Generated module data** — `THERMOCHEM_DATA` (formula and name blobs + offsets, ΔH_f°, ΔG_f°, S° columns,
  `KM_THERMOCHEM_COUNT`) and `SHOMATE_DATA` (formula blob, a … t_max columns, `KM_SHOMATE_COUNT`). The `//` comments
  inside both literals (sections, sources, the C₂H₆ / C₃H₈ / HF caveats) travel into each table's comment block.
  `cyr::string_column` now splits a blob line past 120 columns ("Phosphorus pentachloride") over 16-byte lines;
  the M1/M2a files are unchanged by it. thermochem has no `pub const` f64.
- **Handles** — `lookup_thermochem(formula)` → THERMOCHEM_DATA row + 1 (0 for none), read with
  `thermochem_data_formula` / `_name` (`""` for an invalid handle) and `thermochem_data_delta_hf_kj` / `_delta_gf_kj` /
  `_s_standard_j` (NaN); `lookup_shomate` / `lookup_shomate_at` → SHOMATE_DATA row + 1, read with
  `shomate_data_formula` and `shomate_data_coeffs` (a `KmShomate`, NaN fields for an invalid handle).
- **Species lists** — `&[(&str, f64)]` is a pointer to n 16-byte records (C-string pointer, f64 bits) and n, as
  `KmMolecule_new`'s pairs. A negative n or a null pointer with n > 0 is `Err(KM_ERR_INVALID_INPUT)`; a null formula
  pointer is an unknown formula. `reaction_*(products, np, reactants, nr)`, `adiabatic_flame_temperature(dh, products,
  n, t0)`.
- **ADR 0002 #1** — `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, computed as Rust's kJ difference then
  `× 1000.0` (one multiply); the golden rows hold the oracle's kJ result × 1000.0 (the generator multiplies; the row
  note says so), class B. `reaction_entropy` and the table fields are unchanged.
- **Rust's semantics kept** — `sum_property` from +0.0 with `+=` (the empty sum is +0.0, unlike `.sum()`'s −0.0), no
  coefficient guard (a NaN coefficient gives Ok(NaN)); the flame closures and the total moles fold from −0.0 (Rust's
  `.sum()`); `Segments` is a (first row, count) pair with `_km_segments_lookup` / `_cp` / `_enthalpy_change` in Rust's
  order (`_km_max` / `_km_min`, ±∞ bounds).
- **The flame solver** — hisab 1.4.0's `newton_raphson` and `bisection` transcribed as `_km_flame_newton` /
  `_km_flame_bisection` (ADR 0003 amendment, §6.2): class B, 0 ulp. `FLAME_NEWTON_TOL_J` (10⁻³), `FLAME_BISECTION_TOL`
  (10⁻⁶) and the 40 J/(mol·K) guess are typed formula constants inline, pinned by golden rows (`validate.py`'s
  pinning check) and by `tests/thermochem.tcyr` on the private helpers. Rust's `resolved` Vec is a stack table of 32
  products' segments (the bump allocator can't free, and the function is called in loops); a longer list resolves
  the rest again per evaluation (same values; golden rows with 40 products).

Decided in M3a (2026-11-09; solution, electrochemistry):
- **Constants** — the ten `solution::` constants (K_w, Debye–Hückel A and B, the seven Λ°) are generated into
  `src/constants.cyr` with their Rust doc comments, as M2a's `KM_WATER_SPECIFIC_HEAT` (solution has no data file, and
  constants.cyr precedes it); `KM_KW_WATER_25C` is written positionally (`0.00000000000001`). 35 of the 45 `constants`
  rows are compared; the 10 inorganic ones wait for M4.
- **Generated module data** — `STANDARD_POTENTIALS` in `src/electrochemistry_data.cyr` (`t_electrochemistry.rs`): the
  equation as the first **UTF-8** string blob (`cyr::string_column_utf8`, the bytes copied exactly: "Cu²⁺ + 2e⁻ → Cu" is
  22 bytes, and `strlen` counts bytes, as Rust's `&str` length does), the couple as an ASCII blob, `electrons` a `u8`
  column, `standard_potential` an f64 column, `KM_STANDARD_POTENTIALS_COUNT`; `HalfReaction`'s and the table's doc
  comments and the literal's two trailing "CRC Handbook, Electrochemical Series" notes travel into the comment block.
  `FARADAY` / `ELEMENTARY_CHARGE` stay `src/constants.cyr`'s.
- **Handles** — `lookup_half_reaction(couple)` → STANDARD_POTENTIALS row + 1, 0 for none (a null pointer is none);
  `HalfReaction`'s fields are named after the struct, as M1's `element_<field>(h)`: `half_reaction_equation(h)` and
  `half_reaction_couple(h)` (`""` for an invalid handle), `half_reaction_electrons(h)` (0: no row has 0) and
  `half_reaction_standard_potential(h)` (NaN). `cell_potential_from_couples` → `Err(KM_ERR_INVALID_REACTION)` for an
  unknown couple, the cathode's first.
- **Out-of-type integers** (as M2a's `phase_rule`, never wrapped) — Rust's `charge: i32` is a word: outside the i32 range
  the activity models return `Err(KM_ERR_INVALID_INPUT)`, checked first (Rust's type would reject it at compile time),
  and `ionic_strength` returns NaN (it has no error channel). Rust's `n_electrons: u8` is a word: outside 0 ..= 255 is
  `Err(KM_ERR_INVALID_INPUT)` at the position of Rust's `== 0` check, one test `n < 1 || n > 255` (so the temperature
  check still comes first in `nernst_potential`).
- **Ion records** — `ionic_strength(ions, n)`: Rust's `&[(f64, i32)]` is a pointer to n 16-byte records (the
  concentration's f64 bits, then the charge as an i64 word) and n, as M2b's species lists; Rust's `.sum()` folds from
  −0.0 (no ions: 0.5 × −0.0 = −0.0); a negative n or a null pointer with n > 0 is NaN.
- **The charge-balance Newton** — Rust's private `charge_balance_root` / `checked_root` / `validate_ionic_strength_and_a`
  are `_km_charge_balance_root` / `_km_checked_root` / `_km_validate_ionic_strength_and_a`, operation for operation
  (`_km_min` for `.min`, the step tested against the updated x). The private `CHARGE_BALANCE_REL_TOL` (10⁻¹²) and
  `CHARGE_BALANCE_MAX_ITER` (50) are typed formula constants inline: the tolerance is pinned from both sides by golden rows
  and by `tests/solution.tcyr` on the helper; the cap is unobservable (no input needs more than 7 updates), as the flame
  solver's limits (ADR 0003 amendment). The root is + − × ÷ and √, bit-identical; `weak_acid_ph` / `weak_base_ph` are
  class C through their one log10.
- **Math** — log10 is `_km_log10` (fdlibm) except in `henderson_hasselbalch`, which calls `_km_log10_ln`; the Nernst
  functions' ln is the correctly rounded `_km_ln`; 10^x is `ganita_f64_pow(10, x)` as the Antoine equation's. A
  correctly rounded log10 was tried first and failed the near-cancellation rows: glibc's log10 is not correctly rounded
  on 0.3–29 % of arguments (ADR 0005 amendment M3a, port-parity.md §7). `solution.cyr` includes `lib/ganita.cyr`,
  `electrochemistry.cyr` `lib/string.cyr`. (Superseded in M3c: `_km_cr_log10`, `_km_cr_exp10`, `_km_cr_ln`; solution no
  longer includes ganita.)

Decided in M3b (2026-11-09; kinetics):
- **Parity is functional (the user, 2026-11-09).** The port must give the oracle's results, errors and branches, with
  the best accuracy available, and never reproduce glibc's rounding errors: no primitive is chosen or built because it
  matches glibc's bits. A class-C function whose K would pass 64 keeps its rows *measured* (replayed unbounded, not
  gated) and its grid unchanged, until a later step (ADR 0006) makes the primitives correctly rounded and compares such
  rows against the oracle's formula evaluated with correctly rounded results. M3a's `henderson_hasselbalch` /
  `_km_log10_ln` and ADR 0005 are left for that step.
- **Math** — exp is `f64_exp`, exp_m1 `_km_expm1`, powf `ganita_f64_pow`; `kinetics.cyr` includes `lib/math.cyr` and
  `lib/ganita.cyr` (superseded in M3c: `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_pow`, no includes). All six class-C functions (`arrhenius_rate`, `first_order_concentration`, `nth_order_half_life`,
  `collision_theory_rate`, `eyring_rate`, `eyring_rate_from_activation`) are measured: a subnormal exp or pow result
  scaled into the normal range by a prefactor puts S at 10¹⁴ … 4.4 × 10¹⁸ (port-parity.md §8). The nine class-B functions
  are gated, 0 ulp.
- **Typed constants** — Rust's `std::f64::consts::LN_2` and private `FIRST_ORDER_TOLERANCE` (10⁻¹²) are typed inline, LN_2
  by its bits (`0x3FE62E42FEFA39EF`, as `src/mathx.cyr` writes ln 2's high part); `KM_GAS_CONSTANT`, `KM_BOLTZMANN`
  (Rust's `kinetics::BOLTZMANN`) and `KM_PLANCK` stay `src/constants.cyr`'s. Both typed constants are pinned by golden rows
  (`half_life_first_order(1)` is LN_2 exactly; the orders 1 ± 10⁻¹² and the ulps either side) and `validate.py`'s pinning
  check; `tests/kinetics.tcyr` pins the four boundary orders.
- **Tuple returns** — `lineweaver_burk`'s `Result<(f64, f64)>` follows `peng_robinson_z_roots`: `lineweaver_burk(v_max, km,
  out)` writes the slope (Km/Vmax) at `out` and the intercept (1/Vmax) at `out + 8` of the caller's two-f64 buffer
  (`var lb[16];`) and returns `Ok(2)`, the count written; Rust's errors are unchanged, and on `Err` the buffer is untouched.
  The golden row's output is the tuple's two floats (`f:slope f:intercept`, the format's tuple rule). As for the PR roots,
  `out` is the caller's contract: a null pointer is not checked.
- **Rust's semantics kept** — `zero_order_concentration` is `if remaining < 0.0 { 0.0 } else { remaining }` (NaN
  propagates, −0.0 stays −0.0); `collision_theory_rate`'s steric check is Rust's positive form `p > 0 && p <= 1` (NaN
  fails); `nth_order_half_life` checks `!is_finite(n) || |n − 1| <= 10⁻¹²`, then k, then [A]₀, and returns
  `Err(KM_ERR_COMPUTATION_ERROR)` unless t½ is positive and finite; 0·∞ inside a rate law (`arrhenius_rate(0, −1, 10⁻³⁰⁰)`,
  `eyring_rate_from_activation` with one exp overflowing and the other underflowing) is `Ok(NaN)`, as in Rust.

Decided in M3c (2026-11-09; ADR 0006, the user's decision of that day: "parity was in functionality … we still want to
offer the best accuracy as possible"):
- **The reference is the corrected oracle and the primitives are correctly rounded** (§5). M3c-1: `mpfr_ref.py` (MPFR 4.2
  through ctypes, development only), `validate.py --write-corrections` writes `tests/golden/<module>.cr.tsv`, and the
  harness compares a corrected row against the corrected outputs. M3c-2: `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`,
  `_km_cr_ln`, `_km_cr_log10` and `_km_cr_pow` in `src/mathx.cyr` (Ziv: a fast path with a proven error bound and a
  rounding test, a triple-double accurate path, the exact and midpoint cases explicit), with `_acc` / `_fast` entry
  points for the tests; their tables are generated into `src/mathx_data.cyr` by `scripts/oracle/mathx_tables.py` from
  MPFR (never typed), and `scripts/oracle/cr_verify.py` checks over 1,000,000 arguments per primitive against MPFR on
  every target. M3c-3: every module on them, every class-C function at bound 0.
- **One primitive per Rust call, in Rust's order.** `10.0_f64.powf(x)` is `_km_cr_exp10(x)` (a correctly rounded
  pow(10, x) is the same value; checked on 1.1M arguments); `tc.powf(2.5)` is `_km_cr_pow(tc, 2.5)` — ADR 0005's
  double-double `_km_pow_2_5` (85 ns against 130, but correctly rounded only from 2⁻³⁸⁷, with no rounding test) was
  removed. `_km_acos` and `_km_cos` (ADR 0005) stay. Modules that no longer call ganita don't include it; `gas` does, for
  cbrt.
- **glibc-matching fails.** A port value equal to the oracle's misrounded value on a corrected row means the port's
  primitive is not correctly rounded there: the harness fails it (M3c-1 had counted it as measured).
- **The grids keep and gain the hard rows** (ADR 0006 item 6): subnormal exp results scaled by a prefactor
  (`vant_hoff_k`, `clausius_clapeyron`, `fugacity_pr` at a few kelvin, `equilibrium_constant`'s subnormal end), on the
  shared targets `grids::EXP_SUBNORMAL_X` and seeded sweeps with fresh seeds.
- **`--sensitivity` is a diagnostic of conditioning**: every primitive result nudged, S and the retired K printed for
  comparison, no bound derived from it.

Language traps already probed on 6.7.18 by tanmatra (its plan §0, §5): no exponent literals;
`(0 - x)` on an f64 is wrong; a hex literal in a `const` is an integer; `f64_round` ties to
even; `+?` aborts the process; `cyrius test` must be warning-free.

---

## 7. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Initialize | ✅ done 2026-11-09 | constants 45 + math 18,647 (2.0 MB): 13,340 ok, 0 failed, 5,351 pending (M1 helpers, module-data constants), 1 accepted (ADR 0004) — identical on x86_64 ELF, x86_64 PE, aarch64 | — | scaffold only | Oracle frozen at tag `1.2.0`; scaffold + repairs; hisab 3.4.1; generator + `constants.cyr` (14 `KM_`) + harness; K = 2 ulp for the math primitives ([`port-parity.md`](port-parity.md)); ADRs 0001–0004 |
| M1 Foundation | ✅ done 2026-11-09 | element 975 + molecule 2,472 + constants 14 (of 45) + math 18,647 (now 0 pending): 22,107 ok, 0 failed, 31 pending (module-data constants), 1 accepted (ADR 0004) — identical on x86_64 ELF, x86_64 PE, aarch64; `molecular_weight` 0 ulp; log10 K = 2 (max 1), exp_m1 K = 2 (max 0) | 38 Rust tests ported (error 2, validate 3, element 22, molecule 11; molecule's `serde_roundtrip` doesn't port), plus 4 Cyrius-only groups: 179 assertions | 3 (`lookup_by_symbol` 1.5× Rust, `lookup_by_number` 5.4×, `molecular_weight_glucose` 7.4×; [benchmarks](../benchmarks-rust-v-cyrius.md)) | `error`, `mathx`, `validate`, generated `element_data` + `element`, `molecule`; `[lib]` + `dist/kimiya.cyr`; surface parity 55/55; settled points in §6.4 |
| M2 Thermochemistry & gases | ✅ **done 2026-11-09** (M2a: reaction, thermo, phase, gas; M2b: thermochem) | reaction 1,188 + thermo 1,107 + thermochem 9,487 + phase 3,276 + gas 19,717 + math 26,147 (7,500 new in M2a) + constants 25 of 45: 63,953 ok, 0 failed, 0 measured, 20 pending (M3–M4 constants), 441 accepted (ADR 0004 + 440 ADR 0002 #4 mapped Antoine rows) — identical on x86_64 ELF, x86_64 PE, aarch64; class B all 0 ulp (the flame temperature too: 1,682 Newton + 154 bisection rows); class C K = 2–10 by the rule (ADR 0005; vant_hoff_k 6) | 152 Rust tests ported (reaction 10, thermo 10, thermochem 62, phase 31, gas 39; none skipped; 4 thermochem tests assert J/mol for ADR 0002 #1, `antoine_pole_and_non_finite_give_nan` asserts Err for #4), plus 12 Cyrius-only groups: 729 assertions | 9 (M2a: `gibbs_free_energy` 9.7× Rust, `heat_transfer` 5.4×, `ideal_gas_pressure` 7.7×, `peng_robinson_z_n2_1mpa` 4.4×; M2b: `reaction_enthalpy_ch4` 6.5×, `vant_hoff_k` 6.2×, `enthalpy_change_cp` 6.9×, `shomate_cp` 6.2×, `adiabatic_flame_ch4_air` 10.7×) | generated `thermochem_data`, `phase_data`, `gas_data`; ADR 0005 (correctly rounded acos/cos/ln/x^2.5); ADR 0003 amendment (flame solver transcribed); ADR 0002 #1 (J/mol); surface parity 140/140; settled points in §6.4 |
| M3 Solutions & electrochemistry | ✅ **done 2026-11-09** (M3a: solution, electrochemistry; M3b: kinetics; its six kinetics class-C functions measured until M3c, now gated at 0) | solution 28,566 + electrochemistry 18,574 + kinetics 30,146 + constants 35 of 45: at M3b 119,044 ok, 0 failed, 22,205 measured (M3b's six class-C kinetics functions, K past 64 through subnormal exp / pow results: left for ADR 0006 by the user's decision, §6.4), 10 pending (the M4 inorganic constants), 441 accepted — identical on x86_64 ELF, x86_64 PE, aarch64; class B all 0 ulp; M3a class C K = 2–4 (Henderson–Hasselbalch and Nernst S ≈ 10¹⁸ before the remedy, 0 after, all 0 ulp); the weak acid / base root bit-exact (11,740 f ≤ 0, 6,356 relative-step, 74 non-finite exits; ≤ 7 updates); kinetics S 10¹⁴ … 4.4 × 10¹⁸ (2–4 without the subnormal rows) | 131 Rust tests ported (solution 49, electrochemistry 35, kinetics 47; 2 kinetics tests adapted for the Lineweaver–Burk buffer), plus 12 Cyrius-only groups: 1,667 assertions | 7 (M3a: `ph_from_h` 6.8× Rust, `weak_acid_ph` 5.2×, `lookup_half_reaction` 2.0×, `nernst_potential` 302× — the correctly rounded `_km_ln`; M3b: `arrhenius_rate` 7.3×, `eyring_rate` 6.3×, `michaelis_menten` 5.7×) | generated `electrochemistry_data` (the first UTF-8 blob) and 10 more constants; ADR 0005 amendment M3a (`_km_log10_ln`, `_km_ln` for Nernst; superseded by ADR 0006); the functional-parity decision (§6.4); surface parity 203/203; settled points in §6.4 |
| M3c Correct rounding (ADR 0006) | ✅ **done 2026-11-09** (M3c-1: the corrected oracle; M3c-2: the correctly rounded primitives; M3c-3: the switch) | 12 modules; gas 21,093 + phase 3,646 + reaction 1,558 + thermochem 9,857 (2,486 new hard rows: subnormal exp × prefactor) and 665 corrected rows in `tests/golden/<module>.cr.tsv`: **143,071 ok, 0 failed, 0 measured, 10 pending (the M4 inorganic constants), 1,105 accepted (665 under ADR 6, 440 ADR 0002 #4), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every class-C function at bound 0 and 0 ulp (sin, atan2 and 247 cos rows past \|x\| = 8 on the stdlib at 2 ulp until M4); ADR 0004's row a plain ok | `tests/mathx.tcyr` (the primitives' special, exact, tie, threshold, subnormal and hard cases against MPFR); the Cyrius-only groups that pinned glibc's bits pin the correctly rounded ones; no ported Rust test changed; 15 test files, 0 warnings; `cr_verify.py`: ≥ 1,000,000 arguments per primitive, 0 mismatches on every target | `vant_hoff_k` 9.51 ns (Rust) → 77 ns, `ph_from_h` → 87 ns, `weak_acid_ph` → 219 ns, `nernst_potential` 1.94 µs → 84 ns, `arrhenius_rate` → 66 ns, `eyring_rate` → 59 ns; `mathx/cr_*` benches ([benchmarks](../benchmarks-rust-v-cyrius.md)) | ADR 0006; `mpfr_ref.py`, `mathx_tables.py` (generated `mathx_data`), `cr_verify.py`; `_km_log10`, `_km_expm1`, `_km_log10_ln`, `_km_ln`, `_km_pow_2_5` removed; ADR 0005 amended (its M3a amendment superseded); port-parity.md §9–§10; settled points in §6.4 |
| M4 Structure, nuclear & spectroscopy | ☐ | | | | |
| M5 Numerics | ☐ | | | | ADR 0003 decides the hisab gaps |
| M7 Release 2.0.0 | ☐ | | | | |
