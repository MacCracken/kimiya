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
| Oracle's math library | Rust **std**: exp, exp2, expm1, ln, log10, pow, sin, cos, sincos, acos, atan2 → host **glibc 2.44**; cbrt, floor, trunc, round, powi → Rust's own built-ins. Release build (LLVM turns `2.0.powf` into `exp2`, and hisab 1.4.0's FFT `angle.cos()` / `angle.sin()` into `sincos`, M4e) | [`port-parity.md`](port-parity.md) §3, §15 |
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
| M4 | `organic.rs` | 238 | 7 | 21 | `src/organic_data.cyr` (generated, M4b), `src/organic.cyr` |
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
| 3 | `nuclear::semi_empirical_binding_energy` → NaN for z > a or a = 0 | `Result` (done in M4c) | finite rows class B; oracle NaN ↔ `Err` (23 rows, accepted under ADR 2 by name) |
| 4 | `phase::antoine_pressure_*` → NaN at the pole, silent extrapolation | `Result`, `Err` outside `t_min_c..t_max_c` | in-range rows class C; oracle NaN or out-of-range ↔ `Err` |
| 6 | `FunctionalGroup::hydrogen_bonding`, `Geometry::TShapedThreeLonePairs` | `is_hydrogen_bond_donor` (+ new `is_hydrogen_bond_acceptor`), `TShaped` | renamed fns class A; acceptor has reference tests only |
| 7 | `assign_oxidation_states` → `i32`, `Err` for average states; `balance_equation(&[Vec<f64>])` | exact rational state (`num`/`den`); integer matrix input | oracle `Ok(n)` ↔ `n/1`; oracle "non-integer" `Err` ↔ `Ok(p/q)`, checked by reference tests (Fe₃O₄ 8/3, C₃H₈ −8/3); balancing class A |
| 12 *(2026-11-09, M4b)* | `inorganic::cfse_octahedral` = `t2g·(−0.4) + eg·0.6`, off the exact (3·e_g − 2·t₂g)/5 on 13 of 22 cases (high-spin d⁵, d¹⁰ ≠ 0) | the exact value: integer numerator, one correctly rounded ÷ 5 | class A; equal rows plain; the 13 others mapped (q = round(5 × oracle), port = q/5), accepted under ADR 2 by name. Policy: other inexact-formula artifacts in M4 get their own row |
| 13 *(2026-11-09, M4c; completed in M4d)* | `nuclear::atoms_after_half_lives` = `initial * 2.0.powf(-n)`: ±∞, 0, NaN or up to 2⁵² ulp off where 2⁻ⁿ is not a normal double (n > 1022 or n ≤ −1024) though N₀/2ⁿ is representable | there N₀/2ⁿ correctly rounded (`_km_cr_exp2_scaled`: one rounding, 2⁻ᶠ past 53 bits; M4c's was faithful for a fractional n); ±0, ±∞, NaN N₀ give N₀ | class C; equal rows plain; 1,810 noted with the port's value, mapped, accepted under ADR 2, counted |
| 14 *(2026-11-09, M4c)* | `nuclear::is_secular_equilibrium` = `λ_p < 0.01·λ_d` (0.01 ≠ 1/100, the product rounds or underflows) | the exact λ_p < λ_d/100 (exponents, then `_km_two_prod`), Rust's NaN / ±∞ / ±0 semantics | class A; 208 noted rows mapped, by name |
| 15 *(2026-11-09, M4b's inorganic)* | `inorganic::delta_tet_from_oct` = `delta_oct * 4.0 / 9.0`: ±∞ where Δ·4 overflows (\|Δ\| ≥ 2¹⁰²²) | (Δ/9)·4 there: the correctly rounded 4Δ/9 everywhere | class B; 59 noted rows mapped, by name |
| 16 *(2026-11-09, M4d)* | `potential::lj_equilibrium_distance` = `sigma * 2.0.powf(1.0 / 6.0)`: 1/6 rounded, then the product rounded again (one ulp off on about a fifth of all σ) | the correctly rounded σ·2^(1/6) (`_km_cr_sixth_root2_mul`: a triple-double product rounded once; 2^(1/6)'s continued fraction proves no test is needed) | class C; 1,041 noted rows mapped, counted |
| 17 *(2026-11-09, M4d)* | `potential::lennard_jones` = `(4.0 * epsilon) * d`: ±∞ or NaN where 4·ε overflows (ε ≥ 2¹⁰²²) | (ε·d)·4 there: RN(4·ε·d), Rust's value wherever 4·ε is finite | class B; 180 noted rows mapped, counted |
| 18 *(2026-11-09, M4d)* | `potential::morse_width_parameter`'s `reduced_mass / (2.0 * d_e)`: 0 where 2·D_e overflows (D_e ≥ 2¹⁰²³) | (μ·½)/D_e there: RN(μ/(2·D_e)), Rust's value wherever 2·D_e is finite | class B; 181 noted rows mapped, counted |
| 19 *(2026-11-09, M4e)* | `spectroscopy::wavenumber_from_wavelength` / `wavelength_from_wavenumber` = `0.01 / x`: 1/100 rounded, the quotient rounded again (one ulp off on about 14 % of all x) | the correctly rounded 1/(100·x) (`_km_cr_hundredth_div`: a double-double residual and a rounding test, an exact integer division as the fallback) | class B; 827 noted rows each mapped, counted (as accepted although within B's 2 ulp) |
| 20 *(2026-11-09, M4e-2)* | `potential::distance_3d`, `bond_angle`, `dihedral_angle` square raw coordinates: an Err, 0, ±∞ or a degraded value past about 10^±154 (10^±77 for the dihedral), where the distance or angle is representable | Rust's sequence evaluated with an unbounded exponent (every + − × ÷ √ to 53 bits, no over- or underflow; comparisons and Err checks on those values), the final √ / atan2 rounded once; Rust's bits wherever its chain stays normal, exact scale invariance (`_km_ux_*`, a fast-path range test) | class B / C; 854 / 1,205 / 1,546 noted rows mapped (91 Err), counted |
| 21 *(2026-11-10, M4e-3)* | `spectroscopy::bohr_energy_level` ((−R·Z)·Z)/(n·n), `bohr_transition_energy` E(n₂) − E(n₁), `rydberg_wavelength` / `SpectralSeries::wavelength` 1/(((R·Z)·Z)·(1/n₁² − 1/n₂²)): three to eight roundings, −R one ulp off at 86 of 255 Z = n, cancellation for close large n (to 9.0·10⁷ ulp for λ and 1.8·10⁹ for E on the corpus) | the exact rational, the constant as its double, rounded once (a double-double fast path with a rounding test, the exact `_km_cr_ratio` slow path, ties to even) | class B; 3,823 / 1,351 / 1,376 / 67 noted rows mapped, counted |
| 22 *(2026-11-10, M4e-3)* | `spectroscopy::reduced_mass` = (m₁·m₂)/(m₁ + m₂): m/2 one ulp off on 8.6 % of masses; NaN, 0 or ∞ where m₁·m₂ or m₁ + m₂ over- or underflows | the correctly rounded m₁·m₂/(m₁ + m₂) for every pair (the smaller mass past a gap of 55, else 2^e₂·s₁s₂/(s₁ + s₂·2⁻ᵍ); subnormal results and ties exact) | class B; 2,837 noted rows mapped, counted |
| 23 *(2026-11-10, M4e-4)* | `spectroscopy::isotope_abundances`' table: AME2020 masses to five decimals (¹H 1.00783), no ³⁶S (sulfur sums to 0.9999); `isotope_pattern` convolves it | the cited table (`scripts/oracle/isotopes.tsv`): AME2020 masses at full published precision, the TICE 2009 abundances (Rust's), ³⁶S added; `isotope_pattern`'s algorithm unchanged | class A / B; 7 / 347 rows noted with the port's whole result `(v:n f:… …)`, mapped, counted |
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
     solution → electrochemistry_data → electrochemistry → kinetics → stoichiometry → organic_data →
     organic → inorganic_data → inorganic → potential → nuclear_data → nuclear → spectroscopy → numerics →
     reaction_dynamics → network → fitting → sensitivity. (`gas_data` was added in M2a: the
     `PengRobinsonParams` struct and the `PR_*` presets are generated, so they precede `gas`;
     `organic_data` in M4b: the three enums and their property columns are generated.)
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
     isotope abundances, PR parameters and Λ° inline in their modules. (M4b generated the organic
     bond energies too, into `organic_data.cyr` with the enums, and M2a / M3a generated the PR presets
     and Λ° after all.)
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
| **M4 Structure, nuclear & spectroscopy** | stoichiometry, organic, inorganic, potential, nuclear, spectroscopy | Exact integer balancing (checked i64/i128 → explicit overflow checks; class A) with integer input (#7); rational oxidation states (#7); VSEPR/bond tables; CFSE tables; Born–Landé; SEMF (#3) with `cbrt`; decay chains incl. equal-λ limit (`exp_m1`); Morse/LJ/Coulomb, atan2 dihedral; isotope pattern (pruned convolution, custom `total_cmp` sort); Rydberg/R_H; `find_peaks`; FFT via hisab `num_fft` (M4e: hisab 1.4.0's `num::fft` transcribed instead, ADR 0003 amendment) | balancing corpus (CH₄, KMnO₄/HCl, CuSO₄ both orders, seeded consistent matrices, every error class); oxidation corpus from the 1.2 tests; all 18 isotopes; SEMF Z ≤ A ≤ 3Z for A ≤ 300; isotope patterns up to C₁₆₈/H₈₈; angles at Å and metre scale; FFT n = 2..2²⁰ (M4e: the golden rows to 2¹³, every length to 16; `tests/spectroscopy.tcyr` to 2¹⁶) |
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
| **B — arithmetic** | ≤ 2 ulp, expected 0 | + − × ÷ √ and `cbrt` only: ideal gas, Van der Waals, molarity, ionic strength, Debye–Hückel limiting law, SEMF, mass defect, Shomate Cp and enthalpy, Born–Landé, LJ, Coulomb, Rydberg and isotope patterns (M4e: no transcendental in either) |
| **C — transcendental** | ulp ≤ K against the corrected oracle; K = 0 once the function is on correctly rounded primitives (ADR 0006) | anything through exp/ln/pow/log10/acos/atan2/cos: Arrhenius, Eyring, Nernst, pH, Clausius–Clapeyron, Antoine, PR roots, Morse, decay, transmittance, the FFT (its twiddles) |
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
  (port-parity.md §2–§8); `validate.py --sensitivity` stays as a diagnostic of conditioning. Since M4d every transcendental
  primitive has a correctly rounded port version (sin, cos, atan2 and acos from M4a), so no golden row keeps the stdlib's
  bound of 2 ulp.
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
| `cos` (4), `acos` (1) | **local** `_km_cr_cos`, `_km_cr_acos` (M4a: two-phase, correctly rounded on the whole double range; gas.cyr since M4d, which removed ADR 0005's double-double `_km_cos` / `_km_acos`) | class C, bound 0 (PR trig roots; the `math` rows past \|x\| = 8 too) |
| `sin` (1), `atan2` (2) | **local** `_km_cr_sin` (M4a; the `math` rows; spectroscopy's soorat-side `sin` is f32 and deferred), `_km_cr_atan2` (M4a; potential's `bond_angle` and `dihedral_angle` since M4d) | class C, bound 0 |
| hisab 1.4.0's FFT twiddle `(angle.cos(), angle.sin())` (glibc's `sincos` in the release oracle) | **local** `_km_cr_sincos` (M4a; `fft_spectrum`'s `_km_fft` since M4e) | class C, bound 0 |
| `0.01 / x` (2, spectroscopy's wavenumbers) | `_km_cr_hundredth_div(x)`, the correctly rounded 1/(100·x) (ADR 0002 #19) | class B, 1,654 rows mapped by note |
| `2.0.powf(1.0 / 6.0)` (1, potential) | `_km_cr_sixth_root2_mul(σ)`, the correctly rounded σ·2^(1/6) (ADR 0002 #16) | class C, bound 0, 1,041 rows mapped by note |
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
| `num::fft`, `Complex` (spectroscopy) | `num_fft`, `HComplex` — twiddles from the stdlib's `f64_cos` / `f64_sin` (not correctly rounded), reseeded every 64 steps: bits move from n ≥ 256 vs 3.3.x and vs 1.4.0 | **transcribed in M4e** (ADR 0003 amendment): `_km_fft` in `src/spectroscopy.cyr`, 1.4.0's loop operation for operation with `_km_cr_sincos` twiddles; class C at bound 0, every golden row 0 ulp |
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

Decided in M4b (2026-11-09; organic, inorganic):
- **Rust enums with property methods** — `BondType`, `Geometry`, `FunctionalGroup` and `Ligand` are generated enums
  (`enum KmBondType`, …, `KM_<ENUM>_COUNT`) with explicit append-only values in Rust variant order, under the Rust docs
  and the enum body's `//` section headings, and each `const fn` property is a generated column indexed by the value
  (`src/organic_data.cyr`, `src/inorganic_data.cyr`). `#[non_exhaustive]` forbids an exhaustive `match` outside the
  crate, so the generator proves its variant lists against the source declarations (`cyr::enum_decl`) and each match's
  arms against the variants (`cyr::match_arms`), and `validate.py` parses the enums again. A `matches!` set
  (`is_polar`, `hydrogen_bonding`) is a bool column; a formula over a column (`is_strong_field`, `field_strength() >= 8`)
  stays code, its constant typed inline. Member names are tanmatra's snake rule (`CCDouble` → `KM_BOND_TYPE_CC_DOUBLE`,
  `HCl` → `KM_BOND_TYPE_H_CL`); `surface_parity.py` and `probe_tables.py` now use that rule too (their regex gave
  `CCDOUBLE`).
- **Names** — the methods are free functions of the value, `<enum>_<method>(v)`: `bond_type_energy_kj(b)`,
  `geometry_ideal_bond_angle(g)`, `functional_group_pattern(g)`, `functional_group_is_polar(g)`,
  `ligand_field_strength(l)`, `ligand_is_strong_field(l)`; a word that is no value of the enum gets the property's
  "none": NaN, `""` or 0. **ADR 0002 #6**: `Geometry::TShapedThreeLonePairs` is `KM_GEOMETRY_T_SHAPED` (value 7, its
  Rust position); `FunctionalGroup::hydrogen_bonding` is `functional_group_is_hydrogen_bond_donor`; the new
  `functional_group_is_hydrogen_bond_acceptor` reads a column the generator classifies by IUPAC 2011's criterion
  (Arunan et al., *Pure Appl. Chem.* 83, 1637: an electron-rich region such as a lone pair or a π-bonded pair), the weak
  acceptors included by the user's decision: every group but the alkane accepts — O / N lone pairs (McMurry, *Organic
  Chemistry*, OpenStax 2023, §2.12), X–H···π to alkenes, alkynes and arenes (Steiner, *Angew. Chem. Int. Ed.* 41, 48,
  2002), divalent sulfur (Allen et al., *Acta Cryst.* B53, 696, 1997), C–X halogens (Brammer et al., *Cryst. Growth
  Des.* 1, 277, 2001; organic F the weakest, Dunitz & Taylor, *Chem. Eur. J.* 3, 89, 1997) — with the sources in the
  generated comment, reference tests only, and `validate.py` recomputing the rule (a lone-pair atom or a π bond in the
  pattern).
  The harness compares the renamed items' golden rows under the Rust names, class A, and summarises them.
- **Enum outputs and inputs** — `predict_geometry`'s `Option<Geometry>` is the value-form `Option` of a `KmGeometry`
  value; golden `e:Enum.Variant` tokens map to values through the harness's name tables.
- **Records** — `enthalpy_from_bonds(broken, nb, formed, nf)`: Rust's `&[(BondType, u32)]` is a pointer to n 16-byte
  records (the enum value as a word, then the count as a word), as M2b's species lists. Rust's two `.sum()`s fold from
  −0.0 (two empty lists give +0.0) and `energy * (n as f64)` per term. No error channel: a bad list, an invalid value
  or a count outside u32 is NaN.
- **Out-of-type integers** (as M2a / M3a, never wrapped) — Rust's `u8` (d electrons, unpaired electrons, the VSEPR
  counts), `u32` (ion charges) and `i8` (an ionic radius's charge) are words; outside the type's range a Result function
  returns `Err(KM_ERR_INVALID_INPUT)` (`cfse_octahedral`, `unpaired_electrons_oct`, `born_lande_lattice_energy`, where
  the u32 test sits at Rust's `== 0` check), a bare f64 is NaN (`spin_only_moment`), and an `Option` / lookup is None / 0
  (`predict_geometry`, `lookup_ionic_radius`).
- **Bool inputs** — a Rust `bool` is a word: 0 is false, 1 is true. In a Result function any other word is
  `Err(KM_ERR_INVALID_INPUT)`, checked where Rust branches on it (`cfse_octahedral` / `unpaired_electrons_oct` pick the
  table first, so before the count). A non-Result function with a bool argument returns its documented "none" for any
  other word (NaN, 0, None), as for an out-of-type integer (none in M4b).
- **Private tables through the API** — Rust's private `OCT_UNPAIRED_*` are `unpaired_electrons_oct(n, spin)`; the
  private `OCT_OCCUPANCY_*` (t₂g, e_g) pairs are the unique pair with t₂g + e_g = n whose `t2g · (−0.4) + eg · 0.6` gives
  `cfse_octahedral`'s bits, as M2a's `PA_PER_MMHG`; each is cross-checked against its literal and generated as `_KM_`
  u8 columns. Δ_tet's 4 / 9, the moment's and Born–Haber's 2.0, Born–Landé's 1.0, the strong-field 8 and the exact
  CFSE's divisor 5 are typed inline and pinned by golden rows (`validate.py` perturbs each; the 8 is pinned from above
  only: no ligand has rank 7).
- **ADR 0002 #12, the exact CFSE (the user's decision)** — `cfse_octahedral` returns (3·e_g − 2·t₂g)/5, the integer
  numerator converted to f64 and divided once (correctly rounded), not Rust's `t2g · (−0.4) + eg · 0.6`, which misses it
  on 13 of the 22 cases (high-spin d⁵ is −2⁻⁵² and d¹⁰ −2⁻⁵¹ instead of 0). The golden rows keep the oracle's value, the 13
  differing ones name the exact value in their note, and the harness maps them (q = round(5 × oracle), port = q/5) and
  reports each by name as accepted under ADR 2; the rows are class A. **Policy:** other M4 formulas whose Rust rounding
  misses a simple exact rational or a chemically exact value get the same treatment, each as its own ADR 0002 row.
  Reviewed in M4b and left as Rust's (questions for the user): `delta_tet_from_oct` is already correctly rounded (Δ·4 is
  exact, so (Δ·4)/9 is CR(4Δ/9)) except where Δ·4 overflows, |Δ| > DBL_MAX/4 (5 golden rows, ±∞) — decided in M4c: fixed
  exactly, ADR 0002 #15, below; Born–Landé's
  (1 − 1/n) misses CR((n − 1)/n) on 332 of 3,263 Ok rows (n = 7, Pauling's Ne exponent, among them), but the energy is a
  product of rounded physical constants, so no exact value exists to reach; Born–Haber is a sum of arbitrary inputs; the
  spin-only moment √(n(n + 2)) is already correctly rounded; organic's energies sum exact integers.
- **Constants and handles** — the ten `MADELUNG_*` / `BORN_EXP_*` are generated into `src/inorganic_data.cyr`
  (`t_constants::consts_where`, each group under its source section comment): the last pending `constants` rows, so
  `Home::Module` is gone and nothing is pending. `lookup_ionic_radius(symbol, charge)` → IONIC_RADII row + 1, 0 for None;
  `IonicRadius`'s fields are `ionic_radius_symbol(h)` (`""`), `ionic_radius_charge(h)` and
  `ionic_radius_coordination(h)` (0: no row has either) and `ionic_radius_radius_pm(h)` (NaN), named after the struct as
  `half_reaction_*` are.
- **Classes** — the enum properties, tables, lookups, `predict_geometry`, `unpaired_electrons_oct` and (exact since ADR
  0002 #12) `cfse_octahedral` are class A; the arithmetic (`enthalpy_from_bonds`, Δ_tet, the moment, Born–Landé,
  Born–Haber) is class B, all 0 ulp.

Decided in M4c (2026-11-09; nuclear):
- **Generated module data** — `ISOTOPES` in `src/nuclear_data.cyr` (`t_nuclear.rs`): the label as an ASCII blob, `z` a `u8`
  column, `a` the first `u16` data column, `half_life_s` f64, `decay_mode` a **UTF-8** blob ("β⁻/β⁺/EC" is 14 bytes), and
  `KM_ISOTOPES_COUNT` (18 rows: §4's "all 18 isotopes"). `Isotope`'s and the table's docs and the literal's three NUBASE2020
  notes (C-14, I-131, Pu-239) travel into the comment block, keyed by row and field. `AMU_MEV`, `PROTON_MASS_U`,
  `NEUTRON_MASS_U` and `HYDROGEN_ATOM_MASS_U` stay `src/constants.cyr`'s.
- **Handles** — `lookup_isotope(symbol)` → ISOTOPES row + 1, 0 for None (a null pointer or `""` too); `Isotope`'s fields
  are `isotope_symbol(h)` and `isotope_decay_mode(h)` (`""` for an invalid handle), `isotope_z(h)` / `isotope_a(h)` (0: no
  row has either) and `isotope_half_life_s(h)` (NaN), named after the struct as `ionic_radius_*` are.
- **Typed formula constants** — Rust's `LN_2` (by its bits, as `src/kinetics.cyr`), the private `EQUAL_RATE_REL_TOL`
  (10⁻¹²) and `SEMF_VOLUME` / `_SURFACE` / `_COULOMB` / `_ASYMMETRY` / `_PAIRING` (15.67, 17.23, 0.714, 93.15, 11.2) and the
  asymmetry term's 4 are typed inline, not generated through the API: no single call isolates a SEMF coefficient
  (`PA_PER_MMHG` had one). (`is_secular_equilibrium`'s 0.01 was one until ADR 0002 #14 made the comparison exact.) Each is pinned by golden rows — `validate.py` perturbs it one ulp
  (the tolerance also × ½ and × 2) and requires some row to change (the SEMF coefficients 35–19,404 rows, LN_2 about
  2,000) — and the chain's rows *on* the threshold (δ = fl(10⁻¹²·λ_max) exactly, and δ where 10⁻¹² ± 1 ulp moves the
  branch: 738 rows) pin the tolerance to the ulp, which M3b's first-order tolerance could not be.
- **ADR 0002 #3** — `semi_empirical_binding_energy(z, a)` returns a `Result`: `Err(KM_ERR_INVALID_INPUT)` for a word
  outside u32, then z > a (Rust's `checked_sub`), then a = 0; finite rows class B (cbrt is ganita's, correctly rounded
  as Rust's compiler-builtins' is; √ is IEEE). The golden rows keep the oracle's NaN on the 23 Err rows, noted
  `[ADR 0002 #3: the port returns Err(InvalidInput)]`; the harness maps oracle NaN ↔ `Err(InvalidInput)` (any other
  outcome fails) and prints each row by name as accepted under ADR 2, then a summary after the module. Its accepted
  line no longer prints a got / want bit pair when the port's output is not a float (`parity_core.cyr`).
- **Out-of-type integers** (as M2a / M3a / M4b, never wrapped) — Rust's u32 `z`, `n`, `a` are words: outside
  0 ..= u32::MAX `binding_energy_per_nucleon` and `semi_empirical_binding_energy` return `Err(KM_ERR_INVALID_INPUT)` at the
  position of Rust's first check (the checked add / sub, the same code), `mass_defect` NaN. Z + N past u32::MAX is
  Rust's `checked_add` overflow, `Err` (2³¹ + 2³¹; 2³¹ − 1 + 2³¹ is Ok).
- **The decay chain** — `decay_chain_two_step(n0, l1, l2, t, out)` writes (N_A, N_B, N_C) to the caller's 24-byte buffer
  and returns `Ok(3)`, as `lineweaver_burk`; on `Err` the buffer is untouched. Rust's equal-rate branch and its threshold
  are kept exactly (`delta <= 1e-12 · max(λ₁, λ₂)`, `_km_max` / `_km_min` for Rust's `max` / `min`), and Rust's two
  computations of −λ₁t are one value.
- **Slices** — `q_value(r, nr, p, np)`: each `&[f64]` a pointer and a count, summed by `src/thermo.cyr`'s `_km_sum`
  (Rust's `.sum()` from −0.0; NaN for a negative count or a null pointer with values, M2a's rule) — the first module to
  call another module's private helper, which the `[lib]` order allows (thermo precedes nuclear); `tests/nuclear.tcyr`
  includes `src/thermo.cyr`.
- **Math** — exp, exp_m1 and ln are `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_ln`; `2.0_f64.powf(−n)` is `_km_cr_pow(2, −n)`
  (the release oracle calls glibc's exp2, LLVM's rewrite; `validate.py`'s replica calls Python's `math.exp2` and the
  corrected oracle MPFR's exp2, a new `_prim` kind `"exp2"`: the same real 2^y); cbrt `ganita_f64_cbrt`. `nuclear.cyr`
  includes `lib/ganita.cyr` and `lib/string.cyr`; `dist/kimiya.deps` is unchanged.
- **Performance** — `binding_energy_per_nucleon` calls the private `_km_mass_defect` (no second u32 check) and multiplies by
  `KM_AMU_MEV` inline instead of calling `mass_defect` and `binding_energy_mev`: 41 → 29 ns back to back on the same host
  (27 ns in the final run; Rust 2.02 ns).
- **ADR 0002's inexact-formula policy, reviewed — and decided** (the user, 2026-11-09: "fix exactly"). The review found
  `atoms_after_half_lives` (20 golden rows where 2⁻ⁿ over- or underflows though N₀/2ⁿ is representable: 1000 after 1075
  half-lives gave 0, not 2.47 × 10⁻³²¹) and `is_secular_equilibrium` (0.01·λ_d against the exact λ_d/100, 66 rows); with
  M4b's Δ_tet overflow they became:
  - **ADR 0002 #13** — `atoms_after_half_lives` keeps Rust's N₀·CR(2⁻ⁿ) where 2⁻ⁿ is a normal double (−1024 < n ≤ 1022),
    which is the correctly rounded N₀/2ⁿ for an integer n. For a finite n beyond (`_km_half_lives_edge`): ±0, ±∞, NaN N₀
    return N₀ (scaleB); |n| > 2200 is ±0 or ±∞ with no arithmetic; N₀ = ±2^j with j − n exact (a `_km_two_sum` check) is
    ±`_km_cr_pow(2, j − n)`, correctly rounded; otherwise N₀ = s·2^j, n = k + f, s·CR(2⁻ᶠ) = p + e (`_km_two_prod`, exact),
    then `_km_cr_finish(j − k, p, e, 0)` (`src/mathx.cyr`'s single rounding: exact scaling for a normal result, one
    rounding at the 2⁻¹⁰⁷⁴ grid for a subnormal one, ties to even, ±∞ on overflow). That is the correctly rounded N₀/2ⁿ for
    an integer n (f = 0) and **faithful** for a fractional one — CR(2⁻ᶠ)'s relative 2⁻⁵³ is the one error before the final
    rounding, so it is within one ulp, correctly rounded on 743 of the 1,121 golden rows of that path and 82 % of a 4,000-pair
    sweep; true correct rounding needs 2⁻ᶠ past 53 bits (a correctly rounded scaled exp2, a mathx follow-up). An infinite n
    keeps Rust's limits. 1,812 golden rows differ from the oracle: 564 integer-n, 173 power-of-two, 1,069 faithful, 6 ±0 / ±∞
    N₀.
  - **ADR 0002 #14** — `is_secular_equilibrium` is the exact λ_p < λ_d/100 (`_km_cmp_hundred`: zeros and signs first, then
    the binades, then 100·s_p by `_km_two_prod` against s_d·2⁶ or 2⁷, exactly); NaN false and ±∞ as in Rust. 208 golden
    rows differ, among them (0.03, 3), (−0.02, −2), (0, 5 × 10⁻³²⁴) and (RN(DBL_MAX/100), DBL_MAX).
  - **ADR 0002 #15** (M4b's inorganic) — `delta_tet_from_oct` takes (Δ/9)·4 where Δ·4 overflows: CR(4Δ/9) everywhere.
    59 golden rows differ (the 5 of M4b and a new sweep over [2¹⁰²², DBL_MAX]).
  - **The mapping** — the corpus keeps the oracle's values; a differing row's note names the port's (`… = x (f:<bits>)]` or
    `(i:<0|1>)]`), computed by the generator in integer arithmetic (#13, #14) or as (Δ/9)·4 (#15), proven by `validate.py`
    with `Fraction` and MPFR. The harness reads the value from the note (`_par_note_*`, set per row by `parity_core.cyr`),
    requires the port's bit for bit, and reports the row by name as accepted under ADR 2 (without a bit pair: the note
    has both values).
  Kept as Rust's under the policy: `decay_constant`'s LN_2/t½ (1 ulp from CR(ln 2/t½) on 476 of 2,175 rows; ln 2/t½ is
  irrational, a product of rounded constants). Not misses: `time_for_fraction(1, λ)` is −0.0 (Rust's −ln 1), the exact
  value 0; the SEMF's N = Z asymmetry and odd-A pairing are already exactly 0; the mass defect, binding energies and Q-values
  are sums and products of rounded constants and inputs.
- **Classes** — the table, `lookup_isotope` and `is_secular_equilibrium` are class A; `decay_constant`, `activity`,
  `mass_defect`, `binding_energy_mev`, `binding_energy_per_nucleon`, `semi_empirical_binding_energy` and `q_value` class B
  (0 ulp); `atoms_remaining`, `atoms_after_half_lives`, `time_for_fraction` and `decay_chain_two_step` class C at bound 0
  (`atoms_after_half_lives`' 1,812 #13 rows compared bit for bit with the note's value).

Decided in M4d (2026-11-09; potential, the trig switch, ADR 0002 #13's completion and #16–#18):
- **Points** — Rust's `[f64; 3]` arguments are pointers to three contiguous f64 (24 bytes, `var p: f64[3];`):
  `distance_3d(p1, p2)`, `bond_angle(a, b, c)`, `dihedral_angle(a, b, c, d)`. A null pointer (0) gives NaN from
  `distance_3d` and `Err(KM_ERR_INVALID_INPUT)` from the two Result functions, checked first; Rust's own errors (a
  non-finite coordinate, coinciding atoms, collinear atoms) are `Err(KM_ERR_INVALID_INPUT)` as in Rust. The golden rows
  flatten each point into three `f:` tokens.
- **Private helpers** — Rust's `sub3`, `cross`, `dot3` are `_km_sub3` and `_km_cross` (three results, a multi-value
  return) and `_km_dot3`, on scalars (a point is loaded once into three f64 locals), in Rust's operation order: `dot3` is
  (a0·b0 + a1·b1) + a2·b2.
- **Typed constants** — the private `COINCIDENT_REL_TOL` (10⁻¹², squared at run time as Rust does: `tol_sq = tol·tol`)
  and `COLLINEAR_SIN_SQ_TOL` (10⁻²⁰) are typed inline, `.max` is `_km_max`, and π (the dihedral's −π → π) is
  `std::f64::consts::PI` by its bits. The tolerances are numerical resolution thresholds and stay Rust's (as
  `EQUAL_RATE_REL_TOL`, unlike #14's chemical 1/100); golden rows on each threshold pin them to the ulp (`validate.py`:
  ± 1 ulp changes 544 / 454 bond-angle rows and 108 / 140 dihedral rows, × ½ and × 2 more).
- **Math** — exp is `_km_cr_exp` (Morse), atan2 `_km_cr_atan2` (`bond_angle`'s atan2(|ba × bc|, ba · bc),
  `dihedral_angle`'s atan2(|b₂|·b₁·n₂, n₁·n₂) with Rust's (−π, π] mapping), sqrt `f64_sqrt`. `validate.py`'s replica
  gains the `_prim` kind `"atan2"` (glibc's, corrected by MPFR's). The release oracle's `2.0.powf(1.0 / 6.0)` is modelled
  as exp2(RN(1/6)) (glibc's and MPFR's agree: RN(2^(1/6)) = 0x3FF1F59AC3C7D6C0).
- **The trig switch** — gas.cyr calls M4a's `_km_cr_acos` and `_km_cr_cos`; every one of its 21,093 golden rows keeps its
  bits (the Peng–Robinson z, roots and fugacities all 0 ulp, as before). ADR 0005's double-double `_km_cos`, `_km_acos`
  and their helpers `_km_dd_cos_sin`, `_km_dd_asin`, `_km_dd_add`, `_km_dd_mul`, `_km_dd_div_d`, `_km_dd_sqrt` are
  removed (nothing called them); `src/mathx.cyr` no longer includes `lib/math.cyr` (gas.cyr does, for `f64_clamp`;
  `dist/kimiya.deps` lists the same six leaves, `math` now after `str`, where gas first needs it). The `math` corpus's `f64::sin`, `f64::cos` (all 2,910 rows, the 247 past |x| = 8
  included), `f64::acos` and `f64::atan2` go to `_km_cr_sin`, `_km_cr_cos`, `_km_cr_acos` and `_km_cr_atan2` at bound 0;
  `f64_sin`, `f64_cos` and `ganita_f64_atan2` are shadow lines. No golden transcendental row is on the stdlib any more.
- **ADR 0002 #13 completed** — `atoms_after_half_lives` past the normal range is `_km_cr_exp2_scaled(N₀, −n)`, the
  correctly rounded N₀·2⁻ⁿ (src/mathx.cyr, M4d: exp's fast core with a rounding test at the final scale, exp's accurate
  core in triple-double as the fallback; verified against MPFR by `cr_verify.py`, `validate.py` and `tests/mathx.tcyr`).
  The generator's notes are the correctly rounded value (`scripts/oracle/src/exact.rs`: 2⁻ᶠ to 320 bits in integers).
- **ADR 0002's inexact-formula policy, reviewed for potential** (the user's policy: fix exactly (a)–(c), leave chains of
  rounded operations): **#16** `lj_equilibrium_distance` (a: 1/6 rounded, the product rounded twice), **#17**
  `lennard_jones`' 4·ε overflow and **#18** `morse_width_parameter`'s 2·D_e overflow (b, #15's shape: an exact constant
  multiple of an input overflowing at the top of the range) — fixed. Kept as Rust's and reported to the user (no single
  exact value to restore): (σ/r)⁶ / (σ/r)¹² over- and underflow, k·q₁·q₂ and r·r over- and underflow, the squared
  coordinates' over- and underflow in the geometry (an Err or a degraded value where the angle or distance is
  representable; Rust documents the angles as scale-invariant; fixed in M4e-2 as ADR 0002 #20, the user's decision),
  `morse`'s (r − r_e) overflow and subnormal D_e·x, and the two tolerances.
- **Classes** — `lennard_jones`, `morse_width_parameter`, `coulomb`, `coulomb_force`, `distance_3d` class B (0 ulp);
  `lj_equilibrium_distance` (Rust's powf), `morse`, `bond_angle`, `dihedral_angle` class C at bound 0.
- **Performance** — Rust's one potential bench, `potential/dihedral_angle` (the gauche geometry, whose atan2 is the
  x = 0 special case): 130 ns against Rust's 20.35 (6.4×); Cyrius-only `potential/bond_angle` 215 ns (a full atan2),
  `lennard_jones` 35, `lj_equilibrium_distance` 109, `morse` 66, `mathx/cr_exp2_scaled` 151.

Decided in M4e (2026-11-09; spectroscopy, ADR 0002 #19, ADR 0003's FFT amendment):
- **Enum and methods** — `SpectralSeries` is the generated `enum KmSpectralSeries` (`KM_SPECTRAL_SERIES_LYMAN` = 0 …
  `PFUND` = 4, append-only, proven against the source declaration) with `KM_SPECTRAL_SERIES_N_LOWER` (u8) through the API;
  the methods are `spectral_series_n_lower(s)` (0 for a word that is no value) and `spectral_series_wavelength(s, line)`
  (`Err(KM_ERR_INVALID_INPUT)` for such a word, checked first).
- **Out-of-type integers** (as M2a–M4d, never wrapped) — the u8 Z and u32 n / line at Rust's `== 0` checks
  (`z < 1 || z > 255` in one test), so a Z of 256 is an `Err`, not Rust's Z = 0; NaN for the bare-f64
  `harmonic_oscillator_energy` / `rigid_rotor_energy`; the empty slice for `isotope_abundances`; in `isotope_pattern`, per
  atom in list order, a Z outside 0 ..= 255 is `Err(KM_ERR_INVALID_ELEMENT)` (as an unknown Z: an atomic number is never
  aliased, M1's rule), then a count outside u32 `Err(KM_ERR_INVALID_INPUT)`.
- **Vec results and slices** — a Rust `Vec` result goes to the caller's buffer `out` of `cap` elements (f64, 16-byte
  (mass, intensity) records, or i64 indices): `Ok(count)`; `Err(KM_ERR_INVALID_INPUT)` when it doesn't fit, checked after
  Rust's own errors, with the buffer untouched (`find_peaks` counts first only when cap < (n − 1)/2, the most strict
  maxima n values can have); a negative cap or a null `out` holds nothing. A slice argument is a pointer and a count;
  a negative count or a null pointer with values is `Err(KM_ERR_INVALID_INPUT)` (every slice function here has an error
  channel). `isotope_abundances(z)` returns its slice as `(pointer, count)` into the baked table ((0, 0) for none).
- **Allocation** — `isotope_pattern`'s working pattern and scratch start in two stack tables of 192 records and move to
  the bump allocator (which never frees) only past that (glucose's scratch peaks at 30 records, the Rust test's
  C₁₀₀₀H₂₀₀₀O₅₀₀S₅₀Br₂₀ at 189; C₅₀₀₀ or a thousand Br, Cl or S atoms need 200–546 and allocate); `fft_spectrum` transforms in `out` when cap ≥ 2n (`signal` may then be `out`: it is read from the end) and
  otherwise allocates 16·n bytes; an allocation failure is `Err(KM_ERR_COMPUTATION_ERROR)` (Rust's FFT error maps there,
  unreachable).
- **The stable sort** — Rust's `sort_by(|a, b| a.0.total_cmp(&b.0))` is stable: the port sorts by an insertion sort on the
  total-order key (Rust's `bits ^ ((bits >> 63) >>> 1)`), which keeps equal masses in scratch order and is near linear on
  the convolution's interleaved runs. No natural molecule's scratch has equal masses with unequal intensities (the
  commutative ties are identical records), so `tests/spectroscopy.tcyr` pins the order on the helper with −NaN, ±0, +NaN
  and equal masses.
- **Typed constants** — `ISOTOPE_PRUNE_REL` (10⁻¹²), `ISOTOPE_REPORT_REL` (0.001) and `ISOTOPE_BIN_WIDTH_U` (0.5) are typed
  inline in the helpers that apply them (`_km_isotope_normalize`, `_km_isotope_report`, `_km_isotope_bin`); × ½ / × 2 change
  110–173 golden rows each (the bin width's × ½ none), ± 1 ulp none (no peak lands on a threshold), so
  `tests/spectroscopy.tcyr` pins each to the ulp on its helper. π (ħ's 2π, B's 8π², the FFT's −2π) is
  `std::f64::consts::PI` by its bits; the 2.0, 8.0, 0.5 and the FFT's −2.0 are pinned by golden rows (`validate.py`).
  Rust's private R_H (`RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO)`, const-evaluated) is generated **through the API**:
  `SpectralSeries::Lyman.wavelength(u32::MAX − 1)` is RN(1/R_H) (1/1 − 1/u32::MAX² rounds to 1), which has exactly one
  double reciprocal; it is cross-checked against the expression on the source's literal ratio, also generated.
- **Classes** — `isotope_abundances`, `SpectralSeries::n_lower` and `find_peaks` (indices) are class A; everything + − × ÷
  and √ is class B, `isotope_pattern` and the Rydberg formula included (no transcendental: §5's examples were wrong);
  `transmittance` (exp10), `absorbance_from_transmittance` (log10) and `fft_spectrum` (the twiddles) class C at bound 0.
- **The FFT** — hisab 1.4.0's `num::fft` transcribed as `_km_fft` (ADR 0003 amendment M4e), not hisab 3.4.1's `num_fft`
  (§6.2); the magnitude is Complex::abs's √(re·re + im·im), not hypot (its spurious overflow for |re| past 10¹⁵⁴ is
  Rust's, kept).
- **ADR 0002's inexact-formula policy, reviewed for spectroscopy** — **#19** (case (a): the wavenumbers' 0.01 = 1/100)
  fixed. Kept as Rust's and reported to the user (ADR 0002's M4e addition has the list): the Bohr levels' (−R·Z)·Z/(n·n)
  where Z/n gives −R/k², `reduced_mass(m, m)` against m/2, the Rydberg factor's cancellation for close large n, the
  over- and underflow of products and quotients of inputs, `find_peaks`' interpolated threshold, and the isotope
  thresholds 0.001 / 10⁻¹².
- **Harness** — `_PAR_MAX` is 8448 tokens (the 2¹³-sample FFT rows); a row the dispatcher maps to an ADR's exact value is
  accepted even where its class's tolerance passes it (`_par_force_map`, #19's one-ulp class-B rows).
- **Performance** — `spectroscopy/absorbance` 11 ns (Rust 0.98), `bohr_energy_level` 16 ns (Rust 1.21),
  `isotope_pattern_glucose` 16.8 µs (Rust 2.80, 6.0×); Cyrius-only `fft_spectrum_1024` 134 µs, `find_peaks_1000` 22.5 µs,
  `wavenumber_from_wavelength` 58 ns.

Decided in M4e-2 (2026-11-10; potential's geometry with an unbounded exponent, ADR 0002 #20, the user's decision of
2026-11-09, "Yes, prescale"):
- **The target** — `distance_3d`, `bond_angle` and `dihedral_angle` are Rust's own operation sequence (every subtraction,
  product, sum, comparison, `max` and the Err checks, in Rust's order) evaluated in IEEE binary64 round-to-nearest-even
  with an unbounded exponent range; the final √ or atan2 is correctly rounded to binary64 itself (its subnormal grid, or
  ±∞); comparisons and Err branches use the unbounded values; non-finite inputs keep Rust's behaviour. Rows whose chain
  stays in the normal range keep Rust's bits (`validate.py` proves it on every such golden row); the angles are exactly
  scale-invariant and the distance scales by exactly 2^k where normal.
- **Fast path and test** — Rust's chain as before, plus a cheap proof that it stayed in the normal range: s (`distance_3d`)
  and |ba|², |bc|², max(|a|², |b|²), max(|c|², |b|²), |ba × bc|² (`bond_angle`, in place of Rust's two `is_positive`
  calls) in [2⁻⁹⁰⁰, DBL_MAX] with ba · bc finite; every component of b₁, b₂, b₃ zero or in [2⁻²²⁹, 2²⁵⁴) for the dihedral
  (its squares inline). The proofs are in ADR 0002's M4e-2 addition and the source comments. Where the test fails and the
  coordinates are finite, `_km_distance_3d_ux`, `_km_bond_angle_ux`, `_km_dihedral_angle_ux` (Rust's sub3 / cross / dot3
  on unbounded values); a non-finite coordinate gives Rust's result (NaN or +∞ for the distance, Err for the angles,
  whose first check then always fails).
- **The helpers** — `src/mathx.cyr`'s `_km_ux_*`: (significand, exponent) pairs, each operation one IEEE operation on the
  significands (exact 53-bit rounding), the final operations rounded once with the exact error breaking subnormal
  midpoints (`_km_ux_round` over `_km_cr_finish`), atan2 on `_km_cr_atan2`'s (significand, exponent difference) core.
  General on purpose: the machinery of the library-wide audit after M4.
- **The library-wide audit** (the user's decision, the same day) — after M4, every ported module is checked for
  intermediate over- or underflow of products and quotients of inputs, and each edge is fixed the same way, in-range rows
  keeping Rust's bits; not done in M4e-2 (`state.md`'s follow-up).
- **Golden comparison** — notes `(f:…)` or `(E:InvalidInput)` on exactly the rows where the port's result differs from
  the oracle's and the unbounded chain from Rust's (3,605: 854 / 1,205 / 1,546); the generator's integers
  (`exact::Ux`, `exact::ux_atan2_f64`), `validate.py`'s MPFR evaluator (`mpfr_ref.UX`, `potential_*`), the harness's
  `_par_got_noted20`; 4,823 new edge rows; `cr_verify.py`'s eight new families on all three targets.
- **Performance** — `potential/dihedral_angle` 131 → 133 ns, `bond_angle` 213 → 213, `distance_3d` 14 → 16; the slow path
  314 ns, 1.13 µs and 1.31 µs on the 2⁶⁰⁰-scaled geometries.

Decided in M4e-3 (2026-11-10; the exact Rydberg rationals and the exact reduced mass, ADR 0002 #21 and #22, the user's
decisions of 2026-11-10, "Compute exactly" and "Round once", and the whole expression rounded once, to which the user did
not object, "Committed, continue"):
- **The target** — `bohr_energy_level` RN(−R·Z²/n²), `bohr_transition_energy` RN(R·Z²·(n₂² − n₁²)/(n₁²·n₂²)),
  `rydberg_wavelength` and `spectral_series_wavelength` RN(n₁²·n₂²/(R·Z²·(n₂² − n₁²))), with the constant as its double
  (`KM_RYDBERG_EV`, `KM_RYDBERG`, `_KM_RYDBERG_HYDROGEN`); `reduced_mass` RN(m₁·m₂/(m₁ + m₂)) for every pair of positive
  finite doubles, subnormal results and ties included, no intermediate over- or underflowing. Rust's validation, Err codes
  and check order are kept; an integer outside its Rust type keeps the port's Err.
- **Fast paths** — the exact numerator and denominator as double-doubles (exact for the levels, the reduced mass and
  n₂ < 2¹³; within 2⁻¹⁰²·⁹⁹ each otherwise, `_km_dd_mul`), the quotient and its rounding test by `_km_dd_div_cr`
  (y0 + d within 13.12u² of the double-double quotient, ε = 2⁻¹⁰⁰·y0; fails with probability about 2⁻⁴⁶ and on every
  midpoint); the reduced mass is the smaller mass past an exponent gap of 55 (proven), else 2^e₂·s₁s₂/(s₁ + s₂·2⁻ᵍ).
- **Exact slow path** — `src/mathx.cyr`'s multi-word unsigned integers (`_km_bn_*`, 8 limbs of 32 bits) and the general
  `_km_cr_ratio(p, q, e)` = RN(p/q·2^e) (restoring long division to a 55-bit quotient, the remainder as a sticky bit,
  `_km_cr_from_int`'s single rounding), written for reuse by M4f and the overflow audit.
- **Ties** — none for the levels (n/Z a power of two gives −R/4^k exactly; otherwise at least 56 significant bits) and
  none for the wavelengths (R∞'s and R_H's odd significands carry a prime above 2³²: never dyadic); the transition
  energy's exhaustive search finds 13, (255; 2^k·(58,055, 524,544)), and 5,657 exact values; the reduced mass has whole
  families (subnormal m = m, (6c, 2c)·2⁻¹⁰⁷⁴, normal (2²⁷ ∓ c) and others). Every tie found is a golden row.
- **Golden comparison** — notes `(f:…)` where the port's value differs (#21: 3,823 / 1,351 / 1,376 / 67; #22: 2,837),
  computed by `exact::cr_ratio`; `validate.py` proves each with `Fraction` and every unnoted row equal to the correctly
  rounded value, and re-derives the structure; `cr_verify.py`'s `ratio`, `bohr`, `trans`, `ryd`, `series`, `rmass`
  against `Fraction` on all three targets; 13,706 new rows.
- **The audit's list** — `reduced_mass`' m₁·m₂ and m₁ + m₂ leave the library-wide overflow audit's candidates (#22 covers
  them).
- **Performance** — `spectroscopy/bohr_energy_level` 9 → 49 ns (Rust 1.21), `bohr_transition_energy` 27 → 73,
  `rydberg_wavelength` 24 → 73, `series_wavelength_balmer_alpha` 28 → 78, `reduced_mass_hcl` 14 → 73; the exact path
  2.42 µs and 2.20 µs on a tie.

Decided in M4e-4 (2026-11-10; the cited isotope table, ADR 0002 #23, the user's decision, "Upgrade the data"):
- **The data** — `isotope_abundances` returns the AME2020 masses at full published precision (Wang et al., Chin. Phys. C
  45, 030003 (2021): `mass_1.mas20.txt`, read by its Fortran format, its digits kept) and the TICE 2009 representative
  abundances (Berglund and Wieser, Pure Appl. Chem. 83, 397–410 (2011), Table 1, column 9: Rust's values), ³⁶S added; each
  the double nearest its decimal. CIAAW's current table gives these seven elements intervals only, so TICE 2009's single
  values remain the representative compositions. `isotope_pattern` is Rust's algorithm, unchanged.
- **The source** — `scripts/oracle/isotopes.tsv`, one record per isotope with the AME2020 line verbatim and the TICE value
  as printed, the URLs, sha256s and references in its header; the generator checks it against the oracle (the same
  isotopes plus ³⁶S, Rust's abundances as decimal strings, Rust's masses the AME2020 ones to five decimals, every element
  summing to exactly 1) before writing `src/spectroscopy_data.cyr`; `scripts/oracle/isotope_sources.py` re-verifies it
  against fresh downloads.
- **Golden comparison** — the rows keep the oracle's values; where the port's result differs the note names the whole
  result, `(v:n f:… …)]` (the note's trailing token list generalised from one value), from the generator's
  `pattern_replica` (checked against the oracle on every row with Rust's table); the harness compares those tokens bit
  for bit; `validate.py` proves each note with its own replica and every unnoted row unchanged (7 + 347 rows).
- **Performance** — `spectroscopy/isotope_pattern_glucose` 15.8 µs either way; the new sulfur benches
  `isotope_pattern_cystine` 16.9 → 18.5 µs and `isotope_pattern_s8` 8.4 → 18.8 µs (sulfur's fourth isotope; S₈'s scratch at
  most 68 records, within the stack tables).

Decided in M4f (2026-11-10; stoichiometry, ADR 0002 #7 realised; each point the option closest to Rust's observable
behaviour, listed for the user in M4f's report):
- **The rational states (#7)** — `assign_oxidation_states(atoms, n, charge, out, cap)`: n 16-byte (Z, count) records in,
  n 24-byte (Z, numerator, denominator) records out, the denominator positive and the fraction reduced, `Ok(n)`. Rust's
  private `exact_oxidation_state` is `_km_exact_oxidation_state(num, den)` → (code, p, q): Rust's checks in order, the
  exact average where Rust returned Err "non-integer average oxidation state", and Rust's i32 range check applied to the
  rational's **value** (outside i32 → Err(KM_ERR_COMPUTATION_ERROR): Rust's observable code for those inputs, whose
  non-integer Err came first; they need counts near 2³¹). The peroxide branch's −2 < os < 0 test is made on the rational
  (p < 0 < p + 2q, formed without overflow). Each entry's state is a function of its Z once the solved element is known,
  so the port needs no work table: the records are written in one pass at the end (`out` untouched on Err).
- **The integer matrix (#7)** — `balance_equation(matrix, m, n, out, cap)`: flat row-major i64 words (entry (i, j) at
  i·n + j), the n u32 coefficients to `out` as words, `Ok(n)`. Rust's |entry| ≤ 2⁵³ (MAX_COMPOSITION_ENTRY) is kept
  (Err(KM_ERR_INVALID_INPUT) past it): every integer matrix Rust accepts is accepted, every one it refuses refused. A null
  matrix, m < 1, n < 2 or m·n past 2⁵⁶ is Err(KM_ERR_INVALID_INPUT). Rust's INTEGER_TOLERANCE and ragged-row check have
  no counterpart; the harness converts an oracle row's f64 entries as Rust's test does (finite, within 10⁻⁹ of
  round(v), that integer an i64) and replays the integers class A (41 rows had a near-integer entry); a row without a
  counterpart (55 non-integer, 11 non-finite, 7 beyond i64, 5 ragged) keeps the oracle's Err(InvalidInput), is noted
  `[ADR 0002 #7: no integer-matrix counterpart (…)]` and counted as accepted under ADR 2 (the harness requires the note
  and the Err).
- **Checked arithmetic** — Rust's i64 `checked_add` / `_sub` / `_mul` / `_neg` are `_km_i64_*_ck` (a (1, value) or
  (0, 0) pair; Cyrius's + − × wrap, so each tests its wrapped result; the trapping i64::MIN / −1 is refused first, and
  Rust's plain `lcm / piv` division, which would panic there, is an Err instead: unreachable, the pivots are positive).
  `reduce_row`'s gcd is over unsigned words to 2⁶³ (`_km_gcd_u64`, `_km_umod`); a row of gcd 2⁶³ is Rust's `try_from`
  overflow (a corpus row: [2³¹, 1, 0], [2³¹, −(2³² − 1), 0]). The i128 verification is a signed three-word sum
  (`_km_i128_mul_u32`, `_km_i128_add_ck`: h·2⁶⁴ + m·2³² + l with h an i64), since `src/mathx.cyr`'s `_km_bn_*` are
  unsigned; Rust's `checked_add` edge at −2¹²⁷ is kept, its `checked_mul` can't fail (|term| < 2⁸⁵). Multi-value
  returns take at most three values, so the accumulator is a 24-byte buffer updated in place.
- **Proven fast paths** (no check can fail there, so Rust's results are unchanged): an elimination step with all four
  operands below 2³¹ in magnitude is formed directly; each verification row is summed in i64 first (entries below 2³¹,
  no partial sum overflowing) and goes to the i128 sum otherwise; `reduce_row` stops folding at a gcd of 1; the
  oxidation sums' checks are inline (|os| ≤ 2 and a count below 2³² make Rust's products exact).
- **Out-of-type integers** (as M2a–M4e, never wrapped) — a Z outside 1 ..= 118 is Err(KM_ERR_INVALID_ELEMENT) (256 is no
  hydrogen); a count outside 1 ..= u32::MAX Err(KM_ERR_INVALID_INPUT) at Rust's `count == 0` check (one test); a charge
  outside i32 Err(KM_ERR_INVALID_INPUT), checked first. `stoichiometry` precedes `inorganic` and `nuclear` in the `[lib]`
  order, so it tests these ranges inline rather than calling their `_km_is_u8` / `_km_is_u32`.
- **Slices and buffers** — `limiting_reagent(moles, nm, coeffs, nc)` and `theoretical_yield_moles(…, p)` take two f64
  slices (a negative count or a null pointer with values Err(KM_ERR_INVALID_INPUT), before Rust's length checks, the same
  code); `percent_composition(molecule, out, cap)` takes a `KmMolecule` by value and writes 16-byte (Z, percent) records;
  `empirical_formula(elements, n, out, cap)` reads 16-byte (Z, f64 bits) records and writes (Z, count) records. Every Vec
  result: the cap checked after Rust's errors, `out` untouched (M4e's rule).
- **Work tables** — `balance_equation`'s matrix copy and null vector and `empirical_formula`'s mole ratios are
  function-local tables (512, 64 and 64 words) and move to the bump allocator only past them (the corpus's 30- and
  70-species chains take that path); an allocation failure is Err(KM_ERR_COMPUTATION_ERROR).
- **Typed constants** — Rust's private EMPIRICAL_TOLERANCE (0.1) is typed inline in `_km_empirical_close`, MAX_COMPOSITION_ENTRY
  (2⁵³) in `_KM_MAX_COMPOSITION_ENTRY`; `validate.py` requires a row to change when the tolerance is halved or doubled and
  when the bound moves by one (both pinned); ± 1 ulp of 0.1 changes no golden row (|r·k − round(r·k)| is a multiple of
  the ratio's ulp, at least 2⁻⁵²), so `tests/stoichiometry.tcyr` pins 0.1 to the ulp on the helper.
- **Classes** — `balance_equation`, `assign_oxidation_states`, `limiting_reagent` (an index) and `empirical_formula`
  (counts) class A; `theoretical_yield_moles`, `percent_yield`, `percent_composition` class B, 0 ulp (+ − × ÷ only).
- **Harness** — `_par_map_toks(adr, tag)` takes the note's tag (#23's and #7's whole-result notes); `_par_in_matrix`
  converts the matrix. 9,090 oxidation rows and 78 balancing rows mapped under #7.
- **The five float functions, measured, not fixed** (the inexact-formula policy; the user decides): `percent_yield`
  (one ulp from RN(100a/t) on about a quarter of inputs, 0 for some subnormal results), `theoretical_yield_moles` (one
  ulp on 15–26 %, ±∞ or 0 where m/c over- or underflows though m·p/c is finite), `limiting_reagent` (another index than
  the exact comparison on near-ties and where both ratios over- or underflow), `percent_composition` (up to 2 ulp from
  RN(mass·count·100/mw) on a quarter of entries), `empirical_formula` (exact ratios and an exact 1/10 flip Ok ↔ Err on 14
  corpus rows); and `balance_equation`'s intermediate i64 overflow where a valid u32 solution exists (for the overflow
  audit). ADR 0002's M4f addition has the numbers.
- **Performance** — `stoichiometry/balance_ch4_combustion` 1.20 µs → 0.77 µs with the fast paths (Rust 233.90 ns),
  `oxidation_states_kmno4` 156 → 122 ns (Rust 22.38); Cyrius-only `balance_kmno4_hcl` 3.84 → 2.16 µs,
  `oxidation_states_fe3o4` 132 → 110 ns, `limiting_reagent` 31, `theoretical_yield_moles` 46, `percent_yield` 14,
  `percent_composition_glucose` 70, `empirical_formula_glucose` 217 ns.

Decided in M4f-2 (2026-11-10; the five float functions of M4f made exact, ADR 0002 #24, #25 and #26, the user's
decisions "Round once", "Compare exactly" and "Evaluate exactly"; the M4f points on the typed tolerance and the five
measured functions are superseded by these):
- **The targets** — `percent_yield` RN(100·a/t), `theoretical_yield_moles` RN(m·p/c) at #25's index, `percent_composition`
  RN(mass·count·100/mw) with `molecular_weight`'s own mw (#24); `limiting_reagent` with the exact mᵢ·c_{i₀} < m_{i₀}·cᵢ,
  the first index on an exact tie (#25); `empirical_formula` with exact moles, minimum, ratios, tolerance 1/10 and counts
  (#26). Rust's checks, Err codes and order unchanged; a zero operand keeps Rust's ±0.
- **One helper for #24** — `_km_cr_mul_div(x, y, z)` = RN(x·y/z) for positive finite doubles, the third use of the pattern
  (100·count < 2³⁹ is an exact double, so the composition is x = mass, y = 100·count, z = mw). Its fast path is
  `_km_dd_div_cr`'s with exact operands, plus an exact midpoint resolution where the rounding test fails
  (`_km_mul_div_tie`: four 64-bit integers), since midpoints are common with integer coefficients; `_km_cr_ratio` only
  for subnormal operands or an exponent past ±1021. Kept in `src/stoichiometry.cyr` (its only user), as spectroscopy
  keeps the reduced mass's paths.
- **#25's fast path is Rust's own** — the rounded ratio RN(m/c) decides unequal ones (RN is monotone, over- and underflow
  included); only equal rounded ratios call `_km_ratio_lt`. `limiting_reagent` costs what Rust's did (33 ns).
- **#26** — a ratio of 2³² − ½ or more is Err(KM_ERR_COMPUTATION_ERROR) whichever way the search goes (no multiplier, or
  a count of at least 2³²), so the port returns it once the ratios are known and keeps only ratios below 2³²; a fast
  double-double tolerance test that proves its pass or failure or defers to the exact integer test; EMPIRICAL_TOLERANCE
  is the exact 1/10 (Rust's 0.1 stays pinned on its replica in `validate.py`; the port's 1/10 by the rows exactly on the
  edge); `_km_empirical_close` is removed; `_km_empirical_formula(…, force)` forces the exact tests for `cr_verify.py`.
  The work table holds a double-double per element (128 words local).
- **Rust's empirical Err paths** — "non-positive min" is gone (its inputs Ok or the same ComputationError), "ratio not
  finite" is gone without changing any result (every such ratio is past 2³²); ADR 0002's M4f-2 addition lists the rest.
- **Harness** — a #24 yield row maps through `_par_got_noted`, a #25 row through `_par_map_exact` with an `i:` value, the
  composition and empirical rows through `_par_map_toks`; each counted per function (1,126, 240, 598, 920, 270).
- **Performance** — `percent_yield` 14 → 81 ns, `theoretical_yield_moles` 46 → 115, `percent_composition_glucose`
  72 → 276, `empirical_formula_glucose` 217 → 586, `limiting_reagent` 32 → 33; the exact paths
  `theoretical_yield_moles_tie` 133 ns, `percent_yield_subnormal` 2.15 µs, `limiting_reagent_tie` 94 ns,
  `empirical_formula_edge` 7.46 µs.

Decided in the overflow audit's B1 (2026-11-10; ADR 0002 #27 and #28, the user's decisions "Library-wide audit" of
2026-11-09 and "Round once" of 2026-11-10, quoted in ADR 0002's B1 addition):
- **The treatment rule** — a function whose whole formula is one product or quotient of inputs and constants returns that
  exact value rounded once (a physical constant as its double, an integer literal of Rust's formula exact); a formula with
  a sum inside gets ADR 0002 #20's semantics (Rust's sequence with an unbounded exponent, the final operation rounded
  once, Rust's bits wherever its chain stays normal). B1: `heat_transfer`, `specific_heat_from_calorimetry`, the two
  colligative functions, `osmotic_pressure` (#27), `zero_order_half_life`, `steady_state_intermediate`,
  `pre_equilibrium_rate` rounded once and `second_order_concentration`, `michaelis_menten` under #20's semantics (#28).
  Rust's checks, Err codes and order unchanged everywhere; where Rust gives NaN but the value is out of range the port
  gives the correctly rounded ±∞ or 0 (none in B1).
- **One general helper** — M4f-2's `_km_cr_mul_div` moved from `stoichiometry.cyr` to `src/mathx.cyr` as the 2/1 shape of
  `_km_cr_mono(num, nn, den, nd)` (at most 6 and 4 finite factors, the denominators nonzero; the exact path's capacity
  253 bits of odd significands a side; NaN outside), with the wrappers `_km_cr_mul3`, `_km_cr_div2`, `_km_cr_mul3_div`,
  `_km_cr_mul5`; signs, zeros (the factors' xor), ties to even, subnormal results and overflow handled inside. The fast
  path takes subnormal operands (normalised exactly) and resolves midpoints for up to three numerator and two
  denominator factors by an exact expansion sign (`_km_mono_tie`, `_km_exp_sign`: Shewchuk's Grow-Expansion), so only
  possibly subnormal results and the rare failed test past those shapes reach `_km_cr_ratio`. The general fast path
  is slow in this compiler (its factor tables and loops: `heat_transfer` 148 ns a call through it, 65 through its own
  path), so each shape the modules call has a
  straight-line fast path of its own for normal operands — `_km_mul_div_fast` (M4f-2's, with its integer midpoint
  resolution `_km_mul_div_tie`, now for both signs), `_km_mul3_fast`, `_km_div2_fast`, `_km_mul3_div_fast`,
  `_km_mul5_fast`, on the shared `_km_dd_round_ends` and `_km_mono_fin` — and falls back to `_km_cr_mono`.
- **`zero_order_half_life` keeps Rust's division as its fast path** — wherever 2k is finite it is exact, so the one IEEE
  division is the correctly rounded [A]₀/(2k); the monomial serves k ≥ 2¹⁰²³.
- **The two #20-style fast paths are functions** (`_km_second_order_fast`, `_km_michaelis_menten_fast`: Rust's chain and
  its range test, thresholds 2⁻¹⁰²¹ … DBL_MAX for products and quotients, finiteness for sums of non-negatives), so
  `cr_verify.py` and the tests check them against the slow paths (`_km_second_order_ux`, `_km_michaelis_menten_ux`).
- **Unvalidated functions** (`heat_transfer`, the colligative pair): a non-finite input keeps Rust's expression (IEEE's ±∞
  and NaN), as #20 keeps non-finite coordinates.
- **The corpus** — the oracle's values kept; each module's B1 edge rows appended by `scripts/oracle/src/audit.rs` with
  their own seed (thermo 537, phase 528, kinetics 1,064); 682 #27 and 636 #28 rows noted and mapped.
- **Performance** — `heat_transfer` 5.4 → 65 ns, `specific_heat_from_calorimetry` 18 → 86, the colligative pair 6 → 66–67,
  `osmotic_pressure` 20 → 139, `steady_state_intermediate` 21 → 88, `pre_equilibrium_rate` 27 → 125,
  `zero_order_half_life` 16 → 16, `second_order_concentration` 23 → 33, `michaelis_menten` 17 → 24; stoichiometry as in
  M4f-2.

Decided in the overflow audit's B2 (2026-11-10; ADR 0002 #29 and #30, B1's decisions unchanged; the survey's F3 and F5,
its proposed B2 and B3 in one batch):
- **The treatment, per function** — the monomials rounded once: `dilution`, `solution_conductivity` (1000 exact),
  `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` (Rust's two equivalent weights are quotients
  inside one quotient: RN(m_A·M_B·n_A/(n_B·M_A))) (#29), `absorbance`, `concentration_from_absorbance`,
  `rigid_rotor_energy` (J and J + 1 exact), `rotational_constant` (π as its double each time Rust uses it, 8 exact),
  `diatomic_moment_of_inertia` (#30); `nernst_potential`, with a sum inside, under #20's semantics, the correctly
  rounded ln Q (`_km_cr_ln`) its primitive (#29); `vibrational_frequency` the correctly rounded √ of the exact k/μ (#30,
  the planned treatment; its fallback, #20's semantics with `_km_ux_sqrt_f64`, was not needed). Rust's checks, Err codes
  and order unchanged; no NaN of these functions stands for an out-of-range value.
- **Four more straight-line shapes** on B1's `_km_cr_mono`: `_km_mul2_div2_fast`, `_km_mul3_div2_fast`, `_km_mul4_fast`,
  `_km_div3_fast` (wrappers `_km_cr_mul2_div2`, `_km_cr_mul3_div2`, `_km_cr_mul4`, `_km_cr_div3`); the 2/2 and 3/2
  midpoints fall back to `_km_cr_mono`'s fast path, which resolves them. `rotational_constant` calls `_km_cr_div3` with
  h/8 (exact).
- **The correctly rounded √ of a quotient** — `_km_cr_sqrt_div(x, y)` in `src/mathx.cyr`: a double-double fast path
  (`_km_sqrt_div_fast`, proven within 5.7u², `_km_dd_div_cr`'s rounding test, `_km_mono_fin`'s exact scale) and an exact
  path (`_km_sqrt_div_acc`: R = ⌊√(A·2¹¹⁰/B)⌋ by a guess corrected with multi-word comparisons, then
  `_km_cr_from_int(2R + sticky)`); √(x/y) is never a rounding boundary of binary64 (proven), so the exact path never
  breaks a tie.
- **The Nernst fast path is a function** (`_km_nernst_fast`: R·T in [2⁻¹⁰²², DBL_MAX], (R·T)/(n·F) and the product with
  ln Q at least 2⁻¹⁰²² or ln Q = 0), checked against the slow path `_km_nernst_ux` by `cr_verify.py` and the tests;
  `nernst_potential_25c` always takes it.
- **Unvalidated inputs** (`solution_conductivity`, `absorbance`, `diatomic_moment_of_inertia`, the rotor's B,
  `concentration_from_absorbance`'s A): a non-finite one keeps Rust's expression, as in B1.
- **The corpus** — the oracle's values kept; each module's B2 edge rows appended by `scripts/oracle/src/audit.rs` with
  their own seed (solution 458, electrochemistry 823, spectroscopy 1,178); 2,756 #29 and 2,106 #30 rows noted and
  mapped. The generator's Nernst notes use glibc's ln, which `validate.py` proves correctly rounded on every row.
- **Performance** — the monomials 6–33 → 68–148 ns (`absorbance` 5 → 68, `faraday_mass_ratio` 33 → 148),
  `vibrational_frequency` 18 → 87, `nernst_potential` 83 → 92 (its range test); the edges `nernst_potential_edge` 202
  ns, `vibrational_frequency_subnormal` 883 ns (the exact root), `faraday_mass_ratio_tie` 402 ns.

Decided in the overflow audit's B3 (2026-11-11; ADR 0002 #31 and #32, B1's decisions unchanged; the survey's F2 and F6,
its proposed B4 and B5 in one batch):
- **The treatment, per function** — the monomials rounded once: `ideal_gas_pressure`, `ideal_gas_volume`, `gas_density`
  (#31), `coulomb`, `coulomb_force` (#32); `morse_width_parameter` the correctly rounded ω·√(μ/(2·D_e)) (#32: a monomial
  under a square root, as #30's √(k/μ)); with a sum inside, under #20's semantics: `van_der_waals_pressure`,
  `peng_robinson_pressure`, `redlich_kwong_pressure` (#31), `lennard_jones`, `born_lande_lattice_energy` (its (1 − 1/n)
  Rust's, as M4b decided: in-range rows keep Rust's bits), `binding_energy_per_nucleon` (the mass defect Z·m_H + N·m_n − M
  is a sum of inputs) (#32). Rust's checks, Err codes and order unchanged; the co-volume tests V ≤ n·b and V ≤ b are
  decided on the unbounded co-volume (the semantics' comparisons), so the port can refuse where Rust accepts and the
  reverse; Rust's NaN where the value is out of range (`lennard_jones` past σ/r = 2^170.7, the equations of state) is the
  correctly rounded ±∞ or 0, and the unbounded evaluation gives the correctly signed zero where Rust's is +0.
- **#17 and #18 folded into #32** — #17's (ε·d)·4 is the unbounded chain's value wherever Rust's d is normal (all 180 rows
  kept), and stays the fast path's branch for ε ≥ 2¹⁰²²; #18's (μ·½)/D_e gives way to the whole correctly rounded Morse
  width (128 of its 181 rows kept, 53 moved). Their ADR rows record the fold; the harness counts #32's rows in the two
  regions.
- **The correctly rounded w·√((x/y)·2^e)** — `_km_cr_mul_sqrt_ratio(w, x, y, e)` in `src/mathx.cyr` (B2's
  `_km_cr_sqrt_div` generalised by a factor and a power of two): the fast path B2's root times w (within 9.8u²) with
  `_km_dd_round_ends`' test, the exact path R = ⌊√(W²·A·2⁶/B)⌋ by multi-word comparisons and `_km_cr_from_int(2R + sticky)`;
  unlike √(x/y) the value can be a midpoint (ω·√(9/(2·½)) = 3ω), which the exact path rounds to even.
- **The #20-style fast paths are functions** (`_km_vdw_fast`, `_km_pr_pressure_fast`, `_km_rk_pressure_fast`,
  `_km_lj_fast`, `_km_born_lande_fast`; the binding energy's one test inline), checked against the slow paths by
  `cr_verify.py` and the tests. Their range test is `src/mathx.cyr`'s new `_km_ux_normal` (|x| in [2⁻¹⁰²¹, DBL_MAX]:
  2⁻¹⁰²² excluded, since an exact value just below it rounds there at the subnormal grid but not at 53 bits), skipped
  inside an input region where every intermediate is proven normal (Tc, Pc, T, V within 2^±60, every physical state:
  `_km_pr_box` and inline siblings), which halved their cost.
- **Redlich–Kwong's Tc^2.5 with an unbounded exponent** is `_km_cr_pow(s, 2.5)·2^(5j)` for Tc = s·2^(2j), s in [1, 4)
  (`_km_rk_tc25_ux`): s^2.5 is in binary64's normal range, where the correctly rounded pow is RN₅₃, and RN₅₃ commutes
  with scaling by a power of two. No new pow primitive.
- **Unvalidated inputs** (Coulomb's charges): a non-finite one keeps Rust's expression.
- **The Peng–Robinson Z, roots and fugacity** keep `_km_pr_coefficients`, Rust's chain (the survey's F16, later).
- **The corpus** — the oracle's values kept; each module's B3 edge rows appended by `scripts/oracle/src/audit.rs` with
  their own seed (gas 1,680, potential 1,703, inorganic 254, nuclear 167); 1,157 #31 and 4,182 #32 rows noted and mapped
  (the generator's Redlich–Kwong notes use glibc's pow, which `validate.py` proves correctly rounded on every row; Tc^2.5
  midpoints in range, which glibc misrounds, are left to ADR 0006's corrected oracle and not added).
- **Performance** — the monomials 16–20 → 108–145 ns (`ideal_gas_pressure` 18 → 119, Rust 2.18), the Morse width
  33 → 135, the #20-style functions 35–182 → 39–187 on physical inputs (`peng_robinson_pressure` 78 → 100,
  `binding_energy_per_nucleon` 28 → 29, Rust 2.02); the slow paths 117–618 ns, the Morse width's exact root 1.27 µs.

Decided in the overflow audit's B4 (2026-11-11; ADR 0002 #33 and #34, B1's decisions unchanged; the survey's F7 and F8,
its proposed B6):
- **The treatment** — every function has a sum inside, so all fourteen get #20's semantics: Rust's operation sequence,
  its summation order included (`.sum()` from −0.0, thermochem's `total += coeff * prop` from +0.0), with an unbounded
  exponent, the final operation rounded once (the last addition, the difference, ½·Σ, × 931.494…, ADR 0002 #1's × 1000,
  the Shomate antiderivative's × 1000), comparisons on the unbounded values (`zero_order_concentration`'s clamp).
  Rust's NaN for ∞ − ∞ becomes the correctly rounded ±∞, or the finite value where the two sums cancel; Rust's checks,
  Err codes and order unchanged; a non-finite input keeps Rust's expression.
- **The pure-sum fast test** is "result finite": a sum of two doubles is exact below 2⁻¹⁰²¹ (both are multiples of
  2⁻¹⁰⁷⁴) and the 53-bit rounding above, so only overflow leaves the range, and ±∞ or NaN persists to the result
  (proven in the ADR, checked by `validate.py` partial sum by partial sum). Products inside take `_km_ux_normal` or an
  exact zero (inline on the hot paths).
- **The Shomate polynomials** test a box first (|T| in [2⁻³⁰, 2³⁰], the coefficients 0 or in [2⁻²⁰⁰, 2²⁰⁰]: every
  intermediate provably normal), then each operation; T = ±0 keeps Rust's chain (every product an exact zero, E/0 IEEE's
  special in both). The data-backed callers share the evaluation; every SHOMATE_DATA segment's coefficients are in the
  box, so the table-row path tests T alone, inline (`enthalpy_change_cp` and the flame keep their bits).
- **#29 corrected** — `_km_nernst_fast`'s bound is 2⁻¹⁰²¹ (`_km_ux_normal`'s), not 2⁻¹⁰²²: an IEEE result of exactly
  2⁻¹⁰²² can be an exact value just below it, rounded differently at 53 bits. B1–B3's other fast paths were audited for
  the same bound; none has it.
- **The benches split** (the user's decision of 2026-11-11, "Split the bench file") — `tests/audit.bcyr` holds the
  audit's Cyrius-only benches (B1–B4, typical inputs and edges), `tests/kimiya.bcyr` the rest; names unchanged; B3's
  back at the standard counts.
- **The corpus** — the oracle's values kept; each module's B4 edge rows appended by `scripts/oracle/src/audit.rs` with
  their own seed (thermo 184, reaction 172, nuclear 88, solution 430, inorganic 170, kinetics 171, thermochem 1,292;
  electrochemistry 23 for #29); 697 #33, 210 #34 and 17 more #29 rows noted and mapped.
- **Performance** — the fast tests cost 1–7 ns on physical inputs for the scalar functions and the pure sums
  (`gibbs_free_energy` 17 → 18, Rust 1.75; `limiting_molar_conductivity` 10 → 24, its fast path a function), about 9 ns
  per record for the reaction sums (`reaction_enthalpy_ch4` 259 → 293, Rust 38.58), the Shomate box 17 ns per
  `KmShomate.cp` (15 → 32, Rust 2.34); the table-row path tests T alone (`enthalpy_change_cp` 171 → 171, the flame 3.40
  → 3.29 µs); the slow paths 117–587 ns. The bench files take 177 s (`tests/kimiya.bcyr`) and 189 s
  (`tests/audit.bcyr`).

Decided in the overflow audit's B5 (2026-11-11; ADR 0002 #35, B1's decisions unchanged; the survey's F9, its proposed
B7):
- **The treatment** — `equilibrium_constant`'s exp and the three activity models' 10^x are no products or quotients of
  inputs (the round-once decision's scope), so their *arguments* get #20's edge fix, as `nernst_potential`'s chain did
  around its ln: Rust's sequence for the argument with an unbounded exponent, its final quotient or product rounded once
  into binary64 (`_km_keq_arg_ux`, `_km_dh_limiting_arg_ux`, `_km_dh_extended_arg_ux`, `_km_davies_arg_ux`), then the
  unchanged correctly rounded `_km_cr_exp` / `_km_cr_exp10`. No new primitive: exp / 10^x of the binary64 argument is
  the correctly rounded exp / 10^x of the 53-bit unbounded one (the ADR's lemma: the same double in the normal range,
  ±∞ against ±2¹⁰²⁴ or more, and 1 for both below 2⁻¹⁰²¹). `equilibrium_constant`'s ComputationError is judged on the
  port's K (a subnormal T can move a row across it either way); Rust's checks, Err codes and order unchanged.
- **The fast tests** (inline compares on Rust's own intermediates, proven in the ADR): R·T in [2⁻¹⁰²¹, DBL_MAX]
  (2⁻¹⁰²² excluded: T = 541,658,533,363,376·2⁻¹⁰⁷⁴ rounds R·T up to it from below); A·z² finite for the limiting law
  and Davies (a product with an integer factor is exact below 2⁻¹⁰²¹; √I is 0 or normal; Davies' bracket is always the
  unbounded one); the numerator and (B·å)·√I finite for the extended law (it proves Rust's γ: a subnormal numerator
  gives γ = 1 in both).
- **The corpus** — the oracle's values kept; none of the 2,627 earlier rows changes; edge rows appended by
  `scripts/oracle/src/audit.rs` with their own seed (reaction 676, solution 1,190); 845 #35 rows noted and mapped (326
  `equilibrium_constant`, 27 of them the port's ComputationError and 33 its K where Rust's overflows; 55
  `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity`). The notes' values are glibc's exp and
  pow(10, x) of the port's argument in the generator, each checked against MPFR by `validate.py`.
- **Performance** — the typical inputs keep their cost (`equilibrium_constant` 64 → 64 ns, the activity models 85–100
  ns); the slow paths 93–277 ns (`bench_audit_b5` in a new third bench file, `tests/audit2.bcyr`, 31 s: in
  `tests/audit.bcyr` they took it to 219 s of the 300 s deadline, so it splits as B4's did; `tests/audit.bcyr` 189 s).

Decided in the overflow audit's B6 (2026-11-11; ADR 0002 #36, #37 and #38, B1's decisions unchanged; the survey's F13,
F14 and F16's flame, its proposed B8):
- **The treatment** — `fft_spectrum` (#36, the user's "FFT: both"): #20's semantics over the whole transform, every
  butterfly and twiddle operation of hisab 1.4.0's `num::fft` with an unbounded exponent and each magnitude's √ rounded
  once (`_km_fft_ux`, `_km_fft_mag_ux`, `_km_fft_spectrum_ux`), a NaN or ±∞ sample keeping Rust's IEEE spectrum (B1's
  rule). `find_peaks` (#37, settled under the policy): the exact predicate w > min·(1 − t) + max·t, t the given double,
  min and max Rust's fold, the neighbour tests IEEE (`_km_peak_above_exact` on `src/mathx.cyr`'s new
  `_km_sign_scaled`: the exact sign of a sum of scaled doubles at any exponents). `adiabatic_flame_temperature` (#38, the
  user's "fix all three iterative edges", its flame part): the bracket's and the transcribed bisection's sign tests exact
  (proven equal to Rust's product tests except where the product underflows to a zero; the bisection's can't change a
  result); Newton's absolute 10⁻³ J/mol tolerance kept (at tiny scales it accepts the initial guess, Rust's own answer
  where its product does not underflow). Rust's checks, Err codes and order unchanged.
- **The fast tests** (proven in the ADR): for the FFT, operand bounds on Rust's own transform — every twiddle component 0
  or at least 2⁻¹⁰⁰ (checked in a stage's first block: the recurrence is input-independent), every butterfly's d 0 or at
  least 2⁻⁹²¹, so each product is an exact zero or its exact value at least 2⁻¹⁰²¹ (never a product's IEEE value, which
  can round up to 2⁻¹⁰²²); the retained outputs finite (an overflow anywhere reaches one as ±∞ or NaN); per bin the
  square sum in [2⁻⁹⁰⁰, DBL_MAX] or from re = im = 0, else that bin's unbounded √. For `find_peaks`, Rust's rounded
  threshold T decides when |RN(w − T)| exceeds B = RN(RN(|min|·2⁻⁵¹) + RN(|max|·2⁻⁵¹)) + 2⁻¹⁰⁷², above the proven error
  2.5u(1 + u)(|min| + |max|) + (1 + u)·2⁻¹⁰⁷⁴.
- **The work tables** — the IEEE transform runs in `out` only when cap ≥ 2n and `signal` does not overlap out's first 2n
  values (the slow path needs the signal), so an aliased signal now gets a 16·n-byte heap table; the slow path takes
  32·n heap bytes (the bump allocator never frees either).
- **The corpus** — the oracle's values kept; of the earlier rows 41 of 213 `fft_spectrum`, 3 of 994 `find_peaks` and 16
  of 1,937 flame rows change, the expected counts; edge rows appended with their own seeds (spectroscopy 2,393
  `find_peaks` and 515 `fft_spectrum`, thermochem 443 flame rows); 913 rows noted and mapped (418 #36, 222 #37, 273 #38),
  the two vector results with the port's whole result in the note.
- **Performance** — `fft_spectrum` +11–12 % (256: 27.3 → 30.5 µs; the per-butterfly checks about 6 %), `find_peaks_1000`
  19.0 → 20.3 µs, the flame unchanged (3.37 → 3.40 µs); the unbounded transform about 11× Rust's chain (318.1 µs at
  n = 256), a bin's unbounded √ about 0.1 µs; `bench_audit_b6` in `tests/audit2.bcyr`.

Decided in the overflow audit's B7 (2026-11-12; ADR 0002 #39, B1's decisions unchanged; the survey's F15, its proposed
B9):
- **The treatment** — `balance_equation` (#39, "Balancing: multi-word integers, no size limit", the survey's remedy (a),
  settled by the main session under the user's policy): Rust's algorithm with its i64 and i128 replaced by mathematical
  integers — the pivot search, the swap, reduce_row, the sign flip, the fraction-free row update then reduce_row, the
  lcm of the pivots, x, reduce_row(x), the zero, sign and u32 checks, the i128 verification, as written — so every
  overflow Err disappears and every other Err (ambiguous, only the trivial solution, a zero coefficient, mixed signs,
  past u32, the verification, each InvalidInput, the buffer check after) stays at its place. Rust's i64 code runs first,
  unchanged (`_km_balance_i64`: each overflow check returns a restart code), and at its first overflow the same steps
  run from the original matrix on signed multi-word integers (`_km_balance_big`); `_km_balance_equation(…, force)` skips
  the fast path for the tests and `cr_verify.py`.
- **The integers** — `_km_sb_*` in `src/stoichiometry.cyr` (its only user; `src/mathx.cyr`'s `_km_bn_*` are unsigned and
  fixed at 8 limbs): a signed length and 31-bit limbs, so no intermediate reaches 2⁶³; the schoolbook product, Knuth's
  Algorithm D (its add-back step included), the exact quotient, and Lehmer's gcd (Cohen's Algorithm 1.3.7 on the top 62
  bits, cofactors below 2³⁰, a unimodular combination that keeps the gcd whatever happens).
- **The size bound** (proven in the ADR) — every row after reduce_row is ± the primitive vector of a line spanned by
  minors of order at most min(m, n) (Schur's formula for the non-pivot rows, Cramer's rule for the pivot rows), so by
  Hadamard every stored entry, the lcm (a divisor of the pivot block's determinant) and x are below 2^B, B the smaller of
  the sums of the q largest row and column b = bitlen(max |entry|) + ⌈bitlen(nonzeros)/2⌉, and each product or
  difference before reduce_row below 2^(2B + 1): `_km_balance_bits` sizes the tables once; only allocation limits the
  size (a 4,096-word local table, the bump allocator past it; Err(ComputationError) when it fails).
- **Which overflows can become Ok** (proven) — at the end every pivot row is p_k at c_k and r_k at f, primitive, so x's
  gcd is 1 and |x| ≥ the lcm: an overflow of the lcm or of x is an Err past u32 exactly too, and a reduced row of gcd 2⁶³
  (all its nonzero entries −2⁶³) is orthogonal to no positive x; only the elimination's products, differences and the
  pivot row's negation can turn into Ok. The lcm's gcd conversion, reduce_row(x)'s gcd and the final flip are
  unreachable.
- **The corpus** — the oracle's values kept; of the earlier rows exactly the survey's 104 of 3,329 change (Err → Ok);
  521 edge rows appended with their own seed; 413 rows noted and mapped with the port's whole result.
- **Performance** — the fast path unchanged (`balance_ch4_combustion` 805 → 778 ns, `balance_kmno4_hcl` 2.24 →
  2.18 µs); the multi-word path where Rust gives its overflow Err: the survey's 2 × 3 3.6 µs, a 10 × 9 golden row 57 µs,
  a dense 30 × 31 of 2⁴⁵ entries 158 ms (products past 3,000 bits); `bench_audit_b7` in `tests/audit2.bcyr`.

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
| M4 Structure, nuclear & spectroscopy | ✅ **done 2026-11-10** (M4b organic, inorganic: ✅ done 2026-11-09; M4c nuclear: ✅ done 2026-11-09; M4a, the correctly rounded sin / cos / sincos / atan2 / acos: ✅ done 2026-11-09; M4d potential, the trig switch, ADR 0002 #13 completed and #16–#18: ✅ done 2026-11-09; M4e spectroscopy, ADR 0002 #19, hisab 1.4.0's FFT transcribed: ✅ done 2026-11-09; M4e-2 potential's geometry with an unbounded exponent, ADR 0002 #20: ✅ done 2026-11-10; M4e-3 the exact Rydberg rationals and the exact reduced mass, ADR 0002 #21, #22: ✅ done 2026-11-10; M4e-4 the cited AME2020 / TICE 2009 isotope table, ADR 0002 #23: ✅ done 2026-11-10; M4f stoichiometry, ADR 0002 #7 realised: ✅ done 2026-11-10; M4f-2 exact stoichiometry, ADR 0002 #24–#26: ✅ done 2026-11-10) | M4f-2: stoichiometry 47,776 rows (6,341 hard rows appended: `percent_yield`'s ties at RN, subnormal ones included, and range ends; the yield's over- and underflowing m/c and its midpoints; `limiting_reagent`'s exact ties, near-ties, both ratios over- or underflowing; `empirical_formula`'s knife edges at ±1/10, subnormal percentages and the 2³² edge), the oracle's values kept; ADR 0002 #24: 1,126 `percent_yield`, 240 `theoretical_yield_moles`, 598 `percent_composition` rows noted with the port's correctly rounded result; #25: 920 `limiting_reagent` rows with its exact index; #26: 270 `empirical_formula` rows with its exact result; all mapped: **336,584 ok, 0 failed, 0 measured, 0 pending, 32,153 accepted (ADR 6: 809; ADR 2 as M4f plus 3,154 #24–#26), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64. M4f: stoichiometry 41,435 rows (A: `balance_equation` 3,329 — the Rust tests, CH₄ / CuSO₄ + Fe / KMnO₄ + HCl in every row order, textbook and ionic equations, combustions, chains to 70 species, 820 seeded consistent and 700 random matrices, every error class, the 10⁻⁹ tolerance to the ulp, ±2⁵³, intermediate overflow; `assign_oxidation_states` 23,498; `limiting_reagent` 2,897; `empirical_formula` 3,383; B: `theoretical_yield_moles` 2,965, `percent_yield` 4,071, `percent_composition` 1,292); ADR 0002 #7: 9,090 oxidation rows noted with the port's exact rational states and 78 balancing rows with no integer-matrix counterpart, mapped: **333,397 ok, 0 failed, 0 measured, 0 pending, 28,999 accepted (ADR 6: 809; ADR 2 as M4e-4 plus 9,168 #7), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every float row 0 ulp. M4e-4: spectroscopy 66,691 rows unchanged in number, the oracle's values kept; ADR 0002 #23: 354 rows noted with the port's whole result and mapped (7 `isotope_abundances`, 347 `isotope_pattern`): **301,130 ok, 0 failed, 0 measured, 0 pending, 19,831 accepted (ADR 6: 809; ADR 2 as M4e-3 plus 354 #23), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference. M4e-3: spectroscopy 66,691 (13,706 new rows: every exact Bohr level, the 13 transition ties and the exact-value families of an exhaustive search, n₁ = n₂ − 1 at every scale, the series' large lines, the reduced mass's tie families, exponent gaps and range ends; ADR 0002 #21: 6,617 rows mapped (3,823 levels, 1,351 transition energies, 1,376 wavelengths, 67 series lines), #22: 2,837): **301,484 ok, 0 failed, 0 measured, 0 pending, 19,477 accepted (ADR 6: 809; ADR 2 as M4e-2 plus 6,617 #21 and 2,837 #22), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference. M4e-2: potential 31,761 (4,823 new edge rows; `distance_3d` 1,914, `bond_angle` 4,838, `dihedral_angle` 4,994; ADR 0002 #20: 3,605 rows mapped to the unbounded-exponent evaluation, 91 of them Err) and 3 corrected rows (`potential.cr.tsv`): **297,232 ok, 0 failed, 0 measured, 0 pending, 10,023 accepted (ADR 6: 809; ADR 2 as M4e plus 3,605 #20), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference; no row whose Rust chain stays in the normal range changes (`validate.py`). M4e: spectroscopy 52,985 (A 1,255: the isotope table at every u8, the series' n_lower, `find_peaks`; B 46,415: Beer–Lambert, the photon conversions, the Bohr levels at every u8 Z, transitions and the Rydberg formula, the series to the `checked_add` overflow, the wavenumbers (ADR 0002 #19: 1,654 rows mapped to the correctly rounded 1/(100·x)), the oscillator and rotor, `isotope_pattern` (every u8 Z, counts to 5,000 and u32::MAX); C 5,315: `transmittance`, `absorbance_from_transmittance`, `fft_spectrum` at every length to 16 and every power of two to 2¹³) and 118 corrected rows (`spectroscopy.cr.tsv`: glibc's pow(10, ·) and log10): **296,016 ok, 0 failed, 0 measured, 0 pending, 6,416 accepted (ADR 6: 807; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,810 N₀/2ⁿ + 208 secular + 59 Δ_tet + 1,041 r_min + 180 LJ + 181 Morse width + 1,654 1/(100x)), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference. M4d: potential 26,938 (B 11,320: Lennard-Jones (ADR 0002 #17: 180 rows mapped), the Morse width (#18: 181), Coulomb energy and force, distances; C 15,618: `lj_equilibrium_distance` (#16: 1,041 rows mapped to the correctly rounded σ·2^(1/6)), Morse (exp, near-cancellation rows), bond and dihedral angles (atan2) at Å, metre and femtometre scale, the tolerances to the ulp, the −π → π mapping) and 1 corrected row (`potential.cr.tsv`, glibc's atan2); the `math` trig rows at bound 0; ADR 0002 #13 correctly rounded (1,810 rows): **244,803 ok, 0 failed, 0 measured, 0 pending, 4,644 accepted (ADR 6: 689; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,810 N₀/2ⁿ + 208 secular + 59 Δ_tet + 1,041 r_min + 180 LJ + 181 Morse width), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference. M4c: nuclear 54,057 (A 1,898: ISOTOPES, `lookup_isotope`, `is_secular_equilibrium` (ADR 0002 #14: 208 rows mapped); B 35,370: the SEMF over every Z ≤ A ≤ 3Z with A ≤ 300 (ADR 0002 #3: 23 NaN rows mapped to Err), decay constants, mass defects, binding energies, Q-values; C 16,789: decay, half-lives (ADR 0002 #13: 1,812 rows mapped to the exact N₀/2ⁿ), `time_for_fraction`, the decay chain incl. λ₁ = λ₂, the threshold to the ulp and 554 subnormal-exp rows) and 23 corrected rows (`nuclear.cr.tsv`, glibc's expm1), plus inorganic's Δ_tet (ADR 0002 #15: 59 rows mapped, 55 new): **219,266 ok, 0 failed, 0 measured, 0 pending, 3,243 accepted (ADR 6: 688; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,812 N₀/2ⁿ + 208 secular + 59 Δ_tet, every mapped row by name but the Antoine ones), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every row 0 ulp from its reference. M4b: organic 9,104 (A 7,596: every enum value × property, `predict_geometry` over 7,491 (bonding, lone) pairs; B 1,508: `enthalpy_from_bonds`) + inorganic 15,107 (A 9,397: CFSE (exact, ADR 0002 #12: 13 rows mapped, accepted by name), ligands, `unpaired_electrons_oct`, IONIC_RADII, `lookup_ionic_radius` over every table symbol × i8; B 5,710: Δ_tet, moments, Born–Landé, Born–Haber) + the 10 last `constants` rows: **167,279 ok, 0 failed, 0 measured, 0 pending, 1,118 accepted (ADR 6: 665; ADR 2: 440 Antoine + 13 CFSE), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64; every M4b row 0 ulp from its reference | M4f-2: `tests/stoichiometry.tcyr`'s `cyrius_empirical_tolerance` becomes `cyrius_adr26_exact_tolerance`, with `cyrius_adr24_rounded_once` and `cyrius_adr25_exact_comparison` (no Rust test's assertion changes): 317 assertions; `cr_verify.py`: `muldiv` 640,000, `pyield` 397,518, `tyield` 340,006, `pcomp` 223,604, `limit` 400,008, `empir` 223,899 against `Fraction`, the exact paths forced, 0 mismatches on every target. M4f: 63 Rust tests ported (`tests/stoichiometry.tcyr`; 7 adapted for ADR 0002 #7; no integration test uses the module), plus 6 Cyrius-only groups (the rational states and their i32 range, the integer matrix and its 2⁵³ bound and heap tables, the checked i64 / i128 helpers at their edges, out-of-type words, slices and buffers, the empirical tolerance to the ulp): 265 assertions; `tests/kimiya.tcyr` 47; 21 test files, 0 warnings. M4e-4: `tests/spectroscopy.tcyr`'s `hydrogen_abundances_iupac_2009` adapted (the table's pairs; the Σ = 1 loop gains sulfur), `cyrius_adr23_table` and `cyrius_adr23_s8` (the table cell by cell, S₈'s seven peaks bit for bit against the Python replica, ³⁶S raising M+4 … M+6): 350 assertions; 20 test files, 0 warnings. M4e-3: `tests/spectroscopy.tcyr` gains `cyrius_adr21` and `cyrius_adr22` (exact values, every Z = n, the 13 ties, the reduced mass's ties and range ends, the fast paths against the exact ones on 6,000 seeded transitions and 20,000 seeded mass pairs): 311 assertions; `tests/mathx.tcyr` gains `cr_ratio_values` and `cr_ratio_agrees` (the multi-word integers, `_km_cr_ratio`'s ties and range ends, 10,000 seeded quotients against the double-double): 1,566 assertions; 20 test files, 0 warnings; `cr_verify.py`: `ratio` 320,000, `bohr` 260,429, `trans` and `ryd` 555,944 each, `series` 55,345, `rmass` 502,000 against `Fraction`, 0 mismatches on every target. M4e-2: `tests/potential.tcyr` gains 3 Cyrius-only groups (the edges against MPFR, the scale invariance on 600 random geometries × 75 scales bit for bit, the forced slow path against the fast one on 5,000 in-range geometries): 159 assertions; `tests/mathx.tcyr` gains `ux_values` (48 MPFR values) and `ux_ieee` (20,000 seeded pairs) and the specials: 1,546 assertions; 20 test files, 0 warnings; `cr_verify.py`: 1,040,039 arguments for the five final operations, 114,000 edge geometries, 0 mismatches on every target. M4e: 62 Rust tests ported + the 2 spectroscopy integration tests (`tests/spectroscopy.tcyr`, none adapted), plus 6 Cyrius-only groups (out-of-type integers, the caller buffers and slices, the isotope helpers' stable sort and thresholds to the ulp, the FFT transcription to 2¹⁶, #19's exact values, R_H): 289 assertions; `tests/potential.tcyr`'s `morse_width_parameter_matches_force_constant` calls `vibrational_frequency`; `tests/mathx.tcyr` gains `_km_cr_hundredth_div`'s groups (28 MPFR values, 4,000 seeded fast / exact pairs); 20 test files, 0 warnings. M4d: 33 Rust tests ported (`tests/potential.tcyr`; `morse_width_parameter_matches_force_constant` computes spectroscopy's √(k/μ) inline until M4e), plus 5 Cyrius-only groups (null pointers and the 24-byte point layout, the helpers' operation order, the tolerances to the ulp and the −π → π mapping, #16's exact r_min, #17 / #18's overflow edges): 126 assertions; `tests/mathx.tcyr` gains `_km_cr_exp2_scaled`'s groups (57 MPFR values, 4,000 seeded fast / accurate pairs) and loses `_km_acos` / `_km_cos`; `tests/nuclear.tcyr`'s #13 group asserts the correctly rounded values (243 assertions); `tests/gas.tcyr`'s primitive group calls `_km_cr_acos` / `_km_cr_cos`; 19 test files, 0 warnings. M4c: 36 Rust tests ported (1 adapted for ADR 0002 #3: `semi_empirical_invalid_is_nan` asserts Err; #13 / #14 change none), plus 8 Cyrius-only groups (the SEMF Result, out-of-u32 words and the checked add, the chain buffer and threshold, the q_value slices, the isotope handles and UTF-8 modes, the typed constants and Rust's edges, ADR 0002 #13's exact N₀/2ⁿ, #14's exact comparison): 236 assertions; inorganic's #15 group (332 assertions there); 18 test files, 0 warnings. M4b: 51 Rust tests ported (organic 21, inorganic 30; 5 organic tests adapted for ADR 0002 #6's names) + the 2 organic integration tests, plus 10 Cyrius-only groups (ADR 0002 #6's acceptor reference tests, invalid enum values and handles, out-of-range u8 / u32 / i8 / bool words, the bond records, the exact CFSE): 464 assertions; 17 test files, 0 warnings | M4f-2: `stoichiometry/percent_yield` 14 → 81 ns, `theoretical_yield_moles` 46 → 115, `percent_composition_glucose` 72 → 276, `empirical_formula_glucose` 217 → 586, `limiting_reagent` 32 → 33; the exact paths `theoretical_yield_moles_tie` 133 ns, `percent_yield_subnormal` 2.15 µs, `limiting_reagent_tie` 94 ns, `empirical_formula_edge` 7.46 µs. M4f: `stoichiometry/balance_ch4_combustion` 1.20 µs → 770 ns with the proven fast paths (Rust 233.90 ns), `oxidation_states_kmno4` 156 → 122 ns (Rust 22.38); Cyrius-only `balance_kmno4_hcl` 3.84 → 2.16 µs, `oxidation_states_fe3o4` 132 → 110 ns, `limiting_reagent` 31, `theoretical_yield_moles` 46, `percent_yield` 14, `percent_composition_glucose` 71, `empirical_formula_glucose` 217. M4e-4: `spectroscopy/isotope_pattern_glucose` 15.8 µs either way (Rust 2.80); Cyrius-only `isotope_pattern_cystine` 16.9 → 18.5 µs, `isotope_pattern_s8` 8.4 → 18.8 µs. M4e-3: `spectroscopy/bohr_energy_level` 9 → 49 ns (Rust 1.21); Cyrius-only `bohr_transition_energy` 27 → 73, `rydberg_wavelength` 24 → 73, `series_wavelength_balmer_alpha` 28 → 78, `reduced_mass_hcl` 14 → 73; the exact path `bohr_transition_energy_tie` 2.42 µs, `reduced_mass_tie` 2.20 µs. M4e-2: `potential/dihedral_angle` 131 → 133 ns (Rust 20.35), `potential/bond_angle` 213 → 213, `potential/distance_3d` 14 → 16; the slow path `distance_3d_2p600` 314 ns, `bond_angle_2p600` 1.13 µs, `dihedral_angle_2p600` 1.31 µs. M4e: `spectroscopy/absorbance` 11 ns (Rust 0.98), `spectroscopy/bohr_energy_level` 16 ns (Rust 1.21), `spectroscopy/isotope_pattern_glucose` 16.8 µs (Rust 2.80, 6.0×); Cyrius-only `fft_spectrum_256` / `_1024` / `_4096` 28.2 / 134 / 628 µs, `find_peaks_1000` 22.5 µs, `wavenumber_from_wavelength` 58 ns. M4d: `potential/dihedral_angle` 130 ns (Rust 20.35, 6.4×); Cyrius-only `potential/bond_angle` 215 ns, `potential/lennard_jones` 35, `potential/lj_equilibrium_distance` 109, `potential/morse` 66, `mathx/cr_exp2_scaled` 151; the trig benches lose `_km_cos` (1.68 µs) and `_km_acos` (5.8 µs). M4c: `nuclear/binding_energy_per_nucleon_fe56` 27 ns (Rust 2.02, 13×), `nuclear/semi_empirical_fe56` 131 ns (Rust 18.50, 7.1×), `nuclear/decay_chain_two_step` 256 ns (Cyrius-only). M4b: `organic/enthalpy_from_bonds` 49.5 ns, `inorganic/lookup_ionic_radius` 107 ns (no Rust baseline: rust-old/benches has neither) | M4f-2: ADR 0002 #24–#26 (its M4f-2 addition quotes the decisions: `_km_cr_mul_div` with its midpoint resolution, `_km_ratio_lt`, the exact empirical formula and its 2³² cut); `exact::cr_mul_div`, `ratio_lt`, `empirical_exact` and `m_stoichiometry::emit_exact` in the generator; `validate.py` §18's proofs of the five; the harness's #24–#26 mapping; settled points in §6.4. M4f: ADR 0002 #7 realised (its M4f addition: the settled points and the five float functions' measured candidates); `m_stoichiometry.rs` with `port_oxidation` and `no_counterpart`; `validate.py` §18; the harness's `_par_in_matrix` and `_par_map_toks(adr, tag)`; surface parity 401/401; settled points in §6.4. M4e-4: ADR 0002 #23 (the cited source file `scripts/oracle/isotopes.tsv`, the generator's checks against the oracle, `pattern_replica`, the note's whole-result token list and the harness's `_par_map_toks`, `isotope_sources.py`); settled points in §6.4. M4e-3: ADR 0002 #21 and #22 (`src/mathx.cyr`'s multi-word integers `_km_bn_*`, `_km_cr_ratio`, `_km_dd_mul`, `_km_dd_div_cr`; spectroscopy's fast and exact paths); `exact::cr_ratio` and the transition search in the generator, `cr_verify.py`'s `ratio` … `rmass`; the harness's #21 / #22 mapping; settled points in §6.4. M4e-2: ADR 0002 #20 (`src/mathx.cyr`'s unbounded-exponent helpers `_km_ux_*`, the geometry's fast-path tests and `_km_*_ux` slow paths); `exact::Ux` / `ux_atan2_f64`, `mpfr_ref.UX`, `cr_verify.py`'s `uxsqrt` … `dih`; the harness's `_par_got_noted20`; settled points in §6.4. M4e: generated `spectroscopy_data` (`KmSpectralSeries` and n_lower, the isotope table from every u8 z, R_H through the API); ADR 0002 #19 (`_km_cr_hundredth_div`, `cr_verify.py`'s `h100`); ADR 0003 amendment (hisab 1.4.0's FFT transcribed with correctly rounded twiddles); the harness's `_PAR_MAX` 8448 and `_par_force_map`; surface parity 394/394; settled points in §6.4. M4d: no data file (COULOMB_CONSTANT is `src/constants.cyr`'s); 2^(1/6) as a triple-double generated into `mathx_data` with its continued-fraction proof; ADR 0002 #13 correctly rounded (`_km_cr_exp2_scaled`, `cr_verify.py`'s new `exp2s`), the new #16, #17, #18; gas on M4a's trig, ADR 0005's double-double primitives removed; the generator's exact arithmetic in `exact.rs`; surface parity 357/357; settled points in §6.4. M4c: generated `nuclear_data` (18 ISOTOPES rows, the first `u16` data column, a UTF-8 decay-mode blob); ADR 0002 #3 (SEMF → `Result`) and, by the user's decision on the inexact-formula review, #13 (the exact N₀/2ⁿ), #14 (the exact secular test), #15 (Δ_tet past Δ·4's overflow); the harness reads a mapped row's exact value from its note; surface parity 347/347; settled points in §6.4. M4b: generated `organic_data` (3 enums, the first generated enums with property columns) and `inorganic_data` (KmLigand, the private d-electron tables through the API, the 10 constants, 30 ionic radii); ADR 0002 #6 (`KM_GEOMETRY_T_SHAPED`, `…_is_hydrogen_bond_donor`, the new acceptor by IUPAC 2011's criterion, weak acceptors included) and the new #12 (exact CFSE); surface parity 323/323; settled points in §6.4 |
| Overflow audit (after M4) | **B1 ✅ done 2026-11-10** (F1: thermo / phase products, ADR 0002 #27; F4: kinetics rationals, #28); **B2 ✅ done 2026-11-10** (F3: solution / electrochemistry, #29; F5: spectroscopy, #30; the survey's B2 and B3); **B3 ✅ done 2026-11-11** (F2: gas laws, #31; F6: potentials, Born–Landé, BE per nucleon, #32, folding #17 and #18; the survey's B4 and B5); **B4 ✅ done 2026-11-11** (F7: sums with cancellation, #33; F8: the Shomate polynomials, #34; the survey's B6; #29 corrected; the benches split); **B5 ✅ done 2026-11-11** (F9: the exp / exp10 arguments, #35; the survey's B7); **B6 ✅ done 2026-11-11** (F13: the FFT, #36; F14: `find_peaks`, #37; F16's flame, #38; the survey's B8); **B7 ✅ done 2026-11-12** (F15: `balance_equation` on unbounded integers, #39; the survey's B9); the survey's B10–B12 ☐ | thermo 1,644 (+537), phase 4,174 (+528), kinetics 31,210 (+1,064), the oracle's values kept; #27: 151 `heat_transfer`, 171 `specific_heat_from_calorimetry`, 82 + 82 colligative, 196 `osmotic_pressure`; #28: 89 `second_order_concentration`, 90 `zero_order_half_life`, 118 `steady_state_intermediate`, 179 `pre_equilibrium_rate`, 160 `michaelis_menten` rows noted and mapped; B2: solution 29,024 (+458), electrochemistry 19,397 (+823), spectroscopy 67,869 (+1,178), the oracle's values kept; #29: 94 `dilution`, 148 `solution_conductivity`, 904 `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 119 `nernst_potential`; #30: 195 `absorbance`, 219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524 `rotational_constant`, 118 `diatomic_moment_of_inertia` rows noted and mapped: 334,992 ok, 0 failed, 0 measured, 0 pending, 38,333 accepted (ADR 6: 809), 0 glibc-matching; B3: gas 22,773 (+1,680), potential 33,464 (+1,703), inorganic 15,416 (+254), nuclear 54,224 (+167), the oracle's values kept; #31: 146 `ideal_gas_pressure`, 148 `ideal_gas_volume`, 222 `gas_density`, 206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234 `redlich_kwong_pressure` (98 the port's Err); #32: 582 `lennard_jones`, 1,860 `morse_width_parameter`, 730 `coulomb`, 679 `coulomb_force`, 175 `born_lande_lattice_energy`, 156 `binding_energy_per_nucleon` rows noted and mapped (#17's 180 and #18's 181 among them): 333,818 ok, 0 failed, 0 measured, 0 pending, 43,311 accepted (ADR 6: 809), 0 glibc-matching; B4: thermo 1,828 (+184), reaction 1,730 (+172), nuclear 54,312 (+88), solution 29,454 (+430), inorganic 15,586 (+170), kinetics 31,381 (+171), thermochem 11,149 (+1,292), electrochemistry 19,420 (+23), the oracle's values kept; #33: 28 `hess_law`, 27 `enthalpy_from_formation`, 115 `reaction_enthalpy`, 113 `reaction_gibbs_energy`, 66 `reaction_entropy`, 58 `q_value`, 13 `ionic_strength`, 74 `born_haber_lattice_energy`, 57 `gibbs_free_energy`, 54 `kohlrausch_molar_conductivity`, 31 `limiting_molar_conductivity`, 61 `zero_order_concentration`; #34: 185 `cp`, 25 `enthalpy_antiderivative`; #29: 17 more `nernst_potential` rows noted and mapped: 335,424 ok, 0 failed, 0 measured, 0 pending, 44,235 accepted (ADR 6: 809), 0 glibc-matching; B5: reaction 2,406 (+676), solution 30,644 (+1,190), the oracle's values kept; #35: 326 `equilibrium_constant` (27 the port's ComputationError, 33 its K where Rust's overflows), 55 `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity` rows noted and mapped: 336,445 ok, 0 failed, 0 measured, 0 pending, 45,080 accepted (ADR 6: 809), 0 glibc-matching; B6: spectroscopy 70,777 (+2,908), thermochem 11,592 (+443), the oracle's values kept; #36: 418 `fft_spectrum`, #37: 222 `find_peaks`, #38: 273 `adiabatic_flame_temperature` rows noted and mapped: 338,883 ok, 0 failed, 0 measured, 0 pending, 45,993 accepted (ADR 6: 809), 0 glibc-matching; B7: stoichiometry 48,297 (+521), the oracle's values kept; #39: 413 `balance_equation` rows noted and mapped (104 of the 3,329 earlier ones): **338,991 ok, 0 failed, 0 measured, 0 pending, 46,406 accepted (ADR 6: 809), 0 glibc-matching** — identical on x86_64 ELF, x86_64 PE, aarch64 | `tests/mathx.tcyr` `cr_mono_values`, `cr_mono_sign_and_tie`, `cr_mono_agrees` (1,588 assertions); `cyrius_adr27_rounded_once` (thermo 36, phase 126) and `cyrius_adr28` (kinetics 153); no Rust test changes; `cr_verify.py`: `mono` 430,010, `monotie` 240,000 and the ten functions' families against `Fraction` and `mpfr_ref.UX`, 0 mismatches on every target; B2: `b2_values`, `b2_agrees` (1,610 assertions), `cyrius_adr29` (solution 144, electrochemistry 1,409) and spectroscopy's #30 group (371), no Rust test changes, `cr_verify.py` `sqrtdiv` 326,845 and the twelve functions' families (3,085,916 in all), 0 mismatches on every target; B3: `b3_values`, `b3_agrees` (1,626 assertions), `cyrius_adr31` (gas 181), `cyrius_adr32` (potential 192), `cyrius_adr32_born_lande` (inorganic 337), `cyrius_adr32_bea` (nuclear 247), no Rust test changes, `cr_verify.py` `wsqrt` 316,971 and the twelve functions' families (2,815,142 in all), 0 mismatches on every target; B4: `cyrius_adr33_*` (thermo 44, reaction 34, nuclear 252, solution 151, inorganic 342, kinetics 159), `cyrius_adr33_reaction_sums` and `cyrius_adr34_shomate` (thermochem 415), `cyrius_adr29_2m1022` (electrochemistry 1,413), no Rust test changes, `cr_verify.py` the twelve families (1,250,069) and `nernst`'s 2,800 products at 2⁻¹⁰²², 0 mismatches on every target; B5: `cyrius_adr35_keq` (reaction 49), `cyrius_adr35_activity` (solution 165), no Rust test changes, `cr_verify.py` `keq`, `dhl`, `dhe`, `dav` (505,982), 0 mismatches on every target; B6: `cyrius_adr38_flame_sign` (thermochem 421), the #36 and #37 groups (spectroscopy 403), no Rust test changes, `cr_verify.py` `flame`, `pk3` … `pk16`, `fft1` … `fft256`, `signsc` (644,838), 0 mismatches on every target; B7: `cyrius_adr39_unbounded_integers`, `cyrius_adr39_size_and_forced`, `cyrius_adr39_multiword` (stoichiometry 350), no Rust test changes, `cr_verify.py` `bal23` … `bal2425` (54,306 matrices), 0 mismatches on every target | `heat_transfer` 5.4 → 65 ns (Rust 0.95), `osmotic_pressure` 20 → 139, `pre_equilibrium_rate` 27 → 125, `michaelis_menten` 17 → 24 (Rust 2.92); stoichiometry as in M4f-2; B2: `absorbance` 5 → 68 ns (Rust 0.98), `nernst_potential` 83 → 92 (Rust 6.42), the monomials 6–33 → 72–148, `vibrational_frequency` 18 → 87; B3: `ideal_gas_pressure` 18 → 119 ns (Rust 2.18), `binding_energy_per_nucleon_fe56` 28 → 29 (Rust 2.02), the monomials 16–20 → 108–145, `morse_width_parameter` 33 → 135, `peng_robinson_pressure` 78 → 100; B4: `gibbs_free_energy` 17 → 18 ns (Rust 1.75), `reaction_enthalpy_ch4` 259 → 293 (Rust 38.58), `shomate_cp` 15 → 32 (Rust 2.34), `enthalpy_change_cp` 171 → 171, the flame 3.40 → 3.29 µs, the slow paths 117–587 ns; the bench files 177 s and 189 s; B5: `equilibrium_constant` 64 → 64 ns, `debye_huckel_extended` 105 → 100, the slow paths 93–277 ns, `tests/audit2.bcyr` 31 s (`tests/audit.bcyr` 189 s); B6: `fft_spectrum_256` 27.3 → 30.5 µs (+12 %), `find_peaks_1000` 19.0 → 20.3, the flame 3.37 → 3.40, the unbounded transform 318.1 µs at n = 256; B7: `balance_ch4_combustion` 805 → 778 ns (Rust 233.90), the multi-word path 3.6 µs (2 × 3) … 158 ms (30 × 31) | ADR 0002 #27, #28 (its B1 addition quotes the decisions); `src/mathx.cyr`'s correctly rounded monomial `_km_cr_mono` with straight-line fast paths per shape; `scripts/oracle/src/audit.rs`; `validate.py` §19; settled points in §6.4; B2: ADR 0002 #29, #30, four more shapes and the correctly rounded √ of a quotient `_km_cr_sqrt_div`, `validate.py` §20; B3: ADR 0002 #31, #32 (#17, #18 folded in), the correctly rounded w·√((x/y)·2^e) `_km_cr_mul_sqrt_ratio`, `_km_ux_normal`, `validate.py` §21; B4: ADR 0002 #33, #34, #29's correction, `mpfr_ref.ux_sum`, `validate.py` §22, `tests/audit.bcyr`; B5: ADR 0002 #35 and its lemma, `mpfr_ref.ux_exp_final` / `ux_exp10_final`, `validate.py` §23, `tests/audit2.bcyr`; B6: ADR 0002 #36–#38, ADR 0003's B6 amendment, `_km_sign_scaled`, `exact::peak_above`, `validate.py` §24, `cr_verify.py`'s vector-result families; B7: ADR 0002 #39 and its proofs, `_km_sb_*`, `_km_balance_bits`, `exact::Int`, `audit::port_balance`, `validate.py` §25 |
| M5 Numerics | ☐ | | | | ADR 0003 decides the hisab gaps |
| M7 Release 2.0.0 | ☐ | | | | |
