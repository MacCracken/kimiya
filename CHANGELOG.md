# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.1] - 2026-11-13

The repairs of 2.0.0's P(-1) review sweep. No golden row changes: the corpus is still 342,485 rows ok, 0 failed and
46,634 accepted. Every performance change is bit-identical in results and Err codes, checked by the corpus and by
randomized probes against the 2.0.0 code. The one API addition is `KM_PI`.

### Fixed
- **Null caller buffers** (kinetics, nuclear, gas, fitting, sensitivity). `lineweaver_burk`, `decay_chain_two_step`,
  `peng_robinson_z_roots`, `fit_arrhenius`, `equilibrium_constant_gradient`, `arrhenius_gradient` and `nernst_gradient`
  wrote through a null `out` and faulted. They now return `Err(KM_ERR_INVALID_INPUT)`, checked after Rust's own errors.
  For `peng_robinson_z_roots`, Rust's no-physical-root ComputationError still comes first. The functions that take a
  capacity already treated a null `out` as holding nothing. The port plan's note that a null `out` is the caller's
  contract is updated.
- **Aliased buffers** (spectroscopy, reaction_dynamics, stoichiometry, network, sensitivity). Rust can't alias a slice
  with its output; here each is a raw pointer. Every public function that reads one caller buffer and writes another
  was checked, by reading and by probes of each overlap: `out` = an input, `out` partly inside one, an input partly
  inside `out`, two outputs overlapping. Where a write could reach a value still to be read, the overlap is now
  `Err(KM_ERR_INVALID_INPUT)`, after Rust's own errors; exact in place works wherever the loop allows it at no cost.
  Each function's doc comment gives its rule:
  - spectroscopy `find_peaks`: an `out` inside `data` past its first value (data < out < data + 8n) overwrote values
    still to be scanned. It is now `Err(KM_ERR_INVALID_INPUT)`, checked after Rust's three Ok(0) cases (n < 3, a
    threshold outside [0, 1], no finite value), which write nothing. `out = data` and `out < data` work.
  - reaction_dynamics `simulate_stiff`: in place (`out = y0`), point 0's t overwrote y0[0] before the copy. Point 1 is
    now written first, then point 0's state from its last value down, then its t. Any other overlap with y0 is
    `Err(KM_ERR_INVALID_INPUT)`, after the integration.
  - stoichiometry `assign_oxidation_states`: in place (`out = atoms`), it read records it had already overwritten: wrong
    states, or a fault on a garbage Z. The records are now written from the last down, so in place and any `out` above
    `atoms` work. An `out` below `atoms` that overlaps them is `Err(KM_ERR_INVALID_INPUT)`.
  - stoichiometry `empirical_formula` and `percent_composition`: an `out` inside the input records past their start
    (for `percent_composition`, the molecule's own atom records, its vec's buffer) overwrote atomic numbers still to be
    read: wrong records, or for `percent_composition` a fault on a garbage Z. It is now `Err(KM_ERR_INVALID_INPUT)`. In
    place and an `out` below the records work: record k is written after input record k is read, and below record
    k + 1.
  - network `stoichiometric_matrix` zeroes the matrix before it reads the reaction records, so a matrix overlapping the
    records or any term list, in place included, gave a wrong matrix or faulted on a garbage term pointer. Every such
    overlap is now `Err(KM_ERR_INVALID_INPUT)`. The test is folded into the records' validation pass, and the fill's
    loops now walk pointers (the same stores in the same order, checked against the old fill on 200,000 random
    networks), which pays for it: `network/stoichiometric_matrix_h2_br2` 268 → 252 ns (the checks alone made it 33 ns
    slower in A/B probes, and a call to the helper per span 96 ns).
  - network `rate_from_stoichiometry` reads every rate for each species, and row i + 1 after writing species i's rate:
    an `out` overlapping the rates, or starting inside the matrix past its first entry, read values it had overwritten.
    Both are now `Err(KM_ERR_INVALID_INPUT)`. `out = matrix` and an `out` below it work.
  - sensitivity `compute_gradient`: a tape built through hisab's setters can hold the caller's arrays. `parameters` over
    the tape's header or nodes read values the nodes had overwritten; `out` over the header faulted, and over the nodes,
    or over the adjoints other than exactly, left the tape's records overwritten. These are now
    `Err(KM_ERR_INVALID_INPUT)`, before the tape is reset. `out` = the adjoint array works (the n variables are nodes
    0 … n − 1, whose adjoints are the gradient), and `out` may overlap `parameters`. The check's three loads and five
    span tests cost about 10 ns: `sensitivity/compute_gradient_arrhenius` 500 → 508–516 ns (+2–3 %).
  - Checked and unchanged: `fit_arrhenius`, `fit_polynomial`, `balance_equation`, `isotope_pattern`, `fft_spectrum`,
    `simulate_adaptive` and `simulate_stochastic` read every input before they write `out` (or, for `fft_spectrum`,
    transform in a pool table when the signal lies in out's first 2n values), so every overlap gives the separate
    buffer's result. `peng_robinson_z_roots` takes its parameters by value. Callback contexts are the caller's, and
    kimiya never reads them.
  - The overlap test is src/validate.cyr's `_km_overlaps(a, a_bytes, b, b_bytes)`, which forms no pointer + size and so
    can't wrap. `stoichiometric_matrix` and `compute_gradient` inline its signed form, where no end can wrap (a call
    costs about 5 ns). Cyrius-only test groups check each rule: in place gives the separate buffer's result, and each
    refused overlap gives its Err with the inputs untouched.
  - The timings are `cyrius bench` runs of the tree before and after on the same host, two of each for
    tests/numerics.bcyr. The other benches on a changed path are unchanged: `percent_composition_glucose` 266 → 266 ns,
    `empirical_formula_glucose` 484 → 488 ns, `empirical_formula_edge` 3.90 → 3.92 µs, `find_peaks_1000`
    11.31 → 11.34 µs.
- stoichiometry `empirical_formula`: its 16·n heap table was taken before Rust's per-element validation, so a huge n
  gave ComputationError instead of Rust's InvalidElement or InvalidInput. The table is now taken after the validation,
  and its size is guarded.
- **Memory growth per call** (validate, numerics, reaction_dynamics, fitting, network, sensitivity, spectroscopy,
  stoichiometry). A work table larger than its function's stack table came from the bump allocator, which never frees.
  A call in a loop therefore grew the heap by the table every time, until the allocator's 2 GiB limit turned it into a
  spurious `Err(KM_ERR_COMPUTATION_ERROR)`. Examples: 64 KiB per `fft_spectrum(4096)` with cap n/2 + 1, 12 KiB per
  `isotope_pattern` of Br₃₀₀, 2 KB per `simulate_stiff` with n = 10, and 4.9 KB per `fit_polynomial` of 200 points. The
  tables now come from a pool of eight reusable, grow-only buffers in `src/validate.cyr` (`_km_scratch_take` /
  `_km_scratch_give`; [ADR 0007](docs/adr/0007-shared-work-table-pool.md)):
  - A call claims a slot by an atomic compare-and-swap, so neither two threads nor a re-entrant call (a callback into
    kimiya) ever share one. A slot grows to at least twice its size, and the call gives it back before it returns.
  - After the first call at a given size, a call allocates nothing. With every slot claimed, a fresh bump block serves
    the call, as in 2.0.0.
  - Every site is covered: the transcribed routines' tables (`_km_num_table`), `simulate_stiff`, `simulate_stochastic`,
    `fit_arrhenius`, `fit_polynomial`, `sparse_stoichiometric_matrix`, `compute_gradient`, `isotope_pattern`,
    `fft_spectrum` and its unbounded slow path, `balance_equation`'s three tables and `empirical_formula`'s.
  - The only per-call allocations left are results the caller keeps: KmMolecule's vec and formula string, and the CSR
    handle `sparse_stoichiometric_matrix` returns.
  - Cyrius-only `cyrius_201_no_growth` groups in seven test files check that a loop allocates nothing after its first
    call, on Ok and on every Err after the table, and that every slot is given back. `tests/parity.tcyr` checks that no
    slot is still claimed after the whole corpus.
  - There is no API or result change.
- reaction_dynamics `simulate_stochastic`: `pcg32_new` took 16 bytes from the bump allocator on every call. The
  generator now lives on the stack, seeded exactly as `pcg32_new` seeds it (the same stream, which a test pins).
- **Latent guards** (mathx, numerics). No caller reached any of these today:
  - `_km_cr_finish` and `_km_exp_fin` gave −0 for an exponent above 1024, because 2^(e − 1) wrapped. They also gave +0
    for a normal result at e = −1023, where 2⁻¹⁰²³ encoded as +0. They now give +∞, 2⁻¹⁰²²·(h + m + l)/2 (scaled in
    two exact steps) and, below e = −2044, +0.
  - `_km_f64_int(+0)` looped forever; it now returns (0, 0).
  - `_km_num_table` took a size that had overflowed to a negative word as fitting its stack table; it now fails.
  - The least squares guards n·m + n·n before multiplying it out.
- Stack tables are sized by their constants (`_KM_NUM_STACK`, `_KM_RD_STIFF_STACK`, `_KM_FLAME_CACHED * 2`, …) instead
  of a second literal.

### Changed
- validate: the shared helpers live in `src/validate.cyr`. They are the slice test `_km_is_slice` (formerly
  reaction_dynamics' `_km_rd_slice` and fitting's `_km_fit_slice`, now also used by `isotope_pattern`, `find_peaks` and
  `empirical_formula`) and the Rust integer-type predicates `_km_is_i32`, `_km_is_u8`, `_km_is_u32`, `_km_is_bool` and
  `_km_is_enum_value`. The predicates came from solution, inorganic, nuclear and organic, which spectroscopy had reached
  across for. The buffer-overlap test `_km_overlaps` is new (Fixed, aliased buffers).
- mathx:
  - one special-value ladder, `_km_ln_special`, for ln's and log10's fast and accurate paths;
  - `_km_cr_log10` takes the layout of `_km_cr_ln`: the normal range inline, the rest through
    `_km_cr_log10_fast` and the accurate path;
  - `_km_mul3_rest`, the shared tail of the five x·y·z shapes (`heat_transfer`, `boiling_point_elevation`,
    `freezing_point_depression`, `absorbance`, `diatomic_moment_of_inertia`), goes straight to `_km_cr_mono` once the
    fast path declines (2.0.0 retried the fast path through `_km_cr_mul3`);
  - `_km_ux_sqrt_f64`'s unreachable r ≥ 2 branch is marked defensive, with the proof.
- solution: `_km_dh_az2_ux`, Rust's ((−A)·z)·z with an unbounded exponent, is now the prologue of the three
  activity-coefficient `_ux` helpers.
- error: `km_error_name` and `km_error_message` are now a `match` with an `_ =>` arm, CLAUDE.md's enum rule.
- constants: **`KM_PI`** (Rust's `std::f64::consts::PI`) is generated into `src/constants.cyr` by `scripts/oracle`.
  `t_constants.rs` emits Rust's standard-library constants the port uses, and `validate.py` and `probe_tables.py` check
  each against Python's math value; it has no golden row. It replaces π's bits, typed by hand at six sites in
  spectroscopy and potential. This is the new public constant.
- Comments corrected:
  - spectroscopy: 2⁻⁹²¹ ≈ 5.64·10⁻²⁷⁸;
  - nuclear: only half-lives below ln 2/DBL_MAX ≈ 3.86·10⁻³⁰⁹ overflow;
  - electrochemistry: 1 ..= 8 electrons;
  - inorganic: four functions differ from Rust (ADR 0002 #12, #15, #32, #33), and the Born–Landé bound is ≈ 2¹¹²;
  - numerics: dopri45 takes 9·n words, so it is on the stack for n ≤ 10;
  - molecule: 2²⁶ is a quarter of VEC_CAP_MAX;
  - gas: van der Waals' quotient is at least 2⁻²³⁸;
  - mathx: `_km_exp_core` is valid for |xh| < 1180, by the K < 2²⁰ argument, and `_km_ln_core`'s call costs about
    14 ns.
- Not done, measured:
  - `_km_shomate_cp5` / `_km_shomate_h5` calling `_km_shomate_cp_fast` / `_h_fast`: the extra call costs 7 ns on
    `KmShomate.cp` (34 → 41 ns) and 3 ns on `enthalpy_antiderivative`.
  - `_km_exp_acc_core` ending with `_km_exp_acc_r`: the code moved so that `ideal_gas_pressure`,
    `faraday_charge_required` and `pre_equilibrium_rate` ran 6–7 % slower, an alignment effect, in A/B runs against
    2.0.0.
  - One shared gas box test: the three tests differ, and a call costs about 5 ns on paths of 120–300 ns.
  - `_km_mul5_fast` with its helpers inline: `osmotic_pressure` was slower (139 → 147 ns). `_km_mul4_fast`'s inlining
    measured a win and is in.

### Performance
All bit-identical. Before is the 2.0.0 tree's `cyrius bench` run of 2026-11-13, taken before any change; after is 2.0.1's
run on the same host. Run-to-run noise on this host is about ±5 %, and this 2.0.1 run was a slow one: untouched stdlib
benches (`f64_sin`, `f64_cos`, `ganita_f64_atan2`) ran 5–8 % slower than in the 2.0.0 run. Code-layout shifts move
untouched functions by as much, so a change under about 8 % is within the noise. For those, the A/B probes quoted
above (built from both trees, run back to back) are the better measure. Rust is
`rust-old/benchmarks.md`'s criterion median, for the benchmarks that have a Rust counterpart.

- **Name lookups** (thermochem, phase, element). Each matching first byte used to cost a strlen and a memeq. The scans
  now compare bytes from byte 1 to the first difference or the shared NUL. `lookup_antoine` used to run a `streq` on
  every row; it now has the same first-byte test, so chloroform takes 393 → 43 ns (an A/B probe). A shared compare
  helper measured 11 ns slower per lookup, so each scan keeps its own inline loop.
- **adiabatic_flame_temperature** (thermochem). Once per call, a prefix table stores each product's running total
  before each segment and the antiderivative at each segment's lower end. f(T) for T ≥ T₀ then adds one antiderivative
  per product (the loop's additions, in its order). T < T₀, NaN, or a product past the table runs the full segment
  loop. f(T₀) is −ΔH exactly, so that evaluation is skipped.
- `empirical_formula` (stoichiometry) keeps the counts from the search pass instead of re-running the tolerance test.
- `molecular_weight` (molecule) inlines `element_lookup_by_number`'s two checks.
- The least squares (numerics) no longer zeroes its Q and R table: Q is written whole, and only R's written entries are
  read. Columns 0 and 1 are 1.0 and x, without `_km_powi`. `_km_num_point` inlines `_km_num_copy`.
- **The FFT** (spectroscopy). `_km_fft` and `_km_fft_ux` run the twiddle loop outside the block loop, so each twiddle is
  computed and tested once a stage. The data's range test works on the bits.
- spectroscopy: `find_peaks`' min/max pass, the isotope sort's key and the isotope normalisation's max are inline.
- nuclear: `decay_chain_two_step` reuses e^(−λ₁t) when λ_min is λ₁.
- potential: `bond_angle`'s and `dihedral_angle`'s fast paths inline sub3, cross, dot3 and max.
- electrochemistry: `nernst_potential`'s guards and fast-path test are inline.
- solution: `ionic_strength`'s i32 test is inline, and its finiteness scan runs only in the rare non-finite case.
  `solution_conductivity` goes straight to `_km_cr_mono` when the fast path fails.
- mathx:
  - `_km_cr_log10` inlines the ln core;
  - `_km_cr_expm1` inlines its fast path and exp's core, as `_km_cr_exp` does;
  - `_km_cr_acos` inlines `_km_acos_core` and one copy of `_km_atan_core` (A/B: 147 → 140 ns);
  - `_km_pow_pos_fast` inlines the ln and exp cores (A/B: 136 → 131 ns);
  - `_km_mul3_fast` and `_km_mul4_fast` inline their helpers (`rigid_rotor_energy`, A/B: 103 → 98 ns).

| Bench | Rust 1.2.0 | Cyrius 2.0.0 | Cyrius 2.0.1 | Change |
|---|---|---|---|---|
| `element/lookup_by_symbol` | 44.73 ns | 71 ns | 57 ns | −20 % |
| `molecule/molecular_weight_glucose` | 3.63 ns | 29 ns | 24 ns | −17 % |
| `thermochem/reaction_enthalpy_ch4` | 38.58 ns | 303 ns | 228 ns | −25 % |
| `thermochem/reaction_gibbs_energy` | — | 298 ns | 221 ns | −26 % |
| `thermochem/reaction_entropy` | — | 286 ns | 201 ns | −30 % |
| `thermochem/reaction_enthalpy_edge` | — | 592 ns | 492 ns | −17 % |
| `thermochem/enthalpy_change_cp` | 24.00 ns | 184 ns | 146 ns | −21 % |
| `thermochem/adiabatic_flame_ch4_air` | 301.48 ns | 3.55 µs | 2.04 µs | −42 % |
| `thermochem/adiabatic_flame_tiny` | — | 1.47 µs | 1.15 µs | −22 % |
| `stoichiometry/empirical_formula_glucose` | — | 612 ns | 487 ns | −20 % |
| `stoichiometry/empirical_formula_edge` | — | 7.81 µs | 3.90 µs | −50 % |
| `dynamics/first_order_100steps` | 4.41 µs | 10.32 µs | 9.35 µs | −9 % |
| `dynamics/stochastic_100steps` | — | 22.32 µs | 21.60 µs | −3 % |
| `fitting/arrhenius_6pts` | 233.40 ns | 987 ns | 944 ns | −4 % |
| `fitting/polynomial_deg2_20pts` | 500.05 ns | 2.28 µs | 2.08 µs | −9 % |
| `spectroscopy/fft_spectrum_256` | — | 30.37 µs | 25.21 µs | −17 % |
| `spectroscopy/fft_spectrum_1024` | — | 140.80 µs | 114.02 µs | −19 % |
| `spectroscopy/fft_spectrum_4096` | — | 656.29 µs | 528.29 µs | −20 % |
| `spectroscopy/fft_spectrum_256_huge` | — | 319.07 µs | 246.27 µs | −23 % |
| `spectroscopy/fft_spectrum_256_subnormal` | — | 341.28 µs | 267.35 µs | −22 % |
| `spectroscopy/fft_spectrum_256_squares` | — | 45.02 µs | 39.19 µs | −13 % |
| `spectroscopy/find_peaks_1000` | — | 20.42 µs | 11.63 µs | −43 % |
| `spectroscopy/find_peaks_1000_subnormal` | — | 38.59 µs | 29.30 µs | −24 % |
| `spectroscopy/isotope_pattern_glucose` | 2.80 µs | 15.91 µs | 13.32 µs | −16 % |
| `spectroscopy/isotope_pattern_cystine` | — | 18.68 µs | 15.49 µs | −17 % |
| `spectroscopy/isotope_pattern_s8` | — | 18.76 µs | 14.79 µs | −21 % |
| `nuclear/decay_chain_two_step` | — | 261 ns | 238 ns | −9 % |
| `potential/dihedral_angle` | 20.35 ns | 135 ns | 109 ns | −19 % |
| `potential/bond_angle` | — | 217 ns | 187 ns | −14 % |
| `potential/bond_angle_linear` | — | 140 ns | 96 ns | −31 % |
| `electrochemistry/nernst_potential` | 6.42 ns | 92 ns | 71 ns | −23 % |
| `electrochemistry/nernst_potential_edge` | — | 214 ns | 187 ns | −13 % |
| `solution/ionic_strength` | — | 32 ns | 26 ns | −19 % |
| `solution/solution_conductivity` | — | 74 ns | 71 ns | −4 % |
| `mathx/cr_log10` | — | 80 ns | 67 ns | −16 % |
| `mathx/cr_expm1` | — | 51 ns | 46 ns | −10 % |
| `mathx/cr_expm1_ux` | — | 64 ns | 56 ns | −12 % |
| `mathx/cr_acos` | — | 139 ns | 140 ns | +1 % |
| `mathx/cr_pow` | — | 138 ns | 138 ns | +0 % |
| `thermo/heat_transfer` | 0.95 ns | 68 ns | 62 ns | −9 % |
| `phase/boiling_point_elevation` | — | 69 ns | 61 ns | −12 % |
| `phase/freezing_point_depression` | — | 69 ns | 62 ns | −10 % |
| `spectroscopy/absorbance` | 0.98 ns | 67 ns | 63 ns | −6 % |
| `spectroscopy/diatomic_moment_of_inertia` | — | 80 ns | 68 ns | −15 % |
| `spectroscopy/rigid_rotor_energy` | — | 112 ns | 108 ns | −4 % |
| `thermo/heat_transfer_tie` | — | 377 ns | 311 ns | −18 % |

## [2.0.0] - 2026-11-13

**The Cyrius port.** kimiya 2.0.0 is a Cyrius library (toolchain 6.7.18, hisab 3.4.1), ported from the Rust crate
1.2.0, which stays frozen in `rust-old/` as the parity oracle (ADR 0001; plan: `docs/development/kimiya-port-plan.md`,
milestones M0–M7). All 22 in-scope modules are ported. Their golden corpus is 342,485 rows ok and 46,634 accepted, with
0 failed, on x86_64 ELF, x86_64 PE and aarch64. The accepted rows are ADR 0002's API changes and ADR 0006's correctly
rounded values. Of Rust's 688 tests, 652 are ported, 1 is n/a (serde) and 35 are deferred with their modules. Every one
of `rust-old/src`'s 519 public items is ported, renamed by ADR 0002, n/a or deferred
(`docs/development/rust-parity.md`). The bridges, the soorat integration and the `logging` feature wait on a Cyrius
consumer (port plan §2.2).

### Breaking (relative to Rust 1.2.0)
- The Cyrius API is not source-compatible with the Rust crate. It keeps Rust's names and semantics in Cyrius's form
  (port plan §6.3–§6.4: handles for table lookups, value-form `Result` / `Option`, caller buffers for `Vec` results,
  function pointers with a context for closures), and it changes Rust's API or results on purpose where ADR 0002 says
  so. *Migration:* [ADR 0002](docs/adr/0002-port-api-divergences.md) lists every change with how the parity harness
  compares it; the "Breaking (relative to Rust 1.2.0, ADR 0002 #N)" sections below give each one's migration note.
  Rust 1.x stays on crates.io for the Rust dependents, with 1.2.0's documented behaviour (there is no Rust 2.0)

### Added
- **M7 release** (port plan §4):
  - the Rust-test ledger `docs/development/rust-parity.md`, written and gated (`--check`) by
    `scripts/rust_test_ledger.py`: each of the 688 Rust tests (662 unit, 24 integration, 2 doc) is ported (652, 25 of
    them adapted to an ADR 0002 change or the flat API), n/a (1: molecule's serde round-trip) or deferred (35: bridge 22,
    soorat 12, logging 1); every port's `test_group` is the Rust name and its file's `main` calls it
  - `tests/integration.tcyr`: the 16 integration tests not yet ported and the 2 doc tests (the README's Quick Start,
    `adiabatic_flame_temperature`'s example), against the bundle `dist/kimiya.cyr`; three adapted to ADR 0002 #1's J/mol
  - `scripts/surface_parity.py` covers the whole surface: all 519 public items of `rust-old/src` (429 ported, 2 renamed
    by ADR 0002, 39 n/a, 49 deferred), lib.rs's module declarations and re-exports, every `pub` line accounted for, and
    the 124 derived trait impls (`--markdown` writes the ledger page's table)
  - `docs/benchmarks-rust-v-cyrius.md`: every Rust benchmark with its 2.0.0 Cyrius number and the 160 Cyrius-only ones
  - CI: the ledger's `--check` beside surface parity
  - `CONTRIBUTING.md` and `SECURITY.md` describe the Cyrius library (the cleanliness check, 2.x supported, the Rust
    1.x crate frozen at 1.2.0); the stray probe files `atoms.txt` and `probe.cyr` removed from the root;
    `tests/audit2.bcyr`'s subnormal set by its bits (its one f64 warning in the bench build)
- **M0 initialize** — `cyrius port` scaffold on Cyrius 6.7.18; `cyrius.cyml` with hisab 3.4.1; golden-value
  generator `scripts/oracle/` (constants 45 + math 18,647 rows, glibc 2.44-pinned, `--check`, independent
  `validate.py`); generated `src/constants.cyr` (14 `KM_` physical constants); parity harness
  `tests/parity.tcyr` — 13,340 rows ok, 0 failed on x86_64 ELF, x86_64 PE and aarch64 (`docs/development/port-parity.md`)
- ADRs 0001 (port with Rust oracle), 0002 (API changes in the Cyrius 2.0 API), 0003 (transcribe the hisab 1.4.0
  routines hisab 3.4.1 lacks), 0004 (accept ganita `pow(-1, DBL_MAX)`)
- **M1 foundation** (port plan §4, §6.4):
  - **error** — `enum KmError` (`KM_ERR_*`, Rust variant order, append-only), `km_error_name`, `km_error_message`;
    fallible functions return the value-form `Result` with the code as `Err` payload
  - **mathx** — private `_km_min` / `_km_max` / `_km_signum` / `_km_round_half_away` / `_km_powi` /
    `_km_total_cmp` / `_km_is_finite` with Rust's semantics, and `_km_log10` / `_km_expm1` transcribed from
    Rust libm 0.2.16 (fdlibm): the `math` corpus has 0 pending rows; exp_m1 is bit-identical to glibc on all 810
    rows, log10 within 1 ulp (48 rows, where glibc is the one off correct rounding); class-C K = 2
  - **validate** — `_km_is_positive`, `_km_is_non_negative` (NaN/±∞-rejecting guards)
  - **element** — generated `src/element_data.cyr` (`enum KmElementCategory`, 118 element rows with
    NUL-separated symbol/name blobs, 32 property rows with `i8` oxidation states); `element_lookup_by_number`
    (O(1)), `element_lookup_by_symbol`, `element_lookup_properties` returning handles (0 for none) and
    `element_*` / `element_properties_*` accessors that never read out of bounds
  - **molecule** — `struct KmMolecule` (`new`, `water`, `carbon_dioxide`, `methane`, `glucose`,
    `molecular_weight`, `formula`, `total_atoms`), bit-identical to Rust
  - golden modules `element` (975 rows) and `molecule` (2,472 rows): 22,107 rows ok, 0 failed on x86_64 ELF,
    x86_64 PE and aarch64; `probe_tables.py` compiles all 1,055 element_data cells; `validate.py` recomputes
    both modules from the Rust source literals
  - 38 Rust tests ported (`tests/{error,validate,element,molecule}.tcyr`); `tests/kimiya.tcyr` now smoke-tests
    the bundle
  - `[lib]` and `dist/kimiya.cyr` (`cyrius distlib`); `scripts/surface_parity.py` (55/55); CI runs
    `distlib --check` and surface parity
  - benchmarks `element/lookup_by_symbol`, `element/lookup_by_number`, `molecule/molecular_weight_glucose`
    against Rust ([`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md))
- **M2a reaction, thermo, phase, gas** (port plan §4, §6.4; M2b, thermochem, follows):
  - **reaction** — `gibbs_free_energy`, `equilibrium_constant`, `is_spontaneous`
  - **thermo** — `heat_transfer`, `hess_law(xs, n)`, `enthalpy_from_formation(p, np, r, nr)` (a Rust slice is a
    pointer and a count; the empty sum is Rust's −0.0), `specific_heat_from_calorimetry` (hisab's `EPSILON_F64`
    threshold); `KM_WATER_SPECIFIC_HEAT` generated into `src/constants.cyr`
  - **phase** — generated `src/phase_data.cyr` (`struct KmAntoineCoeffs`, the 8 `ANTOINE_DATA` sets as columns
    with a name blob, `KM_ANTOINE_COUNT`, the mmHg → Pa factor, `KM_KB_*` / `KM_KF_*`); Raoult, Henry,
    colligative properties, osmotic pressure, Clausius–Clapeyron, Antoine (ADR 0002 #4, below),
    `lookup_antoine` (a handle, 0 for none) with `antoine_data_name` / `antoine_data_coeffs`, `phase_rule`
  - **gas** — generated `src/gas_data.cyr` (`struct KmPengRobinsonParams`, `KM_STP_*`, the `KM_PR_*` presets as
    struct values); ideal gas, van der Waals, Dalton, density, Peng–Robinson pressure / Z / roots (a caller's
    three-f64 buffer, `Ok(n)`) / fugacity coefficient / fugacity with Rust's analytic cubic transcribed,
    Redlich–Kwong, `joule_thomson_ideal`
  - **mathx** — correctly rounded `_km_acos`, `_km_cos`, `_km_ln`, `_km_pow_2_5` (double-double, the same bits on
    every target) for the Peng–Robinson and Redlich–Kwong functions, which are ill-conditioned near a double
    root, at extreme states and near P = 0
    ([ADR 0005](docs/adr/0005-correctly-rounded-primitives-for-ill-conditioned-functions.md))
  - golden modules `reaction` (1,188 rows), `thermo` (1,107), `phase` (3,276), `gas` (19,717), and 7,500 more
    `math` rows: 54,466 ok, 0 failed, 0 measured on x86_64 ELF, x86_64 PE and aarch64; every class-B row and every
    PR root bit-identical; class-C bounds K = 2–10 set by the port-parity rule
    ([`port-parity.md`](docs/development/port-parity.md) §5); `validate.py` recomputes all 25,288 rows and
    measures each class-C function's sensitivity (`--sensitivity`); `probe_tables.py` compiles the 65 new table cells
    and 11 new constants
  - 90 Rust tests ported (`tests/{reaction,thermo,phase,gas}.tcyr`, none skipped) plus 6 Cyrius-only groups and a
    `validate` group for the faster guards' edges (DBL_MAX, the smallest subnormal, −NaN, ±∞, −0);
    surface parity 113/113; benches `reaction/gibbs_free_energy`, `thermo/heat_transfer`,
    `gas/ideal_gas_pressure`, `gas/peng_robinson_z_n2_1mpa`
  - ADR 0005 (correctly rounded primitives for the ill-conditioned gas functions)
- **M2b thermochem** (port plan §4, §6.4; M2 is complete):
  - generated `src/thermochem_data.cyr`: `struct KmShomate` (with its Rust field docs), the 76 `THERMOCHEM_DATA` rows
    (formula / name blobs, ΔH_f°, ΔG_f°, S°) and the 34 `SHOMATE_DATA` segments as columns, `KM_THERMOCHEM_COUNT`,
    `KM_SHOMATE_COUNT`, with the literals' sources and caveats (the unsourced C₂H₆ / C₃H₈ / HF fits) as comments
  - **thermochem** — `lookup_thermochem` (a handle, 0 for none) with `thermochem_data_formula` / `_name` /
    `_delta_hf_kj` / `_delta_gf_kj` / `_s_standard_j`; `reaction_enthalpy`, `reaction_gibbs_energy` (J/mol, below),
    `reaction_entropy` over species lists (a pointer to n 16-byte (formula, coefficient) records and n);
    `vant_hoff_k`; `KmShomate.cp` / `.enthalpy_antiderivative`; `lookup_shomate`, `lookup_shomate_at` with
    `shomate_data_formula` / `shomate_data_coeffs`; `enthalpy_change_cp` (analytic, piecewise across segments);
    `adiabatic_flame_temperature` with hisab 1.4.0's `newton_raphson` and `bisection` transcribed as private helpers
    ([ADR 0003, amendment M2b](docs/adr/0003-transcribe-missing-hisab-routines.md)): the whole path is + − × ÷, so the
    flame temperature is bit-identical to Rust's, on the Newton and on the bisection fallback path
  - golden module `thermochem` (9,487 rows): 63,953 ok, 0 failed, 0 measured on x86_64 ELF, x86_64 PE and aarch64;
    every class-B row 0 ulp (1,836 flame temperatures: 1,682 from Newton, 154 from the bisection fallback);
    `vant_hoff_k` class C, K = 6 by the port-parity rule ([`port-parity.md`](docs/development/port-parity.md) §6);
    each flame row's note names its solver path, and `validate.py` recomputes all 9,487 rows and the paths, and checks
    that Rust's typed flame constants are pinned by the corpus; `probe_tables.py` compiles the 654 new table cells
  - 62 Rust tests ported (`tests/thermochem.tcyr`, none skipped; 4 assert J/mol for ADR 0002 #1) plus 6 Cyrius-only
    groups (species lists, handles, J/mol, the solver's tolerances and derivative threshold, > 32 products, the struct
    methods): 400 assertions; surface parity 140/140; benches `thermochem/reaction_enthalpy_ch4`,
    `thermochem/vant_hoff_k`, `thermochem/enthalpy_change_cp`, `thermochem/shomate_cp`,
    `thermochem/adiabatic_flame_ch4_air`
  - generator: `Module::row_noted` (a bracketed note suffix), `cyr::string_column` splits blob lines past 120 columns
- **M3a solution, electrochemistry** (port plan §4, §6.4; M3b, kinetics, follows):
  - **solution** — `molarity`, `molality`, `dilution`, `ph_from_h_concentration`, `poh_from_oh_concentration`,
    `ph_from_poh`, `henderson_hasselbalch`, `h_concentration_from_ph`, `weak_acid_ph` / `weak_base_ph` (Rust's
    charge-balance Newton transcribed as the private `_km_charge_balance_root` / `_km_checked_root`: the root is
    bit-identical), `ionic_strength(ions, n)` over 16-byte (concentration, charge) records, `debye_huckel_limiting`,
    `debye_huckel_extended`, `davies_activity`, `kohlrausch_molar_conductivity`, `limiting_molar_conductivity`,
    `solution_conductivity`; its ten constants (`KM_KW_WATER_25C`, `KM_DH_A_WATER_25C`, `KM_DH_B_WATER_25C`, the seven
    `KM_LAMBDA_*`) generated into `src/constants.cyr`
  - **electrochemistry** — generated `src/electrochemistry_data.cyr` (the 49 `STANDARD_POTENTIALS` rows: a UTF-8
    equation blob, the couple blob, electrons, E°; `KM_STANDARD_POTENTIALS_COUNT`; the literal's CRC Handbook notes);
    `lookup_half_reaction` (a handle, 0 for none) with `half_reaction_equation` / `_couple` / `_electrons` /
    `_standard_potential`; `nernst_potential`, `nernst_potential_25c`, `standard_cell_potential`,
    `cell_potential_from_couples`, `is_spontaneous_cell`, `cell_gibbs_energy`, `faraday_mass_deposited`,
    `faraday_charge_required`, `faraday_mass_ratio`, `charge_to_moles_electrons`, `moles_electrons_to_charge`
  - Rust's integer types kept as words and never wrapped: an i32 charge outside its range is
    `Err(KM_ERR_INVALID_INPUT)` (NaN from `ionic_strength`), a u8 electron count outside 0 ..= 255 is
    `Err(KM_ERR_INVALID_INPUT)` where Rust checks `n == 0`
  - **mathx** — `_km_log10_ln`, glibc 2.44's log10 composition (Sun's e_log10.c) over the correctly rounded `_km_ln`,
    for `henderson_hasselbalch`; the Nernst functions use `_km_ln`. Both cancel by construction near pH 0 / E ≈ 0
    (S ≈ 10¹⁸ with the stdlib primitives); a correctly rounded log10 was tried and failed the near-cancellation rows,
    because glibc's log10 is not correctly rounded on 0.3–29 % of arguments
    ([ADR 0005, amendment M3a](docs/adr/0005-correctly-rounded-primitives-for-ill-conditioned-functions.md))
  - golden modules `solution` (28,566 rows) and `electrochemistry` (18,574), and the 10 solution `constants` rows:
    111,103 ok, 0 failed, 0 measured, 10 pending on x86_64 ELF, x86_64 PE and aarch64; every class-B row 0 ulp;
    class C K = 2–4 by the rule ([`port-parity.md`](docs/development/port-parity.md) §7), Henderson–Hasselbalch and
    Nernst 0 ulp on every row, the near-cancellation rows included; each weak acid / base row's note names the
    Newton's exit and update count (at most 7), and `validate.py` recomputes all 47,140 rows and the paths and checks
    that `CHARGE_BALANCE_REL_TOL` is pinned from both sides (the 50-iteration cap is unobservable);
    `probe_tables.py` compiles the 197 new table cells and 10 new constants
  - 84 Rust tests ported (`tests/solution.tcyr` 49, `tests/electrochemistry.tcyr` 35, none skipped) plus 8 Cyrius-only
    groups (ion records, out-of-i32 charges, the private Newton and its tolerance, the log10 composition, the
    constants, handles, out-of-u8 electron counts, the Nernst ln): 1,529 assertions; surface parity 187/187; benches
    `solution/ph_from_h`, `solution/weak_acid_ph`, `electrochemistry/lookup_half_reaction`,
    `electrochemistry/nernst_potential`
  - generator: `cyr::string_column_utf8` (a UTF-8 blob; the ASCII tables are unchanged)
- **M3b kinetics** (port plan §4, §6.4; M3 is complete):
  - **kinetics** — `arrhenius_rate`, `half_life_first_order`, `first_order_concentration`,
    `second_order_concentration`, `zero_order_concentration`, `zero_order_half_life`, `nth_order_half_life`,
    `steady_state_intermediate`, `steady_state_rate`, `pre_equilibrium_rate`, `michaelis_menten`, `lineweaver_burk`,
    `collision_theory_rate`, `eyring_rate`, `eyring_rate_from_activation`, in Rust's check and operation order, on
    `f64_exp`, `_km_expm1` and `ganita_f64_pow`; `KM_BOLTZMANN` and `KM_PLANCK` are the generated constants; Rust's
    `LN_2` and private 10⁻¹² first-order tolerance are typed inline and pinned by the corpus
  - `lineweaver_burk(v_max, km, out)` writes the (slope, intercept) tuple to the caller's two-f64 buffer and returns
    `Ok(2)`, as `peng_robinson_z_roots` does its roots; on `Err` the buffer is untouched
  - golden module `kinetics` (30,146 rows, T from the smallest subnormal to DBL_MAX, wide signed energies, orders
    around 1 and out to ±300, and on purpose the rows where a subnormal exp or pow result is scaled into the normal
    range): 119,044 ok, 0 failed, 22,205 measured, 10 pending on x86_64 ELF, x86_64 PE and aarch64; every class-B row
    0 ulp; `validate.py` recomputes all 30,146 rows and checks that LN_2 and the tolerance are pinned
  - **the six class-C kinetics functions are measured, not gated**: a subnormal exp or pow result times a prefactor
    puts the rule's S at 10¹⁴ … 4.4 × 10¹⁸ (2–4 without those rows). Under the user's decision of 2026-11-09 —
    functional parity with correctly rounded primitives, never glibc's rounding errors — no glibc-matching primitive
    is built; their bound waits for ADR 0006 ([`port-parity.md`](docs/development/port-parity.md) §8). On two exact
    half-way rows ((3·2⁻²¹⁵)⁵, (9·2⁻⁴³⁰)^2.5) ganita's pow is correctly rounded and glibc's is not
  - `validate.py --sensitivity` lists, for each function past K = 64, the rows that put it there and S without them
  - 47 Rust tests ported (`tests/kinetics.tcyr`; `lineweaver_burk_basic` and `non_finite_inputs_are_rejected` adapted
    for the buffer) plus 4 Cyrius-only groups (the buffer, the inline constants and the tolerance's boundary orders,
    the check order, signed zeros and NaN): 138 assertions; surface parity 203/203; benches `kinetics/arrhenius_rate`,
    `kinetics/eyring_rate`, `kinetics/michaelis_menten`
- **M3c correct rounding** ([ADR 0006](docs/adr/0006-correct-rounding-is-the-accuracy-target.md): correct rounding is the
  accuracy target, parity with the oracle is functional; M3 is complete):
  - **the corrected oracle** (M3c-1) — `scripts/oracle/mpfr_ref.py`, correctly rounded binary64 results from the system
    MPFR 4.2 through ctypes (exact ties, subnormal results and overflow right; development only). `validate.py` decides
    "is glibc correctly rounded here" with it everywhere (the 60-digit `decimal` references, which misjudge exact ties,
    are gone), and evaluates every golden row whose replica calls a transcendental primitive a second time with each
    primitive result correctly rounded, through the same transcribed code (a changed branch or an Ok ↔ Err flip is
    followed; none occurs on the corpus). The 665 rows that change are written by `validate.py . --write-corrections` to
    `tests/golden/<module>.cr.tsv` (`gas` 11, `kinetics` 323, `math` 96, `phase` 4, `reaction` 2, `solution` 229), keyed by
    the golden row's exact fn and inputs, each note naming the oracle's outputs and the calls glibc misrounded; a plain
    `validate.py .` fails on a missing, stale or extra correction and rechecks every listed call against MPFR
  - **correctly rounded primitives** (M3c-2) — `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`
    and `_km_cr_pow` in `src/mathx.cyr`: Ziv's two phases (a fast path with a proven error bound and a rounding test, a
    triple-double accurate path), the exact and midpoint cases explicit (10⁰ … 10²² and the tie 10²³, log10 of the powers
    of 10, pow's exact and half-way powers, subnormal ties included), + − × ÷ and integer bit operations only, so the same
    bits on every target; `_acc` / `_fast` entry points for the tests. Their tables and constants are generated into
    `src/mathx_data.cyr` by `scripts/oracle/mathx_tables.py` from MPFR (asserting the bit-count properties the error
    analysis uses); `scripts/oracle/cr_verify.py` checks ≥ 1,000,000 seeded arguments per primitive (≥ 1,300,000 pairs for
    pow) against MPFR on x86_64 ELF, x86_64 PE and aarch64: 0 mismatches, identical bits and flags on all three;
    `tests/mathx.tcyr` (11 groups, 518 assertions) pins the special, exact, tie, threshold, subnormal and hard cases
  - **every module on them** (M3c-3), in Rust's operation order: exp → `_km_cr_exp` (reaction, phase, thermochem, gas,
    kinetics), exp_m1 → `_km_cr_expm1` and powf → `_km_cr_pow` (`nth_order_half_life`), `10.0.powf(x)` → `_km_cr_exp10`
    (Antoine, `h_concentration_from_ph`, the activity models), log10 → `_km_cr_log10` (solution, Henderson–Hasselbalch
    included), ln → `_km_cr_ln` (Nernst, fugacity), `tc.powf(2.5)` → `_km_cr_pow(tc, 2.5)` (Redlich–Kwong); ADR 0005's
    `_km_acos` / `_km_cos` stay. `phase`, `solution` and `kinetics` no longer include `lib/ganita.cyr` (`gas` keeps it for
    cbrt); `dist/kimiya.deps` is unchanged
  - **the harness** — every class-C function at **bound 0** against the corrected oracle (27 kimiya functions, the six
    kinetics ones gated for the first time, and the `math` primitives through the port's correctly rounded ones); a row
    that matches only the oracle's misrounded value now **fails** (`glibc-matching`, M3c-1's transitional *measured*
    category); ADR 0004's `pow(-1, DBL_MAX)` row is a plain ok under `_km_cr_pow`. `f64::sin`, `f64::atan2` and the 247
    `f64::cos` rows past |x| = 8 stay on the stdlib at 2 ulp until M4 brings correctly rounded sin and atan2. Shadow lines
    show the retired stdlib primitives (`f64_exp`, `f64_ln`, `ganita_f64_acos`, `ganita_f64_pow`) against the corrected
    oracle
  - **hard rows** (ADR 0006 item 6): 2,486 new golden rows where a subnormal exp result is scaled back by a prefactor —
    `vant_hoff_k`, `clausius_clapeyron`, `equilibrium_constant`'s subnormal end (370 each) and the Peng–Robinson fugacity
    at a few kelvin (344 states: water near 10 K, CO₂ near 4 K, methane near 1.7 K), hand-picked on
    `grids::EXP_SUBNORMAL_X` plus seeded sweeps. glibc's exp is correctly rounded on all of them (no new correction);
    with the stdlib's `f64_exp` instead of `_km_cr_exp` they would be up to 27 ulp off on aarch64 (its polyfill rounds
    twice); with it they are 0 ulp on every target
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `143071 ok, 0 failed, 0 measured, 10 pending, 1105 accepted (ADR 6:
    665), 0 glibc-matching (failed)` — every class-C row 0 ulp from its reference
    ([`port-parity.md`](docs/development/port-parity.md) §9–§10); 15 test files, 0 warnings; surface parity 203/203
  - `validate.py --sensitivity` is a diagnostic of conditioning: every primitive result is nudged (S for
    `peng_robinson_z_roots` 2²⁶, the Nernst functions and `henderson_hasselbalch` ~10¹⁸, the subnormal-prefactor rows
    4.4 × 10¹⁸), the retired rule's K is printed for comparison, and no bound is derived from it
- **M4b organic and inorganic** (port plan §4, §6.4 "Decided in M4b"; M4 in progress):
  - **organic** — generated `src/organic_data.cyr`: `enum KmBondType` (45), `enum KmGeometry` (12), `enum
    KmFunctionalGroup` (16), explicit append-only values in Rust variant order under the Rust docs and section headings
    (the generator proves its variant lists against the source declarations, `validate.py` parses them again), and their
    property columns through the oracle's API — bond energies (with the `// Carbon bonds` headings of Rust's match),
    ideal angles (with the arms' notes), the UTF-8 pattern blob, the polar and donor flags, and ADR 0002 #6's acceptor
    column. `src/organic.cyr`: `bond_type_energy_kj(b)`, `geometry_ideal_bond_angle(g)`, `functional_group_pattern(g)`,
    `functional_group_is_polar(g)`, `functional_group_is_hydrogen_bond_donor(g)`, the new
    `functional_group_is_hydrogen_bond_acceptor(g)` (NaN, `""` or 0 for a word that is no value of the enum);
    `predict_geometry(bonding, lone)` → the value-form `Option` of a `KmGeometry` (None outside u8, never wrapped);
    `enthalpy_from_bonds(broken, nb, formed, nf)` over 16-byte (value, count) records, Rust's two `.sum()`s from −0.0
    (NaN for a bad list, an invalid value or a count outside u32)
  - **the hydrogen-bond acceptor classification** (ADR 0002 #6; no oracle data, so reference tests only): IUPAC 2011's
    criterion (E. Arunan et al., *Pure Appl. Chem.* 83, 1637–1641), an electron-rich region such as a lone pair or a
    π-bonded pair, the weak acceptors included (the user's decision): every group but the alkane accepts — the O and N
    lone pairs (McMurry, *Organic Chemistry*, OpenStax 2023, §2.12), X–H···π to alkenes, alkynes and aromatic rings
    (T. Steiner, *Angew. Chem. Int. Ed.* 41, 48 (2002)), the thiol's divalent sulfur (F. H. Allen et al., *Acta Cryst.*
    B53, 696 (1997)) and C–X halogens (L. Brammer et al., *Cryst. Growth Des.* 1, 277 (2001); organic fluorine the
    weakest, J. D. Dunitz & R. Taylor, *Chem. Eur. J.* 3, 89 (1997)). No strength is implied. The sources and each row's
    reason are in the generated comment; `validate.py` and `probe_tables.py` recompute the column from the rule (a
    lone-pair atom or a π bond in the group's pattern)
  - **inorganic** — generated `src/inorganic_data.cyr`: `enum KmLigand` (11) and `KM_LIGAND_FIELD_STRENGTH`; Rust's
    private octahedral occupancy and unpaired-electron tables, read through `cfse_octahedral` / `unpaired_electrons_oct`
    (as M2a's `PA_PER_MMHG`) and cross-checked against their literals; the ten `KM_MADELUNG_*` / `KM_BORN_EXP_*` with their
    section comments (**the last 10 pending `constants` golden rows: 0 pending now**); the 30 Shannon (1976) `IONIC_RADII`
    rows (symbol blob, `i8` charge, coordination, radius) and `KM_IONIC_RADII_COUNT`. `src/inorganic.cyr`:
    `cfse_octahedral` (exact, ADR 0002 #12), `delta_tet_from_oct`, `ligand_field_strength`, `ligand_is_strong_field`, `spin_only_moment`,
    `unpaired_electrons_oct`, `born_lande_lattice_energy`, `born_haber_lattice_energy`, `lookup_ionic_radius` (handle =
    row + 1) and the `ionic_radius_*` accessors. Rust's u8 / u32 / i8 are words never wrapped (Err, NaN or None outside
    the type), and a Rust `bool` is 0 / 1, any other word `Err(KM_ERR_INVALID_INPUT)` (port plan §6.4)
  - golden modules `organic` (9,104 rows: class A 7,596, B 1,508) and `inorganic` (15,107: A 9,397, B 5,710), every enum
    value × property, `predict_geometry` over 7,491 (bonding, lone) pairs, every d count × spin, every u8 moment, every
    table symbol × i8 charge, Born–Landé / Born–Haber over wide grids and every Err path; `validate.py` recomputes every
    row from the Rust source (and perturbs each typed formula constant: each changes some row, the strong-field 8 from
    above only); `probe_tables.py` compiles all 197 organic_data and 220 inorganic_data cells
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `167279 ok, 0 failed, 0 measured, 0 pending, 1118 accepted (ADR 6:
    665), 0 glibc-matching (failed)` — the 13 new accepted rows are ADR 0002 #12's, printed by name — every M4b row 0 ulp
    from its reference ([`port-parity.md`](docs/development/port-parity.md) §11)
  - 51 Rust tests ported one for one (`tests/organic.tcyr` 21, `tests/inorganic.tcyr` 30; 5 organic tests use ADR 0002
    #6's names) plus the two organic integration tests, and 10 Cyrius-only groups (the acceptor reference tests, invalid
    enum values and handles, out-of-range u8 / u32 / i8 / bool words, the bond records, the exact CFSE): 464
    assertions; 17 test files, 0 warnings; surface parity 323/323 (`surface_parity.py` maps enum methods to
    `<enum>_<method>` and now uses the generator's snake rule)
- **M4c nuclear** (port plan §4, §6.4 "Decided in M4c"; M4 in progress):
  - generated `src/nuclear_data.cyr`: the 18 `ISOTOPES` rows (NUBASE2020 / ENSDF half-lives in s) as columns — the label
    blob, `u8` Z, `u16` A, the half-life, a UTF-8 decay-mode blob ("β⁻/β⁺/EC") — and `KM_ISOTOPES_COUNT`, with `Isotope`'s
    docs and the literal's three NUBASE2020 notes; the four nuclear constants stay `src/constants.cyr`'s
  - **nuclear** — `decay_constant`, `atoms_remaining`, `activity`, `atoms_after_half_lives`, `time_for_fraction`,
    `mass_defect`, `binding_energy_mev`, `binding_energy_per_nucleon`, `semi_empirical_binding_energy` (a `Result`, below),
    `q_value(r, nr, p, np)` (two slices, Rust's `.sum()` from −0.0 through thermo's `_km_sum`), `decay_chain_two_step(n0,
    λ₁, λ₂, t, out)` (the (N_A, N_B, N_C) tuple in the caller's 24-byte buffer, `Ok(3)`; Rust's cancellation-free Bateman
    form and its equal-rate limit at |Δλ| ≤ 10⁻¹²·max λ kept exactly), `is_secular_equilibrium`, `lookup_isotope` (a handle,
    0 for none) with `isotope_symbol` / `_z` / `_a` / `_half_life_s` / `_decay_mode`. Rust's check and operation order;
    exp, exp_m1, ln and `2.0.powf(−n)` are the correctly rounded `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_ln` and
    `_km_cr_pow(2, −n)`, cbrt ganita's (correctly rounded); Rust's u32 counts are words, never wrapped (`Err`, or NaN for
    `mass_defect`, outside u32); `LN_2`, the private 10⁻¹² tolerance and SEMF coefficients and the asymmetry's 4 are typed
    inline and pinned by the corpus. `atoms_after_half_lives` and `is_secular_equilibrium` follow ADR 0002 #13 / #14
    (below, Breaking)
  - golden module `nuclear` (54,057 rows: class A 1,898, B 35,370, C 16,789): the SEMF over every Z ≤ A ≤ 3Z with A ≤ 300
    (port plan §4), Z = 0 and the u32 extremes; decay laws and chains over wide λ, t and N₀ grids with every Err path,
    ±0, NaN, ±∞, DBL_MAX and subnormals, λ₁ = λ₂ exactly, λ₁ ≈ λ₂ either side of the threshold and **pairs on it to the ulp
    of 10⁻¹²**, and the subnormal-exp-and-pow-times-prefactor rows (ADR 0006 item 6); every isotope and isotope pair.
    23 rows have an ADR 0006 correction (`tests/golden/nuclear.cr.tsv`, glibc's expm1). `validate.py` recomputes every
    row (glibc's exp2 for the release oracle's `2.0.powf`, a new `_prim` kind), checks the coverage and ADR 0002 #3's
    notes, perturbs each typed constant, checks `src/nuclear_data.cyr` against the literal, and proves every ADR 0002
    #13 / #14 / #15 note with exact rational arithmetic (`Fraction`, MPFR), with seeded exact sweeps beyond the corpus;
    `probe_tables.py` compiles the 91 new cells
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `219266 ok, 0 failed, 0 measured, 0 pending, 3243 accepted (ADR 6:
    688), 0 glibc-matching (failed)` — the new accepted rows are 23 under ADR 6 and, under ADR 2 and printed by name,
    #3's 23, #13's 1,812, #14's 208 and #15's 59 — every row 0 ulp from its reference (the oracle's, the corrected
    oracle's or the note's exact value), every class-C function at bound 0
    ([`port-parity.md`](docs/development/port-parity.md) §12)
  - **ADR 0002's inexact-formula policy, reviewed and decided** (the user, 2026-11-09: "fix exactly"):
    `atoms_after_half_lives` (#13), `is_secular_equilibrium` (#14) and M4b's `delta_tet_from_oct` overflow (#15) now return
    the exact values (Breaking, below); `decay_constant`'s LN_2/t½ is kept (ln 2/t½ is irrational). The harness reads a
    mapped row's exact value from its note (`parity_core.cyr`'s `_par_note_*`)
  - 36 Rust tests ported one for one (`tests/nuclear.tcyr`; `semi_empirical_invalid_is_nan` asserts ADR 0002 #3's `Err`;
    #13 / #14 change none) plus 8 Cyrius-only groups (the SEMF `Result`, out-of-u32 words and the checked add, the chain
    buffer and threshold, the `q_value` slices, the isotope handles and UTF-8 modes, the typed constants and Rust's edges,
    #13's exact N₀/2ⁿ, #14's exact comparison): 236 assertions; `tests/inorganic.tcyr` gains #15's group; the bundle smoke
    test gains 6; 18 test files, 0 warnings; surface parity 347/347; benches
    `nuclear/binding_energy_per_nucleon_fe56`, `nuclear/semi_empirical_fe56` (against Rust) and
    `nuclear/decay_chain_two_step`
  - harness: a mapped `Err` ↔ NaN row's `accepted` line prints no bit pair (`parity_core.cyr`)
- **M4a correctly rounded trigonometry** ([ADR 0006](docs/adr/0006-correct-rounding-is-the-accuracy-target.md) item 1;
  M4 in progress):
  - **mathx** — `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`, `_km_cr_atan2` and `_km_cr_acos` in `src/mathx.cyr`, each
    two-phase like M3c-2's primitives (a fast path with a proven error bound and a rounding test, a triple-double
    accurate path, `_fast` / `_acc` entry points), + − × ÷, √ and integer bit operations only, the same bits on every
    target. sin and cos over the whole double range: Cody–Waite reduction by a five-part π/2 below 2²⁰, Payne–Hanek past it
    (1,344 bits of 2/π in exact 32-bit limbs, the fraction centred and summed into a triple-double), then sin/cos(i/256)
    from a four-word table and a short polynomial in the remainder (fast path within 2⁻⁷⁰ relative), or a triple-double
    Taylor series (below 2⁻¹⁴⁸ by the operation count; 2⁻¹³⁸ at worst, the reduction); `_km_cr_sincos` gives both. atan2
    reduces with atan(j/256) (triple-double table) and δ = (n − c·d)/(d + c·n) built from exact partial products, never
    rounding y/x first; a degree-9 polynomial in δ on the fast path (2⁻⁶⁹·⁴⁴), a triple-double long division and series on
    the accurate one (2⁻¹⁴⁹); a quotient below 2⁻⁵⁹ is rounded directly (its exact subnormal midpoints too); C99 Annex F's
    special values. acos is atan2(√(1 − x²), x) with 1 − x² nearly exact (double-double, then an exact four-term sum) and
    √ by Newton steps. The tables and constants (sin/cos(i/256), atan(j/256), the bits of 2/π checked against Machin's
    formula in integers, π/2's Cody–Waite parts, the polynomials) are generated into `src/mathx_data.cyr` by
    `scripts/oracle/mathx_tables.py` from MPFR, which also proves the reduction bound (no double within 4.687 × 10⁻¹⁹ of a
    nonzero multiple of π/2, by continued fractions, attained at Lefèvre and Muller's worst case 6381956970095103·2⁷⁹⁷)
  - **`_km_td_norm`** — the accurate atan and acos paths divide, and take the square root of, triple-doubles formed after a
    cancellation (n − c·d; the residual w − s²). `_km_td_add`'s last renormalisation does not see its rounded
    third-order term, so after a cancellation its leading word can be far from the sum, while the long division's and
    Newton step's corrections (each q = r₀/b₀) assume the leading word is the sum to 2⁻⁵². `_km_td_norm` renormalises
    (two bottom-up TwoSum passes, exact) before each division and correction, which is what makes `_km_td_div`'s 2⁻¹⁵¹
    and `_km_td_sqrt`'s 2⁻¹⁵⁰ hold: without it the accurate path's quotient could be wrong far above its bound, and with
    it MPFR agrees on every argument tried
  - **verification** — `scripts/oracle/cr_verify.py` gains sin, cos, sincos, atan2 and acos: about 1.0–1.1 million seeded
    arguments each (every binade up to 2¹⁰²³, the Peng–Robinson range, the doubles nearest multiples of π/2 in every
    binade, the Cody–Waite / Payne–Hanek edge, the FFT angles −2π/n, tiny |x|; atan2 over every quadrant, exponent
    difference and the bond / dihedral geometries at Å and metre scale, subnormal and midpoint quotients; acos dense near
    ±1, 0 and the branch points; Lefèvre–Muller's and Hanrot et al.'s hard cases, every C99 special pair, glibc's atan2
    on the specials), plus 100,000 each forced through the accurate path: **0 mismatches against MPFR on x86_64 ELF,
    x86_64 PE and aarch64**, the same bits and flags on all three; the fast path falls back on 7.0 × 10⁻⁵ (sin),
    8.7 × 10⁻⁵ (cos), 1.9 × 10⁻⁵ (atan2) and 1.4 × 10⁻⁵ (acos) of the arguments. `tests/mathx.tcyr` gains the special
    values, hard cases, Payne–Hanek arguments, sincos, atan2's sign and zero matrix and acos (every expected value
    MPFR's) ([`port-parity.md`](docs/development/port-parity.md) §13)
- **M4d potential, the trig switch, ADR 0002 #13 completed** (port plan §4, §6.4 "Decided in M4d"; M4 in progress):
  - **potential** — `lennard_jones`, `lj_equilibrium_distance`, `morse`, `morse_width_parameter`, `coulomb`,
    `coulomb_force`, `distance_3d(p1, p2)`, `bond_angle(a, b, c)`, `dihedral_angle(a, b, c, d)` in Rust's check and operation
    order; a Rust `[f64; 3]` is a pointer to three contiguous f64 (NaN or `Err(KM_ERR_INVALID_INPUT)` for a null
    pointer); Rust's private `sub3` / `cross` / `dot3` are `_km_sub3` / `_km_cross` (multi-value returns) and `_km_dot3`;
    the private 10⁻¹² and 10⁻²⁰ tolerances are typed inline and pinned to the ulp by the corpus; exp and atan2 are the
    correctly rounded `_km_cr_exp` and `_km_cr_atan2`; `KM_COULOMB_CONSTANT` is the generated constant (no data file).
    `lj_equilibrium_distance`, `lennard_jones` and `morse_width_parameter` follow ADR 0002 #16 / #17 / #18 (Breaking,
    below)
  - **the trig switch** — `gas.cyr`'s Peng–Robinson cubic calls M4a's `_km_cr_acos` and `_km_cr_cos`: all 21,093 gas
    rows keep their bits. The `math` corpus's sin, cos (all 2,910 rows, the 247 past |x| = 8 included), atan2 and acos
    rows are gated at bound 0 on `_km_cr_sin` / `_km_cr_cos` / `_km_cr_atan2` / `_km_cr_acos`: no golden transcendental row
    is on the stdlib any more (`f64_sin`, `f64_cos` and `ganita_f64_atan2` are never-gated shadow lines)
  - **mathx** — `_km_cr_exp2_scaled(x, y)`, x·2^y rounded once (exp's fast core and rounding test at the final scale,
    exp's accurate core in triple-double as the fallback; subnormal results at the 2⁻¹⁰⁷⁴ grid), for ADR 0002 #13; and
    `_km_cr_sixth_root2_mul(x)`, x·2^(1/6) rounded once (a triple-double product and one `_km_cr_finish`), for #16, with
    2^(1/6) generated as a triple-double by `mathx_tables.py` from an exact integer root and its continued fraction
    proving that no rounding test is needed (every x·2^(1/6) is at least 2⁻⁵⁵·¹ ulp from a rounding boundary)
  - **ADR 0002 #13 completed** — `atoms_after_half_lives` past the normal range is the correctly rounded N₀/2ⁿ for every
    N₀ and n (M4c's was faithful for a fractional n with N₀ not a power of two): `_km_cr_exp2_scaled(N₀, −n)`, checked
    against MPFR (`mpfr_ref.scaled_exp2`, 2^f to 400 bits bounded from both sides) on 611,518 seeded pairs and 100,228
    through the accurate path alone on x86_64 ELF, x86_64 PE and aarch64 (`cr_verify.py`'s `exp2s`), on every golden row
    (`validate.py`: MPFR and a 60-digit 2⁻ᶠ), and on 4,057 pairs in `tests/mathx.tcyr`; the generator's notes are the
    correctly rounded values (`scripts/oracle/src/exact.rs`: 2⁻ᶠ to 320 bits in integer arithmetic, Ziv's two-sided test)
  - golden module `potential` (26,938 rows: class B 11,320, C 15,618): Lennard-Jones over r/σ 0.5–100 at Å and metre
    scale and the over- and underflow of (σ/r)⁶ and (σ/r)¹², r_min over every binade, subnormals and the overflow end,
    Morse across exp's range with **near-cancellation rows** (1 − e^(−a(r − r_e)) for |a(r − r_e)| from 10⁻¹⁷ to 1, S up to
    8.4 × 10¹⁸), the Morse width and Coulomb's law over every decade and their intermediates' over- and underflow,
    distances, bond and dihedral angles at Å, metre and femtometre scale, right / linear / tiny angles, trans / cis /
    gauche / perpendicular dihedrals with every sign of their zero coordinates (12 rows where atan2's −π becomes π),
    **atoms on the coincidence and collinearity tolerances to the ulp**, every Err path, ±0, NaN, ±∞, DBL_MAX and
    subnormal coordinates; 1 row has an ADR 0006 correction (`tests/golden/potential.cr.tsv`, glibc's atan2). `validate.py`
    recomputes every row, proves #16 (2·σ⁶ against each midpoint's sixth power), #17 and #18 exactly, and pins both
    tolerances (± 1 ulp, × ½, × 2), Rust's 4.0 and 2.0 and the π of the −π mapping
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `244803 ok, 0 failed, 0 measured, 0 pending, 4644 accepted (ADR 6:
    689), 0 glibc-matching (failed)` on 16 modules — under ADR 2, 1,041 #16, 180 #17 and 181 #18 rows besides the earlier
    ones (#13 now 1,810) — every row 0 ulp from its reference, every class-C function at bound 0
    ([`port-parity.md`](docs/development/port-parity.md) §14)
  - 33 Rust tests ported one for one (`tests/potential.tcyr`; `morse_width_parameter_matches_force_constant` computes
    spectroscopy's √(k/μ) inline until M4e) plus 5 Cyrius-only groups (null pointers and the point layout, the helpers'
    operation order, the tolerances to the ulp and the −π → π mapping, #16's exact r_min, #17 / #18's overflow edges):
    126 assertions; `tests/mathx.tcyr` gains `_km_cr_exp2_scaled`'s groups; the bundle smoke test gains potential; 19 test
    files, 0 warnings; surface parity 357/357; benches `potential/dihedral_angle` (against Rust) and the Cyrius-only
    `potential/bond_angle`, `lennard_jones`, `lj_equilibrium_distance`, `morse`, `mathx/cr_exp2_scaled`
  - generator: the exact integer arithmetic moved to `scripts/oracle/src/exact.rs` (with #16's sixth-power rounding);
    `cr_verify.py` gains `exp2s`; `mpfr_ref.py` gains `scaled_exp2`; `surface_parity.py` counts potential as ported
- **M4e spectroscopy** (port plan §4, §6.4 "Decided in M4e"; M4 in progress, stoichiometry next):
  - **spectroscopy** — the Beer–Lambert law (`absorbance`, `transmittance` on `_km_cr_exp10`,
    `absorbance_from_transmittance` on `_km_cr_log10`, `concentration_from_absorbance`), the photon conversions and eV ↔ J,
    `bohr_energy_level`, `bohr_transition_energy`, `rydberg_wavelength`, the spectral series (`enum KmSpectralSeries`,
    `spectral_series_n_lower`, `spectral_series_wavelength` with the reduced-mass R_H), the wavenumber conversions
    (ADR 0002 #19, Breaking, below), `reduced_mass`, `harmonic_oscillator_energy`, `vibrational_frequency`,
    `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia`, `isotope_abundances(z)` (a (pointer, count)
    slice into the baked table), `isotope_pattern(atoms, n, out, cap)`, `fft_spectrum(signal, n, out, cap)` and
    `find_peaks(data, n, t, out, cap)`, in Rust's check and operation order. Rust's u8 / u32 arguments are words, never
    wrapped (Err, NaN or the empty slice outside the type; an isotope-pattern Z outside u8 is `InvalidElement`); Vec results
    go to the caller's buffer of `cap` elements with `Ok(count)`, or `Err(KM_ERR_INVALID_INPUT)` when they don't fit (the
    buffer untouched)
  - **isotope_pattern** — Rust's pruned convolution operation for operation: a stable insertion sort on the masses'
    total order (Rust's `sort_by(total_cmp)` is stable, and the binned means depend on the order of equal masses), the
    binning, `fold(0.0, max)`, the renormalisation, `retain` in order, the `count` loop as Rust's; Rust's private
    thresholds 10⁻¹², 0.001 and 0.5 u typed inline in the helpers that apply them and pinned to the ulp by
    `tests/spectroscopy.tcyr`; the working tables start on the stack and move to the heap only for large molecules
    (`Err(KM_ERR_COMPUTATION_ERROR)` if that allocation fails)
  - **fft_spectrum** — hisab 1.4.0's `num::fft` transcribed as the private `_km_fft` (ADR 0003 amendment M4e): its
    bit-reversal loop, butterflies and `w *= wn` recurrence with Complex's operations in their order, no FMA, the twiddles
    from the correctly rounded `_km_cr_sincos` (hisab 3.4.1's `num_fft` takes them from the stdlib's `f64_cos` / `f64_sin`
    and reseeds them every 64 steps), and Complex::abs's √(re·re + im·im); the transform runs in `out` when cap ≥ 2n
    (nothing allocated; `signal` may be `out`), else in 16·n heap bytes
  - generated `src/spectroscopy_data.cyr`: `enum KmSpectralSeries` (`KM_SPECTRAL_SERIES_*`, append-only) and its n_lower
    column, the isotope abundance table from the oracle's `isotope_abundances(z)` for every u8 z (16 pairs over H, C, N,
    O, S, Cl, Br, with the literal's comments), Rust's private R_H read through the API (the Lyman line u32::MAX − 1 is
    RN(1/R_H), which fixes it uniquely) and its private mₑ/mₚ
  - **mathx** — `_km_cr_hundredth_div(x)`, 1/(100·x) rounded once (a double-double residual and a rounding test, an exact
    64-bit long division as its fallback; 1/100 generated as H + L by `mathx_tables.py`), for ADR 0002 #19, with
    `_fast` / `_acc` entry points: `cr_verify.py`'s `h100` (1,049,032 arguments, 100,018 through the integer path, against
    `mpfr_ref.hundredth_div`) and `tests/mathx.tcyr`'s two new groups
  - golden module `spectroscopy` (52,985 rows: class A 1,255, B 46,415, C 5,315): every u8 Z for the Bohr levels and the
    isotope table, transitions and the Rydberg formula over Z × (n₁, n₂) grids and close large n, every series to
    u32::MAX's `checked_add` overflow, the wavenumbers over every decade, power of two and binade, the isotope patterns of
    the Rust tests' molecules, every u8 Z alone, every element with data × counts to 5,000, elements without data to
    u32::MAX, mixtures and the Err position, the FFT at every length 0 … 16 and every power of two to 2¹³ (impulses,
    constants, ±0, Nyquist, sinusoids, random, huge, subnormal, NaN, ±∞), `find_peaks` over thresholds in and out of
    [0, 1] (NaN, ±∞, −0), plateaus and non-finite data; 118 rows have ADR 0006 corrections
    (`tests/golden/spectroscopy.cr.tsv`: glibc's pow(10, ·) on 3 transmittance rows, 10²³ among them, and its log10 on
    115 absorbance rows); `validate.py` recomputes every row (with its own FFT transcription; the release oracle's
    `sincos` gives glibc's sin and cos bits), proves #19 with `Fraction`, checks the generated data and pins the typed
    constants (`--sensitivity` too); `probe_tables.py` compiles the 81 spectroscopy_data cells
  - the harness: `_PAR_MAX` 1024 → 8448 tokens (the 2¹³-sample FFT rows); a row mapped to an ADR's exact value is counted
    as accepted even within its class's tolerance (`_par_force_map`: #19's one-ulp class-B rows); `fft_spectrum` replays
    through the heap work table up to n = 64 and in place above
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `296016 ok, 0 failed, 0 measured, 0 pending, 6416 accepted (ADR 6:
    807), 0 glibc-matching (failed)` on 17 modules — under ADR 2, 1,654 #19 rows besides the earlier ones — every row 0
    ulp from its reference, every class-C function at bound 0 ([`port-parity.md`](docs/development/port-parity.md) §15)
  - 62 Rust tests ported one for one plus the 2 spectroscopy integration tests (`tests/spectroscopy.tcyr`) and 6
    Cyrius-only groups (out-of-type integers, the caller buffers and slices, the isotope helpers' stable sort and
    thresholds to the ulp, the FFT transcription, #19's exact values, R_H): 289 assertions;
    `tests/potential.tcyr`'s `morse_width_parameter_matches_force_constant` now calls `vibrational_frequency`, as Rust's;
    20 test files, 0 warnings; surface parity 394/394; benches `spectroscopy/absorbance`, `bohr_energy_level`,
    `isotope_pattern_glucose` (against Rust) and the Cyrius-only `fft_spectrum_256` / `_1024` / `_4096`,
    `find_peaks_1000`, `wavenumber_from_wavelength`
  - generator: `m_spectroscopy.rs`, `t_spectroscopy.rs`, `exact::cr_hundredth_div` (128-bit integer division);
    `mpfr_ref.py` gains `hundredth_div`; `surface_parity.py` counts spectroscopy as ported (the SpectralSeries methods as
    `spectral_series_*`); `dist/kimiya.deps` lists `alloc` (spectroscopy calls it)
- **M4e-2 potential's geometry with an unbounded exponent** (ADR 0002 #20, the user's decision of 2026-11-09, "Yes,
  prescale"; port plan §6.4 "Decided in M4e-2"):
  - **potential** — `distance_3d`, `bond_angle` and `dihedral_angle` are Rust's operation sequence evaluated in binary64
    with an unbounded exponent range, the final √ or atan2 rounded once (Breaking, below). Each runs Rust's chain with a
    cheap test that proves it stayed in the normal range (a sum-of-squares window for `distance_3d` and `bond_angle`, which
    replaces `bond_angle`'s two `is_positive` calls; a per-component window on b₁, b₂, b₃ for the dihedral), and only
    where the test fails, for finite coordinates, re-evaluates with the private `_km_distance_3d_ux`,
    `_km_bond_angle_ux` and `_km_dihedral_angle_ux` (Rust's sub3 / cross / dot3 order on unbounded values)
  - **mathx** — binary64 with an unbounded exponent: (significand, exponent) pairs with `_km_ux_from`, `_km_ux_neg`,
    `_km_ux_add`, `_km_ux_sub`, `_km_ux_mul`, `_km_ux_div`, `_km_ux_sqrt` (each one IEEE operation on the significands,
    hence the exact 53-bit rounding), `_km_ux_cmp`, `_km_ux_le`, `_km_ux_max` (Rust's max), `_km_ux_is_positive`, the
    round-once `_km_ux_round` and the final `_km_ux_add_f64`, `_km_ux_sub_f64`, `_km_ux_mul_f64`, `_km_ux_div_f64`,
    `_km_ux_sqrt_f64` and `_km_ux_atan2_f64` (`_fast`, `_acc`), each rounded once into binary64 with its exact error
    breaking subnormal midpoints: the machinery of the library-wide audit of intermediate over- and underflow after M4
  - golden module `potential`: 31,761 rows (4,823 new edge rows: coordinates near 2^±(500 … 1074), mixed magnitudes and
    subnormal coordinates, differences that overflow, the coincidence and collinearity tolerances at extreme scales with
    both crossings found by bisection, b₁ tiny against b₂ huge, tiny angles with subnormal results, trans and cis
    dihedrals at scale for the −π → π mapping, seeded geometries at random scales); 3,605 rows noted with the port's
    result (854 `distance_3d`, 1,205 `bond_angle`, 1,546 `dihedral_angle`, 91 of them the port's Err), 314 of them among
    the 26,938 earlier rows; 2 more ADR 6 rows (glibc's atan2 rounds the subnormal midpoint quotient 1.5·2⁻¹⁰⁷⁴ up)
  - generator: `exact::Ux` (53-bit rounding of u128 mantissas at any exponent) and `exact::ux_atan2_f64` (a 400-bit
    fixed-point atan by Euler's series, Ziv's two-sided test), the f64 replica of Rust's chains checked on every row;
    `mpfr_ref.py` gains `UX`, `unbounded`, the final `ux_sqrt_final`, `ux_atan2_final`, `ux_mul_final`, `ux_div_final`,
    `ux_add_final` and `potential_*` (Rust's sequences on UX); `validate.py` proves every noted row, every unnoted row's
    expected result and the in-range rows' Rust bits against MPFR, and the scale invariance; `cr_verify.py` gains
    `uxsqrt`, `uxmul`, `uxdiv`, `uxadd`, `uxatan2` (1,040,039 arguments) and `dist`, `bond`, `dih` (114,000 edge
    geometries, the slow path forced on all)
  - the harness maps a #20 row's `(f:…)` or `(E:InvalidInput)` note (`_par_got_noted20`)
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `297232 ok, 0 failed, 0 measured, 0 pending, 10023 accepted (ADR 6:
    809), 0 glibc-matching (failed)` on 17 modules — under ADR 2, 3,605 #20 rows besides the earlier ones
    ([`port-parity.md`](docs/development/port-parity.md) §16)
  - tests: `tests/mathx.tcyr` gains `ux_values` (48 MPFR values past binary64's range: subnormal ties, sticky terms,
    overflow, exponent differences past 2^±2000; the signed zeros, IEEE's specials, order and max) and `ux_ieee` (20,000
    seeded pairs: the final operations against IEEE's own, shifted past the range and back, atan2 against
    `_km_cr_atan2`): 1,546 assertions;
    `tests/potential.tcyr` gains `cyrius_unbounded_edges` (MPFR's values), `cyrius_scale_invariance` (600 random
    geometries × 75 scales, bit for bit) and `cyrius_slow_path_agrees` (5,000 in-range geometries, the slow path forced):
    159 assertions; 20 test files, 0 warnings; surface parity 394/394;
    benches `potential/distance_3d` and the slow path's `potential/distance_3d_2p600`, `bond_angle_2p600`,
    `dihedral_angle_2p600`
- **M4e-3 exact Rydberg rationals and the exact reduced mass** (ADR 0002 #21 and #22, the user's decisions of 2026-11-10,
  "Compute exactly" and "Round once"; port plan §6.4 "Decided in M4e-3"):
  - **spectroscopy** — `bohr_energy_level`, `bohr_transition_energy`, `rydberg_wavelength` and
    `spectral_series_wavelength` return their exact rationals, the constant as its double, rounded once; `reduced_mass`
    the correctly rounded m₁·m₂/(m₁ + m₂) for every pair of positive finite masses (Breaking, below). Each has a
    double-double fast path with a proven error bound and a rounding test (`_km_bohr_level_fast`,
    `_km_bohr_transition_fast`, `_km_rydberg_fast`, `_km_reduced_mass_fast`; the operands exact for n₂ < 2¹³ and for the
    reduced mass) and an exact slow path (`_km_bohr_level_acc`, `_km_bohr_transition_acc`, `_km_rydberg_acc`,
    `_km_reduced_mass_acc`), which decides the ties; the checks, Err codes and their order are Rust's
  - **mathx** — multi-word unsigned integers (8 limbs of 32 bits: `_km_bn_set`, `_copy`, `_bitlen`, `_cmp`, `_add`, `_sub`,
    `_mul_small`, `_mul`, `_shl`), `_km_f64_int` and the general `_km_cr_ratio(p, q, e)` = RN(p/q·2^e) (restoring long
    division to a 55-bit quotient and a sticky remainder, rounded once by `_km_cr_from_int`: exact, subnormal results and
    ties included), for reuse by M4f and the overflow audit; the double-double product `_km_dd_mul` (to 2⁻¹⁰²·⁹⁹) and
    quotient with a rounding test `_km_dd_div_cr` (its error budget in the source)
  - golden module `spectroscopy`: 66,691 rows (13,706 new: every exact Bohr level, n = 2^k·Z; n₁ = n₂ − 1 at every binary
    scale and next to u32::MAX; the fast path's 2¹³ edge; the 13 transition ties and each exact-value family's smallest,
    middle and largest member from an exhaustive search of the domain; the series' large lines; for the reduced mass, m = m
    at every power of two and subnormal multiple, three families of normal ties and two of subnormal ones at scales across
    the range, exponent gaps 45 … 64 around the proven 55, the 2⁻¹⁰²¹ edge, both range ends, binade-uniform pairs);
    9,454 rows noted with the port's correctly rounded value (#21: 3,823 `bohr_energy_level`, 1,351
    `bohr_transition_energy`, 1,376 `rydberg_wavelength`, 67 `SpectralSeries::wavelength`; #22: 2,837 `reduced_mass`)
  - generator: `exact::cr_ratio` and `cr_bohr_level`, `cr_transition`, `cr_rydberg`, `cr_reduced_mass` (Big integers, no
    cutoff); the exhaustive search of `bohr_transition_energy`'s exact values (`m_spectroscopy::transition_exact_values`,
    trial division and Miller–Rabin); `validate.py` proves every Ok row of the five functions with `Fraction`, re-derives
    the structure (R's odd significands factored, the search with its own factorisation) and checks every exact level
    and tie is in the corpus; `cr_verify.py` gains `ratio`, `bohr`, `trans`, `ryd`, `series` and `rmass` (2.25 million
    arguments; the public function, the fast path and the exact path forced) against `Fraction`
  - the harness maps a #21 or #22 row's `(f:…)` note and counts #21 per function
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `301484 ok, 0 failed, 0 measured, 0 pending, 19477 accepted (ADR 6:
    809), 0 glibc-matching (failed)` on 17 modules — under ADR 2, 6,617 #21 and 2,837 #22 rows besides the earlier ones
    ([`port-parity.md`](docs/development/port-parity.md) §17)
  - tests: `tests/spectroscopy.tcyr` gains `cyrius_adr21` (−R at every Z = n and −R/4^k at n = 2^k·Z, the exact 16R and 32R,
    the 13 ties declined by the fast path and rounded to even, Lyman α, Hα, Rust's worst cancellation, 6,000 seeded
    (Z, n₁, n₂) fast against exact) and `cyrius_adr22` (m/2 on 10,000 seeded masses, the normal and subnormal ties,
    (DBL_MAX, DBL_MAX), the gap of 55, 20,000 seeded pairs fast against exact): 311 assertions; `tests/mathx.tcyr` gains
    `cr_ratio_values` (the multi-word operations, ratios with ties, subnormal results and the overflow threshold, the
    double-double helpers) and `cr_ratio_agrees` (10,000 seeded quotients, exact against double-double): 1,566 assertions;
    20 test files, 0 warnings; benches `spectroscopy/bohr_transition_energy`, `rydberg_wavelength`,
    `series_wavelength_balmer_alpha`, `reduced_mass_hcl` and the exact path's `bohr_transition_energy_tie`,
    `reduced_mass_tie`
- **M4e-4 the cited isotope table** ([ADR 0002 #23](docs/adr/0002-port-api-divergences.md), the user's decision of
  2026-11-10, "Upgrade the data"; the data change itself is under *Changed*):
  - `scripts/oracle/isotopes.tsv`: the source records of the isotope table — per isotope Z, A, symbol, the AME2020 mass
    in u with its published digits, the TICE 2009 representative abundance and its printed value, and the AME2020 line
    verbatim — under a header with the URLs, sha256s and references; and CIAAW's standard-atomic-weight intervals of H
    and S
  - generator: `t_spectroscopy::isotope_table` reads it (each mass re-read from its line by the file's Fortran format)
    and fails unless its isotopes are the oracle's plus ³⁶S alone, every shared abundance is Rust's literal, every Rust
    mass is the AME2020 one to five decimals and each element sums to exactly 1 (decimal arithmetic);
    `m_spectroscopy::pattern_replica` (Rust's `isotope_pattern` with the table a parameter, checked against the oracle on
    every row) gives the port's results for the notes; `tsv::tokens_note` writes a whole result as a note
  - golden module `spectroscopy`: the rows keep the oracle's values; 7 `isotope_abundances` and 347 `isotope_pattern`
    rows noted with the port's result, `[ADR 0002 #23: … (v:n f:… …)]` (the note's token list extended from one value to
    a whole result); the harness compares a noted row's result with those tokens bit for bit (`_par_map_toks`) and
    counts it per function; `validate.py` repeats the four checks on its own parse, checks the generated table and the
    mean masses against CIAAW's intervals, proves every note with its replica and that every unnoted row is unchanged,
    and pins the thresholds on the port's results too; `probe_tables.py` probes the cells against the notes
  - `scripts/oracle/isotope_sources.py`: re-verifies the file against fresh downloads (the sha256s, every AME2020 line
    verbatim, every TICE value a cell of the PDF's Table 1 in order, the CIAAW intervals)
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `301130 ok, 0 failed, 0 measured, 0 pending, 19831 accepted (ADR 6:
    809), 0 glibc-matching (failed)` on 17 modules — 354 rows move from ok to accepted under ADR 2 (#23)
    ([`port-parity.md`](docs/development/port-parity.md) §18)
  - tests: `tests/spectroscopy.tcyr`'s `hydrogen_abundances_iupac_2009` asserts the table's pairs and sums sulfur too;
    `cyrius_adr23_table` (every element's count, ¹H's full mass, sulfur's four pairs, their sum) and `cyrius_adr23_s8`
    (S₈'s seven peaks bit for bit against `validate.py`'s replica, ³⁶S raising M+4 … M+6): 350 assertions; benches
    `spectroscopy/isotope_pattern_cystine` and `isotope_pattern_s8`

- **M4f stoichiometry** (port plan §4, §6.4 "Decided in M4f"; [ADR 0002 #7](docs/adr/0002-port-api-divergences.md),
  realised; M4 complete):
  - **stoichiometry** (`src/stoichiometry.cyr`, in the `[lib]` order between kinetics and organic_data) —
    `balance_equation(matrix, m, n, out, cap)` on an integer matrix (flat row-major; Rust's fraction-free Gauss–Jordan
    elimination with checked i64 arithmetic, its gcd reductions, the null vector from the RREF, the sign and u32 checks
    and the i128 verification, every Rust overflow an Err at the same step), `assign_oxidation_states(atoms, n, charge,
    out, cap)` with exact rational states (24-byte (Z, numerator, denominator) records), `limiting_reagent`,
    `theoretical_yield_moles`, `percent_yield`, `percent_composition(molecule, out, cap)`, `empirical_formula`; Rust's
    errors, codes and order; the private helpers `_km_i64_{add,sub,mul,neg}_ck`, `_km_gcd_u64`, `_km_reduce_row`, the
    signed 128-bit `_km_i128_mul_u32` / `_km_i128_add_ck`, `_km_exact_oxidation_state`, `_km_rule_oxidation_state`,
    `_km_is_metal`, `_km_empirical_close`
  - generator `m_stoichiometry.rs`: golden module `stoichiometry` (41,435 rows: `balance_equation` 3,329, A —
    the Rust tests, CH₄ / CuSO₄ + Fe / KMnO₄ + HCl in every row order, textbook and ionic equations, every combustion
    C_aH_bO_c with a ≤ 10, chains to 70 species, 820 seeded consistent and 700 random matrices, every error class, the
    10⁻⁹ tolerance to the ulp, ±2⁵³, intermediate overflow; `assign_oxidation_states` 23,498, A — the 1.2 tests, every u8
    Z alone, every element's oxides, hydrides and fluorides, (su)peroxides and ozonides, the u32 / i32 extremes, seeded
    species; `limiting_reagent` 2,897, A; `theoretical_yield_moles` 2,965, B; `percent_yield` 4,071, B;
    `percent_composition` 1,292, B; `empirical_formula` 3,383, A); `m_stoichiometry::port_oxidation`, Rust's algorithm on
    exact rationals, checked against the oracle on every row, writes ADR 0002 #7's notes
  - harness: the stoichiometry dispatcher; the matrix converted as Rust's test converts it (`_par_in_matrix`) and replayed
    on integers; 9,090 oxidation rows noted with the port's exact states and 78 balancing rows with no integer
    counterpart mapped under ADR 2, counted; `_par_map_toks(adr, tag)` takes the note's tag
  - `validate.py` §18: every row recomputed from Rust's semantics (checked i64 / i128 on Python integers, IEEE doubles in
    Rust's order), the port's states with `Fraction` and every note proven, every balancing row checked by exact linear
    algebra (1,976 Ok rows the primitive null vector; 104 Err rows where Rust's i64 elimination overflows although a valid
    u32 solution exists), EMPIRICAL_TOLERANCE, INTEGER_TOLERANCE and MAX_COMPOSITION_ENTRY pinned
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `333397 ok, 0 failed, 0 measured, 0 pending, 28999 accepted (ADR 6:
    809), 0 glibc-matching (failed)` on 18 modules ([`port-parity.md`](docs/development/port-parity.md) §19)
  - tests: `tests/stoichiometry.tcyr`, the 63 Rust tests (7 adapted for ADR 0002 #7) and 6 Cyrius-only groups (the
    rational states and their i32 range, the integer matrix and its 2⁵³ bound, the checked i64 and i128 helpers at their
    edges, out-of-type words, slices and buffers, the tolerance to the ulp): 265 assertions; `tests/kimiya.tcyr` balances
    CH₄ and assigns Fe₃O₄ through the bundle; surface parity 401/401; benches `stoichiometry/balance_ch4_combustion`,
    `oxidation_states_kmno4` (against Rust) and seven Cyrius-only
  - the five float functions are Rust's arithmetic unchanged; their measured misses of the exact values are listed in
    ADR 0002's M4f addition for the user's decision
- **M4f-2 exact stoichiometry** ([ADR 0002 #24, #25 and #26](docs/adr/0002-port-api-divergences.md), the user's
  decisions of 2026-11-10, "Round once", "Compare exactly" and "Evaluate exactly"; port plan §6.4 "Decided in M4f-2"; the
  behaviour changes are under *Breaking*):
  - **stoichiometry** — `percent_yield`, `theoretical_yield_moles` and `percent_composition` return their exact rationals
    rounded once through `_km_cr_mul_div` = RN(x·y/z): a double-double fast path on the significands with
    `_km_dd_div_cr`'s rounding test and an exact resolution of the midpoint where the test fails (`_km_mul_div_tie`: the
    residual's sign as four 64-bit integers), and the exact `_km_mul_div_acc` (`_km_cr_ratio`) for subnormal operands and
    results that may leave the normal range; `limiting_reagent` keeps Rust's rounded ratios where they differ (RN is
    monotone) and compares equal ones exactly (`_km_ratio_lt`: the two cross products, exponents first, then exact
    TwoProducts); `empirical_formula` with exact moles, minimum, ratios (double-doubles with `_km_dd_quot`), the proven
    2³² − ½ cut, a tolerance test that is proven or defers to the exact one (`_km_emp_close_fast`, `_km_emp_close_exact`:
    integers below 2¹⁴³, a 36-bit restoring division) and exact counts; `_km_empirical_formula(…, force)` forces the exact
    tests (`cr_verify.py`). Rust's checks, Err codes and their order are kept
  - generator: `exact::cr_mul_div`, `ratio_lt` and `empirical_exact` (`Big` integers, no cutoff),
    `m_stoichiometry::port_limiting` and `port_empirical` (Rust's validation, asserted against the oracle), the notes;
    `emit_exact` appends 6,341 hard rows with their own seed (`percent_yield`'s ties at RN, subnormal ones included, and
    range ends; the yield's over- and underflowing m/c and its midpoints; `limiting_reagent`'s exact ties, near-ties, both
    ratios over- or underflowing; `empirical_formula`'s knife edges at ±1/10, subnormal percentages and the 2³² edge);
    `tsv::same_toks`
  - golden module `stoichiometry`: 47,776 rows; the oracle's values kept, 3,154 rows noted with the port's result —
    #24: 1,126 `percent_yield`, 240 `theoretical_yield_moles`, 598 `percent_composition`; #25: 920 `limiting_reagent`;
    #26: 270 `empirical_formula` (of M4f's 3,383 rows, 14 flip, the decision's count); the harness maps them
    (`_par_got_noted`, `_par_map_exact` with an index, `_par_map_toks`) and counts them per function
  - `validate.py` §18: the port's result of every row of the five recomputed with `Fraction` and each note proven; the 2³²
    claim and the kept ratios' range checked on the corpus; Rust's empirical Err paths tabulated against the port's results
    (ADR 0002's M4f-2 addition); the hard rows counted; the exact 1/10 pinned
  - `cr_verify.py`: `muldiv` (640,000 triples), `pyield`, `tyield`, `pcomp`, `limit`, `empir` — 2.23 million arguments
    against `Fraction`, the exact paths forced, 0 mismatches and the same bits on x86_64 ELF, x86_64 PE and aarch64
  - **result**, x86_64 ELF, x86_64 PE and aarch64: `336584 ok, 0 failed, 0 measured, 0 pending, 32153 accepted (ADR 6:
    809), 0 glibc-matching (failed)` on 18 modules ([`port-parity.md`](docs/development/port-parity.md) §20)
  - tests: `tests/stoichiometry.tcyr`'s `cyrius_empirical_tolerance` (Rust's 0.1 on `_km_empirical_close`) becomes
    `cyrius_adr26_exact_tolerance` (the exact edges 11/10 and 19/10 and one ulp past them, 2.55·2, C 12.011 % with
    O 17.5989 %, the count's 2³² edge, subnormal percentages, the fast test's three outcomes, the exact test), with
    `cyrius_adr24_rounded_once` (the subnormal and normal ties, 100/3, the overflow end, Rust's ±0, the yield's over- and
    underflowing m/c, water's and methane's percentages, the midpoint resolution above and below) and
    `cyrius_adr25_exact_comparison` (both ratios over- or underflowing, 0 against an underflow, ties and near-ties,
    `_km_ratio_lt`'s branches); no Rust test's assertion changes (`empirical_formula_nan_inf_are_errors`' subnormal case
    is still the ComputationError, its ratio about 3.75·10³²⁰): 317 assertions; benches `stoichiometry/theoretical_yield_moles_tie`,
    `percent_yield_subnormal`, `limiting_reagent_tie`, `empirical_formula_edge`
- **The library-wide overflow audit, B1** ([ADR 0002 #27 and #28](docs/adr/0002-port-api-divergences.md), the user's
  decisions "Library-wide audit" of 2026-11-09 and "Round once" of 2026-11-10; port plan §6.4 "Decided in the overflow
  audit's B1"; the behaviour changes are under *Breaking*):
  - **mathx** — the correctly rounded monomial `_km_cr_mono(num, nn, den, nd)` = RN(Π num / Π den) (up to 6 and 4
    finite factors; signs, IEEE's signed zeros, ties to even, subnormal results, ±∞ past DBL_MAX): M4f-2's `_km_cr_mul_div`
    moved from `stoichiometry.cyr` and generalised — a double-double fast path (`_km_dd_mul_d`, `_km_dd_div_ends` and
    `_km_dd_div_cr`'s rounding test) that takes subnormal operands and resolves midpoints exactly for up to three
    numerator and two denominator factors (`_km_mono_tie`, the sign of an exact sum by Shewchuk's Grow-Expansion,
    `_km_exp_sign`), and the exact `_km_cr_mono_acc` on `_km_bn_*` / `_km_cr_ratio`; the shapes the modules call
    (`_km_cr_mul_div` — the stoichiometry calls, the same results —, `_km_cr_mul3`, `_km_cr_div2`, `_km_cr_mul3_div`,
    `_km_cr_mul5`) each have a straight-line fast path (`_km_mul_div_fast` with M4f-2's `_km_mul_div_tie`,
    `_km_mul3_fast`, `_km_div2_fast`, `_km_mul3_div_fast`, `_km_mul5_fast`; `_km_dd_round_ends`, `_km_mono_fin`) and
    fall back to the general one
  - **thermo**, **phase** — `heat_transfer`, `specific_heat_from_calorimetry`, `boiling_point_elevation`,
    `freezing_point_depression` and `osmotic_pressure` rounded once (#27)
  - **kinetics** — `zero_order_half_life`, `steady_state_intermediate`, `pre_equilibrium_rate` rounded once;
    `second_order_concentration` and `michaelis_menten` under #20's semantics, Rust's chain behind a range test
    (`_km_second_order_fast`, `_km_michaelis_menten_fast`) and an unbounded-exponent slow path (`_km_second_order_ux`,
    `_km_michaelis_menten_ux`) (#28)
  - generator: `scripts/oracle/src/audit.rs` (the port's values independently: `exact::cr_monomial`, `exact::Ux`'s new
    `div` and `div_f64`; Rust's validation replicated and asserted; the notes; each module's B1 edge rows appended with
    their own seed: thermo +537, phase +528, kinetics +1,064)
  - golden modules `thermo` 1,644, `phase` 4,174, `kinetics` 31,210 rows; the oracle's values kept, 1,318 rows noted
    with the port's value — #27: 151 `heat_transfer`, 171 `specific_heat_from_calorimetry`, 82 + 82 colligative, 196
    `osmotic_pressure`; #28: 89 `second_order_concentration`, 90 `zero_order_half_life`, 118
    `steady_state_intermediate`, 179 `pre_equilibrium_rate`, 160 `michaelis_menten`; the harness maps them
    (`_par_got_noted`, the new `_par_got_noted_f`) and counts them per function
  - `validate.py` §19: every row of the ten functions recomputed (`Fraction`; `mpfr_ref.UX` for the two #20-style ones),
    each note proven, no in-range row of the #20-style ones noted, `zero_order_half_life` noted only past 2k's overflow,
    the hard rows counted
  - `cr_verify.py`: `mono`, `monotie` and the ten functions (`heat`, `shc`, `boil`, `frz`, `osm`, `zhl`, `ssi`, `pre`,
    `soc`, `mm`: the public function, the fast path, the slow path forced) against `Fraction` and MPFR on x86_64 ELF,
    x86_64 PE and aarch64
  - tests: `tests/mathx.tcyr` `cr_mono_values` (36 exact cases), `cr_mono_sign_and_tie`, `cr_mono_agrees` (20,000 seeded
    arguments, the fast and exact paths); `tests/thermo.tcyr` `cyrius_adr27_rounded_once`, `tests/phase.tcyr`
    `cyrius_adr27_rounded_once`, `tests/kinetics.tcyr` `cyrius_adr28` (the survey's edges, a 1-ulp in-range case each,
    5,000 seeded inputs against the exact or slow paths); no Rust test changes; one Cyrius-only kinetics assertion of
    Rust's `second_order_concentration(2⁻¹⁰⁷⁴, 0, 1)` = +0 now pins [A]₀; benches for the eight functions without a Rust
    bench and four edge benches
- **The library-wide overflow audit, B2** ([ADR 0002 #29 and #30](docs/adr/0002-port-api-divergences.md), B1's decisions
  unchanged; port plan §6.4 "Decided in the overflow audit's B2"; the behaviour changes are under *Breaking*):
  - **mathx** — four more straight-line monomial shapes on B1's `_km_cr_mono` (`_km_mul2_div2_fast`,
    `_km_mul3_div2_fast`, `_km_mul4_fast`, `_km_div3_fast` with their wrappers `_km_cr_mul2_div2`, `_km_cr_mul3_div2`,
    `_km_cr_mul4`, `_km_cr_div3`), and the correctly rounded √ of a quotient `_km_cr_sqrt_div(x, y)` = RN(√(x/y)): a
    double-double fast path (`_km_sqrt_div_fast`, proven within 5.7u² with `_km_dd_div_cr`'s rounding test) and an
    exact integer square root on `_km_bn_*` (`_km_sqrt_div_acc`); √(x/y) is never a rounding boundary of binary64
    (proven)
  - **solution** — `dilution` and `solution_conductivity` rounded once (#29)
  - **electrochemistry** — `faraday_mass_deposited`, `faraday_charge_required` and `faraday_mass_ratio` rounded once;
    `nernst_potential` under #20's semantics, Rust's chain behind a range test (`_km_nernst_fast`) and an
    unbounded-exponent slow path (`_km_nernst_ux`), the correctly rounded ln Q its primitive; `nernst_potential_25c`
    unchanged (#29)
  - **spectroscopy** — `absorbance`, `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant` and
    `diatomic_moment_of_inertia` rounded once; `vibrational_frequency` the correctly rounded √(k/μ) (#30)
  - generator: `scripts/oracle/src/audit.rs`'s B2 (the port's values independently: `exact::cr_monomial`, the new
    `exact::cr_sqrt_ratio` and `Ux::add_f64` / `sub_f64`; Rust's validation replicated and asserted; the notes; each
    module's B2 edge rows appended with their own seed: solution +458, electrochemistry +823, spectroscopy +1,178)
  - golden modules `solution` 29,024, `electrochemistry` 19,397, `spectroscopy` 67,869 rows; the oracle's values kept,
    4,862 rows noted with the port's value — #29: 94 `dilution`, 148 `solution_conductivity`, 904
    `faraday_mass_deposited`, 644 `faraday_charge_required`, 847 `faraday_mass_ratio`, 119 `nernst_potential`; #30: 195
    `absorbance`, 219 `concentration_from_absorbance`, 706 `vibrational_frequency`, 344 `rigid_rotor_energy`, 524
    `rotational_constant`, 118 `diatomic_moment_of_inertia`; the harness maps and counts them per function
  - `validate.py` §20: every row of the twelve functions recomputed (`Fraction`; `math.isqrt` and MPFR's √ of a 256-bit
    quotient for the √; `mpfr_ref.UX` with MPFR's log and the new `ux_sub_final` for Nernst), each note proven, no
    in-range Nernst row noted, every `nernst_potential_25c` row's chain normal, no √ row a rounding boundary, the hard
    rows counted; `mpfr_ref.py` gains `ux_sub_final` and `ux_sqrt_div_final`, with self-checks
  - `cr_verify.py`: `sqrtdiv` and the twelve functions (`dil`, `cond`, `fmass`, `fcharge`, `fratio`, `nernst`, `abs`,
    `conc`, `vib`, `rotor`, `rot`, `moment`: the public function, the shape's fast path, the slow path forced, roots
    within 2⁻¹⁰⁰ of a midpoint included) on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/mathx.tcyr` `b2_values` (the √'s exact roots, range ends, subnormal result, overflow, domain and
    near-midpoint roots; the shapes' ties and range ends), `b2_agrees` (20,000 seeded arguments: the four shapes against
    `_km_cr_mono_acc`, the √ against its exact path); `tests/solution.tcyr`, `tests/electrochemistry.tcyr`
    `cyrius_adr29` and `tests/spectroscopy.tcyr`'s ADR 0002 #30 group (the survey's edges, a 1-ulp in-range case per
    function, HCl's ω, 5,000 seeded inputs against the exact or slow paths); no Rust test changes; Cyrius-only benches
    for the ten functions without a Rust bench and four edge benches
- **The library-wide overflow audit, B3** ([ADR 0002 #31 and #32](docs/adr/0002-port-api-divergences.md), B1's decisions
  unchanged; port plan §6.4 "Decided in the overflow audit's B3"; the behaviour changes are under *Breaking*):
  - **mathx** — the correctly rounded `_km_cr_mul_sqrt_ratio(w, x, y, e)` = RN(w·√((x/y)·2^e)): B2's double-double root
    times w with a rounding test (`_km_mul_sqrt_ratio_fast`, within 9.8u²) and an exact integer root on `_km_bn_*`
    (`_km_mul_sqrt_ratio_acc`), which rounds the midpoints the value can be (3w·2^k) to even; `_km_ux_normal`, the
    #20-style range test (|x| in [2⁻¹⁰²¹, DBL_MAX])
  - **gas** — `ideal_gas_pressure`, `ideal_gas_volume` and `gas_density` rounded once; `van_der_waals_pressure`,
    `peng_robinson_pressure` and `redlich_kwong_pressure` under #20's semantics, Rust's chain behind a range test
    (`_km_vdw_fast`, `_km_pr_pressure_fast`, `_km_rk_pressure_fast`) and an unbounded-exponent slow path with the
    co-volume test on the unbounded co-volume (`_km_vdw_ux`, `_km_pr_pressure_ux`, `_km_rk_pressure_ux`; Tc^2.5 exactly
    rescaled, `_km_rk_tc25_ux`) (#31)
  - **potential** — `coulomb` and `coulomb_force` rounded once; `morse_width_parameter` the correctly rounded
    ω·√(μ/(2·D_e)); `lennard_jones` under #20's semantics (`_km_lj_fast`, `_km_lj_ux`); M4d's #17 and #18 folded in (#32)
  - **inorganic** — `born_lande_lattice_energy` under #20's semantics (`_km_born_lande_fast`, `_km_born_lande_ux`) (#32)
  - **nuclear** — `binding_energy_per_nucleon` under #20's semantics (`_km_bea_ux` past Δm·931.494…'s overflow) (#32)
  - generator: `scripts/oracle/src/audit.rs`'s B3 (the port's values independently: `exact::cr_monomial`, the new
    `exact::cr_sqrt_monomial`, `Ux::mul_f64` and `Ux::pow25` (Tc^2.5 by an integer square root); Rust's validation
    replicated and asserted, the co-volume tests the only allowed difference; the notes, `(E:InvalidInput)` where the
    port refuses; each module's B3 edge rows appended with their own seed: gas +1,680, potential +1,703, inorganic +254,
    nuclear +167); `m_potential.rs`'s #17 and #18 notes become #32's
  - golden modules `gas` 22,773, `potential` 33,464, `inorganic` 15,416, `nuclear` 54,224 rows; the oracle's values kept,
    5,339 rows noted with the port's result — #31: 146 `ideal_gas_pressure`, 148 `ideal_gas_volume`, 222 `gas_density`,
    206 `van_der_waals_pressure`, 201 `peng_robinson_pressure`, 234 `redlich_kwong_pressure` (98 of them the port's
    Err); #32: 582 `lennard_jones`, 1,860 `morse_width_parameter`, 730 `coulomb`, 679 `coulomb_force`, 175
    `born_lande_lattice_energy`, 156 `binding_energy_per_nucleon`; the harness maps and counts them per function (the new
    `_par_got_noted_r` for a noted Err), #17's and #18's regions among #32's
  - `validate.py` §21: every row of the twelve functions recomputed (`Fraction`; `math.isqrt` and MPFR's √ of a 320-bit
    quotient for the Morse width; `mpfr_ref.UX` with MPFR's own pow for Tc^2.5), each note proven, no in-range row of the
    six #20-style functions noted, #17 kept bit for bit and #18's rows counted, Tc^2.5's rescale checked on every row,
    glibc's pow(Tc, 2.5) correctly rounded on every Redlich–Kwong row, the hard rows counted; §16's #17 / #18 checks fold
    into it; `mpfr_ref.py` gains `UX.pow` and `ux_mul_sqrt_div_final`, with self-checks
  - `cr_verify.py`: `wsqrt` and the twelve functions (`igp`, `igv`, `dens`, `vdw`, `prp`, `rkp`, `lj`, `width`, `coul`,
    `force`, `lande`, `bea`: the public function, the fast path, the slow path forced, the midpoints and co-volume flips
    included) on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/mathx.tcyr` `b3_values` (the root's ties at both grids, the Morse width's edges, the domain,
    `_km_ux_normal`'s edges) and `b3_agrees` (20,000 seeded arguments against the exact path); `tests/gas.tcyr`
    `cyrius_adr31`, `tests/potential.tcyr` `cyrius_adr32` (and `cyrius_exact_overflow_edges` updated for the fold: the
    Morse width at D_e = DBL_MAX and μ = 2⁻¹⁰⁷⁴ is now correctly rounded), `tests/inorganic.tcyr`
    `cyrius_adr32_born_lande`, `tests/nuclear.tcyr` `cyrius_adr32_bea` (the survey's edges, 1-ulp in-range cases, the
    co-volume flips, 5,000 seeded inputs against the slow paths); no Rust test changes; Cyrius-only benches for the nine
    functions without a Rust bench and eight edge and tie benches
- **The library-wide overflow audit, B4** ([ADR 0002 #33 and #34](docs/adr/0002-port-api-divergences.md), B1's decisions
  unchanged; port plan §6.4 "Decided in the overflow audit's B4"; the behaviour changes are under *Breaking*, #29's
  correction under *Fixed*):
  - **thermo** — `hess_law` and `enthalpy_from_formation` under #20's semantics: Rust's sums in Rust's order, from −0.0,
    with an unbounded exponent (`_km_sum_ux`, `_km_sum_ux_f64`, `_km_formation_ux`), the final operation rounded once,
    behind the pure-sum fast test "result finite" (a sum of doubles is exact below 2⁻¹⁰²¹, so only overflow leaves the
    range; proven in the ADR) (#33)
  - **thermochem** — `reaction_enthalpy`, `reaction_gibbs_energy` (#1's × 1000 the last operation) and
    `reaction_entropy` likewise (`_km_sum_property` reports whether every product agreed; `_km_sum_property_ux`,
    `_km_reaction_j_ux`, `_km_reaction_s_ux`) (#33); `KmShomate.cp` and `enthalpy_antiderivative` under #20's semantics
    on Rust's evaluation order, a box test first (`_km_shomate_box`: |T| in [2⁻³⁰, 2³⁰], the coefficients in [2⁻²⁰⁰,
    2²⁰⁰], every tabulated segment), then per-operation tests (`_km_shomate_cp_fast`, `_km_shomate_h_fast`) and the slow
    paths (`_km_shomate_cp_ux`, `_km_shomate_h_ux`); the data-backed callers share the evaluation and stay in the box
    (#34)
  - **nuclear** — `q_value` (`_km_q_value_ux`); **solution** — `ionic_strength` (`_km_ionic_strength_ux`),
    `kohlrausch_molar_conductivity` (`_km_kohlrausch_fast`, `_km_kohlrausch_ux`), `limiting_molar_conductivity`
    (`_km_limiting_fast`, `_km_limiting_ux`); **inorganic** — `born_haber_lattice_energy` (`_km_born_haber_ux`);
    **reaction** — `gibbs_free_energy` (`_km_gibbs_ux`); **kinetics** — `zero_order_concentration` (`_km_zero_order_ux`,
    the clamp on the unbounded difference) (#33)
  - generator: `scripts/oracle/src/audit.rs`'s B4 (the port's values independently: `exact::Ux` in Rust's summation
    orders; Rust's validation replicated and asserted; #33's notes after #1's for the two J/mol sums; each module's B4
    edge rows appended with their own seed: thermo +184, reaction +172, nuclear +88, solution +430, inorganic +170,
    kinetics +171, thermochem +1,292; electrochemistry +23 for #29's correction)
  - golden modules `thermo` 1,828, `reaction` 1,730, `nuclear` 54,312, `solution` 29,454, `inorganic` 15,586, `kinetics`
    31,381, `thermochem` 11,149, `electrochemistry` 19,420 rows; the oracle's values kept, 924 rows noted with the
    port's value — #33: 28 `hess_law`, 27 `enthalpy_from_formation`, 115 `reaction_enthalpy`, 113
    `reaction_gibbs_energy`, 66 `reaction_entropy`, 58 `q_value`, 13 `ionic_strength`, 74 `born_haber_lattice_energy`,
    57 `gibbs_free_energy`, 54 `kohlrausch_molar_conductivity`, 31 `limiting_molar_conductivity`, 61
    `zero_order_concentration`; #34: 185 `cp`, 25 `enthalpy_antiderivative`; #29: 17 more `nernst_potential`; the
    harness maps and counts them per function
  - `validate.py` §22: every row of the fourteen functions recomputed through `mpfr_ref.UX` and the new `ux_sum` /
    `ux_sum_final` (Rust's summation orders; `mpfr_ref.py` self-checks them), each note proven, no in-range row noted,
    the pure-sum lemma checked partial sum by partial sum on 786 rows, the data-backed Shomate callers in the port's box
    (every SHOMATE_DATA segment; the flame rows' 324,925 evaluations replayed); the shared tracer `_t_chk` (§16, §19,
    §20) exact at 2⁻¹⁰²² (#29's correction)
  - `cr_verify.py`: the twelve families `hess`, `eff`, `rsum`, `qv`, `ionic`, `haber`, `gibbs`, `kohl`, `lim`, `zoc`,
    `shcp`, `shh` (1,250,069 arguments: the public function, the fast path, the slow path forced) and `nernst`'s
    products at 2⁻¹⁰²² (2,800), on x86_64 ELF, x86_64 PE and aarch64; the probe includes `src/reaction.cyr` and
    thermochem
  - tests: `tests/thermo.tcyr` `cyrius_adr33_sums`, `tests/reaction.tcyr` `cyrius_adr33_gibbs`, `tests/nuclear.tcyr`
    `cyrius_adr33_q_value`, `tests/solution.tcyr` `cyrius_adr33_solution`, `tests/inorganic.tcyr`
    `cyrius_adr33_born_haber`, `tests/kinetics.tcyr` `cyrius_adr33_zero_order`, `tests/thermochem.tcyr`
    `cyrius_adr33_reaction_sums` and `cyrius_adr34_shomate` (the survey's edges, Rust's results for non-finite inputs,
    3,000 seeded inputs against the slow paths, every SHOMATE_DATA segment in the fast path),
    `tests/electrochemistry.tcyr` `cyrius_adr29_2m1022`; no Rust test changes
  - benches: the audit's Cyrius-only benches (B1–B4) move to the new `tests/audit.bcyr` (the user's decision of
    2026-11-11, "Split the bench file"), names unchanged, B3's back at the standard iteration counts; B4 adds eleven
    typical and twelve edge benches there
- **The library-wide overflow audit, B5** ([ADR 0002 #35](docs/adr/0002-port-api-divergences.md), B1's decisions
  unchanged; port plan §6.4 "Decided in the overflow audit's B5"; the behaviour changes are under *Breaking*):
  - **reaction** — `equilibrium_constant`'s argument −ΔG/(R·T) under #20's semantics: R·T at 53 bits with an unbounded
    exponent and the quotient rounded once into binary64 (`_km_keq_arg_ux`) where R·T leaves [2⁻¹⁰²¹, DBL_MAX] (the
    inline fast test; 2⁻¹⁰²² excluded, B4's lesson), then the unchanged `_km_cr_exp`; the ComputationError on that K
  - **solution** — the activity models' log₁₀ γ likewise (`_km_dh_limiting_arg_ux`, `_km_dh_extended_arg_ux`,
    `_km_davies_arg_ux`), behind inline "finite" tests (A·z²; the extended law's numerator and (B·å)·√I), then the
    unchanged `_km_cr_exp10`; the ADR proves the tests and the lemma that exp / 10^x of the binary64 argument is the
    correctly rounded exp / 10^x of the 53-bit unbounded one (no new primitive)
  - generator: `scripts/oracle/src/audit.rs`'s B5 (the port's values independently: the argument on `exact::Ux`, glibc's
    exp and `powf(10, x)` of it, checked against MPFR by `validate.py`; Rust's validation replicated and asserted;
    `noted_k` for the ComputationError flips), `m_reaction.rs` and `m_solution.rs` through it; edge rows appended with
    their own seeds: reaction +676, solution +1,190
  - golden modules `reaction` 2,406 and `solution` 30,644 rows; the oracle's values kept, 845 rows noted with the port's
    result — 326 `equilibrium_constant` (27 the port's ComputationError where Rust's K is finite, 33 its K where Rust's
    overflows), 55 `debye_huckel_limiting`, 385 `debye_huckel_extended`, 79 `davies_activity`; the harness maps and
    counts them per function (`_par_got_noted_k`); none of the 2,627 earlier rows changes
  - `validate.py` §23: every row of the four functions recomputed as #20's target itself (the 53-bit argument on
    `mpfr_ref.UX` and MPFR's exp / exp10 of it rounded once: the new `mpfr_ref.ux_exp_final`, `ux_exp10_final`, with
    self-checks), the lemma checked on every row, each note proven (so every noted glibc value is MPFR's), the port's fast
    tests proven on the corpus, no in-range row noted (B3's strict tracers), the noted rows counted by class
  - `cr_verify.py`: the families `keq`, `dhl`, `dhe`, `dav` (505,982 arguments: the public function, the fast path, the
    slow path forced on every argument) on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/reaction.tcyr` `cyrius_adr35_keq`, `tests/solution.tcyr` `cyrius_adr35_activity` (the survey's examples
    with MPFR's values, the ComputationError flips, the 2⁻¹⁰²² boundaries, Davies' zero bracket, 3,000 seeded inputs
    against the slow paths); no Rust test changes
  - benches: `bench_audit_b5`, four typical and six edge benches, in a new third bench file `tests/audit2.bcyr` (31 s;
    in `tests/audit.bcyr` they took it to 219 s of the 300 s per-file deadline, so it splits as B4's did)
- **The library-wide overflow audit, B6** (the roadmap's B8; [ADR 0002 #36, #37, #38](docs/adr/0002-port-api-divergences.md),
  B1's decisions unchanged, the user's "FFT: both" and "fix all three iterative edges", find_peaks' exact threshold
  settled under the policy; port plan §6.4 "Decided in the overflow audit's B6"; the behaviour changes are under
  *Breaking*):
  - **spectroscopy** — `fft_spectrum` under #20's semantics over the whole transform (#36): hisab 1.4.0's transcription
    (`_km_fft`) with per-operation operand checks (twiddle components 0 or ≥ 2⁻¹⁰⁰, each butterfly's d 0 or ≥ 2⁻⁹²¹) and
    the retained outputs' finiteness as the fast test, a bin whose square sum leaves [2⁻⁹⁰⁰, DBL_MAX] its unbounded √
    (`_km_fft_mag_ux`), and otherwise the transform with an unbounded exponent (`_km_fft_ux`, `_km_fft_spectrum_ux`; 32·n
    heap bytes); a NaN or ±∞ sample keeps 1.2's IEEE spectrum. `find_peaks`' threshold test exact (#37): Rust's rounded
    threshold behind a proven error bound, else `_km_peak_above_exact`
  - **mathx** — `_km_sign_scaled`: the exact sign of a sum of up to eight terms cᵢ·2^kᵢ at any exponents (rounds of
    Grow-Expansion, the deferred tiny terms carried down), for #37's predicate
  - **thermochem** — `adiabatic_flame_temperature`'s bracket and its transcribed bisection compare the signs of f(T₀) and
    f(T_max) exactly (#38)
  - generator: `scripts/oracle/src/audit.rs`'s B6 (`port_fft` on `exact::Ux`, `port_peaks` on the new
    `exact::peak_above` in `Big` integers), `m_thermochem.rs`'s flame replica with the exact tests (asserting the only
    change is Rust's underflowing bracket); edge rows appended with their own seeds: spectroscopy +2,393 `find_peaks` and
    +515 `fft_spectrum`, thermochem +443 flame rows
  - golden modules `spectroscopy` 70,777 and `thermochem` 11,592 rows; the oracle's values kept, 913 rows noted with the
    port's result — 418 `fft_spectrum` (41 of the 213 earlier rows), 222 `find_peaks` (3 of 994), 273
    `adiabatic_flame_temperature` (16 of 1,937); the harness maps and counts them (`_par_map_toks`, `_par_got_noted`)
  - `validate.py` §24: the flame recomputed with the exact tests (and the sign tests checked against Rust's products over
    1,444 special pairs and 20,000 seeded ones), `find_peaks` with `Fraction`s (the fast test's bound and decisions
    checked on every candidate), the FFT on `mpfr_ref.UX` (MPFR's twiddles, each √ rounded once); each note proven, no
    in-range row noted (B3's strict tracers), the noted rows counted by class
  - `cr_verify.py`: the families `flame` (14,010 arguments), `pk3` … `pk16` (380,000 arrays), `fft1` … `fft256` (110,718
    signals; a vector result per argument) and `signsc` (`_km_sign_scaled`, 140,110 sums against `Fraction`), the slow and
    exact paths forced, on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/thermochem.tcyr` `cyrius_adr38_flame_sign`, `tests/spectroscopy.tcyr` `cyrius: ADR 0002 #36` and
    `#37` (the survey's examples, the fast tests' bounds ulp by ulp, n = 4's product below 2⁻¹⁰²², ±DBL_MAX thresholds,
    `_km_sign_scaled` across rounds, 3,000 and 20,000 seeded inputs against the forced slow and exact paths); no Rust
    test changes
  - benches: `bench_audit_b6`, five edge benches, in `tests/audit2.bcyr`
- **The library-wide overflow audit, B7** (the roadmap's B9; [ADR 0002 #39](docs/adr/0002-port-api-divergences.md), B1's
  decisions unchanged, "Balancing: multi-word integers, no size limit" (the survey's remedy (a), settled by the main
  session under the user's policy, not asked); port plan §6.4 "Decided in the overflow audit's B7"; the behaviour change
  is under *Breaking*):
  - **stoichiometry** — `balance_equation` runs Rust's i64 elimination first (`_km_balance_i64`, the old body with each
    overflow check returning a restart code), and at its first overflow the same steps from the original matrix on
    signed multi-word integers (`_km_balance_big`): the pivot search, the swap, reduce_row, the sign flip, the
    fraction-free row update, the lcm, x, reduce_row(x), the zero, sign and u32 checks, then Rust's i128 verification;
    `_km_balance_equation(…, force)` takes the multi-word path from the start (the tests', `cr_verify.py`'s)
  - **stoichiometry** — `_km_sb_*`: signed integers of any size in slots of 31-bit limbs (from a word, copy, compare,
    the signed difference, the product, Knuth's Algorithm D, the exact quotient, Lehmer's gcd with its unimodular
    combination); `_km_balance_bits`: the Hadamard bound B that sizes the tables (every stored entry, the lcm and x below
    2^B, every product and difference below 2^(2B + 1), proven in the ADR)
  - generator: `scripts/oracle/src/exact.rs`'s `Int` (a sign and a `Big`, with Knuth's division in base 2⁶⁴ and Euclid's
    gcd), `audit::port_balance` (Rust's algorithm on `Int`, with a tracer of Rust's checked operations: the port differs
    from Rust only where Rust overflows, asserted); 521 balancing edge rows appended with their own seed
    (`emit_balance_b7`)
  - golden module `stoichiometry` 48,297 rows (3,850 `balance_equation`); the oracle's values kept, 413 rows noted with
    the port's whole result — 104 of the 3,329 earlier rows, exactly the survey's, and 309 edge rows; the harness maps and
    counts them (`_par_map_toks`)
  - `validate.py` §25: every row recomputed with the port's semantics on Python ints and checked against the exact null
    space; each note proven, a strict tracer of Rust's first overflow site (no row where Rust does not overflow is noted;
    207 where it does stay the ComputationError, counted); the size bound and x's gcd of 1 checked on every row
  - `cr_verify.py`: the families `bal23`, `bal34`, `bal56`, `bal89`, `bal1213` and `bal2425` (54,306 matrices of
    2 × 3 … 24 × 25; the public result, the fast path's overflow report against a Python tracer of Rust's checked
    arithmetic, the multi-word path forced), on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/stoichiometry.tcyr` `cyrius_adr39_unbounded_integers`, `cyrius_adr39_size_and_forced`,
    `cyrius_adr39_multiword`; no Rust test changes (two of them still get the ComputationError, now as past u32 rather
    than by the overflow: their comments say so)
  - benches: `bench_audit_b7`, three edge benches, in `tests/audit2.bcyr`
- **The library-wide overflow audit, B10a** (the first part of the roadmap's B10; for
  [ADR 0002 #40](docs/adr/0002-port-api-divergences.md), whose functions change in the later sub-batches: its
  "Addition, 2026-11-12 (audit B10a)"; the user's "new unbounded-exponent primitives … Rust's chain with an unbounded
  exponent, last step rounded once", the semantics settled by the main session under the user's policy, not asked;
  port plan §6.4). No function's result and no golden row changes:
  - **mathx** — `_km_cr_exp_ux`, `_km_cr_expm1_ux` and `_km_cr_exp10_ux`: e^x, e^x − 1 and 10^x of an unbounded-exponent
    value (a `_km_ux_*` pair), returned as a pair with the significand correctly rounded to 53 bits at any exponent —
    no subnormal grid, no overflow at binary64's limits (e^-800 = 1.795…·2⁻¹¹⁵⁵) — exact up to MPFR's widest exponent
    range (|x| to about 3.2·10¹⁸ for e^x, 1.4·10¹⁸ for 10^x; exponents in [−2⁶² + 154, 2⁶² − 155]) and saturating past it
    as MPFR does (+∞; +0 for e^x and 10^x); IEEE's specials; with `_fast` and `_acc` entry points. Three tiers: binary64's
    own range is `_km_cr_exp` / `_km_cr_expm1` / `_km_cr_exp10` itself; up to |x| < 1024 (10^x: 512) the same cores
    without the final scaling; beyond, a new reduction (a double-double one for the fast path while |x·512/ln 2| < 2⁵¹,
    and an integer one, Payne–Hanek's technique on 384 bits of 512/ln 2 or 512·log₂ 10, for the accurate path and past
    it) before the cores' tails (`_km_exp_core_r`, `_km_exp_acc_r`). The fast path is proven (measured 2^-71.07 against
    its 0.9·2⁻⁷⁰ bound); the accurate path is below 2⁻¹⁵⁰ relative (measured 2^-155.15), correct unless a result lies
    within 2⁻¹⁵⁰ of a midpoint (probability about 2⁻⁹⁶ per argument outside binary64's range, as for pow and
    `_km_cr_exp2_scaled`); no new exact or tie case (10²³ stays the only tie)
  - generator: `scripts/oracle/mathx_tables.py` writes `_KM_CR_EXPUX_C` and `_KM_CR_EXP10UX_C` (the reduction's 384-bit
    tables, checked against integer series of ln 2 and ln 10), ln2/512, 512/ln 2 and 512·log₂ 10 as triple-doubles, and
    the saturation edges bisected with MPFR at its widest range (asserting that no argument rounds to MPFR's smallest
    number instead of +0); `scripts/oracle/mpfr_ref.py`'s `UX.exp`, `UX.expm1`, `UX.exp10` and `UX.parts`
  - `cr_verify.py`: the families `exp_ux`, `expm1_ux` and `exp10_ux` (448,188 arguments against MPFR at its widest
    exponent range: specials and exponents to ±2⁶³, every binade of both signs, binary64's range, each tier's edges, near
    multiples of ln2/512 with K to 2⁷¹, Lefèvre–Muller's hard cases, the saturation edges, expm1's −1 region, exp10's
    integers and 10²³; the function, its fast path, the accurate path forced, and the binary64 function's significand
    wherever it is normal), on x86_64 ELF, x86_64 PE and aarch64
  - tests: `tests/mathx.tcyr` `cyrius_b10a_exp_ux`, `cyrius_b10a_expm1_ux`, `cyrius_b10a_exp10_ux` (96 MPFR rows as bit
    patterns, the specials, 4,000 seeded arguments each)
  - benches: `bench_audit_b10a`, fifteen benches beside the binary64 functions, in `tests/audit2.bcyr`
- **M5a numerics, reaction_dynamics** (port plan §4, §6.2–§6.4; [ADR 0003](docs/adr/0003-transcribe-missing-hisab-routines.md)'s
  M5a amendment). The overflow audit is paused after B10a (the user's decision of 2026-11-12: it resumes after M5 / M7
  over the whole library, M5's modules included), so M5a is ported at the plan's parity classes; the roadmap lists its
  candidate spots for the audit:
  - **numerics** — `src/numerics.cyr`, private: hisab 1.4.0's `rk4_trajectory` (`_km_rk4_trajectory`, with `rk4_step`'s
    Neumaier compensation), `dopri45` (`_km_dopri45`: the Dormand–Prince tableau as Rust's const quotients, the step
    factor 0.9·err^(−0.2) on `_km_cr_pow`, 100 000 attempted steps), `gaussian_elimination` (`_km_gaussian_elimination`:
    partial pivoting, the 10⁻¹² pivot threshold) and `euler_maruyama` (`_km_euler_maruyama_system`), operation for
    operation; PCG32 from hisab 3.4.1 (`pcg32_new` / `pcg32_next`, equal to 1.4.0's `Pcg32::new(seed, 1)` bit for bit),
    1.4.0's `next_f64` and Box–Muller `next_normal` (`_km_cr_ln`, `_km_cr_cos`) transcribed on it. Closures are a function
    pointer and a context pointer (`f(ctx, t, y, dy, n)`, through `callptr`); a trajectory is a flat buffer of (t, y[n])
    points, stored while k < cap, with the final state in `fin`; work tables of up to 96 words on the stack
  - **reaction_dynamics** — `src/reaction_dynamics.cyr`: `simulate_first_order`, `simulate_second_order`,
    `simulate_consecutive`, `simulate_reversible`, `simulate_michaelis_menten` (RK4), `simulate_adaptive` (Dormand–Prince),
    `simulate_stiff` (kimiya's backward-Euler Newton loop on the elimination, its 20 iterations and 10⁻¹⁰ / 10⁻⁶
    tolerances) and `simulate_stochastic` (Euler–Maruyama, PCG32 stream 1), with Rust's checks, Err codes and order and
    finish_trajectory's final-state test. The trajectory goes to the caller's buffer `out` of `cap` points of 1 + n f64
    (Err(KM_ERR_INVALID_INPUT) when it doesn't fit, after Rust's own errors); the Jacobian is a flat n × n buffer, zeroed
    once and reused as Rust's; a null function pointer or a bad state slice is Err(KM_ERR_INVALID_INPUT) (port plan §6.4)
  - generator: golden modules `numerics` (243 rows: hisab 1.4.0's routines called directly, the crate's new direct `hisab`
    dependency, the oracle's own locked version) and `reaction_dynamics` (954 rows: every `simulate_*` on the mechanism
    presets at 10, 100 and 1000 steps, the Robertson stiff problem at 1 … 4000 steps, seeded SDE runs on seeds 0 …
    u64::MAX, a linearly converging Newton across the 20-iteration limit and both tolerances, every Err the 1.2 tests pin
    and the steps / time / state edges); the test systems in `scripts/oracle/src/systems.rs` (the README's table);
    `Module::row_summarised` (a trajectory row's note is its point count and last point)
  - corpus: every row bit-identical — RK4, backward Euler and the elimination class B, dopri45, Euler–Maruyama and the
    normal draws class C at bound 0, PCG32's u32 class A; `tests/golden/numerics.cr.tsv` (3 rows) and
    `reaction_dynamics.cr.tsv` (1 row: two glibc pows in a 1,766-point Dormand–Prince run) are the corrected oracle where
    glibc's misrounding changes a trajectory; parity line 339,944 ok, 0 failed, 0 measured, 0 pending, 46,407 accepted
    (ADR 6: 810), 0 glibc-matching on x86_64 ELF, x86_64 PE and aarch64; `numerics` 240 ok, 3 accepted (ADR 6)
  - `validate.py` §26: every row recomputed with its own transcription of hisab 1.4.0 and `reaction_dynamics.rs` (the
    corrected oracle re-runs a trajectory with pow, ln and cos correctly rounded); `simulate_stiff`'s six Newton exits
    each covered, the 100 000-step limit reached, each typed stiff constant pinned (2–30 rows move when it is perturbed),
    MAX_STEPS pinned at 10⁷ / 10⁷ + 1
  - tests: `tests/numerics.tcyr` (the `numerics` rows through `tests/parity_core.cyr`, `tests/parity.tcyr` skipping the
    module, then 5 Cyrius-only groups: heap work tables, capacity, scratch reuse, the pivot edge, the PCG32 handle);
    `tests/reaction_dynamics.tcyr` (the 26 Rust tests, `stiff_resized_jacobian_is_error_not_panic` adapted, and the 2
    integration tests, plus 5 Cyrius-only groups; 74 assertions); `tests/rd_systems.cyr` and `tests/rd_rows.cyr`, the
    systems as callbacks and the trajectory rows
  - benches: `tests/numerics.bcyr` (a new bench file)
  - `[lib]` gains `src/numerics.cyr` and `src/reaction_dynamics.cyr` (after spectroscopy); `dist/kimiya.cyr` regenerated;
    surface parity 410/410 (`reaction_dynamics::type TrajectoryPoint` mapped to the flat point)
- **M5b network, fitting, sensitivity** (port plan §4, §6.2–§6.4; [ADR 0003](docs/adr/0003-transcribe-missing-hisab-routines.md)'s
  M5b amendment, [ADR 0002](docs/adr/0002-port-api-divergences.md)'s M5b addition). Ported at the plan's parity classes
  (the overflow audit stays paused); the roadmap lists the candidate spots for the audit, and the oracle defect below:
  - **numerics** — `src/numerics.cyr` gains hisab 1.4.0's `least_squares_poly` (`_km_least_squares_poly`: the Vandermonde
    by `powi`, `qr_decompose`'s modified Gram–Schmidt with its 10⁻¹² threshold, back substitution), `matrix_rank` on
    `svd`'s one-sided Jacobi (`_km_matrix_rank`, the singular values only) and `CsrMatrix::from_dense`
    (`_km_csr_from_dense`, |v| > 10⁻¹² into a hisab 3.4.1 `csr_new` handle), operation for operation. hisab 3.4.1's
    `svd_compute_checked` disagrees with 1.4.0 on 191 of the 563 golden rank rows and its `csr_from_dense` keeps the
    entries 1.4.0 drops, so both are transcribed, as ADR 0003 provides
  - **network** — `src/network.cyr`: `stoichiometric_matrix`, `sparse_stoichiometric_matrix`, `reaction_network_rank`,
    `conservation_law_count`, `rate_from_stoichiometry`, with Rust's checks, Err codes and order. A matrix is a flat
    row-major buffer and its dimensions; a reaction a 32-byte record (the reactants' and the products' term lists, each
    term a 16-byte species / coefficient record); the dense results go to the caller's buffer (Err(KM_ERR_INVALID_INPUT)
    when they don't fit, after Rust's errors); the CSR form is hisab 3.4.1's handle; the counts are usize words
  - **fitting** — `src/fitting.cyr`: `fit_arrhenius` (into a `struct KmArrheniusFit`, through a pointer), `fit_polynomial`
    (the coefficients to the caller's buffer), `fit_beer_lambert`, and `struct KmPolynomialFit` (a view of a coefficient
    buffer) with `evaluate` and `degree`; the ln per point and the exp of the intercept correctly rounded
  - **sensitivity** — `src/sensitivity.cyr`: the four gradients (ADR 0002 #2, below) as closed forms — each fixed tape's
    reverse sweep unrolled, with hisab 1.4.0's zero-adjoint skip and the correctly rounded exp and ln — and
    `compute_gradient(f, ctx, tape, parameters, n, out)` on hisab 3.4.1's reverse-mode tape (`ad_tape_new`, `ad_*`,
    `ad_grad`), the callback f(ctx, tape, vars, n) returning the output node through `callptr`
  - generator: golden modules `network` (1,296 rows: the mechanism presets — the Rust tests', the rate laws,
    Michaelis–Menten, Chapman, the H₂ + Br₂ chain, Robertson, Lotka–Volterra, the Brusselator, the Oregonator, H₂/O₂ —
    their transposes and scalings, rank-deficient, noisy, Hilbert, empty and non-finite matrices, the CSR and rank
    thresholds to the ulp, every Err), `fitting` (578: exact synthetic Arrhenius data from k = 2^e, seeded noise, the
    QR's threshold swept, exact integer polynomials of degree 0 … 7, ill-conditioned and degenerate Vandermonde systems,
    every Err) and `sensitivity` (894: the four gradients over each argument's edges and exp's ends, `compute_gradient`
    on 20 test functions built on hisab 1.4.0's `Tape` directly, the crate's `hisab` dependency gaining `autodiff`,
    kimiya's own feature); the test functions in `scripts/oracle/src/m_sensitivity.rs` (the README's table)
  - corpus: every row 0 ulp but `compute_gradient`'s (≤ 40 ulp, class D: hisab 3.4.1's tape against 1.4.0's) — the
    builders and the rank class A, the rates, fits and `evaluate` class B, `fit_arrhenius` class C at bound 0, the four
    gradients class D at 0 ulp; no ADR 0006 correction (glibc 2.44 rounds every call of these rows correctly); parity line
    342,485 ok, 0 failed, 0 measured, 0 pending, 46,634 accepted (ADR 6: 810), 0 glibc-matching on x86_64 ELF, x86_64 PE
    and aarch64 (22 modules)
  - `validate.py` §27: every row recomputed with its own transcription of hisab 1.4.0 (`from_dense`, the Jacobi SVD,
    the QR least squares, the reverse-mode tape) and the three Rust modules; ADR 0002 #2's notes checked against the
    validation rule; the rank's and the CSR's thresholds (one ulp down) and the QR's (halved, doubled) each pinned by
    golden rows
  - tests: `tests/network.tcyr`, `tests/fitting.tcyr`, `tests/sensitivity.tcyr` (the 34 Rust tests, two adapted, and
    fitting's 2 integration tests, plus 11 Cyrius-only groups; 155 assertions); `tests/sens_functions.cyr`, the test
    functions as callbacks
  - benches: 7 in `tests/numerics.bcyr` — `fitting/arrhenius_6pts` 233.40 ns → 1.01 µs, `fitting/polynomial_deg2_20pts`
    500.05 ns → 2.32 µs (Rust's benches), and 5 Cyrius-only ([benchmarks](docs/benchmarks-rust-v-cyrius.md))
  - `[lib]` gains `src/network.cyr`, `src/fitting.cyr` and `src/sensitivity.cyr` (after reaction_dynamics);
    `dist/kimiya.cyr` regenerated; surface parity 431/431 (`network::type ReactionDef` mapped to the record)

### Breaking (relative to Rust 1.2.0, ADR 0002 #1)
- **thermochem** — `reaction_enthalpy` and `reaction_gibbs_energy` return **J/mol**, not kJ/mol: the value is 1.2's
  kJ/mol result × 1000 exactly (Rust's difference, then one multiply). `reaction_entropy` (J/(mol·K)) and the table
  fields `thermochem_data_delta_hf_kj` / `_delta_gf_kj` (kJ/mol, as tabulated) are unchanged.
  *Migration:* drop the `* 1000.0` before `vant_hoff_k` / `adiabatic_flame_temperature` —
  `adiabatic_flame_temperature(dh_kj * 1000.0, …)` becomes `adiabatic_flame_temperature(dh_j, products, n, t0)` with
  `var t, dh_j = reaction_enthalpy(p, np, r, nr);`; code that displays or compares the result in kJ/mol divides by
  1000; in ΔG = ΔH − TΔS with `reaction_entropy`, the `/ 1000.0` on TΔS goes away.

### Breaking (relative to Rust 1.2.0, ADR 0002 #2)
- **sensitivity** — `equilibrium_constant_gradient`, `arrhenius_gradient`, `nernst_gradient` and
  `gibbs_temperature_sensitivity` return a `Result` instead of a bare vector or value: the gradient goes to the caller's
  buffer (two or three f64), `gibbs_temperature_sensitivity` is `Ok(−ΔS)`. They validate as `equilibrium_constant`,
  `arrhenius_rate`, `nernst_potential` and `gibbs_free_energy` do — `Err(KM_ERR_INVALID_TEMPERATURE)` for T not positive
  and finite, `Err(KM_ERR_INVALID_INPUT)` for a non-finite energy, A, E° or ΔS, n outside 1 ..= 255 or Q not positive
  and finite — and return `Err(KM_ERR_COMPUTATION_ERROR)` where an entry is not finite, where 1.2 returned NaN, ±∞ or
  a finite but meaningless entry (T ≤ 0). Every valid input gives 1.2's values (0 ulp on the corpus).
  `compute_gradient` takes the closure as a function pointer and a context, the caller's hisab tape and the output
  buffer, and its `Result` is `Err` only for a null pointer, a negative count or a tape too small (1.2 never failed).
  *Migration:* `let g = arrhenius_gradient(a, ea, t);` becomes `var gr: f64[3]; var t, v = arrhenius_gradient(a, ea, t,
  &gr);` and a check of `is_ok(t)`; code that tested `g[i].is_finite()` tests the `Result`; a `compute_gradient`
  closure becomes `fn f(ctx, tape, vars, n)` building on hisab 3.4.1's `ad_*` operations, called as
  `compute_gradient(&f, ctx, ad_tape_new(cap), params, n, out)` (reuse one tape across calls).

### Breaking (relative to Rust 1.2.0, ADR 0002 #4)
- **phase** — `antoine_pressure_mmhg` / `antoine_pressure_pa` return a `Result` instead of a bare value, and
  return `Err` where 1.2 extrapolated or returned NaN: `Err(KM_ERR_INVALID_TEMPERATURE)` for T outside the
  coefficient set's `t_min_c ..= t_max_c` (NaN included) or at/below the pole (C + T ≤ 0),
  `Err(KM_ERR_INVALID_INPUT)` for a non-finite A, B, C or range end. There is no extrapolating variant.
  *Migration:* `let p = antoine_pressure_mmhg(c, t);` becomes
  `var tag, v = antoine_pressure_mmhg(c, t); if (is_ok(tag) == 1) { var p: f64 = v; … }`; code that relied on
  extrapolation must use a coefficient set fitted over its temperatures (or `clausius_clapeyron`); code that
  tested `p.is_nan()` tests `is_err_result(tag)`. In range, the value is 1.2's to within the class-C bound.

### Breaking (relative to Rust 1.2.0, ADR 0002 #6)
- **organic** — `Geometry::TShapedThreeLonePairs` (three bonding groups and **two** lone pairs, ClF₃: the Rust name is a
  misnomer) is `KM_GEOMETRY_T_SHAPED`, at its Rust position (value 7); `FunctionalGroup::hydrogen_bonding` is
  `functional_group_is_hydrogen_bond_donor`, with the same result for every group (thiols are donors, ethers, aldehydes,
  ketones, esters and nitriles are not). New beside it: `functional_group_is_hydrogen_bond_acceptor`, which 1.2 had no
  counterpart of: by IUPAC 2011's criterion every group but the alkane accepts, weak π, sulfur and halogen acceptors
  included.
  *Migration:* `Geometry::TShapedThreeLonePairs` → `KM_GEOMETRY_T_SHAPED` (`predict_geometry(3, 2)` returns it);
  `g.hydrogen_bonding()` → `functional_group_is_hydrogen_bond_donor(g)`; code that treated "not `hydrogen_bonding`" as
  "doesn't hydrogen-bond" asks `functional_group_is_hydrogen_bond_acceptor(g)` too (an ether accepts, an alkane
  neither donates nor accepts). Code that wants only the strong O / N acceptors can't use the acceptor flag alone: the
  π systems, the thiol and the halide answer 1 as weak acceptors.

### Breaking (relative to Rust 1.2.0, ADR 0002 #12)
- **inorganic** — `cfse_octahedral` returns the exact CFSE, (3·n(e_g) − 2·n(t₂g))/5 correctly rounded (the integer
  numerator, one division), not 1.2's `t2g * (-0.4) + eg * 0.6`. They differ on 13 of the 22 Ok cases: high-spin d⁵
  and d¹⁰ (both spins) are now exactly 0 (1.2: −2.2 × 10⁻¹⁶ and −4.4 × 10⁻¹⁶), and the 10 others by 1–4 ulp (d³: −1.2
  where 1.2 gave −1.2000000000000002; low-spin d⁶: −2.4, not −2.4000000000000004). Errors are unchanged.
  *Migration:* nothing for code that compares with a tolerance, as the 1.2 tests do; code that tested 1.2's exact bits
  (or worked around them, e.g. `if cfse.abs() < 1e-15 { 0.0 }`) can compare with the exact values directly
  (`cfse == 0.0` for d⁵ high-spin and d¹⁰, `cfse == -1.2` for d³).

### Breaking (relative to Rust 1.2.0, ADR 0002 #3)
- **nuclear** — `semi_empirical_binding_energy` returns a `Result` instead of a bare value: `Err(KM_ERR_INVALID_INPUT)`
  where 1.2 returned NaN — Z > A or A = 0, no such nucleus — and for a Z or A outside Rust's u32. Every finite value is
  1.2's, bit for bit.
  *Migration:* `let be = semi_empirical_binding_energy(z, a);` becomes
  `var t, v = semi_empirical_binding_energy(z, a); if (is_ok(t) == 1) { var be: f64 = v; … }`; code that tested
  `be.is_nan()` tests `is_err_result(t)` (the payload is `KM_ERR_INVALID_INPUT`).

### Breaking (relative to Rust 1.2.0, ADR 0002 #13)
- **nuclear** — `atoms_after_half_lives(n0, n)` returns the correctly rounded N₀/2ⁿ (rounded once) where Rust's
  intermediate 2⁻ⁿ is not a normal double (n > 1022 or n ≤ −1024), instead of 1.2's `n0 * 2.0.powf(-n)`, which there
  overflows to ±∞, underflows to 0 or loses bits as a subnormal: 1000 after 1075 half-lives is 500·2⁻¹⁰⁷⁴ ≈ 2.47 × 10⁻³²¹
  (1.2: 0), 2¹⁰⁰ after 1,100 is 2⁻¹⁰⁰⁰ (1.2: 0), the smallest subnormal after −1080 is 64 (1.2: +∞), ties go to even
  (2⁻¹⁰⁷⁵ → 0, 1.5·2⁻¹⁰⁷⁴ → 2⁻¹⁰⁷³), and a fractional n too (DBL_MAX after 1022.25 is 3.363585661014858; M4c's port was
  faithful there, one ulp off on a third of the corpus's fractional rows, where 1.2 could be ±∞, 0 or up to 2⁵² ulp off).
  ±0, ±∞ and NaN N₀ give N₀ there (0 after −1100 half-lives is 0; 1.2: NaN). Everywhere else — every n with 2⁻ⁿ normal,
  an infinite n — the value is 1.2's.
  *Migration:* nothing for code within the normal range (|n| ≲ 1022 half-lives); code that relied on 1.2's 0, ±∞ or NaN
  past it gets the representable N₀/2ⁿ (`r == 0.0` still holds wherever |N₀/2ⁿ| is at most 2⁻¹⁰⁷⁵).

### Breaking (relative to Rust 1.2.0, ADR 0002 #14)
- **nuclear** — `is_secular_equilibrium(λ_p, λ_d)` tests λ_p < λ_d/100 exactly, not 1.2's `λ_p < 0.01 * λ_d`, whose 0.01
  is 0.0100000000000000002 and whose product rounds or underflows. They differ only within an ulp or so of the boundary,
  and where 0.01·λ_d underflows: (0.03, 3) is now true (0.03 is 0.0299999…), (−0.02, −2) true, (0, 5 × 10⁻³²⁴) true
  (1.2: false for all three). NaN still gives false; ±∞ and ±0 compare as in 1.2.
  *Migration:* none in practice — the threshold is a heuristic for "≪"; code that tested the exact 1.2 boundary uses the
  exact one.

### Breaking (relative to Rust 1.2.0, ADR 0002 #15)
- **inorganic** — `delta_tet_from_oct(Δ)` is the correctly rounded 4Δ/9 for every finite Δ: where 1.2's Δ·4 overflows
  (|Δ| ≥ 2¹⁰²², above DBL_MAX/4) it returned ±∞; the port returns (Δ/9)·4 there (`delta_tet_from_oct(DBL_MAX)` ≈
  7.99 × 10³⁰⁷). Every other value is 1.2's, bit for bit.
  *Migration:* none for physical Δ_oct values; code that relied on ±∞ above 4.5 × 10³⁰⁷ gets the finite value.

### Breaking (relative to Rust 1.2.0, ADR 0002 #16)
- **potential** — `lj_equilibrium_distance(σ)` returns the correctly rounded σ·2^(1/6), not 1.2's `sigma * 2.0.powf(1.0 /
  6.0)`, which rounds twice (2^RN(1/6) is RN(2^(1/6)) = 1.122462048309373, then the product) and is one ulp off on about a
  fifth of all σ (4,117 of 20,000 seeded; Ar's 3.4 × 10⁻¹⁰ m is not among them). NaN, ±0 and ±∞ are unchanged.
  *Migration:* none for code that compares with a tolerance; code that tested 1.2's exact bits sees at most one ulp.

### Breaking (relative to Rust 1.2.0, ADR 0002 #17)
- **potential** — `lennard_jones(ε, σ, r)` for ε ≥ 2¹⁰²² J returns (ε·d)·4, the correctly rounded 4·ε·d for the same
  d = (σ/r)¹² − (σ/r)⁶, where 1.2's (4·ε)·d overflowed: ±∞, or NaN at r = σ (`lennard_jones(DBL_MAX, 1, 1)` is 0, 1.2
  NaN). Every other value is 1.2's, bit for bit. (Folded into #32 by the overflow audit's B3, with the same values.)
  *Migration:* none for physical well depths.

### Breaking (relative to Rust 1.2.0, ADR 0002 #18)
- **potential** — `morse_width_parameter(ω, μ, D_e)` for D_e ≥ 2¹⁰²³ J computes μ/(2·D_e) as (μ·½)/D_e, the correctly rounded
  value, where 1.2's 2·D_e overflowed and gave a = 0 (`morse_width_parameter(4, 1, 2¹⁰²³)` is 2⁻⁵¹⁰). Every other value is
  1.2's, bit for bit. (Folded into #32 by the overflow audit's B3, which rounds the whole Morse width once: see there.)
  *Migration:* none for physical well depths.

### Breaking (relative to Rust 1.2.0, ADR 0002 #19)
- **spectroscopy** — `wavenumber_from_wavelength(λ)` and `wavelength_from_wavenumber(ν̃)` return the correctly rounded
  1/(100·x), not 1.2's `0.01 / x`, whose 0.01 is 1/100 + 2.08 × 10⁻¹⁹ and whose quotient rounds again: they differ by one
  ulp on about 14 % of inputs (`wavenumber_from_wavelength(17)` is 1/1700 = 0.000588235294117647, 1.2 the next double
  up). The 1.2 tests' values (500 nm ↔ 20 000 cm⁻¹, 10 µm ↔ 1000 cm⁻¹) are unchanged, and so are the errors.
  *Migration:* none for code that compares with a tolerance; code that tested 1.2's exact bits sees at most one ulp.

### Breaking (relative to Rust 1.2.0, ADR 0002 #20)
- **potential** — `distance_3d`, `bond_angle` and `dihedral_angle` no longer over- or underflow in their squares, cross
  products and dot products: each is Rust's own operation sequence evaluated with an unbounded exponent range (every
  + − × ÷ √ rounded to 53 bits, the comparisons and Err checks on those values), the final √ or atan2 rounded once into
  binary64. Where 1.2's chain stays in the normal range (every physical geometry: coordinates and bonds between about
  10⁻¹⁵⁰ and 10¹⁵⁰, 10⁻⁷⁵ and 10⁷⁵ for dihedrals, or zero) the results are 1.2's, bit for bit; elsewhere 1.2 returned an
  Err, 0, ±∞ or a degraded value.
  Examples: `distance_3d([DBL_MAX, 0, 0], [0, 0, 0])` is DBL_MAX (1.2: +∞) and `distance_3d([5e-324, 0, 0], [0; 3])` is
  5e-324 (1.2: 0); `bond_angle([1, 0, 0], [0; 3], [1, 2⁻¹⁰⁶⁰, 0])` is 2⁻¹⁰⁶⁰ (1.2: 0); the gauche dihedral scaled by 2⁷⁰⁰ is
  π/2, as at 1 Å (1.2: Err). The angles are now exactly scale-invariant (2^k times every coordinate gives the same bits
  whenever the scaled coordinates are exact), as Rust documents them, and the coincidence and collinearity tolerances
  hold at every scale; non-finite coordinates still give 1.2's Err, NaN or ±∞.
  *Migration:* none for physical geometries; code that relied on an Err or ±∞ for extreme coordinates gets the
  geometry's value.

### Breaking (relative to Rust 1.2.0, ADR 0002 #21)
- **spectroscopy** — `bohr_energy_level`, `bohr_transition_energy`, `rydberg_wavelength` and
  `spectral_series_wavelength` return their exact values rounded once: −R·Z²/n², R·Z²·(n₂² − n₁²)/(n₁²·n₂²) and
  n₁²·n₂²/(R·Z²·(n₂² − n₁²)), with R the constant's double (RYDBERG_EV, RYDBERG, or 1.2's const-evaluated R_H). 1.2 rounded
  three to eight times and cancelled for close large n: `bohr_energy_level(5, 5)` is −R exactly (1.2: one ulp off, as for
  86 of the 255 Z = n), ordinary lines move by 0–2 ulp, and close large n by much more (`rydberg_wavelength(255, 2²⁸ + 1,
  2²⁸ + 2)` by 9.0·10⁷ ulp, `bohr_transition_energy(26, 2³² − 10, 2³² − 9)` by 1.8·10⁹ ulp, where 1.2's two levels cancel).
  The errors and their order are 1.2's.
  *Migration:* none for code that compares with a tolerance at ordinary n; at n beyond about 10³ the new values are the
  accurate ones.

### Breaking (relative to Rust 1.2.0, ADR 0002 #22)
- **spectroscopy** — `reduced_mass(m₁, m₂)` returns the correctly rounded m₁·m₂/(m₁ + m₂) for every pair of positive finite
  masses: m/2 exactly for a homonuclear diatomic (1.2 was one ulp off on 8.6 % of masses), and no over- or underflow at the
  range ends (`reduced_mass(DBL_MAX, DBL_MAX)` is DBL_MAX/2, 1.2's NaN; `reduced_mass(1e-160, 1e-160)` is 5e-161, where
  1.2's product underflowed to a few subnormal bits). Subnormal results are rounded once, ties to even; the errors are
  1.2's.
  *Migration:* none for physical masses compared with a tolerance; code that tested 1.2's exact bits sees at most one ulp
  in the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #7)
- **stoichiometry** — `assign_oxidation_states` returns each state as an exact rational, a (Z, numerator, denominator)
  record with the fraction reduced and the denominator positive: where 1.2 returned `Err(ComputationError)` "non-integer
  average oxidation state", 2.0 returns the average (C₃H₈'s C −8/3, Fe₃O₄'s Fe +8/3, KO₂'s O −1/2, O₂⁻ −1/2); an
  integer state n is n/1. An average outside i32 stays `Err(ComputationError)`, as in 1.2. `balance_equation` takes an
  integer matrix (flat row-major, m rows by n columns) instead of `&[Vec<f64>]`; Rust's bound |entry| ≤ 2⁵³ is kept.
  Every other error, its kind and its order are 1.2's.
  *Migration:* read the denominator (1 for every state 1.2 returned); pass the composition counts as integers (1.2's
  near-integer f64 entries are their rounded integers; there is no non-integer entry to pass, and no ragged row).

### Breaking (relative to Rust 1.2.0, ADR 0002 #24)
- **stoichiometry** — `percent_yield`, `theoretical_yield_moles` and `percent_composition` return the correctly rounded
  100·a/t, m·p/c (at #25's limiting reagent) and mass·count·100/mw (mw `molecular_weight`'s own value), rounded once,
  subnormal results and ties included, where 1.2 rounded two to four times: about a quarter of results move by 1–2 ulp
  (`percent_yield(1, 3)` is 33.333333333333336, 1.2's 33.33333333333333; water's O 88.80932556203165 %, 1.2's …163), and
  the range ends are right (`percent_yield(5e-324, 10)` is 5e-323, 1.2's 0; `theoretical_yield_moles([1e300], [1e-10],
  1e-10)` is 1e300, 1.2's +∞). A zero actual yield or limiting amount, or a count of 0, gives 1.2's ±0; the errors and
  their order are 1.2's.
  *Migration:* none for code that compares with a tolerance.

### Breaking (relative to Rust 1.2.0, ADR 0002 #25)
- **stoichiometry** — `limiting_reagent` compares the ratios moles/coefficient exactly (mᵢ·c_j against m_j·cᵢ, the first
  index on an exact tie), so `theoretical_yield_moles` follows: where two rounded ratios tie within an ulp, or both
  overflow to +∞ or underflow to 0, 1.2 could name the wrong reactant (`[DBL_MAX, DBL_MAX]` over `[5e-324, 1e-300]` is 1,
  1.2's 0). The errors and their order are 1.2's.
  *Migration:* none for reactant amounts that differ by more than an ulp in ratio.

### Breaking (relative to Rust 1.2.0, ADR 0002 #26)
- **stoichiometry** — `empirical_formula` evaluates its mole ratios and the tolerance |r·k − round(r·k)| ≤ 1/10 exactly.
  Inputs at the ±1/10 edge change between Ok and Err or change their multiplier (C 12.011 % with O 17.5989 % is CO, 1.2's
  Err; C 12.011 % with O 20.7987 % is Err, 1.2's C₃O₄; C 50 % with C 59 % is 5 : 6, 1.2's 6 : 7), and the Errs that came
  from rounding are gone: percentages whose moles underflow are evaluated (He 2⁻¹⁰⁷⁴ % with He 2·2⁻¹⁰⁷⁴ % is 1 : 2, 1.2's
  Err), an overflowing ratio is still the ComputationError (every such ratio is past 2³²). Every error is still
  `ComputationError` or 1.2's validation error, in 1.2's order.
  *Migration:* none for measured compositions away from the knife edge.

### Breaking (relative to Rust 1.2.0, ADR 0002 #27)
- **thermo**, **phase** — `heat_transfer`, `specific_heat_from_calorimetry`, `boiling_point_elevation`,
  `freezing_point_depression` and `osmotic_pressure` return the exact product or quotient rounded once — RN(m·c·ΔT),
  RN(q/(m·ΔT)), RN(K·m·i), RN(i·M·R·1000·T) with R as `KM_GAS_CONSTANT`'s double — where 1.2 rounded two to four times:
  about a quarter of ordinary results move by 1–2 ulp (`heat_transfer(100, 4.184, −300)` is −125520, 1.2's
  −125520.00000000001; `boiling_point_elevation(1.22, 2.5, 1.9)` 5.795, 1.2's 5.794999999999999), and the range ends are
  right (`specific_heat_from_calorimetry(4184, 1000, DBL_MAX)` is 2.33·10⁻³⁰⁸, 1.2's 0; `osmotic_pressure(5e-324, 0.1,
  298.15)` 1.22·10⁻³¹⁸, 1.2's 0; `heat_transfer(1e200, 1e200, 1e-200)` 1e200, 1.2's +∞; `heat_transfer(1e300, 1e300, 0)`
  +0, 1.2's NaN). A zero factor gives the product's signed zero; a non-finite input of the three unvalidated functions
  gives 1.2's IEEE result; the errors and their order are 1.2's.
  *Migration:* none for code that compares with a tolerance; code that tested 1.2's NaN for a zero ΔT, K or i against an
  overflowing product gets 0.

### Breaking (relative to Rust 1.2.0, ADR 0002 #28)
- **kinetics** — `zero_order_half_life`, `steady_state_intermediate` and `pre_equilibrium_rate` return RN([A]₀/(2k)),
  RN(k₁·[A]/k₂) and RN(k₂·k₁·[A]/k₋₁), rounded once (`zero_order_half_life(DBL_MAX, DBL_MAX)` is ½, 1.2's 0;
  `pre_equilibrium_rate(1e300, 1, 1e300, 0)` 0, 1.2's NaN; `pre_equilibrium_rate(10, 2⁻¹⁰²², 0.1, 1)` 2¹⁰²², 1.2's +∞; about a
  quarter of ordinary `steady_state_intermediate` and `pre_equilibrium_rate` results move by 1–2 ulp;
  `zero_order_half_life` changes only for k ≥ 2¹⁰²³). `second_order_concentration` and `michaelis_menten` evaluate 1.2's
  operation sequence with an unbounded exponent, the final division rounded once: every input whose 1.2 chain stays in
  the normal range keeps its bits, and the edges are right (`second_order_concentration(1, DBL_MAX, 10)` 5.56·10⁻³¹⁰, 1.2's
  0; `michaelis_menten(100, DBL_MAX, DBL_MAX)` 50, 1.2's NaN; `michaelis_menten(DBL_MAX, 5, 10)` 1.198·10³⁰⁸, 1.2's +∞;
  `second_order_concentration(DBL_MAX, 0, 0)` DBL_MAX's predecessor, 1/RN₅₃(1/DBL_MAX), 1.2's +∞). The errors and their
  order are 1.2's.
  *Migration:* none for code that compares with a tolerance or stays within the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #29)
- **solution**, **electrochemistry** — `dilution`, `solution_conductivity`, `faraday_mass_deposited`,
  `faraday_charge_required` and `faraday_mass_ratio` return the exact quotient rounded once — RN(M₁·V₁/V₂),
  RN(Λ·c/1000), RN(Q·M/(n·F)), RN(m·n·F/M), RN(m_A·M_B·n_A/(n_B·M_A)) with F as `KM_FARADAY`'s double — where 1.2
  rounded two to four times: about a quarter of ordinary results move by 1–2 ulp (`dilution(1.5, 0.2, 0.6)` is 0.5,
  1.2's 0.5000000000000001; `solution_conductivity(126.4, 0.01)` 0.0012640000000000001, 1.2's 0.001264), and the range
  ends are right (`faraday_mass_deposited(DBL_MAX, 63.546, 2)` is 5.92·10³⁰⁴, 1.2's +∞;
  `solution_conductivity(126.4, DBL_MAX)` 2.27·10³⁰⁷, 1.2's +∞; `faraday_mass_ratio(1, 5e-324, 3, 1e-300, 1)` 6.07·10²³,
  1.2's +∞). A zero factor gives the quotient's signed zero; a non-finite input of `solution_conductivity` gives 1.2's
  IEEE result; the errors and their order are 1.2's. `nernst_potential` evaluates 1.2's operation sequence with an
  unbounded exponent, the final subtraction rounded once: every input whose 1.2 chain stays in the normal range keeps
  its bits (every `nernst_potential_25c` input does), and the edges are right
  (`nernst_potential(0.342, 2, DBL_MAX, 0.1)` 1.78·10³⁰⁴, 1.2's +∞; at Q = 1 E°, 1.2's NaN; T = 10⁻³¹⁰ K's subnormal
  product no longer degraded).
  *Migration:* none for code that compares with a tolerance or stays within the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #30)
- **spectroscopy** — `absorbance`, `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant` and
  `diatomic_moment_of_inertia` return the exact product or quotient rounded once — RN(ε·l·c), RN(A/(ε·l)),
  RN(h·B·J·(J + 1)), RN(h/(8·π·π·I)) with h as `KM_PLANCK`'s double and π its double, RN(μ·r·r) — and
  `vibrational_frequency` the correctly rounded √ of the exact k/μ, where 1.2 rounded two to four times: a fraction of
  ordinary results move by 1 ulp (HCl's `vibrational_frequency(516, 1.6266e-27)` is 563228298138880, 1.2's
  563228298138879.9; `absorbance(18716, 0.1, 0.0007)` 1.31012, 1.2's 1.3101200000000002), and the range ends are right
  (`vibrational_frequency(500, 2⁻¹⁰²²)` 1.50·10¹⁵⁵, 1.2's +∞; `vibrational_frequency(1e-300, 1e30)` 10⁻¹⁶⁵, 1.2's 0;
  `rigid_rotor_energy(2³² − 1, 1e-300)` 1.22·10⁻³¹⁴, 1.2's 0; `rotational_constant(1e-320)` 8.392109·10²⁸⁴, 1.2's
  8.392090·10²⁸⁴; `concentration_from_absorbance(1e300, 1e200, 1e200)` 10⁻¹⁰⁰, 1.2's 0; `absorbance(1e300, 1e300, 0)`
  +0, 1.2's NaN). A zero factor gives the product's signed zero; a non-finite unvalidated input (absorbance's, the
  moment of inertia's, the rotor's B, `concentration_from_absorbance`'s A) gives 1.2's IEEE result; the errors and their
  order are 1.2's.
  *Migration:* none for code that compares with a tolerance; code that tested 1.2's NaN for a zero concentration against
  an overflowing ε·l gets 0.

### Breaking (relative to Rust 1.2.0, ADR 0002 #31)
- **gas** — `ideal_gas_pressure`, `ideal_gas_volume` and `gas_density` return the exact quotient rounded once —
  RN(n·R·T/V), RN(n·R·T/P), RN(P·M/(R·T)) with R as `KM_GAS_CONSTANT`'s double — where 1.2 rounded three times: about a
  third of ordinary results move by one ulp (air at STP, `gas_density(101325, 0.029, 273.15)`, is 1.2938359687824854,
  1.2's …852), and the range ends are right (`ideal_gas_volume(DBL_MAX, 300, 101325)` 4.43·10³⁰⁶, 1.2's +∞;
  `gas_density(101325, 0.029, DBL_MAX)` 1.97·10⁻³⁰⁶, 1.2's 0). `van_der_waals_pressure`, `peng_robinson_pressure` and
  `redlich_kwong_pressure` evaluate 1.2's operation sequence with an unbounded exponent, the final subtraction rounded
  once: every input whose 1.2 chain stays in the normal range keeps its bits (every physical state), and the edges are
  right (`peng_robinson_pressure` with Tc = 2⁻¹⁰²² 99,773.55 Pa, 1.2's NaN; at V = 2·10¹⁵⁴ with Tc = 10¹⁵⁰ −1.49·10⁻⁷,
  1.2's +1.25·10⁻¹⁵¹; `redlich_kwong_pressure(1e130, 1e300, 300, 1)` −1.71·10²⁵, 1.2's −∞). The co-volume tests V ≤ n·b
  and V ≤ b are decided on the unbounded co-volume, so a state whose 1.2 co-volume rounded at the subnormal grid (or a
  subnormal Tc's rounding scaled back up by a tiny Pc) can be the port's Ok where 1.2 refused, or the reverse; every
  other error and the order of the checks are 1.2's.
  *Migration:* none for code that compares with a tolerance or stays within the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #32)
- **potential** — `coulomb` and `coulomb_force` return the exact quotient rounded once (RN(k·q₁·q₂/r), RN(k·q₁·q₂/(r·r)),
  k as `KM_COULOMB_CONSTANT`'s double), where 1.2 rounded three or four times: a fraction of ordinary results move by one
  ulp, and the range ends are right (`coulomb_force(5e-324, 1.602176634e-19, 1e-10)` 7.11·10⁻³¹³, 1.2's 0;
  `coulomb_force(e, −e, 1e-160)` −2.31·10²⁹², 1.2's −∞); a non-finite charge gives 1.2's IEEE result.
  `morse_width_parameter` returns the correctly rounded ω·√(μ/(2·D_e)), where 1.2 rounded the quotient, its √ and the
  product (a third of ordinary results move by one ulp; `morse_width_parameter(5.634e14, DBL_MAX, 7.399e-19)` is
  6.21·10¹⁷⁷, 1.2's +∞; `(5.634e14, 1.6267e-27, 1e300)` 1.61·10⁻¹⁴⁹, 1.2's 0); since #18 is folded in, its D_e ≥ 2¹⁰²³
  values are now correctly rounded too (`morse_width_parameter(1, 2⁻¹⁰⁷⁴, 2¹⁰²³)` is 2⁻¹⁰⁴⁹ where #18 gave 0, and 53 of
  #18's 181 golden values move, 42 by one ulp). `lennard_jones` evaluates 1.2's sequence with an unbounded exponent and
  the final product rounded once (#17's (ε·d)·4 folded in, the same values): every input whose 1.2 chain stays in the
  normal range keeps its bits, and the edges are right (`lennard_jones(1e-300, 1, 1.72e-26)` 5.83·10⁹, 1.2's +∞; past
  σ/r = 2^170.7 the value or +∞, 1.2's ∞ − ∞ = NaN; below σ/r = 2^-179 with a huge ε a negative value or −0, 1.2's +0).
- **inorganic** — `born_lande_lattice_energy` evaluates 1.2's sequence with an unbounded exponent, the final product
  rounded once: physical salts keep their bits (the (1 − 1/n) factor is 1.2's), and the edges are right
  (`born_lande_lattice_energy(DBL_MAX, 1, 1, DBL_MAX, DBL_MAX)` −1.39·10⁻⁴, 1.2's −∞; a subnormal M −2.13·10⁻³¹⁸, 1.2's
  −0).
- **nuclear** — `binding_energy_per_nucleon` likewise, the final division rounded once: every mass below about
  1.93·10³⁰⁵ u keeps its bits, and past it A ≥ 932 brings the value back (`binding_energy_per_nucleon(1000, 1000, DBL_MAX)`
  −8.37·10³⁰⁷, 1.2's −∞). The errors and their order are 1.2's throughout.
  *Migration:* none for code that compares with a tolerance or stays within the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #33)
- **thermo, thermochem, nuclear, solution, inorganic, reaction, kinetics** — `hess_law`, `enthalpy_from_formation`,
  `reaction_enthalpy`, `reaction_gibbs_energy`, `reaction_entropy`, `q_value`, `ionic_strength`,
  `born_haber_lattice_energy`, `gibbs_free_energy`, `kohlrausch_molar_conductivity`, `limiting_molar_conductivity` and
  `zero_order_concentration` evaluate 1.2's operation sequence (its summation order, `.sum()` from −0.0 and
  thermochem's `+=` from +0.0) with an unbounded exponent, the final operation rounded once: every input whose 1.2 chain
  stays in the normal range keeps its bits (every physical input), and a partial sum or product past DBL_MAX that the
  rest brings back gives the value (`hess_law([DBL_MAX, DBL_MAX, −DBL_MAX])` DBL_MAX, 1.2's +∞;
  `reaction_enthalpy([Al₂O₃ ×10³⁰⁸] → [Al₂O₃ ×10³⁰⁸])` +0, 1.2's NaN; `q_value([M, M], [M, 0.999M])` 1.67·10³⁰⁸, 1.2's
  NaN; `ionic_strength([(0.4·DBL_MAX, 2)])` 1.44·10³⁰⁸, 1.2's +∞; `gibbs_free_energy(DBL_MAX, 2, 0.75·DBL_MAX)`
  −8.99·10³⁰⁷, 1.2's −∞; `zero_order_concentration(−DBL_MAX, −DBL_MAX, 1.5)` 8.99·10³⁰⁷, 1.2's +∞); a subnormal product
  is not rounded twice (`reaction_enthalpy` of 5·10⁻³²⁴ CO₂ is 5·10⁻³²⁴·(−393.51)·1000 rounded once), and a zero keeps
  the unbounded value's sign (`born_haber_lattice_energy(0, 0, 0, 2⁻¹⁰⁷⁴, 0)` −0, 1.2's +0; `zero_order_concentration`'s
  clamp decides on the unbounded difference). A non-finite input keeps 1.2's result; the errors and their order are
  1.2's.
  *Migration:* none for code that compares with a tolerance or stays within the normal range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #34)
- **thermochem** — `KmShomate.cp` and `enthalpy_antiderivative` evaluate 1.2's polynomials in 1.2's order with an
  unbounded exponent, the final addition (Cp) or product by 1000 (H) rounded once: every tabulated segment over its
  range and every physical use keeps its bits (so `enthalpy_change_cp` and `adiabatic_flame_temperature` are unchanged),
  and the edges of user coefficients are right — T → 0 (`Shomate{20.786, 0, 0, 0, 0}.cp(10⁻³⁰⁰)` 20.786, 1.2's NaN; its
  antiderivative at 5·10⁻³²⁴ K 1.04·10⁻³²², 1.2's NaN), large T (H₂O(g)'s Cp at 10³⁰⁰ K −∞, 1.2's ∞ − ∞ = NaN). T = ±0
  and non-finite inputs keep 1.2's results. *Migration:* none for temperatures within a fit's range.

### Breaking (relative to Rust 1.2.0, ADR 0002 #35)
- **reaction, solution** — `equilibrium_constant`, `debye_huckel_limiting`, `debye_huckel_extended` and
  `davies_activity` evaluate the argument of their exp or 10^x in 1.2's operation order with an unbounded exponent, its
  final quotient or product rounded once, then the correctly rounded exp or 10^x: every input whose 1.2 chain stays in
  the normal range keeps its result (every physical input), and an argument whose R·T, A·z², B·å, (B·å)·√I or numerator
  passes DBL_MAX or underflows gives the value (`equilibrium_constant(−DBL_MAX, DBL_MAX)` e^0.12 = 1.1278, 1.2's 1;
  `(−10⁻³²⁰, 10⁻³²²)` 193,216, 1.2's 197,357; `debye_huckel_limiting(100000, 0, 10³⁰⁰)` and `davies_activity` there 1,
  1.2's NaN; `debye_huckel_extended(1, 1, 10³⁰⁸, 10¹⁵⁵, 10¹⁵⁵)` 0.977, 1.2's 1; `(100000, 1, 10³⁰⁰, 10¹⁵⁴, 10¹⁵⁴)`
  10⁻¹⁰⁰, 1.2's 0). `equilibrium_constant`'s ComputationError is decided on the port's K, so at a subnormal T it can be
  the Err where 1.2's K is finite, or the reverse. The errors and their order are 1.2's.
  *Migration:* none for code that stays within the normal range (T above about 5·10⁻³⁰⁹ K and below 2·10³⁰⁷ K; |A|·z²,
  B·å and their products with √I below 10³⁰⁸).

### Breaking (relative to Rust 1.2.0, ADR 0002 #36)
- **spectroscopy** — `fft_spectrum` runs 1.2's transform (hisab 1.4.0's `num::fft`) and magnitudes with an unbounded
  exponent, each magnitude's √ rounded once: every signal whose 1.2 transform stays in the normal range keeps its
  spectrum, and a signal near DBL_MAX or subnormal has its spectrum where 1.2's sums, products or squares over- or
  underflow (`fft_spectrum([1.75·10³⁰⁴])` [1.75·10³⁰⁴], 1.2's [∞]; `([−6.96·10⁻³¹⁶])` [6.96·10⁻³¹⁶], 1.2's [0];
  `([3.3·10³⁰², 3.0·10³⁰⁴])` [3.03·10³⁰⁴, 2.96·10³⁰⁴], 1.2's [∞, ∞]). A signal with a NaN or ±∞ keeps 1.2's spectrum;
  the length check and its error are 1.2's. An aliased `signal` (inside out's first 2n values) is no longer transformed
  in `out`: it gets a 16·n-byte heap table, and the slow path 32·n bytes, from the bump allocator.
  *Migration:* none for samples between about 10⁻²⁷⁷ and DBL_MAX/n in magnitude (or exact zeros); a caller that
  relied on no allocation with `signal == out` passes separate buffers.

### Breaking (relative to Rust 1.2.0, ADR 0002 #37)
- **spectroscopy** — `find_peaks` compares a candidate with the exact threshold min + t·(max − min) (1 − t not
  rounded), where 1.2 rounds min·(1 − t) + max·t first: a candidate within a few ulps of the threshold, or at the
  subnormal grid, can change (`find_peaks([−0, 5·10⁻³²⁴, 0], 0.75)` [1], 1.2's []; `([10, 10.2, 10, 11, 10],
  1 − 2⁻⁵³)` [3], 1.2's []), as can one within 2⁻⁵³·max(|min|, |max|) of it where min and max nearly cancel (±DBL_MAX
  at t = ½ ± an ulp). The checks, the empty results, min / max over the finite values and the neighbour tests are 1.2's.
  *Migration:* none unless a threshold is meant to be 1.2's rounded value.

### Breaking (relative to Rust 1.2.0, ADR 0002 #38)
- **thermochem** — `adiabatic_flame_temperature` tests its bracket by the signs of f(T₀) and f(T_max), not their
  product, which underflows to −0 when |ΔH| and the moles are tiny: such a system (a real one scaled by 10⁻²⁰⁰, or
  (−10⁻³²⁰, [CO₂ ×10⁻³⁰⁰], 298.15)) is Ok — Newton's first guess, accepted by the absolute 10⁻³ J/mol tolerance, 1.2's own
  answer at any scale where the product does not underflow — where 1.2 gives Err(ComputationError). Every other input
  is unchanged (proven: the tests agree wherever the product is not an underflowed zero).
  *Migration:* none; code that relied on the ComputationError for tiny systems checks the scale itself.

### Breaking (relative to Rust 1.2.0, ADR 0002 #39)
- **stoichiometry** — `balance_equation` balances on unbounded integers where 1.2's checked i64 elimination overflows:
  a matrix whose null space is one-dimensional with a positive u32 vector is Ok with it, where 1.2 returns
  Err(ComputationError) "integer overflow while balancing equation" (`[[3037000500, 1, −3037000501], [1, 3037000500,
  −3037000501]]` [1, 1, 1]: 1.2's first product, 3037000500², passes i64::MAX). Every other result is 1.2's: where the
  exact null space has another dimension, or its vector a zero, mixed signs or an entry past u32::MAX, the result is
  still the ComputationError (now for that reason), and every row 1.2 computes without overflow is unchanged. Such an
  input with too small a `cap`, or a null `out`, now reaches the buffer check: Err(InvalidInput). The multi-word path
  allocates its tables (about 8·(m·n·B/31 + 4n·B/31) bytes, B the matrix's Hadamard bound in bits) past a 4,096-word
  local table, from the bump allocator.
  *Migration:* none for matrices 1.2 balances; code that relied on the overflow Err to reject huge coefficients checks
  the matrix's entries itself.

### Changed
- `VERSION` 2.0.0 (`cyrius.cyml` reads it through `${file:VERSION}`; `dist/kimiya.cyr` regenerated); `README.md`
  describes the Cyrius library: its modules, its use from a `cyrius.cyml`, a Quick Start that compiles, and the API's
  mapping from Rust
- Rust crate, CI, cargo scripts, benchmarks and architecture overview moved to `rust-old/`; CI is now Cyrius
  (`test`, `test --bounds`, `deps --locked`, `audit`)
- **Results differ from Rust 1.2.0 in the last bit wherever glibc misrounds** (M3c, ADR 0006). 1.2.0's transcendentals
  are the host glibc's (2.44 for the oracle), which is not correctly rounded; 2.0's are correctly rounded on every
  target, so a Rust 1.x user comparing results will see a 1-ulp difference in the primitive on the calls where glibc is
  off, and whatever the function's conditioning makes of it. Over the golden calls (port-parity.md §9) glibc misrounds
  expm1 on 4.2 % of the `math` rows (109 of 1,267 distinct arguments in the module rows), log10 on 3.7 % (215 of 11,027;
  up to 29 % of arguments in [½, 1)), exp on 10 of 12,506 distinct arguments, log on 5, cos on 5 (0.16 % of random
  arguments), acos 0.07 %, pow(10, y) on 3, pow on 2 — about 0.1–5 % depending on the primitive. Two kinds of example:
  `10.0_f64.powf(23.0)` is exactly the midpoint 10²³, which 2.0 rounds to even, `1e23` (`h_concentration_from_ph(-23)`),
  where glibc returns the next double up; and the exact subnormal ties (3·2⁻²¹⁵)⁵ and (9·2⁻⁴³⁰)^2.5 = 121.5·2⁻¹⁰⁷⁴,
  which 2.0 rounds to 122·2⁻¹⁰⁷⁴ where glibc gives 121·2⁻¹⁰⁷⁴, so `nth_order_half_life` on those inputs differs by
  4.9 × 10¹³ ulp (2.0 is the correct one). On ill-conditioned rows a 1-ulp primitive difference is amplified: the
  Henderson–Hasselbalch row with base/acid = 0.8068…/2.6858… gives pH −2⁻⁵³ where 1.2.0 gave +0.0. Results no longer
  depend on the platform's libm: x86_64 and aarch64 give the same bits. Since M4d sin, cos, atan2 and acos are correctly
  rounded too (glibc misrounds 2 of the 710 `math` sin rows, 4 of the 2,910 cos rows, 1 of 2,620 acos rows and 1 of the
  4,439 potential atan2 calls); a bond angle on such a call differs from 1.2's in the last bit. Since M4e the
  Beer–Lambert `transmittance` (pow(10, ·): 3 golden rows, 10²³ among them) and `absorbance_from_transmittance` (log10:
  115 rows) differ there too; the FFT's twiddles are correctly rounded (glibc's are at the 13 golden angles, so no
  `fft_spectrum` row changes)
- **spectroscopy** — the isotope table is the cited one ([ADR 0002 #23](docs/adr/0002-port-api-divergences.md), M4e-4):
  `isotope_abundances` returns the AME2020 masses at full published precision (¹H 1.007825031898 u where 1.2 had
  1.00783, ²H 2.014101777844 for 2.01410, ³⁴S 33.967867011 for 33.96787) and sulfur's fourth isotope, ³⁶S (35.967080692 u,
  0.0001), so every element's abundances sum to 1 (1.2's sulfur: 0.9999); the abundances are 1.2's, the IUPAC 2009
  (TICE 2009) representative compositions. `isotope_pattern`'s algorithm is unchanged, so its results change wherever
  the molecule has H, C, N, O, S, Cl or Br: the masses by the rounding 1.2's table dropped (Cl₂'s M is 69.937705388, 1.2
  69.9377), and with sulfur the intensities from M+4 on (S₈'s M+4 0.0573616, 0.0565194 without ³⁶S). Errors are 1.2's.
  *Migration:* none for code that compares masses to 10⁻⁴ u; code that tested 1.2's exact bits or sulfur's three pairs
  sees the new table

### Removed
- **mathx** (private, M3c-3): `_km_log10` and `_km_expm1` (M1's Rust libm 0.2.16 / fdlibm transcriptions),
  `_km_log10_ln` (M3a's glibc log10 composition, ADR 0005 amendment M3a, superseded by ADR 0006), the `_km_ln` wrapper,
  `_km_pow_2_5` (ADR 0005's double-double x^2.5, correctly rounded only from 2⁻³⁸⁷) and the unused `_km_dd_div`; their
  Cyrius-only tests now pin the correctly rounded bits (`tests/solution.tcyr`'s Henderson–Hasselbalch group: log10
  `…69DB`, glibc's `…69DA`)
- `validate.py`'s `_km_log10_ln` sensitivity model, its "S before" mode and the e_log10.c composition check
- **mathx** (private, M4d): ADR 0005's double-double `_km_cos` and `_km_acos` (M2a) and their helpers `_km_dd_cos_sin`,
  `_km_dd_asin`, `_km_dd_add`, `_km_dd_mul`, `_km_dd_div_d` and `_km_dd_sqrt`, once `gas.cyr` moved to M4a's
  `_km_cr_cos` / `_km_cr_acos` (the same bits on every golden row); `src/mathx.cyr` no longer includes `lib/math.cyr`.
  Their tests moved to the replacements (`tests/gas.tcyr`), their benches and `cr_verify.py`'s comparison went with them
- **stoichiometry** (private, M4f-2): `_km_empirical_close`, Rust's |scaled − round(scaled)| ≤ 0.1 on a rounded ratio,
  replaced by ADR 0002 #26's exact tolerance test; its Cyrius-only test group became `cyrius_adr26_exact_tolerance`

### Fixed
- **electrochemistry** — `nernst_potential`'s fast path (ADR 0002 #29, the overflow audit's B2) accepted an intermediate
  of exactly 2⁻¹⁰²², which can be the IEEE rounding of an exact value just below it: for a product c·ln Q in
  [2⁻¹⁰²² − 2⁻¹⁰⁷⁵, 2⁻¹⁰²² − 2⁻¹⁰⁷⁶) it returned 1.2's value where #20's semantics differ (with E° an odd multiple of
  2⁻¹⁰⁷⁴: `nernst_potential(5·10⁻³²⁴, 1, 5.1625987955063826·10⁻³⁰³, 1.0512872158898858)` is −0x0.ffffffffffffep-1022, B2
  returned −0x0.fffffffffffffp-1022). The test's bound is now 2⁻¹⁰²¹, `_km_ux_normal`'s; 23 golden edge rows (17 noted)
  and `tests/electrochemistry.tcyr`'s `cyrius_adr29_2m1022` cover it, and `validate.py`'s shared tracer is exact there.
  B1's to B3's other fast paths were audited for the same bound: none has it (ADR 0002's B4 addition)
- **mathx** — a comment on `osmotic_pressure`'s factors said R's odd part has 53 bits; it has 51 (the conclusion, no
  midpoint, holds: with 1000's 7 bits the product has at least 57)

### Performance
- `molecule/molecular_weight_glucose`: 121 ns → 26.7 ns (−78%; Rust 3.63 ns) by walking the atom buffer and
  inlining the element range check; `element/lookup_by_number`: 7.03 → 4.54 ns (−35%; Rust 0.84 ns);
  `element/lookup_by_symbol`: 92.4 → 67.0 ns (−27%, first-byte compare; Rust 44.73 ns)
- validate / mathx guards without nested calls (`_km_is_finite` as `x − x == 0`, `_km_is_positive` /
  `_km_is_non_negative` against +∞): `reaction/gibbs_free_energy` 23.3 → 16.9 ns (−27%; Rust 1.75 ns),
  `gas/ideal_gas_pressure` 26.7 → 16.7 ns (−37%; Rust 2.18 ns), `gas/peng_robinson_z_n2_1mpa` 289 → 270 ns
  (−7%; Rust 61.39 ns)
- M2b thermochem, first measurement (Rust 1.2.0 → Cyrius): `thermochem/reaction_enthalpy_ch4` 38.58 ns → 250 ns,
  `thermochem/vant_hoff_k` 9.51 → 58.9 ns, `thermochem/enthalpy_change_cp` 24.00 → 166 ns, `thermochem/shomate_cp`
  2.34 → 14.5 ns, `thermochem/adiabatic_flame_ch4_air` 301.48 ns → 3.24 µs (6.2–10.7×; two optimizations were
  tried and not kept, [`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md))
- M3a solution / electrochemistry, first measurement (Rust 1.2.0 → Cyrius): `solution/ph_from_h` 9.73 → 66.2 ns,
  `solution/weak_acid_ph` 34.85 → 182 ns, `electrochemistry/lookup_half_reaction` 57.06 → 116 ns (2.0–6.8×);
  `electrochemistry/nernst_potential` 6.42 ns → 1.94 µs (302×): the correctly rounded `_km_ln` (1.90 µs, against 23 ns
  for `f64_ln`) that its cancellation needs (ADR 0005 amendment M3a)
- M3b kinetics, first measurement (Rust 1.2.0 → Cyrius): `kinetics/arrhenius_rate` 6.78 → 49.4 ns (7.3×),
  `kinetics/eyring_rate` 6.73 → 42.4 ns (6.3×), `kinetics/michaelis_menten` 2.92 → 16.5 ns (5.7×)
- M3c correct rounding (Cyrius before → after on the same host and day; Rust 1.2.0 in brackets;
  [`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md)): `electrochemistry/nernst_potential`
  1.94 µs → 84 ns (−96 %, M3c-2a's two-phase `_km_cr_ln` in place of ADR 0005's double-double ln; Rust 6.42 ns);
  the cost of correct rounding elsewhere: `thermochem/vant_hoff_k` 65.5 → 77.1 ns (Rust 9.51), `kinetics/arrhenius_rate`
  52.0 → 65.5 ns (Rust 6.78), `kinetics/eyring_rate` 44.6 → 59.0 ns (Rust 6.73) — `_km_cr_exp` 38.7 ns against about
  11 ns for `f64_exp`; `solution/ph_from_h` 69.9 → 87.3 ns (Rust 9.73), `solution/weak_acid_ph` 196 → 219 ns (Rust
  34.85) — `_km_cr_log10` 78.7 ns against 59 ns for fdlibm's; pow faster: `_km_cr_pow` 138 / 181 ns against
  `ganita_f64_pow` 198 / 218 ns on the two benched pairs. Every other bench unchanged within noise
- M4b organic / inorganic, first measurement (Cyrius only: rust-old/benches has no organic or inorganic bench):
  `organic/enthalpy_from_bonds` 49.5 ns (CH₄ + 2O₂, two lists of two records), `inorganic/lookup_ionic_radius` 107 ns
  (Sn²⁺, a scan to the last of the 30 rows)
- M4c nuclear, first measurement (Rust 1.2.0 → Cyrius): `nuclear/binding_energy_per_nucleon_fe56` 2.02 → 27 ns (13×;
  41 → 29 ns from calling a private, unchecked `_km_mass_defect` and multiplying by m_u·c² inline, same host back to
  back), `nuclear/semi_empirical_fe56` 18.50 → 131 ns (7.1×; ganita's correctly rounded cbrt); Cyrius-only
  `nuclear/decay_chain_two_step` 256 ns (two correctly rounded exps and two expm1s)
- M4a correctly rounded trigonometry (Cyrius-only; [`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md)):
  `mathx/cr_sin` 74 ns against `f64_sin` 37 ns (the stdlib's, not correctly rounded), `mathx/cr_cos` 82 against `f64_cos`
  40, `mathx/cr_atan2` 113 against `ganita_f64_atan2` 53, `mathx/cr_acos` 140 ns against ADR 0005's double-double
  `_km_acos` 5.8 µs (41× faster), and `_km_cos` was 1.68 µs (`cr_cos` 20× faster); `mathx/cr_sincos` 115 ns for both
- M4d potential, first measurement (Rust 1.2.0 → Cyrius): `potential/dihedral_angle` 20.35 → 130 ns (6.4×; Rust's gauche
  geometry, whose atan2 is the x = 0 special case); Cyrius-only `potential/bond_angle` 215 ns (a full correctly rounded
  atan2), `potential/lennard_jones` 35 ns, `potential/lj_equilibrium_distance` 109 ns (the correctly rounded σ·2^(1/6) of
  ADR 0002 #16), `potential/morse` 66 ns, `mathx/cr_exp2_scaled` 151 ns (ADR 0002 #13's N₀·2⁻ⁿ, 1000 after 1075.3
  half-lives). The Peng–Robinson trig roots' acos and cos now cost 142 + 3 × 81 ns instead of 5.8 µs + 3 × 1.68 µs

- M4e spectroscopy, first measurement (Rust 1.2.0 → Cyrius): `spectroscopy/absorbance` 0.98 → 11 ns (11×),
  `spectroscopy/bohr_energy_level` 1.21 → 16 ns (13×; Rust inlines both into a few instructions),
  `spectroscopy/isotope_pattern_glucose` 2.80 → 16.8 µs (6.0×; 24 convolution steps on the stack tables, nothing
  allocated); Cyrius-only `spectroscopy/fft_spectrum_256` 28.2 µs, `_1024` 134 µs, `_4096` 628 µs (in place),
  `spectroscopy/find_peaks_1000` 22.5 µs, `spectroscopy/wavenumber_from_wavelength` 58 ns (ADR 0002 #19's correctly
  rounded 1/(100·λ), the fast path)
- M4e-2, ADR 0002 #20's fast-path test (Cyrius before → after, the same host, back to back): `potential/dihedral_angle`
  131 → 133 ns (Rust 20.35; the per-component window on b₁, b₂, b₃, with the three |b|² computed inline),
  `potential/bond_angle` 213 → 213 ns (the sum-of-squares window replaces Rust's two `is_positive` calls),
  `potential/distance_3d` 14 → 16 ns (new; Cyrius-only); the slow path on the same geometries scaled by 2⁶⁰⁰, where Rust's
  squares overflow: `potential/distance_3d_2p600` 314 ns, `bond_angle_2p600` 1.13 µs, `dihedral_angle_2p600` 1.31 µs
  (1.2's chain gave +∞ or an Err there in 14, 52 and 80 ns); `potential/bond_angle_linear` (CO₂ on an axis, ba × bc
  exactly 0) 1.06 µs → 135 ns, by a proven shortcut for an exactly linear or zero-angle geometry (it had taken the slow
  path)
- M4e-3, ADR 0002 #21 and #22 (Cyrius before → after, the same host, back to back): `spectroscopy/bohr_energy_level`
  9 → 49 ns (Rust 1.21; R·Z² split exactly, the quotient by n² inline with its rounding test), and the Cyrius-only
  `spectroscopy/bohr_transition_energy` 27 → 73 ns, `rydberg_wavelength` 24 → 73, `series_wavelength_balmer_alpha`
  28 → 78, `reduced_mass_hcl` 14 → 73 (double-double operands and `_km_dd_div_cr`); the exact path, on a tie:
  `bohr_transition_energy_tie` 2.42 µs, `reduced_mass_tie` 2.20 µs (1.2's chain took 28 and 14 ns there)
- M4e-4, ADR 0002 #23's table (Cyrius, the old table against the new, the same host, interleaved):
  `spectroscopy/isotope_pattern_glucose` 15.8 µs either way (Rust 2.80); the new Cyrius-only sulfur benches
  `isotope_pattern_cystine` 16.9 → 18.5 µs (+10 %) and `isotope_pattern_s8` 8.4 → 18.8 µs (2.2×): ³⁶S's fourth isotope
  makes each sulfur step's scratch a third longer and its pattern wider (S₈ up to 17 peaks a step instead of 15; the
  scratch at most 68 records, within the 192-record stack tables)

- M4f stoichiometry (Cyrius before → after the proven fast paths, the same host, back to back):
  `stoichiometry/balance_ch4_combustion` 1.20 µs → 770 ns (Rust 233.90 ns), `oxidation_states_kmno4` 156 → 122 ns (Rust
  22.38); the Cyrius-only `balance_kmno4_hcl` 3.84 → 2.16 µs, `oxidation_states_fe3o4` 132 → 110 ns, `limiting_reagent`
  31 ns, `theoretical_yield_moles` 46, `percent_yield` 14, `percent_composition_glucose` 71, `empirical_formula_glucose`
  217 ([`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md))
- M4f-2, ADR 0002 #24–#26 (Cyrius before → after, the same host, back to back; Cyrius-only): `stoichiometry/percent_yield`
  14 → 81 ns, `theoretical_yield_moles` 46 → 115, `percent_composition_glucose` 72 → 276 (three entries),
  `empirical_formula_glucose` 217 → 586, `limiting_reagent` 32 → 33 (the rounded ratios still decide unequal ones); the
  exact paths: `theoretical_yield_moles_tie` 133 ns (a midpoint resolved in the fast path), `percent_yield_subnormal`
  2.15 µs (`_km_cr_ratio`), `limiting_reagent_tie` 94 ns, `empirical_formula_edge` 7.46 µs (two exact tolerance tests).
  Inlining the fast path's helpers was measured (81 → 73 ns) and not kept: the cost is its ~45 floating-point
  operations, as for M4e-3's reduced mass
- The overflow audit's B1, ADR 0002 #27 and #28 (Cyrius before → after, the same host, back to back): `thermo/heat_transfer`
  5.4 → 65 ns (Rust 0.95), `kinetics/michaelis_menten` 17 → 24 ns (Rust 2.92); Cyrius-only
  `specific_heat_from_calorimetry` 18 → 86, `boiling_point_elevation` 6.1 → 66, `freezing_point_depression` 6.4 → 67,
  `osmotic_pressure` 20 → 139, `steady_state_intermediate` 21 → 88, `pre_equilibrium_rate` 27 → 125,
  `zero_order_half_life` 16 → 16, `second_order_concentration` 23 → 33; the other paths `heat_transfer_tie` 365 ns (a
  midpoint, `_km_mono_tie`), `osmotic_pressure_subnormal` 4.04 µs (the exact path), `second_order_concentration_edge`
  230 ns and `michaelis_menten_edge` 126 ns (the unbounded-exponent slow paths); stoichiometry within noise of M4f-2
  (`percent_yield` 82, `theoretical_yield_moles` 115, `percent_composition_glucose` 266, `percent_yield_subnormal`
  2.19 → 1.62 µs). Routing every shape through the general `_km_cr_mono` was measured (`heat_transfer` 148 ns,
  `percent_yield` 108, `percent_composition_glucose` 422) and replaced by straight-line fast paths per shape
- The overflow audit's B2, ADR 0002 #29 and #30 (Cyrius before → after, the same host, back to back):
  `spectroscopy/absorbance` 5 → 68 ns (Rust 0.98), `electrochemistry/nernst_potential` 83 → 92 ns (Rust 6.42: the range
  test); Cyrius-only `dilution` 17 → 86, `solution_conductivity` 6 → 72, `faraday_mass_deposited` 19 → 109,
  `faraday_charge_required` 21 → 122, `faraday_mass_ratio` 33 → 148, `concentration_from_absorbance` 16 → 90,
  `vibrational_frequency` 18 → 87, `rigid_rotor_energy` 14 → 114, `rotational_constant` 18 → 117,
  `diatomic_moment_of_inertia` 11 → 74; the other paths `nernst_potential_edge` 90 → 202 ns (R·T past DBL_MAX: the
  unbounded-exponent slow path), `vibrational_frequency_edge` 19 → 88 ns (k/μ past DBL_MAX: the √ fast path),
  `vibrational_frequency_subnormal` 26 → 883 ns (a subnormal μ: the exact root), `faraday_mass_ratio_tie` 33 → 402 ns (a
  midpoint, `_km_mono_tie`). Inlining the √ fast path's two TwoProducts was measured (87 → 86 ns) and not kept
- The overflow audit's B3, ADR 0002 #31 and #32 (Cyrius before → after, the same host, back to back):
  `gas/ideal_gas_pressure` 18 → 119 ns (Rust 2.18), `nuclear/binding_energy_per_nucleon_fe56` 28 → 29 ns (Rust 2.02),
  `potential/lennard_jones` 36 → 39 ns; Cyrius-only `ideal_gas_volume` 20 → 119, `gas_density` 20 → 108,
  `van_der_waals_pressure` 35 → 55, `peng_robinson_pressure` 78 → 100, `redlich_kwong_pressure` 182 → 187, `coulomb`
  16 → 126, `coulomb_force` 18 → 145, `morse_width_parameter` 33 → 135, `born_lande_lattice_energy` 46 → 50; the slow
  paths `van_der_waals_pressure_edge` 36 → 260 ns, `peng_robinson_pressure_edge` 115 → 618 ns,
  `redlich_kwong_pressure_edge` 134 → 584 ns, `lennard_jones_edge` 46 → 208 ns, `born_lande_lattice_energy_edge`
  66 → 287 ns, `binding_energy_per_nucleon_edge` 39 → 117 ns, and the Morse width's exact root
  (`morse_width_parameter_subnormal` 42 ns → 1.27 µs, `morse_width_parameter_tie` 35 ns → 1.27 µs). A region test on the
  inputs (Tc, Pc, T, V within 2^±60: every intermediate provably normal) lets the equations of state and Born–Landé skip
  their per-operation range tests on physical states: Peng–Robinson 147 → 100 ns, van der Waals 73 → 55, Redlich–Kwong
  217 → 187, Born–Landé 67 → 50; inlining `_km_ux_normal`'s bit test was measured (147 → 156 ns) and not kept. B3's
  seventeen benches run 100,000 × 20 calls in `tests/kimiya.bcyr`: at 1,000,000 × 20 the file passed `cyrius bench`'s
  default 300 s deadline (255 s before B3, 279 s now)
- The overflow audit's B4, ADR 0002 #33 and #34 (Cyrius before → after, the same host): `reaction/gibbs_free_energy` 17
  → 18 ns (Rust 1.75), `thermochem/reaction_enthalpy_ch4` 259 → 293 ns (Rust 38.58), `thermochem/shomate_cp` 15 → 32 ns
  (Rust 2.34: the box test, inlining it measured 34 and 39 ns and not kept), `thermochem/enthalpy_change_cp` 171 → 171
  ns (Rust 24.00) and `thermochem/adiabatic_flame_ch4_air` 3.40 → 3.29 µs (Rust 301.48 ns: the table-row path tests T
  alone, inline), `electrochemistry/nernst_potential` 90 → 89 ns (Rust 6.42; #29's bound); Cyrius-only `hess_law` 13 →
  15, `enthalpy_from_formation` 20 → 24, `reaction_gibbs_energy` 259 → 299, `reaction_entropy` 259 → 286, `q_value` 24 →
  26, `ionic_strength` 24 → 31, `born_haber_lattice_energy` 10 → 14, `kohlrausch_molar_conductivity` 22 → 28,
  `limiting_molar_conductivity` 10 → 24, `zero_order_concentration` 12 → 17, `shomate_enthalpy_antiderivative` 37 → 54
  ns; the slow paths 117–587 ns (`hess_law_edge` 162, `reaction_enthalpy_edge` 587, `gibbs_free_energy_edge` 117), the
  Shomate polynomials' 311 and 474 ns. The bench files split (the user's decision): `tests/kimiya.bcyr` 177 s (261 s
  before), the new `tests/audit.bcyr` 189 s, B3's seventeen back at 1,000,000 × 20 calls with the same per-call figures
- The overflow audit's B5, ADR 0002 #35 (Cyrius before → after, the same host; `bench_audit_b5` in `tests/audit2.bcyr`,
  no Rust bench): the fast tests are inline compares, so the typical inputs keep their cost —
  `reaction/equilibrium_constant` 64 → 64 ns, `solution/debye_huckel_limiting` 86 → 85, `debye_huckel_extended` 105 →
  100, `davies_activity` 98 → 96 ns; the slow paths `equilibrium_constant_edge` 26 → 151 ns (Rust's exp(−0) is trivial),
  `equilibrium_constant_subnormal` 80 → 166, `debye_huckel_limiting_edge` 30 → 93, `debye_huckel_extended_edge` 45 → 277,
  `debye_huckel_extended_overflow` 47 → 276, `davies_activity_edge` 39 → 162 ns; `tests/audit2.bcyr` 31 s, `tests/audit.bcyr` unchanged at 189 s (219 s with them)
- The overflow audit's B6, ADR 0002 #36 … #38 (Cyrius before → after, the same host, alternating runs; no Rust bench):
  `spectroscopy/fft_spectrum_256` 27.3 → 30.5 µs, `fft_spectrum_1024` 129.2 → 143.1, `fft_spectrum_4096` 605.6 → 669.4
  (+11–12 %: the per-butterfly operand checks about 6 %, the restructured magnitude pass the rest; skipping the checks at
  k = 0 gained under 1 % and was not kept), `find_peaks_1000` 19.0 → 20.3 µs, `thermochem/adiabatic_flame_ch4_air` 3.37 →
  3.40 µs (Rust 301.48 ns); `bench_audit_b6`: `fft_spectrum_256_huge` 27.7 → 318.1 µs and `_subnormal` 35.4 → 349.8 µs
  (the unbounded transform), `_squares` 27.4 → 45.1 µs (the per-bin unbounded √), `find_peaks_1000_subnormal` 19.5 →
  39.0 µs, `adiabatic_flame_tiny` 0.96 → 1.48 µs (Rust's Err there); `tests/audit2.bcyr` 45 s (31 s before)
- The overflow audit's B7, ADR 0002 #39 (Cyrius before → after, the same host, alternating runs): the i64 fast path
  unchanged — `stoichiometry/balance_ch4_combustion` 805 → 778 ns (Rust 233.90 ns), `balance_kmno4_hcl` 2.24 →
  2.18 µs; `bench_audit_b7` (Rust's
  ComputationError at its first overflow, then the port's exact result): `balance_survey_2x3` 0.12 → 3.6 µs,
  `balance_overflow_10x9` 6.6 → 57 µs, `balance_dense_30x31` 5.4 µs → 158 ms (products past 3,000 bits; Lehmer's gcd
  took a 39 × 40 matrix from 1.8 to 0.7 s); `tests/audit2.bcyr` 49 s (45 s before)
- The overflow audit's B10a (Cyrius-only; `bench_audit_b10a`, the binary64 function → the unbounded one at the same
  argument): e^x at −50000/(8.314·300) `_km_cr_exp` 46 ns → `_km_cr_exp_ux` 54 ns, e^x − 1 at ln 2 54 → 65 ns, 10^x at
  −4.75 60 → 69 ns (the first tier: the binary64 function and its pair); e^-720.5 (a binary64 subnormal, the cores
  unscaled) 84 ns; e^(−10⁴), e^(10¹²), e^(10⁴) − 1 and 10^(10⁴) 158–165 ns (the double-double reduction); e^(10¹⁵) 614 ns
  (the integer reduction); the accurate path alone 1.86 µs at −20.05, 2.14 µs at −10⁴ and 2.17 µs for 10^(−10¹²);
  `tests/audit2.bcyr` 66 s (49 s before)
- M5a, first measurement (Rust 1.2.0 → Cyrius, `tests/numerics.bcyr`): `dynamics/first_order_100steps` 4.41 → 11.2 µs
  (2.5×), `dynamics/consecutive_100steps` 5.04 → 21.4 µs (4.3×); Cyrius-only against an ad hoc Rust timing of the same
  bodies: `dynamics/adaptive_first_order` 1.80 → 6.19 µs, `dynamics/stiff_robertson_400steps` 121 → 498 µs,
  `dynamics/stochastic_100steps` 4.85 → 24.8 µs, `numerics/gaussian_elimination_8x8` 353 ns → 2.48 µs (the row-pointer addressing
  took it from 3.15 to 2.56 µs in one side-by-side run), `numerics/pcg32_next_normal` 22.6 → 200 ns (the correctly rounded ln and cos)
  ([`docs/benchmarks-rust-v-cyrius.md`](docs/benchmarks-rust-v-cyrius.md))

## [1.2.0] - 2026-11-08

Refresh + audit release: current toolchain and dependencies, then a full-crate audit
and repair sweep in preparation for the Cyrius port. No public signature changed; new
items are additive. Many functions now return corrected values or reject input they
used to accept — see **Behaviour changes** before upgrading.

### Changed — toolchain and dependencies
- MSRV 1.89 → **1.99** (edition 2024); CI MSRV job updated
- hisab 1.1 → **1.4**, now `default-features = false, features = ["num", "autodiff"]` (the `transforms`, `geo` and `calc` modules were unused)
- criterion 0.5 → **0.8** (benchmarks use `std::hint::black_box`); serde 1.0.229, thiserror 2.0.21, tracing-subscriber 0.3.23, serde_json 1.0.154
- CI: actions/checkout v7, codecov-action v7, action-gh-release v3; cargo-audit and cargo-llvm-cov installed as prebuilt binaries (taiki-e/install-action)
- deny.toml: `unused-allowed-license = "allow"` (clean under cargo-deny 0.20)
- `#![warn(missing_docs)]` enabled; every public item documented with units
- README rewritten (all modules, features, dependency stack); its example now runs as a doctest

### Removed
- Unused direct `tracing` dependency (no tracing calls in kimiya; `logging` still uses tracing-subscriber)

### Added
- **gas** — `peng_robinson_z_roots` (all physical Z roots, ascending)
- **thermochem** — `lookup_shomate_at(formula, T)` and `Shomate::enthalpy_antiderivative`; multi-segment NIST Shomate data (N₂, O₂, CO₂, H₂O, H₂, CO, CH₄ to 6000 K); Zn(s), Ti(s), Ba(s), Ag(s), F₂(g) standard-state rows
- **nuclear** — `HYDROGEN_ATOM_MASS_U` (AME2020)
- **solution** — `KW_WATER_25C`
- `validate` (crate-internal) — NaN/±∞-safe guard predicates used by every module
- Benchmarks: balance_equation, oxidation states, isotope_pattern, binding energy, SEMF, Peng–Robinson Z, adiabatic flame, dihedral angle

### Fixed — physics and data
- **nuclear** — `mass_defect` mixed the bare proton mass with neutral-atom masses, so every binding energy was low by Z·mₑc² (He-4 27.27 → 28.296 MeV; Fe-56 BE/A 8.554 → 8.790 MeV); now uses m(¹H)
- **nuclear** — `AMU_MEV` mistyped (931.494 102 22 → 931.494 102 42, CODATA 2018); `PROTON_MASS_U` was CODATA 2014 → 2018
- **nuclear** — semi-empirical mass formula used Povh coefficients with a Z(Z−1) Coulomb term and δ = 12; now self-consistent (Z², δ = 11.2): Fe-56 498.97 → 494.0 MeV (exp. 492.26)
- **nuclear** — half-lives updated to NUBASE2020: I-131 8.0252 d, Pu-239 24 110 y, C-14 5700 y
- **thermochem** — N₂ and O₂ Shomate coefficients were the low-temperature NIST sets mislabelled to 1400/1200 K (N₂ Cp(2000 K) was 127 J/(mol·K); O₂ went negative); CO₂/H₂O lacked high-T segments. Adiabatic flame temperatures were ~550 K low: H₂/air 1974 → 2525 K, CH₄/air 1865 → 2327 K
- **thermochem** — C₂H₄ and Br₂(g) Shomate sets replaced with NIST's; K₂O S° 102.00 → 94.10 (CRC); n-butane row → NBS (ΔHf −126.15, ΔGf −17.03, S° 310.23)
- **thermochem** — `enthalpy_change_cp` / `adiabatic_flame_temperature` silently extrapolated outside the data range; now integrate piecewise (analytic) and return `Err` outside it. A NaN ΔH used to return `Ok(≈299.5 K)`
- **solution** — `weak_acid_ph` / `weak_base_ph` ignored water autoionisation (1e-8 M acid gave pH 8.0, now 6.978) and Newton's absolute tolerance could return the initial guess (Ka 5e-7, C 2e-6: 6.00 → 6.103); now solves the full charge-balance cubic
- **stoichiometry** — `balance_equation` truncated the system to n−1 element rows, so results depended on row order (CuSO₄ + Fe → FeSO₄ + Cu failed in one order), and hid wrong-sign solutions. Rewritten as exact integer null-space elimination with verification
- **stoichiometry** — `assign_oxidation_states` truncated fractional states (propane C −2, Fe₃O₄ +2; now `Err`) and lacked the documented peroxide, OF₂ and metal-hydride rules (H₂O₂, Na₂O₂, OF₂, NaH, LiAlH₄ now correct)
- **stoichiometry** — `empirical_formula` could pick a larger multiplier than needed or return a rounded guess; now the smallest multiplier within 0.1, else `Err`
- **inorganic** — low-spin d⁹ had 3 unpaired electrons (now 1: strong-field Cu²⁺ 3.87 → 1.73 BM); Sn²⁺ radius → Shannon CN 8 (122 pm)
- **gas** — Peng–Robinson Newton from Z = 1 converged to negative roots at high pressure (N₂ 300 K / 90 MPa: −1.65 vs 1.66) and fugacity returned NaN; now analytic cubic roots with Z > B
- **potential** — `morse_width_parameter` used the ordinary-frequency formula with an angular frequency (2π too large; HCl 11.7 → 1.87 Å⁻¹); dihedral/bond-angle degeneracy checks were absolute and rejected SI-metre coordinates
- **spectroscopy** — `isotope_pattern` panicked on large molecules (underflow → NaN → `partial_cmp().unwrap()`; e.g. C₁₆₈); hydrogen `SpectralSeries` now uses R_H (Hα 656.112 → 656.470 nm); H isotope abundances → IUPAC 2009
- **electrochemistry** — table sorted (it was two concatenated runs); K⁺/K −2.924 → −2.931, V²⁺/V −1.130 → −1.175 (CRC); `FARADAY` → exact N_A·e (96 485.332 123 310 01)
- **element** — Li 6.941 → 6.94, Si 28.086 → 28.085, S 32.065 → 32.06 (IUPAC conventional), Zr 91.224 → 91.222 (CIAAW 2024)
- **phase** — mmHg factor 133.322 → 133.322 387 415 Pa; methanol Antoine coefficients and range made consistent (NIST, 15–84 °C)
- **bridge** — used R = 8.314 and F = 96 485.332 instead of the crate constants; `atomic_number_to_valence` returned 2 for every Z ≥ 19
- **reaction_dynamics** — `simulate_stiff` (via hisab `backward_euler`) accepted unconverged Newton iterates and returned `Ok`; now its own backward-Euler step that errors on non-convergence
- **examples/basic.rs** — Nernst example used Q = [Cu²⁺] instead of 1/[Cu²⁺] (0.401 → 0.283 V)

### Fixed — robustness
- Guards written `x <= 0.0` let NaN and ±∞ through and returned `Ok(NaN)`; every validated input in every module now rejects non-finite values
- Removed every `unwrap`/`expect`/`unreachable!` from library code; integer overflow paths (u32 counts, i32 charges, `phase_rule`, SEMF `a − z`, series line index, degree + 1) use checked/saturating/widened arithmetic
- `logging::init()` panicked if a global subscriber was already installed; now uses `try_init`
- `network::rate_from_stoichiometry` silently truncated ragged rows; `reaction_dynamics` caps `steps` at 10⁷ and validates `tol`
- Many documentation errors corrected (units, return semantics, "sparse" vs dense, DOPRI45 "for stiff systems", steady-state condition, 1 atm vs 1 bar)

### Performance
Criterion medians, v1.2 refresh (criterion 0.8, pre-repair) → v1.2 post-repair; full table in `benchmarks.md`.
- thermochem/enthalpy_change_cp: 473 ns → 24.0 ns (−95%) — analytic Shomate enthalpy replaces Simpson integration
- thermochem/adiabatic_flame_ch4_air: 2.76 µs → 301 ns (−89%) — segments resolved once, analytic H(T)
- spectroscopy/isotope_pattern_glucose: 12.7 µs → 2.80 µs (−78%) — pruning + buffer reuse
- stoichiometry/oxidation_states_kmno4: 47.0 ns → 22.4 ns (−52%)
- nuclear/semi_empirical_fe56: 25.1 ns → 18.5 ns (−26%); potential/dihedral_angle: 27.2 ns → 20.4 ns (−25%)
- Cost of correctness (slower, by design):
  - solution/weak_acid_ph: 19.6 → 34.9 ns — full charge-balance cubic instead of the quadratic approximation
  - gas/peng_robinson_z_n2_1mpa: 33.2 → 61.4 ns — analytic cubic roots instead of Newton from Z = 1 (which could return an unphysical root)
  - electrochemistry/lookup_half_reaction: 29.7 → 57.1 ns — linear scan; sorting the table moved Cu²⁺/Cu later
  - Scalar formulas gain ~0.5–1.8 ns from NaN/±∞ guards (e.g. ideal_gas_pressure 1.25 → 2.18 ns, arrhenius_rate 5.80 → 6.78 ns, michaelis_menten 1.16 → 2.92 ns)
- The criterion 0.5 → 0.8 harness change alone shifted several baselines (e.g. michaelis_menten 3.04 → 1.16 ns with no code change); no claims are made for that step

### Behaviour changes (no signature change)
- Corrected values: everything under **Fixed — physics and data** above
- New `Err` cases for input that used to return `Ok`: NaN/±∞ anywhere; `balance_equation` non-integer, ambiguous or mixed-sign input; `assign_oxidation_states` non-integer average states, unknown Z, zero count; Shomate functions outside the data range; `bohr_transition_energy` with n_upper ≤ n_lower; `isotope_pattern` unknown Z; `lennard_jones` ε ≤ 0; negative vdW a/b, Henry x > 1, `phase_rule(_, 0)`, `percent_yield` negative actual
- `decay_chain_two_step` now returns the analytic limit for equal rates (was `Err`)
- `nth_order_half_life` accepts any order ≠ 1 (was order > 1 only)
- Sentinels: `bridge` functions return `NaN` (was `0.0`) for invalid input; Antoine functions return `NaN` at the pole; `semi_empirical_binding_energy` returns `NaN` for z > a or a = 0 (was a panic/wrap)
- `lookup_shomate("N2(g)")` / `("O2(g)")` return the 100 K segment; `SHOMATE_DATA` repeats formulas, one entry per segment
- `predict_geometry(1, n)` returns `Linear` (was `None`); `Molecule::formula` omits zero-count atoms; `find_peaks` threshold is relative to the data's min..max
- soorat `MolecularStructure` bonds form a star around the least electronegative non-H atom (was a chain in input order); radii derive from covalent radii

## [1.1.1] - 2026-03-29

### Fixed
- **license** — corrected SPDX identifier from `GPL-3.0` to `GPL-3.0-only` across Cargo.toml, deny.toml, README, and CLAUDE.md

## [1.1.0]

### Added
- **bridge** — cross-crate primitive-value bridges for ushma (reaction heat release, Arrhenius rate factor, equilibrium to Gibbs), tanmatra (atomic number to valence, ionization to bond energy), bijli (electrochemical current density, ionic conductivity)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `MolecularStructure` (ball-and-stick from `Molecule`), `ReactionNetworkVisualization` (species/reaction graph), `SpectrumVisualization` (x/y data with peak detection), `PhaseDiagramVisualization` (boundary curves with critical/triple points)

### Updated
- hisab 1.1.0 -> 1.3.0, zerocopy 0.8.47 -> 0.8.48

## [1.0.0] - 2026-03-25

Stable release. API frozen. 21 modules, 118 elements, 444 tests.

### Added

- **stoichiometry** module — chemical equation balancing (via hisab gaussian_elimination), oxidation state assignment, limiting reagent, theoretical yield, percent yield, percent composition, empirical formula
- **phase** module — Raoult's law, Henry's law, colligative properties (bp elevation, fp depression, osmotic pressure), Clausius-Clapeyron, Antoine equation with 8 solvents, Gibbs phase rule
- **nuclear** module — radioactive decay law, activity, mass defect, binding energy per nucleon, semi-empirical mass formula (Weizsacker), Q-value, 18 isotopes (H-3 through Np-237)
- **inorganic** module — crystal field theory (CFSE oct/tet), spectrochemical series (11 ligands), spin-only magnetic moment, Born-Landé lattice energy, 5 Madelung constants, 30 Shannon ionic radii
- **potential** module — Lennard-Jones 6-12 potential, Morse potential, Coulomb potential and force
- **reaction_dynamics** module — ODE-powered kinetics simulation (1st/2nd order, consecutive A→B→C, reversible A⇌B, Michaelis-Menten substrate depletion, adaptive DOPRI45) via hisab rk4_trajectory
- **fitting** module — Arrhenius parameter fitting, polynomial calibration curves, Beer-Lambert extinction coefficient fitting via hisab least_squares_poly
- **electrochemistry** — expanded to 49 half-reactions (added MnO₄⁻, Cr₂O₇²⁻, O₂/H₂O, Fe³⁺/Fe²⁺, Ce⁴⁺, H₂O₂, PbO₂, and 23 more)
- **thermochem** — expanded to 71 substances (added 10 elements, 23 inorganic compounds, 7 organic compounds); Shomate Cp(T) expanded to 24 datasets
- **spectroscopy** — wavenumber conversions, reduced mass, quantum harmonic oscillator energy, vibrational frequency from force constant
- **organic** — bond energies expanded to 48 types (added P-H/O/Cl, Si-H/O/C/Cl, B-H/O/N/F, S-S, S=O)
- **kinetics** — zero-order concentration and half-life, nth-order half-life, Michaelis-Menten, Lineweaver-Burk, collision theory, Eyring equation
- **gas** — Peng-Robinson equation of state (pressure, Z-factor), fugacity coefficient, fugacity, 5 substance parameter sets
- **solution** — weak acid/base pH (Newton-Raphson), ionic strength, Debye-Hückel (limiting, extended), Davies equation, Kohlrausch molar conductivity, 7 limiting ionic conductivities
- **element** — full periodic table: 118 elements (H–Og) with Lanthanide and Actinide categories
- **hisab integration** — newton_raphson (8.3x speedup for pH), rk4_trajectory, dopri45, gaussian_elimination, least_squares_poly, integral_simpson, EPSILON_F64

### Changed

- All public functions that can fail now return `Result<T, KimiyaError>` instead of silent sentinel values
- `KimiyaError` now derives `Clone`, `PartialEq`, `Eq`; added `InvalidInput` variant
- `Element` no longer derives `Deserialize` — `&'static str` fields cannot be deserialized
- `lookup_by_number` is O(1) via direct array indexing (was O(n) linear scan)
- Temperature validation added to `ideal_gas_pressure`, `ideal_gas_volume`, `gibbs_free_energy`
- Nuclear constants upgraded to CODATA 2018 full precision
- `hisab` dependency updated from 0.24 to 1.0

### Removed

- `Bond` enum from molecule module (unused dead code)
- `hess_law_enthalpy` from reaction module (duplicate of `thermo::hess_law`)

## [0.1.0] - 2026-03-24

Initial scaffold with real chemistry implementations.

### Added

- **element** — 36 elements (H–Kr), `ElementCategory`, lookup by number/symbol, Avogadro constant, gas constant
- **molecule** — `Molecule` struct, molecular weight, formula string, presets (water, CO₂, methane, glucose)
- **reaction** — Gibbs free energy, equilibrium constant, spontaneity check
- **kinetics** — Arrhenius rate constant, half-life, first/second-order concentration decay
- **gas** — ideal gas law (P and V), Van der Waals equation, partial pressure, gas density, STP constants
- **solution** — molarity, molality, dilution, pH/pOH, Henderson-Hasselbalch equation
- **thermo** — heat transfer (q = mcΔT), Hess's law, formation enthalpy, specific heat from calorimetry
- **error** — `KimiyaError` with `#[non_exhaustive]` and `thiserror` integration

[2.0.0]: https://github.com/MacCracken/kimiya/releases/tag/2.0.0
[1.2.0]: https://github.com/MacCracken/kimiya/releases/tag/1.2.0
[1.0.0]: https://github.com/MacCracken/kimiya/releases/tag/v1.0.0
