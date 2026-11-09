# kimiya — Current State

> Refreshed every release and at every port milestone. `CLAUDE.md` is
> preferences/process/procedures (durable); this file is **state** (volatile).

## Version

**1.2.0** (`VERSION`) — the last Rust release, tag `1.2.0` (`077e2e2`), frozen as the oracle in
`rust-old/`. The Cyrius port is in progress and ships as **2.0.0**; `VERSION` moves to 2.0.0 at M7.
`cyrius.cyml` reads the version through `version = "${file:VERSION}"`.

## Port progress

Plan of record: [`kimiya-port-plan.md`](kimiya-port-plan.md).

| Milestone | Status |
|---|---|
| M0 Initialize | ✅ **done 2026-11-09.** Steps 0–4: oracle frozen green at tag `1.2.0` (662 + 24 + 2 Rust tests, rustc 1.99.1, glibc 2.44; `rust-old/README.md` restored from the tag for its doctest); `cyrius port` scaffold; CLAUDE.md merged, roadmap restored, Cyrius CI extended, Rust CI/scripts/overview under `rust-old/`, `Cargo.lock` tracked; `cyrius.cyml` with hisab 3.4.1 (`d6703c5`); math semantics probed. **Step 5**: `scripts/oracle/` writes `tests/golden/` (constants 45 + math 18,647 rows, 2.0 MB), deterministic, platform-pinned, `--check`ed, independently validated (`validate.py`, 60-digit reference). **Step 6**: `src/constants.cyr` (14 `KM_` constants) generated and compiler-probed. **Step 7**: `tests/parity.tcyr` + `parity_core.cyr` + `ulp.cyr`. **Step 8**: [`port-parity.md`](port-parity.md) — K = 2 ulp for the math primitives on x86_64 ELF, x86_64 PE and aarch64; 13,340 ok, 0 failed, 5,351 pending, 1 accepted on all three. **Step 9**: ADR 0003 (transcribe the hisab 1.4.0 routines hisab 3.4.1 lacks); ADR 0004 (ganita `pow(-1, DBL_MAX)`) |
| M1 Foundation | ✅ **done 2026-11-09.** `src/error.cyr` (`KmError`), `src/mathx.cyr` (the `_km_` math helpers; `_km_log10` / `_km_expm1` from Rust libm 0.2.16), `src/validate.cyr`, generated `src/element_data.cyr` + `src/element.cyr` (handles, 0 for none), `src/molecule.cyr` (`KmMolecule`); golden `element` 975 + `molecule` 2,472 rows; **22,107 ok, 0 failed, 31 pending, 1 accepted** on x86_64 ELF, x86_64 PE and aarch64 (math 0 pending; log10 / exp_m1 K = 2); 38 Rust tests ported; `[lib]` + `dist/kimiya.cyr`; surface parity 55/55; benches vs Rust in [`benchmarks-rust-v-cyrius.md`](../benchmarks-rust-v-cyrius.md) |
| M2 Thermochemistry & gases | ✅ **done 2026-11-09.** M2a: `src/reaction.cyr`, `src/thermo.cyr`, generated `src/phase_data.cyr` + `src/phase.cyr` (Antoine → `Result`, ADR 0002 #4), generated `src/gas_data.cyr` + `src/gas.cyr` (`KmPengRobinsonParams`, the analytic PR cubic); correctly rounded acos/cos/ln/x^2.5 in `src/mathx.cyr` (ADR 0005). M2b: generated `src/thermochem_data.cyr` (`KmShomate`, 76 THERMOCHEM_DATA rows, 34 Shomate segments) + `src/thermochem.cyr` (reaction sums in J/mol, ADR 0002 #1; van 't Hoff; Shomate Cp and enthalpy; the flame temperature with hisab 1.4.0's Newton + bisection transcribed, ADR 0003 amendment). Golden `reaction` 1,188 + `thermo` 1,107 + `thermochem` 9,487 + `phase` 3,276 + `gas` 19,717 rows (math +7,500); **63,953 ok, 0 failed, 0 measured, 20 pending, 441 accepted** on x86_64 ELF, x86_64 PE and aarch64 (class B 0 ulp, the flame rows too; class C K = 2–10); 152 Rust tests ported; surface parity 140/140; 9 benches vs Rust |
| M3 Solutions & electrochemistry | ✅ **done 2026-11-09.** M3a: `src/solution.cyr` (concentration units, pH, Henderson–Hasselbalch, the weak acid / base charge-balance Newton transcribed as `_km_charge_balance_root`, ionic strength over ion records, Debye–Hückel, Davies, Kohlrausch) with its ten constants generated into `src/constants.cyr`; generated `src/electrochemistry_data.cyr` (49 STANDARD_POTENTIALS rows, the first UTF-8 string blob) + `src/electrochemistry.cyr` (handles, Nernst, cell potentials, ΔG, Faraday's laws). M3b: `src/kinetics.cyr` (Arrhenius, the zero / first / second / nth-order laws and half-lives, steady state and pre-equilibrium, Michaelis–Menten, Lineweaver–Burk into the caller's two-f64 buffer, collision theory, Eyring). Golden `solution` 28,566 + `electrochemistry` 18,574 + `kinetics` 30,146 rows; at M3b 119,044 ok, 0 failed, 22,205 measured (the six kinetics class-C functions, gated since M3c), 10 pending, 441 accepted on all three targets; 131 Rust tests ported; surface parity 203/203; 7 benches vs Rust. Its primitives moved to the correctly rounded ones in M3c |
| M3c Correct rounding (ADR 0006) | ✅ **done 2026-11-09.** M3c-1: the corrected oracle — `scripts/oracle/mpfr_ref.py` (MPFR 4.2 through ctypes), `validate.py --write-corrections` → `tests/golden/<module>.cr.tsv` (665 rows), the harness compares corrected rows against them. M3c-2: correctly rounded `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`, `_km_cr_pow` in `src/mathx.cyr` (two-phase, Ziv), tables generated into `src/mathx_data.cyr` by `scripts/oracle/mathx_tables.py`, verified by `scripts/oracle/cr_verify.py` (≥ 1,000,000 arguments per primitive against MPFR, 0 mismatches, the same bits on every target) and `tests/mathx.tcyr`. M3c-3: every module on them in Rust's operation order (`_km_acos` / `_km_cos` kept; x^2.5 is `_km_cr_pow`); `_km_log10`, `_km_expm1`, `_km_log10_ln`, the `_km_ln` wrapper, `_km_pow_2_5` and `_km_dd_div` removed; every class-C function at bound 0; glibc-matching is a failure; 2,486 hard rows (subnormal exp × prefactor) added to `reaction`, `phase`, `thermochem`, `gas`. **143,071 ok, 0 failed, 0 measured, 10 pending, 1,105 accepted (ADR 6: 665; ADR 0002 #4: 440), 0 glibc-matching** on x86_64 ELF, x86_64 PE and aarch64, every class-C row 0 ulp (sin, atan2 and the 247 cos rows past \|x\| = 8 on the stdlib at 2 ulp); ADR 0004's row a plain ok; 15 test files, 0 warnings; surface parity 203/203 ([`port-parity.md`](port-parity.md) §9–§10) |
| M4 Structure, nuclear & spectroscopy | ☐ |
| M5 Numerics | ☐ |
| M7 Release 2.0.0 | ☐ |

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml` `[package].cyrius`)
- **ganita**: 1.2.15, from the 6.7.18 stdlib bundle (`lib/ganita.cyr`)
- **hisab**: 3.4.1 (`[deps.hisab]`, commit-pinned in `cyrius.lock`), with its sakshi 2.5.10
- **Rust oracle**: kimiya 1.2.0 at `rust-old/`, green on rustc 1.99.1; its transcendental results come
  from the host **glibc 2.44** (Rust std)

## Source

- **Rust reference:** `rust-old/` (16,581 lines incl. tests), frozen; do not edit. `rust-old/docs/`,
  `rust-old/.github/` and `rust-old/scripts/` keep the 1.2.0 architecture overview, CI and cargo scripts;
  `rust-old/benchmarks.md` is the Rust benchmark baseline.
- **Cyrius port:** `[lib]` (bundle `dist/kimiya.cyr`, `cyrius distlib`) = `src/error.cyr`, `src/constants.cyr`
  (generated, 25 `KM_` constants: the 14 physical ones, `KM_WATER_SPECIFIC_HEAT` and the ten solution constants — K_w,
  Debye–Hückel A / B, the seven Λ°), `src/mathx_data.cyr` (generated by `scripts/oracle/mathx_tables.py` from MPFR: the
  correctly rounded primitives' tables and constants), `src/mathx.cyr` (the Rust-semantics helpers and the correctly
  rounded primitives, ADR 0006: `_km_acos` / `_km_cos` in double-double (ADR 0005), and the two-phase `_km_cr_exp`,
  `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`, `_km_cr_pow`, each with `_acc` / `_fast` test entry
  points), `src/validate.cyr`, `src/element_data.cyr`
  (generated: 118 elements, 32 property rows, `KmElementCategory`), `src/element.cyr`, `src/molecule.cyr`,
  `src/reaction.cyr`, `src/thermo.cyr`, `src/thermochem_data.cyr` (generated: `KmShomate`, 76 THERMOCHEM_DATA rows
  with formula / name blobs, 34 SHOMATE_DATA segments), `src/thermochem.cyr`, `src/phase_data.cyr` (generated:
  `KmAntoineCoeffs`, 8 Antoine sets, the
  mmHg → Pa factor, K_b/K_f), `src/phase.cyr`, `src/gas_data.cyr` (generated: `KmPengRobinsonParams`, STP, the 5
  `KM_PR_*` presets), `src/gas.cyr`, `src/solution.cyr`, `src/electrochemistry_data.cyr` (generated: the 49
  STANDARD_POTENTIALS rows — UTF-8 equation and couple blobs, electrons, E° — and `KM_STANDARD_POTENTIALS_COUNT`),
  `src/electrochemistry.cyr`, `src/kinetics.cyr`; plus the `src/main.cyr` stub. Generated files are written by
  `scripts/oracle/` and never edited. Every transcendental the modules call is correctly rounded; only `gas` still
  includes `lib/ganita.cyr` (cbrt). `dist/kimiya.deps`: math, string, vec, str, ganita, tagged; the bundle also
  includes `lib/hisab.cyr` (thermo's and thermochem's `EPSILON_F64`), so a consumer needs kimiya's `[deps.hisab]`.

## Tests

- `tests/parity.tcyr` (+ `tests/parity_core.cyr`, `tests/ulp.cyr`): the golden-corpus harness — 12 modules
  (constants, electrochemistry, element, gas, kinetics, math, molecule, phase, reaction, solution, thermo, thermochem),
  143,071 ok, 0 failed, 0 measured, 10 pending (the M4 inorganic constants), 1,105 accepted (665 under ADR 6, against the
  corrected oracle, and the 440 ADR 0002 #4 mapped Antoine rows), 0 glibc-matching, on x86_64 ELF, x86_64 PE and
  aarch64. Every class-C function at bound 0 (ADR 0006); `f64::sin`, `f64::atan2` and the cos rows past |x| = 8 at the
  stdlib's 2 until M4
- `tests/golden/`: `constants.tsv` (45), `electrochemistry.tsv` (18,574), `element.tsv` (975), `gas.tsv` (21,093),
  `kinetics.tsv` (30,146), `math.tsv` (26,147), `molecule.tsv` (2,472), `phase.tsv` (3,646), `reaction.tsv` (1,558),
  `solution.tsv` (28,566), `thermo.tsv` (1,107), `thermochem.tsv` (9,857), `modules.txt` — generated by
  `scripts/oracle/`; the corrected oracle (ADR 0006) `gas.cr.tsv` (11), `kinetics.cr.tsv` (323), `math.cr.tsv` (96),
  `phase.cr.tsv` (4), `reaction.cr.tsv` (2), `solution.cr.tsv` (229) — written by `validate.py --write-corrections`
- `tests/mathx.tcyr`: the correctly rounded primitives (M3c-2) — special values, exact and tie cases (10²³, the
  subnormal ties), thresholds, subnormal results, the known hardest-to-round arguments, the fast path's fallbacks, and
  seeded fast / accurate agreement, every expected value MPFR's (11 groups, 518 assertions); `scripts/oracle/cr_verify.py` checks
  ≥ 1,000,000 arguments per primitive on every target
- Reference tests, the Rust tests ported one for one: `tests/error.tcyr` (2 Rust + 1 Cyrius-only group,
  24 assertions), `tests/validate.tcyr` (3 + 1, 25), `tests/element.tcyr` (22 + 1, 96), `tests/molecule.tcyr`
  (11 + 1, 39; `serde_roundtrip` doesn't port), `tests/reaction.tcyr` (10 + 1, 29), `tests/thermo.tcyr`
  (10 + 1, 24), `tests/thermochem.tcyr` (62 + 6, 400; 4 tests assert J/mol for ADR 0002 #1), `tests/phase.tcyr`
  (31 + 2, 117), `tests/gas.tcyr` (39 + 2, 159), `tests/solution.tcyr` (49 + 5, 135), `tests/electrochemistry.tcyr`
  (35 + 3, 1,394), `tests/kinetics.tcyr` (47 + 4, 138; 2 tests adapted for the Lineweaver–Burk buffer)
- `tests/kimiya.tcyr`: smoke test of `dist/kimiya.cyr` as a consumer includes it (24 assertions)
- `tests/kimiya.bcyr`: `element/lookup_by_symbol`, `element/lookup_by_number`, `molecule/molecular_weight_glucose`,
  `reaction/gibbs_free_energy`, `thermo/heat_transfer`, `gas/ideal_gas_pressure`, `gas/peng_robinson_z_n2_1mpa`,
  `thermochem/reaction_enthalpy_ch4`, `thermochem/vant_hoff_k`, `thermochem/enthalpy_change_cp`,
  `thermochem/shomate_cp`, `thermochem/adiabatic_flame_ch4_air`, `solution/ph_from_h`, `solution/weak_acid_ph`,
  `electrochemistry/lookup_half_reaction`, `electrochemistry/nernst_potential`, `kinetics/arrhenius_rate`,
  `kinetics/eyring_rate`, `kinetics/michaelis_menten`, and (Cyrius-only) `mathx/cr_exp`, `cr_expm1`, `cr_exp10`,
  `cr_ln`, `cr_log10`, `cr_pow` beside `ganita_f64_pow` (two arguments each)
- `tests/kimiya.fcyr`: scaffold fuzz harness

## Dependencies

- stdlib — syscalls, string, alloc, str, fmt, io, vec, args, math, ganita, overflow, tagged, result, fnptr
- hisab 3.4.1 (git tag) → sakshi 2.5.10
- test stdlib — assert, bench
- tooling only (not the shipped library): the system MPFR 4.2 (`libmpfr.so.6`, through ctypes) for
  `scripts/oracle/mpfr_ref.py`, `validate.py`, `mathx_tables.py` and `cr_verify.py`; wine and qemu-aarch64 for the PE and
  aarch64 test runs

## Consumers

- **Rust (crates.io, as of 2026-11-08), all on the 1.x crate:** jivanu, khanij (Arrhenius, pH, first-order
  kinetics, Born–Landé, ionic radii, thermochem fields, `Molecule::new`), jnana (`ELEMENTS`, three constants),
  kiran (re-exports whole modules), tara (declared, unused).
- **Cyrius:** none yet. tanmatra 2.0.0 deferred its kimiya bridges until kimiya ports.

## Next

**M4** (structure, nuclear, spectroscopy, with the 10 pending inorganic constants). It brings correctly rounded **sin**
and **atan2** with their first callers (the dihedral, spectroscopy; ADR 0006 item 1), built and verified as M3c-2's were
(Ziv, generated tables, `cr_verify.py`), so that the `f64::sin` / `f64::atan2` rows and every M4 class-C function are
gated at 0 against the corrected oracle; `_km_cos` would gain its range past |x| = 8 then if an M4 caller needs it.

Follow-ups:
- **pow's accurate path has no third phase.** It is proven correct only for results not within 2⁻¹³⁰ (relative) of a
  midpoint; CORE-MATH's pow has a 256-bit third phase that reports rather than misround (src/mathx.cyr's pow header;
  every published hard case and every verified pair are well clear). A random pair is that close with probability
  about 2⁻⁷⁶.
- **log10 got about 10 ns slower** (to about 77 ns) when pow's absolute ln bound was shared with it in M3c-2b; a separate,
  tighter bound for log10 could win it back.
- **The Cyrius issue** `cyrius/docs/development/issues/2026-11-09-exp-polyfill-double-rounds-subnormal-results.md`
  (filed): the aarch64 exp polyfill rounds subnormal results twice; until a Cyrius library ships correctly rounded
  primitives, kimiya owns them (ADR 0006).
- A hisab feature request (ADR 0003, now also noting that `num_newton` / `num_bisection` take no callback context).
- A note to tanmatra about its `parity_core.cyr` per-function table sizing (768 bytes for a 192-slot limit; overflows past
  96 functions).
- A cyrius note that `distlib`'s sidecar verify misses stdlib leaves a bundle only calls from unreached functions (kimiya
  works around it with explicit `include`s, port plan §6.4).
- From M2a: a note to ganita that its acos differs from glibc by up to 2 ulp and gives different bits on aarch64 than on
  x86_64 (port-parity.md §5); kimiya no longer calls it.
