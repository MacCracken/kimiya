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
| M3 Solutions & electrochemistry | ✅ **ported 2026-11-09** (six kinetics class-C functions measured, not gated, until M3c). M3a: `src/solution.cyr` (concentration units, pH, Henderson–Hasselbalch, the weak acid / base charge-balance Newton transcribed as `_km_charge_balance_root`, ionic strength over ion records, Debye–Hückel, Davies, Kohlrausch) with its ten constants generated into `src/constants.cyr`; generated `src/electrochemistry_data.cyr` (49 STANDARD_POTENTIALS rows, the first UTF-8 string blob) + `src/electrochemistry.cyr` (handles, Nernst, cell potentials, ΔG, Faraday's laws); `_km_log10_ln` in `src/mathx.cyr` (glibc's log10 composition over the correctly rounded `_km_ln`) and `_km_ln` for Nernst (ADR 0005 amendment M3a). M3b: `src/kinetics.cyr` (Arrhenius, the zero / first / second / nth-order laws and half-lives, steady state and pre-equilibrium, Michaelis–Menten, Lineweaver–Burk into the caller's two-f64 buffer, collision theory, Eyring), on `f64_exp` / `_km_expm1` / `ganita_f64_pow`. Golden `solution` 28,566 + `electrochemistry` 18,574 + `kinetics` 30,146 rows; **119,044 ok, 0 failed, 22,205 measured, 10 pending, 441 accepted** on x86_64 ELF, x86_64 PE and aarch64 (class B 0 ulp; M3a class C K = 2–4; the six kinetics class-C functions measured, not gated: subnormal exp / pow results scaled by a prefactor put S at 10¹⁴ … 4.4 × 10¹⁸, left for ADR 0006 under the user's functional-parity decision, port plan §6.4); 131 Rust tests ported; surface parity 203/203; 7 benches vs Rust |
| M3c Correct rounding (ADR 0006) | ☐ next: correctly rounded primitives everywhere; class-C rows compared against the oracle's formula with correctly rounded results |
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
  Debye–Hückel A / B, the seven Λ°), `src/mathx.cyr` (incl. the correctly rounded `_km_acos` / `_km_cos` / `_km_ln` /
  `_km_pow_2_5` and `_km_log10_ln`, glibc's log10 composition over `_km_ln`), `src/validate.cyr`, `src/element_data.cyr`
  (generated: 118 elements, 32 property rows, `KmElementCategory`), `src/element.cyr`, `src/molecule.cyr`,
  `src/reaction.cyr`, `src/thermo.cyr`, `src/thermochem_data.cyr` (generated: `KmShomate`, 76 THERMOCHEM_DATA rows
  with formula / name blobs, 34 SHOMATE_DATA segments), `src/thermochem.cyr`, `src/phase_data.cyr` (generated:
  `KmAntoineCoeffs`, 8 Antoine sets, the
  mmHg → Pa factor, K_b/K_f), `src/phase.cyr`, `src/gas_data.cyr` (generated: `KmPengRobinsonParams`, STP, the 5
  `KM_PR_*` presets), `src/gas.cyr`, `src/solution.cyr`, `src/electrochemistry_data.cyr` (generated: the 49
  STANDARD_POTENTIALS rows — UTF-8 equation and couple blobs, electrons, E° — and `KM_STANDARD_POTENTIALS_COUNT`),
  `src/electrochemistry.cyr`, `src/kinetics.cyr`; plus the `src/main.cyr` stub. Generated files are written by
  `scripts/oracle/` and never edited. `dist/kimiya.deps`: math, string, vec, str, ganita, tagged; the bundle also
  includes `lib/hisab.cyr` (thermo's and thermochem's `EPSILON_F64`), so a consumer needs kimiya's `[deps.hisab]`.

## Tests

- `tests/parity.tcyr` (+ `tests/parity_core.cyr`, `tests/ulp.cyr`): the golden-corpus harness — 12 modules
  (constants, electrochemistry, element, gas, kinetics, math, molecule, phase, reaction, solution, thermo, thermochem),
  119,044 ok, 0 failed, 22,205 measured (the six kinetics class-C functions, no bound until ADR 0006), 10 pending (the
  M4 inorganic constants), 441 accepted (ADR 0004's row and 440 ADR 0002 #4 mapped Antoine rows), on x86_64 ELF, x86_64
  PE and aarch64
- `tests/golden/`: `constants.tsv` (45), `electrochemistry.tsv` (18,574), `element.tsv` (975), `gas.tsv` (19,717),
  `kinetics.tsv` (30,146), `math.tsv` (26,147), `molecule.tsv` (2,472), `phase.tsv` (3,276), `reaction.tsv` (1,188),
  `solution.tsv` (28,566), `thermo.tsv` (1,107), `thermochem.tsv` (9,487), `modules.txt` — generated by `scripts/oracle/`
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
  `kinetics/eyring_rate`, `kinetics/michaelis_menten`
- `tests/kimiya.fcyr`: scaffold fuzz harness

## Dependencies

- stdlib — syscalls, string, alloc, str, fmt, io, vec, args, math, ganita, overflow, tagged, result, fnptr
- hisab 3.4.1 (git tag) → sakshi 2.5.10
- test stdlib — assert, bench

## Consumers

- **Rust (crates.io, as of 2026-11-08), all on the 1.x crate:** jivanu, khanij (Arrhenius, pH, first-order
  kinetics, Born–Landé, ionic radii, thermochem fields, `Molecule::new`), jnana (`ELEMENTS`, three constants),
  kiran (re-exports whole modules), tara (declared, unused).
- **Cyrius:** none yet. tanmatra 2.0.0 deferred its kimiya bridges until kimiya ports.

## Next

**M3c, before M4** (structure, nuclear, spectroscopy, with the 10 pending inorganic constants): the step the
user decided on 2026-11-09 (functional parity, port plan §6.4): **ADR 0006**, correctly rounded primitives in place of
glibc-matching ones, with the class-C rows compared against the oracle's formula evaluated with correctly rounded
results. It takes in M3b's six measured kinetics functions (subnormal exp / pow results scaled into the normal range,
port-parity.md §8; a correctly rounded pow must round exact half-way cases to even, where glibc 2.44's does not), M3a's
`henderson_hasselbalch` / `_km_log10_ln` and ADR 0005, and the M2 functions with the same latent issue outside their
grids (`phase::clausius_clapeyron`, `thermochem::vant_hoff_k`, `gas::fugacity_pr`: a subnormal exp result times a
prefactor, up to 12–15 ulp from the oracle on hand-picked inputs, past their K of 6–10). Also a two-phase (Ziv) correctly
rounded ln to recover `nernst_potential`'s speed: `_km_ln` is 1.9 µs a call (`nernst_potential` is 302× Rust; ADR 0005
amendment M3a). Follow-ups outside this repo: a hisab feature request (ADR 0003, now also noting that `num_newton` /
`num_bisection` take no callback context); a note to tanmatra about its `parity_core.cyr` per-function table sizing (768
bytes for a 192-slot limit; overflows past 96 functions); and a cyrius note that `distlib`'s sidecar verify misses stdlib
leaves a bundle only calls from unreached functions (kimiya works around it with explicit `include`s, port plan §6.4).
From M2a: a note to ganita that its acos differs from glibc by up to 2 ulp and gives different bits on aarch64 than on
x86_64 (port-parity.md §5).
