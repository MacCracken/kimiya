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
| M3 Solutions & electrochemistry | ☐ |
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
  (generated, 15 `KM_` constants: the 14 physical ones and `KM_WATER_SPECIFIC_HEAT`), `src/mathx.cyr` (incl. the
  correctly rounded `_km_acos` / `_km_cos` / `_km_ln` / `_km_pow_2_5`), `src/validate.cyr`, `src/element_data.cyr`
  (generated: 118 elements, 32 property rows, `KmElementCategory`), `src/element.cyr`, `src/molecule.cyr`,
  `src/reaction.cyr`, `src/thermo.cyr`, `src/thermochem_data.cyr` (generated: `KmShomate`, 76 THERMOCHEM_DATA rows
  with formula / name blobs, 34 SHOMATE_DATA segments), `src/thermochem.cyr`, `src/phase_data.cyr` (generated:
  `KmAntoineCoeffs`, 8 Antoine sets, the
  mmHg → Pa factor, K_b/K_f), `src/phase.cyr`, `src/gas_data.cyr` (generated: `KmPengRobinsonParams`, STP, the 5
  `KM_PR_*` presets), `src/gas.cyr`; plus the `src/main.cyr` stub. Generated files are written by
  `scripts/oracle/` and never edited. `dist/kimiya.deps`: math, string, vec, str, ganita, tagged; the bundle also
  includes `lib/hisab.cyr` (thermo's and thermochem's `EPSILON_F64`), so a consumer needs kimiya's `[deps.hisab]`.

## Tests

- `tests/parity.tcyr` (+ `tests/parity_core.cyr`, `tests/ulp.cyr`): the golden-corpus harness — 9 modules
  (constants, element, gas, math, molecule, phase, reaction, thermo, thermochem), 63,953 ok, 0 failed, 0 measured,
  20 pending (module-data constants, M3–M4), 441 accepted (ADR 0004's row and 440 ADR 0002 #4 mapped Antoine
  rows), on x86_64 ELF, x86_64 PE and aarch64
- `tests/golden/`: `constants.tsv` (45), `element.tsv` (975), `gas.tsv` (19,717), `math.tsv` (26,147),
  `molecule.tsv` (2,472), `phase.tsv` (3,276), `reaction.tsv` (1,188), `thermo.tsv` (1,107), `thermochem.tsv`
  (9,487), `modules.txt` — generated by `scripts/oracle/`
- Reference tests, the Rust tests ported one for one: `tests/error.tcyr` (2 Rust + 1 Cyrius-only group,
  24 assertions), `tests/validate.tcyr` (3 + 1, 25), `tests/element.tcyr` (22 + 1, 96), `tests/molecule.tcyr`
  (11 + 1, 39; `serde_roundtrip` doesn't port), `tests/reaction.tcyr` (10 + 1, 29), `tests/thermo.tcyr`
  (10 + 1, 24), `tests/thermochem.tcyr` (62 + 6, 400; 4 tests assert J/mol for ADR 0002 #1), `tests/phase.tcyr`
  (31 + 2, 117), `tests/gas.tcyr` (39 + 2, 159)
- `tests/kimiya.tcyr`: smoke test of `dist/kimiya.cyr` as a consumer includes it (16 assertions)
- `tests/kimiya.bcyr`: `element/lookup_by_symbol`, `element/lookup_by_number`, `molecule/molecular_weight_glucose`,
  `reaction/gibbs_free_energy`, `thermo/heat_transfer`, `gas/ideal_gas_pressure`, `gas/peng_robinson_z_n2_1mpa`,
  `thermochem/reaction_enthalpy_ch4`, `thermochem/vant_hoff_k`, `thermochem/enthalpy_change_cp`,
  `thermochem/shomate_cp`, `thermochem/adiabatic_flame_ch4_air`
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

**M3**: solution, electrochemistry (+ generated `electrochemistry_data`: the 49 potentials), kinetics — the weak
acid/base charge-balance Newton (kimiya's own loop, transcribed), ionic strength, Debye–Hückel, Davies, Kohlrausch,
Nernst, Arrhenius / Eyring / Michaelis–Menten / nth-order (exp, `exp_m1` → class C), and their module-data constants
(K_w, the Debye–Hückel A / B, Λ°: 10 of the 20 pending `constants` rows). Follow-ups outside this repo: a hisab feature
request (ADR 0003, now also noting that `num_newton` / `num_bisection` take no callback context); a note to tanmatra
about its `parity_core.cyr` per-function table sizing (768 bytes for a 192-slot limit; overflows past 96 functions);
and a cyrius note that `distlib`'s sidecar verify misses stdlib leaves a bundle only calls from unreached functions
(kimiya works around it with explicit `include`s, port plan §6.4). From M2a: a note to ganita that its acos differs
from glibc by up to 2 ulp and gives different bits on aarch64 than on x86_64 (port-parity.md §5).
