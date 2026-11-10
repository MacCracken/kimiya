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
| M4 Structure, nuclear & spectroscopy | ◐ **in progress.** **M4b ✅ done 2026-11-09** — generated `src/organic_data.cyr` (`KmBondType`, `KmGeometry`, `KmFunctionalGroup` with their property columns: bond energies, ideal angles, the UTF-8 pattern blob, the polar / donor flags and ADR 0002 #6's new acceptor column: IUPAC 2011's criterion, weak acceptors included, every group but the alkane, one source per family) + `src/organic.cyr` (`bond_type_energy_kj`, `enthalpy_from_bonds` over bond records, `predict_geometry` → `Option`, `geometry_ideal_bond_angle`, the `functional_group_*` properties; ADR 0002 #6: `KM_GEOMETRY_T_SHAPED`, `functional_group_is_hydrogen_bond_donor`, `…_acceptor`); generated `src/inorganic_data.cyr` (`KmLigand`, the private octahedral occupancy and unpaired-electron tables read through the API, the 10 Madelung / Born constants, the 30 Shannon IONIC_RADII rows) + `src/inorganic.cyr` (CFSE exact by ADR 0002 #12, (3·e_g − 2·t₂g)/5 in one correctly rounded division, Δ_tet, the spectrochemical series, spin-only moments, Born–Landé, Born–Haber, `lookup_ionic_radius` handles). Golden `organic` 9,104 + `inorganic` 15,107 rows and the last 10 `constants` rows: **167,279 ok, 0 failed, 0 measured, 0 pending, 1,118 accepted (ADR 6: 665; ADR 2: 440 Antoine + 13 exact-CFSE rows, by name), 0 glibc-matching** on x86_64 ELF, x86_64 PE and aarch64, every M4b row 0 ulp from its reference; 51 Rust tests ported (+ 2 integration); surface parity 323/323. **M4c ✅ done 2026-11-09** — generated `src/nuclear_data.cyr` (the 18 ISOTOPES rows: label blob, `u8` Z, `u16` A, half-life, a UTF-8 decay-mode blob) + `src/nuclear.cyr` (decay constants and laws, mass defect and binding energies, the SEMF → `Result` by ADR 0002 #3, Q-values over slices, the two-step decay chain into the caller's buffer with Rust's equal-rate limit, the secular-equilibrium test, `lookup_isotope` handles), every transcendental correctly rounded; by the user's decision on the review, ADR 0002 #13 (N₀/2ⁿ rounded once where Rust's 2⁻ⁿ is not a normal double: correctly rounded for an integer n or N₀ = ±2^j, faithful otherwise), #14 (the exact λ_p < λ_d/100) and, for M4b's inorganic, #15 (the correctly rounded 4Δ/9 past Δ·4's overflow). Golden `nuclear` 54,057 rows (the SEMF over every Z ≤ A ≤ 3Z with A ≤ 300) + 23 corrected rows, `inorganic` 15,162: **219,266 ok, 0 failed, 0 measured, 0 pending, 3,243 accepted (ADR 6: 688; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,812 N₀/2ⁿ + 208 secular + 59 Δ_tet, all but the Antoine ones by name), 0 glibc-matching** on all three targets, every row 0 ulp from its reference; 36 Rust tests ported; surface parity 347/347. **M4a ✅ done 2026-11-09** — correctly rounded `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`, `_km_cr_atan2`, `_km_cr_acos` (two-phase; Cody–Waite and Payne–Hanek reduction, generated tables, `_km_td_norm` for the accurate division and square root), verified against MPFR by `cr_verify.py` on about 1.0–1.1 million arguments each (100,000 through the accurate path) on all three targets. **M4d ✅ done 2026-11-09** — `src/potential.cyr` (Lennard-Jones, Morse, Coulomb, distances, bond and dihedral angles through three-f64 point pointers, on `_km_cr_exp` and `_km_cr_atan2`); gas.cyr on M4a's `_km_cr_acos` / `_km_cr_cos` (every gas row keeps its bits) and ADR 0005's double-double `_km_cos` / `_km_acos` removed; the `math` trig rows at bound 0 (no golden transcendental row on the stdlib); ADR 0002 #13 correctly rounded (`_km_cr_exp2_scaled`), the new #16 (σ·2^(1/6) correctly rounded), #17 (Lennard-Jones past 4·ε's overflow), #18 (the Morse width past 2·D_e's overflow). Golden `potential` 26,938 rows + 1 corrected: **244,803 ok, 0 failed, 0 measured, 0 pending, 4,644 accepted (ADR 6: 689; ADR 2: 440 Antoine + 13 CFSE + 23 SEMF + 1,810 N₀/2ⁿ + 208 secular + 59 Δ_tet + 1,041 r_min + 180 Lennard-Jones + 181 Morse width), 0 glibc-matching** on all three targets, every row 0 ulp from its reference; 33 Rust tests ported; surface parity 357/357. **M4e** (spectroscopy) is next |
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
  rounded primitives, ADR 0006: the two-phase `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`, `_km_cr_ln`, `_km_cr_log10`,
  `_km_cr_pow` (M3c-2), `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`, `_km_cr_atan2`, `_km_cr_acos` (M4a) and
  `_km_cr_exp2_scaled` (M4d), each with `_acc` / `_fast` test entry points, and `_km_cr_sixth_root2_mul` (M4d, one
  proven rounding); ADR 0005's double-double `_km_acos` / `_km_cos` were removed in M4d), `src/validate.cyr`, `src/element_data.cyr`
  (generated: 118 elements, 32 property rows, `KmElementCategory`), `src/element.cyr`, `src/molecule.cyr`,
  `src/reaction.cyr`, `src/thermo.cyr`, `src/thermochem_data.cyr` (generated: `KmShomate`, 76 THERMOCHEM_DATA rows
  with formula / name blobs, 34 SHOMATE_DATA segments), `src/thermochem.cyr`, `src/phase_data.cyr` (generated:
  `KmAntoineCoeffs`, 8 Antoine sets, the
  mmHg → Pa factor, K_b/K_f), `src/phase.cyr`, `src/gas_data.cyr` (generated: `KmPengRobinsonParams`, STP, the 5
  `KM_PR_*` presets), `src/gas.cyr`, `src/solution.cyr`, `src/electrochemistry_data.cyr` (generated: the 49
  STANDARD_POTENTIALS rows — UTF-8 equation and couple blobs, electrons, E° — and `KM_STANDARD_POTENTIALS_COUNT`),
  `src/electrochemistry.cyr`, `src/kinetics.cyr`, `src/organic_data.cyr` (generated, M4b: `KmBondType`, `KmGeometry`,
  `KmFunctionalGroup` and their property columns, ADR 0002 #6's acceptor column with its sources), `src/organic.cyr`,
  `src/inorganic_data.cyr` (generated, M4b: `KmLigand` and its field strengths, Rust's private octahedral d-electron
  tables, the 10 `KM_MADELUNG_*` / `KM_BORN_EXP_*`, the 30 IONIC_RADII rows), `src/inorganic.cyr`,
  `src/potential.cyr` (M4d: Lennard-Jones, Morse, Coulomb and the 3D geometry; no data file), `src/nuclear_data.cyr`
  (generated, M4c: the 18 ISOTOPES rows, a UTF-8 decay-mode blob, `KM_ISOTOPES_COUNT`),
  `src/nuclear.cyr` (its `q_value` sums with thermo's `_km_sum`); plus the
  `src/main.cyr` stub. Generated files are written by
  `scripts/oracle/` and never edited. Every transcendental the modules call is correctly rounded; `gas` and `nuclear`
  include `lib/ganita.cyr` (cbrt), `gas` `lib/math.cyr` (`f64_clamp`; `mathx.cyr` no longer needs it since M4d). `dist/kimiya.deps`: math, string, vec, str, ganita, tagged; the bundle also
  includes `lib/hisab.cyr` (thermo's and thermochem's `EPSILON_F64`), so a consumer needs kimiya's `[deps.hisab]`.

## Tests

- `tests/parity.tcyr` (+ `tests/parity_core.cyr`, `tests/ulp.cyr`): the golden-corpus harness — 16 modules
  (constants, electrochemistry, element, gas, inorganic, kinetics, math, molecule, nuclear, organic, phase, potential,
  reaction, solution, thermo, thermochem), 244,803 ok, 0 failed, 0 measured, 0 pending, 4,644 accepted (689 under ADR 6,
  against the corrected oracle; under ADR 2 the 440 ADR 0002 #4 mapped Antoine rows, the 13 #12 exact-CFSE and 23 #3 SEMF
  rows printed by name, and, counted, the 1,810 #13, 208 #14, 59 #15, 1,041 #16, 180 #17 and 181 #18 rows whose note names
  the port's exact value), 0 glibc-matching, on x86_64 ELF, x86_64 PE and aarch64; ADR 0002 #6's renamed items compared
  under the Rust names (class A, summarised). Every class-C function and every `math` primitive at bound 0 (ADR 0006;
  sin, cos, atan2 and acos on M4a's primitives since M4d)
- `tests/golden/`: `constants.tsv` (45, all compared since M4b), `electrochemistry.tsv` (18,574), `element.tsv` (975),
  `gas.tsv` (21,093), `inorganic.tsv` (15,162), `kinetics.tsv` (30,146), `math.tsv` (26,147), `molecule.tsv` (2,472),
  `nuclear.tsv` (54,057), `organic.tsv` (9,104), `phase.tsv` (3,646), `potential.tsv` (26,938), `reaction.tsv` (1,558),
  `solution.tsv` (28,566), `thermo.tsv` (1,107),
  `thermochem.tsv` (9,857), `modules.txt` — generated by
  `scripts/oracle/`; the corrected oracle (ADR 0006) `gas.cr.tsv` (11), `kinetics.cr.tsv` (323), `math.cr.tsv` (96),
  `nuclear.cr.tsv` (23), `phase.cr.tsv` (4), `potential.cr.tsv` (1), `reaction.cr.tsv` (2), `solution.cr.tsv` (229) —
  written by `validate.py --write-corrections`
- `tests/mathx.tcyr`: the correctly rounded primitives (M3c-2, M4a, M4d) — special values, exact and tie cases (10²³, the
  subnormal ties), thresholds, subnormal results, the known hardest-to-round arguments, Payne–Hanek arguments, sincos,
  atan2's sign and zero matrix, acos, x·2^y, the fast path's fallbacks, and seeded fast / accurate agreement, every
  expected value MPFR's (20 groups, 1,384 assertions); `scripts/oracle/cr_verify.py` checks ≥ 1,000,000 arguments per
  primitive (611,518 pairs for x·2^y) on every target
- Reference tests, the Rust tests ported one for one: `tests/error.tcyr` (2 Rust + 1 Cyrius-only group,
  24 assertions), `tests/validate.tcyr` (3 + 1, 25), `tests/element.tcyr` (22 + 1, 96), `tests/molecule.tcyr`
  (11 + 1, 39; `serde_roundtrip` doesn't port), `tests/reaction.tcyr` (10 + 1, 29), `tests/thermo.tcyr`
  (10 + 1, 24), `tests/thermochem.tcyr` (62 + 6, 400; 4 tests assert J/mol for ADR 0002 #1), `tests/phase.tcyr`
  (31 + 2, 117), `tests/gas.tcyr` (39 + 2, 159), `tests/solution.tcyr` (49 + 5, 135), `tests/electrochemistry.tcyr`
  (35 + 3, 1,394), `tests/kinetics.tcyr` (47 + 4, 138; 2 tests adapted for the Lineweaver–Burk buffer),
  `tests/organic.tcyr` (21 + the 2 organic integration tests + 4, 136; 5 tests use ADR 0002 #6's names; the acceptor
  property's reference tests), `tests/inorganic.tcyr` (30 + 7, 332; the exact CFSE of ADR 0002 #12, the exact Δ_tet of #15), `tests/nuclear.tcyr`
  (36 + 8, 243; `semi_empirical_invalid_is_nan` asserts ADR 0002 #3's Err; #13's correctly rounded and #14's exact values),
  `tests/potential.tcyr` (33 + 5, 126; `morse_width_parameter_matches_force_constant` computes spectroscopy's √(k/μ)
  inline until M4e; null pointers, the helpers, the tolerances to the ulp and the −π → π mapping, #16's r_min, #17 / #18)
- `tests/kimiya.tcyr`: smoke test of `dist/kimiya.cyr` as a consumer includes it (41 assertions)
- `tests/kimiya.bcyr`: `element/lookup_by_symbol`, `element/lookup_by_number`, `molecule/molecular_weight_glucose`,
  `reaction/gibbs_free_energy`, `thermo/heat_transfer`, `gas/ideal_gas_pressure`, `gas/peng_robinson_z_n2_1mpa`,
  `thermochem/reaction_enthalpy_ch4`, `thermochem/vant_hoff_k`, `thermochem/enthalpy_change_cp`,
  `thermochem/shomate_cp`, `thermochem/adiabatic_flame_ch4_air`, `solution/ph_from_h`, `solution/weak_acid_ph`,
  `electrochemistry/lookup_half_reaction`, `electrochemistry/nernst_potential`, `kinetics/arrhenius_rate`,
  `kinetics/eyring_rate`, `kinetics/michaelis_menten`, and (Cyrius-only) `mathx/cr_exp`, `cr_expm1`, `cr_exp10`,
  `cr_ln`, `cr_log10`, `cr_pow` beside `ganita_f64_pow` (two arguments each), `organic/enthalpy_from_bonds`,
  `inorganic/lookup_ionic_radius` (M4b; rust-old/benches has no organic or inorganic bench),
  `nuclear/binding_energy_per_nucleon_fe56`, `nuclear/semi_empirical_fe56` (M4c, against Rust) and the Cyrius-only
  `nuclear/decay_chain_two_step`, the M4a trig benches (`mathx/cr_sin`, `cr_cos`, `cr_sincos`, `cr_atan2`, `cr_acos`
  beside `f64_sin`, `f64_cos`, `ganita_f64_atan2`), `potential/dihedral_angle` (M4d, against Rust) and the Cyrius-only
  `potential/bond_angle`, `lennard_jones`, `lj_equilibrium_distance`, `morse`, `mathx/cr_exp2_scaled`
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

**M4e**: `spectroscopy` (Rydberg / R_H, Beer–Lambert, `find_peaks`, the FFT via hisab's `num_fft`, isotope patterns with
their pruned convolution; its trigonometry is M4a's), then `stoichiometry` (exact integer balancing, rational oxidation
states, ADR 0002 #7) to complete M4.

Follow-ups:
- **CI: a truncated parity run.** Concurrent `cyrius test --win` runs under wine can print a truncated parity run and
  still exit 0, so the CI gate must check that the parity summary line (`=== parity summary: … ===`) is present (and
  that `cyrius test` reports every test file). The harness cannot detect this itself: whether wine cuts the program's
  output or ends it early and still reports 0, the exit status and the output are what wine hands back, so nothing the
  program prints or returns can mark a cut-off run. Running the three targets one after another (never two wine runs at
  once) avoids it; M4d's gates did.
- **M4d, reported to the user (kept as Rust's under the policy; no single exact value to restore):** in `potential`, the
  over- and underflow of (σ/r)⁶ and (σ/r)¹² (`lennard_jones` for σ/r past 2^85.3 or below 2^-89.5, though 4ε·d may be
  representable), of k·q₁·q₂ and of r·r (`coulomb`, `coulomb_force` for r past 1.3 × 10¹⁵⁴ or below 1.5 × 10⁻¹⁵⁴), of the
  squared coordinates in `distance_3d`, `bond_angle` and `dihedral_angle` (beyond about 10¹⁵⁴, or 10⁷⁷ for the dihedral's
  fourth powers, and below the reciprocals: an Err, 0, ±∞ or a degraded value where the distance or angle is
  representable, though Rust documents the angles as scale-invariant: an exact power-of-two prescaling would fix it and
  keep Rust's bits wherever it doesn't overflow, if the user wants it as an ADR 0002 row), `morse`'s r − r_e overflow and
  its subnormal D_e·x. Fixed in M4d: #16, #17, #18. The two numerical tolerances (10⁻¹², 10⁻²⁰) are kept, as
  `EQUAL_RATE_REL_TOL`.
- **ADR 0002 #13 is correctly rounded since M4d** (`_km_cr_exp2_scaled`; M4c's faithful fractional-n path is gone).
  `decay_constant`'s LN_2/t½ (1 ulp from CR(ln 2/t½) on 476 of 2,175 rows) is kept under the policy.
- **Two-argument primitives' hard cases.** `_km_cr_exp2_scaled`, like pow and atan2, is proven correct for results not
  within 2⁻¹⁵⁰ (relative) of a midpoint; no exhaustive worst-case search exists for two arguments (a random pair is that
  close with probability about 2⁻⁹⁶). `_km_cr_sixth_root2_mul` needs no such caveat (one argument: its continued-fraction
  bound is proven).
- **M4b, open questions for the user** (port plan §6.4): under ADR 0002 #12's policy, `delta_tet_from_oct` differs from
  CR(4Δ/9) only where Δ·4 overflows (|Δ| > DBL_MAX/4, 5 golden rows: decided in M4c, ADR 0002 #15), and Born–Landé's
  (1 − 1/n) from CR((n − 1)/n) on 332 rows inside a product with no exact value (left as Rust's). `organic.rs`'s bond energies carry no source
  citation in the oracle (roadmap, *Data to re-check*). The IUPAC 2011 acceptor wording in `src/organic_data.cyr` is
  quoted from a secondary copy; the published PDF was unreachable.
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
