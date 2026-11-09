# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Cyrius port, toward **2.0.0** (plan: `docs/development/kimiya-port-plan.md`). Rust 1.2.0 is frozen in
`rust-old/` as the parity oracle (ADR 0001).

### Added
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

### Breaking (relative to Rust 1.2.0, ADR 0002 #1)
- **thermochem** — `reaction_enthalpy` and `reaction_gibbs_energy` return **J/mol**, not kJ/mol: the value is 1.2's
  kJ/mol result × 1000 exactly (Rust's difference, then one multiply). `reaction_entropy` (J/(mol·K)) and the table
  fields `thermochem_data_delta_hf_kj` / `_delta_gf_kj` (kJ/mol, as tabulated) are unchanged.
  *Migration:* drop the `* 1000.0` before `vant_hoff_k` / `adiabatic_flame_temperature` —
  `adiabatic_flame_temperature(dh_kj * 1000.0, …)` becomes `adiabatic_flame_temperature(dh_j, products, n, t0)` with
  `var t, dh_j = reaction_enthalpy(p, np, r, nr);`; code that displays or compares the result in kJ/mol divides by
  1000; in ΔG = ΔH − TΔS with `reaction_entropy`, the `/ 1000.0` on TΔS goes away.

### Breaking (relative to Rust 1.2.0, ADR 0002 #4)
- **phase** — `antoine_pressure_mmhg` / `antoine_pressure_pa` return a `Result` instead of a bare value, and
  return `Err` where 1.2 extrapolated or returned NaN: `Err(KM_ERR_INVALID_TEMPERATURE)` for T outside the
  coefficient set's `t_min_c ..= t_max_c` (NaN included) or at/below the pole (C + T ≤ 0),
  `Err(KM_ERR_INVALID_INPUT)` for a non-finite A, B, C or range end. There is no extrapolating variant.
  *Migration:* `let p = antoine_pressure_mmhg(c, t);` becomes
  `var tag, v = antoine_pressure_mmhg(c, t); if (is_ok(tag) == 1) { var p: f64 = v; … }`; code that relied on
  extrapolation must use a coefficient set fitted over its temperatures (or `clausius_clapeyron`); code that
  tested `p.is_nan()` tests `is_err_result(tag)`. In range, the value is 1.2's to within the class-C bound.

### Changed
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
  depend on the platform's libm: x86_64 and aarch64 give the same bits

### Removed
- **mathx** (private, M3c-3): `_km_log10` and `_km_expm1` (M1's Rust libm 0.2.16 / fdlibm transcriptions),
  `_km_log10_ln` (M3a's glibc log10 composition, ADR 0005 amendment M3a, superseded by ADR 0006), the `_km_ln` wrapper,
  `_km_pow_2_5` (ADR 0005's double-double x^2.5, correctly rounded only from 2⁻³⁸⁷) and the unused `_km_dd_div`; their
  Cyrius-only tests now pin the correctly rounded bits (`tests/solution.tcyr`'s Henderson–Hasselbalch group: log10
  `…69DB`, glibc's `…69DA`)
- `validate.py`'s `_km_log10_ln` sensitivity model, its "S before" mode and the e_log10.c composition check

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

[1.2.0]: https://github.com/MacCracken/kimiya/releases/tag/1.2.0
[1.0.0]: https://github.com/MacCracken/kimiya/releases/tag/v1.0.0
