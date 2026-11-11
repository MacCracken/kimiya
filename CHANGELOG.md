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
  NaN). Every other value is 1.2's, bit for bit.
  *Migration:* none for physical well depths.

### Breaking (relative to Rust 1.2.0, ADR 0002 #18)
- **potential** — `morse_width_parameter(ω, μ, D_e)` for D_e ≥ 2¹⁰²³ J computes μ/(2·D_e) as (μ·½)/D_e, the correctly rounded
  value, where 1.2's 2·D_e overflowed and gave a = 0 (`morse_width_parameter(4, 1, 2¹⁰²³)` is 2⁻⁵¹⁰). Every other value is
  1.2's, bit for bit.
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
