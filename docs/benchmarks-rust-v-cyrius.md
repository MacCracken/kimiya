# Benchmarks — Rust 1.2.0 vs the Cyrius port

Each row times one benchmark body: the same work in Rust (`rust-old/`, release build) and in
Cyrius (`tests/kimiya.bcyr`, `cyrius bench`). Port milestones append rows as they land. Numbers
are host-specific, so compare within a column set, never across hosts.

## Method

- **Rust, criterion:** the *latest* column of [`rust-old/benchmarks.md`](../rust-old/benchmarks.md)
  (criterion 0.8 medians, v1.2 post-audit, rustc 1.99.1 release, 2026-11-08), from
  `cargo bench --manifest-path rust-old/Cargo.toml --bench benchmarks -- <name>`.
- **Cyrius:** `cyrius bench`, mean of 20 batches. Each batch is a tight loop of `per_batch` calls
  (as criterion's `b.iter` loop), inputs read from globals (Rust's `black_box`) and results written
  to a sink. lib/bench.cyr measures its clock-read floor and subtracts it from every sample
  (1.31 µs per read on this host); the per-op figure is lib/bench.cyr's ps mean.
- **Host:** x86_64 Linux, HPET clocksource. Toolchains: rustc 1.99.1 (release); Cyrius 6.7.18
  (ganita 1.2.15).

## M1 Foundation (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `element/lookup_by_symbol` | `lookup_by_symbol("Fe")`: a linear scan to row 26 | 44.73 ns | 67.0 ns | 1.5× |
| `element/lookup_by_number` | `lookup_by_number(26)`: range check, one table read | 0.84 ns | 4.54 ns | 5.4× |
| `molecule/molecular_weight_glucose` | `molecular_weight()` of C₆H₁₂O₆, built outside the loop: 3 lookups, 3 multiply-adds | 3.63 ns | 26.7 ns | 7.4× |

All three bodies are bit-identical to Rust (golden rows: `element` 975, `molecule` 2,472, 0 ulp).

**Two optimizations, measured.** The first version of each body measured 92.4 ns, 7.03 ns and
121 ns:
- `molecular_weight` read each atom through `atom_element_number(k)` / `atom_count(k)`, each a
  bounds check, a `vec_get` with its own two checks, and calls. It now walks the vec's buffer
  directly (lib/vec.cyr documents the `{ data, len, cap }` layout) and reads the mass column at
  the checked row: 121 → 57 ns with one `vec_get` per word, → 33 ns walking the buffer.
- `element_lookup_by_number` called `_km_element_row` for its range check; the check is inline
  now: 7.03 → 4.54 ns, and `molecular_weight`, which makes three of them per call, 33 → 26.7 ns.
- `element_lookup_by_symbol` compares the first byte before calling `memeq`, with the probe's
  first byte read once: 92.4 → 67.0 ns.

**Reading the gap.** Rust inlines `lookup_by_number` into a compare and a load, and unrolls the
glucose sum; Cyrius compiles each call as a call (the lookup inside `molecular_weight`, the
`Result` pair return). The symbol scan is the closest, because both sides do the same linear
string search and Rust's per-row cost is a real comparison too. Parity, not speed, is M1's bar.

## M2a reaction, thermo, phase, gas (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `reaction/gibbs_free_energy` | `gibbs_free_energy(-100000, 298, -200)`: three guards, ΔH − T·ΔS, a `Result` | 1.75 ns | 16.9 ns | 9.7× |
| `thermo/heat_transfer` | `heat_transfer(1000, 4.184, 10)`: two multiplies | 0.95 ns | 5.12 ns | 5.4× |
| `gas/ideal_gas_pressure` | `ideal_gas_pressure(1, 300, 0.025)`: three guards, nRT/V, a `Result` | 2.18 ns | 16.7 ns | 7.7× |
| `gas/peng_robinson_z_n2_1mpa` | `peng_robinson_z(PR_NITROGEN, 300 K, 1 MPa)`: validation, a·α and b, A and B, the one-root (Cardano) branch, two Newton steps | 61.39 ns | 270 ns | 4.4× |

All four bodies are bit-identical to Rust (golden rows: `reaction` 1,188, `thermo` 1,107, `gas` 19,717; the PR
roots 0 ulp on every target).

**One optimization, measured.** The NaN/±∞ guards were three nested calls each (`_km_is_positive` →
`_km_is_finite` → `_km_bits`). `_km_is_finite(x)` is now `x − x == 0` and `_km_is_positive` /
`_km_is_non_negative` compare against +∞ directly: `gibbs_free_energy` 23.3 → 16.9 ns (−27%),
`ideal_gas_pressure` 26.7 → 16.7 ns (−37%), `peng_robinson_z_n2_1mpa` 289 → 270 ns (−7%). The M1 benches don't
use the guards and are unchanged (65.3, 4.48, 26.7 ns).

**Reading the gap.** As in M1, Rust inlines each function into a few instructions; Cyrius makes the calls and
returns the `Result` as a pair. In `peng_robinson_z`, ganita's correctly rounded `cbrt` alone takes about 99 ns
(Rust links CORE-MATH's). The three-root branch, not benched in Rust, adds the port's correctly rounded acos
and three cos calls (ADR 0005, double-double): the price of bit-identical roots at the spinodals.

## M2b thermochem (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `thermochem/reaction_enthalpy_ch4` | `reaction_enthalpy([CO2(g) 1, H2O(l) 2], [CH4(g) 1, O2(g) 2])`: four table scans, two sums (and, ADR 0002 #1, one × 1000) | 38.58 ns | 250 ns | 6.5× |
| `thermochem/vant_hoff_k` | `vant_hoff_k(100, -50000, 298.15, 400)`: four guards, one exp | 9.51 ns | 58.9 ns | 6.2× |
| `thermochem/enthalpy_change_cp` | `enthalpy_change_cp("CO2(g)", 300, 600)`: the segment lookup, guards, one segment's antiderivative twice | 24.00 ns | 166 ns | 6.9× |
| `thermochem/shomate_cp` | `s.cp(500)`, s = CO2's first segment (looked up once, outside the loop) | 2.34 ns | 14.5 ns | 6.2× |
| `thermochem/adiabatic_flame_ch4_air` | `adiabatic_flame_temperature(-802300, [CO2(g) 1, H2O(g) 2, N2(g) 7.52], 298.15)`: three segment lookups, the bracket, Newton (3 updates on this input) | 301.48 ns | 3.24 µs | 10.7× |

All five bodies are bit-identical to Rust (golden rows: `thermochem` 9,487, 0 ulp on every class-B row; `vant_hoff_k`
within K = 6).

**Two optimizations tried, neither kept.** Testing the first byte inline before each Shomate formula compare (as
`element_lookup_by_symbol` does) left `enthalpy_change_cp` at 164–166 ns: CO2 is the third row, and the cost is the
`strlen` (16 ns) and the two `memeq` calls (13 ns each) that the run's length needs. Evaluating the antiderivative
straight from the columns, without the six-argument call to the shared `_km_shomate_h5`, took it from 22 to 20 ns
(the flame makes about 90 per call): not worth a second copy of the formula. Of `enthalpy_change_cp`'s 166 ns, the
segment lookup is 70 and the two-antiderivative integral 69.

**Reading the gap.** As in M1 and M2a: Rust inlines the Shomate polynomial, the `Segments` walk and the closures
into straight-line code; Cyrius makes each a call and returns the `Result` as a pair. The flame temperature is the
widest gap because it multiplies the per-call cost: its 3 µs are about six balance evaluations (3 products × up to
3 segments × 2 antiderivatives each) and three derivative evaluations. The port keeps Rust's resolved-segments table
on the stack (32 products), since the bump allocator can't free and the function runs in loops.

## M3a solution, electrochemistry (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `solution/ph_from_h` | `ph_from_h_concentration(1e-4)`: one guard, one log10 (`_km_log10`, fdlibm) | 9.73 ns | 66.2 ns | 6.8× |
| `solution/weak_acid_ph` | `weak_acid_ph(1.8e-5, 0.1)`: two guards, the charge-balance Newton (4 updates on this input), one log10 | 34.85 ns | 182 ns | 5.2× |
| `electrochemistry/lookup_half_reaction` | `lookup_half_reaction("Cu2+/Cu")`: a linear scan to row 27 | 57.06 ns | 116 ns | 2.0× |
| `electrochemistry/nernst_potential` | `nernst_potential(0.342, 2, 298.15, 0.1)`: four guards, one ln — the correctly rounded `_km_ln` | 6.42 ns | 1.94 µs | 302× |

All four bodies are bit-identical to Rust (golden rows: `solution` 28,566, `electrochemistry` 18,574; the Nernst rows 0 ulp
on every target, the pH rows within K = 4).

**The price of the correctly rounded ln.** `nernst_potential` cancels where E° ≈ (RT/nF)·ln Q (S ≈ 9 × 10¹⁸ with the
stdlib ln, port-parity.md §7), so it calls `_km_ln` (ADR 0005, amendment M3a), the double-double atanh series: measured
alone, `_km_ln(0.1)` takes 1.90 µs against 23 ns for `f64_ln`, so the Nernst body would be about 45–60 ns (7–9× Rust)
with the stdlib ln, and is 1.94 µs with the correctly rounded one. `henderson_hasselbalch` (not benched in Rust) pays the
same through `_km_log10_ln` (1.85 µs against 59 ns for `_km_log10`); the four S = 1 log10 functions keep the fast
`_km_log10`. A two-phase (Ziv) correctly rounded log — a fast evaluation with a rounding test, the double-double series
only when the test fails — would recover most of the cost; it is a possible follow-up, not in M3a.

**Reading the gap.** As before, Rust inlines each body; Cyrius makes the guard, primitive and `Result` calls. The lookup is
the closest, because both sides scan the same couple strings (the port tests the first byte before `memeq`, as the
thermochem and element lookups do). `ph_from_h`'s 66 ns are mostly `_km_log10` (59 ns measured alone).

## M3b kinetics (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `kinetics/arrhenius_rate` | `arrhenius_rate(1e13, 50000, 300)`: three guards, −Ea/(RT), one exp (`f64_exp`), A·e | 6.78 ns | 49.4 ns | 7.3× |
| `kinetics/eyring_rate` | `eyring_rate(80000, 298)`: two guards, k_B·T/h, one exp, a multiply | 6.73 ns | 42.4 ns | 6.3× |
| `kinetics/michaelis_menten` | `michaelis_menten(100, 5, 10)`: three guards, Vmax·[S]/(Km + [S]) | 2.92 ns | 16.5 ns | 5.7× |

`michaelis_menten` is class B and bit-identical to Rust (0 ulp on all 370 golden rows). `arrhenius_rate` and
`eyring_rate` are class C and measured, not yet gated (port-parity.md §8: a subnormal exp result scaled by a prefactor
puts their S past the rule's limit; ADR 0006 will set their comparison); on these benched inputs the exp results are
normal, and both bodies return Rust's bits on x86_64 and aarch64.

**Reading the gap.** As before, Rust inlines each body into a few instructions; Cyrius makes the guard calls and returns
the `Result` as a pair. The exp is `f64_exp` (x87 on x86_64); `michaelis_menten`'s three guards and arithmetic take
16.5 ns, so the exp and its call are most of the rest of the Arrhenius and Eyring bodies.

## M3c correct rounding (2026-11-09)

ADR 0006: every transcendental the modules call is now correctly rounded (`_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`,
`_km_cr_ln`, `_km_cr_log10`, `_km_cr_pow`; `_km_acos` / `_km_cos` unchanged). Every bench, before and after the switch,
with Rust's alongside. *Doc* is the number recorded at the bench's milestone above; *before* and *after* were measured
the same day on the same host, `cyrius bench` just before and just after M3c-3's switch (the host ran about 5–10 % slower
than at the earlier milestones, so compare *before* with *after*). Per-op means; ns unless marked.

| Bench | Rust | Cyrius, doc | Cyrius, before | Cyrius, after | after / Rust | What changed |
|---|---|---|---|---|---|---|
| `element/lookup_by_symbol` | 44.73 | 67.0 | 85.9 | 70.2 | 1.6× | nothing (noise) |
| `element/lookup_by_number` | 0.84 | 4.54 | 4.8 | 4.7 | 5.6× | nothing |
| `molecule/molecular_weight_glucose` | 3.63 | 26.7 | 28.8 | 29.1 | 8.0× | nothing |
| `reaction/gibbs_free_energy` | 1.75 | 16.9 | 18.6 | 18.8 | 10.7× | nothing |
| `thermo/heat_transfer` | 0.95 | 5.12 | 5.6 | 5.6 | 5.9× | nothing |
| `gas/ideal_gas_pressure` | 2.18 | 16.7 | 18.4 | 18.1 | 8.3× | nothing |
| `gas/peng_robinson_z_n2_1mpa` | 61.39 | 270 | 297 | 287 | 4.7× | nothing (the one-root branch: cbrt and sqrt) |
| `thermochem/reaction_enthalpy_ch4` | 38.58 | 250 | 276 | 267 | 6.9× | nothing |
| `thermochem/vant_hoff_k` | 9.51 | 58.9 | 65.5 | **77.1** | 8.1× | exp: `f64_exp` → `_km_cr_exp` |
| `thermochem/enthalpy_change_cp` | 24.00 | 166 | 183 | 171 | 7.1× | nothing |
| `thermochem/shomate_cp` | 2.34 | 14.5 | 15.6 | 15.5 | 6.6× | nothing |
| `thermochem/adiabatic_flame_ch4_air` | 301.48 | 3.24 µs | 3.49 µs | 3.43 µs | 11.4× | nothing |
| `solution/ph_from_h` | 9.73 | 66.2 | 69.9 | **87.3** | 9.0× | log10: `_km_log10` (fdlibm) → `_km_cr_log10` |
| `solution/weak_acid_ph` | 34.85 | 182 | 196 | **219** | 6.3× | the same log10 |
| `electrochemistry/lookup_half_reaction` | 57.06 | 116 | 124 | 120 | 2.1× | nothing |
| `electrochemistry/nernst_potential` | 6.42 | 1.94 µs | 86.7 | 83.9 | 13.1× | ln: already `_km_cr_ln` through the `_km_ln` wrapper since M3c-2a (1.94 µs → 87 ns then); now called directly |
| `kinetics/arrhenius_rate` | 6.78 | 49.4 | 52.0 | **65.5** | 9.7× | exp: `f64_exp` → `_km_cr_exp` |
| `kinetics/eyring_rate` | 6.73 | 42.4 | 44.6 | **59.0** | 8.8× | the same exp |
| `kinetics/michaelis_menten` | 2.92 | 16.5 | 17.4 | 17.3 | 5.9× | nothing |
| `mathx/cr_exp` (−50000/(8.314·300)) | — | — | 39.8 | 38.7 | — | Cyrius-only (M3c-2a) |
| `mathx/cr_expm1` (ln 2) | — | — | 51.3 | 49.6 | — | |
| `mathx/cr_exp10` (−4.75) | — | — | 54.4 | 53.5 | — | |
| `mathx/cr_ln` (0.1) | — | — | 49.0 | 49.1 | — | |
| `mathx/cr_log10` (10⁻⁴) | — | — | 77.4 | 78.7 | — | |
| `mathx/cr_pow` (0.1, 1.5) | — | — | 136 | 138 | — | Cyrius-only (M3c-2b) |
| `mathx/ganita_f64_pow` (0.1, 1.5) | — | — | 192 | 198 | — | what the modules called until M3c-3 |
| `mathx/cr_pow_1.409e-132` (1.409·10⁻¹³², 2.451) | — | — | 175 | 181 | — | |
| `mathx/ganita_f64_pow_1.409e-132` | — | — | 204 | 218 | — | |

**The price of correct rounding.** The switch costs 12–15 ns per exp call in a function body (`vant_hoff_k` +11.6,
`arrhenius_rate` +13.5, `eyring_rate` +14.4) and about 20 ns per log10 (`ph_from_h` +17.4, `weak_acid_ph` +23): measured
alone, `_km_cr_exp` takes 38–40 ns against about 11 ns for `f64_exp` (x87) and `_km_cr_log10` 77–79 ns against 59 ns for
`_km_log10` (M3a). pow got faster: `_km_cr_pow` is 136–181 ns against ganita's 192–218 (`nth_order_half_life` is not
benched). Redlich–Kwong's Tc^2.5 went from `_km_pow_2_5` (85 ns, double-double) to `_km_cr_pow(tc, 2.5)` (130 ns); it is
not benched in Rust. `nernst_potential`'s 1.94 µs → 84 ns (23× faster) is M3c-2a's, ADR 0005's double-double ln
replaced by the two-phase one. The benches whose bodies have no transcendental are unchanged within the host's noise.
ADR 0006 asks the fast paths to stay "within a small factor of the stdlib primitives": exp is 3.5×, log10 1.3×, pow 0.7×.
