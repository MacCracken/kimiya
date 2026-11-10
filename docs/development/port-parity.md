# kimiya port parity — measured baselines

Measurements that set the parity classes (port plan §5). The oracle is Rust 1.2.0 on
**glibc 2.44** (Rust std f64 methods call the host libm); the port is Cyrius 6.7.18 with
ganita 1.2.15. Since M3c (§9, §10; ADR 0006) a class-C row's reference is the corrected oracle — the
oracle's operations with every primitive result correctly rounded — and every class-C function is on
correctly rounded primitives and gated at 0 ulp; §2–§8 are the history of the glibc-relative bounds.

## 1. Semantics probe (M0 step 4, 2026-11-08, x86_64 Linux ELF)

A compiled probe of each Cyrius counterpart against the glibc value (via Python's `math`,
which is glibc). Bit patterns shown where they differ.

| Rust (std) | Cyrius | Result | Consequence |
|---|---|---|---|
| `f64::min(NaN, 1.0)`, `min(1.0, NaN)` → `1.0` both | `f64_min` | `min(NaN, 1)` = 1, but **`min(1, NaN)` = NaN** | asymmetric: Rust semantics need `_km_min` / `_km_max` (20 call sites) |
| `f64::max` | `f64_max` | same asymmetry | `_km_max` |
| `f64::round(2.5)` = 3, `round(-2.5)` = −3 (half away from zero) | `f64_round` | 2, −2 (**ties to even**) | `_km_round_half_away` (4 sites) |
| `f64::signum(0.0)` = 1.0, `signum(NaN)` = NaN | `f64_sign` | 0.0 and 0.0 | `_km_signum` |
| `log10`, `exp_m1` | — | no builtin, no ganita fn | `_km_log10`, `_km_expm1` (class C or transcribed) |
| `cbrt(2)` | `ganita_f64_cbrt` | **bit-identical** (`0x3ff428a2f98d728b`) | class A candidate |
| `atan2(1, 2)` | `ganita_f64_atan2` | **bit-identical** (`0x3fddac670561bb4f`) at this point; 1 ulp on 203 of 1,000 grid rows (§2) | class C |
| `ln(2)` | `f64_ln` and `_f64_ln_polyfill` | **bit-identical** (`0x3fe62e42fefa39ef`) | class C, measure |
| `exp(1)` | `f64_exp` (x87 on x86_64) | **bit-identical** (`0x4005bf0a8b145769`) at this point; 1 ulp on 45 of 1,125 grid rows (§2) | class C, per target |
| `exp(1)` | `_f64_exp_polyfill` (same bits on every target) | `…576a`, **1 ulp** | — |
| `acos(0.5)` | `ganita_f64_acos` | `…7365` vs glibc `…7366`, **1 ulp** | class C |
| `powf(2, 0.5)` | `ganita_f64_pow` | `…3bcc` vs glibc `…3bcd` (correctly rounded √2), **1 ulp** | class C |

## 2. Corpus baselines (M0 step 8, 2026-11-08)

`tests/golden/math.tsv` (18,647 rows, 24 functions, glibc 2.44 x86_64 with FMA) replayed by
`tests/parity.tcyr` through each Cyrius counterpart. Histogram = ulp distance from the oracle.

| Rust → Cyrius | Rows | x86_64 ELF | x86_64 PE (wine) | aarch64 (qemu) |
|---|---|---|---|---|
| exp → `f64_exp` | 1,125 | 0: 1,080 · 1: 45 | = ELF | 0: 1,026 · 1: 99 (the polyfill) |
| ln → `f64_ln` | 1,357 | 0: 1,357 | = ELF | 0: 1,341 · 1: 16 (the polyfill) |
| powf → `ganita_f64_pow` | 750 | 0: 697 · 1: 52 · ∞: 1 (ADR 0004) | = ELF | = ELF |
| powf, base 10 | 571 | 0: 506 · 1: 65 | = ELF | = ELF |
| powf, base 2 (glibc `exp2`, see §3) | 369 | 0: 339 · 1: 30 | = ELF | = ELF |
| cos → `f64_cos` | 910 | 0: 892 · 1: 18 | = ELF | = ELF |
| sin → `f64_sin` | 710 | 0: 690 · 1: 20 | = ELF | = ELF |
| acos → `ganita_f64_acos` | 620 | 0: 463 · 1: 157 | = ELF | = ELF |
| atan2 → `ganita_f64_atan2` | 1,000 | 0: 797 · 1: 203 | = ELF | 0: 798 · 1: 202 |
| sqrt, abs, floor, trunc, cbrt, clamp (class A) | 5,915 | all 0 | all 0 | all 0 |

- **Class-C bound for the math primitives: K = 2 ulp** on every target (each side documents
  < 1 ulp; the observed maximum is 1). Set in `_par_c_bound`; gated from now on.
- **One accepted row:** `powf(-1, DBL_MAX)`, ganita's |x| == 1 defect (ADR 0004).
- **Pending until the M1 helpers:** `round`, `min`, `max`, `signum`, `total_cmp`, `is_finite`,
  `powi`, `exp_m1`, `log10` (5,320 rows) and the 31 module-data constants. Done in M1 (§4).
- **Result, all three targets (M0):** 13,340 ok, 0 failed, 0 measured, 5,351 pending, 1 accepted.
- **glibc vs correct rounding** (validate.py, 60-digit `decimal` reference): glibc is correctly
  rounded on every row except exp_m1 (34 rows), log10 (48), powf (1), exp2 (1) and sin (2), all
  1 ulp. Where the Cyrius result is 1 ulp off on the other rows, the Cyrius side is the one off.
- **Polyfill shadows** (reported, never gated): `_f64_exp_polyfill` 99/1,125 at 1 ulp,
  `_f64_ln_polyfill` 16/1,357 at 1 ulp on x86_64 — the "same bits on every target" variants are
  not closer to glibc than the builtins.

Rule for kimiya's own class-C functions (from M2): K_F = 2 · S · (class-C calls on F's path) + 2,
where S is F's measured sensitivity, set before F's first comparison. If K_F > 64, make the
primitive exact instead (transcribe it, as tanmatra did for sin/cos/log). glibc can't be transcribed,
so M2a made four primitives correctly rounded instead, and S counts only the calls whose result can
differ between the port and the oracle ([ADR 0005](../adr/0005-correctly-rounded-primitives-for-ill-conditioned-functions.md), §5).

## 3. Oracle provenance findings (M0 step 5)

- **Not every f64 method is glibc.** With rustc 1.99.1, `cbrt`, `floor`, `trunc`, `round` and
  `powi` come from Rust's own statically linked math (cbrt is correctly rounded). Only exp,
  exp2, expm1, log, log10, pow, sin, cos, acos and atan2 call glibc.
- **The oracle is the release build.** In release, LLVM rewrites `2.0_f64.powf(y)` into glibc
  `exp2(y)` (the generator binary imports `exp2`); in debug, where the Rust tests ran, it stays
  `pow`. kimiya calls `2.0_f64.powf` in `nuclear::atoms_after_half_lives` and
  `potential::lj_equilibrium_distance`. The corpus records release values (what consumers ship),
  and keeps the literal-base forms as separate functions (`powf_base10`, `powf_base2`).
- **Platform pin.** Every TSV header records `libm: glibc 2.44 x86_64 fma`; `--check` refuses a
  different glibc or CPU class, since glibc selects FMA code paths per CPU.

## 4. M1: the local math helpers, element and molecule (2026-11-09)

`src/mathx.cyr` closes the 5,320 pending `math` rows (port plan §6.4). Histogram = ulp distance
from the oracle; identical on x86_64 ELF, x86_64 PE (wine) and aarch64 (qemu).

| Rust → Cyrius | Rows | Class | Result, all three targets |
|---|---|---|---|
| round → `_km_round_half_away` (`trunc(x + copysign(0.5 − ε/4, x))`, Rust's libm generic round) | 661 | A | 661 bit-identical |
| max / min → `_km_max` / `_km_min` (NaN operand ignored; equal operands, ±0 included, give the second) | 500 / 500 | A | all bit-identical |
| signum → `_km_signum` | 374 | A | 374 bit-identical |
| total_cmp → `_km_total_cmp` | 256 | A | 256 equal |
| is_finite → `_km_is_finite` | 33 | A | 33 equal |
| powi → `_km_powi` (compiler-builtins `__powidf2` loop) | 890 | B | 890 at 0 ulp |
| exp_m1 → `_km_expm1` (Rust libm 0.2.16 `expm1.rs`, FreeBSD s_expm1.c) | 810 | C, K = 2 | 0: 810 |
| log10 → `_km_log10` (Rust libm 0.2.16 `log10.rs`, FreeBSD e_log10.c) | 1,296 | C, K = 2 | 0: 1,248 · 1: 48 |

- **K = 2 for exp_m1 and log10**, the math-primitive rule of §2 (each side documents < 1 ulp; log10's
  glibc bound is 2), set after the rows were first replayed as *measured* and seen to stay within 1 ulp
  on every target. Not widened: the observed maximum is 1.
- **The 48 log10 rows are glibc's.** A Python transcription of the same libm code reproduces the 48
  differing rows, and on all 48 the transcription is the correctly rounded result (60-digit `decimal`
  reference); they are exactly the 48 rows where `validate.py` finds glibc 1 ulp off. glibc's expm1 is
  bit-identical to fdlibm's on every row, including the 34 where both are 1 ulp off correct rounding.
- The transcriptions use only + − × ÷, bit manipulation and `f64_trunc` / `f64_from` / `f64_to`, so the
  same bits come out on every target (no x87, no polyfill).
- **math: 0 pending.** The `f64_round` / `f64_min` / `f64_max` / `f64_sign` shadow lines stay, as the
  record of why the helpers exist.

| Module | Rows | Class | Result, all three targets |
|---|---|---|---|
| `element` (`lookup_by_number` 256, `lookup_by_symbol` 463, `lookup_properties` 256) | 975 | A | 975 ok, 0 ulp |
| `molecule` (presets, `new`, `molecular_weight`, `formula`, `total_atoms`) | 2,472 | A, B | 2,472 ok; `molecular_weight` 0 ulp on all 507 Ok rows |
| `constants` (the 14 physical) | 14 of 45 | A | 14 ok; 31 module-data rows pending (M2–M4) |
| `math` | 18,647 | A, B, C | 18,646 ok, 1 accepted (ADR 0004) |

Summary, every target: **4 modules, 0 failing; 22,107 ok, 0 failed, 0 measured, 31 pending, 1 accepted.**

## 5. M2a: reaction, thermo, phase, gas (2026-11-09)

Golden rows: `reaction` 1,188, `thermo` 1,107, `phase` 3,276, `gas` 19,717 (`m_reaction.rs`, `m_thermo.rs`,
`m_phase.rs`, `m_gas.rs`). Every class-A and class-B row (the arithmetic functions, the tables, the presets,
`phase_rule`, `peng_robinson_pressure`) is **0 ulp on all three targets**. `validate.py` recomputes all
25,288 rows from Rust's code transcribed in Python (Python's `math` is glibc, so the class-C rows too)
and finds the corpus's bits on every one.

### Setting the class-C bounds

> **Superseded by [ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md) (M3c-3, §10).** Every class-C function here
> is now on correctly rounded primitives and gated at 0 ulp against the corrected oracle. The S, c and K below are the
> history of the M2a rule; `--sensitivity` still prints S, now as a diagnostic of conditioning.

S was measured by the rule's perturbation method before the first gated comparison:
`python3 -I scripts/oracle/validate.py . --sensitivity` re-evaluates every golden row of each class-C
function with one primitive result nudged one ulp up, then down, one call at a time, and reports the largest
change of any output. (The rows were first replayed unbounded, as *measured*.)

**With the oracle's primitives, four functions were past K = 64** (S per primitive kind, every row):

| Function | S | Primitive | Cause |
|---|---|---|---|
| `peng_robinson_z_roots` | 2²⁶, 14 root-count flips | cos, acos | near a double root (a spinodal) the two close roots are conditioned like √ε; the flips are acos(−1) at P = DBL_MIN, where a double root at 0 sits on Z = B |
| `fugacity_coefficient_pr` / `fugacity_pr` | 318 / 492 | ln | φ = exp(ln φ), so a 1-ulp change in ln(Z − B) moves φ by about \|ln φ\| ulp (N₂ at 2 K) |
| `redlich_kwong_pressure` | 64 | pow(Tc, 2.5) | RT/(V − b) − a/(√T V(V + b)) cancels near P = 0 |
| `peng_robinson_z` | 2 | cos | the vapour-like root is well-conditioned |

glibc 2.44 can't be transcribed (no source in the oracle's dependency graph; FMA ifunc variants), and it is
not correctly rounded (20,000-argument sweeps: acos 0.07 %, cos 0.16 %, ln 0.04 %, pow(x, 2.5) 0.08 % of
arguments 1 ulp off). **ADR 0005**: the port's acos, cos (|x| ≤ 8), ln and x^2.5 are correctly rounded
(`_km_acos`, `_km_cos`, `_km_ln`, `_km_pow_2_5` in `src/mathx.cyr`, double-double, + − × ÷ and sqrt only), so
they return glibc's bits wherever glibc's are correctly rounded, and S counts only the calls that can differ:
every exp and 10^x result (stdlib counterparts), and the acos / cos / ln / x^2.5 results where glibc is not
correctly rounded (`validate.py` decides each golden argument against its 60-digit reference; 5 cos and 5 ln
calls in the grid). Results that are 0, ±∞ or NaN, and acos(±1), are exact on both sides and never nudged.

| Function | S (primitive) | c | K_F | x86_64 ELF and PE | aarch64 |
|---|---|---|---|---|---|
| `reaction::equilibrium_constant` | 1 (exp: K is the exp result) | 1 | **4** | 0: 629 · 1: 31 | 0: 601 · 1: 59 |
| `phase::clausius_clapeyron` | 2 (exp, through P₁·e) | 1 | **6** | 0: 445 | 0: 420 · 1: 22 · 2: 3 |
| `phase::antoine_pressure_mmhg` | 1 (pow: P is the pow result) | 1 | **4** | 0: 304 · 1: 34 | = x86_64 |
| `phase::antoine_pressure_pa` | 2 (pow, through × PA_PER_MMHG) | 1 | **6** | 0: 304 · 1: 33 · 2: 1 | = x86_64 |
| `gas::peng_robinson_z` | 0 (the 5 differing cos calls don't reach the largest root) | 1 | **2** | 0: 4,259 | = x86_64 |
| `gas::peng_robinson_z_roots` | 2 (cos, CO₂ at 280 K, 2.2 MPa) | 1 | **6** | 0: 5,189 roots | = x86_64 |
| `gas::fugacity_coefficient_pr` | 1 (exp) | 2 | **6** | 0: 4,241 · 1: 7 | 0: 4,008 · 1: 240 |
| `gas::fugacity_pr` | 2 (exp, through φ·P) | 2 | **10** | 0: 4,242 · 1: 5 · 2: 1 | 0: 4,035 · 1: 184 · 2: 29 |
| `gas::redlich_kwong_pressure` | 0 (glibc's Tc^2.5 is correctly rounded on every grid Tc) | 0 | **2** | 0: 547 | = x86_64 |

(Histograms count float outputs; Err rows have none. c is the most differing calls in one row.)

- **The Peng–Robinson roots and Redlich–Kwong pressures are bit-identical on every target**, the
  ill-conditioned rows included: the ill-conditioning is real, but both sides now evaluate the same
  correctly rounded primitives in the same order. The fugacities differ only through exp (x87 on x86_64, the
  fdlibm polyfill on aarch64), as do the equilibrium constant and Clausius–Clapeyron.
- **The four primitives against glibc** (shadow lines on the `math` corpus, never gated; identical on all
  three targets): `_km_cos` 2,660 rows with |x| ≤ 8 (4 differ), `_km_acos` 2,620 (1), `_km_ln` 3,857 (3),
  `_km_pow_2_5` 1,100 rows with y = 2.5 (2). Those are exactly the counts of rows where `validate.py` finds
  glibc 1 ulp off correct rounding.
- **The `math` corpus grew by 7,500 rows** (cos on [−2.2, 3.2], acos on [−1, 1] and toward ±1, ln on (0, 4],
  near 1 and wide, powf(x, 2.5) for x 0.1–2000): 26,147 rows. The stdlib counterparts stay within K = 2:
  `f64_cos` 0: 2,830 · 1: 80; `ganita_f64_acos` 0: 1,893 · 1: 726 · **2: 1** (aarch64 0: 1,876 · 1: 742 ·
  2: 2 — ganita's acos is not target-independent); `f64_ln` 0: 3,854 · 1: 3 (aarch64 polyfill 1: 153);
  `ganita_f64_pow` 0: 1,597 · 1: 152 and the ADR 0004 row.

### The Antoine mapping (ADR 0002 #4)

The port's `antoine_pressure_mmhg` / `_pa` return `Err` outside a set's `t_min_c ..= t_max_c`, at or below the
pole and for non-finite input. Of the 1,116 Antoine rows: 676 (338 per function) are in range with a finite
oracle value and are compared as values (class C; the unit set's two rows class A); 174 where the oracle is
NaN (the pole, non-finite fields or T) and 266 where T is out of range (the oracle extrapolated, or a range end
is infinite) are **mapped**: the port must return `Err`, and does, and they count as accepted under ADR 2.
An `Ok` on a mapped row, or an `Err` on an in-range row, would fail.

Summary, every target: **8 modules, 0 failing; 54,466 ok, 0 failed, 0 measured, 20 pending (the M3–M4
module-data constants), 441 accepted (ADR 0004's row and 440 ADR 0002 #4 mapped rows).**

## 6. M2b: thermochem (2026-11-09)

Golden rows: `thermochem` 9,487 (`m_thermochem.rs`). Every class-A and class-B row is **0 ulp on all three
targets**; `validate.py` recomputes all 9,487 from Rust's code transcribed in Python (the tables parsed from the
`thermochem.rs` literals; the flame solver with hisab 1.4.0's loops) and finds the corpus's bits on every one, and
`probe_tables.py` compiles the 654 `thermochem_data.cyr` cells.

| Function | Rows | Class | Result, all three targets |
|---|---|---|---|
| `THERMOCHEM_DATA` / `[]`, `SHOMATE_DATA` / `[]` | 1 + 76, 1 + 34 | A | all bit-identical |
| `lookup_thermochem`, `lookup_shomate`, `lookup_shomate_at` (shared boundaries: the lower segment) | 489, 489, 511 | A | all equal |
| `Shomate::cp`, `Shomate::enthalpy_antiderivative` (every segment; extrapolation, T = 0, ±∞, NaN) | 1,188 each | B | 0 ulp |
| `reaction_enthalpy`, `reaction_gibbs_energy` (J/mol, ADR 0002 #1), `reaction_entropy` | 429 each (422 Ok) | B | 0 ulp |
| `enthalpy_change_cp` (within and across segments, reversed, range ends ± 1 ulp) | 1,601 (619 Ok) | B | 0 ulp |
| `adiabatic_flame_temperature` | 1,937 (1,836 Ok) | B | 0 ulp |
| `vant_hoff_k` | 685 (649 Ok) | C, **K = 6** | x86_64 ELF and PE 0: 649 · aarch64 0: 606 · 1: 36 · 2: 7 |

**ADR 0002 #1.** The reaction-enthalpy and Gibbs-energy rows hold J/mol: the generator multiplies the oracle's
kJ/mol result by 1000.0 (the row note says so), and the port computes Rust's kJ difference then the same one
multiply, so the rows are compared directly as class B — 844 Ok rows at 0 ulp, 14 Err rows. The harness reports the
count after the module (`ADR 0002 #1 kJ -> J`); no row is mapped or accepted.

**vant_hoff_k's bound** (superseded by ADR 0006, §10: bound 0 on `_km_cr_exp`). One exp: K₂ = K₁·exp(−(ΔH/R)(1/T₂ − 1/T₁)). Measured before its first gated comparison
(the 649 Ok rows were replayed as *measured* first, 0 ulp on x86_64 ELF): `validate.py --sensitivity` gives
**S = 2** (the exp result, through K₁·e, as `clausius_clapeyron`), c = 1, so **K = 2·2·1 + 2 = 6**. Observed max 0
(x86_64) and 2 (aarch64, the exp polyfill).

**The flame solver** (ADR 0003, amendment M2b): hisab 1.4.0's `newton_raphson` and `bisection`, transcribed, so the
path is + − × ÷ and the rows class B. The generator rebuilds the solver from the public API (`flame_replica`), must
reproduce the oracle's bits on every row, and writes the solver path into each row's note; `validate.py` recomputes
the path independently and checks the note.

| Path | Rows |
|---|---|
| Newton (0 / 1 / 2 / 3 / 4 / 5 / 6 updates) | 1,682 (32 / 26 / 553 / 954 / 91 / 22 / 4) |
| Bisection after Newton's 50 iterations failed to reach \|f\| < 10⁻³ J/mol (each 33 halvings) | 154 |
| `Err` (66 not bracketable, 20 T₀, 11 input, 4 unknown product) | 101 |

The bisection rows are the mixtures with ΔH and moles scaled by 10⁶ … 10³⁰⁰: at |f| ~ 10¹⁵ J/mol the closest
representable T leaves |f| far above 10⁻³, so Newton exhausts its iterations unless f happens to round to near zero
(about half the scaled rows do; both paths are covered at every scale). The tiny scales (10⁻²⁰⁰ … 10⁻⁶) cover Newton
accepting its first guess and the bracket product underflowing to −0 (`ComputationError`, as Rust). No golden row
reaches Newton's "root outside the range" or "derivative is zero" exits (both replicas report only "no
convergence"): f decreases on the bracket and Newton starts inside it from the clamped guess, and where |f| ≥ 10⁻³
the products total ≳ 10⁻⁹ mol, so |f′| = Σn·Cp ≳ 10⁻⁸ ≫ 10⁻¹² (Cp > 20 J/(mol·K) in every segment's range).
`tests/thermochem.tcyr` drives the private Newton helper into the derivative exit directly.

**Typed constants pinned.** Each of Rust's inline constants changes golden bits when perturbed (`validate.py`): the
Newton tolerance 10⁻³ × ½ or × 2 changes 98 / 104 of the 1,836 Ok rows, the bisection tolerance 10⁻⁶ × ½ or × 2 all
154 bisection rows, the 40 J/(mol·K) guess as 39 or 41 changes 1,209 / 1,199. `tests/thermochem.tcyr` pins 10⁻³ and
10⁻⁶ exactly and the 10⁻¹² derivative threshold within 1 %.

Summary, every target: **9 modules, 0 failing; 63,953 ok, 0 failed, 0 measured, 20 pending (the M3–M4
module-data constants), 441 accepted (ADR 0004's row and the 440 ADR 0002 #4 mapped Antoine rows).**

## 7. M3a: solution, electrochemistry (2026-11-09)

Golden rows: `solution` 28,566, `electrochemistry` 18,574 (`m_solution.rs`, `m_electrochemistry.rs`), and 10 more
`constants` rows (K_w, the Debye–Hückel A and B, the seven Λ°, now in `src/constants.cyr`): 35 of 45, the other 10 (the
inorganic Madelung and Born constants) pending until M4. Every class-A and class-B row is **0 ulp on all three
targets**. `validate.py` recomputes all 47,140 rows from Rust's code transcribed in Python (the potentials parsed from
the `electrochemistry.rs` literal; the weak acid / base Newton with its exit path) and finds the corpus's bits on every
one, and `probe_tables.py` compiles the 197 `electrochemistry_data.cyr` cells (the UTF-8 equations byte for byte) and
the 10 new constants.

| Function | Rows | Class | Result, all three targets |
|---|---|---|---|
| `molarity`, `molality`, `dilution`, `ph_from_poh` | 191, 191, 208, 176 | B | 0 ulp |
| `ionic_strength` (`v:n` + `f:c i:z` pairs; the empty list's −0.0, i32::MIN / MAX charges) | 320 | B | 0 ulp |
| `kohlrausch_molar_conductivity`, `limiting_molar_conductivity`, `solution_conductivity` | 258, 270, 186 | B | 0 ulp |
| `STANDARD_POTENTIALS` / `[]` (49 rows: equation, couple, electrons, E°) | 1 + 49 | A | all bit-identical |
| `lookup_half_reaction` (every couple, case / space variants, every prefix, the equations, `""`) | 529 | A | all equal |
| `standard_cell_potential`, `cell_potential_from_couples` (all 49 × 49 ordered pairs, unknown couples) | 2,509, 2,408 | B | 0 ulp |
| `is_spontaneous_cell` | 74 | A | all equal |
| `cell_gibbs_energy`, `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` | 339, 3,975, 2,168, 2,199 | B | 0 ulp |
| `charge_to_moles_electrons`, `moles_electrons_to_charge` | 853 each | B | 0 ulp |

### Setting the class-C bounds

> **Superseded by [ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md) (M3c-3, §10).** The K column below is the M3a
> rule's; `_km_log10`, `_km_log10_ln` and ADR 0005's amendment M3a are retired (`_km_cr_log10`, `_km_cr_exp10` and
> `_km_cr_ln` now, every row at 0 ulp against the corrected oracle).

Every class-C row was replayed unbounded (*measured*) first; then `validate.py --sensitivity` measured S, the rule set
K = 2·S·c + 2, and the rows were gated. Three functions cancel by construction, and the grids include
**near-cancellation rows on purpose** (pKa = −log₁₀(base/acid), E° = (RT/nF)·ln Q, each to within a few ulp and placed
without a transcendental): `henderson_hasselbalch` 398 of its 831 value rows (90 exactly zero, 152 positive, 156
negative), `nernst_potential` 401 of 1,430 (81, 150, 170), `nernst_potential_25c` 201 of 1,106 (41, 80, 80). With the
stdlib-like primitive their K was far past 64, so ADR 0005's remedy applies (its **amendment M3a**):

| Function | S before (primitive) | Primitive in the port | S | c | K_F | x86_64 ELF, PE and aarch64 |
|---|---|---|---|---|---|---|
| `solution::ph_from_h_concentration` | 1 | `_km_log10` (fdlibm) | 1 | 1 | **4** | 0: 2,836 · 1: 2 |
| `solution::poh_from_oh_concentration` | 1 | `_km_log10` | 1 | 1 | **4** | 0: 2,836 · 1: 2 |
| `solution::henderson_hasselbalch` | 4.4 × 10¹⁸ (`_km_log10`) | `_km_log10_ln`: glibc's log10 composition over `_km_ln` | 0 | 0 | **2** | 0: 831 |
| `solution::h_concentration_from_ph` | 1 | ganita pow(10, ·) | 1 | 1 | **4** | 0: 392 · 1: 47 |
| `solution::weak_acid_ph` | 1 | the Newton root (bit-exact), `_km_log10` | 1 | 1 | **4** | 0: 8,865 · 1: 183 |
| `solution::weak_base_ph` | 1 | the Newton root, `_km_log10`, 14 − pOH | 1 | 1 | **4** | 0: 9,009 · 1: 39 |
| `solution::debye_huckel_limiting` | 1 | ganita pow(10, ·) | 1 | 1 | **4** | 0: 453 · 1: 35 |
| `solution::debye_huckel_extended` | 1 | ganita pow(10, ·) | 1 | 1 | **4** | 0: 475 · 1: 30 |
| `solution::davies_activity` | 1 | ganita pow(10, ·) | 1 | 1 | **4** | 0: 460 · 1: 28 |
| `electrochemistry::nernst_potential` | 8.8 × 10¹⁸ (`f64_ln`) | `_km_ln` (correctly rounded, M2a) | 0 | 0 | **2** | 0: 1,430 |
| `electrochemistry::nernst_potential_25c` | 8.7 × 10¹⁸ (`f64_ln`) | `_km_ln` | 0 | 0 | **2** | 0: 1,106 |

(Histograms count float outputs; Err rows have none. "S before" nudges every result of the primitive the remedy
replaced; S nudges only the calls whose result can differ between port and oracle, ADR 0005: every `_km_log10` and
pow(10, ·) result, and the `_km_ln` results where glibc's log is not correctly rounded — none on these grids.)

- **The Nernst and Henderson–Hasselbalch values are bit-identical to the oracle on every target**, the
  near-cancellation rows included. Before the remedy the replay showed the problem: Nernst through `f64_ln` was 0 ulp
  on x86_64 (x87) but had 6 rows up to 4.3 × 10¹⁸ ulp off on aarch64 (the polyfill).
- **A correctly rounded log10 was tried for `henderson_hasselbalch` first** (`_km_log10_cr`, double-double over
  `_km_ln`): 6 near-cancellation rows ended up to 4.4 × 10¹⁸ ulp off, because glibc's log10 is not correctly rounded at
  base/acid = 0.8068…/2.6858… — and is not on 0.3–29 % of arguments, depending on the range (29 % in [½, 1)). glibc's
  log10 is Sun's e_log10.c composition over glibc's log (bit-identical on 500,000 seeded arguments and every golden one),
  so the port composes it the same way over the correctly rounded `_km_ln`: `_km_log10_ln` matches glibc except where
  glibc's log at the reduced argument is not correctly rounded (13 of 23,710 seeded arguments, none golden; identical
  bits on all three targets). The `math` shadow line `f64::log10 ~ _km_log10_ln` is 0 ulp on all 1,296 rows.
- **`_km_log10` stays** on the four S = 1 functions (and the benched `ph_from_h`): `_km_log10_ln` costs 1.85 µs a
  call, `_km_log10` 59 ns.
- The weak-acid / base and activity rows differ from glibc only through `_km_log10` and ganita's pow (1 ulp at most),
  the same on every target (both are + − × ÷ transcriptions); the Newton root itself is + − × ÷ and √.

### The weak acid / base Newton

Rust's private `charge_balance_root` / `checked_root` are transcribed operation for operation
(`_km_charge_balance_root` / `_km_checked_root`). The generator rebuilds the loop (`root_replica`), must reproduce the
oracle's pH bits or error on every row, and writes each row's exit into its note; `validate.py` recomputes the path
independently and checks the note. Over the 18,222 weak acid and base rows (Ka × C every decade 10⁻⁴⁰ … 10⁴⁰ ×
10⁻⁴⁰ … 10²⁰, seeded pairs over that box and over Ka 10⁻¹⁴ … 10, C 10⁻¹² … 10, the overflow cases, the edges and the
Rust tests' cases):

| Exit | Rows |
|---|---|
| `f(x) ≤ 0` (the iterates met the root; 9,327 of them at the starting bound, 0 updates) | 11,740 |
| relative step `step ≤ 10⁻¹²·x` | 6,356 |
| non-finite step (K·C or x³ overflows) → `ComputationError` | 74 |
| 50 iterations without acceptance | 0 |
| `checked_root` rejection | 0 |
| guard `Err` before the iteration (Ka or C not positive and finite) | 52 |

Updates per row: 0: 9,327 · 1: 3,432 · 2: 2,264 · 3: 1,547 · 4: 824 · 5: 483 · 6: 247 · 7: 46 — **at most 7 updates,
7 iterations entered**, as Rust's comment says. The private constants are typed inline: `CHARGE_BALANCE_REL_TOL`
(10⁻¹²) × ½ changes 11 of the 18,096 Ok rows and × 2 changes 8 (the seeded sweep over the realistic box was added so
that the corpus pins the tolerance from above too: before it, × 2 changed none; the change is one ulp, from accepting
one step earlier); `tests/solution.tcyr` pins it from both sides on the private helper with two of those decade rows.
**`CHARGE_BALANCE_MAX_ITER` (50) is not observable**: no row enters more than 7 iterations, so any cap ≥ 7 gives the
same corpus (a cap of 6 changes 88 rows); it is transcribed, as the flame solver's limits are (ADR 0003 amendment).
No golden row reaches `checked_root`'s rejection either (the start is a positive bound and the iterates decrease onto
the positive root); `tests/solution.tcyr` checks the helper directly.

Summary, every target: **11 modules, 0 failing; 111,103 ok, 0 failed, 0 measured, 10 pending (the M4 inorganic
constants), 441 accepted (ADR 0004's row and the 440 ADR 0002 #4 mapped Antoine rows).**

## 8. M3b: kinetics (2026-11-09)

Golden rows: `kinetics` 30,146 (`m_kinetics.rs`). Every class-B row is **0 ulp on all three targets**. `validate.py`
recomputes all 30,146 rows from Rust's code transcribed in Python (glibc's exp, expm1 and pow through Python's `math`)
and finds the corpus's bits on every one.

| Function | Rows | Class | Result, all three targets |
|---|---|---|---|
| `half_life_first_order` (every decade of k from the smallest subnormal to DBL_MAX; k = 1 is LN_2 exactly) | 1,491 | B | 0 ulp |
| `second_order_concentration`, `zero_order_concentration` (NaN propagates, −0.0 kept), `zero_order_half_life` | 987, 502, 262 | B | 0 ulp |
| `steady_state_intermediate`, `steady_state_rate`, `pre_equilibrium_rate`, `michaelis_menten` | 264, 190, 283, 370 | B | 0 ulp |
| `lineweaver_burk` (output `f:slope f:intercept`, the port's caller buffer, port plan §6.4) | 244 | B | 0 ulp |

**Typed constants pinned.** Rust's `std::f64::consts::LN_2` and private `FIRST_ORDER_TOLERANCE` (10⁻¹²) are typed inline
in `src/kinetics.cyr` (port plan §6.4), and `validate.py` requires each to change golden rows when perturbed: LN_2 ± 1 ulp
changes 1,355 / 1,327 of the 1,491 `half_life_first_order` rows and 2,511 / 3,239 of the 7,822 `nth_order_half_life` rows;
the tolerance × ½ / × 2 changes 20 / 26 rows (the orders 1 ± 10⁻¹² and the ulps either side). The tolerance ± 1 ulp changes
none: near 1 ± 10⁻¹², |n − 1| is a multiple of 2⁻⁵² (above 1) or 2⁻⁵³ (below), so only a change across a multiple can
show; `tests/kinetics.tcyr` pins the four boundary orders (1 + 4503·2⁻⁵² rejected, 4504·2⁻⁵² accepted; 1 − 9007·2⁻⁵³ /
9008·2⁻⁵³).

### The class-C functions: measured, not gated

> **Superseded by [ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md) (M3c-3, §10).** The six functions moved to
> `_km_cr_exp`, `_km_cr_expm1` and `_km_cr_pow` and are gated at 0 ulp against the corrected oracle; nothing is measured.
> The S and K tables below are the conditioning that made the rule fail, and the histograms the stdlib primitives'
> distances, kept as history.

The six class-C functions take exp (`f64_exp`), exp_m1 (`_km_expm1`) and powf (`ganita_f64_pow`), as briefed. Their rows
were replayed unbounded (*measured*), then `validate.py --sensitivity` measured S by the rule (every exp, exp_m1 and pow
result nudged ±1 ulp, one call at a time). The grids include **on purpose the rows where an exp or pow result is subnormal
and a prefactor scales it back into the normal range**: A·e^x with x between −745.2 and −708 and A up to DBL_MAX,
[A]₀·e^(−kt), Z·p·e^x, Eyring at low T, either Eyring factor subnormal and the other one huge, and [A]₀^(n−1) in the
nth-order denominator. There one ulp of the subnormal is a large relative change (all of it at 2⁻¹⁰⁷⁴), and S goes far
past the rule's limit:

| Function | Float rows | S (c) | K by the rule | Rows past K = 64 (with a subnormal exp / pow result) | S without those (no subnormal result, no flip) |
|---|---|---|---|---|---|
| `arrhenius_rate` | 9,232 | 4.4 × 10¹⁸ (1) | 8.8 × 10¹⁸ | 858 (858) | 2 (K 6, 8,128 rows) |
| `first_order_concentration` | 2,529 | 4.4 × 10¹⁸ (1) | 8.8 × 10¹⁸ | 855 (855) | 2 (K 6, 1,378 rows) |
| `nth_order_half_life` | 4,605 of 7,822 | 4.5 × 10¹⁵ (2), 10 branch flips | 1.8 × 10¹⁶ | 384 (382; 10 flips) | 4 (K 18, 7,294 rows) |
| `collision_theory_rate` | 2,228 | 4.4 × 10¹⁸ (1) | 8.8 × 10¹⁸ | 391 (391) | 2 (K 6, 1,669 rows) |
| `eyring_rate` | 1,556 | 1.0 × 10¹⁴ (1) | 2.1 × 10¹⁴ | 334 (334) | 2 (K 6, 1,183 rows) |
| `eyring_rate_from_activation` | 2,055 | ∞ (2): 5·10⁻³²⁴ nudged to 0 turns ∞·0 into NaN | ∞ | 676 (676) | 3 (K 14, 1,213 rows) |

(Worst rows: `arrhenius_rate(DBL_MAX, 6195.3555…, 1)` — e^−745.13 = 2⁻¹⁰⁷⁴ times DBL_MAX; `eyring_rate(30976777.65…,
5000)`; `nth_order_half_life(1.409…e−132, 7.381…e268, 3.451…)`, whose [A]₀^(n−1) ≈ 10⁻³²³. `validate.py --sensitivity`
lists them under "functions past K = 64". "Rows past K = 64" are those whose own S puts 2·S·c + 2 past 64; the 2 of
`nth_order_half_life`'s flips without a subnormal result are pow(DBL_MAX, 1) and pow(1, 1)·DBL_MAX at the overflow edge,
exact on both sides.)

**Decision (the user, 2026-11-09): functional parity, correctly rounded primitives.** The port's goal is the oracle's
results, errors and branches, with the best accuracy available — never glibc's rounding errors. So no primitive is built
to reproduce glibc's bits, and **these six functions stay measured, not gated**: the rule forbids a wider K and a thinner
grid, and their bound waits for a later step (ADR 0006) that makes the primitives correctly rounded and compares such rows
against the oracle's formula evaluated with correctly rounded results. The other kinetics functions are gated as usual.
Measured histograms (float outputs; Err rows have none):

| Function | x86_64 ELF and PE | aarch64 |
|---|---|---|
| `arrhenius_rate` | 0: 9,218 · 1: 13 · 2: 1 | 0: 8,922 · 1: 274 · 2: 34 · 3–16: 2 (max 5) |
| `first_order_concentration` | 0: 2,525 · 1: 2 · 2: 1 · 3–16: 1 (max 7) | 0: 2,429 · 1: 80 · 2: 11 · 3–16: 7 · >16: 2 (max 39) |
| `nth_order_half_life` | 0: 4,221 · 1: 308 · 2: 67 · 3–16: 7 · >16: 2 (max 4.9 × 10¹³) | = x86_64 |
| `collision_theory_rate` | 0: 2,223 · 1: 4 · 2: 1 | 0: 2,137 · 1: 75 · 2: 14 · 3–16: 1 · >16: 1 (max 23) |
| `eyring_rate` | 0: 1,551 · 1: 3 · 2: 1 · 3–16: 1 (max 5) | 0: 1,484 · 1: 59 · 2: 10 · 3–16: 3 (max 11) |
| `eyring_rate_from_activation` | 0: 2,006 · 1: 41 · 2: 6 · 3–16: 2 (max 4) | 0: 1,877 · 1: 136 · 2: 36 · 3–16: 6 (max 4) |

- The large distances are the subnormal rows: x87's exp (x86_64) is close to correctly rounded on subnormal results, the
  aarch64 polyfill rounds twice (`y·2^(k+1000)` then `× 2⁻¹⁰⁰⁰`), so aarch64 has more of them. `ganita_f64_pow` and
  `_km_expm1` are + − × ÷ transcriptions, so `nth_order_half_life` is the same on every target.
- **The two `nth_order_half_life` rows past 10¹³ ulp are exact half-way cases, and there the port is right and the oracle
  is not.** (3·2⁻²¹⁵)⁵ and (9·2⁻⁴³⁰)^2.5 are exactly 121.5·2⁻¹⁰⁷⁴: glibc 2.44's pow returns 121·2⁻¹⁰⁷⁴, ganita's the
  correctly rounded (ties-to-even) 122·2⁻¹⁰⁷⁴, and the half-life divides by it. (2⁻¹⁰⁷⁵ from pow(2, −1075), also a tie,
  rounds to 0 on both sides: Err.) A 60-digit decimal reference misjudges exact ties (it rounds pow(0.5, 1075) up to
  2⁻¹⁰⁷⁴); ADR 0006's reference needs exact rational arithmetic there.
- glibc 2.44's exp is correctly rounded on all 6,801 golden calls with x < −708 (2,118 distinct arguments) and on all but
  2 of 200,000 seeded arguments between −745.2 and −708.39; its pow on all 1,845 golden calls with a result under 2⁻¹⁰²¹.
  Recorded for ADR 0006; under the decision above it decides nothing. (MPFR corrects the pow figure, §9: the two exact
  ties are misrounded.)

Summary, every target: **12 modules, 0 failing; 119,044 ok, 0 failed, 22,205 measured (the six kinetics class-C
functions), 10 pending (the M4 inorganic constants), 441 accepted (ADR 0004's row and the 440 ADR 0002 #4 mapped Antoine
rows).**

## 9. M3c-1, ADR 0006: the corrected oracle (2026-11-09)

[ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md) makes correct rounding the accuracy target. A golden
row's reference is now the **corrected oracle**: Rust 1.2.0's operations as `validate.py` transcribes them, evaluated with
every transcendental primitive result correctly rounded by MPFR 4.2.2 (`scripts/oracle/mpfr_ref.py`: 53 bits, binary64's
exponent range, subnormalized, so exact ties and subnormal results are exact). Where glibc rounded every call of a row
correctly, that is the golden value. Where it didn't and the outputs change, `validate.py . --write-corrections` writes
the corrected outputs to `tests/golden/<module>.cr.tsv` (format: `scripts/oracle/README.md`), and the harness compares the
port against them. This step changes the reference and the harness only: the modules still call the primitives of §2–§8.
The correctly rounded primitives and the switch to them are M3c-2 / M3c-3.

**The harness rule** (`tests/parity_core.cyr`), for a row with a correction, under its class and its current bound:
- it passes against the corrected outputs: **accepted (ADR 6)**, counted quietly;
- it fails there but passes against the golden (misrounded) outputs: **measured, glibc-matching**. Transitional, M3c-1
  only: the port still uses a glibc-like primitive on that path, and M3c-3 makes this a failure;
- it matches neither: FAIL. A correction that no golden row takes fails its module too.

Rows without a correction are compared as before. The histograms are the distance from each row's reference (the corrected
outputs where there are some), and so are the shadow lines.

**The report format.** A function line gains `N accepted (ADR 6: a)`, `N measured (g match glibc's misrounding)` and, at
its end, `| ADR 6 rows: n, port = corrected c, = glibc's g` (the rows bit-identical to the corrected outputs and to the
golden ones). A module with a `.cr.tsv` gains a line `ADR 0006 corrected oracle, tests/golden/<m>.cr.tsv: n rows: …`
with every outcome. The summary line keeps its prefix and appends the two counts:

```text
=== parity summary: M modules, F failing | rows: O ok, X failed, Me measured, P pending, A accepted (ADR 6: a), g glibc-matching (measured) ===
```

`A` includes the `a` ADR 6 rows and `Me` the `g` glibc-matching rows, so a reader of the old line still gets the totals.
(Since M3c-3, §10, a glibc-matching row is a failure: the line ends `g glibc-matching (failed)` and `X` includes `g`.)

### The corrections on the current corpus

665 rows have a correction: `gas` 11, `kinetics` 323, `math` 96, `phase` 4, `reaction` 2, `solution` 229. None changes
shape (Ok ↔ Err, the Err variant, the root count); `validate.py` would follow and report such a row. The outcome is the
same on x86_64 ELF, x86_64 PE and aarch64. Where the bit-identical counts differ, aarch64's follow the `·`.

| Function | Rows | glibc misrounded | Outcome | Port = corrected | Port = glibc's |
|---|---|---|---|---|---|
| `f64::exp_m1` (`_km_expm1`) | 34 | expm1 | accepted (ADR 6) | 0 | 34 |
| `f64::log10` (`_km_log10`) | 48 | log10 | accepted | 48 | 0 |
| `f64::cos` (`f64_cos`) | 4 | cos | accepted | 1 | 3 |
| `f64::ln` (`f64_ln`) | 3 | log | accepted | 3 · 2 | 0 · 1 |
| `f64::powf` (`ganita_f64_pow`) | 3 | pow (2 of them x^2.5) | accepted | 1 | 2 |
| `f64::sin` (`f64_sin`) | 2 | sin | accepted | 1 | 1 |
| `f64::acos` (`ganita_f64_acos`) | 1 | acos | accepted | 0 | 1 |
| `f64::powf_base2` (`ganita_f64_pow(2, ·)`) | 1 | exp2 | accepted | 0 | 1 |
| `reaction::equilibrium_constant` | 2 | exp | accepted | 1 | 1 |
| `phase::antoine_pressure_mmhg`, `_pa` | 2 + 2 | pow(10, ·) | accepted | 1 + 1 | 1 + 1 |
| `gas::fugacity_coefficient_pr` | 6 | exp 5, log 1 | accepted | 6 · 5 | 0 · 1 |
| `gas::fugacity_pr` | 5 | exp 4, log 1 | accepted | 5 · 4 | 0 · 1 |
| `solution::ph_from_h_concentration`, `poh_from_oh_concentration` | 2 + 2 | log10 | accepted | 2 + 2 | 0 |
| `solution::weak_acid_ph` | 178 | log10 | accepted | 178 | 0 |
| `solution::weak_base_ph` | 39 | log10 | accepted | 39 | 0 |
| `solution::debye_huckel_extended` | 1 | pow(10, ·) | accepted | 1 | 0 |
| `solution::henderson_hasselbalch` | 7 | log10 | **glibc-matching (measured)** | 0 | 7 |
| `kinetics::nth_order_half_life` | 320 | expm1 318, pow 1, both 1 | measured (no bound) | 3 | 270 |
| `kinetics::eyring_rate` | 1 | exp | measured | 1 · 0 | 0 · 1 |
| `kinetics::eyring_rate_from_activation` | 2 | exp | measured | 1 · 0 | 1 · 2 |

So **335 rows are accepted under ADR 6, 7 are glibc-matching, and the 323 kinetics rows stay measured** (with the rest of
those six functions' rows, §8). 389 more rows make a misrounded call whose error doesn't reach the outputs:
`peng_robinson_z` 5, `peng_robinson_z_roots` 5, `fugacity_coefficient_pr` 9, `fugacity_pr` 10, `henderson_hasselbalch`
18, `weak_base_ph` 111 (14 − pOH absorbs it), `nth_order_half_life` 231.

- **"Accepted" is within the row's current bound, not bit-identical.** Where the port's primitive is glibc-like it
  passes with glibc's misrounded value: `f64::exp_m1` on all 34 rows (`_km_expm1` is fdlibm, which glibc's expm1 equals bit
  for bit), and 1–3 rows each of cos, ln (aarch64), powf, sin, acos, exp2, the equilibrium constant, Antoine and the
  fugacities on aarch64. At ADR 0006 §3's bound of 0 (M3c-3) these would fail; they are what M3c-2's primitives must fix.
- **ADR 0005's correctly rounded primitives.** glibc's cos misrounds 5 distinct golden arguments (20 calls), and the
  Peng–Robinson Newton polish absorbs every one: no root changes. Its log misrounds 5 (10 calls), and one of them (ln of
  1.0198…) changes both of that state's fugacity rows (accepted, ADR 6); acos and x^2.5 are correctly rounded on every
  golden call. The other 9 fugacity corrections are glibc's exp (1 ulp), where the port's exp (x87) is correctly rounded on
  x86_64; the aarch64 polyfill keeps glibc's value on one row of each function.
- **`henderson_hasselbalch`'s 7 rows** are where glibc's log10 misrounds base/acid: at 0.3004… under 5 near-cancellation
  rows (pKa = −log₁₀(base/acid) to within 2 ulp), and at 0.0660… and 8.3 × 10⁻⁶. The correctly rounded log10 is one ulp
  lower, so the 5 pH values each move down one step of 2⁻⁵³ (the exactly-zero row becomes −2⁻⁵³, the +2⁻⁵³ row 0).
  `_km_log10_ln` reproduces glibc's log10, so the port is the oracle's value there: up to 4.4 × 10¹⁸ ulp from the
  corrected one (0 against −2⁻⁵³). ADR 0006 replaces it with a correctly rounded log10 in M3c-3.
- **`nth_order_half_life`'s two exact half-way rows** ((3·2⁻²¹⁵)⁵ and (9·2⁻⁴³⁰)^2.5, §8), 4.9 × 10¹³ ulp from the oracle,
  are 0 ulp from the corrected oracle: the function's histogram maximum falls from 4.9 × 10¹³ to 13 ulp on every target.
- **The shadow lines** now show the correctly rounded primitives at 0 ulp on every `math` row (`_km_cos` 2,660,
  `_km_acos` 2,620, `_km_ln` 3,857, `_km_pow_2_5` 1,100) and `_km_log10_ln` 1 ulp off on exactly glibc's 48 log10 rows.
- **cbrt**: the oracle's (Rust's compiler-builtins, CORE-MATH) is correctly rounded on all 4,877 distinct golden arguments
  (the `math` rows and the Peng–Robinson cubic), so cbrt corrects nothing. glibc's own cbrt, which the oracle never calls,
  misrounds 2,445 of them (up to 3 ulp; about 45 % of random arguments, cbrt(2) included). ganita's is bit-identical to
  the oracle on every `math` row (class A).

**glibc 2.44's misroundings over the golden calls** (`validate.py`'s report; a call of the M2–M3b replicas is one
primitive evaluation inside a function's row, and the `math` rows are the primitives themselves):

| Primitive | `math` rows (misrounded) | Replica calls (distinct arguments) | Misrounded calls (distinct) |
|---|---|---|---|
| exp | 1,125 (0) | 29,966 (11,257) | 15 (10) |
| expm1 | 810 (34) | 7,741 (1,267) | 550 (109) |
| log | 3,857 (3) | 19,572 (8,877) | 10 (5) |
| log10 | 1,296 (48) | 24,603 (11,027) | 357 (215) |
| pow(x, y) | 1,750 (3) | 7,364 (1,757) | 2 (2: the exact ties) |
| pow(10, y) | 571 (0) | 3,239 (1,906) | 5 (3) |
| pow(x, 2.5) | (in pow's 1,750: 1,100 rows, 2 misrounded) | 593 (13) | 0 |
| exp2 (`powf_base2`) | 369 (1) | — | — |
| cos | 2,910 (4) | 18,444 (4,209) | 20 (5) |
| acos | 2,620 (1) | 6,148 (1,403) | 0 |
| sin | 710 (2) | — | — |
| atan2 | 1,000 (0) | — | — |

Every misrounding is 1 ulp. The `math` counts are the 60-digit `decimal` references' (§2, §4, §5), now decided by MPFR.
One earlier figure changes: §8's "glibc's pow is correctly rounded on all 1,845 golden calls with a result under 2⁻¹⁰²¹"
used the decimal reference, which misjudges ties. With MPFR it misrounds 2 of the 662 distinct arguments (1,851 calls),
the two ties above. Its exp is correctly rounded on all 6,802 calls with x < −708 (2,115 distinct arguments), as §8
says.

**`--sensitivity`** now decides "glibc is not correctly rounded here" with MPFR. Every S, c, K, worst row and count it
prints is unchanged, so no gated K moves.

Summary, every target: **12 modules, 0 failing; 118,702 ok, 0 failed, 22,212 measured (the six kinetics class-C functions
and 7 glibc-matching rows), 10 pending (the M4 inorganic constants), 776 accepted (335 under ADR 6, ADR 0004's row and the
440 ADR 0002 #4 mapped Antoine rows).**

## 10. M3c: every class-C row exact (2026-11-09)

[ADR 0006](../adr/0006-correct-rounding-is-the-accuracy-target.md), completed. M3c-1 made the reference the corrected oracle
(§9); M3c-2 built the correctly rounded primitives of `src/mathx.cyr` — `_km_cr_exp`, `_km_cr_expm1`, `_km_cr_exp10`,
`_km_cr_ln`, `_km_cr_log10` and `_km_cr_pow`, two-phase (Ziv: a fast path with a proven error bound and a rounding test,
a triple-double accurate path), on the MPFR-generated tables of `src/mathx_data.cyr` (`scripts/oracle/mathx_tables.py`),
verified against MPFR by `scripts/oracle/cr_verify.py` and `tests/mathx.tcyr`; M3c-3 moved every module onto them and
compares every class-C row exactly. **Every class-C function the port has is at bound 0, and every row is at 0 ulp from
its reference, on x86_64 ELF, x86_64 PE and aarch64.** Nothing is measured; no row is glibc-matching.

### The switch

Rust's operation order is kept everywhere; only the primitive changed.

| Module (functions) | Rust | Until M3c-3 | Now |
|---|---|---|---|
| `reaction::equilibrium_constant`, `phase::clausius_clapeyron`, `thermochem::vant_hoff_k`, `gas::fugacity_coefficient_pr` / `fugacity_pr`, the five exp functions of `kinetics` | `exp` | `f64_exp` (x87 on x86_64, the fdlibm polyfill on aarch64) | `_km_cr_exp` |
| `kinetics::nth_order_half_life` | `exp_m1`, `powf` | `_km_expm1` (fdlibm), `ganita_f64_pow` | `_km_cr_expm1`, `_km_cr_pow` |
| `phase::antoine_pressure_mmhg` / `_pa`, `solution::h_concentration_from_ph`, `debye_huckel_limiting`, `debye_huckel_extended`, `davies_activity` | `10.0_f64.powf(x)` | `ganita_f64_pow(10, x)` | `_km_cr_exp10(x)` (the same value as a correctly rounded pow(10, x)) |
| `solution::ph_from_h_concentration`, `poh_from_oh_concentration`, `weak_acid_ph`, `weak_base_ph` | `log10` | `_km_log10` (fdlibm) | `_km_cr_log10` |
| `solution::henderson_hasselbalch` | `log10` | `_km_log10_ln` (glibc's composition, ADR 0005 amendment M3a) | `_km_cr_log10` |
| `electrochemistry::nernst_potential` / `_25c`, `gas::fugacity_coefficient_pr` / `fugacity_pr` | `ln` | `_km_ln` (since M3c-2a a wrapper of `_km_cr_ln`) | `_km_cr_ln` |
| `gas::redlich_kwong_pressure` | `tc.powf(2.5)` | `_km_pow_2_5` (double-double, ADR 0005) | `_km_cr_pow(tc, 2.5)` |
| `gas::peng_robinson_z` / `_roots` (and the fugacities' cubic) | `acos`, `cos` | `_km_acos`, `_km_cos` (ADR 0005) | unchanged |

- **x^2.5.** `_km_pow_2_5` was correctly rounded only from 2⁻³⁸⁷ (below, its double-double product's error term underflows)
  and had no rounding test. Measured per call: `_km_pow_2_5` 85 ns, `_km_cr_pow(x, 2.5)` 130 ns (x = 126.2 and 647.1;
  `ganita_f64_pow` 185 ns). `redlich_kwong_pressure` now calls `_km_cr_pow(tc, 2.5)` on its whole range — one primitive
  with a proven bound, a rounding test and the exact / midpoint cases — and `_km_pow_2_5` is gone; the function is not
  benched in Rust, and its 45 ns are the price.
- **Removed, nothing calls them:** `_km_log10` and `_km_expm1` (M1's Rust libm transcriptions), `_km_log10_ln` (M3a), the
  `_km_ln` wrapper, `_km_pow_2_5` and the unused `_km_dd_div`. `phase`, `solution` and `kinetics` no longer include
  `lib/ganita.cyr`; `gas` does, for `ganita_f64_cbrt` (the oracle's cbrt, compiler-builtins', and ganita's are both
  correctly rounded on every golden argument, §9). `dist/kimiya.deps` is unchanged (math, string, vec, str, ganita,
  tagged): `mathx.cyr` keeps `lib/math.cyr` for `_km_cos`'s fallback.
- The Cyrius-only test groups that pinned glibc's bits now pin the correctly rounded ones: `tests/solution.tcyr`'s
  Henderson–Hasselbalch row (base/acid = 0.8068…/2.6858…) asserts log10 = `…69DB` (glibc `…69DA`) and the corrected
  oracle's pH −2⁻⁵³ (the oracle's +0.0). No ported Rust test changed.

### The harness

- **Bound 0** in `_par_c_bound` for all 27 kimiya class-C functions and for the `math` primitives `f64::exp`, `exp_m1`,
  `ln`, `log10`, `powf`, `powf_base10`, `powf_base2` and `acos`, each now dispatched to the port's correctly rounded
  primitive (`powf_base10` to `_km_cr_exp10`, what the modules call; `powf_base2`, glibc's exp2 in the oracle, to
  `_km_cr_pow(2, y)`); `f64::cos` goes to `_km_cos` at bound 0 where `_km_cos` is correctly rounded (|x| ≤ 8, and NaN / ±∞).
- **What stays on the stdlib, at its bound of 2:** `f64::sin` (`f64_sin`) and `f64::atan2` (`ganita_f64_atan2`), which have
  no correctly rounded port primitive yet — M4 adds them with their first callers (the dihedral, spectroscopy; ADR 0006
  item 1) — and the **247 `f64::cos` rows with a finite |x| > 8**, outside `_km_cos`'s correctly rounded range, where it is
  `f64_cos` (6 of them 1 ulp from the correctly rounded value). No kimiya function calls any of them on a golden row: the
  Peng–Robinson cos arguments are in [−2π/3, π].
- **glibc-matching is a failure.** A row whose port value fails against the corrected outputs but matches the oracle's
  misrounded ones fails its module (printed as `FAIL: (matches only the oracle's misrounded outputs, ADR 0006)`) and is
  counted on its own as well; the summary line ends `g glibc-matching (failed)` and `X failed` includes them. 0 on every
  target.
- **ADR 0004's row**, `f64::powf(−1, DBL_MAX)`, goes through `_km_cr_pow`, which gives 1: a plain ok, and its
  `_par_accepted` entry is gone. ADR 0004 stays the record of ganita's defect, and so does the new shadow line
  `f64::powf ~ ganita_f64_pow` (its one nan/inf mismatch).
- **Shadow lines** (never gated) now show the stdlib primitives the modules stopped calling, against the corrected oracle:
  `f64_exp` 45 rows at 1 ulp (aarch64: 99), `f64_ln` 0 (aarch64: 152), `ganita_f64_acos` 727 at 1 ulp and 1 at 2 (aarch64:
  743 and 2), `ganita_f64_pow` 153 at 1 ulp and ADR 0004's row; `_km_cr_pow(10, y)` is bit-identical to `_km_cr_exp10` on all
  571 `powf_base10` rows; the cos rows past 8 on their own line. The ADR 0005 shadows (`_km_cos`, `_km_acos`, `_km_ln`,
  `_km_pow_2_5`, `_km_log10_ln`) and the polyfill shadows are retired: those rows are the gated ones now.

### Every class-C row, all three targets

Bound 0 everywhere but the four `math` lines marked. "ADR 6" counts the rows compared against `tests/golden/<m>.cr.tsv`
(accepted under ADR 6); "floats" the compared float outputs (Err rows have none; the roots count each root). Every
function's histogram is the same on x86_64 ELF, x86_64 PE and aarch64, except `f64::atan2` on aarch64 (0: 798 · 1: 202).

| Function | Rows | Floats | ADR 6 | Primitives | ulp (every target) |
|---|---|---|---|---|---|
| `reaction::equilibrium_constant` | 1,083 | 1,030 | 2 | exp | 0 |
| `phase::clausius_clapeyron` | 849 | 815 | 0 | exp | 0 |
| `phase::antoine_pressure_mmhg` / `_pa` | 558 each (222 mapped, ADR 2) | 338 each | 2 each | exp10 | 0 |
| `gas::peng_robinson_z` | 4,663 | 4,603 | 0 | acos, cos | 0 |
| `gas::peng_robinson_z_roots` | 4,663 | 5,539 roots | 0 | acos, cos | 0 |
| `gas::fugacity_coefficient_pr` | 4,663 | 4,544 | 6 | acos, cos, ln, exp | 0 |
| `gas::fugacity_pr` | 4,663 | 4,544 | 5 | acos, cos, ln, exp | 0 |
| `gas::redlich_kwong_pressure` | 623 | 547 | 0 | pow(x, 2.5) | 0 |
| `thermochem::vant_hoff_k` | 1,055 | 1,019 | 0 | exp | 0 |
| `solution::ph_from_h_concentration` / `poh_from_oh_concentration` | 2,851 each | 2,838 each | 2 each | log10 | 0 |
| `solution::henderson_hasselbalch` | 859 | 831 | 7 | log10 | 0 |
| `solution::h_concentration_from_ph` | 439 | 439 | 0 | exp10 | 0 |
| `solution::weak_acid_ph` | 9,111 | 9,048 | 178 | log10 | 0 |
| `solution::weak_base_ph` | 9,111 | 9,048 | 39 | log10 | 0 |
| `solution::debye_huckel_limiting` / `davies_activity` | 503 each | 488 each | 0 | exp10 | 0 |
| `solution::debye_huckel_extended` | 538 | 505 | 1 | exp10 | 0 |
| `electrochemistry::nernst_potential` | 1,485 | 1,430 | 0 | ln | 0 |
| `electrochemistry::nernst_potential_25c` | 1,132 | 1,106 | 0 | ln | 0 |
| `kinetics::arrhenius_rate` | 9,256 | 9,232 | 0 | exp | 0 |
| `kinetics::first_order_concentration` | 2,529 | 2,529 | 0 | exp | 0 |
| `kinetics::nth_order_half_life` | 7,822 | 4,605 | 320 | expm1, pow | 0 |
| `kinetics::collision_theory_rate` | 2,297 | 2,228 | 0 | exp | 0 |
| `kinetics::eyring_rate` | 1,573 | 1,556 | 1 | exp | 0 |
| `kinetics::eyring_rate_from_activation` | 2,076 | 2,055 | 2 | exp ×2 | 0 |
| `f64::exp` | 1,125 | 1,125 | 0 | `_km_cr_exp` | 0 |
| `f64::exp_m1` | 810 | 810 | 34 | `_km_cr_expm1` | 0 |
| `f64::ln` | 3,857 | 3,857 | 3 | `_km_cr_ln` | 0 |
| `f64::log10` | 1,296 | 1,296 | 48 | `_km_cr_log10` | 0 |
| `f64::powf` | 1,750 | 1,750 | 3 | `_km_cr_pow` (ADR 0004's row ok) | 0 |
| `f64::powf_base10` | 571 | 571 | 0 | `_km_cr_exp10` | 0 |
| `f64::powf_base2` | 369 | 369 | 1 | `_km_cr_pow(2, y)` | 0 |
| `f64::acos` | 2,620 | 2,620 | 1 | `_km_acos` | 0 |
| `f64::cos`, \|x\| ≤ 8 or non-finite | 2,663 | 2,663 | 4 | `_km_cos` | 0 |
| `f64::cos`, finite \|x\| > 8 (**bound 2**) | 247 | 247 | 0 | `f64_cos` | 0: 241 · 1: 6 |
| `f64::sin` (**bound 2**) | 710 | 710 | 2 | `f64_sin` | 0: 690 · 1: 20 |
| `f64::atan2` (**bound 2**) | 1,000 | 1,000 | 0 | `ganita_f64_atan2` | 0: 797 · 1: 203 |

The 665 ADR 6 rows are the corrected rows of §9 (no new correction, below); on every one the port is bit-identical to the
corrected outputs except `f64::sin`'s, where `f64_sin` keeps glibc's value on 1 of its 2 rows (within its bound 2).

### The hard rows (ADR 0006 item 6)

The M2 functions gain the subnormal-exp-times-prefactor inputs M3b found outside their grids, hand-picked on targets
plus a seeded sweep (`scripts/oracle/src/m_*.rs`, fresh seeds, so no existing row moved; `grids::EXP_SUBNORMAL_X`, 14
exponents from −745.2 to −708.39, and `grids::subnormal_prefactors`, DBL_MAX, 10³⁰⁰, 2¹⁰⁰⁰, 10²⁵⁰ and 1; the exponent placed
with + − × ÷ only):

| Function | New rows | Hand-picked | Seeded | New outputs: normal · subnormal · 0 · Err |
|---|---|---|---|---|
| `reaction::equilibrium_constant` (K is the exp result) | 370 | 70: 14 exponents at T = 2, 298.15, 1500, 5000 and 10⁻³ K | 300: x in [−745.2, −708], T 10⁻¹ … 10⁴ K | 13 · 347 · 10 · 0 |
| `thermochem::vant_hoff_k` | 370 | 70: 14 exponents × 5 K₁ at 300 → 299 K | 300: K₁ 10²⁰⁰ … 10³⁰⁸, T₁, T₂ 100–2000 K | 348 · 11 · 11 · 0 |
| `phase::clausius_clapeyron` | 370 | 70: 14 exponents × 5 P₁ at 373.15 → 363.15 K | 300: P₁ 10²⁰⁰ … 10³⁰⁸, T₁, T₂ 150–700 K | 349 · 11 · 10 · 0 |
| `gas::fugacity_coefficient_pr` / `fugacity_pr` (and the state's Z and roots rows) | 344 states × 4 | 44: water 9.76–10.25 K, CO₂ 3.9–4.05 K, methane 1.7–1.75 K, P 10³ … 10⁷ Pa (water at 10.25 K and 5·10⁵ Pa: φ ≈ 2.4·10⁻³⁰⁹) | 300: 100 per gas over the same boxes, P 10³ … 10^7.5 Pa | φ: 51 · 245 · 0 · 48; f = φ·P: 133 · 163 · 0 · 48 |

(The 48 Err rows are states with ln φ below −745.13, where φ rounds to 0: `ComputationError`, on both sides.)

- `validate.py` recomputes every new row bit for bit from the transcription and through the corrected oracle. **None needs
  a correction**: glibc's exp and log are correctly rounded on every new call (exp: 31,764 golden calls now, 12,506
  distinct arguments, still the same 15 misrounded calls at 10 arguments; log 20,948 calls, 9,565 arguments, the same 10 at
  5). So these rows test the port's primitive, not glibc's: they are where a double-rounding exp shows.
- **What they showed.** With the stdlib's `f64_exp` in place of `_km_cr_exp` in the four modules (a scratch build, not
  committed), the functions are at most 1 ulp from the oracle on x86_64 (x87's exp is close to correctly rounded on
  subnormal results), but on aarch64 the polyfill rounds twice (`y·2^(k+1000)`, then `× 2⁻¹⁰⁰⁰`) and they reach
  **27 ulp** (`vant_hoff_k`), **15** (`fugacity_pr`) and **9** (`clausius_clapeyron`); `equilibrium_constant`, whose output
  is the exp result itself, at most 1. With `_km_cr_exp` every row is 0 ulp on every target.
- **Their conditioning** (`--sensitivity`, below): `vant_hoff_k` and `clausius_clapeyron` gain 306 and 303 rows past
  K = 64 (S up to 4.4 × 10¹⁸, all with a subnormal exp result), `fugacity_pr` 245 (S up to 1.6 × 10⁷) and
  `fugacity_coefficient_pr` S = 1,539 at the low-temperature states, where |ln φ| is about 700. Under the old rule none of
  them could have been gated; under ADR 0006 they are gated at 0.

### validate.py and `--sensitivity`

- The port model is the oracle's own code with correctly rounded primitives, which is the corrected oracle, so
  `SENS_RULES` (the `_km_log10_ln` model), the "S before" mode and the e_log10.c composition check are gone.
- `--sensitivity` is a diagnostic of conditioning (ADR 0006 item 3): it nudges **every** primitive result (but 0, ±∞, NaN
  and acos(±1)) one ulp, and prints S, c, the retired rule's K for comparison, and the calls where glibc is not correctly
  rounded. It reproduces the conditioning §5 and §7 measured "before": `peng_robinson_z_roots` S = 2²⁶,
  `redlich_kwong_pressure` 64, `nernst_potential` 8.8 × 10¹⁸, `henderson_hasselbalch` 4.4 × 10¹⁸, and §8's kinetics values.

### The corpus and the summary

Golden rows: `gas` 21,093 (+1,376), `phase` 3,646 (+370), `reaction` 1,558 (+370), `thermochem` 9,857 (+370); the
corrections are unchanged (665 rows: `gas` 11, `kinetics` 323, `math` 96, `phase` 4, `reaction` 2, `solution` 229), and a
second `--write-corrections` is byte-identical.

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **12 modules, 0 failing; 143,071 ok, 0 failed, 0 measured,
10 pending (the M4 inorganic constants), 1,105 accepted (665 under ADR 6 and the 440 ADR 0002 #4 mapped Antoine rows),
0 glibc-matching.**

```text
=== parity summary: 12 modules, 0 failing | rows: 143071 ok, 0 failed, 0 measured, 10 pending, 1105 accepted (ADR 6: 665), 0 glibc-matching (failed) ===
```

## 11. M4b: organic, inorganic (2026-11-09)

Two new golden modules, every row 0 ulp from its reference on x86_64 ELF, x86_64 PE and aarch64 (the oracle's value, or
for 13 CFSE rows the exact value of ADR 0002 #12); neither module calls a
transcendental (inorganic's one primitive is √, IEEE), so no row has an ADR 0006 correction (`validate.py
--write-corrections` writes no new `.cr.tsv`, and the six existing ones are byte-identical).

| Function | Class | Rows | Grid | ulp (every target) |
|---|---|---|---|---|
| `BondType::energy_kj` | A | 45 | every variant | 0 |
| `Geometry::ideal_bond_angle` | A | 12 | every variant (`TShapedThreeLonePairs` → `KM_GEOMETRY_T_SHAPED`, ADR 0002 #6) | 0 |
| `FunctionalGroup::pattern` / `is_polar` / `hydrogen_bonding` | A | 16 each | every variant (`hydrogen_bonding` → `functional_group_is_hydrogen_bond_donor`, ADR 0002 #6) | — (text, flags) |
| `organic::predict_geometry` | A | 7,491 | the 0 … 15 square and every u8 of each argument against the other's edges (0 … 8, 15, 16, 127, 128, 254, 255) | — (`e:` / `n:`) |
| `organic::enthalpy_from_bonds` | B | 1,508 | empty lists, the Rust and integration tests' lists, every bond alone on each side × 9 counts to u32::MAX, self-cancellation, all 45 at once, 600 seeded lists | 0 |
| `inorganic::cfse_octahedral` | A | 512 | every u8 × both spins (Ok 0 … 10, Err beyond); **ADR 0002 #12**: 499 ok, 13 accepted (ADR 2) by name, where the port's exact (3·e_g − 2·t₂g)/5 differs from the oracle's | 0 from the exact value |
| `inorganic::unpaired_electrons_oct` | A | 512 | every u8 × both spins | — (integers) |
| `inorganic::delta_tet_from_oct` | B | 1,644 | the edges, two mantissas over every decade, DBL_MAX / 4 ± 1 ulp, 300 seeded; since M4c's follow-up 2¹⁰²² and its neighbours and 50 seeded Δ in [2¹⁰²², 2¹⁰²⁴) (**ADR 0002 #15**: 59 rows where Rust's Δ·4 overflows, mapped to the port's correctly rounded 4Δ/9 and accepted (ADR 2) by name) | 0 from CR(4Δ/9) |
| `Ligand::field_strength` / `is_strong_field` | A | 11 each | every variant | — |
| `inorganic::spin_only_moment` | B | 256 | every u8 | 0 |
| `inorganic::born_lande_lattice_energy` | B | 3,369 | each argument over the edges, both charges over 11 u32 probes, the Born exponent at 1 and the ulps above, every Madelung constant × Born exponent and mean × charges 1–3 × r₀ 1.5–4 Å, r₀ and M over every tenth decade (−∞ overflow, −0.0 underflow), 600 seeded salts | 0 |
| `inorganic::born_haber_lattice_energy` | B | 496 | each argument over the edges, the Rust tests' NaCl cycle, ∞ − ∞, signed zeros, 400 seeded | 0 |
| `inorganic::IONIC_RADII` / `IONIC_RADII[]` | A | 1 + 30 | every row | 0 |
| `inorganic::lookup_ionic_radius` | A | 8,320 | every table symbol × every i8; case, space and prefix variants, invalid symbols × −4 … 4 | 0 |

- **The constants resolved.** The ten `MADELUNG_*` / `BORN_EXP_*` are generated into `src/inorganic_data.cyr`, so the last
  10 pending `constants` rows are compared (class A, bit-identical): **0 pending** for the first time.
- **Class B on integers.** Every bond energy is an integer and every count a u32, so `enthalpy_from_bonds`'s terms and
  sums are exact (below 2⁵³): its rows pin the record handling, the −0.0 folds and broken − formed, not rounding.
- **ADR 0002 #12, the exact CFSE.** Rust's `t2g·(−0.4) + eg·0.6` is not exact: high-spin d⁵ is 3·(−0.4) + 2·0.6 = −2⁻⁵²
  (3 × 0.4 is a tie that rounds up), d¹⁰ −2⁻⁵¹, and 11 more cases miss the double nearest k/5 by 1–4 ulp. The port returns
  CR((3·e_g − 2·t₂g)/5), one division of an exact integer. The 13 rows keep the oracle's value in the corpus, with the
  exact value in their note; the harness recovers q = round(5 × oracle) (within 10⁻⁹), requires the port's value to be
  q/5 bit for bit, and prints each row as `accepted (ADR 2)`; the rows are class A, so no tolerance lets a wrong value
  through (a port computing q × 0.2 fails 5 of them, tried). The 9 rows where Rust's value is already exact are plain ok.
  The mapped rows' 0-ulp distance from the exact value is kept out of the histogram, which shows the 9.
- **Other formulas, reviewed under the same policy** (port plan §6.4): `delta_tet_from_oct` is already correctly rounded
  wherever its result is finite (Δ·4 is exact; 1,580 rows checked against exact rational arithmetic), and differs from
  CR(4Δ/9) only on the 5 rows where Δ·4 overflows — decided in M4c's follow-up: ADR 0002 #15 takes (Δ/9)·4 there, the
  correctly rounded 4Δ/9 (`validate.py` proves the port's value on all 1,644 rows and 40,000 seeded ones with `Fraction`;
  §12); Born–Landé's (1 − 1/n) differs from CR((n − 1)/n) on 332 of the 3,263 Ok
  rows (n = 7 among them), but the energy is a product of rounded constants with no exact value; it is left as Rust's,
  pending the user's call. The spin-only moment is correctly rounded already (n(n + 2) is an exact integer).
- **ADR 0002 #6.** The harness compares the renamed items' rows under the Rust names (16 `hydrogen_bonding` rows through
  the donor function, 2 rows naming `Geometry.TShapedThreeLonePairs` through `KM_GEOMETRY_T_SHAPED`) as plain class A
  and prints the count after the organic module. The new acceptor property (IUPAC 2011's criterion, weak acceptors
  included: every group but the alkane) has no golden rows (reference tests only); `validate.py` recomputes its column
  from the stated rule and `probe_tables.py` compiles it against the same rule.

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **14 modules, 0 failing; 167,279 ok, 0 failed, 0 measured,
0 pending, 1,118 accepted (665 under ADR 6; under ADR 2, the 440 ADR 0002 #4 mapped Antoine rows and the 13 ADR 0002 #12
CFSE rows), 0 glibc-matching.**

```text
=== parity summary: 14 modules, 0 failing | rows: 167279 ok, 0 failed, 0 measured, 0 pending, 1118 accepted (ADR 6: 665), 0 glibc-matching (failed) ===
```

## 12. M4c: nuclear (2026-11-09)

One new golden module, 54,057 rows, every row 0 ulp from its reference on x86_64 ELF, x86_64 PE and aarch64: the oracle's
value, the corrected oracle's on 23 rows (ADR 6), the port's `Err` on 23 SEMF rows where the oracle returns NaN (ADR 0002
#3), or the exact value the note names on 1,812 `atoms_after_half_lives` and 208 `is_secular_equilibrium` rows (ADR 0002
#13, #14; with inorganic's 59 `delta_tet_from_oct` rows of #15, the user's decision of 2026-11-09 on M4c's review). The four class-C functions are on correctly rounded primitives (exp, expm1, ln and pow(2, ·)), in Rust's
operation order, at **bound 0**; the class-B ones use + − × ÷, √ and cbrt only.

| Function | Class | Rows | Grid | ADR | ulp (every target) |
|---|---|---|---|---|---|
| `nuclear::decay_constant` | B | 2,237 | the edges, every isotope, every decade from the smallest subnormal to DBL_MAX (× 1, 3, ln 2), t½ = 1 (LN_2 exactly), +∞ (λ = 0), seeded | — | 0 |
| `nuclear::atoms_remaining` | C (exp) | 1,575 | each argument over the edges, every isotope over 0.5–1000 half-lives, λt on 28 targets across exp's range × 8 N₀, **e^(−λt) subnormal × N₀ up to DBL_MAX** (`EXP_SUBNORMAL_X` × the prefactors; 300 seeded), realistic and wide sweeps | — | 0 |
| `nuclear::activity` | B | 543 | the edges against each other, every isotope, seeded | — | 0 |
| `nuclear::atoms_after_half_lives` | C (pow(2, ·)) | 3,350 | the edges, integer n over −1100 … −1000, −40 … 40, 1000 … 1100 and every 30th between, fractional n, **2⁻ⁿ subnormal × N₀ up to DBL_MAX** (n = 1022 … 1075.5 in quarters × 5 prefactors; 300 seeded), the tie at n = 1075, 2⁻ⁿ's over / underflow, seeded; for #13 N₀ = ±2^j, ±0, ±∞, NaN, exact subnormal ties, and 900 seeded n past both ends | **#13**: 1,812 rows mapped to the noted N₀/2ⁿ (564 integer n and 173 N₀ = ±2^j correctly rounded, 1,069 faithful, 6 scaleB) | 0 (from the oracle or the note) |
| `nuclear::time_for_fraction` | C (ln) | 5,693 | the edges, the Rust tests' Err paths in order, fractions at 1 and its ulps below, every 4th decade from 10⁻³²³ × λ over every 25th decade 10⁻³⁰⁰ … 10³⁰⁰, every isotope, 1 − 10⁻¹⁶ … 1 − 10⁻¹, seeded | — | 0 |
| `nuclear::mass_defect` | B | 594 | 16 u32 probes (0 … u32::MAX) on each side, the edges of the mass, the AME2020 He-4 and Fe-56 rows, seeded nuclei | — | 0 |
| `nuclear::binding_energy_mev` | B | 641 | the edges, every 3rd decade, seeded | — | 0 |
| `nuclear::binding_energy_per_nucleon` | B | 618 | the u32 probes on each side (the `checked_add` overflow, A = 0), the edges of the mass with A = 0 and the overflow first, seeded | — | 0 |
| `nuclear::semi_empirical_binding_energy` | B (cbrt, √) | 30,351 | **every Z ≤ A ≤ 3Z with A ≤ 300** (30,300 pairs, port plan §4), Z = 0, the u32 extremes, the Err region (z > a, a = 0) | **#3**: 23 rows, oracle NaN ↔ `Err(InvalidInput)`, accepted (ADR 2) and printed by name | 0 (30,328 finite rows) |
| `nuclear::q_value` | B | 386 | empty slices and the −0.0 folds, the Rust tests' D–T fusion both ways, each edge alone, ∞ − ∞, overflow, cancellation, U-235 fission, 300 seeded lists of 0–6 masses | — | 0 |
| `nuclear::decay_chain_two_step` | C (exp ×2, expm1 ×1–2) | 6,171 | each argument over the edges, every Err path in Rust's order, 16 λ × 16 λ × 12 t, N₀ over 6 values, λ₂ = λ₁(1 ± r) for r = 10⁻¹⁵ … 10⁻⁶ either side of 10⁻¹², **pairs on the threshold to the ulp of 10⁻¹²**, **e^(−λ₁t) and e^(−λ_min·t) subnormal × N₀ and N₀·λ₁**, every isotope pair at the parent's t½ and a day, realistic and wide sweeps | **ADR 6: 23** (glibc's expm1) | 0 |
| `nuclear::is_secular_equilibrium` | A | 1,698 | the edges against each other, λ_p at 0.01·λ_d and ± 1–2 ulp for 11 λ_d, the Rust tests' values, seeded; for #14 λ_p at RN(λ_d/100) ± 1 ulp over 300 signed magnitudes, the underflow of 0.01·λ_d, the overflow end | **#14**: 208 rows mapped to the exact λ_p < λ_d/100 | — (bool) |
| `nuclear::ISOTOPES` / `ISOTOPES[]` | A | 1 + 18 | every row | — | 0 |
| `nuclear::lookup_isotope` | A | 181 | every label; its lower- / upper-case forms, a leading / trailing space, without the hyphen, every proper prefix; `""`, `Xx-999`, `C-12`, `Tc-99`, `Tc-99M`, `Carbon-14`, … | — | 0 |

By class: A 1,898, B 35,370, C 16,789 rows.

- **The corrected oracle.** glibc's expm1 misrounds 31 of the module's calls; 23 rows change (all `decay_chain_two_step`,
  through N_C = −N₀·expm1(−λ₁t) − N_B or the Bateman growth term) and go to `tests/golden/nuclear.cr.tsv`; 8 absorb the
  error. glibc's exp, log and exp2 are correctly rounded on every nuclear call (exp2: 2,421 calls, 1,458 distinct
  arguments), so the hard rows below test the port's primitives, not glibc's. A second `--write-corrections` leaves the
  six existing `.cr.tsv` files byte-identical. On all 23 the port is bit-identical to the corrected outputs; 0 match
  glibc's.
- **The hard rows** (ADR 0006 item 6): 370 `atoms_remaining` rows, 1,105 `atoms_after_half_lives` rows and 554
  `decay_chain_two_step` rows have a subnormal primitive result scaled back by a prefactor (`--sensitivity`: S up to
  4.4 × 10¹⁸; without them `atoms_remaining` and `atoms_after_half_lives` have S = 2, K = 6 under the retired rule).
  `decay_chain_two_step` has 2,538 rows past K = 64, the 554 subnormal ones and rows where a nudged result overflows
  N₀·λ₁·e^(−λ_min·t) to ∞ (S = ∞), and S = 6 over the other 3,588. All are gated at 0 like every other row; since
  ADR 0002 #13 the `atoms_after_half_lives` ones are compared with the noted N₀/2ⁿ (the port no longer multiplies by a
  subnormal 2⁻ⁿ there), so the subnormal pow they exercised is tested by the `math` rows, `nth_order_half_life` and
  `cr_verify.py`.
- **ADR 0002 #3.** The 23 Err-region rows keep the oracle's NaN with the port's `Err(InvalidInput)` in their note; the
  harness maps NaN ↔ that `Err` (an `Ok` or another code fails), prints each row as `accepted (ADR 2)` and then
  `ADR 0002 #3 SEMF Result: 23 … rows …, 30328 finite rows compared as values (class B)`. `validate.py` checks that the
  note sits on exactly the NaN rows and that those are exactly z > a or a = 0.
- **The equal-rate threshold.** 946 Ok chain rows take Rust's limit branch (360 with λ₁ = λ₂ exactly) and 5,180 the
  Bateman form; on 738 the branch moves when 10⁻¹² moves one ulp either way, so the typed tolerance is pinned to the ulp
  (`validate.py`: 108 and 158 rows change for ± 1 ulp, 198 and 148 for × ½ and × 2).
- **The typed constants** (`validate.py`'s pinning): LN_2 ± 1 ulp changes 2,044 / 2,005 `decay_constant` rows; the SEMF
  coefficients 8,116 – 19,404 rows each (the pairing δ 35 / 39), the asymmetry's 4 16,798. (The secular 0.01 left with
  ADR 0002 #14.)
- **ADR 0002's inexact-formula policy** (port plan §6.4, "Decided in M4c"). The review found `atoms_after_half_lives`
  missing the exact N₀/2ⁿ on 20 rows of the first corpus (all where 2⁻ⁿ over- or underflows), `is_secular_equilibrium`
  disagreeing with the exact λ_p < λ_d/100 on 66, and M4b's Δ_tet overflow; the user decided "fix exactly" (ADR 0002
  #13–#15). `decay_constant`'s LN_2/t½ (1 ulp from CR(ln 2/t½) on 476 of 2,175 rows) keeps Rust's operations: ln 2/t½ is
  irrational.
- **ADR 0002 #13, #14, #15 — the mapping.** A row where the port's value differs from the oracle's keeps the oracle's,
  and its note names the port's (`[ADR 0002 #13: the port returns <how> = <x> (f:<bits>)]`, `… #14 … (i:1)]`), written by
  the generator in integer arithmetic. `validate.py` proves every note with `Fraction` (and MPFR for CR(2⁻ᶠ) and
  CR(2^(j − n))) and checks that the notes sit exactly where the values differ; the harness reads the value from the note,
  requires the port's bit for bit (anything else fails) and prints the row by name as `accepted (ADR 2)`. On every target
  all 1,812 + 208 + 59 match. Integer-n and power-of-two rows are correctly rounded; the 1,069 noted faithful-path rows (of
  1,121 on that path) are within one ulp of N₀/2ⁿ, correctly rounded on 743 and nearer it than Rust's value on 1,087, as
  near on 28, farther on 6 (by one ulp: N₀·CR(2⁻ⁿ) with a 51–52-bit subnormal 2⁻ⁿ landing right by chance, at n = 1022.25 …
  1023.5); a 4,000-pair sweep beyond the corpus is correctly rounded on 82 % and faithful on all. (M4d made this path
  correctly rounded too, `_km_cr_exp2_scaled`: §14. The figures here are M4c's.) #14's comparison is
  replayed against `Fraction` on 99,940 seeded boundary pairs, #15's (Δ/9)·4 on 20,000 Δ in [2¹⁰²², 2¹⁰²⁴), with
  Rust's form checked on 20,000 more below.

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **15 modules, 0 failing; 219,266 ok, 0 failed, 0 measured,
0 pending, 3,243 accepted (688 under ADR 6; under ADR 2, the 440 ADR 0002 #4 Antoine rows, the 13 #12 CFSE rows, the 23 #3
SEMF rows, and the 1,812 #13, 208 #14 and 59 #15 exact-value rows), 0 glibc-matching.** (Before the follow-up:
219,444 ok and 1,164 accepted, on 52,211 nuclear and 15,107 inorganic rows.)

```text
=== parity summary: 15 modules, 0 failing | rows: 219266 ok, 0 failed, 0 measured, 0 pending, 3243 accepted (ADR 6: 688), 0 glibc-matching (failed) ===
```

## 13. M4a: correctly rounded sin, cos, atan2 and acos (2026-11-09)

ADR 0006 item 1's trigonometry, built and verified as M3c-2's primitives were: `_km_cr_sin`, `_km_cr_cos`, `_km_cr_sincos`,
`_km_cr_atan2` and `_km_cr_acos` in `src/mathx.cyr`, each two-phase (Ziv: a fast path with a proven error bound and a
rounding test, a triple-double accurate path), + − × ÷, √ and integer bit operations only, on tables `mathx_tables.py`
generates from MPFR. No module called them in M4a; M4d moved gas.cyr and the `math` corpus onto them (§14).

- **sin, cos, sincos** — the whole double range. |x| < 2²⁰: Cody–Waite reduction by π/2 = P1 + … + P5 (P1, P2 at 33 bits, so
  k·P1 and k·P2 are exact); |x| ≥ 2²⁰: Payne–Hanek, the 288 bits of 2/π from x's exponent on (1,344 bits generated from
  MPFR and checked against Machin's formula in integers) multiplied by the significand in exact 32-bit columns, the
  fraction centred and summed into a triple-double. The reduced r is split at i/256 and sin/cos(i/256) read from a
  four-word table; the fast path's short polynomials leave |hi + lo − f| < 2⁻⁶⁸·|U_A| + 2⁻⁶⁹·|hi| + the reduction's error
  (measured at most 0.17 of the bound over 750,000 arguments against 400-bit MPFR); the accurate path is a triple-double
  Taylor series and the table's four words (below 2⁻¹⁴⁸ relative by the operation count, 2⁻¹³⁸ at worst through the
  reduction; measured 2⁻¹⁵⁶·⁷). The reduction rests on `mathx_tables.py`'s continued-fraction bound: no double is within
  4.687 × 10⁻¹⁹ of a nonzero multiple of π/2, attained by Lefèvre and Muller's 6381956970095103·2⁷⁹⁷. Lefèvre and Muller's
  worst cases (sin on [2⁻²⁴, 2 + 4675/8192], cos on [2⁻²⁵, 12867/8192]) and Hanrot et al.'s large-argument ones are well
  inside the accurate path's bound.
- **atan2** — never rounds y/x: with c = j/256 near n/d, atan(n/d) = atan(c) + atan(δ), δ = (n − c·d)/(d + c·n) from exact
  partial products (c has 8 bits, the inputs are split on a 2⁻⁴⁴ grid); atan(j/256) is a triple-double table, and the octant
  offsets 0, π, π/2 triple-doubles. Fast path: δ by Dekker's division and a degree-9 odd polynomial (2⁻⁶⁹·⁴⁴ relative;
  measured 2⁻⁷⁰·⁸); accurate path: a triple-double long division and the series to δ¹⁹ (2⁻¹⁴⁹; measured 2⁻¹⁵⁶·⁹). A quotient
  below 2⁻⁵⁹ is rounded directly (π/2, π, or the quotient itself, its exact subnormal midpoints broken toward zero), and C99
  Annex F's special values are glibc's and MPFR's on all 361 special pairs.
- **acos** — atan2(√(1 − x²), x): 1 − x² as a double-double ((1 − |x|)·(1 + |x|) above ½, exact by Sterbenz) on the fast path,
  an exact four-term sum on the accurate one, √ by one and two Newton steps; acos(±1), acos(±0) and the tiny |x| exact.
- **`_km_td_norm`** — the accurate atan and acos paths divide, and take the square root of, triple-doubles formed after a
  cancellation (n − c·d; the residual w − s²). `_km_td_add`'s last renormalisation does not see its rounded third-order term,
  so after a cancellation its leading word can be far from the sum, while each correction of the long division and the
  Newton step is q = r₀/b₀ from the leading word and assumes it is the sum to 2⁻⁵². `_km_td_norm` (two bottom-up TwoSum
  passes, exact) renormalises before each division and correction: it is what makes `_km_td_div`'s 2⁻¹⁵¹ and
  `_km_td_sqrt`'s 2⁻¹⁵⁰ hold, and with it the accurate paths agree with MPFR on every argument below.

**Verification** (`scripts/oracle/cr_verify.py`, every result bit for bit against MPFR 4.2.2; the same run in M4d):

| Function | Arguments | Accurate path forced | Mismatches (ELF · PE · aarch64) | Fast path falls back |
|---|---|---|---|---|
| sin | 1,129,805 | 100,225 | 0 · 0 · 0 | 79 (7.0 × 10⁻⁵) |
| cos | 1,129,805 | 100,225 | 0 · 0 · 0 | 98 (8.7 × 10⁻⁵) |
| sincos (both, against sin and cos) | 1,129,805 | 100,000 | 0 · 0 · 0 | 177 (1.6 × 10⁻⁴) |
| atan2 | 1,014,417 pairs | 100,411 | 0 · 0 · 0 | 19 (1.9 × 10⁻⁵) |
| acos | 1,026,077 | 100,065 | 0 · 0 · 0 | 14 (1.4 × 10⁻⁵) |

The arguments: every binade up to 2¹⁰²³, the Peng–Robinson range [−2.2, 3.2], RN(kπ/2) ± ulps for k < 2⁴⁰ and in every
binade the doubles nearest a multiple of π/2 (from the continued fraction of 2^(e − 52)·2/π), the FFT angles −2π/n, the
Cody–Waite / Payne–Hanek edge, the table edges, tiny |x|; atan2 over every quadrant and exponent difference, results near
±π and ±π/2, subnormal and midpoint quotients, |y/x| near the table points and 1, bond and dihedral angles at Å and metre
scale; acos dense near ±1 and 0, over every binade and at the branch points ½ and 1/√2; the published hard cases. The three
targets give the same bits and fast-path flags on every argument. `tests/mathx.tcyr` holds the special values, hard cases,
Payne–Hanek arguments, sincos, atan2's sign and zero matrix and acos values (MPFR's bits), and seeded fast / accurate
agreement.

**Cost** (`tests/kimiya.bcyr`, M4a): `cr_sin` 74 ns (`f64_sin` 37), `cr_cos` 82 (`f64_cos` 40), `cr_sincos` 115, `cr_atan2`
113 (`ganita_f64_atan2` 53), `cr_acos` 140 — against ADR 0005's double-double `_km_acos` 5.8 µs and `_km_cos` 1.68 µs, which
gave the same bits as the new functions on all 1,026,077 acos and 494,663 cos (|x| ≤ 8) arguments.

## 14. M4d: potential, the trig switch, ADR 0002 #13 completed, #16–#18 (2026-11-09)

### The trig switch

gas.cyr's Peng–Robinson cubic calls `_km_cr_acos` and `_km_cr_cos` (ADR 0005's `_km_acos` and `_km_cos` before): **every gas
row keeps its bits** — `peng_robinson_z` 4,603 values, `peng_robinson_z_roots` 5,539 roots, `fugacity_coefficient_pr` and
`fugacity_pr` 4,544 values each (6 and 5 of them ADR 6 rows), all 0 ulp on every target, as before. `_km_cos`, `_km_acos`
and their helpers (`_km_dd_cos_sin`, `_km_dd_asin`, `_km_dd_add`, `_km_dd_mul`, `_km_dd_div_d`, `_km_dd_sqrt`) are removed;
`src/mathx.cyr` no longer includes `lib/math.cyr`.

The `math` corpus's trigonometry is gated on the correctly rounded primitives at bound 0:

| Row | Port primitive | Rows | ADR 6 | ulp (every target) | Shadow (never gated), against the corrected oracle |
|---|---|---|---|---|---|
| `f64::sin` | `_km_cr_sin` | 710 | 2 (both = corrected) | 0 | `f64_sin`: 0: 690 · 1: 20 |
| `f64::cos` (\|x\| > 8 included) | `_km_cr_cos` | 2,910 | 4 | 0 | `f64_cos`: 0: 2,828 · 1: 82 |
| `f64::acos` | `_km_cr_acos` | 2,620 | 1 | 0 | `ganita_f64_acos`: 0: 1,892 · 1: 727 · 2: 1 (aarch64 1,875 · 743 · 2) |
| `f64::atan2` | `_km_cr_atan2` | 1,000 | 0 | 0 | `ganita_f64_atan2`: 0: 797 · 1: 203 (aarch64 798 · 202) |

Until M4d `f64_sin` (20 rows at 1 ulp), `ganita_f64_atan2` (203) and `f64_cos` past \|x\| = 8 (6) were gated at the stdlib's
2 ulp. **No golden transcendental row is on the stdlib any more.** What still is: `f64::cbrt` (`ganita_f64_cbrt`, class A:
correctly rounded, as the oracle's is on every golden argument, §9), and the exact operations `sqrt` (IEEE), `abs`,
`floor`, `trunc` and `clamp` (`f64_clamp`), which need no port version.

### potential

One new golden module, 26,938 rows, every row 0 ulp from its reference on x86_64 ELF, x86_64 PE and aarch64:

| Function | Class | Rows | Grid | ADR | ulp |
|---|---|---|---|---|---|
| `potential::lennard_jones` | B | 2,105 | each argument over the edges (Ar's values, unit values), the Rust tests' values, r/σ 0.5–100 at Å and metre scale, σ/r at (σ/r)⁶'s and (σ/r)¹²'s over- and underflow, 600 seeded; ε ≥ 2¹⁰²² (6 edges × 10 ratios × 2 scales, 300 seeded) | **#17**: 180 rows | 0 |
| `potential::lj_equilibrium_distance` | C (Rust's powf) | 6,961 | the edges, every decade × 1 and 3.4, every power of two 2⁻¹⁰⁷⁴ … 2¹⁰²³, the 81 doubles around DBL_MAX/2^(1/6), 3,000 seeded over every binade (both signs), 500 Å-scale σ | **#16**: 1,041 rows | 0 |
| `potential::morse` | C (exp) | 2,580 | each argument over the edges, the Rust tests' values, HCl at metre scale over 0.5–5 r_e, a(r − r_e) on targets across exp's range, **1,600 near-cancellation rows** (\|a(r − r_e)\| from 10⁻¹⁷ to 1, both signs), 600 seeded | — | 0 |
| `potential::morse_width_parameter` | B | 4,509 | the edges, the HCl values, every decade of each argument, 400 seeded; D_e ≥ 2¹⁰²³ (4 edges × 9 μ × 4 ω, 300 seeded) | **#18**: 181 rows | 0 |
| `potential::coulomb`, `coulomb_force` | B | 1,930 each | the edges, the Rust tests' charges, every decade of r, the intermediates' over- and underflow, 600 seeded | — | 0 |
| `potential::distance_3d` | B | 846 | the Rust tests' points, each coordinate over the edges, 12 scales from 2⁻¹⁰⁵⁰ to 10³⁰⁰ (squares over- and underflowing), 700 seeded | — | 0 |
| `potential::bond_angle` | C (atan2) | 3,156 | the Rust tests' geometries, tiny angles near 0 and π, each coordinate over the edges, 14 scales, −0.0 coordinates, **80 scans across the 10⁻¹² coincidence threshold to the ulp** (A − B = (0, 0, δ) exactly), 1,000 seeded | ADR 6: 1 | 0 |
| `potential::dihedral_angle` | C (atan2) | 2,921 | the Rust tests' geometries, trans / cis / gauche ±60° / ±90° / ±120° with ±0.0 z-terms at four scales, every sign of the 7 zero coordinates of the Rust test's trans and cis geometries (384 rows; 12 where atan2 gives −π and Rust reports π, 16 of −0.0), each coordinate over the edges, 10 scales from 2⁻¹⁰⁵⁰ to 10¹⁵⁵, consecutive atoms coinciding, **80 scans across the 10⁻²⁰ collinearity threshold to the ulp**, 1,000 seeded | — | 0 |

By class: B 11,320, C 15,618 rows.

- **The corrected oracle.** glibc's atan2 misrounds 1 of the module's 4,439 atan2 calls (2,836 distinct arguments): one
  bond-angle row changes and goes to `tests/golden/potential.cr.tsv`; the port is bit-identical to the corrected value.
  glibc's exp misrounds one Morse call without changing its row, and its exp2 at RN(1/6) is correctly rounded
  (RN(2^(1/6)) = 0x3FF1F59AC3C7D6C0, which is why #16's misses are the second rounding).
- **The hard rows** (ADR 0006 item 6). The Morse near-cancellation rows put `morse`'s S at 8.4 × 10¹⁸ on 1,844 rows
  (`--sensitivity`; S = 30 elsewhere): 1 − e^y cancels, so a 1-ulp exp result changes x = 1 − e^y in its leading bits; only
  a correctly rounded exp meets them. `dihedral_angle` has S = 9.2 × 10¹⁸ on the 12 −π → π rows (a nudged −π is no longer
  ≤ −π, and the result flips sign) and 1 elsewhere; `bond_angle` 1, `lj_equilibrium_distance` 2. All are gated at 0.
- **The tolerances** — Rust's private 10⁻¹² (squared at run time) and 10⁻²⁰, typed inline: `validate.py` changes 544 / 454
  bond-angle rows with 10⁻¹² ± 1 ulp (926 / 1,004 with × ½ / × 2) and 108 / 140 dihedral rows with 10⁻²⁰ ± 1 ulp (528 / 642),
  so the corpus pins both to the ulp; the π of the −π mapping (12 rows), Lennard-Jones' 4.0 and the width's 2.0 too.
- **ADR 0002 #16, #17, #18** (ADR 0002's M4d addition). #16: the port returns σ·2^(1/6) rounded once
  (`_km_cr_sixth_root2_mul`); Rust's σ·RN(2^(1/6)) is correctly rounded on 5,914 of the 6,955 finite nonzero σ rows and one
  ulp off on 1,041 (and on 4,117 of 20,000 seeded σ); `validate.py` proves every value exactly (x·2^(1/6) > t ⟺ 2x⁶ > t⁶,
  each candidate's midpoints compared in `Fraction`), and `mathx_tables.py`'s continued fraction of 2^(1/6) (min over
  q ≤ 2⁵⁴ of \|q·2^(1/6) − p\| = 2⁻⁵⁴·¹⁴) proves that the triple-double product, within 2⁻¹⁰²·⁶ ulp, needs no rounding test.
  #17: 180 of the 382 Ok rows with ε ≥ 2¹⁰²² differ ((ε·d)·4 = RN(4·ε·d); Rust ±∞ or NaN); Rust's (4·ε)·d is proven
  RN(4·ε·d) on every other finite row. #18: 181 of 409 Ok rows with D_e ≥ 2¹⁰²³ differ ((μ·½)/D_e = RN(μ/(2·D_e)); Rust's
  a = 0), proven on them and on 20,000 seeded pairs.
- **Reviewed and kept as Rust's** (no single exact value; reported to the user): (σ/r)⁶ and (σ/r)¹²'s over- and
  underflow, k·q₁·q₂ and r·r's, the squared coordinates' in the geometry, `morse`'s r − r_e overflow and subnormal D_e·x.
  The grids include them, so a later decision can be measured. (The geometry's: fixed in M4e-2 as ADR 0002 #20, §16.)

### ADR 0002 #13 completed

`atoms_after_half_lives` past the normal range is `_km_cr_exp2_scaled(N₀, −n)`, N₀·2⁻ⁿ rounded once for every N₀ and n:
with N₀ = s·2^j and −n = k + f (k = trunc(−n)), the exact s·2^(j + k) when f = 0, and otherwise 2^f from exp's fast core
(f·ln 2 as a double-double, within 2⁻⁹⁶·⁹), s·2^f as an exact TwoProduct plus a rounded tail (within 1.8·2⁻⁷⁰), exp's
rounding test at the final scale (normal, or the subnormal grid) with err = 2⁻⁶⁸, and the triple-double accurate path when
it fails (exp's accurate core, `_km_td_mul_d`, `_km_cr_finish`: below 2⁻¹⁵⁰). 2^f is irrational for a non-zero double f in
(−1, 1), so there are no ties to break.

| Check | Arguments | Result |
|---|---|---|
| `cr_verify.py`'s `exp2s` against `mpfr_ref.scaled_exp2` (2^f to 400 bits, bounded from both sides) | 611,518 seeded pairs (nuclear's edge 151,290, subnormal results 120,000, the overflow end 50,000, exact scalings and ties 60,000, small and near-2200 \|y\| 80,000, every binade × y in [−2300, 2300] 150,000, 228 special pairs), 100,228 through the accurate path alone | 0 mismatches on x86_64 ELF, x86_64 PE and aarch64, the same bits and flags; the fast path falls back on 10 (1.6 × 10⁻⁵) |
| `validate.py` | every one of the 2,505 golden rows where 2⁻ⁿ is not a normal double, and 4,000 seeded pairs | each CR(N₀/2ⁿ) by MPFR and by a 60-digit 2⁻ᶠ with exact rational rounding |
| the generator (`exact::cr_exp2_frac`) | every golden row with a fractional n | 2⁻ᶠ to 320 bits in integers, two-sided; agrees with the port on every row |
| `tests/mathx.tcyr` | 57 MPFR values (specials, ties, subnormal and overflow edges, eight found fast-test failures) and 4,000 seeded fast / accurate pairs | all equal |

1,810 golden rows now carry the port's value in their note (564 integer n, 173 N₀ = ±2^j, 1,067 fractional n, 6 scaleB);
1,812 did with M4c's faithful values, which were one ulp off on 378 of the 1,121 fractional-n rows. Rust's own value is
correctly rounded on 686 of the 2,505 edge rows.

### The corpus and the summary

Golden rows: `potential` 26,938 (new) and its 1 corrected row; `nuclear`'s notes rewritten (1,810 #13 rows); every other
file byte-identical, and a second `--write-corrections` leaves the eight `.cr.tsv` files unchanged.

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **16 modules, 0 failing; 244,803 ok, 0 failed, 0 measured,
0 pending, 4,644 accepted (689 under ADR 6; under ADR 2, 440 #4 Antoine, 13 #12 CFSE, 23 #3 SEMF, 1,810 #13, 208 #14,
59 #15, 1,041 #16, 180 #17 and 181 #18 rows), 0 glibc-matching.**

```text
=== parity summary: 16 modules, 0 failing | rows: 244803 ok, 0 failed, 0 measured, 0 pending, 4644 accepted (ADR 6: 689), 0 glibc-matching (failed) ===
```

**A note for CI.** Two `cyrius test --win` runs at once under wine have printed a truncated parity run and still exited 0.
The harness can't tell (its exit status and output are what wine hands back), so a gate must check that the summary line
above is present; M4d's runs were made one target at a time.

## 15. M4e: spectroscopy, the FFT transcription, ADR 0002 #19 (2026-11-09)

One new golden module, 52,985 rows, every row 0 ulp from its reference on x86_64 ELF, x86_64 PE and aarch64: the oracle's
value, the corrected oracle's on 118 rows (ADR 6), or the exact value the note names on 1,654 rows (ADR 0002 #19). The
three class-C functions are on correctly rounded primitives (exp10, log10, and sin / cos for the FFT's twiddles), in
Rust's and hisab 1.4.0's operation order, at **bound 0**; the class-B ones use + − × ÷ and √ only (the Rydberg formula
and the isotope pattern included: no transcendental, §5's examples notwithstanding).

| Function | Class | Rows | Grid | ADR | ulp (every target) |
|---|---|---|---|---|---|
| `spectroscopy::absorbance` | B | 560 | each argument over the edges, the Rust tests' ε = 100, l = 1, c = 0.01, the products' over- and underflow, seeded | — | 0 |
| `spectroscopy::transmittance` | C (pow(10, ·)) | 2,587 | the edges, every integer −330 … 330 (10⁰ … 10²² exact, 10²³ a tie), A over exp10's range in steps of 0.8, **subnormal results** (A in 307.6 … 323.6), the Rust tests' values, seeded | **ADR 6: 3** (glibc's pow) | 0 |
| `spectroscopy::absorbance_from_transmittance` | C (log10) | 2,515 | the edges, every decade × 4, 1 and the 40 ulps below it (A ≈ 0), the round trips T(A) for A = −10 … 165, (0, 1] seeded uniformly and log-uniformly | **ADR 6: 115** (glibc's log10) | 0 |
| `spectroscopy::concentration_from_absorbance` | B | 564 | each argument over the edges (the guards' order), the Rust tests' values, ε·l over- and underflowing, seeded | — | 0 |
| `spectroscopy::photon_energy_from_frequency`, `photon_energy_from_wavelength`, `wavelength_from_energy`, `frequency_from_wavelength`, `joules_to_ev`, `ev_to_joules` | B | 1,890 each | the edges, every decade × 1 and 5 from 10⁻³²³ to 10³⁰⁸, the visible 380–780 nm, the Rust tests' values, seeded | — | 0 |
| `spectroscopy::bohr_energy_level` | B | 5,120 | **every u8 Z × 20 u32 n** (0 … 12, 100 … 10⁶, 2³¹, u32::MAX − 2 … u32::MAX) | — | 0 |
| `spectroscopy::bohr_transition_energy`, `rydberg_wavelength` | B | 1,816 each | 8 Z (0, 1, 2, 3, 26, 92, 118, 255) × 13 × 13 (n₁, n₂) with every Err path in Rust's order, the Rust tests' values, **close large n** (n₂ = n₁ + 1 … 7 for n₁ from 10 to 2³¹), seeded | — | 0 |
| `SpectralSeries::n_lower`, `::wavelength` | A, B | 5, 195 | every series × lines 0 … 24, 50 … 2³¹ and u32::MAX − 8 … u32::MAX (**the `checked_add` overflow**; Lyman at u32::MAX − 1 is 1/R_H, which pins R_H) | — | 0 |
| `spectroscopy::wavenumber_from_wavelength`, `wavelength_from_wavenumber` | B | 8,119 each | the edges, every decade × 1, 2.5, 5, **every power of two 2⁻¹⁰⁷⁴ … 2¹⁰²³**, the subnormal-result end (x past 2¹⁰¹⁵) and the overflow end (subnormal x), 3,000 seeded over every binade, IR and radio wavelengths | **#19: 827 each** | 0 (from the oracle or the note) |
| `spectroscopy::reduced_mass` | B | 966 | each argument over the edges, both and **equal masses** (300 seeded m), the Rust tests' values, H³⁵Cl in u and kg, the product's and sum's over- and underflow, seeded | — | 0 |
| `spectroscopy::harmonic_oscillator_energy`, `rigid_rotor_energy` | B | 1,015 each | the 20 u32 probes × the Rust tests' ω / B, the edges, every decade, seeded | — | 0 |
| `spectroscopy::vibrational_frequency` | B | 3,047 | each argument over the edges and every decade, the Rust tests' (500, 1.66 × 10⁻²⁷) and potential's HCl (516, 1.6267 × 10⁻²⁷), k/μ over- and underflowing, seeded | — | 0 |
| `spectroscopy::rotational_constant`, `diatomic_moment_of_inertia` | B | 1,587, 440 | the edges, every decade, the Rust tests' values, HCl and CO, seeded | — | 0 |
| `spectroscopy::isotope_abundances` | A | 256 | **every u8** | — | — (32 floats, 0) |
| `spectroscopy::isotope_pattern` | B | 696 | the Rust tests' molecules (CH₄, Cl₂, H₈₈, C₁₆₈, C₁₀₀₀H₂₀₀₀O₅₀₀S₅₀Br₂₀, glucose, the unknown and data-less elements), the empty molecule, **every u8 Z alone**, every element with data × 15 counts to 2,000 and singles to 5,000 (the pruning keeps the peak count near 7.4·σ), elements without data to u32::MAX, count 0, repeated entries, H₂ … (commutative ties), the Err position in a list, 300 seeded molecules | — | 0 (12,526 floats) |
| `spectroscopy::fft_spectrum` | C (cos, sin) | 213 | **every length 0 … 16**, the Rust test's [1; 10], 2^k ± 1 for k = 5 … 9 and 2¹² + 1 (Err); **every power of two 2⁰ … 2¹³** with impulses, constants, ±0, alternating (Nyquist), sinusoids on a bin and off, random, huge (the sums overflow), subnormal, NaN and ±∞ samples (20 kinds to 2⁶, 9 to 2¹⁰, 3 to 2¹², 2 at 2¹³) | — | 0 (27,512 floats) |
| `spectroscopy::find_peaks` | A | 994 | the Rust tests' data × 16 thresholds (0, 1, −0.0, the ulps inside [0, 1] and every way out, NaN, ±∞), plateaus, flat and monotone data, NaN and ±∞ in the data, all non-finite, ±0, DBL_MAX ranges, **thresholds landing on a data value to the ulp**, 500 seeded sequences | — | — (integers) |

By class: A 1,255, B 46,415, C 5,315 rows.

- **The corrected oracle.** glibc's pow(10, ·) misrounds 3 transmittance rows (10^23, the exact midpoint, among them) and
  its log10 115 absorbance rows; they go to `tests/golden/spectroscopy.cr.tsv`, and the port is bit-identical to the
  corrected outputs on all 118 (0 match glibc's). glibc's cos and sin are correctly rounded at the 13 FFT angles
  −2π/2^k (and to k = 20), so no FFT row changes. The release oracle calls glibc's `sincos` for hisab's twiddles (LLVM
  joins `angle.cos()` and `angle.sin()`; the binary imports it once `fft_spectrum` is called): `validate.py`'s replica,
  with Python's `math.cos` and `math.sin`, gives every golden FFT row's bits, so `sincos` is glibc's sin and cos there.
- **The FFT** (ADR 0003 amendment M4e). `_km_fft` is hisab 1.4.0's `num::fft` transcribed — the bit reversal, the
  butterflies with Complex's mul / add / sub, the `w *= wn` recurrence, no FMA — with `_km_cr_sincos` twiddles; hisab
  3.4.1's `num_fft` (stdlib twiddles, reseeded every 64 steps) would differ from n = 256 on. `fft_spectrum` replays the
  rows through the heap work table up to n = 64 and in place above; `tests/spectroscopy.tcyr` checks the two paths
  bit-identical to 2¹⁴, impulses flat to 2¹⁶ and Parseval to 10⁻¹⁰. `--sensitivity`: S = 1 on 129 of the 168 rows with a
  nudged twiddle, ∞ on 39 (a nudged twiddle turns an infinite or overflowing sum into another: the huge, ±∞ and DBL_MAX
  signals); gated at 0 like every row.
- **ADR 0002 #19.** Rust's `0.01 / x` is one ulp from RN(1/(100·x)) on 827 of the 8,107 Ok rows of each function (and
  on 2,748 of 20,000 seeded x over every binade, 14 %); `validate.py` proves every note with `Fraction`, and Rust's value
  is within one ulp everywhere the result is normal. The rows differ by one ulp, inside class B, so the harness marks a
  noted row as mapped (`_par_force_map`) and counts it as accepted under ADR 2, never as a plain ok. `_km_cr_hundredth_div`
  is verified by `cr_verify.py` (`h100`, below), `tests/mathx.tcyr` (28 MPFR values, 4,000 seeded fast / exact pairs) and
  `tests/spectroscopy.tcyr` (20,000 seeded).
- **The typed constants** (`validate.py`'s pinning): the isotope thresholds × ½ / × 2 change 110 / 117 (10⁻¹²), 153 / 173
  (0.001) and 173 (0.5 × 2) of the 693 pattern rows with at most 3,000 steps, ± 1 ulp none (no peak lands on a
  threshold, no mass difference is near 0.5), so `tests/spectroscopy.tcyr` pins each to the ulp on its helper; ħ's 2.0
  + 1 ulp changes 942 oscillator rows, the 0.5 12, B's 8.0 1,331 rotational-constant rows, R_H ± 1 ulp 157 / 166 series
  rows, the FFT's −2.0 − 1 ulp 41 FFT rows; the rotor's (J + 1.0) is unobservable on every input (J + 1 + 2⁻⁵² rounds back
  for J ≥ 1, and J = 0 multiplies by 0).
- **The stable sort.** No natural molecule's scratch holds equal masses with unequal intensities (the ties from
  m₁ + m₂ = m₂ + m₁ are identical records), so reversing the ties changes no golden row; `tests/spectroscopy.tcyr` pins the
  order on `_km_isotope_sort` with −NaN, ±0, +NaN and equal masses of different intensities, and the binned mean's order.
- **Reviewed and kept as Rust's** (ADR 0002's M4e addition; reported to the user): `reduced_mass(m, m)` misses m/2 on 33 of
  the 315 equal-mass rows; `bohr_energy_level(Z, Z)` misses −R on 3 of 12; computing the Rydberg factor 1/n₁² − 1/n₂²
  exactly would move `rydberg_wavelength` by up to 38,777,186 ulp ((255, 2³¹, 2³¹ + 7)); no isotope peak sits exactly on
  RN(0.001). The grids include them, so a later decision can be measured. (The user decided the first three in M4e-3:
  ADR 0002 #21 and #22, §17.)

**Verification of `_km_cr_hundredth_div`** (`scripts/oracle/cr_verify.py`'s `h100`, against `mpfr_ref.hundredth_div`):

| Arguments | Integer path forced | Mismatches (ELF · PE · aarch64) | Fast path declines |
|---|---|---|---|
| 1,049,032 (specials; the fast path's edges 2^±900, the overflow end and the subnormal-result end, ±3,000 ulps each; the integers 1 … 150,000; 150,000 decimal wavelengths k·10⁻ʲ; 700,000 binade-uniform of both signs) | 100,018 | 0 · 0 · 0 (the same bits and flags) | 135,071 (12.9 %: outside [2⁻⁹⁰⁰, 2⁹⁰⁰], by design; inside, none of the 300,000 integers and decimals) |

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **17 modules, 0 failing; 296,016 ok, 0 failed, 0 measured,
0 pending, 6,416 accepted (807 under ADR 6; under ADR 2, 440 #4 Antoine, 13 #12 CFSE, 23 #3 SEMF, 1,810 #13, 208 #14,
59 #15, 1,041 #16, 180 #17, 181 #18 and 1,654 #19 rows), 0 glibc-matching.**

```text
=== parity summary: 17 modules, 0 failing | rows: 296016 ok, 0 failed, 0 measured, 0 pending, 6416 accepted (ADR 6: 807), 0 glibc-matching (failed) ===
```

The harness's token lists hold 8,448 tokens since M4e (`_PAR_MAX`, for the 2¹³-sample FFT rows: 8,193 input tokens).
Every target ran one at a time (§14's note for CI).

## 16. M4e-2: potential's geometry with an unbounded exponent, ADR 0002 #20 (2026-11-10)

The user's decision of 2026-11-09 ("Yes, prescale") on M4d's report: `distance_3d`, `bond_angle` and `dihedral_angle`
are Rust's own operation sequence evaluated in binary64 with an unbounded exponent range (every + − × ÷ √ rounded to 53
bits, nothing over- or underflowing, comparisons and Err checks on those values), the final √ or atan2 rounded once into
binary64 (ADR 0002's M4e-2 addition has the target, the fast-path tests and the helpers' exactness argument). Non-finite
coordinates keep Rust's results.

| Function | Class | Rows (M4d → M4e-2) | Finite coordinates | Rust's chain in the normal range | #20 noted (Err) | ulp (every target) |
|---|---|---|---|---|---|---|
| `potential::distance_3d` | B | 846 → 1,914 | 1,896 | 1,000 | 854 (0) | 0 |
| `potential::bond_angle` | C (atan2) | 3,156 → 4,838 | 4,809 | 3,172 | 1,205 (0) | 0 |
| `potential::dihedral_angle` | C (atan2) | 2,921 → 4,994 | 4,957 | 2,892 | 1,546 (91) | 0 |

- **Only edge rows change.** 314 of M4d's 26,938 rows take a note (the scales 2^±511 and beyond, subnormal and DBL_MAX
  coordinates); `validate.py` proves that no row whose f64 chain stays in the normal range (no product or sum non-finite
  or below 2⁻¹⁰²², no zero product of nonzero factors) is noted, and that every unnoted row's expected result — the
  corrected oracle's where glibc's atan2 misrounded, else the oracle's — is MPFR's unbounded evaluation.
- **The new edge grids** (4,823 rows, `m_potential.rs`'s `edges20`): each geometry family at 2^k for k in ±(230 … 1074),
  one distance at every seventh binary scale, differences that overflow (DBL_MAX against −DBL_MAX), mixed magnitudes and
  subnormal coordinates (seeded, binade-uniform), tiny angles (1, 2^t) for t down to −1500 at six scales (subnormal results,
  their midpoints, and π − t near π), the coincidence and collinearity tolerances at extreme scales with the port's and
  Rust's crossings found by bisection (Rust's squares rounded at the subnormal grid near 2⁻⁴⁸⁰ and 2⁻²⁴⁵, overflowing near
  2⁵¹² and 2²⁷⁶), b₁ tiny (2⁻⁵⁰⁰) against b₂ huge (2³⁰⁰), where Rust's threshold rounds 10⁻²⁰·|b₁|² to a few subnormal
  bits before |b₂|² scales it up (the 91 rows where the port says collinear and Rust does not), trans and cis dihedrals
  with every sign of their zeros at 2^±400 … 2^±1000, and 300 seeded geometries per function at random scales.
- **Two new ADR 6 rows** (`potential.cr.tsv`, 1 → 3): tiny bond angles whose quotient is the subnormal midpoint
  1.5·2⁻¹⁰⁷⁴; glibc's atan2 rounds it to even (2·2⁻¹⁰⁷⁴), the correctly rounded atan2 lies just below it (1·2⁻¹⁰⁷⁴). In
  both rows Rust's chain is in range, so the port's value is the corrected oracle's (no #20 note).
- **The generator** computes the port's result in exact integer arithmetic (`exact::Ux`; `exact::ux_atan2_f64`, a 400-bit
  fixed-point atan with Ziv's test) and notes a row where the unbounded chain differs from its f64 replica of Rust's
  (an Err decision or atan2's arguments) and the result from the oracle's; the replica is checked against the oracle on
  every row.
- **Scale invariance**: `validate.py` evaluates 6,853 rows again at an exact random scale 2^k through MPFR (the angles the
  same bits, distances 2^k times); `tests/potential.tcyr` checks the port bit for bit on 600 random geometries × 75 scales
  from 2⁻¹¹⁰⁰ to 2¹¹⁰⁰ (unit-scale, planar, near-coincident, near-collinear and mixed-scale ones), and the forced slow path
  against the fast one on 5,000 in-range geometries.

**Verification of the helpers and the slow path** (`scripts/oracle/cr_verify.py`, against `mpfr_ref.UX`: 53 bits at MPFR's
widest exponent range, the final operation by `mpfr_check_range` and `mpfr_subnormalize`):

| Family | Arguments | Mismatches (ELF · PE · aarch64) |
|---|---|---|
| `uxsqrt` (`_km_ux_sqrt_f64`) | 200,009 (subnormal results and exact-square ties, the overflow end, every exponent to ±5000) | 0 · 0 · 0 |
| `uxmul`, `uxdiv` | 180,006, 180,005 (subnormal results and short-significand ties, the overflow end, results below 2⁻¹¹³⁰) | 0 · 0 · 0 |
| `uxadd` | 200,011 (near the subnormal range, sticky exponent gaps 56 … 200, cancellations, midpoints and a breaking term) | 0 · 0 · 0 |
| `uxatan2` (accurate path forced on 100,008) | 280,008 (exponent differences to ±2300 at common exponents to ±5000, subnormal and exact quotients) | 0 · 0 · 0 |
| `dist`, `bond`, `dih` (public function and forced slow path) | 38,000 each (random scales 2⁻¹⁰⁷⁴ … 2¹⁰¹⁸, mixed magnitudes, subnormal coordinates, near coincidence and collinearity at scale, overflowing differences) | 0 · 0 · 0 |

`tests/mathx.tcyr` adds 48 MPFR values past binary64's range and 20,000 seeded pairs against IEEE's own operations (at the
doubles' exponents and shifted past the range and back), and atan2 against `_km_cr_atan2` at any common shift.

**Cost** (`tests/kimiya.bcyr`, same host, before → after): `potential/dihedral_angle` 131 → 133 ns (Rust 20.35), `potential/
bond_angle` 213 → 213, `potential/distance_3d` 14 → 16; the slow path on the same geometries scaled by 2⁶⁰⁰ (Rust's
squares overflow): `distance_3d_2p600` 314 ns, `bond_angle_2p600` 1.13 µs, `dihedral_angle_2p600` 1.31 µs (Rust's chain gave
+∞ or an Err there in 14, 52 and 80 ns).

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **17 modules, 0 failing; 297,232 ok, 0 failed, 0 measured,
0 pending, 10,023 accepted (809 under ADR 6; under ADR 2, 440 #4 Antoine, 13 #12 CFSE, 23 #3 SEMF, 1,810 #13, 208 #14,
59 #15, 1,041 #16, 180 #17, 181 #18, 1,654 #19 and 3,605 #20 rows), 0 glibc-matching.**

```text
=== parity summary: 17 modules, 0 failing | rows: 297232 ok, 0 failed, 0 measured, 0 pending, 10023 accepted (ADR 6: 809), 0 glibc-matching (failed) ===
```

## 17. M4e-3: the exact Rydberg rationals and the exact reduced mass, ADR 0002 #21 and #22 (2026-11-10)

The user's decisions of 2026-11-10 ("Compute exactly", "Round once"; ADR 0002's M4e-3 addition quotes them): the Bohr
levels, the transition energy, the Rydberg formula and the series return their exact rationals with the constant as its
double, rounded once, and `reduced_mass` the correctly rounded m₁·m₂/(m₁ + m₂) for every pair of positive finite masses.
Each has a double-double fast path with a proven bound and a rounding test, and an exact slow path on `src/mathx.cyr`'s
multi-word integers (`_km_cr_ratio`, RN(p/q·2^e), ties to even, subnormal results on the 2⁻¹⁰⁷⁴ grid).

| Function | Class | Rows (M4e → M4e-3) | Ok rows | #21 / #22 noted | Exact ties | Exact values | Rust's worst (ulp) | ulp (every target) |
|---|---|---|---|---|---|---|---|---|
| `spectroscopy::bohr_energy_level` | B | 5,120 → 11,899 | 11,624 | 3,823 | 0 (none exist) | 6,869 (every n = 2^k·Z) | 2 | 0 |
| `spectroscopy::bohr_transition_energy` | B | 1,816 → 2,959 | 2,065 | 1,351 | 13 (all that exist) | 665 | 1,814,015,630 | 0 |
| `spectroscopy::rydberg_wavelength` | B | 1,816 → 2,959 | 2,065 | 1,376 | 0 (none exist) | 0 (none exist) | 90,480,107 | 0 |
| `SpectralSeries::wavelength` | B | 195 → 420 | 400 | 67 | 0 (none exist) | 0 (none exist) | 2 | 0 |
| `spectroscopy::reduced_mass` | B | 966 → 5,382 | 5,350 | 2,837 | 195 | 2,575 | 7.0·10¹⁰ (an underflowed product) | 0 |

- **Where the exact values and ties are.** R's odd significands decide (`validate.py` factors them): `KM_RYDBERG_EV`'s is
  683·11214237642539, squarefree with a prime above 2³², so a level is a double exactly when n/Z is a power of two and has
  at least 56 significant bits otherwise; `KM_RYDBERG`'s and R_H's carry a prime above 2³², so no wavelength is dyadic.
  The transition energy's domain is searched exhaustively over its structure (n₁ = g·x, n₂ = g·y, odd(x)·odd(y) | odd(Z),
  odd(g)² | 683·P·j²·w): 5,657 exact values and 13 ties, (255; 2^k·58,055, 2^k·524,544), k = 0 … 12, every one in the
  corpus, and `validate.py` repeats the search with its own factorisation. The reduced mass has tie families: m = m at odd
  multiples of 2⁻¹⁰⁷⁴, (6c, 2c)·2⁻¹⁰⁷⁴, and normal ties such as (2²⁷ − c, 2²⁷ + c) for odd c ≤ 94,906,265, (3·(3·2²⁵ ∓ c))
  and (3·(2²⁶ ∓ c)), at scales across the range and moved onto the subnormal grid.
- **The new rows** (13,706): every exact level; n₁ = n₂ − 1 at 2^k − 1, 2^k, 2^k + 1 for k = 13 … 31 and next to
  u32::MAX; the fast path's n₂ = 2¹³ edge; the 13 ties and each of the 221 exact-value families' smallest, middle and
  largest member; the series at lines 2⁵ … 2³⁰ and u32::MAX − 24 … − 9, and at n₂ = 2¹³; for the reduced mass, m = m at
  every power of two 2⁻¹⁰⁷⁴ … 2¹⁰²³ and at k·2⁻¹⁰⁷⁴, the tie families, exponent gaps 45 … 64 at eight base exponents
  around the proven 55, the 2⁻¹⁰²¹ edge, DBL_MAX pairs and 1e±300 pairs, 600 seeded pairs over every binade (their own
  seed, so every other row keeps its inputs).
- **The proofs.** `validate.py` evaluates every Ok row's exact rational with `Fraction` and checks that a row is noted
  exactly where the rounded value differs from the oracle's and that the note names it bit for bit, so every unnoted row
  is the correctly rounded value too; Rust's own replica still recomputes every row's oracle value.

**Verification** (`scripts/oracle/cr_verify.py`, against `Fraction`; the fast path's declines in brackets):

| Family | Arguments | Exact path forced | Mismatches (ELF · PE · aarch64) |
|---|---|---|---|
| `ratio` (`_km_cr_ratio`) | 320,000 (random sizes to 250 bits and exponents to ±1300; results at the subnormal and overflow ends; exact quotients and midpoints) | (no fast path) | 0 · 0 · 0 |
| `bohr` | 260,429 (every n = 2^k·Z, Z × n ≤ 200, n at 2¹³, 2²⁶ and 2³², random n of every bit length) [0] | 106,869 | 0 · 0 · 0 |
| `trans`, `ryd` | 555,944 each (the ties, the exact values with n₂ ≤ 120 and their 4^k multiples, random pairs at every scale, n₂ = n₁ + 1 at every scale, the 2¹³ edge, next to 2³²) [13: the ties; 0] | 100,013 each | 0 · 0 · 0 |
| `series` | 55,345 (every series × lines 1 … 2,000, the powers of two, the largest lines, 9,000 random) [0] | 110,690 | 0 · 0 · 0 |
| `rmass` | 502,000 (binade-uniform pairs, equal masses, subnormal pairs and ties, the tie families, gaps 45 … 64, the 2⁻¹⁰²¹ edge, the overflow end) [84,878: subnormal results by design, and the ties] | 350,000 | 0 · 0 · 0 |

`tests/spectroscopy.tcyr` checks −R at every Z = n and −R/4^k at n = 2^k·Z, the 13 ties, the exact 16R and 32R, m/2 on 10,000
seeded masses, the reduced mass's ties and range ends, and the fast paths against the exact ones on 6,000 seeded (Z, n₁, n₂)
and 20,000 seeded mass pairs; `tests/mathx.tcyr` checks the multi-word operations, `_km_cr_ratio` on ties, subnormal results
and the overflow threshold, and it against `_km_dd_div_cr` on 10,000 seeded quotients.

**Cost** (`tests/kimiya.bcyr`, before → after, back to back): `spectroscopy/bohr_energy_level` 9 → 49 ns (Rust 1.21),
`bohr_transition_energy` 27 → 73, `rydberg_wavelength` 24 → 73, `series_wavelength_balmer_alpha` 28 → 78,
`reduced_mass_hcl` 14 → 73; the exact path on a tie: `bohr_transition_energy_tie` 2.42 µs, `reduced_mass_tie` 2.20 µs.

Summary, every target (x86_64 ELF, x86_64 PE and aarch64): **17 modules, 0 failing; 301,484 ok, 0 failed, 0 measured,
0 pending, 19,477 accepted (809 under ADR 6; under ADR 2, 440 #4 Antoine, 13 #12 CFSE, 23 #3 SEMF, 1,810 #13, 208 #14,
59 #15, 1,041 #16, 180 #17, 181 #18, 1,654 #19, 3,605 #20, 6,617 #21 and 2,837 #22 rows), 0 glibc-matching.**

```text
=== parity summary: 17 modules, 0 failing | rows: 301484 ok, 0 failed, 0 measured, 0 pending, 19477 accepted (ADR 6: 809), 0 glibc-matching (failed) ===
```

## How to reproduce

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
python3 -I scripts/oracle/validate.py .                  # the corpus, and the ADR 0006 corrections (§9)
python3 -I scripts/oracle/validate.py . --write-corrections   # regenerate tests/golden/<module>.cr.tsv
python3 -I scripts/oracle/validate.py . --sensitivity   # S of the class-C functions: conditioning (§5–§8, §10)
python3 -I scripts/oracle/cr_verify.py                  # the correctly rounded primitives against MPFR, every target
python3 -I scripts/oracle/mathx_tables.py --check       # the generated tables and their proofs (2^(1/6)'s bound included)
cyrius test tests/parity.tcyr                  # x86_64 ELF
cyrius test --win tests/parity.tcyr            # x86_64 PE, needs wine
cyrius test --aarch64 tests/parity.tcyr        # aarch64, needs qemu-aarch64
```
