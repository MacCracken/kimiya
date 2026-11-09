# kimiya port parity — measured baselines

Measurements that set the parity classes (port plan §5). The oracle is Rust 1.2.0 on
**glibc 2.44** (Rust std f64 methods call the host libm); the port is Cyrius 6.7.18 with
ganita 1.2.15.

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

**vant_hoff_k's bound.** One exp: K₂ = K₁·exp(−(ΔH/R)(1/T₂ − 1/T₁)). Measured before its first gated comparison
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
  Recorded for ADR 0006; under the decision above it decides nothing.

Summary, every target: **12 modules, 0 failing; 119,044 ok, 0 failed, 22,205 measured (the six kinetics class-C
functions), 10 pending (the M4 inorganic constants), 441 accepted (ADR 0004's row and the 440 ADR 0002 #4 mapped Antoine
rows).**

## How to reproduce

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
python3 -I scripts/oracle/validate.py . --sensitivity   # S, c and K of the class-C functions (§5–§8)
cyrius test tests/parity.tcyr                  # x86_64 ELF
cyrius test --win tests/parity.tcyr            # x86_64 PE, needs wine
cyrius test --aarch64 tests/parity.tcyr        # aarch64, needs qemu-aarch64
```
