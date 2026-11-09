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
primitive exact instead (transcribe it, as tanmatra did for sin/cos/log).

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

## How to reproduce

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
cyrius test tests/parity.tcyr                  # x86_64 ELF
cyrius test --win tests/parity.tcyr            # x86_64 PE, needs wine
cyrius test --aarch64 tests/parity.tcyr        # aarch64, needs qemu-aarch64
```
