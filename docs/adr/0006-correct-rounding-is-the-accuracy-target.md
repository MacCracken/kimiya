# 0006 — Correct rounding is the accuracy target; parity with the oracle is functional

**Status**: Accepted
**Date**: 2026-11-09

> **Status note (2026-11-09, M4d).** Item 1's list is complete for what the port calls today: M4a added sin, cos, sincos,
> atan2 and acos (two-phase, as below), M4d moved gas.cyr and the `math` corpus onto them (ADR 0005's double-double acos and
> cos removed) and added x·2^y rounded once (ADR 0002 #13) and x·2^(1/6) (ADR 0002 #16). Every class-C golden row,
> the `math` primitives included, is gated at 0 ulp; none uses a stdlib transcendental (port-parity.md §13–§14). M4e added
> 1/(100·x) rounded once (ADR 0002 #19) and put spectroscopy's FFT twiddles on `_km_cr_sincos` (hisab 1.4.0's FFT
> transcribed, ADR 0003's M4e amendment, port-parity.md §15). M4e-3 added exact rationals rounded once in the same
> two-phase style (ADR 0002 #21, #22: a double-double quotient with a rounding test, `_km_dd_div_cr`, and the exact
> multi-word `_km_cr_ratio` = RN(p/q·2^e) behind it; port-parity.md §17).

## Context

The port plan's parity policy (§5) made glibc 2.44, the libm behind the Rust oracle, the reference for every
class-C (transcendental) golden row. A function passed if it stayed within K = 2·S·c + 2 ulp of glibc's result,
with S measured by nudging each primitive result one ulp. Three things followed from that:

- **ADR 0005** made acos, cos, ln and x^2.5 correctly rounded for the ill-conditioned Peng–Robinson and
  Redlich–Kwong functions. Its S still counted only the calls where *glibc* was not correctly rounded, so glibc's
  bits remained the target.
- **ADR 0005's M3a amendment** carried that to its end. glibc's log10 misrounds 4.6 % of the golden arguments, so a
  correctly rounded log10 missed the oracle on the near-cancellation `henderson_hasselbalch` rows by 10¹⁸ ulp. M3a
  therefore built `_km_log10_ln`, which reproduces glibc's own composition and keeps glibc's rounding errors on
  purpose.
- **M3b** found six kinetics functions whose S is 10¹⁴–10¹⁸. A subnormal exp or pow result is scaled back into the
  normal range, and one ulp of a subnormal becomes thousands of ulp. On some of those rows the port is right and the
  oracle is wrong: glibc's pow misrounds exact halfway cases such as (3·2⁻²¹⁵)⁵. The M2 functions `vant_hoff_k`,
  `clausius_clapeyron` and `fugacity_pr` have the same latent problem on inputs outside their grids.

The user's decision (2026-11-09): "WE WANT TO BE BETTER THAT glibc", and "parity was in functionality to bring - we
still want to offer the best accuracy as possible".

## Decision

1. **Every transcendental primitive the port calls is correctly rounded** and target-independent. That means
   round-to-nearest-even, subnormal results and exact ties included. Today that is exp, expm1, ln, log10, pow (x^y;
   10^x as exp10; x^2.5), cos and acos; sin, atan2 and any others follow as M4/M5 need them. sqrt is IEEE already,
   and ganita's cbrt is correctly rounded. The primitives live in `src/mathx.cyr` as private `_km_*` functions, built
   from + − × ÷ √ and integer bit operations only, so they give the same bits on every target. Each follows Ziv's
   strategy:
   - a fast path with a proven error bound and a rounding test;
   - an accurate path (triple-double or better) when the test fails;
   - explicit handling of the exact and halfway cases: exp10 of the integers 0 … 23 (10²³ is an exact midpoint),
     log10 of 10⁰ … 10²², and pow's exact and midpoint powers.
   Their tables are generated, like every other table.
2. **The reference is the "corrected oracle".** A class-C golden row's expected value is the oracle's own sequence
   of operations evaluated with every primitive result correctly rounded:
   - The sequence is Rust 1.2.0's code, transcribed in `scripts/oracle/validate.py`. Evaluated with glibc's
     primitives, that transcription already gives the oracle's bits on every row.
   - The correctly rounded values come from MPFR 4.2 (`scripts/oracle/mpfr_ref.py`, called through ctypes, so no
     Python package is needed) at 53 bits in binary64's exponent range, subnormalized.
   - Where glibc rounded every call of a row correctly, the corrected value is the oracle's and nothing changes.
   - Where it didn't, the corrected outputs are generated into `tests/golden/<module>.cr.tsv`, together with the calls
     that glibc misrounded. The harness then compares the port against the **corrected** value only. It counts such
     rows as **accepted under ADR 6** ("the oracle's libm is not correctly rounded here") and reports the count per
     function.
   - A port value that equals the oracle's misrounded value on such a row is a failure, because it means the port's
     primitive is not correctly rounded.
3. **Class C is compared exactly once a function is on the correctly rounded primitives.** Both sides then evaluate
   the same correctly rounded primitives in the same operation order, so the expected distance is 0 ulp, the bound is
   0, and any difference is a port defect. The sensitivity rule (S, c, K) retires for those functions;
   `validate.py --sensitivity` stays as a diagnostic of conditioning.
4. **Functional parity is unchanged.** The port keeps the oracle's results, `Err` codes, check order and branches.
   A corrected primitive result can change a branch or an `Ok`/`Err` outcome; the corrected oracle follows the
   transcribed code down the other branch, and so must the port.
5. **This supersedes** port plan §5's "the class-C reference is glibc", ADR 0005's glibc-relative S, and ADR 0005's
   M3a amendment: `_km_log10_ln` is replaced by a correctly rounded log10. ADR 0005's correctly rounded primitives
   stay and gain fast paths.
6. **Grids keep the hard rows.** The rows that made the old rule fail — subnormal-times-prefactor results,
   near-cancellation, exact ties — stay in the grids on purpose, and the M2 functions gain them.

## Consequences

- **Positive**
  - The port is at least as accurate as the oracle on every row, and more accurate wherever glibc misrounds.
  - Results are bit-identical on x86_64 ELF, x86_64 PE and aarch64; the platform libm no longer matters.
  - The ill-conditioned rows (cancellation, subnormal scaling) are gated exactly instead of left measured.
- **Negative**
  - kimiya owns correctly rounded math code, with its error analysis and tests, until a Cyrius library provides it
    (`cyrius` issue `2026-11-09-exp-polyfill-double-rounds-subnormal-results`).
  - The transcendental paths cost more than glibc's. The fast paths must keep this within a small factor of the
    stdlib primitives; the double-double `_km_ln` costs 1.9 µs a call today.
  - The tooling now needs the system MPFR library (development only, not the shipped library).
- **Neutral**
  - A Rust 1.x user who compares results will see the last bit differ on the small fraction of calls where glibc
    misrounds (about 0.1–5 % depending on the primitive). The CHANGELOG says so.

## Alternatives considered

- **Keep glibc as the reference** (the M3a amendment's road): rejected by the user. It copies glibc's errors.
- **A wider K or thinner grids:** forbidden by the parity rules, and a K of 10¹⁸ ulp means nothing.
- **Transcribe glibc's libm:** its source is not in the oracle's dependency graph, it has FMA-dispatched variants, and
  it would copy glibc's errors.
- **A fixed-precision decimal reference** (validate.py's former 60-digit one): it misrounds exact ties. MPFR is exact.
