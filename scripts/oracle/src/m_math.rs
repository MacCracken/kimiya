//! `math.tsv` — what Rust **std** returns for every f64 method kimiya 1.2.0's library code
//! calls, the reference for the M0 step 8 baselines (port plan §3 step 5, §6.1).
//!
//! The calls, counted in `rust-old/src` outside the test modules (2026-11-09): sqrt 18,
//! exp 18, max 12, abs 11, powf 10, min 7, ln 7, is_finite 5, round 4, log10 4, cos 3,
//! exp_m1 3, cbrt 2, atan2 2, signum 1, powi 1, clamp 1, acos 1, total_cmp 1. (The `sin`, one
//! more `cos` and one more `max` are `f32`/`usize` calls in the deferred soorat integration;
//! f64 `sin` is recorded anyway, as the port's trigonometry baseline.) `floor` and `trunc`
//! are not called, but they are the building blocks of the port's round-half-away helper,
//! so they are recorded as well.
//!
//! Where each result comes from (rustc 1.99.1, x86_64-unknown-linux-gnu, read off the oracle
//! binary's imports with `nm -D`):
//! - **glibc's libm** (dynamic): exp, exp2, expm1, log, log10, pow, sin, cos, acos, atan2.
//!   These rows are glibc's bits (`glibc.rs`): the header records the glibc, and `--check`
//!   refuses another.
//! - **Rust's compiler-builtins** (static, the pure-Rust `libm` code tanmatra's oracle used):
//!   cbrt (CORE-MATH, correctly rounded), floor, trunc, round, and `__powidf2` for powi.
//! - **inline code**: sqrt (`sqrtsd`), abs, min, max, signum, clamp, total_cmp, is_finite.
//!   min/max lower to `maxsd`/`minsd`, so on (+0, −0) they return the SECOND operand; Rust
//!   leaves that case unspecified, and the rows record what 1.99.1 does.
//!
//! Every input passes through `black_box`, so each row is a run-time library call, never a
//! compile-time fold. kimiya also writes `10.0_f64.powf(y)` (pH, Antoine, activity) and
//! `2.0_f64.powf(y)` (half-lives) with a literal base. LLVM rewrites the second one to
//! `exp2(y)` at opt-level ≥ 1 (the oracle binary imports `exp2`), so the two forms are recorded
//! as written, under `f64::powf_base10` (still glibc `pow`) and `f64::powf_base2` (glibc
//! `exp2`), with only y black-boxed. The corpus comes from the release-built oracle; Rust's own
//! debug-built tests call `pow` for both.
//!
//! Classes:
//! - **A**: exact or correctly rounded on both sides — sqrt, abs, floor, trunc, round
//!   (half away from zero), min, max, signum, clamp, total_cmp, is_finite, and cbrt
//!   (`validate.py` proves every cbrt row correctly rounded with exact arithmetic; ganita's
//!   cbrt is correctly rounded too).
//! - **B**: powi (repeated multiplication, Rust's `__powidf2` order).
//! - **C**: exp, exp_m1, ln, log10, powf (all three forms), cos, sin, acos, atan2.
//!
//! Inputs are the special values, the ranges kimiya feeds each function, and seeded sweeps.

use std::cmp::Ordering;
use std::hint::black_box;

use crate::grids::{Rng, decades, pow2, specials};
use crate::tsv::{Class, Module, b, f, i};

const PI: f64 = core::f64::consts::PI;
const LN_2: f64 = core::f64::consts::LN_2;

fn unary(m: &mut Module, name: &str, class: Class, fun: fn(f64) -> f64, xs: &[f64]) {
    for &x in xs {
        let y = fun(black_box(x));
        m.row(name, class, &[f(x)], &[f(y)]);
    }
}

fn binary(m: &mut Module, name: &str, class: Class, fun: fn(f64, f64) -> f64, xy: &[(f64, f64)]) {
    for &(x, y) in xy {
        let r = fun(black_box(x), black_box(y));
        m.row(name, class, &[f(x), f(y)], &[f(r)]);
    }
}

/// Every ordered pair of `xs`.
fn grid(xs: &[f64]) -> Vec<(f64, f64)> {
    xs.iter()
        .flat_map(|&x| xs.iter().map(move |&y| (x, y)))
        .collect()
}

/// Specials, ±2^k at every 16th binade, and a wide signed log-uniform sweep.
fn wide(seed: u64, n: usize) -> Vec<f64> {
    let mut v = specials();
    for k in (-1074..=1023).step_by(16) {
        let x = pow2(k);
        v.push(x);
        v.push(-x);
    }
    let mut rng = Rng::new(seed);
    v.extend((0..n).map(|_| rng.signed_log_uniform(-300.0, 300.0)));
    v
}

fn sweep(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

fn log_sweep(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| rng.log_uniform(lo_exp10, hi_exp10))
        .collect()
}

fn pairs(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<(f64, f64)> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let x = rng.signed_log_uniform(lo_exp10, hi_exp10);
            (x, rng.signed_log_uniform(lo_exp10, hi_exp10))
        })
        .collect()
}

fn pow10(y: f64) -> f64 {
    10.0_f64.powf(y)
}

fn pow2_of(y: f64) -> f64 {
    2.0_f64.powf(y)
}

/// Rust's `f64::total_cmp` as -1 / 0 / 1.
fn total_cmp(x: f64, y: f64) -> i64 {
    match black_box(x).total_cmp(&black_box(y)) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

#[allow(clippy::too_many_lines)] // one grid per function, in the module comment's order
pub fn emit() -> Module {
    let mut m = Module::new("math");

    // ── exact / correctly rounded (class A) ──
    let mut xs = wide(0x5157, 300);
    xs.extend(log_sweep(0x5158, 300, -40.0, 20.0));
    xs.extend((1..=100).map(|k| f64::from(k * k)));
    xs.extend(sweep(0x5159, 100, 0.0, 1e4));
    unary(&mut m, "f64::sqrt", Class::A, f64::sqrt, &xs);
    unary(&mut m, "f64::abs", Class::A, f64::abs, &wide(0xab5, 100));

    // Ties both ways (±0.5, ±1.5, ±2.5, …), the doubles either side of a tie, the largest
    // double below 0.5, 2^52 − 0.5, then sweeps.
    let mut halves: Vec<f64> = (-20..=20).map(|k| f64::from(k) + 0.5).collect();
    halves.extend((-20..=20).map(|k| f64::from(k) - 0.5 + f64::EPSILON));
    halves.extend([
        0.499_999_999_999_999_94,
        -0.499_999_999_999_999_94,
        4_503_599_627_370_495.5,
        -4_503_599_627_370_495.5,
        4_503_599_627_370_497.0,
    ]);
    halves.extend(sweep(0x7a11, 200, -1e3, 1e3));
    halves.extend(wide(0x7a12, 100));
    unary(&mut m, "f64::floor", Class::A, f64::floor, &halves);
    unary(&mut m, "f64::trunc", Class::A, f64::trunc, &halves);
    unary(&mut m, "f64::round", Class::A, f64::round, &halves);

    // min / max: every special pair (NaN on either side, ±0 both ways, ±∞), equal pairs,
    // a value against 0.0 as kimiya writes `x.max(0.0)`, and random pairs.
    let mut mm = grid(&specials());
    let mut rng = Rng::new(0x3a3a);
    for _ in 0..50 {
        let x = rng.signed_log_uniform(-30.0, 30.0);
        mm.push((x, x));
        mm.push((x, 0.0));
        mm.push((0.0, x));
        mm.push((x, f64::NAN));
    }
    mm.extend(pairs(0x3a3b, 200, -300.0, 300.0));
    binary(&mut m, "f64::max", Class::A, f64::max, &mm);
    binary(&mut m, "f64::min", Class::A, f64::min, &mm);

    unary(
        &mut m,
        "f64::signum",
        Class::A,
        f64::signum,
        &wide(0x5160, 100),
    );

    // clamp: kimiya's one call is `x.clamp(-1.0, 1.0)` (Peng–Robinson's trig root), then
    // other valid bounds (min ≤ max, neither NaN; Rust panics otherwise).
    let mut cl: Vec<f64> = specials();
    cl.extend(sweep(0xc1a0, 300, -2.0, 2.0));
    cl.extend(wide(0xc1a1, 50));
    cl.extend([
        1.0 + f64::EPSILON,
        -1.0 - f64::EPSILON,
        0.999_999_999_999_999_9,
    ]);
    for &x in &cl {
        let r = black_box(x).clamp(-1.0, 1.0);
        m.row("f64::clamp", Class::A, &[f(x), f(-1.0), f(1.0)], &[f(r)]);
    }
    let mut rng = Rng::new(0xc1a2);
    let mut bounds = vec![(0.0, 0.0), (-0.0, 0.0), (-0.0, -0.0), (0.0, f64::INFINITY)];
    for _ in 0..60 {
        let a = rng.signed_log_uniform(-10.0, 10.0);
        let c = rng.signed_log_uniform(-10.0, 10.0);
        bounds.push(if a <= c { (a, c) } else { (c, a) });
    }
    for (k, &(lo, hi)) in bounds.iter().enumerate() {
        let probes = [
            lo,
            hi,
            0.0,
            -0.0,
            f64::NAN,
            rng.signed_log_uniform(-10.0, 10.0),
        ];
        for &x in probes.iter().take(if k < 4 { 6 } else { 2 + k % 4 }) {
            let r = black_box(x).clamp(black_box(lo), black_box(hi));
            m.row("f64::clamp", Class::A, &[f(x), f(lo), f(hi)], &[f(r)]);
        }
    }

    // total_cmp: the specials plus NaNs of both signs and a payload, every pair.
    let mut tc = specials();
    tc.extend([
        f64::from_bits(0xfff8_0000_0000_0000), // -NaN
        f64::from_bits(0x7ff0_0000_0000_0001), // signalling-pattern NaN
        -f64::MIN_POSITIVE,
        -f64::from_bits(1),
    ]);
    for (x, y) in grid(&tc) {
        m.row(
            "f64::total_cmp",
            Class::A,
            &[f(x), f(y)],
            &[i(total_cmp(x, y))],
        );
    }
    for (x, y) in pairs(0x7c70, 60, -300.0, 300.0) {
        m.row(
            "f64::total_cmp",
            Class::A,
            &[f(x), f(y)],
            &[i(total_cmp(x, y))],
        );
    }

    let mut fin = specials();
    fin.extend([
        f64::from_bits(0xfff8_0000_0000_0000),
        -f64::MAX,
        f64::from_bits(0x7fef_ffff_ffff_ffff),
    ]);
    fin.extend(log_sweep(0xf1f1, 20, -300.0, 300.0));
    for &x in &fin {
        m.row(
            "f64::is_finite",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_finite())],
        );
    }

    // cbrt: the mass numbers A = 1..300 (SEMF, A^(1/3)), perfect cubes of both signs, the
    // Peng–Robinson cubic's range, and a wide sweep.
    let mut cb: Vec<f64> = (1..=300).map(f64::from).collect();
    cb.extend((1..=1000).map(|k| f64::from(k * k * k)));
    cb.extend((1..=100).map(|k| -f64::from(k * k * k)));
    cb.extend(log_sweep(0xcb26, 200, -30.0, 30.0));
    cb.extend(wide(0xcb27, 400));
    unary(&mut m, "f64::cbrt", Class::A, f64::cbrt, &cb);

    // ── arithmetic (class B) ──
    // kimiya's one call is `.powi(2)` (Peng–Robinson α). n is black-boxed too, so every row
    // is a `__powidf2` call; for n = 2 that is x·x, as the inline form.
    let mut pw = specials();
    pw.extend(sweep(0x9101, 100, -10.0, 10.0));
    pw.extend(log_sweep(0x9102, 60, -150.0, 150.0));
    for n in [2, -3, -2, -1, 0, 1, 3, 4, 7, 10] {
        let set: &[f64] = if n == 2 { &pw } else { &pw[..80] };
        for &x in set {
            let r = black_box(x).powi(black_box(n));
            m.row("f64::powi", Class::B, &[f(x), i(n)], &[f(r)]);
        }
    }

    // ── transcendental (class C) ──
    // exp: the Arrhenius / Boltzmann factor −Ea/RT in [−700, 0] above all, the full finite
    // range [−745.13, 709.78] and its edges, integers across it, and near 0.
    let mut ex = specials();
    ex.extend(sweep(0xe001, 400, -700.0, 0.0));
    ex.extend(sweep(0xe002, 200, -745.2, 709.8));
    ex.extend(sweep(0xe003, 200, -50.0, 50.0));
    ex.extend(sweep(0xe004, 100, -1e-6, 1e-6));
    ex.extend((-745..=709).step_by(7).map(f64::from));
    ex.extend([
        709.782_712_893_384,
        709.782_712_893_384_1,
        -708.396_418_532_264_1,
        -745.133_219_101_941_1,
        -745.133_219_101_941_2,
        1e-300,
        -1e-300,
    ]);
    unary(&mut m, "f64::exp", Class::C, f64::exp, &ex);

    // exp_m1: kimiya's three calls are (n − 1)·ln 2 (nth-order half-life) and −Δλ·t, −λ·t
    // (decay chains: tiny to large negative arguments), plus the near-0 region it exists for.
    let mut em = specials();
    em.extend(sweep(0xe101, 100, -1e-6, 1e-6));
    em.extend(sweep(0xe102, 200, -1.0, 1.0));
    em.extend(sweep(0xe103, 100, -50.0, 50.0));
    em.extend(sweep(0xe104, 100, -745.2, 709.8));
    em.extend(log_sweep(0xe105, 200, -20.0, 3.0).into_iter().map(|x| -x));
    em.extend(
        sweep(0xe106, 100, -10.0, 10.0)
            .into_iter()
            .map(|n| n * LN_2),
    );
    unary(&mut m, "f64::exp_m1", Class::C, f64::exp_m1, &em);

    // ln: concentrations, rate constants and pressures 1e-40 .. 1e20, both as sweeps and as
    // decades of 1, 2 and 5, plus the region around 1 and a wide sweep.
    let mut lg = wide(0x1060, 300);
    lg.extend(sweep(0x1061, 200, 0.5, 2.0));
    lg.extend(sweep(0x1062, 100, 1.0 - 1e-9, 1.0 + 1e-9));
    lg.extend(log_sweep(0x1063, 300, -40.0, 20.0));
    lg.extend(decades(-40, 20, &["1", "2", "5"]));
    unary(&mut m, "f64::ln", Class::C, f64::ln, &lg);

    // log10: pH and pOH ([H⁺] 1e-15 .. 10), Henderson–Hasselbalch ratios, absorbance, every
    // exact decade 1e-40 .. 1e20, and the same wide sweep as ln.
    let mut l10 = wide(0x1070, 300);
    l10.extend(log_sweep(0x1071, 300, -15.0, 1.0));
    l10.extend(log_sweep(0x1072, 200, -40.0, 20.0));
    l10.extend(sweep(0x1073, 100, 1.0 - 1e-9, 1.0 + 1e-9));
    l10.extend(decades(-40, 20, &["1", "3"]));
    unary(&mut m, "f64::log10", Class::C, f64::log10, &l10);

    // powf, both operands run-time: Tc^2.5 (Peng–Robinson / Redlich–Kwong a), [A]0^(n−1)
    // (nth-order half-life), random positive pairs, the special grid, and negative bases with
    // integer exponents.
    let mut pf: Vec<(f64, f64)> = sweep(0x9001, 100, 5.0, 1000.0)
        .into_iter()
        .map(|t| (t, 2.5))
        .collect();
    let mut rng = Rng::new(0x9002);
    for _ in 0..200 {
        let x = rng.log_uniform(-10.0, 3.0);
        pf.push((x, rng.uniform(-5.0, 5.0)));
    }
    for _ in 0..300 {
        let x = rng.log_uniform(-30.0, 30.0);
        pf.push((x, rng.uniform(-10.0, 10.0)));
    }
    pf.extend(grid(&specials()));
    for _ in 0..50 {
        let x = -rng.uniform(0.1, 10.0);
        #[allow(clippy::cast_precision_loss)]
        let y = rng.int(-5, 5) as f64;
        pf.push((x, y));
    }
    binary(&mut m, "f64::powf", Class::C, f64::powf, &pf);

    // 10.0_f64.powf(y), as kimiya writes it: pH → [H⁺] (−ph), activity coefficients
    // (log γ in [−5, 0.5]), Antoine (A − B/(C + T)), transmittance (−A); every integer
    // exponent −40 ..= 20; and the full range.
    let mut p10 = specials();
    p10.extend((-40..=20).map(f64::from));
    p10.extend(sweep(0x9010, 300, -16.0, 2.0));
    p10.extend(sweep(0x9011, 100, -5.0, 0.5));
    p10.extend(sweep(0x9012, 100, -325.0, 310.0));
    unary(&mut m, "f64::powf_base10", Class::C, pow10, &p10);

    // 2.0_f64.powf(y), as kimiya writes it: the remaining fraction 2^(−n) after n half-lives,
    // the LJ minimum 2^(1/6), integer exponents across the range, and sweeps.
    let mut p2 = specials();
    p2.push(1.0 / 6.0);
    p2.extend((-1080..=1030).step_by(37).map(f64::from));
    p2.extend(sweep(0x9020, 200, -100.0, 10.0));
    p2.extend(sweep(0x9021, 100, -60.0, 0.0));
    unary(&mut m, "f64::powf_base2", Class::C, pow2_of, &p2);

    // cos: Peng–Robinson's three trig roots take cos(θ/3) and cos((θ ± 2π)/3) for θ in [0, π],
    // i.e. [−2π/3, π]; then the general sweeps sin shares.
    let mut tr = specials();
    tr.extend(sweep(0x5101, 400, -2.0 * PI, 2.0 * PI));
    tr.extend(sweep(0x5102, 200, 0.0, 1e3));
    tr.extend(log_sweep(0x5103, 100, -12.0, 15.0));
    let mut co = tr.clone();
    co.extend(sweep(0x5104, 200, -2.0 * PI / 3.0, PI));
    unary(&mut m, "f64::cos", Class::C, f64::cos, &co);
    unary(&mut m, "f64::sin", Class::C, f64::sin, &tr);

    // acos: kimiya's argument is clamped to [−1, 1]; the endpoints, 1 − 2^−k toward each end,
    // just outside (NaN), and tiny values.
    let mut ac = specials();
    ac.extend(sweep(0xac01, 400, -1.0, 1.0));
    for k in 1..=53 {
        let d = pow2(-k);
        ac.push(1.0 - d);
        ac.push(d - 1.0);
    }
    ac.extend([1.0 + f64::EPSILON, -1.0 - f64::EPSILON, 2.0, -2.0]);
    ac.extend(log_sweep(0xac02, 100, -300.0, -1.0));
    unary(&mut m, "f64::acos", Class::C, f64::acos, &ac);

    // atan2(y, x): bond angles (y = |a × b| ≥ 0) and dihedrals (any sign), at ångström and
    // metre scale (coordinates 1e-10 m, so products near 1e-20); every quadrant; the special
    // grid (±0 and ±∞ in each slot).
    let mut at = grid(&specials());
    at.extend(pairs(0xa7a1, 400, -10.0, 10.0));
    let mut rng = Rng::new(0xa7a2);
    for _ in 0..200 {
        let y = rng.uniform(-1.0, 1.0);
        at.push((y, rng.uniform(-1.0, 1.0)));
    }
    at.extend(pairs(0xa7a3, 200, -30.0, 5.0));
    for _ in 0..100 {
        let y = rng.log_uniform(-25.0, 2.0);
        at.push((y, rng.signed_log_uniform(-25.0, 2.0)));
    }
    binary(&mut m, "f64::atan2", Class::C, f64::atan2, &at);

    m
}
