//! `fitting.tsv` — parameter estimation from measurements (port plan §4, M5b): Arrhenius fits, polynomial least squares,
//! Beer–Lambert calibration lines, and `PolynomialFit`'s two methods.
//!
//! - Class C (`ln` per point and `exp` of the intercept, glibc's in the oracle; ADR 0006): `fit_arrhenius` (input the
//!   temperatures `v:n` + `f:`, then the rate constants `v:n` + `f:`; output `ArrheniusFit`'s fields
//!   `f:pre_exponential f:activation_energy_j`, or `E:`): exact synthetic data (k = 2^e exactly, 1/T linear in e, so
//!   ln k is exactly linear in the 1/T the generator rounded once), the same with seeded multiplicative noise, two
//!   points, many points, temperatures across hisab's 10⁻¹² rank threshold (two temperatures δ apart: the second QR
//!   column's norm δ/(√2·T²)), equal temperatures (ComputationError), subnormal and huge temperatures (1/T infinite or
//!   subnormal), rate constants from the smallest subnormal to DBL_MAX, every Err in Rust's order (the lengths, fewer
//!   than 2 points, then per pair the temperature before the rate constant).
//! - Class B (+ − × ÷, √ and `powi`, hisab 1.4.0's `least_squares_poly` transcribed): `fit_polynomial` (input `v:n`
//!   x, `v:n` y, `i:degree` as a usize word; output `v:degree+1` + `f:` coefficients, or `E:`): the Rust tests' lines,
//!   exact integer polynomials of degree 0 … 7 at integer and dyadic x (fitted at their degree, one above and one
//!   below), seeded noise, ill-conditioned Vandermonde systems (x near 1000, near 10⁶, many nodes in [0, 1] at degree
//!   up to 15), degenerate ones (equal x, fewer distinct x than coefficients, x so small the columns vanish or so large
//!   their powers overflow), every Err (the lengths, too few points for the degree, usize::MAX, NaN and ±∞ in x or
//!   y); `fit_beer_lambert` (input `v:n` concentrations, `v:n` absorbances, `f:path_length`; output `f:ε`, or `E:`): the
//!   Rust tests, the blank's intercept, path lengths over the edges, seeded calibration lines; `PolynomialFit::evaluate`
//!   (input `v:n` coefficients, `f:x`; output `f:`): the Rust test, empty and seeded coefficient lists over x's edges
//!   (x^k overflowing against a zero coefficient: ∞ · 0 = NaN).
//! - Class A: `PolynomialFit::degree` (input `v:n` coefficients; output `i:`): lengths 0 … 9.
//!
//! Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::grids::{Rng, pow2};
use crate::systems::vec_toks;
use crate::tsv::{Class, Module, Tok, err, f, i};
use kimiya::fitting::{PolynomialFit, fit_arrhenius, fit_beer_lambert, fit_polynomial};

/// The harness's token list holds 8448 tokens (`tests/parity_core.cyr`): every row must fit.
const MAX_TOKENS: usize = 8448;

fn word(x: usize) -> Tok {
    i(i64::from_ne_bytes((x as u64).to_ne_bytes()))
}

fn arrhenius(m: &mut Module, temps: &[f64], rates: &[f64]) {
    let mut ins = vec_toks(temps);
    ins.extend(vec_toks(rates));
    let outs = match fit_arrhenius(bb(temps), bb(rates)) {
        Ok(fit) => vec![f(fit.pre_exponential), f(fit.activation_energy_j)],
        Err(e) => vec![err(&e)],
    };
    assert!(ins.len() <= MAX_TOKENS, "fit_arrhenius: a row too long");
    m.row("fitting::fit_arrhenius", Class::C, &ins, &outs);
}

fn poly(m: &mut Module, x: &[f64], y: &[f64], degree: usize) {
    let mut ins = vec_toks(x);
    ins.extend(vec_toks(y));
    ins.push(word(degree));
    let outs = match fit_polynomial(bb(x), bb(y), bb(degree)) {
        Ok(fit) => vec_toks(&fit.coefficients),
        Err(e) => vec![err(&e)],
    };
    assert!(ins.len() <= MAX_TOKENS, "fit_polynomial: a row too long");
    m.row("fitting::fit_polynomial", Class::B, &ins, &outs);
}

fn beer(m: &mut Module, c: &[f64], a: &[f64], l: f64) {
    let mut ins = vec_toks(c);
    ins.extend(vec_toks(a));
    ins.push(f(l));
    let outs = match fit_beer_lambert(bb(c), bb(a), bb(l)) {
        Ok(e) => vec![f(e)],
        Err(e) => vec![err(&e)],
    };
    m.row("fitting::fit_beer_lambert", Class::B, &ins, &outs);
}

fn evaluate(m: &mut Module, coeffs: &[f64], x: f64) {
    let fit = PolynomialFit { coefficients: coeffs.to_vec() };
    let mut ins = vec_toks(coeffs);
    ins.push(f(x));
    m.row("fitting::PolynomialFit::evaluate", Class::B, &ins, &[f(bb(&fit).evaluate(bb(x)))]);
}

fn degree(m: &mut Module, coeffs: &[f64]) {
    let fit = PolynomialFit { coefficients: coeffs.to_vec() };
    m.row("fitting::PolynomialFit::degree", Class::A, &vec_toks(coeffs), &[word(bb(&fit).degree())]);
}

/// Integer-coefficient polynomial evaluation in Horner form (exact for the small integers used).
fn horner(c: &[f64], x: f64) -> f64 {
    c.iter().rev().fold(0.0, |acc, &k| acc * x + k)
}

fn usz(x: i64) -> usize {
    usize::try_from(x).unwrap_or(0)
}

fn small(x: i64) -> f64 {
    f64::from(i32::try_from(x).unwrap_or(0))
}

#[allow(clippy::too_many_lines)] // one function after another
pub fn emit() -> Module {
    let mut m = Module::new("fitting");
    let mut rng = Rng::new(0x5b02);

    // ── fit_arrhenius ──
    // Exact synthetic data: k_i = 2^(e0 − i·de) exactly, 1/T_i = x0 + i·dx rounded once (T_i = 1/(x0 + i·dx)).
    for &(e0, de, x0, dx, n) in &[
        (40_i32, 4_i32, 1.0 / 600.0, 1e-4, 5_usize),
        (20, 1, 1.0 / 300.0, 2e-5, 6),
        (-10, 2, 1.0 / 1000.0, 5e-5, 8),
        (0, 3, 1.0 / 250.0, 1e-4, 2),
        (100, 10, 1.0 / 2000.0, 1e-5, 12),
        (-900, 5, 1e-3, 1e-5, 4),
        (1000, 20, 1e-4, 2e-6, 10),
        (10, 1, 0.01, 0.001, 40),
        (60, 1, 1.0 / 400.0, 1e-6, 120),
    ] {
        let temps: Vec<f64> = (0..n).map(|k| 1.0 / (x0 + dx * f64::from(u32::try_from(k).unwrap_or(0)))).collect();
        let rates: Vec<f64> = (0..n).map(|k| pow2(e0 - de * i32::try_from(k).unwrap_or(0))).collect();
        arrhenius(&mut m, &temps, &rates);
        // Seeded multiplicative noise of ±5 % and ±0.1 %.
        for amp in [0.05, 0.001] {
            let noisy: Vec<f64> = rates.iter().map(|&k| k * (1.0 + rng.uniform(-amp, amp))).collect();
            arrhenius(&mut m, &temps, &noisy);
        }
        // Reversed order.
        let rt: Vec<f64> = temps.iter().rev().copied().collect();
        let rr: Vec<f64> = rates.iter().rev().copied().collect();
        arrhenius(&mut m, &rt, &rr);
    }
    // Decimal temperatures and rates (each parsed, correctly rounded).
    let temps = [300.0, 320.0, 340.0, 360.0, 380.0, 400.0];
    arrhenius(&mut m, &temps, &[1e-5, 5e-5, 2.2e-4, 8.4e-4, 2.9e-3, 9.1e-3]);
    arrhenius(&mut m, &[300.0, 350.0, 400.0, 450.0, 500.0], &[1.0, 10.0, 60.0, 250.0, 800.0]);
    arrhenius(&mut m, &[300.0, 400.0], &[1.0, 2.0]);
    arrhenius(&mut m, &[300.0, 400.0], &[2.0, 1.0]); // a negative activation energy
    arrhenius(&mut m, &[300.0, 400.0], &[1.0, 1.0]); // Ea = 0
    // The Rust tests' errors and every other Err, in Rust's order.
    arrhenius(&mut m, &[300.0], &[1.0]);
    arrhenius(&mut m, &[], &[]);
    arrhenius(&mut m, &[300.0, 400.0], &[1.0]);
    arrhenius(&mut m, &[300.0], &[1.0, 2.0]);
    arrhenius(&mut m, &[300.0, -100.0], &[1.0, 2.0]);
    arrhenius(&mut m, &[300.0, 400.0], &[1.0, 0.0]);
    for bad in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        arrhenius(&mut m, &[bad, 400.0, 500.0], &[1.0, 2.0, 3.0]);
        arrhenius(&mut m, &[300.0, 400.0, bad], &[1.0, 2.0, 3.0]);
        arrhenius(&mut m, &[300.0, 400.0, 500.0], &[bad, 2.0, 3.0]);
        arrhenius(&mut m, &[300.0, 400.0, 500.0], &[1.0, 2.0, bad]);
        arrhenius(&mut m, &[300.0, bad, 500.0], &[bad, 2.0, 3.0]); // the rate constant of pair 0 first
        arrhenius(&mut m, &[300.0, bad, 500.0], &[1.0, bad, 3.0]); // the temperature of pair 1 first
    }
    // Equal temperatures (a rank-deficient QR), and two temperatures δ apart across the threshold: the second
    // column's norm after orthogonalization is |1/T₁ − 1/T₂|/√2 ≈ δ/(√2·T²), below 10⁻¹² for δ ≲ 1.27·10⁻⁷ at 300 K.
    arrhenius(&mut m, &[300.0, 300.0], &[1.0, 2.0]);
    arrhenius(&mut m, &[300.0, 300.0, 300.0], &[1.0, 2.0, 3.0]);
    for d in [1e-9, 1e-8, 5e-8, 1e-7, 1.2e-7, 1.25e-7, 1.3e-7, 1.5e-7, 2e-7, 1e-6, 1e-5, 1e-3] {
        arrhenius(&mut m, &[300.0, 300.0 + d], &[1.0, 1.001]);
        arrhenius(&mut m, &[300.0, 300.0 + d, 300.0 + 2.0 * d], &[1.0, 1.001, 1.002]);
    }
    // Extreme temperatures and rate constants.
    for tt in [pow2(-1074), f64::MIN_POSITIVE, 1e-300, 1e-10, 1e10, 1e300, f64::MAX] {
        arrhenius(&mut m, &[300.0, tt], &[1.0, 2.0]);
        arrhenius(&mut m, &[tt, 2.0 * tt], &[1.0, 2.0]);
    }
    for k in [pow2(-1074), f64::MIN_POSITIVE, 1e-300, 1e300, f64::MAX] {
        arrhenius(&mut m, &[300.0, 400.0, 500.0], &[k, 1.0, 2.0]);
        arrhenius(&mut m, &[300.0, 400.0], &[k, k]);
    }
    arrhenius(&mut m, &[300.0, 400.0], &[pow2(-1074), f64::MAX]); // A overflows
    for _ in 0..40 {
        let n = usz(rng.int(2, 30));
        let temps: Vec<f64> = (0..n).map(|_| rng.uniform(200.0, 2500.0)).collect();
        let rates: Vec<f64> = (0..n).map(|_| rng.log_uniform(-12.0, 12.0)).collect();
        arrhenius(&mut m, &temps, &rates);
    }

    // ── fit_polynomial ──
    // The Rust tests.
    let x5 = [1.0, 2.0, 3.0, 4.0, 5.0];
    poly(&mut m, &x5, &x5.map(|xi| 2.0 + 3.0 * xi), 1);
    let x10: Vec<f64> = (0..10).map(f64::from).collect();
    let y10: Vec<f64> = x10.iter().map(|&xi| 1.0 + 2.0 * xi + 0.5 * xi * xi).collect();
    poly(&mut m, &x10, &y10, 2);
    poly(&mut m, &[1.0, 2.0], &[1.0, 2.0], 2);
    poly(&mut m, &[1.0, 2.0], &[1.0, 2.0], usize::MAX);
    poly(&mut m, &[0.0, 1.0, 2.0], &[1.0, f64::NAN, 3.0], 1);
    poly(&mut m, &[0.0, f64::INFINITY, 2.0], &[1.0, 2.0, 3.0], 1);
    poly(&mut m, &[f64::NAN, 1.0, 2.0], &[1.0, 2.0, 3.0], 1);
    // The bench's data: degree 2 at 20 points 0.0 … 1.9.
    let xb: Vec<f64> = (0..20).map(|k| f64::from(k) * 0.1).collect();
    let yb: Vec<f64> = xb.iter().map(|&xi| 1.0 + 2.0 * xi + 0.5 * xi * xi).collect();
    poly(&mut m, &xb, &yb, 2);
    // Every Err: the lengths, too few points, non-finite values anywhere.
    poly(&mut m, &[1.0, 2.0, 3.0], &[1.0, 2.0], 1);
    poly(&mut m, &[], &[], 0);
    poly(&mut m, &[], &[1.0], 0);
    poly(&mut m, &[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0], 3);
    poly(&mut m, &[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0], usize::MAX - 1);
    poly(&mut m, &[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0], 1 << 63);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        poly(&mut m, &[bad, 1.0, 2.0], &[1.0, 2.0, 3.0], 0);
        poly(&mut m, &[0.0, 1.0, 2.0], &[1.0, 2.0, bad], 2);
        poly(&mut m, &[0.0, 1.0, bad], &[bad, 2.0, 3.0], 5); // too few points comes first
    }
    // Exact integer polynomials, fitted at their degree, one above and one below.
    for d in 0..=7_usize {
        for dyadic in [false, true] {
            let coeffs: Vec<f64> = (0..=d).map(|_| small(rng.int(-5, 5))).collect();
            let n = d + 1 + usz(rng.int(0, 12));
            let x: Vec<f64> = (0..n)
                .map(|k| {
                    let v = small(i64::try_from(k).unwrap_or(0)) - small(i64::try_from(n / 2).unwrap_or(0));
                    if dyadic { v * 0.25 } else { v }
                })
                .collect();
            let y: Vec<f64> = x.iter().map(|&xi| horner(&coeffs, xi)).collect();
            poly(&mut m, &x, &y, d);
            if n > d + 1 {
                poly(&mut m, &x, &y, d + 1);
            }
            if d > 0 {
                poly(&mut m, &x, &y, d - 1);
            }
            // Seeded noise.
            let yn: Vec<f64> = y.iter().map(|&v| v + rng.uniform(-0.01, 0.01)).collect();
            poly(&mut m, &x, &yn, d);
        }
    }
    // The interpolating degree (n points, degree n − 1) and degree 0 (the QR's mean).
    for n in 1..=8_usize {
        let x: Vec<f64> = (0..n).map(|_| rng.uniform(-3.0, 3.0)).collect();
        let y: Vec<f64> = (0..n).map(|_| rng.uniform(-10.0, 10.0)).collect();
        poly(&mut m, &x, &y, n - 1);
        poly(&mut m, &x, &y, 0);
    }
    // Ill-conditioned Vandermonde systems.
    for d in [2_usize, 3, 4, 5, 6] {
        let x: Vec<f64> = (0..21).map(|k| 1000.0 + f64::from(k) * 0.5).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 1.0 + 0.001 * xi).collect();
        poly(&mut m, &x, &y, d);
    }
    for d in [1_usize, 2, 3] {
        let x: Vec<f64> = (0..10).map(|k| 1e6 + f64::from(k)).collect();
        let y: Vec<f64> = (0..10).map(|k| f64::from(k) * 2.0).collect();
        poly(&mut m, &x, &y, d);
    }
    for d in [5_usize, 8, 10, 12, 15] {
        let x: Vec<f64> = (0..25).map(|k| f64::from(k) / 24.0).collect();
        let y: Vec<f64> = x.iter().map(|&xi| xi * xi - 0.5).collect();
        poly(&mut m, &x, &y, d);
    }
    // Degenerate: equal x, two distinct x at degree 2, vanishing and overflowing columns, zero data.
    poly(&mut m, &[2.0, 2.0, 2.0], &[1.0, 2.0, 3.0], 1);
    poly(&mut m, &[2.0, 2.0, 2.0], &[1.0, 2.0, 3.0], 0);
    poly(&mut m, &[1.0, 1.0, 1.0, 2.0], &[1.0, 2.0, 3.0, 4.0], 2);
    poly(&mut m, &[1e6, 1e6, 1e6, 2e6], &[1.0, 2.0, 3.0, 4.0], 2);
    for s in [1e-100, 1e-160, 1e-200, 1e100, 1e154, 1e160, 1e200, 1e300] {
        let x: Vec<f64> = (1..=6).map(|k| f64::from(k) * s).collect();
        let y: Vec<f64> = (1..=6).map(f64::from).collect();
        poly(&mut m, &x, &y, 1);
        poly(&mut m, &x, &y, 2);
        poly(&mut m, &x, &y, 4);
    }
    poly(&mut m, &[1.0, 2.0, 3.0], &[0.0, 0.0, 0.0], 1);
    poly(&mut m, &[1.0, 2.0, 3.0], &[-0.0, -0.0, -0.0], 2);
    poly(&mut m, &[-1.0, 0.0, 1.0], &[1e300, -1e300, 1e300], 2);
    poly(&mut m, &[-1.0, 0.0, 1.0], &[f64::MAX, f64::MAX, f64::MAX], 1);
    poly(&mut m, &[0.0, -0.0, 1.0], &[1.0, 2.0, 3.0], 1);
    for _ in 0..30 {
        let n = usz(rng.int(2, 60));
        let d = usz(rng.int(0, 4)).min(n - 1);
        let x: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(-2.0, 3.0)).collect();
        let y: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect();
        poly(&mut m, &x, &y, d);
    }

    // ── fit_beer_lambert ──
    let conc = [0.01, 0.02, 0.03, 0.04, 0.05];
    beer(&mut m, &conc, &conc.map(|c| 100.0 * c), 1.0);
    let conc4 = [0.001, 0.002, 0.003, 0.004];
    beer(&mut m, &conc4, &conc4.map(|c| 200.0 * 2.0 * c), 2.0);
    beer(&mut m, &[0.01], &[1.0], 0.0);
    let conc6 = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05];
    beer(&mut m, &conc6, &conc6.map(|c| 150.0 * c + 0.005), 1.0);
    let conc3 = [0.01, 0.02, 0.03];
    for l in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -0.0, -1.0, pow2(-1074), 1e-300, 0.1, 10.0, 1e300] {
        beer(&mut m, &conc3, &[1.0, 2.0, 3.0], l);
    }
    beer(&mut m, &conc3, &[1.0, f64::NAN, 3.0], 1.0);
    beer(&mut m, &conc3, &[1.0, 2.0], 1.0);
    beer(&mut m, &[0.01], &[1.0], 1.0);
    beer(&mut m, &[0.02, 0.02, 0.02], &[1.0, 1.1, 1.2], 1.0);
    beer(&mut m, &[f64::NAN, 0.02], &[1.0, 1.1], f64::NAN); // the path length first
    for _ in 0..30 {
        let n = usz(rng.int(2, 12));
        let eps = rng.log_uniform(0.0, 5.0);
        let l = rng.log_uniform(-1.0, 1.0);
        let blank = rng.uniform(-0.01, 0.01);
        let c: Vec<f64> = (0..n).map(|_| rng.log_uniform(-6.0, -2.0)).collect();
        let a: Vec<f64> = c.iter().map(|&ci| eps * l * ci + blank + rng.uniform(-1e-4, 1e-4)).collect();
        beer(&mut m, &c, &a, l);
    }

    // ── PolynomialFit::evaluate, ::degree ──
    for x in [0.0, 1.0, 2.0] {
        evaluate(&mut m, &[1.0, 2.0, 3.0], x); // the Rust test
    }
    for x in [0.0, -0.0, 1.0, -1.0, 0.5, 2.0, -3.0, 1e100, 1e200, 1e308, f64::MAX, pow2(-1074), f64::INFINITY,
              f64::NEG_INFINITY, f64::NAN] {
        evaluate(&mut m, &[], x);
        evaluate(&mut m, &[1.0], x);
        evaluate(&mut m, &[1.0, 0.0], x);
        evaluate(&mut m, &[0.0, 0.0, 1.0], x);
        evaluate(&mut m, &[-0.0], x);
        evaluate(&mut m, &[1.0, -2.0, 0.5, 3.0], x);
        evaluate(&mut m, &[f64::NAN, 1.0], x);
    }
    for _ in 0..60 {
        let n = usz(rng.int(0, 9));
        let c: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect();
        let x = rng.signed_log_uniform(-2.0, 2.0);
        evaluate(&mut m, &c, x);
    }
    for n in 0..=9 {
        degree(&mut m, &vec![1.0; n]);
    }
    m
}
