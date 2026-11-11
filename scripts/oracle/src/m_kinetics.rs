//! `kinetics.tsv` — rate constants, integrated rate laws, half-lives, the steady-state and pre-equilibrium
//! approximations, Michaelis–Menten, collision theory and the Eyring equation (port plan §4, M3b).
//!
//! - Class B (+ − × ÷ only, expected 0 ulp): `half_life_first_order` (Rust's `LN_2 / k`: every decade of k from
//!   the smallest subnormal to DBL_MAX, k = 1 giving LN_2 exactly), `second_order_concentration`,
//!   `zero_order_concentration` (NaN propagates, −0.0 stays −0.0), `zero_order_half_life`,
//!   `steady_state_intermediate`, `steady_state_rate`, `pre_equilibrium_rate`, `michaelis_menten` and
//!   `lineweaver_burk` (output `f:slope f:intercept`): each argument over the validation edges with the others
//!   typical, every Err path, decade grids, overflow, the Rust tests' values and seeded sweeps.
//! - Class C, one `exp`: `arrhenius_rate`, `first_order_concentration`, `collision_theory_rate`, `eyring_rate`;
//!   two: `eyring_rate_from_activation`. T over 1–5000 K (port plan §4) and beyond (the smallest subnormal,
//!   10⁻³⁰⁰ … 10³⁰⁰, DBL_MAX); Ea, ΔG‡, ΔH‡, ΔS‡ over wide signed ranges; the exponent placed on targets x across
//!   exp's whole range (Ea = −x·R·T); and **on purpose the rows where the exp result is subnormal or 0 and is
//!   then multiplied by a prefactor** (A = 10¹³ … DBL_MAX with x between −745.2 and −708, Eyring at low T, [A]₀
//!   up to 10³⁰⁸ with k·t in that range, Z·p), and exp's overflow (A·∞, and 0·∞ or ∞·0 → NaN, as Rust gives it).
//! - ADR 0002 #28 (the overflow audit's B1): `second_order_concentration`, `zero_order_half_life`,
//!   `steady_state_intermediate`, `pre_equilibrium_rate` and `michaelis_menten` rows are noted where the port's value
//!   (the correctly rounded monomial, or Rust's sequence with an unbounded exponent; `audit.rs`) differs from Rust's,
//!   and `audit::emit_kinetics` appends their edge rows with their own seed.
//! - Class C, `exp_m1` and `powf`: `nth_order_half_life` — orders on a fine grid around 1 (the 10⁻¹² tolerance
//!   and the ulps either side of it, from both sides), integer and fractional orders out to ±300, [A]₀ and k over
//!   many decades, the extreme-order overflow and underflow Err rows, **rows where [A]₀^(n−1) is subnormal and the
//!   half-life is still finite** (the pow result divided into the normal range), the Rust tests' cases and seeded
//!   sweeps.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::audit;
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, err, f};
use kimiya::element::GAS_CONSTANT;
use kimiya::kinetics::{
    arrhenius_rate, collision_theory_rate, eyring_rate, eyring_rate_from_activation, first_order_concentration,
    half_life_first_order, lineweaver_burk, nth_order_half_life, steady_state_rate, zero_order_concentration,
};

/// `x` moved `d` ulps (up for d > 0).
fn ulps(x: f64, d: i32) -> f64 {
    let mut y = x;
    for _ in 0..d.unsigned_abs() {
        y = if d > 0 { y.next_up() } else { y.next_down() };
    }
    y
}

/// Temperatures: 1–5000 K (port plan §4) and beyond.
fn temperatures() -> Vec<f64> {
    vec![
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1e-300,
        1e-10,
        0.01,
        1.0,
        2.0,
        5.0,
        10.0,
        20.0,
        50.0,
        77.0,
        100.0,
        200.0,
        273.15,
        298.0,
        298.15,
        300.0,
        400.0,
        500.0,
        750.0,
        1000.0,
        1500.0,
        2000.0,
        3000.0,
        4000.0,
        5000.0,
        1e10,
        1e300,
        f64::MAX,
    ]
}

/// Target exponents x = −Ea/(RT): across exp's underflow (x ≈ −745.13) and overflow (x ≈ 709.78) ends, through
/// the subnormal results (−745.13 < x < −708.40) and zero.
const TARGETS: [f64; 31] = [
    -800.0, -746.0, -745.2, -745.14, -745.13, -745.1, -745.0, -744.0, -740.0, -735.5, -730.0, -725.0, -720.25,
    -715.0, -710.0, -708.5, -708.4, -708.39, -708.0, -700.0, -100.0, -1.0, -1e-10, 0.0, 1e-10, 1.0, 100.0, 700.0,
    709.7, 709.79, 710.0,
];

/// A Result<(f64, f64)> as tokens: `f:slope f:intercept`, or `E:Variant`.
fn pair(r: Result<(f64, f64), kimiya::error::KimiyaError>) -> Vec<crate::tsv::Tok> {
    match r {
        Ok((a, b)) => vec![f(a), f(b)],
        Err(e) => vec![err(&e)],
    }
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("kinetics");
    let mut rng = Rng::new(0x3c01);
    let r = GAS_CONSTANT;

    // ── arrhenius_rate ──
    let arr = |m: &mut Module, a: f64, ea: f64, t: f64| {
        let out = res(arrhenius_rate(bb(a), bb(ea), bb(t)));
        m.row("kinetics::arrhenius_rate", Class::C, &[f(a), f(ea), f(t)], &[out]);
    };
    for x in edges() {
        arr(&mut m, x, 50_000.0, 300.0);
        arr(&mut m, 1e13, x, 300.0);
        arr(&mut m, 1e13, 50_000.0, x);
    }
    // the Rust tests' values (and the order of the checks: T, then A, then Ea)
    for (a, ea, t) in [
        (1e13, 100_000.0, 300.0),
        (1e13, 50_000.0, 300.0),
        (1e13, 50_000.0, 400.0),
        (1e13, 30_000.0, 300.0),
        (1e13, 50_000.0, 0.0),
        (1e13, 50_000.0, -100.0),
        (1e13, 0.0, 300.0),
        (1e13, 50_000.0, f64::NAN),
        (1e13, 50_000.0, f64::INFINITY),
        (f64::NAN, 50_000.0, 300.0),
        (1e13, f64::NAN, 300.0),
        (f64::NAN, f64::NAN, 0.0),
        (f64::INFINITY, f64::NAN, 300.0),
        (0.0, -1.0, 1e-300),       // A·exp(+∞) = 0·∞ = NaN
        (-0.0, 50_000.0, 300.0),
        (0.0, 1e300, 1e-300),      // 0·exp(−∞) = 0·0
        (f64::MAX, -1.0, 1.0),     // overflows to +∞
        (-f64::MAX, 1.0, 1e300),   // −DBL_MAX·exp(−0) = −DBL_MAX
    ] {
        arr(&mut m, a, ea, t);
    }
    let prefactors = [0.0, -0.0, 1e-300, 1e-10, 1.0, 1e6, 1e13, 1e17, 1e300, f64::MAX, -1e13, f64::from_bits(1)];
    let energies = [
        -1e300, -1e7, -1e6, -1e5, -5e4, -1.0, -0.0, 0.0, 1.0, 1e3, 1e4, 5e4, 8e4, 1e5, 2e5, 5e5, 1e6, 1e7, 1e300,
    ];
    for a in prefactors {
        for ea in energies {
            for t in temperatures() {
                arr(&mut m, a, ea, t);
            }
        }
    }
    // the exponent on targets, at four temperatures and six prefactors
    for t in [1.0, 298.15, 1000.0, 5000.0] {
        for a in [1.0, 1e13, 1e300, 1e-300, f64::MAX, 0.5] {
            for x in TARGETS {
                arr(&mut m, a, -x * r * t, t);
            }
        }
    }
    // subnormal exp results scaled into the normal range: x between −745.2 and −708, A = 1e13 and up to DBL_MAX
    for _ in 0..400 {
        let x = rng.uniform(-745.2, -708.0);
        let t = rng.log_uniform(0.0, 3.7);
        arr(&mut m, 1e13, -x * r * t, t);
        let a = rng.log_uniform(0.0, 308.0);
        arr(&mut m, a, -x * r * t, t);
    }
    // realistic and wide sweeps
    for _ in 0..500 {
        let a = rng.log_uniform(5.0, 17.0);
        let ea = rng.uniform(0.0, 3e5);
        let t = rng.uniform(200.0, 3000.0);
        arr(&mut m, a, ea, t);
    }
    for _ in 0..300 {
        let a = rng.signed_log_uniform(-300.0, 300.0);
        let ea = rng.signed_log_uniform(-3.0, 8.0);
        let t = rng.log_uniform(-3.0, 5.0);
        arr(&mut m, a, ea, t);
    }

    // ── half_life_first_order ──
    let hl = |m: &mut Module, k: f64| {
        let out = res(half_life_first_order(bb(k)));
        m.row("kinetics::half_life_first_order", Class::B, &[f(k)], &[out]);
    };
    let mut ks = edges();
    ks.extend([0.1, -0.1, 1.0, 1.0_f64.next_up(), 1.0_f64.next_down(), 2.0, 0.5, 1e-3, 1e10]);
    ks.extend(decades(-323, 308, &["1", "3"]));
    for _ in 0..200 {
        ks.push(rng.log_uniform(-320.0, 308.0));
    }
    for k in ks {
        hl(&mut m, k);
    }

    // ── first_order_concentration ──
    let foc = |m: &mut Module, c0: f64, k: f64, t: f64| {
        let out = first_order_concentration(bb(c0), bb(k), bb(t));
        m.row("kinetics::first_order_concentration", Class::C, &[f(c0), f(k), f(t)], &[f(out)]);
    };
    for x in edges() {
        foc(&mut m, x, 0.1, 5.0);
        foc(&mut m, 1.0, x, 5.0);
        foc(&mut m, 1.0, 0.1, x);
    }
    // the Rust tests' values: at t½ and at t = 0
    for (c0, k, t) in [
        (1.0, 0.1, std::f64::consts::LN_2 / 0.1),
        (1.0, 0.1, 0.0),
        (1.0, 0.1, -0.0),
        (f64::INFINITY, 1e300, 1e300), // ∞·exp(−∞) = ∞·0 = NaN
        (0.0, -1e300, 1e300),          // 0·exp(+∞) = NaN
        (f64::MAX, -1.0, 1.0),         // overflows to +∞
        (1e300, 1.0, 745.0),           // a subnormal exp result scaled into the normal range
    ] {
        foc(&mut m, c0, k, t);
    }
    // k·t on targets, with [A]₀ over many decades (k = 1, t = −x; k = 0.25, t = −4x)
    for c0 in [1.0, 1e-300, 1e-10, 0.5, 1e13, 1e300, f64::MAX, -2.0, f64::from_bits(1)] {
        for x in TARGETS {
            foc(&mut m, c0, 1.0, -x);
            foc(&mut m, c0, 0.25, -4.0 * x);
        }
    }
    for c0 in decades(-300, 300, &["1"]).into_iter().step_by(10) {
        for kt in [0.0, 1e-3, 0.5, 1.0, 10.0, 100.0, 700.0, 720.0, 740.0, 745.2] {
            foc(&mut m, c0, kt, 1.0);
        }
    }
    // subnormal exp results scaled into the normal range: k·t between 708 and 745.2, [A]₀ up to 1e308
    for _ in 0..400 {
        let x = rng.uniform(708.0, 745.2);
        let k = rng.log_uniform(-3.0, 3.0);
        let c0 = rng.log_uniform(0.0, 308.0);
        foc(&mut m, c0, k, x / k);
        foc(&mut m, c0, x, 1.0);
    }
    for _ in 0..300 {
        let c0 = rng.log_uniform(-6.0, 1.0);
        let k = rng.log_uniform(-5.0, 2.0);
        let t = rng.log_uniform(-2.0, 4.0);
        foc(&mut m, c0, k, t);
    }
    for _ in 0..200 {
        let c0 = rng.signed_log_uniform(-300.0, 300.0);
        let k = rng.signed_log_uniform(-10.0, 10.0);
        let t = rng.log_uniform(-10.0, 10.0);
        foc(&mut m, c0, k, t);
    }

    // ── second_order_concentration ──
    let soc = audit::second;
    for x in edges() {
        soc(&mut m, x, 0.1, 10.0);
        soc(&mut m, 1.0, x, 10.0);
        soc(&mut m, 1.0, 0.1, x);
    }
    for (c0, k, t) in [
        (1.0, 0.1, 10.0),
        (0.0, 0.1, 10.0),
        (1.0, -0.1, 10.0),
        (1.0, -0.2, 10.0),
        (1.0, 0.1, -10.0),
        (1.0, f64::NAN, 10.0),
        (1.0, 0.1, f64::INFINITY),
        (f64::NAN, 0.1, 10.0),
        (f64::NAN, f64::NAN, f64::NAN), // the initial concentration is checked first
        (1.0, f64::NAN, -1.0),          // then k, then t
        (1e-300, 1e300, 1e300),         // k·t overflows: 1/∞ = 0
        (f64::MAX, 0.0, 0.0),
        (f64::from_bits(1), 0.0, 1.0),  // 1/[A]₀ overflows to ∞: 1/∞ = 0
    ] {
        soc(&mut m, c0, k, t);
    }
    for c0 in decades(-300, 300, &["1"]).into_iter().step_by(20) {
        for k in [0.0, 1e-10, 0.1, 1.0, 1e10] {
            for t in [0.0, 1.0, 1e3, 1e10] {
                soc(&mut m, c0, k, t);
            }
        }
    }
    for _ in 0..300 {
        let c0 = rng.log_uniform(-6.0, 1.0);
        let k = rng.log_uniform(-5.0, 5.0);
        let t = rng.log_uniform(-2.0, 5.0);
        soc(&mut m, c0, k, t);
    }

    // ── zero_order_concentration ──
    let zoc = |m: &mut Module, c0: f64, k: f64, t: f64| {
        let out = zero_order_concentration(bb(c0), bb(k), bb(t));
        m.row("kinetics::zero_order_concentration", Class::B, &[f(c0), f(k), f(t)], &[f(out)]);
    };
    for x in edges() {
        zoc(&mut m, x, 0.1, 5.0);
        zoc(&mut m, 1.0, x, 5.0);
        zoc(&mut m, 1.0, 0.1, x);
        zoc(&mut m, x, 0.0, 5.0);
        zoc(&mut m, x, -0.0, 5.0);
    }
    for (c0, k, t) in [
        (1.0, 0.1, 5.0),
        (1.0, 0.1, 20.0),
        (1.0, 0.1, 10.0),               // exactly complete: 1 − 1 = +0
        (f64::NAN, 0.1, 1.0),
        (1.0, f64::NAN, 1.0),
        (-0.0, 0.0, 1.0),               // −0 − 0 = −0, kept
        (-0.0, -0.0, 1.0),              // −0 − (−0) = +0
        (0.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),               // negative: floored to +0
        (f64::INFINITY, f64::INFINITY, 1.0), // ∞ − ∞ = NaN
        (1.0, f64::INFINITY, 0.0),           // ∞·0 = NaN
        (f64::MAX, -f64::MAX, 2.0),          // overflows to +∞
    ] {
        zoc(&mut m, c0, k, t);
    }
    for _ in 0..300 {
        let c0 = rng.log_uniform(-6.0, 2.0);
        let k = rng.log_uniform(-6.0, 1.0);
        let t = rng.log_uniform(-2.0, 4.0);
        zoc(&mut m, c0, k, t);
    }
    for _ in 0..100 {
        let c0 = rng.signed_log_uniform(-300.0, 300.0);
        let k = rng.signed_log_uniform(-300.0, 300.0);
        let t = rng.signed_log_uniform(-10.0, 10.0);
        zoc(&mut m, c0, k, t);
    }

    // ── zero_order_half_life ──
    let zhl = audit::zero_half;
    for x in edges() {
        zhl(&mut m, x, 0.1);
        zhl(&mut m, 1.0, x);
        zhl(&mut m, x, x);
    }
    for (c0, k) in [
        (1.0, 0.1),
        (0.0, 0.1),
        (-1.0, 0.1),
        (f64::NAN, 0.1),
        (1.0, f64::NAN),
        (f64::NAN, 0.0), // k is checked first
        (f64::MAX, f64::MIN_POSITIVE),
        (f64::MIN_POSITIVE, f64::MAX),
    ] {
        zhl(&mut m, c0, k);
    }
    for _ in 0..200 {
        let c0 = rng.log_uniform(-12.0, 2.0);
        let k = rng.log_uniform(-12.0, 3.0);
        zhl(&mut m, c0, k);
    }

    // ── nth_order_half_life ──
    let nth = |m: &mut Module, c0: f64, k: f64, n: f64| {
        let out = res(nth_order_half_life(bb(c0), bb(k), bb(n)));
        m.row("kinetics::nth_order_half_life", Class::C, &[f(c0), f(k), f(n)], &[out]);
    };
    for x in edges() {
        nth(&mut m, x, 0.1, 2.0);
        nth(&mut m, 1.0, x, 2.0);
        nth(&mut m, 1.0, 0.1, x);
    }
    // the Rust tests' values (and the order of the checks: the order, then k, then [A]₀)
    for (c0, k, n) in [
        (1.0, 0.1, 2.0),
        (1.0, 0.1, 3.0),
        (1.0, 0.1, 1.0),
        (1.0, 0.1, 1.0 + 1e-13),
        (1.0, 0.1, f64::NAN),
        (1.0, 0.1, f64::INFINITY),
        (1.0, 0.1, f64::NEG_INFINITY),
        (1.0, 0.1, 0.0),
        (1.0, 0.1, 0.5),
        (2.0, 0.1, 3.0),
        (1.0, 0.1, 1.0 + 1e-9),
        (1e-300, 0.1, 10.0),
        (f64::NAN, 0.1, 2.0),
        (1.0, f64::NAN, 2.0),
        (f64::NAN, f64::NAN, 1.0),
        (f64::NAN, -1.0, 2.0),
        (1e300, 0.1, 3.0),   // [A]₀^(n−1) overflows: t½ = 0, an Err
        (1e-300, 0.1, -1.0), // [A]₀^(n−1) overflows for n < 1 as well
    ] {
        nth(&mut m, c0, k, n);
    }
    // the order around 1: FIRST_ORDER_TOLERANCE (1e-12) and the ulps either side, from both sides
    let tol = 1e-12_f64;
    let hi = 1.0 + tol; // the doubles near 1 + 1e-12 are 2^-52 apart, those near 1 − 1e-12 2^-53
    let lo = 1.0 - tol;
    for base in [hi, lo] {
        for d in -3..=3 {
            let n = ulps(base, d);
            nth(&mut m, 1.0, 0.1, n);
            nth(&mut m, 0.5, 2.0, n);
        }
    }
    for d in [1e-15, 1e-13, 5e-13, 1.5e-12, 2e-12, 1e-11, 1e-10, 1e-9, 1e-7, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.5] {
        for n in [1.0 + d, 1.0 - d] {
            nth(&mut m, 1.0, 0.1, n);
            nth(&mut m, 1e-3, 5.0, n);
            nth(&mut m, 1e3, 1e-3, n);
        }
    }
    for k in 1..=8 {
        nth(&mut m, 1.0, 0.1, ulps(1.0, k));
        nth(&mut m, 1.0, 0.1, ulps(1.0, -k));
    }
    // integer and fractional orders out to ±300, [A]₀ and k over many decades
    let orders = [
        -300.0, -200.0, -150.5, -100.0, -50.0, -31.0, -10.0, -5.0, -2.0, -1.0, -0.5, 0.0, 0.25, 0.5, 1.5, 2.0, 2.5,
        3.0, 4.0, 5.0, 10.0, 20.0, 32.0, 50.0, 100.0, 105.0, 150.5, 200.0, 300.0,
    ];
    let inits = [
        1e-300, 1e-200, 1e-100, 1e-30, 1e-10, 1e-3, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3, 1e10, 1e30, 1e100, 1e200, 1e300,
    ];
    let rates = [1e-300, 1e-30, 1e-10, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e10, 1e20, 1e30, 1e100, 1e300];
    for n in orders {
        for c0 in inits {
            for k in rates {
                nth(&mut m, c0, k, n);
            }
        }
    }
    // [A]₀^(n−1) subnormal and the half-life still finite: the pow result divided into the normal range.
    // [A]₀ = 2^j·(1 + u/1024) and n − 1 = L/j put [A]₀^(n−1) near 2^L, L between −1074.5 and −1022.5 (n > 1 for
    // [A]₀ < 1, n < 1 for [A]₀ > 1); k from 10⁸ to 10³⁰⁰ keeps (n − 1)·k·[A]₀^(n−1) normal for most of them.
    for _ in 0..400 {
        let j = rng.int(9, 1000);
        let j = i32::try_from(if rng.int(0, 1) == 0 { -j } else { j }).unwrap_or(9);
        let l2 = rng.uniform(-1074.5, -1022.5);
        let n = 1.0 + l2 / f64::from(j);
        let c0 = (1.0 + rng.unit() * (1.0 / 1024.0)) * pow2(j);
        let k = rng.log_uniform(8.0, 300.0);
        nth(&mut m, c0, k, n);
    }
    // exact powers on the subnormal grid, at and next to its half-way points: 2^-1074, 2^-1075 (a tie, rounded to
    // even: 0, so Err), 2^-1073.5, 243·2^-1075 = 121.5·2^-1074 from (3·2^-215)^5 and from (9·2^-430)^2.5 (ties)
    for (c0, k, n) in [
        (2.0, 1e20, -1073.0),
        (2.0, 1e20, -1074.0),
        (2.0, 1e20, -1072.5),
        (0.5, 1e20, 1075.0),
        (3.0 * pow2(-215), 1e300, 6.0),
        (3.0 * pow2(-215), 1e300, 5.0),
        (9.0 * pow2(-430), 1e300, 3.5),
    ] {
        nth(&mut m, c0, k, n);
    }
    // realistic and wide sweeps
    for _ in 0..500 {
        let c0 = rng.log_uniform(-12.0, 2.0);
        let k = rng.log_uniform(-6.0, 6.0);
        let n = rng.uniform(-5.0, 5.0);
        nth(&mut m, c0, k, n);
    }
    for _ in 0..300 {
        let c0 = rng.log_uniform(-300.0, 300.0);
        let k = rng.log_uniform(-300.0, 300.0);
        let n = rng.uniform(-300.0, 300.0);
        nth(&mut m, c0, k, n);
    }

    // ── steady_state_intermediate, steady_state_rate ──
    let ssi = audit::steady;
    for x in edges() {
        ssi(&mut m, x, 1.0, 1.0);
        ssi(&mut m, 0.1, x, 1.0);
        ssi(&mut m, 0.1, 1.0, x);
    }
    for (k1, a, k2) in [
        (0.1, 1.0, 1.0),
        (0.1, 1.0, f64::NAN),
        (f64::NAN, 1.0, 1.0),
        (0.1, f64::NAN, 1.0),
        (f64::NAN, f64::NAN, 0.0), // k2 is checked first, then k1, then [A]
        (-1.0, f64::NAN, 1.0),
        (0.0, 0.0, 1.0),
        (-0.0, 1.0, 1.0),
        (1e300, 1e300, 1.0),       // overflows to +∞
        (1e-300, 1e-300, 1e300),   // underflows to 0
    ] {
        ssi(&mut m, k1, a, k2);
    }
    for _ in 0..200 {
        let k1 = rng.log_uniform(-6.0, 3.0);
        let a = rng.log_uniform(-6.0, 1.0);
        let k2 = rng.log_uniform(-3.0, 6.0);
        ssi(&mut m, k1, a, k2);
    }
    let ssr = |m: &mut Module, k: f64, c: f64| {
        let out = steady_state_rate(bb(k), bb(c));
        m.row("kinetics::steady_state_rate", Class::B, &[f(k), f(c)], &[f(out)]);
    };
    for x in edges() {
        ssr(&mut m, x, 2.0);
        ssr(&mut m, 0.05, x);
    }
    for (k, c) in [(0.05, 2.0), (f64::INFINITY, 0.0), (1e300, 1e300), (-0.0, 1.0)] {
        ssr(&mut m, k, c);
    }
    for _ in 0..150 {
        let k = rng.signed_log_uniform(-6.0, 6.0);
        let c = rng.signed_log_uniform(-6.0, 2.0);
        ssr(&mut m, k, c);
    }

    // ── pre_equilibrium_rate ──
    let pre = audit::pre;
    for x in edges() {
        pre(&mut m, x, 2.0, 0.1, 1.0);
        pre(&mut m, 10.0, x, 0.1, 1.0);
        pre(&mut m, 10.0, 2.0, x, 1.0);
        pre(&mut m, 10.0, 2.0, 0.1, x);
    }
    for (kf, kr, ks, a) in [
        (10.0, 2.0, 0.1, 1.0),
        (1.0, 0.0, 0.1, 1.0),
        (1.0, f64::NAN, 0.1, 1.0),
        (f64::NAN, 2.0, 0.1, 1.0),
        (1.0, 2.0, f64::NAN, 1.0),
        (1.0, 2.0, 0.1, f64::NAN),
        (f64::NAN, 0.0, f64::NAN, f64::NAN), // k_reverse first, then k_forward / k_slow, then [A]
        (1.0, 2.0, -1.0, f64::NAN),
        (1e300, 1e-300, 1.0, 1.0),           // K_eq overflows to +∞
        (1e300, 1.0, 1e300, 0.0),            // ∞·0 = NaN
        (1e-300, 1e300, 1e-10, 1.0),         // underflows to 0
    ] {
        pre(&mut m, kf, kr, ks, a);
    }
    for _ in 0..200 {
        let kf = rng.log_uniform(-3.0, 6.0);
        let kr = rng.log_uniform(-3.0, 6.0);
        let ks = rng.log_uniform(-6.0, 3.0);
        let a = rng.log_uniform(-6.0, 1.0);
        pre(&mut m, kf, kr, ks, a);
    }

    // ── michaelis_menten, lineweaver_burk ──
    let mm = audit::mm;
    for x in edges() {
        mm(&mut m, x, 5.0, 10.0);
        mm(&mut m, 100.0, x, 10.0);
        mm(&mut m, 100.0, 5.0, x);
    }
    for (v, km, s) in [
        (100.0, 5.0, 5.0),
        (100.0, 5.0, 1e6),
        (100.0, 5.0, 0.0),
        (100.0, 5.0, 10.0),
        (0.0, 5.0, 1.0),
        (100.0, 0.0, 1.0),
        (f64::NAN, 5.0, 1.0),
        (100.0, f64::NAN, 1.0),
        (100.0, 5.0, f64::NAN),
        (100.0, 5.0, f64::INFINITY),
        (f64::NAN, f64::NAN, f64::NAN), // Vmax is checked first, then Km, then [S]
        (1.0, f64::NAN, -1.0),
        (f64::MAX, 1.0, f64::MAX),      // Vmax·[S] overflows: ∞/∞ = NaN
        (1e300, 1e-300, 1e300),         // ∞/1e300 = ∞
        (1e-300, 1e300, 1e-300),        // underflows to 0
        (100.0, f64::MAX, f64::MAX),    // Km + [S] overflows: (Vmax·[S])/∞ = 0
    ] {
        mm(&mut m, v, km, s);
    }
    for _ in 0..300 {
        let v = rng.log_uniform(-6.0, 3.0);
        let km = rng.log_uniform(-6.0, 1.0);
        let s = rng.log_uniform(-8.0, 3.0);
        mm(&mut m, v, km, s);
    }
    let lb = |m: &mut Module, v: f64, km: f64| {
        let out = pair(lineweaver_burk(bb(v), bb(km)));
        m.row("kinetics::lineweaver_burk", Class::B, &[f(v), f(km)], &out);
    };
    for x in edges() {
        lb(&mut m, x, 5.0);
        lb(&mut m, 100.0, x);
    }
    for (v, km) in [
        (100.0, 5.0),
        (f64::NAN, 5.0),
        (100.0, f64::INFINITY),
        (f64::NAN, f64::NAN),
        (0.0, f64::NAN), // Vmax is checked first
        (f64::MIN_POSITIVE, f64::MAX),       // Km/Vmax overflows to +∞
        (f64::from_bits(1), 1.0),            // 1/Vmax overflows to +∞
        (f64::MAX, f64::MIN_POSITIVE),       // both underflow
    ] {
        lb(&mut m, v, km);
    }
    for _ in 0..200 {
        let v = rng.log_uniform(-6.0, 3.0);
        let km = rng.log_uniform(-6.0, 1.0);
        lb(&mut m, v, km);
    }

    // ── collision_theory_rate ──
    let col = |m: &mut Module, z: f64, p: f64, ea: f64, t: f64| {
        let out = res(collision_theory_rate(bb(z), bb(p), bb(ea), bb(t)));
        m.row("kinetics::collision_theory_rate", Class::C, &[f(z), f(p), f(ea), f(t)], &[out]);
    };
    for x in edges() {
        col(&mut m, x, 0.1, 50_000.0, 300.0);
        col(&mut m, 1e10, x, 50_000.0, 300.0);
        col(&mut m, 1e10, 0.1, x, 300.0);
        col(&mut m, 1e10, 0.1, 50_000.0, x);
    }
    let sterics = [
        0.0,
        -0.0,
        f64::from_bits(1),
        1e-300,
        1e-6,
        0.01,
        0.1,
        0.5,
        1.0_f64.next_down(),
        1.0,
        1.0_f64.next_up(),
        1.5,
        -0.1,
        f64::NAN,
        f64::INFINITY,
    ];
    for p in sterics {
        col(&mut m, 1e10, p, 50_000.0, 300.0);
        col(&mut m, f64::NAN, p, f64::NAN, 300.0); // the steric factor is checked before Z and Ea
    }
    for (z, p, ea, t) in [
        (1e10, 0.1, 50_000.0, 300.0),
        (1e10, 1.0, 50_000.0, 300.0),
        (1e10, 0.01, 50_000.0, 300.0),
        (1e10, 0.0, 50_000.0, 300.0),
        (1e10, 1.5, 50_000.0, 300.0),
        (-1e10, 0.1, 50_000.0, 300.0),
        (0.0, 0.1, 50_000.0, 300.0),
        (f64::NAN, 0.1, 50_000.0, 300.0),
        (1e10, f64::NAN, 50_000.0, 300.0),
        (1e10, 0.1, f64::NAN, 300.0),
        (1e10, 0.1, 50_000.0, f64::NAN),
        (f64::NAN, f64::NAN, f64::NAN, 0.0), // the temperature is checked first
        (f64::NAN, 0.5, f64::NAN, 300.0),    // then the steric factor, then Z, then Ea
        (1e10, 0.5, f64::NAN, 300.0),
        (f64::MAX, 1.0, -1.0, 1.0),          // overflows to +∞
    ] {
        col(&mut m, z, p, ea, t);
    }
    for t in [1.0, 298.15, 1000.0, 5000.0] {
        for (z, p) in [(1e10, 0.1), (1e300, 1.0), (f64::MAX, 0.5), (1e-300, 1e-6), (1.0, 1.0)] {
            for x in TARGETS {
                col(&mut m, z, p, -x * r * t, t);
            }
        }
    }
    for z in [1e-300, 1.0, 1e10, 1e300] {
        for ea in [-1e6, -1e3, 0.0, 1e3, 5e4, 1e5, 1e6, 1e300] {
            for t in temperatures() {
                col(&mut m, z, 0.1, ea, t);
            }
        }
    }
    // subnormal exp results scaled into the normal range by Z·p
    for _ in 0..300 {
        let x = rng.uniform(-745.2, -708.0);
        let t = rng.log_uniform(0.0, 3.7);
        let z = rng.log_uniform(0.0, 308.0);
        let p = rng.uniform(0.05, 1.0);
        col(&mut m, z, p, -x * r * t, t);
    }
    for _ in 0..300 {
        let z = rng.log_uniform(8.0, 12.0);
        let p = rng.log_uniform(-6.0, 0.0);
        let ea = rng.uniform(0.0, 3e5);
        let t = rng.uniform(200.0, 3000.0);
        col(&mut m, z, p, ea, t);
    }

    // ── eyring_rate ──
    let ey = |m: &mut Module, dg: f64, t: f64| {
        let out = res(eyring_rate(bb(dg), bb(t)));
        m.row("kinetics::eyring_rate", Class::C, &[f(dg), f(t)], &[out]);
    };
    for x in edges() {
        ey(&mut m, x, 298.15);
        ey(&mut m, 80_000.0, x);
    }
    for (dg, t) in [
        (80_000.0, 298.15),
        (0.0, 298.15),
        (80_000.0, 298.0),
        (80_000.0, 400.0),
        (80_000.0, 0.0),
        (80_000.0, f64::NAN),
        (f64::NAN, 298.15),
        (f64::NAN, 0.0), // the temperature is checked first
        (80_000.0 - 298.0 * -50.0, 298.0),
        (-1.0, 1e-300), // exp(+∞)·prefactor: ∞
        (1.0, f64::MAX), // the prefactor overflows: ∞·exp(−0) = ∞
    ] {
        ey(&mut m, dg, t);
    }
    for dg in [-1e7, -1e6, -1e5, -1e4, -1.0, -0.0, 0.0, 1.0, 1e3, 1e4, 5e4, 8e4, 1e5, 2e5, 1e6, 1e7, 1e300, -1e300] {
        for t in temperatures() {
            ey(&mut m, dg, t);
        }
    }
    for t in [1.0, 2.0, 10.0, 298.15, 1000.0, 5000.0] {
        for x in TARGETS {
            ey(&mut m, -x * r * t, t);
        }
    }
    // at low T the prefactor is ~2e10·T, and ΔG‡ of a few kJ/mol already puts exp in the subnormal range
    for _ in 0..300 {
        let t = rng.uniform(0.5, 30.0);
        let x = rng.uniform(-745.2, -708.0);
        ey(&mut m, -x * r * t, t);
    }
    for _ in 0..300 {
        let dg = rng.uniform(0.0, 2e5);
        let t = rng.uniform(200.0, 3000.0);
        ey(&mut m, dg, t);
    }
    for _ in 0..200 {
        let dg = rng.signed_log_uniform(-3.0, 8.0);
        let t = rng.log_uniform(-3.0, 5.0);
        ey(&mut m, dg, t);
    }

    // ── eyring_rate_from_activation ──
    let eya = |m: &mut Module, dh: f64, ds: f64, t: f64| {
        let out = res(eyring_rate_from_activation(bb(dh), bb(ds), bb(t)));
        m.row("kinetics::eyring_rate_from_activation", Class::C, &[f(dh), f(ds), f(t)], &[out]);
    };
    for x in edges() {
        eya(&mut m, x, -50.0, 298.15);
        eya(&mut m, 80_000.0, x, 298.15);
        eya(&mut m, 80_000.0, -50.0, x);
    }
    for (dh, ds, t) in [
        (80_000.0, -50.0, 298.0),
        (80_000.0, -50.0, 298.15),
        (f64::NAN, -50.0, 298.15),
        (80_000.0, f64::NAN, 298.15),
        (80_000.0, -50.0, f64::INFINITY),
        (f64::NAN, f64::NAN, 0.0),     // the temperature is checked first, then ΔH‡, then ΔS‡
        (f64::NAN, f64::NAN, 298.15),
        (0.0, 0.0, 298.15),
        (1e7, 6000.0, 298.15),         // exp(ΔS‡/R) overflows, exp(−ΔH‡/RT) underflows: ∞·0 = NaN
        (-1e7, -6000.0, 298.15),       // 0·∞ = NaN
        (-3e6, -6000.0, 298.15),       // exp(−ΔH‡/RT) overflows: (prefactor·subnormal)·∞ = ∞
        (-1.7e6, -6000.0, 298.15),     // a subnormal exp(ΔS‡/R) times a large exp(−ΔH‡/RT): a normal k
    ] {
        eya(&mut m, dh, ds, t);
    }
    let dss = [-1e300, -7000.0, -6000.0, -5000.0, -500.0, -50.0, -0.0, 0.0, 50.0, 500.0, 5000.0, 6000.0, 1e300];
    let dhs = [-1e7, -1e6, -1e5, -0.0, 0.0, 1e4, 8e4, 1e5, 1e6, 1e7];
    for ds in dss {
        for dh in dhs {
            for t in [1.0, 77.0, 298.15, 1000.0, 5000.0, 1e300] {
                eya(&mut m, dh, ds, t);
            }
        }
    }
    // either exp result subnormal and the other one scaling it back: x_S = ΔS‡/R and x_H = −ΔH‡/RT on targets
    for t in [1.0, 298.15, 3000.0] {
        for xs in [-745.0, -730.0, -720.25, -710.0, -700.0, -1.0, 0.0, 1.0, 300.0, 700.0, 709.0] {
            for xh in [-745.0, -730.0, -710.0, -1.0, 0.0, 1.0, 300.0, 700.0, 709.0, 710.0] {
                eya(&mut m, -xh * r * t, xs * r, t);
            }
        }
    }
    for _ in 0..300 {
        let t = rng.log_uniform(0.0, 3.7);
        let xs = rng.uniform(-745.2, -708.0);
        let xh = rng.uniform(0.0, 700.0);
        eya(&mut m, -xh * r * t, xs * r, t);
        eya(&mut m, -xs * r * t, xh * r, t);
    }
    for _ in 0..300 {
        let dh = rng.uniform(0.0, 2e5);
        let ds = rng.uniform(-300.0, 100.0);
        let t = rng.uniform(200.0, 3000.0);
        eya(&mut m, dh, ds, t);
    }
    audit::emit_kinetics(&mut m);
    m
}
