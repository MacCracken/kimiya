//! `thermo.tsv` — calorimetry and Hess's law (port plan §4, M2). Every function is + − × ÷ only,
//! so every row is class B (expected 0 ulp).
//!
//! - `thermo::heat_transfer` (m·c·ΔT, no validation): each argument over the edges with the others
//!   typical, the Rust tests' values, overflow and NaN propagation, and a seeded sweep.
//! - `thermo::hess_law` (Σ of a slice): the input is `v:n` then n `f:` values. Rust's f64 `sum` folds
//!   from −0.0, so the empty slice, `[0.0]`, `[-0.0]` and the zero mixes pin the sign of a zero sum;
//!   then ∞ − ∞, NaN, overflow, cancellation (1e16 + 1 − 1e16) and seeded slices of 1–12 values.
//! - `thermo::enthalpy_from_formation` (Σ products − Σ reactants): two slices, each `v:n` then its
//!   values: the Rust test, empty sides, and seeded pairs.
//! - `thermo::specific_heat_from_calorimetry` (q / (m·ΔT), `Result`): the edges per argument,
//!   |ΔT| at the hisab `EPSILON_F64` = 10⁻¹² threshold and one ulp either side, both signs, the Rust
//!   tests' values and a seeded sweep.

use std::hint::black_box as bb;

use crate::grids::Rng;
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, f};
use kimiya::thermo::{
    WATER_SPECIFIC_HEAT, enthalpy_from_formation, heat_transfer, hess_law,
    specific_heat_from_calorimetry,
};

fn seq(xs: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(xs.len())];
    v.extend(xs.iter().map(|&x| f(x)));
    v
}

fn slice(rng: &mut Rng, lo: i64, hi: i64) -> Vec<f64> {
    let n = rng.int(lo, hi);
    (0..n).map(|_| rng.signed_log_uniform(-3.0, 7.0)).collect()
}

pub fn emit() -> Module {
    let mut m = Module::new("thermo");
    let heat = |m: &mut Module, ma: f64, c: f64, dt: f64| {
        let q = heat_transfer(bb(ma), bb(c), bb(dt));
        m.row("thermo::heat_transfer", Class::B, &[f(ma), f(c), f(dt)], &[f(q)]);
    };
    for x in edges() {
        heat(&mut m, x, WATER_SPECIFIC_HEAT, 10.0);
        heat(&mut m, 100.0, x, 10.0);
        heat(&mut m, 100.0, WATER_SPECIFIC_HEAT, x);
    }
    for (ma, c, dt) in [
        (1000.0, WATER_SPECIFIC_HEAT, 1.0),
        (100.0, WATER_SPECIFIC_HEAT, 10.0),
        (100.0, WATER_SPECIFIC_HEAT, -5.0),
        (1e300, 1e300, 1.0),
        (f64::INFINITY, 0.0, 1.0),
    ] {
        heat(&mut m, ma, c, dt);
    }
    let mut rng = Rng::new(0x2f01);
    for _ in 0..300 {
        let ma = rng.log_uniform(-3.0, 6.0);
        let c = rng.log_uniform(-1.0, 1.5);
        let dt = rng.signed_log_uniform(-3.0, 3.0);
        heat(&mut m, ma, c, dt);
    }

    let mut slices: Vec<Vec<f64>> = vec![
        vec![],
        vec![0.0],
        vec![-0.0],
        vec![-0.0, -0.0],
        vec![0.0, -0.0],
        vec![-0.0, 0.0],
        vec![-100.0, 50.0, -30.0],
        vec![1e308, 1e308],
        vec![f64::INFINITY, f64::NEG_INFINITY],
        vec![f64::INFINITY, 1.0],
        vec![f64::NAN],
        vec![1.0, f64::NAN, 2.0],
        vec![1e16, 1.0, -1e16],
        vec![1.0, 1e16, -1e16],
        vec![f64::MIN_POSITIVE, f64::from_bits(1)],
    ];
    for _ in 0..200 {
        slices.push(slice(&mut rng, 1, 12));
    }
    for xs in &slices {
        let s = hess_law(bb(xs));
        m.row("thermo::hess_law", Class::B, &seq(xs), &[f(s)]);
    }

    let mut pairs: Vec<(Vec<f64>, Vec<f64>)> = vec![
        (vec![-400.0], vec![-200.0, -100.0]),
        (vec![], vec![]),
        (vec![], vec![-0.0]),
        (vec![-0.0], vec![]),
        (vec![-0.0], vec![0.0]),
        (vec![5.0], vec![]),
        (vec![], vec![5.0]),
        (vec![f64::INFINITY], vec![f64::INFINITY]),
        (vec![f64::NAN], vec![1.0]),
    ];
    for _ in 0..150 {
        let p = slice(&mut rng, 0, 6);
        pairs.push((p, slice(&mut rng, 0, 6)));
    }
    for (p, r) in &pairs {
        let dh = enthalpy_from_formation(bb(p), bb(r));
        let mut ins = seq(p);
        ins.extend(seq(r));
        m.row("thermo::enthalpy_from_formation", Class::B, &ins, &[f(dh)]);
    }

    let cal = |m: &mut Module, e: f64, ma: f64, dt: f64| {
        let out = res(specific_heat_from_calorimetry(bb(e), bb(ma), bb(dt)));
        m.row("thermo::specific_heat_from_calorimetry", Class::B, &[f(e), f(ma), f(dt)], &[out]);
    };
    for x in edges() {
        cal(&mut m, 4184.0, x, 1.0);
        cal(&mut m, x, 1000.0, 1.0);
        cal(&mut m, 4184.0, 1000.0, x);
    }
    // |ΔT| against hisab's EPSILON_F64 (1e-12): the threshold and one ulp either side.
    let eps = 1e-12_f64;
    for dt in [eps.next_down(), eps, eps.next_up(), 2.0 * eps, 0.5 * eps] {
        cal(&mut m, 4184.0, 1000.0, dt);
        cal(&mut m, 4184.0, 1000.0, -dt);
    }
    for (e, ma, dt) in [
        (4184.0, 1000.0, 1.0),
        (100.0, 0.0, 1.0),
        (100.0, 100.0, 0.0),
        (100.0, f64::NAN, 1.0),
        (100.0, 100.0, f64::NAN),
        (f64::NAN, 100.0, 1.0),
        (100.0, f64::INFINITY, 1.0),
        (100.0, 100.0, f64::INFINITY),
        (1e300, 1e-300, 1e-12), // q / (m·ΔT) overflows to +∞ on valid inputs
        (-0.0, 100.0, 1.0),
    ] {
        cal(&mut m, e, ma, dt);
    }
    for _ in 0..300 {
        let e = rng.signed_log_uniform(-3.0, 7.0);
        let ma = rng.log_uniform(-3.0, 6.0);
        let dt = rng.signed_log_uniform(-3.0, 3.0);
        cal(&mut m, e, ma, dt);
    }
    m
}
