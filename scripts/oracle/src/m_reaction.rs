//! `reaction.tsv` — Gibbs free energy, the equilibrium constant and spontaneity (port plan §4, M2).
//!
//! - `reaction::gibbs_free_energy` (class B: ΔH − T·ΔS, one multiply and one subtract, expected
//!   0 ulp): each argument over the validation edges with the others typical, then the Rust tests'
//!   values and a seeded sweep (ΔH ±10⁻²…10⁷ J/mol, T 10⁻³…10⁵ K, ΔS ±10⁻²…10⁴ J/(mol·K)).
//! - `reaction::equilibrium_constant` (class C: one `exp`): the edges, then ΔG built from a target
//!   x = −ΔG/RT at three temperatures, across the underflow (x ≈ −745) and overflow (x ≈ 709.78)
//!   ends and zero, the Rust tests' values, and a seeded sweep of x over −760…720; then (ADR 0006
//!   item 6, M3c-3) the subnormal results densely, on targets at five temperatures and seeded.
//! - `reaction::is_spontaneous` (class A, `i:0/1`): the special values and both signs.
//!
//! Every input goes through `black_box`, so each row is a run-time call (the transcendental bits
//! are glibc's, never a compile-time fold).

use std::hint::black_box as bb;

use crate::grids::{EXP_SUBNORMAL_X, Rng, specials};
use crate::tsv::{Class, Module, Tok, b, err, f};
use kimiya::element::GAS_CONSTANT;
use kimiya::error::KimiyaError;
use kimiya::reaction::{equilibrium_constant, gibbs_free_energy, is_spontaneous};

pub fn res(r: Result<f64, KimiyaError>) -> Tok {
    r.map_or_else(|e| err(&e), f)
}

/// The special values plus the magnitudes the validation guards meet: negative and tiny
/// positive, the extremes and a half.
pub fn edges() -> Vec<f64> {
    let mut v = specials();
    v.extend([-0.5, -300.0, 0.5, 2.0, 1e-300, -1e-300, 1e300, -1e300]);
    v
}

pub fn emit() -> Module {
    let mut m = Module::new("reaction");
    let gibbs = |m: &mut Module, h: f64, t: f64, s: f64| {
        let out = res(gibbs_free_energy(bb(h), bb(t), bb(s)));
        m.row("reaction::gibbs_free_energy", Class::B, &[f(h), f(t), f(s)], &[out]);
    };
    for x in edges() {
        gibbs(&mut m, -100_000.0, x, -200.0);
        gibbs(&mut m, x, 298.0, -200.0);
        gibbs(&mut m, -100_000.0, 298.0, x);
    }
    for (h, t, s) in [
        (-100_000.0, 298.0, -200.0),
        (-100_000.0, 0.0, -200.0),
        (1e308, 1e308, -1e308), // ΔH − TΔS overflows to +∞ on finite inputs
        (-1e308, 1e308, 1e308),
        (0.0, 298.15, 0.0),
        (-0.0, 298.15, 0.0),
    ] {
        gibbs(&mut m, h, t, s);
    }
    let mut rng = Rng::new(0x2e01);
    for _ in 0..400 {
        let h = rng.signed_log_uniform(-2.0, 7.0);
        let t = rng.log_uniform(-3.0, 5.0);
        let s = rng.signed_log_uniform(-2.0, 4.0);
        gibbs(&mut m, h, t, s);
    }

    let keq = |m: &mut Module, dg: f64, t: f64| {
        let out = res(equilibrium_constant(bb(dg), bb(t)));
        m.row("reaction::equilibrium_constant", Class::C, &[f(dg), f(t)], &[out]);
    };
    for x in edges() {
        keq(&mut m, -5000.0, x);
        keq(&mut m, x, 298.0);
    }
    // ΔG = −x·R·T for a target exponent x = −ΔG/RT, across both ends of exp's range.
    let xs = [
        -800.0, -746.0, -745.2, -745.1, -745.0, -744.0, -709.0, -708.4, -700.0, -1.0, -1e-10,
        0.0, 1e-10, 1.0, 700.0, 708.4, 709.0, 709.7, 709.78, 709.79, 709.8, 710.0, 800.0,
    ];
    for t in [1.0, 298.0, 1000.0] {
        for &x in &xs {
            keq(&mut m, -x * GAS_CONSTANT * t, t);
        }
    }
    for (dg, t) in [
        (-5000.0, 298.0),
        (5000.0, 298.0),
        (0.0, 298.0),
        (-0.0, 298.0),
        (0.0, 0.0),
        (-2.0e6, 298.0),
        (2.0e6, 298.0),
        (-700.0 * GAS_CONSTANT * 298.0, 298.0),
    ] {
        keq(&mut m, dg, t);
    }
    let mut rng = Rng::new(0x2e02);
    for _ in 0..500 {
        let t = rng.log_uniform(0.0, 4.0);
        let x = rng.uniform(-760.0, 720.0);
        keq(&mut m, -x * GAS_CONSTANT * t, t);
    }
    for _ in 0..100 {
        let t = rng.log_uniform(-3.0, 5.0);
        let dg = rng.signed_log_uniform(-3.0, 7.0);
        keq(&mut m, dg, t);
    }
    // ADR 0006's hard rows (M3c-3): K is the exp result itself, so the subnormal end of exp's range, densely:
    // x = −ΔG/RT on the subnormal targets at five temperatures (not the three above, so no row repeats), then
    // seeded.
    for t in [2.0, 298.15, 1500.0, 5000.0, 1e-3] {
        for x in EXP_SUBNORMAL_X {
            keq(&mut m, -x * GAS_CONSTANT * t, t);
        }
    }
    let mut hard = Rng::new(0x2e03);
    for _ in 0..300 {
        let t = hard.log_uniform(-1.0, 4.0);
        let x = hard.uniform(-745.2, -708.0);
        keq(&mut m, -x * GAS_CONSTANT * t, t);
    }

    let mut probes = specials();
    probes.extend([-100.0, 100.0, -1e-300, 1e-300, -5e-324]);
    for x in probes {
        m.row("reaction::is_spontaneous", Class::A, &[f(x)], &[b(is_spontaneous(bb(x)))]);
    }
    m
}
