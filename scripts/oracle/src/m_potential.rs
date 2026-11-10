//! `potential.tsv` — intermolecular potentials (Lennard-Jones, Morse, Coulomb) and molecular geometry from 3D
//! coordinates (distance, bond and dihedral angles) (port plan §4, M4d).
//!
//! Rust's `[f64; 3]` points are flattened into their three `f:` tokens, in argument order (`distance_3d` has 6 input
//! floats, `bond_angle` 9, `dihedral_angle` 12).
//! - Class B (+ − × ÷ and √, expected 0 ulp): `lennard_jones` (each argument over the edges, the Rust tests' values,
//!   r/σ from 0.5 to 100 at Å and metre scale, (σ/r)¹² and (σ/r)⁶ over- and underflowing, seeded; **ADR 0002 #17**: ε ≥
//!   2¹⁰²², where Rust's 4·ε overflows, the rows noted with the port's (ε·d)·4), `morse_width_parameter` (the edges, the
//!   Rust tests' HCl values, every decade, seeded; **ADR 0002 #18**: D_e ≥ 2¹⁰²³, where Rust's 2·D_e overflows, noted
//!   with the port's ω·√((μ·½)/D_e)), `coulomb`, `coulomb_force` (the edges, the Rust tests' charges, every decade of
//!   r, signed charges, the intermediates' over- and underflow, seeded), `distance_3d` (each coordinate over the edges,
//!   the Rust tests' points, Å and metre scale, over- and underflowing squares, seeded).
//! - Class C (exp, atan2, exp2; bound 0 against the corrected oracle): `lj_equilibrium_distance` (`σ·2.0.powf(1/6)`: σ
//!   over every binade, subnormals, the overflow end; **ADR 0002 #16**: the rows where the port's correctly rounded
//!   σ·2^(1/6) differs from Rust's σ·RN(2^(1/6)), noted), `morse` (each argument over the edges, the Rust tests' values,
//!   r = r_e exactly, a(r − r_e) on targets across exp's range, **near-cancellation rows** where 1 − e^(−a(r − r_e))
//!   cancels (|a(r − r_e)| from 10⁻¹⁷ to 1), seeded), `bond_angle` and `dihedral_angle` (atan2: the Rust tests'
//!   geometries, right / linear / tiny angles, trans / cis / gauche / perpendicular dihedrals and their mirror images
//!   (the −0.0 y-term that Rust maps from −π to π), every coordinate slot over the edges, Å, metre and femtometre
//!   scale, near-coincident atoms **on Rust's 10⁻¹² tolerance to the ulp** and near-collinear atoms **on its 10⁻²⁰
//!   tolerance to the ulp**, every Err path, ±0, NaN, ±∞, DBL_MAX and subnormal coordinates, geometries scaled until
//!   their squares over- or underflow, seeded).
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::exact::cr_root6_2_mul;
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, f};
use kimiya::potential::{
    bond_angle, coulomb, coulomb_force, dihedral_angle, distance_3d, lennard_jones, lj_equilibrium_distance, morse,
    morse_width_parameter,
};

type P3 = [f64; 3];

fn toks(ps: &[P3]) -> Vec<Tok> {
    ps.iter().flat_map(|p| p.iter().map(|&x| f(x))).collect()
}

fn scale(p: P3, s: f64) -> P3 {
    [p[0] * s, p[1] * s, p[2] * s]
}

/// ADR 0002 #17: the port's ε·d·4 where Rust's 4·ε overflows (ε ≥ 2¹⁰²²), d = (σ/r)¹² − (σ/r)⁶ in Rust's operation order;
/// `None` elsewhere (Rust's `4.0 * epsilon * (sr12 - sr6)` stays).
fn lj_exact(eps: f64, sigma: f64, r: f64) -> Option<f64> {
    let ok = |x: f64| x > 0.0 && x.is_finite();
    if !(ok(eps) && ok(sigma) && ok(r)) || eps < pow2(1022) {
        return None;
    }
    let sr = sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    Some(eps * (sr12 - sr6) * 4.0)
}

/// ADR 0002 #18: the port's ω·√((μ·½)/D_e) where Rust's 2·D_e overflows (D_e ≥ 2¹⁰²³): (μ·½)/D_e is the correctly rounded
/// μ/(2·D_e) there (`validate.py` proves it); `None` elsewhere.
fn width_exact(w: f64, mu: f64, d: f64) -> Option<f64> {
    let ok = |x: f64| x > 0.0 && x.is_finite();
    if !(ok(w) && ok(mu) && ok(d)) || d < pow2(1023) {
        return None;
    }
    Some(w * ((mu * 0.5) / d).sqrt())
}

fn lj(m: &mut Module, eps: f64, sigma: f64, r: f64) {
    let out = lennard_jones(bb(eps), bb(sigma), bb(r));
    let name = "potential::lennard_jones";
    let ins = [f(eps), f(sigma), f(r)];
    match (out.as_ref().ok(), lj_exact(eps, sigma, r)) {
        (Some(&v), Some(x)) if v.to_bits() != x.to_bits() && !(v.is_nan() && x.is_nan()) => {
            let note = format!(
                "ADR 0002 #17: the port returns (ε·d)·4, the correctly rounded 4·ε·d for Rust's d = (σ/r)¹² − (σ/r)⁶ = \
                 {x:?} (f:{:016x})",
                x.to_bits()
            );
            m.row_noted(name, Class::B, &ins, &[f(v)], &note);
        }
        _ => m.row(name, Class::B, &ins, &[res(out)]),
    }
}

fn rmin(m: &mut Module, sigma: f64) {
    let v = lj_equilibrium_distance(bb(sigma));
    let name = "potential::lj_equilibrium_distance";
    if sigma.is_finite() && sigma != 0.0 {
        let x = cr_root6_2_mul(sigma, v);
        if x.to_bits() != v.to_bits() {
            let note =
                format!("ADR 0002 #16: the port returns the correctly rounded σ·2^(1/6) = {x:?} (f:{:016x})", x.to_bits());
            m.row_noted(name, Class::C, &[f(sigma)], &[f(v)], &note);
            return;
        }
    }
    m.row(name, Class::C, &[f(sigma)], &[f(v)]);
}

fn morse_row(m: &mut Module, d: f64, a: f64, re: f64, r: f64) {
    let out = res(morse(bb(d), bb(a), bb(re), bb(r)));
    m.row("potential::morse", Class::C, &[f(d), f(a), f(re), f(r)], &[out]);
}

fn width(m: &mut Module, w: f64, mu: f64, d: f64) {
    let out = morse_width_parameter(bb(w), bb(mu), bb(d));
    let name = "potential::morse_width_parameter";
    let ins = [f(w), f(mu), f(d)];
    match (out.as_ref().ok(), width_exact(w, mu, d)) {
        (Some(&v), Some(x)) if v.to_bits() != x.to_bits() => {
            let note = format!(
                "ADR 0002 #18: the port returns ω·√((μ·½)/D_e), with the correctly rounded μ/(2·D_e) = {x:?} (f:{:016x})",
                x.to_bits()
            );
            m.row_noted(name, Class::B, &ins, &[f(v)], &note);
        }
        _ => m.row(name, Class::B, &ins, &[res(out)]),
    }
}

fn coul(m: &mut Module, q1: f64, q2: f64, r: f64) {
    let v = res(coulomb(bb(q1), bb(q2), bb(r)));
    m.row("potential::coulomb", Class::B, &[f(q1), f(q2), f(r)], &[v]);
    let fv = res(coulomb_force(bb(q1), bb(q2), bb(r)));
    m.row("potential::coulomb_force", Class::B, &[f(q1), f(q2), f(r)], &[fv]);
}

fn dist(m: &mut Module, p1: P3, p2: P3) {
    let d = distance_3d(bb(p1), bb(p2));
    m.row("potential::distance_3d", Class::B, &toks(&[p1, p2]), &[f(d)]);
}

fn angle(m: &mut Module, a: P3, b: P3, c: P3) {
    let out = res(bond_angle(bb(a), bb(b), bb(c)));
    m.row("potential::bond_angle", Class::C, &toks(&[a, b, c]), &[out]);
}

fn dihedral(m: &mut Module, a: P3, b: P3, c: P3, d: P3) {
    let out = res(dihedral_angle(bb(a), bb(b), bb(c), bb(d)));
    m.row("potential::dihedral_angle", Class::C, &toks(&[a, b, c, d]), &[out]);
}

/// A point with coordinates uniform in [−s, s].
fn rpt(rng: &mut Rng, s: f64) -> P3 {
    [rng.uniform(-s, s), rng.uniform(-s, s), rng.uniform(-s, s)]
}

fn add(a: P3, b: P3) -> P3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Rust's bond-angle coincidence test for A against B (`ba_sq <= tol² · max(|a|², |b|²)`, the `is_positive` checks
/// aside), with the tolerance `tol`: the replica that finds the rows on the threshold.
fn coincident(a: P3, b: P3, tol: f64) -> bool {
    let dot = |x: P3, y: P3| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let ba_sq = dot(ba, ba);
    let b_sq = dot(b, b);
    ba_sq <= tol * tol * dot(a, a).max(b_sq)
}

/// Rust's dihedral collinearity test for A-B-C (`n1_sq <= tol · |b1|² · |b2|²`), with the tolerance `tol`.
fn collinear(a: P3, b: P3, c: P3, tol: f64) -> bool {
    let dot = |x: P3, y: P3| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let sub = |x: P3, y: P3| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let cross = |x: P3, y: P3| [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
    let b1 = sub(b, a);
    let b2 = sub(c, b);
    let n1 = cross(b1, b2);
    dot(n1, n1) <= tol * dot(b1, b1) * dot(b2, b2)
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("potential");
    let mut rng = Rng::new(0x9071);
    let big = f64::MAX;

    // ── lennard_jones ──
    for x in edges() {
        lj(&mut m, x, 3.4e-10, 3.8e-10);
        lj(&mut m, 1.65e-21, x, 3.8e-10);
        lj(&mut m, 1.65e-21, 3.4e-10, x);
        lj(&mut m, x, 1.0, 1.0);
        lj(&mut m, 1.0, x, 1.0);
        lj(&mut m, 1.0, 1.0, x);
    }
    let rmin_ar = lj_equilibrium_distance(3.4e-10);
    for (e, s, r) in [
        (1.65e-21, 3.4e-10, rmin_ar), // the Rust tests' values
        (1.0, 1.0, 1.0),
        (1.0, 1.0, 0.5),
        (1.0, 1.0, 100.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ] {
        lj(&mut m, e, s, r);
    }
    for k in 0..=400 {
        // r/σ from 0.5 to 100 in steps of 1/4 below 10, at Å scale (σ = 3.4 Å) and metre scale (3.4e-10 m)
        let ratio = if k <= 38 { 0.5 + 0.25 * f64::from(k) } else { 10.0 + 0.25 * f64::from(k - 38) };
        lj(&mut m, 1.65e-21, 3.4e-10, 3.4e-10 * ratio);
        lj(&mut m, 0.997, 3.4, 3.4 * ratio);
    }
    // (σ/r)¹² overflows past σ/r = 2^85.33, (σ/r)⁶ past 2^170.67; (σ/r)⁶ underflows below 2^-179, (σ/r)¹² below 2^-89.5
    for &e in &[84, 85, 86, 170, 171, 172, -89, -90, -91, -178, -179, -180, -500, -1000] {
        for &mant in &[1.0, 1.25, 1.5, 1.75] {
            let sr = mant * pow2(e);
            lj(&mut m, 1.65e-21, sr, 1.0);
            lj(&mut m, 1e300, sr, 1.0);
            lj(&mut m, 1e-300, 1.0, 1.0 / sr);
        }
    }
    for _ in 0..600 {
        let e = rng.log_uniform(-300.0, 300.0);
        let s = rng.log_uniform(-200.0, 200.0);
        let r = s * rng.log_uniform(-3.0, 3.0);
        lj(&mut m, e, s, r);
    }
    // ADR 0002 #17: ε ≥ 2¹⁰²², where Rust's 4·ε overflows (σ = r gives d = 0, Rust's ∞·0 = NaN)
    let eps_edge = [pow2(1022), pow2(1022).next_down(), (big / 4.0).next_up(), big / 4.0, big, 1.5 * pow2(1022)];
    for &e in &eps_edge {
        for &ratio in &[1.0, 0.5, 2.0, 1.122_462_048_309_373, 0.9, 1.1, 10.0, 0.1, 1e-30, 1e30] {
            lj(&mut m, e, 1.0, ratio);
            lj(&mut m, e, 3.4e-10, 3.4e-10 * ratio);
        }
    }
    for _ in 0..300 {
        let e = pow2(1022) * rng.uniform(1.0, 4.0);
        let s = rng.log_uniform(-12.0, 3.0);
        let r = s * rng.log_uniform(-1.0, 1.0);
        lj(&mut m, e, s, r);
    }

    // ── lj_equilibrium_distance ──
    for x in edges() {
        rmin(&mut m, x);
    }
    for x in decades(-323, 308, &["1", "3.4"]) {
        rmin(&mut m, x);
    }
    for k in -1074..=1023 {
        rmin(&mut m, pow2(k));
    }
    // the overflow end: σ·2^(1/6) passes DBL_MAX near DBL_MAX/1.1225
    let top = big / 1.122_462_048_309_373;
    let mut t = top;
    for _ in 0..40 {
        t = t.next_down();
    }
    for _ in 0..81 {
        rmin(&mut m, t);
        t = t.next_up();
    }
    for _ in 0..3000 {
        let x = f64::from_bits(rng.int(1, 0x7FEF_FFFF_FFFF_FFFF).cast_unsigned());
        rmin(&mut m, if rng.unit() < 0.5 { -x } else { x });
    }
    for _ in 0..500 {
        rmin(&mut m, rng.uniform(1e-10, 5e-10)); // Å-scale collision diameters in metres
    }

    // ── morse ──
    for x in edges() {
        morse_row(&mut m, x, 1.0, 1.5, 1.6);
        morse_row(&mut m, 5.0, x, 1.5, 1.6);
        morse_row(&mut m, 5.0, 1.0, x, 1.6);
        morse_row(&mut m, 5.0, 1.0, 1.5, x);
    }
    for (d, a, re, r) in [
        (5.0, 1.0, 1.5, 1.5), // the Rust tests' values
        (5.0, 1.0, 1.5, 100.0),
        (5.0, 2.0, 1.5, 1.6),
        (5.0, 2.0, 1.5, 1.4),
        (f64::NAN, 1.0, 1.5, 1.5),
        (5.0, f64::NAN, 1.5, 1.5),
        (f64::INFINITY, 1.0, 1.5, 1.5),
    ] {
        morse_row(&mut m, d, a, re, r);
    }
    // HCl at metre scale: a = 1.868e10 m⁻¹, r_e = 1.2746e-10 m, D_e = 7.399e-19 J, r from 0.5 r_e to 5 r_e
    for k in 0..=180 {
        let r = 1.2746e-10 * (0.5 + 0.025 * f64::from(k));
        morse_row(&mut m, 7.399e-19, 1.868e10, 1.2746e-10, r);
    }
    // the exponent −a(r − r_e) on targets across exp's range (subnormal, under- and overflowing e^x)
    for &x in &[-1000.0, -745.2, -745.1, -744.0, -720.0, -708.5, -708.39, -100.0, -1.0, 1.0, 100.0, 700.0, 709.78, 709.8, 1000.0] {
        for &a in &[1.0, 1e10, 1e-10, 3.0] {
            let r = 1.0 - x / a; // −a(r − 1) ≈ x
            morse_row(&mut m, 7.399e-19, a, 1.0, r);
            morse_row(&mut m, 1e300, a, 1.0, r);
        }
    }
    // near-cancellation: 1 − e^y with |y| = |a(r − r_e)| from 10⁻¹⁷ to 1 (both signs), where only a correctly rounded exp
    // gives the corrected oracle's x and so V ≈ D_e·y²
    for _ in 0..800 {
        let y = rng.log_uniform(-17.0, 0.0);
        let dr = if rng.unit() < 0.5 { y } else { -y };
        let re = 1.2746e-10;
        let a = 1.868e10;
        morse_row(&mut m, 7.399e-19, a, re, re + dr / a);
        morse_row(&mut m, 5.0, 1.0, 1.5, 1.5 + dr);
    }
    for _ in 0..600 {
        let d = rng.log_uniform(-300.0, 300.0);
        let a = rng.log_uniform(-15.0, 15.0);
        let re = rng.log_uniform(-12.0, 3.0);
        let r = re * rng.uniform(0.0, 5.0);
        morse_row(&mut m, d, a, re, r);
    }

    // ── morse_width_parameter ──
    for x in edges() {
        width(&mut m, x, 1.6267e-27, 7.399e-19);
        width(&mut m, 5.634e14, x, 7.399e-19);
        width(&mut m, 5.634e14, 1.6267e-27, x);
    }
    for (w, mu, d) in [
        (5.634e14, 1.6267e-27, 7.399e-19), // the Rust tests' HCl values
        ((516.0_f64 / 1.6267e-27).sqrt(), 1.6267e-27, 7.399e-19),
        (f64::NAN, 1.0, 1.0),
        (1.0, f64::NAN, 1.0),
        (1.0, 1.0, f64::INFINITY),
    ] {
        width(&mut m, w, mu, d);
    }
    for x in decades(-300, 300, &["1", "5.634"]) {
        width(&mut m, x, 1.6267e-27, 7.399e-19);
        width(&mut m, 5.634e14, x, 7.399e-19);
        width(&mut m, 5.634e14, 1.6267e-27, x);
    }
    for _ in 0..400 {
        width(&mut m, rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0));
    }
    // ADR 0002 #18: D_e ≥ 2¹⁰²³, where Rust's 2·D_e overflows (μ small enough that μ·½ is inexact included)
    let d_edge = [pow2(1023), pow2(1023).next_down(), big, 1.5 * pow2(1023)];
    for &d in &d_edge {
        for &mu in &[1.0, 1.6267e-27, 1e300, big, pow2(-1021), pow2(-1022), f64::from_bits(3), f64::from_bits(1), 1e-300] {
            for &w in &[1.0, big, 5.634e14, 1e-300] {
                width(&mut m, w, mu, d);
            }
        }
    }
    for _ in 0..300 {
        let d = pow2(1023) * rng.uniform(1.0, 2.0);
        width(&mut m, rng.log_uniform(-300.0, 308.0), rng.log_uniform(-320.0, 308.0), d);
    }

    // ── coulomb, coulomb_force ──
    let e = 1.602_176_634e-19;
    for x in edges() {
        coul(&mut m, x, e, 1e-10);
        coul(&mut m, e, x, 1e-10);
        coul(&mut m, e, e, x);
    }
    for (q1, q2, r) in [
        (e, -e, 1e-10), // the Rust tests' values
        (e, e, 1e-10),
        (e, e, 2e-10),
        (1.0, 1.0, 0.0),
        (1.0, 1.0, f64::NAN),
        (1.0, 1.0, f64::NEG_INFINITY),
        (0.0, -0.0, 1.0),
        (-0.0, -0.0, 1.0),
    ] {
        coul(&mut m, q1, q2, r);
    }
    for r in decades(-320, 308, &["1", "2.5"]) {
        coul(&mut m, e, -e, r);
    }
    // the intermediates: k·q₁ overflow (|q₁| > 2e298), k·q₁·q₂ underflow, r·r over- (r > 1.34e154) and underflow
    for &(q1, q2, r) in &[
        (1e299, 1e-299, 1.0),
        (1e299, 1e-10, 1e10),
        (1e-200, 1e-200, 1e-100),
        (1e-160, 1e-160, 1e-10),
        (e, e, 1e155),
        (1e100, 1e100, 1e155),
        (e, e, 1e-160),
        (1e-150, 1e-150, 1e-160),
        (big, 1.0, big),
        (big, -big, big),
    ] {
        coul(&mut m, q1, q2, r);
    }
    for _ in 0..600 {
        let q1 = rng.signed_log_uniform(-200.0, 200.0);
        let q2 = rng.signed_log_uniform(-200.0, 200.0);
        coul(&mut m, q1, q2, rng.log_uniform(-200.0, 200.0));
    }

    // ── distance_3d ──
    let o: P3 = [0.0; 3];
    dist(&mut m, o, [3.0, 4.0, 0.0]); // the Rust tests' points
    dist(&mut m, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
    for x in edges() {
        for k in 0..6 {
            let mut v = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0];
            v[k] = x;
            dist(&mut m, [v[0], v[1], v[2]], [v[3], v[4], v[5]]);
        }
    }
    for &s in &[1.0, 1e-10, 1e-15, 1e10, 1e150, 1e155, 1e-160, 1e-170, pow2(-1050), 1e300, pow2(511), pow2(512)] {
        dist(&mut m, scale([1.0, 2.0, 2.0], s), o);
        dist(&mut m, scale([-1.0, 0.5, 2.0], s), scale([2.0, -3.0, 0.25], s));
        dist(&mut m, scale([3.0, 4.0, 12.0], s), o); // 13: exact
    }
    for _ in 0..500 {
        let s = rng.log_uniform(-12.0, 3.0);
        dist(&mut m, rpt(&mut rng, s), rpt(&mut rng, s));
    }
    for _ in 0..200 {
        let s = rng.log_uniform(-300.0, 300.0);
        dist(&mut m, rpt(&mut rng, s), rpt(&mut rng, s));
    }

    // ── bond_angle ──
    let x1: P3 = [1.0, 0.0, 0.0];
    let y1: P3 = [0.0, 1.0, 0.0];
    angle(&mut m, x1, o, y1); // the Rust tests' geometries: 90°, 180°, coincident, NaN / ∞
    angle(&mut m, x1, o, [-1.0, 0.0, 0.0]);
    angle(&mut m, o, o, x1);
    angle(&mut m, [1.0 + 1e-15, 1.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]);
    angle(&mut m, [f64::NAN, 0.0, 0.0], o, x1);
    angle(&mut m, x1, o, [f64::INFINITY, 0.0, 0.0]);
    // water in metres (the Rust test builds H₂ from cos θ and sin θ for θ = 104.52°; here they are typed as decimals,
    // since a grid input never calls a transcendental and only needs to be a fixed double), femtometre scale, the
    // small-angle test
    let r_oh = 0.9572e-10;
    let (c104, s104) = (-0.250_380_004_394_034_9_f64, 0.968_147_640_378_108_1_f64);
    angle(&mut m, [r_oh, 0.0, 0.0], o, [r_oh * c104, r_oh * s104, 0.0]);
    angle(&mut m, [1e-15, 0.0, 0.0], o, [0.0, 1e-15, 0.0]);
    angle(&mut m, x1, o, [0.999_999_999_999_999_9, 1e-9, 0.0]);
    for k in 0..12 {
        let t = pow2(-4 * k - 1); // tiny angles: the point (1, t) at angle ≈ t
        angle(&mut m, x1, o, [1.0, t, 0.0]);
        angle(&mut m, x1, o, [-1.0, t, 0.0]); // near π
    }
    for x in edges() {
        for k in 0..9 {
            let mut v = [1.0, 0.25, 0.0, 0.0, 0.0, 0.0, -0.5, 1.0, 0.0];
            v[k] = x;
            angle(&mut m, [v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]);
        }
    }
    for &s in &[1.0, 1e-10, 1e-15, 1e10, 1e100, 1e150, 1e155, 1e160, 1e-150, 1e-160, 1e-165, pow2(-1050), pow2(511), pow2(512)] {
        angle(&mut m, scale([1.0, 0.25, 0.0], s), o, scale([-0.5, 1.0, 0.0], s));
        angle(&mut m, scale([2.0, 1.0, 1.0], s), scale([1.0, 1.0, 1.0], s), scale([1.0, 2.0, -1.0], s));
    }
    // −0.0 coordinates (the dot product's sign of zero, atan2's)
    angle(&mut m, [1.0, -0.0, -0.0], [-0.0, -0.0, -0.0], [-0.0, 1.0, -0.0]);
    angle(&mut m, [-1.0, -0.0, 0.0], [0.0, -0.0, -0.0], [1.0, 0.0, -0.0]);
    angle(&mut m, [-0.0, -1.0, 0.0], o, [-0.0, 1.0, -0.0]);
    // on the coincidence tolerance to the ulp: B = (B₁, B₂, 0) of a random size and direction, A = B + (0, 0, δ) (so
    // A − B = (0, 0, δ) exactly), δ scanned ulp by ulp across Rust's threshold δ² = RN(10⁻²⁴·|A|²): the rows on either
    // side of the crossing, and every row where 10⁻¹² ± 1 ulp moves Rust's test (the rows just below it move with 10⁻¹² × ½
    // and × 2 too); A coincident with B, and C coincident with B instead
    let tol = 1e-12_f64;
    for _ in 0..80 {
        let sz = rng.log_uniform(-12.0, 4.0);
        let b: P3 = [rng.uniform(-sz, sz), rng.uniform(-sz, sz), 0.0];
        let bn = (b[0] * b[0] + b[1] * b[1]).sqrt();
        let mut d = bn * tol;
        for _ in 0..200 {
            d = d.next_down();
        }
        let far: P3 = [b[0] + sz, b[1] - 0.5 * sz, sz];
        let mut prev = None;
        let mut emitted = 0;
        for _ in 0..400 {
            let a: P3 = [b[0], b[1], d];
            let base = coincident(a, b, tol);
            let crossing = prev.is_some_and(|p| p != base);
            if crossing || base != coincident(a, b, tol.next_up()) || base != coincident(a, b, tol.next_down()) {
                let mut dd = d;
                for _ in 0..2 {
                    dd = dd.next_down();
                }
                for _ in 0..5 {
                    angle(&mut m, [b[0], b[1], dd], b, far);
                    angle(&mut m, far, b, [b[0], b[1], dd]);
                    dd = dd.next_up();
                }
                emitted += 1;
            }
            prev = Some(base);
            d = d.next_up();
            if emitted >= 3 {
                break;
            }
        }
    }
    for _ in 0..800 {
        let s = rng.log_uniform(-15.0, 2.0);
        let b = rpt(&mut rng, s);
        angle(&mut m, add(b, rpt(&mut rng, s)), b, add(b, rpt(&mut rng, s)));
    }
    for _ in 0..200 {
        let s = rng.log_uniform(-170.0, 170.0);
        let b = rpt(&mut rng, s);
        angle(&mut m, add(b, rpt(&mut rng, s)), b, add(b, rpt(&mut rng, s)));
    }

    // ── dihedral_angle ──
    let b0: P3 = [1.0, 0.0, 0.0];
    dihedral(&mut m, o, b0, [2.0, 1.0, 0.0], [3.0, 1.0, 0.0]); // the Rust tests' geometries
    dihedral(&mut m, o, b0, [2.0, -1.0, 0.0], [3.0, -1.0, 0.0]);
    dihedral(&mut m, [0.0, 1.0, 0.0], o, b0, [1.0, 1.0, 0.0]);
    dihedral(&mut m, b0, o, y1, [0.0, 1.0, 1.0]);
    dihedral(&mut m, b0, o, y1, [0.0, 1.0, -1.0]);
    let gauche = [[1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.5, 0.0, 1.0]];
    dihedral(&mut m, gauche[0], gauche[1], gauche[2], gauche[3]);
    let gm = gauche.map(|p| p.map(|x| x * 1e-10));
    dihedral(&mut m, gm[0], gm[1], gm[2], gm[3]);
    dihedral(&mut m, o, b0, [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]);
    dihedral(&mut m, o, [1e-10, 0.0, 0.0], [2e-10, 1e-10, 0.0], [3e-10, 2e-10, 0.0]);
    dihedral(&mut m, [f64::NAN; 3], o, b0, [1.0, 1.0, 0.0]);
    // trans, cis, gauche ±60°, ±90°, ±120° about the x axis (the 1, √3 points), their mirror images, with +0.0 and −0.0
    // z terms (the y-term −0.0 that Rust maps from −π to π)
    let h = 3.0_f64.sqrt();
    let ends: [(f64, f64); 8] = [(1.0, 0.0), (-1.0, 0.0), (0.5, 0.5 * h), (0.5, -0.5 * h), (0.0, 1.0), (0.0, -1.0), (-0.5, 0.5 * h), (-0.5, -0.5 * h)];
    for &(cy, cz) in &ends {
        for &zs in &[0.0, -0.0] {
            for &s in &[1.0, 1e-10, 1e-15, 1e10] {
                let a = scale([-1.0, 1.0, zs], s);
                let d = scale([2.5, cy, if cz == 0.0 { zs } else { cz }], s);
                dihedral(&mut m, a, o, scale([1.5, 0.0, zs], s), d);
                dihedral(&mut m, scale([-1.0, -1.0, zs], s), o, scale([1.5, -0.0, -zs], s), scale([2.5, -cy, -cz], s));
            }
        }
    }
    // every sign of the zero coordinates of the Rust test's two planar trans geometries and of a cis one (7 zeros
    // each): where b₁·n₂ comes out −0.0, atan2 gives −π for trans (Rust reports π) and −0.0 for cis
    for mask in 0..128_u32 {
        let z = |k: u32| if mask >> k & 1 == 1 { -0.0 } else { 0.0 };
        for &yc in &[1.0, -1.0] {
            dihedral(&mut m, [z(0), z(1), z(2)], [1.0, z(3), z(4)], [2.0, yc, z(5)], [3.0, yc, z(6)]);
        }
        dihedral(&mut m, [z(0), 1.0, z(1)], [z(2), z(3), z(4)], [1.0, z(5), z(6)], [1.0, 1.0, 0.0]);
    }
    for x in edges() {
        for k in 0..12 {
            let mut v = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 2.5, 0.0, 1.0];
            v[k] = x;
            dihedral(
                &mut m,
                [v[0], v[1], v[2]],
                [v[3], v[4], v[5]],
                [v[6], v[7], v[8]],
                [v[9], v[10], v[11]],
            );
        }
    }
    for &s in &[1e-77, 1e-80, 1e-100, 1e70, 1e77, 1e80, 1e100, pow2(-1050), 1e-160, 1e155] {
        let g = gauche.map(|p| scale(p, s));
        dihedral(&mut m, g[0], g[1], g[2], g[3]);
    }
    // consecutive atoms coinciding (b₁, b₂ or b₃ zero)
    dihedral(&mut m, o, o, b0, [1.0, 1.0, 0.0]);
    dihedral(&mut m, [0.0, 1.0, 0.0], b0, b0, [2.0, 1.0, 0.0]);
    dihedral(&mut m, [0.0, 1.0, 0.0], o, b0, b0);
    // on the collinearity tolerance to the ulp: A = (−L, h, 0), B = 0, C = (L, 0, 0) (so |b₁ × b₂|² = (hL)² rounded,
    // monotone in h), h scanned ulp by ulp across Rust's threshold near 10⁻¹⁰·L: the rows on either side of the
    // crossing and every row where 10⁻²⁰ ± 1 ulp moves the test (× ½ and × 2 too); A-B-C bent, and B-C-D instead
    let ctol = 1e-20_f64;
    for _ in 0..80 {
        let len = rng.log_uniform(-11.0, 4.0);
        let mut hgt = len * 1e-10;
        for _ in 0..200 {
            hgt = hgt.next_down();
        }
        let c: P3 = [len, 0.0, 0.0];
        let mut prev = None;
        let mut emitted = 0;
        for _ in 0..400 {
            let a: P3 = [-len, hgt, 0.0];
            let base = collinear(a, o, c, ctol);
            let crossing = prev.is_some_and(|p| p != base);
            if crossing || base != collinear(a, o, c, ctol.next_up()) || base != collinear(a, o, c, ctol.next_down()) {
                let mut hh = hgt;
                for _ in 0..2 {
                    hh = hh.next_down();
                }
                for _ in 0..5 {
                    dihedral(&mut m, [-len, hh, 0.0], o, c, [2.0 * len, 0.0, len]);
                    // B-C-D bent instead: A = (−L, 0, L) and D = C + (L, h, 0), so b₂ × b₃ = (0, 0, L·h)
                    dihedral(&mut m, [-len, 0.0, len], o, c, [2.0 * len, hh, 0.0]);
                    hh = hh.next_up();
                }
                emitted += 1;
            }
            prev = Some(base);
            hgt = hgt.next_up();
            if emitted >= 3 {
                break;
            }
        }
    }
    for _ in 0..800 {
        let s = rng.log_uniform(-15.0, 2.0);
        let b = rpt(&mut rng, s);
        let c = add(b, rpt(&mut rng, s));
        dihedral(&mut m, add(b, rpt(&mut rng, s)), b, c, add(c, rpt(&mut rng, s)));
    }
    for _ in 0..200 {
        let s = rng.log_uniform(-80.0, 80.0);
        let b = rpt(&mut rng, s);
        let c = add(b, rpt(&mut rng, s));
        dihedral(&mut m, add(b, rpt(&mut rng, s)), b, c, add(c, rpt(&mut rng, s)));
    }
    m
}
