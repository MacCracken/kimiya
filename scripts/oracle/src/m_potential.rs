//! `potential.tsv` — intermolecular potentials (Lennard-Jones, Morse, Coulomb) and molecular geometry from 3D
//! coordinates (distance, bond and dihedral angles) (port plan §4, M4d).
//!
//! Rust's `[f64; 3]` points are flattened into their three `f:` tokens, in argument order (`distance_3d` has 6 input
//! floats, `bond_angle` 9, `dihedral_angle` 12).
//! - Class B (+ − × ÷ and √, expected 0 ulp): `lennard_jones` (each argument over the edges, the Rust tests' values,
//!   r/σ from 0.5 to 100 at Å and metre scale, (σ/r)¹² and (σ/r)⁶ over- and underflowing, seeded; ε ≥ 2¹⁰²², where
//!   Rust's 4·ε overflows: ADR 0002 #17's rows), `morse_width_parameter` (the edges, the Rust tests' HCl values, every
//!   decade, seeded; D_e ≥ 2¹⁰²³, where Rust's 2·D_e overflows: #18's rows), `coulomb`, `coulomb_force` (the edges, the
//!   Rust tests' charges, every decade of r, signed charges, the intermediates' over- and underflow, seeded),
//!   `distance_3d` (each coordinate over the edges, the Rust tests' points, Å and metre scale, over- and underflowing
//!   squares, seeded). **ADR 0002 #32** (the overflow audit's B3, which folds #17 and #18 in): every `lennard_jones`,
//!   `morse_width_parameter`, `coulomb` and `coulomb_force` row is noted where the port's value differs from Rust's —
//!   `lennard_jones` Rust's sequence with an unbounded exponent and the final product rounded once (#17's (ε·d)·4 where
//!   ε ≥ 2¹⁰²², the same bits), the others rounded once (`audit.rs`) — and `audit::emit_potential` appends their edge
//!   rows with its own seed.
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
//! - **ADR 0002 #20** (M4e-2): `distance_3d`, `bond_angle` and `dihedral_angle` are Rust's sequence evaluated in binary64
//!   with an unbounded exponent range (every operation rounded to 53 bits, nothing over- or underflowing; the final √ or
//!   atan2 rounded once into binary64). Where Rust's chain is not that evaluation (an intermediate left the normal range)
//!   and the result differs, the row keeps the oracle's value and its note names the port's, computed here in exact
//!   integer arithmetic (`exact::Ux`, `exact::ux_atan2_f64`); every row's replica of Rust's chain is checked against the
//!   oracle. The edge grids: coordinates near 2^±(500…1074), mixed magnitudes and subnormal coordinates, differences that
//!   overflow, the tolerance scans at extreme scales, tiny angles with subnormal results, trans and cis dihedrals at
//!   scale, seeded geometries at random scales.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::audit;
use crate::exact::{Ux, cr_root6_2_mul, ux_atan2_f64};
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, f};
use kimiya::potential::{bond_angle, dihedral_angle, distance_3d, lj_equilibrium_distance, morse};

type P3 = [f64; 3];

fn toks(ps: &[P3]) -> Vec<Tok> {
    ps.iter().flat_map(|p| p.iter().map(|&x| f(x))).collect()
}

fn scale(p: P3, s: f64) -> P3 {
    [p[0] * s, p[1] * s, p[2] * s]
}

fn lj(m: &mut Module, eps: f64, sigma: f64, r: f64) {
    audit::lj(m, eps, sigma, r);
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
    audit::width(m, w, mu, d);
}

fn coul(m: &mut Module, q1: f64, q2: f64, r: f64) {
    audit::coul(m, q1, q2, r);
}

// ── ADR 0002 #20: Rust's geometry sequences, on f64 (the replica, checked against the oracle on every row) and with an
// unbounded exponent (exact integer arithmetic) ──

fn vsub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vdot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vcross(a: P3, b: P3) -> P3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

type U3 = [Ux; 3];

fn u3(p: P3) -> U3 {
    p.map(Ux::from_f64)
}

fn usub(a: U3, b: U3) -> U3 {
    [a[0].sub(b[0]), a[1].sub(b[1]), a[2].sub(b[2])]
}

fn udot(a: U3, b: U3) -> Ux {
    a[0].mul(b[0]).add(a[1].mul(b[1])).add(a[2].mul(b[2]))
}

fn ucross(a: U3, b: U3) -> U3 {
    [
        a[1].mul(b[2]).sub(a[2].mul(b[1])),
        a[2].mul(b[0]).sub(a[0].mul(b[2])),
        a[0].mul(b[1]).sub(a[1].mul(b[0])),
    ]
}

fn finite(ps: &[P3]) -> bool {
    ps.iter().flatten().all(|x| x.is_finite())
}

const COINCIDENT_REL_TOL: f64 = 1e-12;
const COLLINEAR_SIN_SQ_TOL: f64 = 1e-20;

/// Rust's `bond_angle` up to atan2's arguments: `None` for its Err, `Some((y, x))` otherwise.
fn bond_parts(a: P3, b: P3, c: P3, tol: f64) -> Option<(f64, f64)> {
    let ok = |x: f64| x > 0.0 && x.is_finite();
    let (ba, bc) = (vsub(a, b), vsub(c, b));
    let (ba_sq, bc_sq) = (vdot(ba, ba), vdot(bc, bc));
    let tol_sq = tol * tol;
    let b_sq = vdot(b, b);
    if !ok(ba_sq) || !ok(bc_sq) || ba_sq <= tol_sq * vdot(a, a).max(b_sq) || bc_sq <= tol_sq * vdot(c, c).max(b_sq) {
        return None;
    }
    let cr = vcross(ba, bc);
    Some((vdot(cr, cr).sqrt(), vdot(ba, bc)))
}

/// The same sequence with an unbounded exponent.
fn bond_parts_ux(a: P3, b: P3, c: P3, tol: f64) -> Option<(Ux, Ux)> {
    let (ua, ub, uc) = (u3(a), u3(b), u3(c));
    let (ba, bc) = (usub(ua, ub), usub(uc, ub));
    let (ba_sq, bc_sq) = (udot(ba, ba), udot(bc, bc));
    let t = Ux::from_f64(tol);
    let tol_sq = t.mul(t);
    let b_sq = udot(ub, ub);
    if !ba_sq.is_positive()
        || !bc_sq.is_positive()
        || ba_sq.le(tol_sq.mul(udot(ua, ua).max(b_sq)))
        || bc_sq.le(tol_sq.mul(udot(uc, uc).max(b_sq)))
    {
        return None;
    }
    let cr = ucross(ba, bc);
    Some((udot(cr, cr).sqrt(), udot(ba, bc)))
}

/// Rust's `dihedral_angle` up to atan2's arguments: `None` for its Err, `Some((y, x))` otherwise.
fn dihedral_parts(a: P3, b: P3, c: P3, d: P3, tol: f64) -> Option<(f64, f64)> {
    let ok = |x: f64| x > 0.0 && x.is_finite();
    let (b1, b2, b3) = (vsub(b, a), vsub(c, b), vsub(d, c));
    let (n1, n2) = (vcross(b1, b2), vcross(b2, b3));
    let (n1_sq, n2_sq, b2_sq) = (vdot(n1, n1), vdot(n2, n2), vdot(b2, b2));
    if !ok(n1_sq) || !ok(n2_sq) || n1_sq <= tol * vdot(b1, b1) * b2_sq || n2_sq <= tol * b2_sq * vdot(b3, b3) {
        return None;
    }
    Some((b2_sq.sqrt() * vdot(b1, n2), vdot(n1, n2)))
}

/// The same sequence with an unbounded exponent.
fn dihedral_parts_ux(a: P3, b: P3, c: P3, d: P3, tol: f64) -> Option<(Ux, Ux)> {
    let (ua, ub, uc, ud) = (u3(a), u3(b), u3(c), u3(d));
    let (b1, b2, b3) = (usub(ub, ua), usub(uc, ub), usub(ud, uc));
    let (n1, n2) = (ucross(b1, b2), ucross(b2, b3));
    let (n1_sq, n2_sq, b2_sq) = (udot(n1, n1), udot(n2, n2), udot(b2, b2));
    let t = Ux::from_f64(tol);
    if !n1_sq.is_positive()
        || !n2_sq.is_positive()
        || n1_sq.le(t.mul(udot(b1, b1)).mul(b2_sq))
        || n2_sq.le(t.mul(b2_sq).mul(udot(b3, b3)))
    {
        return None;
    }
    Some((b2_sq.sqrt().mul(udot(b1, n2)), udot(n1, n2)))
}

/// Rust's −π → π (dihedral_angle reports (−π, π]).
fn dihedral_map(phi: f64) -> f64 {
    if phi <= -std::f64::consts::PI { std::f64::consts::PI } else { phi }
}

/// ADR 0002 #20: where Rust's chain is not the unbounded-exponent evaluation (an Err decision or atan2's arguments
/// differ), the port's result (`None` for Err), else `None` from the outer Option (the port's value is Rust's sequence
/// with a correctly rounded atan2, the corrected oracle's).
fn port_angle(rust: Option<(f64, f64)>, ux: Option<(Ux, Ux)>, map: fn(f64) -> f64) -> Option<Option<f64>> {
    let same = match (rust, ux) {
        (None, None) => true,
        (Some((y, x)), Some((uy, ux_))) => uy.is(y) && ux_.is(x),
        _ => false,
    };
    if same {
        return None;
    }
    Some(ux.map(|(uy, ux_)| map(ux_atan2_f64(uy, ux_))))
}

const P20: &str = "ADR 0002 #20: the port returns Rust's sequence evaluated with an unbounded exponent";

/// The #20 note for an angle row where the port's result differs from the oracle's `out`, or `None`.
fn angle_note(out: &Result<f64, kimiya::error::KimiyaError>, port: Option<Option<f64>>) -> Option<String> {
    let port = port?;
    match (out.as_ref().ok(), port) {
        (None, None) => None,
        (Some(&v), Some(x)) if v.to_bits() == x.to_bits() => None,
        (_, Some(x)) => Some(format!("{P20} (atan2 rounded once) = Ok({x:?}) (f:{:016x})", x.to_bits())),
        (Some(_), None) => Some(format!("{P20} = Err(InvalidInput) (E:InvalidInput)")),
    }
}

fn dist(m: &mut Module, p1: P3, p2: P3) {
    let d = distance_3d(bb(p1), bb(p2));
    let name = "potential::distance_3d";
    let ins = toks(&[p1, p2]);
    let dd = vsub(p2, p1);
    let s = vdot(dd, dd);
    assert!(s.sqrt().to_bits() == d.to_bits() || d.is_nan(), "distance_3d replica");
    if finite(&[p1, p2]) {
        let du = usub(u3(p2), u3(p1));
        let us = udot(du, du);
        let x = us.sqrt_f64();
        if !us.is(s) && x.to_bits() != d.to_bits() {
            let note = format!("{P20} (√ rounded once) = {x:?} (f:{:016x})", x.to_bits());
            m.row_noted(name, Class::B, &ins, &[f(d)], &note);
            return;
        }
    }
    m.row(name, Class::B, &ins, &[f(d)]);
}

fn angle(m: &mut Module, a: P3, b: P3, c: P3) {
    let out = bond_angle(bb(a), bb(b), bb(c));
    let name = "potential::bond_angle";
    let ins = toks(&[a, b, c]);
    let rust = bond_parts(a, b, c, COINCIDENT_REL_TOL);
    match (&out, rust) {
        (Err(_), None) => {}
        (Ok(v), Some((y, x))) => assert!(v.to_bits() == y.atan2(x).to_bits(), "bond_angle replica"),
        _ => panic!("bond_angle replica: {a:?} {b:?} {c:?}"),
    }
    if finite(&[a, b, c]) {
        let port = port_angle(rust, bond_parts_ux(a, b, c, COINCIDENT_REL_TOL), |x| x);
        if let Some(note) = angle_note(&out, port) {
            m.row_noted(name, Class::C, &ins, &[res(out)], &note);
            return;
        }
    }
    m.row(name, Class::C, &ins, &[res(out)]);
}

fn dihedral(m: &mut Module, a: P3, b: P3, c: P3, d: P3) {
    let out = dihedral_angle(bb(a), bb(b), bb(c), bb(d));
    let name = "potential::dihedral_angle";
    let ins = toks(&[a, b, c, d]);
    let rust = dihedral_parts(a, b, c, d, COLLINEAR_SIN_SQ_TOL);
    match (&out, rust) {
        (Err(_), None) => {}
        (Ok(v), Some((y, x))) => assert!(v.to_bits() == dihedral_map(y.atan2(x)).to_bits(), "dihedral replica"),
        _ => panic!("dihedral_angle replica: {a:?} {b:?} {c:?} {d:?}"),
    }
    if finite(&[a, b, c, d]) {
        let port = port_angle(rust, dihedral_parts_ux(a, b, c, d, COLLINEAR_SIN_SQ_TOL), dihedral_map);
        if let Some(note) = angle_note(&out, port) {
            m.row_noted(name, Class::C, &ins, &[res(out)], &note);
            return;
        }
    }
    m.row(name, Class::C, &ins, &[res(out)]);
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
    let ba = vsub(a, b);
    vdot(ba, ba) <= tol * tol * vdot(a, a).max(vdot(b, b))
}

/// The same test with an unbounded exponent (ADR 0002 #20: the port's, at any scale).
fn coincident_ux(a: P3, b: P3, tol: f64) -> bool {
    let (ua, ub) = (u3(a), u3(b));
    let ba = usub(ua, ub);
    let t = Ux::from_f64(tol);
    udot(ba, ba).le(t.mul(t).mul(udot(ua, ua).max(udot(ub, ub))))
}

/// Rust's dihedral collinearity test for A-B-C (`n1_sq <= tol · |b1|² · |b2|²`), with the tolerance `tol`.
fn collinear(a: P3, b: P3, c: P3, tol: f64) -> bool {
    let (b1, b2) = (vsub(b, a), vsub(c, b));
    let n1 = vcross(b1, b2);
    vdot(n1, n1) <= tol * vdot(b1, b1) * vdot(b2, b2)
}

/// The same test with an unbounded exponent (ADR 0002 #20).
fn collinear_ux(a: P3, b: P3, c: P3, tol: f64) -> bool {
    let (b1, b2) = (usub(u3(b), u3(a)), usub(u3(c), u3(b)));
    let n1 = ucross(b1, b2);
    udot(n1, n1).le(Ux::from_f64(tol).mul(udot(b1, b1)).mul(udot(b2, b2)))
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
    edges20(&mut m);
    audit::emit_potential(&mut m);
    m
}

/// x·2^k for any k, in at most two exact scalings (the second one rounds once where the result is subnormal).
fn sc2(x: f64, k: i32) -> f64 {
    if k > 1023 {
        x * pow2(1023) * pow2(k - 1023)
    } else if k < -1022 {
        x * pow2(-1022) * pow2((k + 1022).max(-1074))
    } else {
        x * pow2(k)
    }
}

fn sc(p: P3, k: i32) -> P3 {
    p.map(|x| sc2(x, k))
}

/// A double with a uniformly random binade in [2^lo, 2^hi] (lo ≥ −1074, the subnormal binades included), a random
/// significand and sign.
fn binade(rng: &mut Rng, lo: i32, hi: i32) -> f64 {
    let e = i32::try_from(rng.int(i64::from(lo), i64::from(hi))).unwrap_or(0);
    let sig = 1.0 + rng.unit();
    let x = if e < -1022 { f64::from_bits(u64::try_from(rng.int(1, (1_i64 << (e + 1075)) - 1)).unwrap_or(1)) } else { sig * pow2(e) };
    if rng.unit() < 0.5 { -x } else { x }
}

/// The least positive double x in (lo, hi] with `test(x)` false, for a test true at lo, false at hi and monotone
/// between (bisection on the bit patterns); `None` when the ends don't bracket a turn.
fn crossing(test: impl Fn(f64) -> bool, lo: f64, hi: f64) -> Option<f64> {
    if !test(lo) || test(hi) {
        return None;
    }
    let (mut a, mut b) = (lo.to_bits(), hi.to_bits());
    while b - a > 1 {
        let mid = a + (b - a) / 2;
        if test(f64::from_bits(mid)) { a = mid } else { b = mid }
    }
    Some(f64::from_bits(b))
}

/// The doubles two ulps either side of each crossing, and three between them (none if they coincide).
fn around_crossings(port: Option<f64>, rust: Option<f64>) -> Vec<f64> {
    let mut out = Vec::new();
    for x in [port, rust].into_iter().flatten() {
        let mut y = x.next_down().next_down();
        for _ in 0..4 {
            out.push(y);
            y = y.next_up();
        }
    }
    if let (Some(p), Some(r)) = (port, rust) {
        let (a, b) = if p.to_bits() < r.to_bits() { (p.to_bits(), r.to_bits()) } else { (r.to_bits(), p.to_bits()) };
        if b - a > 4 {
            for q in 1..4_u64 {
                out.push(f64::from_bits(a + (b - a) * q / 4));
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out.dedup_by(|x, y| x.to_bits() == y.to_bits());
    out
}

/// ADR 0002 #20 (M4e-2): the geometry where Rust's squares, cross products and dot products leave the normal range.
#[allow(clippy::too_many_lines)] // one grid after another
fn edges20(m: &mut Module) {
    let mut rng = Rng::new(0x9020);
    let big = f64::MAX;
    let o: P3 = [0.0; 3];
    let x1: P3 = [1.0, 0.0, 0.0];
    let y1: P3 = [0.0, 1.0, 0.0];
    // binary exponents near both ends: the squares (2^±511), the fourth powers (2^±255), the subnormal range
    let ks: [i32; 36] = [
        -1074, -1070, -1060, -1050, -1040, -1030, -1023, -1022, -1000, -950, -900, -800, -700, -600, -540, -538, -537,
        -520, -512, -511, -510, -500, -300, -256, -255, -230, 230, 254, 255, 256, 500, 511, 512, 600, 800, 1000,
    ];
    let ks_top: [i32; 4] = [1010, 1016, 1019, 1020];

    // ── distance_3d ──
    let pairs: [(P3, P3); 4] = [
        ([1.0, 2.0, 2.0], o),
        ([-1.0, 0.5, 2.0], [2.0, -3.0, 0.25]),
        ([3.0, 4.0, 12.0], o),
        (x1, y1),
    ];
    for &k in ks.iter().chain(&ks_top) {
        for &(p, q) in &pairs {
            dist(m, sc(p, k), sc(q, k));
        }
    }
    // one geometry at every seventh binary scale (the distance scales by exactly 2^k wherever it is normal)
    let mut k = -1074;
    while k <= 1019 {
        dist(m, sc([-1.0, 0.5, 2.0], k), sc([2.0, -3.0, 0.25], k));
        k += 7;
    }
    // differences that overflow, and their neighbours
    for (p, q) in [
        ([big, 0.0, 0.0], [-big, 0.0, 0.0]),
        ([big, big, big], [-big, -big, -big]),
        ([1e308, -1e308, 1e308], [-1e308, 1e308, -1e308]),
        ([big, 0.0, 0.0], [-big / 2.0, 0.0, 0.0]),
        ([big / 2.0, big / 2.0, 0.0], [-big / 2.0, -big / 2.0, 0.0]),
        ([big, 1.0, -1.0], [0.0, -big, big]),
        ([big, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ([big, big, 0.0], [0.0, 0.0, 0.0]),
    ] {
        dist(m, p, q);
    }
    // mixed magnitudes and subnormal coordinates, seeded
    for _ in 0..300 {
        let p = [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)];
        let q = [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)];
        dist(m, p, q);
    }
    for _ in 0..150 {
        let p = [binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023)];
        let q = [binade(&mut rng, -1074, -1023), 0.0, binade(&mut rng, -1074, -1023)];
        dist(m, p, q);
    }
    for _ in 0..150 {
        let s = rng.int(-1074, 1018);
        let k = i32::try_from(s).unwrap_or(0);
        let p = rpt(&mut rng, 1.0);
        dist(m, sc(p, k), sc(add(p, rpt(&mut rng, 1e-6)), k)); // near points at a random scale
    }

    // ── bond_angle ──
    let tris: [(P3, P3, P3); 5] = [
        ([1.0, 0.25, 0.0], o, [-0.5, 1.0, 0.0]),
        ([2.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 2.0, -1.0]),
        (x1, o, y1),
        (x1, o, [-1.0, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0009765625, 0.0]),
    ];
    for &k in ks.iter().chain(&ks_top) {
        for &(a, b, c) in &tris {
            angle(m, sc(a, k), sc(b, k), sc(c, k));
        }
    }
    // mixed magnitudes: the two bonds, or the coordinates of one point, at very different scales
    for (a, b, c) in [
        ([pow2(600), 0.0, 0.0], o, [pow2(-600), pow2(-600), 0.0]),
        (x1, o, [pow2(-700), pow2(-800), 0.0]),
        ([pow2(1000), pow2(1000), 0.0], o, [pow2(-1000), 0.0, 0.0]),
        ([pow2(1000), 1.0, 0.0], o, [1.0, pow2(-1000), 0.0]),
        ([pow2(-1074), 0.0, 0.0], o, [0.0, pow2(-1074), 0.0]),
        ([pow2(-1074), 0.0, 0.0], o, [-pow2(-1074), pow2(-1074), 0.0]),
        ([3.0 * pow2(-1074), pow2(-1074), 0.0], o, [pow2(-1073), 5.0 * pow2(-1074), 0.0]),
        ([pow2(-1074), pow2(-1074), pow2(-1074)], o, [pow2(-1074), 0.0, 0.0]),
        ([big, 0.0, 0.0], [-big, 0.0, 0.0], [-big, big, 0.0]),
        ([big, big, 0.0], o, [-big, big, 0.0]),
        ([big, 0.0, 0.0], o, [0.0, big, 0.0]),
        ([big, 0.0, 0.0], [big / 2.0, 0.0, 0.0], [big, big / 4.0, 0.0]),
        ([1e300, 1e-300, 0.0], [0.0, 0.0, 0.0], [1e-300, 1e300, 0.0]),
    ] {
        angle(m, a, b, c);
        angle(m, c, b, a);
    }
    // tiny angles (results near and in the subnormal range) and angles near π: A = (1, 0, 0)·2^s, C = (±1, t, 0)·2^s
    for &s in &[0, -300, -500, 300, 500, 1000] {
        for &t in &[-1022, -1030, -1040, -1050, -1060, -1070, -1074, -1075, -1076, -1100, -1500] {
            let tv = if t + s >= -1074 { sc2(1.0, t + s) } else { 0.0 };
            let tv3 = if t + s >= -1072 { sc2(3.0, t + s) } else { tv };
            for &cx in &[1.0, -1.0] {
                angle(m, sc(x1, s), o, [sc2(cx, s), tv, 0.0]);
                angle(m, sc(x1, s), o, [sc2(cx, s), tv3, 0.0]);
            }
        }
    }
    for _ in 0..120 {
        let s = i32::try_from(rng.int(-500, 1000)).unwrap_or(0);
        let t = binade(&mut rng, -1074, -1000).abs();
        angle(m, sc(x1, s), o, [sc2(1.0, s), sc2(t, s.min(0)), 0.0]);
    }
    // the coincidence tolerance to the ulp at extreme scales: B and A − B = (0, 0, δ) scaled by 2^k exactly, δ found by
    // bisection where the port's exact test turns (the unit-scale one) and where Rust's turns (its squares over- or
    // underflow, or are rounded at the subnormal grid for k near −480 … −520, or |b|² overflows near 2^512), the rows
    // around both and between them, where the two disagree
    for &k in &[-1000, -900, -700, -600, -540, -520, -500, -490, -480, -470, 508, 511, 512, 514, 600, 700, 900, 1000] {
        for _ in 0..3 {
            let sz = rng.log_uniform(-1.0, 1.0);
            let b0: P3 = [rng.uniform(-sz, sz), rng.uniform(-sz, sz), 0.0];
            let far0: P3 = [b0[0] + sz, b0[1] - 0.5 * sz, sz];
            let bn = (b0[0] * b0[0] + b0[1] * b0[1]).sqrt();
            let (b, far) = (sc(b0, k), sc(far0, k));
            let at = |d: f64| sc([b0[0], b0[1], d], k);
            let (lo, hi) = (bn * COINCIDENT_REL_TOL * 0.25, bn * COINCIDENT_REL_TOL * 4.0);
            let port = crossing(|d| coincident_ux(at(d), b, COINCIDENT_REL_TOL), lo, hi);
            let rust = crossing(|d| bond_parts(at(d), b, far, COINCIDENT_REL_TOL).is_none(), lo, hi);
            for d in around_crossings(port, rust) {
                angle(m, at(d), b, far);
                angle(m, far, b, at(d));
            }
        }
    }
    // seeded geometries at random scales, and with subnormal coordinates
    for _ in 0..300 {
        let k = i32::try_from(rng.int(-1074, 1018)).unwrap_or(0);
        let b = rpt(&mut rng, 1.0);
        angle(m, sc(add(b, rpt(&mut rng, 1.0)), k), sc(b, k), sc(add(b, rpt(&mut rng, 1.0)), k));
    }
    for _ in 0..100 {
        let a = [binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023)];
        let c = [binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023), 0.0];
        angle(m, a, o, c);
    }
    for _ in 0..100 {
        let a = [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)];
        let b = [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)];
        let c = [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)];
        angle(m, a, b, c);
    }

    // ── dihedral_angle ──
    let b0: P3 = [1.0, 0.0, 0.0];
    let quads: [[P3; 4]; 5] = [
        [[1.0, 1.0, 0.0], o, [1.5, 0.0, 0.0], [2.5, 0.0, 1.0]], // gauche (Rust's bench)
        [o, b0, [2.0, 1.0, 0.0], [3.0, 1.0, 0.0]],              // trans
        [o, b0, [2.0, -1.0, 0.0], [3.0, -1.0, 0.0]],            // trans, mirrored (atan2's −π)
        [[0.0, 1.0, 0.0], o, b0, [1.0, 1.0, 0.0]],              // cis
        [b0, o, y1, [0.0, 1.0, 1.0]],                           // −90°
    ];
    for &k in ks.iter().chain(&ks_top[..3]) {
        for q in &quads {
            dihedral(m, sc(q[0], k), sc(q[1], k), sc(q[2], k), sc(q[3], k));
        }
    }
    // trans and cis with every sign of their zero coordinates at extreme scales (Rust's −π → π)
    for &k in &[-1000, -600, -400, 400, 600, 1000] {
        for mask in (0..128_u32).step_by(5) {
            let z = |j: u32| if mask >> j & 1 == 1 { -0.0 } else { 0.0 };
            for &yc in &[1.0, -1.0] {
                dihedral(m, sc([z(0), z(1), z(2)], k), sc([1.0, z(3), z(4)], k), sc([2.0, yc, z(5)], k), sc([3.0, yc, z(6)], k));
            }
            dihedral(m, sc([z(0), 1.0, z(1)], k), sc([z(2), z(3), z(4)], k), sc([1.0, z(5), z(6)], k), sc([1.0, 1.0, 0.0], k));
        }
    }
    // tiny dihedrals with results near and in the subnormal range: A = (0, 1, 0), B = 0, C = (1, 0, 0), D = (1, 1, t),
    // all scaled by 2^s (φ ≈ t; Rust's b₁ · n₂ underflows)
    for &s in &[-500, -400, -300, -200, 200, 300, 400] {
        for &t in &[-1000, -1022, -1030, -1050, -1070, -1074] {
            for &mant in &[1.0, 3.0, -5.0] {
                let tv = sc2(mant, t);
                dihedral(m, sc(y1, s), o, sc(b0, s), [sc2(1.0, s), sc2(1.0, s), sc2(tv, s)]);
            }
        }
    }
    // the collinearity tolerance to the ulp at extreme scales (A = (−L, h, 0), B = 0, C = (L, 0, 0), scaled by 2^k), h
    // found by bisection where the port's exact test turns and where Rust's does (|n₁|² and its threshold subnormal for
    // k near −240 … −255, overflowing near 272 … 280), the rows around both and between them; B-C-D bent instead too
    for &k in &[-900, -700, -500, -300, -255, -250, -245, -240, 272, 274, 276, 280, 300, 500, 700, 900] {
        for _ in 0..3 {
            let len = rng.log_uniform(-1.0, 1.0);
            let c = sc([len, 0.0, 0.0], k);
            let d1 = sc([2.0 * len, 0.0, len], k);
            let at = |h: f64| sc([-len, h, 0.0], k);
            let (lo, hi) = (len * 1e-10 * 0.25, len * 1e-10 * 4.0);
            let port = crossing(|h| collinear_ux(at(h), o, c, COLLINEAR_SIN_SQ_TOL), lo, hi);
            let rust = crossing(|h| dihedral_parts(at(h), o, c, d1, COLLINEAR_SIN_SQ_TOL).is_none(), lo, hi);
            for h in around_crossings(port, rust) {
                dihedral(m, at(h), o, c, d1);
                // B-C-D bent instead: A = (−L, 0, L) and D = C + (L, h, 0)
                dihedral(m, sc([-len, 0.0, len], k), o, c, sc([2.0 * len, h, 0.0], k));
            }
        }
    }
    // the collinearity tolerance with b₁ tiny and b₂ huge: Rust's threshold (10⁻²⁰·|b₁|²)·|b₂|² rounds 10⁻²⁰·|b₁|² at the
    // subnormal grid (a few bits) before |b₂|² scales it back up, so its turn is far from the exact one, on either side:
    // between the two crossings Rust says Ok where the port says Err, or the reverse. A = (−L₁, h, 0), B = 0,
    // C = (L₂, 0, 0), D = (2L₂, 0, L₂); and the mirror image for B-C-D (A = (−L₁, 0, L₁), D = C + (L₁, h, 0))
    for &(k1, k2) in &[(-500, 300), (-505, 250), (-510, 200), (-515, 150), (-520, 100), (-498, 400)] {
        for _ in 0..4 {
            let l1 = sc2(1.0 + rng.unit(), k1);
            let l2 = sc2(1.0 + rng.unit(), k2);
            let c = [l2, 0.0, 0.0];
            let d1 = [2.0 * l2, 0.0, l2];
            let (lo, hi) = (l1 * 1e-10 * 0.25, l1 * 1e-10 * 4.0);
            let at = |h: f64| [-l1, h, 0.0];
            let port = crossing(|h| collinear_ux(at(h), o, c, COLLINEAR_SIN_SQ_TOL), lo, hi);
            let rust = crossing(|h| dihedral_parts(at(h), o, c, d1, COLLINEAR_SIN_SQ_TOL).is_none(), lo, hi);
            for h in around_crossings(port, rust) {
                dihedral(m, at(h), o, c, d1);
            }
            // B-C-D: b₂ = C − B tiny, b₃ = D − C huge: A = (−L₂, 0, L₂), B = 0, C = (L₁, 0, 0), D = (L₁ + L₂, h·L₂/L₁, 0)
            let c2 = [l1, 0.0, 0.0];
            let a2 = [-l2, 0.0, l2];
            let at2 = |h: f64| [l1 + l2, h * (l2 / l1), 0.0];
            let port2 = crossing(|h| collinear_ux(o, c2, at2(h), COLLINEAR_SIN_SQ_TOL), lo, hi);
            let rust2 = crossing(|h| dihedral_parts(a2, o, c2, at2(h), COLLINEAR_SIN_SQ_TOL).is_none(), lo, hi);
            for h in around_crossings(port2, rust2) {
                dihedral(m, a2, o, c2, at2(h));
            }
        }
    }
    // mixed magnitudes, overflowing differences, subnormal coordinates, seeded scales
    for q in [
        [[big, 0.0, 0.0], [-big, 0.0, 0.0], [-big, big, 0.0], [big, big, big]],
        [[big, big, 0.0], o, [big, 0.0, 0.0], [big, 0.0, big]],
        [[pow2(600), pow2(600), 0.0], o, [pow2(-600), 0.0, 0.0], [pow2(-600), 0.0, pow2(-600)]],
        [[1.0, 1.0, 0.0], o, [pow2(-700), 0.0, 0.0], [1.0, 0.0, 1.0]],
        [[pow2(-1074), pow2(-1074), 0.0], o, [pow2(-1073), 0.0, 0.0], [pow2(-1072), 0.0, pow2(-1074)]],
        [[1e300, 1e300, 0.0], o, [1e-300, 0.0, 0.0], [1e300, 0.0, 1e300]],
    ] {
        dihedral(m, q[0], q[1], q[2], q[3]);
        dihedral(m, q[3], q[2], q[1], q[0]);
    }
    for _ in 0..300 {
        let k = i32::try_from(rng.int(-1074, 1018)).unwrap_or(0);
        let b = rpt(&mut rng, 1.0);
        let c = add(b, rpt(&mut rng, 1.0));
        dihedral(m, sc(add(b, rpt(&mut rng, 1.0)), k), sc(b, k), sc(c, k), sc(add(c, rpt(&mut rng, 1.0)), k));
    }
    for _ in 0..100 {
        let p: [P3; 4] = std::array::from_fn(|_| [binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023), binade(&mut rng, -1074, -1023)]);
        dihedral(m, p[0], p[1], p[2], p[3]);
    }
    for _ in 0..100 {
        let p: [P3; 4] = std::array::from_fn(|_| [binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020), binade(&mut rng, -1074, 1020)]);
        dihedral(m, p[0], p[1], p[2], p[3]);
    }
}
