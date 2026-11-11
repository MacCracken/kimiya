//! The library-wide overflow audit (the user's decisions of 2026-11-09 and 2026-11-10): the port's values of the
//! functions each batch fixes, independent of the port's code, the noted rows, and the edge rows each batch appends.
//!
//! **B1** (ADR 0002 #27, #28): thermo's `heat_transfer` and `specific_heat_from_calorimetry`, phase's
//! `boiling_point_elevation`, `freezing_point_depression` and `osmotic_pressure` (#27), and kinetics'
//! `zero_order_half_life`, `steady_state_intermediate`, `pre_equilibrium_rate` (#27's treatment) and
//! `second_order_concentration`, `michaelis_menten` (#20's) under #28:
//! - a formula that is one product or quotient of inputs and constants is **rounded once**: the exact rational with each
//!   constant as its double and each integer literal of Rust's formula exact, by [`exact::cr_monomial`] (Rust's own
//!   expression for a non-finite input of the two unvalidated functions);
//! - a formula with a sum inside is **Rust's sequence evaluated with an unbounded exponent**, every operation rounded to
//!   53 bits and the final division once into binary64 ([`exact::Ux`]).
//!
//! Every row of these functions goes through [`noted`]: where the port's value differs from the oracle's (by its bits;
//! two NaNs are the same) the row is noted `[ADR 0002 #n: the port returns <how> = <x> (f:<bits>)]`. Rust's validation is
//! replicated here and must agree with the oracle's (asserted on every row). `validate.py` proves every note again,
//! with `Fraction` and MPFR.

use std::hint::black_box as bb;

use crate::exact::{Ux, cr_monomial};
use crate::grids::{Rng, pow2};
use crate::m_reaction::res;
use crate::tsv::{Class, Module, Tok, f};
use kimiya::element::GAS_CONSTANT;
use kimiya::error::KimiyaError;
use kimiya::kinetics::{
    michaelis_menten, pre_equilibrium_rate, second_order_concentration, steady_state_intermediate,
    zero_order_half_life,
};
use kimiya::phase::{boiling_point_elevation, freezing_point_depression, osmotic_pressure};
use kimiya::thermo::{WATER_SPECIFIC_HEAT, heat_transfer, specific_heat_from_calorimetry};

/// How the port obtains its value under #20's semantics (the note's `<how>`).
pub const UX_HOW: &str = "Rust's sequence evaluated with an unbounded exponent (the final division rounded once)";

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

fn note(adr: u32, how: &str, x: f64) -> String {
    format!("ADR 0002 #{adr}: the port returns {how} = {x:?} (f:{:016x})", x.to_bits())
}

/// A class-B row of an f64-valued function, noted where the port's value differs from the oracle's.
pub fn noted_f(m: &mut Module, name: &str, ins: &[Tok], oracle: f64, port: f64, adr: u32, how: &str) {
    if same(oracle, port) {
        m.row(name, Class::B, ins, &[f(oracle)]);
    } else {
        m.row_noted(name, Class::B, ins, &[f(oracle)], &note(adr, how, port));
    }
}

/// A class-B row of a `Result<f64>` function: the port's value (`None`: its replica of Rust's validation refuses, which
/// must be exactly where the oracle errs), noted where it differs from the oracle's.
pub fn noted(
    m: &mut Module,
    name: &str,
    ins: &[Tok],
    oracle: Result<f64, KimiyaError>,
    port: Option<f64>,
    adr: u32,
    how: &str,
) {
    assert_eq!(oracle.is_ok(), port.is_some(), "{name}: the port's validation differs from the oracle's: {ins:?}");
    if let (Ok(v), Some(p)) = (&oracle, port)
        && !same(*v, p)
    {
        m.row_noted(name, Class::B, ins, &[f(*v)], &note(adr, how, p));
        return;
    }
    m.row(name, Class::B, ins, &[res(oracle)]);
}

fn pos(x: f64) -> bool {
    x > 0.0 && x.is_finite()
}

fn nonneg(x: f64) -> bool {
    x >= 0.0 && x.is_finite()
}

/// ADR 0002 #27: `heat_transfer` and the two colligative functions, x·y·z rounded once (Rust's (x·y)·z for a
/// non-finite input).
pub fn port_mul3(x: f64, y: f64, z: f64) -> f64 {
    if x.is_finite() && y.is_finite() && z.is_finite() { cr_monomial(&[x, y, z], &[]) } else { x * y * z }
}

/// ADR 0002 #27: `specific_heat_from_calorimetry`, q/(m·ΔT) rounded once, after Rust's checks.
pub fn port_specific_heat(e: f64, m: f64, dt: f64) -> Option<f64> {
    (e.is_finite() && pos(m) && dt.is_finite() && dt.abs() >= 1e-12).then(|| cr_monomial(&[e], &[m, dt]))
}

/// ADR 0002 #27: `osmotic_pressure`, i·M·R·1000·T rounded once (R the double GAS_CONSTANT, 1000 exact).
pub fn port_osmotic(i: f64, c: f64, t: f64) -> Option<f64> {
    (pos(t) && pos(i) && nonneg(c)).then(|| cr_monomial(&[i, c, GAS_CONSTANT, 1000.0, t], &[]))
}

/// ADR 0002 #28: `zero_order_half_life`, [A]₀/(2k) rounded once.
pub fn port_zero_half(c0: f64, k: f64) -> Option<f64> {
    (pos(k) && pos(c0)).then(|| cr_monomial(&[c0], &[2.0, k]))
}

/// ADR 0002 #28: `steady_state_intermediate`, k₁·[A]/k₂ rounded once.
pub fn port_steady(k1: f64, a: f64, k2: f64) -> Option<f64> {
    (pos(k2) && nonneg(k1) && nonneg(a)).then(|| cr_monomial(&[k1, a], &[k2]))
}

/// ADR 0002 #28: `pre_equilibrium_rate`, k₂·k₁·[A]/k₋₁ rounded once.
pub fn port_pre(kf: f64, kr: f64, ks: f64, a: f64) -> Option<f64> {
    (pos(kr) && nonneg(kf) && nonneg(ks) && nonneg(a)).then(|| cr_monomial(&[ks, kf, a], &[kr]))
}

/// ADR 0002 #28: `second_order_concentration`, Rust's 1/(1/[A]₀ + k·t) with an unbounded exponent, the final division
/// rounded once.
pub fn port_second(c0: f64, k: f64, t: f64) -> Option<f64> {
    (pos(c0) && nonneg(k) && nonneg(t)).then(|| {
        let one = Ux::from_f64(1.0);
        let s = one.div(Ux::from_f64(c0)).add(Ux::from_f64(k).mul(Ux::from_f64(t)));
        one.div_f64(s)
    })
}

/// ADR 0002 #28: `michaelis_menten`, Rust's (Vmax·[S])/(Km + [S]) with an unbounded exponent, the final division
/// rounded once.
pub fn port_mm(v: f64, km: f64, s: f64) -> Option<f64> {
    (pos(v) && pos(km) && nonneg(s)).then(|| {
        let p = Ux::from_f64(v).mul(Ux::from_f64(s));
        p.div_f64(Ux::from_f64(km).add(Ux::from_f64(s)))
    })
}

// ── the rows ──────────────────────────────────────────────────────────────────────────────────────────────────────

pub fn heat(m: &mut Module, ma: f64, c: f64, dt: f64) {
    let q = heat_transfer(bb(ma), bb(c), bb(dt));
    let how = "the correctly rounded m·c·ΔT";
    noted_f(m, "thermo::heat_transfer", &[f(ma), f(c), f(dt)], q, port_mul3(ma, c, dt), 27, how);
}

pub fn calorimetry(m: &mut Module, e: f64, ma: f64, dt: f64) {
    let out = specific_heat_from_calorimetry(bb(e), bb(ma), bb(dt));
    let (name, how) = ("thermo::specific_heat_from_calorimetry", "the correctly rounded q/(m·ΔT)");
    noted(m, name, &[f(e), f(ma), f(dt)], out, port_specific_heat(e, ma, dt), 27, how);
}

pub fn colligative(m: &mut Module, k: f64, mol: f64, i: f64) {
    let how = "the correctly rounded K·m·i";
    let bp = boiling_point_elevation(bb(k), bb(mol), bb(i));
    noted_f(m, "phase::boiling_point_elevation", &[f(k), f(mol), f(i)], bp, port_mul3(k, mol, i), 27, how);
    let fp = freezing_point_depression(bb(k), bb(mol), bb(i));
    noted_f(m, "phase::freezing_point_depression", &[f(k), f(mol), f(i)], fp, port_mul3(k, mol, i), 27, how);
}

pub fn osmotic(m: &mut Module, i: f64, c: f64, t: f64) {
    let out = osmotic_pressure(bb(i), bb(c), bb(t));
    let how = "the correctly rounded i·M·R·1000·T";
    noted(m, "phase::osmotic_pressure", &[f(i), f(c), f(t)], out, port_osmotic(i, c, t), 27, how);
}

pub fn second(m: &mut Module, c0: f64, k: f64, t: f64) {
    let out = second_order_concentration(bb(c0), bb(k), bb(t));
    noted(m, "kinetics::second_order_concentration", &[f(c0), f(k), f(t)], out, port_second(c0, k, t), 28, UX_HOW);
}

pub fn zero_half(m: &mut Module, c0: f64, k: f64) {
    let out = zero_order_half_life(bb(c0), bb(k));
    let how = "the correctly rounded [A]₀/(2k)";
    noted(m, "kinetics::zero_order_half_life", &[f(c0), f(k)], out, port_zero_half(c0, k), 28, how);
}

pub fn steady(m: &mut Module, k1: f64, a: f64, k2: f64) {
    let out = steady_state_intermediate(bb(k1), bb(a), bb(k2));
    let how = "the correctly rounded k₁·[A]/k₂";
    noted(m, "kinetics::steady_state_intermediate", &[f(k1), f(a), f(k2)], out, port_steady(k1, a, k2), 28, how);
}

pub fn pre(m: &mut Module, kf: f64, kr: f64, ks: f64, a: f64) {
    let out = pre_equilibrium_rate(bb(kf), bb(kr), bb(ks), bb(a));
    let how = "the correctly rounded k₂·k₁·[A]/k₋₁";
    let ins = [f(kf), f(kr), f(ks), f(a)];
    noted(m, "kinetics::pre_equilibrium_rate", &ins, out, port_pre(kf, kr, ks, a), 28, how);
}

pub fn mm(m: &mut Module, v: f64, km: f64, s: f64) {
    let out = michaelis_menten(bb(v), bb(km), bb(s));
    noted(m, "kinetics::michaelis_menten", &[f(v), f(km), f(s)], out, port_mm(v, km, s), 28, UX_HOW);
}

// ── B1's edge rows (appended to each module with their own seed, so the rows before keep theirs) ─────────────────

const DMAX: f64 = f64::MAX;

fn sub(k: u64) -> f64 {
    f64::from_bits(k)
}

/// An odd integer y in [lo, hi) with 3y of 54 bits and odd (a midpoint at 53 bits), as a double.
fn odd3(rng: &mut Rng) -> f64 {
    let lo = ((1_i64 << 53) / 3 + 1) / 2;
    let hi = ((1_i64 << 54) / 3 - 1) / 2;
    #[allow(clippy::cast_precision_loss)]
    let y = (rng.int(lo, hi) * 2 + 1) as f64; // below 2⁵³: exact
    y
}

fn scale(rng: &mut Rng, lo: i64, hi: i64) -> f64 {
    pow2(i32::try_from(rng.int(lo, hi)).unwrap_or(0))
}

/// thermo's B1 rows: heat_transfer's partial product m·c past DBL_MAX or below 2⁻¹⁰²² with ΔT bringing it back, the true
/// overflow and underflow ends, 0 against an overflowing m·c (Rust's NaN, the exact ±0), the signs, midpoints of the
/// 53-bit grid (3·y with 3y of 54 bits, at scales across the range) and of the subnormal grid (q·2⁻¹⁰⁷⁴ × ½); the
/// calorimetry's m·ΔT over- and underflowing while q/(m·ΔT) is representable, and its subnormal results and ties.
pub fn emit_thermo(m: &mut Module) {
    let mut rng = Rng::new(0x2f27);
    for (ma, c, dt) in [
        (1e200, 1e200, 1e-200),
        (1e-200, 1e-200, 1e200),
        (-1e200, 1e200, 1e-200),
        (1e300, 1e300, 0.0),
        (1e300, -1e300, 0.0),
        (1e300, 1e300, -0.0),
        (DMAX, 2.0, 0.25),
        (DMAX, DMAX, sub(1)),
        (sub(1), sub(1), DMAX),
        (sub(1), 0.5, 1.0),
        (sub(3), 0.5, 1.0),
        (sub(5), 0.5, -1.0),
        (sub(1), 0.25, 3.0),
        (f64::MIN_POSITIVE, 0.5, 1.0),
        (DMAX, 1.0, 1.0),
        (DMAX, 1.0000000000000002, 1.0),
        (1e308, 10.0, 1.0),
        (1e-308, 0.1, 1.0),
        (1000.0, WATER_SPECIFIC_HEAT, 1e-320),
        (3.0, 3_002_399_751_580_331.0, 1.0),
        (f64::INFINITY, 1e-300, 1e-300),
        (1e-300, 1e-300, f64::INFINITY),
        (f64::NAN, 0.0, 1.0),
    ] {
        heat(m, ma, c, dt);
    }
    for _ in 0..40 {
        // m·c past DBL_MAX, ΔT bringing it back; m·c below 2⁻¹⁰²², ΔT bringing it up
        let (a, b) = (rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        heat(m, a, b, rng.signed_log_uniform(-300.0, -1.0));
        let (a, b) = (rng.signed_log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0));
        heat(m, a, b, rng.log_uniform(1.0, 300.0));
        // the result past DBL_MAX or below 2⁻¹⁰⁷⁴ (true overflow and underflow), and subnormal
        heat(m, rng.log_uniform(150.0, 308.0), rng.log_uniform(150.0, 308.0), rng.log_uniform(0.0, 100.0));
        heat(m, rng.log_uniform(-320.0, -100.0), rng.log_uniform(-300.0, -100.0), rng.log_uniform(-100.0, 100.0));
        heat(m, rng.log_uniform(-300.0, -150.0), rng.log_uniform(-200.0, -150.0), rng.log_uniform(0.0, 50.0));
    }
    for _ in 0..60 {
        let y = odd3(&mut rng);
        let sc = scale(&mut rng, -900, 900);
        heat(m, 3.0, y * sc, 1.0);
        heat(m, -1.5, y, scale(&mut rng, -60, 60));
        // q·2⁻¹⁰⁷⁵ for an odd q (a tie of the subnormal grid): q·2^(a − 1074) times 2^-(a + 1)
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        let q = sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a);
        heat(m, q, pow2(-a - 1), if rng.int(0, 1) == 0 { 1.0 } else { -1.0 });
    }
    for (e, ma, dt) in [
        (1e300, 1e200, 1e200),
        (-1e300, 1e200, 1e200),
        (1e-300, 1e-300, 1e-12),
        (1e-300, 1e-200, 1e-100),
        (sub(1), 1.0, 1.0),
        (sub(3), 2.0, 1.0),
        (sub(5), 2.0, -1.0),
        (sub(1), 2.0, 1.0),
        (sub(3), 4.0, 1.0),
        (f64::MIN_POSITIVE, 3.0, 1.0),
        (DMAX, 0.5, 1.0),
        (DMAX, f64::MIN_POSITIVE, 1e-12),
        (1.0, DMAX, DMAX),
        (-0.0, DMAX, DMAX),
    ] {
        calorimetry(m, e, ma, dt);
    }
    for _ in 0..40 {
        // m·ΔT past DBL_MAX while q is large; m·ΔT subnormal while q is tiny
        calorimetry(m, rng.signed_log_uniform(200.0, 308.0), rng.log_uniform(154.0, 308.0), rng.signed_log_uniform(154.0, 300.0));
        calorimetry(m, rng.signed_log_uniform(-320.0, -250.0), rng.log_uniform(-320.0, -300.0), rng.signed_log_uniform(-12.0, -5.0));
        calorimetry(m, rng.signed_log_uniform(-320.0, -200.0), rng.log_uniform(-1.0, 100.0), rng.signed_log_uniform(-1.0, 100.0));
    }
}

/// phase's B1 rows: the colligative K·m·i over the same families as heat_transfer, and osmotic_pressure's (i·M)·R·…
/// with i·M past DBL_MAX or below 2⁻¹⁰²² while T brings it back, subnormal results and the overflow end.
pub fn emit_phase(m: &mut Module) {
    let mut rng = Rng::new(0x3127);
    for (k, mol, i) in [
        (1e200, 1e200, 1e-200),
        (1e-200, 1e-200, 1e200),
        (1e300, 1e300, 0.0),
        (-1e300, 1e300, 0.0),
        (DMAX, 2.0, 0.25),
        (sub(1), sub(1), DMAX),
        (sub(3), 0.5, 1.0),
        (sub(1), 0.5, 1.0),
        (0.512, 3_002_399_751_580_331.0, 3.0),
        (f64::INFINITY, 0.0, 1.0),
    ] {
        colligative(m, k, mol, i);
    }
    for _ in 0..40 {
        let (a, b) = (rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        colligative(m, a, b, rng.signed_log_uniform(-300.0, -1.0));
        let (a, b) = (rng.signed_log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0));
        colligative(m, a, b, rng.log_uniform(1.0, 300.0));
        colligative(m, rng.log_uniform(-300.0, -150.0), rng.log_uniform(-200.0, -150.0), rng.log_uniform(0.0, 50.0));
    }
    for _ in 0..40 {
        let y = odd3(&mut rng);
        colligative(m, 3.0, y * scale(&mut rng, -900, 900), 1.0);
    }
    for (i, c, t) in [
        (1.0, sub(1), 298.15),
        (1e300, 1e10, 1e-20),
        (1e200, 1e200, 1e-100),
        (1.0, 1e-310, 1e300),
        (1.0, 1e-320, 1e10),
        (DMAX, 1.0, 1e-4),
        (DMAX, DMAX, sub(1)),
        (1.0, 1e300, 1e10),
    ] {
        osmotic(m, i, c, t);
    }
    for _ in 0..60 {
        osmotic(m, rng.log_uniform(150.0, 308.0), rng.log_uniform(150.0, 308.0), rng.log_uniform(-300.0, -10.0));
        osmotic(m, rng.log_uniform(-320.0, -150.0), rng.log_uniform(-320.0, -10.0), rng.log_uniform(10.0, 300.0));
        osmotic(m, rng.log_uniform(-320.0, -300.0), rng.log_uniform(-10.0, 1.0), rng.log_uniform(0.0, 3.0));
    }
}

/// kinetics' B1 rows: second_order_concentration's 1/[A]₀ past DBL_MAX or below 2⁻¹⁰²², k·t over- or underflowing, the
/// sum past DBL_MAX, the final division's subnormal results and its tie 2⁻¹⁰⁷⁵; zero_order_half_life past 2·k's overflow
/// (k ≥ 2¹⁰²³); steady_state_intermediate's and pre_equilibrium_rate's partial products past DBL_MAX or below 2⁻¹⁰²²,
/// their midpoints and subnormal ties, 0 against an overflowing product; michaelis_menten's Vmax·[S] and Km + [S]
/// overflowing (∞/∞) or Vmax·[S] underflowing.
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_kinetics(m: &mut Module) {
    let mut rng = Rng::new(0x3c28);
    for (c0, k, t) in [
        (sub(1), 1.0, 1.0),
        (sub(3), 1e-300, 1e-10),
        (f64::MIN_POSITIVE, 0.0, 0.0),
        (1e308, 0.0, 1.0),
        (1e308, 1e-320, 1.0),
        (DMAX, 1e-300, 1e-300),
        (1.0, 1e300, 1e300),
        (1.0, DMAX, DMAX),
        (1.0, 1e-200, 1e-200),
        (1.0, pow2(1000), pow2(75)),
        (1.0, pow2(1000), pow2(74)),
        (1e-10, 1e300, 1e10),
        (1e300, 1e-300, 1e-300),
    ] {
        second(m, c0, k, t);
    }
    for _ in 0..60 {
        second(m, rng.log_uniform(-324.0, -308.0), rng.log_uniform(-10.0, 10.0), rng.log_uniform(-10.0, 10.0));
        second(m, rng.log_uniform(307.0, 308.0), rng.log_uniform(-320.0, -290.0), rng.log_uniform(-10.0, 10.0));
        second(m, rng.log_uniform(-10.0, 10.0), rng.log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        second(m, rng.log_uniform(-10.0, 10.0), rng.log_uniform(-320.0, -154.0), rng.log_uniform(-200.0, -100.0));
        second(m, rng.log_uniform(-300.0, -200.0), rng.log_uniform(250.0, 308.0), rng.log_uniform(1.0, 60.0));
    }
    for (c0, k) in [
        (1.0, pow2(1023)),
        (DMAX, pow2(1023)),
        (sub(1), pow2(1023)),
        (sub(3), pow2(1023)),
        (f64::MIN_POSITIVE, pow2(1023)),
        (1.0, f64::from_bits(0x7FDF_FFFF_FFFF_FFFF)),
        (3.0, pow2(1023) * 1.5),
    ] {
        zero_half(m, c0, k);
    }
    for _ in 0..60 {
        let k = pow2(1023) * rng.uniform(1.0, 2.0);
        zero_half(m, rng.log_uniform(-324.0, 308.0), k);
        zero_half(m, rng.log_uniform(280.0, 308.0), k);
    }
    for (k1, a, k2) in [
        (1e300, 1e300, 1e300),
        (1e-300, 1e-300, 1e-300),
        (DMAX, 2.0, 4.0),
        (sub(1), sub(1), sub(1)),
        (sub(3), 1.0, 2.0),
        (sub(1), 1.0, 2.0),
        (3.0, 3_002_399_751_580_331.0, 1.0),
        (0.0, DMAX, sub(1)),
        (-0.0, 1e300, 1e-300),
    ] {
        steady(m, k1, a, k2);
    }
    for _ in 0..60 {
        steady(m, rng.log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0), rng.log_uniform(1.0, 308.0));
        steady(m, rng.log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0), rng.log_uniform(-320.0, -1.0));
        let y = odd3(&mut rng);
        steady(m, y * scale(&mut rng, -900, 900), 3.0, scale(&mut rng, -60, 60));
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        steady(m, sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a), 1.0, pow2(a + 1));
    }
    for (kf, kr, ks, a) in [
        (1e300, 1.0, 1e300, -0.0),
        (1e300, 1e-300, 0.0, 1.0),
        (1e300, 1e-300, 1e-300, 1e-300),
        (1e-300, 1e300, 1e300, 1e300),
        (sub(1), 1.0, 1.0, 1.0),
        (sub(3), 2.0, 1.0, 1.0),
        (3.0, 1.0, 3_002_399_751_580_331.0, 1.0),
        (DMAX, DMAX, DMAX, 1.0),
    ] {
        pre(m, kf, kr, ks, a);
    }
    for _ in 0..60 {
        pre(m, rng.log_uniform(150.0, 308.0), rng.log_uniform(-308.0, -150.0), rng.log_uniform(-308.0, -10.0), rng.log_uniform(-300.0, 0.0));
        pre(m, rng.log_uniform(-320.0, -150.0), rng.log_uniform(150.0, 308.0), rng.log_uniform(10.0, 308.0), rng.log_uniform(1.0, 300.0));
        let y = odd3(&mut rng);
        pre(m, 3.0, scale(&mut rng, -60, 60), y * scale(&mut rng, -800, 800), 1.0);
    }
    for (v, km, s) in [
        (DMAX, DMAX, DMAX),
        (1e-300, 1e-300, 1e-300),
        (sub(1), 1.0, 1.0),
        (sub(3), 1.0, 1.0),
        (1e300, 1e300, 1e300),
        (1.0, 1e-300, sub(1)),
        (DMAX, 1e-300, 1e-300),
    ] {
        mm(m, v, km, s);
    }
    for _ in 0..60 {
        mm(m, rng.log_uniform(1.0, 308.0), rng.log_uniform(300.0, 308.0), rng.log_uniform(300.0, 308.0));
        mm(m, rng.log_uniform(150.0, 308.0), rng.log_uniform(-300.0, 300.0), rng.log_uniform(154.0, 308.0));
        mm(m, rng.log_uniform(-320.0, -150.0), rng.log_uniform(-320.0, 0.0), rng.log_uniform(-320.0, -150.0));
    }
}
