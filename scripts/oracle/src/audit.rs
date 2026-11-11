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
//!
//! **B2** (ADR 0002 #29, #30): solution's `dilution` and `solution_conductivity`, electrochemistry's
//! `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` (#29, rounded once) and `nernst_potential`
//! (#29, #20's semantics: Rust's E° − ((R·T)/(n·F))·ln Q with an unbounded exponent and the final subtraction rounded
//! once, ln Q the correctly rounded ln of the input, as the port's), and spectroscopy's `absorbance`,
//! `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia` (#30,
//! rounded once; π as its double each time Rust uses it, 8, 1000 and J + 1 exact) and `vibrational_frequency` (#30, the
//! correctly rounded √ of the exact k/μ, [`cr_sqrt_ratio`]).
//!
//! Every row of these functions goes through [`noted`]: where the port's value differs from the oracle's (by its bits;
//! two NaNs are the same) the row is noted `[ADR 0002 #n: the port returns <how> = <x> (f:<bits>)]`. Rust's validation is
//! replicated here and must agree with the oracle's (asserted on every row). `validate.py` proves every note again,
//! with `Fraction` and MPFR.

use std::hint::black_box as bb;

use crate::exact::{Ux, cr_monomial, cr_sqrt_ratio};
use crate::grids::{Rng, pow2};
use crate::m_reaction::res;
use crate::tsv::{Class, Module, Tok, f, i};
use kimiya::electrochemistry::{
    FARADAY, faraday_charge_required, faraday_mass_deposited, faraday_mass_ratio, nernst_potential,
};
use kimiya::element::GAS_CONSTANT;
use kimiya::error::KimiyaError;
use kimiya::kinetics::{
    michaelis_menten, pre_equilibrium_rate, second_order_concentration, steady_state_intermediate,
    zero_order_half_life,
};
use kimiya::phase::{boiling_point_elevation, freezing_point_depression, osmotic_pressure};
use kimiya::solution::{dilution, solution_conductivity};
use kimiya::spectroscopy::{
    PLANCK, absorbance, concentration_from_absorbance, diatomic_moment_of_inertia, rigid_rotor_energy,
    rotational_constant, vibrational_frequency,
};
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
    noted_cls(m, name, Class::B, ins, oracle, port, adr, how);
}

/// [`noted`] for a row of class `cls` (B2's `nernst_potential` is class C: one ln).
#[allow(clippy::too_many_arguments)] // noted's, and the class
pub fn noted_cls(
    m: &mut Module,
    name: &str,
    cls: Class,
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
        m.row_noted(name, cls, ins, &[f(*v)], &note(adr, how, p));
        return;
    }
    m.row(name, cls, ins, &[res(oracle)]);
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


// ── The overflow audit's B2 (ADR 0002 #29, #30) ─────────────────────────────────────────────────────────────────

/// How the port obtains `nernst_potential`'s value under #20's semantics (the note's `<how>`).
pub const UX_HOW_SUB: &str =
    "Rust's sequence evaluated with an unbounded exponent (the final subtraction rounded once)";

/// ADR 0002 #29: `dilution`, M₁·V₁/V₂ rounded once, after Rust's checks.
pub fn port_dilution(m1: f64, v1: f64, v2: f64) -> Option<f64> {
    (pos(v2) && m1.is_finite() && v1.is_finite()).then(|| cr_monomial(&[m1, v1], &[v2]))
}

/// ADR 0002 #29: `solution_conductivity`, Λ·c/1000 rounded once (Rust's (Λ·c)/1000 for a non-finite input).
pub fn port_conductivity(l: f64, c: f64) -> f64 {
    if l.is_finite() && c.is_finite() { cr_monomial(&[l, c], &[1000.0]) } else { l * c / 1000.0 }
}

/// ADR 0002 #29: `faraday_mass_deposited`, Q·M/(n·F) rounded once (F as its double).
pub fn port_fmass(q: f64, mm: f64, n: u8) -> Option<f64> {
    (n != 0 && pos(mm) && q.is_finite()).then(|| cr_monomial(&[q, mm], &[f64::from(n), FARADAY]))
}

/// ADR 0002 #29: `faraday_charge_required`, m·n·F/M rounded once.
pub fn port_fcharge(ma: f64, mm: f64, n: u8) -> Option<f64> {
    (n != 0 && pos(mm) && ma.is_finite()).then(|| cr_monomial(&[ma, f64::from(n), FARADAY], &[mm]))
}

/// ADR 0002 #29: `faraday_mass_ratio`, Rust's (m_A·(M_B/n_B))/(M_A/n_A) exactly, m_A·M_B·n_A/(n_B·M_A), rounded once.
pub fn port_fratio(ma: f64, mma: f64, na: u8, mmb: f64, nb: u8) -> Option<f64> {
    (na != 0 && nb != 0 && pos(mma) && pos(mmb) && ma.is_finite())
        .then(|| cr_monomial(&[ma, mmb, f64::from(na)], &[f64::from(nb), mma]))
}

/// ADR 0002 #29: `nernst_potential`, Rust's E° − ((R·T)/(n·F))·ln Q with an unbounded exponent, the final subtraction
/// rounded once into binary64, after Rust's checks. ln Q is the correctly rounded ln of Q (the port's `_km_cr_ln`):
/// glibc's here, which `validate.py` checks against MPFR on every row.
pub fn port_nernst(e0: f64, n: u8, t: f64, q: f64) -> Option<f64> {
    (pos(t) && n != 0 && pos(q) && e0.is_finite()).then(|| {
        let a = Ux::from_f64(GAS_CONSTANT).mul(Ux::from_f64(t));
        let c = a.div(Ux::from_f64(f64::from(n)).mul(Ux::from_f64(FARADAY)));
        let p = c.mul(Ux::from_f64(bb(q).ln()));
        Ux::from_f64(e0).sub_f64(p)
    })
}

/// ADR 0002 #30: `absorbance` and `diatomic_moment_of_inertia`, x·y·z rounded once (Rust's (x·y)·z for a non-finite
/// input).
pub fn port_mul3_30(x: f64, y: f64, z: f64) -> f64 {
    port_mul3(x, y, z)
}

/// ADR 0002 #30: `concentration_from_absorbance`, A/(ε·l) rounded once for a finite A (Rust's A/(ε·l) otherwise).
pub fn port_conc(a: f64, e: f64, l: f64) -> Option<f64> {
    (pos(e) && pos(l)).then(|| if a.is_finite() { cr_monomial(&[a], &[e, l]) } else { a / (e * l) })
}

/// ADR 0002 #30: `vibrational_frequency`, the correctly rounded √ of the exact k/μ.
pub fn port_vib(k: f64, mu: f64) -> Option<f64> {
    (pos(k) && pos(mu)).then(|| cr_sqrt_ratio(k, mu))
}

/// ADR 0002 #30: `rigid_rotor_energy`, h·B·J·(J + 1) rounded once for a finite B (Rust's ((h·B)·J)·(J + 1) otherwise).
pub fn port_rotor(j: u32, b: f64) -> f64 {
    let (jf, j1) = (f64::from(j), f64::from(j) + 1.0);
    if b.is_finite() { cr_monomial(&[PLANCK, b, jf, j1], &[]) } else { PLANCK * b * jf * j1 }
}

/// ADR 0002 #30: `rotational_constant`, h/(8·π·π·I) rounded once (π as its double).
pub fn port_rot_const(i_: f64) -> Option<f64> {
    pos(i_).then(|| cr_monomial(&[PLANCK], &[8.0, std::f64::consts::PI, std::f64::consts::PI, i_]))
}

pub fn dilution_row(m: &mut Module, a: f64, b: f64, c: f64) {
    let out = dilution(bb(a), bb(b), bb(c));
    noted(m, "solution::dilution", &[f(a), f(b), f(c)], out, port_dilution(a, b, c), 29, "the correctly rounded M₁·V₁/V₂");
}

pub fn conductivity(m: &mut Module, a: f64, c: f64) {
    let out = solution_conductivity(bb(a), bb(c));
    let how = "the correctly rounded Λ·c/1000";
    noted_f(m, "solution::solution_conductivity", &[f(a), f(c)], out, port_conductivity(a, c), 29, how);
}

pub fn fmass(m: &mut Module, q: f64, mm: f64, n: u8) {
    let out = faraday_mass_deposited(bb(q), bb(mm), bb(n));
    let (name, how) = ("electrochemistry::faraday_mass_deposited", "the correctly rounded Q·M/(n·F)");
    noted(m, name, &[f(q), f(mm), i(n)], out, port_fmass(q, mm, n), 29, how);
}

pub fn fcharge(m: &mut Module, ma: f64, mm: f64, n: u8) {
    let out = faraday_charge_required(bb(ma), bb(mm), bb(n));
    let (name, how) = ("electrochemistry::faraday_charge_required", "the correctly rounded m·n·F/M");
    noted(m, name, &[f(ma), f(mm), i(n)], out, port_fcharge(ma, mm, n), 29, how);
}

pub fn fratio(m: &mut Module, ma: f64, mma: f64, na: u8, mmb: f64, nb: u8) {
    let out = faraday_mass_ratio(bb(ma), bb(mma), bb(na), bb(mmb), bb(nb));
    let (name, how) = ("electrochemistry::faraday_mass_ratio", "the correctly rounded m_A·M_B·n_A/(n_B·M_A)");
    let ins = [f(ma), f(mma), i(na), f(mmb), i(nb)];
    noted(m, name, &ins, out, port_fratio(ma, mma, na, mmb, nb), 29, how);
}

pub fn nernst(m: &mut Module, e0: f64, n: u8, t: f64, q: f64) {
    let out = nernst_potential(bb(e0), bb(n), bb(t), bb(q));
    let name = "electrochemistry::nernst_potential";
    noted_cls(m, name, Class::C, &[f(e0), i(n), f(t), f(q)], out, port_nernst(e0, n, t, q), 29, UX_HOW_SUB);
}

pub fn absorbance_row(m: &mut Module, e: f64, l: f64, c: f64) {
    let out = absorbance(bb(e), bb(l), bb(c));
    noted_f(m, "spectroscopy::absorbance", &[f(e), f(l), f(c)], out, port_mul3_30(e, l, c), 30, "the correctly rounded ε·l·c");
}

pub fn conc(m: &mut Module, a: f64, e: f64, l: f64) {
    let out = concentration_from_absorbance(bb(a), bb(e), bb(l));
    let (name, how) = ("spectroscopy::concentration_from_absorbance", "the correctly rounded A/(ε·l)");
    noted(m, name, &[f(a), f(e), f(l)], out, port_conc(a, e, l), 30, how);
}

pub fn vib(m: &mut Module, k: f64, mu: f64) {
    let out = vibrational_frequency(bb(k), bb(mu));
    let (name, how) = ("spectroscopy::vibrational_frequency", "the correctly rounded √(k/μ)");
    noted(m, name, &[f(k), f(mu)], out, port_vib(k, mu), 30, how);
}

pub fn rotor(m: &mut Module, j: u32, b: f64) {
    let out = rigid_rotor_energy(bb(j), bb(b));
    let how = "the correctly rounded h·B·J·(J + 1)";
    noted_f(m, "spectroscopy::rigid_rotor_energy", &[i(j), f(b)], out, port_rotor(j, b), 30, how);
}

pub fn rot_const(m: &mut Module, x: f64) {
    let out = rotational_constant(bb(x));
    let (name, how) = ("spectroscopy::rotational_constant", "the correctly rounded h/(8·π·π·I)");
    noted(m, name, &[f(x)], out, port_rot_const(x), 30, how);
}

pub fn moment(m: &mut Module, mu: f64, r: f64) {
    let out = diatomic_moment_of_inertia(bb(mu), bb(r));
    let how = "the correctly rounded μ·r·r";
    noted_f(m, "spectroscopy::diatomic_moment_of_inertia", &[f(mu), f(r)], out, port_mul3_30(mu, r, r), 30, how);
}

// ── B2's edge rows (appended to each module with their own seed, so the rows before keep theirs) ─────────────────

/// An odd integer r with 3r² odd and of 54 bits (a midpoint at 53 bits) and r² exact (below 2⁵³), as a double:
/// `diatomic_moment_of_inertia(3, r)`'s 3·r·r.
fn odd_sq3(rng: &mut Rng) -> f64 {
    // 3r² in [2^53, 2^54) and r² < 2^53: r in [⌈√(2^53/3)⌉, ⌊√(2^54/3)⌋]
    let lo: i64 = 54_794_159;
    let hi: i64 = 77_490_641;
    #[allow(clippy::cast_precision_loss)]
    let r = (rng.int(lo / 2, hi / 2 - 1) * 2 + 1) as f64; // below 2^27: exact
    r
}

/// solution's B2 rows: dilution's M₁·V₁ past DBL_MAX or below 2⁻¹⁰²² with V₂ bringing it back, the true overflow and
/// underflow ends, the signs and signed zeros, 53-bit midpoints at scales across the range and the subnormal grid's
/// ties; solution_conductivity's Λ·c in (DBL_MAX, 1000·DBL_MAX] (Rust's ∞), its subnormal results, midpoints
/// (3000·y/1000 = 3y) and non-finite inputs (Rust's expression).
pub fn emit_solution(m: &mut Module) {
    let mut rng = Rng::new(0x5029);
    for (a, b, c) in [
        (1e200, 1e200, 1e200),
        (1e-200, 1e-200, 1e-200),
        (-1e200, 1e200, 1e200),
        (DMAX, 2.0, 4.0),
        (DMAX, DMAX, DMAX),
        (sub(1), sub(1), sub(1)),
        (sub(3), 1.0, 2.0),
        (sub(1), 1.0, 2.0),
        (sub(5), -1.0, 2.0),
        (3.0, 3_002_399_751_580_331.0, 1.0),
        (0.0, DMAX, sub(1)),
        (-0.0, 1e300, 1e-300),
        (1.0, -1e300, 1e-300),
        (f64::MIN_POSITIVE, 0.5, 1.0),
        (1e300, 1e300, 0.0),
    ] {
        dilution_row(m, a, b, c);
    }
    for _ in 0..50 {
        let (a, b) = (rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        dilution_row(m, a, b, rng.log_uniform(1.0, 308.0));
        let (a, b) = (rng.signed_log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0));
        dilution_row(m, a, b, rng.log_uniform(-320.0, -1.0));
        dilution_row(m, rng.log_uniform(-320.0, -150.0), rng.log_uniform(-200.0, -100.0), rng.log_uniform(0.0, 30.0));
        let y = odd3(&mut rng);
        dilution_row(m, y * scale(&mut rng, -900, 900), 3.0, scale(&mut rng, -60, 60));
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        dilution_row(m, sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a), 1.0, pow2(a + 1));
    }
    for (a, c) in [
        (126.4, DMAX),
        (DMAX, DMAX),
        (1e300, 1e10),
        (-1e300, 1e10),
        (sub(1), 0.5),
        (sub(1000), 1.0),
        (1e-300, 1e-20),
        (3000.0, 3_002_399_751_580_331.0),
        (f64::INFINITY, 0.0),
        (f64::NAN, 1.0),
        (f64::INFINITY, 1e-300),
        (0.0, DMAX),
        (-0.0, 1.0),
    ] {
        conductivity(m, a, c);
    }
    for _ in 0..60 {
        // Λ·c in [2¹⁰²⁴, 2¹⁰³⁴): past DBL_MAX (Rust's ∞), the value finite below 1000·DBL_MAX
        let (l, c) = (rng.signed_log_uniform(154.13, 155.62), rng.log_uniform(154.13, 155.62));
        conductivity(m, l, c);
        conductivity(m, rng.signed_log_uniform(-320.0, -150.0), rng.log_uniform(-170.0, 0.0));
        let y = odd3(&mut rng);
        conductivity(m, 3000.0, y * scale(&mut rng, -900, 900));
    }
}

/// electrochemistry's B2 rows: nernst_potential's R·T past DBL_MAX (T > 2.16·10³⁰⁷ K: Rust's ±∞, or NaN at Q = 1) and
/// its subnormal (R·T)/(n·F) for tiny T (Rust's degraded or zero product), E° against an overflowing or underflowing
/// product, and the true overflow ends; Faraday's partial products past either end, an equivalent weight M/n
/// underflowing, the true ends, and the 53-bit midpoints 3F·y/F and y·3·F/F (3F is exact: F's odd part has 51 bits).
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_electrochemistry(m: &mut Module) {
    let mut rng = Rng::new(0xec29);
    for (e0, n, t, q) in [
        (0.342, 2, DMAX, 0.1),
        (0.342, 2, DMAX, 1.0),
        (-0.762, 2, DMAX, 1e10),
        (0.0, 1, DMAX, DMAX),
        (0.0, 255, DMAX, sub(1)),
        (DMAX, 1, DMAX, sub(1)),
        (-DMAX, 1, DMAX, DMAX),
        (-DMAX, 1, DMAX, sub(1)),
        (1e308, 255, 3e307, 0.5),
        (0.0, 1, 1e-310, 1e-300),
        (0.0, 255, 1e-305, 0.1),
        (1.0, 1, 1e-310, 1e-300),
        (0.0, 1, sub(1), 1e-300),
        (sub(5), 1, 1e-320, 0.5),
        (-0.0, 1, sub(1), 2.0),
        (0.0, 1, sub(1), 1.0),
    ] {
        nernst(m, e0, n, t, q);
    }
    for _ in 0..80 {
        let n = u8::try_from(rng.int(1, 255)).unwrap_or(1);
        let q = rng.log_uniform(-320.0, 308.0);
        let t = rng.log_uniform(307.34, 308.25);
        nernst(m, rng.signed_log_uniform(-3.0, 308.0), n, t, q);
        let t = rng.log_uniform(-324.0, -300.0);
        nernst(m, if rng.int(0, 2) == 0 { 0.0 } else { rng.signed_log_uniform(-324.0, -300.0) }, n, t, q);
        nernst(m, rng.uniform(-3.0, 3.0), n, rng.log_uniform(-324.0, -302.0), q);
    }
    let f3 = 3.0 * FARADAY; // exact
    for (q, mm, n) in [
        (DMAX, 63.546, 2),
        (1e300, 1e10, 1),
        (-1e300, 1e10, 255),
        (sub(1), 1e10, 1),
        (1e-300, 1e-20, 255),
        (-DMAX, 1e10, 255),
        (0.0, DMAX, 1),
        (-0.0, DMAX, 1),
        (DMAX, DMAX, 1),
        (f3, 3_002_399_751_580_331.0, 1),
    ] {
        fmass(m, q, mm, n);
    }
    for (ma, mm, n) in [
        (DMAX, 1.0, 2),
        (DMAX, 1e10, 255),
        (1e300, 1e-10, 1),
        (sub(1), 1e10, 1),
        (1e-300, 1e10, 255),
        (3_002_399_751_580_331.0, FARADAY, 3),
        (-0.0, DMAX, 7),
        (1e-320, 1e-300, 1),
    ] {
        fcharge(m, ma, mm, n);
    }
    for (ma, mma, na, mmb, nb) in [
        (1.0, sub(1), 3, 1e-300, 1),
        (DMAX, 1.0, 1, DMAX, 1),
        (DMAX, DMAX, 1, 2.0, 1),
        (1e-300, 1e300, 255, 1e-300, 1),
        (1e300, sub(1), 1, 1e-300, 255),
        (3_002_399_751_580_331.0, 1.0, 3, 1.0, 1),
        (3_002_399_751_580_331.0, 7.0, 3, 14.0, 2),
        (-0.0, sub(1), 255, DMAX, 1),
        (sub(3), 1.0, 1, 1.0, 2),
    ] {
        fratio(m, ma, mma, na, mmb, nb);
    }
    for _ in 0..60 {
        let n = u8::try_from(rng.int(1, 255)).unwrap_or(1);
        // Q·M past DBL_MAX, n·F bringing it back; Q·M below 2⁻¹⁰²²
        fmass(m, rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0), n);
        fmass(m, rng.signed_log_uniform(-320.0, -150.0), rng.log_uniform(-170.0, 10.0), n);
        let y = odd3(&mut rng);
        fmass(m, f3 * scale(&mut rng, -600, 600), y, 1);
        fcharge(m, rng.signed_log_uniform(290.0, 308.0), rng.log_uniform(1.0, 308.0), n);
        fcharge(m, rng.signed_log_uniform(-324.0, -290.0), rng.log_uniform(-300.0, 0.0), n);
        fcharge(m, y * scale(&mut rng, -600, 600), FARADAY * scale(&mut rng, -60, 60), 3);
        let nb = u8::try_from(rng.int(1, 255)).unwrap_or(1);
        // M_A/n_A subnormal or 0; m·(M_B/n_B) past DBL_MAX
        fratio(m, rng.signed_log_uniform(-300.0, 300.0), rng.log_uniform(-324.0, -305.0), n, rng.log_uniform(-300.0, 0.0), nb);
        fratio(m, rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0), n, rng.log_uniform(154.0, 308.0), nb);
        fratio(m, y * scale(&mut rng, -600, 600), 7.0, 3, 7.0 * f64::from(nb), nb);
    }
}

/// spectroscopy's B2 rows: absorbance's and the moment of inertia's partial products past either end, the true ends, 0
/// against an overflowing product, signs, midpoints (3·y·1, 3·r·r with 3r² of 54 bits) and subnormal ties;
/// concentration_from_absorbance's ε·l past either end and its subnormal ties; vibrational_frequency's k/μ past
/// DBL_MAX or below 2⁻¹⁰⁷⁴ (Rust's ∞ and 0), subnormal operands and results, exact roots, roots within 2⁻¹⁰⁷ of a
/// 53-bit midpoint; rigid_rotor_energy's h·B
/// underflowing (|B| below 3.4·10⁻²⁷⁵) over every J magnitude; rotational_constant's 8π²·I subnormal (I below
/// 2.8·10⁻³¹⁰).
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_spectroscopy(m: &mut Module) {
    let mut rng = Rng::new(0x5c30);
    for (e, l, c) in [
        (1e300, 1e300, 0.0),
        (1e300, -1e300, 0.0),
        (DMAX, 2.0, 0.25),
        (DMAX, DMAX, sub(1)),
        (sub(1), sub(1), DMAX),
        (sub(3), 0.5, 1.0),
        (sub(1), 0.5, 1.0),
        (3.0, 3_002_399_751_580_331.0, 1.0),
        (f64::INFINITY, 1e-300, 1e-300),
        (f64::NAN, 0.0, 1.0),
        (1e-300, 1e-10, 1e-10),
    ] {
        absorbance_row(m, e, l, c);
    }
    for _ in 0..50 {
        let (a, b) = (rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        absorbance_row(m, a, b, rng.signed_log_uniform(-300.0, -1.0));
        let (a, b) = (rng.signed_log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0));
        absorbance_row(m, a, b, rng.log_uniform(1.0, 300.0));
        absorbance_row(m, rng.log_uniform(-300.0, -150.0), rng.log_uniform(-200.0, -150.0), rng.log_uniform(0.0, 50.0));
        let y = odd3(&mut rng);
        absorbance_row(m, 3.0, y * scale(&mut rng, -900, 900), 1.0);
    }
    for (a, e, l) in [
        (sub(3), 1.0, 2.0),
        (sub(1), 1.0, 2.0),
        (sub(5), 2.0, 1.0),
        (1.0, DMAX, DMAX),
        (DMAX, DMAX, 1e-300),
        (1e-300, 1e-300, 1e-300),
        (f64::INFINITY, 1e200, 1e200),
        (f64::NEG_INFINITY, 1.0, 1.0),
        (-0.0, DMAX, DMAX),
        (f64::MIN_POSITIVE, 3.0, 1.0),
    ] {
        conc(m, a, e, l);
    }
    for _ in 0..50 {
        // ε·l past DBL_MAX with A large; ε·l subnormal with A tiny; subnormal results
        conc(m, rng.signed_log_uniform(200.0, 308.0), rng.log_uniform(154.0, 308.0), rng.log_uniform(154.0, 300.0));
        conc(m, rng.signed_log_uniform(-320.0, -250.0), rng.log_uniform(-320.0, -160.0), rng.log_uniform(-160.0, -5.0));
        conc(m, rng.signed_log_uniform(-320.0, -200.0), rng.log_uniform(-1.0, 100.0), rng.log_uniform(-1.0, 100.0));
    }
    let nmin = f64::MIN_POSITIVE;
    for (k, mu) in [
        (500.0, nmin),
        (1e-300, 1e30),
        (DMAX, sub(1)),
        (sub(1), DMAX),
        (DMAX, nmin),
        (nmin, DMAX),
        (sub(1), sub(1)),
        (4.0, 1.0),
        (9.0, 0.25),
        (2.0, 0.5),
        (sub(4), sub(1)),
        (1e300, 1e-100),
        (1.0, 1e-308),
        (DMAX, 1e-300),
    ] {
        vib(m, k, mu);
    }
    // √(k/μ) within about 2⁻¹⁰⁷ of a 53-bit midpoint m (the fast path's rounding test fails, the exact root decides):
    // k/μ = p/q·2^2j for m²'s last continued-fraction convergent p/q with p, q below 2⁵³ (found offline: cr_verify.py's
    // `_cf_best`, seed 30), relative distance from m² about 10⁻³²
    for (k, mu) in [
        (0x3E57_E304_6A1F_3250, 0x3DD3_D9DC_063D_8070),
        (0x2A6C_395D_F951_EF97, 0x2941_8CDE_D749_A2CF),
        (0x26DB_419D_7C34_BE51, 0x28C8_B543_1F13_0CE2),
        (0x5F51_FE81_6A7A_0BD4, 0x5D04_A920_083D_AC32),
        (0x3ECA_2DCF_E781_B324, 0x4124_4FF5_0297_1110),
        (0x4C78_EE5B_6FCA_4E99, 0x4B10_A790_EB0F_DE66),
        (0x2D72_43A9_D721_C14C, 0x2D7B_E4A6_EFCF_E06C),
        (0x3537_2495_0259_4498, 0x33A1_8FD5_99CF_A4E8),
        (0x34A9_33E8_5C5E_8E94, 0x349E_ECED_0426_B670),
        (0x5FD6_67A8_7077_B619, 0x5E3B_1275_06A6_F204),
        (0x469A_B94F_98E9_0FF4, 0x47AC_F96F_C17E_B0E0),
        (0x2673_D0B0_F3CA_A493, 0x282E_AE8B_6508_5268),
    ] {
        vib(m, f64::from_bits(k), f64::from_bits(mu));
    }
    for _ in 0..80 {
        // k/μ past DBL_MAX (Rust's ∞), below 2⁻¹⁰⁷⁴ (Rust's 0), subnormal operands, and the exact roots s²·μ/μ
        vib(m, rng.log_uniform(154.2, 308.0), rng.log_uniform(-308.0, -154.2));
        vib(m, rng.log_uniform(-308.0, -162.2), rng.log_uniform(162.2, 308.0));
        vib(m, rng.log_uniform(-324.0, -308.0), rng.log_uniform(-300.0, 308.0));
        vib(m, rng.log_uniform(-300.0, 308.0), rng.log_uniform(-324.0, -308.0));
        #[allow(clippy::cast_precision_loss)]
        let s = rng.int(1, 1 << 26) as f64;
        vib(m, s * s * pow2(i32::try_from(rng.int(-500, 500)).unwrap_or(0) * 2), 1.0);
    }
    for (j, b) in [(u32::MAX, 1e-300), (2, 1e-280), (1, 3e-275), (u32::MAX, DMAX), (0, sub(1)), (1, -1e-300), (7, sub(1))] {
        rotor(m, j, b);
    }
    for _ in 0..60 {
        let j = u32::try_from(rng.int(0, i64::from(u32::MAX))).unwrap_or(0);
        rotor(m, j, rng.signed_log_uniform(-324.0, -270.0));
        rotor(m, u32::try_from(rng.int(1, 100)).unwrap_or(1), rng.log_uniform(-324.0, -270.0));
    }
    for x in [sub(1), sub(2), sub(79), 1e-320, 2.8e-310, nmin, 1e-315] {
        rot_const(m, x);
    }
    for _ in 0..60 {
        rot_const(m, rng.log_uniform(-324.0, -309.5));
    }
    for (mu, r) in [
        (1e-300, 1e-5),
        (1e-310, 0.1),
        (sub(3), 0.5),
        (3.0, 67_108_865.0),
        (f64::INFINITY, 0.0),
        (-1e-300, 1e-6),
        (1e300, 1e-200),
    ] {
        moment(m, mu, r);
    }
    for _ in 0..60 {
        moment(m, rng.signed_log_uniform(-300.0, -100.0), rng.log_uniform(-110.0, -1.0));
        moment(m, rng.signed_log_uniform(-324.0, -300.0), rng.log_uniform(-10.0, 0.0));
        let r = odd_sq3(&mut rng);
        moment(m, 3.0, r * scale(&mut rng, -400, 400));
    }
}
