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

use crate::exact::{Ux, cr_monomial, cr_sqrt_monomial, cr_sqrt_ratio};
use crate::grids::{Rng, pow2};
use crate::m_reaction::res;
use crate::tsv::{Class, Module, Tok, f, i};
use kimiya::electrochemistry::{
    ELEMENTARY_CHARGE, FARADAY, faraday_charge_required, faraday_mass_deposited, faraday_mass_ratio, nernst_potential,
};
use kimiya::element::{AVOGADRO, GAS_CONSTANT};
use kimiya::error::KimiyaError;
use kimiya::gas::{
    PengRobinsonParams, gas_density, ideal_gas_pressure, ideal_gas_volume, peng_robinson_pressure, redlich_kwong_pressure,
    van_der_waals_pressure,
};
use kimiya::inorganic::born_lande_lattice_energy;
use kimiya::kinetics::{
    michaelis_menten, pre_equilibrium_rate, second_order_concentration, steady_state_intermediate,
    zero_order_half_life,
};
use kimiya::nuclear::{AMU_MEV, HYDROGEN_ATOM_MASS_U, NEUTRON_MASS_U, binding_energy_per_nucleon};
use kimiya::phase::{boiling_point_elevation, freezing_point_depression, osmotic_pressure};
use kimiya::potential::{COULOMB_CONSTANT, coulomb, coulomb_force, lennard_jones, morse_width_parameter};
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

// ── The overflow audit's B3 (ADR 0002 #31, #32) ─────────────────────────────────────────────────────────────────

/// How the port obtains a #20-style value whose final operation is a product (the note's `<how>`).
pub const UX_HOW_MUL: &str = "Rust's sequence evaluated with an unbounded exponent (the final product rounded once)";

/// A B3 port result: `Ok` with the port's value, or the name of its Err — Rust's validation's (`InvalidTemperature`,
/// `InvalidInput`), or [`CO_VOLUME`] for the three equations of state's V ≤ n·b or V ≤ b on the unbounded co-volume,
/// which can refuse where Rust's rounded one accepts or the reverse (an `InvalidInput` either way).
pub type Port = Result<f64, &'static str>;

/// The co-volume test's Err in a [`Port`] (Rust's `InvalidInput`).
pub const CO_VOLUME: &str = "co-volume";

fn ux(x: f64) -> Ux {
    Ux::from_f64(x)
}

fn err_name(e: &KimiyaError) -> &'static str {
    match e {
        KimiyaError::InvalidTemperature(_) => "InvalidTemperature",
        KimiyaError::InvalidInput(_) => "InvalidInput",
        _ => "other",
    }
}

/// A row of a `Result<f64>` function under B3: plain where the port's result is the oracle's (the same bits, or the same
/// Err), noted otherwise — `= <x> (f:<bits>)` for the port's value, `= Err(InvalidInput) (E:InvalidInput)` for its Err.
/// The port's validation is Rust's, so only the co-volume test may differ: asserted (the oracle's Err against the port's
/// value is then the co-volume test's, since the same checks passed in the port).
#[allow(clippy::too_many_arguments)] // noted_cls's
pub fn noted_r(
    m: &mut Module,
    name: &str,
    cls: Class,
    ins: &[Tok],
    oracle: Result<f64, KimiyaError>,
    port: Port,
    adr: u32,
    how: &str,
) {
    let pe = port.map_err(|n| if n == CO_VOLUME { "InvalidInput" } else { n });
    match (&oracle, pe) {
        (Ok(v), Ok(p)) if same(*v, p) => m.row(name, cls, ins, &[f(*v)]),
        (Err(e), Err(n)) => {
            assert_eq!(err_name(e), n, "{name}: the port's Err is not the oracle's: {ins:?}");
            m.row(name, cls, ins, &[res(oracle)]);
        }
        (Ok(_), Ok(p)) => m.row_noted(name, cls, ins, &[res(oracle)], &note(adr, how, p)),
        (Err(e), Ok(p)) => {
            assert_eq!(err_name(e), "InvalidInput", "{name}: the port is Ok where Rust refuses: {ins:?}");
            m.row_noted(name, cls, ins, &[res(oracle)], &note(adr, how, p));
        }
        (Ok(_), Err(_)) => {
            assert_eq!(port, Err(CO_VOLUME), "{name}: the port's validation refuses where the oracle's accepts: {ins:?}");
            let extra = format!("ADR 0002 #{adr}: the port returns {how} = Err(InvalidInput) (E:InvalidInput)");
            m.row_noted(name, cls, ins, &[res(oracle)], &extra);
        }
    }
}

/// ADR 0002 #31: `ideal_gas_pressure` and `ideal_gas_volume`, n·R·T/x rounded once (R as its double), after Rust's checks.
pub fn port_ideal(n: f64, t: f64, x: f64) -> Option<f64> {
    (n.is_finite() && pos(t) && pos(x)).then(|| cr_monomial(&[n, GAS_CONSTANT, t], &[x]))
}

/// ADR 0002 #31: `gas_density`, P·M/(R·T) rounded once, after Rust's checks.
pub fn port_density(p: f64, mm: f64, t: f64) -> Option<f64> {
    (pos(t) && nonneg(p) && pos(mm)).then(|| cr_monomial(&[p, mm], &[GAS_CONSTANT, t]))
}

/// ADR 0002 #31: `van_der_waals_pressure`, Rust's (((n·R)·T)/(V − n·b)) − ((a·(n/V))·(n/V)) with an unbounded exponent,
/// V ≤ n·b decided on the unbounded n·b, the final subtraction rounded once; Rust's checks in Rust's order.
pub fn port_vdw(n: f64, t: f64, v: f64, a: f64, b: f64) -> Port {
    if !pos(t) {
        return Err("InvalidTemperature");
    }
    if !n.is_finite() || !nonneg(a) || !nonneg(b) || !pos(v) {
        return Err("InvalidInput");
    }
    let nb = ux(n).mul(ux(b));
    if ux(v).le(nb) {
        return Err(CO_VOLUME);
    }
    let t1 = ux(n).mul(ux(GAS_CONSTANT)).mul(ux(t)).div(ux(v).sub(nb));
    let nv = ux(n).div(ux(v));
    Ok(t1.sub_f64(ux(a).mul(nv).mul(nv)))
}

/// Peng–Robinson's a·α and b with an unbounded exponent: Rust's `pr_coefficients`, operation for operation (α = u·u, as
/// Rust's `powi(2)`).
fn pr_coefficients_ux(tc: f64, pc: f64, omega: f64, t: f64) -> (Ux, Ux) {
    let (r, tcu, pcu) = (ux(GAS_CONSTANT), ux(tc), ux(pc));
    let a = ux(0.45724).mul(r).mul(r).mul(tcu).mul(tcu).div(pcu);
    let b = ux(0.07780).mul(r).mul(tcu).div(pcu);
    let om = ux(omega);
    let kappa = ux(0.37464).add(ux(1.54226).mul(om)).sub(ux(0.26992).mul(om).mul(om));
    let sq = ux(t).div(tcu).sqrt();
    let u = ux(1.0).add(kappa.mul(ux(1.0).sub(sq)));
    (a.mul(u.mul(u)), b)
}

/// ADR 0002 #31: `peng_robinson_pressure`, Rust's RT/(V − b) − a·α/(V² + 2bV − b²) with an unbounded exponent, V ≤ b on
/// the unbounded b, the final subtraction rounded once; Rust's checks in Rust's order.
pub fn port_pr(tc: f64, pc: f64, omega: f64, t: f64, v: f64) -> Port {
    if !pos(t) || !pos(tc) {
        return Err("InvalidTemperature");
    }
    if !pos(pc) || !omega.is_finite() {
        return Err("InvalidInput");
    }
    let (aa, b) = pr_coefficients_ux(tc, pc, omega, t);
    if !v.is_finite() {
        return Err("InvalidInput");
    }
    if ux(v).le(b) {
        return Err(CO_VOLUME);
    }
    let vu = ux(v);
    let t1 = ux(GAS_CONSTANT).mul(ux(t)).div(vu.sub(b));
    let den = vu.mul(vu).add(ux(2.0).mul(b).mul(vu)).sub(b.mul(b));
    Ok(t1.sub_f64(aa.div(den)))
}

/// ADR 0002 #31: `redlich_kwong_pressure`, Rust's RT/(V − b) − a/((√T·V)·(V + b)) with an unbounded exponent (Tc^2.5 its
/// 53-bit rounding at any exponent, [`Ux::pow25`]), V ≤ b on the unbounded b, the final subtraction rounded once; Rust's
/// checks in Rust's order.
pub fn port_rk(tc: f64, pc: f64, t: f64, v: f64) -> Port {
    if !pos(t) || !pos(tc) {
        return Err("InvalidTemperature");
    }
    if !pos(pc) {
        return Err("InvalidInput");
    }
    let (r, tcu, pcu) = (ux(GAS_CONSTANT), ux(tc), ux(pc));
    let a = ux(0.42748).mul(r).mul(r).mul(tcu.pow25()).div(pcu);
    let b = ux(0.08664).mul(r).mul(tcu).div(pcu);
    if !v.is_finite() {
        return Err("InvalidInput");
    }
    if ux(v).le(b) {
        return Err(CO_VOLUME);
    }
    let vu = ux(v);
    let t1 = r.mul(ux(t)).div(vu.sub(b));
    let den = ux(t).sqrt().mul(vu).mul(vu.add(b));
    Ok(t1.sub_f64(a.div(den)))
}

/// ADR 0002 #32: `lennard_jones`, Rust's 4·ε·((σ/r)¹² − (σ/r)⁶) with an unbounded exponent (4·ε exact), the final product
/// rounded once (#17's (ε·d)·4 where ε ≥ 2¹⁰²² and Rust's d is normal: the same value); Rust's checks.
pub fn port_lj(eps: f64, sigma: f64, r: f64) -> Option<f64> {
    (pos(eps) && pos(sigma) && pos(r)).then(|| {
        let sr = ux(sigma).div(ux(r));
        let sr6 = sr.mul(sr).mul(sr).mul(sr).mul(sr).mul(sr);
        let d = sr6.mul(sr6).sub(sr6);
        ux(4.0).mul(ux(eps)).mul_f64(d)
    })
}

/// ADR 0002 #32: `morse_width_parameter`, the correctly rounded ω·√(μ/(2·D_e)) = √(ω·ω·μ/(2·D_e)) (#18's edge
/// included), after Rust's checks.
pub fn port_width(w: f64, mu: f64, d: f64) -> Option<f64> {
    (pos(d) && pos(mu) && pos(w)).then(|| cr_sqrt_monomial(&[w, w, mu], &[2.0, d]))
}

/// ADR 0002 #32: `coulomb`, k_e·q₁·q₂/r rounded once (k_e as its double) for finite charges (Rust's ((k_e·q₁)·q₂)/r
/// otherwise), after Rust's check of r.
pub fn port_coulomb(q1: f64, q2: f64, r: f64) -> Option<f64> {
    pos(r).then(|| {
        if q1.is_finite() && q2.is_finite() { cr_monomial(&[COULOMB_CONSTANT, q1, q2], &[r]) } else { COULOMB_CONSTANT * q1 * q2 / r }
    })
}

/// ADR 0002 #32: `coulomb_force`, k_e·q₁·q₂/(r·r) rounded once for finite charges (Rust's ((k_e·q₁)·q₂)/(r·r) otherwise).
pub fn port_coulomb_force(q1: f64, q2: f64, r: f64) -> Option<f64> {
    pos(r).then(|| {
        if q1.is_finite() && q2.is_finite() {
            cr_monomial(&[COULOMB_CONSTANT, q1, q2], &[r, r])
        } else {
            COULOMB_CONSTANT * q1 * q2 / (r * r)
        }
    })
}

/// ADR 0002 #32: `born_lande_lattice_energy`, Rust's −k_e·N_A·M·z⁺·z⁻·e·e/r₀·(1 − 1/n) left to right with an unbounded
/// exponent, the final product rounded once; Rust's checks in Rust's order.
pub fn port_born_lande(madelung: f64, zc: u32, za: u32, r0: f64, n: f64) -> Option<f64> {
    (pos(r0) && n.is_finite() && n > 1.0 && pos(madelung) && zc != 0 && za != 0).then(|| {
        let e = ux(ELEMENTARY_CHARGE);
        let u = ux(-COULOMB_CONSTANT).mul(ux(AVOGADRO)).mul(ux(madelung)).mul(ux(f64::from(zc))).mul(ux(f64::from(za)));
        let u = u.mul(e).mul(e).div(ux(r0));
        u.mul_f64(ux(1.0).sub(ux(1.0).div(ux(n))))
    })
}

/// ADR 0002 #32: `binding_energy_per_nucleon`, Rust's ((Z·m_H + N·m_n − M)·931.494…)/A with an unbounded exponent, the
/// final division rounded once; Rust's checks in Rust's order.
pub fn port_bea(z: u32, n: u32, mass: f64) -> Option<f64> {
    let a = z.checked_add(n)?;
    (a != 0 && pos(mass)).then(|| {
        let md = ux(f64::from(z)).mul(ux(HYDROGEN_ATOM_MASS_U)).add(ux(f64::from(n)).mul(ux(NEUTRON_MASS_U))).sub(ux(mass));
        md.mul(ux(AMU_MEV)).div_f64(ux(f64::from(a)))
    })
}

pub fn igp(m: &mut Module, n: f64, t: f64, v: f64) {
    let out = ideal_gas_pressure(bb(n), bb(t), bb(v));
    let how = "the correctly rounded n·R·T/V";
    noted(m, "gas::ideal_gas_pressure", &[f(n), f(t), f(v)], out, port_ideal(n, t, v), 31, how);
}

pub fn igv(m: &mut Module, n: f64, t: f64, p: f64) {
    let out = ideal_gas_volume(bb(n), bb(t), bb(p));
    let how = "the correctly rounded n·R·T/P";
    noted(m, "gas::ideal_gas_volume", &[f(n), f(t), f(p)], out, port_ideal(n, t, p), 31, how);
}

pub fn dens(m: &mut Module, p: f64, mm: f64, t: f64) {
    let out = gas_density(bb(p), bb(mm), bb(t));
    let how = "the correctly rounded P·M/(R·T)";
    noted(m, "gas::gas_density", &[f(p), f(mm), f(t)], out, port_density(p, mm, t), 31, how);
}

pub fn vdw(m: &mut Module, n: f64, t: f64, v: f64, a: f64, b: f64) {
    let out = van_der_waals_pressure(bb(n), bb(t), bb(v), bb(a), bb(b));
    let ins = [f(n), f(t), f(v), f(a), f(b)];
    noted_r(m, "gas::van_der_waals_pressure", Class::B, &ins, out, port_vdw(n, t, v, a, b), 31, UX_HOW_SUB);
}

pub fn prp(m: &mut Module, p: &PengRobinsonParams, t: f64, v: f64) {
    let out = peng_robinson_pressure(bb(p), bb(t), bb(v));
    let ins = [f(p.tc), f(p.pc), f(p.omega), f(t), f(v)];
    let port = port_pr(p.tc, p.pc, p.omega, t, v);
    noted_r(m, "gas::peng_robinson_pressure", Class::B, &ins, out, port, 31, UX_HOW_SUB);
}

pub fn rk(m: &mut Module, tc: f64, pc: f64, t: f64, v: f64) {
    let out = redlich_kwong_pressure(bb(tc), bb(pc), bb(t), bb(v));
    let ins = [f(tc), f(pc), f(t), f(v)];
    noted_r(m, "gas::redlich_kwong_pressure", Class::C, &ins, out, port_rk(tc, pc, t, v), 31, UX_HOW_SUB);
}

pub fn lj(m: &mut Module, eps: f64, sigma: f64, r: f64) {
    let out = lennard_jones(bb(eps), bb(sigma), bb(r));
    noted(m, "potential::lennard_jones", &[f(eps), f(sigma), f(r)], out, port_lj(eps, sigma, r), 32, UX_HOW_MUL);
}

pub fn width(m: &mut Module, w: f64, mu: f64, d: f64) {
    let out = morse_width_parameter(bb(w), bb(mu), bb(d));
    let (name, how) = ("potential::morse_width_parameter", "the correctly rounded ω·√(μ/(2·D_e))");
    noted(m, name, &[f(w), f(mu), f(d)], out, port_width(w, mu, d), 32, how);
}

pub fn coul(m: &mut Module, q1: f64, q2: f64, r: f64) {
    let ins = [f(q1), f(q2), f(r)];
    let out = coulomb(bb(q1), bb(q2), bb(r));
    noted(m, "potential::coulomb", &ins, out, port_coulomb(q1, q2, r), 32, "the correctly rounded k_e·q₁·q₂/r");
    let out = coulomb_force(bb(q1), bb(q2), bb(r));
    let how = "the correctly rounded k_e·q₁·q₂/(r·r)";
    noted(m, "potential::coulomb_force", &ins, out, port_coulomb_force(q1, q2, r), 32, how);
}

pub fn lande(m: &mut Module, madelung: f64, zc: u32, za: u32, r0: f64, n: f64) {
    let out = born_lande_lattice_energy(bb(madelung), bb(zc), bb(za), bb(r0), bb(n));
    let ins = [f(madelung), i(zc), i(za), f(r0), f(n)];
    let port = port_born_lande(madelung, zc, za, r0, n);
    noted(m, "inorganic::born_lande_lattice_energy", &ins, out, port, 32, UX_HOW_MUL);
}

pub fn bea(m: &mut Module, z: u32, n: u32, mass: f64) {
    let out = binding_energy_per_nucleon(bb(z), bb(n), bb(mass));
    let name = "nuclear::binding_energy_per_nucleon";
    noted(m, name, &[i(z), i(n), f(mass)], out, port_bea(z, n, mass), 32, UX_HOW);
}

// ── B3's edge rows (appended to each module with their own seed, so the rows before keep theirs) ─────────────────

/// gas's B3 rows: the ideal gas law's n·R·T past DBL_MAX or below 2⁻¹⁰²² with V or P bringing it back, the true ends,
/// signed zeros, 53-bit midpoints (q·R for q = 9 … 15: R's odd part has 51 bits) and the subnormal grid's ties
/// (q·2⁻¹⁰⁷⁴·2^a·R/(R·2^(a + 1))); gas_density's R·T past DBL_MAX or subnormal, P·M past either end, its subnormal ties;
/// van der Waals' n·R·T, a·(n/V)² and n/V past either end (n/V = ∞ with a = 0: Rust's 0·∞ = NaN), V − n·b past DBL_MAX
/// (n·b = −∞), a subnormal n·b against a subnormal V (Rust's rounded n·b refuses where the exact one does not);
/// Peng–Robinson's subnormal Tc (T/Tc = ∞, Rust's NaN), Tc·Tc, ω², V·V and R·T past DBL_MAX, κ = −∞ at T = Tc (Rust's
/// −∞·0 = NaN), a subnormal b against V on its ulps; Redlich–Kwong's Tc^2.5 past DBL_MAX or below 2⁻¹⁰²² (exactly
/// rescaled), its exact midpoints there (Tc = t²·4^j with t⁵ of 54 bits), (√T·V)·(V + b) past DBL_MAX, a subnormal b
/// against V between Rust's b and the unbounded one.
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_gas(m: &mut Module) {
    let mut rng = Rng::new(0x6a31);
    let r = GAS_CONSTANT;
    for (n, t, v) in [
        (DMAX, 300.0, 1e10),
        (1e300, 1e10, 1e300),
        (-1e300, 1e10, 1e300),
        (DMAX, DMAX, 1.0),
        (1e-300, 1e-20, 1e-300),
        (sub(1), 1.0, 1.0),
        (sub(3), 0.25, 1.0),
        (-sub(5), 1.0, 2.0),
        (1e-310, 1.0, 1e-10),
        (sub(1), 1e-300, 1e10),
        (0.0, DMAX, sub(1)),
        (-0.0, 1e300, 1e-300),
        (9.0, 1.0, 1.0),
        (11.0, 1.0, 1.0),
        (-13.0, 1.0, 1.0),
        (15.0, 1.0, 1.0),
        (9.0, 3.0, 3.0),
        (sub(3), 1.0, r * 2.0),
    ] {
        igp(m, n, t, v);
        igv(m, n, t, v);
    }
    for _ in 0..40 {
        // n·R·T past DBL_MAX, V bringing it back; n·R·T below 2⁻¹⁰²², V bringing it up; subnormal results
        let (n, t) = (rng.signed_log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0));
        igp(m, n, t, rng.log_uniform(1.0, 308.0));
        igv(m, n, t, rng.log_uniform(1.0, 308.0));
        let (n, t) = (rng.signed_log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0));
        igp(m, n, t, rng.log_uniform(-320.0, -1.0));
        igv(m, n, t, rng.log_uniform(-320.0, -1.0));
        igp(m, rng.log_uniform(-320.0, -150.0), rng.log_uniform(-200.0, -100.0), rng.log_uniform(0.0, 30.0));
        // midpoints q·R·2^k (q = 9 … 15 odd) at scales across the range, and the subnormal grid's ties
        let q = [9.0, 11.0, 13.0, 15.0][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        igp(m, q * scale(&mut rng, -900, 900), scale(&mut rng, -60, 60), scale(&mut rng, -60, 60));
        igv(m, q * scale(&mut rng, -900, 900), 3.0 * scale(&mut rng, -60, 60), 3.0 * scale(&mut rng, -60, 60));
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        let b = i32::try_from(rng.int(-30, 30)).unwrap_or(0);
        let qs = sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a);
        igp(m, qs, pow2(b), r * pow2(a + b + 1));
    }
    for (p, mm, t) in [
        (101_325.0, 0.029, DMAX),
        (1e300, 1e10, DMAX),
        (DMAX, 10.0, 1e300),
        (1e-300, 1e-20, 1e-300),
        (1e-310, 1.0, sub(1)),
        (sub(1), 1e-10, 1e10),
        (0.0, DMAX, sub(1)),
        (-0.0, 1.0, 1.0),
        (sub(3), r, 2.0),
    ] {
        dens(m, p, mm, t);
    }
    for _ in 0..40 {
        dens(m, rng.log_uniform(154.0, 308.0), rng.log_uniform(154.0, 308.0), rng.log_uniform(300.0, 308.0));
        dens(m, rng.log_uniform(-1.0, 308.0), rng.log_uniform(-3.0, 3.0), rng.log_uniform(307.4, 308.25));
        dens(m, rng.log_uniform(-320.0, -154.0), rng.log_uniform(-300.0, -154.0), rng.log_uniform(-320.0, -10.0));
        dens(m, rng.log_uniform(-310.0, -300.0), rng.log_uniform(-3.0, 3.0), rng.log_uniform(-324.0, -308.0));
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        let b = i32::try_from(rng.int(-30, 30)).unwrap_or(0);
        let qs = sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a);
        dens(m, qs, r * pow2(b), pow2(a + b + 1));
    }
    for (n, t, v, a, b) in [
        (1e306, 300.0, 1e300, 0.0, 1e-10),
        (1e306, 300.0, 1e300, 1e-300, 1e-10),
        (1e300, 1e-20, 1e-10, 0.0, 0.0),
        (1e300, 1e-20, 1e-10, 1e-300, 0.0),
        (1e-300, 1e-30, 1e-300, 0.0, 0.0),
        (1e-300, 1e-30, 1e-300, 1e-300, 1e-20),
        (-1e300, 300.0, 1e300, 0.0, 1e10),
        (-1e300, 300.0, 1e300, 0.364, 1e10),
        (sub(3), 300.0, sub(1), 0.0, 0.2),
        (sub(3), 300.0, sub(1), 0.0, 0.4),
        (sub(3), 300.0, sub(2), 0.0, 0.4),
        (1e200, 300.0, 1e-100, 1e10, 0.0),
        (DMAX, DMAX, DMAX, 0.0, 0.0),
    ] {
        vdw(m, n, t, v, a, b);
    }
    for _ in 0..60 {
        // n·R·T past DBL_MAX with V and a small b keeping V > n·b; below 2⁻¹⁰²² with V tiny
        let n = rng.signed_log_uniform(300.0, 308.0);
        vdw(m, n, rng.log_uniform(1.0, 10.0), rng.log_uniform(300.0, 308.0), rng.unit(), rng.log_uniform(-12.0, -8.0));
        let n = rng.signed_log_uniform(-320.0, -290.0);
        vdw(m, n, rng.log_uniform(-30.0, 3.0), rng.log_uniform(-320.0, -290.0), rng.log_uniform(-320.0, 0.0), 0.0);
        // n/V past DBL_MAX (a·(n/V)² too) with a = 0 or tiny; V − n·b past DBL_MAX for a negative n
        let a = if rng.int(0, 1) == 0 { 0.0 } else { rng.log_uniform(-320.0, -300.0) };
        vdw(m, rng.log_uniform(250.0, 308.0), rng.log_uniform(-30.0, -10.0), rng.log_uniform(-100.0, -40.0), a, 0.0);
        vdw(m, -rng.log_uniform(290.0, 308.0), rng.log_uniform(1.0, 3.0), rng.log_uniform(280.0, 308.0), rng.unit(),
            rng.log_uniform(1.0, 30.0));
        // a subnormal n·b against a subnormal V on its grid
        let k = rng.int(1, 40);
        let n = sub(u64::try_from(rng.int(1, 1 << 20)).unwrap_or(1));
        let b = rng.uniform(0.0, 1.0);
        let nb = (n * b) / sub(1);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let j = (nb.round() as i64 + k - 20).max(1);
        vdw(m, n, 300.0, sub(u64::try_from(j).unwrap_or(1)), 0.0, b);
    }
    let prs = |tc: f64, pc: f64, om: f64| PengRobinsonParams { tc, pc, omega: om };
    for (tc, pc, om, t, v) in [
        (pow2(-1022), 3_394_000.0, 0.0377, 300.0, 0.025),
        (sub(1), 3_394_000.0, 0.0377, 300.0, 0.025),
        (sub(12_345), 3_394_000.0, 0.0377, 300.0, 1.0),
        (1e150, 1.0, 0.0, 300.0, 2e154),
        (1e200, 1e300, 0.0377, 300.0, 1e200),
        (126.2, 3_394_000.0, 1e200, 126.2, 0.025),
        (126.2, 3_394_000.0, -1e160, 126.2, 0.025),
        (126.2, 3_394_000.0, 1e155, 300.0, 0.025),
        (126.2, 3_394_000.0, 0.0377, DMAX, 1e300),
        (126.2, 3_394_000.0, 0.0377, 1e308, 1e306),
        (126.2, 3_394_000.0, 0.0377, 300.0, 1e160),
        (126.2, 3_394_000.0, 0.0377, 1e-300, 1e-300),
        (1e-300, 1e10, 0.0377, 300.0, 1e-300),
    ] {
        prp(m, &prs(tc, pc, om), t, v);
    }
    for _ in 0..60 {
        let tc = rng.log_uniform(-323.0, 308.0);
        let pc = rng.log_uniform(-320.0, 308.0);
        let om = if rng.int(0, 3) == 0 { rng.signed_log_uniform(100.0, 308.0) } else { rng.uniform(-0.5, 1.0) };
        let t = if rng.int(0, 3) == 0 { tc } else { rng.log_uniform(-320.0, 308.0) };
        let b = 0.07780 * GAS_CONSTANT * tc / pc;
        let v = if b.is_finite() && b > 0.0 { b * (1.0 + rng.log_uniform(-3.0, 5.0)) } else { rng.log_uniform(-300.0, 308.0) };
        prp(m, &prs(tc, pc, om), t, v);
        prp(m, &prs(126.2, 3_394_000.0, 0.0377), rng.log_uniform(-320.0, 308.0), rng.log_uniform(150.0, 308.0));
        // a subnormal b against V on its ulps
        let tc = rng.log_uniform(-320.0, -300.0);
        let pc = rng.log_uniform(1.0, 15.0);
        let b = 0.07780 * GAS_CONSTANT * tc / pc;
        let v = if b > 0.0 { f64::from_bits((b.to_bits() + 1).saturating_sub(u64::try_from(rng.int(0, 2)).unwrap_or(0))) } else { sub(1) };
        prp(m, &prs(tc, pc, 0.0377), 300.0, v);
    }
    // a subnormal Tc and a tiny Pc: Rust's (0.07780·R)·Tc rounds at the subnormal grid and /Pc scales its error back into
    // the normal range, so Rust's b and the unbounded one differ by up to a third; V between them, or on either, decides
    // V ≤ b differently (the port's Err where Rust accepts, and the reverse); Redlich–Kwong's 0.08664 likewise
    for _ in 0..30 {
        let tc = sub(u64::try_from(rng.int(1, 300)).unwrap_or(1));
        let pc = pow2(-i32::try_from(rng.int(30, 60)).unwrap_or(30));
        for (c, is_pr) in [(0.07780, true), (0.08664, false)] {
            let bf = c * GAS_CONSTANT * tc / pc;
            let bu = ux(c).mul(ux(GAS_CONSTANT)).mul(ux(tc)).div(ux(pc)).mul_f64(ux(1.0));
            for v in [bf, bu, 0.5 * (bf + bu), f64::from_bits(bf.to_bits() + 1), f64::from_bits(bu.to_bits() - 1)] {
                if is_pr {
                    prp(m, &prs(tc, pc, 0.0377), 300.0, v);
                } else {
                    rk(m, tc, pc, 300.0, v);
                }
            }
        }
    }
    for (tc, pc, t, v) in [
        (1e130, 1e300, 300.0, 1.0),
        (1e124, 1e10, 300.0, 1e120),
        (1e200, 1e300, 1e10, 1e100),
        (1e-125, 1e-300, 300.0, 1e175),
        (1e-200, 1e-300, 1e-10, 1e110),
        (126.2, 3_394_000.0, 1e200, 1e200),
        (126.2, 3_394_000.0, DMAX, 1e300),
        (sub(1), 1e-300, 300.0, 1e-50),
        (1e-300, 1e10, 300.0, 1e-300),
    ] {
        rk(m, tc, pc, t, v);
    }
    // t = 1553/1024 … 1781/1024 odd: t⁵ has 54 bits, so Tc = t²·4^j gives Tc^2.5 = t⁵·2^(5j) a 53-bit midpoint
    let mid = |rng: &mut Rng| {
        #[allow(clippy::cast_precision_loss)]
        let t = (rng.int(776, 890) * 2 + 1) as f64 / 1024.0;
        t * t
    };
    for _ in 0..60 {
        // Tc^2.5 past DBL_MAX with Pc keeping a = 0.42748·R²·Tc^2.5/Pc finite (log₁₀ Pc ≥ 2.5·log₁₀ Tc − 306), V above b
        // and large enough that a/(√T·V²) is finite (the decades bound both: log₁₀ a < 2.5·(k + 1) + 1.5 − kp)
        let k = rng.int(124, 240);
        #[allow(clippy::cast_precision_loss)]
        let (kf, kp) = (k as f64, rng.int((5 * (k + 1)) / 2 - 305, 307));
        #[allow(clippy::cast_precision_loss)]
        let kpf = kp as f64;
        let tc = rng.log_uniform(kf, kf + 1.0);
        let pc = rng.log_uniform(kpf, kpf + 1.0);
        let lv = (kf + 1.0 - kpf).max((2.5 * (kf + 1.0) + 1.5 - kpf - 300.0) / 2.0) + 0.5;
        rk(m, tc, pc, rng.log_uniform(-3.0, 4.0), rng.log_uniform(lv, lv + 20.0));
        // Tc^2.5 below 2⁻¹⁰²¹ with Pc keeping a normal (log₁₀ Pc ≤ 2.5·log₁₀ Tc + 308), V above b
        #[allow(clippy::cast_precision_loss)]
        let k = rng.int(-320, -124) as f64;
        let tc = rng.log_uniform(k, k + 1.0);
        let hi = (2.5 * k + 308.0).max(-319.0);
        let pc = rng.log_uniform((hi - 300.0).max(-320.0), hi);
        let b = 0.08664 * GAS_CONSTANT * tc / pc;
        let v = if b.is_finite() && b > 0.0 { b * rng.log_uniform(0.0005, 12.0) } else { rng.log_uniform(-300.0, 308.0) };
        rk(m, tc, pc, rng.log_uniform(-3.0, 4.0), v);
        // Tc^2.5's exact midpoints with Tc^2.5 past either end (the rescale's ties; in range glibc's pow, the oracle's,
        // misrounds such a tie, and those rows are ADR 0006's, not B3's)
        let j = i32::try_from(rng.int(-500, -205)).unwrap_or(0);
        let tc = mid(&mut rng) * pow2(2 * j);
        let pc = pow2(-1000);
        rk(m, tc, pc, 300.0, 0.08664 * GAS_CONSTANT * tc / pc * 4.0);
        let j = i32::try_from(rng.int(205, 400)).unwrap_or(0);
        let tc = mid(&mut rng) * pow2(2 * j);
        let pc = pow2(1000);
        rk(m, tc, pc, 300.0, pow2(((5 * j - 1000) / 2 + 20).max(2 * j - 996)));
        // (√T·V)·(V + b) past DBL_MAX
        rk(m, 126.2, 3_394_000.0, rng.log_uniform(150.0, 308.0), rng.log_uniform(150.0, 308.0));
    }
}

/// potential's B3 rows: lennard_jones' σ/r past 2^85.3 (Rust's +∞ where 4ε·d is finite) and past 2^170.7 (Rust's ∞ − ∞ =
/// NaN, the value finite for a tiny ε, else +∞), below 2^-179 with a huge ε (Rust's 0), a subnormal σ/r, ε ≥ 2¹⁰²²
/// (#17's) with σ/r out of range, the final product's subnormal ties; morse_width_parameter's μ/(2·D_e) past either end,
/// subnormal operands and results, exact roots, its 53-bit midpoints 3ω·2^k and the subnormal grid's ties, roots within
/// about 2⁻¹⁰⁷ of a midpoint (B2's twelve convergent pairs, as ω·√(k/(2·μ/2))); coulomb's and coulomb_force's k·q₁·q₂ and
/// r·r past either end, subnormal charges and results, their midpoints q·k_e (q = 5, 7: k_e's odd part has 52 bits),
/// the subnormal ties, non-finite charges (Rust's expression).
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_potential(m: &mut Module) {
    let mut rng = Rng::new(0x9032);
    for (e, s, r) in [
        (1e-300, 1.0, 1.723_292_942_818_807_6e-26),
        (sub(1), pow2(171), 1.0),
        (sub(1), pow2(200), 1.0),
        (1.0, pow2(171), 1.0),
        (1e300, 3e-54, 1.0),
        (1e300, 1e-60, 1.0),
        (pow2(1023), pow2(-350), 1.0),
        (3.0 * pow2(1017), pow2(-349), 1.0),
        (5.0 * pow2(1016), pow2(-349), 1.0),
        (1.0, sub(1), 1.0),
        (DMAX, sub(1), 1.0),
        (DMAX, 1e-60, 1.0),
        (DMAX, 1e30, 1.0),
        (pow2(1022), pow2(-200), 1.0),
        (1.0, 1.0, pow2(1000)),
        (1.0, pow2(-500), pow2(500)),
    ] {
        lj(m, e, s, r);
    }
    for _ in 0..60 {
        // (σ/r)¹² past DBL_MAX with ε tiny; σ/r past 2^170.7 (Rust's NaN); (σ/r)⁶ below 2⁻¹⁰²² with ε huge
        let s = rng.log_uniform(-200.0, 200.0);
        lj(m, rng.log_uniform(-323.0, -200.0), s, s / rng.log_uniform(25.7, 51.4));
        lj(m, rng.log_uniform(-323.0, 0.0), rng.log_uniform(51.5, 80.0), 1.0);
        lj(m, rng.log_uniform(280.0, 308.0), s, s / rng.log_uniform(-60.0, -53.9));
        lj(m, pow2(1022) * rng.uniform(1.0, 4.0), s, s / rng.log_uniform(-120.0, -50.0));
        lj(m, rng.log_uniform(-323.0, 308.0), rng.log_uniform(-323.0, -300.0), rng.log_uniform(-10.0, 10.0));
    }
    let y = |rng: &mut Rng| odd3(rng);
    for (w, mu, d) in [
        (563_400_000_000_000.0, DMAX, 7.399e-19),
        (563_400_000_000_000.0, 1.6267e-27, DMAX),
        (563_400_000_000_000.0, 1e300, 7.399e-19),
        (563_400_000_000_000.0, 1.6267e-27, 1e300),
        (1.0, sub(1), pow2(1023)),
        (1.0, sub(1), DMAX),
        (DMAX, DMAX, sub(1)),
        (sub(1), sub(1), DMAX),
        (sub(1), 1.0, 0.5),
        (1.0, 4.0, 0.5),
        (3.0, 18.0, 1.0),
        (3_002_399_751_580_331.0, 9.0, 0.5),
        (pow2(-975), 1.0, pow2(199)),
        (3.0 * pow2(-975), 1.0, pow2(199)),
        (5.0 * pow2(-976), 4.0, pow2(201)),
        (DMAX, 1.0, 0.5),
        (DMAX, 1.0000000000000002, 0.5),
    ] {
        width(m, w, mu, d);
    }
    // B2's √(k/μ) pairs within about 2⁻¹⁰⁷ of a 53-bit midpoint, as ω·√(k/(2·(μ/2))) with ω = 2^j (exact scalings)
    for (k, mu) in [
        (0x3E57_E304_6A1F_3250_u64, 0x3DD3_D9DC_063D_8070_u64),
        (0x2A6C_395D_F951_EF97, 0x2941_8CDE_D749_A2CF),
        (0x26DB_419D_7C34_BE51, 0x28C8_B543_1F13_0CE2),
        (0x5F51_FE81_6A7A_0BD4, 0x5D04_A920_083D_AC32),
        (0x3ECA_2DCF_E781_B324, 0x4124_4FF5_0297_1110),
        (0x4C78_EE5B_6FCA_4E99, 0x4B10_A790_EB0F_DE66),
    ] {
        let (kf, mf) = (f64::from_bits(k), f64::from_bits(mu));
        width(m, 1.0, kf, mf * 0.5);
        width(m, pow2(i32::try_from(rng.int(-200, 200)).unwrap_or(0)), kf, mf * 0.5);
    }
    for _ in 0..60 {
        // μ/(2·D_e) past DBL_MAX (Rust's ∞) and below 2⁻¹⁰⁷⁴ (Rust's 0); subnormal operands; subnormal results
        width(m, rng.log_uniform(-100.0, 100.0), rng.log_uniform(150.0, 308.0), rng.log_uniform(-320.0, -150.0));
        width(m, rng.log_uniform(-100.0, 300.0), rng.log_uniform(-320.0, -150.0), rng.log_uniform(150.0, 308.0));
        width(m, rng.log_uniform(-324.0, -308.0), rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0));
        width(m, rng.log_uniform(-100.0, 100.0), rng.log_uniform(-324.0, -308.0), rng.log_uniform(-300.0, 300.0));
        width(m, rng.log_uniform(-300.0, -200.0), rng.log_uniform(-300.0, -150.0), rng.log_uniform(150.0, 300.0));
        // the 53-bit midpoints 3ω·2^(k − m) and the subnormal grid's ties q·2⁻¹⁰⁷⁵, exact roots
        let (k, mm) = (i32::try_from(rng.int(-200, 200)).unwrap_or(0), i32::try_from(rng.int(-200, 200)).unwrap_or(0));
        width(m, y(&mut rng), 9.0 * pow2(2 * k), 0.5 * pow2(2 * mm));
        let s = i32::try_from(rng.int(60, 120)).unwrap_or(60);
        #[allow(clippy::cast_precision_loss)]
        let q = (rng.int(0, 1 << 40) * 2 + 1) as f64;
        width(m, q * pow2(-1075 + s), 1.0, 0.5 * pow2(2 * s));
        #[allow(clippy::cast_precision_loss)]
        let sq = rng.int(1, 1 << 26) as f64;
        width(m, rng.log_uniform(-100.0, 100.0), sq * sq * pow2(2 * k), 0.5 * pow2(2 * mm));
    }
    let e = ELEMENTARY_CHARGE;
    let kc = COULOMB_CONSTANT;
    for (q1, q2, r) in [
        (sub(1), e, 1e-10),
        (e, -e, 1e-160),
        (e, -e, 2.5e-168),
        (1e300, e, 1e-10),
        (-1e300, e, 1e-10),
        (1e300, 1e-300, 1.0),
        (1e-300, 1e-20, 1e-300),
        (1e100, 1e100, 1e155),
        (e, e, 1e155),
        (DMAX, DMAX, DMAX),
        (5.0, 1.0, 1.0),
        (-7.0, 1.0, 1.0),
        (5.0, 3.0, 3.0),
        (sub(3), 1.0, kc * 2.0),
        (f64::INFINITY, 0.0, 1.0),
        (f64::NAN, 1.0, 1.0),
        (f64::INFINITY, 1e-300, 1e-300),
        (0.0, -0.0, sub(1)),
        (-0.0, DMAX, DMAX),
    ] {
        coul(m, q1, q2, r);
    }
    for _ in 0..60 {
        // k·q₁·q₂ past DBL_MAX or below 2⁻¹⁰²² with r bringing it back; r·r past either end; subnormal charges
        coul(m, rng.signed_log_uniform(150.0, 308.0), rng.signed_log_uniform(0.0, 308.0), rng.log_uniform(100.0, 308.0));
        coul(m, rng.signed_log_uniform(-320.0, -150.0), rng.signed_log_uniform(-170.0, -20.0), rng.log_uniform(-320.0, -100.0));
        coul(m, rng.signed_log_uniform(100.0, 200.0), rng.signed_log_uniform(100.0, 200.0), rng.log_uniform(154.2, 250.0));
        coul(m, rng.signed_log_uniform(-200.0, -100.0), rng.signed_log_uniform(-100.0, -19.0), rng.log_uniform(-320.0, -154.2));
        coul(m, rng.signed_log_uniform(-324.0, -308.0), rng.signed_log_uniform(-30.0, 30.0), rng.log_uniform(-30.0, 30.0));
        // midpoints q·k_e·2^j, and the subnormal grid's ties q·2⁻¹⁰⁷⁴·2^a·k_e/(k_e·2^(a + 1))
        let q = [5.0, 7.0][usize::try_from(rng.int(0, 1)).unwrap_or(0)];
        coul(m, q * scale(&mut rng, -600, 600), scale(&mut rng, -60, 60), scale(&mut rng, -60, 60));
        let a = i32::try_from(rng.int(0, 50)).unwrap_or(0);
        let b = i32::try_from(rng.int(-30, 30)).unwrap_or(0);
        let qs = sub(u64::try_from(rng.int(1, 1 << 40)).unwrap_or(1) * 2 + 1) * pow2(a);
        coul(m, qs, pow2(b), kc * pow2(a + b + 1));
    }
}

/// inorganic's B3 rows: born_lande_lattice_energy's prefix −k_e·N_A·M·z⁺·z⁻ past DBL_MAX (M above about
/// 3·10²⁷⁴/(z⁺·z⁻)) with e² and r₀ bringing it back, a subnormal or tiny M (·e·e below 2⁻¹⁰²²) with a tiny r₀, the true
/// ends, n past 2¹⁰²¹ (1/n subnormal: 1 − 1/n = 1 in both) and n near 1.
pub fn emit_inorganic(m: &mut Module) {
    let mut rng = Rng::new(0x1c32);
    for (mm, zc, za, r0, n) in [
        (DMAX, 1, 1, DMAX, DMAX),
        (DMAX, 1, 1, 1e300, 8.0),
        (1e300, 1, 1, 1e300, 1e300),
        (1e280, 1000, 1000, 1e300, 8.0),
        (1e280, u32::MAX, u32::MAX, DMAX, 12.0),
        (sub(1), 1, 1, 2.82e-10, 8.0),
        (sub(1), u32::MAX, u32::MAX, 1e-300, 8.0),
        (1e-310, 1, 1, 1e-300, 8.0),
        (1e-320, 3, 2, 1e-305, 9.0),
        (1.748, 1, 1, 2.82e-10, DMAX),
        (1.748, 1, 1, 2.82e-10, 1e300),
        (DMAX, 1, 1, DMAX, 1.000_000_000_000_000_2),
        (DMAX, 2, 2, 1e-10, 8.0),
        (sub(1), 1, 1, 1.0, 8.0),
    ] {
        lande(m, mm, zc, za, r0, n);
    }
    for _ in 0..80 {
        let zc = u32::try_from(rng.int(1, 4)).unwrap_or(1);
        let za = if rng.int(0, 3) == 0 { u32::try_from(rng.int(1, i64::from(u32::MAX))).unwrap_or(1) } else { 1 };
        let n = rng.uniform(5.0, 12.0);
        lande(m, rng.log_uniform(274.0, 308.0), zc, za, rng.log_uniform(250.0, 308.0), n);
        lande(m, rng.log_uniform(-323.0, -300.0), zc, za, rng.log_uniform(-320.0, -290.0), n);
        lande(m, rng.log_uniform(-323.0, 308.0), zc, za, rng.log_uniform(-320.0, 308.0), rng.log_uniform(0.0, 308.0));
    }
}

/// nuclear's B3 rows: binding_energy_per_nucleon's Δm·931.494… past DBL_MAX (masses above about 1.93·10³⁰⁵ u) with
/// A ≥ 932 bringing it back (Rust's −∞), and A below that (the value itself past DBL_MAX: −∞ in both).
pub fn emit_nuclear(m: &mut Module) {
    let mut rng = Rng::new(0x2e32);
    for (z, n, mass) in [
        (1000, 1000, DMAX),
        (466, 466, DMAX),
        (465, 465, DMAX),
        (0, 932, DMAX),
        (u32::MAX - 1, 1, DMAX),
        (26, 30, DMAX),
        (1 << 31, (1 << 31) - 1, 1e306),
    ] {
        bea(m, z, n, mass);
    }
    for _ in 0..80 {
        let z = u32::try_from(rng.int(0, 1 << 31)).unwrap_or(0);
        let n = u32::try_from(rng.int(1, 1 << 31)).unwrap_or(1);
        bea(m, z, n, rng.log_uniform(305.3, 308.25));
        let z = u32::try_from(rng.int(0, 600)).unwrap_or(0);
        let n = u32::try_from(rng.int(1, 600)).unwrap_or(1);
        bea(m, z, n, rng.log_uniform(305.0, 308.25));
    }
}
