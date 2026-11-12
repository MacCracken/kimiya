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
//! **B2** (ADR 0002 #29, #30): solution's `dilution` and `solution_conductivity`, electrochemistry's
//! `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` (#29, rounded once) and `nernst_potential`
//! (#29, #20's semantics: Rust's E° − ((R·T)/(n·F))·ln Q with an unbounded exponent and the final subtraction rounded
//! once, ln Q the correctly rounded ln of the input, as the port's), and spectroscopy's `absorbance`,
//! `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia` (#30,
//! rounded once; π as its double each time Rust uses it, 8, 1000 and J + 1 exact) and `vibrational_frequency` (#30, the
//! correctly rounded √ of the exact k/μ, [`cr_sqrt_ratio`]).
//!
//! **B3** (ADR 0002 #31, #32): gas's ideal gas law, `gas_density` and the three equations of state (#31), and
//! potential's `lennard_jones`, `coulomb`, `coulomb_force`, `morse_width_parameter`, inorganic's
//! `born_lande_lattice_energy` and nuclear's `binding_energy_per_nucleon` (#32), with their own section below.
//!
//! **B4** (ADR 0002 #33, #34, and a correction of #29): #20's semantics, Rust's sequence (its order of operations and
//! summation) with an unbounded exponent and the final operation rounded once, for thermo's `hess_law` and
//! `enthalpy_from_formation`, thermochem's `reaction_enthalpy`, `reaction_gibbs_energy` (ADR 0002 #1's (Σ kJ)·1000) and
//! `reaction_entropy`, nuclear's `q_value`, solution's `ionic_strength`, `kohlrausch_molar_conductivity` and
//! `limiting_molar_conductivity`, inorganic's `born_haber_lattice_energy`, reaction's `gibbs_free_energy` and kinetics'
//! `zero_order_concentration` (#33), and `Shomate::cp` and `Shomate::enthalpy_antiderivative` on Rust's exact
//! polynomial order (#34); and #29's `nernst_potential` rows whose product (R·T)/(n·F)·ln Q IEEE rounds up to 2⁻¹⁰²²
//! from below ([`NERNST_2M1022`]).
//!
//! **B5** (ADR 0002 #35): the exp / exp10 *arguments* of reaction's `equilibrium_constant` (−ΔG/(R·T)) and solution's
//! `debye_huckel_limiting`, `debye_huckel_extended` and `davies_activity` (log₁₀ γ): Rust's sequence for the argument
//! with an unbounded exponent, its final quotient or product rounded once into binary64, then the correctly rounded
//! exp / 10^x of that argument (glibc's here, as B2's ln: `validate.py` checks every noted value against MPFR's exp /
//! exp10 of the 53-bit unbounded argument), and `equilibrium_constant`'s ComputationError judged on that K.
//!
//! Every row of these functions goes through [`noted`]: where the port's value differs from the oracle's (by its bits;
//! two NaNs are the same) the row is noted `[ADR 0002 #n: the port returns <how> = <x> (f:<bits>)]`. Rust's validation is
//! replicated here and must agree with the oracle's (asserted on every row). `validate.py` proves every note again,
//! with `Fraction` and MPFR.

use std::hint::black_box as bb;

use crate::exact::{Ux, cr_monomial, cr_sqrt_monomial, cr_sqrt_ratio};
use crate::grids::{Rng, pow2};
use crate::m_reaction::res;
use crate::tsv::{Class, Module, Tok, f, i, s};
use kimiya::electrochemistry::{
    ELEMENTARY_CHARGE, FARADAY, faraday_charge_required, faraday_mass_deposited, faraday_mass_ratio, nernst_potential,
};
use crate::t_thermochem::shomate_fields;
use kimiya::inorganic::born_haber_lattice_energy;
use kimiya::kinetics::zero_order_concentration;
use kimiya::nuclear::q_value;
use kimiya::reaction::{equilibrium_constant, gibbs_free_energy};
use kimiya::solution::{
    davies_activity, debye_huckel_extended, debye_huckel_limiting, ionic_strength, kohlrausch_molar_conductivity,
    limiting_molar_conductivity,
};
use kimiya::thermo::{enthalpy_from_formation, hess_law};
use kimiya::thermochem::{
    Shomate, THERMOCHEM_DATA, lookup_thermochem, reaction_enthalpy, reaction_entropy, reaction_gibbs_energy,
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
/// underflowing, the true ends, and the 53-bit midpoints 3F·y/F and y·3·F/F (3F is exact: F's odd part has 51 bits);
/// then B4's correction of #29 ([`NERNST_2M1022`]: a product c·l that IEEE rounds up to 2⁻¹⁰²² from below).
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
    // B4's correction of #29 (appended after B2's rows, so they keep theirs): nernst_potential where the product
    // c·l = ((R·T)/(n·F))·ln Q lies in [2⁻¹⁰²² − 2⁻¹⁰⁷⁵, 2⁻¹⁰²² − 2⁻¹⁰⁷⁶), which IEEE rounds up to 2⁻¹⁰²² and 53-bit
    // rounding to 2⁻¹⁰²² − 2⁻¹⁰⁷⁵. B2's fast path took 2⁻¹⁰²² itself as normal and kept Rust's value; with E° an odd
    // multiple of 2⁻¹⁰⁷⁴ (±1, 3) #20's value differs (a subnormal tie, to even), with an even one (0, 2) it does not.
    for &(e0, n, t, q) in &NERNST_2M1022 {
        nernst(m, f64::from_bits(e0), n, f64::from_bits(t), f64::from_bits(q));
    }
}

/// ADR 0002 #29's B4 correction: (E°, n, T, Q) bits with c·l in [2⁻¹⁰²² − 2⁻¹⁰⁷⁵, 2⁻¹⁰²² − 2⁻¹⁰⁷⁶) (found by a seeded
/// search with MPFR's log, which glibc's equals on each; the first three are the regression the main session found).
const NERNST_2M1022: [(u64, u8, u64, u64); 23] = [
    (0x0000_0000_0000_0001, 1, 0x012c_5299_aebb_56a1, 0x3ff0_d212_8b2f_3300),
    (0x0000_0000_0000_0001, 1, 0x015d_a605_365d_ec16, 0x3ff0_1889_2f90_2b00),
    (0x0000_0000_0000_0001, 1, 0x013b_5b45_c314_13cc, 0x3ff0_6b6f_0367_5a00),
    (0x8000_0000_0000_0001, 7, 0x012d_4ef9_5f21_0a80, 0x3fe6_d090_03db_9e8a),
    (0x8000_0000_0000_0001, 1, 0x0136_0e4f_6050_3912, 0x3fee_fd1b_c883_8237),
    (0x8000_0000_0000_0001, 255, 0x0184_59b4_eaac_df85, 0x3fe8_3fd3_5aa3_b4ff),
    (0x8000_0000_0000_0001, 255, 0x0191_b54b_d246_5e98, 0x3feb_4932_0184_af3e),
    (0x0000_0000_0000_0001, 2, 0x011f_cd71_eb55_3fde, 0x3ff3_1ed6_b393_7b45),
    (0x0000_0000_0000_0001, 3, 0x015a_f517_7370_7b06, 0x3ff0_5183_cbac_8628),
    (0x8000_0000_0000_0001, 7, 0x012b_48b8_639e_6baa, 0x3fe6_3fc4_c95c_fefc),
    (0x0000_0000_0000_0001, 1, 0x0111_b58b_7e7f_1aa1, 0x3ff2_c69c_5171_ec09),
    (0x8000_0000_0000_0001, 1, 0x0122_3118_d417_935f, 0x3fed_9a4e_14f8_84e7),
    (0x0000_0000_0000_0001, 3, 0x012c_c67d_7a4b_bdf0, 0x3ff2_8bdf_0a3b_5122),
    (0x0000_0000_0000_0001, 255, 0x0179_1b07_966e_9d3b, 0x3ff9_1571_2636_d1db),
    (0x0000_0000_0000_0003, 1, 0x012c_5299_aebb_56a1, 0x3ff0_d212_8b2f_3300),
    (0x0000_0000_0000_0002, 1, 0x012c_5299_aebb_56a1, 0x3ff0_d212_8b2f_3300),
    (0x0000_0000_0000_0000, 1, 0x012c_5299_aebb_56a1, 0x3ff0_d212_8b2f_3300),
    (0x0000_0000_0000_0003, 1, 0x015d_a605_365d_ec16, 0x3ff0_1889_2f90_2b00),
    (0x0000_0000_0000_0002, 1, 0x015d_a605_365d_ec16, 0x3ff0_1889_2f90_2b00),
    (0x0000_0000_0000_0000, 1, 0x015d_a605_365d_ec16, 0x3ff0_1889_2f90_2b00),
    (0x0000_0000_0000_0003, 1, 0x013b_5b45_c314_13cc, 0x3ff0_6b6f_0367_5a00),
    (0x0000_0000_0000_0002, 1, 0x013b_5b45_c314_13cc, 0x3ff0_6b6f_0367_5a00),
    (0x0000_0000_0000_0000, 1, 0x013b_5b45_c314_13cc, 0x3ff0_6b6f_0367_5a00),
];

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
        KimiyaError::ComputationError(_) => "ComputationError",
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

// ── The overflow audit's B4 (ADR 0002 #33, #34) ─────────────────────────────────────────────────────────────────

/// How the port obtains a #20-style value whose final operation is an addition (the note's `<how>`).
pub const UX_HOW_ADD: &str = "Rust's sequence evaluated with an unbounded exponent (the final addition rounded once)";

const ADR1: &str = "J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002 #1)";

fn finite_all(xs: &[f64]) -> bool {
    xs.iter().all(|x| x.is_finite())
}

/// Rust's `iter().sum()` with an unbounded exponent: from −0.0 in order, every addition rounded to 53 bits.
fn sum_ux(xs: &[f64]) -> Ux {
    xs.iter().fold(ux(-0.0), |s, &x| s.add(ux(x)))
}

/// ADR 0002 #33: `hess_law`, Rust's Σ from −0.0 with an unbounded exponent, the last addition rounded once (−0.0 for no
/// steps); Rust's sum for a non-finite step.
pub fn port_hess(xs: &[f64]) -> f64 {
    if !finite_all(xs) || xs.is_empty() {
        return hess_law(xs);
    }
    let (head, last) = xs.split_at(xs.len() - 1);
    sum_ux(head).add_f64(ux(last[0]))
}

/// ADR 0002 #33: `enthalpy_from_formation`, Rust's Σp − Σr with an unbounded exponent, the difference rounded once;
/// Rust's for a non-finite value.
pub fn port_formation(p: &[f64], r: &[f64]) -> f64 {
    if !finite_all(p) || !finite_all(r) {
        return enthalpy_from_formation(p, r);
    }
    sum_ux(p).sub_f64(sum_ux(r))
}

/// ADR 0002 #33: `q_value`, Rust's (Σr − Σp)·931.494… with an unbounded exponent, the product rounded once; Rust's for
/// a non-finite mass.
pub fn port_q(r: &[f64], p: &[f64]) -> f64 {
    if !finite_all(r) || !finite_all(p) {
        return q_value(r, p);
    }
    sum_ux(r).sub(sum_ux(p)).mul_f64(ux(AMU_MEV))
}

/// ADR 0002 #33: `ionic_strength`, Rust's ½·Σ (c·z)·z from −0.0 with an unbounded exponent, the product by ½ rounded
/// once; Rust's for a non-finite concentration.
pub fn port_ionic(xs: &[(f64, i32)]) -> f64 {
    if !xs.iter().all(|(c, _)| c.is_finite()) {
        return ionic_strength(xs);
    }
    let s = xs.iter().fold(ux(-0.0), |s, &(c, z)| {
        let zf = ux(f64::from(z));
        s.add(ux(c).mul(zf).mul(zf))
    });
    ux(0.5).mul_f64(s)
}

/// ADR 0002 #33: `born_haber_lattice_energy`, Rust's ΔH_f − ΔH_sub − IE − ΔH_diss/2 − EA left to right with an unbounded
/// exponent, the last difference rounded once; Rust's for a non-finite input.
pub fn port_haber(x: [f64; 5]) -> f64 {
    if !finite_all(&x) {
        return born_haber_lattice_energy(x[0], x[1], x[2], x[3], x[4]);
    }
    ux(x[0]).sub(ux(x[1])).sub(ux(x[2])).sub(ux(x[3]).div(ux(2.0))).sub_f64(ux(x[4]))
}

/// ADR 0002 #33: `gibbs_free_energy`, Rust's ΔH − T·ΔS with an unbounded exponent, the difference rounded once, after
/// Rust's checks.
pub fn port_gibbs(h: f64, t: f64, s: f64) -> Option<f64> {
    (pos(t) && h.is_finite() && s.is_finite()).then(|| ux(h).sub_f64(ux(t).mul(ux(s))))
}

/// ADR 0002 #33: `kohlrausch_molar_conductivity`, Rust's Λ° − K·√c with an unbounded exponent, the difference rounded
/// once, after Rust's checks.
pub fn port_kohlrausch(l0: f64, k: f64, c: f64) -> Option<f64> {
    (nonneg(c) && l0.is_finite() && k.is_finite()).then(|| ux(l0).sub_f64(ux(k).mul(ux(c).sqrt())))
}

/// ADR 0002 #33: `limiting_molar_conductivity`, Rust's ν₊·λ₊ + ν₋·λ₋ with an unbounded exponent, the sum rounded once;
/// Rust's for a non-finite input.
pub fn port_limiting(cat: f64, nc: f64, an: f64, na: f64) -> f64 {
    if !finite_all(&[cat, nc, an, na]) {
        return limiting_molar_conductivity(cat, nc, an, na);
    }
    ux(nc).mul(ux(cat)).add_f64(ux(na).mul(ux(an)))
}

/// ADR 0002 #33: `zero_order_concentration`, Rust's [A]₀ − k·t with an unbounded exponent, 0 where that difference is
/// negative (the comparison on the unbounded value), otherwise the difference rounded once; Rust's for a non-finite
/// input.
pub fn port_zero_conc(c0: f64, k: f64, t: f64) -> f64 {
    if !finite_all(&[c0, k, t]) {
        return zero_order_concentration(c0, k, t);
    }
    let p = ux(k).mul(ux(t));
    if ux(c0).sub(p).cmp(ux(0.0)) == std::cmp::Ordering::Less {
        return 0.0;
    }
    ux(c0).sub_f64(p)
}

/// The thermochem column a reaction sum reads: 0 ΔH_f°, 1 ΔG_f°, 2 S°.
fn tc_prop(formula: &str, col: u8) -> Option<f64> {
    lookup_thermochem(formula).map(|d| match col {
        0 => d.delta_hf_kj,
        1 => d.delta_gf_kj,
        _ => d.s_standard_j,
    })
}

/// Rust's `sum_property` with an unbounded exponent: from +0.0, `total += coeff * prop`, every operation rounded to 53
/// bits; `None` for an unknown formula (Rust's Err).
fn sum_property_ux(species: &[(&str, f64)], col: u8) -> Option<Ux> {
    let mut total = ux(0.0);
    for &(fo, c) in species {
        total = total.add(ux(c).mul(ux(tc_prop(fo, col)?)));
    }
    Some(total)
}

/// ADR 0002 #33: thermochem's reaction sums, the port's J/mol (`col` 0, 1: ADR 0002 #1's Σp − Σr × 1000) or J/(mol·K)
/// (`col` 2: Σp − Σr) with an unbounded exponent, the final operation rounded once; Rust's value (× 1000) for a
/// non-finite coefficient; `None` where Rust errs (an unknown formula, the products first).
pub fn port_reaction(p: &[(&str, f64)], r: &[(&str, f64)], col: u8) -> Option<f64> {
    if p.iter().chain(r).any(|(fo, _)| tc_prop(fo, col).is_none()) {
        return None;
    }
    if !p.iter().chain(r).all(|(_, c)| c.is_finite()) {
        let kj = match col {
            0 => reaction_enthalpy(p, r),
            1 => reaction_gibbs_energy(p, r),
            _ => return reaction_entropy(p, r).ok(),
        };
        return kj.ok().map(|v| v * 1000.0);
    }
    let sp = sum_property_ux(p, col)?;
    let sr = sum_property_ux(r, col)?;
    Some(if col == 2 { sp.sub_f64(sr) } else { sp.sub(sr).mul_f64(ux(1000.0)) })
}

/// ADR 0002 #34: `Shomate::cp`, Rust's (((A + B·t) + (C·t)·t) + ((D·t)·t)·t) + E/(t·t), t = T/1000, with an unbounded
/// exponent, the last addition rounded once; Rust's for a non-finite coefficient or T, and for T = ±0 (every product an
/// exact zero and E/0 IEEE's special in both evaluations).
pub fn port_cp(x: &Shomate, temp: f64) -> f64 {
    if !finite_all(&[x.a, x.b, x.c, x.d, x.e, temp]) || temp == 0.0 {
        return x.cp(temp);
    }
    let t = ux(temp).div(ux(1000.0));
    let s = ux(x.a).add(ux(x.b).mul(t)).add(ux(x.c).mul(t).mul(t)).add(ux(x.d).mul(t).mul(t).mul(t));
    s.add_f64(ux(x.e).div(t.mul(t)))
}

/// ADR 0002 #34: `Shomate::enthalpy_antiderivative`, Rust's 1000·(t·(A + t·(B/2 + t·(C/3 + t·(D/4)))) − E/t) with an
/// unbounded exponent, the product by 1000 rounded once; Rust's as [`port_cp`].
pub fn port_h(x: &Shomate, temp: f64) -> f64 {
    if !finite_all(&[x.a, x.b, x.c, x.d, x.e, temp]) || temp == 0.0 {
        return x.enthalpy_antiderivative(temp);
    }
    let t = ux(temp).div(ux(1000.0));
    let inner3 = ux(x.c).div(ux(3.0)).add(t.mul(ux(x.d).div(ux(4.0))));
    let inner2 = ux(x.b).div(ux(2.0)).add(t.mul(inner3));
    let inner = t.mul(ux(x.a).add(t.mul(inner2)));
    ux(1000.0).mul_f64(inner.sub(ux(x.e).div(t)))
}

pub fn hess(m: &mut Module, ins: &[Tok], xs: &[f64]) {
    let out = hess_law(bb(xs));
    noted_f(m, "thermo::hess_law", ins, out, port_hess(xs), 33, UX_HOW_ADD);
}

pub fn formation(m: &mut Module, ins: &[Tok], p: &[f64], r: &[f64]) {
    let out = enthalpy_from_formation(bb(p), bb(r));
    noted_f(m, "thermo::enthalpy_from_formation", ins, out, port_formation(p, r), 33, UX_HOW_SUB);
}

pub fn qv(m: &mut Module, ins: &[Tok], r: &[f64], p: &[f64]) {
    let out = q_value(bb(r), bb(p));
    noted_f(m, "nuclear::q_value", ins, out, port_q(r, p), 33, UX_HOW_MUL);
}

pub fn ionic(m: &mut Module, ins: &[Tok], xs: &[(f64, i32)]) {
    let out = ionic_strength(bb(xs));
    noted_f(m, "solution::ionic_strength", ins, out, port_ionic(xs), 33, UX_HOW_MUL);
}

pub fn haber(m: &mut Module, x: [f64; 5]) {
    let out = born_haber_lattice_energy(bb(x[0]), bb(x[1]), bb(x[2]), bb(x[3]), bb(x[4]));
    let ins: Vec<Tok> = x.iter().map(|&v| f(v)).collect();
    noted_f(m, "inorganic::born_haber_lattice_energy", &ins, out, port_haber(x), 33, UX_HOW_SUB);
}

pub fn gibbs(m: &mut Module, h: f64, t: f64, s: f64) {
    let out = gibbs_free_energy(bb(h), bb(t), bb(s));
    noted(m, "reaction::gibbs_free_energy", &[f(h), f(t), f(s)], out, port_gibbs(h, t, s), 33, UX_HOW_SUB);
}

pub fn kohl(m: &mut Module, l0: f64, k: f64, c: f64) {
    let out = kohlrausch_molar_conductivity(bb(l0), bb(k), bb(c));
    let (name, ins) = ("solution::kohlrausch_molar_conductivity", [f(l0), f(k), f(c)]);
    noted(m, name, &ins, out, port_kohlrausch(l0, k, c), 33, UX_HOW_SUB);
}

pub fn limiting(m: &mut Module, cat: f64, nc: f64, an: f64, na: f64) {
    let out = limiting_molar_conductivity(bb(cat), bb(nc), bb(an), bb(na));
    let (name, ins) = ("solution::limiting_molar_conductivity", [f(cat), f(nc), f(an), f(na)]);
    noted_f(m, name, &ins, out, port_limiting(cat, nc, an, na), 33, UX_HOW_ADD);
}

pub fn zero_conc(m: &mut Module, c0: f64, k: f64, t: f64) {
    let out = zero_order_concentration(bb(c0), bb(k), bb(t));
    let (name, ins) = ("kinetics::zero_order_concentration", [f(c0), f(k), f(t)]);
    noted_f(m, name, &ins, out, port_zero_conc(c0, k, t), 33, UX_HOW_SUB);
}

/// thermochem's three reaction-sum rows of one reaction (`ins` its tokens): reaction_enthalpy and reaction_gibbs_energy
/// with ADR 0002 #1's note (the oracle's kJ/mol × 1000.0), and #33's after it where the port's value differs; and
/// reaction_entropy.
pub fn sums(m: &mut Module, ins: &[Tok], p: &[(&str, f64)], r: &[(&str, f64)]) {
    for (col, name) in [(0_u8, "thermochem::reaction_enthalpy"), (1, "thermochem::reaction_gibbs_energy")] {
        let oracle = if col == 0 { reaction_enthalpy(bb(p), bb(r)) } else { reaction_gibbs_energy(bb(p), bb(r)) };
        let oracle = oracle.map(|kj| kj * 1000.0);
        let port = port_reaction(p, r, col);
        assert_eq!(oracle.is_ok(), port.is_some(), "{name}: the port's validation differs from the oracle's: {ins:?}");
        match (&oracle, port) {
            (Ok(v), Some(pv)) if !same(*v, pv) => {
                let extra = format!("{ADR1}] [{}", note(33, UX_HOW_MUL, pv));
                m.row_noted(name, Class::B, ins, &[f(*v)], &extra);
            }
            _ => m.row_noted(name, Class::B, ins, &[res(oracle)], ADR1),
        }
    }
    let out = reaction_entropy(bb(p), bb(r));
    noted(m, "thermochem::reaction_entropy", ins, out, port_reaction(p, r, 2), 33, UX_HOW_SUB);
}

/// thermochem's Shomate rows at T (`ins` the seven fields and T).
pub fn shomate(m: &mut Module, ins: &[Tok], x: &Shomate, t: f64) {
    let cp = bb(x).cp(bb(t));
    noted_f(m, "Shomate::cp", ins, cp, port_cp(x, t), 34, UX_HOW_ADD);
    let h = bb(x).enthalpy_antiderivative(bb(t));
    noted_f(m, "Shomate::enthalpy_antiderivative", ins, h, port_h(x, t), 34, UX_HOW_MUL);
}

// ── B4's edge rows (appended to each module with their own seed, so the rows before keep theirs) ─────────────────

fn seq_toks(xs: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(xs.len())];
    v.extend(xs.iter().map(|&x| f(x)));
    v
}

/// A signed value near DBL_MAX: ±DBL_MAX·u, u in [lo, hi].
fn big(rng: &mut Rng, lo: f64, hi: f64) -> f64 {
    let s = if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
    s * DMAX * rng.uniform(lo, hi)
}

/// thermo's B4 rows: hess_law's partial sums past ±DBL_MAX that later steps bring back (Rust's ±∞, or ∞ − ∞ = NaN), the
/// rounding past DBL_MAX itself (DBL_MAX + 2⁹⁷⁰), sums that are themselves out of range (±∞ in both), subnormal sums
/// (exact: Rust's bits); enthalpy_from_formation's two sums past DBL_MAX that cancel, one past it against a finite one,
/// a difference past DBL_MAX.
pub fn emit_thermo_b4(m: &mut Module) {
    let mut rng = Rng::new(0x7433);
    let h970 = pow2(970);
    for xs in [
        vec![DMAX, DMAX, -DMAX],
        vec![DMAX, DMAX, -DMAX, -DMAX],
        vec![-DMAX, -DMAX, DMAX],
        vec![DMAX, 1e308, -1e308],
        vec![1e308, 1e308, -1e308, 1.0],
        vec![DMAX, DMAX],
        vec![DMAX, DMAX, -DMAX, -DMAX, sub(1)],
        vec![DMAX, h970],
        vec![DMAX, h970, -DMAX],
        vec![DMAX, h970, -h970],
        vec![-DMAX, -h970, DMAX, 1.0],
        vec![sub(1), sub(1), -sub(3)],
        vec![sub(5), -sub(5)],
        vec![DMAX, -DMAX],
        vec![0.9 * DMAX, 0.9 * DMAX, -0.9 * DMAX, -0.9 * DMAX, 0.5],
    ] {
        hess(m, &seq_toks(&xs), &xs);
    }
    for _ in 0..80 {
        let n = usize::try_from(rng.int(2, 8)).unwrap_or(2);
        let xs: Vec<f64> = (0..n).map(|_| big(&mut rng, 0.3, 1.0)).collect();
        hess(m, &seq_toks(&xs), &xs);
    }
    for (p, r) in [
        (vec![DMAX, DMAX], vec![DMAX, DMAX]),
        (vec![DMAX, DMAX], vec![DMAX]),
        (vec![DMAX, DMAX, -DMAX], vec![]),
        (vec![], vec![DMAX, DMAX, -DMAX]),
        (vec![1e308, 1e308], vec![1e308, 0.999e308]),
        (vec![DMAX], vec![-DMAX]),
        (vec![-DMAX, -DMAX], vec![-DMAX, -DMAX, 1.0]),
        (vec![DMAX, h970], vec![DMAX]),
        (vec![sub(3)], vec![sub(1), sub(1)]),
    ] {
        let mut ins = seq_toks(&p);
        ins.extend(seq_toks(&r));
        formation(m, &ins, &p, &r);
    }
    for _ in 0..80 {
        let np = usize::try_from(rng.int(0, 5)).unwrap_or(0);
        let nr = usize::try_from(rng.int(0, 5)).unwrap_or(0);
        let p: Vec<f64> = (0..np).map(|_| big(&mut rng, 0.3, 1.0)).collect();
        let r: Vec<f64> = (0..nr).map(|_| big(&mut rng, 0.3, 1.0)).collect();
        let mut ins = seq_toks(&p);
        ins.extend(seq_toks(&r));
        formation(m, &ins, &p, &r);
    }
}

/// reaction's B4 rows: gibbs_free_energy's T·ΔS past DBL_MAX with ΔH bringing it back (Rust's ±∞), the true ends, T·ΔS
/// subnormal against a tiny ΔH (Rust's rounding at the subnormal grid twice), and T·ΔS below 2⁻¹⁰⁷⁵ (Rust's 0).
pub fn emit_reaction(m: &mut Module) {
    let mut rng = Rng::new(0x6b33);
    for (h, t, s) in [
        (DMAX, 2.0, 0.75 * DMAX),
        (-DMAX, 2.0, -0.75 * DMAX),
        (DMAX, DMAX, 1.0),
        (0.0, 1e-300, 1e-20),
        (sub(3), 1.5, sub(1)),
        (sub(1), 0.75, sub(1)),
        (-sub(3), 1.5, -sub(1)),
        (0.0, 1e-200, 1e-200),
        (-0.0, 1e-200, 1e-200),
        (1.0, 1e-200, -1e-200),
        (1e300, 1e300, 1e10),
        (DMAX, 1.5, DMAX),
    ] {
        gibbs(m, h, t, s);
    }
    for _ in 0..80 {
        // T·ΔS in (DBL_MAX, 2·DBL_MAX), ΔH of the same sign: the value ±DBL_MAX·(v − u) is finite
        let t = rng.log_uniform(0.31, 10.0);
        let sg = if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let s = sg * (DMAX / t) * rng.uniform(1.0, 1.99);
        gibbs(m, sg * DMAX * rng.uniform(0.2, 1.0), t, s);
        // T·ΔS subnormal, ΔH tiny or zero
        let t = rng.log_uniform(-3.0, 5.0);
        let s = rng.signed_log_uniform(-323.0, -305.0) / t;
        let h = if rng.int(0, 2) == 0 { 0.0 } else { rng.signed_log_uniform(-323.0, -300.0) };
        gibbs(m, h, t, s);
    }
}

/// nuclear's B4 rows: q_value's two sums past DBL_MAX that cancel (Rust's ∞ − ∞ = NaN), one past it against a finite one,
/// the difference past DBL_MAX with the product beyond it (±∞ in both), and subnormal sums (exact: Rust's bits).
pub fn emit_nuclear_b4(m: &mut Module) {
    let mut rng = Rng::new(0x2e33);
    let toks = |r: &[f64], p: &[f64]| {
        let mut ins = seq_toks(r);
        ins.extend(seq_toks(p));
        ins
    };
    let mm1 = f64::from_bits(DMAX.to_bits() - 1);
    for (r, p) in [
        (vec![DMAX, DMAX], vec![DMAX, 0.999 * DMAX]),
        (vec![DMAX, DMAX], vec![DMAX, DMAX]),
        (vec![DMAX, DMAX], vec![DMAX, mm1]),
        (vec![DMAX, DMAX, -DMAX], vec![]),
        (vec![], vec![DMAX, DMAX, -DMAX]),
        (vec![1e308, 1e308], vec![1e308, 1e308, 1e-10]),
        (vec![DMAX, 1e308], vec![DMAX]),
        (vec![sub(3)], vec![sub(1)]),
    ] {
        qv(m, &toks(&r, &p), &r, &p);
    }
    for _ in 0..80 {
        let nr = usize::try_from(rng.int(1, 4)).unwrap_or(1);
        let np = usize::try_from(rng.int(1, 4)).unwrap_or(1);
        let r: Vec<f64> = (0..nr).map(|_| DMAX * rng.uniform(0.5, 1.0)).collect();
        let mut p: Vec<f64> = (0..np).map(|_| DMAX * rng.uniform(0.5, 1.0)).collect();
        // the products' sum close to the reactants' (cancellation of two overflowing sums)
        let sr: f64 = r.iter().map(|x| x / 8.0).sum::<f64>() * 8.0;
        let spp: f64 = p.iter().map(|x| x / 8.0).sum::<f64>() * 8.0;
        if sr.is_finite() && spp.is_finite() && np > 0 {
            p[0] += (sr - spp) * rng.uniform(0.999, 1.0);
        }
        qv(m, &toks(&r, &p), &r, &p);
    }
}

/// solution's B4 rows: ionic_strength's term (c·z)·z or its sum past DBL_MAX while ½·Σ is finite (Rust's ∞), terms of
/// opposite infinities (Rust's NaN, the value 0), the true overflow ends, subnormal terms (exact); kohlrausch's K·√c past
/// DBL_MAX with Λ° bringing it back, K·√c subnormal (Rust's rounding twice); limiting_molar_conductivity's products past
/// DBL_MAX that cancel or that the other brings back, subnormal products.
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_solution_b4(m: &mut Module) {
    let mut rng = Rng::new(0x5033);
    let ions_toks = |xs: &[(f64, i32)]| {
        let mut v = vec![Tok::V(xs.len())];
        for &(c, z) in xs {
            v.push(f(c));
            v.push(i(z));
        }
        v
    };
    for xs in [
        vec![(0.4 * DMAX, 2)],
        vec![(0.4 * DMAX, -2)],
        vec![(DMAX, 1), (DMAX, 1)],
        vec![(DMAX, 1), (DMAX, -1), (-DMAX, 1)],
        vec![(DMAX, 2), (-DMAX, 2)],
        vec![(0.3 * DMAX, 2), (-0.3 * DMAX, 2), (1.0, 1)],
        vec![(0.6 * DMAX, 2)],
        vec![(1e300, 46341)],
        vec![(sub(3), 3), (-sub(1), 1)],
        vec![(DMAX, i32::MAX), (-DMAX, i32::MIN)],
        vec![(0.3 * DMAX, 1), (0.3 * DMAX, 1), (0.3 * DMAX, 1), (0.3 * DMAX, 1)],
    ] {
        ionic(m, &ions_toks(&xs), &xs);
    }
    for _ in 0..80 {
        let n = usize::try_from(rng.int(1, 5)).unwrap_or(1);
        let xs: Vec<(f64, i32)> = (0..n)
            .map(|_| {
                let z = i32::try_from(rng.int(1, 4)).unwrap_or(1) * if rng.int(0, 1) == 0 { 1 } else { -1 };
                (big(&mut rng, 0.05, 0.5) / f64::from(z * z), z)
            })
            .collect();
        ionic(m, &ions_toks(&xs), &xs);
    }
    for (l0, k, c) in [
        (DMAX, DMAX, 1.21),
        (-DMAX, -DMAX, 1.21),
        (DMAX, 1e300, 1e20),
        (1e300, 1e300, 1e10),
        (sub(3), sub(3), 0.25),
        (sub(2), sub(3), 0.25),
        (0.0, sub(1), 0.25),
        (sub(1), 1e-300, 1e-30),
        (0.0, 1e-200, 1e-250),
        (DMAX, DMAX, 4.0),
    ] {
        kohl(m, l0, k, c);
    }
    for _ in 0..80 {
        // K·√c in (DBL_MAX, 2·DBL_MAX), Λ° of the same sign
        let c = rng.log_uniform(0.61, 10.0);
        let sg = if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let k = sg * (DMAX / c.sqrt()) * rng.uniform(1.0, 1.99);
        kohl(m, sg * DMAX * rng.uniform(0.2, 1.0), k, c);
        // K·√c subnormal, Λ° tiny or zero
        let c = rng.log_uniform(-6.0, 6.0);
        let k = rng.signed_log_uniform(-323.0, -306.0) / c.sqrt();
        let l0 = if rng.int(0, 2) == 0 { 0.0 } else { rng.signed_log_uniform(-323.0, -305.0) };
        kohl(m, l0, k, c);
    }
    for (cat, nc, an, na) in [
        (DMAX, 1.5, -DMAX, 1.0),
        (DMAX, 2.0, -DMAX, 1.5),
        (-DMAX, 1.5, DMAX, 1.0),
        (DMAX, 2.0, DMAX, -2.0),
        (DMAX, 2.0, 1.0, 1.0),
        (sub(3), 0.2, sub(3), 0.2),
        (sub(1), 0.75, sub(1), 0.5),
        (1e-200, 1e-200, 1e-200, 1e-200),
        (1e-200, 1e-200, -1e-200, 1e-200),
    ] {
        limiting(m, cat, nc, an, na);
    }
    for _ in 0..80 {
        let (cat, an) = (big(&mut rng, 0.3, 1.0), big(&mut rng, 0.3, 1.0));
        limiting(m, cat, rng.uniform(0.5, 3.0), an, rng.uniform(0.5, 3.0));
        let (cat, an) = (rng.signed_log_uniform(-323.0, -300.0), rng.signed_log_uniform(-323.0, -300.0));
        limiting(m, cat, rng.log_uniform(-10.0, 1.0), an, rng.log_uniform(-10.0, 1.0));
    }
}

/// inorganic's B4 rows: born_haber_lattice_energy's partial differences past ±DBL_MAX that later terms bring back (Rust's
/// ±∞ or NaN), the true ends, and ΔH_diss/2 inexact (an odd subnormal ΔH_diss: Rust's half rounds to even).
pub fn emit_inorganic_b4(m: &mut Module) {
    let mut rng = Rng::new(0x1033);
    for x in [
        [DMAX, -DMAX, DMAX, 0.0, 0.0],
        [-DMAX, DMAX, -DMAX, 0.0, 0.0],
        [DMAX, -DMAX, 0.0, 0.0, DMAX],
        [DMAX, -DMAX, 0.0, DMAX, 0.5 * DMAX],
        [DMAX, -DMAX, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, sub(1), 0.0],
        [sub(1), 0.0, 0.0, sub(5), 0.0],
        [sub(2), 0.0, 0.0, sub(7), -sub(1)],
        [0.0, 0.0, 0.0, sub(3), sub(1)],
        [-DMAX, 0.5 * DMAX, 0.6 * DMAX, -2.0 * 0.9 * DMAX, 0.0],
    ] {
        haber(m, x);
    }
    for _ in 0..80 {
        let x = [big(&mut rng, 0.3, 1.0), big(&mut rng, 0.3, 1.0), big(&mut rng, 0.3, 1.0), big(&mut rng, 0.3, 1.0), big(&mut rng, 0.3, 1.0)];
        haber(m, x);
        let k = u64::try_from(rng.int(0, 1 << 20)).unwrap_or(0) * 2 + 1;
        let y = [sub(u64::try_from(rng.int(0, 1 << 20)).unwrap_or(0)), 0.0, 0.0, sub(k), -sub(u64::try_from(rng.int(0, 8)).unwrap_or(0))];
        haber(m, y);
    }
}

/// kinetics' B4 rows: zero_order_concentration's k·t past DBL_MAX with a negative [A]₀ bringing it back (Rust's +∞), the
/// true ends, k·t subnormal (Rust's rounding twice) and below 2⁻¹⁰⁷⁵ against ±0 (the sign of the clamp's zero).
pub fn emit_kinetics_b4(m: &mut Module) {
    let mut rng = Rng::new(0x4b33);
    for (c0, k, t) in [
        (-DMAX, -DMAX, 1.5),
        (DMAX, DMAX, 1.5),
        (-DMAX, DMAX, -1.5),
        (sub(3), 1.5, sub(1)),
        (sub(1), 0.75, sub(1)),
        (sub(5), 0.5, sub(5)),
        (0.0, 1e-200, -1e-200),
        (-0.0, 1e-200, 1e-200),
        (-0.0, 1e-200, -1e-200),
        (0.0, 1e-200, 1e-200),
        (-DMAX, -2.0, DMAX),
    ] {
        zero_conc(m, c0, k, t);
    }
    for _ in 0..80 {
        // k·t in (−2·DBL_MAX, −DBL_MAX), [A]₀ negative: the value DBL_MAX·(u − v) can be positive and finite
        let t = rng.log_uniform(0.31, 10.0);
        let k = -(DMAX / t) * rng.uniform(1.0, 1.99);
        zero_conc(m, -DMAX * rng.uniform(0.2, 1.0), k, t);
        let t = rng.log_uniform(-3.0, 5.0);
        let k = rng.signed_log_uniform(-323.0, -305.0) / t;
        let c0 = if rng.int(0, 2) == 0 { 0.0 } else { rng.signed_log_uniform(-323.0, -305.0) };
        zero_conc(m, c0, k, t);
    }
}

/// A reaction side, (formula, coefficient) pairs.
type Side<'a> = Vec<(&'a str, f64)>;

/// thermochem's B4 rows: the reaction sums' products past DBL_MAX whose sums cancel (Rust's ∞ − ∞ = NaN), subnormal
/// products coeff·prop (Rust's subnormal rounding, then × 1000); and Shomate polynomials with user coefficients at any
/// T: T → 0 (t·t underflowing: E/(t·t) = ±∞ or 0/0 = NaN where Cp is finite, the antiderivative's t subnormal or 0),
/// large T (C·t², D·t³ and t⁴ overflowing: ∞ − ∞ = NaN where the value is ±∞, and an exact cancellation of two such
/// terms at t = 2⁶⁰⁰), subnormal coefficients, and seeded coefficient sets over every binade.
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_thermochem(m: &mut Module) {
    let mut rng = Rng::new(0x3c33);
    let sp_toks = |p: &[(&str, f64)], r: &[(&str, f64)]| {
        let mut v = vec![Tok::V(p.len())];
        for (fo, c) in p {
            v.push(s(*fo));
            v.push(f(*c));
        }
        v.push(Tok::V(r.len()));
        for (fo, c) in r {
            v.push(s(*fo));
            v.push(f(*c));
        }
        v
    };
    let a2 = "Al2O3(s)";
    let cases: Vec<(Side, Side)> = vec![
        (vec![(a2, 1e305), (a2, 1e305)], vec![(a2, 1e305), (a2, 0.999e305)]),
        (vec![(a2, 1e306), (a2, -1e306)], vec![]),
        (vec![(a2, 2e305), ("CO2(g)", 1.0)], vec![(a2, 2e305)]),
        (vec![("CO2(g)", 5e-324)], vec![]),
        (vec![("CO2(g)", 1e-320)], vec![("H2O(l)", 1e-321)]),
        (vec![("H2(g)", 1.0), ("CO2(g)", sub(3))], vec![("H2(g)", 1.0)]),
        (vec![("CO2(g)", 4e306), ("CO2(g)", 4e306)], vec![("CO2(g)", 4e306), ("CO2(g)", 3.99e306)]),
        (vec![("CH4(g)", 1e-310)], vec![("CO2(g)", 1e-311)]),
    ];
    for (p, r) in &cases {
        sums(m, &sp_toks(p, r), p, r);
    }
    let formulas: Vec<&str> = THERMOCHEM_DATA.iter().map(|d| d.formula).collect();
    let pick = |rng: &mut Rng| formulas[usize::try_from(rng.int(0, i64::try_from(formulas.len()).unwrap_or(1) - 1)).unwrap_or(0)];
    for _ in 0..80 {
        // products past DBL_MAX: coefficients near DBL_MAX/|prop|, the two sides' sums close
        let f1 = pick(&mut rng);
        let pr = lookup_thermochem(f1).map_or(1.0, |d| d.delta_hf_kj.abs().max(1.0));
        let c1 = DMAX / pr * rng.uniform(0.5, 1.0);
        let c2 = c1 * rng.uniform(0.999, 1.001);
        let p = vec![(f1, c1), (f1, c1)];
        let r = vec![(f1, c2), (f1, c1)];
        sums(m, &sp_toks(&p, &r), &p, &r);
        // subnormal products
        let p = vec![(pick(&mut rng), rng.signed_log_uniform(-323.0, -308.0))];
        let r = vec![(pick(&mut rng), rng.signed_log_uniform(-323.0, -308.0))];
        sums(m, &sp_toks(&p, &r), &p, &r);
    }

    let sh = |m: &mut Module, x: &Shomate, t: f64| {
        let mut ins: Vec<Tok> = shomate_fields(x).iter().map(|&v| f(v)).collect();
        ins.push(f(t));
        shomate(m, &ins, x, t);
    };
    let mk = |a: f64, b: f64, c: f64, d: f64, e: f64| Shomate { a, b, c, d, e, t_min: 298.0, t_max: 6000.0 };
    let t600 = 125.0 * pow2(603); // T/1000 = 2⁶⁰⁰ exactly
    let sets = [
        mk(20.786, 0.0, 0.0, 0.0, 0.0),
        mk(1.0, 0.0, 0.0, 0.0, 1e-310),
        mk(1.0, 2.0, 3.0, 4.0, 1e-300),
        mk(0.0, 0.0, 1.0, -1.0, 0.0),
        mk(0.0, 0.0, 1.0, -pow2(-600), 0.0),
        mk(5.0, 0.0, 1.0, -pow2(-600), 0.0),
        mk(1e300, -1e300, 1e300, -1e300, 0.0),
        mk(sub(3), sub(1), 0.0, 0.0, 0.0),
        mk(30.092, 6.832514, 6.793435, -2.53448, 0.082139),
        mk(0.0, 0.0, 0.0, 1.0, 0.0),
        mk(1e-300, 1e-300, 1e-300, 1e-300, 1e-300),
    ];
    for x in &sets {
        for t in [1e-300, 5e-324, 1e-160, 1e-163, 1e-155, 2.2e-305, 1e200, 1e300, -1e300, DMAX, t600, -t600, 1e-10, 1e10] {
            sh(m, x, t);
        }
    }
    for _ in 0..120 {
        let co = |rng: &mut Rng| if rng.int(0, 4) == 0 { 0.0 } else { rng.signed_log_uniform(-320.0, 308.0) };
        let x = mk(co(&mut rng), co(&mut rng), co(&mut rng), co(&mut rng), co(&mut rng));
        let t = rng.signed_log_uniform(-323.0, 308.0);
        sh(m, &x, t);
        // physical coefficients, extreme T
        let y = mk(rng.uniform(-50.0, 50.0), rng.uniform(-50.0, 50.0), rng.uniform(-20.0, 20.0), rng.uniform(-5.0, 5.0),
                   rng.uniform(-2.0, 2.0));
        let t = if rng.int(0, 1) == 0 { rng.log_uniform(-323.0, -150.0) } else { rng.signed_log_uniform(75.0, 308.0) };
        sh(m, &y, t);
    }
}

// ── The overflow audit's B5 (ADR 0002 #35): the exp / exp10 arguments ───────────────────────────────────────────────

/// How the port obtains `equilibrium_constant`'s K (the note's `<how>`).
pub const KEQ_HOW: &str =
    "the argument evaluated with an unbounded exponent (its final quotient rounded once), then the correctly rounded exp";

/// How the port obtains the limiting law's and Davies' γ (the note's `<how>`).
pub const DH_HOW_MUL: &str =
    "the argument evaluated with an unbounded exponent (its final product rounded once), then the correctly rounded 10^x";

/// How the port obtains the extended law's γ (the note's `<how>`).
pub const DH_HOW_DIV: &str =
    "the argument evaluated with an unbounded exponent (its final quotient rounded once), then the correctly rounded 10^x";

/// ADR 0002 #35: `equilibrium_constant`, Rust's checks in Rust's order, then x = −ΔG/(R·T) with R·T rounded to 53 bits
/// and the quotient rounded once into binary64, K = e^x (glibc's exp, which `validate.py` checks is correctly rounded on
/// every noted row: the port's `_km_cr_exp`), and Rust's ComputationError where that K is not finite.
pub fn port_keq(dg: f64, t: f64) -> Port {
    if !pos(t) {
        return Err("InvalidTemperature");
    }
    if !dg.is_finite() {
        return Err("InvalidInput");
    }
    let x = ux(dg).neg().div_f64(ux(GAS_CONSTANT).mul(ux(t)));
    let k = bb(x).exp();
    if k.is_finite() { Ok(k) } else { Err("ComputationError") }
}

/// Rust's A·z² for the activity models, ((−A)·z)·z, rounded to 53 bits at each product.
fn az2(z: i32, a: f64) -> Ux {
    let zf = ux(f64::from(z));
    ux(a).neg().mul(zf).mul(zf)
}

/// 10^x as Rust writes it (`10.0_f64.powf(x)`: glibc's, checked against MPFR's exp10 by `validate.py` on every noted row).
fn exp10(x: f64) -> f64 {
    10.0_f64.powf(bb(x))
}

/// ADR 0002 #35: `debye_huckel_limiting`, Rust's checks, then log₁₀ γ = (((−A)·z)·z)·√I with the product A·z² and √I
/// rounded to 53 bits and the final product rounded once into binary64, γ = 10^x.
pub fn port_dh_lim(z: i32, ii: f64, a: f64) -> Option<f64> {
    (nonneg(ii) && a.is_finite()).then(|| exp10(az2(z, a).mul_f64(ux(ii).sqrt())))
}

/// ADR 0002 #35: `debye_huckel_extended`, Rust's checks in Rust's order, then log₁₀ γ =
/// ((((−A)·z)·z)·√I)/(1 + (B·å)·√I) with every operation rounded to 53 bits but the final quotient, rounded once.
pub fn port_dh_ext(z: i32, ii: f64, a: f64, b: f64, s: f64) -> Option<f64> {
    (nonneg(ii) && a.is_finite() && nonneg(b) && nonneg(s)).then(|| {
        let sq = ux(ii).sqrt();
        let num = az2(z, a).mul(sq);
        let den = ux(1.0).add(ux(b).mul(ux(s)).mul(sq));
        exp10(num.div_f64(den))
    })
}

/// ADR 0002 #35: `davies_activity`, Rust's checks, then log₁₀ γ = (((−A)·z)·z)·(√I/(1 + √I) − 0.3·I) with every
/// operation rounded to 53 bits but the final product, rounded once.
pub fn port_davies(z: i32, ii: f64, a: f64) -> Option<f64> {
    (nonneg(ii) && a.is_finite()).then(|| {
        let sq = ux(ii).sqrt();
        let br = sq.div(ux(1.0).add(sq)).sub(ux(0.3).mul(ux(ii)));
        exp10(az2(z, a).mul_f64(br))
    })
}

/// A row of `equilibrium_constant` under B5: plain where the port's result is the oracle's (the same bits, or the same
/// Err), noted otherwise — `= <x> (f:<bits>)` for the port's K, `= Err(ComputationError) (E:ComputationError)` where the
/// port's K overflows and Rust's does not. The validation is Rust's, so only ComputationError may differ (asserted).
#[allow(clippy::too_many_arguments)] // noted_r's
pub fn noted_k(
    m: &mut Module,
    name: &str,
    cls: Class,
    ins: &[Tok],
    oracle: Result<f64, KimiyaError>,
    port: Port,
    adr: u32,
    how: &str,
) {
    match (&oracle, port) {
        (Ok(v), Ok(p)) if same(*v, p) => m.row(name, cls, ins, &[f(*v)]),
        (Err(e), Err(n)) => {
            assert_eq!(err_name(e), n, "{name}: the port's Err is not the oracle's: {ins:?}");
            m.row(name, cls, ins, &[res(oracle)]);
        }
        (Ok(_), Ok(p)) => m.row_noted(name, cls, ins, &[res(oracle)], &note(adr, how, p)),
        (Err(e), Ok(p)) => {
            assert_eq!(err_name(e), "ComputationError", "{name}: the port is Ok where Rust's validation refuses: {ins:?}");
            m.row_noted(name, cls, ins, &[res(oracle)], &note(adr, how, p));
        }
        (Ok(_), Err(n)) => {
            assert_eq!(n, "ComputationError", "{name}: the port's validation refuses where the oracle's accepts: {ins:?}");
            let extra = format!("ADR 0002 #{adr}: the port returns {how} = Err(ComputationError) (E:ComputationError)");
            m.row_noted(name, cls, ins, &[res(oracle)], &extra);
        }
    }
}

pub fn keq(m: &mut Module, dg: f64, t: f64) {
    let out = equilibrium_constant(bb(dg), bb(t));
    noted_k(m, "reaction::equilibrium_constant", Class::C, &[f(dg), f(t)], out, port_keq(dg, t), 35, KEQ_HOW);
}

pub fn dh_lim(m: &mut Module, z: i32, ii: f64, a: f64) {
    let out = debye_huckel_limiting(bb(z), bb(ii), bb(a));
    let name = "solution::debye_huckel_limiting";
    noted_cls(m, name, Class::C, &[i(z), f(ii), f(a)], out, port_dh_lim(z, ii, a), 35, DH_HOW_MUL);
}

pub fn dh_ext(m: &mut Module, z: i32, ii: f64, a: f64, b: f64, s: f64) {
    let out = debye_huckel_extended(bb(z), bb(ii), bb(a), bb(b), bb(s));
    let (name, ins) = ("solution::debye_huckel_extended", [i(z), f(ii), f(a), f(b), f(s)]);
    noted_cls(m, name, Class::C, &ins, out, port_dh_ext(z, ii, a, b, s), 35, DH_HOW_DIV);
}

pub fn davies(m: &mut Module, z: i32, ii: f64, a: f64) {
    let out = davies_activity(bb(z), bb(ii), bb(a));
    let name = "solution::davies_activity";
    noted_cls(m, name, Class::C, &[i(z), f(ii), f(a)], out, port_davies(z, ii, a), 35, DH_HOW_MUL);
}

// ── B5's edge rows (appended to each module with their own seed, so the rows before keep theirs) ─────────────────

/// T = k·2⁻¹⁰⁷⁴ with R·k = 2⁵² − 0.27…: R·T just below 2⁻¹⁰²², which IEEE rounds up to 2⁻¹⁰²² at the subnormal grid and
/// 53-bit rounding to 2⁻¹⁰²² − 2⁻¹⁰⁷⁵ (the B4 pitfall of #29: a fast path must not take 2⁻¹⁰²² as normal).
pub const KEQ_T_2M1022: u64 = 541_658_533_363_376;

/// The two doubles I near 1.9404 where √I/(1 + √I) and 0.3·I round to the same double: Davies' bracket exactly 0.
pub const DAVIES_ZERO_I: [f64; 2] = [1.940_363_884_733_242, 1.940_363_884_733_242_5];

/// A charge for the activity models' edges: |z| ≥ 2 where `big` (A·z² can pass DBL_MAX), any i32 otherwise.
fn charge(rng: &mut Rng, big: bool) -> i32 {
    let mag = match rng.int(0, 3) {
        0 => rng.int(2, 4),
        1 => rng.int(5, 1000),
        2 => rng.int(1001, 100_000),
        _ => rng.int(100_001, i64::from(i32::MAX)),
    };
    let mag = if big { mag } else { mag - 1 };
    let z = i32::try_from(mag).unwrap_or(2);
    if rng.int(0, 1) == 0 { z } else { -z }
}

/// reaction's B5 rows: equilibrium_constant's R·T past DBL_MAX (T > 2.16·10³⁰⁷ K: Rust's x = ∓0, K = 1, where e^x with
/// |x| < 1/R is the value) and subnormal (T < 2⁻¹⁰²¹/R: Rust's degraded x, K off, or Ok where the port's K overflows
/// (Err) and the reverse), R·T at 2⁻¹⁰²² from below, subnormal and huge ΔG, x past ±DBL_MAX (K = ∞: Err, K = 0 in both),
/// x subnormal (K = 1 in both), the boundaries T = DBL_MAX/R and 2⁻¹⁰²¹/R, and every binade of both inputs.
pub fn emit_reaction_b5(m: &mut Module) {
    let mut rng = Rng::new(0x6b35);
    let tk = sub(KEQ_T_2M1022);
    for (dg, t) in [
        (-DMAX, DMAX),
        (DMAX, DMAX),
        (-1e308, DMAX),
        (-DMAX, 1e308),
        (1e308, 1e308),
        (-DMAX, 2.2e307),
        (-1e-320, 1e-322),
        (1e-320, 1e-322),
        (-1e-320, sub(1)),
        (0.0, DMAX),
        (-0.0, DMAX),
        (0.0, sub(1)),
        (-0.0, sub(1)),
        (-DMAX, sub(1)),
        (DMAX, sub(1)),
        (-1.0, sub(1)),
        (sub(1), DMAX),
        (-sub(1), 1e-300),
        (-sub(3), 1.0),
        (-700.0 * pow2(-1022), tk),
        (-pow2(-1022), tk),
        (pow2(-1022), tk),
        (-300.0 * pow2(-1022), tk),
        (-3.0 * pow2(-1022), tk),
    ] {
        keq(m, dg, t);
    }
    // the boundaries: T around DBL_MAX/R (R·T at DBL_MAX) and around 2⁻¹⁰²¹/R (R·T at 2⁻¹⁰²¹), ulp by ulp
    let t_ov = DMAX / GAS_CONSTANT;
    let t_un = pow2(-1021) / GAS_CONSTANT;
    for k in -3_i64..=3 {
        let tv = f64::from_bits(t_ov.to_bits().wrapping_add_signed(k));
        keq(m, -DMAX, tv);
        keq(m, -0.5 * DMAX, tv);
        let tu = f64::from_bits(t_un.to_bits().wrapping_add_signed(k));
        for x in [1.0, 100.0, 700.0, -700.0] {
            keq(m, -x * GAS_CONSTANT * tu, tu);
        }
    }
    // R·T at 2⁻¹⁰²² from below, and its neighbours
    for d in 0_u64..6 {
        let t = sub(KEQ_T_2M1022 - 2 + d);
        for x in [0.5, 2.0, 50.0, 709.0, -744.0] {
            keq(m, -x * pow2(-1022), t);
        }
    }
    for _ in 0..80 {
        // R·T past DBL_MAX: |x| < 1/R, K in (e^(−1/R), e^(1/R)), Rust's K = 1
        let t = DMAX * rng.uniform(0.121, 1.0);
        let sg = if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        keq(m, sg * DMAX * rng.uniform(1e-3, 1.0), t);
        keq(m, rng.signed_log_uniform(-323.0, 308.0), t);
        // R·T subnormal: ΔG = −x·R·T for a target x across exp's range
        let t = rng.log_uniform(-323.3, -308.3);
        let x = rng.uniform(-750.0, 715.0);
        keq(m, -x * GAS_CONSTANT * t, t);
        // near the overflow and underflow ends at a subnormal T: Rust's K Ok where the port's overflows, and the reverse
        let t = rng.log_uniform(-323.3, -308.3);
        let x = if rng.int(0, 1) == 0 { rng.uniform(709.0, 710.5) } else { rng.uniform(-746.0, -744.0) };
        keq(m, -x * GAS_CONSTANT * t, t);
    }
    // a short subnormal T = k·2⁻¹⁰⁷⁴ and ΔG = −n·2⁻¹⁰⁷⁴ between the two ends' thresholds as Rust's R·T (k·R rounded to an
    // integer) and the port's (rounded to 53 bits) place them: K overflows (Err) in one and not the other, or is 0 in
    // one and subnormal in the other
    for _ in 0..60 {
        let k = u64::try_from(rng.int(1, 1 << 24)).unwrap_or(1);
        #[allow(clippy::cast_precision_loss)] // k < 2²⁴, exact
        let rk53 = GAS_CONSTANT * k as f64;
        #[allow(clippy::cast_precision_loss)] // the subnormal's integer, below 2²⁸: exact
        let rki = (GAS_CONSTANT * sub(k)).to_bits() as f64;
        let mid = 0.5 * (rk53 + rki);
        for (xe, sg) in [(709.782_712_893_384, -1.0), (745.133_219_101_941_1, 1.0)] {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)] // below 2⁴⁰, non-negative
            let n = (xe * mid).floor() as u64 + u64::try_from(rng.int(0, 1)).unwrap_or(0);
            keq(m, sg * sub(n), sub(k));
        }
    }
    for _ in 0..40 {
        // subnormal T against a huge or subnormal ΔG: x past ±DBL_MAX (K = ∞: Err, or 0), or ordinary
        let t = rng.log_uniform(-323.3, -308.3);
        keq(m, rng.signed_log_uniform(-323.0, 308.0), t);
        // x subnormal: K = 1 in both
        keq(m, rng.signed_log_uniform(-323.3, -300.0), rng.log_uniform(1.0, 300.0));
    }
    for _ in 0..60 {
        // every binade of both inputs
        keq(m, rng.signed_log_uniform(-323.3, 308.2), rng.log_uniform(-323.3, 308.2));
    }
}

/// solution's B5 rows: the limiting law's and Davies' A·z² past DBL_MAX at I = ±0 (Rust's ∞·0 = NaN, γ = 1) and with √I
/// or the bracket bringing it back (beyond ±2⁴⁸⁴: Rust's γ = ∞ or 0 in both), A·z² at DBL_MAX ulp by ulp, A·z²
/// subnormal (exact in both), Davies' bracket exactly 0 (I = 1.9404…, Rust's NaN against an overflowing A·z²) and
/// negative; the extended law's B·å and (B·å)·√I past DBL_MAX (Rust's γ = 1 where log₁₀ γ is finite, or NaN at I = 0),
/// A·z² and the numerator past DBL_MAX with the denominator bringing them back (Rust's γ = 0 or ∞, or ∞/∞ = NaN), the
/// numerator and B·å subnormal (γ = 1, or the denominator 1, in both), B and å ±0; and every binade of every input.
#[allow(clippy::too_many_lines)] // one function after another
pub fn emit_solution_b5(m: &mut Module) {
    let mut rng = Rng::new(0x5035);
    for (z, ii, a) in [
        (100_000, 0.0, 1e300),
        (100_000, -0.0, 1e300),
        (100_000, 0.0, -1e300),
        (i32::MAX, 0.0, 1e290),
        (i32::MIN, -0.0, -DMAX),
        (1, 0.0, DMAX),
        (-1, 0.0, -DMAX),
        (2, 0.0, DMAX),
        (-3, 0.0, -0.25 * DMAX),
        (2, 1e-300, DMAX),
        (2, sub(1), DMAX),
        (2, 1.0, DMAX),
        (i32::MAX, sub(1), DMAX),
        (1, 0.01, sub(3)),
        (3, 0.5, sub(1)),
        (46_341, 1.0, 1e-320),
        (0, 0.0, DMAX),
        (0, 1.0, DMAX),
    ] {
        dh_lim(m, z, ii, a);
        davies(m, z, ii, a);
    }
    for ii in DAVIES_ZERO_I {
        for (z, a) in [(2, DMAX), (-2, -DMAX), (100_000, 1e300), (1, 0.509), (i32::MIN, 1.0)] {
            davies(m, z, ii, a);
            dh_lim(m, z, ii, a);
        }
    }
    // A·z² at DBL_MAX, ulp by ulp, at I = 0 (Rust's NaN only past it)
    for z in [2_i32, 3, 46_341, 100_000] {
        let zf = f64::from(z);
        let a0 = DMAX / zf / zf;
        for k in -2_i64..=2 {
            let a = f64::from_bits(a0.to_bits().wrapping_add_signed(k));
            dh_lim(m, z, 0.0, a);
            davies(m, -z, -0.0, -a);
        }
    }
    for _ in 0..80 {
        // A·z² past DBL_MAX: at I = ±0 (Rust's NaN), past it (the argument beyond ±2⁴⁸⁴: Rust's γ in both)
        let z = charge(&mut rng, true);
        let zf = f64::from(z);
        let lo = (DMAX / zf / zf).min(DMAX);
        let sg = if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let a = sg * (lo * rng.uniform(1.0001, 4.0)).min(DMAX);
        let ii = match rng.int(0, 3) {
            0 => 0.0,
            1 => -0.0,
            2 => DAVIES_ZERO_I[usize::try_from(rng.int(0, 1)).unwrap_or(0)],
            _ => rng.log_uniform(-323.3, 308.2),
        };
        dh_lim(m, z, ii, a);
        davies(m, z, ii, a);
        // Davies' bracket negative (I past 1.94) or tiny, against an ordinary or overflowing A·z²
        let z = charge(&mut rng, false);
        let ii = if rng.int(0, 1) == 0 { rng.uniform(1.9, 2.0) } else { rng.log_uniform(0.3, 308.2) };
        davies(m, z, ii, rng.signed_log_uniform(-1.0, 308.2));
        // A·z² subnormal (exact in both)
        let z = charge(&mut rng, false);
        dh_lim(m, z, rng.log_uniform(-10.0, 10.0), rng.signed_log_uniform(-323.3, -310.0));
        davies(m, z, rng.log_uniform(-10.0, 10.0), rng.signed_log_uniform(-323.3, -310.0));
    }
    for _ in 0..60 {
        // every binade of every input
        let z = charge(&mut rng, false);
        let (ii, a) = (rng.log_uniform(-323.3, 308.2), rng.signed_log_uniform(-323.3, 308.2));
        dh_lim(m, z, ii, a);
        davies(m, z, ii, a);
    }

    for (z, ii, a, b, s) in [
        (1, 1.0, 1e308, 1e155, 1e155),
        (100_000, 1.0, 1e300, 1e154, 1e154),
        (1, 0.0, 1.0, DMAX, 2.0),
        (1, -0.0, 0.509, 1e200, 1e200),
        (1, 0.0, 0.509, 1e200, 1e200),
        (1, 1e20, 1e300, 1e298, 1.0),
        (2, 1.0, DMAX, DMAX, DMAX),
        (1, 1e-300, 1e300, 1e300, 1e300),
        (1, 1.0, 1e300, 1e155, 1e155),
        (100_000, 0.0, 1e300, 0.328, 3.0),
        (100_000, -0.0, -1e300, 1e200, 1e200),
        (1, 1.0, 1e308, 1e300, 1e300),
        (1, 4.0, DMAX, DMAX, 1.0),
        (1, 0.01, 0.509, -0.0, 3.0),
        (1, 0.01, 0.509, 0.328, -0.0),
        (1, -0.0, 0.509, -0.0, -0.0),
        (1, 0.01, sub(1), 0.328, 3.0),
        (2, sub(1), 1.0, 0.328, 3.0),
        (1, 0.01, 0.509, sub(1), sub(1)),
        (1, 1e300, 0.509, sub(5), 1e-300),
        (0, 0.0, DMAX, DMAX, DMAX),
        (0, 1.0, DMAX, DMAX, DMAX),
    ] {
        dh_ext(m, z, ii, a, b, s);
    }
    for _ in 0..80 {
        // B·å past DBL_MAX, the numerator finite: log₁₀ γ ≈ ∓A·z²/(B·å) = ∓L (Rust's 10^(x/∞) = 1): A·z² = P in
        // [10⁻⁵·DBL_MAX, DBL_MAX), L below P/DBL_MAX (so B·å = P/L passes DBL_MAX), split as r·u and r/u, √I ≤ 1
        let z = charge(&mut rng, false);
        let zf = f64::from(z);
        let p = DMAX * rng.log_uniform(-5.0, 0.0);
        let a = p / zf / zf * if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let l = p / DMAX * rng.log_uniform(-6.0, 0.0);
        let r = p.sqrt() / l.sqrt();
        let u = rng.log_uniform(-3.0, 3.0);
        dh_ext(m, z, rng.log_uniform(-100.0, 0.0), a, r * u, r / u);
        // A·z² past DBL_MAX by a factor w below 100, the denominator finite and bringing it back (Rust's 10^(±∞/y) =
        // ∞ or 0): B·å = |A|·z²/L = DBL_MAX·w/L with L over (w, 320), and √I in [1, L/w) so that (B·å)·√I is finite
        let z = i32::try_from(rng.int(2, 10)).unwrap_or(2) * if rng.int(0, 1) == 0 { 1 } else { -1 };
        let zf = f64::from(z);
        let w = rng.uniform(1.01, 100.0);
        let a = (DMAX / zf / zf * w).min(DMAX) * if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let w = a.abs() / DMAX * zf * zf;
        let l = rng.uniform(w * 1.01, 320.0);
        let sq = rng.uniform(1.0, l / (1.01 * w));
        let r = DMAX.sqrt() * (w / l).sqrt();
        let u = rng.log_uniform(-3.0, 3.0);
        dh_ext(m, z, sq * sq, a, r * u, r / u);
        // A·z² past DBL_MAX, the denominator bringing it back: B·å·√I ≈ |A|·z²/L
        let z = charge(&mut rng, true);
        let zf = f64::from(z);
        let a = DMAX * rng.uniform(0.3, 1.0) * if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let ii = rng.log_uniform(-50.0, 50.0);
        let l = rng.log_uniform(-3.0, 2.5);
        let r = (a.abs().sqrt() * zf.abs()) / (l * ii.sqrt()).sqrt();
        let u = rng.log_uniform(-3.0, 3.0);
        if (r * u).is_finite() && (r / u).is_finite() {
            dh_ext(m, z, ii, a, r * u, r / u);
        }
        // the numerator past DBL_MAX with A·z² finite: (A·z²)·√I overflows, B·å·√I brings it back
        let z = charge(&mut rng, false);
        let zf = f64::from(z);
        let a = (DMAX / zf / zf) * rng.uniform(0.01, 1.0) * if rng.int(0, 1) == 0 { 1.0 } else { -1.0 };
        let ii = rng.log_uniform(1.0, 300.0);
        let b = rng.log_uniform(-1.0, 308.0);
        dh_ext(m, z, ii, a, b, rng.log_uniform(-1.0, 308.0));
        // I = ±0 with B·å or A·z² past DBL_MAX (Rust's ∞·0 = NaN, γ = 1)
        let z = charge(&mut rng, true);
        let ii = if rng.int(0, 1) == 0 { 0.0 } else { -0.0 };
        let a = rng.signed_log_uniform(250.0, 308.2);
        dh_ext(m, z, ii, a, rng.log_uniform(154.5, 308.2), rng.log_uniform(154.5, 308.2));
    }
    for _ in 0..40 {
        // the numerator or B·å subnormal (γ = 1, or the denominator 1, in both)
        let z = charge(&mut rng, false);
        let ii = rng.log_uniform(-323.3, 0.0);
        dh_ext(m, z, ii, rng.signed_log_uniform(-323.3, -300.0), rng.log_uniform(-10.0, 1.0), rng.log_uniform(0.0, 1.0));
        dh_ext(m, z, rng.log_uniform(-5.0, 300.0), 0.509, rng.log_uniform(-323.3, -310.0), rng.log_uniform(-10.0, 0.0));
    }
    for _ in 0..60 {
        // every binade of every input
        let z = charge(&mut rng, false);
        let ii = rng.log_uniform(-323.3, 308.2);
        let a = rng.signed_log_uniform(-323.3, 308.2);
        dh_ext(m, z, ii, a, rng.log_uniform(-323.3, 308.2), rng.log_uniform(-323.3, 308.2));
    }
    // the argument's final product a 53-bit midpoint (A·z²·√I = (2⁵³ + 1)·2^k with √I = 3·2^j, A = y·2^i, 3y = 2⁵³ + 1):
    // in range, Rust's own tie (rounded to even in both evaluations); past DBL_MAX, the slow path's tie (γ = ∞ or 0)
    let y = 3_002_399_751_580_331.0;
    for (z, i, j) in [(1, -50, 0), (1, -52, 1), (-1, -49, -1), (2, -54, 2), (2, 971, -400), (-2, 970, -500)] {
        let (a, sq) = (y * pow2(i), 3.0 * pow2(j));
        dh_lim(m, z, sq * sq, a);
        dh_lim(m, z, sq * sq, -a);
    }
}

// ── The overflow audit's B6 (ADR 0002 #36, #37): the FFT and find_peaks' threshold ──────────────────────────────────

/// How the port obtains `find_peaks`' result (the note's `<how>`).
pub const PEAKS_HOW: &str = "the threshold test evaluated exactly (w > min·(1 − t) + max·t in rational arithmetic)";

/// ADR 0002 #37: `find_peaks` with Rust's checks, Rust's fold for min and max over the finite values, the neighbour
/// tests as IEEE, and the threshold test exact ([`exact::peak_above`]); a non-finite candidate compared as IEEE would
/// against a finite threshold (+∞ above it, NaN never).
pub fn port_peaks(data: &[f64], t: f64) -> Vec<usize> {
    if data.len() < 3 || !(0.0..=1.0).contains(&t) {
        return Vec::new();
    }
    let (lo, hi) = data
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)));
    if lo > hi {
        return Vec::new();
    }
    data.windows(3)
        .enumerate()
        .filter(|(_, w)| {
            let above = if w[1].is_finite() { crate::exact::peak_above(w[1], lo, hi, t) } else { w[1] > 0.0 };
            w[1] > w[0] && w[1] > w[2] && above
        })
        .map(|(i, _)| i + 1)
        .collect()
}

/// How the port obtains `fft_spectrum`'s result (the note's `<how>`).
pub const FFT_HOW: &str =
    "hisab 1.4.0's transform and the magnitudes evaluated with an unbounded exponent (each √ rounded once)";

/// ADR 0002 #36: `fft_spectrum` of a finite signal (a non-power-of-two length, or a NaN or ±∞ sample, keeps Rust's
/// result: `None`): hisab 1.4.0's `num::fft` — the bit reversal, per stage wn = (cos θ, sin θ) of θ = (−2·π)/len (glibc's
/// cos and sin, the oracle's, which `validate.py` checks are correctly rounded at every FFT angle: the port's
/// `_km_cr_sincos`), w = (1, 0), t = w·d (re = w.re·d.re − w.im·d.im, im = w.re·d.im + w.im·d.re), u ± t, w *= wn — on
/// [`Ux`], every operation rounded to 53 bits with an unbounded exponent, then Complex::abs's √(re·re + im·im) of bins
/// 0 … n/2 with the squares and their sum on [`Ux`] and the √ rounded once into binary64.
pub fn port_fft(signal: &[f64]) -> Option<Vec<f64>> {
    let n = signal.len();
    if n == 0 || !n.is_power_of_two() || signal.iter().any(|x| !x.is_finite()) {
        return None;
    }
    let mut re: Vec<Ux> = signal.iter().map(|&x| ux(x)).collect();
    let mut im: Vec<Ux> = vec![ux(0.0); n];
    if n > 1 {
        let mut j = 0_usize;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;
            if i < j {
                re.swap(i, j);
                im.swap(i, j);
            }
        }
        let mut len = 2_usize;
        while len <= n {
            let half = len / 2;
            #[allow(clippy::cast_precision_loss)] // a power of two, exact
            let angle = -2.0 * std::f64::consts::PI / len as f64;
            let (wnr, wni) = (ux(bb(angle).cos()), ux(bb(angle).sin()));
            let mut start = 0;
            while start < n {
                let (mut wr, mut wi) = (ux(1.0), ux(0.0));
                for k in 0..half {
                    let (a, b) = (start + k, start + k + half);
                    let tr = wr.mul(re[b]).sub(wi.mul(im[b]));
                    let ti = wr.mul(im[b]).add(wi.mul(re[b]));
                    let (ur, ui) = (re[a], im[a]);
                    re[a] = ur.add(tr);
                    im[a] = ui.add(ti);
                    re[b] = ur.sub(tr);
                    im[b] = ui.sub(ti);
                    let nr = wr.mul(wnr).sub(wi.mul(wni));
                    let ni = wr.mul(wni).add(wi.mul(wnr));
                    wr = nr;
                    wi = ni;
                }
                start += len;
            }
            len <<= 1;
        }
    }
    Some((0..=n / 2).map(|k| re[k].mul(re[k]).add(im[k].mul(im[k])).sqrt_f64()).collect())
}
