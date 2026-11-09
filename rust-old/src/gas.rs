//! Gas laws — ideal gas, van der Waals, Peng–Robinson and Redlich–Kwong
//! equations of state, Dalton's law, gas density and PR fugacity.
//!
//! All quantities are SI: pressure in Pa, temperature in K, volume in m³
//! (molar volume in m³/mol), amount in mol. Fallible functions reject NaN and
//! ±∞ inputs with an error instead of returning `Ok(NaN)`.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};

/// Ideal gas law: P = nRT/V
///
/// n = moles, T = temperature (K), V = volume (m³). Returns pressure in Pascals.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if moles is not finite or volume is
/// not positive and finite, or [`KimiyaError::InvalidTemperature`] if
/// temperature is not positive and finite.
#[inline]
pub fn ideal_gas_pressure(moles: f64, temperature_k: f64, volume_m3: f64) -> Result<f64> {
    if !is_finite(moles) {
        return Err(KimiyaError::InvalidInput("moles must be finite".into()));
    }
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_positive(volume_m3) {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    Ok(moles * GAS_CONSTANT * temperature_k / volume_m3)
}

/// Ideal gas law: V = nRT/P
///
/// n = moles, T = temperature (K), P = pressure (Pa). Returns volume in m³.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if moles is not finite or pressure is
/// not positive and finite, or [`KimiyaError::InvalidTemperature`] if
/// temperature is not positive and finite.
#[inline]
pub fn ideal_gas_volume(moles: f64, temperature_k: f64, pressure_pa: f64) -> Result<f64> {
    if !is_finite(moles) {
        return Err(KimiyaError::InvalidInput("moles must be finite".into()));
    }
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_positive(pressure_pa) {
        return Err(KimiyaError::InvalidInput(
            "pressure must be positive".into(),
        ));
    }
    Ok(moles * GAS_CONSTANT * temperature_k / pressure_pa)
}

/// Van der Waals equation: (P + a·n²/V²)(V - n·b) = nRT
///
/// - `moles`: n (mol)
/// - `temperature_k`: T (K)
/// - `volume_m3`: V (m³)
/// - `a`: attraction constant (Pa·m⁶/mol²); 1 L²·bar/mol² = 0.1 Pa·m⁶/mol²
/// - `b`: excluded volume (m³/mol); 1 L/mol = 10⁻³ m³/mol
///
/// Returns pressure in Pascals.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, or [`KimiyaError::InvalidInput`] if moles is not finite, `a` or
/// `b` is negative or non-finite, volume is not positive and finite, or volume
/// is too small (V ≤ n·b).
#[inline]
pub fn van_der_waals_pressure(
    moles: f64,
    temperature_k: f64,
    volume_m3: f64,
    a: f64,
    b: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_finite(moles) {
        return Err(KimiyaError::InvalidInput("moles must be finite".into()));
    }
    if !is_non_negative(a) || !is_non_negative(b) {
        return Err(KimiyaError::InvalidInput(
            "van der Waals constants a and b must be non-negative".into(),
        ));
    }
    if !is_positive(volume_m3) {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    if volume_m3 <= moles * b {
        return Err(KimiyaError::InvalidInput(
            "volume must be greater than n·b".into(),
        ));
    }
    let n_over_v = moles / volume_m3;
    Ok((moles * GAS_CONSTANT * temperature_k / (volume_m3 - moles * b)) - a * n_over_v * n_over_v)
}

/// Dalton's law: partial pressure = total pressure × mole fraction.
///
/// Returns the partial pressure in the unit of `total_pressure`.
#[must_use]
#[inline]
pub fn partial_pressure(total_pressure: f64, mole_fraction: f64) -> f64 {
    total_pressure * mole_fraction
}

/// Gas density: ρ = PM / (RT)
///
/// P in Pa, M in kg/mol, T in K. Returns density in kg/m³.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, or [`KimiyaError::InvalidInput`] if pressure is negative or
/// non-finite, or molar mass is not positive and finite.
#[inline]
pub fn gas_density(pressure_pa: f64, molar_mass_kg: f64, temperature_k: f64) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_non_negative(pressure_pa) {
        return Err(KimiyaError::InvalidInput(
            "pressure must be non-negative".into(),
        ));
    }
    if !is_positive(molar_mass_kg) {
        return Err(KimiyaError::InvalidInput(
            "molar mass must be positive".into(),
        ));
    }
    Ok(pressure_pa * molar_mass_kg / (GAS_CONSTANT * temperature_k))
}

/// Standard temperature for [`STP_PRESSURE`]: 273.15 K (0 °C).
pub const STP_TEMPERATURE: f64 = 273.15;

/// "Standard pressure" of the traditional STP convention: 1 atm = 101 325 Pa.
///
/// This is the pre-1982 convention (with [`STP_TEMPERATURE`] it gives the
/// familiar 22.414 L/mol ideal-gas molar volume). IUPAC has defined STP as
/// 273.15 K and 100 kPa (1 bar) since 1982, which gives 22.711 L/mol; the
/// 1 atm value is kept here for compatibility.
pub const STP_PRESSURE: f64 = 101_325.0;

// ── Peng-Robinson equation of state ──────────────────────────────────

/// Peng-Robinson EOS parameters for a substance.
///
/// The α(T) function uses the original 1976 κ polynomial,
/// κ = 0.37464 + 1.54226ω − 0.26992ω², which was fitted for
/// ω ≤ 0.49. Heavier components (ω > 0.49) are accepted but are less
/// accurate; the 1978 cubic κ(ω) correlation is not implemented.
#[derive(Debug, Clone)]
pub struct PengRobinsonParams {
    /// Critical temperature Tc (K). Must be positive and finite.
    pub tc: f64,
    /// Critical pressure Pc (Pa). Must be positive and finite.
    pub pc: f64,
    /// Acentric factor ω (dimensionless). Must be finite; the κ polynomial
    /// is valid for ω ≤ 0.49.
    pub omega: f64,
}

/// Peng-Robinson pressure: P = RT/(V-b) - a·α(T) / (V² + 2bV - b²)
///
/// - `params`: critical properties and acentric factor
/// - `temperature_k`: T (K)
/// - `molar_volume`: V (m³/mol)
///
/// Returns pressure in Pa.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature or the critical
/// temperature is not positive and finite, or [`KimiyaError::InvalidInput`]
/// if the critical pressure is not positive and finite, ω is not finite, or
/// the molar volume is non-finite or not greater than b.
#[inline]
pub fn peng_robinson_pressure(
    params: &PengRobinsonParams,
    temperature_k: f64,
    molar_volume: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    validate_pr_params(params)?;
    let (a_alpha, b) = pr_coefficients(params, temperature_k);
    if !is_finite(molar_volume) || molar_volume <= b {
        return Err(KimiyaError::InvalidInput(
            "molar volume must be greater than b".into(),
        ));
    }
    let term1 = GAS_CONSTANT * temperature_k / (molar_volume - b);
    let term2 = a_alpha / (molar_volume * molar_volume + 2.0 * b * molar_volume - b * b);
    Ok(term1 - term2)
}

/// Compressibility factor Z from Peng-Robinson.
///
/// Solves the cubic Z³ − (1−B)Z² + (A−3B²−2B)Z − (AB−B²−B³) = 0 analytically
/// (see [`peng_robinson_z_roots`]) and returns the largest physical root
/// (Z > B), i.e. the vapor-like root.
///
/// The vapor root is not necessarily the stable phase: below Tc the cubic
/// can have a liquid root too, and the stable phase is the one with the lower
/// fugacity. For example water at 300 K and 1 atm has a vapor root
/// Z ≈ 0.985 while the liquid root Z ≈ 0.00086 is the stable one. Use
/// [`peng_robinson_z_roots`] to get every root.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature or the critical
/// temperature is not positive and finite, [`KimiyaError::InvalidInput`] if
/// pressure or the critical pressure is not positive and finite or ω is not
/// finite, or [`KimiyaError::ComputationError`] if no root with Z > B is found
/// (only possible when intermediate values overflow).
pub fn peng_robinson_z(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<f64> {
    let (cap_a, cap_b) = pr_dimensionless(params, temperature_k, pressure_pa)?;
    let (roots, n) = pr_physical_roots(cap_a, cap_b);
    match n.checked_sub(1).and_then(|last| roots.get(last)) {
        Some(&z) => Ok(z),
        None => Err(no_physical_root()),
    }
}

/// All physical Peng-Robinson compressibility-factor roots, ascending.
///
/// The cubic Z³ − (1−B)Z² + (A−3B²−2B)Z − (AB−B²−B³) = 0, with
/// A = a·α·P/(RT)² and B = bP/(RT), is solved analytically (Cardano, or the
/// trigonometric form when there are three real roots). Only roots with
/// Z > B (molar volume greater than the co-volume b) are returned, so the
/// result has one to three entries; a repeated root (at a spinodal or the
/// critical point) is listed once per multiplicity.
///
/// With three roots, the smallest is liquid-like, the largest vapor-like,
/// and the middle one is mechanically unstable (∂P/∂V > 0) with no physical
/// meaning. Which of liquid and vapor is *stable* is decided by comparing
/// fugacities (lower wins), not by root size: for water at 300 K and 1 atm
/// the roots are ≈ \[0.00086, 0.0138, 0.985\] and the liquid root is stable.
///
/// # Errors
///
/// Same conditions as [`peng_robinson_z`].
pub fn peng_robinson_z_roots(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<Vec<f64>> {
    let (cap_a, cap_b) = pr_dimensionless(params, temperature_k, pressure_pa)?;
    let (roots, n) = pr_physical_roots(cap_a, cap_b);
    match roots.get(..n) {
        Some(physical) if !physical.is_empty() => Ok(physical.to_vec()),
        _ => Err(no_physical_root()),
    }
}

#[cold]
fn no_physical_root() -> KimiyaError {
    KimiyaError::ComputationError("PR Z-factor: no root with Z > B".into())
}

/// Reject non-physical critical constants (NaN/±∞ included).
#[inline]
fn validate_pr_params(params: &PengRobinsonParams) -> Result<()> {
    if !is_positive(params.tc) {
        return Err(KimiyaError::InvalidTemperature(
            "critical temperature must be positive".into(),
        ));
    }
    if !is_positive(params.pc) {
        return Err(KimiyaError::InvalidInput(
            "critical pressure must be positive".into(),
        ));
    }
    if !is_finite(params.omega) {
        return Err(KimiyaError::InvalidInput(
            "acentric factor must be finite".into(),
        ));
    }
    Ok(())
}

/// Validate inputs and return the dimensionless PR parameters (A, B).
#[inline]
fn pr_dimensionless(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<(f64, f64)> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_positive(pressure_pa) {
        return Err(KimiyaError::InvalidInput(
            "pressure must be positive".into(),
        ));
    }
    validate_pr_params(params)?;
    let (a_alpha, b) = pr_coefficients(params, temperature_k);
    let rt = GAS_CONSTANT * temperature_k;
    Ok((a_alpha * pressure_pa / (rt * rt), b * pressure_pa / rt))
}

/// Compute Peng-Robinson a·α(T) (Pa·m⁶/mol²) and b (m³/mol) from critical
/// properties. Callers must have run [`validate_pr_params`].
#[inline]
fn pr_coefficients(params: &PengRobinsonParams, temperature_k: f64) -> (f64, f64) {
    let tc = params.tc;
    let pc = params.pc;
    let omega = params.omega;

    let a = 0.45724 * GAS_CONSTANT * GAS_CONSTANT * tc * tc / pc;
    let b = 0.07780 * GAS_CONSTANT * tc / pc;
    // Original 1976 κ(ω) polynomial, fitted for ω ≤ 0.49.
    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let tr = temperature_k / tc;
    let alpha = (1.0 + kappa * (1.0 - tr.sqrt())).powi(2);
    (a * alpha, b)
}

/// Roots of the PR Z-cubic with Z > B (and finite), ascending, packed at the
/// front of the array; the second element is how many there are.
#[inline]
fn pr_physical_roots(cap_a: f64, cap_b: f64) -> ([f64; 3], usize) {
    // Z³ + c2·Z² + c1·Z + c0 = 0
    let c2 = cap_b - 1.0;
    let c1 = cap_a - 3.0 * cap_b * cap_b - 2.0 * cap_b;
    let c0 = -(cap_a * cap_b - cap_b * cap_b - cap_b * cap_b * cap_b);
    let (all, n_real) = cubic_real_roots(c2, c1, c0);
    let mut out = [0.0; 3];
    let mut n = 0;
    for &z in all.iter().take(n_real) {
        if is_finite(z)
            && z > cap_b
            && let Some(slot) = out.get_mut(n)
        {
            *slot = z;
            n += 1;
        }
    }
    (out, n)
}

/// Relative tolerance on the cubic discriminant below which a computed
/// "one real root" case is treated as a (near-)double root, so a double root
/// lost to rounding is still reported.
const CUBIC_DEGENERACY_TOL: f64 = 1e-12;

/// Real roots of z³ + c2·z² + c1·z + c0 = 0, ascending, packed at the front of
/// the array; the second element is how many (1 or 3) there are.
///
/// Uses the trigonometric form when there are three real roots and Cardano's
/// formula (in the cancellation-free form of Numerical Recipes §5.6)
/// otherwise, then polishes each root with a guarded Newton step.
#[inline]
fn cubic_real_roots(c2: f64, c1: f64, c0: f64) -> ([f64; 3], usize) {
    let shift = c2 / 3.0;
    let q = (c2 * c2 - 3.0 * c1) / 9.0;
    let r = (c2 * (2.0 * c2 * c2 - 9.0 * c1) + 27.0 * c0) / 54.0;
    let q3 = q * q * q;
    let r2 = r * r;

    let (mut roots, n) = if q > 0.0 && r2 <= q3 * (1.0 + CUBIC_DEGENERACY_TOL) {
        // Three real roots (possibly coincident). Clamp: rounding can push
        // |r|/√q³ slightly above 1 near a double root.
        let sqrt_q = q.sqrt();
        let theta = (r / (sqrt_q * q)).clamp(-1.0, 1.0).acos();
        let m = -2.0 * sqrt_q;
        let tau = std::f64::consts::TAU;
        (
            [
                m * (theta / 3.0).cos() - shift,
                m * ((theta - tau) / 3.0).cos() - shift,
                m * ((theta + tau) / 3.0).cos() - shift,
            ],
            3,
        )
    } else {
        // One real root.
        let big_a = -r.signum() * (r.abs() + (r2 - q3).max(0.0).sqrt()).cbrt();
        let big_b = if big_a == 0.0 { 0.0 } else { q / big_a };
        ([big_a + big_b - shift, 0.0, 0.0], 1)
    };

    let poly = |z: f64| ((z + c2) * z + c1) * z + c0;
    for z in roots.iter_mut().take(n) {
        for _ in 0..2 {
            let f = poly(*z);
            let df = (3.0 * *z + 2.0 * c2) * *z + c1;
            if f == 0.0 || df == 0.0 {
                break;
            }
            let next = *z - f / df;
            if poly(next).abs() < f.abs() {
                *z = next;
            } else {
                break;
            }
        }
    }
    if n == 3 {
        // Already ascending by construction; re-sort in case polishing
        // swapped two nearly coincident roots.
        let [a, b, c] = &mut roots;
        if *a > *b {
            std::mem::swap(a, b);
        }
        if *b > *c {
            std::mem::swap(b, c);
        }
        if *a > *b {
            std::mem::swap(a, b);
        }
    }
    (roots, n)
}

/// Water (H₂O): Tc = 647.1 K, Pc = 22.064 MPa, ω = 0.3443.
pub const PR_WATER: PengRobinsonParams = PengRobinsonParams {
    tc: 647.1,
    pc: 22_064_000.0,
    omega: 0.3443,
};
/// Nitrogen (N₂): Tc = 126.2 K, Pc = 3.394 MPa, ω = 0.0377.
pub const PR_NITROGEN: PengRobinsonParams = PengRobinsonParams {
    tc: 126.2,
    pc: 3_394_000.0,
    omega: 0.0377,
};
/// Oxygen (O₂): Tc = 154.6 K, Pc = 5.043 MPa, ω = 0.0222.
pub const PR_OXYGEN: PengRobinsonParams = PengRobinsonParams {
    tc: 154.6,
    pc: 5_043_000.0,
    omega: 0.0222,
};
/// Carbon dioxide (CO₂): Tc = 304.1 K, Pc = 7.375 MPa, ω = 0.2236.
pub const PR_CO2: PengRobinsonParams = PengRobinsonParams {
    tc: 304.1,
    pc: 7_375_000.0,
    omega: 0.2236,
};
/// Methane (CH₄): Tc = 190.6 K, Pc = 4.599 MPa, ω = 0.0115.
pub const PR_METHANE: PengRobinsonParams = PengRobinsonParams {
    tc: 190.6,
    pc: 4_599_000.0,
    omega: 0.0115,
};

/// Fugacity coefficient from Peng-Robinson compressibility factor:
/// ln(φ) = (Z-1) - ln(Z-B) - A/(2√2·B) × ln((Z + (1+√2)B) / (Z + (1-√2)B))
///
/// Evaluated at the vapor (largest) root returned by [`peng_robinson_z`];
/// see that function for why this is not always the stable phase.
/// Dimensionless.
///
/// # Errors
///
/// Same input conditions as [`peng_robinson_z`]; additionally returns
/// [`KimiyaError::ComputationError`] if φ is not finite and positive.
pub fn fugacity_coefficient_pr(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<f64> {
    let (cap_a, cap_b) = pr_dimensionless(params, temperature_k, pressure_pa)?;
    let (roots, n) = pr_physical_roots(cap_a, cap_b);
    let z = match n.checked_sub(1).and_then(|last| roots.get(last)) {
        Some(&z) => z,
        None => return Err(no_physical_root()),
    };

    let sqrt2 = std::f64::consts::SQRT_2;
    let ln_phi = (z - 1.0)
        - (z - cap_b).ln()
        - cap_a / (2.0 * sqrt2 * cap_b)
            * ((z + (1.0 + sqrt2) * cap_b) / (z + (1.0 - sqrt2) * cap_b)).ln();
    let phi = ln_phi.exp();
    if !is_positive(phi) {
        return Err(KimiyaError::ComputationError(
            "PR fugacity coefficient is not finite and positive".into(),
        ));
    }
    Ok(phi)
}

/// Fugacity: f = φ × P
///
/// Returns fugacity in Pa, using φ from [`fugacity_coefficient_pr`] (vapor
/// root).
///
/// # Errors
///
/// Returns error if fugacity coefficient calculation fails.
#[inline]
pub fn fugacity_pr(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<f64> {
    let phi = fugacity_coefficient_pr(params, temperature_k, pressure_pa)?;
    Ok(phi * pressure_pa)
}

// ── Redlich-Kwong equation of state ──────────────────────────────────

/// Redlich-Kwong pressure: P = RT/(V-b) - a / (T^½ × V(V+b))
///
/// - `tc`: critical temperature (K)
/// - `pc`: critical pressure (Pa)
/// - `temperature_k`: T (K)
/// - `molar_volume`: V (m³/mol)
///
/// Returns pressure in Pa.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature or `tc` is not
/// positive and finite, or [`KimiyaError::InvalidInput`] if `pc` is not
/// positive and finite or the molar volume is non-finite or not greater
/// than b.
pub fn redlich_kwong_pressure(
    tc: f64,
    pc: f64,
    temperature_k: f64,
    molar_volume: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_positive(tc) {
        return Err(KimiyaError::InvalidTemperature(
            "critical temperature must be positive".into(),
        ));
    }
    if !is_positive(pc) {
        return Err(KimiyaError::InvalidInput(
            "critical pressure must be positive".into(),
        ));
    }
    let a = 0.42748 * GAS_CONSTANT * GAS_CONSTANT * tc.powf(2.5) / pc;
    let b = 0.08664 * GAS_CONSTANT * tc / pc;
    if !is_finite(molar_volume) || molar_volume <= b {
        return Err(KimiyaError::InvalidInput(
            "molar volume must be greater than b".into(),
        ));
    }
    let term1 = GAS_CONSTANT * temperature_k / (molar_volume - b);
    let term2 = a / (temperature_k.sqrt() * molar_volume * (molar_volume + b));
    Ok(term1 - term2)
}

/// Joule-Thomson coefficient for an ideal gas (always zero, K/Pa).
///
/// Real gas JT requires EOS derivatives — use with Peng-Robinson for real gases.
#[must_use]
#[inline]
pub fn joule_thomson_ideal() -> f64 {
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ideal_gas_stp_volume() {
        let v = ideal_gas_volume(1.0, STP_TEMPERATURE, STP_PRESSURE).unwrap();
        let v_liters = v * 1000.0;
        // 1 mol at 273.15 K / 101 325 Pa = 22.413 97 L (CODATA molar volume).
        assert!(
            (v_liters - 22.413_969).abs() < 1e-5,
            "1 mol at STP should be 22.414 L, got {v_liters}"
        );
    }

    #[test]
    fn ideal_gas_roundtrip() {
        let p = ideal_gas_pressure(2.0, 300.0, 0.05).unwrap();
        let v = ideal_gas_volume(2.0, 300.0, p).unwrap();
        assert!((v - 0.05).abs() < 1e-10, "roundtrip should recover volume");
    }

    #[test]
    fn van_der_waals_approaches_ideal() {
        let p_ideal = ideal_gas_pressure(1.0, 300.0, 0.025).unwrap();
        let p_vdw = van_der_waals_pressure(1.0, 300.0, 0.025, 0.0, 0.0).unwrap();
        assert!(
            (p_ideal - p_vdw).abs() < 1e-9,
            "a=0 b=0 should match ideal gas"
        );
    }

    #[test]
    fn van_der_waals_co2_value() {
        // CO₂: a = 0.3640 Pa·m⁶/mol², b = 4.267e-5 m³/mol, 1 mol in 1 L at 300 K.
        // P = RT/(V−b) − a/V² = 2494.34/9.5733e-4 − 0.364/1e-6 = 2.2416 MPa.
        let p = van_der_waals_pressure(1.0, 300.0, 1e-3, 0.3640, 4.267e-5).unwrap();
        let expected = GAS_CONSTANT * 300.0 / (1e-3 - 4.267e-5) - 0.3640 / (1e-3 * 1e-3);
        assert!((p - expected).abs() < 1e-6);
        assert!((p - 2.2416e6).abs() < 1e3, "got {p}");
    }

    #[test]
    fn van_der_waals_rejects_bad_temperature() {
        for t in [0.0, -10.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    van_der_waals_pressure(1.0, t, 0.025, 0.1, 1e-5),
                    Err(KimiyaError::InvalidTemperature(_))
                ),
                "T = {t} should be rejected"
            );
        }
    }

    #[test]
    fn van_der_waals_rejects_nan_inputs() {
        assert!(van_der_waals_pressure(f64::NAN, 300.0, 0.025, 0.1, 1e-5).is_err());
        assert!(van_der_waals_pressure(1.0, 300.0, f64::NAN, 0.1, 1e-5).is_err());
        assert!(van_der_waals_pressure(1.0, 300.0, 0.025, f64::NAN, 1e-5).is_err());
        assert!(van_der_waals_pressure(1.0, 300.0, 0.025, 0.1, f64::NAN).is_err());
        assert!(van_der_waals_pressure(1.0, 300.0, 0.025, -0.1, 1e-5).is_err());
    }

    #[test]
    fn partial_pressure_basic() {
        let pp = partial_pressure(101325.0, 0.21);
        assert!((pp - 21278.25).abs() < 1e-9);
    }

    #[test]
    fn gas_density_air_at_stp() {
        let rho = gas_density(STP_PRESSURE, 0.029, STP_TEMPERATURE).unwrap();
        // 101325 × 0.029 / (8.314462618 × 273.15) = 1.293 86 kg/m³
        assert!(
            (rho - 1.293_86).abs() < 1e-4,
            "air density at STP should be ~1.294, got {rho}"
        );
    }

    #[test]
    fn gas_density_rejects_nan() {
        assert!(gas_density(f64::NAN, 0.029, 300.0).is_err());
        assert!(gas_density(101325.0, f64::NAN, 300.0).is_err());
        assert!(gas_density(101325.0, 0.029, f64::NAN).is_err());
        assert!(gas_density(-1.0, 0.029, 300.0).is_err());
    }

    #[test]
    fn stp_pressure_is_one_atm() {
        assert!((STP_PRESSURE - 101_325.0).abs() < f64::EPSILON);
        assert!((STP_TEMPERATURE - 273.15).abs() < f64::EPSILON);
    }

    #[test]
    fn zero_volume_is_error() {
        assert!(ideal_gas_pressure(1.0, 300.0, 0.0).is_err());
    }

    #[test]
    fn zero_pressure_is_error() {
        assert!(ideal_gas_volume(1.0, 300.0, 0.0).is_err());
    }

    #[test]
    fn zero_temp_density_is_error() {
        assert!(gas_density(101325.0, 0.029, 0.0).is_err());
    }

    #[test]
    fn negative_temp_pressure_is_error() {
        assert!(ideal_gas_pressure(1.0, -100.0, 0.025).is_err());
    }

    #[test]
    fn negative_temp_volume_is_error() {
        assert!(ideal_gas_volume(1.0, -100.0, 101325.0).is_err());
    }

    #[test]
    fn van_der_waals_volume_too_small_is_error() {
        // V = 0.001 m³ < n*b = 1.0 * 0.01 = 0.01
        assert!(van_der_waals_pressure(1.0, 300.0, 0.001, 0.0, 0.01).is_err());
    }

    #[test]
    fn nan_inputs_are_rejected() {
        assert!(ideal_gas_pressure(f64::NAN, 300.0, 0.025).is_err());
        assert!(ideal_gas_pressure(1.0, f64::NAN, 0.025).is_err());
        assert!(ideal_gas_pressure(1.0, 300.0, f64::NAN).is_err());
        assert!(ideal_gas_pressure(1.0, f64::INFINITY, 0.025).is_err());
        assert!(ideal_gas_volume(f64::NAN, 300.0, 101325.0).is_err());
        assert!(ideal_gas_volume(1.0, f64::NAN, 101325.0).is_err());
        assert!(ideal_gas_volume(1.0, 300.0, f64::NAN).is_err());
        assert!(ideal_gas_volume(1.0, 300.0, f64::INFINITY).is_err());
    }

    // ── Peng-Robinson ────────────────────────────────────────────────

    /// Residual of the PR cubic at Z, for checking roots.
    fn pr_residual(params: &PengRobinsonParams, t: f64, p: f64, z: f64) -> f64 {
        let (a, b) = pr_dimensionless(params, t, p).unwrap();
        z * z * z - (1.0 - b) * z * z + (a - 3.0 * b * b - 2.0 * b) * z
            - (a * b - b * b - b * b * b)
    }

    #[test]
    fn pr_z_ideal_at_low_pressure() {
        // At very low pressure, Z → 1 (ideal gas behavior)
        let z = peng_robinson_z(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (z - 1.0).abs() < 1e-4,
            "Z at low P should approach 1.0, got {z}"
        );
    }

    #[test]
    fn pr_z_nitrogen_at_stp() {
        // N₂ at 273.15 K / 1 atm: PR gives Z = 0.999 230 (reference: numpy.roots
        // on the same cubic).
        let z = peng_robinson_z(&PR_NITROGEN, 273.15, 101325.0).unwrap();
        assert!(
            (z - 0.999_230_19).abs() < 1e-7,
            "N₂ Z at STP should be 0.99923, got {z}"
        );
    }

    #[test]
    fn pr_z_nitrogen_1mpa() {
        // Benchmark case. Reference root: 0.996 007 28.
        let z = peng_robinson_z(&PR_NITROGEN, 300.0, 1.0e6).unwrap();
        assert!((z - 0.996_007_28).abs() < 1e-7, "got {z}");
    }

    #[test]
    fn pr_z_high_pressure_takes_physical_root() {
        // Newton from Z = 1 used to converge to the negative roots
        // (−1.65, −1.87, −1.49). Reference physical roots from numpy.roots.
        let cases = [
            (&PR_NITROGEN, 300.0, 90.0e6, 1.663_635_53),
            (&PR_NITROGEN, 300.0, 100.0e6, 1.764_087_16),
            (&PR_METHANE, 1000.0, 200.0e6, 1.604_907_83),
        ];
        for (params, t, p, expected) in cases {
            let z = peng_robinson_z(params, t, p).unwrap();
            assert!(
                (z - expected).abs() < 1e-7,
                "T={t} P={p}: expected Z={expected}, got {z}"
            );
            let roots = peng_robinson_z_roots(params, t, p).unwrap();
            assert_eq!(roots.len(), 1, "only one root has Z > B: {roots:?}");
            assert!(pr_residual(params, t, p, z).abs() < 1e-12);
        }
    }

    #[test]
    fn pr_z_consistent_with_pr_pressure() {
        // V = ZRT/P fed back into the explicit PR pressure must give P back.
        let cases = [
            (&PR_NITROGEN, 300.0, 1.0e6),
            (&PR_NITROGEN, 300.0, 90.0e6),
            (&PR_METHANE, 1000.0, 200.0e6),
            (&PR_CO2, 304.1, 7.375e6),
            (&PR_WATER, 300.0, 101_325.0),
            (&PR_OXYGEN, 200.0, 10.0e6),
        ];
        for (params, t, p) in cases {
            for z in peng_robinson_z_roots(params, t, p).unwrap() {
                let v = z * GAS_CONSTANT * t / p;
                let p_back = peng_robinson_pressure(params, t, v).unwrap();
                assert!(
                    ((p_back - p) / p).abs() < 1e-8,
                    "T={t} P={p} Z={z}: P(V) = {p_back}"
                );
            }
        }
    }

    #[test]
    fn pr_z_roots_water_three_roots() {
        // Water at 300 K / 1 atm: liquid 0.000 863 44, unstable 0.013 783 68,
        // vapor 0.984 582 22 (numpy.roots).
        let roots = peng_robinson_z_roots(&PR_WATER, 300.0, 101_325.0).unwrap();
        assert_eq!(roots.len(), 3, "{roots:?}");
        let expected = [0.000_863_440_8, 0.013_783_678, 0.984_582_222];
        for (z, e) in roots.iter().zip(expected) {
            assert!(((z - e) / e).abs() < 1e-7, "{roots:?}");
        }
        assert!(roots.windows(2).all(|w| w[0] < w[1]), "ascending");
        let z = peng_robinson_z(&PR_WATER, 300.0, 101_325.0).unwrap();
        assert!((z - expected[2]).abs() < 1e-8, "vapor root returned");
    }

    #[test]
    fn pr_z_roots_co2_subcritical_three_roots() {
        // CO₂ at 280 K / 4 MPa: 0.089 225 14, 0.202 907 08, 0.662 039 25.
        let roots = peng_robinson_z_roots(&PR_CO2, 280.0, 4.0e6).unwrap();
        let expected = [0.089_225_14, 0.202_907_08, 0.662_039_25];
        assert_eq!(roots.len(), 3, "{roots:?}");
        for (z, e) in roots.iter().zip(expected) {
            assert!((z - e).abs() < 1e-7, "{roots:?}");
        }
    }

    #[test]
    fn pr_z_near_critical_point() {
        // At (Tc, Pc) the cubic is nearly a triple root (A = 0.45724,
        // B = 0.0778); the solver must still return a valid root.
        let z = peng_robinson_z(&PR_CO2, PR_CO2.tc, PR_CO2.pc).unwrap();
        assert!((z - 0.321_379_03).abs() < 1e-6, "got {z}");
        assert!(pr_residual(&PR_CO2, PR_CO2.tc, PR_CO2.pc, z).abs() < 1e-14);
    }

    #[test]
    fn cubic_double_and_triple_roots() {
        // (z − 1)²(z − 2) = z³ − 4z² + 5z − 2
        let (r, n) = cubic_real_roots(-4.0, 5.0, -2.0);
        assert_eq!(n, 3);
        assert!((r[0] - 1.0).abs() < 1e-7 && (r[1] - 1.0).abs() < 1e-7);
        assert!((r[2] - 2.0).abs() < 1e-12);
        // (z − 0.5)³ = z³ − 1.5z² + 0.75z − 0.125
        let (r, n) = cubic_real_roots(-1.5, 0.75, -0.125);
        assert!(n >= 1);
        for z in r.iter().take(n) {
            assert!((z - 0.5).abs() < 1e-5, "{r:?}");
        }
        // z³ + z = 0 → single real root 0
        let (r, n) = cubic_real_roots(0.0, 1.0, 0.0);
        assert_eq!(n, 1);
        assert!(r[0].abs() < 1e-15);
    }

    #[test]
    fn pr_pressure_positive() {
        // V = 0.025 m³/mol at 300K for N₂
        let p = peng_robinson_pressure(&PR_NITROGEN, 300.0, 0.025).unwrap();
        assert!(p > 0.0, "pressure should be positive, got {p}");
    }

    #[test]
    fn pr_approaches_ideal_at_high_v() {
        // At very large molar volume, PR should approach ideal gas
        let p_ideal = ideal_gas_pressure(1.0, 300.0, 1.0).unwrap(); // 1 m³/mol — very dilute
        let p_pr = peng_robinson_pressure(&PR_NITROGEN, 300.0, 1.0).unwrap();
        assert!(
            (p_pr - p_ideal).abs() / p_ideal < 1e-4,
            "PR should approach ideal at large V"
        );
    }

    #[test]
    fn pr_zero_temp_is_error() {
        assert!(peng_robinson_z(&PR_NITROGEN, 0.0, 101325.0).is_err());
        assert!(peng_robinson_pressure(&PR_NITROGEN, 0.0, 0.025).is_err());
    }

    #[test]
    fn pr_rejects_nan_state() {
        assert!(peng_robinson_z(&PR_NITROGEN, f64::NAN, 101325.0).is_err());
        assert!(peng_robinson_z(&PR_NITROGEN, 300.0, f64::NAN).is_err());
        assert!(peng_robinson_z(&PR_NITROGEN, 300.0, f64::INFINITY).is_err());
        assert!(peng_robinson_z_roots(&PR_NITROGEN, 300.0, -1.0).is_err());
        assert!(peng_robinson_pressure(&PR_NITROGEN, f64::NAN, 0.025).is_err());
        assert!(peng_robinson_pressure(&PR_NITROGEN, 300.0, f64::NAN).is_err());
        assert!(peng_robinson_pressure(&PR_NITROGEN, 300.0, f64::INFINITY).is_err());
        assert!(fugacity_coefficient_pr(&PR_NITROGEN, f64::NAN, 1e5).is_err());
    }

    #[test]
    fn pr_rejects_bad_critical_constants() {
        let bad = [
            PengRobinsonParams {
                tc: 0.0,
                ..PR_NITROGEN
            },
            PengRobinsonParams {
                tc: f64::NAN,
                ..PR_NITROGEN
            },
            PengRobinsonParams {
                pc: -1.0,
                ..PR_NITROGEN
            },
            PengRobinsonParams {
                pc: f64::INFINITY,
                ..PR_NITROGEN
            },
            PengRobinsonParams {
                omega: f64::NAN,
                ..PR_NITROGEN
            },
        ];
        for params in &bad {
            assert!(peng_robinson_z(params, 300.0, 1e5).is_err(), "{params:?}");
            assert!(peng_robinson_z_roots(params, 300.0, 1e5).is_err());
            assert!(peng_robinson_pressure(params, 300.0, 0.025).is_err());
            assert!(fugacity_coefficient_pr(params, 300.0, 1e5).is_err());
        }
    }

    // ── Fugacity ─────────────────────────────────────────────────────

    #[test]
    fn fugacity_coefficient_ideal_at_low_p() {
        // At low pressure, φ → 1
        let phi = fugacity_coefficient_pr(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (phi - 1.0).abs() < 1e-4,
            "φ at low P should be ~1.0, got {phi}"
        );
    }

    #[test]
    fn fugacity_equals_pressure_at_low_p() {
        let f = fugacity_pr(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (f - 1000.0).abs() < 0.1,
            "fugacity should ≈ P at low P, got {f}"
        );
    }

    #[test]
    fn fugacity_finite_positive_at_high_pressure() {
        // Used to be Ok(NaN): ln(Z − B) of the negative root. Reference φ
        // from the physical root: 1.418 438, 1.529 092, 1.755 352.
        let cases = [
            (&PR_NITROGEN, 300.0, 90.0e6, 1.418_438_3),
            (&PR_NITROGEN, 300.0, 100.0e6, 1.529_091_7),
            (&PR_METHANE, 1000.0, 200.0e6, 1.755_352_2),
        ];
        for (params, t, p, expected) in cases {
            let phi = fugacity_coefficient_pr(params, t, p).unwrap();
            assert!(phi.is_finite() && phi > 0.0);
            assert!((phi - expected).abs() < 1e-6, "T={t} P={p}: φ={phi}");
            let f = fugacity_pr(params, t, p).unwrap();
            assert!((f - phi * p).abs() <= 1e-9 * f);
        }
    }

    #[test]
    fn fugacity_water_vapor_root() {
        // φ at the vapor root of water, 300 K / 1 atm: 0.984 806.
        let phi = fugacity_coefficient_pr(&PR_WATER, 300.0, 101_325.0).unwrap();
        assert!((phi - 0.984_806).abs() < 1e-5, "got {phi}");
    }

    #[test]
    fn joule_thomson_ideal_is_zero() {
        assert!((joule_thomson_ideal()).abs() < f64::EPSILON);
    }

    // ── Redlich-Kwong ────────────────────────────────────────────────

    #[test]
    fn rk_approaches_ideal_at_high_v() {
        // N₂: Tc=126.2, Pc=3394000
        let p_ideal = ideal_gas_pressure(1.0, 300.0, 1.0).unwrap();
        let p_rk = redlich_kwong_pressure(126.2, 3_394_000.0, 300.0, 1.0).unwrap();
        assert!(
            (p_rk - p_ideal).abs() / p_ideal < 1e-4,
            "RK should approach ideal at large V"
        );
    }

    #[test]
    fn rk_zero_temp_is_error() {
        assert!(redlich_kwong_pressure(126.2, 3_394_000.0, 0.0, 0.025).is_err());
    }

    #[test]
    fn rk_rejects_nan_inputs() {
        assert!(redlich_kwong_pressure(126.2, 3_394_000.0, f64::NAN, 0.025).is_err());
        assert!(redlich_kwong_pressure(f64::NAN, 3_394_000.0, 300.0, 0.025).is_err());
        assert!(redlich_kwong_pressure(126.2, f64::NAN, 300.0, 0.025).is_err());
        assert!(redlich_kwong_pressure(126.2, 0.0, 300.0, 0.025).is_err());
        assert!(redlich_kwong_pressure(126.2, 3_394_000.0, 300.0, f64::NAN).is_err());
        assert!(redlich_kwong_pressure(126.2, 3_394_000.0, 300.0, f64::INFINITY).is_err());
    }
}
