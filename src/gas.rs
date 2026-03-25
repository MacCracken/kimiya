use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};

/// Ideal gas law: P = nRT/V
///
/// n = moles, T = temperature (K), V = volume (m³). Returns pressure in Pascals.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if volume is not positive, or
/// [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn ideal_gas_pressure(moles: f64, temperature_k: f64, volume_m3: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if volume_m3 <= 0.0 {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    Ok(moles * GAS_CONSTANT * temperature_k / volume_m3)
}

/// Ideal gas law: V = nRT/P
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if pressure is not positive, or
/// [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn ideal_gas_volume(moles: f64, temperature_k: f64, pressure_pa: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if pressure_pa <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "pressure must be positive".into(),
        ));
    }
    Ok(moles * GAS_CONSTANT * temperature_k / pressure_pa)
}

/// Van der Waals equation: (P + a·n²/V²)(V - n·b) = nRT
///
/// Returns pressure in Pascals. a and b are Van der Waals constants.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if volume is too small (V ≤ n·b).
pub fn van_der_waals_pressure(
    moles: f64,
    temperature_k: f64,
    volume_m3: f64,
    a: f64,
    b: f64,
) -> Result<f64> {
    if volume_m3 <= moles * b {
        return Err(KimiyaError::InvalidInput(
            "volume must be greater than n·b".into(),
        ));
    }
    let n_over_v = moles / volume_m3;
    Ok((moles * GAS_CONSTANT * temperature_k / (volume_m3 - moles * b)) - a * n_over_v * n_over_v)
}

/// Dalton's law: partial pressure = total pressure × mole fraction.
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
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn gas_density(pressure_pa: f64, molar_mass_kg: f64, temperature_k: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    Ok(pressure_pa * molar_mass_kg / (GAS_CONSTANT * temperature_k))
}

/// STP (Standard Temperature and Pressure): 273.15 K, 101325 Pa.
pub const STP_TEMPERATURE: f64 = 273.15;
pub const STP_PRESSURE: f64 = 101325.0;

// ── Peng-Robinson equation of state ──────────────────────────────────

/// Peng-Robinson EOS parameters for a substance.
#[derive(Debug, Clone)]
pub struct PengRobinsonParams {
    /// Critical temperature Tc (K).
    pub tc: f64,
    /// Critical pressure Pc (Pa).
    pub pc: f64,
    /// Acentric factor ω.
    pub omega: f64,
}

/// Peng-Robinson pressure: P = RT/(V-b) - a·α(T) / (V² + 2bV - b²)
///
/// - `params`: critical properties and acentric factor
/// - `temperature_k`: T (K)
/// - `molar_volume`: V (m³/mol)
///
/// # Errors
///
/// Returns error if temperature or molar volume is invalid.
pub fn peng_robinson_pressure(
    params: &PengRobinsonParams,
    temperature_k: f64,
    molar_volume: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    let (a_alpha, b) = pr_coefficients(params, temperature_k);
    if molar_volume <= b {
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
/// Solves the cubic: Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
///
/// Returns the vapor-phase root (largest real root).
///
/// # Errors
///
/// Returns error if temperature or pressure is invalid.
pub fn peng_robinson_z(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if pressure_pa <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "pressure must be positive".into(),
        ));
    }
    let (a_alpha, b) = pr_coefficients(params, temperature_k);
    let rt = GAS_CONSTANT * temperature_k;
    let cap_a = a_alpha * pressure_pa / (rt * rt);
    let cap_b = b * pressure_pa / rt;

    // Cubic: Z³ + p₂Z² + p₁Z + p₀ = 0
    let p2 = -(1.0 - cap_b);
    let p1 = cap_a - 3.0 * cap_b * cap_b - 2.0 * cap_b;
    let p0 = -(cap_a * cap_b - cap_b * cap_b - cap_b * cap_b * cap_b);

    // Newton-Raphson starting from Z=1 (vapor-like)
    let f = |z: f64| z * z * z + p2 * z * z + p1 * z + p0;
    let df = |z: f64| 3.0 * z * z + 2.0 * p2 * z + p1;
    hisab::num::newton_raphson(f, df, 1.0, 1e-10, 100)
        .map_err(|e| KimiyaError::ComputationError(format!("PR Z-factor solver failed: {e}")))
}

/// Compute Peng-Robinson a·α(T) and b from critical properties.
fn pr_coefficients(params: &PengRobinsonParams, temperature_k: f64) -> (f64, f64) {
    let tc = params.tc;
    let pc = params.pc;
    let omega = params.omega;

    let a = 0.45724 * GAS_CONSTANT * GAS_CONSTANT * tc * tc / pc;
    let b = 0.07780 * GAS_CONSTANT * tc / pc;
    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let tr = temperature_k / tc;
    let alpha = (1.0 + kappa * (1.0 - tr.sqrt())).powi(2);
    (a * alpha, b)
}

/// Common Peng-Robinson parameters.
pub const PR_WATER: PengRobinsonParams = PengRobinsonParams {
    tc: 647.1,
    pc: 22_064_000.0,
    omega: 0.3443,
};
pub const PR_NITROGEN: PengRobinsonParams = PengRobinsonParams {
    tc: 126.2,
    pc: 3_394_000.0,
    omega: 0.0377,
};
pub const PR_OXYGEN: PengRobinsonParams = PengRobinsonParams {
    tc: 154.6,
    pc: 5_043_000.0,
    omega: 0.0222,
};
pub const PR_CO2: PengRobinsonParams = PengRobinsonParams {
    tc: 304.1,
    pc: 7_375_000.0,
    omega: 0.2236,
};
pub const PR_METHANE: PengRobinsonParams = PengRobinsonParams {
    tc: 190.6,
    pc: 4_599_000.0,
    omega: 0.0115,
};

/// Fugacity coefficient from Peng-Robinson compressibility factor:
/// ln(φ) = (Z-1) - ln(Z-B) - A/(2√2·B) × ln((Z + (1+√2)B) / (Z + (1-√2)B))
///
/// # Errors
///
/// Returns error if Z-factor calculation fails.
pub fn fugacity_coefficient_pr(
    params: &PengRobinsonParams,
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<f64> {
    let z = peng_robinson_z(params, temperature_k, pressure_pa)?;
    let (a_alpha, b) = pr_coefficients(params, temperature_k);
    let rt = GAS_CONSTANT * temperature_k;
    let cap_a = a_alpha * pressure_pa / (rt * rt);
    let cap_b = b * pressure_pa / rt;

    let sqrt2 = std::f64::consts::SQRT_2;
    let ln_phi = (z - 1.0)
        - (z - cap_b).ln()
        - cap_a / (2.0 * sqrt2 * cap_b)
            * ((z + (1.0 + sqrt2) * cap_b) / (z + (1.0 - sqrt2) * cap_b)).ln();
    Ok(ln_phi.exp())
}

/// Fugacity: f = φ × P
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

/// Joule-Thomson coefficient for an ideal gas (always zero).
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
        assert!(
            (v_liters - 22.414).abs() < 0.1,
            "1 mol at STP should be ~22.4 L, got {v_liters}"
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
            (p_ideal - p_vdw).abs() < 1.0,
            "a=0 b=0 should match ideal gas"
        );
    }

    #[test]
    fn partial_pressure_basic() {
        let pp = partial_pressure(101325.0, 0.21);
        assert!((pp - 21278.25).abs() < 1.0);
    }

    #[test]
    fn gas_density_air_at_stp() {
        let rho = gas_density(STP_PRESSURE, 0.029, STP_TEMPERATURE).unwrap();
        assert!(
            (rho - 1.29).abs() < 0.05,
            "air density at STP should be ~1.29, got {rho}"
        );
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
    fn nan_inputs_propagate() {
        let result = ideal_gas_pressure(f64::NAN, 300.0, 0.025).unwrap();
        assert!(result.is_nan());
    }

    // ── Peng-Robinson ────────────────────────────────────────────────

    #[test]
    fn pr_z_ideal_at_low_pressure() {
        // At very low pressure, Z → 1 (ideal gas behavior)
        let z = peng_robinson_z(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (z - 1.0).abs() < 0.01,
            "Z at low P should approach 1.0, got {z}"
        );
    }

    #[test]
    fn pr_z_nitrogen_at_stp() {
        // N₂ at STP: Z should be close to 1 (~0.9997)
        let z = peng_robinson_z(&PR_NITROGEN, 273.15, 101325.0).unwrap();
        assert!(
            (z - 1.0).abs() < 0.01,
            "N₂ Z at STP should be ~1.0, got {z}"
        );
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
            (p_pr - p_ideal).abs() / p_ideal < 0.01,
            "PR should approach ideal at large V"
        );
    }

    #[test]
    fn pr_zero_temp_is_error() {
        assert!(peng_robinson_z(&PR_NITROGEN, 0.0, 101325.0).is_err());
        assert!(peng_robinson_pressure(&PR_NITROGEN, 0.0, 0.025).is_err());
    }

    // ── Fugacity ─────────────────────────────────────────────────────

    #[test]
    fn fugacity_coefficient_ideal_at_low_p() {
        // At low pressure, φ → 1
        let phi = fugacity_coefficient_pr(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (phi - 1.0).abs() < 0.01,
            "φ at low P should be ~1.0, got {phi}"
        );
    }

    #[test]
    fn fugacity_equals_pressure_at_low_p() {
        let f = fugacity_pr(&PR_NITROGEN, 300.0, 1000.0).unwrap();
        assert!(
            (f - 1000.0).abs() < 20.0,
            "fugacity should ≈ P at low P, got {f}"
        );
    }

    #[test]
    fn joule_thomson_ideal_is_zero() {
        assert!((joule_thomson_ideal()).abs() < f64::EPSILON);
    }
}
