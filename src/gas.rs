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
        // NaN inputs should produce NaN outputs, not silently succeed
        let result = ideal_gas_pressure(f64::NAN, 300.0, 0.025).unwrap();
        assert!(result.is_nan());
    }
}
