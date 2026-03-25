use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};

/// Gibbs free energy: ΔG = ΔH - TΔS
///
/// All values in SI: ΔH in J/mol, T in Kelvin, ΔS in J/(mol·K).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn gibbs_free_energy(enthalpy_j: f64, temperature_k: f64, entropy_j_per_k: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    Ok(enthalpy_j - temperature_k * entropy_j_per_k)
}

/// Equilibrium constant from Gibbs free energy: K = exp(-ΔG / (R·T))
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn equilibrium_constant(delta_g_j: f64, temperature_k: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    Ok((-delta_g_j / (GAS_CONSTANT * temperature_k)).exp())
}

/// Check if a reaction is spontaneous (ΔG < 0).
#[must_use]
#[inline]
pub fn is_spontaneous(delta_g: f64) -> bool {
    delta_g < 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gibbs_basic() {
        let dg = gibbs_free_energy(-100_000.0, 298.0, -200.0).unwrap();
        assert!(
            (dg - (-40400.0)).abs() < 1.0,
            "ΔG should be ~-40400, got {dg}"
        );
    }

    #[test]
    fn gibbs_zero_temp_is_error() {
        assert!(gibbs_free_energy(-100_000.0, 0.0, -200.0).is_err());
    }

    #[test]
    fn equilibrium_constant_negative_dg() {
        let k = equilibrium_constant(-5000.0, 298.0).unwrap();
        assert!(k > 1.0, "negative ΔG should give K > 1, got {k}");
    }

    #[test]
    fn equilibrium_constant_positive_dg() {
        let k = equilibrium_constant(5000.0, 298.0).unwrap();
        assert!(k < 1.0, "positive ΔG should give K < 1, got {k}");
    }

    #[test]
    fn equilibrium_constant_zero_dg() {
        let k = equilibrium_constant(0.0, 298.0).unwrap();
        assert!(
            (k - 1.0).abs() < 0.001,
            "zero ΔG should give K = 1, got {k}"
        );
    }

    #[test]
    fn equilibrium_constant_zero_temp_is_error() {
        assert!(equilibrium_constant(0.0, 0.0).is_err());
    }

    #[test]
    fn spontaneous_check() {
        assert!(is_spontaneous(-100.0));
        assert!(!is_spontaneous(100.0));
        assert!(!is_spontaneous(0.0));
    }
}
