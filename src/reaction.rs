//! Reaction thermodynamics — Gibbs free energy (ΔG = ΔH − TΔS), the
//! equilibrium constant K = exp(−ΔG/RT), and spontaneity.
//!
//! SI units throughout: J/mol, J/(mol·K), K. Fallible functions reject NaN
//! and ±∞ inputs instead of returning `Ok(NaN)`.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_positive};

/// Gibbs free energy: ΔG = ΔH - TΔS
///
/// All values in SI: ΔH in J/mol, T in Kelvin, ΔS in J/(mol·K).
///
/// Returns ΔG in J/mol.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, or [`KimiyaError::InvalidInput`] if ΔH or ΔS is not finite.
#[inline]
pub fn gibbs_free_energy(enthalpy_j: f64, temperature_k: f64, entropy_j_per_k: f64) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_finite(enthalpy_j) || !is_finite(entropy_j_per_k) {
        return Err(KimiyaError::InvalidInput(
            "enthalpy and entropy must be finite".into(),
        ));
    }
    Ok(enthalpy_j - temperature_k * entropy_j_per_k)
}

/// Equilibrium constant from Gibbs free energy: K = exp(-ΔG / (R·T))
///
/// `delta_g_j` in J/mol, `temperature_k` in K; K is dimensionless.
/// For ΔG/RT above ≈ 745 the result underflows to 0.0.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, [`KimiyaError::InvalidInput`] if ΔG is not finite, or
/// [`KimiyaError::ComputationError`] if K overflows `f64` (−ΔG/RT > ≈ 709.8,
/// e.g. ΔG = −2 MJ/mol at 298 K); work with ln K = −ΔG/RT instead.
#[inline]
pub fn equilibrium_constant(delta_g_j: f64, temperature_k: f64) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_finite(delta_g_j) {
        return Err(KimiyaError::InvalidInput(
            "Gibbs energy must be finite".into(),
        ));
    }
    let k = (-delta_g_j / (GAS_CONSTANT * temperature_k)).exp();
    if !is_finite(k) {
        return Err(KimiyaError::ComputationError(
            "equilibrium constant overflows f64 (−ΔG/RT too large)".into(),
        ));
    }
    Ok(k)
}

/// Check if a reaction is spontaneous (ΔG < 0).
///
/// `delta_g` in any energy unit; NaN is not spontaneous.
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
            (dg - (-40400.0)).abs() < 1e-9,
            "ΔG should be ~-40400, got {dg}"
        );
    }

    #[test]
    fn gibbs_zero_temp_is_error() {
        assert!(gibbs_free_energy(-100_000.0, 0.0, -200.0).is_err());
    }

    #[test]
    fn gibbs_nan_is_error() {
        assert!(matches!(
            gibbs_free_energy(-100_000.0, f64::NAN, -200.0),
            Err(KimiyaError::InvalidTemperature(_))
        ));
        assert!(gibbs_free_energy(-100_000.0, f64::INFINITY, -200.0).is_err());
        assert!(gibbs_free_energy(f64::NAN, 298.0, -200.0).is_err());
        assert!(gibbs_free_energy(-100_000.0, 298.0, f64::NAN).is_err());
    }

    #[test]
    fn equilibrium_constant_nan_is_error() {
        assert!(matches!(
            equilibrium_constant(-5000.0, f64::NAN),
            Err(KimiyaError::InvalidTemperature(_))
        ));
        assert!(equilibrium_constant(f64::NAN, 298.0).is_err());
        assert!(equilibrium_constant(f64::NEG_INFINITY, 298.0).is_err());
    }

    #[test]
    fn equilibrium_constant_overflow_is_error() {
        // −ΔG/RT = 2e6 / (8.314 × 298) ≈ 807 > ln(f64::MAX) ≈ 709.8
        assert!(matches!(
            equilibrium_constant(-2.0e6, 298.0),
            Err(KimiyaError::ComputationError(_))
        ));
        // Large positive ΔG underflows to 0, which is a valid limit.
        let k = equilibrium_constant(2.0e6, 298.0).unwrap();
        assert!((0.0..1e-300).contains(&k));
        // Just inside the range: −ΔG/RT = 700.
        let dg = -700.0 * GAS_CONSTANT * 298.0;
        let k = equilibrium_constant(dg, 298.0).unwrap();
        assert!(((k.ln() - 700.0) / 700.0).abs() < 1e-12);
    }

    #[test]
    fn equilibrium_constant_negative_dg() {
        let k = equilibrium_constant(-5000.0, 298.0).unwrap();
        // exp(5000 / (8.314462618 × 298)) = 7.5232
        let expected = (5000.0 / (GAS_CONSTANT * 298.0)).exp();
        assert!(
            (k - expected).abs() < 1e-12 && k > 7.52 && k < 7.53,
            "got {k}"
        );
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
            (k - 1.0).abs() < f64::EPSILON,
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
