//! Basic calorimetry and thermochemistry — heat transfer (q = mcΔT),
//! Hess's law, reaction enthalpy from formation enthalpies, and specific heat
//! from calorimetry.
//!
//! Calorimetry functions use grams and J/(g·°C) (a temperature *difference*
//! in °C equals one in K); enthalpies are in J or J/mol as passed in. For
//! tabulated standard-state data and Cp(T) integration see
//! [`crate::thermochem`].

use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_positive};
use hisab::EPSILON_F64;

/// Specific heat capacity of water in J/(g·°C).
pub const WATER_SPECIFIC_HEAT: f64 = 4.184;

/// Heat transfer: q = m × c × ΔT
///
/// mass in grams, specific_heat in J/(g·°C), delta_t in °C. Returns energy in Joules.
#[must_use]
#[inline]
pub fn heat_transfer(mass_g: f64, specific_heat: f64, delta_t: f64) -> f64 {
    mass_g * specific_heat * delta_t
}

/// Hess's law: total enthalpy change = sum of step enthalpies.
///
/// Returns the sum in the unit of the inputs (e.g. J/mol).
#[must_use]
#[inline]
pub fn hess_law(enthalpies_j: &[f64]) -> f64 {
    enthalpies_j.iter().sum()
}

/// Enthalpy of reaction from formation enthalpies:
/// ΔH_rxn = Σ ΔH_f(products) - Σ ΔH_f(reactants)
///
/// Each entry is ΔH_f already multiplied by its stoichiometric coefficient
/// (e.g. J/mol of reaction); returns ΔH_rxn in the same unit.
#[must_use]
#[inline]
pub fn enthalpy_from_formation(products_hf: &[f64], reactants_hf: &[f64]) -> f64 {
    let sum_products: f64 = products_hf.iter().sum();
    let sum_reactants: f64 = reactants_hf.iter().sum();
    sum_products - sum_reactants
}

/// Specific heat at constant pressure from energy and temperature change.
/// c_p = q / (m × ΔT)
///
/// `energy_j` in J, `mass_g` in g, `delta_t` in °C (or K). Returns J/(g·°C).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if energy is not finite, mass is not
/// positive and finite, or the temperature change is zero (|ΔT| < 10⁻¹²) or
/// non-finite.
#[inline]
pub fn specific_heat_from_calorimetry(energy_j: f64, mass_g: f64, delta_t: f64) -> Result<f64> {
    if !is_finite(energy_j) {
        return Err(KimiyaError::InvalidInput("energy must be finite".into()));
    }
    if !is_positive(mass_g) {
        return Err(KimiyaError::InvalidInput("mass must be positive".into()));
    }
    if !is_finite(delta_t) || delta_t.abs() < EPSILON_F64 {
        return Err(KimiyaError::InvalidInput(
            "temperature change must be non-zero".into(),
        ));
    }
    Ok(energy_j / (mass_g * delta_t))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heat_1kg_water_1c() {
        let q = heat_transfer(1000.0, WATER_SPECIFIC_HEAT, 1.0);
        assert!(
            (q - 4184.0).abs() < 1e-9,
            "1kg water +1°C = 4184 J, got {q}"
        );
    }

    #[test]
    fn heat_100g_water_10c() {
        let q = heat_transfer(100.0, WATER_SPECIFIC_HEAT, 10.0);
        assert!(
            (q - 4184.0).abs() < 1e-9,
            "100g water +10°C = 4184 J, got {q}"
        );
    }

    #[test]
    fn heat_cooling_is_negative() {
        let q = heat_transfer(100.0, WATER_SPECIFIC_HEAT, -5.0);
        assert!(q < 0.0, "cooling should give negative heat");
    }

    #[test]
    fn hess_law_basic() {
        let total = hess_law(&[-100.0, 50.0, -30.0]);
        assert!((total - (-80.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn formation_enthalpy() {
        let dh = enthalpy_from_formation(&[-400.0], &[-200.0, -100.0]);
        assert!((dh - (-100.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn specific_heat_roundtrip() {
        let c = specific_heat_from_calorimetry(4184.0, 1000.0, 1.0).unwrap();
        assert!((c - WATER_SPECIFIC_HEAT).abs() < 1e-12);
    }

    #[test]
    fn water_specific_heat_value() {
        assert!((WATER_SPECIFIC_HEAT - 4.184).abs() < 0.001);
    }

    #[test]
    fn zero_mass_calorimetry_is_error() {
        assert!(specific_heat_from_calorimetry(100.0, 0.0, 1.0).is_err());
    }

    #[test]
    fn zero_delta_t_calorimetry_is_error() {
        assert!(specific_heat_from_calorimetry(100.0, 100.0, 0.0).is_err());
    }

    #[test]
    fn nan_calorimetry_is_error() {
        assert!(specific_heat_from_calorimetry(100.0, f64::NAN, 1.0).is_err());
        assert!(specific_heat_from_calorimetry(100.0, 100.0, f64::NAN).is_err());
        assert!(specific_heat_from_calorimetry(f64::NAN, 100.0, 1.0).is_err());
        assert!(specific_heat_from_calorimetry(100.0, f64::INFINITY, 1.0).is_err());
        assert!(specific_heat_from_calorimetry(100.0, 100.0, f64::INFINITY).is_err());
    }
}
