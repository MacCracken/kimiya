use crate::error::{KimiyaError, Result};

/// Molarity: M = moles / volume_liters
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if volume is not positive.
#[inline]
pub fn molarity(moles: f64, volume_liters: f64) -> Result<f64> {
    if volume_liters <= 0.0 {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    Ok(moles / volume_liters)
}

/// Molality: m = moles / solvent_kg
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if solvent mass is not positive.
#[inline]
pub fn molality(moles: f64, solvent_kg: f64) -> Result<f64> {
    if solvent_kg <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "solvent mass must be positive".into(),
        ));
    }
    Ok(moles / solvent_kg)
}

/// Dilution: M1×V1 = M2×V2 → M2 = M1×V1/V2
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if final volume is not positive.
#[inline]
pub fn dilution(m1: f64, v1: f64, v2: f64) -> Result<f64> {
    if v2 <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "final volume must be positive".into(),
        ));
    }
    Ok(m1 * v1 / v2)
}

/// pH from hydrogen ion concentration: pH = -log₁₀([H⁺])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not positive.
#[inline]
pub fn ph_from_h_concentration(h_molar: f64) -> Result<f64> {
    if h_molar <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "[H⁺] must be positive".into(),
        ));
    }
    Ok(-h_molar.log10())
}

/// pOH from hydroxide ion concentration: pOH = -log₁₀([OH⁻])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not positive.
#[inline]
pub fn poh_from_oh_concentration(oh_molar: f64) -> Result<f64> {
    if oh_molar <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "[OH⁻] must be positive".into(),
        ));
    }
    Ok(-oh_molar.log10())
}

/// pH + pOH = 14 (at 25°C)
#[must_use]
#[inline]
pub fn ph_from_poh(poh: f64) -> f64 {
    14.0 - poh
}

/// Henderson-Hasselbalch equation: pH = pKa + log₁₀([A⁻]/[HA])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if acid or base concentration is not positive.
#[inline]
pub fn henderson_hasselbalch(pka: f64, conjugate_base: f64, acid: f64) -> Result<f64> {
    if acid <= 0.0 || conjugate_base <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "acid and conjugate base concentrations must be positive".into(),
        ));
    }
    Ok(pka + (conjugate_base / acid).log10())
}

/// H⁺ concentration from pH: [H⁺] = 10^(-pH)
#[must_use]
#[inline]
pub fn h_concentration_from_ph(ph: f64) -> f64 {
    10.0_f64.powf(-ph)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_water_ph() {
        let ph = ph_from_h_concentration(1e-7).unwrap();
        assert!(
            (ph - 7.0).abs() < 0.01,
            "pure water pH should be ~7, got {ph}"
        );
    }

    #[test]
    fn strong_acid_ph() {
        let ph = ph_from_h_concentration(0.1).unwrap();
        assert!((ph - 1.0).abs() < 0.01);
    }

    #[test]
    fn ph_poh_sum_14() {
        let ph = 3.0;
        let poh = 14.0 - ph;
        assert!((ph_from_poh(poh) - ph).abs() < f64::EPSILON);
    }

    #[test]
    fn ph_h_roundtrip() {
        let ph = 4.5;
        let h = h_concentration_from_ph(ph);
        let back = ph_from_h_concentration(h).unwrap();
        assert!((back - ph).abs() < 0.001);
    }

    #[test]
    fn molarity_basic() {
        let m = molarity(0.5, 0.25).unwrap();
        assert!((m - 2.0).abs() < f64::EPSILON, "0.5 mol in 0.25 L = 2M");
    }

    #[test]
    fn dilution_basic() {
        let m2 = dilution(1.0, 0.1, 0.5).unwrap();
        assert!((m2 - 0.2).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_equal_ratio() {
        let ph = henderson_hasselbalch(4.75, 0.1, 0.1).unwrap();
        assert!((ph - 4.75).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_more_base() {
        let ph = henderson_hasselbalch(4.75, 1.0, 0.1).unwrap();
        assert!(ph > 4.75);
    }

    #[test]
    fn zero_volume_molarity_is_error() {
        assert!(molarity(1.0, 0.0).is_err());
    }

    #[test]
    fn zero_concentration_ph_is_error() {
        assert!(ph_from_h_concentration(0.0).is_err());
    }

    #[test]
    fn negative_concentration_ph_is_error() {
        assert!(ph_from_h_concentration(-1.0).is_err());
    }

    #[test]
    fn molality_zero_solvent_is_error() {
        assert!(molality(1.0, 0.0).is_err());
    }

    #[test]
    fn dilution_zero_volume_is_error() {
        assert!(dilution(1.0, 0.1, 0.0).is_err());
    }

    #[test]
    fn henderson_hasselbalch_zero_acid_is_error() {
        assert!(henderson_hasselbalch(4.75, 0.1, 0.0).is_err());
    }

    #[test]
    fn henderson_hasselbalch_zero_base_is_error() {
        assert!(henderson_hasselbalch(4.75, 0.0, 0.1).is_err());
    }

    #[test]
    fn poh_zero_concentration_is_error() {
        assert!(poh_from_oh_concentration(0.0).is_err());
    }

    #[test]
    fn h_concentration_from_ph_roundtrip_strong_acid() {
        // pH 0 → [H⁺] = 1.0 M
        let h = h_concentration_from_ph(0.0);
        assert!((h - 1.0).abs() < f64::EPSILON);
    }
}
