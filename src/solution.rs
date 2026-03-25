/// Molarity: M = moles / volume_liters
#[must_use]
#[inline]
pub fn molarity(moles: f64, volume_liters: f64) -> f64 {
    if volume_liters <= 0.0 { return 0.0; }
    moles / volume_liters
}

/// Molality: m = moles / solvent_kg
#[must_use]
#[inline]
pub fn molality(moles: f64, solvent_kg: f64) -> f64 {
    if solvent_kg <= 0.0 { return 0.0; }
    moles / solvent_kg
}

/// Dilution: M1×V1 = M2×V2 → M2 = M1×V1/V2
#[must_use]
#[inline]
pub fn dilution(m1: f64, v1: f64, v2: f64) -> f64 {
    if v2 <= 0.0 { return 0.0; }
    m1 * v1 / v2
}

/// pH from hydrogen ion concentration: pH = -log₁₀([H⁺])
#[must_use]
#[inline]
pub fn ph_from_h_concentration(h_molar: f64) -> f64 {
    if h_molar <= 0.0 { return 14.0; }
    -h_molar.log10()
}

/// pOH from hydroxide ion concentration: pOH = -log₁₀([OH⁻])
#[must_use]
#[inline]
pub fn poh_from_oh_concentration(oh_molar: f64) -> f64 {
    if oh_molar <= 0.0 { return 14.0; }
    -oh_molar.log10()
}

/// pH + pOH = 14 (at 25°C)
#[must_use]
#[inline]
pub fn ph_from_poh(poh: f64) -> f64 {
    14.0 - poh
}

/// Henderson-Hasselbalch equation: pH = pKa + log₁₀([A⁻]/[HA])
#[must_use]
#[inline]
pub fn henderson_hasselbalch(pka: f64, conjugate_base: f64, acid: f64) -> f64 {
    if acid <= 0.0 || conjugate_base <= 0.0 { return pka; }
    pka + (conjugate_base / acid).log10()
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
        // [H+] = 1e-7 → pH = 7
        let ph = ph_from_h_concentration(1e-7);
        assert!((ph - 7.0).abs() < 0.01, "pure water pH should be ~7, got {ph}");
    }

    #[test]
    fn strong_acid_ph() {
        // 0.1M HCl → pH ≈ 1
        let ph = ph_from_h_concentration(0.1);
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
        let back = ph_from_h_concentration(h);
        assert!((back - ph).abs() < 0.001);
    }

    #[test]
    fn molarity_basic() {
        let m = molarity(0.5, 0.25);
        assert!((m - 2.0).abs() < f64::EPSILON, "0.5 mol in 0.25 L = 2M");
    }

    #[test]
    fn dilution_basic() {
        // 1M × 100mL diluted to 500mL → 0.2M
        let m2 = dilution(1.0, 0.1, 0.5);
        assert!((m2 - 0.2).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_equal_ratio() {
        // When [A-] = [HA], pH = pKa
        let ph = henderson_hasselbalch(4.75, 0.1, 0.1);
        assert!((ph - 4.75).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_more_base() {
        // More conjugate base → pH > pKa
        let ph = henderson_hasselbalch(4.75, 1.0, 0.1);
        assert!(ph > 4.75);
    }
}
