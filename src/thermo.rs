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
#[must_use]
pub fn hess_law(enthalpies_j: &[f64]) -> f64 {
    enthalpies_j.iter().sum()
}

/// Enthalpy of reaction from formation enthalpies:
/// ΔH_rxn = Σ ΔH_f(products) - Σ ΔH_f(reactants)
#[must_use]
pub fn enthalpy_from_formation(products_hf: &[f64], reactants_hf: &[f64]) -> f64 {
    let sum_products: f64 = products_hf.iter().sum();
    let sum_reactants: f64 = reactants_hf.iter().sum();
    sum_products - sum_reactants
}

/// Specific heat at constant pressure from energy and temperature change.
/// c_p = q / (m × ΔT)
#[must_use]
#[inline]
pub fn specific_heat_from_calorimetry(energy_j: f64, mass_g: f64, delta_t: f64) -> f64 {
    if mass_g <= 0.0 || delta_t.abs() < f64::EPSILON { return 0.0; }
    energy_j / (mass_g * delta_t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heat_1kg_water_1c() {
        // 1 kg = 1000 g of water heated by 1°C
        let q = heat_transfer(1000.0, WATER_SPECIFIC_HEAT, 1.0);
        assert!((q - 4184.0).abs() < 1.0, "1kg water +1°C = 4184 J, got {q}");
    }

    #[test]
    fn heat_100g_water_10c() {
        let q = heat_transfer(100.0, WATER_SPECIFIC_HEAT, 10.0);
        assert!((q - 4184.0).abs() < 1.0, "100g water +10°C = 4184 J, got {q}");
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
        // Simple: products = [-400], reactants = [-200, -100] → ΔH = -400 - (-300) = -100
        let dh = enthalpy_from_formation(&[-400.0], &[-200.0, -100.0]);
        assert!((dh - (-100.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn specific_heat_roundtrip() {
        let c = specific_heat_from_calorimetry(4184.0, 1000.0, 1.0);
        assert!((c - WATER_SPECIFIC_HEAT).abs() < 0.001);
    }

    #[test]
    fn water_specific_heat_value() {
        assert!((WATER_SPECIFIC_HEAT - 4.184).abs() < 0.001);
    }
}
