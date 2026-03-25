use crate::element::GAS_CONSTANT;

/// Gibbs free energy: ΔG = ΔH - TΔS
///
/// All values in SI: ΔH in J/mol, T in Kelvin, ΔS in J/(mol·K).
#[must_use]
#[inline]
pub fn gibbs_free_energy(enthalpy_j: f64, temperature_k: f64, entropy_j_per_k: f64) -> f64 {
    enthalpy_j - temperature_k * entropy_j_per_k
}

/// Equilibrium constant from Gibbs free energy: K = exp(-ΔG / (R·T))
#[must_use]
#[inline]
pub fn equilibrium_constant(delta_g_j: f64, temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    (-delta_g_j / (GAS_CONSTANT * temperature_k)).exp()
}

/// Check if a reaction is spontaneous (ΔG < 0).
#[must_use]
#[inline]
pub fn is_spontaneous(delta_g: f64) -> bool {
    delta_g < 0.0
}

/// Hess's law: total enthalpy = sum of individual step enthalpies.
#[must_use]
pub fn hess_law_enthalpy(step_enthalpies: &[f64]) -> f64 {
    step_enthalpies.iter().sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gibbs_basic() {
        // ΔH = -100 kJ, T = 298K, ΔS = -200 J/K
        let dg = gibbs_free_energy(-100_000.0, 298.0, -200.0);
        // ΔG = -100000 - 298*(-200) = -100000 + 59600 = -40400
        assert!((dg - (-40400.0)).abs() < 1.0, "ΔG should be ~-40400, got {dg}");
    }

    #[test]
    fn equilibrium_constant_negative_dg() {
        // Negative ΔG → K > 1 (products favored)
        let k = equilibrium_constant(-5000.0, 298.0);
        assert!(k > 1.0, "negative ΔG should give K > 1, got {k}");
    }

    #[test]
    fn equilibrium_constant_positive_dg() {
        // Positive ΔG → K < 1 (reactants favored)
        let k = equilibrium_constant(5000.0, 298.0);
        assert!(k < 1.0, "positive ΔG should give K < 1, got {k}");
    }

    #[test]
    fn equilibrium_constant_zero_dg() {
        let k = equilibrium_constant(0.0, 298.0);
        assert!((k - 1.0).abs() < 0.001, "zero ΔG should give K = 1, got {k}");
    }

    #[test]
    fn spontaneous_check() {
        assert!(is_spontaneous(-100.0));
        assert!(!is_spontaneous(100.0));
        assert!(!is_spontaneous(0.0));
    }

    #[test]
    fn hess_law() {
        let total = hess_law_enthalpy(&[-100.0, 50.0, -30.0]);
        assert!((total - (-80.0)).abs() < f64::EPSILON);
    }
}
