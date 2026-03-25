use crate::element::GAS_CONSTANT;

/// Arrhenius rate constant: k = A × exp(-Ea / (R × T))
///
/// A = pre-exponential factor (s⁻¹), Ea = activation energy (J/mol), T = temperature (K).
#[must_use]
#[inline]
pub fn arrhenius_rate(pre_exponential: f64, activation_energy_j: f64, temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    pre_exponential * (-activation_energy_j / (GAS_CONSTANT * temperature_k)).exp()
}

/// Half-life for a first-order reaction: t½ = ln(2) / k
#[must_use]
#[inline]
pub fn half_life_first_order(rate_constant: f64) -> f64 {
    if rate_constant <= 0.0 {
        return f64::INFINITY;
    }
    std::f64::consts::LN_2 / rate_constant
}

/// Concentration remaining for first-order decay: [A] = [A]₀ × exp(-k × t)
#[must_use]
#[inline]
pub fn first_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    initial * (-rate_constant * time).exp()
}

/// Concentration remaining for second-order decay: 1/[A] = 1/[A]₀ + k × t
#[must_use]
#[inline]
pub fn second_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    if initial <= 0.0 {
        return 0.0;
    }
    1.0 / (1.0 / initial + rate_constant * time)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_basic() {
        // A = 1e13, Ea = 100 kJ/mol, T = 300K
        let k = arrhenius_rate(1e13, 100_000.0, 300.0);
        assert!(k > 0.0 && k < 1e13, "rate should be positive and less than A, got {k}");
    }

    #[test]
    fn arrhenius_higher_temp_faster() {
        let k_low = arrhenius_rate(1e13, 50_000.0, 300.0);
        let k_high = arrhenius_rate(1e13, 50_000.0, 400.0);
        assert!(k_high > k_low, "higher temperature should give faster rate");
    }

    #[test]
    fn arrhenius_higher_ea_slower() {
        let k_low_ea = arrhenius_rate(1e13, 30_000.0, 300.0);
        let k_high_ea = arrhenius_rate(1e13, 100_000.0, 300.0);
        assert!(k_low_ea > k_high_ea, "higher activation energy should give slower rate");
    }

    #[test]
    fn half_life_basic() {
        let t = half_life_first_order(0.1);
        assert!((t - 6.931).abs() < 0.01, "t½ for k=0.1 should be ~6.93s, got {t}");
    }

    #[test]
    fn half_life_zero_rate() {
        assert!(half_life_first_order(0.0).is_infinite());
    }

    #[test]
    fn first_order_at_half_life() {
        let k = 0.1;
        let t_half = half_life_first_order(k);
        let remaining = first_order_concentration(1.0, k, t_half);
        assert!((remaining - 0.5).abs() < 0.001, "should be half at t½, got {remaining}");
    }

    #[test]
    fn second_order_decreases() {
        let c0 = 1.0;
        let c1 = second_order_concentration(c0, 0.1, 10.0);
        assert!(c1 < c0, "concentration should decrease");
        assert!(c1 > 0.0, "concentration should stay positive");
    }
}
