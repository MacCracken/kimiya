use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};

/// Arrhenius rate constant: k = A × exp(-Ea / (R × T))
///
/// A = pre-exponential factor (s⁻¹), Ea = activation energy (J/mol), T = temperature (K).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive.
#[inline]
pub fn arrhenius_rate(
    pre_exponential: f64,
    activation_energy_j: f64,
    temperature_k: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    Ok(pre_exponential * (-activation_energy_j / (GAS_CONSTANT * temperature_k)).exp())
}

/// Half-life for a first-order reaction: t½ = ln(2) / k
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if rate constant is not positive.
#[inline]
pub fn half_life_first_order(rate_constant: f64) -> Result<f64> {
    if rate_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    Ok(std::f64::consts::LN_2 / rate_constant)
}

/// Concentration remaining for first-order decay: \[A\] = \[A\]₀ × exp(-k × t)
#[must_use]
#[inline]
pub fn first_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    initial * (-rate_constant * time).exp()
}

/// Concentration remaining for second-order decay: 1/\[A\] = 1/\[A\]₀ + k × t
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if initial concentration is not positive.
#[inline]
pub fn second_order_concentration(initial: f64, rate_constant: f64, time: f64) -> Result<f64> {
    if initial <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be positive".into(),
        ));
    }
    Ok(1.0 / (1.0 / initial + rate_constant * time))
}

/// Zero-order concentration: \[A\] = \[A\]₀ - k × t
///
/// Returns zero if the result would be negative (reaction complete).
#[must_use]
#[inline]
pub fn zero_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    (initial - rate_constant * time).max(0.0)
}

/// Zero-order half-life: t½ = \[A\]₀ / (2k)
///
/// # Errors
///
/// Returns error if rate constant is not positive.
#[inline]
pub fn zero_order_half_life(initial: f64, rate_constant: f64) -> Result<f64> {
    if rate_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    Ok(initial / (2.0 * rate_constant))
}

/// General nth-order half-life: t½ = (2^(n-1) - 1) / ((n-1) × k × \[A\]₀^(n-1))
///
/// Valid for n > 1 (integer or fractional). For n=1, use [`half_life_first_order`].
/// For n=0, use [`zero_order_half_life`].
///
/// # Errors
///
/// Returns error if order ≤ 1 or parameters are invalid.
#[inline]
pub fn nth_order_half_life(initial: f64, rate_constant: f64, order: f64) -> Result<f64> {
    if order <= 1.0 {
        return Err(KimiyaError::InvalidInput(
            "order must be > 1 (use half_life_first_order for n=1)".into(),
        ));
    }
    if rate_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    if initial <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be positive".into(),
        ));
    }
    let n_minus_1 = order - 1.0;
    Ok((2.0_f64.powf(n_minus_1) - 1.0) / (n_minus_1 * rate_constant * initial.powf(n_minus_1)))
}

// ── Steady-state approximation ────────────────────────────────────────

/// Steady-state concentration of an intermediate in A →(k₁) B →(k₂) C.
///
/// At steady state: d\[B\]/dt = 0 → k₁\[A\] = k₂\[B\] → \[B\]_ss = k₁\[A\]/k₂
///
/// # Errors
///
/// Returns error if k₂ is not positive.
#[inline]
pub fn steady_state_intermediate(
    k1: f64,
    concentration_a: f64,
    k2: f64,
) -> Result<f64> {
    if k2 <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "k2 must be positive".into(),
        ));
    }
    Ok(k1 * concentration_a / k2)
}

/// Rate of product formation under steady-state approximation.
///
/// For A →(k₁) B →(k₂) C: rate = k₂\[B\]_ss = k₁\[A\]
///
/// The rate-determining step determines the overall rate.
/// If k₁ << k₂, rate ≈ k₁\[A\] (first step is rate-determining).
/// If k₂ << k₁, rate ≈ k₂\[B\] (second step is rate-determining).
#[must_use]
#[inline]
pub fn steady_state_rate(k_slow: f64, concentration: f64) -> f64 {
    k_slow * concentration
}

/// Pre-equilibrium approximation: fast equilibrium followed by slow step.
///
/// A ⇌(K_eq) B →(k₂) C where K_eq = k₁/k₋₁
///
/// Rate = k₂ × K_eq × \[A\] = k₂ × (k₁/k₋₁) × \[A\]
///
/// # Errors
///
/// Returns error if k_reverse is not positive.
#[inline]
pub fn pre_equilibrium_rate(
    k_forward: f64,
    k_reverse: f64,
    k_slow: f64,
    concentration_a: f64,
) -> Result<f64> {
    if k_reverse <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "reverse rate constant must be positive".into(),
        ));
    }
    Ok(k_slow * (k_forward / k_reverse) * concentration_a)
}

// ── Michaelis-Menten enzyme kinetics ─────────────────────────────────

/// Michaelis-Menten rate: v = Vmax × \[S\] / (Km + \[S\])
///
/// - `v_max`: maximum reaction rate
/// - `km`: Michaelis constant (substrate concentration at half Vmax)
/// - `substrate`: substrate concentration \[S\]
///
/// # Errors
///
/// Returns error if Vmax, Km, or substrate is not positive.
#[inline]
pub fn michaelis_menten(v_max: f64, km: f64, substrate: f64) -> Result<f64> {
    if v_max <= 0.0 {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if km <= 0.0 {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    if substrate < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "substrate concentration must be non-negative".into(),
        ));
    }
    Ok(v_max * substrate / (km + substrate))
}

/// Lineweaver-Burk linearization: 1/v = (Km/Vmax)(1/\[S\]) + 1/Vmax
///
/// Returns (slope, intercept) for a Lineweaver-Burk plot.
///
/// - slope = Km / Vmax
/// - intercept = 1 / Vmax
///
/// # Errors
///
/// Returns error if Vmax or Km is not positive.
#[inline]
pub fn lineweaver_burk(v_max: f64, km: f64) -> Result<(f64, f64)> {
    if v_max <= 0.0 {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if km <= 0.0 {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    Ok((km / v_max, 1.0 / v_max))
}

// ── Collision theory ─────────────────────────────────────────────────

/// Collision theory rate constant:
/// k = Z × p × exp(-Ea / (RT))
///
/// - `collision_frequency`: Z (collisions per second per unit concentration²)
/// - `steric_factor`: p (probability of correct orientation, 0 < p ≤ 1)
/// - `activation_energy_j`: Ea in J/mol
/// - `temperature_k`: T in Kelvin
///
/// # Errors
///
/// Returns error if temperature or steric factor is invalid.
#[inline]
pub fn collision_theory_rate(
    collision_frequency: f64,
    steric_factor: f64,
    activation_energy_j: f64,
    temperature_k: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if steric_factor <= 0.0 || steric_factor > 1.0 {
        return Err(KimiyaError::InvalidInput(
            "steric factor must be in (0, 1]".into(),
        ));
    }
    Ok(collision_frequency
        * steric_factor
        * (-activation_energy_j / (GAS_CONSTANT * temperature_k)).exp())
}

// ── Transition state theory (Eyring equation) ────────────────────────

/// Boltzmann constant (J/K).
pub const BOLTZMANN: f64 = 1.380649e-23;

/// Eyring equation (transition state theory):
/// k = (kB × T / h) × exp(-ΔG‡ / (RT))
///
/// - `delta_g_activation_j`: ΔG‡ (Gibbs energy of activation) in J/mol
/// - `temperature_k`: T in Kelvin
///
/// Returns rate constant in s⁻¹.
///
/// # Errors
///
/// Returns error if temperature is not positive.
#[inline]
pub fn eyring_rate(delta_g_activation_j: f64, temperature_k: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    let prefactor = BOLTZMANN * temperature_k / crate::spectroscopy::PLANCK;
    Ok(prefactor * (-delta_g_activation_j / (GAS_CONSTANT * temperature_k)).exp())
}

/// Eyring equation decomposed into enthalpy and entropy of activation:
/// k = (kB × T / h) × exp(ΔS‡/R) × exp(-ΔH‡/(RT))
///
/// - `delta_h_activation_j`: ΔH‡ in J/mol
/// - `delta_s_activation_j_per_k`: ΔS‡ in J/(mol·K)
/// - `temperature_k`: T in Kelvin
///
/// # Errors
///
/// Returns error if temperature is not positive.
#[inline]
pub fn eyring_rate_from_activation(
    delta_h_activation_j: f64,
    delta_s_activation_j_per_k: f64,
    temperature_k: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    let prefactor = BOLTZMANN * temperature_k / crate::spectroscopy::PLANCK;
    Ok(prefactor
        * (delta_s_activation_j_per_k / GAS_CONSTANT).exp()
        * (-delta_h_activation_j / (GAS_CONSTANT * temperature_k)).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_basic() {
        let k = arrhenius_rate(1e13, 100_000.0, 300.0).unwrap();
        assert!(
            k > 0.0 && k < 1e13,
            "rate should be positive and less than A, got {k}"
        );
    }

    #[test]
    fn arrhenius_higher_temp_faster() {
        let k_low = arrhenius_rate(1e13, 50_000.0, 300.0).unwrap();
        let k_high = arrhenius_rate(1e13, 50_000.0, 400.0).unwrap();
        assert!(k_high > k_low, "higher temperature should give faster rate");
    }

    #[test]
    fn arrhenius_higher_ea_slower() {
        let k_low_ea = arrhenius_rate(1e13, 30_000.0, 300.0).unwrap();
        let k_high_ea = arrhenius_rate(1e13, 100_000.0, 300.0).unwrap();
        assert!(
            k_low_ea > k_high_ea,
            "higher activation energy should give slower rate"
        );
    }

    #[test]
    fn arrhenius_zero_temp_is_error() {
        assert!(arrhenius_rate(1e13, 50_000.0, 0.0).is_err());
    }

    #[test]
    fn half_life_basic() {
        let t = half_life_first_order(0.1).unwrap();
        assert!(
            (t - 6.931).abs() < 0.01,
            "t½ for k=0.1 should be ~6.93s, got {t}"
        );
    }

    #[test]
    fn half_life_zero_rate_is_error() {
        assert!(half_life_first_order(0.0).is_err());
    }

    #[test]
    fn first_order_at_half_life() {
        let k = 0.1;
        let t_half = half_life_first_order(k).unwrap();
        let remaining = first_order_concentration(1.0, k, t_half);
        assert!(
            (remaining - 0.5).abs() < 0.001,
            "should be half at t½, got {remaining}"
        );
    }

    #[test]
    fn second_order_decreases() {
        let c0 = 1.0;
        let c1 = second_order_concentration(c0, 0.1, 10.0).unwrap();
        assert!(c1 < c0, "concentration should decrease");
        assert!(c1 > 0.0, "concentration should stay positive");
    }

    #[test]
    fn second_order_zero_initial_is_error() {
        assert!(second_order_concentration(0.0, 0.1, 10.0).is_err());
    }

    #[test]
    fn arrhenius_negative_temp_is_error() {
        assert!(arrhenius_rate(1e13, 50_000.0, -100.0).is_err());
    }

    #[test]
    fn half_life_negative_rate_is_error() {
        assert!(half_life_first_order(-0.1).is_err());
    }

    #[test]
    fn first_order_at_time_zero() {
        let remaining = first_order_concentration(1.0, 0.1, 0.0);
        assert!((remaining - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn arrhenius_zero_ea_equals_pre_exponential() {
        // With Ea=0, k = A * exp(0) = A
        let k = arrhenius_rate(1e13, 0.0, 300.0).unwrap();
        assert!((k - 1e13).abs() < 1.0, "zero Ea should give k=A, got {k}");
    }

    // ── Michaelis-Menten ─────────────────────────────────────────────

    #[test]
    fn michaelis_menten_at_km() {
        // When [S] = Km, v = Vmax/2
        let v = michaelis_menten(100.0, 5.0, 5.0).unwrap();
        assert!((v - 50.0).abs() < f64::EPSILON);
    }

    #[test]
    fn michaelis_menten_saturation() {
        // At very high [S], v → Vmax
        let v = michaelis_menten(100.0, 5.0, 1e6).unwrap();
        assert!((v - 100.0).abs() < 0.01);
    }

    #[test]
    fn michaelis_menten_zero_substrate() {
        let v = michaelis_menten(100.0, 5.0, 0.0).unwrap();
        assert!(v.abs() < f64::EPSILON);
    }

    #[test]
    fn michaelis_menten_zero_vmax_is_error() {
        assert!(michaelis_menten(0.0, 5.0, 1.0).is_err());
    }

    #[test]
    fn michaelis_menten_zero_km_is_error() {
        assert!(michaelis_menten(100.0, 0.0, 1.0).is_err());
    }

    #[test]
    fn lineweaver_burk_basic() {
        let (slope, intercept) = lineweaver_burk(100.0, 5.0).unwrap();
        assert!((slope - 0.05).abs() < 1e-10); // Km/Vmax = 5/100
        assert!((intercept - 0.01).abs() < 1e-10); // 1/Vmax = 1/100
    }

    // ── Collision theory ─────────────────────────────────────────────

    #[test]
    fn collision_theory_basic() {
        let k = collision_theory_rate(1e10, 0.1, 50_000.0, 300.0).unwrap();
        assert!(k > 0.0);
        assert!(k < 1e10, "rate should be less than Z×p");
    }

    #[test]
    fn collision_theory_steric_reduces_rate() {
        let k_full = collision_theory_rate(1e10, 1.0, 50_000.0, 300.0).unwrap();
        let k_steric = collision_theory_rate(1e10, 0.01, 50_000.0, 300.0).unwrap();
        assert!(k_steric < k_full);
    }

    #[test]
    fn collision_theory_zero_steric_is_error() {
        assert!(collision_theory_rate(1e10, 0.0, 50_000.0, 300.0).is_err());
    }

    #[test]
    fn collision_theory_steric_over_one_is_error() {
        assert!(collision_theory_rate(1e10, 1.5, 50_000.0, 300.0).is_err());
    }

    // ── Eyring equation ──────────────────────────────────────────────

    #[test]
    fn eyring_basic() {
        // ΔG‡ = 80 kJ/mol at 298K should give a reasonable rate
        let k = eyring_rate(80_000.0, 298.0).unwrap();
        assert!(k > 0.0 && k < 1e15);
    }

    #[test]
    fn eyring_higher_temp_faster() {
        let k_low = eyring_rate(80_000.0, 298.0).unwrap();
        let k_high = eyring_rate(80_000.0, 400.0).unwrap();
        assert!(k_high > k_low);
    }

    #[test]
    fn eyring_zero_temp_is_error() {
        assert!(eyring_rate(80_000.0, 0.0).is_err());
    }

    #[test]
    fn eyring_decomposed_matches_combined() {
        let dh = 80_000.0;
        let ds = -50.0;
        let t = 298.0;
        let dg = dh - t * ds;
        let k_combined = eyring_rate(dg, t).unwrap();
        let k_decomposed = eyring_rate_from_activation(dh, ds, t).unwrap();
        assert!(
            (k_combined - k_decomposed).abs() / k_combined < 1e-6,
            "combined ({k_combined}) should match decomposed ({k_decomposed})"
        );
    }

    // ── Zero-order + nth-order ───────────────────────────────────────

    #[test]
    fn zero_order_linear_decrease() {
        let c = zero_order_concentration(1.0, 0.1, 5.0);
        assert!((c - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn zero_order_floors_at_zero() {
        let c = zero_order_concentration(1.0, 0.1, 20.0);
        assert!((c).abs() < f64::EPSILON, "should floor at 0, got {c}");
    }

    #[test]
    fn zero_order_half_life_basic() {
        let t = zero_order_half_life(1.0, 0.1).unwrap();
        assert!((t - 5.0).abs() < f64::EPSILON);
    }

    #[test]
    fn nth_order_second_order_matches() {
        // For n=2: t½ = 1/(k·[A]₀)
        let t_nth = nth_order_half_life(1.0, 0.1, 2.0).unwrap();
        let t_direct = 1.0 / (0.1 * 1.0);
        assert!(
            (t_nth - t_direct).abs() < 0.001,
            "n=2 half-life should be {t_direct}, got {t_nth}"
        );
    }

    #[test]
    fn nth_order_increases_with_order() {
        let t2 = nth_order_half_life(1.0, 0.1, 2.0).unwrap();
        let t3 = nth_order_half_life(1.0, 0.1, 3.0).unwrap();
        assert!(
            t3 > t2,
            "higher order should have longer half-life at same [A]₀"
        );
    }

    #[test]
    fn nth_order_below_1_is_error() {
        assert!(nth_order_half_life(1.0, 0.1, 1.0).is_err());
        assert!(nth_order_half_life(1.0, 0.1, 0.5).is_err());
    }
}
