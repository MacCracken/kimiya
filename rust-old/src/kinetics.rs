//! Chemical kinetics — rate constants, integrated rate laws, half-lives,
//! steady-state and pre-equilibrium approximations, Michaelis–Menten enzyme
//! kinetics, collision theory and transition-state (Eyring) theory.
//!
//! Units: concentrations in mol/L, time in s, energies in J/mol,
//! temperatures in K. Rate constants carry the units of their reaction
//! order n, (mol/L)^(1−n)·s⁻¹: zero order mol·L⁻¹·s⁻¹, first order s⁻¹,
//! second order L·mol⁻¹·s⁻¹. Any consistent time unit may be used instead of
//! seconds, except in the Eyring functions, whose prefactor k_B·T/h is in s⁻¹.
//!
//! For numerical time integration of reaction systems see
//! [`crate::reaction_dynamics`].

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};

/// Reject a non-finite energy argument.
#[inline]
fn require_finite(value: f64, what: &str) -> Result<()> {
    if is_finite(value) {
        Ok(())
    } else {
        Err(KimiyaError::InvalidInput(format!("{what} must be finite")))
    }
}

/// Arrhenius rate constant: k = A × exp(-Ea / (R × T))
///
/// - `pre_exponential`: A, in the units of the returned rate constant
/// - `activation_energy_j`: Ea (J/mol)
/// - `temperature_k`: T (K)
///
/// Returns k in the units of A.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, and [`KimiyaError::InvalidInput`] if A or Ea is not finite.
#[inline]
pub fn arrhenius_rate(
    pre_exponential: f64,
    activation_energy_j: f64,
    temperature_k: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    require_finite(pre_exponential, "pre-exponential factor")?;
    require_finite(activation_energy_j, "activation energy")?;
    Ok(pre_exponential * (-activation_energy_j / (GAS_CONSTANT * temperature_k)).exp())
}

/// Half-life for a first-order reaction: t½ = ln(2) / k
///
/// - `rate_constant`: k (s⁻¹)
///
/// Returns t½ in s.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if rate constant is not positive and
/// finite.
#[inline]
pub fn half_life_first_order(rate_constant: f64) -> Result<f64> {
    if !is_positive(rate_constant) {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    Ok(std::f64::consts::LN_2 / rate_constant)
}

/// Concentration remaining for first-order decay: \[A\] = \[A\]₀ × exp(-k × t)
///
/// - `initial`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (s⁻¹)
/// - `time`: t (s)
///
/// Returns \[A\] in mol/L. Inputs are not validated; non-finite input
/// propagates to the result.
#[must_use]
#[inline]
pub fn first_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    initial * (-rate_constant * time).exp()
}

/// Concentration remaining for second-order decay: 1/\[A\] = 1/\[A\]₀ + k × t
///
/// - `initial`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (L·mol⁻¹·s⁻¹)
/// - `time`: t (s)
///
/// Returns \[A\] in mol/L, which lies in (0, \[A\]₀\].
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if the initial concentration
/// is not positive and finite, and [`KimiyaError::InvalidInput`] if the rate
/// constant or time is negative or not finite (k·t ≤ −1/\[A\]₀ would give an
/// infinite or negative concentration).
#[inline]
pub fn second_order_concentration(initial: f64, rate_constant: f64, time: f64) -> Result<f64> {
    if !is_positive(initial) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be positive".into(),
        ));
    }
    if !is_non_negative(rate_constant) {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be non-negative".into(),
        ));
    }
    if !is_non_negative(time) {
        return Err(KimiyaError::InvalidInput(
            "time must be non-negative".into(),
        ));
    }
    Ok(1.0 / (1.0 / initial + rate_constant * time))
}

/// Zero-order concentration: \[A\] = \[A\]₀ - k × t
///
/// - `initial`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (mol·L⁻¹·s⁻¹)
/// - `time`: t (s)
///
/// Returns \[A\] in mol/L, or zero if the result would be negative (reaction
/// complete). Inputs are not validated; NaN input propagates to the result.
#[must_use]
#[inline]
pub fn zero_order_concentration(initial: f64, rate_constant: f64, time: f64) -> f64 {
    let remaining = initial - rate_constant * time;
    // Not `f64::max`, which would turn NaN into 0.0.
    if remaining < 0.0 { 0.0 } else { remaining }
}

/// Zero-order half-life: t½ = \[A\]₀ / (2k)
///
/// - `initial`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (mol·L⁻¹·s⁻¹)
///
/// Returns t½ in s.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if the rate constant is not positive
/// and finite, and [`KimiyaError::InvalidConcentration`] if the initial
/// concentration is not positive and finite.
#[inline]
pub fn zero_order_half_life(initial: f64, rate_constant: f64) -> Result<f64> {
    if !is_positive(rate_constant) {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    if !is_positive(initial) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be positive".into(),
        ));
    }
    Ok(initial / (2.0 * rate_constant))
}

/// Orders within this distance of 1 are rejected by [`nth_order_half_life`].
const FIRST_ORDER_TOLERANCE: f64 = 1e-12;

/// General nth-order half-life: t½ = (2^(n-1) - 1) / ((n-1) × k × \[A\]₀^(n-1))
///
/// Valid for any finite order n ≠ 1, integer or fractional (n = 0 gives
/// \[A\]₀/(2k), n = 2 gives 1/(k\[A\]₀)). For n = 1, use
/// [`half_life_first_order`]. 2^(n-1) − 1 is evaluated with `exp_m1`, so
/// orders close to 1 keep full precision.
///
/// - `initial`: \[A\]₀ (mol/L)
/// - `rate_constant`: k ((mol/L)^(1−n)·s⁻¹)
/// - `order`: n
///
/// Returns t½ in s.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if the order is not finite or is
/// within 1e-12 of 1, or the rate constant is not positive and finite;
/// [`KimiyaError::InvalidConcentration`] if the initial concentration is not
/// positive and finite; and [`KimiyaError::ComputationError`] if the result
/// overflows or underflows (extreme orders).
#[inline]
pub fn nth_order_half_life(initial: f64, rate_constant: f64, order: f64) -> Result<f64> {
    if !is_finite(order) || (order - 1.0).abs() <= FIRST_ORDER_TOLERANCE {
        return Err(KimiyaError::InvalidInput(
            "order must be finite and ≠ 1 (use half_life_first_order for n=1)".into(),
        ));
    }
    if !is_positive(rate_constant) {
        return Err(KimiyaError::InvalidInput(
            "rate constant must be positive".into(),
        ));
    }
    if !is_positive(initial) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be positive".into(),
        ));
    }
    let n_minus_1 = order - 1.0;
    // 2^(n−1) − 1 = expm1((n−1)·ln 2); same sign as n − 1, so t½ > 0.
    let numerator = (n_minus_1 * std::f64::consts::LN_2).exp_m1();
    let t_half = numerator / (n_minus_1 * rate_constant * initial.powf(n_minus_1));
    if is_positive(t_half) {
        Ok(t_half)
    } else {
        Err(KimiyaError::ComputationError(
            "half-life is not representable as a finite positive f64".into(),
        ))
    }
}

// ── Steady-state approximation ────────────────────────────────────────

/// Steady-state concentration of an intermediate in A →(k₁) B →(k₂) C.
///
/// At steady state: d\[B\]/dt = 0 → k₁\[A\] = k₂\[B\] → \[B\]_ss = k₁\[A\]/k₂
///
/// The approximation is only valid when the intermediate is consumed much
/// faster than it is formed (k₂ ≫ k₁).
///
/// - `k1`, `k2`: first-order rate constants (s⁻¹)
/// - `concentration_a`: \[A\] (mol/L)
///
/// Returns \[B\]_ss in mol/L.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if k₂ is not positive and finite or
/// k₁ is negative or not finite, and [`KimiyaError::InvalidConcentration`] if
/// \[A\] is negative or not finite.
#[inline]
pub fn steady_state_intermediate(k1: f64, concentration_a: f64, k2: f64) -> Result<f64> {
    if !is_positive(k2) {
        return Err(KimiyaError::InvalidInput("k2 must be positive".into()));
    }
    if !is_non_negative(k1) {
        return Err(KimiyaError::InvalidInput("k1 must be non-negative".into()));
    }
    if !is_non_negative(concentration_a) {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be non-negative".into(),
        ));
    }
    Ok(k1 * concentration_a / k2)
}

/// Rate of the rate-limiting step: rate = k_slow × \[X\]
///
/// For A →(k₁) B →(k₂) C:
///
/// - If k₂ ≫ k₁ (intermediate consumed much faster than formed), the
///   steady-state approximation d\[B\]/dt ≈ 0 holds and the first step is
///   rate-determining: rate = k₂\[B\]_ss = k₁\[A\]. Pass `k_slow` = k₁ and
///   `concentration` = \[A\].
/// - If k₂ ≪ k₁, B accumulates, the steady-state approximation does *not*
///   apply, and the second step is rate-determining: rate = k₂\[B\] with the
///   actual (not steady-state) \[B\]. Pass `k_slow` = k₂ and
///   `concentration` = \[B\].
///
/// Units: k_slow in s⁻¹, concentration in mol/L → rate in mol·L⁻¹·s⁻¹.
/// Inputs are not validated; non-finite input propagates to the result.
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
/// - `k_forward`, `k_reverse`: k₁, k₋₁ of the fast equilibrium (same units)
/// - `k_slow`: k₂ (s⁻¹)
/// - `concentration_a`: \[A\] (mol/L)
///
/// Returns the rate in mol·L⁻¹·s⁻¹.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if k_reverse is not positive and
/// finite or another rate constant is negative or not finite, and
/// [`KimiyaError::InvalidConcentration`] if \[A\] is negative or not finite.
#[inline]
pub fn pre_equilibrium_rate(
    k_forward: f64,
    k_reverse: f64,
    k_slow: f64,
    concentration_a: f64,
) -> Result<f64> {
    if !is_positive(k_reverse) {
        return Err(KimiyaError::InvalidInput(
            "reverse rate constant must be positive".into(),
        ));
    }
    if !is_non_negative(k_forward) || !is_non_negative(k_slow) {
        return Err(KimiyaError::InvalidInput(
            "rate constants must be non-negative".into(),
        ));
    }
    if !is_non_negative(concentration_a) {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be non-negative".into(),
        ));
    }
    Ok(k_slow * (k_forward / k_reverse) * concentration_a)
}

// ── Michaelis-Menten enzyme kinetics ─────────────────────────────────

/// Michaelis-Menten rate: v = Vmax × \[S\] / (Km + \[S\])
///
/// - `v_max`: maximum reaction rate (mol·L⁻¹·s⁻¹)
/// - `km`: Michaelis constant, the substrate concentration at half Vmax (mol/L)
/// - `substrate`: substrate concentration \[S\] (mol/L)
///
/// Returns v in the units of `v_max`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if Vmax or Km is not positive and
/// finite, and [`KimiyaError::InvalidConcentration`] if \[S\] is negative or
/// not finite.
#[inline]
pub fn michaelis_menten(v_max: f64, km: f64, substrate: f64) -> Result<f64> {
    if !is_positive(v_max) {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if !is_positive(km) {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    if !is_non_negative(substrate) {
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
/// - slope = Km / Vmax (s)
/// - intercept = 1 / Vmax (L·s/mol)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if Vmax or Km is not positive and
/// finite.
#[inline]
pub fn lineweaver_burk(v_max: f64, km: f64) -> Result<(f64, f64)> {
    if !is_positive(v_max) {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if !is_positive(km) {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    Ok((km / v_max, 1.0 / v_max))
}

// ── Collision theory ─────────────────────────────────────────────────

/// Collision theory rate constant:
/// k = Z × p × exp(-Ea / (RT))
///
/// - `collision_frequency`: Z, the collision frequency factor, in the units of
///   the returned rate constant (e.g. L·mol⁻¹·s⁻¹ for a bimolecular reaction)
/// - `steric_factor`: p (probability of correct orientation, 0 < p ≤ 1)
/// - `activation_energy_j`: Ea in J/mol
/// - `temperature_k`: T in Kelvin
///
/// Returns k in the units of Z.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if the temperature is not
/// positive and finite, and [`KimiyaError::InvalidInput`] if the steric factor
/// is outside (0, 1], Z is not positive and finite, or Ea is not finite.
#[inline]
pub fn collision_theory_rate(
    collision_frequency: f64,
    steric_factor: f64,
    activation_energy_j: f64,
    temperature_k: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    // Written positively so that NaN fails the range check.
    if !(steric_factor > 0.0 && steric_factor <= 1.0) {
        return Err(KimiyaError::InvalidInput(
            "steric factor must be in (0, 1]".into(),
        ));
    }
    if !is_positive(collision_frequency) {
        return Err(KimiyaError::InvalidInput(
            "collision frequency must be positive".into(),
        ));
    }
    require_finite(activation_energy_j, "activation energy")?;
    Ok(collision_frequency
        * steric_factor
        * (-activation_energy_j / (GAS_CONSTANT * temperature_k)).exp())
}

// ── Transition state theory (Eyring equation) ────────────────────────

/// Boltzmann constant k_B (J/K), exact in the 2019 SI.
pub const BOLTZMANN: f64 = 1.380649e-23;

/// Eyring equation (transition state theory):
/// k = (kB × T / h) × exp(-ΔG‡ / (RT))
///
/// - `delta_g_activation_j`: ΔG‡ (Gibbs energy of activation) in J/mol
/// - `temperature_k`: T in Kelvin
///
/// Returns rate constant in s⁻¹ (transmission coefficient κ = 1).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, and [`KimiyaError::InvalidInput`] if ΔG‡ is not finite.
#[inline]
pub fn eyring_rate(delta_g_activation_j: f64, temperature_k: f64) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    require_finite(delta_g_activation_j, "Gibbs energy of activation")?;
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
/// Returns rate constant in s⁻¹ (transmission coefficient κ = 1).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, and [`KimiyaError::InvalidInput`] if ΔH‡ or ΔS‡ is not finite.
#[inline]
pub fn eyring_rate_from_activation(
    delta_h_activation_j: f64,
    delta_s_activation_j_per_k: f64,
    temperature_k: f64,
) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    require_finite(delta_h_activation_j, "enthalpy of activation")?;
    require_finite(delta_s_activation_j_per_k, "entropy of activation")?;
    let prefactor = BOLTZMANN * temperature_k / crate::spectroscopy::PLANCK;
    Ok(prefactor
        * (delta_s_activation_j_per_k / GAS_CONSTANT).exp()
        * (-delta_h_activation_j / (GAS_CONSTANT * temperature_k)).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert |actual/expected − 1| ≤ rel.
    fn assert_rel(actual: f64, expected: f64, rel: f64, what: &str) {
        assert!(
            ((actual - expected) / expected).abs() <= rel,
            "{what}: expected {expected}, got {actual}"
        );
    }

    #[test]
    fn arrhenius_basic() {
        // 1e13 · exp(−100000 / (8.314462618 · 300)), 50-digit decimal
        let k = arrhenius_rate(1e13, 100_000.0, 300.0).unwrap();
        assert_rel(k, 3.879_656_655_666_024e-5, 1e-12, "Arrhenius k");
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
    fn second_order_value() {
        // 1/(1/1.0 + 0.1·10) = 0.5
        let c = second_order_concentration(1.0, 0.1, 10.0).unwrap();
        assert_rel(c, 0.5, 1e-15, "second order");
    }

    #[test]
    fn second_order_zero_initial_is_error() {
        assert!(second_order_concentration(0.0, 0.1, 10.0).is_err());
    }

    #[test]
    fn second_order_rejects_negative_or_nan_k_and_t() {
        // k·t = −1/[A]₀ gave Ok(inf); k·t < −1/[A]₀ gave a negative concentration.
        assert!(second_order_concentration(1.0, -0.1, 10.0).is_err());
        assert!(second_order_concentration(1.0, -0.2, 10.0).is_err());
        assert!(second_order_concentration(1.0, 0.1, -10.0).is_err());
        assert!(second_order_concentration(1.0, f64::NAN, 10.0).is_err());
        assert!(second_order_concentration(1.0, 0.1, f64::INFINITY).is_err());
        assert!(second_order_concentration(f64::NAN, 0.1, 10.0).is_err());
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
        assert_eq!(k, 1e13, "zero Ea should give k=A");
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
        // 1e10 · 0.1 · exp(−50000 / (8.314462618 · 300)), 50-digit decimal
        let k = collision_theory_rate(1e10, 0.1, 50_000.0, 300.0).unwrap();
        assert_rel(k, 1.969_684_405_092_862_7, 1e-12, "collision k");
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

    #[test]
    fn collision_theory_rejects_bad_z_and_nan_steric() {
        // A negative Z used to give a negative rate constant.
        assert!(collision_theory_rate(-1e10, 0.1, 50_000.0, 300.0).is_err());
        assert!(collision_theory_rate(0.0, 0.1, 50_000.0, 300.0).is_err());
        assert!(collision_theory_rate(f64::NAN, 0.1, 50_000.0, 300.0).is_err());
        // NaN slipped through `p <= 0 || p > 1`.
        assert!(collision_theory_rate(1e10, f64::NAN, 50_000.0, 300.0).is_err());
        assert!(collision_theory_rate(1e10, 0.1, f64::NAN, 300.0).is_err());
        assert!(collision_theory_rate(1e10, 0.1, 50_000.0, f64::NAN).is_err());
    }

    // ── Eyring equation ──────────────────────────────────────────────

    #[test]
    fn eyring_basic() {
        // k_B·T/h = 6.212437991620e12 s⁻¹ at 298.15 K;
        // × exp(−80000 / (R·298.15)) = 0.059960941088 s⁻¹ (50-digit decimal)
        let k = eyring_rate(80_000.0, 298.15).unwrap();
        assert_rel(k, 5.996_094_108_818_953e-2, 1e-12, "Eyring k");
        let prefactor = eyring_rate(0.0, 298.15).unwrap();
        assert_rel(prefactor, 6.212_437_991_620_116e12, 1e-12, "k_B·T/h");
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
            (k_combined - k_decomposed).abs() / k_combined < 1e-12,
            "combined ({k_combined}) should match decomposed ({k_decomposed})"
        );
        assert_rel(
            k_decomposed,
            1.441_830_328_253_839_6e-4,
            1e-12,
            "Eyring ΔH‡/ΔS‡",
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
    fn zero_order_half_life_rejects_bad_initial() {
        assert!(zero_order_half_life(0.0, 0.1).is_err());
        assert!(zero_order_half_life(-1.0, 0.1).is_err());
        assert!(zero_order_half_life(f64::NAN, 0.1).is_err());
        assert!(zero_order_half_life(1.0, f64::NAN).is_err());
    }

    #[test]
    fn zero_order_concentration_propagates_nan() {
        // f64::max(NaN, 0.0) silently returned 0.0
        assert!(zero_order_concentration(f64::NAN, 0.1, 1.0).is_nan());
        assert!(zero_order_concentration(1.0, f64::NAN, 1.0).is_nan());
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
    fn nth_order_one_or_non_finite_is_error() {
        assert!(nth_order_half_life(1.0, 0.1, 1.0).is_err());
        assert!(nth_order_half_life(1.0, 0.1, 1.0 + 1e-13).is_err());
        assert!(nth_order_half_life(1.0, 0.1, f64::NAN).is_err());
        assert!(nth_order_half_life(1.0, 0.1, f64::INFINITY).is_err());
        assert!(nth_order_half_life(1.0, 0.1, f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn nth_order_below_one_is_valid() {
        // n = 0 reproduces the zero-order half-life [A]₀/(2k)
        let t0 = nth_order_half_life(1.0, 0.1, 0.0).unwrap();
        assert_rel(t0, zero_order_half_life(1.0, 0.1).unwrap(), 1e-15, "n = 0");
        // n = ½: t½ = 2(1 − 1/√2)·√[A]₀ / k
        let t_half = nth_order_half_life(1.0, 0.1, 0.5).unwrap();
        assert_rel(t_half, 5.85786437626905, 1e-14, "n = 1/2");
    }

    #[test]
    fn nth_order_third_order_value() {
        // n = 3: t½ = 3 / (2k[A]₀²) = 3 / (2 · 0.1 · 4) = 3.75
        let t = nth_order_half_life(2.0, 0.1, 3.0).unwrap();
        assert_rel(t, 3.75, 1e-15, "n = 3");
    }

    #[test]
    fn nth_order_near_one_approaches_first_order() {
        // (2^δ − 1)/δ → ln 2; powf(2, δ) − 1 lost ~7 digits at δ = 1e-9.
        let t = nth_order_half_life(1.0, 0.1, 1.0 + 1e-9).unwrap();
        let t1 = half_life_first_order(0.1).unwrap();
        assert_rel(t, t1, 1e-9, "n → 1");
    }

    #[test]
    fn nth_order_overflow_is_error() {
        assert!(nth_order_half_life(1e-300, 0.1, 10.0).is_err());
    }

    // ── Steady-state ─────────────────────────────────────────────────

    #[test]
    fn steady_state_basic() {
        // A →(0.1) B →(1.0) C, [A]=1.0 → [B]_ss = 0.1*1.0/1.0 = 0.1
        let b_ss = steady_state_intermediate(0.1, 1.0, 1.0).unwrap();
        assert!((b_ss - 0.1).abs() < f64::EPSILON);
    }

    #[test]
    fn steady_state_rate_basic() {
        let rate = steady_state_rate(0.05, 2.0);
        assert!((rate - 0.1).abs() < f64::EPSILON);
    }

    #[test]
    fn pre_equilibrium_basic() {
        // K_eq = k1/k-1 = 10/2 = 5, rate = k2 × 5 × [A] = 0.1 × 5 × 1.0 = 0.5
        let rate = pre_equilibrium_rate(10.0, 2.0, 0.1, 1.0).unwrap();
        assert!((rate - 0.5).abs() < 1e-10);
    }

    #[test]
    fn pre_equilibrium_zero_kreverse_is_error() {
        assert!(pre_equilibrium_rate(1.0, 0.0, 0.1, 1.0).is_err());
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        assert!(arrhenius_rate(1e13, 50_000.0, nan).is_err());
        assert!(arrhenius_rate(1e13, 50_000.0, inf).is_err());
        assert!(arrhenius_rate(nan, 50_000.0, 300.0).is_err());
        assert!(arrhenius_rate(1e13, nan, 300.0).is_err());
        assert!(half_life_first_order(nan).is_err());
        assert!(half_life_first_order(inf).is_err());
        assert!(nth_order_half_life(nan, 0.1, 2.0).is_err());
        assert!(nth_order_half_life(1.0, nan, 2.0).is_err());
        assert!(steady_state_intermediate(0.1, 1.0, nan).is_err());
        assert!(steady_state_intermediate(nan, 1.0, 1.0).is_err());
        assert!(steady_state_intermediate(0.1, nan, 1.0).is_err());
        assert!(pre_equilibrium_rate(1.0, nan, 0.1, 1.0).is_err());
        assert!(pre_equilibrium_rate(nan, 2.0, 0.1, 1.0).is_err());
        assert!(pre_equilibrium_rate(1.0, 2.0, nan, 1.0).is_err());
        assert!(pre_equilibrium_rate(1.0, 2.0, 0.1, nan).is_err());
        assert!(michaelis_menten(nan, 5.0, 1.0).is_err());
        assert!(michaelis_menten(100.0, nan, 1.0).is_err());
        assert!(michaelis_menten(100.0, 5.0, nan).is_err());
        assert!(michaelis_menten(100.0, 5.0, inf).is_err());
        assert!(lineweaver_burk(nan, 5.0).is_err());
        assert!(lineweaver_burk(100.0, inf).is_err());
        assert!(eyring_rate(80_000.0, nan).is_err());
        assert!(eyring_rate(nan, 298.15).is_err());
        assert!(eyring_rate_from_activation(nan, -50.0, 298.15).is_err());
        assert!(eyring_rate_from_activation(80_000.0, nan, 298.15).is_err());
        assert!(eyring_rate_from_activation(80_000.0, -50.0, inf).is_err());
    }
}
