//! Intermolecular potentials — Lennard-Jones, Morse, Coulomb.

use crate::error::{KimiyaError, Result};

/// Coulomb constant k_e (N·m²/C²).
pub const COULOMB_CONSTANT: f64 = 8.987_551_792_3e9;

/// Lennard-Jones 6-12 potential: V(r) = 4ε\[(σ/r)¹² - (σ/r)⁶\]
///
/// - `epsilon`: well depth ε (J)
/// - `sigma`: collision diameter σ (m)
/// - `r`: separation distance (m)
///
/// # Errors
///
/// Returns error if sigma or r is not positive.
#[inline]
pub fn lennard_jones(epsilon: f64, sigma: f64, r: f64) -> Result<f64> {
    if sigma <= 0.0 {
        return Err(KimiyaError::InvalidInput("sigma must be positive".into()));
    }
    if r <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "distance must be positive".into(),
        ));
    }
    let sr = sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    Ok(4.0 * epsilon * (sr12 - sr6))
}

/// Lennard-Jones equilibrium distance: r_min = 2^(1/6) × σ ≈ 1.122σ
#[must_use]
#[inline]
pub fn lj_equilibrium_distance(sigma: f64) -> f64 {
    sigma * 2.0_f64.powf(1.0 / 6.0)
}

/// Morse potential: V(r) = D_e × (1 - exp(-a(r - r_e)))²
///
/// - `d_e`: well depth (J)
/// - `a`: width parameter (m⁻¹), controls potential well width
/// - `r_e`: equilibrium bond distance (m)
/// - `r`: separation distance (m)
///
/// # Errors
///
/// Returns error if d_e or a is not positive.
#[inline]
pub fn morse(d_e: f64, a: f64, r_e: f64, r: f64) -> Result<f64> {
    if d_e <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "well depth must be positive".into(),
        ));
    }
    if a <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "width parameter must be positive".into(),
        ));
    }
    let x = 1.0 - (-a * (r - r_e)).exp();
    Ok(d_e * x * x)
}

/// Morse potential width parameter from spectroscopic data:
/// a = ω_e × √(2π²μ / D_e)
///
/// - `omega_e`: harmonic vibrational frequency (rad/s)
/// - `reduced_mass`: μ (kg)
/// - `d_e`: well depth (J)
///
/// # Errors
///
/// Returns error if parameters are not positive.
#[inline]
pub fn morse_width_parameter(omega_e: f64, reduced_mass: f64, d_e: f64) -> Result<f64> {
    if d_e <= 0.0 || reduced_mass <= 0.0 || omega_e <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "all parameters must be positive".into(),
        ));
    }
    Ok(omega_e * (2.0 * std::f64::consts::PI * std::f64::consts::PI * reduced_mass / d_e).sqrt())
}

/// Coulomb potential between two point charges: V(r) = k_e × q₁ × q₂ / r
///
/// - `q1`, `q2`: charges (C)
/// - `r`: separation distance (m)
///
/// Returns energy in Joules.
///
/// # Errors
///
/// Returns error if distance is not positive.
#[inline]
pub fn coulomb(q1: f64, q2: f64, r: f64) -> Result<f64> {
    if r <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "distance must be positive".into(),
        ));
    }
    Ok(COULOMB_CONSTANT * q1 * q2 / r)
}

/// Coulomb force between two point charges: F = k_e × q₁ × q₂ / r²
///
/// Positive = repulsion, negative = attraction.
///
/// # Errors
///
/// Returns error if distance is not positive.
#[inline]
pub fn coulomb_force(q1: f64, q2: f64, r: f64) -> Result<f64> {
    if r <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "distance must be positive".into(),
        ));
    }
    Ok(COULOMB_CONSTANT * q1 * q2 / (r * r))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lj_at_equilibrium_is_negative() {
        let sigma = 3.4e-10; // Ar
        let r = lj_equilibrium_distance(sigma);
        let v = lennard_jones(1.65e-21, sigma, r).unwrap();
        assert!(v < 0.0, "V at r_min should be negative (bound), got {v}");
    }

    #[test]
    fn lj_at_sigma_is_zero() {
        let v = lennard_jones(1.0, 1.0, 1.0).unwrap();
        assert!(v.abs() < f64::EPSILON, "V(σ) should be 0, got {v}");
    }

    #[test]
    fn lj_repulsive_at_short_range() {
        let v = lennard_jones(1.0, 1.0, 0.5).unwrap();
        assert!(v > 0.0, "V at r<σ should be repulsive");
    }

    #[test]
    fn lj_approaches_zero_at_long_range() {
        let v = lennard_jones(1.0, 1.0, 100.0).unwrap();
        assert!(v.abs() < 1e-10, "V at large r should approach 0");
    }

    #[test]
    fn lj_zero_distance_is_error() {
        assert!(lennard_jones(1.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn morse_at_equilibrium_is_zero() {
        let v = morse(5.0, 1.0, 1.5, 1.5).unwrap();
        assert!(v.abs() < f64::EPSILON, "V(r_e) should be 0, got {v}");
    }

    #[test]
    fn morse_at_infinity_is_de() {
        let v = morse(5.0, 1.0, 1.5, 100.0).unwrap();
        assert!((v - 5.0).abs() < 0.01, "V(∞) should approach D_e, got {v}");
    }

    #[test]
    fn morse_symmetric_near_equilibrium() {
        let d = 5.0;
        let a = 2.0;
        let re = 1.5;
        let dr = 0.1;
        let v_plus = morse(d, a, re, re + dr).unwrap();
        let v_minus = morse(d, a, re, re - dr).unwrap();
        // Near equilibrium, Morse is approximately harmonic → nearly symmetric
        assert!(
            (v_plus - v_minus).abs() < 0.5 * v_plus.max(v_minus),
            "should be roughly symmetric near r_e"
        );
    }

    #[test]
    fn coulomb_opposite_charges_attractive() {
        let e = 1.602176634e-19;
        let v = coulomb(e, -e, 1e-10).unwrap();
        assert!(v < 0.0, "opposite charges should attract");
    }

    #[test]
    fn coulomb_same_charges_repulsive() {
        let e = 1.602176634e-19;
        let v = coulomb(e, e, 1e-10).unwrap();
        assert!(v > 0.0, "same charges should repel");
    }

    #[test]
    fn coulomb_inverse_r() {
        let e = 1.602176634e-19;
        let v1 = coulomb(e, e, 1e-10).unwrap();
        let v2 = coulomb(e, e, 2e-10).unwrap();
        assert!((v1 / v2 - 2.0).abs() < 0.001, "V should scale as 1/r");
    }

    #[test]
    fn coulomb_zero_distance_is_error() {
        assert!(coulomb(1.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn coulomb_force_basic() {
        let e = 1.602176634e-19;
        let f = coulomb_force(e, e, 1e-10).unwrap();
        // F ≈ 8.99e9 × (1.6e-19)² / (1e-10)² ≈ 2.3e-8 N
        assert!(
            f > 2e-8 && f < 3e-8,
            "proton-proton force should be ~2.3e-8 N, got {f}"
        );
    }
}
