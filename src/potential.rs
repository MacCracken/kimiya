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

// ── Molecular geometry from 3D coordinates ───────────────────────────

/// Distance between two 3D points: d = √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
#[must_use]
#[inline]
pub fn distance_3d(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dz = p2[2] - p1[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Bond angle between three atoms (A-B-C) in radians.
///
/// B is the central atom.
///
/// # Errors
///
/// Returns error if any two atoms coincide.
pub fn bond_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Result<f64> {
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
    let mag_ba = (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]).sqrt();
    let mag_bc = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]).sqrt();
    if mag_ba < 1e-15 || mag_bc < 1e-15 {
        return Err(KimiyaError::InvalidInput("atoms must not coincide".into()));
    }
    let cos_theta = (dot / (mag_ba * mag_bc)).clamp(-1.0, 1.0);
    Ok(cos_theta.acos())
}

/// Dihedral angle between four atoms (A-B-C-D) in radians.
///
/// The dihedral is the angle between the ABC and BCD planes.
///
/// # Errors
///
/// Returns error if any three consecutive atoms are collinear.
pub fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Result<f64> {
    let b1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let b2 = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let b3 = [d[0] - c[0], d[1] - c[1], d[2] - c[2]];

    // n1 = b1 × b2, n2 = b2 × b3
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let mag_n1 = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
    let mag_n2 = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
    if mag_n1 < 1e-15 || mag_n2 < 1e-15 {
        return Err(KimiyaError::InvalidInput(
            "collinear atoms — dihedral undefined".into(),
        ));
    }
    let cos_phi = (dot3(n1, n2) / (mag_n1 * mag_n2)).clamp(-1.0, 1.0);
    // Sign from (n1 × n2) · b2
    let sign = dot3(cross(n1, n2), b2);
    let angle = cos_phi.acos();
    if sign < 0.0 { Ok(-angle) } else { Ok(angle) }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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
        assert!(
            f > 2e-8 && f < 3e-8,
            "proton-proton force should be ~2.3e-8 N, got {f}"
        );
    }

    // ── 3D geometry ──────────────────────────────────────────────────

    #[test]
    fn distance_3d_basic() {
        let d = distance_3d([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < f64::EPSILON);
    }

    #[test]
    fn distance_3d_same_point() {
        let d = distance_3d([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(d.abs() < f64::EPSILON);
    }

    #[test]
    fn bond_angle_right_angle() {
        // 90° angle
        let angle = bond_angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        assert!(
            (angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "should be 90°, got {:.1}°",
            angle.to_degrees()
        );
    }

    #[test]
    fn bond_angle_linear() {
        let angle = bond_angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]).unwrap();
        assert!(
            (angle - std::f64::consts::PI).abs() < 1e-10,
            "should be 180°, got {:.1}°",
            angle.to_degrees()
        );
    }

    #[test]
    fn bond_angle_coincident_is_error() {
        assert!(bond_angle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn dihedral_angle_planar() {
        // All in xy-plane → dihedral ≈ 0 or π
        let angle = dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(
            angle.abs() < 0.1 || (angle.abs() - std::f64::consts::PI).abs() < 0.1,
            "planar dihedral should be ~0 or ~180°, got {:.1}°",
            angle.to_degrees()
        );
    }

    #[test]
    fn dihedral_angle_perpendicular() {
        // 90° dihedral
        let angle = dihedral_angle(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        )
        .unwrap();
        assert!(
            (angle.abs() - std::f64::consts::FRAC_PI_2).abs() < 0.1,
            "should be ~90°, got {:.1}°",
            angle.to_degrees()
        );
    }
}
