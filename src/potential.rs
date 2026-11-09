//! Intermolecular potentials — Lennard-Jones, Morse, Coulomb — and molecular
//! geometry (distances, bond and dihedral angles) from 3D coordinates.

use crate::error::{KimiyaError, Result};
use crate::validate::is_positive;

/// Coulomb constant k_e = 1/(4πε₀) (N·m²/C²), CODATA 2018.
pub const COULOMB_CONSTANT: f64 = 8.987_551_792_3e9;

/// Relative resolution below which a bond vector is treated as zero
/// (atoms coincide): |b| ≤ 10⁻¹² × the coordinate magnitude.
const COINCIDENT_REL_TOL: f64 = 1e-12;

/// Squared sine of the bond angle below which three atoms count as
/// collinear in [`dihedral_angle`]: |b₁×b₂|² ≤ (10⁻¹⁰)²·|b₁|²·|b₂|².
const COLLINEAR_SIN_SQ_TOL: f64 = 1e-20;

/// Lennard-Jones 6-12 potential: V(r) = 4ε\[(σ/r)¹² - (σ/r)⁶\]
///
/// - `epsilon`: well depth ε (J)
/// - `sigma`: collision diameter σ (m)
/// - `r`: separation distance (m)
///
/// Returns energy in the unit of `epsilon` (J); any consistent length unit
/// may be used for `sigma` and `r`.
///
/// # Errors
///
/// Returns error if epsilon, sigma or r is not positive and finite.
#[inline]
pub fn lennard_jones(epsilon: f64, sigma: f64, r: f64) -> Result<f64> {
    if !is_positive(epsilon) {
        return Err(KimiyaError::InvalidInput(
            "well depth must be positive".into(),
        ));
    }
    if !is_positive(sigma) {
        return Err(KimiyaError::InvalidInput("sigma must be positive".into()));
    }
    if !is_positive(r) {
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
///
/// Returns r_min in the unit of `sigma`.
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
/// Returns energy in the unit of `d_e` (J), zero at r = r_e.
///
/// # Errors
///
/// Returns error if d_e or a is not positive and finite.
#[inline]
pub fn morse(d_e: f64, a: f64, r_e: f64, r: f64) -> Result<f64> {
    if !is_positive(d_e) {
        return Err(KimiyaError::InvalidInput(
            "well depth must be positive".into(),
        ));
    }
    if !is_positive(a) {
        return Err(KimiyaError::InvalidInput(
            "width parameter must be positive".into(),
        ));
    }
    let x = 1.0 - (-a * (r - r_e)).exp();
    Ok(d_e * x * x)
}

/// Morse potential width parameter from spectroscopic data:
/// a = ω_e × √(μ / (2·D_e))
///
/// Follows from matching the Morse curvature at r_e, k = 2·D_e·a², to the
/// harmonic force constant k = μ·ω_e².
///
/// - `omega_e`: harmonic vibrational **angular** frequency ω_e (rad/s), as
///   returned by [`crate::spectroscopy::vibrational_frequency`]
///   (ω = 2π·ν = 2π·c·ω̃ for a wavenumber ω̃)
/// - `reduced_mass`: μ (kg)
/// - `d_e`: well depth (J)
///
/// Returns a in m⁻¹.
///
/// # Errors
///
/// Returns error if parameters are not positive and finite.
#[inline]
pub fn morse_width_parameter(omega_e: f64, reduced_mass: f64, d_e: f64) -> Result<f64> {
    if !is_positive(d_e) || !is_positive(reduced_mass) || !is_positive(omega_e) {
        return Err(KimiyaError::InvalidInput(
            "all parameters must be positive".into(),
        ));
    }
    Ok(omega_e * (reduced_mass / (2.0 * d_e)).sqrt())
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
/// Returns error if distance is not positive and finite.
#[inline]
pub fn coulomb(q1: f64, q2: f64, r: f64) -> Result<f64> {
    if !is_positive(r) {
        return Err(KimiyaError::InvalidInput(
            "distance must be positive".into(),
        ));
    }
    Ok(COULOMB_CONSTANT * q1 * q2 / r)
}

/// Coulomb force between two point charges: F = k_e × q₁ × q₂ / r²
///
/// - `q1`, `q2`: charges (C)
/// - `r`: separation distance (m)
///
/// Returns force in newtons. Positive = repulsion, negative = attraction.
///
/// # Errors
///
/// Returns error if distance is not positive and finite.
#[inline]
pub fn coulomb_force(q1: f64, q2: f64, r: f64) -> Result<f64> {
    if !is_positive(r) {
        return Err(KimiyaError::InvalidInput(
            "distance must be positive".into(),
        ));
    }
    Ok(COULOMB_CONSTANT * q1 * q2 / (r * r))
}

// ── Molecular geometry from 3D coordinates ───────────────────────────

/// Distance between two 3D points: d = √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
///
/// Returns the distance in the unit of the coordinates.
#[must_use]
#[inline]
pub fn distance_3d(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dz = p2[2] - p1[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Bond angle between three atoms (A-B-C) in radians, in \[0, π\].
///
/// B is the central atom. Coordinates may be in any length unit (Å, m, …);
/// the result is scale-invariant.
///
/// # Errors
///
/// Returns error if any coordinate is not finite or if A or C coincides
/// with B (bond length zero, or below floating-point resolution relative to
/// the coordinates).
pub fn bond_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Result<f64> {
    let ba = sub3(a, b);
    let bc = sub3(c, b);
    let ba_sq = dot3(ba, ba);
    let bc_sq = dot3(bc, bc);
    let tol_sq = COINCIDENT_REL_TOL * COINCIDENT_REL_TOL;
    let b_sq = dot3(b, b);
    if !is_positive(ba_sq)
        || !is_positive(bc_sq)
        || ba_sq <= tol_sq * dot3(a, a).max(b_sq)
        || bc_sq <= tol_sq * dot3(c, c).max(b_sq)
    {
        return Err(KimiyaError::InvalidInput("atoms must not coincide".into()));
    }
    // atan2(|ba × bc|, ba · bc) is accurate near 0 and π, unlike acos.
    let cr = cross(ba, bc);
    Ok(dot3(cr, cr).sqrt().atan2(dot3(ba, bc)))
}

/// Dihedral (torsion) angle between four atoms (A-B-C-D) in radians,
/// in (−π, π\].
///
/// The dihedral is the angle between the ABC and BCD planes, with the IUPAC
/// sign convention: positive when, viewed along B→C, the bond B–A must be
/// rotated clockwise to eclipse C–D. Planar trans (anti) is π, cis is 0.
/// Coordinates may be in any length unit; the result is scale-invariant.
///
/// # Errors
///
/// Returns error if any coordinate is not finite, any two consecutive atoms
/// coincide, or A-B-C or B-C-D is collinear (bond angle within ~10⁻¹⁰ rad
/// of 0 or π).
pub fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Result<f64> {
    let b1 = sub3(b, a);
    let b2 = sub3(c, b);
    let b3 = sub3(d, c);

    // n1 = b1 × b2, n2 = b2 × b3
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let n1_sq = dot3(n1, n1);
    let n2_sq = dot3(n2, n2);
    let b2_sq = dot3(b2, b2);
    // |b1 × b2| = |b1|·|b2|·sin θ: compare the sine, not an absolute length,
    // so the check works for Å and SI-metre coordinates alike.
    if !is_positive(n1_sq)
        || !is_positive(n2_sq)
        || n1_sq <= COLLINEAR_SIN_SQ_TOL * dot3(b1, b1) * b2_sq
        || n2_sq <= COLLINEAR_SIN_SQ_TOL * b2_sq * dot3(b3, b3)
    {
        return Err(KimiyaError::InvalidInput(
            "collinear atoms — dihedral undefined".into(),
        ));
    }
    // φ = atan2(|b2|·b1·n2, n1·n2) (IUPAC sign).
    let phi = (b2_sq.sqrt() * dot3(b1, n2)).atan2(dot3(n1, n2));
    // atan2 yields −π for a trans geometry whose y-term is −0.0; report (−π, π].
    if phi <= -std::f64::consts::PI {
        Ok(std::f64::consts::PI)
    } else {
        Ok(phi)
    }
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
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
    fn lj_invalid_epsilon_is_error() {
        // Consistent with `morse` rejecting D_e ≤ 0.
        assert!(lennard_jones(0.0, 1.0, 1.0).is_err());
        assert!(lennard_jones(-1.0, 1.0, 1.0).is_err());
        assert!(lennard_jones(f64::NAN, 1.0, 1.0).is_err());
        assert!(lennard_jones(f64::INFINITY, 1.0, 1.0).is_err());
    }

    #[test]
    fn lj_nan_inputs_are_errors() {
        assert!(lennard_jones(1.0, f64::NAN, 1.0).is_err());
        assert!(lennard_jones(1.0, 1.0, f64::NAN).is_err());
        assert!(lennard_jones(1.0, 1.0, f64::INFINITY).is_err());
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
    fn morse_nan_inputs_are_errors() {
        assert!(morse(f64::NAN, 1.0, 1.5, 1.5).is_err());
        assert!(morse(5.0, f64::NAN, 1.5, 1.5).is_err());
        assert!(morse(f64::INFINITY, 1.0, 1.5, 1.5).is_err());
    }

    #[test]
    fn morse_width_parameter_hcl() {
        // H³⁵Cl: ω̃_e = 2990.9 cm⁻¹ → ω_e = 2πcω̃ = 5.634e14 rad/s,
        // μ = 0.9796 u = 1.6267e-27 kg, D_e = 4.618 eV = 7.399e-19 J.
        // a = ω√(μ/2D_e) = 1.868e10 m⁻¹ (literature ≈ 1.87 Å⁻¹).
        let a = morse_width_parameter(5.634e14, 1.6267e-27, 7.399e-19).unwrap();
        assert!(
            (a - 1.868e10).abs() < 0.005e10,
            "HCl a ≈ 1.87e10 m⁻¹, got {a}"
        );
    }

    #[test]
    fn morse_width_parameter_matches_force_constant() {
        // Morse curvature 2·D_e·a² must equal the harmonic k = μω² that
        // `spectroscopy::vibrational_frequency` inverts.
        let (k, mu, d_e) = (516.0, 1.6267e-27, 7.399e-19);
        let omega = crate::spectroscopy::vibrational_frequency(k, mu).unwrap();
        let a = morse_width_parameter(omega, mu, d_e).unwrap();
        assert!(((2.0 * d_e * a * a) / k - 1.0).abs() < 1e-12);
    }

    #[test]
    fn morse_width_parameter_nan_is_error() {
        assert!(morse_width_parameter(f64::NAN, 1.0, 1.0).is_err());
        assert!(morse_width_parameter(1.0, f64::NAN, 1.0).is_err());
        assert!(morse_width_parameter(1.0, 1.0, f64::INFINITY).is_err());
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
    fn coulomb_nan_distance_is_error() {
        assert!(coulomb(1.0, 1.0, f64::NAN).is_err());
        assert!(coulomb_force(1.0, 1.0, f64::NAN).is_err());
        assert!(coulomb_force(1.0, 1.0, f64::NEG_INFINITY).is_err());
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
        // Differs from B only by rounding noise relative to the coordinates.
        assert!(bond_angle([1.0 + 1e-15, 1.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]).is_err());
    }

    #[test]
    fn bond_angle_nan_is_error() {
        assert!(bond_angle([f64::NAN, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
        assert!(bond_angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]).is_err());
    }

    #[test]
    fn bond_angle_metre_scale() {
        // Water in SI metres: O–H 0.9572 Å, H–O–H 104.52°.
        let r = 0.9572e-10;
        let theta = 104.52_f64.to_radians();
        let h1 = [r, 0.0, 0.0];
        let h2 = [r * theta.cos(), r * theta.sin(), 0.0];
        let angle = bond_angle(h1, [0.0; 3], h2).unwrap();
        assert!((angle - theta).abs() < 1e-12, "got {}", angle.to_degrees());
        // Femtometre scale works too.
        let angle = bond_angle([1e-15, 0.0, 0.0], [0.0; 3], [0.0, 1e-15, 0.0]).unwrap();
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn bond_angle_small_angle_is_accurate() {
        // acos(cos θ) collapses to 0 for θ ≲ 1e-8; atan2 does not.
        let theta = 1e-9_f64;
        let angle = bond_angle(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [theta.cos(), theta.sin(), 0.0],
        )
        .unwrap();
        assert!((angle / theta - 1.0).abs() < 1e-6, "got {angle}");
    }

    #[test]
    fn dihedral_angle_planar_trans_is_pi() {
        // All in xy-plane, A and D on opposite sides of B–C → trans, φ = π.
        let angle = dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
        )
        .unwrap();
        assert_eq!(angle, std::f64::consts::PI, "got {}°", angle.to_degrees());
        // Mirror image through the B–C axis: atan2 sees −0.0 → still +π.
        let angle = dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, -1.0, 0.0],
            [3.0, -1.0, 0.0],
        )
        .unwrap();
        assert_eq!(angle, std::f64::consts::PI, "got {}°", angle.to_degrees());
    }

    #[test]
    fn dihedral_angle_planar_cis_is_zero() {
        let angle = dihedral_angle(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(angle.abs() < 1e-15, "got {}°", angle.to_degrees());
    }

    #[test]
    fn dihedral_angle_perpendicular() {
        // Viewed along B→C (+y), B–A points along +x and C–D along +z:
        // eclipsing needs a counter-clockwise turn → IUPAC −90°.
        let angle = dihedral_angle(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        )
        .unwrap();
        assert!(
            (angle + std::f64::consts::FRAC_PI_2).abs() < 1e-12,
            "should be −90°, got {:.1}°",
            angle.to_degrees()
        );
        // Mirror image (D along −z) → +90°.
        let angle = dihedral_angle(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, -1.0],
        )
        .unwrap();
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn dihedral_angle_metre_scale() {
        // Butane-like gauche geometry: Å vs SI metres must agree. The old
        // absolute 1e-15 check rejected every metre-scale input as collinear.
        let pts = [
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [2.5, 0.0, 1.0],
        ];
        let angstrom = dihedral_angle(pts[0], pts[1], pts[2], pts[3]).unwrap();
        let m = pts.map(|p| p.map(|x| x * 1e-10));
        let metre = dihedral_angle(m[0], m[1], m[2], m[3]).unwrap();
        assert!((angstrom - metre).abs() < 1e-12, "{angstrom} vs {metre}");
        // Viewed along B→C (+x): B–A at 3 o'clock, C–D at 6 o'clock → +90°.
        assert!((angstrom - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn dihedral_angle_collinear_is_error() {
        assert!(
            dihedral_angle(
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
            )
            .is_err()
        );
        // Collinear in metres too.
        assert!(
            dihedral_angle(
                [0.0, 0.0, 0.0],
                [1e-10, 0.0, 0.0],
                [2e-10, 1e-10, 0.0],
                [3e-10, 2e-10, 0.0],
            )
            .is_err()
        );
        assert!(dihedral_angle([f64::NAN; 3], [0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]).is_err());
    }
}
