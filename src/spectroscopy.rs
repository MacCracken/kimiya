//! Spectroscopy — Beer-Lambert law, Bohr model, photon energy conversions.

use crate::error::{KimiyaError, Result};

// ── Physical constants ───────────────────────────────────────────────

/// Planck's constant (J·s).
pub const PLANCK: f64 = 6.62607015e-34;

/// Speed of light in vacuum (m/s).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Rydberg constant (m⁻¹).
pub const RYDBERG: f64 = 1.097_373_156_816_0e7;

/// Rydberg energy (eV) — ionization energy of hydrogen.
pub const RYDBERG_EV: f64 = 13.605693122994;

// ── Beer-Lambert law ─────────────────────────────────────────────────

/// Beer-Lambert law: A = ε × l × c
///
/// - `molar_absorptivity`: ε in L/(mol·cm)
/// - `path_length`: l in cm
/// - `concentration`: c in mol/L
///
/// Returns absorbance (dimensionless).
#[must_use]
#[inline]
pub fn absorbance(molar_absorptivity: f64, path_length: f64, concentration: f64) -> f64 {
    molar_absorptivity * path_length * concentration
}

/// Transmittance from absorbance: T = 10^(-A)
#[must_use]
#[inline]
pub fn transmittance(absorbance: f64) -> f64 {
    10.0_f64.powf(-absorbance)
}

/// Absorbance from transmittance: A = -log₁₀(T)
///
/// # Errors
///
/// Returns error if transmittance is not in (0, 1].
#[inline]
pub fn absorbance_from_transmittance(transmittance: f64) -> Result<f64> {
    if transmittance <= 0.0 || transmittance > 1.0 {
        return Err(KimiyaError::InvalidInput(
            "transmittance must be in (0, 1]".into(),
        ));
    }
    Ok(-transmittance.log10())
}

/// Concentration from absorbance: c = A / (ε × l)
///
/// # Errors
///
/// Returns error if molar absorptivity or path length is not positive.
#[inline]
pub fn concentration_from_absorbance(
    absorbance: f64,
    molar_absorptivity: f64,
    path_length: f64,
) -> Result<f64> {
    if molar_absorptivity <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "molar absorptivity must be positive".into(),
        ));
    }
    if path_length <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "path length must be positive".into(),
        ));
    }
    Ok(absorbance / (molar_absorptivity * path_length))
}

// ── Photon energy/wavelength/frequency ───────────────────────────────

/// Photon energy from frequency: E = hν
#[must_use]
#[inline]
pub fn photon_energy_from_frequency(frequency_hz: f64) -> f64 {
    PLANCK * frequency_hz
}

/// Photon energy from wavelength: E = hc/λ
///
/// # Errors
///
/// Returns error if wavelength is not positive.
#[inline]
pub fn photon_energy_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if wavelength_m <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "wavelength must be positive".into(),
        ));
    }
    Ok(PLANCK * SPEED_OF_LIGHT / wavelength_m)
}

/// Wavelength from photon energy: λ = hc/E
///
/// # Errors
///
/// Returns error if energy is not positive.
#[inline]
pub fn wavelength_from_energy(energy_j: f64) -> Result<f64> {
    if energy_j <= 0.0 {
        return Err(KimiyaError::InvalidInput("energy must be positive".into()));
    }
    Ok(PLANCK * SPEED_OF_LIGHT / energy_j)
}

/// Frequency from wavelength: ν = c/λ
///
/// # Errors
///
/// Returns error if wavelength is not positive.
#[inline]
pub fn frequency_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if wavelength_m <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "wavelength must be positive".into(),
        ));
    }
    Ok(SPEED_OF_LIGHT / wavelength_m)
}

/// Convert energy in Joules to electron-volts.
#[must_use]
#[inline]
pub fn joules_to_ev(energy_j: f64) -> f64 {
    energy_j / 1.602176634e-19
}

/// Convert energy in electron-volts to Joules.
#[must_use]
#[inline]
pub fn ev_to_joules(energy_ev: f64) -> f64 {
    energy_ev * 1.602176634e-19
}

// ── Bohr model ───────────────────────────────────────────────────────

/// Bohr model energy level for hydrogen-like atoms: E_n = -13.6 × Z² / n² (eV)
///
/// - `z`: nuclear charge (1 for hydrogen)
/// - `n`: principal quantum number (≥ 1)
///
/// # Errors
///
/// Returns error if n is zero or Z is zero.
#[inline]
pub fn bohr_energy_level(z: u8, n: u32) -> Result<f64> {
    if z == 0 {
        return Err(KimiyaError::InvalidInput(
            "nuclear charge Z must be positive".into(),
        ));
    }
    if n == 0 {
        return Err(KimiyaError::InvalidInput(
            "quantum number n must be at least 1".into(),
        ));
    }
    Ok(-RYDBERG_EV * (z as f64) * (z as f64) / (n as f64 * n as f64))
}

/// Energy of a photon emitted in a transition from n_upper to n_lower.
///
/// Returns positive energy (the photon energy released).
///
/// # Errors
///
/// Returns error if n_lower or n_upper is zero.
#[inline]
pub fn bohr_transition_energy(z: u8, n_lower: u32, n_upper: u32) -> Result<f64> {
    let e_lower = bohr_energy_level(z, n_lower)?;
    let e_upper = bohr_energy_level(z, n_upper)?;
    Ok(e_upper - e_lower) // upper is less negative, so this is positive for emission
}

/// Wavelength of a photon from a hydrogen-like transition (Rydberg formula).
///
/// 1/λ = R × Z² × (1/n₁² - 1/n₂²) where n₂ > n₁
///
/// Returns wavelength in meters.
///
/// # Errors
///
/// Returns error if quantum numbers are invalid or n_upper ≤ n_lower.
pub fn rydberg_wavelength(z: u8, n_lower: u32, n_upper: u32) -> Result<f64> {
    if z == 0 {
        return Err(KimiyaError::InvalidInput("Z must be positive".into()));
    }
    if n_lower == 0 || n_upper == 0 {
        return Err(KimiyaError::InvalidInput(
            "quantum numbers must be at least 1".into(),
        ));
    }
    if n_upper <= n_lower {
        return Err(KimiyaError::InvalidInput(
            "n_upper must be greater than n_lower".into(),
        ));
    }
    let inv_lambda = RYDBERG
        * (z as f64)
        * (z as f64)
        * (1.0 / (n_lower as f64 * n_lower as f64) - 1.0 / (n_upper as f64 * n_upper as f64));
    Ok(1.0 / inv_lambda)
}

/// Named hydrogen spectral series.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SpectralSeries {
    /// n_lower = 1 (UV)
    Lyman,
    /// n_lower = 2 (visible)
    Balmer,
    /// n_lower = 3 (IR)
    Paschen,
    /// n_lower = 4 (IR)
    Brackett,
    /// n_lower = 5 (far IR)
    Pfund,
}

impl SpectralSeries {
    /// Lower quantum number for this series.
    #[must_use]
    pub const fn n_lower(self) -> u32 {
        match self {
            Self::Lyman => 1,
            Self::Balmer => 2,
            Self::Paschen => 3,
            Self::Brackett => 4,
            Self::Pfund => 5,
        }
    }

    /// Wavelength of a specific line in this series.
    ///
    /// `line` is the line number (1 = first line, α; 2 = β; etc.)
    pub fn wavelength(self, line: u32) -> Result<f64> {
        if line == 0 {
            return Err(KimiyaError::InvalidInput(
                "line number must be at least 1".into(),
            ));
        }
        rydberg_wavelength(1, self.n_lower(), self.n_lower() + line)
    }
}

// ── Wavenumber conversions ───────────────────────────────────────────

/// Convert wavelength (m) to wavenumber (cm⁻¹): ν̃ = 1/λ (in cm)
///
/// # Errors
///
/// Returns error if wavelength is not positive.
#[inline]
pub fn wavenumber_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if wavelength_m <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "wavelength must be positive".into(),
        ));
    }
    Ok(0.01 / wavelength_m) // 1/(λ in cm) = 0.01/λ_m
}

/// Convert wavenumber (cm⁻¹) to wavelength (m): λ = 1/ν̃ (in cm) → m
///
/// # Errors
///
/// Returns error if wavenumber is not positive.
#[inline]
pub fn wavelength_from_wavenumber(wavenumber_cm_inv: f64) -> Result<f64> {
    if wavenumber_cm_inv <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "wavenumber must be positive".into(),
        ));
    }
    Ok(0.01 / wavenumber_cm_inv)
}

// ── Harmonic oscillator ──────────────────────────────────────────────

/// Reduced mass of a diatomic: μ = m₁·m₂ / (m₁ + m₂)
///
/// Masses in any consistent unit (kg, u, etc.)
///
/// # Errors
///
/// Returns error if either mass is not positive.
#[inline]
pub fn reduced_mass(m1: f64, m2: f64) -> Result<f64> {
    if m1 <= 0.0 || m2 <= 0.0 {
        return Err(KimiyaError::InvalidInput("masses must be positive".into()));
    }
    Ok(m1 * m2 / (m1 + m2))
}

/// Quantum harmonic oscillator energy: E_n = (n + ½)ℏω
///
/// - `n`: vibrational quantum number (0, 1, 2, ...)
/// - `angular_frequency`: ω (rad/s) = √(k/μ)
///
/// Returns energy in Joules.
#[must_use]
#[inline]
pub fn harmonic_oscillator_energy(n: u32, angular_frequency: f64) -> f64 {
    let hbar = PLANCK / (2.0 * std::f64::consts::PI);
    (n as f64 + 0.5) * hbar * angular_frequency
}

/// Vibrational frequency from force constant and reduced mass:
/// ω = √(k/μ)
///
/// - `force_constant`: k (N/m)
/// - `reduced_mass_kg`: μ (kg)
///
/// Returns angular frequency in rad/s.
///
/// # Errors
///
/// Returns error if parameters are not positive.
#[inline]
pub fn vibrational_frequency(force_constant: f64, reduced_mass_kg: f64) -> Result<f64> {
    if force_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "force constant must be positive".into(),
        ));
    }
    if reduced_mass_kg <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "reduced mass must be positive".into(),
        ));
    }
    Ok((force_constant / reduced_mass_kg).sqrt())
}

// ── Rotational spectroscopy ──────────────────────────────────────────

/// Rigid rotor rotational energy: E_J = B × J(J+1) × h
///
/// where B = h / (8π²I) is the rotational constant in Hz.
///
/// - `j`: rotational quantum number
/// - `rotational_constant_hz`: B in Hz
///
/// Returns energy in Joules.
#[must_use]
#[inline]
pub fn rigid_rotor_energy(j: u32, rotational_constant_hz: f64) -> f64 {
    PLANCK * rotational_constant_hz * j as f64 * (j as f64 + 1.0)
}

/// Rotational constant B from moment of inertia:
/// B = h / (8π²I)
///
/// - `moment_of_inertia`: I in kg·m²
///
/// Returns B in Hz.
///
/// # Errors
///
/// Returns error if moment of inertia is not positive.
#[inline]
pub fn rotational_constant(moment_of_inertia: f64) -> Result<f64> {
    if moment_of_inertia <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "moment of inertia must be positive".into(),
        ));
    }
    Ok(PLANCK / (8.0 * std::f64::consts::PI * std::f64::consts::PI * moment_of_inertia))
}

/// Moment of inertia for a diatomic: I = μ × r²
///
/// - `reduced_mass_kg`: μ in kg
/// - `bond_length_m`: r in meters
///
/// Returns I in kg·m².
#[must_use]
#[inline]
pub fn diatomic_moment_of_inertia(reduced_mass_kg: f64, bond_length_m: f64) -> f64 {
    reduced_mass_kg * bond_length_m * bond_length_m
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Beer-Lambert ─────────────────────────────────────────────────

    #[test]
    fn beer_lambert_basic() {
        // ε=100, l=1cm, c=0.01M → A=1.0
        let a = absorbance(100.0, 1.0, 0.01);
        assert!((a - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn transmittance_at_a1() {
        // A=1 → T=0.1 (10%)
        let t = transmittance(1.0);
        assert!((t - 0.1).abs() < 0.001);
    }

    #[test]
    fn transmittance_at_a0() {
        // A=0 → T=1.0 (100%)
        let t = transmittance(0.0);
        assert!((t - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn absorbance_transmittance_roundtrip() {
        let a = 0.75;
        let t = transmittance(a);
        let back = absorbance_from_transmittance(t).unwrap();
        assert!((back - a).abs() < 1e-10);
    }

    #[test]
    fn absorbance_from_transmittance_zero_is_error() {
        assert!(absorbance_from_transmittance(0.0).is_err());
    }

    #[test]
    fn absorbance_from_transmittance_over_one_is_error() {
        assert!(absorbance_from_transmittance(1.5).is_err());
    }

    #[test]
    fn concentration_from_absorbance_basic() {
        let c = concentration_from_absorbance(1.0, 100.0, 1.0).unwrap();
        assert!((c - 0.01).abs() < 1e-10);
    }

    // ── Photon conversions ───────────────────────────────────────────

    #[test]
    fn photon_energy_wavelength_roundtrip() {
        let lambda = 500e-9; // 500 nm (green light)
        let e = photon_energy_from_wavelength(lambda).unwrap();
        let back = wavelength_from_energy(e).unwrap();
        assert!((back - lambda).abs() / lambda < 1e-10);
    }

    #[test]
    fn visible_light_energy_range() {
        // Visible: 400–700 nm → ~1.77–3.10 eV
        let e_red = joules_to_ev(photon_energy_from_wavelength(700e-9).unwrap());
        let e_violet = joules_to_ev(photon_energy_from_wavelength(400e-9).unwrap());
        assert!(
            e_red > 1.7 && e_red < 1.85,
            "red should be ~1.77 eV, got {e_red}"
        );
        assert!(
            e_violet > 3.0 && e_violet < 3.2,
            "violet should be ~3.10 eV, got {e_violet}"
        );
    }

    #[test]
    fn ev_joules_roundtrip() {
        let ev = 13.6;
        let j = ev_to_joules(ev);
        let back = joules_to_ev(j);
        assert!((back - ev).abs() < 1e-10);
    }

    #[test]
    fn frequency_wavelength_roundtrip() {
        let lambda = 500e-9;
        let freq = frequency_from_wavelength(lambda).unwrap();
        let back = SPEED_OF_LIGHT / freq;
        assert!((back - lambda).abs() / lambda < 1e-10);
    }

    #[test]
    fn zero_wavelength_is_error() {
        assert!(photon_energy_from_wavelength(0.0).is_err());
        assert!(frequency_from_wavelength(0.0).is_err());
    }

    #[test]
    fn zero_energy_is_error() {
        assert!(wavelength_from_energy(0.0).is_err());
    }

    // ── Bohr model ───────────────────────────────────────────────────

    #[test]
    fn hydrogen_ground_state() {
        let e = bohr_energy_level(1, 1).unwrap();
        assert!(
            (e - (-13.6)).abs() < 0.01,
            "H ground state should be ~-13.6 eV, got {e}"
        );
    }

    #[test]
    fn hydrogen_n2() {
        let e = bohr_energy_level(1, 2).unwrap();
        assert!(
            (e - (-3.4)).abs() < 0.01,
            "H n=2 should be ~-3.4 eV, got {e}"
        );
    }

    #[test]
    fn helium_ion_ground_state() {
        // He⁺ (Z=2): E₁ = -13.6 × 4 = -54.4 eV
        let e = bohr_energy_level(2, 1).unwrap();
        assert!(
            (e - (-54.4)).abs() < 0.1,
            "He⁺ ground state should be ~-54.4 eV, got {e}"
        );
    }

    #[test]
    fn bohr_zero_n_is_error() {
        assert!(bohr_energy_level(1, 0).is_err());
    }

    #[test]
    fn bohr_zero_z_is_error() {
        assert!(bohr_energy_level(0, 1).is_err());
    }

    #[test]
    fn hydrogen_lyman_alpha() {
        // Lyman-α: n=2→1, λ ≈ 121.6 nm
        let lambda = rydberg_wavelength(1, 1, 2).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 121.6).abs() < 0.5,
            "Lyman-α should be ~121.6 nm, got {lambda_nm}"
        );
    }

    #[test]
    fn hydrogen_balmer_alpha() {
        // Balmer-α (Hα): n=3→2, λ ≈ 656.3 nm
        let lambda = rydberg_wavelength(1, 2, 3).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 656.3).abs() < 0.5,
            "Hα should be ~656.3 nm, got {lambda_nm}"
        );
    }

    #[test]
    fn rydberg_invalid_quantum_numbers() {
        assert!(rydberg_wavelength(1, 0, 2).is_err());
        assert!(rydberg_wavelength(1, 2, 2).is_err()); // n_upper must be > n_lower
        assert!(rydberg_wavelength(1, 3, 2).is_err());
    }

    #[test]
    fn transition_energy_emission_positive() {
        // Emission: higher → lower gives positive energy
        let e = bohr_transition_energy(1, 1, 3).unwrap();
        assert!(e > 0.0);
    }

    // ── Spectral series ──────────────────────────────────────────────

    #[test]
    fn lyman_alpha_from_series() {
        let lambda = SpectralSeries::Lyman.wavelength(1).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!((lambda_nm - 121.6).abs() < 0.5);
    }

    #[test]
    fn balmer_alpha_from_series() {
        let lambda = SpectralSeries::Balmer.wavelength(1).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!((lambda_nm - 656.3).abs() < 0.5);
    }

    #[test]
    fn series_zero_line_is_error() {
        assert!(SpectralSeries::Lyman.wavelength(0).is_err());
    }

    #[test]
    fn balmer_series_is_visible() {
        for line in 1..=4 {
            let lambda_nm = SpectralSeries::Balmer.wavelength(line).unwrap() * 1e9;
            assert!(
                lambda_nm > 380.0 && lambda_nm < 700.0,
                "Balmer line {line} should be visible, got {lambda_nm} nm"
            );
        }
    }

    // ── Wavenumber ───────────────────────────────────────────────────

    #[test]
    fn wavenumber_roundtrip() {
        let lambda = 500e-9; // 500 nm
        let wn = wavenumber_from_wavelength(lambda).unwrap();
        let back = wavelength_from_wavenumber(wn).unwrap();
        assert!((back - lambda).abs() / lambda < 1e-10);
    }

    #[test]
    fn wavenumber_ir_range() {
        // 10 μm → 1000 cm⁻¹ (mid-IR)
        let wn = wavenumber_from_wavelength(10e-6).unwrap();
        assert!((wn - 1000.0).abs() < 0.1);
    }

    // ── Harmonic oscillator ──────────────────────────────────────────

    #[test]
    fn reduced_mass_equal_masses() {
        let mu = reduced_mass(1.0, 1.0).unwrap();
        assert!((mu - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn reduced_mass_heavy_light() {
        // H + very heavy → μ ≈ m_H
        let mu = reduced_mass(1.0, 1e6).unwrap();
        assert!((mu - 1.0).abs() < 0.001);
    }

    #[test]
    fn harmonic_oscillator_ground_state() {
        // E₀ = ½ℏω
        let hbar = PLANCK / (2.0 * std::f64::consts::PI);
        let omega = 1e14; // typical molecular vibration
        let e0 = harmonic_oscillator_energy(0, omega);
        assert!((e0 - 0.5 * hbar * omega).abs() < 1e-30);
    }

    #[test]
    fn harmonic_oscillator_spacing() {
        let omega = 1e14;
        let e0 = harmonic_oscillator_energy(0, omega);
        let e1 = harmonic_oscillator_energy(1, omega);
        let e2 = harmonic_oscillator_energy(2, omega);
        // Equal spacing: E₁ - E₀ = E₂ - E₁
        assert!(((e1 - e0) - (e2 - e1)).abs() < 1e-30);
    }

    #[test]
    fn vibrational_frequency_basic() {
        let omega = vibrational_frequency(500.0, 1.66e-27).unwrap();
        assert!(omega > 5e14 && omega < 6e14);
    }

    // ── Rotational spectroscopy ──────────────────────────────────────

    #[test]
    fn rigid_rotor_j0_is_zero() {
        let e = rigid_rotor_energy(0, 1e10);
        assert!(e.abs() < f64::EPSILON);
    }

    #[test]
    fn rigid_rotor_increases_with_j() {
        let e1 = rigid_rotor_energy(1, 1e10);
        let e2 = rigid_rotor_energy(2, 1e10);
        assert!(e2 > e1);
    }

    #[test]
    fn rotational_constant_positive() {
        let b = rotational_constant(1e-46).unwrap();
        assert!(b > 0.0);
    }

    #[test]
    fn rotational_constant_zero_inertia_is_error() {
        assert!(rotational_constant(0.0).is_err());
    }

    #[test]
    fn diatomic_moment_of_inertia_basic() {
        let mu = 1e-26; // ~10 u
        let r = 1e-10; // 1 Å
        let i = diatomic_moment_of_inertia(mu, r);
        assert!((i - 1e-46).abs() < 1e-48);
    }
}
