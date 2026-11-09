//! Spectroscopy — Beer-Lambert law, Bohr model, photon energy conversions,
//! vibrational/rotational levels, isotope patterns and peak detection.

use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_positive};

// ── Physical constants ───────────────────────────────────────────────

/// Planck's constant h (J·s), exact (SI 2019).
pub const PLANCK: f64 = 6.62607015e-34;

/// Speed of light in vacuum c (m/s), exact.
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Rydberg constant R∞ for infinite nuclear mass (m⁻¹), CODATA 2018.
pub const RYDBERG: f64 = 1.097_373_156_816_0e7;

/// Rydberg energy R∞·h·c (eV), CODATA 2018.
///
/// This is the Bohr-model ground-state binding energy for an infinitely
/// heavy nucleus, not hydrogen's ionization energy (13.598 eV, lower because
/// of the reduced-mass correction).
pub const RYDBERG_EV: f64 = 13.605693122994;

/// Electron-to-proton mass ratio mₑ/mₚ (dimensionless), CODATA 2018.
const ELECTRON_PROTON_MASS_RATIO: f64 = 5.446_170_214_87e-4;

/// Rydberg constant for ¹H with the reduced-mass correction,
/// R_H = R∞ / (1 + mₑ/mₚ) ≈ 10 967 758.34 m⁻¹.
const RYDBERG_HYDROGEN: f64 = RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO);

/// Peaks below this fraction of the current base peak are discarded after
/// each convolution step in [`isotope_pattern`].
const ISOTOPE_PRUNE_REL: f64 = 1e-12;

/// Peaks below this fraction of the base peak are omitted from the
/// [`isotope_pattern`] result.
const ISOTOPE_REPORT_REL: f64 = 0.001;

/// Peaks closer than this (u) are merged in [`isotope_pattern`].
const ISOTOPE_BIN_WIDTH_U: f64 = 0.5;

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
/// Returns error if transmittance is not in (0, 1] (NaN is rejected).
#[inline]
pub fn absorbance_from_transmittance(transmittance: f64) -> Result<f64> {
    if !is_positive(transmittance) || transmittance > 1.0 {
        return Err(KimiyaError::InvalidInput(
            "transmittance must be in (0, 1]".into(),
        ));
    }
    Ok(-transmittance.log10())
}

/// Concentration from absorbance: c = A / (ε × l)
///
/// - `absorbance`: A (dimensionless)
/// - `molar_absorptivity`: ε in L/(mol·cm)
/// - `path_length`: l in cm
///
/// Returns concentration in mol/L.
///
/// # Errors
///
/// Returns error if molar absorptivity or path length is not positive and finite.
#[inline]
pub fn concentration_from_absorbance(
    absorbance: f64,
    molar_absorptivity: f64,
    path_length: f64,
) -> Result<f64> {
    if !is_positive(molar_absorptivity) {
        return Err(KimiyaError::InvalidInput(
            "molar absorptivity must be positive".into(),
        ));
    }
    if !is_positive(path_length) {
        return Err(KimiyaError::InvalidInput(
            "path length must be positive".into(),
        ));
    }
    Ok(absorbance / (molar_absorptivity * path_length))
}

// ── Photon energy/wavelength/frequency ───────────────────────────────

/// Photon energy from frequency: E = hν
///
/// `frequency_hz` in Hz; returns energy in J.
#[must_use]
#[inline]
pub fn photon_energy_from_frequency(frequency_hz: f64) -> f64 {
    PLANCK * frequency_hz
}

/// Photon energy from wavelength: E = hc/λ
///
/// `wavelength_m` in m; returns energy in J.
///
/// # Errors
///
/// Returns error if wavelength is not positive and finite.
#[inline]
pub fn photon_energy_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if !is_positive(wavelength_m) {
        return Err(KimiyaError::InvalidInput(
            "wavelength must be positive".into(),
        ));
    }
    Ok(PLANCK * SPEED_OF_LIGHT / wavelength_m)
}

/// Wavelength from photon energy: λ = hc/E
///
/// `energy_j` in J; returns wavelength in m.
///
/// # Errors
///
/// Returns error if energy is not positive and finite.
#[inline]
pub fn wavelength_from_energy(energy_j: f64) -> Result<f64> {
    if !is_positive(energy_j) {
        return Err(KimiyaError::InvalidInput("energy must be positive".into()));
    }
    Ok(PLANCK * SPEED_OF_LIGHT / energy_j)
}

/// Frequency from wavelength: ν = c/λ
///
/// `wavelength_m` in m; returns frequency in Hz.
///
/// # Errors
///
/// Returns error if wavelength is not positive and finite.
#[inline]
pub fn frequency_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if !is_positive(wavelength_m) {
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
    energy_j / crate::electrochemistry::ELEMENTARY_CHARGE
}

/// Convert energy in electron-volts to Joules.
#[must_use]
#[inline]
pub fn ev_to_joules(energy_ev: f64) -> f64 {
    energy_ev * crate::electrochemistry::ELEMENTARY_CHARGE
}

// ── Bohr model ───────────────────────────────────────────────────────

/// Bohr model energy level for hydrogen-like atoms: E_n = −R∞hc × Z² / n² (eV)
///
/// - `z`: nuclear charge (1 for hydrogen)
/// - `n`: principal quantum number (≥ 1)
///
/// Uses [`RYDBERG_EV`] (13.6057 eV), i.e. assumes an infinitely heavy nucleus.
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
/// Returns positive energy in eV (the photon energy released), Bohr model
/// with infinite nuclear mass (see [`bohr_energy_level`]).
///
/// # Errors
///
/// Returns error if Z, n_lower or n_upper is zero, or n_upper ≤ n_lower.
#[inline]
pub fn bohr_transition_energy(z: u8, n_lower: u32, n_upper: u32) -> Result<f64> {
    let e_lower = bohr_energy_level(z, n_lower)?;
    let e_upper = bohr_energy_level(z, n_upper)?;
    if n_upper <= n_lower {
        return Err(KimiyaError::InvalidInput(
            "n_upper must be greater than n_lower".into(),
        ));
    }
    Ok(e_upper - e_lower) // upper is less negative, so this is positive for emission
}

/// Wavelength of a photon from a hydrogen-like transition (Rydberg formula).
///
/// 1/λ = R∞ × Z² × (1/n₁² - 1/n₂²) where n₂ > n₁
///
/// Uses [`RYDBERG`] (R∞), i.e. assumes an infinitely heavy nucleus; for ¹H
/// this gives wavelengths about 0.05 % short (Hα 656.112 nm vs 656.470 nm).
/// [`SpectralSeries::wavelength`] applies the hydrogen reduced-mass
/// correction.
///
/// Returns wavelength in meters (vacuum).
///
/// # Errors
///
/// Returns error if quantum numbers are invalid or n_upper ≤ n_lower.
pub fn rydberg_wavelength(z: u8, n_lower: u32, n_upper: u32) -> Result<f64> {
    rydberg_formula(RYDBERG, z, n_lower, n_upper)
}

/// Rydberg formula with an explicit Rydberg constant `rydberg` (m⁻¹).
#[inline]
fn rydberg_formula(rydberg: f64, z: u8, n_lower: u32, n_upper: u32) -> Result<f64> {
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
    let inv_lambda = rydberg
        * (z as f64)
        * (z as f64)
        * (1.0 / (n_lower as f64 * n_lower as f64) - 1.0 / (n_upper as f64 * n_upper as f64));
    Ok(1.0 / inv_lambda)
}

/// Named hydrogen (¹H) spectral series.
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

    /// Vacuum wavelength (m) of a specific line in this ¹H series.
    ///
    /// `line` is the line number (1 = first line, α; 2 = β; etc.), so
    /// n_upper = n_lower + line.
    ///
    /// Uses the reduced-mass-corrected hydrogen Rydberg constant
    /// R_H = R∞/(1 + mₑ/mₚ) = 10 967 758.34 m⁻¹ (Bohr model, no fine
    /// structure): Lyα 121.568 nm, Hα 656.470 nm.
    ///
    /// # Errors
    ///
    /// Returns error if `line` is zero or n_lower + line overflows `u32`.
    pub fn wavelength(self, line: u32) -> Result<f64> {
        if line == 0 {
            return Err(KimiyaError::InvalidInput(
                "line number must be at least 1".into(),
            ));
        }
        let n_upper = self
            .n_lower()
            .checked_add(line)
            .ok_or_else(|| KimiyaError::InvalidInput("line number too large".into()))?;
        rydberg_formula(RYDBERG_HYDROGEN, 1, self.n_lower(), n_upper)
    }
}

// ── Wavenumber conversions ───────────────────────────────────────────

/// Convert wavelength (m) to wavenumber (cm⁻¹): ν̃ = 1/λ (in cm)
///
/// # Errors
///
/// Returns error if wavelength is not positive and finite.
#[inline]
pub fn wavenumber_from_wavelength(wavelength_m: f64) -> Result<f64> {
    if !is_positive(wavelength_m) {
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
/// Returns error if wavenumber is not positive and finite.
#[inline]
pub fn wavelength_from_wavenumber(wavenumber_cm_inv: f64) -> Result<f64> {
    if !is_positive(wavenumber_cm_inv) {
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
/// Returns error if either mass is not positive and finite.
#[inline]
pub fn reduced_mass(m1: f64, m2: f64) -> Result<f64> {
    if !is_positive(m1) || !is_positive(m2) {
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
/// Returns error if parameters are not positive and finite.
#[inline]
pub fn vibrational_frequency(force_constant: f64, reduced_mass_kg: f64) -> Result<f64> {
    if !is_positive(force_constant) {
        return Err(KimiyaError::InvalidInput(
            "force constant must be positive".into(),
        ));
    }
    if !is_positive(reduced_mass_kg) {
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
/// Returns error if moment of inertia is not positive and finite.
#[inline]
pub fn rotational_constant(moment_of_inertia: f64) -> Result<f64> {
    if !is_positive(moment_of_inertia) {
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

// ── Mass spectrometry isotope patterns ───────────────────────────────

/// Natural abundance data for common isotopes.
///
/// Returns (isotope mass in u, abundance as a mole fraction) pairs for a
/// given element; abundances are IUPAC 2009 representative isotopic
/// compositions. Only elements with significant isotope patterns are
/// included; all others return an empty slice.
#[must_use]
pub fn isotope_abundances(atomic_number: u8) -> &'static [(f64, f64)] {
    match atomic_number {
        1 => &[(1.00783, 0.999_885), (2.01410, 0.000_115)], // H, D (IUPAC 2009)
        6 => &[(12.0, 0.9893), (13.00335, 0.0107)],         // 12C, 13C
        7 => &[(14.00307, 0.99636), (15.00011, 0.00364)],   // 14N, 15N
        8 => &[
            (15.99491, 0.99757),
            (16.99913, 0.00038),
            (17.99916, 0.00205),
        ], // 16O, 17O, 18O
        16 => &[(31.97207, 0.9499), (32.97146, 0.0075), (33.96787, 0.0425)], // 32S, 33S, 34S
        17 => &[(34.96885, 0.7576), (36.96590, 0.2424)],    // 35Cl, 37Cl
        35 => &[(78.91834, 0.5069), (80.91629, 0.4931)],    // 79Br, 81Br
        _ => &[],
    }
}

/// Compute the isotope pattern (mass spectrum envelope) of a molecule.
///
/// Takes a list of (atomic_number, count) pairs and returns
/// (mass in u, relative intensity) peaks sorted by mass.
///
/// The isotope distribution is convolved atom by atom; peaks closer than
/// 0.5 u are merged (intensity-weighted mean mass) and intermediate peaks
/// below 10⁻¹² of the base peak are discarded, which keeps large molecules
/// fast and free of underflow. Intensities are normalized to the base peak
/// (most abundant = 1.0), and peaks below 0.1 % of it are omitted.
///
/// Elements without data in [`isotope_abundances`] contribute a single peak
/// at their *standard atomic weight* (not their monoisotopic mass) — a known
/// limitation: their own isotopic structure is ignored.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if the molecule is empty, or
/// [`KimiyaError::InvalidElement`] for an unknown atomic number (0 or > 118).
pub fn isotope_pattern(atoms: &[(u8, u32)]) -> Result<Vec<(f64, f64)>> {
    if atoms.is_empty() {
        return Err(KimiyaError::InvalidInput("need at least one atom".into()));
    }

    // Start with a single peak at mass=0, intensity=1
    let mut pattern: Vec<(f64, f64)> = vec![(0.0, 1.0)];
    let mut scratch: Vec<(f64, f64)> = Vec::new();

    for &(z, count) in atoms {
        let abundances = isotope_abundances(z);
        if abundances.is_empty() {
            // No isotope data — shift every peak by count × standard atomic weight.
            let elem = crate::element::lookup_by_number(z)
                .ok_or_else(|| KimiyaError::InvalidElement(format!("unknown atomic number {z}")))?;
            let shift = f64::from(count) * elem.atomic_mass;
            for peak in &mut pattern {
                peak.0 += shift;
            }
            continue;
        }

        // Convolve pattern with this element's isotope distribution, `count` times
        for _ in 0..count {
            scratch.clear();
            for &(m1, i1) in &pattern {
                for &(m2, i2) in abundances {
                    scratch.push((m1 + m2, i1 * i2));
                }
            }
            scratch.sort_by(|a, b| a.0.total_cmp(&b.0));

            // Bin peaks within 0.5 u
            pattern.clear();
            for &(m, i) in &scratch {
                if let Some(last) = pattern.last_mut()
                    && (m - last.0).abs() < ISOTOPE_BIN_WIDTH_U
                {
                    let total = last.1 + i;
                    if total > 0.0 {
                        last.0 = (last.0 * last.1 + m * i) / total;
                    }
                    last.1 = total;
                    continue;
                }
                pattern.push((m, i));
            }

            // Renormalize to the base peak and prune the negligible tail, so
            // intensities never underflow and the peak count stays bounded.
            let max_intensity = pattern.iter().fold(0.0_f64, |acc, p| acc.max(p.1));
            if max_intensity > 0.0 {
                for peak in &mut pattern {
                    peak.1 /= max_intensity;
                }
            }
            pattern.retain(|p| p.1 >= ISOTOPE_PRUNE_REL);
        }
    }

    // Filter out very small peaks (intensities are already relative to the base peak)
    pattern.retain(|(_, i)| *i > ISOTOPE_REPORT_REL);

    Ok(pattern)
}

// ── Fourier transform spectroscopy ───────────────────────────────────

/// Apply FFT to a real-valued spectrum (e.g., interferogram → spectrum).
///
/// Input: real-valued signal of length N (must be power of 2).
/// Returns: magnitude spectrum (N/2 + 1 points).
///
/// # Errors
///
/// Returns error if input length is not a power of 2.
pub fn fft_spectrum(signal: &[f64]) -> Result<Vec<f64>> {
    let n = signal.len();
    if n == 0 || (n & (n - 1)) != 0 {
        return Err(KimiyaError::InvalidInput(
            "signal length must be a non-zero power of 2".into(),
        ));
    }
    let mut data: Vec<hisab::Complex> = signal
        .iter()
        .map(|&x| hisab::Complex::from_real(x))
        .collect();
    hisab::num::fft(&mut data)
        .map_err(|e| KimiyaError::ComputationError(format!("FFT failed: {e}")))?;

    // Return magnitude spectrum (first half + DC)
    Ok(data.iter().take(n / 2 + 1).map(|c| c.abs()).collect())
}

/// Simple peak detection: find strict local maxima above a threshold.
///
/// - `data`: spectrum values
/// - `threshold`: minimum peak height as a fraction (0–1) of the data range;
///   a peak must exceed min + threshold × (max − min), where min and max are
///   taken over the finite values
///
/// Returns indices of detected peaks (end points are never peaks). Returns
/// an empty `Vec` if `data` has fewer than 3 points, has no finite values,
/// or `threshold` is not in \[0, 1\] (including NaN).
#[must_use]
pub fn find_peaks(data: &[f64], threshold: f64) -> Vec<usize> {
    if data.len() < 3 || !(0.0..=1.0).contains(&threshold) {
        return Vec::new();
    }
    let (min_val, max_val) = data
        .iter()
        .copied()
        .filter(|v| is_finite(*v))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        });
    if min_val > max_val {
        return Vec::new(); // no finite values
    }
    // Interpolate without forming (max − min), which can overflow.
    let abs_threshold = min_val * (1.0 - threshold) + max_val * threshold;
    data.windows(3)
        .enumerate()
        .filter(|(_, w)| w[1] > w[0] && w[1] > w[2] && w[1] > abs_threshold)
        .map(|(i, _)| i + 1)
        .collect()
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

    #[test]
    fn beer_lambert_nan_is_error() {
        assert!(absorbance_from_transmittance(f64::NAN).is_err());
        assert!(concentration_from_absorbance(1.0, f64::NAN, 1.0).is_err());
        assert!(concentration_from_absorbance(1.0, 100.0, f64::NAN).is_err());
        assert!(concentration_from_absorbance(1.0, f64::INFINITY, 1.0).is_err());
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

    #[test]
    fn photon_nan_inf_is_error() {
        assert!(photon_energy_from_wavelength(f64::NAN).is_err());
        assert!(photon_energy_from_wavelength(f64::INFINITY).is_err());
        assert!(wavelength_from_energy(f64::NAN).is_err());
        assert!(frequency_from_wavelength(f64::NAN).is_err());
        assert!(wavenumber_from_wavelength(f64::NAN).is_err());
        assert!(wavelength_from_wavenumber(f64::NAN).is_err());
        assert!(wavelength_from_wavenumber(f64::INFINITY).is_err());
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
        // Lyman-α with R∞ (infinite nuclear mass): 4/(3R∞) = 121.502 nm
        let lambda = rydberg_wavelength(1, 1, 2).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 121.502).abs() < 0.01,
            "R∞ Lyman-α should be 121.502 nm, got {lambda_nm}"
        );
    }

    #[test]
    fn hydrogen_balmer_alpha() {
        // Balmer-α with R∞ (infinite nuclear mass): 36/(5R∞) = 656.112 nm
        let lambda = rydberg_wavelength(1, 2, 3).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 656.112).abs() < 0.01,
            "R∞ Hα should be 656.112 nm, got {lambda_nm}"
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

    #[test]
    fn transition_energy_requires_upper_above_lower() {
        // Used to return a negative (or zero) "emitted" energy.
        assert!(bohr_transition_energy(1, 3, 1).is_err());
        assert!(bohr_transition_energy(1, 2, 2).is_err());
        assert!(bohr_transition_energy(1, 0, 2).is_err());
    }

    // ── Spectral series ──────────────────────────────────────────────

    #[test]
    fn lyman_alpha_from_series() {
        // ¹H with R_H = R∞/(1 + mₑ/mₚ): 4/(3R_H) = 121.568 nm (vacuum)
        let lambda = SpectralSeries::Lyman.wavelength(1).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 121.568).abs() < 0.01,
            "Lyα should be 121.568 nm, got {lambda_nm}"
        );
    }

    #[test]
    fn balmer_alpha_from_series() {
        // ¹H with R_H: 36/(5R_H) = 656.470 nm (vacuum)
        let lambda = SpectralSeries::Balmer.wavelength(1).unwrap();
        let lambda_nm = lambda * 1e9;
        assert!(
            (lambda_nm - 656.470).abs() < 0.01,
            "Hα should be 656.470 nm, got {lambda_nm}"
        );
    }

    #[test]
    fn balmer_beta_and_paschen_alpha_from_series() {
        let hb = SpectralSeries::Balmer.wavelength(2).unwrap() * 1e9;
        assert!((hb - 486.274).abs() < 0.01, "Hβ: {hb}");
        let pa = SpectralSeries::Paschen.wavelength(1).unwrap() * 1e9;
        assert!((pa - 1875.627).abs() < 0.01, "Paα: {pa}");
    }

    #[test]
    fn rydberg_hydrogen_value() {
        // R_H = 10 967 758.34 m⁻¹
        assert!((RYDBERG_HYDROGEN - 10_967_758.34).abs() < 0.01);
    }

    #[test]
    fn series_zero_line_is_error() {
        assert!(SpectralSeries::Lyman.wavelength(0).is_err());
    }

    #[test]
    fn series_line_overflow_is_error() {
        // n_lower + line used to overflow u32 (panic in debug builds).
        assert!(SpectralSeries::Pfund.wavelength(u32::MAX).is_err());
        assert!(SpectralSeries::Lyman.wavelength(u32::MAX - 1).is_ok());
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
    fn oscillator_rotor_nan_is_error() {
        assert!(reduced_mass(f64::NAN, 1.0).is_err());
        assert!(reduced_mass(1.0, f64::INFINITY).is_err());
        assert!(vibrational_frequency(f64::NAN, 1.0).is_err());
        assert!(vibrational_frequency(1.0, f64::NAN).is_err());
        assert!(rotational_constant(f64::NAN).is_err());
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
        let mu = 1e-26;
        let r = 1e-10;
        let i = diatomic_moment_of_inertia(mu, r);
        assert!((i - 1e-46).abs() < 1e-48);
    }

    // ── Isotope patterns ─────────────────────────────────────────────

    #[test]
    fn isotope_pattern_ch4() {
        // CH₄: M=16, M+1 peak from 13C
        let pattern = isotope_pattern(&[(6, 1), (1, 4)]).unwrap();
        assert!(!pattern.is_empty());
        // Base peak should be near 16
        assert!((pattern[0].0 - 16.0).abs() < 0.1);
        assert!((pattern[0].1 - 1.0).abs() < f64::EPSILON); // normalized to 1
        // M+1 / M = 0.0107/0.9893 + 4 × 0.000115/0.999885 = 0.011276 (IUPAC 2009)
        assert_eq!(pattern.len(), 2);
        assert!(
            (pattern[1].1 - 0.011_276).abs() < 1e-6,
            "M+1 = {}",
            pattern[1].1
        );
    }

    #[test]
    fn hydrogen_abundances_iupac_2009() {
        let h = isotope_abundances(1);
        assert_eq!(h, &[(1.00783, 0.999_885), (2.01410, 0.000_115)]);
        for z in [1, 6, 7, 8, 17, 35] {
            let total: f64 = isotope_abundances(z).iter().map(|(_, a)| a).sum();
            assert!((total - 1.0).abs() < 1e-9, "Z={z}: Σ = {total}");
        }
    }

    #[test]
    fn isotope_pattern_chlorine_doublet() {
        // Cl₂: M : M+2 : M+4 = 0.7576² : 2·0.7576·0.2424 : 0.2424²
        //                     = 1 : 0.6399 : 0.1024
        let pattern = isotope_pattern(&[(17, 2)]).unwrap();
        assert_eq!(pattern.len(), 3, "{pattern:?}");
        let expected = [(69.9377, 1.0), (71.9348, 0.6399), (73.9318, 0.1024)];
        for (&(m, i), (em, ei)) in pattern.iter().zip(expected) {
            assert!((m - em).abs() < 1e-3, "mass {m} vs {em}");
            assert!((i - ei).abs() < 1e-4, "intensity {i} vs {ei}");
        }
    }

    #[test]
    fn isotope_pattern_empty_is_error() {
        assert!(isotope_pattern(&[]).is_err());
    }

    #[test]
    fn isotope_pattern_large_molecules_do_not_panic() {
        // Tail intensities used to underflow to 0, merging two zero peaks
        // gave a NaN mass (0/0) and the next sort's unwrap() panicked.
        for atoms in [
            &[(1u8, 88u32)][..],
            &[(6, 168)][..],
            &[(6, 1000), (1, 2000), (8, 500), (16, 50), (35, 20)][..],
        ] {
            let pattern = isotope_pattern(atoms).unwrap();
            assert!(!pattern.is_empty());
            assert!(pattern.iter().all(|(m, i)| m.is_finite() && i.is_finite()));
            assert!(pattern.windows(2).all(|w| w[0].0 < w[1].0), "sorted");
            let base = pattern.iter().fold(0.0_f64, |acc, p| acc.max(p.1));
            assert!((base - 1.0).abs() < f64::EPSILON);
        }
    }

    #[test]
    fn isotope_pattern_c168_matches_binomial() {
        // 168 C: P(k ¹³C) ∝ C(168, k)·0.0107^k·0.9893^(168−k).
        // M+1/M = 168 × 0.0107/0.9893 = 1.8170 → M+1 is the base peak.
        let pattern = isotope_pattern(&[(6, 168)]).unwrap();
        let r = 0.0107 / 0.9893;
        let m0 = 1.0 / (168.0 * r);
        assert!((pattern[0].1 - m0).abs() < 1e-9, "M = {}", pattern[0].1);
        assert!((pattern[1].1 - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn isotope_pattern_unknown_element_is_error() {
        assert!(matches!(
            isotope_pattern(&[(0, 1)]),
            Err(KimiyaError::InvalidElement(_))
        ));
        assert!(matches!(
            isotope_pattern(&[(6, 1), (119, 2)]),
            Err(KimiyaError::InvalidElement(_))
        ));
    }

    #[test]
    fn isotope_pattern_element_without_data_uses_standard_weight() {
        // Fe has no isotope table: one peak at 2 × 55.845 (standard atomic weight).
        let fe = crate::element::lookup_by_number(26).unwrap().atomic_mass;
        let pattern = isotope_pattern(&[(26, 2)]).unwrap();
        assert_eq!(pattern.len(), 1);
        assert!((pattern[0].0 - 2.0 * fe).abs() < 1e-12);
    }

    // ── FFT spectroscopy ─────────────────────────────────────────────

    #[test]
    fn fft_spectrum_sine_wave() {
        // Pure sine wave → single peak in FFT
        let n = 64;
        let freq = 4.0; // 4 cycles in N samples
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 / n as f64).sin())
            .collect();
        let spectrum = fft_spectrum(&signal).unwrap();
        // Peak should be at index 4 (frequency bin 4)
        let peak_idx = spectrum
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        assert_eq!(peak_idx, 4, "peak should be at frequency bin 4");
    }

    #[test]
    fn fft_spectrum_non_power_of_2_is_error() {
        assert!(fft_spectrum(&[1.0; 10]).is_err());
    }

    #[test]
    fn find_peaks_basic() {
        let data = [0.0, 1.0, 3.0, 1.0, 0.5, 2.0, 0.5, 0.0];
        let peaks = find_peaks(&data, 0.1);
        assert!(peaks.contains(&2)); // peak at index 2 (value 3.0)
        assert!(peaks.contains(&5)); // peak at index 5 (value 2.0)
    }

    #[test]
    fn find_peaks_threshold_filters() {
        let data = [0.0, 1.0, 3.0, 1.0, 0.5, 0.6, 0.5, 0.0];
        let peaks = find_peaks(&data, 0.5); // threshold = 50% of max (3.0) = 1.5
        assert!(peaks.contains(&2)); // 3.0 > 1.5
        assert!(!peaks.contains(&5)); // 0.6 < 1.5
    }

    #[test]
    fn find_peaks_empty() {
        assert!(find_peaks(&[1.0, 2.0], 0.1).is_empty());
        assert!(find_peaks(&[], 0.1).is_empty());
    }

    #[test]
    fn find_peaks_all_negative() {
        // The max fold used to start at 0.0, so negative spectra found nothing.
        let data = [-5.0, -3.0, -1.0, -3.0, -5.0, -4.0, -5.0];
        assert_eq!(find_peaks(&data, 0.5), vec![2]); // −4 < −5 + 0.5·4 = −3
        assert_eq!(find_peaks(&data, 0.0), vec![2, 5]);
    }

    #[test]
    fn find_peaks_threshold_relative_to_range() {
        // Baseline offset: threshold is measured from the minimum.
        let data = [10.0, 10.2, 10.0, 11.0, 10.0];
        assert_eq!(find_peaks(&data, 0.5), vec![3]); // 10.2 < 10.5
    }

    #[test]
    fn find_peaks_invalid_threshold_is_empty() {
        let data = [0.0, 1.0, 3.0, 1.0, 0.0];
        assert!(find_peaks(&data, f64::NAN).is_empty());
        assert!(find_peaks(&data, -0.1).is_empty());
        assert!(find_peaks(&data, 1.5).is_empty());
        assert!(find_peaks(&data, f64::INFINITY).is_empty());
        assert_eq!(find_peaks(&data, 1.0), Vec::<usize>::new()); // nothing exceeds max
    }

    #[test]
    fn find_peaks_ignores_non_finite() {
        let data = [0.0, f64::NAN, 0.0, 2.0, 0.0, 1.0, 0.0];
        assert_eq!(find_peaks(&data, 0.25), vec![3, 5]);
        assert!(find_peaks(&[f64::NAN; 4], 0.5).is_empty());
    }
}
