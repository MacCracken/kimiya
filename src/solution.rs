use crate::error::{KimiyaError, Result};
use hisab::EPSILON_F64;

/// Molarity: M = moles / volume_liters
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if volume is not positive.
#[inline]
pub fn molarity(moles: f64, volume_liters: f64) -> Result<f64> {
    if volume_liters <= 0.0 {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    Ok(moles / volume_liters)
}

/// Molality: m = moles / solvent_kg
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if solvent mass is not positive.
#[inline]
pub fn molality(moles: f64, solvent_kg: f64) -> Result<f64> {
    if solvent_kg <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "solvent mass must be positive".into(),
        ));
    }
    Ok(moles / solvent_kg)
}

/// Dilution: M1×V1 = M2×V2 → M2 = M1×V1/V2
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if final volume is not positive.
#[inline]
pub fn dilution(m1: f64, v1: f64, v2: f64) -> Result<f64> {
    if v2 <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "final volume must be positive".into(),
        ));
    }
    Ok(m1 * v1 / v2)
}

/// pH from hydrogen ion concentration: pH = -log₁₀([H⁺])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not positive.
#[inline]
pub fn ph_from_h_concentration(h_molar: f64) -> Result<f64> {
    if h_molar <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "[H⁺] must be positive".into(),
        ));
    }
    Ok(-h_molar.log10())
}

/// pOH from hydroxide ion concentration: pOH = -log₁₀([OH⁻])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not positive.
#[inline]
pub fn poh_from_oh_concentration(oh_molar: f64) -> Result<f64> {
    if oh_molar <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "[OH⁻] must be positive".into(),
        ));
    }
    Ok(-oh_molar.log10())
}

/// pH + pOH = 14 (at 25°C)
#[must_use]
#[inline]
pub fn ph_from_poh(poh: f64) -> f64 {
    14.0 - poh
}

/// Henderson-Hasselbalch equation: pH = pKa + log₁₀(\[A⁻\]/\[HA\])
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if acid or base concentration is not positive.
#[inline]
pub fn henderson_hasselbalch(pka: f64, conjugate_base: f64, acid: f64) -> Result<f64> {
    if acid <= 0.0 || conjugate_base <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "acid and conjugate base concentrations must be positive".into(),
        ));
    }
    Ok(pka + (conjugate_base / acid).log10())
}

/// H⁺ concentration from pH: [H⁺] = 10^(-pH)
#[must_use]
#[inline]
pub fn h_concentration_from_ph(ph: f64) -> f64 {
    10.0_f64.powf(-ph)
}

/// pH of a weak acid solution.
///
/// Solves the equilibrium equation: Ka = x² / (C - x) where x = \[H⁺\].
/// Uses [`hisab::num::newton_raphson`] with the analytical derivative f'(x) = 2x + Ka.
///
/// - `ka`: acid dissociation constant
/// - `concentration`: initial acid concentration (M)
///
/// # Errors
///
/// Returns error if Ka or concentration is not positive, or if root finding fails.
pub fn weak_acid_ph(ka: f64, concentration: f64) -> Result<f64> {
    if ka <= 0.0 {
        return Err(KimiyaError::InvalidInput("Ka must be positive".into()));
    }
    if concentration <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be positive".into(),
        ));
    }
    // Solve: x² + Ka·x - Ka·C = 0, f'(x) = 2x + Ka
    let f = |x: f64| x * x + ka * x - ka * concentration;
    let df = |x: f64| 2.0 * x + ka;
    // Initial guess: x₀ = √(Ka·C) (weak acid approximation)
    let x0 = (ka * concentration).sqrt();
    let root = hisab::num::newton_raphson(f, df, x0, EPSILON_F64, 50)
        .map_err(|e| KimiyaError::ComputationError(format!("newton_raphson failed: {e}")))?;
    ph_from_h_concentration(root.abs())
}

/// pH of a weak base solution.
///
/// Solves the equilibrium equation: Kb = x² / (C - x) where x = \[OH⁻\].
/// Uses [`hisab::num::newton_raphson`] with analytical derivative.
/// Returns pH = 14 - pOH.
///
/// - `kb`: base dissociation constant
/// - `concentration`: initial base concentration (M)
///
/// # Errors
///
/// Returns error if Kb or concentration is not positive.
pub fn weak_base_ph(kb: f64, concentration: f64) -> Result<f64> {
    if kb <= 0.0 {
        return Err(KimiyaError::InvalidInput("Kb must be positive".into()));
    }
    if concentration <= 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be positive".into(),
        ));
    }
    let f = |x: f64| x * x + kb * x - kb * concentration;
    let df = |x: f64| 2.0 * x + kb;
    let x0 = (kb * concentration).sqrt();
    let root = hisab::num::newton_raphson(f, df, x0, EPSILON_F64, 50)
        .map_err(|e| KimiyaError::ComputationError(format!("newton_raphson failed: {e}")))?;
    let poh = poh_from_oh_concentration(root.abs())?;
    Ok(ph_from_poh(poh))
}

// ── Ionic strength & activity coefficients ───────────────────────────

/// Ionic strength: I = ½ Σ cᵢ·zᵢ²
///
/// - `ions`: slice of (concentration_mol_per_L, charge) pairs
#[must_use]
pub fn ionic_strength(ions: &[(f64, i32)]) -> f64 {
    0.5 * ions.iter().map(|&(c, z)| c * (z * z) as f64).sum::<f64>()
}

/// Debye-Huckel limiting law: log₁₀(γ) = -A·z²·√I
///
/// Valid for I < 0.01 M.
///
/// - `charge`: ion charge z
/// - `ionic_str`: ionic strength I (mol/L)
/// - `a_constant`: Debye-Huckel A parameter (0.509 for water at 25°C)
///
/// Returns activity coefficient γ.
///
/// # Errors
///
/// Returns error if ionic strength is negative.
#[inline]
pub fn debye_huckel_limiting(charge: i32, ionic_str: f64, a_constant: f64) -> Result<f64> {
    if ionic_str < 0.0 {
        return Err(KimiyaError::InvalidInput(
            "ionic strength must be non-negative".into(),
        ));
    }
    let log_gamma = -a_constant * (charge * charge) as f64 * ionic_str.sqrt();
    Ok(10.0_f64.powf(log_gamma))
}

/// Extended Debye-Huckel: log₁₀(γ) = -A·z²·√I / (1 + B·a·√I)
///
/// Valid for I < 0.1 M.
///
/// - `charge`: ion charge z
/// - `ionic_str`: ionic strength I (mol/L)
/// - `a_constant`: Debye-Huckel A (0.509 for water at 25°C)
/// - `b_constant`: Debye-Huckel B (0.328 × 10¹⁰ m⁻¹ for water at 25°C)
/// - `ion_size`: effective ion diameter å (Angstroms, typically 3–9)
///
/// # Errors
///
/// Returns error if ionic strength is negative.
#[inline]
pub fn debye_huckel_extended(
    charge: i32,
    ionic_str: f64,
    a_constant: f64,
    b_constant: f64,
    ion_size: f64,
) -> Result<f64> {
    if ionic_str < 0.0 {
        return Err(KimiyaError::InvalidInput(
            "ionic strength must be non-negative".into(),
        ));
    }
    let sqrt_i = ionic_str.sqrt();
    let log_gamma =
        -a_constant * (charge * charge) as f64 * sqrt_i / (1.0 + b_constant * ion_size * sqrt_i);
    Ok(10.0_f64.powf(log_gamma))
}

/// Davies equation: log₁₀(γ) = -A·z²·(√I/(1+√I) - 0.3·I)
///
/// Valid for I < 0.5 M. No ion-size parameter needed.
///
/// # Errors
///
/// Returns error if ionic strength is negative.
#[inline]
pub fn davies_activity(charge: i32, ionic_str: f64, a_constant: f64) -> Result<f64> {
    if ionic_str < 0.0 {
        return Err(KimiyaError::InvalidInput(
            "ionic strength must be non-negative".into(),
        ));
    }
    let sqrt_i = ionic_str.sqrt();
    let log_gamma =
        -a_constant * (charge * charge) as f64 * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_str);
    Ok(10.0_f64.powf(log_gamma))
}

/// Debye-Huckel A parameter for water at 25°C.
pub const DH_A_WATER_25C: f64 = 0.509;

/// Debye-Huckel B parameter for water at 25°C (in Å⁻¹·M⁻¹/²).
pub const DH_B_WATER_25C: f64 = 0.328;

// ── Conductivity ─────────────────────────────────────────────────────

/// Kohlrausch's law: Λ_m = Λ°_m - K√c
///
/// Molar conductivity at concentration c from limiting molar conductivity.
///
/// - `lambda_zero`: Λ°_m (S·cm²/mol) — limiting molar conductivity at infinite dilution
/// - `k_kohlrausch`: K — Kohlrausch constant (S·cm²·L^½ / mol^(3/2))
/// - `concentration`: c (mol/L)
///
/// # Errors
///
/// Returns error if concentration is negative.
#[inline]
pub fn kohlrausch_molar_conductivity(
    lambda_zero: f64,
    k_kohlrausch: f64,
    concentration: f64,
) -> Result<f64> {
    if concentration < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be non-negative".into(),
        ));
    }
    Ok(lambda_zero - k_kohlrausch * concentration.sqrt())
}

/// Kohlrausch's law of independent ion migration:
/// Λ°_m = ν₊·λ°₊ + ν₋·λ°₋
///
/// - `cation_conductivity`: λ°₊ (S·cm²/mol)
/// - `cation_stoich`: ν₊
/// - `anion_conductivity`: λ°₋ (S·cm²/mol)
/// - `anion_stoich`: ν₋
#[must_use]
#[inline]
pub fn limiting_molar_conductivity(
    cation_conductivity: f64,
    cation_stoich: f64,
    anion_conductivity: f64,
    anion_stoich: f64,
) -> f64 {
    cation_stoich * cation_conductivity + anion_stoich * anion_conductivity
}

/// Solution conductivity: κ = Λ_m × c / 1000
///
/// - `molar_conductivity`: Λ_m (S·cm²/mol)
/// - `concentration`: c (mol/L)
///
/// Returns conductivity κ in S/cm.
#[must_use]
#[inline]
pub fn solution_conductivity(molar_conductivity: f64, concentration: f64) -> f64 {
    molar_conductivity * concentration / 1000.0
}

/// Common limiting ionic conductivities λ° at 25°C (S·cm²/mol).
pub const LAMBDA_H_PLUS: f64 = 349.8;
pub const LAMBDA_OH_MINUS: f64 = 198.0;
pub const LAMBDA_NA_PLUS: f64 = 50.1;
pub const LAMBDA_K_PLUS: f64 = 73.5;
pub const LAMBDA_CL_MINUS: f64 = 76.3;
pub const LAMBDA_NO3_MINUS: f64 = 71.4;
pub const LAMBDA_SO4_2MINUS: f64 = 160.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_water_ph() {
        let ph = ph_from_h_concentration(1e-7).unwrap();
        assert!(
            (ph - 7.0).abs() < 0.01,
            "pure water pH should be ~7, got {ph}"
        );
    }

    #[test]
    fn strong_acid_ph() {
        let ph = ph_from_h_concentration(0.1).unwrap();
        assert!((ph - 1.0).abs() < 0.01);
    }

    #[test]
    fn ph_poh_sum_14() {
        let ph = 3.0;
        let poh = 14.0 - ph;
        assert!((ph_from_poh(poh) - ph).abs() < f64::EPSILON);
    }

    #[test]
    fn ph_h_roundtrip() {
        let ph = 4.5;
        let h = h_concentration_from_ph(ph);
        let back = ph_from_h_concentration(h).unwrap();
        assert!((back - ph).abs() < 0.001);
    }

    #[test]
    fn molarity_basic() {
        let m = molarity(0.5, 0.25).unwrap();
        assert!((m - 2.0).abs() < f64::EPSILON, "0.5 mol in 0.25 L = 2M");
    }

    #[test]
    fn dilution_basic() {
        let m2 = dilution(1.0, 0.1, 0.5).unwrap();
        assert!((m2 - 0.2).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_equal_ratio() {
        let ph = henderson_hasselbalch(4.75, 0.1, 0.1).unwrap();
        assert!((ph - 4.75).abs() < 0.001);
    }

    #[test]
    fn henderson_hasselbalch_more_base() {
        let ph = henderson_hasselbalch(4.75, 1.0, 0.1).unwrap();
        assert!(ph > 4.75);
    }

    #[test]
    fn zero_volume_molarity_is_error() {
        assert!(molarity(1.0, 0.0).is_err());
    }

    #[test]
    fn zero_concentration_ph_is_error() {
        assert!(ph_from_h_concentration(0.0).is_err());
    }

    #[test]
    fn negative_concentration_ph_is_error() {
        assert!(ph_from_h_concentration(-1.0).is_err());
    }

    #[test]
    fn molality_zero_solvent_is_error() {
        assert!(molality(1.0, 0.0).is_err());
    }

    #[test]
    fn dilution_zero_volume_is_error() {
        assert!(dilution(1.0, 0.1, 0.0).is_err());
    }

    #[test]
    fn henderson_hasselbalch_zero_acid_is_error() {
        assert!(henderson_hasselbalch(4.75, 0.1, 0.0).is_err());
    }

    #[test]
    fn henderson_hasselbalch_zero_base_is_error() {
        assert!(henderson_hasselbalch(4.75, 0.0, 0.1).is_err());
    }

    #[test]
    fn poh_zero_concentration_is_error() {
        assert!(poh_from_oh_concentration(0.0).is_err());
    }

    #[test]
    fn h_concentration_from_ph_roundtrip_strong_acid() {
        // pH 0 → [H⁺] = 1.0 M
        let h = h_concentration_from_ph(0.0);
        assert!((h - 1.0).abs() < f64::EPSILON);
    }

    // ── Weak acid/base tests (hisab-powered) ─────────────────────────

    #[test]
    fn acetic_acid_ph() {
        // 0.1 M acetic acid, Ka = 1.8e-5 → pH ≈ 2.87
        let ph = weak_acid_ph(1.8e-5, 0.1).unwrap();
        assert!(
            (ph - 2.87).abs() < 0.05,
            "0.1M acetic acid should have pH ~2.87, got {ph}"
        );
    }

    #[test]
    fn hf_weak_acid_ph() {
        // 0.1 M HF, Ka = 6.8e-4
        // x² + 6.8e-4·x - 6.8e-5 = 0 → x ≈ 0.00791 → pH ≈ 2.10
        let ph = weak_acid_ph(6.8e-4, 0.1).unwrap();
        assert!(
            (ph - 2.10).abs() < 0.05,
            "0.1M HF should have pH ~2.10, got {ph}"
        );
    }

    #[test]
    fn weak_acid_dilute() {
        // Very dilute: 1e-5 M acetic acid → pH closer to 7
        let ph = weak_acid_ph(1.8e-5, 1e-5).unwrap();
        assert!(
            ph > 5.0 && ph < 7.0,
            "very dilute weak acid pH should be 5-7, got {ph}"
        );
    }

    #[test]
    fn weak_acid_zero_ka_is_error() {
        assert!(weak_acid_ph(0.0, 0.1).is_err());
    }

    #[test]
    fn weak_acid_zero_concentration_is_error() {
        assert!(weak_acid_ph(1.8e-5, 0.0).is_err());
    }

    #[test]
    fn ammonia_weak_base_ph() {
        // 0.1 M NH₃, Kb = 1.8e-5 → pOH ≈ 2.87 → pH ≈ 11.13
        let ph = weak_base_ph(1.8e-5, 0.1).unwrap();
        assert!(
            (ph - 11.13).abs() < 0.05,
            "0.1M ammonia should have pH ~11.13, got {ph}"
        );
    }

    #[test]
    fn weak_base_zero_kb_is_error() {
        assert!(weak_base_ph(0.0, 0.1).is_err());
    }

    #[test]
    fn weak_acid_large_ka_approaches_strong() {
        // Very large Ka (essentially strong acid) → pH approaches -log₁₀(C)
        let ph = weak_acid_ph(1.0, 0.1).unwrap();
        let strong_ph = 1.0; // -log10(0.1)
        assert!(
            (ph - strong_ph).abs() < 0.2,
            "large Ka should approach strong acid pH, got {ph}"
        );
    }

    #[test]
    fn weak_acid_base_symmetry() {
        let ka = 1.8e-5;
        let c = 0.1;
        let ph_acid = weak_acid_ph(ka, c).unwrap();
        let ph_base = weak_base_ph(ka, c).unwrap();
        assert!(
            (ph_acid + ph_base - 14.0).abs() < 0.1,
            "acid pH + base pH should ≈ 14, got {} + {} = {}",
            ph_acid,
            ph_base,
            ph_acid + ph_base
        );
    }

    // ── Activity coefficients ────────────────────────────────────────

    #[test]
    fn ionic_strength_nacl() {
        // 0.1 M NaCl: I = 0.5 × (0.1×1² + 0.1×1²) = 0.1
        let i = ionic_strength(&[(0.1, 1), (0.1, -1)]);
        assert!((i - 0.1).abs() < f64::EPSILON);
    }

    #[test]
    fn ionic_strength_cacl2() {
        // 0.1 M CaCl₂: I = 0.5 × (0.1×4 + 0.2×1) = 0.3
        let i = ionic_strength(&[(0.1, 2), (0.2, -1)]);
        assert!((i - 0.3).abs() < f64::EPSILON);
    }

    #[test]
    fn debye_huckel_limiting_monovalent() {
        // I=0.01, z=1, A=0.509 → log(γ) = -0.509×1×0.1 = -0.0509 → γ ≈ 0.889
        let gamma = debye_huckel_limiting(1, 0.01, DH_A_WATER_25C).unwrap();
        assert!(
            (gamma - 0.889).abs() < 0.01,
            "γ for z=1 at I=0.01 should be ~0.889, got {gamma}"
        );
    }

    #[test]
    fn debye_huckel_limiting_divalent() {
        // z=2 at same I should give much lower γ
        let gamma_1 = debye_huckel_limiting(1, 0.01, DH_A_WATER_25C).unwrap();
        let gamma_2 = debye_huckel_limiting(2, 0.01, DH_A_WATER_25C).unwrap();
        assert!(gamma_2 < gamma_1, "divalent should have lower γ");
    }

    #[test]
    fn debye_huckel_zero_ionic_strength() {
        // At I=0, γ = 1 (ideal)
        let gamma = debye_huckel_limiting(1, 0.0, DH_A_WATER_25C).unwrap();
        assert!((gamma - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn extended_dh_larger_than_limiting() {
        // Extended DH gives larger γ (closer to 1) than limiting law at same I
        let gamma_lim = debye_huckel_limiting(1, 0.05, DH_A_WATER_25C).unwrap();
        let gamma_ext =
            debye_huckel_extended(1, 0.05, DH_A_WATER_25C, DH_B_WATER_25C, 3.0).unwrap();
        assert!(
            gamma_ext > gamma_lim,
            "extended DH should give larger γ: {gamma_ext} vs {gamma_lim}"
        );
    }

    #[test]
    fn davies_between_0_and_1() {
        let gamma = davies_activity(1, 0.1, DH_A_WATER_25C).unwrap();
        assert!(gamma > 0.0 && gamma < 1.0, "γ should be 0-1, got {gamma}");
    }

    #[test]
    fn negative_ionic_strength_is_error() {
        assert!(debye_huckel_limiting(1, -0.1, 0.509).is_err());
        assert!(debye_huckel_extended(1, -0.1, 0.509, 0.328, 3.0).is_err());
        assert!(davies_activity(1, -0.1, 0.509).is_err());
    }

    // ── Conductivity ─────────────────────────────────────────────────

    #[test]
    fn limiting_molar_conductivity_nacl() {
        // NaCl: Λ° = λ°(Na⁺) + λ°(Cl⁻) = 50.1 + 76.3 = 126.4
        let lambda = limiting_molar_conductivity(LAMBDA_NA_PLUS, 1.0, LAMBDA_CL_MINUS, 1.0);
        assert!((lambda - 126.4).abs() < 0.1);
    }

    #[test]
    fn limiting_molar_conductivity_hcl() {
        // HCl: H⁺ has very high mobility
        let lambda = limiting_molar_conductivity(LAMBDA_H_PLUS, 1.0, LAMBDA_CL_MINUS, 1.0);
        assert!(
            lambda > 400.0,
            "HCl should have high conductivity, got {lambda}"
        );
    }

    #[test]
    fn kohlrausch_decreases_with_concentration() {
        let lm1 = kohlrausch_molar_conductivity(126.4, 10.0, 0.01).unwrap();
        let lm2 = kohlrausch_molar_conductivity(126.4, 10.0, 0.1).unwrap();
        assert!(lm2 < lm1, "higher c should give lower Λ_m");
    }

    #[test]
    fn solution_conductivity_basic() {
        // κ = Λ_m × c / 1000
        let kappa = solution_conductivity(126.4, 0.1);
        assert!((kappa - 0.01264).abs() < 0.001);
    }
}
