//! Solution chemistry — concentration units, pH, weak acid/base equilibria,
//! ionic strength and activity coefficients, and electrolytic conductivity.
//!
//! Concentrations are in mol/L (M) unless stated otherwise. Acid–base
//! functions assume aqueous solution at 25 °C (Kw = [`KW_WATER_25C`],
//! pH + pOH = 14).

use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};

/// Ion product of water at 25 °C: Kw = \[H⁺\]·\[OH⁻\] = 1.0 × 10⁻¹⁴ (mol/L)².
///
/// This is the conventional value (pKw = 14) that [`ph_from_poh`] also
/// assumes; the measured value at 25 °C is 1.008 × 10⁻¹⁴ (pKw = 13.995).
pub const KW_WATER_25C: f64 = 1.0e-14;

/// Molarity: M = moles / volume_liters
///
/// - `moles`: amount of solute (mol)
/// - `volume_liters`: solution volume (L)
///
/// Returns concentration in mol/L.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if volume is not positive and finite,
/// or if `moles` is not finite.
#[inline]
pub fn molarity(moles: f64, volume_liters: f64) -> Result<f64> {
    if !is_positive(volume_liters) {
        return Err(KimiyaError::InvalidInput("volume must be positive".into()));
    }
    if !is_finite(moles) {
        return Err(KimiyaError::InvalidInput("moles must be finite".into()));
    }
    Ok(moles / volume_liters)
}

/// Molality: m = moles / solvent_kg
///
/// - `moles`: amount of solute (mol)
/// - `solvent_kg`: mass of solvent (kg)
///
/// Returns molality in mol/kg.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if solvent mass is not positive and
/// finite, or if `moles` is not finite.
#[inline]
pub fn molality(moles: f64, solvent_kg: f64) -> Result<f64> {
    if !is_positive(solvent_kg) {
        return Err(KimiyaError::InvalidInput(
            "solvent mass must be positive".into(),
        ));
    }
    if !is_finite(moles) {
        return Err(KimiyaError::InvalidInput("moles must be finite".into()));
    }
    Ok(moles / solvent_kg)
}

/// Dilution: M1×V1 = M2×V2 → M2 = M1×V1/V2
///
/// - `m1`: initial concentration (mol/L)
/// - `v1`: initial volume (any volume unit)
/// - `v2`: final volume (same unit as `v1`)
///
/// Returns the final concentration M2 in mol/L.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if final volume is not positive and
/// finite, or if `m1` or `v1` is not finite.
#[inline]
pub fn dilution(m1: f64, v1: f64, v2: f64) -> Result<f64> {
    if !is_positive(v2) {
        return Err(KimiyaError::InvalidInput(
            "final volume must be positive".into(),
        ));
    }
    if !is_finite(m1) || !is_finite(v1) {
        return Err(KimiyaError::InvalidInput(
            "initial concentration and volume must be finite".into(),
        ));
    }
    Ok(m1 * v1 / v2)
}

/// pH from hydrogen ion concentration: pH = -log₁₀(\[H⁺\])
///
/// - `h_molar`: \[H⁺\] (mol/L)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not
/// positive and finite.
#[inline]
pub fn ph_from_h_concentration(h_molar: f64) -> Result<f64> {
    if !is_positive(h_molar) {
        return Err(KimiyaError::InvalidConcentration(
            "[H⁺] must be positive".into(),
        ));
    }
    Ok(-h_molar.log10())
}

/// pOH from hydroxide ion concentration: pOH = -log₁₀(\[OH⁻\])
///
/// - `oh_molar`: \[OH⁻\] (mol/L)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is not
/// positive and finite.
#[inline]
pub fn poh_from_oh_concentration(oh_molar: f64) -> Result<f64> {
    if !is_positive(oh_molar) {
        return Err(KimiyaError::InvalidConcentration(
            "[OH⁻] must be positive".into(),
        ));
    }
    Ok(-oh_molar.log10())
}

/// pH + pOH = 14 (at 25°C)
///
/// Non-finite input propagates to the result.
#[must_use]
#[inline]
pub fn ph_from_poh(poh: f64) -> f64 {
    14.0 - poh
}

/// Henderson-Hasselbalch equation: pH = pKa + log₁₀(\[A⁻\]/\[HA\])
///
/// - `pka`: −log₁₀(Ka) of the acid
/// - `conjugate_base`: \[A⁻\] (mol/L)
/// - `acid`: \[HA\] (mol/L)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if acid or base concentration
/// is not positive and finite, and [`KimiyaError::InvalidInput`] if `pka` is
/// not finite.
#[inline]
pub fn henderson_hasselbalch(pka: f64, conjugate_base: f64, acid: f64) -> Result<f64> {
    if !is_positive(acid) || !is_positive(conjugate_base) {
        return Err(KimiyaError::InvalidConcentration(
            "acid and conjugate base concentrations must be positive".into(),
        ));
    }
    if !is_finite(pka) {
        return Err(KimiyaError::InvalidInput("pKa must be finite".into()));
    }
    Ok(pka + (conjugate_base / acid).log10())
}

/// H⁺ concentration from pH: \[H⁺\] = 10^(-pH)
///
/// Returns \[H⁺\] in mol/L. Non-finite input propagates to the result.
#[must_use]
#[inline]
pub fn h_concentration_from_ph(ph: f64) -> f64 {
    10.0_f64.powf(-ph)
}

/// Iteration cap for [`charge_balance_root`]. Newton from the upper bound
/// converges in ≤ 7 iterations for Ka ∈ \[1e-40, 1e40\], C ∈ \[1e-40, 1e20\].
const CHARGE_BALANCE_MAX_ITER: usize = 50;

/// Relative Newton step below which [`charge_balance_root`] stops. Newton
/// converges quadratically here, so the remaining error is ~(1e-12)².
const CHARGE_BALANCE_REL_TOL: f64 = 1e-12;

/// Positive root x of the monoprotic weak-electrolyte charge balance
///
/// x³ + K·x² − (K·C + Kw)·x − K·Kw = 0
///
/// where x = \[H⁺\] for an acid (K = Ka) or \[OH⁻\] for a base (K = Kb).
/// It follows from x = K·C/(K + x) + Kw/x (ionised electrolyte + water).
///
/// The cubic has exactly one positive root, is convex for x > 0 and
/// increasing above the root, and both √(K·C + Kw) and
/// (C + √(C² + 4·Kw))/2 are upper bounds of the root. Newton started from
/// the smaller bound therefore decreases monotonically onto the root.
/// The stop test is on the relative step, so it is scale-free.
#[inline]
fn charge_balance_root(k: f64, c: f64) -> Result<f64> {
    let kc_kw = k * c + KW_WATER_25C;
    let weak_bound = kc_kw.sqrt();
    let strong_bound = 0.5 * (c + (c * c + 4.0 * KW_WATER_25C).sqrt());
    let mut x = weak_bound.min(strong_bound);
    for _ in 0..CHARGE_BALANCE_MAX_ITER {
        let fx = x * x * (x + k) - k * (c * x + KW_WATER_25C) - KW_WATER_25C * x;
        if fx <= 0.0 {
            // At the root to within rounding (iterates approach from above).
            return checked_root(x);
        }
        let dfx = x * (3.0 * x + 2.0 * k) - kc_kw;
        let step = fx / dfx;
        if !step.is_finite() {
            // Only reachable when K·C or x³ overflows.
            break;
        }
        x -= step;
        if step <= CHARGE_BALANCE_REL_TOL * x {
            return checked_root(x);
        }
    }
    Err(KimiyaError::ComputationError(
        "weak electrolyte charge balance did not converge".into(),
    ))
}

/// Accept a charge-balance root only if it is a usable concentration.
#[inline]
fn checked_root(x: f64) -> Result<f64> {
    if is_positive(x) {
        Ok(x)
    } else {
        Err(KimiyaError::ComputationError(
            "weak electrolyte charge balance did not converge".into(),
        ))
    }
}

/// pH of a monoprotic weak acid solution at 25 °C.
///
/// Solves the full charge balance \[H⁺\] = \[A⁻\] + \[OH⁻\], including water
/// autoionisation, for h = \[H⁺\]:
///
/// h³ + Ka·h² − (Ka·C + Kw)·h − Ka·Kw = 0, Kw = [`KW_WATER_25C`]
///
/// with a Newton iteration that converges monotonically from an upper bound
/// of the root. Dilute solutions therefore approach pH 7 from below (e.g.
/// 10⁻⁸ M acetic acid → pH 6.98), and very large Ka reproduces the strong
/// acid limit.
///
/// - `ka`: acid dissociation constant Ka (mol/L)
/// - `concentration`: analytical acid concentration C (mol/L)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if Ka is not positive and finite,
/// [`KimiyaError::InvalidConcentration`] if the concentration is not positive
/// and finite, and [`KimiyaError::ComputationError`] if the root cannot be
/// represented (e.g. Ka·C overflows).
#[inline]
pub fn weak_acid_ph(ka: f64, concentration: f64) -> Result<f64> {
    if !is_positive(ka) {
        return Err(KimiyaError::InvalidInput("Ka must be positive".into()));
    }
    if !is_positive(concentration) {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be positive".into(),
        ));
    }
    let h = charge_balance_root(ka, concentration)?;
    ph_from_h_concentration(h)
}

/// pH of a monoprotic weak base solution at 25 °C.
///
/// Mirror of [`weak_acid_ph`]: solves the full charge balance
/// \[OH⁻\] = \[BH⁺\] + \[H⁺\] for x = \[OH⁻\],
///
/// x³ + Kb·x² − (Kb·C + Kw)·x − Kb·Kw = 0, Kw = [`KW_WATER_25C`]
///
/// and returns pH = 14 − pOH. Dilute solutions approach pH 7 from above.
///
/// - `kb`: base dissociation constant Kb (mol/L)
/// - `concentration`: analytical base concentration C (mol/L)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if Kb is not positive and finite,
/// [`KimiyaError::InvalidConcentration`] if the concentration is not positive
/// and finite, and [`KimiyaError::ComputationError`] if the root cannot be
/// represented (e.g. Kb·C overflows).
#[inline]
pub fn weak_base_ph(kb: f64, concentration: f64) -> Result<f64> {
    if !is_positive(kb) {
        return Err(KimiyaError::InvalidInput("Kb must be positive".into()));
    }
    if !is_positive(concentration) {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be positive".into(),
        ));
    }
    let oh = charge_balance_root(kb, concentration)?;
    let poh = poh_from_oh_concentration(oh)?;
    Ok(ph_from_poh(poh))
}

// ── Ionic strength & activity coefficients ───────────────────────────

/// Ionic strength: I = ½ Σ cᵢ·zᵢ²
///
/// - `ions`: slice of (concentration_mol_per_L, charge) pairs
///
/// Returns I in mol/L. The charge is squared in `f64`, so any `i32` charge is
/// accepted. Inputs are not validated: negative or non-finite concentrations
/// propagate to the result.
#[must_use]
#[inline]
pub fn ionic_strength(ions: &[(f64, i32)]) -> f64 {
    0.5 * ions
        .iter()
        .map(|&(c, z)| {
            let z = f64::from(z);
            c * z * z
        })
        .sum::<f64>()
}

/// Validate the ionic strength and Debye–Hückel A parameter shared by the
/// activity-coefficient models.
#[inline]
fn validate_ionic_strength_and_a(ionic_str: f64, a_constant: f64) -> Result<()> {
    if !is_non_negative(ionic_str) {
        return Err(KimiyaError::InvalidInput(
            "ionic strength must be non-negative".into(),
        ));
    }
    if !is_finite(a_constant) {
        return Err(KimiyaError::InvalidInput(
            "Debye-Huckel A parameter must be finite".into(),
        ));
    }
    Ok(())
}

/// Debye-Huckel limiting law: log₁₀(γ) = -A·z²·√I
///
/// Valid for I < 0.01 M.
///
/// - `charge`: ion charge z
/// - `ionic_str`: ionic strength I (mol/L)
/// - `a_constant`: Debye-Huckel A parameter in (L/mol)^½
///   ([`DH_A_WATER_25C`] = 0.509 for water at 25°C)
///
/// Returns the (dimensionless) activity coefficient γ.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if ionic strength is negative or not
/// finite, or if `a_constant` is not finite.
#[inline]
pub fn debye_huckel_limiting(charge: i32, ionic_str: f64, a_constant: f64) -> Result<f64> {
    validate_ionic_strength_and_a(ionic_str, a_constant)?;
    let z = f64::from(charge);
    let log_gamma = -a_constant * z * z * ionic_str.sqrt();
    Ok(10.0_f64.powf(log_gamma))
}

/// Extended Debye-Huckel: log₁₀(γ) = -A·z²·√I / (1 + B·å·√I)
///
/// Valid for I < 0.1 M.
///
/// - `charge`: ion charge z
/// - `ionic_str`: ionic strength I (mol/L)
/// - `a_constant`: Debye-Huckel A in (L/mol)^½ ([`DH_A_WATER_25C`] = 0.509
///   for water at 25°C)
/// - `b_constant`: Debye-Huckel B in Å⁻¹·(mol/L)^−½ ([`DH_B_WATER_25C`] =
///   0.328 for water at 25°C), so that B·å·√I is dimensionless
/// - `ion_size`: effective ion diameter å in Å (typically 3–9)
///
/// Returns the (dimensionless) activity coefficient γ.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if ionic strength is negative or not
/// finite, if `a_constant` is not finite, or if `b_constant` or `ion_size` is
/// negative or not finite (a negative B·å can make 1 + B·å·√I zero).
#[inline]
pub fn debye_huckel_extended(
    charge: i32,
    ionic_str: f64,
    a_constant: f64,
    b_constant: f64,
    ion_size: f64,
) -> Result<f64> {
    validate_ionic_strength_and_a(ionic_str, a_constant)?;
    if !is_non_negative(b_constant) {
        return Err(KimiyaError::InvalidInput(
            "Debye-Huckel B parameter must be finite and non-negative".into(),
        ));
    }
    if !is_non_negative(ion_size) {
        return Err(KimiyaError::InvalidInput(
            "ion size must be finite and non-negative".into(),
        ));
    }
    let sqrt_i = ionic_str.sqrt();
    let z = f64::from(charge);
    let log_gamma = -a_constant * z * z * sqrt_i / (1.0 + b_constant * ion_size * sqrt_i);
    Ok(10.0_f64.powf(log_gamma))
}

/// Davies equation: log₁₀(γ) = -A·z²·(√I/(1+√I) - 0.3·I)
///
/// Valid for I < 0.5 M. No ion-size parameter needed.
///
/// - `charge`: ion charge z
/// - `ionic_str`: ionic strength I (mol/L)
/// - `a_constant`: Debye-Huckel A in (L/mol)^½ ([`DH_A_WATER_25C`] for water
///   at 25°C)
///
/// Returns the (dimensionless) activity coefficient γ.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if ionic strength is negative or not
/// finite, or if `a_constant` is not finite.
#[inline]
pub fn davies_activity(charge: i32, ionic_str: f64, a_constant: f64) -> Result<f64> {
    validate_ionic_strength_and_a(ionic_str, a_constant)?;
    let sqrt_i = ionic_str.sqrt();
    let z = f64::from(charge);
    let log_gamma = -a_constant * z * z * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_str);
    Ok(10.0_f64.powf(log_gamma))
}

/// Debye-Huckel A parameter for water at 25°C, in (L/mol)^½.
pub const DH_A_WATER_25C: f64 = 0.509;

/// Debye-Huckel B parameter for water at 25°C, in Å⁻¹·(mol/L)^−½.
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
/// Returns Λ_m in S·cm²/mol.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if concentration is negative
/// or not finite, and [`KimiyaError::InvalidInput`] if `lambda_zero` or
/// `k_kohlrausch` is not finite.
#[inline]
pub fn kohlrausch_molar_conductivity(
    lambda_zero: f64,
    k_kohlrausch: f64,
    concentration: f64,
) -> Result<f64> {
    if !is_non_negative(concentration) {
        return Err(KimiyaError::InvalidConcentration(
            "concentration must be non-negative".into(),
        ));
    }
    if !is_finite(lambda_zero) || !is_finite(k_kohlrausch) {
        return Err(KimiyaError::InvalidInput(
            "conductivity parameters must be finite".into(),
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
///
/// Returns Λ°_m in S·cm²/mol. Non-finite input propagates to the result.
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
/// Returns conductivity κ in S/cm. Non-finite input propagates to the result.
#[must_use]
#[inline]
pub fn solution_conductivity(molar_conductivity: f64, concentration: f64) -> f64 {
    molar_conductivity * concentration / 1000.0
}

// Common limiting ionic conductivities λ° at infinite dilution, 25°C.

/// Limiting ionic conductivity λ° of H⁺ at 25°C (S·cm²/mol).
pub const LAMBDA_H_PLUS: f64 = 349.8;
/// Limiting ionic conductivity λ° of OH⁻ at 25°C (S·cm²/mol).
pub const LAMBDA_OH_MINUS: f64 = 198.0;
/// Limiting ionic conductivity λ° of Na⁺ at 25°C (S·cm²/mol).
pub const LAMBDA_NA_PLUS: f64 = 50.1;
/// Limiting ionic conductivity λ° of K⁺ at 25°C (S·cm²/mol).
pub const LAMBDA_K_PLUS: f64 = 73.5;
/// Limiting ionic conductivity λ° of Cl⁻ at 25°C (S·cm²/mol).
pub const LAMBDA_CL_MINUS: f64 = 76.3;
/// Limiting ionic conductivity λ° of NO₃⁻ at 25°C (S·cm²/mol).
pub const LAMBDA_NO3_MINUS: f64 = 71.4;
/// Limiting ionic conductivity λ° of SO₄²⁻ at 25°C (S·cm²/mol), per mole of
/// SO₄²⁻ ions (i.e. 2 × λ°(½SO₄²⁻) = 2 × 80.0).
pub const LAMBDA_SO4_2MINUS: f64 = 160.0;

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert |actual − expected| ≤ tol.
    fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol,
            "{what}: expected {expected}, got {actual} (|Δ| = {:e})",
            (actual - expected).abs()
        );
    }

    #[test]
    fn pure_water_ph() {
        let ph = ph_from_h_concentration(1e-7).unwrap();
        assert_close(ph, 7.0, 1e-12, "pure water pH");
    }

    #[test]
    fn strong_acid_ph() {
        let ph = ph_from_h_concentration(0.1).unwrap();
        assert_close(ph, 1.0, 1e-12, "0.1 M strong acid pH");
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
        assert_close(back, ph, 1e-12, "pH round trip");
    }

    #[test]
    fn molarity_basic() {
        let m = molarity(0.5, 0.25).unwrap();
        assert!((m - 2.0).abs() < f64::EPSILON, "0.5 mol in 0.25 L = 2M");
    }

    #[test]
    fn dilution_basic() {
        let m2 = dilution(1.0, 0.1, 0.5).unwrap();
        assert_close(m2, 0.2, 1e-15, "dilution");
    }

    #[test]
    fn henderson_hasselbalch_equal_ratio() {
        let ph = henderson_hasselbalch(4.75, 0.1, 0.1).unwrap();
        assert_close(ph, 4.75, 1e-15, "equal ratio");
    }

    #[test]
    fn henderson_hasselbalch_more_base() {
        let ph = henderson_hasselbalch(4.75, 1.0, 0.1).unwrap();
        assert_close(ph, 5.75, 1e-12, "10:1 base:acid");
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

    #[test]
    fn non_finite_inputs_are_rejected() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        assert!(molarity(1.0, nan).is_err());
        assert!(molarity(1.0, inf).is_err());
        assert!(molarity(nan, 1.0).is_err());
        assert!(molality(1.0, nan).is_err());
        assert!(molality(inf, 1.0).is_err());
        assert!(dilution(1.0, 0.1, nan).is_err());
        assert!(dilution(nan, 0.1, 1.0).is_err());
        assert!(dilution(1.0, inf, 1.0).is_err());
        assert!(ph_from_h_concentration(nan).is_err());
        assert!(ph_from_h_concentration(inf).is_err());
        assert!(poh_from_oh_concentration(nan).is_err());
        assert!(henderson_hasselbalch(nan, 0.1, 0.1).is_err());
        assert!(henderson_hasselbalch(4.75, nan, 0.1).is_err());
        assert!(henderson_hasselbalch(4.75, 0.1, inf).is_err());
        assert!(kohlrausch_molar_conductivity(126.4, 10.0, nan).is_err());
        assert!(kohlrausch_molar_conductivity(nan, 10.0, 0.1).is_err());
        assert!(kohlrausch_molar_conductivity(126.4, inf, 0.1).is_err());
    }

    // ── Weak acid/base: full charge balance incl. water ─────────────
    //
    // Reference values: positive root of h³ + K·h² − (K·C + Kw)·h − K·Kw = 0
    // with Kw = 1e-14, solved by 400-step bisection in 50-digit decimal
    // arithmetic (Python `decimal`).

    #[test]
    fn acetic_acid_ph() {
        // 0.1 M acetic acid, Ka = 1.8e-5 → pH 2.8753 (textbook "≈ 2.87")
        let ph = weak_acid_ph(1.8e-5, 0.1).unwrap();
        assert_close(ph, 2.875_277_060_319_266_5, 1e-12, "0.1 M acetic acid");
    }

    #[test]
    fn hf_weak_acid_ph() {
        // 0.1 M HF, Ka = 6.8e-4 → [H⁺] ≈ 7.913e-3 → pH ≈ 2.1016
        let ph = weak_acid_ph(6.8e-4, 0.1).unwrap();
        assert_close(ph, 2.101646894217442, 1e-12, "0.1 M HF");
    }

    #[test]
    fn weak_acid_dilute() {
        // 1e-5 M acetic acid: ~72 % dissociated
        let ph = weak_acid_ph(1.8e-5, 1e-5).unwrap();
        assert_close(ph, 5.145294329911299, 1e-12, "1e-5 M acetic acid");
    }

    #[test]
    fn weak_acid_very_dilute_stays_acidic() {
        // 1e-8 M acetic acid: water autoionisation dominates. Without the Kw
        // term the old solver returned pH 8.0 (basic!).
        let ph = weak_acid_ph(1.8e-5, 1e-8).unwrap();
        assert!(ph < 7.0, "an acid can never be basic, got {ph}");
        assert_close(ph, 6.978420204421044, 1e-12, "1e-8 M acetic acid");
    }

    #[test]
    fn weak_acid_tiny_ka_times_c_converges() {
        // Ka·C = 1e-12: the old absolute |f| < 1e-12 stop returned the
        // initial guess (pH 6.00).
        let ph = weak_acid_ph(5e-7, 2e-6).unwrap();
        assert_close(ph, 6.103103512160274, 1e-12, "Ka 5e-7, C 2e-6");
    }

    #[test]
    fn weak_acid_very_weak() {
        let ph = weak_acid_ph(1e-12, 0.01).unwrap();
        assert_close(ph, 6.849485769896728, 1e-12, "Ka 1e-12, C 0.01");
        let ph = weak_acid_ph(1e-20, 1.0).unwrap();
        assert_close(ph, 6.999999782852868, 1e-12, "Ka 1e-20, C 1");
        let ph = weak_acid_ph(1e-30, 1e-30).unwrap();
        assert_close(ph, 7.0, 1e-12, "negligible acid → neutral water");
    }

    #[test]
    fn weak_acid_large_ka_approaches_strong() {
        // Ka = 1, C = 0.1: 91.6 % dissociated
        let ph = weak_acid_ph(1.0, 0.1).unwrap();
        assert_close(ph, 1.038_066_701_166_792_4, 1e-12, "Ka 1, C 0.1");
        let ph = weak_acid_ph(1e3, 0.1).unwrap();
        assert_close(ph, 1.000043422934786, 1e-12, "Ka 1e3, C 0.1");
        let ph = weak_acid_ph(1e10, 1e-3).unwrap();
        assert_close(ph, 2.999999995657099, 1e-12, "Ka 1e10, C 1e-3");
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
    fn weak_acid_base_non_finite_is_error() {
        for (k, c) in [
            (f64::NAN, 0.1),
            (1.8e-5, f64::NAN),
            (f64::INFINITY, 0.1),
            (1.8e-5, f64::INFINITY),
        ] {
            assert!(weak_acid_ph(k, c).is_err(), "acid ({k}, {c})");
            assert!(weak_base_ph(k, c).is_err(), "base ({k}, {c})");
        }
    }

    #[test]
    fn weak_acid_overflow_is_error_not_garbage() {
        assert!(matches!(
            weak_acid_ph(1e300, 1e300),
            Err(KimiyaError::ComputationError(_))
        ));
    }

    #[test]
    fn ammonia_weak_base_ph() {
        // 0.1 M NH₃, Kb = 1.8e-5 → pOH 2.8753 → pH 11.1247
        let ph = weak_base_ph(1.8e-5, 0.1).unwrap();
        assert_close(ph, 11.124_722_939_680_733, 1e-12, "0.1 M ammonia");
    }

    #[test]
    fn weak_base_very_dilute_stays_basic() {
        let ph = weak_base_ph(1.8e-5, 1e-8).unwrap();
        assert!(ph > 7.0, "a base can never be acidic, got {ph}");
        assert_close(ph, 7.021579795578956, 1e-12, "1e-8 M base");
    }

    #[test]
    fn weak_base_tiny_kb_times_c_converges() {
        let ph = weak_base_ph(5e-7, 2e-6).unwrap();
        assert_close(ph, 7.896896487839726, 1e-12, "Kb 5e-7, C 2e-6");
    }

    #[test]
    fn weak_base_zero_kb_is_error() {
        assert!(weak_base_ph(0.0, 0.1).is_err());
    }

    #[test]
    fn weak_acid_base_symmetry() {
        // Same cubic → acid pH + base pH = 14 exactly (pKw = 14).
        for (k, c) in [(1.8e-5, 0.1), (1.8e-5, 1e-8), (1.0, 0.1), (1e-12, 0.01)] {
            let ph_acid = weak_acid_ph(k, c).unwrap();
            let ph_base = weak_base_ph(k, c).unwrap();
            assert_close(ph_acid + ph_base, 14.0, 1e-12, "acid + base pH");
        }
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
    fn large_charges_do_not_overflow() {
        // 50_000² overflows i32; the square must be taken in f64.
        let i = ionic_strength(&[(1.0, 50_000), (1.0, i32::MIN)]);
        let expected = 0.5 * (2.5e9 + (i32::MIN as f64).powi(2));
        assert_close(i, expected, expected * 1e-15, "ionic strength");
        // γ underflows to 0 for absurd charges, but nothing panics.
        assert_eq!(debye_huckel_limiting(i32::MAX, 0.01, 0.509).unwrap(), 0.0);
        assert_eq!(
            debye_huckel_extended(i32::MIN, 0.01, 0.509, 0.328, 3.0).unwrap(),
            0.0
        );
        assert_eq!(davies_activity(50_000, 0.01, 0.509).unwrap(), 0.0);
    }

    #[test]
    fn debye_huckel_limiting_monovalent() {
        // I=0.01, z=1, A=0.509 → log(γ) = -0.0509 → γ = 0.889406
        let gamma = debye_huckel_limiting(1, 0.01, DH_A_WATER_25C).unwrap();
        assert_close(gamma, 0.889_405_887_555_891_7, 1e-12, "γ(z=1, I=0.01)");
    }

    #[test]
    fn debye_huckel_limiting_divalent() {
        // z=2 at I=0.01 → log(γ) = -0.2036 → γ = 0.625749
        let gamma = debye_huckel_limiting(2, 0.01, DH_A_WATER_25C).unwrap();
        assert_close(gamma, 0.625_748_763_354_440_6, 1e-12, "γ(z=2, I=0.01)");
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
        assert_close(gamma_ext, 0.80669717443588, 1e-12, "extended DH");
    }

    #[test]
    fn extended_dh_rejects_bad_parameters() {
        let (a, b) = (DH_A_WATER_25C, DH_B_WATER_25C);
        // 1 + B·å·√I = 0 for å = −1/(0.328·0.1)
        let singular_size = -1.0 / (b * 0.1);
        assert!(debye_huckel_extended(1, 0.01, a, b, singular_size).is_err());
        assert!(debye_huckel_extended(1, 0.01, a, b, -3.0).is_err());
        assert!(debye_huckel_extended(1, 0.01, a, b, f64::NAN).is_err());
        assert!(debye_huckel_extended(1, 0.01, a, b, f64::INFINITY).is_err());
        assert!(debye_huckel_extended(1, 0.01, a, f64::NAN, 3.0).is_err());
        assert!(debye_huckel_extended(1, 0.01, a, -b, 3.0).is_err());
        assert!(debye_huckel_extended(1, 0.01, f64::NAN, b, 3.0).is_err());
        // å = 0 reduces to the limiting law
        let ext = debye_huckel_extended(1, 0.01, a, b, 0.0).unwrap();
        let lim = debye_huckel_limiting(1, 0.01, a).unwrap();
        assert_close(ext, lim, 1e-15, "å = 0");
    }

    #[test]
    fn davies_value() {
        let gamma = davies_activity(1, 0.1, DH_A_WATER_25C).unwrap();
        assert_close(gamma, 0.781_593_943_946_833_4, 1e-12, "Davies z=1 I=0.1");
    }

    #[test]
    fn negative_ionic_strength_is_error() {
        assert!(debye_huckel_limiting(1, -0.1, 0.509).is_err());
        assert!(debye_huckel_extended(1, -0.1, 0.509, 0.328, 3.0).is_err());
        assert!(davies_activity(1, -0.1, 0.509).is_err());
    }

    #[test]
    fn non_finite_ionic_strength_or_a_is_error() {
        for i in [f64::NAN, f64::INFINITY] {
            assert!(debye_huckel_limiting(1, i, 0.509).is_err());
            assert!(debye_huckel_extended(1, i, 0.509, 0.328, 3.0).is_err());
            assert!(davies_activity(1, i, 0.509).is_err());
        }
        assert!(debye_huckel_limiting(1, 0.01, f64::NAN).is_err());
        assert!(davies_activity(1, 0.01, f64::INFINITY).is_err());
    }

    // ── Conductivity ─────────────────────────────────────────────────

    #[test]
    fn limiting_molar_conductivity_nacl() {
        // NaCl: Λ° = λ°(Na⁺) + λ°(Cl⁻) = 50.1 + 76.3 = 126.4
        let lambda = limiting_molar_conductivity(LAMBDA_NA_PLUS, 1.0, LAMBDA_CL_MINUS, 1.0);
        assert_close(lambda, 126.4, 1e-12, "Λ°(NaCl)");
    }

    #[test]
    fn limiting_molar_conductivity_hcl() {
        // HCl: H⁺ has very high mobility, Λ° = 349.8 + 76.3
        let lambda = limiting_molar_conductivity(LAMBDA_H_PLUS, 1.0, LAMBDA_CL_MINUS, 1.0);
        assert_close(lambda, 426.1, 1e-12, "Λ°(HCl)");
    }

    #[test]
    fn kohlrausch_decreases_with_concentration() {
        let lm1 = kohlrausch_molar_conductivity(126.4, 10.0, 0.01).unwrap();
        let lm2 = kohlrausch_molar_conductivity(126.4, 10.0, 0.1).unwrap();
        assert!(lm2 < lm1, "higher c should give lower Λ_m");
        assert_close(lm1, 125.4, 1e-12, "Λ_m at 0.01 M");
    }

    #[test]
    fn solution_conductivity_basic() {
        // κ = Λ_m × c / 1000
        let kappa = solution_conductivity(126.4, 0.1);
        assert_close(kappa, 0.01264, 1e-15, "κ");
    }
}
