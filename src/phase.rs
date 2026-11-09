//! Phase equilibria — Raoult's law, Henry's law, colligative properties,
//! Clausius-Clapeyron, Antoine equation.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};
use serde::Serialize;

// ── Raoult's law ─────────────────────────────────────────────────────

/// Raoult's law: P = x × P°
///
/// Vapor pressure of a component in an ideal solution.
///
/// - `mole_fraction`: x (0 to 1)
/// - `pure_vapor_pressure`: P° of pure component (Pa)
///
/// Returns the partial vapor pressure in the unit of `pure_vapor_pressure`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if mole fraction is not in \[0, 1\]
/// or pressure is negative or non-finite (NaN included).
#[inline]
pub fn raoult_pressure(mole_fraction: f64, pure_vapor_pressure: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&mole_fraction) {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be between 0 and 1".into(),
        ));
    }
    if !is_non_negative(pure_vapor_pressure) {
        return Err(KimiyaError::InvalidInput(
            "vapor pressure must be non-negative".into(),
        ));
    }
    Ok(mole_fraction * pure_vapor_pressure)
}

/// Total vapor pressure over a binary ideal solution (Raoult's law).
///
/// P_total = x_A × P°_A + (1 - x_A) × P°_B
///
/// Pressures in any consistent unit (e.g. Pa); the result is in that unit.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `x_a` is not in \[0, 1\] or
/// either pure-component vapor pressure is negative or non-finite.
#[inline]
pub fn raoult_total_binary(x_a: f64, p_pure_a: f64, p_pure_b: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&x_a) {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be between 0 and 1".into(),
        ));
    }
    if !is_non_negative(p_pure_a) || !is_non_negative(p_pure_b) {
        return Err(KimiyaError::InvalidInput(
            "vapor pressure must be non-negative".into(),
        ));
    }
    Ok(x_a * p_pure_a + (1.0 - x_a) * p_pure_b)
}

// ── Henry's law ──────────────────────────────────────────────────────

/// Henry's law: P = K_H × x
///
/// Dissolved gas partial pressure from mole fraction.
///
/// - `henry_constant`: K_H in Pa (or same pressure unit as output)
/// - `mole_fraction`: x of dissolved gas (0 to 1)
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if K_H is not positive and finite,
/// or the mole fraction is not in \[0, 1\] (NaN included).
#[inline]
pub fn henry_pressure(henry_constant: f64, mole_fraction: f64) -> Result<f64> {
    if !is_positive(henry_constant) {
        return Err(KimiyaError::InvalidInput(
            "Henry's constant must be positive".into(),
        ));
    }
    if !(0.0..=1.0).contains(&mole_fraction) {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be between 0 and 1".into(),
        ));
    }
    Ok(henry_constant * mole_fraction)
}

/// Solubility (mole fraction) of a gas from Henry's law: x = P / K_H
///
/// `partial_pressure` and `henry_constant` in the same pressure unit (e.g.
/// Pa). Henry's law is a dilute-solution limit; results approaching 1 are
/// outside its range of validity.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if K_H is not positive and finite,
/// or the partial pressure is negative or non-finite.
#[inline]
pub fn henry_solubility(partial_pressure: f64, henry_constant: f64) -> Result<f64> {
    if !is_positive(henry_constant) {
        return Err(KimiyaError::InvalidInput(
            "Henry's constant must be positive".into(),
        ));
    }
    if !is_non_negative(partial_pressure) {
        return Err(KimiyaError::InvalidInput(
            "partial pressure must be non-negative".into(),
        ));
    }
    Ok(partial_pressure / henry_constant)
}

// ── Colligative properties ───────────────────────────────────────────

/// Boiling point elevation: ΔT_b = K_b × m × i
///
/// - `kb`: ebullioscopic constant (K·kg/mol)
/// - `molality`: m (mol solute / kg solvent)
/// - `van_t_hoff_factor`: i (1 for non-electrolytes, >1 for electrolytes)
///
/// Returns ΔT_b in K.
#[must_use]
#[inline]
pub fn boiling_point_elevation(kb: f64, molality: f64, van_t_hoff_factor: f64) -> f64 {
    kb * molality * van_t_hoff_factor
}

/// Freezing point depression: ΔT_f = K_f × m × i
///
/// - `kf`: cryoscopic constant (K·kg/mol)
/// - `molality`: m (mol solute / kg solvent)
/// - `van_t_hoff_factor`: i
///
/// Returns the (positive) depression ΔT_f in K.
#[must_use]
#[inline]
pub fn freezing_point_depression(kf: f64, molality: f64, van_t_hoff_factor: f64) -> f64 {
    kf * molality * van_t_hoff_factor
}

/// Osmotic pressure: π = i × M × R × T
///
/// - `van_t_hoff_factor`: i
/// - `molarity`: M (mol/L)
/// - `temperature_k`: T in Kelvin
///
/// Returns pressure in Pa.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidTemperature`] if temperature is not positive
/// and finite, or [`KimiyaError::InvalidInput`] if the van 't Hoff factor is
/// not positive and finite or the molarity is negative or non-finite.
#[inline]
pub fn osmotic_pressure(van_t_hoff_factor: f64, molarity: f64, temperature_k: f64) -> Result<f64> {
    if !is_positive(temperature_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if !is_positive(van_t_hoff_factor) {
        return Err(KimiyaError::InvalidInput(
            "van 't Hoff factor must be positive".into(),
        ));
    }
    if !is_non_negative(molarity) {
        return Err(KimiyaError::InvalidInput(
            "molarity must be non-negative".into(),
        ));
    }
    // R in L·Pa/(mol·K) = 8314.462618 (GAS_CONSTANT × 1000 for L→m³)
    Ok(van_t_hoff_factor * molarity * GAS_CONSTANT * 1000.0 * temperature_k)
}

// Common ebullioscopic constants K_b (K·kg/mol).

/// Ebullioscopic constant K_b of water (K·kg/mol).
pub const KB_WATER: f64 = 0.512;
/// Ebullioscopic constant K_b of benzene (K·kg/mol).
pub const KB_BENZENE: f64 = 2.53;
/// Ebullioscopic constant K_b of ethanol (K·kg/mol).
pub const KB_ETHANOL: f64 = 1.22;
/// Ebullioscopic constant K_b of chloroform (K·kg/mol).
pub const KB_CHLOROFORM: f64 = 3.63;

// Common cryoscopic constants K_f (K·kg/mol).

/// Cryoscopic constant K_f of water (K·kg/mol).
pub const KF_WATER: f64 = 1.86;
/// Cryoscopic constant K_f of benzene (K·kg/mol).
pub const KF_BENZENE: f64 = 5.12;
/// Cryoscopic constant K_f of cyclohexane (K·kg/mol).
pub const KF_CYCLOHEXANE: f64 = 20.0;
/// Cryoscopic constant K_f of camphor (K·kg/mol).
pub const KF_CAMPHOR: f64 = 40.0;

// ── Clausius-Clapeyron ───────────────────────────────────────────────

/// Clausius-Clapeyron equation: ln(P₂/P₁) = -(ΔH_vap/R)(1/T₂ - 1/T₁)
///
/// Calculate vapor pressure at T₂ given P₁ at T₁.
///
/// - `p1`: vapor pressure at T₁ (any consistent unit)
/// - `t1_k`: temperature T₁ (K)
/// - `t2_k`: temperature T₂ (K)
/// - `delta_h_vap_j`: enthalpy of vaporization (J/mol)
///
/// Returns P₂ in the unit of `p1`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `p1` is not positive and finite or
/// ΔH_vap is not finite, or [`KimiyaError::InvalidTemperature`] if either
/// temperature is not positive and finite.
pub fn clausius_clapeyron(p1: f64, t1_k: f64, t2_k: f64, delta_h_vap_j: f64) -> Result<f64> {
    if !is_positive(p1) {
        return Err(KimiyaError::InvalidInput(
            "reference pressure must be positive".into(),
        ));
    }
    if !is_positive(t1_k) || !is_positive(t2_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
        ));
    }
    if !is_finite(delta_h_vap_j) {
        return Err(KimiyaError::InvalidInput(
            "enthalpy of vaporization must be finite".into(),
        ));
    }
    let exponent = -(delta_h_vap_j / GAS_CONSTANT) * (1.0 / t2_k - 1.0 / t1_k);
    Ok(p1 * exponent.exp())
}

// ── Antoine equation ─────────────────────────────────────────────────

/// Antoine equation coefficients: log₁₀(P) = A - B/(C + T)
///
/// P in mmHg, T in °C (traditional Antoine convention).
#[derive(Debug, Clone, Serialize)]
pub struct AntoineCoeffs {
    /// A (dimensionless; log₁₀ of pressure in mmHg).
    pub a: f64,
    /// B (°C; log₁₀ mmHg × °C).
    pub b: f64,
    /// C (°C). The equation has a pole at T = −C.
    pub c: f64,
    /// Lower end of the temperature range the coefficients were fitted
    /// over (°C).
    pub t_min_c: f64,
    /// Upper end of the temperature range the coefficients were fitted
    /// over (°C).
    pub t_max_c: f64,
}

/// Pascals per conventional millimetre of mercury (13.5951 g/cm³ ×
/// 9.80665 m/s² × 1 mm; NIST SP 811). The torr (101 325/760 =
/// 133.322 368 Pa) differs by 1.4 × 10⁻⁷ relative.
const PA_PER_MMHG: f64 = 133.322_387_415;

/// Vapor pressure from Antoine equation.
///
/// `temperature_c` in °C; returns pressure in mmHg.
///
/// The coefficients are only fitted over
/// \[[`t_min_c`](AntoineCoeffs::t_min_c), [`t_max_c`](AntoineCoeffs::t_max_c)\];
/// the range is **not** enforced, and temperatures outside it extrapolate
/// (accuracy degrades quickly). Returns `f64::NAN` when the equation is
/// undefined: `C + T ≤ 0` (at and below the pole) or any non-finite input.
#[must_use]
#[inline]
pub fn antoine_pressure_mmhg(coeffs: &AntoineCoeffs, temperature_c: f64) -> f64 {
    // `is_positive` also rejects a NaN/±∞ `c` or `temperature_c`.
    let denom = coeffs.c + temperature_c;
    if !is_positive(denom) || !is_finite(coeffs.a) || !is_finite(coeffs.b) {
        return f64::NAN;
    }
    10.0_f64.powf(coeffs.a - coeffs.b / denom)
}

/// Vapor pressure from Antoine equation in Pa.
///
/// Same as [`antoine_pressure_mmhg`] (including the `f64::NAN` cases and the
/// unenforced validity range), converted with 1 mmHg = 133.322 387 415 Pa.
#[must_use]
#[inline]
pub fn antoine_pressure_pa(coeffs: &AntoineCoeffs, temperature_c: f64) -> f64 {
    antoine_pressure_mmhg(coeffs, temperature_c) * PA_PER_MMHG
}

/// Built-in Antoine coefficients for common solvents (P in mmHg, T in °C).
pub static ANTOINE_DATA: &[(&str, AntoineCoeffs)] = &[
    (
        "water",
        AntoineCoeffs {
            a: 8.07131,
            b: 1730.63,
            c: 233.426,
            t_min_c: 1.0,
            t_max_c: 100.0,
        },
    ),
    (
        "ethanol",
        AntoineCoeffs {
            a: 8.20417,
            b: 1642.89,
            c: 230.300,
            t_min_c: -57.0,
            t_max_c: 80.0,
        },
    ),
    (
        "methanol",
        // Lange's Handbook, 15–84 °C set (consistent with NIST WebBook,
        // Ambrose & Sprake 1970, 288.1–356.8 K). The −14 to 65 °C range
        // belongs to a different set (7.89750 / 1474.08 / 229.13).
        AntoineCoeffs {
            a: 8.08097,
            b: 1582.271,
            c: 239.726,
            t_min_c: 15.0,
            t_max_c: 84.0,
        },
    ),
    (
        "benzene",
        AntoineCoeffs {
            a: 6.90565,
            b: 1211.03,
            c: 220.790,
            t_min_c: 8.0,
            t_max_c: 80.0,
        },
    ),
    (
        "acetone",
        AntoineCoeffs {
            a: 7.02447,
            b: 1161.00,
            c: 224.000,
            t_min_c: -20.0,
            t_max_c: 77.0,
        },
    ),
    (
        "toluene",
        AntoineCoeffs {
            a: 6.95464,
            b: 1344.80,
            c: 219.480,
            t_min_c: 6.0,
            t_max_c: 137.0,
        },
    ),
    (
        "diethyl_ether",
        AntoineCoeffs {
            a: 6.92032,
            b: 1064.07,
            c: 228.800,
            t_min_c: -60.0,
            t_max_c: 35.0,
        },
    ),
    (
        "chloroform",
        AntoineCoeffs {
            a: 6.95465,
            b: 1170.97,
            c: 226.232,
            t_min_c: -10.0,
            t_max_c: 60.0,
        },
    ),
];

/// Look up Antoine coefficients by name.
#[must_use]
#[inline]
pub fn lookup_antoine(name: &str) -> Option<&'static AntoineCoeffs> {
    ANTOINE_DATA
        .iter()
        .find(|(n, _)| *n == name)
        .map(|(_, c)| c)
}

// ── Gibbs phase rule ─────────────────────────────────────────────────

/// Gibbs phase rule: F = C - P + 2
///
/// - `components`: number of independent chemical components
/// - `phases`: number of phases present
///
/// Returns degrees of freedom.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `phases` is 0, if
/// phases > components + 2 (F would be negative), or if F does not fit in a
/// `u32` (components = `u32::MAX` with one phase).
#[inline]
pub fn phase_rule(components: u32, phases: u32) -> Result<u32> {
    if phases == 0 {
        return Err(KimiyaError::InvalidInput(
            "at least one phase must be present".into(),
        ));
    }
    // i64 cannot overflow here: |components − phases + 2| ≤ 2³² + 1.
    let f = i64::from(components) - i64::from(phases) + 2;
    if f < 0 {
        return Err(KimiyaError::InvalidInput(
            "too many phases for the number of components".into(),
        ));
    }
    u32::try_from(f)
        .map_err(|_| KimiyaError::InvalidInput("degrees of freedom exceed u32::MAX".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Raoult's law ─────────────────────────────────────────────────

    #[test]
    fn raoult_basic() {
        let p = raoult_pressure(0.8, 101325.0).unwrap();
        assert!((p - 81060.0).abs() < 1.0);
    }

    #[test]
    fn raoult_pure_component() {
        let p = raoult_pressure(1.0, 101325.0).unwrap();
        assert!((p - 101325.0).abs() < f64::EPSILON);
    }

    #[test]
    fn raoult_invalid_fraction() {
        assert!(raoult_pressure(1.5, 101325.0).is_err());
        assert!(raoult_pressure(-0.1, 101325.0).is_err());
    }

    #[test]
    fn raoult_binary_total() {
        // x_A=0.6, P°_A=100, P°_B=50 → P = 60 + 20 = 80
        let p = raoult_total_binary(0.6, 100.0, 50.0).unwrap();
        assert!((p - 80.0).abs() < 1e-12);
    }

    #[test]
    fn raoult_rejects_bad_pressures_and_nan() {
        assert!(raoult_pressure(0.5, f64::NAN).is_err());
        assert!(raoult_pressure(0.5, f64::INFINITY).is_err());
        assert!(raoult_pressure(f64::NAN, 100.0).is_err());
        assert!(raoult_total_binary(0.5, -1.0, 50.0).is_err());
        assert!(raoult_total_binary(0.5, 100.0, -1.0).is_err());
        assert!(raoult_total_binary(0.5, f64::NAN, 50.0).is_err());
        assert!(raoult_total_binary(0.5, 100.0, f64::INFINITY).is_err());
        assert!(raoult_total_binary(f64::NAN, 100.0, 50.0).is_err());
    }

    // ── Henry's law ──────────────────────────────────────────────────

    #[test]
    fn henry_basic() {
        // O₂ in water: K_H ≈ 4.3e9 Pa, at P_O2 = 21278 Pa → x ≈ 4.95e-6
        let x = henry_solubility(21278.0, 4.3e9).unwrap();
        assert!(x > 4e-6 && x < 6e-6);
    }

    #[test]
    fn henry_roundtrip() {
        let kh = 1e8;
        let x = 0.001;
        let p = henry_pressure(kh, x).unwrap();
        let back = henry_solubility(p, kh).unwrap();
        assert!((back - x).abs() < 1e-15);
    }

    #[test]
    fn henry_rejects_invalid_inputs() {
        // Negative partial pressure, x > 1, NaN P and K_H used to pass.
        assert!(henry_solubility(-1.0, 4.3e9).is_err());
        assert!(henry_solubility(f64::NAN, 4.3e9).is_err());
        assert!(henry_solubility(21278.0, f64::NAN).is_err());
        assert!(henry_solubility(21278.0, f64::INFINITY).is_err());
        assert!(henry_pressure(1e8, 1.5).is_err());
        assert!(henry_pressure(1e8, f64::NAN).is_err());
        assert!(henry_pressure(f64::NAN, 0.001).is_err());
        assert!(henry_pressure(1e8, 1.0).is_ok());
        assert!(henry_solubility(0.0, 4.3e9).is_ok());
    }

    // ── Colligative properties ───────────────────────────────────────

    #[test]
    fn bp_elevation_water() {
        // 1 molal NaCl (i=2) in water: ΔT = 0.512 × 1.0 × 2 = 1.024 K
        let dt = boiling_point_elevation(KB_WATER, 1.0, 2.0);
        assert!((dt - 1.024).abs() < 0.001);
    }

    #[test]
    fn fp_depression_water() {
        // 1 molal glucose (i=1): ΔT = 1.86 × 1.0 × 1 = 1.86 K
        let dt = freezing_point_depression(KF_WATER, 1.0, 1.0);
        assert!((dt - 1.86).abs() < f64::EPSILON);
    }

    #[test]
    fn osmotic_pressure_basic() {
        // 0.1 M NaCl (i=2) at 298K
        let pi = osmotic_pressure(2.0, 0.1, 298.15).unwrap();
        // π = 2 × 0.1 mol/L × 8314.462618 L·Pa/(mol·K) × 298.15 K = 495 791.4 Pa
        assert!((pi - 495_791.406).abs() < 0.01, "got {pi}");
    }

    #[test]
    fn osmotic_pressure_zero_temp_is_error() {
        assert!(osmotic_pressure(1.0, 0.1, 0.0).is_err());
    }

    #[test]
    fn osmotic_pressure_rejects_invalid_inputs() {
        assert!(osmotic_pressure(1.0, -0.1, 298.15).is_err());
        assert!(osmotic_pressure(1.0, f64::NAN, 298.15).is_err());
        assert!(osmotic_pressure(1.0, 0.1, f64::NAN).is_err());
        assert!(osmotic_pressure(f64::NAN, 0.1, 298.15).is_err());
        assert!(osmotic_pressure(0.0, 0.1, 298.15).is_err());
        assert!(osmotic_pressure(1.0, 0.0, 298.15).is_ok());
    }

    // ── Clausius-Clapeyron ───────────────────────────────────────────

    #[test]
    fn clausius_clapeyron_water() {
        // Water: P₁ = 101325 Pa at T₁ = 373.15 K, ΔH_vap = 40670 J/mol
        // At T₂ = 363.15 K (90°C), P₂ should be ~70 kPa
        let p2 = clausius_clapeyron(101325.0, 373.15, 363.15, 40670.0).unwrap();
        // 101325 × exp(−(40670/R)(1/363.15 − 1/373.15)) = 70 623.6 Pa
        // (measured: 70.1 kPa).
        assert!(
            (p2 - 70_623.55).abs() < 0.1,
            "water at 90°C should be ~70.6 kPa, got {p2}"
        );
    }

    #[test]
    fn clausius_clapeyron_rejects_nan() {
        assert!(clausius_clapeyron(f64::NAN, 373.15, 363.15, 40670.0).is_err());
        assert!(clausius_clapeyron(101325.0, f64::NAN, 363.15, 40670.0).is_err());
        assert!(clausius_clapeyron(101325.0, 373.15, f64::NAN, 40670.0).is_err());
        assert!(clausius_clapeyron(101325.0, 373.15, 363.15, f64::NAN).is_err());
    }

    #[test]
    fn clausius_clapeyron_same_temp() {
        let p2 = clausius_clapeyron(101325.0, 373.15, 373.15, 40670.0).unwrap();
        assert!((p2 - 101325.0).abs() < 1.0);
    }

    #[test]
    fn clausius_clapeyron_zero_temp_is_error() {
        assert!(clausius_clapeyron(101325.0, 0.0, 373.15, 40670.0).is_err());
    }

    // ── Antoine equation ─────────────────────────────────────────────

    #[test]
    fn antoine_water_100c() {
        let c = lookup_antoine("water").unwrap();
        let p = antoine_pressure_mmhg(c, 100.0);
        // Water at 100°C ≈ 760 mmHg
        assert!(
            (p - 760.0).abs() < 10.0,
            "water at 100°C should be ~760 mmHg, got {p}"
        );
    }

    #[test]
    fn antoine_water_pa() {
        let c = lookup_antoine("water").unwrap();
        let p = antoine_pressure_pa(c, 100.0);
        // 760.086 mmHg × 133.322 387 Pa/mmHg = 101 336.5 Pa
        assert!(
            (p - 101_336.53).abs() < 0.05,
            "water at 100°C should be ~101.3 kPa, got {p}"
        );
    }

    #[test]
    fn antoine_mmhg_to_pa_factor() {
        // 1 mmHg = 133.322 387 415 Pa (was truncated to 133.322).
        let c = lookup_antoine("water").unwrap();
        let ratio = antoine_pressure_pa(c, 50.0) / antoine_pressure_mmhg(c, 50.0);
        assert!((ratio - 133.322_387_415).abs() < 1e-9, "got {ratio}");
    }

    #[test]
    fn antoine_pole_and_non_finite_give_nan() {
        let c = lookup_antoine("water").unwrap();
        // T = −C is the pole; below it the formula is meaningless.
        assert!(antoine_pressure_mmhg(c, -c.c).is_nan());
        assert!(antoine_pressure_mmhg(c, -c.c - 10.0).is_nan());
        assert!(antoine_pressure_pa(c, -300.0).is_nan());
        assert!(antoine_pressure_mmhg(c, f64::NAN).is_nan());
        assert!(antoine_pressure_mmhg(c, f64::INFINITY).is_nan());
        let bad = AntoineCoeffs {
            a: f64::NAN,
            ..c.clone()
        };
        assert!(antoine_pressure_mmhg(&bad, 50.0).is_nan());
        // Just above the pole: finite and tiny.
        let p = antoine_pressure_mmhg(c, -c.c + 1.0);
        assert!(p.is_finite() && p >= 0.0);
    }

    #[test]
    fn antoine_methanol_matches_nist() {
        let c = lookup_antoine("methanol").unwrap();
        assert!((c.c - 239.726).abs() < 1e-9);
        assert!((c.t_min_c - 15.0).abs() < 1e-9 && (c.t_max_c - 84.0).abs() < 1e-9);
        // NIST WebBook (Ambrose & Sprake 1970): log10(P/bar) =
        // 5.20409 − 1581.341/(T/K − 33.50) → 16 938.45 Pa at 25 °C.
        let p = antoine_pressure_pa(c, 25.0);
        assert!((p - 16_938.45).abs() < 2.0, "got {p}");
        // Normal boiling point 64.7 °C → ~760 mmHg.
        let p_bp = antoine_pressure_mmhg(c, 64.7);
        assert!((p_bp - 760.0).abs() < 6.0, "got {p_bp}");
    }

    #[test]
    fn antoine_data_ranges_are_ordered() {
        for (name, c) in ANTOINE_DATA {
            assert!(c.t_min_c < c.t_max_c, "{name}");
            assert!(c.c + c.t_min_c > 0.0, "{name}: pole inside range");
        }
    }

    #[test]
    fn antoine_ethanol_78c() {
        let c = lookup_antoine("ethanol").unwrap();
        let p = antoine_pressure_mmhg(c, 78.4);
        // Ethanol boils at 78.4°C → ~760 mmHg
        assert!(
            (p - 760.0).abs() < 30.0,
            "ethanol at bp should be ~760 mmHg, got {p}"
        );
    }

    #[test]
    fn antoine_lookup_nonexistent() {
        assert!(lookup_antoine("unobtainium").is_none());
    }

    // ── Gibbs phase rule ─────────────────────────────────────────────

    #[test]
    fn phase_rule_water_triple_point() {
        // C=1, P=3 → F = 0 (invariant)
        let f = phase_rule(1, 3).unwrap();
        assert_eq!(f, 0);
    }

    #[test]
    fn phase_rule_water_liquid_vapor() {
        // C=1, P=2 → F = 1
        let f = phase_rule(1, 2).unwrap();
        assert_eq!(f, 1);
    }

    #[test]
    fn phase_rule_single_phase() {
        // C=1, P=1 → F = 2
        let f = phase_rule(1, 1).unwrap();
        assert_eq!(f, 2);
    }

    #[test]
    fn phase_rule_too_many_phases() {
        assert!(phase_rule(1, 4).is_err());
    }

    #[test]
    fn phase_rule_rejects_zero_phases() {
        assert!(phase_rule(2, 0).is_err());
        // Used to overflow i32 (panic in debug builds).
        assert!(phase_rule(i32::MAX as u32, 0).is_err());
    }

    #[test]
    fn phase_rule_large_counts() {
        // `as i32` used to wrap: u32::MAX → −1 gave Ok(0) here.
        assert!(phase_rule(u32::MAX, 1).is_err());
        assert_eq!(phase_rule(u32::MAX, 2).unwrap(), u32::MAX);
        assert_eq!(phase_rule(u32::MAX, u32::MAX).unwrap(), 2);
        assert_eq!(phase_rule(3_000_000_000, 1).unwrap(), 3_000_000_001);
        assert!(phase_rule(0, u32::MAX).is_err());
    }
}
