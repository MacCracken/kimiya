//! Phase equilibria — Raoult's law, Henry's law, colligative properties,
//! Clausius-Clapeyron, Antoine equation.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use serde::Serialize;

// ── Raoult's law ─────────────────────────────────────────────────────

/// Raoult's law: P = x × P°
///
/// Vapor pressure of a component in an ideal solution.
///
/// - `mole_fraction`: x (0 to 1)
/// - `pure_vapor_pressure`: P° of pure component (Pa)
///
/// # Errors
///
/// Returns error if mole fraction is not in \[0, 1\] or pressure is negative.
#[inline]
pub fn raoult_pressure(mole_fraction: f64, pure_vapor_pressure: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&mole_fraction) {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be between 0 and 1".into(),
        ));
    }
    if pure_vapor_pressure < 0.0 {
        return Err(KimiyaError::InvalidInput(
            "vapor pressure must be non-negative".into(),
        ));
    }
    Ok(mole_fraction * pure_vapor_pressure)
}

/// Total vapor pressure over a binary ideal solution (Raoult's law).
///
/// P_total = x_A × P°_A + (1 - x_A) × P°_B
#[inline]
pub fn raoult_total_binary(x_a: f64, p_pure_a: f64, p_pure_b: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&x_a) {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be between 0 and 1".into(),
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
/// - `mole_fraction`: x of dissolved gas
///
/// # Errors
///
/// Returns error if inputs are invalid.
#[inline]
pub fn henry_pressure(henry_constant: f64, mole_fraction: f64) -> Result<f64> {
    if henry_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "Henry's constant must be positive".into(),
        ));
    }
    if mole_fraction < 0.0 {
        return Err(KimiyaError::InvalidInput(
            "mole fraction must be non-negative".into(),
        ));
    }
    Ok(henry_constant * mole_fraction)
}

/// Solubility (mole fraction) of a gas from Henry's law: x = P / K_H
#[inline]
pub fn henry_solubility(partial_pressure: f64, henry_constant: f64) -> Result<f64> {
    if henry_constant <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "Henry's constant must be positive".into(),
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
#[must_use]
#[inline]
pub fn boiling_point_elevation(kb: f64, molality: f64, van_t_hoff_factor: f64) -> f64 {
    kb * molality * van_t_hoff_factor
}

/// Freezing point depression: ΔT_f = K_f × m × i
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
#[inline]
pub fn osmotic_pressure(van_t_hoff_factor: f64, molarity: f64, temperature_k: f64) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    // R in L·Pa/(mol·K) = 8314.462618 (GAS_CONSTANT × 1000 for L→m³)
    Ok(van_t_hoff_factor * molarity * GAS_CONSTANT * 1000.0 * temperature_k)
}

/// Common ebullioscopic constants K_b (K·kg/mol).
pub const KB_WATER: f64 = 0.512;
pub const KB_BENZENE: f64 = 2.53;
pub const KB_ETHANOL: f64 = 1.22;
pub const KB_CHLOROFORM: f64 = 3.63;

/// Common cryoscopic constants K_f (K·kg/mol).
pub const KF_WATER: f64 = 1.86;
pub const KF_BENZENE: f64 = 5.12;
pub const KF_CYCLOHEXANE: f64 = 20.0;
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
/// # Errors
///
/// Returns error if temperatures or pressure are not positive.
pub fn clausius_clapeyron(p1: f64, t1_k: f64, t2_k: f64, delta_h_vap_j: f64) -> Result<f64> {
    if p1 <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "reference pressure must be positive".into(),
        ));
    }
    if t1_k <= 0.0 || t2_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
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
    pub a: f64,
    pub b: f64,
    pub c: f64,
    /// Valid temperature range (°C).
    pub t_min_c: f64,
    pub t_max_c: f64,
}

/// Vapor pressure from Antoine equation.
///
/// Returns pressure in mmHg.
#[must_use]
#[inline]
pub fn antoine_pressure_mmhg(coeffs: &AntoineCoeffs, temperature_c: f64) -> f64 {
    10.0_f64.powf(coeffs.a - coeffs.b / (coeffs.c + temperature_c))
}

/// Vapor pressure from Antoine equation in Pa.
#[must_use]
#[inline]
pub fn antoine_pressure_pa(coeffs: &AntoineCoeffs, temperature_c: f64) -> f64 {
    antoine_pressure_mmhg(coeffs, temperature_c) * 133.322
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
        AntoineCoeffs {
            a: 8.08097,
            b: 1582.27,
            c: 239.700,
            t_min_c: -14.0,
            t_max_c: 65.0,
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
/// Returns error if phases > components + 2 (F would be negative).
#[inline]
pub fn phase_rule(components: u32, phases: u32) -> Result<u32> {
    let f = components as i32 - phases as i32 + 2;
    if f < 0 {
        return Err(KimiyaError::InvalidInput(
            "too many phases for the number of components".into(),
        ));
    }
    Ok(f as u32)
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
        assert!((p - 80.0).abs() < f64::EPSILON);
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
        // π ≈ 2 × 0.1 × 8314.46 × 298.15 ≈ 495800 Pa ≈ 4.9 atm
        assert!(pi > 400_000.0 && pi < 600_000.0);
    }

    #[test]
    fn osmotic_pressure_zero_temp_is_error() {
        assert!(osmotic_pressure(1.0, 0.1, 0.0).is_err());
    }

    // ── Clausius-Clapeyron ───────────────────────────────────────────

    #[test]
    fn clausius_clapeyron_water() {
        // Water: P₁ = 101325 Pa at T₁ = 373.15 K, ΔH_vap = 40670 J/mol
        // At T₂ = 363.15 K (90°C), P₂ should be ~70 kPa
        let p2 = clausius_clapeyron(101325.0, 373.15, 363.15, 40670.0).unwrap();
        assert!(
            p2 > 60_000.0 && p2 < 80_000.0,
            "water at 90°C should be ~70 kPa, got {p2}"
        );
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
        // ~101325 Pa
        assert!(
            (p - 101325.0).abs() < 2000.0,
            "water at 100°C should be ~101 kPa, got {p}"
        );
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
}
