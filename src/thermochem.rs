//! Thermochemical data — standard formation enthalpies, Gibbs energies, entropies.
//!
//! All values at standard conditions: 25°C (298.15 K), 1 atm.
//! Reference state: elements in their most stable form have ΔH_f° = 0, ΔG_f° = 0.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use serde::Serialize;

// ── Thermochemical data table ────────────────────────────────────────

/// Standard thermochemical data for a substance at 25°C, 1 atm.
#[derive(Debug, Clone, Serialize)]
pub struct ThermochemData {
    /// Chemical formula (e.g., "H2O(l)").
    pub formula: &'static str,
    /// Common name.
    pub name: &'static str,
    /// Standard enthalpy of formation ΔH_f° (kJ/mol).
    pub delta_hf_kj: f64,
    /// Standard Gibbs free energy of formation ΔG_f° (kJ/mol).
    pub delta_gf_kj: f64,
    /// Standard molar entropy S° (J/(mol·K)).
    pub s_standard_j: f64,
}

/// Built-in thermochemical data for common substances.
///
/// Values from NIST/JANAF and CRC Handbook.
/// Elements in standard state have ΔH_f° = 0 and ΔG_f° = 0.
pub static THERMOCHEM_DATA: &[ThermochemData] = &[
    // ── Elements in standard state ───────────────────────────────
    ThermochemData {
        formula: "H2(g)",
        name: "Hydrogen gas",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 130.68,
    },
    ThermochemData {
        formula: "O2(g)",
        name: "Oxygen gas",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 205.15,
    },
    ThermochemData {
        formula: "N2(g)",
        name: "Nitrogen gas",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 191.61,
    },
    ThermochemData {
        formula: "C(graphite)",
        name: "Graphite",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 5.74,
    },
    ThermochemData {
        formula: "Fe(s)",
        name: "Iron",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 27.28,
    },
    ThermochemData {
        formula: "Cu(s)",
        name: "Copper",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 33.15,
    },
    ThermochemData {
        formula: "Al(s)",
        name: "Aluminum",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 28.30,
    },
    // ── Common compounds ─────────────────────────────────────────
    ThermochemData {
        formula: "H2O(l)",
        name: "Water (liquid)",
        delta_hf_kj: -285.83,
        delta_gf_kj: -237.13,
        s_standard_j: 69.91,
    },
    ThermochemData {
        formula: "H2O(g)",
        name: "Water (gas)",
        delta_hf_kj: -241.82,
        delta_gf_kj: -228.57,
        s_standard_j: 188.83,
    },
    ThermochemData {
        formula: "CO2(g)",
        name: "Carbon dioxide",
        delta_hf_kj: -393.51,
        delta_gf_kj: -394.36,
        s_standard_j: 213.79,
    },
    ThermochemData {
        formula: "CO(g)",
        name: "Carbon monoxide",
        delta_hf_kj: -110.53,
        delta_gf_kj: -137.17,
        s_standard_j: 197.67,
    },
    ThermochemData {
        formula: "CH4(g)",
        name: "Methane",
        delta_hf_kj: -74.81,
        delta_gf_kj: -50.72,
        s_standard_j: 186.26,
    },
    ThermochemData {
        formula: "C2H6(g)",
        name: "Ethane",
        delta_hf_kj: -84.68,
        delta_gf_kj: -32.82,
        s_standard_j: 229.60,
    },
    ThermochemData {
        formula: "C2H4(g)",
        name: "Ethylene",
        delta_hf_kj: 52.26,
        delta_gf_kj: 68.15,
        s_standard_j: 219.56,
    },
    ThermochemData {
        formula: "C2H2(g)",
        name: "Acetylene",
        delta_hf_kj: 226.73,
        delta_gf_kj: 209.20,
        s_standard_j: 200.94,
    },
    ThermochemData {
        formula: "C6H12O6(s)",
        name: "Glucose",
        delta_hf_kj: -1274.4,
        delta_gf_kj: -910.3,
        s_standard_j: 212.1,
    },
    ThermochemData {
        formula: "C2H5OH(l)",
        name: "Ethanol (liquid)",
        delta_hf_kj: -277.69,
        delta_gf_kj: -174.78,
        s_standard_j: 160.70,
    },
    ThermochemData {
        formula: "NH3(g)",
        name: "Ammonia",
        delta_hf_kj: -46.11,
        delta_gf_kj: -16.45,
        s_standard_j: 192.45,
    },
    ThermochemData {
        formula: "NO(g)",
        name: "Nitric oxide",
        delta_hf_kj: 90.25,
        delta_gf_kj: 86.55,
        s_standard_j: 210.76,
    },
    ThermochemData {
        formula: "NO2(g)",
        name: "Nitrogen dioxide",
        delta_hf_kj: 33.18,
        delta_gf_kj: 51.31,
        s_standard_j: 240.06,
    },
    ThermochemData {
        formula: "SO2(g)",
        name: "Sulfur dioxide",
        delta_hf_kj: -296.83,
        delta_gf_kj: -300.19,
        s_standard_j: 248.22,
    },
    ThermochemData {
        formula: "SO3(g)",
        name: "Sulfur trioxide",
        delta_hf_kj: -395.72,
        delta_gf_kj: -371.06,
        s_standard_j: 256.76,
    },
    ThermochemData {
        formula: "HCl(g)",
        name: "Hydrogen chloride",
        delta_hf_kj: -92.31,
        delta_gf_kj: -95.30,
        s_standard_j: 186.91,
    },
    ThermochemData {
        formula: "NaCl(s)",
        name: "Sodium chloride",
        delta_hf_kj: -411.15,
        delta_gf_kj: -384.14,
        s_standard_j: 72.13,
    },
    ThermochemData {
        formula: "CaCO3(s)",
        name: "Calcium carbonate",
        delta_hf_kj: -1206.9,
        delta_gf_kj: -1128.8,
        s_standard_j: 92.9,
    },
    ThermochemData {
        formula: "Fe2O3(s)",
        name: "Iron(III) oxide",
        delta_hf_kj: -824.2,
        delta_gf_kj: -742.2,
        s_standard_j: 87.40,
    },
    ThermochemData {
        formula: "Al2O3(s)",
        name: "Aluminum oxide",
        delta_hf_kj: -1675.7,
        delta_gf_kj: -1582.3,
        s_standard_j: 50.92,
    },
    ThermochemData {
        formula: "CaO(s)",
        name: "Calcium oxide",
        delta_hf_kj: -635.09,
        delta_gf_kj: -604.03,
        s_standard_j: 39.75,
    },
    ThermochemData {
        formula: "H2SO4(l)",
        name: "Sulfuric acid",
        delta_hf_kj: -813.99,
        delta_gf_kj: -690.00,
        s_standard_j: 156.90,
    },
    ThermochemData {
        formula: "HNO3(l)",
        name: "Nitric acid",
        delta_hf_kj: -174.10,
        delta_gf_kj: -80.71,
        s_standard_j: 155.60,
    },
    ThermochemData {
        formula: "NaOH(s)",
        name: "Sodium hydroxide",
        delta_hf_kj: -425.61,
        delta_gf_kj: -379.49,
        s_standard_j: 64.46,
    },
];

/// Look up thermochemical data by formula (e.g., "H2O(l)", "CO2(g)").
#[must_use]
#[inline]
pub fn lookup_thermochem(formula: &str) -> Option<&'static ThermochemData> {
    THERMOCHEM_DATA.iter().find(|d| d.formula == formula)
}

// ── Reaction thermochemistry from tables ─────────────────────────────

/// Standard reaction enthalpy from formation enthalpies:
/// ΔH°_rxn = Σ n_i·ΔH_f°(products) − Σ n_j·ΔH_f°(reactants)
///
/// Each entry is (formula, stoichiometric coefficient).
///
/// # Errors
///
/// Returns error if any formula is not found in the thermochemical data.
pub fn reaction_enthalpy(products: &[(&str, f64)], reactants: &[(&str, f64)]) -> Result<f64> {
    let sum_products = sum_property(products, |d| d.delta_hf_kj)?;
    let sum_reactants = sum_property(reactants, |d| d.delta_hf_kj)?;
    Ok(sum_products - sum_reactants)
}

/// Standard reaction Gibbs energy:
/// ΔG°_rxn = Σ n_i·ΔG_f°(products) − Σ n_j·ΔG_f°(reactants)
///
/// # Errors
///
/// Returns error if any formula is not found in the thermochemical data.
pub fn reaction_gibbs_energy(products: &[(&str, f64)], reactants: &[(&str, f64)]) -> Result<f64> {
    let sum_products = sum_property(products, |d| d.delta_gf_kj)?;
    let sum_reactants = sum_property(reactants, |d| d.delta_gf_kj)?;
    Ok(sum_products - sum_reactants)
}

/// Standard reaction entropy change:
/// ΔS°_rxn = Σ n_i·S°(products) − Σ n_j·S°(reactants)
///
/// Returns entropy in J/(mol·K).
///
/// # Errors
///
/// Returns error if any formula is not found in the thermochemical data.
pub fn reaction_entropy(products: &[(&str, f64)], reactants: &[(&str, f64)]) -> Result<f64> {
    let sum_products = sum_property(products, |d| d.s_standard_j)?;
    let sum_reactants = sum_property(reactants, |d| d.s_standard_j)?;
    Ok(sum_products - sum_reactants)
}

/// Helper: sum n_i × property(formula_i) for a list of (formula, coefficient) pairs.
fn sum_property(species: &[(&str, f64)], prop: impl Fn(&ThermochemData) -> f64) -> Result<f64> {
    let mut total = 0.0;
    for &(formula, coeff) in species {
        let data = lookup_thermochem(formula).ok_or_else(|| {
            KimiyaError::InvalidReaction(format!("no thermochemical data for '{formula}'"))
        })?;
        total += coeff * prop(data);
    }
    Ok(total)
}

// ── Van't Hoff equation ──────────────────────────────────────────────

/// Van't Hoff equation: ln(K₂/K₁) = −(ΔH°/R)(1/T₂ − 1/T₁)
///
/// Calculates equilibrium constant at a new temperature given K at a reference
/// temperature, assuming ΔH° is constant over the temperature range.
///
/// - `k_ref`: equilibrium constant at `t_ref_k`
/// - `delta_h_j`: standard reaction enthalpy in J/mol (not kJ!)
/// - `t_ref_k`: reference temperature (K)
/// - `t_new_k`: target temperature (K)
///
/// # Errors
///
/// Returns error if temperatures or K_ref are not positive.
pub fn vant_hoff_k(k_ref: f64, delta_h_j: f64, t_ref_k: f64, t_new_k: f64) -> Result<f64> {
    if t_ref_k <= 0.0 || t_new_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
        ));
    }
    if k_ref <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "reference K must be positive".into(),
        ));
    }
    let exponent = -(delta_h_j / GAS_CONSTANT) * (1.0 / t_new_k - 1.0 / t_ref_k);
    Ok(k_ref * exponent.exp())
}

// ── Heat capacity integration ────────────────────────────────────────

/// Shomate equation coefficients for Cp(T) in J/(mol·K).
///
/// Cp(T) = A + B·t + C·t² + D·t³ + E/t²
///
/// where t = T(K) / 1000.
#[derive(Debug, Clone, Serialize)]
pub struct Shomate {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
    /// Valid temperature range (K).
    pub t_min: f64,
    pub t_max: f64,
}

impl Shomate {
    /// Evaluate Cp at temperature T (K) in J/(mol·K).
    #[must_use]
    #[inline]
    pub fn cp(&self, temperature_k: f64) -> f64 {
        let t = temperature_k / 1000.0;
        self.a + self.b * t + self.c * t * t + self.d * t * t * t + self.e / (t * t)
    }
}

/// Built-in Shomate coefficients for common gases (NIST, 298–1200 K range).
pub static SHOMATE_DATA: &[(&str, Shomate)] = &[
    (
        "H2O(g)",
        Shomate {
            a: 30.092,
            b: 6.832514,
            c: 6.793435,
            d: -2.534480,
            e: 0.082139,
            t_min: 500.0,
            t_max: 1700.0,
        },
    ),
    (
        "CO2(g)",
        Shomate {
            a: 24.99735,
            b: 55.18696,
            c: -33.69137,
            d: 7.948387,
            e: -0.136638,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "N2(g)",
        Shomate {
            a: 28.98641,
            b: 1.853978,
            c: -9.647459,
            d: 16.63537,
            e: 0.000117,
            t_min: 298.0,
            t_max: 1400.0,
        },
    ),
    (
        "O2(g)",
        Shomate {
            a: 31.32234,
            b: -20.23531,
            c: 57.86644,
            d: -36.50624,
            e: -0.007374,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "CH4(g)",
        Shomate {
            a: -0.703029,
            b: 108.4773,
            c: -42.52157,
            d: 5.862788,
            e: 0.678565,
            t_min: 298.0,
            t_max: 1300.0,
        },
    ),
    (
        "CO(g)",
        Shomate {
            a: 25.56759,
            b: 6.096130,
            c: 4.054656,
            d: -2.671301,
            e: 0.131021,
            t_min: 298.0,
            t_max: 1300.0,
        },
    ),
];

/// Look up Shomate coefficients by formula.
#[must_use]
#[inline]
pub fn lookup_shomate(formula: &str) -> Option<&'static Shomate> {
    SHOMATE_DATA
        .iter()
        .find(|(f, _)| *f == formula)
        .map(|(_, s)| s)
}

/// Enthalpy change from T₁ to T₂ via Cp(T) integration:
/// ΔH = ∫[T₁→T₂] Cp(T) dT
///
/// Uses [`hisab::calc::integral_simpson`] for numerical integration of
/// the Shomate polynomial.
///
/// Returns enthalpy change in J/mol.
///
/// # Errors
///
/// Returns error if formula not found, temperatures invalid, or integration fails.
pub fn enthalpy_change_cp(formula: &str, t1_k: f64, t2_k: f64) -> Result<f64> {
    let shomate = lookup_shomate(formula)
        .ok_or_else(|| KimiyaError::InvalidReaction(format!("no Shomate data for '{formula}'")))?;
    if t1_k <= 0.0 || t2_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
        ));
    }
    let cp_fn = |t: f64| shomate.cp(t);
    let result = hisab::calc::integral_simpson(cp_fn, t1_k, t2_k, 100)
        .map_err(|e| KimiyaError::ComputationError(format!("integration failed: {e}")))?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Data lookup ──────────────────────────────────────────────────

    #[test]
    fn lookup_water_liquid() {
        let d = lookup_thermochem("H2O(l)").unwrap();
        assert!((d.delta_hf_kj - (-285.83)).abs() < 0.01);
        assert!((d.delta_gf_kj - (-237.13)).abs() < 0.01);
        assert!((d.s_standard_j - 69.91).abs() < 0.01);
    }

    #[test]
    fn lookup_co2() {
        let d = lookup_thermochem("CO2(g)").unwrap();
        assert!((d.delta_hf_kj - (-393.51)).abs() < 0.01);
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_thermochem("XeF6(g)").is_none());
    }

    #[test]
    fn elements_have_zero_formation() {
        for formula in &["H2(g)", "O2(g)", "N2(g)", "C(graphite)", "Fe(s)"] {
            let d = lookup_thermochem(formula).unwrap();
            assert!(
                d.delta_hf_kj.abs() < f64::EPSILON,
                "{formula} should have ΔH_f° = 0"
            );
            assert!(
                d.delta_gf_kj.abs() < f64::EPSILON,
                "{formula} should have ΔG_f° = 0"
            );
        }
    }

    #[test]
    fn entropy_always_positive() {
        for d in THERMOCHEM_DATA.iter() {
            assert!(
                d.s_standard_j > 0.0,
                "{} has non-positive entropy: {}",
                d.formula,
                d.s_standard_j
            );
        }
    }

    // ── Reaction calculations ────────────────────────────────────────

    #[test]
    fn combustion_of_methane_enthalpy() {
        // CH₄(g) + 2O₂(g) → CO₂(g) + 2H₂O(l)
        // ΔH° = [-393.51 + 2(-285.83)] - [-74.81 + 0] = -890.36 kJ/mol
        let dh = reaction_enthalpy(
            &[("CO2(g)", 1.0), ("H2O(l)", 2.0)],
            &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
        )
        .unwrap();
        assert!(
            (dh - (-890.36)).abs() < 0.1,
            "CH₄ combustion should be ~-890 kJ/mol, got {dh}"
        );
    }

    #[test]
    fn combustion_of_methane_gibbs() {
        // ΔG° = [-394.36 + 2(-237.13)] - [-50.72 + 0] = -817.90 kJ/mol
        let dg = reaction_gibbs_energy(
            &[("CO2(g)", 1.0), ("H2O(l)", 2.0)],
            &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
        )
        .unwrap();
        assert!(
            (dg - (-817.90)).abs() < 0.1,
            "CH₄ combustion ΔG° should be ~-818 kJ/mol, got {dg}"
        );
    }

    #[test]
    fn combustion_of_methane_entropy() {
        // ΔS° = [213.79 + 2(69.91)] - [186.26 + 2(205.15)]
        //      = 353.61 - 596.56 = -242.95 J/(mol·K)
        let ds = reaction_entropy(
            &[("CO2(g)", 1.0), ("H2O(l)", 2.0)],
            &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
        )
        .unwrap();
        assert!(
            (ds - (-242.95)).abs() < 0.1,
            "CH₄ combustion ΔS° should be ~-243 J/(mol·K), got {ds}"
        );
    }

    #[test]
    fn reaction_unknown_formula_is_error() {
        assert!(reaction_enthalpy(&[("XeF6(g)", 1.0)], &[("H2(g)", 1.0)]).is_err());
    }

    #[test]
    fn haber_process_enthalpy() {
        // N₂(g) + 3H₂(g) → 2NH₃(g)
        // ΔH° = 2(-46.11) - [0 + 0] = -92.22 kJ/mol
        let dh = reaction_enthalpy(&[("NH3(g)", 2.0)], &[("N2(g)", 1.0), ("H2(g)", 3.0)]).unwrap();
        assert!(
            (dh - (-92.22)).abs() < 0.1,
            "Haber process should be ~-92 kJ/mol, got {dh}"
        );
    }

    // ── Gibbs-enthalpy-entropy consistency ────────────────────────────

    #[test]
    fn gibbs_enthalpy_entropy_consistency() {
        // For any compound: ΔG_f° ≈ ΔH_f° - T·ΔS_f°
        // where ΔS_f° = S°(compound) - Σ S°(elements)
        // Check water: ΔG = ΔH - TΔS
        // ΔH_f° = -285.83, ΔG_f° = -237.13
        // At 298.15 K: TΔS_f° = ΔH_f° - ΔG_f° = -285.83 - (-237.13) = -48.70 kJ
        // ΔS_f° = -48700 / 298.15 = -163.3 J/(mol·K)
        // S°(H2O) - S°(H2) - 0.5·S°(O2) = 69.91 - 130.68 - 0.5(205.15) = -163.34
        let d = lookup_thermochem("H2O(l)").unwrap();
        let t = 298.15;
        let t_delta_s = d.delta_hf_kj - d.delta_gf_kj; // kJ
        let delta_s_f = t_delta_s * 1000.0 / t; // J/(mol·K)

        // Calculate ΔS_f° from element entropies
        let s_h2 = lookup_thermochem("H2(g)").unwrap().s_standard_j;
        let s_o2 = lookup_thermochem("O2(g)").unwrap().s_standard_j;
        let delta_s_calc = d.s_standard_j - s_h2 - 0.5 * s_o2;

        assert!(
            (delta_s_f - delta_s_calc).abs() < 0.5,
            "ΔS_f° from ΔG/ΔH ({delta_s_f:.1}) should match element entropies ({delta_s_calc:.1})"
        );
    }

    // ── Van't Hoff equation ──────────────────────────────────────────

    #[test]
    fn vant_hoff_same_temperature() {
        let k2 = vant_hoff_k(1.0, -50_000.0, 298.15, 298.15).unwrap();
        assert!(
            (k2 - 1.0).abs() < 1e-10,
            "same temp should give same K, got {k2}"
        );
    }

    #[test]
    fn vant_hoff_exothermic_higher_temp_decreases_k() {
        // Exothermic (ΔH < 0): raising T decreases K
        let k1 = 100.0;
        let k2 = vant_hoff_k(k1, -50_000.0, 298.15, 400.0).unwrap();
        assert!(k2 < k1, "exothermic: higher T should decrease K, got {k2}");
    }

    #[test]
    fn vant_hoff_endothermic_higher_temp_increases_k() {
        // Endothermic (ΔH > 0): raising T increases K
        let k1 = 0.01;
        let k2 = vant_hoff_k(k1, 50_000.0, 298.15, 400.0).unwrap();
        assert!(k2 > k1, "endothermic: higher T should increase K, got {k2}");
    }

    #[test]
    fn vant_hoff_zero_temp_is_error() {
        assert!(vant_hoff_k(1.0, -50_000.0, 0.0, 300.0).is_err());
        assert!(vant_hoff_k(1.0, -50_000.0, 300.0, 0.0).is_err());
    }

    #[test]
    fn vant_hoff_zero_k_is_error() {
        assert!(vant_hoff_k(0.0, -50_000.0, 298.15, 400.0).is_err());
    }

    // ── Audit-driven tests ───────────────────────────────────────────

    #[test]
    fn empty_reaction_returns_zero() {
        let dh = reaction_enthalpy(&[], &[]).unwrap();
        assert!(dh.abs() < f64::EPSILON);
    }

    #[test]
    fn vant_hoff_extreme_enthalpy_doesnt_panic() {
        // Very large ΔH can cause overflow to infinity — should not panic
        let result = vant_hoff_k(1.0, 1e8, 298.15, 1000.0);
        assert!(result.is_ok());
        let k = result.unwrap();
        assert!(k.is_finite() || k.is_infinite()); // either is acceptable, no NaN
        assert!(!k.is_nan());
    }

    #[test]
    fn thermodynamic_consistency_co2() {
        // Verify ΔG ≈ ΔH - TΔS for CO2 formation
        // C(graphite) + O2(g) → CO2(g)
        let d = lookup_thermochem("CO2(g)").unwrap();
        let s_c = lookup_thermochem("C(graphite)").unwrap().s_standard_j;
        let s_o2 = lookup_thermochem("O2(g)").unwrap().s_standard_j;
        let delta_s = d.s_standard_j - s_c - s_o2; // J/(mol·K)
        let t = 298.15;
        let dg_calc = d.delta_hf_kj - t * delta_s / 1000.0; // kJ
        assert!(
            (dg_calc - d.delta_gf_kj).abs() < 1.0,
            "CO2: ΔG_calc={dg_calc:.1}, ΔG_tab={:.1}",
            d.delta_gf_kj
        );
    }

    #[test]
    fn unique_formulas() {
        // No duplicate formula entries
        for (i, a) in THERMOCHEM_DATA.iter().enumerate() {
            for b in THERMOCHEM_DATA.iter().skip(i + 1) {
                assert_ne!(a.formula, b.formula, "duplicate formula: {}", a.formula);
            }
        }
    }

    #[test]
    fn water_vaporization_enthalpy() {
        // ΔH_vap = ΔH_f°(H2O,g) - ΔH_f°(H2O,l)
        let liquid = lookup_thermochem("H2O(l)").unwrap();
        let gas = lookup_thermochem("H2O(g)").unwrap();
        let dh_vap = gas.delta_hf_kj - liquid.delta_hf_kj;
        // Known: ΔH_vap(H2O) ≈ 44.01 kJ/mol at 25°C
        assert!(
            (dh_vap - 44.01).abs() < 0.1,
            "water vaporization should be ~44 kJ/mol, got {dh_vap}"
        );
    }

    #[test]
    fn gases_have_higher_entropy_than_solids() {
        // S°(H2O,g) > S°(H2O,l) — fundamental thermodynamic principle
        let liquid = lookup_thermochem("H2O(l)").unwrap();
        let gas = lookup_thermochem("H2O(g)").unwrap();
        assert!(
            gas.s_standard_j > liquid.s_standard_j,
            "gas entropy should exceed liquid"
        );
    }

    // ── Shomate / Cp(T) tests (hisab-powered) ───────────────────────

    #[test]
    fn shomate_co2_cp_at_298() {
        // CO2 Cp at 298 K ≈ 37.1 J/(mol·K)
        let s = lookup_shomate("CO2(g)").unwrap();
        let cp = s.cp(298.0);
        assert!(
            (cp - 37.1).abs() < 1.0,
            "CO2 Cp at 298K should be ~37 J/(mol·K), got {cp}"
        );
    }

    #[test]
    fn shomate_n2_cp_at_298() {
        // N2 Cp at 298 K ≈ 29.1 J/(mol·K)
        let s = lookup_shomate("N2(g)").unwrap();
        let cp = s.cp(298.0);
        assert!(
            (cp - 29.1).abs() < 0.5,
            "N2 Cp at 298K should be ~29.1 J/(mol·K), got {cp}"
        );
    }

    #[test]
    fn shomate_cp_increases_with_temperature() {
        // For polyatomic molecules, Cp generally increases with T
        let s = lookup_shomate("CO2(g)").unwrap();
        let cp_300 = s.cp(300.0);
        let cp_1000 = s.cp(1000.0);
        assert!(
            cp_1000 > cp_300,
            "CO2 Cp should increase with T: {cp_300} at 300K vs {cp_1000} at 1000K"
        );
    }

    #[test]
    fn enthalpy_change_co2_heating() {
        // Heat CO2 from 300K to 600K
        let dh = enthalpy_change_cp("CO2(g)", 300.0, 600.0).unwrap();
        // Rough estimate: ~37 J/(mol·K) × 300 K ≈ 11,100 J/mol
        assert!(
            dh > 10_000.0 && dh < 15_000.0,
            "CO2 300→600K should be ~11-13 kJ/mol, got {dh}"
        );
    }

    #[test]
    fn enthalpy_change_zero_range() {
        let dh = enthalpy_change_cp("CO2(g)", 500.0, 500.0).unwrap();
        assert!(dh.abs() < 0.01, "same T should give ΔH≈0, got {dh}");
    }

    #[test]
    fn enthalpy_change_cooling_is_negative() {
        let dh = enthalpy_change_cp("N2(g)", 600.0, 300.0).unwrap();
        assert!(dh < 0.0, "cooling should give negative ΔH, got {dh}");
    }

    #[test]
    fn enthalpy_change_unknown_formula_is_error() {
        assert!(enthalpy_change_cp("XeF6(g)", 300.0, 600.0).is_err());
    }

    #[test]
    fn enthalpy_change_zero_temp_is_error() {
        assert!(enthalpy_change_cp("CO2(g)", 0.0, 600.0).is_err());
    }

    #[test]
    fn shomate_lookup_nonexistent() {
        assert!(lookup_shomate("XeF6(g)").is_none());
    }

    #[test]
    fn shomate_cp_positive_in_valid_range() {
        // Cp should always be positive in the valid temperature range
        for (formula, s) in SHOMATE_DATA.iter() {
            for t in [s.t_min, (s.t_min + s.t_max) / 2.0, s.t_max] {
                let cp = s.cp(t);
                assert!(
                    cp > 0.0,
                    "{formula}: Cp should be positive at {t}K, got {cp}"
                );
            }
        }
    }

    #[test]
    fn enthalpy_change_cp_n2_matches_constant_cp_approx() {
        // N2 has nearly constant Cp ≈ 29.1 J/(mol·K)
        // ΔH from 300→600K ≈ 29.1 × 300 = 8730 J/mol
        let dh = enthalpy_change_cp("N2(g)", 300.0, 600.0).unwrap();
        assert!(
            (dh - 8730.0).abs() < 500.0,
            "N2 300→600K should be ~8730 J/mol, got {dh}"
        );
    }

    #[test]
    fn all_shomate_data_has_valid_range() {
        for (formula, s) in SHOMATE_DATA.iter() {
            assert!(
                s.t_min < s.t_max,
                "{formula}: t_min ({}) must be < t_max ({})",
                s.t_min,
                s.t_max
            );
            assert!(s.t_min > 0.0, "{formula}: t_min must be positive");
        }
    }
}
