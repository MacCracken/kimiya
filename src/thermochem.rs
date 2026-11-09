//! Thermochemical data — standard formation enthalpies, Gibbs energies,
//! entropies, and Shomate heat-capacity fits.
//!
//! Tabulated values are at 298.15 K and the IUPAC standard pressure
//! p° = 1 bar (100 kPa), the convention of the NBS, CODATA, and CRC tables.
//! Reference state: elements in their most stable form have ΔH_f° = 0, ΔG_f° = 0.
//!
//! # Units
//!
//! Energy units differ between functions; check before chaining them:
//!
//! | Function | Energy unit |
//! |---|---|
//! | [`reaction_enthalpy`], [`reaction_gibbs_energy`] | returns **kJ/mol** |
//! | [`reaction_entropy`] | returns J/(mol·K) |
//! | [`enthalpy_change_cp`] | returns **J/mol** |
//! | [`vant_hoff_k`], [`adiabatic_flame_temperature`] | take **J/mol** |
//!
//! Multiply a kJ/mol result by 1000 before passing it to a function that
//! takes J/mol. "Per mol" means per mole of reaction as written.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_positive};
use serde::Serialize;

// ── Thermochemical data table ────────────────────────────────────────

/// Standard thermochemical data for a substance at 298.15 K and 1 bar.
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

/// Built-in thermochemical data for common substances at 298.15 K, 1 bar.
///
/// ΔH_f° and ΔG_f° in kJ/mol, S° in J/(mol·K). Values from the NBS tables
/// (Wagman et al. 1982), NIST/JANAF, CODATA, and the CRC Handbook.
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
    // ── Elements in standard state (additional) ─────────────────
    ThermochemData {
        formula: "S(rhombic)",
        name: "Sulfur (rhombic)",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 31.80,
    },
    ThermochemData {
        formula: "Cl2(g)",
        name: "Chlorine gas",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 223.08,
    },
    ThermochemData {
        formula: "Br2(l)",
        name: "Bromine liquid",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 152.23,
    },
    ThermochemData {
        formula: "I2(s)",
        name: "Iodine solid",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 116.14,
    },
    ThermochemData {
        formula: "Na(s)",
        name: "Sodium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 51.21,
    },
    ThermochemData {
        formula: "K(s)",
        name: "Potassium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 64.18,
    },
    ThermochemData {
        formula: "Ca(s)",
        name: "Calcium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 41.42,
    },
    ThermochemData {
        formula: "Mg(s)",
        name: "Magnesium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 32.68,
    },
    ThermochemData {
        formula: "P(white)",
        name: "Phosphorus (white)",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 41.09,
    },
    ThermochemData {
        formula: "Si(s)",
        name: "Silicon",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 18.83,
    },
    // S° from CODATA key values (Zn, Ti, Ag, F2) and the CRC Handbook (Ba).
    ThermochemData {
        formula: "Zn(s)",
        name: "Zinc",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 41.63,
    },
    ThermochemData {
        formula: "Ti(s)",
        name: "Titanium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 30.72,
    },
    ThermochemData {
        formula: "Ba(s)",
        name: "Barium",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 62.5,
    },
    ThermochemData {
        formula: "Ag(s)",
        name: "Silver",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 42.55,
    },
    ThermochemData {
        formula: "F2(g)",
        name: "Fluorine gas",
        delta_hf_kj: 0.0,
        delta_gf_kj: 0.0,
        s_standard_j: 202.79,
    },
    // ── Inorganic compounds (additional) ──────────────────────────
    ThermochemData {
        formula: "MgO(s)",
        name: "Magnesium oxide",
        delta_hf_kj: -601.60,
        delta_gf_kj: -569.31,
        s_standard_j: 26.95,
    },
    ThermochemData {
        formula: "Na2O(s)",
        name: "Sodium oxide",
        delta_hf_kj: -414.22,
        delta_gf_kj: -375.46,
        s_standard_j: 75.06,
    },
    ThermochemData {
        formula: "K2O(s)",
        name: "Potassium oxide",
        delta_hf_kj: -361.50,
        delta_gf_kj: -322.09,
        // CRC Handbook / NBS (94.14); 102.00 was not a tabulated value.
        s_standard_j: 94.10,
    },
    ThermochemData {
        formula: "ZnO(s)",
        name: "Zinc oxide",
        delta_hf_kj: -350.46,
        delta_gf_kj: -320.48,
        s_standard_j: 43.64,
    },
    ThermochemData {
        formula: "CuO(s)",
        name: "Copper(II) oxide",
        delta_hf_kj: -157.30,
        delta_gf_kj: -129.70,
        s_standard_j: 42.63,
    },
    ThermochemData {
        formula: "TiO2(s)",
        name: "Titanium dioxide (rutile)",
        delta_hf_kj: -944.00,
        delta_gf_kj: -888.80,
        s_standard_j: 50.33,
    },
    ThermochemData {
        formula: "SiO2(s)",
        name: "Silicon dioxide (quartz)",
        delta_hf_kj: -910.70,
        delta_gf_kj: -856.30,
        s_standard_j: 41.46,
    },
    ThermochemData {
        formula: "MgCO3(s)",
        name: "Magnesium carbonate",
        delta_hf_kj: -1095.80,
        delta_gf_kj: -1012.10,
        s_standard_j: 65.70,
    },
    ThermochemData {
        formula: "BaCO3(s)",
        name: "Barium carbonate",
        delta_hf_kj: -1216.30,
        delta_gf_kj: -1137.60,
        s_standard_j: 112.10,
    },
    ThermochemData {
        formula: "AgCl(s)",
        name: "Silver chloride",
        delta_hf_kj: -127.01,
        delta_gf_kj: -109.79,
        s_standard_j: 96.25,
    },
    ThermochemData {
        formula: "BaSO4(s)",
        name: "Barium sulfate",
        delta_hf_kj: -1473.20,
        delta_gf_kj: -1362.20,
        s_standard_j: 132.20,
    },
    ThermochemData {
        formula: "KCl(s)",
        name: "Potassium chloride",
        delta_hf_kj: -436.75,
        delta_gf_kj: -409.14,
        s_standard_j: 82.59,
    },
    ThermochemData {
        formula: "NaHCO3(s)",
        name: "Sodium bicarbonate",
        delta_hf_kj: -950.81,
        delta_gf_kj: -851.00,
        s_standard_j: 101.70,
    },
    ThermochemData {
        formula: "Na2CO3(s)",
        name: "Sodium carbonate",
        delta_hf_kj: -1130.70,
        delta_gf_kj: -1044.40,
        s_standard_j: 135.00,
    },
    ThermochemData {
        formula: "KNO3(s)",
        name: "Potassium nitrate",
        delta_hf_kj: -494.63,
        delta_gf_kj: -394.86,
        s_standard_j: 133.10,
    },
    ThermochemData {
        formula: "NH4Cl(s)",
        name: "Ammonium chloride",
        delta_hf_kj: -314.43,
        delta_gf_kj: -202.87,
        s_standard_j: 94.60,
    },
    ThermochemData {
        formula: "NH4NO3(s)",
        name: "Ammonium nitrate",
        delta_hf_kj: -365.56,
        delta_gf_kj: -183.87,
        s_standard_j: 151.08,
    },
    ThermochemData {
        formula: "HF(g)",
        name: "Hydrogen fluoride",
        delta_hf_kj: -273.30,
        delta_gf_kj: -275.40,
        s_standard_j: 173.78,
    },
    ThermochemData {
        formula: "HBr(g)",
        name: "Hydrogen bromide",
        delta_hf_kj: -36.29,
        delta_gf_kj: -53.45,
        s_standard_j: 198.70,
    },
    ThermochemData {
        formula: "HI(g)",
        name: "Hydrogen iodide",
        delta_hf_kj: 26.48,
        delta_gf_kj: 1.70,
        s_standard_j: 206.59,
    },
    ThermochemData {
        formula: "H2S(g)",
        name: "Hydrogen sulfide",
        delta_hf_kj: -20.60,
        delta_gf_kj: -33.40,
        s_standard_j: 205.79,
    },
    ThermochemData {
        formula: "PCl3(g)",
        name: "Phosphorus trichloride",
        delta_hf_kj: -287.00,
        delta_gf_kj: -267.80,
        s_standard_j: 311.78,
    },
    ThermochemData {
        formula: "PCl5(g)",
        name: "Phosphorus pentachloride",
        delta_hf_kj: -374.90,
        delta_gf_kj: -305.00,
        s_standard_j: 364.58,
    },
    // ── Organic compounds (additional) ────────────────────────────
    ThermochemData {
        formula: "CH3OH(l)",
        name: "Methanol",
        delta_hf_kj: -238.40,
        delta_gf_kj: -166.27,
        s_standard_j: 126.80,
    },
    ThermochemData {
        formula: "CH3COOH(l)",
        name: "Acetic acid",
        delta_hf_kj: -484.50,
        delta_gf_kj: -389.90,
        s_standard_j: 159.80,
    },
    ThermochemData {
        formula: "C6H6(l)",
        name: "Benzene",
        delta_hf_kj: 49.04,
        delta_gf_kj: 124.50,
        s_standard_j: 173.40,
    },
    ThermochemData {
        formula: "C3H8(g)",
        name: "Propane",
        delta_hf_kj: -103.85,
        delta_gf_kj: -23.40,
        s_standard_j: 270.02,
    },
    // n-Butane: all three values from the NBS tables (Wagman et al. 1982),
    // which are mutually consistent (ΔG = ΔH − TΔS within 0.05 kJ/mol).
    ThermochemData {
        formula: "C4H10(g)",
        name: "Butane",
        delta_hf_kj: -126.15,
        delta_gf_kj: -17.03,
        s_standard_j: 310.23,
    },
    ThermochemData {
        formula: "C8H18(l)",
        name: "Octane",
        delta_hf_kj: -250.10,
        delta_gf_kj: 6.40,
        s_standard_j: 361.20,
    },
    ThermochemData {
        formula: "HCHO(g)",
        name: "Formaldehyde",
        delta_hf_kj: -108.57,
        delta_gf_kj: -102.53,
        s_standard_j: 218.77,
    },
];

/// Look up thermochemical data by formula (e.g., "H2O(l)", "CO2(g)").
#[must_use]
#[inline]
pub fn lookup_thermochem(formula: &str) -> Option<&'static ThermochemData> {
    THERMOCHEM_DATA.iter().find(|d| d.formula == formula)
}

// ── Reaction thermochemistry from tables ─────────────────────────────

/// Standard reaction enthalpy from formation enthalpies at 298.15 K:
/// ΔH°_rxn = Σ n_i·ΔH_f°(products) − Σ n_j·ΔH_f°(reactants)
///
/// Each entry is (formula, stoichiometric coefficient).
///
/// Returns ΔH°_rxn in **kJ/mol** (per mole of reaction as written).
/// [`vant_hoff_k`] and [`adiabatic_flame_temperature`] take J/mol:
/// multiply this result by 1000 before passing it to them.
///
/// # Errors
///
/// Returns error if any formula is not found in the thermochemical data.
pub fn reaction_enthalpy(products: &[(&str, f64)], reactants: &[(&str, f64)]) -> Result<f64> {
    let sum_products = sum_property(products, |d| d.delta_hf_kj)?;
    let sum_reactants = sum_property(reactants, |d| d.delta_hf_kj)?;
    Ok(sum_products - sum_reactants)
}

/// Standard reaction Gibbs energy at 298.15 K:
/// ΔG°_rxn = Σ n_i·ΔG_f°(products) − Σ n_j·ΔG_f°(reactants)
///
/// Each entry is (formula, stoichiometric coefficient).
///
/// Returns ΔG°_rxn in **kJ/mol** (per mole of reaction as written).
/// Multiply by 1000 for functions that take J/mol.
///
/// # Errors
///
/// Returns error if any formula is not found in the thermochemical data.
pub fn reaction_gibbs_energy(products: &[(&str, f64)], reactants: &[(&str, f64)]) -> Result<f64> {
    let sum_products = sum_property(products, |d| d.delta_gf_kj)?;
    let sum_reactants = sum_property(reactants, |d| d.delta_gf_kj)?;
    Ok(sum_products - sum_reactants)
}

/// Standard reaction entropy change at 298.15 K:
/// ΔS°_rxn = Σ n_i·S°(products) − Σ n_j·S°(reactants)
///
/// Each entry is (formula, stoichiometric coefficient).
///
/// Returns ΔS°_rxn in **J/(mol·K)** (per mole of reaction as written) —
/// note J, not kJ: divide by 1000 before combining with
/// [`reaction_enthalpy`] in ΔG = ΔH − TΔS.
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
/// - `k_ref`: equilibrium constant at `t_ref_k` (dimensionless, > 0)
/// - `delta_h_j`: standard reaction enthalpy in **J/mol** (not kJ!) —
///   multiply a [`reaction_enthalpy`] result by 1000
/// - `t_ref_k`: reference temperature (K)
/// - `t_new_k`: target temperature (K)
///
/// Returns K at `t_new_k` (dimensionless). May overflow to `+∞` for
/// extreme ΔH°.
///
/// # Errors
///
/// Returns error if temperatures or K_ref are not positive and finite, or if
/// ΔH° is not finite.
pub fn vant_hoff_k(k_ref: f64, delta_h_j: f64, t_ref_k: f64, t_new_k: f64) -> Result<f64> {
    if !is_positive(t_ref_k) || !is_positive(t_new_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
        ));
    }
    if !is_positive(k_ref) {
        return Err(KimiyaError::InvalidInput(
            "reference K must be positive".into(),
        ));
    }
    if !is_finite(delta_h_j) {
        return Err(KimiyaError::InvalidInput(
            "reaction enthalpy must be finite".into(),
        ));
    }
    let exponent = -(delta_h_j / GAS_CONSTANT) * (1.0 / t_new_k - 1.0 / t_ref_k);
    Ok(k_ref * exponent.exp())
}

// ── Heat capacity integration ────────────────────────────────────────

/// Shomate equation coefficients for one temperature segment.
///
/// Cp(T) = A + B·t + C·t² + D·t³ + E/t²  in J/(mol·K), where t = T(K) / 1000.
///
/// The fit is valid only for `t_min ≤ T ≤ t_max`.
#[derive(Debug, Clone, Serialize)]
pub struct Shomate {
    /// Coefficient A, J/(mol·K).
    pub a: f64,
    /// Coefficient B (multiplies t), J/(mol·K) per kK.
    pub b: f64,
    /// Coefficient C (multiplies t²), J/(mol·K) per kK².
    pub c: f64,
    /// Coefficient D (multiplies t³), J/(mol·K) per kK³.
    pub d: f64,
    /// Coefficient E (multiplies 1/t²), J/(mol·K)·kK².
    pub e: f64,
    /// Lower bound of the fit's valid temperature range (K).
    pub t_min: f64,
    /// Upper bound of the fit's valid temperature range (K).
    pub t_max: f64,
}

impl Shomate {
    /// Evaluate Cp at temperature T (K) in J/(mol·K).
    ///
    /// Does not check `t_min`/`t_max`: outside that range the polynomial is
    /// an extrapolation. Use [`lookup_shomate_at`] to pick the right segment.
    #[must_use]
    #[inline]
    pub fn cp(&self, temperature_k: f64) -> f64 {
        let t = temperature_k / 1000.0;
        self.a + self.b * t + self.c * t * t + self.d * t * t * t + self.e / (t * t)
    }

    /// Antiderivative of [`Shomate::cp`] with respect to T, in J/mol:
    ///
    /// ∫Cp dT = 1000·(A·t + B·t²/2 + C·t³/3 + D·t⁴/4 − E/t),  t = T(K)/1000
    ///
    /// The integration constant (NIST's F − H terms) is omitted, so the value
    /// itself is not H(T) − H(298.15). Only differences are meaningful:
    /// within one segment, H(T₂) − H(T₁) =
    /// `enthalpy_antiderivative(T₂) − enthalpy_antiderivative(T₁)`.
    /// Use [`enthalpy_change_cp`] for ranges that span several segments.
    #[must_use]
    #[inline]
    pub fn enthalpy_antiderivative(&self, temperature_k: f64) -> f64 {
        let t = temperature_k / 1000.0;
        1000.0
            * (t * (self.a + t * (self.b / 2.0 + t * (self.c / 3.0 + t * (self.d / 4.0))))
                - self.e / t)
    }

    /// `t_min ≤ T ≤ t_max` (false for NaN).
    #[inline]
    #[must_use]
    fn covers(&self, temperature_k: f64) -> bool {
        temperature_k >= self.t_min && temperature_k <= self.t_max
    }
}

/// Built-in Shomate coefficients for common gases.
///
/// Coefficients give Cp in J/(mol·K) with t = T(K)/1000. Source is the NIST
/// Chemistry WebBook unless an entry says otherwise. Valid ranges differ per
/// entry: read each entry's `t_min`/`t_max`.
///
/// Species that need several fits to cover their range (N₂, O₂, CO₂, H₂O)
/// appear once per temperature segment under the same formula string. The
/// segments of one species are adjacent, in increasing temperature order,
/// and contiguous (each `t_max` equals the next segment's `t_min`). Use
/// [`lookup_shomate_at`] to select the segment for a temperature;
/// [`lookup_shomate`] returns the lowest-temperature segment.
pub static SHOMATE_DATA: &[(&str, Shomate)] = &[
    // NIST WebBook gives no gas-phase set below 500 K. Its 500–1700 K set is
    // used down to 298 K as a documented extrapolation so flame calculations
    // can start at 298.15 K: Cp(298.15) = 33.59 and Cp(400) = 34.26 J/(mol·K),
    // matching JANAF (33.590, 34.262) to < 0.1 %.
    (
        "H2O(g)",
        Shomate {
            a: 30.09200,
            b: 6.832514,
            c: 6.793435,
            d: -2.534480,
            e: 0.082139,
            t_min: 298.0,
            t_max: 1700.0,
        },
    ),
    // NIST WebBook Shomate, 1700–6000 K.
    (
        "H2O(g)",
        Shomate {
            a: 41.96426,
            b: 8.622053,
            c: -1.499780,
            d: 0.098119,
            e: -11.15764,
            t_min: 1700.0,
            t_max: 6000.0,
        },
    ),
    // NIST WebBook Shomate, 298–1200 K.
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
    // NIST WebBook Shomate, 1200–6000 K.
    (
        "CO2(g)",
        Shomate {
            a: 58.16639,
            b: 2.720074,
            c: -0.492289,
            d: 0.038844,
            e: -6.447293,
            t_min: 1200.0,
            t_max: 6000.0,
        },
    ),
    // NIST WebBook Shomate, 100–500 K.
    (
        "N2(g)",
        Shomate {
            a: 28.98641,
            b: 1.853978,
            c: -9.647459,
            d: 16.63537,
            e: 0.000117,
            t_min: 100.0,
            t_max: 500.0,
        },
    ),
    // NIST WebBook Shomate, 500–2000 K.
    (
        "N2(g)",
        Shomate {
            a: 19.50583,
            b: 19.88705,
            c: -8.598535,
            d: 1.369784,
            e: 0.527601,
            t_min: 500.0,
            t_max: 2000.0,
        },
    ),
    // NIST WebBook Shomate, 2000–6000 K.
    (
        "N2(g)",
        Shomate {
            a: 35.51872,
            b: 1.128728,
            c: -0.196103,
            d: 0.014662,
            e: -4.553760,
            t_min: 2000.0,
            t_max: 6000.0,
        },
    ),
    // NIST WebBook Shomate, 100–700 K.
    (
        "O2(g)",
        Shomate {
            a: 31.32234,
            b: -20.23531,
            c: 57.86644,
            d: -36.50624,
            e: -0.007374,
            t_min: 100.0,
            t_max: 700.0,
        },
    ),
    // NIST WebBook Shomate, 700–2000 K.
    (
        "O2(g)",
        Shomate {
            a: 30.03235,
            b: 8.772972,
            c: -3.988133,
            d: 0.788313,
            e: -0.741599,
            t_min: 700.0,
            t_max: 2000.0,
        },
    ),
    // NIST WebBook Shomate, 2000–6000 K.
    (
        "O2(g)",
        Shomate {
            a: 20.91111,
            b: 10.72071,
            c: -2.020498,
            d: 0.146449,
            e: 9.245722,
            t_min: 2000.0,
            t_max: 6000.0,
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
    // NIST WebBook Shomate, 1300–6000 K
    (
        "CH4(g)",
        Shomate {
            a: 85.81217,
            b: 11.26467,
            c: -2.114146,
            d: 0.13819,
            e: -26.42221,
            t_min: 1300.0,
            t_max: 6000.0,
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
    // NIST WebBook Shomate, 1300–6000 K
    (
        "CO(g)",
        Shomate {
            a: 35.1507,
            b: 1.300095,
            c: -0.205921,
            d: 0.01355,
            e: -3.28278,
            t_min: 1300.0,
            t_max: 6000.0,
        },
    ),
    (
        "H2(g)",
        Shomate {
            a: 33.066178,
            b: -11.363417,
            c: 11.432816,
            d: -2.772874,
            e: -0.158558,
            t_min: 298.0,
            t_max: 1000.0,
        },
    ),
    // NIST WebBook Shomate, 1000–2500 K
    (
        "H2(g)",
        Shomate {
            a: 18.563083,
            b: 12.257357,
            c: -2.859786,
            d: 0.268238,
            e: 1.97799,
            t_min: 1000.0,
            t_max: 2500.0,
        },
    ),
    // NIST WebBook Shomate, 2500–6000 K
    (
        "H2(g)",
        Shomate {
            a: 43.41356,
            b: -4.293079,
            c: 1.272428,
            d: -0.096876,
            e: -20.533862,
            t_min: 2500.0,
            t_max: 6000.0,
        },
    ),
    (
        "NH3(g)",
        Shomate {
            a: 19.99563,
            b: 49.77119,
            c: -15.37599,
            d: 1.921168,
            e: 0.189174,
            t_min: 298.0,
            t_max: 1400.0,
        },
    ),
    (
        "SO2(g)",
        Shomate {
            a: 21.43049,
            b: 74.35094,
            c: -57.75217,
            d: 16.35534,
            e: 0.086731,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "HCl(g)",
        Shomate {
            a: 32.12392,
            b: -13.45805,
            c: 19.86852,
            d: -6.853936,
            e: -0.049672,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "NO(g)",
        Shomate {
            a: 23.83491,
            b: 12.58878,
            c: -1.139011,
            d: -1.497459,
            e: 0.214194,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "NO2(g)",
        Shomate {
            a: 16.10857,
            b: 75.89525,
            c: -54.38740,
            d: 14.30777,
            e: 0.239423,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    // Approximate: not a NIST WebBook Shomate set (NIST publishes none for
    // ethane) and the source is unknown. Known error: Cp(298.15) = 54.3 vs
    // reference 52.5 J/(mol·K) (+3.4 %). Use with caution.
    (
        "C2H6(g)",
        Shomate {
            a: 6.16331,
            b: 173.5507,
            c: -64.14486,
            d: 7.523975,
            e: 0.168137,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    // NIST WebBook Shomate, 298–1200 K (Cp(298.15) = 42.85 J/(mol·K)).
    (
        "C2H4(g)",
        Shomate {
            a: -6.387880,
            b: 184.4019,
            c: -112.9718,
            d: 28.49593,
            e: 0.315540,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "C2H2(g)",
        Shomate {
            a: 40.68697,
            b: 40.73279,
            c: -16.17840,
            d: 3.669741,
            e: -0.658411,
            t_min: 298.0,
            t_max: 1100.0,
        },
    ),
    // ── Additional Shomate data ──────────────────────────────────
    // Approximate: not a NIST WebBook Shomate set (NIST publishes none for
    // propane) and the source is unknown. Known error: Cp(298.15) = 79.8 vs
    // reference 73.6 J/(mol·K) (+8.4 %). Use with caution.
    (
        "C3H8(g)",
        Shomate {
            a: -4.224394,
            b: 306.0191,
            c: -158.5991,
            d: 32.13038,
            e: 0.534834,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "H2S(g)",
        Shomate {
            a: 26.88412,
            b: 18.67809,
            c: 3.434203,
            d: -3.378702,
            e: 0.135882,
            t_min: 298.0,
            t_max: 1400.0,
        },
    ),
    (
        "Cl2(g)",
        Shomate {
            a: 33.0506,
            b: 12.2294,
            c: -12.0651,
            d: 4.38533,
            e: -0.159494,
            t_min: 298.0,
            t_max: 1000.0,
        },
    ),
    // NIST WebBook Shomate, 332.5–3400 K.
    (
        "Br2(g)",
        Shomate {
            a: 38.52723,
            b: -1.976835,
            c: 1.526107,
            d: -0.198398,
            e: -0.185815,
            t_min: 332.5,
            t_max: 3400.0,
        },
    ),
    // KNOWN DEFECT — not the NIST set and of unknown source. This fit gives
    // Cp(298.15) = 28.77 J/(mol·K) (JANAF: 29.14), falls to ~27.8 near 700 K,
    // and stays below the 7R/2 = 29.10 floor of a diatomic gas across its
    // whole range, which is physically impossible. Do not rely on HF heat
    // capacities until this is replaced with a sourced fit.
    (
        "HF(g)",
        Shomate {
            a: 29.14000,
            b: -6.56728,
            c: 8.27622,
            d: -2.91489,
            e: 0.08296,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "HBr(g)",
        Shomate {
            a: 29.14100,
            b: -3.42838,
            c: 8.16490,
            d: -3.05563,
            e: 0.02633,
            t_min: 298.0,
            t_max: 1200.0,
        },
    ),
    (
        "Ar(g)",
        Shomate {
            a: 20.786,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 0.0,
            t_min: 298.0,
            t_max: 6000.0,
        },
    ),
    (
        "He(g)",
        Shomate {
            a: 20.786,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 0.0,
            t_min: 298.0,
            t_max: 6000.0,
        },
    ),
    (
        "Ne(g)",
        Shomate {
            a: 20.786,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 0.0,
            t_min: 298.0,
            t_max: 6000.0,
        },
    ),
];

/// Look up Shomate coefficients by formula.
///
/// Returns the first (lowest-temperature) segment. For species with several
/// segments (see [`SHOMATE_DATA`]) use [`lookup_shomate_at`] instead.
#[must_use]
#[inline]
pub fn lookup_shomate(formula: &str) -> Option<&'static Shomate> {
    SHOMATE_DATA
        .iter()
        .find(|(f, _)| *f == formula)
        .map(|(_, s)| s)
}

/// Look up the Shomate segment for `formula` whose range covers
/// `temperature_k` (K), i.e. `t_min ≤ T ≤ t_max`.
///
/// At a boundary shared by two segments the lower-temperature segment is
/// returned (both are valid there). Returns `None` if the formula is unknown
/// or no segment covers the temperature (including NaN).
#[must_use]
#[inline]
pub fn lookup_shomate_at(formula: &str, temperature_k: f64) -> Option<&'static Shomate> {
    SHOMATE_DATA
        .iter()
        .find(|(f, s)| *f == formula && s.covers(temperature_k))
        .map(|(_, s)| s)
}

/// The contiguous run of Shomate segments for one species.
struct Segments {
    segs: &'static [(&'static str, Shomate)],
    /// Lowest temperature covered (K).
    t_min: f64,
    /// Highest temperature covered (K).
    t_max: f64,
}

impl Segments {
    /// Resolve the segments of `formula` (adjacent entries in [`SHOMATE_DATA`]).
    #[must_use]
    fn lookup(formula: &str) -> Option<Self> {
        let start = SHOMATE_DATA.iter().position(|(f, _)| *f == formula)?;
        let rest = SHOMATE_DATA.get(start..)?;
        let len = rest.iter().take_while(|(f, _)| *f == formula).count();
        let segs = rest.get(..len)?;
        let (_, first) = segs.first()?;
        let (_, last) = segs.last()?;
        Some(Self {
            segs,
            t_min: first.t_min,
            t_max: last.t_max,
        })
    }

    /// Cp(T) in J/(mol·K) from the segment covering T. Below the first or
    /// above the last segment, that end segment is extrapolated.
    #[inline]
    #[must_use]
    fn cp(&self, temperature_k: f64) -> f64 {
        self.segs
            .iter()
            .find(|(_, s)| temperature_k <= s.t_max)
            .or_else(|| self.segs.last())
            .map_or(f64::NAN, |(_, s)| s.cp(temperature_k))
    }

    /// ∫[T₁→T₂] Cp dT in J/mol, summed analytically segment by segment.
    /// Segment boundaries are the interior `t_max` values; the first and
    /// last segments are extrapolated beyond the covered range.
    #[inline]
    #[must_use]
    fn enthalpy_change(&self, t1_k: f64, t2_k: f64) -> f64 {
        let (lo, hi, sign) = if t1_k <= t2_k {
            (t1_k, t2_k, 1.0)
        } else {
            (t2_k, t1_k, -1.0)
        };
        let n = self.segs.len();
        let mut total = 0.0;
        let mut lower = f64::NEG_INFINITY;
        for (i, (_, s)) in self.segs.iter().enumerate() {
            let upper = if i + 1 == n { f64::INFINITY } else { s.t_max };
            let a = lo.max(lower);
            let b = hi.min(upper);
            if b > a {
                total += s.enthalpy_antiderivative(b) - s.enthalpy_antiderivative(a);
            }
            if upper >= hi {
                break;
            }
            lower = upper;
        }
        sign * total
    }
}

/// Enthalpy change of one species from T₁ to T₂ via Cp(T) integration:
/// ΔH = ∫[T₁→T₂] Cp(T) dT
///
/// Integrates the Shomate polynomial analytically, piecewise across the
/// species' temperature segments (see [`SHOMATE_DATA`]).
///
/// Returns the sensible-enthalpy change in **J/mol** (not kJ); negative when
/// T₂ < T₁.
///
/// # Errors
///
/// Returns error if the formula has no Shomate data, if a temperature is not
/// positive and finite, or if any part of [T₁, T₂] lies outside the range
/// covered by the species' Shomate segments.
pub fn enthalpy_change_cp(formula: &str, t1_k: f64, t2_k: f64) -> Result<f64> {
    let segments = Segments::lookup(formula)
        .ok_or_else(|| KimiyaError::InvalidReaction(format!("no Shomate data for '{formula}'")))?;
    if !is_positive(t1_k) || !is_positive(t2_k) {
        return Err(KimiyaError::InvalidTemperature(
            "temperatures must be positive".into(),
        ));
    }
    let (lo, hi) = (t1_k.min(t2_k), t1_k.max(t2_k));
    if lo < segments.t_min || hi > segments.t_max {
        return Err(KimiyaError::InvalidTemperature(format!(
            "{lo}–{hi} K is outside the Shomate range {}–{} K for '{formula}'",
            segments.t_min, segments.t_max
        )));
    }
    Ok(segments.enthalpy_change(t1_k, t2_k))
}

/// Newton–Raphson stopping tolerance on the energy balance (J/mol).
const FLAME_NEWTON_TOL_J: f64 = 1e-3;
/// Bisection stopping tolerance (K, and J/mol on the residual).
const FLAME_BISECTION_TOL: f64 = 1e-6;

/// Adiabatic flame temperature (no dissociation, constant pressure).
///
/// Finds T where: ΔH_rxn + Σ nᵢ·∫\[T₀→T\] Cpᵢ(T) dT = 0
///
/// Sensible enthalpies are integrated analytically and piecewise across each
/// product's Shomate segments. Uses [`hisab::num::newton_raphson`] with
/// derivative f'(T) = −Σ nᵢ·Cpᵢ(T), falling back to [`hisab::num::bisection`]
/// on the bracket [T₀, T_max].
///
/// - `reaction_enthalpy_j`: ΔH_rxn in **J/mol** (not kJ; exothermic =
///   negative), evaluated at `t_initial_k`. [`reaction_enthalpy`] returns
///   kJ/mol at 298.15 K — multiply it by 1000, and correct it to T₀ if
///   T₀ ≠ 298.15 K.
/// - `products`: list of (formula, moles) for product gases, per mole of
///   reaction as written; every product needs Shomate data.
/// - `t_initial_k`: initial (reactant) temperature T₀ (K).
///
/// Returns adiabatic flame temperature in Kelvin. The search is restricted
/// to the range covered by every product's Shomate data, T_max being the
/// lowest upper bound among them.
///
/// # Examples
///
/// ```
/// use kimiya::thermochem::{adiabatic_flame_temperature, reaction_enthalpy};
///
/// // CH₄ + 2 O₂ → CO₂ + 2 H₂O(g), burned in air (7.52 mol N₂ per mol CH₄).
/// let dh_kj = reaction_enthalpy(
///     &[("CO2(g)", 1.0), ("H2O(g)", 2.0)],
///     &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
/// )?;
/// let t = adiabatic_flame_temperature(
///     dh_kj * 1000.0, // kJ/mol → J/mol
///     &[("CO2(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 7.52)],
///     298.15,
/// )?;
/// assert!((t - 2327.0).abs() < 10.0);
/// # Ok::<(), kimiya::KimiyaError>(())
/// ```
///
/// # Errors
///
/// Returns error if `t_initial_k` is not positive and finite or lies outside
/// the products' common Shomate range, if ΔH_rxn is not finite, if
/// `products` is empty, if any product has no Shomate data or non-positive /
/// non-finite moles, if the root cannot be bracketed (ΔH_rxn ≥ 0, or the
/// flame temperature exceeds the products' common Shomate range), or if
/// root finding fails.
pub fn adiabatic_flame_temperature(
    reaction_enthalpy_j: f64,
    products: &[(&str, f64)],
    t_initial_k: f64,
) -> Result<f64> {
    if !is_positive(t_initial_k) {
        return Err(KimiyaError::InvalidTemperature(
            "initial temperature must be positive".into(),
        ));
    }
    if !is_finite(reaction_enthalpy_j) {
        return Err(KimiyaError::InvalidInput(
            "reaction enthalpy must be finite".into(),
        ));
    }
    if products.is_empty() {
        return Err(KimiyaError::InvalidInput(
            "at least one product is required".into(),
        ));
    }

    // Resolve each product's Shomate segments once, and the temperature
    // range covered by all of them.
    let mut resolved = Vec::with_capacity(products.len());
    let mut range_min = f64::NEG_INFINITY;
    let mut range_max = f64::INFINITY;
    for &(formula, moles) in products {
        let segments = Segments::lookup(formula).ok_or_else(|| {
            KimiyaError::InvalidReaction(format!("no Shomate data for '{formula}'"))
        })?;
        if !is_positive(moles) {
            return Err(KimiyaError::InvalidInput(format!(
                "moles of '{formula}' must be positive and finite"
            )));
        }
        range_min = range_min.max(segments.t_min);
        range_max = range_max.min(segments.t_max);
        resolved.push((segments, moles));
    }
    if t_initial_k < range_min || t_initial_k >= range_max {
        return Err(KimiyaError::InvalidTemperature(format!(
            "initial temperature {t_initial_k} K is outside the Shomate range \
             {range_min}–{range_max} K shared by the products"
        )));
    }

    // f(T) = −ΔH_rxn − Σ nᵢ × ∫[T₀→T] Cpᵢ dT
    // f'(T) = −Σ nᵢ × Cpᵢ(T)
    let f = |t: f64| -> f64 {
        let heat_absorbed: f64 = resolved
            .iter()
            .map(|(s, moles)| moles * s.enthalpy_change(t_initial_k, t))
            .sum();
        -reaction_enthalpy_j - heat_absorbed
    };
    let df = |t: f64| -> f64 {
        let total_cp: f64 = resolved.iter().map(|(s, moles)| moles * s.cp(t)).sum();
        -total_cp
    };

    let t_lo = t_initial_k;
    let t_hi = range_max;
    // f(T₀) = −ΔH_rxn and f is decreasing, so a sign change on [T₀, T_max]
    // means the reaction is exothermic and the root lies in the data range.
    // Written as a product test so NaN fails it.
    let bracketed = f(t_lo) * f(t_hi) < 0.0;
    if !bracketed {
        return Err(KimiyaError::ComputationError(format!(
            "cannot bracket flame temperature in {t_lo}–{t_hi} K — reaction enthalpy \
             must be negative (J/mol) and the flame temperature within the \
             products' Shomate range"
        )));
    }

    // Initial guess: assume average Cp ~ 40 J/(mol·K) per mol of product.
    let total_moles: f64 = resolved.iter().map(|(_, moles)| moles).sum();
    let t_guess = (t_initial_k - reaction_enthalpy_j / (total_moles * 40.0))
        .max(t_lo)
        .min(t_hi);

    // Newton–Raphson first (fast convergence); accept only roots in range.
    if let Ok(t) = hisab::num::newton_raphson(f, df, t_guess, FLAME_NEWTON_TOL_J, 50)
        && (t_lo..=t_hi).contains(&t)
    {
        return Ok(t);
    }

    // Fallback: bisection on the verified bracket.
    hisab::num::bisection(f, t_lo, t_hi, FLAME_BISECTION_TOL, 200)
        .map_err(|e| KimiyaError::ComputationError(format!("flame temperature solver failed: {e}")))
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
        // CO2 Cp at 298.15 K = 37.12 J/(mol·K) (NIST)
        let s = lookup_shomate("CO2(g)").unwrap();
        let cp = s.cp(298.15);
        assert!(
            (cp - 37.12).abs() < 0.1,
            "CO2 Cp at 298K should be ~37.12 J/(mol·K), got {cp}"
        );
    }

    #[test]
    fn shomate_n2_cp_at_298() {
        // N2 Cp at 298.15 K = 29.12 J/(mol·K) (NIST)
        let s = lookup_shomate("N2(g)").unwrap();
        let cp = s.cp(298.15);
        assert!(
            (cp - 29.12).abs() < 0.1,
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

    #[test]
    fn expanded_shomate_h2_cp_at_298() {
        let s = lookup_shomate("H2(g)").unwrap();
        let cp = s.cp(298.15);
        // H2 Cp at 298.15 K = 28.84 J/(mol·K) (NIST)
        assert!(
            (cp - 28.84).abs() < 0.1,
            "H2 Cp at 298K should be ~28.8, got {cp}"
        );
    }

    #[test]
    fn shomate_data_count() {
        // 24 species; N2 and O2 have 3 segments, CO2 and H2O have 2.
        assert_eq!(SHOMATE_DATA.len(), 34, "should have 34 Shomate entries");
    }

    // ── Adiabatic flame temperature ──────────────────────────────────

    #[test]
    fn adiabatic_flame_h2_air() {
        // H₂ + ½O₂ + 1.88 N₂ → H₂O(g) + 1.88 N₂ (stoichiometric in air)
        // ΔH = -241.82 kJ/mol = -241820 J/mol
        // N₂ diluent absorbs heat, bringing flame temp to a solvable range
        // No-dissociation literature value ≈ 2520 K; this data gives 2524.7 K.
        let t =
            adiabatic_flame_temperature(-241_820.0, &[("H2O(g)", 1.0), ("N2(g)", 1.88)], 298.15)
                .unwrap();
        assert!((t - 2520.0).abs() < 30.0, "H₂/air vs literature: {t}");
        assert!((t - 2524.7).abs() < 10.0, "H₂/air flame: {t}");
    }

    #[test]
    fn adiabatic_flame_ch4_air() {
        // CH₄ + 2O₂ → CO₂ + 2H₂O(g), ΔH = -802.3 kJ/mol = -802300 J/mol
        // In air with excess N₂ (~7.52 mol N₂ per mol CH₄)
        // No-dissociation literature value ≈ 2325 K; this data gives 2327.2 K.
        let t = adiabatic_flame_temperature(
            -802_300.0,
            &[("CO2(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 7.52)],
            298.15,
        )
        .unwrap();
        assert!((t - 2325.0).abs() < 30.0, "CH₄/air vs literature: {t}");
        assert!((t - 2327.2).abs() < 10.0, "CH₄/air flame: {t}");
    }

    #[test]
    fn adiabatic_flame_unknown_product_is_error() {
        assert!(adiabatic_flame_temperature(-100_000.0, &[("XeF6(g)", 1.0)], 298.15).is_err());
    }

    #[test]
    fn adiabatic_flame_zero_temp_is_error() {
        assert!(adiabatic_flame_temperature(-100_000.0, &[("CO2(g)", 1.0)], 0.0).is_err());
    }

    // ── Multi-segment Shomate data (audit fixes) ─────────────────────

    /// Segment-aware Cp, J/(mol·K).
    fn cp_at(formula: &str, t: f64) -> f64 {
        lookup_shomate_at(formula, t)
            .unwrap_or_else(|| panic!("no {formula} segment at {t} K"))
            .cp(t)
    }

    #[test]
    fn cp_matches_reference_at_298_1000_2000() {
        // JANAF / NIST WebBook Cp, J/(mol·K).
        let cases = [
            ("N2(g)", 298.15, 29.12),
            ("N2(g)", 1000.0, 32.70),
            ("N2(g)", 2000.0, 35.97),
            ("O2(g)", 298.15, 29.38),
            ("O2(g)", 1000.0, 34.87),
            ("O2(g)", 2000.0, 37.78),
            ("CO2(g)", 298.15, 37.12),
            ("CO2(g)", 1000.0, 54.31),
            ("CO2(g)", 2000.0, 60.35),
            ("H2O(g)", 298.15, 33.59),
            ("H2O(g)", 1000.0, 41.27),
            ("H2O(g)", 2000.0, 51.18),
        ];
        for (formula, t, expected) in cases {
            let cp = cp_at(formula, t);
            assert!(
                (cp - expected).abs() < 0.1,
                "{formula} Cp({t} K) = {cp}, expected {expected}"
            );
        }
    }

    #[test]
    fn shomate_segments_are_adjacent_ordered_and_contiguous() {
        for (i, (formula, s)) in SHOMATE_DATA.iter().enumerate() {
            // All segments of a formula form one adjacent run.
            let first = SHOMATE_DATA.iter().position(|(f, _)| f == formula).unwrap();
            let run = SHOMATE_DATA[first..]
                .iter()
                .take_while(|(f, _)| f == formula)
                .count();
            assert!(
                i < first + run,
                "{formula}: segment at index {i} is not adjacent to its other segments"
            );
            if let Some((next_formula, next)) = SHOMATE_DATA.get(i + 1)
                && next_formula == formula
            {
                assert!(
                    (s.t_max - next.t_min).abs() < f64::EPSILON,
                    "{formula}: gap/overlap {}–{} K",
                    s.t_max,
                    next.t_min
                );
            }
        }
    }

    #[test]
    fn shomate_segment_continuity() {
        let mut boundaries = 0;
        for pair in SHOMATE_DATA.windows(2) {
            let [(fa, a), (fb, b)] = pair else { continue };
            if fa != fb {
                continue;
            }
            boundaries += 1;
            let t = a.t_max;
            let dcp = (a.cp(t) - b.cp(t)).abs();
            assert!(dcp < 0.5, "{fa}: Cp jump {dcp} J/(mol·K) at {t} K");
        }
        assert_eq!(
            boundaries, 10,
            "N2 ×2, O2 ×2, CO2 ×1, H2O ×1, CH4 ×1, CO ×1, H2 ×2"
        );
    }

    #[test]
    fn lookup_shomate_at_selects_covering_segment() {
        assert_eq!(lookup_shomate_at("N2(g)", 300.0).unwrap().t_max, 500.0);
        assert_eq!(lookup_shomate_at("N2(g)", 1000.0).unwrap().t_min, 500.0);
        assert_eq!(lookup_shomate_at("N2(g)", 3000.0).unwrap().t_min, 2000.0);
        // Shared boundary → lower segment.
        assert_eq!(lookup_shomate_at("CO2(g)", 1200.0).unwrap().t_max, 1200.0);
        assert!(lookup_shomate_at("N2(g)", 50.0).is_none());
        assert!(lookup_shomate_at("N2(g)", 7000.0).is_none());
        assert!(lookup_shomate_at("N2(g)", f64::NAN).is_none());
        assert!(lookup_shomate_at("XeF6(g)", 500.0).is_none());
        // lookup_shomate keeps returning the first segment.
        assert_eq!(lookup_shomate("N2(g)").unwrap().t_min, 100.0);
    }

    #[test]
    fn enthalpy_antiderivative_differentiates_to_cp() {
        for (formula, s) in SHOMATE_DATA.iter() {
            let t = 0.5 * (s.t_min + s.t_max);
            let h = 1e-3;
            let numeric =
                (s.enthalpy_antiderivative(t + h) - s.enthalpy_antiderivative(t - h)) / (2.0 * h);
            assert!(
                (numeric - s.cp(t)).abs() < 1e-4,
                "{formula}: d/dT antiderivative {numeric} vs Cp {}",
                s.cp(t)
            );
        }
    }

    #[test]
    fn enthalpy_change_matches_janaf() {
        // JANAF H°(T) − H°(298.15), kJ/mol.
        let cases = [
            ("N2(g)", 1000.0, 21.463),
            ("N2(g)", 2000.0, 56.137),
            ("N2(g)", 3000.0, 92.715),
            ("O2(g)", 1000.0, 22.707),
            ("O2(g)", 2000.0, 59.199),
            ("CO2(g)", 1000.0, 33.397),
            ("CO2(g)", 2000.0, 91.439),
            ("CO2(g)", 3000.0, 152.852),
            ("H2O(g)", 1000.0, 26.000),
            ("H2O(g)", 2000.0, 72.790),
            ("H2O(g)", 3000.0, 126.549),
        ];
        for (formula, t, expected_kj) in cases {
            let dh = enthalpy_change_cp(formula, 298.15, t).unwrap() / 1000.0;
            assert!(
                (dh - expected_kj).abs() < 0.05,
                "{formula} H({t})−H(298.15) = {dh} kJ/mol, JANAF {expected_kj}"
            );
        }
    }

    #[test]
    fn enthalpy_change_across_segments_matches_numeric_integral() {
        // CO2 crosses 1200 K; N2 crosses 500 and 2000 K; O2 crosses 700 and 2000 K.
        for (formula, t1, t2) in [
            ("CO2(g)", 1000.0, 1500.0),
            ("N2(g)", 400.0, 2500.0),
            ("O2(g)", 298.15, 3000.0),
            ("H2O(g)", 298.15, 2500.0),
        ] {
            let analytic = enthalpy_change_cp(formula, t1, t2).unwrap();
            // Composite midpoint rule on the segment-aware Cp(T).
            let steps = 20_000;
            let h = (t2 - t1) / f64::from(steps);
            let numeric: f64 = (0..steps)
                .map(|i| h * cp_at(formula, t1 + (f64::from(i) + 0.5) * h))
                .sum();
            assert!(
                (analytic - numeric).abs() < 0.5,
                "{formula} {t1}→{t2}: analytic {analytic} vs numeric {numeric}"
            );
            let reverse = enthalpy_change_cp(formula, t2, t1).unwrap();
            assert!(
                (analytic + reverse).abs() < 1e-6,
                "{formula}: not antisymmetric"
            );
        }
    }

    #[test]
    fn enthalpy_change_outside_data_range_is_error() {
        let out_of_range = [
            ("CO2(g)", 250.0, 600.0),  // below 298 K
            ("NH3(g)", 300.0, 1500.0), // above 1400 K
            ("N2(g)", 300.0, 6500.0),  // above 6000 K
            ("Br2(g)", 298.15, 400.0), // below 332.5 K
        ];
        for (formula, t1, t2) in out_of_range {
            assert!(
                matches!(
                    enthalpy_change_cp(formula, t1, t2),
                    Err(KimiyaError::InvalidTemperature(_))
                ),
                "{formula} {t1}→{t2} should be out of range"
            );
        }
        assert!(enthalpy_change_cp("N2(g)", 100.0, 6000.0).is_ok());
    }

    #[test]
    fn enthalpy_change_nan_inf_is_error() {
        for (t1, t2) in [
            (f64::NAN, 600.0),
            (300.0, f64::NAN),
            (f64::INFINITY, 600.0),
            (300.0, f64::INFINITY),
            (-300.0, 600.0),
        ] {
            assert!(
                matches!(
                    enthalpy_change_cp("N2(g)", t1, t2),
                    Err(KimiyaError::InvalidTemperature(_))
                ),
                "{t1}→{t2} should be rejected"
            );
        }
    }

    #[test]
    fn c2h4_cp_matches_nist() {
        // NIST: Cp(C2H4, 298.15 K) = 42.90 J/(mol·K)
        let cp = lookup_shomate("C2H4(g)").unwrap().cp(298.15);
        assert!((cp - 42.90).abs() < 0.1, "C2H4 Cp(298.15) = {cp}");
    }

    #[test]
    fn br2_cp_at_lower_bound() {
        // JANAF Cp(Br2, g, 298.15 K) = 36.06 J/(mol·K) and Cp rises with T,
        // so Cp at the 332.5 K lower bound must be ≥ 36.06.
        let s = lookup_shomate("Br2(g)").unwrap();
        assert!((s.t_min - 332.5).abs() < f64::EPSILON);
        let cp = s.cp(332.5);
        assert!(cp > 36.06 && cp < 36.6, "Br2 Cp(332.5) = {cp}");
    }

    // ── Van't Hoff NaN guards ────────────────────────────────────────

    #[test]
    fn vant_hoff_nan_inf_is_error() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        assert!(vant_hoff_k(1.0, -50_000.0, nan, 400.0).is_err());
        assert!(vant_hoff_k(1.0, -50_000.0, 298.15, nan).is_err());
        assert!(vant_hoff_k(1.0, -50_000.0, inf, 400.0).is_err());
        assert!(vant_hoff_k(1.0, -50_000.0, 298.15, inf).is_err());
        assert!(vant_hoff_k(nan, -50_000.0, 298.15, 400.0).is_err());
        assert!(vant_hoff_k(inf, -50_000.0, 298.15, 400.0).is_err());
        assert!(vant_hoff_k(1.0, nan, 298.15, 400.0).is_err());
        assert!(vant_hoff_k(1.0, -inf, 298.15, 400.0).is_err());
    }

    // ── Adiabatic flame (audit fixes) ────────────────────────────────

    #[test]
    fn adiabatic_flame_satisfies_energy_balance() {
        let products = [("CO2(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 7.52)];
        let t = adiabatic_flame_temperature(-802_300.0, &products, 298.15).unwrap();
        let absorbed: f64 = products
            .iter()
            .map(|&(f, n)| n * enthalpy_change_cp(f, 298.15, t).unwrap())
            .sum();
        assert!((absorbed - 802_300.0).abs() < 1.0, "residual {absorbed}");
    }

    #[test]
    fn adiabatic_flame_nan_enthalpy_is_error() {
        // Previously returned Ok(≈299.5 K): NaN.signum() == NaN.signum() is false.
        let p = [("CO2(g)", 1.0), ("N2(g)", 3.76)];
        assert!(adiabatic_flame_temperature(f64::NAN, &p, 298.15).is_err());
        assert!(adiabatic_flame_temperature(f64::NEG_INFINITY, &p, 298.15).is_err());
    }

    #[test]
    fn adiabatic_flame_invalid_inputs_are_errors() {
        let p = [("CO2(g)", 1.0)];
        assert!(adiabatic_flame_temperature(-100_000.0, &p, f64::NAN).is_err());
        assert!(adiabatic_flame_temperature(-100_000.0, &p, f64::INFINITY).is_err());
        assert!(adiabatic_flame_temperature(-100_000.0, &[], 298.15).is_err());
        for moles in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                adiabatic_flame_temperature(-100_000.0, &[("CO2(g)", moles)], 298.15).is_err(),
                "moles = {moles} should be rejected"
            );
        }
        // Endothermic / zero ΔH: no flame temperature above T₀.
        assert!(adiabatic_flame_temperature(100_000.0, &p, 298.15).is_err());
        assert!(adiabatic_flame_temperature(0.0, &p, 298.15).is_err());
    }

    #[test]
    fn adiabatic_flame_root_outside_data_range_is_error() {
        // NH3(g) data end at 1400 K; the root here lies far above that.
        assert!(adiabatic_flame_temperature(-200_000.0, &[("NH3(g)", 1.0)], 298.15).is_err());
        // Mixed products: the common range ends at SO2's 1200 K.
        assert!(
            adiabatic_flame_temperature(-241_820.0, &[("H2O(g)", 1.0), ("SO2(g)", 1.0)], 298.15)
                .is_err()
        );
        // Below the products' data range (CO2 starts at 298 K).
        assert!(adiabatic_flame_temperature(-100_000.0, &[("CO2(g)", 1.0)], 250.0).is_err());
    }

    // ── Table data (audit fixes) ─────────────────────────────────────

    #[test]
    fn k2o_entropy_matches_crc() {
        let d = lookup_thermochem("K2O(s)").unwrap();
        assert!(
            (d.s_standard_j - 94.1).abs() < 0.05,
            "K2O S° = {}",
            d.s_standard_j
        );
    }

    #[test]
    fn butane_row_is_self_consistent() {
        // NBS: ΔHf −126.15, ΔGf −17.03, S° 310.23.
        let d = lookup_thermochem("C4H10(g)").unwrap();
        let elements = 4.0 * lookup_thermochem("C(graphite)").unwrap().s_standard_j
            + 5.0 * lookup_thermochem("H2(g)").unwrap().s_standard_j;
        let dg_calc = d.delta_hf_kj - 298.15 * (d.s_standard_j - elements) / 1000.0;
        assert!(
            (dg_calc - d.delta_gf_kj).abs() < 0.2,
            "butane: ΔG_calc = {dg_calc}, table {}",
            d.delta_gf_kj
        );
    }

    #[test]
    fn added_standard_state_elements() {
        for (formula, s) in [
            ("Zn(s)", 41.63),
            ("Ti(s)", 30.72),
            ("Ba(s)", 62.5),
            ("Ag(s)", 42.55),
            ("F2(g)", 202.79),
        ] {
            let d = lookup_thermochem(formula).unwrap();
            assert!(d.delta_hf_kj.abs() < f64::EPSILON && d.delta_gf_kj.abs() < f64::EPSILON);
            assert!((d.s_standard_j - s).abs() < 0.01, "{formula}");
        }
        // Formation reactions that needed them now resolve.
        let dg =
            reaction_gibbs_energy(&[("ZnO(s)", 1.0)], &[("Zn(s)", 1.0), ("O2(g)", 0.5)]).unwrap();
        assert!((dg - (-320.48)).abs() < 0.01);
    }

    /// Reference standard state of each element: (formula, atoms per unit).
    fn reference_state(symbol: &str) -> Option<(&'static str, f64)> {
        Some(match symbol {
            "H" => ("H2(g)", 2.0),
            "O" => ("O2(g)", 2.0),
            "N" => ("N2(g)", 2.0),
            "F" => ("F2(g)", 2.0),
            "Cl" => ("Cl2(g)", 2.0),
            "Br" => ("Br2(l)", 2.0),
            "I" => ("I2(s)", 2.0),
            "C" => ("C(graphite)", 1.0),
            "S" => ("S(rhombic)", 1.0),
            "P" => ("P(white)", 1.0),
            "Fe" => ("Fe(s)", 1.0),
            "Cu" => ("Cu(s)", 1.0),
            "Al" => ("Al(s)", 1.0),
            "Na" => ("Na(s)", 1.0),
            "K" => ("K(s)", 1.0),
            "Ca" => ("Ca(s)", 1.0),
            "Mg" => ("Mg(s)", 1.0),
            "Si" => ("Si(s)", 1.0),
            "Zn" => ("Zn(s)", 1.0),
            "Ti" => ("Ti(s)", 1.0),
            "Ba" => ("Ba(s)", 1.0),
            "Ag" => ("Ag(s)", 1.0),
            _ => return None,
        })
    }

    /// Parse a parenthesis-free formula such as "CH3COOH(l)" into
    /// (element symbol, count) pairs; the phase suffix is ignored.
    fn parse_simple_formula(formula: &str) -> Vec<(String, f64)> {
        let body = formula.split('(').next().unwrap();
        let mut out: Vec<(String, f64)> = Vec::new();
        let chars: Vec<char> = body.chars().collect();
        let mut i = 0;
        while i < chars.len() {
            let mut symbol = chars[i].to_string();
            i += 1;
            while i < chars.len() && chars[i].is_ascii_lowercase() {
                symbol.push(chars[i]);
                i += 1;
            }
            let mut digits = String::new();
            while i < chars.len() && chars[i].is_ascii_digit() {
                digits.push(chars[i]);
                i += 1;
            }
            let n: f64 = if digits.is_empty() {
                1.0
            } else {
                digits.parse().unwrap()
            };
            match out.iter_mut().find(|(s, _)| *s == symbol) {
                Some((_, count)) => *count += n,
                None => out.push((symbol, n)),
            }
        }
        out
    }

    /// Rows whose source values (NBS / CRC) are themselves inconsistent
    /// beyond 2 kJ/mol, with the observed ΔG_calc − ΔG_table. Pinned so a
    /// data change forces a re-check.
    const KNOWN_INCONSISTENT: &[(&str, f64)] = &[("Na2CO3(s)", -2.55)];

    #[test]
    fn table_gibbs_enthalpy_entropy_consistency() {
        // ΔG_f° ≈ ΔH_f° − 298.15·ΔS_f°, ΔS_f° = S°(compound) − Σ S°(elements).
        let mut checked = 0;
        for d in THERMOCHEM_DATA.iter() {
            let composition = parse_simple_formula(d.formula);
            // Skip reference-state elements themselves.
            if composition.len() == 1
                && reference_state(&composition[0].0).map(|(f, _)| f) == Some(d.formula)
            {
                continue;
            }
            let mut s_elements = 0.0;
            let mut all_present = true;
            for (symbol, count) in &composition {
                let ref_data = reference_state(symbol)
                    .and_then(|(f, atoms)| lookup_thermochem(f).map(|r| (r, atoms)));
                match ref_data {
                    Some((r, atoms)) => s_elements += count / atoms * r.s_standard_j,
                    None => all_present = false,
                }
            }
            if !all_present {
                continue;
            }
            checked += 1;
            let dg_calc = d.delta_hf_kj - 298.15 * (d.s_standard_j - s_elements) / 1000.0;
            let diff = dg_calc - d.delta_gf_kj;
            if let Some(&(_, pinned)) = KNOWN_INCONSISTENT.iter().find(|(f, _)| *f == d.formula) {
                assert!(
                    (diff - pinned).abs() < 0.05,
                    "{}: deviation changed to {diff:.2} — re-check data",
                    d.formula
                );
                continue;
            }
            assert!(
                diff.abs() < 2.0,
                "{}: ΔG_calc = {dg_calc:.2}, table {:.2} (diff {diff:.2} kJ/mol)",
                d.formula,
                d.delta_gf_kj
            );
        }
        // Every non-element row has all its elements in the table.
        assert_eq!(checked, THERMOCHEM_DATA.len() - 22, "rows checked");
    }

    #[test]
    fn high_temperature_segments_match_janaf() {
        // JANAF Cp° at 2000 K (J/(mol·K)): H2 34.280, CO 36.250, CH4 94.399; H2 at 3000 K 37.087.
        let cases = [
            ("H2(g)", 2000.0, 34.280),
            ("H2(g)", 3000.0, 37.087),
            ("CO(g)", 2000.0, 36.250),
            ("CH4(g)", 2000.0, 94.399),
        ];
        for (formula, t, expected) in cases {
            let cp = lookup_shomate_at(formula, t).unwrap().cp(t);
            assert!(
                (cp - expected).abs() < 0.1,
                "{formula} Cp({t}) = {cp}, JANAF {expected}"
            );
        }
    }

    #[test]
    fn adiabatic_flame_rich_methane_with_co_product() {
        // CH4 + 1.5 O2 (+ 5.64 N2) -> CO + 2 H2O (+ 5.64 N2):
        // ΔH = -110.53 - 2(241.83) + 74.87 = -519.32 kJ/mol.
        // CO must be usable above 1300 K for this to solve.
        let products = [("CO(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 5.64)];
        let t = adiabatic_flame_temperature(-519_320.0, &products, 298.15).unwrap();
        assert!(t > 1300.0 && t < 3000.0, "rich CH4 flame T = {t}");
        let absorbed: f64 = products
            .iter()
            .map(|&(f, n)| n * enthalpy_change_cp(f, 298.15, t).unwrap())
            .sum();
        assert!(
            (absorbed - 519_320.0).abs() < 1.0,
            "energy balance off: {absorbed}"
        );
    }
}
