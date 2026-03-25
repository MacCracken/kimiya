//! Electrochemistry — electrode potentials, Nernst equation, Faraday's laws.
//!
//! Standard electrode potentials are given for reduction half-reactions at 25°C,
//! 1 atm, 1 M concentrations, referenced to the Standard Hydrogen Electrode (SHE).

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};
use serde::Serialize;

// ── Constants ────────────────────────────────────────────────────────

/// Faraday constant (C/mol) — charge of one mole of electrons.
pub const FARADAY: f64 = 96_485.332_12;

/// Elementary charge (C).
pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

// ── Standard electrode potentials ────────────────────────────────────

/// A standard reduction half-reaction with its electrode potential.
#[derive(Debug, Clone, Serialize)]
pub struct HalfReaction {
    /// The reduction half-reaction equation (e.g., "Cu²⁺ + 2e⁻ → Cu").
    pub equation: &'static str,
    /// Species being reduced (e.g., "Cu2+/Cu").
    pub couple: &'static str,
    /// Number of electrons transferred.
    pub electrons: u8,
    /// Standard reduction potential E° (V vs SHE).
    pub standard_potential: f64,
}

/// Built-in standard reduction potentials (V vs SHE at 25°C).
///
/// Ordered from most negative (strongest reductant) to most positive (strongest oxidant).
pub static STANDARD_POTENTIALS: &[HalfReaction] = &[
    HalfReaction {
        equation: "Li⁺ + e⁻ → Li",
        couple: "Li+/Li",
        electrons: 1,
        standard_potential: -3.04,
    },
    HalfReaction {
        equation: "K⁺ + e⁻ → K",
        couple: "K+/K",
        electrons: 1,
        standard_potential: -2.924,
    },
    HalfReaction {
        equation: "Ca²⁺ + 2e⁻ → Ca",
        couple: "Ca2+/Ca",
        electrons: 2,
        standard_potential: -2.868,
    },
    HalfReaction {
        equation: "Na⁺ + e⁻ → Na",
        couple: "Na+/Na",
        electrons: 1,
        standard_potential: -2.71,
    },
    HalfReaction {
        equation: "Mg²⁺ + 2e⁻ → Mg",
        couple: "Mg2+/Mg",
        electrons: 2,
        standard_potential: -2.372,
    },
    HalfReaction {
        equation: "Al³⁺ + 3e⁻ → Al",
        couple: "Al3+/Al",
        electrons: 3,
        standard_potential: -1.662,
    },
    HalfReaction {
        equation: "Zn²⁺ + 2e⁻ → Zn",
        couple: "Zn2+/Zn",
        electrons: 2,
        standard_potential: -0.762,
    },
    HalfReaction {
        equation: "Fe²⁺ + 2e⁻ → Fe",
        couple: "Fe2+/Fe",
        electrons: 2,
        standard_potential: -0.447,
    },
    HalfReaction {
        equation: "Ni²⁺ + 2e⁻ → Ni",
        couple: "Ni2+/Ni",
        electrons: 2,
        standard_potential: -0.257,
    },
    HalfReaction {
        equation: "Sn²⁺ + 2e⁻ → Sn",
        couple: "Sn2+/Sn",
        electrons: 2,
        standard_potential: -0.136,
    },
    HalfReaction {
        equation: "Pb²⁺ + 2e⁻ → Pb",
        couple: "Pb2+/Pb",
        electrons: 2,
        standard_potential: -0.126,
    },
    HalfReaction {
        equation: "2H⁺ + 2e⁻ → H₂",
        couple: "H+/H2",
        electrons: 2,
        standard_potential: 0.000,
    },
    HalfReaction {
        equation: "Cu²⁺ + 2e⁻ → Cu",
        couple: "Cu2+/Cu",
        electrons: 2,
        standard_potential: 0.342,
    },
    HalfReaction {
        equation: "I₂ + 2e⁻ → 2I⁻",
        couple: "I2/I-",
        electrons: 2,
        standard_potential: 0.536,
    },
    HalfReaction {
        equation: "Ag⁺ + e⁻ → Ag",
        couple: "Ag+/Ag",
        electrons: 1,
        standard_potential: 0.7996,
    },
    HalfReaction {
        equation: "Br₂ + 2e⁻ → 2Br⁻",
        couple: "Br2/Br-",
        electrons: 2,
        standard_potential: 1.066,
    },
    HalfReaction {
        equation: "Cl₂ + 2e⁻ → 2Cl⁻",
        couple: "Cl2/Cl-",
        electrons: 2,
        standard_potential: 1.358,
    },
    HalfReaction {
        equation: "Au³⁺ + 3e⁻ → Au",
        couple: "Au3+/Au",
        electrons: 3,
        standard_potential: 1.498,
    },
    HalfReaction {
        equation: "F₂ + 2e⁻ → 2F⁻",
        couple: "F2/F-",
        electrons: 2,
        standard_potential: 2.866,
    },
];

/// Look up a half-reaction by its couple string (e.g., "Cu2+/Cu").
#[must_use]
pub fn lookup_half_reaction(couple: &str) -> Option<&'static HalfReaction> {
    STANDARD_POTENTIALS.iter().find(|hr| hr.couple == couple)
}

// ── Nernst equation ──────────────────────────────────────────────────

/// Nernst equation: E = E° − (RT / nF) ln(Q)
///
/// Calculates the electrode potential under non-standard conditions.
///
/// - `standard_potential`: E° in volts
/// - `n_electrons`: number of electrons transferred
/// - `temperature_k`: temperature in Kelvin
/// - `reaction_quotient`: Q (ratio of product activities to reactant activities)
///
/// # Errors
///
/// Returns error if temperature, electron count, or reaction quotient is invalid.
#[inline]
pub fn nernst_potential(
    standard_potential: f64,
    n_electrons: u8,
    temperature_k: f64,
    reaction_quotient: f64,
) -> Result<f64> {
    if temperature_k <= 0.0 {
        return Err(KimiyaError::InvalidTemperature(
            "temperature must be positive".into(),
        ));
    }
    if n_electrons == 0 {
        return Err(KimiyaError::InvalidInput(
            "electron count must be at least 1".into(),
        ));
    }
    if reaction_quotient <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "reaction quotient must be positive".into(),
        ));
    }
    Ok(standard_potential
        - (GAS_CONSTANT * temperature_k / (n_electrons as f64 * FARADAY)) * reaction_quotient.ln())
}

/// Simplified Nernst equation at 25°C (298.15 K) using log₁₀:
/// E = E° − (0.05916 / n) log₁₀(Q)
///
/// # Errors
///
/// Returns error if electron count or reaction quotient is invalid.
#[inline]
pub fn nernst_potential_25c(
    standard_potential: f64,
    n_electrons: u8,
    reaction_quotient: f64,
) -> Result<f64> {
    nernst_potential(standard_potential, n_electrons, 298.15, reaction_quotient)
}

// ── Cell potential ───────────────────────────────────────────────────

/// Standard cell potential: E°_cell = E°_cathode − E°_anode.
///
/// Both potentials are standard reduction potentials.
#[must_use]
#[inline]
pub fn standard_cell_potential(e_cathode: f64, e_anode: f64) -> f64 {
    e_cathode - e_anode
}

/// Standard cell potential from couple strings.
///
/// `cathode_couple` is reduced (gains electrons), `anode_couple` is oxidized.
///
/// # Errors
///
/// Returns error if either couple is not found in the standard potentials table.
pub fn cell_potential_from_couples(cathode_couple: &str, anode_couple: &str) -> Result<f64> {
    let cathode = lookup_half_reaction(cathode_couple).ok_or_else(|| {
        KimiyaError::InvalidReaction(format!("unknown cathode couple: {cathode_couple}"))
    })?;
    let anode = lookup_half_reaction(anode_couple).ok_or_else(|| {
        KimiyaError::InvalidReaction(format!("unknown anode couple: {anode_couple}"))
    })?;
    Ok(standard_cell_potential(
        cathode.standard_potential,
        anode.standard_potential,
    ))
}

/// Check if a galvanic cell is spontaneous (E°_cell > 0).
#[must_use]
#[inline]
pub fn is_spontaneous_cell(e_cell: f64) -> bool {
    e_cell > 0.0
}

/// Gibbs free energy of an electrochemical cell: ΔG = −nFE
///
/// Returns energy in Joules per mole of reaction.
///
/// # Errors
///
/// Returns error if electron count is zero.
#[inline]
pub fn cell_gibbs_energy(n_electrons: u8, cell_potential: f64) -> Result<f64> {
    if n_electrons == 0 {
        return Err(KimiyaError::InvalidInput(
            "electron count must be at least 1".into(),
        ));
    }
    Ok(-(n_electrons as f64) * FARADAY * cell_potential)
}

// ── Faraday's laws ───────────────────────────────────────────────────

/// Faraday's first law: mass deposited m = (Q × M) / (n × F)
///
/// - `charge_c`: total charge passed (Coulombs)
/// - `molar_mass_g`: molar mass of substance (g/mol)
/// - `n_electrons`: electrons per ion (valence)
///
/// Returns mass in grams.
///
/// # Errors
///
/// Returns error if electron count is zero.
#[inline]
pub fn faraday_mass_deposited(charge_c: f64, molar_mass_g: f64, n_electrons: u8) -> Result<f64> {
    if n_electrons == 0 {
        return Err(KimiyaError::InvalidInput(
            "electron count must be at least 1".into(),
        ));
    }
    Ok(charge_c * molar_mass_g / (n_electrons as f64 * FARADAY))
}

/// Charge required to deposit a given mass: Q = (m × n × F) / M
///
/// Returns charge in Coulombs.
///
/// # Errors
///
/// Returns error if molar mass is not positive or electron count is zero.
#[inline]
pub fn faraday_charge_required(mass_g: f64, molar_mass_g: f64, n_electrons: u8) -> Result<f64> {
    if n_electrons == 0 {
        return Err(KimiyaError::InvalidInput(
            "electron count must be at least 1".into(),
        ));
    }
    if molar_mass_g <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "molar mass must be positive".into(),
        ));
    }
    Ok(mass_g * n_electrons as f64 * FARADAY / molar_mass_g)
}

/// Faraday's second law: mass ratio of substances deposited by the same charge
/// equals the ratio of their equivalent weights (M/n).
///
/// Returns the mass of substance B deposited when `mass_a_g` of substance A
/// is deposited by the same charge.
///
/// # Errors
///
/// Returns error if electron counts are zero or molar mass A is not positive.
#[inline]
pub fn faraday_mass_ratio(
    mass_a_g: f64,
    molar_mass_a: f64,
    n_electrons_a: u8,
    molar_mass_b: f64,
    n_electrons_b: u8,
) -> Result<f64> {
    if n_electrons_a == 0 || n_electrons_b == 0 {
        return Err(KimiyaError::InvalidInput(
            "electron counts must be at least 1".into(),
        ));
    }
    if molar_mass_a <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "molar mass A must be positive".into(),
        ));
    }
    let equiv_a = molar_mass_a / n_electrons_a as f64;
    let equiv_b = molar_mass_b / n_electrons_b as f64;
    Ok(mass_a_g * equiv_b / equiv_a)
}

// ── Utility ──────────────────────────────────────────────────────────

/// Convert charge to moles of electrons: n_e = Q / F
#[must_use]
#[inline]
pub fn charge_to_moles_electrons(charge_c: f64) -> f64 {
    charge_c / FARADAY
}

/// Convert moles of electrons to charge: Q = n_e × F
#[must_use]
#[inline]
pub fn moles_electrons_to_charge(moles: f64) -> f64 {
    moles * FARADAY
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Constants ────────────────────────────────────────────────────

    #[test]
    fn faraday_constant_value() {
        // F = N_A × e
        let computed = crate::element::AVOGADRO * ELEMENTARY_CHARGE;
        assert!(
            (computed - FARADAY).abs() / FARADAY < 1e-6,
            "F should equal N_A × e, got {computed}"
        );
    }

    // ── Half-reaction lookup ─────────────────────────────────────────

    #[test]
    fn lookup_copper() {
        let hr = lookup_half_reaction("Cu2+/Cu").unwrap();
        assert_eq!(hr.electrons, 2);
        assert!((hr.standard_potential - 0.342).abs() < 0.001);
    }

    #[test]
    fn lookup_hydrogen() {
        let hr = lookup_half_reaction("H+/H2").unwrap();
        assert!((hr.standard_potential).abs() < f64::EPSILON);
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_half_reaction("Xx/Xx").is_none());
    }

    #[test]
    fn potentials_ordered() {
        for window in STANDARD_POTENTIALS.windows(2) {
            assert!(
                window[0].standard_potential <= window[1].standard_potential,
                "{} ({}) should be ≤ {} ({})",
                window[0].couple,
                window[0].standard_potential,
                window[1].couple,
                window[1].standard_potential
            );
        }
    }

    // ── Nernst equation ──────────────────────────────────────────────

    #[test]
    fn nernst_at_standard_conditions() {
        // Q=1 → E = E°
        let e = nernst_potential(0.342, 2, 298.15, 1.0).unwrap();
        assert!((e - 0.342).abs() < 1e-10, "Q=1 should give E=E°, got {e}");
    }

    #[test]
    fn nernst_q_greater_than_1() {
        // Q > 1 → E < E° for reduction
        let e = nernst_potential(0.342, 2, 298.15, 100.0).unwrap();
        assert!(e < 0.342, "Q>1 should reduce E below E°, got {e}");
    }

    #[test]
    fn nernst_q_less_than_1() {
        // Q < 1 → E > E°
        let e = nernst_potential(0.342, 2, 298.15, 0.01).unwrap();
        assert!(e > 0.342, "Q<1 should increase E above E°, got {e}");
    }

    #[test]
    fn nernst_25c_shorthand() {
        let e_full = nernst_potential(0.342, 2, 298.15, 0.1).unwrap();
        let e_short = nernst_potential_25c(0.342, 2, 0.1).unwrap();
        assert!((e_full - e_short).abs() < 1e-10);
    }

    #[test]
    fn nernst_zero_electrons_is_error() {
        assert!(nernst_potential(0.342, 0, 298.15, 1.0).is_err());
    }

    #[test]
    fn nernst_zero_q_is_error() {
        assert!(nernst_potential(0.342, 2, 298.15, 0.0).is_err());
    }

    #[test]
    fn nernst_zero_temp_is_error() {
        assert!(nernst_potential(0.342, 2, 0.0, 1.0).is_err());
    }

    // ── Cell potential ───────────────────────────────────────────────

    #[test]
    fn daniell_cell() {
        // Cu²⁺/Cu (cathode, +0.342) vs Zn²⁺/Zn (anode, -0.762)
        let e = standard_cell_potential(0.342, -0.762);
        assert!(
            (e - 1.104).abs() < 0.001,
            "Daniell cell should be ~1.10V, got {e}"
        );
    }

    #[test]
    fn daniell_cell_from_couples() {
        let e = cell_potential_from_couples("Cu2+/Cu", "Zn2+/Zn").unwrap();
        assert!((e - 1.104).abs() < 0.001);
    }

    #[test]
    fn unknown_couple_is_error() {
        assert!(cell_potential_from_couples("Xx/Xx", "Zn2+/Zn").is_err());
    }

    #[test]
    fn spontaneous_cell_positive() {
        assert!(is_spontaneous_cell(1.1));
        assert!(!is_spontaneous_cell(-0.5));
        assert!(!is_spontaneous_cell(0.0));
    }

    // ── Gibbs energy ─────────────────────────────────────────────────

    #[test]
    fn cell_gibbs_daniell() {
        // ΔG = -nFE = -2 × 96485 × 1.104 ≈ -213,038 J/mol
        let dg = cell_gibbs_energy(2, 1.104).unwrap();
        assert!(dg < 0.0, "spontaneous cell should have negative ΔG");
        assert!(
            (dg - (-213_038.0)).abs() < 100.0,
            "ΔG should be ~-213 kJ/mol, got {dg}"
        );
    }

    #[test]
    fn cell_gibbs_zero_electrons_is_error() {
        assert!(cell_gibbs_energy(0, 1.0).is_err());
    }

    // ── Faraday's laws ───────────────────────────────────────────────

    #[test]
    fn faraday_copper_deposition() {
        // 1 Faraday of charge deposits M/(n) = 63.546/2 = 31.773 g of Cu
        let mass = faraday_mass_deposited(FARADAY, 63.546, 2).unwrap();
        assert!(
            (mass - 31.773).abs() < 0.01,
            "1F should deposit ~31.77g Cu, got {mass}"
        );
    }

    #[test]
    fn faraday_charge_roundtrip() {
        let mass = 10.0; // 10g Cu
        let q = faraday_charge_required(mass, 63.546, 2).unwrap();
        let back = faraday_mass_deposited(q, 63.546, 2).unwrap();
        assert!(
            (back - mass).abs() < 0.001,
            "roundtrip failed: {mass} → {q}C → {back}g"
        );
    }

    #[test]
    fn faraday_zero_electrons_is_error() {
        assert!(faraday_mass_deposited(1000.0, 63.546, 0).is_err());
        assert!(faraday_charge_required(10.0, 63.546, 0).is_err());
    }

    #[test]
    fn faraday_second_law() {
        // Same charge deposits Cu and Ag in ratio of equivalent weights
        // Cu: M=63.546, n=2, equiv=31.773
        // Ag: M=107.868, n=1, equiv=107.868
        // If 31.773g Cu deposited, should get 107.868g Ag
        let mass_ag = faraday_mass_ratio(31.773, 63.546, 2, 107.868, 1).unwrap();
        assert!(
            (mass_ag - 107.868).abs() < 0.1,
            "should deposit ~107.87g Ag, got {mass_ag}"
        );
    }

    // ── Utility ──────────────────────────────────────────────────────

    #[test]
    fn charge_moles_roundtrip() {
        let moles = 2.5;
        let q = moles_electrons_to_charge(moles);
        let back = charge_to_moles_electrons(q);
        assert!((back - moles).abs() < 1e-10);
    }

    #[test]
    fn one_faraday_is_one_mole() {
        let moles = charge_to_moles_electrons(FARADAY);
        assert!((moles - 1.0).abs() < 1e-10);
    }
}
