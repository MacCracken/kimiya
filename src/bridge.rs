//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into kimiya chemistry parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! ushma    (thermodynamics) ──┐
//! tanmatra (atomic physics)   ┼──> bridge ──> kimiya chemistry parameters
//! bijli    (electromagnetism)┘
//! ```
//!
//! Functions here return a bare `f64`. Where an input makes the formula
//! undefined (e.g. a non-positive temperature) they return `f64::NAN`, never
//! a plausible-looking number such as 0.0, so invalid input stays detectable
//! with `is_nan()`.

use crate::electrochemistry::FARADAY;
use crate::element::GAS_CONSTANT;
use crate::validate::is_positive;

// ── Ushma bridges (thermodynamics) ─────────────────────────────────────────

/// Convert reaction enthalpy (J/mol) and molar flow (mol/s) to heat
/// release rate (W).
///
/// Q̇ = -ΔH × ṅ (negative for exothermic releases positive heat).
#[must_use]
#[inline]
pub fn reaction_heat_release(enthalpy_j_per_mol: f64, molar_rate_mol_per_s: f64) -> f64 {
    -enthalpy_j_per_mol * molar_rate_mol_per_s
}

/// Convert activation energy (J/mol) and temperature (K) to Arrhenius
/// rate scaling factor.
///
/// k/A = exp(-Ea / (R × T)) (dimensionless)
///
/// Returns `f64::NAN` if the temperature is not positive and finite.
#[must_use]
#[inline]
pub fn activation_energy_to_rate_factor(
    activation_energy_j_per_mol: f64,
    temperature_k: f64,
) -> f64 {
    if !is_positive(temperature_k) {
        return f64::NAN;
    }
    (-activation_energy_j_per_mol / (GAS_CONSTANT * temperature_k)).exp()
}

/// Convert equilibrium constant K and temperature (K) to standard
/// Gibbs free energy of reaction (J/mol).
///
/// ΔG° = -R × T × ln(K)
///
/// Returns `f64::NAN` if K or the temperature is not positive and finite
/// (previously 0.0, which is indistinguishable from the valid result for
/// K = 1).
#[must_use]
#[inline]
pub fn equilibrium_to_gibbs(equilibrium_constant: f64, temperature_k: f64) -> f64 {
    if !is_positive(equilibrium_constant) || !is_positive(temperature_k) {
        return f64::NAN;
    }
    -GAS_CONSTANT * temperature_k * equilibrium_constant.ln()
}

// ── Tanmatra bridges (atomic/nuclear physics) ──────────────────────────────

/// Convert atomic number to number of valence electrons (simplified).
///
/// Convention, for Z = 1–118:
///
/// - **s and p block (main group):** the outer-shell s + p electron count,
///   i.e. the group number for groups 1–2 and group − 10 for groups 13–18
///   (H 1, He 2, C 4, Cl 7, Ar 8, K 1, Ga 3, Br 7, I 7, Cs 1, Pb 4, Rn 8,
///   Og 8).
/// - **d and f block** (transition metals, lanthanides, actinides:
///   Z = 21–30, 39–48, 57–80, 89–112): the conventional value 2, the
///   nominal outer ns² pair. (n−1)d / (n−2)f electrons and ground-state
///   anomalies such as Cr/Cu (4s¹) or Pd (5s⁰) are ignored.
/// - **Invalid Z** (0 or > 118): returns 0.
///
/// This counts valence *electrons*, not bonding valence (combining
/// capacity): oxygen gives 6, not 2.
#[must_use]
#[inline]
pub fn atomic_number_to_valence(atomic_number: u8) -> u8 {
    match atomic_number {
        1 => 1,                           // H
        2 => 2,                           // He (1s², noble gas)
        3..=10 => atomic_number - 2,      // Li 1, Be 2, B 3 … Ne 8
        11..=18 => atomic_number - 10,    // Na 1, Mg 2, Al 3 … Ar 8
        19..=20 => atomic_number - 18,    // K 1, Ca 2
        31..=36 => atomic_number - 28,    // Ga 3 … Kr 8
        37..=38 => atomic_number - 36,    // Rb 1, Sr 2
        49..=54 => atomic_number - 46,    // In 3 … Xe 8
        55..=56 => atomic_number - 54,    // Cs 1, Ba 2
        81..=86 => atomic_number - 78,    // Tl 3 … Rn 8
        87..=88 => atomic_number - 86,    // Fr 1, Ra 2
        113..=118 => atomic_number - 110, // Nh 3 … Og 8
        // d and f block: conventional ns² = 2
        21..=30 | 39..=48 | 57..=80 | 89..=112 => 2,
        // Z = 0 or Z > 118: not an element
        _ => 0,
    }
}

/// Crude heuristic for the bond dissociation energy of a homonuclear
/// diatomic from the atomic ionization energy: D ≈ 0.5 × IE (both in eV).
///
/// **This is not a reliable correlation.** It is off by a factor of up to
/// ~5 in either direction (estimate vs. experimental bond energy, eV):
///
/// | Molecule | 0.5 × IE | Experimental |
/// |----------|----------|--------------|
/// | H₂       | 6.80     | 4.48         |
/// | N₂       | 7.27     | 9.79         |
/// | F₂       | 8.71     | 1.60         |
/// | Cl₂      | 6.48     | 2.48         |
///
/// Use only as an order-of-magnitude placeholder when no tabulated bond
/// energy is available. Negative input is clamped to 0.0; non-finite input
/// (NaN, ±∞) returns `f64::NAN`.
#[must_use]
#[inline]
pub fn ionization_to_bond_energy_estimate(ionization_energy_ev: f64) -> f64 {
    if !ionization_energy_ev.is_finite() {
        return f64::NAN;
    }
    (0.5 * ionization_energy_ev).max(0.0)
}

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert electrochemical cell potential (V) and electron count to
/// current density estimate (A/m²).
///
/// Uses the linearised (low-overpotential) Butler-Volmer form:
/// j ≈ j₀ × n × F × η / (R × T), with η = E − E_eq.
///
/// - `cell_potential_v`: electrode potential E (V).
/// - `equilibrium_potential_v`: equilibrium (Nernst) potential E_eq (V).
/// - `exchange_current_density_a_m2`: j₀ (A/m²; typically 1e-3 to 1e3 for
///   metals).
/// - `temperature_k`: electrolyte temperature (K).
/// - `n_electrons`: number of electrons transferred.
///
/// Returns current density in A/m² (positive = anodic).
///
/// The linearisation is only valid for |η| ≪ RT/F ≈ 25.7 mV at 298.15 K
/// (more strictly |η| ≪ RT/(α·n·F)). At larger overpotentials the full
/// Butler-Volmer equation grows exponentially and this underestimates |j|
/// (Tafel regime). Returns `f64::NAN` if the temperature is not positive and
/// finite.
#[must_use]
#[inline]
pub fn electrochemical_current_density(
    cell_potential_v: f64,
    equilibrium_potential_v: f64,
    exchange_current_density_a_m2: f64,
    temperature_k: f64,
    n_electrons: u32,
) -> f64 {
    if !is_positive(temperature_k) {
        return f64::NAN;
    }
    let eta = cell_potential_v - equilibrium_potential_v; // overpotential
    let n = f64::from(n_electrons);
    // Low-overpotential Butler-Volmer
    exchange_current_density_a_m2 * n * FARADAY * eta / (GAS_CONSTANT * temperature_k)
}

/// Convert ionic concentration (mol/L) to electrical conductivity (S/m).
///
/// Simplified Kohlrausch: κ ≈ c × Λ_m where Λ_m ≈ 0.01 S·m²/mol for strong electrolytes.
#[must_use]
#[inline]
pub fn ionic_concentration_to_conductivity(concentration_mol_per_l: f64) -> f64 {
    // Convert mol/L to mol/m³ (* 1000), then κ = c × Λ_m
    let c_mol_m3 = concentration_mol_per_l * 1000.0;
    let lambda_m = 0.01; // S·m²/mol (typical strong electrolyte)
    c_mol_m3 * lambda_m
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ELEMENTS, ElementCategory};

    // ── Ushma ──────────────────────────────────────────────────────────

    #[test]
    fn reaction_heat_exothermic() {
        // ΔH = -100 kJ/mol (exo), 0.1 mol/s → 10 kW
        let q = reaction_heat_release(-100_000.0, 0.1);
        assert!((q - 10_000.0).abs() < 0.1);
    }

    #[test]
    fn reaction_heat_endothermic() {
        // ΔH = +50 kJ/mol (endo), 0.1 mol/s → -5 kW (absorbs heat)
        let q = reaction_heat_release(50_000.0, 0.1);
        assert!((q - (-5_000.0)).abs() < 0.1);
    }

    #[test]
    fn arrhenius_factor_room_temp() {
        // Ea = 50 kJ/mol at 298K → exp(-50000/(8.314462618×298)) = 1.721 75e-9
        // (R = 8.314 gave 1.7236e-9).
        let f = activation_energy_to_rate_factor(50_000.0, 298.0);
        assert!(((f - 1.721_748_757_6e-9) / f).abs() < 1e-9, "got {f}");
    }

    #[test]
    fn arrhenius_factor_invalid_temp_is_nan() {
        assert!(activation_energy_to_rate_factor(50_000.0, 0.0).is_nan());
        assert!(activation_energy_to_rate_factor(50_000.0, -10.0).is_nan());
        assert!(activation_energy_to_rate_factor(50_000.0, f64::NAN).is_nan());
    }

    #[test]
    fn gibbs_from_equilibrium() {
        // K = 1.0 → ΔG° = 0
        let g = equilibrium_to_gibbs(1.0, 298.0);
        assert!(g.abs() < f64::EPSILON);
    }

    #[test]
    fn gibbs_favorable() {
        // K = 100 at 298 K → ΔG° = −8.314462618 × 298 × ln 100 = −11 410.2756 J/mol
        // (R = 8.314 gave −11 409.64).
        let g = equilibrium_to_gibbs(100.0, 298.0);
        assert!((g - (-11_410.275_577_6)).abs() < 1e-6, "got {g}");
    }

    #[test]
    fn gibbs_from_invalid_equilibrium_is_nan() {
        // Used to return 0.0 — the valid answer for K = 1.
        assert!(equilibrium_to_gibbs(0.0, 298.0).is_nan());
        assert!(equilibrium_to_gibbs(-1.0, 298.0).is_nan());
        assert!(equilibrium_to_gibbs(f64::NAN, 298.0).is_nan());
        assert!(equilibrium_to_gibbs(f64::INFINITY, 298.0).is_nan());
        assert!(equilibrium_to_gibbs(100.0, 0.0).is_nan());
        assert!(equilibrium_to_gibbs(100.0, f64::NAN).is_nan());
    }

    // ── Tanmatra ───────────────────────────────────────────────────────

    #[test]
    fn valence_carbon() {
        assert_eq!(atomic_number_to_valence(6), 4);
    }

    #[test]
    fn valence_sodium() {
        assert_eq!(atomic_number_to_valence(11), 1);
    }

    #[test]
    fn valence_chlorine() {
        assert_eq!(atomic_number_to_valence(17), 7);
    }

    #[test]
    fn valence_transition_metal() {
        // Iron (26) → conventional 2
        assert_eq!(atomic_number_to_valence(26), 2);
    }

    #[test]
    fn valence_heavier_main_group() {
        // All of these used to return 2.
        assert_eq!(atomic_number_to_valence(19), 1); // K
        assert_eq!(atomic_number_to_valence(31), 3); // Ga
        assert_eq!(atomic_number_to_valence(35), 7); // Br
        assert_eq!(atomic_number_to_valence(37), 1); // Rb
        assert_eq!(atomic_number_to_valence(53), 7); // I
        assert_eq!(atomic_number_to_valence(55), 1); // Cs
        assert_eq!(atomic_number_to_valence(82), 4); // Pb
        assert_eq!(atomic_number_to_valence(86), 8); // Rn
        assert_eq!(atomic_number_to_valence(118), 8); // Og
    }

    #[test]
    fn valence_d_and_f_block_conventional() {
        for z in [21, 30, 39, 48, 57, 64, 71, 72, 80, 89, 92, 103, 104, 112] {
            assert_eq!(atomic_number_to_valence(z), 2, "Z = {z}");
        }
    }

    #[test]
    fn valence_invalid_atomic_number_is_zero() {
        assert_eq!(atomic_number_to_valence(0), 0);
        assert_eq!(atomic_number_to_valence(119), 0);
        assert_eq!(atomic_number_to_valence(u8::MAX), 0);
    }

    #[test]
    fn valence_matches_element_categories() {
        for el in ELEMENTS {
            let v = atomic_number_to_valence(el.atomic_number);
            assert!((1..=8).contains(&v), "{}", el.symbol);
            let expected = match el.category {
                ElementCategory::AlkaliMetal => Some(1),
                ElementCategory::AlkalineEarth => Some(2),
                ElementCategory::Halogen => Some(7),
                ElementCategory::NobleGas if el.atomic_number != 2 => Some(8),
                ElementCategory::TransitionMetal
                | ElementCategory::Lanthanide
                | ElementCategory::Actinide => Some(2),
                _ => None,
            };
            if let Some(e) = expected {
                assert_eq!(v, e, "{} ({:?})", el.symbol, el.category);
            }
        }
    }

    #[test]
    fn bond_energy_estimate() {
        // H ionization ≈ 13.6 eV → bond ≈ 6.8 eV (actual H₂ ≈ 4.5 eV, rough)
        let d = ionization_to_bond_energy_estimate(13.6);
        assert!((d - 6.8).abs() < 1e-12);
        assert!(ionization_to_bond_energy_estimate(-1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn bond_energy_estimate_non_finite_is_nan() {
        // `NaN.max(0.0)` used to turn NaN into 0.0.
        assert!(ionization_to_bond_energy_estimate(f64::NAN).is_nan());
        assert!(ionization_to_bond_energy_estimate(f64::INFINITY).is_nan());
        assert!(ionization_to_bond_energy_estimate(f64::NEG_INFINITY).is_nan());
    }

    // ── Bijli ──────────────────────────────────────────────────────────

    #[test]
    fn electrochemical_current_basic() {
        // 0.1 V is ~4 RT/F, well outside the linear regime (|η| ≪ 25.7 mV):
        // this only checks the formula, not physical accuracy.
        // j = 1 × 2 × F × 0.1 / (R × 298) = 7.788 267 A/m²
        // (R = 8.314, F = 96 485.332 gave 7.788 70).
        let j = electrochemical_current_density(0.1, 0.0, 1.0, 298.0, 2);
        assert!(((j - 7.788_267_195_5) / j).abs() < 1e-10, "got {j}");
    }

    #[test]
    fn electrochemical_current_small_overpotential() {
        // η = −5 mV (cathodic), inside the linear regime.
        let j = electrochemical_current_density(-0.005, 0.0, 2.0, 298.15, 1);
        let expected = 2.0 * FARADAY * -0.005 / (GAS_CONSTANT * 298.15);
        assert!((j - expected).abs() < 1e-12 && j < 0.0);
    }

    #[test]
    fn electrochemical_current_zero_overpotential() {
        let j = electrochemical_current_density(0.5, 0.5, 1.0, 298.0, 1);
        assert!(j.abs() < 0.01);
    }

    #[test]
    fn electrochemical_current_invalid_temp_is_nan() {
        assert!(electrochemical_current_density(0.1, 0.0, 1.0, 0.0, 1).is_nan());
        assert!(electrochemical_current_density(0.1, 0.0, 1.0, -5.0, 1).is_nan());
        assert!(electrochemical_current_density(0.1, 0.0, 1.0, f64::NAN, 1).is_nan());
    }

    #[test]
    fn conductivity_basic() {
        // 1 mol/L → κ = 1000 × 0.01 = 10 S/m
        let k = ionic_concentration_to_conductivity(1.0);
        assert!((k - 10.0).abs() < 0.01);
    }
}
