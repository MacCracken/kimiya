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
/// k/A = exp(-Ea / (R × T))
#[must_use]
#[inline]
pub fn activation_energy_to_rate_factor(
    activation_energy_j_per_mol: f64,
    temperature_k: f64,
) -> f64 {
    const R: f64 = 8.314; // J/(mol·K)
    if temperature_k <= 0.0 {
        return 0.0;
    }
    (-activation_energy_j_per_mol / (R * temperature_k)).exp()
}

/// Convert equilibrium constant K and temperature (K) to standard
/// Gibbs free energy of reaction (J/mol).
///
/// ΔG° = -R × T × ln(K)
#[must_use]
#[inline]
pub fn equilibrium_to_gibbs(equilibrium_constant: f64, temperature_k: f64) -> f64 {
    const R: f64 = 8.314;
    if equilibrium_constant <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    -R * temperature_k * equilibrium_constant.ln()
}

// ── Tanmatra bridges (atomic/nuclear physics) ──────────────────────────────

/// Convert atomic number to number of valence electrons (simplified).
///
/// Uses the main-group rule: valence = group number (s and p block).
/// For transition metals, returns a conventional value.
#[must_use]
pub fn atomic_number_to_valence(atomic_number: u8) -> u8 {
    // Period and group from atomic number (simplified main-group only)
    match atomic_number {
        1 => 1,
        2 => 2,                        // He (noble gas, but 2 electrons)
        3..=4 => atomic_number - 2,    // Li=1, Be=2
        5..=10 => atomic_number - 2,   // B=3, C=4, N=5, O=6, F=7, Ne=8
        11..=12 => atomic_number - 10, // Na=1, Mg=2
        13..=18 => atomic_number - 10, // Al=3, Si=4, P=5, S=6, Cl=7, Ar=8
        _ => {
            // Transition metals and beyond — conventional 2
            2
        }
    }
}

/// Convert ionization energy (eV) to approximate bond dissociation energy
/// scaling for homonuclear diatomics (empirical correlation).
///
/// D_e ≈ 0.5 × IE for many light diatomics. Returns eV.
#[must_use]
#[inline]
pub fn ionization_to_bond_energy_estimate(ionization_energy_ev: f64) -> f64 {
    // Very rough: bond energy scales with IE for light elements
    (0.5 * ionization_energy_ev).max(0.0)
}

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert electrochemical cell potential (V) and electron count to
/// current density estimate (A/m²).
///
/// Uses Butler-Volmer at low overpotential: j ≈ j₀ × n × F × η / (R × T)
/// `cell_potential_v`: electrode potential.
/// `equilibrium_potential_v`: equilibrium (Nernst) potential.
/// `exchange_current_density_a_m2`: j₀ (typically 1e-3 to 1e3 for metals).
/// `temperature_k`: electrolyte temperature.
/// `n_electrons`: number of electrons transferred.
#[must_use]
pub fn electrochemical_current_density(
    cell_potential_v: f64,
    equilibrium_potential_v: f64,
    exchange_current_density_a_m2: f64,
    temperature_k: f64,
    n_electrons: u32,
) -> f64 {
    const F: f64 = 96_485.332; // C/mol
    const R: f64 = 8.314; // J/(mol·K)
    if temperature_k <= 0.0 {
        return 0.0;
    }
    let eta = cell_potential_v - equilibrium_potential_v; // overpotential
    let n = n_electrons as f64;
    // Low-overpotential Butler-Volmer
    exchange_current_density_a_m2 * n * F * eta / (R * temperature_k)
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
        // Ea = 50 kJ/mol at 298K → exp(-50000/(8.314×298)) ≈ 1.7e-9
        let f = activation_energy_to_rate_factor(50_000.0, 298.0);
        assert!(f > 1e-10 && f < 1e-8);
    }

    #[test]
    fn arrhenius_factor_zero_temp() {
        assert_eq!(activation_energy_to_rate_factor(50_000.0, 0.0), 0.0);
    }

    #[test]
    fn gibbs_from_equilibrium() {
        // K = 1.0 → ΔG° = 0
        let g = equilibrium_to_gibbs(1.0, 298.0);
        assert!(g.abs() < 0.01);
    }

    #[test]
    fn gibbs_favorable() {
        // K > 1 → ΔG° < 0
        let g = equilibrium_to_gibbs(100.0, 298.0);
        assert!(g < 0.0);
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
    fn bond_energy_estimate() {
        // H ionization ≈ 13.6 eV → bond ≈ 6.8 eV (actual H₂ ≈ 4.5 eV, rough)
        let d = ionization_to_bond_energy_estimate(13.6);
        assert!((d - 6.8).abs() < 0.01);
    }

    // ── Bijli ──────────────────────────────────────────────────────────

    #[test]
    fn electrochemical_current_basic() {
        // 0.1V overpotential, j₀ = 1.0 A/m², n=2, T=298K
        let j = electrochemical_current_density(0.1, 0.0, 1.0, 298.0, 2);
        assert!(j > 0.0);
    }

    #[test]
    fn electrochemical_current_zero_overpotential() {
        let j = electrochemical_current_density(0.5, 0.5, 1.0, 298.0, 1);
        assert!(j.abs() < 0.01);
    }

    #[test]
    fn electrochemical_current_zero_temp() {
        assert_eq!(electrochemical_current_density(0.1, 0.0, 1.0, 0.0, 1), 0.0);
    }

    #[test]
    fn conductivity_basic() {
        // 1 mol/L → κ = 1000 × 0.01 = 10 S/m
        let k = ionic_concentration_to_conductivity(1.0);
        assert!((k - 10.0).abs() < 0.01);
    }
}
