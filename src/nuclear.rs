//! Nuclear chemistry — radioactive decay, binding energy, isotope data.

use crate::error::{KimiyaError, Result};
use serde::Serialize;

// ── Constants ────────────────────────────────────────────────────────

/// Atomic mass unit in MeV/c² (CODATA 2018).
pub const AMU_MEV: f64 = 931.494_102_22;

/// Proton mass in atomic mass units (CODATA 2018).
pub const PROTON_MASS_U: f64 = 1.007_276_466_88;

/// Neutron mass in atomic mass units (CODATA 2018).
pub const NEUTRON_MASS_U: f64 = 1.008_664_915_95;

// ── Radioactive decay ────────────────────────────────────────────────

/// Decay constant from half-life: λ = ln(2) / t½
///
/// # Errors
///
/// Returns error if half-life is not positive.
#[inline]
pub fn decay_constant(half_life: f64) -> Result<f64> {
    if half_life <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "half-life must be positive".into(),
        ));
    }
    Ok(std::f64::consts::LN_2 / half_life)
}

/// Number of atoms remaining: N(t) = N₀ × exp(-λt)
#[must_use]
#[inline]
pub fn atoms_remaining(initial: f64, decay_const: f64, time: f64) -> f64 {
    initial * (-decay_const * time).exp()
}

/// Activity: A = λN (decays per second, Bq)
#[must_use]
#[inline]
pub fn activity(decay_const: f64, num_atoms: f64) -> f64 {
    decay_const * num_atoms
}

/// Number of atoms remaining after n half-lives: N = N₀ / 2ⁿ
#[must_use]
#[inline]
pub fn atoms_after_half_lives(initial: f64, n_half_lives: f64) -> f64 {
    initial * 2.0_f64.powf(-n_half_lives)
}

/// Time for a given fraction to remain: t = -ln(fraction) / λ
///
/// # Errors
///
/// Returns error if fraction is not in (0, 1] or decay constant is not positive.
#[inline]
pub fn time_for_fraction(fraction_remaining: f64, decay_const: f64) -> Result<f64> {
    if fraction_remaining <= 0.0 || fraction_remaining > 1.0 {
        return Err(KimiyaError::InvalidInput(
            "fraction must be in (0, 1]".into(),
        ));
    }
    if decay_const <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "decay constant must be positive".into(),
        ));
    }
    Ok(-fraction_remaining.ln() / decay_const)
}

// ── Binding energy ───────────────────────────────────────────────────

/// Mass defect: Δm = Z·m_p + N·m_n - M_atom (in atomic mass units)
///
/// - `z`: number of protons
/// - `n`: number of neutrons
/// - `atomic_mass_u`: measured atomic mass (u)
#[must_use]
#[inline]
pub fn mass_defect(z: u32, n: u32, atomic_mass_u: f64) -> f64 {
    z as f64 * PROTON_MASS_U + n as f64 * NEUTRON_MASS_U - atomic_mass_u
}

/// Binding energy from mass defect: BE = Δm × 931.494 MeV
#[must_use]
#[inline]
pub fn binding_energy_mev(mass_defect_u: f64) -> f64 {
    mass_defect_u * AMU_MEV
}

/// Binding energy per nucleon: BE/A
///
/// # Errors
///
/// Returns error if mass number A is zero.
#[inline]
pub fn binding_energy_per_nucleon(z: u32, n: u32, atomic_mass_u: f64) -> Result<f64> {
    let a = z + n;
    if a == 0 {
        return Err(KimiyaError::InvalidInput(
            "mass number must be positive".into(),
        ));
    }
    if atomic_mass_u <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "atomic mass must be positive".into(),
        ));
    }
    let md = mass_defect(z, n, atomic_mass_u);
    Ok(binding_energy_mev(md) / a as f64)
}

/// Semi-empirical mass formula (Weizsacker / liquid drop model):
/// BE = a_v·A - a_s·A^(2/3) - a_c·Z(Z-1)/A^(1/3) - a_a·(A-2Z)²/A + δ(A,Z)
///
/// Returns binding energy in MeV.
#[must_use]
#[inline]
pub fn semi_empirical_binding_energy(z: u32, a: u32) -> f64 {
    let af = a as f64;
    let zf = z as f64;
    let n = (a - z) as f64;

    // Coefficients (MeV)
    let a_v = 15.67; // volume
    let a_s = 17.23; // surface
    let a_c = 0.714; // Coulomb
    let a_a = 93.15; // asymmetry (uses (N-Z)²/(4A) form: 93.15/4 ≈ 23.3)

    let volume = a_v * af;
    let surface = -a_s * af.powf(2.0 / 3.0);
    let coulomb = -a_c * zf * (zf - 1.0) / af.powf(1.0 / 3.0);
    let asymmetry = -a_a * (n - zf) * (n - zf) / (4.0 * af);

    // Pairing term
    let delta = if !a.is_multiple_of(2) {
        0.0 // odd A
    } else if z.is_multiple_of(2) {
        12.0 / af.sqrt() // even-even
    } else {
        -12.0 / af.sqrt() // odd-odd
    };

    volume + surface + coulomb + asymmetry + delta
}

/// Q-value of a nuclear reaction: Q = (Σ m_reactants - Σ m_products) × c²
///
/// - `reactant_masses_u`: atomic masses of reactants (u)
/// - `product_masses_u`: atomic masses of products (u)
///
/// Returns Q in MeV. Positive Q = exothermic (energy released).
#[must_use]
pub fn q_value(reactant_masses_u: &[f64], product_masses_u: &[f64]) -> f64 {
    let sum_r: f64 = reactant_masses_u.iter().sum();
    let sum_p: f64 = product_masses_u.iter().sum();
    (sum_r - sum_p) * AMU_MEV
}

// ── Isotope data ─────────────────────────────────────────────────────

/// A radioactive isotope with its half-life.
#[derive(Debug, Clone, Serialize)]
pub struct Isotope {
    pub symbol: &'static str,
    pub z: u8,
    pub a: u16,
    /// Half-life in seconds. Use f64::INFINITY for stable isotopes.
    pub half_life_s: f64,
    pub decay_mode: &'static str,
}

/// Built-in isotope data for common radioisotopes.
pub static ISOTOPES: &[Isotope] = &[
    Isotope {
        symbol: "H-3",
        z: 1,
        a: 3,
        half_life_s: 3.888e8,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "C-14",
        z: 6,
        a: 14,
        half_life_s: 1.808e11,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "K-40",
        z: 19,
        a: 40,
        half_life_s: 3.938e16,
        decay_mode: "β⁻/β⁺/EC",
    },
    Isotope {
        symbol: "Co-60",
        z: 27,
        a: 60,
        half_life_s: 1.663e8,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "Sr-90",
        z: 38,
        a: 90,
        half_life_s: 9.085e8,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "Tc-99m",
        z: 43,
        a: 99,
        half_life_s: 21_624.0,
        decay_mode: "IT",
    },
    Isotope {
        symbol: "I-131",
        z: 53,
        a: 131,
        half_life_s: 693_792.0,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "Cs-137",
        z: 55,
        a: 137,
        half_life_s: 9.496e8,
        decay_mode: "β⁻",
    },
    Isotope {
        symbol: "Ra-226",
        z: 88,
        a: 226,
        half_life_s: 5.049e10,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Rn-222",
        z: 86,
        a: 222,
        half_life_s: 330_350.0,
        decay_mode: "α",
    },
    Isotope {
        symbol: "U-235",
        z: 92,
        a: 235,
        half_life_s: 2.221e16,
        decay_mode: "α",
    },
    Isotope {
        symbol: "U-238",
        z: 92,
        a: 238,
        half_life_s: 1.409e17,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Pu-239",
        z: 94,
        a: 239,
        half_life_s: 7.594e11,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Am-241",
        z: 95,
        a: 241,
        half_life_s: 1.364e10,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Po-210",
        z: 84,
        a: 210,
        half_life_s: 1.196e7,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Th-232",
        z: 90,
        a: 232,
        half_life_s: 4.434e17,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Pa-231",
        z: 91,
        a: 231,
        half_life_s: 1.034e12,
        decay_mode: "α",
    },
    Isotope {
        symbol: "Np-237",
        z: 93,
        a: 237,
        half_life_s: 6.765e13,
        decay_mode: "α",
    },
];

/// Look up an isotope by symbol (e.g., "C-14", "U-235").
#[must_use]
#[inline]
pub fn lookup_isotope(symbol: &str) -> Option<&'static Isotope> {
    ISOTOPES.iter().find(|i| i.symbol == symbol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decay_constant_c14() {
        // C-14: t½ = 5730 years = 1.808e11 s → λ ≈ 3.83e-12 s⁻¹
        let lambda = decay_constant(1.808e11).unwrap();
        assert!(
            (lambda - 3.83e-12).abs() < 1e-13,
            "C-14 λ should be ~3.83e-12, got {lambda}"
        );
    }

    #[test]
    fn decay_constant_zero_is_error() {
        assert!(decay_constant(0.0).is_err());
    }

    #[test]
    fn atoms_remaining_one_half_life() {
        let lambda = decay_constant(100.0).unwrap();
        let n = atoms_remaining(1000.0, lambda, 100.0);
        assert!(
            (n - 500.0).abs() < 0.1,
            "should have half after one half-life, got {n}"
        );
    }

    #[test]
    fn atoms_remaining_two_half_lives() {
        let n = atoms_after_half_lives(1000.0, 2.0);
        assert!((n - 250.0).abs() < f64::EPSILON);
    }

    #[test]
    fn activity_basic() {
        let a = activity(0.01, 1e6);
        assert!((a - 1e4).abs() < f64::EPSILON);
    }

    #[test]
    fn time_for_half() {
        let lambda = decay_constant(100.0).unwrap();
        let t = time_for_fraction(0.5, lambda).unwrap();
        assert!((t - 100.0).abs() < 0.01);
    }

    #[test]
    fn time_for_fraction_invalid() {
        assert!(time_for_fraction(0.0, 0.01).is_err());
        assert!(time_for_fraction(1.5, 0.01).is_err());
        assert!(time_for_fraction(0.5, 0.0).is_err());
    }

    #[test]
    fn mass_defect_he4() {
        // He-4: Z=2, N=2, M=4.002602 u
        // Δm = 2×1.007276 + 2×1.008665 - 4.002602 = 0.03028 u
        let md = mass_defect(2, 2, 4.002602);
        assert!(
            (md - 0.0293).abs() < 0.002,
            "He-4 mass defect should be ~0.029 u, got {md}"
        );
    }

    #[test]
    fn binding_energy_he4() {
        let md = mass_defect(2, 2, 4.002602);
        let be = binding_energy_mev(md);
        // ~27.3 MeV (using simplified nucleon masses)
        assert!(
            (be - 27.3).abs() < 1.5,
            "He-4 BE should be ~27 MeV, got {be}"
        );
    }

    #[test]
    fn be_per_nucleon_fe56() {
        // Fe-56: Z=26, N=30, M=55.9349 u, BE/A ≈ 8.79 MeV
        let bea = binding_energy_per_nucleon(26, 30, 55.9349).unwrap();
        assert!(
            (bea - 8.55).abs() < 0.5,
            "Fe-56 BE/A should be ~8.5-8.8 MeV, got {bea}"
        );
    }

    #[test]
    fn semi_empirical_fe56() {
        // SEMF should give ~490 MeV for Fe-56 (actual ~492 MeV)
        let be = semi_empirical_binding_energy(26, 56);
        assert!(
            (be - 492.0).abs() < 15.0,
            "SEMF Fe-56 should be ~492 MeV, got {be}"
        );
    }

    #[test]
    fn semi_empirical_he4() {
        let be = semi_empirical_binding_energy(2, 4);
        // SEMF is less accurate for light nuclei, but should be > 0
        assert!(be > 0.0, "He-4 BE should be positive, got {be}");
    }

    #[test]
    fn lookup_c14() {
        let iso = lookup_isotope("C-14").unwrap();
        assert_eq!(iso.z, 6);
        assert_eq!(iso.a, 14);
        assert!(iso.half_life_s > 1e11);
    }

    #[test]
    fn lookup_u235() {
        let iso = lookup_isotope("U-235").unwrap();
        assert_eq!(iso.z, 92);
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_isotope("Xx-999").is_none());
    }

    #[test]
    fn isotope_count() {
        assert!(ISOTOPES.len() >= 18);
    }

    // ── Q-value ──────────────────────────────────────────────────────

    #[test]
    fn q_value_deuterium_tritium_fusion() {
        // D + T → He-4 + n
        // m_D=2.01410, m_T=3.01605, m_He4=4.00260, m_n=1.00866
        let q = q_value(&[2.01410, 3.01605], &[4.00260, 1.00866]);
        // Q ≈ 17.6 MeV
        assert!(
            (q - 17.6).abs() < 0.5,
            "D-T fusion Q should be ~17.6 MeV, got {q}"
        );
    }

    #[test]
    fn q_value_endothermic_is_negative() {
        // Reverse reaction has negative Q
        let q = q_value(&[4.00260, 1.00866], &[2.01410, 3.01605]);
        assert!(q < 0.0, "reverse fusion should be endothermic");
    }
}
