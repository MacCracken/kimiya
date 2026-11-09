//! Nuclear chemistry — radioactive decay, binding energy, isotope data.

use crate::error::{KimiyaError, Result};
use crate::validate::{is_non_negative, is_positive};
use serde::Serialize;

// ── Constants ────────────────────────────────────────────────────────

/// Atomic mass unit energy equivalent m_u·c² in MeV (CODATA 2018).
pub const AMU_MEV: f64 = 931.494_102_42;

/// Proton mass in atomic mass units, u (CODATA 2018).
pub const PROTON_MASS_U: f64 = 1.007_276_466_621;

/// Neutron mass in atomic mass units, u (CODATA 2018).
pub const NEUTRON_MASS_U: f64 = 1.008_664_915_95;

/// Mass of the neutral ¹H atom (proton + electron, less the 13.6 eV binding
/// energy) in atomic mass units, u (AME2020).
///
/// Used with neutral-atom masses in [`mass_defect`] so the electron masses cancel.
pub const HYDROGEN_ATOM_MASS_U: f64 = 1.007_825_031_898;

/// Relative tolerance below which two decay constants are treated as equal
/// in [`decay_chain_two_step`].
const EQUAL_RATE_REL_TOL: f64 = 1e-12;

// Semi-empirical mass formula coefficients (MeV), Povh, Rith, Scholz &
// Zetsche, *Particles and Nuclei*, eq. 2.8 — fitted with a Z² Coulomb term,
// an (N − Z)²/(4A) asymmetry term and a pairing term δ/√A.
/// SEMF volume coefficient a_v (MeV).
const SEMF_VOLUME: f64 = 15.67;
/// SEMF surface coefficient a_s (MeV).
const SEMF_SURFACE: f64 = 17.23;
/// SEMF Coulomb coefficient a_c (MeV), multiplies Z²/A^(1/3).
const SEMF_COULOMB: f64 = 0.714;
/// SEMF asymmetry coefficient a_a (MeV), multiplies (N − Z)²/(4A).
const SEMF_ASYMMETRY: f64 = 93.15;
/// SEMF pairing coefficient δ (MeV), multiplies 1/√A.
const SEMF_PAIRING: f64 = 11.2;

// ── Radioactive decay ────────────────────────────────────────────────

/// Decay constant from half-life: λ = ln(2) / t½
///
/// `half_life` is in any time unit; λ is returned in the reciprocal unit.
/// A half-life of `f64::INFINITY` (the stable-isotope convention used by
/// [`Isotope::half_life_s`]) gives λ = 0.
///
/// # Errors
///
/// Returns error if half-life is not positive (zero, negative, NaN or −∞).
#[inline]
pub fn decay_constant(half_life: f64) -> Result<f64> {
    let stable = half_life == f64::INFINITY;
    if !is_positive(half_life) && !stable {
        return Err(KimiyaError::InvalidInput(
            "half-life must be positive".into(),
        ));
    }
    Ok(std::f64::consts::LN_2 / half_life)
}

/// Number of atoms remaining: N(t) = N₀ × exp(-λt)
///
/// `decay_const` and `time` must use reciprocal units (e.g. s⁻¹ and s).
#[must_use]
#[inline]
pub fn atoms_remaining(initial: f64, decay_const: f64, time: f64) -> f64 {
    initial * (-decay_const * time).exp()
}

/// Activity: A = λN (decays per unit time; Bq when λ is in s⁻¹)
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
/// Returns time in the reciprocal unit of `decay_const`.
///
/// # Errors
///
/// Returns error if fraction is not in (0, 1] or decay constant is not
/// positive and finite (NaN is rejected).
#[inline]
pub fn time_for_fraction(fraction_remaining: f64, decay_const: f64) -> Result<f64> {
    if !is_positive(fraction_remaining) || fraction_remaining > 1.0 {
        return Err(KimiyaError::InvalidInput(
            "fraction must be in (0, 1]".into(),
        ));
    }
    if !is_positive(decay_const) {
        return Err(KimiyaError::InvalidInput(
            "decay constant must be positive".into(),
        ));
    }
    Ok(-fraction_remaining.ln() / decay_const)
}

// ── Binding energy ───────────────────────────────────────────────────

/// Mass defect: Δm = Z·m(¹H) + N·m_n − M_atom (in atomic mass units, u)
///
/// - `z`: number of protons
/// - `n`: number of neutrons
/// - `atomic_mass_u`: **neutral-atom** mass (u), as tabulated in AME2020
///
/// Uses the ¹H atom mass ([`HYDROGEN_ATOM_MASS_U`]) rather than the bare
/// proton mass, so the Z electrons included in `atomic_mass_u` cancel
/// (electron binding energies, at most a few hundred keV in total, are
/// neglected).
#[must_use]
#[inline]
pub fn mass_defect(z: u32, n: u32, atomic_mass_u: f64) -> f64 {
    z as f64 * HYDROGEN_ATOM_MASS_U + n as f64 * NEUTRON_MASS_U - atomic_mass_u
}

/// Binding energy from mass defect: BE = Δm × m_u·c² ([`AMU_MEV`])
///
/// `mass_defect_u` is in u; returns MeV.
#[must_use]
#[inline]
pub fn binding_energy_mev(mass_defect_u: f64) -> f64 {
    mass_defect_u * AMU_MEV
}

/// Binding energy per nucleon: BE/A (MeV)
///
/// `atomic_mass_u` is the neutral-atom mass (u); see [`mass_defect`].
///
/// # Errors
///
/// Returns error if mass number A = Z + N is zero or overflows `u32`, or if
/// the atomic mass is not positive and finite.
#[inline]
pub fn binding_energy_per_nucleon(z: u32, n: u32, atomic_mass_u: f64) -> Result<f64> {
    let a = z
        .checked_add(n)
        .ok_or_else(|| KimiyaError::InvalidInput("mass number overflows u32".into()))?;
    if a == 0 {
        return Err(KimiyaError::InvalidInput(
            "mass number must be positive".into(),
        ));
    }
    if !is_positive(atomic_mass_u) {
        return Err(KimiyaError::InvalidInput(
            "atomic mass must be positive".into(),
        ));
    }
    let md = mass_defect(z, n, atomic_mass_u);
    Ok(binding_energy_mev(md) / a as f64)
}

/// Semi-empirical mass formula (Weizsäcker / liquid drop model):
/// BE = a_v·A − a_s·A^(2/3) − a_c·Z²/A^(1/3) − a_a·(N − Z)²/(4A) + δ/√A
///
/// Coefficients are Povh et al. (*Particles and Nuclei*): a_v = 15.67,
/// a_s = 17.23, a_c = 0.714, a_a = 93.15, δ = +11.2 (even-even),
/// 0 (odd A), −11.2 (odd-odd), all MeV. Accuracy is about 0.4 % for
/// medium and heavy nuclei (Fe-56: 494.0 vs 492.26 MeV; U-238: 1808.0 vs
/// 1801.69 MeV) and poor for very light ones.
///
/// Returns binding energy in MeV, or `f64::NAN` if `a == 0` or `z > a`
/// (no such nucleus).
#[must_use]
#[inline]
pub fn semi_empirical_binding_energy(z: u32, a: u32) -> f64 {
    let Some(n) = a.checked_sub(z) else {
        return f64::NAN;
    };
    if a == 0 {
        return f64::NAN;
    }
    let af = a as f64;
    let zf = z as f64;
    let asym = n as f64 - zf;
    let cbrt_a = af.cbrt();

    let volume = SEMF_VOLUME * af;
    let surface = SEMF_SURFACE * cbrt_a * cbrt_a;
    let coulomb = SEMF_COULOMB * zf * zf / cbrt_a;
    let asymmetry = SEMF_ASYMMETRY * asym * asym / (4.0 * af);

    let pairing = if !a.is_multiple_of(2) {
        0.0 // odd A
    } else if z.is_multiple_of(2) {
        SEMF_PAIRING / af.sqrt() // even-even
    } else {
        -SEMF_PAIRING / af.sqrt() // odd-odd
    };

    volume - surface - coulomb - asymmetry + pairing
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

// ── Decay chains ─────────────────────────────────────────────────────

/// Two-step decay chain: A →(λ₁) B →(λ₂) C (Bateman equation)
///
/// Returns (N_A(t), N_B(t), N_C(t)) given initial N_A₀ = `n0` (and no B or C
/// at t = 0). Decay constants and `time` must use reciprocal units.
///
/// N_B(t) = N₀ × λ₁/(λ₂-λ₁) × (e^(-λ₁t) - e^(-λ₂t))
///
/// evaluated in a cancellation-free form. When λ₁ and λ₂ agree to within a
/// relative 10⁻¹², the analytic limit N_B(t) = N₀·λ·t·e^(−λt) is used.
///
/// # Errors
///
/// Returns error if a decay constant is not positive and finite, or if `n0`
/// or `time` is negative or not finite.
pub fn decay_chain_two_step(
    n0: f64,
    lambda1: f64,
    lambda2: f64,
    time: f64,
) -> Result<(f64, f64, f64)> {
    if !is_positive(lambda1) || !is_positive(lambda2) {
        return Err(KimiyaError::InvalidInput(
            "decay constants must be positive".into(),
        ));
    }
    if !is_non_negative(n0) {
        return Err(KimiyaError::InvalidInput(
            "initial amount must be non-negative".into(),
        ));
    }
    if !is_non_negative(time) {
        return Err(KimiyaError::InvalidInput(
            "time must be non-negative".into(),
        ));
    }
    let na = n0 * (-lambda1 * time).exp();

    // (e^(−λ₁t) − e^(−λ₂t)) / (λ₂ − λ₁) = e^(−λ_min·t) · (1 − e^(−|Δλ|t)) / |Δλ|
    // The right-hand side has no catastrophic cancellation and tends to
    // e^(−λt)·t as Δλ → 0.
    let delta = (lambda2 - lambda1).abs();
    let lambda_min = lambda1.min(lambda2);
    let growth = if delta <= EQUAL_RATE_REL_TOL * lambda1.max(lambda2) {
        time
    } else {
        -(-delta * time).exp_m1() / delta
    };
    let nb = n0 * lambda1 * (-lambda_min * time).exp() * growth;
    // N_C = N₀(1 − e^(−λ₁t)) − N_B, with expm1 for accuracy at small t.
    let nc = -n0 * (-lambda1 * time).exp_m1() - nb;
    Ok((na, nb, nc))
}

/// Secular equilibrium condition: when λ₁ << λ₂ (parent much longer-lived than daughter).
///
/// At secular equilibrium: A_daughter = A_parent (activities equal).
///
/// Returns true if λ₁ < 0.01 × λ₂.
#[must_use]
#[inline]
pub fn is_secular_equilibrium(lambda_parent: f64, lambda_daughter: f64) -> bool {
    lambda_parent < 0.01 * lambda_daughter
}

// ── Isotope data ─────────────────────────────────────────────────────

/// A radioactive isotope with its half-life.
#[derive(Debug, Clone, Serialize)]
pub struct Isotope {
    /// Nuclide label, element symbol and mass number (e.g. `"C-14"`,
    /// `"Tc-99m"` for a metastable state).
    pub symbol: &'static str,
    /// Atomic number Z (protons).
    pub z: u8,
    /// Mass number A (protons + neutrons).
    pub a: u16,
    /// Half-life in seconds. Use f64::INFINITY for stable isotopes.
    pub half_life_s: f64,
    /// Principal decay mode(s): `"α"`, `"β⁻"`, `"β⁺"`, `"EC"` (electron
    /// capture), `"IT"` (isomeric transition); multiple modes are `/`-separated.
    pub decay_mode: &'static str,
}

/// Built-in isotope data for common radioisotopes.
///
/// Half-lives are in seconds (NUBASE2020 / ENSDF evaluations to 3–5
/// significant figures), converted with the Julian year
/// (1 a = 365.25 d = 3.15576 × 10⁷ s).
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
        half_life_s: 1.7988e11, // 5700 a (NUBASE2020)
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
        half_life_s: 693_377.28, // 8.0252 d (NUBASE2020)
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
        half_life_s: 7.6085e11, // 24 110 a (NUBASE2020)
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
        // C-14: t½ = 5700 a (NUBASE2020) = 1.7988e11 s → λ = 3.8534e-12 s⁻¹
        let iso = lookup_isotope("C-14").unwrap();
        let lambda = decay_constant(iso.half_life_s).unwrap();
        assert!(
            (lambda - 3.8534e-12).abs() < 1e-16,
            "C-14 λ should be 3.8534e-12, got {lambda}"
        );
    }

    #[test]
    fn decay_constant_zero_is_error() {
        assert!(decay_constant(0.0).is_err());
    }

    #[test]
    fn decay_constant_nan_is_error() {
        assert!(decay_constant(f64::NAN).is_err());
        assert!(decay_constant(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn decay_constant_stable_is_zero() {
        // f64::INFINITY is the documented stable-isotope half-life.
        assert_eq!(decay_constant(f64::INFINITY).unwrap(), 0.0);
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
    fn time_for_fraction_nan_is_error() {
        assert!(time_for_fraction(f64::NAN, 0.01).is_err());
        assert!(time_for_fraction(0.5, f64::NAN).is_err());
        assert!(time_for_fraction(0.5, f64::INFINITY).is_err());
    }

    // ── Constants ────────────────────────────────────────────────────

    #[test]
    fn constants_codata_2018() {
        assert_eq!(AMU_MEV, 931.494_102_42);
        assert_eq!(PROTON_MASS_U, 1.007_276_466_621);
        assert_eq!(NEUTRON_MASS_U, 1.008_664_915_95);
        assert_eq!(HYDROGEN_ATOM_MASS_U, 1.007_825_031_898); // AME2020
        // Cross-check: CODATA 2018 m_p·c² = 938.272 088 16 MeV,
        // m_n·c² = 939.565 420 52 MeV.
        assert!((PROTON_MASS_U * AMU_MEV - 938.272_088_16).abs() < 1e-8);
        assert!((NEUTRON_MASS_U * AMU_MEV - 939.565_420_52).abs() < 1e-8);
    }

    #[test]
    fn hydrogen_atom_mass_is_proton_plus_electron() {
        // m(¹H) = m_p + m_e − 13.6 eV/c²; CODATA 2018 m_e = 5.485 799 090 65e-4 u.
        // AME2020 and CODATA 2018 agree to ~3e-11 u (within their uncertainties).
        let m_e = 5.485_799_090_65e-4;
        let binding_u = 13.598_434e-6 / AMU_MEV;
        assert!((HYDROGEN_ATOM_MASS_U - (PROTON_MASS_U + m_e - binding_u)).abs() < 1e-10);
    }

    #[test]
    fn mass_defect_he4() {
        // He-4: Z=2, N=2, neutral-atom M = 4.002 603 254 u (AME2020)
        // Δm = 2×1.007 825 032 + 2×1.008 664 916 − 4.002 603 254 = 0.030 377 u
        let md = mass_defect(2, 2, 4.002_603_254);
        assert!(
            (md - 0.030_377).abs() < 1e-5,
            "He-4 mass defect should be 0.030377 u, got {md}"
        );
    }

    #[test]
    fn binding_energy_he4() {
        // He-4 BE = 28.296 MeV (AME2020)
        let md = mass_defect(2, 2, 4.002_603_254);
        let be = binding_energy_mev(md);
        assert!(
            (be - 28.296).abs() < 0.01,
            "He-4 BE should be 28.296 MeV, got {be}"
        );
    }

    #[test]
    fn be_per_nucleon_fe56() {
        // Fe-56: Z=26, N=30, M = 55.934 936 33 u (AME2020), BE/A = 8.790 MeV
        let bea = binding_energy_per_nucleon(26, 30, 55.934_936_33).unwrap();
        assert!(
            (bea - 8.790).abs() < 0.01,
            "Fe-56 BE/A should be 8.790 MeV, got {bea}"
        );
    }

    #[test]
    fn be_per_nucleon_overflow_is_error() {
        assert!(binding_energy_per_nucleon(u32::MAX, 1, 1.0).is_err());
    }

    #[test]
    fn be_per_nucleon_nan_mass_is_error() {
        assert!(binding_energy_per_nucleon(26, 30, f64::NAN).is_err());
        assert!(binding_energy_per_nucleon(26, 30, f64::INFINITY).is_err());
    }

    #[test]
    fn semi_empirical_fe56() {
        // Experimental BE(Fe-56) = 492.26 MeV (AME2020); Povh SEMF gives 494.0.
        let be = semi_empirical_binding_energy(26, 56);
        assert!(
            (be - 492.26).abs() < 2.0,
            "SEMF Fe-56 should be ~492.3 MeV, got {be}"
        );
    }

    #[test]
    fn semi_empirical_u238() {
        // Experimental BE(U-238) = 1801.69 MeV (AME2020); Povh SEMF gives 1808.0.
        let be = semi_empirical_binding_energy(92, 238);
        assert!(
            (be - 1801.69).abs() < 7.0,
            "SEMF U-238 should be ~1801.7 MeV, got {be}"
        );
    }

    #[test]
    fn semi_empirical_invalid_is_nan() {
        // z > a used to panic (debug) / wrap (release); a = 0 gave NaN by accident.
        assert!(semi_empirical_binding_energy(3, 2).is_nan());
        assert!(semi_empirical_binding_energy(0, 0).is_nan());
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
    fn half_lives_nubase2020() {
        let julian_year = 3.155_76e7;
        let i131 = lookup_isotope("I-131").unwrap().half_life_s;
        assert!((i131 - 8.0252 * 86_400.0).abs() < 1.0, "I-131: {i131}");
        let pu239 = lookup_isotope("Pu-239").unwrap().half_life_s;
        assert!(
            (pu239 / (24_110.0 * julian_year) - 1.0).abs() < 1e-4,
            "Pu-239: {pu239}"
        );
        let c14 = lookup_isotope("C-14").unwrap().half_life_s;
        assert!(
            (c14 / (5_700.0 * julian_year) - 1.0).abs() < 1e-4,
            "C-14: {c14}"
        );
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
        let q = q_value(&[4.00260, 1.00866], &[2.01410, 3.01605]);
        assert!(q < 0.0, "reverse fusion should be endothermic");
    }

    // ── Decay chains ─────────────────────────────────────────────────

    #[test]
    fn decay_chain_mass_conservation() {
        let (na, nb, nc) = decay_chain_two_step(1000.0, 0.1, 0.05, 10.0).unwrap();
        let total = na + nb + nc;
        assert!(
            (total - 1000.0).abs() < 0.1,
            "mass should be conserved: {total}"
        );
    }

    #[test]
    fn decay_chain_at_t0() {
        let (na, nb, nc) = decay_chain_two_step(1000.0, 0.1, 0.05, 0.0).unwrap();
        assert!((na - 1000.0).abs() < f64::EPSILON);
        assert!(nb.abs() < f64::EPSILON);
        assert!(nc.abs() < f64::EPSILON);
    }

    #[test]
    fn decay_chain_long_time() {
        // After very long time, almost all should be C
        let (na, nb, nc) = decay_chain_two_step(1000.0, 0.1, 0.05, 200.0).unwrap();
        assert!(nc > 990.0, "most should be C after long time, got {nc}");
        assert!(na < 1.0);
        assert!(nb < 10.0);
    }

    #[test]
    fn decay_chain_equal_rates_uses_limit() {
        // λ₁ = λ₂ = λ: N_B = N₀·λ·t·e^(−λt) = 1000 × 0.1 × 10 × e⁻¹
        let (na, nb, nc) = decay_chain_two_step(1000.0, 0.1, 0.1, 10.0).unwrap();
        let expected = 1000.0 * (-1.0_f64).exp();
        assert!((nb - expected).abs() < 1e-9, "N_B = {nb}, want {expected}");
        assert!((na - expected).abs() < 1e-9);
        assert!((na + nb + nc - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn decay_chain_nearly_equal_rates_is_accurate() {
        // Δλ/λ = 1e-10 is outside the equal-rate tolerance; the naive Bateman
        // form loses ~6 digits here to cancellation.
        let (_, nb, _) = decay_chain_two_step(1000.0, 0.1, 0.1 * (1.0 + 1e-10), 10.0).unwrap();
        let limit = 1000.0 * (-1.0_f64).exp();
        assert!(
            ((nb - limit) / limit).abs() < 1e-9,
            "N_B = {nb}, limit {limit}"
        );
        // λ₂ < λ₁ branch, far from equal: compare with the textbook formula.
        let (l1, l2, t) = (2.0, 0.5, 3.0);
        let (_, nb, _) = decay_chain_two_step(1.0, l1, l2, t).unwrap();
        let textbook = l1 / (l2 - l1) * ((-l1 * t).exp() - (-l2 * t).exp());
        assert!((nb - textbook).abs() < 1e-14, "{nb} vs {textbook}");
    }

    #[test]
    fn decay_chain_no_overflow_for_fast_parent() {
        // λ₁ ≫ λ₂ at long times must not produce inf·0 = NaN.
        let (na, nb, nc) = decay_chain_two_step(1.0, 10.0, 0.1, 100.0).unwrap();
        assert!(na == 0.0 && nb.is_finite() && nc.is_finite());
        let textbook = 10.0 / (0.1 - 10.0) * (0.0 - (-10.0_f64).exp());
        assert!((nb - textbook).abs() < 1e-15);
    }

    #[test]
    fn decay_chain_invalid_inputs_are_errors() {
        assert!(decay_chain_two_step(1000.0, 0.1, 0.05, -1.0).is_err());
        assert!(decay_chain_two_step(-1.0, 0.1, 0.05, 1.0).is_err());
        assert!(decay_chain_two_step(f64::NAN, 0.1, 0.05, 1.0).is_err());
        assert!(decay_chain_two_step(1000.0, f64::NAN, 0.05, 1.0).is_err());
        assert!(decay_chain_two_step(1000.0, 0.1, f64::INFINITY, 1.0).is_err());
        assert!(decay_chain_two_step(1000.0, 0.1, 0.05, f64::NAN).is_err());
    }

    #[test]
    fn secular_equilibrium_check() {
        // U-238 (λ ~ 5e-18) → Th-234 (λ ~ 3.3e-7) → very secular
        assert!(is_secular_equilibrium(5e-18, 3.3e-7));
        // Two similar rates → not secular
        assert!(!is_secular_equilibrium(0.1, 0.2));
    }
}
