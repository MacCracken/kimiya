//! Inorganic chemistry — crystal field theory, lattice energy, ionic radii.

use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_positive};
use serde::Serialize;

// ── Crystal field theory ─────────────────────────────────────────────

/// Octahedral (t₂g, e_g) occupancies for high-spin dⁿ, n = 0..=10
/// (Hund's rule: singly occupy all five d orbitals before pairing).
const OCT_OCCUPANCY_HIGH_SPIN: [(u8, u8); 11] = [
    (0, 0),
    (1, 0),
    (2, 0),
    (3, 0),
    (3, 1),
    (3, 2),
    (4, 2),
    (5, 2),
    (6, 2),
    (6, 3),
    (6, 4),
];

/// Octahedral (t₂g, e_g) occupancies for low-spin dⁿ, n = 0..=10
/// (t₂g filled completely before e_g).
const OCT_OCCUPANCY_LOW_SPIN: [(u8, u8); 11] = [
    (0, 0),
    (1, 0),
    (2, 0),
    (3, 0),
    (4, 0),
    (5, 0),
    (6, 0),
    (6, 1),
    (6, 2),
    (6, 3),
    (6, 4),
];

/// Unpaired electrons for high-spin octahedral dⁿ, n = 0..=10.
const OCT_UNPAIRED_HIGH_SPIN: [u8; 11] = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0];

/// Unpaired electrons for low-spin octahedral dⁿ, n = 0..=10.
const OCT_UNPAIRED_LOW_SPIN: [u8; 11] = [0, 1, 2, 3, 2, 1, 0, 1, 2, 1, 0];

/// Error for a d-electron count outside 0–10.
fn invalid_d_count() -> KimiyaError {
    KimiyaError::InvalidInput("d electron count must be 0-10".into())
}

/// Crystal field stabilization energy (CFSE) of an octahedral dⁿ complex,
/// in units of Δ_oct (dimensionless multiplier; negative = stabilizing).
///
/// CFSE = −0.4·n(t₂g) + 0.6·n(e_g). Electron-pairing energy is **not**
/// included, so for d⁴–d⁷ the low-spin value is more negative than the
/// high-spin one even when the high-spin state is the ground state.
///
/// - `d_electrons`: number of d electrons (0–10)
/// - `high_spin`: true for high-spin (weak field), false for low-spin (strong field)
///
/// # Errors
///
/// Returns error if d_electrons > 10.
pub fn cfse_octahedral(d_electrons: u8, high_spin: bool) -> Result<f64> {
    let table = if high_spin {
        &OCT_OCCUPANCY_HIGH_SPIN
    } else {
        &OCT_OCCUPANCY_LOW_SPIN
    };
    let &(t2g, eg) = table
        .get(usize::from(d_electrons))
        .ok_or_else(invalid_d_count)?;
    Ok(f64::from(t2g) * (-0.4) + f64::from(eg) * 0.6)
}

/// Crystal field splitting energy Δ_tet for tetrahedral complexes, in the
/// same units as `delta_oct`.
///
/// Δ_tet ≈ (4/9) × Δ_oct (same metal, ligands and metal–ligand distance).
#[must_use]
#[inline]
pub fn delta_tet_from_oct(delta_oct: f64) -> f64 {
    delta_oct * 4.0 / 9.0
}

/// Spectrochemical series — common ligands, ordered from weakest to
/// strongest field.
///
/// [`Ligand::field_strength`] gives a dimensionless ordinal rank (1–12), not
/// a Δ_oct value; only the ordering is meaningful.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[non_exhaustive]
pub enum Ligand {
    /// I⁻ — weakest-field ligand in this list.
    Iodide,
    /// Br⁻.
    Bromide,
    /// Cl⁻.
    Chloride,
    /// F⁻.
    Fluoride,
    /// OH⁻.
    Hydroxide,
    /// H₂O.
    Water,
    /// NH₃ — intermediate/borderline field strength.
    Ammonia,
    /// Ethylenediamine (en), bidentate.
    Ethylenediamine,
    /// NO₂⁻ (N-bound nitro).
    Nitrite,
    /// CN⁻.
    Cyanide,
    /// CO — strongest-field ligand in this list.
    CarbonMonoxide,
}

impl Ligand {
    /// Relative field strength as an ordinal rank 1–12 (dimensionless;
    /// higher = stronger field). Gaps between ranks carry no meaning.
    #[must_use]
    #[inline]
    pub const fn field_strength(self) -> u8 {
        match self {
            Self::Iodide => 1,
            Self::Bromide => 2,
            Self::Chloride => 3,
            Self::Fluoride => 4,
            Self::Hydroxide => 5,
            Self::Water => 6,
            Self::Ammonia => 8,
            Self::Ethylenediamine => 9,
            Self::Nitrite => 10,
            Self::Cyanide => 11,
            Self::CarbonMonoxide => 12,
        }
    }

    /// Whether this is a strong-field ligand (typically causes low-spin),
    /// i.e. rank ≥ 8 ([`Ligand::Ammonia`] and stronger).
    ///
    /// NH₃ is borderline/intermediate: whether its complexes are low-spin
    /// depends on the metal ion and its oxidation state (e.g. \[Co(NH₃)₆\]³⁺
    /// is low-spin, \[Co(NH₃)₆\]²⁺ is high-spin).
    #[must_use]
    pub const fn is_strong_field(self) -> bool {
        self.field_strength() >= 8
    }
}

// ── Spin-only magnetic moment ────────────────────────────────────────

/// Spin-only magnetic moment in Bohr magnetons (μ_B): μ = √(n(n+2))
///
/// - `unpaired_electrons`: n
#[must_use]
#[inline]
pub fn spin_only_moment(unpaired_electrons: u8) -> f64 {
    let n = unpaired_electrons as f64;
    (n * (n + 2.0)).sqrt()
}

/// Number of unpaired electrons for octahedral dⁿ.
///
/// - `d_electrons`: number of d electrons (0–10)
/// - `high_spin`: true for high-spin (weak field), false for low-spin (strong field)
///
/// # Errors
///
/// Returns error if d_electrons > 10.
pub fn unpaired_electrons_oct(d_electrons: u8, high_spin: bool) -> Result<u8> {
    let table = if high_spin {
        &OCT_UNPAIRED_HIGH_SPIN
    } else {
        &OCT_UNPAIRED_LOW_SPIN
    };
    table
        .get(usize::from(d_electrons))
        .copied()
        .ok_or_else(invalid_d_count)
}

// ── Lattice energy ───────────────────────────────────────────────────

/// Born-Landé equation for lattice energy:
/// U = -(N_A × M × e² × z⁺ × z⁻) / (4πε₀ × r₀) × (1 - 1/n)
///
/// Simplified form using Coulomb constant:
/// U = -(k_e × N_A × M × z⁺ × z⁻ × e²) / r₀ × (1 - 1/n)
///
/// - `madelung`: Madelung constant M (dimensionless, positive)
/// - `z_cation`: cation charge magnitude (positive integer)
/// - `z_anion`: anion charge magnitude (positive integer, sign handled internally)
/// - `r0_m`: nearest-neighbor distance in meters
/// - `born_exponent`: Born exponent n (dimensionless, typically 5–12)
///
/// Returns lattice energy in J/mol (negative = stable).
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `r0_m` or `madelung` is not
/// positive and finite, if `born_exponent` is not finite and > 1, or if
/// either charge is zero.
pub fn born_lande_lattice_energy(
    madelung: f64,
    z_cation: u32,
    z_anion: u32,
    r0_m: f64,
    born_exponent: f64,
) -> Result<f64> {
    if !is_positive(r0_m) {
        return Err(KimiyaError::InvalidInput(
            "interionic distance must be positive".into(),
        ));
    }
    if !is_finite(born_exponent) || born_exponent <= 1.0 {
        return Err(KimiyaError::InvalidInput(
            "Born exponent must be > 1".into(),
        ));
    }
    if !is_positive(madelung) {
        return Err(KimiyaError::InvalidInput(
            "Madelung constant must be positive".into(),
        ));
    }
    if z_cation == 0 || z_anion == 0 {
        return Err(KimiyaError::InvalidInput(
            "ion charges must be non-zero".into(),
        ));
    }
    let na = crate::element::AVOGADRO;
    let e = crate::electrochemistry::ELEMENTARY_CHARGE;
    let ke = crate::potential::COULOMB_CONSTANT;

    Ok(
        -ke * na * madelung * f64::from(z_cation) * f64::from(z_anion) * e * e / r0_m
            * (1.0 - 1.0 / born_exponent),
    )
}

/// Born-Haber cycle: calculate lattice energy from thermodynamic cycle.
///
/// ΔH_f° = ΔH_sub + IE + ΔH_diss/2 + EA + U_lattice
///
/// Rearranged: U_lattice = ΔH_f° - ΔH_sub - IE - ΔH_diss/2 - EA
///
/// - `formation_enthalpy_kj`: ΔH_f° of the ionic compound (kJ/mol)
/// - `sublimation_enthalpy_kj`: ΔH_sub of the metal (kJ/mol)
/// - `ionization_energy_kj`: IE of the metal → cation (kJ/mol)
/// - `dissociation_enthalpy_kj`: bond dissociation of X₂ → 2X (kJ/mol, full value — halved internally for MX)
/// - `electron_affinity_kj`: EA of the nonmetal (kJ/mol, negative = energy released)
///
/// Returns lattice energy in kJ/mol (negative = stable).
#[must_use]
pub fn born_haber_lattice_energy(
    formation_enthalpy_kj: f64,
    sublimation_enthalpy_kj: f64,
    ionization_energy_kj: f64,
    dissociation_enthalpy_kj: f64,
    electron_affinity_kj: f64,
) -> f64 {
    formation_enthalpy_kj
        - sublimation_enthalpy_kj
        - ionization_energy_kj
        - dissociation_enthalpy_kj / 2.0
        - electron_affinity_kj
}

// Madelung constants (dimensionless), referred to the nearest cation–anion
// distance r₀ and used with |z₊·z₋| in [`born_lande_lattice_energy`].

/// Madelung constant of the rock-salt (NaCl) structure, 6:6 coordination.
pub const MADELUNG_NACL: f64 = 1.747565;
/// Madelung constant of the caesium chloride (CsCl) structure, 8:8 coordination.
pub const MADELUNG_CSCL: f64 = 1.762675;
/// Madelung constant of the zinc blende (sphalerite, cubic ZnS) structure,
/// 4:4 coordination (wurtzite: 1.6413).
pub const MADELUNG_ZNS: f64 = 1.6381;
/// Madelung constant of the fluorite (CaF₂) structure, 8:4 coordination.
pub const MADELUNG_FLUORITE: f64 = 2.51939;
/// Madelung constant of the rutile (TiO₂) structure, 6:3 coordination.
pub const MADELUNG_RUTILE: f64 = 2.408;

// Born exponents n (dimensionless) by the ion's noble-gas electron
// configuration (Pauling). For a salt whose ions have different
// configurations use the mean, e.g. NaCl: (7 + 9) / 2 = 8.

/// Born exponent for ions with the He configuration (Li⁺, Be²⁺, H⁻).
pub const BORN_EXP_HE: f64 = 5.0;
/// Born exponent for ions with the Ne configuration (Na⁺, Mg²⁺, F⁻, O²⁻).
pub const BORN_EXP_NE: f64 = 7.0;
/// Born exponent for ions with the Ar configuration or \[Ar\]3d¹⁰ (K⁺, Ca²⁺,
/// Cl⁻, S²⁻, Cu⁺).
pub const BORN_EXP_AR: f64 = 9.0;
/// Born exponent for ions with the Kr configuration or \[Kr\]4d¹⁰ (Rb⁺, Sr²⁺,
/// Br⁻, Ag⁺).
pub const BORN_EXP_KR: f64 = 10.0;
/// Born exponent for ions with the Xe configuration or \[Xe\]4f¹⁴5d¹⁰ (Cs⁺,
/// Ba²⁺, I⁻, Au⁺).
pub const BORN_EXP_XE: f64 = 12.0;

// ── Ionic radii ──────────────────────────────────────────────────────

/// Shannon (1976) effective ionic radius entry.
#[derive(Debug, Clone, Serialize)]
pub struct IonicRadius {
    /// Element symbol, e.g. `"Fe"`.
    pub symbol: &'static str,
    /// Ion charge in units of e (negative for anions).
    pub charge: i8,
    /// Coordination number the radius applies to.
    pub coordination: u8,
    /// Effective ionic radius in pm.
    pub radius_pm: f64,
}

/// Common Shannon (1976) effective ionic radii in pm, one entry per ion.
///
/// All entries are 6-coordinate except Sn²⁺, which Shannon lists only for
/// coordination number 8 — check [`IonicRadius::coordination`]. The Fe²⁺,
/// Fe³⁺, Mn²⁺ and Co²⁺ values are the high-spin radii.
pub static IONIC_RADII: &[IonicRadius] = &[
    IonicRadius {
        symbol: "Li",
        charge: 1,
        coordination: 6,
        radius_pm: 76.0,
    },
    IonicRadius {
        symbol: "Na",
        charge: 1,
        coordination: 6,
        radius_pm: 102.0,
    },
    IonicRadius {
        symbol: "K",
        charge: 1,
        coordination: 6,
        radius_pm: 138.0,
    },
    IonicRadius {
        symbol: "Rb",
        charge: 1,
        coordination: 6,
        radius_pm: 152.0,
    },
    IonicRadius {
        symbol: "Cs",
        charge: 1,
        coordination: 6,
        radius_pm: 167.0,
    },
    IonicRadius {
        symbol: "Be",
        charge: 2,
        coordination: 6,
        radius_pm: 45.0,
    },
    IonicRadius {
        symbol: "Mg",
        charge: 2,
        coordination: 6,
        radius_pm: 72.0,
    },
    IonicRadius {
        symbol: "Ca",
        charge: 2,
        coordination: 6,
        radius_pm: 100.0,
    },
    IonicRadius {
        symbol: "Sr",
        charge: 2,
        coordination: 6,
        radius_pm: 118.0,
    },
    IonicRadius {
        symbol: "Ba",
        charge: 2,
        coordination: 6,
        radius_pm: 135.0,
    },
    IonicRadius {
        symbol: "Al",
        charge: 3,
        coordination: 6,
        radius_pm: 53.5,
    },
    IonicRadius {
        symbol: "Fe",
        charge: 2,
        coordination: 6,
        radius_pm: 78.0,
    },
    IonicRadius {
        symbol: "Fe",
        charge: 3,
        coordination: 6,
        radius_pm: 64.5,
    },
    IonicRadius {
        symbol: "Cu",
        charge: 1,
        coordination: 6,
        radius_pm: 77.0,
    },
    IonicRadius {
        symbol: "Cu",
        charge: 2,
        coordination: 6,
        radius_pm: 73.0,
    },
    IonicRadius {
        symbol: "Zn",
        charge: 2,
        coordination: 6,
        radius_pm: 74.0,
    },
    IonicRadius {
        symbol: "Ag",
        charge: 1,
        coordination: 6,
        radius_pm: 115.0,
    },
    IonicRadius {
        symbol: "F",
        charge: -1,
        coordination: 6,
        radius_pm: 133.0,
    },
    IonicRadius {
        symbol: "Cl",
        charge: -1,
        coordination: 6,
        radius_pm: 181.0,
    },
    IonicRadius {
        symbol: "Br",
        charge: -1,
        coordination: 6,
        radius_pm: 196.0,
    },
    IonicRadius {
        symbol: "I",
        charge: -1,
        coordination: 6,
        radius_pm: 220.0,
    },
    IonicRadius {
        symbol: "O",
        charge: -2,
        coordination: 6,
        radius_pm: 140.0,
    },
    IonicRadius {
        symbol: "S",
        charge: -2,
        coordination: 6,
        radius_pm: 184.0,
    },
    IonicRadius {
        symbol: "Mn",
        charge: 2,
        coordination: 6,
        radius_pm: 83.0,
    },
    IonicRadius {
        symbol: "Ni",
        charge: 2,
        coordination: 6,
        radius_pm: 69.0,
    },
    IonicRadius {
        symbol: "Co",
        charge: 2,
        coordination: 6,
        radius_pm: 74.5,
    },
    IonicRadius {
        symbol: "Cr",
        charge: 3,
        coordination: 6,
        radius_pm: 61.5,
    },
    IonicRadius {
        symbol: "Ti",
        charge: 4,
        coordination: 6,
        radius_pm: 60.5,
    },
    IonicRadius {
        symbol: "Pb",
        charge: 2,
        coordination: 6,
        radius_pm: 119.0,
    },
    IonicRadius {
        symbol: "Sn",
        charge: 2,
        // Shannon (1976) lists Sn²⁺ only at CN 8.
        coordination: 8,
        radius_pm: 122.0,
    },
];

/// Look up an ionic radius by symbol and charge (case-sensitive symbol).
///
/// Returns the single [`IONIC_RADII`] entry for that ion, whatever its
/// coordination number.
#[must_use]
#[inline]
pub fn lookup_ionic_radius(symbol: &str, charge: i8) -> Option<&'static IonicRadius> {
    IONIC_RADII
        .iter()
        .find(|r| r.symbol == symbol && r.charge == charge)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Crystal field theory ─────────────────────────────────────────

    #[test]
    fn cfse_d0() {
        assert!((cfse_octahedral(0, true).unwrap()).abs() < f64::EPSILON);
    }

    #[test]
    fn cfse_d3_high_spin() {
        // d³: 3 in t2g → -1.2 Δ_oct
        let cfse = cfse_octahedral(3, true).unwrap();
        assert!((cfse - (-1.2)).abs() < 1e-10);
    }

    #[test]
    fn cfse_d6_high_vs_low() {
        let hs = cfse_octahedral(6, true).unwrap();
        let ls = cfse_octahedral(6, false).unwrap();
        // Low-spin d⁶ (all in t2g) → -2.4, high-spin d⁶ → -0.4
        assert!(ls < hs, "low-spin should have more negative CFSE");
    }

    #[test]
    fn cfse_d10() {
        // d¹⁰: full → CFSE = 0
        let cfse = cfse_octahedral(10, true).unwrap();
        assert!(cfse.abs() < 1e-10);
    }

    #[test]
    fn cfse_invalid_d_electrons() {
        assert!(cfse_octahedral(11, true).is_err());
        assert!(cfse_octahedral(11, false).is_err());
        assert!(cfse_octahedral(u8::MAX, false).is_err());
    }

    #[test]
    fn cfse_all_d_counts_both_spins() {
        // CFSE in units of Δ_oct, pairing energy excluded.
        let high = [
            0.0, -0.4, -0.8, -1.2, -0.6, 0.0, -0.4, -0.8, -1.2, -0.6, 0.0,
        ];
        let low = [
            0.0, -0.4, -0.8, -1.2, -1.6, -2.0, -2.4, -1.8, -1.2, -0.6, 0.0,
        ];
        for n in 0..=10u8 {
            let i = usize::from(n);
            let hs = cfse_octahedral(n, true).unwrap();
            let ls = cfse_octahedral(n, false).unwrap();
            assert!((hs - high[i]).abs() < 1e-12, "HS d{n}: {hs}");
            assert!((ls - low[i]).abs() < 1e-12, "LS d{n}: {ls}");
        }
    }

    #[test]
    fn unpaired_all_d_counts_both_spins() {
        let high = [0u8, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0];
        let low = [0u8, 1, 2, 3, 2, 1, 0, 1, 2, 1, 0];
        for n in 0..=10u8 {
            let i = usize::from(n);
            assert_eq!(unpaired_electrons_oct(n, true).unwrap(), high[i], "HS d{n}");
            assert_eq!(unpaired_electrons_oct(n, false).unwrap(), low[i], "LS d{n}");
        }
    }

    #[test]
    fn unpaired_low_spin_d9_is_one() {
        // Cu²⁺ (d⁹) has one unpaired electron in any octahedral field.
        assert_eq!(unpaired_electrons_oct(9, false).unwrap(), 1);
    }

    #[test]
    fn unpaired_invalid_d_electrons() {
        assert!(unpaired_electrons_oct(11, true).is_err());
        assert!(unpaired_electrons_oct(11, false).is_err());
    }

    #[test]
    fn delta_tet_ratio() {
        let tet = delta_tet_from_oct(10000.0);
        assert!((tet - 4444.4).abs() < 1.0);
    }

    // ── Spectrochemical series ───────────────────────────────────────

    #[test]
    fn ligand_ordering() {
        assert!(Ligand::Cyanide.field_strength() > Ligand::Water.field_strength());
        assert!(Ligand::Water.field_strength() > Ligand::Chloride.field_strength());
    }

    #[test]
    fn strong_field_ligands() {
        assert!(Ligand::Cyanide.is_strong_field());
        assert!(Ligand::CarbonMonoxide.is_strong_field());
        assert!(!Ligand::Water.is_strong_field());
        assert!(!Ligand::Chloride.is_strong_field());
    }

    // ── Magnetic moment ──────────────────────────────────────────────

    #[test]
    fn spin_only_moment_fe3_high_spin() {
        // Fe³⁺ high-spin d⁵: 5 unpaired → μ = √35 ≈ 5.92 BM
        let n = unpaired_electrons_oct(5, true).unwrap();
        assert_eq!(n, 5);
        let mu = spin_only_moment(n);
        assert!((mu - 5.92).abs() < 0.01);
    }

    #[test]
    fn spin_only_moment_zero_unpaired() {
        assert!((spin_only_moment(0)).abs() < f64::EPSILON);
    }

    // ── Lattice energy ───────────────────────────────────────────────

    #[test]
    fn nacl_lattice_energy() {
        // NaCl: M=1.747565, z+=1, z-=1, r0=2.82e-10 m, n=8 (mean of Ne 7, Ar 9)
        // Analytic Born–Landé: −k_e·N_A·M·e²/r₀·(1 − 1/8) = −753.365 kJ/mol
        // (the Born–Haber value, −786 kJ/mol, differs by model error).
        let u = born_lande_lattice_energy(MADELUNG_NACL, 1, 1, 2.82e-10, 8.0).unwrap();
        let u_kj = u / 1000.0;
        assert!(
            (u_kj - (-753.365)).abs() < 0.5,
            "NaCl Born–Landé lattice energy should be -753.4 kJ/mol, got {u_kj}"
        );
    }

    #[test]
    fn lattice_energy_scales_with_charges() {
        let u1 = born_lande_lattice_energy(MADELUNG_NACL, 1, 1, 2.82e-10, 8.0).unwrap();
        let u4 = born_lande_lattice_energy(MADELUNG_NACL, 2, 2, 2.82e-10, 8.0).unwrap();
        assert!((u4 / u1 - 4.0).abs() < 1e-12);
    }

    #[test]
    fn lattice_energy_zero_distance_is_error() {
        assert!(born_lande_lattice_energy(1.748, 1, 1, 0.0, 8.0).is_err());
    }

    #[test]
    fn lattice_energy_born_exp_too_low_is_error() {
        assert!(born_lande_lattice_energy(1.748, 1, 1, 2.82e-10, 1.0).is_err());
    }

    #[test]
    fn lattice_energy_non_finite_inputs_are_errors() {
        let ok = born_lande_lattice_energy(1.748, 1, 1, 2.82e-10, 8.0);
        assert!(ok.is_ok());
        for r0 in [f64::NAN, f64::INFINITY, -2.82e-10] {
            assert!(
                born_lande_lattice_energy(1.748, 1, 1, r0, 8.0).is_err(),
                "r0={r0}"
            );
        }
        for n in [f64::NAN, f64::INFINITY, 0.5] {
            assert!(
                born_lande_lattice_energy(1.748, 1, 1, 2.82e-10, n).is_err(),
                "n={n}"
            );
        }
        for m in [f64::NAN, f64::INFINITY, 0.0, -1.748] {
            assert!(
                born_lande_lattice_energy(m, 1, 1, 2.82e-10, 8.0).is_err(),
                "M={m}"
            );
        }
    }

    #[test]
    fn lattice_energy_zero_charge_is_error() {
        assert!(born_lande_lattice_energy(1.748, 0, 1, 2.82e-10, 8.0).is_err());
        assert!(born_lande_lattice_energy(1.748, 1, 0, 2.82e-10, 8.0).is_err());
    }

    // ── Ionic radii ──────────────────────────────────────────────────

    #[test]
    fn lookup_na_plus() {
        let r = lookup_ionic_radius("Na", 1).unwrap();
        assert!((r.radius_pm - 102.0).abs() < 0.1);
    }

    #[test]
    fn lookup_cl_minus() {
        let r = lookup_ionic_radius("Cl", -1).unwrap();
        assert!((r.radius_pm - 181.0).abs() < 0.1);
    }

    #[test]
    fn lookup_fe_two_oxidation_states() {
        let fe2 = lookup_ionic_radius("Fe", 2).unwrap();
        let fe3 = lookup_ionic_radius("Fe", 3).unwrap();
        assert!(
            fe2.radius_pm > fe3.radius_pm,
            "Fe²⁺ should be larger than Fe³⁺"
        );
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_ionic_radius("Xx", 1).is_none());
    }

    #[test]
    fn cation_smaller_than_anion() {
        let na = lookup_ionic_radius("Na", 1).unwrap();
        let cl = lookup_ionic_radius("Cl", -1).unwrap();
        assert!(na.radius_pm < cl.radius_pm, "cation should be smaller");
    }

    #[test]
    fn lookup_sn2_is_shannon_cn8() {
        let sn = lookup_ionic_radius("Sn", 2).unwrap();
        assert_eq!(sn.coordination, 8);
        assert!((sn.radius_pm - 122.0).abs() < 1e-9);
    }

    #[test]
    fn radii_coordination_matches_doc() {
        // Doc: every entry is 6-coordinate except Sn²⁺ (CN 8).
        for r in IONIC_RADII {
            let expected = if r.symbol == "Sn" && r.charge == 2 {
                8
            } else {
                6
            };
            assert_eq!(r.coordination, expected, "{}{:+}", r.symbol, r.charge);
        }
    }

    #[test]
    fn radii_count() {
        assert!(IONIC_RADII.len() >= 30);
    }

    // ── Born-Haber cycle ─────────────────────────────────────────────

    #[test]
    fn born_haber_nacl() {
        // NaCl Born-Haber cycle:
        // ΔH_f° = -411.15, ΔH_sub(Na) = 107.3, IE(Na) = 495.8
        // ΔH_diss(Cl₂) = 242.0, EA(Cl) = -349.0
        let u = born_haber_lattice_energy(-411.15, 107.3, 495.8, 242.0, -349.0);
        // U = -411.15 - 107.3 - 495.8 - 121.0 - (-349.0) = -786.25
        assert!(
            (u - (-786.25)).abs() < 1.0,
            "NaCl Born-Haber U should be ~-786 kJ/mol, got {u}"
        );
    }

    #[test]
    fn born_haber_matches_born_lande() {
        // Born-Haber (−786) and Born-Landé (−753) agree only to within the
        // ~4 % model error of Born-Landé, hence the loose tolerance.
        let u_bh = born_haber_lattice_energy(-411.15, 107.3, 495.8, 242.0, -349.0);
        let u_bl = born_lande_lattice_energy(MADELUNG_NACL, 1, 1, 2.82e-10, 8.0).unwrap() / 1000.0;
        assert!(
            (u_bh - u_bl).abs() < 50.0,
            "Born-Haber ({u_bh:.0}) and Born-Landé ({u_bl:.0}) should roughly agree"
        );
    }
}
