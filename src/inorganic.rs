//! Inorganic chemistry — crystal field theory, lattice energy, ionic radii.

use crate::error::{KimiyaError, Result};
use serde::Serialize;

// ── Crystal field theory ─────────────────────────────────────────────

/// Crystal field splitting energy Δ_oct for octahedral complexes.
///
/// Returns CFSE (Crystal Field Stabilization Energy) in units of Δ_oct.
///
/// - `d_electrons`: number of d electrons (0–10)
/// - `high_spin`: true for high-spin (weak field), false for low-spin (strong field)
///
/// # Errors
///
/// Returns error if d_electrons > 10.
pub fn cfse_octahedral(d_electrons: u8, high_spin: bool) -> Result<f64> {
    if d_electrons > 10 {
        return Err(KimiyaError::InvalidInput(
            "d electron count must be 0-10".into(),
        ));
    }
    // CFSE = (# in t2g × -0.4 + # in eg × +0.6) × Δ_oct
    // High-spin filling: t2g fills 1 each, then eg 1 each, then pair
    // Low-spin filling: t2g fills completely (6) before eg
    let (t2g, eg) = if high_spin {
        match d_electrons {
            0 => (0, 0),
            1 => (1, 0),
            2 => (2, 0),
            3 => (3, 0),
            4 => (3, 1),
            5 => (3, 2),
            6 => (4, 2),
            7 => (5, 2),
            8 => (6, 2),
            9 => (6, 3),
            10 => (6, 4),
            _ => unreachable!(),
        }
    } else {
        match d_electrons {
            0 => (0, 0),
            1 => (1, 0),
            2 => (2, 0),
            3 => (3, 0),
            4 => (4, 0),
            5 => (5, 0),
            6 => (6, 0),
            7 => (6, 1),
            8 => (6, 2),
            9 => (6, 3),
            10 => (6, 4),
            _ => unreachable!(),
        }
    };
    Ok(t2g as f64 * (-0.4) + eg as f64 * 0.6)
}

/// Crystal field splitting energy Δ_tet for tetrahedral complexes.
///
/// Δ_tet ≈ (4/9) × Δ_oct
#[must_use]
#[inline]
pub fn delta_tet_from_oct(delta_oct: f64) -> f64 {
    delta_oct * 4.0 / 9.0
}

/// Spectrochemical series — common ligands ordered by field strength.
///
/// Returns approximate Δ_oct in cm⁻¹ for a given ligand with a generic M²⁺ ion.
/// These are qualitative/relative values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[non_exhaustive]
pub enum Ligand {
    Iodide,
    Bromide,
    Chloride,
    Fluoride,
    Hydroxide,
    Water,
    Ammonia,
    Ethylenediamine,
    Nitrite,
    Cyanide,
    CarbonMonoxide,
}

impl Ligand {
    /// Relative field strength (higher = stronger field).
    #[must_use]
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

    /// Whether this is a strong-field ligand (typically causes low-spin).
    #[must_use]
    pub const fn is_strong_field(self) -> bool {
        self.field_strength() >= 8
    }
}

// ── Spin-only magnetic moment ────────────────────────────────────────

/// Spin-only magnetic moment: μ = √(n(n+2)) Bohr magnetons
///
/// - `unpaired_electrons`: n
#[must_use]
#[inline]
pub fn spin_only_moment(unpaired_electrons: u8) -> f64 {
    let n = unpaired_electrons as f64;
    (n * (n + 2.0)).sqrt()
}

/// Number of unpaired electrons for octahedral d^n.
///
/// # Errors
///
/// Returns error if d_electrons > 10.
pub fn unpaired_electrons_oct(d_electrons: u8, high_spin: bool) -> Result<u8> {
    if d_electrons > 10 {
        return Err(KimiyaError::InvalidInput(
            "d electron count must be 0-10".into(),
        ));
    }
    let n = if high_spin {
        match d_electrons {
            0 | 10 => 0,
            1 | 9 => 1,
            2 | 8 => 2,
            3 | 7 => 3,
            4 | 6 => 4,
            5 => 5,
            _ => unreachable!(),
        }
    } else {
        match d_electrons {
            0 | 6 | 10 => 0,
            1 | 5 | 7 => 1,
            2 | 4 | 8 => 2,
            3 | 9 => 3,
            _ => unreachable!(),
        }
    };
    Ok(n)
}

// ── Lattice energy ───────────────────────────────────────────────────

/// Born-Landé equation for lattice energy:
/// U = -(N_A × M × e² × z⁺ × z⁻) / (4πε₀ × r₀) × (1 - 1/n)
///
/// Simplified form using Coulomb constant:
/// U = -(k_e × N_A × M × z⁺ × z⁻ × e²) / r₀ × (1 - 1/n)
///
/// - `madelung`: Madelung constant M
/// - `z_cation`: cation charge (positive integer)
/// - `z_anion`: anion charge (positive integer, sign handled internally)
/// - `r0_m`: nearest-neighbor distance in meters
/// - `born_exponent`: Born exponent n (typically 5-12)
///
/// Returns lattice energy in J/mol (negative = stable).
///
/// # Errors
///
/// Returns error if parameters are invalid.
pub fn born_lande_lattice_energy(
    madelung: f64,
    z_cation: u32,
    z_anion: u32,
    r0_m: f64,
    born_exponent: f64,
) -> Result<f64> {
    if r0_m <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "interionic distance must be positive".into(),
        ));
    }
    if born_exponent <= 1.0 {
        return Err(KimiyaError::InvalidInput(
            "Born exponent must be > 1".into(),
        ));
    }
    let na = crate::element::AVOGADRO;
    let e = crate::electrochemistry::ELEMENTARY_CHARGE;
    let ke = crate::potential::COULOMB_CONSTANT;

    Ok(
        -ke * na * madelung * z_cation as f64 * z_anion as f64 * e * e / r0_m
            * (1.0 - 1.0 / born_exponent),
    )
}

/// Common Madelung constants.
pub const MADELUNG_NACL: f64 = 1.747565;
pub const MADELUNG_CSCL: f64 = 1.762675;
pub const MADELUNG_ZNS: f64 = 1.6381;
pub const MADELUNG_FLUORITE: f64 = 2.51939;
pub const MADELUNG_RUTILE: f64 = 2.408;

/// Common Born exponents (based on electron configuration).
pub const BORN_EXP_HE: f64 = 5.0;
pub const BORN_EXP_NE: f64 = 7.0;
pub const BORN_EXP_AR: f64 = 9.0;
pub const BORN_EXP_KR: f64 = 10.0;
pub const BORN_EXP_XE: f64 = 12.0;

// ── Ionic radii ──────────────────────────────────────────────────────

/// Shannon ionic radius entry (pm).
#[derive(Debug, Clone, Serialize)]
pub struct IonicRadius {
    pub symbol: &'static str,
    pub charge: i8,
    pub coordination: u8,
    pub radius_pm: f64,
}

/// Common Shannon ionic radii (6-coordinate, pm).
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
        coordination: 6,
        radius_pm: 93.0,
    },
];

/// Look up ionic radius by symbol and charge.
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
        // NaCl: M=1.7476, z+=1, z-=1, r0=2.82e-10 m, n≈8 (avg Ne+Ar)
        let u = born_lande_lattice_energy(MADELUNG_NACL, 1, 1, 2.82e-10, 8.0).unwrap();
        // Known: ~-786 kJ/mol
        let u_kj = u / 1000.0;
        assert!(
            (u_kj - (-786.0)).abs() < 50.0,
            "NaCl lattice energy should be ~-786 kJ/mol, got {u_kj}"
        );
    }

    #[test]
    fn lattice_energy_zero_distance_is_error() {
        assert!(born_lande_lattice_energy(1.748, 1, 1, 0.0, 8.0).is_err());
    }

    #[test]
    fn lattice_energy_born_exp_too_low_is_error() {
        assert!(born_lande_lattice_energy(1.748, 1, 1, 2.82e-10, 1.0).is_err());
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
    fn radii_count() {
        assert!(IONIC_RADII.len() >= 30);
    }
}
