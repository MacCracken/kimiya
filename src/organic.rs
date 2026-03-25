//! Organic chemistry — bond energies, VSEPR geometry, functional groups.

use serde::Serialize;

// ── Bond energies ────────────────────────────────────────────────────

/// A chemical bond type with its average dissociation energy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[non_exhaustive]
pub enum BondType {
    // ── Carbon bonds ─────────────────────────────────────
    CH,
    CC,
    CCDouble,
    CCTriple,
    CO,
    CODouble,
    CN,
    CNDouble,
    CNTriple,
    CCl,
    CBr,
    CF,
    CS,
    // ── Hydrogen bonds ───────────────────────────────────
    HH,
    OH,
    NH,
    SH,
    HF,
    HCl,
    HBr,
    HI,
    // ── Oxygen bonds ─────────────────────────────────────
    OO,
    OODouble,
    ON,
    // ── Nitrogen bonds ───────────────────────────────────
    NN,
    NNDouble,
    NNTriple,
    // ── Phosphorus bonds ─────────────────────────────────
    PH,
    PO,
    PODouble,
    PCl,
    // ── Silicon bonds ──────────────────────────────────────
    SiH,
    SiO,
    SiC,
    SiCl,
    // ── Boron bonds ────────────────────────────────────────
    BH,
    BO,
    BN,
    BF,
    // ── Sulfur bonds ───────────────────────────────────────
    SS,
    SODouble,
    // ── Halogen bonds ────────────────────────────────────
    FF,
    ClCl,
    BrBr,
    II,
}

impl BondType {
    /// Average bond dissociation energy in kJ/mol.
    #[must_use]
    pub const fn energy_kj(self) -> f64 {
        match self {
            // Carbon bonds
            Self::CH => 413.0,
            Self::CC => 347.0,
            Self::CCDouble => 614.0,
            Self::CCTriple => 839.0,
            Self::CO => 358.0,
            Self::CODouble => 745.0,
            Self::CN => 305.0,
            Self::CNDouble => 615.0,
            Self::CNTriple => 891.0,
            Self::CCl => 339.0,
            Self::CBr => 276.0,
            Self::CF => 485.0,
            Self::CS => 259.0,
            // Hydrogen bonds
            Self::HH => 436.0,
            Self::OH => 463.0,
            Self::NH => 391.0,
            Self::SH => 363.0,
            Self::HF => 567.0,
            Self::HCl => 431.0,
            Self::HBr => 366.0,
            Self::HI => 298.0,
            // Oxygen bonds
            Self::OO => 146.0,
            Self::OODouble => 498.0,
            Self::ON => 201.0,
            // Nitrogen bonds
            Self::NN => 160.0,
            Self::NNDouble => 418.0,
            Self::NNTriple => 945.0,
            // Phosphorus bonds
            Self::PH => 322.0,
            Self::PO => 335.0,
            Self::PODouble => 544.0,
            Self::PCl => 326.0,
            // Silicon bonds
            Self::SiH => 318.0,
            Self::SiO => 452.0,
            Self::SiC => 301.0,
            Self::SiCl => 381.0,
            // Boron bonds
            Self::BH => 389.0,
            Self::BO => 536.0,
            Self::BN => 389.0,
            Self::BF => 613.0,
            // Sulfur bonds
            Self::SS => 266.0,
            Self::SODouble => 522.0,
            // Halogen bonds
            Self::FF => 155.0,
            Self::ClCl => 242.0,
            Self::BrBr => 193.0,
            Self::II => 151.0,
        }
    }
}

/// Estimate reaction enthalpy from bond energies:
/// ΔH ≈ Σ(bonds broken) − Σ(bonds formed)
///
/// Each entry is (BondType, count).
#[must_use]
pub fn enthalpy_from_bonds(
    bonds_broken: &[(BondType, u32)],
    bonds_formed: &[(BondType, u32)],
) -> f64 {
    let broken: f64 = bonds_broken
        .iter()
        .map(|(b, n)| b.energy_kj() * *n as f64)
        .sum();
    let formed: f64 = bonds_formed
        .iter()
        .map(|(b, n)| b.energy_kj() * *n as f64)
        .sum();
    broken - formed
}

// ── VSEPR geometry ───────────────────────────────────────────────────

/// Molecular geometry predicted by VSEPR theory.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum Geometry {
    Linear,
    BentTwoLonePairs,
    BentOneLonePair,
    TrigonalPlanar,
    TrigonalPyramidal,
    Tetrahedral,
    SeesawShaped,
    TShapedThreeLonePairs,
    TrigonalBipyramidal,
    SquarePlanar,
    SquarePyramidal,
    Octahedral,
}

impl Geometry {
    /// Ideal bond angle in degrees (where applicable).
    #[must_use]
    pub const fn ideal_bond_angle(self) -> f64 {
        match self {
            Self::Linear => 180.0,
            Self::BentOneLonePair => 117.0, // ~117° due to lone pair compression
            Self::BentTwoLonePairs => 104.5,
            Self::TrigonalPlanar => 120.0,
            Self::TrigonalPyramidal => 107.0,
            Self::Tetrahedral => 109.5,
            Self::SeesawShaped => 120.0, // equatorial
            Self::TShapedThreeLonePairs => 90.0,
            Self::TrigonalBipyramidal => 120.0, // equatorial
            Self::SquarePlanar => 90.0,
            Self::SquarePyramidal => 90.0,
            Self::Octahedral => 90.0,
        }
    }
}

/// Predict molecular geometry from VSEPR theory.
///
/// - `bonding_groups`: number of bonding electron groups (atoms bonded to central atom)
/// - `lone_pairs`: number of lone pairs on the central atom
///
/// Returns `None` for unsupported electron group arrangements.
#[must_use]
pub fn predict_geometry(bonding_groups: u8, lone_pairs: u8) -> Option<Geometry> {
    match (bonding_groups, lone_pairs) {
        (2, 0) => Some(Geometry::Linear),
        (2, 1) => Some(Geometry::BentOneLonePair),
        (2, 2) => Some(Geometry::BentTwoLonePairs),
        (2, 3) => Some(Geometry::Linear), // XeF2-like
        (3, 0) => Some(Geometry::TrigonalPlanar),
        (3, 1) => Some(Geometry::TrigonalPyramidal),
        (3, 2) => Some(Geometry::TShapedThreeLonePairs),
        (4, 0) => Some(Geometry::Tetrahedral),
        (4, 1) => Some(Geometry::SeesawShaped),
        (4, 2) => Some(Geometry::SquarePlanar),
        (5, 0) => Some(Geometry::TrigonalBipyramidal),
        (5, 1) => Some(Geometry::SquarePyramidal),
        (6, 0) => Some(Geometry::Octahedral),
        _ => None,
    }
}

// ── Functional groups ────────────────────────────────────────────────

/// Common organic functional groups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[non_exhaustive]
pub enum FunctionalGroup {
    Alkane,
    Alkene,
    Alkyne,
    Alcohol,
    Aldehyde,
    Ketone,
    CarboxylicAcid,
    Ester,
    Ether,
    Amine,
    Amide,
    Nitrile,
    Thiol,
    Halide,
    Phenol,
    AromaticRing,
}

impl FunctionalGroup {
    /// General formula or identifying pattern.
    #[must_use]
    pub const fn pattern(self) -> &'static str {
        match self {
            Self::Alkane => "C-C (sp3)",
            Self::Alkene => "C=C (sp2)",
            Self::Alkyne => "C≡C (sp)",
            Self::Alcohol => "-OH",
            Self::Aldehyde => "-CHO",
            Self::Ketone => "C(=O)C",
            Self::CarboxylicAcid => "-COOH",
            Self::Ester => "-COO-",
            Self::Ether => "C-O-C",
            Self::Amine => "-NH₂ / -NHR / -NR₂",
            Self::Amide => "-CONH₂",
            Self::Nitrile => "-C≡N",
            Self::Thiol => "-SH",
            Self::Halide => "-X (F, Cl, Br, I)",
            Self::Phenol => "Ar-OH",
            Self::AromaticRing => "C₆H₅-",
        }
    }

    /// Typical polarity (qualitative).
    #[must_use]
    pub const fn is_polar(self) -> bool {
        !matches!(
            self,
            Self::Alkane | Self::Alkene | Self::Alkyne | Self::AromaticRing
        )
    }

    /// Whether the group can participate in hydrogen bonding.
    #[must_use]
    pub const fn hydrogen_bonding(self) -> bool {
        matches!(
            self,
            Self::Alcohol
                | Self::CarboxylicAcid
                | Self::Amine
                | Self::Amide
                | Self::Phenol
                | Self::Thiol
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Bond energies ────────────────────────────────────────────────

    #[test]
    fn bond_energy_ordering() {
        // Triple > double > single for carbon-carbon
        assert!(BondType::CCTriple.energy_kj() > BondType::CCDouble.energy_kj());
        assert!(BondType::CCDouble.energy_kj() > BondType::CC.energy_kj());
    }

    #[test]
    fn bond_energy_halogens_decrease_down_group() {
        // F-F is anomalously weak, but H-X decreases: HF > HCl > HBr > HI
        assert!(BondType::HF.energy_kj() > BondType::HCl.energy_kj());
        assert!(BondType::HCl.energy_kj() > BondType::HBr.energy_kj());
        assert!(BondType::HBr.energy_kj() > BondType::HI.energy_kj());
    }

    #[test]
    fn bond_energy_all_positive() {
        let bonds = [
            BondType::CH,
            BondType::CC,
            BondType::CCDouble,
            BondType::CCTriple,
            BondType::CO,
            BondType::CODouble,
            BondType::CN,
            BondType::CNDouble,
            BondType::CNTriple,
            BondType::HH,
            BondType::OH,
            BondType::NH,
            BondType::OO,
            BondType::OODouble,
            BondType::NN,
            BondType::NNTriple,
            BondType::FF,
            BondType::ClCl,
            BondType::BrBr,
            BondType::II,
        ];
        for b in &bonds {
            assert!(b.energy_kj() > 0.0, "{b:?} should have positive energy");
        }
    }

    #[test]
    fn enthalpy_from_bonds_h2_combustion() {
        // H₂ + ½O₂ → H₂O
        // Broken: 1×H-H + 0.5×O=O = 436 + 249 = 685
        // Formed: 2×O-H = 926
        // ΔH ≈ 685 - 926 = -241 kJ/mol (actual: -242 kJ/mol)
        let dh = enthalpy_from_bonds(&[(BondType::HH, 1)], &[(BondType::OH, 2)]);
        // We'd also need to account for 0.5 O=O broken, but this tests the API
        // Full calculation:
        let dh_full = 1.0 * BondType::HH.energy_kj() + 0.5 * BondType::OODouble.energy_kj()
            - 2.0 * BondType::OH.energy_kj();
        assert!(
            (dh_full - (-241.0)).abs() < 10.0,
            "H₂ combustion should be ~-241 kJ/mol, got {dh_full}"
        );
        // Basic API test with just H-H broken and O-H formed
        assert!(dh < 0.0, "forming stronger bonds should be exothermic");
    }

    // ── VSEPR ────────────────────────────────────────────────────────

    #[test]
    fn vsepr_water() {
        // H₂O: 2 bonding, 2 lone pairs → bent
        assert_eq!(predict_geometry(2, 2), Some(Geometry::BentTwoLonePairs));
    }

    #[test]
    fn vsepr_methane() {
        // CH₄: 4 bonding, 0 lone pairs → tetrahedral
        assert_eq!(predict_geometry(4, 0), Some(Geometry::Tetrahedral));
    }

    #[test]
    fn vsepr_ammonia() {
        // NH₃: 3 bonding, 1 lone pair → trigonal pyramidal
        assert_eq!(predict_geometry(3, 1), Some(Geometry::TrigonalPyramidal));
    }

    #[test]
    fn vsepr_co2() {
        // CO₂: 2 bonding, 0 lone pairs → linear
        assert_eq!(predict_geometry(2, 0), Some(Geometry::Linear));
    }

    #[test]
    fn vsepr_bf3() {
        // BF₃: 3 bonding, 0 lone pairs → trigonal planar
        assert_eq!(predict_geometry(3, 0), Some(Geometry::TrigonalPlanar));
    }

    #[test]
    fn vsepr_sf6() {
        // SF₆: 6 bonding, 0 lone pairs → octahedral
        assert_eq!(predict_geometry(6, 0), Some(Geometry::Octahedral));
    }

    #[test]
    fn vsepr_xef2() {
        // XeF₂: 2 bonding, 3 lone pairs → linear
        assert_eq!(predict_geometry(2, 3), Some(Geometry::Linear));
    }

    #[test]
    fn vsepr_unsupported() {
        assert!(predict_geometry(7, 0).is_none());
    }

    #[test]
    fn tetrahedral_bond_angle() {
        assert!((Geometry::Tetrahedral.ideal_bond_angle() - 109.5).abs() < f64::EPSILON);
    }

    // ── Functional groups ────────────────────────────────────────────

    #[test]
    fn alkane_nonpolar() {
        assert!(!FunctionalGroup::Alkane.is_polar());
    }

    #[test]
    fn alcohol_polar_and_hbonds() {
        assert!(FunctionalGroup::Alcohol.is_polar());
        assert!(FunctionalGroup::Alcohol.hydrogen_bonding());
    }

    #[test]
    fn ether_polar_no_hbonds() {
        assert!(FunctionalGroup::Ether.is_polar());
        assert!(!FunctionalGroup::Ether.hydrogen_bonding());
    }

    #[test]
    fn carboxylic_acid_hbonds() {
        assert!(FunctionalGroup::CarboxylicAcid.hydrogen_bonding());
    }
}
