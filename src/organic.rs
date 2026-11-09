//! Organic chemistry — bond energies, VSEPR geometry, functional groups.

use serde::Serialize;

// ── Bond energies ────────────────────────────────────────────────────

/// A covalent bond type with a tabulated average bond dissociation energy
/// (see [`BondType::energy_kj`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[non_exhaustive]
pub enum BondType {
    // ── Carbon bonds ─────────────────────────────────────
    /// C–H single bond.
    CH,
    /// C–C single bond.
    CC,
    /// C=C double bond.
    CCDouble,
    /// C≡C triple bond.
    CCTriple,
    /// C–O single bond.
    CO,
    /// C=O double bond (carbonyl; average value).
    CODouble,
    /// C–N single bond.
    CN,
    /// C=N double bond.
    CNDouble,
    /// C≡N triple bond.
    CNTriple,
    /// C–Cl single bond.
    CCl,
    /// C–Br single bond.
    CBr,
    /// C–F single bond.
    CF,
    /// C–S single bond.
    CS,
    // ── Hydrogen bonds ───────────────────────────────────
    /// H–H bond.
    HH,
    /// O–H bond.
    OH,
    /// N–H bond.
    NH,
    /// S–H bond.
    SH,
    /// H–F bond.
    HF,
    /// H–Cl bond.
    HCl,
    /// H–Br bond.
    HBr,
    /// H–I bond.
    HI,
    // ── Oxygen bonds ─────────────────────────────────────
    /// O–O single bond (peroxide).
    OO,
    /// O=O double bond (O₂).
    OODouble,
    /// O–N single bond.
    ON,
    // ── Nitrogen bonds ───────────────────────────────────
    /// N–N single bond.
    NN,
    /// N=N double bond.
    NNDouble,
    /// N≡N triple bond (N₂).
    NNTriple,
    // ── Phosphorus bonds ─────────────────────────────────
    /// P–H bond.
    PH,
    /// P–O single bond.
    PO,
    /// P=O double bond.
    PODouble,
    /// P–Cl bond.
    PCl,
    // ── Silicon bonds ──────────────────────────────────────
    /// Si–H bond.
    SiH,
    /// Si–O single bond.
    SiO,
    /// Si–C bond.
    SiC,
    /// Si–Cl bond.
    SiCl,
    // ── Boron bonds ────────────────────────────────────────
    /// B–H bond.
    BH,
    /// B–O bond.
    BO,
    /// B–N bond.
    BN,
    /// B–F bond.
    BF,
    // ── Sulfur bonds ───────────────────────────────────────
    /// S–S single bond.
    SS,
    /// S=O double bond.
    SODouble,
    // ── Halogen bonds ────────────────────────────────────
    /// F–F bond (F₂).
    FF,
    /// Cl–Cl bond (Cl₂).
    ClCl,
    /// Br–Br bond (Br₂).
    BrBr,
    /// I–I bond (I₂).
    II,
}

impl BondType {
    /// Average bond dissociation energy in kJ/mol (gas phase, averaged over
    /// typical compounds).
    #[must_use]
    #[inline]
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

/// Estimate reaction enthalpy from average bond energies, in kJ (per mole of
/// reaction as written):
/// ΔH ≈ Σ(bonds broken) − Σ(bonds formed)
///
/// Each entry is (BondType, count). Negative = exothermic.
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
    /// Linear (180°): 2 bonding groups with 0 or 3 lone pairs (CO₂, XeF₂), or a diatomic.
    Linear,
    /// Bent with 2 lone pairs (~104.5°): 2 bonding groups + 2 lone pairs (H₂O).
    BentTwoLonePairs,
    /// Bent with 1 lone pair (<120°): 2 bonding groups + 1 lone pair (SO₂).
    BentOneLonePair,
    /// Trigonal planar (120°): 3 bonding groups, no lone pairs (BF₃).
    TrigonalPlanar,
    /// Trigonal pyramidal (~107°): 3 bonding groups + 1 lone pair (NH₃).
    TrigonalPyramidal,
    /// Tetrahedral (109.5°): 4 bonding groups, no lone pairs (CH₄).
    Tetrahedral,
    /// Seesaw: 4 bonding groups + 1 lone pair (SF₄).
    SeesawShaped,
    /// T-shaped (90°): 3 bonding groups + **2** lone pairs (ClF₃). The
    /// variant name is a misnomer — there are two lone pairs, not three.
    TShapedThreeLonePairs,
    /// Trigonal bipyramidal (120° equatorial, 90° axial): 5 bonding groups (PCl₅).
    TrigonalBipyramidal,
    /// Square planar (90°): 4 bonding groups + 2 lone pairs (XeF₄).
    SquarePlanar,
    /// Square pyramidal (~90°): 5 bonding groups + 1 lone pair (BrF₅).
    SquarePyramidal,
    /// Octahedral (90°): 6 bonding groups, no lone pairs (SF₆).
    Octahedral,
}

impl Geometry {
    /// Ideal bond angle in degrees (where applicable; equatorial angle for
    /// the trigonal-bipyramidal family).
    #[must_use]
    #[inline]
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
/// A single bonding group (a diatomic such as HCl, or a terminal atom) is
/// always [`Geometry::Linear`] for 0–3 lone pairs.
///
/// Returns `None` for unsupported electron group arrangements.
#[must_use]
#[inline]
pub fn predict_geometry(bonding_groups: u8, lone_pairs: u8) -> Option<Geometry> {
    match (bonding_groups, lone_pairs) {
        (1, 0..=3) => Some(Geometry::Linear), // diatomic: HCl, N₂, F₂
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
    /// Alkane: saturated C–C/C–H framework (sp³).
    Alkane,
    /// Alkene: C=C double bond (sp²).
    Alkene,
    /// Alkyne: C≡C triple bond (sp).
    Alkyne,
    /// Alcohol: –OH on an sp³ carbon.
    Alcohol,
    /// Aldehyde: terminal carbonyl, –CHO.
    Aldehyde,
    /// Ketone: internal carbonyl, C(=O)C.
    Ketone,
    /// Carboxylic acid: –COOH.
    CarboxylicAcid,
    /// Ester: –COO–R.
    Ester,
    /// Ether: C–O–C.
    Ether,
    /// Amine: primary, secondary or tertiary (–NH₂ / –NHR / –NR₂).
    Amine,
    /// Amide: –CONH₂ (and N-substituted).
    Amide,
    /// Nitrile: –C≡N.
    Nitrile,
    /// Thiol: –SH.
    Thiol,
    /// Alkyl/aryl halide: –X (F, Cl, Br, I).
    Halide,
    /// Phenol: –OH on an aromatic ring.
    Phenol,
    /// Aromatic ring: phenyl, C₆H₅–.
    AromaticRing,
}

impl FunctionalGroup {
    /// General formula or identifying pattern.
    #[must_use]
    #[inline]
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

    /// Whether the group can act as a hydrogen-bond **donor** (has an O–H,
    /// N–H or S–H).
    ///
    /// Acceptor-only groups (ethers, aldehydes, ketones, esters, nitriles)
    /// return `false` even though they accept hydrogen bonds. [`Self::Amine`]
    /// returns `true` for primary/secondary amines; tertiary amines (–NR₂)
    /// have no N–H and cannot donate. S–H is only a weak donor.
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
        // 2H₂ + O₂ → 2H₂O
        // Broken: 2×H–H + 1×O=O = 2×436 + 498 = 1370 kJ
        // Formed: 4×O–H = 4×463 = 1852 kJ
        // ΔH = 1370 − 1852 = −482 kJ (≈ 2 × −241.8 kJ/mol experimental, gas)
        let dh = enthalpy_from_bonds(
            &[(BondType::HH, 2), (BondType::OODouble, 1)],
            &[(BondType::OH, 4)],
        );
        assert!(
            (dh - (-482.0)).abs() < 1e-9,
            "2H₂ + O₂ → 2H₂O should be -482 kJ from the table, got {dh}"
        );
    }

    #[test]
    fn enthalpy_from_bonds_empty_is_zero() {
        assert!(enthalpy_from_bonds(&[], &[]).abs() < f64::EPSILON);
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
    fn vsepr_diatomic_linear() {
        // HCl, N₂, F₂: one bonding group → linear regardless of lone pairs
        for lp in 0..=3 {
            assert_eq!(predict_geometry(1, lp), Some(Geometry::Linear), "lp={lp}");
        }
    }

    #[test]
    fn vsepr_clf3_t_shaped() {
        // ClF₃: 3 bonding groups + 2 lone pairs → T-shaped
        assert_eq!(
            predict_geometry(3, 2),
            Some(Geometry::TShapedThreeLonePairs)
        );
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
    fn ether_polar_not_hbond_donor() {
        // Ethers accept hydrogen bonds but have no O–H to donate.
        assert!(FunctionalGroup::Ether.is_polar());
        assert!(!FunctionalGroup::Ether.hydrogen_bonding());
    }

    #[test]
    fn acceptor_only_groups_are_not_donors() {
        for g in [
            FunctionalGroup::Ether,
            FunctionalGroup::Ketone,
            FunctionalGroup::Aldehyde,
            FunctionalGroup::Ester,
            FunctionalGroup::Nitrile,
        ] {
            assert!(!g.hydrogen_bonding(), "{g:?} is an acceptor, not a donor");
        }
    }

    #[test]
    fn carboxylic_acid_hbonds() {
        assert!(FunctionalGroup::CarboxylicAcid.hydrogen_bonding());
    }
}
