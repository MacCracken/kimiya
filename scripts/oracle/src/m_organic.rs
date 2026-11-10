//! `organic.tsv` — bond energies, VSEPR geometry and functional groups (port plan §4, M4b).
//!
//! - Class A, every enum value × every property: `BondType::energy_kj` (all 45 variants),
//!   `Geometry::ideal_bond_angle` (all 12), `FunctionalGroup::pattern`, `is_polar` and `hydrogen_bonding` (all 16;
//!   the port's `functional_group_is_hydrogen_bond_donor`, ADR 0002 #6, compared under the Rust name). The input is
//!   `e:<Enum>.<Variant>`; `pattern`'s output an `s:` token, the flags `i:0/1`.
//! - Class A, `organic::predict_geometry` (two Rust u8s, `i:` tokens): the whole square 0 ..= 15 × 0 ..= 15 and every
//!   u8 of each argument against the other's edges (0 … 8, 15, 16, 127, 128, 254, 255), each pair once, in order. Output `e:Geometry.<Variant>`
//!   (`Geometry.TShapedThreeLonePairs` is the port's `KM_GEOMETRY_T_SHAPED`, ADR 0002 #6) or `n:`.
//! - Class B, `organic::enthalpy_from_bonds` (Σ energy·count over each list from −0.0, then broken − formed): both lists
//!   as `v:n` + n pairs `e:BondType.<Variant> i:count` (Rust's `&[(BondType, u32)]`). Empty lists (−0.0 − −0.0 = +0.0),
//!   the Rust tests' lists (2H₂ + O₂, and the integration test's CH₄ + 2O₂), every bond alone on either side with the
//!   counts 0, 1, 2, 3, 7, 1000, 65 535, 2³¹ and u32::MAX, a list against itself (exact cancellation), all 45 bonds at
//!   once, and seeded lists of 1–8 entries with small and large counts. The terms are integers times integers below
//!   2⁵³ or exact products, so no row is NaN: a NaN needs an input Rust's types can't hold (the port's bad lists,
//!   invalid enum values and out-of-u32 counts are Cyrius-only tests).
//!
//! Every input goes through `black_box`.

use std::collections::BTreeSet;
use std::hint::black_box as bb;

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, b, e, f, i, s};
use kimiya::organic::{BondType, FunctionalGroup, Geometry, enthalpy_from_bonds, predict_geometry};

/// Every `BondType`, in Rust declaration order. `t_organic` proves the list against the enum in
/// `rust-old/src/organic.rs` (the enum is `#[non_exhaustive]`, so no match here could catch a missing variant).
pub const BOND_TYPES: [BondType; 45] = [
    BondType::CH,
    BondType::CC,
    BondType::CCDouble,
    BondType::CCTriple,
    BondType::CO,
    BondType::CODouble,
    BondType::CN,
    BondType::CNDouble,
    BondType::CNTriple,
    BondType::CCl,
    BondType::CBr,
    BondType::CF,
    BondType::CS,
    BondType::HH,
    BondType::OH,
    BondType::NH,
    BondType::SH,
    BondType::HF,
    BondType::HCl,
    BondType::HBr,
    BondType::HI,
    BondType::OO,
    BondType::OODouble,
    BondType::ON,
    BondType::NN,
    BondType::NNDouble,
    BondType::NNTriple,
    BondType::PH,
    BondType::PO,
    BondType::PODouble,
    BondType::PCl,
    BondType::SiH,
    BondType::SiO,
    BondType::SiC,
    BondType::SiCl,
    BondType::BH,
    BondType::BO,
    BondType::BN,
    BondType::BF,
    BondType::SS,
    BondType::SODouble,
    BondType::FF,
    BondType::ClCl,
    BondType::BrBr,
    BondType::II,
];

/// Every `Geometry`, in Rust declaration order (proved by `t_organic`).
pub const GEOMETRIES: [Geometry; 12] = [
    Geometry::Linear,
    Geometry::BentTwoLonePairs,
    Geometry::BentOneLonePair,
    Geometry::TrigonalPlanar,
    Geometry::TrigonalPyramidal,
    Geometry::Tetrahedral,
    Geometry::SeesawShaped,
    Geometry::TShapedThreeLonePairs,
    Geometry::TrigonalBipyramidal,
    Geometry::SquarePlanar,
    Geometry::SquarePyramidal,
    Geometry::Octahedral,
];

/// Every `FunctionalGroup`, in Rust declaration order (proved by `t_organic`).
pub const FUNCTIONAL_GROUPS: [FunctionalGroup; 16] = [
    FunctionalGroup::Alkane,
    FunctionalGroup::Alkene,
    FunctionalGroup::Alkyne,
    FunctionalGroup::Alcohol,
    FunctionalGroup::Aldehyde,
    FunctionalGroup::Ketone,
    FunctionalGroup::CarboxylicAcid,
    FunctionalGroup::Ester,
    FunctionalGroup::Ether,
    FunctionalGroup::Amine,
    FunctionalGroup::Amide,
    FunctionalGroup::Nitrile,
    FunctionalGroup::Thiol,
    FunctionalGroup::Halide,
    FunctionalGroup::Phenol,
    FunctionalGroup::AromaticRing,
];

/// A bond list as tokens: `v:n`, then each pair `e:BondType.<Variant> i:count`.
fn bond_toks(list: &[(BondType, u32)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(list.len())];
    for &(bt, n) in list {
        v.push(e("BondType", bt));
        v.push(i(n));
    }
    v
}

fn enthalpy(m: &mut Module, broken: &[(BondType, u32)], formed: &[(BondType, u32)]) {
    let dh = enthalpy_from_bonds(bb(broken), bb(formed));
    let mut ins = bond_toks(broken);
    ins.extend(bond_toks(formed));
    m.row("organic::enthalpy_from_bonds", Class::B, &ins, &[f(dh)]);
}

/// The counts every bond is probed with alone: small, 2³¹ and u32::MAX.
const COUNTS: [u32; 9] = [0, 1, 2, 3, 7, 1000, 65_535, 1 << 31, u32::MAX];

/// The (bonding groups, lone pairs) edges each argument's full u8 sweep is crossed with.
const GEOMETRY_EDGES: [u8; 15] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 127, 128, 254, 255];

fn geometry(m: &mut Module, bonding: u8, lone: u8) {
    let out = predict_geometry(bb(bonding), bb(lone)).map_or(Tok::None, |g| e("Geometry", g));
    m.row("organic::predict_geometry", Class::A, &[i(bonding), i(lone)], &[out]);
}

pub fn emit() -> Module {
    let mut m = Module::new("organic");

    // ── the enums' properties ──
    for bt in BOND_TYPES {
        m.row("BondType::energy_kj", Class::A, &[e("BondType", bt)], &[f(bb(bt).energy_kj())]);
    }
    for g in GEOMETRIES {
        m.row("Geometry::ideal_bond_angle", Class::A, &[e("Geometry", g)], &[f(bb(g).ideal_bond_angle())]);
    }
    for g in FUNCTIONAL_GROUPS {
        m.row("FunctionalGroup::pattern", Class::A, &[e("FunctionalGroup", g)], &[s(bb(g).pattern())]);
    }
    for g in FUNCTIONAL_GROUPS {
        m.row("FunctionalGroup::is_polar", Class::A, &[e("FunctionalGroup", g)], &[b(bb(g).is_polar())]);
    }
    for g in FUNCTIONAL_GROUPS {
        let donor = bb(g).hydrogen_bonding();
        m.row("FunctionalGroup::hydrogen_bonding", Class::A, &[e("FunctionalGroup", g)], &[b(donor)]);
    }

    // ── VSEPR ──
    let mut pairs = BTreeSet::new();
    for bonding in 0..=15u8 {
        for lone in 0..=15u8 {
            pairs.insert((bonding, lone));
        }
    }
    for edge in GEOMETRY_EDGES {
        for x in 0..=u8::MAX {
            pairs.insert((x, edge));
            pairs.insert((edge, x));
        }
    }
    for (bonding, lone) in pairs {
        geometry(&mut m, bonding, lone);
    }

    // ── bond enthalpies ──
    use BondType as T;
    enthalpy(&mut m, &[], &[]);
    enthalpy(&mut m, &[(T::HH, 2), (T::OODouble, 1)], &[(T::OH, 4)]); // 2H₂ + O₂ → 2H₂O, −482 kJ
    enthalpy(&mut m, &[(T::CH, 4), (T::OODouble, 2)], &[(T::CODouble, 2), (T::OH, 4)]); // CH₄ + 2O₂
    enthalpy(&mut m, &[(T::OH, 4)], &[(T::HH, 2), (T::OODouble, 1)]);
    for bt in BOND_TYPES {
        for n in COUNTS {
            enthalpy(&mut m, &[(bt, n)], &[]);
            enthalpy(&mut m, &[], &[(bt, n)]);
        }
        enthalpy(&mut m, &[(bt, 1)], &[(bt, 1)]);
        enthalpy(&mut m, &[(bt, u32::MAX)], &[(bt, u32::MAX)]);
    }
    let all: Vec<(BondType, u32)> = BOND_TYPES.iter().map(|&bt| (bt, 1)).collect();
    let all_max: Vec<(BondType, u32)> = BOND_TYPES.iter().map(|&bt| (bt, u32::MAX)).collect();
    enthalpy(&mut m, &all, &[]);
    enthalpy(&mut m, &all, &all);
    enthalpy(&mut m, &all_max, &all);
    enthalpy(&mut m, &[], &all_max);
    let mut rng = Rng::new(0x4b01);
    let pick = |rng: &mut Rng, large: bool| -> Vec<(BondType, u32)> {
        let len = usize::try_from(rng.int(1, 8)).unwrap_or(1);
        (0..len)
            .map(|_| {
                let k = usize::try_from(rng.int(0, 44)).unwrap_or(0);
                let n = if large { rng.int(0, i64::from(u32::MAX)) } else { rng.int(0, 12) };
                (BOND_TYPES[k], u32::try_from(n).unwrap_or(0))
            })
            .collect()
    };
    for k in 0..600 {
        let large = k % 3 == 2;
        let broken = pick(&mut rng, large);
        let formed = if k % 5 == 4 { Vec::new() } else { pick(&mut rng, large) };
        enthalpy(&mut m, &broken, &formed);
    }
    m
}
