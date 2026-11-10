//! `src/organic_data.cyr` — the organic-chemistry data of `rust-old/src/organic.rs`: the enums `BondType`, `Geometry`
//! and `FunctionalGroup` and their per-variant properties as baked columns.
//!
//! - `enum KmBondType`, `enum KmGeometry`, `enum KmFunctionalGroup` (with `KM_<ENUM>_COUNT`): the Rust variants in
//!   declaration order, value k for the k-th (append-only), `KM_` + snake(enum) + `_` + snake(variant), under the Rust
//!   doc comments and the enum body's `//` section headings (`── Carbon bonds ──`). `m_organic`'s variant lists are
//!   proved against the source declarations (`cyr::enum_decl`): the enums are `#[non_exhaustive]`, so no match here
//!   could catch a missing variant, and `validate.py` parses the Rust enums again and checks the generated ones.
//!   **ADR 0002 #6:** `Geometry::TShapedThreeLonePairs` (a misnomer: two lone pairs) is `KM_GEOMETRY_T_SHAPED`, at its
//!   Rust position.
//! - The `const fn` properties, read through the API for every variant: `BondType::energy_kj` (f64, with the `//`
//!   section comments of its `match` arms, `// Carbon bonds`, inside the initializer), `Geometry::ideal_bond_angle`
//!   (f64, each arm's trailing note as its cell's label), `FunctionalGroup::pattern` (a **UTF-8** blob, "C≡C (sp)"),
//!   `is_polar` and `hydrogen_bonding` (bool columns; the latter is the port's
//!   `functional_group_is_hydrogen_bond_donor`, ADR 0002 #6). The match arms are checked to cover every variant once.
//! - **ADR 0002 #6's new `functional_group_is_hydrogen_bond_acceptor`:** the oracle has no such data, so its column
//!   is [`acceptor`]'s classification, with its sources, below: IUPAC 2011's criterion, an electron-rich region (a
//!   lone pair or a π-bonded pair), the weak acceptors included (the user's decision). It has no golden rows (ADR
//!   0002 #6: reference tests only, `tests/organic.tcyr`); `validate.py` recomputes it from the patterns by the rule.
//!
//! Values come from the oracle's API; the Rust source is read for comments only (and the variant and arm names, as a
//! completeness check).

use crate::cyr::{
    Elem, Member, doc_comment, enum_decl, header, match_arms, snake, string_column_utf8, table,
    table_sectioned,
};
use crate::m_organic::{BOND_TYPES, FUNCTIONAL_GROUPS, GEOMETRIES};
use kimiya::organic::FunctionalGroup;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// The sources of the acceptor classification, as the column's doc lines.
const ACCEPTOR_DOC: &str = "\
ADR 0002 #6 adds this property: whether the group can act as a hydrogen-bond acceptor. The oracle has no
such data, so this column is scripts/oracle/src/t_organic.rs's classification, not the oracle's API.
- The criterion is IUPAC's (Recommendations 2011: E. Arunan et al., \"Definition of the hydrogen bond\",
  Pure Appl. Chem. 83 (2011) 1637–1641, doi:10.1351/PAC-REC-10-01-02). In X–H···Y–Z, X–H is the donor and
  \"the acceptor is an electron-rich region such as, but not limited to, a lone pair in Y or a π-bonded
  pair in Y–Z\" (the wording as reproduced from the Recommendation in S. M. Bachrach's Computational
  Organic Chemistry blog, comporgchem.com/blog/archives/1402: the published PDF was not reachable to
  check it against). The definition also asks for evidence of bond formation; each weak family has it:
- O and N lone pairs, the conventional acceptors (J. McMurry, Organic Chemistry, OpenStax 2023, §2.12;
  the oracle's own doc of `hydrogen_bonding` names ethers, aldehydes, ketones, esters and nitriles as
  acceptors that cannot donate).
- π-bonded pairs, X–H···π to C=C, C≡C and aromatic rings: T. Steiner, \"The hydrogen bond in the solid
  state\", Angew. Chem. Int. Ed. 41 (2002) 48–76, doi:10.1002/1521-3773(20020104)41:1<48::AID-ANIE48>3.0.CO;2-U.
- Divalent sulfur, the thiol's S: a weak acceptor (F. H. Allen, C. M. Bird, R. S. Rowland, P. R. Raithby,
  \"Hydrogen-bond acceptor and donor properties of divalent sulfur (Y–S–Z and R–S–H)\", Acta Cryst. B53
  (1997) 696–701, doi:10.1107/S0108768197002644).
- Carbon-bound halogens: weak acceptors (L. Brammer, E. A. Bruton, P. Sherwood, \"Understanding the
  behavior of halogens as hydrogen bond acceptors\", Cryst. Growth Des. 1 (2001) 277–290,
  doi:10.1021/cg015522k), organic fluorine the weakest (J. D. Dunitz, R. Taylor, \"Organic fluorine hardly
  ever accepts hydrogen bonds\", Chem. Eur. J. 3 (1997) 89–98, doi:10.1002/chem.19970030115).
- An alkane has no lone pair and no π bond, only C–C and C–H σ bonds, so it is the one group that does
  not accept. No strength is implied: a 1 here can be a carbonyl O or the π face of a benzene ring.";

/// ADR 0002 #6's `is_hydrogen_bond_acceptor` for one group, with the reason its cell's label carries. The rule and
/// its sources are [`ACCEPTOR_DOC`]: IUPAC 2011's electron-rich region, a lone pair (O, N, S, halogen) or a
/// π-bonded pair, with a source for each weak family (the user's decision, 2026-11-09: the weak acceptors count).
fn acceptor(g: FunctionalGroup) -> Result<(bool, &'static str), String> {
    use FunctionalGroup as G;
    Ok(match g {
        G::Alkane => (false, "no lone pair and no π bond: C–C and C–H σ bonds only"),
        G::Alkene => (true, "the C=C π-bonded pair, a weak X–H···π acceptor (Steiner 2002)"),
        G::Alkyne => (true, "the C≡C π-bonded pairs, weak X–H···π acceptors (Steiner 2002)"),
        G::Alcohol => (true, "the O–H oxygen's two lone pairs"),
        G::Aldehyde => (true, "the carbonyl O's lone pairs (Rust: acceptor-only)"),
        G::Ketone => (true, "the carbonyl O's lone pairs (Rust: acceptor-only)"),
        G::CarboxylicAcid => (true, "the C=O and O–H oxygens' lone pairs"),
        G::Ester => (true, "the carbonyl and ether-type O's lone pairs (Rust: acceptor-only)"),
        G::Ether => (true, "the C–O–C oxygen's two lone pairs (Rust: acceptor-only)"),
        G::Amine => (true, "the N's lone pair, tertiary amines included"),
        G::Amide => (true, "chiefly the carbonyl O's lone pairs (the N pair is delocalised)"),
        G::Nitrile => (true, "the sp N's lone pair (Rust: acceptor-only)"),
        G::Thiol => (true, "the divalent S's lone pairs, a weak acceptor (Allen et al. 1997)"),
        G::Halide => (true, "the C–X halogen's lone pairs, weak; C–F the weakest (Brammer 2001, Dunitz 1997)"),
        G::Phenol => (true, "the O–H oxygen's lone pairs, and the ring's π face"),
        G::AromaticRing => (true, "the ring's π face, a weak X–H···π acceptor (Steiner 2002)"),
        _ => return Err(format!("{g:?}: no hydrogen-bond acceptor classification in t_organic::acceptor")),
    })
}

/// The members of a generated enum: `KM_<SNAKE(enum)>_<SNAKE(variant)>`, but for the ADR renames in `renames`
/// (variant → member and its note).
fn members(prefix: &str, variants: &[String], renames: &[(&str, &str, &str)]) -> Vec<Member> {
    variants
        .iter()
        .map(|v| match renames.iter().find(|(rv, _, _)| rv == v) {
            Some((_, name, note)) => Member {
                variant: v.clone(),
                name: (*name).to_string(),
                extra: vec![(*note).to_string()],
            },
            None => Member {
                variant: v.clone(),
                name: format!("{prefix}{}", snake(v)),
                extra: Vec::new(),
            },
        })
        .collect()
}

/// A method's Rust doc as the plain lines of a column doc, then which Rust item the column mirrors.
fn method_doc(src: &str, item: &str, mirrors: &str, index: &str) -> Result<String, String> {
    let docs = doc_comment(src, item)?;
    let mut out = String::new();
    for l in docs.lines() {
        let _ = writeln!(out, "{}", l.trim_start_matches('#').trim_start());
    }
    let _ = write!(out, "{index}; mirrors `{mirrors}`.");
    Ok(out)
}

/// The arms of `item`'s `match self`, checked to name every variant once: (sections before each variant's row, in
/// declaration order; each variant's trailing note). Sections are kept only when the arms are in declaration order,
/// so a heading stays above the rows it heads.
type Arms = (Vec<Vec<String>>, BTreeMap<String, String>);
fn arms(src: &str, item: &str, variants: &[String]) -> Result<Arms, String> {
    let arms = match_arms(src, item)?;
    let mut named: Vec<&str> = arms.iter().map(|a| a.variant.as_str()).collect();
    named.sort_unstable();
    let mut want: Vec<&str> = variants.iter().map(String::as_str).collect();
    want.sort_unstable();
    if named != want {
        return Err(format!("{item}: the match arms {named:?} are not the variants {want:?}"));
    }
    let in_order = arms.iter().map(|a| a.variant.as_str()).eq(variants.iter().map(String::as_str));
    let mut sections = vec![Vec::new(); variants.len()];
    let mut notes = BTreeMap::new();
    for a in &arms {
        if !a.sections.is_empty() {
            if !in_order {
                return Err(format!("{item}: section comments in a match out of declaration order"));
            }
            let k = variants.iter().position(|v| *v == a.variant).unwrap_or(0);
            sections[k].clone_from(&a.sections);
        }
        if let Some(n) = &a.note {
            notes.insert(a.variant.clone(), n.clone());
        }
    }
    Ok((sections, notes))
}

#[allow(clippy::too_many_lines)] // one enum after another, each with its columns
pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let bonds: Vec<String> = BOND_TYPES.iter().map(|v| format!("{v:?}")).collect();
    let geoms: Vec<String> = GEOMETRIES.iter().map(|v| format!("{v:?}")).collect();
    let groups: Vec<String> = FUNCTIONAL_GROUPS.iter().map(|v| format!("{v:?}")).collect();

    let mut out = header(
        "organic-chemistry data: bond energies, VSEPR geometries, functional groups",
        provenance,
        &[
            "Each enum follows the Rust declaration order and is append-only: value k is the k-th variant, and a\n\
             property column is indexed by the value. f64 cells are the oracle's values as round-trip decimals;\n\
             tests/golden/organic.tsv pins every cell. ADR 0002 #6: Geometry::TShapedThreeLonePairs is\n\
             KM_GEOMETRY_T_SHAPED, FunctionalGroup::hydrogen_bonding is the donor column, and the acceptor column\n\
             is new (its sources are in its comment; it has no golden rows).",
        ],
    );

    // ── BondType ──
    out.push_str(&enum_decl(
        src,
        "BondType",
        "KmBondType",
        "KM_BOND_TYPE_COUNT",
        &members("KM_BOND_TYPE_", &bonds, &[]),
    )?);
    let (sections, notes) = arms(src, "pub const fn energy_kj(", &bonds)?;
    if !notes.is_empty() {
        return Err("energy_kj: a trailing arm note, which this generator does not carry".to_string());
    }
    let energy: Vec<f64> = BOND_TYPES.iter().map(|b| b.energy_kj()).collect();
    out.push_str(&table_sectioned(
        "KM_BOND_TYPE_ENERGY_KJ",
        &method_doc(src, "pub const fn energy_kj(", "BondType::energy_kj", "Indexed by KmBondType value")?,
        &Elem::F64(energy),
        &bonds,
        &sections,
    )?);

    // ── Geometry ──
    out.push_str(&enum_decl(
        src,
        "Geometry",
        "KmGeometry",
        "KM_GEOMETRY_COUNT",
        &members(
            "KM_GEOMETRY_",
            &geoms,
            &[(
                "TShapedThreeLonePairs",
                "KM_GEOMETRY_T_SHAPED",
                "Rust: `Geometry::TShapedThreeLonePairs`, renamed by ADR 0002 #6 (docs/adr/0002-port-api-divergences.md).",
            )],
        ),
    )?);
    let (sections, notes) = arms(src, "pub const fn ideal_bond_angle(", &geoms)?;
    if sections.iter().any(|s| !s.is_empty()) {
        return Err("ideal_bond_angle: section comments, which this generator does not carry".to_string());
    }
    let labels: Vec<String> = geoms
        .iter()
        .map(|g| notes.get(g).map_or_else(|| g.clone(), |n| format!("{g}: {n}")))
        .collect();
    let angles: Vec<f64> = GEOMETRIES.iter().map(|g| g.ideal_bond_angle()).collect();
    out.push_str(&table(
        "KM_GEOMETRY_IDEAL_BOND_ANGLE",
        &method_doc(
            src,
            "pub const fn ideal_bond_angle(",
            "Geometry::ideal_bond_angle",
            "Indexed by KmGeometry value",
        )?,
        &Elem::F64(angles),
        &labels,
    )?);

    // ── FunctionalGroup ──
    out.push_str(&enum_decl(
        src,
        "FunctionalGroup",
        "KmFunctionalGroup",
        "KM_FUNCTIONAL_GROUP_COUNT",
        &members("KM_FUNCTIONAL_GROUP_", &groups, &[]),
    )?);
    let (sections, notes) = arms(src, "pub const fn pattern(", &groups)?;
    if sections.iter().any(|s| !s.is_empty()) || !notes.is_empty() {
        return Err("pattern: comments in the match, which this generator does not carry".to_string());
    }
    let patterns: Vec<&str> = FUNCTIONAL_GROUPS.iter().map(|g| g.pattern()).collect();
    let idx = "Indexed by KmFunctionalGroup value";
    out.push_str(&string_column_utf8(
        "_KM_FUNCTIONAL_GROUP_PATTERN_BYTES",
        "_KM_FUNCTIONAL_GROUP_PATTERN_OFF",
        &method_doc(src, "pub const fn pattern(", "FunctionalGroup::pattern", idx)?,
        &patterns,
    )?);
    let polar: Vec<bool> = FUNCTIONAL_GROUPS.iter().map(|g| g.is_polar()).collect();
    out.push_str(&table(
        "KM_FUNCTIONAL_GROUP_IS_POLAR",
        &method_doc(src, "pub const fn is_polar(", "FunctionalGroup::is_polar", idx)?,
        &Elem::Bool(polar),
        &groups,
    )?);
    let donor: Vec<bool> = FUNCTIONAL_GROUPS.iter().map(|g| g.hydrogen_bonding()).collect();
    out.push_str(&table(
        "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_DONOR",
        &format!(
            "{}\nADR 0002 #6: Rust's `hydrogen_bonding` is functional_group_is_hydrogen_bond_donor here.",
            method_doc(src, "pub const fn hydrogen_bonding(", "FunctionalGroup::hydrogen_bonding", idx)?
        ),
        &Elem::Bool(donor),
        &groups,
    )?);
    let mut flags = Vec::new();
    let mut reasons = Vec::new();
    for (g, name) in FUNCTIONAL_GROUPS.iter().zip(&groups) {
        let (flag, why) = acceptor(*g)?;
        flags.push(flag);
        reasons.push(format!("{name}: {why}"));
    }
    out.push_str(&table(
        "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_ACCEPTOR",
        &format!("{ACCEPTOR_DOC}\n{idx}; no Rust counterpart."),
        &Elem::Bool(flags),
        &reasons,
    )?);
    Ok(out.trim_end().to_string() + "\n")
}
