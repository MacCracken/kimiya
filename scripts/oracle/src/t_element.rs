//! `src/element_data.cyr` — the periodic table (`ELEMENTS`, Z = 1–118) and the extended
//! properties of 32 elements (`ELEMENT_PROPERTIES`), as baked column arrays.
//!
//! - `enum KmElementCategory`: the Rust `ElementCategory` variants in declaration order,
//!   `KM_ELEMENT_CATEGORY_` + snake(variant). `e:ElementCategory.NobleGas` in the golden corpus
//!   is `KM_ELEMENT_CATEGORY_NOBLE_GAS`.
//! - `ELEMENTS` row i (Z = i + 1, which `emit` asserts): atomic number, symbol and name (a
//!   NUL-separated `u8` blob plus a `u16` offset column each: baked arrays can't hold string
//!   literals, port plan §6.4), atomic mass, electronegativity (a `bool` presence column plus the
//!   value, 0.0 where Rust has `None`: a compile-time f64 can't be NaN) and category.
//! - `ELEMENT_PROPERTIES` row k: atomic number, IE₁, EA, covalent radius, and the oxidation
//!   states flattened into start / count columns over one `i8` member table (probed: `i8`
//!   arrays sign-extend on read).
//!
//! Every value is read through the oracle's public API (`ELEMENTS`, `ELEMENT_PROPERTIES`). The
//! Rust source is read only for comments: the doc comments of the enum, its variants, the
//! statics and the struct fields, and the `//` notes on individual cells (the IUPAC conventional
//! and CIAAW 2024 atomic weights), which travel as the cells' label comments.

use crate::cyr::{Elem, doc_comment, doc_comment_after, header, snake, string_column, table};
use crate::m_element::CATEGORIES;
use kimiya::element::{ELEMENT_PROPERTIES, ELEMENTS, lookup_by_number};
use std::collections::BTreeMap;
use std::fmt::Write as _;

fn member(variant: &str) -> String {
    format!("KM_ELEMENT_CATEGORY_{}", snake(variant))
}

/// The variant names of `pub enum ElementCategory` in the Rust source, in order. Names only:
/// this checks that `CATEGORIES` is complete and in declaration order.
fn source_variants(src: &str) -> Result<Vec<String>, String> {
    let anchor = "pub enum ElementCategory {";
    let at = src.find(anchor).ok_or("no `pub enum ElementCategory {`")?;
    let body = &src[at + anchor.len()..];
    let body = &body[..body
        .find("\n}")
        .ok_or("ElementCategory: no closing brace")?];
    Ok(body
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with("//") && !l.starts_with("#["))
        .map(|l| l.trim_end_matches(',').to_string())
        .collect())
}

/// The `//` notes inside the `ELEMENTS` literal: (atomic number, field) → the note, for the
/// field line the note precedes.
fn cell_notes(src: &str) -> Result<BTreeMap<(u8, String), String>, String> {
    let at = src
        .find("pub static ELEMENTS:")
        .ok_or("no `pub static ELEMENTS:`")?;
    let body = &src[at..];
    let body = &body[..body.find("\n];").ok_or("ELEMENTS: no closing `];`")?];
    let mut notes = BTreeMap::new();
    let mut z: Option<u8> = None;
    let mut pending: Vec<String> = Vec::new();
    for line in body.lines().map(str::trim) {
        if let Some(note) = line.strip_prefix("// ") {
            pending.push(note.to_string());
            continue;
        }
        let Some((field, value)) = line.split_once(": ") else {
            continue;
        };
        if field == "atomic_number" {
            let n = value.trim_end_matches(',');
            z = Some(
                n.parse()
                    .map_err(|e| format!("ELEMENTS: atomic_number {n}: {e}"))?,
            );
        }
        if !pending.is_empty() {
            let zz = z.ok_or("ELEMENTS: a note before the first atomic_number")?;
            notes.insert((zz, field.to_string()), pending.join(" "));
            pending.clear();
        }
    }
    Ok(notes)
}

/// The doc line(s) of `pub struct Element`'s field, then which Rust field the column mirrors.
fn field_doc(src: &str, strukt: &str, field: &str, index: &str) -> Result<String, String> {
    let docs = doc_comment_after(
        src,
        &format!("pub struct {strukt} {{"),
        &format!("pub {field}:"),
    )?;
    let mut out = String::new();
    for l in docs.lines() {
        let _ = writeln!(out, "{}", l.trim_start_matches('#').trim_start());
    }
    let _ = write!(out, "{index}; mirrors `{strukt}::{field}`.");
    Ok(out)
}

#[allow(clippy::too_many_lines)] // one column after another, in field order
pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let variants: Vec<String> = CATEGORIES.iter().map(|c| format!("{c:?}")).collect();
    let in_source = source_variants(src)?;
    if variants != in_source {
        return Err(format!(
            "m_element::CATEGORIES {variants:?} is not ElementCategory in source order {in_source:?}"
        ));
    }
    for (k, el) in ELEMENTS.iter().enumerate() {
        if usize::from(el.atomic_number) != k + 1 {
            return Err(format!(
                "ELEMENTS[{k}] has atomic number {}",
                el.atomic_number
            ));
        }
    }
    let notes = cell_notes(src)?;
    let sym: Vec<String> = ELEMENTS.iter().map(|el| el.symbol.to_string()).collect();

    let mut out = header(
        "periodic-table data (Z = 1–118) and extended properties of 32 elements",
        provenance,
        &[
            "Element rows are indexed by Z − 1: a handle (src/element.cyr) is row + 1, which is Z.\n\
             Property rows are in Rust's order; a property handle is row + 1. KmElementCategory\n\
             follows the Rust declaration order and is append-only. f64 cells are the oracle's\n\
             values as round-trip decimals; tests/golden/element.tsv pins every cell.",
        ],
    );

    // ── ElementCategory ──
    out.push_str(&doc_comment(src, "pub enum ElementCategory {")?);
    out.push_str("enum KmElementCategory {\n");
    for (k, v) in variants.iter().enumerate() {
        let docs = doc_comment_after(src, "pub enum ElementCategory {", &format!("{v},"))?;
        for l in docs.lines() {
            let _ = writeln!(out, "    {l}");
        }
        let _ = writeln!(out, "    {} = {k};", member(v));
    }
    out.push_str("}\n");
    let _ = writeln!(
        out,
        "const KM_ELEMENT_CATEGORY_COUNT = {};\n",
        variants.len()
    );

    // ── ELEMENTS ──
    out.push_str(&doc_comment(src, "pub static ELEMENTS:")?);
    out.push_str("# Rust: kimiya::element::ELEMENTS (its length).\n");
    let _ = writeln!(out, "const KM_ELEMENT_COUNT = {};\n", ELEMENTS.len());

    let idx = "Indexed by Z − 1";
    let v = ELEMENTS.iter().map(|el| el.atomic_number).collect();
    out.push_str(&table(
        "KM_ELEMENT_ATOMIC_NUMBER",
        &field_doc(src, "Element", "atomic_number", idx)?,
        &Elem::U8(v),
        &sym,
    )?);

    let symbols: Vec<&str> = ELEMENTS.iter().map(|el| el.symbol).collect();
    out.push_str(&string_column(
        "_KM_ELEMENT_SYMBOL_BYTES",
        "_KM_ELEMENT_SYMBOL_OFF",
        &field_doc(src, "Element", "symbol", idx)?,
        &symbols,
    )?);
    let names: Vec<&str> = ELEMENTS.iter().map(|el| el.name).collect();
    out.push_str(&string_column(
        "_KM_ELEMENT_NAME_BYTES",
        "_KM_ELEMENT_NAME_OFF",
        &field_doc(src, "Element", "name", idx)?,
        &names,
    )?);

    let v = ELEMENTS.iter().map(|el| el.atomic_mass).collect();
    let mass_labels: Vec<String> = ELEMENTS
        .iter()
        .map(
            |el| match notes.get(&(el.atomic_number, "atomic_mass".to_string())) {
                Some(n) => format!("{}: {n}", el.symbol),
                None => el.symbol.to_string(),
            },
        )
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_ATOMIC_MASS",
        &field_doc(src, "Element", "atomic_mass", idx)?,
        &Elem::F64(v),
        &mass_labels,
    )?);

    let en_doc = field_doc(src, "Element", "electronegativity", idx)?;
    let v = ELEMENTS
        .iter()
        .map(|el| el.electronegativity.is_some())
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_HAS_ELECTRONEGATIVITY",
        &format!("{en_doc}\ntrue where Rust has Some(value), false for None."),
        &Elem::Bool(v),
        &sym,
    )?);
    let v = ELEMENTS
        .iter()
        .map(|el| el.electronegativity.unwrap_or(0.0))
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_ELECTRONEGATIVITY",
        &format!("{en_doc}\nThe value where KM_ELEMENT_HAS_ELECTRONEGATIVITY is true, else 0.0."),
        &Elem::F64(v),
        &sym,
    )?);

    let v = ELEMENTS
        .iter()
        .map(|el| member(&format!("{:?}", el.category)))
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_CATEGORY",
        &format!(
            "{} (KmElementCategory values).",
            field_doc(src, "Element", "category", idx)?.trim_end_matches('.')
        ),
        &Elem::Names(v),
        &sym,
    )?);

    // ── ELEMENT_PROPERTIES ──
    out.push_str(&doc_comment(src, "pub static ELEMENT_PROPERTIES:")?);
    out.push_str("# Rust: kimiya::element::ELEMENT_PROPERTIES (its length).\n");
    let _ = writeln!(
        out,
        "const KM_ELEMENT_PROPERTIES_COUNT = {};\n",
        ELEMENT_PROPERTIES.len()
    );
    let plab: Vec<String> = ELEMENT_PROPERTIES
        .iter()
        .map(|p| {
            lookup_by_number(p.atomic_number).map_or_else(
                || format!("Z = {}", p.atomic_number),
                |el| el.symbol.to_string(),
            )
        })
        .collect();
    let pidx = "Indexed by property row";
    let v = ELEMENT_PROPERTIES.iter().map(|p| p.atomic_number).collect();
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_ATOMIC_NUMBER",
        &field_doc(src, "ElementProperties", "atomic_number", pidx)?,
        &Elem::U8(v),
        &plab,
    )?);
    let v = ELEMENT_PROPERTIES
        .iter()
        .map(|p| p.ionization_energy_1)
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_IONIZATION_ENERGY_1",
        &field_doc(src, "ElementProperties", "ionization_energy_1", pidx)?,
        &Elem::F64(v),
        &plab,
    )?);
    let v = ELEMENT_PROPERTIES
        .iter()
        .map(|p| p.electron_affinity)
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_ELECTRON_AFFINITY",
        &field_doc(src, "ElementProperties", "electron_affinity", pidx)?,
        &Elem::F64(v),
        &plab,
    )?);
    let v = ELEMENT_PROPERTIES
        .iter()
        .map(|p| p.covalent_radius_pm)
        .collect();
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_COVALENT_RADIUS_PM",
        &field_doc(src, "ElementProperties", "covalent_radius_pm", pidx)?,
        &Elem::F64(v),
        &plab,
    )?);

    // oxidation_states is a slice per row: flatten it into start / count / members.
    let mut start = Vec::new();
    let mut count = Vec::new();
    let mut members = Vec::new();
    let mut member_labels = Vec::new();
    for (p, lab) in ELEMENT_PROPERTIES.iter().zip(&plab) {
        start.push(u8::try_from(members.len()).map_err(|e| e.to_string())?);
        count.push(u8::try_from(p.oxidation_states.len()).map_err(|e| e.to_string())?);
        for &state in p.oxidation_states {
            members.push(state);
            member_labels.push(lab.clone());
        }
    }
    let ox = field_doc(
        src,
        "ElementProperties",
        "oxidation_states",
        "Flattened by property row",
    )?;
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_OXIDATION_STATE_START",
        &format!(
            "{ox} Row k's states are\nKM_ELEMENT_PROPERTIES_OXIDATION_STATES[START[k] .. START[k] + COUNT[k]].\n\
             Start index into KM_ELEMENT_PROPERTIES_OXIDATION_STATES, by property row."
        ),
        &Elem::U8(start),
        &plab,
    )?);
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_OXIDATION_STATE_COUNT",
        "Number of oxidation states, by property row.",
        &Elem::U8(count),
        &plab,
    )?);
    out.push_str(&table(
        "KM_ELEMENT_PROPERTIES_OXIDATION_STATES",
        "Every row's oxidation states, in row order (each row's most common first).",
        &Elem::I8(members),
        &member_labels,
    )?);

    Ok(out.trim_end().to_string() + "\n")
}
