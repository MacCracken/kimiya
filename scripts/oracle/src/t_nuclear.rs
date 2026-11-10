//! `src/nuclear_data.cyr` — the isotope data of `rust-old/src/nuclear.rs`: the `ISOTOPES` table as baked columns.
//!
//! - `ISOTOPES` row k (a handle is k + 1, port plan §6.4): the nuclide label (ASCII blob + `u16` offsets, "C-14",
//!   "Tc-99m"), `z` (a `u8` column), `a` (a `u16` column), `half_life_s` (f64) and `decay_mode` (a NUL-separated
//!   **UTF-8** blob plus offsets: "β⁻/β⁺/EC", its bytes copied exactly), read through the public static, and
//!   `KM_ISOTOPES_COUNT`. `Isotope` is no Cyrius struct: its fields are the `isotope_<field>(h)` accessors of
//!   `src/nuclear.cyr`, and its Rust docs (the struct's and each field's) travel into the table's comment block and
//!   the column docs.
//! - The `//` comments inside the literal (C-14's 5700 a, I-131's 8.0252 d and Pu-239's 24 110 a, NUBASE2020) travel
//!   into the table's comment block, keyed by row and field; the literal's symbols are cross-checked against the API's
//!   rows, in order.
//! - `AMU_MEV`, `PROTON_MASS_U`, `NEUTRON_MASS_U` and `HYDROGEN_ATOM_MASS_U` are physical constants (`src/constants.cyr`,
//!   `t_constants.rs`), not here. Rust's private `EQUAL_RATE_REL_TOL` and `SEMF_*` coefficients are typed formula
//!   constants inline in `src/nuclear.cyr` (port plan §6.4), pinned by golden rows (`validate.py` perturbs each).
//!
//! Values come from the oracle's API; the Rust source is read for comments only.

use crate::cyr::{Elem, doc_comment, doc_comment_after, header, string_column, string_column_utf8, table};
use kimiya::nuclear::ISOTOPES;
use std::fmt::Write as _;

/// The `Isotope` fields, in declaration order.
pub const FIELDS: [&str; 5] = ["symbol", "z", "a", "half_life_s", "decay_mode"];

/// One row of the literal: its symbol, and the comments inside it as `field: text` (a full-line comment belongs to the
/// field line after it).
struct RowNotes {
    symbol: String,
    notes: Vec<String>,
}

/// The rows of `pub static ISOTOPES:` … `];`, with every `//` comment inside an `Isotope { … }`.
fn literal_notes(src: &str) -> Result<Vec<RowNotes>, String> {
    let anchor = "pub static ISOTOPES:";
    let at = src.find(anchor).ok_or_else(|| format!("no `{anchor}`"))?;
    let body = &src[at..];
    let body = &body[..body.find("\n];").ok_or("ISOTOPES: no closing `];`")?];
    let mut rows: Vec<RowNotes> = Vec::new();
    let mut pending: Vec<String> = Vec::new();
    let mut inside = false;
    for line in body.lines().skip(1).map(str::trim) {
        if line == "Isotope {" {
            if inside {
                return Err("ISOTOPES: a row opens inside a row".to_string());
            }
            inside = true;
            rows.push(RowNotes { symbol: String::new(), notes: Vec::new() });
            continue;
        }
        if line == "}," {
            inside = false;
            continue;
        }
        let (code, comment) = match line.split_once("//") {
            Some((c, n)) => (c.trim(), Some(n.trim().to_string())),
            None => (line, None),
        };
        if !inside {
            if comment.is_some() || !code.is_empty() {
                return Err(format!("ISOTOPES: {line:?} outside a row"));
            }
            continue;
        }
        let row = rows.last_mut().ok_or("ISOTOPES: no row")?;
        if code.is_empty() {
            pending.extend(comment);
            continue;
        }
        let field = code.split(':').next().unwrap_or_default().trim().to_string();
        if !FIELDS.contains(&field.as_str()) {
            return Err(format!("ISOTOPES: {field:?} is no Isotope field"));
        }
        if field == "symbol" {
            let value = code.split_once(':').map_or("", |(_, v)| v).trim().trim_end_matches(',');
            row.symbol = value.trim_matches('"').to_string();
        }
        for n in pending.drain(..).chain(comment) {
            row.notes.push(format!("{field}: {n}"));
        }
    }
    if !pending.is_empty() {
        return Err("ISOTOPES: a comment after the last field".to_string());
    }
    Ok(rows)
}

/// A field's Rust doc as the plain lines of a column doc.
fn field_doc(src: &str, field: &str) -> Result<String, String> {
    let docs = doc_comment_after(src, "pub struct Isotope {", &format!("pub {field}:"))?;
    let mut doc = String::new();
    for l in docs.lines() {
        let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
    }
    Ok(doc)
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let rows = literal_notes(src)?;
    let lit: Vec<&str> = rows.iter().map(|r| r.symbol.as_str()).collect();
    let api: Vec<&str> = ISOTOPES.iter().map(|i| i.symbol).collect();
    if lit != api {
        return Err(format!("ISOTOPES: the literal's symbols {lit:?} are not the API's {api:?}"));
    }
    let mut seen = std::collections::BTreeSet::new();
    if !api.iter().all(|s| seen.insert(*s)) {
        return Err(format!("ISOTOPES: a symbol is listed twice: {api:?}"));
    }

    let mut out = header(
        "nuclear data: the ISOTOPES table",
        provenance,
        &[
            "ISOTOPES rows are in Rust's order; a handle (src/nuclear.cyr) is row + 1, and 0 is None. f64 cells\n\
             are the oracle's values as round-trip decimals; tests/golden/nuclear.tsv pins every cell (the\n\
             `nuclear::ISOTOPES[]` rows) and tests/golden/constants.tsv the four nuclear constants, which are\n\
             src/constants.cyr's (KM_AMU_MEV, KM_PROTON_MASS_U, KM_NEUTRON_MASS_U, KM_HYDROGEN_ATOM_MASS_U).",
        ],
    );

    out.push_str(&doc_comment(src, "pub struct Isotope {")?);
    out.push_str("#\n# Rust's `Isotope` is no Cyrius struct: a row is a handle, its fields the accessors\n");
    out.push_str("# isotope_<field>(h) of src/nuclear.cyr:\n");
    for field in FIELDS {
        for (k, line) in field_doc(src, field)?.lines().enumerate() {
            let lead = if k == 0 { format!("- {field}: ") } else { "  ".to_string() };
            let _ = writeln!(out, "# {lead}{line}");
        }
    }
    out.push('\n');
    out.push_str(&doc_comment(src, "pub static ISOTOPES:")?);
    if rows.iter().any(|r| !r.notes.is_empty()) {
        out.push_str("#\n# Notes on rows, from the Rust literal:\n");
        for (k, r) in rows.iter().enumerate() {
            for n in &r.notes {
                let _ = writeln!(out, "# - row {k}, {}, {n}", r.symbol);
            }
        }
    }
    out.push_str("# Rust: kimiya::nuclear::ISOTOPES (its length).\n");
    let _ = writeln!(out, "const KM_ISOTOPES_COUNT = {};\n", ISOTOPES.len());

    let labels: Vec<String> = api.iter().map(|s| (*s).to_string()).collect();
    out.push_str(&string_column(
        "_KM_ISOTOPES_SYMBOL_BYTES",
        "_KM_ISOTOPES_SYMBOL_OFF",
        &format!("{}Indexed by ISOTOPES row; mirrors `Isotope::symbol`.", field_doc(src, "symbol")?),
        &api,
    )?);
    let z: Vec<u8> = ISOTOPES.iter().map(|i| i.z).collect();
    out.push_str(&table(
        "KM_ISOTOPES_Z",
        &format!("{}Indexed by ISOTOPES row; mirrors `Isotope::z` (a u8).", field_doc(src, "z")?),
        &Elem::U8(z),
        &labels,
    )?);
    let a: Vec<u16> = ISOTOPES.iter().map(|i| i.a).collect();
    out.push_str(&table(
        "KM_ISOTOPES_A",
        &format!("{}Indexed by ISOTOPES row; mirrors `Isotope::a` (a u16).", field_doc(src, "a")?),
        &Elem::U16(a),
        &labels,
    )?);
    let hl: Vec<f64> = ISOTOPES.iter().map(|i| i.half_life_s).collect();
    out.push_str(&table(
        "KM_ISOTOPES_HALF_LIFE_S",
        &format!("{}Indexed by ISOTOPES row; mirrors `Isotope::half_life_s`.", field_doc(src, "half_life_s")?),
        &Elem::F64(hl),
        &labels,
    )?);
    let modes: Vec<&str> = ISOTOPES.iter().map(|i| i.decay_mode).collect();
    out.push_str(&string_column_utf8(
        "_KM_ISOTOPES_DECAY_MODE_BYTES",
        "_KM_ISOTOPES_DECAY_MODE_OFF",
        &format!("{}Indexed by ISOTOPES row; mirrors `Isotope::decay_mode`.", field_doc(src, "decay_mode")?),
        &modes,
    )?);
    Ok(out.trim_end().to_string() + "\n")
}
