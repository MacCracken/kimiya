//! `src/electrochemistry_data.cyr` — the standard reduction potentials of `rust-old/src/electrochemistry.rs`: the
//! 49 `STANDARD_POTENTIALS` rows as baked columns.
//!
//! - `STANDARD_POTENTIALS` row k (a handle is k + 1, port plan §6.4): the equation (a NUL-separated **UTF-8** blob
//!   plus a `u16` offset column: "Cu²⁺ + 2e⁻ → Cu", its bytes copied exactly) and the couple (ASCII blob +
//!   offsets, "Cu2+/Cu"), `electrons` (a `u8` column) and `standard_potential` (an f64 column), read through the
//!   public static. `HalfReaction` is no Cyrius struct: its fields are the `half_reaction_<field>(h)` accessors of
//!   `src/electrochemistry.cyr`, and its Rust docs (the struct's and each field's) travel into the table's comment
//!   block and the column docs.
//! - The `//` comments inside the literal — the two "CRC Handbook, Electrochemical Series" sources, trailing the
//!   `standard_potential` of K⁺/K and V²⁺/V — travel with the table's comment block, keyed by row and field. The
//!   literal's couples are cross-checked against the API's rows, in order.
//! - `FARADAY` and `ELEMENTARY_CHARGE` are physical constants (`src/constants.cyr`, `t_constants.rs`), not here.
//!
//! Values come from the oracle's API; the Rust source is read for comments only.

use crate::cyr::{Elem, doc_comment, doc_comment_after, header, string_column, string_column_utf8, table};
use kimiya::electrochemistry::STANDARD_POTENTIALS;
use std::fmt::Write as _;

/// The `HalfReaction` fields, in declaration order.
pub const FIELDS: [&str; 4] = ["equation", "couple", "electrons", "standard_potential"];

/// One row of the literal: its couple, and the comments inside it as `field: text` (or the text alone for a
/// full-line comment).
struct RowNotes {
    couple: String,
    notes: Vec<String>,
}

/// The rows of `pub static STANDARD_POTENTIALS:` … `];`, with every `//` comment inside a `HalfReaction { … }`.
fn literal_notes(src: &str) -> Result<Vec<RowNotes>, String> {
    let anchor = "pub static STANDARD_POTENTIALS:";
    let at = src.find(anchor).ok_or_else(|| format!("no `{anchor}`"))?;
    let body = &src[at..];
    let body = &body[..body.find("\n];").ok_or("STANDARD_POTENTIALS: no closing `];`")?];
    let mut rows: Vec<RowNotes> = Vec::new();
    let mut inside = false;
    for line in body.lines().skip(1).map(str::trim) {
        if line == "HalfReaction {" {
            if inside {
                return Err("STANDARD_POTENTIALS: a row opens inside a row".to_string());
            }
            inside = true;
            rows.push(RowNotes { couple: String::new(), notes: Vec::new() });
            continue;
        }
        if line == "}," {
            inside = false;
            continue;
        }
        let (code, comment) = match line.split_once("//") {
            Some((c, n)) => (c.trim(), Some(n.trim())),
            None => (line, None),
        };
        if !inside {
            if comment.is_some() || !code.is_empty() {
                return Err(format!("STANDARD_POTENTIALS: {line:?} outside a row"));
            }
            continue;
        }
        let row = rows.last_mut().ok_or("STANDARD_POTENTIALS: no row")?;
        if let Some(rest) = code.strip_prefix("couple: \"") {
            row.couple = rest.split('"').next().unwrap_or_default().to_string();
        }
        if let Some(n) = comment {
            let field = code.split(':').next().unwrap_or_default().trim();
            row.notes.push(if field.is_empty() { n.to_string() } else { format!("{field}: {n}") });
        }
    }
    Ok(rows)
}

/// A field's Rust doc as the plain lines of a column doc.
fn field_doc(src: &str, field: &str) -> Result<String, String> {
    let docs = doc_comment_after(src, "pub struct HalfReaction {", &format!("pub {field}:"))?;
    let mut doc = String::new();
    for l in docs.lines() {
        let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
    }
    Ok(doc)
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let rows = literal_notes(src)?;
    let lit: Vec<&str> = rows.iter().map(|r| r.couple.as_str()).collect();
    let api: Vec<&str> = STANDARD_POTENTIALS.iter().map(|h| h.couple).collect();
    if lit != api {
        return Err(format!("STANDARD_POTENTIALS: the literal's couples {lit:?} are not the API's {api:?}"));
    }

    let mut out = header(
        "electrochemistry data: the standard reduction potentials",
        provenance,
        &[
            "STANDARD_POTENTIALS rows are in Rust's order (by E°, most negative first); a handle\n\
             (src/electrochemistry.cyr) is row + 1, and 0 is None. f64 cells are the oracle's values as\n\
             round-trip decimals; tests/golden/electrochemistry.tsv pins every cell (the\n\
             `electrochemistry::STANDARD_POTENTIALS[]` rows). FARADAY and ELEMENTARY_CHARGE are\n\
             src/constants.cyr's KM_FARADAY and KM_ELEMENTARY_CHARGE.",
        ],
    );

    // ── HalfReaction ──
    out.push_str(&doc_comment(src, "pub struct HalfReaction {")?);
    out.push_str("#\n# Rust's `HalfReaction` is no Cyrius struct: a row is a handle, its fields the accessors\n");
    out.push_str("# half_reaction_<field>(h) of src/electrochemistry.cyr:\n");
    for field in FIELDS {
        let doc = field_doc(src, field)?;
        for (k, line) in doc.lines().enumerate() {
            let lead = if k == 0 { format!("- {field}: ") } else { "  ".to_string() };
            let _ = writeln!(out, "# {lead}{line}");
        }
    }
    out.push('\n');

    // ── STANDARD_POTENTIALS ──
    out.push_str(&doc_comment(src, "pub static STANDARD_POTENTIALS:")?);
    if rows.iter().any(|r| !r.notes.is_empty()) {
        out.push_str("#\n# Notes on rows, from the Rust literal:\n");
        for (k, r) in rows.iter().enumerate() {
            for n in &r.notes {
                let _ = writeln!(out, "# - row {k}, {}: {n}", r.couple);
            }
        }
    }
    out.push_str("# Rust: kimiya::electrochemistry::STANDARD_POTENTIALS (its length).\n");
    let _ = writeln!(out, "const KM_STANDARD_POTENTIALS_COUNT = {};\n", STANDARD_POTENTIALS.len());
    let labels: Vec<String> = api.iter().map(|c| (*c).to_string()).collect();
    let equations: Vec<&str> = STANDARD_POTENTIALS.iter().map(|h| h.equation).collect();
    out.push_str(&string_column_utf8(
        "_KM_STANDARD_POTENTIALS_EQUATION_BYTES",
        "_KM_STANDARD_POTENTIALS_EQUATION_OFF",
        &format!(
            "{}Indexed by STANDARD_POTENTIALS row; mirrors `HalfReaction::equation`.",
            field_doc(src, "equation")?
        ),
        &equations,
    )?);
    out.push_str(&string_column(
        "_KM_STANDARD_POTENTIALS_COUPLE_BYTES",
        "_KM_STANDARD_POTENTIALS_COUPLE_OFF",
        &format!(
            "{}Indexed by STANDARD_POTENTIALS row; mirrors `HalfReaction::couple`.",
            field_doc(src, "couple")?
        ),
        &api,
    )?);
    let mut doc = field_doc(src, "electrons")?;
    let _ = write!(doc, "Indexed by STANDARD_POTENTIALS row; mirrors `HalfReaction::electrons` (a u8).");
    let electrons: Vec<u8> = STANDARD_POTENTIALS.iter().map(|h| h.electrons).collect();
    out.push_str(&table("KM_STANDARD_POTENTIALS_ELECTRONS", &doc, &Elem::U8(electrons), &labels)?);
    let mut doc = field_doc(src, "standard_potential")?;
    let _ = write!(doc, "Indexed by STANDARD_POTENTIALS row; mirrors `HalfReaction::standard_potential`.");
    let e0: Vec<f64> = STANDARD_POTENTIALS.iter().map(|h| h.standard_potential).collect();
    out.push_str(&table("KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL", &doc, &Elem::F64(e0), &labels)?);
    Ok(out.trim_end().to_string() + "\n")
}
