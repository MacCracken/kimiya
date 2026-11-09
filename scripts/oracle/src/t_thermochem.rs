//! `src/thermochem_data.cyr` — the thermochemical data of `rust-old/src/thermochem.rs`: the
//! `Shomate` struct, the 76 `THERMOCHEM_DATA` rows and the 34 `SHOMATE_DATA` segments, as baked
//! columns.
//!
//! - `struct KmShomate`: `Shomate`'s fields (a, b, c, d, e, t_min, t_max) in declaration order, all
//!   f64, with their Rust field docs. A by-value struct, as M2a's `KmAntoineCoeffs`.
//! - `THERMOCHEM_DATA` row k (a handle is k + 1, port plan §6.4): formula and name (a NUL-separated
//!   `u8` blob plus a `u16` offset column each) and ΔH_f°, ΔG_f°, S° (one f64 column each), read
//!   through the public static. `ThermochemData` is no Cyrius struct: its fields are the
//!   `thermochem_data_<field>(h)` accessors of `src/thermochem.cyr`.
//! - `SHOMATE_DATA` row k (a handle is k + 1): the formula (blob + offsets; a species with several
//!   segments repeats it, as the Rust tuples do) and the seven `Shomate` fields (one f64 column each).
//! - The `//` comments inside both literals (section headers, sources, the known defects) travel with
//!   each table's comment block, keyed by row, word-wrapped: they carry the citations and caveats
//!   (NIST WebBook segment ranges, the unsourced C₂H₆ / C₃H₈ / HF fits).
//!
//! thermochem has no `pub const` f64 (t_constants proves the list), so this file has no constant but
//! the two table lengths. Values come from the oracle's API; the Rust source is read for comments
//! only, and each literal's formulas are cross-checked against the API's rows in order.

use crate::cyr::{Elem, doc_comment, doc_comment_after, header, string_column, table};
use kimiya::thermochem::{SHOMATE_DATA, Shomate, THERMOCHEM_DATA};
use std::fmt::Write as _;

/// The `Shomate` fields, in declaration order.
pub const SHOMATE_FIELDS: [&str; 7] = ["a", "b", "c", "d", "e", "t_min", "t_max"];

/// The `ThermochemData` f64 fields, in declaration order (after `formula` and `name`).
pub const THERMOCHEM_FIELDS: [&str; 3] = ["delta_hf_kj", "delta_gf_kj", "s_standard_j"];

/// A segment's fields, in `SHOMATE_FIELDS` order.
pub fn shomate_fields(s: &Shomate) -> [f64; 7] {
    [s.a, s.b, s.c, s.d, s.e, s.t_min, s.t_max]
}

/// One row of a Rust table literal as the comment parser sees it.
#[derive(Default)]
struct RowComments {
    /// The row's key (its formula), read from the literal.
    key: String,
    /// A `// ── X ──` section header that starts at this row.
    section: Option<String>,
    /// Every other comment before the row (after the previous one) and inside it.
    notes: Vec<String>,
}

/// The comments of a table literal `pub static <NAME>:` … `];`, by row. A row opens on a line equal
/// to `open` and closes on one equal to `close` (both trimmed); its key is the first string literal
/// after `key_prefix` inside it. A comment outside a row belongs to the next row.
fn literal_comments(
    src: &str,
    name: &str,
    open: &str,
    close: &str,
    key_prefix: &str,
) -> Result<Vec<RowComments>, String> {
    let anchor = format!("pub static {name}:");
    let at = src.find(&anchor).ok_or_else(|| format!("no `{anchor}`"))?;
    let body = &src[at..];
    let body = &body[..body
        .find("\n];")
        .ok_or_else(|| format!("{name}: no closing `];`"))?];
    let mut rows: Vec<RowComments> = Vec::new();
    let mut pending = RowComments::default();
    let mut inside = false;
    for line in body.lines().skip(1).map(str::trim) {
        if let Some(c) = line.strip_prefix("//") {
            let c = c.trim();
            let target = if inside {
                rows.last_mut().ok_or_else(|| format!("{name}: a comment in no row"))?
            } else {
                &mut pending
            };
            if let Some(h) = c.strip_prefix("── ") {
                let title = h.trim_end_matches(['─', ' ']).to_string();
                if inside || target.section.is_some() {
                    return Err(format!("{name}: section header {title:?} inside a row or twice"));
                }
                target.section = Some(title);
            } else {
                target.notes.push(c.to_string());
            }
        } else if line == open {
            if inside {
                return Err(format!("{name}: a row opens inside a row"));
            }
            inside = true;
            rows.push(std::mem::take(&mut pending));
        } else if line == close {
            inside = false;
        } else if inside
            && let Some(r) = rows.last_mut()
            && r.key.is_empty()
            && let Some(rest) = line.strip_prefix(key_prefix)
            && let Some(k) = rest.strip_prefix('"').and_then(|x| x.split('"').next())
        {
            r.key = k.to_string();
        }
    }
    if !pending.notes.is_empty() || pending.section.is_some() {
        return Err(format!("{name}: comments after the last row"));
    }
    Ok(rows)
}

/// `text` greedily word-wrapped to lines of at most `width` characters.
fn wrap(text: &str, width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut cur = String::new();
    for word in text.split_whitespace() {
        if !cur.is_empty() && cur.chars().count() + 1 + word.chars().count() > width {
            lines.push(std::mem::take(&mut cur));
        }
        if !cur.is_empty() {
            cur.push(' ');
        }
        cur.push_str(word);
    }
    if !cur.is_empty() {
        lines.push(cur);
    }
    lines
}

/// One `# - …` list item, wrapped, continuation lines indented.
fn item(out: &mut String, text: &str) {
    for (k, line) in wrap(text, 92).iter().enumerate() {
        let lead = if k == 0 { "- " } else { "  " };
        let _ = writeln!(out, "# {lead}{line}");
    }
}

/// The sections and row notes of a literal, as `#` comment lines (empty if it has none).
fn comment_block(rows: &[RowComments], label: impl Fn(usize) -> String) -> String {
    let mut out = String::new();
    let starts: Vec<(usize, &str)> = rows
        .iter()
        .enumerate()
        .filter_map(|(k, r)| r.section.as_deref().map(|s| (k, s)))
        .collect();
    if !starts.is_empty() {
        out.push_str("#\n# Sections of the Rust literal:\n");
        for (j, (k, title)) in starts.iter().enumerate() {
            let end = starts.get(j + 1).map_or(rows.len(), |(n, _)| *n) - 1;
            item(&mut out, &format!("rows {k}–{end} ({} … {}): {title}", label(*k), label(end)));
        }
    }
    if rows.iter().any(|r| !r.notes.is_empty()) {
        out.push_str("#\n# Notes on rows, from the Rust literal:\n");
        for (k, r) in rows.iter().enumerate() {
            if !r.notes.is_empty() {
                item(&mut out, &format!("row {k}, {}: {}", label(k), r.notes.join(" ")));
            }
        }
    }
    out
}

/// A field's Rust doc as the plain lines of a column doc.
fn field_doc(src: &str, anchor: &str, field: &str) -> Result<String, String> {
    let docs = doc_comment_after(src, anchor, &format!("pub {field}:"))?;
    let mut doc = String::new();
    for l in docs.lines() {
        let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
    }
    Ok(doc)
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let tc = literal_comments(src, "THERMOCHEM_DATA", "ThermochemData {", "},", "formula: ")?;
    let tc_keys: Vec<&str> = tc.iter().map(|r| r.key.as_str()).collect();
    let api_keys: Vec<&str> = THERMOCHEM_DATA.iter().map(|d| d.formula).collect();
    if tc_keys != api_keys {
        return Err(format!(
            "THERMOCHEM_DATA: the literal's formulas {tc_keys:?} are not the API's {api_keys:?}"
        ));
    }
    let sh = literal_comments(src, "SHOMATE_DATA", "(", "),", "")?;
    let sh_keys: Vec<&str> = sh.iter().map(|r| r.key.as_str()).collect();
    let api_sh: Vec<&str> = SHOMATE_DATA.iter().map(|(f, _)| *f).collect();
    if sh_keys != api_sh {
        return Err(format!(
            "SHOMATE_DATA: the literal's formulas {sh_keys:?} are not the API's {api_sh:?}"
        ));
    }

    let mut out = header(
        "thermochemical data: formation enthalpies, Gibbs energies, entropies, Shomate segments",
        provenance,
        &[
            "THERMOCHEM_DATA and SHOMATE_DATA rows are in Rust's order; a handle (src/thermochem.cyr)\n\
             is row + 1, and 0 is None. f64 cells are the oracle's values as round-trip decimals;\n\
             tests/golden/thermochem.tsv pins every cell (the `thermochem::THERMOCHEM_DATA[]` and\n\
             `thermochem::SHOMATE_DATA[]` rows).",
        ],
    );

    // ── Shomate ──
    let anchor = "pub struct Shomate {";
    out.push_str(&doc_comment(src, anchor)?);
    out.push_str("struct KmShomate {\n");
    for field in SHOMATE_FIELDS {
        let docs = doc_comment_after(src, anchor, &format!("pub {field}:"))?;
        for l in docs.lines() {
            let _ = writeln!(out, "    {l}");
        }
        let _ = writeln!(out, "    {field}: f64;");
    }
    out.push_str("}\n\n");

    // ── THERMOCHEM_DATA ──
    out.push_str(&doc_comment(src, "pub static THERMOCHEM_DATA:")?);
    out.push_str(&comment_block(&tc, |k| THERMOCHEM_DATA[k].formula.to_string()));
    out.push_str("# Rust: kimiya::thermochem::THERMOCHEM_DATA (its length).\n");
    let _ = writeln!(out, "const KM_THERMOCHEM_COUNT = {};\n", THERMOCHEM_DATA.len());
    let labels: Vec<String> = api_keys.iter().map(|f| (*f).to_string()).collect();
    let anchor_t = "pub struct ThermochemData {";
    for (field, col) in [("formula", "FORMULA"), ("name", "NAME")] {
        let strings: Vec<&str> = THERMOCHEM_DATA
            .iter()
            .map(|d| if field == "formula" { d.formula } else { d.name })
            .collect();
        let doc = field_doc(src, anchor_t, field)?;
        out.push_str(&string_column(
            &format!("_KM_THERMOCHEM_{col}_BYTES"),
            &format!("_KM_THERMOCHEM_{col}_OFF"),
            &format!("{}Indexed by THERMOCHEM_DATA row; mirrors `ThermochemData::{field}`.", doc),
            &strings,
        )?);
    }
    for field in THERMOCHEM_FIELDS {
        let mut doc = field_doc(src, anchor_t, field)?;
        let _ = write!(doc, "Indexed by THERMOCHEM_DATA row; mirrors `ThermochemData::{field}`.");
        let col: Vec<f64> = THERMOCHEM_DATA
            .iter()
            .map(|d| match field {
                "delta_hf_kj" => d.delta_hf_kj,
                "delta_gf_kj" => d.delta_gf_kj,
                _ => d.s_standard_j,
            })
            .collect();
        out.push_str(&table(
            &format!("KM_THERMOCHEM_{}", field.to_ascii_uppercase()),
            &doc,
            &Elem::F64(col),
            &labels,
        )?);
    }

    // ── SHOMATE_DATA ──
    let seg_label = |k: usize| {
        let (f, s) = &SHOMATE_DATA[k];
        format!("{f} {}–{} K", s.t_min, s.t_max)
    };
    out.push_str(&doc_comment(src, "pub static SHOMATE_DATA:")?);
    out.push_str(&comment_block(&sh, seg_label));
    out.push_str("# Rust: kimiya::thermochem::SHOMATE_DATA (its length).\n");
    let _ = writeln!(out, "const KM_SHOMATE_COUNT = {};\n", SHOMATE_DATA.len());
    out.push_str(&string_column(
        "_KM_SHOMATE_FORMULA_BYTES",
        "_KM_SHOMATE_FORMULA_OFF",
        "The formula of each SHOMATE_DATA segment (`\"H2O(g)\"`, once per segment); mirrors the \
         tuple's `&str`.",
        &api_sh,
    )?);
    let seg_labels: Vec<String> = (0..SHOMATE_DATA.len()).map(seg_label).collect();
    for (k, field) in SHOMATE_FIELDS.iter().enumerate() {
        let mut doc = field_doc(src, anchor, field)?;
        let _ = write!(doc, "Indexed by SHOMATE_DATA row; mirrors `Shomate::{field}`.");
        let col: Vec<f64> = SHOMATE_DATA.iter().map(|(_, s)| shomate_fields(s)[k]).collect();
        out.push_str(&table(
            &format!("KM_SHOMATE_{}", field.to_ascii_uppercase()),
            &doc,
            &Elem::F64(col),
            &seg_labels,
        )?);
    }
    Ok(out.trim_end().to_string() + "\n")
}
