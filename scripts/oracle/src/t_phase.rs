//! `src/phase_data.cyr` — the phase-equilibrium data of `rust-old/src/phase.rs`: the
//! `AntoineCoeffs` struct, the 8 `ANTOINE_DATA` sets as baked columns, the mmHg → Pa factor, and
//! the ebullioscopic / cryoscopic constants `KB_*` / `KF_*`.
//!
//! - `struct KmAntoineCoeffs`: `AntoineCoeffs`'s fields in declaration order, all f64.
//! - `ANTOINE_DATA` row k (a handle is k + 1, port plan §6.4): the name (a NUL-separated `u8` blob
//!   plus a `u16` offset column) and A, B, C, t_min_c, t_max_c (one f64 column each), read through
//!   the public static. The `//` note inside a row's literal (methanol's source) travels with the
//!   table's comment block, word-wrapped (a label that long would break the 120-column lint).
//! - `_KM_PA_PER_MMHG`: Rust's private `PA_PER_MMHG`, which the API exposes only as the factor
//!   `antoine_pressure_pa` applies. With A = B = 0, C = 1 at T = 0, `antoine_pressure_mmhg` is
//!   10⁰ = 1.0 exactly (checked here), so `antoine_pressure_pa` returns the factor itself. The golden
//!   row of that call is class A (`m_phase.rs`), so the factor is pinned bit for bit.
//! - `KB_*` / `KF_*`: `t_constants::consts` for `Home::Data("phase_data")`.
//!
//! Values come from the oracle's API; the Rust source is read for comments only (and the
//! factor's literal, as a cross-check of the API value).

use crate::cyr::{Elem, comment_line, doc_comment, doc_comment_after, header, string_column, table};
use crate::t_constants::{Home, consts};
use kimiya::phase::{ANTOINE_DATA, AntoineCoeffs, antoine_pressure_mmhg, antoine_pressure_pa};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

/// The `AntoineCoeffs` fields, in declaration order.
const FIELDS: [&str; 5] = ["a", "b", "c", "t_min_c", "t_max_c"];

/// The coefficient set's fields, in `FIELDS` order.
pub fn fields(c: &AntoineCoeffs) -> [f64; 5] {
    [c.a, c.b, c.c, c.t_min_c, c.t_max_c]
}

/// The `//` notes inside the `ANTOINE_DATA` literal, by the row name they follow.
fn row_notes(src: &str) -> Result<BTreeMap<String, String>, String> {
    let at = src
        .find("pub static ANTOINE_DATA:")
        .ok_or("no `pub static ANTOINE_DATA:`")?;
    let body = &src[at..];
    let body = &body[..body.find("\n];").ok_or("ANTOINE_DATA: no closing `];`")?];
    let mut notes: BTreeMap<String, String> = BTreeMap::new();
    let mut name: Option<String> = None;
    for line in body.lines().map(str::trim) {
        if let Some(n) = line.strip_prefix('"').and_then(|l| l.strip_suffix("\",")) {
            name = Some(n.to_string());
        } else if let Some(note) = line.strip_prefix("// ") {
            let n = name.clone().ok_or("ANTOINE_DATA: a note before the first name")?;
            let e = notes.entry(n).or_default();
            if !e.is_empty() {
                e.push(' ');
            }
            e.push_str(note);
        }
    }
    Ok(notes)
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

/// `PA_PER_MMHG` through the API: `antoine_pressure_pa` of the unit set (see the module comment).
pub fn pa_per_mmhg() -> Result<f64, String> {
    let unit = AntoineCoeffs {
        a: 0.0,
        b: 0.0,
        c: 1.0,
        t_min_c: 0.0,
        t_max_c: 0.0,
    };
    let one = antoine_pressure_mmhg(&unit, 0.0);
    if one.to_bits() != 1.0_f64.to_bits() {
        return Err(format!("antoine_pressure_mmhg(unit, 0) = {one:?}, not 1.0"));
    }
    Ok(antoine_pressure_pa(&unit, 0.0))
}

/// The literal of `const PA_PER_MMHG: f64 = …;` in the source, parsed (a cross-check only).
fn pa_literal(src: &str) -> Result<f64, String> {
    let line = src
        .lines()
        .find(|l| l.starts_with("const PA_PER_MMHG: f64 = "))
        .ok_or("no `const PA_PER_MMHG: f64 =`")?;
    let lit = line
        .trim_start_matches("const PA_PER_MMHG: f64 = ")
        .trim_end_matches(';')
        .replace('_', "");
    lit.parse().map_err(|e| format!("PA_PER_MMHG {lit}: {e}"))
}

pub fn emit(root: &Path, src: &str, provenance: &str) -> Result<String, String> {
    let names: Vec<&str> = ANTOINE_DATA.iter().map(|(n, _)| *n).collect();
    let notes = row_notes(src)?;
    for n in notes.keys() {
        if !names.contains(&n.as_str()) {
            return Err(format!("ANTOINE_DATA: a note for {n:?}, which is no row"));
        }
    }
    let pa = pa_per_mmhg()?;
    let lit = pa_literal(src)?;
    if pa.to_bits() != lit.to_bits() {
        return Err(format!(
            "PA_PER_MMHG: the API gives {pa:?}, the source says {lit:?}"
        ));
    }

    let mut out = header(
        "phase-equilibrium data: Antoine coefficients, the mmHg → Pa factor, K_b and K_f",
        provenance,
        &[
            "ANTOINE_DATA rows are in Rust's order; a handle (src/phase.cyr) is row + 1, and 0 is None.\n\
             f64 cells are the oracle's values as round-trip decimals; tests/golden/phase.tsv pins\n\
             every cell (the `phase::ANTOINE_DATA[]` rows) and tests/golden/constants.tsv the K_b / K_f\n\
             constants.",
        ],
    );

    // ── AntoineCoeffs ──
    out.push_str(&doc_comment(src, "pub struct AntoineCoeffs {")?);
    out.push_str("struct KmAntoineCoeffs {\n");
    for field in FIELDS {
        let docs = doc_comment_after(src, "pub struct AntoineCoeffs {", &format!("pub {field}:"))?;
        for l in docs.lines() {
            let _ = writeln!(out, "    {l}");
        }
        let _ = writeln!(out, "    {field}: f64;");
    }
    out.push_str("}\n\n");

    // ── PA_PER_MMHG (private) ──
    out.push_str(&doc_comment(src, "const PA_PER_MMHG:")?);
    out.push_str(
        "# Rust: kimiya::phase::PA_PER_MMHG (private), read as antoine_pressure_pa(A = B = 0, C = 1, T = 0):\n\
         # 10⁰ × the factor, the factor itself.\n",
    );
    let _ = writeln!(
        out,
        "const _KM_PA_PER_MMHG = {};\n",
        crate::cyr::f64_lit(pa)?
    );

    // ── ANTOINE_DATA ──
    out.push_str(&doc_comment(src, "pub static ANTOINE_DATA:")?);
    if !notes.is_empty() {
        out.push_str("#\n# Notes on rows, from the Rust literal:\n");
        for (name, note) in &notes {
            for (k, line) in wrap(&format!("{name}: {note}"), 92).iter().enumerate() {
                let lead = if k == 0 { "- " } else { "  " };
                let _ = writeln!(out, "# {lead}{line}");
            }
        }
    }
    out.push_str("# Rust: kimiya::phase::ANTOINE_DATA (its length).\n");
    let _ = writeln!(out, "const KM_ANTOINE_COUNT = {};\n", ANTOINE_DATA.len());
    let labels: Vec<String> = names.iter().map(|n| (*n).to_string()).collect();
    out.push_str(&string_column(
        "_KM_ANTOINE_NAME_BYTES",
        "_KM_ANTOINE_NAME_OFF",
        "The solvent name of each ANTOINE_DATA row (`\"water\"`); mirrors the tuple's `&str`.",
        &names,
    )?);
    for (k, field) in FIELDS.iter().enumerate() {
        let docs = doc_comment_after(src, "pub struct AntoineCoeffs {", &format!("pub {field}:"))?;
        let mut doc = String::new();
        for l in docs.lines() {
            let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
        }
        let _ = write!(
            doc,
            "Indexed by ANTOINE_DATA row; mirrors `AntoineCoeffs::{field}`."
        );
        let col: Vec<f64> = ANTOINE_DATA.iter().map(|(_, c)| fields(c)[k]).collect();
        out.push_str(&table(
            &format!("KM_ANTOINE_{}", field.to_ascii_uppercase()),
            &doc,
            &Elem::F64(col),
            &labels,
        )?);
    }

    // ── K_b, K_f ──
    let _ = writeln!(
        out,
        "{}",
        comment_line("Common ebullioscopic constants K_b and cryoscopic constants K_f (K·kg/mol).")
    );
    out.push('\n');
    out.push_str(&consts(root, Home::Data("phase_data"))?);
    Ok(out.trim_end().to_string() + "\n")
}
