//! `src/inorganic_data.cyr` — the inorganic-chemistry data of `rust-old/src/inorganic.rs`: the spectrochemical series,
//! the octahedral occupancy and unpaired-electron tables, the Madelung constants and Born exponents, and the
//! `IONIC_RADII` table.
//!
//! - `enum KmLigand` (`KM_LIGAND_COUNT`): the Rust variants in declaration order, value k for the k-th (append-only),
//!   proved against the source declaration (`cyr::enum_decl`; `#[non_exhaustive]` keeps `m_inorganic::LIGANDS` from
//!   an exhaustive match), and `KM_LIGAND_FIELD_STRENGTH` (u8), `Ligand::field_strength` through the API.
//!   `is_strong_field` is a formula (`field_strength() >= 8`), so it stays code in `src/inorganic.cyr`.
//! - Rust's private `OCT_UNPAIRED_HIGH_SPIN` / `_LOW_SPIN` are `unpaired_electrons_oct(n, spin)` for n = 0 ..= 10
//!   (u8 columns `_KM_OCT_UNPAIRED_*`); the API returns Err from n = 11, which is checked, so the tables have 11 rows.
//!   The private `OCT_OCCUPANCY_HIGH_SPIN` / `_LOW_SPIN` (t₂g, e_g) pairs are read through `cfse_octahedral`: for each
//!   n and spin, exactly one pair with t₂g + e_g = n, t₂g ≤ 6 and e_g ≤ 4 gives `f64::from(t2g) * (-0.4) +
//!   f64::from(eg) * 0.6` with the API's bits (checked unique), as `PA_PER_MMHG` is read through
//!   `antoine_pressure_pa` (`t_phase`). The source literals of all four tables are cross-checked.
//! - `MADELUNG_*` / `BORN_EXP_*`: `t_constants::consts_where` for `Home::Data("inorganic_data")`, each group under the
//!   source's `//` section comment (the Madelung reference distance; Pauling's Born exponents by configuration).
//! - `IONIC_RADII` row k (a handle is k + 1, port plan §6.4): the symbol (ASCII blob + `u16` offsets), `charge` (`i8`,
//!   sign-extended on read), `coordination` (`u8`) and `radius_pm` (f64), read through the public static, and
//!   `KM_IONIC_RADII_COUNT`. `IonicRadius` is no Cyrius struct: its fields are the `ionic_radius_<field>(h)` accessors
//!   of `src/inorganic.cyr`, and its Rust docs travel into the table's comment block. The `//` comment inside the
//!   literal (Sn²⁺'s CN 8) travels too, keyed by row and field; the literal's (symbol, charge) keys are cross-checked
//!   against the API's rows, in order.
//!
//! Values come from the oracle's API; the Rust source is read for comments only (and the private tables' literals,
//! as a cross-check of the API values).

use crate::cyr::{Elem, Member, comment_line, doc_comment, doc_comment_after, enum_decl, header, snake, string_column, table};
use crate::m_inorganic::LIGANDS;
use crate::t_constants::{Home, consts_where};
use kimiya::inorganic::{IONIC_RADII, cfse_octahedral, unpaired_electrons_oct};
use std::fmt::Write as _;
use std::path::Path;

/// The `IonicRadius` fields, in declaration order.
const FIELDS: [&str; 4] = ["symbol", "charge", "coordination", "radius_pm"];

/// The rows of an occupancy table, (t₂g, e_g) per n = 0 ..= 10, read through `cfse_octahedral` (module comment).
fn occupancy(high_spin: bool) -> Result<Vec<(u8, u8)>, String> {
    let mut out = Vec::new();
    for n in 0..=10u8 {
        let want = cfse_octahedral(n, high_spin).map_err(|e| format!("cfse_octahedral({n}, {high_spin}): {e:?}"))?;
        let hits: Vec<(u8, u8)> = (0..=6u8)
            .filter_map(|t2g| n.checked_sub(t2g).map(|eg| (t2g, eg)))
            .filter(|&(_, eg)| eg <= 4)
            .filter(|&(t2g, eg)| (f64::from(t2g) * (-0.4) + f64::from(eg) * 0.6).to_bits() == want.to_bits())
            .collect();
        match hits.as_slice() {
            [one] => out.push(*one),
            _ => return Err(format!("cfse_octahedral({n}, {high_spin}) = {want:?}: occupancies {hits:?}")),
        }
    }
    if cfse_octahedral(11, high_spin).is_ok() {
        return Err(format!("cfse_octahedral(11, {high_spin}) is Ok: the table is longer than 11 rows"));
    }
    Ok(out)
}

/// The rows of an unpaired-electron table, n = 0 ..= 10, through `unpaired_electrons_oct`.
fn unpaired(high_spin: bool) -> Result<Vec<u8>, String> {
    let out = (0..=10u8)
        .map(|n| unpaired_electrons_oct(n, high_spin).map_err(|e| format!("unpaired_electrons_oct({n}): {e:?}")))
        .collect::<Result<Vec<u8>, String>>()?;
    if unpaired_electrons_oct(11, high_spin).is_ok() {
        return Err(format!("unpaired_electrons_oct(11, {high_spin}) is Ok: the table is longer than 11 rows"));
    }
    Ok(out)
}

/// The integers of `const <name>: … = [ … ];` in the source, in order (a cross-check only).
fn literal_ints(src: &str, name: &str) -> Result<Vec<u8>, String> {
    let at = src.find(&format!("const {name}:")).ok_or_else(|| format!("no `const {name}:`"))?;
    let body = &src[at..];
    let body = &body[body.find("= [").ok_or_else(|| format!("{name}: no `= [`"))? + 3..];
    let body = &body[..body.find("];").ok_or_else(|| format!("{name}: no `];`"))?];
    body.split(|c: char| !c.is_ascii_digit())
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<u8>().map_err(|e| format!("{name}: {t}: {e}")))
        .collect()
}

/// The `//` lines (not `///`) of the comment block above the doc comment of the item whose line starts with `item`,
/// separated from it by blank lines: a section heading such as `// Madelung constants (dimensionless), …`.
fn section_comment(src: &str, item: &str) -> Result<Vec<String>, String> {
    let lines: Vec<&str> = src.lines().collect();
    let mut k = lines
        .iter()
        .position(|l| l.starts_with(item))
        .ok_or_else(|| format!("no `{item}`"))?;
    while k > 0 && lines[k - 1].trim_start().starts_with("///") {
        k -= 1;
    }
    while k > 0 && lines[k - 1].trim().is_empty() {
        k -= 1;
    }
    let mut out = Vec::new();
    while k > 0 {
        let l = lines[k - 1].trim_start();
        match l.strip_prefix("//") {
            Some(rest) if !l.starts_with("///") => out.push(rest.trim().to_string()),
            _ => break,
        }
        k -= 1;
    }
    if out.is_empty() {
        return Err(format!("no `//` section comment above `{item}`"));
    }
    out.reverse();
    Ok(out)
}

/// One row of the IONIC_RADII literal: its (symbol, charge) key, and the comments inside it as `field: text` (a
/// full-line comment belongs to the field line after it).
struct RowNotes {
    symbol: String,
    charge: i8,
    notes: Vec<String>,
}

fn literal_notes(src: &str) -> Result<Vec<RowNotes>, String> {
    let anchor = "pub static IONIC_RADII:";
    let at = src.find(anchor).ok_or_else(|| format!("no `{anchor}`"))?;
    let body = &src[at..];
    let body = &body[..body.find("\n];").ok_or("IONIC_RADII: no closing `];`")?];
    let mut rows: Vec<RowNotes> = Vec::new();
    let mut pending: Vec<String> = Vec::new();
    let mut inside = false;
    for line in body.lines().skip(1).map(str::trim) {
        if line == "IonicRadius {" {
            inside = true;
            rows.push(RowNotes { symbol: String::new(), charge: 0, notes: Vec::new() });
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
                return Err(format!("IONIC_RADII: {line:?} outside a row"));
            }
            continue;
        }
        let row = rows.last_mut().ok_or("IONIC_RADII: no row")?;
        if code.is_empty() {
            pending.extend(comment);
            continue;
        }
        let field = code.split(':').next().unwrap_or_default().trim().to_string();
        let value = code.split_once(':').map_or("", |(_, v)| v).trim().trim_end_matches(',');
        if field == "symbol" {
            row.symbol = value.trim_matches('"').to_string();
        }
        if field == "charge" {
            row.charge = value.parse().map_err(|e| format!("IONIC_RADII: charge {value}: {e}"))?;
        }
        for n in pending.drain(..).chain(comment) {
            row.notes.push(format!("{field}: {n}"));
        }
    }
    if !pending.is_empty() {
        return Err("IONIC_RADII: a comment after the last field".to_string());
    }
    Ok(rows)
}

/// A field's Rust doc as the plain lines of a column doc.
fn field_doc(src: &str, field: &str) -> Result<String, String> {
    let docs = doc_comment_after(src, "pub struct IonicRadius {", &format!("pub {field}:"))?;
    let mut doc = String::new();
    for l in docs.lines() {
        let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
    }
    Ok(doc)
}

#[allow(clippy::too_many_lines)] // one table after another, in source order
pub fn emit(root: &Path, src: &str, provenance: &str) -> Result<String, String> {
    let ligands: Vec<String> = LIGANDS.iter().map(|l| format!("{l:?}")).collect();
    let rows = literal_notes(src)?;
    let lit: Vec<(String, i8)> = rows.iter().map(|r| (r.symbol.clone(), r.charge)).collect();
    let api: Vec<(String, i8)> = IONIC_RADII.iter().map(|r| (r.symbol.to_string(), r.charge)).collect();
    if lit != api {
        return Err(format!("IONIC_RADII: the literal's keys {lit:?} are not the API's {api:?}"));
    }

    let mut out = header(
        "inorganic-chemistry data: ligands, d-electron tables, Madelung and Born constants, ionic radii",
        provenance,
        &[
            "KmLigand follows the Rust declaration order and is append-only; KM_LIGAND_FIELD_STRENGTH is indexed by\n\
             its value. The d-electron tables are indexed by n = 0 ..= 10. IONIC_RADII rows are in Rust's order; a\n\
             handle (src/inorganic.cyr) is row + 1, and 0 is None. f64 cells are the oracle's values as round-trip\n\
             decimals; tests/golden/inorganic.tsv pins every cell (the tables through cfse_octahedral,\n\
             unpaired_electrons_oct and the `inorganic::IONIC_RADII[]` rows) and tests/golden/constants.tsv the\n\
             Madelung constants and Born exponents.",
        ],
    );

    // ── Ligand ──
    let members: Vec<Member> = ligands
        .iter()
        .map(|v| Member { variant: v.clone(), name: format!("KM_LIGAND_{}", snake(v)), extra: Vec::new() })
        .collect();
    out.push_str(&enum_decl(src, "Ligand", "KmLigand", "KM_LIGAND_COUNT", &members)?);
    let strength: Vec<u8> = LIGANDS.iter().map(|l| l.field_strength()).collect();
    let mut doc = String::new();
    for l in doc_comment(src, "pub const fn field_strength(")?.lines() {
        let _ = writeln!(doc, "{}", l.trim_start_matches('#').trim_start());
    }
    let _ = write!(doc, "Indexed by KmLigand value; mirrors `Ligand::field_strength` (a u8).");
    out.push_str(&table("KM_LIGAND_FIELD_STRENGTH", &doc, &Elem::U8(strength), &ligands)?);

    // ── the octahedral d-electron tables (private in Rust) ──
    let nlab: Vec<String> = (0..=10).map(|n| format!("d{n}")).collect();
    for (rust, high) in [("OCT_OCCUPANCY_HIGH_SPIN", true), ("OCT_OCCUPANCY_LOW_SPIN", false)] {
        let occ = occupancy(high)?;
        let flat: Vec<u8> = occ.iter().flat_map(|&(t, e)| [t, e]).collect();
        let literal = literal_ints(src, rust)?;
        if flat != literal {
            return Err(format!("{rust}: through cfse_octahedral {occ:?}, the source says {literal:?}"));
        }
        let doc = doc_comment(src, &format!("const {rust}:"))?;
        for (k, part) in ["t₂g", "e_g"].iter().enumerate() {
            let col: Vec<u8> = occ.iter().map(|p| if k == 0 { p.0 } else { p.1 }).collect();
            let suffix = if k == 0 { "T2G" } else { "EG" };
            let mut d = String::new();
            for l in doc.lines() {
                let _ = writeln!(d, "{}", l.trim_start_matches('#').trim_start());
            }
            let _ = write!(
                d,
                "The {part} count, indexed by n. Rust: kimiya::inorganic::{rust} (private), read through cfse_octahedral."
            );
            out.push_str(&table(&format!("_KM_{rust}_{suffix}"), &d, &Elem::U8(col), &nlab)?);
        }
    }
    for (rust, high) in [("OCT_UNPAIRED_HIGH_SPIN", true), ("OCT_UNPAIRED_LOW_SPIN", false)] {
        let col = unpaired(high)?;
        let literal = literal_ints(src, rust)?;
        if col != literal {
            return Err(format!("{rust}: through unpaired_electrons_oct {col:?}, the source says {literal:?}"));
        }
        let mut d = String::new();
        for l in doc_comment(src, &format!("const {rust}:"))?.lines() {
            let _ = writeln!(d, "{}", l.trim_start_matches('#').trim_start());
        }
        let _ = write!(d, "Rust: kimiya::inorganic::{rust} (private), read through unpaired_electrons_oct.");
        out.push_str(&table(&format!("_KM_{rust}"), &d, &Elem::U8(col), &nlab)?);
    }

    // ── Madelung constants, Born exponents ──
    let madelung = consts_where(root, Home::Data("inorganic_data"), |n| n.starts_with("MADELUNG_"))?;
    let born = consts_where(root, Home::Data("inorganic_data"), |n| n.starts_with("BORN_EXP_"))?;
    let all = consts_where(root, Home::Data("inorganic_data"), |_| true)?;
    if madelung.len() + born.len() != all.len() {
        return Err("inorganic_data: a constant that is neither MADELUNG_* nor BORN_EXP_*".to_string());
    }
    for (group, text) in [("pub const MADELUNG_NACL:", madelung), ("pub const BORN_EXP_HE:", born)] {
        for l in section_comment(src, group)? {
            let _ = writeln!(out, "{}", comment_line(&l));
        }
        out.push('\n');
        out.push_str(&text);
    }

    // ── IonicRadius, IONIC_RADII ──
    out.push_str(&doc_comment(src, "pub struct IonicRadius {")?);
    out.push_str("#\n# Rust's `IonicRadius` is no Cyrius struct: a row is a handle, its fields the accessors\n");
    out.push_str("# ionic_radius_<field>(h) of src/inorganic.cyr:\n");
    for field in FIELDS {
        for (k, line) in field_doc(src, field)?.lines().enumerate() {
            let lead = if k == 0 { format!("- {field}: ") } else { "  ".to_string() };
            let _ = writeln!(out, "# {lead}{line}");
        }
    }
    out.push('\n');
    out.push_str(&doc_comment(src, "pub static IONIC_RADII:")?);
    if rows.iter().any(|r| !r.notes.is_empty()) {
        out.push_str("#\n# Notes on rows, from the Rust literal:\n");
        for (k, r) in rows.iter().enumerate() {
            for n in &r.notes {
                let _ = writeln!(out, "# - row {k}, {}{:+}, {n}", r.symbol, r.charge);
            }
        }
    }
    out.push_str("# Rust: kimiya::inorganic::IONIC_RADII (its length).\n");
    let _ = writeln!(out, "const KM_IONIC_RADII_COUNT = {};\n", IONIC_RADII.len());
    let labels: Vec<String> = IONIC_RADII.iter().map(|r| format!("{}{:+}", r.symbol, r.charge)).collect();
    let symbols: Vec<&str> = IONIC_RADII.iter().map(|r| r.symbol).collect();
    out.push_str(&string_column(
        "_KM_IONIC_RADII_SYMBOL_BYTES",
        "_KM_IONIC_RADII_SYMBOL_OFF",
        &format!("{}Indexed by IONIC_RADII row; mirrors `IonicRadius::symbol`.", field_doc(src, "symbol")?),
        &symbols,
    )?);
    let charge: Vec<i8> = IONIC_RADII.iter().map(|r| r.charge).collect();
    out.push_str(&table(
        "KM_IONIC_RADII_CHARGE",
        &format!("{}Indexed by IONIC_RADII row; mirrors `IonicRadius::charge` (an i8).", field_doc(src, "charge")?),
        &Elem::I8(charge),
        &labels,
    )?);
    let cn: Vec<u8> = IONIC_RADII.iter().map(|r| r.coordination).collect();
    out.push_str(&table(
        "KM_IONIC_RADII_COORDINATION",
        &format!(
            "{}Indexed by IONIC_RADII row; mirrors `IonicRadius::coordination` (a u8).",
            field_doc(src, "coordination")?
        ),
        &Elem::U8(cn),
        &labels,
    )?);
    let radius: Vec<f64> = IONIC_RADII.iter().map(|r| r.radius_pm).collect();
    out.push_str(&table(
        "KM_IONIC_RADII_RADIUS_PM",
        &format!("{}Indexed by IONIC_RADII row; mirrors `IonicRadius::radius_pm`.", field_doc(src, "radius_pm")?),
        &Elem::F64(radius),
        &labels,
    )?);
    Ok(out.trim_end().to_string() + "\n")
}
