//! Cyrius source emission for the generated modules (port plan §3 step 6). Adapted from
//! tanmatra's `scripts/oracle/src/cyr.rs`.
//!
//! Rules the output must keep. Each was probed on the 6.7.18 toolchain (by tanmatra, and
//! re-checked here by `probe_tables.py`).
//! - **f64 values are decimal literals, written out positionally.** The lexer has no
//!   exponent notation. Decimals have been correctly rounded since 6.6.10, so the shortest
//!   round-trip digits give back the oracle's exact bits.
//! - **Hex is never used in a `const`**, where it is an integer and `CONST * x` then does
//!   integer arithmetic.
//! - **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** NaN can't be a compile-time f64 (the
//!   compiler refuses it: its bits depend on the target), so a NaN is an error here.
//! - **Every global gets the `KM_` prefix.** A bundle is concatenated into its consumers, so
//!   names share one namespace.

use std::fmt::Write as _;

/// The Cyrius literal for `x`, checked to parse back to the same bits.
pub fn f64_lit(x: f64) -> Result<String, String> {
    if x.is_nan() {
        return Err("NaN cannot be a compile-time f64 in Cyrius".to_string());
    }
    if x.is_infinite() {
        return Ok(if x > 0.0 {
            "(1.0 / 0.0)"
        } else {
            "-(1.0 / 0.0)"
        }
        .to_string());
    }
    let lit = positional(x);
    if lit.contains(['e', 'E', 'x', 'X']) {
        return Err(format!("{lit}: not a positional decimal"));
    }
    let back: f64 = lit.parse().map_err(|e| format!("{lit}: {e}"))?;
    if back.to_bits() != x.to_bits() {
        return Err(format!("{lit} does not round-trip to {:016x}", x.to_bits()));
    }
    Ok(lit)
}

/// The shortest round-trip digits of `x` (Rust's `{:e}`), written without an exponent.
fn positional(x: f64) -> String {
    let sci = format!("{:e}", x.abs());
    let (mant, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i64 = exp.parse().unwrap_or(0);
    let digits: String = mant.chars().filter(char::is_ascii_digit).collect();
    let point = exp + 1; // digits before the decimal point
    let n = i64::try_from(digits.len()).unwrap_or(i64::MAX);
    let body = if point <= 0 {
        let zeros = "0".repeat(usize::try_from(-point).unwrap_or(0));
        format!("0.{zeros}{digits}")
    } else if point >= n {
        let zeros = "0".repeat(usize::try_from(point - n).unwrap_or(0));
        format!("{digits}{zeros}.0")
    } else {
        let at = usize::try_from(point).unwrap_or(0);
        format!("{}.{}", &digits[..at], &digits[at..])
    };
    if x.is_sign_negative() {
        format!("-{body}")
    } else {
        body
    }
}

/// The `///` doc block directly above the first line of `src` that starts with `item`
/// (attributes in between are skipped), as `# ` comment lines. Empty if there is none.
/// Only comments are read from the Rust source; values always come from the oracle's API.
pub fn doc_comment(src: &str, item: &str) -> Result<String, String> {
    doc_comment_after(src, "", item)
}

/// As [`doc_comment`], for the first line starting with `item` (after leading whitespace) that
/// follows the first line starting with `anchor`: a field inside `pub struct Element {`, a
/// variant inside `pub enum ElementCategory {`.
pub fn doc_comment_after(src: &str, anchor: &str, item: &str) -> Result<String, String> {
    let lines: Vec<&str> = src.lines().collect();
    let from = lines
        .iter()
        .position(|l| l.trim_start().starts_with(anchor))
        .ok_or_else(|| format!("no `{anchor}` in the Rust source"))?;
    let at = from
        + lines[from..]
            .iter()
            .position(|l| l.trim_start().starts_with(item))
            .ok_or_else(|| format!("no `{item}` after `{anchor}` in the Rust source"))?;
    let mut docs = Vec::new();
    let mut k = at;
    while k > 0 {
        k -= 1;
        let l = lines[k].trim_start();
        if l.starts_with("#[") {
            continue;
        }
        match l.strip_prefix("///") {
            Some(rest) => docs.push(rest.strip_prefix(' ').unwrap_or(rest).to_string()),
            None => break,
        }
    }
    docs.reverse();
    let mut out = String::new();
    for d in docs {
        let _ = writeln!(out, "{}", comment_line(d.trim_end()));
    }
    Ok(out)
}

/// `line` as a `#` comment, with no trailing space when it is empty (lint-clean).
pub fn comment_line(line: &str) -> String {
    if line.is_empty() {
        "#".to_string()
    } else {
        format!("# {line}")
    }
}

/// `AlkaliMetal` → `ALKALI_METAL`, `NobleGas` → `NOBLE_GAS`, `TAI` → `TAI`. An underscore goes
/// before a capital that follows a lowercase letter or digit, or that starts a new capitalised
/// word after a run of capitals (tanmatra's rule).
pub fn snake(camel: &str) -> String {
    let c: Vec<char> = camel.chars().collect();
    let mut out = String::new();
    for (k, &ch) in c.iter().enumerate() {
        if ch.is_ascii_uppercase() && k > 0 {
            let prev = c[k - 1];
            let next_lower = c.get(k + 1).is_some_and(char::is_ascii_lowercase);
            if prev.is_ascii_lowercase()
                || prev.is_ascii_digit()
                || (prev.is_ascii_uppercase() && next_lower)
            {
                out.push('_');
            }
        }
        out.push(ch.to_ascii_uppercase());
    }
    out
}

/// A Cyrius array element type and how one value is written.
pub enum Elem {
    F64(Vec<f64>),
    U8(Vec<u8>),
    U16(Vec<u16>),
    I8(Vec<i8>),
    Bool(Vec<bool>),
    /// Names of constants (enum members), stored in a `u8` table.
    Names(Vec<String>),
}

/// `var NAME: T[N] = { .. };`, one element per line with a trailing label comment, under the
/// `doc` lines as `#` comments.
pub fn table(name: &str, doc: &str, elem: &Elem, labels: &[String]) -> Result<String, String> {
    table_sectioned(name, doc, elem, labels, &[])
}

/// As [`table`], with comment lines inside the initializer: `sections[i]` (when present and not
/// empty) is written as `# <line>` lines before element i. M4b: the `// Carbon bonds` headings of a
/// Rust `match` travel to the rows they head.
pub fn table_sectioned(
    name: &str,
    doc: &str,
    elem: &Elem,
    labels: &[String],
    sections: &[Vec<String>],
) -> Result<String, String> {
    let (ty, vals): (&str, Vec<String>) = match elem {
        Elem::F64(v) => (
            "f64",
            v.iter().map(|&x| f64_lit(x)).collect::<Result<_, _>>()?,
        ),
        Elem::U8(v) => ("u8", v.iter().map(u8::to_string).collect()),
        Elem::U16(v) => ("u16", v.iter().map(u16::to_string).collect()),
        Elem::I8(v) => ("i8", v.iter().map(i8::to_string).collect()),
        Elem::Bool(v) => ("bool", v.iter().map(bool::to_string).collect()),
        Elem::Names(v) => ("u8", v.clone()),
    };
    if vals.len() != labels.len() {
        return Err(format!(
            "{name}: {} values, {} labels",
            vals.len(),
            labels.len()
        ));
    }
    let mut out = String::new();
    for line in doc.lines() {
        let _ = writeln!(out, "{}", comment_line(line));
    }
    let _ = writeln!(out, "var {name}: {ty}[{}] = {{", vals.len());
    for (k, (v, label)) in vals.iter().zip(labels).enumerate() {
        for line in sections.get(k).map_or(&[][..], Vec::as_slice) {
            let _ = writeln!(out, "    {}", comment_line(line));
        }
        let _ = writeln!(out, "    {v},  # {label}");
    }
    out.push_str("};\n\n");
    Ok(out)
}

/// A string column as one NUL-terminated `u8` blob plus a `u16` offset column (port plan §6.4:
/// baked arrays can't hold string literals). `&BYTES + OFF[i]` is row i's C string. One blob
/// line per row, its bytes then the text as the label. ASCII only, and no NUL inside a string.
pub fn string_column(
    bytes_name: &str,
    off_name: &str,
    doc: &str,
    strings: &[&str],
) -> Result<String, String> {
    blob_column(bytes_name, off_name, doc, strings, false)
}

/// As [`string_column`], for UTF-8 text (M3a: the half-reaction equations, "Cu²⁺ + 2e⁻ → Cu"): the blob
/// holds each string's UTF-8 bytes exactly, so `strlen` counts bytes and a C-string compare is a byte
/// compare, as Rust's `&str ==` is. No NUL inside a string.
pub fn string_column_utf8(
    bytes_name: &str,
    off_name: &str,
    doc: &str,
    strings: &[&str],
) -> Result<String, String> {
    blob_column(bytes_name, off_name, doc, strings, true)
}

fn blob_column(
    bytes_name: &str,
    off_name: &str,
    doc: &str,
    strings: &[&str],
    utf8: bool,
) -> Result<String, String> {
    let mut offsets = Vec::new();
    let mut lines = Vec::new();
    let mut at = 0usize;
    for s in strings {
        if (!utf8 && !s.is_ascii()) || s.contains('\0') {
            let what = if utf8 { "NUL-free" } else { "NUL-free ASCII" };
            return Err(format!("{bytes_name}: {s:?} is not {what}"));
        }
        offsets.push(u16::try_from(at).map_err(|_| format!("{off_name}: offset {at} > u16"))?);
        let codes: Vec<String> = s.bytes().chain([0]).map(|b| b.to_string()).collect();
        let one = format!("    {},  # {s}", codes.join(", "));
        if one.chars().count() <= 120 {
            lines.push(one);
        } else {
            // A long string (thermochem's "Phosphorus pentachloride") would break the 120-column
            // lint: 16 bytes per line, the text as the first line's label.
            for (k, chunk) in codes.chunks(16).enumerate() {
                let label = if k == 0 { format!("  # {s}") } else { String::new() };
                lines.push(format!("    {},{label}", chunk.join(", ")));
            }
        }
        at += s.len() + 1;
    }
    let labels: Vec<String> = strings.iter().map(|s| (*s).to_string()).collect();
    let mut out = table(
        off_name,
        &format!("{doc}\nOffset of row i's C string in {bytes_name}."),
        &Elem::U16(offsets),
        &labels,
    )?;
    let enc = if utf8 { "UTF-8 bytes" } else { "ASCII" };
    let _ = writeln!(
        out,
        "{}\nvar {bytes_name}: u8[{at}] = {{",
        comment_line(&format!("Every string, NUL-terminated, in row order ({enc}).")),
    );
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("};\n\n");
    Ok(out)
}

/// The fixed header of every generated `.cyr` file.
pub fn header(title: &str, provenance: &str, notes: &[&str]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# kimiya — {title}");
    out.push_str("#\n");
    let _ = writeln!(
        out,
        "# GENERATED by scripts/oracle from {provenance}. Do not edit."
    );
    out.push_str(
        "# Regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml\n",
    );
    for n in notes {
        out.push_str("#\n");
        for line in n.lines() {
            let _ = writeln!(out, "# {line}");
        }
    }
    out.push('\n');
    out
}

/// One variant of a Rust enum as the source declares it (M4b).
pub struct SourceVariant {
    /// The `//` comment lines between the previous variant and this one (a section heading such as
    /// `── Carbon bonds ───`), without the `//`.
    pub sections: Vec<String>,
    /// The variant's name.
    pub name: String,
}

/// The unit variants of `pub enum <name> {` in the Rust source, in declaration order, each with the
/// `//` comments before it (doc comments and attributes are skipped). Names only: a generator proves
/// its variant list against this, since a `#[non_exhaustive]` enum can't be matched exhaustively
/// from outside its crate.
pub fn enum_variants(src: &str, name: &str) -> Result<Vec<SourceVariant>, String> {
    let anchor = format!("pub enum {name} {{");
    let at = src.find(&anchor).ok_or_else(|| format!("no `{anchor}`"))?;
    let body = &src[at + anchor.len()..];
    let body = &body[..body
        .find("\n}")
        .ok_or_else(|| format!("{name}: no closing brace"))?];
    let mut out = Vec::new();
    let mut sections = Vec::new();
    for l in body.lines().map(str::trim) {
        if l.is_empty() || l.starts_with("///") || l.starts_with("#[") {
            continue;
        }
        if let Some(c) = l.strip_prefix("//") {
            sections.push(c.trim().to_string());
            continue;
        }
        let v = l
            .strip_suffix(',')
            .filter(|v| !v.is_empty() && v.chars().all(|c| c.is_ascii_alphanumeric()))
            .ok_or_else(|| format!("{name}: {l:?} is no unit variant"))?;
        out.push(SourceVariant {
            sections: std::mem::take(&mut sections),
            name: v.to_string(),
        });
    }
    Ok(out)
}

/// One arm `Self::<Variant> => <value>,` of a `match self` in the Rust source (M4b).
pub struct SourceArm {
    /// The `//` comment lines between the previous arm and this one (`Carbon bonds`).
    pub sections: Vec<String>,
    /// The variant the arm matches.
    pub variant: String,
    /// The arm's trailing `// note`, if any (`~117° due to lone pair compression`).
    pub note: Option<String>,
}

/// The arms of the `match self {` in the body of the first item whose line starts with `item` (after
/// leading whitespace), e.g. `pub const fn energy_kj(`. Comments only: the values come from the API.
pub fn match_arms(src: &str, item: &str) -> Result<Vec<SourceArm>, String> {
    let at = src
        .lines()
        .position(|l| l.trim_start().starts_with(item))
        .ok_or_else(|| format!("no `{item}` in the Rust source"))?;
    let lines: Vec<&str> = src.lines().skip(at).collect();
    let open = lines
        .iter()
        .position(|l| l.trim() == "match self {")
        .ok_or_else(|| format!("{item}: no `match self {{`"))?;
    let mut arms = Vec::new();
    let mut sections = Vec::new();
    for l in lines[open + 1..].iter().map(|l| l.trim()) {
        if l == "}" {
            return Ok(arms);
        }
        if l.is_empty() {
            continue;
        }
        if let Some(c) = l.strip_prefix("//") {
            sections.push(c.trim().to_string());
            continue;
        }
        let rest = l
            .strip_prefix("Self::")
            .ok_or_else(|| format!("{item}: {l:?} is no `Self::X => …` arm"))?;
        let (variant, after) = rest
            .split_once(" =>")
            .ok_or_else(|| format!("{item}: {l:?} has no `=>`"))?;
        arms.push(SourceArm {
            sections: std::mem::take(&mut sections),
            variant: variant.to_string(),
            note: after.split_once("//").map(|(_, n)| n.trim().to_string()),
        });
    }
    Err(format!("{item}: the match has no closing brace"))
}

/// One member of a generated Cyrius enum: the Rust variant it mirrors, its Cyrius name, and doc lines
/// added below the variant's Rust doc (an ADR rename's note).
pub struct Member {
    pub variant: String,
    pub name: String,
    pub extra: Vec<String>,
}

/// `enum <cyr> { … }` then `const <count> = N;` for the Rust `pub enum <rust>`, under the enum's doc
/// comment: member k is `= k`, in declaration order (append-only), under the variant's `//` section
/// comments and doc comment. `members` must be the source's variants, in order.
pub fn enum_decl(
    src: &str,
    rust: &str,
    cyr: &str,
    count: &str,
    members: &[Member],
) -> Result<String, String> {
    let source = enum_variants(src, rust)?;
    let names: Vec<&str> = source.iter().map(|v| v.name.as_str()).collect();
    let listed: Vec<&str> = members.iter().map(|m| m.variant.as_str()).collect();
    if names != listed {
        return Err(format!(
            "{rust}: the generator lists {listed:?}, the source declares {names:?}"
        ));
    }
    let anchor = format!("pub enum {rust} {{");
    let mut out = doc_comment(src, &anchor)?;
    let _ = writeln!(out, "enum {cyr} {{");
    for (k, (m, sv)) in members.iter().zip(&source).enumerate() {
        for s in &sv.sections {
            let _ = writeln!(out, "    {}", comment_line(s));
        }
        let docs = doc_comment_after(src, &anchor, &format!("{},", m.variant))?;
        for l in docs.lines() {
            let _ = writeln!(out, "    {l}");
        }
        for l in &m.extra {
            let _ = writeln!(out, "    {}", comment_line(l));
        }
        let _ = writeln!(out, "    {} = {k};", m.name);
    }
    out.push_str("}\n");
    let _ = writeln!(out, "const {count} = {};\n", members.len());
    Ok(out)
}
