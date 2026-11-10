//! `src/spectroscopy_data.cyr` — the spectroscopy data of `rust-old/src/spectroscopy.rs`: the `SpectralSeries` enum and
//! its `n_lower` column, the isotope abundance table, and Rust's private Rydberg constant for ¹H.
//!
//! - `enum KmSpectralSeries` (`KM_SPECTRAL_SERIES_COUNT`): the Rust variants in declaration order, value k for the k-th
//!   (append-only), proved against the source declaration (`cyr::enum_decl`; `#[non_exhaustive]` keeps
//!   `m_spectroscopy::SERIES` from an exhaustive match), and `KM_SPECTRAL_SERIES_N_LOWER` (u8), `SpectralSeries::n_lower`
//!   through the API (its match arms checked to name every variant once).
//! - The isotope abundances: `isotope_abundances(z)` through the API for **every u8 z**, the non-empty slices
//!   concatenated as (mass, abundance) pairs into `_KM_ISOTOPE_ABUNDANCES` (f64, interleaved), with `_KM_ISOTOPE_
//!   ABUNDANCE_START` / `_COUNT` (u8, pair index and count) indexed by z for z below `_KM_ISOTOPE_ABUNDANCE_ZS` (one past
//!   the last z with data; every z from there to 255 is checked empty). The comments of the source's match arms (`// H,
//!   D (IUPAC 2009)`, `// 12C, 13C`) travel to their rows, the arms' z cross-checked against the API's non-empty set.
//! - Rust's private `RYDBERG_HYDROGEN` (= `RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO)`, const-evaluated) **through the
//!   API**, as `PA_PER_MMHG` is (`t_phase`): `SpectralSeries::Lyman.wavelength(u32::MAX − 1)` is 1/R_H (n₂ = u32::MAX, so
//!   1/1 − 1/n₂² is 1.0 exactly, checked), and R_H is the one double whose reciprocal rounds to it (checked unique: its
//!   neighbours' reciprocals differ). It is cross-checked against the source's const expression evaluated here (IEEE
//!   f64, as Rust's const evaluation is) with the private `ELECTRON_PROTON_MASS_RATIO` literal, which is generated too
//!   (`_KM_ELECTRON_PROTON_MASS_RATIO`, for the record and `tests/spectroscopy.tcyr`; the port reads R_H).
//! - `PLANCK`, `SPEED_OF_LIGHT`, `RYDBERG` and `RYDBERG_EV` are physical constants (`src/constants.cyr`, `t_constants.rs`).
//!   Rust's private `ISOTOPE_PRUNE_REL`, `ISOTOPE_REPORT_REL` and `ISOTOPE_BIN_WIDTH_U` are typed formula constants inline
//!   in `src/spectroscopy.cyr` (port plan §6.4), pinned by golden rows where a row can reach them and by
//!   `tests/spectroscopy.tcyr` on the private helpers.
//!
//! Values come from the oracle's API; the Rust source is read for comments only (and the private ratio's literal, as a
//! cross-check of R_H).

use crate::cyr::{Elem, Member, comment_line, doc_comment, enum_decl, f64_lit, header, match_arms, snake, table};
use crate::m_spectroscopy::SERIES;
use kimiya::spectroscopy::{RYDBERG, SpectralSeries, isotope_abundances};
use std::fmt::Write as _;

/// The value of the private `const <name>: f64 = <literal>;` in the source (underscores removed).
fn private_literal(src: &str, name: &str) -> Result<f64, String> {
    let anchor = format!("const {name}: f64 = ");
    let line = src
        .lines()
        .find(|l| l.trim_start().starts_with(&anchor))
        .ok_or_else(|| format!("no `{anchor}`"))?;
    let lit = line.trim_start()[anchor.len()..].trim_end().trim_end_matches(';').replace('_', "");
    lit.parse::<f64>().map_err(|e| format!("{name}: {lit:?}: {e}"))
}

/// Rust's private RYDBERG_HYDROGEN through the API (module comment), cross-checked against the const expression.
fn rydberg_hydrogen(ratio: f64) -> Result<f64, String> {
    let n1 = f64::from(SpectralSeries::Lyman.n_lower());
    let n2 = f64::from(u32::MAX);
    if 1.0 / (n1 * n1) - 1.0 / (n2 * n2) != 1.0 {
        return Err("1/1 − 1/u32::MAX² is not 1.0: the Lyman line does not isolate R_H".to_string());
    }
    let lam = SpectralSeries::Lyman.wavelength(u32::MAX - 1).map_err(|e| format!("Lyman line u32::MAX − 1: {e:?}"))?;
    let guess = 1.0 / lam;
    let mut x = guess;
    for _ in 0..4 {
        x = x.next_down();
    }
    let mut hits = Vec::new();
    for _ in 0..9 {
        if (1.0 / x).to_bits() == lam.to_bits() {
            hits.push(x);
        }
        x = x.next_up();
    }
    let [rh] = hits.as_slice() else {
        return Err(format!("RYDBERG_HYDROGEN: {} doubles have 1/x = {lam:?}: {hits:?}", hits.len()));
    };
    let lit = RYDBERG / (1.0 + ratio);
    if lit.to_bits() != rh.to_bits() {
        return Err(format!("RYDBERG_HYDROGEN: the API gives {rh:?}, RYDBERG / (1 + {ratio:?}) is {lit:?}"));
    }
    Ok(*rh)
}

/// The arms `<z> => &[ … ],  // note` of `isotope_abundances`' match: (z, note), in source order.
fn abundance_arms(src: &str) -> Result<Vec<(u8, String)>, String> {
    let at = src.find("pub fn isotope_abundances(").ok_or("no `pub fn isotope_abundances(`")?;
    let body = &src[at..];
    let body = &body[body.find("match atomic_number {").ok_or("isotope_abundances: no `match atomic_number {`")?..];
    let mut out = Vec::new();
    let mut open: Option<u8> = None;
    for l in body.lines().skip(1).map(str::trim) {
        if l.starts_with("_ =>") {
            if open.is_some() {
                return Err("isotope_abundances: the `_` arm inside an arm".to_string());
            }
            return Ok(out);
        }
        if open.is_none() {
            let (z, _) = l.split_once(" =>").ok_or_else(|| format!("isotope_abundances: {l:?} is no `<z> =>` arm"))?;
            open = Some(z.parse::<u8>().map_err(|e| format!("isotope_abundances: {z:?}: {e}"))?);
        }
        let code = l.split("//").next().unwrap_or_default();
        if code.trim_end().ends_with("],") {
            let note = l.split_once("//").map_or(String::new(), |(_, n)| n.trim().to_string());
            out.push((open.take().unwrap_or_default(), note));
        }
    }
    Err("isotope_abundances: the match has no `_` arm".to_string())
}

/// A function's or method's Rust doc as the plain lines of a column doc.
fn plain_doc(src: &str, item: &str) -> Result<String, String> {
    let mut out = String::new();
    for l in doc_comment(src, item)?.lines() {
        let _ = writeln!(out, "{}", l.trim_start_matches('#').trim_start());
    }
    Ok(out)
}

#[allow(clippy::too_many_lines)] // the enum, then the table, then the constants
pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let names: Vec<String> = SERIES.iter().map(|s| format!("{s:?}")).collect();
    let mut out = header(
        "spectroscopy data: the spectral series, the isotope abundances, R_H",
        provenance,
        &[
            "KmSpectralSeries follows the Rust declaration order and is append-only: value k is the k-th variant,\n\
             and KM_SPECTRAL_SERIES_N_LOWER is indexed by the value. The isotope abundances are\n\
             isotope_abundances(z) for every u8 z, the oracle's values as round-trip decimals; tests/golden/\n\
             spectroscopy.tsv pins every cell (the `spectroscopy::isotope_abundances` rows, every u8) and every\n\
             series. PLANCK, SPEED_OF_LIGHT, RYDBERG and RYDBERG_EV are src/constants.cyr's.",
        ],
    );

    // ── SpectralSeries ──
    let members: Vec<Member> = names
        .iter()
        .map(|v| Member { variant: v.clone(), name: format!("KM_SPECTRAL_SERIES_{}", snake(v)), extra: Vec::new() })
        .collect();
    out.push_str(&enum_decl(src, "SpectralSeries", "KmSpectralSeries", "KM_SPECTRAL_SERIES_COUNT", &members)?);
    let arms = match_arms(src, "pub const fn n_lower(")?;
    let mut armed: Vec<&str> = arms.iter().map(|a| a.variant.as_str()).collect();
    armed.sort_unstable();
    let mut want: Vec<&str> = names.iter().map(String::as_str).collect();
    want.sort_unstable();
    if armed != want || arms.iter().any(|a| !a.sections.is_empty() || a.note.is_some()) {
        return Err(format!("n_lower: the match arms {armed:?} are not the variants {want:?}, or carry comments"));
    }
    let n_lower = SERIES
        .iter()
        .map(|s| u8::try_from(s.n_lower()).map_err(|_| format!("{s:?}.n_lower() > 255")))
        .collect::<Result<Vec<u8>, String>>()?;
    out.push_str(&table(
        "KM_SPECTRAL_SERIES_N_LOWER",
        &format!(
            "{}Indexed by KmSpectralSeries value; mirrors `SpectralSeries::n_lower` (a u32 in Rust, at most 5).",
            plain_doc(src, "pub const fn n_lower(")?
        ),
        &Elem::U8(n_lower),
        &names,
    )?);

    // ── the isotope abundances ──
    let arms = abundance_arms(src)?;
    let mut with_data = Vec::new();
    for z in 0..=255u8 {
        if !isotope_abundances(z).is_empty() {
            with_data.push(z);
        }
    }
    let armed: Vec<u8> = arms.iter().map(|(z, _)| *z).collect();
    if armed != with_data {
        return Err(format!("isotope_abundances: the match arms {armed:?} are not the API's non-empty z {with_data:?}"));
    }
    let zs = usize::from(*with_data.last().ok_or("isotope_abundances: no element has data")?) + 1;
    let mut pairs: Vec<f64> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    let mut start = vec![0u8; zs];
    let mut count = vec![0u8; zs];
    let mut zlabels: Vec<String> = Vec::new();
    for z in 0..zs {
        let zu = u8::try_from(z).map_err(|_| "z > 255".to_string())?;
        let a = isotope_abundances(zu);
        let note = arms.iter().find(|(az, _)| *az == zu).map_or("", |(_, n)| n.as_str());
        start[z] = u8::try_from(pairs.len() / 2).map_err(|_| "more than 255 isotope pairs".to_string())?;
        count[z] = u8::try_from(a.len()).map_err(|_| format!("z = {z}: more than 255 isotopes"))?;
        zlabels.push(if a.is_empty() { format!("Z {z}: none") } else { format!("Z {z}: {note}") });
        for (k, &(mass, ab)) in a.iter().enumerate() {
            pairs.push(mass);
            pairs.push(ab);
            labels.push(format!("Z {z} isotope {k} mass (u){}", if k == 0 { format!(": {note}") } else { String::new() }));
            labels.push(format!("Z {z} isotope {k} abundance"));
        }
    }
    out.push_str(&plain_doc(src, "pub fn isotope_abundances(")?.lines().map(comment_line).collect::<Vec<_>>().join("\n"));
    out.push_str("\n#\n# The pairs of every z with data, concatenated in z order: z's slice is the _KM_ISOTOPE_ABUNDANCE_COUNT[z]\n");
    out.push_str("# pairs from pair _KM_ISOTOPE_ABUNDANCE_START[z] on (16 bytes each, mass then abundance); src/spectroscopy.cyr's\n");
    out.push_str("# isotope_abundances(z) returns (pointer, count). Every z from _KM_ISOTOPE_ABUNDANCE_ZS to 255 has none.\n");
    let _ = writeln!(out, "const _KM_ISOTOPE_ABUNDANCE_ZS = {zs};\n");
    out.push_str(&table(
        "_KM_ISOTOPE_ABUNDANCES",
        "The (isotope mass in u, abundance as a mole fraction) pairs, interleaved; mirrors `isotope_abundances(z)`.",
        &Elem::F64(pairs),
        &labels,
    )?);
    out.push_str(&table(
        "_KM_ISOTOPE_ABUNDANCE_START",
        "The first pair of z's slice in _KM_ISOTOPE_ABUNDANCES (pair index), indexed by z.",
        &Elem::U8(start),
        &zlabels,
    )?);
    out.push_str(&table(
        "_KM_ISOTOPE_ABUNDANCE_COUNT",
        "The length of z's slice (0: no isotope data), indexed by z.",
        &Elem::U8(count),
        &zlabels,
    )?);

    // ── the private constants ──
    let ratio = private_literal(src, "ELECTRON_PROTON_MASS_RATIO")?;
    let rh = rydberg_hydrogen(ratio)?;
    out.push_str(&doc_comment(src, "const ELECTRON_PROTON_MASS_RATIO:")?);
    out.push_str("# Rust's private constant, from its literal: R_H's cross-check (tests/spectroscopy.tcyr evaluates\n");
    out.push_str("# KM_RYDBERG / (1.0 + this) and requires _KM_RYDBERG_HYDROGEN).\n");
    out.push_str("# Rust: kimiya::spectroscopy::ELECTRON_PROTON_MASS_RATIO (private)\n");
    let _ = writeln!(out, "const _KM_ELECTRON_PROTON_MASS_RATIO = {};\n", f64_lit(ratio)?);
    out.push_str(&doc_comment(src, "const RYDBERG_HYDROGEN:")?);
    out.push_str("# Rust's private constant, RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO) as its const evaluation gives\n");
    out.push_str("# it, read through the API: SpectralSeries::Lyman.wavelength(u32::MAX − 1) is RN(1/R_H), and this is\n");
    out.push_str("# the one double whose reciprocal rounds to it (golden row `SpectralSeries::wavelength` Lyman 4294967294).\n");
    out.push_str("# Rust: kimiya::spectroscopy::RYDBERG_HYDROGEN (private)\n");
    let _ = writeln!(out, "const _KM_RYDBERG_HYDROGEN = {};", f64_lit(rh)?);
    Ok(out.trim_end().to_string() + "\n")
}
