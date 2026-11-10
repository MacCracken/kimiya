//! `src/spectroscopy_data.cyr` — the spectroscopy data of `rust-old/src/spectroscopy.rs`: the `SpectralSeries` enum and
//! its `n_lower` column, the isotope abundance table, and Rust's private Rydberg constant for ¹H.
//!
//! - `enum KmSpectralSeries` (`KM_SPECTRAL_SERIES_COUNT`): the Rust variants in declaration order, value k for the k-th
//!   (append-only), proved against the source declaration (`cyr::enum_decl`; `#[non_exhaustive]` keeps
//!   `m_spectroscopy::SERIES` from an exhaustive match), and `KM_SPECTRAL_SERIES_N_LOWER` (u8), `SpectralSeries::n_lower`
//!   through the API (its match arms checked to name every variant once).
//! - The isotope abundances (ADR 0002 #23, M4e-4): **not the oracle's values** but the cited source table
//!   `scripts/oracle/isotopes.tsv` (AME2020 masses at full published precision, the TICE 2009 representative abundances,
//!   ³⁶S included), read by [`isotope_table`] and checked against the oracle before anything is written: the same elements
//!   and isotopes as `isotope_abundances(z)` at **every u8 z** plus ³⁶S as the only addition, every shared abundance Rust's
//!   literal as a decimal string, every Rust mass the AME2020 mass rounded to five decimals, each element's abundances
//!   summing to exactly 1 (exact decimal arithmetic). The (mass, abundance) pairs of every z with data, concatenated, are
//!   `_KM_ISOTOPE_ABUNDANCES` (f64, interleaved: each the double nearest its decimal), with `_KM_ISOTOPE_ABUNDANCE_START` /
//!   `_COUNT` (u8, pair index and count) indexed by z for z below `_KM_ISOTOPE_ABUNDANCE_ZS` (one past the last z with
//!   data; every z from there to 255 has none, as the oracle's).
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
//! Values come from the oracle's API, but for the isotope table (ADR 0002 #23: from its cited source file, checked against
//! the oracle's API and the literal); the Rust source is read for comments, the private ratio's literal (a cross-check of
//! R_H) and the isotope literal's decimal strings.

use crate::cyr::{Elem, Member, doc_comment, enum_decl, f64_lit, header, match_arms, snake, table};
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

/// The source file of ADR 0002 #23's isotope table, from the repository root.
pub const ISOTOPES_TSV: &str = "scripts/oracle/isotopes.tsv";

/// One `isotope` record of [`ISOTOPES_TSV`] (ADR 0002 #23, M4e-4).
pub struct Iso {
    pub z: u8,
    /// The mass number.
    pub a: u16,
    pub symbol: String,
    /// The AME2020 atomic mass in u, its published decimal digits.
    pub mass: String,
    /// The TICE 2009 representative abundance (mole fraction), its printed digits.
    pub abundance: String,
}

impl Iso {
    /// `³⁶S`: the mass number in superscript digits and the symbol.
    pub fn label(&self) -> String {
        const SUP: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
        let a: String = self.a.to_string().bytes().map(|d| SUP[usize::from(d - b'0')]).collect();
        format!("{a}{}", self.symbol)
    }
}

/// A decimal `d…d.d…d` (no sign, exponent or spaces) as (digits, scale): the value is digits / 10^scale.
fn decimal(s: &str) -> Result<(u128, u32), String> {
    let (i, f) = s.split_once('.').ok_or_else(|| format!("{s:?}: no decimal point"))?;
    if i.is_empty() || f.is_empty() || !i.bytes().chain(f.bytes()).all(|c| c.is_ascii_digit()) || i.len() + f.len() > 30 {
        return Err(format!("{s:?}: not a decimal d.d"));
    }
    let digits = format!("{i}{f}").parse::<u128>().map_err(|e| format!("{s:?}: {e}"))?;
    Ok((digits, u32::try_from(f.len()).map_err(|_| format!("{s:?}: too long"))?))
}

/// x (digits, scale) at `scale` digits, rounded to nearest; an exact tie is refused (it would make "rounded" ambiguous).
fn decimal_at(x: (u128, u32), scale: u32) -> Result<u128, String> {
    let (d, s) = x;
    if s <= scale {
        return Ok(d * 10u128.pow(scale - s));
    }
    let unit = 10u128.pow(s - scale);
    let (q, r) = (d / unit, d % unit);
    match (2 * r).cmp(&unit) {
        std::cmp::Ordering::Less => Ok(q),
        std::cmp::Ordering::Greater => Ok(q + 1),
        std::cmp::Ordering::Equal => Err(format!("{d}e-{s} is a tie at {scale} decimals")),
    }
}

/// Each z's (mass, abundance) literal strings, by z.
type Literals = Vec<(u8, Vec<(String, String)>)>;

/// The (mass, abundance) literals of `isotope_abundances`' match arms, by z in source order, underscores removed (Rust's
/// `0.999_885` is `0.999885`); comments are not read.
fn rust_literals(src: &str) -> Result<Literals, String> {
    let at = src.find("pub fn isotope_abundances(").ok_or("no `pub fn isotope_abundances(`")?;
    let body = &src[at..];
    let body = &body[body.find("match atomic_number {").ok_or("isotope_abundances: no `match atomic_number {`")?..];
    let mut code = String::new();
    let mut closed = false;
    for l in body.lines().skip(1) {
        let c = l.split("//").next().unwrap_or_default().trim();
        if c.starts_with("_ =>") {
            closed = true;
            break;
        }
        code.push_str(c);
        code.push(' ');
    }
    if !closed {
        return Err("isotope_abundances: the match has no `_` arm".to_string());
    }
    let mut out = Vec::new();
    let mut rest = code.as_str();
    while let Some(arrow) = rest.find(" => &[") {
        let z = rest[..arrow].trim().parse::<u8>().map_err(|e| format!("isotope_abundances: arm {:?}: {e}", &rest[..arrow]))?;
        let close = arrow + rest[arrow..].find(']').ok_or("isotope_abundances: an arm without `]`")?;
        let mut pairs = Vec::new();
        for part in rest[arrow + 6..close].split('(').skip(1) {
            let tuple = part.split(')').next().unwrap_or_default();
            let (m, a) = tuple.split_once(',').ok_or_else(|| format!("isotope_abundances: z = {z}: {tuple:?}"))?;
            pairs.push((m.trim().replace('_', ""), a.trim().replace('_', "")));
        }
        out.push((z, pairs));
        rest = rest[close + 1..].trim_start().trim_start_matches(',');
    }
    if !rest.trim().is_empty() {
        return Err(format!("isotope_abundances: unparsed {:?}", rest.trim()));
    }
    Ok(out)
}

/// One `isotope` line of [`ISOTOPES_TSV`], its mass checked against the verbatim AME2020 line and its abundance against
/// the printed TICE 2009 value.
fn iso_record(f: &[&str]) -> Result<Iso, String> {
    let [_, z, a, sym, mass, ab, tice, line] = f else {
        return Err(format!("{ISOTOPES_TSV}: an isotope record of {} fields, not 8", f.len()));
    };
    let iso = Iso {
        z: z.parse().map_err(|e| format!("{ISOTOPES_TSV}: Z {z:?}: {e}"))?,
        a: a.parse().map_err(|e| format!("{ISOTOPES_TSV}: A {a:?}: {e}"))?,
        symbol: (*sym).to_string(),
        mass: (*mass).to_string(),
        abundance: (*ab).to_string(),
    };
    let who = format!("{ISOTOPES_TSV}: {}", iso.label());
    // mass_1.mas20.txt's format (its header): a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,
    // i3,1x,f13.6,f12.6; Z at columns 10–14, A 15–19, the element 21–23, the atomic mass's u 107–109 and micro-u 111–123.
    if !line.is_ascii() || line.len() < 123 {
        return Err(format!("{who}: the mas20 line is not an AME2020 mass line"));
    }
    let col = |a: usize, b: usize| line[a - 1..b].trim();
    if col(10, 14) != *z || col(15, 19) != *a || col(21, 23) != *sym {
        return Err(format!("{who}: the mas20 line is Z {}, A {}, {}", col(10, 14), col(15, 19), col(21, 23)));
    }
    let (u, micro) = (col(107, 109), col(111, 123));
    let (mi, mf) = micro.split_once('.').ok_or_else(|| format!("{who}: micro-u {micro:?}"))?;
    if u.parse::<u32>().is_err() || mi.len() != 6 || format!("{}.{mi}{mf}", u.parse::<u32>().unwrap_or(0)) != *mass {
        return Err(format!("{who}: mass_u {mass} is not the mas20 line's {u} {micro} micro-u"));
    }
    let printed = tice.split('(').next().unwrap_or_default().replace(' ', "");
    if printed != *ab || !tice.ends_with(')') {
        return Err(format!("{who}: abundance {ab} is not the printed {tice:?}"));
    }
    decimal(mass)?;
    decimal(ab)?;
    Ok(iso)
}

/// ADR 0002 #23's isotope table: the `isotope` records of [`ISOTOPES_TSV`] (each mass checked against its verbatim AME2020
/// line, each abundance against its printed TICE 2009 value), in Z then A order. It fails unless
/// (a) its elements and isotopes are the oracle's — `isotope_abundances(z)` at every u8 z, an isotope's A its mass rounded
///     (Rust's literal is that, checked against the API's bits) — plus ³⁶S, the only addition;
/// (b) every shared abundance is Rust's literal as a decimal string;
/// (c) every Rust mass is the AME2020 mass rounded to five decimals (exact decimal arithmetic; a tie is refused);
/// (d) each element's abundances sum to exactly 1 in decimal arithmetic.
pub fn isotope_table(root: &std::path::Path, src: &str) -> Result<Vec<Iso>, String> {
    let path = root.join(ISOTOPES_TSV);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut table: Vec<Iso> = Vec::new();
    for l in text.lines().filter(|l| l.starts_with("isotope\t")) {
        let iso = iso_record(&l.split('\t').collect::<Vec<_>>())?;
        if let Some(p) = table.last()
            && (p.z, p.a) >= (iso.z, iso.a)
        {
            return Err(format!("{ISOTOPES_TSV}: {} after {}: not in Z then A order", iso.label(), p.label()));
        }
        table.push(iso);
    }
    // (a) The oracle's isotopes, each by its mass number, against the source's without ³⁶S.
    let lits = rust_literals(src)?;
    let mut oracle: Vec<(u8, u16, String, String)> = Vec::new();
    for z in 0..=255u8 {
        let api = isotope_abundances(z);
        let lit = lits.iter().find(|(lz, _)| *lz == z).map_or(&[][..], |(_, p)| p.as_slice());
        if api.len() != lit.len() {
            return Err(format!("isotope_abundances({z}): {} pairs, the literal {}", api.len(), lit.len()));
        }
        for (&(m, ab), (lm, la)) in api.iter().zip(lit) {
            let (pm, pa) = (lm.parse::<f64>(), la.parse::<f64>());
            if pm.map(f64::to_bits) != Ok(m.to_bits()) || pa.map(f64::to_bits) != Ok(ab.to_bits()) {
                return Err(format!("isotope_abundances({z}): ({m:?}, {ab:?}) is not the literal ({lm}, {la})"));
            }
            let a = m.round();
            if !(1.0..=300.0).contains(&a) || (m - a).abs() > 0.1 {
                return Err(format!("isotope_abundances({z}): mass {m:?} names no mass number"));
            }
            oracle.push((z, u16::try_from(a as i64).map_err(|e| e.to_string())?, lm.clone(), la.clone()));
        }
    }
    let shared: Vec<&Iso> = table.iter().filter(|i| (i.z, i.a) != (16, 36)).collect();
    let file_za: Vec<(u8, u16)> = shared.iter().map(|i| (i.z, i.a)).collect();
    let oracle_za: Vec<(u8, u16)> = oracle.iter().map(|o| (o.0, o.1)).collect();
    if file_za != oracle_za || table.len() != shared.len() + 1 {
        return Err(format!("(a) {ISOTOPES_TSV}'s isotopes less ³⁶S {file_za:?} are not the oracle's {oracle_za:?} (or no ³⁶S)"));
    }
    for (iso, (_, _, lm, la)) in shared.iter().zip(&oracle) {
        // (b) The abundance as Rust wrote it.
        if iso.abundance != *la {
            return Err(format!("(b) {}: abundance {} is not Rust's literal {la}", iso.label(), iso.abundance));
        }
        // (c) Rust's mass, the AME2020 one to five decimals.
        if decimal_at(decimal(&iso.mass)?, 5)? != decimal_at(decimal(lm)?, 5)? || decimal(lm)?.1 > 5 {
            return Err(format!("(c) {}: Rust's mass {lm} is not AME2020's {} to five decimals", iso.label(), iso.mass));
        }
    }
    // (d) Each element's abundances sum to 1, exactly.
    for z in table.iter().map(|i| i.z).collect::<std::collections::BTreeSet<_>>() {
        let abs = table.iter().filter(|i| i.z == z).map(|i| decimal(&i.abundance)).collect::<Result<Vec<_>, _>>()?;
        let scale = abs.iter().map(|d| d.1).max().unwrap_or(0);
        let mut sum = 0u128;
        for d in abs {
            sum += decimal_at(d, scale)?; // exact: no abundance has more digits than `scale`
        }
        if sum != 10u128.pow(scale) {
            return Err(format!("(d) Z = {z}: the abundances sum to {sum}e-{scale}, not 1"));
        }
    }
    Ok(table)
}

/// z's (mass, abundance) pairs in ADR 0002 #23's table, each the double nearest its decimal (Rust's `parse` is correctly
/// rounded); empty for a z without data.
pub fn new_abundances(table: &[Iso], z: u8) -> Vec<(f64, f64)> {
    table
        .iter()
        .filter(|i| i.z == z)
        .map(|i| (i.mass.parse::<f64>().unwrap_or(f64::NAN), i.abundance.parse::<f64>().unwrap_or(f64::NAN)))
        .collect()
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
pub fn emit(src: &str, provenance: &str, isos: &[Iso]) -> Result<String, String> {
    let names: Vec<String> = SERIES.iter().map(|s| format!("{s:?}")).collect();
    let mut out = header(
        "spectroscopy data: the spectral series, the isotope abundances, R_H",
        provenance,
        &[
            "KmSpectralSeries follows the Rust declaration order and is append-only: value k is the k-th variant,\n\
             and KM_SPECTRAL_SERIES_N_LOWER is indexed by the value. The isotope abundances are ADR 0002 #23's\n\
             table, from scripts/oracle/isotopes.tsv (AME2020 masses, TICE 2009 abundances, ³⁶S included), not the\n\
             oracle's: checked against isotope_abundances(z) at every u8 z (the same elements and isotopes, ³⁶S the\n\
             only addition; the abundances Rust's literals; Rust's masses these to five decimals), written as the\n\
             round-trip decimals of the doubles nearest the published values. tests/golden/spectroscopy.tsv pins\n\
             every cell (the `spectroscopy::isotope_abundances` rows, every u8: the oracle's values, the port's in\n\
             their ADR 0002 #23 notes) and every series. PLANCK, SPEED_OF_LIGHT, RYDBERG and RYDBERG_EV are\n\
             src/constants.cyr's.",
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

    // ── the isotope abundances (ADR 0002 #23: the cited table, checked against the oracle by `isotope_table`) ──
    let zs = usize::from(isos.last().ok_or("isotope_abundances: no element has data")?.z) + 1;
    if (zs..=255).any(|z| u8::try_from(z).is_ok_and(|z| !isotope_abundances(z).is_empty())) {
        return Err(format!("isotope_abundances: the oracle has data at or past z = {zs}"));
    }
    let mut pairs: Vec<f64> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    let mut start = vec![0u8; zs];
    let mut count = vec![0u8; zs];
    let mut zlabels: Vec<String> = Vec::new();
    for z in 0..zs {
        let zu = u8::try_from(z).map_err(|_| "z > 255".to_string())?;
        let here: Vec<&Iso> = isos.iter().filter(|i| i.z == zu).collect();
        start[z] = u8::try_from(pairs.len() / 2).map_err(|_| "more than 255 isotope pairs".to_string())?;
        count[z] = u8::try_from(here.len()).map_err(|_| format!("z = {z}: more than 255 isotopes"))?;
        let names: Vec<String> = here.iter().map(|i| i.label()).collect();
        zlabels.push(if here.is_empty() { format!("Z {z}: none") } else { format!("Z {z}: {}", names.join(", ")) });
        for (iso, (mass, ab)) in here.iter().zip(new_abundances(isos, zu)) {
            pairs.push(mass);
            pairs.push(ab);
            labels.push(format!("{} mass (u), AME2020: {}", iso.label(), iso.mass));
            labels.push(format!("{} abundance, TICE 2009: {}", iso.label(), iso.abundance));
        }
    }
    out.push_str(
        "# Natural abundance data for common isotopes (ADR 0002 #23, M4e-4): (isotope mass in u, abundance as a mole\n\
         # fraction) pairs for H, C, N, O, S, Cl and Br; every other element has none. The masses are AME2020's at full\n\
         # published precision (M. Wang, W.J. Huang, F.G. Kondev, G. Audi and S. Naimi, Chinese Physics C 45, 030003\n\
         # (2021): mass_1.mas20.txt), the abundances the IUPAC 2009 representative isotopic compositions (M. Berglund and\n\
         # M.E. Wieser, Pure Appl. Chem. 83, 397-410 (2011): TICE 2009, Table 1, column 9), ³⁶S included; each value is the\n\
         # double nearest its published decimal (the label). The records, with the sources' URLs and sha256s and each\n\
         # mass's verbatim AME2020 line, are scripts/oracle/isotopes.tsv. Rust 1.2.0's table rounds the masses to five\n\
         # decimals and omits ³⁶S (its sulfur sums to 0.9999): its values are the golden rows', these their ADR 0002 #23\n\
         # notes.\n\
         #\n\
         # The pairs of every z with data, concatenated in z order: z's slice is the _KM_ISOTOPE_ABUNDANCE_COUNT[z]\n\
         # pairs from pair _KM_ISOTOPE_ABUNDANCE_START[z] on (16 bytes each, mass then abundance); src/spectroscopy.cyr's\n\
         # isotope_abundances(z) returns (pointer, count). Every z from _KM_ISOTOPE_ABUNDANCE_ZS to 255 has none.\n",
    );
    let _ = writeln!(out, "const _KM_ISOTOPE_ABUNDANCE_ZS = {zs};\n");
    out.push_str(&table(
        "_KM_ISOTOPE_ABUNDANCES",
        "The (isotope mass in u, abundance as a mole fraction) pairs, interleaved; mirrors `isotope_abundances(z)` (ADR\n\
         0002 #23's values).",
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
