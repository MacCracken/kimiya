//! `src/gas_data.cyr` — the gas data of `rust-old/src/gas.rs`: the `PengRobinsonParams` struct,
//! the five `PR_*` presets, and the STP constants.
//!
//! - `struct KmPengRobinsonParams`: `PengRobinsonParams`'s fields (tc, pc, omega) in declaration
//!   order, all f64. A by-value struct, as tanmatra's `TmFourMomentum` (port plan §6.4).
//! - `PR_*`: each preset is a top-level struct value, `var KM_PR_<NAME> = KmPengRobinsonParams {
//!   tc: …, pc: …, omega: … };` (probed on 6.7.18: a top-level struct literal of f64 fields is
//!   initialised data, and passes to a `p: KmPengRobinsonParams` parameter). `emit` proves the list
//!   is every `pub const PR_*: PengRobinsonParams` in the source. tests/golden/gas.tsv pins each
//!   field (`gas::PR_<NAME>` rows, `m_gas.rs`).
//! - `STP_TEMPERATURE`, `STP_PRESSURE`: `t_constants::consts` for `Home::Data("gas_data")`.
//!
//! Values come from the oracle's API; the Rust source is read for comments only.

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::t_constants::{Home, consts};
use kimiya::gas::{PR_CO2, PR_METHANE, PR_NITROGEN, PR_OXYGEN, PR_WATER, PengRobinsonParams};
use std::fmt::Write as _;
use std::path::Path;

/// Every `PR_*` preset, in source order.
pub const PRESETS: [(&str, PengRobinsonParams); 5] = [
    ("PR_WATER", PR_WATER),
    ("PR_NITROGEN", PR_NITROGEN),
    ("PR_OXYGEN", PR_OXYGEN),
    ("PR_CO2", PR_CO2),
    ("PR_METHANE", PR_METHANE),
];

/// The `PengRobinsonParams` fields, in declaration order.
const FIELDS: [&str; 3] = ["tc", "pc", "omega"];

/// The `pub const PR_*: PengRobinsonParams` names in the source, in order.
fn source_presets(src: &str) -> Vec<String> {
    src.lines()
        .filter_map(|l| l.strip_prefix("pub const "))
        .filter_map(|l| l.split_once(": PengRobinsonParams"))
        .map(|(n, _)| n.to_string())
        .collect()
}

pub fn emit(root: &Path, src: &str, provenance: &str) -> Result<String, String> {
    let listed: Vec<String> = PRESETS.iter().map(|(n, _)| (*n).to_string()).collect();
    let in_source = source_presets(src);
    if listed != in_source {
        return Err(format!(
            "t_gas::PRESETS {listed:?} is not the source's PengRobinsonParams consts {in_source:?}"
        ));
    }
    let mut out = header(
        "gas data: STP and the Peng–Robinson parameters of five gases",
        provenance,
        &[
            "f64 values are the oracle's as round-trip decimals; tests/golden/constants.tsv pins the\n\
             STP constants and tests/golden/gas.tsv every preset field (`gas::PR_*` rows).",
        ],
    );

    out.push_str(&consts(root, Home::Data("gas_data"))?);

    out.push_str(&doc_comment(src, "pub struct PengRobinsonParams {")?);
    out.push_str("struct KmPengRobinsonParams {\n");
    for field in FIELDS {
        let docs =
            doc_comment_after(src, "pub struct PengRobinsonParams {", &format!("pub {field}:"))?;
        for l in docs.lines() {
            let _ = writeln!(out, "    {l}");
        }
        let _ = writeln!(out, "    {field}: f64;");
    }
    out.push_str("}\n\n");

    for (name, p) in &PRESETS {
        out.push_str(&doc_comment(src, &format!("pub const {name}:"))?);
        let _ = writeln!(out, "# Rust: kimiya::gas::{name}");
        let _ = writeln!(
            out,
            "var KM_{name} = KmPengRobinsonParams {{ tc: {}, pc: {}, omega: {} }};\n",
            f64_lit(p.tc)?,
            f64_lit(p.pc)?,
            f64_lit(p.omega)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}
