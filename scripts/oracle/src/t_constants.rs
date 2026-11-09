//! `src/constants.cyr` — kimiya's physical and universal constants, as `KM_` consts.
//!
//! kimiya 1.2.0 has 45 public scalar `f64` consts spread over 11 modules. `ALL` lists every
//! one (`complete` proves it against `rust-old/src`), and `tests/golden/constants.tsv` pins
//! every bit (`m_constants.rs`).
//! - **Physical constants** (`Home::Constants`) are generated here, into `src/constants.cyr`,
//!   which comes first in `[lib]` because 14 modules reach across for them.
//! - **Module data** (`Home::Module`: K_b/K_f, Madelung, Born exponents, Λ°, Debye–Hückel,
//!   K_w, STP, c_water) stays with its module and is generated with that module's data at
//!   its milestone. Until then its golden rows are *pending* in the harness.
//!
//! Values come from the oracle's API; the comments are the Rust doc comments, which carry
//! the citations (CODATA 2018, SI 2019, AME2020, …). Struct-valued consts (`gas::PR_*`)
//! arrive with their module (M2).

use crate::cyr::{doc_comment, f64_lit, header};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// Where a constant is generated.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Home {
    /// `src/constants.cyr` (this generator).
    Constants,
    /// With its module's data, at the named milestone.
    Module(&'static str),
}

/// One public `f64` const of the oracle.
pub struct Constant {
    pub module: &'static str,
    pub name: &'static str,
    pub value: f64,
    pub home: Home,
}

impl Constant {
    /// `element::AVOGADRO`, the golden row's input.
    pub fn path(&self) -> String {
        format!("{}::{}", self.module, self.name)
    }
}

macro_rules! constants {
    ($($module:ident :: $name:ident => $home:expr),* $(,)?) => {
        &[ $( Constant {
            module: stringify!($module),
            name: stringify!($name),
            value: kimiya::$module::$name,
            home: $home,
        } ),* ]
    };
}

const PHYS: Home = Home::Constants;
const M2: Home = Home::Module("M2");
const M3: Home = Home::Module("M3");
const M4: Home = Home::Module("M4");

/// Every public scalar `f64` const of kimiya 1.2.0: modules in file-name order, then
/// source order.
const ALL: &[Constant] = constants!(
    electrochemistry::FARADAY => PHYS,
    electrochemistry::ELEMENTARY_CHARGE => PHYS,
    element::AVOGADRO => PHYS,
    element::GAS_CONSTANT => PHYS,
    gas::STP_TEMPERATURE => M2,
    gas::STP_PRESSURE => M2,
    inorganic::MADELUNG_NACL => M4,
    inorganic::MADELUNG_CSCL => M4,
    inorganic::MADELUNG_ZNS => M4,
    inorganic::MADELUNG_FLUORITE => M4,
    inorganic::MADELUNG_RUTILE => M4,
    inorganic::BORN_EXP_HE => M4,
    inorganic::BORN_EXP_NE => M4,
    inorganic::BORN_EXP_AR => M4,
    inorganic::BORN_EXP_KR => M4,
    inorganic::BORN_EXP_XE => M4,
    kinetics::BOLTZMANN => PHYS,
    nuclear::AMU_MEV => PHYS,
    nuclear::PROTON_MASS_U => PHYS,
    nuclear::NEUTRON_MASS_U => PHYS,
    nuclear::HYDROGEN_ATOM_MASS_U => PHYS,
    phase::KB_WATER => M2,
    phase::KB_BENZENE => M2,
    phase::KB_ETHANOL => M2,
    phase::KB_CHLOROFORM => M2,
    phase::KF_WATER => M2,
    phase::KF_BENZENE => M2,
    phase::KF_CYCLOHEXANE => M2,
    phase::KF_CAMPHOR => M2,
    potential::COULOMB_CONSTANT => PHYS,
    solution::KW_WATER_25C => M3,
    solution::DH_A_WATER_25C => M3,
    solution::DH_B_WATER_25C => M3,
    solution::LAMBDA_H_PLUS => M3,
    solution::LAMBDA_OH_MINUS => M3,
    solution::LAMBDA_NA_PLUS => M3,
    solution::LAMBDA_K_PLUS => M3,
    solution::LAMBDA_CL_MINUS => M3,
    solution::LAMBDA_NO3_MINUS => M3,
    solution::LAMBDA_SO4_2MINUS => M3,
    spectroscopy::PLANCK => PHYS,
    spectroscopy::SPEED_OF_LIGHT => PHYS,
    spectroscopy::RYDBERG => PHYS,
    spectroscopy::RYDBERG_EV => PHYS,
    thermo::WATER_SPECIFIC_HEAT => M2,
);

/// Every constant, for the golden emitter too.
pub fn all() -> &'static [Constant] {
    ALL
}

/// Proves `ALL` is every `pub const NAME: f64` item at the top level of `rust-old/src/*.rs`
/// (no more, no fewer), and that the Cyrius names are unique.
pub fn complete(root: &Path) -> Result<(), String> {
    let dir = root.join("rust-old/src");
    let mut found = BTreeSet::new();
    let entries = std::fs::read_dir(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().is_none_or(|x| x != "rs") {
            continue;
        }
        let module = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string();
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("pub const ")
                && let Some((name, ty)) = rest.split_once(':')
                && ty.trim_start().starts_with("f64 ")
            {
                found.insert(format!("{module}::{name}"));
            }
        }
    }
    let listed: BTreeSet<String> = ALL.iter().map(Constant::path).collect();
    if listed.len() != ALL.len() {
        return Err("t_constants: a constant is listed twice".to_string());
    }
    let missing: Vec<_> = found.difference(&listed).collect();
    let extra: Vec<_> = listed.difference(&found).collect();
    if !missing.is_empty() || !extra.is_empty() {
        return Err(format!(
            "t_constants: not listed {missing:?}; listed but not in rust-old/src {extra:?}"
        ));
    }
    let mut names = BTreeMap::new();
    for c in ALL {
        if let Some(prev) = names.insert(c.name, c.module) {
            return Err(format!(
                "KM_{} would name both {prev}::{} and {}",
                c.name,
                c.name,
                c.path()
            ));
        }
    }
    Ok(())
}

pub fn emit(root: &Path, provenance: &str) -> Result<String, String> {
    complete(root)?;
    let mut out = header(
        "physical and universal constants",
        provenance,
        &[
            "Each value is the shortest decimal that round-trips to the oracle's f64, written out\n\
             positionally: the Cyrius lexer has no exponent notation, and decimal literals are\n\
             correctly rounded (6.6.10+). tests/golden/constants.tsv pins every bit.",
            "Module data constants (K_b/K_f, Madelung, Born exponents, Λ°, Debye–Hückel, K_w, STP,\n\
             water's specific heat) are generated with their module at its milestone, not here.",
        ],
    );
    let mut sources = BTreeMap::new();
    for c in ALL.iter().filter(|c| c.home == Home::Constants) {
        if !sources.contains_key(c.module) {
            let p = root.join(format!("rust-old/src/{}.rs", c.module));
            let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
            sources.insert(c.module, text);
        }
        let src = &sources[c.module];
        out.push_str(&doc_comment(src, &format!("pub const {}:", c.name))?);
        out.push_str(&format!("# Rust: kimiya::{}\n", c.path()));
        out.push_str(&format!("const KM_{} = {};\n\n", c.name, f64_lit(c.value)?));
    }
    Ok(out.trim_end().to_string() + "\n")
}
