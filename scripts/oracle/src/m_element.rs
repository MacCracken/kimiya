//! `element.tsv` — every element and every extended-property row, through the oracle's three
//! lookups (port plan §4, M1 grid: all 118 elements × every field; lookups incl. invalid symbol
//! and Z). Class A throughout: the rows are table cells.
//!
//! - `element::lookup_by_number`: every u8, 0 ..= 255 (Rust's whole input domain).
//! - `element::lookup_by_symbol`: every symbol; its lower- and upper-case forms and the
//!   element's name (all `None`: the match is exact and case-sensitive); and the invalid inputs
//!   `""`, `"Xx"`, `"FE"`, `"Fe "`, `" Fe"` and a 300-byte string.
//! - `element::lookup_properties`: every u8, 0 ..= 255.
//!
//! A found `Element` is flattened in declaration order — `i:Z s:symbol s:name f:mass` then the
//! electronegativity (`f:` or `n:` for `None`) and `e:ElementCategory.<Variant>`. An
//! `ElementProperties` is `i:Z f:IE₁ f:EA f:radius v:n` and its n oxidation states as `i:`.
//! `None` is `n:`.

use std::collections::BTreeSet;

use crate::tsv::{Class, Module, Tok, e, f, i, s};
use kimiya::element::{
    ELEMENTS, Element, ElementCategory, ElementProperties, lookup_by_number, lookup_by_symbol,
    lookup_properties,
};

/// Every `ElementCategory`, in Rust declaration order. `t_element` proves the list against
/// the enum in `rust-old/src/element.rs` (the enum is `#[non_exhaustive]`, so no match here
/// could catch a missing variant).
pub const CATEGORIES: [ElementCategory; 10] = [
    ElementCategory::AlkaliMetal,
    ElementCategory::AlkalineEarth,
    ElementCategory::TransitionMetal,
    ElementCategory::PostTransitionMetal,
    ElementCategory::Metalloid,
    ElementCategory::Nonmetal,
    ElementCategory::Halogen,
    ElementCategory::NobleGas,
    ElementCategory::Lanthanide,
    ElementCategory::Actinide,
];

fn element_toks(el: Option<&Element>) -> Vec<Tok> {
    match el {
        None => vec![Tok::None],
        Some(el) => vec![
            i(el.atomic_number),
            s(el.symbol),
            s(el.name),
            f(el.atomic_mass),
            el.electronegativity.map_or(Tok::None, f),
            e("ElementCategory", el.category),
        ],
    }
}

fn properties_toks(p: Option<&ElementProperties>) -> Vec<Tok> {
    match p {
        None => vec![Tok::None],
        Some(p) => {
            let mut v = vec![
                i(p.atomic_number),
                f(p.ionization_energy_1),
                f(p.electron_affinity),
                f(p.covalent_radius_pm),
                Tok::V(p.oxidation_states.len()),
            ];
            v.extend(p.oxidation_states.iter().map(|&x| i(x)));
            v
        }
    }
}

/// The symbol probes, in emission order, each once.
fn symbol_probes() -> Vec<String> {
    let mut probes: Vec<String> = ELEMENTS.iter().map(|el| el.symbol.to_string()).collect();
    for el in ELEMENTS {
        probes.push(el.symbol.to_lowercase());
        probes.push(el.symbol.to_uppercase());
        probes.push(el.name.to_string());
    }
    probes.extend(
        ["", "Xx", "FE", "Fe ", " Fe"]
            .iter()
            .map(|x| (*x).to_string()),
    );
    probes.push("Fe".repeat(150));
    let mut seen = BTreeSet::new();
    probes.retain(|p| seen.insert(p.clone()));
    probes
}

pub fn emit() -> Module {
    let mut m = Module::new("element");
    for z in 0..=u8::MAX {
        m.row(
            "element::lookup_by_number",
            Class::A,
            &[i(z)],
            &element_toks(lookup_by_number(z)),
        );
    }
    for sym in symbol_probes() {
        let out = element_toks(lookup_by_symbol(&sym));
        m.row("element::lookup_by_symbol", Class::A, &[s(sym)], &out);
    }
    for z in 0..=u8::MAX {
        m.row(
            "element::lookup_properties",
            Class::A,
            &[i(z)],
            &properties_toks(lookup_properties(z)),
        );
    }
    m
}
