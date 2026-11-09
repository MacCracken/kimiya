//! `molecule.tsv` — every `Molecule` function over the four presets and seeded random
//! compositions (port plan §4, M1 grid: molecular weight and formula for the presets and seeded
//! random compositions).
//!
//! A molecule is written as its atoms, `v:n` then `i:Z i:count` per atom in stored order.
//! - `Molecule::water` / `carbon_dioxide` / `methane` / `glucose`: no input, the atoms.
//! - `Molecule::new`: the (Z, count) pairs in, the atoms out.
//! - `Molecule::molecular_weight` (class B: Σ mass·count in atom order, + and × only, expected
//!   0 ulp), `Molecule::formula` and `Molecule::total_atoms` (class A): the receiver's atoms in,
//!   `f:` / `s:` / `i:` out, or `E:InvalidElement`.
//!
//! The compositions: the cases the 1.2 tests pin (zero counts before and after a valid atom, an
//! unknown element with count 0 and 1, input order kept, u32::MAX saturation, the empty
//! molecule), then 600 seeded ones of 1–8 atoms with counts 0..=200, about 1 in 7 atoms
//! repeating an earlier element, 1 in 25 an invalid Z (0, 119 or 255) and 1 in 33 a count of
//! u32::MAX.

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use kimiya::molecule::Molecule;

fn atoms(m: &Molecule) -> Vec<Tok> {
    let mut v = vec![Tok::V(m.atoms.len())];
    for a in &m.atoms {
        v.push(i(a.element_number));
        v.push(i(a.count));
    }
    v
}

fn pairs_toks(pairs: &[(u8, u32)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(pairs.len())];
    for &(z, c) in pairs {
        v.push(i(z));
        v.push(i(c));
    }
    v
}

/// The three methods on one receiver.
fn methods(md: &mut Module, mol: &Molecule) {
    let recv = atoms(mol);
    let mw = mol.molecular_weight().map_or_else(|x| err(&x), f);
    md.row("Molecule::molecular_weight", Class::B, &recv, &[mw]);
    let fo = mol.formula().map_or_else(|x| err(&x), s);
    md.row("Molecule::formula", Class::A, &recv, &[fo]);
    md.row(
        "Molecule::total_atoms",
        Class::A,
        &recv,
        &[i(mol.total_atoms())],
    );
}

fn small(rng: &mut Rng, lo: i64, hi: i64) -> usize {
    usize::try_from(rng.int(lo, hi)).unwrap_or(0)
}

fn compositions() -> Vec<Vec<(u8, u32)>> {
    let mut out: Vec<Vec<(u8, u32)>> = vec![
        vec![],
        vec![(1, 0), (8, 1)],
        vec![(6, 1), (1, 0)],
        vec![(1, 0)],
        vec![(8, 1), (1, 2)],
        vec![(200, 0)],
        vec![(200, 1)],
        vec![(6, 1), (200, 0)],
        vec![(1, u32::MAX), (8, 2)],
        vec![(1, u32::MAX), (1, u32::MAX)],
        vec![(118, u32::MAX)],
        vec![(1, 1), (1, 1), (8, 1)],
        vec![(6, 6), (1, 12), (8, 6)],
        vec![(255, 0), (0, 0)],
    ];
    let mut rng = Rng::new(0x3a01);
    for _ in 0..600 {
        let n = small(&mut rng, 1, 8);
        let mut mol: Vec<(u8, u32)> = Vec::with_capacity(n);
        for _ in 0..n {
            let r = rng.int(0, 99);
            let z = if r < 4 {
                [0_u8, 119, 255][small(&mut rng, 0, 2)]
            } else if r < 18 && !mol.is_empty() {
                mol[small(&mut rng, 0, i64::try_from(mol.len()).unwrap_or(1) - 1)].0
            } else {
                u8::try_from(rng.int(1, 118)).unwrap_or(1)
            };
            let c = if rng.int(0, 99) < 3 {
                u32::MAX
            } else {
                u32::try_from(rng.int(0, 200)).unwrap_or(0)
            };
            mol.push((z, c));
        }
        out.push(mol);
    }
    out
}

pub fn emit() -> Module {
    let mut m = Module::new("molecule");
    let presets: [(&str, Molecule); 4] = [
        ("Molecule::water", Molecule::water()),
        ("Molecule::carbon_dioxide", Molecule::carbon_dioxide()),
        ("Molecule::methane", Molecule::methane()),
        ("Molecule::glucose", Molecule::glucose()),
    ];
    for (name, mol) in &presets {
        m.row(name, Class::A, &[], &atoms(mol));
        methods(&mut m, mol);
    }
    for pairs in compositions() {
        let mol = Molecule::new(&pairs);
        m.row("Molecule::new", Class::A, &pairs_toks(&pairs), &atoms(&mol));
        methods(&mut m, &mol);
    }
    m
}
