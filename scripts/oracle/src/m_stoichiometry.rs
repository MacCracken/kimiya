//! `stoichiometry.tsv` — equation balancing, oxidation states, the limiting reagent, the yields, percent composition and
//! empirical formulas (port plan §4, M4f).
//!
//! - Class A: `balance_equation` (input the matrix as Rust's `&[Vec<f64>]`, `v:m` then each row `v:len` and its `f:`
//!   entries; output `v:n` and n `i:` coefficients, or `E:`): the Rust tests' matrices, CH₄ in every row order and
//!   negated, CuSO₄ + Fe in every row order (port plan §4: both orders), KMnO₄ + HCl in every row order and with its
//!   columns permuted, a list of textbook equations (ionic ones with a charge row), every combustion
//!   C_aH_bO_c + O₂ → CO₂ + H₂O with a ≤ 10, chains of 30 and 70 species (the port's heap tables), seeded consistent
//!   matrices (rows orthogonal to a random positive x, so x is the answer), seeded random matrices (every error class),
//!   and every error class by construction: empty, one species, ragged rows, non-integer and non-finite entries,
//!   entries on Rust's 10⁻⁹ integer tolerance to the ulp at several magnitudes, the 2⁵³ bound and integers past it (and
//!   past i64), the zero coefficient, mixed signs, the trivial and ambiguous null spaces, a coefficient at and past
//!   u32::MAX, intermediate i64 overflow (also where a valid u32 solution exists), and the reduced row whose gcd is 2⁶³.
//!   **ADR 0002 #7**: the port takes an integer matrix, so a row whose matrix has no integer counterpart (an entry that
//!   is not within 10⁻⁹ of an integer, non-finite, an integer beyond i64, or ragged rows) keeps the oracle's
//!   `E:InvalidInput` and is noted `[ADR 0002 #7: no integer-matrix counterpart (…)]`; the harness counts it.
//!   `assign_oxidation_states` (input `v:n` and n `i:z i:count` pairs, then `i:charge`; output `v:n` and n `i:z i:os`
//!   pairs, or `E:`): the 1.2 tests' species (port plan §4's oxidation corpus), every element alone over counts and
//!   charges, every element's binary oxides, hydrides and fluorides, oxyanions, peroxides, superoxides and ozonides,
//!   repeated entries, every error path, the u32 / i32 extremes and the i32 range's edges, seeded species. **ADR 0002
//!   #7**: where the oracle returns Err "non-integer average oxidation state" the port returns the exact rational
//!   states; those rows are noted with the port's whole result `(v:n i:z i:p i:q …)` (a reduced fraction, q > 0), from
//!   `port_oxidation`, this module's exact replica, which must give the oracle's result on every other row.
//!   `limiting_reagent` (two `v:n` slices; output `i:index`): the edges in every slot, the Rust tests, exact ties, ratios
//!   that overflow or underflow, seeded lists. `empirical_formula` (input `v:n` and n `i:z f:percent` pairs; output
//!   `v:n` and n `i:z i:count` pairs): the Rust tests, every multiplier 1 … 6, real compounds' percentages (from the
//!   oracle's own `percent_composition`), ratios on both sides of the 0.1 tolerance, every Err path, seeded formulas.
//! - Class B (+ − × ÷ only, expected 0 ulp): `theoretical_yield_moles` (the slices, then `f:product`), `percent_yield`
//!   (`f:actual f:theoretical`): the edges, the Rust tests, overflow and underflow, seeded values; `percent_composition`
//!   (input the molecule's atoms, `v:n` and `i:z i:count` pairs; output `v:n` and n `i:z f:percent` pairs): the presets,
//!   every element alone, binary compounds, zero and u32::MAX counts, invalid elements, the empty molecule, seeded
//!   compositions.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::exact::{cr_mul_div, empirical_exact, ratio_lt};
use crate::grids::{Rng, pow2, specials};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, err, f, i, same_toks, tokens_note};
use kimiya::element::{ElementCategory, lookup_by_number};
use kimiya::error::KimiyaError;
use kimiya::molecule::Molecule;
use kimiya::stoichiometry::{
    assign_oxidation_states, balance_equation, empirical_formula, limiting_reagent, percent_composition,
    percent_yield, theoretical_yield_moles,
};

// ── balance_equation ────────────────────────────────────────────────────────────────

fn matrix_toks(mat: &[Vec<f64>]) -> Vec<Tok> {
    let mut v = vec![Tok::V(mat.len())];
    for row in mat {
        v.push(Tok::V(row.len()));
        v.extend(row.iter().map(|&x| f(x)));
    }
    v
}

/// Why the matrix has no counterpart in the port's integer API (ADR 0002 #7), or `None` when every entry is Rust's
/// rounded integer of an i64 and the rows have one length. Rust's own acceptance test, in its order: finite, then
/// within 10⁻⁹ of `round(v)` (halves away from zero); an integer beyond i64 (|r| ≥ 2⁶³) has no i64 to pass.
pub fn no_counterpart(mat: &[Vec<f64>]) -> Option<&'static str> {
    if let Some(first) = mat.first()
        && mat.iter().any(|r| r.len() != first.len())
    {
        return Some("ragged rows");
    }
    for &v in mat.iter().flatten() {
        if !v.is_finite() {
            return Some("a non-finite entry");
        }
        let r = v.round();
        if (v - r).abs() > 1e-9 {
            return Some("a non-integer entry");
        }
        if r.abs() >= 9_223_372_036_854_775_808.0 {
            return Some("an integer beyond i64");
        }
    }
    None
}

fn balance(m: &mut Module, mat: &[Vec<f64>]) {
    let outs = match balance_equation(bb(mat)) {
        Ok(x) => {
            let mut v = vec![Tok::V(x.len())];
            v.extend(x.iter().map(|&k| i(k)));
            v
        }
        Err(e_) => vec![err(&e_)],
    };
    let ins = matrix_toks(mat);
    if let Some(why) = no_counterpart(mat) {
        assert!(
            matches!(&outs[..], [Tok::Err(v)] if v == "InvalidInput"),
            "a matrix without an integer counterpart that Rust does not refuse: {mat:?}"
        );
        let note = format!(
            "ADR 0002 #7: no integer-matrix counterpart ({why}): the port's balance_equation takes integers, and \
             Rust's Err(InvalidInput) is this input's result"
        );
        m.row_noted("stoichiometry::balance_equation", Class::A, &ins, &outs, &note);
    } else {
        m.row("stoichiometry::balance_equation", Class::A, &ins, &outs);
    }
}

fn mat(rows: &[&[f64]]) -> Vec<Vec<f64>> {
    rows.iter().map(|r| r.to_vec()).collect()
}

/// Every permutation of 0 .. n (Heap's algorithm order is not needed: lexicographic, so the rows are stable).
fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut out = vec![];
    let mut p: Vec<usize> = (0..n).collect();
    loop {
        out.push(p.clone());
        let Some(k) = (0..n.saturating_sub(1)).rev().find(|&k| p[k] < p[k + 1]) else {
            return out;
        };
        let l = (k + 1..n).rev().find(|&l| p[k] < p[l]).unwrap_or(k);
        p.swap(k, l);
        p[k + 1..].reverse();
    }
}

fn rows_in(m0: &[Vec<f64>], order: &[usize]) -> Vec<Vec<f64>> {
    order.iter().map(|&k| m0[k].clone()).collect()
}

fn cols_in(m0: &[Vec<f64>], order: &[usize]) -> Vec<Vec<f64>> {
    m0.iter().map(|r| order.iter().map(|&k| r[k]).collect()).collect()
}

fn negated(m0: &[Vec<f64>]) -> Vec<Vec<f64>> {
    m0.iter().map(|r| r.iter().map(|&x| -x).collect()).collect()
}

/// A species as (element index, count) pairs; `charge` is the last element row when the equation is ionic.
type Species = &'static [(usize, i32)];

/// The composition matrix of reactants → products over `elements` element rows (reactants positive, products negative).
fn equation(elements: usize, reactants: &[Species], products: &[Species]) -> Vec<Vec<f64>> {
    let n = reactants.len() + products.len();
    let mut out = vec![vec![0.0; n]; elements];
    for (j, sp) in reactants.iter().chain(products).enumerate() {
        let sign = if j < reactants.len() { 1.0 } else { -1.0 };
        for &(e, c) in *sp {
            out[e][j] += sign * f64::from(c);
        }
    }
    out
}

/// Textbook equations, each with its element rows (a charge row last where ionic).
fn textbook() -> Vec<Vec<Vec<f64>>> {
    // element indices per equation, in the order listed in each comment
    vec![
        // Fe + O₂ → Fe₂O₃ (Fe, O)
        equation(2, &[&[(0, 1)], &[(1, 2)]], &[&[(0, 2), (1, 3)]]),
        // Al + O₂ → Al₂O₃
        equation(2, &[&[(0, 1)], &[(1, 2)]], &[&[(0, 2), (1, 3)]]),
        // N₂ + H₂ → NH₃ (N, H)
        equation(2, &[&[(0, 2)], &[(1, 2)]], &[&[(0, 1), (1, 3)]]),
        // KClO₃ → KCl + O₂ (K, Cl, O)
        equation(3, &[&[(0, 1), (1, 1), (2, 3)]], &[&[(0, 1), (1, 1)], &[(2, 2)]]),
        // C₆H₁₂O₆ + O₂ → CO₂ + H₂O (C, H, O)
        equation(3, &[&[(0, 6), (1, 12), (2, 6)], &[(2, 2)]], &[&[(0, 1), (2, 2)], &[(1, 2), (2, 1)]]),
        // CO₂ + H₂O → C₆H₁₂O₆ + O₂ (photosynthesis)
        equation(3, &[&[(0, 1), (2, 2)], &[(1, 2), (2, 1)]], &[&[(0, 6), (1, 12), (2, 6)], &[(2, 2)]]),
        // Fe₂O₃ + Al → Al₂O₃ + Fe (Fe, O, Al: thermite)
        equation(3, &[&[(0, 2), (1, 3)], &[(2, 1)]], &[&[(2, 2), (1, 3)], &[(0, 1)]]),
        // Cu + HNO₃ → Cu(NO₃)₂ + NO + H₂O (Cu, H, N, O)
        equation(
            4,
            &[&[(0, 1)], &[(1, 1), (2, 1), (3, 3)]],
            &[&[(0, 1), (2, 2), (3, 6)], &[(2, 1), (3, 1)], &[(1, 2), (3, 1)]],
        ),
        // Ca₃(PO₄)₂ + SiO₂ + C → CaSiO₃ + P₄ + CO (Ca, P, O, Si, C)
        equation(
            5,
            &[&[(0, 3), (1, 2), (2, 8)], &[(3, 1), (2, 2)], &[(4, 1)]],
            &[&[(0, 1), (3, 1), (2, 3)], &[(1, 4)], &[(4, 1), (2, 1)]],
        ),
        // MnO₄⁻ + Fe²⁺ + H⁺ → Mn²⁺ + Fe³⁺ + H₂O (Mn, O, Fe, H, charge)
        equation(
            5,
            &[&[(0, 1), (1, 4), (4, -1)], &[(2, 1), (4, 2)], &[(3, 1), (4, 1)]],
            &[&[(0, 1), (4, 2)], &[(2, 1), (4, 3)], &[(3, 2), (1, 1)]],
        ),
        // Cr₂O₇²⁻ + I⁻ + H⁺ → Cr³⁺ + I₂ + H₂O (Cr, O, I, H, charge)
        equation(
            5,
            &[&[(0, 2), (1, 7), (4, -2)], &[(2, 1), (4, -1)], &[(3, 1), (4, 1)]],
            &[&[(0, 1), (4, 3)], &[(2, 2)], &[(3, 2), (1, 1)]],
        ),
        // K₄Fe(CN)₆ + KMnO₄ + H₂SO₄ → KHSO₄ + Fe₂(SO₄)₃ + MnSO₄ + HNO₃ + CO₂ + H₂O (K, Fe, C, N, Mn, O, H, S)
        equation(
            8,
            &[&[(0, 4), (1, 1), (2, 6), (3, 6)], &[(0, 1), (4, 1), (5, 4)], &[(6, 2), (7, 1), (5, 4)]],
            &[
                &[(0, 1), (6, 1), (7, 1), (5, 4)],
                &[(1, 2), (7, 3), (5, 12)],
                &[(4, 1), (7, 1), (5, 4)],
                &[(6, 1), (3, 1), (5, 3)],
                &[(2, 1), (5, 2)],
                &[(6, 2), (5, 1)],
            ],
        ),
        // NaOH + H₂SO₄ → Na₂SO₄ + H₂O (Na, O, H, S)
        equation(4, &[&[(0, 1), (1, 1), (2, 1)], &[(2, 2), (3, 1), (1, 4)]], &[&[(0, 2), (3, 1), (1, 4)], &[(2, 2), (1, 1)]]),
        // C₃H₈ + O₂ → CO₂ + H₂O
        equation(3, &[&[(0, 3), (1, 8)], &[(2, 2)]], &[&[(0, 1), (2, 2)], &[(1, 2), (2, 1)]]),
        // P₄O₁₀ + H₂O → H₃PO₄ (P, O, H)
        equation(3, &[&[(0, 4), (1, 10)], &[(2, 2), (1, 1)]], &[&[(2, 3), (0, 1), (1, 4)]]),
        // Pb(NO₃)₂ → PbO + NO₂ + O₂ (Pb, N, O)
        equation(3, &[&[(0, 1), (1, 2), (2, 6)]], &[&[(0, 1), (2, 1)], &[(1, 1), (2, 2)], &[(2, 2)]]),
        // H₂ + O₂ → H₂O + H₂O₂: two independent reactions (ambiguous)
        equation(2, &[&[(0, 2)], &[(1, 2)]], &[&[(0, 2), (1, 1)], &[(0, 2), (1, 2)]]),
        // C + O₂ → CO + CO₂: ambiguous
        equation(2, &[&[(0, 1)], &[(1, 2)]], &[&[(0, 1), (1, 1)], &[(0, 1), (1, 2)]]),
        // H₂O → H₂O₂: only the trivial solution
        equation(2, &[&[(0, 2), (1, 1)]], &[&[(0, 2), (1, 2)]]),
    ]
}

/// The combustion C_aH_bO_c + O₂ → CO₂ + H₂O (rows C, H, O; species fuel, O₂, CO₂, H₂O).
fn combustion(a: u32, b: u32, c: u32) -> Vec<Vec<f64>> {
    let (a, b, c) = (f64::from(a), f64::from(b), f64::from(c));
    vec![vec![a, 0.0, -1.0, 0.0], vec![b, 0.0, 0.0, -2.0], vec![c, 2.0, -2.0, -1.0]]
}

/// A chain A₀ ⇌ … of n species: row k balances species k against k + 1 with the coefficients x[k], x[k + 1]
/// (x[k + 1]·A_k − x[k]·A_{k+1}), so the answer is x / gcd(x).
fn chain(x: &[i64]) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut out = vec![];
    for k in 0..n - 1 {
        let mut row = vec![0.0; n];
        #[allow(clippy::cast_precision_loss)]
        {
            row[k] = x[k + 1] as f64;
            row[k + 1] = -(x[k] as f64);
        }
        out.push(row);
    }
    out
}

/// A consistent matrix for a positive x: each row a random small combination of x[k]·e_j − x[j]·e_k, so M·x = 0; with
/// m ≥ n − 1 rows the null space is generically x alone.
fn consistent(rng: &mut Rng, x: &[i64], m: usize, terms: usize, coef: i64) -> Vec<Vec<f64>> {
    let n = x.len();
    let mut out = vec![];
    for _ in 0..m {
        let mut row = vec![0_i64; n];
        for _ in 0..terms {
            let j = usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)).unwrap_or(0);
            let k = usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)).unwrap_or(0);
            if j == k {
                continue;
            }
            let t = rng.int(-coef, coef);
            row[j] += t * x[k];
            row[k] -= t * x[j];
        }
        #[allow(clippy::cast_precision_loss)]
        out.push(row.iter().map(|&v| v as f64).collect());
    }
    out
}

/// Rust's acceptance of one entry: finite and within 10⁻⁹ of `round(v)`.
fn accepted(v: f64) -> bool {
    v.is_finite() && (v - v.round()).abs() <= 1e-9
}

/// Probes on Rust's 10⁻⁹ integer tolerance for an integer r ≥ 0: in each direction from ±r, the farthest entry Rust
/// accepts and the next double, which it refuses (found by bisection over the bit patterns: the deviation is a
/// multiple of the entry's ulp, so the boundary is exact).
fn tolerance_probes(r: f64) -> Vec<f64> {
    let mut out = vec![];
    for up in [true, false] {
        if r == 0.0 && !up {
            continue; // below 0 is −0's side: the sign flip below covers it
        }
        let b = r.to_bits();
        let at = |k: u64| f64::from_bits(if up { b + k } else { b - k });
        let far = if up { r + 0.25 } else { r - 0.25 };
        let (mut lo, mut hi) = (0_u64, if up { far.to_bits() - b } else { b - far.to_bits() });
        while hi - lo > 1 {
            let mid = lo + (hi - lo) / 2;
            if accepted(at(mid)) { lo = mid } else { hi = mid }
        }
        for v in [at(lo), at(hi)] {
            out.push(v);
            out.push(-v);
        }
    }
    out
}

fn emit_balance(m: &mut Module) {
    // The Rust tests (rust-old/src/stoichiometry.rs), in order.
    let big = 9_007_199_254_740_991.0; // 2⁵³ − 1
    let ch4 = mat(&[&[1.0, 0.0, -1.0, 0.0], &[4.0, 0.0, 0.0, -2.0], &[0.0, 2.0, -2.0, -1.0]]);
    let cu: &[f64] = &[1.0, 0.0, 0.0, -1.0];
    let s_: &[f64] = &[1.0, 0.0, -1.0, 0.0];
    let o: &[f64] = &[4.0, 0.0, -4.0, 0.0];
    let fe: &[f64] = &[0.0, 1.0, -1.0, 0.0];
    let cuso4 = mat(&[cu, s_, o, fe]);
    let kmno4 = mat(&[
        &[1.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        &[1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
        &[4.0, 0.0, 0.0, 0.0, -1.0, 0.0],
        &[0.0, 1.0, 0.0, 0.0, -2.0, 0.0],
        &[0.0, 1.0, -1.0, -2.0, 0.0, -2.0],
    ]);
    let tests: Vec<Vec<Vec<f64>>> = vec![
        mat(&[&[2.0, 0.0, -2.0], &[0.0, 2.0, -1.0]]),
        ch4.clone(),
        negated(&ch4),
        cuso4.clone(),
        mat(&[fe, cu, s_, o]),
        kmno4.clone(),
        mat(&[&[1.0, 0.0, 0.0, -1.0], &[1.0, 0.0, -1.0, 0.0], &[0.0, 1.0, -1.0, 0.0], &[0.0, 1.0, 0.0, -1.0], &[
            0.0, 3.0, 0.0, -3.0,
        ]]),
        mat(&[&[1.0, 0.0, -1.0, 0.0], &[4.0, 0.0, 0.0, 2.0], &[0.0, 2.0, -2.0, 1.0]]),
        mat(&[&[2.0, 0.0, -2.0, -2.0], &[0.0, 2.0, -1.0, -2.0]]),
        mat(&[&[2.0, 0.0], &[0.0, -2.0]]),
        mat(&[&[2.0, 0.0, -2.0, 0.0], &[0.0, 2.0, -1.0, 0.0], &[0.0, 0.0, 0.0, 1.0]]),
        mat(&[&[1.5, -1.0]]),
        mat(&[&[f64::NAN, -1.0]]),
        mat(&[&[f64::INFINITY, -1.0]]),
        mat(&[&[1e300, -1.0]]),
        mat(&[&[2.0 + 1e-12, 0.0, -2.0], &[0.0, 2.0, -1.0 - 1e-12]]),
        mat(&[&[1.0, -8_589_934_592.0]]),
        mat(&[&[big, 1.0, -1.0], &[1.0, big, -1.0]]),
        mat(&[&[big, -(big - 2.0), 3.0, -7.0], &[big - 4.0, 5.0, -big, 11.0], &[13.0, big - 6.0, -17.0, -big]]),
        vec![],
        mat(&[&[1.0]]),
    ];
    for t in &tests {
        balance(m, t);
    }
    // CH₄ in every row order, negated, and with its columns permuted; CuSO₄ + Fe in every row order (port plan §4:
    // both orders); KMnO₄ + HCl in every row order and with its columns permuted.
    for p in permutations(3) {
        balance(m, &rows_in(&ch4, &p));
        balance(m, &negated(&rows_in(&ch4, &p)));
    }
    for p in permutations(4) {
        balance(m, &cols_in(&ch4, &p));
        balance(m, &rows_in(&cuso4, &p));
    }
    for p in permutations(5) {
        balance(m, &rows_in(&kmno4, &p));
    }
    for p in permutations(6).iter().step_by(7) {
        balance(m, &cols_in(&kmno4, p));
    }
    for t in textbook() {
        balance(m, &t);
        balance(m, &negated(&t));
        let rev: Vec<usize> = (0..t.len()).rev().collect();
        balance(m, &rows_in(&t, &rev));
    }
    // Every combustion C_aH_bO_c + O₂ → CO₂ + H₂O with a ≤ 10, b ≤ 2a + 2, c ≤ a + 1: balanced, mixed sign (too much
    // fuel oxygen), the zero O₂ coefficient (c = 2a + b/2), and odd b.
    for a in 1..=10 {
        for b in 0..=2 * a + 2 {
            for c in 0..=a + 1 {
                balance(m, &combustion(a, b, c));
            }
        }
    }
    // Chains: 30 species (29 × 30 = 870 cells, past the port's 512-word table) and 70 (the 64-word null vector).
    for n in [2_i64, 3, 10, 30, 64, 65, 70] {
        let x: Vec<i64> = (0..n).map(|k| 1 + (k * 7) % 5).collect();
        balance(m, &chain(&x));
        let y: Vec<i64> = (0..n).map(|k| 1 + k % 3).collect();
        balance(m, &chain(&y));
    }
    // Seeded consistent matrices: x random positive, rows orthogonal to it.
    let mut rng = Rng::new(0x5f01);
    for _ in 0..700 {
        let n = usize::try_from(rng.int(2, 9)).unwrap_or(2);
        let hi = [3, 12, 40, 1000][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        let x: Vec<i64> = (0..n).map(|_| rng.int(1, hi)).collect();
        let extra = usize::try_from(rng.int(0, 2)).unwrap_or(0);
        let mm = n - 1 + extra;
        let terms = usize::try_from(rng.int(1, 4)).unwrap_or(1);
        let t = consistent(&mut rng, &x, mm, terms, 3);
        balance(m, &t);
    }
    // Large coefficients: x up to 2³² − 1 and past it (the u32 conversion), and entries near 2⁵³ (the elimination's
    // checked i64 products).
    for _ in 0..120 {
        let n = usize::try_from(rng.int(2, 4)).unwrap_or(2);
        let hi = [65_535_i64, 4_294_967_295, 8_589_934_592, 1 << 40][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        let x: Vec<i64> = (0..n).map(|_| rng.int(1, hi)).collect();
        let t = consistent(&mut rng, &x, n - 1, 1, 1);
        balance(m, &t);
    }
    for (x0, x1) in [(4_294_967_295_i64, 1_i64), (4_294_967_296, 1), (4_294_967_295, 4_294_967_294), (1, 4_294_967_297)] {
        #[allow(clippy::cast_precision_loss)]
        balance(m, &[vec![x1 as f64, -(x0 as f64)]]);
    }
    // Intermediate i64 overflow although a valid u32 solution exists: rows [a, b, −(a + b)] with a, b near 2⁴⁰ … 2⁵²
    // (x = [1, 1, 1]), and rows near 2³² (products past 2⁶³ only after the gcd reduction fails to shrink them).
    for _ in 0..60 {
        let e = rng.int(30, 51);
        let lo = 1_i64 << e;
        let a = rng.int(lo, 2 * lo);
        let b = rng.int(lo, 2 * lo);
        let c = rng.int(lo, 2 * lo);
        let d = rng.int(lo, 2 * lo);
        #[allow(clippy::cast_precision_loss)]
        balance(m, &[vec![a as f64, b as f64, -((a + b) as f64)], vec![c as f64, d as f64, -((c + d) as f64)]]);
    }
    // The reduced row whose gcd is 2⁶³ (Rust's `i64::try_from(g)` overflow): row 1 becomes [0, −2⁶³, 0].
    let p31 = 2_147_483_648.0;
    balance(m, &[vec![p31, 1.0, 0.0], vec![p31, -4_294_967_295.0, 0.0]]);
    balance(m, &[vec![p31, 1.0, 0.0], vec![p31, -4_294_967_295.0, 0.0], vec![0.0, 0.0, 1.0]]);
    // Seeded random matrices: whatever Rust returns (the trivial, ambiguous, zero and mixed-sign classes mostly).
    for _ in 0..700 {
        let mm = usize::try_from(rng.int(1, 6)).unwrap_or(1);
        let n = usize::try_from(rng.int(2, 7)).unwrap_or(2);
        let r = rng.int(1, 5);
        #[allow(clippy::cast_precision_loss)]
        let t: Vec<Vec<f64>> = (0..mm).map(|_| (0..n).map(|_| rng.int(-r, r) as f64).collect()).collect();
        balance(m, &t);
    }
    // Error classes by construction.
    let e2 = pow2(53);
    for t in [
        vec![vec![]],
        vec![vec![], vec![]],
        vec![vec![0.0, 0.0]],
        vec![vec![0.0, 0.0], vec![0.0, 0.0]],
        vec![vec![1.0, -1.0], vec![1.0]],
        vec![vec![1.0], vec![1.0, -1.0]],
        vec![vec![1.0, -1.0], vec![1.0, -1.0, 0.0]],
        vec![vec![1.0, -1.0, 0.5], vec![1.0, -1.0]],
        vec![vec![1.0, -1.0], vec![f64::NAN]],
        vec![vec![e2, -e2]],
        vec![vec![e2, -1.0], vec![1.0, -e2]],
        vec![vec![-e2, e2, 0.0], vec![1.0, 0.0, -1.0]],
        vec![vec![e2 + 2.0, -1.0]],
        vec![vec![1.0, -(e2 + 2.0)]],
        vec![vec![pow2(62), -1.0]],
        vec![vec![-pow2(62), 1.0]],
        vec![vec![pow2(63), -1.0]],
        vec![vec![-pow2(63), 1.0]],
        vec![vec![f64::MAX, -1.0]],
        vec![vec![1e18, -1e18]],
        vec![vec![1.0, -1.0, 2.0, -2.0]],
        vec![vec![1.0, 1.0]],
        vec![vec![1.0, 0.0]],
        vec![vec![0.0, 1.0]],
        vec![vec![1.0, -1.0], vec![1.0, -2.0]],
        vec![vec![2.0, -4.0]],
        vec![vec![-0.0, 0.0, 1.0, -1.0], vec![1.0, -1.0, 0.0, 0.0]],
    ] {
        balance(m, &t);
    }
    for x in specials().into_iter().chain([0.5, -0.5, 2.5, -2.5, 1e-9, 1e-10, 1.0 + 1e-6, 3.000_000_001]) {
        balance(m, &[vec![x, -1.0], vec![1.0, -1.0]]);
        balance(m, &[vec![1.0, -1.0, 0.0], vec![0.0, 1.0, x]]);
        balance(m, &[vec![2.0, 0.0, -2.0], vec![0.0, 2.0, x]]);
    }
    // Entries on Rust's 10⁻⁹ integer tolerance, to the ulp, at 0, 1, 2, 2²⁰ and 2²⁹ (the deviation is a multiple of
    // the entry's ulp): the largest deviation accepted and the next one, both signs and directions, in H₂ + O₂.
    for r in [0.0, 1.0, 2.0, pow2(20), pow2(29)] {
        for v in tolerance_probes(r) {
            balance(m, &[vec![2.0, 0.0, -2.0], vec![0.0, 2.0, v]]);
            balance(m, &[vec![v, -1.0]]);
        }
    }
}

// ── assign_oxidation_states ─────────────────────────────────────────────────────────

fn is_metal(c: ElementCategory) -> bool {
    matches!(
        c,
        ElementCategory::AlkaliMetal
            | ElementCategory::AlkalineEarth
            | ElementCategory::TransitionMetal
            | ElementCategory::PostTransitionMetal
            | ElementCategory::Lanthanide
            | ElementCategory::Actinide
    )
}

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// ADR 0002 #7: `exact_oxidation_state` with rationals — num/den reduced (den > 0), or Rust's Err variant. Rust's
/// checks and order; where Rust returns Err "non-integer average oxidation state", the exact average; a value outside
/// i32 is Err(ComputationError) (Rust's "out of range", applied to the rational).
fn exact_state(num: i64, den: i64) -> Result<(i64, i64), &'static str> {
    if den == 0 {
        return Err("InvalidInput");
    }
    let (n, d) = (i128::from(num), i128::from(den));
    let g = i128::try_from(gcd(n.unsigned_abs(), d.unsigned_abs())).map_err(|_| "ComputationError")?;
    let (mut p, mut q) = (n / g, d / g);
    if q < 0 {
        p = -p;
        q = -q;
    }
    let fl = p.div_euclid(q);
    let out_of_range = if p % q == 0 { fl < i128::from(i32::MIN) || fl > i128::from(i32::MAX) } else {
        fl >= i128::from(i32::MAX) || fl < i128::from(i32::MIN)
    };
    if out_of_range {
        return Err("ComputationError");
    }
    Ok((i64::try_from(p).map_err(|_| "ComputationError")?, i64::try_from(q).map_err(|_| "ComputationError")?))
}

fn rule_state(z: u8, c: ElementCategory, has_f: bool, hydride: bool) -> Option<i64> {
    match z {
        9 => Some(-1),
        1 if hydride => Some(-1),
        1 => Some(1),
        8 if has_f => None,
        8 => Some(-2),
        _ => match c {
            ElementCategory::AlkaliMetal => Some(1),
            ElementCategory::AlkalineEarth => Some(2),
            _ => None,
        },
    }
}

/// The port's `assign_oxidation_states` (ADR 0002 #7), Rust 1.2.0's algorithm with exact rational states: one
/// (z, p, q) per entry, or Rust's Err variant. The integer arithmetic is exact (i128); Rust's checked-i64 overflows
/// can't occur on these inputs (fewer than 2³⁰ entries), which the oracle comparison confirms on every row.
pub fn port_oxidation(atoms: &[(u8, u32)], charge: i32) -> Result<Vec<(u8, i64, i64)>, &'static str> {
    let Some(&(first_z, _)) = atoms.first() else {
        return Err("InvalidInput");
    };
    let mut single = true;
    let mut has_f = false;
    let mut others_metal = true;
    let mut total: i64 = 0;
    for &(z, count) in atoms {
        if count == 0 {
            return Err("InvalidInput");
        }
        let Some(e) = lookup_by_number(z) else {
            return Err("InvalidElement");
        };
        single &= z == first_z;
        has_f |= z == 9;
        if z != 1 && !is_metal(e.category) {
            others_metal = false;
        }
        total += i64::from(count);
    }
    let cat = |z: u8| lookup_by_number(z).map_or(ElementCategory::Nonmetal, |e| e.category);
    if single {
        let (p, q) = exact_state(i64::from(charge), total)?;
        return Ok(atoms.iter().map(|&(z, _)| (z, p, q)).collect());
    }
    let mut assigned: i64 = 0;
    let mut oxygen: i64 = 0;
    let mut unknown: Option<(u8, i64)> = None;
    for &(z, count) in atoms {
        match rule_state(z, cat(z), has_f, others_metal) {
            Some(os) => {
                assigned += os * i64::from(count);
                if z == 8 {
                    oxygen += i64::from(count);
                }
            }
            None => {
                unknown = match unknown {
                    None => Some((z, i64::from(count))),
                    Some((uz, uc)) if uz == z => Some((uz, uc + i64::from(count))),
                    Some(_) => return Err("ComputationError"),
                };
            }
        }
    }
    let remaining = i64::from(charge) - assigned;
    let solved = match unknown {
        Some((uz, uc)) => Some((uz, exact_state(remaining, uc)?)),
        None if remaining == 0 => None,
        None if oxygen > 0 => {
            let (p, q) = exact_state(remaining - 2 * oxygen, oxygen)?;
            if !(-2 * q < p && p < 0) {
                return Err("ComputationError");
            }
            Some((8, (p, q)))
        }
        None => return Err("ComputationError"),
    };
    Ok(atoms
        .iter()
        .map(|&(z, _)| match solved {
            Some((sz, (p, q))) if sz == z => (z, p, q),
            _ => (z, rule_state(z, cat(z), has_f, others_metal).unwrap_or(0), 1),
        })
        .collect())
}

fn oxidation(m: &mut Module, atoms: &[(u8, u32)], charge: i32) {
    let oracle = assign_oxidation_states(bb(atoms), bb(charge));
    let mut ins = vec![Tok::V(atoms.len())];
    for &(z, c) in atoms {
        ins.push(i(z));
        ins.push(i(c));
    }
    ins.push(i(charge));
    let port = port_oxidation(atoms, charge);
    let outs = match &oracle {
        Ok(v) => {
            let mut t = vec![Tok::V(v.len())];
            for &(z, os) in v {
                t.push(i(z));
                t.push(i(os));
            }
            t
        }
        Err(e_) => vec![err(e_)],
    };
    let non_integer = matches!(&oracle, Err(KimiyaError::ComputationError(msg)) if msg == "non-integer average oxidation state");
    match (&oracle, &port) {
        (Ok(v), Ok(p)) => assert!(
            v.len() == p.len() && v.iter().zip(p).all(|(&(z, os), &(pz, pp, pq))| z == pz && i64::from(os) == pp && pq == 1),
            "port_oxidation disagrees with the oracle: {atoms:?} {charge}"
        ),
        (Err(e_), Err(pe)) => assert!(
            matches!(&err(e_), Tok::Err(v) if v == pe),
            "port_oxidation's Err differs: {atoms:?} {charge}"
        ),
        (Err(_), Ok(_)) => assert!(non_integer, "port_oxidation is Ok where the oracle's Err is not the non-integer one: {atoms:?} {charge}"),
        (Ok(_), Err(_)) => panic!("port_oxidation errs where the oracle is Ok: {atoms:?} {charge}"),
    }
    if let (true, Ok(p)) = (non_integer, &port) {
        let mut toks = vec![Tok::V(p.len())];
        for &(z, pp, pq) in p {
            toks.push(i(z));
            toks.push(i(pp));
            toks.push(i(pq));
        }
        let note = format!(
            "ADR 0002 #7: the port returns the exact average oxidation states (Z, numerator, denominator) {}",
            tokens_note(&toks)
        );
        m.row_noted("stoichiometry::assign_oxidation_states", Class::A, &ins, &outs, &note);
    } else {
        m.row("stoichiometry::assign_oxidation_states", Class::A, &ins, &outs);
    }
}

fn emit_oxidation(m: &mut Module) {
    // The 1.2 tests' species (port plan §4's oxidation corpus), in order.
    let big = 1_u32 << 31;
    let tests: &[(&[(u8, u32)], i32)] = &[
        (&[(1, 2), (8, 1)], 0),
        (&[(11, 1), (17, 1)], 0),
        (&[(16, 1), (8, 4)], -2),
        (&[(25, 1), (8, 4)], -1),
        (&[(19, 1), (25, 1), (8, 4)], 0),
        (&[(24, 2), (8, 7)], -2),
        (&[(7, 1), (1, 4)], 1),
        (&[(6, 1), (8, 2)], 0),
        (&[(1, 2), (8, 2)], 0),
        (&[(11, 2), (8, 2)], 0),
        (&[(19, 1), (8, 2)], 0),
        (&[(8, 1), (9, 2)], 0),
        (&[(11, 1), (1, 1)], 0),
        (&[(20, 1), (1, 2)], 0),
        (&[(3, 1), (13, 1), (1, 4)], 0),
        (&[(8, 2)], 0),
        (&[(9, 2)], 0),
        (&[(17, 1)], -1),
        (&[(8, 2)], -2),
        (&[(8, 2)], -1),
        (&[(6, 3), (1, 8)], 0),
        (&[(26, 3), (8, 4)], 0),
        (&[(6, 1), (1, 6), (6, 1)], 0),
        (&[(11, 1), (9, 1)], 1),
        (&[(1, 2), (8, 1)], -1),
        (&[(1, 2), (8, 1)], 1),
        (&[(8, 1), (1, 1)], -1),
        (&[(1, 1), (8, 2)], -1),
        (&[(6, 1), (7, 1)], -1),
        (&[(25, 1), (8, big)], 0),
        (&[(25, 1), (8, u32::MAX)], 0),
        (&[(11, u32::MAX), (17, u32::MAX)], 0),
        (&[(8, u32::MAX), (8, u32::MAX)], i32::MIN),
        (&[(25, 0), (8, 4)], -1),
        (&[(200, 1), (8, 1)], 0),
        (&[], 0),
    ];
    for &(a, q) in tests {
        oxidation(m, a, q);
    }
    // Every element alone (rule 1): counts and charges, the averages included.
    for z in 0..=255_u8 {
        let counts: &[u32] = if (1..=118).contains(&z) { &[1, 2, 3, 4, 6, 8] } else { &[1] };
        for &c in counts {
            for q in -3..=3 {
                oxidation(m, &[(z, c)], q);
            }
        }
    }
    // Binary oxides, hydrides and fluorides of every element, both entry orders for the oxides.
    for z in 1..=118_u8 {
        if z != 8 {
            for zc in 1..=3 {
                for oc in 1..=7 {
                    for q in [-2, -1, 0, 1] {
                        oxidation(m, &[(z, zc), (8, oc)], q);
                    }
                    oxidation(m, &[(8, oc), (z, zc)], 0);
                }
            }
        }
        if z != 1 {
            for zc in 1..=2 {
                for hc in 1..=4 {
                    oxidation(m, &[(z, zc), (1, hc)], 0);
                }
            }
            oxidation(m, &[(z, 1), (1, 4)], -1);
        }
        if z != 9 {
            for fc in 1..=7 {
                oxidation(m, &[(z, 1), (9, fc)], 0);
            }
        }
    }
    // Oxyanions, salts and peroxides / superoxides / ozonides (the re-solve branch and its −2 < os < 0 test).
    let more: &[(&[(u8, u32)], i32)] = &[
        (&[(7, 1), (8, 3)], -1),
        (&[(7, 1), (8, 2)], -1),
        (&[(15, 1), (8, 4)], -3),
        (&[(17, 1), (8, 4)], -1),
        (&[(17, 1), (8, 1)], -1),
        (&[(24, 1), (8, 4)], -2),
        (&[(16, 2), (8, 3)], -2),
        (&[(16, 4), (8, 6)], -2),
        (&[(6, 1), (8, 3)], -2),
        (&[(6, 2), (8, 4)], -2),
        (&[(11, 1), (16, 1), (8, 4), (1, 1)], 0),
        (&[(19, 2), (24, 2), (8, 7)], 0),
        (&[(11, 2), (16, 2), (8, 3)], 0),
        (&[(11, 2), (16, 4), (8, 6)], 0),
        (&[(26, 2), (8, 3)], 0),
        (&[(26, 1), (8, 1)], 0),
        (&[(29, 2), (8, 1)], 0),
        (&[(1, 1), (8, 1)], 0),
        (&[(1, 1), (8, 2)], 0),
        (&[(1, 2), (8, 3)], 0),
        (&[(1, 1), (8, 3)], -1),
        (&[(11, 1), (8, 1), (1, 1)], 0),
        (&[(20, 1), (8, 2), (1, 2)], 0),
        (&[(56, 1), (8, 2)], 0),
        (&[(20, 1), (8, 4)], 0),
        (&[(8, 2), (9, 2)], 0),
        (&[(8, 1), (9, 1)], 0),
        (&[(8, 3)], 0),
        (&[(8, 3)], -1),
        (&[(8, 3)], -2),
        (&[(6, 1), (1, 4), (8, 1)], 0),
        (&[(6, 2), (1, 6), (8, 1)], 0),
        (&[(6, 6), (1, 12), (8, 6)], 0),
        (&[(6, 4), (1, 10)], 0),
        (&[(6, 8), (1, 18)], 0),
        (&[(5, 1), (1, 4), (11, 1)], 0),
        (&[(14, 1), (1, 4)], 0),
        (&[(3, 1), (1, 1)], 0),
        (&[(13, 1), (1, 3)], 0),
        (&[(1, 1), (9, 1), (8, 1)], 0),
        (&[(1, 2), (8, 1), (1, 2)], 0),
        (&[(8, 1), (1, 2), (8, 1)], 0),
        (&[(26, 1), (26, 2), (8, 4)], 0),
        (&[(6, 1), (6, 1), (6, 1), (1, 8)], 0),
        (&[(6, 1), (7, 1), (8, 1)], -1),
        (&[(6, 1), (1, 1), (7, 1)], 0),
        (&[(11, 1), (19, 1), (8, 1)], 0),
        (&[(11, 1), (17, 1), (1, 1)], 0),
    ];
    for &(a, q) in more {
        oxidation(m, a, q);
    }
    for z in [3_u8, 11, 19, 37, 55, 87, 4, 12, 20, 38, 56, 88] {
        for mc in 1..=2 {
            for oc in 1..=4 {
                for q in -1..=1 {
                    oxidation(m, &[(z, mc), (8, oc)], q);
                }
            }
        }
    }
    // The u32 and i32 extremes and the i32 range's edges: integer states at i32::MAX and past it, averages just inside
    // and just outside (ADR 0002 #7: the range applies to the rational's value), the counts' sum.
    let umax = u32::MAX;
    let edges_: &[(&[(u8, u32)], i32)] = &[
        (&[(25, 1), (8, 1 << 30)], -1),
        (&[(25, 1), (8, 1 << 30)], 0),
        (&[(25, 1), (8, 1 << 30)], 1),
        (&[(7, 1), (1, umax)], 0),
        (&[(7, 1), (1, 1 << 31)], 0),
        (&[(7, 1), (1, (1 << 31) - 1)], -1),
        (&[(7, 1), (1, (1 << 31) - 1)], -2),
        (&[(25, 2), (8, (1 << 31) - 1)], -1),
        (&[(25, 2), (8, (1 << 31) - 1)], 1),
        (&[(25, 2), (8, (1 << 31) - 1)], 0),
        (&[(7, 2), (1, umax)], 0),
        (&[(7, 2), (1, umax), (1, 2)], 0),
        (&[(7, 2), (1, umax), (1, 1)], 0),
        (&[(25, 3), (8, big)], 0),
        (&[(25, 3), (8, umax)], 0),
        (&[(25, 3), (8, umax)], 1),
        (&[(25, 3), (8, umax), (8, 1)], 0),
        (&[(25, 7), (8, umax)], 0),
        (&[(7, 3), (1, umax), (1, umax)], 0),
        (&[(7, 3), (1, umax), (1, umax), (1, 1)], 0),
        (&[(25, umax), (8, 1)], 0),
        (&[(25, umax), (8, umax)], 0),
        (&[(25, umax), (25, umax), (8, 3)], 0),
        (&[(8, umax)], i32::MAX),
        (&[(8, umax)], i32::MIN),
        (&[(8, 1)], i32::MAX),
        (&[(8, 1)], i32::MIN),
        (&[(8, 2)], i32::MAX),
        (&[(8, 2)], i32::MIN),
        (&[(11, umax), (8, umax)], 0),
        (&[(11, umax), (8, umax)], i32::MAX),
        (&[(19, 1), (8, umax)], 0),
        (&[(1, umax), (8, umax)], i32::MIN),
        (&[(1, 1), (8, 1)], i32::MAX),
        (&[(1, 1), (8, 1)], i32::MIN),
        (&[(9, umax), (8, 1)], 0),
        (&[(9, 1), (8, umax)], 0),
        (&[(0, 1)], 0),
        (&[(119, 1)], 0),
        (&[(255, 1), (8, 1)], 0),
        (&[(8, 1), (255, 1)], 0),
        (&[(8, 1), (0, 0)], 0),
        (&[(0, 0), (8, 1)], 0),
        (&[(8, 0)], 0),
        (&[(6, 1), (7, 1), (16, 1)], 0),
        (&[(6, 1), (7, 1), (0, 1)], 0),
        (&[(6, 1), (7, 1), (8, 0)], 0),
    ];
    for &(a, q) in edges_ {
        oxidation(m, a, q);
    }
    // Seeded species: 1–5 entries, common elements weighted, counts 1–12, charges −4 … 4.
    let pool: [u8; 30] = [1, 1, 1, 6, 6, 7, 8, 8, 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 22, 24, 25, 26, 29, 30, 35, 53, 56, 92, 2];
    let mut rng = Rng::new(0x5f02);
    for _ in 0..3000 {
        let n = usize::try_from(rng.int(1, 5)).unwrap_or(1);
        let atoms: Vec<(u8, u32)> = (0..n)
            .map(|_| {
                let z = if rng.int(0, 9) == 0 {
                    u8::try_from(rng.int(1, 118)).unwrap_or(1)
                } else {
                    pool[usize::try_from(rng.int(0, 29)).unwrap_or(0)]
                };
                (z, u32::try_from(rng.int(1, 12)).unwrap_or(1))
            })
            .collect();
        let q = i32::try_from(rng.int(-4, 4)).unwrap_or(0);
        oxidation(m, &atoms, q);
    }
}

// ── limiting reagent and yields ─────────────────────────────────────────────────────

fn slice_toks(xs: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(xs.len())];
    v.extend(xs.iter().map(|&x| f(x)));
    v
}

/// ADR 0002 #25: the port's limiting reagent: Rust's validation, in its order, with the ratios compared exactly
/// (mᵢ·c_{i₀} < m_{i₀}·cᵢ by `exact::ratio_lt`; the first index on an exact tie); None where Rust errs.
fn port_limiting(mo: &[f64], co: &[f64]) -> Option<usize> {
    if mo.len() != co.len() || mo.is_empty() {
        return None;
    }
    let mut best = 0;
    for k in 0..mo.len() {
        if !(mo[k].is_finite() && mo[k] >= 0.0) || !(co[k].is_finite() && co[k] > 0.0) {
            return None;
        }
        if k > 0 && ratio_lt(mo[k], co[k], mo[best], co[best]) {
            best = k;
        }
    }
    Some(best)
}

/// A class-B row of an f64 result, noted `[ADR 0002 #24: the port returns the correctly rounded <what> = <x> (f:…)]`
/// where the port's value (None: the oracle's Err, asserted) differs from the oracle's.
fn noted24(m: &mut Module, name: &str, ins: &[Tok], out: Result<f64, KimiyaError>, port: Option<f64>, what: &str) {
    assert_eq!(out.is_ok(), port.is_some(), "{name}: the port's validation differs from the oracle's: {ins:?}");
    if let (Ok(v), Some(c)) = (&out, port)
        && c.to_bits() != v.to_bits()
    {
        let note = format!("ADR 0002 #24: the port returns the correctly rounded {what} = {c:?} (f:{:016x})", c.to_bits());
        m.row_noted(name, Class::B, ins, &[f(*v)], &note);
        return;
    }
    m.row(name, Class::B, ins, &[res(out)]);
}

fn limiting(m: &mut Module, mo: &[f64], co: &[f64]) {
    let oracle = limiting_reagent(bb(mo), bb(co));
    let out = oracle.as_ref().map_or_else(err, |&k| i(i64::try_from(k).unwrap_or(-1)));
    let port = port_limiting(mo, co);
    assert_eq!(port.is_some(), oracle.is_ok(), "port_limiting's validation differs: {mo:?} {co:?}");
    let mut ins = slice_toks(mo);
    ins.extend(slice_toks(co));
    if let (Ok(k), Some(p)) = (&oracle, port)
        && *k != p
    {
        let note = format!(
            "ADR 0002 #25: the port returns the exact limiting reagent (mᵢ·c_j < m_j·cᵢ, the first on a tie) = {p} (i:{p})"
        );
        m.row_noted("stoichiometry::limiting_reagent", Class::A, &ins, &[out], &note);
        return;
    }
    m.row("stoichiometry::limiting_reagent", Class::A, &ins, &[out]);
}

fn yield_moles(m: &mut Module, mo: &[f64], co: &[f64], p: f64) {
    let out = theoretical_yield_moles(bb(mo), bb(co), bb(p));
    let port = (p.is_finite() && p > 0.0).then(|| port_limiting(mo, co)).flatten().map(|k| {
        let (mm, c) = (mo[k], co[k]);
        if mm == 0.0 { mm / c * p } else { cr_mul_div(mm, p, c) }
    });
    let mut ins = slice_toks(mo);
    ins.extend(slice_toks(co));
    ins.push(f(p));
    noted24(m, "stoichiometry::theoretical_yield_moles", &ins, out, port, "m·p/c at the exact limiting reagent");
}

/// ADR 0002 #24: percent_yield's row; the port's value is RN(100·a/t) (an actual yield of ±0: Rust's ±0).
fn pyield(m: &mut Module, a: f64, t: f64) {
    let out = percent_yield(bb(a), bb(t));
    let ok = t.is_finite() && t > 0.0 && a.is_finite() && a >= 0.0;
    let port = ok.then(|| if a == 0.0 { a / t * 100.0 } else { cr_mul_div(a, 100.0, t) });
    noted24(m, "stoichiometry::percent_yield", &[f(a), f(t)], out, port, "100·a/t");
}

fn emit_yields(m: &mut Module) {
    // limiting_reagent and theoretical_yield_moles: the Rust tests, each slot over the edges, ties, over- and underflow.
    let lists: &[(&[f64], &[f64])] = &[
        (&[3.0, 2.0], &[2.0, 1.0]),
        (&[10.0, 1.0], &[2.0, 1.0]),
        (&[], &[]),
        (&[1.0], &[1.0, 2.0]),
        (&[1.0, 2.0], &[1.0]),
        (&[f64::NAN, f64::NAN], &[1.0, 1.0]),
        (&[1.0, 1.0], &[f64::NAN, 1.0]),
        (&[1.0, 1.0], &[1.0, f64::INFINITY]),
        (&[f64::INFINITY, 1.0], &[1.0, 1.0]),
        (&[-1.0, 1.0], &[1.0, 1.0]),
        (&[1.0, 1.0], &[0.0, 1.0]),
        (&[5.0, 0.0], &[1.0, 1.0]),
        (&[1.0, 1.0], &[1.0, 1.0]),
        (&[2.0, 1.0], &[2.0, 1.0]),
        (&[0.3, 0.1], &[3.0, 1.0]),
        (&[0.1, 0.3], &[1.0, 3.0]),
        (&[0.0, 0.0, 0.0], &[1.0, 2.0, 3.0]),
        (&[-0.0, 0.0], &[1.0, 1.0]),
        (&[f64::MAX, f64::MAX], &[f64::from_bits(1), 1e-300]),
        (&[f64::MAX, 1e300], &[f64::from_bits(1), f64::from_bits(1)]),
        (&[f64::from_bits(1), f64::from_bits(2)], &[f64::MAX, f64::MAX]),
        (&[1e-300, 2e-300], &[1e300, 1e300]),
        (&[f64::from_bits(3), f64::from_bits(2)], &[3.0, 2.0]),
        (&[1.0, 2.0, 3.0, 0.5], &[2.0, 4.0, 6.0, 1.0]),
        (&[1e308, 1e308], &[1e-10, 2e-10]),
    ];
    for &(mo, co) in lists {
        limiting(m, mo, co);
        for p in [2.0, 1.0, 0.5] {
            yield_moles(m, mo, co, p);
        }
    }
    for x in edges() {
        limiting(m, &[x, 2.0], &[2.0, 1.0]);
        limiting(m, &[3.0, x], &[2.0, 1.0]);
        limiting(m, &[3.0, 2.0], &[x, 1.0]);
        limiting(m, &[3.0, 2.0], &[2.0, x]);
        yield_moles(m, &[x, 2.0], &[2.0, 1.0], 2.0);
        yield_moles(m, &[3.0, 2.0], &[2.0, x], 2.0);
        yield_moles(m, &[3.0, 2.0], &[2.0, 1.0], x);
        yield_moles(m, &[x], &[1.0], 3.0);
        yield_moles(m, &[1.0], &[x], 3.0);
    }
    // Seeded lists: 1–8 reactants, moles and coefficients log-uniform; near-ties (one ratio a few ulps from another).
    let mut rng = Rng::new(0x5f03);
    for _ in 0..1200 {
        let n = usize::try_from(rng.int(1, 8)).unwrap_or(1);
        let (lo, hi) = [(-3.0, 3.0), (-300.0, 300.0), (-320.0, 308.0)][usize::try_from(rng.int(0, 2)).unwrap_or(0)];
        let mo: Vec<f64> = (0..n).map(|_| rng.log_uniform(lo, hi)).collect();
        let co: Vec<f64> = (0..n).map(|_| rng.log_uniform(-1.0, 1.0).ceil()).collect();
        limiting(m, &mo, &co);
        let p = rng.log_uniform(-1.0, 1.0).ceil();
        yield_moles(m, &mo, &co, p);
        let pw = rng.log_uniform(lo, hi);
        yield_moles(m, &mo, &co, pw);
    }
    for _ in 0..400 {
        // ties: m₂/c₂ = m₁/c₁ exactly in the reals, rounded (k·m, k·c), and its neighbours
        let a = rng.log_uniform(-5.0, 5.0);
        let c = rng.int(1, 20);
        let k = rng.int(2, 20);
        #[allow(clippy::cast_precision_loss)]
        let (cf, kf) = (c as f64, k as f64);
        let am = a * kf;
        let cm = cf * kf;
        limiting(m, &[a, am], &[cf, cm]);
        limiting(m, &[am, a], &[cm, cf]);
        limiting(m, &[a, am.next_down()], &[cf, cm]);
        limiting(m, &[a, am.next_up()], &[cf, cm]);
        yield_moles(m, &[a, am], &[cf, cm], 2.0);
    }

    // percent_yield: each argument over the edges with the other typical, the Rust tests, over- and underflow, seeded.
    let py = pyield;
    for x in edges() {
        py(m, x, 10.0);
        py(m, 8.0, x);
        py(m, x, x);
    }
    for (a, t) in [
        (8.0, 10.0),
        (11.0, 10.0),
        (5.0, 0.0),
        (0.0, 10.0),
        (1.0, 3.0),
        (2.0, 3.0),
        (1.0, 7.0),
        (0.7, 1.0),
        (f64::MAX, 1.0),
        (f64::MAX, 2.0),
        (1e307, 1e-10),
        (f64::from_bits(1), 1.0),
        (f64::from_bits(1), 1e300),
        (1e-300, 1e10),
        (f64::MIN_POSITIVE, 3.0),
        (1.0, f64::from_bits(1)),
        (1.0, f64::MAX),
    ] {
        py(m, a, t);
    }
    let mut rng2 = Rng::new(0x5f04);
    for _ in 0..1500 {
        let t = rng2.log_uniform(-300.0, 300.0);
        let a = t * rng2.uniform(0.0, 1.2);
        py(m, a, t);
        let a2 = rng2.log_uniform(-320.0, 308.0);
        let t2 = rng2.log_uniform(-320.0, 308.0);
        py(m, a2, t2);
    }
    for k in 1..=1000 {
        // the integer percentages: k/1000 of a whole yield
        py(m, f64::from(k), 1000.0);
    }
}

// ── percent composition and the empirical formula ───────────────────────────────────

fn pairs_in(pairs: &[(u8, u32)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(pairs.len())];
    for &(z, c) in pairs {
        v.push(i(z));
        v.push(i(c));
    }
    v
}

fn composition(m: &mut Module, pairs: &[(u8, u32)]) {
    let mol = Molecule::new(pairs);
    let oracle = percent_composition(bb(&mol));
    let entries = |v: &[(u8, f64)]| {
        let mut t = vec![Tok::V(v.len())];
        for &(z, p) in v {
            t.push(i(z));
            t.push(f(p));
        }
        t
    };
    let outs = match &oracle {
        Ok(v) => entries(v),
        Err(e_) => vec![err(e_)],
    };
    // ADR 0002 #24: each entry RN(mass·(100·count)/mw) with mw the oracle's molecular_weight (a count of 0: +0, Rust's)
    if oracle.is_ok() {
        let mw = mol.molecular_weight().unwrap_or(f64::NAN);
        let port: Vec<(u8, f64)> = mol
            .atoms
            .iter()
            .map(|a| {
                let mass = lookup_by_number(a.element_number).map_or(f64::NAN, |e| e.atomic_mass);
                (a.element_number, if a.count == 0 { 0.0 } else { cr_mul_div(mass, 100.0 * f64::from(a.count), mw) })
            })
            .collect();
        let pt = entries(&port);
        if !same_toks(&pt, &outs) {
            let note = format!(
                "ADR 0002 #24: the port returns the correctly rounded mass·count·100/mw per entry {}",
                tokens_note(&pt)
            );
            m.row_noted("stoichiometry::percent_composition", Class::B, &pairs_in(pairs), &outs, &note);
            return;
        }
    }
    m.row("stoichiometry::percent_composition", Class::B, &pairs_in(pairs), &outs);
}

/// ADR 0002 #26: the port's empirical formula: Rust's validation, in its order, then `exact::empirical_exact` (exact
/// moles, minimum, ratios, tolerance 1/10 and counts). The Err is the variant's name.
fn port_empirical(elements: &[(u8, f64)]) -> Result<Vec<(u8, u32)>, &'static str> {
    if elements.is_empty() {
        return Err("InvalidInput");
    }
    let mut mp = Vec::with_capacity(elements.len());
    for &(z, p) in elements {
        let e = lookup_by_number(z).ok_or("InvalidElement")?;
        if !(p.is_finite() && p > 0.0) {
            return Err("InvalidInput");
        }
        mp.push((e.atomic_mass, p));
    }
    let counts = empirical_exact(&mp).ok_or("ComputationError")?;
    Ok(elements.iter().zip(counts).map(|(&(z, _), c)| (z, c)).collect())
}

fn empirical(m: &mut Module, elements: &[(u8, f64)]) {
    let counts = |v: &[(u8, u32)]| {
        let mut t = vec![Tok::V(v.len())];
        for &(z, c) in v {
            t.push(i(z));
            t.push(i(c));
        }
        t
    };
    let oracle = empirical_formula(bb(elements));
    let outs = match &oracle {
        Ok(v) => counts(v),
        Err(e_) => vec![err(e_)],
    };
    let mut ins = vec![Tok::V(elements.len())];
    for &(z, p) in elements {
        ins.push(i(z));
        ins.push(f(p));
    }
    let port = match port_empirical(elements) {
        Ok(v) => counts(&v),
        Err(name) => vec![Tok::Err(name.to_string())],
    };
    if let Err(KimiyaError::InvalidInput(_) | KimiyaError::InvalidElement(_)) = &oracle {
        assert!(same_toks(&port, &outs), "port_empirical's validation differs: {elements:?}");
    }
    if !same_toks(&port, &outs) {
        let note = format!(
            "ADR 0002 #26: with exact ratios and an exact 1/10 tolerance the port returns {}",
            tokens_note(&port)
        );
        m.row_noted("stoichiometry::empirical_formula", Class::A, &ins, &outs, &note);
        return;
    }
    m.row("stoichiometry::empirical_formula", Class::A, &ins, &outs);
}

/// The mass percentages of a formula, from the oracle's own `percent_composition` (Err: none).
fn percents(pairs: &[(u8, u32)]) -> Option<Vec<(u8, f64)>> {
    percent_composition(&Molecule::new(pairs)).ok()
}

fn emit_composition(m: &mut Module) {
    let mut formulas: Vec<Vec<(u8, u32)>> = vec![
        vec![(1, 2), (8, 1)],
        vec![(6, 1), (8, 2)],
        vec![(6, 1), (1, 4)],
        vec![(6, 6), (1, 12), (8, 6)],
        vec![],
        vec![(1, 0)],
        vec![(1, 0), (8, 0)],
        vec![(1, 0), (8, 1)],
        vec![(200, 1)],
        vec![(0, 1), (8, 1)],
        vec![(6, 1), (200, 0)],
        vec![(1, u32::MAX), (8, 1)],
        vec![(1, u32::MAX), (1, u32::MAX)],
        vec![(118, u32::MAX)],
        vec![(26, 2), (8, 3)],
        vec![(26, 3), (8, 4)],
        vec![(1, 1), (1, 1), (8, 1)],
        vec![(19, 1), (25, 1), (8, 4)],
        vec![(6, 2), (1, 6), (8, 1)],
        vec![(20, 3), (15, 2), (8, 8)],
    ];
    for z in 1..=118_u8 {
        formulas.push(vec![(z, 1)]);
        formulas.push(vec![(z, 7)]);
        formulas.push(vec![(z, 1), (8, 2)]);
        formulas.push(vec![(z, 2), (17, 3)]);
    }
    let mut rng = Rng::new(0x5f05);
    for _ in 0..800 {
        let n = usize::try_from(rng.int(1, 6)).unwrap_or(1);
        let f_: Vec<(u8, u32)> = (0..n)
            .map(|_| {
                let z = u8::try_from(rng.int(1, 118)).unwrap_or(1);
                let c = match rng.int(0, 20) {
                    0 => 0,
                    1 => u32::MAX,
                    _ => u32::try_from(rng.int(1, 30)).unwrap_or(1),
                };
                (z, c)
            })
            .collect();
        formulas.push(f_);
    }
    for pairs in &formulas {
        composition(m, pairs);
    }

    // empirical_formula: the Rust tests, in order.
    let o108 = 1.08 * 30.0 / 12.011 * 15.999;
    let o143 = 1.43 * 30.0 / 12.011 * 15.999;
    let tests: &[&[(u8, f64)]] = &[
        &[(1, 11.19), (8, 88.81)],
        &[(6, 40.00), (1, 6.71), (8, 53.29)],
        &[(26, 69.94), (8, 30.06)],
        &[(26, 72.36), (8, 27.64)],
        &[(6, 30.0), (8, o108)],
        &[(6, 30.0), (8, o143)],
        &[(6, f64::NAN), (8, 50.0)],
        &[(6, f64::INFINITY), (8, 50.0)],
        &[(6, 40.0), (8, -5.0)],
        &[(6, 1e-320), (8, 50.0)],
        &[],
    ];
    for t in tests {
        empirical(m, t);
    }
    // Each multiplier 1 … 6 and none: O/C ratios k/6 … over a fine ladder around each.
    for num in 0..=420_u32 {
        let r = 1.0 + f64::from(num) / 120.0;
        empirical(m, &[(6, 12.011), (8, r * 15.999)]);
    }
    // Real compounds: their exact formulas' percentages (from the oracle's percent_composition) give the formula back,
    // or its reduction, or an Err where a ratio is not within 0.1 at any k ≤ 6.
    for pairs in formulas.iter().filter(|p| !p.is_empty() && p.iter().all(|&(z, c)| (1..=118).contains(&z) && c > 0 && c <= 40)) {
        if let Some(pc) = percents(pairs) {
            empirical(m, &pc);
        }
    }
    for pairs in [
        vec![(6_u8, 1_u32), (1, 2), (8, 1)],
        vec![(26, 2), (8, 3)],
        vec![(26, 3), (8, 4)],
        vec![(6, 4), (1, 5)],
        vec![(6, 5), (1, 6)],
        vec![(6, 6), (1, 7)],
        vec![(6, 5), (1, 7)],
        vec![(6, 7), (1, 6)],
        vec![(15, 4), (8, 10)],
        vec![(13, 2), (16, 3), (8, 12)],
        vec![(20, 3), (15, 2), (8, 8)],
        vec![(6, 8), (1, 10), (7, 4), (8, 2)],
    ] {
        if let Some(pc) = percents(&pairs) {
            empirical(m, &pc);
            let rev: Vec<(u8, f64)> = pc.iter().rev().copied().collect();
            empirical(m, &rev);
        }
    }
    // The 0.1 tolerance from both sides: moles of C = 1 exactly (12.011 g), O/C = n ± 0.1 to a few ulps either side
    // (the deviation is a multiple of the ratio's ulp, so 0.1 lies between two such multiples).
    for n in [1.0_f64, 2.0, 3.0, 5.0] {
        for dev in [0.1_f64, -0.1] {
            let target = n + dev;
            let mut o = target * 15.999;
            for _ in 0..40 {
                o = o.next_down();
            }
            for _ in 0..80 {
                empirical(m, &[(6, 12.011), (8, o)]);
                o = o.next_up();
            }
        }
    }
    // Err paths and edges: invalid elements, each percentage edge, a subnormal or huge ratio, counts past u32::MAX.
    for x in edges() {
        empirical(m, &[(6, x), (8, 50.0)]);
        empirical(m, &[(6, 40.0), (8, x)]);
    }
    for z in [0_u8, 119, 200, 255] {
        empirical(m, &[(z, 10.0), (8, 50.0)]);
        empirical(m, &[(6, f64::NAN), (z, 50.0)]);
    }
    for (a, b) in [
        (1e-310, 50.0),
        (5e-324, 1.0),
        (1.0, 1e300),
        (1.0, 4.0e10),
        (1.0, 4.294_967_295e9 * 15.999 / 12.011),
        (12.011, 4_294_967_295.0 * 15.999),
        (12.011, 4_294_967_295.4 * 15.999),
        (12.011, 4_294_967_295.6 * 15.999),
        (12.011, 4_294_967_296.0 * 15.999),
        (f64::MAX, 1.0),
        (f64::MAX, f64::MAX),
        (1e-300, 1e-300),
    ] {
        empirical(m, &[(6, a), (8, b)]);
    }
    empirical(m, &[(6, 40.0)]);
    empirical(m, &[(6, 40.0), (6, 40.0)]);
    empirical(m, &[(6, 40.0), (6, 20.0)]);
    // Seeded formulas: random small counts, their percentages scaled by a random factor (only ratios matter).
    let mut rng3 = Rng::new(0x5f06);
    for _ in 0..1200 {
        let n = usize::try_from(rng3.int(1, 4)).unwrap_or(1);
        let pairs: Vec<(u8, u32)> =
            (0..n).map(|_| (u8::try_from(rng3.int(1, 92)).unwrap_or(1), u32::try_from(rng3.int(1, 12)).unwrap_or(1))).collect();
        if let Some(pc) = percents(&pairs) {
            let s = rng3.log_uniform(-3.0, 3.0);
            let scaled: Vec<(u8, f64)> = pc.iter().map(|&(z, p)| (z, p * s)).collect();
            empirical(m, &scaled);
        }
    }
}

/// ADR 0002 #24–#26's hard rows (M4f-2; appended, with their own seed, so the rows before keep theirs): percent_yield's
/// exact ties at RN — subnormal (q·2⁻¹⁰⁷⁴ over 200·2^j, q odd: (q/2)·2⁻¹⁰⁷⁴) and normal (m·2^e over 5·2^j with 5m a 54-bit
/// odd number; 3m·2^e over 3·2^j with 25m one) at scales across the range —, its exact values (over 25 and 100·2^j), the
/// subnormal-result and overflow ends to the ulp; theoretical_yield_moles where m/c over- or underflows though m·p/c is
/// finite, and its midpoints (m·3/1, m·3/2); limiting_reagent's exact ties, the near-ties the rounded ratios cannot tell
/// apart (both orders), both ratios overflowing or underflowing, 0 against an underflowed ratio; empirical_formula's
/// knife edges at ±1/10 — the same element twice with r = (10n ± 1)/(10k) exactly and one ulp either side, other elements
/// with r·k within a few ulps of n ± 1/10 —, the subnormal percentages whose moles Rust rounds to 0, and the count's
/// 2³² edge.
fn emit_exact(m: &mut Module) {
    let mut rng = Rng::new(0x5f24);
    let dmax = f64::MAX;
    let sub = f64::from_bits;
    // percent_yield
    pyield(m, sub(1), 10.0);
    for q in [1_u64, 3, 5, 7, 9, 99, (1 << 52) - 1, (1 << 52) - 3, 1 << 51] {
        for j in [0, 1, 2] {
            pyield(m, sub(q), 200.0 * pow2(j));
        }
        pyield(m, sub(q), 300.0);
        pyield(m, sub(q), 100.0);
    }
    for (lo, hi) in [((1_u64 << 53) / 5 + 1, (1_u64 << 54) / 5), ((1_u64 << 53) / 25 + 1, (1_u64 << 54) / 25)] {
        for c in [lo | 1, ((lo + hi) / 2) | 1, (hi - 1) | 1] {
            #[allow(clippy::cast_precision_loss)]
            let cf = c as f64; // c < 2⁵³: exact
            for e in [-1000, -500, -60, 0, 60, 500, 900] {
                if lo < (1 << 50) {
                    pyield(m, 3.0 * cf * pow2(e), 3.0 * pow2(e / 2)); // 25c a 54-bit odd number
                } else {
                    pyield(m, cf * pow2(e), 5.0 * pow2(e / 2)); // 5c a 54-bit odd number
                }
            }
        }
    }
    for _ in 0..200 {
        let a = rng.log_uniform(-300.0, 300.0);
        pyield(m, a, 25.0);
        pyield(m, a, 100.0 * pow2(i32::try_from(rng.int(-60, 60)).unwrap_or(0)));
    }
    for t in [1.0, 2.0, 0.5, 3.0, 7.0, 99.999, 1e-10, pow2(-60), pow2(60)] {
        let a0 = dmax / 100.0 * t; // near 100·a/t = DBL_MAX
        let mut a = a0;
        for _ in 0..3 {
            a = a.next_down();
        }
        for _ in 0..7 {
            if a.is_finite() {
                pyield(m, a, t);
            }
            a = a.next_up();
        }
    }
    for k in [1_u64, 2, 3, 4, 5, 1 << 20] {
        for j in [0, 10, 52] {
            pyield(m, sub(k) * pow2(j), 100.0 * pow2(j)); // exactly k·2⁻¹⁰⁷⁴
            pyield(m, sub(k) * pow2(j), 200.0 * pow2(j)); // (k/2)·2⁻¹⁰⁷⁴: a tie for odd k
            pyield(m, sub(k) * pow2(j), 300.0 * pow2(j));
        }
    }
    pyield(m, f64::MIN_POSITIVE, 100.0 * pow2(60));
    pyield(m, f64::MIN_POSITIVE, 3.0);
    // theoretical_yield_moles
    for (mo, c, p) in [
        (1e300, 1e-10, 1e-10),
        (dmax, 0.5, 0.25),
        (dmax, 0.5, 0.5),
        (1e-300, 1e10, 1e10),
        (sub(1), 2.0, 2.0),
        (sub(3), 2.0, 2.0),
        (sub(5), 4.0, 3.0),
        (1e-310, 1e10, 1e20),
        (1e308, 1e-300, 1e-300),
    ] {
        yield_moles(m, &[mo], &[c], p);
    }
    for _ in 0..40 {
        // m/c past DBL_MAX or below 2⁻¹⁰²² while m·p/c is a normal double
        let (mo, c) = (rng.log_uniform(300.0, 308.0), rng.log_uniform(-12.0, -1.0));
        yield_moles(m, &[mo], &[c], c * rng.uniform(0.01, 1.0));
        let (mu, cu) = (rng.log_uniform(-320.0, -300.0), rng.log_uniform(1.0, 12.0));
        yield_moles(m, &[mu], &[cu], cu * rng.uniform(0.5, 2.0));
    }
    for _ in 0..100 {
        #[allow(clippy::cast_precision_loss)]
        let y = (rng.int(((1 << 53) / 3 + 1) / 2, ((1 << 54) / 3 - 1) / 2) * 2 + 1) as f64; // 3y a 54-bit odd number
        let sc = pow2(i32::try_from(rng.int(-900, 900)).unwrap_or(0));
        yield_moles(m, &[y * sc], &[1.0], 3.0);
        yield_moles(m, &[y * sc], &[2.0], 3.0);
        yield_moles(m, &[y * sc, 1e300], &[2.0, 1.0], 3.0);
    }
    // limiting_reagent: exact ties, near-ties (both orders), over- and underflowing ratios, 0 against an underflow
    let lists: &[(&[f64], &[f64])] = &[
        (&[1.0, 2.0], &[1.0, 2.0]),
        (&[0.1, 0.2], &[1.0, 2.0]),
        (&[0.2, 0.1], &[2.0, 1.0]),
        (&[3.0, 6.0, 1.5], &[1.0, 2.0, 0.5]),
        (&[f64::MAX, f64::MAX], &[f64::from_bits(1), 1e-300]),
        (&[f64::MAX, f64::MAX], &[1e-300, f64::from_bits(1)]),
        (&[f64::MAX, 1e308], &[1e-300, 1e-300]),
        (&[1e308, f64::MAX], &[1e-300, 1e-300]),
        (&[f64::from_bits(2), f64::from_bits(1)], &[f64::MAX, f64::MAX]),
        (&[f64::from_bits(1), f64::from_bits(2)], &[f64::MAX, f64::MAX]),
        (&[f64::from_bits(1), 0.0], &[f64::MAX, 1.0]),
        (&[0.0, f64::from_bits(1)], &[1.0, f64::MAX]),
        (&[-0.0, 0.0], &[1.0, 1.0]),
        (&[9.051_597_531_392_973, 54.309_585_188_357_836], &[1.0, 6.0]),
        (&[54.309_585_188_357_836, 9.051_597_531_392_973], &[6.0, 1.0]),
    ];
    for &(mo, co) in lists {
        limiting(m, mo, co);
        yield_moles(m, mo, co, 2.0);
    }
    for _ in 0..400 {
        let a = rng.log_uniform(-300.0, 300.0);
        #[allow(clippy::cast_precision_loss)]
        let (c, k) = (rng.int(1, 30) as f64, rng.int(2, 30) as f64);
        let am = a * k;
        for b in [am.next_down(), am, am.next_up()] {
            if (a / c).to_bits() == (b / (c * k)).to_bits() {
                limiting(m, &[a, b], &[c, c * k]); // equal rounded ratios: the exact comparison decides
                limiting(m, &[b, a], &[c * k, c]);
                yield_moles(m, &[a, b], &[c, c * k], 3.0);
            }
        }
        let big = |r: &mut Rng| r.log_uniform(300.0, 308.0);
        let small = |r: &mut Rng| r.log_uniform(-320.0, -300.0);
        let (m1, m2, c1, c2) = (big(&mut rng), big(&mut rng), small(&mut rng), small(&mut rng));
        limiting(m, &[m1, m2], &[c1, c2]); // both ratios +∞ in Rust
        limiting(m, &[c1, c2], &[m1, m2]); // both 0 in Rust
    }
    // empirical_formula: knife edges at ±1/10
    let o_edge = 17.5989; // C 12.011 %, O 17.5989 %: r − 1 = 1/10 − 1/(10·9006636304787571)
    let mut o = o_edge;
    for _ in 0..3 {
        o = f64::next_down(o);
    }
    for _ in 0..7 {
        empirical(m, &[(6, 12.011), (8, o)]);
        o = o.next_up();
    }
    for z in [6_u8, 8, 26, 92] {
        for k in 1..=6_u32 {
            for n in [k, k + 1, 2 * k + 1, 5 * k] {
                for sgn in [1_i64, -1] {
                    let num = 10 * i64::from(n) + sgn; // r = (10n ± 1)/(10k), so r·k − n = ±1/10
                    #[allow(clippy::cast_precision_loss)]
                    let (pa, pb) = (f64::from(10 * k), num as f64);
                    for w in [1.0, pow2(-1000), pow2(900)] {
                        let (a, b) = (pa * w, pb * w);
                        empirical(m, &[(z, a), (z, b)]);
                        empirical(m, &[(z, a), (z, b.next_up())]);
                        empirical(m, &[(z, a), (z, b.next_down())]);
                    }
                }
            }
        }
    }
    let mass = |z: u8| lookup_by_number(z).map_or(1.0, |e| e.atomic_mass);
    for _ in 0..300 {
        let (zi, zj) = (u8::try_from(rng.int(1, 118)).unwrap_or(1), u8::try_from(rng.int(1, 118)).unwrap_or(1));
        let k = rng.int(1, 6);
        let n = rng.int(k, 6 * k);
        let sgn = if rng.int(0, 1) == 0 { 1 } else { -1 };
        #[allow(clippy::cast_precision_loss)]
        let target = (10 * n + sgn) as f64 / (10 * k) as f64;
        let pj = rng.log_uniform(-2.0, 2.0);
        let pi = target * mass(zi) / mass(zj) * pj;
        let mut x = pi;
        for _ in 0..2 {
            x = x.next_down();
        }
        for _ in 0..5 {
            empirical(m, &[(zj, pj), (zi, x)]);
            x = x.next_up();
        }
    }
    // subnormal percentages whose moles Rust rounds to 0, or rounds coarsely
    for (a, b) in [(sub(1), sub(2)), (sub(1), sub(1)), (sub(2), sub(4)), (sub(3), sub(6)), (sub(1), sub(3))] {
        empirical(m, &[(2, a), (2, b)]);
        empirical(m, &[(1, a), (1, b)]);
        empirical(m, &[(118, a), (118, b)]);
    }
    for e in [-1074, -1070, -1060, -1050, -1040, -1023] {
        let s = pow2(e);
        empirical(m, &[(6, 12.011 * s), (8, 15.999 * s)]);
        empirical(m, &[(6, 12.011 * s), (8, 2.0 * 15.999 * s)]);
        empirical(m, &[(1, 2.0 * 1.008 * s), (8, 15.999 * s)]);
    }
    // the count's 2³² edge (one element twice: the ratio is the percentages' quotient)
    for b in [4_294_967_295.0, 4_294_967_295.25, 4_294_967_295.5, 4_294_967_295.75, 4_294_967_294.5, 4_294_967_296.0] {
        empirical(m, &[(6, 1.0), (6, b)]);
        empirical(m, &[(6, 1.0), (6, f64::next_up(b))]);
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("stoichiometry");
    emit_balance(&mut m);
    emit_oxidation(&mut m);
    emit_yields(&mut m);
    emit_composition(&mut m);
    emit_exact(&mut m);
    m
}
