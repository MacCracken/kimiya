//! `network.tsv` — reaction network analysis (port plan §4, M5b): stoichiometric matrices, their CSR form, the rank and
//! the conservation-law count, species rates of change.
//!
//! - Class A: `stoichiometric_matrix` (input `i:n_species v:n_reactions`, then per reaction its reactants `v:k` + k
//!   pairs `i:species f:coefficient` and its products the same; output the dense matrix as Rust's `Vec<Vec<f64>>`,
//!   `v:rows` then each row `v:len` + `f:` entries, or `E:`) and `sparse_stoichiometric_matrix` (the same input; output
//!   hisab 1.4.0's `CsrMatrix`: `i:nrows i:ncols v:nnz`, each stored entry `i:col f:value` row by row, then the row
//!   offsets `v:nrows+1` + `i:`): the mechanism presets below (the Rust tests', the reaction_dynamics laws, Chapman's
//!   ozone cycle, the H₂ + Br₂ chain, Robertson, Lotka–Volterra, the Brusselator and the Oregonator, H₂/O₂ combustion,
//!   Michaelis–Menten, a catalyst, sources and sinks), repeated species (summed), coefficients on hisab's 10⁻¹² CSR
//!   threshold to the ulp, ±0, NaN and ±∞ coefficients, every Err (no species, no reactions, a species index at and past
//!   n_species in a reactant or a product of any reaction), seeded networks. `reaction_network_rank` (input the matrix
//!   as Rust's `&[Vec<f64>]`; output `i:rank` or `E:ComputationError`): the presets' matrices and their transposes,
//!   scaled across the 10⁻¹⁰ singular-value threshold and to the over- and underflow of the column sums of squares,
//!   rank-deficient ones (duplicated, reversed and combined reactions; seeded exact products U·V of rank r), noisy
//!   ones, Hilbert matrices, the empty matrix and rows of length 0, NaN and ±∞ entries, single rows and columns, seeded
//!   uniform and integer matrices. hisab 1.4.0's one-sided Jacobi SVD tests its rotations against EPSILON_F64² =
//!   10⁻²⁴, which rounding rarely reaches, so most dense matrices (Chapman's, Robertson's, the H₂ + Br₂ chain's, about
//!   six in ten seeded ones) end in NoConvergence after 100·max(m, n) sweeps: the oracle's ComputationError, which the
//!   port's transcription (`src/numerics.cyr`, ADR 0003's M5b amendment) reproduces. `conservation_law_count` (input
//!   `i:n_species i:rank`, usize words; output `i:`): the edges and every saturating case.
//! - Class B: `rate_from_stoichiometry` (the matrix, then the rates `v:n` + `f:`; output `v:rows` + `f:`): the presets'
//!   matrices at seeded rates, the −0.0 sum of an empty row, ±0 terms, overflow and ∞ − ∞, NaN rates, every Err
//!   (an empty matrix, a rate vector longer or shorter than the rows). Rust's ragged matrices (a row of another length)
//!   have no flat counterpart in the port (port plan §6.3) and are not emitted.
//!
//! Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::grids::{Rng, pow2};
use crate::tsv::{Class, Module, Tok, err, f, i};
use hisab::CsrMatrix;
use kimiya::network::{
    ReactionDef, conservation_law_count, rate_from_stoichiometry, reaction_network_rank, sparse_stoichiometric_matrix,
    stoichiometric_matrix,
};

/// The harness's token list holds 8448 tokens (`tests/parity_core.cyr`): every row must fit.
const MAX_TOKENS: usize = 8448;

/// A usize as its i64 bits, the Cyrius word.
fn word(x: usize) -> Tok {
    i(i64::from_ne_bytes((x as u64).to_ne_bytes()))
}

fn terms_toks(v: &mut Vec<Tok>, terms: &[(usize, f64)]) {
    v.push(Tok::V(terms.len()));
    for &(s, c) in terms {
        v.push(word(s));
        v.push(f(c));
    }
}

/// `i:n_species v:n_reactions`, then each reaction's reactants and products.
fn net_toks(n: usize, rx: &[ReactionDef]) -> Vec<Tok> {
    let mut v = vec![word(n), Tok::V(rx.len())];
    for (r, p) in rx {
        terms_toks(&mut v, r);
        terms_toks(&mut v, p);
    }
    v
}

/// Rust's `&[Vec<f64>]`: `v:rows`, then each row `v:len` and its entries.
pub fn matrix_toks(mat: &[Vec<f64>]) -> Vec<Tok> {
    let mut v = vec![Tok::V(mat.len())];
    for row in mat {
        v.push(Tok::V(row.len()));
        v.extend(row.iter().map(|&x| f(x)));
    }
    v
}

fn vec_toks(x: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(x.len())];
    v.extend(x.iter().map(|&y| f(y)));
    v
}

/// hisab 1.4.0's `CsrMatrix` as tokens. Its arrays are private, so they are read back through `to_dense` (an entry is
/// stored exactly where it is nonzero there: `from_dense` keeps |v| > 10⁻¹², never a zero or a NaN) and proven to be the
/// matrix's own by rebuilding it with `CsrMatrix::new` and comparing (`PartialEq` compares every field).
fn csr_toks(csr: &CsrMatrix) -> Vec<Tok> {
    let dense = csr.to_dense();
    let (mut values, mut cols, mut offsets) = (Vec::new(), Vec::new(), vec![0_usize]);
    for row in &dense {
        for (j, &x) in row.iter().enumerate() {
            if x != 0.0 {
                values.push(x);
                cols.push(j);
            }
        }
        offsets.push(values.len());
    }
    assert_eq!(values.len(), csr.nnz(), "CSR: the stored entries are not to_dense's nonzero ones");
    let rebuilt = CsrMatrix::new(csr.nrows, csr.ncols, values.clone(), cols.clone(), offsets.clone());
    assert!(matches!(rebuilt, Ok(ref r) if r == csr), "CSR: the rebuilt matrix differs from the oracle's");
    let mut v = vec![word(csr.nrows), word(csr.ncols), Tok::V(values.len())];
    for (&c, &x) in cols.iter().zip(&values) {
        v.push(word(c));
        v.push(f(x));
    }
    v.push(Tok::V(offsets.len()));
    v.extend(offsets.iter().map(|&o| word(o)));
    v
}

fn row(m: &mut Module, func: &str, class: Class, ins: &[Tok], outs: &[Tok]) {
    assert!(ins.len() <= MAX_TOKENS && outs.len() <= MAX_TOKENS, "{func}: a row of {} / {} tokens", ins.len(), outs.len());
    m.row(func, class, ins, outs);
}

/// One network through both matrix builders.
fn builders(m: &mut Module, n: usize, rx: &[ReactionDef]) {
    let ins = net_toks(n, rx);
    let outs = match stoichiometric_matrix(bb(n), bb(rx)) {
        Ok(mat) => matrix_toks(&mat),
        Err(e) => vec![err(&e)],
    };
    row(m, "network::stoichiometric_matrix", Class::A, &ins, &outs);
    let outs = match sparse_stoichiometric_matrix(bb(n), bb(rx)) {
        Ok(csr) => csr_toks(&csr),
        Err(e) => vec![err(&e)],
    };
    row(m, "network::sparse_stoichiometric_matrix", Class::A, &ins, &outs);
}

fn rank(m: &mut Module, mat: &[Vec<f64>]) {
    let outs = match reaction_network_rank(bb(mat)) {
        Ok(r) => vec![word(r)],
        Err(e) => vec![err(&e)],
    };
    row(m, "network::reaction_network_rank", Class::A, &matrix_toks(mat), &outs);
}

fn rate(m: &mut Module, mat: &[Vec<f64>], rates: &[f64]) {
    let mut ins = matrix_toks(mat);
    ins.extend(vec_toks(rates));
    let outs = match rate_from_stoichiometry(bb(mat), bb(rates)) {
        Ok(dx) => vec_toks(&dx),
        Err(e) => vec![err(&e)],
    };
    row(m, "network::rate_from_stoichiometry", Class::B, &ins, &outs);
}

fn claws(m: &mut Module, n: usize, r: usize) {
    let c = conservation_law_count(bb(n), bb(r));
    m.row("network::conservation_law_count", Class::A, &[word(n), word(r)], &[word(c)]);
}

fn t(terms: &[(usize, f64)]) -> Vec<(usize, f64)> {
    terms.to_vec()
}

/// The mechanism presets: (name, species count, reactions). Species are numbered in the comment's order.
#[allow(clippy::too_many_lines)] // one mechanism per entry
fn presets() -> Vec<(&'static str, usize, Vec<ReactionDef>)> {
    vec![
        // The Rust tests: A → B; A → B → C.
        ("first_order", 2, vec![(t(&[(0, 1.0)]), t(&[(1, 1.0)]))]),
        ("consecutive", 3, vec![(t(&[(0, 1.0)]), t(&[(1, 1.0)])), (t(&[(1, 1.0)]), t(&[(2, 1.0)]))]),
        // A ⇌ B (rank 1: the reverse reaction is the forward one negated).
        ("reversible", 2, vec![(t(&[(0, 1.0)]), t(&[(1, 1.0)])), (t(&[(1, 1.0)]), t(&[(0, 1.0)]))]),
        // 2A → B (second order).
        ("second_order", 2, vec![(t(&[(0, 2.0)]), t(&[(1, 1.0)]))]),
        // E, S, ES, P: E + S ⇌ ES → E + P (rank 2: E and S + ES + P conserved).
        ("michaelis_menten", 4, vec![
            (t(&[(0, 1.0), (1, 1.0)]), t(&[(2, 1.0)])),
            (t(&[(2, 1.0)]), t(&[(0, 1.0), (1, 1.0)])),
            (t(&[(2, 1.0)]), t(&[(0, 1.0), (3, 1.0)])),
        ]),
        // O, O₂, O₃: Chapman's ozone cycle (rank 2: oxygen atoms conserved).
        ("chapman", 3, vec![
            (t(&[(1, 1.0)]), t(&[(0, 2.0)])),
            (t(&[(0, 1.0), (1, 1.0)]), t(&[(2, 1.0)])),
            (t(&[(2, 1.0)]), t(&[(0, 1.0), (1, 1.0)])),
            (t(&[(0, 1.0), (2, 1.0)]), t(&[(1, 2.0)])),
        ]),
        // H₂, Br₂, HBr, H, Br: the H₂ + Br₂ chain (rank 3: H and Br conserved).
        ("h2_br2", 5, vec![
            (t(&[(1, 1.0)]), t(&[(4, 2.0)])),
            (t(&[(4, 1.0), (0, 1.0)]), t(&[(2, 1.0), (3, 1.0)])),
            (t(&[(3, 1.0), (1, 1.0)]), t(&[(2, 1.0), (4, 1.0)])),
            (t(&[(3, 1.0), (2, 1.0)]), t(&[(0, 1.0), (4, 1.0)])),
            (t(&[(4, 2.0)]), t(&[(1, 1.0)])),
        ]),
        // A, B, C: Robertson's stiff kinetics, A → B, 2B → B + C, B + C → A + C (rank 2).
        ("robertson", 3, vec![
            (t(&[(0, 1.0)]), t(&[(1, 1.0)])),
            (t(&[(1, 2.0)]), t(&[(1, 1.0), (2, 1.0)])),
            (t(&[(1, 1.0), (2, 1.0)]), t(&[(0, 1.0), (2, 1.0)])),
        ]),
        // X, Y: Lotka–Volterra, X → 2X, X + Y → 2Y, Y → ∅ (a sink).
        ("lotka_volterra", 2, vec![
            (t(&[(0, 1.0)]), t(&[(0, 2.0)])),
            (t(&[(0, 1.0), (1, 1.0)]), t(&[(1, 2.0)])),
            (t(&[(1, 1.0)]), vec![]),
        ]),
        // X, Y: the Brusselator with A, B held constant, ∅ → X, 2X + Y → 3X, X → Y, X → ∅.
        ("brusselator", 2, vec![
            (vec![], t(&[(0, 1.0)])),
            (t(&[(0, 2.0), (1, 1.0)]), t(&[(0, 3.0)])),
            (t(&[(0, 1.0)]), t(&[(1, 1.0)])),
            (t(&[(0, 1.0)]), vec![]),
        ]),
        // X, Y, Z, P: the Oregonator with A, B held constant, f = 1/2: Y → X + P, X + Y → 2P, X → 2X + 2Z, 2X → P,
        // Z → ½Y (a fractional coefficient).
        ("oregonator", 4, vec![
            (t(&[(1, 1.0)]), t(&[(0, 1.0), (3, 1.0)])),
            (t(&[(0, 1.0), (1, 1.0)]), t(&[(3, 2.0)])),
            (t(&[(0, 1.0)]), t(&[(0, 2.0), (2, 2.0)])),
            (t(&[(0, 2.0)]), t(&[(3, 1.0)])),
            (t(&[(2, 1.0)]), t(&[(1, 0.5)])),
        ]),
        // H₂, O₂, OH, H₂O, H, O: H₂/O₂ combustion (rank 4: H and O conserved).
        ("h2_o2", 6, vec![
            (t(&[(0, 1.0), (1, 1.0)]), t(&[(2, 2.0)])),
            (t(&[(2, 1.0), (0, 1.0)]), t(&[(3, 1.0), (4, 1.0)])),
            (t(&[(4, 1.0), (1, 1.0)]), t(&[(2, 1.0), (5, 1.0)])),
            (t(&[(5, 1.0), (0, 1.0)]), t(&[(2, 1.0), (4, 1.0)])),
            (t(&[(4, 1.0), (4, 1.0)]), t(&[(0, 1.0)])),
            (t(&[(2, 2.0)]), t(&[(3, 1.0), (5, 1.0)])),
        ]),
        // A, B, E: a catalyst on both sides (its entry cancels to +0), A + E → B + E.
        ("catalysis", 3, vec![(t(&[(0, 1.0), (2, 1.0)]), t(&[(1, 1.0), (2, 1.0)]))]),
        // A, B: A + B → 2B, autocatalytic.
        ("autocatalytic", 2, vec![(t(&[(0, 1.0), (1, 1.0)]), t(&[(1, 2.0)]))]),
        // A, B, C: A → B, B → C, A → C (rank 2: the third is the sum of the first two).
        ("combined", 3, vec![
            (t(&[(0, 1.0)]), t(&[(1, 1.0)])),
            (t(&[(1, 1.0)]), t(&[(2, 1.0)])),
            (t(&[(0, 1.0)]), t(&[(2, 1.0)])),
        ]),
        // A → B twice (rank 1).
        ("duplicate", 2, vec![(t(&[(0, 1.0)]), t(&[(1, 1.0)])), (t(&[(0, 1.0)]), t(&[(1, 1.0)]))]),
        // A: ∅ → A, A → ∅.
        ("source_sink", 1, vec![(vec![], t(&[(0, 1.0)])), (t(&[(0, 1.0)]), vec![])]),
        // H₂, O₂, H₂O: H₂ + ½O₂ → H₂O, written with a repeated reactant (summed).
        ("water", 3, vec![(t(&[(0, 1.0), (1, 0.25), (1, 0.25)]), t(&[(2, 1.0)]))]),
        // An empty reaction (no reactants, no products) among others.
        ("empty_reaction", 2, vec![(vec![], vec![]), (t(&[(0, 1.0)]), t(&[(1, 1.0)]))]),
    ]
}

fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.first().map_or(0, Vec::len);
    (0..n).map(|j| a.iter().map(|r| r[j]).collect()).collect()
}

fn scaled(a: &[Vec<f64>], s: f64) -> Vec<Vec<f64>> {
    a.iter().map(|r| r.iter().map(|&x| x * s).collect()).collect()
}

fn mat(rows: &[&[f64]]) -> Vec<Vec<f64>> {
    rows.iter().map(|r| r.to_vec()).collect()
}

/// The n × n Hilbert matrix.
fn hilbert(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|r| (0..n).map(|c| 1.0 / f64::from(u32::try_from(r + c + 1).unwrap_or(1))).collect())
        .collect()
}

fn usz(x: i64) -> usize {
    usize::try_from(x).unwrap_or(0)
}

#[allow(clippy::too_many_lines)] // one function after another
pub fn emit() -> Module {
    let mut m = Module::new("network");
    let mut rng = Rng::new(0x5b01);
    let eps = 1e-12_f64; // hisab 1.4.0's EPSILON_F64, CsrMatrix::from_dense's threshold
    let thr = 1e-10_f64; // reaction_network_rank's singular-value threshold

    // ── stoichiometric_matrix, sparse_stoichiometric_matrix ──
    let presets = presets();
    for (_, n, rx) in &presets {
        builders(&mut m, *n, rx);
    }
    // The errors, in Rust's order: no species, no reactions, then a species index at or past n_species, reaction by
    // reaction, reactants before products.
    builders(&mut m, 0, &[(t(&[(0, 1.0)]), t(&[(1, 1.0)]))]);
    builders(&mut m, 0, &[]);
    builders(&mut m, 2, &[]);
    builders(&mut m, 2, &[(t(&[(2, 1.0)]), t(&[(1, 1.0)]))]);
    builders(&mut m, 2, &[(t(&[(0, 1.0)]), t(&[(2, 1.0)]))]);
    builders(&mut m, 2, &[(t(&[(5, 1.0)]), t(&[(1, 1.0)]))]); // the Rust test
    builders(&mut m, 2, &[(t(&[(usize::MAX, 1.0)]), vec![])]);
    builders(&mut m, 3, &[(t(&[(0, 1.0)]), t(&[(1, 1.0)])), (t(&[(1, 1.0)]), t(&[(3, 1.0)]))]);
    builders(&mut m, 3, &[(t(&[(0, 1.0)]), t(&[(1, 1.0)])), (t(&[(1, 1.0), (7, 2.0)]), t(&[(2, 1.0)]))]);
    builders(&mut m, 3, &[(t(&[(0, 1.0)]), t(&[(9, 1.0)])), (t(&[(8, 1.0)]), t(&[(2, 1.0)]))]);
    builders(&mut m, 1, &[(vec![], vec![])]);
    // Coefficients: hisab's CSR threshold to the ulp (|v| > 10⁻¹² is stored), ±0, the specials, extremes; repeated
    // species summed in list order.
    for c in [
        eps,
        eps.next_up(),
        eps.next_down(),
        -eps,
        -eps.next_up(),
        2e-12,
        1e-13,
        0.0,
        -0.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
        pow2(-1074),
        0.5,
        1.5,
        1e300,
    ] {
        builders(&mut m, 2, &[(t(&[(0, c)]), t(&[(1, c)]))]);
        builders(&mut m, 3, &[(t(&[(0, 1.0), (0, c)]), t(&[(1, 1.0), (2, c)])), (t(&[(2, c)]), t(&[(0, c)]))]);
    }
    builders(&mut m, 2, &[(t(&[(0, 1.0)]), t(&[(0, 1.0)]))]); // a species on both sides: 0 − 1 + 1 = +0
    builders(&mut m, 2, &[(t(&[(0, 1.0), (0, -1.0)]), t(&[(1, -0.0)]))]);
    builders(&mut m, 2, &[(t(&[(0, f64::MAX), (0, f64::MAX)]), t(&[(1, f64::INFINITY), (1, f64::NEG_INFINITY)]))]);
    builders(&mut m, 2, &[(t(&[(0, 1e-12), (0, 1e-12)]), t(&[(1, 5e-13), (1, 5e-13), (1, 5e-13)]))]);
    for _ in 0..120 {
        let n = usz(rng.int(1, 12));
        let k = usz(rng.int(1, 14));
        let rx: Vec<ReactionDef> = (0..k)
            .map(|_| {
                let side = |rng: &mut Rng| -> Vec<(usize, f64)> {
                    let len = usz(rng.int(0, 3));
                    (0..len).map(|_| (usz(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)), f64::from(u32::try_from(rng.int(1, 4)).unwrap_or(1)))).collect()
                };
                (side(&mut rng), side(&mut rng))
            })
            .collect();
        builders(&mut m, n, &rx);
    }
    for _ in 0..20 {
        let n = usz(rng.int(1, 6));
        let k = usz(rng.int(1, 6));
        let rx: Vec<ReactionDef> = (0..k)
            .map(|_| {
                let r = vec![(usz(rng.int(0, i64::try_from(n).unwrap_or(1))), rng.signed_log_uniform(-14.0, 3.0))];
                let p = vec![(usz(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)), rng.signed_log_uniform(-14.0, 3.0))];
                (r, p)
            })
            .collect();
        builders(&mut m, n, &rx); // some index n: an Err
    }

    // ── reaction_network_rank ──
    let mut mats = Vec::new();
    for (_, n, rx) in &presets {
        if let Ok(a) = stoichiometric_matrix(*n, rx) {
            mats.push(a);
        }
    }
    for a in &mats {
        rank(&mut m, a);
        rank(&mut m, &transpose(a));
        for s in [1e-11, 1e-10, 1e-9, 1e100, 1e154, 1e155, 1e160, 1e200, 1e-150, 1e-154, 1e-160, 1e-200, -1.0] {
            rank(&mut m, &scaled(a, s));
        }
    }
    // The threshold to the ulp: [[s]] has the singular value |s|.
    for s in [thr, thr.next_up(), thr.next_down(), -thr, -thr.next_up(), 0.0, -0.0, 1.0, f64::MIN_POSITIVE, pow2(-1074)] {
        rank(&mut m, &[vec![s]]);
        rank(&mut m, &[vec![s, 0.0], vec![0.0, 1.0]]);
        rank(&mut m, &[vec![1.0, 0.0, 0.0], vec![0.0, 2.0, 0.0], vec![0.0, 0.0, s]]);
    }
    for s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX, -f64::MAX, 1e308, 1e155, 1e-155, 1e-170] {
        rank(&mut m, &[vec![s]]);
        rank(&mut m, &[vec![s], vec![0.0], vec![1.0]]);
        rank(&mut m, &[vec![s, 0.0, 1.0]]);
        rank(&mut m, &[vec![s, 0.0], vec![0.0, 1.0]]);
        rank(&mut m, &[vec![1.0, 2.0], vec![3.0, s]]);
        rank(&mut m, &[vec![s, s], vec![s, s]]);
    }
    rank(&mut m, &[]); // the empty matrix
    rank(&mut m, &[vec![]]);
    rank(&mut m, &[vec![], vec![]]);
    rank(&mut m, &mat(&[&[0.0, 0.0], &[0.0, 0.0]]));
    rank(&mut m, &mat(&[&[-0.0, -0.0, -0.0]]));
    rank(&mut m, &mat(&[&[1.0, 2.0], &[2.0, 4.0]]));
    rank(&mut m, &mat(&[&[1.0, 2.0], &[2.0, 4.000000000000001]]));
    rank(&mut m, &mat(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], &[7.0, 8.0, 9.0]]));
    rank(&mut m, &mat(&[&[1.0, 0.0], &[0.0, 1.0]]));
    rank(&mut m, &mat(&[&[3.0, 4.0]]));
    rank(&mut m, &mat(&[&[3.0], &[4.0]]));
    rank(&mut m, &mat(&[&[1.0, 1.0], &[1.0, -1.0]]));
    rank(&mut m, &mat(&[&[2.0, 0.0, 0.0], &[0.0, 0.0, 3.0]]));
    for n in 2..=12 {
        rank(&mut m, &hilbert(n));
    }
    for _ in 0..4 {
        let row_: Vec<f64> = (0..6).map(|_| rng.uniform(-10.0, 10.0)).collect();
        rank(&mut m, std::slice::from_ref(&row_));
        rank(&mut m, &transpose(&[row_]));
    }
    // Seeded matrices, small ones of each shape: uniform reals, small integers (stoichiometry-like), exact rank-r
    // products U·V of integer factors, and the same with noise across the threshold.
    for rows in 1..=6 {
        for cols in 1..=6 {
            let a: Vec<Vec<f64>> = (0..rows).map(|_| (0..cols).map(|_| rng.uniform(-10.0, 10.0)).collect()).collect();
            rank(&mut m, &a);
            let b: Vec<Vec<f64>> =
                (0..rows).map(|_| (0..cols).map(|_| f64::from(i32::try_from(rng.int(-3, 3)).unwrap_or(0))).collect()).collect();
            rank(&mut m, &b);
        }
    }
    let int = |rng: &mut Rng, lo: i64, hi: i64| f64::from(i32::try_from(rng.int(lo, hi)).unwrap_or(0));
    for _ in 0..60 {
        let rows = usz(rng.int(1, 7));
        let cols = usz(rng.int(1, 7));
        let r = usz(rng.int(1, 3));
        let u: Vec<Vec<f64>> = (0..rows).map(|_| (0..r).map(|_| int(&mut rng, -2, 2)).collect()).collect();
        let v: Vec<Vec<f64>> = (0..r).map(|_| (0..cols).map(|_| int(&mut rng, -2, 2)).collect()).collect();
        let a: Vec<Vec<f64>> = (0..rows)
            .map(|i_| (0..cols).map(|j| (0..r).map(|k| u[i_][k] * v[k][j]).sum::<f64>()).collect())
            .collect();
        rank(&mut m, &a);
    }
    for noise in [1e-14, 1e-12, 1e-11, 1e-10, 1e-9, 1e-6] {
        for _ in 0..4 {
            let base = [[1.0, -1.0, 0.0], [0.0, 1.0, -1.0], [1.0, 0.0, -1.0]]; // rank 2
            let a: Vec<Vec<f64>> = base.iter().map(|r| r.iter().map(|&x| x + noise * rng.uniform(-1.0, 1.0)).collect()).collect();
            rank(&mut m, &a);
        }
    }
    // Larger sparse stoichiometric matrices: chains A₀ → A₁ → … (rank n − 1), which converge.
    for n in [5_usize, 10, 16] {
        let rx: Vec<ReactionDef> = (0..n - 1).map(|k| (vec![(k, 1.0)], vec![(k + 1, 1.0)])).collect();
        let a = stoichiometric_matrix(n, &rx).unwrap_or_default();
        rank(&mut m, &a);
        rank(&mut m, &transpose(&a));
    }

    // ── conservation_law_count ──
    for n in [0, 1, 2, 3, 5, 10, 1 << 32, usize::MAX - 1, usize::MAX, 1 << 63] {
        for r in [0, 1, 2, 3, 5, 10, 1 << 32, usize::MAX - 1, usize::MAX, 1 << 63] {
            claws(&mut m, n, r);
        }
    }

    // ── rate_from_stoichiometry ──
    for a in &mats {
        let k = a.first().map_or(0, Vec::len);
        for _ in 0..3 {
            let rates: Vec<f64> = (0..k).map(|_| rng.log_uniform(-6.0, 2.0)).collect();
            rate(&mut m, a, &rates);
        }
        rate(&mut m, a, &vec![0.5; k]);
        rate(&mut m, a, &vec![0.0; k]);
        rate(&mut m, a, &vec![-0.0; k]);
        rate(&mut m, a, &vec![1.0; k + 1]); // longer than a row
        if k > 0 {
            rate(&mut m, a, &vec![1.0; k - 1]); // shorter
        }
    }
    let ab = mat(&[&[-1.0], &[1.0]]);
    for v in [0.5, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -0.0, f64::MAX, pow2(-1074)] {
        rate(&mut m, &ab, &[v]);
    }
    rate(&mut m, &ab, &[0.5, 0.3]); // the Rust test
    rate(&mut m, &[], &[]);
    rate(&mut m, &[], &[1.0]);
    rate(&mut m, &[vec![]], &[]);
    rate(&mut m, &[vec![], vec![]], &[]);
    rate(&mut m, &[vec![], vec![]], &[1.0]);
    rate(&mut m, &mat(&[&[-0.0]]), &[1.0]);
    rate(&mut m, &mat(&[&[-0.0, -0.0]]), &[1.0, 2.0]);
    rate(&mut m, &mat(&[&[-1.0, 1.0]]), &[1.0, 1.0]);
    rate(&mut m, &mat(&[&[1e308, 1e308]]), &[10.0, 10.0]);
    rate(&mut m, &mat(&[&[1e308, -1e308]]), &[10.0, 10.0]);
    rate(&mut m, &mat(&[&[f64::INFINITY, 1.0]]), &[0.0, 1.0]);
    rate(&mut m, &mat(&[&[1.0, 1e-16, -1.0]]), &[1.0, 1.0, 1.0]);
    for _ in 0..40 {
        let rows = usz(rng.int(1, 8));
        let cols = usz(rng.int(1, 8));
        let a: Vec<Vec<f64>> = (0..rows).map(|_| (0..cols).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect()).collect();
        let rates: Vec<f64> = (0..cols).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect();
        rate(&mut m, &a, &rates);
    }
    m
}
