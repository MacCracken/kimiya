//! `numerics.tsv` — hisab 1.4.0's routines that `src/numerics.cyr` transcribes (ADR 0003, port plan §6.2; M5a),
//! called directly, so the edges kimiya's own guards never pass on are checked too. Replayed by `tests/numerics.tcyr`
//! (`tests/parity.tcyr` skips the module).
//!
//! - `num::Pcg32::next_u32` (A): `Pcg32::new(seed, 1)` (stream 1, as `simulate_stochastic`), then `count` draws: input
//!   `i:seed i:count` (a u64 seed as its i64 bits), output `v:count` + `i:` values. hisab 3.4.1's `pcg32_new` /
//!   `pcg32_next` must give these bits.
//! - `num::Pcg32::next_normal` (C: Box–Muller's ln and cos): the same seeds, `v:count` + `f:` values.
//! - `num::rk4_trajectory` (B): `s:system v:n f:params… f:t0 v:n f:y0… f:t_end i:steps` on the systems of
//!   `systems.rs`: a start time other than 0, backward integration (t_end < t0), t_end = t0, the empty state, an
//!   explicitly time-dependent rate, a derivative never written or NaN, steps 0 (ZeroSteps), 1, 2, 3, 10, 100, 1000.
//! - `num::dopri45` (C: the step factor's powf): `s:system v:n f:params… f:t0 v:n f:y0… f:t_end f:tol f:h_init`:
//!   h_init 0, negative, NaN, tiny, larger than the span; the empty state; backward integration; t_end = t0; a NaN
//!   derivative (the 100 000-step limit); tolerances 10⁻³ … 10⁻¹²; extreme scales.
//! - `num::gaussian_elimination` (B): `v:n` then n rows of `v:n+1 f:…`; output the solution `v:n f:x…`, then the
//!   eliminated matrix `v:n(n+1) f:…` (row-major), or `E:Variant`: the empty matrix, pivots exactly at, just above and
//!   just below hisab's EPSILON_F64 (10⁻¹²), zero diagonals that need a swap, ties in the pivot search, NaN and ±∞
//!   entries, singular, Hilbert and seeded random matrices of size 1 … 12.
//! - `num::euler_maruyama` (C): `s:drift v:n f:… s:diffusion v:n f:… f:t0 v:n f:y0… f:t_end i:steps i:seed`, the
//!   generator `Pcg32::new(seed, 1)`: a start time other than 0, backward integration (dt < 0, √dt NaN), the empty
//!   state, steps 0 (ZeroSteps), every diffusion, seeds 0 … u64::MAX.
//!
//! Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::grids::{Rng, pow2};
use crate::systems::{diffusion, rhs, sys_toks, traj_note, traj_toks, vec_toks};
use crate::tsv::{Class, Module, Tok, f, i};
use hisab::HisabError;
use hisab::num::{Pcg32, dopri45, euler_maruyama, gaussian_elimination, rk4_trajectory};

/// The harness's token list holds 8448 tokens (`tests/parity_core.cyr`): every row must fit.
const MAX_TOKENS: usize = 8448;

/// A u64 (`seed`) or usize (`steps`) as its i64 bits, the Cyrius word.
fn word(x: u64) -> Tok {
    i(i64::from_ne_bytes(x.to_ne_bytes()))
}

/// `E:<Variant>` for a HisabError: its `Debug` text up to the payload.
fn herr(e: &HisabError) -> Tok {
    Tok::Err(format!("{e:?}").split('(').next().unwrap_or_default().to_string())
}

fn herr_name(e: &HisabError) -> String {
    format!("Err({})", format!("{e:?}").split('(').next().unwrap_or_default())
}

fn traj_row(m: &mut Module, func: &str, class: Class, ins: &[Tok], r: &Result<Vec<(f64, Vec<f64>)>, HisabError>) {
    let (outs, summary) = match r {
        Ok(tr) => (traj_toks(tr), traj_note(tr)),
        Err(e) => (vec![herr(e)], herr_name(e)),
    };
    assert!(outs.len() <= MAX_TOKENS && ins.len() <= MAX_TOKENS, "{func}: a row of {} tokens", outs.len());
    m.row_summarised(func, class, ins, &outs, &summary);
}

#[allow(clippy::too_many_arguments)]
fn rk4(m: &mut Module, sys: &str, p: &[f64], t0: f64, y0: &[f64], t_end: f64, n: usize) {
    let fs = |t: f64, y: &[f64], dy: &mut [f64]| rhs(sys, p, t, y, dy);
    let r = rk4_trajectory(fs, bb(t0), bb(y0), bb(t_end), bb(n));
    let mut ins = sys_toks(sys, p);
    ins.push(f(t0));
    ins.extend(vec_toks(y0));
    ins.extend([f(t_end), word(n as u64)]);
    traj_row(m, "num::rk4_trajectory", Class::B, &ins, &r);
}

#[allow(clippy::too_many_arguments)]
fn dp(m: &mut Module, sys: &str, p: &[f64], t0: f64, y0: &[f64], t_end: f64, tol: f64, h: f64) {
    let fs = |t: f64, y: &[f64], dy: &mut [f64]| rhs(sys, p, t, y, dy);
    let r = dopri45(fs, bb(t0), bb(y0), bb(t_end), bb(tol), bb(h));
    let mut ins = sys_toks(sys, p);
    ins.push(f(t0));
    ins.extend(vec_toks(y0));
    ins.extend([f(t_end), f(tol), f(h)]);
    traj_row(m, "num::dopri45", Class::C, &ins, &r);
}

#[allow(clippy::too_many_arguments)]
fn em(m: &mut Module, dr: &str, p: &[f64], df: &str, q: &[f64], t0: f64, y0: &[f64], t_end: f64, n: usize, seed: u64) {
    let fd = |t: f64, y: &[f64], dy: &mut [f64]| rhs(dr, p, t, y, dy);
    let fg = |t: f64, y: &[f64], g: &mut [f64]| diffusion(df, q, t, y, g);
    let mut rng = Pcg32::new(bb(seed), 1);
    let r = euler_maruyama(fd, fg, bb(t0), bb(y0), bb(t_end), bb(n), &mut rng);
    let mut ins = sys_toks(dr, p);
    ins.extend(sys_toks(df, q));
    ins.push(f(t0));
    ins.extend(vec_toks(y0));
    ins.extend([f(t_end), word(n as u64), word(seed)]);
    traj_row(m, "num::euler_maruyama", Class::C, &ins, &r);
}

fn gauss(m: &mut Module, a: &[Vec<f64>]) {
    let mut ins = vec![Tok::V(a.len())];
    for row in a {
        ins.extend(vec_toks(row));
    }
    let mut work = a.to_vec();
    let outs = match gaussian_elimination(bb(&mut work)) {
        Ok(x) => {
            let mut v = vec_toks(&x);
            let flat: Vec<f64> = work.iter().flatten().copied().collect();
            v.extend(vec_toks(&flat));
            v
        }
        Err(e) => vec![herr(&e)],
    };
    assert!(outs.len() <= MAX_TOKENS && ins.len() <= MAX_TOKENS, "gaussian_elimination: a row too long");
    m.row("num::gaussian_elimination", Class::B, &ins, &outs);
}

fn pcg(m: &mut Module, seed: u64, count: usize) {
    let mut rng = Pcg32::new(bb(seed), 1);
    let mut outs = vec![Tok::V(count)];
    outs.extend((0..count).map(|_| i(rng.next_u32())));
    m.row("num::Pcg32::next_u32", Class::A, &[word(seed), word(count as u64)], &outs);
    let mut rng = Pcg32::new(bb(seed), 1);
    let mut outs = vec![Tok::V(count)];
    outs.extend((0..count).map(|_| f(rng.next_normal())));
    m.row("num::Pcg32::next_normal", Class::C, &[word(seed), word(count as u64)], &outs);
}

/// The n × n Hilbert matrix with right-hand side 1 (badly conditioned from n ≈ 8).
fn hilbert(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|r| {
            let mut row: Vec<f64> = (0..n).map(|c| 1.0 / f64::from(u32::try_from(r + c + 1).unwrap_or(1))).collect();
            row.push(1.0);
            row
        })
        .collect()
}

#[allow(clippy::too_many_lines)] // one routine after another
pub fn emit() -> Module {
    let mut m = Module::new("numerics");
    let mut rng = Rng::new(0x5a02);

    // ── Pcg32 (stream 1) ──
    for seed in [0, 1, 2, 7, 42, 2024, 1 << 32, 1 << 63, u64::MAX, u64::MAX - 1] {
        pcg(&mut m, seed, 64);
    }
    for _ in 0..10 {
        let seed = u64::try_from(rng.int(0, 1 << 62)).unwrap_or(0) * 3 + 1;
        pcg(&mut m, seed, 32);
    }
    pcg(&mut m, 5, 0);
    pcg(&mut m, 9, 2000);

    // ── rk4_trajectory ──
    for n in [0, 1, 2, 3, 10, 100, 1000] {
        rk4(&mut m, "decay", &[0.1], 0.0, &[1.0], 10.0, n);
    }
    for n in [1, 10, 100] {
        rk4(&mut m, "forced", &[2.0, 1.0], 1.0, &[1.0], 6.0, n); // t0 = 1
        rk4(&mut m, "decay", &[0.5], 5.0, &[2.0], -5.0, n); // backward
        rk4(&mut m, "brusselator", &[1.0, 3.0], -2.0, &[1.5, 3.0], 18.0, n);
        rk4(&mut m, "robertson", &[], 0.0, &[1.0, 0.0, 0.0], 0.001, n);
        rk4(&mut m, "consecutive", &[0.5, 0.1], 0.0, &[1.0, 0.0, 0.0], 20.0, n);
    }
    rk4(&mut m, "decay", &[0.1], 3.0, &[1.0], 3.0, 5); // t_end = t0: h = 0
    rk4(&mut m, "zero", &[], 0.0, &[], 1.0, 10); // the empty state
    rk4(&mut m, "zero", &[], 0.0, &[1.0, -0.0, 2.5], 1.0, 4); // dy never written: stays zero
    rk4(&mut m, "nan", &[], 0.0, &[1.0], 1.0, 3);
    rk4(&mut m, "decay", &[1.0], 0.0, &[f64::INFINITY], 1.0, 2);
    rk4(&mut m, "decay", &[1.0], f64::NAN, &[1.0], 1.0, 2);
    rk4(&mut m, "decay", &[1.0], 0.0, &[1.0], f64::INFINITY, 2);
    rk4(&mut m, "decay", &[1.0], 0.0, &[1.0], f64::MAX, 1);
    rk4(&mut m, "decay", &[1e3], 0.0, &[1.0], 10.0, 100); // unstable: overflow
    rk4(&mut m, "decay", &[1.0], 0.0, &[1.0], pow2(-1074), 3); // h underflows to 0
    rk4(&mut m, "decay", &[0.3], 0.0, &[1e308, 1e-308, -1.0, 0.0], 1.0, 7);
    for _ in 0..15 {
        let k = rng.log_uniform(-2.0, 1.0);
        let t0 = rng.uniform(-10.0, 10.0);
        let t1 = rng.uniform(-10.0, 10.0);
        let y0 = rng.signed_log_uniform(-3.0, 3.0);
        let n = usize::try_from(rng.int(1, 120)).unwrap_or(1);
        rk4(&mut m, "decay", &[k], t0, &[y0], t1, n);
    }

    // ── dopri45 ──
    for h in [0.0, -0.0, -1.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY, 1e-300, 1e-3, 0.1, 100.0] {
        dp(&mut m, "decay", &[0.1], 0.0, &[1.0], 10.0, 1e-8, h);
    }
    dp(&mut m, "zero", &[], 0.0, &[], 1.0, 1e-8, 0.01); // the empty state
    dp(&mut m, "decay", &[0.5], 5.0, &[2.0], -5.0, 1e-8, 0.1); // backward
    dp(&mut m, "decay", &[0.5], 5.0, &[2.0], 5.0, 1e-8, 0.1); // t_end = t0
    dp(&mut m, "forced", &[2.0, 1.0], 1.0, &[1.0], 6.0, 1e-8, 0.05);
    dp(&mut m, "nan", &[], 0.0, &[1.0], 1.0, 1e-6, 0.01); // the 100 000-step limit
    dp(&mut m, "zero", &[], 0.0, &[1.0], 1e6, 1e-8, 1e-6); // ×5 steps to the end
    for tol in [1e-3, 1e-6, 1e-9, 1e-12] {
        dp(&mut m, "brusselator", &[1.0, 3.0], 0.0, &[1.5, 3.0], 5.0, tol, 0.05);
        dp(&mut m, "decay2", &[3.0], -1.0, &[2.0], 1.0, tol, 0.02);
    }
    dp(&mut m, "decay", &[1.0], 0.0, &[1.0], 1.0, f64::MIN_POSITIVE, 0.01); // tol tiny: tiny steps
    dp(&mut m, "decay", &[1.0], 0.0, &[1.0], 1.0, f64::NAN, 0.01); // NaN tol: never accepted
    dp(&mut m, "decay", &[1.0], 0.0, &[1.0], 1.0, -1e-6, 0.01);
    dp(&mut m, "decay", &[1.0], 0.0, &[1e-300, 1e300], 1.0, 1e-8, 0.01);
    for _ in 0..12 {
        let k = rng.log_uniform(-2.0, 1.0);
        let t0 = rng.uniform(-5.0, 5.0);
        let t1 = rng.uniform(-5.0, 5.0);
        let y0 = rng.signed_log_uniform(-3.0, 3.0);
        let tol = rng.log_uniform(-10.0, -4.0);
        let h = rng.log_uniform(-4.0, 0.0);
        dp(&mut m, "decay", &[k], t0, &[y0], t1, tol, h);
    }

    // ── gaussian_elimination ──
    let eps = 1e-12_f64; // hisab 1.4.0's EPSILON_F64
    gauss(&mut m, &[]);
    gauss(&mut m, &[vec![2.0, 4.0]]);
    gauss(&mut m, &[vec![0.0, 1.0]]);
    for piv in [eps, eps.next_up(), eps.next_down(), -eps, -eps.next_down(), 2e-12, 0.0, -0.0, f64::NAN, f64::INFINITY] {
        gauss(&mut m, &[vec![piv, 1.0]]);
        gauss(&mut m, &[vec![piv, 1.0, 2.0], vec![0.0, 3.0, 4.0]]);
    }
    gauss(&mut m, &[vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 2.0]]); // a swap
    gauss(&mut m, &[vec![1.0, 2.0, 3.0], vec![-1.0, 5.0, 6.0]]); // |−1| = |1|: no swap (strictly larger only)
    gauss(&mut m, &[vec![1.0, 2.0, 3.0], vec![2.0, 4.0, 6.0]]); // singular
    gauss(&mut m, &[vec![f64::NAN, 2.0, 3.0], vec![1.0, 4.0, 6.0]]); // a NaN pivot is never replaced
    gauss(&mut m, &[vec![1.0, f64::NAN, 3.0], vec![2.0, 4.0, 6.0]]);
    gauss(&mut m, &[vec![1.0, 2.0, f64::INFINITY], vec![3.0, 4.0, 6.0]]);
    gauss(&mut m, &[vec![1e300, 1e300, 1.0], vec![1e-300, 1.0, 2.0]]);
    gauss(&mut m, &[vec![0.0, 0.0, 1.0, 1.0], vec![0.0, 1.0, 0.0, 2.0], vec![1.0, 0.0, 0.0, 3.0]]);
    gauss(&mut m, &[vec![4.0, -2.0, 1.0, 11.0], vec![-2.0, 4.0, -2.0, -16.0], vec![1.0, -2.0, 4.0, 17.0]]);
    for n in [2, 4, 6, 8, 10, 12] {
        gauss(&mut m, &hilbert(n));
    }
    for n in 1..=12 {
        for _ in 0..4 {
            let a: Vec<Vec<f64>> = (0..n).map(|_| (0..=n).map(|_| rng.uniform(-10.0, 10.0)).collect()).collect();
            gauss(&mut m, &a);
        }
    }
    for _ in 0..10 {
        let n = usize::try_from(rng.int(2, 6)).unwrap_or(2);
        let a: Vec<Vec<f64>> = (0..n).map(|_| (0..=n).map(|_| rng.signed_log_uniform(-200.0, 200.0)).collect()).collect();
        gauss(&mut m, &a);
    }

    // ── euler_maruyama ──
    for n in [0, 1, 10, 100, 1000] {
        em(&mut m, "decay", &[1.0], "const", &[0.1], 0.0, &[1.0], 1.0, n, 42);
    }
    em(&mut m, "decay", &[1.0], "const", &[0.1], 2.0, &[1.0], 3.0, 20, 7); // t0 = 2
    em(&mut m, "decay", &[1.0], "const", &[0.1], 1.0, &[1.0], 0.0, 5, 7); // backward: √dt is NaN
    em(&mut m, "zero", &[], "none", &[], 0.0, &[], 1.0, 10, 1); // the empty state
    em(&mut m, "zero", &[], "none", &[], 0.0, &[1.0, 2.0], 1.0, 10, 1);
    em(&mut m, "reversible", &[1.0, 0.5], "mult", &[0.3], 0.0, &[1.0, 0.0], 5.0, 50, 3);
    em(&mut m, "consecutive", &[0.5, 0.1], "zero", &[], 0.0, &[1.0, 0.0, 0.0], 5.0, 50, 3);
    for seed in [0, 1, 1 << 63, u64::MAX] {
        em(&mut m, "decay", &[0.5], "const", &[1.0], 0.0, &[0.0], 1.0, 64, seed);
    }
    m
}
