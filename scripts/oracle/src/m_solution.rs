//! `solution.tsv` — concentration units, pH, weak acids and bases, ionic strength, activity coefficients and
//! conductivity (port plan §4, M3a).
//!
//! - Class B (+ − × ÷ and √ only, expected 0 ulp): `molarity`, `molality`, `dilution`, `ph_from_poh`,
//!   `ionic_strength`, `kohlrausch_molar_conductivity`, `limiting_molar_conductivity`, `solution_conductivity`.
//!   Each argument over the validation edges with the others typical, the Rust tests' values, overflow and a
//!   seeded sweep. `ionic_strength`'s input is `v:n` then n pairs `f:c i:z` (the charge an i32): the Rust tests'
//!   lists, the empty list and zero / −0 mixes (Rust's `.sum()` folds from −0.0), i32::MIN / MAX charges,
//!   NaN / ±∞ / negative concentrations, and seeded lists of 1–8 ions.
//! - Class C, one `log10`: `ph_from_h_concentration`, `poh_from_oh_concentration` (the edges, every decade from
//!   the smallest subnormal to DBL_MAX with four mantissas, the Rust tests' values, a seeded sweep) and
//!   `henderson_hasselbalch` (the edges per argument, the Rust tests, a seeded sweep, and **near-cancellation
//!   rows on purpose**: pKa = −log₁₀(base/acid) to within a few ulp, so pH ≈ 0 of both signs, over ratios from
//!   10⁻³⁰⁰ to 10³⁰⁰, and the exact powers of ten, where the cancellation is exact).
//! - Class C, one `10^x`: `h_concentration_from_ph` (the edges, pH −310 … 330 across 10^x's overflow and
//!   underflow, a seeded sweep) and `debye_huckel_limiting`, `debye_huckel_extended`, `davies_activity` (input
//!   `i:z` first: every charge from −4 to 4, ±50 000 and i32::MIN / MAX; I, A, B and å over their edges; the Rust
//!   tests' values incl. the singular ion size; a seeded sweep).
//! - Class C, the charge-balance Newton then one `log10`: `weak_acid_ph`, `weak_base_ph` over the **Ka (or Kb) ×
//!   C grid of every decade 10⁻⁴⁰ … 10⁴⁰ × 10⁻⁴⁰ … 10²⁰** (port plan §4), seeded log-uniform pairs over the same
//!   box and over the chemically realistic one (Ka 10⁻¹⁴ … 10, C 10⁻¹² … 10), the edges, the overflow cases
//!   (Ka·C or x³ past DBL_MAX) and the Rust tests' cases.
//!
//! **Newton paths.** `root_replica` transcribes Rust's private `charge_balance_root` (and `checked_root`) and
//! must give the oracle's pH bits, or its error, on every weak acid / base row, or generation fails. Every row
//! that reaches the iteration carries its exit in the note — `[root: f(x) <= 0, k updates]`, `[root: relative
//! step, k updates]`, `[root: non-finite step, k updates]`, `[root: 50 iterations]`, with `, checked_root
//! rejected` when the accepted x was not positive and finite — and the generator prints the tally and the
//! largest update count.
//!
//! - ADR 0002 #29 (the overflow audit's B2): `dilution` and `solution_conductivity` rows are noted where the port's
//!   correctly rounded M₁·V₁/V₂ or Λ·c/1000 differs from Rust's value (`audit.rs`), and `audit::emit_solution` appends
//!   their edge rows (partial products past DBL_MAX or below 2⁻¹⁰²², midpoints and subnormal ties) with their own seed.
//! - ADR 0002 #33 (the overflow audit's B4): `ionic_strength`, `kohlrausch_molar_conductivity` and
//!   `limiting_molar_conductivity` rows are noted where Rust's sequence with an unbounded exponent (the final product
//!   or subtraction rounded once) differs from Rust's value (`audit.rs`), and `audit::emit_solution_b4` appends their
//!   edge rows (terms, sums or products past DBL_MAX or below 2⁻¹⁰²² that the rest brings back) with their own seed.
//! - ADR 0002 #35 (the overflow audit's B5): `debye_huckel_limiting`, `debye_huckel_extended` and `davies_activity` rows
//!   are noted where the port's γ — log₁₀ γ with an unbounded exponent (its final product or quotient rounded once),
//!   then the correctly rounded 10^x — differs from Rust's value (`audit.rs`), and `audit::emit_solution_b5` appends
//!   their edge rows (A·z², B·å, (B·å)·√I or the numerator past DBL_MAX, at I = 0 or brought back, Davies' bracket 0,
//!   subnormal products) with their own seed.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations
//! (`grids.rs`). The near-cancellation pKa values come from `log10_arith`, + − × ÷ only.

use std::hint::black_box as bb;

use crate::audit;
use crate::grids::{Rng, decades};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, f, i};
use kimiya::error::KimiyaError;
use kimiya::solution::{
    DH_A_WATER_25C, DH_B_WATER_25C, KW_WATER_25C, LAMBDA_CL_MINUS, LAMBDA_H_PLUS, LAMBDA_K_PLUS,
    LAMBDA_NA_PLUS, LAMBDA_NO3_MINUS, LAMBDA_OH_MINUS, LAMBDA_SO4_2MINUS, h_concentration_from_ph,
    henderson_hasselbalch,
    molality, molarity, ph_from_h_concentration, ph_from_poh, poh_from_oh_concentration,
    weak_acid_ph, weak_base_ph,
};

// ── the charge-balance replica ──────────────────────────────────────────────────────

/// How Rust's `charge_balance_root` ended.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Exit {
    /// `fx <= 0.0`: x accepted where the iterates met the root.
    FxNonPositive,
    /// `step <= CHARGE_BALANCE_REL_TOL * x`, after the update.
    RelStep,
    /// `!step.is_finite()`: the loop breaks (K·C or x³ overflowed).
    NonFinite,
    /// `CHARGE_BALANCE_MAX_ITER` iterations without acceptance.
    Exhausted,
}

/// The exit, the number of `x -= step` updates, and whether `checked_root` refused the accepted x.
#[derive(Clone, Copy)]
pub struct RootPath {
    pub exit: Exit,
    pub updates: usize,
    pub rejected: bool,
}

/// Rust's private `charge_balance_root(k, c)` and `checked_root`, transcribed operation for operation from
/// `rust-old/src/solution.rs`, with the path it took. `Err(())` is its `ComputationError`.
pub fn root_replica(k: f64, c: f64, tol: f64, max_iter: usize) -> (Result<f64, ()>, RootPath) {
    let checked = |x: f64| if x > 0.0 && x.is_finite() { Ok(x) } else { Err(()) };
    let kc_kw = k * c + KW_WATER_25C;
    let weak_bound = kc_kw.sqrt();
    let strong_bound = 0.5 * (c + (c * c + 4.0 * KW_WATER_25C).sqrt());
    let mut x = weak_bound.min(strong_bound);
    for it in 0..max_iter {
        let fx = x * x * (x + k) - k * (c * x + KW_WATER_25C) - KW_WATER_25C * x;
        if fx <= 0.0 {
            let r = checked(x);
            return (r, RootPath { exit: Exit::FxNonPositive, updates: it, rejected: r.is_err() });
        }
        let dfx = x * (3.0 * x + 2.0 * k) - kc_kw;
        let step = fx / dfx;
        if !step.is_finite() {
            return (Err(()), RootPath { exit: Exit::NonFinite, updates: it, rejected: false });
        }
        x -= step;
        if step <= tol * x {
            let r = checked(x);
            return (r, RootPath { exit: Exit::RelStep, updates: it + 1, rejected: r.is_err() });
        }
    }
    (Err(()), RootPath { exit: Exit::Exhausted, updates: max_iter, rejected: false })
}

/// Rust's `CHARGE_BALANCE_REL_TOL` and `CHARGE_BALANCE_MAX_ITER` (private; transcribed).
const REL_TOL: f64 = 1e-12;
const MAX_ITER: usize = 50;

/// The note suffix for a path.
fn path_note(p: RootPath) -> String {
    let exit = match p.exit {
        Exit::FxNonPositive => format!("f(x) <= 0, {} updates", p.updates),
        Exit::RelStep => format!("relative step, {} updates", p.updates),
        Exit::NonFinite => format!("non-finite step, {} updates", p.updates),
        Exit::Exhausted => format!("{} iterations", p.updates),
    };
    let rej = if p.rejected { ", checked_root rejected" } else { "" };
    format!("root: {exit}{rej}")
}

/// Path tallies for the module comment and port-parity.md.
#[derive(Default)]
struct Paths {
    by_exit: [usize; 4],
    rejected: usize,
    max_updates: usize,
    hist: [usize; 12],
    guard_err: usize,
}

fn variant(e: &KimiyaError) -> String {
    format!("{e:?}").split('(').next().unwrap_or_default().to_string()
}

fn pos(x: f64) -> bool {
    x > 0.0 && x.is_finite()
}

/// One weak acid or base row: the oracle's value, checked against the replica, with the path in the note.
fn weak(m: &mut Module, paths: &mut Paths, base: bool, k: f64, c: f64) {
    let got = if base { weak_base_ph(bb(k), bb(c)) } else { weak_acid_ph(bb(k), bb(c)) };
    let name = if base { "solution::weak_base_ph" } else { "solution::weak_acid_ph" };
    let (rep, path): (Result<f64, &str>, Option<RootPath>) = if !pos(k) {
        (Err("InvalidInput"), None)
    } else if !pos(c) {
        (Err("InvalidConcentration"), None)
    } else {
        let (r, p) = root_replica(k, c, REL_TOL, MAX_ITER);
        let v = r.map(|x| {
            let l = -x.log10();
            if base { 14.0 - l } else { l }
        });
        (v.map_err(|()| "ComputationError"), Some(p))
    };
    let same = match (&got, &rep) {
        (Ok(a), Ok(b)) => a.to_bits() == b.to_bits(),
        (Err(a), Err(b)) => variant(a) == *b,
        _ => false,
    };
    assert!(same, "root replica disagrees with the oracle: {name}({k:?}, {c:?}): {got:?} vs {rep:?}");
    let ins = [f(k), f(c)];
    let out = res(got);
    match path {
        Some(p) => {
            let e = match p.exit {
                Exit::FxNonPositive => 0,
                Exit::RelStep => 1,
                Exit::NonFinite => 2,
                Exit::Exhausted => 3,
            };
            paths.by_exit[e] += 1;
            paths.rejected += usize::from(p.rejected);
            paths.max_updates = paths.max_updates.max(p.updates);
            paths.hist[p.updates.min(11)] += 1;
            m.row_noted(name, Class::C, &ins, &[out], &path_note(p));
        }
        None => {
            paths.guard_err += 1;
            m.row(name, Class::C, &ins, &[out]);
        }
    }
}

/// log₁₀(x) for a normal x > 0 from + − × ÷ and bit manipulation only (grids may not call a transcendental
/// std method, `grids.rs`): x = 2ᵉ·m, m in [√½, √2), ln m = 2·atanh((m − 1)/(m + 1)) by its series. Within a
/// few ulp: it only places the near-cancellation pKa values, whose rows are near cancellation either way.
fn log10_arith(x: f64) -> f64 {
    let bits = x.to_bits();
    let biased = i32::try_from((bits >> 52) & 0x7ff).unwrap_or(1023);
    let mut e = f64::from(biased - 1023);
    let mut m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    if m > std::f64::consts::SQRT_2 {
        m *= 0.5;
        e += 1.0;
    }
    let s = (m - 1.0) / (m + 1.0);
    let s2 = s * s;
    let (mut term, mut sum, mut k) = (s, s, 1.0);
    for _ in 0..40 {
        term *= s2;
        k += 2.0;
        sum += term / k;
    }
    (2.0 * sum + e * std::f64::consts::LN_2) / std::f64::consts::LN_10
}

/// `x` moved `d` ulps (up for d > 0).
fn ulps(x: f64, d: i32) -> f64 {
    let mut y = x;
    for _ in 0..d.unsigned_abs() {
        y = if d > 0 { y.next_up() } else { y.next_down() };
    }
    y
}

fn ions_toks(xs: &[(f64, i32)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(xs.len())];
    for &(c, z) in xs {
        v.push(f(c));
        v.push(i(z));
    }
    v
}

/// The charges every activity model is probed with.
const CHARGES: [i32; 13] = [i32::MIN, -50_000, -4, -3, -2, -1, 0, 1, 2, 3, 4, 50_000, i32::MAX];

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("solution");
    let mut rng = Rng::new(0x3a01);

    // ── molarity, molality, dilution ──
    let two = |m: &mut Module, name: &str, x: f64, y: f64| {
        let out = match name {
            "solution::molarity" => res(molarity(bb(x), bb(y))),
            _ => res(molality(bb(x), bb(y))),
        };
        m.row(name, Class::B, &[f(x), f(y)], &[out]);
    };
    for name in ["solution::molarity", "solution::molality"] {
        for x in edges() {
            two(&mut m, name, x, 0.25);
            two(&mut m, name, 0.5, x);
        }
        for (x, y) in [(0.5, 0.25), (1.0, 0.0), (1e308, 1e-308), (-2.0, 3.0), (0.0, 1e-300)] {
            two(&mut m, name, x, y);
        }
        for _ in 0..150 {
            let x = rng.signed_log_uniform(-6.0, 4.0);
            let y = rng.log_uniform(-4.0, 3.0);
            two(&mut m, name, x, y);
        }
    }
    let dil = audit::dilution_row;
    for x in edges() {
        dil(&mut m, x, 0.1, 0.5);
        dil(&mut m, 1.0, x, 0.5);
        dil(&mut m, 1.0, 0.1, x);
    }
    for (a, b, c) in [(1.0, 0.1, 0.5), (1.0, 0.1, 0.0), (1e300, 1e300, 1.0), (1e-300, 1e-300, 1e300)] {
        dil(&mut m, a, b, c);
    }
    for _ in 0..150 {
        let a = rng.log_uniform(-6.0, 2.0);
        let b = rng.log_uniform(-3.0, 3.0);
        let c = rng.log_uniform(-3.0, 3.0);
        dil(&mut m, a, b, c);
    }

    // ── pH, pOH ──
    let mut conc = edges();
    conc.extend(decades(-323, 308, &["1", "2", "3.16", "5"]));
    conc.extend([1e-7, 0.1, 1e-4, h_concentration_from_ph(4.5), f64::MIN_POSITIVE.next_down()]);
    for _ in 0..300 {
        conc.push(rng.log_uniform(-320.0, 308.0));
    }
    for &x in &conc {
        let ph = res(ph_from_h_concentration(bb(x)));
        m.row("solution::ph_from_h_concentration", Class::C, &[f(x)], &[ph]);
        let poh = res(poh_from_oh_concentration(bb(x)));
        m.row("solution::poh_from_oh_concentration", Class::C, &[f(x)], &[poh]);
    }
    let mut pohs = edges();
    pohs.extend([11.0, 7.0, 14.0, -1.0, 15.0, 14.0_f64.next_up(), 1e16, -1e16]);
    for _ in 0..150 {
        pohs.push(rng.uniform(-5.0, 20.0));
    }
    for x in pohs {
        let out = ph_from_poh(bb(x));
        m.row("solution::ph_from_poh", Class::B, &[f(x)], &[f(out)]);
    }

    // ── Henderson–Hasselbalch ──
    let hh = |m: &mut Module, pka: f64, b: f64, a: f64| {
        let out = res(henderson_hasselbalch(bb(pka), bb(b), bb(a)));
        m.row("solution::henderson_hasselbalch", Class::C, &[f(pka), f(b), f(a)], &[out]);
    };
    for x in edges() {
        hh(&mut m, x, 0.1, 0.1);
        hh(&mut m, 4.75, x, 0.1);
        hh(&mut m, 4.75, 0.1, x);
    }
    for (pka, b, a) in [
        (4.75, 0.1, 0.1),
        (4.75, 1.0, 0.1),
        (4.75, 0.1, 0.0),
        (4.75, 0.0, 0.1),
        (4.75, f64::NAN, 0.1),
        (4.75, 0.1, f64::INFINITY),
        (f64::NAN, 0.1, 0.1),
        (4.75, f64::MAX, f64::MIN_POSITIVE), // the ratio overflows to +∞: log10 = +∞
        (4.75, f64::MIN_POSITIVE, f64::MAX), // the ratio underflows to 0: log10 = −∞
    ] {
        hh(&mut m, pka, b, a);
    }
    for _ in 0..400 {
        let pka = rng.uniform(-5.0, 20.0);
        let b = rng.log_uniform(-10.0, 2.0);
        let a = rng.log_uniform(-10.0, 2.0);
        hh(&mut m, pka, b, a);
    }
    // Near cancellation, on purpose: pKa ≈ −log10(base/acid), so the pH is a few ulp of the log either side
    // of 0 — one ulp of the log is many of the pH (the reason for _km_log10_cr, ADR 0005 amendment M3a).
    for _ in 0..60 {
        let a = rng.log_uniform(-12.0, 2.0);
        let b = a * rng.log_uniform(-15.0, 15.0);
        if !pos(b) || !pos(b / a) {
            continue;
        }
        let pka0 = -log10_arith(b / a);
        for d in -2..=2 {
            hh(&mut m, ulps(pka0, d), b, a);
        }
    }
    for _ in 0..20 {
        let a = rng.log_uniform(-150.0, 150.0);
        let b = rng.log_uniform(-150.0, 150.0);
        if !pos(b / a) || b / a < f64::MIN_POSITIVE {
            continue;
        }
        let pka0 = -log10_arith(b / a);
        for d in [-1, 0, 1] {
            hh(&mut m, ulps(pka0, d), b, a);
        }
    }
    // Exact cancellation: base/acid an exact power of ten, log10 an exact integer.
    for k in [1_i32, 2, 5, 10, 15, 22] {
        let b: f64 = format!("1e{k}").parse().unwrap_or(f64::NAN);
        for d in [-1, 0, 1] {
            hh(&mut m, ulps(-f64::from(k), d), b, 1.0);
            hh(&mut m, ulps(f64::from(k), d), 1.0, b);
        }
    }

    // ── 10^(−pH) ──
    let mut phs = edges();
    phs.extend([4.5, 0.0, 7.0, 14.0, -1.0, 2.875, 307.0, 308.0, 308.3, 308.26, 309.0, 323.3, 323.31, 324.0]);
    phs.extend([-307.0, -308.0, -308.25, -308.26, -309.0, -310.0, 330.0]);
    for _ in 0..300 {
        phs.push(rng.uniform(-310.0, 330.0));
    }
    for _ in 0..100 {
        phs.push(rng.uniform(-1.0, 15.0));
    }
    for x in phs {
        let out = h_concentration_from_ph(bb(x));
        m.row("solution::h_concentration_from_ph", Class::C, &[f(x)], &[f(out)]);
    }

    // ── weak acids and bases ──
    let mut paths = Paths::default();
    let ks = decades(-40, 40, &["1"]);
    let cs = decades(-40, 20, &["1"]);
    for base in [false, true] {
        for &k in &ks {
            for &c in &cs {
                weak(&mut m, &mut paths, base, k, c);
            }
        }
        for _ in 0..2000 {
            let k = rng.log_uniform(-40.0, 40.0);
            let c = rng.log_uniform(-40.0, 20.0);
            weak(&mut m, &mut paths, base, k, c);
        }
        // the chemically realistic box (Ka 1e-14 … 10, C 1e-12 … 10), where water and the acid compete
        for _ in 0..2000 {
            let k = rng.log_uniform(-14.0, 1.0);
            let c = rng.log_uniform(-12.0, 1.0);
            weak(&mut m, &mut paths, base, k, c);
        }
        for x in edges() {
            weak(&mut m, &mut paths, base, x, 0.1);
            weak(&mut m, &mut paths, base, 1.8e-5, x);
        }
        // the Rust tests' cases
        for (k, c) in [
            (1.8e-5, 0.1),
            (6.8e-4, 0.1),
            (1.8e-5, 1e-5),
            (1.8e-5, 1e-8),
            (5e-7, 2e-6),
            (1e-12, 0.01),
            (1e-20, 1.0),
            (1e-30, 1e-30),
            (1.0, 0.1),
            (1e3, 0.1),
            (1e10, 1e-3),
            (0.0, 0.1),
            (1.8e-5, 0.0),
            (f64::NAN, 0.1),
            (1.8e-5, f64::NAN),
            (f64::INFINITY, 0.1),
            (1.8e-5, f64::INFINITY),
            (1e300, 1e300),
        ] {
            weak(&mut m, &mut paths, base, k, c);
        }
        // overflow: K·C past DBL_MAX, and x³ (x ≈ C for K ≫ C) past it, across both ends
        for &(k, c) in &[(1e200, 1e200), (1e154, 1e154), (1e155, 1e154), (f64::MAX, 2.0), (f64::MAX, 1.0), (2.0, f64::MAX)] {
            weak(&mut m, &mut paths, base, k, c);
        }
        for e in [90, 100, 102, 103, 104, 105, 110, 150, 154, 155, 160, 200, 300] {
            let c: f64 = format!("1e{e}").parse().unwrap_or(f64::NAN);
            for k in [1e-10, 1.0, 1e10, 1e100, 1e300] {
                weak(&mut m, &mut paths, base, k, c);
            }
        }
        for e in [100, 150, 200, 250, 290, 300, 305, 308] {
            let k: f64 = format!("1e{e}").parse().unwrap_or(f64::NAN);
            for c in [1e-300, 1e-100, 1e-10, 1.0, 1e10] {
                weak(&mut m, &mut paths, base, k, c);
            }
        }
        // the subnormal and tiny corner: K·C and C² underflow, the root is water's 1e-7
        for (k, c) in [(5e-324, 5e-324), (1e-300, 1e-300), (f64::MIN_POSITIVE, 1.0), (1.0, 5e-324), (1e-200, 1e200)] {
            weak(&mut m, &mut paths, base, k, c);
        }
    }
    eprintln!(
        "solution: weak acid/base rows reaching the Newton: {} f(x) <= 0, {} relative step, {} non-finite step, \
         {} exhausted, {} checked_root rejections; max {} updates; updates histogram {:?}; {} guard Err rows",
        paths.by_exit[0],
        paths.by_exit[1],
        paths.by_exit[2],
        paths.by_exit[3],
        paths.rejected,
        paths.max_updates,
        paths.hist,
        paths.guard_err
    );

    // ── ionic strength ──
    let ion = |m: &mut Module, xs: &[(f64, i32)]| audit::ionic(m, &ions_toks(xs), xs);
    let lists: Vec<Vec<(f64, i32)>> = vec![
        vec![],
        vec![(0.1, 1), (0.1, -1)],
        vec![(0.1, 2), (0.2, -1)],
        vec![(1.0, 50_000), (1.0, i32::MIN)],
        vec![(1.0, i32::MAX)],
        vec![(0.0, 1)],
        vec![(-0.0, 1)],
        vec![(-0.0, 1), (-0.0, -1)],
        vec![(0.5, 0)],
        vec![(-0.5, 0)],
        vec![(-0.1, 2)],
        vec![(f64::NAN, 1)],
        vec![(f64::INFINITY, 1), (1.0, 1)],
        vec![(f64::INFINITY, 1), (f64::NEG_INFINITY, 1)],
        vec![(f64::INFINITY, 0)],
        vec![(1e308, 4), (1e308, 4)],
        vec![(f64::MIN_POSITIVE, 1), (5e-324, 3)],
        vec![(1e16, 1), (1.0, 1), (-1e16, 1)],
        vec![(0.01, 3), (0.03, -1)],
        vec![(0.05, 2), (0.05, -2)],
    ];
    for xs in &lists {
        ion(&mut m, xs);
    }
    for _ in 0..300 {
        let n = rng.int(1, 8);
        let xs: Vec<(f64, i32)> = (0..n)
            .map(|_| (rng.log_uniform(-6.0, 1.0), i32::try_from(rng.int(-4, 4)).unwrap_or(1)))
            .collect();
        ion(&mut m, &xs);
    }

    // ── activity coefficients ──
    let lim = audit::dh_lim;
    let dav = audit::davies;
    let ext = audit::dh_ext;
    let ionic = [0.0, -0.0, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 10.0, 1e300];
    for z in CHARGES {
        for ii in ionic {
            lim(&mut m, z, ii, DH_A_WATER_25C);
            dav(&mut m, z, ii, DH_A_WATER_25C);
            ext(&mut m, z, ii, DH_A_WATER_25C, DH_B_WATER_25C, 3.0);
        }
    }
    for x in edges() {
        lim(&mut m, 1, x, DH_A_WATER_25C);
        lim(&mut m, 2, 0.01, x);
        dav(&mut m, 1, x, DH_A_WATER_25C);
        dav(&mut m, 2, 0.01, x);
        ext(&mut m, 1, x, DH_A_WATER_25C, DH_B_WATER_25C, 3.0);
        ext(&mut m, 2, 0.01, x, DH_B_WATER_25C, 3.0);
        ext(&mut m, 1, 0.01, DH_A_WATER_25C, x, 3.0);
        ext(&mut m, 1, 0.01, DH_A_WATER_25C, DH_B_WATER_25C, x);
    }
    let singular = -1.0 / (DH_B_WATER_25C * 0.1);
    for (ii, a, b, s) in [
        (0.01, DH_A_WATER_25C, DH_B_WATER_25C, singular),
        (0.01, DH_A_WATER_25C, DH_B_WATER_25C, -3.0),
        (0.01, DH_A_WATER_25C, -DH_B_WATER_25C, 3.0),
        (0.01, DH_A_WATER_25C, DH_B_WATER_25C, 0.0),
        (0.05, DH_A_WATER_25C, DH_B_WATER_25C, 3.0),
        (0.01, DH_A_WATER_25C, 1e300, 1e300), // 1 + B·å·√I overflows: γ = 10^-0 = 1
    ] {
        ext(&mut m, 1, ii, a, b, s);
    }
    for (z, ii) in [(1, 0.01), (2, 0.01), (1, 0.0), (1, 0.05), (1, 0.1), (i32::MAX, 0.01), (50_000, 0.01)] {
        lim(&mut m, z, ii, DH_A_WATER_25C);
        dav(&mut m, z, ii, DH_A_WATER_25C);
    }
    ext(&mut m, i32::MIN, 0.01, 0.509, 0.328, 3.0);
    for (z, ii) in [(1, -0.1), (1, f64::NAN), (1, f64::INFINITY)] {
        lim(&mut m, z, ii, 0.509);
        dav(&mut m, z, ii, 0.509);
        ext(&mut m, z, ii, 0.509, 0.328, 3.0);
    }
    lim(&mut m, 1, 0.01, f64::NAN);
    dav(&mut m, 1, 0.01, f64::INFINITY);
    for _ in 0..300 {
        let z = i32::try_from(rng.int(-4, 4)).unwrap_or(1);
        let ii = rng.log_uniform(-6.0, 0.0);
        let a = rng.uniform(0.4, 0.6);
        let b = rng.uniform(0.3, 0.35);
        let s = rng.uniform(2.0, 10.0);
        lim(&mut m, z, ii, a);
        dav(&mut m, z, ii, a);
        ext(&mut m, z, ii, a, b, s);
    }

    // ── conductivity ──
    let kohl = audit::kohl;
    for x in edges() {
        kohl(&mut m, x, 10.0, 0.01);
        kohl(&mut m, 126.4, x, 0.01);
        kohl(&mut m, 126.4, 10.0, x);
    }
    for (l0, k, c) in [(126.4, 10.0, 0.01), (126.4, 10.0, 0.1), (126.4, 10.0, f64::NAN), (f64::NAN, 10.0, 0.1)] {
        kohl(&mut m, l0, k, c);
    }
    for _ in 0..200 {
        let l0 = rng.uniform(50.0, 450.0);
        let k = rng.uniform(0.0, 200.0);
        let c = rng.log_uniform(-6.0, 0.5);
        kohl(&mut m, l0, k, c);
    }
    let lambdas = [
        LAMBDA_H_PLUS,
        LAMBDA_OH_MINUS,
        LAMBDA_NA_PLUS,
        LAMBDA_K_PLUS,
        LAMBDA_CL_MINUS,
        LAMBDA_NO3_MINUS,
        LAMBDA_SO4_2MINUS,
    ];
    let lmc = audit::limiting;
    for &cat in &lambdas[..4] {
        for &an in &[lambdas[1], lambdas[4], lambdas[5], lambdas[6]] {
            for (nc, na) in [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0)] {
                lmc(&mut m, cat, nc, an, na);
            }
        }
    }
    for x in edges() {
        lmc(&mut m, x, 1.0, 76.3, 1.0);
        lmc(&mut m, 50.1, x, 76.3, 1.0);
        lmc(&mut m, 50.1, 1.0, x, 1.0);
        lmc(&mut m, 50.1, 1.0, 76.3, x);
    }
    for _ in 0..150 {
        let a = rng.uniform(30.0, 360.0);
        let b = f64::from(i32::try_from(rng.int(1, 3)).unwrap_or(1));
        let c = rng.uniform(30.0, 200.0);
        let d = f64::from(i32::try_from(rng.int(1, 3)).unwrap_or(1));
        lmc(&mut m, a, b, c, d);
    }
    let cond = audit::conductivity;
    for x in edges() {
        cond(&mut m, x, 0.1);
        cond(&mut m, 126.4, x);
    }
    for _ in 0..150 {
        let a = rng.uniform(50.0, 450.0);
        let c = rng.log_uniform(-6.0, 1.0);
        cond(&mut m, a, c);
    }
    audit::emit_solution(&mut m);
    audit::emit_solution_b4(&mut m);
    audit::emit_solution_b5(&mut m);
    m
}
