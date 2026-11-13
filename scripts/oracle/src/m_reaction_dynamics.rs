//! `reaction_dynamics.tsv` — concentration trajectories by ODE integration (port plan §4, M5a).
//!
//! - Class B (+ − × ÷ and √ only, RK4 and backward Euler, expected 0 ulp): the five built-in rate laws,
//!   `simulate_first_order`, `simulate_second_order`, `simulate_consecutive`, `simulate_reversible` and
//!   `simulate_michaelis_menten` (inputs their floats, then `i:steps`), on mechanism presets at 10, 100 and 1000 steps
//!   (the Rust tests' and benches' values among them), each argument over the validation edges, the steps edges (0, 1,
//!   10⁷ + 1, usize::MAX), RK4's stability boundary and overflow (the unstable-RK4 Err), extreme magnitudes and seeded
//!   sweeps; `simulate_stiff` on the systems of `systems.rs` (input `s:system v:n f:params… v:n f:y0… f:t_end
//!   i:steps`): the Rust tests' (linear and nonlinear recurrences, Robertson, the singular and the wrong-sign Newton
//!   systems, a NaN derivative), the mechanisms at 10, 100 and 1000 steps, Robertson at 10 … 4000 steps, a Newton that
//!   converges linearly (`slowjac`) swept across the 20-iteration limit and both tolerances, the empty system at 10⁷
//!   steps, the Err edges.
//! - Class C (`powf` in dopri45's step factor; `ln` and `cos` in the Box–Muller draws, glibc's in the oracle):
//!   `simulate_adaptive` (`s:system v:n f:params… v:n f:y0… f:t_end f:tol`) on the mechanisms and the edges (tolerance
//!   and state Err, the empty state, a vanishing first step, a derivative that is NaN or never written, a blow-up and
//!   the 100 000-step limit, tolerances 10⁻³ … 10⁻¹²); `simulate_stochastic` (`s:drift v:n f:… s:diffusion v:n f:… v:n
//!   f:y0… f:t_end i:steps i:seed`) with additive, multiplicative, zero and unwritten noise, at 10, 100 and 1000 steps,
//!   on seeds 0, 1, 7, 42, 2⁶³ and u64::MAX (a u64 is its i64 bits), the Rust tests' runs and the Err edges.
//!
//! Every output is the trajectory, `v:points` then each point's `f:t v:n f:y…`, or `E:Variant`; a row's note gives its
//! point count and last point. Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::grids::Rng;
use crate::m_reaction::edges;
use crate::systems::{diffusion, jac, rhs, sys_toks, traj_note, traj_toks, vec_toks};
use crate::tsv::{Class, Module, Tok, err, f, i};
use kimiya::error::KimiyaError;
use kimiya::reaction_dynamics::{
    TrajectoryPoint, simulate_adaptive, simulate_consecutive, simulate_first_order, simulate_michaelis_menten,
    simulate_reversible, simulate_second_order, simulate_stiff, simulate_stochastic,
};

/// The harness's token list holds 8448 tokens (`tests/parity_core.cyr`): every row must fit.
const MAX_TOKENS: usize = 8448;

/// A usize (`steps`) or u64 (`seed`) as its i64 bits, the Cyrius word.
fn word(x: u64) -> Tok {
    i(i64::from_ne_bytes(x.to_ne_bytes()))
}

/// One trajectory row: the outputs, their summary in the note.
fn traj_row(m: &mut Module, func: &str, class: Class, ins: &[Tok], r: &Result<Vec<TrajectoryPoint>, KimiyaError>) {
    let (outs, summary) = match r {
        Ok(tr) => (traj_toks(tr), traj_note(tr)),
        Err(e) => (vec![err(e)], format!("Err({})", format!("{e:?}").split('(').next().unwrap_or_default())),
    };
    assert!(outs.len() <= MAX_TOKENS && ins.len() <= MAX_TOKENS, "{func}: a row of {} tokens", outs.len());
    m.row_summarised(func, class, ins, &outs, &summary);
}

fn first(m: &mut Module, c: f64, k: f64, t: f64, n: usize) {
    let r = simulate_first_order(bb(c), bb(k), bb(t), bb(n));
    traj_row(m, "reaction_dynamics::simulate_first_order", Class::B, &[f(c), f(k), f(t), word(n as u64)], &r);
}

fn second(m: &mut Module, c: f64, k: f64, t: f64, n: usize) {
    let r = simulate_second_order(bb(c), bb(k), bb(t), bb(n));
    traj_row(m, "reaction_dynamics::simulate_second_order", Class::B, &[f(c), f(k), f(t), word(n as u64)], &r);
}

fn consecutive(m: &mut Module, a: f64, k1: f64, k2: f64, t: f64, n: usize) {
    let r = simulate_consecutive(bb(a), bb(k1), bb(k2), bb(t), bb(n));
    let ins = [f(a), f(k1), f(k2), f(t), word(n as u64)];
    traj_row(m, "reaction_dynamics::simulate_consecutive", Class::B, &ins, &r);
}

fn reversible(m: &mut Module, a: f64, b: f64, kf: f64, kr: f64, t: f64, n: usize) {
    let r = simulate_reversible(bb(a), bb(b), bb(kf), bb(kr), bb(t), bb(n));
    let ins = [f(a), f(b), f(kf), f(kr), f(t), word(n as u64)];
    traj_row(m, "reaction_dynamics::simulate_reversible", Class::B, &ins, &r);
}

fn mm(m: &mut Module, s0: f64, v: f64, km: f64, t: f64, n: usize) {
    let r = simulate_michaelis_menten(bb(s0), bb(v), bb(km), bb(t), bb(n));
    let ins = [f(s0), f(v), f(km), f(t), word(n as u64)];
    traj_row(m, "reaction_dynamics::simulate_michaelis_menten", Class::B, &ins, &r);
}

fn adaptive(m: &mut Module, sys: &str, p: &[f64], y0: &[f64], t: f64, tol: f64) {
    let fs = |t: f64, y: &[f64], dy: &mut [f64]| rhs(sys, p, t, y, dy);
    let r = simulate_adaptive(fs, bb(y0), bb(t), bb(tol));
    let mut ins = sys_toks(sys, p);
    ins.extend(vec_toks(y0));
    ins.extend([f(t), f(tol)]);
    traj_row(m, "reaction_dynamics::simulate_adaptive", Class::C, &ins, &r);
}

fn stiff(m: &mut Module, sys: &str, p: &[f64], y0: &[f64], t: f64, n: usize) {
    let fs = |t: f64, y: &[f64], dy: &mut [f64]| rhs(sys, p, t, y, dy);
    let js = |t: f64, y: &[f64], j: &mut Vec<Vec<f64>>| jac(sys, p, t, y, j);
    let r = simulate_stiff(fs, js, bb(y0), bb(t), bb(n));
    let mut ins = sys_toks(sys, p);
    ins.extend(vec_toks(y0));
    ins.extend([f(t), word(n as u64)]);
    traj_row(m, "reaction_dynamics::simulate_stiff", Class::B, &ins, &r);
}

#[allow(clippy::too_many_arguments)]
fn stochastic(m: &mut Module, dr: &str, p: &[f64], df: &str, q: &[f64], y0: &[f64], t: f64, n: usize, seed: u64) {
    let fd = |t: f64, y: &[f64], dy: &mut [f64]| rhs(dr, p, t, y, dy);
    let fg = |t: f64, y: &[f64], g: &mut [f64]| diffusion(df, q, t, y, g);
    let r = simulate_stochastic(fd, fg, bb(y0), bb(t), bb(n), bb(seed));
    let mut ins = sys_toks(dr, p);
    ins.extend(sys_toks(df, q));
    ins.extend(vec_toks(y0));
    ins.extend([f(t), word(n as u64), word(seed)]);
    traj_row(m, "reaction_dynamics::simulate_stochastic", Class::C, &ins, &r);
}

/// The steps edges: 0, 1, 2, 3, 10⁷ + 1 and usize::MAX are Err or tiny; 10⁷ itself is accepted (simulate_stiff's
/// empty system shows it without a 10⁷-point output).
const BAD_STEPS: [usize; 3] = [0, 10_000_001, usize::MAX];

/// The step counts every preset runs at (port plan §4: 10, 100 and 1000).
const PRESET_STEPS: [usize; 3] = [10, 100, 1000];

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("reaction_dynamics");
    let mut rng = Rng::new(0x5a01);

    // ── simulate_first_order ──
    // the Rust tests' and benches' values
    for (c, k, t, n) in [
        (1.0, 0.1, 10.0, 100),
        (1.0, 0.5, 5.0, 50),
        (1.0, 0.1, 10.0, 1000),
        (-1.0, 0.1, 10.0, 100),
        (1.0, 0.1, 0.0, 100),
        (1.0, 0.0, 10.0, 10),
        (1.0, 1e3, 10.0, 100), // k·h = 100: the unstable-RK4 Err
    ] {
        first(&mut m, c, k, t, n);
    }
    // mechanism presets at 10, 100, 1000 steps
    for (c, k, t) in [(1.0, 0.1, 10.0), (2.5, 1.0, 5.0), (1e-3, 50.0, 1.0), (0.0, 1.0, 1.0), (7.0, 0.0, 3.0)] {
        for n in PRESET_STEPS {
            first(&mut m, c, k, t, n);
        }
    }
    // each argument over the edges, the others typical
    for x in edges() {
        first(&mut m, x, 0.1, 10.0, 10);
        first(&mut m, 1.0, x, 10.0, 10);
        first(&mut m, 1.0, 0.1, x, 10);
    }
    for n in BAD_STEPS.into_iter().chain([1, 2, 3]) {
        first(&mut m, 1.0, 0.1, 10.0, n);
    }
    // RK4's stability boundary (k·h ≈ 2.785) and past it: finite growth, then overflow to the Err
    for kh in [2.0, 2.7, 2.78, 2.785, 2.79, 2.8, 3.0, 5.0, 10.0, 100.0] {
        first(&mut m, 1.0, kh, 10.0, 10);
        first(&mut m, 1.0, kh * 100.0, 10.0, 1000);
    }
    // extreme magnitudes: [A]₀ subnormal and near DBL_MAX, t_end vanishing (h underflows) or near DBL_MAX
    for (c, k, t, n) in [
        (f64::MAX, 0.1, 10.0, 10),
        (1e300, 1e-300, 1e300, 10),
        (f64::from_bits(1), 0.1, 10.0, 10),
        (1.0, 0.1, f64::from_bits(1), 10),
        (1.0, 0.1, f64::from_bits(1), 1),
        (1.0, 0.1, f64::MAX, 1),
        (1.0, 0.0, f64::MAX, 3),
        (1.0, f64::MAX, 1.0, 1),
        (1.0, 1e-320, 1e-10, 10),
    ] {
        first(&mut m, c, k, t, n);
    }
    for _ in 0..40 {
        let c = rng.log_uniform(-6.0, 3.0);
        let k = rng.log_uniform(-4.0, 1.0);
        let t = rng.log_uniform(-1.0, 2.0);
        let n = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        first(&mut m, c, k, t, n);
    }

    // ── simulate_second_order ──
    for (c, k, t, n) in [(1.0, 0.1, 10.0, 100), (1.0, 0.1, 10.0, 10), (-1.0, 0.1, 10.0, 10), (1.0, 1e3, 10.0, 100)] {
        second(&mut m, c, k, t, n);
    }
    for (c, k, t) in [(1.0, 0.1, 10.0), (2.0, 3.0, 2.0), (0.01, 100.0, 50.0), (5.0, 0.0, 1.0)] {
        for n in PRESET_STEPS {
            second(&mut m, c, k, t, n);
        }
    }
    for x in edges() {
        second(&mut m, x, 0.1, 10.0, 10);
        second(&mut m, 1.0, x, 10.0, 10);
        second(&mut m, 1.0, 0.1, x, 10);
    }
    for n in BAD_STEPS.into_iter().chain([1, 2]) {
        second(&mut m, 1.0, 0.1, 10.0, n);
    }
    for (c, k, t, n) in [(1e200, 1e200, 1.0, 10), (1e-200, 1e-200, 1e300, 10), (1.0, 1e10, 1.0, 10)] {
        second(&mut m, c, k, t, n);
    }
    for _ in 0..30 {
        let c = rng.log_uniform(-4.0, 2.0);
        let k = rng.log_uniform(-3.0, 1.0);
        let t = rng.log_uniform(-1.0, 2.0);
        let n = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        second(&mut m, c, k, t, n);
    }

    // ── simulate_consecutive ──
    for (a, k1, k2, t, n) in [
        (1.0, 0.5, 0.1, 20.0, 200),
        (1.0, 1.0, 0.1, 10.0, 100),
        (1.0, 0.5, 0.1, 20.0, 100), // the bench
    ] {
        consecutive(&mut m, a, k1, k2, t, n);
    }
    for (a, k1, k2, t) in [(1.0, 0.5, 0.1, 20.0), (1.0, 0.3, 0.3, 30.0), (2.0, 10.0, 1.0, 5.0), (1.0, 0.0, 1.0, 1.0)] {
        for n in PRESET_STEPS {
            consecutive(&mut m, a, k1, k2, t, n);
        }
    }
    for x in edges() {
        consecutive(&mut m, x, 0.1, 0.1, 10.0, 10);
        consecutive(&mut m, 1.0, x, 0.1, 10.0, 10);
        consecutive(&mut m, 1.0, 0.1, x, 10.0, 10);
        consecutive(&mut m, 1.0, 0.1, 0.1, x, 10);
    }
    for n in BAD_STEPS.into_iter().chain([1]) {
        consecutive(&mut m, 1.0, 0.1, 0.1, 10.0, n);
    }
    for _ in 0..30 {
        let a = rng.log_uniform(-3.0, 1.0);
        let k1 = rng.log_uniform(-2.0, 1.0);
        let k2 = rng.log_uniform(-2.0, 1.0);
        let t = rng.log_uniform(0.0, 1.5);
        let n = usize::try_from(rng.int(5, 150)).unwrap_or(5);
        consecutive(&mut m, a, k1, k2, t, n);
    }

    // ── simulate_reversible ──
    for (a, b, kf, kr, t, n) in [
        (1.0, 0.0, 1.0, 0.5, 20.0, 200),
        (0.8, 0.2, 0.5, 0.3, 10.0, 100),
        (1.0, 0.0, 2.0, 1.0, 10.0, 1000), // the integration test
    ] {
        reversible(&mut m, a, b, kf, kr, t, n);
    }
    for (a, b, kf, kr, t) in [(1.0, 0.0, 1.0, 0.5, 20.0), (0.0, 1.0, 0.2, 2.0, 10.0), (0.5, 0.5, 0.0, 0.0, 1.0)] {
        for n in PRESET_STEPS {
            reversible(&mut m, a, b, kf, kr, t, n);
        }
    }
    for x in edges() {
        reversible(&mut m, x, 0.0, 0.1, 0.1, 10.0, 10);
        reversible(&mut m, 1.0, x, 0.1, 0.1, 10.0, 10);
        reversible(&mut m, 1.0, 0.0, x, 0.1, 10.0, 10);
        reversible(&mut m, 1.0, 0.0, 0.1, x, 10.0, 10);
        reversible(&mut m, 1.0, 0.0, 0.1, 0.1, x, 10);
    }
    for n in BAD_STEPS.into_iter().chain([1]) {
        reversible(&mut m, 1.0, 0.0, 0.1, 0.1, 10.0, n);
    }
    for _ in 0..30 {
        let a = rng.log_uniform(-3.0, 1.0);
        let b = rng.log_uniform(-3.0, 1.0);
        let kf = rng.log_uniform(-2.0, 1.0);
        let kr = rng.log_uniform(-2.0, 1.0);
        let t = rng.log_uniform(0.0, 1.5);
        let n = usize::try_from(rng.int(5, 150)).unwrap_or(5);
        reversible(&mut m, a, b, kf, kr, t, n);
    }

    // ── simulate_michaelis_menten ──
    for (s0, v, km, t, n) in [(10.0, 1.0, 2.0, 30.0, 300), (5.0, 2.0, 1.0, 10.0, 100)] {
        mm(&mut m, s0, v, km, t, n);
    }
    for (s0, v, km, t) in [(10.0, 1.0, 2.0, 30.0), (1.0, 1.0, 1e-6, 2.0), (1e-3, 1.0, 1e3, 1e3), (0.0, 1.0, 1.0, 1.0)] {
        for n in PRESET_STEPS {
            mm(&mut m, s0, v, km, t, n);
        }
    }
    for x in edges() {
        mm(&mut m, x, 1.0, 1.0, 10.0, 10);
        mm(&mut m, 1.0, x, 1.0, 10.0, 10);
        mm(&mut m, 1.0, 1.0, x, 10.0, 10);
        mm(&mut m, 1.0, 1.0, 1.0, x, 10);
    }
    for n in BAD_STEPS.into_iter().chain([1]) {
        mm(&mut m, 1.0, 1.0, 1.0, 10.0, n);
    }
    for _ in 0..30 {
        let s0 = rng.log_uniform(-3.0, 2.0);
        let v = rng.log_uniform(-2.0, 1.0);
        let km = rng.log_uniform(-3.0, 2.0);
        let t = rng.log_uniform(0.0, 2.0);
        let n = usize::try_from(rng.int(5, 150)).unwrap_or(5);
        mm(&mut m, s0, v, km, t, n);
    }

    // ── simulate_adaptive ──
    adaptive(&mut m, "decay", &[0.1], &[1.0], 10.0, 1e-8); // the Rust test
    for tol in [0.0, -1e-8, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MIN_POSITIVE, f64::MAX] {
        adaptive(&mut m, "decay", &[1.0], &[1.0], 1.0, tol);
    }
    for y in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        adaptive(&mut m, "decay", &[1.0], &[y], 1.0, 1e-8);
        adaptive(&mut m, "reversible", &[1.0, 0.5], &[1.0, y], 1.0, 1e-8);
    }
    for t in edges() {
        adaptive(&mut m, "zero", &[], &[1.0], t, 1e-8);
    }
    adaptive(&mut m, "zero", &[], &[], 1.0, 1e-8); // the empty state: hisab's InvalidInput
    adaptive(&mut m, "decay", &[1.0], &[1.0], f64::from_bits(1), 1e-8); // h_init = t_end/100 = 0
    adaptive(&mut m, "decay", &[1.0], &[1.0], 1e-300, 1e-8);
    adaptive(&mut m, "decay", &[1.0], &[1.0], 1e300, 1e-6);
    adaptive(&mut m, "zero", &[], &[1.0, 2.0], 100.0, 1e-8); // dy never written: zero error, ×5 steps
    adaptive(&mut m, "nan", &[], &[1.0], 1.0, 1e-8); // never accepted: the 100 000-step limit
    adaptive(&mut m, "blowup", &[1e3], &[1.0], 1.0, 1e-6); // the solution blows up at t = 10⁻³
    for tol in [1e-3, 1e-6, 1e-9, 1e-12] {
        adaptive(&mut m, "decay", &[0.1], &[1.0], 10.0, tol);
        adaptive(&mut m, "decay2", &[3.0], &[2.0], 2.0, tol);
        adaptive(&mut m, "consecutive", &[0.5, 0.1], &[1.0, 0.0, 0.0], 20.0, tol);
        adaptive(&mut m, "reversible", &[1.0, 0.5], &[1.0, 0.0], 20.0, tol);
        adaptive(&mut m, "michaelis_menten", &[1.0, 2.0], &[10.0, 0.0], 30.0, tol);
        adaptive(&mut m, "brusselator", &[1.0, 3.0], &[1.5, 3.0], 20.0, tol);
        adaptive(&mut m, "forced", &[2.0, 1.0], &[1.0], 5.0, tol);
    }
    adaptive(&mut m, "decay", &[0.1, 0.2], &[1.0, 2.0, 3.0, 4.0], 10.0, 1e-8);
    adaptive(&mut m, "decay", &[-0.5], &[1e-300], 100.0, 1e-8); // growth from a tiny state
    adaptive(&mut m, "robertson", &[], &[1.0, 0.0, 0.0], 0.1, 1e-4); // stiff: small steps throughout
    for _ in 0..20 {
        let k = rng.log_uniform(-2.0, 1.0);
        let y0 = rng.log_uniform(-3.0, 2.0);
        let t = rng.log_uniform(-1.0, 1.5);
        let tol = rng.log_uniform(-10.0, -3.0);
        adaptive(&mut m, "decay", &[k], &[y0], t, tol);
    }

    // ── simulate_stiff ──
    stiff(&mut m, "decay", &[50.0], &[1.0], 1.0, 100); // the Rust tests
    stiff(&mut m, "decay2", &[3.0], &[2.0], 0.05 * 40.0, 40);
    stiff(&mut m, "robertson", &[], &[1.0, 0.0, 0.0], 40.0, 400);
    stiff(&mut m, "affine", &[10.0], &[1.0], 1.0, 10); // I − h·J singular
    stiff(&mut m, "badjac", &[], &[1.0], 1.0, 1); // Newton diverges
    stiff(&mut m, "decay", &[1.0], &[f64::NAN], 1.0, 1);
    stiff(&mut m, "nan", &[], &[1.0], 1.0, 1);
    for t in [f64::NAN, f64::INFINITY, -1.0, 0.0, -0.0] {
        stiff(&mut m, "zero", &[], &[1.0], t, 10);
    }
    for n in BAD_STEPS {
        stiff(&mut m, "zero", &[], &[1.0], 1.0, n);
        stiff(&mut m, "zero", &[], &[], 1.0, n);
    }
    stiff(&mut m, "zero", &[], &[], 1.0, 10_000_000); // the step limit itself: nothing to integrate
    stiff(&mut m, "zero", &[], &[1.0, -2.0], 1.0, 5);
    for y in [f64::INFINITY, f64::NEG_INFINITY] {
        stiff(&mut m, "decay", &[1.0], &[y], 1.0, 1);
    }
    for n in PRESET_STEPS {
        stiff(&mut m, "decay", &[50.0], &[1.0], 1.0, n);
        stiff(&mut m, "decay2", &[3.0], &[2.0], 2.0, n);
        stiff(&mut m, "consecutive", &[0.5, 0.1], &[1.0, 0.0, 0.0], 20.0, n);
        stiff(&mut m, "reversible", &[1.0, 0.5], &[1.0, 0.0], 20.0, n);
        stiff(&mut m, "michaelis_menten", &[1.0, 2.0], &[10.0, 0.0], 30.0, n);
        stiff(&mut m, "brusselator", &[1.0, 3.0], &[1.5, 3.0], 20.0, n);
        stiff(&mut m, "forced", &[2.0, 1.0], &[1.0], 5.0, n);
        stiff(&mut m, "robertson", &[], &[1.0, 0.0, 0.0], 40.0, n);
    }
    for n in [1, 2, 4, 40, 4000] {
        stiff(&mut m, "robertson", &[], &[1.0, 0.0, 0.0], 40.0, n);
    }
    stiff(&mut m, "robertson", &[], &[1.0, 0.0, 0.0], 1e5, 100);
    stiff(&mut m, "affine", &[-10.0], &[1.0], 1.0, 10);
    stiff(&mut m, "blowup", &[1e3], &[1.0], 1.0, 100);
    // a Newton that converges linearly, at r = |1 − h·k|/2 per iteration (h·c = 1): converged within 20 iterations,
    // accepted at the 1e-6 floor, or Err, across the sweep
    for k10 in 0..=40 {
        let k = f64::from(k10) * 0.5 - 5.0;
        stiff(&mut m, "slowjac", &[k, 10.0], &[1.0], 0.1, 1);
        stiff(&mut m, "slowjac", &[k, 10.0], &[1.0], 1.0, 10);
    }
    for _ in 0..20 {
        let k = rng.log_uniform(-1.0, 4.0);
        let y0 = rng.log_uniform(-3.0, 2.0);
        let t = rng.log_uniform(-1.0, 1.0);
        let n = usize::try_from(rng.int(1, 300)).unwrap_or(1);
        stiff(&mut m, "decay2", &[k], &[y0], t, n);
    }

    // ── simulate_stochastic ──
    stochastic(&mut m, "decay", &[1.0], "zero", &[], &[1.0], 1.0, 10, 7); // the Rust tests
    stochastic(&mut m, "decay", &[1.0], "const", &[0.1], &[1.0], 1.0, 50, 42);
    stochastic(&mut m, "blowup", &[1e3], "zero", &[], &[1.0], 1.0, 100, 1);
    stochastic(&mut m, "blowup", &[1e3], "zero", &[], &[f64::NAN], 1.0, 10, 1);
    for t in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        stochastic(&mut m, "zero", &[], "none", &[], &[1.0], t, 10, 1);
    }
    for n in BAD_STEPS {
        stochastic(&mut m, "zero", &[], "none", &[], &[1.0], 1.0, n, 1);
    }
    stochastic(&mut m, "zero", &[], "none", &[], &[], 1.0, 10, 1); // the empty state
    stochastic(&mut m, "zero", &[], "none", &[], &[1.0, 2.0], 1.0, 10, 3); // nothing written: y stays
    for seed in [0, 1, 7, 42, 1 << 63, u64::MAX] {
        for n in PRESET_STEPS {
            stochastic(&mut m, "decay", &[1.0], "const", &[0.1], &[1.0], 1.0, n, seed);
        }
        stochastic(&mut m, "reversible", &[1.0, 0.5], "const", &[0.05], &[1.0, 0.0], 5.0, 100, seed);
        stochastic(&mut m, "consecutive", &[0.5, 0.1], "mult", &[0.2], &[1.0, 0.0, 0.0], 10.0, 100, seed);
    }
    for n in PRESET_STEPS {
        stochastic(&mut m, "brusselator", &[1.0, 3.0], "const", &[0.1], &[1.5, 3.0], 10.0, n, 2024);
        stochastic(&mut m, "michaelis_menten", &[1.0, 2.0], "mult", &[0.1], &[10.0, 0.0], 30.0, n, 5);
    }
    stochastic(&mut m, "decay", &[1.0], "const", &[1e300], &[1.0], 1.0, 10, 1); // the noise overflows the state
    stochastic(&mut m, "decay", &[1.0], "const", &[0.1], &[1.0], f64::from_bits(1), 10, 1);
    stochastic(&mut m, "decay", &[1.0], "const", &[0.1], &[1.0], f64::MAX, 2, 1);
    for _ in 0..20 {
        let seed = u64::try_from(rng.int(0, 1 << 62)).unwrap_or(0) * 4 + 3;
        let k = rng.log_uniform(-2.0, 1.0);
        let sg = rng.log_uniform(-3.0, 0.0);
        let n = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        stochastic(&mut m, "decay", &[k], "mult", &[sg], &[1.0], 2.0, n, seed);
    }
    m
}
