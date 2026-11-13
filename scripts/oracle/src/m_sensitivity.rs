//! `sensitivity.tsv` — gradients of chemical quantities (port plan §4, M5b), all class D (relative 10⁻¹², port plan §5:
//! autodiff gradients).
//!
//! - The four fixed-formula gradients, `equilibrium_constant_gradient` (input `f:ΔG f:T`), `arrhenius_gradient`
//!   (`f:A f:Ea f:T`), `nernst_gradient` (`f:E° i:n f:T f:Q`, n a u8) and `gibbs_temperature_sensitivity` (`f:ΔH f:T
//!   f:ΔS`); output the oracle's `Vec<f64>` (`v:n` + `f:`) or `f64`, hisab 1.4.0's reverse sweep of each function's tape.
//!   **ADR 0002 #2**: the port returns `Result`s from closed-form derivatives (the same sweep unrolled), validated as
//!   `reaction::equilibrium_constant`, `kinetics::arrhenius_rate`, `electrochemistry::nernst_potential` and
//!   `reaction::gibbs_free_energy` validate (codes and order), with Err(ComputationError) where an entry is not finite.
//!   A row where the port returns Err is noted with it, `[ADR 0002 #2: the port's Err(V) (E:V)]` (the oracle's own
//!   output, finite or not, stays the row's): every non-finite oracle row, and the finite ones at inputs the port
//!   refuses (T ≤ 0 with finite entries, n = 0 …). The grids: the Rust tests' values, each argument over the edges
//!   around a physical point, the exponent −ΔG/(RT) across exp's overflow and underflow, Q across ln's range, n over
//!   every u8 class, seeded points.
//! - `compute_gradient` on the test functions of [`build`] (input `s:name v:k f:constants… v:n f:parameters…`; output
//!   the gradient `v:n` + `f:`): the Rust test's x² + 3y, the four fixed formulas' tapes, van 't Hoff through ΔG = ΔH −
//!   TΔS, Eyring, the ideal gas, Clausius–Clapeyron, Lennard-Jones (powf), Morse, Rosenbrock, sin·cos + x/y, ln sums,
//!   sums of squares and products of up to 16 parameters, a node reached along several paths, an unused parameter, the
//!   identity and a constant of no parameters; at physical and seeded points. The port's compute_gradient builds on
//!   hisab 3.4.1's tape (`ad_*`), so these rows compare two autodiff implementations at class D.
//!
//! Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::grids::{Rng, pow2};
use crate::systems::{sys_toks, vec_toks};
use crate::tsv::{Class, Module, Tok, f, i, tokens_note};
use hisab::autodiff::{Tape, Var};
use kimiya::electrochemistry::{FARADAY, nernst_potential};
use kimiya::element::GAS_CONSTANT;
use kimiya::error::KimiyaError;
use kimiya::kinetics::{BOLTZMANN, arrhenius_rate};
use kimiya::reaction::{equilibrium_constant, gibbs_free_energy};
use kimiya::sensitivity::{
    arrhenius_gradient, compute_gradient, equilibrium_constant_gradient, gibbs_temperature_sensitivity, nernst_gradient,
};
use kimiya::spectroscopy::PLANCK;

/// The port's Err for a row (ADR 0002 #2): the validating function's Err (its code and order), else ComputationError
/// where an entry is not finite; `None` when the port returns the gradient.
fn port_err<T>(v: Result<T, KimiyaError>, grads: &[f64]) -> Option<String> {
    let finite = grads.iter().all(|x| x.is_finite());
    match v {
        Err(KimiyaError::ComputationError(_)) => {
            assert!(!finite, "a ComputationError of the validator with a finite gradient");
            Some("ComputationError".to_string())
        }
        Err(e) => Some(format!("{e:?}").split('(').next().unwrap_or_default().to_string()),
        Ok(_) if finite => None,
        Ok(_) => Some("ComputationError".to_string()),
    }
}

fn noted(m: &mut Module, func: &str, ins: &[Tok], outs: &[Tok], port: Option<String>) {
    match port {
        None => m.row(func, Class::D, ins, outs),
        Some(e) => {
            let extra = format!("ADR 0002 #2: the port's {}", tokens_note(&[Tok::Err(e)]));
            m.row_noted(func, Class::D, ins, outs, &extra);
        }
    }
}

fn keq(m: &mut Module, dg: f64, t: f64) {
    let g = equilibrium_constant_gradient(bb(dg), bb(t));
    let port = port_err(equilibrium_constant(dg, t), &g);
    noted(m, "sensitivity::equilibrium_constant_gradient", &[f(dg), f(t)], &vec_toks(&g), port);
}

fn arr(m: &mut Module, a: f64, ea: f64, t: f64) {
    let g = arrhenius_gradient(bb(a), bb(ea), bb(t));
    let port = port_err(arrhenius_rate(a, ea, t), &g);
    noted(m, "sensitivity::arrhenius_gradient", &[f(a), f(ea), f(t)], &vec_toks(&g), port);
}

fn nernst(m: &mut Module, e0: f64, n: u8, t: f64, q: f64) {
    let g = nernst_gradient(bb(e0), bb(n), bb(t), bb(q));
    let port = port_err(nernst_potential(e0, n, t, q), &g);
    noted(m, "sensitivity::nernst_gradient", &[f(e0), i(n), f(t), f(q)], &vec_toks(&g), port);
}

fn gibbs(m: &mut Module, h: f64, t: f64, s: f64) {
    let g = gibbs_temperature_sensitivity(bb(h), bb(t), bb(s));
    let port = port_err(gibbs_free_energy(h, t, s), &[g]);
    noted(m, "sensitivity::gibbs_temperature_sensitivity", &[f(h), f(t), f(s)], &[f(g)], port);
}

/// The test functions of the `compute_gradient` rows: `name` with constants `p` on `tape` over the parameters `v`. The
/// parity harness (`tests/sens_functions.cyr`, on hisab 3.4.1's tape) and `validate.py` (on its replica of hisab
/// 1.4.0's) build the same operations in the same order; the formulas are in `scripts/oracle/README.md`.
pub fn build(name: &str, p: &[f64], tape: &mut Tape, v: &[Var]) -> Var {
    match name {
        // x² + 3y (the Rust test custom_gradient)
        "quad" => {
            let x2 = tape.mul(v[0], v[0]);
            let three = tape.constant(3.0);
            let ty = tape.mul(three, v[1]);
            tape.add(x2, ty)
        }
        // K = exp(−ΔG/(R·T)), parameters [ΔG, T], p = [R] (equilibrium_constant_gradient's tape)
        "keq" => {
            let r = tape.constant(p[0]);
            let n = tape.neg(v[0]);
            let rt = tape.mul(r, v[1]);
            let e = tape.div(n, rt);
            tape.exp(e)
        }
        // k = A·exp(−Ea/(R·T)), [A, Ea, T], p = [R] (arrhenius_gradient's tape)
        "arrhenius" => {
            let r = tape.constant(p[0]);
            let n = tape.neg(v[1]);
            let rt = tape.mul(r, v[2]);
            let e = tape.div(n, rt);
            let b = tape.exp(e);
            tape.mul(v[0], b)
        }
        // E = E° − (R·T/nF)·ln Q, [E°, T, Q], p = [R, nF] (nernst_gradient's tape)
        "nernst" => {
            let r = tape.constant(p[0]);
            let nf = tape.constant(p[1]);
            let rt = tape.mul(r, v[1]);
            let lq = tape.ln(v[2]);
            let q = tape.div(rt, nf);
            let c = tape.mul(q, lq);
            tape.sub(v[0], c)
        }
        // ΔG = ΔH − T·ΔS, [ΔH, T, ΔS] (gibbs_temperature_sensitivity's tape)
        "gibbs" => {
            let ts = tape.mul(v[1], v[2]);
            tape.sub(v[0], ts)
        }
        // K = exp(−(ΔH − T·ΔS)/(R·T)), [ΔH, ΔS, T], p = [R]: van 't Hoff, as the module doc suggests
        "vant_hoff" => {
            let ts = tape.mul(v[2], v[1]);
            let dg = tape.sub(v[0], ts);
            let n = tape.neg(dg);
            let r = tape.constant(p[0]);
            let rt = tape.mul(r, v[2]);
            let e = tape.div(n, rt);
            tape.exp(e)
        }
        // k = (k_B·T/h)·exp(−ΔG‡/(R·T)), [ΔG‡, T], p = [k_B, h, R]
        "eyring" => {
            let kb = tape.constant(p[0]);
            let h = tape.constant(p[1]);
            let kbt = tape.mul(kb, v[1]);
            let pre = tape.div(kbt, h);
            let r = tape.constant(p[2]);
            let rt = tape.mul(r, v[1]);
            let n = tape.neg(v[0]);
            let e = tape.div(n, rt);
            let b = tape.exp(e);
            tape.mul(pre, b)
        }
        // P = n·R·T/V, [n, T, V], p = [R]
        "ideal_gas" => {
            let r = tape.constant(p[0]);
            let nr = tape.mul(v[0], r);
            let nrt = tape.mul(nr, v[1]);
            tape.div(nrt, v[2])
        }
        // p₂ = p₁·exp(−(ΔH/R)·(1/T₂ − 1/T₁)), [p₁, ΔH, T₁, T₂], p = [R]
        "clausius" => {
            let one = tape.constant(1.0);
            let i2 = tape.div(one, v[3]);
            let i1 = tape.div(one, v[2]);
            let d = tape.sub(i2, i1);
            let hd = tape.mul(v[1], d);
            let r = tape.constant(p[0]);
            let x = tape.div(hd, r);
            let nx = tape.neg(x);
            let e = tape.exp(nx);
            tape.mul(v[0], e)
        }
        // V = 4ε((σ/r)¹² − (σ/r)⁶), [ε, σ, r]
        "lj" => {
            let q = tape.div(v[1], v[2]);
            let q12 = tape.powf(q, 12.0);
            let q6 = tape.powf(q, 6.0);
            let d = tape.sub(q12, q6);
            let four = tape.constant(4.0);
            let fe = tape.mul(four, v[0]);
            tape.mul(fe, d)
        }
        // V = D_e·(1 − exp(−a·(r − r_e)))², [D_e, a, r, r_e]
        "morse" => {
            let dr = tape.sub(v[2], v[3]);
            let adr = tape.mul(v[1], dr);
            let n = tape.neg(adr);
            let e = tape.exp(n);
            let one = tape.constant(1.0);
            let om = tape.sub(one, e);
            let sq = tape.mul(om, om);
            tape.mul(v[0], sq)
        }
        // (1 − x)² + 100·(y − x²)², [x, y]
        "rosenbrock" => {
            let one = tape.constant(1.0);
            let a = tape.sub(one, v[0]);
            let a2 = tape.mul(a, a);
            let x2 = tape.mul(v[0], v[0]);
            let b = tape.sub(v[1], x2);
            let b2 = tape.mul(b, b);
            let h = tape.constant(100.0);
            let hb = tape.mul(h, b2);
            tape.add(a2, hb)
        }
        // sin(x)·cos(y) + x/y, [x, y]
        "trig" => {
            let s = tape.sin(v[0]);
            let c = tape.cos(v[1]);
            let sc = tape.mul(s, c);
            let q = tape.div(v[0], v[1]);
            tape.add(sc, q)
        }
        // ln x + z·ln y, [x, y, z]
        "logs" => {
            let lx = tape.ln(v[0]);
            let ly = tape.ln(v[1]);
            let zl = tape.mul(v[2], ly);
            tape.add(lx, zl)
        }
        // Σ x_k², [x_0 … x_{n−1}], n ≥ 1
        "sum_squares" => {
            let mut acc = tape.mul(v[0], v[0]);
            for &x in &v[1..] {
                let sq = tape.mul(x, x);
                acc = tape.add(acc, sq);
            }
            acc
        }
        // Π x_k, n ≥ 1
        "product" => {
            let mut acc = v[0];
            for &x in &v[1..] {
                acc = tape.mul(acc, x);
            }
            acc
        }
        // (x² + x)·x² − x: x reached along several paths, [x]
        "shared" => {
            let a = tape.mul(v[0], v[0]);
            let b = tape.add(a, v[0]);
            let c = tape.mul(b, a);
            tape.sub(c, v[0])
        }
        // x², [x, y]: y unused
        "unused" => tape.mul(v[0], v[0]),
        // x, [x]
        "ident" => v[0],
        // p[0], no parameters
        "const" => tape.constant(p[0]),
        _ => panic!("no gradient function {name}"),
    }
}

fn grad(m: &mut Module, name: &str, p: &[f64], x: &[f64]) {
    let g = compute_gradient(|tape, v| build(name, p, tape, v), bb(x)).unwrap_or_default();
    assert_eq!(g.len(), x.len(), "compute_gradient: one entry per parameter");
    let mut ins = sys_toks(name, p);
    ins.extend(vec_toks(x));
    m.row("sensitivity::compute_gradient", Class::D, &ins, &vec_toks(&g));
}

/// The temperatures every fixed formula is probed at: physical ones, the edges the validation meets.
fn temps() -> Vec<f64> {
    vec![
        298.15,
        298.0,
        1.0,
        1000.0,
        5000.0,
        1e-3,
        1e-150,
        1e-300,
        pow2(-1074),
        f64::MIN_POSITIVE,
        1e150,
        1e300,
        f64::MAX,
        0.0,
        -0.0,
        -1.0,
        -298.15,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ]
}

/// Signed magnitudes for the energies (J/mol) and the edges.
fn energies() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        1.0,
        -1.0,
        1e3,
        -5e3,
        8e4,
        -1e5,
        1e6,
        -1e7,
        1e300,
        -1e300,
        f64::MAX,
        -f64::MAX,
        pow2(-1074),
        -pow2(-1074),
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ]
}

#[allow(clippy::too_many_lines)] // one function after another
pub fn emit() -> Module {
    let mut m = Module::new("sensitivity");
    let mut rng = Rng::new(0x5b03);
    let r = GAS_CONSTANT;

    // ── equilibrium_constant_gradient ──
    keq(&mut m, -5000.0, 298.0); // the Rust tests
    keq(&mut m, -5000.0, 0.0);
    keq(&mut m, f64::NAN, 298.0);
    for t in temps() {
        keq(&mut m, -5000.0, t);
        keq(&mut m, 5000.0, t);
    }
    for dg in energies() {
        keq(&mut m, dg, 298.15);
        keq(&mut m, dg, 1e-300);
    }
    // −ΔG/(RT) across exp's overflow (709.78) and underflow (−745.13) ends, and where K is finite but K/(RT) or
    // K·ΔG/(RT)² overflows (large K, small T) or vanishes.
    for x in [700.0, 709.0, 709.7, 709.78, 709.79, 710.0, -700.0, -708.4, -740.0, -745.0, -745.13, -745.2, -800.0] {
        for t in [298.15, 1e-3, 1e-200, 1e200] {
            keq(&mut m, -x * (r * t), t);
        }
    }
    for _ in 0..80 {
        let dg = rng.signed_log_uniform(0.0, 7.0);
        let t = rng.log_uniform(0.0, 4.0);
        keq(&mut m, dg, t);
    }

    // ── arrhenius_gradient ──
    arr(&mut m, 1e13, 80_000.0, 300.0); // the Rust tests
    arr(&mut m, 1e13, 80_000.0, f64::NAN);
    for t in temps() {
        arr(&mut m, 1e13, 80_000.0, t);
        arr(&mut m, 0.0, 80_000.0, t);
        arr(&mut m, 1e13, -1e4, t);
    }
    for a in [1e13, 1.0, 0.0, -0.0, -1e5, 1e300, f64::MAX, pow2(-1074), f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        arr(&mut m, a, 80_000.0, 300.0);
        arr(&mut m, a, 0.0, 300.0);
        arr(&mut m, a, 1e6, 1e-300);
    }
    for ea in energies() {
        arr(&mut m, 1e13, ea, 300.0);
        arr(&mut m, 1.0, ea, 1e-170);
    }
    for x in [700.0, 709.78, 709.79, -745.13, -745.2] {
        arr(&mut m, 1.0, -x * (r * 300.0), 300.0);
        arr(&mut m, 1e-300, -x * (r * 300.0), 300.0);
        arr(&mut m, 1e300, -x * (r * 300.0), 300.0);
    }
    for _ in 0..80 {
        let a = rng.log_uniform(-3.0, 16.0);
        let ea = rng.signed_log_uniform(2.0, 6.0);
        let t = rng.log_uniform(1.0, 4.0);
        arr(&mut m, a, ea, t);
    }

    // ── nernst_gradient ──
    for (e0, n, t, q) in
        [(0.342, 0_u8, 298.0, 0.1), (0.342, 2, 298.0, 0.0), (0.342, 2, 298.0, -1.0), (0.342, 2, 298.0, 1.0), (0.342, 2, 298.0, 0.1)]
    {
        nernst(&mut m, e0, n, t, q); // the Rust tests
    }
    for n in [0_u8, 1, 2, 3, 4, 6, 8, 127, 128, 254, 255] {
        nernst(&mut m, 0.342, n, 298.15, 0.1);
        nernst(&mut m, -0.76, n, 1e300, 1e-300);
    }
    for t in temps() {
        nernst(&mut m, 0.342, 2, t, 0.1);
        nernst(&mut m, 0.342, 2, t, 1.0);
        nernst(&mut m, 0.342, 1, t, pow2(-1074));
    }
    for q in [pow2(-1074), f64::MIN_POSITIVE, 1e-300, 1e-10, 0.1, 1.0, 1.0_f64.next_up(), 10.0, 1e10, 1e300, f64::MAX, 0.0,
              -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        nernst(&mut m, 0.342, 2, 298.15, q);
        nernst(&mut m, 0.342, 2, 1e300, q);
    }
    for e0 in [0.342, -0.76, 0.0, -0.0, 1e300, -f64::MAX, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        nernst(&mut m, e0, 2, 298.15, 0.1);
    }
    for _ in 0..80 {
        let e0 = rng.uniform(-3.0, 3.0);
        let n = u8::try_from(rng.int(1, 8)).unwrap_or(1);
        let t = rng.log_uniform(1.0, 4.0);
        let q = rng.log_uniform(-20.0, 20.0);
        nernst(&mut m, e0, n, t, q);
    }

    // ── gibbs_temperature_sensitivity ──
    gibbs(&mut m, -100_000.0, 298.0, -200.0); // the Rust test
    for x in energies() {
        gibbs(&mut m, x, 298.0, -200.0);
        gibbs(&mut m, -100_000.0, 298.0, x);
    }
    for t in temps() {
        gibbs(&mut m, -100_000.0, t, -200.0);
        gibbs(&mut m, -100_000.0, t, 0.0);
    }
    for _ in 0..40 {
        let h = rng.signed_log_uniform(2.0, 7.0);
        let t = rng.log_uniform(1.0, 4.0);
        let s = rng.signed_log_uniform(-1.0, 3.0);
        gibbs(&mut m, h, t, s);
    }

    // ── compute_gradient ──
    let nf = 2.0 * FARADAY;
    grad(&mut m, "quad", &[], &[2.0, 5.0]); // the Rust test
    for &(dg, t) in &[(-5000.0, 298.0), (5000.0, 298.15), (-1e5, 1000.0), (2e4, 500.0)] {
        grad(&mut m, "keq", &[r], &[dg, t]);
    }
    for &(a, ea, t) in &[(1e13, 80_000.0, 300.0), (1e10, 50_000.0, 350.0), (1.0, -1e4, 1000.0)] {
        grad(&mut m, "arrhenius", &[r], &[a, ea, t]);
    }
    for &(e0, t, q) in &[(0.342, 298.0, 0.1), (0.342, 298.15, 10.0), (-0.76, 310.0, 1e-5)] {
        grad(&mut m, "nernst", &[r, nf], &[e0, t, q]);
    }
    grad(&mut m, "gibbs", &[], &[-100_000.0, 298.0, -200.0]);
    for &(dh, ds, t) in &[(-92_000.0, -199.0, 298.15), (57_200.0, 175.8, 400.0), (-41_000.0, -42.0, 700.0)] {
        grad(&mut m, "vant_hoff", &[r], &[dh, ds, t]);
    }
    for &(dg, t) in &[(80_000.0, 298.15), (50_000.0, 400.0)] {
        grad(&mut m, "eyring", &[BOLTZMANN, PLANCK, r], &[dg, t]);
    }
    grad(&mut m, "ideal_gas", &[r], &[1.0, 273.15, 0.0224]);
    grad(&mut m, "clausius", &[r], &[101_325.0, 40_700.0, 373.15, 363.15]);
    for &(eps, sigma, rr) in &[(1.65e-21, 3.4e-10, 3.8e-10), (1.0, 1.0, 1.122), (1.0, 1.0, 0.9)] {
        grad(&mut m, "lj", &[], &[eps, sigma, rr]);
    }
    grad(&mut m, "morse", &[], &[7.4e-19, 1.9e10, 1.5e-10, 1.27e-10]);
    grad(&mut m, "morse", &[], &[1.0, 1.0, 2.0, 1.0]);
    grad(&mut m, "rosenbrock", &[], &[-1.2, 1.0]);
    grad(&mut m, "rosenbrock", &[], &[1.0, 1.0]);
    grad(&mut m, "trig", &[], &[0.5, 1.5]);
    grad(&mut m, "trig", &[], &[-2.0, 3.0]);
    grad(&mut m, "logs", &[], &[2.0, 3.0, 0.5]);
    grad(&mut m, "shared", &[], &[1.5]);
    grad(&mut m, "shared", &[], &[-0.25]);
    grad(&mut m, "unused", &[], &[3.0, 7.0]);
    grad(&mut m, "ident", &[], &[4.0]);
    grad(&mut m, "const", &[5.0], &[]);
    for n in [1_usize, 2, 3, 8, 16] {
        let x: Vec<f64> = (0..n).map(|_| rng.uniform(-3.0, 3.0)).collect();
        grad(&mut m, "sum_squares", &[], &x);
        grad(&mut m, "product", &[], &x);
    }
    for _ in 0..20 {
        let dg = rng.signed_log_uniform(2.0, 5.0);
        let t = rng.log_uniform(2.0, 3.5);
        grad(&mut m, "keq", &[r], &[dg, t]);
        let x = rng.uniform(-2.0, 2.0);
        let y = rng.uniform(-1.0, 3.0);
        grad(&mut m, "rosenbrock", &[], &[x, y]);
        let lx = rng.log_uniform(-3.0, 3.0);
        let ly = rng.log_uniform(-3.0, 3.0);
        let z = rng.uniform(-2.0, 2.0);
        grad(&mut m, "logs", &[], &[lx, ly, z]);
        let tx = rng.uniform(-3.0, 3.0);
        let ty = rng.uniform(0.5, 3.0);
        grad(&mut m, "trig", &[], &[tx, ty]);
    }
    m
}
