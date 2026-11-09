//! `gas.tsv` — ideal gas, van der Waals, Dalton, gas density, Peng–Robinson (pressure, Z, its
//! roots, fugacity), Redlich–Kwong and the `PR_*` presets (port plan §4, M2).
//!
//! A `PengRobinsonParams` input is flattened as `f:tc f:pc f:omega`.
//! - Class A: `gas::PR_<NAME>` (no input; the preset's three fields), `joule_thomson_ideal`.
//! - Class B (+ − × ÷ and `sqrt` / `powi` only, expected 0 ulp): `ideal_gas_pressure`,
//!   `ideal_gas_volume`, `van_der_waals_pressure` (also at V = n·b and one ulp above),
//!   `partial_pressure`, `gas_density`, `peng_robinson_pressure` (at V = b, one ulp above, and
//!   b·(1.001 … 10⁶)). Each argument over the validation edges, the Rust tests' values, a seeded sweep.
//! - Class C: `peng_robinson_z`, `peng_robinson_z_roots` (`v:n` then the n ascending roots),
//!   `fugacity_coefficient_pr`, `fugacity_pr` and `redlich_kwong_pressure` (`tc.powf(2.5)`). The
//!   cubic's three-root (trigonometric) branch goes through `acos` and `cos`; its one-root branch
//!   is `cbrt` and `sqrt`, both correctly rounded, so those rows are expected at 0 ulp.
//!
//! The Peng–Robinson state grid, for each preset: T = T_r·Tc over 20 reduced temperatures from 0.3
//! to 10 (the near-critical 0.999, 1, 1.001 included) and the Rust tests' 200, 273.15, 280, 300,
//! 304.1 and 1000 K; P = P_r·Pc over 20 reduced pressures from 10⁻⁶ to 60 and the absolute
//! 1 kPa … 200 MPa of the Rust tests (101 325 Pa, 4 MPa, 7.375 MPa, 90 and 100 MPa). That covers the
//! one-root (supercritical, low and very high pressure) and the three-root (subcritical, between the
//! spinodals: water at 300 K / 1 atm) regions and the critical point. Then (Tc, Pc) and one part in
//! 10⁶ around it, a seeded sweep, every validation edge for T, P and each parameter (ω = 5 and −0.5
//! are accepted), and overflow (P up to 10²⁰ Pa and DBL_MAX, T = 10⁻³⁰⁰ K). Then (ADR 0006 item 6, M3c-3)
//! the hard fugacity rows: states a few kelvin above 0 where φ = exp(ln φ) is subnormal and φ·P is normal.

use std::hint::black_box as bb;

use crate::grids::Rng;
use crate::m_reaction::{edges, res};
use crate::t_gas::PRESETS;
use crate::tsv::{Class, Module, Tok, err, f};
use kimiya::element::GAS_CONSTANT;
use kimiya::gas::{
    PR_CO2, PR_METHANE, PR_NITROGEN, PR_WATER, PengRobinsonParams, fugacity_coefficient_pr, fugacity_pr, gas_density,
    ideal_gas_pressure, ideal_gas_volume, joule_thomson_ideal, partial_pressure,
    peng_robinson_pressure, peng_robinson_z, peng_robinson_z_roots, redlich_kwong_pressure,
    van_der_waals_pressure,
};

fn params(p: &PengRobinsonParams) -> Vec<Tok> {
    vec![f(p.tc), f(p.pc), f(p.omega)]
}

/// The four Peng–Robinson functions of (T, P) at one state.
fn pr_state(m: &mut Module, p: &PengRobinsonParams, t: f64, pr: f64) {
    let mut ins = params(p);
    ins.extend([f(t), f(pr)]);
    let z = res(peng_robinson_z(bb(p), bb(t), bb(pr)));
    m.row("gas::peng_robinson_z", Class::C, &ins, &[z]);
    let roots = match peng_robinson_z_roots(bb(p), bb(t), bb(pr)) {
        Ok(v) => {
            let mut out = vec![Tok::V(v.len())];
            out.extend(v.iter().map(|&x| f(x)));
            out
        }
        Err(e) => vec![err(&e)],
    };
    m.row("gas::peng_robinson_z_roots", Class::C, &ins, &roots);
    let phi = res(fugacity_coefficient_pr(bb(p), bb(t), bb(pr)));
    m.row("gas::fugacity_coefficient_pr", Class::C, &ins, &[phi]);
    let fu = res(fugacity_pr(bb(p), bb(t), bb(pr)));
    m.row("gas::fugacity_pr", Class::C, &ins, &[fu]);
}

fn with(base: &PengRobinsonParams, k: usize, x: f64) -> PengRobinsonParams {
    let mut v = [base.tc, base.pc, base.omega];
    v[k] = x;
    PengRobinsonParams {
        tc: v[0],
        pc: v[1],
        omega: v[2],
    }
}

const T_REDUCED: [f64; 20] = [
    0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1, 1.3, 1.5, 2.0,
    3.0, 5.0, 10.0,
];
const T_ABS: [f64; 6] = [200.0, 273.15, 280.0, 300.0, 304.1, 1000.0];
const P_REDUCED: [f64; 20] = [
    1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.01, 1.1, 1.5, 2.0, 5.0,
    10.0, 30.0, 60.0,
];
const P_ABS: [f64; 9] = [
    1000.0, 101_325.0, 1e6, 4e6, 7.375e6, 1e7, 90e6, 100e6, 200e6,
];

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("gas");
    let mut rng = Rng::new(0x3201);

    for (name, p) in &PRESETS {
        m.row(&format!("gas::{name}"), Class::A, &[], &params(p));
    }
    m.row(
        "gas::joule_thomson_ideal",
        Class::A,
        &[],
        &[f(joule_thomson_ideal())],
    );

    // ── ideal gas ──
    let igp = |m: &mut Module, n: f64, t: f64, v: f64| {
        let out = res(ideal_gas_pressure(bb(n), bb(t), bb(v)));
        m.row("gas::ideal_gas_pressure", Class::B, &[f(n), f(t), f(v)], &[out]);
    };
    let igv = |m: &mut Module, n: f64, t: f64, p: f64| {
        let out = res(ideal_gas_volume(bb(n), bb(t), bb(p)));
        m.row("gas::ideal_gas_volume", Class::B, &[f(n), f(t), f(p)], &[out]);
    };
    for x in edges() {
        igp(&mut m, x, 300.0, 0.025);
        igp(&mut m, 1.0, x, 0.025);
        igp(&mut m, 1.0, 300.0, x);
        igv(&mut m, x, 300.0, 101_325.0);
        igv(&mut m, 1.0, x, 101_325.0);
        igv(&mut m, 1.0, 300.0, x);
    }
    igp(&mut m, 2.0, 300.0, 0.05);
    igp(&mut m, 1.0, 300.0, 1.0);
    igv(&mut m, 1.0, 273.15, 101_325.0);
    for _ in 0..200 {
        let n = rng.signed_log_uniform(-3.0, 3.0);
        let t = rng.log_uniform(0.0, 4.0);
        igp(&mut m, n, t, rng.log_uniform(-6.0, 1.0));
        igv(&mut m, n, t, rng.log_uniform(0.0, 9.0));
    }

    // ── van der Waals ──
    let vdw = |m: &mut Module, n: f64, t: f64, v: f64, a: f64, b: f64| {
        let out = res(van_der_waals_pressure(bb(n), bb(t), bb(v), bb(a), bb(b)));
        m.row(
            "gas::van_der_waals_pressure",
            Class::B,
            &[f(n), f(t), f(v), f(a), f(b)],
            &[out],
        );
    };
    for x in edges() {
        vdw(&mut m, x, 300.0, 0.025, 0.1, 1e-5);
        vdw(&mut m, 1.0, x, 0.025, 0.1, 1e-5);
        vdw(&mut m, 1.0, 300.0, x, 0.1, 1e-5);
        vdw(&mut m, 1.0, 300.0, 0.025, x, 1e-5);
        vdw(&mut m, 1.0, 300.0, 0.025, 0.1, x);
    }
    for (n, t, v, a, b) in [
        (1.0, 300.0, 0.025, 0.0, 0.0),
        (1.0, 300.0, 1e-3, 0.3640, 4.267e-5),
        (1.0, 300.0, 0.001, 0.0, 0.01),
        (1.0, 300.0, 0.01, 0.0, 0.01),            // V = n·b exactly
        (1.0, 300.0, 0.01_f64.next_up(), 0.0, 0.01), // one ulp above
        (2.0, 300.0, 0.02, 0.1, 0.01),
        (-1.0, 300.0, 0.025, 0.1, 1e-5),          // negative n passes V > n·b
    ] {
        vdw(&mut m, n, t, v, a, b);
    }
    for _ in 0..200 {
        let n = rng.log_uniform(-3.0, 3.0);
        let t = rng.log_uniform(1.0, 3.5);
        let b = rng.log_uniform(-6.0, -4.0);
        let v = n * b * (1.0 + rng.log_uniform(-3.0, 4.0));
        vdw(&mut m, n, t, v, rng.uniform(0.0, 1.0), b);
    }

    // ── Dalton, density ──
    for (p, x) in [(101_325.0, 0.21), (0.0, 0.5), (-1.0, 0.5)] {
        let pp = partial_pressure(bb(p), bb(x));
        m.row("gas::partial_pressure", Class::B, &[f(p), f(x)], &[f(pp)]);
    }
    for x in edges() {
        for (p, y) in [(x, 0.21), (101_325.0, x)] {
            let pp = partial_pressure(bb(p), bb(y));
            m.row("gas::partial_pressure", Class::B, &[f(p), f(y)], &[f(pp)]);
        }
    }
    for _ in 0..100 {
        let p = rng.log_uniform(0.0, 8.0);
        let x = rng.unit();
        let pp = partial_pressure(bb(p), bb(x));
        m.row("gas::partial_pressure", Class::B, &[f(p), f(x)], &[f(pp)]);
    }
    let dens = |m: &mut Module, p: f64, mm: f64, t: f64| {
        let out = res(gas_density(bb(p), bb(mm), bb(t)));
        m.row("gas::gas_density", Class::B, &[f(p), f(mm), f(t)], &[out]);
    };
    for x in edges() {
        dens(&mut m, x, 0.029, 300.0);
        dens(&mut m, 101_325.0, x, 300.0);
        dens(&mut m, 101_325.0, 0.029, x);
    }
    dens(&mut m, 101_325.0, 0.029, 273.15);
    for _ in 0..150 {
        let p = rng.log_uniform(0.0, 8.0);
        let mm = rng.log_uniform(-3.0, 0.0);
        dens(&mut m, p, mm, rng.log_uniform(1.0, 3.5));
    }

    // ── Peng–Robinson pressure ──
    let prp = |m: &mut Module, p: &PengRobinsonParams, t: f64, v: f64| {
        let mut ins = params(p);
        ins.extend([f(t), f(v)]);
        let out = res(peng_robinson_pressure(bb(p), bb(t), bb(v)));
        m.row("gas::peng_robinson_pressure", Class::B, &ins, &[out]);
    };
    for (_, p) in &PRESETS {
        let b = 0.07780 * GAS_CONSTANT * p.tc / p.pc; // the co-volume, as pr_coefficients forms it
        for tr in [0.5, 0.9, 1.0, 1.5, 3.0] {
            let t = tr * p.tc;
            for k in [1.0, 1.001, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0, 1e3, 1e4, 1e6] {
                prp(&mut m, p, t, b * k);
            }
            prp(&mut m, p, t, b.next_up());
        }
        for v in [1e-4, 1e-3, 0.025, 1.0] {
            prp(&mut m, p, 300.0, v);
        }
    }
    for x in edges() {
        prp(&mut m, &PR_NITROGEN, x, 0.025);
        prp(&mut m, &PR_NITROGEN, 300.0, x);
        for k in 0..3 {
            prp(&mut m, &with(&PR_NITROGEN, k, x), 300.0, 0.025);
        }
    }
    for _ in 0..200 {
        let p = &PRESETS[usize::try_from(rng.int(0, 4)).unwrap_or(0)].1;
        let t = p.tc * rng.log_uniform(-0.6, 1.0);
        let b = 0.07780 * GAS_CONSTANT * p.tc / p.pc;
        prp(&mut m, p, t, b * (1.0 + rng.log_uniform(-3.0, 5.0)));
    }

    // ── Peng–Robinson Z, roots, fugacity ──
    for (_, p) in &PRESETS {
        let temps: Vec<f64> = T_REDUCED
            .iter()
            .map(|tr| tr * p.tc)
            .chain(T_ABS)
            .collect();
        let pressures: Vec<f64> = P_REDUCED
            .iter()
            .map(|pr| pr * p.pc)
            .chain(P_ABS)
            .collect();
        for &t in &temps {
            for &pr in &pressures {
                pr_state(&mut m, p, t, pr);
            }
        }
        for dt in [-1e-6, 0.0, 1e-6] {
            for dp in [-1e-6, 0.0, 1e-6] {
                pr_state(&mut m, p, p.tc * (1.0 + dt), p.pc * (1.0 + dp));
            }
        }
    }
    for _ in 0..400 {
        let p = &PRESETS[usize::try_from(rng.int(0, 4)).unwrap_or(0)].1;
        let t = p.tc * rng.log_uniform(-0.6, 1.3);
        let pr = rng.log_uniform(0.0, 9.0);
        pr_state(&mut m, p, t, pr);
    }
    for x in edges() {
        pr_state(&mut m, &PR_NITROGEN, x, 101_325.0);
        pr_state(&mut m, &PR_NITROGEN, 300.0, x);
        for k in 0..3 {
            pr_state(&mut m, &with(&PR_NITROGEN, k, x), 300.0, 1e5);
        }
    }
    for (t, pr) in [
        (300.0, 1e10),
        (300.0, 1e12),
        (300.0, 1e15),
        (300.0, 1e20),
        (1e-300, 101_325.0),
        (1e-10, 101_325.0),
    ] {
        pr_state(&mut m, &PR_NITROGEN, t, pr);
    }
    for omega in [5.0, -0.5, 0.49, 1.0] {
        pr_state(&mut m, &with(&PR_NITROGEN, 2, omega), 100.0, 1e6);
        pr_state(&mut m, &with(&PR_NITROGEN, 2, omega), 300.0, 1e6);
    }
    // ADR 0006's hard rows (M3c-3): liquid-like roots at a few kelvin, where ln φ is between −745.2 and −708, so
    // φ = exp(ln φ) is subnormal (or rounds to 0: Err) and fugacity_pr's φ·P scales it back into the normal
    // range: water near 10 K (at 10.25 K and 5·10⁵ Pa, φ ≈ 2.4·10⁻³⁰⁹), CO₂ near 4 K, methane near 1.7 K.
    // Hand-picked states, then a seeded sweep over the same boxes (about half of its states land in the range).
    for t in [9.8, 9.85, 9.9, 9.95, 10.0, 10.05, 10.1, 10.15, 10.2, 10.25] {
        for pr in [1e4, 5e5, 1e7] {
            pr_state(&mut m, &PR_WATER, t, pr);
        }
    }
    pr_state(&mut m, &PR_WATER, 9.76, 1e3);
    pr_state(&mut m, &PR_WATER, 9.76, 1e4); // ln φ ≈ −745.2: φ rounds to 0, Err
    for t in [3.9, 3.95, 4.0, 4.05] {
        for pr in [1e4, 1e6] {
            pr_state(&mut m, &PR_CO2, t, pr);
        }
    }
    for t in [1.7, 1.75] {
        for pr in [1e4, 1e7] {
            pr_state(&mut m, &PR_METHANE, t, pr);
        }
    }
    let mut hard = Rng::new(0x3203);
    for (p, lo, hi) in [(&PR_WATER, 9.75, 10.3), (&PR_CO2, 3.85, 4.1), (&PR_METHANE, 1.68, 1.82)] {
        for _ in 0..100 {
            let t = hard.uniform(lo, hi);
            let pr = hard.log_uniform(3.0, 7.5);
            pr_state(&mut m, p, t, pr);
        }
    }

    // ── Redlich–Kwong ──
    let rk = |m: &mut Module, tc: f64, pc: f64, t: f64, v: f64| {
        let out = res(redlich_kwong_pressure(bb(tc), bb(pc), bb(t), bb(v)));
        m.row(
            "gas::redlich_kwong_pressure",
            Class::C,
            &[f(tc), f(pc), f(t), f(v)],
            &[out],
        );
    };
    for (_, p) in &PRESETS {
        let b = 0.08664 * GAS_CONSTANT * p.tc / p.pc; // the co-volume, as redlich_kwong_pressure forms it
        for tr in [0.5, 0.9, 1.0, 1.5, 3.0] {
            let t = tr * p.tc;
            for k in [1.0, 1.001, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0, 1e3, 1e4, 1e6] {
                rk(&mut m, p.tc, p.pc, t, b * k);
            }
            rk(&mut m, p.tc, p.pc, t, b.next_up());
        }
    }
    for x in edges() {
        rk(&mut m, 126.2, 3_394_000.0, x, 0.025);
        rk(&mut m, x, 3_394_000.0, 300.0, 0.025);
        rk(&mut m, 126.2, x, 300.0, 0.025);
        rk(&mut m, 126.2, 3_394_000.0, 300.0, x);
    }
    rk(&mut m, 126.2, 3_394_000.0, 300.0, 1.0);
    for _ in 0..200 {
        let p = &PRESETS[usize::try_from(rng.int(0, 4)).unwrap_or(0)].1;
        let t = p.tc * rng.log_uniform(-0.6, 1.0);
        let b = 0.08664 * GAS_CONSTANT * p.tc / p.pc;
        rk(&mut m, p.tc, p.pc, t, b * (1.0 + rng.log_uniform(-3.0, 5.0)));
    }
    m
}
