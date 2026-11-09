//! `phase.tsv` — Raoult, Henry, colligative properties, Clausius–Clapeyron, Antoine and the phase
//! rule (port plan §4, M2), and the `ANTOINE_DATA` table.
//!
//! - Class B (+ − × ÷ only, expected 0 ulp): `raoult_pressure`, `raoult_total_binary`,
//!   `henry_pressure`, `henry_solubility`, `boiling_point_elevation`, `freezing_point_depression`,
//!   `osmotic_pressure`. Each argument over the validation edges (the mole fraction also at 0, 1
//!   and one ulp outside), the others typical; the Rust tests' values; a seeded sweep.
//! - Class C: `clausius_clapeyron` (one `exp`) over the edges, the Rust tests' values, exp's
//!   overflow and underflow, and a seeded sweep; `antoine_pressure_mmhg` / `antoine_pressure_pa`
//!   (one `pow`), see below.
//! - Class A: `lookup_antoine` (every name, case and spacing variants, unknown names; the coefficient
//!   set flattened as `f:a f:b f:c f:t_min_c f:t_max_c`, or `n:`), `phase_rule` (every pair of 13
//!   u32 probes from 0 to u32::MAX and a seeded sweep; `i:F` or `E:InvalidInput`), and the table:
//!   `phase::ANTOINE_DATA` (`v:8`) and one `phase::ANTOINE_DATA[]` row per entry (input `i:k`,
//!   output `s:name` and the five fields).
//!
//! **Antoine** (ADR 0002 #4: the Cyrius functions return `Err` outside a set's `t_min_c..t_max_c`,
//! at or below the pole and for non-finite input; Rust 1.2 extrapolates or returns NaN). The input
//! is the coefficient set flattened, then T. For each built-in set: both range ends, one ulp inside
//! and outside each, 0.5, 1 and 10 °C outside, a 5 °C ladder and seeded points across the range,
//! the pole −C and one ulp either side, −C + 1, −C − 10, −300, ±0, NaN and ±∞, and the Rust tests'
//! temperatures. Then custom sets: each field NaN or infinite, an inverted range, a range that
//! contains the pole, exponents that overflow and underflow, and the unit set A = B = 0, C = 1 on
//! [0, 0] at T = 0, whose rows are class A: 10⁰ = 1 exactly, so `antoine_pressure_pa` returns the
//! private `PA_PER_MMHG` itself (`t_phase.rs` generates `_KM_PA_PER_MMHG` from that call).

use std::hint::black_box as bb;

use crate::grids::Rng;
use crate::m_reaction::{edges, res};
use crate::t_phase::fields;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use kimiya::phase::{
    ANTOINE_DATA, AntoineCoeffs, KB_BENZENE, KB_CHLOROFORM, KB_ETHANOL, KB_WATER, KF_BENZENE,
    KF_CAMPHOR, KF_CYCLOHEXANE, KF_WATER, antoine_pressure_mmhg, antoine_pressure_pa,
    boiling_point_elevation, clausius_clapeyron, freezing_point_depression, henry_pressure,
    henry_solubility, lookup_antoine, osmotic_pressure, phase_rule, raoult_pressure,
    raoult_total_binary,
};

/// The mole-fraction probes: the edges, both ends of [0, 1] and one ulp outside each.
fn fractions() -> Vec<f64> {
    let mut v = edges();
    v.extend([1.0_f64.next_up(), 1.0_f64.next_down(), 0.21, 0.8, 0.6]);
    v
}

fn coeff_toks(c: &AntoineCoeffs) -> Vec<Tok> {
    fields(c).iter().map(|&x| f(x)).collect()
}

fn antoine(m: &mut Module, c: &AntoineCoeffs, t: f64, class: Class) {
    let mut ins = coeff_toks(c);
    ins.push(f(t));
    let mm = antoine_pressure_mmhg(bb(c), bb(t));
    m.row("phase::antoine_pressure_mmhg", class, &ins, &[f(mm)]);
    let pa = antoine_pressure_pa(bb(c), bb(t));
    m.row("phase::antoine_pressure_pa", class, &ins, &[f(pa)]);
}

fn antoine_temperatures(c: &AntoineCoeffs, rng: &mut Rng) -> Vec<f64> {
    let (lo, hi) = (c.t_min_c, c.t_max_c);
    let mut t = vec![
        lo,
        hi,
        lo.next_up(),
        lo.next_down(),
        hi.next_up(),
        hi.next_down(),
        lo - 0.5,
        lo - 1.0,
        lo - 10.0,
        hi + 0.5,
        hi + 1.0,
        hi + 10.0,
        -c.c,
        (-c.c).next_up(),
        (-c.c).next_down(),
        -c.c + 1.0,
        -c.c - 10.0,
        -300.0,
        0.0,
        -0.0,
        25.0,
        50.0,
        64.7,
        78.4,
        100.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    let mut x = lo;
    while x < hi {
        t.push(x);
        x += 5.0;
    }
    for _ in 0..12 {
        t.push(rng.uniform(lo, hi));
    }
    t
}

pub fn emit() -> Module {
    let mut m = Module::new("phase");
    let mut rng = Rng::new(0x3101);

    // ── Raoult ──
    let raoult = |m: &mut Module, x: f64, p: f64| {
        let out = res(raoult_pressure(bb(x), bb(p)));
        m.row("phase::raoult_pressure", Class::B, &[f(x), f(p)], &[out]);
    };
    for x in fractions() {
        raoult(&mut m, x, 101_325.0);
    }
    for p in edges() {
        raoult(&mut m, 0.5, p);
    }
    for (x, p) in [(0.8, 101_325.0), (1.0, 101_325.0), (1.5, 101_325.0), (-0.1, 101_325.0)] {
        raoult(&mut m, x, p);
    }
    for _ in 0..150 {
        let x = rng.unit();
        let p = rng.log_uniform(0.0, 7.0);
        raoult(&mut m, x, p);
    }
    let binary = |m: &mut Module, x: f64, pa: f64, pb: f64| {
        let out = res(raoult_total_binary(bb(x), bb(pa), bb(pb)));
        m.row("phase::raoult_total_binary", Class::B, &[f(x), f(pa), f(pb)], &[out]);
    };
    for x in fractions() {
        binary(&mut m, x, 100.0, 50.0);
    }
    for p in edges() {
        binary(&mut m, 0.5, p, 50.0);
        binary(&mut m, 0.5, 100.0, p);
    }
    binary(&mut m, 0.6, 100.0, 50.0);
    for _ in 0..150 {
        let x = rng.unit();
        let pa = rng.log_uniform(0.0, 7.0);
        let pb = rng.log_uniform(0.0, 7.0);
        binary(&mut m, x, pa, pb);
    }

    // ── Henry ──
    let hp = |m: &mut Module, kh: f64, x: f64| {
        let out = res(henry_pressure(bb(kh), bb(x)));
        m.row("phase::henry_pressure", Class::B, &[f(kh), f(x)], &[out]);
    };
    let hs = |m: &mut Module, p: f64, kh: f64| {
        let out = res(henry_solubility(bb(p), bb(kh)));
        m.row("phase::henry_solubility", Class::B, &[f(p), f(kh)], &[out]);
    };
    for x in fractions() {
        hp(&mut m, 1e8, x);
    }
    for k in edges() {
        hp(&mut m, k, 0.001);
        hs(&mut m, 21_278.0, k);
        hs(&mut m, k, 4.3e9);
    }
    for (p, kh) in [(21_278.0, 4.3e9), (0.0, 4.3e9), (-1.0, 4.3e9), (1e8 * 0.001, 1e8)] {
        hs(&mut m, p, kh);
    }
    for _ in 0..150 {
        let kh = rng.log_uniform(3.0, 11.0);
        let x = rng.unit();
        hp(&mut m, kh, x);
        let p = rng.log_uniform(-2.0, 7.0);
        hs(&mut m, p, kh);
    }

    // ── colligative ──
    let consts = [
        KB_WATER,
        KB_BENZENE,
        KB_ETHANOL,
        KB_CHLOROFORM,
        KF_WATER,
        KF_BENZENE,
        KF_CYCLOHEXANE,
        KF_CAMPHOR,
    ];
    let coll = |m: &mut Module, k: f64, mol: f64, i: f64| {
        let bp = boiling_point_elevation(bb(k), bb(mol), bb(i));
        m.row("phase::boiling_point_elevation", Class::B, &[f(k), f(mol), f(i)], &[f(bp)]);
        let fp = freezing_point_depression(bb(k), bb(mol), bb(i));
        m.row("phase::freezing_point_depression", Class::B, &[f(k), f(mol), f(i)], &[f(fp)]);
    };
    for k in consts {
        for (mol, i) in [(1.0, 2.0), (1.0, 1.0), (0.1, 3.0), (2.5, 1.9)] {
            coll(&mut m, k, mol, i);
        }
    }
    for x in edges() {
        coll(&mut m, x, 1.0, 2.0);
        coll(&mut m, KB_WATER, x, 2.0);
        coll(&mut m, KF_WATER, 1.0, x);
    }
    for _ in 0..100 {
        let k = consts[usize::try_from(rng.int(0, 7)).unwrap_or(0)];
        let mol = rng.log_uniform(-4.0, 1.0);
        let i = rng.uniform(1.0, 4.0);
        coll(&mut m, k, mol, i);
    }

    // ── osmotic pressure ──
    let osm = |m: &mut Module, i: f64, c: f64, t: f64| {
        let out = res(osmotic_pressure(bb(i), bb(c), bb(t)));
        m.row("phase::osmotic_pressure", Class::B, &[f(i), f(c), f(t)], &[out]);
    };
    for x in edges() {
        osm(&mut m, x, 0.1, 298.15);
        osm(&mut m, 2.0, x, 298.15);
        osm(&mut m, 2.0, 0.1, x);
    }
    for (i, c, t) in [(2.0, 0.1, 298.15), (1.0, 0.0, 298.15), (1e300, 1e300, 1e300)] {
        osm(&mut m, i, c, t);
    }
    for _ in 0..150 {
        let i = rng.uniform(1.0, 4.0);
        let c = rng.log_uniform(-6.0, 1.0);
        let t = rng.log_uniform(0.0, 3.5);
        osm(&mut m, i, c, t);
    }

    // ── Clausius–Clapeyron ──
    let cc = |m: &mut Module, p1: f64, t1: f64, t2: f64, dh: f64| {
        let out = res(clausius_clapeyron(bb(p1), bb(t1), bb(t2), bb(dh)));
        m.row("phase::clausius_clapeyron", Class::C, &[f(p1), f(t1), f(t2), f(dh)], &[out]);
    };
    for x in edges() {
        cc(&mut m, x, 373.15, 363.15, 40_670.0);
        cc(&mut m, 101_325.0, x, 363.15, 40_670.0);
        cc(&mut m, 101_325.0, 373.15, x, 40_670.0);
        cc(&mut m, 101_325.0, 373.15, 363.15, x);
    }
    for (p1, t1, t2, dh) in [
        (101_325.0, 373.15, 363.15, 40_670.0),
        (101_325.0, 373.15, 373.15, 40_670.0),
        (101_325.0, 0.0, 373.15, 40_670.0),
        (101_325.0, 1000.0, 1.0, -1e7),    // exponent ≫ 709: overflows to +∞
        (101_325.0, 1.0, 1000.0, -1e7),    // exponent ≪ −745: underflows to 0
        (1e300, 300.0, 400.0, 1e5),
        (5e-324, 300.0, 299.0, 1e5),
    ] {
        cc(&mut m, p1, t1, t2, dh);
    }
    for _ in 0..400 {
        let p1 = rng.log_uniform(0.0, 7.0);
        let t1 = rng.uniform(150.0, 700.0);
        let t2 = rng.uniform(150.0, 700.0);
        let dh = rng.signed_log_uniform(3.0, 5.3);
        cc(&mut m, p1, t1, t2, dh);
    }

    // ── Antoine ──
    for (_, c) in ANTOINE_DATA {
        for t in antoine_temperatures(c, &mut rng) {
            antoine(&mut m, c, t, Class::C);
        }
    }
    let water = &ANTOINE_DATA[0].1;
    let with = |k: usize, x: f64| {
        let mut v = fields(water);
        v[k] = x;
        AntoineCoeffs {
            a: v[0],
            b: v[1],
            c: v[2],
            t_min_c: v[3],
            t_max_c: v[4],
        }
    };
    for k in 0..5 {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for t in [1.0, 50.0, 100.0] {
                antoine(&mut m, &with(k, x), t, Class::C);
            }
        }
    }
    let custom = [
        // inverted range
        AntoineCoeffs {
            t_min_c: 100.0,
            t_max_c: 1.0,
            ..water.clone()
        },
        // the pole (T = 10) inside the range [0, 100]
        AntoineCoeffs {
            c: -10.0,
            t_min_c: 0.0,
            t_max_c: 100.0,
            ..water.clone()
        },
        // 10^(A − B/(C + T)) overflows, and underflows
        AntoineCoeffs {
            a: 400.0,
            ..water.clone()
        },
        AntoineCoeffs {
            a: -400.0,
            ..water.clone()
        },
        // a degenerate one-point range
        AntoineCoeffs {
            t_min_c: 50.0,
            t_max_c: 50.0,
            ..water.clone()
        },
    ];
    for c in &custom {
        for t in [0.0, 1.0, 5.0, 10.0, 10.5, 50.0, 100.0] {
            antoine(&mut m, c, t, Class::C);
        }
    }
    let unit = AntoineCoeffs {
        a: 0.0,
        b: 0.0,
        c: 1.0,
        t_min_c: 0.0,
        t_max_c: 0.0,
    };
    antoine(&mut m, &unit, 0.0, Class::A);

    // ── lookup_antoine and the table ──
    let mut names: Vec<String> = Vec::new();
    for (n, _) in ANTOINE_DATA {
        names.push((*n).to_string());
        names.push(n.to_uppercase());
        names.push(format!(" {n}"));
        names.push(format!("{n} "));
    }
    names.extend(
        ["", "unobtainium", "Water", "diethyl ether", "wate", "waterr"]
            .iter()
            .map(|x| (*x).to_string()),
    );
    for n in &names {
        let out = lookup_antoine(bb(n)).map_or_else(|| vec![Tok::None], coeff_toks);
        m.row("phase::lookup_antoine", Class::A, &[s(n.clone())], &out);
    }
    m.row(
        "phase::ANTOINE_DATA",
        Class::A,
        &[],
        &[Tok::V(ANTOINE_DATA.len())],
    );
    for (k, (n, c)) in ANTOINE_DATA.iter().enumerate() {
        let mut out = vec![s(*n)];
        out.extend(coeff_toks(c));
        m.row("phase::ANTOINE_DATA[]", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &out);
    }

    // ── phase rule ──
    let probes: [u32; 13] = [
        0,
        1,
        2,
        3,
        4,
        10,
        100,
        0x7fff_ffff,
        0x8000_0000,
        3_000_000_000,
        u32::MAX - 2,
        u32::MAX - 1,
        u32::MAX,
    ];
    let rule = |m: &mut Module, c: u32, p: u32| {
        let out = phase_rule(bb(c), bb(p)).map_or_else(|e| err(&e), i);
        m.row("phase::phase_rule", Class::A, &[i(c), i(p)], &[out]);
    };
    for &c in &probes {
        for &p in &probes {
            rule(&mut m, c, p);
        }
    }
    for _ in 0..100 {
        let c = u32::try_from(rng.int(0, 20)).unwrap_or(0);
        let p = u32::try_from(rng.int(0, 25)).unwrap_or(0);
        rule(&mut m, c, p);
    }
    m
}
