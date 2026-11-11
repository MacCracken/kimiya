//! `electrochemistry.tsv` — the standard potentials, the Nernst equation, cell potentials, cell Gibbs energy and
//! Faraday's laws (port plan §4, M3a).
//!
//! The electron count (Rust's `u8`) is an `i:` token.
//! - Class A: `electrochemistry::STANDARD_POTENTIALS` (`v:49`) and one `STANDARD_POTENTIALS[]` row per entry (input
//!   `i:k`; output `s:equation s:couple i:electrons f:standard_potential`), `lookup_half_reaction` (every couple,
//!   its lower- and upper-case forms, a leading / trailing space, every proper prefix, the equation text, `""`,
//!   `"Xx/Xx"`; the four fields or `n:`), `is_spontaneous_cell` (the specials, ± tiny, ± the table's potentials).
//! - Class B (+ − × ÷ only, expected 0 ulp): `standard_cell_potential` and `cell_potential_from_couples` (every
//!   ordered pair of the 49 couples; unknown couples on either side, the cathode first; the edges),
//!   `cell_gibbs_energy`, `faraday_mass_deposited`, `faraday_charge_required`, `faraday_mass_ratio` (every
//!   argument over the edges, the electron counts 0, 1 … 6, 127, 128, 254, 255, charges and masses over every
//!   decade from 10⁻³⁰⁰ to 10³⁰⁰, the Rust tests' values, seeded sweeps), `charge_to_moles_electrons`,
//!   `moles_electrons_to_charge` (the edges, every decade from the smallest subnormal to DBL_MAX, seeded).
//! - Class C, one `ln`: `nernst_potential` and `nernst_potential_25c` — the edges per argument, the electron
//!   counts, T from 10⁻³⁰⁰ to 10³⁰⁰ K, Q over every tenth decade from 10⁻³²⁰ to 10³⁰⁰, every couple's E° and n at
//!   seven quotients (the standard rows), the Rust tests' values, a seeded sweep, and **near-cancellation rows on
//!   purpose**: E° = (RT/nF)·ln Q to within a few ulp, so E ≈ 0 of both signs (E° placed with `ln_arith`, + − × ÷
//!   only, as grids may not call a transcendental).
//!
//! - ADR 0002 #29 (the overflow audit's B2): the three Faraday functions' rows are noted where the port's correctly
//!   rounded value differs from Rust's, and `nernst_potential`'s where Rust's sequence with an unbounded exponent (the
//!   final subtraction rounded once) does (`audit.rs`; in-range rows keep Rust's bits); `audit::emit_electrochemistry`
//!   appends their edge rows (R·T past DBL_MAX, a subnormal (R·T)/(n·F), partial products past either end, an
//!   equivalent weight underflowing, midpoints) with their own seed.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::audit;
use crate::grids::{Rng, decades};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, b, f, i, s};
use kimiya::electrochemistry::{
    FARADAY, HalfReaction, STANDARD_POTENTIALS, cell_gibbs_energy, cell_potential_from_couples,
    charge_to_moles_electrons, is_spontaneous_cell, lookup_half_reaction, moles_electrons_to_charge,
    nernst_potential_25c, standard_cell_potential,
};
use kimiya::element::GAS_CONSTANT;

/// A half-reaction's four fields, in declaration order.
pub fn hr_toks(h: &HalfReaction) -> Vec<Tok> {
    vec![s(h.equation), s(h.couple), i(h.electrons), f(h.standard_potential)]
}

/// The electron counts every n-taking function is probed with (Rust's u8: 0 is the Err case).
const COUNTS: [u8; 11] = [0, 1, 2, 3, 4, 5, 6, 127, 128, 254, 255];

/// ln(x) for a normal x > 0 from + − × ÷ and bit manipulation only (grids may not call a transcendental std
/// method, `grids.rs`): x = 2ᵉ·m, m in [√½, √2), ln m = 2·atanh((m − 1)/(m + 1)) by its series. Within a few ulp:
/// it only places the near-cancellation E° values, whose rows are near cancellation either way.
fn ln_arith(x: f64) -> f64 {
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
    2.0 * sum + e * std::f64::consts::LN_2
}

/// `x` moved `d` ulps (up for d > 0).
fn ulps(x: f64, d: i32) -> f64 {
    let mut y = x;
    for _ in 0..d.unsigned_abs() {
        y = if d > 0 { y.next_up() } else { y.next_down() };
    }
    y
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("electrochemistry");
    let mut rng = Rng::new(0x3b01);

    // ── the table and the lookup ──
    m.row(
        "electrochemistry::STANDARD_POTENTIALS",
        Class::A,
        &[],
        &[Tok::V(STANDARD_POTENTIALS.len())],
    );
    for (k, h) in STANDARD_POTENTIALS.iter().enumerate() {
        m.row(
            "electrochemistry::STANDARD_POTENTIALS[]",
            Class::A,
            &[i(i64::try_from(k).unwrap_or(-1))],
            &hr_toks(h),
        );
    }
    let mut names: Vec<String> = Vec::new();
    for h in STANDARD_POTENTIALS {
        let c = h.couple;
        names.push(c.to_string());
        names.push(c.to_lowercase());
        names.push(c.to_uppercase());
        names.push(format!(" {c}"));
        names.push(format!("{c} "));
        for n in 1..c.len() {
            names.push(c[..n].to_string());
        }
        names.push(h.equation.to_string());
    }
    names.extend(["", "Xx/Xx", "Cu2+/Cu2+", "Cu/Cu2+", "H+/H2 ", "cu2+/cu", "Cu2+ /Cu"].iter().map(|x| (*x).to_string()));
    names.sort();
    names.dedup();
    for n in &names {
        let out = lookup_half_reaction(bb(n)).map_or_else(|| vec![Tok::None], hr_toks);
        m.row("electrochemistry::lookup_half_reaction", Class::A, &[s(n.clone())], &out);
    }

    // ── the Nernst equation ──
    let nernst = audit::nernst;
    let n25 = |m: &mut Module, e0: f64, n: u8, q: f64| {
        let out = res(nernst_potential_25c(bb(e0), bb(n), bb(q)));
        m.row("electrochemistry::nernst_potential_25c", Class::C, &[f(e0), i(n), f(q)], &[out]);
    };
    for x in edges() {
        nernst(&mut m, x, 2, 298.15, 0.1);
        nernst(&mut m, 0.342, 2, x, 0.1);
        nernst(&mut m, 0.342, 2, 298.15, x);
        n25(&mut m, x, 2, 0.1);
        n25(&mut m, 0.342, 2, x);
    }
    for n in COUNTS {
        nernst(&mut m, 0.342, n, 298.15, 0.1);
        nernst(&mut m, 0.342, n, 0.0, 0.0); // the temperature is checked first
        nernst(&mut m, f64::NAN, n, 298.15, 0.0); // then n, then Q, then E°
        n25(&mut m, 0.342, n, 100.0);
        n25(&mut m, f64::NAN, n, -1.0);
    }
    for t in decades(-300, 300, &["1"]).into_iter().step_by(25).chain([1.0, 77.0, 273.15, 298.15, 373.15, 1000.0, 5000.0]) {
        nernst(&mut m, 0.342, 2, t, 0.1);
        nernst(&mut m, -0.762, 2, t, 1e10);
    }
    for q in decades(-320, 300, &["1", "3"]).into_iter().step_by(10).chain([1.0, 1.0_f64.next_up(), 1.0_f64.next_down()]) {
        nernst(&mut m, 0.342, 2, 298.15, q);
        n25(&mut m, 1.229, 4, q);
    }
    // the standard rows: every couple's E° and n, at seven quotients
    for h in STANDARD_POTENTIALS {
        for q in [1e-10, 0.01, 0.1, 1.0, 10.0, 100.0, 1e10] {
            nernst(&mut m, h.standard_potential, h.electrons, 298.15, q);
            n25(&mut m, h.standard_potential, h.electrons, q);
        }
        nernst(&mut m, h.standard_potential, h.electrons, 310.15, 0.5);
    }
    // the Rust tests' values
    for (e0, n, t, q) in [
        (0.342, 2, 298.15, 1.0),
        (0.342, 2, 298.15, 100.0),
        (0.342, 2, 298.15, 0.01),
        (0.342, 2, 298.15, 0.1),
        (0.342, 0, 298.15, 1.0),
        (0.342, 2, 298.15, 0.0),
        (0.342, 2, 0.0, 1.0),
        (0.342, 2, 298.15, -1.0),
        (0.342, 2, f64::NAN, 1.0),
        (0.342, 2, f64::INFINITY, 1.0),
        (0.342, 2, 298.15, f64::NAN),
        (0.342, 2, 298.15, f64::INFINITY),
        (f64::NAN, 2, 298.15, 1.0),
        (1e308, 1, 1e308, 1e300), // RT overflows: E° − ∞
    ] {
        nernst(&mut m, e0, n, t, q);
    }
    for (e0, n, q) in [(0.342, 2, 0.1), (0.342, 2, 0.01), (0.342, 2, f64::NAN)] {
        n25(&mut m, e0, n, q);
    }
    for _ in 0..400 {
        let e0 = rng.uniform(-3.1, 2.9);
        let n = u8::try_from(rng.int(1, 6)).unwrap_or(1);
        let t = rng.uniform(200.0, 1500.0);
        let q = rng.log_uniform(-30.0, 30.0);
        nernst(&mut m, e0, n, t, q);
        n25(&mut m, e0, n, q);
    }
    // Near cancellation, on purpose: E° = (RT/nF)·ln Q to within a few ulp, so E is a few ulp of the log term
    // either side of 0 — one ulp of ln Q is many of E (the reason for _km_ln, ADR 0005 amendment M3a).
    let near = |m: &mut Module, rng: &mut Rng, t: f64, at25: bool| {
        let n = u8::try_from(rng.int(1, 5)).unwrap_or(1);
        let q = rng.log_uniform(-60.0, 60.0);
        if !(q > 0.0 && q.is_finite()) || q == 1.0 {
            return;
        }
        let factor = GAS_CONSTANT * t / (f64::from(n) * FARADAY);
        let e00 = factor * ln_arith(q);
        for d in -2..=2 {
            let e0 = ulps(e00, d);
            if at25 {
                let out = res(nernst_potential_25c(bb(e0), bb(n), bb(q)));
                m.row("electrochemistry::nernst_potential_25c", Class::C, &[f(e0), i(n), f(q)], &[out]);
            } else {
                audit::nernst(m, e0, n, t, q);
            }
        }
    };
    for _ in 0..40 {
        let t = rng.uniform(250.0, 1000.0);
        near(&mut m, &mut rng, t, false);
        near(&mut m, &mut rng, 298.15, false);
        near(&mut m, &mut rng, 298.15, true);
    }

    // ── cell potentials ──
    for c in STANDARD_POTENTIALS {
        for a in STANDARD_POTENTIALS {
            let e = standard_cell_potential(bb(c.standard_potential), bb(a.standard_potential));
            m.row(
                "electrochemistry::standard_cell_potential",
                Class::B,
                &[f(c.standard_potential), f(a.standard_potential)],
                &[f(e)],
            );
            let out = res(cell_potential_from_couples(bb(c.couple), bb(a.couple)));
            m.row("electrochemistry::cell_potential_from_couples", Class::B, &[s(c.couple), s(a.couple)], &[out]);
        }
    }
    for (c, a) in [
        ("Xx/Xx", "Zn2+/Zn"),
        ("Cu2+/Cu", "Xx/Xx"),
        ("Xx/Xx", "Yy/Yy"),
        ("", "Zn2+/Zn"),
        ("Cu2+/Cu", ""),
        ("cu2+/cu", "Zn2+/Zn"),
        ("Cu2+/Cu", "Zn2+/Zn"),
    ] {
        let out = res(cell_potential_from_couples(bb(c), bb(a)));
        m.row("electrochemistry::cell_potential_from_couples", Class::B, &[s(c), s(a)], &[out]);
    }
    for x in edges() {
        for y in [0.342, -0.762, x] {
            let e = standard_cell_potential(bb(x), bb(y));
            m.row("electrochemistry::standard_cell_potential", Class::B, &[f(x), f(y)], &[f(e)]);
            let e = standard_cell_potential(bb(y), bb(x));
            m.row("electrochemistry::standard_cell_potential", Class::B, &[f(y), f(x)], &[f(e)]);
        }
    }
    let mut probes = edges();
    probes.extend([1.1, -0.5, 1e-300, -1e-300, 5e-324, -5e-324, -0.0]);
    probes.extend(STANDARD_POTENTIALS.iter().map(|h| h.standard_potential));
    for x in probes {
        m.row("electrochemistry::is_spontaneous_cell", Class::A, &[f(x)], &[b(is_spontaneous_cell(bb(x)))]);
    }

    // ── Gibbs energy ──
    let gibbs = |m: &mut Module, n: u8, e: f64| {
        let out = res(cell_gibbs_energy(bb(n), bb(e)));
        m.row("electrochemistry::cell_gibbs_energy", Class::B, &[i(n), f(e)], &[out]);
    };
    for n in COUNTS {
        for e in [1.104, 1.1, -0.5, 0.0, -0.0, 1e300, f64::MAX, -f64::MAX, 5e-324, f64::NAN, f64::INFINITY] {
            gibbs(&mut m, n, e);
        }
    }
    for x in edges() {
        gibbs(&mut m, 2, x);
    }
    for _ in 0..200 {
        let n = u8::try_from(rng.int(1, 255)).unwrap_or(1);
        let e = rng.signed_log_uniform(-6.0, 1.0);
        gibbs(&mut m, n, e);
    }

    // ── Faraday's laws ──
    let dep = audit::fmass;
    let req = audit::fcharge;
    let ratio = audit::fratio;
    for x in edges() {
        dep(&mut m, x, 63.546, 2);
        dep(&mut m, FARADAY, x, 2);
        req(&mut m, x, 63.546, 2);
        req(&mut m, 10.0, x, 2);
        ratio(&mut m, x, 63.546, 2, 107.868, 1);
        ratio(&mut m, 31.773, x, 2, 107.868, 1);
        ratio(&mut m, 31.773, 63.546, 2, x, 1);
    }
    for n in COUNTS {
        dep(&mut m, FARADAY, 63.546, n);
        req(&mut m, 10.0, 63.546, n);
        ratio(&mut m, 31.773, 63.546, n, 107.868, 1);
        ratio(&mut m, 31.773, 63.546, 2, 107.868, n);
        dep(&mut m, f64::NAN, f64::NAN, n); // n is checked first, then M, then the charge
        req(&mut m, f64::NAN, -1.0, n);
        ratio(&mut m, f64::NAN, 0.0, n, 0.0, n);
    }
    // wide magnitude grids: charge and mass over every decade, molar masses from light to absurd
    for q in decades(-300, 300, &["1", "3"]).into_iter().step_by(4) {
        for (mm, n) in [(63.546, 2), (1.008, 1), (107.868, 1), (238.03, 3), (1e-300, 255), (1e300, 1)] {
            dep(&mut m, q, mm, n);
            dep(&mut m, -q, mm, n);
            req(&mut m, q, mm, n);
            ratio(&mut m, q, mm, n, 26.98, 3);
        }
    }
    for (q, mm, n) in [(FARADAY, 63.546, 2), (-FARADAY, 63.546, 2), (1000.0, 63.546, 0), (FARADAY, 0.0, 2), (FARADAY, -63.546, 2)] {
        dep(&mut m, q, mm, n);
    }
    for (ma, mm, n) in [(10.0, 63.546, 2), (10.0, 63.546, 0), (f64::NAN, 63.546, 2), (10.0, f64::NAN, 2)] {
        req(&mut m, ma, mm, n);
    }
    for (ma, mma, na, mmb, nb) in [
        (31.773, 63.546, 2, 107.868, 1),
        (10.0, 63.546, 2, 0.0, 1),
        (f64::NAN, 63.546, 2, 107.868, 1),
        (10.0, f64::NAN, 2, 107.868, 1),
        (10.0, 63.546, 2, f64::INFINITY, 1),
        (1e300, 1e-300, 1, 1e300, 255), // the equivalents overflow and underflow
    ] {
        ratio(&mut m, ma, mma, na, mmb, nb);
    }
    for _ in 0..300 {
        let q = rng.signed_log_uniform(-3.0, 8.0);
        let mm = rng.uniform(1.0, 250.0);
        let n = u8::try_from(rng.int(1, 6)).unwrap_or(1);
        dep(&mut m, q, mm, n);
        let ma = rng.log_uniform(-6.0, 4.0);
        req(&mut m, ma, mm, n);
        let mmb = rng.uniform(1.0, 250.0);
        let nb = u8::try_from(rng.int(1, 6)).unwrap_or(1);
        ratio(&mut m, ma, mm, n, mmb, nb);
    }

    // ── utility ──
    let mut xs = edges();
    xs.extend(decades(-323, 308, &["1", "2.5", "9.6485"]).into_iter().step_by(3));
    xs.extend([FARADAY, -FARADAY, 2.5]);
    for _ in 0..200 {
        xs.push(rng.signed_log_uniform(-10.0, 10.0));
    }
    for x in xs {
        let n = charge_to_moles_electrons(bb(x));
        m.row("electrochemistry::charge_to_moles_electrons", Class::B, &[f(x)], &[f(n)]);
        let q = moles_electrons_to_charge(bb(x));
        m.row("electrochemistry::moles_electrons_to_charge", Class::B, &[f(x)], &[f(q)]);
    }
    audit::emit_electrochemistry(&mut m);
    m
}
