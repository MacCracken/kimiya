//! `thermochem.tsv` — the thermochemical tables, reaction sums, van 't Hoff, Shomate heat
//! capacities and enthalpies, and the adiabatic flame temperature (port plan §4, M2b).
//!
//! Species lists (`&[(&str, f64)]`) are `v:n` then n pairs `s:formula f:coefficient`. A `Shomate`
//! is flattened as its seven fields `f:a f:b f:c f:d f:e f:t_min f:t_max`.
//! - Class A: `thermochem::THERMOCHEM_DATA` (`v:76`) and one `THERMOCHEM_DATA[]` row per entry
//!   (input `i:k`; output `s:formula s:name` and the three f64 fields), `SHOMATE_DATA` / `SHOMATE_DATA[]`
//!   (`s:formula` and the seven fields), `lookup_thermochem` (every formula, case / spacing / phase
//!   variants, unknown ones; the five fields or `n:`), `lookup_shomate` and `lookup_shomate_at` (every
//!   species at each segment's ends, midpoints and shared boundaries — the lower segment wins — one
//!   ulp either side, outside the range, specials; the seven fields or `n:`).
//! - Class B (+ − × ÷ only, expected 0 ulp):
//!   - `reaction_enthalpy`, `reaction_gibbs_energy` (**ADR 0002 #1**: the Cyrius functions return
//!     J/mol, so the expected value is the oracle's kJ/mol result × 1000.0, computed here; the row
//!     note says so) and `reaction_entropy` (J/(mol·K), unchanged): the Rust tests' reactions, empty
//!     sides, products only, unknown formulas first / last on either side, zero / −0 / negative /
//!     NaN / ±∞ / huge coefficients (∞ − ∞, ∞ × 0), and seeded lists of 0–6 species per side;
//!   - `Shomate::cp`, `Shomate::enthalpy_antiderivative`: every segment at its range ends, one ulp
//!     either side, midpoint, extrapolation below and above, the reference temperatures, T = 0, −0,
//!     negative, tiny, huge, NaN, ±∞, seeded points; and custom sets (all zero, E only);
//!   - `enthalpy_change_cp`: every species within and across segments (each shared boundary, from
//!     one ulp either side), reversed, equal temperatures, the range ends and one ulp outside them,
//!     the Rust tests' cases, NaN / ±∞ / 0 / negative temperatures, unknown formulas, seeded pairs;
//!   - `adiabatic_flame_temperature` (input `f:ΔH`, the products, `f:T₀`): the transcribed hisab
//!     1.4.0 Newton + bisection make the whole path + − × ÷. Fuels whose products have Shomate data
//!     (H₂, CH₄, CO, C₂H₆, C₂H₄, C₂H₂, C₃H₈, NH₃, and rich CH₄ → CO), stoichiometric, lean (excess O₂ as
//!     a product) and rich (excess H₂), in air (3.76 N₂ per O₂, with and without Ar) and in pure O₂,
//!     at several T₀; seeded mixtures whose ΔH puts the root at a seeded temperature; the moles and
//!     ΔH scaled by 10⁻²⁰⁰ … 10³⁰⁰ (small scales: Newton accepts its guess at once, the bracket product
//!     underflows; large ones: |f| can't get below the 10⁻³ J/mol tolerance, so the bisection
//!     fallback answers); and every Err path in Rust's check order.
//! - Class C: `vant_hoff_k` (one `exp`): each argument over the edges, the Rust tests' values, exp's
//!   overflow and underflow, and a seeded sweep; and (ADR 0006 item 6, M3c-3) the hard rows where the exp
//!   result is subnormal and K₁ scales it back into the normal range, hand-picked and seeded.
//! - ADR 0002 #33 and #34 (the overflow audit's B4): `reaction_enthalpy`, `reaction_gibbs_energy` and
//!   `reaction_entropy` rows (#33: Rust's sums with an unbounded exponent, the final subtraction rounded once, then
//!   ADR 0002 #1's × 1000) and `Shomate::cp` and `Shomate::enthalpy_antiderivative` rows (#34: Rust's polynomial order
//!   with an unbounded exponent, the final addition or product rounded once) are noted where the port's value differs
//!   from Rust's (`audit.rs`), and `audit::emit_thermochem` appends their edge rows with their own seed.
//!
//! **Solver paths.** `flame_replica` re-implements the flame solver from the public API (`SHOMATE_DATA`,
//! `Shomate::cp`, `Shomate::enthalpy_antiderivative`, and hisab 1.4.0's `newton_raphson` and `bisection`
//! loops transcribed) and must give the oracle's bits on every row, or generation fails. It reports
//! which solver produced each result, and the row note carries it (`[newton, k steps]`,
//! `[bisection after newton: <reason>, k halvings]`).

use std::hint::black_box as bb;

use crate::audit;
use crate::grids::{EXP_SUBNORMAL_X, Rng, specials, subnormal_prefactors};
use crate::m_reaction::{edges, res};
use crate::t_thermochem::shomate_fields;
use crate::tsv::{Class, Module, Tok, f, i, s};
use kimiya::element::GAS_CONSTANT;
use kimiya::error::KimiyaError;
use kimiya::thermochem::{
    SHOMATE_DATA, Shomate, THERMOCHEM_DATA, adiabatic_flame_temperature, enthalpy_change_cp,
    lookup_shomate, lookup_shomate_at, lookup_thermochem, reaction_enthalpy, vant_hoff_k,
};

type Species = Vec<(String, f64)>;

fn sp(xs: &[(&str, f64)]) -> Species {
    xs.iter().map(|(f, c)| ((*f).to_string(), *c)).collect()
}

fn borrow(xs: &Species) -> Vec<(&str, f64)> {
    xs.iter().map(|(f, c)| (f.as_str(), *c)).collect()
}

fn species_toks(xs: &Species) -> Vec<Tok> {
    let mut v = vec![Tok::V(xs.len())];
    for (fo, c) in xs {
        v.push(s(fo.clone()));
        v.push(f(*c));
    }
    v
}

fn shomate_toks(x: &Shomate) -> Vec<Tok> {
    shomate_fields(x).iter().map(|&v| f(v)).collect()
}

fn opt_shomate(x: Option<&Shomate>) -> Vec<Tok> {
    x.map_or_else(|| vec![Tok::None], shomate_toks)
}

// ── the flame-solver replica ────────────────────────────────────────────────────────

/// The contiguous run of one species' segments in `SHOMATE_DATA` (Rust's private `Segments`).
struct Segs {
    run: &'static [(&'static str, Shomate)],
    t_min: f64,
    t_max: f64,
}

fn segs(formula: &str) -> Option<Segs> {
    let start = SHOMATE_DATA.iter().position(|(f, _)| *f == formula)?;
    let rest = &SHOMATE_DATA[start..];
    let len = rest.iter().take_while(|(f, _)| *f == formula).count();
    let run = &rest[..len];
    Some(Segs {
        run,
        t_min: run[0].1.t_min,
        t_max: run[len - 1].1.t_max,
    })
}

impl Segs {
    fn cp(&self, t: f64) -> f64 {
        self.run
            .iter()
            .find(|(_, s)| t <= s.t_max)
            .or_else(|| self.run.last())
            .map_or(f64::NAN, |(_, s)| s.cp(t))
    }

    fn enthalpy_change(&self, t1: f64, t2: f64) -> f64 {
        let (lo, hi, sign) = if t1 <= t2 { (t1, t2, 1.0) } else { (t2, t1, -1.0) };
        let n = self.run.len();
        let mut total = 0.0;
        let mut lower = f64::NEG_INFINITY;
        for (k, (_, s)) in self.run.iter().enumerate() {
            let upper = if k + 1 == n { f64::INFINITY } else { s.t_max };
            let a = lo.max(lower);
            let b = hi.min(upper);
            if b > a {
                total += s.enthalpy_antiderivative(b) - s.enthalpy_antiderivative(a);
            }
            if upper >= hi {
                break;
            }
            lower = upper;
        }
        sign * total
    }
}

/// Which solver gave a flame temperature.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Path {
    /// Newton–Raphson, after this many updates.
    Newton(usize),
    /// Bisection, after Newton failed for this reason; this many halvings.
    Bisection(&'static str, usize),
}

/// The flame solver rebuilt from the public API, with hisab 1.4.0's loops transcribed: the result
/// (or the error variant) and, for an `Ok`, the solver path. With `exact`, the port's (ADR 0002 #38, the overflow
/// audit's B6): the bracket test and the bisection's same-sign test compare the signs of f(T₀) and f(T_max) exactly
/// (f_lo < 0 < f_hi or f_hi < 0 < f_lo; both > 0 or both < 0), where Rust tests their product, which underflows to
/// ±0 for a tiny |ΔH| and moles; everything else is Rust's.
fn flame_replica(
    dh: f64,
    products: &[(&str, f64)],
    t0: f64,
    exact: bool,
) -> (Result<f64, &'static str>, Option<Path>) {
    let pos = |x: f64| x > 0.0 && x.is_finite();
    if !pos(t0) {
        return (Err("InvalidTemperature"), None);
    }
    if !dh.is_finite() {
        return (Err("InvalidInput"), None);
    }
    if products.is_empty() {
        return (Err("InvalidInput"), None);
    }
    let mut resolved = Vec::new();
    let mut range_min = f64::NEG_INFINITY;
    let mut range_max = f64::INFINITY;
    for &(formula, moles) in products {
        let Some(sg) = segs(formula) else {
            return (Err("InvalidReaction"), None);
        };
        if !pos(moles) {
            return (Err("InvalidInput"), None);
        }
        range_min = range_min.max(sg.t_min);
        range_max = range_max.min(sg.t_max);
        resolved.push((sg, moles));
    }
    if t0 < range_min || t0 >= range_max {
        return (Err("InvalidTemperature"), None);
    }
    let fun = |t: f64| -> f64 {
        let heat: f64 = resolved.iter().map(|(s, n)| n * s.enthalpy_change(t0, t)).sum();
        -dh - heat
    };
    let dfun = |t: f64| -> f64 {
        let cp: f64 = resolved.iter().map(|(s, n)| n * s.cp(t)).sum();
        -cp
    };
    let (t_lo, t_hi) = (t0, range_max);
    let (fa, fb) = (fun(t_lo), fun(t_hi));
    let bracketed = if exact { (fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0) } else { fa * fb < 0.0 };
    if !bracketed {
        return (Err("ComputationError"), None);
    }
    let total: f64 = resolved.iter().map(|(_, n)| n).sum();
    let guess = (t0 - dh / (total * 40.0)).max(t_lo).min(t_hi);

    // hisab 1.4.0 newton_raphson(f, df, guess, 1e-3, 50)
    let mut x = guess;
    let mut why = "no convergence in 50 iterations";
    for k in 0..50 {
        let fx = fun(x);
        if fx.abs() < 1e-3 {
            if (t_lo..=t_hi).contains(&x) {
                return (Ok(x), Some(Path::Newton(k)));
            }
            why = "root outside the range";
            break;
        }
        let dfx = dfun(x);
        if dfx.abs() < 1e-12 {
            why = "zero derivative";
            break;
        }
        x -= fx / dfx;
    }

    // hisab 1.4.0 bisection(f, t_lo, t_hi, 1e-6, 200)
    let (mut lo, mut hi) = (t_lo, t_hi);
    let (f_lo, f_hi) = (fun(lo), fun(hi));
    let same_sign = if exact { (f_lo > 0.0 && f_hi > 0.0) || (f_lo < 0.0 && f_hi < 0.0) } else { f_lo * f_hi > 0.0 };
    if same_sign {
        return (Err("ComputationError"), None);
    }
    if f_lo > 0.0 {
        std::mem::swap(&mut lo, &mut hi);
    }
    for k in 0..200 {
        let mid = (lo + hi) * 0.5;
        let fm = fun(mid);
        if fm.abs() < 1e-6 || (hi - lo).abs() < 1e-6 {
            return (Ok(mid), Some(Path::Bisection(why, k)));
        }
        if fm < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (Ok((lo + hi) * 0.5), Some(Path::Bisection(why, 200)))
}

/// f(T₀) and f(T_max), the flame balance at the bracket's ends (for inputs Rust's validation accepts): ADR 0002 #38's
/// rows, where their product underflows to a zero.
fn flame_ends(dh: f64, products: &[(&str, f64)], t0: f64) -> (f64, f64) {
    let resolved: Vec<(Segs, f64)> = products.iter().filter_map(|&(fo, n)| segs(fo).map(|sg| (sg, n))).collect();
    let t_max = resolved.iter().map(|(sg, _)| sg.t_max).fold(f64::INFINITY, f64::min);
    let fun = |t: f64| -> f64 {
        let heat: f64 = resolved.iter().map(|(sg, n)| n * sg.enthalpy_change(t0, t)).sum();
        -dh - heat
    };
    (fun(t0), fun(t_max))
}

/// Solver-path tallies for the module comment and port-parity.md.
#[derive(Default)]
pub struct Paths {
    pub newton: usize,
    pub bisection: usize,
    pub err: usize,
}

fn variant(e: &KimiyaError) -> String {
    format!("{e:?}").split('(').next().unwrap_or_default().to_string()
}

/// How the port obtains a flame temperature under ADR 0002 #38 (the note's `<how>`).
pub const FLAME_HOW: &str = "Rust's solver with an exact bracket sign test";

/// A solver path as the row notes write it.
fn path_text(path: Option<Path>) -> String {
    match path {
        Some(Path::Newton(k)) => format!("newton, {k} steps"),
        Some(Path::Bisection(why, k)) => format!("bisection after newton: {why}, {k} halvings"),
        None => String::new(),
    }
}

fn flame(m: &mut Module, paths: &mut Paths, dh: f64, products: &Species, t0: f64) {
    let p = borrow(products);
    let got = adiabatic_flame_temperature(bb(dh), bb(&p), bb(t0));
    let (rep, path) = flame_replica(dh, &p, t0, false);
    let same = match (&got, &rep) {
        (Ok(a), Ok(b)) => a.to_bits() == b.to_bits(),
        (Err(a), Err(b)) => variant(a) == *b,
        _ => false,
    };
    assert!(
        same,
        "flame replica disagrees with the oracle: ΔH {dh:?}, {products:?}, T0 {t0:?}: {got:?} vs {rep:?}"
    );
    let mut ins = vec![f(dh)];
    ins.extend(species_toks(products));
    ins.push(f(t0));
    match path {
        Some(Path::Newton(_)) => paths.newton += 1,
        Some(Path::Bisection(..)) => paths.bisection += 1,
        None => paths.err += 1,
    }
    let extra = path_text(path);
    // ADR 0002 #38: the port's result, Rust's solver with the exact sign tests. It can differ only where Rust's bracket
    // product underflows to −0 although the signs differ: Rust's ComputationError, the port's Ok (the rest of the
    // solver is the same, and the bisection's test cannot refuse after the bracket passed).
    let (port, ppath) = flame_replica(dh, &p, t0, true);
    let changed = match (&rep, &port) {
        (Ok(a), Ok(b)) => a.to_bits() != b.to_bits(),
        (Err(a), Err(b)) => a != b,
        _ => true,
    };
    if changed {
        let (fa, fb) = flame_ends(dh, &p, t0);
        assert!(
            rep == Err("ComputationError") && port.is_ok() && fa * fb == 0.0 && fa != 0.0 && fb != 0.0,
            "ADR 0002 #38: the port's flame differs from Rust's beyond an underflowing bracket: \
             {dh:?} {products:?} {t0:?}"
        );
        let x = port.unwrap_or(f64::NAN);
        let note = format!(
            "ADR 0002 #38: the port returns {FLAME_HOW} ({}) = {x:?} (f:{:016x})",
            path_text(ppath),
            x.to_bits()
        );
        m.row_noted("thermochem::adiabatic_flame_temperature", Class::B, &ins, &[res(got)], &note);
        return;
    }
    let out = res(got);
    if extra.is_empty() {
        m.row("thermochem::adiabatic_flame_temperature", Class::B, &ins, &[out]);
    } else {
        m.row_noted("thermochem::adiabatic_flame_temperature", Class::B, &ins, &[out], &extra);
    }
}

// ── reaction sums ───────────────────────────────────────────────────────────────────

fn sums(m: &mut Module, prod: &Species, reac: &Species) {
    let (p, r) = (borrow(prod), borrow(reac));
    let mut ins = species_toks(prod);
    ins.extend(species_toks(reac));
    audit::sums(m, &ins, &p, &r);
}

fn coefficient(rng: &mut Rng) -> f64 {
    const SMALL: [f64; 10] = [1.0, 2.0, 3.0, 0.5, 1.5, 4.0, 6.0, 2.5, 7.52, 0.25];
    if rng.int(0, 9) < 6 {
        SMALL[usize::try_from(rng.int(0, 9)).unwrap_or(0)]
    } else {
        rng.log_uniform(-2.0, 2.0)
    }
}

fn random_side(rng: &mut Rng, formulas: &[&str], max: i64) -> Species {
    let n = rng.int(0, max);
    (0..n)
        .map(|_| {
            let k = usize::try_from(rng.int(0, i64::try_from(formulas.len()).unwrap_or(1) - 1)).unwrap_or(0);
            (formulas[k].to_string(), coefficient(rng))
        })
        .collect()
}

// ── Shomate ─────────────────────────────────────────────────────────────────────────

fn shomate_temps(x: &Shomate, rng: &mut Rng) -> Vec<f64> {
    let (lo, hi) = (x.t_min, x.t_max);
    let mut t = vec![
        lo,
        hi,
        lo.next_up(),
        lo.next_down(),
        hi.next_up(),
        hi.next_down(),
        0.5 * (lo + hi),
        0.5 * lo,
        1.5 * hi,
        298.15,
        1000.0,
        2000.0,
        3000.0,
        10_000.0,
        50.0,
        1e-300,
        1e300,
        -300.0,
    ];
    t.extend(specials());
    for _ in 0..6 {
        t.push(rng.uniform(lo, hi));
    }
    t
}

fn shomate_rows(m: &mut Module, x: &Shomate, temps: &[f64]) {
    for &t in temps {
        let mut ins = shomate_toks(x);
        ins.push(f(t));
        audit::shomate(m, &ins, x, t);
    }
}

/// The distinct formulas of `SHOMATE_DATA`, in first-appearance order.
fn shomate_species() -> Vec<&'static str> {
    let mut v: Vec<&'static str> = Vec::new();
    for (fo, _) in SHOMATE_DATA {
        if !v.contains(fo) {
            v.push(fo);
        }
    }
    v
}

fn dh_cp(m: &mut Module, formula: &str, t1: f64, t2: f64) {
    let out = res(enthalpy_change_cp(bb(formula), bb(t1), bb(t2)));
    m.row("thermochem::enthalpy_change_cp", Class::B, &[s(formula), f(t1), f(t2)], &[out]);
}

// ── flame mixtures ──────────────────────────────────────────────────────────────────

/// A fuel's combustion: products per mole of fuel at stoichiometric O₂, the O₂ that takes, and
/// the reaction written with H₂O(g) for the oracle's `reaction_enthalpy`.
struct Fuel {
    fuel: &'static str,
    o2: f64,
    products: &'static [(&'static str, f64)],
}

const FUELS: [Fuel; 9] = [
    Fuel { fuel: "H2(g)", o2: 0.5, products: &[("H2O(g)", 1.0)] },
    Fuel { fuel: "CH4(g)", o2: 2.0, products: &[("CO2(g)", 1.0), ("H2O(g)", 2.0)] },
    Fuel { fuel: "CO(g)", o2: 0.5, products: &[("CO2(g)", 1.0)] },
    Fuel { fuel: "C2H6(g)", o2: 3.5, products: &[("CO2(g)", 2.0), ("H2O(g)", 3.0)] },
    Fuel { fuel: "C2H4(g)", o2: 3.0, products: &[("CO2(g)", 2.0), ("H2O(g)", 2.0)] },
    Fuel { fuel: "C2H2(g)", o2: 2.5, products: &[("CO2(g)", 2.0), ("H2O(g)", 1.0)] },
    Fuel { fuel: "C3H8(g)", o2: 5.0, products: &[("CO2(g)", 3.0), ("H2O(g)", 4.0)] },
    Fuel { fuel: "NH3(g)", o2: 0.75, products: &[("N2(g)", 0.5), ("H2O(g)", 1.5)] },
    // rich methane to CO (the Rust test's reaction)
    Fuel { fuel: "CH4(g)", o2: 1.5, products: &[("CO(g)", 1.0), ("H2O(g)", 2.0)] },
];

/// ΔH (J/mol) of `fuel + o2 O₂ → products` from the oracle's table, kJ × 1000.
fn fuel_dh(fu: &Fuel) -> f64 {
    let reac = [(fu.fuel, 1.0), ("O2(g)", fu.o2)];
    reaction_enthalpy(fu.products, &reac).map_or(f64::NAN, |kj| kj * 1000.0)
}

/// The products of `fu` burned at equivalence ratio `phi` (≤ 1: excess O₂ is a product; > 1 is
/// used for H₂ only: excess H₂ is a product), with `n2_per_o2` N₂ and `ar_per_o2` Ar per O₂ fed.
/// Returns (ΔH J/mol of fuel, products).
fn mixture(fu: &Fuel, phi: f64, n2_per_o2: f64, ar_per_o2: f64) -> (f64, Species) {
    let dh = fuel_dh(fu);
    let mut prods: Species = Vec::new();
    let (scale, o2_fed) = if phi <= 1.0 { (1.0, fu.o2 / phi) } else { (1.0 / phi, fu.o2 / phi) };
    for (pf, n) in fu.products {
        prods.push(((*pf).to_string(), n * scale));
    }
    if phi < 1.0 {
        prods.push(("O2(g)".to_string(), o2_fed - fu.o2));
    }
    if phi > 1.0 {
        prods.push((fu.fuel.to_string(), 1.0 - scale));
    }
    if n2_per_o2 > 0.0 {
        // NH₃ already makes N₂: add the air's N₂ to it.
        if let Some(e) = prods.iter_mut().find(|(pf, _)| pf == "N2(g)") {
            e.1 += n2_per_o2 * o2_fed;
        } else {
            prods.push(("N2(g)".to_string(), n2_per_o2 * o2_fed));
        }
    }
    if ar_per_o2 > 0.0 {
        prods.push(("Ar(g)".to_string(), ar_per_o2 * o2_fed));
    }
    (dh * scale, prods)
}

fn scaled(products: &Species, k: f64) -> Species {
    products.iter().map(|(pf, n)| (pf.clone(), n * k)).collect()
}

pub fn emit() -> Module {
    let mut m = Module::new("thermochem");
    let mut rng = Rng::new(0x3201);

    // ── the tables ──
    m.row("thermochem::THERMOCHEM_DATA", Class::A, &[], &[Tok::V(THERMOCHEM_DATA.len())]);
    let row_toks = |d: &kimiya::thermochem::ThermochemData| {
        vec![s(d.formula), s(d.name), f(d.delta_hf_kj), f(d.delta_gf_kj), f(d.s_standard_j)]
    };
    for (k, d) in THERMOCHEM_DATA.iter().enumerate() {
        m.row("thermochem::THERMOCHEM_DATA[]", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &row_toks(d));
    }
    m.row("thermochem::SHOMATE_DATA", Class::A, &[], &[Tok::V(SHOMATE_DATA.len())]);
    for (k, (fo, x)) in SHOMATE_DATA.iter().enumerate() {
        let mut out = vec![s(*fo)];
        out.extend(shomate_toks(x));
        m.row("thermochem::SHOMATE_DATA[]", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &out);
    }

    // ── lookups ──
    let mut names: Vec<String> = Vec::new();
    for d in THERMOCHEM_DATA {
        names.push(d.formula.to_string());
        names.push(d.formula.to_uppercase());
        names.push(d.formula.to_lowercase());
        names.push(format!(" {}", d.formula));
        names.push(format!("{} ", d.formula));
        names.push(d.formula.split('(').next().unwrap_or_default().to_string());
    }
    for fo in shomate_species() {
        names.push(fo.to_string());
    }
    names.extend(
        ["", "XeF6(g)", "H2O(aq)", "H2O", "h2o(l)", "H2O(l)H2O(l)", "C(diamond)", "Ar(s)"]
            .iter()
            .map(|x| (*x).to_string()),
    );
    names.push("H".repeat(300));
    for n in &names {
        let out = lookup_thermochem(bb(n)).map_or_else(|| vec![Tok::None], row_toks);
        m.row("thermochem::lookup_thermochem", Class::A, &[s(n.clone())], &out);
        let out = opt_shomate(lookup_shomate(bb(n)));
        m.row("thermochem::lookup_shomate", Class::A, &[s(n.clone())], &out);
    }
    let species = shomate_species();
    let at = |m: &mut Module, fo: &str, t: f64| {
        let out = opt_shomate(lookup_shomate_at(bb(fo), bb(t)));
        m.row("thermochem::lookup_shomate_at", Class::A, &[s(fo), f(t)], &out);
    };
    for &fo in &species {
        let mut temps = vec![50.0, 7000.0, 298.15, 1000.0, 2000.0, 0.0, -0.0, -300.0];
        temps.extend([f64::NAN, f64::INFINITY, f64::NEG_INFINITY]);
        for (sf, x) in SHOMATE_DATA {
            if *sf == fo {
                let (lo, hi) = (x.t_min, x.t_max);
                temps.extend([lo, hi, lo.next_up(), lo.next_down(), hi.next_up(), hi.next_down()]);
                temps.push(0.5 * (lo + hi));
            }
        }
        for t in temps {
            at(&mut m, fo, t);
        }
    }
    for t in [300.0, 500.0, f64::NAN] {
        at(&mut m, "XeF6(g)", t);
        at(&mut m, "H2O(l)", t);
        at(&mut m, "", t);
    }

    // ── Shomate::cp, Shomate::enthalpy_antiderivative ──
    for (_, x) in SHOMATE_DATA {
        let temps = shomate_temps(x, &mut rng);
        shomate_rows(&mut m, x, &temps);
    }
    let custom = [
        Shomate { a: 0.0, b: 0.0, c: 0.0, d: 0.0, e: 0.0, t_min: 0.0, t_max: 0.0 },
        Shomate { a: 0.0, b: 0.0, c: 0.0, d: 0.0, e: 1.0, t_min: 1.0, t_max: 2.0 },
        Shomate { a: 1e300, b: 1e300, c: 1e300, d: 1e300, e: 1e300, t_min: 298.0, t_max: 6000.0 },
        Shomate { a: f64::NAN, b: 1.0, c: 1.0, d: 1.0, e: 1.0, t_min: 298.0, t_max: 6000.0 },
    ];
    for x in &custom {
        shomate_rows(&mut m, x, &[0.0, -0.0, 1.0, 500.0, 1000.0, 1e-300, f64::INFINITY, f64::NAN]);
    }

    // ── reaction sums ──
    let rxn = |p: &[(&str, f64)], r: &[(&str, f64)]| (sp(p), sp(r));
    let mut cases = vec![
        // the Rust tests
        rxn(&[("CO2(g)", 1.0), ("H2O(l)", 2.0)], &[("CH4(g)", 1.0), ("O2(g)", 2.0)]),
        rxn(&[("XeF6(g)", 1.0)], &[("H2(g)", 1.0)]),
        rxn(&[("NH3(g)", 2.0)], &[("N2(g)", 1.0), ("H2(g)", 3.0)]),
        rxn(&[], &[]),
        rxn(&[("ZnO(s)", 1.0)], &[("Zn(s)", 1.0), ("O2(g)", 0.5)]),
        rxn(&[("CO2(g)", 1.0), ("H2O(g)", 2.0)], &[("CH4(g)", 1.0), ("O2(g)", 2.0)]),
        // one side empty
        rxn(&[("CO2(g)", 1.0)], &[]),
        rxn(&[], &[("CO2(g)", 1.0)]),
        rxn(&[("H2(g)", 1.0)], &[]),
        // unknown formulas: first or last, on either side, both
        rxn(&[("XeF6(g)", 1.0), ("CO2(g)", 1.0)], &[("CH4(g)", 1.0)]),
        rxn(&[("CO2(g)", 1.0), ("XeF6(g)", 1.0)], &[("CH4(g)", 1.0)]),
        rxn(&[("CO2(g)", 1.0)], &[("CH4(g)", 1.0), ("co2(g)", 1.0)]),
        rxn(&[("CO2(g)", 1.0)], &[("", 1.0), ("CH4(g)", 1.0)]),
        rxn(&[("XeF6(g)", 1.0)], &[("XeF6(g)", 1.0)]),
        rxn(&[("CO2(g)", f64::NAN), ("XeF6(g)", 1.0)], &[]),
        // coefficients
        rxn(&[("CO2(g)", 0.0)], &[("CH4(g)", -0.0)]),
        rxn(&[("CO2(g)", -1.0)], &[("CH4(g)", -2.0)]),
        rxn(&[("CO2(g)", f64::NAN)], &[("CH4(g)", 1.0)]),
        rxn(&[("CO2(g)", 1.0)], &[("CH4(g)", f64::NAN)]),
        rxn(&[("CO2(g)", f64::INFINITY)], &[("CH4(g)", 1.0)]),
        rxn(&[("CO2(g)", f64::NEG_INFINITY)], &[("CH4(g)", 1.0)]),
        rxn(&[("CO2(g)", f64::INFINITY)], &[("CO2(g)", f64::INFINITY)]),
        rxn(&[("H2(g)", f64::INFINITY)], &[]),
        rxn(&[("Al2O3(s)", 1e308)], &[]),
        rxn(&[("Al2O3(s)", 1e308)], &[("Al2O3(s)", 1e308)]),
        rxn(&[("Al2O3(s)", 1e308), ("Al2O3(s)", -1e308)], &[]),
        rxn(&[("CO2(g)", 5e-324)], &[("CH4(g)", 1e-300)]),
        rxn(&[("CO2(g)", 1e16), ("CO2(g)", 1.0), ("CO2(g)", -1e16)], &[]),
        rxn(&[("H2(g)", 1.0), ("O2(g)", 1.0), ("N2(g)", -1.0)], &[("Fe(s)", 3.0)]),
    ];
    let formulas: Vec<&str> = THERMOCHEM_DATA.iter().map(|d| d.formula).collect();
    for _ in 0..400 {
        let p = random_side(&mut rng, &formulas, 6);
        let r = random_side(&mut rng, &formulas, 6);
        cases.push((p, r));
    }
    for (p, r) in &cases {
        sums(&mut m, p, r);
    }

    // ── van 't Hoff ──
    let vh = |m: &mut Module, k: f64, dh: f64, t1: f64, t2: f64| {
        let out = res(vant_hoff_k(bb(k), bb(dh), bb(t1), bb(t2)));
        m.row("thermochem::vant_hoff_k", Class::C, &[f(k), f(dh), f(t1), f(t2)], &[out]);
    };
    for x in edges() {
        vh(&mut m, x, -50_000.0, 298.15, 400.0);
        vh(&mut m, 100.0, x, 298.15, 400.0);
        vh(&mut m, 100.0, -50_000.0, x, 400.0);
        vh(&mut m, 100.0, -50_000.0, 298.15, x);
    }
    for (k, dh, t1, t2) in [
        (1.0, -50_000.0, 298.15, 298.15),
        (100.0, -50_000.0, 298.15, 400.0),
        (0.01, 50_000.0, 298.15, 400.0),
        (1.0, -50_000.0, 0.0, 300.0),
        (1.0, -50_000.0, 300.0, 0.0),
        (0.0, -50_000.0, 298.15, 400.0),
        (1.0, 1e8, 298.15, 1000.0), // exponent ≈ 9.6e3: overflows to +∞
        (1.0, -1e8, 298.15, 1000.0), // underflows to 0
        (1.0, 0.0, 298.15, 1000.0),
        (1.0, -0.0, 298.15, 1000.0),
        (1e300, -1e5, 300.0, 299.0),
        (5e-324, 1e5, 300.0, 299.0),
        (2.5, -1e5, 1e-300, 1e300),
    ] {
        vh(&mut m, k, dh, t1, t2);
    }
    for _ in 0..600 {
        let k = rng.log_uniform(-10.0, 10.0);
        let dh = rng.signed_log_uniform(2.0, 6.5);
        let t1 = rng.uniform(100.0, 2000.0);
        let t2 = rng.uniform(100.0, 2000.0);
        vh(&mut m, k, dh, t1, t2);
    }
    // ADR 0006's hard rows (M3c-3): the exponent −(ΔH/R)(1/T₂ − 1/T₁) placed on subnormal-exp targets x
    // (ΔH = −x·R/(1/T₂ − 1/T₁)) and K₁ scaling e^x back into the normal range. With the stdlib's exp these rows
    // were up to 12 ulp from the oracle on x86_64 (M3b); only a correctly rounded exp matches the corrected oracle.
    let (t1, t2) = (300.0, 299.0);
    let d = 1.0 / t2 - 1.0 / t1;
    for k in subnormal_prefactors() {
        for x in EXP_SUBNORMAL_X {
            vh(&mut m, k, -x * GAS_CONSTANT / d, t1, t2);
        }
    }
    let mut hard = Rng::new(0x3c3c);
    for _ in 0..300 {
        let x = hard.uniform(-745.2, -708.0);
        let k = hard.log_uniform(200.0, 308.0);
        let t1 = hard.uniform(100.0, 2000.0);
        let t2 = hard.uniform(100.0, 2000.0);
        let d = 1.0 / t2 - 1.0 / t1;
        vh(&mut m, k, -x * GAS_CONSTANT / d, t1, t2);
    }

    // ── enthalpy_change_cp ──
    for &fo in &species {
        let sg = segs(fo).unwrap_or_else(|| panic!("{fo}: no segments"));
        let (lo, hi) = (sg.t_min, sg.t_max);
        let mut pairs = vec![
            (lo, hi),
            (hi, lo),
            (lo, lo),
            (hi, hi),
            (lo.next_down(), hi),
            (lo, hi.next_up()),
            (hi.next_up(), lo),
            (lo.next_up(), hi.next_down()),
            (298.15, 1000.0),
            (298.15, 2000.0),
            (298.15, 3000.0),
            (300.0, 600.0),
            (600.0, 300.0),
            (500.0, 500.0),
        ];
        for (_, x) in sg.run.iter().take(sg.run.len() - 1) {
            let b = x.t_max; // a shared boundary
            pairs.extend([
                (b, b),
                (b.next_down(), b.next_up()),
                (b.next_up(), b.next_down()),
                (lo, b),
                (b, hi),
                (b.next_down(), b),
                (b, b.next_up()),
                (0.5 * (lo + b), 0.5 * (b + hi)),
                (0.5 * (b + hi), 0.5 * (lo + b)),
            ]);
        }
        for x in edges() {
            pairs.push((x, 0.5 * (lo + hi)));
            pairs.push((0.5 * (lo + hi), x));
        }
        for _ in 0..12 {
            pairs.push((rng.uniform(lo, hi), rng.uniform(lo, hi)));
        }
        for (t1, t2) in pairs {
            dh_cp(&mut m, fo, t1, t2);
        }
    }
    for (fo, t1, t2) in [
        ("XeF6(g)", 300.0, 600.0),
        ("XeF6(g)", f64::NAN, 600.0),
        ("XeF6(g)", 0.0, 1e9),
        ("H2O(l)", 300.0, 350.0),
        ("", 300.0, 350.0),
        ("CO2(g)", 0.0, 600.0),
        ("CO2(g)", 250.0, 600.0),
        ("NH3(g)", 300.0, 1500.0),
        ("N2(g)", 300.0, 6500.0),
        ("Br2(g)", 298.15, 400.0),
        ("N2(g)", 100.0, 6000.0),
        ("N2(g)", f64::NAN, 600.0),
        ("N2(g)", 300.0, f64::NAN),
        ("N2(g)", f64::INFINITY, 600.0),
        ("N2(g)", 300.0, f64::INFINITY),
        ("N2(g)", -300.0, 600.0),
        ("CO2(g)", 1000.0, 1500.0),
        ("N2(g)", 400.0, 2500.0),
        ("O2(g)", 298.15, 3000.0),
        ("H2O(g)", 298.15, 2500.0),
        ("CO2(g)", 300.0, 600.0),
        ("CO2(g)", 500.0, 500.0),
        ("N2(g)", 600.0, 300.0),
    ] {
        dh_cp(&mut m, fo, t1, t2);
    }

    // ── adiabatic flame temperature ──
    let mut paths = Paths::default();
    // the Rust tests
    let air = |x: &[(&str, f64)]| sp(x);
    flame(&mut m, &mut paths, -241_820.0, &air(&[("H2O(g)", 1.0), ("N2(g)", 1.88)]), 298.15);
    let ch4 = air(&[("CO2(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 7.52)]);
    flame(&mut m, &mut paths, -802_300.0, &ch4, 298.15);
    flame(&mut m, &mut paths, -100_000.0, &air(&[("XeF6(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -100_000.0, &air(&[("CO2(g)", 1.0)]), 0.0);
    let co2n2 = air(&[("CO2(g)", 1.0), ("N2(g)", 3.76)]);
    flame(&mut m, &mut paths, f64::NAN, &co2n2, 298.15);
    flame(&mut m, &mut paths, f64::NEG_INFINITY, &co2n2, 298.15);
    let co2 = air(&[("CO2(g)", 1.0)]);
    flame(&mut m, &mut paths, -100_000.0, &co2, f64::NAN);
    flame(&mut m, &mut paths, -100_000.0, &co2, f64::INFINITY);
    flame(&mut m, &mut paths, -100_000.0, &Vec::new(), 298.15);
    for moles in [0.0, -1.0, f64::NAN, f64::INFINITY, -0.0, 5e-324] {
        flame(&mut m, &mut paths, -100_000.0, &air(&[("CO2(g)", moles)]), 298.15);
    }
    flame(&mut m, &mut paths, 100_000.0, &co2, 298.15);
    flame(&mut m, &mut paths, 0.0, &co2, 298.15);
    flame(&mut m, &mut paths, -0.0, &co2, 298.15);
    flame(&mut m, &mut paths, -200_000.0, &air(&[("NH3(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -241_820.0, &air(&[("H2O(g)", 1.0), ("SO2(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -100_000.0, &co2, 250.0);
    let rich = air(&[("CO(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 5.64)]);
    flame(&mut m, &mut paths, -519_320.0, &rich, 298.15);
    // Err order: the first bad product in list order decides; every check before the list.
    flame(&mut m, &mut paths, -1e5, &air(&[("XeF6(g)", 1.0), ("CO2(g)", -1.0)]), 298.15);
    flame(&mut m, &mut paths, -1e5, &air(&[("CO2(g)", -1.0), ("XeF6(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -1e5, &air(&[("CO2(g)", 1.0), ("N2(g)", 1.0), ("XeF6(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -1e5, &air(&[("H2O(l)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, f64::NAN, &Vec::new(), f64::NAN);
    flame(&mut m, &mut paths, f64::NAN, &Vec::new(), 298.15);
    flame(&mut m, &mut paths, f64::NAN, &air(&[("XeF6(g)", 1.0)]), 298.15);
    flame(&mut m, &mut paths, -1e5, &air(&[("XeF6(g)", 1.0)]), -1.0);
    // T₀ against the products' common range: [298, 6000) for CO2 + N2; [298, 1000) with H2/Cl2.
    for t0 in [298.0, 298.0_f64.next_down(), 297.0, 5999.0, 6000.0, 6000.0_f64.next_down(), 6500.0, 1e-300] {
        flame(&mut m, &mut paths, -1e5, &co2n2, t0);
    }
    let cl = air(&[("HCl(g)", 2.0), ("Cl2(g)", 0.1)]);
    for t0 in [298.15, 999.0, 1000.0, 1000.0_f64.next_down(), 1200.0] {
        flame(&mut m, &mut paths, -184_620.0, &cl, t0);
        flame(&mut m, &mut paths, -1000.0, &cl, t0);
    }

    // fuels × equivalence ratio × oxidiser × T₀
    let t0s = [298.15, 300.0, 400.0, 500.0, 700.0, 1000.0];
    for fu in &FUELS {
        let phis: &[f64] = if fu.fuel == "H2(g)" {
            &[0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0]
        } else {
            &[0.5, 0.7, 0.8, 0.9, 1.0]
        };
        for &phi in phis {
            for (n2, ar) in [(3.76, 0.0), (3.727, 0.0446), (0.0, 0.0), (1.88, 0.0)] {
                let (dh, prods) = mixture(fu, phi, n2, ar);
                for &t0 in &t0s {
                    flame(&mut m, &mut paths, dh, &prods, t0);
                }
            }
        }
    }

    // scaled: tiny (Newton accepts its guess; the bracket product underflows) and huge (|f| can't
    // reach the 1e-3 J/mol tolerance: the bisection fallback)
    let scales = [
        1e-200, 1e-170, 1e-160, 1e-100, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 10.0, 1e3, 1e5, 1e6,
        3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 1e10, 1e12, 1e15, 1e20, 1e50, 1e100, 1e200, 1e300,
    ];
    let base: [(f64, Species); 4] = [
        (-802_300.0, ch4.clone()),
        (-241_820.0, sp(&[("H2O(g)", 1.0), ("N2(g)", 1.88)])),
        (-519_320.0, rich.clone()),
        (-282_980.0, sp(&[("CO2(g)", 1.0), ("N2(g)", 1.88)])),
    ];
    for (dh, prods) in &base {
        for &k in &scales {
            flame(&mut m, &mut paths, dh * k, &scaled(prods, k), 298.15);
            flame(&mut m, &mut paths, dh * k, &scaled(prods, k), 500.0);
        }
    }

    // the bisection sweep: each fuel at φ = 0.8 and 1 in air, scaled by 10¹¹ … 10²⁴, where Newton's
    // |f| < 1e-3 J/mol is reached only when f happens to round to (near) zero
    for fu in &FUELS {
        for phi in [0.8, 1.0] {
            let (dh, prods) = mixture(fu, phi, 3.76, 0.0);
            for k in [1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e21, 1e24] {
                flame(&mut m, &mut paths, dh * k, &scaled(&prods, k), 298.15);
            }
        }
    }

    // seeded mixtures: 1–5 products, ΔH chosen so that the root sits at a seeded temperature
    // (heat to T* from the oracle's enthalpy_change_cp), then nudged
    let wide: Vec<&str> = species
        .iter()
        .copied()
        .filter(|fo| segs(fo).is_some_and(|sg| sg.t_max >= 3400.0 && sg.t_min <= 298.15))
        .collect();
    for _ in 0..300 {
        let n = rng.int(1, 5);
        let mut prods: Species = Vec::new();
        for _ in 0..n {
            let k = usize::try_from(rng.int(0, i64::try_from(wide.len()).unwrap_or(1) - 1)).unwrap_or(0);
            prods.push((wide[k].to_string(), rng.log_uniform(-1.0, 1.5)));
        }
        let t0 = if rng.int(0, 3) == 0 { rng.uniform(298.15, 1500.0) } else { 298.15 };
        let t_star = rng.uniform(t0 + 1.0, 3300.0);
        let heat: f64 = prods
            .iter()
            .map(|(pf, nn)| nn * enthalpy_change_cp(pf, t0, t_star).unwrap_or(f64::NAN))
            .sum();
        let dh = -heat * (1.0 + 0.01 * (rng.unit() - 0.5));
        let k = if rng.int(0, 2) == 0 { rng.log_uniform(6.0, 20.0) } else { 1.0 };
        flame(&mut m, &mut paths, dh * k, &scaled(&prods, k), t0);
    }
    // a long product list (more than any fixed per-call table): every Shomate species once
    let all: Species = species.iter().map(|fo| ((*fo).to_string(), 0.5)).collect();
    let many: Species = (0..40_u32)
        .map(|k| {
            let w = wide[usize::try_from(k).unwrap_or(0) % wide.len()];
            (w.to_string(), 0.1 + 0.05 * f64::from(k))
        })
        .collect();
    for dh in [-1e5, -1e6, -2e6, -1e7, -1e12] {
        flame(&mut m, &mut paths, dh, &all, 298.15);
        flame(&mut m, &mut paths, dh, &many, 298.15);
        flame(&mut m, &mut paths, dh, &many, 600.0);
    }
    eprintln!(
        "thermochem: adiabatic_flame_temperature rows: {} newton, {} bisection, {} Err",
        paths.newton, paths.bisection, paths.err
    );
    audit::emit_thermochem(&mut m);
    emit_flame_b6(&mut m, &mut paths);
    m
}

/// 10^k as a double (the decimal literal's correctly rounded value, by parsing).
fn ten(k: i32) -> f64 {
    format!("1e{k}").parse().unwrap_or(f64::NAN)
}

/// The overflow audit's B6 (ADR 0002 #38), appended with their own seed so the rows before keep theirs: flame rows whose
/// bracket product f(T₀)·f(T_max) underflows to −0 while the signs differ (Rust's ComputationError, the port's Ok) and
/// their neighbours — the survey's example and its endothermic mirror (Err in both); real systems scaled across every
/// binade from 10⁻³²⁰ to 10⁻¹⁵⁰ (below about 10⁻¹⁶⁹ the product underflows), at two T₀; for each, the scale where
/// |f(T₀)·f(T_max)| crosses 2⁻¹⁰⁷⁵ (found by bisection on the scale), ulp by ulp around it; a subnormal ΔH against moles
/// from 10⁻⁹ to 10², both signs; an endothermic system at the same scales (both ends negative, their product +0: Err in
/// both); seeded fuel mixtures at seeded tiny scales. (The bisection fallback can't be reached with an underflowing
/// bracket: |f| at the guess is then far below Newton's 10⁻³ J/mol, which accepts it at once.)
fn emit_flame_b6(m: &mut Module, paths: &mut Paths) {
    let mut rng = Rng::new(0x3806);
    let co2 = sp(&[("CO2(g)", 1.0)]);
    flame(m, paths, -1e-320, &scaled(&co2, 1e-300), 298.15);
    flame(m, paths, 1e-320, &scaled(&co2, 1e-300), 298.15);
    let base: [(f64, Species); 6] = [
        (-802_300.0, sp(&[("CO2(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 7.52)])),
        (-241_820.0, sp(&[("H2O(g)", 1.0), ("N2(g)", 1.88)])),
        (-519_320.0, sp(&[("CO(g)", 1.0), ("H2O(g)", 2.0), ("N2(g)", 5.64)])),
        (-282_980.0, sp(&[("CO2(g)", 1.0), ("N2(g)", 1.88)])),
        (-393_510.0, co2.clone()),
        (-2_043_900.0, sp(&[("CO2(g)", 3.0), ("H2O(g)", 4.0), ("N2(g)", 18.8)])),
    ];
    let exps = [-320, -310, -300, -290, -280, -260, -240, -220, -200, -190, -180, -175, -172, -170, -169, -168, -165,
        -160, -155, -150];
    for (dh, prods) in &base {
        for &k in &exps {
            let sc = ten(k);
            for t0 in [298.15, 500.0] {
                flame(m, paths, dh * sc, &scaled(prods, sc), t0);
            }
        }
        // the scale where |f(T₀)·f(T_max)| crosses 2⁻¹⁰⁷⁵: below it the product is −0 (Rust's Err), above −2⁻¹⁰⁷⁴ or
        // more (Rust's Ok); the rows at the two sides, ulp by ulp
        let pv = |sc: f64| {
            let sp_ = scaled(prods, sc);
            let (fa, fb) = flame_ends(dh * sc, &borrow(&sp_), 298.15);
            fa * fb
        };
        let (mut lo, mut hi) = (ten(-200), ten(-150));
        for _ in 0..200 {
            let mid = lo.sqrt() * hi.sqrt();
            if mid <= lo || mid >= hi {
                break;
            }
            if pv(mid) == 0.0 { lo = mid } else { hi = mid }
        }
        for sc in [lo, hi] {
            for d in -2_i64..=2 {
                let x = f64::from_bits(sc.to_bits().wrapping_add_signed(d));
                flame(m, paths, dh * x, &scaled(prods, x), 298.15);
            }
        }
    }
    // a subnormal ΔH against moles from 10⁻⁹ to 10², both signs
    for kk in [1_u64, 3, 1 << 20, (1 << 52) - 1] {
        let dh = f64::from_bits(kk);
        for e in [-9, -7, -6, -5, -3, 0, 2] {
            flame(m, paths, -dh, &scaled(&co2, ten(e)), 298.15);
            flame(m, paths, dh, &scaled(&co2, ten(e)), 298.15);
        }
    }
    // an endothermic system at tiny scales: both ends negative, their product +0 (Err in both)
    for k in [-300, -250, -200, -170, -160] {
        let sc = ten(k);
        flame(m, paths, 241_820.0 * sc, &scaled(&base[1].1, sc), 298.15);
    }
    // seeded: a fuel mixture at a seeded tiny scale and T₀
    for _ in 0..80 {
        let fu = &FUELS[usize::try_from(rng.int(0, 8)).unwrap_or(0)];
        let phi = [0.7, 0.9, 1.0, 1.5][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        let phi = if phi > 1.0 && fu.fuel != "H2(g)" { 1.0 } else { phi }; // rich only for H₂ (its excess has data)
        let (dh, prods) = mixture(fu, phi, 3.76, 0.0);
        let sc = rng.log_uniform(-320.0, -150.0);
        let t0 = [298.15, 300.0, 500.0, 1000.0][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        flame(m, paths, dh * sc, &scaled(&prods, sc), t0);
    }
}
