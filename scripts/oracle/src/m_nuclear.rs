//! `nuclear.tsv` — radioactive decay, binding energy, the semi-empirical mass formula, Q-values, the two-step decay
//! chain and the isotope table (port plan §4, M4c).
//!
//! Rust's `u32` (z, n, a) are `i:` tokens; a slice is `v:n` and n `f:` tokens; `decay_chain_two_step`'s tuple is its
//! three floats `f:N_A f:N_B f:N_C`.
//! - Class A: `nuclear::ISOTOPES` (`v:18`) and one `ISOTOPES[]` row per entry (input `i:k`; output `s:symbol i:z i:a
//!   f:half_life_s s:decay_mode`); `nuclear::lookup_isotope` (every table symbol, its lower- and upper-case forms, a
//!   leading / trailing space, every proper prefix, `""`, `"Xx-999"` and near-miss labels; the five fields or `n:`);
//!   `is_secular_equilibrium` (the edges against each other, λ_p at and one and two ulps either side of Rust's
//!   0.01·λ_d, the Rust tests' values, seeded).
//! - Class B (+ − × ÷, √ and cbrt, expected 0 ulp): `decay_constant` (Rust's `LN_2 / t½`: the edges, every isotope,
//!   every decade from the smallest subnormal to DBL_MAX, t½ = 1 giving LN_2 exactly, +∞ the stable isotope's λ = 0),
//!   `activity`, `mass_defect` and `binding_energy_per_nucleon` (u32 probes from 0 to u32::MAX on each side, the
//!   `checked_add` overflow and A = 0 Err paths in Rust's order, the edges of the mass), `binding_energy_mev`,
//!   `semi_empirical_binding_energy` (every Z ≤ A ≤ 3Z with A ≤ 300 (port plan §4), Z = 0, the u32 extremes; **ADR 0002
//!   #3**: the oracle's NaN for z > a or a = 0, the rows noted with the port's `Err(InvalidInput)`, which the harness
//!   maps), `q_value` (empty slices and the −0.0 folds, the Rust tests' D–T fusion, ∞ − ∞, overflow, cancellation,
//!   seeded lists).
//! - Class C: `atoms_remaining` (one `exp`), `atoms_after_half_lives` (`2.0.powf(−n)`, glibc's exp2 in the release
//!   oracle; the port's correctly rounded pow(2, −n)), `time_for_fraction` (one `ln`), `decay_chain_two_step` (`exp` ×2,
//!   `exp_m1` ×1 or ×2): each argument over the edges and every Err path in Rust's order, the Rust tests' values, the
//!   exponent across exp's range, **on purpose the rows where an exp or pow result is subnormal and a prefactor scales it
//!   back into the normal range** (N₀ up to DBL_MAX; 2⁻ⁿ for n from 1022 to 1075.5; the chain's e^(−λ_min·t) times
//!   N₀·λ₁), 2⁻ⁿ's tie at n = 1075 and its over / underflow; for the chain the equal-λ branch (λ₁ = λ₂ exactly), λ₁ ≈ λ₂
//!   at relative 10⁻¹⁵ … 10⁻⁶ either side of Rust's 10⁻¹² threshold, and **pairs on the threshold itself** (δ equal to
//!   fl(10⁻¹²·λ_max), and where 10⁻¹² ± 1 ulp would move the branch), every isotope pair, seeded sweeps.
//! - ADR 0002 #32 (the overflow audit's B3): `binding_energy_per_nucleon` rows are noted where the port's Rust sequence
//!   with an unbounded exponent (the final division rounded once) differs from Rust's value (`audit.rs`: none of the
//!   rows above), and `audit::emit_nuclear` appends its edge rows (Δm·931.494… past DBL_MAX, A bringing it back) with its
//!   own seed.
//! - ADR 0002 #33 (the overflow audit's B4): `q_value` rows are noted where Rust's two sums and their difference with an
//!   unbounded exponent (the final subtraction rounded once) differ from Rust's value (`audit.rs`), and
//!   `audit::emit_nuclear_b4` appends its edge rows (sums past DBL_MAX that cancel) with its own seed.
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`).

use std::hint::black_box as bb;

use crate::audit;
use crate::exact::{cr_exp2_frac, decompose, round_scaled};
use crate::grids::{EXP_SUBNORMAL_X, Rng, decades, pow2, subnormal_prefactors};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, b, err, f, i, s};
use kimiya::error::KimiyaError;
use kimiya::nuclear::{
    ISOTOPES, Isotope, activity, atoms_after_half_lives, atoms_remaining, binding_energy_mev,
    decay_chain_two_step, decay_constant, is_secular_equilibrium, lookup_isotope,
    mass_defect, semi_empirical_binding_energy, time_for_fraction,
};

/// An isotope's five fields, in declaration order.
pub fn isotope_toks(r: &Isotope) -> Vec<Tok> {
    vec![s(r.symbol), i(r.z), i(r.a), f(r.half_life_s), s(r.decay_mode)]
}

/// The u32 probes of the integer arguments: the small values, real nuclei's Z and N, the u8 / u16 / i32 edges and
/// u32::MAX.
const U32S: [u32; 16] = [0, 1, 2, 3, 6, 8, 26, 30, 92, 118, 146, 255, 256, 65_535, 1 << 31, u32::MAX];

/// A Result<(f64, f64, f64)> as tokens: `f:N_A f:N_B f:N_C`, or `E:Variant`.
fn triple(r: Result<(f64, f64, f64), KimiyaError>) -> Vec<Tok> {
    match r {
        Ok((a, b_, c)) => vec![f(a), f(b_), f(c)],
        Err(e) => vec![err(&e)],
    }
}

/// `x` moved `d` ulps (up for d > 0).
fn ulps(x: f64, d: i32) -> f64 {
    let mut y = x;
    for _ in 0..d.unsigned_abs() {
        y = if d > 0 { y.next_up() } else { y.next_down() };
    }
    y
}

/// (λ₁, λ₂) pairs on Rust's equal-rate threshold `|λ₂ − λ₁| <= 1e-12 · max(λ₁, λ₂)`, both orientations. For each target
/// δ, the λ_max within 96 ulps of δ / 10⁻¹² whose λ_min = λ_max − δ is exact (so |λ₂ − λ₁| is δ) and where δ is at the
/// threshold: fl(10⁻¹²·λ_max) = δ (the equal branch, by `<=`), or 10⁻¹² one ulp up or down would move the branch. These
/// rows pin the tolerance to the ulp.
fn threshold_pairs() -> Vec<(f64, f64)> {
    let tol = 1e-12_f64;
    let (up, down) = (tol.next_up(), tol.next_down());
    let mut deltas = vec![pow2(-40), pow2(-20), pow2(-60), 1.5 * pow2(-30), 1.25 * pow2(-45), 3.0 * pow2(-50)];
    // seeded δ = m·2^e with m below 2^13, so λ_max − δ can be exact (λ_max ≈ 2^40·δ: its ulp is about 2^e)
    let mut rng = Rng::new(0x4c0d);
    for _ in 0..40 {
        let mant = rng.int(2049, 8191);
        let e = i32::try_from(rng.int(-900, 900)).unwrap_or(0);
        #[allow(clippy::cast_precision_loss)]
        deltas.push(mant as f64 * pow2(e));
    }
    let mut out = Vec::new();
    for delta in deltas {
        let mut lam = ulps(delta / tol, -96);
        for _ in 0..193 {
            let lo = lam - delta;
            if lo > 0.0 && lam - lo == delta {
                let thr = tol * lam;
                let on = thr == delta;
                let by_up = thr < delta && delta <= up * lam;
                let by_down = down * lam < delta && delta <= thr;
                if on || by_up || by_down {
                    out.push((lo, lam));
                    out.push((lam, lo));
                }
            }
            lam = lam.next_up();
        }
    }
    out
}

/// ADR 0002 #13: the port's N₀/2ⁿ where Rust's 2⁻ⁿ is not a normal double (finite n > 1022 or n ≤ −1024), with how it is
/// obtained; `None` elsewhere (Rust's `initial * 2.0.powf(-n)` stays). ±0, ±∞ and NaN N₀ give N₀ (scaleB). Otherwise the
/// correctly rounded N₀/2ⁿ: for an integer n exactly in integers, for N₀ = ±2^j with j − n exact CR(2^(j − n)) (glibc's
/// exp2 here; `validate.py` proves it against MPFR), and since M4d for a fractional n (n = k + f) by `cr_exp2_frac`,
/// 2⁻ᶠ to 320 bits (M4c's was N₀·CR(2⁻ᶠ)·2⁻ᵏ rounded once: faithful).
fn half_lives_exact(n0: f64, n: f64) -> Option<(f64, &'static str)> {
    if !n.is_finite() || !(n > 1022.0 || n <= -1024.0) {
        return None;
    }
    if n0.is_nan() || n0 == 0.0 || n0.is_infinite() {
        return Some((n0, "N₀ itself (scaleB)"));
    }
    let neg = n0 < 0.0;
    if n > 2200.0 {
        return Some((if neg { -0.0 } else { 0.0 }, "the correctly rounded N₀/2ⁿ"));
    }
    if n < -2200.0 {
        return Some((if neg { f64::NEG_INFINITY } else { f64::INFINITY }, "the correctly rounded N₀/2ⁿ"));
    }
    let (mant, e0) = decompose(n0);
    let k = n.floor();
    let fr = n - k;
    #[allow(clippy::cast_possible_truncation)]
    let ki = k as i64;
    if mant.is_power_of_two() {
        let j = e0 + i64::from(mant.trailing_zeros());
        #[allow(clippy::cast_precision_loss)]
        let jf = j as f64;
        let y = jf - n;
        let bv = y - jf; // TwoSum(jf, −n): j − n is exact when its error term is 0
        let err = (jf - (y - bv)) + (-n - bv);
        if err == 0.0 {
            let v = y.exp2();
            return Some((if neg { -v } else { v }, "the correctly rounded N₀/2ⁿ (2^(j − n))"));
        }
    }
    if fr == 0.0 {
        return Some((round_scaled(u128::from(mant), e0 - ki, neg), "the correctly rounded N₀/2ⁿ"));
    }
    Some((cr_exp2_frac(mant, e0 - ki, fr, neg), "the correctly rounded N₀/2ⁿ (2⁻ᶠ to 320 bits)"))
}

/// ADR 0002 #14: λ_p < λ_d/100 exactly (100·λ_p against λ_d in integers), Rust's NaN / ±∞ / ±0 semantics kept.
fn secular_exact(lp: f64, ld: f64) -> bool {
    if lp.is_nan() || ld.is_nan() {
        return false;
    }
    if lp.is_infinite() || ld.is_infinite() {
        return lp < ld;
    }
    let sgn = |x: f64| if x == 0.0 { 0 } else if x < 0.0 { -1 } else { 1 };
    let (sp, sd) = (sgn(lp), sgn(ld));
    if sp != sd || sp == 0 {
        return sp < sd;
    }
    let (mp, ep) = decompose(lp);
    let (md, ed) = decompose(ld);
    let a = u128::from(mp) * 100;
    let b_ = u128::from(md);
    let la = 128 - i64::from(a.leading_zeros()) + ep;
    let lb = 128 - i64::from(b_.leading_zeros()) + ed;
    let ord = if la == lb {
        if ep >= ed {
            (a << u32::try_from(ep - ed).unwrap_or(0)).cmp(&b_)
        } else {
            a.cmp(&(b_ << u32::try_from(ed - ep).unwrap_or(0)))
        }
    } else {
        la.cmp(&lb)
    };
    if sp > 0 { ord == std::cmp::Ordering::Less } else { ord == std::cmp::Ordering::Greater }
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("nuclear");
    let mut rng = Rng::new(0x4c01);
    let lambdas: Vec<f64> = ISOTOPES.iter().map(|r| decay_constant(r.half_life_s).unwrap_or(f64::NAN)).collect();

    // ── decay_constant ──
    let dc = |m: &mut Module, t: f64| {
        let out = res(decay_constant(bb(t)));
        m.row("nuclear::decay_constant", Class::B, &[f(t)], &[out]);
    };
    let mut hs = edges();
    hs.extend([std::f64::consts::LN_2, 100.0, 1.0_f64.next_up(), 1.0_f64.next_down(), f64::MIN_POSITIVE / 2.0]);
    hs.extend(ISOTOPES.iter().map(|r| r.half_life_s));
    hs.extend(decades(-323, 308, &["1", "3", "6.931471805599453"]));
    for _ in 0..300 {
        hs.push(rng.log_uniform(-320.0, 308.0));
    }
    for t in hs {
        dc(&mut m, t);
    }

    // ── atoms_remaining ──
    let ar = |m: &mut Module, n0: f64, l: f64, t: f64| {
        let out = atoms_remaining(bb(n0), bb(l), bb(t));
        m.row("nuclear::atoms_remaining", Class::C, &[f(n0), f(l), f(t)], &[f(out)]);
    };
    for x in edges() {
        ar(&mut m, x, 0.1, 5.0);
        ar(&mut m, 1000.0, x, 5.0);
        ar(&mut m, 1000.0, 0.1, x);
    }
    // the Rust tests' values: one half-life of t½ = 100
    let l100 = decay_constant(100.0).unwrap_or(f64::NAN);
    for (n0, l, t) in [
        (1000.0, l100, 100.0),
        (1000.0, l100, 0.0),
        (1000.0, 0.0, f64::INFINITY), // 0·∞ = NaN in the exponent
        (f64::INFINITY, 1e300, 1e300), // ∞·exp(−∞) = ∞·0 = NaN
        (0.0, -1e300, 1e300),          // 0·exp(+∞) = NaN
        (f64::MAX, -1.0, 1.0),         // overflows to +∞
        (1e300, 1.0, 745.0),           // a subnormal exp result scaled into the normal range
    ] {
        ar(&mut m, n0, l, t);
    }
    // every isotope's λ, over a few of its half-lives
    for (r, &l) in ISOTOPES.iter().zip(&lambdas) {
        for k in [0.5, 1.0, 2.0, 10.0, 100.0, 1000.0] {
            ar(&mut m, 6.022e23, l, k * r.half_life_s);
        }
        ar(&mut m, 1.0, l, 3.155_76e7);
    }
    // λ·t on targets across exp's range, N₀ over many decades
    let targets = [
        -800.0, -746.0, -745.2, -745.14, -745.13, -745.1, -745.0, -744.0, -740.0, -730.0, -720.25, -710.0, -708.5,
        -708.4, -708.39, -708.0, -700.0, -100.0, -1.0, -1e-10, 0.0, 1e-10, 1.0, 100.0, 700.0, 709.7, 709.79, 710.0,
    ];
    for n0 in [1.0, 1e-300, 0.5, 1e13, 1e300, f64::MAX, -2.0, f64::from_bits(1)] {
        for x in targets {
            ar(&mut m, n0, 1.0, -x);
            ar(&mut m, n0, 0.25, -4.0 * x);
        }
    }
    // ADR 0006 item 6: e^(−λt) subnormal, N₀ scaling it back
    for x in EXP_SUBNORMAL_X {
        for p in subnormal_prefactors() {
            ar(&mut m, p, 1.0, -x);
            ar(&mut m, p, 2.0, -x / 2.0);
        }
    }
    for _ in 0..300 {
        let x = rng.uniform(-745.2, -708.0);
        let l = rng.log_uniform(-3.0, 3.0);
        let n0 = rng.log_uniform(0.0, 308.0);
        ar(&mut m, n0, l, -x / l);
    }
    for _ in 0..300 {
        let n0 = rng.log_uniform(0.0, 25.0);
        let l = rng.log_uniform(-18.0, 0.0);
        let t = rng.log_uniform(0.0, 18.0);
        ar(&mut m, n0, l, t);
    }
    for _ in 0..200 {
        let n0 = rng.signed_log_uniform(-300.0, 300.0);
        let l = rng.signed_log_uniform(-10.0, 10.0);
        let t = rng.log_uniform(-10.0, 10.0);
        ar(&mut m, n0, l, t);
    }

    // ── activity ──
    let act = |m: &mut Module, l: f64, n: f64| {
        let out = activity(bb(l), bb(n));
        m.row("nuclear::activity", Class::B, &[f(l), f(n)], &[f(out)]);
    };
    for x in edges() {
        for y in edges() {
            act(&mut m, x, y);
        }
    }
    act(&mut m, 0.01, 1e6);
    for &l in &lambdas {
        act(&mut m, l, 6.022e23);
    }
    for _ in 0..200 {
        let l = rng.log_uniform(-20.0, 5.0);
        let n = rng.log_uniform(0.0, 26.0);
        act(&mut m, l, n);
    }

    // ── atoms_after_half_lives ──
    let aah = |m: &mut Module, n0: f64, n: f64| {
        let out = atoms_after_half_lives(bb(n0), bb(n));
        let ins = [f(n0), f(n)];
        match half_lives_exact(n0, n) {
            Some((x, how)) if x.to_bits() != out.to_bits() && !(x.is_nan() && out.is_nan()) => {
                let note = format!("ADR 0002 #13: the port returns {how} = {x:?} (f:{:016x})", x.to_bits());
                m.row_noted("nuclear::atoms_after_half_lives", Class::C, &ins, &[f(out)], &note);
            }
            _ => m.row("nuclear::atoms_after_half_lives", Class::C, &ins, &[f(out)]),
        }
    };
    for x in edges() {
        aah(&mut m, x, 2.0);
        aah(&mut m, 1000.0, x);
    }
    aah(&mut m, 1000.0, 2.0); // the Rust test
    // integer n: 2⁻ⁿ is exact from n = −1023 to 1074 (subnormal from 1023), so N₀·2⁻ⁿ is correctly rounded there; past
    // the ends 2⁻ⁿ is ∞ or 0 (at n = 1075 a tie, rounded to even: 0)
    let mut ints: Vec<i32> = (-1100..=-1000).chain(-40..=40).chain(1000..=1100).collect();
    ints.extend((-990..=990).step_by(30));
    for k in ints {
        aah(&mut m, 1000.0, f64::from(k));
    }
    // fractional n, and 2⁻ⁿ over / underflowing where N₀ would bring N₀/2ⁿ back into range
    for k in [-1100.0, -1024.0, -1023.5, -100.0, -1.0, 0.0, 1.0, 10.0, 100.0, 1022.0, 1074.0, 1075.0] {
        for frac in [0.25, 0.5, 0.75, 1.0 / 3.0] {
            aah(&mut m, 1000.0, k + frac);
        }
    }
    for (n0, n) in [
        (pow2(100), 1100.0),
        (pow2(-100), -1100.0),
        (f64::MAX, 1080.0),
        (f64::from_bits(1), -1080.0),
        (pow2(52), 1126.0),
        (3.0, 1074.0),
        (3.0, 1075.0),
        (3.0, 1074.5),
        (0.0, -1100.0), // 0·∞ = NaN
        (f64::INFINITY, 1100.0), // ∞·0 = NaN
        (-1.0, -1100.0),
    ] {
        aah(&mut m, n0, n);
    }
    // ADR 0006 item 6: 2⁻ⁿ subnormal, N₀ scaling it back into the normal range
    for k in 0..=214 {
        let n = 1022.0 + f64::from(k) * 0.25;
        for p in subnormal_prefactors() {
            aah(&mut m, p, n);
        }
    }
    for _ in 0..300 {
        let n = rng.uniform(1022.0, 1075.0);
        let n0 = rng.log_uniform(0.0, 308.0);
        aah(&mut m, n0, n);
    }
    for _ in 0..300 {
        let n0 = rng.log_uniform(0.0, 25.0);
        let n = rng.uniform(0.0, 60.0);
        aah(&mut m, n0, n);
    }
    for _ in 0..300 {
        let n0 = rng.signed_log_uniform(-300.0, 300.0);
        let n = rng.uniform(-1100.0, 1100.0);
        aah(&mut m, n0, n);
    }

    // ADR 0002 #13 (M4c follow-up, a fresh seed so no earlier row moves): N₀ = ±2^j, ±0, ±∞, NaN and negative N₀ at the
    // edges; exact subnormal ties for an integer n; integer and fractional n past both ends with N₀ bringing N₀/2ⁿ into
    // range
    let mut rng13 = Rng::new(0x4c13);
    for (n0, n) in [
        (pow2(100), 1100.5),
        (pow2(-1074), -1100.25),
        (pow2(1023), 2097.5),
        (pow2(-1000), -2000.75),
        (-pow2(60), 1100.125),
        (0.0, -1100.0),
        (-0.0, -1100.0),
        (0.0, 1100.5),
        (-0.0, 1100.5),
        (f64::INFINITY, 1100.0),
        (f64::NEG_INFINITY, 1100.0),
        (f64::NEG_INFINITY, -1100.0),
        (f64::NAN, 1100.0),
        (-1000.0, 1075.0),
        (-3.0, 1075.0),
        (5.0, 1075.0),
        (7.0, 1075.0),
        (1.0, 1075.0),
        (1.0, 1074.5),
        (1.0, 1023.5),
        (3.0, 1076.0),
        (f64::MAX, 2098.0),
        (f64::MAX, 2099.0),
        (f64::from_bits(1), -2097.0),
        (f64::from_bits(1), -2098.0),
        (f64::MIN_POSITIVE, -2046.0),
        (f64::MIN_POSITIVE, -2045.5),
        (1.5, 2300.0),
        (1.5, -2300.0),
    ] {
        aah(&mut m, n0, n);
    }
    for _ in 0..300 {
        let n = f64::from(i32::try_from(rng13.int(1075, 1200)).unwrap_or(1075));
        let n0 = rng13.log_uniform(0.0, 308.0);
        aah(&mut m, n0, n);
    }
    for _ in 0..300 {
        let n = f64::from(i32::try_from(rng13.int(-1200, -1024)).unwrap_or(-1024));
        let n0 = rng13.log_uniform(-320.0, 0.0);
        aah(&mut m, n0, n);
    }
    for _ in 0..300 {
        let n = if rng13.int(0, 1) == 0 { rng13.uniform(1022.0, 1200.0) } else { rng13.uniform(-1200.0, -1024.0) };
        let n0 = if n > 0.0 { rng13.log_uniform(-10.0, 308.0) } else { rng13.log_uniform(-320.0, 10.0) };
        aah(&mut m, n0, n);
    }

    // ── time_for_fraction ──
    let tff = |m: &mut Module, fr: f64, l: f64| {
        let out = res(time_for_fraction(bb(fr), bb(l)));
        m.row("nuclear::time_for_fraction", Class::C, &[f(fr), f(l)], &[out]);
    };
    for x in edges() {
        tff(&mut m, x, 0.01);
        tff(&mut m, 0.5, x);
    }
    // the Rust tests' values (and the order of the checks: the fraction, then λ)
    for (fr, l) in [
        (0.5, l100),
        (0.0, 0.01),
        (1.5, 0.01),
        (0.5, 0.0),
        (f64::NAN, 0.01),
        (0.5, f64::NAN),
        (0.5, f64::INFINITY),
        (f64::NAN, f64::NAN),
        (1.5, -1.0),
        (1.0, 0.01), // −ln 1 = −0.0
    ] {
        tff(&mut m, fr, l);
    }
    let mut fracs = vec![1.0, 1.0_f64.next_up(), 0.5, 0.25, 0.1, 0.01, 1e-300, f64::MIN_POSITIVE, f64::from_bits(1)];
    for k in 1..=4 {
        fracs.push(ulps(1.0, -k));
    }
    fracs.extend(decades(-323, -1, &["1", "5"]).into_iter().step_by(4));
    let mut ls: Vec<f64> = decades(-300, 300, &["1"]).into_iter().step_by(25).collect();
    ls.extend([f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1), l100]);
    for &fr in &fracs {
        for &l in &ls {
            tff(&mut m, fr, l);
        }
    }
    for &l in &lambdas {
        for fr in [0.5, 0.1, 0.01, 0.999] {
            tff(&mut m, fr, l);
        }
    }
    for _ in 0..300 {
        let fr = rng.unit();
        let l = rng.log_uniform(-20.0, 5.0);
        tff(&mut m, fr, l);
    }
    for _ in 0..200 {
        let fr = 1.0 - rng.log_uniform(-16.0, -1.0);
        let l = rng.log_uniform(-300.0, 300.0);
        tff(&mut m, fr, l);
    }

    // ── mass_defect, binding_energy_mev, binding_energy_per_nucleon ──
    let md = |m: &mut Module, z: u32, n: u32, mass: f64| {
        let out = mass_defect(bb(z), bb(n), bb(mass));
        m.row("nuclear::mass_defect", Class::B, &[i(z), i(n), f(mass)], &[f(out)]);
    };
    for z in U32S {
        for n in U32S {
            md(&mut m, z, n, 4.002_603_254);
        }
    }
    for x in edges() {
        md(&mut m, 2, 2, x);
        md(&mut m, 0, 0, x);
    }
    md(&mut m, 2, 2, 4.002_603_254); // the Rust tests' He-4 (AME2020)
    md(&mut m, 26, 30, 55.934_936_33); // Fe-56 (AME2020)
    for _ in 0..300 {
        let z = u32::try_from(rng.int(1, 118)).unwrap_or(1);
        let n = u32::try_from(rng.int(0, 180)).unwrap_or(0);
        let mass = f64::from(z + n) * rng.uniform(0.99, 1.01);
        md(&mut m, z, n, mass);
    }
    let bem = |m: &mut Module, x: f64| {
        let out = binding_energy_mev(bb(x));
        m.row("nuclear::binding_energy_mev", Class::B, &[f(x)], &[f(out)]);
    };
    let mut xs = edges();
    xs.push(mass_defect(2, 2, 4.002_603_254));
    xs.extend(decades(-323, 308, &["1", "3.0377"]).into_iter().step_by(3));
    for _ in 0..200 {
        xs.push(rng.signed_log_uniform(-6.0, 1.0));
    }
    for x in xs {
        bem(&mut m, x);
    }
    let bea = audit::bea;
    for z in U32S {
        for n in U32S {
            bea(&mut m, z, n, 55.934_936_33);
        }
    }
    for x in edges() {
        bea(&mut m, 26, 30, x);
        bea(&mut m, 0, 0, x); // A = 0 is checked before the mass
        bea(&mut m, u32::MAX, 1, x); // the overflow first
    }
    for (z, n, mass) in [
        (26, 30, 55.934_936_33), // the Rust tests: Fe-56 (AME2020)
        (u32::MAX, 1, 1.0),
        (26, 30, f64::NAN),
        (26, 30, f64::INFINITY),
        (2, 2, 4.002_603_254),
        (u32::MAX, 0, 1.0),
        (1 << 31, 1 << 31, 1.0), // 2³² overflows
        ((1 << 31) - 1, 1 << 31, 1.0), // u32::MAX: no overflow
    ] {
        bea(&mut m, z, n, mass);
    }
    for _ in 0..300 {
        let z = u32::try_from(rng.int(1, 118)).unwrap_or(1);
        let n = u32::try_from(rng.int(0, 180)).unwrap_or(0);
        let mass = f64::from(z + n) * rng.uniform(0.99, 1.01);
        bea(&mut m, z, n, mass);
    }

    // ── semi_empirical_binding_energy (ADR 0002 #3: the oracle's NaN is the port's Err) ──
    let semf = |m: &mut Module, z: u32, a: u32| {
        let be = semi_empirical_binding_energy(bb(z), bb(a));
        if be.is_nan() {
            m.row_noted(
                "nuclear::semi_empirical_binding_energy",
                Class::B,
                &[i(z), i(a)],
                &[f(be)],
                "ADR 0002 #3: the port returns Err(InvalidInput)",
            );
        } else {
            m.row("nuclear::semi_empirical_binding_energy", Class::B, &[i(z), i(a)], &[f(be)]);
        }
    };
    for z in 1..=300u32 {
        for a in z..=(3 * z).min(300) {
            semf(&mut m, z, a);
        }
    }
    for a in [1, 2, 3, 4, 5, 8, 10, 50, 100, 300, 65_535, u32::MAX] {
        semf(&mut m, 0, a);
    }
    for (z, a) in [
        (26, 56), // the Rust tests: Fe-56, U-238, He-4
        (92, 238),
        (2, 4),
        (82, 208),
        (1, 1000),
        (92, 1000),
        (65_535, 65_535),
        (1 << 16, 1 << 17),
        (100_000, 300_000),
        (u32::MAX, u32::MAX),
        (1, u32::MAX),
        (u32::MAX / 3, u32::MAX),
        (u32::MAX / 2, u32::MAX),
        (1 << 31, u32::MAX),
        (1 << 31, 1 << 31),
        ((1 << 31) - 1, u32::MAX - 1),
    ] {
        semf(&mut m, z, a);
    }
    // the Err region: z > a (Rust's checked_sub) and a = 0
    let mut errs: Vec<(u32, u32)> = vec![(3, 2), (0, 0), (1, 0), (2, 0), (300, 0), (u32::MAX, 0), (10, 5), (300, 1)];
    for a in [1, 2, 3, 4, 5, 10, 56, 237, 299, 300, 65_535, u32::MAX - 1] {
        errs.push((a + 1, a));
    }
    errs.extend([(u32::MAX, 1), (u32::MAX, 1 << 31), (1 << 31, (1 << 31) - 1), (1000, 999)]);
    let mut seen = std::collections::BTreeSet::new();
    for (z, a) in errs {
        if seen.insert((z, a)) {
            semf(&mut m, z, a);
        }
    }

    // ── q_value ──
    let qv = |m: &mut Module, r: &[f64], p: &[f64]| {
        let mut ins = vec![Tok::V(r.len())];
        ins.extend(r.iter().map(|&x| f(x)));
        ins.push(Tok::V(p.len()));
        ins.extend(p.iter().map(|&x| f(x)));
        audit::qv(m, &ins, r, p);
    };
    qv(&mut m, &[2.014_10, 3.016_05], &[4.002_60, 1.008_66]); // the Rust tests: D + T → He-4 + n, and back
    qv(&mut m, &[4.002_60, 1.008_66], &[2.014_10, 3.016_05]);
    qv(&mut m, &[], &[]);
    for x in edges() {
        qv(&mut m, &[x], &[]);
        qv(&mut m, &[], &[x]);
        qv(&mut m, &[x], &[x]);
        qv(&mut m, &[x, 1.0], &[1.0]);
    }
    for (r, p) in [
        (vec![-0.0], vec![-0.0]),
        (vec![-0.0, -0.0], vec![]),
        (vec![], vec![-0.0]),
        (vec![0.0], vec![-0.0]),
        (vec![f64::INFINITY], vec![f64::INFINITY]),
        (vec![f64::INFINITY, f64::NEG_INFINITY], vec![]),
        (vec![f64::MAX, f64::MAX], vec![]),
        (vec![f64::MAX], vec![-f64::MAX]),
        (vec![1e16, 1.0, -1e16], vec![]),
        (vec![1.0, 1e16, -1e16], vec![]),
        (vec![235.043_928, 1.008_665], vec![140.914_41, 91.926_156, 1.008_665, 1.008_665, 1.008_665]),
    ] {
        qv(&mut m, &r, &p);
    }
    for _ in 0..300 {
        let nr = usize::try_from(rng.int(0, 6)).unwrap_or(0);
        let np = usize::try_from(rng.int(0, 6)).unwrap_or(0);
        let r: Vec<f64> = (0..nr).map(|_| rng.log_uniform(0.0, 2.5)).collect();
        let p: Vec<f64> = (0..np).map(|_| rng.log_uniform(0.0, 2.5)).collect();
        qv(&mut m, &r, &p);
    }

    // ── decay_chain_two_step ──
    let chain = |m: &mut Module, n0: f64, l1: f64, l2: f64, t: f64| {
        let out = triple(decay_chain_two_step(bb(n0), bb(l1), bb(l2), bb(t)));
        m.row("nuclear::decay_chain_two_step", Class::C, &[f(n0), f(l1), f(l2), f(t)], &out);
    };
    for x in edges() {
        chain(&mut m, x, 0.1, 0.05, 10.0);
        chain(&mut m, 1000.0, x, 0.05, 10.0);
        chain(&mut m, 1000.0, 0.1, x, 10.0);
        chain(&mut m, 1000.0, 0.1, 0.05, x);
    }
    // the Rust tests' values (and the order of the checks: λ₁ and λ₂, then N₀, then t)
    for (n0, l1, l2, t) in [
        (1000.0, 0.1, 0.05, 10.0),
        (1000.0, 0.1, 0.05, 0.0),
        (1000.0, 0.1, 0.05, 200.0),
        (1000.0, 0.1, 0.1, 10.0),
        (1000.0, 0.1, 0.1 * (1.0 + 1e-10), 10.0),
        (1.0, 2.0, 0.5, 3.0),
        (1.0, 10.0, 0.1, 100.0),
        (1000.0, 0.1, 0.05, -1.0),
        (-1.0, 0.1, 0.05, 1.0),
        (f64::NAN, 0.1, 0.05, 1.0),
        (1000.0, f64::NAN, 0.05, 1.0),
        (1000.0, 0.1, f64::INFINITY, 1.0),
        (1000.0, 0.1, 0.05, f64::NAN),
        (f64::NAN, f64::NAN, 0.05, f64::NAN),
        (f64::NAN, 0.1, f64::NAN, f64::NAN),
        (-1.0, 0.1, 0.05, f64::NAN),
        (1000.0, 0.1, 0.05, -0.0),
        (-0.0, 0.1, 0.05, 1.0),
        (-0.0, 0.1, 0.1, -0.0),
    ] {
        chain(&mut m, n0, l1, l2, t);
    }
    let lgrid = [
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1e-300,
        1e-20,
        1e-10,
        1e-3,
        0.05,
        0.1,
        0.5,
        1.0,
        2.0,
        10.0,
        1e3,
        1e10,
        1e300,
        f64::MAX,
    ];
    let tgrid = [0.0, -0.0, 1e-300, 1e-10, 0.1, 1.0, 10.0, 100.0, 1e3, 1e10, 1e300, f64::MAX];
    for l1 in lgrid {
        for l2 in lgrid {
            for t in tgrid {
                chain(&mut m, 1000.0, l1, l2, t);
            }
        }
    }
    for n0 in [0.0, -0.0, 1.0, 1e300, f64::MAX, f64::from_bits(1)] {
        for (l1, l2) in [(0.1, 0.05), (0.05, 0.1), (0.1, 0.1), (1e10, 1.0), (1e-300, 1e-300), (f64::MAX, 1.0)] {
            for t in [0.0, 1.0, 10.0, 1e3, 1e300] {
                chain(&mut m, n0, l1, l2, t);
            }
        }
    }
    // λ₁ ≈ λ₂: relative differences either side of Rust's 1e-12 threshold, both orders
    for l1 in [1e-10, 0.1, 1.0, 1e10] {
        for r in [1e-15, 1e-13, 5e-13, 9e-13, 1e-12, 1.1e-12, 2e-12, 1e-11, 1e-10, 1e-6] {
            for sign in [1.0, -1.0] {
                let l2 = l1 * (1.0 + sign * r);
                for t in [1.0, 10.0, 1e3, 1e9] {
                    chain(&mut m, 1000.0, l1, l2, t / l1);
                }
            }
        }
    }
    // on the threshold itself (to the ulp of 1e-12)
    for (l1, l2) in threshold_pairs() {
        for t in [10.0, 1e6] {
            chain(&mut m, 1000.0, l1, l2, t);
        }
        let lmax = l1.max(l2);
        chain(&mut m, 1000.0, l1, l2, 10.0 / lmax);
    }
    // ADR 0006 item 6: e^(−λ₁t) (N_A) and e^(−λ_min·t) (N_B) subnormal, N₀ and N₀·λ₁ scaling them back
    for x in EXP_SUBNORMAL_X {
        for p in subnormal_prefactors() {
            chain(&mut m, p, 1.0, 2.0, -x); // N_A = N₀·e^x
            chain(&mut m, p, 1e10, 1.0, -x); // N_B = N₀·λ₁·e^x·growth, growth ≈ 1/(λ₁ − λ₂)
            chain(&mut m, p, 1.0, 1.0, -x); // the equal-rate limit: N_B = N₀·λ·e^x·t
        }
    }
    for _ in 0..300 {
        let x = rng.uniform(-745.2, -708.0);
        let lmin = rng.log_uniform(-3.0, 3.0);
        let lmax = lmin * rng.uniform(1.0, 1e3);
        let n0 = rng.log_uniform(0.0, 308.0);
        if rng.int(0, 1) == 0 {
            chain(&mut m, n0, lmax, lmin, -x / lmin);
        } else {
            chain(&mut m, n0, lmin, lmax, -x / lmin);
        }
    }
    // every isotope pair as parent → daughter, at the parent's half-life and a day
    for (rp, &lp) in ISOTOPES.iter().zip(&lambdas) {
        for &ld in &lambdas {
            for t in [rp.half_life_s, 86_400.0] {
                chain(&mut m, 6.022e23, lp, ld, t);
            }
        }
    }
    for _ in 0..400 {
        let n0 = rng.log_uniform(0.0, 24.0);
        let l1 = rng.log_uniform(-20.0, 0.0);
        let l2 = rng.log_uniform(-20.0, 0.0);
        let t = rng.log_uniform(0.0, 12.0);
        chain(&mut m, n0, l1, l2, t);
    }
    for _ in 0..200 {
        let n0 = rng.log_uniform(-300.0, 300.0);
        let l1 = rng.log_uniform(-300.0, 300.0);
        let l2 = rng.log_uniform(-300.0, 300.0);
        let t = rng.log_uniform(-300.0, 300.0);
        chain(&mut m, n0, l1, l2, t);
    }

    // ── is_secular_equilibrium ──
    let sec = |m: &mut Module, lp: f64, ld: f64| {
        let out = is_secular_equilibrium(bb(lp), bb(ld));
        let ex = secular_exact(lp, ld);
        if ex == out {
            m.row("nuclear::is_secular_equilibrium", Class::A, &[f(lp), f(ld)], &[b(out)]);
        } else {
            let note = format!("ADR 0002 #14: the port returns the exact λ_p < λ_d/100 = {ex} (i:{})", i64::from(ex));
            m.row_noted("nuclear::is_secular_equilibrium", Class::A, &[f(lp), f(ld)], &[b(out)], &note);
        }
    };
    for x in edges() {
        for y in edges() {
            sec(&mut m, x, y);
        }
    }
    sec(&mut m, 5e-18, 3.3e-7); // the Rust tests: U-238 → Th-234, and two similar rates
    sec(&mut m, 0.1, 0.2);
    for ld in [1.0, 3.0, 7.0, 0.1, 100.0, 1e-300, 1e300, f64::MAX, f64::from_bits(1), 3.3e-7, -1.0] {
        let at = 0.01 * ld;
        for d in -2..=2 {
            sec(&mut m, ulps(at, d), ld);
        }
    }
    for _ in 0..200 {
        let lp = rng.log_uniform(-30.0, 0.0);
        let ld = rng.log_uniform(-30.0, 0.0);
        sec(&mut m, lp, ld);
        let at = 0.01 * ld;
        sec(&mut m, at, ld);
    }

    // ADR 0002 #14 (a fresh seed): λ_p at the rounded λ_d/100 and one ulp either side over every magnitude and sign, the
    // underflow of 0.01·λ_d, the overflow end, signed zeros
    let mut rng14 = Rng::new(0x4c14);
    for _ in 0..300 {
        let ld = rng14.signed_log_uniform(-300.0, 300.0);
        let at = ld / 100.0;
        for d in -1..=1 {
            sec(&mut m, ulps(at, d), ld);
        }
    }
    for (lp, ld) in [
        (0.03, 3.0),
        (-0.02, -2.0),
        (-0.03, -3.0),
        (0.0, f64::from_bits(1)),
        (-0.0, f64::from_bits(1)),
        (f64::from_bits(1), f64::from_bits(100)),
        (f64::from_bits(1), f64::from_bits(101)),
        (-f64::from_bits(1), -f64::from_bits(99)),
        (f64::from_bits(1), f64::MIN_POSITIVE),
        (f64::MAX / 100.0, f64::MAX),
        ((f64::MAX / 100.0).next_up(), f64::MAX),
        (f64::MAX, f64::MAX),
        (1e306, f64::INFINITY),
        (f64::NEG_INFINITY, -1e300),
        (f64::INFINITY, f64::INFINITY),
        (-0.0, -0.0),
        (-1e-300, 0.0),
    ] {
        sec(&mut m, lp, ld);
    }

    // ── isotope data ──
    m.row("nuclear::ISOTOPES", Class::A, &[], &[Tok::V(ISOTOPES.len())]);
    for (k, r) in ISOTOPES.iter().enumerate() {
        m.row("nuclear::ISOTOPES[]", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &isotope_toks(r));
    }
    let lookup = |m: &mut Module, sym: &str| {
        let out = lookup_isotope(bb(sym)).map_or_else(|| vec![Tok::None], isotope_toks);
        m.row("nuclear::lookup_isotope", Class::A, &[s(sym)], &out);
    };
    let syms: Vec<&str> = ISOTOPES.iter().map(|r| r.symbol).collect();
    for sym in &syms {
        lookup(&mut m, sym);
    }
    let mut others: Vec<String> = Vec::new();
    for sym in &syms {
        others.push(sym.to_lowercase());
        others.push(sym.to_uppercase());
        others.push(format!(" {sym}"));
        others.push(format!("{sym} "));
        others.push(sym.replace('-', ""));
        for n in 1..sym.len() {
            others.push(sym[..n].to_string());
        }
    }
    others.extend(
        ["", "Xx-999", "C-12", "Tc-99", "Tc-99M", "U-235m", "U-2350", "Carbon-14", "14C", "C 14", "H-2", "Pu-240"]
            .iter()
            .map(|x| (*x).to_string()),
    );
    others.sort();
    others.dedup();
    others.retain(|x| !syms.contains(&x.as_str()));
    for sym in &others {
        lookup(&mut m, sym);
    }
    audit::emit_nuclear(&mut m);
    audit::emit_nuclear_b4(&mut m);
    m
}
