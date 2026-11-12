//! `inorganic.tsv` — crystal field theory, the spectrochemical series, spin-only moments, lattice energies and ionic
//! radii (port plan §4, M4b).
//!
//! Rust's `u8` (d electrons, unpaired electrons), `u32` (ion charges) and `i8` (an ionic radius's charge) are `i:`
//! tokens; a `bool` is `i:0` / `i:1`.
//! - Class A: `Ligand::field_strength` and `Ligand::is_strong_field` (all 11 variants, input `e:Ligand.<Variant>`);
//!   `inorganic::cfse_octahedral` (every u8 × both spins: Ok for 0 ..= 10, the Err beyond; the oracle's value. **ADR 0002
//!   #12** makes the port's exact, (3·e_g − 2·t2g)/5 in one correctly rounded division, so the rows compare bit for bit:
//!   the 13 where the oracle's `t2g·(−0.4) + eg·0.6` misses the exact value carry it in their note, and the harness
//!   maps them, accepted under ADR 2, by name);
//!   `inorganic::unpaired_electrons_oct` (every u8 × both spins: Ok for 0 ..= 10, the Err beyond);
//!   `inorganic::IONIC_RADII` (`v:30`) and one `IONIC_RADII[]` row per entry (input `i:k`; output
//!   `s:symbol i:charge i:coordination f:radius_pm`); `inorganic::lookup_ionic_radius` (every table symbol × every i8
//!   charge; lower- and upper-case forms, a leading / trailing space, every proper prefix, `""`, `"Xx"`, `"Fe2+"`,
//!   element names × the charges −4 … 4; the four fields or `n:`).
//! - Class B (+ − × ÷ and √, expected 0 ulp): `inorganic::cfse_octahedral` (every u8 × both spins; the occupancy
//!   tables times −0.4 and 0.6, the oracle's value; **ADR 0002 #12**: where it misses the exact (3·e_g − 2·t2g)/5, on
//!   13 of the 22 Ok rows, the row is noted with the exact value, which the port returns and the harness maps), `delta_tet_from_oct` (the edges, every decade from the smallest subnormal to DBL_MAX,
//!   DBL_MAX / 4 ± 1 ulp where 4x overflows, seeded), `spin_only_moment` (every u8), `born_lande_lattice_energy`
//!   (each argument over the edges with the others NaCl's, both charges over 0, 1, 2, 3, 4, 7, 255, 256, 65 535, 2³¹,
//!   u32::MAX, the Born exponent at 1 and the ulps above, every Madelung constant × every Born exponent and every mean
//!   of two × charges 1–3 × r₀ 1.5–4 Å, r₀ and M over every tenth decade (overflow to −∞, underflow to −0.0), seeded
//!   salts; every Err path), `born_haber_lattice_energy` (each argument over the edges, the Rust tests' NaCl cycle,
//!   ∞ − ∞, seeded cycles).
//! - ADR 0002 #32 (the overflow audit's B3): `born_lande_lattice_energy` rows are noted where the port's Rust sequence
//!   with an unbounded exponent (the final product rounded once) differs from Rust's value (`audit.rs`), and
//!   `audit::emit_inorganic` appends its edge rows (the prefix past DBL_MAX or below 2⁻¹⁰²², brought back by e² and r₀)
//!   with its own seed.
//! - ADR 0002 #33 (the overflow audit's B4): `born_haber_lattice_energy` rows are noted where Rust's differences with
//!   an unbounded exponent (the final subtraction rounded once) differ from Rust's value (`audit.rs`), and
//!   `audit::emit_inorganic_b4` appends its edge rows (partial differences past ±DBL_MAX that later terms bring back)
//!   with their own seed.
//!
//! Every input goes through `black_box`.

use std::hint::black_box as bb;

use crate::audit;
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, b, e, err, f, i, s};
use kimiya::inorganic::{
    BORN_EXP_AR, BORN_EXP_HE, BORN_EXP_KR, BORN_EXP_NE, BORN_EXP_XE, IONIC_RADII, IonicRadius, Ligand, MADELUNG_CSCL,
    MADELUNG_FLUORITE, MADELUNG_NACL, MADELUNG_RUTILE, MADELUNG_ZNS,
    cfse_octahedral, delta_tet_from_oct, lookup_ionic_radius, spin_only_moment,
    unpaired_electrons_oct,
};

/// Every `Ligand`, in Rust declaration order. `t_inorganic` proves the list against the enum in
/// `rust-old/src/inorganic.rs` (the enum is `#[non_exhaustive]`, so no match here could catch a missing variant).
pub const LIGANDS: [Ligand; 11] = [
    Ligand::Iodide,
    Ligand::Bromide,
    Ligand::Chloride,
    Ligand::Fluoride,
    Ligand::Hydroxide,
    Ligand::Water,
    Ligand::Ammonia,
    Ligand::Ethylenediamine,
    Ligand::Nitrite,
    Ligand::Cyanide,
    Ligand::CarbonMonoxide,
];

/// An ionic radius's four fields, in declaration order.
pub fn radius_toks(r: &IonicRadius) -> Vec<Tok> {
    vec![s(r.symbol), i(r.charge), i(r.coordination), f(r.radius_pm)]
}

/// The charge probes of `born_lande_lattice_energy` (Rust's u32: 0 is the Err case).
const CHARGES: [u32; 11] = [0, 1, 2, 3, 4, 7, 255, 256, 65_535, 1 << 31, u32::MAX];

/// The Madelung constants and the Born exponents.
const MADELUNG: [f64; 5] = [MADELUNG_NACL, MADELUNG_CSCL, MADELUNG_ZNS, MADELUNG_FLUORITE, MADELUNG_RUTILE];
const BORN: [f64; 5] = [BORN_EXP_HE, BORN_EXP_NE, BORN_EXP_AR, BORN_EXP_KR, BORN_EXP_XE];

/// ADR 0002 #12: the exact CFSE of an oracle value v ≈ q/5 (Rust's `t2g * (-0.4) + eg * 0.6` is within a few ulps of
/// (3·e_g − 2·t2g)/5): q = round(5v), checked to be that close, then q / 5 in one correctly rounded division, the
/// port's value. `+ 0.0` makes a zero numerator +0 (the port converts the integer 0). Its rows keep the oracle's
/// value; the note names the exact one, and the harness maps the two.
fn exact_cfse(v: f64) -> f64 {
    let q = (v * 5.0).round() + 0.0;
    assert!((v * 5.0 - q).abs() < 1e-9 && q.abs() <= 12.0, "cfse {v:?} is no fifth");
    q / 5.0
}

fn lande(m: &mut Module, madelung: f64, zc: u32, za: u32, r0: f64, n: f64) {
    audit::lande(m, madelung, zc, za, r0, n);
}

fn haber(m: &mut Module, x: [f64; 5]) {
    audit::haber(m, x);
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("inorganic");
    let mut rng = Rng::new(0x4b02);

    // ── crystal field theory ──
    for d in 0..=u8::MAX {
        for hs in [false, true] {
            let r = cfse_octahedral(bb(d), bb(hs));
            let ins = [i(d), b(hs)];
            match r {
                Ok(v) if exact_cfse(v).to_bits() != v.to_bits() => {
                    let x = exact_cfse(v);
                    let note = format!(
                        "ADR 0002 #12: the port returns the exact (3·e_g − 2·t2g)/5 = {x:?} (f:{:016x})",
                        x.to_bits()
                    );
                    m.row_noted("inorganic::cfse_octahedral", Class::A, &ins, &[f(v)], &note);
                }
                _ => m.row("inorganic::cfse_octahedral", Class::A, &ins, &[res(r)]),
            }
        }
    }
    let mut tet = edges();
    tet.extend(decades(-323, 308, &["1", "4.5"]));
    let quarter = f64::MAX / 4.0;
    tet.extend([quarter, quarter.next_up(), quarter.next_down(), -quarter.next_up(), 10_000.0, 1.0 / 9.0, 9.0]);
    for _ in 0..300 {
        tet.push(rng.signed_log_uniform(-300.0, 300.0));
    }
    // ADR 0002 #15 (M4c follow-up; a fresh seed, so no earlier row moves): Δ ≥ 2¹⁰²², where Rust's Δ·4 overflows
    let mut rng15 = Rng::new(0x4b15);
    tet.extend([pow2(1022), -pow2(1022), pow2(1022).next_down(), f64::MAX.next_down(), 1.5 * pow2(1022)]);
    for _ in 0..50 {
        tet.push(pow2(1022) * rng15.uniform(1.0, 4.0));
    }
    for x in tet {
        let out = delta_tet_from_oct(bb(x));
        // ADR 0002 #15: where Δ·4 overflows, the port's (Δ/9)·4, the correctly rounded 4Δ/9 (validate.py proves it)
        let exact = if out.is_infinite() && x.is_finite() { (x / 9.0) * 4.0 } else { out };
        if exact.to_bits() == out.to_bits() || (exact.is_nan() && out.is_nan()) {
            m.row("inorganic::delta_tet_from_oct", Class::B, &[f(x)], &[f(out)]);
        } else {
            let note = format!(
                "ADR 0002 #15: the port returns the correctly rounded 4Δ/9 = {exact:?} (f:{:016x})",
                exact.to_bits()
            );
            m.row_noted("inorganic::delta_tet_from_oct", Class::B, &[f(x)], &[f(out)], &note);
        }
    }

    // ── spectrochemical series ──
    for l in LIGANDS {
        m.row("Ligand::field_strength", Class::A, &[e("Ligand", l)], &[i(bb(l).field_strength())]);
    }
    for l in LIGANDS {
        m.row("Ligand::is_strong_field", Class::A, &[e("Ligand", l)], &[b(bb(l).is_strong_field())]);
    }

    // ── magnetic moment and unpaired electrons ──
    for n in 0..=u8::MAX {
        let mu = spin_only_moment(bb(n));
        m.row("inorganic::spin_only_moment", Class::B, &[i(n)], &[f(mu)]);
    }
    for d in 0..=u8::MAX {
        for hs in [false, true] {
            let out = unpaired_electrons_oct(bb(d), bb(hs)).map_or_else(|e| err(&e), i);
            m.row("inorganic::unpaired_electrons_oct", Class::A, &[i(d), b(hs)], &[out]);
        }
    }

    // ── Born–Landé ──
    let (nacl_r0, nacl_n) = (2.82e-10, 8.0);
    for x in edges() {
        lande(&mut m, x, 1, 1, nacl_r0, nacl_n);
        lande(&mut m, MADELUNG_NACL, 1, 1, x, nacl_n);
        lande(&mut m, MADELUNG_NACL, 1, 1, nacl_r0, x);
        lande(&mut m, x, 1, 1, x, x); // r₀ is checked first, then n, then M
        lande(&mut m, x, 0, 0, nacl_r0, nacl_n); // M before the charges
    }
    for zc in CHARGES {
        for za in CHARGES {
            lande(&mut m, MADELUNG_NACL, zc, za, nacl_r0, nacl_n);
        }
        lande(&mut m, f64::NAN, zc, 1, nacl_r0, nacl_n);
    }
    let one = 1.0_f64;
    let mut ns = vec![one, one.next_up(), one.next_up().next_up(), one.next_down(), 1.0 + 1e-15, 1.5, 0.5, -1.0, 1e300];
    ns.extend(BORN);
    for &n in &ns {
        lande(&mut m, MADELUNG_NACL, 1, 1, nacl_r0, n);
    }
    let mut exps: Vec<f64> = BORN.to_vec();
    for (k, &a) in BORN.iter().enumerate() {
        for &c in &BORN[k + 1..] {
            exps.push((a + c) / 2.0); // the mean of two configurations (NaCl: (7 + 9) / 2 = 8), exact
        }
    }
    for madelung in MADELUNG {
        for &n in &exps {
            for z in 1..=3u32 {
                for r0_pm in [150.0, 200.0, 236.0, 282.0, 314.0, 356.0, 400.0] {
                    lande(&mut m, madelung, z, z, r0_pm * 1e-12, n);
                }
                lande(&mut m, madelung, z, 1, nacl_r0, n);
                lande(&mut m, madelung, 1, z, nacl_r0, n);
            }
        }
    }
    for x in decades(-320, 308, &["1", "2.82"]).into_iter().step_by(10) {
        lande(&mut m, MADELUNG_NACL, 1, 1, x, nacl_n);
        lande(&mut m, x, 1, 1, nacl_r0, nacl_n);
        lande(&mut m, x, u32::MAX, u32::MAX, 1e-300, nacl_n);
        lande(&mut m, x, 1, 1, f64::MAX, 1e300);
    }
    lande(&mut m, 1.748, 1, 1, 2.82e-10, 8.0);
    lande(&mut m, MADELUNG_NACL, 2, 2, 2.82e-10, 8.0);
    lande(&mut m, -1.748, 1, 1, 2.82e-10, 8.0);
    lande(&mut m, 1.748, 1, 1, -2.82e-10, 8.0);
    for _ in 0..600 {
        let madelung = rng.uniform(1.0, 3.0);
        let zc = u32::try_from(rng.int(1, 4)).unwrap_or(1);
        let za = u32::try_from(rng.int(1, 4)).unwrap_or(1);
        let r0 = rng.log_uniform(-11.0, -9.0);
        let n = rng.uniform(5.0, 12.0);
        lande(&mut m, madelung, zc, za, r0, n);
    }

    // ── Born–Haber ──
    let nacl = [-411.15, 107.3, 495.8, 242.0, -349.0];
    haber(&mut m, nacl);
    for x in edges() {
        for k in 0..5 {
            let mut v = nacl;
            v[k] = x;
            haber(&mut m, v);
        }
    }
    haber(&mut m, [f64::INFINITY, f64::INFINITY, 0.0, 0.0, 0.0]);
    haber(&mut m, [f64::MAX, -f64::MAX, -f64::MAX, -f64::MAX, -f64::MAX]);
    haber(&mut m, [0.0, 0.0, 0.0, 0.0, 0.0]);
    haber(&mut m, [-0.0, 0.0, 0.0, 0.0, 0.0]);
    haber(&mut m, [-0.0, 0.0, 0.0, -0.0, 0.0]);
    for _ in 0..400 {
        let mut v = [0.0; 5];
        for x in &mut v {
            *x = rng.signed_log_uniform(-3.0, 4.0);
        }
        haber(&mut m, v);
    }

    // ── ionic radii ──
    m.row("inorganic::IONIC_RADII", Class::A, &[], &[Tok::V(IONIC_RADII.len())]);
    for (k, r) in IONIC_RADII.iter().enumerate() {
        m.row("inorganic::IONIC_RADII[]", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &radius_toks(r));
    }
    let mut table_syms: Vec<&str> = IONIC_RADII.iter().map(|r| r.symbol).collect();
    table_syms.dedup();
    let mut seen = std::collections::BTreeSet::new();
    table_syms.retain(|x| seen.insert(*x));
    let lookup = |m: &mut Module, sym: &str, q: i8| {
        let out = lookup_ionic_radius(bb(sym), bb(q)).map_or_else(|| vec![Tok::None], radius_toks);
        m.row("inorganic::lookup_ionic_radius", Class::A, &[s(sym), i(q)], &out);
    };
    for sym in &table_syms {
        for q in i8::MIN..=i8::MAX {
            lookup(&mut m, sym, q);
        }
    }
    let mut others: Vec<String> = Vec::new();
    for sym in &table_syms {
        others.push(sym.to_lowercase());
        others.push(sym.to_uppercase());
        others.push(format!(" {sym}"));
        others.push(format!("{sym} "));
        for n in 1..sym.len() {
            others.push(sym[..n].to_string());
        }
    }
    others.extend(
        ["", "Xx", "Fe2+", "Fe3+", "Cl-", "Na+", "Iron", "Sodium", "Chlorine", "S2-", "Sn2"].iter().map(|x| (*x).to_string()),
    );
    others.sort();
    others.dedup();
    others.retain(|x| !table_syms.contains(&x.as_str()));
    for sym in &others {
        for q in -4..=4i8 {
            lookup(&mut m, sym, q);
        }
    }
    audit::emit_inorganic(&mut m);
    audit::emit_inorganic_b4(&mut m);
    m
}
