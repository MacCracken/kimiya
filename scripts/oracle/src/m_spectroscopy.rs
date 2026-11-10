//! `spectroscopy.tsv` — the Beer–Lambert law, photon energy conversions, the Bohr model and the Rydberg formula,
//! vibrational and rotational levels, isotope patterns, the FFT magnitude spectrum and peak detection (port plan §4,
//! M4e).
//!
//! Rust's u8 / u32 arguments are `i:` tokens, `SpectralSeries` an `e:SpectralSeries.<Variant>` token. A slice is `v:n`
//! and its n elements: `isotope_pattern`'s `&[(u8, u32)]` as n `i:z i:count` pairs, `fft_spectrum`'s and `find_peaks`'
//! `&[f64]` as n `f:` tokens (`find_peaks`' threshold follows). A Vec result is `v:n` and its elements: the isotope
//! pattern's (mass, intensity) pairs, the spectrum's floats, the peaks' `i:` indices; `isotope_abundances`' slice the
//! same way.
//! - Class A: `isotope_abundances` (every u8), `SpectralSeries::n_lower` (every variant), `find_peaks` (integers: the
//!   Rust tests' data, fewer than 3 points, thresholds 0, 1, −0.0, the ulps inside [0, 1] and every way out of it, NaN
//!   and ±∞; plateaus, flat and monotone data, NaN and ±∞ in the data, all non-finite, ±0, DBL_MAX ranges where the
//!   interpolated threshold over- or underflows, seeded data).
//! - Class B (+ − × ÷ and √, expected 0 ulp): `absorbance`, `concentration_from_absorbance`, the photon conversions,
//!   `joules_to_ev` / `ev_to_joules`, `reduced_mass` (equal masses included), `harmonic_oscillator_energy`,
//!   `vibrational_frequency`, `rigid_rotor_energy`, `rotational_constant`, `diatomic_moment_of_inertia` (each argument
//!   over the edges, the Rust tests' values, decades, the intermediates' over- and underflow, seeded);
//!   `bohr_energy_level` (every u8 Z × 20 n), `bohr_transition_energy` and `rydberg_wavelength` (Z × (n₁, n₂) grids,
//!   every Err path in Rust's order, adjacent large n), `SpectralSeries::wavelength` (every series × lines to u32::MAX,
//!   Rust's `checked_add` overflow); `isotope_pattern` (the Rust tests' molecules, the empty molecule, every u8 Z alone,
//!   every element with data × counts to 5000 (the pruning keeps the peak count bounded), elements without data with
//!   counts to u32::MAX, H₂… (equal scratch masses: the stable sort's order decides the binned mean), mixtures, the Err
//!   position in a list, seeded molecules); and **ADR 0002 #19**: `wavenumber_from_wavelength` and
//!   `wavelength_from_wavenumber` (the edges, every decade, every power of two, the subnormal-result and overflow ends,
//!   seeded over every binade), the rows where the port's correctly rounded 1/(100·x) differs from Rust's `0.01 / x`
//!   noted with it.
//! - Class C (bound 0 against the corrected oracle): `transmittance` (`10.0.powf(−A)`: A over exp10's whole range, the
//!   integers, subnormal results, seeded), `absorbance_from_transmittance` (log10: every decade of T, T near 1, the
//!   round trips, seeded), `fft_spectrum` (cos and sin twiddles: every length 0 … 16 and every power of two to 2¹³,
//!   2^k ± 1, with impulses, constants, ±0, alternating and Nyquist signals, sinusoids, random, huge and subnormal
//!   values, NaN and ±∞).
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`): the
//! sinusoids are a Taylor series in + − × ÷ (deterministic, not std's sin).

use std::hint::black_box as bb;

use crate::exact::cr_hundredth_div;
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::tsv::{Class, Module, Tok, e, err, f, i};
use kimiya::spectroscopy::{
    SpectralSeries, absorbance, absorbance_from_transmittance, bohr_energy_level, bohr_transition_energy,
    concentration_from_absorbance, diatomic_moment_of_inertia, ev_to_joules, fft_spectrum, find_peaks,
    frequency_from_wavelength, harmonic_oscillator_energy, isotope_abundances, isotope_pattern, joules_to_ev,
    photon_energy_from_frequency, photon_energy_from_wavelength, reduced_mass, rigid_rotor_energy,
    rotational_constant, rydberg_wavelength, transmittance, vibrational_frequency, wavelength_from_energy,
    wavelength_from_wavenumber, wavenumber_from_wavelength,
};

/// Every `SpectralSeries`, in Rust declaration order. `t_spectroscopy` proves the list against the enum in
/// `rust-old/src/spectroscopy.rs` (the enum is `#[non_exhaustive]`, so no match here could catch a missing variant).
pub const SERIES: [SpectralSeries; 5] = [
    SpectralSeries::Lyman,
    SpectralSeries::Balmer,
    SpectralSeries::Paschen,
    SpectralSeries::Brackett,
    SpectralSeries::Pfund,
];

/// The u32 probes of the quantum numbers: the small values, the u8 / u16 / i32 edges, near u32::MAX.
const U32S: [u32; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 100, 1000, 65_535, 1_000_000, 1 << 31, u32::MAX - 2, u32::MAX - 1, u32::MAX,
];

/// The u8 nuclear charges of the transition grids.
const ZS: [u8; 8] = [0, 1, 2, 3, 26, 92, 118, 255];

fn ff(m: &mut Module, name: &str, cls: Class, ins: &[f64], out: f64) {
    let toks: Vec<Tok> = ins.iter().map(|&x| f(x)).collect();
    m.row(name, cls, &toks, &[f(out)]);
}

fn fr(m: &mut Module, name: &str, cls: Class, ins: &[f64], out: Result<f64, kimiya::error::KimiyaError>) {
    let toks: Vec<Tok> = ins.iter().map(|&x| f(x)).collect();
    m.row(name, cls, &toks, &[res(out)]);
}

/// ADR 0002 #19: a wavenumber conversion's row, noted with the port's correctly rounded 1/(100·x) where it differs
/// from Rust's `0.01 / x`.
fn wn(m: &mut Module, name: &str, x: f64, out: Result<f64, kimiya::error::KimiyaError>, var: &str) {
    if let Ok(v) = out {
        let c = cr_hundredth_div(x);
        if c.to_bits() != v.to_bits() {
            let note = format!(
                "ADR 0002 #19: the port returns the correctly rounded 1/(100·{var}) = {c:?} (f:{:016x})",
                c.to_bits()
            );
            m.row_noted(name, Class::B, &[f(x)], &[f(v)], &note);
            return;
        }
    }
    m.row(name, Class::B, &[f(x)], &[res(out)]);
}

fn wavenumbers(m: &mut Module, x: f64) {
    wn(m, "spectroscopy::wavenumber_from_wavelength", x, wavenumber_from_wavelength(bb(x)), "λ");
    wn(m, "spectroscopy::wavelength_from_wavenumber", x, wavelength_from_wavenumber(bb(x)), "ν̃");
}

fn photon(m: &mut Module, x: f64) {
    ff(m, "spectroscopy::photon_energy_from_frequency", Class::B, &[x], photon_energy_from_frequency(bb(x)));
    fr(m, "spectroscopy::photon_energy_from_wavelength", Class::B, &[x], photon_energy_from_wavelength(bb(x)));
    fr(m, "spectroscopy::wavelength_from_energy", Class::B, &[x], wavelength_from_energy(bb(x)));
    fr(m, "spectroscopy::frequency_from_wavelength", Class::B, &[x], frequency_from_wavelength(bb(x)));
    ff(m, "spectroscopy::joules_to_ev", Class::B, &[x], joules_to_ev(bb(x)));
    ff(m, "spectroscopy::ev_to_joules", Class::B, &[x], ev_to_joules(bb(x)));
}

fn bohr(m: &mut Module, z: u8, n: u32) {
    let out = res(bohr_energy_level(bb(z), bb(n)));
    m.row("spectroscopy::bohr_energy_level", Class::B, &[i(z), i(n)], &[out]);
}

fn transition(m: &mut Module, z: u8, lo: u32, hi: u32) {
    let t = res(bohr_transition_energy(bb(z), bb(lo), bb(hi)));
    m.row("spectroscopy::bohr_transition_energy", Class::B, &[i(z), i(lo), i(hi)], &[t]);
    let w = res(rydberg_wavelength(bb(z), bb(lo), bb(hi)));
    m.row("spectroscopy::rydberg_wavelength", Class::B, &[i(z), i(lo), i(hi)], &[w]);
}

fn series_line(m: &mut Module, s: SpectralSeries, line: u32) {
    let out = res(bb(s).wavelength(bb(line)));
    m.row("SpectralSeries::wavelength", Class::B, &[e("SpectralSeries", s), i(line)], &[out]);
}

fn oscillator(m: &mut Module, n: u32, w: f64) {
    let out = harmonic_oscillator_energy(bb(n), bb(w));
    m.row("spectroscopy::harmonic_oscillator_energy", Class::B, &[i(n), f(w)], &[f(out)]);
}

fn rotor(m: &mut Module, j: u32, b: f64) {
    let out = rigid_rotor_energy(bb(j), bb(b));
    m.row("spectroscopy::rigid_rotor_energy", Class::B, &[i(j), f(b)], &[f(out)]);
}

fn pattern(m: &mut Module, atoms: &[(u8, u32)]) {
    let mut ins = vec![Tok::V(atoms.len())];
    for &(z, c) in atoms {
        ins.push(i(z));
        ins.push(i(c));
    }
    let outs = match isotope_pattern(bb(atoms)) {
        Ok(p) => {
            let mut v = vec![Tok::V(p.len())];
            for (mass, int) in p {
                v.push(f(mass));
                v.push(f(int));
            }
            v
        }
        Err(e_) => vec![err(&e_)],
    };
    m.row("spectroscopy::isotope_pattern", Class::B, &ins, &outs);
}

fn spectrum(m: &mut Module, signal: &[f64]) {
    let mut ins = vec![Tok::V(signal.len())];
    ins.extend(signal.iter().map(|&x| f(x)));
    let outs = match fft_spectrum(bb(signal)) {
        Ok(s) => {
            let mut v = vec![Tok::V(s.len())];
            v.extend(s.iter().map(|&x| f(x)));
            v
        }
        Err(e_) => vec![err(&e_)],
    };
    m.row("spectroscopy::fft_spectrum", Class::C, &ins, &outs);
}

fn peaks(m: &mut Module, data: &[f64], threshold: f64) {
    let mut ins = vec![Tok::V(data.len())];
    ins.extend(data.iter().map(|&x| f(x)));
    ins.push(f(threshold));
    let p = find_peaks(bb(data), bb(threshold));
    let mut outs = vec![Tok::V(p.len())];
    for k in p {
        outs.push(i(i64::try_from(k).unwrap_or(-1)));
    }
    m.row("spectroscopy::find_peaks", Class::A, &ins, &outs);
}

/// sin(x) for |x| ≤ π + 1, a Taylor series in + − × ÷ only: a fixed double on every machine (a grid input, not std's
/// sin, which is glibc's).
fn taylor_sin(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for k in 1..30 {
        let kf = f64::from(2 * k);
        term = -term * x2 / (kf * (kf + 1.0));
        sum += term;
    }
    sum
}

/// sin(2π·c·k/n) for the signal sample k of a sinusoid of c cycles over n samples (n a power of two: the fraction
/// (c·k mod n)/n is exact), reduced to [−π, π] before the series.
fn sine(c: u64, k: u64, n: u64) -> f64 {
    let pi = std::f64::consts::PI;
    #[allow(clippy::cast_precision_loss)]
    let frac = ((c * k) % n) as f64 / n as f64;
    let mut x = 2.0 * pi * frac;
    if x > pi {
        x -= 2.0 * pi;
    }
    taylor_sin(x)
}

#[allow(clippy::too_many_lines)] // one function after another, in source order
pub fn emit() -> Module {
    let mut m = Module::new("spectroscopy");
    let mut rng = Rng::new(0x5ec7);
    let big = f64::MAX;
    let tiny = f64::from_bits(1);

    // ── absorbance (B) ──
    let ab = |m: &mut Module, e_: f64, l: f64, c: f64| {
        ff(m, "spectroscopy::absorbance", Class::B, &[e_, l, c], absorbance(bb(e_), bb(l), bb(c)));
    };
    for x in edges() {
        ab(&mut m, x, 1.0, 0.01);
        ab(&mut m, 100.0, x, 0.01);
        ab(&mut m, 100.0, 1.0, x);
    }
    for (e_, l, c) in [
        (100.0, 1.0, 0.01), // the Rust tests' values (A = 1)
        (1e200, 1e200, 1e-200),
        (1e-200, 1e-200, 1e200),
        (big, 2.0, 0.5),
        (tiny, 0.5, 2.0),
        (0.0, f64::INFINITY, 1.0),
    ] {
        ab(&mut m, e_, l, c);
    }
    for _ in 0..400 {
        ab(&mut m, rng.log_uniform(-3.0, 6.0), rng.log_uniform(-3.0, 2.0), rng.log_uniform(-12.0, 1.0));
    }
    for _ in 0..100 {
        ab(&mut m, rng.signed_log_uniform(-300.0, 300.0), rng.signed_log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0));
    }

    // ── transmittance (C: 10^(−A)) ──
    let tr = |m: &mut Module, a: f64| ff(m, "spectroscopy::transmittance", Class::C, &[a], transmittance(bb(a)));
    for x in edges() {
        tr(&mut m, x);
    }
    for k in -330..=330 {
        tr(&mut m, f64::from(k)); // the integers: 10⁰ … 10²² exact, 10⁻ᵏ inexact
    }
    for k in 0..=800 {
        tr(&mut m, -309.0 + 0.8 * f64::from(k)); // A over exp10's range, the overflow at A < −308.25 included
    }
    // subnormal results: A from 307.65 to 323.31 (10^(−A) below 2^-1022 and down to the smallest subnormal)
    for k in 0..=400 {
        tr(&mut m, 307.6 + 0.04 * f64::from(k));
    }
    for a in [0.75, 1.0, 0.0, -0.0, 2.0, std::f64::consts::LOG10_2] {
        tr(&mut m, a); // the Rust tests' absorbances; log₁₀ 2
    }
    for _ in 0..500 {
        tr(&mut m, rng.uniform(-5.0, 10.0));
    }
    for _ in 0..200 {
        tr(&mut m, rng.signed_log_uniform(-320.0, 3.0));
    }

    // ── absorbance_from_transmittance (C: log10) ──
    let aft = |m: &mut Module, t: f64| {
        fr(m, "spectroscopy::absorbance_from_transmittance", Class::C, &[t], absorbance_from_transmittance(bb(t)));
    };
    for x in edges() {
        aft(&mut m, x);
    }
    for x in decades(-323, 0, &["1", "2", "5", "7.5"]) {
        aft(&mut m, x);
    }
    let mut t1 = 1.0_f64;
    for _ in 0..40 {
        aft(&mut m, t1); // 1 and the ulps below (A ≈ 0, its cancellation)
        t1 = t1.next_down();
    }
    for x in [1.0_f64.next_up(), 1.5, 2.0, 0.1, 0.01, 0.001, pow2(-1074), pow2(-1022), 0.5, 0.25] {
        aft(&mut m, x);
    }
    for k in -20..=330 {
        let a = 0.5 * f64::from(k);
        aft(&mut m, transmittance(a)); // the round trips T(A)
    }
    for _ in 0..500 {
        aft(&mut m, rng.uniform(0.0, 1.0));
    }
    for _ in 0..300 {
        aft(&mut m, rng.log_uniform(-323.0, 0.0));
    }

    // ── concentration_from_absorbance (B) ──
    let conc = |m: &mut Module, a: f64, e_: f64, l: f64| {
        fr(m, "spectroscopy::concentration_from_absorbance", Class::B, &[a, e_, l], concentration_from_absorbance(bb(a), bb(e_), bb(l)));
    };
    for x in edges() {
        conc(&mut m, x, 100.0, 1.0);
        conc(&mut m, 1.0, x, 1.0);
        conc(&mut m, 1.0, 100.0, x);
    }
    for (a, e_, l) in [
        (1.0, 100.0, 1.0), // the Rust tests' values
        (1.0, f64::NAN, 1.0),
        (1.0, 100.0, f64::NAN),
        (1.0, f64::INFINITY, 1.0),
        (1e300, 1e200, 1e200), // ε·l overflows: Rust 0
        (1e-300, 1e-200, 1e-200), // ε·l underflows to 0: Rust ∞
        (1.0, 1e-170, 1e-170),
        (-1.0, 1e-160, 1e-160),
        (0.0, 1e-200, 1e-200),
        (big, 0.5, 0.5),
    ] {
        conc(&mut m, a, e_, l);
    }
    for _ in 0..400 {
        conc(&mut m, rng.signed_log_uniform(-4.0, 1.0), rng.log_uniform(-1.0, 6.0), rng.log_uniform(-2.0, 2.0));
    }
    for _ in 0..100 {
        conc(&mut m, rng.signed_log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0));
    }

    // ── photon energy, wavelength, frequency; eV ↔ J (B) ──
    for x in edges() {
        photon(&mut m, x);
    }
    for x in decades(-323, 308, &["1", "5"]) {
        photon(&mut m, x);
    }
    for k in 0..=100 {
        photon(&mut m, (380.0 + 4.0 * f64::from(k)) * 1e-9); // the visible range, 380–780 nm
    }
    for x in [500e-9, 700e-9, 400e-9, 13.6, 5e14, 1.0, 2.0 * 1.602_176_634e-19] {
        photon(&mut m, x); // the Rust tests' values
    }
    for _ in 0..300 {
        photon(&mut m, rng.log_uniform(-15.0, 20.0));
    }
    for _ in 0..200 {
        photon(&mut m, rng.signed_log_uniform(-323.0, 308.0));
    }

    // ── Bohr levels and transitions (B) ──
    for z in 0..=255u8 {
        for n in U32S {
            bohr(&mut m, z, n);
        }
    }
    let pairs: Vec<u32> = vec![0, 1, 2, 3, 4, 5, 6, 10, 100, 1000, 1_000_000, u32::MAX - 1, u32::MAX];
    for &z in &ZS {
        for &lo in &pairs {
            for &hi in &pairs {
                transition(&mut m, z, lo, hi);
            }
        }
    }
    for (z, lo, hi) in [(1, 1, 3), (1, 3, 1), (1, 2, 2), (1, 0, 2), (1, 1, 2), (1, 2, 3), (1, 3, 2), (2, 1, 2)] {
        transition(&mut m, z, lo, hi); // the Rust tests' values
    }
    // adjacent and close large n (1/n₁² − 1/n₂² and E₂ − E₁ cancel)
    for n in [10_u32, 50, 100, 137, 500, 1000, 10_000, 100_000, 1 << 20, 1 << 26, 94_906_265, 94_906_266, 1 << 31] {
        for d in [1_u32, 2, 3, 7] {
            for z in [1_u8, 2, 255] {
                transition(&mut m, z, n, n + d);
            }
        }
    }
    for _ in 0..300 {
        let z = u8::try_from(rng.int(1, 255)).unwrap_or(1);
        let lo = u32::try_from(rng.int(1, 10_000)).unwrap_or(1);
        let hi = lo + u32::try_from(rng.int(1, 10_000)).unwrap_or(1);
        transition(&mut m, z, lo, hi);
    }

    // ── SpectralSeries (A n_lower, B wavelength) ──
    for s in SERIES {
        m.row("SpectralSeries::n_lower", Class::A, &[e("SpectralSeries", s)], &[i(s.n_lower())]);
    }
    for s in SERIES {
        for line in (0..=24).chain([50, 100, 1000, 1_000_000, 1 << 31]) {
            series_line(&mut m, s, line);
        }
        for d in 0..=8 {
            series_line(&mut m, s, u32::MAX - d); // the checked_add overflow: n_lower + line past u32::MAX
        }
    }

    // ── wavenumber conversions (B; ADR 0002 #19) ──
    for x in edges() {
        wavenumbers(&mut m, x);
    }
    for x in decades(-323, 308, &["1", "2.5", "5"]) {
        wavenumbers(&mut m, x);
    }
    for k in -1074..=1023 {
        wavenumbers(&mut m, pow2(k));
    }
    for x in [500e-9, 10e-6, 20_000.0, 1000.0, 0.01, 0.1, 100.0] {
        wavenumbers(&mut m, x); // the Rust tests' values
    }
    // the subnormal-result end (x past 2^1018) and the overflow end (x below 2^-1029)
    for _ in 0..300 {
        wavenumbers(&mut m, pow2(1015) * rng.uniform(1.0, 256.0));
        wavenumbers(&mut m, f64::from_bits(u64::try_from(rng.int(1, 1 << 52)).unwrap_or(1)));
    }
    for _ in 0..3000 {
        let x = f64::from_bits(rng.int(1, 0x7FEF_FFFF_FFFF_FFFF).cast_unsigned());
        wavenumbers(&mut m, x);
    }
    for _ in 0..500 {
        wavenumbers(&mut m, rng.log_uniform(-8.0, 5.0)); // IR to radio wavelengths, and wavenumbers
    }

    // ── reduced_mass (B) ──
    let rm = |m: &mut Module, a: f64, b_: f64| fr(m, "spectroscopy::reduced_mass", Class::B, &[a, b_], reduced_mass(bb(a), bb(b_)));
    for x in edges() {
        rm(&mut m, x, 1.0);
        rm(&mut m, 1.0, x);
        rm(&mut m, x, x);
    }
    for (a, b_) in [
        (1.0, 1.0), // the Rust tests' values
        (1.0, 1e6),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (1.007_825, 34.968_853), // H³⁵Cl in u
        (1.6735e-27, 5.8067e-26),
        (1e200, 1e200),
        (1e-200, 1e-200),
        (big, big),
        (big, 1.0),
        (tiny, tiny),
        (tiny, 1.0),
    ] {
        rm(&mut m, a, b_);
    }
    for _ in 0..300 {
        let x = rng.log_uniform(-30.0, 3.0);
        rm(&mut m, x, x); // a homonuclear diatomic: Rust's m·m/(m + m), one ulp from m/2 on some masses
    }
    for _ in 0..400 {
        rm(&mut m, rng.log_uniform(-28.0, 3.0), rng.log_uniform(-28.0, 3.0));
    }
    for _ in 0..200 {
        rm(&mut m, rng.log_uniform(-320.0, 308.0), rng.log_uniform(-320.0, 308.0));
    }

    // ── harmonic_oscillator_energy (B) ──
    for x in edges() {
        oscillator(&mut m, 0, x);
        oscillator(&mut m, 1, x);
        oscillator(&mut m, u32::MAX, x);
    }
    for n in U32S {
        oscillator(&mut m, n, 1e14); // the Rust tests' ω
        oscillator(&mut m, n, 5.634e14);
    }
    for x in decades(-320, 300, &["1"]) {
        oscillator(&mut m, 3, x);
    }
    for _ in 0..300 {
        oscillator(&mut m, u32::try_from(rng.int(0, 1000)).unwrap_or(0), rng.log_uniform(10.0, 16.0));
    }

    // ── vibrational_frequency (B) ──
    let vf = |m: &mut Module, k: f64, mu: f64| fr(m, "spectroscopy::vibrational_frequency", Class::B, &[k, mu], vibrational_frequency(bb(k), bb(mu)));
    for x in edges() {
        vf(&mut m, x, 1.66e-27);
        vf(&mut m, 500.0, x);
    }
    for (k, mu) in [(500.0, 1.66e-27), (516.0, 1.6267e-27), (f64::NAN, 1.0), (1.0, f64::NAN), (1e300, 1e-100), (1e-300, 1e100), (big, tiny)] {
        vf(&mut m, k, mu); // the Rust tests' values (spectroscopy's and potential's), k/μ over- and underflowing
    }
    for x in decades(-300, 300, &["1", "5"]) {
        vf(&mut m, x, 1.66e-27);
        vf(&mut m, 500.0, x);
    }
    for _ in 0..400 {
        vf(&mut m, rng.log_uniform(0.0, 4.0), rng.log_uniform(-28.0, -24.0));
    }
    for _ in 0..200 {
        vf(&mut m, rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0));
    }

    // ── rigid_rotor_energy (B) ──
    for x in edges() {
        rotor(&mut m, 0, x);
        rotor(&mut m, 1, x);
        rotor(&mut m, u32::MAX, x);
    }
    for j in U32S {
        rotor(&mut m, j, 1e10); // the Rust tests' B
        rotor(&mut m, j, 3.17e11);
    }
    for x in decades(-320, 300, &["1"]) {
        rotor(&mut m, 2, x);
    }
    for _ in 0..300 {
        rotor(&mut m, u32::try_from(rng.int(0, 100)).unwrap_or(0), rng.log_uniform(8.0, 13.0));
    }

    // ── rotational_constant (B) ──
    let rc = |m: &mut Module, x: f64| fr(m, "spectroscopy::rotational_constant", Class::B, &[x], rotational_constant(bb(x)));
    for x in edges() {
        rc(&mut m, x);
    }
    for x in decades(-323, 308, &["1", "2.6"]) {
        rc(&mut m, x);
    }
    for x in [1e-46, 0.0, f64::NAN, 2.6e-47, 1.45e-46] {
        rc(&mut m, x); // the Rust tests' values; HCl's and CO's moments
    }
    for _ in 0..300 {
        rc(&mut m, rng.log_uniform(-48.0, -44.0));
    }

    // ── diatomic_moment_of_inertia (B) ──
    let mi = |m: &mut Module, mu: f64, r: f64| {
        ff(m, "spectroscopy::diatomic_moment_of_inertia", Class::B, &[mu, r], diatomic_moment_of_inertia(bb(mu), bb(r)));
    };
    for x in edges() {
        mi(&mut m, x, 1e-10);
        mi(&mut m, 1e-26, x);
    }
    for (mu, r) in [(1e-26, 1e-10), (1.6267e-27, 1.2746e-10), (1e200, 1e100), (1e-200, 1e-100)] {
        mi(&mut m, mu, r); // the Rust tests' values; HCl
    }
    for _ in 0..300 {
        mi(&mut m, rng.log_uniform(-28.0, -24.0), rng.log_uniform(-11.0, -9.0));
    }
    for _ in 0..100 {
        mi(&mut m, rng.signed_log_uniform(-300.0, 300.0), rng.signed_log_uniform(-200.0, 200.0));
    }

    // ── isotope_abundances (A): every u8 ──
    for z in 0..=255u8 {
        let a = isotope_abundances(bb(z));
        let mut outs = vec![Tok::V(a.len())];
        for &(mass, ab_) in a {
            outs.push(f(mass));
            outs.push(f(ab_));
        }
        m.row("spectroscopy::isotope_abundances", Class::A, &[i(z)], &outs);
    }

    // ── isotope_pattern (B) ──
    pattern(&mut m, &[]);
    for atoms in [
        &[(6u8, 1u32), (1, 4)][..], // the Rust tests' molecules: CH₄, Cl₂, H₈₈, C₁₆₈, the large one
        &[(17, 2)][..],
        &[(1, 88)][..],
        &[(6, 168)][..],
        &[(6, 1000), (1, 2000), (8, 500), (16, 50), (35, 20)][..],
        &[(0, 1)][..],
        &[(6, 1), (119, 2)][..],
        &[(26, 2)][..],
        &[(6, 6), (1, 12), (8, 6)][..], // glucose (the bench)
    ] {
        pattern(&mut m, atoms);
    }
    for z in 0..=255u8 {
        pattern(&mut m, &[(z, 1)]);
    }
    for z in [1u8, 6, 7, 8, 16, 17, 35] {
        for c in [0u32, 1, 2, 3, 4, 5, 7, 10, 20, 50, 100, 200, 500, 1000, 2000] {
            pattern(&mut m, &[(z, c)]);
        }
    }
    for atoms in [
        &[(1u8, 5000u32)][..], // the pruning bounds the peak count: about 7.4·σ peaks survive
        &[(6, 5000)][..],
        &[(35, 3000)][..],
        &[(17, 3000)][..],
        &[(16, 3000)][..],
        &[(8, 3000)][..],
        &[(26, u32::MAX)][..], // elements without data: one shift by count × the standard atomic weight
        &[(92, u32::MAX)][..],
        &[(118, 1000)][..],
        &[(26, 0)][..],
        &[(6, 0)][..],
        &[(6, 0), (1, 0)][..],
        &[(26, 1), (6, 2)][..], // a shift before and after a convolution
        &[(6, 2), (26, 1)][..],
        &[(6, 1), (6, 1)][..], // repeated entries against one entry of count 2
        &[(6, 2)][..],
        &[(1, 2)][..], // equal scratch masses (1H + 2H, 2H + 1H): the stable order decides the bin's mean
        &[(1, 3), (8, 1)][..],
        &[(17, 1), (35, 1)][..],
        &[(35, 1), (17, 1)][..],
        &[(6, 1), (0, 1)][..], // the Err position: the first bad atom in list order
        &[(0, 1), (6, 1)][..],
        &[(200, 1), (119, 1)][..],
        &[(6, 1), (1, 4), (8, 1), (7, 1), (16, 1), (17, 1), (35, 1), (26, 1), (53, 1), (11, 1)][..],
        &[(6, 60), (1, 98), (7, 16), (8, 18), (16, 1)][..], // a peptide-sized molecule
    ] {
        pattern(&mut m, atoms);
    }
    let zs = [1u8, 6, 7, 8, 16, 17, 35, 26, 11, 53, 92, 0, 119, 200];
    for _ in 0..300 {
        let len = usize::try_from(rng.int(1, 5)).unwrap_or(1);
        let atoms: Vec<(u8, u32)> = (0..len)
            .map(|_| {
                let z = zs[usize::try_from(rng.int(0, 10)).unwrap_or(0)]; // mostly valid; 0 / 119 / 200 now and then
                let z = if rng.unit() < 0.04 { zs[usize::try_from(rng.int(11, 13)).unwrap_or(11)] } else { z };
                (z, u32::try_from(rng.int(0, 60)).unwrap_or(0))
            })
            .collect();
        pattern(&mut m, &atoms);
    }

    // ── fft_spectrum (C) ──
    for n in 0..=16 {
        let sig: Vec<f64> = (0..n).map(|k| 1.0 + f64::from(k)).collect();
        spectrum(&mut m, &sig); // every length 0 … 16: Err off the powers of two
    }
    spectrum(&mut m, &[1.0; 10]); // the Rust test's non-power of two
    for k in 5..=9 {
        let n = 1usize << k;
        for len in [n - 1, n + 1] {
            spectrum(&mut m, &vec![0.5; len]); // 2^k ± 1
        }
    }
    spectrum(&mut m, &vec![0.5; (1 << 12) + 1]);
    for k in 0..=13u32 {
        let n = 1u64 << k;
        let nu = usize::try_from(n).unwrap_or(0);
        let mut kinds: Vec<Vec<f64>> = Vec::new();
        let sin_c = |c: u64| -> Vec<f64> { (0..n).map(|j| sine(c, j, n)).collect() };
        let random = |rng: &mut Rng| -> Vec<f64> { (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect() };
        if k <= 10 {
            let mut imp = vec![0.0; nu];
            imp[0] = 1.0;
            kinds.push(imp); // an impulse at 0: a flat spectrum
            let mut imp2 = vec![0.0; nu];
            imp2[nu / 2] = 2.5;
            kinds.push(imp2);
            kinds.push(vec![1.0; nu]); // a constant: DC only
            kinds.push(sin_c((n / 8).max(1))); // a sinusoid on a bin
            kinds.push(random(&mut rng));
            kinds.push((0..n).map(|_| rng.signed_log_uniform(299.0, 305.0)).collect()); // huge: the sums overflow
            kinds.push((0..n).map(|_| rng.signed_log_uniform(-323.0, -309.0)).collect()); // subnormal values
            let mut nan = random(&mut rng);
            nan[nu - 1] = f64::NAN;
            kinds.push(nan);
            let mut inf = random(&mut rng);
            inf[0] = f64::INFINITY;
            kinds.push(inf);
        } else if k <= 12 {
            let mut imp = vec![0.0; nu];
            imp[1] = 1.0;
            kinds.push(imp);
            kinds.push(sin_c(n / 8));
            kinds.push(random(&mut rng));
        } else {
            kinds.push(sin_c(n / 8));
            kinds.push(random(&mut rng));
        }
        if k <= 6 {
            let mut imp3 = vec![0.0; nu];
            imp3[usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)).unwrap_or(0)] = -3.0;
            kinds.push(imp3);
            kinds.push(vec![3.7; nu]);
            kinds.push(vec![0.0; nu]);
            kinds.push(vec![-0.0; nu]); // ±0: the signs of the zero magnitudes' sums
            kinds.push((0..n).map(|j| if j % 2 == 0 { 1.0 } else { -1.0 }).collect()); // the Nyquist bin
            kinds.push((0..n).map(|j| sine(1, j, n) + 0.5 * sine(3, j, n)).collect()); // two tones, off the axes
            kinds.push((0..n).map(|j| sine(1, j, n) * 1e300).collect());
            kinds.push((0..n).map(|j| if j == 0 { f64::NEG_INFINITY } else if j == 1 { f64::INFINITY } else { 0.0 }).collect());
            kinds.push(vec![f64::NAN; nu]);
            kinds.push((0..n).map(|j| if j % 3 == 0 { big } else { -big }).collect());
            kinds.push((0..n).map(|_| tiny).collect());
        }
        for sig in kinds {
            spectrum(&mut m, &sig);
        }
    }

    // ── find_peaks (A) ──
    let th = [0.0, 1.0, -0.0, 0.5, 0.1, 0.25, 0.75, 1.0_f64.next_down(), tiny, -0.1, 1.5, f64::NAN, f64::INFINITY,
              f64::NEG_INFINITY, -tiny, 1.0_f64.next_up()];
    let basic: [&[f64]; 10] = [
        &[0.0, 1.0, 3.0, 1.0, 0.5, 2.0, 0.5, 0.0], // the Rust tests' data
        &[0.0, 1.0, 3.0, 1.0, 0.5, 0.6, 0.5, 0.0],
        &[1.0, 2.0],
        &[],
        &[-5.0, -3.0, -1.0, -3.0, -5.0, -4.0, -5.0],
        &[10.0, 10.2, 10.0, 11.0, 10.0],
        &[0.0, 1.0, 3.0, 1.0, 0.0],
        &[0.0, f64::NAN, 0.0, 2.0, 0.0, 1.0, 0.0],
        &[f64::NAN; 4],
        &[5.0],
    ];
    for d in basic {
        for &t in &th {
            peaks(&mut m, d, t);
        }
    }
    let shapes: [&[f64]; 18] = [
        &[0.0, 1.0, 1.0, 0.0], // plateaus: never a strict maximum
        &[0.0, 1.0, 1.0, 1.0, 0.0],
        &[0.0, 1.0, 2.0, 2.0, 1.0, 3.0, 0.0],
        &[2.0, 2.0, 2.0, 2.0],
        &[1.0, 2.0, 3.0, 4.0, 5.0],
        &[5.0, 4.0, 3.0, 2.0, 1.0],
        &[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        &[0.0, 1.0, 0.0],
        &[1.0, 0.0, 1.0],
        &[0.0, f64::INFINITY, 0.0, 1.0, 0.0], // ±∞ in the data: compared, but not in min / max
        &[0.0, f64::NEG_INFINITY, 0.0, 1.0, 0.0],
        &[f64::INFINITY, f64::INFINITY, f64::INFINITY],
        &[f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY],
        &[f64::NAN, 1.0, f64::NAN],
        &[-0.0, 0.0, -0.0, 0.0, -0.0], // +0 > −0 is false
        &[-0.0, 5e-324, 0.0],
        &[-1.797_693_134_862_315_7e308, 1.797_693_134_862_315_7e308, -1.797_693_134_862_315_7e308, 0.0],
        &[1.797_693_134_862_315_7e308, 1.797_693_134_862_315_7e308, f64::INFINITY, 1.797_693_134_862_315_7e308],
    ];
    for d in shapes {
        for &t in &th {
            peaks(&mut m, d, t);
        }
    }
    // thresholds where min·(1 − t) + max·t lands on a data value, to the ulp
    let ladder = [0.0, 1.0, 0.25, 0.5, 0.75, 1.0, 0.0, 0.5000000000000001, 0.0, 0.4999999999999999, 0.0];
    for k in 0..=40 {
        peaks(&mut m, &ladder, f64::from(k) / 40.0);
    }
    for t in [0.5_f64.next_down(), 0.5, 0.5_f64.next_up(), 0.25_f64.next_down(), 0.25_f64.next_up()] {
        peaks(&mut m, &ladder, t);
    }
    for _ in 0..400 {
        let len = usize::try_from(rng.int(0, 40)).unwrap_or(0);
        let data: Vec<f64> = (0..len)
            .map(|_| {
                let u = rng.unit();
                if u < 0.02 {
                    f64::NAN
                } else if u < 0.03 {
                    f64::INFINITY
                } else if u < 0.04 {
                    f64::NEG_INFINITY
                } else if u < 0.2 {
                    f64::from(u32::try_from(rng.int(0, 4)).unwrap_or(0)) // small integers: ties and plateaus
                } else {
                    rng.uniform(-10.0, 10.0)
                }
            })
            .collect();
        let t = if rng.unit() < 0.1 { th[usize::try_from(rng.int(0, 15)).unwrap_or(0)] } else { rng.unit() };
        peaks(&mut m, &data, t);
    }
    for _ in 0..100 {
        let len = usize::try_from(rng.int(3, 30)).unwrap_or(3);
        let data: Vec<f64> = (0..len).map(|_| rng.signed_log_uniform(-300.0, 308.0)).collect();
        peaks(&mut m, &data, rng.unit());
    }
    m
}
