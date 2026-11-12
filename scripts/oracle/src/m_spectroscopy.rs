//! `spectroscopy.tsv` — the Beer–Lambert law, photon energy conversions, the Bohr model and the Rydberg formula,
//! vibrational and rotational levels, isotope patterns, the FFT magnitude spectrum and peak detection (port plan §4,
//! M4e).
//!
//! Rust's u8 / u32 arguments are `i:` tokens, `SpectralSeries` an `e:SpectralSeries.<Variant>` token. A slice is `v:n`
//! and its n elements: `isotope_pattern`'s `&[(u8, u32)]` as n `i:z i:count` pairs, `fft_spectrum`'s and `find_peaks`'
//! `&[f64]` as n `f:` tokens (`find_peaks`' threshold follows). A Vec result is `v:n` and its elements: the isotope
//! pattern's (mass, intensity) pairs, the spectrum's floats, the peaks' `i:` indices; `isotope_abundances`' slice the
//! same way.
//! - Class A: `isotope_abundances` (every u8; **ADR 0002 #23**, M4e-4: the 7 elements' rows noted with the port's slice from
//!   the cited AME2020 / TICE 2009 table), `SpectralSeries::n_lower` (every variant), `find_peaks` (integers: the
//!   Rust tests' data, fewer than 3 points, thresholds 0, 1, −0.0, the ulps inside [0, 1] and every way out of it, NaN
//!   and ±∞; plateaus, flat and monotone data, NaN and ±∞ in the data, all non-finite, ±0, DBL_MAX ranges where the
//!   interpolated threshold over- or underflows, seeded data).
//! - Class B (+ − × ÷ and √, expected 0 ulp): `absorbance`, `concentration_from_absorbance`, the photon conversions,
//!   `joules_to_ev` / `ev_to_joules`, `harmonic_oscillator_energy`, `vibrational_frequency`, `rigid_rotor_energy`,
//!   `rotational_constant`, `diatomic_moment_of_inertia` (each argument over the edges, the Rust tests' values, decades,
//!   the intermediates' over- and underflow, seeded); **ADR 0002 #21** (M4e-3): `bohr_energy_level` (every u8 Z × 20 n,
//!   and every exact level, n = 2^k·Z), `bohr_transition_energy` and `rydberg_wavelength` (Z × (n₁, n₂) grids, every Err
//!   path in Rust's order, adjacent large n, n₁ = n₂ − 1 at every binary scale and next to u32::MAX, the fast path's
//!   2¹³ edge, the exhaustive search's 13 ties and its exact values), `SpectralSeries::wavelength` (every series × lines
//!   to u32::MAX, Rust's `checked_add` overflow, large lines), each row noted with the port's correctly rounded exact
//!   rational where it differs from Rust's; **ADR 0002 #22**: `reduced_mass` (equal masses at every power of two and
//!   subnormal multiple, three families of normal ties and two of subnormal ones, exponent gaps around the proven 55,
//!   the overflow and underflow ends, seeded over every binade), noted likewise; `isotope_pattern` (the Rust tests'
//!   molecules, the empty molecule, every u8 Z alone, every element with data × counts to 5000 (the pruning keeps the
//!   peak count bounded), elements without data with counts to u32::MAX, H₂… (equal scratch masses: the stable sort's
//!   order decides the binned mean), mixtures, the Err position in a list, seeded molecules; **ADR 0002 #23**: each row
//!   whose result changes with the AME2020 / TICE 2009 table noted with the port's, from [`pattern_replica`], Rust's
//!   function transcribed with the table a parameter and checked against the oracle on every row); and **ADR 0002 #19**:
//!   `wavenumber_from_wavelength` and `wavelength_from_wavenumber` (the edges, every decade, every power of two, the
//!   subnormal-result and overflow ends, seeded over every binade), the rows where the port's correctly rounded
//!   1/(100·x) differs from Rust's `0.01 / x` noted with it; **ADR 0002 #30** (the overflow audit's B2): `absorbance`,
//!   `concentration_from_absorbance`, `rigid_rotor_energy`, `rotational_constant` and `diatomic_moment_of_inertia` rows
//!   noted where the port's correctly rounded monomial differs from Rust's value, `vibrational_frequency`'s where its
//!   correctly rounded √(k/μ) does (`audit.rs`), and `audit::emit_spectroscopy`'s edge rows appended with their own seed
//!   (k/μ past DBL_MAX or below 2⁻¹⁰⁷⁴, h·B and 8π²·I subnormal, partial products past either end, midpoints).
//! - Class C (bound 0 against the corrected oracle): `transmittance` (`10.0.powf(−A)`: A over exp10's whole range, the
//!   integers, subnormal results, seeded), `absorbance_from_transmittance` (log10: every decade of T, T near 1, the
//!   round trips, seeded), `fft_spectrum` (cos and sin twiddles: every length 0 … 16 and every power of two to 2¹³,
//!   2^k ± 1, with impulses, constants, ±0, alternating and Nyquist signals, sinusoids, random, huge and subnormal
//!   values, NaN and ±∞).
//!
//! Every input goes through `black_box`; grid inputs use only exact or correctly rounded operations (`grids.rs`): the
//! sinusoids are a Taylor series in + − × ÷ (deterministic, not std's sin).

use std::hint::black_box as bb;

use crate::audit;
use crate::exact::{cr_bohr_level, cr_hundredth_div, cr_reduced_mass, cr_rydberg, cr_transition};
use crate::grids::{Rng, decades, pow2};
use crate::m_reaction::{edges, res};
use crate::t_spectroscopy::{Iso, new_abundances};
use crate::tsv::{Class, Module, Tok, e, err, f, i, tokens_note};
use kimiya::spectroscopy::{
    SpectralSeries, absorbance_from_transmittance, bohr_energy_level, bohr_transition_energy, ev_to_joules,
    fft_spectrum, find_peaks, frequency_from_wavelength, harmonic_oscillator_energy, isotope_abundances,
    isotope_pattern, joules_to_ev, photon_energy_from_frequency, photon_energy_from_wavelength, reduced_mass,
    rydberg_wavelength, transmittance, wavelength_from_energy, wavelength_from_wavenumber,
    wavenumber_from_wavelength, RYDBERG, RYDBERG_EV,
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

/// ADR 0002 #21 / #22: a Result row, noted with the port's correctly rounded value where it differs from the oracle's
/// (`what` names the exact value: "−R·Z²/n²", …).
fn noted(m: &mut Module, name: &str, ins: &[Tok], out: Result<f64, kimiya::error::KimiyaError>, exact: Option<f64>,
         adr: u32, what: &str) {
    if let (Ok(v), Some(c)) = (&out, exact)
        && c.to_bits() != v.to_bits()
    {
        let note = format!("ADR 0002 #{adr}: the port returns the correctly rounded {what} = {c:?} (f:{:016x})", c.to_bits());
        m.row_noted(name, Class::B, ins, &[f(*v)], &note);
        return;
    }
    m.row(name, Class::B, ins, &[res(out)]);
}

fn bohr(m: &mut Module, z: u8, n: u32) {
    let out = bohr_energy_level(bb(z), bb(n));
    let exact = (z > 0 && n > 0).then(|| -cr_bohr_level(RYDBERG_EV, z, n));
    noted(m, "spectroscopy::bohr_energy_level", &[i(z), i(n)], out, exact, 21, "−R·Z²/n²");
}

fn transition(m: &mut Module, z: u8, lo: u32, hi: u32) {
    let valid = z > 0 && lo > 0 && hi > lo;
    let t = bohr_transition_energy(bb(z), bb(lo), bb(hi));
    let et = valid.then(|| cr_transition(RYDBERG_EV, z, lo, hi));
    noted(m, "spectroscopy::bohr_transition_energy", &[i(z), i(lo), i(hi)], t, et, 21, "R·Z²·(n₂² − n₁²)/(n₁²·n₂²)");
    let w = rydberg_wavelength(bb(z), bb(lo), bb(hi));
    let ew = valid.then(|| cr_rydberg(RYDBERG, z, lo, hi));
    noted(m, "spectroscopy::rydberg_wavelength", &[i(z), i(lo), i(hi)], w, ew, 21, "n₁²·n₂²/(R∞·Z²·(n₂² − n₁²))");
}

fn series_line(m: &mut Module, s: SpectralSeries, line: u32, rh: f64) {
    let out = bb(s).wavelength(bb(line));
    let exact = s.n_lower().checked_add(line).filter(|_| line > 0).map(|hi| cr_rydberg(rh, 1, s.n_lower(), hi));
    noted(m, "SpectralSeries::wavelength", &[e("SpectralSeries", s), i(line)], out, exact, 21, "n₁²·n₂²/(R_H·(n₂² − n₁²))");
}

/// Rust's private RYDBERG_HYDROGEN through the API, as `t_spectroscopy` derives it: `Lyman.wavelength(u32::MAX − 1)` is
/// 1/R_H (1/1 − 1/u32::MAX² is 1.0 exactly), and R_H is the one double near 1/λ whose reciprocal rounds to λ.
fn rydberg_hydrogen() -> f64 {
    let lam = SpectralSeries::Lyman.wavelength(u32::MAX - 1).unwrap_or(f64::NAN);
    let mut x = 1.0 / lam;
    for _ in 0..4 {
        x = x.next_down();
    }
    let mut hits = Vec::new();
    for _ in 0..9 {
        if (1.0 / x).to_bits() == lam.to_bits() {
            hits.push(x);
        }
        x = x.next_up();
    }
    assert!(hits.len() == 1, "RYDBERG_HYDROGEN: {} doubles have 1/x = {lam:?}", hits.len());
    hits[0]
}

/// The prime factorisation of n (trial division by every prime below 2¹⁷, then a Miller–Rabin prime cofactor, which
/// must be below 2⁶⁴): (prime, exponent) pairs. For the tie search's odd significands and n₂ ∓ n₁ < 2³³.
fn factor(mut n: u64, primes: &[u64]) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    for &p in primes {
        if p * p > n {
            break;
        }
        let mut k = 0;
        while n.is_multiple_of(p) {
            n /= p;
            k += 1;
        }
        if k > 0 {
            out.push((p, k));
        }
    }
    if n > 1 {
        assert!(n < 1 << 34 || is_prime(n), "factor: composite cofactor {n} past the trial division");
        out.push((n, 1));
    }
    out
}

/// Deterministic Miller–Rabin for n < 2⁶⁴ (the first twelve primes as bases).
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let mulmod = |a: u64, b: u64| u64::try_from(u128::from(a) * u128::from(b) % u128::from(n)).unwrap_or(0);
    let powmod = |mut b: u64, mut e: u64| {
        let mut r = 1_u64;
        while e > 0 {
            if e & 1 == 1 {
                r = mulmod(r, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        r
    };
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'base: for a in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}

/// The odd part of a binary64's significand (a positive finite x).
fn odd_significand(x: f64) -> u64 {
    let (m, _) = crate::exact::decompose(x);
    m >> m.trailing_zeros()
}

/// A (Z, n₁, n₂) argument triple of the transition functions.
type Transition = (u8, u32, u32);

/// ADR 0002 #21's exhaustive search of `bohr_transition_energy`'s domain for exact values: R·Z²·(n₂² − n₁²)/(n₁²·n₂²)
/// with n₁ = g·x, n₂ = g·y, gcd(x, y) = 1, is dyadic iff the odd part of g²x²y² divides K·w, K = odd(R)·odd(Z)² and
/// w = odd(y² − x²); w is coprime to x and y, so odd(x)·odd(y) | odd(Z) (R's odd part, 683·11214237642539, is squarefree),
/// and with j = odd(Z)/(odd(x)·odd(y)) the odd significand is T/og², T = 683·11214237642539·j²·w, og = odd(g), og² | T.
/// So: every coprime odd (ox, oy) with ox·oy ≤ 255, every x = ox·2^a, y = oy·2^b (min(a, b) = 0, x < y < 2³²), every
/// Z with ox·oy | odd(Z), every og with og² | T and og·y < 2³², every g = og·2^c with g·y < 2³². Returns (ties, exact
/// values): (Z, n₁, n₂) whose odd significand has 54 bits (a midpoint) or at most 53 (a double).
fn transition_exact_values() -> (Vec<Transition>, Vec<Transition>) {
    let mut sieve = vec![true; 1 << 17];
    let mut primes = Vec::new();
    for k in 2..sieve.len() {
        if sieve[k] {
            primes.push(u64::try_from(k).unwrap_or(0));
            let mut j = k * k;
            while j < sieve.len() {
                sieve[j] = false;
                j += k;
            }
        }
    }
    let kr = odd_significand(RYDBERG_EV);
    let fr = factor(kr, &primes);
    assert!(fr.iter().all(|&(_, e)| e == 1), "RYDBERG_EV's odd significand {kr} is not squarefree: {fr:?}");
    let gcd = |mut a: u64, mut b: u64| {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    };
    let odd = |n: u64| n >> n.trailing_zeros();
    let (mut ties, mut exact) = (Vec::new(), Vec::new());
    let max = u64::from(u32::MAX);
    for ox in (1..=255_u64).step_by(2) {
        for oy in (1..=255_u64).step_by(2) {
            if ox * oy > 255 || gcd(ox, oy) != 1 {
                continue;
            }
            for a in 0..33_u32 {
                for b in 0..33_u32 {
                    if a > 0 && b > 0 {
                        continue;
                    }
                    let (x, y) = (ox << a, oy << b);
                    if x >= y || y > max {
                        continue;
                    }
                    let mut fw: Vec<(u64, u32)> = Vec::new();
                    for (p, e) in factor(y - x, &primes).into_iter().chain(factor(y + x, &primes)) {
                        if p == 2 {
                            continue;
                        }
                        if let Some(t) = fw.iter_mut().find(|t| t.0 == p) { t.1 += e } else { fw.push((p, e)) }
                    }
                    for z in 1..=255_u8 {
                        let zo = odd(u64::from(z));
                        if zo % (ox * oy) != 0 {
                            continue;
                        }
                        let j = zo / (ox * oy);
                        let mut ft = fr.clone();
                        for (p, e) in fw.iter().copied().chain(factor(j, &primes).into_iter().map(|(p, e)| (p, 2 * e))) {
                            if let Some(t) = ft.iter_mut().find(|t| t.0 == p) { t.1 += e } else { ft.push((p, e)) }
                        }
                        // og over the square divisors of T (exponent k ≤ e/2 per prime), with T/og² by exponents
                        let mut cands: Vec<(u64, Vec<u32>)> = vec![(1, vec![0; ft.len()])];
                        for (idx, &(p, e)) in ft.iter().enumerate() {
                            let mut more = Vec::new();
                            for (d, ks) in &cands {
                                let mut v = *d;
                                for k in 1..=e / 2 {
                                    v = v.saturating_mul(p);
                                    let mut ks2 = ks.clone();
                                    ks2[idx] = k;
                                    more.push((v, ks2));
                                }
                            }
                            cands.extend(more);
                        }
                        for (og, ks) in cands {
                            if og.saturating_mul(y) > max {
                                continue;
                            }
                            let mut sig: Option<u128> = Some(1);
                            for (idx, &(p, e)) in ft.iter().enumerate() {
                                let pe = u128::from(p).checked_pow(e - 2 * ks[idx]);
                                sig = sig.and_then(|v| pe.and_then(|q| v.checked_mul(q)));
                            }
                            let Some(sig) = sig else { continue };
                            let bits = 128 - sig.leading_zeros();
                            if bits > 54 {
                                continue;
                            }
                            let mut g = og;
                            while g * y <= max {
                                let row = (z, u32::try_from(g * x).unwrap_or(0), u32::try_from(g * y).unwrap_or(0));
                                if bits == 54 { ties.push(row) } else { exact.push(row) }
                                g <<= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    ties.sort_unstable();
    exact.sort_unstable();
    (ties, exact)
}

fn oscillator(m: &mut Module, n: u32, w: f64) {
    let out = harmonic_oscillator_energy(bb(n), bb(w));
    m.row("spectroscopy::harmonic_oscillator_energy", Class::B, &[i(n), f(w)], &[f(out)]);
}

fn rotor(m: &mut Module, j: u32, b: f64) {
    audit::rotor(m, j, b);
}

/// Rust 1.2.0's `isotope_pattern` (rust-old/src/spectroscopy.rs), transcribed operation for operation with the isotope
/// table a parameter (`table(z)` is z's (mass, abundance) pairs, empty for none) and its private `ISOTOPE_BIN_WIDTH_U`
/// (0.5), `ISOTOPE_PRUNE_REL` (10⁻¹²) and `ISOTOPE_REPORT_REL` (0.001) typed: the pruned convolution, the stable
/// `sort_by(total_cmp)`, the binning, `fold(0.0, max)` and the two `retain`s. With the oracle's `isotope_abundances` it
/// must give the oracle's result on every row (`pattern` checks it), so with ADR 0002 #23's table it gives the port's.
/// The Err is the variant's name.
fn pattern_replica(atoms: &[(u8, u32)], table: &dyn Fn(u8) -> Vec<(f64, f64)>) -> Result<Vec<(f64, f64)>, &'static str> {
    if atoms.is_empty() {
        return Err("InvalidInput");
    }
    let mut pattern: Vec<(f64, f64)> = vec![(0.0, 1.0)];
    let mut scratch: Vec<(f64, f64)> = Vec::new();
    for &(z, count) in atoms {
        let abundances = table(z);
        if abundances.is_empty() {
            let Some(elem) = kimiya::element::lookup_by_number(z) else {
                return Err("InvalidElement");
            };
            let shift = f64::from(count) * elem.atomic_mass;
            for peak in &mut pattern {
                peak.0 += shift;
            }
            continue;
        }
        for _ in 0..count {
            scratch.clear();
            for &(m1, i1) in &pattern {
                for &(m2, i2) in &abundances {
                    scratch.push((m1 + m2, i1 * i2));
                }
            }
            scratch.sort_by(|a, b| a.0.total_cmp(&b.0));
            pattern.clear();
            for &(m, i) in &scratch {
                if let Some(last) = pattern.last_mut()
                    && (m - last.0).abs() < 0.5
                {
                    let total = last.1 + i;
                    if total > 0.0 {
                        last.0 = (last.0 * last.1 + m * i) / total;
                    }
                    last.1 = total;
                    continue;
                }
                pattern.push((m, i));
            }
            let max_intensity = pattern.iter().fold(0.0_f64, |acc, p| acc.max(p.1));
            if max_intensity > 0.0 {
                for peak in &mut pattern {
                    peak.1 /= max_intensity;
                }
            }
            pattern.retain(|p| p.1 >= 1e-12);
        }
    }
    pattern.retain(|(_, i)| *i > 0.001);
    Ok(pattern)
}

/// A `Vec<(f64, f64)>` result's golden tokens: `v:n`, then each pair's two floats; an Err its `E:` token.
fn pairs_toks(r: Result<&[(f64, f64)], &str>) -> Vec<Tok> {
    match r {
        Ok(p) => {
            let mut v = vec![Tok::V(p.len())];
            for &(a, b) in p {
                v.push(f(a));
                v.push(f(b));
            }
            v
        }
        Err(name) => vec![Tok::Err(name.to_string())],
    }
}

/// Whether two token lists are the same bits (every float bit for bit; a NaN is no isotope value, so none is special).
fn same_toks(a: &[Tok], b: &[Tok]) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(x, y)| match (x, y) {
            (Tok::F(p), Tok::F(q)) => p.to_bits() == q.to_bits(),
            (Tok::V(p), Tok::V(q)) => p == q,
            (Tok::Err(p), Tok::Err(q)) => p == q,
            _ => false,
        })
}

/// ADR 0002 #23 (M4e-4): a row whose result is a `Vec<(f64, f64)>`, the oracle's outputs, noted with the port's result
/// from the cited table where the two differ.
fn row23(m: &mut Module, name: &str, cls: Class, ins: &[Tok], oracle: &[Tok], port: &[Tok]) {
    if same_toks(oracle, port) {
        m.row(name, cls, ins, oracle);
        return;
    }
    let note = format!(
        "ADR 0002 #23: with the AME2020 masses and the TICE 2009 abundances, ³⁶S included, the port returns {}",
        tokens_note(port)
    );
    m.row_noted(name, cls, ins, oracle, &note);
}

fn pattern(m: &mut Module, isos: &[Iso], atoms: &[(u8, u32)]) {
    let mut ins = vec![Tok::V(atoms.len())];
    for &(z, c) in atoms {
        ins.push(i(z));
        ins.push(i(c));
    }
    let outs = match isotope_pattern(bb(atoms)) {
        Ok(p) => pairs_toks(Ok(&p)),
        Err(e_) => vec![err(&e_)],
    };
    let rust = pattern_replica(atoms, &|z| isotope_abundances(z).to_vec());
    assert!(same_toks(&pairs_toks(rust.as_deref().map_err(|e| *e)), &outs), "isotope_pattern replica disagrees with the oracle: {atoms:?}");
    let port = pattern_replica(atoms, &|z| new_abundances(isos, z));
    row23(m, "spectroscopy::isotope_pattern", Class::B, &ins, &outs, &pairs_toks(port.as_deref().map_err(|e| *e)));
}

fn spectrum(m: &mut Module, signal: &[f64]) {
    let mut ins = vec![Tok::V(signal.len())];
    ins.extend(signal.iter().map(|&x| f(x)));
    let vals = |s: &[f64]| {
        let mut v = vec![Tok::V(s.len())];
        v.extend(s.iter().map(|&x| f(x)));
        v
    };
    let rust = fft_spectrum(bb(signal));
    let outs = match &rust {
        Ok(s) => vals(s),
        Err(e_) => vec![err(e_)],
    };
    // ADR 0002 #36 (the overflow audit's B6): the port's spectrum of a finite signal, the unbounded-exponent transform,
    // noted where it differs from Rust's (by its bits; NaN never occurs for a finite signal)
    if let (Ok(r), Some(p)) = (&rust, audit::port_fft(signal))
        && (r.len() != p.len() || r.iter().zip(&p).any(|(a, b)| a.to_bits() != b.to_bits()))
    {
        let note = format!("ADR 0002 #36: the port returns {} = {}", audit::FFT_HOW, tokens_note(&vals(&p)));
        m.row_noted("spectroscopy::fft_spectrum", Class::C, &ins, &outs, &note);
        return;
    }
    m.row("spectroscopy::fft_spectrum", Class::C, &ins, &outs);
}

fn peaks(m: &mut Module, data: &[f64], threshold: f64) {
    let mut ins = vec![Tok::V(data.len())];
    ins.extend(data.iter().map(|&x| f(x)));
    ins.push(f(threshold));
    let idx = |p: &[usize]| {
        let mut v = vec![Tok::V(p.len())];
        v.extend(p.iter().map(|&k| i(i64::try_from(k).unwrap_or(-1))));
        v
    };
    let rust = find_peaks(bb(data), bb(threshold));
    let outs = idx(&rust);
    // ADR 0002 #37 (the overflow audit's B6): the port's indices with the exact threshold test, noted where they differ
    let port = audit::port_peaks(data, threshold);
    if port == rust {
        m.row("spectroscopy::find_peaks", Class::A, &ins, &outs);
    } else {
        let note = format!("ADR 0002 #37: the port returns {} = {}", audit::PEAKS_HOW, tokens_note(&idx(&port)));
        m.row_noted("spectroscopy::find_peaks", Class::A, &ins, &outs, &note);
    }
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
pub fn emit(isos: &[Iso]) -> Module {
    let mut m = Module::new("spectroscopy");
    let mut rng = Rng::new(0x5ec7);
    let big = f64::MAX;
    let tiny = f64::from_bits(1);

    // ── absorbance (B) ──
    let ab = audit::absorbance_row;
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
    let conc = audit::conc;
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

    // ── Bohr levels and transitions (B; ADR 0002 #21) ──
    for z in 0..=255u8 {
        for n in U32S {
            bohr(&mut m, z, n);
        }
    }
    // every exact level: −R·Z²/n² is a double exactly when n/Z is a power of two (R's odd significand, 683 times a prime
    // above 2³², is squarefree, so odd(n)² must divide odd(Z)², and a square factor j² ≥ 9 left over puts the significand
    // at 56 bits or more: never a midpoint either)
    for z in 1..=255u8 {
        let mut n = u32::from(z);
        while n % 2 == 0 {
            n /= 2;
        }
        loop {
            if !U32S.contains(&n) {
                bohr(&mut m, z, n);
            }
            match n.checked_mul(2) {
                Some(v) => n = v,
                None => break,
            }
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
    // ADR 0002 #21's hard rows: n₁ = n₂ − 1 at every binary scale and next to u32::MAX, the fast path's n₂ = 2¹³ edge,
    // and the exhaustive search's exact values (each family's smallest, middle and largest g) and every tie
    for k in 13..=31_u32 {
        for n in [(1_u32 << k) - 1, 1 << k, (1 << k) + 1] {
            for z in [1_u8, 2, 3, 92, 255] {
                transition(&mut m, z, n, n + 1);
            }
        }
    }
    for k in 0..=9_u32 {
        for d in [1_u32, 2, 3, 7, 1000] {
            for z in [1_u8, 26, 255] {
                transition(&mut m, z, u32::MAX - k - d, u32::MAX - k);
            }
        }
    }
    for hi in [8190_u32, 8191, 8192, 8193] {
        for lo in [1_u32, hi / 2, hi - 2, hi - 1] {
            for z in [1_u8, 255] {
                transition(&mut m, z, lo, hi);
            }
        }
    }
    let (ties, exact) = transition_exact_values();
    assert!(ties.len() == 13 && exact.len() == 5657, "the #21 search found {} ties and {} exact values", ties.len(), exact.len());
    for &(z, lo, hi) in &ties {
        transition(&mut m, z, lo, hi);
    }
    let mut fams: Vec<(Transition, Vec<(u32, u32)>)> = Vec::new();
    for &(z, lo, hi) in &exact {
        let g = {
            let (mut a, mut b) = (lo, hi);
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a
        };
        match fams.iter_mut().find(|f| f.0 == (z, lo / g, hi / g)) {
            Some(f) => f.1.push((lo, hi)),
            None => fams.push(((z, lo / g, hi / g), vec![(lo, hi)])),
        }
    }
    for ((z, _, _), mut members) in fams {
        members.sort_unstable();
        let last = members.len() - 1;
        let mut picks = vec![members[0], members[last / 2], members[last]];
        picks.dedup();
        for (lo, hi) in picks {
            transition(&mut m, z, lo, hi);
        }
    }

    // ── SpectralSeries (A n_lower, B wavelength; ADR 0002 #21) ──
    let rh = rydberg_hydrogen();
    for s in SERIES {
        m.row("SpectralSeries::n_lower", Class::A, &[e("SpectralSeries", s)], &[i(s.n_lower())]);
    }
    for s in SERIES {
        for line in (0..=24).chain([50, 100, 1000, 1_000_000, 1 << 31]) {
            series_line(&mut m, s, line, rh);
        }
        for d in 0..=8 {
            series_line(&mut m, s, u32::MAX - d, rh); // the checked_add overflow: n_lower + line past u32::MAX
        }
        // large lines (n₂ = n₁ + line: 1/n₁² − 1/n₂² near 1/n₁²), the largest, and the fast path's n₂ = 2¹³ edge
        for k in 5..=30_u32 {
            series_line(&mut m, s, 1 << k, rh);
        }
        for d in 9..=24 {
            series_line(&mut m, s, u32::MAX - d, rh);
        }
        for n2 in [8191_u32, 8192, 8193] {
            series_line(&mut m, s, n2 - s.n_lower(), rh);
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

    // ── reduced_mass (B; ADR 0002 #22) ──
    let rm = |m: &mut Module, a: f64, b_: f64| {
        let ok = a > 0.0 && a.is_finite() && b_ > 0.0 && b_.is_finite();
        let exact = ok.then(|| cr_reduced_mass(a, b_));
        noted(m, "spectroscopy::reduced_mass", &[f(a), f(b_)], reduced_mass(bb(a), bb(b_)), exact, 22, "m₁·m₂/(m₁ + m₂)");
    };
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
    // ADR 0002 #22's hard rows (their own seed, so the rows after them keep theirs): m = m at every power of two (m/2,
    // and at 2⁻¹⁰⁷⁴ the tie 2⁻¹⁰⁷⁵, to 0) and at subnormal odd multiples of 2⁻¹⁰⁷⁴ (k·2⁻¹⁰⁷⁵: ties); (6c, 2c)·2⁻¹⁰⁷⁴
    // (3c/2·2⁻¹⁰⁷⁴: ties); three families of normal ties — (2²⁷ ∓ c) (μ = (2⁵⁴ − c²)·2⁻²⁸), 3·(3·2²⁵ ∓ c) ((9·2⁵⁰ − c²)·2⁻²⁶)
    // and 3·(2²⁶ ∓ c) (3·(2⁵² − c²)·2⁻²⁷), each a 54-bit odd significand — at scales across the range, and the first
    // two moved to the subnormal grid (2⁻¹⁰⁷⁵ times an odd significand below 2⁵³); exponent gaps 45 … 64 around the proven
    // 55 (m₂ − μ < 2⁻⁵⁴·m₂ there), with the largest and smallest significands; the fast path's subnormal edge (m₂ near
    // 2⁻¹⁰²¹); m₁ + m₂ and m₁·m₂ overflowing and underflowing; binade-uniform pairs over the whole range.
    let mut r22 = Rng::new(0x5e22);
    for k in -1074..=1023 {
        rm(&mut m, pow2(k), pow2(k));
    }
    for k in (1_u64..=64).chain([(1 << 52) - 3, (1 << 52) - 1, (1 << 52) + 1, (1 << 53) - 1]) {
        let x = f64::from_bits(k);
        rm(&mut m, x, x);
    }
    for c in (1_u64..=31).step_by(2).chain([(1 << 48) + 1, (1 << 49) - 1]) {
        rm(&mut m, f64::from_bits(6 * c), f64::from_bits(2 * c));
    }
    let fam = |h: f64, base: f64, c: f64| (h * (base - c), h * (base + c));
    let fams = [
        (1.0, pow2(27), vec![1.0, 3.0, 5.0, 7.0, 99.0, 12_345.0, 1_048_577.0, 94_906_265.0]),
        (3.0, 3.0 * pow2(25), vec![1.0, 3.0, 1001.0, 33_554_431.0]),
        (3.0, pow2(26), vec![1.0, 3.0, 7.0, 12_345.0]),
    ];
    for (h, base, cs) in &fams {
        for &c in cs {
            let (a, b_) = fam(*h, *base, c);
            for t in [-1000, -500, -60, 0, 60, 500, 900, 994] {
                rm(&mut m, a * pow2(t), b_ * pow2(t));
            }
        }
    }
    for c in [94_906_267.0, pow2(26) + 1.0, pow2(27) - 1.0] {
        let (a, b_) = fam(1.0, pow2(27), c);
        rm(&mut m, a * pow2(-1047), b_ * pow2(-1047)); // (2⁵⁴ − c²)·2⁻¹⁰⁷⁵ with 2⁵⁴ − c² < 2⁵³: a subnormal tie
    }
    for c in [pow2(25) + 1.0, 3.0 * pow2(25) - 1.0] {
        let (a, b_) = fam(3.0, 3.0 * pow2(25), c);
        rm(&mut m, a * pow2(-1049), b_ * pow2(-1049));
    }
    for g in 45..=64_i32 {
        for e2 in [-1074, -1060, -1022, -1021, -1, 0, 400, 1023 - 64] {
            for (s1, s2) in [(1.0, 1.0), (1.0, 2.0 - pow2(-52)), (2.0 - pow2(-52), 1.0), (1.5, 1.25)] {
                let m2 = if e2 == -1074 { pow2(-1074) } else { s2 * pow2(e2) };
                let m1 = s1 * pow2(e2.max(-1022) + g);
                rm(&mut m, m1, m2);
                rm(&mut m, m2, m1);
            }
            let s1 = r22.uniform(1.0, 2.0);
            let s2 = r22.uniform(1.0, 2.0);
            rm(&mut m, s1 * pow2(e2.max(-1022) + g), s2 * pow2(e2.max(-1022)));
        }
    }
    for m2 in [pow2(-1021), pow2(-1021).next_down(), pow2(-1021).next_up(), pow2(-1022), pow2(-1022).next_up(),
               pow2(-1022).next_down(), pow2(-1020)] {
        for m1 in [m2, 2.0 * m2, 3.0 * m2, 1.0, pow2(-1000), pow2(-970)] {
            rm(&mut m, m1, m2);
        }
    }
    for (a, b_) in [
        (big, big),
        (big, pow2(1023)),
        (1.5 * pow2(1023), 1.5 * pow2(1023)),
        (big, 1e300),
        (1e308, 1e308),
        (big, big.next_down()),
        (big, pow2(1023 - 55)),
        (big, pow2(1023 - 54)),
        (1e300, 1e300),
        (pow2(1000), pow2(1000)),
        (pow2(-1074), pow2(-1073)),
        (pow2(-1074), 3.0 * pow2(-1074)),
        (pow2(-1074), pow2(-1022)),
        (1e-320, 1e-320),
        (1e-320, 3e-320),
        (1e-160, 1e-160),
        (pow2(-1022), pow2(-1022)),
    ] {
        rm(&mut m, a, b_);
    }
    for _ in 0..400 {
        let a = f64::from_bits(r22.int(1, 0x7FEF_FFFF_FFFF_FFFF).cast_unsigned());
        let b_ = if r22.unit() < 0.25 { a } else { f64::from_bits(r22.int(1, 0x7FEF_FFFF_FFFF_FFFF).cast_unsigned()) };
        rm(&mut m, a, b_);
    }
    for _ in 0..200 {
        let a = f64::from_bits(r22.int(1, 1 << 54).cast_unsigned());
        let b_ = f64::from_bits(r22.int(1, 1 << 54).cast_unsigned());
        rm(&mut m, a, b_); // subnormal and low-binade masses: subnormal results
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
    let vf = audit::vib;
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
    let rc = audit::rot_const;
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
    let mi = audit::moment;
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

    // ── isotope_abundances (A): every u8; ADR 0002 #23's notes ──
    for z in 0..=255u8 {
        let outs = pairs_toks(Ok(isotope_abundances(bb(z))));
        let port = pairs_toks(Ok(&new_abundances(isos, z)));
        row23(&mut m, "spectroscopy::isotope_abundances", Class::A, &[i(z)], &outs, &port);
    }

    // ── isotope_pattern (B; ADR 0002 #23's notes) ──
    pattern(&mut m, isos, &[]);
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
        pattern(&mut m, isos, atoms);
    }
    for z in 0..=255u8 {
        pattern(&mut m, isos, &[(z, 1)]);
    }
    for z in [1u8, 6, 7, 8, 16, 17, 35] {
        for c in [0u32, 1, 2, 3, 4, 5, 7, 10, 20, 50, 100, 200, 500, 1000, 2000] {
            pattern(&mut m, isos, &[(z, c)]);
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
        pattern(&mut m, isos, atoms);
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
        pattern(&mut m, isos, &atoms);
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
    audit::emit_spectroscopy(&mut m);
    emit_peaks_b6(&mut m);
    emit_fft_b6(&mut m);
    m
}

/// The overflow audit's B6 (ADR 0002 #36), appended with their own seed so the rows before keep theirs: fft_spectrum rows
/// at every power of two to 256 whose transform or magnitudes leave the normal range — huge samples (the sums and the
/// squares past DBL_MAX, the magnitude finite or not), DBL_MAX impulses, constants and alternations (the DC and Nyquist
/// sums past DBL_MAX), samples near √DBL_MAX and 2⁻⁵¹¹ (a square over- or underflows alone), subnormal samples, huge
/// and tiny samples mixed, near-cancelling samples (x and x·(1 + k·2⁻⁵²): the butterflies' differences tiny), the data
/// bound 2⁻⁹²¹ of the port's fast test ulp by ulp, a product w·d below 2⁻¹⁰²² (n = 4, w = (6.1·10⁻¹⁷, −1)), ±0 among
/// them, and NaN and ±∞ beside huge samples (Rust's IEEE spectrum kept); and seeded signals over every binade.
#[allow(clippy::too_many_lines)] // one family after another
fn emit_fft_b6(m: &mut Module) {
    let mut rng = Rng::new(0x3606);
    let (big, tiny) = (f64::MAX, f64::from_bits(1));
    let sq = big.sqrt(); // 1.34·10¹⁵⁴: a square past DBL_MAX
    let lo_sq = pow2(-511); // a square at 2⁻¹⁰²²
    let bound = pow2(-921);
    for k in 0..=8_u32 {
        let n = 1_usize << k;
        let mut imp = vec![0.0; n];
        imp[0] = big;
        spectrum(m, &imp);
        spectrum(m, &vec![big; n]);
        spectrum(m, &vec![-big * 0.75; n]);
        spectrum(m, &(0..n).map(|j| if j % 2 == 0 { big } else { -big }).collect::<Vec<_>>());
        spectrum(m, &vec![sq; n]);
        spectrum(m, &vec![sq * 0.5; n]);
        spectrum(m, &vec![lo_sq; n]);
        spectrum(m, &vec![tiny; n]);
        spectrum(m, &(0..n).map(|j| f64::from_bits(1 + j as u64)).collect::<Vec<_>>());
        for &b in &[bound, bound.next_down(), bound.next_up(), -bound] {
            let mut v = vec![0.0; n];
            v[n - 1] = b;
            spectrum(m, &v);
            spectrum(m, &vec![b; n]);
        }
        // huge, subnormal, mixed, near-cancelling, with ±0
        for _ in 0..4 {
            spectrum(m, &(0..n).map(|_| rng.signed_log_uniform(300.0, 308.25)).collect::<Vec<_>>());
            spectrum(m, &(0..n).map(|_| rng.signed_log_uniform(-323.3, -300.0)).collect::<Vec<_>>());
            spectrum(m, &(0..n).map(|_| rng.signed_log_uniform(150.0, 160.0)).collect::<Vec<_>>());
            spectrum(m, &(0..n).map(|_| rng.signed_log_uniform(-170.0, -150.0)).collect::<Vec<_>>());
            spectrum(
                m,
                &(0..n)
                    .map(|_| match rng.int(0, 3) {
                        0 => rng.signed_log_uniform(290.0, 308.25),
                        1 => rng.signed_log_uniform(-323.3, -290.0),
                        2 => 0.0,
                        _ => -0.0,
                    })
                    .collect::<Vec<_>>(),
            );
            let x = rng.signed_log_uniform(-300.0, 300.0);
            #[allow(clippy::cast_precision_loss)] // small integers
            let near: Vec<f64> = (0..n).map(|j| x * (1.0 + (j % 3) as f64 * f64::EPSILON)).collect();
            spectrum(m, &near);
        }
    }
    // a product w·d of the n = 4 stage just below 2⁻¹⁰²² (w = (cos(−π/2), −1) correctly rounded, d = x₁ − x₃): IEEE's
    // rounds it at the subnormal grid; and its magnitude below the squares' range
    let w4 = (-std::f64::consts::PI / 2.0).cos();
    let c = pow2(-1022) / w4;
    for d in -3_i64..=3 {
        let x1 = f64::from_bits(c.to_bits().wrapping_add_signed(d));
        spectrum(m, &[0.0, x1, 0.0, 0.0]);
        spectrum(m, &[x1, x1, -x1, 0.0]);
    }
    // NaN and ±∞ beside huge samples: Rust's IEEE spectrum, unchanged
    for k in 1..=4_u32 {
        let n = 1_usize << k;
        let mut v: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(300.0, 308.25)).collect();
        v[0] = f64::INFINITY;
        spectrum(m, &v);
        v[n - 1] = f64::NAN;
        spectrum(m, &v);
        v[0] = tiny;
        spectrum(m, &v);
    }
    // seeded: every binade, lengths to 64
    for _ in 0..120 {
        let n = 1_usize << rng.int(0, 6);
        let e0 = rng.uniform(-323.3, 308.0);
        let spread = rng.uniform(0.0, 40.0);
        spectrum(
            m,
            &(0..n).map(|_| rng.signed_log_uniform(e0, (e0 + spread).min(308.25))).collect::<Vec<_>>(),
        );
    }
}

/// The smallest double w > lo·(1 − t) + hi·t exactly (the exact threshold's successor on the double grid), +∞ when no
/// finite double lies above it: a bisection on the doubles' order (±0 one point), the exact predicate monotone in w.
/// (Rust's rounded threshold can be 2⁵² ulps off it: lo = −DBL_MAX, hi = DBL_MAX, t = ½ − 2⁻⁵⁴.)
fn first_above(lo: f64, hi: f64, t: f64) -> f64 {
    let key = |x: f64| -> i128 {
        let m = i128::from(x.to_bits() & 0x7FFF_FFFF_FFFF_FFFF);
        if x.is_sign_negative() { -m } else { m }
    };
    let val = |k: i128| -> f64 {
        let b = u64::try_from(k.unsigned_abs()).unwrap_or(0);
        if k < 0 { -f64::from_bits(b) } else { f64::from_bits(b) }
    };
    let (mut a, mut b) = (key(-f64::MAX), key(f64::MAX));
    if !crate::exact::peak_above(f64::MAX, lo, hi, t) {
        return f64::INFINITY;
    }
    if crate::exact::peak_above(-f64::MAX, lo, hi, t) {
        return -f64::MAX;
    }
    // invariant: val(a) is not above, val(b) is
    while b - a > 1 {
        let mid = a + (b - a) / 2;
        if crate::exact::peak_above(val(mid), lo, hi, t) { b = mid } else { a = mid }
    }
    val(b)
}

/// The overflow audit's B6 (ADR 0002 #37), appended with their own seed so the rows before keep theirs: find_peaks rows
/// with a candidate at the exact threshold's successor and its neighbours, ulp by ulp (min and max normal, subnormal,
/// huge, of mixed signs, ±DBL_MAX; t over 0, 1, 1 − 2⁻⁵³, 2⁻¹⁰⁷⁴, ½ ± an ulp and seeded), subnormal data on every t of a
/// grid, ±0 mixes, a max of DBL_MAX with a +∞ candidate, and seeded near-ties. (With −DBL_MAX and DBL_MAX, Rust's
/// threshold is 2⁵² ulps off the exact one at t = ½ ± an ulp: its 1 − t rounds by 2⁻⁵⁴, times DBL_MAX.)
#[allow(clippy::too_many_lines)] // one family after another
fn emit_peaks_b6(m: &mut Module) {
    let mut rng = Rng::new(0x3707);
    let (big, tiny) = (f64::MAX, f64::from_bits(1));
    let ts = [0.0, 1.0, 1.0_f64.next_down(), tiny, 0.5, 0.5_f64.next_down(), 0.5_f64.next_up(), 0.1, 0.75, 0.3,
        0.9, 0.25_f64.next_up(), 2.0_f64.powi(-60), 0.999];
    let ends: [(f64, f64); 14] = [
        (0.0, 1.0),
        (1.0, 3.0),
        (10.0, 11.0),
        (-5.0, 7.0),
        (-3.0, -1.0),
        (0.1, 0.7),
        (tiny, 5.0 * tiny),
        (-0.0, tiny),
        (-7.0 * tiny, 3.0 * tiny),
        (pow2(-1022), 3.0 * pow2(-1022)),
        (-big, big),
        (big * 0.5, big),
        (-big, -big * 0.25),
        (1e-300, 1e300),
    ];
    for &(lo, hi) in &ends {
        for &t in &ts {
            let w0 = first_above(lo, hi, t);
            for d in -2_i64..=2 {
                let w = f64::from_bits(w0.to_bits().wrapping_add_signed(d));
                if !(lo..=hi).contains(&w) {
                    continue;
                }
                peaks(m, &[lo, w, lo, hi, lo], t);
                peaks(m, &[hi, lo, w, lo], t);
            }
        }
    }
    // subnormal data at every t of a grid: k·2⁻¹⁰⁷⁴ around a few units
    for k in 1_u64..=6 {
        for j in 0..=16 {
            let t = f64::from(j) / 16.0;
            peaks(m, &[-0.0, f64::from_bits(k), 0.0], t);
            peaks(m, &[0.0, f64::from_bits(k), f64::from_bits(1), f64::from_bits(k + 1), -0.0], t);
            peaks(m, &[-f64::from_bits(k), -0.0, -f64::from_bits(k + 2), f64::from_bits(k), -f64::from_bits(1)], t);
        }
    }
    // ±0 mixes
    for t in [0.0, 0.5, 1.0, tiny, 1.0_f64.next_down()] {
        let mixes: [&[f64]; 4] =
            [&[-0.0, 0.0, -0.0], &[0.0, -0.0, 0.0], &[-0.0, tiny, -0.0], &[0.0, -0.0, tiny, -0.0, 0.0]];
        for d in mixes {
            peaks(m, d, t);
        }
        peaks(m, &[-0.0, 0.0, -0.0, tiny, 0.0, -0.0], t);
    }
    // a max of DBL_MAX with a +∞ candidate (above every finite threshold), ±DBL_MAX ends, a −∞ neighbour
    for t in [0.0, 1.0, 0.5, 1.0_f64.next_down(), tiny, 0.3] {
        peaks(m, &[big, f64::INFINITY, big], t);
        peaks(m, &[-big, big, -big, 0.0, f64::from_bits(1), 0.0], t);
        peaks(m, &[f64::NEG_INFINITY, -big, f64::NEG_INFINITY, big], t);
    }
    // seeded near-ties: a candidate at the exact threshold's successor or a neighbour, over every binade of the ends
    for _ in 0..400 {
        let a = rng.signed_log_uniform(-323.3, 308.2);
        let b = rng.signed_log_uniform(-323.3, 308.2);
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let t = match rng.int(0, 3) {
            0 => rng.unit(),
            1 => ts[usize::try_from(rng.int(0, 13)).unwrap_or(0)],
            2 => rng.log_uniform(-323.3, 0.0).min(1.0),
            _ => 1.0 - rng.log_uniform(-16.0, -0.3),
        };
        let w0 = first_above(lo, hi, t);
        let w = f64::from_bits(w0.to_bits().wrapping_add_signed(rng.int(-2, 1)));
        if (lo..=hi).contains(&w) {
            peaks(m, &[lo, w, lo, hi], t);
        }
    }
}
