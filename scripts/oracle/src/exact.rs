//! Exact integer arithmetic for the exact values ADR 0002 #13, #16, #19 … #22 and #24 … #32 name (M4c … the overflow
//! audit's B3): binary64 rounding of m·2^e for a mantissa of any size, the two irrational factors the port rounds once,
//! 2⁻ᶠ (N₀/2ⁿ for a fractional n) and 2^(1/6) (σ·2^(1/6)), the rational 1/(100·x) (the wavenumber conversions), binary64
//! with an unbounded exponent (the geometry, B1's kinetics, B2's Nernst equation and B3's gas laws, potentials, lattice
//! energy and binding energy, with Redlich–Kwong's Tc^2.5), the exact rationals p/q·2^e of the Bohr model, the Rydberg
//! formula and the reduced mass, the correctly rounded monomials of stoichiometry, B1, B2 and B3, B2's correctly rounded
//! √ of a quotient and B3's of a monomial (the Morse width).
//! Integer arithmetic only, independent of the port's double-word and multi-word code (`validate.py` checks every value
//! again, with `Fraction` and MPFR).

/// A finite nonzero |x| as (mant, exp) with |x| = mant·2^exp exactly, mant < 2^53.
pub fn decompose(x: f64) -> (u64, i64) {
    let bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let be = i64::try_from(bits >> 52).unwrap_or(0);
    let frac = bits & ((1_u64 << 52) - 1);
    if be == 0 { (frac, -1074) } else { (frac | (1_u64 << 52), be - 1075) }
}

/// mant·2^exp rounded once to binary64, to nearest with ties to even (subnormal results on the 2⁻¹⁰⁷⁴ grid, ±∞ past
/// DBL_MAX), negated when `neg`. Integer arithmetic only: the generator's reference for ADR 0002 #13, independent of the
/// port's double-double path (`validate.py` checks both against `Fraction`).
pub fn round_scaled(mant: u128, exp: i64, neg: bool) -> f64 {
    let sign = if neg { -1.0 } else { 1.0 };
    if mant == 0 {
        return sign * 0.0;
    }
    let top = 127 - i64::from(mant.leading_zeros()) + exp; // mant·2^exp is in [2^top, 2^(top + 1))
    if top >= 1024 {
        return sign * f64::INFINITY;
    }
    let mut lsb = if top >= -1022 { top - 52 } else { -1074 }; // the result's ulp
    let shift = lsb - exp;
    let mut q: u128 = if shift <= 0 {
        mant << u32::try_from(-shift).unwrap_or(0)
    } else if shift >= 127 {
        0
    } else {
        let s = u32::try_from(shift).unwrap_or(0);
        let q0 = mant >> s;
        let rem = mant & ((1_u128 << s) - 1);
        let half = 1_u128 << (s - 1);
        if rem > half || (rem == half && q0 & 1 == 1) { q0 + 1 } else { q0 }
    };
    if q == 1_u128 << 53 {
        q >>= 1;
        lsb += 1;
    }
    if q == 0 {
        return sign * 0.0;
    }
    let bits = if q < 1_u128 << 52 {
        u64::try_from(q).unwrap_or(0) // subnormal: q units of 2⁻¹⁰⁷⁴ (lsb is −1074)
    } else {
        let e = lsb + 52;
        if e >= 1024 {
            return sign * f64::INFINITY;
        }
        (u64::try_from(e + 1023).unwrap_or(0) << 52) | (u64::try_from(q).unwrap_or(0) & ((1_u64 << 52) - 1))
    };
    sign * f64::from_bits(bits)
}

/// A non-negative integer of any size, little-endian 64-bit limbs: the exact arithmetic of ADR 0002 #13's reference for
/// a fractional n (M4d), where N₀/2ⁿ needs 2⁻ᶠ well past 53 bits. Only what that needs: shifts, +, −, ×, ÷ by a word.
#[derive(Clone)]
pub struct Big(Vec<u64>);

impl Big {
    fn new(x: u128) -> Self {
        #[allow(clippy::cast_possible_truncation)]
        Self(vec![x as u64, (x >> 64) as u64]).norm()
    }

    fn norm(mut self) -> Self {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
        self
    }

    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The bit length (0 for zero).
    fn bits(&self) -> u64 {
        self.0.last().map_or(0, |&t| 64 * (self.0.len() as u64 - 1) + u64::from(64 - t.leading_zeros()))
    }

    fn shl(&self, n: u64) -> Self {
        let (w, b) = (usize::try_from(n / 64).unwrap_or(0), u32::try_from(n % 64).unwrap_or(0));
        let mut out = vec![0_u64; w + self.0.len() + 1];
        for (k, &x) in self.0.iter().enumerate() {
            out[w + k] |= x << b;
            if b > 0 {
                out[w + k + 1] |= x >> (64 - b);
            }
        }
        Self(out).norm()
    }

    fn shr(&self, n: u64) -> Self {
        let (w, b) = (usize::try_from(n / 64).unwrap_or(usize::MAX), u32::try_from(n % 64).unwrap_or(0));
        if w >= self.0.len() {
            return Self(Vec::new());
        }
        let mut out = vec![0_u64; self.0.len() - w];
        for (k, o) in out.iter_mut().enumerate() {
            *o = self.0[k + w] >> b;
            if b > 0 && k + w + 1 < self.0.len() {
                *o |= self.0[k + w + 1] << (64 - b);
            }
        }
        Self(out).norm()
    }

    fn add(&self, o: &Self) -> Self {
        let n = self.0.len().max(o.0.len());
        let mut out = Vec::with_capacity(n + 1);
        let mut c = 0_u128;
        for k in 0..n {
            let t = u128::from(self.0.get(k).copied().unwrap_or(0)) + u128::from(o.0.get(k).copied().unwrap_or(0)) + c;
            #[allow(clippy::cast_possible_truncation)]
            out.push(t as u64);
            c = t >> 64;
        }
        #[allow(clippy::cast_possible_truncation)]
        out.push(c as u64);
        Self(out).norm()
    }

    /// self − o, for self ≥ o.
    fn sub(&self, o: &Self) -> Self {
        assert!(self.cmp(o) != std::cmp::Ordering::Less, "Big::sub would be negative");
        let mut out = Vec::with_capacity(self.0.len());
        let mut borrow = false;
        for (k, &a) in self.0.iter().enumerate() {
            let (d1, b1) = a.overflowing_sub(o.0.get(k).copied().unwrap_or(0));
            let (d2, b2) = d1.overflowing_sub(u64::from(borrow));
            out.push(d2);
            borrow = b1 || b2;
        }
        Self(out).norm()
    }

    fn mul(&self, o: &Self) -> Self {
        let mut out = vec![0_u64; self.0.len() + o.0.len() + 1];
        for (i, &a) in self.0.iter().enumerate() {
            let mut c = 0_u128;
            for (j, &b_) in o.0.iter().enumerate() {
                let t = u128::from(a) * u128::from(b_) + u128::from(out[i + j]) + c;
                #[allow(clippy::cast_possible_truncation)]
                let lo = t as u64;
                out[i + j] = lo;
                c = t >> 64;
            }
            let mut k = i + o.0.len();
            while c > 0 {
                let t = u128::from(out[k]) + c;
                #[allow(clippy::cast_possible_truncation)]
                let lo = t as u64;
                out[k] = lo;
                c = t >> 64;
                k += 1;
            }
        }
        Self(out).norm()
    }

    /// ⌊self / d⌋ for a word d > 0.
    fn div_word(&self, d: u64) -> Self {
        let mut out = vec![0_u64; self.0.len()];
        let mut r = 0_u128;
        for k in (0..self.0.len()).rev() {
            let cur = (r << 64) | u128::from(self.0[k]);
            #[allow(clippy::cast_possible_truncation)]
            let q = (cur / u128::from(d)) as u64;
            out[k] = q;
            r = cur % u128::from(d);
        }
        Self(out).norm()
    }

    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.0.len().cmp(&o.0.len()).then_with(|| self.0.iter().rev().cmp(o.0.iter().rev()))
    }
}

/// The fixed-point precision of the fractional-n reference: 2⁻ᶠ is computed to 2⁻³²⁰.
const PREC: u64 = 320;

/// m·2^exp (a non-negative `Big`) rounded once to binary64, to nearest with ties to even (subnormal results on the
/// 2⁻¹⁰⁷⁴ grid, ±∞ past DBL_MAX), negated when `neg`: [`round_scaled`] for a mantissa of any size.
fn round_big(m: &Big, exp: i64, neg: bool) -> f64 {
    if m.bits() <= 127 {
        let lo = u128::from(m.0.first().copied().unwrap_or(0)) | (u128::from(m.0.get(1).copied().unwrap_or(0)) << 64);
        return round_scaled(lo, exp, neg);
    }
    // Drop all but the top 120 bits, keeping a sticky bit below them: the rounding to 53 bits (or fewer, at the
    // subnormal grid) is unchanged, since a nonzero dropped part only breaks a tie that the sticky bit already breaks.
    let drop = m.bits() - 120;
    let top = m.shr(drop);
    let sticky = u128::from(top.shl(drop).cmp(m) != std::cmp::Ordering::Equal);
    let t = u128::from(top.0.first().copied().unwrap_or(0)) | (u128::from(top.0.get(1).copied().unwrap_or(0)) << 64);
    round_scaled((t << 1) | sticky, exp + i64::try_from(drop).unwrap_or(0) - 1, neg)
}

/// ln 2·2^PREC to within 2 units: Σ 2^(PREC + 16 − k)/k for k ≥ 1, each term floored, with 16 guard bits and the tail
/// past k = PREC + 32 (below 2^(−PREC − 32)) dropped.
fn ln2_fixed() -> Big {
    let one = Big::new(1).shl(PREC + 16);
    let mut sum = Big::new(0);
    for k in 1..=(PREC + 32) {
        sum = sum.add(&one.shr(k).div_word(k));
    }
    sum.shr(16)
}

/// e^(−x)·2^PREC for a fixed-point x = x_fix·2^(−PREC) in [0, ln 2]: the alternating Taylor series, each term
/// ⌊⌊t·x⌋/n⌋ (its error below 3 units: δtₙ ≤ x·δtₙ₋₁/n + 2; the terms fall below one unit by n = 70), so within 3 units a
/// term, 240 in all.
fn exp_neg_fixed(x_fix: &Big) -> Big {
    let mut term = Big::new(1).shl(PREC);
    let (mut pos, mut neg) = (term.clone(), Big::new(0));
    for n in 1..200_u64 {
        term = term.mul(x_fix).shr(PREC).div_word(n);
        if term.is_zero() {
            break;
        }
        if n % 2 == 1 { neg = neg.add(&term) } else { pos = pos.add(&term) }
    }
    pos.sub(&neg)
}

/// ADR 0002 #13 for a fractional n (M4d): mant·2^exp·2⁻ᶠ (f a double in (0, 1)) correctly rounded, negated when `neg`.
/// 2⁻ᶠ = e^(−f·ln 2) to within 2^12 units of 2⁻³²⁰ (ln 2 to 2 units, f·ln 2 to 3, which moves e^(−x) by at most 3,
/// the series' 240), so the value lies
/// between mant·(V − 2^12) and mant·(V + 2^12) (V the computed 2⁻ᶠ·2^PREC); both bounds rounded, and they must agree
/// (Ziv's test; 2⁻ᶠ is irrational for f in (0, 1), so no value is a midpoint, and a disagreement would mean one within
/// 2⁻³⁰⁶ of it). Integer arithmetic only, independent of the port's exp machinery (`validate.py` checks it against
/// MPFR).
pub fn cr_exp2_frac(mant: u64, exp: i64, f: f64, neg: bool) -> f64 {
    let (fm, fe) = decompose(f); // f = fm·2^fe exactly, fe < 0
    let x_fix = ln2_fixed().mul(&Big::new(u128::from(fm))).shr(u64::try_from(-fe).unwrap_or(0));
    let v = exp_neg_fixed(&x_fix);
    let slack = Big::new(1 << 12);
    let m = Big::new(u128::from(mant));
    let p = i64::try_from(PREC).unwrap_or(0);
    let lo = round_big(&m.mul(&v.sub(&slack)), exp - p, neg);
    let hi = round_big(&m.mul(&v.add(&slack)), exp - p, neg);
    assert!(lo.to_bits() == hi.to_bits(), "ADR 0002 #13: N₀/2ⁿ undecided at 2⁻³⁰⁶ (mant {mant}, exp {exp}, f {f:?})");
    lo
}


/// ADR 0002 #16 (M4d): RN(x·2^(1/6)) for a finite x ≠ 0, decided exactly: 2^(1/6) is irrational with (2^(1/6))⁶ = 2, so
/// x·2^(1/6) > t ⟺ 2·x⁶ > t⁶ for t > 0, an integer comparison. Starting from Rust's x·RN(2^(1/6)) (within an ulp),
/// the candidate moves by an ulp until both midpoints around it bracket |x|·2^(1/6); no value is a midpoint.
pub fn cr_root6_2_mul(x: f64, rust: f64) -> f64 {
    assert!(x.is_finite() && x != 0.0, "cr_root6_2_mul({x:?})");
    let neg = x < 0.0;
    let (sm, se) = decompose(x);
    // 2·|x|⁶ as an integer times 2^(6·se + 1)
    let sm_b = Big::new(u128::from(sm));
    let x6 = sm_b.mul(&sm_b).mul(&sm_b).mul(&sm_b).mul(&sm_b).mul(&sm_b);
    let x6e = 6 * se + 1;
    // |x|·2^(1/6) > m·2^me (m > 0) ⟺ x6·2^x6e > m⁶·2^(6·me)
    let above = |m: u128, me: i64| -> bool {
        let mb = Big::new(m);
        let m6 = mb.mul(&mb).mul(&mb).mul(&mb).mul(&mb).mul(&mb);
        let me6 = 6 * me;
        let (l, r) = if x6e >= me6 {
            (x6.shl(u64::try_from(x6e - me6).unwrap_or(0)), m6)
        } else {
            (x6.clone(), m6.shl(u64::try_from(me6 - x6e).unwrap_or(0)))
        };
        l.cmp(&r) == std::cmp::Ordering::Greater
    };
    let mut v = rust.abs();
    if v.is_infinite() {
        v = f64::MAX;
    }
    loop {
        let (vm, ve) = decompose(v);
        let vm = u128::from(vm);
        // the midpoints around v: above at (2·vm + 1)·2^(ve − 1); below at (2·vm − 1)·2^(ve − 1), or (4·vm − 1)·2^(ve − 2)
        // when v is a normal power of two (the binade below has half the ulp)
        let up = above(2 * vm + 1, ve - 1);
        let down_ok = if vm == 1 << 52 && ve > -1074 { above(4 * vm - 1, ve - 2) } else { above(2 * vm - 1, ve - 1) };
        if up {
            if v == f64::MAX {
                return if neg { f64::NEG_INFINITY } else { f64::INFINITY }; // past DBL_MAX + ulp/2: overflow
            }
            v = v.next_up();
        } else if !down_ok {
            v = v.next_down();
            assert!(v > 0.0, "cr_root6_2_mul({x:?}) below the smallest subnormal");
        } else {
            return if neg { -v } else { v };
        }
    }
}

/// ADR 0002 #19 (M4e): RN(1/(100·x)), the exact (1/100)/x rounded once, for any x: ±0 gives ±∞, ±∞ ±0, NaN NaN (Rust's
/// 0.01/x there). |x| = m·2^e; q = ⌊2^120/(100·m)⌋ (at least 2^60, so at least 8 bits below the result's last) and the
/// remainder as a sticky bit, (2q + sticky)·2^(−121 − e) rounded by `round_scaled`. Integer arithmetic only, independent of
/// the port's double-double fast path (`validate.py` checks every value again with `Fraction`).
pub fn cr_hundredth_div(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return if x > 0.0 { 0.0 } else { -0.0 };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { f64::NEG_INFINITY } else { f64::INFINITY };
    }
    let (m, e) = decompose(x);
    let d = u128::from(m) * 100;
    let num = 1_u128 << 120;
    let q = num / d;
    let sticky = u128::from(!num.is_multiple_of(d));
    round_scaled(2 * q + sticky, -121 - e, x < 0.0)
}

// ── ADR 0002 #20 (M4e-2): binary64 with an unbounded exponent ──────────────────────────────────────────────

/// A binary64 value with an unbounded exponent range: ±m·2^e, with m = 0 (a signed zero) or 2^52 ≤ m < 2^53. The
/// generator's reference for ADR 0002 #20 (potential's geometry): Rust's operation sequence with every operation's exact
/// result rounded to 53 bits, nearest-even, and no overflow or underflow. Integer arithmetic only, independent of the
/// port's (significand, exponent) doubles (`validate.py` checks every value again with MPFR at an unbounded exponent).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ux {
    neg: bool,
    m: u64,
    e: i64,
}

/// ⌊√n⌋ (Newton's iteration from above).
fn isqrt(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x: u128 = 1 << (128 - n.leading_zeros()).div_ceil(2);
    loop {
        let y = (x + n / x) / 2;
        if y >= x {
            return x;
        }
        x = y;
    }
}

impl Ux {
    /// The exact value of a finite double.
    pub fn from_f64(x: f64) -> Self {
        assert!(x.is_finite(), "Ux::from_f64({x:?})");
        let neg = x.is_sign_negative();
        if x == 0.0 {
            return Self { neg, m: 0, e: 0 };
        }
        let (m, e) = decompose(x);
        Self::round(neg, u128::from(m), e)
    }

    pub fn is_zero(self) -> bool {
        self.m == 0
    }

    /// ±mant·2^exp rounded to 53 bits, to nearest with ties to even, at any exponent. A caller with a nonzero part below
    /// `mant` passes it as a sticky low bit, with at least two bits between it and the 53rd.
    fn round(neg: bool, mant: u128, exp: i64) -> Self {
        if mant == 0 {
            return Self { neg, m: 0, e: 0 };
        }
        let b = 128 - mant.leading_zeros();
        if b <= 53 {
            let s = 53 - b;
            return Self { neg, m: u64::try_from(mant << s).unwrap_or(0), e: exp - i64::from(s) };
        }
        let s = b - 53;
        let mut q = mant >> s;
        let rem = mant & ((1_u128 << s) - 1);
        let half = 1_u128 << (s - 1);
        if rem > half || (rem == half && q & 1 == 1) {
            q += 1;
        }
        let mut e = exp + i64::from(s);
        if q == 1 << 53 {
            q >>= 1;
            e += 1;
        }
        Self { neg, m: u64::try_from(q).unwrap_or(0), e }
    }

    pub fn neg(self) -> Self {
        Self { neg: !self.neg, ..self }
    }

    /// RN53(x·y).
    pub fn mul(self, o: Self) -> Self {
        let neg = self.neg != o.neg;
        if self.is_zero() || o.is_zero() {
            return Self { neg, m: 0, e: 0 };
        }
        Self::round(neg, u128::from(self.m) * u128::from(o.m), self.e + o.e)
    }

    /// RN53(x + y), with IEEE's signed zeros (an exact cancellation is +0).
    pub fn add(self, o: Self) -> Self {
        if self.is_zero() && o.is_zero() {
            return Self { neg: self.neg && o.neg, m: 0, e: 0 };
        }
        if self.is_zero() {
            return o;
        }
        if o.is_zero() {
            return self;
        }
        let (h, l) = if self.e >= o.e { (self, o) } else { (o, self) };
        let d = h.e - l.e;
        if d > 72 {
            // l is below 2^-19 of h's last place: a sticky bit three places below it decides every rounding
            let base = u128::from(h.m) << 3;
            let mant = if h.neg == l.neg { base | 1 } else { base - 1 };
            return Self::round(h.neg, mant, h.e - 3);
        }
        let hm = u128::from(h.m) << u32::try_from(d).unwrap_or(0);
        let lm = u128::from(l.m);
        if h.neg == l.neg {
            return Self::round(h.neg, hm + lm, l.e);
        }
        match hm.cmp(&lm) {
            std::cmp::Ordering::Equal => Self { neg: false, m: 0, e: 0 },
            std::cmp::Ordering::Greater => Self::round(h.neg, hm - lm, l.e),
            std::cmp::Ordering::Less => Self::round(l.neg, lm - hm, l.e),
        }
    }

    pub fn sub(self, o: Self) -> Self {
        self.add(o.neg())
    }

    /// (M, e') with √(m·2^e) = √M·2^e', M ≥ 2^112 (so ⌊√M⌋ has at least 56 bits).
    fn sqrt_parts(self) -> (u128, i64) {
        assert!(!self.neg || self.is_zero(), "Ux::sqrt of a negative value");
        let k: u32 = if (self.e - 60).rem_euclid(2) == 0 { 60 } else { 61 };
        let mm = u128::from(self.m) << k;
        (mm, (self.e - i64::from(k)) / 2)
    }

    /// RN53(√x) for x ≥ 0 (√±0 = ±0).
    pub fn sqrt(self) -> Self {
        if self.is_zero() {
            return self;
        }
        let (mm, e) = self.sqrt_parts();
        let r = isqrt(mm);
        let sticky = u128::from(r * r != mm);
        Self::round(false, (r << 1) | sticky, e - 1)
    }

    /// The exact order of two values (zeros equal whatever their signs).
    pub fn cmp(self, o: Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let (sz, oz) = (self.is_zero(), o.is_zero());
        if sz && oz {
            return Ordering::Equal;
        }
        let sgn = |v: Self, z: bool| if z { 0 } else if v.neg { -1 } else { 1 };
        let (ss, os) = (sgn(self, sz), sgn(o, oz));
        if ss != os {
            return ss.cmp(&os);
        }
        let mag = (self.e, self.m).cmp(&(o.e, o.m));
        if ss < 0 { mag.reverse() } else { mag }
    }

    pub fn le(self, o: Self) -> bool {
        self.cmp(o) != std::cmp::Ordering::Greater
    }

    /// Rust's `x.max(y)` on the exact values: x where x > y, y otherwise.
    pub fn max(self, o: Self) -> Self {
        if self.cmp(o) == std::cmp::Ordering::Greater { self } else { o }
    }

    pub fn is_positive(self) -> bool {
        !self.neg && !self.is_zero()
    }

    /// √x rounded once into binary64, for x ≥ 0: ⌊√M⌋ and a sticky bit, then `round_scaled` (which drops at least
    /// five bits, so the sticky bit decides every midpoint).
    pub fn sqrt_f64(self) -> f64 {
        if self.is_zero() {
            return if self.neg { -0.0 } else { 0.0 };
        }
        let (mm, e) = self.sqrt_parts();
        let r = isqrt(mm);
        let sticky = u128::from(r * r != mm);
        round_scaled((r << 1) | sticky, e - 1, false)
    }

    /// (2Q + sticky, exponent) for |x/y| = (Q + f)·2^(e − 64), Q = ⌊m₁·2⁶⁴/m₂⌋ (at least 64 bits: m₁ ≥ 2⁵², m₂ < 2⁵³) and
    /// sticky = [f > 0]; x nonzero, y nonzero.
    fn div_parts(self, o: Self) -> (u128, i64) {
        let num = u128::from(self.m) << 64;
        let d = u128::from(o.m);
        ((num / d) << 1 | u128::from(num % d != 0), self.e - o.e - 65)
    }

    /// RN53(x/y) for y ≠ 0 (ADR 0002 #28, the overflow audit's B1: kinetics' 1/[A]₀ and 1/(…), Vmax·[S]/(Km + [S])), a
    /// zero x giving IEEE's signed zero.
    pub fn div(self, o: Self) -> Self {
        assert!(!o.is_zero(), "Ux::div by zero");
        let neg = self.neg != o.neg;
        if self.is_zero() {
            return Self { neg, m: 0, e: 0 };
        }
        let (mant, e) = self.div_parts(o);
        Self::round(neg, mant, e)
    }

    /// x/y rounded once into binary64 (its subnormal grid, ±∞ past DBL_MAX), y ≠ 0: `div_parts` and `round_scaled`
    /// (which drops at least twelve bits, so the sticky bit decides every midpoint).
    pub fn div_f64(self, o: Self) -> f64 {
        assert!(!o.is_zero(), "Ux::div_f64 by zero");
        let neg = self.neg != o.neg;
        if self.is_zero() {
            return if neg { -0.0 } else { 0.0 };
        }
        let (mant, e) = self.div_parts(o);
        round_scaled(mant, e, neg)
    }

    /// x + y rounded once into binary64 (its subnormal grid, ±∞ past DBL_MAX), with IEEE's signed zeros (an exact
    /// cancellation is +0; −0 + −0 is −0): ADR 0002 #29 (the overflow audit's B2), the Nernst equation's final
    /// subtraction. The exact sum for an exponent gap up to 72 (below 2¹²⁶), otherwise the larger term with a sticky bit
    /// three places below its last (the smaller is below 2⁻¹⁹ of that place, and no binary64 rounding boundary lies
    /// strictly between the two: `add`'s argument, and at the subnormal grid the larger term's last place is at most
    /// 2⁻¹⁰⁷⁵, so its boundaries are multiples of it); then `round_scaled`.
    pub fn add_f64(self, o: Self) -> f64 {
        if self.is_zero() && o.is_zero() {
            return if self.neg && o.neg { -0.0 } else { 0.0 };
        }
        if self.is_zero() {
            return round_scaled(u128::from(o.m), o.e, o.neg);
        }
        if o.is_zero() {
            return round_scaled(u128::from(self.m), self.e, self.neg);
        }
        let (h, l) = if self.e >= o.e { (self, o) } else { (o, self) };
        let d = h.e - l.e;
        if d > 72 {
            let base = u128::from(h.m) << 3;
            let mant = if h.neg == l.neg { base | 1 } else { base - 1 };
            return round_scaled(mant, h.e - 3, h.neg);
        }
        let hm = u128::from(h.m) << u32::try_from(d).unwrap_or(0);
        let lm = u128::from(l.m);
        if h.neg == l.neg {
            return round_scaled(hm + lm, l.e, h.neg);
        }
        match hm.cmp(&lm) {
            std::cmp::Ordering::Equal => 0.0,
            std::cmp::Ordering::Greater => round_scaled(hm - lm, l.e, h.neg),
            std::cmp::Ordering::Less => round_scaled(lm - hm, l.e, l.neg),
        }
    }

    /// x − y rounded once into binary64 (`add_f64` of −y).
    pub fn sub_f64(self, o: Self) -> f64 {
        self.add_f64(o.neg())
    }

    /// x·y rounded once into binary64 (its subnormal grid, ±∞ past DBL_MAX), with IEEE's signed zeros: ADR 0002 #32 (the
    /// overflow audit's B3), the Lennard-Jones and Born–Landé final products. The exact product of the two 53-bit
    /// mantissas (below 2¹⁰⁶), then `round_scaled`.
    pub fn mul_f64(self, o: Self) -> f64 {
        let neg = self.neg != o.neg;
        if self.is_zero() || o.is_zero() {
            return if neg { -0.0 } else { 0.0 };
        }
        round_scaled(u128::from(self.m) * u128::from(o.m), self.e + o.e, neg)
    }

    /// RN53(x^2.5) for x > 0, at any exponent: ADR 0002 #31 (the overflow audit's B3), Redlich–Kwong's Tc^2.5 with an
    /// unbounded exponent. x^2.5 = √(x⁵) = √(m⁵·2^(5e)) exactly, by [`sqrt_parts_big`] (⌊√⌋ of a quotient of at least 121
    /// bits and a sticky bit, at least eight bits below the 53rd) and `round`. Integer arithmetic only, independent of the
    /// port's rescaled `_km_cr_pow(s, 2.5)` (`validate.py` checks every value again with MPFR's pow at 53 bits and an
    /// unbounded exponent).
    pub fn pow25(self) -> Self {
        assert!(!self.neg && !self.is_zero(), "Ux::pow25 of a value that is not positive");
        let m = Big::new(u128::from(self.m));
        let m5 = m.mul(&m).mul(&m).mul(&m).mul(&m);
        let (mant, h) = sqrt_parts_big(&m5, &Big::new(1), 5 * self.e);
        Self::round(false, mant, h)
    }

    /// Whether this is exactly the double x (a finite x; zeros by sign too).
    pub fn is(self, x: f64) -> bool {
        x.is_finite() && Self::from_f64(x) == self
    }
}

impl Big {
    fn bit(&self, k: u64) -> bool {
        let w = usize::try_from(k / 64).unwrap_or(usize::MAX);
        self.0.get(w).is_some_and(|&x| (x >> (k % 64)) & 1 == 1)
    }

    /// (⌊self / d⌋, self mod d) for d > 0: binary long division.
    fn divrem(&self, d: &Self) -> (Self, Self) {
        assert!(!d.is_zero(), "Big::divrem by zero");
        let mut q = vec![0_u64; self.0.len()];
        let mut r = Self(Vec::new());
        for k in (0..self.bits()).rev() {
            r = r.shl(1);
            if self.bit(k) {
                r = r.add(&Self::new(1));
            }
            if r.cmp(d) != std::cmp::Ordering::Less {
                r = r.sub(d);
                let w = usize::try_from(k / 64).unwrap_or(0);
                q[w] |= 1 << (k % 64);
            }
        }
        (Self(q).norm(), r)
    }
}

/// The fixed-point precision of the atan2 reference: θ to 2⁻⁴⁰⁰.
const APREC: u64 = 400;

/// atan(z)·2^APREC for z = zf·2^-APREC in [0, 1], within 2^11 units: Euler's series atan z = Σₖ cₖ·u·ρᵏ with
/// u = z/(1 + z²), ρ = z²/(1 + z²) ≤ ½, c₀ = 1, cₖ = cₖ₋₁·2k/(2k + 1), each from zf exactly (one division each) and
/// floored; every term's error stays below 4 units, the terms fall below one unit within about 420, and zf's own unit
/// moves the sum by at most 6.
fn atan_fixed(zf: &Big) -> Big {
    let one2 = Big::new(1).shl(2 * APREC);
    let den = one2.add(&zf.mul(zf));
    let (rho, _) = zf.mul(zf).shl(APREC).divrem(&den);
    let (u, _) = zf.shl(2 * APREC).divrem(&den);
    let mut term = u.clone();
    let mut sum = u;
    let mut k = 1_u64;
    while !term.is_zero() {
        term = term.mul(&rho).shr(APREC).mul(&Big::new(u128::from(2 * k))).div_word(2 * k + 1);
        sum = sum.add(&term);
        k += 1;
    }
    sum
}

/// atan2(y, x), correctly rounded in binary64, of two unbounded-exponent values (the angle of (x, y) in [−π, π]).
/// A zero argument has C99's special values (glibc's atan2 of ±0 and ±1 gives them exactly). Otherwise the octant
/// reduction: θ = atan(n/d) with n ≤ d, the result θ, π − θ, π/2 − θ or π/2 + θ with y's sign, in fixed point to 2⁻⁴⁰⁰
/// with a proven error below 2^14 units; both ends of ± 2^20 must round to the same double (Ziv), which they do unless
/// the result is within 2⁻³⁷⁰ of a midpoint (atan of a nonzero rational is transcendental, and so are π ± it). For
/// n/d < 2⁻⁶⁰ with x > 0 the result is atan(z) ∈ (z − z³/3, z), bracketed from z itself at 200 bits (z is never within
/// 2⁻¹⁰⁷ of a midpoint unless it is one, when atan(z) lies just below it).
pub fn ux_atan2_f64(y: Ux, x: Ux) -> f64 {
    if y.is_zero() || x.is_zero() {
        let py = if y.neg { -1.0 } else { 1.0 } * if y.is_zero() { 0.0 } else { 1.0 };
        let px = if x.neg { -1.0 } else { 1.0 } * if x.is_zero() { 0.0 } else { 1.0 };
        return f64::atan2(py, px);
    }
    let (ym, xm) = ((y.e, y.m), (x.e, x.m));
    let swap = ym > xm;
    let (n, d) = if swap { (x, y) } else { (y, x) };
    let oct = u32::from(x.neg) + 2 * u32::from(swap);
    let k = n.e - d.e;
    let (mn, md) = (Big::new(u128::from(n.m)), Big::new(u128::from(d.m)));
    let p = i64::try_from(APREC).unwrap_or(0);
    if oct == 0 && k <= -61 {
        // atan(z) ∈ (z − z³/3, z), z = n/d < 2^-60: Z = ⌊mn·2^200/md⌋, so z·2^(200 − k) ∈ [Z, Z + 1); in units of
        // 2^(k − 202): atan(z) ∈ (4Z − 2^83, 4Z) (z³/3 < z·2^-119), and below 4Z + 4 when Z is not exact
        let (zq, rem) = mn.shl(200).divrem(&md);
        let base = zq.shl(2);
        let lo = base.sub(&Big::new(1 << 83));
        let hi = if rem.is_zero() { base.sub(&Big::new(1)) } else { base.add(&Big::new(4)) };
        let (a, b) = (round_big(&lo, k - 202, y.neg), round_big(&hi, k - 202, y.neg));
        assert!(a.to_bits() == b.to_bits(), "ADR 0002 #20: tiny atan2 undecided ({y:?}, {x:?})");
        return a;
    }
    let zf = if p + k >= 0 { mn.shl(u64::try_from(p + k).unwrap_or(0)).divrem(&md).0 } else { Big::new(0) };
    let theta = atan_fixed(&zf);
    static PI4: std::sync::OnceLock<Big> = std::sync::OnceLock::new();
    let pi4 = PI4.get_or_init(|| atan_fixed(&Big::new(1).shl(APREC)));
    let v = match oct {
        0 => theta,
        1 => pi4.shl(2).sub(&theta),
        2 => pi4.shl(1).sub(&theta),
        _ => pi4.shl(1).add(&theta),
    };
    let slack = Big::new(1 << 20);
    let lo = if v.cmp(&slack) == std::cmp::Ordering::Greater { v.sub(&slack) } else { Big::new(0) };
    let hi = v.add(&slack);
    let (a, b) = (round_big(&lo, -p, y.neg), round_big(&hi, -p, y.neg));
    assert!(a.to_bits() == b.to_bits(), "ADR 0002 #20: atan2 undecided at 2^-380 ({y:?}, {x:?})");
    a
}

// ── ADR 0002 #21, #22 (M4e-3): exact rationals rounded once ─────────────────────────────────────────────────

/// p/q·2^e rounded once to binary64 (nearest, ties to even; subnormal results on the 2⁻¹⁰⁷⁴ grid, ±∞ past DBL_MAX), for
/// integers p ≥ 0 and q ≥ 1 of any size: Q = ⌊p·2^k/q⌋ with at least 60 bits and the remainder as a sticky bit
/// (`round_big` drops at least seven bits below Q's, so the sticky bit decides every midpoint). Integer arithmetic only,
/// independent of the port's multi-word long division (`validate.py` checks every value again with `Fraction`).
fn cr_ratio(p: &Big, q: &Big, e: i64) -> f64 {
    assert!(!q.is_zero(), "cr_ratio: q = 0");
    if p.is_zero() {
        return 0.0;
    }
    let k = (60 + i64::try_from(q.bits()).unwrap_or(0) - i64::try_from(p.bits()).unwrap_or(0)).max(0);
    let (quo, rem) = p.shl(u64::try_from(k).unwrap_or(0)).divrem(q);
    let mant = quo.shl(1).add(&Big::new(u128::from(!rem.is_zero())));
    round_big(&mant, e - k - 1, false)
}

/// ADR 0002 #21: RN(r·Z²/n²) for a positive finite r (the Bohr level's magnitude), Z and n ≥ 1.
pub fn cr_bohr_level(r: f64, z: u8, n: u32) -> f64 {
    let (rm, re) = decompose(r);
    let p = Big::new(u128::from(rm) * u128::from(z) * u128::from(z));
    let q = Big::new(u128::from(n) * u128::from(n));
    cr_ratio(&p, &q, re)
}

/// The exact Rydberg integers for n₁ < n₂ and r = M·2^E: (M·Z²·(n₂² − n₁²), n₁²·n₂², E).
fn rydberg_ints(r: f64, z: u8, n_lower: u32, n_upper: u32) -> (Big, Big, i64) {
    assert!(n_lower >= 1 && n_lower < n_upper, "rydberg_ints({n_lower}, {n_upper})");
    let (rm, re) = decompose(r);
    let (l, u) = (u128::from(n_lower), u128::from(n_upper));
    let num = Big::new(u128::from(rm) * u128::from(z) * u128::from(z)).mul(&Big::new(u * u - l * l));
    let w = Big::new(l * u);
    (num, w.mul(&w), re)
}

/// ADR 0002 #21: RN(r·Z²·(n₂² − n₁²)/(n₁²·n₂²)) (the transition energy with r = RYDBERG_EV).
pub fn cr_transition(r: f64, z: u8, n_lower: u32, n_upper: u32) -> f64 {
    let (num, den, re) = rydberg_ints(r, z, n_lower, n_upper);
    cr_ratio(&num, &den, re)
}

/// ADR 0002 #21: RN(n₁²·n₂²/(r·Z²·(n₂² − n₁²))) (the Rydberg wavelength with r = RYDBERG or R_H).
pub fn cr_rydberg(r: f64, z: u8, n_lower: u32, n_upper: u32) -> f64 {
    let (num, den, re) = rydberg_ints(r, z, n_lower, n_upper);
    cr_ratio(&den, &num, -re)
}

/// ADR 0002 #22: RN(m₁·m₂/(m₁ + m₂)) for positive finite doubles, exactly at any exponent gap: m = M·2^E, and with
/// E₁ ≥ E₂, μ = M₁·M₂·2^E₁/(M₁·2^(E₁ − E₂) + M₂) (no cutoff: the port's proven gap of 55 is checked, not assumed).
pub fn cr_reduced_mass(m1: f64, m2: f64) -> f64 {
    assert!(m1 > 0.0 && m1.is_finite() && m2 > 0.0 && m2.is_finite(), "cr_reduced_mass({m1:?}, {m2:?})");
    let ((a, ea), (b, eb)) = (decompose(m1), decompose(m2));
    let ((ma, e1), (mb, e2)) = if ea >= eb { ((a, ea), (b, eb)) } else { ((b, eb), (a, ea)) };
    let p = Big::new(u128::from(ma) * u128::from(mb));
    let q = Big::new(u128::from(ma)).shl(u64::try_from(e1 - e2).unwrap_or(0)).add(&Big::new(u128::from(mb)));
    cr_ratio(&p, &q, e1)
}

// ── ADR 0002 #24, #25, #26 (M4f-2): stoichiometry's exact rationals and comparisons ─────────────────────────────

/// ADR 0002 #24: RN(x·y/z) for positive finite doubles x, y, z, at any exponents: x·y/z = (Mₓ·M_y)/M_z·2^(Eₓ + E_y − E_z)
/// rounded once by [`cr_ratio`] (percent_yield's 100·a/t, theoretical_yield_moles' m·p/c, percent_composition's
/// mass·(100·count)/mw). Integer arithmetic only, independent of the port's double-double fast path and its multi-word
/// division (`validate.py` checks every value again with `Fraction`).
pub fn cr_mul_div(x: f64, y: f64, z: f64) -> f64 {
    assert!(
        [x, y, z].iter().all(|v| v.is_finite() && *v > 0.0),
        "cr_mul_div({x:?}, {y:?}, {z:?})"
    );
    let ((mx, ex), (my, ey), (mz, ez)) = (decompose(x), decompose(y), decompose(z));
    cr_ratio(&Big::new(u128::from(mx) * u128::from(my)), &Big::new(u128::from(mz)), ex + ey - ez)
}

/// The overflow audit's B1 (ADR 0002 #27, #28): RN(Π num / Π den) for finite doubles with nonzero denominator factors,
/// rounded once (nearest, ties to even; subnormal results on the 2⁻¹⁰⁷⁴ grid; ±∞ past DBL_MAX), its sign the xor of the
/// factors' signs (a zero factor, or a value that rounds to zero, gives that signed zero): Π Mᵢ over Π Mⱼ·2^(ΣEᵢ − ΣEⱼ)
/// by [`cr_ratio`]. Integer arithmetic only, independent of the port's double-double fast path and its multi-word
/// division (`validate.py` checks every value again with `Fraction`).
pub fn cr_monomial(num: &[f64], den: &[f64]) -> f64 {
    assert!(
        num.iter().chain(den).all(|v| v.is_finite()) && den.iter().all(|v| *v != 0.0),
        "cr_monomial({num:?}, {den:?})"
    );
    let neg = num.iter().chain(den).fold(false, |a, v| a != v.is_sign_negative());
    let sign = if neg { -1.0 } else { 1.0 };
    if num.contains(&0.0) {
        return sign * 0.0;
    }
    let (mut p, mut q, mut e) = (Big::new(1), Big::new(1), 0_i64);
    for &x in num {
        let (m, ex) = decompose(x);
        p = p.mul(&Big::new(u128::from(m)));
        e += ex;
    }
    for &y in den {
        let (m, ey) = decompose(y);
        q = q.mul(&Big::new(u128::from(m)));
        e -= ey;
    }
    sign * cr_ratio(&p, &q, e)
}

/// The overflow audit's B2 (ADR 0002 #30): RN(√(x/y)) for positive finite doubles, rounded once (nearest; subnormal
/// results on the 2⁻¹⁰⁷⁴ grid; +∞ past DBL_MAX): with x/y = mₓ/m_y·2^(eₓ − e_y), Q = ⌊mₓ·2^k/m_y⌋ of 121 or 122 bits
/// (k ≥ 0 with eₓ − e_y − k even), R = ⌊√Q⌋ (of at least 60 bits; ⌊√⌊z⌋⌋ = ⌊√z⌋), sticky = [R² ≠ Q or the division
/// left a remainder], then `round_scaled(2R + sticky)` (which drops at least eight bits, so the sticky bit decides every
/// rounding: v·2^-(h − 1) lies in [2R, 2R + 2), strictly inside when sticky). Integer arithmetic only, independent of
/// the port's double-double fast path and its multi-word root (`validate.py` checks every value again with Python's
/// integers and MPFR).
pub fn cr_sqrt_ratio(x: f64, y: f64) -> f64 {
    assert!(x > 0.0 && x.is_finite() && y > 0.0 && y.is_finite(), "cr_sqrt_ratio({x:?}, {y:?})");
    let (mx, ex) = decompose(x);
    let (my, ey) = decompose(y);
    let (p, q) = (Big::new(u128::from(mx)), Big::new(u128::from(my)));
    let mut k = (121 + i64::try_from(q.bits()).unwrap_or(0) - i64::try_from(p.bits()).unwrap_or(0)).max(0);
    if (ex - ey - k).rem_euclid(2) == 1 {
        k += 1;
    }
    let (quo, rem) = p.shl(u64::try_from(k).unwrap_or(0)).divrem(&q);
    assert!(quo.bits() <= 124, "cr_sqrt_ratio: Q has {} bits", quo.bits());
    let qq = u128::from(quo.0.first().copied().unwrap_or(0)) | (u128::from(quo.0.get(1).copied().unwrap_or(0)) << 64);
    let r = isqrt(qq);
    let sticky = u128::from(r * r != qq || !rem.is_zero());
    round_scaled((r << 1) | sticky, (ex - ey - k) / 2 - 1, false)
}

/// (2R + sticky, h − 1) for v = √(p/q·2^e), p ≥ 1, q ≥ 1 integers of any size: Q = ⌊p·2^k/q⌋ of 121 or 122 bits (k of
/// either sign, e − k even), R = ⌊√Q⌋ (of at least 60 bits; ⌊√⌊z⌋⌋ = ⌊√z⌋), sticky = [R² ≠ Q or the division left a
/// remainder], h = (e − k)/2: v·2^(1 − h) lies in [2R, 2R + 2), strictly inside when sticky, so a rounding of
/// (2R + sticky)·2^(h − 1) that drops at least two bits is v's (the overflow audit's B3: the Morse width's √ of a
/// monomial, Redlich–Kwong's Tc^2.5).
fn sqrt_parts_big(p: &Big, q: &Big, e: i64) -> (u128, i64) {
    let mut k = 121 + i64::try_from(q.bits()).unwrap_or(0) - i64::try_from(p.bits()).unwrap_or(0);
    if (e - k).rem_euclid(2) == 1 {
        k += 1;
    }
    let (num, den) = if k >= 0 {
        (p.shl(u64::try_from(k).unwrap_or(0)), q.clone())
    } else {
        (p.clone(), q.shl(u64::try_from(-k).unwrap_or(0)))
    };
    let (quo, rem) = num.divrem(&den);
    assert!((121..=123).contains(&quo.bits()), "sqrt_parts_big: Q has {} bits", quo.bits());
    let qq = u128::from(quo.0.first().copied().unwrap_or(0)) | (u128::from(quo.0.get(1).copied().unwrap_or(0)) << 64);
    let r = isqrt(qq);
    let sticky = u128::from(r * r != qq || !rem.is_zero());
    ((r << 1) | sticky, (e - k) / 2 - 1)
}

/// The overflow audit's B3 (ADR 0002 #32): RN(√(Π num / Π den)) for positive finite doubles, rounded once (nearest, ties
/// to even: unlike #30's √(x/y) the value can be a midpoint, ω·√(9/(2·½)) = 3ω; subnormal results on the 2⁻¹⁰⁷⁴ grid;
/// +∞ past DBL_MAX): the Morse width ω·√(μ/(2·D_e)) as √(ω·ω·μ/(2·D_e)). The mantissas' products over each other by
/// [`sqrt_parts_big`], then `round_scaled` (which drops at least eight bits, so the sticky bit decides every rounding,
/// an exact midpoint to even). Integer arithmetic only, independent of the port's double-double fast path and its
/// multi-word root (`validate.py` checks every value again with Python's integers and MPFR).
pub fn cr_sqrt_monomial(num: &[f64], den: &[f64]) -> f64 {
    assert!(num.iter().chain(den).all(|v| v.is_finite() && *v > 0.0), "cr_sqrt_monomial({num:?}, {den:?})");
    let (mut p, mut q, mut e) = (Big::new(1), Big::new(1), 0_i64);
    for &x in num {
        let (m, ex) = decompose(x);
        p = p.mul(&Big::new(u128::from(m)));
        e += ex;
    }
    for &y in den {
        let (m, ey) = decompose(y);
        q = q.mul(&Big::new(u128::from(m)));
        e -= ey;
    }
    let (mant, h) = sqrt_parts_big(&p, &q, e);
    round_scaled(mant, h, false)
}

/// a·2^ea against b·2^eb, for non-negative integers a and b, exactly.
fn cmp_scaled(a: &Big, ea: i64, b: &Big, eb: i64) -> std::cmp::Ordering {
    let e = ea.min(eb);
    a.shl(u64::try_from(ea - e).unwrap_or(0)).cmp(&b.shl(u64::try_from(eb - e).unwrap_or(0)))
}

/// ADR 0002 #25: m₁/c₁ < m₂/c₂ exactly, for finite m₁, m₂ ≥ 0 and finite c₁, c₂ > 0 (m₁·c₂ < m₂·c₁; a zero m is the
/// least, two zeros tie).
pub fn ratio_lt(m1: f64, c1: f64, m2: f64, c2: f64) -> bool {
    if m1 == 0.0 {
        return m2 != 0.0;
    }
    if m2 == 0.0 {
        return false;
    }
    let ((a1, e1), (b1, f1), (a2, e2), (b2, f2)) = (decompose(m1), decompose(c1), decompose(m2), decompose(c2));
    let l = Big::new(u128::from(a1) * u128::from(b2));
    let r = Big::new(u128::from(a2) * u128::from(b1));
    cmp_scaled(&l, e1 + f2, &r, e2 + f1) == std::cmp::Ordering::Less
}

/// ADR 0002 #26: the empirical formula of (atomic mass, mass percent) pairs, each positive and finite, evaluated exactly:
/// the moles p/m and their minimum (element j, by [`ratio_lt`]), the ratios rᵢ = pᵢ·m_j/(mᵢ·p_j) as integer fractions
/// A/B (the exponents folded in, at any size: the port's 2³² − ½ cut is checked by `validate.py`, not assumed here), the
/// smallest k in 1 ..= 6 with |rᵢ·k − round(rᵢ·k)| ≤ 1/10 for every i (k·A = q·B + R: 10·R ≤ B or 10·R ≥ 9·B), and the
/// counts round(rᵢ·k), halves away from zero (q, or q + 1 when 2·R ≥ B). None for Rust's ComputationError: no such k,
/// or a count past u32::MAX.
pub fn empirical_exact(mp: &[(f64, f64)]) -> Option<Vec<u32>> {
    use std::cmp::Ordering::{Greater, Less};
    let mut j = 0;
    for k in 1..mp.len() {
        if ratio_lt(mp[k].1, mp[k].0, mp[j].1, mp[j].0) {
            j = k;
        }
    }
    let (pj, mj) = (decompose(mp[j].1), decompose(mp[j].0));
    let fr: Vec<(Big, Big)> = mp
        .iter()
        .map(|&(m, p)| {
            let (pi, mi) = (decompose(p), decompose(m));
            let s = pi.1 + mj.1 - mi.1 - pj.1;
            let a = Big::new(u128::from(pi.0) * u128::from(mj.0));
            let b = Big::new(u128::from(mi.0) * u128::from(pj.0));
            if s >= 0 {
                (a.shl(u64::try_from(s).unwrap_or(0)), b)
            } else {
                (a, b.shl(u64::try_from(-s).unwrap_or(0)))
            }
        })
        .collect();
    let (ten, nine, two) = (Big::new(10), Big::new(9), Big::new(2));
    let k = (1..=6_u128).find(|&k| {
        fr.iter().all(|(a, b)| {
            let (_, r) = a.mul(&Big::new(k)).divrem(b);
            let r10 = r.mul(&ten);
            r10.cmp(b) != Greater || r10.cmp(&b.mul(&nine)) != Less
        })
    })?;
    fr.iter()
        .map(|(a, b)| {
            let (q, r) = a.mul(&Big::new(k)).divrem(b);
            let c = if r.mul(&two).cmp(b) == Less { q } else { q.add(&Big::new(1)) };
            if c.bits() > 32 { None } else { u32::try_from(c.0.first().copied().unwrap_or(0)).ok() }
        })
        .collect()
}

/// The overflow audit's B6 (ADR 0002 #37): whether w > lo·(1 − t) + hi·t exactly, for finite doubles w, lo, hi and t
/// (`find_peaks`' threshold with 1 − t not rounded). Every double is an integer multiple of 2⁻¹⁰⁷⁴, so with X = x·2¹⁰⁷⁴
/// the test w + t·lo > lo + t·hi is W·2¹⁰⁷⁴ + T·L > L·2¹⁰⁷⁴ + T·H in integers, each side a sum of signed products of
/// non-negative `Big`s. Integer arithmetic only, independent of the port's scaled double expansions (`validate.py`
/// checks every value again with `Fraction`).
pub fn peak_above(w: f64, lo: f64, hi: f64, t: f64) -> bool {
    let fin = w.is_finite() && lo.is_finite() && hi.is_finite() && t.is_finite();
    assert!(fin, "peak_above({w:?}, {lo:?}, {hi:?}, {t:?})");
    // |x|·2¹⁰⁷⁴ and x's sign (zero as non-negative)
    let int = |x: f64| -> (bool, Big) {
        if x == 0.0 {
            return (false, Big::new(0));
        }
        let (m, e) = decompose(x);
        (x < 0.0, Big::new(u128::from(m)).shl(u64::try_from(e + 1074).unwrap_or(0)))
    };
    let (wn, wm) = int(w);
    let (ln, lm) = int(lo);
    let (hn, hm) = int(hi);
    let (tn, tm) = int(t);
    let p1074 = Big::new(1).shl(1074);
    // the signed terms of the left side W·2¹⁰⁷⁴ + T·L and the right side L·2¹⁰⁷⁴ + T·H
    let left = [(wn, wm.mul(&p1074)), (tn != ln, tm.mul(&lm))];
    let right = [(ln, lm.mul(&p1074)), (tn != hn, tm.mul(&hm))];
    // Σ positive − Σ negative of (left − right), compared as two non-negative sums
    let (mut pos, mut neg) = (Big::new(0), Big::new(0));
    for (n, v) in left {
        if n { neg = neg.add(&v) } else { pos = pos.add(&v) }
    }
    for (n, v) in right {
        if n { pos = pos.add(&v) } else { neg = neg.add(&v) }
    }
    pos.cmp(&neg) == std::cmp::Ordering::Greater
}

// ── The overflow audit's B7 (ADR 0002 #39): signed integers of any size, for balance_equation's exact replica ──────────

impl Big {
    /// (⌊self / d⌋, self mod d) for d > 0 by Knuth's Algorithm D in base 2⁶⁴ (u128 intermediates): [`Big::divrem`]'s
    /// binary long division is too slow for the thousand-bit operands of ADR 0002 #39's replica.
    fn divrem_d(&self, d: &Self) -> (Self, Self) {
        assert!(!d.is_zero(), "Big::divrem_d by zero");
        if self.cmp(d) == std::cmp::Ordering::Less {
            return (Self(Vec::new()), self.clone());
        }
        let n = d.0.len();
        if n == 1 {
            let dv = u128::from(d.0[0]);
            let mut q = vec![0_u64; self.0.len()];
            let mut r = 0_u128;
            for k in (0..self.0.len()).rev() {
                let cur = (r << 64) | u128::from(self.0[k]);
                #[allow(clippy::cast_possible_truncation)] // cur / dv < 2⁶⁴: r < dv
                let qd = (cur / dv) as u64;
                q[k] = qd;
                r = cur % dv;
            }
            return (Self(q).norm(), Self::new(r));
        }
        let s = d.0[n - 1].leading_zeros();
        let vn = d.shl(u64::from(s)).0;
        let mut un = self.shl(u64::from(s)).0;
        un.resize(self.0.len() + 1, 0);
        let m = self.0.len() - n;
        let mut q = vec![0_u64; m + 1];
        let b = 1_u128 << 64;
        for j in (0..=m).rev() {
            let num = (u128::from(un[j + n]) << 64) | u128::from(un[j + n - 1]);
            let mut qh = num / u128::from(vn[n - 1]);
            let mut rh = num % u128::from(vn[n - 1]);
            while qh >= b || qh * u128::from(vn[n - 2]) > ((rh << 64) | u128::from(un[j + n - 2])) {
                qh -= 1;
                rh += u128::from(vn[n - 1]);
                if rh >= b {
                    break;
                }
            }
            // un[j ..= j + n] −= q̂·vn, with the borrow as a signed accumulator
            let mut borrow = 0_i128;
            let mut carry = 0_u128;
            for i in 0..n {
                let p = qh * u128::from(vn[i]) + carry;
                carry = p >> 64;
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let t = i128::from(un[i + j]) - i128::from(p as u64) + borrow;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let lo = t as u64;
                un[i + j] = lo;
                borrow = t >> 64;
            }
            #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
            let t = i128::from(un[j + n]) - (carry as i128) + borrow;
            if t < 0 {
                // q̂ one too large: add vn back
                qh -= 1;
                let mut c = 0_u128;
                for i in 0..n {
                    let sm = u128::from(un[i + j]) + u128::from(vn[i]) + c;
                    #[allow(clippy::cast_possible_truncation)]
                    let lo = sm as u64;
                    un[i + j] = lo;
                    c = sm >> 64;
                }
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let top = (t + i128::from(c != 0)) as u64;
                un[j + n] = top;
            } else {
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let top = t as u64;
                un[j + n] = top;
            }
            #[allow(clippy::cast_possible_truncation)]
            let qd = qh as u64;
            q[j] = qd;
        }
        un.truncate(n);
        (Self(q).norm(), Self(un).norm().shr(u64::from(s)))
    }

    /// gcd(self, o) by Euclid's algorithm on [`Big::divrem_d`]; gcd(x, 0) = x.
    fn gcd(&self, o: &Self) -> Self {
        let (mut a, mut b) = (self.clone(), o.clone());
        while !b.is_zero() {
            let (_, r) = a.divrem_d(&b);
            a = b;
            b = r;
        }
        a
    }
}

/// A signed integer of any size: a sign and a magnitude (zero is non-negative). ADR 0002 #39's replica of the port's
/// multi-word balance_equation (`audit::port_balance`): only what Rust's elimination needs. Integer arithmetic only,
/// independent of the port's 31-bit limbs and Lehmer gcd (`validate.py` checks every result again with Python ints).
#[derive(Clone)]
pub struct Int {
    neg: bool,
    mag: Big,
}

impl Int {
    pub fn from_i64(v: i64) -> Self {
        Self { neg: v < 0, mag: Big::new(u128::from(v.unsigned_abs())) }
    }

    fn of(neg: bool, mag: Big) -> Self {
        let neg = neg && !mag.is_zero();
        Self { neg, mag }
    }

    pub fn is_zero(&self) -> bool {
        self.mag.is_zero()
    }

    pub fn is_negative(&self) -> bool {
        self.neg
    }

    pub fn is_positive(&self) -> bool {
        !self.neg && !self.mag.is_zero()
    }

    #[must_use]
    pub fn neg(&self) -> Self {
        Self::of(!self.neg, self.mag.clone())
    }

    #[must_use]
    pub fn mul(&self, o: &Self) -> Self {
        Self::of(self.neg != o.neg, self.mag.mul(&o.mag))
    }

    #[must_use]
    pub fn sub(&self, o: &Self) -> Self {
        use std::cmp::Ordering::Less;
        if self.neg != o.neg {
            return Self::of(self.neg, self.mag.add(&o.mag)); // a − b = sign(a)·(|a| + |b|)
        }
        if self.mag.cmp(&o.mag) == Less {
            Self::of(!self.neg, o.mag.sub(&self.mag))
        } else {
            Self::of(self.neg, self.mag.sub(&o.mag))
        }
    }

    /// gcd(|self|, |o|), non-negative.
    #[must_use]
    pub fn gcd(&self, o: &Self) -> Self {
        Self::of(false, self.mag.gcd(&o.mag))
    }

    /// self / d for a d ≠ 0 that divides self (asserted): Rust's truncating `/`, exact here.
    #[must_use]
    pub fn div_exact(&self, d: &Self) -> Self {
        let (q, r) = self.mag.divrem_d(&d.mag);
        assert!(r.is_zero(), "Int::div_exact: not a divisor");
        Self::of(self.neg != d.neg, q)
    }

    /// Whether this is 1 exactly.
    pub fn is_one(&self) -> bool {
        !self.neg && self.mag.bits() == 1
    }

    /// Whether this lies in i64's range.
    pub fn fits_i64(&self) -> bool {
        let b = self.mag.bits();
        b <= 63 || (self.neg && b == 64 && self.mag.0.len() == 1 && self.mag.0[0] == 1 << 63)
    }

    /// The value as a u32, if it is one (0 ..= u32::MAX).
    pub fn to_u32(&self) -> Option<u32> {
        if self.neg || self.mag.bits() > 32 {
            return None;
        }
        u32::try_from(self.mag.0.first().copied().unwrap_or(0)).ok()
    }
}
