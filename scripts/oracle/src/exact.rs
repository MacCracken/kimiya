//! Exact integer arithmetic for the exact values ADR 0002 #13 and #16 name (M4c, M4d): binary64 rounding of
//! m·2^e for a mantissa of any size, and the two irrational factors the port rounds once, 2⁻ᶠ (N₀/2ⁿ for a fractional
//! n) and 2^(1/6) (σ·2^(1/6)). Integer arithmetic only, independent of the port's double-word code (`validate.py`
//! checks every value again, with `Fraction` and MPFR).

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
