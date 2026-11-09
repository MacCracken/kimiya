//! Crate-internal input-validation predicates.
//!
//! Comparisons such as `x <= 0.0` are `false` for NaN, so a guard written
//! `if x <= 0.0 { return Err(..) }` lets NaN (and ±∞) through and the
//! function silently returns `Ok(NaN)`. Guards use these predicates in the
//! negated form instead: `if !is_positive(x) { return Err(..) }`.

/// `x` is finite and strictly greater than zero.
#[must_use]
#[inline]
pub(crate) const fn is_positive(x: f64) -> bool {
    x > 0.0 && x.is_finite()
}

/// `x` is finite and greater than or equal to zero.
#[must_use]
#[inline]
pub(crate) const fn is_non_negative(x: f64) -> bool {
    x >= 0.0 && x.is_finite()
}

/// `x` is neither NaN nor ±∞.
#[must_use]
#[inline]
pub(crate) const fn is_finite(x: f64) -> bool {
    x.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_rejects_nan_inf_zero_negative() {
        assert!(is_positive(1e-300));
        assert!(!is_positive(0.0));
        assert!(!is_positive(-1.0));
        assert!(!is_positive(f64::NAN));
        assert!(!is_positive(f64::INFINITY));
    }

    #[test]
    fn non_negative_accepts_zero() {
        assert!(is_non_negative(0.0));
        assert!(is_non_negative(-0.0));
        assert!(!is_non_negative(-1e-300));
        assert!(!is_non_negative(f64::NAN));
        assert!(!is_non_negative(f64::INFINITY));
    }

    #[test]
    fn finite_rejects_nan_inf() {
        assert!(is_finite(-1.0));
        assert!(!is_finite(f64::NAN));
        assert!(!is_finite(f64::NEG_INFINITY));
    }
}
