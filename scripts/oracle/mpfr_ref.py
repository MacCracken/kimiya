"""Correctly rounded binary64 results from the system MPFR library (ADR 0006).

Each function returns the IEEE-754 binary64 round-to-nearest-even value of the exact mathematical
result: MPFR at 53 bits of precision with binary64's exponent range (emin = -1073, emax = 1024) and
`mpfr_subnormalize`, so subnormal results, overflow to infinity and exact ties (pow(3*2^-215, 5),
exp10(23)) round exactly as IEEE round-to-nearest does. A decimal reference at fixed precision can't
do that: it misjudges exact ties (it rounds pow(0.5, 1075) = 2^-1075 up instead of to 0).

UX (M4e-2, ADR 0002 #20) is binary64 with an unbounded exponent range: 53-bit MPFR numbers operated on at MPFR's widest
exponent range, the final operation rounded once into binary64 (ux_*_final), and potential.rs's geometry on them
(potential_*), the reference for the port's unbounded-exponent helpers and geometry; since the overflow audit's B2 also
the final subtraction (ux_sub_final) and the √ of a quotient (ux_sqrt_div_final).

No Python package is needed: MPFR (libmpfr.so.6, MPFR 4.x) is called through ctypes. Run
`python3 -I scripts/oracle/mpfr_ref.py` for the self-checks.
"""
import ctypes
import ctypes.util
import math
import struct

_lib = ctypes.util.find_library("mpfr") or "libmpfr.so.6"
_m = ctypes.CDLL(_lib)

# mpfr_t on LP64: _mpfr_prec (long), _mpfr_sign (int, padded), _mpfr_exp (long), _mpfr_d (limb *).
_MPFR_T = ctypes.c_byte * 32
_RNDN = 0
_P = ctypes.c_void_p

_m.mpfr_get_version.restype = ctypes.c_char_p
_m.mpfr_set_emin.argtypes = [ctypes.c_long]
_m.mpfr_set_emax.argtypes = [ctypes.c_long]
_m.mpfr_init2.argtypes = [_P, ctypes.c_long]
_m.mpfr_set_d.argtypes = [_P, ctypes.c_double, ctypes.c_int]
_m.mpfr_get_d.argtypes = [_P, ctypes.c_int]
_m.mpfr_get_d.restype = ctypes.c_double
_m.mpfr_subnormalize.argtypes = [_P, ctypes.c_int, ctypes.c_int]

_UNARY = ("exp", "expm1", "exp2", "exp10", "log", "log2", "log10", "log1p",
          "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "cbrt", "sqrt")
_BINARY = ("pow", "atan2", "hypot")
for _name in _UNARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, ctypes.c_int]
for _name in _BINARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, _P, ctypes.c_int]

_m.mpfr_set_emin(-1073)
_m.mpfr_set_emax(1024)
_a = _MPFR_T()
_b = _MPFR_T()
_r = _MPFR_T()
for _v in (_a, _b, _r):
    _m.mpfr_init2(ctypes.byref(_v), 53)

VERSION = _m.mpfr_get_version().decode()


def _finish(t):
    t = _m.mpfr_subnormalize(ctypes.byref(_r), t, _RNDN)
    return _m.mpfr_get_d(ctypes.byref(_r), _RNDN)


def cr(name, x):
    """The correctly rounded binary64 value of the unary MPFR function `name` at x (a float)."""
    if name not in _UNARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), _RNDN)
    return _finish(t)


def cr2(name, x, y):
    """The correctly rounded binary64 value of the binary MPFR function `name` (pow, atan2, hypot)."""
    if name not in _BINARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    _m.mpfr_set_d(ctypes.byref(_b), y, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), ctypes.byref(_b), _RNDN)
    return _finish(t)


# x·2^y rounded once (ADR 0002 #13, M4d): the port's _km_cr_exp2_scaled. 2^f (f = y − trunc(y), |f| < 1) to _HP bits,
# bounded from below and from above (RNDD, RNDU); each bound times x's significand exactly (_HP + 53 bits), then scaled
# by 2^(trunc(y) + x's exponent) into a 53-bit result with binary64's exponent range and subnormalized, as cr() rounds. The two must agree, which they
# do unless x·2^y is within 2^-_HP (relative) of a rounding boundary: for f ≠ 0 it is irrational (2^q is, for every
# rational q that is not an integer), so a disagreement means only that _HP is too small, and it raises.
_HP = 400
_RNDU, _RNDD = 2, 3
_m.mpfr_mul_d.argtypes = [_P, _P, ctypes.c_double, ctypes.c_int]
_m.mpfr_mul_2si.argtypes = [_P, _P, ctypes.c_long, ctypes.c_int]
_m.mpfr_exp2.argtypes = [_P, _P, ctypes.c_int]
_f = _MPFR_T()
_e2 = _MPFR_T()
_xe = _MPFR_T()
_m.mpfr_init2(ctypes.byref(_f), 53)
_m.mpfr_init2(ctypes.byref(_e2), _HP)
_m.mpfr_init2(ctypes.byref(_xe), _HP + 60)


def _scaled_bound(x, f, k, rnd):
    _m.mpfr_set_d(ctypes.byref(_f), f, _RNDN)
    _m.mpfr_exp2(ctypes.byref(_e2), ctypes.byref(_f), rnd)                    # 2^f in [1/2, 2], bounded
    _m.mpfr_mul_d(ctypes.byref(_xe), ctypes.byref(_e2), x, _RNDN)             # exact: _HP + 53 bits fit
    t = _m.mpfr_mul_2si(ctypes.byref(_r), ctypes.byref(_xe), k, _RNDN)        # one rounding to 53 bits, binary64 range
    return _finish(t)


def scaled_exp2(x, y):
    """RN(x·2^y), rounded once (subnormal results at the 2^-1074 grid, ±∞ past DBL_MAX), for finite x and y; x = ±0
    gives x. The correctly rounded binary64 value of the exact product."""
    import math
    if x == 0.0 or not math.isfinite(x) or not math.isfinite(y):
        raise ValueError((x, y))
    k = math.trunc(y)
    f = y - k                                       # exact
    m, e = math.frexp(x)                            # x = m·2^e, |m| in [1/2, 1): the intermediates stay in range
    if f == 0.0:
        return _scaled_bound(m, 0.0, k + e, _RNDN)  # 2^0 = 1: exact, one rounding
    lo = _scaled_bound(m, f, k + e, _RNDD)
    hi = _scaled_bound(m, f, k + e, _RNDU)
    if bits(lo) != bits(hi):
        raise ArithmeticError(f"scaled_exp2({x!r}, {y!r}) undecided at {_HP} bits")
    return lo


# 1/(100·x) rounded once (ADR 0002 #19, M4e): the port's _km_cr_hundredth_div, spectroscopy's wavenumber conversions.
# 100·x is exact at 64 bits (x has 53 and 100 seven), then one division of 1 by it to 53 bits. 100·x can pass binary64's
# exponent range (x near DBL_MAX), so both run in MPFR's widest range; then mpfr_check_range and mpfr_subnormalize take
# the result (with its ternary value) into binary64's range, MPFR's documented emulation of IEEE rounding, as cr()
# rounds. ±0 gives ±∞, ±∞ gives ±0, NaN NaN (IEEE's 1/(100·x)).
_m.mpfr_ui_div.argtypes = [_P, ctypes.c_ulong, _P, ctypes.c_int]
_m.mpfr_check_range.argtypes = [_P, ctypes.c_int, ctypes.c_int]
_m.mpfr_get_emin_min.restype = ctypes.c_long
_m.mpfr_get_emax_max.restype = ctypes.c_long
_h100 = _MPFR_T()
_m.mpfr_init2(ctypes.byref(_h100), 64)


def hundredth_div(x):
    """RN(1/(100·x)), the correctly rounded (1/100)/x."""
    _m.mpfr_set_emin(_m.mpfr_get_emin_min())
    _m.mpfr_set_emax(_m.mpfr_get_emax_max())
    try:
        _m.mpfr_set_d(ctypes.byref(_h100), x, _RNDN)
        _m.mpfr_mul_d(ctypes.byref(_h100), ctypes.byref(_h100), 100.0, _RNDN)      # exact: 60 bits
        t = _m.mpfr_ui_div(ctypes.byref(_r), 1, ctypes.byref(_h100), _RNDN)
    finally:
        _m.mpfr_set_emin(-1073)
        _m.mpfr_set_emax(1024)
    t = _m.mpfr_check_range(ctypes.byref(_r), t, _RNDN)
    return _finish(t)


# Binary64 with an unbounded exponent range (ADR 0002 #20, M4e-2): the port's potential geometry evaluates Rust's operation
# sequence with every + − × ÷ √ rounded to 53 bits, to nearest-even, and no overflow or underflow, then rounds the final
# operation (√ or atan2) once into binary64. UX is that value: an MPFR number at 53 bits, each operation done while MPFR's
# exponent range is its widest (emin and emax at mpfr_get_emin_min / _emax_max, about ±2^62: nothing a polynomial of
# doubles reaches can over- or underflow there), so it is exactly the 53-bit rounding of the exact result; MPFR keeps
# IEEE's signed zeros (x − x = +0 under round-to-nearest, a product's sign is the xor, √−0 = −0). The final operations
# (ux_sqrt_final, ux_atan2_final) are computed in the same range with their ternary value, then mpfr_check_range and
# mpfr_subnormalize take the result into binary64's range, the IEEE emulation hundredth_div uses: one rounding, the
# subnormal grid and ±∞ included. Independent of the port's (significand, exponent) doubles and of the generator's
# integer arithmetic.
_UXOPS = ("add", "sub", "mul", "div")
for _name in _UXOPS:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_clear.argtypes = [_P]
_m.mpfr_neg.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_cmp.argtypes = [_P, _P]
_m.mpfr_sgn.argtypes = [_P]
_m.mpfr_zero_p.argtypes = [_P]
_m.mpfr_signbit.argtypes = [_P]
_m.mpfr_number_p.argtypes = [_P]
_m.mpfr_get_exp.argtypes = [_P]
_m.mpfr_get_exp.restype = ctypes.c_long
_UX_EMIN = _m.mpfr_get_emin_min()
_UX_EMAX = _m.mpfr_get_emax_max()


_WIDE = [0]     # the depth of `with unbounded():` (MPFR's exponent range is the widest inside it)


def _wide():
    _m.mpfr_set_emin(_UX_EMIN)
    _m.mpfr_set_emax(_UX_EMAX)


def _narrow():
    _m.mpfr_set_emin(-1073)
    _m.mpfr_set_emax(1024)


class unbounded:
    """`with unbounded():` keeps MPFR's widest exponent range for a whole evaluation (the UX operations otherwise switch
    to it and back around each call; cr() and the other functions need binary64's range outside)."""
    def __enter__(self):
        if _WIDE[0] == 0:
            _wide()
        _WIDE[0] += 1
        return self

    def __exit__(self, *exc):
        _WIDE[0] -= 1
        if _WIDE[0] == 0:
            _narrow()
        return False


_UX_FN = {name: getattr(_m, "mpfr_" + name) for name in _UXOPS}


class UX:
    """A binary64 value with an unbounded exponent range (an MPFR number at 53 bits). UX(x) is the double x exactly."""
    __slots__ = ("v",)

    def __init__(self, x=None):
        self.v = _MPFR_T()
        _m.mpfr_init2(self.v, 53)
        if x is not None:
            if not math.isfinite(x):
                raise ValueError(x)
            _m.mpfr_set_d(self.v, x, _RNDN)                    # exact: 53 bits, a double's exponent

    def __del__(self):
        _m.mpfr_clear(self.v)

    def _op(self, name, o):
        r = UX()
        with unbounded():
            _UX_FN[name](r.v, self.v, o.v, _RNDN)
        return r

    def __add__(self, o): return self._op("add", o)
    def __sub__(self, o): return self._op("sub", o)
    def __mul__(self, o): return self._op("mul", o)
    def __truediv__(self, o): return self._op("div", o)

    def sqrt(self):
        r = UX()
        with unbounded():
            _m.mpfr_sqrt(r.v, self.v, _RNDN)
        return r

    def cmp(self, o):
        """−1, 0 or 1 by exact value (zeros equal)."""
        with unbounded():
            c = _m.mpfr_cmp(self.v, o.v)
        return (c > 0) - (c < 0)

    def __le__(self, o): return self.cmp(o) <= 0
    def __lt__(self, o): return self.cmp(o) < 0
    def __gt__(self, o): return self.cmp(o) > 0

    def max(self, o):
        """Rust's x.max(y) on non-NaN values: x where x > y, y otherwise."""
        return self if self.cmp(o) > 0 else o

    def is_positive(self):
        return _m.mpfr_sgn(self.v) > 0

    def is_zero(self):
        return _m.mpfr_zero_p(self.v) != 0

    def neg_sign(self):
        return _m.mpfr_signbit(self.v) != 0

    def is_f64(self, x):
        """True if this is exactly the double x (a finite x, zeros by sign)."""
        if not math.isfinite(x):
            return False
        if self.is_zero() or x == 0.0:
            return self.is_zero() and x == 0.0 and self.neg_sign() == (math.copysign(1.0, x) < 0)
        e = _m.mpfr_get_exp(self.v)                            # value in [2^(e − 1), 2^e)
        if e < -1073 or e > 1024:
            return False
        return self.cmp(UX(x)) == 0


def _final(fn, *args):
    r = UX()
    with unbounded():
        t = fn(r.v, *[a.v for a in args], _RNDN)
        _narrow()
        try:
            t = _m.mpfr_check_range(r.v, t, _RNDN)
            t = _m.mpfr_subnormalize(r.v, t, _RNDN)
            return _m.mpfr_get_d(r.v, _RNDN)
        finally:
            _wide()


def ux_sqrt_final(a):
    """√a rounded once into binary64 (subnormal results at the 2^-1074 grid, ±∞ past DBL_MAX), a ≥ 0."""
    return _final(_m.mpfr_sqrt, a)


def ux_atan2_final(y, x):
    """atan2(y, x) rounded once into binary64."""
    return _final(_m.mpfr_atan2, y, x)


def ux_mul_final(a, b):
    """a·b rounded once into binary64."""
    return _final(_UX_FN["mul"], a, b)


def ux_div_final(a, b):
    """a/b rounded once into binary64."""
    return _final(_UX_FN["div"], a, b)


def ux_add_final(a, b):
    """a + b rounded once into binary64."""
    return _final(_UX_FN["add"], a, b)


def ux_sub_final(a, b):
    """a − b rounded once into binary64 (MPFR's subtraction: IEEE's signed zeros, x − x = +0): ADR 0002 #29's Nernst
    equation, E° − ((R·T)/(n·F))·ln Q (the overflow audit's B2)."""
    return _final(_UX_FN["sub"], a, b)


class _Wide:
    """An MPFR number of `prec` bits (the quotient of ux_sqrt_div_final)."""
    __slots__ = ("v",)

    def __init__(self, prec):
        self.v = _MPFR_T()
        _m.mpfr_init2(self.v, prec)

    def __del__(self):
        _m.mpfr_clear(self.v)


def ux_sqrt_div_final(a, b, prec=256):
    """√(a/b) rounded once into binary64 for positive a, b (ADR 0002 #30, the overflow audit's B2: vibrational_frequency's
    √(k/μ)): the quotient at `prec` bits (relative error at most 2^-prec, MPFR's widest exponent range), then its square
    root rounded once (_final). For doubles a, b this is RN(√(a/b)) itself when prec ≥ 166: with a/b = (p/q)·2^E (p, q
    odd, below 2^53) and a rounding boundary m = M·2^f (M odd, at most 54 bits), p·2^E − q·M²·2^(2f) is a nonzero
    multiple of 2^min(E, 2f) (validate.py's §20 proves it is never 0), so a/b is at least 2^-161 (relative) from every m²
    and √(a/b) at least 2^-162 from every m, far beyond the quotient's error."""
    q = _Wide(prec)
    with unbounded():
        _UX_FN["div"](q.v, a.v, b.v, _RNDN)
    return _final(_m.mpfr_sqrt, q)


def ux_scaled(s, e):
    """s·2^e as a UX, exactly, for a double s and any integer e."""
    r = UX(s)
    with unbounded():
        _m.mpfr_mul_2si(r.v, r.v, e, _RNDN)     # exact at the widest range (53 bits kept)
    return r


# potential.rs's geometry (rust-old/src/potential.rs, Rust 1.2.0's operation order: sub3, cross, dot3 as
# (a0·b0 + a1·b1) + a2·b2) on UX values, for ADR 0002 #20's reference: the coordinates (finite doubles, flattened) in, the
# result rounded once into binary64 out, None for Rust's Err. `tol` is Rust's COINCIDENT_REL_TOL / COLLINEAR_SIN_SQ_TOL
# (validate.py passes the values it parses from the source).
def _ux3(p): return [UX(x) for x in p]
def _ux_sub(a, b): return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
def _ux_dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def _ux_cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def potential_distance_3d(v):
    with unbounded():
        d = _ux_sub(_ux3(v[3:6]), _ux3(v[0:3]))
        return ux_sqrt_final(_ux_dot(d, d))


def potential_bond_angle(v, tol):
    with unbounded():
        a, b, c = _ux3(v[0:3]), _ux3(v[3:6]), _ux3(v[6:9])
        ba, bc = _ux_sub(a, b), _ux_sub(c, b)
        ba_sq, bc_sq = _ux_dot(ba, ba), _ux_dot(bc, bc)
        t = UX(tol)
        tol_sq = t * t
        b_sq = _ux_dot(b, b)
        if (not ba_sq.is_positive() or not bc_sq.is_positive() or ba_sq <= tol_sq * _ux_dot(a, a).max(b_sq)
                or bc_sq <= tol_sq * _ux_dot(c, c).max(b_sq)):
            return None
        cr = _ux_cross(ba, bc)
        return ux_atan2_final(_ux_dot(cr, cr).sqrt(), _ux_dot(ba, bc))


def potential_dihedral_angle(v, tol, pi=math.pi):
    with unbounded():
        a, b, c, d = _ux3(v[0:3]), _ux3(v[3:6]), _ux3(v[6:9]), _ux3(v[9:12])
        b1, b2, b3 = _ux_sub(b, a), _ux_sub(c, b), _ux_sub(d, c)
        n1, n2 = _ux_cross(b1, b2), _ux_cross(b2, b3)
        n1_sq, n2_sq, b2_sq = _ux_dot(n1, n1), _ux_dot(n2, n2), _ux_dot(b2, b2)
        t = UX(tol)
        if (not n1_sq.is_positive() or not n2_sq.is_positive() or n1_sq <= t * _ux_dot(b1, b1) * b2_sq
                or n2_sq <= t * b2_sq * _ux_dot(b3, b3)):
            return None
        phi = ux_atan2_final(b2_sq.sqrt() * _ux_dot(b1, n2), _ux_dot(n1, n2))
        return pi if phi <= -pi else phi


def exp(x): return cr("exp", x)
def expm1(x): return cr("expm1", x)
def exp10(x): return cr("exp10", x)
def log(x): return cr("log", x)
def log10(x): return cr("log10", x)
def cos(x): return cr("cos", x)
def sin(x): return cr("sin", x)
def acos(x): return cr("acos", x)
def cbrt(x): return cr("cbrt", x)
def pow(x, y): return cr2("pow", x, y)
def atan2(y, x): return cr2("atan2", y, x)


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


if __name__ == "__main__":
    import math
    # Exact ties, rounded to even (glibc 2.44's pow returns 121 here: port-parity.md §8).
    assert bits(pow(3 * 2.0 ** -215, 5)) == 122
    assert bits(pow(9 * 2.0 ** -430, 2.5)) == 122
    assert pow(0.5, 1075) == 0.0 and bits(pow(0.5, 1074)) == 1
    assert exp10(23.0) == float("1e23")            # 10^23 is exactly a midpoint; the parser rounds it to even too
    assert exp10(22.0) == 1e22 and log10(1e22) == 22.0
    # Subnormal and overflow ends.
    assert bits(exp(-745.1332191019411)) == 1 and exp(-745.2) == 0.0
    assert exp(709.78) < math.inf and exp(709.8) == math.inf
    # Signed zeros, specials.
    assert bits(expm1(-0.0)) == bits(-0.0) and log(1.0) == 0.0 and math.isnan(log(-1.0))
    assert log(0.0) == -math.inf and exp(-math.inf) == 0.0
    # Agrees with glibc where glibc is correctly rounded.
    assert exp(1.0) == math.e and log(10.0) == math.log(10.0) and acos(-1.0) == math.pi
    # x·2^y rounded once (ADR 0002 #13): exact scalings and their ties, a subnormal result from a fractional y, overflow.
    assert bits(scaled_exp2(1000.0, -1075.0)) == 500 and scaled_exp2(1.0, -1075.0) == 0.0
    assert bits(scaled_exp2(1.5, -1074.0)) == 2 and bits(scaled_exp2(3.0, -1076.0)) == 1
    assert scaled_exp2(1.0, 0.5) == math.sqrt(2.0) and scaled_exp2(-2.0, -0.5) == -math.sqrt(2.0)
    assert scaled_exp2(fb := from_bits(0x7FEFFFFFFFFFFFFF), 1.0) == math.inf and scaled_exp2(fb, -1.0) == fb / 2
    assert bits(scaled_exp2(1.0, -1074.5)) == 1 and scaled_exp2(1.0, -1075.5) == 0.0   # 2^-0.5 above / below the tie
    # 1/(100·x) rounded once (ADR 0002 #19): exact quotients, one Rust misses, a subnormal and an infinite result, signs.
    from fractions import Fraction
    assert hundredth_div(0.0000005) == 20000.0 and hundredth_div(0.25) == 0.04
    assert bits(hundredth_div(17.0)) == 0x3F434679ACE01346 != bits(0.01 / 17.0)
    assert bits(hundredth_div(from_bits(0x7FEFFFFFFFFFFFFF))) == 0xA3D70A3D70A
    assert hundredth_div(from_bits(1)) == math.inf and bits(hundredth_div(-0.0)) == bits(-math.inf)
    assert bits(hundredth_div(math.inf)) == 0 and bits(hundredth_div(-17.0)) == 0xBF434679ACE01346
    for k in range(1, 2000):
        assert hundredth_div(float(k)) == float(Fraction(1, 100 * k)), k
    # The unbounded exponent (ADR 0002 #20): squares that overflow and underflow binary64 stay exact; the final √ and
    # atan2 round once (a subnormal result, ties to even, overflow), and Rust's own operations where they stay normal.
    big = UX(2.0 ** 1000)
    assert ux_sqrt_final(big * big) == 2.0 ** 1000 and ux_sqrt_final(big * big * UX(2.0 ** 50)) == math.inf
    tiny = UX(2.0 ** -1074)
    assert bits(ux_sqrt_final(tiny * tiny)) == 1 and ux_sqrt_final(tiny * tiny * UX(0.25)) == 0.0   # 2^-1075: a tie
    assert bits(ux_sqrt_final(tiny * tiny * UX(0.25) + tiny * tiny * UX(2.0 ** -52))) == 1   # just above the tie
    assert bits(ux_atan2_final(UX(3.0 * 2.0 ** -774), UX(2.0 ** 300))) == 3           # just below 3·2^-1074
    assert bits(ux_atan2_final(UX(2.0 ** -1000) * UX(2.0 ** -1000), UX(1.0))) == 0
    assert ux_atan2_final(UX(1.0), UX(2.0 ** -1000) * UX(2.0 ** -1000)) == math.pi / 2
    assert (UX(1.0) - UX(1.0)).is_zero() and not (UX(1.0) - UX(1.0)).neg_sign() and (UX(-0.0) * UX(1.0)).neg_sign()
    assert (big * big).is_f64(2.0 ** 1000) is False and (UX(3.0) * UX(5.0)).is_f64(15.0) and UX(-0.0).is_f64(-0.0)
    # The overflow audit's B2: the final subtraction's signed zeros and range ends; √ of a 256-bit quotient.
    assert bits(ux_sub_final(UX(-0.0), UX(0.0))) == bits(-0.0) and bits(ux_sub_final(UX(1.5), UX(1.5))) == 0
    assert ux_sub_final(UX(from_bits(0x7FEFFFFFFFFFFFFF)), UX(-2.0 ** 1000) * UX(2.0 ** 24)) == math.inf
    assert bits(ux_sub_final(UX(0.0), UX(from_bits(1)) * UX(1.5))) == bits(-from_bits(2))       # −1.5·2^-1074: a tie
    assert ux_sqrt_div_final(UX(4.0), UX(1.0)) == 2.0 and ux_sqrt_div_final(UX(2.0), UX(1.0)) == math.sqrt(2.0)
    assert bits(ux_sqrt_div_final(UX(500.0), UX(2.0 ** -1022))) == 0x60265C55827DF1D2       # Rust's (500/2^-1022).sqrt() = ∞
    assert bits(ux_sqrt_div_final(UX(from_bits(1)), UX(from_bits(0x7FEFFFFFFFFFFFFF)))) == 0x2000000
    assert ux_sqrt_div_final(UX(from_bits(0x7FEFFFFFFFFFFFFF)), UX(from_bits(1))) == math.inf
    rng = __import__("random").Random(20)
    for _ in range(2000):
        a, b = rng.uniform(-1e3, 1e3), rng.uniform(-1e3, 1e3)
        assert (UX(a) * UX(b)).is_f64(a * b) and (UX(a) + UX(b)).is_f64(a + b) and (UX(a) - UX(b)).is_f64(a - b)
        assert (UX(a) / UX(b)).is_f64(a / b) and (UX(abs(a)).sqrt()).is_f64(math.sqrt(abs(a)))
        assert ux_sub_final(UX(a), UX(b)) == a - b
        if a > 0 and b > 0:
            q = Fraction(a) / Fraction(b)
            r = ux_sqrt_div_final(UX(a), UX(b))
            assert Fraction(r) ** 2 * (1 - Fraction(1, 2 ** 52)) < q < Fraction(r) ** 2 * (1 + Fraction(1, 2 ** 52))
    print(f"mpfr_ref: MPFR {VERSION}, self-checks ok")
