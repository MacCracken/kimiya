"""Correctly rounded binary64 results from the system MPFR library (ADR 0006).

Each function returns the IEEE-754 binary64 round-to-nearest-even value of the exact mathematical
result: MPFR at 53 bits of precision with binary64's exponent range (emin = -1073, emax = 1024) and
`mpfr_subnormalize`, so subnormal results, overflow to infinity and exact ties (pow(3*2^-215, 5),
exp10(23)) round exactly as IEEE round-to-nearest does. A decimal reference at fixed precision can't
do that: it misjudges exact ties (it rounds pow(0.5, 1075) = 2^-1075 up instead of to 0).

No Python package is needed: MPFR (libmpfr.so.6, MPFR 4.x) is called through ctypes. Run
`python3 -I scripts/oracle/mpfr_ref.py` for the self-checks.
"""
import ctypes
import ctypes.util
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
    print(f"mpfr_ref: MPFR {VERSION}, self-checks ok")
