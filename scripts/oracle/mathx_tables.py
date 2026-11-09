"""Generate src/mathx_data.cyr: the tables and constants of the correctly rounded exp, expm1, exp10, ln
and log10 in src/mathx.cyr (ADR 0006, M3c-2a).

Every value is computed from MPFR (through scripts/oracle/mpfr_ref.py's ctypes binding, at 400 bits)
or from exact rational arithmetic (Python's Fraction), then split into doubles, each the round-to-
nearest of what is left. So a k-word entry is the value to about 2^(-53k) relative. The output is
deterministic: no randomness, no time stamps, no dict ordering.

Usage, from the repo root:
    python3 -I scripts/oracle/mathx_tables.py            # write src/mathx_data.cyr
    python3 -I scripts/oracle/mathx_tables.py --check    # exit 1 if src/mathx_data.cyr is stale

Besides writing the values, the generator checks the properties the code's error analysis relies on
(asserting, so a table that would break a proof is never written):
- exp: A_j has at most 26 significant bits and A_j + B_j + C_j is 2^(j/512) to 2^-128 relative;
  L1 has at most 33 bits (K·L1 is exact for |K| < 2^20); the hi part of ln 10 splits 26 + 27 bits.
- ln: for each of the 512 intervals of m, r = m·c_j − 1 is exactly representable for every double m
  in the interval and |r| < 2^-9 (a bit-count argument on the interval's end points); c = 1 and
  L = 0 on the interval holding 1; |L_j| ≥ max |r| elsewhere (the Fast2Sum precondition).
"""
from fractions import Fraction
import ctypes
import struct
import sys
from pathlib import Path

sys.dont_write_bytecode = True   # no __pycache__ in scripts/oracle/
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpfr_ref  # noqa: E402  (the binding: emin/emax already set for binary64)

_m = mpfr_ref._m
_P = ctypes.c_void_p
_RNDN = 0
PREC = 400

_m.mpfr_exp2.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_log.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_const_log2.argtypes = [_P, ctypes.c_int]
_m.mpfr_sub.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_add.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_mul.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_div.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_set.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_neg.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_clear.argtypes = [_P]


class _V:
    """An MPFR number at `prec` bits (400 by default)."""

    def __init__(self, prec=PREC):
        self.t = (ctypes.c_byte * 32)()
        _m.mpfr_init2(ctypes.byref(self.t), prec)

    @property
    def p(self):
        return ctypes.byref(self.t)

    def __del__(self):
        _m.mpfr_clear(ctypes.byref(self.t))


def _d(x):
    v = _V()
    _m.mpfr_set_d(v.p, x, _RNDN)
    return v


def _op(fn, a, b=None):
    r = _V()
    if b is None:
        getattr(_m, fn)(r.p, a.p, _RNDN)
    else:
        getattr(_m, fn)(r.p, a.p, b.p, _RNDN)
    return r


def _ln2():
    r = _V()
    _m.mpfr_const_log2(r.p, _RNDN)
    return r


def _get(v):
    return _m.mpfr_get_d(v.p, _RNDN)


def _round_bits(v, bits):
    """v rounded to nearest at `bits` significant bits, as a double (exact for bits <= 53)."""
    r = _V(bits)
    _m.mpfr_set(r.p, v.p, _RNDN)
    return _m.mpfr_get_d(r.p, _RNDN)


def _split(v, n, first_bits=53):
    """n doubles: the first rounded to `first_bits` bits, each next the nearest double to what is left."""
    out = []
    rest = _op("mpfr_set", v)
    for i in range(n):
        h = _round_bits(rest, first_bits) if i == 0 else _get(rest)
        out.append(h)
        rest = _op("mpfr_sub", rest, _d(h))
    return out


def _fsplit(q, n, first_bits=53):
    """The same split of an exact Fraction."""
    out = []
    rest = q
    for i in range(n):
        if i == 0 and first_bits < 53:
            h = _fround_bits(rest, first_bits)
        else:
            h = float(rest)       # Fraction -> float is correctly rounded
        out.append(h)
        rest = rest - Fraction(h)
    return out


def _fround_bits(q, bits):
    """A Fraction rounded to nearest-even at `bits` significant bits."""
    if q == 0:
        return 0.0
    s = -1 if q < 0 else 1
    a = abs(q)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    while Fraction(2) ** e > a:
        e -= 1
    while Fraction(2) ** (e + 1) <= a:
        e += 1
    scale = Fraction(2) ** (bits - 1 - e)
    n = a * scale
    fl = n.numerator // n.denominator
    rem = n - fl
    if rem > Fraction(1, 2) or (rem == Fraction(1, 2) and fl % 2 == 1):
        fl += 1
    return s * float(Fraction(fl) / scale)


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def _sigbits(x):
    """The number of significant bits of a double (0 for 0)."""
    if x == 0.0:
        return 0
    f = Fraction(abs(x))
    n, d = f.numerator, f.denominator
    while n % 2 == 0:
        n //= 2
    return n.bit_length()


def _exact(xs):
    return sum((Fraction(x) for x in xs), Fraction(0))


def _rel(approx, v):
    """|approx - v| / |v| as a float, for an exact Fraction approx and an MPFR v."""
    d = _op("mpfr_sub", _d_frac(approx), v)
    return abs(_get(d) / _get(v))


def _d_frac(q):
    """An MPFR number holding the Fraction q (a sum of a few doubles) exactly, at 400 bits."""
    acc = _d(0.0)
    rest = q
    for _ in range(12):
        if rest == 0:
            break
        h = float(rest)
        acc = _op("mpfr_add", acc, _d(h))
        rest -= Fraction(h)
    assert rest == 0, "Fraction is not a short sum of doubles"
    return acc


# ── exp: 2^(j/512) = A + B + C, the reduction and the polynomials ────────────────────────────────

N_EXP = 512


def exp_table():
    rows = []
    for j in range(N_EXP):
        v = _op("mpfr_exp2", _d(j / N_EXP))
        a, b, c, d = _split(v, 4, first_bits=26)
        assert _sigbits(a) <= 26
        assert _rel(_exact([a, b, c]), v) < 2.0 ** -128, j
        assert _rel(_exact([a, b, c, d]), v) < 2.0 ** -180, j
        if j == 0:
            assert (a, b, c, d) == (1.0, 0.0, 0.0, 0.0)
        rows.append((a, b, c, d))
    return rows


def exp_consts():
    ln2 = _ln2()
    ln2n = _op("mpfr_div", ln2, _d(float(N_EXP)))
    inv = _get(_op("mpfr_div", _d(float(N_EXP)), ln2))
    l1, l2, l3, l4 = _split(ln2n, 4, first_bits=33)
    assert _sigbits(l1) <= 33
    ln10 = _op("mpfr_log", _d(10.0))
    t10 = _split(ln10, 3)
    h1, h2 = _fsplit(Fraction(t10[0]), 2, first_bits=26)
    assert _sigbits(h1) <= 26 and _sigbits(h2) <= 27 and h1 + h2 == t10[0]
    return {
        "INVLN2N": inv, "L1": l1, "L2": l2, "L3": l3, "L4": l4,
        "P3": float(Fraction(1, 6)), "P4": float(Fraction(1, 24)), "P5": float(Fraction(1, 120)),
        "LN10H": t10[0], "LN10M": t10[1], "LN10L": t10[2], "LN10H1": h1, "LN10H2": h2,
    }


def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def expm1_poly():
    """1/3!, 1/4!, ..., 1/9!: the coefficients of expm1's small-|x| fast path, x³·S(x)."""
    return [float(Fraction(1, _fact(k))) for k in range(3, 10)]


def exp_acc():
    """1/k! for k = 2 ... 17 as triple-doubles: the accurate paths' Taylor series."""
    return [_fsplit(Fraction(1, _fact(k)), 3) for k in range(2, 18)]


# ── ln: m·c_j ≈ 1 over 512 intervals of m in [0.7495, 1.4990), −ln c_j as a triple-double ───────

N_LN = 512
LN_OFF = 0x3FE7FC0000000000   # bits(1) − 256.5·2^43: interval 256 is centred (in bits) on 1


def _m_of_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def ln_table():
    rows = []
    for j in range(N_LN):
        lo_b = LN_OFF + j * (1 << 43)
        hi_b = lo_b + (1 << 43) - 1                 # the last double of the interval
        m_lo, m_hi = _m_of_bits(lo_b), _m_of_bits(hi_b)
        if j == 256:
            c = 1.0
        else:
            mid = (Fraction(m_lo) + Fraction(m_hi)) / 2
            inv = 1 / mid
            # c on the grid 2^-9 (c >= 1, m < 1) or 2^-10 (c < 1, m >= 1): r = m·c − 1 then fits 53 bits.
            g = 9 if m_hi < 1.0 else 10
            c = float(Fraction(round(inv * 2 ** g), 2 ** g))
        cf = Fraction(c)
        r_lo = Fraction(m_lo) * cf - 1
        r_hi = Fraction(m_hi) * cf - 1
        rmax = max(abs(r_lo), abs(r_hi))
        assert rmax < Fraction(1, 2 ** 9), (j, float(rmax))
        # Exactness of r = m·c − 1 for every double m of the interval: r is a multiple of
        # ulp(m)·ulp(c), and |r| < 2^-9, so it has at most log2(2^-9 / (ulp(m)·ulp(c))) bits.
        for mm in (m_lo, m_hi):
            if j == 256:
                continue
            um = Fraction(2) ** (_exp2(mm) - 52)
            uc = _ulp_grid(c)
            nbits = (Fraction(1, 2 ** 9) / (um * uc)).numerator.bit_length() - 1
            assert nbits <= 53, (j, nbits)
            assert _exp2(mm) == _exp2(m_lo), j    # the interval stays in one binade
        if j == 256:
            lh, ll, lt = 0.0, 0.0, 0.0
            assert m_lo < 1.0 < m_hi
        else:
            lv = _op("mpfr_neg", _op("mpfr_log", _d(c)))
            lh, ll, lt = _split(lv, 3)
            assert abs(lh) >= float(rmax), j      # Fast2Sum(L, r) precondition
        rows.append((c, lh, ll, lt))
    return rows


def _exp2(x):
    """The binary exponent of a positive double."""
    return ((bits(x) >> 52) & 0x7FF) - 1023


def _ulp_grid(c):
    """The coarsest power of two c is a multiple of."""
    f = Fraction(c)
    g = Fraction(1)
    while (f / g).denominator != 1:
        g /= 2
    while (f / (g * 2)).denominator == 1:
        g *= 2
    return g


def ln_consts():
    ln2 = _ln2()
    h, mid, lo = _split(ln2, 3, first_bits=42)
    assert _sigbits(h) <= 42
    ln10 = _op("mpfr_log", _d(10.0))
    iv = _op("mpfr_div", _d(1.0), ln10)
    ivh, ivm, ivl = _split(iv, 3)
    i1, i2 = _fsplit(Fraction(ivh), 2, first_bits=26)
    assert _sigbits(i1) <= 26 and _sigbits(i2) <= 27 and i1 + i2 == ivh
    return {
        "LN2H": h, "LN2M": mid, "LN2L": lo,
        "C3": float(Fraction(1, 3)), "C4": -0.25, "C5": float(Fraction(1, 5)),
        "C6": float(Fraction(-1, 6)), "C7": float(Fraction(1, 7)), "C8": -0.125,
        "IVH": ivh, "IVM": ivm, "IVL": ivl, "IVH1": i1, "IVH2": i2,
    }


def ln_acc():
    """(−1)^(k+1)/k for k = 2 ... 16 as triple-doubles: the accurate log1p series."""
    return [_fsplit(Fraction((-1) ** (k + 1), k), 3) for k in range(2, 17)]


# ── thresholds, from the correctly rounded MPFR functions ───────────────────────────────────────

def _last_true(pred, lo_b, hi_b):
    """The largest bit pattern b in [lo_b, hi_b] with pred(b), for pred monotone (true then false)."""
    assert pred(lo_b) and not pred(hi_b)
    while hi_b - lo_b > 1:
        mid = (lo_b + hi_b) // 2
        if pred(mid):
            lo_b = mid
        else:
            hi_b = mid
    return lo_b


def limits():
    fb = _m_of_bits
    inf = float("inf")
    # Positive arguments: bit patterns increase with x.
    exp_ovf = _last_true(lambda b: mpfr_ref.exp(fb(b)) < inf, bits(700.0), bits(720.0))
    exp10_ovf = _last_true(lambda b: mpfr_ref.exp10(fb(b)) < inf, bits(300.0), bits(320.0))
    expm1_ovf = _last_true(lambda b: mpfr_ref.expm1(fb(b)) < inf, bits(700.0), bits(720.0))
    # Negative arguments: bit patterns increase with |x|.
    exp_unf = _last_true(lambda b: mpfr_ref.exp(fb(b)) > 0.0, bits(-740.0), bits(-750.0))
    exp10_unf = _last_true(lambda b: mpfr_ref.exp10(fb(b)) > 0.0, bits(-320.0), bits(-330.0))
    expm1_m1 = _last_true(lambda b: mpfr_ref.expm1(fb(b)) > -1.0, bits(-30.0), bits(-40.0))
    return {
        "EXP_OVF": exp_ovf, "EXP_UNF": exp_unf, "EXP10_OVF": exp10_ovf, "EXP10_UNF": exp10_unf,
        "EXPM1_OVF": expm1_ovf, "EXPM1_M1": expm1_m1,
    }


# ── output ───────────────────────────────────────────────────────────────────────────────────────

def _hex(x):
    return f"0x{bits(x):016X}"


def _table(name, words, comment, per_line=4):
    out = [f"# {line}" if line else "#" for line in comment]
    out.append(f"var {name}: i64[{len(words)}] = {{")
    for i in range(0, len(words), per_line):
        chunk = words[i:i + per_line]
        sep = "," if i + per_line < len(words) else ""
        out.append("    " + ", ".join(chunk) + sep)
    out.append("};")
    out.append("")
    return out


def _scalars(prefix, d, order, comment):
    """One i64 global per value (`var _KM_CR_<prefix>_<NAME> = 0x…;`): a scalar global reads in one load,
    a table element in an address computation and a load, which the fast paths feel."""
    lines = [f"# {line}" if line else "#" for line in comment]
    for k in order:
        lines.append(f"var _KM_CR_{prefix}{k} = {_hex(d[k])};   # {d[k]!r}")
    lines.append("")
    return lines


EXP_K_ORDER = ["INVLN2N", "L1", "L2", "L3", "L4", "P3", "P4", "P5",
               "LN10H", "LN10M", "LN10L", "LN10H1", "LN10H2"]
EXPM1_ORDER = ["Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"]
LN_K_ORDER = ["LN2H", "LN2M", "LN2L", "C3", "C4", "C5", "C6", "C7", "C8",
              "IVH", "IVM", "IVL", "IVH1", "IVH2"]
LIM_ORDER = ["EXP_OVF", "EXP_UNF", "EXP10_OVF", "EXP10_UNF", "EXPM1_OVF", "EXPM1_M1"]


def render():
    et = exp_table()
    ek = exp_consts()
    lt = ln_table()
    lk = ln_consts()
    lim = limits()
    L = []
    L += [
        "# kimiya — tables of the correctly rounded exp, expm1, exp10, ln and log10 (src/mathx.cyr)",
        "#",
        f"# GENERATED by scripts/oracle/mathx_tables.py from MPFR {mpfr_ref.VERSION} (400-bit values) and exact",
        "# rationals. Do not edit.",
        "# Regenerate: python3 -I scripts/oracle/mathx_tables.py   (--check verifies this file is current)",
        "#",
        "# Every value is an IEEE-754 binary64 bit pattern, in an i64 table (read with load64) or an i64 global,",
        "# assigned to an f64 local (a hex literal in an f64 initialiser would be an integer). A multi-word entry",
        "# is a value split into doubles, each the nearest double to what the previous ones leave: a + b + c is",
        "# the value to about 2^(-53·k) relative. ADR 0006; the error analysis is in src/mathx.cyr.",
        "",
    ]
    L += _table("_KM_CR_EXP_T", [_hex(v) for row in et for v in row], [
        "2^(j/512) for j = 0 … 511, four words per j at index 4j: A_j (rounded to 26 significant bits,",
        "so A_j·r1 is exact for a 27-bit r1), then B_j, C_j and D_j (each the nearest double to the rest).",
        "The fast path reads A_j + B_j (2^(j/512) to 2^-79), the accurate path all four (to 2^-180).",
        "Row 0 is (1, 0, 0, 0).",
    ], per_line=4)
    L += _scalars("EXP_", ek, EXP_K_ORDER, [
        "exp / exp10 reduction: x = (512·E + j)·ln2/512 + r. INVLN2N = 512/ln 2 rounded; ln2/512 =",
        "L1 + L2 + L3 + L4 (L1 at 33 bits, so K·L1 is exact for |K| < 2^20). P3 … P5: 1/6, 1/24, 1/120,",
        "the fast path's Taylor coefficients (1/2 is exact). ln 10 = LN10H + LN10M + LN10L, and",
        "LN10H = LN10H1 + LN10H2 (26 + 27 bits, Dekker's product of x·LN10H without a split of LN10H).",
    ])
    L += _scalars("EXPM1_", dict(zip(EXPM1_ORDER, expm1_poly())), EXPM1_ORDER, [
        "Q3 … Q9 = 1/3!, …, 1/9!: expm1's small-|x| fast path, x + x²/2 + x³·(Q3 + x·Q4 + … + x⁶·Q9).",
    ])
    L += _table("_KM_CR_EXP_ACC", [_hex(v) for row in exp_acc() for v in row], [
        "1/k! for k = 2 … 17 as triple-doubles (index 3·(k − 2)): the accurate paths' Taylor series of",
        "exp(r) − 1 and expm1(x).",
    ], per_line=3)
    L += _table("_KM_CR_LN_T", [_hex(v) for row in lt for v in row], [
        "ln: x = 2^E·m, m in [0.74951171875, 1.4990234375), j = the 9 bits of bits(x) − 0x3FE7FC0000000000",
        "below the exponent field. Four words per j at index 4j: c_j ≈ 1/m on interval j (a multiple of",
        "2^-9 if c_j > 1, of 2^-10 if c_j < 1, so r = m·c_j − 1 is exact and |r| < 2^-9), then",
        "−ln(c_j) as a triple-double. Interval 256 holds m = 1: c = 1, −ln c = 0.",
    ], per_line=4)
    L += _scalars("LN_", lk, LN_K_ORDER, [
        "ln 2 = LN2H + LN2M + LN2L (LN2H at 42 bits, so E·LN2H is exact for |E| < 2^11). C3 … C8: the",
        "Taylor coefficients 1/3, −1/4, 1/5, −1/6, 1/7, −1/8 of log1p(r) − r + r²/2. 1/ln 10 = IVH + IVM +",
        "IVL, and IVH = IVH1 + IVH2 (26 + 27 bits).",
    ])
    L += _table("_KM_CR_LN_ACC", [_hex(v) for row in ln_acc() for v in row], [
        "(−1)^(k+1)/k for k = 2 … 16 as triple-doubles (index 3·(k − 2)): the accurate log1p series.",
    ], per_line=3)
    lim_d = {k: _m_of_bits(v) for k, v in lim.items()}
    L += _scalars("", lim_d, LIM_ORDER, [
        "Thresholds, each the last double before the result changes (bisection over bit patterns with",
        "MPFR's correctly rounded function): EXP_OVF the largest x with exp(x) finite, EXP_UNF the",
        "most negative x with exp(x) > 0 (below it exp(x) < 2^-1075 rounds to +0), the same for exp10,",
        "EXPM1_OVF the largest x with expm1(x) finite, EXPM1_M1 the most negative x with expm1(x) > −1.",
    ])
    return "\n".join(L).rstrip("\n") + "\n"


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "src/mathx_data.cyr"
    text = render()
    if "--check" in sys.argv[1:]:
        old = out.read_text() if out.exists() else ""
        if old != text:
            print("STALE    src/mathx_data.cyr (run python3 -I scripts/oracle/mathx_tables.py)")
            sys.exit(1)
        print(f"mathx_tables: src/mathx_data.cyr current (MPFR {mpfr_ref.VERSION})")
        return
    out.write_text(text)
    print(f"mathx_tables: wrote src/mathx_data.cyr ({text.count(chr(10))} lines, MPFR {mpfr_ref.VERSION})")


if __name__ == "__main__":
    main()
