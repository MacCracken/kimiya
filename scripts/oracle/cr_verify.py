"""Verify the correctly rounded exp, exp10, expm1, ln, log10, pow, sin, cos, atan2, acos and x·2^y of src/mathx.cyr
against MPFR (ADR 0006).

For each primitive it generates seeded arguments (deterministic: random.Random with fixed seeds), builds
a probe program that calls the Cyrius function on every argument, runs it on each target, and compares
every result bit for bit with the correctly rounded MPFR value (scripts/oracle/mpfr_ref.py). Each
primitive gets:
- >= 1,000,000 arguments: binade-uniform over the whole domain (random sign, binary exponent and
  significand), plus the dense regions where the code changes path or the rounding is delicate
  (subnormal results, the overflow and underflow ends, the thresholds +- a few thousand ulps, the
  reduction boundaries, the tiny-argument ends, ln/log10 near 1 and at powers of 10, expm1 near 0 and at
  multiples of ln 2), every special and exact case, and the known hard-to-round arguments;
- the accurate path forced (the _km_cr_<f>_acc entry point) on >= 100,000 of them.
pow gets >= 1,000,000 (x, y) pairs: x over every binade with y in [-1100, 1100], and with |y ln x| <= 745;
results near overflow, near underflow and subnormal; x near 1 with huge y; subnormal x; negative x with
integer y; tiny |y|; the exact and midpoint families (x = m^n 2^e, perfect squares, cubes, 4th and 8th
powers with dyadic y, powers of two around the tie 2^-1075); the published hard cases; every C99 special
pair. It also checks pow(10, x) against exp10(x), _km_cr_pow(x, 2.5) (Redlich–Kwong's Tc^2.5; ADR 0005's
double-double _km_pow_2_5 until M3c-3) against MPFR over every binade, and the special values against glibc's pow
(ctypes, libm.so.6) as well as MPFR.
sin and cos (M4a) get >= 1,000,000 arguments each: every binade up to 2^1023, the Peng–Robinson range [-2.2, 3.2],
doubles near multiples of pi/2 (RN(k pi/2) +- ulps for k < 2^40, and in every binade the doubles closest to a
multiple of pi/2, found from the continued fraction of 2^(e-52)*2/pi, with Lefevre and Muller's worst case for
the reduction), the FFT angles -2pi/n, the Cody–Waite / Payne–Hanek edge, the table edges, tiny |x|, and
Lefevre–Muller's and Hanrot et al.'s hard cases; _km_cr_sincos is checked against both on the sin arguments.
atan2 gets >= 1,000,000 (y, x) pairs: every quadrant over every binade, every exponent difference, results near
+-pi and +-pi/2, subnormal and exactly-midpoint quotients, |y/x| near the table points and 1, bond and dihedral
angles at Angstrom and metre scales (near-collinear and perpendicular included), arctan's hard cases, and every
C99 special pair (also against glibc's atan2). acos gets >= 1,000,000: [-1, 1], dense near +-1 and 0, every
binade, the branch points 1/2 and 1/sqrt(2), and arccos's and arcsin's hard cases. (Until M4d it also compared M2a's
double-double _km_acos and _km_cos with them: bit-identical on every argument; M4d removed the two.)
exp2s, _km_cr_exp2_scaled(x, y) = x*2^y rounded once (M4d, ADR 0002 #13: nuclear's N0/2^n past the normal range), gets
>= 600,000 (x, y) pairs: nuclear's edge (N0 = m*2^e, y = -n with n beyond 1022 or -1024), results across the subnormal
range and at the overflow end, exact scalings and their ties (integer y), tiny and subnormal |y|, |y| near 2200, every
binade of x with y over [-2300, 2300], and the special pairs; the reference is mpfr_ref.scaled_exp2 (2^f to 400 bits,
bounded from both sides).
Through each function's _km_cr_<f>_fast entry it reports the fast path's fallback rate per region, and
checks the fast path's result equals the function's wherever its test passes. It checks that the three
targets give identical bits and flags.

Usage, from the repo root:
    python3 -I scripts/oracle/cr_verify.py                 # x86_64 ELF, x86_64 PE (wine), aarch64 (qemu)
    python3 -I scripts/oracle/cr_verify.py --targets elf   # one target
    python3 -I scripts/oracle/cr_verify.py --big           # also 10,000,000 more per primitive on ELF
    python3 -I scripts/oracle/cr_verify.py --only pow,exp  # some functions (exp exp10 expm1 ln log10 pow sin cos
                                                           # atan2 acos exp2s)
    python3 -I scripts/oracle/cr_verify.py --work DIR      # probe sources and binaries in DIR
                                                           # (default build/cr-verify, gitignored)
Exit status 1 on any mismatch or cross-target difference, 2 on a build or run failure.
"""
import math
import random
import struct
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True   # no __pycache__ in scripts/oracle/
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpfr_ref  # noqa: E402

FUNCS = [  # (id, name, fast entry, accurate entry, MPFR reference)
    (0, "exp", "_km_cr_exp", "_km_cr_exp_acc", mpfr_ref.exp),
    (1, "exp10", "_km_cr_exp10", "_km_cr_exp10_acc", mpfr_ref.exp10),
    (2, "expm1", "_km_cr_expm1", "_km_cr_expm1_acc", mpfr_ref.expm1),
    (3, "ln", "_km_cr_ln", "_km_cr_ln_acc", mpfr_ref.log),
    (4, "log10", "_km_cr_log10", "_km_cr_log10_acc", mpfr_ref.log10),
    (5, "pow", "_km_cr_pow", "_km_cr_pow_acc", mpfr_ref.pow),
    (6, "sin", "_km_cr_sin", "_km_cr_sin_acc", mpfr_ref.sin),
    (7, "cos", "_km_cr_cos", "_km_cr_cos_acc", mpfr_ref.cos),
    (8, "atan2", "_km_cr_atan2", "_km_cr_atan2_acc", mpfr_ref.atan2),
    (9, "acos", "_km_cr_acos", "_km_cr_acos_acc", mpfr_ref.acos),
    (11, "exp2s", "_km_cr_exp2_scaled", "_km_cr_exp2_scaled_acc", None),   # the reference is ref_exp2s, below
]
PAIRS = (5, 8, 11)   # the functions of two arguments
FID_SINCOS = 10      # probe-only id: _km_cr_sincos
N_UNIFORM = 700_000
N_ACC = 100_000

fb = mpfr_ref.from_bits
bits = mpfr_ref.bits
INF = float("inf")
NAN = float("nan")


def ulps(x, k):
    """x moved k ulps along the ordered bit line (k may be negative), staying finite."""
    b = bits(x)
    o = b if b < (1 << 63) else -(b & ((1 << 63) - 1))
    o += k
    nb = o if o >= 0 else ((-o) | (1 << 63))
    y = fb(nb)
    return y if math.isfinite(y) else x


def binade_uniform(rng, n, emin, emax, signs=(1, -1)):
    """n doubles with a uniformly random binade 2^e (emin <= e <= emax; the subnormal binades too) and a
    uniformly random significand in it."""
    out = []
    for _ in range(n):
        e = rng.randint(emin, emax)
        if e < -1022:   # a subnormal binade: the leading bit at 2^e
            k = e + 1074
            x = fb((1 << k) | rng.randrange(0, 1 << k))
        else:
            x = fb(((e + 1023) << 52) | rng.randrange(0, 1 << 52))
        out.append(rng.choice(signs) * x)
    return out


def around(x, k, rng=None, n=None):
    """x and k ulps either side (all of them), or n random offsets within +-k."""
    if n is None:
        return [ulps(x, d) for d in range(-k, k + 1)]
    return [ulps(x, rng.randint(-k, k)) for _ in range(n)]


SPECIALS = [0.0, -0.0, INF, -INF, NAN, fb(1), -fb(1), fb(0x0010000000000000), -fb(0x0010000000000000),
            fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 1.0, -1.0, 2.0, 0.5, 10.0, 0.1, -0.1,
            fb(0x000FFFFFFFFFFFFF), 2.0 ** -54, -2.0 ** -54, 2.0 ** -53, -2.0 ** -53, 2.0 ** -30, -2.0 ** -30,
            2.0 ** -6, -2.0 ** -6, 708.0, -708.0, 307.0, -307.0]


def _b2f(mant, e, neg=False):
    ip, fp = mant.split(".")
    from fractions import Fraction
    v = Fraction(int(ip + fp, 2), 2 ** len(fp)) * Fraction(2) ** e
    x = float(v)
    return -x if neg else x


# Lefèvre–Muller 2003 worst-case database (Tables 4 and 5): exp and ln hardest-to-round arguments.
HARD_EXP = [
    _b2f("1.1110110100110001100011101111101101100010011111101010", -27, True),
    _b2f("1.0100000000000000000000000000000000000000000000110010", -46, True),
    _b2f("1.0000000000000000000000000000000000000000000000000001", -51, True),
    _b2f("1.1111111111111111111111111111111111111111111111111111", -53),
    _b2f("1.1111111111111111111111111111111111111111111110000000", -46),
    _b2f("1.0001111111111111111111111111111111111111111110101111", -45),
    _b2f("110.00001111010100101111001101111010111011001111110100", 0),
]
HARD_LN = [
    _b2f("1.1110101001110001110110000101110011101110000000100000", -509),
    _b2f("1.1001010001110110111000110000010011001101011111000111", -384),
    _b2f("1.0010011011101001110001001101001100100111100101100000", -232),
    _b2f("1.0110000100111001010101011101110010000000001011111000", -35),
    _b2f("1.0110001010101000100001100001001101100010100110110110", 678),
]
LN2 = 0.6931471805599453
LN10 = 2.302585092994046


def args_exp(rng, lim):
    ovf, unf = lim["EXP_OVF"], lim["EXP_UNF"]
    red = []
    for _ in range(30_000):  # reduction boundaries: (K + 1/2)·ln2/512 (K's rounding flips), K·ln2/512 (r ≈ 0)
        k = rng.randint(-550_000, 524_000)
        red.append(ulps((k + 0.5) * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-700,700]", [rng.uniform(-700.0, 700.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("subnormal results", [rng.uniform(unf - 0.01, -708.0) for _ in range(120_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(40_000)]),
        ("thresholds+-ulps", around(ovf, 2000) + around(unf, 2000) + around(708.0, 500) + around(-708.0, 500)
         + around(-708.3964185322641, 2000)),                               # the last: results near DBL_MIN
        ("tiny |x|", [math.copysign(2.0 ** rng.uniform(-60, -25), rng.random() - 0.5) for _ in range(60_000)]
         + around(2.0 ** -54, 300) + around(-2.0 ** -54, 300) + around(2.0 ** -30, 300) + around(-2.0 ** -30, 300)),
        ("reduction edges", red),
    ]


def args_exp10(rng, lim):
    ovf, unf = lim["EXP10_OVF"], lim["EXP10_UNF"]
    red = []
    for _ in range(30_000):
        k = rng.randint(-550_000, 524_000)
        red.append(ulps((k + 0.5) * LN2 / 512 / LN10, rng.randint(-3, 3)))
        red.append(ulps(k * LN2 / 512 / LN10, rng.randint(-3, 3)))
    return [
        ("special+integers", list(SPECIALS) + [float(n) for n in range(-330, 331)]
         + [ulps(float(n), d) for n in range(-30, 31) for d in (-3, -2, -1, 1, 2, 3)]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 8)),
        ("typical [-300,300]", [rng.uniform(-300.0, 300.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("subnormal results", [rng.uniform(unf - 0.01, -307.65) for _ in range(120_000)]),
        ("overflow end", [rng.uniform(300.0, ovf + 0.01) for _ in range(40_000)]),
        ("thresholds+-ulps", around(ovf, 2000) + around(unf, 2000) + around(307.0, 500) + around(-307.0, 500)
         + around(-307.6526555685888, 2000)),
        ("tiny |x|", [math.copysign(2.0 ** rng.uniform(-62, -25), rng.random() - 0.5) for _ in range(60_000)]
         + around(2.0 ** -56, 300) + around(-2.0 ** -56, 300)),
        ("reduction edges", red),
    ]


def args_expm1(rng, lim):
    ovf, m1 = lim["EXPM1_OVF"], lim["EXPM1_M1"]
    red = []
    for _ in range(40_000):     # ± multiples of ln 2 / 512 and of ln 2 (j = 0, r ≈ 0; E changes)
        k = rng.randint(-30_000, 30_000)
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(rng.randint(-60, 60) * LN2, rng.randint(-200, 200)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-40,40]", [rng.uniform(-40.0, 40.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("near 0 (|x|<2^-6)", [math.copysign(2.0 ** rng.uniform(-60, -6), rng.random() - 0.5) for _ in range(150_000)]
         + around(2.0 ** -6, 2000) + around(-2.0 ** -6, 2000) + around(2.0 ** -53, 500) + around(-2.0 ** -53, 500)),
        ("the -1 end", around(m1, 3000) + [rng.uniform(-40.0, -30.0) for _ in range(30_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(20_000)] + around(ovf, 1000) + around(708.0, 500)),
        ("k·ln2/512, k·ln2", red),
    ]


def args_log(rng, lim, log10=False):
    edges = []
    for _ in range(40_000):     # interval edges of the 512-entry table, in random binades
        j = rng.randrange(0, 512)
        e = rng.randint(-1021, 1022)
        b = 0x3FE7FC0000000000 + j * (1 << 43) + rng.randint(-4, 4)
        b = (b & 0x000FFFFFFFFFFFFF) | (((b >> 52) + e) << 52)
        if 0 < b < 0x7FF0000000000000:
            edges.append(fb(b))
    return [
        ("special+hard", list(SPECIALS) + HARD_LN + [-1.0, -INF, -2.5]),
        ("powers of 10", [10.0 ** k for k in range(-30, 31)]
         + [ulps(10.0 ** k, d) for k in range(-25, 26) for d in range(-50, 51) if d]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 1023, signs=(1,))),
        ("negative", binade_uniform(rng, 20_000, -1074, 1023, signs=(-1,))),
        ("near 1", [rng.uniform(1 - 2.0 ** -6, 1 + 2.0 ** -6) for _ in range(150_000)] + around(1.0, 3000)),
        ("[0.5,2]", [rng.uniform(0.5, 2.0) for _ in range(50_000)]),
        ("subnormal x", [fb(rng.randrange(1, 1 << 52)) for _ in range(30_000)]),
        ("table edges", edges),
    ]


# Published hard pow cases (Hubrecht, Jeannerod and Zimmermann, "Towards a correctly-rounded and fast power
# function in binary64 arithmetic", ARITH 2023, slides 2, 5 and 6): libraries that misrounded them.
HARD_POW = [
    (float.fromhex("0x1.0fbcd29b56829p-1"), float.fromhex("0x1.57347b643eac2p-1")),
    (float.fromhex("0x1.30b3e414e3d3bp-1"), float.fromhex("0x1.a6c0a38da8066p-1")),
    (float.fromhex("0x1.470574d68e0afp+1"), float.fromhex("0x1.02e0706205c0ep+1")),
    (float.fromhex("0x1.f80b060553772p-1"), float.fromhex("0x1.99cp+13")),      # 49 identical bits
    (float.fromhex("0x1.524ebae943097p+1"), float.fromhex("0x1.ep-2")),         # 64 identical bits
]


def pow_specials():
    """Every C99 Annex F special pair, over the special values of x and y."""
    xs = [0.0, -0.0, 1.0, -1.0, INF, -INF, NAN, 0.5, -0.5, 2.0, -2.0, fb(1), -fb(1), 1.5, -1.5,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 10.0, -10.0]
    ys = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5, 2.5, -2.5, INF, -INF, NAN,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 2.0 ** 53, 2.0 ** 53 + 2, 2.0 ** 53 - 1,
          -(2.0 ** 53 - 1), 2.0 ** 64, 2.0 ** 63 + 2 ** 11, fb(1), -fb(1), 1e-300, 1075.0, -1075.0]
    return [(x, y) for x in xs for y in ys]


def args_pow(rng, lim):
    """(x, y) pairs over the regions that matter, by section."""
    def t_over(lo, hi, n, xgen):
        out = []
        for _ in range(n):
            x = xgen()
            l = math.log(x)
            if l == 0.0:
                continue
            out.append((x, rng.uniform(lo, hi) / l))
        return out
    pos = lambda: binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0]
    near1 = lambda: 1.0 + math.copysign(2.0 ** rng.uniform(-52, -1), rng.random() - 0.5)
    exact = []
    for m in (1, 3, 5, 7, 9, 11, 13, 15, 17, 21, 25, 27, 49, 81, 125, 243, 625, 2187, 3**10, 3**20, 3**33, 5**22, 7**18):
        for e in range(-1100, 1050, 37):
            x = math.ldexp(float(m), e) if -1200 < e + m.bit_length() < 1024 else 0.0
            if not (0 < x < INF) or x != x:
                continue
            for n in range(1, 41):
                exact.append((x, float(n)))
                exact.append((x, -float(n)))
    for a in (1, 3, 5, 7, 9, 11, 13, 15, 27, 31, 101, 1001, 3**6, 65535, 2**26 - 1):
        for e in range(-540, 520, 23):
            for k, root in ((1, 2), (2, 4), (3, 8), (1, 2)):
                if a ** root >= 2 ** 53 or not (-1074 < root * e + (a ** root).bit_length() < 1024):
                    continue
                x = math.ldexp(float(a ** root), root * e)
                if not (0 < x < INF):
                    continue
                for n in range(1, 36, 2):
                    exact.append((x, n / 2 ** k))
                    exact.append((x, -n / 2 ** k))
    for e in range(-1080, 1024):                       # powers of two: exact, the tie 2^-1075, out of range
        for y in (0.5, 1.5, 2.0, 3.0, -0.5, -1.0, 0.25, 1.25):
            exact.append((math.ldexp(1.0, max(e, -1074)), y))
    exact += [(3 * 2.0 ** -215, 5.0), (9 * 2.0 ** -430, 2.5), (0.5, 1075.0), (0.5, 1074.0), (2.0, -1075.0),
              (4.0, -537.5), (0.25, 537.5), (10.0, 23.0), (10.0, 22.0), (27.0, float.fromhex("0x1.5555555555555p-2"))]
    sq = []                                            # perfect squares and cubes with y = 0.5, 1.5, 2.5, 1/3
    for _ in range(20_000):
        a = rng.randrange(1, 1 << 26)
        e = rng.randint(-330, 330)
        xs = math.ldexp(float(a * a), 2 * e)
        xc = math.ldexp(float(a ** 3), 3 * e) if a < (1 << 17) else xs
        for y in (0.5, 1.5, 2.5, -0.5, float.fromhex("0x1.5555555555555p-2")):
            sq.append((xs, y))
            sq.append((xc, y))
    negint = []
    for _ in range(60_000):
        x = -binade_uniform(rng, 1, -60, 60, signs=(1,))[0]
        n = rng.choice([rng.randint(-300, 300), rng.randint(-40, 40), 2 ** 53 + 2 * rng.randint(0, 99),
                        2 ** 53 - 1 - 2 * rng.randint(0, 99), rng.randint(-2000, 2000)])
        negint.append((x, float(n)))
    return [
        ("special", pow_specials() + HARD_POW),
        ("x binade-uniform, y in [-1100,1100]", [(pos(), rng.uniform(-1100.0, 1100.0)) for _ in range(250_000)]),
        ("x binade-uniform, |y ln x| <= 745", t_over(-745.0, 709.0, 200_000, pos)),
        ("typical x in (0,100], y in [-10,10]", [(rng.uniform(1e-3, 100.0), rng.uniform(-10.0, 10.0)) for _ in range(50_000)]),
        ("subnormal results", t_over(-745.2, -708.3, 80_000, pos)),
        ("overflow end", t_over(700.0, 709.8, 60_000, pos)),
        ("x near 1, huge y", t_over(-745.0, 709.0, 120_000, near1)),
        ("x subnormal", [(fb(rng.randrange(1, 1 << 52)), rng.uniform(-1.5, 1.5)) for _ in range(30_000)]),
        ("negative x, integer y", negint),
        ("exact and midpoint families", exact),
        ("squares/cubes, y = 1/2, 3/2, 5/2, 1/3", sq),
        ("tiny |y|", [(pos(), math.copysign(2.0 ** rng.uniform(-80, -20), rng.random() - 0.5)) for _ in range(30_000)]),
        ("y = 2.5 (x^2.5)", [(binade_uniform(rng, 1, -500, 409, signs=(1,))[0], 2.5) for _ in range(40_000)]),
    ]


# Lefevre–Muller 2003 worst cases (Tables 8–11 and 13; the run after the round bit in the comment), checked with
# 400-bit MPFR: sin on [2^-24, 2 + 4675/8192], arcsin, cos on [2^-25, 12867/8192], arccos, arctan.
HARD_SIN = [float.fromhex(h) for h in (
    "0x1.e0000000001c2p-20",   # 0 0^72
    "0x1.598bae9e632f6p-7",    # 0 1^59
    "0x1.fe767739d0f6dp-2",    # 1 1^65
    "0x1.921fb54442d18p+0",    # 1 1^54 (RN(pi/2))
    "0x1.6756745770a51p+1")]   # 0 0^54
HARD_ASIN = [float.fromhex(h) for h in ("0x1.dfffffffff02ep-20", "0x1.dffffffffc0b8p-19", "0x1.e9950730c4696p-2")]
HARD_COS = [float.fromhex(h) for h in (
    "0x1.8000000000009p-23",   # 0 0^88
    "0x1.8000000000024p-22",   # 0 0^82
    "0x1.2000000000f30p-18",   # 0 0^60
    "0x1.06b505550e6b2p-9",    # 0 0^58
    "0x1.97ccd3d2c438fp-6",    # 1 1^55
    "0x1.6b8a6273d7c21p+0")]   # 1 0^54
HARD_ACOS = [float.fromhex("0x1.fd737be914578p-11")]   # 1 1^62
HARD_ATAN = [float.fromhex(h) for h in (
    "0x1.e000000000546p-21", "0x1.22e8d75e2bc7fp-11", "0x1.a933fe176b375p-3", "0x1.7ba49f739829fp-1",
    "0x1.06b8e2c3cc000p+52")]  # the last: 4621447055448553/2^11, arctan's worst case above 2.25e12 (1 0^45)
# Hanrot, Lefevre, Stehle and Zimmermann (ARITH 18, 2007): sin in [2^1023, 2^1024) (1^43 and 0^45 after the 53
# bits), and Lefevre and Muller's (2003, sect. 2.3) double closest to a nonzero multiple of pi/2 (4.687e-19 away).
HARD_SC_BIG = [4621478864517314 * 2.0 ** 971, 5501214608935005 * 2.0 ** 971, 6381956970095103 * 2.0 ** 797]
PI = math.pi


def _twopi_int():
    """2/pi to 1,344 bits as an integer T (2/pi ~ T/2^1344), from the generator's table."""
    import mathx_tables
    t = 0
    for w in mathx_tables.twopi_words()[2:]:
        t = (t << 32) | w
    return t, 32 * (mathx_tables.N_TWOPI - 2)


def closest_to_pio2_multiples():
    """In every binade [2^e, 2^(e+1)), e = 0 ... 1023, the doubles m*2^(e-52) closest to a multiple of pi/2: m among
    the small multiples of the continued-fraction convergent denominators q_n of alpha = 2^(e-52)*2/pi with
    2^52 <= m < 2^53 (|m alpha - round(m alpha)| is then within a few times its minimum over the binade)."""
    t, nb = _twopi_int()
    out = []
    for e in range(0, 1024):
        num, den = t, 1 << (nb - (e - 52))          # alpha = num/den (den may be fractional for big e: rescale)
        if nb - (e - 52) < 0:
            break
        # continued fraction of num/den
        p0, q0, p1, q1 = 0, 1, 1, 0
        a, b = num, den
        conv = []
        while b and q1 < (1 << 53):
            c = a // b
            a, b = b, a - c * b
            p0, q0, p1, q1 = p1, q1, c * p1 + p0, c * q1 + q0
            conv.append(q1)
        for q in conv:
            if q < (1 << 40) or q >= (1 << 53):
                continue
            tlo = -(-(1 << 52) // q)
            for tt in range(tlo, tlo + 4):
                m = tt * q
                if (1 << 52) <= m < (1 << 53):
                    out.append(math.ldexp(float(m), e - 52))
    return out


def args_sc(rng, lim):
    tiny_edges = [2.0 ** -26, 2.0 ** -27, 2.0 ** -25]
    sp = list(SPECIALS) + [ulps(v, d) for v in tiny_edges for d in range(-3, 4)] + \
        [k * PI / 2 for k in range(-40, 41)] + [ulps(k * PI / 2, d) for k in range(1, 9) for d in (-2, -1, 1, 2)] + \
        [2 * PI, -2 * PI, PI, -PI, 2.0 ** 20, ulps(2.0 ** 20, -1), 2.0 ** 1023, fb(0x7FEFFFFFFFFFFFFF)]
    hard = HARD_SIN + HARD_COS + HARD_SC_BIG + [ulps(h, d) for h in HARD_SC_BIG for d in (-2, -1, 1, 2)]
    hard += [-h for h in hard]
    near = []
    import mathx_tables
    pio2 = mathx_tables._pio2()
    for _ in range(25_000):                        # RN(k pi/2) +- ulps, k < 2^40
        k = rng.choice([rng.randint(1, 64), rng.randint(1, 1 << 12), rng.randint(1, 1 << 20), rng.randint(1, 1 << 40)])
        v = mathx_tables._V()
        mathx_tables._m.mpfr_mul_ui(v.p, pio2.p, k, 0)
        x = mathx_tables._get(v)
        for d in (rng.randint(-40, 40), rng.randint(-3, 3), 0, rng.randint(-2000, 2000)):
            near.append(math.copysign(ulps(x, d), rng.random() - 0.5))
    closest = closest_to_pio2_multiples()
    closest = [ulps(x, d) for x in closest for d in (-1, 0, 1)]
    closest += [-x for x in closest]
    fft = [-(2.0 * PI) / (1 << k) for k in range(1, 21)] + [-(2.0 * PI) / n for n in range(2, 4097)] + \
        [-(2.0 * PI) / rng.randint(2, 1 << 20) for _ in range(40_000)]
    edges = []
    for _ in range(15_000):                        # r near (i + 1/2)/256: the table index's rounding
        i = rng.randint(-201, 201)
        k = rng.choice([0, 0, 1, 2, 3, rng.randint(1, 1 << 19)])
        x0 = k * PI / 2 + (i + 0.5) / 256
        edges += [ulps(x0, rng.randint(-3, 3)), ulps(k * PI / 2 + i / 256, rng.randint(-3, 3))]
    return [
        ("special+hard", sp + hard),
        ("binade-uniform", binade_uniform(rng, 470_000, -1074, 1023)),
        ("PR range [-2.2, 3.2]", [rng.uniform(-2.2, 3.2) for _ in range(100_000)]),
        ("[-8, 8]", [rng.uniform(-8.0, 8.0) for _ in range(40_000)]),
        ("near k*pi/2 (k < 2^40)", near),
        ("closest to k*pi/2 per binade", closest),
        ("FFT angles -2pi/n", fft),
        ("huge [2^20, 2^1024)", binade_uniform(rng, 90_000, 20, 1023)),
        ("[0, 2^20] and the 2^20 edge", [rng.uniform(-2.0 ** 20, 2.0 ** 20) for _ in range(40_000)]
         + around(2.0 ** 20, 2000) + around(-2.0 ** 20, 2000)),
        ("tiny |x|", [math.copysign(2.0 ** rng.uniform(-60, -18), rng.random() - 0.5) for _ in range(40_000)]
         + [ulps(v, d) for v in tiny_edges for d in range(-500, 501)]),
        ("table edges (i + 1/2)/256", edges),
    ]


def atan2_specials():
    vals = [0.0, -0.0, fb(1), -fb(1), 1.0, -1.0, 2.0, -0.5, fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF),
            INF, -INF, NAN, 3.0, -3.0, 2.0 ** -1022, -(2.0 ** -1022), 1e300, -1e-300]
    return [(y, x) for y in vals for x in vals]


def _vec(rng, scale):
    return [rng.uniform(-1.0, 1.0) * scale for _ in range(3)]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def bond_and_dihedral(rng, n):
    """atan2(|ba x bc|, ba.bc) (potential::bond_angle) and atan2(|b2| b1.n2, n1.n2) (dihedral_angle), as the
    oracle computes them in doubles, at Angstrom (1) and metre (1e-10) scales, with near-collinear and
    perpendicular geometries."""
    bond, dih = [], []
    for _ in range(n):
        sc = rng.choice([1.0, 1e-10, 1.5, 1e-9, 3.7e-11])
        ba = _vec(rng, sc)
        bc = _vec(rng, sc)
        shape = rng.random()
        if shape < 0.25:                         # near-collinear (angle near 0 or pi)
            eps = 2.0 ** rng.uniform(-50, -5)
            sgn = rng.choice([1.0, -1.0])
            bc = [sgn * v * rng.uniform(0.5, 2.0) + eps * sc * rng.uniform(-1, 1) for v in ba]
        elif shape < 0.4:                        # near-perpendicular
            c = _cross(ba, _vec(rng, sc))
            nn = math.sqrt(_dot(c, c)) or 1.0
            bc = [v / nn * sc + 2.0 ** rng.uniform(-50, -5) * sc * rng.uniform(-1, 1) for v in c]
        cr = _cross(ba, bc)
        bond.append((math.sqrt(_dot(cr, cr)), _dot(ba, bc)))
        b1, b2, b3 = _vec(rng, sc), _vec(rng, sc), _vec(rng, sc)
        if shape < 0.2:
            b3 = [v + 2.0 ** rng.uniform(-40, -5) * sc * rng.uniform(-1, 1) for v in b1]
        n1 = _cross(b1, b2)
        n2 = _cross(b2, b3)
        dih.append((math.sqrt(_dot(b2, b2)) * _dot(b1, n2), _dot(n1, n2)))
    return bond, dih


def args_atan2(rng, lim):
    sgn = lambda: rng.choice((1.0, -1.0))
    quad = [(sgn() * binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0],
             sgn() * binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0]) for _ in range(300_000)]
    diffs = []
    for _ in range(200_000):                     # every exponent difference, -2097 ... 2097
        d = rng.randint(-2097, 2097)
        lo, hi = max(-1074, -1074 - d), min(1023, 1023 - d)
        e = rng.randint(lo, hi)
        y = binade_uniform(rng, 1, e + d, e + d, signs=(1,))[0]
        x = binade_uniform(rng, 1, e, e, signs=(1,))[0]
        diffs.append((sgn() * y, sgn() * x))
    near_pi = [(sgn() * x * 2.0 ** -rng.uniform(0, 80), -x) for x in
               (binade_uniform(rng, 1, -900, 900, signs=(1,))[0] for _ in range(60_000))]
    near_pio2 = [(sgn() * y, sgn() * y * 2.0 ** -rng.uniform(0, 80)) for y in
                 (binade_uniform(rng, 1, -900, 900, signs=(1,))[0] for _ in range(60_000))]
    sub = []
    for _ in range(60_000):                      # quotients from 2^-1100 to 2^-1015 (x > 0): subnormal, and below
        eq = rng.randint(-1100, -1015)
        ex = rng.randint(-1074 - eq, 1023)
        x = binade_uniform(rng, 1, ex, ex, signs=(1,))[0]
        y = binade_uniform(rng, 1, eq + ex, eq + ex, signs=(1,))[0]
        sub.append((sgn() * y, x))
    for _ in range(20_000):                      # exact quotients, midpoints (2K+1)*2^-1075 among them
        k = rng.randint(1, 1 << rng.randint(1, 52)) | 1
        b = rng.randint(1, 1021)
        a = rng.choice([-1075, -1075, -1076, -1074, -1060, rng.randint(-1100, -1000)])
        y = math.ldexp(float(k), a + b)
        if y != 0.0 and math.isfinite(y):
            sub.append((sgn() * y, math.ldexp(1.0, b)))
            sub.append((sgn() * y * 3.0, math.ldexp(3.0, b)))
    sub += [(3 * 2.0 ** -60, 2.0 ** 1015), (-3 * 2.0 ** -60, 2.0 ** 1015), (5 * 2.0 ** -61, 2.0 ** 1015),
            (2.0 ** -1074, 1.5), (2.0 ** -1074, 2.0 ** 1023), (fb(3), 2.0)]
    jpts = []
    for _ in range(30_000):                      # |y/x| near (j + 1/2)/256, near j/256, near 1
        j = rng.randint(0, 255)
        x = sgn() * binade_uniform(rng, 1, -500, 500, signs=(1,))[0]
        for z in ((j + 0.5) / 256, j / 256 + 2.0 ** -rng.uniform(20, 52), 1.0 - 2.0 ** -rng.uniform(1, 53)):
            jpts.append((sgn() * abs(x) * z, x))
    jpts += [(v, v) for v in binade_uniform(rng, 2000, -1074, 1023)] + \
        [(v, -v) for v in binade_uniform(rng, 2000, -1074, 1023)]
    bond, dih = bond_and_dihedral(rng, 100_000)
    hard = [(h, 1.0) for h in HARD_ATAN] + [(-h, 1.0) for h in HARD_ATAN] + \
        [(h * 2.0 ** k, 2.0 ** k) for h in HARD_ATAN for k in (-1000, -60, -1, 1, 60, 900)] + \
        [(1.0, h) for h in HARD_ATAN] + [(-1.0, -h) for h in HARD_ATAN]
    return [
        ("special", atan2_specials() + hard),
        ("quadrants, binade-uniform", quad),
        ("every exponent difference", diffs),
        ("near +-pi (|y| << |x|, x < 0)", near_pi),
        ("near +-pi/2 (|x| << |y|)", near_pio2),
        ("subnormal and exact quotients", sub),
        ("|y/x| near j/256 and 1, |y| = |x|", jpts),
        ("bond angles", bond),
        ("dihedral angles", dih),
    ]


def args_acos(rng, lim):
    sp = list(SPECIALS) + [ulps(2.0 ** -56, d) for d in range(-3, 4)] + [-ulps(2.0 ** -56, d) for d in range(-3, 4)] + \
        [1.0000000000000002, -1.0000000000000002, 2.0, -2.0] + HARD_ACOS + [-h for h in HARD_ACOS] + \
        [math.sin(h) for h in HARD_ASIN] + [ulps(h, d) for h in HARD_ACOS for d in range(-5, 6)]
    near1 = [1.0 - 2.0 ** -rng.uniform(1, 53) for _ in range(100_000)] + \
        [-1.0 + 2.0 ** -rng.uniform(1, 53) for _ in range(80_000)] + around(1.0, 3000) + around(-1.0, 3000)
    near1 = [v for v in near1 if abs(v) <= 1.0]
    branch = []
    for v in (0.5, 0.7071067811865476, 0.7071067811865475, 0.25, 0.1):
        branch += around(v, 2000) + around(-v, 2000)
    return [
        ("special+hard", sp),
        ("[-1, 1]", [rng.uniform(-1.0, 1.0) for _ in range(330_000)]),
        ("near +-1", near1),
        ("near 0", [math.copysign(2.0 ** -rng.uniform(1, 60), rng.random() - 0.5) for _ in range(150_000)]),
        ("binade-uniform |x| < 1", [v for v in binade_uniform(rng, 220_000, -1074, -1)]),
        ("branch points +-1/2, +-1/sqrt2, ...", branch),
        ("PR-like (-q/2)/sqrt(-p^3/27)", [math.copysign(min(1.0, rng.uniform(0, 1) ** 0.25), rng.random() - 0.5)
                                         for _ in range(100_000)]),
    ]


def ref_exp2s(x, y):
    """x*2^y rounded once, with _km_cr_exp2_scaled's special values (src/mathx.cyr): x NaN, +-0 or +-inf gives x for
    a finite y; y NaN gives NaN; y = +inf gives +-inf (NaN for x = 0), y = -inf gives +-0 (NaN for x = +-inf)."""
    if math.isnan(x):
        return x
    if math.isnan(y):
        return y
    if math.isinf(y):
        if y > 0:
            return NAN if x == 0 else math.copysign(INF, x)
        return NAN if math.isinf(x) else math.copysign(0.0, x)
    if x == 0 or math.isinf(x):
        return x
    return mpfr_ref.scaled_exp2(x, y)


def args_exp2s(rng):
    """(x, y) pairs for x*2^y (M4d)."""
    def sig(neg_ok=True):
        x = fb((0x3FF << 52) | rng.randrange(0, 1 << 52))
        return -x if neg_ok and rng.random() < 0.5 else x

    def target(lo, hi, n):
        out = []
        for _ in range(n):   # x over every binade, y placing x*2^y near 2^t (t uniform in [lo, hi])
            x = binade_uniform(rng, 1, -1074, 1023)[0]
            t = rng.uniform(lo, hi)
            out.append((x, t - math.log2(abs(x))))
        return out
    edge = []
    for _ in range(150_000):    # nuclear's ADR 0002 #13 rows: N0 = m*2^e, n past 1022 or -1024, y = -n
        n = rng.uniform(1022.0, 1200.0) if rng.random() < 0.5 else rng.uniform(-1200.0, -1024.0)
        e = rng.randint(-1074, 1023) if n > 0 else rng.randint(-1074, 0)
        edge.append((math.ldexp(sig(), e) if e >= -1022 else sig() * 2.0 ** e, -n))
    quarters = [(m, -(1022.0 + q / 4)) for m in (fb(0x7FEFFFFFFFFFFFFF), 1e300, 2.0 ** 1000, 1e250, 1.0, 1000.0)
                for q in range(0, 215)]
    ties = []
    for _ in range(30_000):     # exact scalings: an odd integer times 2^y near and in the subnormal range (ties)
        o = 2 * rng.randint(0, 1 << 20) + 1
        ties.append((float(o) * rng.choice((1.0, -1.0)), float(rng.randint(-1100, -1040))))
        x = binade_uniform(rng, 1, -1074, 1023)[0]
        ties.append((x, float(rng.randint(-2150, 2150))))
    small = []
    for _ in range(40_000):     # |y| < 1, tiny and subnormal |y| included, and |y| near the 2200 cut
        x = binade_uniform(rng, 1, -1074, 1023)[0]
        small.append((x, math.copysign(2.0 ** rng.uniform(-1074, 0), rng.random() - 0.5)))
        small.append((sig(), rng.choice((1.0, -1.0)) * ulps(2200.0, rng.randint(-1000, 1000))))
    sp = [0.0, -0.0, INF, -INF, NAN, fb(1), -fb(1), fb(0x7FEFFFFFFFFFFFFF), 1.0, -1.0, 2.0 ** -1022, 1.5]
    spy = [0.0, -0.0, INF, -INF, NAN, 0.5, -0.5, 2200.0, -2200.0, 2200.5, -2200.5, 1074.0, -1074.0, -1075.0, 1024.0,
           -1022.25, 1e-300, -fb(1), 1023.0]
    return [
        ("special", [(x, y) for x in sp for y in spy]),
        ("nuclear edge", edge + quarters),
        ("subnormal results", target(-1080.0, -1018.0, 120_000)),
        ("overflow end", target(1020.0, 1025.0, 50_000)),
        ("exact and ties", ties),
        ("small |y|, |y| ~ 2200", small),
        ("every binade, y in [-2300, 2300]", [(binade_uniform(rng, 1, -1074, 1023)[0], rng.uniform(-2300.0, 2300.0))
                                              for _ in range(150_000)]),
    ]


def gen_args(name, lim):
    rng = random.Random({"exp": 1001, "exp10": 1002, "expm1": 1003, "ln": 1004, "log10": 1005, "pow": 1006,
                         "sin": 1007, "cos": 1008, "atan2": 1009, "acos": 1010, "exp2s": 1011}[name])
    if name == "exp":
        sections = args_exp(rng, lim)
    elif name == "exp10":
        sections = args_exp10(rng, lim)
    elif name == "expm1":
        sections = args_expm1(rng, lim)
    elif name == "pow":
        sections = args_pow(rng, lim)
    elif name in ("sin", "cos"):
        sections = args_sc(rng, lim)
    elif name == "atan2":
        sections = args_atan2(rng, lim)
    elif name == "acos":
        sections = args_acos(rng, lim)
    elif name == "exp2s":
        sections = args_exp2s(rng)
    else:
        sections = args_log(rng, lim, name == "log10")
    xs = [x for _, sec in sections for x in sec]
    spans = []
    i = 0
    for tag, sec in sections:
        spans.append((tag, i, i + len(sec)))
        i += len(sec)
    head = len(sections[0][1])
    acc = xs[:head] + [xs[i] for i in random.Random(77).sample(range(len(xs)), N_ACC)]
    return xs, acc, spans


PROBE = r'''# GENERATED by scripts/oracle/cr_verify.py. Reads job.bin (u64 words: function id, mode, n, then n arguments,
# two words each for pow and atan2) and writes res.bin. Mode 0: for each argument the function's result, the fast
# path's result and one byte, 1 if the fast path's test passed. Mode 1: the accurate path's result alone.
# Id 10 (sincos) gives sin and cos in the two slots, and the byte 1 or 0 as above, 2 if the fast pair differs
# from the full pair, 3 if the pair differs from _km_cr_sin and _km_cr_cos; id 11 is x·2^y (M4d), two words.
include "src/mathx_data.cyr"
include "src/mathx.cyr"

#must_use
fn run(fid, mode, x: f64, y: f64): (f64, f64, i64) {
    if (mode == 1) {
        if (fid == 0) { var a: f64 = _km_cr_exp_acc(x); return (a, a, 0); }
        if (fid == 1) { var b: f64 = _km_cr_exp10_acc(x); return (b, b, 0); }
        if (fid == 2) { var c: f64 = _km_cr_expm1_acc(x); return (c, c, 0); }
        if (fid == 3) { var d: f64 = _km_cr_ln_acc(x); return (d, d, 0); }
        if (fid == 4) { var e: f64 = _km_cr_log10_acc(x); return (e, e, 0); }
        if (fid == 6) { var sa: f64 = _km_cr_sin_acc(x); return (sa, sa, 0); }
        if (fid == 7) { var ca: f64 = _km_cr_cos_acc(x); return (ca, ca, 0); }
        if (fid == 8) { var ta: f64 = _km_cr_atan2_acc(x, y); return (ta, ta, 0); }
        if (fid == 9) { var aa: f64 = _km_cr_acos_acc(x); return (aa, aa, 0); }
        if (fid == 10) { var ps, pc = _km_cr_sincos_acc(x); return (ps, pc, 0); }
        if (fid == 11) { var za: f64 = _km_cr_exp2_scaled_acc(x, y); return (za, za, 0); }
        var f: f64 = _km_cr_pow_acc(x, y);
        return (f, f, 0);
    }
    if (fid == 6) { var m: f64 = _km_cr_sin(x); var mf, mo = _km_cr_sin_fast(x); return (m, mf, mo); }
    if (fid == 7) { var n: f64 = _km_cr_cos(x); var nf, no = _km_cr_cos_fast(x); return (n, nf, no); }
    if (fid == 8) { var o: f64 = _km_cr_atan2(x, y); var of, oo = _km_cr_atan2_fast(x, y); return (o, of, oo); }
    if (fid == 9) { var q: f64 = _km_cr_acos(x); var qf, qo = _km_cr_acos_fast(x); return (q, qf, qo); }
    if (fid == 10) {
        var s, c = _km_cr_sincos(x);
        var fs, fc, fo = _km_cr_sincos_fast(x);
        var flag = fo;
        if (fo == 1 && (load64(&fs) != load64(&s) || load64(&fc) != load64(&c))) { flag = 2; }
        var s1: f64 = _km_cr_sin(x);
        var c1: f64 = _km_cr_cos(x);
        if (load64(&s1) != load64(&s) || load64(&c1) != load64(&c)) { flag = 3; }
        return (s, c, flag);
    }
    if (fid == 11) {
        var z: f64 = _km_cr_exp2_scaled(x, y);
        var zf, zo = _km_cr_exp2_scaled_fast(x, y);
        return (z, zf, zo);
    }
    if (fid == 0) { var g: f64 = _km_cr_exp(x); var gf, go = _km_cr_exp_fast(x); return (g, gf, go); }
    if (fid == 1) { var h: f64 = _km_cr_exp10(x); var hf, ho = _km_cr_exp10_fast(x); return (h, hf, ho); }
    if (fid == 2) { var i: f64 = _km_cr_expm1(x); var if1, io = _km_cr_expm1_fast(x); return (i, if1, io); }
    if (fid == 3) { var j: f64 = _km_cr_ln(x); var jf, jo = _km_cr_ln_fast(x); return (j, jf, jo); }
    if (fid == 4) { var k: f64 = _km_cr_log10(x); var kf, ko = _km_cr_log10_fast(x); return (k, kf, ko); }
    var l: f64 = _km_cr_pow(x, y);
    var lf, lo = _km_cr_pow_fast(x, y);
    return (l, lf, lo);
}

fn main(): i64 {
    alloc_init();
    var cap = 8 * (3 + 2 * NMAX);
    var buf = alloc(cap);
    var got = file_read_all("job.bin", buf, cap);
    if (got < 24) { return 2; }
    var fid = load64(buf);
    var mode = load64(buf + 8);
    var n = load64(buf + 16);
    var ar = 1;
    if (fid == 5 || fid == 8 || fid == 11) { ar = 2; }
    if (got != 8 * (3 + ar * n)) { return 3; }
    var out = alloc(17 * n + 8);
    var i = 0;
    while (i < n) {
        var x: f64 = load64(buf + 24 + 8 * ar * i);
        var y: f64 = 0.0;
        if (ar == 2) { y = load64(buf + 32 + 16 * i); }
        var r, rf, ok = run(fid, mode, x, y);
        store64(out + 8 * i, r);
        store64(out + 8 * n + 8 * i, rf);
        store8(out + 16 * n + i, ok);
        i = i + 1;
    }
    var w = file_write_all("res.bin", out, 17 * n);
    if (w != 17 * n) { return 4; }
    return 0;
}

var r = main();
sys_exit_group(r);
'''

TARGETS = {
    "elf": ([], "probe", []),
    "pe": (["--win"], "probe.exe", ["wine"]),
    "aarch64": (["--aarch64"], "probe-aarch64", ["qemu-aarch64"]),
}


def build(root, work, nmax, targets):
    src = work / "probe.cyr"
    src.write_text(PROBE.replace("NMAX", str(nmax)))
    for t in targets:
        flags, out, _ = TARGETS[t]
        r = subprocess.run(["cyrius", "build", *flags, str(src), str(work / out)], cwd=root,
                           capture_output=True, text=True)
        log = r.stdout + r.stderr
        if r.returncode != 0 or "warning" in log:
            print(log)
            print(f"cr_verify: probe build failed for {t}")
            sys.exit(2)


def run_probe(work, target, fid, mode, args):
    """(results, fast results, fast flags) for the arguments (floats, or (x, y) pairs for pow)."""
    _, out, runner = TARGETS[target]
    n = len(args)
    flat = [bits(v) for a in args for v in (a if isinstance(a, tuple) else (a,))]
    if fid in PAIRS and args and not isinstance(args[0], tuple):
        raise ValueError("pairs expected")
    (work / "job.bin").write_bytes(struct.pack(f"<3Q{len(flat)}Q", fid, mode, n, *flat))
    res = work / "res.bin"
    for attempt in range(2):            # wine has been seen to exit 0 without writing the file, once
        if res.exists():
            res.unlink()
        r = subprocess.run([*runner, str(work / out)], cwd=work, capture_output=True)
        if r.returncode == 0 and res.exists():
            break
    if r.returncode != 0 or not res.exists():
        print(r.stdout.decode(errors="replace")[-2000:], r.stderr.decode(errors="replace")[-2000:])
        print(f"cr_verify: probe failed on {target} (exit {r.returncode})")
        sys.exit(2)
    data = res.read_bytes()
    words = struct.unpack(f"<{2 * n}Q", data[:16 * n])
    return list(words[:n]), list(words[n:]), data[16 * n:17 * n]


def same(a, b):
    """Bit-identical, or both NaN (the results compared against MPFR, whose NaN has no sign)."""
    if a == b:
        return True
    return math.isnan(fb(a)) and math.isnan(fb(b))


def ref_of(name, ref, a):
    if name == "exp2s":
        return ref_exp2s(*a)
    return ref(*a) if isinstance(a, tuple) else ref(a)


def glibc_atan2():
    import ctypes
    lm = ctypes.CDLL("libm.so.6")
    lm.atan2.argtypes = [ctypes.c_double, ctypes.c_double]
    lm.atan2.restype = ctypes.c_double
    return lm.atan2


def m4a_checks(work, targets, plan):
    """sincos against MPFR's sin and cos on every target; atan2's special values against glibc. True if anything
    failed."""
    failed = False
    if "sin" in plan:
        xs = plan["sin"][0]
        ws = [bits(mpfr_ref.sin(x)) for x in xs]
        wc = [bits(mpfr_ref.cos(x)) for x in xs]
        first = None
        for t in targets:
            gs, gc, fl = run_probe(work, t, FID_SINCOS, 0, xs)
            gas, gac, _ = run_probe(work, t, FID_SINCOS, 1, xs[:N_ACC])
            bad = sum(1 for i in range(len(xs)) if not (same(gs[i], ws[i]) and same(gc[i], wc[i])))
            bad_acc = sum(1 for i in range(len(gas)) if not (same(gas[i], ws[i]) and same(gac[i], wc[i])))
            f2 = sum(1 for f in fl if f >= 2)
            slow = sum(1 for f in fl if f == 0)
            xd = 0 if first is None else sum(1 for a, b in zip(gs + gc + list(fl), first) if a != b)
            if first is None:
                first = gs + gc + list(fl)
            print(f"sincos {t:8} {len(xs):9,} args: {bad} mismatches (sin or cos); accurate pair on {len(gas):,}: "
                  f"{bad_acc}; fast/full or sin/cos pair differences {f2}; fell back {slow:,} "
                  f"({slow / len(xs):.2e}); bits or flags differing from {targets[0]}: {xd}")
            failed |= bool(bad or bad_acc or f2 or xd)
    if "atan2" in plan:
        ga = glibc_atan2()
        sp = atan2_specials()
        gs, _, _ = run_probe(work, targets[0], 8, 0, sp)
        n_om = n_og = n_mg = 0
        for i, (y, x) in enumerate(sp):
            m = bits(mpfr_ref.atan2(y, x))
            g = bits(ga(y, x))
            if not same(m, g):
                n_mg += 1
                print(f"  special: MPFR and glibc differ at atan2({y!r}, {x!r}): MPFR {m:016x}, glibc {g:016x}")
            n_om += not same(gs[i], m)
            n_og += not same(gs[i], g)
        print(f"atan2 special values: {len(sp)} pairs; ours vs MPFR {n_om} differ, ours vs glibc {n_og}, "
              f"MPFR vs glibc {n_mg}")
        failed |= bool(n_om or n_og)
    return failed


def glibc_pow():
    import ctypes
    lm = ctypes.CDLL("libm.so.6")
    lm.pow.argtypes = [ctypes.c_double, ctypes.c_double]
    lm.pow.restype = ctypes.c_double
    return lm.pow


def main():
    argv = sys.argv[1:]
    root = Path(__file__).resolve().parents[2]
    work = root / "build/cr-verify"
    targets = ["elf", "pe", "aarch64"]
    big = "--big" in argv
    if "--work" in argv:
        work = Path(argv[argv.index("--work") + 1]).resolve()
    if "--targets" in argv:
        targets = argv[argv.index("--targets") + 1].split(",")
    only = None
    if "--only" in argv:
        only = argv[argv.index("--only") + 1].split(",")
    work.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(root / "scripts/oracle"))
    import mathx_tables
    lim = {k: fb(v) for k, v in mathx_tables.limits().items()}

    funcs = [f for f in FUNCS if only is None or f[1] in only]
    plan = {}
    nmax = 0
    for fid, name, _, _, ref in funcs:
        xs, acc, spans = gen_args(name, lim)
        plan[name] = (xs, acc, spans)
        nmax = max(nmax, len(xs))
    if big:
        nmax = max(nmax, 10_000_000)
    build(root, work, nmax, targets)
    print(f"cr_verify: MPFR {mpfr_ref.VERSION}; targets {', '.join(targets)}; probe in {work}")

    failed = False
    results = {}
    for fid, name, fast, accname, ref in funcs:
        xs, acc, spans = plan[name]
        t0 = time.time()
        want = [bits(ref_of(name, ref, a)) for a in xs]
        want_acc = [bits(ref_of(name, ref, a)) for a in acc]
        first = None
        for t in targets:
            got, gotf, flags = run_probe(work, t, fid, 0, xs)
            got_acc, _, _ = run_probe(work, t, fid, 1, acc)
            bad = [i for i in range(len(xs)) if not same(got[i], want[i])]
            bad_acc = [i for i in range(len(acc)) if not same(got_acc[i], want_acc[i])]
            bad_fast = [i for i in range(len(xs)) if flags[i] == 1 and gotf[i] != got[i]]
            for i in bad[:5]:
                print(f"  MISMATCH {fast}{xs[i]!r} on {t}: got {got[i]:016x}, want {want[i]:016x}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {accname}{acc[i]!r} on {t}: got {got_acc[i]:016x}, want {want_acc[i]:016x}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {fast}_fast{xs[i]!r} on {t}: {gotf[i]:016x} vs {got[i]:016x}")
            slow = len(xs) - sum(flags)
            if first is None:
                first = (got, got_acc, flags)
                xdiff = 0
                if t == targets[0]:
                    results[name] = (got, flags)
            else:
                xdiff = sum(1 for a, b in zip(got, first[0]) if a != b) + \
                    sum(1 for a, b in zip(got_acc, first[1]) if a != b) + \
                    sum(1 for a, b in zip(flags, first[2]) if a != b)
            print(f"{name:6} {t:8} {len(xs):9,} args: {len(bad)} mismatches; accurate path forced on "
                  f"{len(acc):,}: {len(bad_acc)} mismatches; fast path: {len(bad_fast)} differ, fell back "
                  f"{slow:,} ({slow / len(xs):.2e}); bits or flags differing from {targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0]:
                rates = [f"{tag} {(b - a) - sum(flags[a:b])}/{b - a}" for tag, a, b in spans]
                print("         fallbacks by region: " + "; ".join(rates))
        if big and "elf" in targets and name not in ("pow", "atan2", "exp2s"):
            rng = random.Random(9000 + fid)
            if name in ("ln", "log10"):
                ys = binade_uniform(rng, 10_000_000, -1074, 1023, signs=(1,))
            elif name in ("sin", "cos"):
                ys = [rng.uniform(-8.0, 8.0) for _ in range(5_000_000)] + binade_uniform(rng, 5_000_000, -1074, 1023)
            elif name == "acos":
                ys = [rng.uniform(-1.0, 1.0) for _ in range(10_000_000)]
            elif name == "exp10":
                ys = [rng.uniform(lim["EXP10_UNF"] - 1, lim["EXP10_OVF"] + 1) for _ in range(5_000_000)] + \
                    binade_uniform(rng, 5_000_000, -1074, 8)
            else:
                ys = [rng.uniform(-746.0, 710.0) for _ in range(5_000_000)] + binade_uniform(rng, 5_000_000, -1074, 9)
            wy = [bits(ref(y)) for y in ys]
            gy, _, fy = run_probe(work, "elf", fid, 0, ys)
            bad = [i for i in range(len(ys)) if not same(gy[i], wy[i])]
            for i in bad[:5]:
                print(f"  MISMATCH {fast}({ys[i]!r}) big run: got {gy[i]:016x}, want {wy[i]:016x}")
            slow = len(ys) - sum(fy)
            print(f"{name:6} elf-big  {len(ys):9,} args: {len(bad)} mismatches; fell back {slow:,} "
                  f"({slow / len(ys):.2e})")
            failed |= bool(bad)
        if big and "elf" in targets and name == "pow":
            rng = random.Random(9005)
            ys = [(binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0], 0.0) for _ in range(4_000_000)]
            ys = [(x, rng.uniform(-745.0, 709.0) / math.log(x)) for x, _ in ys if math.log(x) != 0.0]
            ys += [(rng.uniform(1e-3, 100.0), rng.uniform(-10.0, 10.0)) for _ in range(2_000_000)]
            wy = [bits(ref(*a)) for a in ys]
            gy, _, fy = run_probe(work, "elf", fid, 0, ys)
            bad = [i for i in range(len(ys)) if not same(gy[i], wy[i])]
            for i in bad[:5]:
                print(f"  MISMATCH {fast}{ys[i]!r} big run: got {gy[i]:016x}, want {wy[i]:016x}")
            slow = len(ys) - sum(fy)
            print(f"{name:6} elf-big  {len(ys):9,} args: {len(bad)} mismatches; fell back {slow:,} "
                  f"({slow / len(ys):.2e})")
            failed |= bool(bad)
        print(f"       [{time.time() - t0:.0f} s]")

    if "pow" in plan:
        # 10^x: pow(10, x) and exp10(x) are the same function: every exp10 argument through the full pow, on ELF.
        if "exp10" in plan:
            xs10 = plan["exp10"][0]
            g10, _, f10 = run_probe(work, targets[0], 5, 0, [(10.0, x) for x in xs10])
            e10 = results["exp10"][0]
            diff = sum(1 for a, b in zip(g10, e10) if not same(a, b))
            print(f"10^x: pow(10, x) vs exp10(x) on {len(xs10):,} arguments: {diff} differ; pow's fast path fell "
                  f"back on {len(xs10) - sum(f10):,}")
            failed |= bool(diff)
        # x^2.5: _km_cr_pow(x, 2.5), Redlich–Kwong's Tc^2.5 (gas.cyr), against MPFR over every binade, around
        # 2^-387 (where ADR 0005's double-double x^2.5, replaced in M3c-3, stopped being correctly rounded) and at
        # the underflow and overflow ends.
        rng = random.Random(2525)
        x25 = binade_uniform(rng, 200_000, -1074, 1023, signs=(1,)) + \
            [rng.uniform(1e-3, 1e3) for _ in range(100_000)] + \
            [9 * 2.0 ** -430, 2.0 ** -387, 2.0 ** -386, ulps(2.0 ** -387, -1), 4.0, 2.0 ** 409.6, 2.0 ** 410]
        c25, _, _ = run_probe(work, targets[0], 5, 0, [(x, 2.5) for x in x25])
        w25 = [bits(mpfr_ref.pow(x, 2.5)) for x in x25]
        cr_bad = sum(1 for i in range(len(x25)) if not same(c25[i], w25[i]))
        print(f"x^2.5 on {len(x25):,} arguments: _km_cr_pow(x, 2.5) vs MPFR {cr_bad} mismatches")
        failed |= bool(cr_bad)
        # Special values: ours, MPFR and glibc (ctypes libm pow, what Rust's f64::powf calls).
        gp = glibc_pow()
        sp = pow_specials() + HARD_POW
        gs, _, _ = run_probe(work, targets[0], 5, 0, sp)
        n_mg = n_om = n_og = 0
        for i, (x, y) in enumerate(sp):
            m = bits(mpfr_ref.pow(x, y))
            g = bits(gp(x, y))
            if not same(m, g):
                n_mg += 1
                print(f"  special: MPFR and glibc differ at pow({x!r}, {y!r}): MPFR {m:016x}, glibc {g:016x}")
            n_om += not same(gs[i], m)
            n_og += not same(gs[i], g)
        print(f"pow special values: {len(sp)} pairs; ours vs MPFR {n_om} differ, ours vs glibc {n_og}, "
              f"MPFR vs glibc {n_mg}")
        failed |= bool(n_om)
    failed |= m4a_checks(work, targets, plan)
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
