"""Verify the correctly rounded exp, exp10, expm1, ln, log10 and pow of src/mathx.cyr against MPFR (ADR 0006).

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
Through each function's _km_cr_<f>_fast entry it reports the fast path's fallback rate per region, and
checks the fast path's result equals the function's wherever its test passes. It checks that the three
targets give identical bits and flags.

Usage, from the repo root:
    python3 -I scripts/oracle/cr_verify.py                 # x86_64 ELF, x86_64 PE (wine), aarch64 (qemu)
    python3 -I scripts/oracle/cr_verify.py --targets elf   # one target
    python3 -I scripts/oracle/cr_verify.py --big           # also 10,000,000 more per primitive on ELF
    python3 -I scripts/oracle/cr_verify.py --only pow,exp  # some functions
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
]
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


def gen_args(name, lim):
    rng = random.Random({"exp": 1001, "exp10": 1002, "expm1": 1003, "ln": 1004, "log10": 1005, "pow": 1006}[name])
    if name == "exp":
        sections = args_exp(rng, lim)
    elif name == "exp10":
        sections = args_exp10(rng, lim)
    elif name == "expm1":
        sections = args_expm1(rng, lim)
    elif name == "pow":
        sections = args_pow(rng, lim)
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
# two words each for pow) and writes res.bin. Mode 0: for each argument the function's result, the fast path's
# result and one byte, 1 if the fast path's test passed. Mode 1: the accurate path's result alone.
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
        var f: f64 = _km_cr_pow_acc(x, y);
        return (f, f, 0);
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
    if (fid == 5) { ar = 2; }
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
    return ref(*a) if isinstance(a, tuple) else ref(a)


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
        if big and "elf" in targets and name != "pow":
            rng = random.Random(9000 + fid)
            if name in ("ln", "log10"):
                ys = binade_uniform(rng, 10_000_000, -1074, 1023, signs=(1,))
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
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
