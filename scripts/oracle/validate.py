"""Independent check of tests/golden/*.tsv (format 1) and src/constants.cyr, using no Rust
and no Python packages. Usage, from the repo root:

    python3 -I scripts/oracle/validate.py .

1. Format: every file has the five-line header (with the right row count), every row has five
   tab fields, a known class and well-formed tokens; modules.txt lists exactly the .tsv files.
2. constants.tsv matches the `pub const NAME: f64 = <literal>;` items of rust-old/src/*.rs,
   all of them (Python's float() is correctly rounded, like rustc's literal conversion).
3. src/constants.cyr: every `const KM_X = <decimal>;` is a positional decimal (no exponent,
   no hex), names its Rust path, and parses to that constant's golden bits.
4. math.tsv, the exact functions, recomputed in Python with Rust's semantics: sqrt
   (math.sqrt is correctly rounded), abs, floor, trunc, round (half away from zero), min and
   max (a NaN operand is ignored; the sign of a zero result is unspecified in Rust, so only
   the value is checked and the observed rule printed), signum, clamp, total_cmp, is_finite,
   and powi (compiler-builtins' __powidf2 loop). cbrt is proven correctly rounded with exact
   rational arithmetic, which is why it is class A.
5. math.tsv, the class-C functions (glibc's bits): the IEEE / C99 Annex F special cases, and
   every finite row against a correctly rounded reference computed with the decimal module
   at 60+ digits. Each function must be within glibc's documented bound (1 ulp; log10 2) of
   it; the histogram is printed. This proves the reference is the function it names, on
   every row, independently of the platform that recorded it.
Every function in math.tsv must have a rule here: a new one fails until it gets one.
"""
import math, re, struct, sys
from decimal import Decimal, localcontext, Overflow, Underflow, InvalidOperation
from fractions import Fraction
from pathlib import Path

root = Path(sys.argv[1])
gold = root / "tests/golden"

def b2f(b): return struct.unpack("<d", struct.pack("<Q", b))[0]
def f2b(x): return struct.unpack("<Q", struct.pack("<d", x))[0]
def isnan_b(b): return (b & 0x7FF0000000000000) == 0x7FF0000000000000 and (b & 0xFFFFFFFFFFFFF) != 0
def ordered(b): return -(b & 0x7FFFFFFFFFFFFFFF) if b >> 63 else b
def ulp_dist(g, w):
    if isnan_b(g) or isnan_b(w):
        return 0 if isnan_b(g) and isnan_b(w) else math.inf
    return abs(ordered(g) - ordered(w))

fails = 0
def check(ok, msg):
    global fails
    if not ok:
        fails += 1
        if fails <= 30:
            print("FAIL", msg)

# ── 1. format ───────────────────────────────────────────────────────────────────────
TOK = {
    "f": re.compile(r"[0-9a-f]{16}"),
    "i": re.compile(r"-?[0-9]+"),
    "e": re.compile(r"[A-Za-z][A-Za-z0-9_]*\.[A-Za-z][A-Za-z0-9_]*"),
    "s": re.compile(r"(?:[\x21-\x24\x26-\x7e]|%[0-9A-F]{2})*"),
    "n": re.compile(r""),
    "E": re.compile(r"[A-Za-z][A-Za-z0-9_]*"),
    "v": re.compile(r"[0-9]+"),
}

def parse_tok(t, where):
    kind, sep, val = t.partition(":")
    check(sep == ":" and kind in TOK and TOK[kind].fullmatch(val), f"{where}: bad token {t!r}")
    if kind == "f":
        b = int(val, 16)
        return ("f", b2f(b), b)
    if kind in ("i", "v"):
        return (kind, int(val), None)
    if kind == "s":
        raw = re.sub(rb"%([0-9A-F]{2})", lambda m: bytes([int(m.group(1), 16)]), val.encode())
        try:
            return ("s", raw.decode("utf-8"), None)
        except UnicodeDecodeError:
            check(False, f"{where}: s: token is not UTF-8")
            return ("s", "", None)
    return (kind, val, None)

def load(name):
    lines = (gold / f"{name}.tsv").read_text().split("\n")
    check(lines[-1] == "", f"{name}.tsv: no final newline")
    lines = lines[:-1]
    head = lines[:5]
    check(len(head) == 5 and all(l.startswith("# ") for l in head), f"{name}.tsv: header")
    check(head[0] == f"# kimiya golden corpus, format 1: {name}", f"{name}.tsv: header line 1")
    check(head[1].startswith("# oracle: kimiya ") and "; libm: glibc " in head[1], f"{name}.tsv: provenance")
    check(head[2].startswith("# regenerate: cargo run --release --locked"), f"{name}.tsv: header line 3")
    check(head[3].startswith("# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note"), f"{name}.tsv: header line 4")
    out = []
    for n, line in enumerate(lines[5:], 6):
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        where = f"{name}.tsv:{n}"
        check(len(cols) == 5, f"{where}: {len(cols)} fields")
        if len(cols) != 5:
            continue
        fn, cls, ins, outs, note = cols
        check(cls in ("A", "B", "C", "D"), f"{where}: class {cls!r}")
        check(fn != "" and note.startswith(fn + "("), f"{where}: fn / note")
        pi = [] if ins == "-" else [parse_tok(t, where) for t in ins.split(" ")]
        po = [] if outs == "-" else [parse_tok(t, where) for t in outs.split(" ")]
        out.append((fn, cls, pi, po, where))
    check(head[4] == f"# rows: {len(out)}", f"{name}.tsv: header says {head[4]!r}, has {len(out)} rows")
    return out

tsvs = sorted(p.stem for p in gold.glob("*.tsv"))
index = [l for l in (gold / "modules.txt").read_text().splitlines() if l and not l.startswith("#")]
check(index == tsvs, f"modules.txt lists {index}, tests/golden has {tsvs}")
R = {name: load(name) for name in tsvs}
print("parsed:", {k: len(v) for k, v in R.items()})

def same(x, want_bits):
    if isinstance(x, float) and math.isnan(x):
        return isnan_b(want_bits)
    return f2b(x) == want_bits

# ── 2. constants vs the Rust source literals ───────────────────────────────────────
lits = {}
for p in sorted((root / "rust-old/src").glob("*.rs")):
    for m in re.finditer(r"^pub const ([A-Z0-9_]+): f64 = ([^;]+);", p.read_text(), re.M):
        lit = m.group(2).replace("_", "")
        check(re.fullmatch(r"[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?", lit), f"{p.name}: {m.group(1)} = {m.group(2)} is not a literal")
        lits[f"{p.stem}::{m.group(1)}"] = float(lit)
golden_consts = {}
for fn, cls, ins, outs, where in R.get("constants", []):
    check(fn == "const" and cls == "A" and len(ins) == 1 and ins[0][0] == "s" and len(outs) == 1 and outs[0][0] == "f", f"{where}: shape")
    path = ins[0][1]
    check(path not in golden_consts, f"{where}: {path} twice")
    golden_consts[path] = outs[0][2]
    check(path in lits, f"{where}: {path} is not a pub const f64 in rust-old/src")
    if path in lits:
        check(same(lits[path], outs[0][2]), f"{where}: {path} = {lits[path]!r}, golden {outs[0][1]!r}")
missing = sorted(set(lits) - set(golden_consts))
check(not missing, f"constants.tsv lacks {missing}")
print(f"constants checked: {len(golden_consts)} rows against {len(lits)} Rust literals")

# ── 3. src/constants.cyr ────────────────────────────────────────────────────────────
cyr = (root / "src/constants.cyr").read_text()
n = 0
for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = (.+);$", cyr, re.M):
    path, name, lit = m.groups()
    check(name == "KM_" + path.split("::")[1], f"constants.cyr: {name} names kimiya::{path}")
    check(re.fullmatch(r"-?[0-9]+\.[0-9]+", lit), f"constants.cyr: {name} = {lit} is not a positional decimal")
    check(path in golden_consts and f2b(float(lit)) == golden_consts.get(path), f"constants.cyr: {name} = {lit} is not the golden {path}")
    n += 1
check(n == len(re.findall(r"^const ", cyr, re.M)) and n > 0, "constants.cyr: a const without its `# Rust:` line")
print("constants.cyr literals checked:", n)

# ── 4. math.tsv: exact functions ───────────────────────────────────────────────────
def r_sqrt(x):
    if math.isnan(x) or x < 0: return math.nan
    return math.sqrt(x) if x != 0 else x
def r_floor(x):
    if not math.isfinite(x) or abs(x) >= 2.0**52: return x
    return math.copysign(float(math.floor(x)), x)
def r_trunc(x):
    if not math.isfinite(x) or abs(x) >= 2.0**52: return x
    return math.copysign(float(math.trunc(x)), x)
def r_round(x):  # half away from zero, exactly
    if not math.isfinite(x) or abs(x) >= 2.0**52: return x
    t = math.trunc(x)
    if abs(Fraction(x) - t) >= Fraction(1, 2):
        t += 1 if x > 0 else -1
    return math.copysign(float(t), x)
def r_signum(x):
    if math.isnan(x): return math.nan
    return math.copysign(1.0, x)
def r_clamp(x, lo, hi):
    if x < lo: x = lo
    if x > hi: x = hi
    return x
def r_total_cmp(x, y):
    def key(b):
        b = b - (1 << 64) if b >> 63 else b
        return b ^ ((((b >> 63) & 0xFFFFFFFFFFFFFFFF) >> 1))
    a, c = key(f2b(x)), key(f2b(y))
    return (a > c) - (a < c)
def r_powi(a, b):  # compiler-builtins float/pow.rs
    recip = b < 0
    p = abs(b)
    mul = 1.0
    while True:
        if p & 1: mul *= a
        p >>= 1
        if p == 0: break
        a *= a
    if not recip: return mul
    if mul == 0.0: return math.copysign(math.inf, mul)
    return 1.0 / mul

def nextup(y): return math.nextafter(y, math.inf)
def nextdown(y): return math.nextafter(y, -math.inf)
def cbrt_correctly_rounded(x, y):
    """y == cbrt(x) correctly rounded (to nearest even), by exact rational arithmetic."""
    if math.isnan(x): return math.isnan(y)
    if x == 0 or math.isinf(x): return f2b(y) == f2b(x)
    if math.copysign(1.0, y) != math.copysign(1.0, x) or not math.isfinite(y): return False
    ax, ay = Fraction(abs(x)), abs(y)
    lo = (Fraction(nextdown(ay)) + Fraction(ay)) / 2
    hi = (Fraction(ay) + Fraction(nextup(ay))) / 2
    lo3, hi3 = lo ** 3, hi ** 3
    if lo3 < ax < hi3: return True
    even = f2b(ay) & 1 == 0
    return even and (ax == lo3 or ax == hi3)

minmax_zero_rule = set()
def check_minmax(fn, x, y, got_bits, where):
    if math.isnan(x) and math.isnan(y): ok = isnan_b(got_bits)
    elif math.isnan(x): ok = f2b(y) == got_bits
    elif math.isnan(y): ok = f2b(x) == got_bits
    else:
        want = max(x, y) if fn == "f64::max" else min(x, y)
        g = b2f(got_bits)
        ok = g == want and (want != 0 or got_bits in (f2b(x), f2b(y)))
        if x == 0 and y == 0 and f2b(x) != f2b(y):
            minmax_zero_rule.add(f"{fn}({x}, {y}) = {g} ({'second' if got_bits == f2b(y) else 'first'} operand)")
    check(ok, f"{where}: {fn}({x!r}, {y!r}) -> {b2f(got_bits)!r}")

exact = {
    "f64::sqrt": lambda a: r_sqrt(a[0]),
    "f64::abs": lambda a: abs(a[0]) if not math.isnan(a[0]) else math.nan,
    "f64::floor": lambda a: r_floor(a[0]),
    "f64::trunc": lambda a: r_trunc(a[0]),
    "f64::round": lambda a: r_round(a[0]),
    "f64::signum": lambda a: r_signum(a[0]),
    "f64::clamp": lambda a: r_clamp(*a),
    "f64::powi": lambda a: r_powi(*a),
}

# ── 5. math.tsv: class C against a correctly rounded decimal reference ──────────────
def pi_dec(prec):
    with localcontext() as ctx:  # the decimal module documentation's recipe
        ctx.prec = prec + 2
        three = Decimal(3)
        lasts, t, s, k, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            k, na = k + na, na + 8
            d, da = d + da, da + 32
            t = (t * k) / d
            s += t
    with localcontext() as ctx:
        ctx.prec = prec
        return +s

PREC = 60
PI = pi_dec(450)
def prec_for(x):  # digits needed to reduce x by 2π and keep PREC significant digits
    return PREC + max(0, int(abs(x)).bit_length() * 31 // 100 if math.isfinite(x) else 0) + 10

def d_sin_cos(x):
    """(sin x, cos x) as Decimals with PREC good digits (Taylor after reduction mod 2π)."""
    with localcontext() as ctx:
        ctx.prec = prec_for(x)
        two_pi = 2 * PI
        r = Decimal(x)
        r = r - two_pi * (r / two_pi).to_integral_value()
        ctx.prec = PREC + 10
        r = +r
        s, c, term_s, term_c, k = r, Decimal(1), r, Decimal(1), 0
        r2 = r * r
        while True:
            k += 2
            term_s = -term_s * r2 / (k * (k + 1))
            term_c = -term_c * r2 / ((k - 1) * k)
            if abs(term_s) < Decimal(10) ** (-(PREC + 20)) and abs(term_c) < Decimal(10) ** (-(PREC + 20)):
                break
            s += term_s
            c += term_c
        return s, c

def d_atan(z):
    """atan of a Decimal z, PREC good digits."""
    with localcontext() as ctx:
        ctx.prec = PREC + 20
        if z == 0: return Decimal(0)
        if abs(z) > 1:
            half = PI / 2 if z > 0 else -PI / 2
            return half - d_atan(1 / z)
        halvings = 0
        while abs(z) > Decimal("0.05"):
            z = z / (1 + (1 + z * z).sqrt())
            halvings += 1
        s, term, k, z2 = z, z, 1, z * z
        while True:
            term = -term * z2
            k += 2
            t = term / k
            if abs(t) < Decimal(10) ** (-(PREC + 25)): break
            s += t
        return s * (2 ** halvings)

def d_ln(x): return Decimal(x).ln()
def d_exp_or_inf(z):
    """exp of a Decimal, saturating far outside the double range."""
    if z > 710: return math.inf
    if z < -746: return 0.0
    return z.exp()

def to_double(d):
    if isinstance(d, float): return d
    return float(d)

def ref_exp(x):
    if math.isnan(x): return math.nan
    if math.isinf(x): return math.inf if x > 0 else 0.0
    with localcontext() as ctx:
        ctx.prec = PREC
        return to_double(d_exp_or_inf(Decimal(x)))

def ref_expm1(x):
    if math.isnan(x): return math.nan
    if math.isinf(x): return math.inf if x > 0 else -1.0
    if x == 0: return x
    with localcontext() as ctx:
        ctx.prec = PREC + 20
        d = Decimal(x)
        if abs(x) < 1e-3:  # the series keeps relative precision for tiny x
            s, term, k = d, d, 1
            while abs(term) > abs(s) * Decimal(10) ** (-(PREC + 10)):
                k += 1
                term = term * d / k
                s += term
            return to_double(s)
        e = d_exp_or_inf(d)
        return math.inf if e == math.inf else to_double(e - 1)

def ref_log(x, base10):
    if math.isnan(x) or x < 0: return math.nan
    if x == 0: return -math.inf
    if math.isinf(x): return math.inf
    with localcontext() as ctx:
        ctx.prec = PREC
        return to_double(Decimal(x).log10() if base10 else Decimal(x).ln())

def is_int(y): return math.isfinite(y) and y == math.floor(y)
def is_odd_int(y): return is_int(y) and abs(y) < 2.0**53 and int(y) % 2 == 1

def ref_pow(x, y):
    """C99 Annex F.9.4.4, then a correctly rounded exp(y·ln|x|)."""
    if y == 0: return 1.0
    if x == 1: return 1.0
    if math.isnan(x) or math.isnan(y): return math.nan
    if x == 0:
        neg = math.copysign(1.0, x) < 0 and is_odd_int(y)
        if y < 0: return -math.inf if neg else math.inf
        return -0.0 if neg else 0.0
    if math.isinf(y):
        if x == -1: return 1.0
        if (abs(x) < 1) == (y < 0): return math.inf
        return 0.0
    if math.isinf(x):
        if x > 0: return 0.0 if y < 0 else math.inf
        odd = is_odd_int(y)
        if y < 0: return -0.0 if odd else 0.0
        return -math.inf if odd else math.inf
    if x < 0 and not is_int(y): return math.nan
    sign = -1.0 if x < 0 and is_odd_int(y) else 1.0
    with localcontext() as ctx:
        ctx.prec = PREC + 30
        z = Decimal(y) * Decimal(abs(x)).ln()
        e = d_exp_or_inf(z)
        ctx.prec = PREC
        return sign * to_double(+e if isinstance(e, Decimal) else e)

def ref_sin(x):
    if not math.isfinite(x): return math.nan
    if x == 0: return x
    with localcontext() as ctx:
        ctx.prec = PREC
        return to_double(+d_sin_cos(x)[0])

def ref_cos(x):
    if not math.isfinite(x): return math.nan
    with localcontext() as ctx:
        ctx.prec = PREC
        return to_double(+d_sin_cos(x)[1])

def ref_acos(x):
    if math.isnan(x) or abs(x) > 1: return math.nan
    if x == 1: return 0.0
    with localcontext() as ctx:
        ctx.prec = PREC + 20
        if x == -1: r = PI
        else:
            d = Decimal(x)
            r = PI / 2 - d_atan(d / (1 - d * d).sqrt())
        ctx.prec = PREC
        return to_double(+r)

def ref_atan2(y, x):
    """C99 Annex F.9.1.4, then the quadrant-corrected atan(y/x)."""
    if math.isnan(x) or math.isnan(y): return math.nan
    sy = -1 if math.copysign(1.0, y) < 0 else 1
    with localcontext() as ctx:
        ctx.prec = PREC + 20
        if y == 0:
            if math.copysign(1.0, x) < 0: r = sy * PI
            else: return math.copysign(0.0, y)
        elif math.isinf(y):
            if math.isinf(x): r = sy * (3 * PI / 4 if x < 0 else PI / 4)
            else: r = sy * PI / 2
        elif x == 0: r = sy * PI / 2
        elif math.isinf(x):
            if x > 0: return math.copysign(0.0, y)
            r = sy * PI
        else:
            r = d_atan(Decimal(y) / Decimal(x))
            if x < 0: r = r + sy * PI
        ctx.prec = PREC
        return to_double(+r)

transcendental = {
    "f64::exp": (lambda a: ref_exp(a[0]), 1),
    "f64::exp_m1": (lambda a: ref_expm1(a[0]), 1),
    "f64::ln": (lambda a: ref_log(a[0], False), 1),
    "f64::log10": (lambda a: ref_log(a[0], True), 2),
    "f64::powf": (lambda a: ref_pow(a[0], a[1]), 1),
    "f64::powf_base10": (lambda a: ref_pow(10.0, a[0]), 1),
    "f64::powf_base2": (lambda a: ref_pow(2.0, a[0]), 1),
    "f64::sin": (lambda a: ref_sin(a[0]), 1),
    "f64::cos": (lambda a: ref_cos(a[0]), 1),
    "f64::acos": (lambda a: ref_acos(a[0]), 1),
    "f64::atan2": (lambda a: ref_atan2(a[0], a[1]), 1),
}

hist = {}
counts = {}
unknown = set()
for fn, cls, ins, outs, where in R.get("math", []):
    counts[fn] = counts.get(fn, 0) + 1
    a = [t[1] for t in ins]
    if fn in exact:
        check(cls == ("B" if fn == "f64::powi" else "A"), f"{where}: {fn} class {cls}")
        check(same(exact[fn](a), outs[0][2]), f"{where}: {fn}{tuple(a)} -> {outs[0][1]!r}, want {exact[fn](a)!r}")
    elif fn in ("f64::max", "f64::min"):
        check(cls == "A", f"{where}: {fn} class {cls}")
        check_minmax(fn, a[0], a[1], outs[0][2], where)
    elif fn == "f64::total_cmp":
        check(cls == "A" and outs[0][0] == "i", f"{where}: shape")
        check(r_total_cmp(a[0], a[1]) == outs[0][1], f"{where}: total_cmp{tuple(a)} -> {outs[0][1]}")
    elif fn == "f64::is_finite":
        check(cls == "A" and outs[0][0] == "i", f"{where}: shape")
        check(int(math.isfinite(a[0])) == outs[0][1], f"{where}: is_finite({a[0]!r}) -> {outs[0][1]}")
    elif fn == "f64::cbrt":
        check(cls == "A", f"{where}: {fn} class {cls}")
        check(cbrt_correctly_rounded(a[0], outs[0][1]), f"{where}: cbrt({a[0]!r}) = {outs[0][1]!r} is not correctly rounded")
    elif fn in transcendental:
        check(cls == "C", f"{where}: {fn} class {cls}")
        ref, bound = transcendental[fn]
        try:
            want = ref(a)
        except (Overflow, Underflow, InvalidOperation, ZeroDivisionError) as e:
            check(False, f"{where}: reference failed: {e!r}")
            continue
        d = ulp_dist(outs[0][2], f2b(want))
        h = hist.setdefault(fn, {0: 0, 1: 0, 2: 0, "more": 0, "max": 0})
        h[d if d in (0, 1, 2) else "more"] += 1
        if d != math.inf: h["max"] = max(h["max"], d)
        check(d <= bound, f"{where}: {fn}{tuple(a)} = {outs[0][1]!r}, correctly rounded {want!r}: {d} ulp > {bound}")
    else:
        unknown.add(fn)
check(not unknown, f"math.tsv functions with no rule here: {sorted(unknown)}")
print("math rows checked:", sum(counts.values()), "over", len(counts), "functions")
for rule in sorted(minmax_zero_rule):
    print("  observed:", rule)
print("glibc vs correctly rounded (decimal reference), ulp histogram:")
for fn in sorted(hist):
    h = hist[fn]
    print(f"  {fn:20} 0:{h[0]} 1:{h[1]} 2:{h[2]} more:{h['more']} max:{h['max']} (bound {transcendental[fn][1]})")
print("FAILURES:", fails)
sys.exit(1 if fails else 0)
