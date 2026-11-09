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
5. math.tsv, the class-C functions (glibc's bits): every row against its correctly rounded value from
   MPFR (scripts/oracle/mpfr_ref.py: 53 bits, binary64's exponent range, subnormalized, so exact ties,
   subnormal results and the C99 Annex F special cases are exact). Each function must be within glibc's
   documented bound (1 ulp; log10 2) of it; the histogram is printed. This proves the reference is the
   function it names, on every row, independently of the platform that recorded it. A row where glibc is
   not correctly rounded gets an ADR 0006 correction (item 12).
6. element.tsv against the ELEMENTS / ELEMENT_PROPERTIES literals parsed from
   rust-old/src/element.rs: every lookup_by_number and lookup_properties row (all 256 u8), every
   lookup_by_symbol row (exact, case-sensitive match), and the generated KmElementCategory order.
7. molecule.tsv recomputed from the Rust semantics with those masses: molecular weight summed in
   atom order (bit-exact; Python floats are IEEE doubles), formula, saturating total_atoms, the
   InvalidElement errors, and the presets parsed from rust-old/src/molecule.rs.
8. reaction.tsv, thermo.tsv, phase.tsv and gas.tsv (M2a), every row recomputed from Rust 1.2.0's
   code transcribed here in its operation order (Python floats are IEEE doubles), with the
   validation guards in Rust's order. The tables and presets (ANTOINE_DATA, PR_*, PA_PER_MMHG) are
   parsed from rust-old/src literals, and the generated src/phase_data.cyr and src/gas_data.cyr
   literals are checked against them. Class B rows must be bit-identical. Class C rows call
   Python's math (glibc's exp, log, pow, acos, cos; cbrt made correctly rounded, as Rust's
   compiler-builtins cbrt is): bit-identical when this machine is the corpus's recorded libm
   platform, else only the Ok/Err shape is checked (printed).
9. thermochem.tsv (M2b), every row recomputed the same way, from THERMOCHEM_DATA and SHOMATE_DATA parsed from
   rust-old/src/thermochem.rs: the reaction sums (with ADR 0002 #1's J/mol × 1000.0), van 't Hoff, the Shomate
   polynomial, Rust's private Segments, and the flame temperature with hisab 1.4.0's newton_raphson and bisection
   transcribed. Each flame row's noted solver path must match the recomputation; each of Rust's typed flame
   constants must change some row when perturbed (the corpus pins it). src/thermochem_data.cyr is checked
   against the literals.
10. solution.tsv and electrochemistry.tsv (M3a), every row recomputed the same way: the weak acid / base
   charge-balance Newton transcribed (each row's noted exit path checked; CHARGE_BALANCE_REL_TOL must change some
   row when halved and when doubled; the 50-iteration cap is reported unobservable), ionic strength with Rust's
   .sum() from -0.0, STANDARD_POTENTIALS parsed from rust-old/src/electrochemistry.rs, and
   src/electrochemistry_data.cyr checked against the literal.
11. kinetics.tsv (M3b), every row recomputed the same way (glibc's exp, expm1 and pow via Python's math;
   lineweaver_burk's tuple as two floats); Rust's LN_2 and private FIRST_ORDER_TOLERANCE (typed inline in
   src/kinetics.cyr) must each change some golden row when perturbed.
12. The corrected oracle (ADR 0006). Every golden row whose replica calls a transcendental primitive (the
   M2–M3b modules' rows; the math rows of the primitives themselves) is evaluated a second time with every
   primitive result replaced by its correctly rounded value from MPFR, through the same transcribed code (so
   a corrected result may take another branch, or turn Ok into Err). Where the corrected outputs differ from
   the oracle's, the row gets a correction in tests/golden/<module>.cr.tsv, in golden-file order, with the
   oracle's outputs and each call glibc misrounded in its note. This run fails when a .cr.tsv is missing,
   stale or has an extra row; it recomputes every correction (the outputs through the replica, each listed
   call against MPFR). Per-primitive counts of glibc's misroundings over every golden call are printed, and
   glibc's cbrt (the oracle's is Rust's correctly rounded one) is measured on the golden cbrt arguments.
Every function in every .tsv must have a rule here: a new one fails until it gets one.

    python3 -I scripts/oracle/validate.py . --write-corrections

writes the tests/golden/<module>.cr.tsv files (deterministic: one per module with at least one correction,
rows in golden-file order; a module without corrections has none, and a stale file is removed), then runs
every check. It refuses to write unless this machine is the corpus's recorded libm platform, since the notes
hold glibc's results.

    python3 -I scripts/oracle/validate.py . --sensitivity

also measures the conditioning of each class-C function (M2a–M3b) on its golden grid: every golden row is
re-evaluated with each primitive result nudged one ulp up and down, one call at a time, and S is the largest
change of any output, in ulp. A result that is 0, ±∞ or NaN, or acos(±1), is not nudged; a nudge that flips
Ok/Err or the root count is reported as a branch flip. It prints S, c (the most nudged calls in one row), the
retired rule's K = 2·S·c + 2 for comparison, and the calls where glibc's result is not correctly rounded (MPFR
decides it): exactly the calls the corrected oracle changes. For each function whose K passes 64 it lists the
rows that put it there (how many nudge a subnormal exp / pow result, the worst three) and S over the other rows.
Since ADR 0006 (M3c-3) it is a diagnostic only: the port evaluates the oracle's operations with correctly rounded
primitives, as the corrected oracle does, so every class-C row is gated at 0 ulp in tests/parity.tcyr whatever
its S (port-parity.md §10; §5–§8 hold the bounds the rule set before).

MPFR (libmpfr.so.6, 4.x) is the one system library this needs, through ctypes (scripts/oracle/mpfr_ref.py).
"""
import importlib.util, math, platform, re, struct, sys
from fractions import Fraction
from pathlib import Path

root = Path(sys.argv[1])
SENSITIVITY = "--sensitivity" in sys.argv[2:]
WRITE_CORRECTIONS = "--write-corrections" in sys.argv[2:]
gold = root / "tests/golden"

# mpfr_ref.py sits beside this script; `-I` keeps the script's directory off sys.path, so load it by path
# (and leave no __pycache__ behind).
sys.dont_write_bytecode = True
_mpfr_spec = importlib.util.spec_from_file_location("mpfr_ref", Path(__file__).resolve().parent / "mpfr_ref.py")
MPFR = importlib.util.module_from_spec(_mpfr_spec)
_mpfr_spec.loader.exec_module(MPFR)

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

NOTES = {}       # "<module>.tsv:<line>" -> the row's note (the flame rows' solver path, M2b)
RAW = {}         # "<module>.tsv:<line>" -> the row's fn, inputs and outputs fields, verbatim (the .cr.tsv key, ADR 0006)
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
        NOTES[where] = note
        RAW[where] = (fn, ins, outs)
    check(head[4] == f"# rows: {len(out)}", f"{name}.tsv: header says {head[4]!r}, has {len(out)} rows")
    return out

CR_SUFFIX = ".cr.tsv"   # tests/golden/<module>.cr.tsv: the corrected oracle (ADR 0006, item 12), not a module
tsvs = sorted(p.name[:-4] for p in gold.glob("*.tsv") if not p.name.endswith(CR_SUFFIX))
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

# ── 5. math.tsv: class C against MPFR (ADR 0006) ───────────────────────────────────────
# The correctly rounded reference is MPFR (scripts/oracle/mpfr_ref.py): exact, so it decides exact ties and
# subnormal results as IEEE round-to-nearest-even does, where a fixed-precision decimal reference misjudged ties
# (port-parity.md §8). The C99 Annex F special cases are MPFR's too.
def same_f(x, y):
    """x and y are the same binary64 (any NaN is any NaN)."""
    if math.isnan(x) or math.isnan(y): return math.isnan(x) and math.isnan(y)
    return f2b(x) == f2b(y)

_cr_vals = {}
def cr_of(name, args):
    """The correctly rounded binary64 value of MPFR's `name` at `args` (floats), cached by their bits."""
    key = (name,) + tuple(f2b(x) for x in args)
    v = _cr_vals.get(key)
    if v is None:
        v = MPFR.cr(name, args[0]) if len(args) == 1 else MPFR.cr2(name, args[0], args[1])
        if name == "pow" and f2b(args[0]) == f2b(10.0):   # 10.0.powf(y) is pow(10, y), correctly rounded exp10(y)
            check(same_f(v, MPFR.exp10(args[1])), f"MPFR: pow(10, {args[1]!r}) is not exp10")
        _cr_vals[key] = v
    return v

# The ADR 0006 corrections (item 12): module -> [(where, corrected output tokens, misrounded calls)], in golden
# order. A misrounded call is (MPFR name, arguments, glibc's result, the correctly rounded result).
CORRECTIONS = {}
def add_correction(mod, where, ctoks, mis):
    CORRECTIONS.setdefault(mod, []).append((where, ctoks, mis))

# The math functions that are transcendental primitives: the MPFR function the oracle's call rounds, its
# arguments from the row's inputs, and glibc's documented bound in ulp. powf_base10 is pow(10, y);
# powf_base2 is glibc's exp2 (LLVM's rewrite of 2.0_f64.powf(y), README); atan2's inputs are (y, x).
CR_MATH = {
    "f64::exp": ("exp", lambda a: (a[0],), 1),
    "f64::exp_m1": ("expm1", lambda a: (a[0],), 1),
    "f64::ln": ("log", lambda a: (a[0],), 1),
    "f64::log10": ("log10", lambda a: (a[0],), 2),
    "f64::powf": ("pow", lambda a: (a[0], a[1]), 1),
    "f64::powf_base10": ("pow", lambda a: (10.0, a[0]), 1),
    "f64::powf_base2": ("exp2", lambda a: (a[0],), 1),
    "f64::sin": ("sin", lambda a: (a[0],), 1),
    "f64::cos": ("cos", lambda a: (a[0],), 1),
    "f64::acos": ("acos", lambda a: (a[0],), 1),
    "f64::atan2": ("atan2", lambda a: (a[0], a[1]), 1),
}
CBRT_ARGS = []   # every golden cbrt argument (math rows, the Peng–Robinson cubic): glibc's cbrt is measured on them

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
        check(same_f(cr_of("cbrt", (a[0],)), outs[0][1]), f"{where}: cbrt({a[0]!r}) = {outs[0][1]!r} is not MPFR's")
        CBRT_ARGS.append(a[0])
    elif fn in CR_MATH:
        check(cls == "C", f"{where}: {fn} class {cls}")
        name, argf, bound = CR_MATH[fn]
        args = argf(a)
        want = cr_of(name, args)
        d = ulp_dist(outs[0][2], f2b(want))
        h = hist.setdefault(fn, {0: 0, 1: 0, 2: 0, "more": 0, "max": 0})
        h[d if d in (0, 1, 2) else "more"] += 1
        if d != math.inf: h["max"] = max(h["max"], d)
        check(d <= bound, f"{where}: {fn}{tuple(a)} = {outs[0][1]!r}, correctly rounded {want!r}: {d} ulp > {bound}")
        if not same_f(outs[0][1], want):
            add_correction("math", where, [("f", want, f2b(want))], [(name, args, outs[0][1], want)])
    else:
        unknown.add(fn)
check(not unknown, f"math.tsv functions with no rule here: {sorted(unknown)}")
print("math rows checked:", sum(counts.values()), "over", len(counts), "functions")
for rule in sorted(minmax_zero_rule):
    print("  observed:", rule)
print(f"glibc vs correctly rounded (MPFR {MPFR.VERSION}), ulp histogram:")
for fn in sorted(hist):
    h = hist[fn]
    print(f"  {fn:20} 0:{h[0]} 1:{h[1]} 2:{h[2]} more:{h['more']} max:{h['max']} (bound {CR_MATH[fn][2]})")
# ── 6. element.tsv against the Rust source literals ──────────────────────────────────
esrc = (root / "rust-old/src/element.rs").read_text()
def literal_blocks(name, struct):
    body = esrc[esrc.index(f"pub static {name}:"):]
    body = body[: body.index("\n];")]
    return re.findall(rf"{struct} \{{(.*?)\n    \}}", body, re.S)
ELEMS = {}
for b in literal_blocks("ELEMENTS", "Element"):
    f_ = dict(re.findall(r"^\s*(\w+): (.+?),\s*$", b, re.M))
    z = int(f_["atomic_number"])
    en = f_["electronegativity"]
    ELEMS[z] = (z, f_["symbol"].strip('"'), f_["name"].strip('"'), float(f_["atomic_mass"]),
                None if en == "None" else float(re.fullmatch(r"Some\((.+)\)", en).group(1)),
                f_["category"].split("::")[1])
PROPS = {}
for b in literal_blocks("ELEMENT_PROPERTIES", "ElementProperties"):
    f_ = dict(re.findall(r"^\s*(\w+): (.+?),\s*$", b, re.M))
    states = [int(x) for x in re.fullmatch(r"&\[(.*)\]", f_["oxidation_states"]).group(1).split(",")]
    PROPS[int(f_["atomic_number"])] = (int(f_["atomic_number"]), float(f_["ionization_energy_1"]),
                                       float(f_["electron_affinity"]), float(f_["covalent_radius_pm"]), states)
check(sorted(ELEMS) == list(range(1, 119)), f"element.rs: ELEMENTS parsed as Z {sorted(ELEMS)[:3]}.. ({len(ELEMS)})")
check(len(PROPS) == 32, f"element.rs: ELEMENT_PROPERTIES parsed as {len(PROPS)} rows")
cats = re.findall(r"^    ([A-Z]\w+),$", esrc[esrc.index("pub enum ElementCategory {"):esrc.index("pub struct Element {")], re.M)
def element_out(e):
    if e is None: return [("n", "", None)]
    z, sym, name, mass, en, cat = e
    return [("i", z, None), ("s", sym, None), ("s", name, None), ("f", mass, f2b(mass)),
            ("n", "", None) if en is None else ("f", en, f2b(en)), ("e", f"ElementCategory.{cat}", None)]
def props_out(pr):
    if pr is None: return [("n", "", None)]
    z, ie, ea, rad, st = pr
    return ([("i", z, None), ("f", ie, f2b(ie)), ("f", ea, f2b(ea)), ("f", rad, f2b(rad)), ("v", len(st), None)]
            + [("i", x, None) for x in st])
def same_toks(got, want):
    return [(k, v if k != "f" else None, b) for k, v, b in got] == [(k, v if k != "f" else None, b) for k, v, b in want]
ecount = {}
for fn, cls, ins, outs, where in R.get("element", []):
    ecount[fn] = ecount.get(fn, 0) + 1
    check(cls == "A", f"{where}: {fn} class {cls}")
    if fn == "element::lookup_by_number":
        z = ins[0][1]
        check(same_toks(outs, element_out(ELEMS.get(z))), f"{where}: lookup_by_number({z})")
    elif fn == "element::lookup_by_symbol":
        sym = ins[0][1]
        hit = next((e for e in ELEMS.values() if e[1] == sym), None)  # exact, case-sensitive
        check(same_toks(outs, element_out(hit)), f"{where}: lookup_by_symbol({sym!r})")
    elif fn == "element::lookup_properties":
        z = ins[0][1]
        check(same_toks(outs, props_out(PROPS.get(z))), f"{where}: lookup_properties({z})")
    else:
        check(False, f"{where}: element.tsv function {fn} has no rule here")
check(ecount.get("element::lookup_by_number") == 256 and ecount.get("element::lookup_properties") == 256,
      "element.tsv: lookup_by_number / lookup_properties must cover every u8")
syms = {i[0][1] for fn, _, i, _, _ in R.get("element", []) if fn == "element::lookup_by_symbol"}
check(all(e[1] in syms for e in ELEMS.values()), "element.tsv: a symbol is not probed")
print(f"element rows checked: {sum(ecount.values())} against {len(ELEMS)} ELEMENTS and {len(PROPS)} "
      f"ELEMENT_PROPERTIES literals; categories in source order: {len(cats)}")
# The generated enum follows the Rust declaration order.
edata = (root / "src/element_data.cyr").read_text()
gen = re.findall(r"^    (KM_ELEMENT_CATEGORY_[A-Z_]+) = (\d+);", edata, re.M)
check([g[0] for g in gen] == ["KM_ELEMENT_CATEGORY_" + re.sub(r"(?<=[a-z])([A-Z])", r"_\1", c).upper() for c in cats]
      and [int(g[1]) for g in gen] == list(range(len(cats))), "element_data.cyr: KmElementCategory is not ElementCategory in order")

# ── 7. molecule.tsv, recomputed from the Rust source's semantics ────────────────────
msrc = (root / "rust-old/src/molecule.rs").read_text()
presets = {}
for name, pairs in re.findall(r"pub fn (\w+)\(\) -> Self \{\s*Self::new\(&\[(.*?)\]\)", msrc, re.S):
    presets[f"Molecule::{name}"] = [(int(a), int(b)) for a, b in re.findall(r"\((\d+), (\d+)\)", pairs)]
check(sorted(presets) == ["Molecule::carbon_dioxide", "Molecule::glucose", "Molecule::methane", "Molecule::water"],
      f"molecule.rs: presets parsed as {sorted(presets)}")
def mol_in(toks):
    n = toks[0][1]
    check(toks[0][0] == "v" and len(toks) == 1 + 2 * n, "molecule tokens")
    return [(toks[1 + 2 * k][1], toks[2 + 2 * k][1]) for k in range(n)]
def mol_toks(pairs):
    return [("v", len(pairs), None)] + [t for z, c in pairs for t in (("i", z, None), ("i", c, None))]
mcount = {}
for fn, cls, ins, outs, where in R.get("molecule", []):
    mcount[fn] = mcount.get(fn, 0) + 1
    if fn in presets:
        check(cls == "A" and ins == [] and same_toks(outs, mol_toks(presets[fn])), f"{where}: {fn}")
        continue
    pairs = mol_in(ins)
    check(all(0 <= z <= 255 and 0 <= c <= 0xFFFFFFFF for z, c in pairs), f"{where}: (u8, u32) inputs")
    bad = any(z not in ELEMS for z, _ in pairs)
    if fn == "Molecule::new":
        check(cls == "A" and same_toks(outs, mol_toks(pairs)), f"{where}: new")
    elif fn == "Molecule::molecular_weight":
        check(cls == "B", f"{where}: class {cls}")
        if bad:
            check(outs == [("E", "InvalidElement", None)], f"{where}: molecular_weight error")
        else:
            total = 0.0
            for z, c in pairs:
                total += ELEMS[z][3] * float(c)   # Rust: total += mass * count as f64, in atom order
            check(outs[0][0] == "f" and outs[0][2] == f2b(total), f"{where}: molecular_weight = {outs[0][1]!r}, recomputed {total!r}")
    elif fn == "Molecule::formula":
        check(cls == "A", f"{where}: class {cls}")
        if bad:
            check(outs == [("E", "InvalidElement", None)], f"{where}: formula error")
        else:
            want = "".join(ELEMS[z][1] + (str(c) if c > 1 else "") for z, c in pairs if c != 0)
            check(outs == [("s", want, None)], f"{where}: formula {outs[0][1]!r}, want {want!r}")
    elif fn == "Molecule::total_atoms":
        check(cls == "A" and outs == [("i", min(sum(c for _, c in pairs), 0xFFFFFFFF), None)], f"{where}: total_atoms")
    else:
        check(False, f"{where}: molecule.tsv function {fn} has no rule here")
print(f"molecule rows checked: {sum(mcount.values())} ({len(presets)} presets parsed from molecule.rs)")

# ── 8. M2a: reaction, thermo, phase, gas ────────────────────────────────────────────
def recorded_platform(name):
    head = (gold / f"{name}.tsv").read_text().split("\n", 2)[1]
    return head.split("libm: ", 1)[1].split(". ", 1)[0]

def running_platform():
    lib, ver = platform.libc_ver()
    arch = platform.machine()
    cpu = "generic"
    if arch == "x86_64":
        flags = set()
        try:
            for line in open("/proc/cpuinfo"):
                if line.startswith("flags"):
                    flags = set(line.split(":", 1)[1].split())
                    break
        except OSError:
            pass
        cpu = "fma" if {"fma", "avx2"} <= flags else "no-fma"
    return f"{lib} {ver} {arch} {cpu}"

LIBM_EXACT = all(recorded_platform(n) == running_platform() for n in ("reaction", "phase", "gas", "thermochem", "solution", "electrochemistry",
                                                                     "kinetics"))
print(f"M2a class-C recompute: {'bit-exact' if LIBM_EXACT else 'shape only'} (corpus libm "
      f"{recorded_platform('gas')}, running {running_platform()})")

class Err:
    def __init__(self, name): self.name = name
    def __eq__(self, o): return isinstance(o, Err) and o.name == self.name
    def __repr__(self): return f"Err({self.name})"

def dv(a, b):
    """IEEE division: Python raises on a zero divisor, where Rust gives ±∞ or NaN."""
    if b == 0.0:
        if a == 0.0 or math.isnan(a): return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b
def is_finite(x): return math.isfinite(x)
def is_positive(x): return x > 0.0 and math.isfinite(x)
def is_non_negative(x): return x >= 0.0 and math.isfinite(x)
def r_max(x, y):
    if math.isnan(x): return y
    if math.isnan(y): return x
    return x if x > y else y

def cbrt_cr(x):
    """cbrt(x), correctly rounded (Rust's compiler-builtins cbrt is CORE-MATH)."""
    if math.isnan(x) or x == 0 or math.isinf(x): return x
    if NUDGE[0] is None and not CORRECT[0]:
        CBRT_ARGS.append(x)
    y = math.cbrt(x)           # glibc's cbrt is not correctly rounded (seen 2 ulp off): search out
    up = down = y
    for _ in range(16):
        for cand in (up, down):
            if cbrt_correctly_rounded(x, cand): return cand
        up, down = nextup(up), nextdown(down)
    raise AssertionError(f"cbrt({x!r}): no correctly rounded value within 16 ulp of {y!r}")

# The class-C primitives (glibc through Python's math). Each call is counted and logged (its kind,
# arguments and result), and the call NUDGE[0] gets its result moved one ulp in direction NUDGE[1]
# (the sensitivity mode). The kinds name the port's correctly rounded counterpart (ADR 0006, M3c-3): "exp"
# _km_cr_exp, "expm1" _km_cr_expm1, "ln" _km_cr_ln, "log10" _km_cr_log10, "pow" _km_cr_pow, "pow10"
# _km_cr_exp10 (= a correctly rounded pow(10, y)), "pow25" _km_cr_pow(x, 2.5), "cos" _km_cos (correctly
# rounded for |x| ≤ 8, which every Peng–Robinson argument is) and "acos" _km_acos.
# With CORRECT[0] set (the corrected oracle, ADR 0006), every result is replaced by its correctly rounded
# value from MPFR before it is logged or used. cbrt is not logged: the oracle's (Rust's) is correctly rounded.
CALLS = [0]
NUDGE = [None]
CORRECT = [False]
CALL_LOG = []
# The MPFR function each kind rounds; pow10's arguments are (10.0, y) and pow25's (x, 2.5).
KIND_MPFR = {"exp": "exp", "expm1": "expm1", "ln": "log", "log10": "log10", "pow": "pow", "pow10": "pow",
             "pow25": "pow", "cos": "cos", "acos": "acos"}
def _prim(y, kind, args):
    k = CALLS[0]
    CALLS[0] += 1
    if CORRECT[0]:
        y = cr_of(KIND_MPFR[kind], args)
    CALL_LOG.append((kind, args, y))
    if NUDGE[0] is not None and NUDGE[0][0] == k and math.isfinite(y) and y != 0:
        return math.nextafter(y, math.inf if NUDGE[0][1] > 0 else -math.inf)
    return y
def p_exp(x):
    if math.isnan(x): return _prim(x, "exp", (x,))
    try: return _prim(math.exp(x), "exp", (x,))
    except OverflowError: return _prim(math.inf, "exp", (x,))
def p_ln(x):
    if math.isnan(x) or x < 0: return _prim(math.nan, "ln", (x,))
    if x == 0: return _prim(-math.inf, "ln", (x,))
    return _prim(math.log(x), "ln", (x,))
def _pow(x, y, kind):
    try: return _prim(math.pow(x, y), kind, (x, y))
    except OverflowError: return _prim(math.inf, kind, (x, y))
    except ValueError: return _prim(math.nan, kind, (x, y))
def p_pow(x, y): return _pow(x, y, "pow10" if x == 10.0 else "pow")
def p_pow25(x): return _pow(x, 2.5, "pow25")
def p_acos(x):
    if math.isnan(x): return _prim(x, "acos", (x,))
    return _prim(math.acos(x), "acos", (x,))
def p_cos(x):
    if not math.isfinite(x): return _prim(math.nan, "cos", (x,))
    return _prim(math.cos(x), "cos", (x,))

# The sensitivity mode's nudges (a diagnostic of conditioning since ADR 0006): every logged result but one that is
# 0, ±∞ or NaN, or acos(±1) (0 and π rounded, by C99 and in the port). The port's primitives are all correctly
# rounded, so its result differs from the oracle's exactly where glibc's is not correctly rounded (MPFR decides;
# glibc_not_cr), and from the corrected oracle's nowhere.
def glibc_not_cr(kind, args, y):
    return not same_f(cr_of(KIND_MPFR[kind], args), y)
def nudged(kind, args, y):
    if not math.isfinite(y) or y == 0: return False
    if kind == "acos" and abs(args[0]) == 1.0: return False
    return True

R_GAS = b2f(golden_consts["element::GAS_CONSTANT"])
TAU = 2 * math.pi              # std::f64::consts::TAU (doubling is exact)
SQRT2 = math.sqrt(2.0)         # std::f64::consts::SQRT_2 (sqrt is correctly rounded)

# reaction.rs
def r_gibbs(h, t, s_):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(h) or not is_finite(s_): return Err("InvalidInput")
    return h - t * s_
def r_keq(dg, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(dg): return Err("InvalidInput")
    k = p_exp(dv(-dg, R_GAS * t))
    if not is_finite(k): return Err("ComputationError")
    return k

# thermo.rs
def r_sum(xs):
    acc = -0.0                 # Rust's f64 Sum folds from -0.0
    for x in xs: acc = acc + x
    return acc
def r_calorimetry(e, m, dt):
    if not is_finite(e): return Err("InvalidInput")
    if not is_positive(m): return Err("InvalidInput")
    if not is_finite(dt) or abs(dt) < 1e-12: return Err("InvalidInput")   # hisab 1.4.0 EPSILON_F64
    return dv(e, m * dt)

# phase.rs
psrc = (root / "rust-old/src/phase.rs").read_text()
PA_PER_MMHG = float(re.search(r"^const PA_PER_MMHG: f64 = ([0-9_.]+);", psrc, re.M).group(1).replace("_", ""))
ANTOINE = []
for name, body in re.findall(r'\(\s*"(\w+)",\s*(?://[^\n]*\s*)*AntoineCoeffs \{(.*?)\}', psrc[psrc.index("pub static ANTOINE_DATA"):], re.S):
    fl = dict(re.findall(r"(\w+): (-?[0-9.]+),", body))
    ANTOINE.append((name, [float(fl[k]) for k in ("a", "b", "c", "t_min_c", "t_max_c")]))
check(len(ANTOINE) == 8, f"phase.rs: ANTOINE_DATA parsed as {len(ANTOINE)} rows")
def unit_interval(x): return 0.0 <= x <= 1.0
def r_raoult(x, p):
    if not unit_interval(x): return Err("InvalidInput")
    if not is_non_negative(p): return Err("InvalidInput")
    return x * p
def r_binary(x, pa, pb):
    if not unit_interval(x): return Err("InvalidInput")
    if not is_non_negative(pa) or not is_non_negative(pb): return Err("InvalidInput")
    return x * pa + (1.0 - x) * pb
def r_henry_p(kh, x):
    if not is_positive(kh): return Err("InvalidInput")
    if not unit_interval(x): return Err("InvalidInput")
    return kh * x
def r_henry_s(p, kh):
    if not is_positive(kh): return Err("InvalidInput")
    if not is_non_negative(p): return Err("InvalidInput")
    return dv(p, kh)
def r_osmotic(i, c, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_positive(i): return Err("InvalidInput")
    if not is_non_negative(c): return Err("InvalidInput")
    return i * c * R_GAS * 1000.0 * t
def r_clausius(p1, t1, t2, dh):
    if not is_positive(p1): return Err("InvalidInput")
    if not is_positive(t1) or not is_positive(t2): return Err("InvalidTemperature")
    if not is_finite(dh): return Err("InvalidInput")
    e = -dv(dh, R_GAS) * (dv(1.0, t2) - dv(1.0, t1))
    return p1 * p_exp(e)
def r_antoine_mmhg(a, b, c, tmin, tmax, t):
    denom = c + t
    if not is_positive(denom) or not is_finite(a) or not is_finite(b): return math.nan
    return p_pow(10.0, a - dv(b, denom))
def r_antoine_pa(a, b, c, tmin, tmax, t): return r_antoine_mmhg(a, b, c, tmin, tmax, t) * PA_PER_MMHG
def r_phase_rule(c, p):
    if p == 0: return Err("InvalidInput")
    f_ = c - p + 2
    if f_ < 0 or f_ > 0xFFFFFFFF: return Err("InvalidInput")
    return f_

# gas.rs
gsrc = (root / "rust-old/src/gas.rs").read_text()
PRESETS = {}
for name, tc, pc, om in re.findall(r"pub const (PR_\w+): PengRobinsonParams = PengRobinsonParams \{\s*tc: ([0-9_.]+),\s*pc: ([0-9_.]+),\s*omega: ([0-9_.]+),\s*\};", gsrc):
    PRESETS[name] = [float(v.replace("_", "")) for v in (tc, pc, om)]
check(len(PRESETS) == 5, f"gas.rs: PR_* parsed as {sorted(PRESETS)}")
def r_ideal_p(n, t, v):
    if not is_finite(n): return Err("InvalidInput")
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_positive(v): return Err("InvalidInput")
    return dv(n * R_GAS * t, v)
def r_ideal_v(n, t, p):
    if not is_finite(n): return Err("InvalidInput")
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_positive(p): return Err("InvalidInput")
    return dv(n * R_GAS * t, p)
def r_vdw(n, t, v, a, b):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(n): return Err("InvalidInput")
    if not is_non_negative(a) or not is_non_negative(b): return Err("InvalidInput")
    if not is_positive(v): return Err("InvalidInput")
    if v <= n * b: return Err("InvalidInput")
    nv = dv(n, v)
    return dv(n * R_GAS * t, v - n * b) - a * nv * nv
def r_density(p, mm, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_non_negative(p): return Err("InvalidInput")
    if not is_positive(mm): return Err("InvalidInput")
    return dv(p * mm, R_GAS * t)
def pr_params_err(tc, pc, om):
    if not is_positive(tc): return Err("InvalidTemperature")
    if not is_positive(pc): return Err("InvalidInput")
    if not is_finite(om): return Err("InvalidInput")
    return None
def pr_coefficients(tc, pc, om, t):
    a = dv(0.45724 * R_GAS * R_GAS * tc * tc, pc)
    b = dv(0.07780 * R_GAS * tc, pc)
    kappa = 0.37464 + 1.54226 * om - 0.26992 * om * om
    tr = dv(t, tc)
    u = 1.0 + kappa * (1.0 - math.sqrt(tr))
    return a * (u * u), b              # powi(u, 2) is u·u in __powidf2
def r_pr_pressure(tc, pc, om, t, v):
    if not is_positive(t): return Err("InvalidTemperature")
    e = pr_params_err(tc, pc, om)
    if e: return e
    aa, b = pr_coefficients(tc, pc, om, t)
    if not is_finite(v) or v <= b: return Err("InvalidInput")
    return dv(R_GAS * t, v - b) - dv(aa, v * v + 2.0 * b * v - b * b)
def pr_dimensionless(tc, pc, om, t, p):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_positive(p): return Err("InvalidInput")
    e = pr_params_err(tc, pc, om)
    if e: return e
    aa, b = pr_coefficients(tc, pc, om, t)
    rt = R_GAS * t
    return (dv(aa * p, rt * rt), dv(b * p, rt))
def cubic_real_roots(c2, c1, c0):
    shift = c2 / 3.0
    q = (c2 * c2 - 3.0 * c1) / 9.0
    r = (c2 * (2.0 * c2 * c2 - 9.0 * c1) + 27.0 * c0) / 54.0   # (non-zero divisors)
    q3 = q * q * q
    r2 = r * r
    if q > 0.0 and r2 <= q3 * (1.0 + 1e-12):
        sq = math.sqrt(q)
        x = dv(r, sq * q)
        x = -1.0 if x < -1.0 else (1.0 if x > 1.0 else x)   # clamp(-1, 1); NaN passes
        theta = p_acos(x)
        m = -2.0 * sq
        roots = [m * p_cos(theta / 3.0) - shift, m * p_cos((theta - TAU) / 3.0) - shift,
                 m * p_cos((theta + TAU) / 3.0) - shift]
        n = 3
    else:
        sig = math.nan if math.isnan(r) else math.copysign(1.0, r)
        big_a = -sig * cbrt_cr(abs(r) + math.sqrt(r_max(r2 - q3, 0.0)))
        big_b = 0.0 if big_a == 0.0 else dv(q, big_a)
        roots = [big_a + big_b - shift, 0.0, 0.0]
        n = 1
    poly = lambda z: ((z + c2) * z + c1) * z + c0
    for k in range(n):
        z = roots[k]
        for _ in range(2):
            fz = poly(z)
            df = (3.0 * z + 2.0 * c2) * z + c1
            if fz == 0.0 or df == 0.0: break
            nxt = z - dv(fz, df)
            if abs(poly(nxt)) < abs(fz): z = nxt
            else: break
        roots[k] = z
    if n == 3:
        if roots[0] > roots[1]: roots[0], roots[1] = roots[1], roots[0]
        if roots[1] > roots[2]: roots[1], roots[2] = roots[2], roots[1]
        if roots[0] > roots[1]: roots[0], roots[1] = roots[1], roots[0]
    return roots, n
def pr_physical(a_, b_):
    c2 = b_ - 1.0
    c1 = a_ - 3.0 * b_ * b_ - 2.0 * b_
    c0 = -(a_ * b_ - b_ * b_ - b_ * b_ * b_)
    roots, n = cubic_real_roots(c2, c1, c0)
    return [z for z in roots[:n] if math.isfinite(z) and z > b_]
def r_pr_roots(tc, pc, om, t, p):
    d = pr_dimensionless(tc, pc, om, t, p)
    if isinstance(d, Err): return d
    phys = pr_physical(*d)
    return phys if phys else Err("ComputationError")
def r_pr_z(tc, pc, om, t, p):
    rr = r_pr_roots(tc, pc, om, t, p)
    return rr if isinstance(rr, Err) else rr[-1]
def r_fug_coeff(tc, pc, om, t, p):
    d = pr_dimensionless(tc, pc, om, t, p)
    if isinstance(d, Err): return d
    a_, b_ = d
    phys = pr_physical(a_, b_)
    if not phys: return Err("ComputationError")
    z = phys[-1]
    ln_phi = (z - 1.0) - p_ln(z - b_) - dv(a_, 2.0 * SQRT2 * b_) * p_ln(dv(z + (1.0 + SQRT2) * b_, z + (1.0 - SQRT2) * b_))
    phi = p_exp(ln_phi)
    if not is_positive(phi): return Err("ComputationError")
    return phi
def r_fugacity(tc, pc, om, t, p):
    phi = r_fug_coeff(tc, pc, om, t, p)
    return phi if isinstance(phi, Err) else phi * p
def r_rk(tc, pc, t, v):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_positive(tc): return Err("InvalidTemperature")
    if not is_positive(pc): return Err("InvalidInput")
    a = dv(0.42748 * R_GAS * R_GAS * p_pow25(tc), pc)
    b = dv(0.08664 * R_GAS * tc, pc)
    if not is_finite(v) or v <= b: return Err("InvalidInput")
    return dv(R_GAS * t, v - b) - dv(a, math.sqrt(t) * v * (v + b))

def out_toks(x):
    """The golden tokens of a Rust result: Err, a root list (`v:n` + floats), an int, or a float."""
    if isinstance(x, Err): return [("E", x.name, None)]
    if isinstance(x, list): return [("v", len(x), None)] + [("f", z, f2b(z)) for z in x]
    if isinstance(x, tuple): return [("f", z, f2b(z)) for z in x]     # a Rust tuple: its fields in order (M3b)
    if isinstance(x, bool): return [("i", int(x), None)]
    if isinstance(x, int): return [("i", x, None)]
    return [("f", x, f2b(x))]
def toks_equal(got, want):
    if len(got) != len(want): return False
    for (gk, gv, gb), (wk, wv, wb) in zip(got, want):
        if gk != wk: return False
        if gk == "f":
            if not (gb == wb or (isnan_b(gb) and isnan_b(wb))): return False
        elif gv != wv: return False
    return True
def toks_shape(got, want):
    return len(got) == len(want) and all(g[0] == w[0] and (g[0] == "f" or g[1] == w[1]) for g, w in zip(got, want))

def fl(ins): return [t[1] for t in ins]
def slice_args(ins):
    """[v:n f..][v:m f..] → the lists."""
    out, k = [], 0
    while k < len(ins):
        n = ins[k][1]
        out.append([t[1] for t in ins[k + 1:k + 1 + n]])
        k += 1 + n
    return out

# thermochem.rs (M2b)
tsrc = (root / "rust-old/src/thermochem.rs").read_text()
def _literal(name):
    body = tsrc[tsrc.index(f"pub static {name}:"):]
    return body[:body.index("\n];")]
TC = [(f, n, [float(hf), float(gf), float(s0)]) for f, n, hf, gf, s0 in re.findall(
    r'ThermochemData \{\s*formula: "([^"]*)",\s*name: "([^"]*)",\s*delta_hf_kj: (-?[0-9.]+),\s*'
    r'delta_gf_kj: (-?[0-9.]+),\s*(?://[^\n]*\s*)*s_standard_j: (-?[0-9.]+),\s*\}', _literal("THERMOCHEM_DATA"))]
check(len(TC) == 76, f"thermochem.rs: THERMOCHEM_DATA parsed as {len(TC)} rows")
SH_FIELDS = ("a", "b", "c", "d", "e", "t_min", "t_max")
SH = []
for f_, body in re.findall(r'\(\s*"([^"]+)",\s*Shomate \{(.*?)\},\s*\)', _literal("SHOMATE_DATA"), re.S):
    fd = dict(re.findall(r"(\w+): (-?[0-9.]+),", body))
    SH.append((f_, [float(fd[k]) for k in SH_FIELDS]))
check(len(SH) == 34, f"thermochem.rs: SHOMATE_DATA parsed as {len(SH)} rows")
def tc_lookup(formula):
    return next((r for r in TC if r[0] == formula), None)
def sh_lookup(formula):
    return next((c for f_, c in SH if f_ == formula), None)
def sh_lookup_at(formula, t):
    return next((c for f_, c in SH if f_ == formula and c[5] <= t <= c[6]), None)
def r_min(x, y):
    if math.isnan(x): return y
    if math.isnan(y): return x
    return x if x < y else y
def sh_cp(c, temp):
    a, b, c_, d, e = c[:5]
    t = temp / 1000.0
    return a + b * t + c_ * t * t + d * t * t * t + dv(e, t * t)
def sh_h(c, temp):
    a, b, c_, d, e = c[:5]
    t = temp / 1000.0
    return 1000.0 * (t * (a + t * (b / 2.0 + t * (c_ / 3.0 + t * (d / 4.0)))) - dv(e, t))
def segs(formula):
    """Rust's private Segments::lookup: the adjacent run of `formula`'s rows, as [coefficients]."""
    k = next((k for k, (f_, _) in enumerate(SH) if f_ == formula), None)
    if k is None: return None
    run = []
    while k < len(SH) and SH[k][0] == formula:
        run.append(SH[k][1])
        k += 1
    return run
def seg_cp(run, t):
    for c in run:
        if t <= c[6]: return sh_cp(c, t)
    return sh_cp(run[-1], t)
def seg_dh(run, t1, t2):
    lo, hi, sign = (t1, t2, 1.0) if t1 <= t2 else (t2, t1, -1.0)
    total, lower = 0.0, -math.inf
    for k, c in enumerate(run):
        upper = math.inf if k + 1 == len(run) else c[6]
        a, b = r_max(lo, lower), r_min(hi, upper)
        if b > a: total += sh_h(c, b) - sh_h(c, a)
        if upper >= hi: break
        lower = upper
    return sign * total
def r_sum_property(species, col):
    total = 0.0                # sum_property starts from +0.0 and uses `+=` (not .sum())
    for f_, c in species:
        row = tc_lookup(f_)
        if row is None: return Err("InvalidReaction")
        total += c * row[2][col]
    return total
def r_reaction(prod, reac, col, factor):
    sp_, sr = r_sum_property(prod, col), None
    if isinstance(sp_, Err): return sp_
    sr = r_sum_property(reac, col)
    if isinstance(sr, Err): return sr
    d = sp_ - sr
    return d * 1000.0 if factor else d     # ADR 0002 #1: the port's J/mol, the oracle's kJ × 1000.0
def r_vant_hoff(k, dh, t1, t2):
    if not is_positive(t1) or not is_positive(t2): return Err("InvalidTemperature")
    if not is_positive(k): return Err("InvalidInput")
    if not is_finite(dh): return Err("InvalidInput")
    return k * p_exp(-dv(dh, R_GAS) * (dv(1.0, t2) - dv(1.0, t1)))
def r_dh_cp(formula, t1, t2):
    run = segs(formula)
    if run is None: return Err("InvalidReaction")
    if not is_positive(t1) or not is_positive(t2): return Err("InvalidTemperature")
    lo, hi = r_min(t1, t2), r_max(t1, t2)
    if lo < run[0][5] or hi > run[-1][6]: return Err("InvalidTemperature")
    return seg_dh(run, t1, t2)
FLAME_PATH = [None]
def r_flame(dh, prods, t0, ntol=1e-3, btol=1e-6, gcp=40.0, eps=1e-12):
    """adiabatic_flame_temperature with hisab 1.4.0's newton_raphson (tol ntol, 50 iterations,
    derivative threshold eps) and bisection (tol btol, 200 iterations) transcribed; FLAME_PATH[0] is
    the solver that answered. The keyword arguments are Rust's typed constants (the pinning check)."""
    FLAME_PATH[0] = None
    if not is_positive(t0): return Err("InvalidTemperature")
    if not is_finite(dh): return Err("InvalidInput")
    if not prods: return Err("InvalidInput")
    res, rmin, rmax = [], -math.inf, math.inf
    for f_, n in prods:
        run = segs(f_)
        if run is None: return Err("InvalidReaction")
        if not is_positive(n): return Err("InvalidInput")
        rmin, rmax = r_max(rmin, run[0][5]), r_min(rmax, run[-1][6])
        res.append((run, n))
    if t0 < rmin or t0 >= rmax: return Err("InvalidTemperature")
    def fun(t):
        heat = -0.0                # Rust's .sum() folds from -0.0
        for run, n in res: heat += n * seg_dh(run, t0, t)
        return -dh - heat
    def dfun(t):
        cp = -0.0
        for run, n in res: cp += n * seg_cp(run, t)
        return -cp
    t_lo, t_hi = t0, rmax
    if not fun(t_lo) * fun(t_hi) < 0.0: return Err("ComputationError")
    total = -0.0
    for _, n in res: total += n
    x = r_min(r_max(t0 - dv(dh, total * gcp), t_lo), t_hi)
    why = "no convergence in 50 iterations"
    for k in range(50):
        fx = fun(x)
        if abs(fx) < ntol:
            if t_lo <= x <= t_hi:
                FLAME_PATH[0] = f"newton, {k} steps"
                return x
            why = "root outside the range"
            break
        dfx = dfun(x)
        if abs(dfx) < eps:
            why = "zero derivative"
            break
        x -= dv(fx, dfx)
    lo, hi = t_lo, t_hi
    f_lo, f_hi = fun(lo), fun(hi)
    if f_lo * f_hi > 0.0: return Err("ComputationError")
    if f_lo > 0.0: lo, hi = hi, lo
    for k in range(200):
        mid = (lo + hi) * 0.5
        fm = fun(mid)
        if abs(fm) < btol or abs(hi - lo) < btol:
            FLAME_PATH[0] = f"bisection after newton: {why}, {k} halvings"
            return mid
        if fm < 0.0: lo = mid
        else: hi = mid
    FLAME_PATH[0] = f"bisection after newton: {why}, 200 halvings"
    return (lo + hi) * 0.5
def species_at(ins, k):
    """The species list at token k (`v:n`, then n `s:` / `f:` pairs) and the index after it."""
    n = ins[k][1]
    return [(ins[k + 1 + 2 * j][1], ins[k + 2 + 2 * j][1]) for j in range(n)], k + 1 + 2 * n
def reaction_rule(col, factor):
    def rule(a):
        prod, k = species_at(a, 0)
        reac, _ = species_at(a, k)
        return r_reaction(prod, reac, col, factor)
    return rule
def flame_rule(a):
    prods, k = species_at(a, 1)
    return r_flame(a[0][1], prods, a[k][1])
def shomate_rule(fn):
    def rule(a):
        v = fl(a)
        return sh_cp(v[:7], v[7]) if fn == "cp" else sh_h(v[:7], v[7])
    return rule


# solution.rs (M3a)
KW = b2f(golden_consts["solution::KW_WATER_25C"])
def p_log10(x):
    """glibc's log10 (Python's math); the port's is the correctly rounded _km_cr_log10, the corrected oracle's."""
    if math.isnan(x) or x < 0: return _prim(math.nan, "log10", (x,))
    if x == 0: return _prim(-math.inf, "log10", (x,))
    return _prim(math.log10(x), "log10", (x,))
def r_molarity(n, v):
    if not is_positive(v): return Err("InvalidInput")
    if not is_finite(n): return Err("InvalidInput")
    return dv(n, v)
def r_dilution(m1, v1, v2):
    if not is_positive(v2): return Err("InvalidInput")
    if not is_finite(m1) or not is_finite(v1): return Err("InvalidInput")
    return dv(m1 * v1, v2)
def r_ph_h(h):
    if not is_positive(h): return Err("InvalidConcentration")
    return -p_log10(h)
def r_hh(pka, base, acid):
    if not is_positive(acid) or not is_positive(base): return Err("InvalidConcentration")
    if not is_finite(pka): return Err("InvalidInput")
    return pka + p_log10(dv(base, acid))
ROOT_PATH = [None]
def r_root(k, c, tol=1e-12, max_iter=50):
    """Rust's private charge_balance_root and checked_root, transcribed; ROOT_PATH[0] is the note's path."""
    def checked(x, path):
        ok = is_positive(x)
        ROOT_PATH[0] = path + ("" if ok else ", checked_root rejected")
        return x if ok else Err("ComputationError")
    kc_kw = k * c + KW
    weak_b = math.sqrt(kc_kw)
    strong = 0.5 * (c + math.sqrt(c * c + 4.0 * KW))
    x = r_min(weak_b, strong)
    for it in range(max_iter):
        fx = x * x * (x + k) - k * (c * x + KW) - KW * x
        if fx <= 0.0: return checked(x, f"root: f(x) <= 0, {it} updates")
        dfx = x * (3.0 * x + 2.0 * k) - kc_kw
        step = dv(fx, dfx)
        if not math.isfinite(step):
            ROOT_PATH[0] = f"root: non-finite step, {it} updates"
            return Err("ComputationError")
        x -= step
        if step <= tol * x: return checked(x, f"root: relative step, {it + 1} updates")
    ROOT_PATH[0] = f"root: {max_iter} iterations"
    return Err("ComputationError")
def r_weak(k, c, base, tol=1e-12, max_iter=50):
    ROOT_PATH[0] = None
    if not is_positive(k): return Err("InvalidInput")
    if not is_positive(c): return Err("InvalidConcentration")
    x = r_root(k, c, tol, max_iter)
    if isinstance(x, Err): return x
    p = -p_log10(x)
    return 14.0 - p if base else p
def r_ionic(ins):
    n = ins[0][1]
    acc = -0.0                 # Rust's f64 Sum folds from -0.0
    for j in range(n):
        c, z = ins[1 + 2 * j][1], float(ins[2 + 2 * j][1])
        acc += c * z * z
    return 0.5 * acc
def r_ia(i_, a):
    if not is_non_negative(i_): return Err("InvalidInput")
    if not is_finite(a): return Err("InvalidInput")
    return None
def r_dh_lim(z, i_, a):
    e = r_ia(i_, a)
    if e: return e
    z = float(z)
    return p_pow(10.0, -a * z * z * math.sqrt(i_))
def r_dh_ext(z, i_, a, b, size):
    e = r_ia(i_, a)
    if e: return e
    if not is_non_negative(b) or not is_non_negative(size): return Err("InvalidInput")
    sq, z = math.sqrt(i_), float(z)
    return p_pow(10.0, dv(-a * z * z * sq, 1.0 + b * size * sq))
def r_davies(z, i_, a):
    e = r_ia(i_, a)
    if e: return e
    sq, z = math.sqrt(i_), float(z)
    return p_pow(10.0, -a * z * z * (dv(sq, 1.0 + sq) - 0.3 * i_))
def r_kohlrausch(l0, k, c):
    if not is_non_negative(c): return Err("InvalidConcentration")
    if not is_finite(l0) or not is_finite(k): return Err("InvalidInput")
    return l0 - k * math.sqrt(c)
def weak_rule(base):
    return lambda a: r_weak(*fl(a), base)

# electrochemistry.rs (M3a)
ecsrc = (root / "rust-old/src/electrochemistry.rs").read_text()
_ec_body = ecsrc[ecsrc.index("pub static STANDARD_POTENTIALS:"):]
_ec_body = _ec_body[:_ec_body.index("\n];")]
SP = [(eq, cp, int(n), float(e)) for eq, cp, n, e in re.findall(
    r'HalfReaction \{\s*equation: "([^"]*)",\s*couple: "([^"]*)",\s*electrons: (\d+),\s*'
    r'standard_potential: (-?[0-9.]+),(?:\s*//[^\n]*)?\s*\}', _ec_body)]
check(len(SP) == 49, f"electrochemistry.rs: STANDARD_POTENTIALS parsed as {len(SP)} rows")
FARADAY = b2f(golden_consts["electrochemistry::FARADAY"])
def sp_lookup(couple):
    return next((r for r in SP if r[1] == couple), None)
def sp_toks(r):
    return [("s", r[0], None), ("s", r[1], None), ("i", r[2], None), ("f", r[3], f2b(r[3]))]
def r_nernst(e0, n, t, q):
    if not is_positive(t): return Err("InvalidTemperature")
    if n == 0: return Err("InvalidInput")
    if not is_positive(q): return Err("InvalidInput")
    if not is_finite(e0): return Err("InvalidInput")
    return e0 - dv(R_GAS * t, float(n) * FARADAY) * p_ln(q)
def r_cell_couples(c, a):
    hc = sp_lookup(c)
    if hc is None: return Err("InvalidReaction")
    ha = sp_lookup(a)
    if ha is None: return Err("InvalidReaction")
    return hc[3] - ha[3]
def r_cell_gibbs(n, e):
    if n == 0: return Err("InvalidInput")
    if not is_finite(e): return Err("InvalidInput")
    return -float(n) * FARADAY * e
def r_faraday_mass(q, mm, n):
    if n == 0: return Err("InvalidInput")
    if not is_positive(mm): return Err("InvalidInput")
    if not is_finite(q): return Err("InvalidInput")
    return dv(q * mm, float(n) * FARADAY)
def r_faraday_charge(ma, mm, n):
    if n == 0: return Err("InvalidInput")
    if not is_positive(mm): return Err("InvalidInput")
    if not is_finite(ma): return Err("InvalidInput")
    return dv(ma * float(n) * FARADAY, mm)
def r_faraday_ratio(ma, mma, na, mmb, nb):
    if na == 0 or nb == 0: return Err("InvalidInput")
    if not is_positive(mma): return Err("InvalidInput")
    if not is_positive(mmb): return Err("InvalidInput")
    if not is_finite(ma): return Err("InvalidInput")
    return dv(ma * dv(mmb, float(nb)), dv(mma, float(na)))

# kinetics.rs (M3b)
LN2 = b2f(0x3FE62E42FEFA39EF)        # std::f64::consts::LN_2, typed inline in src/kinetics.cyr
BOLTZMANN = b2f(golden_consts["kinetics::BOLTZMANN"])
PLANCK = b2f(golden_consts["spectroscopy::PLANCK"])
def p_expm1(x):
    """glibc's expm1 (Python's math); the port's is the correctly rounded _km_cr_expm1, the corrected oracle's."""
    if math.isnan(x): return _prim(x, "expm1", (x,))
    try: return _prim(math.expm1(x), "expm1", (x,))
    except OverflowError: return _prim(math.inf, "expm1", (x,))
def r_arrhenius(a, ea, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(a): return Err("InvalidInput")
    if not is_finite(ea): return Err("InvalidInput")
    return a * p_exp(dv(-ea, R_GAS * t))
def r_half_first(k, ln2=LN2):
    if not is_positive(k): return Err("InvalidInput")
    return dv(ln2, k)
def r_second(c0, k, t):
    if not is_positive(c0): return Err("InvalidConcentration")
    if not is_non_negative(k): return Err("InvalidInput")
    if not is_non_negative(t): return Err("InvalidInput")
    return dv(1.0, dv(1.0, c0) + k * t)
def r_zero(c0, k, t):
    remaining = c0 - k * t
    return 0.0 if remaining < 0.0 else remaining
def r_zero_half(c0, k):
    if not is_positive(k): return Err("InvalidInput")
    if not is_positive(c0): return Err("InvalidConcentration")
    return dv(c0, 2.0 * k)
def r_nth(c0, k, n, tol=1e-12, ln2=LN2):
    if not is_finite(n) or abs(n - 1.0) <= tol: return Err("InvalidInput")
    if not is_positive(k): return Err("InvalidInput")
    if not is_positive(c0): return Err("InvalidConcentration")
    nm1 = n - 1.0
    num = p_expm1(nm1 * ln2)
    t = dv(num, nm1 * k * p_pow(c0, nm1))
    return t if is_positive(t) else Err("ComputationError")
def r_ssi(k1, a, k2):
    if not is_positive(k2): return Err("InvalidInput")
    if not is_non_negative(k1): return Err("InvalidInput")
    if not is_non_negative(a): return Err("InvalidConcentration")
    return dv(k1 * a, k2)
def r_pre(kf, kr, ks, a):
    if not is_positive(kr): return Err("InvalidInput")
    if not is_non_negative(kf) or not is_non_negative(ks): return Err("InvalidInput")
    if not is_non_negative(a): return Err("InvalidConcentration")
    return ks * dv(kf, kr) * a
def r_mm(v, km, s_):
    if not is_positive(v): return Err("InvalidInput")
    if not is_positive(km): return Err("InvalidInput")
    if not is_non_negative(s_): return Err("InvalidConcentration")
    return dv(v * s_, km + s_)
def r_lb(v, km):
    if not is_positive(v): return Err("InvalidInput")
    if not is_positive(km): return Err("InvalidInput")
    return (dv(km, v), dv(1.0, v))
def r_first(c0, k, t):
    return c0 * p_exp(-k * t)
def r_collision(z, p, ea, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not (p > 0.0 and p <= 1.0): return Err("InvalidInput")
    if not is_positive(z): return Err("InvalidInput")
    if not is_finite(ea): return Err("InvalidInput")
    return z * p * p_exp(dv(-ea, R_GAS * t))
def r_eyring(dg, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(dg): return Err("InvalidInput")
    prefactor = dv(BOLTZMANN * t, PLANCK)
    return prefactor * p_exp(dv(-dg, R_GAS * t))
def r_eyring_act(dh, ds, t):
    if not is_positive(t): return Err("InvalidTemperature")
    if not is_finite(dh): return Err("InvalidInput")
    if not is_finite(ds): return Err("InvalidInput")
    prefactor = dv(BOLTZMANN * t, PLANCK)
    es = p_exp(dv(ds, R_GAS))
    return prefactor * es * p_exp(dv(-dh, R_GAS * t))


RULES = {
    "reaction::gibbs_free_energy": ("B", lambda a: r_gibbs(*fl(a))),
    "reaction::equilibrium_constant": ("C", lambda a: r_keq(*fl(a))),
    "reaction::is_spontaneous": ("A", lambda a: fl(a)[0] < 0.0),
    "thermo::heat_transfer": ("B", lambda a: fl(a)[0] * fl(a)[1] * fl(a)[2]),
    "thermo::hess_law": ("B", lambda a: r_sum(slice_args(a)[0])),
    "thermo::enthalpy_from_formation": ("B", lambda a: r_sum(slice_args(a)[0]) - r_sum(slice_args(a)[1])),
    "thermo::specific_heat_from_calorimetry": ("B", lambda a: r_calorimetry(*fl(a))),
    "phase::raoult_pressure": ("B", lambda a: r_raoult(*fl(a))),
    "phase::raoult_total_binary": ("B", lambda a: r_binary(*fl(a))),
    "phase::henry_pressure": ("B", lambda a: r_henry_p(*fl(a))),
    "phase::henry_solubility": ("B", lambda a: r_henry_s(*fl(a))),
    "phase::boiling_point_elevation": ("B", lambda a: fl(a)[0] * fl(a)[1] * fl(a)[2]),
    "phase::freezing_point_depression": ("B", lambda a: fl(a)[0] * fl(a)[1] * fl(a)[2]),
    "phase::osmotic_pressure": ("B", lambda a: r_osmotic(*fl(a))),
    "phase::clausius_clapeyron": ("C", lambda a: r_clausius(*fl(a))),
    "phase::antoine_pressure_mmhg": ("C", lambda a: r_antoine_mmhg(*fl(a))),
    "phase::antoine_pressure_pa": ("C", lambda a: r_antoine_pa(*fl(a))),
    "phase::phase_rule": ("A", lambda a: r_phase_rule(*fl(a))),
    "gas::ideal_gas_pressure": ("B", lambda a: r_ideal_p(*fl(a))),
    "gas::ideal_gas_volume": ("B", lambda a: r_ideal_v(*fl(a))),
    "gas::van_der_waals_pressure": ("B", lambda a: r_vdw(*fl(a))),
    "gas::partial_pressure": ("B", lambda a: fl(a)[0] * fl(a)[1]),
    "gas::gas_density": ("B", lambda a: r_density(*fl(a))),
    "gas::peng_robinson_pressure": ("B", lambda a: r_pr_pressure(*fl(a))),
    "gas::peng_robinson_z": ("C", lambda a: r_pr_z(*fl(a))),
    "gas::peng_robinson_z_roots": ("C", lambda a: r_pr_roots(*fl(a))),
    "gas::fugacity_coefficient_pr": ("C", lambda a: r_fug_coeff(*fl(a))),
    "gas::fugacity_pr": ("C", lambda a: r_fugacity(*fl(a))),
    "gas::redlich_kwong_pressure": ("C", lambda a: r_rk(*fl(a))),
    "gas::joule_thomson_ideal": ("A", lambda a: 0.0),
    "thermochem::reaction_enthalpy": ("B", reaction_rule(0, True)),
    "thermochem::reaction_gibbs_energy": ("B", reaction_rule(1, True)),
    "thermochem::reaction_entropy": ("B", reaction_rule(2, False)),
    "thermochem::vant_hoff_k": ("C", lambda a: r_vant_hoff(*fl(a))),
    "thermochem::enthalpy_change_cp": ("B", lambda a: r_dh_cp(*fl(a))),
    "thermochem::adiabatic_flame_temperature": ("B", flame_rule),
    "Shomate::cp": ("B", shomate_rule("cp")),
    "Shomate::enthalpy_antiderivative": ("B", shomate_rule("h")),
    "solution::molarity": ("B", lambda a: r_molarity(*fl(a))),
    "solution::molality": ("B", lambda a: r_molarity(*fl(a))),
    "solution::dilution": ("B", lambda a: r_dilution(*fl(a))),
    "solution::ph_from_h_concentration": ("C", lambda a: r_ph_h(*fl(a))),
    "solution::poh_from_oh_concentration": ("C", lambda a: r_ph_h(*fl(a))),
    "solution::ph_from_poh": ("B", lambda a: 14.0 - fl(a)[0]),
    "solution::henderson_hasselbalch": ("C", lambda a: r_hh(*fl(a))),
    "solution::h_concentration_from_ph": ("C", lambda a: p_pow(10.0, -fl(a)[0])),
    "solution::weak_acid_ph": ("C", weak_rule(False)),
    "solution::weak_base_ph": ("C", weak_rule(True)),
    "solution::ionic_strength": ("B", r_ionic),
    "solution::debye_huckel_limiting": ("C", lambda a: r_dh_lim(*fl(a))),
    "solution::debye_huckel_extended": ("C", lambda a: r_dh_ext(*fl(a))),
    "solution::davies_activity": ("C", lambda a: r_davies(*fl(a))),
    "solution::kohlrausch_molar_conductivity": ("B", lambda a: r_kohlrausch(*fl(a))),
    "solution::limiting_molar_conductivity": ("B", lambda a: fl(a)[1] * fl(a)[0] + fl(a)[3] * fl(a)[2]),
    "solution::solution_conductivity": ("B", lambda a: dv(fl(a)[0] * fl(a)[1], 1000.0)),
    "electrochemistry::nernst_potential": ("C", lambda a: r_nernst(*fl(a))),
    "electrochemistry::nernst_potential_25c": ("C", lambda a: r_nernst(fl(a)[0], fl(a)[1], 298.15, fl(a)[2])),
    "electrochemistry::standard_cell_potential": ("B", lambda a: fl(a)[0] - fl(a)[1]),
    "electrochemistry::cell_potential_from_couples": ("B", lambda a: r_cell_couples(*fl(a))),
    "electrochemistry::is_spontaneous_cell": ("A", lambda a: fl(a)[0] > 0.0),
    "electrochemistry::cell_gibbs_energy": ("B", lambda a: r_cell_gibbs(*fl(a))),
    "electrochemistry::faraday_mass_deposited": ("B", lambda a: r_faraday_mass(*fl(a))),
    "electrochemistry::faraday_charge_required": ("B", lambda a: r_faraday_charge(*fl(a))),
    "electrochemistry::faraday_mass_ratio": ("B", lambda a: r_faraday_ratio(*fl(a))),
    "electrochemistry::charge_to_moles_electrons": ("B", lambda a: dv(fl(a)[0], FARADAY)),
    "electrochemistry::moles_electrons_to_charge": ("B", lambda a: fl(a)[0] * FARADAY),
    "kinetics::arrhenius_rate": ("C", lambda a: r_arrhenius(*fl(a))),
    "kinetics::half_life_first_order": ("B", lambda a: r_half_first(*fl(a))),
    "kinetics::first_order_concentration": ("C", lambda a: r_first(*fl(a))),
    "kinetics::second_order_concentration": ("B", lambda a: r_second(*fl(a))),
    "kinetics::zero_order_concentration": ("B", lambda a: r_zero(*fl(a))),
    "kinetics::zero_order_half_life": ("B", lambda a: r_zero_half(*fl(a))),
    "kinetics::nth_order_half_life": ("C", lambda a: r_nth(*fl(a))),
    "kinetics::steady_state_intermediate": ("B", lambda a: r_ssi(*fl(a))),
    "kinetics::steady_state_rate": ("B", lambda a: fl(a)[0] * fl(a)[1]),
    "kinetics::pre_equilibrium_rate": ("B", lambda a: r_pre(*fl(a))),
    "kinetics::michaelis_menten": ("B", lambda a: r_mm(*fl(a))),
    "kinetics::lineweaver_burk": ("B", lambda a: r_lb(*fl(a))),
    "kinetics::collision_theory_rate": ("C", lambda a: r_collision(*fl(a))),
    "kinetics::eyring_rate": ("C", lambda a: r_eyring(*fl(a))),
    "kinetics::eyring_rate_from_activation": ("C", lambda a: r_eyring_act(*fl(a))),
}
TABLE_RULES = {}
for name, v in PRESETS.items():
    TABLE_RULES[f"gas::{name}"] = [("f", x, f2b(x)) for x in v]
TABLE_RULES["phase::ANTOINE_DATA"] = [("v", len(ANTOINE), None)]
TABLE_RULES["thermochem::THERMOCHEM_DATA"] = [("v", len(TC), None)]
TABLE_RULES["thermochem::SHOMATE_DATA"] = [("v", len(SH), None)]
TABLE_RULES["electrochemistry::STANDARD_POTENTIALS"] = [("v", len(SP), None)]
def tc_toks(row):
    return [("s", row[0], None), ("s", row[1], None)] + [("f", x, f2b(x)) for x in row[2]]
def sh_toks(c):
    return [("f", x, f2b(x)) for x in c]

def antoine_toks(row):
    return [("f", x, f2b(x)) for x in row]

def nudge_s(fn, rule, ins, want, a, md, log):
    """S over one row: each call in md nudged ±1 ulp in turn; (largest output change in ulp, branch flips, where)."""
    S, flips, worst = 0, 0, None
    for k in md:
        for d in (1, -1):
            CALLS[0] = 0
            NUDGE[0] = (k, d)
            got = out_toks(rule(ins))
            if not toks_shape(got, want):
                flips += 1
                continue
            for (gk, gv, gb), (wk, wv, wb) in zip(got, want):
                if gk == "f":
                    dist = ulp_dist(gb, wb)
                    if dist > S:
                        S = dist
                        worst = f"{fn}{tuple(a)}: {log[k][0]} call {k} nudged {'+' if d > 0 else '-'}1 ulp"
    NUDGE[0] = None
    return S, flips, worst

root_paths = {}   # weak acid / base Newton exits, recomputed
root_iters = [0, 0]  # the most iterations entered, the most updates

# ADR 0006: the corrected oracle of one row, and the tallies of glibc's misroundings over the golden calls.
PRIM_CALLS = {}      # kind -> [calls, misrounded calls, {arguments}, {misrounded arguments}, max ulp]
CORR_FLIPS = {}      # fn -> corrected rows whose outputs change shape (Ok <-> Err, an Err variant, the root count)
CORR_ABSORBED = {}   # fn -> rows with a misrounded call whose outputs don't change
def misrounded(log):
    """The calls of a row's log, (kind, args, glibc's result), whose result is not correctly rounded: each distinct
    call once, in call order, as (MPFR name, args, glibc's result, the correctly rounded result)."""
    out, seen = [], set()
    for kind, args, y in log:
        name = KIND_MPFR[kind]
        cr = cr_of(name, args)
        key = (name,) + tuple(f2b(x) for x in args)
        st = PRIM_CALLS.setdefault(kind, [0, 0, set(), set(), 0])
        st[0] += 1
        st[2].add(key)
        if same_f(y, cr):
            continue
        st[1] += 1
        st[3].add(key)
        st[4] = max(st[4], ulp_dist(f2b(y), f2b(cr)))
        if key not in seen:
            seen.add(key)
            out.append((name, args, y, cr))
    return out
def corrected_row(mod, fn, rule, ins, outs, where, log):
    """Evaluates the row's replica again with every primitive result correctly rounded (MPFR), down whatever branch
    that takes, and records a correction where the outputs differ from the oracle's."""
    mis = misrounded(log)
    CALLS[0] = 0
    CALL_LOG.clear()
    CORRECT[0] = True
    try:
        ctoks = out_toks(rule(ins))
    finally:
        CORRECT[0] = False
    if toks_equal(ctoks, outs):
        if mis: CORR_ABSORBED[fn] = CORR_ABSORBED.get(fn, 0) + 1
        return
    check(mis or not LIBM_EXACT, f"{where}: the corrected outputs differ, but glibc rounded every call correctly")
    if not toks_shape(ctoks, outs): CORR_FLIPS[fn] = CORR_FLIPS.get(fn, 0) + 1
    add_correction(mod, where, ctoks, mis)

m2a_counts = {}
shape_only = 0
sens = {}       # fn -> [S, rows nudged, branch flips, worst row]
sens_rows = {}  # fn -> [(row S, branch flips, a nudged result is subnormal, inputs)] (M3b: the rows past K = 64)
flame_paths = {}  # "newton" / "bisection" / "err" -> rows
for mod in ("reaction", "thermo", "phase", "gas", "thermochem", "solution", "electrochemistry", "kinetics"):
    for fn, cls, ins, outs, where in R.get(mod, []):
        m2a_counts[fn] = m2a_counts.get(fn, 0) + 1
        if fn in TABLE_RULES:
            check(cls == "A" and ins == [] and toks_equal(outs, TABLE_RULES[fn]), f"{where}: {fn}")
            continue
        if fn == "thermochem::THERMOCHEM_DATA[]":
            k = ins[0][1]
            want = tc_toks(TC[k]) if 0 <= k < len(TC) else []
            check(cls == "A" and toks_equal(outs, want), f"{where}: THERMOCHEM_DATA[{k}]")
            continue
        if fn == "thermochem::SHOMATE_DATA[]":
            k = ins[0][1]
            want = [("s", SH[k][0], None)] + sh_toks(SH[k][1]) if 0 <= k < len(SH) else []
            check(cls == "A" and toks_equal(outs, want), f"{where}: SHOMATE_DATA[{k}]")
            continue
        if fn in ("thermochem::lookup_thermochem", "thermochem::lookup_shomate", "thermochem::lookup_shomate_at"):
            name_ = ins[0][1]
            if fn == "thermochem::lookup_thermochem":
                hit = tc_lookup(name_)
                want = [("n", "", None)] if hit is None else tc_toks(hit)
            else:
                hit = sh_lookup(name_) if fn == "thermochem::lookup_shomate" else sh_lookup_at(name_, ins[1][1])
                want = [("n", "", None)] if hit is None else sh_toks(hit)
            check(cls == "A" and toks_equal(outs, want), f"{where}: {fn}({name_!r})")
            continue
        if fn == "phase::ANTOINE_DATA[]":
            k = ins[0][1]
            want = [("s", ANTOINE[k][0], None)] + antoine_toks(ANTOINE[k][1]) if 0 <= k < len(ANTOINE) else []
            check(cls == "A" and toks_equal(outs, want), f"{where}: ANTOINE_DATA[{k}]")
            continue
        if fn == "electrochemistry::STANDARD_POTENTIALS[]":
            k = ins[0][1]
            want = sp_toks(SP[k]) if 0 <= k < len(SP) else []
            check(cls == "A" and toks_equal(outs, want), f"{where}: STANDARD_POTENTIALS[{k}]")
            continue
        if fn == "electrochemistry::lookup_half_reaction":
            hit = sp_lookup(ins[0][1])            # exact, case-sensitive, first in table order
            want = [("n", "", None)] if hit is None else sp_toks(hit)
            check(cls == "A" and toks_equal(outs, want), f"{where}: lookup_half_reaction({ins[0][1]!r})")
            continue
        if fn == "phase::lookup_antoine":
            hit = next((c for n, c in ANTOINE if n == ins[0][1]), None)
            want = [("n", "", None)] if hit is None else antoine_toks(hit)
            check(cls == "A" and toks_equal(outs, want), f"{where}: lookup_antoine({ins[0][1]!r})")
            continue
        if fn not in RULES:
            check(False, f"{where}: {mod}.tsv function {fn} has no rule here")
            continue
        want_cls, rule = RULES[fn]
        a = fl(ins)
        unit_row = fn.startswith("phase::antoine_pressure") and a == [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        check(cls == ("A" if unit_row else want_cls), f"{where}: {fn} class {cls}")
        CALLS[0] = 0
        NUDGE[0] = None
        CALL_LOG.clear()
        base = rule(ins)
        ncalls = CALLS[0]
        log = list(CALL_LOG)
        want = out_toks(base)
        if fn == "thermochem::adiabatic_flame_temperature":
            path = FLAME_PATH[0]
            note = NOTES.get(where, "")
            noted = note.rsplit(" [", 1)[1][:-1] if note.endswith("]") else None
            check(noted == path, f"{where}: solver path noted {noted!r}, recomputed {path!r}")
            kind = "err" if path is None else path.split(" ", 1)[0].rstrip(",")
            flame_paths[kind] = flame_paths.get(kind, 0) + 1
        if fn in ("thermochem::reaction_enthalpy", "thermochem::reaction_gibbs_energy"):
            check(NOTES.get(where, "").endswith("[J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002 #1)]"),
                  f"{where}: {fn} row note does not say ADR 0002 #1's × 1000")
        if cls == "C" and not LIBM_EXACT:
            check(toks_shape(outs, want), f"{where}: {fn}{tuple(a)} -> {outs}, recomputed {base!r}")
            shape_only += 1
        else:
            check(toks_equal(outs, want), f"{where}: {fn}{tuple(a)} -> {[t[1] for t in outs]}, recomputed {base!r}")
        if fn in ("solution::weak_acid_ph", "solution::weak_base_ph"):
            path = ROOT_PATH[0]
            note = NOTES.get(where, "")
            noted = note.rsplit(" [", 1)[1][:-1] if note.endswith("]") else None
            check(noted == path, f"{where}: Newton path noted {noted!r}, recomputed {path!r}")
            if path is not None:
                root_paths[path.split(",")[0]] = root_paths.get(path.split(",")[0], 0) + 1
                ups = int(re.search(r"(\d+) (?:updates|iterations)", path).group(1))
                root_iters[0] = max(root_iters[0], ups if "relative step" in path else ups + 1)
                root_iters[1] = max(root_iters[1], ups)
        if log:
            corrected_row(mod, fn, rule, ins, outs, where, log)
        if SENSITIVITY and cls == "C" and ncalls > 0:
            st = sens.setdefault(fn, [0, 0, 0, None, 0, 0])
            st[1] += 1
            md = [k for k, (kind, args, y) in enumerate(log) if nudged(kind, args, y)]
            st[4] = max(st[4], len(md))
            st[5] += sum(1 for kind, args, y in log if glibc_not_cr(kind, args, y))
            S, flips, worst = nudge_s(fn, rule, ins, want, a, md, log)
            sub = any(0 < abs(log[k][2]) < 2.0 ** -1022 for k in md)
            sens_rows.setdefault(fn, []).append((S, flips, sub, tuple(a)))
            st[2] += flips
            if S > st[0]: st[0], st[3] = S, worst
known = set(RULES) | set(TABLE_RULES) | {"phase::ANTOINE_DATA[]", "phase::lookup_antoine",
    "thermochem::THERMOCHEM_DATA[]", "thermochem::SHOMATE_DATA[]", "thermochem::lookup_thermochem",
    "thermochem::lookup_shomate", "thermochem::lookup_shomate_at", "electrochemistry::STANDARD_POTENTIALS[]",
    "electrochemistry::lookup_half_reaction"}
for fn in sorted(known):
    check(fn in m2a_counts, f"M2/M3: no golden rows for {fn}")
print("M2/M3 rows checked:", sum(m2a_counts.values()), "over", len(m2a_counts), "functions"
      + (f" ({shape_only} class-C rows by shape only)" if shape_only else ""))
print(f"thermochem: adiabatic_flame_temperature solver paths (recomputed, and as noted): "
      f"{flame_paths.get('newton', 0)} newton, {flame_paths.get('bisection', 0)} bisection, "
      f"{flame_paths.get('err', 0)} Err")
check(flame_paths.get("bisection", 0) >= 100, "thermochem: fewer than 100 bisection-path flame rows")

# Rust's typed flame constants (FLAME_NEWTON_TOL_J = 1e-3, FLAME_BISECTION_TOL = 1e-6, the 40 J/(mol·K)
# guess) are inline in src/thermochem.cyr: each must change some golden row's bits when perturbed, so
# the corpus pins it. (The 1e-12 derivative threshold can't bind on a bracketed input: tests/thermochem.tcyr
# pins it on the private Newton helper.)
flame_rows = [(ins, outs) for fn, cls, ins, outs, where in R.get("thermochem", [])
              if fn == "thermochem::adiabatic_flame_temperature" and outs[0][0] == "f"]
for label, kw in (("FLAME_NEWTON_TOL_J × 0.5", {"ntol": 5e-4}), ("FLAME_NEWTON_TOL_J × 2", {"ntol": 2e-3}),
                  ("FLAME_BISECTION_TOL × 0.5", {"btol": 5e-7}), ("FLAME_BISECTION_TOL × 2", {"btol": 2e-6}),
                  ("guess Cp 39", {"gcp": 39.0}), ("guess Cp 41", {"gcp": 41.0})):
    changed = 0
    for ins, outs in flame_rows:
        prods, k = species_at(ins, 1)
        got = r_flame(ins[0][1], prods, ins[k][1], **kw)
        if isinstance(got, Err) or f2b(got) != outs[0][2]: changed += 1
    print(f"  pinning: with {label}, {changed} of {len(flame_rows)} Ok flame rows change")
    check(changed > 0, f"thermochem: no golden flame row pins {label}")

# The weak acid / base Newton (M3a): its exits as recomputed (and as noted), and Rust's private constants, inline in
# src/solution.cyr: CHARGE_BALANCE_REL_TOL (1e-12) must change some golden row's bits when perturbed;
# CHARGE_BALANCE_MAX_ITER (50) is unobservable when every row ends in fewer iterations (any cap at or above the
# most iterations entered gives the same corpus), which is reported, not failed.
print(f"solution: weak acid / base Newton exits (recomputed, and as noted): "
      + ", ".join(f"{root_paths[k]} {k[len('root: '):]}" for k in sorted(root_paths))
      + f"; at most {root_iters[1]} updates, {root_iters[0]} iterations entered")
weak_rows = [(fn == "solution::weak_base_ph", ins, outs) for fn, cls, ins, outs, where in R.get("solution", [])
             if fn in ("solution::weak_acid_ph", "solution::weak_base_ph") and outs[0][0] == "f"]
for label, kw in (("CHARGE_BALANCE_REL_TOL × 0.5", {"tol": 5e-13}), ("CHARGE_BALANCE_REL_TOL × 2", {"tol": 2e-12}),
                  (f"CHARGE_BALANCE_MAX_ITER = {root_iters[0]}", {"max_iter": root_iters[0]}),
                  (f"CHARGE_BALANCE_MAX_ITER = {root_iters[0] - 1}", {"max_iter": root_iters[0] - 1})):
    changed = 0
    for base, ins, outs in weak_rows:
        got = r_weak(ins[0][1], ins[1][1], base, **kw)
        if isinstance(got, Err) or f2b(got) != outs[0][2]: changed += 1
    print(f"  pinning: with {label}, {changed} of {len(weak_rows)} Ok weak acid / base rows change")
    if "REL_TOL" in label:
        check(changed > 0, f"solution: no golden weak acid / base row pins {label}")
    elif label.endswith(f"= {root_iters[0]}"):
        check(changed == 0, f"solution: {label} changes rows: the cap is observable")
print(f"  CHARGE_BALANCE_MAX_ITER = 50 is not observable on the corpus (no row enters more than {root_iters[0]} "
      f"iterations); tests/solution.tcyr covers it on the private helper as far as an input can reach")

# kinetics (M3b): Rust's std::f64::consts::LN_2 and private FIRST_ORDER_TOLERANCE (1e-12) are typed inline in
# src/kinetics.cyr: each must change some golden row (its bits or Ok/Err) when perturbed, so the corpus pins it.
kin_rows = {fn: [(ins, outs) for f_, cls, ins, outs, where in R.get("kinetics", []) if f_ == fn]
            for fn in ("kinetics::half_life_first_order", "kinetics::nth_order_half_life")}
for label, fn, rule in (
        ("LN_2 + 1 ulp", "kinetics::half_life_first_order", lambda a: r_half_first(*a, ln2=nextup(LN2))),
        ("LN_2 − 1 ulp", "kinetics::half_life_first_order", lambda a: r_half_first(*a, ln2=nextdown(LN2))),
        ("LN_2 + 1 ulp", "kinetics::nth_order_half_life", lambda a: r_nth(*a, ln2=nextup(LN2))),
        ("LN_2 − 1 ulp", "kinetics::nth_order_half_life", lambda a: r_nth(*a, ln2=nextdown(LN2))),
        ("FIRST_ORDER_TOLERANCE × 0.5", "kinetics::nth_order_half_life", lambda a: r_nth(*a, tol=5e-13)),
        ("FIRST_ORDER_TOLERANCE × 2", "kinetics::nth_order_half_life", lambda a: r_nth(*a, tol=2e-12))):
    rows_ = kin_rows[fn]
    changed = sum(1 for ins, outs in rows_ if not toks_equal(out_toks(rule(fl(ins))), outs))
    print(f"  pinning: with {label}, {changed} of {len(rows_)} {fn.split('::')[1]} rows change")
    check(changed > 0, f"kinetics: no golden {fn} row pins {label}")
for label, tol in (("+ 1 ulp", nextup(1e-12)), ("− 1 ulp", nextdown(1e-12))):
    changed = sum(1 for ins, outs in kin_rows["kinetics::nth_order_half_life"]
                  if not toks_equal(out_toks(r_nth(*fl(ins), tol=tol)), outs))
    print(f"  FIRST_ORDER_TOLERANCE {label}: {changed} rows change (|n − 1| is a multiple of 2^-52 or 2^-53 near "
          f"1e-12, so only a change across a multiple can show)")

# The generated data files against the Rust literals.
for rel in ("src/phase_data.cyr", "src/gas_data.cyr"):
    text = (root / rel).read_text()
    n = 0
    for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = (.+);$", text, re.M):
        path, name, lit = m.groups()
        check(name == "KM_" + path.split("::")[1], f"{rel}: {name} names kimiya::{path}")
        check(re.fullmatch(r"-?[0-9]+\.[0-9]+", lit), f"{rel}: {name} = {lit} is not a positional decimal")
        check(path in golden_consts and f2b(float(lit)) == golden_consts.get(path), f"{rel}: {name} = {lit} is not the golden {path}")
        n += 1
    lengths = re.findall(r"^const (KM_\w+_COUNT) = (\d+);", text, re.M)   # a table's length, not a constant
    check(n + len(lengths) == len(re.findall(r"^const KM_", text, re.M)), f"{rel}: a const without its `# Rust:` line")
    for name, val in lengths:
        check(name == "KM_ANTOINE_COUNT" and int(val) == len(ANTOINE), f"{rel}: {name} = {val}")
    print(f"{rel} constant literals checked: {n}")
gtext = (root / "src/gas_data.cyr").read_text()
gen = {name: [float(x) for x in vals] for name, *vals in re.findall(
    r"^var KM_(PR_\w+) = KmPengRobinsonParams \{ tc: ([0-9.]+), pc: ([0-9.]+), omega: ([0-9.]+) \};$", gtext, re.M)}
check(gen == PRESETS, f"gas_data.cyr: PR presets {gen} are not gas.rs's {PRESETS}")
ptext = (root / "src/phase_data.cyr").read_text()
check(re.search(r"^const _KM_PA_PER_MMHG = ([0-9.]+);$", ptext, re.M) and
      float(re.search(r"^const _KM_PA_PER_MMHG = ([0-9.]+);$", ptext, re.M).group(1)) == PA_PER_MMHG,
      "phase_data.cyr: _KM_PA_PER_MMHG is not phase.rs's PA_PER_MMHG")
cols = {}
for col in ("A", "B", "C", "T_MIN_C", "T_MAX_C"):
    body = re.search(rf"^var KM_ANTOINE_{col}: f64\[8\] = \{{(.*?)\}};", ptext, re.M | re.S).group(1)
    cols[col] = [float(x) for x in re.findall(r"^\s+(-?[0-9.]+),", body, re.M)]
gen_rows = [[cols[c][k] for c in ("A", "B", "C", "T_MIN_C", "T_MAX_C")] for k in range(8)]
check(gen_rows == [c for _, c in ANTOINE], "phase_data.cyr: the Antoine columns are not ANTOINE_DATA")
print("phase_data.cyr / gas_data.cyr tables checked against phase.rs / gas.rs")
ttext = (root / "src/thermochem_data.cyr").read_text()
def cyr_col(name, kind=float):
    m = re.search(rf"^var {name}: \w+\[(\d+)\] = \{{(.*?)\}};", ttext, re.M | re.S)
    check(m is not None, f"thermochem_data.cyr: no {name}")
    if m is None: return []
    vals = [kind(x) for line in m.group(2).splitlines()
            for x in re.findall(r"(-?[0-9.]+),", line.split("#", 1)[0])]
    check(len(vals) == int(m.group(1)), f"thermochem_data.cyr: {name} declares {m.group(1)}, holds {len(vals)}")
    return vals
def cyr_strings(prefix):
    raw = bytes(cyr_col(f"_{prefix}_BYTES", int))
    offs = cyr_col(f"_{prefix}_OFF", int)
    return [raw[o:raw.index(0, o)].decode("ascii") for o in offs], raw.count(0) == len(offs)
check(re.search(r"^const KM_THERMOCHEM_COUNT = 76;$", ttext, re.M), "thermochem_data.cyr: KM_THERMOCHEM_COUNT")
check(re.search(r"^const KM_SHOMATE_COUNT = 34;$", ttext, re.M), "thermochem_data.cyr: KM_SHOMATE_COUNT")
for prefix, want in (("KM_THERMOCHEM_FORMULA", [r[0] for r in TC]), ("KM_THERMOCHEM_NAME", [r[1] for r in TC]),
                     ("KM_SHOMATE_FORMULA", [f_ for f_, _ in SH])):
    got_s, one_each = cyr_strings(prefix)
    check(got_s == want and one_each, f"thermochem_data.cyr: the {prefix} blob is not thermochem.rs's strings")
for k, col in enumerate(("DELTA_HF_KJ", "DELTA_GF_KJ", "S_STANDARD_J")):
    check(cyr_col(f"KM_THERMOCHEM_{col}") == [r[2][k] for r in TC], f"thermochem_data.cyr: KM_THERMOCHEM_{col}")
for k, fld in enumerate(SH_FIELDS):
    check(cyr_col(f"KM_SHOMATE_{fld.upper()}") == [c[k] for _, c in SH], f"thermochem_data.cyr: KM_SHOMATE_{fld.upper()}")
sm = re.search(r"^struct KmShomate \{(.*?)\n\}", ttext, re.M | re.S)
check(sm is not None and re.findall(r"^\s+(\w+): f64;", sm.group(1), re.M) == list(SH_FIELDS),
      "thermochem_data.cyr: struct KmShomate is not Shomate's fields in order")
print("thermochem_data.cyr tables checked against thermochem.rs (76 + 34 rows, 3 string blobs, 10 f64 columns)")

# src/electrochemistry_data.cyr against the Rust literal (M3a).
etext = (root / "src/electrochemistry_data.cyr").read_text()
def cyr_col_in(text, rel, name, kind=float):
    m = re.search(rf"^var {name}: \w+\[(\d+)\] = \{{(.*?)\}};", text, re.M | re.S)
    check(m is not None, f"{rel}: no {name}")
    if m is None: return []
    vals = [kind(x) for line in m.group(2).splitlines()
            for x in re.findall(r"(-?[0-9.]+),", line.split("#", 1)[0])]
    check(len(vals) == int(m.group(1)), f"{rel}: {name} declares {m.group(1)}, holds {len(vals)}")
    return vals
def cyr_strings_in(text, rel, prefix, enc):
    raw = bytes(cyr_col_in(text, rel, f"_{prefix}_BYTES", int))
    offs = cyr_col_in(text, rel, f"_{prefix}_OFF", int)
    return [raw[o:raw.index(0, o)].decode(enc) for o in offs], raw.count(0) == len(offs)
rel_e = "src/electrochemistry_data.cyr"
check(re.search(r"^const KM_STANDARD_POTENTIALS_COUNT = 49;$", etext, re.M), f"{rel_e}: KM_STANDARD_POTENTIALS_COUNT")
for prefix, enc, k in (("KM_STANDARD_POTENTIALS_EQUATION", "utf-8", 0), ("KM_STANDARD_POTENTIALS_COUPLE", "ascii", 1)):
    got_s, one_each = cyr_strings_in(etext, rel_e, prefix, enc)
    check(got_s == [r[k] for r in SP] and one_each, f"{rel_e}: the {prefix} blob is not electrochemistry.rs's strings")
check(cyr_col_in(etext, rel_e, "KM_STANDARD_POTENTIALS_ELECTRONS", int) == [r[2] for r in SP], f"{rel_e}: electrons")
check(cyr_col_in(etext, rel_e, "KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL") == [r[3] for r in SP], f"{rel_e}: E°")
check(len(re.findall(r"CRC Handbook, Electrochemical Series", etext)) == 2 and
      len(re.findall(r"CRC Handbook, Electrochemical Series", _ec_body)) == 2,
      f"{rel_e}: the literal's two CRC Handbook notes are not carried over")
print(f"{rel_e} checked against electrochemistry.rs (49 rows: a UTF-8 equation blob, a couple blob, electrons, E°)")

# ── 12. ADR 0006: the corrected oracle, tests/golden/<module>.cr.tsv ─────────────────────
def tok_text(t):
    k, v, b = t
    if k == "f": return f"f:{b:016x}"
    if k in ("i", "v"): return f"{k}:{v}"
    if k == "E": return f"E:{v}"
    raise ValueError(f"no corrected-output token for {t!r}")
def toks_text(ts): return " ".join(tok_text(t) for t in ts) if ts else "-"
def fbits(x): return f"f:{f2b(x):016x}"
def call_text(c):
    name, args, g, cr = c
    return f"{name}({', '.join(fbits(x) for x in args)}) glibc {fbits(g)} correctly rounded {fbits(cr)}"
def cr_header(mod):
    """The five header lines of <mod>.cr.tsv, in the golden files' style; the provenance is the golden file's."""
    prov = (gold / f"{mod}.tsv").read_text().split("\n", 2)[1].split(". Generated by", 1)[0]
    rows = len(CORRECTIONS.get(mod, []))
    return [f"# kimiya corrected oracle (ADR 0006), format 1: {mod}",
            f"{prov}. Correctly rounded: MPFR {MPFR.VERSION}. Generated by scripts/oracle/validate.py; do not edit.",
            "# regenerate: python3 -I scripts/oracle/validate.py . --write-corrections",
            "# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: "
            "the oracle's operations with every primitive result correctly rounded; note: the oracle's outputs, then "
            "each call glibc misrounded)",
            f"# rows: {rows}"]
def cr_rows(mod):
    rows = []
    for where, ctoks, mis in CORRECTIONS.get(mod, []):
        fn_t, ins_t, outs_t = RAW[where]
        rows.append((fn_t, ins_t, toks_text(ctoks), "; ".join([f"oracle {outs_t}"] + [call_text(c) for c in mis])))
    return rows
def cr_text(mod):
    return "\n".join(cr_header(mod) + ["\t".join(r) for r in cr_rows(mod)]) + "\n"

# The serializer writes the generator's tokens: every golden output field that gets a correction is rewritten
# from its parsed tokens, byte for byte.
for mod in CORRECTIONS:
    by_where = {w: o for f_, c_, i_, o, w in R[mod]}
    for where, ctoks, mis in CORRECTIONS[mod]:
        check(toks_text(by_where[where]) == RAW[where][2], f"{where}: the outputs don't re-serialize as written")

GLIBC_OF = {"exp": math.exp, "expm1": math.expm1, "log": math.log, "log10": math.log10, "pow": math.pow,
            "cos": math.cos, "sin": math.sin, "acos": math.acos, "atan2": math.atan2, "exp2": math.exp2}
CR_CALL = re.compile(r"(\w+)\((f:[0-9a-f]{16}(?:, f:[0-9a-f]{16})?)\) glibc f:([0-9a-f]{16}) correctly rounded "
                     r"f:([0-9a-f]{16})")
def check_listed_call(rel, k, text):
    """One call a note lists, against MPFR directly (not the replica): its correctly rounded value is MPFR's, glibc's
    differs from it, and (on the recorded platform) glibc's is this machine's libm result."""
    m = CR_CALL.fullmatch(text)
    check(m is not None and m.group(1) in GLIBC_OF, f"{rel}:{k}: malformed call {text!r}")
    if m is None or m.group(1) not in GLIBC_OF: return
    name = m.group(1)
    args = tuple(b2f(int(x[2:], 16)) for x in m.group(2).split(", "))
    g, cr = b2f(int(m.group(3), 16)), b2f(int(m.group(4), 16))
    want = MPFR.cr(name, args[0]) if len(args) == 1 else MPFR.cr2(name, args[0], args[1])
    check(same_f(cr, want), f"{rel}:{k}: {name}{args} correctly rounded is {want!r}, the note says {cr!r}")
    check(not same_f(g, cr), f"{rel}:{k}: {name}{args}: glibc's result is correctly rounded")
    if LIBM_EXACT:
        try: lib = GLIBC_OF[name](*args)
        except (OverflowError, ValueError): lib = None
        check(lib is not None and same_f(g, lib), f"{rel}:{k}: {name}{args}: glibc gives {lib!r}, the note says {g!r}")

cr_on_disk = sorted(p.name[:-len(CR_SUFFIX)] for p in gold.glob("*" + CR_SUFFIX))
if WRITE_CORRECTIONS:
    check(LIBM_EXACT, f"--write-corrections: this machine ({running_platform()}) is not the corpus's libm platform "
                      f"({recorded_platform('gas')}); the notes hold glibc's results, so nothing is written")
    if LIBM_EXACT:
        for mod in tsvs:
            path = gold / f"{mod}{CR_SUFFIX}"
            if CORRECTIONS.get(mod):
                path.write_text(cr_text(mod))
                print(f"wrote tests/golden/{mod}{CR_SUFFIX}: {len(CORRECTIONS[mod])} rows")
            elif path.exists():
                path.unlink()
                print(f"removed tests/golden/{mod}{CR_SUFFIX}: no row to correct")
        cr_on_disk = sorted(p.name[:-len(CR_SUFFIX)] for p in gold.glob("*" + CR_SUFFIX))
for mod in cr_on_disk:
    check(mod in tsvs, f"tests/golden/{mod}{CR_SUFFIX}: no golden module {mod}")
for mod in sorted(set(CORRECTIONS) | set(m for m in cr_on_disk if m in tsvs)):
    rel = f"tests/golden/{mod}{CR_SUFFIX}"
    want_rows = cr_rows(mod)
    path = gold / f"{mod}{CR_SUFFIX}"
    if not path.exists():
        check(False, f"{rel}: missing; {len(want_rows)} rows need a correction (run with --write-corrections)")
        continue
    check(want_rows, f"{rel}: stale: no {mod} row needs a correction")
    text = path.read_text()
    lines = text.split("\n")
    check(lines[-1] == "", f"{rel}: no final newline")
    lines = lines[:-1]
    head, want_head = lines[:5], cr_header(mod)
    for k in range(5):
        if k == 1 and not LIBM_EXACT: continue
        check(k < len(head) and head[k] == want_head[k], f"{rel}: header line {k + 1} is "
              f"{head[k] if k < len(head) else None!r}, want {want_head[k]!r}")
    got_rows = []
    for k, line in enumerate(lines[5:], 6):
        if line.startswith("#"): continue
        cols = tuple(line.split("\t"))
        check(len(cols) == 4, f"{rel}:{k}: {len(cols)} fields")
        if len(cols) == 4: got_rows.append((k, cols))
    golden = {}
    for f_, c_, i_, o_, w in R[mod]:
        golden.setdefault(RAW[w][:2], RAW[w][2])
    want_by = {}
    for r in want_rows: want_by.setdefault(r[:2], []).append(r)
    got_by = {}
    for k, r in got_rows: got_by.setdefault(r[:2], []).append((k, r))
    for key, rs in got_by.items():
        for k, r in rs[len(want_by.get(key, [])):]:
            why = "no golden row has this fn and inputs" if key not in golden else "the oracle's outputs need no correction"
            check(False, f"{rel}:{k}: extra row ({why}): {r[0]} {r[1]}")
    for key, rs in want_by.items():
        have = got_by.get(key, [])
        for r in rs[len(have):]:
            check(False, f"{rel}: missing row: {r[0]} {r[1]} -> {r[2]}")
        for (k, g), w in zip(have, rs):
            check(g[2] == w[2], f"{rel}:{k}: stale outputs {g[2]!r}, the corrected oracle gives {w[2]!r}")
            if LIBM_EXACT:
                check(g[3] == w[3], f"{rel}:{k}: stale note {g[3]!r}, want {w[3]!r}")
    check([r[:2] for k, r in got_rows] == [r[:2] for r in want_rows] or len(got_rows) != len(want_rows),
          f"{rel}: rows are not in golden-file order")
    for k, r in got_rows:
        parts = r[3].split("; ")
        check(parts[0] == f"oracle {golden.get(r[:2])}", f"{rel}:{k}: the note's oracle outputs {parts[0]!r} are not the golden row's")
        check(len(parts) > 1, f"{rel}:{k}: the note lists no misrounded call")
        for call in parts[1:]:
            check_listed_call(rel, k, call)
    if LIBM_EXACT:
        check(text == cr_text(mod), f"{rel}: differs from the corrected oracle (regenerate with --write-corrections)")
    print(f"{rel}: {len(got_rows)} rows checked against the corrected oracle (recomputed) and MPFR (each listed call)")

# The report: corrected rows per module and function, the calls glibc misrounds, glibc's cbrt.
print(f"ADR 0006 corrected oracle (MPFR {MPFR.VERSION}): rows whose outputs change when every primitive result is "
      f"correctly rounded")
for mod in sorted(CORRECTIONS):
    per = {}
    for where, ctoks, mis in CORRECTIONS[mod]:
        per[RAW[where][0]] = per.get(RAW[where][0], 0) + 1
    print(f"  {mod}: {len(CORRECTIONS[mod])} rows")
    for fn in sorted(per):
        extra = f" ({CORR_FLIPS[fn]} change shape: Ok <-> Err, the Err, or the root count)" if CORR_FLIPS.get(fn) else ""
        print(f"    {fn:44} {per[fn]}{extra}")
if CORR_ABSORBED:
    print("  rows with a misrounded call whose outputs don't change: "
          + ", ".join(f"{fn} {n}" for fn, n in sorted(CORR_ABSORBED.items())))
print("glibc vs correctly rounded over the M2–M3b replicas' golden calls (the math rows are above):")
for kind in sorted(PRIM_CALLS):
    c, m_, args, margs, mx = PRIM_CALLS[kind]
    print(f"  {kind:6} {c} calls ({len(args)} distinct arguments): {m_} misrounded ({len(margs)} distinct), max {mx} ulp")
cb = {f2b(x): x for x in CBRT_ARGS}
cb_or = sum(1 for x in cb.values() if not same_f(cbrt_cr(x), cr_of("cbrt", (x,))))
cb_gl = [x for x in cb.values() if not same_f(math.cbrt(x), cr_of("cbrt", (x,)))]
check(cb_or == 0, f"cbrt: the oracle's (Rust's) cbrt is not MPFR's on {cb_or} golden arguments")
print(f"  cbrt: the oracle's (Rust's compiler-builtins) is correctly rounded on all {len(cb)} distinct golden arguments "
      f"(math rows and the Peng–Robinson cubic); glibc's own cbrt (not the oracle's) misrounds {len(cb_gl)} of them"
      + (f", max {max(ulp_dist(f2b(math.cbrt(x)), f2b(cr_of('cbrt', (x,)))) for x in cb_gl)} ulp" if cb_gl else ""))

if SENSITIVITY:
    print("sensitivity, a diagnostic of conditioning (ADR 0006): S = the largest output change, in ulp, from a 1-ulp")
    print("nudge of one primitive result; c = the most nudged calls in one row; K = 2·S·c + 2 is the retired rule's")
    print("bound, for comparison (the gate is 0 ulp against the corrected oracle, port-parity.md §10)")
    for fn in sorted(sens):
        S, rows, flips, worst, calls, crcalls = sens[fn]
        print(f"  {fn:40} S = {S:<3} c = {calls}  K = {2 * S * calls + 2:<4} ({rows} rows; {crcalls} calls where "
              f"glibc is not correctly rounded; {flips} branch flips) worst: {worst}")
    print("ill-conditioned on their grids (K past 64; gated exactly since M3c-3): the rows whose own S puts K past 64")
    print("(or flips a branch), how many of them nudge a subnormal primitive result, and S over the other rows:")
    for fn in sorted(sens):
        S, rows, flips, worst, calls, crcalls = sens[fn]
        if 2 * S * calls + 2 <= 64 and flips == 0:
            continue
        over = [r for r in sens_rows[fn] if 2 * r[0] * calls + 2 > 64 or r[1] > 0]
        rest = max((r[0] for r in sens_rows[fn] if not (2 * r[0] * calls + 2 > 64 or r[1] > 0)), default=0)
        nsub = sum(1 for r in over if r[2])
        normal = [r for r in sens_rows[fn] if not r[2] and r[1] == 0]
        snorm = max((r[0] for r in normal), default=0)
        print(f"  {fn:40} {len(over)} rows past K = 64 ({nsub} with a subnormal exp / pow result, "
              f"{sum(1 for r in over if r[1])} with branch flips); S over the other {rows - len(over)} rows = {rest} "
              f"(K = {2 * rest * calls + 2}); S over the {len(normal)} rows with no subnormal primitive result and no "
              f"flip = {snorm} (K = {2 * snorm * calls + 2})")
        for r in sorted(over, key=lambda r: -r[0])[:3]:
            print(f"      e.g. {fn}{r[3]}: S = {r[0]}")

print("FAILURES:", fails)
sys.exit(1 if fails else 0)
