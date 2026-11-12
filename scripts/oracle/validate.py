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
13. organic.tsv (M4b), every row against rust-old/src/organic.rs: the enums BondType, Geometry and FunctionalGroup
   parsed from their declarations (which proves the generated enums complete and in order, ADR 0002 #6's
   KM_GEOMETRY_T_SHAPED included: the generator's lists can't be matched exhaustively, the enums being
   #[non_exhaustive]), each const fn's match arms or matches! set, predict_geometry's match parsed and replayed over
   its grid (the 0..15 square and every u8 against the edges), enthalpy_from_bonds recomputed as two .sum()s from
   -0.0. src/organic_data.cyr is checked against the source (the section comments and notes of the matches carried),
   and ADR 0002 #6's acceptor column, which has no oracle, against its stated rule recomputed from the patterns (IUPAC
   2011's electron-rich region: a lone-pair atom O, N, S or a halogen, or a π bond).
14. inorganic.tsv (M4b), every row recomputed from rust-old/src/inorganic.rs: Ligand and field_strength, the private
   occupancy and unpaired-electron literals, IONIC_RADII, the formulas in Rust's order with the guards in Rust's
   order; every d count and spin, every u8 moment, every table symbol × every i8 charge covered. Each typed formula
   constant (the oracle's CFSE -0.4 and 0.6, 4/9, the 2.0s, Born-Lande's 1.0, the strong-field threshold 8 from above,
   the exact CFSE's divisor 5) must change some golden row when perturbed. ADR 0002 #12: every Ok CFSE row is checked
   both as the oracle's value and against the exact (3·eg − 2·t2g)/5 (Fraction arithmetic): the oracle within 1e-9 of
   the same fifth, and the note naming the exact value on exactly the 13 rows where the two differ. src/inorganic_data.cyr is checked against the literals, its ten constants
   against constants.tsv. Every golden module must have a section here (VALIDATED).
15. nuclear.tsv (M4c), every row recomputed from rust-old/src/nuclear.rs (glibc's exp, expm1, log and exp2), the grids'
   coverage, ADR 0002 #3's notes, the typed constants' pinning, src/nuclear_data.cyr against the literal, and ADR 0002
   #13 / #14 / #15's notes proven exactly (since M4d #13's correctly rounded N₀/2ⁿ by MPFR's scaled_exp2 and a 60-digit
   2⁻ᶠ).
16. potential.tsv (M4d), every row recomputed from rust-old/src/potential.rs (glibc's exp, atan2 and exp2), ADR 0002
   #16's correctly rounded σ·2^(1/6) proven with exact sixth powers, Rust's (4·ε)·d proven RN(4·ε·d) where 4·ε is finite
   (#17 and #18 folded into the overflow audit's #32 since B3, proven with it in §21), and both private tolerances pinned
   to the ulp;
   since M4e-2, ADR 0002 #20: distance_3d, bond_angle and dihedral_angle evaluated through MPFR at an unbounded exponent
   (mpfr_ref.UX, mpfr_ref.potential_*) on every row with finite coordinates — each #20 note is that value, each unnoted
   row's expected result equals it, no row whose f64 chain stays in the normal range is noted — and their scale
   invariance through the same evaluator.
17. spectroscopy.tsv (M4e), every row recomputed from rust-old/src/spectroscopy.rs (glibc's pow(10, ·), log10, and the
   FFT twiddles' cos and sin, which the release oracle calls as sincos: the rows prove its bits are glibc's sin and cos),
   with hisab 1.4.0's num::fft transcribed and the isotope pattern's stable sort; ADR 0002 #19's correctly rounded
   1/(100·x) proven with Fraction on every wavenumber row; src/spectroscopy_data.cyr against the literal (the
   SpectralSeries enum and n_lower, the isotope table) and R_H against RYDBERG / (1 + the private mass ratio) and the
   API; the typed thresholds' and constants' pinning where a row can reach them (the rest, and the stable order of
   equal masses, which no natural molecule's scratch shows with unequal intensities, are pinned by tests/spectroscopy.tcyr
   on the private helpers); the grids' coverage; and, since M4e-3, ADR 0002 #21's and #22's exact rationals rounded
   once, proven with Fraction on every Bohr, transition, Rydberg, series and reduced-mass row, with the structure that
   decides their exact values and ties (R's odd significands factored, the transition energy's domain searched
   exhaustively for midpoints); since M4e-4, ADR 0002 #23: scripts/oracle/isotopes.tsv (the cited AME2020 masses and TICE
   2009 abundances, ³⁶S included) parsed independently and checked against Rust's literal (the same isotopes plus ³⁶S,
   the shared abundances equal, Rust's masses the AME2020 ones to five decimals, each element summing to exactly 1), the
   mean masses against CIAAW's intervals, src/spectroscopy_data.cyr against the file, and every isotope_abundances and
   isotope_pattern row's port result (r_pattern on the new table): noted exactly where it differs from the oracle's, the
   note's tokens it bit for bit; the isotope thresholds pinned on the port's results too.
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
import decimal, importlib.util, math, platform, random, re, struct, sys
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
                                                                     "kinetics", "nuclear", "potential", "spectroscopy"))
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
# _km_cr_exp10 (= a correctly rounded pow(10, y)), "pow25" _km_cr_pow(x, 2.5), "cos" _km_cr_cos and "acos" _km_cr_acos
# (since M4d; ADR 0005's double-double _km_cos and _km_acos until then), (M4c) "exp2", glibc's exp2 in the release
# oracle (LLVM's rewrite of 2.0_f64.powf(y)), _km_cr_pow(2, y) in the port: the same real 2^y, (M4d) "atan2"
# _km_cr_atan2 (potential's angles; its args are (y, x)), and (M4e) "sin" with "cos", the FFT's twiddles: the release
# oracle calls glibc's sincos (LLVM joins hisab's angle.cos() and angle.sin()), whose bits are glibc's sin and cos (the
# golden rows prove it), and the port _km_cr_sincos.
# With CORRECT[0] set (the corrected oracle, ADR 0006), every result is replaced by its correctly rounded
# value from MPFR before it is logged or used. cbrt is not logged: the oracle's (Rust's) is correctly rounded.
CALLS = [0]
NUDGE = [None]
CORRECT = [False]
CALL_LOG = []
# The MPFR function each kind rounds; pow10's arguments are (10.0, y) and pow25's (x, 2.5).
KIND_MPFR = {"exp": "exp", "expm1": "expm1", "ln": "log", "log10": "log10", "pow": "pow", "pow10": "pow",
             "pow25": "pow", "cos": "cos", "acos": "acos", "exp2": "exp2", "atan2": "atan2", "sin": "sin"}
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

class Pairs(list):
    """A Vec<(f64, f64)> (M4e: the isotope pattern and abundances): `v:n`, then each pair's two floats."""
class Ints(list):
    """A Vec<usize> (M4e: find_peaks' indices): `v:n`, then n `i:` tokens."""
def out_toks(x):
    """The golden tokens of a Rust result: Err, a root list (`v:n` + floats), an int, or a float."""
    if isinstance(x, Err): return [("E", x.name, None)]
    if isinstance(x, Pairs): return [("v", len(x), None)] + [("f", z, f2b(z)) for pr in x for z in pr]
    if isinstance(x, Ints): return [("v", len(x), None)] + [("i", z, None) for z in x]
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
FLAME_ENDS = [None]              # (f(T₀), f(T_max)) of the last flame evaluated past validation (ADR 0002 #38)
def r_flame(dh, prods, t0, ntol=1e-3, btol=1e-6, gcp=40.0, eps=1e-12, exact=False):
    """adiabatic_flame_temperature with hisab 1.4.0's newton_raphson (tol ntol, 50 iterations,
    derivative threshold eps) and bisection (tol btol, 200 iterations) transcribed; FLAME_PATH[0] is
    the solver that answered. The keyword arguments are Rust's typed constants (the pinning check), and
    with `exact` the port's ADR 0002 #38 (the overflow audit's B6): the bracket's and the bisection's sign
    tests compare the signs of f(T₀) and f(T_max) exactly instead of testing their product."""
    FLAME_PATH[0] = None
    FLAME_ENDS[0] = None
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
    fa, fb = fun(t_lo), fun(t_hi)
    FLAME_ENDS[0] = (fa, fb)
    if exact:
        if not ((fa < 0.0 < fb) or (fb < 0.0 < fa)): return Err("ComputationError")
    elif not fa * fb < 0.0: return Err("ComputationError")
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
    if (((f_lo > 0.0 and f_hi > 0.0) or (f_lo < 0.0 and f_hi < 0.0)) if exact else f_lo * f_hi > 0.0):
        return Err("ComputationError")
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

# nuclear.rs (M4c): the constants from the golden corpus, Rust's private ones from the source literals (they are typed
# inline in src/nuclear.cyr; the pinning checks below perturb each), the ISOTOPES table from its literal.
nsrc = (root / "rust-old/src/nuclear.rs").read_text()
def n_private(name):
    m = re.search(rf"^const {name}: f64 = ([0-9_.e-]+);", nsrc, re.M)
    check(m is not None, f"nuclear.rs: no private const {name}")
    return float(m.group(1).replace("_", "")) if m else math.nan
N_TOL = n_private("EQUAL_RATE_REL_TOL")
SEMF = tuple(n_private(f"SEMF_{k}") for k in ("VOLUME", "SURFACE", "COULOMB", "ASYMMETRY", "PAIRING"))
check(N_TOL == 1e-12 and SEMF == (15.67, 17.23, 0.714, 93.15, 11.2), f"nuclear.rs: private consts {N_TOL}, {SEMF}")
AMU_MEV = b2f(golden_consts["nuclear::AMU_MEV"])
M_H = b2f(golden_consts["nuclear::HYDROGEN_ATOM_MASS_U"])
M_N = b2f(golden_consts["nuclear::NEUTRON_MASS_U"])
N_ISO = [(sym, int(z), int(a), float(hl.replace("_", "")), mode) for sym, z, a, hl, mode in re.findall(
    r'Isotope \{\s*symbol: "([^"]+)",\s*z: (\d+),\s*a: (\d+),\s*half_life_s: ([0-9_.e]+),(?:\s*//[^\n]*)?\s*'
    r'decay_mode: "([^"]+)",\s*\}', nsrc[nsrc.index("pub static ISOTOPES:"):])]
check(len(N_ISO) == 18 and len({r[0] for r in N_ISO}) == 18, f"nuclear.rs: ISOTOPES parsed as {len(N_ISO)} rows")
def iso_toks(r):
    return [("s", r[0], None), ("i", r[1], None), ("i", r[2], None), ("f", r[3], f2b(r[3])), ("s", r[4], None)]
def p_exp2(x):
    """glibc's exp2 (Python's math): the release oracle's 2.0_f64.powf(y), which LLVM rewrites to exp2(y)."""
    if math.isnan(x): return _prim(x, "exp2", (x,))
    try: return _prim(math.exp2(x), "exp2", (x,))
    except OverflowError: return _prim(math.inf, "exp2", (x,))
def r_decay_constant(t, ln2=LN2):
    stable = t == math.inf
    if not is_positive(t) and not stable: return Err("InvalidInput")
    return dv(ln2, t)
def r_atoms_remaining(n0, l, t): return n0 * p_exp(-l * t)
def r_after_half_lives(n0, n): return n0 * p_exp2(-n)
def r_time_for_fraction(fr, l):
    if not is_positive(fr) or fr > 1.0: return Err("InvalidInput")
    if not is_positive(l): return Err("InvalidInput")
    return dv(-p_ln(fr), l)
def r_mass_defect(z, n, m): return float(z) * M_H + float(n) * M_N - m
def r_bea(z, n, m):
    a = z + n
    if a > 0xFFFFFFFF: return Err("InvalidInput")        # Rust's checked_add
    if a == 0: return Err("InvalidInput")
    if not is_positive(m): return Err("InvalidInput")
    return dv(r_mass_defect(z, n, m) * AMU_MEV, float(a))
_semf_cbrt = {}
def r_semf(z, a, c=SEMF, four=4.0):
    """Rust 1.2's semi_empirical_binding_energy: NaN for z > a or a = 0 (the port's Err, ADR 0002 #3)."""
    if z > a or a == 0: return math.nan
    af, zf = float(a), float(z)
    asym = float(a - z) - zf
    cb = _semf_cbrt.get(a)
    if cb is None:
        cb = _semf_cbrt[a] = cbrt_cr(af)
    volume = c[0] * af
    surface = c[1] * cb * cb
    coulomb = dv(c[2] * zf * zf, cb)
    asymmetry = dv(c[3] * asym * asym, four * af)
    if a % 2: pairing = 0.0
    elif z % 2 == 0: pairing = dv(c[4], math.sqrt(af))
    else: pairing = dv(-c[4], math.sqrt(af))
    return volume - surface - coulomb - asymmetry + pairing
def r_q_value(a):
    r, p_ = slice_args(a)
    return (r_sum(r) - r_sum(p_)) * AMU_MEV
def r_chain(n0, l1, l2, t, tol=N_TOL):
    if not is_positive(l1) or not is_positive(l2): return Err("InvalidInput")
    if not is_non_negative(n0): return Err("InvalidInput")
    if not is_non_negative(t): return Err("InvalidInput")
    na = n0 * p_exp(-l1 * t)
    delta = abs(l2 - l1)
    lmin = r_min(l1, l2)
    growth = t if delta <= tol * r_max(l1, l2) else dv(-p_expm1(-delta * t), delta)
    nb = n0 * l1 * p_exp(-lmin * t) * growth
    nc = -n0 * p_expm1(-l1 * t) - nb
    return (na, nb, nc)
def r_secular(lp, ld, hundredth=0.01): return lp < hundredth * ld

# M4d: potential.rs, transcribed in its operation order (Rust's [f64; 3] points flattened, three floats each). Its typed
# private tolerances are parsed from the source; the replicas take them as parameters for the pinning checks.
psrc_pot = (root / "rust-old/src/potential.rs").read_text()
def pot_private(name):
    return float(re.search(rf"^const {name}: f64 = ([0-9_.e-]+);", psrc_pot, re.M).group(1).replace("_", ""))
P_TOL = pot_private("COINCIDENT_REL_TOL")
P_CTOL = pot_private("COLLINEAR_SIN_SQ_TOL")
KE_POT = b2f(golden_consts["potential::COULOMB_CONSTANT"])
def p_atan2(y, x):
    """glibc's atan2 (Python's math), Rust's y.atan2(x)."""
    return _prim(math.atan2(y, x), "atan2", (y, x))
def v_sub(a, b): return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
def v_cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def v_dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def r_lj(e, sg, r, four=4.0):
    if not is_positive(e) or not is_positive(sg) or not is_positive(r): return Err("InvalidInput")
    sr = dv(sg, r)
    sr6 = sr * sr * sr * sr * sr * sr
    sr12 = sr6 * sr6
    return four * e * (sr12 - sr6)
def r_lj_rmin(sg, six=6.0): return sg * p_exp2(1.0 / six)
def r_morse(d, a, re_, r):
    if not is_positive(d) or not is_positive(a): return Err("InvalidInput")
    x = 1.0 - p_exp(-a * (r - re_))
    return d * x * x
def r_width(w, mu, d, two=2.0):
    if not is_positive(d) or not is_positive(mu) or not is_positive(w): return Err("InvalidInput")
    return w * r_sqrt(dv(mu, two * d))
def r_coulomb(q1, q2, r):
    if not is_positive(r): return Err("InvalidInput")
    return dv(KE_POT * q1 * q2, r)
def r_coulomb_force(q1, q2, r):
    if not is_positive(r): return Err("InvalidInput")
    return dv(KE_POT * q1 * q2, r * r)
def r_distance(v):
    dx, dy, dz = v[3] - v[0], v[4] - v[1], v[5] - v[2]
    return r_sqrt(dx * dx + dy * dy + dz * dz)
def r_bond_angle(v, tol=None):
    a, b, c = v[0:3], v[3:6], v[6:9]
    t = P_TOL if tol is None else tol
    ba, bc = v_sub(a, b), v_sub(c, b)
    ba_sq, bc_sq = v_dot(ba, ba), v_dot(bc, bc)
    tol_sq = t * t
    b_sq = v_dot(b, b)
    if (not is_positive(ba_sq) or not is_positive(bc_sq) or ba_sq <= tol_sq * r_max(v_dot(a, a), b_sq)
            or bc_sq <= tol_sq * r_max(v_dot(c, c), b_sq)):
        return Err("InvalidInput")
    cr = v_cross(ba, bc)
    return p_atan2(r_sqrt(v_dot(cr, cr)), v_dot(ba, bc))
POT_MAPPED = [0]   # dihedral rows where atan2's −π became π (recomputed)
def r_dihedral(v, tol=None, pi=math.pi):
    a, b, c, d = v[0:3], v[3:6], v[6:9], v[9:12]
    t = P_CTOL if tol is None else tol
    b1, b2, b3 = v_sub(b, a), v_sub(c, b), v_sub(d, c)
    n1, n2 = v_cross(b1, b2), v_cross(b2, b3)
    n1_sq, n2_sq, b2_sq = v_dot(n1, n1), v_dot(n2, n2), v_dot(b2, b2)
    if (not is_positive(n1_sq) or not is_positive(n2_sq) or n1_sq <= t * v_dot(b1, b1) * b2_sq
            or n2_sq <= t * b2_sq * v_dot(b3, b3)):
        return Err("InvalidInput")
    phi = p_atan2(r_sqrt(b2_sq) * v_dot(b1, n2), v_dot(n1, n2))
    if phi <= -pi:
        POT_MAPPED[0] += 1
        return pi
    return phi

# M4e: spectroscopy.rs, transcribed in its operation order. The private constants are parsed from the source (the
# typed thresholds, which the replicas take as parameters for the pinning checks, and ELECTRON_PROTON_MASS_RATIO, whose
# const expression RYDBERG / (1.0 + ratio) is evaluated here as Rust's const evaluation does, in IEEE f64); the isotope
# table is the source's `match atomic_number` literal; hisab 1.4.0's fft (src/num/fft.rs) is transcribed for
# fft_spectrum, its twiddles glibc's cos and sin (the release oracle's sincos) through _prim.
ssrc = (root / "rust-old/src/spectroscopy.rs").read_text()
def s_private(name):
    return float(re.search(rf"^const {name}: f64 = ([0-9_.e-]+);", ssrc, re.M).group(1).replace("_", ""))
S_PRUNE = s_private("ISOTOPE_PRUNE_REL")
S_REPORT = s_private("ISOTOPE_REPORT_REL")
S_BIN = s_private("ISOTOPE_BIN_WIDTH_U")
S_RATIO = s_private("ELECTRON_PROTON_MASS_RATIO")
check(re.search(r"^const RYDBERG_HYDROGEN: f64 = RYDBERG / \(1\.0 \+ ELECTRON_PROTON_MASS_RATIO\);$", ssrc, re.M),
      "spectroscopy.rs: RYDBERG_HYDROGEN is not RYDBERG / (1.0 + ELECTRON_PROTON_MASS_RATIO)")
S_C = b2f(golden_consts["spectroscopy::SPEED_OF_LIGHT"])
S_RYD = b2f(golden_consts["spectroscopy::RYDBERG"])
S_RYD_EV = b2f(golden_consts["spectroscopy::RYDBERG_EV"])
S_QE = b2f(golden_consts["electrochemistry::ELEMENTARY_CHARGE"])
S_RYD_H = S_RYD / (1.0 + S_RATIO)
S_ISO, S_ISO_LIT = {}, {}
_iso_body = ssrc[ssrc.index("pub fn isotope_abundances("):]
_iso_body = _iso_body[_iso_body.index("match atomic_number {"):_iso_body.index("_ => &[],")]
for _z, _arm in re.findall(r"(\d+) => &\[(.*?)\]", _iso_body, re.S):
    S_ISO_LIT[int(_z)] = [(m_.replace("_", ""), a_.replace("_", ""))
                          for m_, a_ in re.findall(r"\(([0-9_.]+), ([0-9_.]+)\)", _arm)]
    S_ISO[int(_z)] = [(float(m_), float(a_)) for m_, a_ in S_ISO_LIT[int(_z)]]
# ADR 0002 #23 (M4e-4): the port's isotope table, scripts/oracle/isotopes.tsv (AME2020 masses, TICE 2009 abundances, ³⁶S
# included), parsed here independently of the generator: each `isotope` record's mass re-read from its verbatim AME2020
# line by the file's Fortran format (Z at columns 10-14, A 15-19, the element 21-23, the atomic mass's u 107-109 and
# micro-u 111-123), its abundance from the printed TICE value. S_ISO_NEW[z] is z's (mass, abundance) pairs, each the
# double nearest its decimal (Python's float() is correctly rounded); section 17 checks it against the oracle.
ISO23_PATH = root / "scripts/oracle/isotopes.tsv"
ISO23, ISO23_IV = [], []
for _l in ISO23_PATH.read_text(encoding="utf-8").splitlines():
    _f = _l.split("\t")
    if _f[0] == "isotope":
        check(len(_f) == 8, f"isotopes.tsv: {_l[:40]!r}: {len(_f)} fields")
        _z, _a, _sym, _mass, _ab, _tice, _line = int(_f[1]), int(_f[2]), _f[3], _f[4], _f[5], _f[6], _f[7]
        check(_line.isascii() and len(_line) >= 123 and int(_line[9:14]) == _z and int(_line[14:19]) == _a
              and _line[20:23].strip() == _sym, f"isotopes.tsv: {_a}{_sym}: the mas20 line is not Z {_z}, A {_a}, {_sym}")
        _u, _mu = _line[106:109].strip(), _line[110:123].strip()
        check(re.fullmatch(r"[0-9]{6}\.[0-9]+", _mu) is not None and _mass == f"{int(_u)}.{_mu.replace('.', '')}",
              f"isotopes.tsv: {_a}{_sym}: mass_u {_mass} is not the mas20 line's {_u} {_mu} micro-u")
        check(re.fullmatch(r"[0-9. ]+\([0-9]+\)", _tice) is not None and _ab == _tice.split("(")[0].replace(" ", ""),
              f"isotopes.tsv: {_a}{_sym}: abundance {_ab} is not the printed {_tice!r}")
        check(re.fullmatch(r"[0-9]+\.[0-9]+", _mass) is not None and re.fullmatch(r"0\.[0-9]+", _ab) is not None,
              f"isotopes.tsv: {_a}{_sym}: {_mass} / {_ab} not decimals")
        ISO23.append((_z, _a, _sym, _mass, _ab))
    elif _f[0] == "interval":
        check(len(_f) == 6, f"isotopes.tsv: {_l[:40]!r}: {len(_f)} fields")
        ISO23_IV.append((_f[1], Fraction(_f[2]), Fraction(_f[3])))
    else:
        check(_l.startswith("#"), f"isotopes.tsv: {_l[:40]!r} is no record or comment")
check(len(ISO23) > 0 and [(r[0], r[1]) for r in ISO23] == sorted({(r[0], r[1]) for r in ISO23}),
      "isotopes.tsv: the isotope records are not in Z then A order, or repeat one")
S_ISO_NEW = {}
for _z, _a, _sym, _mass, _ab in ISO23:
    S_ISO_NEW.setdefault(_z, []).append((float(_mass), float(_ab)))
S_SERIES_SRC = ssrc[ssrc.index("pub enum SpectralSeries {"):]
S_SERIES = re.findall(r"^    ([A-Z]\w*),$", S_SERIES_SRC[:S_SERIES_SRC.index("\n}")], re.M)
_nl = ssrc[ssrc.index("pub const fn n_lower("):]
S_N_LOWER = {v: int(x) for v, x in re.findall(r"Self::(\w+) => (\d+),", _nl[:_nl.index("\n    }\n")])}
def p_sin(x):
    if not math.isfinite(x): return _prim(math.nan, "sin", (x,))
    return _prim(math.sin(x), "sin", (x,))
def r_transmittance(a): return p_pow(10.0, -a)
def r_aft(tt):
    if not is_positive(tt) or tt > 1.0: return Err("InvalidInput")
    return -p_log10(tt)
def r_conc(a, e_, l_):
    if not is_positive(e_) or not is_positive(l_): return Err("InvalidInput")
    return dv(a, e_ * l_)
def r_by_wavelength(x, num):
    if not is_positive(x): return Err("InvalidInput")
    return dv(num, x)
def r_bohr(z, n):
    if z == 0 or n == 0: return Err("InvalidInput")
    return dv(-S_RYD_EV * float(z) * float(z), float(n) * float(n))
def r_transition(z, lo, hi):
    el = r_bohr(z, lo)
    if isinstance(el, Err): return el
    eu = r_bohr(z, hi)
    if isinstance(eu, Err): return eu
    if hi <= lo: return Err("InvalidInput")
    return eu - el
def r_rydberg_formula(r, z, lo, hi):
    if z == 0 or lo == 0 or hi == 0 or hi <= lo: return Err("InvalidInput")
    inv = r * float(z) * float(z) * (dv(1.0, float(lo) * float(lo)) - dv(1.0, float(hi) * float(hi)))
    return dv(1.0, inv)
def r_series(name, line, rh=None):
    nl = S_N_LOWER[name.split(".", 1)[1]]
    if line == 0 or nl + line > 0xFFFFFFFF: return Err("InvalidInput")
    return r_rydberg_formula(S_RYD_H if rh is None else rh, 1, nl, nl + line)
def r_wavenumber(x, hundredth=0.01): return r_by_wavelength(x, hundredth)
def r_reduced(m1, m2):
    if not is_positive(m1) or not is_positive(m2): return Err("InvalidInput")
    return dv(m1 * m2, m1 + m2)
def r_oscillator(n, w, two=2.0, half=0.5): return (float(n) + half) * dv(PLANCK, two * math.pi) * w
def r_vibrational(k, mu):
    if not is_positive(k) or not is_positive(mu): return Err("InvalidInput")
    return r_sqrt(dv(k, mu))
def r_rotor(j, b_, one=1.0): return PLANCK * b_ * float(j) * (float(j) + one)
def r_rot_const(i_, eight=8.0):
    if not is_positive(i_): return Err("InvalidInput")
    return dv(PLANCK, eight * math.pi * math.pi * i_)
def total_key(x):
    """Rust's total_cmp key of x: the bits as i64, the magnitude bits of a negative pattern flipped."""
    b = f2b(x)
    s = b - (1 << 64) if b >> 63 else b
    return s ^ 0x7FFFFFFFFFFFFFFF if s < 0 else s
def r_pattern(atoms, prune=None, report=None, width=None, ties_reversed=False, table=None):
    """Rust's isotope_pattern: the pruned convolution with a stable sort (Python's is stable) on total_cmp's key.
    ties_reversed puts equal masses in the reverse order (an unstable sort's possible result), for the check that the
    corpus pins the stable order. table is the isotope table, z → [(mass, abundance)]: Rust's literal by default,
    S_ISO_NEW for the port's (ADR 0002 #23)."""
    table = S_ISO if table is None else table
    prune = S_PRUNE if prune is None else prune
    report = S_REPORT if report is None else report
    width = S_BIN if width is None else width
    if not atoms: return Err("InvalidInput")
    pattern = [(0.0, 1.0)]
    for z, count in atoms:
        ab = table.get(z, [])
        if not ab:
            if z not in ELEMS: return Err("InvalidElement")
            shift = float(count) * ELEMS[z][3]
            pattern = [(m_ + shift, i_) for m_, i_ in pattern]
            continue
        for _ in range(count):
            scratch = [(m1 + m2, i1 * i2) for m1, i1 in pattern for m2, i2 in ab]
            if ties_reversed:
                scratch = [pr for _, pr in sorted(enumerate(scratch), key=lambda kp: (total_key(kp[1][0]), -kp[0]))]
            else:
                scratch.sort(key=lambda pr: total_key(pr[0]))
            pattern = []
            for m_, i_ in scratch:
                if pattern and abs(m_ - pattern[-1][0]) < width:
                    lm, li = pattern[-1]
                    total = li + i_
                    if total > 0.0: lm = (lm * li + m_ * i_) / total
                    pattern[-1] = (lm, total)
                    continue
                pattern.append((m_, i_))
            mx = 0.0
            for pr in pattern: mx = r_max(mx, pr[1])
            if mx > 0.0: pattern = [(m_, i_ / mx) for m_, i_ in pattern]
            pattern = [pr for pr in pattern if pr[1] >= prune]
    return Pairs([pr for pr in pattern if pr[1] > report])
def pattern_atoms(ins):
    n = ins[0][1]
    return [(ins[1 + 2 * k][1], ins[2 + 2 * k][1]) for k in range(n)]
def r_fft(sig, m2=-2.0):
    """Rust's fft_spectrum: hisab 1.4.0's num::fft transcribed (bit reversal, butterflies, w *= wn), then
    Complex::abs = sqrt(re·re + im·im) of the first n/2 + 1 bins."""
    n = len(sig)
    if n == 0 or n & (n - 1): return Err("InvalidInput")
    re_, im_ = list(sig), [0.0] * n
    if n > 1:
        j = 0
        for i_ in range(1, n):
            bit = n >> 1
            while j & bit:
                j ^= bit
                bit >>= 1
            j ^= bit
            if i_ < j:
                re_[i_], re_[j] = re_[j], re_[i_]
                im_[i_], im_[j] = im_[j], im_[i_]
        length = 2
        while length <= n:
            half = length >> 1
            angle = m2 * math.pi / float(length)
            wr, wi = p_cos(angle), p_sin(angle)
            for start in range(0, n, length):
                w_re, w_im = 1.0, 0.0
                for k in range(half):
                    a_, b_ = start + k, start + k + half
                    t_re = w_re * re_[b_] - w_im * im_[b_]
                    t_im = w_re * im_[b_] + w_im * re_[b_]
                    u_re, u_im = re_[a_], im_[a_]
                    re_[a_], im_[a_] = u_re + t_re, u_im + t_im
                    re_[b_], im_[b_] = u_re - t_re, u_im - t_im
                    w_re, w_im = w_re * wr - w_im * wi, w_re * wi + w_im * wr
            length <<= 1
    return [r_sqrt(re_[k] * re_[k] + im_[k] * im_[k]) for k in range(n // 2 + 1)]
def r_find_peaks(data, th, one=1.0):
    if len(data) < 3 or not (0.0 <= th <= 1.0): return Ints([])
    lo, hi = math.inf, -math.inf
    for v in data:
        if math.isfinite(v):
            lo, hi = r_min(lo, v), r_max(hi, v)
    if lo > hi: return Ints([])
    thr = lo * (one - th) + hi * th
    return Ints([k + 1 for k in range(len(data) - 2) if data[k + 1] > data[k] and data[k + 1] > data[k + 2]
                 and data[k + 1] > thr])
def peaks_args(ins):
    n = ins[0][1]
    return [t_[1] for t_ in ins[1:1 + n]], ins[1 + n][1]

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
    "nuclear::decay_constant": ("B", lambda a: r_decay_constant(*fl(a))),
    "nuclear::atoms_remaining": ("C", lambda a: r_atoms_remaining(*fl(a))),
    "nuclear::activity": ("B", lambda a: fl(a)[0] * fl(a)[1]),
    "nuclear::atoms_after_half_lives": ("C", lambda a: r_after_half_lives(*fl(a))),
    "nuclear::time_for_fraction": ("C", lambda a: r_time_for_fraction(*fl(a))),
    "nuclear::mass_defect": ("B", lambda a: r_mass_defect(*fl(a))),
    "nuclear::binding_energy_mev": ("B", lambda a: fl(a)[0] * AMU_MEV),
    "nuclear::binding_energy_per_nucleon": ("B", lambda a: r_bea(*fl(a))),
    "nuclear::semi_empirical_binding_energy": ("B", lambda a: r_semf(*fl(a))),
    "nuclear::q_value": ("B", r_q_value),
    "nuclear::decay_chain_two_step": ("C", lambda a: r_chain(*fl(a))),
    "nuclear::is_secular_equilibrium": ("A", lambda a: r_secular(*fl(a))),
    "potential::lennard_jones": ("B", lambda a: r_lj(*fl(a))),
    "potential::lj_equilibrium_distance": ("C", lambda a: r_lj_rmin(*fl(a))),
    "potential::morse": ("C", lambda a: r_morse(*fl(a))),
    "potential::morse_width_parameter": ("B", lambda a: r_width(*fl(a))),
    "potential::coulomb": ("B", lambda a: r_coulomb(*fl(a))),
    "potential::coulomb_force": ("B", lambda a: r_coulomb_force(*fl(a))),
    "potential::distance_3d": ("B", lambda a: r_distance(fl(a))),
    "potential::bond_angle": ("C", lambda a: r_bond_angle(fl(a))),
    "potential::dihedral_angle": ("C", lambda a: r_dihedral(fl(a))),
    "spectroscopy::absorbance": ("B", lambda a: fl(a)[0] * fl(a)[1] * fl(a)[2]),
    "spectroscopy::transmittance": ("C", lambda a: r_transmittance(fl(a)[0])),
    "spectroscopy::absorbance_from_transmittance": ("C", lambda a: r_aft(fl(a)[0])),
    "spectroscopy::concentration_from_absorbance": ("B", lambda a: r_conc(*fl(a))),
    "spectroscopy::photon_energy_from_frequency": ("B", lambda a: PLANCK * fl(a)[0]),
    "spectroscopy::photon_energy_from_wavelength": ("B", lambda a: r_by_wavelength(fl(a)[0], PLANCK * S_C)),
    "spectroscopy::wavelength_from_energy": ("B", lambda a: r_by_wavelength(fl(a)[0], PLANCK * S_C)),
    "spectroscopy::frequency_from_wavelength": ("B", lambda a: r_by_wavelength(fl(a)[0], S_C)),
    "spectroscopy::joules_to_ev": ("B", lambda a: dv(fl(a)[0], S_QE)),
    "spectroscopy::ev_to_joules": ("B", lambda a: fl(a)[0] * S_QE),
    "spectroscopy::bohr_energy_level": ("B", lambda a: r_bohr(*fl(a))),
    "spectroscopy::bohr_transition_energy": ("B", lambda a: r_transition(*fl(a))),
    "spectroscopy::rydberg_wavelength": ("B", lambda a: r_rydberg_formula(S_RYD, *fl(a))),
    "SpectralSeries::n_lower": ("A", lambda a: S_N_LOWER[fl(a)[0].split(".", 1)[1]]),
    "SpectralSeries::wavelength": ("B", lambda a: r_series(*fl(a))),
    "spectroscopy::wavenumber_from_wavelength": ("B", lambda a: r_wavenumber(fl(a)[0])),
    "spectroscopy::wavelength_from_wavenumber": ("B", lambda a: r_wavenumber(fl(a)[0])),
    "spectroscopy::reduced_mass": ("B", lambda a: r_reduced(*fl(a))),
    "spectroscopy::harmonic_oscillator_energy": ("B", lambda a: r_oscillator(*fl(a))),
    "spectroscopy::vibrational_frequency": ("B", lambda a: r_vibrational(*fl(a))),
    "spectroscopy::rigid_rotor_energy": ("B", lambda a: r_rotor(*fl(a))),
    "spectroscopy::rotational_constant": ("B", lambda a: r_rot_const(fl(a)[0])),
    "spectroscopy::diatomic_moment_of_inertia": ("B", lambda a: fl(a)[0] * fl(a)[1] * fl(a)[1]),
    "spectroscopy::isotope_abundances": ("A", lambda a: Pairs(S_ISO.get(fl(a)[0], []))),
    "spectroscopy::isotope_pattern": ("B", lambda a: r_pattern(pattern_atoms(a))),
    "spectroscopy::fft_spectrum": ("C", lambda a: r_fft(fl(a)[1:])),
    "spectroscopy::find_peaks": ("A", lambda a: r_find_peaks(*peaks_args(a))),
}
TABLE_RULES = {}
for name, v in PRESETS.items():
    TABLE_RULES[f"gas::{name}"] = [("f", x, f2b(x)) for x in v]
TABLE_RULES["phase::ANTOINE_DATA"] = [("v", len(ANTOINE), None)]
TABLE_RULES["thermochem::THERMOCHEM_DATA"] = [("v", len(TC), None)]
TABLE_RULES["thermochem::SHOMATE_DATA"] = [("v", len(SH), None)]
TABLE_RULES["electrochemistry::STANDARD_POTENTIALS"] = [("v", len(SP), None)]
TABLE_RULES["nuclear::ISOTOPES"] = [("v", len(N_ISO), None)]
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
for mod in ("reaction", "thermo", "phase", "gas", "thermochem", "solution", "electrochemistry", "kinetics", "nuclear",
            "potential", "spectroscopy"):
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
        if fn == "nuclear::ISOTOPES[]":
            k = ins[0][1]
            want = iso_toks(N_ISO[k]) if 0 <= k < len(N_ISO) else []
            check(cls == "A" and toks_equal(outs, want), f"{where}: ISOTOPES[{k}]")
            continue
        if fn == "nuclear::lookup_isotope":
            hit = next((r for r in N_ISO if r[0] == ins[0][1]), None)    # exact, case-sensitive, first in table order
            want = [("n", "", None)] if hit is None else iso_toks(hit)
            check(cls == "A" and toks_equal(outs, want), f"{where}: lookup_isotope({ins[0][1]!r})")
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
            if "[ADR 0002 #38:" in note: noted = None     # the oracle's Err; §24 checks the port's path in the note
            check(noted == path, f"{where}: solver path noted {noted!r}, recomputed {path!r}")
            kind = "err" if path is None else path.split(" ", 1)[0].rstrip(",")
            flame_paths[kind] = flame_paths.get(kind, 0) + 1
        if fn in ("thermochem::reaction_enthalpy", "thermochem::reaction_gibbs_energy"):
            adr1 = "[J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002 #1)]"
            nt_ = NOTES.get(where, "")
            check(nt_.endswith(adr1) or (adr1 + " [ADR 0002 #33: " in nt_ and nt_.endswith(")]")),
                  f"{where}: {fn} row note does not say ADR 0002 #1's × 1000 (then, at most, ADR 0002 #33's value)")
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
    "electrochemistry::lookup_half_reaction", "nuclear::ISOTOPES[]", "nuclear::lookup_isotope"}
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

# ── 13. organic.tsv (M4b) against the Rust source ───────────────────────────────────────
# The enums are parsed from rust-old/src/organic.rs (which proves the generated enums complete and in order: the
# generator's variant lists can't be matched exhaustively, the enums being #[non_exhaustive]), the const fns' match
# arms and matches! sets give each property, and predict_geometry's match is parsed and replayed; enthalpy_from_bonds
# is recomputed as Rust's two .sum()s from -0.0 and one subtraction.
osrc = (root / "rust-old/src/organic.rs").read_text()
def rust_enum(src, name):
    body = src[src.index(f"pub enum {name} {{"):]
    return re.findall(r"^    ([A-Z]\w*),$", body[:body.index("\n}")], re.M)
def rust_fn_body(src, item):
    body = src[src.index(item):]
    return body[:body.index("\n    }\n")]
def rust_arms(src, item):
    body = rust_fn_body(src, item)
    return {v: x for v, x in re.findall(r'Self::(\w+) => ("(?:[^"\\]|\\.)*"|[^,\n]+),', body)}
def rust_matches(src, item):
    body = rust_fn_body(src, item)
    return set(re.findall(r"Self::(\w+)", body[body.index("matches!("):]))
def snake_t(camel):
    """cyr.rs's snake (tanmatra's rule): `_` before a capital after a lowercase letter or digit, or before a capital
    that starts a capitalised word after a run of capitals (`CCDouble` → `CC_DOUBLE`, `HCl` → `H_CL`)."""
    out = []
    for k, ch in enumerate(camel):
        if ch.isupper() and k > 0:
            prev, nxt = camel[k - 1], camel[k + 1] if k + 1 < len(camel) else ""
            if prev.islower() or prev.isdigit() or (prev.isupper() and nxt.islower()):
                out.append("_")
        out.append(ch.upper())
    return "".join(out)
O_BONDS = rust_enum(osrc, "BondType")
O_GEOMS = rust_enum(osrc, "Geometry")
O_GROUPS = rust_enum(osrc, "FunctionalGroup")
check((len(O_BONDS), len(O_GEOMS), len(O_GROUPS)) == (45, 12, 16),
      f"organic.rs: enums parsed as {len(O_BONDS)}, {len(O_GEOMS)}, {len(O_GROUPS)} variants")
O_ENERGY = {v: float(x) for v, x in rust_arms(osrc, "pub const fn energy_kj(").items()}
O_ANGLE = {v: float(x) for v, x in rust_arms(osrc, "pub const fn ideal_bond_angle(").items()}
O_PATTERN = {v: x[1:-1] for v, x in rust_arms(osrc, "pub const fn pattern(").items()}
O_NONPOLAR = rust_matches(osrc, "pub const fn is_polar(")
O_DONOR = rust_matches(osrc, "pub const fn hydrogen_bonding(")
check(set(O_ENERGY) == set(O_BONDS) and set(O_ANGLE) == set(O_GEOMS) and set(O_PATTERN) == set(O_GROUPS)
      and O_NONPOLAR <= set(O_GROUPS) and O_DONOR <= set(O_GROUPS), "organic.rs: a property's arms are not its enum's variants")
check(all('\\' not in p for p in O_PATTERN.values()), "organic.rs: an escape in a pattern literal")
O_VSEPR = []   # (bonding, lone lo, lone hi, geometry) from predict_geometry's match
for b_, lo, hi, g in re.findall(r"\((\d+), (\d+)(?:\.\.=(\d+))?\) => Some\(Geometry::(\w+)\)",
                                rust_fn_body(osrc, "pub fn predict_geometry(")):
    O_VSEPR.append((int(b_), int(lo), int(hi or lo), g))
check(len(O_VSEPR) == 14 and all(g in O_GEOMS for *_, g in O_VSEPR), f"organic.rs: predict_geometry parsed as {O_VSEPR}")
def r_predict(bonding, lone):
    for b_, lo, hi, g in O_VSEPR:          # first matching arm, as Rust's match
        if bonding == b_ and lo <= lone <= hi:
            return g
    return None
def r_bond_sum(ins, k):
    """Rust's .map(|(b, n)| b.energy_kj() * *n as f64).sum() over the list at token k: (sum, tokens used)."""
    n = ins[k][1]
    check(ins[k][0] == "v", "organic: a bond list without its v: length")
    s_ = -0.0
    for j in range(n):
        (ek, ev, _), (ik, iv, _) = ins[k + 1 + 2 * j], ins[k + 2 + 2 * j]
        check(ek == "e" and ev.startswith("BondType.") and ik == "i" and 0 <= iv <= 0xFFFFFFFF,
              f"organic: bond list token {ev!r} {iv!r}")
        s_ += O_ENERGY[ev.split(".", 1)[1]] * float(iv)
    return s_, 1 + 2 * n
def o_enum(tok, enum, variants):
    check(tok[0] == "e" and tok[1].startswith(enum + ".") and tok[1].split(".", 1)[1] in variants,
          f"organic: {tok!r} is no {enum} variant")
    return tok[1].split(".", 1)[1]
ocount, oseen, opairs = {}, {}, set()
for fn, cls, ins, outs, where in R.get("organic", []):
    ocount[fn] = ocount.get(fn, 0) + 1
    if fn == "organic::enthalpy_from_bonds":
        check(cls == "B", f"{where}: {fn} class {cls}")
        sb, used = r_bond_sum(ins, 0)
        sf, used2 = r_bond_sum(ins, used)
        check(used + used2 == len(ins), f"{where}: {fn} has {len(ins)} input tokens")
        check(toks_equal(outs, out_toks(sb - sf)), f"{where}: {fn} -> {outs}, recomputed {sb - sf!r}")
        continue
    check(cls == "A", f"{where}: {fn} class {cls}")
    if fn == "organic::predict_geometry":
        b_, lp = ins[0][1], ins[1][1]
        check(ins[0][0] == ins[1][0] == "i" and 0 <= b_ <= 255 and 0 <= lp <= 255, f"{where}: inputs {ins}")
        opairs.add((b_, lp))
        g = r_predict(b_, lp)
        check(toks_equal(outs, [("n", "", None)] if g is None else [("e", f"Geometry.{g}", None)]),
              f"{where}: predict_geometry({b_}, {lp}) -> {outs}, recomputed {g}")
        continue
    if fn == "BondType::energy_kj":
        v = o_enum(ins[0], "BondType", O_BONDS)
        want = [("f", O_ENERGY[v], f2b(O_ENERGY[v]))]
    elif fn == "Geometry::ideal_bond_angle":
        v = o_enum(ins[0], "Geometry", O_GEOMS)
        want = [("f", O_ANGLE[v], f2b(O_ANGLE[v]))]
    elif fn == "FunctionalGroup::pattern":
        v = o_enum(ins[0], "FunctionalGroup", O_GROUPS)
        want = [("s", O_PATTERN[v], None)]
    elif fn == "FunctionalGroup::is_polar":
        v = o_enum(ins[0], "FunctionalGroup", O_GROUPS)
        want = [("i", int(v not in O_NONPOLAR), None)]
    elif fn == "FunctionalGroup::hydrogen_bonding":
        v = o_enum(ins[0], "FunctionalGroup", O_GROUPS)
        want = [("i", int(v in O_DONOR), None)]
    else:
        check(False, f"{where}: organic.tsv function {fn} has no rule here")
        continue
    oseen.setdefault(fn, []).append(v)
    check(toks_equal(outs, want), f"{where}: {fn}({v}) -> {outs}, want {want}")
for fn, variants in (("BondType::energy_kj", O_BONDS), ("Geometry::ideal_bond_angle", O_GEOMS),
                     ("FunctionalGroup::pattern", O_GROUPS), ("FunctionalGroup::is_polar", O_GROUPS),
                     ("FunctionalGroup::hydrogen_bonding", O_GROUPS)):
    check(oseen.get(fn) == variants, f"organic.tsv: {fn} rows are not every variant once, in declaration order")
edges_o = [0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 127, 128, 254, 255]
need = {(x, y) for x in range(16) for y in range(16)} | {(x, e_) for e_ in edges_o for x in range(256)} \
    | {(e_, x) for e_ in edges_o for x in range(256)}
check(opairs == need, f"organic.tsv: predict_geometry covers {len(opairs)} pairs, the grid is {len(need)}")
check(ocount.get("organic::enthalpy_from_bonds", 0) >= 1000, "organic.tsv: fewer than 1000 enthalpy_from_bonds rows")
print(f"organic rows checked: {sum(ocount.values())} ({', '.join(f'{k} {v}' for k, v in ocount.items())}) against "
      f"organic.rs's {len(O_BONDS)} + {len(O_GEOMS)} + {len(O_GROUPS)} variants and predict_geometry's {len(O_VSEPR)} arms")

# src/organic_data.cyr against the Rust source: the enums (ADR 0002 #6's rename), the columns, the UTF-8 pattern blob,
# the bool columns, the match's section comments and notes, and the acceptor column against its stated rule.
otext = (root / "src/organic_data.cyr").read_text()
rel_o = "src/organic_data.cyr"
def cyr_enum(text, name):
    m = re.search(rf"^enum {name} \{{(.*?)\n\}}", text, re.M | re.S)
    return [] if m is None else [(n, int(v)) for n, v in re.findall(r"^    (KM_[A-Z0-9_]+) = (\d+);", m.group(1), re.M)]
def cyr_bools(text, rel, name):
    m = re.search(rf"^var {name}: bool\[(\d+)\] = \{{(.*?)\}};", text, re.M | re.S)
    check(m is not None, f"{rel}: no {name}")
    if m is None: return []
    vals = [x == "true" for line in m.group(2).splitlines() for x in re.findall(r"\b(true|false),", line.split("#", 1)[0])]
    check(len(vals) == int(m.group(1)), f"{rel}: {name} declares {m.group(1)}, holds {len(vals)}")
    return vals
for enum, cyr_name, prefix, variants, count in (("BondType", "KmBondType", "KM_BOND_TYPE_", O_BONDS, "KM_BOND_TYPE_COUNT"),
                                                 ("Geometry", "KmGeometry", "KM_GEOMETRY_", O_GEOMS, "KM_GEOMETRY_COUNT"),
                                                 ("FunctionalGroup", "KmFunctionalGroup", "KM_FUNCTIONAL_GROUP_", O_GROUPS,
                                                  "KM_FUNCTIONAL_GROUP_COUNT")):
    want = [("KM_GEOMETRY_T_SHAPED" if v == "TShapedThreeLonePairs" else prefix + snake_t(v), k)
            for k, v in enumerate(variants)]
    check(cyr_enum(otext, cyr_name) == want, f"{rel_o}: enum {cyr_name} is not {enum}'s variants in order")
    check(re.search(rf"^const {count} = {len(variants)};$", otext, re.M), f"{rel_o}: {count}")
check(cyr_col_in(otext, rel_o, "KM_BOND_TYPE_ENERGY_KJ") == [O_ENERGY[v] for v in O_BONDS], f"{rel_o}: energies")
check(cyr_col_in(otext, rel_o, "KM_GEOMETRY_IDEAL_BOND_ANGLE") == [O_ANGLE[v] for v in O_GEOMS], f"{rel_o}: angles")
got_s, one_each = cyr_strings_in(otext, rel_o, "KM_FUNCTIONAL_GROUP_PATTERN", "utf-8")
check(got_s == [O_PATTERN[v] for v in O_GROUPS] and one_each, f"{rel_o}: the pattern blob is not organic.rs's strings")
check(cyr_bools(otext, rel_o, "KM_FUNCTIONAL_GROUP_IS_POLAR") == [v not in O_NONPOLAR for v in O_GROUPS], f"{rel_o}: is_polar")
check(cyr_bools(otext, rel_o, "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_DONOR") == [v in O_DONOR for v in O_GROUPS],
      f"{rel_o}: the donor column is not hydrogen_bonding's set")
# ADR 0002 #6's acceptor column has no oracle: its stated rule, IUPAC 2011's electron-rich region (a lone pair or a
# π-bonded pair, the weak acceptors included: the user's decision), is recomputed from the patterns — an atom that
# carries lone pairs in these groups (O, N, S, a halogen: F, Cl, Br, I or the placeholder X) or a π bond (=, ≡, or an
# aryl ring, Ar / C₆H₅).
def has_acceptor_site(pattern):
    atoms = set(re.findall(r"[A-Z][a-z]?", pattern))
    return bool(atoms & {"O", "N", "S", "F", "Cl", "Br", "I", "X"}) or "=" in pattern or "≡" in pattern \
        or "Ar" in atoms or "C₆H₅" in pattern
acc = [has_acceptor_site(O_PATTERN[v]) for v in O_GROUPS]
check(cyr_bools(otext, rel_o, "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_ACCEPTOR") == acc,
      f"{rel_o}: the acceptor column is not 'the pattern has a lone-pair atom or a π bond'")
check(all(a for v, a in zip(O_GROUPS, acc) if v in O_DONOR), f"{rel_o}: a donor that is no acceptor")
sections = re.findall(r"^\s+// (.+)$", rust_fn_body(osrc, "pub const fn energy_kj("), re.M)
ebody = re.search(r"^var KM_BOND_TYPE_ENERGY_KJ: f64\[45\] = \{(.*?)\};", otext, re.M | re.S).group(1)
check(re.findall(r"^    # (.+)$", ebody, re.M) == sections and len(sections) == 9,
      f"{rel_o}: energy_kj's section comments are not carried into the initializer")
for v, note in re.findall(r"Self::(\w+) => [0-9.]+, // (.+)$", rust_fn_body(osrc, "pub const fn ideal_bond_angle("), re.M):
    check(f"# {v}: {note}" in otext, f"{rel_o}: ideal_bond_angle's note on {v} is not carried")
print(f"{rel_o} checked against organic.rs (3 enums, ADR 0002 #6's rename; energies, angles, the UTF-8 pattern blob, "
      f"polar and donor flags; the acceptor column by its rule: {sum(acc)} of {len(acc)} groups, all but "
      f"{', '.join(v for v, a in zip(O_GROUPS, acc) if not a)})")

# ── 14. inorganic.tsv (M4b) against the Rust source ─────────────────────────────────────
# KmLigand and field_strength from the enum and its match, the four private d-electron tables and IONIC_RADII from
# their literals, the formulas transcribed in Rust's operation order (Python floats are IEEE doubles), the guards in
# Rust's order. Every typed formula constant must change some golden row when perturbed (the corpus pins it).
isrc = (root / "rust-old/src/inorganic.rs").read_text()
I_LIGANDS = rust_enum(isrc, "Ligand")
I_STRENGTH = {v: int(x) for v, x in rust_arms(isrc, "pub const fn field_strength(").items()}
check(len(I_LIGANDS) == 11 and set(I_STRENGTH) == set(I_LIGANDS), f"inorganic.rs: Ligand parsed as {I_LIGANDS}")
I_THRESHOLD = int(re.search(r"self\.field_strength\(\) >= (\d+)", isrc).group(1))
def literal_u8s(name):
    body = isrc[isrc.index(f"const {name}:"):]
    body = body[body.index("= [") + 3:body.index("];")]
    return [int(x) for x in re.findall(r"\d+", body)]
I_OCC = {hs: list(zip(*[iter(literal_u8s(n))] * 2)) for hs, n in ((True, "OCT_OCCUPANCY_HIGH_SPIN"), (False, "OCT_OCCUPANCY_LOW_SPIN"))}
I_UNP = {True: literal_u8s("OCT_UNPAIRED_HIGH_SPIN"), False: literal_u8s("OCT_UNPAIRED_LOW_SPIN")}
check(all(len(I_OCC[h]) == 11 and len(I_UNP[h]) == 11 for h in (True, False))
      and all(t + e_ == n for h in (True, False) for n, (t, e_) in enumerate(I_OCC[h])),
      "inorganic.rs: the d-electron tables are not 11 rows of dⁿ")
I_RADII = [(sym, int(q), int(cn), float(r)) for sym, q, cn, r in re.findall(
    r'IonicRadius \{\s*symbol: "(\w+)",\s*charge: (-?\d+),\s*(?://[^\n]*\s*)*coordination: (\d+),\s*radius_pm: ([0-9.]+),\s*\}',
    isrc[isrc.index("pub static IONIC_RADII:"):])]
check(len(I_RADII) == 30, f"inorganic.rs: IONIC_RADII parsed as {len(I_RADII)} rows")
KE = b2f(golden_consts["potential::COULOMB_CONSTANT"])
NA = b2f(golden_consts["element::AVOGADRO"])
QE = b2f(golden_consts["electrochemistry::ELEMENTARY_CHARGE"])
INVALID = Err("InvalidInput")
def r_cfse(d, hs, a=-0.4, b_=0.6):
    if not 0 <= d <= 10: return INVALID
    t2g, eg = I_OCC[hs][d]
    return float(t2g) * a + float(eg) * b_
def r_cfse_exact(d, hs, five=5):
    """ADR 0002 #12, the port's CFSE: the exact (3·eg − 2·t2g)/5, correctly rounded (Python's Fraction → float is)."""
    t2g, eg = I_OCC[hs][d]
    return float(Fraction(3 * eg - 2 * t2g) / Fraction(five))
def r_unpaired(d, hs):
    return I_UNP[hs][d] if 0 <= d <= 10 else INVALID
def r_tet(x, four=4.0, nine=9.0): return dv(x * four, nine)
def r_spin(n, two=2.0):
    n = float(n)
    return math.sqrt(n * (n + two))
def r_lande(m, zc, za, r0, n, one=1.0):
    if not is_positive(r0): return INVALID
    if not is_finite(n) or n <= 1.0: return INVALID
    if not is_positive(m): return INVALID
    if zc == 0 or za == 0: return INVALID
    u = -KE * NA * m * float(zc) * float(za) * QE * QE
    return dv(u, r0) * (one - dv(one, n))
def r_haber(h, s_, ie, d, ea, two=2.0): return h - s_ - ie - dv(d, two) - ea
def radius_toks(r):
    return [("s", r[0], None), ("i", r[1], None), ("i", r[2], None), ("f", r[3], f2b(r[3]))]
def i_bool(tok):
    check(tok[0] == "i" and tok[1] in (0, 1), f"inorganic: {tok!r} is no bool token")
    return tok[1] == 1
I_RULES = {
    "inorganic::cfse_octahedral": ("A", lambda a: r_cfse(a[0][1], i_bool(a[1]))),
    "inorganic::unpaired_electrons_oct": ("A", lambda a: r_unpaired(a[0][1], i_bool(a[1]))),
    "inorganic::delta_tet_from_oct": ("B", lambda a: r_tet(a[0][1])),
    "inorganic::spin_only_moment": ("B", lambda a: r_spin(a[0][1])),
    "inorganic::born_lande_lattice_energy": ("B", lambda a: r_lande(*fl(a))),
    "inorganic::born_haber_lattice_energy": ("B", lambda a: r_haber(*fl(a))),
}
icount, igrid = {}, {}
cfse_exact_rows = [0, 0]     # Ok cfse rows where the exact value equals / differs from the oracle's (ADR 0002 #12)
for fn, cls, ins, outs, where in R.get("inorganic", []):
    icount[fn] = icount.get(fn, 0) + 1
    if fn == "inorganic::IONIC_RADII":
        check(cls == "A" and ins == [] and toks_equal(outs, [("v", len(I_RADII), None)]), f"{where}: {fn}")
        continue
    if fn == "inorganic::IONIC_RADII[]":
        k = ins[0][1]
        check(cls == "A" and 0 <= k < len(I_RADII) and toks_equal(outs, radius_toks(I_RADII[k])), f"{where}: IONIC_RADII[{k}]")
        igrid.setdefault(fn, set()).add(k)
        continue
    if fn == "inorganic::lookup_ionic_radius":
        sym, q = ins[0][1], ins[1][1]
        check(ins[1][0] == "i" and -128 <= q <= 127, f"{where}: charge {q!r} is no i8")
        hit = next((r for r in I_RADII if r[0] == sym and r[1] == q), None)     # exact, case-sensitive, first match
        check(cls == "A" and toks_equal(outs, [("n", "", None)] if hit is None else radius_toks(hit)),
              f"{where}: lookup_ionic_radius({sym!r}, {q})")
        igrid.setdefault(fn, set()).add((sym, q))
        continue
    if fn in ("Ligand::field_strength", "Ligand::is_strong_field"):
        v = o_enum(ins[0], "Ligand", I_LIGANDS)
        want = I_STRENGTH[v] if fn == "Ligand::field_strength" else I_STRENGTH[v] >= I_THRESHOLD
        check(cls == "A" and toks_equal(outs, out_toks(want)), f"{where}: {fn}({v}) -> {outs}")
        igrid.setdefault(fn, []).append(v)
        continue
    if fn not in I_RULES:
        check(False, f"{where}: inorganic.tsv function {fn} has no rule here")
        continue
    want_cls, rule = I_RULES[fn]
    check(cls == want_cls, f"{where}: {fn} class {cls}")
    if fn in ("inorganic::cfse_octahedral", "inorganic::unpaired_electrons_oct"):
        check(ins[0][0] == "i" and 0 <= ins[0][1] <= 255, f"{where}: d = {ins[0]!r} is no u8")
        igrid.setdefault(fn, set()).add((ins[0][1], ins[1][1]))
    if fn == "inorganic::spin_only_moment":
        check(ins[0][0] == "i" and 0 <= ins[0][1] <= 255, f"{where}: n = {ins[0]!r} is no u8")
        igrid.setdefault(fn, set()).add(ins[0][1])
    if fn == "inorganic::born_lande_lattice_energy":
        check(all(ins[k][0] == "i" and 0 <= ins[k][1] <= 0xFFFFFFFF for k in (1, 2)), f"{where}: a charge is no u32")
    got = rule(ins)
    check(toks_equal(outs, out_toks(got)), f"{where}: {fn}{tuple(fl(ins))} -> {[t[1] for t in outs]}, recomputed {got!r}")
    if fn == "inorganic::cfse_octahedral" and not isinstance(got, Err):
        # ADR 0002 #12: the port returns the exact value; the oracle's is within a few ulps of the same fifth (the
        # harness's mapping, round(5v)), and the row's note names the exact value exactly where the two differ.
        d, hs = ins[0][1], i_bool(ins[1])
        exact = r_cfse_exact(d, hs)
        t2g, eg = I_OCC[hs][d]
        q = round(got * 5.0)
        check(q == 3 * eg - 2 * t2g and abs(got * 5.0 - q) < 1e-9, f"{where}: the oracle's {got!r} is not near (3·eg − 2·t2g)/5")
        tag = f"[ADR 0002 #12: the port returns the exact (3·e_g − 2·t2g)/5 = {exact!r} (f:{f2b(exact):016x})]"
        noted = NOTES.get(where, "").endswith(tag)
        differs = f2b(exact) != f2b(got)
        check(noted == differs and ("ADR 0002 #12" in NOTES.get(where, "")) == differs,
              f"{where}: cfse({d}, {hs}) exact {exact!r}, oracle {got!r}: the ADR 0002 #12 note is {'missing' if differs else 'stray'}")
        cfse_exact_rows[differs] += 1
check(cfse_exact_rows == [9, 13], f"inorganic: ADR 0002 #12 rows (equal, differing) = {cfse_exact_rows}, not [9, 13]")
check(igrid.get("inorganic::cfse_octahedral") == igrid.get("inorganic::unpaired_electrons_oct")
      == {(d, h) for d in range(256) for h in (0, 1)}, "inorganic.tsv: the d-electron functions don't cover every u8 × spin")
check(igrid.get("inorganic::spin_only_moment") == set(range(256)), "inorganic.tsv: spin_only_moment doesn't cover every u8")
check(igrid.get("inorganic::IONIC_RADII[]") == set(range(len(I_RADII))), "inorganic.tsv: an IONIC_RADII row is not probed")
check({(r[0], q) for r in I_RADII for q in range(-128, 128)} <= igrid.get("inorganic::lookup_ionic_radius", set()),
      "inorganic.tsv: lookup_ionic_radius doesn't probe every table symbol × every i8")
check(igrid.get("Ligand::field_strength") == I_LIGANDS == igrid.get("Ligand::is_strong_field"),
      "inorganic.tsv: the Ligand rows are not every variant once, in declaration order")
print(f"inorganic rows checked: {sum(icount.values())} ({', '.join(f'{k} {v}' for k, v in icount.items())}) against "
      f"inorganic.rs's Ligand ({len(I_LIGANDS)}), d-electron tables (2 × 11 + 2 × 11) and IONIC_RADII ({len(I_RADII)})")
# The typed formula constants (src/inorganic.cyr has them inline, as Rust): each must change some golden row.
def changes(fn, rule):
    for f_, cls, ins, outs, where in R.get("inorganic", []):
        if f_ == fn and not toks_equal(outs, out_toks(rule(ins))):
            return True
    return False
for label, fn, rule in (
        # The oracle's −0.4 and 0.6 (this replica's; since ADR 0002 #12 the port computes (3·eg − 2·t2g)/5).
        ("cfse −0.4 + 1 ulp", "inorganic::cfse_octahedral", lambda a: r_cfse(a[0][1], i_bool(a[1]), a=nextup(-0.4))),
        ("cfse 0.6 − 1 ulp", "inorganic::cfse_octahedral", lambda a: r_cfse(a[0][1], i_bool(a[1]), b_=nextdown(0.6))),
        ("delta_tet 4.0 + 1 ulp", "inorganic::delta_tet_from_oct", lambda a: r_tet(a[0][1], four=nextup(4.0))),
        ("delta_tet 9.0 + 1 ulp", "inorganic::delta_tet_from_oct", lambda a: r_tet(a[0][1], nine=nextup(9.0))),
        ("spin_only 2.0 + 1 ulp", "inorganic::spin_only_moment", lambda a: r_spin(a[0][1], two=nextup(2.0))),
        ("Born–Landé 1.0 + 1 ulp", "inorganic::born_lande_lattice_energy", lambda a: r_lande(*fl(a), one=nextup(1.0))),
        ("Born–Haber 2.0 + 1 ulp", "inorganic::born_haber_lattice_energy", lambda a: r_haber(*fl(a), two=nextup(2.0)))):
    check(changes(fn, rule), f"inorganic: {label} changes no golden row: the typed constant is not pinned")
# The port's divisor 5 (ADR 0002 #12, typed inline): another divisor must change some row's exact value (its note).
check(any(r_cfse_exact(d, h, five=6) != r_cfse_exact(d, h) for d in range(11) for h in (True, False)),
      "inorganic: the exact CFSE's divisor 5 is not pinned")
strong = [(f_, ins, outs) for f_, cls, ins, outs, w in R.get("inorganic", []) if f_ == "Ligand::is_strong_field"]
for thr, observable in ((I_THRESHOLD + 1, True), (I_THRESHOLD - 1, False)):
    moved = any(outs[0][1] != int(I_STRENGTH[ins[0][1].split(".")[1]] >= thr) for _, ins, outs in strong)
    check(moved == observable, f"inorganic: is_strong_field's threshold {thr} {'changes' if moved else 'changes no'} row")
print(f"inorganic: ADR 0002 #12 exact CFSE: {cfse_exact_rows[1]} Ok rows where (3·eg − 2·t2g)/5 differs from the "
      f"oracle's t2g·(−0.4) + eg·0.6 (noted; the oracle within 1e-9 of the same fifth), {cfse_exact_rows[0]} equal")
print(f"inorganic: typed constants pinned (the oracle's −0.4, 0.6; 4.0, 9.0, the 2.0s, 1.0; the strong-field threshold "
      f"{I_THRESHOLD} from above: {I_THRESHOLD - 1} is unobservable, no ligand has rank {I_THRESHOLD - 1})")

# src/inorganic_data.cyr against the Rust source.
itext = (root / "src/inorganic_data.cyr").read_text()
rel_i = "src/inorganic_data.cyr"
check(cyr_enum(itext, "KmLigand") == [("KM_LIGAND_" + snake_t(v), k) for k, v in enumerate(I_LIGANDS)],
      f"{rel_i}: enum KmLigand is not Ligand's variants in order")
check(re.search(r"^const KM_LIGAND_COUNT = 11;$", itext, re.M), f"{rel_i}: KM_LIGAND_COUNT")
check(cyr_col_in(itext, rel_i, "KM_LIGAND_FIELD_STRENGTH", int) == [I_STRENGTH[v] for v in I_LIGANDS], f"{rel_i}: strengths")
for hs, name in ((True, "HIGH_SPIN"), (False, "LOW_SPIN")):
    check(cyr_col_in(itext, rel_i, f"_KM_OCT_OCCUPANCY_{name}_T2G", int) == [p[0] for p in I_OCC[hs]], f"{rel_i}: {name} t2g")
    check(cyr_col_in(itext, rel_i, f"_KM_OCT_OCCUPANCY_{name}_EG", int) == [p[1] for p in I_OCC[hs]], f"{rel_i}: {name} eg")
    check(cyr_col_in(itext, rel_i, f"_KM_OCT_UNPAIRED_{name}", int) == I_UNP[hs], f"{rel_i}: {name} unpaired")
nk = 0
for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = (.+);$", itext, re.M):
    path, name, lit = m.groups()
    check(name == "KM_" + path.split("::")[1] and path.startswith("inorganic::"), f"{rel_i}: {name} names kimiya::{path}")
    check(re.fullmatch(r"-?[0-9]+\.[0-9]+", lit) and f2b(float(lit)) == golden_consts.get(path) == f2b(lits.get(path, math.nan)),
          f"{rel_i}: {name} = {lit} is not the golden and literal {path}")
    nk += 1
check(nk == 10, f"{rel_i}: {nk} constants, not the 10 MADELUNG_* / BORN_EXP_*")
for heading in ("Madelung constants (dimensionless), referred to the nearest cation–anion",
                "Born exponents n (dimensionless) by the ion's noble-gas electron"):
    check(f"# {heading}" in itext, f"{rel_i}: the section comment {heading!r} is not carried")
check(re.search(r"^const KM_IONIC_RADII_COUNT = 30;$", itext, re.M), f"{rel_i}: KM_IONIC_RADII_COUNT")
got_s, one_each = cyr_strings_in(itext, rel_i, "KM_IONIC_RADII_SYMBOL", "ascii")
check(got_s == [r[0] for r in I_RADII] and one_each, f"{rel_i}: the symbol blob is not inorganic.rs's symbols")
check(cyr_col_in(itext, rel_i, "KM_IONIC_RADII_CHARGE", int) == [r[1] for r in I_RADII], f"{rel_i}: charges")
check(cyr_col_in(itext, rel_i, "KM_IONIC_RADII_COORDINATION", int) == [r[2] for r in I_RADII], f"{rel_i}: coordination")
check(cyr_col_in(itext, rel_i, "KM_IONIC_RADII_RADIUS_PM") == [r[3] for r in I_RADII], f"{rel_i}: radii")
check("Shannon (1976) lists Sn²⁺ only at CN 8." in itext, f"{rel_i}: the literal's Sn²⁺ note is not carried")
print(f"{rel_i} checked against inorganic.rs (KmLigand, field strengths, the 4 d-electron tables, 10 constants, "
      f"30 IONIC_RADII rows: symbol blob, charge, coordination, radius)")
# ── 15. nuclear.tsv (M4c): beyond the row recomputation ─────────────────────────────────
# Every nuclear row was recomputed in the M2/M3 loop above (Rust 1.2.0's nuclear.rs transcribed in its operation order;
# glibc's exp, expm1, ln and exp2 via Python, and the corrected oracle through _prim). Here: the grids' coverage, ADR 0002
# #3's notes, the typed formula constants' pinning, the generated data file against the literal, and the M4b policy's
# diagnostics of the formulas that miss an exact value (ADR 0002, "the policy for inexact formulas").
nrows = {}
for fn, cls, ins, outs, where in R.get("nuclear", []):
    nrows.setdefault(fn, []).append((ins, outs, where))
semf_pairs, adr3 = set(), 0
SEMF_TAG = "[ADR 0002 #3: the port returns Err(InvalidInput)]"
for ins, outs, where in nrows.get("nuclear::semi_empirical_binding_energy", []):
    z, a = ins[0][1], ins[1][1]
    check(ins[0][0] == ins[1][0] == "i" and 0 <= z <= 0xFFFFFFFF and 0 <= a <= 0xFFFFFFFF, f"{where}: inputs {ins}")
    semf_pairs.add((z, a))
    nan = outs[0][0] == "f" and isnan_b(outs[0][2])
    check(nan == (z > a or a == 0), f"{where}: SEMF({z}, {a}) is {'NaN' if nan else 'finite'}")
    check(NOTES.get(where, "").endswith(SEMF_TAG) == nan, f"{where}: the ADR 0002 #3 note is {'missing' if nan else 'stray'}")
    adr3 += nan
check({(z, a) for z in range(1, 301) for a in range(z, min(3 * z, 300) + 1)} <= semf_pairs,
      "nuclear.tsv: semi_empirical_binding_energy doesn't cover every Z ≤ A ≤ 3Z with A ≤ 300 (port plan §4)")
check({(0xFFFFFFFF, 0xFFFFFFFF), (0, 0), (0xFFFFFFFF, 0)} <= semf_pairs, "nuclear.tsv: SEMF misses a u32 extreme")
check({ins[0][1] for ins, _, _ in nrows.get("nuclear::ISOTOPES[]", [])} == set(range(len(N_ISO))),
      "nuclear.tsv: an ISOTOPES row is not probed")
check({r[0] for r in N_ISO} <= {ins[0][1] for ins, _, _ in nrows.get("nuclear::lookup_isotope", [])},
      "nuclear.tsv: lookup_isotope doesn't probe every table symbol")
for fn in ("nuclear::mass_defect", "nuclear::binding_energy_per_nucleon"):
    for ins, outs, where in nrows.get(fn, []):
        check(all(ins[k][0] == "i" and 0 <= ins[k][1] <= 0xFFFFFFFF for k in (0, 1)), f"{where}: a count is no u32")
# The chain's branches: the equal-rate limit (λ₁ = λ₂ exactly, and within 1e-12), the rows on the threshold itself (where
# 1e-12 ± 1 ulp moves the branch), and the rows whose e^(−λ_min·t) is subnormal (ADR 0006 item 6).
def chain_equal(a, tol=N_TOL):
    """The branch Rust's chain takes on inputs a: True for the equal-rate limit, None for an Err row."""
    n0, l1, l2, t = a
    if isinstance(r_chain(n0, l1, l2, t), Err): return None
    return abs(l2 - l1) <= tol * r_max(l1, l2)
ok_chain = [fl(ins) for ins, outs, w in nrows.get("nuclear::decay_chain_two_step", []) if outs[0][0] == "f"]
CALLS[0] = 0
br = [chain_equal(a) for a in ok_chain]
exact_eq = sum(1 for a in ok_chain if a[1] == a[2])
on_thr = sum(1 for a in ok_chain if chain_equal(a) != chain_equal(a, nextup(N_TOL))
             or chain_equal(a) != chain_equal(a, nextdown(N_TOL)))
def _sub_exponent(a):
    x = -r_min(a[1], a[2]) * a[3]
    return -745.1332 < x < -708.3964
sub_rows = sum(1 for a in ok_chain if _sub_exponent(a))
check(exact_eq >= 100 and on_thr >= 4 and sub_rows >= 100,
      f"nuclear.tsv: the chain has {exact_eq} λ₁ = λ₂ rows, {on_thr} on the threshold, {sub_rows} with a subnormal e^(−λ_min·t)")
print(f"nuclear rows checked: {sum(len(v) for v in nrows.values())} ({', '.join(f'{k} {len(v)}' for k, v in nrows.items())})")
print(f"nuclear: ADR 0002 #3: {adr3} semi_empirical_binding_energy rows where the oracle's NaN (z > a or a = 0) is the "
      f"port's Err (noted), {len(semf_pairs) - adr3} finite pairs; decay_chain_two_step: {sum(1 for x in br if x)} Ok rows "
      f"on the equal-rate branch ({exact_eq} with λ₁ = λ₂ exactly), {sum(1 for x in br if x is False)} on the Bateman "
      f"branch, {on_thr} where 1e-12 ± 1 ulp moves the branch, {sub_rows} with a subnormal e^(−λ_min·t)")
# The typed formula constants (src/nuclear.cyr has them inline): each must change some golden row when perturbed.
def n_changes(fn, rule):
    return sum(1 for ins, outs, where in nrows.get(fn, []) if not toks_equal(out_toks(rule(ins)), outs))
def semf_c(k, d):
    return tuple(math.nextafter(c, math.inf if d > 0 else -math.inf) if j == k else c for j, c in enumerate(SEMF))
N_PINS = [
    ("LN_2 + 1 ulp", "nuclear::decay_constant", lambda a: r_decay_constant(*fl(a), ln2=nextup(LN2))),
    ("LN_2 − 1 ulp", "nuclear::decay_constant", lambda a: r_decay_constant(*fl(a), ln2=nextdown(LN2))),
    ("EQUAL_RATE_REL_TOL × 0.5", "nuclear::decay_chain_two_step", lambda a: r_chain(*fl(a), tol=N_TOL * 0.5)),
    ("EQUAL_RATE_REL_TOL × 2", "nuclear::decay_chain_two_step", lambda a: r_chain(*fl(a), tol=N_TOL * 2.0)),
    ("EQUAL_RATE_REL_TOL + 1 ulp", "nuclear::decay_chain_two_step", lambda a: r_chain(*fl(a), tol=nextup(N_TOL))),
    ("EQUAL_RATE_REL_TOL − 1 ulp", "nuclear::decay_chain_two_step", lambda a: r_chain(*fl(a), tol=nextdown(N_TOL))),
]
for k, name in enumerate(("VOLUME", "SURFACE", "COULOMB", "ASYMMETRY", "PAIRING")):
    for d in (1, -1):
        N_PINS.append((f"SEMF_{name} {'+' if d > 0 else '−'} 1 ulp", "nuclear::semi_empirical_binding_energy",
                       (lambda c: lambda a: r_semf(*fl(a), c=c))(semf_c(k, d))))
N_PINS += [
    ("the asymmetry term's 4.0 + 1 ulp", "nuclear::semi_empirical_binding_energy", lambda a: r_semf(*fl(a), four=nextup(4.0))),
]   # (Rust's secular 0.01 is no port constant since ADR 0002 #14: the comparison is exact, checked below)
for label, fn, rule in N_PINS:
    n_ch = n_changes(fn, rule)
    print(f"  pinning: with {label}, {n_ch} of {len(nrows.get(fn, []))} {fn.split('::')[1]} rows change")
    check(n_ch > 0, f"nuclear: no golden {fn} row pins {label}")
# ADR 0002 #13, #14 (the user's decision, 2026-11-09: fix exactly). The golden rows keep the oracle's values; a row where
# the port's value differs carries it in its note, `[ADR 0002 #N: the port returns <how> = <x> (f:<bits>)]` or `(i:<b>)]`,
# and the harness maps the two. Here every row's port value is recomputed with exact arithmetic (Fraction, MPFR for the
# correctly rounded 2^y) and the notes checked against it: present exactly where it differs from the oracle's, with its
# value. Rust's form is checked to stay wherever it is already correctly rounded.
def rn_frac(q):
    """A Fraction rounded once to binary64 (CPython's int true division is correctly rounded, subnormals and ties
    included); ±∞ past DBL_MAX."""
    try: return float(q)
    except OverflowError: return math.inf if q > 0 else -math.inf
def same_bits(x, y): return (math.isnan(x) and math.isnan(y)) or f2b(x) == f2b(y)
def port_half_lives(n0, n):
    """ADR 0002 #13: the port's N₀/2ⁿ where Rust's 2⁻ⁿ is not a normal double, with how it is obtained; None elsewhere."""
    if not math.isfinite(n) or not (n > 1022.0 or n <= -1024.0): return None
    if math.isnan(n0) or n0 == 0 or math.isinf(n0): return n0, "N₀ itself (scaleB)"
    neg = n0 < 0
    if n > 2200: return (-0.0 if neg else 0.0), "the correctly rounded N₀/2ⁿ"
    if n < -2200: return (-math.inf if neg else math.inf), "the correctly rounded N₀/2ⁿ"
    mfr, efr = math.frexp(abs(n0))
    k = math.floor(n)
    f_ = n - k
    if mfr == 0.5:
        j = efr - 1
        y = float(j) - n
        if Fraction(y) == Fraction(j) - Fraction(n):
            v = cr_of("exp2", (y,))
            return (-v if neg else v), "the correctly rounded N₀/2ⁿ (2^(j − n))"
    if f_ == 0:
        return rn_frac(Fraction(n0) * Fraction(2) ** (-k)), "the correctly rounded N₀/2ⁿ"
    # M4d: N₀·2⁻ⁿ rounded once for a fractional n too (M4c's port was N₀·CR(2⁻ᶠ)·2⁻ᵏ rounded once, faithful): MPFR,
    # 2^f to 400 bits bounded from both sides (mpfr_ref.scaled_exp2), checked below against a 60-digit 2⁻ᶠ as well
    return MPFR.scaled_exp2(n0, -n), "the correctly rounded N₀/2ⁿ (2⁻ᶠ to 320 bits)"
_dec50 = decimal.Context(prec=60)
def true_half_lives(n0, n):
    """CR(N₀·2⁻ⁿ) for a finite nonzero N₀ and finite n, |n| ≤ 2200: 2⁻ᶠ to 60 digits (decimal), the rest exact."""
    k = math.floor(n)
    f_ = n - k
    t = Fraction(_dec50.power(decimal.Decimal(2), decimal.Decimal(-f_))) if f_ else Fraction(1)
    return rn_frac(Fraction(n0) * t * Fraction(2) ** (-k))
def note_value(where, tag):
    note = NOTES.get(where, "")
    m = re.search(rf"\[{re.escape(tag)}: the port returns (.+) = (\S+) \((f|i):([0-9a-f]+)\)\]$", note)
    return m
n13 = {"N₀ itself (scaleB)": 0, "the correctly rounded N₀/2ⁿ": 0, "the correctly rounded N₀/2ⁿ (2^(j − n))": 0,
       "the correctly rounded N₀/2ⁿ (2⁻ᶠ to 320 bits)": 0}
edge_rows, frac_rows, rust_cr = 0, 0, 0   # rows past the normal range; of them, fractional n; Rust correctly rounded
for ins, outs, where in nrows.get("nuclear::atoms_after_half_lives", []):
    n0, n = fl(ins)
    oracle = outs[0][1]
    pv = port_half_lives(n0, n)
    noted = "[ADR 0002 #13:" in NOTES.get(where, "")
    if pv is None:
        check(not noted, f"{where}: an ADR 0002 #13 note where Rust's 2⁻ⁿ is a normal double")
        if math.isfinite(n0) and n0 != 0 and math.isfinite(n) and n == math.floor(n):
            check(same_bits(oracle, rn_frac(Fraction(n0) * Fraction(2) ** int(-n))),
                  f"{where}: Rust's N₀·2⁻ⁿ is not the correctly rounded N₀/2ⁿ for an integer n with 2⁻ⁿ normal")
        continue
    v, how = pv
    edge_rows += 1
    differs = not same_bits(v, oracle)
    check(noted == differs, f"{where}: atoms_after_half_lives{tuple(fl(ins))}: port {v!r}, oracle {oracle!r}: the ADR 0002 #13 "
          f"note is {'missing' if differs else 'stray'}")
    if noted:
        mm = note_value(where, "ADR 0002 #13")
        check(mm is not None and mm.group(1) == how and mm.group(3) == "f" and int(mm.group(4), 16) == f2b(v),
              f"{where}: the ADR 0002 #13 note does not name {how} = {v!r}")
        n13[how] += 1
    if math.isfinite(n0) and n0 != 0 and abs(n) <= 2200:
        tr = true_half_lives(n0, n)
        frac_rows += how.endswith("(2⁻ᶠ to 320 bits)")
        rust_cr += same_bits(oracle, tr)
        check(same_bits(v, tr), f"{where}: {how} {v!r} is not CR(N₀/2ⁿ) = {tr!r} (60-digit 2⁻ᶠ)")
        check(same_bits(v, MPFR.scaled_exp2(n0, -n)), f"{where}: {how} {v!r} is not MPFR's N₀·2⁻ⁿ rounded once")
# an exact-arithmetic sweep beyond the golden rows (seeded): MPFR's N₀·2⁻ⁿ against the 60-digit decimal 2⁻ᶠ, both rounded
# once with exact rational arithmetic (scripts/oracle/cr_verify.py runs the port's _km_cr_exp2_scaled on far more)
_rs = random.Random(0x4c13)
sw_n = 0
for _ in range(4000):
    n = _rs.uniform(1022.0, 1200.0) if _rs.random() < 0.5 else _rs.uniform(-1200.0, -1024.0)
    n0 = math.ldexp(_rs.uniform(1.0, 2.0), _rs.randint(-1074, 1023) if n > 0 else _rs.randint(-1074, 0))
    v, how = port_half_lives(n0, n)
    tr = true_half_lives(n0, n)
    sw_n += 1
    check(same_bits(v, tr), f"ADR 0002 #13 sweep: N₀ = {n0!r}, n = {n!r}: {v!r} is not CR(N₀/2ⁿ) = {tr!r}")
print(f"nuclear: ADR 0002 #13: {sum(n13.values())} atoms_after_half_lives rows noted with the port's value ("
      + ", ".join(f"{k} {c}" for k, c in n13.items()) + f") of {edge_rows} rows where 2⁻ⁿ is not a normal double "
      f"({frac_rows} with a fractional n); every one with a finite N₀ ≠ 0 is CR(N₀/2ⁿ), by MPFR and by a 60-digit 2⁻ᶠ "
      f"(Rust's value is too on {rust_cr}), and so are {sw_n} seeded pairs beyond the corpus")

def exact_secular(lp, ld):
    """ADR 0002 #14: λ_p < λ_d/100 exactly, Rust's NaN (false) and ±∞ (compared as themselves) semantics."""
    if math.isnan(lp) or math.isnan(ld): return False
    if math.isinf(lp) or math.isinf(ld): return lp < ld
    return Fraction(lp) < Fraction(ld) / 100
n14 = 0
for ins, outs, where in nrows.get("nuclear::is_secular_equilibrium", []):
    lp, ld = fl(ins)
    ex = exact_secular(lp, ld)
    differs = ex != bool(outs[0][1])
    noted = "[ADR 0002 #14:" in NOTES.get(where, "")
    check(noted == differs, f"{where}: secular({lp!r}, {ld!r}): exact {ex}, Rust {bool(outs[0][1])}: the ADR 0002 #14 "
          f"note is {'missing' if differs else 'stray'}")
    if noted:
        mm = note_value(where, "ADR 0002 #14")
        check(mm is not None and mm.group(1) == "the exact λ_p < λ_d/100" and mm.group(2) == str(ex).lower()
              and mm.group(3) == "i" and int(mm.group(4)) == int(ex), f"{where}: the ADR 0002 #14 note does not say {ex}")
        n14 += 1
# The port's comparison (src/nuclear.cyr's _km_cmp_hundred, its case analysis replayed with exact arithmetic) against
# Fraction, over a seeded sweep that crowds the boundary and the exponent shortcuts.
def cmp_hundred_replica(a, b):
    ma, ea = math.frexp(a)
    mb, eb = math.frexp(b)
    sa, ea = ma * 2, ea - 1
    sb, eb = mb * 2, eb - 1
    if eb >= ea + 8: return -1
    if eb <= ea + 5: return 1
    p = sa * 100.0
    e = Fraction(sa) * 100 - Fraction(p)
    bs = sb * (128.0 if eb == ea + 7 else 64.0)
    if p != bs: return -1 if p < bs else 1
    return (e > 0) - (e < 0)
sw14 = 0
for _ in range(20000):
    ld = math.ldexp(_rs.uniform(1.0, 2.0), _rs.randint(-1070, 1020))
    lp = ld / 100.0
    for d in (-2, -1, 0, 1, 2):
        x = lp
        for _k in range(abs(d)): x = math.nextafter(x, math.inf if d > 0 else -math.inf)
        if x <= 0: continue
        want = (Fraction(x) * 100 > Fraction(ld)) - (Fraction(x) * 100 < Fraction(ld))
        check(cmp_hundred_replica(x, ld) == want, f"ADR 0002 #14 sweep: 100·{x!r} against {ld!r}")
        sw14 += 1
print(f"nuclear: ADR 0002 #14: {n14} is_secular_equilibrium rows noted where the exact λ_p < λ_d/100 differs from "
      f"Rust's λ_p < 0.01·λ_d; the port's comparison replayed against Fraction on {sw14} seeded pairs at the boundary")
# decay_constant's LN_2 / t½ is kept (ADR 0002's policy: ln 2/t½ is irrational, a product of rounded constants).
_dc = decimal.Context(prec=80)
_ln2_dec = _dc.ln(decimal.Decimal(2))
miss_dc, n_dc = 0, 0
for ins, outs, where in nrows.get("nuclear::decay_constant", []):
    t_ = ins[0][1]
    if outs[0][0] == "f" and math.isfinite(t_) and math.isfinite(outs[0][1]) and outs[0][1] != 0:
        n_dc += 1
        if not same_f(float(_dc.divide(_ln2_dec, decimal.Decimal(t_))), outs[0][1]): miss_dc += 1
print(f"nuclear: decay_constant's LN_2/t½ is 1 ulp from the correctly rounded ln 2/t½ on {miss_dc} of {n_dc} finite "
      f"nonzero rows (kept: LN_2 is itself rounded and ln 2/t½ irrational)")

# ADR 0002 #15 (M4b's inorganic, the user's decision of 2026-11-09): Δ_tet is the correctly rounded 4Δ/9. Rust's
# (Δ·4)/9 already is wherever Δ·4 doesn't overflow (Δ·4 is exact); where it does (|Δ| ≥ 2¹⁰²²) the port takes (Δ/9)·4.
# Every golden row's port value is checked against Fraction, and the notes sit exactly where it differs from Rust's.
n15 = 0
for f_, cls, ins, outs, where in R.get("inorganic", []):
    if f_ != "inorganic::delta_tet_from_oct": continue
    x, oracle = ins[0][1], outs[0][1]
    exact = rn_frac(Fraction(4) * Fraction(x) / 9) if math.isfinite(x) and x != 0 else oracle   # ±0, ±∞, NaN: as Rust
    port = (x / 9.0) * 4.0 if (math.isinf(x * 4.0) and math.isfinite(x)) else oracle
    check(same_bits(port, exact), f"{where}: delta_tet_from_oct({x!r}): the port's {port!r} is not CR(4Δ/9) = {exact!r}")
    differs = not same_bits(exact, oracle)
    noted = "[ADR 0002 #15:" in NOTES.get(where, "")
    check(noted == differs and (not differs or abs(x) >= 2.0 ** 1022),
          f"{where}: delta_tet_from_oct({x!r}): the ADR 0002 #15 note is {'missing' if differs else 'stray'}")
    if noted:
        mm = note_value(where, "ADR 0002 #15")
        check(mm is not None and mm.group(1) == "the correctly rounded 4Δ/9" and mm.group(3) == "f"
              and int(mm.group(4), 16) == f2b(exact), f"{where}: the ADR 0002 #15 note does not name {exact!r}")
        n15 += 1
sw15 = 0
for _ in range(20000):
    big = math.ldexp(_rs.uniform(1.0, 2.0), 1022 + _rs.randint(0, 1)) * (1 if _rs.random() < 0.5 else -1)
    if math.isinf(big): continue
    check(same_bits((big / 9.0) * 4.0, rn_frac(Fraction(4) * Fraction(big) / 9)), f"ADR 0002 #15 sweep: Δ = {big!r}")
    any_ = math.ldexp(_rs.uniform(1.0, 2.0), _rs.randint(-1074, 1021)) * (1 if _rs.random() < 0.5 else -1)
    check(same_bits((any_ * 4.0) / 9.0, rn_frac(Fraction(4) * Fraction(any_) / 9)), f"ADR 0002 #15 sweep: Rust's form, Δ = {any_!r}")
    sw15 += 2
print(f"inorganic: ADR 0002 #15: {n15} delta_tet_from_oct rows noted (Rust's Δ·4 overflows; the port's (Δ/9)·4 is "
      f"CR(4Δ/9)); the port's value is CR(4Δ/9) on every golden row and {sw15} seeded ones (both branches)")

# src/nuclear_data.cyr against the Rust literal.
ntext = (root / "src/nuclear_data.cyr").read_text()
rel_n = "src/nuclear_data.cyr"
check(re.search(rf"^const KM_ISOTOPES_COUNT = {len(N_ISO)};$", ntext, re.M), f"{rel_n}: KM_ISOTOPES_COUNT")
got_s, one_each = cyr_strings_in(ntext, rel_n, "KM_ISOTOPES_SYMBOL", "ascii")
check(got_s == [r[0] for r in N_ISO] and one_each, f"{rel_n}: the symbol blob is not nuclear.rs's labels")
check(cyr_col_in(ntext, rel_n, "KM_ISOTOPES_Z", int) == [r[1] for r in N_ISO], f"{rel_n}: Z")
check(cyr_col_in(ntext, rel_n, "KM_ISOTOPES_A", int) == [r[2] for r in N_ISO], f"{rel_n}: A")
check(cyr_col_in(ntext, rel_n, "KM_ISOTOPES_HALF_LIFE_S") == [r[3] for r in N_ISO], f"{rel_n}: half-lives")
got_s, one_each = cyr_strings_in(ntext, rel_n, "KM_ISOTOPES_DECAY_MODE", "utf-8")
check(got_s == [r[4] for r in N_ISO] and one_each, f"{rel_n}: the UTF-8 decay-mode blob is not nuclear.rs's modes")
n_notes = re.findall(r"half_life_s: [0-9_.e]+, // (.+)$", nsrc, re.M)
check(len(n_notes) == 3 and all(f"half_life_s: {note}" in ntext for note in n_notes),
      f"{rel_n}: the literal's NUBASE2020 notes {n_notes} are not carried")
check(re.findall(r"^const (_?KM_\w+)", ntext, re.M) == ["KM_ISOTOPES_COUNT"], f"{rel_n}: a const other than the count")
print(f"{rel_n} checked against nuclear.rs ({len(N_ISO)} ISOTOPES rows: symbol blob, Z, A, half-life, the UTF-8 "
      f"decay-mode blob, the literal's {len(n_notes)} NUBASE2020 notes)")

# ── 16. potential.tsv (M4d): beyond the row recomputation ──────────────────────────────────────────────
# Every potential row was recomputed in the M2/M3 loop above (Rust 1.2.0's potential.rs transcribed in its operation
# order; glibc's exp, atan2 and exp2 via Python, the corrected oracle through _prim). Here: ADR 0002 #16, #17 and #18's
# notes, proven exactly; Rust's form proven exact where it stays; the typed tolerances' pinning; the grid's coverage.
prows = {}
for fn, cls, ins, outs, where in R.get("potential", []):
    prows.setdefault(fn, []).append((ins, outs, where))
def cr_root6_2(x):
    """RN(x·2^(1/6)) for a finite x ≠ 0, decided exactly: x·2^(1/6) > t ⟺ 2·x⁶ > t⁶ (t > 0), so each candidate's two
    midpoints are compared in rational arithmetic and the candidate moved until they bracket the value (never a
    midpoint: 2^(1/6) is irrational). Past DBL_MAX + ulp/2 = 2^1024 − 2^970, ±∞."""
    two_x6 = 2 * abs(Fraction(x)) ** 6
    v = abs(x) * 1.122462048309373
    if math.isinf(v): v = sys.float_info.max
    while True:
        up, down = nextup(v), nextdown(v)
        hi = (Fraction(v) + Fraction(up)) / 2 if math.isfinite(up) else Fraction(2) ** 1024 - Fraction(2) ** 970
        lo = (Fraction(v) + Fraction(down)) / 2
        if two_x6 > hi ** 6:
            if not math.isfinite(up): return math.copysign(math.inf, x)
            v = up
        elif two_x6 < lo ** 6:
            v = down
        else:
            return math.copysign(v, x)
def note_val(where, tag, how):
    m = re.search(rf"\[{re.escape(tag)}: the port returns {re.escape(how)} = (\S+) \(f:([0-9a-f]+)\)\]$", NOTES.get(where, ""))
    return None if m is None else b2f(int(m.group(2), 16))
# #16: σ·2^(1/6) rounded once, against Rust's σ·RN(2^(1/6)).
P16_HOW = "the correctly rounded σ·2^(1/6)"
n16, n16_rows, rust16 = 0, 0, 0
for ins, outs, where in prows.get("potential::lj_equilibrium_distance", []):
    sg = ins[0][1]
    oracle = outs[0][1]
    noted = "[ADR 0002 #16:" in NOTES.get(where, "")
    if not math.isfinite(sg) or sg == 0:
        check(not noted and same_bits(oracle, sg), f"{where}: lj_equilibrium_distance({sg!r}) is not σ itself, or noted")
        continue
    n16_rows += 1
    v = cr_root6_2(sg)
    differs = not same_bits(v, oracle)
    rust16 += not differs
    check(noted == differs, f"{where}: lj_equilibrium_distance({sg!r}): port {v!r}, oracle {oracle!r}: the #16 note is "
          f"{'missing' if differs else 'stray'}")
    if noted:
        nv = note_val(where, "ADR 0002 #16", P16_HOW)
        check(nv is not None and same_bits(nv, v), f"{where}: the ADR 0002 #16 note does not name {v!r}")
        n16 += 1
check(f2b(MPFR.cr("exp2", 1.0 / 6.0)) == 0x3FF1F59AC3C7D6C0 == f2b(1.122462048309373),
      "Rust's 2.0.powf(1.0 / 6.0) = CR(2^RN(1/6)) is not RN(2^(1/6)) = 0x3FF1F59AC3C7D6C0")
# #17 and #18 (M4d's (ε·d)·4 and (μ·½)/D_e past 4·ε's and 2·D_e's overflow) are folded into ADR 0002 #32 (the overflow
# audit's B3): every lennard_jones and morse_width_parameter row is proven in §21, #17's rows bit for bit there. What
# stays here is M4d's arithmetic fact behind them: Rust's (4·ε)·d is RN(4·ε·d) wherever 4·ε is finite and d is Rust's.
n17 = 0
for ins, outs, where in prows.get("potential::lennard_jones", []):
    e_, sg, r = fl(ins)
    if outs[0][0] != "f" or e_ >= 2.0 ** 1022: continue
    sr = dv(sg, r)
    sr6 = sr * sr * sr * sr * sr * sr
    d = sr6 * sr6 - sr6
    if math.isfinite(d) and math.isfinite(e_):
        check(same_bits(outs[0][1], rn_frac(4 * Fraction(e_) * Fraction(d))), f"{where}: Rust's (4·ε)·d is not RN(4·ε·d)")
        n17 += 1
# seeded exact sweeps beyond the corpus: #16's σ over every binade
_rp = random.Random(0x9016)
sw16 = 0
for _ in range(20000):
    x = b2f(_rp.randrange(1, 0x7FEFFFFFFFFFFFFF))
    v = cr_root6_2(x)
    rust = x * 1.122462048309373
    sw16 += not same_bits(v, rust)
    check(abs(ulp_dist(f2b(v), f2b(rust))) <= 1, f"#16 sweep: RN(σ·2^(1/6)) {v!r} is not within an ulp of Rust's {rust!r}")
print(f"potential rows checked: {sum(len(v) for v in prows.values())} ({', '.join(f'{k} {len(v)}' for k, v in prows.items())})")
print(f"potential: ADR 0002 #16: {n16} lj_equilibrium_distance rows noted with the correctly rounded σ·2^(1/6) (of "
      f"{n16_rows} finite nonzero σ; Rust's σ·RN(2^(1/6)) is correctly rounded on {rust16}), proven exactly (2·σ⁶ against "
      f"each midpoint's sixth power), and Rust's value differs on {sw16} of 20,000 seeded σ (always by one ulp); Rust's "
      f"(4·ε)·d is RN(4·ε·d) on all {n17} finite lennard_jones rows with ε < 2¹⁰²² (#17 and #18 are folded into #32: §21)")
# ADR 0002 #20 (M4e-2): distance_3d, bond_angle and dihedral_angle are Rust's sequence evaluated in binary64 with an
# unbounded exponent range (every operation rounded to 53 bits, nothing over- or underflowing), the final √ or atan2
# rounded once into binary64. The evaluator is MPFR's (mpfr_ref.UX: 53 bits at MPFR's widest exponent range; the final
# operation by mpfr_check_range and mpfr_subnormalize), independent of the port's doubles and of the generator's
# integers, with Rust's sequences transcribed in mpfr_ref.potential_* (sub3 / cross / dot3 in Rust's order; the
# tolerances are the ones parsed from the source above). On every row with finite coordinates: a noted row's port
# result is the evaluator's and differs from the oracle's; an unnoted row's expected result (the corrected oracle's
# where there is one, else the oracle's) is the evaluator's. A row whose f64 chain stays in the normal range (no product
# or sum non-finite or below 2^-1022, no zero product of nonzero factors) is never noted: there Rust's bits are the
# port's.
def _e20(r): return Err("InvalidInput") if r is None else r
def e20_distance(v): return MPFR.potential_distance_3d(v)
def e20_bond(v): return _e20(MPFR.potential_bond_angle(v, P_TOL))
def e20_dihedral(v): return _e20(MPFR.potential_dihedral_angle(v, P_CTOL))
# Rust's f64 chains, traced: OOR[0] is set when a product, quotient or sum leaves the normal range: a result that is not
# finite, 0 from nonzero factors, nonzero and below 2⁻¹⁰²² in magnitude, or a product or quotient of exactly ±2⁻¹⁰²² that
# is not the 53-bit rounding of its exact value (one just below 2⁻¹⁰²², which IEEE rounds up there at the subnormal grid
# and 53-bit rounding to 2⁻¹⁰²² − 2⁻¹⁰⁷⁵: B4's correction of ADR 0002 #29, whose B2 tracer b2_nernst_in_range took every
# 2⁻¹⁰²² as normal, as the port's fast path did; #20's and B1's tracers share these). A sum of two doubles is exact
# wherever it is below 2⁻¹⁰²¹ in magnitude, so ±2⁻¹⁰²² is in range for a sum.
OOR = [False]
N1022 = 2.0 ** -1022
def _t_chk(r, *fs, rn53=None):
    if not math.isfinite(r) or (r != 0.0 and abs(r) < N1022) or (r == 0.0 and fs and all(x != 0.0 for x in fs)):
        OOR[0] = True
    elif abs(r) == N1022 and rn53 is not None and not rn53().is_f64(r):
        OOR[0] = True
    return r
def tm(x, y): return _t_chk(x * y, x, y, rn53=lambda: MPFR.UX(x) * MPFR.UX(y))
def tq(x, y): return _t_chk(dv(x, y), x, rn53=lambda: MPFR.UX(x) / MPFR.UX(y))
def ta(x, y): return _t_chk(x + y)
def ts(x, y): return _t_chk(x - y)
def t_sub(a, b): return [ts(a[0], b[0]), ts(a[1], b[1]), ts(a[2], b[2])]
def t_dot(a, b): return ta(ta(tm(a[0], b[0]), tm(a[1], b[1])), tm(a[2], b[2]))
def t_cross(a, b):
    return [ts(tm(a[1], b[2]), tm(a[2], b[1])), ts(tm(a[2], b[0]), tm(a[0], b[2])), ts(tm(a[0], b[1]), tm(a[1], b[0]))]
def t_sqrt(x): return _t_chk(r_sqrt(x)) if x != 0.0 else x
def in_range_distance(v):
    OOR[0] = False
    d = t_sub(v[3:6], v[0:3])
    t_sqrt(t_dot(d, d))
    return not OOR[0]
def in_range_bond(v):
    OOR[0] = False
    a, b, c = v[0:3], v[3:6], v[6:9]
    ba, bc = t_sub(a, b), t_sub(c, b)
    ba_sq, bc_sq = t_dot(ba, ba), t_dot(bc, bc)
    tol_sq = tm(P_TOL, P_TOL)
    b_sq = t_dot(b, b)
    if (not is_positive(ba_sq) or not is_positive(bc_sq) or ba_sq <= tm(tol_sq, r_max(t_dot(a, a), b_sq))
            or bc_sq <= tm(tol_sq, r_max(t_dot(c, c), b_sq))):
        return not OOR[0]
    cr = t_cross(ba, bc)
    t_sqrt(t_dot(cr, cr))
    t_dot(ba, bc)
    return not OOR[0]
def in_range_dihedral(v):
    OOR[0] = False
    a, b, c, d = v[0:3], v[3:6], v[6:9], v[9:12]
    b1, b2, b3 = t_sub(b, a), t_sub(c, b), t_sub(d, c)
    n1, n2 = t_cross(b1, b2), t_cross(b2, b3)
    n1_sq, n2_sq, b2_sq = t_dot(n1, n1), t_dot(n2, n2), t_dot(b2, b2)
    if (not is_positive(n1_sq) or not is_positive(n2_sq) or n1_sq <= tm(tm(P_CTOL, t_dot(b1, b1)), b2_sq)
            or n2_sq <= tm(tm(P_CTOL, b2_sq), t_dot(b3, b3))):
        return not OOR[0]
    tm(t_sqrt(b2_sq), t_dot(b1, n2))
    t_dot(n1, n2)
    return not OOR[0]
P20_RE = re.compile(r"\[ADR 0002 #20: the port returns Rust's sequence evaluated with an unbounded exponent(?: \((?:√|atan2) "
                    r"rounded once\))? = (?:(Ok\()?(\S+?)\)? \(f:([0-9a-f]{16})\)|Err\(InvalidInput\) \(E:InvalidInput\))\]$")
CORR_POT = {w: ct for w, ct, _ in CORRECTIONS.get("potential", [])}
E20 = {"potential::distance_3d": (e20_distance, in_range_distance, False),
       "potential::bond_angle": (e20_bond, in_range_bond, True),
       "potential::dihedral_angle": (e20_dihedral, in_range_dihedral, True)}
n20 = {fn: [0, 0, 0, 0] for fn in E20}     # rows with finite coordinates, in range, noted, noted Err
e20_vals = {}
for fn, (ev, inr, isres) in E20.items():
    for ins, outs, where in prows.get(fn, []):
        v = fl(ins)
        note = NOTES.get(where, "")
        noted = "[ADR 0002 #20:" in note
        if not all(math.isfinite(x) for x in v):
            check(not noted, f"{where}: an ADR 0002 #20 note on a row with a non-finite coordinate (Rust's result stays)")
            continue
        n20[fn][0] += 1
        with MPFR.unbounded():
            u = ev(v)
        e20_vals[where] = u
        utoks = out_toks(u)
        expect = CORR_POT.get(where, outs)
        in_range = inr(v)
        n20[fn][1] += in_range
        if in_range:
            check(not noted and toks_equal(utoks, expect), f"{where}: {fn}{tuple(v)}: Rust's chain stays in the normal "
                  f"range, so the unbounded evaluation {u!r} must be the expected result {expect}, unnoted")
        if not noted:
            check(toks_equal(utoks, expect), f"{where}: {fn}{tuple(v)}: the unbounded evaluation {u!r} differs from the "
                  f"expected {expect}, and the row has no ADR 0002 #20 note")
            continue
        mm = P20_RE.search(note)
        check(mm is not None, f"{where}: malformed ADR 0002 #20 note")
        if mm is None: continue
        if mm.group(3) is None:
            check(isinstance(u, Err), f"{where}: the #20 note says Err, the unbounded evaluation gives {u!r}")
            n20[fn][3] += 1
        else:
            nv = b2f(int(mm.group(3), 16))
            check(not isinstance(u, Err) and same_bits(nv, u), f"{where}: the #20 note names {nv!r}, the unbounded "
                  f"evaluation gives {u!r}")
            check((mm.group(1) is not None) == isres and same_bits(float(mm.group(2)), nv), f"{where}: the #20 note's "
                  f"value text")
        check(not toks_equal(utoks, outs), f"{where}: an ADR 0002 #20 note where the port's result is the oracle's")
        n20[fn][2] += 1
# Scale invariance through the reference (the property the port's tests check bit for bit): each angle row with
# finite coordinates, scaled by 2^k for three k where every scaled coordinate is exact, gives the same result; a
# distance gives 2^k times it where both are normal.
_r20 = random.Random(0x9020)
inv20 = 0
def exact_scale(v, k):
    try:
        w = [math.ldexp(x, k) for x in v]
    except OverflowError:
        return None
    return w if all(math.isfinite(y) and math.ldexp(y, -k) == x for x, y in zip(v, w)) else None
for fn, (ev, inr, isres) in E20.items():
    for ins, outs, where in prows.get(fn, []):
        if where not in e20_vals or _r20.random() > 0.25: continue
        v, u = fl(ins), e20_vals[where]
        for _ in range(3):
            k = _r20.randint(-1100, 1100)
            w = exact_scale(v, k)
            if w is None: continue
            with MPFR.unbounded():
                us = ev(w)
            if isres:
                check(toks_equal(out_toks(us), out_toks(u)), f"{where}: {fn} is not scale-invariant at 2^{k}: {us!r}, {u!r}")
            elif (u == 0.0 or 2.0 ** -1022 <= u <= sys.float_info.max) and (us == 0.0 or 2.0 ** -1022 <= us < math.inf):
                check(same_bits(us, math.ldexp(u, k)), f"{where}: distance_3d at 2^{k} is {us!r}, not 2^{k}·{u!r}")
            inv20 += 1
for fn, (a_, b_, c_, d_) in n20.items():
    print(f"potential: ADR 0002 #20: {fn.split('::')[1]}: {a_} rows with finite coordinates, {b_} with Rust's chain in "
          f"the normal range (all unnoted, Rust's bits); {c_} noted with the unbounded evaluation ({d_} of them Err), each "
          f"proven against MPFR, and every other row's expected result is MPFR's")
print(f"potential: ADR 0002 #20: scale invariance through MPFR on {inv20} scaled rows")
check(all(n20[fn][2] > 0 for fn in n20), "potential: ADR 0002 #20 has no noted row for some function")
# The typed tolerances (src/potential.cyr has them inline): each must change some golden row when perturbed.
def p_changes(fn, rule):
    return sum(1 for ins, outs, where in prows.get(fn, []) if not toks_equal(out_toks(rule(ins)), outs))
P_PINS = [
    ("COINCIDENT_REL_TOL + 1 ulp", "potential::bond_angle", lambda a: r_bond_angle(fl(a), tol=nextup(P_TOL))),
    ("COINCIDENT_REL_TOL − 1 ulp", "potential::bond_angle", lambda a: r_bond_angle(fl(a), tol=nextdown(P_TOL))),
    ("COINCIDENT_REL_TOL × 0.5", "potential::bond_angle", lambda a: r_bond_angle(fl(a), tol=P_TOL * 0.5)),
    ("COINCIDENT_REL_TOL × 2", "potential::bond_angle", lambda a: r_bond_angle(fl(a), tol=P_TOL * 2.0)),
    ("COLLINEAR_SIN_SQ_TOL + 1 ulp", "potential::dihedral_angle", lambda a: r_dihedral(fl(a), tol=nextup(P_CTOL))),
    ("COLLINEAR_SIN_SQ_TOL − 1 ulp", "potential::dihedral_angle", lambda a: r_dihedral(fl(a), tol=nextdown(P_CTOL))),
    ("COLLINEAR_SIN_SQ_TOL × 0.5", "potential::dihedral_angle", lambda a: r_dihedral(fl(a), tol=P_CTOL * 0.5)),
    ("COLLINEAR_SIN_SQ_TOL × 2", "potential::dihedral_angle", lambda a: r_dihedral(fl(a), tol=P_CTOL * 2.0)),
    ("the −π → π mapping's π − 1 ulp", "potential::dihedral_angle", lambda a: r_dihedral(fl(a), pi=nextdown(math.pi))),
    ("Lennard-Jones' 4.0 + 1 ulp", "potential::lennard_jones", lambda a: r_lj(*fl(a), four=nextup(4.0))),
    ("the width's 2.0 + 1 ulp", "potential::morse_width_parameter", lambda a: r_width(*fl(a), two=nextup(2.0))),
]
for label, fn, rule in P_PINS:
    CALLS[0] = 0
    n_ch = p_changes(fn, rule)
    print(f"  pinning: with {label}, {n_ch} of {len(prows.get(fn, []))} {fn.split('::')[1]} rows change")
    check(n_ch > 0, f"potential: no golden {fn} row pins {label}")
POT_MAPPED[0] = 0
for ins, outs, where in prows.get("potential::dihedral_angle", []):
    r_dihedral(fl(ins))
mapped_pi = POT_MAPPED[0]
neg_zero = sum(1 for ins, outs, w in prows.get("potential::dihedral_angle", []) if outs[0][0] == "f" and outs[0][2] == 0x8000000000000000)
check(mapped_pi >= 8, f"potential: only {mapped_pi} dihedral rows take Rust's −π → π mapping")
for fn, need in (("potential::bond_angle", [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]]),
                 ("potential::dihedral_angle", [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, -1.0, 0.0, 3.0, -1.0, 0.0]]),
                 ("potential::distance_3d", [[0.0, 0.0, 0.0, 3.0, 4.0, 0.0]])):
    have = {tuple(fl(ins)) for ins, _, _ in prows.get(fn, [])}
    check(all(tuple(v) in have for v in need), f"potential: {fn} misses a Rust test's geometry")
print(f"potential: dihedral_angle: {mapped_pi} rows where atan2 gives −π and Rust reports π, {neg_zero} rows of −0.0")

# ── 17. spectroscopy.tsv (M4e): beyond the row recomputation ──────────────────────────────────────────
# Every spectroscopy row was recomputed in the M2/M3 loop above (Rust 1.2.0's spectroscopy.rs transcribed in its operation
# order, with hisab 1.4.0's fft for fft_spectrum; glibc's pow(10, ·), log10, cos and sin via Python, the corrected oracle
# through _prim). Here: ADR 0002 #19's notes, proven exactly; the generated data against the literal; the typed constants'
# pinning; the stable sort's pinning; the grids' coverage; and ADR 0002 #21's and #22's notes (M4e-3), proven exactly,
# with the structure of their exact values and ties.
srows = {}
for fn, cls, ins, outs, where in R.get("spectroscopy", []):
    srows.setdefault(fn, []).append((ins, outs, where))
# ADR 0002 #19: RN(1/(100·x)) against Rust's RN(RN(0.01)/x).
def note_val19(where, var):
    m = re.search(rf"\[ADR 0002 #19: the port returns the correctly rounded 1/\(100·{var}\) = (\S+) \(f:([0-9a-f]+)\)\]$",
                  NOTES.get(where, ""))
    return None if m is None else b2f(int(m.group(2), 16))
n19, ok19 = {}, {}
for fn, var in (("spectroscopy::wavenumber_from_wavelength", "λ"), ("spectroscopy::wavelength_from_wavenumber", "ν̃")):
    n19[fn] = ok19[fn] = 0
    for ins, outs, where in srows.get(fn, []):
        x = ins[0][1]
        noted = "[ADR 0002 #19:" in NOTES.get(where, "")
        if outs[0][0] != "f":
            check(not noted and not is_positive(x), f"{where}: {fn}({x!r}): an Err row on a positive finite input, or noted")
            continue
        ok19[fn] += 1
        exact = rn_frac(Fraction(1, 100) / Fraction(x))
        differs = not same_bits(exact, outs[0][1])
        check(noted == differs, f"{where}: {fn}({x!r}): exact {exact!r}, oracle {outs[0][1]!r}: the #19 note is "
              f"{'missing' if differs else 'stray'}")
        if noted:
            nv = note_val19(where, var)
            check(nv is not None and same_bits(nv, exact), f"{where}: the ADR 0002 #19 note does not name {exact!r}")
            n19[fn] += 1
            check(ulp_dist(f2b(exact), outs[0][2]) <= 1 or not math.isfinite(exact) or exact == 0.0 or abs(exact) < 2.0 ** -1022,
                  f"{where}: Rust's 0.01/x is more than one ulp from RN(1/(100x)) on a normal result")
_r19 = random.Random(0x9019)
sw19 = 0
for _ in range(20000):
    x = b2f(_r19.randrange(1, 0x7FEFFFFFFFFFFFFF))
    sw19 += not same_bits(rn_frac(Fraction(1, 100) / Fraction(x)), dv(0.01, x))
check(f2b(0.01) == 0x3F847AE147AE147B and Fraction(0.01) != Fraction(1, 100), "0.01 is not RN(1/100), or is exact")
print(f"spectroscopy: ADR 0002 #19: {n19['spectroscopy::wavenumber_from_wavelength']} of "
      f"{ok19['spectroscopy::wavenumber_from_wavelength']} wavenumber_from_wavelength and "
      f"{n19['spectroscopy::wavelength_from_wavenumber']} of {ok19['spectroscopy::wavelength_from_wavenumber']} "
      f"wavelength_from_wavenumber Ok rows noted with RN(1/(100·x)) (proven with Fraction); Rust's 0.01/x differs from it "
      f"on {sw19} of 20,000 seeded x over every binade")
# The generated data (src/spectroscopy_data.cyr) against the literal and the API.
sdata = (root / "src/spectroscopy_data.cyr").read_text()
def s_col(name, kind=float):
    m = re.search(rf"^var {name}: \w+\[(\d+)\] = \{{(.*?)\}};", sdata, re.M | re.S)
    check(m is not None, f"spectroscopy_data.cyr: no {name}")
    if m is None: return []
    vals = [kind(x) for line in m.group(2).splitlines() for x in re.findall(r"(-?[0-9.]+),", line.split("#", 1)[0])]
    check(len(vals) == int(m.group(1)), f"spectroscopy_data.cyr: {name} declares {m.group(1)}, holds {len(vals)}")
    return vals
s_members = re.findall(r"^    (KM_SPECTRAL_SERIES_[A-Z_]+) = (\d+);", sdata, re.M)
check(s_members == [(f"KM_SPECTRAL_SERIES_{snake_t(v)}", str(k)) for k, v in enumerate(S_SERIES)],
      f"spectroscopy_data.cyr: KmSpectralSeries {s_members} is not SpectralSeries {S_SERIES} in order")
check(re.search(rf"^const KM_SPECTRAL_SERIES_COUNT = {len(S_SERIES)};$", sdata, re.M), "spectroscopy_data.cyr: the count")
check(s_col("KM_SPECTRAL_SERIES_N_LOWER", int) == [S_N_LOWER[v] for v in S_SERIES], "spectroscopy_data.cyr: n_lower")
# ADR 0002 #23 (M4e-4): the isotope table is scripts/oracle/isotopes.tsv's (parsed above), checked here against the
# oracle's: (a) the same elements and isotopes as Rust's literal (an isotope's A is its mass rounded), plus ³⁶S as the only
# addition; (b) every shared abundance is Rust's literal as a decimal string; (c) every Rust mass is the AME2020 mass
# rounded to five decimals (Fraction; a tie would be ambiguous and fails); (d) each element's abundances sum to exactly 1.
iso_rust = [(z, round(m_)) for z in sorted(S_ISO) for m_, a_ in S_ISO[z]]
check(all(abs(m_ - round(m_)) < 0.1 for z in S_ISO for m_, a_ in S_ISO[z]), "spectroscopy.rs: an isotope mass names no A")
iso_new = [(r[0], r[1]) for r in ISO23]
check(iso_new.count((16, 36)) == 1 and [za for za in iso_new if za != (16, 36)] == iso_rust,
      f"(a) isotopes.tsv's isotopes {iso_new} are not Rust's {iso_rust} plus ³⁶S")
shared23 = [r for r in ISO23 if (r[0], r[1]) != (16, 36)]
lits23 = [pr for z in sorted(S_ISO_LIT) for pr in S_ISO_LIT[z]]
for (z, a, sym, mass, ab), (lm, la) in zip(shared23, lits23):
    check(ab == la, f"(b) isotopes.tsv: {a}{sym}: abundance {ab} is not Rust's literal {la}")
    q = Fraction(mass) * 10 ** 5
    check(q - math.floor(q) != Fraction(1, 2) and Fraction(round(q), 10 ** 5) == Fraction(lm) and len(lm.split(".")[1]) <= 5,
          f"(c) isotopes.tsv: {a}{sym}: Rust's mass {lm} is not AME2020's {mass} to five decimals")
for z in sorted(S_ISO_NEW):
    check(sum(Fraction(r[4]) for r in ISO23 if r[0] == z) == 1, f"(d) isotopes.tsv: Z = {z}'s abundances do not sum to 1")
check(sum(Fraction(la) for la in (a_ for m_, a_ in S_ISO_LIT[16])) == Fraction(9999, 10000),
      "spectroscopy.rs: Rust's sulfur no longer sums to 0.9999 (ADR 0002 #23's premise)")
# CIAAW's standard atomic weight intervals (the `interval` records): each element's abundance-weighted mean mass inside.
iv_text = []
for sym, lo, hi in ISO23_IV:
    recs = [r for r in ISO23 if r[2] == sym]
    mean = sum(Fraction(r[3]) * Fraction(r[4]) for r in recs)
    check(len(recs) > 0 and lo <= mean <= hi, f"isotopes.tsv: {sym}'s mean mass {float(mean)!r} is outside [{lo}, {hi}]")
    iv_text.append(f"{sym} {float(mean):.6f} in [{float(lo)}, {float(hi)}]")
print(f"spectroscopy: ADR 0002 #23: isotopes.tsv's {len(ISO23)} isotopes are Rust's {len(iso_rust)} plus ³⁶S, each mass "
      f"its verbatim AME2020 line's, the shared abundances Rust's literals, Rust's masses the AME2020 ones to five "
      f"decimals, each element summing to 1 exactly (Rust's sulfur 0.9999); mean masses {', '.join(iv_text)} (CIAAW)")
zs_ = max(S_ISO_NEW) + 1
check(zs_ == max(S_ISO) + 1, "isotopes.tsv: its last element is not Rust's")
check(re.search(rf"^const _KM_ISOTOPE_ABUNDANCE_ZS = {zs_};$", sdata, re.M), "spectroscopy_data.cyr: _KM_ISOTOPE_ABUNDANCE_ZS")
flat = [x for z in sorted(S_ISO_NEW) for pr in S_ISO_NEW[z] for x in pr]
check(s_col("_KM_ISOTOPE_ABUNDANCES") == flat, "spectroscopy_data.cyr: the abundance pairs are not isotopes.tsv's (ADR 0002 #23)")
check(all(f2b(x) == f2b(y) for x, y in zip(s_col("_KM_ISOTOPE_ABUNDANCES"), flat)), "spectroscopy_data.cyr: a pair's bits")
starts, counts, at = [], [], 0
for z in range(zs_):
    starts.append(at)
    counts.append(len(S_ISO_NEW.get(z, [])))
    at += counts[-1]
check(s_col("_KM_ISOTOPE_ABUNDANCE_START", int) == starts and s_col("_KM_ISOTOPE_ABUNDANCE_COUNT", int) == counts,
      "spectroscopy_data.cyr: the abundance start / count columns")
def s_const(name):
    m = re.search(rf"^const {name} = ([0-9.]+);$", sdata, re.M)
    check(m is not None and re.fullmatch(r"[0-9]+\.[0-9]+", m.group(1)), f"spectroscopy_data.cyr: {name}")
    return float(m.group(1)) if m else math.nan
check(same_bits(s_const("_KM_ELECTRON_PROTON_MASS_RATIO"), S_RATIO), "spectroscopy_data.cyr: the mass ratio literal")
check(same_bits(s_const("_KM_RYDBERG_HYDROGEN"), S_RYD_H), "spectroscopy_data.cyr: R_H is not RYDBERG / (1 + ratio)")
lyman_top = [outs for ins, outs, w in srows.get("SpectralSeries::wavelength", [])
             if ins[0][1] == "SpectralSeries.Lyman" and ins[1][1] == 0xFFFFFFFE]
check(len(lyman_top) == 1 and same_bits(lyman_top[0][0][1], 1.0 / S_RYD_H)
      and sum(1 for x in (nextdown(S_RYD_H), S_RYD_H, nextup(S_RYD_H)) if same_bits(1.0 / x, 1.0 / S_RYD_H)) == 1,
      "spectroscopy: the Lyman line u32::MAX − 1 does not pin R_H uniquely")
gz = {ins[0][1]: outs for ins, outs, w in srows.get("spectroscopy::isotope_abundances", [])}
check(sorted(gz) == list(range(256)), "spectroscopy.tsv: isotope_abundances is not probed at every u8")
print(f"spectroscopy_data.cyr checked: {len(S_SERIES)} series, {len(flat) // 2} isotope pairs over {len(S_ISO_NEW)} elements "
      f"(isotopes.tsv's; every u8 probed), R_H = {S_RYD_H!r} = RYDBERG / (1 + {S_RATIO!r}), unique through the Lyman line")
# ADR 0002 #23's notes: on every isotope_abundances and isotope_pattern row, the port's result from S_ISO_NEW (the
# replica r_pattern with the table swapped; the rows' own check above proves the replica gives the oracle's bits with
# Rust's table) is computed; a row is noted exactly where it differs from the oracle's outputs, and the note's tokens are
# it, bit for bit. So every unnoted row's port result is the oracle's too.
P23_RE = re.compile(r"\[ADR 0002 #23: with the AME2020 masses and the TICE 2009 abundances, ³⁶S included, the port returns "
                    r"(\[\d+\]|Err\(\w+\))(?:, [^()]*)? \(([^()]*)\)\]$")
def port23(fn, ins):
    if fn == "spectroscopy::isotope_abundances": return Pairs(S_ISO_NEW.get(ins[0][1], []))
    return r_pattern(pattern_atoms(ins), table=S_ISO_NEW)
def expect23(fn, ins, outs, where):
    """The port's expected tokens for an isotope row: the #23 note's where there is one, else the oracle's outputs."""
    m = P23_RE.search(NOTES.get(where, ""))
    return outs if m is None else [parse_tok(t, where) for t in m.group(2).split(" ")]
n23, ok23 = {}, {}
for fn in ("spectroscopy::isotope_abundances", "spectroscopy::isotope_pattern"):
    n23[fn] = ok23[fn] = 0
    for ins, outs, where in srows.get(fn, []):
        note = NOTES.get(where, "")
        noted = "[ADR 0002 #23:" in note
        port = out_toks(port23(fn, ins))
        differs = not toks_equal(port, outs)
        check(noted == differs, f"{where}: {fn}: the #23 note is {'missing' if differs else 'stray'}")
        ok23[fn] += outs[0][0] == "v"
        if noted:
            m = P23_RE.search(note)
            check(m is not None and toks_equal(expect23(fn, ins, outs, where), port)
                  and m.group(1) == (f"[{port[0][1]}]" if port[0][0] == "v" else f"Err({port[0][1]})"),
                  f"{where}: the ADR 0002 #23 note is not the port's result from isotopes.tsv")
            check(port[0][0] == "v" and outs[0][0] == "v", f"{where}: a #23 note on an Err row")
            n23[fn] += 1
print(f"spectroscopy: ADR 0002 #23: {n23['spectroscopy::isotope_abundances']} isotope_abundances rows (of 256) and "
      f"{n23['spectroscopy::isotope_pattern']} isotope_pattern rows (of {ok23['spectroscopy::isotope_pattern']} Ok, "
      f"{len(srows.get('spectroscopy::isotope_pattern', []))} in all) noted with the port's result from isotopes.tsv "
      f"(r_pattern with the table swapped, bit for bit); every other row's port result is the oracle's")
# The typed constants (inline in src/spectroscopy.cyr) and the stable sort: each must change some golden row when
# perturbed, where a row can reach it (isotope rows with at most 3,000 convolution steps, for speed).
small_iso = [(ins, outs, w) for ins, outs, w in srows.get("spectroscopy::isotope_pattern", [])
             if sum(c for z, c in pattern_atoms(ins) if z in S_ISO) <= 3000]
def s_changes(rows_, rule):
    CALLS[0] = 0
    return sum(1 for ins, outs, where in rows_ if not toks_equal(out_toks(rule(ins)), outs))
S_PINS = [
    ("ISOTOPE_PRUNE_REL × 0.5", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), prune=S_PRUNE * 0.5, table=t), True),
    ("ISOTOPE_PRUNE_REL × 2", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), prune=S_PRUNE * 2.0, table=t), True),
    ("ISOTOPE_PRUNE_REL + 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), prune=nextup(S_PRUNE), table=t), False),
    ("ISOTOPE_PRUNE_REL − 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), prune=nextdown(S_PRUNE), table=t), False),
    ("ISOTOPE_REPORT_REL × 0.5", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), report=S_REPORT * 0.5, table=t), True),
    ("ISOTOPE_REPORT_REL × 2", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), report=S_REPORT * 2.0, table=t), True),
    ("ISOTOPE_REPORT_REL + 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), report=nextup(S_REPORT), table=t), False),
    ("ISOTOPE_REPORT_REL − 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), report=nextdown(S_REPORT), table=t), False),
    ("ISOTOPE_BIN_WIDTH_U × 2", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), width=S_BIN * 2.0, table=t), True),
    ("ISOTOPE_BIN_WIDTH_U × 0.5", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), width=S_BIN * 0.5, table=t), False),
    ("ISOTOPE_BIN_WIDTH_U + 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), width=nextup(S_BIN), table=t), False),
    ("ISOTOPE_BIN_WIDTH_U − 1 ulp", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), width=nextdown(S_BIN), table=t), False),
    ("the stable sort (ties reversed)", "iso", lambda a, t=None: r_pattern(pattern_atoms(a), ties_reversed=True, table=t), False),
    ("ħ's 2.0 + 1 ulp", "spectroscopy::harmonic_oscillator_energy", lambda a: r_oscillator(*fl(a), two=nextup(2.0)), True),
    ("the oscillator's 0.5 + 1 ulp", "spectroscopy::harmonic_oscillator_energy",
     lambda a: r_oscillator(*fl(a), half=nextup(0.5)), True),
    ("B's 8.0 + 1 ulp", "spectroscopy::rotational_constant", lambda a: r_rot_const(fl(a)[0], eight=nextup(8.0)), True),
    ("the rotor's 1.0 + 1 ulp", "spectroscopy::rigid_rotor_energy", lambda a: r_rotor(*fl(a), one=nextup(1.0)), None),
    ("R_H + 1 ulp", "SpectralSeries::wavelength", lambda a: r_series(*fl(a), rh=nextup(S_RYD_H)), True),
    ("R_H − 1 ulp", "SpectralSeries::wavelength", lambda a: r_series(*fl(a), rh=nextdown(S_RYD_H)), True),
    ("the FFT's −2.0 − 1 ulp", "spectroscopy::fft_spectrum", lambda a: r_fft(fl(a)[1:], m2=nextdown(-2.0)), True),
    ("find_peaks' 1.0 − 1 ulp", "spectroscopy::find_peaks", lambda a: r_find_peaks(*peaks_args(a), one=nextdown(1.0)),
     False),
]
for label, fn, rule, must in S_PINS:
    rows_ = small_iso if fn == "iso" else srows.get(fn, [])
    n_ch = s_changes(rows_, rule)
    what = "isotope_pattern rows (≤ 3,000 steps)" if fn == "iso" else f"{fn.split('::')[1]} rows"
    why = ""
    if must is None and n_ch == 0:
        why = " (unobservable on every input: j + 1 + 2⁻⁵² rounds back to j + 1 for j ≥ 1, and j = 0 multiplies by 0)"
    elif not must and n_ch == 0:
        why = " (no row can reach it: tests/spectroscopy.tcyr pins it on the private helper)"
    print(f"  pinning: with {label}, {n_ch} of {len(rows_)} {what} change{why}")
    if must:
        check(n_ch > 0, f"spectroscopy: no golden row pins {label}")
# ADR 0002 #23: the harness compares a noted isotope row with its note (the port's result from isotopes.tsv), so the
# port's typed constants must also be pinned there: each isotope perturbation above (its rule's table swapped for the
# #23 one) must change some row's expected port result (the note's tokens, or the oracle's outputs on an unnoted row).
for label, fn, rule, must in S_PINS:
    if fn != "iso": continue
    n_ch = sum(1 for ins, outs, where in small_iso
               if not toks_equal(out_toks(rule(ins, S_ISO_NEW)),
                                 expect23("spectroscopy::isotope_pattern", ins, outs, where)))
    print(f"  pinning (ADR 0002 #23's table, against the port's expected results): with {label}, {n_ch} of "
          f"{len(small_iso)} isotope_pattern rows change")
    if must:
        check(n_ch > 0, f"spectroscopy: no golden row pins {label} on the port's #23 results")
# Coverage.
fft_lens = {ins[0][1] for ins, outs, w in srows.get("spectroscopy::fft_spectrum", [])}
check(set(range(17)) | {1 << k for k in range(14)} <= fft_lens, f"spectroscopy.tsv: fft_spectrum lengths {sorted(fft_lens)}")
fft_ok = [fl(ins)[1:] for ins, outs, w in srows.get("spectroscopy::fft_spectrum", []) if outs[0][0] == "v"]
check(any(any(math.isnan(v) for v in s_) for s_ in fft_ok) and any(any(math.isinf(v) for v in s_) for s_ in fft_ok)
      and any(any(0 < abs(v) < 2.0 ** -1022 for v in s_) for s_ in fft_ok) and any(any(abs(v) > 1e300 for v in s_) for s_ in fft_ok),
      "spectroscopy.tsv: fft_spectrum misses NaN, ±∞, subnormal or huge samples")
th_seen = {f2b(peaks_args(ins)[1]) for ins, outs, w in srows.get("spectroscopy::find_peaks", [])}
check({f2b(0.0), f2b(1.0), f2b(-0.0), f2b(math.inf)} <= th_seen and any(isnan_b(b) for b in th_seen),
      "spectroscopy.tsv: find_peaks misses a threshold edge")
check({(z, n) for z in range(256) for n in (0, 1, 0xFFFFFFFF)} <= {tuple(fl(ins)) for ins, o, w in srows.get("spectroscopy::bohr_energy_level", [])},
      "spectroscopy.tsv: bohr_energy_level misses a u8 Z")
iso_atoms = [pattern_atoms(ins) for ins, outs, w in srows.get("spectroscopy::isotope_pattern", [])]
check([] in iso_atoms and [(0, 1)] in iso_atoms and all([(z, 1)] in iso_atoms for z in range(256)),
      "spectroscopy.tsv: isotope_pattern misses the empty molecule or a u8 Z")
print(f"spectroscopy rows checked: {sum(len(v) for v in srows.values())} ({', '.join(f'{k} {len(v)}' for k, v in srows.items())})")
# ADR 0002 #21 (M4e-3): bohr_energy_level, bohr_transition_energy, rydberg_wavelength and SpectralSeries::wavelength are
# their exact rationals rounded once — −R·Z²/n², R·Z²·(n₂² − n₁²)/(n₁²·n₂²) and n₁²·n₂²/(R·Z²·(n₂² − n₁²)), R the double —
# and #22 reduced_mass is m₁·m₂/(m₁ + m₂) rounded once. Every Ok row's exact value, with Fraction: a row is noted exactly
# where it differs from the oracle's, and the note names it, bit for bit; so every unnoted row is the correctly rounded
# value too. Err rows are never noted.
def note_val_cr(where, adr):
    m = re.search(rf"\[ADR 0002 #{adr}: the port returns the correctly rounded (.+) = (\S+) \(f:([0-9a-f]+)\)\]$",
                  NOTES.get(where, ""))
    return None if m is None else (m.group(1), b2f(int(m.group(3), 16)))
def q_level(z, n): return -Fraction(S_RYD_EV) * z * z / (n * n)
def q_trans(z, lo, hi): return Fraction(S_RYD_EV) * z * z * (hi * hi - lo * lo) / (lo * lo * hi * hi)
def q_ryd(r, z, lo, hi): return Fraction(lo * lo * hi * hi) / (Fraction(r) * z * z * (hi * hi - lo * lo))
def q_series(name, line):
    nl = S_N_LOWER[name.split(".", 1)[1]]
    return q_ryd(S_RYD_H, 1, nl, nl + line)
def q_mass(m1, m2): return Fraction(m1) * Fraction(m2) / (Fraction(m1) + Fraction(m2))
EXACT_CR = {   # fn: (ADR, the expression the note names, the exact value from the row's inputs)
    "spectroscopy::bohr_energy_level": (21, "−R·Z²/n²", lambda a: q_level(*a)),
    "spectroscopy::bohr_transition_energy": (21, "R·Z²·(n₂² − n₁²)/(n₁²·n₂²)", lambda a: q_trans(*a)),
    "spectroscopy::rydberg_wavelength": (21, "n₁²·n₂²/(R∞·Z²·(n₂² − n₁²))", lambda a: q_ryd(S_RYD, *a)),
    "SpectralSeries::wavelength": (21, "n₁²·n₂²/(R_H·(n₂² − n₁²))", lambda a: q_series(*a)),
    "spectroscopy::reduced_mass": (22, "m₁·m₂/(m₁ + m₂)", lambda a: q_mass(*a)),
}
def is_midpoint(q):
    """1 if the positive rational q is a binary64 midpoint: a dyadic with a 54-bit odd significand in the normal range, or
    an odd multiple of 2⁻¹⁰⁷⁵ below 2⁻¹⁰²² (the subnormal grid's)."""
    num, den = q.numerator, q.denominator
    if den & (den - 1): return False
    while num % 2 == 0 and num:
        num //= 2
    if q < Fraction(1, 2 ** 1022):
        return (q * 2 ** 1075).denominator == 1 and (q * 2 ** 1075).numerator % 2 == 1
    return num.bit_length() == 54
ncr, okcr, tiescr, exactcr, worst = {}, {}, {}, {}, {}
for fn, (adr, what, rule) in EXACT_CR.items():
    ncr[fn] = okcr[fn] = tiescr[fn] = exactcr[fn] = 0
    worst[fn] = (0, None)
    for ins, outs, where in srows.get(fn, []):
        a = fl(ins)
        noted = f"[ADR 0002 #{adr}:" in NOTES.get(where, "")
        if outs[0][0] != "f":
            check(not noted, f"{where}: {fn}{tuple(a)}: an Err row noted")
            continue
        okcr[fn] += 1
        q = rule(a)
        exact = rn_frac(q)
        mid = is_midpoint(abs(q))
        tiescr[fn] += mid
        exactcr[fn] += (not mid) and rn_frac(q) == q if math.isfinite(exact) else 0
        differs = not same_bits(exact, outs[0][1])
        check(noted == differs, f"{where}: {fn}{tuple(a)}: exact {exact!r}, oracle {outs[0][1]!r}: the #{adr} note is "
              f"{'missing' if differs else 'stray'}")
        if noted:
            nv = note_val_cr(where, adr)
            check(nv is not None and nv[0] == what and same_bits(nv[1], exact),
                  f"{where}: the ADR 0002 #{adr} note does not name {what} = {exact!r}")
            ncr[fn] += 1
            if math.isfinite(outs[0][1]) and outs[0][1] != 0.0 and math.isfinite(exact):
                d_ = ulp_dist(f2b(exact), outs[0][2])
                if d_ > worst[fn][0]: worst[fn] = (d_, tuple(a))
for fn in EXACT_CR:
    w_ = f"; Rust's chain at most {worst[fn][0]:,} ulp from it on a finite row ({worst[fn][1]})" if worst[fn][1] else ""
    print(f"spectroscopy: ADR 0002 #{EXACT_CR[fn][0]}: {ncr[fn]} of {okcr[fn]} {fn.split('::')[1]} Ok rows noted with "
          f"the correctly rounded {EXACT_CR[fn][1]} (proven with Fraction; {tiescr[fn]} exact ties, {exactcr[fn]} exact "
          f"values){w_}")

# The structure behind #21's hard rows, proven independently of the generator. R's odd significand decides which
# exact values and ties exist: for RYDBERG_EV it is 683·P with P prime above 2³² (so squarefree), for RYDBERG and R_H it
# has a prime factor above 2³², which no product of n₁, n₂ < 2³² can cancel (λ is never dyadic: no exact value, no tie).
def mr_prime(n):
    if n < 2: return False
    for p_ in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p_ == 0: return n == p_
    d, r_ = n - 1, 0
    while d % 2 == 0: d //= 2; r_ += 1
    for a_ in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x_ = pow(a_, d, n)
        if x_ in (1, n - 1): continue
        for _ in range(r_ - 1):
            x_ = x_ * x_ % n
            if x_ == n - 1: break
        else: return False
    return True
_SMALL_PRIMES = [p_ for p_ in range(2, 1 << 17) if all(p_ % d for d in range(2, int(p_ ** 0.5) + 1))]
def factor_trial(n):
    """{prime: exponent} of n < 2⁶⁴ by trial division below 2¹⁷ and a Miller–Rabin prime cofactor."""
    out = {}
    for p_ in _SMALL_PRIMES:
        if p_ * p_ > n: break
        while n % p_ == 0: out[p_] = out.get(p_, 0) + 1; n //= p_
    if n > 1:
        check(n < 1 << 34 or mr_prime(n), f"factor_trial: composite cofactor {n}")
        out[n] = out.get(n, 0) + 1
    return out
def odd_sig(x):
    m_ = Fraction(x).numerator
    while m_ % 2 == 0: m_ //= 2
    return m_
f_ev = factor_trial(odd_sig(S_RYD_EV))
check(all(e_ == 1 for e_ in f_ev.values()) and max(f_ev) > 1 << 32,
      f"RYDBERG_EV's odd significand {odd_sig(S_RYD_EV)} = {f_ev}: not squarefree with a prime above 2³²")
for r_, nm in ((S_RYD, "RYDBERG"), (S_RYD_H, "R_H")):
    f_r = factor_trial(odd_sig(r_))
    check(max(f_r) > 1 << 32, f"{nm}'s odd significand {odd_sig(r_)} = {f_r}: no prime above 2³²")
# Bohr levels: odd(n)² | odd(Z)²·683·P forces odd(n) | odd(Z); the significand 683·P·j² (j = odd(Z)/odd(n)) has 53 bits
# for j = 1 and at least 56 otherwise (683·P·9 > 2⁵⁵): no tie, and an exact value exactly when n/Z is a power of two.
K_EV = odd_sig(S_RYD_EV)
lvl_bits = {(K_EV * j * j).bit_length() for z in range(1, 256) for j in range(1, 256, 2)
            if (z >> ((z & -z).bit_length() - 1)) % j == 0}
check(54 not in lvl_bits and min(b_ for b_ in lvl_bits if b_ > 53) >= 56 and 53 in lvl_bits,
      f"the Bohr levels' significand bit lengths {sorted(lvl_bits)}")
have_lvl = {tuple(fl(ins)) for ins, outs, w in srows.get("spectroscopy::bohr_energy_level", [])}
want_lvl = {(z, (z >> ((z & -z).bit_length() - 1)) << k) for z in range(1, 256) for k in range(33)
            if (z >> ((z & -z).bit_length() - 1)) << k <= 0xFFFFFFFF}
check(want_lvl <= have_lvl, f"spectroscopy.tsv: {len(want_lvl - have_lvl)} exact Bohr levels (n = 2^k·Z) missing")
check(tiescr["spectroscopy::bohr_energy_level"] == 0 and tiescr["spectroscopy::rydberg_wavelength"] == 0
      and tiescr["SpectralSeries::wavelength"] == 0, "a Bohr-level or Rydberg-wavelength row is a midpoint")
# Transition energies: with n₁ = g·x, n₂ = g·y (x, y coprime), the value is dyadic iff odd(x)·odd(y) | odd(Z) and odd(g)²
# divides 683·P·j²·w (w = odd(y² − x²), j = odd(Z)/(odd(x)·odd(y))); the odd significand is that over odd(g)². Every
# such (Z, n₁, n₂) with n₂ < 2³² is enumerated here and the midpoints (54 bits) collected.
def transition_ties():
    ties, exact = set(), 0
    U = 0xFFFFFFFF
    for ox in range(1, 256, 2):
        for oy in range(1, 256, 2):
            if ox * oy > 255 or math.gcd(ox, oy) != 1: continue
            for a_ in range(33):
                for b_ in range(33):
                    if a_ and b_: continue
                    x_, y_ = ox << a_, oy << b_
                    if not x_ < y_ <= U: continue
                    fw = {}
                    for part in (y_ - x_, y_ + x_):
                        for p_, e_ in factor_trial(part).items():
                            if p_ != 2: fw[p_] = fw.get(p_, 0) + e_
                    for z in range(1, 256):
                        zo = z >> ((z & -z).bit_length() - 1)
                        if zo % (ox * oy): continue
                        j = zo // (ox * oy)
                        ft = dict(f_ev)
                        for p_, e_ in list(fw.items()) + [(p_, 2 * e_) for p_, e_ in factor_trial(j).items()]:
                            ft[p_] = ft.get(p_, 0) + e_
                        T = 1
                        for p_, e_ in ft.items(): T *= p_ ** e_
                        ogs = [1]
                        for p_, e_ in ft.items():
                            ogs = [d * p_ ** k for d in ogs for k in range(e_ // 2 + 1)]
                        for og in ogs:
                            if og * y_ > U: continue
                            bl = (T // (og * og)).bit_length()
                            if bl > 54: continue
                            g = og
                            while g * y_ <= U:
                                if bl == 54: ties.add((z, g * x_, g * y_))
                                else: exact += 1
                                g *= 2
    return ties, exact
t_ties, t_exact = transition_ties()
have_tr = {tuple(fl(ins)) for ins, outs, w in srows.get("spectroscopy::bohr_transition_energy", [])}
have_ry = {tuple(fl(ins)) for ins, outs, w in srows.get("spectroscopy::rydberg_wavelength", [])}
check(t_ties <= have_tr and t_ties <= have_ry, f"spectroscopy.tsv: {len(t_ties - have_tr)} transition ties missing")
check(all(is_midpoint(q_trans(*t)) for t in t_ties), "a searched transition tie is not a midpoint")
check(tiescr["spectroscopy::bohr_transition_energy"] == len(t_ties),
      f"the corpus has {tiescr['spectroscopy::bohr_transition_energy']} midpoint transition rows, the search {len(t_ties)}")
print(f"spectroscopy: ADR 0002 #21's structure: R's odd significands {K_EV} = {' · '.join(map(str, sorted(f_ev)))} "
      f"(RYDBERG_EV), {odd_sig(S_RYD)} and {odd_sig(S_RYD_H)} (each with a prime above 2³²): no level is a midpoint, a level "
      f"is exact iff n/Z = 2^k ({len(want_lvl)} such rows, all present), no Rydberg wavelength is dyadic; the exhaustive "
      f"search of bohr_transition_energy's domain finds {t_exact} exact values and {len(t_ties)} ties "
      f"({sorted(t_ties)[0]} … ×2^k), every tie in the corpus")
mass_ties = [tuple(fl(ins)) for ins, outs, w in srows.get("spectroscopy::reduced_mass", [])
             if outs[0][0] == "f" and is_midpoint(q_mass(*fl(ins)))]
check(len(mass_ties) >= 40 and any(q_mass(*t) < Fraction(1, 2 ** 1022) for t in mass_ties)
      and any(q_mass(*t) >= Fraction(1, 2 ** 1022) for t in mass_ties), f"reduced_mass: only {len(mass_ties)} tie rows")
eq_m = [(fl(ins), outs) for ins, outs, w in srows.get("spectroscopy::reduced_mass", []) if outs[0][0] == "f" and fl(ins)[0] == fl(ins)[1]]
eq_miss = sum(1 for (a, b_), outs in eq_m if not same_bits(outs[0][1], rn_frac(Fraction(a) / 2)))
zn = [(fl(ins), outs) for ins, outs, w in srows.get("spectroscopy::bohr_energy_level", []) if outs[0][0] == "f" and fl(ins)[0] == fl(ins)[1]]
zn_miss = sum(1 for a, outs in zn if not same_bits(outs[0][1], -S_RYD_EV))
near_rep = sum(1 for ins, outs, w in srows.get("spectroscopy::isotope_pattern", []) if outs[0][0] == "v"
               for k in range(outs[0][1]) if outs[2 + 2 * k][1] == S_REPORT)
near_rep23 = 0
for ins, outs, w in srows.get("spectroscopy::isotope_pattern", []):
    e23 = expect23("spectroscopy::isotope_pattern", ins, outs, w)
    if e23[0][0] == "v":
        near_rep23 += sum(1 for k in range(e23[0][1]) if e23[2 + 2 * k][1] == S_REPORT)
print(f"spectroscopy: reduced_mass has {len(mass_ties)} tie rows (normal and subnormal); the oracle misses m/2 on {eq_miss} "
      f"of {len(eq_m)} equal-mass rows and −R on {zn_miss} of the {len(zn)} Z = n level rows (each noted); {near_rep} "
      f"isotope peaks exactly at the 0.001 report threshold ({near_rep23} in the port's #23 results)")

# ── 18. stoichiometry.tsv (M4f): Rust's semantics recomputed, ADR 0002 #7's mapping proven ─────────────────────────
# Every row is recomputed from rust-old/src/stoichiometry.rs's operations: balance_equation and assign_oxidation_states
# on Python integers checked against i64 at every Rust checked_* call (i128 in the verification), the float functions
# on IEEE doubles in Rust's order. Then:
# - balance_equation: every row is checked again by exact linear algebra (Fraction RREF): an Ok result is the
#   primitive positive basis vector of a one-dimensional null space; an Err is a null space of another dimension, a
#   basis vector with a zero or mixed signs or past u32, or (counted) Rust's i64 overflow where a valid u32 solution
#   exists (the port's Ok there since the overflow audit's B7, ADR 0002 #39: §25). The `[ADR 0002 #7: no
#   integer-matrix counterpart` notes sit exactly on the rows whose f64 matrix has no
#   integer one (an entry Rust's test refuses, an integer beyond i64, ragged rows), all Err(InvalidInput);
# - assign_oxidation_states: the port's exact rational states (ADR 0002 #7: Rust's algorithm on Fractions, the i32
#   range applied to the value) are computed; a row is noted exactly where Rust's Err is its "non-integer average" and
#   the port is Ok, the note's tokens are the port's result (reduced, the denominator positive), so every unnoted row's
#   port result is the oracle's;
# - the typed constants: EMPIRICAL_TOLERANCE (0.1) × ½ and × 2, INTEGER_TOLERANCE (10⁻⁹: the generator's and the
#   harness's conversion) ± 1 ulp and MAX_COMPOSITION_ENTRY (2⁵³) ± 1 each change some row;
# - ADR 0002 #24–#26 (M4f-2): the port's results recomputed with Fraction, independently of the generator — percent_yield
#   RN(100·a/t), theoretical_yield_moles RN(m·p/c) at the exact limiting reagent, percent_composition RN(mass·count·100/mw)
#   with Rust's mw, limiting_reagent's exact least m/c (the first on a tie), empirical_formula with exact moles, minimum,
#   ratios, tolerance 1/10 and counts — and every row of the five checked: noted exactly where the port's result differs
#   from the oracle's, the note naming it (a value bit for bit, an index, or the whole result's tokens), so every unnoted
#   row is the port's result too; the claims behind the port's algorithm checked on the corpus (a ratio of 2³² − ½ or
#   more is always the ComputationError; the ratios the port keeps lie in [1, 2³²)), Rust's Err paths classified by what
#   the port returns, the hard rows' coverage counted, and the exact 1/10 pinned (1/10 − 2⁻⁸⁰ and the nearest larger
#   deviation each change some row).
stsrc = (root / "rust-old/src/stoichiometry.rs").read_text()
ST_MAX_ENTRY = float(re.search(r"^const MAX_COMPOSITION_ENTRY: f64 = ([0-9_.]+);", stsrc, re.M).group(1).replace("_", ""))
ST_INT_TOL = float(re.search(r"^const INTEGER_TOLERANCE: f64 = ([0-9e.-]+);", stsrc, re.M).group(1))
ST_EMP_TOL = float(re.search(r"^const EMPIRICAL_TOLERANCE: f64 = ([0-9.]+);", stsrc, re.M).group(1))
check((ST_MAX_ENTRY, ST_INT_TOL, ST_EMP_TOL) == (2.0 ** 53, 1e-9, 0.1),
      f"stoichiometry.rs: constants parsed as {ST_MAX_ENTRY}, {ST_INT_TOL}, {ST_EMP_TOL}")
I64_MIN, I64_MAX, I128_MIN, I128_MAX = -2 ** 63, 2 ** 63 - 1, -2 ** 127, 2 ** 127 - 1
I32_MIN, I32_MAX, U32_MAX = -2 ** 31, 2 ** 31 - 1, 2 ** 32 - 1
class Ovf(Exception):
    """A Rust checked_* call returned None: the overflow Err (ComputationError)."""
def c64(v):
    if v < I64_MIN or v > I64_MAX: raise Ovf
    return v
def tdiv(a, b):
    """Rust's i64 `/`: truncation toward zero."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b > 0) else -q
def st_reduce(a, lo, n):
    g = 0
    for v in a[lo:lo + n]: g = math.gcd(g, abs(v))
    if g > 1:
        if g > I64_MAX: raise Ovf
        for k in range(lo, lo + n): a[k] = tdiv(a[k], g)
def st_matrix(ins):
    """Rust's &[Vec<f64>] from `v:m`, then per row `v:len` and its floats."""
    m, k, rows = ins[0][1], 1, []
    for _ in range(m):
        ln = ins[k][1]
        rows.append([t[1] for t in ins[k + 1:k + 1 + ln]])
        k += 1 + ln
    check(k == len(ins), "balance_equation: malformed matrix tokens")
    return rows
def st_counterpart(mat, tol=None):
    """The port's integer matrix (ADR 0002 #7), row-major, or None: equal rows, every entry finite and within 10⁻⁹ of
    round(v) (Rust's test), that integer below 2⁶³ in magnitude."""
    tol = ST_INT_TOL if tol is None else tol
    if mat and any(len(r) != len(mat[0]) for r in mat): return None
    out = []
    for row in mat:
        for v in row:
            if not math.isfinite(v): return None
            r = r_round(v)
            if abs(v - r) > tol or abs(r) >= 2.0 ** 63: return None
            out.append(int(r))
    return out
def st_balance(mat, tol=None, bound=None):
    """Rust 1.2.0's balance_equation on its f64 matrix: the u32 coefficients, or Err."""
    tol = ST_INT_TOL if tol is None else tol
    bound = ST_MAX_ENTRY if bound is None else bound
    if not mat: return Err("InvalidInput")
    n = len(mat[0])
    if n < 2 or any(len(r) != n for r in mat): return Err("InvalidInput")
    a = []
    for row in mat:
        for v in row:
            r = r_round(v) if math.isfinite(v) else v
            if not math.isfinite(v) or abs(v - r) > tol or abs(r) > bound: return Err("InvalidInput")
            a.append(int(r))
    return st_balance_int(a, len(mat), n, bound)
def st_balance_int(a, m, n, bound=None):
    """Rust's balance_equation on ADR 0002 #7's integer matrix a, row-major: Rust's checks with its 2⁵³ bound on the
    integers, then Rust's elimination with its i64 checks; the u32 coefficients, or Err. (The port's own result, on
    unbounded integers where these checks overflow, is §25's b7_port: ADR 0002 #39.)"""
    bound = ST_MAX_ENTRY if bound is None else bound
    if m < 1 or n < 2: return Err("InvalidInput")
    if any(abs(v) > bound for v in a): return Err("InvalidInput")
    a = list(a)
    orig = list(a)
    try:
        rank, free = 0, None
        for col in range(n):
            p = next((r for r in range(rank, m) if a[r * n + col] != 0), None)
            if p is None:
                if free is not None: return Err("ComputationError")       # ambiguous
                free = col
                continue
            if p != rank:
                for j in range(n): a[p * n + j], a[rank * n + j] = a[rank * n + j], a[p * n + j]
            pr = rank * n
            st_reduce(a, pr, n)
            if a[pr + col] < 0:
                for j in range(n): a[pr + j] = c64(-a[pr + j])
            piv = a[pr + col]
            for r in range(m):
                rr, fac = r * n, a[r * n + col]
                if r == rank or fac == 0: continue
                for j in range(n): a[rr + j] = c64(c64(piv * a[rr + j]) - c64(fac * a[pr + j]))
                st_reduce(a, rr, n)
            rank += 1
        if free is None: return Err("ComputationError")                   # only the trivial solution
        lcm = 1
        for k in range(rank):
            pv = next((v for v in a[k * n:k * n + n] if v != 0), None)
            if pv is not None:
                g = math.gcd(abs(lcm), abs(pv))
                if g > I64_MAX: raise Ovf
                lcm = c64(tdiv(lcm, g) * pv)
        x = [0] * n
        x[free] = lcm
        for k in range(rank):
            row = a[k * n:k * n + n]
            c = next((j for j, v in enumerate(row) if v != 0), None)
            if c is not None: x[c] = c64(-c64(row[free] * tdiv(lcm, row[c])))
        st_reduce(x, 0, n)
        if 0 in x: return Err("ComputationError")
        pos = sum(1 for v in x if v > 0)
        if pos == 0: x = [c64(-v) for v in x]
        elif pos != n: return Err("ComputationError")
        if any(v < 0 or v > U32_MAX for v in x): return Err("ComputationError")
        for i in range(m):
            s_ = 0
            for j in range(n):
                s_ += orig[i * n + j] * x[j]                                   # |term| < 2⁹⁵: checked_mul can't fail
                if s_ < I128_MIN or s_ > I128_MAX: return Err("ComputationError")
            if s_ != 0: return Err("ComputationError")
        return x
    except Ovf:
        return Err("ComputationError")
def st_nullspace(a, m, n):
    """The exact null space of the integer matrix: (dimension, its primitive basis vector when the dimension is 1)."""
    M = [[Fraction(a[i * n + j]) for j in range(n)] for i in range(m)]
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = M[r][c]
        M[r] = [v / inv for v in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f_ = M[i][c]
                M[i] = [vi - f_ * vr for vi, vr in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(n) if c not in pivots]
    if len(free) != 1: return len(free), None
    f_ = free[0]
    v = [Fraction(0)] * n
    v[f_] = Fraction(1)
    for i, c in enumerate(pivots): v[c] = -M[i][f_]
    den = 1
    for q in v: den = den * q.denominator // math.gcd(den, q.denominator)
    w = [int(q * den) for q in v]
    g = 0
    for q in w: g = math.gcd(g, abs(q))
    w = [q // g for q in w]
    if all(q <= 0 for q in w): w = [-q for q in w]
    return 1, w

ST_CAT = {z: e[5] for z, e in ELEMS.items()}
ST_METALS = {"AlkaliMetal", "AlkalineEarth", "TransitionMetal", "PostTransitionMetal", "Lanthanide", "Actinide"}
NONINT = "non-integer"
def st_exact(num, den, port):
    """exact_oxidation_state: Rust's integer (port=False; Err(NONINT) for its non-integer average) or the port's exact
    Fraction (port=True; ADR 0002 #7: the i32 range on the value)."""
    if den == 0: return Err("InvalidInput")
    if num == I64_MIN and den == -1: return Err("ComputationError")
    q = Fraction(num, den)
    if q.denominator != 1:
        if not port: return Err(NONINT)
        fl_ = math.floor(q)
        return Err("ComputationError") if (fl_ >= I32_MAX or fl_ < I32_MIN) else q
    return Err("ComputationError") if (q < I32_MIN or q > I32_MAX) else q
def st_rule(z, has_f, hydride):
    if z == 9: return -1
    if z == 1: return -1 if hydride else 1
    if z == 8: return None if has_f else -2
    return {"AlkaliMetal": 1, "AlkalineEarth": 2}.get(ST_CAT[z])
def st_oxidation(atoms, charge, port):
    """Rust 1.2.0's assign_oxidation_states (port=False: integer states, Err(NONINT) for its non-integer Err) or the
    port's (port=True: Fractions)."""
    if not atoms: return Err("InvalidInput")
    first = atoms[0][0]
    single, has_f, others_metal, total = True, False, True, 0
    try:
        for z, c in atoms:
            if c == 0: return Err("InvalidInput")
            if z not in ST_CAT: return Err("InvalidElement")
            single &= z == first
            has_f |= z == 9
            if z != 1 and ST_CAT[z] not in ST_METALS: others_metal = False
            total = c64(total + c)
        if single:
            os_ = st_exact(charge, total, port)
            return os_ if isinstance(os_, Err) else [(z, os_) for z, _ in atoms]
        assigned, oxygen, unknown = 0, 0, None
        for z, c in atoms:
            os_ = st_rule(z, has_f, others_metal)
            if os_ is not None:
                assigned = c64(assigned + c64(os_ * c))
                if z == 8: oxygen = c64(oxygen + c)
            elif unknown is None: unknown = (z, c)
            elif unknown[0] == z: unknown = (z, c64(unknown[1] + c))
            else: return Err("ComputationError")
        remaining = c64(charge - assigned)
        solved = None
        if unknown is not None:
            os_ = st_exact(remaining, unknown[1], port)
            if isinstance(os_, Err): return os_
            solved = (unknown[0], os_)
        elif remaining != 0:
            if oxygen <= 0: return Err("ComputationError")
            os_ = st_exact(c64(remaining + c64(oxygen * -2)), oxygen, port)
            if isinstance(os_, Err): return os_
            if not (-2 < os_ < 0): return Err("ComputationError")
            solved = (8, os_)
        return [(z, solved[1] if solved and solved[0] == z else Fraction(st_rule(z, has_f, others_metal)))
                for z, _ in atoms]
    except Ovf:
        return Err("ComputationError")
def st_pairs(ins):
    n = ins[0][1]
    return [(ins[1 + 2 * k][1], ins[2 + 2 * k][1]) for k in range(n)]
def st_out_rust(r):
    if isinstance(r, Err): return [("E", "ComputationError" if r.name == NONINT else r.name, None)]
    return [("v", len(r), None)] + [t for z, q in r for t in (("i", z, None), ("i", int(q), None))]
def st_out_port(r):
    if isinstance(r, Err): return [("E", r.name, None)]
    return [("v", len(r), None)] + [t for z, q in r for t in (("i", z, None), ("i", q.numerator, None), ("i", q.denominator, None))]
def st_limiting(mo, co):
    if len(mo) != len(co) or not mo: return Err("InvalidInput")
    best, idx = math.inf, 0
    for i, (m_, c_) in enumerate(zip(mo, co)):
        if not is_non_negative(m_) or not is_positive(c_): return Err("InvalidInput")
        r = m_ / c_
        if r < best: best, idx = r, i
    return idx
def st_two_slices(ins):
    n1 = ins[0][1]
    n2 = ins[1 + n1][1]
    return [t[1] for t in ins[1:1 + n1]], [t[1] for t in ins[2 + n1:2 + n1 + n2]], ins[2 + n1 + n2:]
def st_yield(mo, co, p):
    if not is_positive(p): return Err("InvalidInput")
    k = st_limiting(mo, co)
    return k if isinstance(k, Err) else dv(mo[k], co[k]) * p
def st_percent_yield(a, t):
    if not is_positive(t): return Err("InvalidInput")
    if not is_non_negative(a): return Err("InvalidInput")
    return dv(a, t) * 100.0
def st_composition(pairs):
    mw = 0.0
    for z, c in pairs:
        if z not in ELEMS: return Err("InvalidElement")
        mw += ELEMS[z][3] * float(c)
    if not is_positive(mw): return Err("ComputationError")
    return [(z, dv(ELEMS[z][3] * float(c), mw) * 100.0) for z, c in pairs]
def st_empirical(els, tol=None):
    tol = ST_EMP_TOL if tol is None else tol
    if not els: return Err("InvalidInput")
    moles = []
    for z, p in els:
        if z not in ELEMS: return Err("InvalidElement")
        if not is_positive(p): return Err("InvalidInput")
        moles.append(dv(p, ELEMS[z][3]))
    mn = math.inf
    for v in moles: mn = v if v < mn else mn                                 # f64::min: no NaN here
    if not is_positive(mn): return Err("ComputationError")
    rat = []
    for v in moles:
        r = dv(v, mn)
        if not math.isfinite(r): return Err("ComputationError")
        rat.append(r)
    mult = next((float(k) for k in range(1, 7) if all(abs(r * k - r_round(r * k)) <= tol for r in rat)), None)
    if mult is None: return Err("ComputationError")
    out = []
    for (z, _), r in zip(els, rat):
        c = r_round(r * mult)
        if c > 4294967295.0: return Err("ComputationError")
        out.append((z, int(c)))
    return out
def st_toks(fn, r):
    if isinstance(r, Err): return [("E", r.name, None)]
    if fn == "stoichiometry::balance_equation": return [("v", len(r), None)] + [("i", v, None) for v in r]
    if fn == "stoichiometry::percent_composition":
        return [("v", len(r), None)] + [t for z, p in r for t in (("i", z, None), ("f", p, f2b(p)))]
    if fn == "stoichiometry::empirical_formula":
        return [("v", len(r), None)] + [t for z, c in r for t in (("i", z, None), ("i", c, None))]
    if isinstance(r, int): return [("i", r, None)]
    return [("f", r, f2b(r))]
# ADR 0002 #24–#26 (M4f-2): the port's results, exact (Fraction; rn_frac rounds once, correctly, ties to even).
def st_port_limiting(mo, co):
    if len(mo) != len(co) or not mo: return Err("InvalidInput")
    best = None
    for i, (m_, c_) in enumerate(zip(mo, co)):
        if not is_non_negative(m_) or not is_positive(c_): return Err("InvalidInput")
        q = Fraction(m_) / Fraction(c_)
        if best is None or q < best[0]: best = (q, i)
    return best[1]
def st_port_yield(mo, co, p):
    if not is_positive(p): return Err("InvalidInput")
    k = st_port_limiting(mo, co)
    if isinstance(k, Err): return k
    m_, c_ = mo[k], co[k]
    return dv(m_, c_) * p if m_ == 0 else rn_frac(Fraction(m_) * Fraction(p) / Fraction(c_))
def st_port_percent_yield(a, t):
    if not is_positive(t) or not is_non_negative(a): return Err("InvalidInput")
    return dv(a, t) * 100.0 if a == 0 else rn_frac(Fraction(a) * 100 / Fraction(t))
def st_port_composition(pairs):
    r = st_composition(pairs)                   # Rust's checks and its mw
    if isinstance(r, Err): return r
    mw = 0.0
    for z, c in pairs: mw += ELEMS[z][3] * float(c)
    return [(z, 0.0 if c == 0 else rn_frac(Fraction(ELEMS[z][3]) * c * 100 / Fraction(mw))) for z, c in pairs]
def st_rnd(q): return math.floor(q + Fraction(1, 2))       # Rust's round of a q ≥ 0: halves away from zero
def st_port_ratios(els):
    """The exact ratios rᵢ = molesᵢ / min moles (or an Err of Rust's validation)."""
    if not els: return Err("InvalidInput")
    moles = []
    for z, p in els:
        if z not in ELEMS: return Err("InvalidElement")
        if not is_positive(p): return Err("InvalidInput")
        moles.append(Fraction(p) / Fraction(ELEMS[z][3]))
    mn = min(moles)
    return [q / mn for q in moles]
def st_port_empirical(els, tol=Fraction(1, 10)):
    rs = st_port_ratios(els)
    if isinstance(rs, Err): return rs
    mult = next((k for k in range(1, 7) if all(abs(r * k - st_rnd(r * k)) <= tol for r in rs)), None)
    if mult is None: return Err("ComputationError")
    out = []
    for (z, _), r in zip(els, rs):
        c = st_rnd(r * mult)
        if c > U32_MAX: return Err("ComputationError")
        out.append((z, c))
    return out
def st_empirical_path(els):
    """Rust's empirical_formula's outcome by path: 'Ok', 'InvalidInput', 'InvalidElement', 'non-positive min',
    'ratio not finite', 'no multiplier', 'count past u32'."""
    if not els: return "InvalidInput"
    moles = []
    for z, p in els:
        if z not in ELEMS: return "InvalidElement"
        if not is_positive(p): return "InvalidInput"
        moles.append(dv(p, ELEMS[z][3]))
    mn = math.inf
    for v in moles: mn = v if v < mn else mn
    if not is_positive(mn): return "non-positive min"
    rat = [dv(v, mn) for v in moles]
    if not all(math.isfinite(r) for r in rat): return "ratio not finite"
    mult = next((float(k) for k in range(1, 7) if all(abs(r * k - r_round(r * k)) <= ST_EMP_TOL for r in rat)), None)
    if mult is None: return "no multiplier"
    if any(r_round(r * mult) > 4294967295.0 for r in rat): return "count past u32"
    return "Ok"
NOTE24_VAL = re.compile(r"\[ADR 0002 #24: the port returns the correctly rounded (.+) = (\S+) \(f:([0-9a-f]{16})\)\]$")
NOTE25 = re.compile(r"\[ADR 0002 #25: the port returns the exact limiting reagent \(mᵢ·c_j < m_j·cᵢ, the first on a tie\) = "
                    r"(\d+) \(i:(\d+)\)\]$")
NOTE_TOKS = re.compile(r" \(([^()]*)\)\]$")
ST_EXACT = {   # fn: (ADR, what a #24 value note names)
    "stoichiometry::percent_yield": (24, "100·a/t"),
    "stoichiometry::theoretical_yield_moles": (24, "m·p/c at the exact limiting reagent"),
    "stoichiometry::percent_composition": (24, None),
    "stoichiometry::limiting_reagent": (25, None),
    "stoichiometry::empirical_formula": (26, None),
}
st_x_noted, st_x_ok = {}, {}
def st_check_port(fn, outs, port, where, note):
    """The ADR 0002 #24–#26 note of one row: present exactly where the port's result differs from the oracle's, naming
    it."""
    adr, what = ST_EXACT[fn]
    pt = st_toks(fn, port)
    differs = not toks_equal(outs, pt)
    noted = f"[ADR 0002 #{adr}:" in note
    check(noted == differs, f"{where}: {fn}: the port's {pt} vs the oracle's {outs}: the #{adr} note is "
          f"{'missing' if differs else 'stray'}")
    st_x_ok[fn] = st_x_ok.get(fn, 0) + 1
    if not noted: return
    st_x_noted[fn] = st_x_noted.get(fn, 0) + 1
    if what is not None:
        m = NOTE24_VAL.search(note)
        check(m is not None and m.group(1) == what and len(pt) == 1 and pt[0][0] == "f" and int(m.group(3), 16) == pt[0][2],
              f"{where}: the #24 note does not name {what} = {port!r}")
    elif adr == 25:
        m = NOTE25.search(note)
        check(m is not None and int(m.group(1)) == port and int(m.group(2)) == port, f"{where}: the #25 note is not {port}")
    else:
        m = NOTE_TOKS.search(note)
        toks = [parse_tok(t, where) for t in m.group(1).split(" ")] if m else []
        check(toks_equal(toks, pt), f"{where}: the #{adr} note's result is not the port's {pt}")
ST_CLASS = {"stoichiometry::balance_equation": "A", "stoichiometry::assign_oxidation_states": "A",
            "stoichiometry::limiting_reagent": "A", "stoichiometry::theoretical_yield_moles": "B",
            "stoichiometry::percent_yield": "B", "stoichiometry::percent_composition": "B",
            "stoichiometry::empirical_formula": "A"}
NOTE7_NONE = "[ADR 0002 #7: no integer-matrix counterpart ("
NOTE7_OX = "[ADR 0002 #7: the port returns the exact average oxidation states"
st_count, st_none, st_round, st_ox, st_bal_ok, st_bal_err = {}, {}, 0, 0, 0, {}
st_ovf_valid, st_ovf_examples = 0, []
st_rows = {}
st_x_rows = {}
for fn, cls, ins, outs, where in R.get("stoichiometry", []):
    st_count[fn] = st_count.get(fn, 0) + 1
    st_rows.setdefault(fn, []).append((ins, outs, where))
    check(ST_CLASS.get(fn) == cls, f"{where}: {fn} class {cls}")
    note = NOTES.get(where, "")
    if fn == "stoichiometry::balance_equation":
        mat = st_matrix(ins)
        rust = st_balance(mat)
        check(toks_equal(outs, st_toks(fn, rust)), f"{where}: balance_equation -> {outs}, recomputed {rust}")
        cp = st_counterpart(mat)
        noted = NOTE7_NONE in note
        check(noted == (cp is None), f"{where}: the #7 no-counterpart note is {'stray' if noted else 'missing'}")
        if cp is None:
            check(outs == [("E", "InvalidInput", None)], f"{where}: a matrix without an integer counterpart is not Err(InvalidInput)")
            why = re.search(r"no integer-matrix counterpart \(([^)]*)\)", note)
            st_none[why.group(1) if why else "?"] = st_none.get(why.group(1) if why else "?", 0) + 1
            continue
        st_round += any(math.isfinite(v) and v != r_round(v) for row in mat for v in row)
        m_, n_ = len(mat), len(mat[0]) if mat else 0
        # Rust's algorithm on the integer counterpart gives the oracle's result: the f64 entry test is the one step the
        # integer API lacks, and every entry Rust accepts is its rounded integer (the port's result is §25's)
        check(toks_equal(outs, st_toks(fn, st_balance_int(cp, m_, n_))), f"{where}: the integer counterpart's result")
        if m_ == 0 or n_ < 2 or any(abs(v) > ST_MAX_ENTRY for v in cp):
            continue
        dim, w = st_nullspace(cp, m_, n_)
        valid = dim == 1 and all(q > 0 for q in w) and max(w) <= U32_MAX
        if isinstance(rust, Err):
            st_bal_err[dim if dim != 1 else ("valid" if valid else "invalid")] = st_bal_err.get(
                dim if dim != 1 else ("valid" if valid else "invalid"), 0) + 1
            if valid:
                st_ovf_valid += 1
                if len(st_ovf_examples) < 2: st_ovf_examples.append((mat, w))
        else:
            st_bal_ok += 1
            check(valid and rust == w, f"{where}: balance_equation's Ok {rust} is not the exact null vector {w} (dim {dim})")
    elif fn == "stoichiometry::assign_oxidation_states":
        atoms, charge = st_pairs(ins), ins[-1][1]
        check(len(ins) == 2 + 2 * len(atoms) and all(0 <= z <= 255 and 0 <= c <= U32_MAX for z, c in atoms)
              and I32_MIN <= charge <= I32_MAX, f"{where}: (u8, u32) atoms and an i32 charge")
        rust = st_oxidation(atoms, charge, False)
        port = st_oxidation(atoms, charge, True)
        check(toks_equal(outs, st_out_rust(rust)), f"{where}: assign_oxidation_states -> {outs}, recomputed {rust}")
        nonint = isinstance(rust, Err) and rust.name == NONINT
        if not nonint:
            check(st_out_port(port) == ([("E", rust.name, None)] if isinstance(rust, Err) else
                  [("v", len(rust), None)] + [t for z, q in rust for t in (("i", z, None), ("i", int(q), None), ("i", 1, None))]),
                  f"{where}: the port's states differ from Rust's on a row Rust answers")
        expect_note = nonint and not isinstance(port, Err)
        noted = NOTE7_OX in note
        check(noted == expect_note, f"{where}: the #7 oxidation note is {'stray' if noted else 'missing'}")
        if nonint and isinstance(port, Err):
            check(port.name == "ComputationError", f"{where}: the port's Err for a non-integer average is {port.name}")
        if noted:
            m = re.search(r" \(([^()]*)\)\]$", note)
            toks = [parse_tok(t, where) for t in m.group(1).split(" ")] if m else []
            check(toks == st_out_port(port), f"{where}: the #7 note is not the port's exact states {port}")
            check(all(q.denominator > 0 and math.gcd(q.numerator, q.denominator) == 1 for _, q in port),
                  f"{where}: a state not reduced")
            st_ox += 1
    elif fn in ("stoichiometry::limiting_reagent", "stoichiometry::theoretical_yield_moles"):
        mo, co, rest = st_two_slices(ins)
        r = st_limiting(mo, co) if fn.endswith("limiting_reagent") else st_yield(mo, co, rest[0][1])
        check(toks_equal(outs, st_toks(fn, r)), f"{where}: {fn} -> {outs}, recomputed {r}")
        port = st_port_limiting(mo, co) if fn.endswith("limiting_reagent") else st_port_yield(mo, co, rest[0][1])
        st_check_port(fn, outs, port, where, note)
        st_x_rows.setdefault(fn, []).append((mo, co, rest, r, port))
    elif fn == "stoichiometry::percent_yield":
        r = st_percent_yield(ins[0][1], ins[1][1])
        check(toks_equal(outs, st_toks(fn, r)), f"{where}: percent_yield -> {outs}, recomputed {r}")
        port = st_port_percent_yield(ins[0][1], ins[1][1])
        st_check_port(fn, outs, port, where, note)
        st_x_rows.setdefault(fn, []).append((ins[0][1], ins[1][1], r, port))
    elif fn == "stoichiometry::percent_composition":
        r = st_composition(st_pairs(ins))
        check(toks_equal(outs, st_toks(fn, r)), f"{where}: percent_composition -> {outs}, recomputed {r}")
        st_check_port(fn, outs, st_port_composition(st_pairs(ins)), where, note)
    elif fn == "stoichiometry::empirical_formula":
        r = st_empirical(st_pairs(ins))
        check(toks_equal(outs, st_toks(fn, r)), f"{where}: empirical_formula -> {outs}, recomputed {r}")
        port = st_port_empirical(st_pairs(ins))
        st_check_port(fn, outs, port, where, note)
        st_x_rows.setdefault(fn, []).append((st_pairs(ins), r, port))
    else:
        check(False, f"{where}: stoichiometry.tsv function {fn} has no rule here")
check(set(st_count) == set(ST_CLASS), f"stoichiometry: functions {sorted(st_count)}")
print(f"stoichiometry rows checked: {sum(st_count.values())} ({', '.join(f'{k.split(chr(58))[-1]} {v}' for k, v in sorted(st_count.items()))})")
print(f"stoichiometry: ADR 0002 #7: {st_ox} assign_oxidation_states rows noted with the port's exact rational states "
      f"(Fraction, reduced; the oracle's non-integer-average Err), {sum(st_none.values())} balance_equation rows without "
      f"an integer counterpart ({', '.join(f'{v} {k}' for k, v in sorted(st_none.items()))}), "
      f"{st_round} replayed with a near-integer entry rounded as Rust does")
print(f"stoichiometry: balance_equation by exact linear algebra: {st_bal_ok} Ok rows the primitive positive null "
      f"vector; Err rows by null-space dimension / vector: {dict(sorted(st_bal_err.items(), key=str))}; "
      f"{st_ovf_valid} rows where Rust's i64 elimination overflows although a valid u32 solution exists")
for mat, w in st_ovf_examples:
    print(f"  e.g. {[[int(v) for v in row] for row in mat]}: Err(ComputationError), the exact answer {w}")
# The typed constants: each must change some golden row when perturbed.
def st_changes(fn, rule):
    return sum(1 for ins, outs, w in st_rows.get(fn, []) if not toks_equal(st_toks(fn, rule(ins)), outs))
ST_PINS = [
    ("EMPIRICAL_TOLERANCE × 0.5", "stoichiometry::empirical_formula", lambda a: st_empirical(st_pairs(a), tol=0.05)),
    ("EMPIRICAL_TOLERANCE × 2", "stoichiometry::empirical_formula", lambda a: st_empirical(st_pairs(a), tol=0.2)),
    ("INTEGER_TOLERANCE + 1 ulp", "stoichiometry::balance_equation", lambda a: st_balance(st_matrix(a), tol=nextup(1e-9))),
    ("INTEGER_TOLERANCE − 1 ulp", "stoichiometry::balance_equation", lambda a: st_balance(st_matrix(a), tol=nextdown(1e-9))),
    ("MAX_COMPOSITION_ENTRY − 1", "stoichiometry::balance_equation", lambda a: st_balance(st_matrix(a), bound=2.0 ** 53 - 1)),
    ("MAX_COMPOSITION_ENTRY + 2", "stoichiometry::balance_equation", lambda a: st_balance(st_matrix(a), bound=2.0 ** 53 + 2)),
]
for label, fn, rule in ST_PINS:
    ch = st_changes(fn, rule)
    print(f"  pinning: with {label}, {ch} of {len(st_rows.get(fn, []))} {fn.split('::')[1]} rows change")
    check(ch > 0, f"stoichiometry: no golden {fn} row pins {label}")
for label, tol in (("+ 1 ulp", nextup(0.1)), ("− 1 ulp", nextdown(0.1))):
    ch = st_changes("stoichiometry::empirical_formula", lambda a, t=tol: st_empirical(st_pairs(a), tol=t))
    up_ = label[0] == "+"
    k = next((kk for kk in range(1, 33) if st_changes("stoichiometry::empirical_formula", lambda a, kk=kk:
              st_empirical(st_pairs(a), tol=(0.1 + kk * 2.0 ** -56) if up_ else (0.1 - kk * 2.0 ** -56))) > 0), None)
    print(f"  EMPIRICAL_TOLERANCE {label}: {ch} rows change in Rust's replica (|r·k − round(r·k)| is a multiple of the "
          f"ratio's ulp, at least 2⁻⁵² here, so only a change across a multiple shows: the first at {k} ulps; the port's "
          f"exact 1/10 is pinned below)")
# ADR 0002 #24–#26 (M4f-2): the noted rows per function, the port's claims on the corpus, Rust's Err paths, coverage.
print("stoichiometry: ADR 0002 #24–#26: " + ", ".join(f"{st_x_noted.get(fn, 0)} of {st_x_ok.get(fn, 0)} {fn.split('::')[1]}"
      for fn in ST_EXACT) + " rows noted with the port's exact result (Fraction; every unnoted row is the port's result too)")
emp_rows = st_x_rows.get("stoichiometry::empirical_formula", [])
big_err = kept_max = 0
for els, r, port in emp_rows:
    rs = st_port_ratios(els)
    if isinstance(rs, Err): continue
    if max(rs) >= Fraction(2 ** 32) - Fraction(1, 2):
        big_err += 1
        check(port == Err("ComputationError"), f"empirical_formula{els}: a ratio of 2³² − ½ or more is not the ComputationError")
    else:
        kept_max = max(kept_max, max(rs))
        check(min(rs) == 1 and max(rs) < 2 ** 32, f"empirical_formula{els}: a kept ratio outside [1, 2³²)")
paths = {}
for els, r, port in emp_rows:
    pth = st_empirical_path(els)
    if pth in ("InvalidInput", "InvalidElement"):
        check(port == Err(pth), f"empirical_formula{els}: the port's validation is not Rust's")
        continue
    got = "Ok" if not isinstance(port, Err) else port.name
    if pth == "Ok" and got == "Ok":
        got = "Ok, the same" if toks_equal(st_toks("stoichiometry::empirical_formula", r),
                                           st_toks("stoichiometry::empirical_formula", port)) else "Ok, another result"
    key = (pth, got)
    if key not in paths: paths[key] = [0, els]
    paths[key][0] += 1
print(f"stoichiometry: ADR 0002 #26: {big_err} empirical_formula rows with a ratio of 2³² − ½ or more, every one the "
      f"ComputationError; the largest kept ratio {float(kept_max):.6g} (< 2³²); Rust's result by path → the port's:")
for (pth, got), (n, ex) in sorted(paths.items()):
    print(f"  Rust {pth} → port {got}: {n} rows, e.g. {ex}")
# the exact 1/10, pinned: 1/10 − 2⁻⁸⁰ drops the rows exactly on the edge; the smallest tolerance past 1/10 that would
# make some row pass at a smaller multiplier than now (or at all), admitted, changes that row
def st_port_dev(els, k):
    """The tolerance a row needs to pass at multiplier k: the largest |rᵢ·k − round(rᵢ·k)|."""
    rs = st_port_ratios(els)
    return None if isinstance(rs, Err) else max(abs(x * k - st_rnd(x * k)) for x in rs)
def st_port_mult(els):
    return next((k for k in range(1, 7) if st_port_dev(els, k) <= Fraction(1, 10)), 7)
beyond = min((d for els, r, port in emp_rows if not isinstance(st_port_ratios(els), Err)
              for k in range(1, st_port_mult(els)) for d in [st_port_dev(els, k)] if d > Fraction(1, 10)), default=None)
for label, tol in (("1/10 − 2⁻⁸⁰", Fraction(1, 10) - Fraction(1, 2 ** 80)),
                   (f"the nearest larger deviation, 1/10 + {float(beyond - Fraction(1, 10)):.3g}", beyond)):
    ch = sum(1 for els, r, port in emp_rows if st_toks("stoichiometry::empirical_formula", st_port_empirical(els, tol))
             != st_toks("stoichiometry::empirical_formula", port))
    print(f"  pinning the port's tolerance: with {label}, {ch} of {len(emp_rows)} empirical_formula rows change")
    check(ch > 0, f"stoichiometry: no golden empirical_formula row pins the exact 1/10 against {label}")
edge_rows = sum(1 for els, r, port in emp_rows for k in range(1, 7) if st_port_dev(els, k) == Fraction(1, 10))
check(edge_rows >= 50, f"stoichiometry: only {edge_rows} (row, k) with a deviation of exactly 1/10")
# coverage of the hard rows
py_rows = st_x_rows.get("stoichiometry::percent_yield", [])
def st_q100(a, t): return Fraction(a) * 100 / Fraction(t)
py_ties = [(a, t) for a, t, r, port in py_rows if not isinstance(port, Err) and a != 0 and is_midpoint(st_q100(a, t))]
py_sub_ties = [x for x in py_ties if st_q100(*x) < Fraction(1, 2 ** 1022)]
py_ovf = sum(1 for a, t, r, port in py_rows if not isinstance(port, Err) and a != 0
             and Fraction(2 ** 1023) <= st_q100(a, t) < Fraction(2 ** 1025))
py_unf = sum(1 for a, t, r, port in py_rows if not isinstance(r, Err) and r == 0 and a != 0)
check(len(py_ties) >= 40 and len(py_sub_ties) >= 10 and (b2f(3), 200.0) in [(a, t) for a, t in py_sub_ties],
      f"percent_yield: {len(py_ties)} tie rows, {len(py_sub_ties)} subnormal (3·2⁻¹⁰⁷⁴, 200 among them?)")
check(py_ovf >= 20 and py_unf >= 10, f"percent_yield: {py_ovf} rows at the overflow end, {py_unf} where Rust underflows to 0")
ty_rows = st_x_rows.get("stoichiometry::theoretical_yield_moles", [])
ty_mc = {"over": 0, "under": 0, "tie": 0}
for mo, co, rest, r, port in ty_rows:
    if isinstance(port, Err): continue
    k = st_port_limiting(mo, co)
    q = Fraction(mo[k]) / Fraction(co[k])
    if mo[k] != 0 and math.isfinite(port) and q > Fraction(b2f(0x7FEFFFFFFFFFFFFF)): ty_mc["over"] += 1
    if mo[k] != 0 and q < Fraction(1, 2 ** 1022) and port >= 2.0 ** -1022: ty_mc["under"] += 1
    if mo[k] != 0 and is_midpoint(q * Fraction(rest[0][1])): ty_mc["tie"] += 1
check(min(ty_mc.values()) >= 5, f"theoretical_yield_moles: hard rows {ty_mc}")
lim_rows = st_x_rows.get("stoichiometry::limiting_reagent", [])
lim_cov = {"exact tie": 0, "rounded tie": 0, "both overflow": 0, "both underflow": 0}
for mo, co, rest, r, port in lim_rows:
    if isinstance(port, Err) or len(mo) < 2: continue
    qs = [Fraction(m_) / Fraction(c_) for m_, c_ in zip(mo, co)]
    rq = [dv(m_, c_) for m_, c_ in zip(mo, co)]
    if len(set(qs)) < len(qs): lim_cov["exact tie"] += 1
    if any(rq[a_] == rq[b_] and qs[a_] != qs[b_] for a_ in range(len(qs)) for b_ in range(a_)): lim_cov["rounded tie"] += 1
    if sum(1 for x in rq if x == math.inf) >= 2: lim_cov["both overflow"] += 1
    if sum(1 for x, q in zip(rq, qs) if x == 0 and q != 0) >= 2: lim_cov["both underflow"] += 1
check(min(lim_cov.values()) >= 20, f"limiting_reagent: hard rows {lim_cov}")
print(f"stoichiometry: ADR 0002 #24–#26 hard rows: percent_yield {len(py_ties)} ties ({len(py_sub_ties)} subnormal), "
      f"{py_ovf} at the overflow end, {py_unf} where Rust underflows to 0; theoretical_yield_moles {ty_mc}; "
      f"limiting_reagent {lim_cov}; empirical_formula {edge_rows} (row, k) exactly 1/10 off")

# ── 19. The overflow audit's B1 (ADR 0002 #27, #28): thermo's, phase's and kinetics' rationals ──────────────────────────
# Every row of the ten functions is recomputed from the port's semantics, independently of the generator and the port:
# - the monomials (#27: heat_transfer, specific_heat_from_calorimetry, boiling_point_elevation, freezing_point_depression,
#   osmotic_pressure; #28: zero_order_half_life, steady_state_intermediate, pre_equilibrium_rate): the exact rational
#   (Fraction; R as its double, Rust's 1000 and 2 exact) rounded once (rn_frac), its sign the xor of the factors' signs, a
#   zero factor that signed zero; heat_transfer's and the colligative functions' non-finite inputs keep Rust's expression;
# - #28's second_order_concentration and michaelis_menten: Rust's sequence on mpfr_ref.UX (53 bits at MPFR's widest
#   exponent range, nothing over- or underflowing), the final division rounded once into binary64 (ux_div_final).
# Rust's validation is the port's (the replicas above decide the Err rows, which carry no note). Every row: noted exactly
# where the port's value differs from the oracle's (by its bits; NaN is one value), the note naming it bit for bit, so
# every unnoted row's port value is the oracle's. For the two unbounded-exponent functions, no row whose f64 chain stays
# in the normal range is noted (Rust's bits); zero_order_half_life is noted only where 2·k overflows (elsewhere Rust's
# one division of exact operands is the correctly rounded value). The hard rows are counted: 53-bit midpoints, subnormal
# results and their ties, and Rust's 0, ±∞ or NaN where the value is representable (and NaN where it is out of range:
# the correctly rounded ±∞ or 0, the policy's).
B1_HOW = {
    "thermo::heat_transfer": (27, "the correctly rounded m·c·ΔT"),
    "thermo::specific_heat_from_calorimetry": (27, "the correctly rounded q/(m·ΔT)"),
    "phase::boiling_point_elevation": (27, "the correctly rounded K·m·i"),
    "phase::freezing_point_depression": (27, "the correctly rounded K·m·i"),
    "phase::osmotic_pressure": (27, "the correctly rounded i·M·R·1000·T"),
    "kinetics::zero_order_half_life": (28, "the correctly rounded [A]₀/(2k)"),
    "kinetics::steady_state_intermediate": (28, "the correctly rounded k₁·[A]/k₂"),
    "kinetics::pre_equilibrium_rate": (28, "the correctly rounded k₂·k₁·[A]/k₋₁"),
    "kinetics::second_order_concentration": (28, "Rust's sequence evaluated with an unbounded exponent (the final division "
                                             "rounded once)"),
    "kinetics::michaelis_menten": (28, "Rust's sequence evaluated with an unbounded exponent (the final division rounded "
                                   "once)"),
}
def b1_mono(num, den):
    """RN(Π num / Π den), rounded once: the sign the xor of the factors' signs, a zero factor that signed zero."""
    neg = sum(1 for v in num + den if math.copysign(1.0, v) < 0) % 2 == 1
    if any(v == 0.0 for v in num): return -0.0 if neg else 0.0
    q = Fraction(1)
    for v in num: q *= Fraction(v)
    for v in den: q /= Fraction(v)
    r = rn_frac(q if q > 0 else -q)
    return -r if neg else r
def b1_ux_second(c0, k, t):
    with MPFR.unbounded():
        one = MPFR.UX(1.0)
        return MPFR.ux_div_final(one, one / MPFR.UX(c0) + MPFR.UX(k) * MPFR.UX(t))
def b1_ux_mm(v, km, s_):
    with MPFR.unbounded():
        return MPFR.ux_div_final(MPFR.UX(v) * MPFR.UX(s_), MPFR.UX(km) + MPFR.UX(s_))
def b1_factors(fn, a):
    """A monomial's (numerator, denominator) factor lists, in the port's order."""
    if fn == "thermo::specific_heat_from_calorimetry": return [a[0]], [a[1], a[2]]
    if fn == "phase::osmotic_pressure": return [a[0], a[1], R_GAS, 1000.0, a[2]], []
    if fn == "kinetics::zero_order_half_life": return [a[0]], [2.0, a[1]]
    if fn == "kinetics::steady_state_intermediate": return [a[0], a[1]], [a[2]]
    if fn == "kinetics::pre_equilibrium_rate": return [a[2], a[0], a[3]], [a[1]]
    return list(a), []
def b1_port(fn, a):
    rust = RULES[fn][1]
    if fn in ("thermo::heat_transfer", "phase::boiling_point_elevation", "phase::freezing_point_depression"):
        return b1_mono(list(a), []) if all(math.isfinite(v) for v in a) else a[0] * a[1] * a[2]
    r = rust([("f", v, f2b(v)) for v in a])
    if isinstance(r, Err): return r
    if fn == "kinetics::second_order_concentration": return b1_ux_second(*a)
    if fn == "kinetics::michaelis_menten": return b1_ux_mm(*a)
    return b1_mono(*b1_factors(fn, a))
# Rust's f64 chains of the two unbounded-exponent functions, traced (OOR as #20's): an operation's result non-finite, a
# nonzero subnormal, or 0 from nonzero factors leaves the normal range. The final division's own range does not count:
# it is rounded once into binary64 in both.
def b1_in_range(fn, a):
    OOR[0] = False
    if fn == "kinetics::second_order_concentration":
        c0, k, t = a
        ta(tq(1.0, c0), tm(k, t))
    else:
        v, km, s_ = a
        tm(v, s_)
        ta(km, s_)
    return not OOR[0]
B1_RE = re.compile(r"\[ADR 0002 #(2[78]): the port returns (.+) = (\S+) \(f:([0-9a-f]{16})\)\]$")
b1_n = {fn: [0, 0, 0, 0] for fn in B1_HOW}          # Ok rows, noted, in range (the #20-style), qualitative
b1_cov = {"53-bit midpoint": 0, "subnormal result": 0, "subnormal tie": 0, "Rust 0, ±∞ or NaN, the value representable": 0,
          "Rust's special differs, the value out of range": 0}
def b1_midpoint(q):
    """1 if the exact positive rational q is a midpoint of binary64's grid (53-bit or subnormal), else 0."""
    if q <= 0: return False
    e = q.numerator.bit_length() - q.denominator.bit_length()
    if q < Fraction(2) ** e: e -= 1                # q in [2^e, 2^(e + 1))
    u = q / Fraction(2) ** (max(e, -1022) - 53)    # in units of half the ulp (or of the subnormal grid's)
    return u.denominator == 1 and u.numerator % 2 == 1
for mod_ in ("thermo", "phase", "kinetics"):
    for fn, cls, ins, outs, where in R.get(mod_, []):
        if fn not in B1_HOW: continue
        adr, how = B1_HOW[fn]
        a = fl(ins)
        port = b1_port(fn, a)
        note = NOTES.get(where, "")
        noted = f"[ADR 0002 #{adr}:" in note
        pt = out_toks(port)
        differs = not toks_equal(pt, outs)
        check(noted == differs, f"{where}: {fn}{tuple(a)}: the port's {pt} vs the oracle's {outs}: the #{adr} note is "
              f"{'missing' if differs else 'stray'}")
        if isinstance(port, Err): continue
        b1_n[fn][0] += 1
        if noted:
            b1_n[fn][1] += 1
            mm_ = B1_RE.search(note)
            check(mm_ is not None and int(mm_.group(1)) == adr and mm_.group(2) == how and int(mm_.group(4), 16) == f2b(port)
                  and same_bits(float(mm_.group(3)), port), f"{where}: the #{adr} note does not name {how} = {port!r}")
        if fn in ("kinetics::second_order_concentration", "kinetics::michaelis_menten"):
            inr = b1_in_range(fn, a)
            b1_n[fn][2] += inr
            check(not (inr and noted), f"{where}: {fn}{tuple(a)}: Rust's chain stays in the normal range, so its bits are "
                  f"the port's, but the row is noted")
        if fn == "kinetics::zero_order_half_life" and 2.0 * a[1] < math.inf:
            check(not noted, f"{where}: zero_order_half_life{tuple(a)} noted where 2k is finite")
        o = outs[0][1] if outs[0][0] == "f" else None
        if o is not None and (o == 0.0 or not math.isfinite(o)) and not same_bits(o, port):
            if math.isfinite(port) and (port != 0.0 or any(v == 0.0 for v in a)):
                b1_n[fn][3] += 1                   # the value is representable: a nonzero double, or an exact 0
                b1_cov["Rust 0, ±∞ or NaN, the value representable"] += 1
            else:
                b1_cov["Rust's special differs, the value out of range"] += 1
        if fn in ("kinetics::second_order_concentration", "kinetics::michaelis_menten") or not all(math.isfinite(v) for v in a):
            continue
        num, den = b1_factors(fn, a)
        if any(v == 0.0 for v in num): continue
        q = Fraction(1)
        for v in num: q *= abs(Fraction(v))
        for v in den: q /= abs(Fraction(v))
        if q < Fraction(2) ** -1022:
            b1_cov["subnormal result"] += 1
            b1_cov["subnormal tie"] += b1_midpoint(q)
        elif b1_midpoint(q):
            b1_cov["53-bit midpoint"] += 1
for fn, (ok_, nt, inr, ql) in b1_n.items():
    extra = f", {inr} with Rust's chain in the normal range (all unnoted, Rust's bits)" if fn in (
        "kinetics::second_order_concentration", "kinetics::michaelis_menten") else ""
    print(f"{fn.split('::')[0]}: ADR 0002 #{B1_HOW[fn][0]}: {fn.split('::')[1]}: {nt} of {ok_} Ok rows noted with the port's "
          f"value ({ql} where Rust's is 0, ±∞ or NaN and the value representable){extra}")
    check(nt > 0 and ql > 0, f"{fn}: no noted row, or no qualitative edge row")
print(f"ADR 0002 #27, #28 hard rows: {b1_cov}")
check(all(v > 0 for k, v in b1_cov.items() if "out of range" not in k), f"B1: a hard-row class has no golden row: {b1_cov}")

# ── 20. The overflow audit's B2 (ADR 0002 #29, #30): solution's, electrochemistry's and spectroscopy's functions ──────────
# Every row of the twelve functions is recomputed from the port's semantics, independently of the generator and the port:
# - the monomials (#29: dilution, solution_conductivity, faraday_mass_deposited, faraday_charge_required,
#   faraday_mass_ratio; #30: absorbance, concentration_from_absorbance, rigid_rotor_energy, rotational_constant,
#   diatomic_moment_of_inertia): the exact rational (Fraction; F, h and π as their doubles, Rust's 1000, 8, n and J + 1
#   exact) rounded once (b1_mono), the unvalidated inputs' non-finite values keeping Rust's expression;
# - #30's vibrational_frequency: RN(√(k/μ)) of the exact quotient, by Python's integer square root (b2_sqrt, exact) and
#   again by MPFR (mpfr_ref.ux_sqrt_div_final, a 256-bit quotient: its argument needs the gap this section proves, that
#   √(k/μ) is never a rounding boundary of binary64);
# - #29's nernst_potential: Rust's sequence on mpfr_ref.UX (53 bits, MPFR's widest exponent range) with ln Q MPFR's
#   correctly rounded log (the port's _km_cr_ln), the final subtraction rounded once into binary64 (ux_sub_final).
# Rust's validation is replicated here and must agree with the oracle (Err rows carry no note). Every row: noted exactly
# where the port's value differs from the oracle's (by its bits; NaN is one value), the note naming it bit for bit, so
# every unnoted row's port value is the oracle's. nernst_potential: no row whose f64 chain stays in the normal range is
# noted (Rust's bits), and every nernst_potential_25c row's chain stays in it (so the 25 °C function, which calls
# nernst_potential at 298.15 K, is unchanged). The hard rows are counted: 53-bit midpoints, subnormal results and their
# ties, Rust's 0, ±∞ or NaN where the value is representable, √'s exact roots and its quotients past either end.
B2_HOW = {
    "solution::dilution": (29, "the correctly rounded M₁·V₁/V₂"),
    "solution::solution_conductivity": (29, "the correctly rounded Λ·c/1000"),
    "electrochemistry::faraday_mass_deposited": (29, "the correctly rounded Q·M/(n·F)"),
    "electrochemistry::faraday_charge_required": (29, "the correctly rounded m·n·F/M"),
    "electrochemistry::faraday_mass_ratio": (29, "the correctly rounded m_A·M_B·n_A/(n_B·M_A)"),
    "electrochemistry::nernst_potential": (29, "Rust's sequence evaluated with an unbounded exponent (the final "
                                           "subtraction rounded once)"),
    "spectroscopy::absorbance": (30, "the correctly rounded ε·l·c"),
    "spectroscopy::concentration_from_absorbance": (30, "the correctly rounded A/(ε·l)"),
    "spectroscopy::vibrational_frequency": (30, "the correctly rounded √(k/μ)"),
    "spectroscopy::rigid_rotor_energy": (30, "the correctly rounded h·B·J·(J + 1)"),
    "spectroscopy::rotational_constant": (30, "the correctly rounded h/(8·π·π·I)"),
    "spectroscopy::diatomic_moment_of_inertia": (30, "the correctly rounded μ·r·r"),
}
def b2_sqrt(k, mu):
    """RN(√(k/μ)) for positive finite doubles, exactly: Q = ⌊N·2^s/D⌋ with at least 120 bits and the exponent even, R =
    ⌊√Q⌋ (math.isqrt), and RN((2R + sticky)·2^(e − 1)) by CPython's correctly rounded int true division."""
    q = Fraction(k) / Fraction(mu)
    n, d = q.numerator, q.denominator
    s = max(0, 121 + d.bit_length() - n.bit_length())
    if s % 2: s += 1
    quo, rem = divmod(n << s, d)
    r = math.isqrt(quo)
    st = 1 if (r * r != quo or rem) else 0
    return rn_frac(Fraction(2 * r + st, 1 << (s // 2 + 1)))
def b2_no_boundary(k, mu):
    """True if √(k/μ) is no rounding boundary of binary64, exactly: a boundary is dyadic, so k/μ must be the square of a
    dyadic r = M·2^f (M odd), and r a midpoint: M of 54 bits on the 53-bit grid (r ≥ 2⁻¹⁰²²), f = −1075 on the subnormal
    one. (The general argument: k/μ = (p/q)·2^E with p, q odd needs q = 1 and p = M² < 2⁵³, so M has at most 26 bits,
    and E = 2f ≥ −2097 > −2150.)"""
    q = Fraction(k) / Fraction(mu)
    n, d = q.numerator, q.denominator
    if d & (d - 1): return True                     # not dyadic
    j = d.bit_length() - 1
    r = math.isqrt(n)
    if r * r != n or j % 2: return True             # not a dyadic square
    mm, t = r, 0
    while mm % 2 == 0:
        mm //= 2
        t += 1
    if Fraction(r, 1 << (j // 2)) >= Fraction(2) ** -1022: return mm.bit_length() != 54
    return t - j // 2 != -1075
DMAX_V = b2f(0x7FEFFFFFFFFFFFFF)
def b2_nernst(e0, n, t, q):
    with MPFR.unbounded():
        a = MPFR.UX(R_GAS) * MPFR.UX(t)
        c = a / (MPFR.UX(float(n)) * MPFR.UX(FARADAY))
        p_ = c * MPFR.UX(MPFR.log(q))
        return MPFR.ux_sub_final(MPFR.UX(e0), p_)
def b2_nernst_in_range(t, n, l):
    OOR[0] = False
    a = tm(R_GAS, t)
    c = tq(a, float(n) * FARADAY)
    tm(c, l)
    return not OOR[0]
def b2_factors(fn, a):
    """A monomial's (numerator, denominator) factor lists, in the port's order (rotational_constant's 8 included)."""
    if fn == "solution::dilution": return [a[0], a[1]], [a[2]]
    if fn == "solution::solution_conductivity": return [a[0], a[1]], [1000.0]
    if fn == "electrochemistry::faraday_mass_deposited": return [a[0], a[1]], [float(a[2]), FARADAY]
    if fn == "electrochemistry::faraday_charge_required": return [a[0], float(a[2]), FARADAY], [a[1]]
    if fn == "electrochemistry::faraday_mass_ratio": return [a[0], a[3], float(a[2])], [float(a[4]), a[1]]
    if fn == "spectroscopy::concentration_from_absorbance": return [a[0]], [a[1], a[2]]
    if fn == "spectroscopy::rigid_rotor_energy": return [PLANCK, a[1], float(a[0]), float(a[0]) + 1.0], []
    if fn == "spectroscopy::rotational_constant": return [PLANCK], [8.0, math.pi, math.pi, a[0]]
    if fn == "spectroscopy::diatomic_moment_of_inertia": return [a[0], a[1], a[1]], []
    return list(a), []                              # absorbance
def b2_valid(fn, a):
    """Rust's validation (None: Ok), in Rust's order."""
    if fn == "solution::dilution":
        return None if is_positive(a[2]) and is_finite(a[0]) and is_finite(a[1]) else Err("InvalidInput")
    if fn in ("electrochemistry::faraday_mass_deposited", "electrochemistry::faraday_charge_required"):
        return None if a[2] != 0 and is_positive(a[1]) and is_finite(a[0]) else Err("InvalidInput")
    if fn == "electrochemistry::faraday_mass_ratio":
        ok = a[2] != 0 and a[4] != 0 and is_positive(a[1]) and is_positive(a[3]) and is_finite(a[0])
        return None if ok else Err("InvalidInput")
    if fn == "electrochemistry::nernst_potential":
        if not is_positive(a[2]): return Err("InvalidTemperature")
        return None if a[1] != 0 and is_positive(a[3]) and is_finite(a[0]) else Err("InvalidInput")
    if fn == "spectroscopy::concentration_from_absorbance":
        return None if is_positive(a[1]) and is_positive(a[2]) else Err("InvalidInput")
    if fn == "spectroscopy::vibrational_frequency":
        return None if is_positive(a[0]) and is_positive(a[1]) else Err("InvalidInput")
    if fn == "spectroscopy::rotational_constant":
        return None if is_positive(a[0]) else Err("InvalidInput")
    return None
def b2_port(fn, a):
    e = b2_valid(fn, a)
    if e is not None: return e
    if fn == "electrochemistry::nernst_potential": return b2_nernst(*a)
    if fn == "spectroscopy::vibrational_frequency": return b2_sqrt(a[0], a[1])
    if fn == "solution::solution_conductivity" and not all(math.isfinite(v) for v in a):
        return dv(a[0] * a[1], 1000.0)
    if fn in ("spectroscopy::absorbance", "spectroscopy::diatomic_moment_of_inertia") and not all(
            math.isfinite(v) for v in a):
        return a[0] * a[1] * (a[2] if fn == "spectroscopy::absorbance" else a[1])
    if fn == "spectroscopy::concentration_from_absorbance" and not math.isfinite(a[0]):
        return dv(a[0], a[1] * a[2])
    if fn == "spectroscopy::rigid_rotor_energy" and not math.isfinite(a[1]):
        return PLANCK * a[1] * float(a[0]) * (float(a[0]) + 1.0)
    return b1_mono(*b2_factors(fn, a))
B2_RE = re.compile(r"\[ADR 0002 #(29|30): the port returns (.+) = (\S+) \(f:([0-9a-f]{16})\)\]$")
b2_n = {fn: [0, 0, 0, 0] for fn in B2_HOW}          # Ok rows, noted, in range (nernst), qualitative
b2_cov = {"53-bit midpoint": 0, "subnormal result": 0, "subnormal tie": 0, "Rust 0, ±∞ or NaN, the value representable": 0,
          "Rust's special differs, the value out of range": 0, "√ exact": 0, "√: k/μ past DBL_MAX": 0,
          "√: k/μ below 2⁻¹⁰⁷⁴": 0, "√: a subnormal operand": 0, "√: subnormal result": 0,
          "√ within 2⁻⁹⁸ of a 53-bit midpoint": 0}
b2_mpfr = [0, 0]                                     # √ rows checked against MPFR; Nernst rows whose glibc ln misrounds
for mod_ in ("solution", "electrochemistry", "spectroscopy"):
    for fn, cls, ins, outs, where in R.get(mod_, []):
        if fn == "electrochemistry::nernst_potential_25c":
            a = fl(ins)
            if b2_valid("electrochemistry::nernst_potential", [a[0], a[1], 298.15, a[2]]) is None:
                check(b2_nernst_in_range(298.15, a[1], math.log(a[2])), f"{where}: nernst_potential_25c{tuple(a)} leaves "
                      f"the normal range")
            continue
        if fn not in B2_HOW: continue
        adr, how = B2_HOW[fn]
        a = fl(ins)
        port = b2_port(fn, a)
        note = NOTES.get(where, "")
        noted = f"[ADR 0002 #{adr}:" in note
        pt = out_toks(port)
        differs = not toks_equal(pt, outs)
        check(noted == differs, f"{where}: {fn}{tuple(a)}: the port's {pt} vs the oracle's {outs}: the #{adr} note is "
              f"{'missing' if differs else 'stray'}")
        if isinstance(port, Err): continue
        b2_n[fn][0] += 1
        if noted:
            b2_n[fn][1] += 1
            mm_ = B2_RE.search(note)
            check(mm_ is not None and int(mm_.group(1)) == adr and mm_.group(2) == how and int(mm_.group(4), 16) == f2b(port)
                  and same_bits(float(mm_.group(3)), port), f"{where}: the #{adr} note does not name {how} = {port!r}")
        if fn == "electrochemistry::nernst_potential":
            lg = MPFR.log(a[3])
            if not same_bits(lg, math.log(a[3])): b2_mpfr[1] += 1
            inr = b2_nernst_in_range(a[2], a[1], lg)
            b2_n[fn][2] += inr
            check(not (inr and noted), f"{where}: {fn}{tuple(a)}: Rust's chain stays in the normal range, so its bits are "
                  f"the port's, but the row is noted")
        o = outs[0][1] if outs[0][0] == "f" else None
        if o is not None and (o == 0.0 or not math.isfinite(o)) and not same_bits(o, port):
            zero_in = any(v == 0.0 for v in a if isinstance(v, float))
            if math.isfinite(port) and (port != 0.0 or zero_in):
                b2_n[fn][3] += 1                   # the value is representable: a nonzero double, or an exact 0
                b2_cov["Rust 0, ±∞ or NaN, the value representable"] += 1
            else:
                b2_cov["Rust's special differs, the value out of range"] += 1
        if fn == "spectroscopy::vibrational_frequency":
            k, mu = a
            check(b2_no_boundary(k, mu), f"{where}: √({k!r}/{mu!r}) is a rounding boundary")
            wm = MPFR.ux_sqrt_div_final(MPFR.UX(k), MPFR.UX(mu))
            check(same_bits(wm, port), f"{where}: MPFR's √(k/μ) {wm!r} vs the integer root's {port!r}")
            b2_mpfr[0] += 1
            q = Fraction(k) / Fraction(mu)
            if math.isfinite(port) and Fraction(port) ** 2 == q: b2_cov["√ exact"] += 1
            if q > Fraction(DMAX_V): b2_cov["√: k/μ past DBL_MAX"] += 1
            if q < Fraction(2) ** -1074: b2_cov["√: k/μ below 2⁻¹⁰⁷⁴"] += 1
            if min(k, mu) < 2.0 ** -1022: b2_cov["√: a subnormal operand"] += 1
            if port < 2.0 ** -1022: b2_cov["√: subnormal result"] += 1
            if 2.0 ** -1022 <= port < math.inf:
                h = Fraction(math.ulp(port)) / 2
                if any(abs(q / (Fraction(port) + d) ** 2 - 1) < Fraction(1, 1 << 97) for d in (h, -h)):
                    b2_cov["√ within 2⁻⁹⁸ of a 53-bit midpoint"] += 1
            continue
        if fn == "electrochemistry::nernst_potential": continue
        num, den = b2_factors(fn, a)
        if not all(math.isfinite(v) for v in num + den) or any(v == 0.0 for v in num): continue
        q = Fraction(1)
        for v in num: q *= abs(Fraction(v))
        for v in den: q /= abs(Fraction(v))
        if q < Fraction(2) ** -1022:
            b2_cov["subnormal result"] += 1
            b2_cov["subnormal tie"] += b1_midpoint(q)
        elif b1_midpoint(q):
            b2_cov["53-bit midpoint"] += 1
for fn, (ok_, nt, inr, ql) in b2_n.items():
    extra = f", {inr} with Rust's chain in the normal range (all unnoted, Rust's bits)" if fn == (
        "electrochemistry::nernst_potential") else ""
    print(f"{fn.split('::')[0]}: ADR 0002 #{B2_HOW[fn][0]}: {fn.split('::')[1]}: {nt} of {ok_} Ok rows noted with the port's "
          f"value ({ql} where Rust's is 0, ±∞ or NaN and the value representable){extra}")
    check(nt > 0, f"{fn}: no noted row")
print(f"ADR 0002 #29, #30 hard rows: {b2_cov}")
print(f"ADR 0002 #30: {b2_mpfr[0]} vibrational_frequency rows checked against MPFR's √ of a 256-bit quotient as well; "
      f"ADR 0002 #29: {b2_mpfr[1]} nernst_potential rows where glibc's ln is not correctly rounded")
check(all(v > 0 for k, v in b2_cov.items() if "out of range" not in k), f"B2: a hard-row class has no golden row: {b2_cov}")
check(b2_mpfr[1] == 0, "B2: a nernst_potential row whose glibc ln misrounds (its note, from glibc's ln, is not the port's)")

# ── 21. The overflow audit's B3 (ADR 0002 #31, #32): gas's, potential's, inorganic's and nuclear's functions ─────────────
# Every row of the twelve functions is recomputed from the port's semantics, independently of the generator and the port:
# - the monomials (#31: ideal_gas_pressure, ideal_gas_volume, gas_density; #32: coulomb, coulomb_force): the exact rational
#   (Fraction; R and k_e as their doubles) rounded once (b1_mono), the unvalidated charges' non-finite values keeping
#   Rust's expression;
# - #32's morse_width_parameter: RN(ω·√(μ/(2·D_e))) of the exact value, by Python's integer square root (b3_sqrt_frac,
#   exact, ties to even: the value can be a midpoint) and again by MPFR (mpfr_ref.ux_mul_sqrt_div_final, a 320-bit
#   quotient: exact where it is dyadic, else at least 2^-162 from every boundary's square);
# - the #20-style ones (#31: van_der_waals_pressure, peng_robinson_pressure, redlich_kwong_pressure; #32: lennard_jones,
#   born_lande_lattice_energy, binding_energy_per_nucleon): Rust's sequence on mpfr_ref.UX (53 bits, MPFR's widest
#   exponent range), Redlich–Kwong's Tc^2.5 by MPFR's own pow at 53 bits there (UX.pow: not the port's rescaling), the
#   co-volume tests V ≤ n·b and V ≤ b on the UX values, the final operation rounded once into binary64.
# Rust's validation is replicated here (Err rows that both refuse carry no note). Every row: noted exactly where the
# port's result differs from the oracle's (by its bits, NaN one value; an Err against a value), the note naming the port's
# value bit for bit, or its Err (`= Err(InvalidInput) (E:InvalidInput)`), so every unnoted row's port result is the
# oracle's. For the #20-style functions no row whose f64 chain stays in the normal range is noted (Rust's bits): the
# tracers below replay Rust's chain, and a product, quotient or root leaves the range when it is not finite, is 0 from
# nonzero operands, or is nonzero and at most 2⁻¹⁰²² in magnitude (2⁻¹⁰²² itself included: a value just below it
# rounds there at the subnormal grid but not at 53 bits), a sum when it is not finite (a subnormal sum of two doubles is
# exact); the final operation's own range does not count (it is rounded once in both). The folds are proven: #17's
# (ε·d)·4 is the port's value on every lennard_jones row with ε ≥ 2¹⁰²² whose d chain is normal (bit-identical), and
# the #18 rows (D_e ≥ 2¹⁰²³) are counted by whether #18's ω·√(RN(μ/(2·D_e))) is the new correctly rounded value.
# Redlich–Kwong's exact rescale is checked on every row: RN₅₃(Tc^2.5) = RN(s^2.5)·2^(5j) for Tc = s·2^(2j), s in [1, 4)
# (MPFR's pow both sides). The hard rows are counted.
B3_SUB = "Rust's sequence evaluated with an unbounded exponent (the final subtraction rounded once)"
B3_MUL = "Rust's sequence evaluated with an unbounded exponent (the final product rounded once)"
B3_DIV = "Rust's sequence evaluated with an unbounded exponent (the final division rounded once)"
B3_HOW = {
    "gas::ideal_gas_pressure": (31, "the correctly rounded n·R·T/V"),
    "gas::ideal_gas_volume": (31, "the correctly rounded n·R·T/P"),
    "gas::gas_density": (31, "the correctly rounded P·M/(R·T)"),
    "gas::van_der_waals_pressure": (31, B3_SUB),
    "gas::peng_robinson_pressure": (31, B3_SUB),
    "gas::redlich_kwong_pressure": (31, B3_SUB),
    "potential::lennard_jones": (32, B3_MUL),
    "potential::morse_width_parameter": (32, "the correctly rounded ω·√(μ/(2·D_e))"),
    "potential::coulomb": (32, "the correctly rounded k_e·q₁·q₂/r"),
    "potential::coulomb_force": (32, "the correctly rounded k_e·q₁·q₂/(r·r)"),
    "inorganic::born_lande_lattice_energy": (32, B3_MUL),
    "nuclear::binding_energy_per_nucleon": (32, B3_DIV),
}
def p_pow25_glibc(x):
    """The platform libm's pow(x, 2.5) (glibc's on the corpus's platform: the oracle's tc.powf(2.5))."""
    try: return math.pow(x, 2.5)
    except OverflowError: return math.inf
B3_UX = ("gas::van_der_waals_pressure", "gas::peng_robinson_pressure", "gas::redlich_kwong_pressure",
         "potential::lennard_jones", "inorganic::born_lande_lattice_energy", "nuclear::binding_energy_per_nucleon")
INVT = Err("InvalidTemperature")
def b3_sqrt_frac(q):
    """RN(√q) for a positive Fraction q, exactly (ties to even): Q = ⌊N·2^s/D⌋ with at least 120 bits and an even
    exponent, R = ⌊√Q⌋, and RN((2R + sticky)·2^-(s/2 + 1)) by CPython's correctly rounded int true division (a tie only
    with sticky = 0, when √q is the exact midpoint)."""
    n, d = q.numerator, q.denominator
    s = 121 + d.bit_length() - n.bit_length()
    if s % 2: s += 1
    quo, rem = divmod(n << s, d) if s >= 0 else divmod(n, d << -s)
    r = math.isqrt(quo)
    st = 1 if (r * r != quo or rem) else 0
    return rn_frac(Fraction(2 * r + st) * Fraction(2) ** (-(s // 2) - 1))
def b3_width(w, mu, d): return b3_sqrt_frac(Fraction(w) ** 2 * Fraction(mu) / (2 * Fraction(d)))
def b3_valid(fn, a):
    """Rust's validation, in Rust's order (None: it passes; the co-volume tests are the port functions')."""
    if fn in ("gas::ideal_gas_pressure", "gas::ideal_gas_volume"):
        if not is_finite(a[0]): return INVALID
        if not is_positive(a[1]): return INVT
        return None if is_positive(a[2]) else INVALID
    if fn == "gas::gas_density":
        if not is_positive(a[2]): return INVT
        return None if is_non_negative(a[0]) and is_positive(a[1]) else INVALID
    if fn == "gas::van_der_waals_pressure":
        if not is_positive(a[1]): return INVT
        ok = is_finite(a[0]) and is_non_negative(a[3]) and is_non_negative(a[4]) and is_positive(a[2])
        return None if ok else INVALID
    if fn == "gas::peng_robinson_pressure":
        if not is_positive(a[3]) or not is_positive(a[0]): return INVT
        return None if is_positive(a[1]) and is_finite(a[2]) and is_finite(a[4]) else INVALID
    if fn == "gas::redlich_kwong_pressure":
        if not is_positive(a[2]) or not is_positive(a[0]): return INVT
        return None if is_positive(a[1]) and is_finite(a[3]) else INVALID
    if fn == "potential::lennard_jones":
        return None if all(is_positive(v) for v in a) else INVALID
    if fn == "potential::morse_width_parameter":
        return None if all(is_positive(v) for v in a) else INVALID
    if fn in ("potential::coulomb", "potential::coulomb_force"):
        return None if is_positive(a[2]) else INVALID
    if fn == "inorganic::born_lande_lattice_energy":
        mm, zc, za, r0, n = a
        ok = is_positive(r0) and is_finite(n) and n > 1.0 and is_positive(mm) and zc != 0 and za != 0
        return None if ok else INVALID
    z, n, m = a                                          # binding_energy_per_nucleon
    return None if z + n <= 0xFFFFFFFF and z + n != 0 and is_positive(m) else INVALID
U = MPFR.UX
def b3_pr_coef(tc, pc, om, t):
    """Rust's pr_coefficients on UX: (a·α, b)."""
    r, tcu, pcu = U(R_GAS), U(tc), U(pc)
    a = U(0.45724) * r * r * tcu * tcu / pcu
    b = U(0.07780) * r * tcu / pcu
    o = U(om)
    kappa = (U(0.37464) + U(1.54226) * o) - U(0.26992) * o * o
    u = U(1.0) + kappa * (U(1.0) - (U(t) / tcu).sqrt())
    return a * (u * u), b
def b3_port(fn, a):
    e = b3_valid(fn, a)
    if e is not None: return e
    if fn in ("gas::ideal_gas_pressure", "gas::ideal_gas_volume"): return b1_mono([a[0], R_GAS, a[1]], [a[2]])
    if fn == "gas::gas_density": return b1_mono([a[0], a[1]], [R_GAS, a[2]])
    if fn in ("potential::coulomb", "potential::coulomb_force"):
        q1, q2, r = a
        den = [r] if fn == "potential::coulomb" else [r, r]
        if not (is_finite(q1) and is_finite(q2)): return dv(KE_POT * q1 * q2, r if len(den) == 1 else r * r)
        return b1_mono([KE_POT, q1, q2], den)
    if fn == "potential::morse_width_parameter": return b3_width(*a)
    with MPFR.unbounded():
        if fn == "gas::van_der_waals_pressure":
            n, t, v, aa, b = a
            nb = U(n) * U(b)
            if U(v) <= nb: return INVALID
            nv = U(n) / U(v)
            return MPFR.ux_sub_final(U(n) * U(R_GAS) * U(t) / (U(v) - nb), U(aa) * nv * nv)
        if fn == "gas::peng_robinson_pressure":
            tc, pc, om, t, v = a
            aal, b = b3_pr_coef(tc, pc, om, t)
            if U(v) <= b: return INVALID
            vu = U(v)
            den = (vu * vu + U(2.0) * b * vu) - b * b
            return MPFR.ux_sub_final(U(R_GAS) * U(t) / (vu - b), aal / den)
        if fn == "gas::redlich_kwong_pressure":
            tc, pc, t, v = a
            r, tcu, pcu = U(R_GAS), U(tc), U(pc)
            aa = U(0.42748) * r * r * tcu.pow(2.5) / pcu
            b = U(0.08664) * r * tcu / pcu
            if U(v) <= b: return INVALID
            vu = U(v)
            return MPFR.ux_sub_final(r * U(t) / (vu - b), aa / (U(t).sqrt() * vu * (vu + b)))
        if fn == "potential::lennard_jones":
            e_, s, r = a
            sr = U(s) / U(r)
            sr6 = sr * sr * sr * sr * sr * sr
            return MPFR.ux_mul_final(U(4.0) * U(e_), sr6 * sr6 - sr6)
        if fn == "inorganic::born_lande_lattice_energy":
            mm, zc, za, r0, n = a
            u = U(-KE) * U(NA) * U(mm) * U(float(zc)) * U(float(za)) * U(QE) * U(QE) / U(r0)
            return MPFR.ux_mul_final(u, U(1.0) - U(1.0) / U(n))
        z, n, m = a
        md = U(float(z)) * U(M_H) + U(float(n)) * U(M_N) - U(m)
        return MPFR.ux_div_final(md * U(AMU_MEV), U(float(z + n)))
# Rust's f64 chains of the #20-style functions, traced (B3's strict range: OOR[0] set as the section comment says).
def t3m(r, *fs):
    if not math.isfinite(r) or (r == 0.0 and fs and all(x != 0.0 for x in fs)) or (r != 0.0 and abs(r) <= 2.0 ** -1022):
        OOR[0] = True
    return r
def t3a(r):
    if not math.isfinite(r): OOR[0] = True
    return r
def b3_in_range(fn, a):
    """True when Rust's f64 chain of the row (validated) stays in the normal range, so Rust's bits are the port's."""
    OOR[0] = False
    if fn == "gas::van_der_waals_pressure":
        n, t, v, aa, b = a
        nb = t3m(n * b, n, b)
        nr = t3m(n * R_GAS, n, R_GAS)
        nrt = t3m(nr * t, nr, t)
        t3m(dv(nrt, t3a(v - nb)), nrt)
        nv = t3m(dv(n, v), n)
        an = t3m(aa * nv, aa, nv)
        t3m(an * nv, an, nv)
        return not OOR[0]
    if fn == "gas::peng_robinson_pressure":
        tc, pc, om, t, v = a
        r = R_GAS
        c4 = t3m(t3m(0.45724 * r * r * tc, tc) * tc, tc)
        aa = t3m(dv(c4, pc), c4)
        b = t3m(dv(t3m(0.07780 * r * tc, tc), pc), tc)
        k1 = t3m(1.54226 * om, om)
        k4 = t3m(t3m(0.26992 * om, om) * om, om)
        kappa = t3a(t3a(0.37464 + k1) - k4)
        tr = t3m(dv(t, tc), t)
        s1 = t3a(1.0 - math.sqrt(tr))
        s2 = t3m(kappa * s1, kappa, s1)
        u = t3a(1.0 + s2)
        al = t3m(u * u, u)
        aal = t3m(aa * al, aa, al)
        if not OOR[0] and v <= b: return True             # the Err: decided on a b that is the unbounded one
        t3m(dv(t3m(r * t, t), t3a(v - b)), t)
        vv = t3m(v * v, v)
        b2v = t3m(t3m(2.0 * b, b) * v, b, v)
        den = t3a(t3a(vv + b2v) - t3m(b * b, b))
        t3m(dv(aal, den), aal)
        return not OOR[0]
    if fn == "gas::redlich_kwong_pressure":
        tc, pc, t, v = a
        r = R_GAS
        aa = t3m(dv(t3m(0.42748 * r * r * t3m(MPFR.pow(tc, 2.5), tc), tc), pc), tc)
        b = t3m(dv(t3m(0.08664 * r * tc, tc), pc), tc)
        if not OOR[0] and v <= b: return True
        t3m(dv(t3m(r * t, t), t3a(v - b)), t)
        t3m(dv(aa, t3m(t3m(math.sqrt(t) * v, v) * t3a(v + b), v)), aa)
        return not OOR[0]
    if fn == "potential::lennard_jones":
        e_, s, r = a
        sr = t3m(dv(s, r), s)
        p = sr
        for _ in range(5): p = t3m(p * sr, p)
        t3a(t3m(p * p, p) - p)
        t3m(4.0 * e_, e_)                                 # Rust's 4·ε (#17's overflow)
        return not OOR[0]
    if fn == "inorganic::born_lande_lattice_energy":
        mm, zc, za, r0, n = a
        u = -KE * NA
        for f_ in (mm, float(zc), float(za), QE, QE): u = t3m(u * f_, u)
        t3m(dv(u, r0), u)
        t3a(1.0 - t3m(dv(1.0, n), 1.0))
        return not OOR[0]
    z, n, m = a
    md = t3a(t3a(t3m(float(z) * M_H, float(z)) + t3m(float(n) * M_N, float(n))) - m)
    t3m(md * AMU_MEV, md)
    return not OOR[0]
B3_RE = re.compile(r"\[ADR 0002 #(31|32): the port returns (.+) = (\S+) \((f|E):([0-9A-Za-z]+)\)\]$")
b3_n = {fn: [0, 0, 0, 0, 0] for fn in B3_HOW}        # rows (non-Err for either), noted, in range, qualitative, Err flips
b3_cov = {"53-bit midpoint": 0, "subnormal result": 0, "subnormal tie": 0, "Rust 0, ±∞ or NaN, the value representable": 0,
          "Rust's special differs, the value out of range": 0, "Rust's NaN, the value ±∞": 0,
          "the port's Err, Rust's value (co-volume)": 0, "the port's value, Rust's Err (co-volume)": 0,
          "√ exact": 0, "√ a midpoint (tie)": 0, "√: μ/(2·D_e) past DBL_MAX": 0, "√: μ/(2·D_e) below 2⁻¹⁰⁷⁴": 0,
          "√: a subnormal operand": 0, "√: subnormal result": 0, "√ within 2⁻⁹⁸ of a 53-bit midpoint": 0,
          "Tc^2.5 past either end": 0}
b3_fold = {"#17 rows (ε ≥ 2¹⁰²²)": 0, "#17 rows with d normal (proven: (ε·d)·4 the port's)": 0,
           "#17's (ε·d)·4 the port's": 0, "#17's (ε·d)·4 the port's, Rust's differs (#17's notes)": 0,
           "#17's (ε·d)·4 not the port's (d out of range)": 0,
           "#18 rows (D_e ≥ 2¹⁰²³)": 0, "#18 notes: the value kept": 0, "#18 notes: one ulp off": 0,
           "#18 notes: further off": 0, "#18's value Rust's 0, the new one representable": 0,
           "Tc^2.5 rescale checked": 0, "RK rows whose glibc pow(Tc, 2.5) misrounds": 0, "√ checked against MPFR": 0}
for mod_ in ("gas", "potential", "inorganic", "nuclear"):
    for fn, cls, ins, outs, where in R.get(mod_, []):
        if fn not in B3_HOW: continue
        adr, how = B3_HOW[fn]
        a = fl(ins)
        port = b3_port(fn, a)
        note = NOTES.get(where, "")
        noted = f"[ADR 0002 #{adr}:" in note
        pt = out_toks(port)
        differs = not toks_equal(pt, outs)
        check(noted == differs, f"{where}: {fn}{tuple(a)}: the port's {pt} vs the oracle's {outs}: the #{adr} note is "
              f"{'missing' if differs else 'stray'}")
        oerr = outs[0][0] == "E"
        if isinstance(port, Err) and oerr: continue
        b3_n[fn][0] += 1
        if noted:
            b3_n[fn][1] += 1
            mm_ = B3_RE.search(note)
            if isinstance(port, Err):
                ok_ = mm_ is not None and mm_.group(3) == "Err(InvalidInput)" and mm_.group(4) == "E" and mm_.group(5) == \
                    "InvalidInput"
            else:
                ok_ = mm_ is not None and mm_.group(4) == "f" and int(mm_.group(5), 16) == f2b(port) and same_bits(
                    float(mm_.group(3)), port)
            check(ok_ and int(mm_.group(1)) == adr and mm_.group(2) == how, f"{where}: the #{adr} note does not name "
                  f"{how} = {port!r}")
        if isinstance(port, Err) or oerr:
            b3_n[fn][4] += 1
            b3_cov["the port's Err, Rust's value (co-volume)" if isinstance(port, Err) else
                   "the port's value, Rust's Err (co-volume)"] += 1
            check(fn in ("gas::van_der_waals_pressure", "gas::peng_robinson_pressure", "gas::redlich_kwong_pressure"),
                  f"{where}: {fn}: an Err of only one side outside the co-volume tests")
            continue
        if fn in B3_UX:
            inr = b3_in_range(fn, a)
            b3_n[fn][2] += inr
            check(not (inr and noted), f"{where}: {fn}{tuple(a)}: Rust's chain stays in the normal range, so its bits are "
                  f"the port's, but the row is noted")
        o = outs[0][1]
        if (o == 0.0 or not math.isfinite(o)) and not same_bits(o, port):
            zero_in = any(isinstance(v, float) and v == 0.0 for v in a)
            if math.isnan(o) and math.isinf(port): b3_cov["Rust's NaN, the value ±∞"] += 1
            elif math.isfinite(port) and (port != 0.0 or zero_in):
                b3_n[fn][3] += 1
                b3_cov["Rust 0, ±∞ or NaN, the value representable"] += 1
            else:
                b3_cov["Rust's special differs, the value out of range"] += 1
        if fn == "potential::lennard_jones" and a[0] >= 2.0 ** 1022:
            b3_fold["#17 rows (ε ≥ 2¹⁰²²)"] += 1
            e_, s, r = a
            sr = dv(s, r)
            p6 = sr * sr * sr * sr * sr * sr
            d = p6 * p6 - p6
            OOR[0] = False
            p = t3m(sr, s)
            for _ in range(5): p = t3m(p * sr, p)
            t3a(t3m(p * p, p) - p)
            old17 = e_ * d * 4.0
            if not OOR[0]:
                b3_fold["#17 rows with d normal (proven: (ε·d)·4 the port's)"] += 1
                check(same_bits(old17, port), f"{where}: #17's (ε·d)·4 = {old17!r} is not the port's {port!r} although "
                      f"Rust's d is normal")
            if same_bits(old17, port):
                b3_fold["#17's (ε·d)·4 the port's"] += 1
                b3_fold["#17's (ε·d)·4 the port's, Rust's differs (#17's notes)"] += not same_bits(old17, o)
            else:
                b3_fold["#17's (ε·d)·4 not the port's (d out of range)"] += 1
        if fn == "potential::morse_width_parameter":
            w, mu, d = a
            wm = MPFR.ux_mul_sqrt_div_final(U(w), U(mu), U(2.0) * U(d))
            check(same_bits(wm, port), f"{where}: MPFR's ω·√(μ/(2·D_e)) {wm!r} vs the integer root's {port!r}")
            b3_fold["√ checked against MPFR"] += 1
            if d >= 2.0 ** 1023:
                b3_fold["#18 rows (D_e ≥ 2¹⁰²³)"] += 1
                old18 = w * r_sqrt(rn_frac(Fraction(mu) / (2 * Fraction(d))))
                if same_bits(old18, o):
                    b3_fold["#18's value Rust's 0, the new one representable"] += not same_bits(port, o)
                else:
                    k = abs(ulp_dist(f2b(old18), f2b(port))) if math.isfinite(old18) and math.isfinite(port) else 99
                    b3_fold["#18 notes: the value kept" if k == 0 else "#18 notes: one ulp off" if k == 1 else
                            "#18 notes: further off"] += 1
            q = Fraction(w) ** 2 * Fraction(mu) / (2 * Fraction(d))
            if math.isfinite(port) and Fraction(port) ** 2 == q: b3_cov["√ exact"] += 1
            if q / Fraction(w) ** 2 > Fraction(DMAX_V): b3_cov["√: μ/(2·D_e) past DBL_MAX"] += 1
            if q / Fraction(w) ** 2 < Fraction(2) ** -1074: b3_cov["√: μ/(2·D_e) below 2⁻¹⁰⁷⁴"] += 1
            if min(w, mu, d) < 2.0 ** -1022: b3_cov["√: a subnormal operand"] += 1
            if port < 2.0 ** -1022: b3_cov["√: subnormal result"] += 1
            if 0.0 < port < math.inf:
                h = Fraction(math.ulp(port)) / 2 if port >= 2.0 ** -1022 else Fraction(2) ** -1075
                for dd in (h, -h):
                    mid = Fraction(port) + dd
                    if mid * mid == q: b3_cov["√ a midpoint (tie)"] += 1
                    elif port >= 2.0 ** -1022 and abs(q / (mid * mid) - 1) < Fraction(1, 1 << 97):
                        b3_cov["√ within 2⁻⁹⁸ of a 53-bit midpoint"] += 1
            continue
        if fn == "gas::redlich_kwong_pressure":
            tc = a[0]
            m_, ex = math.frexp(tc)                       # tc = m·2^ex, m in [½, 1)
            ex -= 1
            s_, j = (m_ * 2.0, ex // 2) if ex % 2 == 0 else (m_ * 4.0, (ex - 1) // 2)
            with MPFR.unbounded():
                lhs = U(tc).pow(2.5)
                rhs = MPFR.ux_scaled(MPFR.pow(s_, 2.5), 5 * j)
            check(1.0 <= s_ < 4.0 and lhs.cmp(rhs) == 0, f"{where}: RN₅₃(Tc^2.5) is not RN(s^2.5)·2^(5j) for Tc = {tc!r}")
            b3_fold["Tc^2.5 rescale checked"] += 1
            b3_fold["RK rows whose glibc pow(Tc, 2.5) misrounds"] += not same_bits(p_pow25_glibc(tc), MPFR.pow(tc, 2.5))
            if not (2.0 ** -1021 <= MPFR.pow(tc, 2.5) < math.inf): b3_cov["Tc^2.5 past either end"] += 1
        if fn in ("gas::ideal_gas_pressure", "gas::ideal_gas_volume", "gas::gas_density", "potential::coulomb",
                  "potential::coulomb_force"):
            if fn.startswith("potential::") and not all(is_finite(v) for v in a): continue
            if fn in ("gas::ideal_gas_pressure", "gas::ideal_gas_volume"): num, den = [a[0], R_GAS, a[1]], [a[2]]
            elif fn == "gas::gas_density": num, den = [a[0], a[1]], [R_GAS, a[2]]
            else: num, den = [KE_POT, a[0], a[1]], [a[2]] if fn == "potential::coulomb" else [a[2], a[2]]
            if any(v == 0.0 for v in num): continue
            q = Fraction(1)
            for v in num: q *= abs(Fraction(v))
            for v in den: q /= abs(Fraction(v))
            if q < Fraction(2) ** -1022:
                b3_cov["subnormal result"] += 1
                b3_cov["subnormal tie"] += b1_midpoint(q)
            elif b1_midpoint(q):
                b3_cov["53-bit midpoint"] += 1
for fn, (ok_, nt, inr, ql, fl_) in b3_n.items():
    extra = f", {inr} with Rust's chain in the normal range (all unnoted, Rust's bits)" if fn in B3_UX else ""
    extra += f", {fl_} where the co-volume test decides otherwise than Rust's" if fl_ else ""
    print(f"{fn.split('::')[0]}: ADR 0002 #{B3_HOW[fn][0]}: {fn.split('::')[1]}: {nt} of {ok_} rows noted with the port's "
          f"result ({ql} where Rust's is 0, ±∞ or NaN and the value representable){extra}")
    check(nt > 0, f"{fn}: no noted row")
print(f"ADR 0002 #31, #32 hard rows: {b3_cov}")
print(f"ADR 0002 #32's folds and checks: {b3_fold}")
check(all(v > 0 for k, v in b3_cov.items() if "out of range" not in k), f"B3: a hard-row class has no golden row: {b3_cov}")
check(b3_fold["#17's (ε·d)·4 the port's, Rust's differs (#17's notes)"] >= 180, f"B3: fewer #17 rows than M4d's 180: "
      f"{b3_fold}")
check(sum(b3_fold[k] for k in b3_fold if k.startswith("#18 notes")) == 181, f"B3: M4d's #18 rows are not 181: {b3_fold}")
check(b3_fold["RK rows whose glibc pow(Tc, 2.5) misrounds"] == 0, "B3: a redlich_kwong_pressure row whose glibc pow "
      "misrounds (its note, against glibc's oracle, would not be against the corrected one)")

# ── 22. The overflow audit's B4 (ADR 0002 #33, #34): the sums and the Shomate polynomials ──────────────────────────────
# Every row of the fourteen functions is recomputed from the port's semantics, independently of the generator and the port:
# Rust's operation sequence on mpfr_ref.UX (53 bits at MPFR's widest exponent range; Rust's sums by mpfr_ref.ux_sum in
# Rust's order, `.sum()` from −0.0, thermochem's `total += coeff * prop` from +0.0), the final operation rounded once into
# binary64 (ux_add_final, ux_sub_final, ux_mul_final, ux_sum_final); Rust's own expression where an input is not finite,
# and for the Shomate polynomials at T = ±0 (every product is then an exact zero and E/0 IEEE's special in both
# evaluations); zero_order_concentration's clamp decided on the unbounded difference. Rust's validation is replicated
# (rows where it errs carry no note). Every row: noted exactly where the port's value differs from the oracle's (by its
# bits, NaN one value), the note naming it bit for bit. No row whose f64 chain stays in the normal range is noted (B3's
# strict tracers: a product or quotient leaves the range when it is not finite, 0 from nonzero factors, or nonzero and at
# most 2⁻¹⁰²²; a sum only when it is not finite; the final operation's own range does not count).
# The pure-sum lemma (src/thermo.cyr's comment on #33's sums: a sum of two doubles is exact below 2⁻¹⁰²¹ and the 53-bit
# rounding above, unless it overflows) is checked on every row of the three pure sums (hess_law, enthalpy_from_formation,
# q_value's two sums) with finite inputs: where Rust's result is finite, every partial sum is exactly MPFR's.
# The data-backed Shomate callers: every SHOMATE_DATA coefficient is 0 or in [2⁻²⁰⁰, 2²⁰⁰] and every range end in
# [2⁻³⁰, 2³⁰], the port's box (src/thermochem.cyr's _km_shomate_box), so enthalpy_change_cp, which evaluates the
# polynomials only at temperatures within the ranges, always takes Rust's chain; adiabatic_flame_temperature's evaluations
# (Newton's iterates may leave the ranges) are replayed on every golden flame row, each inside the box, traced normal, or
# equal to MPFR's.
B4_ADD = "Rust's sequence evaluated with an unbounded exponent (the final addition rounded once)"
B4_HOW = {
    "thermo::hess_law": (33, B4_ADD),
    "thermo::enthalpy_from_formation": (33, B3_SUB),
    "thermochem::reaction_enthalpy": (33, B3_MUL),
    "thermochem::reaction_gibbs_energy": (33, B3_MUL),
    "thermochem::reaction_entropy": (33, B3_SUB),
    "nuclear::q_value": (33, B3_MUL),
    "solution::ionic_strength": (33, B3_MUL),
    "inorganic::born_haber_lattice_energy": (33, B3_SUB),
    "reaction::gibbs_free_energy": (33, B3_SUB),
    "solution::kohlrausch_molar_conductivity": (33, B3_SUB),
    "solution::limiting_molar_conductivity": (33, B4_ADD),
    "kinetics::zero_order_concentration": (33, B3_SUB),
    "Shomate::cp": (34, B4_ADD),
    "Shomate::enthalpy_antiderivative": (34, B3_MUL),
}
B4_COL = {"thermochem::reaction_enthalpy": 0, "thermochem::reaction_gibbs_energy": 1, "thermochem::reaction_entropy": 2}
def b4_ions(ins):
    n = ins[0][1]
    return [(ins[1 + 2 * j][1], ins[2 + 2 * j][1]) for j in range(n)]
def b4_species(ins):
    prod, k = species_at(ins, 0)
    reac, _ = species_at(ins, k)
    return prod, reac
def b4_port(fn, ins):
    """The port's result of the row (a double, or Err)."""
    if fn == "thermo::hess_law":
        xs = slice_args(ins)[0]
        return MPFR.ux_sum_final(xs) if all(map(math.isfinite, xs)) else r_sum(xs)
    if fn == "thermo::enthalpy_from_formation":
        p_, r_ = slice_args(ins)
        if not all(map(math.isfinite, p_ + r_)): return r_sum(p_) - r_sum(r_)
        with MPFR.unbounded(): return MPFR.ux_sub_final(MPFR.ux_sum(p_), MPFR.ux_sum(r_))
    if fn == "nuclear::q_value":
        r_, p_ = slice_args(ins)
        if not all(map(math.isfinite, r_ + p_)): return r_q_value(ins)
        with MPFR.unbounded(): return MPFR.ux_mul_final(MPFR.ux_sum(r_) - MPFR.ux_sum(p_), U(AMU_MEV))
    if fn == "solution::ionic_strength":
        ions = b4_ions(ins)
        if not all(math.isfinite(c) for c, _ in ions): return r_ionic(ins)
        with MPFR.unbounded():
            return MPFR.ux_mul_final(U(0.5), MPFR.ux_sum([U(c) * U(float(z)) * U(float(z)) for c, z in ions]))
    if fn in B4_COL:
        prod, reac = b4_species(ins)
        col = B4_COL[fn]
        for f_, _ in prod + reac:
            if tc_lookup(f_) is None: return Err("InvalidReaction")
        if not all(math.isfinite(c) for _, c in prod + reac): return r_reaction(prod, reac, col, col != 2)
        with MPFR.unbounded():
            sp_ = MPFR.ux_sum([U(c) * U(tc_lookup(f_)[2][col]) for f_, c in prod], 0.0)
            sr = MPFR.ux_sum([U(c) * U(tc_lookup(f_)[2][col]) for f_, c in reac], 0.0)
            return MPFR.ux_sub_final(sp_, sr) if col == 2 else MPFR.ux_mul_final(sp_ - sr, U(1000.0))
    a = fl(ins)
    if fn == "inorganic::born_haber_lattice_energy":
        if not all(map(math.isfinite, a)): return r_haber(*a)
        with MPFR.unbounded(): return MPFR.ux_sub_final(U(a[0]) - U(a[1]) - U(a[2]) - U(a[3]) / U(2.0), U(a[4]))
    if fn == "reaction::gibbs_free_energy":
        h, t, s_ = a
        if not is_positive(t): return INVT
        if not (is_finite(h) and is_finite(s_)): return INVALID
        with MPFR.unbounded(): return MPFR.ux_sub_final(U(h), U(t) * U(s_))
    if fn == "solution::kohlrausch_molar_conductivity":
        l0, k, c = a
        if not is_non_negative(c): return Err("InvalidConcentration")
        if not (is_finite(l0) and is_finite(k)): return INVALID
        with MPFR.unbounded(): return MPFR.ux_sub_final(U(l0), U(k) * U(c).sqrt())
    if fn == "solution::limiting_molar_conductivity":
        cat, nc, an, na = a
        if not all(map(math.isfinite, a)): return nc * cat + na * an
        with MPFR.unbounded(): return MPFR.ux_add_final(U(nc) * U(cat), U(na) * U(an))
    if fn == "kinetics::zero_order_concentration":
        c0, k, t = a
        if not all(map(math.isfinite, a)): return r_zero(c0, k, t)
        with MPFR.unbounded():
            p_ = U(k) * U(t)
            if (U(c0) - p_).cmp(U(0.0)) < 0: return 0.0
            return MPFR.ux_sub_final(U(c0), p_)
    c, temp = a[:5], a[7]
    if not all(map(math.isfinite, c + [temp])) or temp == 0.0:
        return _sh_cp0(c, temp) if fn == "Shomate::cp" else _sh_h0(c, temp)
    return b4_shomate_ux(fn == "Shomate::cp", c, temp)
def b4_shomate_ux(is_cp, c, temp):
    aa, bb_, cc, dd, ee = (U(x) for x in c)
    with MPFR.unbounded():
        t = U(temp) / U(1000.0)
        if is_cp: return MPFR.ux_add_final(aa + bb_ * t + cc * t * t + dd * t * t * t, ee / (t * t))
        inner = t * (aa + t * (bb_ / U(2.0) + t * (cc / U(3.0) + t * (dd / U(4.0)))))
        return MPFR.ux_mul_final(U(1000.0), inner - ee / t)
def b4_shomate_in_range(is_cp, c, temp):
    """Rust's f64 Shomate chain traced (B3's strict range): True when it stays in the normal range."""
    a_, b_, c_, d_, e_ = c
    OOR[0] = False
    if temp == 0.0: return True
    t = t3m(dv(temp, 1000.0), temp)
    if is_cp:
        s1 = t3a(a_ + t3m(b_ * t, b_, t))
        s2 = t3a(s1 + t3m(t3m(c_ * t, c_, t) * t, c_, t))
        t3a(s2 + t3m(t3m(t3m(d_ * t, d_, t) * t, d_, t) * t, d_, t))
        t3m(dv(e_, t3m(t * t, t)), e_)
        return not OOR[0]
    h4 = t3a(t3m(dv(c_, 3.0), c_) + t3m(t * t3m(dv(d_, 4.0), d_), t, d_))
    h7 = t3a(t3m(dv(b_, 2.0), b_) + t3m(t * h4, t, h4))
    inner = t3m(t * t3a(a_ + t3m(t * h7, t, h7)), t)
    t3a(inner - t3m(dv(e_, t), e_))
    return not OOR[0]
def b4_in_range(fn, ins):
    """True when Rust's f64 chain of the row (inputs finite, validated) stays in the normal range."""
    OOR[0] = False
    if fn == "thermo::hess_law":
        xs = slice_args(ins)[0]
        acc = -0.0
        for x in xs[:-1]: acc = t3a(acc + x)
        return not OOR[0]
    if fn in ("thermo::enthalpy_from_formation", "nuclear::q_value"):
        sums_ = []
        for xs in slice_args(ins):
            acc = -0.0
            for x in xs: acc = t3a(acc + x)
            sums_.append(acc)
        if fn == "nuclear::q_value": t3a(sums_[0] - sums_[1])
        return not OOR[0]
    if fn == "solution::ionic_strength":
        acc = -0.0
        for c, z in b4_ions(ins):
            zf = float(z)
            acc = t3a(acc + t3m(t3m(c * zf, c, zf) * zf, c, zf))
        return not OOR[0]
    if fn in B4_COL:
        prod, reac = b4_species(ins)
        col = B4_COL[fn]
        sums_ = []
        for side in (prod, reac):
            acc = 0.0
            for f_, c in side:
                v = tc_lookup(f_)[2][col]
                acc = t3a(acc + t3m(c * v, c, v))
            sums_.append(acc)
        if col != 2: t3a(sums_[0] - sums_[1])
        return not OOR[0]
    a = fl(ins)
    if fn == "inorganic::born_haber_lattice_energy":
        t3a(t3a(t3a(a[0] - a[1]) - a[2]) - t3m(dv(a[3], 2.0), a[3]))
        return not OOR[0]
    if fn == "reaction::gibbs_free_energy":
        t3m(a[1] * a[2], a[1], a[2])
        return not OOR[0]
    if fn == "solution::kohlrausch_molar_conductivity":
        sq = math.sqrt(a[2])
        t3m(a[1] * sq, a[1], sq)
        return not OOR[0]
    if fn == "solution::limiting_molar_conductivity":
        t3m(a[1] * a[0], a[1], a[0])
        t3m(a[3] * a[2], a[3], a[2])
        return not OOR[0]
    if fn == "kinetics::zero_order_concentration":
        t3m(a[1] * a[2], a[1], a[2])
        return not OOR[0]
    return b4_shomate_in_range(fn == "Shomate::cp", a[:5], a[7])
def b4_inputs_finite(fn, ins):
    if fn in ("thermo::hess_law", "thermo::enthalpy_from_formation", "nuclear::q_value"):
        return all(math.isfinite(x) for xs in slice_args(ins) for x in xs)
    if fn == "solution::ionic_strength": return all(math.isfinite(c) for c, _ in b4_ions(ins))
    if fn in B4_COL:
        prod, reac = b4_species(ins)
        return all(math.isfinite(c) for _, c in prod + reac)
    a = fl(ins)
    if fn.startswith("Shomate::"): return all(map(math.isfinite, a[:5] + [a[7]]))
    return all(map(math.isfinite, a))
def b4_lemma(fn, ins):
    """The pure-sum lemma on the row: (True, checked) where Rust's sums are finite and every partial sum is MPFR's."""
    groups = slice_args(ins)
    for xs in groups:
        acc, ux = -0.0, U(-0.0)
        with MPFR.unbounded():
            for x in xs:
                acc = acc + x
                ux = ux + U(x)
                if not math.isfinite(acc): return True, False
                if not ux.is_f64(acc): return False, True
    return True, True
B4_RE = re.compile(r"\[ADR 0002 #(33|34): the port returns (.+) = (\S+) \(f:([0-9a-f]{16})\)\]$")
B4_ADR1 = "[J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002 #1)]"
b4_n = {fn: [0, 0, 0] for fn in B4_HOW}             # rows with a value (either side), noted, Rust's chain in range
b4_cov = {"Rust 0, ±∞ or NaN, the value representable": 0, "Rust's NaN, the value ±∞": 0,
          "Rust's special differs, the value out of range": 0, "a signed zero differs": 0,
          "degraded: the value normal, Rust's more than 2 ulp off": 0, "within 2 ulp, or the value subnormal": 0}
b4_pure = [0, 0]                                      # rows the lemma checked, rows with a non-finite f64 sum
_sh_cp0, _sh_h0 = sh_cp, sh_h
for mod_ in ("thermo", "thermochem", "nuclear", "solution", "inorganic", "reaction", "kinetics"):
    for fn, cls, ins, outs, where in R.get(mod_, []):
        if fn not in B4_HOW: continue
        adr, how = B4_HOW[fn]
        port = b4_port(fn, ins)
        note = NOTES.get(where, "")
        noted = f"[ADR 0002 #{adr}:" in note
        pt = out_toks(port)
        differs = not toks_equal(pt, outs)
        check(noted == differs, f"{where}: {fn}: the port's {pt} vs the oracle's {outs}: the #{adr} note is "
              f"{'missing' if differs else 'stray'}")
        if fn in ("thermochem::reaction_enthalpy", "thermochem::reaction_gibbs_energy") and noted:
            check(B4_ADR1 + " [ADR 0002 #33: " in note, f"{where}: {fn}'s #33 note does not follow ADR 0002 #1's")
        if isinstance(port, Err): continue
        b4_n[fn][0] += 1
        if noted:
            b4_n[fn][1] += 1
            mm_ = B4_RE.search(note)
            check(mm_ is not None and int(mm_.group(1)) == adr and mm_.group(2) == how and int(mm_.group(4), 16) == f2b(port)
                  and same_bits(float(mm_.group(3)), port), f"{where}: the #{adr} note does not name {how} = {port!r}")
        fin = b4_inputs_finite(fn, ins)
        if fin:
            inr = b4_in_range(fn, ins)
            b4_n[fn][2] += inr
            check(not (inr and noted), f"{where}: {fn}: Rust's chain stays in the normal range, so its bits are the "
                  f"port's, but the row is noted")
        if fin and fn in ("thermo::hess_law", "thermo::enthalpy_from_formation", "nuclear::q_value"):
            ok_, checked = b4_lemma(fn, ins)
            check(ok_, f"{where}: {fn}: a finite partial sum of doubles is not the unbounded one (the lemma)")
            b4_pure[0 if checked else 1] += 1
        if not noted: continue
        o = outs[0][1]
        if o == 0.0 and port == 0.0:
            b4_cov["a signed zero differs"] += 1
        elif o == 0.0 or not math.isfinite(o):
            if math.isnan(o) and math.isinf(port): b4_cov["Rust's NaN, the value ±∞"] += 1
            elif math.isfinite(port): b4_cov["Rust 0, ±∞ or NaN, the value representable"] += 1
            else: b4_cov["Rust's special differs, the value out of range"] += 1
        elif math.isfinite(port) and abs(port) >= 2.0 ** -1022 and abs(ulp_dist(f2b(o), f2b(port))) > 2:
            b4_cov["degraded: the value normal, Rust's more than 2 ulp off"] += 1
        else:
            b4_cov["within 2 ulp, or the value subnormal"] += 1
# The data-backed Shomate callers: the static box over SHOMATE_DATA, then every flame row's evaluations replayed.
B4_LO, B4_HI = 2.0 ** -200, 2.0 ** 200
def b4_box(c, temp):
    return 2.0 ** -30 <= abs(temp) <= 2.0 ** 30 and all(x == 0.0 or B4_LO <= abs(x) <= B4_HI for x in c[:5])
check(all(all(x == 0.0 or B4_LO <= abs(x) <= B4_HI for x in c[:5]) and 2.0 ** -30 <= c[5] <= c[6] <= 2.0 ** 30
          for _, c in SH), "B4: a SHOMATE_DATA segment outside the port's box (enthalpy_change_cp would not be provably "
                           "Rust's chain)")
B4_EVALS = []
def sh_cp(c, temp):
    B4_EVALS.append((True, c, temp))
    return _sh_cp0(c, temp)
def sh_h(c, temp):
    B4_EVALS.append((False, c, temp))
    return _sh_h0(c, temp)
b4_flame = {"rows": 0, "evaluations": 0, "in the box": 0, "traced normal": 0, "checked against MPFR": 0}
for fn, cls, ins, outs, where in R.get("thermochem", []):
    if fn != "thermochem::adiabatic_flame_temperature": continue
    B4_EVALS.clear()
    RULES[fn][1](ins)
    fprods_, fk_ = species_at(ins, 1)
    r_flame(ins[0][1], fprods_, ins[fk_][1], exact=True)   # the port's own evaluations (ADR 0002 #38's rows go on)
    b4_flame["rows"] += 1
    for is_cp, c, temp in B4_EVALS:
        b4_flame["evaluations"] += 1
        if b4_box(c, temp): b4_flame["in the box"] += 1
        elif not all(map(math.isfinite, list(c[:5]) + [temp])) or b4_shomate_in_range(is_cp, list(c[:5]), temp):
            b4_flame["traced normal"] += 1
        else:
            b4_flame["checked against MPFR"] += 1
            v = _sh_cp0(c, temp) if is_cp else _sh_h0(c, temp)
            check(same_bits(v, b4_shomate_ux(is_cp, list(c[:5]), temp)), f"{where}: a flame evaluation of Shomate "
                  f"{'cp' if is_cp else 'H'} at {temp!r} would differ under ADR 0002 #34")
sh_cp, sh_h = _sh_cp0, _sh_h0
for fn, (ok_, nt, inr) in b4_n.items():
    print(f"{fn.split('::')[0]}: ADR 0002 #{B4_HOW[fn][0]}: {fn.split('::')[1]}: {nt} of {ok_} rows with a value noted "
          f"with the port's, {inr} with Rust's chain in the normal range (all unnoted, Rust's bits)")
    check(nt > 0, f"{fn}: no noted row")
print(f"ADR 0002 #33, #34 noted rows: {b4_cov}")
print(f"ADR 0002 #33's lemma: {b4_pure[0]} pure-sum rows with finite inputs and finite sums checked partial sum by "
      f"partial sum against MPFR, {b4_pure[1]} with a sum past DBL_MAX")
print(f"ADR 0002 #34's data-backed callers: every SHOMATE_DATA segment in the port's box; the flame rows' evaluations: "
      f"{b4_flame}")
check(b4_cov["Rust 0, ±∞ or NaN, the value representable"] > 0 and b4_cov["Rust's NaN, the value ±∞"] > 0,
      f"B4: a hard-row class has no golden row: {b4_cov}")

# ── 23. The overflow audit's B5 (ADR 0002 #35): the exp / exp10 arguments ──────────────────────────────────────────────
# Every row of the four functions is recomputed from the port's semantics, independently of the generator and the port:
# Rust's validation in Rust's order, then the argument — equilibrium_constant's −ΔG/(R·T), the activity models' log₁₀ γ
# (((−A)·z)·z)·√I, ((((−A)·z)·z)·√I)/(1 + (B·å)·√I) and (((−A)·z)·z)·(√I/(1 + √I) − 0.3·I) — in Rust's sequence on
# mpfr_ref.UX (53 bits at MPFR's widest exponent range), every operation rounded to 53 bits, the last one too: X, the
# argument of ADR 0002 #20's target; the result is MPFR's exp / exp10 of X rounded once into binary64
# (mpfr_ref.ux_exp_final, ux_exp10_final), and equilibrium_constant's ComputationError where that K is not finite.
# ADR 0002 #35's lemma — the port's exp / 10^x of x, X's last operation rounded once into binary64 (ux_div_final,
# ux_mul_final), is that value — is checked on every row (x normal: the same double; |X| past DBL_MAX: x = ±∞, both
# ∞ or 0; |X| ≤ 2⁻¹⁰²¹: |x| ≤ 2⁻¹⁰²¹, both 1). Every row: noted exactly where the port's result differs from the row's
# reference (the corrected oracle's outputs where tests/golden/<module>.cr.tsv has the row, else the oracle's), the note
# naming it bit for bit (so every noted value, glibc's exp or pow(10, x) in the generator, is MPFR's), or its Err
# (`= Err(ComputationError) (E:ComputationError)`). The port's fast tests (src/reaction.cyr, src/solution.cyr: R·T in
# [2⁻¹⁰²¹, DBL_MAX]; A·z² finite; the extended law's numerator and (B·å)·√I finite) are checked on every row: where one
# passes the row is unnoted (Rust's result). And no row whose f64 chain stays in the normal range is noted (B3's strict
# tracers: a product or quotient leaves the range when it is not finite, 0 from nonzero factors, or nonzero and at most
# 2⁻¹⁰²²; a sum only when it is not finite; the argument's own last operation does not count). The noted rows are
# counted by class.
B5_DIV = "the argument evaluated with an unbounded exponent (its final quotient rounded once), then the correctly rounded"
B5_MUL = "the argument evaluated with an unbounded exponent (its final product rounded once), then the correctly rounded"
B5_HOW = {
    "reaction::equilibrium_constant": B5_DIV + " exp",
    "solution::debye_huckel_limiting": B5_MUL + " 10^x",
    "solution::debye_huckel_extended": B5_DIV + " 10^x",
    "solution::davies_activity": B5_MUL + " 10^x",
}
B5_KEQ, B5_LIM, B5_EXT, B5_DAV = B5_HOW
B5_LO, B5_HI = 2.0 ** -1021, DMAX_V
def b5_valid(fn, a):
    """Rust's validation (None: Ok), in Rust's order."""
    if fn == B5_KEQ:
        if not is_positive(a[1]): return INVT
        return None if is_finite(a[0]) else INVALID
    if fn == B5_EXT:
        if r_ia(a[1], a[2]): return INVALID
        return None if is_non_negative(a[3]) and is_non_negative(a[4]) else INVALID
    return INVALID if r_ia(a[1], a[2]) else None
def b5_arg(fn, a):
    """(X, x): the argument with every operation rounded to 53 bits at MPFR's widest exponent range (X, a UX), and with
    its last operation rounded once into binary64 instead (x, the port's)."""
    with MPFR.unbounded():
        if fn == B5_KEQ:
            n, rt = U(-a[0]), U(R_GAS) * U(a[1])
            return n / rt, MPFR.ux_div_final(n, rt)
        z = U(float(a[0]))
        p_ = U(-a[2]) * z * z
        sq = U(a[1]).sqrt()
        if fn == B5_LIM: return p_ * sq, MPFR.ux_mul_final(p_, sq)
        if fn == B5_EXT:
            num, den = p_ * sq, U(1.0) + U(a[3]) * U(a[4]) * sq
            return num / den, MPFR.ux_div_final(num, den)
        br = sq / (U(1.0) + sq) - U(0.3) * U(a[1])
        return p_ * br, MPFR.ux_mul_final(p_, br)
def b5_lemma_class(xx):
    """The lemma's case for the argument X: its 53-bit value normal, past DBL_MAX, or at most 2⁻¹⁰²¹ (0 included)."""
    if xx.is_zero(): return "|X| ≤ 2⁻¹⁰²¹"
    if U(B5_HI).cmp(xx) < 0 or U(-B5_HI).cmp(xx) > 0: return "|X| past DBL_MAX"
    if U(B5_LO).cmp(xx) > 0 and U(-B5_LO).cmp(xx) < 0: return "|X| ≤ 2⁻¹⁰²¹"
    return "X normal"
def b5_port(fn, a):
    """The port's result of the row (a double, or Err) and the lemma's case."""
    e = b5_valid(fn, a)
    if e is not None: return e, None
    xx, x = b5_arg(fn, a)
    if fn == B5_KEQ:
        y, y64 = MPFR.ux_exp_final(xx), MPFR.exp(x)
    else:
        y, y64 = MPFR.ux_exp10_final(xx), MPFR.exp10(x)
    cls = b5_lemma_class(xx)
    check(same_bits(y, y64), f"{fn}{tuple(a)}: ADR 0002 #35's lemma: the correctly rounded exp of the 53-bit argument "
          f"{y!r} is not that of its binary64 rounding {x!r}, {y64!r} ({cls})")
    if fn == B5_KEQ and not math.isfinite(y): return Err("ComputationError"), cls
    return y, cls
def b5_fast(fn, a):
    """True where the port's fast test passes on the row (validated): Rust's chain is its result."""
    if fn == B5_KEQ: return B5_LO <= R_GAS * a[1] <= B5_HI
    z = float(a[0])
    p_ = -a[2] * z * z
    if fn != B5_EXT: return math.isfinite(p_)
    sq = math.sqrt(a[1])
    return math.isfinite(p_ * sq) and math.isfinite(a[3] * a[4] * sq)
def b5_in_range(fn, a):
    """True when Rust's f64 chain of the row (validated) stays in the normal range, its last operation apart."""
    OOR[0] = False
    if fn == B5_KEQ:
        t3m(R_GAS * a[1], R_GAS, a[1])
        return not OOR[0]
    z, ii, aa = float(a[0]), a[1], a[2]
    p1 = t3m(-aa * z, aa, z)
    p_ = t3m(p1 * z, p1, z)
    sq = math.sqrt(ii)
    if fn == B5_EXT:
        t3m(p_ * sq, p_, sq)
        b1 = t3m(a[3] * a[4], a[3], a[4])
        t3a(1.0 + t3m(b1 * sq, b1, sq))
    elif fn == B5_DAV:
        d1 = t3a(1.0 + sq)
        q = t3m(dv(sq, d1), sq)
        t3a(q - t3m(0.3 * ii, 0.3, ii))
    return not OOR[0]
B5_RE = re.compile(r"\[ADR 0002 #35: the port returns (.+) = (\S+) \((f:([0-9a-f]{16})|E:ComputationError)\)\]$")
b5_n = {fn: [0, 0, 0, 0] for fn in B5_HOW}            # rows past validation, noted, fast test passing, chain in range
b5_lemma = {"X normal": 0, "|X| past DBL_MAX": 0, "|X| ≤ 2⁻¹⁰²¹": 0}
b5_cov = {"Rust's K finite, the port's overflows (its ComputationError)": 0,
          "Rust's K overflows (ComputationError), the port's finite": 0,
          "Rust NaN, the value representable": 0, "Rust 0 or ∞, the value representable": 0,
          "Rust's special differs, the value out of range": 0,
          "degraded: the value normal, Rust's more than 2 ulp off": 0, "within 2 ulp, or the value subnormal": 0}
b5_ref = [0, 0]                                       # rows whose reference is a correction; noted rows among them
B5_CORR = {w: ct for m_ in ("reaction", "solution") for w, ct, _ in CORRECTIONS.get(m_, [])}
for mod_ in ("reaction", "solution"):
    for fn, cls, ins, outs, where in R.get(mod_, []):
        if fn not in B5_HOW: continue
        how = B5_HOW[fn]
        a = fl(ins)
        port, lcls = b5_port(fn, a)
        note = NOTES.get(where, "")
        noted = "[ADR 0002 #35:" in note
        ref = B5_CORR.get(where, outs)
        if where in B5_CORR:
            b5_ref[0] += 1
            b5_ref[1] += noted
        pt = out_toks(port)
        differs = not toks_equal(pt, ref)
        check(noted == differs, f"{where}: {fn}{tuple(a)}: the port's {pt} vs the reference {ref}: the #35 note is "
              f"{'missing' if differs else 'stray'}")
        if lcls is None:
            check(not noted, f"{where}: a row Rust's validation refuses is noted")
            continue
        b5_n[fn][0] += 1
        b5_lemma[lcls] += 1
        fast = b5_fast(fn, a)
        b5_n[fn][2] += fast
        check(not (fast and noted), f"{where}: {fn}{tuple(a)}: the port's fast test passes, so its result is Rust's, but "
              f"the row is noted")
        inr = b5_in_range(fn, a)
        b5_n[fn][3] += inr
        check(not (inr and noted), f"{where}: {fn}{tuple(a)}: Rust's chain stays in the normal range, so its result is "
              f"the port's, but the row is noted")
        if not noted: continue
        b5_n[fn][1] += 1
        mm_ = B5_RE.search(note)
        if isinstance(port, Err):
            ok_ = mm_ is not None and mm_.group(1) == how and mm_.group(2) == "Err(ComputationError)" and mm_.group(4) is None
        else:
            ok_ = (mm_ is not None and mm_.group(1) == how and mm_.group(4) is not None
                   and int(mm_.group(4), 16) == f2b(port) and same_bits(float(mm_.group(2)), port))
        check(ok_, f"{where}: the #35 note does not name {how} = {port!r}")
        o = ref[0]
        if isinstance(port, Err):
            b5_cov["Rust's K finite, the port's overflows (its ComputationError)"] += 1
        elif o[0] == "E":
            b5_cov["Rust's K overflows (ComputationError), the port's finite"] += 1
        elif math.isnan(o[1]) or o[1] == 0.0 or math.isinf(o[1]):
            if math.isfinite(port) and port != 0.0:
                b5_cov["Rust NaN, the value representable" if math.isnan(o[1]) else "Rust 0 or ∞, the value representable"] += 1
            else:
                b5_cov["Rust's special differs, the value out of range"] += 1
        elif port >= 2.0 ** -1022 and abs(ulp_dist(o[2], f2b(port))) > 2:
            b5_cov["degraded: the value normal, Rust's more than 2 ulp off"] += 1
        else:
            b5_cov["within 2 ulp, or the value subnormal"] += 1
for fn, (ok_, nt, fast, inr) in b5_n.items():
    print(f"{fn.split('::')[0]}: ADR 0002 #35: {fn.split('::')[1]}: {nt} of {ok_} rows past Rust's validation noted with "
          f"the port's result, {fast} with the port's fast test passing and {inr} with Rust's chain in the normal range "
          f"(all unnoted, Rust's result)")
    check(nt > 0, f"{fn}: no noted row")
print(f"ADR 0002 #35 noted rows: {b5_cov}")
print(f"ADR 0002 #35's lemma, checked on every row past validation: {b5_lemma}; {b5_ref[0]} rows compared against the "
      f"corrected oracle ({b5_ref[1]} of them noted)")
check(all(v > 0 for k, v in b5_cov.items() if "out of range" not in k), f"B5: a noted-row class has no golden row: "
      f"{b5_cov}")
check(all(v > 0 for v in b5_lemma.values()), f"B5: a case of the lemma has no golden row: {b5_lemma}")

# ── 24. The overflow audit's B6 (ADR 0002 #36, #37, #38): the FFT, find_peaks' threshold, the flame's bracket ─────────
# #38, adiabatic_flame_temperature: every row is recomputed with the port's solver (r_flame with `exact`: the bracket's
# and the bisection's sign tests compare the signs exactly, everything else Rust's) and must be noted exactly where its
# result differs from the oracle's, the note naming the value bit for bit and its solver path. The sign tests are
# checked against Rust's product tests over every pair of a set of special and boundary values (±0, the subnormal and
# normal ends, ±∞, NaN) and seeded pairs: they differ only where the product underflows to a zero from nonzero operands.
# On the rows: a noted row's ends have opposite signs and an underflowing product (Rust's ComputationError, the port's
# Ok); every row whose product is not such a zero is unnoted (strict); the bisection's test never refuses after the
# bracket passed.
B6_FLAME = "thermochem::adiabatic_flame_temperature"
B6_FLAME_HOW = "Rust's solver with an exact bracket sign test"
B6_FLAME_RE = re.compile(r"\[ADR 0002 #38: the port returns (.+) \((newton, \d+ steps|bisection after newton: [^)]*)\) "
                         r"= (\S+) \(f:([0-9a-f]{16})\)\]$")
def b6_lt_rust(a, b): return a * b < 0.0
def b6_lt_exact(a, b): return (a < 0.0 < b) or (b < 0.0 < a)
def b6_gt_rust(a, b): return a * b > 0.0
def b6_gt_exact(a, b): return (a > 0.0 and b > 0.0) or (a < 0.0 and b < 0.0)
B6_SPECIALS = [0.0, 2.0 ** -1074, 3 * 2.0 ** -1074, 2.0 ** -1060, 2.0 ** -1022, 2.0 ** -1000, 2.0 ** -537, 2.0 ** -538,
               1e-162, 1.5e-162, 1e-200, 1e-3, 1.0, 1e200, 2.0 ** 512, 2.0 ** 1000, DMAX_V, math.inf, math.nan]
B6_SPECIALS = B6_SPECIALS + [-x for x in B6_SPECIALS]
b6_pairs = [(a, b) for a in B6_SPECIALS for b in B6_SPECIALS]
_rng6 = random.Random(0x3806)
for _ in range(20000):
    b6_pairs.append(tuple(_rng6.choice([-1.0, 1.0]) * 2.0 ** _rng6.uniform(-1074, 1023) for _ in range(2)))
b6_diff = [0, 0]
for a_, b_ in b6_pairs:
    for k_, (rt, et) in enumerate(((b6_lt_rust, b6_lt_exact), (b6_gt_rust, b6_gt_exact))):
        if rt(a_, b_) != et(a_, b_):
            b6_diff[k_] += 1
            check(a_ * b_ == 0.0 and a_ != 0.0 and b_ != 0.0 and not math.isnan(a_) and not math.isnan(b_),
                  f"ADR 0002 #38: the sign test and Rust's product test differ at ({a_!r}, {b_!r}) beyond an underflow")
print(f"ADR 0002 #38's sign tests against Rust's products: {len(b6_pairs)} pairs; they differ on {b6_diff[0]} "
      f"(bracket) and {b6_diff[1]} (bisection) pairs, each an underflowing product")
check(b6_diff[0] > 0 and b6_diff[1] > 0, "ADR 0002 #38: no pair shows the underflow")
b6_fl = {"rows past validation": 0, "noted": 0, "the product an underflowed zero, signs equal (Err in both)": 0,
         "the product subnormal, nonzero (bracketed in both)": 0}
for fn, cls, ins, outs, where in R.get("thermochem", []):
    if fn != B6_FLAME: continue
    prods6, k6 = species_at(ins, 1)
    rust = r_flame(ins[0][1], prods6, ins[k6][1])
    ends = FLAME_ENDS[0]
    port = r_flame(ins[0][1], prods6, ins[k6][1], exact=True)
    ppath = FLAME_PATH[0]
    note = NOTES.get(where, "")
    noted = "[ADR 0002 #38:" in note
    differs = not toks_equal(out_toks(port), outs)
    check(noted == differs, f"{where}: the port's flame {port!r} vs the oracle's {outs}: the #38 note is "
          f"{'missing' if differs else 'stray'}")
    check(toks_equal(out_toks(rust), outs), f"{where}: Rust's flame recomputed differs from the oracle")
    if ends is None:
        check(not noted, f"{where}: a row Rust's validation refuses is noted")
        continue
    b6_fl["rows past validation"] += 1
    fa, fb = ends
    uz = fa * fb == 0.0 and fa != 0.0 and fb != 0.0
    check(uz or not noted, f"{where}: a noted flame row whose bracket product does not underflow to a zero")
    if uz and b6_gt_exact(fa, fb): b6_fl["the product an underflowed zero, signs equal (Err in both)"] += 1
    if 0.0 < abs(fa * fb) < 2.0 ** -1022: b6_fl["the product subnormal, nonzero (bracketed in both)"] += 1
    if not noted: continue
    b6_fl["noted"] += 1
    check(isinstance(rust, Err) and rust.name == "ComputationError" and not isinstance(port, Err) and b6_lt_exact(fa, fb),
          f"{where}: a #38 row that is not Rust's ComputationError against the port's Ok")
    mm_ = B6_FLAME_RE.search(note)
    check(mm_ is not None and mm_.group(1) == B6_FLAME_HOW and mm_.group(2) == ppath and int(mm_.group(4), 16) == f2b(port)
          and same_bits(float(mm_.group(3)), port), f"{where}: the #38 note does not name {B6_FLAME_HOW} ({ppath}) = "
          f"{port!r}")
print(f"thermochem: ADR 0002 #38: adiabatic_flame_temperature: {b6_fl}")
check(b6_fl["noted"] > 0 and b6_fl["the product an underflowed zero, signs equal (Err in both)"] > 0,
      f"ADR 0002 #38: a class has no golden row: {b6_fl}")

# #37, find_peaks: every row is recomputed with the port's semantics — Rust's checks and empty results, Rust's fold for
# min and max over the finite values, the neighbour tests as IEEE, the threshold test w > min·(1 − t) + max·t in
# Fraction arithmetic with t the given double (1 − t exact), a non-finite candidate as IEEE against a finite threshold
# (+∞ above) — and must be noted exactly where its indices differ from the oracle's, the note's tokens being them. On
# every row past the checks, each candidate is checked against the port's fast test (src/spectroscopy.cyr's
# _km_find_peaks_run): Rust's rounded threshold T within the proven bound 2.5u(1 + u)(|min| + |max|) + (1 + u)·2⁻¹⁰⁷⁴ of the
# exact θ, the port's computed bound B at least that, and where RN(w − T) is beyond ±B Rust's decision is the exact one;
# a row whose every candidate the fast test decides is unnoted (strict). The noted rows and the candidates are counted.
B6_PEAKS = "spectroscopy::find_peaks"
B6_PEAKS_HOW = "the threshold test evaluated exactly (w > min·(1 − t) + max·t in rational arithmetic)"
B6_PEAKS_RE = re.compile(r"\[ADR 0002 #37: the port returns (.+) = (\[\d+\])(?:, [^()]*)? \(([^()]*)\)\]$")
def b6_peaks_port(data, th):
    """The port's find_peaks (Ints), and per candidate (w, T, θ) for the checks."""
    if len(data) < 3 or not (0.0 <= th <= 1.0): return Ints([]), [], None
    lo, hi = math.inf, -math.inf
    for v in data:
        if math.isfinite(v): lo, hi = r_min(lo, v), r_max(hi, v)
    if lo > hi: return Ints([]), [], None
    q = Fraction(lo) * (1 - Fraction(th)) + Fraction(hi) * Fraction(th)
    rt = lo * (1.0 - th) + hi * th
    out, cands = [], []
    for k in range(len(data) - 2):
        w0, w1, w2 = data[k], data[k + 1], data[k + 2]
        if w1 > w0 and w1 > w2:
            above = Fraction(w1) > q if math.isfinite(w1) else w1 > 0.0
            cands.append((w1, rt, above))
            if above: out.append(k + 1)
    return Ints(out), cands, (lo, hi, q, rt)
b6_pk = {"rows": 0, "rows past the checks": 0, "noted": 0, "candidates": 0, "candidates decided by the fast test": 0,
         "candidates where Rust's threshold decides otherwise": 0, "of them subnormal or zero thresholds": 0,
         "of them with min = −DBL_MAX and max = DBL_MAX": 0, "of them other normal near-ties": 0}
U6 = Fraction(1, 2 ** 53)
for ins, outs, where in srows.get(B6_PEAKS, []):
    data, th = peaks_args(ins)
    b6_pk["rows"] += 1
    port, cands, ctx = b6_peaks_port(data, th)
    note = NOTES.get(where, "")
    noted = "[ADR 0002 #37:" in note
    pt = out_toks(port)
    differs = not toks_equal(pt, outs)
    check(noted == differs, f"{where}: find_peaks{(data, th)}: the port's {port} vs the oracle's {outs}: the #37 note is "
          f"{'missing' if differs else 'stray'}")
    if noted:
        b6_pk["noted"] += 1
        mm_ = B6_PEAKS_RE.search(note)
        check(mm_ is not None and mm_.group(1) == B6_PEAKS_HOW and mm_.group(2) == f"[{len(port)}]"
              and toks_equal([parse_tok(t_, where) for t_ in mm_.group(3).split(" ")], pt),
              f"{where}: the #37 note does not name {B6_PEAKS_HOW} = {port}")
    if ctx is None: continue
    b6_pk["rows past the checks"] += 1
    lo, hi, q, rt = ctx
    alo, ahi = abs(lo), abs(hi)
    bound = Fraction(5, 2) * U6 * (1 + U6) * (Fraction(alo) + Fraction(ahi)) + (1 + U6) * Fraction(2) ** -1074
    if math.isfinite(rt):
        check(abs(Fraction(rt) - q) <= bound, f"{where}: Rust's threshold {rt!r} is beyond the proven bound of θ")
    bnd = (alo * 2.0 ** -51 + ahi * 2.0 ** -51) + 2.0 ** -1072
    check(Fraction(bnd) >= bound, f"{where}: the port's bound {bnd!r} is below the proven one")
    all_fast = True
    for w1, _, above in cands:
        b6_pk["candidates"] += 1
        if not math.isfinite(w1):
            b6_pk["candidates decided by the fast test"] += 1
            continue
        fast = None
        if math.isfinite(rt):
            d = w1 - rt
            if d > bnd: fast = True
            elif d < -bnd: fast = False
        if fast is None:
            all_fast = False
        else:
            b6_pk["candidates decided by the fast test"] += 1
            check(fast == above, f"{where}: the fast test decides {w1!r} against {rt!r} wrongly")
        if (w1 > rt) != above:
            b6_pk["candidates where Rust's threshold decides otherwise"] += 1
            if abs(rt) < 2.0 ** -1022 or abs(q) < Fraction(2) ** -1022:
                b6_pk["of them subnormal or zero thresholds"] += 1
            elif lo == -DMAX_V and hi == DMAX_V:
                b6_pk["of them with min = −DBL_MAX and max = DBL_MAX"] += 1
            else:
                b6_pk["of them other normal near-ties"] += 1
    check(not (all_fast and noted), f"{where}: find_peaks{(data, th)}: every candidate is decided by Rust's threshold "
          f"with the proven bound, but the row is noted")
print(f"spectroscopy: ADR 0002 #37: find_peaks: {b6_pk}")
check(all(v > 0 for v in b6_pk.values()), f"ADR 0002 #37: a class has no golden row: {b6_pk}")

# #36, fft_spectrum: every row with a finite signal is recomputed with the port's semantics on mpfr_ref.UX — hisab
# 1.4.0's transform as transcribed (the bit reversal, per stage wn = (cos θ, sin θ) of θ = (−2·π)/len with MPFR's
# correctly rounded cos and sin, w = (1, 0), t = w·d, u ± t, w *= wn), every operation at 53 bits at MPFR's widest
# exponent range, then each bin's √(re·re + im·im) with the squares and the sum at 53 bits and the √ rounded once into
# binary64 (mpfr_ref.ux_sqrt_final) — and must be noted exactly where it differs from the oracle's (its reference: no
# fft row has a correction), the note's tokens being it bit for bit; a signal with a NaN or ±∞ (Rust's IEEE spectrum, B1's
# rule) is never noted. Rows whose f64 chain stays in the normal range are unnoted (B3's strict tracers over every
# product, sum and square of the transform and the magnitudes, the twiddle recurrence included; the √ apart), and so are
# rows where the port's fast test passes (its operand bounds, 2⁻¹⁰⁰ for the twiddles and 2⁻⁹²¹ for the data, the
# retained outputs finite, each bin's square sum in [2⁻⁹⁰⁰, DBL_MAX] or from re = im = 0). The noted rows are counted
# by class, as the survey's.
B6_FFT = "spectroscopy::fft_spectrum"
B6_FFT_HOW = "hisab 1.4.0's transform and the magnitudes evaluated with an unbounded exponent (each √ rounded once)"
B6_FFT_RE = re.compile(r"\[ADR 0002 #36: the port returns (.+) = (\[\d+\])(?:, [^()]*)? \(([^()]*)\)\]$")
def b6_bitrev(n):
    perm, j = list(range(n)), 0
    for i_ in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i_ < j: perm[i_], perm[j] = perm[j], perm[i_]
    return perm
def b6_fft_ux(sig):
    """The port's spectrum of a finite signal (power-of-two length): the transform on UX, the magnitudes rounded once."""
    n = len(sig)
    U_ = MPFR.UX
    with MPFR.unbounded():
        re_ = [U_(x) for x in sig]
        im_ = [U_(0.0) for _ in sig]
        perm = b6_bitrev(n)
        re_, im_ = [re_[k] for k in perm], [im_[k] for k in perm]
        length = 2
        while length <= n:
            half = length >> 1
            angle = -2.0 * math.pi / float(length)
            wr0, wi0 = U_(MPFR.cos(angle)), U_(MPFR.sin(angle))
            for start in range(0, n, length):
                wr, wi = U_(1.0), U_(0.0)
                for k in range(half):
                    a_, b_ = start + k, start + k + half
                    tr = wr * re_[b_] - wi * im_[b_]
                    ti = wr * im_[b_] + wi * re_[b_]
                    ur, ui = re_[a_], im_[a_]
                    re_[a_], im_[a_] = ur + tr, ui + ti
                    re_[b_], im_[b_] = ur - tr, ui - ti
                    wr, wi = wr * wr0 - wi * wi0, wr * wi0 + wi * wr0
            length <<= 1
        return [MPFR.ux_sqrt_final(re_[k] * re_[k] + im_[k] * im_[k]) for k in range(n // 2 + 1)]
def b6_fft_trace(sig):
    """(in range, fast): Rust's IEEE chain stays in the normal range (every product of nonzero factors normal and above
    2⁻¹⁰²², every sum and product finite; the √ apart), and the port's fast test passes (_km_fft's operand bounds, the
    retained outputs finite, each bin's square sum in [2⁻⁹⁰⁰, DBL_MAX] or from re = im = 0)."""
    n = len(sig)
    ok, fast = [True], [True]
    LW, LD, LO, HI = 2.0 ** -100, 2.0 ** -921, 2.0 ** -900, DMAX_V
    def mul(x, y):
        p_ = x * y
        if x != 0.0 and y != 0.0 and not (math.isfinite(p_) and abs(p_) > 2.0 ** -1022): ok[0] = False
        if not math.isfinite(p_): ok[0] = False
        return p_
    def add(x, y):
        z = x + y
        if not math.isfinite(z): ok[0] = False
        return z
    def chk(v, lim):
        if 0.0 < abs(v) < lim: fast[0] = False
    perm = b6_bitrev(n)
    re_, im_ = [sig[k] for k in perm], [0.0] * n
    length = 2
    while length <= n:
        half = length >> 1
        angle = -2.0 * math.pi / float(length)
        wr0, wi0 = MPFR.cos(angle), MPFR.sin(angle)
        chk(wr0, LW)
        chk(wi0, LW)
        for start in range(0, n, length):
            wr, wi = 1.0, 0.0
            for k in range(half):
                a_, b_ = start + k, start + k + half
                chk(re_[b_], LD)
                chk(im_[b_], LD)
                if start == 0:
                    chk(wr, LW)
                    chk(wi, LW)
                tr = add(mul(wr, re_[b_]), -mul(wi, im_[b_]))
                ti = add(mul(wr, im_[b_]), mul(wi, re_[b_]))
                ur, ui = re_[a_], im_[a_]
                re_[a_], im_[a_] = add(ur, tr), add(ui, ti)
                re_[b_], im_[b_] = add(ur, -tr), add(ui, -ti)
                wr, wi = add(mul(wr, wr0), -mul(wi, wi0)), add(mul(wr, wi0), mul(wi, wr0))
        length <<= 1
    for k in range(n // 2 + 1):
        rr, ii = re_[k], im_[k]
        if not (math.isfinite(rr) and math.isfinite(ii)): fast[0] = False
        sq_ = add(mul(rr, rr), mul(ii, ii))
        if not (LO <= sq_ <= HI or (rr == 0.0 and ii == 0.0)): fast[0] = False
    return ok[0], fast[0]
b6_ff = {"rows": 0, "finite signals": 0, "noted": 0, "noted bins": 0, "Rust's chain in range": 0,
         "the port's fast test passing": 0}
b6_ffc = {"qual: Rust 0, ±∞ or NaN where the value is representable": 0,
          "NaN-for-∞: Rust's NaN where the value is ±∞": 0,
          "degraded: the value normal, Rust's more than 2 ulp off": 0, "within 2 ulp, or the value subnormal": 0}
for ins, outs, where in srows.get(B6_FFT, []):
    sig = fl(ins)[1:]
    b6_ff["rows"] += 1
    note = NOTES.get(where, "")
    noted = "[ADR 0002 #36:" in note
    n = len(sig)
    if n == 0 or n & (n - 1) or not all(map(math.isfinite, sig)):
        check(not noted, f"{where}: an fft_spectrum row with a non-finite or invalid signal is noted")
        continue
    b6_ff["finite signals"] += 1
    inr, fast = b6_fft_trace(sig)
    b6_ff["Rust's chain in range"] += inr
    b6_ff["the port's fast test passing"] += fast
    check(not (inr and noted), f"{where}: fft_spectrum of {n}: Rust's chain stays in the normal range, but the row is noted")
    check(not (fast and noted), f"{where}: fft_spectrum of {n}: the port's fast test passes, but the row is noted")
    if (inr or fast) and n > 1024:
        continue                                        # the lemma's rows (Rust's bits), too long for UX here
    port = b6_fft_ux(sig)
    pt = out_toks(port)
    differs = not toks_equal(pt, outs)
    check(noted == differs, f"{where}: fft_spectrum of {n}: the port's spectrum vs the oracle's: the #36 note is "
          f"{'missing' if differs else 'stray'}")
    if not noted: continue
    b6_ff["noted"] += 1
    mm_ = B6_FFT_RE.search(note)
    check(mm_ is not None and mm_.group(1) == B6_FFT_HOW and mm_.group(2) == f"[{len(port)}]"
          and toks_equal([parse_tok(t_, where) for t_ in mm_.group(3).split(" ")], pt),
          f"{where}: the #36 note does not name the port's spectrum")
    worst = 3
    for (gk, gv, gb), x in zip(outs[1:], port):
        if gb == f2b(x): continue
        b6_ff["noted bins"] += 1
        o = gv
        if math.isnan(o) or o == 0.0 or math.isinf(o):
            c_ = 0 if (math.isfinite(x) and (x != 0.0 or o != 0.0)) else (1 if math.isnan(o) and math.isinf(x) else 3)
        elif x >= 2.0 ** -1022 and math.isfinite(x) and abs(ulp_dist(gb, f2b(x))) > 2:
            c_ = 2
        else:
            c_ = 3
        worst = min(worst, c_)
    b6_ffc[list(b6_ffc)[worst]] += 1
print(f"spectroscopy: ADR 0002 #36: fft_spectrum: {b6_ff}")
print(f"ADR 0002 #36 noted rows by their worst bin: {b6_ffc}")
check(b6_ff["noted"] > 0 and b6_ffc["qual: Rust 0, ±∞ or NaN where the value is representable"] > 0,
      f"ADR 0002 #36: a class has no golden row: {b6_ff} {b6_ffc}")

# ── 25. The overflow audit's B7 (ADR 0002 #39): balance_equation on unbounded integers ─────────────────────────────────
# Every balance_equation row with an integer counterpart (ADR 0002 #7) is recomputed with the port's semantics,
# independently of the generator: Rust's algorithm (st_balance_int's steps: the pivot search, the swap, reduce_row, the
# sign flip, row_r ← pivot·row_r − factor·row_pivot then reduce_row, the lcm of the pivots, x[c] = −(row[f]·(lcm/piv)),
# reduce_row(x), the zero, sign and u32 checks, the verification) on Python ints with no i64 check, b7_port; and that
# result checked against the exact null space (st_nullspace, Fraction): an Ok result is the primitive positive basis
# vector of a one-dimensional null space, an Err a null space of another dimension or a basis vector with a zero, mixed
# signs or an entry past u32. A row is noted exactly where the port's result differs from the oracle's, the note's
# tokens being the port's whole result. Strict tracer: a noted row is one where Rust's checked arithmetic overflows
# (b7_rust_site: the site of its first overflow, Rust's ComputationError), and every row where Rust does not overflow is
# unnoted; the rows where Rust overflows and the exact result is still an Err are counted by class. The proofs behind
# the port's tables (ADR 0002 #39) are checked on every row: every entry after reduce_row, every original entry, the lcm
# and x are below 2^B, every product and difference before reduce_row below 2^(2B + 1), B the port's
# _km_balance_bits; before reduce_row(x) the gcd of x is 1 wherever no coefficient is zero (so |x| ≥ lcm, and an
# overflow of the lcm or of x is always an Err past u32); a reduced row of gcd 2⁶³ has every nonzero entry −2⁶³.
B7_HOW = "Rust's fraction-free elimination on unbounded integers"
B7_RE = re.compile(r"\[ADR 0002 #39: the port returns (.+) = (\[\d+\])(?:, [^()]*)? \(([^()]*)\)\]$")
def b7_bits(a, m, n):
    """The port's size bound B (src/stoichiometry.cyr's _km_balance_bits): the smaller of the sums of the q = min(m, n)
    largest row and column b's, b = bitlen(max |entry|) + ⌈bitlen(nonzero count)/2⌉ (0 for a zero line)."""
    def bs(lines):
        out = []
        for ln in lines:
            nz = [abs(v) for v in ln if v]
            out.append(max(nz).bit_length() + (len(nz).bit_length() + 1) // 2 if nz else 0)
        return sum(sorted(out, reverse=True)[:min(m, n)])
    rows_ = [a[i * n:(i + 1) * n] for i in range(m)]
    return min(bs(rows_), bs([[a[i * n + j] for i in range(m)] for j in range(n)]))
def b7_port(a, m, n, B=None, tr=None):
    """The port's balance_equation (ADR 0002 #39) on the integer matrix a: Rust's steps on Python ints. With B, every
    value is checked against the bound; tr (a dict) receives the largest value's bit length, and the gcd of x before its
    reduce_row."""
    if m < 1 or n < 2 or any(abs(v) > ST_MAX_ENTRY for v in a): return Err("InvalidInput")
    A = [list(a[i * n:(i + 1) * n]) for i in range(m)]
    big = [0, 0]
    def stored(v):
        big[0] = max(big[0], abs(v).bit_length())
        if B is not None: check(abs(v) < 2 ** B, f"balance_equation: a stored value of {abs(v).bit_length()} bits past 2^{B}")
        return v
    def red(row):
        g = 0
        for v in row: g = math.gcd(g, abs(v))
        return [stored(v // g) for v in row] if g > 1 else [stored(v) for v in row]
    for row in A: [stored(v) for v in row]
    rank, free = 0, None
    for col in range(n):
        p = next((r for r in range(rank, m) if A[r][col] != 0), None)
        if p is None:
            if free is not None: return Err("ComputationError")
            free = col
            continue
        A[p], A[rank] = A[rank], A[p]
        A[rank] = red(A[rank])
        if A[rank][col] < 0: A[rank] = [-v for v in A[rank]]
        piv = A[rank][col]
        for r in range(m):
            fac = A[r][col]
            if r == rank or fac == 0: continue
            new = []
            for j in range(n):
                d = piv * A[r][j] - fac * A[rank][j]
                big[1] = max(big[1], abs(piv * A[r][j]).bit_length(), abs(fac * A[rank][j]).bit_length(),
                             abs(d).bit_length())
                if B is not None:
                    check(max(abs(piv * A[r][j]), abs(fac * A[rank][j]), abs(d)) < 2 ** (2 * B + 1),
                          f"balance_equation: a product or difference past 2^(2·{B} + 1)")
                new.append(d)
            A[r] = red(new)
        rank += 1
    if free is None: return Err("ComputationError")
    lcm = 1
    for k in range(rank):
        pv = next(v for v in A[k] if v != 0)
        lcm = stored(lcm // math.gcd(lcm, abs(pv)) * pv)
    x = [0] * n
    x[free] = lcm
    for k in range(rank):
        c = next(j for j, v in enumerate(A[k]) if v != 0)
        x[c] = stored(-(A[k][free] * (lcm // A[k][c])))
    g = 0
    for v in x: g = math.gcd(g, abs(v))
    if tr is not None: tr["x gcd"] = g; tr["bits"] = big[0]; tr["pbits"] = big[1]
    x = red(x)
    if 0 in x: return Err("ComputationError")
    pos = sum(1 for v in x if v > 0)
    if pos == 0: x = [-v for v in x]
    elif pos != n: return Err("ComputationError")
    if any(v < 0 or v > U32_MAX for v in x): return Err("ComputationError")
    for i in range(m):
        if sum(a[i * n + j] * x[j] for j in range(n)) != 0: return Err("ComputationError")
    return x
def b7_rust_site(a, m, n):
    """The site of Rust's first i64 overflow (its checked operations in its order), or None."""
    def ok(v): return I64_MIN <= v <= I64_MAX
    A = [list(a[i * n:(i + 1) * n]) for i in range(m)]
    def red(row, site):
        g = 0
        for v in row: g = math.gcd(g, abs(v))
        if g > I64_MAX: return None, site
        return ([tdiv(v, g) for v in row] if g > 1 else row), None
    rank, free = 0, None
    for col in range(n):
        p = next((r for r in range(rank, m) if A[r][col] != 0), None)
        if p is None:
            if free is not None: return None
            free = col
            continue
        A[p], A[rank] = A[rank], A[p]
        A[rank], s_ = red(A[rank], "a reduced row's gcd of 2⁶³")
        if s_: return s_
        if A[rank][col] < 0:
            if I64_MIN in A[rank]: return "the pivot row's negation"
            A[rank] = [-v for v in A[rank]]
        piv = A[rank][col]
        for r in range(m):
            fac = A[r][col]
            if r == rank or fac == 0: continue
            new = []
            for j in range(n):
                if not ok(piv * A[r][j]): return "pivot·a"
                if not ok(fac * A[rank][j]): return "factor·a"
                if not ok(piv * A[r][j] - fac * A[rank][j]): return "the difference"
                new.append(piv * A[r][j] - fac * A[rank][j])
            A[r], s_ = red(new, "a reduced row's gcd of 2⁶³")
            if s_: return s_
        rank += 1
    if free is None: return None
    lcm = 1
    for k in range(rank):
        pv = next(v for v in A[k] if v != 0)
        if math.gcd(lcm, abs(pv)) > I64_MAX: return "the lcm's gcd"
        lcm = tdiv(lcm, math.gcd(lcm, abs(pv))) * pv
        if not ok(lcm): return "the lcm's product"
    x = [0] * n
    x[free] = lcm
    for k in range(rank):
        c = next(j for j, v in enumerate(A[k]) if v != 0)
        px = A[k][free] * tdiv(lcm, A[k][c])
        if not ok(px): return "x's product"
        if not ok(-px): return "x's negation"
        x[c] = -px
    g = 0
    for v in x: g = math.gcd(g, abs(v))
    if g > I64_MAX: return "reduce_row(x)'s gcd"
    if 0 in x: return None
    if sum(1 for v in x if v > 0) == 0 and I64_MIN in x: return "the final sign flip"
    return None
b7 = {"rows": 0, "noted": 0, "Rust overflows": 0}
b7_still_err, b7_sites, b7_shapes, b7_shapes_before, b7_bits_max, b7_pbits_max, b7_gcd1 = {}, {}, {}, {}, 0, 0, 0
b7_still_before, b7_ovf_before = 0, 0
for idx, (ins, outs, where) in enumerate(st_rows.get("stoichiometry::balance_equation", [])):
    mat = st_matrix(ins)
    cp = st_counterpart(mat)
    note = NOTES.get(where, "")
    noted = "[ADR 0002 #39:" in note
    if cp is None:
        check(not noted, f"{where}: a #39 note on a matrix without an integer counterpart")
        continue
    m_, n_ = len(mat), len(mat[0]) if mat else 0
    b7["rows"] += 1
    rust = st_balance_int(cp, m_, n_)
    check(toks_equal(outs, st_toks("stoichiometry::balance_equation", rust)), f"{where}: Rust's balance recomputed")
    if m_ < 1 or n_ < 2 or any(abs(v) > ST_MAX_ENTRY for v in cp):
        check(not noted, f"{where}: a #39 note on a row Rust's checks refuse")
        continue
    tr = {}
    port = b7_port(cp, m_, n_, B=b7_bits(cp, m_, n_), tr=tr)
    b7_bits_max = max(b7_bits_max, tr.get("bits", 0))
    b7_pbits_max = max(b7_pbits_max, tr.get("pbits", 0))
    dim, w = st_nullspace(cp, m_, n_)
    valid = dim == 1 and all(q > 0 for q in w) and max(w) <= U32_MAX
    check(port == w if valid else isinstance(port, Err) and port.name == "ComputationError",
          f"{where}: the port's {port} is not the exact null space's (dimension {dim}, {w})")
    if dim == 1 and w is not None and 0 not in w:
        check(tr.get("x gcd") == 1, f"{where}: x's gcd before reduce_row is {tr.get('x gcd')}, not 1")
        b7_gcd1 += 1
    pt = st_toks("stoichiometry::balance_equation", port)
    differs = not toks_equal(pt, outs)
    check(noted == differs, f"{where}: the port's {port} vs the oracle's {outs}: the #39 note is "
          f"{'missing' if differs else 'stray'}")
    site = b7_rust_site(cp, m_, n_)
    if site is not None:
        b7["Rust overflows"] += 1
        b7_ovf_before += idx < 3329
        b7_still_before += idx < 3329 and isinstance(port, Err)
        b7_sites.setdefault(site, {"rows": 0, "Ok exactly": 0})
        b7_sites[site]["rows"] += 1
        b7_sites[site]["Ok exactly"] += not isinstance(port, Err)
        check(rust == Err("ComputationError"), f"{where}: Rust overflows but its result is {rust}")
        if isinstance(port, Err):
            cls = "ambiguous or trivial" if dim != 1 else ("zero" if 0 in w else ("mixed" if not (all(q > 0 for q in w)) else "past u32"))
            b7_still_err[cls] = b7_still_err.get(cls, 0) + 1
    else:
        check(not noted and toks_equal(pt, outs), f"{where}: Rust does not overflow, yet the port's result differs")
    if site == "a reduced row's gcd of 2⁶³" or site in ("the lcm's product", "x's product", "x's negation"):
        check(isinstance(port, Err), f"{where}: an overflow at {site} with an exact Ok")
    if not noted: continue
    b7["noted"] += 1
    check(site is not None, f"{where}: a noted row where Rust does not overflow")
    mm_ = B7_RE.search(note)
    check(mm_ is not None and mm_.group(1) == B7_HOW and mm_.group(2) == f"[{len(port)}]"
          and toks_equal([parse_tok(t_, where) for t_ in mm_.group(3).split(" ")], pt),
          f"{where}: the #39 note does not name the port's result {port}")
    b7_shapes[(m_, n_)] = b7_shapes.get((m_, n_), 0) + 1
    if idx < 3329: b7_shapes_before[(m_, n_)] = b7_shapes_before.get((m_, n_), 0) + 1
print(f"stoichiometry: ADR 0002 #39: balance_equation: {b7}; the noted rows Rust's ComputationError, the port's Ok "
      f"(every one checked against the Fraction null space); the exact elimination's largest stored value {b7_bits_max} "
      f"bits, its largest product or difference before reduce_row {b7_pbits_max} bits")
print(f"  the rows before B7: {sum(b7_shapes_before.values())} noted, by (m, n): "
      f"{dict(sorted(b7_shapes_before.items(), key=lambda kv: (-kv[1], kv[0])))}")
print(f"  every noted row by (m, n): {dict(sorted(b7_shapes.items()))}")
print(f"  Rust's first overflow by site (rows, of them Ok exactly): "
      + "; ".join(f"{k} {v['rows']} ({v['Ok exactly']})" for k, v in sorted(b7_sites.items())))
print(f"  where Rust overflows and the exact result is still the ComputationError: {dict(sorted(b7_still_err.items()))}; "
      f"x's gcd 1 before reduce_row on all {b7_gcd1} rows with a one-dimensional null space and no zero coefficient")
print(f"  of the rows before B7: Rust overflows on {b7_ovf_before}, the exact result Ok on "
      f"{b7_ovf_before - b7_still_before} (noted) and still the ComputationError on {b7_still_before} (unnoted)")
check(sum(b7_shapes_before.values()) == 104 and b7_shapes_before.get((2, 3)) == 55,
      f"ADR 0002 #39: {sum(b7_shapes_before.values())} rows before B7 noted, not the survey's 104 (55 of them 2 × 3)")
check(all(b7_sites.get(k, {}).get("Ok exactly", 0) > 0 for k in ("pivot·a", "factor·a", "the difference",
      "the pivot row's negation")), f"ADR 0002 #39: a site that can lead to Ok has no golden row: {b7_sites}")
check(all(b7_sites.get(k, {}).get("rows", 0) > 0 for k in ("a reduced row's gcd of 2⁶³", "the lcm's product",
      "x's product", "x's negation")), f"ADR 0002 #39: an overflow site has no golden row: {b7_sites}")
check(all(b7_still_err.get(k, 0) > 0 for k in ("ambiguous or trivial", "zero", "mixed", "past u32")),
      f"ADR 0002 #39: an Err class after Rust's overflow has no golden row: {b7_still_err}")
check(b7_pbits_max >= 3000, f"ADR 0002 #39: the largest product before reduce_row is only {b7_pbits_max} bits")

VALIDATED = {"constants", "math", "element", "molecule", "reaction", "thermo", "phase", "gas", "thermochem", "solution",
             "electrochemistry", "kinetics", "organic", "inorganic", "nuclear", "potential", "spectroscopy", "stoichiometry"}
check(set(tsvs) == VALIDATED, f"golden modules with no section here: {sorted(set(tsvs) - VALIDATED)}; "
      f"sections with no module: {sorted(VALIDATED - set(tsvs))}")

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
