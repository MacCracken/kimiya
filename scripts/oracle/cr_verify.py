"""Verify the correctly rounded exp, exp10, expm1, ln, log10, pow, sin, cos, atan2, acos, x·2^y and 1/(100·x) of
src/mathx.cyr against MPFR (ADR 0006), and its unbounded-exponent final operations with potential's geometry (ADR 0002
#20); and the exact rationals of spectroscopy (#21, #22), stoichiometry (#24 ... #26) and the overflow audit's B1 (#27,
#28: the correctly rounded monomial, thermo's, phase's and kinetics' rationals) and B2 (#29, #30: solution's,
electrochemistry's and spectroscopy's, and the correctly rounded square root of a quotient) against Fraction, exact
integer square roots and MPFR.

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
h100, _km_cr_hundredth_div(x) = 1/(100*x) rounded once (M4e, ADR 0002 #19: spectroscopy's wavenumber conversions), gets
>= 1,000,000 arguments: every binade of both signs, the fast path's edges 2^+-900, the overflow end (the subnormal x
whose 1/(100*x) passes DBL_MAX) and the subnormal-result end (x near 2^1022/100^-1), the integers, decimal wavelengths
k*10^-j, and the special values; the reference is mpfr_ref.hundredth_div.
M4e-2 (ADR 0002 #20) adds the unbounded-exponent final operations of src/mathx.cyr, each on (significand, exponent) pairs
far outside binary64's range: uxsqrt (_km_ux_sqrt_f64), uxmul, uxdiv, uxadd (_km_ux_mul_f64, _km_ux_div_f64,
_km_ux_add_f64) and uxatan2 (_km_ux_atan2_f64, its accurate path forced through _km_ux_atan2_acc), with subnormal results
and their ties (short significands), the overflow end, sticky exponent gaps, cancellations and every exponent difference;
and potential's distance_3d, bond_angle and dihedral_angle (dist, bond, dih) on edge geometries (random points at
scales 2^-1074 ... 2^1018, mixed magnitudes, subnormal coordinates, near coincidence and collinearity at scale,
overflowing differences), the public function and the forced slow path (_km_*_ux) both. The references are mpfr_ref's
UX (53 bits at MPFR's widest exponent range; the final operation by mpfr_check_range and mpfr_subnormalize) and
mpfr_ref.potential_* (Rust's sequence on UX).
M4e-3 (ADR 0002 #21, #22) adds the exact rationals rounded once: ratio, _km_cr_ratio(p, q, e) = RN(p/q·2^e) for multi-word
integers p, q < 2^250 (random sizes, every exponent of the result, subnormal results and their ties, exact quotients and
midpoints, the overflow end); bohr (bohr_energy_level: every Z with n = 2^k·Z, small n exhaustively, n at the 2^13 and
2^32 ends, random n of every bit length), trans and ryd (bohr_transition_energy, rydberg_wavelength: random (Z, n1, n2)
at every scale, n2 = n1 + 1 … up to 2^32 − 1, the fast path's 2^13 edge, every transition tie and the exact values with
n2 <= 120), series (spectral_series_wavelength: every series, lines 1 … 2^32 − 6) and rmass (reduced_mass: binade-uniform
pairs, equal masses, subnormal masses and their ties, three families of normal ties at random scales, exponent gaps
45 … 64 around the proven 55, the fast path's 2^-1021 edge, the overflow end): the public function, its fast path and its
exact path forced. Their reference is exact rational arithmetic (Fraction, rounded by CPython's correctly rounded int
division).
M4f-2 (ADR 0002 #24 ... #26) adds stoichiometry's exact rationals and comparisons, against Fraction as well: muldiv
(_km_cr_mul_div = RN(x*y/z): binade-uniform triples, the fast path's 2^+-1021 edges and past, exact quotients and
midpoints, values within 2^-100 of a midpoint on both sides, subnormal results and ties; the fast path and
_km_mul_div_acc forced), pyield, tyield, pcomp (percent_yield, theoretical_yield_moles with two reactants,
percent_composition's first entry: the public function, its fast path and the exact path forced), limit
(limiting_reagent: the public function and the exact comparisons alone, on ties, near-ties and ratios that both
overflow or underflow) and empir (empirical_formula, and _km_empirical_formula with every tolerance test exact:
formulas' percentages, r*k at n +- 1/10 within ulps and exactly, subnormal percentages, ratios near 2^32).
The overflow audit's B1 (ADR 0002 #27, #28) adds mono (_km_cr_mono = RN(prod num / prod den), every shape up to four
factors a side and five or six with short ones, both signs, zeros, the result at the range ends, 53-bit midpoints of
the 3/0, 2/1, 2/2 and 3/1 shapes and the subnormal grid's ties, the domain's NaNs; its fast path and its exact path
forced), monotie (_km_mono_tie: the exact sign decision for z1 < z2 = nextup(z1) against Fraction, F near, far and
exactly at their midpoint), and the ten functions: heat, shc, boil, frz, osm, zhl, ssi, pre (heat_transfer,
specific_heat_from_calorimetry, boiling_point_elevation, freezing_point_depression, osmotic_pressure,
zero_order_half_life, steady_state_intermediate, pre_equilibrium_rate: the public function, its shape's fast path, the
monomial's exact path forced; references Fraction) and soc, mm (second_order_concentration,
michaelis_menten: the public function, its Rust-chain fast path where the range test passes, the unbounded-exponent slow
path forced; reference mpfr_ref.UX), each on binade-uniform, typical, partial-products-past-the-range, midpoint,
subnormal and Err arguments.
The overflow audit's B2 (ADR 0002 #29, #30) adds sqrtdiv (_km_cr_sqrt_div = RN(sqrt(x/y)): every binade of both, x/y
past either end, subnormal operands and results, exact roots, roots within ulps and within 2^-100 of a 53-bit midpoint
(continued-fraction convergents of m^2), the domain; the fast path and _km_sqrt_div_acc forced; reference an exact
integer square root) and the twelve functions: dil, cond, fmass, fcharge, fratio, abs, conc, rotor, rot, moment
(dilution, solution_conductivity, faraday_mass_deposited, faraday_charge_required, faraday_mass_ratio, absorbance,
concentration_from_absorbance, rigid_rotor_energy, rotational_constant, diatomic_moment_of_inertia: the public function,
its shape's fast path, the monomial's exact path forced; reference Fraction), vib (vibrational_frequency: the root's two
paths) and nernst (nernst_potential: the public function, its Rust-chain fast path where the range test passes, the
unbounded-exponent slow path forced; reference mpfr_ref.UX with MPFR's log and ux_sub_final), each on binade-uniform,
typical, past-the-range, midpoint, subnormal and Err arguments.
Through each function's _km_cr_<f>_fast entry it reports the fast path's fallback rate per region, and
checks the fast path's result equals the function's wherever its test passes. It checks that the three
targets give identical bits and flags.

Usage, from the repo root:
    python3 -I scripts/oracle/cr_verify.py                 # x86_64 ELF, x86_64 PE (wine), aarch64 (qemu)
    python3 -I scripts/oracle/cr_verify.py --targets elf   # one target
    python3 -I scripts/oracle/cr_verify.py --big           # also 10,000,000 more per primitive on ELF
    python3 -I scripts/oracle/cr_verify.py --only pow,exp  # some functions (exp exp10 expm1 ln log10 pow sin cos
                                                           # atan2 acos exp2s h100 uxsqrt uxmul uxdiv uxadd
                                                           # uxatan2 dist bond dih ratio bohr trans ryd series
                                                           # rmass muldiv pyield tyield pcomp limit empir mono
                                                           # heat shc boil frz osm zhl ssi pre soc mm monotie
                                                           # sqrtdiv dil cond fmass fcharge fratio nernst abs
                                                           # conc vib rotor rot moment)
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
    (12, "h100", "_km_cr_hundredth_div", "_km_cr_hundredth_div_acc", mpfr_ref.hundredth_div),
    # ADR 0002 #20 (M4e-2): the unbounded-exponent final operations ((significand, exponent) pairs in, binary64 out) and
    # potential's geometry by its forced slow path; the references are below (ref_ux)
    (13, "uxsqrt", "_km_ux_sqrt_f64", "_km_ux_sqrt_f64", None),
    (14, "uxmul", "_km_ux_mul_f64", "_km_ux_mul_f64", None),
    (15, "uxdiv", "_km_ux_div_f64", "_km_ux_div_f64", None),
    (16, "uxadd", "_km_ux_add_f64", "_km_ux_add_f64", None),
    (17, "uxatan2", "_km_ux_atan2_f64", "_km_ux_atan2_acc", None),
    (18, "dist", "distance_3d", "_km_distance_3d_ux", None),
    (19, "bond", "bond_angle", "_km_bond_angle_ux", None),
    (20, "dih", "dihedral_angle", "_km_dihedral_angle_ux", None),
    # ADR 0002 #21, #22 (M4e-3): the exact rationals rounded once, against Fraction (ref_spec, below)
    (26, "ratio", "_km_cr_ratio", "_km_cr_ratio", None),
    (21, "bohr", "bohr_energy_level", "_km_bohr_level_acc", None),
    (22, "trans", "bohr_transition_energy", "_km_bohr_transition_acc", None),
    (23, "ryd", "rydberg_wavelength", "_km_rydberg_acc", None),
    (24, "series", "spectral_series_wavelength", "_km_rydberg_acc", None),
    (25, "rmass", "reduced_mass", "_km_reduced_mass_acc", None),
]
ARITY = {5: 2, 8: 2, 11: 2, 13: 2, 14: 4, 15: 4, 16: 4, 17: 4, 18: 6, 19: 9, 20: 12,
         21: 2, 22: 3, 23: 3, 24: 2, 25: 2, 26: 9}   # argument words; 1 elsewhere
UX_NAMES = ("uxsqrt", "uxmul", "uxdiv", "uxadd", "uxatan2", "dist", "bond", "dih")
SPEC_NAMES = ("ratio", "bohr", "trans", "ryd", "series", "rmass")
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


def args_h100(rng):
    """Arguments for 1/(100*x) (M4e)."""
    sp = [0.0, -0.0, INF, -INF, NAN, fb(1), -fb(1), fb(0x0010000000000000), fb(0x7FEFFFFFFFFFFFFF),
          -fb(0x7FEFFFFFFFFFFFFF), 1.0, -1.0, 17.0, 29.0, 100.0, 0.01, 0.5e-6, 1e-5]
    edges = []
    for x in (2.0 ** 900, 2.0 ** -900, fb(0x7FEFFFFFFFFFFFFF), 1.0 / (100.0 * fb(0x7FEFFFFFFFFFFFFF)),
              1.0 / (100.0 * 2.0 ** -1022), 2.0 ** -1022, 1.0 / 100.0):
        edges += around(x, 3000)
        edges += [-v for v in around(x, 500)]
    ints = [float(k) for k in range(1, 150_001)]
    dec = [rng.randint(1, 9999) * 10.0 ** -rng.randint(0, 22) for _ in range(150_000)]
    return [
        ("special", sp),
        ("fast-path edges, overflow and subnormal-result ends", edges),
        ("integers 1 .. 150000", ints),
        ("decimal wavelengths k*10^-j", dec),
        ("binade-uniform", binade_uniform(rng, 700_000, -1074, 1023)),
    ]


def _sig(rng, neg_ok=True):
    """A random significand in [1, 2), of either sign."""
    x = fb((0x3FF << 52) | rng.randrange(0, 1 << 52))
    return -x if neg_ok and rng.random() < 0.5 else x


def _short_sig(rng):
    """A significand of at most 12 bits: exact products, quotients and square roots, and their ties."""
    k = rng.randint(1 << 11, (1 << 12) - 1)
    return k / 2048.0


def args_uxsqrt(rng):
    """(s, e) for RN(√(s·2^e)) (M4e-2): subnormal results and their ties, the overflow end, exact squares, every e."""
    sub = [(_sig(rng, False), rng.randint(-2200, -2040)) for _ in range(60_000)]
    top = [(_sig(rng, False), rng.randint(2040, 2050)) for _ in range(20_000)]
    exact = []
    for _ in range(40_000):            # (k/2048)² with k < 2^12, so the root is exact: ties at the subnormal grid
        r = _short_sig(rng)
        e = rng.randint(-2200, -2100) if rng.random() < 0.7 else rng.randint(-2000, 2000)
        q, ee = r * r, 2 * e
        if q >= 2.0:
            q, ee = q / 2.0, ee + 1
        exact.append((q, ee))
    every = [(_sig(rng, False), rng.randint(-5000, 5000)) for _ in range(80_000)]
    sp = [(1.0, -2150), (1.0, -2148), (1.0, 2048), (fb(0x3FFFFFFFFFFFFFFF), 2047), (1.0, 2046), (0.0, 0), (-0.0, 0),
          (1.5, -2149), (1.0 + 2.0 ** -52, -2150)]
    return [("special", sp), ("subnormal results", sub), ("the overflow end", top), ("exact squares", exact),
            ("every exponent", every)]


def _ux_pairs(rng, n, targets, exact=False, combine="mul"):
    """(s1, e1, s2, e2) whose product (combine = mul) or quotient (div) lands at 2^t for t in targets (a range)."""
    out = []
    for _ in range(n):
        t = rng.randint(*targets)
        e1 = rng.randint(-5000, 5000)
        s1 = _short_sig(rng) if exact else _sig(rng)
        s2 = _short_sig(rng) if exact else _sig(rng)
        if exact and rng.random() < 0.5:
            s1 = -s1
        e2 = t - e1 if combine == "mul" else e1 - t
        out.append((s1, e1, s2, e2))
    return out


def args_uxmul(rng):
    return [("special", [(1.5, -1074, 1.0, -1), (1.0, -1075, 1.0, 0), (1.5, 1023, 1.5, 0), (0.0, 0, 1.5, 9000),
                         (-0.0, 0, 1.5, -9000), (1.0, 2000, 1.0, -1000)]),
            ("subnormal results", _ux_pairs(rng, 70_000, (-1130, -1020))),
            ("subnormal ties (short significands)", _ux_pairs(rng, 40_000, (-1090, -1060), exact=True)),
            ("the overflow end", _ux_pairs(rng, 20_000, (1020, 1025))),
            ("normal results", _ux_pairs(rng, 40_000, (-1020, 1020))),
            ("below 2^-1130", _ux_pairs(rng, 10_000, (-3000, -1131)))]


def args_uxdiv(rng):
    return [("special", [(1.0, -1000, 1.5, 74), (1.0, 2000, 1.0, 977), (0.0, 0, 1.5, 5), (1.5, 5, 0.0, 0),
                         (-1.25, -1074, 1.0, 1)]),
            ("subnormal results", _ux_pairs(rng, 70_000, (-1130, -1020), combine="div")),
            ("subnormal ties (short significands)", _ux_pairs(rng, 40_000, (-1090, -1060), exact=True, combine="div")),
            ("the overflow end", _ux_pairs(rng, 20_000, (1020, 1025), combine="div")),
            ("normal results", _ux_pairs(rng, 40_000, (-1020, 1020), combine="div")),
            ("below 2^-1130", _ux_pairs(rng, 10_000, (-3000, -1131), combine="div"))]


def args_uxadd(rng):
    """(s1, e1, s2, e2) for RN(x + y): results near and in the subnormal range, exponent gaps past the 55 where the
    smaller term becomes a sticky bit, cancellations, the overflow end."""
    def gap(lo, hi, n, emin, emax):
        out = []
        for _ in range(n):
            e1 = rng.randint(emin, emax)
            out.append((_sig(rng), e1, _sig(rng), e1 - rng.randint(lo, hi)))
        return out
    canc = []
    for _ in range(30_000):
        s1 = _sig(rng)
        e1 = rng.randint(-1200, 1200)
        s2 = -ulps(abs(s1), rng.randint(-2000, 2000)) * math.copysign(1.0, s1)
        canc.append((s1, e1, s2 if 1.0 <= abs(s2) < 2.0 else -s1, e1 - rng.randint(0, 1)))
    ties = []
    for _ in range(30_000):             # 2^-1075-scale midpoints and a smaller term that breaks them
        e1 = rng.randint(-1080, -1070)
        ties.append((_short_sig(rng), e1, _sig(rng), e1 - rng.randint(50, 200)))
    sp = [(1.0, -1075, 1.0, -1200), (1.0, -1075, -1.0, -1200), (1.0, -1075, 1.0, -1075), (1.0, 1023, 1.0, 1023),
          (fb(0x3FFFFFFFFFFFFFFF), 1023, 1.0, 970), (fb(0x3FFFFFFFFFFFFFFF), 1023, 1.0, 969), (1.0, 0, -1.0, 0),
          (0.0, 0, -0.0, 0), (-0.0, 0, -0.0, 0), (0.0, 0, 1.5, -1100), (-1.5, -1100, 0.0, 0)]
    return [("special", sp), ("near the subnormal range", gap(0, 70, 60_000, -1130, -1015)),
            ("sticky gaps (56 ... 200)", gap(56, 200, 30_000, -1100, 1023)), ("cancellations", canc),
            ("midpoints and a breaking term", ties), ("the overflow end", gap(0, 60, 20_000, 1020, 1024)),
            ("every exponent", gap(0, 120, 30_000, -5000, 5000))]


def args_uxatan2(rng):
    """(ys, ye, xs, xe) for atan2 of two unbounded-exponent values: every exponent difference to ±2300 at any common
    exponent, quotients in and below the subnormal range (x > 0), equal exponents, zeros."""
    diffs = [(_sig(rng), (c := rng.randint(-5000, 5000)) + rng.randint(-2300, 2300), _sig(rng), c)
             for _ in range(120_000)]
    sub = []
    for _ in range(80_000):
        c = rng.randint(-5000, 5000)
        sub.append((_sig(rng), c + rng.randint(-1140, -1010), _sig(rng, False), c))
    exact = []
    for _ in range(20_000):             # short significands: exact quotients, subnormal midpoints among them
        c = rng.randint(-5000, 5000)
        exact.append((_short_sig(rng), c + rng.randint(-1080, -1065), 1.0, c))
    eq = [(_sig(rng), (c := rng.randint(-5000, 5000)), _sig(rng), c) for _ in range(60_000)]
    sp = [(0.0, 0, 1.0, 5000), (-0.0, 0, -1.0, -5000), (1.0, 7000, 0.0, 0), (-1.5, -7000, -0.0, 0), (0.0, 0, -0.0, 0),
          (1.5, -774, 1.0, 300), (1.0, 0, 1.0, -2000), (1.0, -2000, -1.0, 0)]
    return [("special", sp), ("every exponent difference", diffs), ("subnormal quotients", sub),
            ("exact quotients", exact), ("equal exponents", eq)]


def _geom_point(rng, k):
    return [math.ldexp(rng.uniform(-1.0, 1.0), k) for _ in range(3)]


def _binade_any(rng):
    return binade_uniform(rng, 1, -1074, 1020)[0]


def args_geom(rng, n_pts):
    """Edge geometries for potential's distance_3d (2 points), bond_angle (3), dihedral_angle (4), finite coordinates:
    random points at a random scale 2^k (k from -1074 to 1018: Rust's squares, fourth powers and cross products over-
    and underflow), mixed binary magnitudes, subnormal coordinates, points near coincidence and collinearity at random
    scales, and overflowing differences."""
    def scaled(n):
        out = []
        for _ in range(n):
            k = rng.randint(-1074, 1018)
            out.append(tuple(x for _ in range(n_pts) for x in _geom_point(rng, k)))
        return out
    def mixed(n):
        return [tuple(_binade_any(rng) for _ in range(3 * n_pts)) for _ in range(n)]
    def subn(n):
        return [tuple(binade_uniform(rng, 1, -1074, -1023)[0] for _ in range(3 * n_pts)) for _ in range(n)]
    def near(n):
        out = []
        for _ in range(n):
            k = rng.randint(-1074, 1018)
            p = [rng.uniform(-1.0, 1.0) for _ in range(3)]
            pts = [p]
            for j in range(1, n_pts):
                if j == 1:          # B near A (coincidence) or the points near a line (collinearity)
                    off = 2.0 ** rng.uniform(-44, -30)
                    pts.append([p[0] + off * rng.uniform(-1, 1), p[1] + off * rng.uniform(-1, 1), p[2]])
                else:
                    t = rng.uniform(1.0, 2.0)
                    off = 2.0 ** rng.uniform(-40, -26)
                    pts.append([p[0] + t * (pts[1][0] - p[0] + 1.0), p[1] + t * (pts[1][1] - p[1]) + off, p[2]])
            out.append(tuple(math.ldexp(x, k) for q in pts for x in q))
        return out
    big = fb(0x7FEFFFFFFFFFFFFF)
    over = []
    for _ in range(2_000):
        pts = [[rng.choice((big, -big, big / 2.0, -big / 3.0, 0.0, 1.0)) for _ in range(3)] for _ in range(n_pts)]
        over.append(tuple(x for q in pts for x in q))
    return [("scaled", scaled(14_000)), ("mixed magnitudes", mixed(8_000)), ("subnormal coordinates", subn(4_000)),
            ("near coincidence / collinearity at scale", near(10_000)), ("overflowing differences", over)]


def ref_ux(name, a):
    """The MPFR reference for ADR 0002 #20's functions (mpfr_ref.UX at MPFR's widest exponent range, the final
    operation rounded once into binary64); Err (bond and dihedral) encoded as +inf, as the probe does."""
    S = mpfr_ref.ux_scaled
    if name == "uxsqrt":
        return mpfr_ref.ux_sqrt_final(S(a[0], a[1]))
    if name in ("uxmul", "uxdiv", "uxadd"):
        f = {"uxmul": mpfr_ref.ux_mul_final, "uxdiv": mpfr_ref.ux_div_final, "uxadd": mpfr_ref.ux_add_final}[name]
        return f(S(a[0], a[1]), S(a[2], a[3]))
    if name == "uxatan2":
        return mpfr_ref.ux_atan2_final(S(a[0], a[1]), S(a[2], a[3]))
    if name == "dist":
        return mpfr_ref.potential_distance_3d(list(a))
    r = mpfr_ref.potential_bond_angle(list(a), 1e-12) if name == "bond" else \
        mpfr_ref.potential_dihedral_angle(list(a), 1e-20)
    return INF if r is None else r


# ── ADR 0002 #21, #22 (M4e-3): exact rationals rounded once ─────────────────────────────────────────────────────────
import functools  # noqa: E402
from fractions import Fraction  # noqa: E402

_ROOT = Path(__file__).resolve().parents[2]


def _const(path, name):
    """A generated decimal constant (`const NAME = 1.5;`) of a src file, as its double."""
    import re
    m = re.search(rf"^const {name} = ([0-9.]+);$", (_ROOT / path).read_text(), re.M)
    return float(m.group(1))


R_EV = _const("src/constants.cyr", "KM_RYDBERG_EV")
R_INF = _const("src/constants.cyr", "KM_RYDBERG")
R_H = _const("src/spectroscopy_data.cyr", "_KM_RYDBERG_HYDROGEN")
SERIES_NL = (1, 2, 3, 4, 5)        # KmSpectralSeries 0 … 4: Lyman … Pfund (src/spectroscopy_data.cyr)
U32 = (1 << 32) - 1


def rn(q):
    """A Fraction rounded once to binary64 (CPython's int true division is correctly rounded, subnormals and ties
    included); +inf past DBL_MAX."""
    try:
        return float(q)
    except OverflowError:
        return INF


def ref_spec(name, a):
    """The exact reference of the M4e-3 functions (an Err is encoded as +inf, as the probe does)."""
    if name == "ratio":
        p = sum(a[k] << (64 * k) for k in range(4))
        q = sum(a[4 + k] << (64 * k) for k in range(4))
        e = a[8] - (1 << 64) if a[8] >> 63 else a[8]
        return rn(Fraction(p, q) * Fraction(2) ** e)
    if name == "bohr":
        z, n = a
        return -rn(Fraction(R_EV) * z * z / (n * n))
    if name == "trans":
        z, lo, hi = a
        return rn(Fraction(R_EV) * z * z * (hi * hi - lo * lo) / (lo * lo * hi * hi))
    if name in ("ryd", "series"):
        if name == "ryd":
            z, lo, hi, r = a[0], a[1], a[2], R_INF
        else:
            z, lo, r = 1, SERIES_NL[a[0]], R_H
            hi = lo + a[1]
            if hi > U32:
                return INF
        return rn(Fraction(lo * lo * hi * hi) / (Fraction(r) * z * z * (hi * hi - lo * lo)))
    m1, m2 = a
    return rn(Fraction(m1) * Fraction(m2) / (Fraction(m1) + Fraction(m2)))


def _bits_n(rng, lo=1, hi=32):
    """A random integer of a random bit length in lo ... hi (at least 1)."""
    b = rng.randint(lo, hi)
    return rng.randrange(1 << (b - 1), 1 << b)


def args_ratio(rng):
    """(p as 4 words, q as 4 words, e) for RN(p/q·2^e)."""
    def words(v):
        return tuple((v >> (64 * k)) & 0xFFFFFFFFFFFFFFFF for k in range(4))

    def pack(p, q, e):
        return words(p) + words(q) + (e & 0xFFFFFFFFFFFFFFFF,)

    rand = []
    for _ in range(200_000):
        p = _bits_n(rng, 1, 250) if rng.random() < 0.97 else 0
        q = _bits_n(rng, 1, 250)
        e = rng.randint(-1300, 1300)
        rand.append(pack(p, q, e))
    # results at the subnormal grid and near the overflow threshold: e placed so p/q·2^e lands there
    edge = []
    for _ in range(60_000):
        p, q = _bits_n(rng, 1, 250), _bits_n(rng, 1, 250)
        t = rng.choice([rng.randint(-1080, -1018), rng.randint(1018, 1025)])
        e = t - (p.bit_length() - q.bit_length())
        edge.append(pack(p, q, e))
    # exact quotients and midpoints (q | p·2^k), at the 53/54-bit and subnormal boundaries
    ties = []
    for _ in range(60_000):
        q = _bits_n(rng, 1, 120)
        m = rng.randrange(1 << 52, 1 << 54) | rng.choice([0, 1])
        p = m * q
        t = rng.choice([rng.randint(-1130, -1022), rng.randint(-200, 200), rng.randint(970, 1024)])
        ties.append(pack(p, q, t - m.bit_length()))
    return [("random sizes and exponents", rand), ("subnormal and overflow ends", edge),
            ("exact quotients and midpoints", ties)]


def args_bohr(rng):
    """(Z, n) for bohr_energy_level."""
    exact = []
    for z in range(1, 256):
        o = z >> ((z & -z).bit_length() - 1)
        n = o
        while n <= U32:
            exact.append((z, n))
            n *= 2
    small = [(z, n) for z in range(1, 256) for n in range(1, 201)]
    edges = [(rng.randint(1, 255), n) for n in list(range(8180, 8200)) + list(range(U32 - 100, U32 + 1))
             + [(1 << 26) + d for d in range(-3, 4)] for _ in range(20)]
    rand = [(rng.randint(1, 255), _bits_n(rng)) for _ in range(200_000)]
    return [("exact: n = 2^k·Z", exact), ("Z × n <= 200", small), ("n at 2^13, 2^26 and 2^32", edges),
            ("random n of every bit length", rand)]


def _pair(rng):
    lo = _bits_n(rng)
    lo = min(lo, U32 - 1)
    r = rng.random()
    if r < 0.4:
        d = rng.randint(1, 8)
    elif r < 0.7:
        d = _bits_n(rng, 1, 20)
    else:
        d = rng.randint(1, U32 - lo)
    return lo, min(U32, lo + d)


@functools.lru_cache(maxsize=1)
def _exact_transitions():
    """(Z, n1, n2) with n2 <= 120 whose transition energy is exact or a midpoint (by Fraction), and their multiples."""
    out = []
    kr = Fraction(R_EV)
    for z in range(1, 256):
        for lo in range(1, 120):
            for hi in range(lo + 1, 121):
                q = kr * z * z * (hi * hi - lo * lo) / (lo * lo * hi * hi)
                d = q.denominator
                if d & (d - 1) == 0:
                    g = 1
                    while g * hi <= U32 and g <= 1 << 20:
                        out.append((z, g * lo, g * hi))
                        g *= 4
    return out


def args_trans(rng):
    """(Z, n1, n2) for bohr_transition_energy and rydberg_wavelength."""
    ties = [(255, 58055 << k, 524544 << k) for k in range(13)]
    rand = [(rng.randint(1, 255),) + _pair(rng) for _ in range(220_000)]
    adj = []
    for k in range(1, 33):
        for _ in range(400):
            n = min(U32 - 1, rng.randrange(1 << (k - 1), 1 << k))
            adj.append((rng.randint(1, 255), n, n + 1))
    edge = [(rng.randint(1, 255), lo, hi) for hi in range(8186, 8200) for lo in (1, 2, hi // 2, hi - 2, hi - 1)
            for _ in range(10)]
    top = [(rng.randint(1, 255), U32 - k - d, U32 - k) for k in range(200) for d in (1, 2, 3, 1000, 1 << 20)]
    return [("ties", ties), ("exact values (n2 <= 120, and ×4^k)", _exact_transitions()), ("random", rand),
            ("n2 = n1 + 1 at every scale", adj), ("the 2^13 edge", edge), ("next to 2^32 - 1", top)]


def args_series(rng):
    """(series, line) for spectral_series_wavelength."""
    lines = list(range(1, 2001)) + [1 << k for k in range(32)] + [U32 - d for d in range(1, 40)]
    lines += [_bits_n(rng) for _ in range(9000)]
    return [("every series × lines", [(s, ln) for s in range(5) for ln in lines if 1 <= ln <= U32 - SERIES_NL[s]])]


def args_rmass(rng):
    """(m1, m2) for reduced_mass."""
    pos = lambda n, lo=-1074, hi=1023: binade_uniform(rng, n, lo, hi, signs=(1,))
    rand = list(zip(pos(250_000), pos(250_000)))
    eq = [(x, x) for x in pos(80_000)]
    sub = list(zip(pos(40_000, -1074, -1015), pos(40_000, -1074, -1015)))
    sub_eq = [(fb(k), fb(k)) for k in range(1, 2001)] + [(fb(6 * c), fb(2 * c)) for c in range(1, 4001, 2)]
    ties = []
    for _ in range(40_000):
        fam = rng.randint(0, 2)
        if fam == 0:
            c = rng.randrange(1, 94906266, 2)
            a, b = 2.0 ** 27 - c, 2.0 ** 27 + c
        elif fam == 1:
            c = rng.randrange(1, 1 << 25, 2)
            a, b = 3.0 * (3 * 2 ** 25 - c), 3.0 * (3 * 2 ** 25 + c)
        else:
            c = rng.randrange(1, 1 << 26, 2)
            a, b = 3.0 * (2 ** 26 - c), 3.0 * (2 ** 26 + c)
        t = rng.randint(-1049, 994)
        ties.append((a * 2.0 ** t, b * 2.0 ** t) if rng.random() < 0.5 else (b * 2.0 ** t, a * 2.0 ** t))
    gaps = []
    for _ in range(60_000):
        g = rng.randint(45, 64)
        e2 = rng.randint(-1074, 1023 - g)
        m2 = pos(1, e2, e2)[0]
        m1 = pos(1, max(e2, -1022) + g, max(e2, -1022) + g)[0]
        gaps.append((m1, m2) if rng.random() < 0.5 else (m2, m1))
    edge = []
    for m2 in (2.0 ** -1021, 2.0 ** -1022, fb(0x0020000000000000 - 1), fb(0x0020000000000000 + 1)):
        for _ in range(2000):
            edge.append((m2 * rng.uniform(1.0, 2.0 ** rng.randint(0, 60)), m2))
    top = [(fb(0x7FEFFFFFFFFFFFFF - rng.randrange(0, 1 << 52)), fb(0x7FEFFFFFFFFFFFFF - rng.randrange(0, 1 << 54)))
           for _ in range(20_000)]
    return [("binade-uniform pairs", rand), ("equal masses", eq), ("subnormal pairs", sub),
            ("subnormal equal and (6c, 2c) ties", sub_eq), ("normal and subnormal tie families", ties),
            ("exponent gaps 45 ... 64", gaps), ("the 2^-1021 edge", edge), ("the overflow end", top)]


# ── ADR 0002 #24 … #26 (M4f-2): stoichiometry's exact rationals and comparisons ───────────────────────────────────
NANW = 0x7FF8000000000001          # the probe's Err word for the float functions (none of their results is NaN)
ERRW = {1: 0xFFFFFFFFFFFFFF01, 5: 0xFFFFFFFFFFFFFF05, 6: 0xFFFFFFFFFFFFFF06}   # src/error.cyr's codes, encoded
W64 = (1 << 64) - 1
STOICH = [(27, "pyield"), (28, "tyield"), (29, "pcomp"), (30, "limit"), (31, "empir"), (32, "muldiv")]
STOICH_NAMES = tuple(name for _, name in STOICH)
ARITY.update({27: 2, 28: 5, 29: 4, 30: 7, 31: 7, 32: 3})


@functools.lru_cache(maxsize=1)
def masses():
    """src/element_data.cyr's KM_ELEMENT_ATOMIC_MASS (generated from the oracle), by atomic number."""
    import re
    src = (_ROOT / "src/element_data.cyr").read_text()
    block = src[src.index("var KM_ELEMENT_ATOMIC_MASS: f64[118] = {"):]
    vals = re.findall(r"^\s+([0-9.]+),", block[:block.index("};")], re.M)
    assert len(vals) == 118
    return {z + 1: float(v) for z, v in enumerate(vals)}


def _pos(x):
    return math.isfinite(x) and x > 0


def _nonneg(x):
    return math.isfinite(x) and x >= 0


def _rnd(q):
    """Rust's round (halves away from zero) of a rational q >= 0."""
    return math.floor(q + Fraction(1, 2))


def _fw(x):
    return NANW if x is None else bits(x)


def ref_limit(ms, cs):
    """limiting_reagent with exact ratios: the first index of the least m/c; None for Rust's Err."""
    if len(ms) != len(cs) or not ms:
        return None
    best = None
    for i, (m, c) in enumerate(zip(ms, cs)):
        if not _nonneg(m) or not _pos(c):
            return None
        q = Fraction(m) / Fraction(c)
        if best is None or q < best[0]:
            best = (q, i)
    return best[1]


def ref_empir(els):
    """empirical_formula with exact moles, ratios, tolerance 1/10 and counts: ("ok", counts) or ("E", code)."""
    M = masses()
    if not els:
        return ("E", 5)
    moles = []
    for z, p in els:
        if z not in M:
            return ("E", 1)
        if not _pos(p):
            return ("E", 5)
        moles.append(Fraction(p) / Fraction(M[z]))
    mn = min(moles)
    rs = [q / mn for q in moles]
    mult = next((k for k in range(1, 7) if all(abs(r * k - _rnd(r * k)) <= Fraction(1, 10) for r in rs)), None)
    if mult is None:
        return ("E", 6)
    cs = [_rnd(r * mult) for r in rs]
    if any(c > U32 for c in cs):
        return ("E", 6)
    return ("ok", cs)


def ref_stoich(name, a):
    """The exact reference words (w0, w1) of ids 27 … 32 (w1 None: the second slot is the fast path's, compared with the
    first where its flag is 1)."""
    if name == "muldiv":
        x, y, z = a
        return bits(rn(Fraction(x) * Fraction(y) / Fraction(z))), None
    if name == "pyield":
        x, t = a
        if not _pos(t) or not _nonneg(x):
            return NANW, None
        return _fw(x / t * 100.0 if x == 0 else rn(Fraction(x) * 100 / Fraction(t))), None
    if name == "tyield":
        m0, m1, c0, c1, p = a
        if not _pos(p):
            return NANW, None
        i = ref_limit([m0, m1], [c0, c1])
        if i is None:
            return NANW, None
        m, c = (m0, c0) if i == 0 else (m1, c1)
        return _fw(m / c * p if m == 0 else rn(Fraction(m) * Fraction(p) / Fraction(c))), None
    if name == "pcomp":
        z0, n0, z1, n1 = a
        M = masses()
        if z0 not in M or z1 not in M:
            return NANW, None
        mw = 0.0
        for z, c in ((z0, n0), (z1, n1)):
            mw += M[z] * float(c)                  # Rust's molecular_weight: Σ mass·count in f64, in atom order
        if not _pos(mw):
            return NANW, None
        return _fw(0.0 if n0 == 0 else rn(Fraction(M[z0]) * n0 * 100 / Fraction(mw))), None
    if name == "limit":
        n = a[0]
        i = ref_limit(list(a[1:1 + n]), list(a[4:4 + n]))
        return (W64 if i is None else i), None
    n = a[0]
    kind, v = ref_empir([(a[1 + 2 * k], a[2 + 2 * k]) for k in range(n)])
    if kind == "E":
        return ERRW[v], W64
    return v[0] | ((v[1] if n >= 2 else 0) << 32), (v[2] if n == 3 else 0)


def _pp(rng, lo=-1074, hi=1023):
    """A positive double of a random binade in [2^lo, 2^hi]."""
    return binade_uniform(rng, 1, lo, hi, signs=(1,))[0]


def _odd_in(rng, lo, hi):
    """A random odd integer in [lo, hi)."""
    return rng.randrange(lo | 1, hi, 2)


def args_muldiv(rng):
    """(x, y, z) for _km_cr_mul_div = RN(x·y/z)."""
    rand = [(_pp(rng), _pp(rng), _pp(rng)) for _ in range(300_000)]
    mid = [(_pp(rng, -200, 200), _pp(rng, -200, 200), _pp(rng, -200, 200)) for _ in range(100_000)]
    edge = []                                   # k = ex + ey − ez near the fast path's ±1021, and past the range
    while len(edge) < 80_000:
        k = rng.choice([rng.randint(-1085, -1015), rng.randint(1015, 1030)])
        ex, ez = rng.randint(-1074, 1023), rng.randint(-1074, 1023)
        ey = k - ex + ez
        if -1074 <= ey <= 1023:
            edge.append((_pp(rng, ex, ex), _pp(rng, ey, ey), _pp(rng, ez, ez)))
    ties = []                                   # x·y/z = 3y with 3y a 54-bit odd integer (a midpoint), or exact
    for _ in range(80_000):
        z = _odd_in(rng, 1, 1 << 40)
        y = _odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3) if rng.random() < 0.7 else rng.randrange(1, 1 << 53)
        t = rng.randint(-1100, 900)
        sx, sy = rng.randint(-60, 60), t
        x = 3.0 * z * 2.0 ** sx
        yy = float(y) * 2.0 ** sy if sy >= -1074 + 53 else float(Fraction(y) * Fraction(2) ** sy)
        zz = float(z) * 2.0 ** sx
        if _pos(x) and _pos(yy) and _pos(zz):
            ties.append((x, yy, zz))
    near = []                                   # x/z = M/2^53 + c/(z·2^53), M odd with 54 bits: within 2^-100 of a midpoint
    while len(near) < 40_000:
        z = rng.randrange((1 << 52) + 1, 1 << 53, 2)
        c = rng.choice([1, -1, 3, -3, 7, -7])
        x0 = (c * pow(1 << 53, -1, z)) % z      # x0·2^53 ≡ c (mod z)
        x = x0 + z * rng.randint(1, max(1, ((1 << 53) - x0) // z))
        if not z < x < (1 << 53):
            continue
        m = (x * (1 << 53) - c) // z
        if m % 2 == 1 and (1 << 53) <= m < (1 << 54):
            t, j = rng.randint(-900, 900), rng.randint(-60, 60)
            near.append((float(x) * 2.0 ** t, 2.0 ** j, float(z) * 2.0 ** rng.randint(-60, 60)))
    sub = []                                    # subnormal results and their ties: x = q·2^-1074, y/z = 1/2, 3/2, …
    for _ in range(40_000):
        q = rng.randrange(1, 1 << 52)
        y, z = rng.choice([(1.0, 2.0), (3.0, 2.0), (5.0, 2.0), (1.0, 4.0), (3.0, 4.0), (100.0, 200.0), (1.0, 3.0)])
        sub.append((fb(q), y * 2.0 ** rng.randint(0, 40), z * 2.0 ** rng.randint(0, 40)))
    return [("random triples", rand), ("|exponents| <= 200", mid), ("the 2^±1021 edges and past", edge),
            ("exact quotients and midpoints", ties), ("within 2^-100 of a midpoint", near),
            ("subnormal results and ties", sub)]


def args_pyield(rng):
    """(a, t) for percent_yield = RN(100·a/t)."""
    rand = [(_pp(rng), _pp(rng)) for _ in range(200_000)]
    typical = []
    for _ in range(100_000):
        t = _pp(rng, -40, 40)
        typical.append((t * rng.uniform(0.0, 1.2), t))
    ints = [(float(k), 1000.0) for k in range(0, 1001)] + [(float(x), float(t)) for x in range(0, 80) for t in range(1, 80)]
    ties = []
    for _ in range(60_000):
        kind = rng.randint(0, 3)
        e, j = rng.randint(-1000, 900), rng.randint(-60, 60)
        if kind == 0:                           # t = 5·2^j, a = m·2^e with 5m a 54-bit odd integer: a midpoint
            m = _odd_in(rng, (1 << 53) // 5 + 1, (1 << 54) // 5)
            ties.append((m * 2.0 ** e, 5.0 * 2.0 ** j))
        elif kind == 1:                         # t = 3·2^j, a = 3m·2^e with 25m 54 bits and odd
            m = _odd_in(rng, (1 << 53) // 25 + 1, (1 << 54) // 25)
            ties.append((3.0 * m * 2.0 ** e, 3.0 * 2.0 ** j))
        elif kind == 2:                         # a subnormal, t = 200·2^j: (q/2)·2^(-1074 - j), a tie for odd q
            ties.append((fb(rng.randrange(1, 1 << 52)), 200.0 * 2.0 ** rng.randint(0, 3)))
        else:                                   # exact: t = 25·2^j or 100·2^j
            ties.append((_pp(rng), rng.choice([25.0, 100.0]) * 2.0 ** j))
    edge = []
    for _ in range(40_000):
        t = _pp(rng, -60, 60)
        if rng.random() < 0.5:                  # 100·a/t near DBL_MAX
            a = rn(Fraction(fb(0x7FEFFFFFFFFFFFFF)) * Fraction(t) / 100)
        else:                                   # near the subnormal grid's first few units
            a = float(Fraction(rng.randint(1, 1 << 12)) * Fraction(2) ** -1074 * Fraction(t) / 100)
        a = ulps(a, rng.randint(-8, 8))
        if _pos(a):
            edge.append((a, t))
    specials = [(0.0, 10.0), (-0.0, 10.0), (5.0, 0.0), (-1.0, 10.0), (NAN, 10.0), (INF, 10.0), (5.0, INF),
                (5.0, NAN), (fb(1), 10.0), (fb(3), 200.0), (fb(0x7FEFFFFFFFFFFFFF), 1.0), (1.0, fb(1))]
    return [("binade-uniform pairs", rand), ("a = t·U(0, 1.2)", typical), ("integer percentages", ints),
            ("ties and exact values", ties), ("the overflow and subnormal ends", edge), ("specials", specials)]


def args_tyield(rng):
    """(m0, m1, c0, c1, p) for theoretical_yield_moles with two reactants."""
    rand = [(_pp(rng), _pp(rng), _pp(rng), _pp(rng), _pp(rng)) for _ in range(100_000)]
    typical = [(_pp(rng, -10, 10), _pp(rng, -10, 10), float(rng.randint(1, 12)), float(rng.randint(1, 12)),
                float(rng.randint(1, 12))) for _ in range(100_000)]
    ties = []                                   # (a, k·a ± ulps) over (c, k·c): the ratios tie or nearly
    for _ in range(80_000):
        a, c, k = _pp(rng, -1000, 1000), float(rng.randint(1, 50)), float(rng.randint(2, 50))
        ties.append((a, ulps(a * k, rng.randint(-2, 2)), c, c * k, _pp(rng, -100, 100)))
    edge = []                                   # m/c past either end while m·p/c is finite
    for _ in range(60_000):
        if rng.random() < 0.5:
            m, c = _pp(rng, 900, 1023), _pp(rng, -300, -100)
            p = _pp(rng, -1074, -500)
        else:
            m, c = _pp(rng, -1074, -900), _pp(rng, 100, 300)
            p = _pp(rng, 500, 1023)
        edge.append((m, _pp(rng), c, _pp(rng), p))
    specials = [(0.0, 1.0, 1.0, 1.0, 2.0), (-0.0, 0.0, 1.0, 1.0, 2.0), (3.0, 2.0, 2.0, 1.0, 2.0),
                (NAN, 2.0, 2.0, 1.0, 2.0), (3.0, 2.0, 2.0, 1.0, 0.0), (3.0, 2.0, 2.0, 1.0, INF)]
    return [("binade-uniform", rand), ("typical", typical), ("ratio ties and near-ties", ties),
            ("m/c over- and underflowing", edge), ("specials", specials)]


def args_pcomp(rng):
    """(z0, count0, z1, count1) for percent_composition, entry 0."""
    def count():
        r = rng.random()
        if r < 0.05:
            return 0
        if r < 0.1:
            return U32
        if r < 0.6:
            return rng.randint(1, 30)
        return rng.randint(1, U32)
    rand = [(rng.randint(1, 118), count(), rng.randint(1, 118), count()) for _ in range(200_000)]
    same = [(z, count(), z, count()) for z in range(1, 119) for _ in range(200)]
    bad = [(0, 1, 8, 1), (119, 1, 8, 1), (1, 0, 8, 0), (8, 2, 255, 1)]
    return [("random pairs", rand), ("the same element twice", same), ("Err paths", bad)]


def args_limit(rng):
    """(n, m0, m1, m2, c0, c1, c2) for limiting_reagent, n = 1 … 3 (unused slots 0)."""
    def pack(ms, cs):
        n = len(ms)
        return (n,) + tuple(ms) + (0.0,) * (3 - n) + tuple(cs) + (0.0,) * (3 - n)
    rand = []
    for _ in range(150_000):
        n = rng.randint(1, 3)
        rand.append(pack([_pp(rng) for _ in range(n)], [_pp(rng) for _ in range(n)]))
    ties = []
    for _ in range(150_000):
        a, c = _pp(rng, -1000, 1000), float(rng.randint(1, 60))
        k = float(rng.randint(2, 60))
        b = ulps(a * k, rng.randint(-2, 2))
        ms, cs = [a, b], [c, c * k]
        if rng.random() < 0.5:
            ms, cs = ms[::-1], cs[::-1]
        if rng.random() < 0.3:
            ms.append(ulps(a * 3.0, rng.randint(-1, 1)))
            cs.append(c * 3.0)
        ties.append(pack(ms, cs))
    ends = []
    for _ in range(100_000):
        if rng.random() < 0.5:                  # both ratios overflow
            ms, cs = [_pp(rng, 1000, 1023), _pp(rng, 1000, 1023)], [_pp(rng, -1074, -30), _pp(rng, -1074, -30)]
        else:                                   # both underflow
            ms, cs = [_pp(rng, -1074, -1000), _pp(rng, -1074, -1000)], [_pp(rng, 30, 1023), _pp(rng, 30, 1023)]
        if rng.random() < 0.2:
            ms[rng.randint(0, 1)] = 0.0
        ends.append(pack(ms, cs))
    specials = [pack([3.0, 2.0], [2.0, 1.0]), pack([0.0, 0.0], [1.0, 2.0]), pack([-0.0, 0.0], [1.0, 1.0]),
                pack([5.0, NAN], [1.0, 1.0]), pack([5.0, 1.0], [1.0, 0.0]), pack([-1.0, 1.0], [1.0, 1.0]),
                pack([fb(0x7FEFFFFFFFFFFFFF)] * 2, [fb(1), 1e-300]), pack([fb(1), fb(2)], [fb(0x7FEFFFFFFFFFFFFF)] * 2)]
    return [("random lists", rand), ("ties and near-ties", ties), ("both over- or underflowing", ends),
            ("specials", specials)]


def args_empir(rng):
    """(n, z0, p0, z1, p1, z2, p2) for empirical_formula, n = 1 … 3 (unused slots 0)."""
    M = masses()

    def pack(els):
        flat = [v for e in els for v in e]
        return (len(els),) + tuple(flat) + (0, 0.0) * (3 - len(els))
    formulas = []
    for _ in range(60_000):                     # a formula's own percentages (Rust's percent_composition), scaled
        n = rng.randint(1, 3)
        pairs = [(rng.randint(1, 92), rng.randint(1, 12)) for _ in range(n)]
        mw = 0.0
        for z, c in pairs:
            mw += M[z] * float(c)
        s = 10.0 ** rng.uniform(-3, 3)
        formulas.append(pack([(z, M[z] * float(c) / mw * 100.0 * s) for z, c in pairs]))
    edge = []                                   # r·k within a few ulps of n ± 1/10 (the tolerance's knife edge)
    for _ in range(80_000):
        zi, zj = rng.randint(1, 118), rng.randint(1, 118)
        k = rng.randint(1, 6)
        target = Fraction(10 * rng.randint(1, 40) + rng.choice([1, -1]), 10 * k)
        if target < 1:
            continue
        pj = _pp(rng, -20, 20)
        pi = float(target * Fraction(M[zi]) * Fraction(pj) / Fraction(M[zj]))
        edge.append(pack([(zj, pj), (zi, ulps(pi, rng.randint(-40, 40)))]))
    exact = []                                  # the same element twice: r = (10n ± 1)/(10k) exactly, so r·k − n = ±1/10
    for _ in range(40_000):
        z, k = rng.randint(1, 118), rng.randint(1, 6)
        num = 10 * rng.randint(k, 40) + rng.choice([1, -1])
        w = rng.randrange(1, 1 << 40) * 2.0 ** rng.randint(-200, 200)
        els = [(z, float(10 * k) * w), (z, float(num) * w)]
        if rng.random() < 0.5:
            els = els[::-1]
        if rng.random() < 0.3:
            els.append((z, float(10 * k) * w * rng.randint(1, 4)))
        exact.append(pack(els))
    sub = []                                    # subnormal percentages: Rust's moles underflow
    for _ in range(30_000):
        n = rng.randint(1, 3)
        q = _pp(rng, -1074, -1020)
        sub.append(pack([(rng.randint(1, 118), q * rng.choice([1.0, 2.0, 3.0, 0.5])) for _ in range(n)]))
    big = []                                    # ratios near 2^32 − 1/2, 2^32 − 1/4 and past
    for _ in range(20_000):
        zi, zj = rng.randint(1, 118), rng.randint(1, 118)
        target = rng.choice([Fraction(2 ** 32) - Fraction(1, 2), Fraction(2 ** 32) - Fraction(1, 4), Fraction(2 ** 32),
                             Fraction(U32), Fraction(2 ** 33) + 3])
        pj = _pp(rng, -10, 10)
        pi = ulps(float(target * Fraction(M[zi]) * Fraction(pj) / Fraction(M[zj])), rng.randint(-3, 3))
        big.append(pack([(zj, pj), (zi, pi)]))
    specials = [pack([(1, 11.19), (8, 88.81)]), pack([(6, 12.011), (8, 17.5989)]), pack([(6, NAN), (8, 50.0)]),
                pack([(0, 10.0), (8, 50.0)]), pack([(6, 1e-320), (8, 50.0)]), pack([(2, fb(1)), (2, fb(1))]),
                pack([(6, 40.0), (6, 40.0)]), pack([(6, 40.0), (6, 20.0)])]
    return [("formulas' percentages", formulas), ("r·k at n ± 1/10 ± ulps", edge), ("exactly n ± 1/10", exact),
            ("subnormal percentages", sub), ("ratios near 2^32", big), ("specials", specials)]


def stoich_checks(work, targets, which):
    """Ids 27 … 32 against their exact references on every target, the exact path forced; True if anything failed."""
    failed = False
    gens = {"pyield": args_pyield, "tyield": args_tyield, "pcomp": args_pcomp, "limit": args_limit,
            "empir": args_empir, "muldiv": args_muldiv}
    for fid, name in STOICH:
        if name not in which:
            continue
        t0 = time.time()
        rng = random.Random(1027 + fid)
        sections = gens[name](rng)
        xs = [x for _, sec in sections for x in sec]
        want = [ref_stoich(name, a) for a in xs]
        first = None
        for t in targets:
            g0, g1, fl = run_probe(work, t, fid, 0, xs)
            a0, a1, _ = run_probe(work, t, fid, 1, xs)
            two = name == "empir"
            bad = [i for i in range(len(xs)) if g0[i] != want[i][0] or (two and g1[i] != want[i][1])]
            bad_acc = [i for i in range(len(xs)) if a0[i] != want[i][0] or (two and a1[i] != want[i][1])]
            bad_fast = [] if two else [i for i in range(len(xs)) if fl[i] == 1 and g1[i] != g0[i]]
            for i in bad[:5]:
                print(f"  MISMATCH {name}{xs[i]!r} on {t}: got {g0[i]:016x} {g1[i]:016x}, want {want[i]}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {name} (exact path){xs[i]!r} on {t}: got {a0[i]:016x} {a1[i]:016x}, want {want[i]}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {name}{xs[i]!r} on {t}: {g1[i]:016x} vs {g0[i]:016x}")
            cur = (g0, g1, fl, a0, a1)
            xdiff = 0 if first is None else sum(1 for u, v in zip(cur, first) for p, q in zip(u, v) if p != q)
            if first is None:
                first = cur
            slow = len(xs) - sum(1 for f in fl if f == 1)
            print(f"{name:6} {t:8} {len(xs):9,} args: {len(bad)} mismatches; exact path forced: {len(bad_acc)} "
                  f"mismatches; fast path: {len(bad_fast)} differ, fell back {slow:,}; bits or flags differing from "
                  f"{targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0] and not two and name != "limit":
                rates = []
                i = 0
                for tag, sec in sections:
                    rates.append(f"{tag} {sum(1 for f in fl[i:i + len(sec)] if f != 1)}/{len(sec)}")
                    i += len(sec)
                print("         fallbacks by region: " + "; ".join(rates))
        print(f"       [{time.time() - t0:.0f} s]")
    return failed



# ── The overflow audit's B1 (ADR 0002 #27, #28): the correctly rounded monomial and thermo's, phase's, kinetics' ──────
# ids 33 … 44: mono (_km_cr_mono: n_num, n_den, six numerator and four denominator words, the unused 1.0), heat, shc,
# boil, frz, osm, zhl, ssi, pre (heat_transfer, specific_heat_from_calorimetry, boiling_point_elevation,
# freezing_point_depression, osmotic_pressure, zero_order_half_life, steady_state_intermediate, pre_equilibrium_rate:
# the public function, its shape's straight-line fast path and flag; mode 1 the monomial's exact path at its factors),
# soc, mm (second_order_concentration, michaelis_menten: the public function, Rust's chain where the function's fast
# test passes and the flag; mode 1 the unbounded-exponent slow path forced) and monotie (_km_mono_tie: cn, cd, three
# numerator and two denominator significands, z1; the sign decision for z1 < z2 = nextup(z1), against Fraction whatever
# the distance of F from their midpoint). A float function's Err is the NaN word (none of their Ok results is NaN where
# the inputs are finite), a non-finite input of heat / boil / frz Rust's expression (both modes).
AUDIT = [(33, "mono"), (34, "heat"), (35, "shc"), (36, "boil"), (37, "frz"), (38, "osm"), (39, "zhl"), (40, "ssi"),
         (41, "pre"), (42, "soc"), (43, "mm"), (44, "monotie"),
         (45, "sqrtdiv"), (46, "dil"), (47, "cond"), (48, "fmass"), (49, "fcharge"), (50, "fratio"), (51, "nernst"),
         (52, "abs"), (53, "conc"), (54, "vib"), (55, "rotor"), (56, "rot"), (57, "moment")]
AUDIT_NAMES = tuple(name for _, name in AUDIT)
AUDIT2_NAMES = ("sqrtdiv", "dil", "cond", "fmass", "fcharge", "fratio", "nernst", "abs", "conc", "vib", "rotor", "rot",
                "moment")
ARITY.update({33: 12, 34: 3, 35: 3, 36: 3, 37: 3, 38: 3, 39: 2, 40: 3, 41: 4, 42: 3, 43: 3, 44: 8})
ARITY.update({45: 2, 46: 3, 47: 2, 48: 3, 49: 3, 50: 5, 51: 4, 52: 3, 53: 3, 54: 2, 55: 2, 56: 1, 57: 2})
F_CONST = _const("src/constants.cyr", "KM_FARADAY")
H_CONST = _const("src/constants.cyr", "KM_PLANCK")
R_GAS = _const("src/constants.cyr", "KM_GAS_CONSTANT")
DMAX = fb(0x7FEFFFFFFFFFFFFF)


def mono_ref(num, den):
    """RN(Π num / Π den) for finite factors, nonzero denominators: the sign the xor of the factors' signs, a zero factor
    that signed zero, the magnitude rounded once (rn: CPython's correctly rounded int division)."""
    neg = sum(1 for v in num + den if math.copysign(1.0, v) < 0) % 2 == 1
    if any(v == 0.0 for v in num):
        return -0.0 if neg else 0.0
    q = Fraction(1)
    for v in num:
        q *= abs(Fraction(v))
    for v in den:
        q /= abs(Fraction(v))
    r = rn(q)
    return -r if neg else r


def _odd_bits(x):
    """The bit length of the odd part of |x|'s significand (0 for 0)."""
    if x == 0.0:
        return 0
    m = abs(Fraction(x))
    n, d = m.numerator, m.denominator
    while n % 2 == 0:
        n //= 2
    return n.bit_length()


def _ok_side(vals):
    return sum(_odd_bits(v) for v in vals) <= 253


def sqrt_ratio_ref(k, mu):
    """RN(√(k/μ)) for positive finite doubles, exactly (Python's integer square root; B2's vib and sqrtdiv)."""
    q = Fraction(k) / Fraction(mu)
    n, d = q.numerator, q.denominator
    s = max(0, 121 + d.bit_length() - n.bit_length())
    if s % 2:
        s += 1
    quo, rem = divmod(n << s, d)
    r = math.isqrt(quo)
    st = 1 if (r * r != quo or rem) else 0
    return rn(Fraction(2 * r + st, 1 << (s // 2 + 1)))


def _ieee_div(a, b):
    """IEEE division (Python raises on a zero divisor)."""
    if b == 0.0:
        return NAN if (a == 0.0 or math.isnan(a)) else math.copysign(INF, a) * math.copysign(1.0, b)
    return a / b


def _count(n):
    return isinstance(n, int) and 1 <= n <= 255


def ref_audit2(name, a):
    """The reference word of ids 45 … 57 (the overflow audit's B2; NANW for an Err or outside a helper's domain)."""
    if name in ("sqrtdiv", "vib"):
        k, mu = a
        return bits(sqrt_ratio_ref(k, mu)) if _pos(k) and _pos(mu) else NANW
    if name == "dil":
        m1, v1, v2 = a
        if not (_pos(v2) and math.isfinite(m1) and math.isfinite(v1)):
            return NANW
        return bits(mono_ref([m1, v1], [v2]))
    if name == "cond":
        l, c = a
        return bits(mono_ref([l, c], [1000.0]) if math.isfinite(l) and math.isfinite(c) else l * c / 1000.0)
    if name in ("fmass", "fcharge"):
        x, mm, n = a
        if not (_count(n) and _pos(mm) and math.isfinite(x)):
            return NANW
        if name == "fmass":
            return bits(mono_ref([x, mm], [float(n), F_CONST]))
        return bits(mono_ref([x, float(n), F_CONST], [mm]))
    if name == "fratio":
        ma, mma, na, mmb, nb = a
        if not (_count(na) and _count(nb) and _pos(mma) and _pos(mmb) and math.isfinite(ma)):
            return NANW
        return bits(mono_ref([ma, mmb, float(na)], [float(nb), mma]))
    if name == "nernst":
        e0, n, t, q = a
        if not (_pos(t) and _count(n) and _pos(q) and math.isfinite(e0)):
            return NANW
        with mpfr_ref.unbounded():
            c = (mpfr_ref.UX(R_GAS) * mpfr_ref.UX(t)) / (mpfr_ref.UX(float(n)) * mpfr_ref.UX(F_CONST))
            p = c * mpfr_ref.UX(mpfr_ref.log(q))
            return bits(mpfr_ref.ux_sub_final(mpfr_ref.UX(e0), p))
    if name in ("abs", "moment"):
        x, y, z = (a[0], a[1], a[2]) if name == "abs" else (a[0], a[1], a[1])
        return bits(mono_ref([x, y, z], []) if all(math.isfinite(v) for v in (x, y, z)) else x * y * z)
    if name == "conc":
        av, e, l = a
        if not (_pos(e) and _pos(l)):
            return NANW
        return bits(mono_ref([av], [e, l]) if math.isfinite(av) else _ieee_div(av, e * l))
    if name == "rotor":
        j, b = a
        if not (isinstance(j, int) and 0 <= j <= 0xFFFFFFFF):
            return NANW
        jf = float(j)
        return bits(mono_ref([H_CONST, b, jf, jf + 1.0], []) if math.isfinite(b) else H_CONST * b * jf * (jf + 1.0))
    i_ = a[0]
    if not _pos(i_):
        return NANW
    return bits(mono_ref([H_CONST], [8.0, math.pi, math.pi, i_]))


def ref_audit(name, a):
    """The reference word of ids 33 … 44 (NANW for an Err)."""
    if name in AUDIT2_NAMES:
        return ref_audit2(name, a)
    fin = all(math.isfinite(v) for v in a if isinstance(v, float))
    if name == "mono":
        nn, nd = a[0], a[1]
        num, den = list(a[2:2 + nn]), list(a[8:8 + nd])
        if not all(math.isfinite(v) for v in num + den) or any(v == 0.0 for v in den):
            return NANW
        return bits(mono_ref(num, den))
    if name in ("heat", "boil", "frz"):
        x, y, z = a
        return bits(mono_ref([x, y, z], []) if fin else x * y * z)
    if name == "shc":
        e, m, dt = a
        if not (math.isfinite(e) and _pos(m) and math.isfinite(dt) and abs(dt) >= 1e-12):
            return NANW
        return bits(mono_ref([e], [m, dt]))
    if name == "osm":
        i, c, t = a
        if not (_pos(t) and _pos(i) and _nonneg(c)):
            return NANW
        return bits(mono_ref([i, c, R_GAS, 1000.0, t], []))
    if name == "zhl":
        c0, k = a
        if not (_pos(k) and _pos(c0)):
            return NANW
        return bits(mono_ref([c0], [2.0, k]))
    if name == "ssi":
        k1, c, k2 = a
        if not (_pos(k2) and _nonneg(k1) and _nonneg(c)):
            return NANW
        return bits(mono_ref([k1, c], [k2]))
    if name == "pre":
        kf, kr, ks, c = a
        if not (_pos(kr) and _nonneg(kf) and _nonneg(ks) and _nonneg(c)):
            return NANW
        return bits(mono_ref([ks, kf, c], [kr]))
    if name == "soc":
        c0, k, t = a
        if not (_pos(c0) and _nonneg(k) and _nonneg(t)):
            return NANW
        with mpfr_ref.unbounded():
            one = mpfr_ref.UX(1.0)
            return bits(mpfr_ref.ux_div_final(one, one / mpfr_ref.UX(c0) + mpfr_ref.UX(k) * mpfr_ref.UX(t)))
    if name == "mm":
        v, km, sb = a
        if not (_pos(v) and _pos(km) and _nonneg(sb)):
            return NANW
        with mpfr_ref.unbounded():
            p = mpfr_ref.UX(v) * mpfr_ref.UX(sb)
            return bits(mpfr_ref.ux_div_final(p, mpfr_ref.UX(km) + mpfr_ref.UX(sb)))
    cn, cd = a[0], a[1]
    f = Fraction(1)
    for v in a[2:2 + cn]:
        f *= Fraction(v)
    for v in a[5:5 + cd]:
        f /= Fraction(v)
    z1 = a[7]
    z2 = math.nextafter(z1, INF)
    mu = (Fraction(z1) + Fraction(z2)) / 2
    if f > mu:
        return bits(z2)
    if f < mu:
        return bits(z1)
    return bits(z1) if bits(z1) % 2 == 0 else bits(z2)


def _unit_sig(rng):
    """A random double in [1, 2) (B1's monotie)."""
    return 1.0 + rng.randrange(1 << 52) * 2.0 ** -52


def _spp(rng, lo=-1074, hi=1023):
    """A signed double of a random binade in [2^lo, 2^hi]."""
    return _pp(rng, lo, hi) * rng.choice((1.0, -1.0))


def _around(rng, c, w=200):
    """A positive double of a random binade within w of 2^c, clamped to the range (B1's compensating factors)."""
    lo = min(max(-1074, c - w), 1023)
    return _pp(rng, lo, max(lo, min(1023, c + w)))


SHAPES = [(1, 0), (2, 0), (3, 0), (4, 0), (1, 1), (2, 1), (3, 1), (4, 1), (1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (1, 4),
          (4, 2), (2, 4)]


def _mono_pack(num, den):
    return (len(num), len(den)) + tuple(num) + (1.0,) * (6 - len(num)) + tuple(den) + (1.0,) * (4 - len(den))


def args_mono(rng):
    """Arguments of _km_cr_mono: every shape up to four factors a side (five and six with short factors: osmotic's),
    within the exact path's capacity."""
    def factor(lo=-1074, hi=1023):
        r = rng.random()
        if r < 0.03:
            return 2.0 ** rng.randint(max(lo, -1074), min(hi, 1023)) * rng.choice((1.0, -1.0))
        if r < 0.08:
            return float(rng.randint(1, 1000)) * rng.choice((1.0, -1.0)) * 2.0 ** rng.randint(-40, 40)
        return _spp(rng, lo, hi)
    rand, mid, edge, ties, zeros, five, dom = [], [], [], [], [], [], []
    for _ in range(140_000):
        nn, nd = rng.choice(SHAPES)
        rand.append(_mono_pack([factor() for _ in range(nn)], [factor() for _ in range(nd)]))
    for _ in range(100_000):
        nn, nd = rng.choice(SHAPES)
        mid.append(_mono_pack([factor(-200, 200) for _ in range(nn)], [factor(-200, 200) for _ in range(nd)]))
    while len(edge) < 60_000:                       # the total exponent near the subnormal range, 2^-1075 or 2^1024
        nn, nd = rng.choice(SHAPES)
        num, den = [factor(-300, 300) for _ in range(nn)], [factor(-300, 300) for _ in range(nd)]
        q = Fraction(1)
        for v in num:
            q *= abs(Fraction(v))
        for v in den:
            q /= abs(Fraction(v))
        e = q.numerator.bit_length() - q.denominator.bit_length()
        sh = rng.choice([rng.randint(-1080, -1018), rng.randint(1018, 1026)]) - e
        lead = num[0]
        try:
            moved = math.ldexp(lead, sh)
        except OverflowError:
            continue
        if math.isfinite(moved) and moved != 0.0 and math.ldexp(moved, -sh) == lead:
            num[0] = moved
            edge.append(_mono_pack(num, den))
    for _ in range(80_000):                         # exact 53-bit midpoints, and the subnormal grid's ties
        y = float(_odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3))
        k, j = rng.randint(-900, 900), rng.randint(-60, 60)
        z, w = float(rng.randrange(1, 1 << 20) | 1), float(rng.randrange(1, 1 << 20) | 1)
        kind = rng.randint(0, 6)
        if kind == 0:
            num, den = [3.0, y * 2.0 ** k], []
        elif kind == 1:
            num, den = [3.0, y * 2.0 ** k, 2.0 ** j], []
        elif kind == 2:
            num, den = [3.0 * z, y * 2.0 ** k], [z]
        elif kind == 3:
            num, den = [3.0, y * 2.0 ** k, w], [w]
        elif kind == 4:
            num, den = [3.0 * z * w, y * 2.0 ** k], [z, w]
        elif kind == 5:
            num, den = [-1.5, y, 2.0 ** j, w * 2.0 ** -20], [w]
        else:
            q = rng.randrange(1, 1 << 52) | 1
            num, den = [fb(q), 0.5 * 2.0 ** rng.randint(0, 40)], [2.0 ** rng.randint(0, 40)]
        ties.append(_mono_pack(num, den))
    for _ in range(10_000):                         # zeros and signs, overflowing partial products included
        nn, nd = rng.choice(SHAPES)
        num = [factor() for _ in range(nn)]
        num[rng.randrange(nn)] = rng.choice((0.0, -0.0))
        zeros.append(_mono_pack(num, [factor() for _ in range(nd)]))
    while len(five) < 40_000:                       # five and six numerator factors (osmotic's i·M·R·1000·T)
        nn = rng.choice((5, 6))
        num = [factor() for _ in range(3)] + [R_GAS, 1000.0] + ([float(rng.randint(1, 12))] if nn == 6 else [])
        rng.shuffle(num)
        den = [factor() for _ in range(rng.randint(0, 2))]
        if _ok_side(num) and _ok_side(den):
            five.append(_mono_pack(num, den))
    dom = [_mono_pack([1.0, INF], []), _mono_pack([NAN], [1.0]), _mono_pack([1.0], [0.0]), _mono_pack([0.0], [-0.0]),
           _mono_pack([DMAX, DMAX], [fb(1)]), _mono_pack([fb(1), fb(1)], [DMAX]), _mono_pack([], []),
           _mono_pack([fb(1)], [2.0]), _mono_pack([fb(3)], [2.0]), _mono_pack([DMAX], [0.5])]
    return [("random shapes", rand), ("|exponents| <= 200", mid), ("the result at the range ends", edge),
            ("midpoints and subnormal ties", ties), ("zeros", zeros), ("five and six factors", five), ("domain", dom)]


def args_mul3(rng):
    """(x, y, z) for heat_transfer and the colligative functions."""
    rand = [(_spp(rng), _spp(rng), _spp(rng)) for _ in range(80_000)]
    mid = [(_spp(rng, -150, 150), _spp(rng, -150, 150), _spp(rng, -150, 150)) for _ in range(60_000)]
    typical = [(_pp(rng, -10, 20), _pp(rng, -4, 5), _spp(rng, -10, 10)) for _ in range(60_000)]
    comp = []                                       # x·y past DBL_MAX or below 2^-1022, z bringing it back
    for _ in range(40_000):
        if rng.random() < 0.5:
            x, y = _spp(rng, 520, 1023), _spp(rng, 520, 1023)
        else:
            x, y = _spp(rng, -1074, -520), _spp(rng, -1074, -520)
        e = math.frexp(x)[1] + math.frexp(y)[1]
        comp.append((x, y, _around(rng, -e) * rng.choice((1.0, -1.0))))
    ties = []
    for _ in range(40_000):
        y = float(_odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3))
        if rng.random() < 0.6:
            ties.append((3.0 * rng.choice((1.0, -1.0)), y * 2.0 ** rng.randint(-900, 900), 2.0 ** rng.randint(-60, 60)))
        else:
            ties.append((fb(rng.randrange(1, 1 << 52) | 1), 0.5 * 2.0 ** rng.randint(0, 40), 2.0 ** -rng.randint(0, 40)))
    sp = [0.0, -0.0, INF, -INF, NAN, 1.0, fb(1), DMAX, 1e300, 1e-300]
    specials = [(a, b, c) for a in sp for b in sp for c in sp]
    return [("binade-uniform", rand), ("|exponents| <= 150", mid), ("typical", typical),
            ("partial products past the range", comp), ("midpoints and subnormal ties", ties), ("specials", specials)]


def args_shc(rng):
    """(q, m, ΔT) for specific_heat_from_calorimetry."""
    rand = [(_spp(rng), _pp(rng), _spp(rng, -39, 1023)) for _ in range(80_000)]
    typical = [(_spp(rng, -10, 30), _pp(rng, -10, 20), _spp(rng, -10, 10)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            m, dt = _pp(rng, 520, 1023), _spp(rng, 520, 1023)
            q = _spp(rng, 600, 1023)
        else:
            m, dt = _pp(rng, -1074, -600), _spp(rng, -39, 20)
            q = _spp(rng, -1074, -900)
        comp.append((q, m, dt))
    sub = [(fb(rng.randrange(1, 1 << 52)) * rng.choice((1.0, -1.0)), float(rng.choice((1, 2, 4, 3, 5))),
            rng.choice((1.0, -1.0, 2.0, 0.5))) for _ in range(20_000)]
    eps = 1e-12
    bad = [(1.0, 1.0, math.nextafter(eps, 0.0)), (1.0, 1.0, eps), (1.0, 0.0, 1.0), (NAN, 1.0, 1.0), (1.0, 1.0, INF),
           (INF, 1.0, 1.0), (1.0, -1.0, 1.0), (0.0, 1.0, 1.0), (-0.0, 1.0, -1.0), (0.0, DMAX, DMAX)]
    return [("binade-uniform", rand), ("typical", typical), ("m·ΔT past the range", comp),
            ("subnormal results and ties", sub), ("Err and zeros", bad)]


def args_osm(rng):
    """(i, M, T) for osmotic_pressure."""
    rand = [(_pp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(1.0, 4.0), _pp(rng, -20, 4), _pp(rng, 0, 13)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            i, c = _pp(rng, 520, 1023), _pp(rng, 520, 1023)
        else:
            i, c = _pp(rng, -1074, -520), _pp(rng, -1074, -520)
        e = math.frexp(i)[1] + math.frexp(c)[1] + 13
        comp.append((i, c, _around(rng, -e)))
    sub = [(_pp(rng, -1074, -1000), _pp(rng, -40, 2), _pp(rng, -20, 13)) for _ in range(20_000)]
    bad = [(2.0, 0.0, 298.15), (2.0, -0.0, 298.15), (0.0, 0.1, 298.15), (2.0, 0.1, 0.0), (2.0, -1.0, 298.15),
           (NAN, 0.1, 298.15), (2.0, INF, 298.15), (2.0, 0.1, INF)]
    return [("binade-uniform", rand), ("typical", typical), ("partial products past the range", comp),
            ("subnormal results", sub), ("Err and zeros", bad)]


def args_zhl(rng):
    """([A]0, k) for zero_order_half_life."""
    rand = [(_pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 5), _pp(rng, -30, 10)) for _ in range(40_000)]
    big = [(_pp(rng), 2.0 ** 1023 * rng.uniform(1.0, 2.0)) for _ in range(40_000)] + \
        [(_pp(rng), ulps(2.0 ** 1023, rng.randint(-5, 5))) for _ in range(5_000)]
    sub = [(_pp(rng, -1074, -900), _pp(rng, 0, 200)) for _ in range(20_000)]
    bad = [(1.0, 0.0), (0.0, 1.0), (NAN, 1.0), (1.0, INF), (-1.0, 1.0), (1.0, -1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("k >= 2^1023: 2k overflows", big),
            ("subnormal results", sub), ("Err", bad)]


def args_ssi(rng):
    """(k1, [A], k2) for steady_state_intermediate."""
    rand = [(_pp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 10), _pp(rng, -20, 4), _pp(rng, -10, 20)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            a, b = _pp(rng, 520, 1023), _pp(rng, 520, 1023)
        else:
            a, b = _pp(rng, -1074, -520), _pp(rng, -1074, -520)
        e = math.frexp(a)[1] + math.frexp(b)[1]
        comp.append((a, b, _around(rng, e)))
    ties = []
    for _ in range(40_000):
        y = float(_odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3))
        if rng.random() < 0.6:
            ties.append((y * 2.0 ** rng.randint(-900, 900), 3.0, 2.0 ** rng.randint(-60, 60)))
        else:
            ties.append((fb(rng.randrange(1, 1 << 52) | 1), 1.0, 2.0 * 2.0 ** rng.randint(0, 40)))
    bad = [(0.0, 1.0, 1.0), (-0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0), (NAN, 1.0, 1.0), (1.0, -1.0, 1.0),
           (0.0, DMAX, fb(1))]
    return [("binade-uniform", rand), ("typical", typical), ("k1·[A] past the range", comp),
            ("midpoints and subnormal ties", ties), ("Err and zeros", bad)]


def args_pre(rng):
    """(k_forward, k_reverse, k_slow, [A]) for pre_equilibrium_rate."""
    rand = [(_pp(rng), _pp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_pp(rng, -10, 20), _pp(rng, -10, 20), _pp(rng, -20, 10), _pp(rng, -20, 4)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            kf, kr = _pp(rng, 400, 1023), _pp(rng, -1074, -400)
        else:
            kf, kr = _pp(rng, -1074, -400), _pp(rng, 400, 1023)
        e = math.frexp(kf)[1] - math.frexp(kr)[1]
        ks = _around(rng, -e, 300)
        comp.append((kf, kr, ks, _pp(rng, -100, 100)))
    ties = []
    for _ in range(40_000):
        y = float(_odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3))
        ties.append((3.0, 2.0 ** rng.randint(-60, 60), y * 2.0 ** rng.randint(-800, 800), 2.0 ** rng.randint(-60, 60)))
    bad = [(1e300, 1.0, 1e300, 0.0), (1e300, 1e-300, 0.0, 1.0), (0.0, 1.0, 1.0, 1.0), (1.0, 0.0, 1.0, 1.0),
           (NAN, 1.0, 1.0, 1.0), (1.0, 1.0, -1.0, 1.0), (1.0, 1.0, 1.0, -1.0), (-0.0, 1.0, 1.0, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("k_f/k_r past the range", comp), ("midpoints", ties),
            ("Err and zeros", bad)]


def args_soc(rng):
    """([A]0, k, t) for second_order_concentration."""
    def nz(lo=-1074, hi=1023):
        return 0.0 if rng.random() < 0.03 else _pp(rng, lo, hi)
    rand = [(_pp(rng), nz(), nz()) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 4), nz(-20, 20), nz(-10, 20)) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        if kind == 0:
            edge.append((_pp(rng, -1074, -1015), nz(-40, 40), nz(-40, 40)))
        elif kind == 1:
            edge.append((_pp(rng, 1015, 1023), nz(-1074, -900), nz(-40, 40)))
        elif kind == 2:
            edge.append((_pp(rng, -40, 40), _pp(rng, 500, 1023), _pp(rng, 500, 1023)))
        elif kind == 3:
            edge.append((_pp(rng, -40, 40), _pp(rng, -1074, -500), _pp(rng, -1074, -500)))
        else:
            edge.append((_pp(rng, -1074, -900), _pp(rng, 800, 1023), _pp(rng, 0, 200)))
    bad = [(0.0, 1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, -1.0), (NAN, 1.0, 1.0), (1.0, INF, 1.0), (1.0, 1.0, NAN),
           (DMAX, 0.0, 0.0), (fb(1), 0.0, 1.0), (1.0, 2.0 ** 1000, 2.0 ** 75), (1.0, -0.0, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range", edge), ("Err and specials", bad)]


def args_mm(rng):
    """(Vmax, Km, [S]) for michaelis_menten."""
    def nz(lo=-1074, hi=1023):
        return 0.0 if rng.random() < 0.03 else _pp(rng, lo, hi)
    rand = [(_pp(rng), _pp(rng), nz()) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 10), _pp(rng, -20, 4), nz(-27, 10)) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 3)
        if kind == 0:
            edge.append((_pp(rng, 0, 1023), _pp(rng, 1000, 1023), _pp(rng, 1000, 1023)))
        elif kind == 1:
            edge.append((_pp(rng, 500, 1023), _pp(rng, -1074, 1023), _pp(rng, 500, 1023)))
        elif kind == 2:
            edge.append((_pp(rng, -1074, -500), _pp(rng, -1074, 0), _pp(rng, -1074, -500)))
        else:
            edge.append((_pp(rng, -1074, 1023), _pp(rng, -1074, -1015), _pp(rng, -1074, -1015)))
    bad = [(0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, -1.0), (NAN, 1.0, 1.0), (1.0, 1.0, INF), (100.0, DMAX, DMAX),
           (DMAX, DMAX, DMAX), (1.0, 1.0, -0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range", edge), ("Err and specials", bad)]


def args_monotie(rng):
    """(cn, cd, s0, s1, s2, t0, t1, z1) for _km_mono_tie: z1 the double below F or a nearby one, and exact midpoints."""
    near, far, mids = [], [], []
    for _ in range(100_000):
        cn, cd = rng.randint(0, 3), rng.randint(0, 2)
        ss, ts = [_unit_sig(rng) for _ in range(cn)], [_unit_sig(rng) for _ in range(cd)]
        f = Fraction(1)
        for v in ss:
            f *= Fraction(v)
        for v in ts:
            f /= Fraction(v)
        z = rn(f)
        z1 = z if Fraction(z) <= f else math.nextafter(z, 0.0)
        near.append((cn, cd) + tuple(ss + [1.0] * (3 - cn)) + tuple(ts + [1.0] * (2 - cd)) + (z1,))
        far.append((cn, cd) + tuple(ss + [1.0] * (3 - cn)) + tuple(ts + [1.0] * (2 - cd)) +
                   (ulps(z1, rng.choice((-3, -2, -1, 1, 2, 3))),))
    for _ in range(40_000):                         # F = 1.5·s_y = 3y/2^b exactly a midpoint (3y odd, of 54 bits)
        y = _odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3)
        b = y.bit_length()
        sy = y * 2.0 ** -(b - 1)
        f = Fraction(3 * y, 1 << b)
        zr = rn(f)
        z1 = zr if Fraction(zr) < f else math.nextafter(zr, 0.0)
        mids.append((2, 0, 1.5, sy, 1.0, 1.0, 1.0, z1))
    return [("z1 below F", near), ("z1 a few ulps away", far), ("exact midpoints", mids)]



# ── The overflow audit's B2 (ADR 0002 #29, #30): the new shapes, √ of a quotient, and the twelve functions ─────────────
# ids 45 … 57: sqrtdiv (_km_cr_sqrt_div(x, y): the public helper, its fast path and flag; mode 1 _km_sqrt_div_acc), dil,
# cond (dilution, solution_conductivity), fmass, fcharge, fratio (Faraday's three, an electron count as its integer
# word), abs, conc, rotor (J as its integer word), rot, moment (absorbance, concentration_from_absorbance,
# rigid_rotor_energy, rotational_constant, diatomic_moment_of_inertia): the public function, its shape's straight-line fast
# path (_km_mul_div_fast, _km_mul2_div2_fast, _km_mul3_div_fast, _km_mul3_div2_fast, _km_mul3_fast, _km_div2_fast,
# _km_mul4_fast, _km_div3_fast) and flag, mode 1 the monomial's exact path at its factors; vib (vibrational_frequency:
# _km_sqrt_div_fast, mode 1 _km_sqrt_div_acc); nernst (nernst_potential: the public function, _km_nernst_fast where its
# range test passes and the flag, mode 1 _km_nernst_ux). References: Fraction (the monomials), Python's integer square
# root (√), mpfr_ref.UX with MPFR's log and the final subtraction (Nernst). An Err row, a NaN result of a helper outside
# its domain, and a non-finite input of an unvalidated function (Rust's expression) give the public word in both slots.
def _y3(rng):
    """An odd y with 3y odd and of 54 bits (a 53-bit midpoint), as a double."""
    return float(_odd_in(rng, (1 << 53) // 3 + 1, (1 << 54) // 3))


def _sub_tie(rng):
    """(q·2^a, 2^(a + 1)) for an odd q: their quotient a tie of the subnormal grid."""
    a = rng.randint(0, 50)
    return fb(rng.randrange(1, 1 << 41) | 1) * 2.0 ** a, 2.0 ** (a + 1)


def _cf_best(x, lim):
    """The last continued-fraction convergent p/q of the positive Fraction x with p, q below lim."""
    h0, h1, k0, k1 = 0, 1, 1, 0
    xx = x
    while True:
        a = xx.numerator // xx.denominator
        h2, k2 = a * h1 + h0, a * k1 + k0
        if h2 >= lim or k2 >= lim:
            return h1, k1
        h0, h1, k0, k1 = h1, h2, k1, k2
        f = xx - a
        if f == 0:
            return h1, k1
        xx = 1 / f


def args_sqrtdiv(rng, n_rand=120_000):
    """(k, μ) for _km_cr_sqrt_div and vibrational_frequency: every binade of both, k/μ past DBL_MAX and below 2^-1074
    (Rust's ∞ and 0), subnormal operands and results, the overflow end, exact roots s²·2^2j/μ, near-ties (the root
    within a few ulps of a 53-bit midpoint, from the midpoint's square), and the domain."""
    rand = [(_pp(rng), _pp(rng)) for _ in range(n_rand)]
    typical = [(_pp(rng, -10, 12), _pp(rng, -95, -80)) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        if kind == 0:
            edge.append((_pp(rng, 512, 1023), _pp(rng, -1022, -512)))          # k/μ past DBL_MAX
        elif kind == 1:
            edge.append((_pp(rng, -1022, -540), _pp(rng, 538, 1023)))          # k/μ below 2^-1074
        elif kind == 2:
            edge.append((_pp(rng, -1074, -1023), _pp(rng)))                   # a subnormal k
        elif kind == 3:
            edge.append((_pp(rng), _pp(rng, -1074, -1023)))                   # a subnormal μ
        else:
            edge.append((_pp(rng, -1074, -980), _pp(rng, 990, 1023)))          # a subnormal result
    exact = []
    for _ in range(30_000):
        sq = float(rng.randrange(1, 1 << 26))
        mu = _pp(rng, -500, 500)
        j = rng.randint(-200, 200)
        k = sq * sq * 2.0 ** (2 * j) * mu
        if _pos(k) and k / mu == sq * sq * 2.0 ** (2 * j):
            exact.append((k, mu))
    near, close = [], []
    for _ in range(40_000):                         # m = a 53-bit midpoint; k = RN(m²·μ) puts √(k/μ) within ulps of m
        z = 1.0 + rng.randrange(1 << 52) * 2.0 ** -52
        m = Fraction(z) + Fraction(1, 1 << 53)
        mu = _pp(rng, -300, 300)
        k = rn(m * m * Fraction(mu))
        k = ulps(k, rng.randint(-3, 3))
        if _pos(k):
            near.append((k, mu))
    for _ in range(20_000):                         # k/μ = p/q, m²'s best approximation with p, q below 2^53: within
        z = 1.0 + rng.randrange(1 << 52) * 2.0 ** -52   # about 2^-105 of m², inside the fast path's 2^-98 band
        m = Fraction(z) + Fraction(1, 1 << 53)
        p, q = _cf_best(m * m, 1 << 53)
        e = rng.randint(-480, 480)
        close.append((float(p) * 2.0 ** e, float(q) * 2.0 ** (e - 2 * rng.randint(-20, 20))))
    dom = [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (1.0, -1.0), (INF, 1.0), (1.0, INF), (NAN, 1.0), (DMAX, fb(1)),
           (fb(1), DMAX), (fb(1), fb(1)), (4.0, 1.0), (2.0, 1.0)]
    return [("binade-uniform", rand), ("typical (HCl's k, μ)", typical), ("k/μ past either end, subnormals", edge),
            ("exact roots", exact), ("within ulps of a 53-bit midpoint", near),
            ("within 2^-100 of a 53-bit midpoint", close), ("domain", dom)]


def args_dil(rng):
    """(M1, V1, V2) for dilution."""
    rand = [(_spp(rng), _spp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 4), _pp(rng, -10, 10), _pp(rng, -10, 10)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            a, b = _spp(rng, 520, 1023), _spp(rng, 520, 1023)
        else:
            a, b = _spp(rng, -1074, -520), _spp(rng, -1074, -520)
        e = math.frexp(a)[1] + math.frexp(b)[1]
        comp.append((a, b, _around(rng, e)))
    ties = []
    for _ in range(40_000):
        if rng.random() < 0.6:
            ties.append((_y3(rng) * 2.0 ** rng.randint(-900, 900), 3.0, 2.0 ** rng.randint(-60, 60)))
        else:
            q, d = _sub_tie(rng)
            ties.append((q, 1.0, d))
    bad = [(1.0, 1.0, 0.0), (1.0, 1.0, -1.0), (NAN, 1.0, 1.0), (1.0, INF, 1.0), (1.0, 1.0, INF), (0.0, DMAX, fb(1)),
           (-0.0, 1e300, 1e-300), (DMAX, DMAX, DMAX)]
    return [("binade-uniform", rand), ("typical", typical), ("M1·V1 past the range", comp),
            ("midpoints and subnormal ties", ties), ("Err and zeros", bad)]


def args_cond(rng):
    """(Λ, c) for solution_conductivity."""
    rand = [(_spp(rng), _spp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(50.0, 450.0), _pp(rng, -20, 3)) for _ in range(60_000)]
    edge = []
    for _ in range(40_000):
        a = _spp(rng, 1, 1023)
        e = math.frexp(a)[1]
        edge.append((a, _around(rng, 1025 - e, 9)))                           # Λ·c past DBL_MAX, Λ·c/1000 not
        edge.append((_spp(rng, -1074, -900), _pp(rng, -200, 60)))             # subnormal results
    ties = [(3000.0 * rng.choice((1.0, -1.0)), _y3(rng) * 2.0 ** rng.randint(-900, 900)) for _ in range(30_000)]
    sp = [0.0, -0.0, INF, -INF, NAN, 1.0, fb(1), DMAX, 1e300, 1e-300, 126.4]
    specials = [(a, b) for a in sp for b in sp]
    return [("binade-uniform", rand), ("typical", typical), ("Λ·c past the range, subnormal results", edge),
            ("midpoints 3000·y/1000", ties), ("specials", specials)]


def _n(rng):
    return rng.choice((1, 2, 3, 4, 5, 6, 127, 128, 254, 255, rng.randint(1, 255)))


def args_fmass(rng):
    """(Q, M, n) for faraday_mass_deposited."""
    rand = [(_spp(rng), _pp(rng), _n(rng)) for _ in range(80_000)]
    typical = [(_spp(rng, -10, 27), rng.uniform(1.0, 250.0), _n(rng)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            q, mm = _spp(rng, 520, 1023), _pp(rng, 520, 1023)
        else:
            q, mm = _spp(rng, -1074, -520), _pp(rng, -1074, -400)
        comp.append((q, mm, _n(rng)))
    f3 = 3.0 * F_CONST
    ties = [(f3 * 2.0 ** rng.randint(-600, 600) * rng.choice((1.0, -1.0)), _y3(rng), 1) for _ in range(30_000)]
    bad = [(1.0, 1.0, 0), (1.0, 1.0, 256), (1.0, 1.0, -1), (1.0, 0.0, 2), (1.0, -1.0, 2), (NAN, 1.0, 2), (INF, 1.0, 2),
           (0.0, DMAX, 1), (-0.0, DMAX, 255), (DMAX, DMAX, 1)]
    return [("binade-uniform", rand), ("typical", typical), ("Q·M past the range", comp), ("midpoints 3F·y/F", ties),
            ("Err and zeros", bad)]


def args_fcharge(rng):
    """(m, M, n) for faraday_charge_required."""
    rand = [(_spp(rng), _pp(rng), _n(rng)) for _ in range(80_000)]
    typical = [(_spp(rng, -20, 13), rng.uniform(1.0, 250.0), _n(rng)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            ma, mm = _spp(rng, 990, 1023), _pp(rng, 0, 1023)
        else:
            ma, mm = _spp(rng, -1074, -960), _pp(rng, -1022, 0)
        comp.append((ma, mm, _n(rng)))
    ties = [(_y3(rng) * 2.0 ** rng.randint(-600, 600), F_CONST * 2.0 ** rng.randint(-60, 60), 3) for _ in range(30_000)]
    bad = [(1.0, 1.0, 0), (1.0, 1.0, 300), (1.0, 0.0, 2), (NAN, 1.0, 2), (1.0, INF, 2), (-0.0, DMAX, 7), (0.0, fb(1), 1)]
    return [("binade-uniform", rand), ("typical", typical), ("m·n·F past the range", comp), ("midpoints y·3·F/F", ties),
            ("Err and zeros", bad)]


def args_fratio(rng):
    """(m_A, M_A, n_A, M_B, n_B) for faraday_mass_ratio."""
    rand = [(_spp(rng), _pp(rng), _n(rng), _pp(rng), _n(rng)) for _ in range(80_000)]
    typical = [(_pp(rng, -20, 13), rng.uniform(1.0, 250.0), _n(rng), rng.uniform(1.0, 250.0), _n(rng))
               for _ in range(60_000)]
    edge = []
    for _ in range(40_000):
        if rng.random() < 0.5:                      # M_A/n_A subnormal or 0
            edge.append((_spp(rng, -1000, 1000), _pp(rng, -1074, -1010), _n(rng), _pp(rng, -1000, 0), _n(rng)))
        else:                                       # m_A·(M_B/n_B) past DBL_MAX
            edge.append((_spp(rng, 520, 1023), _pp(rng, 520, 1023), _n(rng), _pp(rng, 520, 1023), _n(rng)))
    ties = []
    for _ in range(30_000):
        nb = _n(rng)
        w = float(rng.randrange(1, 1 << 20) | 1)
        ties.append((_y3(rng) * 2.0 ** rng.randint(-600, 600), w, 3, w * nb, nb))
    bad = [(1.0, 1.0, 0, 1.0, 1), (1.0, 1.0, 1, 1.0, 0), (1.0, 0.0, 1, 1.0, 1), (1.0, 1.0, 1, -1.0, 1),
           (NAN, 1.0, 1, 1.0, 1), (1.0, fb(1), 3, 1e-300, 1), (-0.0, fb(1), 255, DMAX, 1), (DMAX, DMAX, 1, 2.0, 1)]
    return [("binade-uniform", rand), ("typical", typical), ("M_A/n_A or m_A·M_B/n_B past the range", edge),
            ("midpoints y·w·3/(n_B·w)·n_B", ties), ("Err and zeros", bad)]


def args_nernst(rng):
    """(E°, n, T, Q) for nernst_potential."""
    rand = [(_spp(rng), _n(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(-3.1, 2.9), _n(rng), rng.uniform(200.0, 1500.0), _pp(rng, -100, 100)) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 3)
        if kind == 0:                               # R·T past DBL_MAX
            edge.append((_spp(rng, -10, 1023), _n(rng), _pp(rng, 1019, 1023), _pp(rng)))
        elif kind == 1:                             # (R·T)/(n·F) subnormal, E° tiny or 0
            e0 = 0.0 if rng.random() < 0.3 else _spp(rng, -1074, -1000)
            edge.append((e0, _n(rng), _pp(rng, -1074, -1000), _pp(rng)))
        elif kind == 2:                             # the product with ln Q subnormal: Q near 1
            edge.append((_spp(rng, -1074, -1000), _n(rng), _pp(rng, -1050, -990), ulps(1.0, rng.randint(-40, 40))))
        else:                                       # the difference past ±DBL_MAX
            edge.append((DMAX * rng.choice((1.0, -1.0)) * rng.uniform(0.5, 1.0), _n(rng), _pp(rng, 1018, 1023),
                         _pp(rng)))
    bad = [(0.342, 2, 0.0, 0.1), (0.342, 0, 298.15, 0.1), (0.342, 2, 298.15, 0.0), (NAN, 2, 298.15, 0.1),
           (0.342, 2, INF, 0.1), (0.342, 2, DMAX, 1.0), (-0.0, 1, fb(1), 1.0), (0.0, 1, DMAX, 1.0), (0.342, 256, 298.15, 0.1)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range", edge), ("Err and specials", bad)]


def args_conc(rng):
    """(A, ε, l) for concentration_from_absorbance."""
    rand = [(_spp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(0.0, 3.0), _pp(rng, 0, 17), _pp(rng, -4, 4)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            comp.append((_spp(rng, 600, 1023), _pp(rng, 520, 1023), _pp(rng, 520, 1023)))
        else:
            comp.append((_spp(rng, -1074, -900), _pp(rng, -1074, -520), _pp(rng, -600, -10)))
    ties = []
    for _ in range(20_000):
        q, d = _sub_tie(rng)
        ties.append((q * rng.choice((1.0, -1.0)), 1.0, d))
    bad = [(1.0, 0.0, 1.0), (1.0, 1.0, 0.0), (1.0, NAN, 1.0), (INF, 1e200, 1e200), (-INF, 1.0, 1.0), (NAN, 1.0, 1.0),
           (-0.0, DMAX, DMAX), (1.0, DMAX, DMAX)]
    return [("binade-uniform", rand), ("typical", typical), ("ε·l past the range", comp), ("subnormal ties", ties),
            ("Err and specials", bad)]


def args_mul3_30(rng, moment=False):
    """(ε, l, c) for absorbance, or (μ, r) for diatomic_moment_of_inertia (moment)."""
    secs = args_mul3(rng)
    if not moment:
        return secs
    out = []
    for tag, sec in secs:
        out.append((tag, [(x, y) for x, y, _ in sec]))
    ties = []
    for _ in range(30_000):                         # 3·r·r with 3r² odd and of 54 bits
        r = float(_odd_in(rng, 54_794_159, 77_490_641))
        ties.append((3.0 * rng.choice((1.0, -1.0)), r * 2.0 ** rng.randint(-400, 400)))
    out.append(("midpoints 3·r·r", ties))
    return out


def args_rotor(rng):
    """(J, B) for rigid_rotor_energy."""
    def j(rng_):
        return rng_.choice((0, 1, 2, 3, 10, 0xFFFFFFFF, rng_.randint(0, 100), rng_.randint(0, 0xFFFFFFFF)))
    rand = [(j(rng), _spp(rng)) for _ in range(80_000)]
    typical = [(rng.randint(0, 100), _pp(rng, 30, 45)) for _ in range(60_000)]
    edge = [(j(rng), _spp(rng, -1074, -900)) for _ in range(60_000)]            # h·B below 2^-1022
    bad = [(-1, 1.0), (1 << 32, 1.0), (0, INF), (1, INF), (0, NAN), (2, -INF), (0xFFFFFFFF, DMAX), (0, -0.0), (5, -0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("h·B past the range", edge), ("domain and specials", bad)]


def args_rot(rng):
    """(I,) for rotational_constant."""
    rand = [(_pp(rng),) for _ in range(80_000)]
    typical = [(_pp(rng, -160, -145),) for _ in range(60_000)]
    edge = [(_pp(rng, -1074, -1020),) for _ in range(60_000)]                  # 8π²·I subnormal
    bad = [(0.0,), (-1.0,), (INF,), (NAN,), (fb(1),), (DMAX,)]
    return [("binade-uniform", rand), ("typical", typical), ("8π²·I subnormal", edge), ("Err and ends", bad)]


def audit_checks(work, targets, which):
    """Ids 33 … 57 against their exact references on every target, the exact path forced; True if anything failed."""
    failed = False
    gens = {"mono": args_mono, "heat": args_mul3, "shc": args_shc, "boil": args_mul3, "frz": args_mul3, "osm": args_osm,
            "zhl": args_zhl, "ssi": args_ssi, "pre": args_pre, "soc": args_soc, "mm": args_mm, "monotie": args_monotie,
            "sqrtdiv": args_sqrtdiv, "dil": args_dil, "cond": args_cond, "fmass": args_fmass, "fcharge": args_fcharge,
            "fratio": args_fratio, "nernst": args_nernst, "abs": args_mul3_30, "conc": args_conc,
            "vib": lambda r: args_sqrtdiv(r, 60_000), "rotor": args_rotor, "rot": args_rot,
            "moment": lambda r: args_mul3_30(r, True)}

    def eq(a, b):
        return a == b or (math.isnan(fb(a)) and math.isnan(fb(b)))
    for fid, name in AUDIT:
        if name not in which:
            continue
        t0 = time.time()
        rng = random.Random(1033 + fid)
        sections = gens[name](rng)
        xs = [x for _, sec in sections for x in sec]
        want = [ref_audit(name, a) for a in xs]
        first = None
        for t in targets:
            g0, g1, fl = run_probe(work, t, fid, 0, xs)
            a0, _, _ = run_probe(work, t, fid, 1, xs)
            bad = [i for i in range(len(xs)) if not eq(g0[i], want[i])]
            bad_acc = [i for i in range(len(xs)) if not eq(a0[i], want[i])]
            bad_fast = [i for i in range(len(xs)) if fl[i] == 1 and not eq(g1[i], g0[i])]
            for i in bad[:5]:
                print(f"  MISMATCH {name}{xs[i]!r} on {t}: got {g0[i]:016x}, want {want[i]:016x}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {name} (slow path){xs[i]!r} on {t}: got {a0[i]:016x}, want {want[i]:016x}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {name}{xs[i]!r} on {t}: {g1[i]:016x} vs {g0[i]:016x}")
            # A NaN word is one value across targets: heat / boil / frz keep Rust's (m·c)·ΔT for a non-finite input, and
            # the hardware's default NaN of ∞·0 is negative on x86_64 and positive on aarch64 (as Rust's own result is;
            # the parity rules match any NaN). Every other word must be identical.
            def canon(ws):
                return [NANW if math.isnan(fb(w)) else w for w in ws]
            cur = (canon(g0), canon(g1), fl, canon(a0))
            xdiff = 0 if first is None else sum(1 for u, v in zip(cur, first) for p, q in zip(u, v) if p != q)
            if first is None:
                first = cur
            slow = len(xs) - sum(1 for f in fl if f == 1)
            print(f"{name:7} {t:8} {len(xs):9,} args: {len(bad)} mismatches; slow path forced: {len(bad_acc)} "
                  f"mismatches; fast path: {len(bad_fast)} differ, fell back {slow:,}; bits or flags differing from "
                  f"{targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0] and name != "monotie":
                rates = []
                i = 0
                for tag, sec in sections:
                    rates.append(f"{tag} {sum(1 for f in fl[i:i + len(sec)] if f != 1)}/{len(sec)}")
                    i += len(sec)
                print("         fallbacks by region: " + "; ".join(rates))
        print(f"       [{time.time() - t0:.0f} s]")
    return failed

def gen_args(name, lim):
    rng = random.Random({"exp": 1001, "exp10": 1002, "expm1": 1003, "ln": 1004, "log10": 1005, "pow": 1006,
                         "sin": 1007, "cos": 1008, "atan2": 1009, "acos": 1010, "exp2s": 1011, "h100": 1012,
                         "uxsqrt": 1013, "uxmul": 1014, "uxdiv": 1015, "uxadd": 1016, "uxatan2": 1017, "dist": 1018,
                         "bond": 1019, "dih": 1020, "bohr": 1021, "trans": 1022, "ryd": 1022, "series": 1024,
                         "rmass": 1025, "ratio": 1026}[name])
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
    elif name == "h100":
        sections = args_h100(rng)
    elif name in ("uxsqrt", "uxmul", "uxdiv", "uxadd", "uxatan2"):
        sections = {"uxsqrt": args_uxsqrt, "uxmul": args_uxmul, "uxdiv": args_uxdiv, "uxadd": args_uxadd,
                    "uxatan2": args_uxatan2}[name](rng)
    elif name in ("dist", "bond", "dih"):
        sections = args_geom(rng, {"dist": 2, "bond": 3, "dih": 4}[name])
    elif name in SPEC_NAMES:
        sections = {"ratio": args_ratio, "bohr": args_bohr, "trans": args_trans, "ryd": args_trans,
                    "series": args_series, "rmass": args_rmass}[name](rng)
    else:
        sections = args_log(rng, lim, name == "log10")
    xs = [x for _, sec in sections for x in sec]
    spans = []
    i = 0
    for tag, sec in sections:
        spans.append((tag, i, i + len(sec)))
        i += len(sec)
    head = len(sections[0][1])
    acc = xs[:head] + [xs[i] for i in random.Random(77).sample(range(len(xs)), min(N_ACC, len(xs)))]
    return xs, acc, spans


PROBE = r'''# GENERATED by scripts/oracle/cr_verify.py. Reads job.bin (u64 words: function id, mode, n, then n arguments,
# two words each for pow and atan2) and writes res.bin. Mode 0: for each argument the function's result, the fast
# path's result and one byte, 1 if the fast path's test passed. Mode 1: the accurate path's result alone.
# Id 10 (sincos) gives sin and cos in the two slots, and the byte 1 or 0 as above, 2 if the fast pair differs
# from the full pair, 3 if the pair differs from _km_cr_sin and _km_cr_cos; id 11 is x·2^y (M4d), two words; id 12
# is 1/(100·x) (M4e).
include "src/error.cyr"
include "src/constants.cyr"
include "src/mathx_data.cyr"
include "src/mathx.cyr"
include "src/validate.cyr"
include "src/element_data.cyr"
include "src/element.cyr"
include "src/molecule.cyr"
include "src/stoichiometry.cyr"
include "src/thermo.cyr"
include "src/phase_data.cyr"
include "src/phase.cyr"
include "src/solution.cyr"
include "src/electrochemistry_data.cyr"
include "src/electrochemistry.cyr"
include "src/kinetics.cyr"
include "src/organic_data.cyr"
include "src/organic.cyr"
include "src/inorganic_data.cyr"
include "src/inorganic.cyr"
include "src/potential.cyr"
include "src/nuclear_data.cyr"
include "src/nuclear.cyr"
include "src/spectroscopy_data.cyr"
include "src/spectroscopy.cyr"

# The argument words of function fid.
#must_use
fn arity(fid): i64 {
    if (fid == 5 || fid == 8 || fid == 11 || fid == 13) { return 2; }
    if (fid >= 14 && fid <= 17) { return 4; }
    if (fid == 18) { return 6; }
    if (fid == 19) { return 9; }
    if (fid == 20) { return 12; }
    if (fid == 21 || fid == 24 || fid == 25) { return 2; }
    if (fid == 22 || fid == 23) { return 3; }
    if (fid == 26) { return 9; }
    if (fid == 27) { return 2; }
    if (fid == 28) { return 5; }
    if (fid == 29) { return 4; }
    if (fid == 30 || fid == 31) { return 7; }
    if (fid == 32) { return 3; }
    if (fid == 33) { return 12; }
    if (fid == 39) { return 2; }
    if (fid == 41) { return 4; }
    if (fid == 44) { return 8; }
    if (fid >= 34 && fid <= 43) { return 3; }
    if (fid == 45 || fid == 47 || fid == 54 || fid == 55 || fid == 57) { return 2; }
    if (fid == 50) { return 5; }
    if (fid == 51) { return 4; }
    if (fid == 56) { return 1; }
    if (fid >= 46 && fid <= 53) { return 3; }
    return 1;
}

# ADR 0002 #21, #22 (M4e-3): ids 21 … 25 are bohr_energy_level (Z, n), bohr_transition_energy and rydberg_wavelength
# (Z, n1, n2), spectral_series_wavelength (series, line) and reduced_mass (m1, m2): mode 0 the public function (an Err as
# +inf), its fast path in the second slot and the flag; mode 1 the exact path alone. Id 26 is _km_cr_ratio(p, q, e): p and
# q four little-endian 64-bit words each, then e.
#must_use
fn run_spec(fid, mode, a): (f64, f64, i64) {
    var inf: f64 = 0x7FF0000000000000;
    var w0 = load64(a);
    var w1 = load64(a + 8);
    var w2 = load64(a + 16);
    var ri: f64 = KM_RYDBERG;
    var rh: f64 = _KM_RYDBERG_HYDROGEN;
    if (fid == 26) {
        var p[64];
        var q[64];
        var t[64];
        var k = 3;
        _km_bn_set(&p, 0);
        _km_bn_set(&q, 0);
        while (k >= 0) {
            _km_bn_shl(&p, 64);
            _km_bn_set(&t, load64(a + 8 * k));
            _km_bn_add(&p, &p, &t);
            _km_bn_shl(&q, 64);
            _km_bn_set(&t, load64(a + 32 + 8 * k));
            _km_bn_add(&q, &q, &t);
            k = k - 1;
        }
        var rr: f64 = _km_cr_ratio(&p, &q, load64(a + 64));
        return (rr, rr, 1);
    }
    if (fid == 21) {
        if (mode == 1) { var la: f64 = _km_bohr_level_acc(w0, w1); var nla: f64 = -la; return (nla, nla, 0); }
        var tl, vl = bohr_energy_level(w0, w1);
        var gl: f64 = inf;
        if (is_ok(tl) == 1) { gl = vl; }
        var fl, ol = _km_bohr_level_fast(w0, w1);
        var nfl: f64 = -fl;
        return (gl, nfl, ol);
    }
    if (fid == 22) {
        if (mode == 1) { var ta: f64 = _km_bohr_transition_acc(w0, w1, w2); return (ta, ta, 0); }
        var tt, vt = bohr_transition_energy(w0, w1, w2);
        var gt: f64 = inf;
        if (is_ok(tt) == 1) { gt = vt; }
        var ft, ot = _km_bohr_transition_fast(w0, w1, w2);
        return (gt, ft, ot);
    }
    if (fid == 23) {
        if (mode == 1) { var ya: f64 = _km_rydberg_acc(ri, w0, w1, w2); return (ya, ya, 0); }
        var ty, vy = rydberg_wavelength(w0, w1, w2);
        var gy: f64 = inf;
        if (is_ok(ty) == 1) { gy = vy; }
        var fy, oy = _km_rydberg_fast(ri, w0, w1, w2);
        return (gy, fy, oy);
    }
    if (fid == 24) {
        var nl = spectral_series_n_lower(w0);
        var nu = nl + w1;
        if (mode == 1) {
            if (nu > 0xFFFFFFFF) { return (inf, inf, 0); }
            var sa: f64 = _km_rydberg_acc(rh, 1, nl, nu);
            return (sa, sa, 0);
        }
        var ts, vs = spectral_series_wavelength(w0, w1);
        var gs: f64 = inf;
        if (is_ok(ts) == 1) { gs = vs; }
        if (nu > 0xFFFFFFFF) { return (gs, gs, 1); }
        var fs, os = _km_rydberg_fast(rh, 1, nl, nu);
        return (gs, fs, os);
    }
    var m1: f64 = load64(a);
    var m2: f64 = load64(a + 8);
    if (mode == 1) { var ma: f64 = _km_reduced_mass_acc(m1, m2); return (ma, ma, 0); }
    var tm, vm = reduced_mass(m1, m2);
    var gm: f64 = inf;
    if (is_ok(tm) == 1) { gm = vm; }
    var fm, om = _km_reduced_mass_fast(m1, m2);
    return (gm, fm, om);
}

# ADR 0002 #24 … #26 (M4f-2): ids 27 … 32 are percent_yield (a, t), theoretical_yield_moles (m0, m1, c0, c1, p: two
# reactants), percent_composition (z0, c0, z1, c1: entry 0's percent), limiting_reagent (n, m0, m1, m2, c0, c1, c2),
# empirical_formula (n, z0, p0, z1, p1, z2, p2) and _km_cr_mul_div (x, y, z), as result words: a float function's Err is
# the NaN word 0x7FF8000000000001 (none of their results is NaN), limiting_reagent's Err −1, empirical_formula's Ok counts
# (c0 | c1·2^32, c2) and its Err (0xFFFFFFFFFFFFFF00 | code, −1). Mode 0: the public function, the fast path's result
# (_km_mul_div_fast; for limiting_reagent the index by exact comparisons alone) and its flag; mode 1: the exact path
# forced (_km_mul_div_acc at the same operands, the exact comparisons alone, _km_empirical_formula's exact tests).
#must_use
fn run_stoich(fid, mode, a): (i64, i64, i64) {
    var nanw = 0x7FF8000000000001;
    if (fid == 32) {
        var x: f64 = load64(a);
        var y: f64 = load64(a + 8);
        var z: f64 = load64(a + 16);
        if (mode == 1) {
            var xa: f64 = _km_mul_div_acc(x, y, z);
            return (_km_bits(xa), _km_bits(xa), 0);
        }
        var xr: f64 = _km_cr_mul_div(x, y, z);
        var xf, xo = _km_mul_div_fast(x, y, z);
        return (_km_bits(xr), _km_bits(xf), xo);
    }
    if (fid == 27) {
        var ya: f64 = load64(a);
        var yt: f64 = load64(a + 8);
        var hundred: f64 = 100.0;
        var tp, vp = percent_yield(ya, yt);
        if (is_ok(tp) != 1) { return (nanw, nanw, 1); }
        if (ya == 0.0) { return (vp, vp, 1); }
        if (mode == 1) {
            var pa: f64 = _km_mul_div_acc(ya, hundred, yt);
            return (_km_bits(pa), _km_bits(pa), 0);
        }
        var pf, po = _km_mul_div_fast(ya, hundred, yt);
        return (vp, _km_bits(pf), po);
    }
    if (fid == 28) {
        var pp: f64 = load64(a + 32);
        var ty, vy = theoretical_yield_moles(a, 2, a + 16, 2, pp);
        if (is_ok(ty) != 1) { return (nanw, nanw, 1); }
        var tl, il = limiting_reagent(a, 2, a + 16, 2);
        var mm: f64 = load64(a + 8 * il);
        var cc: f64 = load64(a + 16 + 8 * il);
        if (mm == 0.0) { return (vy, vy, 1); }
        if (mode == 1) {
            var ta: f64 = _km_mul_div_acc(mm, pp, cc);
            return (_km_bits(ta), _km_bits(ta), 0);
        }
        var tf, tfo = _km_mul_div_fast(mm, pp, cc);
        return (vy, _km_bits(tf), tfo);
    }
    if (fid == 29) {
        var pairs[32];
        store64(&pairs, load64(a));
        store64(&pairs + 8, load64(a + 8));
        store64(&pairs + 16, load64(a + 16));
        store64(&pairs + 24, load64(a + 24));
        var mol = KmMolecule_new(&pairs, 2);
        var ob[64];
        var tc, vc = percent_composition(mol, &ob, 4);
        if (is_err_result(tc) == 1) { return (nanw, nanw, 1); }
        var e0 = load64(&ob + 8);
        var cnt0 = load64(a + 8);
        if (cnt0 == 0) { return (e0, e0, 1); }
        var tw, vw = mol.molecular_weight();
        var mw: f64 = vw;
        var ms: f64 = KM_ELEMENT_ATOMIC_MASS[load64(a) - 1];
        var hc: f64 = f64_from(100 * cnt0);
        if (mode == 1) {
            var ca: f64 = _km_mul_div_acc(ms, hc, mw);
            return (_km_bits(ca), _km_bits(ca), 0);
        }
        var cf, cfo = _km_mul_div_fast(ms, hc, mw);
        return (e0, _km_bits(cf), cfo);
    }
    if (fid == 30) {
        var nl = load64(a);
        var tg, vg = limiting_reagent(a + 8, nl, a + 32, nl);
        var got = vg;
        if (is_err_result(tg) == 1) { got = 0 - 1; }
        var slow = got;
        if (got >= 0) {
            slow = 0;
            var i = 1;
            while (i < nl) {
                var mi: f64 = load64(a + 8 + 8 * i);
                var ci: f64 = load64(a + 32 + 8 * i);
                var ms0: f64 = load64(a + 8 + 8 * slow);
                var cs0: f64 = load64(a + 32 + 8 * slow);
                if (_km_ratio_lt(mi, ci, ms0, cs0) == 1) { slow = i; }
                i = i + 1;
            }
        }
        if (mode == 1) { return (slow, slow, 0); }
        return (got, slow, 1);
    }
    var ne = load64(a);
    var eb[64];
    var k = 0;
    while (k < 8) {
        store64(&eb + 8 * k, 0);
        k = k + 1;
    }
    var force = 0;
    if (mode == 1) { force = 1; }
    var te, ve = _km_empirical_formula(a + 8, ne, &eb, 4, force);
    if (is_err_result(te) == 1) { return ((0 - 256) | ve, 0 - 1, 1); }
    var w0 = load64(&eb + 8) | (load64(&eb + 24) << 32);
    return (w0, load64(&eb + 40), 1);
}

# The overflow audit's B2 (ADR 0002 #29, #30): ids 45 … 57 as result words (cr_verify.py's AUDIT): mode 0 the public
# function (an Err as the NaN word 0x7FF8000000000001), its shape's fast path at the function's factors and its flag
# (_km_sqrt_div_fast for 45 and 54, _km_nernst_fast for 51); mode 1 the slow path forced (_km_cr_mono_acc at the factors,
# _km_sqrt_div_acc, _km_nernst_ux). Id 45 is _km_cr_sqrt_div itself (its NaN outside the domain in both slots). An
# electron count and J are integer words. An Err row, and a non-finite input of an unvalidated function, return the
# public result in both slots and modes.
#must_use
fn run_audit2(fid, mode, a): (i64, i64, i64) {
    var nanw = 0x7FF8000000000001;
    var x: f64 = load64(a);
    var y: f64 = load64(a + 8);
    var z: f64 = load64(a + 16);
    var w: f64 = load64(a + 24);
    var wi1 = load64(a + 8);
    var wi2 = load64(a + 16);
    var fc: f64 = KM_FARADAY;
    var hc: f64 = KM_PLANCK;
    if (fid == 45 || fid == 54) {
        var pw = 0;
        if (fid == 45) {
            var r45: f64 = _km_cr_sqrt_div(x, y);
            pw = _km_bits(r45);
            if (_km_is_nan(r45) == 1) { return (nanw, nanw, 1); }
        } else {
            var t54, v54 = vibrational_frequency(x, y);
            if (is_ok(t54) != 1) { return (nanw, nanw, 1); }
            pw = v54;
        }
        if (mode == 1) {
            var sa: f64 = _km_sqrt_div_acc(x, y);
            return (_km_bits(sa), _km_bits(sa), 0);
        }
        var sf, so = _km_sqrt_div_fast(x, y);
        return (pw, _km_bits(sf), so);
    }
    if (fid == 51) {
        var t51, v51 = nernst_potential(x, wi1, z, w);
        if (is_ok(t51) != 1) { return (nanw, nanw, 1); }
        var nf: f64 = f64_from(wi1);
        var lq: f64 = _km_cr_ln(w);
        if (mode == 1) {
            var nu: f64 = _km_nernst_ux(x, nf, z, lq);
            return (_km_bits(nu), _km_bits(nu), 0);
        }
        var nfa, nok = _km_nernst_fast(x, nf, z, lq);
        return (v51, _km_bits(nfa), nok);
    }
    var nb: f64[4];
    var db: f64[4];
    var nn = 0;
    var nd = 0;
    var pw2 = 0;
    var ok = 1;
    var ff: f64 = 0.0;
    var fo = 0;
    if (fid == 46) {
        var t46, v46 = dilution(x, y, z);
        if (is_ok(t46) != 1) { return (nanw, nanw, 1); }
        pw2 = v46;
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&db, z);
        nn = 2;
        nd = 1;
        var f46, o46 = _km_mul_div_fast(x, y, z);
        ff = f46;
        fo = o46;
    } elif (fid == 47) {
        var r47: f64 = solution_conductivity(x, y);
        pw2 = _km_bits(r47);
        if (_km_is_finite(x) != 1 || _km_is_finite(y) != 1) { ok = 0; }
        var th: f64 = 1000.0;
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&db, th);
        nn = 2;
        nd = 1;
        var f47, o47 = _km_mul_div_fast(x, y, 1000.0);
        ff = f47;
        fo = o47;
    } elif (fid == 48 || fid == 49) {
        var nf2: f64 = f64_from(wi2);
        if (fid == 48) {
            var t48, v48 = faraday_mass_deposited(x, y, wi2);
            if (is_ok(t48) != 1) { return (nanw, nanw, 1); }
            pw2 = v48;
            store64(&nb, x);
            store64(&nb + 8, y);
            store64(&db, nf2);
            store64(&db + 8, fc);
            nn = 2;
            nd = 2;
            var f48, o48 = _km_mul2_div2_fast(x, y, nf2, fc);
            ff = f48;
            fo = o48;
        } else {
            var t49, v49 = faraday_charge_required(x, y, wi2);
            if (is_ok(t49) != 1) { return (nanw, nanw, 1); }
            pw2 = v49;
            store64(&nb, x);
            store64(&nb + 8, nf2);
            store64(&nb + 16, fc);
            store64(&db, y);
            nn = 3;
            nd = 1;
            var f49, o49 = _km_mul3_div_fast(x, nf2, fc, y);
            ff = f49;
            fo = o49;
        }
    } elif (fid == 50) {
        var wi4 = load64(a + 32);
        var t50, v50 = faraday_mass_ratio(x, y, wi2, w, wi4);
        if (is_ok(t50) != 1) { return (nanw, nanw, 1); }
        pw2 = v50;
        var na: f64 = f64_from(wi2);
        var nbf: f64 = f64_from(wi4);
        store64(&nb, x);
        store64(&nb + 8, w);
        store64(&nb + 16, na);
        store64(&db, nbf);
        store64(&db + 8, y);
        nn = 3;
        nd = 2;
        var f50, o50 = _km_mul3_div2_fast(x, w, na, nbf, y);
        ff = f50;
        fo = o50;
    } elif (fid == 52 || fid == 57) {
        var zz: f64 = z;
        var r52: f64 = 0.0;
        if (fid == 52) {
            r52 = absorbance(x, y, z);
        } else {
            zz = y;
            r52 = diatomic_moment_of_inertia(x, y);
        }
        pw2 = _km_bits(r52);
        if (_km_is_finite(x) != 1 || _km_is_finite(y) != 1 || _km_is_finite(zz) != 1) { ok = 0; }
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&nb + 16, zz);
        nn = 3;
        var f52, o52 = _km_mul3_fast(x, y, zz);
        ff = f52;
        fo = o52;
    } elif (fid == 53) {
        var t53, v53 = concentration_from_absorbance(x, y, z);
        if (is_ok(t53) != 1) { return (nanw, nanw, 1); }
        pw2 = v53;
        if (_km_is_finite(x) != 1) { ok = 0; }
        store64(&nb, x);
        store64(&db, y);
        store64(&db + 8, z);
        nn = 1;
        nd = 2;
        var f53, o53 = _km_div2_fast(x, y, z);
        ff = f53;
        fo = o53;
    } elif (fid == 55) {
        var j = load64(a);
        var r55: f64 = rigid_rotor_energy(j, y);
        pw2 = _km_bits(r55);
        if (_km_is_u32(j) != 1 || _km_is_finite(y) != 1) { ok = 0; }
        var jf: f64 = f64_from(j);
        var j1: f64 = jf + 1.0;
        store64(&nb, hc);
        store64(&nb + 8, y);
        store64(&nb + 16, jf);
        store64(&nb + 24, j1);
        nn = 4;
        var f55, o55 = _km_mul4_fast(hc, y, jf, j1);
        ff = f55;
        fo = o55;
    } else {
        var t56, v56 = rotational_constant(x);
        if (is_ok(t56) != 1) { return (nanw, nanw, 1); }
        pw2 = v56;
        var pi: f64 = 0x400921FB54442D18;
        var eight: f64 = 8.0;
        store64(&nb, hc);
        store64(&db, eight);
        store64(&db + 8, pi);
        store64(&db + 16, pi);
        store64(&db + 24, x);
        nn = 1;
        nd = 4;
        var h8: f64 = hc * 0.125;
        var f56, o56 = _km_div3_fast(h8, pi, pi, x);
        ff = f56;
        fo = o56;
    }
    if (ok != 1) { return (pw2, pw2, 1); }
    if (mode == 1) {
        var aa: f64 = _km_cr_mono_acc(&nb, nn, &db, nd);
        return (_km_bits(aa), _km_bits(aa), 0);
    }
    return (pw2, _km_bits(ff), fo);
}

# The overflow audit's B1 (ADR 0002 #27, #28): ids 33 … 44 as result words (cr_verify.py's AUDIT): mode 0 the public
# function (an Err as the NaN word 0x7FF8000000000001), the fast path's result at the function's factors and its flag
# (the monomial's _km_cr_mono_fast; for second_order_concentration and michaelis_menten their Rust-chain fast path);
# mode 1 the slow path forced (_km_cr_mono_acc; _km_second_order_ux, _km_michaelis_menten_ux). Id 33 is _km_cr_mono
# itself (n_num, n_den, 6 + 4 factor words), id 44 _km_mono_tie (cn, cd, s0 … s2, t0, t1, z1; z2 = nextup(z1)). An Err
# row, and a non-finite input of the unvalidated functions, return the public result in both slots and modes.
#must_use
fn run_audit(fid, mode, a): (i64, i64, i64) {
    var nanw = 0x7FF8000000000001;
    if (fid >= 45) {
        var b0, b1, bf = run_audit2(fid, mode, a);
        return (b0, b1, bf);
    }
    if (fid == 33) {
        var nn = load64(a);
        var nd = load64(a + 8);
        if (mode == 1) {
            var ma: f64 = _km_cr_mono_acc(a + 16, nn, a + 64, nd);
            return (_km_bits(ma), _km_bits(ma), 0);
        }
        var mr: f64 = _km_cr_mono(a + 16, nn, a + 64, nd);
        var mf, mo = _km_cr_mono_fast(a + 16, nn, a + 64, nd);
        return (_km_bits(mr), _km_bits(mf), mo);
    }
    if (fid == 44) {
        var z1: f64 = load64(a + 56);
        var zb = _km_bits(z1) + 1;
        var z2: f64 = zb;
        var tr: f64 = _km_mono_tie(a + 16, load64(a), a + 40, load64(a + 8), z1, z2);
        return (_km_bits(tr), _km_bits(tr), 1);
    }
    var x: f64 = load64(a);
    var y: f64 = load64(a + 8);
    var z: f64 = load64(a + 16);
    var w: f64 = load64(a + 24);
    var nb: f64[5];
    var db: f64[2];
    var nn2 = 3;
    var nd2 = 0;
    var pw = 0;
    var ok = 1;
    if (fid == 34 || fid == 36 || fid == 37) {
        var r3: f64 = 0.0;
        if (fid == 34) { r3 = heat_transfer(x, y, z); }
        if (fid == 36) { r3 = boiling_point_elevation(x, y, z); }
        if (fid == 37) { r3 = freezing_point_depression(x, y, z); }
        pw = _km_bits(r3);
        if (_km_is_finite(x) != 1 || _km_is_finite(y) != 1 || _km_is_finite(z) != 1) { ok = 0; }
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&nb + 16, z);
    } elif (fid == 35) {
        var t35, v35 = specific_heat_from_calorimetry(x, y, z);
        pw = v35;
        if (is_ok(t35) != 1) { pw = nanw; ok = 0; }
        store64(&nb, x);
        store64(&db, y);
        store64(&db + 8, z);
        nn2 = 1;
        nd2 = 2;
    } elif (fid == 38) {
        var t38, v38 = osmotic_pressure(x, y, z);
        pw = v38;
        if (is_ok(t38) != 1) { pw = nanw; ok = 0; }
        var rg: f64 = KM_GAS_CONSTANT;
        var th: f64 = 1000.0;
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&nb + 16, rg);
        store64(&nb + 24, th);
        store64(&nb + 32, z);
        nn2 = 5;
    } elif (fid == 39) {
        var t39, v39 = zero_order_half_life(x, y);
        pw = v39;
        if (is_ok(t39) != 1) { pw = nanw; ok = 0; }
        var two: f64 = 2.0;
        store64(&nb, x);
        store64(&db, two);
        store64(&db + 8, y);
        nn2 = 1;
        nd2 = 2;
    } elif (fid == 40) {
        var t40, v40 = steady_state_intermediate(x, y, z);
        pw = v40;
        if (is_ok(t40) != 1) { pw = nanw; ok = 0; }
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&db, z);
        nn2 = 2;
        nd2 = 1;
    } elif (fid == 41) {
        var t41, v41 = pre_equilibrium_rate(x, y, z, w);
        pw = v41;
        if (is_ok(t41) != 1) { pw = nanw; ok = 0; }
        store64(&nb, z);
        store64(&nb + 8, x);
        store64(&nb + 16, w);
        store64(&db, y);
        nd2 = 1;
    } elif (fid == 42) {
        var t42, v42 = second_order_concentration(x, y, z);
        if (is_ok(t42) != 1) { return (nanw, nanw, 1); }
        if (mode == 1) {
            var su: f64 = _km_second_order_ux(x, y, z);
            return (_km_bits(su), _km_bits(su), 0);
        }
        var sf, so = _km_second_order_fast(x, y, z);
        return (v42, _km_bits(sf), so);
    } else {
        var t43, v43 = michaelis_menten(x, y, z);
        if (is_ok(t43) != 1) { return (nanw, nanw, 1); }
        if (mode == 1) {
            var mu: f64 = _km_michaelis_menten_ux(x, y, z);
            return (_km_bits(mu), _km_bits(mu), 0);
        }
        var vf, vo = _km_michaelis_menten_fast(x, y, z);
        return (v43, _km_bits(vf), vo);
    }
    if (ok != 1) { return (pw, pw, 1); }
    if (mode == 1) {
        var aa: f64 = _km_cr_mono_acc(&nb, nn2, &db, nd2);
        return (_km_bits(aa), _km_bits(aa), 0);
    }
    var ff: f64 = 0.0;
    var fo = 0;
    if (fid == 34 || fid == 36 || fid == 37) {
        var f3, o3 = _km_mul3_fast(x, y, z);
        ff = f3;
        fo = o3;
    } elif (fid == 35) {
        var f5, o5 = _km_div2_fast(x, y, z);
        ff = f5;
        fo = o5;
    } elif (fid == 38) {
        var rg5: f64 = KM_GAS_CONSTANT;
        var f8, o8 = _km_mul5_fast(x, y, rg5, 1000.0, z);
        ff = f8;
        fo = o8;
    } elif (fid == 39) {
        var f9, o9 = _km_div2_fast(x, 2.0, y);
        ff = f9;
        fo = o9;
    } elif (fid == 40) {
        var f0, o0 = _km_mul_div_fast(x, y, z);
        ff = f0;
        fo = o0;
    } else {
        var f1, o1 = _km_mul3_div_fast(z, x, w, y);
        ff = f1;
        fo = o1;
    }
    return (pw, _km_bits(ff), fo);
}

# ADR 0002 #20 (M4e-2): ids 13 … 17 take (significand, exponent) words: _km_ux_sqrt_f64(s, e), _km_ux_mul_f64,
# _km_ux_div_f64, _km_ux_add_f64 and _km_ux_atan2_f64 (mode 1: _km_ux_atan2_acc; mode 0's second slot: _km_ux_atan2_fast);
# ids 18 … 20 the geometry's points (6, 9, 12 f64): mode 0 the public function, its second slot the slow path forced
# (_km_distance_3d_ux, _km_bond_angle_ux, _km_dihedral_angle_ux), mode 1 the forced slow path alone. An Err is +inf.
#must_use
fn run_ux(fid, mode, a): (f64, f64, i64) {
    var s1: f64 = load64(a);
    var e1 = load64(a + 8);
    var s2: f64 = load64(a + 16);
    var e2 = load64(a + 24);
    if (fid == 13) { var q: f64 = _km_ux_sqrt_f64(s1, e1); return (q, q, 1); }
    if (fid == 14) { var m: f64 = _km_ux_mul_f64(s1, e1, s2, e2); return (m, m, 1); }
    if (fid == 15) { var d: f64 = _km_ux_div_f64(s1, e1, s2, e2); return (d, d, 1); }
    if (fid == 16) { var p: f64 = _km_ux_add_f64(s1, e1, s2, e2); return (p, p, 1); }
    if (fid == 17) {
        if (mode == 1) { var t1: f64 = _km_ux_atan2_acc(s1, e1, s2, e2); return (t1, t1, 0); }
        var t0: f64 = _km_ux_atan2_f64(s1, e1, s2, e2);
        var tf, tok = _km_ux_atan2_fast(s1, e1, s2, e2);
        return (t0, tf, tok);
    }
    var inf: f64 = 0x7FF0000000000000;
    if (fid == 18) {
        var ds: f64 = _km_distance_3d_ux(a, a + 24);
        if (mode == 1) { return (ds, ds, 0); }
        var dp: f64 = distance_3d(a, a + 24);
        return (dp, ds, 1);
    }
    if (fid == 19) {
        var bs: f64 = inf;
        var tb, vb = _km_bond_angle_ux(a, a + 24, a + 48);
        if (is_ok(tb) == 1) { bs = vb; }
        if (mode == 1) { return (bs, bs, 0); }
        var bp: f64 = inf;
        var tc, vc = bond_angle(a, a + 24, a + 48);
        if (is_ok(tc) == 1) { bp = vc; }
        return (bp, bs, 1);
    }
    var hs: f64 = inf;
    var th, vh = _km_dihedral_angle_ux(a, a + 24, a + 48, a + 72);
    if (is_ok(th) == 1) { hs = vh; }
    if (mode == 1) { return (hs, hs, 0); }
    var hp: f64 = inf;
    var ti, vi = dihedral_angle(a, a + 24, a + 48, a + 72);
    if (is_ok(ti) == 1) { hp = vi; }
    return (hp, hs, 1);
}

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
        if (fid == 12) { var ha: f64 = _km_cr_hundredth_div_acc(x); return (ha, ha, 0); }
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
    if (fid == 12) {
        var hh: f64 = _km_cr_hundredth_div(x);
        var hf, ho = _km_cr_hundredth_div_fast(x);
        return (hh, hf, ho);
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
    var cap = 8 * (3 + WMAX);
    var buf = alloc(cap);
    var got = file_read_all("job.bin", buf, cap);
    if (got < 24) { return 2; }
    var fid = load64(buf);
    var mode = load64(buf + 8);
    var n = load64(buf + 16);
    var ar = arity(fid);
    if (got != 8 * (3 + ar * n)) { return 3; }
    var out = alloc(17 * n + 8);
    var i = 0;
    while (i < n) {
        var at = buf + 24 + 8 * ar * i;
        var x: f64 = load64(at);
        var y: f64 = 0.0;
        if (ar == 2) { y = load64(at + 8); }
        var r: f64 = 0.0;
        var rf: f64 = 0.0;
        var ok = 0;
        if (fid >= 33) {
            var u0, u1, uf = run_audit(fid, mode, at);
            store64(out + 8 * i, u0);
            store64(out + 8 * n + 8 * i, u1);
            store8(out + 16 * n + i, uf);
        } elif (fid >= 27) {
            var w0, w1, wf = run_stoich(fid, mode, at);
            store64(out + 8 * i, w0);
            store64(out + 8 * n + 8 * i, w1);
            store8(out + 16 * n + i, wf);
        } elif (fid >= 21) {
            var r3, rf3, ok3 = run_spec(fid, mode, at);
            r = r3;
            rf = rf3;
            ok = ok3;
        } elif (fid >= 13) {
            var r1, rf1, ok1 = run_ux(fid, mode, at);
            r = r1;
            rf = rf1;
            ok = ok1;
        } else {
            var r2, rf2, ok2 = run(fid, mode, x, y);
            r = r2;
            rf = rf2;
            ok = ok2;
        }
        if (fid < 27) {
            store64(out + 8 * i, r);
            store64(out + 8 * n + 8 * i, rf);
            store8(out + 16 * n + i, ok);
        }
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


def build(root, work, wmax, targets):
    src = work / "probe.cyr"
    src.write_text(PROBE.replace("WMAX", str(wmax)))
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
    flat = [v & 0xFFFFFFFFFFFFFFFF if isinstance(v, int) else bits(v)
            for a in args for v in (a if isinstance(a, tuple) else (a,))]
    if ARITY.get(fid, 1) > 1 and args and not (isinstance(args[0], tuple) and len(args[0]) == ARITY[fid]):
        raise ValueError(f"{ARITY[fid]}-word arguments expected")
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
    if name in UX_NAMES:
        return ref_ux(name, a)
    if name in SPEC_NAMES:
        return ref_spec(name, a)
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
    stoich = [name for _, name in STOICH if only is None or name in only]
    audit = [name for _, name in AUDIT if only is None or name in only]
    plan = {}
    wmax = 7 * 600_000 if stoich else 0
    if audit:
        wmax = max(wmax, 12 * 450_000)
    for fid, name, _, _, ref in funcs:
        xs, acc, spans = gen_args(name, lim)
        plan[name] = (xs, acc, spans)
        wmax = max(wmax, ARITY.get(fid, 1) * len(xs))
    if big:
        wmax = max(wmax, 2 * 10_000_000)
    build(root, work, wmax, targets)
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
        if big and "elf" in targets and name not in ("pow", "atan2", "exp2s") + UX_NAMES + SPEC_NAMES:
            rng = random.Random(9000 + fid)
            if name in ("ln", "log10"):
                ys = binade_uniform(rng, 10_000_000, -1074, 1023, signs=(1,))
            elif name in ("sin", "cos"):
                ys = [rng.uniform(-8.0, 8.0) for _ in range(5_000_000)] + binade_uniform(rng, 5_000_000, -1074, 1023)
            elif name == "acos":
                ys = [rng.uniform(-1.0, 1.0) for _ in range(10_000_000)]
            elif name == "h100":
                ys = binade_uniform(rng, 10_000_000, -1074, 1023)
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
    failed |= stoich_checks(work, targets, stoich)
    failed |= audit_checks(work, targets, audit)
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
