"""Verify the correctly rounded exp, exp10, expm1, ln, log10, pow, sin, cos, atan2, acos, x·2^y and 1/(100·x) of
src/mathx.cyr against MPFR (ADR 0006), and its unbounded-exponent final operations with potential's geometry (ADR 0002
#20); and the exact rationals of spectroscopy (#21, #22), stoichiometry (#24 ... #26) and the overflow audit's B1 (#27,
#28: the correctly rounded monomial, thermo's, phase's and kinetics' rationals), B2 (#29, #30: solution's,
electrochemistry's and spectroscopy's, and the correctly rounded square root of a quotient) and B3 (#31, #32: gas's,
potential's, inorganic's and nuclear's, and the correctly rounded w·√(x/y·2^e)), B4 (#33, #34: the sums with
cancellation and the Shomate polynomials) and B5 (#35: the exp / exp10 arguments) against Fraction, exact integer square
roots and MPFR; and B10a's e^x, e^x − 1 and 10^x with an unbounded exponent against MPFR at its widest exponent range.

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
The overflow audit's B3 (ADR 0002 #31, #32) adds wsqrt (_km_cr_mul_sqrt_ratio = RN(w·sqrt((x/y)·2^e)): every binade of
each, e over [-4, 4] and far, x/y past either end, results past DBL_MAX and subnormal, exact roots, the 53-bit midpoints
3w·2^k (x/y = 9·4^k) and the subnormal grid's ties, roots within ulps of a midpoint, the domain; the fast path and
_km_mul_sqrt_ratio_acc forced; reference an exact integer square root) and the twelve functions: igp, igv, dens, coul,
force (ideal_gas_pressure, ideal_gas_volume, gas_density, coulomb, coulomb_force: the public function, its shape's fast
path, the monomial's exact path forced; reference Fraction), width (morse_width_parameter: the root's two paths) and
vdw, prp, rkp, lj, lande, bea (van_der_waals_pressure, peng_robinson_pressure, redlich_kwong_pressure, lennard_jones,
born_lande_lattice_energy, binding_energy_per_nucleon: the public function, its Rust-chain fast path where the range
test passes, the unbounded-exponent slow path forced on every argument; reference mpfr_ref.UX, Tc^2.5 by UX.pow), each
on binade-uniform, typical, past-the-range, co-volume (V near n·b or b), midpoint, subnormal and Err arguments.
The overflow audit's B4 (ADR 0002 #33, #34) adds the twelve families of the sums with cancellation and the Shomate
polynomials: hess, eff, rsum, qv, ionic, haber, gibbs, kohl, lim, zoc (hess_law, enthalpy_from_formation, thermochem's
reaction_enthalpy / reaction_gibbs_energy / reaction_entropy, q_value, ionic_strength, born_haber_lattice_energy,
gibbs_free_energy, kohlrausch_molar_conductivity, limiting_molar_conductivity, zero_order_concentration) and shcp, shh
(Shomate::cp, enthalpy_antiderivative): the public function, its Rust-chain fast path where the test passes, the
unbounded-exponent slow path forced on every argument; reference mpfr_ref.UX with mpfr_ref.ux_sum in Rust's summation
order; each on binade-uniform, typical, partial-sums-past-DBL_MAX, subnormal and special arguments. B4 also adds to
nernst the products c·l that IEEE rounds up to 2^-1022 from below (ADR 0002 #29's correction).
The overflow audit's B5 (ADR 0002 #35) adds keq, dhl, dhe, dav (equilibrium_constant, debye_huckel_limiting,
debye_huckel_extended, davies_activity: the public function, its Rust-chain fast path where the test passes, the
unbounded-exponent argument forced on every argument, then the correctly rounded exp or 10^x); reference the argument
on mpfr_ref.UX with every operation at 53 bits and MPFR's exp / exp10 of it rounded once (mpfr_ref.ux_exp_final,
ux_exp10_final), so the lemma that the port's binary64 argument gives the same result is checked on every argument; each
on binade-uniform, typical, R·T past DBL_MAX or subnormal, A·z², B·å or the numerator past DBL_MAX (at I = 0, or brought
back), subnormal products, Davies' bracket at 0, the 2^-1022 boundaries, and Err arguments.
The overflow audit's B6 (ADR 0002 #36 … #38) adds fft1 … fft256 (fft_spectrum of n = 1, 2, 4, …, 64 and 256 values: the
public spectrum, Rust's IEEE one where the port's fast test passes, the unbounded-exponent transform forced; reference
hisab 1.4.0's transform on mpfr_ref.UX with MPFR's twiddles, each √ rounded once; 110,718 signals: binade-uniform,
typical, huge, subnormal, mixed, near-cancelling, at the fast test's bounds, non-finite samples, and n = 4's product just
below 2^-1022), pk3 … pk16 (find_peaks of 3, 5, 8 and 16 values: the public function, the fast test's scan, every
candidate through the exact predicate; reference Fraction; 380,000 arrays with near-ties, subnormal data and ±DBL_MAX
ends), flame (adiabatic_flame_temperature and its transcribed bisection with the exact sign tests, 14,010 arguments
at tiny scales, endothermic and subnormal ΔH; reference Rust's solver in Python on the generated SHOMATE table) and
signsc (src/mathx.cyr's _km_sign_scaled, 140,110 sums against Fraction).
The overflow audit's B7 (ADR 0002 #39) adds bal23, bal34, bal56, bal89, bal1213 and bal2425 (balance_equation of 2 × 3,
3 × 4, 5 × 6, 8 × 9, 12 × 13 and 24 × 25 integer matrices, n + 1 words a result: the public result, the i64 fast path's
own outcome and the multi-word path forced; reference Rust's algorithm on Python ints with a tracer of its checked i64
operations, which the fast path's overflow report must match; 54,306 matrices over known x with large and small
entries, x with a zero, negative or u32-edge entry, random small and large entries, ±2^53 entries and dependent rows).
The overflow audit's B10a (ADR 0002 #40's primitives) adds exp_ux, expm1_ux and exp10_ux (_km_cr_exp_ux, _km_cr_expm1_ux
and _km_cr_exp10_ux of (s, e) pairs, the result a pair: the function, its fast path and the accurate path forced, and at
every normal double x the binary64 function's pair; reference MPFR at its widest exponent range, mpfr_ref.UX.exp, .expm1
and .exp10, whose own binary64 rounding is checked equal wherever that is normal; 448,188 arguments: specials and extreme
exponents up to ±2^63, every binade from 2^-1200 to 2^65 of both signs, binary64's range, the tiny end, each tier's edges,
near multiples of ln2/512 (log10(2)/512) with K up to 2^71, Lefevre–Muller's hard cases, large arguments to 3.2·10^18,
the saturation edges, expm1's −1 region and exp10's integers −400 … 400 and the tie 10^23).
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
                                                           # conc vib rotor rot moment wsqrt igp igv dens
                                                           # vdw prp rkp lj width coul force lande bea
                                                           # hess eff rsum qv ionic haber gibbs kohl lim zoc
                                                           # shcp shh keq dhl dhe dav flame pk3 pk5 pk8
                                                           # pk16 signsc fft1 … fft256 bal23 bal34 bal56
                                                           # bal89 bal1213 bal2425 exp_ux expm1_ux
                                                           # exp10_ux)
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
         (52, "abs"), (53, "conc"), (54, "vib"), (55, "rotor"), (56, "rot"), (57, "moment"),
         (58, "wsqrt"), (59, "igp"), (60, "igv"), (61, "dens"), (62, "vdw"), (63, "prp"), (64, "rkp"), (65, "lj"),
         (66, "width"), (67, "coul"), (68, "force"), (69, "lande"), (70, "bea"),
         (71, "hess"), (72, "eff"), (73, "rsum"), (74, "qv"), (75, "ionic"), (76, "haber"), (77, "gibbs"), (78, "kohl"),
         (79, "lim"), (80, "zoc"), (81, "shcp"), (82, "shh"),
         (83, "keq"), (84, "dhl"), (85, "dhe"), (86, "dav"),
         (87, "flame"), (88, "pk3"), (89, "pk5"), (90, "pk8"), (91, "pk16"), (100, "signsc")]
AUDIT_NAMES = tuple(name for _, name in AUDIT)
AUDIT2_NAMES = ("sqrtdiv", "dil", "cond", "fmass", "fcharge", "fratio", "nernst", "abs", "conc", "vib", "rotor", "rot",
                "moment")
ARITY.update({33: 12, 34: 3, 35: 3, 36: 3, 37: 3, 38: 3, 39: 2, 40: 3, 41: 4, 42: 3, 43: 3, 44: 8})
ARITY.update({45: 2, 46: 3, 47: 2, 48: 3, 49: 3, 50: 5, 51: 4, 52: 3, 53: 3, 54: 2, 55: 2, 56: 1, 57: 2})
ARITY.update({58: 4, 59: 3, 60: 3, 61: 3, 62: 5, 63: 5, 64: 4, 65: 3, 66: 3, 67: 3, 68: 3, 69: 5, 70: 3})
ARITY.update({71: 6, 72: 8, 73: 11, 74: 8, 75: 7, 76: 5, 77: 3, 78: 3, 79: 4, 80: 3, 81: 6, 82: 6})
ARITY.update({83: 2, 84: 3, 85: 5, 86: 3})
ARITY.update({87: 8, 88: 4, 89: 6, 90: 9, 91: 17, 100: 17})
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
    if name in AUDIT3_NAMES:
        return ref_audit3(name, a)
    if name in AUDIT4_NAMES:
        return ref_audit4(name, a)
    if name in AUDIT5_NAMES:
        return ref_audit5(name, a)
    if name in AUDIT6_NAMES:
        return ref_audit6(name, a)
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
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range", edge), ("Err and specials", bad),
            ("the product at 2^-1022 (B4)", args_nernst_2m1022(rng))]


# B4's correction of ADR 0002 #29 (the golden rows' NERNST_2M1022 in scripts/oracle/src/audit.rs, then a seeded search):
# (E°, n, T, Q) whose product c·l = ((R·T)/(n·F))·ln Q lies in [2^-1022 − 2^-1075, 2^-1022 − 2^-1076), which IEEE rounds
# up to 2^-1022 and 53-bit rounding to 2^-1022 − 2^-1075; B2's fast path took 2^-1022 as normal. E° an odd or even
# multiple of 2^-1074 (#20's value differs from Rust's for an odd one: a subnormal tie, to even), 0, or normal.
NERNST_2M1022 = [
    (0x0000000000000001, 1, 0x012c5299aebb56a1, 0x3ff0d2128b2f3300),
    (0x0000000000000001, 1, 0x015da605365dec16, 0x3ff018892f902b00),
    (0x0000000000000001, 1, 0x013b5b45c31413cc, 0x3ff06b6f03675a00),
    (0x8000000000000001, 7, 0x012d4ef95f210a80, 0x3fe6d09003db9e8a),
    (0x8000000000000001, 1, 0x01360e4f60503912, 0x3feefd1bc8838237),
    (0x8000000000000001, 255, 0x018459b4eaacdf85, 0x3fe83fd35aa3b4ff),
    (0x8000000000000001, 255, 0x0191b54bd2465e98, 0x3feb49320184af3e),
    (0x0000000000000001, 2, 0x011fcd71eb553fde, 0x3ff31ed6b3937b45),
    (0x0000000000000001, 3, 0x015af51773707b06, 0x3ff05183cbac8628),
    (0x8000000000000001, 7, 0x012b48b8639e6baa, 0x3fe63fc4c95cfefc),
    (0x0000000000000001, 1, 0x0111b58b7e7f1aa1, 0x3ff2c69c5171ec09),
    (0x8000000000000001, 1, 0x01223118d417935f, 0x3fed9a4e14f884e7),
    (0x0000000000000001, 3, 0x012cc67d7a4bbdf0, 0x3ff28bdf0a3b5122),
    (0x0000000000000001, 255, 0x01791b07966e9d3b, 0x3ff915712636d1db),
]


def args_nernst_2m1022(rng):
    """(E°, n, T, Q) with c·l in the window above: the fixed rows, then a seeded search (T over the binades where c is
    about 2^-1018, Q within 60 ulps of e^(l*) for the target l*), each with E° from {±2^-1074, ±3·2^-1074, 0, ±2^-1073,
    a normal one}."""
    lo = Fraction(2) ** -1022 - Fraction(2) ** -1075
    hi = Fraction(2) ** -1022 - Fraction(2) ** -1076
    base = [(n, fb(t), fb(q)) for _, n, t, q in NERNST_2M1022]
    tries = 0
    while len(base) < 400 and tries < 200_000:
        tries += 1
        n = rng.choice((1, 1, 2, 3, 7, 96, 255))
        sg = rng.choice((1.0, -1.0))
        t = math.ldexp(rng.uniform(1.0, 2.0), rng.randint(-1006, -990))
        c = (R_GAS * t) / (float(n) * F_CONST)
        if not c >= 2.0 ** -1021:
            continue
        ls = sg * float(lo / Fraction(c))
        if not abs(ls) < 700.0:
            continue
        q = ulps(math.exp(ls), -60)
        for _ in range(121):
            p = abs(Fraction(c) * Fraction(mpfr_ref.log(q)))
            if lo <= p < hi:
                base.append((n, t, q))
                break
            q = math.nextafter(q, INF)
    out = []
    for n, t, q in base:
        for k in (1, 3, 0, 2, -1, -3):
            out.append((fb(abs(k)) if k >= 0 else -fb(-k), n, t, q))
        out.append((rng.uniform(-3.0, 3.0), n, t, q))
    return out


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


# ── The overflow audit's B3 (ADR 0002 #31, #32): w·√((x/y)·2^e), the gas laws, the potentials, Born–Landé, BE/A ──────
# ids 58 … 70: wsqrt (_km_cr_mul_sqrt_ratio(w, x, y, e), e an integer word: the public helper, its fast path and flag;
# mode 1 _km_mul_sqrt_ratio_acc), igp, igv, dens, coul, force (ideal_gas_pressure, ideal_gas_volume, gas_density, coulomb,
# coulomb_force: the public function, its shape's straight-line fast path and flag, mode 1 the monomial's exact path at
# its factors), width (morse_width_parameter: _km_mul_sqrt_ratio_fast, mode 1 the exact root), vdw, prp, rkp, lj, lande,
# bea (van_der_waals_pressure, peng_robinson_pressure, redlich_kwong_pressure, lennard_jones, born_lande_lattice_energy,
# binding_energy_per_nucleon: the public function, the Rust-chain fast path where its range test passes and the flag,
# mode 1 the unbounded-exponent slow path forced). A charge or a count is an integer word. References: Fraction (the
# monomials), Python's integer square root (the roots, ties included), mpfr_ref.UX with MPFR's pow at 53 bits and the
# final operations (the #20-style ones). An Err (validation, or the co-volume test) is the NaN word in both slots.
K_CONST = _const("src/constants.cyr", "KM_COULOMB_CONSTANT")
NA_CONST = _const("src/constants.cyr", "KM_AVOGADRO")
QE_CONST = _const("src/constants.cyr", "KM_ELEMENTARY_CHARGE")
AMU_CONST = _const("src/constants.cyr", "KM_AMU_MEV")
MH_CONST = _const("src/constants.cyr", "KM_HYDROGEN_ATOM_MASS_U")
MN_CONST = _const("src/constants.cyr", "KM_NEUTRON_MASS_U")
AUDIT3_NAMES = ("wsqrt", "igp", "igv", "dens", "vdw", "prp", "rkp", "lj", "width", "coul", "force", "lande", "bea")


def sqrt_frac_ref(q):
    """RN(√q) for a positive Fraction q, exactly, ties to even (B3's wsqrt and width): Q = ⌊N·2^s/D⌋ of at least 120 bits
    and an even exponent, R = ⌊√Q⌋, RN((2R + sticky)·2^-(s/2 + 1)) by CPython's correctly rounded int division."""
    n, d = q.numerator, q.denominator
    s = 121 + d.bit_length() - n.bit_length()
    if s % 2:
        s += 1
    quo, rem = divmod(n << s, d) if s >= 0 else divmod(n, d << -s)
    r = math.isqrt(quo)
    st = 1 if (r * r != quo or rem) else 0
    return rn(Fraction(2 * r + st) * Fraction(2) ** (-(s // 2) - 1))


def _nn(x):
    return math.isfinite(x) and x >= 0


def ref_audit3(name, a):
    """The reference word of ids 58 … 70 (NANW for an Err or outside a helper's domain)."""
    UX = mpfr_ref.UX
    if name == "wsqrt":
        w, x, y, e = a
        if not (_pos(w) and _pos(x) and _pos(y)):
            return NANW
        return bits(sqrt_frac_ref(Fraction(w) ** 2 * Fraction(x) / Fraction(y) * Fraction(2) ** e))
    if name in ("igp", "igv"):
        n, t, v = a
        if not (math.isfinite(n) and _pos(t) and _pos(v)):
            return NANW
        return bits(mono_ref([n, R_GAS, t], [v]))
    if name == "dens":
        p, mm, t = a
        if not (_pos(t) and _nn(p) and _pos(mm)):
            return NANW
        return bits(mono_ref([p, mm], [R_GAS, t]))
    if name in ("coul", "force"):
        q1, q2, r = a
        if not _pos(r):
            return NANW
        den = [r] if name == "coul" else [r, r]
        if not (math.isfinite(q1) and math.isfinite(q2)):
            return bits(_ieee_div(K_CONST * q1 * q2, r if name == "coul" else r * r))
        return bits(mono_ref([K_CONST, q1, q2], den))
    if name == "width":
        w, mu, d = a
        if not (_pos(w) and _pos(mu) and _pos(d)):
            return NANW
        return bits(sqrt_frac_ref(Fraction(w) ** 2 * Fraction(mu) / (2 * Fraction(d))))
    with mpfr_ref.unbounded():
        if name == "vdw":
            n, t, v, aa, b = a
            if not (_pos(t) and math.isfinite(n) and _nn(aa) and _nn(b) and _pos(v)):
                return NANW
            nb = UX(n) * UX(b)
            if UX(v) <= nb:
                return NANW
            nv = UX(n) / UX(v)
            return bits(mpfr_ref.ux_sub_final(UX(n) * UX(R_GAS) * UX(t) / (UX(v) - nb), UX(aa) * nv * nv))
        if name == "prp":
            tc, pc, om, t, v = a
            if not (_pos(t) and _pos(tc) and _pos(pc) and math.isfinite(om) and math.isfinite(v)):
                return NANW
            r, tcu, pcu = UX(R_GAS), UX(tc), UX(pc)
            aa = UX(0.45724) * r * r * tcu * tcu / pcu
            b = UX(0.07780) * r * tcu / pcu
            if UX(v) <= b:
                return NANW
            o = UX(om)
            kappa = (UX(0.37464) + UX(1.54226) * o) - UX(0.26992) * o * o
            u = UX(1.0) + kappa * (UX(1.0) - (UX(t) / tcu).sqrt())
            vu = UX(v)
            den = (vu * vu + UX(2.0) * b * vu) - b * b
            return bits(mpfr_ref.ux_sub_final(r * UX(t) / (vu - b), aa * (u * u) / den))
        if name == "rkp":
            tc, pc, t, v = a
            if not (_pos(t) and _pos(tc) and _pos(pc) and math.isfinite(v)):
                return NANW
            r, tcu, pcu = UX(R_GAS), UX(tc), UX(pc)
            aa = UX(0.42748) * r * r * tcu.pow(2.5) / pcu
            b = UX(0.08664) * r * tcu / pcu
            if UX(v) <= b:
                return NANW
            vu = UX(v)
            return bits(mpfr_ref.ux_sub_final(r * UX(t) / (vu - b), aa / (UX(t).sqrt() * vu * (vu + b))))
        if name == "lj":
            e, s, r = a
            if not (_pos(e) and _pos(s) and _pos(r)):
                return NANW
            sr = UX(s) / UX(r)
            sr6 = sr * sr * sr * sr * sr * sr
            return bits(mpfr_ref.ux_mul_final(UX(4.0) * UX(e), sr6 * sr6 - sr6))
        if name == "lande":
            mm, zc, za, r0, n = a
            if not (_pos(r0) and math.isfinite(n) and n > 1.0 and _pos(mm) and 1 <= zc <= U32 and 1 <= za <= U32):
                return NANW
            u = UX(-K_CONST) * UX(NA_CONST) * UX(mm) * UX(float(zc)) * UX(float(za)) * UX(QE_CONST) * UX(QE_CONST)
            return bits(mpfr_ref.ux_mul_final(u / UX(r0), UX(1.0) - UX(1.0) / UX(n)))
        z, n, m = a                                  # bea
        if not (0 <= z <= U32 and 0 <= n <= U32 and z + n <= U32 and z + n != 0 and _pos(m)):
            return NANW
        md = UX(float(z)) * UX(MH_CONST) + UX(float(n)) * UX(MN_CONST) - UX(m)
        return bits(mpfr_ref.ux_div_final(md * UX(AMU_CONST), UX(float(z + n))))


def args_wsqrt(rng):
    """(w, x, y, e) for _km_cr_mul_sqrt_ratio: every binade of each, e over [−4, 4] and far, x/y past either end, the
    result past DBL_MAX and subnormal, exact roots, the midpoints 3w·2^k (w with 3w of 54 bits, x/y = 9·4^k) and the
    subnormal grid's ties q·2^-1075, roots within ulps of a midpoint, and the domain."""
    def e_(rng_):
        return rng_.choice((-1, -1, 0, 1, -2, 2, rng_.randint(-4, 4), rng_.randint(-200, 200)))
    rand = [(_pp(rng), _pp(rng), _pp(rng), e_(rng)) for _ in range(100_000)]
    typical = [(_pp(rng, 40, 60), _pp(rng, -95, -80), _pp(rng, -70, -55), -1) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        if kind == 0:
            edge.append((_pp(rng, -100, 100), _pp(rng, 512, 1023), _pp(rng, -1022, -512), e_(rng)))
        elif kind == 1:
            edge.append((_pp(rng, -100, 600), _pp(rng, -1022, -540), _pp(rng, 538, 1023), e_(rng)))
        elif kind == 2:
            edge.append((_pp(rng, -1074, -1023), _pp(rng), _pp(rng), e_(rng)))
        elif kind == 3:
            edge.append((_pp(rng), _pp(rng, -1074, -1023), _pp(rng, -1074, 1023), e_(rng)))
        else:
            edge.append((_pp(rng, 900, 1023), _pp(rng, 0, 1023), _pp(rng, -1074, 0), e_(rng)))
    exact, ties = [], []
    for _ in range(30_000):
        sq = float(rng.randrange(1, 1 << 26))
        k = rng.randint(-200, 200)
        y = _pp(rng, -400, 400)
        x = sq * sq * 2.0 ** (2 * k) * y
        if _pos(x) and x / y == sq * sq * 2.0 ** (2 * k):
            exact.append((_pp(rng, -300, 300), x, y, 0))
    for _ in range(40_000):
        if rng.random() < 0.6:
            k, j = rng.randint(-300, 300), rng.randint(-300, 300)
            ties.append((_y3(rng) * 2.0 ** rng.randint(-100, 100), 9.0 * 2.0 ** (2 * k), 2.0 ** (2 * j), -1 + 1))
        else:
            s = rng.randint(60, 120)
            q = float(rng.randrange(1, 1 << 41) | 1)
            ties.append((q * 2.0 ** (-1075 + s), 1.0, 2.0 ** (2 * s), 0))
    near = []
    for _ in range(30_000):
        z = 1.0 + rng.randrange(1 << 52) * 2.0 ** -52
        mm = Fraction(z) + Fraction(1, 1 << 53)
        w = 2.0 ** rng.randint(-200, 200)
        y = _pp(rng, -300, 300)
        x = rn(mm * mm * Fraction(y))
        near.append((w, ulps(x, rng.randint(-3, 3)), y, 0))
    dom = [(0.0, 1.0, 1.0, 0), (1.0, 0.0, 1.0, 0), (1.0, 1.0, 0.0, 0), (-1.0, 1.0, 1.0, 0), (INF, 1.0, 1.0, 0),
           (1.0, INF, 1.0, 0), (1.0, 1.0, INF, 0), (NAN, 1.0, 1.0, 0), (DMAX, DMAX, fb(1), 0), (fb(1), fb(1), DMAX, 0),
           (DMAX, 4.0, 1.0, 0), (1.0, 4.0, 1.0, 0)]
    return [("binade-uniform", rand), ("typical (HCl's Morse width)", typical), ("x/y and the result past the range",
            edge), ("exact roots", exact), ("midpoints and subnormal ties", ties), ("within ulps of a midpoint", near),
            ("domain", dom)]


def args_ideal(rng):
    """(n, T, x) for ideal_gas_pressure and ideal_gas_volume."""
    rand = [(_spp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_spp(rng, -10, 10), rng.uniform(1.0, 3000.0), _pp(rng, -20, 30)) for _ in range(60_000)]
    comp = []
    for _ in range(40_000):
        if rng.random() < 0.5:
            n, t = _spp(rng, 520, 1023), _pp(rng, 400, 1023)
        else:
            n, t = _spp(rng, -1074, -520), _pp(rng, -1074, -400)
        e = math.frexp(n)[1] + math.frexp(t)[1]
        comp.append((n, t, _around(rng, e)))
    ties = []
    for _ in range(40_000):
        if rng.random() < 0.6:
            q = rng.choice((9.0, 11.0, 13.0, 15.0)) * rng.choice((1.0, -1.0))
            ties.append((q * 2.0 ** rng.randint(-900, 900), 2.0 ** rng.randint(-60, 60), 2.0 ** rng.randint(-60, 60)))
        else:
            qq, d = _sub_tie(rng)
            ties.append((qq, 1.0, R_GAS * d))
    bad = [(NAN, 1.0, 1.0), (INF, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, -1.0, 1.0), (1.0, INF, 1.0), (1.0, 1.0, 0.0),
           (1.0, 1.0, NAN), (0.0, DMAX, fb(1)), (-0.0, 1e300, 1e-300), (DMAX, DMAX, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("n·R·T past the range", comp),
            ("midpoints q·R and subnormal ties", ties), ("Err and zeros", bad)]


def args_dens(rng):
    """(P, M, T) for gas_density."""
    rand = [(_pp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(0.0, 1e8), rng.uniform(0.002, 0.4), rng.uniform(1.0, 3000.0)) for _ in range(60_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 3)
        if kind == 0:                               # R·T past DBL_MAX
            edge.append((_pp(rng, 0, 1023), _pp(rng, -100, 1023), _pp(rng, 1020, 1023)))
        elif kind == 1:                             # R·T subnormal
            edge.append((_pp(rng, -1074, -900), _pp(rng, -100, 100), _pp(rng, -1074, -1020)))
        elif kind == 2:                             # P·M past either end
            edge.append((_pp(rng, 520, 1023), _pp(rng, 520, 1023), _pp(rng, 0, 1023)))
        else:
            edge.append((_pp(rng, -1074, -520), _pp(rng, -1074, -400), _pp(rng, -1074, 0)))
    ties = []
    for _ in range(20_000):
        qq, d = _sub_tie(rng)
        b = rng.randint(-30, 30)
        ties.append((qq, R_GAS * 2.0 ** b, d * 2.0 ** b))
    bad = [(1.0, 1.0, 0.0), (-1.0, 1.0, 1.0), (1.0, 0.0, 1.0), (NAN, 1.0, 1.0), (1.0, 1.0, INF), (-0.0, 1.0, 1.0),
           (0.0, DMAX, fb(1)), (101325.0, 0.029, DMAX)]
    return [("binade-uniform", rand), ("typical", typical), ("R·T or P·M past the range", edge), ("subnormal ties", ties),
            ("Err and zeros", bad)]


def args_coul(rng):
    """(q1, q2, r) for coulomb and coulomb_force."""
    rand = [(_spp(rng), _spp(rng), _pp(rng)) for _ in range(80_000)]
    typical = [(_spp(rng, -70, -55), _spp(rng, -70, -55), _pp(rng, -40, -20)) for _ in range(60_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 3)
        if kind == 0:
            edge.append((_spp(rng, 500, 1023), _spp(rng, 0, 1023), _pp(rng, 300, 1023)))
        elif kind == 1:
            edge.append((_spp(rng, -1074, -500), _spp(rng, -600, -60), _pp(rng, -1074, -300)))
        elif kind == 2:                             # r·r past DBL_MAX
            edge.append((_spp(rng, 300, 700), _spp(rng, 300, 700), _pp(rng, 512, 900)))
        else:                                       # r·r below 2^-1022
            edge.append((_spp(rng, -700, -300), _spp(rng, -400, -60), _pp(rng, -1074, -512)))
    ties = []
    for _ in range(30_000):
        if rng.random() < 0.6:
            q = rng.choice((5.0, 7.0)) * rng.choice((1.0, -1.0))
            ties.append((q * 2.0 ** rng.randint(-600, 600), 2.0 ** rng.randint(-60, 60), 2.0 ** rng.randint(-60, 60)))
        else:
            qq, d = _sub_tie(rng)
            ties.append((qq, 1.0, K_CONST * d))
    sp = [0.0, -0.0, INF, -INF, NAN, 1.0, fb(1), DMAX, 1e300, -1e-300]
    specials = [(a, b, c) for a in sp for b in sp for c in (1.0, fb(1), DMAX, 0.0, -1.0, INF, NAN)]
    return [("binade-uniform", rand), ("typical", typical), ("k·q₁·q₂ or r·r past the range", edge),
            ("midpoints q·k_e and subnormal ties", ties), ("specials and Err", specials)]


def args_width(rng):
    """(ω, μ, D_e) for morse_width_parameter."""
    secs = []
    for tag, sec in args_wsqrt(rng):
        out = []
        for w, x, y, e in sec:
            if e == -1:
                out.append((w, x, y))
            elif e == 0 and math.isfinite(y) and y > 0 and y * 0.5 > 0 and (y * 0.5) * 2.0 == y:
                out.append((w, x, y * 0.5))         # ω·√(x/(2·y/2)) = ω·√(x/y)
            else:
                out.append((w, x, y))
        secs.append((tag, out))
    return secs


def args_vdw(rng):
    """(n, T, V, a, b) for van_der_waals_pressure."""
    rand = [(_spp(rng), _pp(rng), _pp(rng), _pp(rng), _pp(rng)) for _ in range(60_000)]
    typical = []
    for _ in range(60_000):
        n, b = _pp(rng, -10, 10), rng.uniform(1e-5, 2e-4)
        typical.append((n, rng.uniform(100.0, 2000.0), n * b * (1.0 + 10.0 ** rng.uniform(-3, 4)), rng.uniform(0, 4), b))
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        if kind == 0:                               # n·R·T past DBL_MAX
            edge.append((_spp(rng, 990, 1023), _pp(rng, 0, 40), _pp(rng, 990, 1023), rng.random(), _pp(rng, -40, -30)))
        elif kind == 1:                             # n·R·T below 2^-1022, V tiny
            edge.append((_spp(rng, -1074, -960), _pp(rng, -100, 10), _pp(rng, -1074, -960), _pp(rng, -1074, 0), 0.0))
        elif kind == 2:                             # n/V and a·(n/V)² past DBL_MAX
            edge.append((_pp(rng, 800, 1023), _pp(rng, -100, -30), _pp(rng, -330, -130), rng.choice((0.0, _pp(rng, -1074, -1000))), 0.0))
        elif kind == 3:                             # V − n·b past DBL_MAX
            edge.append((-_pp(rng, 960, 1023), _pp(rng, 0, 10), _pp(rng, 930, 1023), rng.random(), _pp(rng, 1, 100)))
        else:                                       # a subnormal n·b against a subnormal V
            n = fb(rng.randrange(1, 1 << 20))
            b = rng.random()
            k = round(n * b / fb(1)) + rng.randint(-3, 3)
            edge.append((n, 300.0, fb(max(1, k)), 0.0, b))
    bad = [(1.0, 0.0, 1.0, 0.0, 0.0), (NAN, 1.0, 1.0, 0.0, 0.0), (1.0, 1.0, 1.0, -1.0, 0.0), (1.0, 1.0, 1.0, 0.0, -1.0),
           (1.0, 1.0, 0.0, 0.0, 0.0), (1.0, 300.0, 0.01, 0.0, 0.01), (1.0, 300.0, 0.001, 0.0, 0.01),
           (DMAX, DMAX, DMAX, 0.0, 0.0), (fb(3), 300.0, fb(1), 0.0, 0.2)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range, the co-volume at the grid", edge),
            ("Err and specials", bad)]


def args_prp(rng):
    """(Tc, Pc, ω, T, V) for peng_robinson_pressure."""
    def pr_b(tc, pc, c=0.07780):
        b = c * R_GAS * tc / pc
        return b if math.isfinite(b) and b > 0 else None
    rand = []
    for _ in range(60_000):
        tc, pc = _pp(rng), _pp(rng)
        b = pr_b(tc, pc)
        v = b * (1.0 + 10.0 ** rng.uniform(-3, 8)) if b else _pp(rng)
        rand.append((tc, pc, rng.choice((rng.uniform(-0.5, 1.0), _spp(rng))), _pp(rng), v))
    typical = []
    for _ in range(60_000):
        tc, pc = rng.uniform(5.0, 700.0), rng.uniform(1e5, 2.5e7)
        b = pr_b(tc, pc)
        typical.append((tc, pc, rng.uniform(-0.3, 1.0), tc * rng.uniform(0.3, 10.0), b * (1.0 + 10.0 ** rng.uniform(-3, 5))))
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        if kind == 0:                               # a subnormal or tiny Tc (T/Tc past DBL_MAX)
            tc, pc = _pp(rng, -1074, -1000), _pp(rng, -100, 30)
            b = pr_b(tc, pc)
            edge.append((tc, pc, rng.uniform(-0.5, 1.0), _pp(rng, -10, 20), (b or fb(1)) * (1.0 + 10.0 ** rng.uniform(-3, 8))))
        elif kind == 1:                             # V·V past DBL_MAX
            edge.append((_pp(rng, 400, 600), _pp(rng, -10, 30), rng.uniform(-0.5, 1.0), _pp(rng, 0, 20), _pp(rng, 520, 1000)))
        elif kind == 2:                             # ω² past DBL_MAX, T = Tc sometimes (κ·0)
            tc = rng.uniform(5.0, 700.0)
            edge.append((tc, 3.4e6, _spp(rng, 510, 1023), rng.choice((tc, _pp(rng, -10, 20))), 0.025))
        elif kind == 3:                             # R·T past DBL_MAX
            edge.append((126.2, 3.394e6, 0.0377, _pp(rng, 1020, 1023), _pp(rng, 990, 1023)))
        else:                                       # a subnormal Tc, a tiny Pc: V between Rust's b and the exact one
            tc, pc = fb(rng.randrange(1, 300)), 2.0 ** -rng.randint(30, 60)
            b = pr_b(tc, pc)
            v = ulps(b, rng.randint(-40, 40)) if b else fb(1)
            edge.append((tc, pc, 0.0377, 300.0, v))
    bad = [(126.2, 3.394e6, 0.0377, 0.0, 0.025), (0.0, 3.394e6, 0.0377, 300.0, 0.025), (126.2, 0.0, 0.0377, 300.0, 0.025),
           (126.2, 3.394e6, NAN, 300.0, 0.025), (126.2, 3.394e6, 0.0377, 300.0, INF), (126.2, 3.394e6, 0.0377, 300.0, 1e-6),
           (fb(1), 3394000.0, 0.0377, 300.0, 0.025)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range, the co-volume", edge),
            ("Err and specials", bad)]


def args_rkp(rng):
    """(Tc, Pc, T, V) for redlich_kwong_pressure."""
    def rk_b(tc, pc):
        b = 0.08664 * R_GAS * tc / pc
        return b if math.isfinite(b) and b > 0 else None
    rand = []
    for _ in range(60_000):
        tc, pc = _pp(rng), _pp(rng)
        b = rk_b(tc, pc)
        rand.append((tc, pc, _pp(rng), b * (1.0 + 10.0 ** rng.uniform(-3, 8)) if b else _pp(rng)))
    typical = []
    for _ in range(60_000):
        tc, pc = rng.uniform(5.0, 700.0), rng.uniform(1e5, 2.5e7)
        typical.append((tc, pc, tc * rng.uniform(0.3, 10.0), rk_b(tc, pc) * (1.0 + 10.0 ** rng.uniform(-3, 5))))
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 3)
        if kind == 0:                               # Tc^2.5 past DBL_MAX, a finite
            e = rng.randint(410, 800)
            tc = _pp(rng, e, e)
            pc = _pp(rng, min(1023, int(2.5 * e) - 1020), 1023)
            b = rk_b(tc, pc) or 1.0
            edge.append((tc, pc, _pp(rng, -10, 12), b * 2.0 ** rng.randint(1, 600)))
        elif kind == 1:                             # Tc^2.5 below 2^-1022, a normal
            e = rng.randint(-1074, -410)
            tc = _pp(rng, e, e)
            pc = _pp(rng, max(-1074, int(2.5 * e) + 20), max(-1074, int(2.5 * e) + 1000))
            b = rk_b(tc, pc)
            edge.append((tc, pc, _pp(rng, -10, 12), b * 2.0 ** rng.randint(1, 600) if b else _pp(rng)))
        elif kind == 2:                             # Tc^2.5's ties at the rescale: Tc = t²·4^j, t⁵ of 54 bits
            t = (rng.randint(776, 890) * 2 + 1) / 1024.0
            j = rng.choice((rng.randint(-500, -205), rng.randint(205, 500)))
            tc = t * t * 2.0 ** (2 * j)
            pc = 2.0 ** (-1000 if j < 0 else 1000)
            b = rk_b(tc, pc) or 1.0
            edge.append((tc, pc, 300.0, b * 4.0))
        else:                                       # a subnormal Tc, a tiny Pc: V around Rust's b
            tc, pc = fb(rng.randrange(1, 300)), 2.0 ** -rng.randint(30, 60)
            b = rk_b(tc, pc)
            edge.append((tc, pc, 300.0, ulps(b, rng.randint(-40, 40)) if b else fb(1)))
    bad = [(126.2, 3.394e6, 0.0, 0.025), (0.0, 3.394e6, 300.0, 0.025), (126.2, 0.0, 300.0, 0.025),
           (126.2, 3.394e6, 300.0, NAN), (126.2, 3.394e6, 300.0, 1e-6), (1e130, 1e300, 300.0, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("Tc^2.5 past the range, its ties, the co-volume", edge),
            ("Err and specials", bad)]


def args_lj(rng):
    """(ε, σ, r) for lennard_jones."""
    rand = [(_pp(rng), _pp(rng), _pp(rng)) for _ in range(80_000)]
    typical = []
    for _ in range(60_000):
        s = _pp(rng, -40, -25)
        typical.append((_pp(rng, -80, -60), s, s * rng.uniform(0.5, 100.0)))
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 4)
        s = _pp(rng, -600, 600)
        if kind == 0:                               # (σ/r)¹² past DBL_MAX, ε tiny
            edge.append((_pp(rng, -1074, -600), s, s / _pp(rng, 85, 170)))
        elif kind == 1:                             # σ/r past 2^170.7: Rust's ∞ − ∞
            edge.append((_pp(rng, -1074, 0), s, s / _pp(rng, 171, 300)))
        elif kind == 2:                             # (σ/r)⁶ below 2^-1022, ε huge
            edge.append((_pp(rng, 900, 1023), s, s / _pp(rng, -200, -170)))
        elif kind == 3:                             # ε ≥ 2^1022 (#17's) with σ/r anywhere
            edge.append((2.0 ** 1022 * rng.uniform(1.0, 4.0) if rng.random() < 0.9 else DMAX, s, s / _pp(rng, -200, 200)))
        else:                                       # the final product's subnormal ties: 4ε·2^(6k) = q·2^-1075
            k = rng.randint(-360, -330)
            q = rng.randrange(1, 1 << 7) | 1
            m = -1075 - 6 * k
            eps = q * 2.0 ** (m - 2) if m - 2 < 1018 else q * 2.0 ** 1010
            edge.append((eps, 2.0 ** k, 1.0))
    bad = [(0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0), (NAN, 1.0, 1.0), (1.0, INF, 1.0), (1.0, 1.0, -1.0),
           (DMAX, 1.0, 1.0), (fb(1), fb(1), 1.0), (1.0, 1.0, fb(1))]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range, #17's region, ties", edge),
            ("Err and specials", bad)]


def args_lande(rng):
    """(M, z⁺, z⁻, r₀, n) for born_lande_lattice_energy."""
    def z(rng_):
        return rng_.choice((1, 2, 3, 4, rng_.randint(1, 1000), rng_.randint(1, U32), U32))
    rand = [(_pp(rng), z(rng), z(rng), _pp(rng), 1.0 + _pp(rng)) for _ in range(80_000)]
    typical = [(rng.uniform(1.6, 2.6), rng.randint(1, 3), rng.randint(1, 3), rng.uniform(1.5e-10, 4e-10),
                rng.uniform(5.0, 12.0)) for _ in range(60_000)]
    edge = []
    for _ in range(60_000):
        kind = rng.randint(0, 2)
        if kind == 0:                               # the prefix past DBL_MAX, e² and r₀ bringing it back
            edge.append((_pp(rng, 900, 1023), z(rng), z(rng), _pp(rng, 800, 1023), rng.uniform(5.0, 12.0)))
        elif kind == 1:                             # ·e·e below 2^-1022, r₀ tiny
            edge.append((_pp(rng, -1074, -1000), z(rng), z(rng), _pp(rng, -1074, -960), rng.uniform(5.0, 12.0)))
        else:                                       # n past 2^1021 or near 1
            edge.append((_pp(rng), z(rng), z(rng), _pp(rng), rng.choice((_pp(rng, 1000, 1023), 1.0 + 2.0 ** -rng.randint(1, 52)))))
    bad = [(1.0, 1, 1, 0.0, 8.0), (1.0, 1, 1, 1.0, 1.0), (1.0, 1, 1, 1.0, NAN), (0.0, 1, 1, 1.0, 8.0), (1.0, 0, 1, 1.0, 8.0),
           (1.0, 1, 1 << 32, 1.0, 8.0), (DMAX, 1, 1, DMAX, DMAX), (fb(1), 1, 1, 2.82e-10, 8.0)]
    return [("binade-uniform", rand), ("typical", typical), ("the chain past the range, n's ends", edge),
            ("Err and specials", bad)]


def args_bea(rng):
    """(Z, N, M) for binding_energy_per_nucleon."""
    rand = [(rng.randint(0, U32), rng.randint(0, U32), _pp(rng)) for _ in range(60_000)]
    typical = []
    for _ in range(60_000):
        zz, nn = rng.randint(1, 118), rng.randint(0, 180)
        typical.append((zz, nn, (zz + nn) * rng.uniform(0.99, 1.01)))
    edge = []
    for _ in range(60_000):
        zz = rng.randint(0, 1 << 31)
        nn = rng.randint(1, 1 << 31)
        if rng.random() < 0.3:
            zz, nn = rng.randint(0, 600), rng.randint(1, 600)
        edge.append((zz, nn, _pp(rng, 1013, 1023)))
    bad = [(0, 0, 1.0), (U32, 1, 1.0), (1 << 31, 1 << 31, 1.0), (26, 30, 0.0), (26, 30, NAN), (26, 30, INF),
           (1000, 1000, DMAX), (465, 465, DMAX), (466, 466, DMAX)]
    return [("binade-uniform", rand), ("typical", typical), ("Δm·931.494 past DBL_MAX", edge), ("Err and specials", bad)]


# ── The overflow audit's B4 (ADR 0002 #33, #34): the sums with cancellation and the Shomate polynomials ─────────────────
# ids 71 … 82: hess (hess_law: n ≤ 5, then five value words), eff (enthalpy_from_formation: np, nr ≤ 3, three product and
# three reactant words), rsum (thermochem's reaction sums: col 0 reaction_enthalpy, 1 reaction_gibbs_energy,
# 2 reaction_entropy; np, nr ≤ 2; then two (formula, coefficient) pairs a side, a formula as its THERMOCHEM_DATA row,
# 76 or more an unknown formula: Rust's Err), qv (q_value: as eff), ionic (ionic_strength: n ≤ 3, then (c, z) pairs, z an
# integer word), haber (born_haber_lattice_energy), gibbs (gibbs_free_energy), kohl (kohlrausch_molar_conductivity), lim
# (limiting_molar_conductivity), zoc (zero_order_concentration), shcp and shh (Shomate::cp and enthalpy_antiderivative:
# A … E, T): mode 0 the public function, its Rust-chain fast path where the test passes and the flag; mode 1 the
# unbounded-exponent slow path forced (the _km_*_ux functions). Reference: mpfr_ref.UX with mpfr_ref.ux_sum in Rust's
# order (from −0.0, thermochem's `+=` from +0.0) and the final operation rounded once; Rust's own expression for a
# non-finite input, and the Shomate polynomials' at T = ±0 (exact in both). An Err is the NaN word in both slots.
AUDIT4_NAMES = ("hess", "eff", "rsum", "qv", "ionic", "haber", "gibbs", "kohl", "lim", "zoc", "shcp", "shh")
AMU4 = _const("src/constants.cyr", "KM_AMU_MEV")


def _tc_column(name):
    """A THERMOCHEM_DATA f64 column of src/thermochem_data.cyr (generated), in row order."""
    import re
    text = (_ROOT / "src/thermochem_data.cyr").read_text()
    body = text[text.index(f"var {name}: f64[76] = {{"):]
    body = body[:body.index("};")]
    return [float(v) for v in re.findall(r"^\s+(-?[0-9.]+),", body, re.M)]


TC_COLS = [_tc_column(n) for n in ("KM_THERMOCHEM_DELTA_HF_KJ", "KM_THERMOCHEM_DELTA_GF_KJ", "KM_THERMOCHEM_S_STANDARD_J")]
assert all(len(c) == 76 for c in TC_COLS)


def _fin(*xs):
    return all(math.isfinite(x) for x in xs)


def _rsum(xs, start=-0.0):
    s = start
    for x in xs:
        s = s + x
    return s


def _shomate_rust(is_cp, a, b, c, d, e, temp):
    t = temp / 1000.0
    if is_cp:
        return a + b * t + c * t * t + d * t * t * t + _ieee_div(e, t * t)
    return 1000.0 * (t * (a + t * (b / 2.0 + t * (c / 3.0 + t * (d / 4.0)))) - _ieee_div(e, t))


def ref_audit4(name, a):
    """The reference word of ids 71 … 82 (NANW for an Err)."""
    UX = mpfr_ref.UX
    with mpfr_ref.unbounded():
        if name == "hess":
            xs = list(a[1:1 + a[0]])
            return bits(mpfr_ref.ux_sum_final(xs) if _fin(*xs) else _rsum(xs))
        if name in ("eff", "qv"):
            p, r = list(a[2:2 + a[0]]), list(a[5:5 + a[1]])
            if not _fin(*p, *r):
                d = _rsum(p) - _rsum(r)
                return bits(d if name == "eff" else d * AMU4)
            if name == "eff":
                return bits(mpfr_ref.ux_sub_final(mpfr_ref.ux_sum(p), mpfr_ref.ux_sum(r)))
            return bits(mpfr_ref.ux_mul_final(mpfr_ref.ux_sum(p) - mpfr_ref.ux_sum(r), UX(AMU4)))
        if name == "rsum":
            col, np_, nr = a[0], a[1], a[2]
            prod = [(a[3 + 2 * j], a[4 + 2 * j]) for j in range(np_)]
            reac = [(a[7 + 2 * j], a[8 + 2 * j]) for j in range(nr)]
            if any(f >= 76 for f, _ in prod + reac):
                return NANW
            v = TC_COLS[col]
            if not _fin(*[c for _, c in prod + reac]):
                d = _rsum([c * v[f] for f, c in prod], 0.0) - _rsum([c * v[f] for f, c in reac], 0.0)
                return bits(d * 1000.0 if col != 2 else d)
            sp = mpfr_ref.ux_sum([UX(c) * UX(v[f]) for f, c in prod], 0.0)
            sr = mpfr_ref.ux_sum([UX(c) * UX(v[f]) for f, c in reac], 0.0)
            return bits(mpfr_ref.ux_sub_final(sp, sr) if col == 2 else mpfr_ref.ux_mul_final(sp - sr, UX(1000.0)))
        if name == "ionic":
            ions = [(a[1 + 2 * j], a[2 + 2 * j]) for j in range(a[0])]
            if not _fin(*[c for c, _ in ions]):
                return bits(0.5 * _rsum([c * float(z) * float(z) for c, z in ions]))
            terms = [UX(c) * UX(float(z)) * UX(float(z)) for c, z in ions]
            return bits(mpfr_ref.ux_mul_final(UX(0.5), mpfr_ref.ux_sum(terms)))
        if name == "haber":
            f_, s, ie, d, ea = a
            if not _fin(*a):
                return bits(f_ - s - ie - d / 2.0 - ea)
            return bits(mpfr_ref.ux_sub_final(UX(f_) - UX(s) - UX(ie) - UX(d) / UX(2.0), UX(ea)))
        if name == "gibbs":
            h, t, s = a
            if not (_pos(t) and _fin(h, s)):
                return NANW
            return bits(mpfr_ref.ux_sub_final(UX(h), UX(t) * UX(s)))
        if name == "kohl":
            l0, k, c = a
            if not (_nonneg(c) and _fin(l0, k)):
                return NANW
            return bits(mpfr_ref.ux_sub_final(UX(l0), UX(k) * UX(c).sqrt()))
        if name == "lim":
            cat, nc, an, na = a
            if not _fin(*a):
                return bits(nc * cat + na * an)
            return bits(mpfr_ref.ux_add_final(UX(nc) * UX(cat), UX(na) * UX(an)))
        if name == "zoc":
            c0, k, t = a
            if not _fin(*a):
                rem = c0 - k * t
                return bits(0.0 if rem < 0.0 else rem)
            p = UX(k) * UX(t)
            if (UX(c0) - p).cmp(UX(0.0)) < 0:
                return bits(0.0)
            return bits(mpfr_ref.ux_sub_final(UX(c0), p))
        aa, bb, cc, dd, ee, temp = a                    # shcp, shh
        if not _fin(*a) or temp == 0.0:
            return bits(_shomate_rust(name == "shcp", aa, bb, cc, dd, ee, temp))
        t = UX(temp) / UX(1000.0)
        if name == "shcp":
            return bits(mpfr_ref.ux_add_final(UX(aa) + UX(bb) * t + UX(cc) * t * t + UX(dd) * t * t * t,
                                              UX(ee) / (t * t)))
        inner = t * (UX(aa) + t * (UX(bb) / UX(2.0) + t * (UX(cc) / UX(3.0) + t * (UX(dd) / UX(4.0)))))
        return bits(mpfr_ref.ux_mul_final(UX(1000.0), inner - UX(ee) / t))


def _bigs(rng, lo=0.3, hi=1.0):
    return DMAX * rng.uniform(lo, hi) * rng.choice((1.0, -1.0))


def _slice_args(rng, n_max, width, gen):
    """(n, x1 … x_width) with n ≤ n_max values from gen, the unused words 0."""
    n = rng.randint(0, n_max)
    xs = [gen() for _ in range(n)]
    return (n, *xs, *([0.0] * (width - n)))


def args_hess(rng):
    rand = [_slice_args(rng, 5, 5, lambda: _spp(rng)) for _ in range(40_000)]
    typical = [_slice_args(rng, 5, 5, lambda: rng.uniform(-1e6, 1e6)) for _ in range(20_000)]
    edge = [_slice_args(rng, 5, 5, lambda: _bigs(rng)) for _ in range(40_000)]
    sub = [_slice_args(rng, 5, 5, lambda: _spp(rng, -1074, -1015)) for _ in range(10_000)]
    bad = [(2, DMAX, DMAX, 0.0, 0.0, 0.0), (3, DMAX, DMAX, -DMAX, 0.0, 0.0), (2, INF, -INF, 0.0, 0.0, 0.0),
           (1, NAN, 0.0, 0.0, 0.0, 0.0), (0, 0.0, 0.0, 0.0, 0.0, 0.0), (2, -0.0, -0.0, 0.0, 0.0, 0.0),
           (3, DMAX, 2.0 ** 970, -DMAX, 0.0, 0.0), (2, DMAX, 2.0 ** 970, 0.0, 0.0, 0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("partial sums past DBL_MAX", edge), ("subnormal", sub),
            ("specials", bad)]


def _two_slices(rng, gen):
    np_, nr = rng.randint(0, 3), rng.randint(0, 3)
    p = [gen() for _ in range(np_)] + [0.0] * (3 - np_)
    r = [gen() for _ in range(nr)] + [0.0] * (3 - nr)
    return (np_, nr, *p, *r)


def args_eff(rng):
    rand = [_two_slices(rng, lambda: _spp(rng)) for _ in range(40_000)]
    typical = [_two_slices(rng, lambda: rng.uniform(-5e6, 5e6)) for _ in range(20_000)]
    edge = [_two_slices(rng, lambda: _bigs(rng)) for _ in range(40_000)]
    sub = [_two_slices(rng, lambda: _spp(rng, -1074, -1015)) for _ in range(10_000)]
    bad = [(2, 2, DMAX, DMAX, 0.0, DMAX, DMAX, 0.0), (1, 1, INF, 0.0, 0.0, INF, 0.0, 0.0), (0, 0, 0, 0, 0, 0, 0, 0),
           (1, 0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0), (1, 1, NAN, 0.0, 0.0, 1.0, 0.0, 0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("sums past DBL_MAX", edge), ("subnormal", sub),
            ("specials", bad)]


def args_qv(rng):
    rand = [_two_slices(rng, lambda: _pp(rng)) for _ in range(40_000)]
    typical = [_two_slices(rng, lambda: rng.uniform(1.0, 250.0)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):
        x = _two_slices(rng, lambda: DMAX * rng.uniform(0.3, 1.0))
        edge.append(x)
    sub = [_two_slices(rng, lambda: _pp(rng, -1074, -1015)) for _ in range(10_000)]
    bad = [(2, 2, DMAX, DMAX, 0.0, DMAX, 0.999 * DMAX, 0.0), (2, 2, DMAX, DMAX, 0.0, DMAX, DMAX, 0.0),
           (1, 1, INF, 0.0, 0.0, INF, 0.0, 0.0), (0, 0, 0, 0, 0, 0, 0, 0), (1, 0, NAN, 0.0, 0.0, 0.0, 0.0, 0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("sums past DBL_MAX", edge), ("subnormal", sub),
            ("specials", bad)]


def args_rsum(rng):
    def side(gen):
        n = rng.randint(0, 2)
        out = []
        for j in range(2):
            out += [rng.randint(0, 75), gen()] if j < n else [0, 0.0]
        return n, out
    def mk(gen):
        np_, p = side(gen)
        nr, r = side(gen)
        return (rng.randint(0, 2), np_, nr, *p, *r)
    rand = [mk(lambda: _spp(rng)) for _ in range(40_000)]
    typical = [mk(lambda: rng.choice((1.0, 2.0, 0.5, 3.0)) * rng.uniform(0.5, 2.0)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):                         # products past DBL_MAX: both sides close, the same formula
        col, f = rng.randint(0, 2), rng.randint(0, 75)
        v = abs(TC_COLS[col][f]) or 1.0
        c = DMAX / max(v, 1.0) * rng.uniform(0.5, 1.0)
        edge.append((col, 2, 2, f, c, f, c, f, c * rng.uniform(0.999, 1.001), f, c))
    sub = [mk(lambda: _spp(rng, -1074, -1010)) for _ in range(10_000)]
    bad = [(0, 1, 0, 76, 1.0, 0, 0.0, 0, 0.0, 0, 0.0), (2, 1, 1, 3, 1.0, 0, 0.0, 80, 1.0, 0, 0.0),
           (0, 1, 1, 3, INF, 0, 0.0, 3, INF, 0, 0.0), (1, 1, 0, 0, NAN, 0, 0.0, 0, 0.0, 0, 0.0),
           (0, 0, 0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("products past DBL_MAX, cancelling", edge),
            ("subnormal products", sub), ("Err and specials", bad)]


def args_ionic(rng):
    def mk(gen):
        n = rng.randint(0, 3)
        out = []
        for j in range(3):
            if j < n:
                z = rng.choice((1, -1, 2, -2, 3, rng.randint(-1000, 1000), rng.choice((2 ** 31 - 1, -2 ** 31))))
                out += [gen(z), z]
            else:
                out += [0.0, 0]
        return (n, *out)
    rand = [mk(lambda z: _spp(rng)) for _ in range(40_000)]
    typical = [mk(lambda z: rng.uniform(0.0, 2.0)) for _ in range(20_000)]
    edge = [mk(lambda z: _bigs(rng, 0.05, 1.0) / max(1.0, float(z) * float(z))) for _ in range(40_000)]
    bad = [(1, 0.4 * DMAX, 2, 0.0, 0, 0.0, 0), (2, DMAX, 2, -DMAX, 2, 0.0, 0), (1, INF, 1, 0.0, 0, 0.0, 0),
           (1, NAN, 1, 0.0, 0, 0.0, 0), (0, 0.0, 0, 0.0, 0, 0.0, 0), (1, fb(3), 3, 0.0, 0, 0.0, 0)]
    return [("binade-uniform", rand), ("typical", typical), ("terms past DBL_MAX", edge), ("specials", bad)]


def args_haber(rng):
    rand = [tuple(_spp(rng) for _ in range(5)) for _ in range(40_000)]
    typical = [(rng.uniform(-1000, 0), rng.uniform(50, 200), rng.uniform(300, 1000), rng.uniform(100, 500),
                rng.uniform(-400, 0)) for _ in range(20_000)]
    edge = [tuple(_bigs(rng) for _ in range(5)) for _ in range(40_000)]
    sub = [(fb(rng.randrange(0, 1 << 20)), 0.0, 0.0, fb(rng.randrange(0, 1 << 20) | 1), -fb(rng.randrange(0, 8)))
           for _ in range(10_000)]
    bad = [(DMAX, -DMAX, DMAX, 0.0, 0.0), (0.0, 0.0, 0.0, fb(1), 0.0), (INF, 0.0, 0.0, 0.0, 0.0), (NAN, 1, 1, 1, 1),
           (1.0, 1.0, 1.0, INF, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("differences past DBL_MAX", edge),
            ("an odd subnormal ΔH_diss", sub), ("specials", bad)]


def args_gibbs(rng):
    rand = [(_spp(rng), _pp(rng), _spp(rng)) for _ in range(40_000)]
    typical = [(rng.uniform(-1e6, 1e6), rng.uniform(1.0, 5000.0), rng.uniform(-500.0, 500.0)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 1)
        if kind == 0:                               # T·ΔS in (DBL_MAX, 2·DBL_MAX), ΔH of the same sign
            t = _pp(rng, 1, 60)
            sg = rng.choice((1.0, -1.0))
            edge.append((sg * DMAX * rng.uniform(0.2, 1.0), t, sg * (DMAX / t) * rng.uniform(1.0, 1.99)))
        else:                                       # T·ΔS subnormal, ΔH tiny or 0
            t = _pp(rng, -10, 20)
            edge.append((rng.choice((0.0, _spp(rng, -1074, -1000))), t, _spp(rng, -1074, -1015) / t))
    bad = [(1.0, 0.0, 1.0), (1.0, -1.0, 1.0), (NAN, 1.0, 1.0), (1.0, 1.0, INF), (1.0, INF, 1.0), (DMAX, DMAX, 1.0)]
    return [("binade-uniform", rand), ("typical", typical), ("T·ΔS past the range", edge), ("Err and specials", bad)]


def args_kohl(rng):
    rand = [(_spp(rng), _spp(rng), _pp(rng)) for _ in range(40_000)]
    typical = [(rng.uniform(50, 450), rng.uniform(0, 200), rng.uniform(0, 1)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 1)
        if kind == 0:
            c = _pp(rng, 2, 60)
            sg = rng.choice((1.0, -1.0))
            edge.append((sg * DMAX * rng.uniform(0.2, 1.0), sg * (DMAX / math.sqrt(c)) * rng.uniform(1.0, 1.99), c))
        else:
            c = _pp(rng, -40, 40)
            edge.append((rng.choice((0.0, _spp(rng, -1074, -1000))), _spp(rng, -1074, -1015) / math.sqrt(c), c))
    bad = [(1.0, 1.0, -1.0), (1.0, 1.0, NAN), (NAN, 1.0, 1.0), (1.0, INF, 1.0), (1.0, 1.0, -0.0), (1.0, 0.0, 0.0)]
    return [("binade-uniform", rand), ("typical", typical), ("K·√c past the range", edge), ("Err and specials", bad)]


def args_lim(rng):
    rand = [tuple(_spp(rng) for _ in range(4)) for _ in range(40_000)]
    typical = [(rng.uniform(10, 350), rng.choice((1.0, 2.0, 3.0)), rng.uniform(10, 200), rng.choice((1.0, 2.0, 3.0)))
               for _ in range(20_000)]
    edge = [(_bigs(rng), rng.uniform(0.5, 3.0), _bigs(rng), rng.uniform(0.5, 3.0)) for _ in range(30_000)]
    edge += [(_spp(rng, -1074, -1000), _pp(rng, -40, 1), _spp(rng, -1074, -1000), _pp(rng, -40, 1))
             for _ in range(10_000)]
    bad = [(INF, 0.0, 1.0, 1.0), (NAN, 1.0, 1.0, 1.0), (1.0, INF, 1.0, 1.0), (DMAX, 2.0, -DMAX, 1.5)]
    return [("binade-uniform", rand), ("typical", typical), ("products past the range", edge), ("specials", bad)]


def args_zoc(rng):
    rand = [(_spp(rng), _spp(rng), _spp(rng)) for _ in range(40_000)]
    typical = [(rng.uniform(0, 2), rng.uniform(0, 0.1), rng.uniform(0, 100)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 1)
        if kind == 0:                               # k·t in (−2·DBL_MAX, −DBL_MAX), [A]₀ negative
            t = _pp(rng, 1, 60)
            edge.append((-DMAX * rng.uniform(0.2, 1.0), -(DMAX / t) * rng.uniform(1.0, 1.99), t))
        else:                                       # k·t subnormal or below 2^-1075, [A]₀ tiny or ±0
            t = _pp(rng, -30, 20)
            edge.append((rng.choice((0.0, -0.0, _spp(rng, -1074, -1000))), _spp(rng, -1074, -1015) / t, t))
    bad = [(NAN, 1.0, 1.0), (1.0, INF, 1.0), (-0.0, 1e-200, 1e-200), (0.0, 1e-200, -1e-200), (-DMAX, -DMAX, 1.5)]
    return [("binade-uniform", rand), ("typical", typical), ("k·t past the range", edge), ("specials", bad)]


def args_shomate(rng):
    def co():
        return 0.0 if rng.random() < 0.25 else _spp(rng)
    rand = [(co(), co(), co(), co(), co(), _spp(rng)) for _ in range(40_000)]
    typical = [(rng.uniform(-50, 50), rng.uniform(-50, 50), rng.uniform(-20, 20), rng.uniform(-5, 5), rng.uniform(-2, 2),
                rng.uniform(100, 6000)) for _ in range(20_000)]
    edge = []
    for _ in range(40_000):
        kind = rng.randint(0, 3)
        phys = (rng.uniform(-50, 50), rng.uniform(-50, 50), rng.uniform(-20, 20), rng.uniform(-5, 5),
                rng.choice((0.0, rng.uniform(-2, 2), _spp(rng, -1074, -900))))
        if kind == 0:                               # T → 0: t·t underflowing
            edge.append((*phys, _spp(rng, -1074, -150)))
        elif kind == 1:                             # large T: C·t², D·t³, t⁴ overflowing
            edge.append((*phys, _spp(rng, 150, 1023)))
        elif kind == 2:                             # an exact cancellation of two overflowing terms at t = 2^k
            k = rng.randint(400, 1000)
            edge.append((rng.uniform(-50, 50), 0.0, 1.0, -(2.0 ** -k), 0.0, 125.0 * 2.0 ** (k + 3)))
        else:                                       # every binade of the coefficients at a tiny or huge T
            edge.append((co(), co(), co(), co(), co(), rng.choice((_spp(rng, -1074, -100), _spp(rng, 100, 1023)))))
    bad = [(1.0, 1.0, 1.0, 1.0, 1.0, 0.0), (20.786, 0.0, 0.0, 0.0, 0.0, -0.0), (1.0, 1.0, 1.0, 1.0, 1.0, INF),
           (NAN, 1.0, 1.0, 1.0, 1.0, 300.0), (20.786, 0.0, 0.0, 0.0, 0.0, 1e-300), (20.786, 0.0, 0.0, 0.0, 0.0, fb(1)),
           (0.0, 0.0, 1.0, -1.0, 0.0, 1e200)]
    return [("binade-uniform", rand), ("typical", typical), ("T → 0, T large, cancellation", edge),
            ("T = ±0 and specials", bad)]


# ── The overflow audit's B5 (ADR 0002 #35): the exp / exp10 arguments ────────────────────────────────────────────────
# ids 83 … 86: keq (equilibrium_constant: ΔG, T), dhl (debye_huckel_limiting: z, I, A), dhe (debye_huckel_extended:
# z, I, A, B, å), dav (davies_activity: z, I, A), the charge z an integer word: mode 0 the public function, Rust's chain
# where the function's fast test passes (R·T in [2^-1021, DBL_MAX]; A·z² finite; the extended law's numerator and
# (B·å)·√I finite) and the flag; mode 1 the unbounded-exponent argument forced (_km_keq_arg_ux, _km_dh_limiting_arg_ux,
# _km_dh_extended_arg_ux, _km_davies_arg_ux), then the correctly rounded exp or 10^x. Reference: Rust's validation in
# Rust's order, then the argument in Rust's sequence on mpfr_ref.UX with every operation at 53 bits, the last one too,
# and MPFR's exp or exp10 of that 53-bit value rounded once (mpfr_ref.ux_exp_final, ux_exp10_final): ADR 0002 #20's
# target itself, so every argument checks ADR 0002 #35's lemma as well (the port rounds the argument's last operation
# into binary64 first). An Err is the word 2^63 | its KM_ERR_* code (K and γ are never negative), in both slots and
# modes; equilibrium_constant's ComputationError where its K is not finite.
AUDIT5_NAMES = ("keq", "dhl", "dhe", "dav")


def _err_code(name):
    import re
    m = re.search(rf"^\s+{name} = ([0-9]+);", (_ROOT / "src/error.cyr").read_text(), re.M)
    return int(m.group(1))


ERR5 = {n: (1 << 63) | _err_code("KM_ERR_" + n) for n in ("INVALID_TEMPERATURE", "INVALID_INPUT", "COMPUTATION_ERROR")}
KEQ_T_2M1022 = 541_658_533_363_376      # T = k·2^-1074 with R·k = 2^52 − 0.27…: IEEE's R·T rounds up to 2^-1022
DAVIES_ZERO_I = (1.940363884733242, 1.9403638847332425)     # √I/(1 + √I) = RN(0.3·I): Davies' bracket exactly 0


def ref_audit5(name, a):
    """The reference word of ids 83 … 86 (an Err as 2^63 | its code)."""
    UX = mpfr_ref.UX
    if name == "keq":
        dg, t = a
        if not _pos(t):
            return ERR5["INVALID_TEMPERATURE"]
        if not math.isfinite(dg):
            return ERR5["INVALID_INPUT"]
        with mpfr_ref.unbounded():
            x = UX(-dg) / (UX(R_GAS) * UX(t))
        k = mpfr_ref.ux_exp_final(x)
        return bits(k) if math.isfinite(k) else ERR5["COMPUTATION_ERROR"]
    z, ii, aa = a[0], a[1], a[2]
    if not (-(1 << 31) <= z < (1 << 31) and _nonneg(ii) and math.isfinite(aa)):
        return ERR5["INVALID_INPUT"]
    if name == "dhe" and not (_nonneg(a[3]) and _nonneg(a[4])):
        return ERR5["INVALID_INPUT"]
    with mpfr_ref.unbounded():
        zf = UX(float(z))
        p = UX(-aa) * zf * zf
        sq = UX(ii).sqrt()
        if name == "dhl":
            x = p * sq
        elif name == "dhe":
            x = (p * sq) / (UX(1.0) + UX(a[3]) * UX(a[4]) * sq)
        else:
            x = p * (sq / (UX(1.0) + sq) - UX(0.3) * UX(ii))
    return bits(mpfr_ref.ux_exp10_final(x))


def _z5(rng, big=False):
    """A charge of a random magnitude (|z| ≥ 2 when big, so that A·z² can pass DBL_MAX), a random sign."""
    k = rng.randint(0, 3)
    m = (rng.randint(2, 4), rng.randint(5, 1000), rng.randint(1001, 100_000), rng.randint(100_001, (1 << 31) - 1))[k]
    m = m if big else m - 1
    return m * rng.choice((1, -1))


def _past(rng, z):
    """A with |A·z²| past DBL_MAX (|z| ≥ 2), A itself finite."""
    lo = DMAX / float(z) / float(z)
    return min(lo * rng.uniform(1.0001, 4.0), DMAX) * rng.choice((1.0, -1.0))


def args_keq(rng):
    rand = [(_spp(rng), _pp(rng)) for _ in range(40_000)]
    typical = [(rng.uniform(-1e6, 1e6), rng.uniform(1.0, 5000.0)) for _ in range(20_000)]
    over = []
    for _ in range(30_000):                         # R·T past DBL_MAX: |x| < 1/R, Rust's K = 1
        t = DMAX * rng.uniform(0.121, 1.0)
        over.append((rng.choice((_spp(rng, 900, 1023), _spp(rng))), t))
    sub = []
    for _ in range(30_000):                         # R·T subnormal: ΔG = −x·R·T for x across exp's range
        t = _pp(rng, -1074, -1025)
        sub.append((-rng.uniform(-760.0, 720.0) * R_GAS * t, t))
    ends = []
    for _ in range(20_000):                         # a short subnormal T and ΔG between the two ends' thresholds
        k = rng.randint(1, 1 << 24)
        rk53, rki = R_GAS * k, float(bits(R_GAS * fb(k)))
        mid = 0.5 * (rk53 + rki)
        xe, sg = rng.choice(((709.782712893384, -1.0), (745.1332191019411, 1.0)))
        n = int(xe * mid) + rng.randint(-2, 2)
        ends.append((sg * fb(max(n, 0)), fb(k)))
    bound = []
    t_ov, t_un = DMAX / R_GAS, 2.0 ** -1021 / R_GAS
    for d in range(-20, 21):                        # the boundaries T = DBL_MAX/R and 2^-1021/R, R·T at 2^-1022
        for t in (fb(bits(t_ov) + d), fb(bits(t_un) + d), fb(KEQ_T_2M1022 + d)):
            for _ in range(40):
                x = rng.uniform(-760.0, 720.0)
                bound.append((-x * 2.0 ** -1022 if t < 1.0 else -DMAX * rng.uniform(0.0, 1.0), t))
    bad = [(1.0, 0.0), (1.0, -1.0), (NAN, 1.0), (1.0, NAN), (INF, 1.0), (1.0, INF), (-2e6, 298.0), (2e6, 298.0),
           (-0.0, fb(1)), (0.0, DMAX), (-DMAX, DMAX), (-1e-320, 1e-322), (-DMAX, fb(1)), (DMAX, fb(1))]
    return [("binade-uniform", rand), ("typical", typical), ("R·T past DBL_MAX", over), ("R·T subnormal", sub),
            ("the ends at a short subnormal T", ends), ("T at the boundaries", bound), ("Err and specials", bad)]


def args_dh(rng, name):
    rand = []
    for _ in range(40_000):
        z = _z5(rng)
        rand.append((z, _pp(rng), _spp(rng)) if name != "dhe" else (z, _pp(rng), _spp(rng), _pp(rng), _pp(rng)))
    typical = []
    for _ in range(20_000):
        z, ii, aa = rng.randint(-4, 4), rng.uniform(0.0, 1.0), rng.uniform(0.4, 0.6)
        typical.append((z, ii, aa) if name != "dhe" else (z, ii, aa, rng.uniform(0.3, 0.35), rng.uniform(2.0, 10.0)))
    past = []
    for _ in range(30_000):                         # A·z² past DBL_MAX at I = ±0, Davies' zero bracket, or brought back
        z = _z5(rng, True)
        ii = rng.choice((0.0, -0.0, rng.choice(DAVIES_ZERO_I), _pp(rng), rng.uniform(1.9, 2.0)))
        aa = _past(rng, z)
        if name != "dhe":
            past.append((z, ii, aa))
        else:
            sq = math.sqrt(ii)
            lv = 10.0 ** rng.uniform(-3.0, 2.5)
            r = math.sqrt(abs(aa)) * abs(z) / math.sqrt(lv * sq) if sq > 0 else _pp(rng)
            u = 10.0 ** rng.uniform(-3.0, 3.0)
            b_, s_ = r * u, r / u
            past.append((z, ii, aa, b_ if math.isfinite(b_) else DMAX, s_ if math.isfinite(s_) else 1.0))
    small = []
    for _ in range(10_000):                         # A·z² subnormal (exact in both)
        z = _z5(rng)
        ii, aa = _pp(rng, -40, 40), _spp(rng, -1074, -1023)
        small.append((z, ii, aa) if name != "dhe" else (z, ii, aa, _pp(rng, -40, 40), _pp(rng, -40, 40)))
    extra = []
    if name == "dhe":
        for _ in range(30_000):
            kind = rng.randint(0, 4)
            z = _z5(rng)
            if kind == 0:                           # B·å past DBL_MAX, the numerator finite: log₁₀ γ ≈ −A·z²/(B·å)
                b_, s_ = _pp(rng, 512, 1023), _pp(rng, 512, 1023)
                aa = 10.0 ** rng.uniform(-6.0, 2.5) * rng.choice((1.0, -1.0)) * (b_ / z) * (s_ / z)
                extra.append((z, _pp(rng, -300, 300), aa if math.isfinite(aa) else 1.0, b_, s_))
            elif kind == 1:                         # the numerator past DBL_MAX with A·z² finite
                aa = (DMAX / z / z) * rng.uniform(0.01, 1.0) * rng.choice((1.0, -1.0))
                extra.append((z, _pp(rng, 2, 1000), aa, _pp(rng, -10, 1023), _pp(rng, -10, 1023)))
            elif kind == 2:                         # I = ±0 with B·å or A·z² past DBL_MAX
                extra.append((_z5(rng, True), rng.choice((0.0, -0.0)), _spp(rng, 800, 1023), _pp(rng, 512, 1023),
                              _pp(rng, 512, 1023)))
            elif kind == 3:                         # the numerator subnormal, or IEEE's rounded up to 2^-1022 from below
                ii = _pp(rng, 0, 40)
                k = int(2.0 ** 52 / math.sqrt(ii)) + rng.randint(-1, 1)
                aa = -fb(max(1, min(k, (1 << 52) - 1))) * rng.choice((1.0, -1.0))
                extra.append((1, ii, aa, rng.choice((0.0, 0.328)), rng.choice((0.0, 3.0))))
            else:                                   # B·å subnormal (the denominator 1 in both), B or å ±0
                extra.append((z, _pp(rng, -100, 1023), _spp(rng, -10, 10), rng.choice((_pp(rng, -1074, -1023), -0.0)),
                              rng.choice((_pp(rng, -40, 1), 0.0, -0.0))))
    elif name == "dav":
        for _ in range(30_000):                     # the bracket near 0 (I ≈ 1.9404), negative, or I subnormal
            ii = rng.choice((rng.uniform(1.9403, 1.9404), rng.uniform(1.9, 2.0), _pp(rng, 1, 1023), _pp(rng, -1074, -1019)))
            extra.append((_z5(rng), ii, _spp(rng, -10, 1023)))
    bad = [(1, -0.1, 0.509), (1, NAN, 0.509), (1, INF, 0.509), (1, 0.01, NAN), (1, 0.01, INF), (1 << 31, 0.01, 0.509),
           (-(1 << 31) - 1, 0.01, 0.509), (100_000, 0.0, 1e300), (100_000, -0.0, -1e300), (0, 0.0, DMAX),
           ((1 << 31) - 1, 0.0, 1e290), (-(1 << 31), 1.0, DMAX)]
    if name == "dhe":
        bad = [(z, ii, aa, 0.328, 3.0) for z, ii, aa in bad] + [
            (1, 0.01, 0.509, -0.328, 3.0), (1, 0.01, 0.509, 0.328, -3.0), (1, 0.01, 0.509, NAN, 3.0),
            (1, 0.01, 0.509, 0.328, INF), (1, 1.0, 1e308, 1e155, 1e155), (100_000, 1.0, 1e300, 1e154, 1e154),
            (1, 0.0, 1.0, DMAX, 2.0), (2, 1.0, DMAX, DMAX, DMAX), (1, 0.01, 0.509, -0.0, -0.0)]
    elif name == "dav":
        bad += [(2, DAVIES_ZERO_I[0], DMAX), (-2, DAVIES_ZERO_I[1], -DMAX), (1, DAVIES_ZERO_I[0], 0.509)]
    return [("binade-uniform", rand), ("typical", typical), ("A·z² past DBL_MAX", past), ("A·z² subnormal", small),
            ("the function's own edges", extra), ("Err and specials", bad)]


# ── The overflow audit's B6 (ADR 0002 #36, #37, #38): the FFT, find_peaks' threshold, the flame's bracket ────────────
# id 87, flame (adiabatic_flame_temperature: ΔH, T₀, then three (species, moles) pairs, a species an index into
# FLAME6_SPECIES, −1 unused): the public result word (an Err as 2^63 | its code), and where Rust's validation passes (Ok
# or the ComputationError) the transcribed bisection over [T₀, T_max] (its exact same-sign test) in the second slot, else
# the public word; both modes. Reference: Rust's solver from src/thermochem_data.cyr's SHOMATE table (generated), the
# Shomate polynomials in Rust's order (every evaluation checked inside ADR 0002 #34's box, where they are the port's),
# hisab 1.4.0's newton_raphson and bisection, with ADR 0002 #38's exact sign tests. Ids 88 … 91, pk3 … pk16 (find_peaks
# of 3, 5, 8, 16 values, then t): the peaks as a bit mask; mode 0 the public function, the scan with the fast test and
# the flag 1 when no candidate needed the exact predicate, mode 1 every candidate through the exact predicate. Reference:
# Rust's checks, fold and neighbour tests with ADR 0002 #37's threshold test in Fraction arithmetic.
AUDIT6_NAMES = ("flame", "pk3", "pk5", "pk8", "pk16", "signsc")
FLAME6_SPECIES = ("CO2(g)", "H2O(g)", "N2(g)", "CO(g)", "O2(g)", "H2(g)", "Ar(g)")


def _shomate_table():
    """SHOMATE_DATA from src/thermochem_data.cyr (generated): (formula, A, B, C, D, E, t_min, t_max) per row."""
    import re
    text = (_ROOT / "src/thermochem_data.cyr").read_text()
    off = text[text.index("var _KM_SHOMATE_FORMULA_OFF"):]
    off = off[:off.index("};")]
    names = re.findall(r"^\s+[0-9]+,\s+#\s+(\S+)$", off, re.M)

    def col(name):
        body = text[text.index(f"var {name}: f64[34] = {{"):]
        body = body[:body.index("};")]
        return [float(v) for v in re.findall(r"^\s+(-?[0-9.]+),", body, re.M)]
    cols = [col("KM_SHOMATE_" + c) for c in ("A", "B", "C", "D", "E", "T_MIN", "T_MAX")]
    assert len(names) == 34 and all(len(c) == 34 for c in cols)
    return [(names[i], *(c[i] for c in cols)) for i in range(34)]


SHOMATE6 = _shomate_table()


def _segs6(formula):
    rows = [r for r in SHOMATE6 if r[0] == formula]
    return rows if rows else None


def _sh6(is_cp, r, temp):
    assert 2.0 ** -30 <= abs(temp) <= 2.0 ** 30, f"a flame evaluation at {temp!r}, outside ADR 0002 #34's box"
    return _shomate_rust(is_cp, *r[1:6], temp)


def _seg_dh6(run, t1, t2):
    lo, hi, sign = (t1, t2, 1.0) if t1 <= t2 else (t2, t1, -1.0)
    total, lower = 0.0, -INF
    for k, r in enumerate(run):
        upper = INF if k + 1 == len(run) else r[7]
        a_, b_ = max(lo, lower), min(hi, upper)
        if b_ > a_:
            total += _sh6(False, r, b_) - _sh6(False, r, a_)
        if upper >= hi:
            break
        lower = upper
    return sign * total


def _seg_cp6(run, t):
    for r in run:
        if t <= r[7]:
            return _sh6(True, r, t)
    return _sh6(True, run[-1], t)


def ref_flame6(a):
    """(the public word, the bisection's word) of id 87."""
    dh, t0 = a[0], a[1]
    prods = [(FLAME6_SPECIES[a[2 + 2 * j]], a[3 + 2 * j]) for j in range(3) if a[2 + 2 * j] >= 0]
    if not _pos(t0):
        w = ERR5["INVALID_TEMPERATURE"]
        return w, w
    if not math.isfinite(dh) or not prods:
        return ERR5["INVALID_INPUT"], ERR5["INVALID_INPUT"]
    res, tmax = [], INF
    for f_, n_ in prods:
        run = _segs6(f_)
        if run is None:
            w = (1 << 63) | _err_code("KM_ERR_INVALID_REACTION")
            return w, w
        if not _pos(n_):
            return ERR5["INVALID_INPUT"], ERR5["INVALID_INPUT"]
        tmax = min(tmax, run[-1][7])
        res.append((run, n_, run[0][6]))
    tmin = max(r_[2] for r_ in res)
    if t0 < tmin or t0 >= tmax:
        w = ERR5["INVALID_TEMPERATURE"]
        return w, w

    def fun(t):
        heat = -0.0
        for run, n_, _ in res:
            heat += n_ * _seg_dh6(run, t0, t)
        return -dh - heat

    def dfun(t):
        cp = -0.0
        for run, n_, _ in res:
            cp += n_ * _seg_cp6(run, t)
        return -cp

    def bisect(lo, hi):
        f_lo, f_hi = fun(lo), fun(hi)
        if (f_lo > 0.0 and f_hi > 0.0) or (f_lo < 0.0 and f_hi < 0.0):
            return (1 << 63) | _err_code("KM_ERR_INVALID_INPUT")
        if f_lo > 0.0:
            lo, hi = hi, lo
        for _ in range(200):
            mid = (lo + hi) * 0.5
            fm = fun(mid)
            if abs(fm) < 1e-6 or abs(hi - lo) < 1e-6:
                return bits(mid)
            if fm < 0.0:
                lo = mid
            else:
                hi = mid
        return bits((lo + hi) * 0.5)
    wb = bisect(t0, tmax)
    fa, fb_ = fun(t0), fun(tmax)
    if not ((fa < 0.0 < fb_) or (fb_ < 0.0 < fa)):
        return ERR5["COMPUTATION_ERROR"], wb
    total = -0.0
    for _, n_, _ in res:
        total += n_
    x = min(max(t0 - dh / (total * 40.0), t0), tmax)
    for _ in range(50):
        fx = fun(x)
        if abs(fx) < 1e-3:
            if t0 <= x <= tmax:
                return bits(x), wb
            break
        dfx = dfun(x)
        if abs(dfx) < 1e-12:
            break
        x -= fx / dfx
    return wb, wb


def ref_peaks6(a):
    """The bit mask of id 88 … 91's peaks (ADR 0002 #37's exact threshold test)."""
    from fractions import Fraction
    data, th = list(a[:-1]), a[-1]
    if len(data) < 3 or not (0.0 <= th <= 1.0):
        return 0
    fin = [v for v in data if math.isfinite(v)]
    if not fin:
        return 0
    lo, hi = INF, -INF
    for v in fin:                                   # Rust's fold: f64::min / max (equal values: either, the same)
        lo, hi = min(lo, v), max(hi, v)
    q = Fraction(lo) * (1 - Fraction(th)) + Fraction(hi) * Fraction(th)
    mask = 0
    for k in range(1, len(data) - 1):
        w = data[k]
        if w > data[k - 1] and w > data[k + 1] and (Fraction(w) > q if math.isfinite(w) else w > 0.0):
            mask |= 1 << k
    return mask


def ref_audit6(name, a):
    if name == "flame":
        return ref_flame6(a)[0]
    if name == "signsc":
        from fractions import Fraction
        n = a[0]
        tot = sum((Fraction(a[1 + i]) * Fraction(2) ** a[9 + i] for i in range(n)), Fraction(0))
        return (tot > 0) - (tot < 0) + 1
    return ref_peaks6(a)


def args_signsc(rng):
    """(n, eight significands, eight exponents) for _km_sign_scaled: random terms over its whole exponent range; exact
    cancellations at the top (c·2^k against −c·2^k or −2c·2^(k − 1)) with small remainders 850 … 2,210 binades below and
    deferred terms around them, so that the sign is decided one, two or three rounds down; a round's remainder of
    2^-1010 … 2^-901 of its scale against a larger deferred term (the bound B6's review corrected); and exact zeros."""
    def sig():
        e = rng.randint(-110, 1)
        return rng.choice((1.0, -1.0)) * rng.uniform(1.0, 2.0) * 2.0 ** e

    def pack(cs, ks):
        cs, ks = list(cs)[:8], list(ks)[:8]
        n = len(cs)
        return (n, *(cs + [0.0] * (8 - n)), *(ks + [0] * (8 - n)))
    rand = [pack([sig() for _ in range(m)], [rng.randint(-2200, 1030) for _ in range(m)])
            for m in (rng.randint(1, 8) for _ in range(40_000))]
    canc = []
    for _ in range(60_000):
        k = rng.randint(-1000, 1000)
        c = sig()
        if rng.random() < 0.5:
            cs, ks = [c, -c], [k, k]
        else:
            cs, ks = [c, -2.0 * c], [k, k - 1] if abs(2.0 * c) < 4.0 else [k, k]
            if ks[1] == k:
                cs[1] = -c
        for _ in range(rng.randint(1, 6)):
            cs.append(sig())
            ks.append(k - rng.randint(850, 1110) - rng.choice((0, 0, rng.randint(0, 1100))))
        canc.append(pack(cs, ks))
    edge = []
    for _ in range(40_000):
        k = rng.randint(-900, 1000)
        c = sig()
        small = rng.choice((1.0, -1.0)) * 2.0 ** rng.randint(-110, -1)
        cs = [c, -c, small, rng.choice((1.0, -1.0)) * rng.uniform(1.0, 3.99)]
        ks = [k, k, k - 900, k - rng.randint(901, 1000)]
        edge.append(pack(cs, ks))
    zeros = [pack([1.0, -1.0], [k, k]) for k in range(-50, 50)] + [pack([0.0, 2.0, -1.0], [5, 3, 4]) for _ in range(10)]
    return [("binade-uniform", rand), ("cancellations", canc), ("the round's bound", edge), ("zeros", zeros)]


def args_flame6(rng):
    """(ΔH, T₀, s₁, n₁, s₂, n₂, s₃, n₃) with real systems at tiny scales (the bracket product underflowing: Rust's Err),
    across the scale where it stops underflowing, endothermic ones (both ends negative: Err), subnormal ΔH against
    moles over many decades, ordinary ones, and every Err path."""
    systems = [(-802_300.0, ((0, 1.0), (1, 2.0), (2, 7.52))), (-241_820.0, ((1, 1.0), (2, 1.88))),
               (-519_320.0, ((3, 1.0), (1, 2.0), (2, 5.64))), (-282_980.0, ((0, 1.0), (2, 1.88))),
               (-393_510.0, ((0, 1.0),)), (-241_820.0, ((1, 1.0), (5, 0.2), (6, 0.5)))]

    def pack(dh, t0, prods, k):
        pr = list(prods) + [(-1, 0.0)] * (3 - len(prods))
        return (dh * k, t0, *(v for s_, n_ in pr for v in (s_, n_ * k)))
    tiny = []
    for _ in range(6000):
        dh, prods = rng.choice(systems)
        tiny.append(pack(dh, rng.choice((298.15, 300.0, 500.0, 1000.0)), prods, 10.0 ** rng.uniform(-320.0, -150.0)))
    edge = []
    for _ in range(3000):                           # around the scale 10^-169 where the product stops underflowing
        dh, prods = rng.choice(systems)
        edge.append(pack(dh, 298.15, prods, 10.0 ** rng.uniform(-172.0, -166.0)))
    endo = []
    for _ in range(2000):
        dh, prods = rng.choice(systems)
        endo.append(pack(-dh, rng.choice((298.15, 500.0)), prods, 10.0 ** rng.uniform(-320.0, 0.0)))
    sub = []
    for _ in range(2000):                           # a subnormal ΔH against moles from 10^-9 to 10^2
        k = rng.randint(1, (1 << 52) - 1)
        sub.append((-fb(k) * rng.choice((1.0, -1.0)), 298.15, 0, 10.0 ** rng.uniform(-9.0, 2.0), -1, 0.0, -1, 0.0))
    typical = []
    for _ in range(1000):
        dh, prods = rng.choice(systems)
        typical.append(pack(dh * rng.uniform(0.2, 2.0), rng.choice((298.15, 400.0, 700.0)), prods, 1.0))
    bad = [(-1e5, 0.0, 0, 1.0, -1, 0.0, -1, 0.0), (NAN, 298.15, 0, 1.0, -1, 0.0, -1, 0.0),
           (-1e5, 298.15, 0, -1.0, -1, 0.0, -1, 0.0), (-1e5, 298.15, -1, 0.0, -1, 0.0, -1, 0.0),
           (-1e5, 7000.0, 0, 1.0, -1, 0.0, -1, 0.0), (0.0, 298.15, 0, 1.0, -1, 0.0, -1, 0.0),
           (-0.0, 298.15, 0, 1.0, -1, 0.0, -1, 0.0), (-1e-320, 298.15, 0, 1e-300, -1, 0.0, -1, 0.0),
           (1e-320, 298.15, 0, 1e-300, -1, 0.0, -1, 0.0), (-1e5, 298.15, 0, 5e-324, -1, 0.0, -1, 0.0)]
    return [("tiny scales", tiny), ("around the underflow's scale", edge), ("endothermic", endo),
            ("a subnormal ΔH", sub), ("typical", typical), ("Err and specials", bad)]


def args_peaks6(rng, n):
    """n values and a threshold: binade-uniform, typical, near-ties (a candidate at the exact threshold's neighbours),
    subnormal, ±DBL_MAX ends, ±0 and non-finite values."""
    from fractions import Fraction
    ts = (0.0, 1.0, nextdown_(1.0), fb(1), 0.5, nextdown_(0.5), nextup_(0.5), 0.1, 0.75)

    def t_of():
        return rng.choice(ts) if rng.random() < 0.3 else rng.random()
    rand = [tuple(_spp(rng) for _ in range(n)) + (t_of(),) for _ in range(20_000)]
    typical = [tuple(rng.uniform(-10.0, 10.0) for _ in range(n)) + (t_of(),) for _ in range(20_000)]
    ties = []
    for _ in range(30_000):
        e0 = rng.randint(-1074, 1000)
        xs = [_spp(rng, e0, min(e0 + rng.randint(0, 60), 1023)) for _ in range(n)]
        t = t_of()
        lo, hi = min(xs), max(xs)
        q = Fraction(lo) * (1 - Fraction(t)) + Fraction(hi) * Fraction(t)
        w = float(q)                                # the double nearest θ, then a neighbour
        w = ulps(w, rng.randint(-2, 2))
        if lo <= w <= hi:
            xs[rng.randint(1, n - 2)] = w
        ties.append(tuple(xs) + (t,))
    sub = [tuple(fb(rng.randint(0, 8)) * rng.choice((1.0, -1.0)) for _ in range(n)) + (t_of(),) for _ in range(10_000)]
    wide = []
    for _ in range(10_000):
        xs = [rng.choice((DMAX, -DMAX, _spp(rng, 1000, 1023), _spp(rng))) for _ in range(n)]
        wide.append(tuple(xs) + (rng.choice((nextdown_(0.5), nextup_(0.5), 0.5, t_of())),))
    special = []
    for _ in range(5_000):
        xs = [rng.choice((0.0, -0.0, INF, -INF, NAN, fb(1), 1.0, DMAX)) for _ in range(n)]
        special.append(tuple(xs) + (rng.choice((t_of(), -0.1, 1.5, NAN, INF)),))
    return [("binade-uniform", rand), ("typical", typical), ("near-ties", ties), ("subnormal", sub),
            ("±DBL_MAX ends", wide), ("±0 and non-finite", special)]


def nextdown_(x):
    return math.nextafter(x, -INF)


def nextup_(x):
    return math.nextafter(x, INF)


# ids 92 … 99, fft1 … fft256 (fft_spectrum of n = 1, 2, 4, …, 64 and 256 values): n/2 + 1 words per slot; mode 0 the
# public spectrum, the fast path's (_km_fft_spectrum_fast: Rust's IEEE transform, each bin's magnitude Rust's or the
# unbounded one) and the flag 1 when its transform test passed, mode 1 the
# unbounded-exponent spectrum forced (_km_fft_spectrum_ux; a non-finite signal's public one). Reference: hisab 1.4.0's
# transform on mpfr_ref.UX (every operation at 53 bits, MPFR's correctly rounded cos and sin for the twiddles) and each
# bin's √ rounded once (mpfr_ref.ux_sqrt_final); for a signal with a NaN or ±∞, Rust's IEEE transform and magnitudes
# in binary64 (B1's rule). Sizes and counts (the reference costs about 16 UX operations a butterfly, so the counts fall
# with n): per region 3,000 signals for n = 1 … 8, 1,500 for 16, 800 for 32, 400 for 64 and 60 for 256 (seven regions
# and half one of non-finite samples: 22,500 … 450 signals a size, about 110,718 in all).
FFT6 = [(92, "fft1", 1, 3000), (93, "fft2", 2, 3000), (94, "fft4", 4, 3000), (95, "fft8", 8, 3000),
        (96, "fft16", 16, 1500), (97, "fft32", 32, 800), (98, "fft64", 64, 400), (99, "fft256", 256, 60)]
for _fid, _nm, _n6, _c6 in FFT6:
    ARITY[_fid] = _n6


def _bitrev6(n):
    perm, j = list(range(n)), 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            perm[i], perm[j] = perm[j], perm[i]
    return perm


def _fft6(sig, U):
    """hisab 1.4.0's transform with the value type U (mpfr_ref.UX, or float for Rust's IEEE one)."""
    n = len(sig)
    perm = _bitrev6(n)
    re_, im_ = [U(sig[k]) for k in perm], [U(0.0) for _ in sig]
    length = 2
    while length <= n:
        half = length >> 1
        ang = -2.0 * math.pi / float(length)
        wr0, wi0 = U(mpfr_ref.cos(ang)), U(mpfr_ref.sin(ang))
        for start in range(0, n, length):
            wr, wi = U(1.0), U(0.0)
            for k in range(half):
                a_, b_ = start + k, start + k + half
                tr = wr * re_[b_] - wi * im_[b_]
                ti = wr * im_[b_] + wi * re_[b_]
                ur, ui = re_[a_], im_[a_]
                re_[a_], im_[a_] = ur + tr, ui + ti
                re_[b_], im_[b_] = ur - tr, ui - ti
                wr, wi = wr * wr0 - wi * wi0, wr * wi0 + wi * wr0
        length <<= 1
    return re_, im_


def ref_fft6(sig):
    """The reference words of an FFT family argument (n/2 + 1 of them)."""
    m = len(sig) // 2 + 1
    if not all(math.isfinite(x) for x in sig):
        re_, im_ = _fft6(sig, float)
        out = []
        for k in range(m):
            sq_ = re_[k] * re_[k] + im_[k] * im_[k]
            out.append(bits(math.sqrt(sq_) if sq_ == sq_ else NAN))
        return out
    with mpfr_ref.unbounded():
        re_, im_ = _fft6(sig, mpfr_ref.UX)
        return [bits(mpfr_ref.ux_sqrt_final(re_[k] * re_[k] + im_[k] * im_[k])) for k in range(m)]


def args_fft6(rng, n, count):
    c = count
    rand = [tuple(_spp(rng) for _ in range(n)) for _ in range(c)]
    typical = [tuple(rng.uniform(-1.0, 1.0) for _ in range(n)) for _ in range(c)]
    huge = [tuple(_spp(rng, 990, 1023) for _ in range(n)) for _ in range(c)]
    subn = [tuple(_spp(rng, -1074, -1000) for _ in range(n)) for _ in range(c)]
    mixed = [tuple(rng.choice((_spp(rng, 950, 1023), _spp(rng, -1074, -950), 0.0, -0.0, _spp(rng)))
                   for _ in range(n)) for _ in range(c)]
    near = []                                       # near-cancelling samples: tiny butterfly differences
    for _ in range(c):
        x = _spp(rng, -1000, 1000)
        near.append(tuple(ulps(x, rng.randint(-3, 3)) for _ in range(n)))
    bound = []                                      # around the data bound 2^-921 and the squares' 2^-511, 2^512
    for _ in range(c):
        b0 = rng.choice((2.0 ** -921, 2.0 ** -511, 2.0 ** 511, math.sqrt(DMAX)))
        bound.append(tuple(rng.choice((0.0, ulps(b0, rng.randint(-2, 2)) * rng.choice((1.0, -1.0)), rng.uniform(-1, 1)))
                           for _ in range(n)))
    special = []
    for _ in range(max(c // 2, 1)):
        xs = [rng.choice((_spp(rng, 990, 1023), rng.uniform(-1.0, 1.0))) for _ in range(n)]
        xs[rng.randint(0, n - 1)] = rng.choice((INF, -INF, NAN))
        special.append(tuple(xs))
    if n == 4:                                      # a product w·d just below 2^-1022 (w = (cos(−π/2), −1))
        cc = 2.0 ** -1022 / math.cos(-math.pi / 2.0)
        special += [(0.0, ulps(cc, d), 0.0, 0.0) for d in range(-4, 5)]
        special += [(1.0, ulps(cc, d), -1.0, 0.0) for d in range(-4, 5)]
    return [("binade-uniform", rand), ("typical", typical), ("huge", huge), ("subnormal", subn), ("mixed", mixed),
            ("near-cancelling", near), ("the bounds", bound), ("non-finite and 2^-1022 products", special)]


def run_probe_vec(work, target, fid, mode, args, vw):
    """(slot 0 word lists, slot 1 word lists, flags) for the FFT families' arguments (floats) and the balancing
    families' (integer words)."""
    _, out, runner = TARGETS[target]
    n = len(args)
    flat = [v & 0xFFFFFFFFFFFFFFFF if isinstance(v, int) else bits(v) for a in args for v in a]
    (work / "job.bin").write_bytes(struct.pack(f"<3Q{len(flat)}Q", fid, mode, n, *flat))
    res = work / "res.bin"
    for attempt in range(2):
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
    words = struct.unpack(f"<{2 * vw * n}Q", data[:16 * vw * n])
    s0 = [list(words[vw * i:vw * i + vw]) for i in range(n)]
    s1 = [list(words[vw * n + vw * i:vw * n + vw * i + vw]) for i in range(n)]
    return s0, s1, data[16 * vw * n:16 * vw * n + n]


def fft6_checks(work, targets, which):
    """Ids 92 … 99 against the unbounded-exponent spectrum on every target, the slow path forced; True if anything
    failed."""
    failed = False

    def eq(a, b):
        return all(x == y or (math.isnan(fb(x)) and math.isnan(fb(y))) for x, y in zip(a, b))
    for fid, name, n, count in FFT6:
        if name not in which:
            continue
        t0 = time.time()
        rng = random.Random(1033 + fid)
        sections = args_fft6(rng, n, count)
        xs = [x for _, sec in sections for x in sec]
        want = [ref_fft6(a) for a in xs]
        vw = n // 2 + 1
        first = None
        for t in targets:
            g0, g1, fl = run_probe_vec(work, t, fid, 0, xs, vw)
            a0, _, _ = run_probe_vec(work, t, fid, 1, xs, vw)
            bad = [i for i in range(len(xs)) if not eq(g0[i], want[i])]
            bad_acc = [i for i in range(len(xs)) if not eq(a0[i], want[i])]
            bad_fast = [i for i in range(len(xs)) if fl[i] == 1 and not eq(g1[i], g0[i])]
            for i in bad[:3]:
                print(f"  MISMATCH {name}{xs[i]!r} on {t}: got {[hex(w) for w in g0[i]]}, "
                      f"want {[hex(w) for w in want[i]]}")
            for i in bad_acc[:3]:
                print(f"  MISMATCH {name} (slow path){xs[i]!r} on {t}: got {[hex(w) for w in a0[i]]}")
            for i in bad_fast[:3]:
                print(f"  FAST/FULL DIFFER {name}{xs[i]!r} on {t}")

            def canon(ws):
                return [[NANW if math.isnan(fb(w)) else w for w in v] for v in ws]
            cur = (canon(g0), canon(g1), list(fl), canon(a0))
            xdiff = 0 if first is None else sum(1 for u, v in zip(cur, first) for p, q in zip(u, v) if p != q)
            if first is None:
                first = cur
            slow = len(xs) - sum(1 for f in fl if f == 1)
            print(f"{name:7} {t:8} {len(xs):9,} args: {len(bad)} mismatches; slow path forced: {len(bad_acc)} "
                  f"mismatches; fast path: {len(bad_fast)} differ, fell back {slow:,}; bits or flags differing from "
                  f"{targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0]:
                rates, i = [], 0
                for tag, sec in sections:
                    rates.append(f"{tag} {sum(1 for f in fl[i:i + len(sec)] if f != 1)}/{len(sec)}")
                    i += len(sec)
                print("         fallbacks by region: " + "; ".join(rates))
        print(f"       [{time.time() - t0:.0f} s]")
    return failed


# ── The overflow audit's B7 (ADR 0002 #39): balance_equation on unbounded integers ─────────────────────────────────
# Ids 101 … 106: balance_equation of m × n integer matrices (the m·n entries as words), n + 1 words a slot: the status
# (n for Ok, 2⁶³ | the Err code) and the n coefficients (0 where Err). Mode 0: the public result (slot 0), the i64 fast
# path's own outcome (_km_balance_i64: slot 1 its code 0, 1 or 2, then its coefficients where 0) and the flag 1 when it
# did not overflow (code 0 or 1); mode 1 the multi-word path forced (_km_balance_equation with force) in both slots.
# Reference: Rust's algorithm on Python ints (ref_bal7, exact) with a tracer of Rust's checked i64 operations, so the
# fast path must report an overflow exactly where Rust's arithmetic overflows, and every result must be the exact one.
BAL7 = [(101, "bal23", 2, 3, 3000), (102, "bal34", 3, 4, 2500), (103, "bal56", 5, 6, 1500), (104, "bal89", 8, 9, 600),
        (105, "bal1213", 12, 13, 150), (106, "bal2425", 24, 25, 8)]
for _fid, _nm, _m7, _n7, _c7 in BAL7:
    ARITY[_fid] = _m7 * _n7
BAL7_ERR = {"InvalidInput": 5, "ComputationError": 6}


def ref_bal7(a, m, n):
    """(the expected n + 1 words of balance_equation, whether Rust's checked i64 arithmetic overflows): Rust 1.2.0's
    algorithm on Python ints, every value Rust forms with a checked operation tested against i64 until it returns."""
    lo, hi = -2 ** 63, 2 ** 63 - 1
    ovf = [False]

    def c(v):
        if not lo <= v <= hi:
            ovf[0] = True
        return v

    def tdiv(p, q):
        r = abs(p) // abs(q)
        return r if (p >= 0) == (q > 0) else -r

    def red(row):
        g = 0
        for v in row:
            g = math.gcd(g, abs(v))
        if g > 1:
            c(g)
            return [tdiv(v, g) for v in row]
        return row

    def err(name):
        return [(1 << 63) | BAL7_ERR[name]] + [0] * n, ovf[0]
    if m < 1 or n < 2 or any(abs(v) > 2 ** 53 for v in a):
        return err("InvalidInput")
    A = [list(a[i * n:(i + 1) * n]) for i in range(m)]
    rank, free = 0, None
    for col in range(n):
        p = next((r for r in range(rank, m) if A[r][col] != 0), None)
        if p is None:
            if free is not None:
                return err("ComputationError")
            free = col
            continue
        A[p], A[rank] = A[rank], A[p]
        A[rank] = red(A[rank])
        if A[rank][col] < 0:
            A[rank] = [c(-v) for v in A[rank]]
        piv = A[rank][col]
        for r in range(m):
            fac = A[r][col]
            if r == rank or fac == 0:
                continue
            A[r] = red([c(c(piv * A[r][j]) - c(fac * A[rank][j])) for j in range(n)])
        rank += 1
    if free is None:
        return err("ComputationError")
    lcm = 1
    for k in range(rank):
        pv = next(v for v in A[k] if v != 0)
        g = c(math.gcd(abs(lcm), abs(pv)))
        lcm = c(tdiv(lcm, g) * pv)
    x = [0] * n
    x[free] = lcm
    for k in range(rank):
        cc = next(j for j, v in enumerate(A[k]) if v != 0)
        x[cc] = c(-c(A[k][free] * tdiv(lcm, A[k][cc])))
    x = red(x)
    if 0 in x:
        return err("ComputationError")
    pos = sum(1 for v in x if v > 0)
    if pos == 0:
        x = [c(-v) for v in x]
    elif pos != n:
        return err("ComputationError")
    if any(v < 0 or v > 2 ** 32 - 1 for v in x):
        return err("ComputationError")
    for i in range(m):
        if sum(a[i * n + j] * x[j] for j in range(n)) != 0:
            return err("ComputationError")
    return [n] + x, ovf[0]


def _bal7_solved(rng, x, rows, shrink=0, zp=0):
    """rows rows orthogonal to x (x's last entry 1): each entry but the last uniform in ±⌊2⁵³/((n − 1)·|xⱼ|)⌋ >> shrink,
    zero with probability zp, the last solved (|last| ≤ 2⁵³)."""
    n = len(x)
    out = []
    for _ in range(rows):
        row, s_ = [], 0
        for xj in x[:-1]:
            mag = max(((2 ** 53) // ((n - 1) * abs(xj)) if xj else (2 ** 53) // n) >> shrink, 1)
            e = 0 if rng.random() < zp else rng.randint(-mag, mag)
            row.append(e)
            s_ += e * xj
        row.append(-s_)
        out += row
    return out


def args_bal7(rng, m, n, count):
    """Regions of m × n matrices: orthogonal to a known x with large entries (Rust overflows; Ok exactly) and with small
    ones (the fast path); x with an entry 0, negative, 2³² − 1 or 2³² (after the overflow); random small and random
    large entries (mostly ambiguous or trivial); ±2⁵³ entries; dependent rows (a duplicate, a sum, a zero row)."""
    c = count

    def xs_(hi):
        return [rng.randint(1, hi) for _ in range(n - 1)] + [1]

    def fit(rows_flat):
        rows_ = [rows_flat[i * n:(i + 1) * n] for i in range(len(rows_flat) // n)]
        rows_ = (rows_ + [[0] * n] * m)[:m]
        return tuple(v for r in rows_ for v in r)
    big = [fit(_bal7_solved(rng, xs_(rng.choice((3, 1000, 2 ** 20, 2 ** 31 - 1))), n - 1 if m >= n - 1 else m,
                            rng.randint(0, 12), rng.choice((0.0, 0.0, 0.2, 0.4)))) for _ in range(c)]
    small = [fit(_bal7_solved(rng, xs_(rng.choice((3, 50))), min(m, n - 1), 40)) for _ in range(c)]
    special = []
    for _ in range(c):
        x = xs_(50)
        x[rng.randint(0, n - 2)] = rng.choice((0, -1, -1000, 2 ** 32 - 1, 2 ** 32, 2 ** 32 + 1))
        special.append(fit(_bal7_solved(rng, x, min(m, n - 1))))
    rsmall = [tuple(rng.randint(-r, r) for _ in range(m * n)) for r in (rng.choice((1, 3, 1000)) for _ in range(c))]
    rbig = [tuple(rng.randint(-2 ** b, 2 ** b) for _ in range(m * n)) for b in (rng.randint(30, 53) for _ in range(c))]
    e53 = []
    for _ in range(c):
        t = []
        for _ in range(min(m, n - 1)):
            while True:
                row = [rng.choice((2 ** 53, -2 ** 53, 2 ** 53 - rng.randint(1, 9), rng.randint(1, 9) - 2 ** 53,
                                   rng.randint(-9, 9))) for _ in range(n - 1)]
                if abs(sum(row)) <= 2 ** 53:
                    t += row + [-sum(row)]
                    break
        e53.append(fit(t))
    dep = []
    for _ in range(c):
        base = _bal7_solved(rng, xs_(100), max(min(m, n - 1) - 1, 1), 1)
        rows_ = [base[i * n:(i + 1) * n] for i in range(len(base) // n)]
        extra = rng.choice((rows_[0], [-v for v in rows_[-1]], [0] * n,
                            [a + b for a, b in zip(rows_[0], rows_[-1])]))
        if any(abs(v) > 2 ** 53 for v in extra):
            extra = rows_[0]
        rows_.insert(rng.randint(0, len(rows_)), extra)
        dep.append(fit([v for r in rows_ for v in r]))
    return [("known x, large entries", big), ("known x, small entries", small), ("x with 0, < 0, 2^32 - 1, 2^32", special),
            ("random small", rsmall), ("random large", rbig), ("+-2^53 entries", e53), ("dependent rows", dep)]


def bal7_checks(work, targets, which):
    """Ids 101 … 106 against the exact algorithm on every target, the multi-word path forced; True if anything
    failed."""
    failed = False
    for fid, name, m, n, count in BAL7:
        if name not in which:
            continue
        t0 = time.time()
        rng = random.Random(1039 + fid)
        sections = args_bal7(rng, m, n, count)
        xs = [x for _, sec in sections for x in sec]
        ref = [ref_bal7(list(a), m, n) for a in xs]
        want = [w for w, _ in ref]
        vw = n + 1
        first = None
        for t in targets:
            g0, g1, fl = run_probe_vec(work, t, fid, 0, xs, vw)
            a0, _, _ = run_probe_vec(work, t, fid, 1, xs, vw)
            bad = [i for i in range(len(xs)) if g0[i] != want[i]]
            bad_acc = [i for i in range(len(xs)) if a0[i] != want[i]]
            # the fast path: an overflow exactly where Rust's arithmetic overflows; elsewhere the exact outcome
            bad_fast = [i for i in range(len(xs)) if fl[i] != (0 if ref[i][1] else 1)
                        or (fl[i] == 1 and (g1[i][0] == 0) != (want[i][0] == n))
                        or (fl[i] == 1 and g1[i][0] == 0 and g1[i][1:] != want[i][1:])]
            for i in bad[:3]:
                print(f"  MISMATCH {name}{xs[i]!r} on {t}: got {g0[i]}, want {want[i]}")
            for i in bad_acc[:3]:
                print(f"  MISMATCH {name} (multi-word path forced){xs[i]!r} on {t}: got {a0[i]}, want {want[i]}")
            for i in bad_fast[:3]:
                print(f"  FAST PATH {name}{xs[i]!r} on {t}: code {g1[i][0]}, flag {fl[i]}, Rust overflows {ref[i][1]}")
            cur = (g0, g1, list(fl), a0)
            xdiff = 0 if first is None else sum(1 for u, v in zip(cur, first) for p, q in zip(u, v) if p != q)
            if first is None:
                first = cur
            slow = len(xs) - sum(1 for f in fl if f == 1)
            oks = sum(1 for w in want if w[0] == n)
            print(f"{name:7} {t:8} {len(xs):9,} matrices ({oks:,} Ok): {len(bad)} mismatches; multi-word path forced: "
                  f"{len(bad_acc)} mismatches; fast path: {len(bad_fast)} differ, Rust overflows on {slow:,}; words or "
                  f"flags differing from {targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0]:
                rates, i = [], 0
                for tag, sec in sections:
                    rates.append(f"{tag} {sum(1 for f in fl[i:i + len(sec)] if f != 1)}/{len(sec)}")
                    i += len(sec)
                print("         Rust overflows by region: " + "; ".join(rates))
        print(f"       [{time.time() - t0:.0f} s]")
    return failed


# The overflow audit's B10a (ADR 0002 #40's primitives): e^x, e^x − 1 and 10^x of an unbounded-exponent value x = s·2^e,
# correctly rounded to 53 bits at any exponent, as (s, e) pairs; the reference is MPFR at its widest exponent range
# (mpfr_ref.UX.exp, .expm1, .exp10), whose saturation (+∞, and +0 past the underflow edge) the port's matches.
B10A = [(107, "exp_ux", "exp"), (108, "expm1_ux", "expm1"), (109, "exp10_ux", "exp10")]
B10A_NAMES = tuple(name for _, name, _ in B10A)
ARITY.update({107: 2, 108: 2, 109: 2})


def _ux_pair(x):
    """(s, e) of a finite nonzero double, 1 ≤ |s| < 2 (subnormals included)."""
    m, e = math.frexp(x)
    return (m * 2.0, e - 1)


def _ux_arg(s, e):
    """The MPFR value of the pair (s, e); a zero, an infinity or NaN as itself."""
    if s == 0.0 or not math.isfinite(s):
        u = mpfr_ref.UX()
        mpfr_ref._m.mpfr_set_d(u.v, s, 0)
        return u
    return mpfr_ref.ux_scaled(s, e)


def ref_b10a(fn, a):
    """The expected (s, e) words: MPFR's correctly rounded value at its widest exponent range, as a pair."""
    s, e = getattr(_ux_arg(*a), fn)().parts()
    return [bits(s), e & W64]


def _rand_sig(rng):
    return 1.0 + rng.randrange(1 << 52) * 2.0 ** -52


def _near_ln2_multiples(rng, n, log10=False):
    """Doubles at and next to K·ln2/512 and (K + 1/2)·ln2/512 (where K's rounding flips), K of every bit length up to
    71, both signs; for 10^x the same multiples of log10(2)/512."""
    import mathx_tables
    from fractions import Fraction
    lb = 400
    l2 = mathx_tables._ln2_int(lb)                          # ⌊2^400·ln 2⌋
    l10 = mathx_tables._ln10_int(lb) if log10 else 1
    out = []
    for _ in range(n):
        b = rng.randint(1, 71)
        k = rng.randrange(1 << (b - 1), 1 << b)
        num = 2 * k + (1 if rng.random() < 0.5 else 0)     # (num/2)·ln2/512
        q = Fraction(num * l2, 1 << (lb + 10))
        if log10:
            q = q * Fraction(1 << lb, l10)
        x = ulps(float(q), rng.randint(-2, 2))
        out.append(_ux_pair(rng.choice((1, -1)) * x))
    return out


def args_b10a(rng, name):
    """Sections of (s, e) arguments: specials and extreme exponents, every binade of both signs, binary64's range,
    the tiny end, the tier edges, near multiples of ln2/512, the known hard cases, large arguments and the saturation
    edges (exp and expm1 at s·2^61, exp10 at s·2^60), and each function's own regions."""
    import mathx_tables
    lb = mathx_tables.ux_limits()
    sg = (1, -1)
    big_e = [62, 63, 64, 100, 1000, 1 << 20, 1 << 40, (1 << 62) - 2, (1 << 62) - 1, 1 << 62]
    small_e = [-55, -56, -57, -100, -1100, -5000, -(1 << 20), -(1 << 40), -(1 << 62) + 2, -(1 << 62)]
    if name != "expm1_ux":
        big_e += [(1 << 63) - 1]
        small_e += [-(1 << 63), -(1 << 62) - 1]
    sigs = [1.0, 1.5, 2.0 - 2.0 ** -52, 1.0 + 2.0 ** -52]
    special = [(0.0, 0), (-0.0, 0), (INF, 0), (-INF, 0), (NAN, 0), (0.0, 77), (-0.0, -5)]
    special += [(sv * s, e) for e in big_e + small_e for s in sigs for sv in sg]
    binade = [(rng.choice(sg) * _rand_sig(rng), e) for e in range(-1200, 66) for _ in range(36)]
    binade += [(sv * _rand_sig(rng), e) for e in (-(1 << 30), -40000, -3000, 3000, 40000, 1 << 30) for sv in sg
               for _ in range(20)]
    tiny = [(sv * s, e) for e in range(-60, -48) for s in sigs + [_rand_sig(rng) for _ in range(60)] for sv in sg]
    if name == "exp10_ux":
        lo, hi, d1, d2, sh = -323.7, 308.3, 307.0, 512.0, 60
        ov, un = lb["EXP10UX_OVF"], lb["EXP10UX_UNF"]
    else:
        lo, hi, d1, d2, sh = -745.2, 709.8, 708.0, 1024.0, 61
        ov, un = lb["EXPUX_OVF"], lb["EXPUX_UNF"]
        if name == "expm1_ux":
            lo = -40.0
    b64 = [_ux_pair(rng.uniform(lo, hi)) for _ in range(30_000)]
    b64 += [_ux_pair(x) for x in SPECIALS if x != 0.0 and math.isfinite(x)]
    edges = []
    for c in (d1, d2, 745.1332191019411, 709.782712893384, 323.60724533877976, 308.2547155599167, 64.0,
              37.42994775023704, 2.0 ** 40, 2.0 ** 41):                # the last two: the fast reduction's end
        for x in around(c, 150):
            edges += [_ux_pair(x), _ux_pair(-x)]
    edges += [_ux_pair(rng.uniform(d1, d2) * rng.choice(sg)) for _ in range(8_000)]
    red = _near_ln2_multiples(rng, 30_000, name == "exp10_ux")
    hard = [_ux_pair(h) for h in HARD_EXP] + [_ux_pair(-h) for h in HARD_EXP]
    large = [_ux_pair(rng.choice(sg) * 10.0 ** rng.uniform(3.0, 18.5)) for _ in range(15_000)]
    large += [_ux_pair(x * sv) for c in (1e4, 1e12) for x in around(c, 300) for sv in sg]
    sat = []
    for c, sv in ((ov, 1), (un, -1)):
        for d in range(-60, 61):
            sat.append((sv * ulps(c, d), sh))
    sat += [(sv * s, e) for e in (sh - 1, sh, sh + 1) for s in sigs + [_rand_sig(rng) for _ in range(30)] for sv in sg]
    sections = [("special", special), ("every binade", binade), ("binary64 range", b64), ("tiny", tiny),
                ("tier edges", edges), ("near K·ln2/512", red), ("hard cases", hard), ("large", large),
                ("saturation", sat)]
    if name == "expm1_ux":
        minus1 = [_ux_pair(rng.uniform(-40.0, 40.0)) for _ in range(20_000)]
        minus1 += [_ux_pair(x) for c in (-37.42994775023704, -0.5, 0.5, 2.0 ** -6, -2.0 ** -6, 83.17766166719343)
                   for x in around(c, 300)]
        minus1 += [_ux_pair(rng.uniform(60.0, 100.0)) for _ in range(5_000)]           # E around 120: the −1 at 2^-120
        sections.append(("the −1 region", minus1))
    if name == "exp10_ux":
        ints = [_ux_pair(float(n)) for n in range(-400, 401) if n != 0]
        ints += [_ux_pair(ulps(float(n), d)) for n in range(-30, 31) if n != 0 for d in (-3, -2, -1, 1, 2, 3)]
        ints += [_ux_pair(ulps(23.0, d)) for d in range(-40, 41)]
        ints += [_ux_pair(float(n)) for n in (10 ** 3, 10 ** 6, 2 ** 40, 10 ** 15, 10 ** 18, -10 ** 18)]
        sections.append(("integers", ints))
    return sections


def b10a_checks(work, targets, which):
    """Ids 107 … 109 against MPFR at its widest exponent range on every target, the accurate path forced; True if
    anything failed. The flag's bit 2 marks a normal double x whose binary64 function has a normal result other than
    the pair (it must never be set)."""
    failed = False
    for fid, name, fn in B10A:
        if name not in which:
            continue
        t0 = time.time()
        rng = random.Random(1040 + fid)
        sections = args_b10a(rng, name)
        xs = [x for _, sec in sections for x in sec]
        want = [ref_b10a(fn, a) for a in xs]
        # the references' own consistency: MPFR's 53-bit value is its binary64 one wherever that is normal
        nref = 0
        for a, w in zip(xs, want):
            s, e = a
            if s == 0.0 or not math.isfinite(s) or not -1022 <= e <= 1023:
                continue
            y = getattr(mpfr_ref, fn)(math.ldexp(s, e))
            if math.isfinite(y) and abs(y) >= 2.0 ** -1022:
                nref += 1
                if [bits(_ux_pair(y)[0]), _ux_pair(y)[1] & W64] != w:
                    print(f"  REFERENCE {name}{a!r}: MPFR's binary64 {fn} {y!r} differs from its 53-bit value")
                    failed = True

        def eq(g, w):
            return g == w or (math.isnan(fb(g[0])) and math.isnan(fb(w[0])))
        first = None
        for t in targets:
            g0, g1, fl = run_probe_vec(work, t, fid, 0, xs, 2)
            a0, _, _ = run_probe_vec(work, t, fid, 1, xs, 2)
            bad = [i for i in range(len(xs)) if not eq(g0[i], want[i])]
            bad_acc = [i for i in range(len(xs)) if not eq(a0[i], want[i])]
            bad_fast = [i for i in range(len(xs)) if fl[i] & 1 and not eq(g1[i], g0[i])]
            bad_d = [i for i in range(len(xs)) if fl[i] & 2]
            for i in bad[:5]:
                print(f"  MISMATCH {name}{xs[i]!r} on {t}: got {g0[i]}, want {want[i]}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {name} (accurate path forced){xs[i]!r} on {t}: got {a0[i]}, want {want[i]}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {name}{xs[i]!r} on {t}: {g1[i]} vs {g0[i]}")
            for i in bad_d[:5]:
                print(f"  BINARY64 DIFFERS {name}{xs[i]!r} on {t}: the pair {g0[i]}")

            def canon(ws):
                return [[NANW, 0] if math.isnan(fb(w[0])) else w for w in ws]
            cur = (canon(g0), canon(g1), list(fl), canon(a0))
            xdiff = 0 if first is None else sum(1 for u, v in zip(cur, first) for p, q in zip(u, v) if p != q)
            if first is None:
                first = cur
            slow = sum(1 for f in fl if not f & 1)
            print(f"{name:8} {t:8} {len(xs):9,} args: {len(bad)} mismatches; accurate path forced: {len(bad_acc)} "
                  f"mismatches; fast path: {len(bad_fast)} differ, fell back {slow:,}; binary64 function differs: "
                  f"{len(bad_d)}; bits or flags differing from {targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or bad_d or xdiff)
            if t == targets[0]:
                rates, i = [], 0
                for tag, sec in sections:
                    rates.append(f"{tag} {sum(1 for f in fl[i:i + len(sec)] if not f & 1)}/{len(sec)}")
                    i += len(sec)
                print("         fallbacks by region: " + "; ".join(rates))
                print(f"         MPFR's binary64 {fn} equals its 53-bit value on all {nref:,} arguments with a normal "
                      f"binary64 result")
        print(f"       [{time.time() - t0:.0f} s]")
    return failed


def flame6_bisection_check(work, targets, which):
    """Id 87's second slot, the transcribed bisection over [T₀, T_max] with ADR 0002 #38's exact same-sign test, against
    its reference on every target; True if anything failed."""
    if "flame" not in which:
        return False
    failed = False
    rng = random.Random(1033 + 87)
    xs = [x for _, sec in args_flame6(rng) for x in sec]
    want = [ref_flame6(a)[1] for a in xs]
    for t in targets:
        _, g1, _ = run_probe(work, t, 87, 0, xs)
        bad = [i for i in range(len(xs)) if g1[i] != want[i]]
        for i in bad[:5]:
            print(f"  MISMATCH flame (bisection){xs[i]!r} on {t}: got {g1[i]:016x}, want {want[i]:016x}")
        refused = sum(1 for w in want if w == (1 << 63) | _err_code("KM_ERR_INVALID_INPUT"))
        print(f"flame   {t:8} {len(xs):9,} args: the bisection's exact same-sign test: {len(bad)} mismatches "
              f"({refused:,} refused)")
        failed |= bool(bad)
    return failed


def audit_checks(work, targets, which):
    """Ids 33 … 82 against their exact references on every target, the exact path forced; True if anything failed."""
    failed = False
    gens = {"mono": args_mono, "heat": args_mul3, "shc": args_shc, "boil": args_mul3, "frz": args_mul3, "osm": args_osm,
            "zhl": args_zhl, "ssi": args_ssi, "pre": args_pre, "soc": args_soc, "mm": args_mm, "monotie": args_monotie,
            "sqrtdiv": args_sqrtdiv, "dil": args_dil, "cond": args_cond, "fmass": args_fmass, "fcharge": args_fcharge,
            "fratio": args_fratio, "nernst": args_nernst, "abs": args_mul3_30, "conc": args_conc,
            "vib": lambda r: args_sqrtdiv(r, 60_000), "rotor": args_rotor, "rot": args_rot,
            "moment": lambda r: args_mul3_30(r, True),
            "wsqrt": args_wsqrt, "igp": args_ideal, "igv": args_ideal, "dens": args_dens, "vdw": args_vdw,
            "prp": args_prp, "rkp": args_rkp, "lj": args_lj, "width": args_width, "coul": args_coul, "force": args_coul,
            "lande": args_lande, "bea": args_bea,
            "hess": args_hess, "eff": args_eff, "rsum": args_rsum, "qv": args_qv, "ionic": args_ionic, "haber": args_haber,
            "gibbs": args_gibbs, "kohl": args_kohl, "lim": args_lim, "zoc": args_zoc, "shcp": args_shomate,
            "shh": args_shomate,
            "keq": args_keq, "dhl": lambda r: args_dh(r, "dhl"), "dhe": lambda r: args_dh(r, "dhe"),
            "dav": lambda r: args_dh(r, "dav"),
            "flame": args_flame6, "pk3": lambda r: args_peaks6(r, 3), "pk5": lambda r: args_peaks6(r, 5),
            "pk8": lambda r: args_peaks6(r, 8), "pk16": lambda r: args_peaks6(r, 16), "signsc": args_signsc}

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
include "src/reaction.cyr"
include "src/stoichiometry.cyr"
include "src/thermo.cyr"
include "src/thermochem_data.cyr"
include "src/thermochem.cyr"
include "src/phase_data.cyr"
include "src/phase.cyr"
include "src/gas_data.cyr"
include "src/gas.cyr"
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

# The overflow audit's B7 (ADR 0002 #39): the (m, n) of ids 101 … 106 (cr_verify.py's BAL7).
#must_use
fn bal7_dims(fid): (i64, i64) {
    if (fid == 101) { return (2, 3); }
    if (fid == 102) { return (3, 4); }
    if (fid == 103) { return (5, 6); }
    if (fid == 104) { return (8, 9); }
    if (fid == 105) { return (12, 13); }
    return (24, 25);
}

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
    if (fid == 58 || fid == 64) { return 4; }
    if (fid == 62 || fid == 63 || fid == 69) { return 5; }
    if (fid >= 59 && fid <= 70) { return 3; }
    if (fid == 71 || fid == 81 || fid == 82) { return 6; }
    if (fid == 72 || fid == 74) { return 8; }
    if (fid == 73) { return 11; }
    if (fid == 75) { return 7; }
    if (fid == 76) { return 5; }
    if (fid == 79) { return 4; }
    if (fid >= 77 && fid <= 80) { return 3; }
    if (fid == 83) { return 2; }
    if (fid == 84 || fid == 86) { return 3; }
    if (fid == 85) { return 5; }
    if (fid == 87) { return 8; }
    if (fid == 88) { return 4; }
    if (fid == 89) { return 6; }
    if (fid == 90) { return 9; }
    if (fid == 91) { return 17; }
    if (fid >= 92 && fid <= 98) { return 1 << (fid - 92); }
    if (fid == 99) { return 256; }
    if (fid == 100) { return 17; }
    if (fid >= 101 && fid <= 106) {
        var m7, n7 = bal7_dims(fid);
        return m7 * n7;
    }
    if (fid >= 107 && fid <= 109) { return 2; }
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

# The overflow audit's B3 (ADR 0002 #31, #32): ids 58 … 70 as result words (cr_verify.py's AUDIT): mode 0 the public
# function (an Err as the NaN word 0x7FF8000000000001), its fast path and flag (the shape's straight-line monomial path,
# _km_mul_sqrt_ratio_fast, or the Rust-chain fast path of a #20-style function where its range test passes; for
# binding_energy_per_nucleon the inline test, Δm·931.494… finite); mode 1 the slow path forced (_km_cr_mono_acc at the
# factors, _km_mul_sqrt_ratio_acc, the _km_*_ux evaluation). Id 58 is _km_cr_mul_sqrt_ratio itself (its NaN outside the
# domain in both slots). A charge, a count or wsqrt's e is an integer word. An Err row, and a non-finite charge (Rust's
# expression), return the public result in both slots and modes.
#must_use
fn run_audit3(fid, mode, a): (i64, i64, i64) {
    var nanw = 0x7FF8000000000001;
    var x: f64 = load64(a);
    var y: f64 = load64(a + 8);
    var z: f64 = load64(a + 16);
    var w: f64 = load64(a + 24);
    var u: f64 = load64(a + 32);
    var rg: f64 = KM_GAS_CONSTANT;
    var kc: f64 = KM_COULOMB_CONSTANT;
    if (fid == 58 || fid == 66) {
        var e = 0 - 1;
        var pw = 0;
        if (fid == 58) {
            e = load64(a + 24);
            var r58: f64 = _km_cr_mul_sqrt_ratio(x, y, z, e);
            if (_km_is_nan(r58) == 1) { return (nanw, nanw, 1); }
            pw = _km_bits(r58);
        } else {
            var t66, v66 = morse_width_parameter(x, y, z);
            if (is_ok(t66) != 1) { return (nanw, nanw, 1); }
            pw = v66;
        }
        if (mode == 1) {
            var sa: f64 = _km_mul_sqrt_ratio_acc(x, y, z, e);
            return (_km_bits(sa), _km_bits(sa), 0);
        }
        var sf, so = _km_mul_sqrt_ratio_fast(x, y, z, e);
        return (pw, _km_bits(sf), so);
    }
    if (fid >= 62 && fid <= 64) {
        var tp = 0;
        var vp = 0;
        if (fid == 62) {
            var t62, v62 = van_der_waals_pressure(x, y, z, w, u);
            tp = t62;
            vp = v62;
        } elif (fid == 63) {
            var pr = KmPengRobinsonParams { tc: x, pc: y, omega: z };
            var t63, v63 = peng_robinson_pressure(pr, w, u);
            tp = t63;
            vp = v63;
        } else {
            var t64, v64 = redlich_kwong_pressure(x, y, z, w);
            tp = t64;
            vp = v64;
        }
        var pubw = nanw;
        if (is_ok(tp) == 1) { pubw = vp; }
        # Rust's validation before the co-volume test: an Err there is the public Err in both slots
        var pre = 1;
        if (fid == 62) {
            if (_km_is_positive(y) != 1 || _km_is_finite(x) != 1 || _km_is_non_negative(w) != 1) { pre = 0; }
            if (_km_is_non_negative(u) != 1 || _km_is_positive(z) != 1) { pre = 0; }
        } elif (fid == 63) {
            if (_km_is_positive(w) != 1 || _km_is_positive(x) != 1 || _km_is_positive(y) != 1) { pre = 0; }
            if (_km_is_finite(z) != 1 || _km_is_finite(u) != 1) { pre = 0; }
        } else {
            if (_km_is_positive(z) != 1 || _km_is_positive(x) != 1 || _km_is_positive(y) != 1) { pre = 0; }
            if (_km_is_finite(w) != 1) { pre = 0; }
        }
        if (pre != 1) { return (pubw, pubw, 1); }
        var sp: f64 = 0.0;
        var se = 0;
        if (mode == 1) {
            if (fid == 62) {
                var p62, e62 = _km_vdw_ux(x, y, z, w, u);
                sp = p62;
                se = e62;
            } elif (fid == 63) {
                var p63, e63 = _km_pr_pressure_ux(x, y, z, w, u);
                sp = p63;
                se = e63;
            } else {
                var p64, e64 = _km_rk_pressure_ux(x, y, z, w);
                sp = p64;
                se = e64;
            }
            var sw = _km_bits(sp);
            if (se != KM_ERR_NONE) { sw = nanw; }
            return (sw, sw, 0);
        }
        var fp: f64 = 0.0;
        var fe = 0;
        var fo = 0;
        if (fid == 62) {
            var q62, f62, o62 = _km_vdw_fast(x, y, z, w, u);
            fp = q62;
            fe = f62;
            fo = o62;
        } elif (fid == 63) {
            var q63, f63, o63 = _km_pr_pressure_fast(x, y, z, w, u);
            fp = q63;
            fe = f63;
            fo = o63;
        } else {
            var q64, f64e, o64 = _km_rk_pressure_fast(x, y, z, w);
            fp = q64;
            fe = f64e;
            fo = o64;
        }
        var fw = _km_bits(fp);
        if (fe != KM_ERR_NONE) { fw = nanw; }
        return (pubw, fw, fo);
    }
    if (fid == 65) {
        var t65, v65 = lennard_jones(x, y, z);
        if (is_ok(t65) != 1) { return (nanw, nanw, 1); }
        if (mode == 1) {
            var l65: f64 = _km_lj_ux(x, y, z);
            return (_km_bits(l65), _km_bits(l65), 0);
        }
        var f65, o65 = _km_lj_fast(x, y, z);
        return (v65, _km_bits(f65), o65);
    }
    if (fid == 69) {
        var zc = load64(a + 8);
        var za = load64(a + 16);
        var t69, v69 = born_lande_lattice_energy(x, zc, za, w, u);
        if (is_ok(t69) != 1) { return (nanw, nanw, 1); }
        var zcf: f64 = f64_from(zc);
        var zaf: f64 = f64_from(za);
        if (mode == 1) {
            var l69: f64 = _km_born_lande_ux(x, zcf, zaf, w, u);
            return (_km_bits(l69), _km_bits(l69), 0);
        }
        var f69, o69 = _km_born_lande_fast(x, zcf, zaf, w, u);
        return (v69, _km_bits(f69), o69);
    }
    if (fid == 70) {
        var zz = load64(a);
        var nn = load64(a + 8);
        var t70, v70 = binding_energy_per_nucleon(zz, nn, z);
        if (is_ok(t70) != 1) { return (nanw, nanw, 1); }
        var zf: f64 = f64_from(zz);
        var nf: f64 = f64_from(nn);
        var af: f64 = f64_from(zz + nn);
        var md: f64 = _km_mass_defect(zf, nf, z);
        if (mode == 1) {
            var l70: f64 = _km_bea_ux(md, af);
            return (_km_bits(l70), _km_bits(l70), 0);
        }
        var be: f64 = md * KM_AMU_MEV;
        var dbe: f64 = be - be;
        if (dbe == 0.0) {
            var q70: f64 = be / af;
            return (v70, _km_bits(q70), 1);
        }
        return (v70, 0, 0);
    }
    # the monomials: 59 igp, 60 igv, 61 dens, 67 coul, 68 force
    var nb: f64[4];
    var db: f64[4];
    var nm = 3;
    var nd = 1;
    var pw2 = 0;
    var ff: f64 = 0.0;
    var fo2 = 0;
    if (fid == 59 || fid == 60) {
        var t59 = 0;
        var v59 = 0;
        if (fid == 59) {
            var ta, va = ideal_gas_pressure(x, y, z);
            t59 = ta;
            v59 = va;
        } else {
            var tb, vb = ideal_gas_volume(x, y, z);
            t59 = tb;
            v59 = vb;
        }
        if (is_ok(t59) != 1) { return (nanw, nanw, 1); }
        pw2 = v59;
        store64(&nb, x);
        store64(&nb + 8, rg);
        store64(&nb + 16, y);
        store64(&db, z);
        var f59, o59 = _km_mul3_div_fast(x, rg, y, z);
        ff = f59;
        fo2 = o59;
    } elif (fid == 61) {
        var t61, v61 = gas_density(x, y, z);
        if (is_ok(t61) != 1) { return (nanw, nanw, 1); }
        pw2 = v61;
        store64(&nb, x);
        store64(&nb + 8, y);
        store64(&db, rg);
        store64(&db + 8, z);
        nm = 2;
        nd = 2;
        var f61, o61 = _km_mul2_div2_fast(x, y, rg, z);
        ff = f61;
        fo2 = o61;
    } else {
        var t67 = 0;
        var v67 = 0;
        if (fid == 67) {
            var tc1, vc1 = coulomb(x, y, z);
            t67 = tc1;
            v67 = vc1;
        } else {
            var tc2, vc2 = coulomb_force(x, y, z);
            t67 = tc2;
            v67 = vc2;
        }
        if (is_ok(t67) != 1) { return (nanw, nanw, 1); }
        pw2 = v67;
        if (_km_is_finite(x) != 1 || _km_is_finite(y) != 1) { return (pw2, pw2, 1); }
        store64(&nb, kc);
        store64(&nb + 8, x);
        store64(&nb + 16, y);
        store64(&db, z);
        if (fid == 67) {
            var f67, o67 = _km_mul3_div_fast(kc, x, y, z);
            ff = f67;
            fo2 = o67;
        } else {
            store64(&db + 8, z);
            nd = 2;
            var f68, o68 = _km_mul3_div2_fast(kc, x, y, z, z);
            ff = f68;
            fo2 = o68;
        }
    }
    if (mode == 1) {
        var aa: f64 = _km_cr_mono_acc(&nb, nm, &db, nd);
        return (_km_bits(aa), _km_bits(aa), 0);
    }
    return (pw2, _km_bits(ff), fo2);
}

# The overflow audit's B4 (ADR 0002 #33, #34): ids 71 … 82 (cr_verify.py's AUDIT4_NAMES has the argument layouts): mode 0
# the public function (an Err as the NaN word), Rust's chain where the function's fast test passes and the flag (the
# pure sums' "finite", _km_reaction_property's flag, the inline normal-range tests, _km_kohlrausch_fast,
# _km_limiting_fast, _km_shomate_cp_fast, _km_shomate_h_fast); mode 1 the unbounded-exponent slow path forced (the
# _km_*_ux functions). A non-finite input keeps Rust's expression: the public result in both slots and modes.
#must_use
fn run_audit4(fid, mode, a): (i64, i64, i64) {
    var nanw = 0x7FF8000000000001;
    if (fid == 71) {
        var n71 = load64(a);
        var r71: f64 = hess_law(a + 8, n71);
        if (_km_finite_slice(a + 8, n71) != 1) { return (_km_bits(r71), _km_bits(r71), 1); }
        if (mode == 1) {
            var u71: f64 = _km_sum_ux_f64(a + 8, n71);
            return (_km_bits(u71), _km_bits(u71), 0);
        }
        var s71: f64 = _km_sum(a + 8, n71);
        var d71: f64 = s71 - s71;
        if (d71 == 0.0) { return (_km_bits(r71), _km_bits(s71), 1); }
        return (_km_bits(r71), 0, 0);
    }
    if (fid == 72 || fid == 74) {
        var np = load64(a);
        var nr = load64(a + 8);
        var pp = a + 16;
        var rp = a + 40;
        var r72: f64 = 0.0;
        if (fid == 72) { r72 = enthalpy_from_formation(pp, np, rp, nr); } else { r72 = q_value(pp, np, rp, nr); }
        if (_km_finite_slice(pp, np) != 1 || _km_finite_slice(rp, nr) != 1) {
            return (_km_bits(r72), _km_bits(r72), 1);
        }
        if (mode == 1) {
            var u72: f64 = 0.0;
            if (fid == 72) { u72 = _km_formation_ux(pp, np, rp, nr); } else { u72 = _km_q_value_ux(pp, np, rp, nr); }
            return (_km_bits(u72), _km_bits(u72), 0);
        }
        var sp: f64 = _km_sum(pp, np);
        var sr: f64 = _km_sum(rp, nr);
        var d72: f64 = sp - sr;
        var dd72: f64 = d72 - d72;
        if (fid == 72) {
            var dp: f64 = sp - sp;
            var dr: f64 = sr - sr;
            if (dp == 0.0 && dr == 0.0) { return (_km_bits(r72), _km_bits(d72), 1); }
            return (_km_bits(r72), 0, 0);
        }
        if (dd72 == 0.0) {
            var q74: f64 = d72 * KM_AMU_MEV;
            return (_km_bits(r72), _km_bits(q74), 1);
        }
        return (_km_bits(r72), 0, 0);
    }
    if (fid == 73) {
        var col = load64(a);
        var np3 = load64(a + 8);
        var nr3 = load64(a + 16);
        var sp3[32];
        var sr3[32];
        var j = 0;
        while (j < 2) {
            var fp = load64(a + 24 + 16 * j);
            var fr = load64(a + 56 + 16 * j);
            var np_ = "XeF6(g)";
            var nr_ = "XeF6(g)";
            if (fp >= 0 && fp < KM_THERMOCHEM_COUNT) { np_ = thermochem_data_formula(fp + 1); }
            if (fr >= 0 && fr < KM_THERMOCHEM_COUNT) { nr_ = thermochem_data_formula(fr + 1); }
            store64(&sp3 + 16 * j, np_);
            store64(&sp3 + 16 * j + 8, load64(a + 32 + 16 * j));
            store64(&sr3 + 16 * j, nr_);
            store64(&sr3 + 16 * j + 8, load64(a + 64 + 16 * j));
            j = j + 1;
        }
        var colp = &KM_THERMOCHEM_DELTA_HF_KJ;
        if (col == 1) { colp = &KM_THERMOCHEM_DELTA_GF_KJ; }
        if (col == 2) { colp = &KM_THERMOCHEM_S_STANDARD_J; }
        var t73 = 0;
        var v73 = 0;
        if (col == 0) {
            var t0, v0 = reaction_enthalpy(&sp3, np3, &sr3, nr3);
            t73 = t0;
            v73 = v0;
        } elif (col == 1) {
            var t1, v1 = reaction_gibbs_energy(&sp3, np3, &sr3, nr3);
            t73 = t1;
            v73 = v1;
        } else {
            var t2, v2 = reaction_entropy(&sp3, np3, &sr3, nr3);
            t73 = t2;
            v73 = v2;
        }
        if (is_ok(t73) != 1) { return (nanw, nanw, 1); }
        var e73, d73, f73 = _km_reaction_property(&sp3, np3, &sr3, nr3, colp);
        if (f73 == 2) { return (v73, v73, 1); }
        if (mode == 1) {
            var u73: f64 = 0.0;
            if (col == 2) { u73 = _km_reaction_s_ux(&sp3, np3, &sr3, nr3); } else {
                u73 = _km_reaction_j_ux(&sp3, np3, &sr3, nr3, colp);
            }
            return (_km_bits(u73), _km_bits(u73), 0);
        }
        var dd73: f64 = d73 - d73;
        if (col == 2) {
            if (f73 == 1) { return (v73, _km_bits(d73), 1); }
            return (v73, 0, 0);
        }
        if (f73 == 1 && dd73 == 0.0) {
            var j73: f64 = d73 * 1000.0;
            return (v73, _km_bits(j73), 1);
        }
        return (v73, 0, 0);
    }
    if (fid == 75) {
        var n75 = load64(a);
        var ions[48];
        var k = 0;
        var fin = 1;
        while (k < n75) {
            var c75: f64 = load64(a + 8 + 16 * k);
            store64(&ions + 16 * k, load64(a + 8 + 16 * k));
            store64(&ions + 16 * k + 8, load64(a + 16 + 16 * k));
            if (_km_is_finite(c75) != 1) { fin = 0; }
            k = k + 1;
        }
        var r75: f64 = ionic_strength(&ions, n75);
        if (fin != 1) { return (_km_bits(r75), _km_bits(r75), 1); }
        if (mode == 1) {
            var u75: f64 = _km_ionic_strength_ux(&ions, n75);
            return (_km_bits(u75), _km_bits(u75), 0);
        }
        var s75: f64 = -0.0;
        k = 0;
        while (k < n75) {
            var c: f64 = load64(&ions + 16 * k);
            var z: f64 = f64_from(load64(&ions + 16 * k + 8));
            var term: f64 = (c * z) * z;
            s75 = s75 + term;
            k = k + 1;
        }
        var ds75: f64 = s75 - s75;
        if (ds75 == 0.0) {
            var i75: f64 = 0.5 * s75;
            return (_km_bits(r75), _km_bits(i75), 1);
        }
        return (_km_bits(r75), 0, 0);
    }
    var x: f64 = load64(a);
    var y: f64 = load64(a + 8);
    var z3: f64 = load64(a + 16);
    var w: f64 = load64(a + 24);
    var v: f64 = load64(a + 32);
    var t6: f64 = load64(a + 40);
    if (fid == 76) {
        var r76: f64 = born_haber_lattice_energy(x, y, z3, w, v);
        if (_km_finite6(x, y, z3, w, v, 1.0) != 1) { return (_km_bits(r76), _km_bits(r76), 1); }
        if (mode == 1) {
            var u76: f64 = _km_born_haber_ux(x, y, z3, w, v);
            return (_km_bits(u76), _km_bits(u76), 0);
        }
        var u: f64 = x - y;
        u = u - z3;
        var half: f64 = w / 2.0;
        u = u - half;
        u = u - v;
        var du: f64 = u - u;
        var h2: f64 = half * 2.0;
        if (du == 0.0 && h2 == w) { return (_km_bits(r76), _km_bits(u), 1); }
        return (_km_bits(r76), 0, 0);
    }
    if (fid == 77) {
        var t77, v77 = gibbs_free_energy(x, y, z3);
        if (is_ok(t77) != 1) { return (nanw, nanw, 1); }
        if (mode == 1) {
            var u77: f64 = _km_gibbs_ux(x, y, z3);
            return (_km_bits(u77), _km_bits(u77), 0);
        }
        var p77: f64 = y * z3;
        if (_km_ux_normal(p77) == 1 || z3 == 0.0) {
            var g77: f64 = x - p77;
            return (v77, _km_bits(g77), 1);
        }
        return (v77, 0, 0);
    }
    if (fid == 78) {
        var t78, v78 = kohlrausch_molar_conductivity(x, y, z3);
        if (is_ok(t78) != 1) { return (nanw, nanw, 1); }
        if (mode == 1) {
            var u78: f64 = _km_kohlrausch_ux(x, y, z3);
            return (_km_bits(u78), _km_bits(u78), 0);
        }
        var f78, o78 = _km_kohlrausch_fast(x, y, z3);
        return (v78, _km_bits(f78), o78);
    }
    if (fid == 79) {
        var r79: f64 = limiting_molar_conductivity(x, y, z3, w);
        if (_km_finite6(x, y, z3, w, 1.0, 1.0) != 1) { return (_km_bits(r79), _km_bits(r79), 1); }
        if (mode == 1) {
            var u79: f64 = _km_limiting_ux(x, y, z3, w);
            return (_km_bits(u79), _km_bits(u79), 0);
        }
        var f79, o79 = _km_limiting_fast(x, y, z3, w);
        return (_km_bits(r79), _km_bits(f79), o79);
    }
    if (fid == 80) {
        var r80: f64 = zero_order_concentration(x, y, z3);
        if (_km_finite6(x, y, z3, 1.0, 1.0, 1.0) != 1) { return (_km_bits(r80), _km_bits(r80), 1); }
        if (mode == 1) {
            var u80: f64 = _km_zero_order_ux(x, y, z3);
            return (_km_bits(u80), _km_bits(u80), 0);
        }
        var p80: f64 = y * z3;
        if (_km_ux_normal(p80) == 1 || y == 0.0 || z3 == 0.0) {
            var rem: f64 = x - p80;
            if (rem < 0.0) { rem = 0.0; }
            return (_km_bits(r80), _km_bits(rem), 1);
        }
        return (_km_bits(r80), 0, 0);
    }
    var sh = KmShomate { a: x, b: y, c: z3, d: w, e: v, t_min: 298.0, t_max: 6000.0 };
    var r81: f64 = 0.0;
    if (fid == 81) { r81 = sh.cp(t6); } else { r81 = sh.enthalpy_antiderivative(t6); }
    if (t6 == 0.0 || _km_finite6(x, y, z3, w, v, t6) != 1) { return (_km_bits(r81), _km_bits(r81), 1); }
    if (mode == 1) {
        var u81: f64 = 0.0;
        if (fid == 81) { u81 = _km_shomate_cp_ux(x, y, z3, w, v, t6); } else { u81 = _km_shomate_h_ux(x, y, z3, w, v, t6); }
        return (_km_bits(u81), _km_bits(u81), 0);
    }
    var f81: f64 = 0.0;
    var o81 = 0;
    if (fid == 81) {
        var fa, oa = _km_shomate_cp_fast(x, y, z3, w, v, t6);
        f81 = fa;
        o81 = oa;
    } else {
        var fb, ob = _km_shomate_h_fast(x, y, z3, w, v, t6);
        f81 = fb;
        o81 = ob;
    }
    return (_km_bits(r81), _km_bits(f81), o81);
}


# The overflow audit's B5 (ADR 0002 #35): an Err as the word 2⁶³ | its KM_ERR_* code (K and γ are never negative).
#must_use
fn _km_res5(t, v): i64 {
    if (is_ok(t) == 1) { return v; }
    return (1 << 63) | v;
}

# equilibrium_constant's result word for a K: Rust's ComputationError where it is not finite.
#must_use
fn _km_k5(k: f64): i64 {
    if (_km_is_finite(k) != 1) { return (1 << 63) | KM_ERR_COMPUTATION_ERROR; }
    return _km_bits(k);
}

# The overflow audit's B5 (ADR 0002 #35): ids 83 … 86 (cr_verify.py's AUDIT5_NAMES has the argument layouts): mode 0
# the public function, Rust's chain where the function's fast test passes (R·T in [2⁻¹⁰²¹, DBL_MAX]; A·z² finite; the
# extended law's numerator and (B·å)·√I finite) and the flag; mode 1 the unbounded-exponent argument forced
# (_km_keq_arg_ux, _km_dh_limiting_arg_ux, _km_dh_extended_arg_ux, _km_davies_arg_ux), then the correctly rounded exp or
# 10^x. A row Rust's validation refuses returns the public Err in both slots and modes.
#must_use
fn run_audit5(fid, mode, a): (i64, i64, i64) {
    if (fid == 83) {
        var dg: f64 = load64(a);
        var tk: f64 = load64(a + 8);
        var t83, v83 = equilibrium_constant(dg, tk);
        var w83 = _km_res5(t83, v83);
        if (_km_is_positive(tk) != 1 || _km_is_finite(dg) != 1) { return (w83, w83, 1); }
        if (mode == 1) {
            var u83: f64 = _km_keq_arg_ux(dg, tk);
            var k83: f64 = _km_cr_exp(u83);
            var s83 = _km_k5(k83);
            return (s83, s83, 0);
        }
        var r: f64 = KM_GAS_CONSTANT;
        var rt: f64 = r * tk;
        var lo: f64 = 0x0020000000000000;
        var hi: f64 = 0x7FEFFFFFFFFFFFFF;
        if (rt >= lo && rt <= hi) {
            var x83: f64 = -dg / rt;
            var f83: f64 = _km_cr_exp(x83);
            return (w83, _km_k5(f83), 1);
        }
        return (w83, 0, 0);
    }
    var z = load64(a);
    var ii: f64 = load64(a + 8);
    var aa: f64 = load64(a + 16);
    var bc: f64 = load64(a + 24);
    var sz: f64 = load64(a + 32);
    var tp = 0;
    var vp = 0;
    if (fid == 84) {
        var t84, v84 = debye_huckel_limiting(z, ii, aa);
        tp = t84;
        vp = v84;
    } elif (fid == 85) {
        var t85, v85 = debye_huckel_extended(z, ii, aa, bc, sz);
        tp = t85;
        vp = v85;
    } else {
        var t86, v86 = davies_activity(z, ii, aa);
        tp = t86;
        vp = v86;
    }
    var w = _km_res5(tp, vp);
    if (is_ok(tp) != 1) { return (w, w, 1); }
    if (mode == 1) {
        var u: f64 = 0.0;
        if (fid == 84) {
            u = _km_dh_limiting_arg_ux(z, ii, aa);
        } elif (fid == 85) {
            u = _km_dh_extended_arg_ux(z, ii, aa, bc, sz);
        } else {
            u = _km_davies_arg_ux(z, ii, aa);
        }
        var g: f64 = _km_cr_exp10(u);
        return (_km_bits(g), _km_bits(g), 0);
    }
    var zf: f64 = f64_from(z);
    var sq: f64 = f64_sqrt(ii);
    var p: f64 = ((-aa) * zf) * zf;
    var dp: f64 = p - p;
    if (fid == 85) {
        var num: f64 = p * sq;
        var yv: f64 = (bc * sz) * sq;
        var dn: f64 = num - num;
        var dy: f64 = yv - yv;
        if (dn == 0.0 && dy == 0.0) {
            var l85: f64 = num / (1.0 + yv);
            var g85: f64 = _km_cr_exp10(l85);
            return (w, _km_bits(g85), 1);
        }
        return (w, 0, 0);
    }
    if (dp != 0.0) { return (w, 0, 0); }
    var lg: f64 = p * sq;
    if (fid == 86) {
        var br: f64 = (sq / (1.0 + sq)) - (0.3 * ii);
        lg = p * br;
    }
    var gf: f64 = _km_cr_exp10(lg);
    return (w, _km_bits(gf), 1);
}

# The overflow audit's B6 (ADR 0002 #38): the flame family's product formulas by index (cr_verify.py's FLAME6_SPECIES).
#must_use
fn _km_species6(k): i64 {
    if (k == 0) { return "CO2(g)"; }
    if (k == 1) { return "H2O(g)"; }
    if (k == 2) { return "N2(g)"; }
    if (k == 3) { return "CO(g)"; }
    if (k == 4) { return "O2(g)"; }
    if (k == 5) { return "H2(g)"; }
    return "Ar(g)";
}

# The overflow audit's B6 (ADR 0002 #37, #38): ids 87 … 91 as result words (cr_verify.py's AUDIT6_NAMES). Id 87, flame
# (adiabatic_flame_temperature: ΔH, T₀, then three (species index, moles) pairs, an index below 0 unused): the public
# result (an Err as 2⁶³ | its code) and, where Rust's validation passes (Ok or the ComputationError), the transcribed
# bisection over [T₀, T_max] in the second slot (its exact same-sign test), else the public word; the flag 0 (no fast
# path: cr_verify.py checks the second slot on its own); both modes.
# Ids 88 … 91, find_peaks of 3, 5, 8 and 16 values (then the threshold): the peaks as a bit mask (bit i for index i);
# mode 0 the public function, the scan with the fast test (_km_find_peaks_run) and the flag 1 when no candidate needed
# the exact predicate; mode 1 every candidate through the exact predicate (force). Id 100, signsc: src/mathx.cyr's
# _km_sign_scaled of n ≤ 8 terms (n, eight significand words, eight exponent words): the sign + 1 in both slots and
# modes.
#must_use
fn run_audit6(fid, mode, a): (i64, i64, i64) {
    if (fid == 100) {
        var sg = _km_sign_scaled(a + 8, a + 72, load64(a));
        return (sg + 1, sg + 1, 1);
    }
    if (fid == 87) {
        var dh: f64 = load64(a);
        var t0: f64 = load64(a + 8);
        var pr: i64[6];
        var n = 0;
        var j = 0;
        while (j < 3) {
            var kk = load64(a + 16 + 16 * j);
            if (kk >= 0) {
                store64(&pr + 16 * n, _km_species6(kk));
                store64(&pr + 16 * n + 8, load64(a + 24 + 16 * j));
                n = n + 1;
            }
            j = j + 1;
        }
        var t, v = adiabatic_flame_temperature(dh, &pr, n, t0);
        var w0 = _km_res5(t, v);
        if (is_ok(t) != 1 && v != KM_ERR_COMPUTATION_ERROR) { return (w0, w0, 0); }
        var cache: i64[6];
        var inf: f64 = _KM_INF;
        var tmax: f64 = inf;
        j = 0;
        while (j < n) {
            var st, ct = _km_segments_lookup(load64(&pr + 16 * j));
            store64(&cache + 16 * j, st);
            store64(&cache + 16 * j + 8, ct);
            var hi: f64 = KM_SHOMATE_T_MAX[st + ct - 1];
            tmax = _km_min(tmax, hi);
            j = j + 1;
        }
        var e, x = _km_flame_bisection(dh, &pr, n, &cache, t0, t0, tmax);
        var xb: f64 = x;
        var w1 = _km_bits(xb);
        if (e != KM_ERR_NONE) { w1 = (1 << 63) | e; }
        return (w0, w1, 0);
    }
    var len = 16;
    if (fid == 88) { len = 3; }
    if (fid == 89) { len = 5; }
    if (fid == 90) { len = 8; }
    var th: f64 = load64(a + 8 * len);
    var ob: i64[16];
    var tp, vp = find_peaks(a, len, th, &ob, 16);
    var mask = 0;
    var q = 0;
    while (q < vp) {
        mask = mask | (1 << load64(&ob + 8 * q));
        q = q + 1;
    }
    var lo: f64 = _KM_INF;
    var hi2: f64 = -lo;
    q = 0;
    while (q < len) {
        var y: f64 = load64(a + 8 * q);
        if (_km_is_finite(y) == 1) {
            lo = _km_min(lo, y);
            hi2 = _km_max(hi2, y);
        }
        q = q + 1;
    }
    var inr = 0;
    if (th >= 0.0 && th <= 1.0) { inr = 1; }
    if (inr != 1 || lo > hi2) { return (mask, mask, 1); }
    var force = 0;
    if (mode == 1) { force = 1; }
    var c, sl = _km_find_peaks_run(a, len, th, lo, hi2, &ob, 16, force);
    var m2 = 0;
    q = 0;
    while (q < c) {
        m2 = m2 | (1 << load64(&ob + 8 * q));
        q = q + 1;
    }
    if (mode == 1) { return (m2, m2, 0); }
    var fl = 0;
    if (sl == 0) { fl = 1; }
    return (mask, m2, fl);
}

# The overflow audit's B6 (ADR 0002 #36): ids 92 … 99, fft_spectrum of n = 1, 2, 4, …, 64 and 256 values, n/2 + 1 words
# per slot: mode 0 the public spectrum (o0), the fast path's (_km_fft_spectrum_fast, o1: Rust's IEEE transform, each
# bin's magnitude Rust's or the unbounded one) and the flag 1 when its transform test passed; mode 1 the unbounded-exponent spectrum forced (_km_fft_spectrum_ux) in both slots for a finite signal,
# the public one otherwise. The work tables are allocated once.
var _fft6_out = 0;
var _fft6_work = 0;
#must_use
fn run_fft6(fid, mode, a, o0, o1): i64 {
    if (_fft6_out == 0) {
        _fft6_out = alloc(8 * 512);
        _fft6_work = alloc(16 * 256);
    }
    var n = arity(fid);
    var m = n / 2 + 1;
    var fin = 1;
    var k = 0;
    while (k < n) {
        var x: f64 = load64(a + 8 * k);
        if (_km_is_finite(x) != 1) { fin = 0; }
        k = k + 1;
    }
    if (mode == 1 && fin == 1) {
        if (_km_fft_spectrum_ux(a, n, o0) != 1) { return 2; }
        k = 0;
        while (k < m) {
            store64(o1 + 8 * k, load64(o0 + 8 * k));
            k = k + 1;
        }
        return 0;
    }
    var t, v = fft_spectrum(a, n, _fft6_out, 2 * n);
    k = 0;
    while (k < m) {
        store64(o0 + 8 * k, load64(_fft6_out + 8 * k));
        k = k + 1;
    }
    if (mode == 1) {
        k = 0;
        while (k < m) {
            store64(o1 + 8 * k, load64(o0 + 8 * k));
            k = k + 1;
        }
        return 0;
    }
    var f = _km_fft_spectrum_fast(a, n, _fft6_work);
    if (f == 2) { f = 1; }
    k = 0;
    while (k < m) {
        store64(o1 + 8 * k, load64(_fft6_work + 8 * k));
        k = k + 1;
    }
    return f;
}

# The overflow audit's B7 (ADR 0002 #39): ids 101 … 106, balance_equation of an m × n integer matrix (bal7_dims), n + 1
# words a slot: the status (n for Ok, 2⁶³ | the Err code) and the n coefficients (0 where Err). Mode 0: the public
# result (o0), the i64 fast path's outcome (_km_balance_i64: its code 0, 1 or 2, then its coefficients where 0; o1) and
# the flag 1 when it did not overflow; mode 1 the multi-word path forced (_km_balance_equation with force) in both.
#must_use
fn run_bal7(fid, mode, a, o0, o1): i64 {
    var m, n = bal7_dims(fid);
    var k = 0;
    while (k <= n) {
        store64(o0 + 8 * k, 0);
        store64(o1 + 8 * k, 0);
        k = k + 1;
    }
    var force = 0;
    if (mode == 1) { force = 1; }
    var t, v = _km_balance_equation(a, m, n, o0 + 8, n, force);
    if (is_ok(t) == 1) {
        store64(o0, v);
    } else {
        store64(o0, (1 << 63) | v);
    }
    if (mode == 1) {
        k = 0;
        while (k <= n) {
            store64(o1 + 8 * k, load64(o0 + 8 * k));
            k = k + 1;
        }
        return 0;
    }
    var f = _km_balance_i64(a, m, n, o1 + 8);
    store64(o1, f);
    if (f != 0) {
        k = 1;
        while (k <= n) {
            store64(o1 + 8 * k, 0);
            k = k + 1;
        }
    }
    if (f == 2) { return 0; }
    return 1;
}

# The overflow audit's B10a (ADR 0002 #40's primitives): ids 107 … 109, _km_cr_exp_ux, _km_cr_expm1_ux and
# _km_cr_exp10_ux of (s, e) (the significand's bits, then the exponent word), two words a slot, the result pair (s, e).
# Mode 0: the function (o0), its fast path (o1: (0, 0) where its test fails) and the flag, 1 when the fast test passed,
# plus 2 when x is a normal double whose binary64 function (_km_cr_exp, _km_cr_expm1, _km_cr_exp10) has a normal result
# other than o0; mode 1: the accurate path alone in both slots.
fn run_b10a(fid, mode, a, o0, o1): i64 {
    var xs: f64 = load64(a);
    var xe = load64(a + 8);
    var s: f64 = 0.0;
    var e = 0;
    if (mode == 1) {
        if (fid == 107) { var s1, e1 = _km_cr_exp_ux_acc(xs, xe); s = s1; e = e1; }
        if (fid == 108) { var s2, e2 = _km_cr_expm1_ux_acc(xs, xe); s = s2; e = e2; }
        if (fid == 109) { var s3, e3 = _km_cr_exp10_ux_acc(xs, xe); s = s3; e = e3; }
        store64(o0, s);
        store64(o0 + 8, e);
        store64(o1, s);
        store64(o1 + 8, e);
        return 0;
    }
    var fs: f64 = 0.0;
    var fe = 0;
    var ok = 0;
    if (fid == 107) {
        var s4, e4 = _km_cr_exp_ux(xs, xe);
        s = s4;
        e = e4;
        var f4, g4, k4 = _km_cr_exp_ux_fast(xs, xe);
        fs = f4;
        fe = g4;
        ok = k4;
    }
    if (fid == 108) {
        var s5, e5 = _km_cr_expm1_ux(xs, xe);
        s = s5;
        e = e5;
        var f5, g5, k5 = _km_cr_expm1_ux_fast(xs, xe);
        fs = f5;
        fe = g5;
        ok = k5;
    }
    if (fid == 109) {
        var s6, e6 = _km_cr_exp10_ux(xs, xe);
        s = s6;
        e = e6;
        var f6, g6, k6 = _km_cr_exp10_ux_fast(xs, xe);
        fs = f6;
        fe = g6;
        ok = k6;
    }
    store64(o0, s);
    store64(o0 + 8, e);
    store64(o1, fs);
    store64(o1 + 8, fe);
    var flag = ok;
    var dx: f64 = xs - xs;
    if (dx == 0.0 && xs != 0.0 && xe >= 0 - 1022 && xe <= 1023) {
        var sc: f64 = _km_cr_pow2(xe);
        var x: f64 = xs * sc;
        var y: f64 = 0.0;
        if (fid == 107) { y = _km_cr_exp(x); }
        if (fid == 108) { y = _km_cr_expm1(x); }
        if (fid == 109) { y = _km_cr_exp10(x); }
        var yb = load64(&y) & 0x7FFFFFFFFFFFFFFF;
        if (yb >= 0x0010000000000000 && yb < 0x7FF0000000000000) {
            var ys, ye = _km_ux_norm(y, 0);
            if (load64(&ys) != load64(&s) || ye != e) { flag = flag | 2; }
        }
    }
    return flag;
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
    if (fid >= 87) {
        var g0, g1, gf = run_audit6(fid, mode, a);
        return (g0, g1, gf);
    }
    if (fid >= 83) {
        var e0, e1, ef = run_audit5(fid, mode, a);
        return (e0, e1, ef);
    }
    if (fid >= 71) {
        var d0, d1, df = run_audit4(fid, mode, a);
        return (d0, d1, df);
    }
    if (fid >= 58) {
        var c0, c1, cf = run_audit3(fid, mode, a);
        return (c0, c1, cf);
    }
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
    if ((fid >= 92 && fid <= 99) || (fid >= 101 && fid <= 109)) {
        # the FFT families: n/2 + 1 words per slot and argument, then the flag bytes; the balancing families n7 + 1;
        # the unbounded exponentials (B10a) 2
        var vw = ar / 2 + 1;
        if (fid >= 101 && fid <= 106) {
            var m7, n7 = bal7_dims(fid);
            vw = n7 + 1;
        }
        if (fid >= 107) { vw = 2; }
        var outv = alloc(16 * vw * n + n + 8);
        var iv = 0;
        while (iv < n) {
            var av = buf + 24 + 8 * ar * iv;
            var fv = 0;
            if (fid >= 107) {
                fv = run_b10a(fid, mode, av, outv + 8 * vw * iv, outv + 8 * vw * n + 8 * vw * iv);
            } elif (fid >= 101) {
                fv = run_bal7(fid, mode, av, outv + 8 * vw * iv, outv + 8 * vw * n + 8 * vw * iv);
            } else {
                fv = run_fft6(fid, mode, av, outv + 8 * vw * iv, outv + 8 * vw * n + 8 * vw * iv);
            }
            store8(outv + 16 * vw * n + iv, fv);
            iv = iv + 1;
        }
        var wv = file_write_all("res.bin", outv, 16 * vw * n + n);
        if (wv != 16 * vw * n + n) { return 4; }
        return 0;
    }
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
    fft6 = [name for _, name, _, _ in FFT6 if only is None or name in only]
    bal7 = [name for _, name, _, _, _ in BAL7 if only is None or name in only]
    b10a = [name for _, name, _ in B10A if only is None or name in only]
    plan = {}
    wmax = 7 * 600_000 if stoich else 0
    if audit or fft6 or bal7 or b10a:
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
    failed |= flame6_bisection_check(work, targets, audit)
    failed |= fft6_checks(work, targets, fft6)
    failed |= bal7_checks(work, targets, bal7)
    failed |= b10a_checks(work, targets, b10a)
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
