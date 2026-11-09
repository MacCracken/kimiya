# scripts/oracle: what the port derives from the Rust oracle

`oracle` runs the frozen Rust kimiya 1.2.0 in [`rust-old/`](../../rust-old/) and writes
two kinds of output:
- **The golden corpus, `tests/golden/<module>.tsv`** (port plan §3 step 5). It has one row
  per call, with every float as its IEEE-754 bit pattern. The Cyrius parity harness
  (`tests/parity.tcyr`, step 7) replays every row against the port and checks it within its
  tolerance class.
- **The generated Cyrius sources, `src/constants.cyr` and (from M1) `src/<module>_data.cyr`**
  (step 6). These are the constants and tables the port reads, so no value is ever typed by
  hand.

The crate, the format and the rules are tanmatra's
(`/home/macro/Repos/tanmatra/scripts/oracle/`), adapted. The one structural difference is
the math library: see [*glibc*](#glibc-not-the-libm-crate).

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write both
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify; exit 1 on any difference
python3 -I scripts/oracle/validate.py .                                            # independent re-check of the corpus
python3 -I scripts/oracle/probe_tables.py .                                        # generated .cyr through the Cyrius compiler
```

## Rules

- **`rust-old/` is never built in place or edited.** This crate depends on it by path, with
  no features (`soorat-compat` and `logging` are deferred, port plan §2.2, and change no value
  the corpus records). `Cargo.lock` was seeded from `rust-old/Cargo.lock` and then pruned by
  cargo, so every shared dependency is the version the oracle was tested with: **hisab 1.4.0**,
  serde 1.0.229, thiserror 2.0.21, … (18 packages, none new but the oracle itself). Build with
  `--locked`.
- **The output is deterministic.** Running it twice gives byte-identical files. Grid inputs use
  seeded xorshift and only exact or correctly rounded operations (+ − × ÷, `floor`,
  `str::parse::<f64>`, bit construction) — never a transcendental std method, timestamps or
  hash-map order. So the *inputs* are the same on any machine; the *outputs* of the
  transcendental rows are glibc's (next section).
- **The corpus is regenerated only from Rust 1.2.0.** A wrong oracle value is fixed in
  Rust: bump to 1.2.x, then regenerate. It is never patched in a `.tsv`.
- **Modules arrive with their milestone.** Each port milestone starts by adding its golden
  emitter (`src/m_<module>.rs`, with the grid from the port plan §4) and its data-module
  emitter (`src/t_<module>.rs`), and registering them in `src/main.rs`.
- **`--check` flags any `.tsv` in `tests/golden/` without an emitter as stale**, and any
  corpus file recorded on another libm platform (next section).
- **`tests/golden/modules.txt`**, also generated, lists every module. `tests/parity.tcyr`
  replays exactly that list, so an emitter added here is replayed at once, as pending until
  its functions get dispatcher entries.

| Module | Milestone | Golden rows | Golden emitter | Generated Cyrius |
|---|---|---|---|---|
| `constants` | M0 | 45 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`): the 14 physical constants, `thermo::WATER_SPECIFIC_HEAT` (M2a) and the ten `solution::` constants (M3a: K_w, Debye–Hückel A / B, Λ°). K_b/K_f and STP are generated with their module's data (M2a); the 10 Madelung and Born rows wait for M4 |
| `electrochemistry` | M3a | 18,574 | `src/m_electrochemistry.rs`: the table and `lookup_half_reaction` (every couple, case / space variants, every prefix, the equations, `""`) (A); the cell potentials over all 49 × 49 couple pairs, ΔG, Faraday's laws over every decade of charge and mass, the electron counts 0 … 6, 127, 128, 254, 255, the edges (B); `nernst_potential` / `_25c` (C) over the edges, T and Q decades, every couple's E° and n, a seeded sweep and **near-cancellation rows** (E° = (RT/nF)·ln Q to a few ulp, E ≈ 0 of both signs) | `src/electrochemistry_data.cyr` (`src/t_electrochemistry.rs`): the 49 `STANDARD_POTENTIALS` rows (a **UTF-8** equation blob, the couple blob, electrons `u8`, E° f64), `KM_STANDARD_POTENTIALS_COUNT`, `HalfReaction`'s docs and the literal's two CRC Handbook notes |
| `element` | M1 | 975 | `src/m_element.rs`: `lookup_by_number` and `lookup_properties` over every u8 (0..=255), `lookup_by_symbol` over every symbol, its lower- and upper-case forms, every name, and `""`, `"Xx"`, `"FE"`, `"Fe "`, `" Fe"`, a 300-byte string | `src/element_data.cyr` (`src/t_element.rs`): `enum KmElementCategory`, the 118 `ELEMENTS` rows (atomic number, symbol and name blobs, mass, electronegativity flag + value, category) and the 32 `ELEMENT_PROPERTIES` rows (oxidation states flattened over an `i8` table) |
| `math` | M0 step 8 (M2a: +7,500) | 26,147 | `src/m_math.rs`; M2a added dense cos / acos / ln / powf(·, 2.5) sweeps over the Peng–Robinson and Redlich–Kwong argument ranges | — |
| `molecule` | M1 | 2,472 | `src/m_molecule.rs`: the four presets, the cases the 1.2 tests pin, 600 seeded compositions (1–8 atoms, counts 0..=200, repeats, invalid Z 0/119/255, u32::MAX counts); `new`, `molecular_weight` (B), `formula`, `total_atoms` on each | — |
| `reaction` | M2a | 1,188 | `src/m_reaction.rs`: `gibbs_free_energy` (B) and `equilibrium_constant` (C) over each argument's edges, ΔG from target exponents across exp's underflow/overflow ends, the Rust tests' values, seeded sweeps; `is_spontaneous` (A) | — |
| `solution` | M3a | 28,566 | `src/m_solution.rs`: molarity, molality, dilution, ionic strength (`v:n` + `f:c i:z` pairs), Kohlrausch, conductivity (B); pH / pOH over every decade, h from pH across 10^x's ends, the activity models over every charge −4 … 4, ±50 000, i32::MIN / MAX (C); `henderson_hasselbalch` (C) with **near-cancellation rows** (pKa = −log₁₀(base/acid) to a few ulp, pH ≈ 0 of both signs); `weak_acid_ph` / `weak_base_ph` (C) over the Ka (Kb) × C grid of every decade 10⁻⁴⁰ … 10⁴⁰ × 10⁻⁴⁰ … 10²⁰, seeded pairs, the overflow cases, **each row's note naming the Newton's exit** (`[root: relative step, 4 updates]`) | `src/constants.cyr` gains the ten solution constants (solution has no data file) |
| `thermo` | M2a | 1,107 | `src/m_thermo.rs`: all class B — `heat_transfer`, `hess_law` and `enthalpy_from_formation` (slices as `v:n` + values; −0.0 sums, ∞ − ∞, cancellation, seeded slices), `specific_heat_from_calorimetry` (the 1e-12 threshold ± 1 ulp) | `src/constants.cyr` gains `KM_WATER_SPECIFIC_HEAT` (thermo has no data file) |
| `phase` | M2a | 3,276 | `src/m_phase.rs`: Raoult, Henry, colligative, osmotic (B); Clausius–Clapeyron (C); Antoine mmHg / Pa (C) for every set at and one ulp past each range end, 0.5/1/10 °C outside, a 5 °C ladder, seeded points, the pole ± 1 ulp, specials, custom sets (non-finite fields, inverted range, pole inside the range, overflow), and the unit set's class-A rows that pin `PA_PER_MMHG`; `lookup_antoine`, `phase_rule` (all pairs of 13 u32 probes), `ANTOINE_DATA` / `ANTOINE_DATA[]` (A) | `src/phase_data.cyr` (`src/t_phase.rs`): `struct KmAntoineCoeffs`, the 8 sets (name blob + A, B, C, t_min_c, t_max_c columns), `_KM_PA_PER_MMHG` (through the API, see below), `KM_KB_*` / `KM_KF_*` |
| `thermochem` | M2b | 9,487 | `src/m_thermochem.rs`: the tables and lookups (A); `Shomate::cp` / `enthalpy_antiderivative` per segment incl. range ends, extrapolation, T = 0 / ±∞ / NaN, `enthalpy_change_cp` within and across segments (B); `reaction_enthalpy` / `_gibbs_energy` (B, **J/mol: the oracle's kJ/mol × 1000.0, computed by the generator**, ADR 0002 #1, said in the row note) and `reaction_entropy` over seeded species lists (`v:n` + `s:formula f:coeff` pairs); `vant_hoff_k` (C); `adiabatic_flame_temperature` (B) over fuels, equivalence ratios, oxidisers, T₀, 10⁻²⁰⁰…10³⁰⁰ scalings and every Err path, **each Ok row's note naming the solver that answered** (`[newton, k steps]`, `[bisection after newton: …, k halvings]`) | `src/thermochem_data.cyr` (`src/t_thermochem.rs`): `struct KmShomate`, the 76 `THERMOCHEM_DATA` rows (formula and name blobs, ΔH_f°, ΔG_f°, S°), the 34 `SHOMATE_DATA` segments (formula blob, a … t_max), `KM_THERMOCHEM_COUNT`, `KM_SHOMATE_COUNT`, and the literals' `//` comments |
| `gas` | M2a | 19,717 | `src/m_gas.rs`: ideal gas, van der Waals, Dalton, density, `peng_robinson_pressure` (B); `peng_robinson_z`, `_z_roots` (`v:n` + roots), `fugacity_coefficient_pr`, `fugacity_pr` over each preset × 26 temperatures × 29 pressures (one- and three-root regions, the critical point ± 10⁻⁶, 90–200 MPa), seeded states, validation edges and overflow, and `redlich_kwong_pressure` (C); the `PR_*` presets and `joule_thomson_ideal` (A) | `src/gas_data.cyr` (`src/t_gas.rs`): `struct KmPengRobinsonParams`, `KM_STP_*`, the 5 `KM_PR_*` presets as top-level struct values |
| `kinetics` | M3b | 30,146 | `src/m_kinetics.rs`: `half_life_first_order` (every decade of k; k = 1 is LN_2 exactly), the second- and zero-order laws, `zero_order_half_life`, steady state, pre-equilibrium, Michaelis–Menten and `lineweaver_burk` (output `f:slope f:intercept`, the tuple's fields) (B), each argument over the edges and every Err path; `arrhenius_rate`, `first_order_concentration`, `collision_theory_rate`, `eyring_rate`, `eyring_rate_from_activation` (C) over T from the smallest subnormal to DBL_MAX, wide signed energies, the exponent on targets across exp's range, and **on purpose the rows where a subnormal exp result is scaled into the normal range** (A, [A]₀, Z·p up to DBL_MAX; Eyring at low T; one Eyring factor subnormal, the other huge), exp overflow (A·∞, 0·∞ → NaN); `nth_order_half_life` (C) on a fine order grid around 1 (the 10⁻¹² tolerance and the ulps either side), orders to ±300, [A]₀ and k over many decades, **[A]₀^(n−1) subnormal with a finite half-life**, and exact powers at the subnormal grid's half-way points | — (kinetics' constants are `src/constants.cyr`'s) |

## glibc, not the `libm` crate

tanmatra's oracle called the pinned pure-Rust `libm` crate, so its corpus was the same on
every machine. kimiya 1.2.0 calls the **std f64 methods**, and on Linux those reach several
different implementations. Read off the release oracle's imports (`nm -D`), rustc 1.99.1,
x86_64-unknown-linux-gnu:

| Implementation | Functions | Bits depend on |
|---|---|---|
| **glibc libm** (dynamic) | exp, exp2, expm1, log, log10, pow, sin, cos, acos, atan2 | the glibc version, and on x86_64 the CPU: glibc picks an FMA variant (ifunc) when the CPU has FMA and AVX2 |
| **Rust compiler-builtins** (static: the `libm` crate's code) | cbrt (CORE-MATH, correctly rounded), floor, trunc, round, `__powidf2` (powi) | nothing: exact or correctly rounded |
| **inline code** | sqrt (`sqrtsd`), abs, min, max, signum, clamp, total_cmp, is_finite | nothing |

So the header of every `.tsv` records the platform, e.g.
`# oracle: kimiya 1.2.0 (rust-old/), hisab 1.4.0; libm: glibc 2.44 x86_64 fma.` (version from
`gnu_get_libc_version()`, `src/glibc.rs`). **`--check` refuses a corpus recorded on another
platform** with a `PLATFORM` line naming both, before comparing anything: its transcendental
rows can only be checked where they were made. Regenerating on another platform is a
deliberate act that changes the reference; record why.

Two more consequences, both recorded rather than assumed:
- **LLVM rewrites `2.0_f64.powf(y)` to `exp2(y)`** at opt-level ≥ 1 (the oracle binary imports
  `exp2`). kimiya writes that form for half-lives, so `math.tsv` records it as written, under
  `f64::powf_base2`; `10.0_f64.powf(y)` stays `pow` (`f64::powf_base10`). The corpus comes from
  the release-built oracle; Rust's own debug-built tests call `pow` for both.
- **`max`/`min` of +0 and −0 return the second operand** (`maxsd`/`minsd`). Rust leaves that
  case unspecified; the rows record what rustc 1.99.1 does, and the port's helpers must match.

`validate.py` checks every glibc row against a correctly rounded reference, so the reference
is known good on every row (the histogram is below).

## Format 1

### Header

Each file starts with five `#` lines:
- the format and module;
- the oracle provenance (kimiya and hisab versions) and the libm platform;
- the regenerate command;
- the column list;
- the row count.

Lines starting with `#` are comments.

### Rows

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **`fn`** is the Rust path of what was called:
  - `module::function` or `Type::method` (the receiver is the first input) from M1 on;
  - `const` for a constant, whose input is its path, `s:element::AVOGADRO`;
  - `f64::<method>` for a std math reference value, `self` first (`f64::atan2` is
    `y.atan2(x)`, inputs `y x`); `f64::powf_base10` / `f64::powf_base2` are
    `10.0_f64.powf(y)` / `2.0_f64.powf(y)` with input `y`.
- **`class`** is `A`, `B`, `C` or `D` (see *Comparison*).
- **`inputs` and `outputs`** are space-separated tokens, or `-` when empty.
  - A struct is flattened into its fields in declaration order; a tuple in order.
- **`note`** is the same call in readable form, e.g. `f64::exp(-1.5) -> 0.22313016014842982`.
  It is for review only: the harness ignores it.

### Tokens

| Token | Meaning |
|---|---|
| `f:<16 hex>` | f64, by IEEE-754 bit pattern (lowercase hex) |
| `i:<decimal>` | integer (u8 / u32 / i64 / usize); booleans are `0` / `1`; `total_cmp` is -1 / 0 / 1 |
| `e:<Enum>.<Variant>` | enum variant by its Rust name. In Cyrius, enum values follow the Rust declaration order and are append-only |
| `s:<text>` | UTF-8 string. Bytes outside printable ASCII (0x21–0x7e), and `%` itself, are `%XX` (uppercase hex) |
| `n:` | `Option::None`. `Some(x)` is just `x` |
| `E:<Variant>` | `Err(KimiyaError::<Variant>)`. Only the variant is compared, not the payload |
| `v:<len>` | a sequence of `len` elements, whose tokens follow in order. The harness knows each function's element shape |

A function returning a `Vec` too long for one row gets a row with its length (`v:<len>`) plus
one row per element under `fn[]`, whose input is the index `i:k`.

### Comparison

**Non-float tokens** (`i`, `e`, `s`, `n`, `E`, `v` lengths) must be equal in every class.

**Floats** are compared by class:

| Class | A float output passes when |
|---|---|
| **A** exact | the bits are identical. Any NaN matches any NaN, because NaN payloads and signs differ between implementations |
| **B** arithmetic | the ulp distance is ≤ 2, measured on the ordered-integer line, where +0 and −0 are the same point. **Expected 0**: a non-zero distance means the operation order differs from Rust |
| **C** transcendental | the ulp distance is ≤ K, a bound per function (`_par_c_bound` in `tests/parity.tcyr`). K is derived and recorded in `docs/development/port-parity.md` (M0 step 8) before the function's first comparison. A function with no K yet is *measured*: its histogram is reported but not gated |
| **D** iterative | \|got − want\| ≤ 1e-12 · \|want\| |

In classes B, C and D:
- NaN matches only NaN, and ±∞ match only the same infinity;
- a zero `want` needs a zero `got`, of either sign (class D; in B and C, +0 and −0 are one point
  on the ulp line).

A tolerance is never widened to make a row pass (port plan §5). A divergence accepted by an
ADR in `docs/adr/` is listed by its exact function and input field in `_par_accepted`. It is
reported as `accepted (ADR n)` on every run.

**Shadow lines.** The harness can also replay a row through a second Cyrius candidate
(`fn ~ candidate`, e.g. `f64::exp ~ _f64_exp_polyfill`, `f64::max ~ f64_max`). A shadow line
reports match / differ counts and the histogram and is never gated: it is evidence for the
step-8 and M1 choices, not a test.

## Generated Cyrius sources

Values come from the oracle's public API. The Rust source is read only for the `///` doc
comments, which carry the citations (SI 2019, CODATA 2018, AME2020, …) into the Cyrius file,
and `t_constants::complete` proves the list is every `pub const NAME: f64` in
`rust-old/src` — no more, no fewer — before anything is written.

### `PA_PER_MMHG`: a private constant through the API

Rust's `phase::PA_PER_MMHG` is private. `t_phase.rs` reads it as `antoine_pressure_pa` of the unit set
(A = B = 0, C = 1, range [0, 0]) at T = 0: `antoine_pressure_mmhg` of that set is 10⁰ = 1.0 exactly (checked),
so the Pa value is the factor itself. Its golden row is class A, and the source literal is a cross-check.

### The flame solver's path, and ADR 0002 #1's × 1000

Two kinds of golden row carry a bracketed suffix in their note (`Module::row_noted`; the harness ignores notes):
- `thermochem::reaction_enthalpy` / `reaction_gibbs_energy`: `[J/mol: the oracle's kJ/mol result × 1000.0 (ADR 0002
  #1)]`. The Cyrius functions return J/mol; the generator multiplies Rust's kJ/mol value by 1000.0, the port's one
  multiply after the same subtraction, so the rows compare directly (class B, 0 ulp).
- `thermochem::adiabatic_flame_temperature`: the solver that produced an `Ok` value. `m_thermochem::flame_replica`
  rebuilds the solver from the public API (`SHOMATE_DATA`, `Shomate::cp`, `Shomate::enthalpy_antiderivative`, Rust's
  private `Segments` logic, hisab 1.4.0's `newton_raphson` and `bisection` loops), and generation fails unless it
  gives the oracle's bits, or the oracle's error, on every row. `validate.py` recomputes the path independently and
  checks the note. The generator prints the tally (1,682 Newton, 154 bisection, 101 Err).

### The weak acid / base Newton's path (M3a)

`solution::weak_acid_ph` / `weak_base_ph` rows that reach the iteration carry its exit in their note:
`[root: f(x) <= 0, k updates]`, `[root: relative step, k updates]`, `[root: non-finite step, k updates]`,
`[root: 50 iterations]`, with `, checked_root rejected` if the accepted x was not positive and finite.
`m_solution::root_replica` transcribes Rust's private `charge_balance_root` and `checked_root`, and generation fails
unless it gives the oracle's pH bits, or its error, on every row. `validate.py` recomputes the path independently and
checks the note. The generator prints the tally (11,740 `f(x) <= 0`, 6,356 relative step, 74 non-finite step, no
exhaustion or rejection; at most 7 updates).

### Emission rules

Probed on Cyrius 6.7.18 (by tanmatra; re-proven here by `probe_tables.py`), and enforced by
the generator:
- **f64s are decimal literals written positionally.** They use the shortest digits that
  round-trip, and each is checked to parse back to the oracle's bits; the lexer has no
  exponent notation. Decimal literals have been correctly rounded since 6.6.10.
- **Hex is never put in a `const`.** There it is an integer, and `CONST * x` does integer
  arithmetic.
- **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** A NaN is a generator error: Cyrius refuses a
  compile-time NaN, whose bits depend on the target.
- **Tables (from M1) are baked column arrays**, `var KM_<X>_<FIELD>: f64|u8|u16|i8|bool[N] = { .. };`,
  one element per line with a label comment (`cyr::table`, tanmatra's). `i8` arrays sign-extend on
  read (probed on x86_64 ELF/PE and aarch64), so oxidation states keep Rust's `i8`.
- **String columns are a NUL-separated `u8` blob plus a `u16` offset column** (`cyr::string_column`):
  an array initializer can't hold a string literal. `&_KM_X_BYTES + _KM_X_OFF[i]` is row i's C string;
  one blob line per row, its bytes then the text as the label (a row past 120 columns, M2b's longest names, is
  split into lines of 16 bytes with the label on the first). ASCII and NUL-free, or a generator error — except a
  column written with `cyr::string_column_utf8` (M3a: the half-reaction equations), whose blob holds each string's UTF-8
  bytes exactly (NUL-free), so `strlen` counts bytes and a C-string compare is Rust's `&str ==`.
- **Comments inside a table literal** (M2a's methanol source, M2b's sections and Shomate sources / caveats, M3a's two
  trailing "CRC Handbook" notes, keyed by row and field) travel into the table's comment block, keyed by row and
  word-wrapped; `t_thermochem` and `t_electrochemistry` cross-check each literal's row keys against the API's rows, in
  order.
- **`Option<f64>` cells are a `bool` presence column plus the value** (0.0 for `None`), since a
  compile-time f64 can't be NaN; the accessor returns the value-form `Option`.
- **Per-cell `//` notes in the Rust literal** (the IUPAC conventional and CIAAW 2024 atomic weights)
  travel as that cell's label comment. Only comments are read from the source.
- **The output is `cyrius fmt`- and `cyrius lint`-clean as generated.**

### Names

A bundle is concatenated into its consumers, so every generated global is prefixed `KM_`:

| Oracle / golden token | Generated Cyrius |
|---|---|
| `const` row `s:element::AVOGADRO` | `const KM_AVOGADRO` (the generator refuses two modules' constants with one name) |
| enum `ElementCategory` (M1) | `enum KmElementCategory`, plus `const KM_ELEMENT_CATEGORY_COUNT`; variant `KM_` + snake(enum) + `_` + snake(variant), as tanmatra (`e:ElementCategory.NobleGas` is `KM_ELEMENT_CATEGORY_NOBLE_GAS`). `t_element` proves the variant list against the Rust declaration |
| `ELEMENTS[i].field` (M1) | `KM_ELEMENT_<FIELD>[i]`, row i = Z − 1; symbol and name `_KM_ELEMENT_{SYMBOL,NAME}_{BYTES,OFF}`; `KM_ELEMENT_COUNT` |
| `ELEMENT_PROPERTIES[k].field` (M1) | `KM_ELEMENT_PROPERTIES_<FIELD>[k]`; oxidation states `..._OXIDATION_STATE_START` / `_COUNT` over `..._OXIDATION_STATES`; `KM_ELEMENT_PROPERTIES_COUNT` |
| `THERMOCHEM_DATA[k].field` (M2b) | `KM_THERMOCHEM_<FIELD>[k]`; formula and name `_KM_THERMOCHEM_{FORMULA,NAME}_{BYTES,OFF}`; `KM_THERMOCHEM_COUNT` |
| `SHOMATE_DATA[k]` (M2b) | `KM_SHOMATE_<FIELD>[k]` (A … T_MAX); formula `_KM_SHOMATE_FORMULA_{BYTES,OFF}`; `KM_SHOMATE_COUNT`; `Shomate` is `struct KmShomate` |
| `STANDARD_POTENTIALS[k].field` (M3a) | `KM_STANDARD_POTENTIALS_ELECTRONS[k]`, `KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL[k]`; equation and couple `_KM_STANDARD_POTENTIALS_{EQUATION,COUPLE}_{BYTES,OFF}`; `KM_STANDARD_POTENTIALS_COUNT`; `HalfReaction` has no struct: its fields are `half_reaction_<field>(h)` |

Each constant's comment block ends with `# Rust: kimiya::<module>::<NAME>`, which
`validate.py` and `probe_tables.py` read to tie it to its golden row.

### `probe_tables.py`

`probe_tables.py` compiles a program that prints every generated constant and every
`element_data.cyr` cell (1,055 in M1: each `ELEMENTS` field, symbol and name read through the blob
offsets, each property and oxidation state, and each category checked against the enum member the
golden variant names), and from M2a every `phase_data.cyr` / `gas_data.cyr` value (65 cells: each Antoine
name and coefficient, `KM_ANTOINE_COUNT`, `_KM_PA_PER_MMHG`, each preset field read through its struct; and
their 10 constants), and from M2b every `thermochem_data.cyr` cell (654: each THERMOCHEM_DATA formula, name and
f64 field, each SHOMATE_DATA formula and the seven columns, both counts), and from M3a every `electrochemistry_data.cyr`
cell (197: each STANDARD_POTENTIALS equation — its UTF-8 bytes — couple, electron count and E°, the count) and the ten
solution constants. It reads the arrays directly, not through `src/element.cyr`, so it tests the
data while the harness tests the accessors. It builds in `build/oracle-probe/` with the project's
manifest, then compares each value bit for bit with the golden rows. It fails on any mismatch, any
compiler warning, a generated constant with no golden row, and a table it doesn't probe. A mass, an
oxidation state or a symbol byte tampered by one digit fails (checked in M1; in M2a an Antoine B, a name
byte, the mmHg factor and a preset's ω).

## `validate.py`

`validate.py` is an independent re-check of the corpus that doesn't use the Rust code. It
needs only Python 3, with no packages (run it with `-I`).
- **Format:** the five-line header (with the true row count), five fields per row, a known
  class, well-formed tokens (UTF-8 inside `s:`), and `modules.txt` listing exactly the `.tsv`
  files.
- **Constants:** `constants.tsv` matches every `pub const NAME: f64` literal in
  `rust-old/src/*.rs` (Python's `float()` is correctly rounded, like rustc), and every
  `const KM_X = <decimal>;` in `src/constants.cyr` is positional and parses to its row's bits.
- **math, exact functions,** recomputed with Rust's semantics: sqrt (`math.sqrt` is correctly
  rounded), abs, floor, trunc, round (half away from zero, exactly), min/max (a NaN operand is
  ignored; for ±0 only the value is checked and the observed operand rule is printed), signum,
  clamp, total_cmp, is_finite, and powi (compiler-builtins' loop). **cbrt is proven correctly
  rounded on every row with exact rational arithmetic**, which is why it is class A.
- **math, class C** (glibc's bits): the IEEE / C99 Annex F special cases, then every finite row
  against a **correctly rounded reference** from the `decimal` module at 60+ digits (exp, ln,
  log10, pow directly; sin, cos by reduction mod 2π with a 450-digit π; acos, atan2 through a
  halving atan series). Each must be within glibc's documented bound (1 ulp; log10 2 ulp).
  Measured on glibc 2.44 x86_64 fma (ulp from correctly rounded):

  | Function | rows | 0 ulp | 1 ulp |
  |---|---|---|---|
  | exp, atan2, powf_base10 | 2,696 | all | 0 |
  | ln | 3,857 | 3,854 | 3 |
  | cos | 2,910 | 2,906 | 4 |
  | acos | 2,620 | 2,619 | 1 |
  | exp_m1 | 810 | 776 | 34 |
  | log10 | 1,296 | 1,248 | 48 |
  | powf | 1,750 | 1,747 | 3 (2 of them powf(x, 2.5)) |
  | powf_base2 (exp2) | 369 | 368 | 1 |
  | sin | 710 | 708 | 2 |

  (Before M2a's sweeps, ln, cos and acos had no 1-ulp row: glibc 2.44 is not correctly rounded on about
  0.1 % of arguments.)

- **element** (M1): every row against the `ELEMENTS` and `ELEMENT_PROPERTIES` literals parsed from
  `rust-old/src/element.rs` (Python's `float()` is correctly rounded, like rustc), with
  `lookup_by_symbol` recomputed as an exact, case-sensitive match; every u8 covered; and the
  generated `KmElementCategory` in the Rust declaration order.
- **molecule** (M1): every row recomputed from those masses with Rust's semantics — the molecular
  weight summed in atom order (`total += mass * count`, bit-exact: Python floats are IEEE doubles),
  the formula (input order, count 1 without a digit, zero counts skipped after the element check),
  saturating `total_atoms`, `InvalidElement` for Z outside 1–118, and the presets parsed from
  `rust-old/src/molecule.rs`.

- **reaction, thermo, phase, gas** (M2a): every row recomputed from Rust 1.2.0's code transcribed in Python
  in its operation order, with the validation guards in Rust's order and IEEE division (Python raises on a zero
  divisor). Class B rows must match bit for bit. The class-C rows call Python's `math`, which is glibc's
  exp, log, pow, acos and cos (and a correctly rounded cbrt, as Rust's compiler-builtins cbrt is), so they too
  match bit for bit when the machine is the corpus's recorded libm platform; elsewhere only their Ok/Err shape
  is checked. The tables and presets (`ANTOINE_DATA`, `PR_*`, `PA_PER_MMHG`) are parsed from the Rust source
  literals, and the generated `src/phase_data.cyr` / `src/gas_data.cyr` literals are checked against them.
- **thermochem** (M2b): every row recomputed from Rust 1.2.0's code transcribed in Python — `sum_property` from
  +0.0 with `+=`, the J/mol × 1000.0 of ADR 0002 #1, `Segments` (lookup, cp, enthalpy_change with Rust's max/min),
  the Shomate polynomial and antiderivative, `vant_hoff_k` (glibc's exp via Python), and the flame temperature with
  hisab 1.4.0's `newton_raphson` and `bisection` loops; the solver path each flame row's note claims is checked
  against the recomputation, at least 100 rows must be bisection rows, and each of Rust's typed flame constants
  (10⁻³, 10⁻⁶, the 40 J/(mol·K) guess) must change some golden row's bits when perturbed. `THERMOCHEM_DATA` and
  `SHOMATE_DATA` are parsed from the Rust literals, and the generated `src/thermochem_data.cyr` columns, string
  blobs and `struct KmShomate` are checked against them.
- **solution, electrochemistry** (M3a): every row recomputed from Rust 1.2.0's code transcribed in Python — the guards in
  Rust's order, Rust's `.sum()` from −0.0 for `ionic_strength`, the weak acid / base Newton (`charge_balance_root`,
  `checked_root`) with its exit path checked against each row's note, and `CHARGE_BALANCE_REL_TOL` (10⁻¹²) required to
  change some golden row's bits when halved and when doubled (the 50-iteration cap is reported as unobservable: no row
  enters more than 7 iterations); `henderson_hasselbalch`'s log10 as glibc's e_log10.c composition over log, checked to
  be glibc's log10 bit for bit on every golden argument; `STANDARD_POTENTIALS` parsed from the `electrochemistry.rs`
  literal, and the generated `src/electrochemistry_data.cyr` (the UTF-8 equation blob, the couple blob, electrons, E°,
  the two CRC notes) checked against it.
- **kinetics** (M3b): every row recomputed from Rust 1.2.0's code transcribed in Python (glibc's exp, expm1 and pow via
  Python's `math`), the guards in Rust's order, `lineweaver_burk`'s tuple as two floats; LN_2 (± 1 ulp) and
  `FIRST_ORDER_TOLERANCE` (× ½, × 2) must each change some golden row (± 1 ulp of the tolerance is reported as
  unobservable: near 1 ± 10⁻¹², |n − 1| moves in steps of 2⁻⁵² or 2⁻⁵³).
- **`--sensitivity`** (M2a; M2b adds `vant_hoff_k`; M3a the solution and Nernst functions; M3b the kinetics functions):
  `python3 -I scripts/oracle/validate.py . --sensitivity` re-evaluates every class-C row with one primitive result
  nudged one ulp up and down, one call at a time, and prints each function's S, c and K = 2·S·c + 2 — the port-parity
  rule's measurement (`docs/development/port-parity.md` §5–§7). Following ADR 0005, only calls whose result can differ
  between port and oracle are nudged: exp, 10^x and (M3a) `_km_log10` results, and the acos / cos / ln / x^2.5 results
  where glibc is not correctly rounded (the port's are; the 60-digit references decide). For the functions the M3a
  remedy changed (`henderson_hasselbalch`, `nernst_potential`, `nernst_potential_25c`) it also prints S "before": every
  result of the stdlib-like primitive the remedy replaced (`_km_log10`, `f64_ln`) nudged — about 10¹⁸ for each. Every
  kinetics exp, exp_m1 and pow result is nudged (stdlib-like counterparts). For each function whose K passes 64 (M3b: the
  six class-C kinetics functions, left measured until ADR 0006) it also lists the rows whose own S puts it there, how
  many of them nudge a subnormal primitive result, the worst three, and S over the remaining rows and over the rows with
  no subnormal primitive result (port-parity.md §8).

Every function in every `.tsv` must have a rule in `validate.py`: a new one fails until it gets one. Extend it
as modules arrive. A check that recomputes a row independently is worth more than one that re-reads the
oracle.
