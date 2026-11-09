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
| `constants` | M0 | 45 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`): the 14 physical constants. The other 31 rows (K_b/K_f, Madelung, Born exponents, Λ°, Debye–Hückel, K_w, STP, c_water) are generated with their module's data at M2–M4 |
| `element` | M1 | 975 | `src/m_element.rs`: `lookup_by_number` and `lookup_properties` over every u8 (0..=255), `lookup_by_symbol` over every symbol, its lower- and upper-case forms, every name, and `""`, `"Xx"`, `"FE"`, `"Fe "`, `" Fe"`, a 300-byte string | `src/element_data.cyr` (`src/t_element.rs`): `enum KmElementCategory`, the 118 `ELEMENTS` rows (atomic number, symbol and name blobs, mass, electronegativity flag + value, category) and the 32 `ELEMENT_PROPERTIES` rows (oxidation states flattened over an `i8` table) |
| `math` | M0 step 8 | 18,647 | `src/m_math.rs` | — |
| `molecule` | M1 | 2,472 | `src/m_molecule.rs`: the four presets, the cases the 1.2 tests pin, 600 seeded compositions (1–8 atoms, counts 0..=200, repeats, invalid Z 0/119/255, u32::MAX counts); `new`, `molecular_weight` (B), `formula`, `total_atoms` on each | — |

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
  one blob line per row, its bytes then the text as the label. ASCII and NUL-free, or a generator error.
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

Each constant's comment block ends with `# Rust: kimiya::<module>::<NAME>`, which
`validate.py` and `probe_tables.py` read to tie it to its golden row.

### `probe_tables.py`

`probe_tables.py` compiles a program that prints every generated constant and every
`element_data.cyr` cell (1,055 in M1: each `ELEMENTS` field, symbol and name read through the blob
offsets, each property and oxidation state, and each category checked against the enum member the
golden variant names). It reads the arrays directly, not through `src/element.cyr`, so it tests the
data while the harness tests the accessors. It builds in `build/oracle-probe/` with the project's
manifest, then compares each value bit for bit with the golden rows. It fails on any mismatch, any
compiler warning, a generated constant with no golden row, and a table it doesn't probe. A mass, an
oxidation state or a symbol byte tampered by one digit fails (checked in M1).

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
  | exp, ln, cos, acos, atan2, powf_base10 | 5,583 | all | 0 |
  | exp_m1 | 810 | 776 | 34 |
  | log10 | 1,296 | 1,248 | 48 |
  | powf | 750 | 749 | 1 |
  | powf_base2 (exp2) | 369 | 368 | 1 |
  | sin | 710 | 708 | 2 |

- **element** (M1): every row against the `ELEMENTS` and `ELEMENT_PROPERTIES` literals parsed from
  `rust-old/src/element.rs` (Python's `float()` is correctly rounded, like rustc), with
  `lookup_by_symbol` recomputed as an exact, case-sensitive match; every u8 covered; and the
  generated `KmElementCategory` in the Rust declaration order.
- **molecule** (M1): every row recomputed from those masses with Rust's semantics — the molecular
  weight summed in atom order (`total += mass * count`, bit-exact: Python floats are IEEE doubles),
  the formula (input order, count 1 without a digit, zero counts skipped after the element check),
  saturating `total_atoms`, `InvalidElement` for Z outside 1–118, and the presets parsed from
  `rust-old/src/molecule.rs`.

Every function in `math.tsv`, `element.tsv` and `molecule.tsv` must have a rule in `validate.py`:
a new one fails until it gets one. Extend it as modules arrive. A check that recomputes a row
independently is worth more than one that re-reads the oracle.
