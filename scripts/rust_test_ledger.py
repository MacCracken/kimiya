"""The Rust-test ledger: every one of the oracle's tests is ported, folded into golden rows, or recorded with a reason.

Port plan M7 (the mela pattern; adapted from tanmatra's scripts/rust_test_ledger.py). Usage, from the repo root:
    python3 -I -B scripts/rust_test_ledger.py .           # write docs/development/rust-parity.md
    python3 -I -B scripts/rust_test_ledger.py . --check   # exit 1 if the ledger is stale or a test is unaccounted

It lists each `#[test]` fn in rust-old/src/**.rs and rust-old/tests/*.rs, and each doc test (a Rust code block in a
`///` or `//!` comment, and the README's, which src/lib.rs includes as ReadmeDoctests), and finds its port: a
`fn _t_<name>(` in the Cyrius test files FILES gives for that Rust file, called from the file's main, whose
`test_group` is the Rust test's name (or the name and a parenthesised remark). Each Rust test is then:
- ported: it has that port; ADAPTED says why a port asserts something other than Rust's assertion (an ADR 0002 API
  change, or the flat API); every other port asserts what Rust asserts in the port's form (port plan §6.4);
- folded: it has no port, and FOLDED names the golden rows that check it instead;
- n/a: REASONS says why it does not apply to the Cyrius library;
- deferred: its whole file waits on a §2.2 module (WHOLE_FILE).
A Rust test with no status is an error, and so is an ADAPTED, FOLDED or REASONS entry that names no Rust test, or a
ported one (FOLDED, REASONS). The page also carries the whole-surface table (`scripts/surface_parity.py --markdown`).
"""
import re
import subprocess
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
check = "--check" in sys.argv[2:]
out_path = root / "docs/development/rust-parity.md"

# Rust file -> the Cyrius test files its tests are ported into (None: its tests are deferred with the file).
MODULES = ["error", "validate", "element", "molecule", "reaction", "thermo", "thermochem", "phase", "gas", "solution",
           "electrochemistry", "kinetics", "stoichiometry", "organic", "inorganic", "potential", "nuclear",
           "spectroscopy", "reaction_dynamics", "network", "fitting", "sensitivity"]
FILES = {f"src/{m}.rs": [f"tests/{m}.tcyr"] for m in MODULES}
FILES.update({
    "src/bridge.rs": None,
    "src/integration/soorat.rs": None,
    "src/logging.rs": None,
    # The integration tests: 16 in tests/integration.tcyr (M7), 8 ported earlier beside their module's tests.
    "tests/integration.rs": ["tests/integration.tcyr", "tests/organic.tcyr", "tests/spectroscopy.tcyr",
                             "tests/reaction_dynamics.tcyr", "tests/fitting.tcyr"],
})
NO_TESTS = ["src/lib.rs", "src/integration/mod.rs"]   # files with no #[test] (lib.rs carries the README doc test)

# The doc tests: (Rust file, the item it documents) -> (Cyrius file, test name).
DOC_TESTS = {
    ("src/lib.rs", "ReadmeDoctests"): ("tests/integration.tcyr", "doc_readme_quick_start"),
    ("src/thermochem.rs", "adiabatic_flame_temperature"): ("tests/integration.tcyr", "doc_adiabatic_flame_temperature"),
}
DOC_WHAT = {
    ("src/lib.rs", "ReadmeDoctests"): "the README's Quick Start (rust-old/README.md, included as `ReadmeDoctests`)",
    ("src/thermochem.rs", "adiabatic_flame_temperature"): "`adiabatic_flame_temperature`'s example",
}

ADR1 = ("ADR 0002 #1: `reaction_enthalpy` / `reaction_gibbs_energy` return J/mol, so Rust's kJ/mol values and "
        "tolerances are asserted × 1000")
ADR6 = "ADR 0002 #6: the renamed `{}`"
ADR7 = "ADR 0002 #7: "
FLAT = "the flat API: "
# (Rust file, test) -> how its port's assertion differs from Rust's.
ADAPTED = {
    ("src/thermochem.rs", "combustion_of_methane_enthalpy"): ADR1,
    ("src/thermochem.rs", "combustion_of_methane_gibbs"): ADR1,
    ("src/thermochem.rs", "haber_process_enthalpy"): ADR1,
    ("src/thermochem.rs", "added_standard_state_elements"): ADR1,
    ("src/phase.rs", "antoine_pole_and_non_finite_give_nan"): (
        "ADR 0002 #4: each case is Err where Rust asserts NaN (\"just above the pole\" lies outside water's fitted "
        "range)"),
    ("src/nuclear.rs", "semi_empirical_invalid_is_nan"): (
        "ADR 0002 #3: `semi_empirical_binding_energy` is a Result, so Err(InvalidInput) where Rust asserts NaN"),
    ("src/organic.rs", "vsepr_clf3_t_shaped"): ADR6.format("KM_GEOMETRY_T_SHAPED") + " (Rust's `TShapedThreeLonePairs`)",
    ("src/organic.rs", "alcohol_polar_and_hbonds"): ADR6.format("functional_group_is_hydrogen_bond_donor"),
    ("src/organic.rs", "ether_polar_not_hbond_donor"): ADR6.format("functional_group_is_hydrogen_bond_donor"),
    ("src/organic.rs", "acceptor_only_groups_are_not_donors"): ADR6.format("functional_group_is_hydrogen_bond_donor"),
    ("src/organic.rs", "carboxylic_acid_hbonds"): ADR6.format("functional_group_is_hydrogen_bond_donor"),
    ("src/stoichiometry.rs", "balance_non_integer_entry_is_error"): ADR7 + (
        "the integer matrix: 1.5, NaN and ∞ have no counterpart (the corpus counts those rows); Err(InvalidInput) for "
        "the integers Rust's 2⁵³ bound refuses, i64::MAX and i64::MIN standing in for 1e300"),
    ("src/stoichiometry.rs", "balance_near_integer_entry_accepted"): ADR7 + (
        "the integer matrix: 2 + 10⁻¹² and −1 − 10⁻¹² are the integers 2 and −1 Rust rounds them to"),
    ("src/stoichiometry.rs", "oxidation_state_superoxide_is_error"): ADR7 + "Ok with KO₂'s exact O −1/2 where Rust "
                                                                           "asserts Err",
    ("src/stoichiometry.rs", "oxidation_state_propane_non_integer_is_error"): ADR7 + "Ok with C₃H₈'s exact C −8/3 where "
                                                                                    "Rust asserts Err",
    ("src/stoichiometry.rs", "oxidation_state_magnetite_non_integer_is_error"): ADR7 + "Ok with Fe₃O₄'s exact Fe +8/3 "
                                                                                      "where Rust asserts Err",
    ("src/stoichiometry.rs", "oxidation_state_single_element"): ADR7 + "O₂⁻'s state is the exact −1/2",
    ("src/stoichiometry.rs", "oxidation_state_huge_count_does_not_panic"): ADR7 + (
        "(O, u32::MAX)×2 at i32::MIN is the exact −2³⁰/(2³² − 1)"),
    ("src/spectroscopy.rs", "hydrogen_abundances_iupac_2009"): (
        "ADR 0002 #23: the exact pairs are the cited table's (¹H 1.007825031898, ²H 2.014101777844), and the Σ = 1 loop "
        "gains sulfur (Rust's ³⁶S-less sulfur sums to 0.9999)"),
    ("src/sensitivity.rs", "invalid_input_gives_non_finite_gradients"): (
        "ADR 0002 #2: the port's Err (with its code) where Rust asserts a non-finite entry"),
    ("src/network.rs", "rate_ragged_matrix_is_error"): FLAT + (
        "a flat matrix has one row length, so the group checks the nearest case, a rate vector of another length than "
        "every row"),
    ("src/reaction_dynamics.rs", "stiff_resized_jacobian_is_error_not_panic"): FLAT + (
        "the flat n × n Jacobian cannot be resized, so the group checks the nearest case, a Jacobian callback that "
        "writes nothing (J = 0, no panic either)"),
    ("tests/integration.rs", "methane_combustion_cross_module"): ADR1 + " (−800 kJ/mol is −800,000 J/mol)",
    ("tests/integration.rs", "thermochem_enthalpy_entropy_gibbs_triangle"): ADR1 + (
        " (ΔH − T·ΔS in J/mol, Rust's 1 kJ/mol tolerance as 1,000 J/mol)"),
}
DOC_ADAPTED = {
    ("src/thermochem.rs", "adiabatic_flame_temperature"): (
        "ADR 0002 #1: `reaction_enthalpy` is J/mol already, so the example's `dh_kj * 1000.0` is the J/mol value as "
        "it is"),
}

# (Rust file, test) -> the golden rows that check it in place of a port. None: every in-scope test is ported.
FOLDED = {}

SERDE = "n/a: a serde round-trip; serde is not carried over and no Cyrius consumer exchanges the type (port plan §2.4)"
# (Rust file, test) -> why it does not apply.
REASONS = {
    ("src/molecule.rs", "serde_roundtrip"): SERDE,
}

BRIDGE = ("deferred: `bridge.rs` (the ushma, tanmatra and bijli bridges) waits on a Cyrius consumer; tanmatra deferred "
          "its kimiya bridges likewise, so the two sets port together (port plan §2.2)")
SOORAT = ("deferred: `integration/soorat.rs` (feature `soorat-compat`) waits on soorat's port; its output is a serde JSON "
          "shape (port plan §2.2)")
LOGGING = ("deferred: `logging.rs` (feature `logging`) wraps tracing-subscriber and kimiya emits no log events; the "
           "Cyrius analogue would be sakshi, if a consumer asks (port plan §2.2)")
WHOLE_FILE = {"src/bridge.rs": BRIDGE, "src/integration/soorat.rs": SOORAT, "src/logging.rs": LOGGING}


def rust_tests(rel):
    lines = (root / "rust-old" / rel).read_text().split("\n")
    names = []
    for i, line in enumerate(lines):
        if line.strip().startswith("#[test]"):
            for nxt in lines[i + 1 : i + 8]:
                m = re.match(r"\s*(?:pub )?fn (\w+)", nxt)
                if m:
                    names.append(m.group(1))
                    break
    return names


def doc_tests():
    """(Rust file, item) for every Rust code block in a doc comment, and the README's (lib.rs's ReadmeDoctests)."""
    found = []
    for path in sorted((root / "rust-old/src").rglob("*.rs")):
        rel = str(path.relative_to(root / "rust-old"))
        lines = path.read_text().split("\n")
        k = 0
        while k < len(lines):
            m = re.match(r"\s*//[/!]\s?```(\w*)", lines[k])
            if m:
                end = k + 1
                while end < len(lines) and not re.match(r"\s*//[/!]\s?```", lines[end]):
                    end += 1
                if m.group(1) in ("", "rust", "no_run", "should_panic"):
                    item = None
                    for nxt in lines[end + 1 :]:
                        t = re.match(r"\s*pub(?:\([\w:]+\))? (?:const )?(?:fn|struct|enum|const|static|mod) (\w+)", nxt)
                        if t:
                            item = t.group(1)
                            break
                    found.append((rel, item or "the module"))
                k = end + 1
                continue
            m = re.match(r'\s*#\[doc = include_str!\("(.*)"\)\]', lines[k])
            if m:
                inc = (path.parent / m.group(1)).resolve()
                blocks = len(re.findall(r"^```rust\s*$", inc.read_text(), re.M))
                ty = None
                for nxt in lines[k + 1 : k + 4]:
                    t = re.match(r"\s*struct (\w+)", nxt)
                    if t:
                        ty = t.group(1)
                        break
                found.extend([(rel, ty)] * blocks)
            k += 1
    return found


def cyrius_ports(rel):
    """name -> (file, test_group) for every `fn _t_<name>(` that its file's main calls."""
    ports = {}
    for f in rel or []:
        text = (root / f).read_text()
        called = set(re.findall(r"^\s+_t_(\w+)\(\);", text, re.M))
        for m in re.finditer(r"^fn _t_(\w+)\(\) \{\n\s+test_group\(\"([^\"]*)\"\);", text, re.M):
            if m.group(1) in called:
                ports[m.group(1)] = (f, m.group(2))
    return ports


def cyrius_only(f, rust_names):
    text = (root / f).read_text()
    called = set(re.findall(r"^\s+_t_(\w+)\(\);", text, re.M))
    return sorted(n for n in re.findall(r"^fn _t_(\w+)\(", text, re.M) if n in called and n not in rust_names)


errors = []
rows = []        # (rust file, cyrius files, tests, counts)
adapted = []     # (rust file, test, cyrius file, reason)
notes = []       # (rust file, test, status text) for the tests without a port
elsewhere = []   # (rust file, test, cyrius file) for integration tests ported beside a module's tests
all_rust = set()
for rel in FILES:
    all_rust |= set(rust_tests(rel))
for rel in NO_TESTS:
    if rust_tests(rel):
        errors.append(f"{rel} has #[test]s; give it a FILES entry")
known = set(FILES) | set(NO_TESTS)
for path in sorted(list((root / "rust-old/src").rglob("*.rs")) + list((root / "rust-old/tests").glob("*.rs"))):
    rel = str(path.relative_to(root / "rust-old"))
    if rel not in known:
        errors.append(f"{rel} is in neither FILES nor NO_TESTS")
for rel, cy in FILES.items():
    tests = rust_tests(rel)
    ports = cyrius_ports(cy)
    counts = {"ported": 0, "adapted": 0, "folded": 0, "n/a": 0, "deferred": 0}
    for t in tests:
        if t in ports:
            f, group = ports[t]
            if group != t and not group.startswith(t + " ("):
                errors.append(f"{rel}: {t}'s port in {f} is test_group(\"{group}\"), not the Rust name")
            counts["ported"] += 1
            if (rel, t) in ADAPTED:
                counts["adapted"] += 1
                adapted.append((rel, t, f, ADAPTED[(rel, t)]))
            for table, what in ((FOLDED, "FOLDED"), (REASONS, "REASONS")):
                if (rel, t) in table:
                    errors.append(f"{rel}: {t} is ported but has a {what} entry")
            if cy and f != cy[0]:
                elsewhere.append((rel, t, f))
            continue
        if (rel, t) in ADAPTED:
            errors.append(f"{rel}: {t} has an ADAPTED entry but no port")
        if (rel, t) in FOLDED:
            counts["folded"] += 1
            notes.append((rel, t, "folded: " + FOLDED[(rel, t)]))
            continue
        reason = REASONS.get((rel, t)) or WHOLE_FILE.get(rel)
        if reason is None:
            errors.append(f"{rel}: {t} is neither ported, folded nor recorded")
            continue
        counts[reason.split(":")[0]] += 1
        if rel not in WHOLE_FILE:
            notes.append((rel, t, reason))
    for table, what in ((ADAPTED, "ADAPTED"), (FOLDED, "FOLDED"), (REASONS, "REASONS")):
        for (r, t) in table:
            if r == rel and t not in tests:
                errors.append(f"{rel}: {what} names {t}, which is not a test")
    rows.append((rel, cy, len(tests), counts))
found_docs = doc_tests()
if sorted(found_docs) != sorted(DOC_TESTS):
    errors.append(f"the doc tests are {sorted(found_docs)}; DOC_TESTS lists {sorted(DOC_TESTS)}")
doc_rows = []
for key, (cy, name) in DOC_TESTS.items():
    ports = cyrius_ports([cy])
    if name not in ports:
        errors.append(f"doc test {key}: {cy} has no _t_{name} called from its main")
    elif ports[name][1] != name:
        errors.append(f"doc test {key}: _t_{name}'s test_group is \"{ports[name][1]}\"")
    doc_rows.append((key, cy, name))

# Cyrius-only groups per Cyrius test file (the ported files only).
cy_files = []
for rel, cy in FILES.items():
    for f in cy or []:
        if f not in cy_files:
            cy_files.append(f)
only = {f: cyrius_only(f, all_rust | {n for _, n in DOC_TESTS.values()}) for f in cy_files}

n_unit = sum(r[2] for r in rows if r[0].startswith("src/"))
n_int = sum(r[2] for r in rows if r[0].startswith("tests/"))
n_doc = len(DOC_TESTS)
total = n_unit + n_int + n_doc
sums = {k: sum(r[3][k] for r in rows) for k in ("ported", "adapted", "folded", "n/a", "deferred")}
sums["ported"] += n_doc
sums["adapted"] += len(DOC_ADAPTED)
in_scope_unit = sum(r[2] for r in rows if r[0].startswith("src/") and r[0] not in WHOLE_FILE)
m7_int = [t for t, (f, _) in cyrius_ports(["tests/integration.tcyr"]).items() if t in all_rust]

surface = subprocess.run([sys.executable, "-I", "-B", str(root / "scripts/surface_parity.py"), str(root), "--markdown"],
                         capture_output=True, text=True)
if surface.returncode != 0:
    errors.append("scripts/surface_parity.py --markdown failed: " + surface.stderr.strip())

md = []
md.append("# Rust-test ledger and surface parity")
md.append("")
md.append("> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan M7, the mela")
md.append("> pattern). Every test of the frozen Rust oracle (`rust-old/`, kimiya 1.2.0: `#[test]` fns, the integration")
md.append("> tests and the doc tests) is ported to a Cyrius test, folded into golden rows, or recorded here with its")
md.append("> reason; every public item of `rust-old/src` is ported, renamed by ADR 0002, n/a or deferred. The golden")
md.append("> corpus (`tests/golden/`, `tests/parity.tcyr`) checks every ported function on its grid besides; this page is")
md.append("> about the oracle's own assertions and its API.")
md.append("")
md.append("## Summary")
md.append("")
md.append(f"**{total} Rust tests** ({n_unit} unit, {n_int} integration, {n_doc} doc): **{sums['ported']} ported** "
          f"({sums['adapted']} of them adapted to an ADR 0002 change or the flat API), {sums['folded']} folded into "
          f"golden rows, {sums['n/a']} not applicable, {sums['deferred']} deferred with their modules (port plan §2.2).")
md.append(f"Of the unit tests, the {in_scope_unit} in scope (port plan §2.1) are {in_scope_unit - sums['n/a']} ported "
          f"and {sums['n/a']} n/a. M7 ported the last {len(m7_int) + n_doc}: {len(m7_int)} integration tests and the "
          f"{n_doc} doc tests, in `tests/integration.tcyr`.")
md.append("")
md.append("A port's `test_group` is the Rust test's name, and the script checks that the Cyrius file's `main` calls it. "
          "A port asserts what Rust asserts, in the port's form (port plan §6.4: `Result::unwrap` is `is_ok` + the "
          "payload, a lookup's `Option<&T>` a handle, a Vec the caller's buffer, a slice a pointer and a count, a "
          "value Rust writes with an exponent written positionally); the adapted ones are listed below.")
md.append("")
md.append("| Rust file | Cyrius tests | Rust tests | Ported | Adapted | Folded | N/A | Deferred | Cyrius-only groups |")
md.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
for rel, cy, n, c in rows:
    cyt = ", ".join(f"`{f}`" for f in cy) if cy else "—"
    extra = len(only[cy[0]]) if cy else "—"
    md.append(f"| `{rel}` | {cyt} | {n} | {c['ported']} | {c['adapted']} | {c['folded']} | {c['n/a']} | "
              f"{c['deferred']} | {extra} |")
md.append(f"| doc tests (`src/lib.rs`, `src/thermochem.rs`) | `tests/integration.tcyr` | {n_doc} | {n_doc} | "
          f"{len(DOC_ADAPTED)} | 0 | 0 | 0 | — |")
md.append(f"| **Total** | | **{total}** | **{sums['ported']}** | **{sums['adapted']}** | **{sums['folded']}** | "
          f"**{sums['n/a']}** | **{sums['deferred']}** | **{sum(len(v) for v in only.values())}** |")
md.append("")
md.append("The Cyrius-only groups (`cyrius_*`, and element's `out_of_range_handles`, validate's `guard_edges`, error's "
          "`codes`) test what no Rust test reaches: the port's conventions (handles, caller buffers, words for Rust's "
          "integer types), the ADR 0002 changes and the correctly rounded primitives. `tests/mathx.tcyr` and "
          "`tests/numerics.tcyr` hold only such tests (the primitives, and hisab 1.4.0's transcribed routines against "
          "their golden rows); `tests/parity.tcyr` replays the corpus and `tests/kimiya.tcyr` smoke-tests the bundle.")
md.append("")
md.append("## Integration and doc tests")
md.append("")
md.append(f"`rust-old/tests/integration.rs` has {n_int} tests. {len(m7_int)} are in `tests/integration.tcyr` (M7), "
          "which includes the bundle `dist/kimiya.cyr` as a consumer does; the other "
          f"{len(elsewhere)} were ported beside their module's tests:")
md.append("")
for rel, t, f in elsewhere:
    md.append(f"- `{t}`: `{f}`")
md.append("")
md.append("The doc tests:")
md.append("")
for key, cy, name in doc_rows:
    note = f" Adapted: {DOC_ADAPTED[key]}." if key in DOC_ADAPTED else ""
    md.append(f"- `rust-old/{key[0]}`, {DOC_WHAT[key]}: `{cy}` `_t_{name}`.{note}")
md.append("")
md.append("## Adapted ports")
md.append("")
md.append("| Rust file | Test | Cyrius file | Adaptation |")
md.append("|---|---|---|---|")
for rel, t, f, reason in adapted:
    md.append(f"| `{rel}` | `{t}` | `{f}` | {reason} |")
md.append("")
md.append("## Tests without a port")
md.append("")
md.append("| Rust file | Test | Status |")
md.append("|---|---|---|")
for rel, t, reason in notes:
    md.append(f"| `{rel}` | `{t}` | {reason} |")
for rel, reason in WHOLE_FILE.items():
    left = rust_tests(rel)
    md.append(f"| `{rel}` | all {len(left)} | {reason} |")
md.append("")
md.append("No test is folded into golden rows: every in-scope Rust test has its own port, and the corpus checks the "
          "functions on their grids besides.")
md.append("")
md.append("## Surface parity")
md.append("")
md.append("`python3 -I -B scripts/surface_parity.py .` checks every public item of `rust-old/src` (before `#[cfg(test)]`; "
          "every `pub` line is one of the items counted) against `src/`: a ported module's item by the Cyrius name the "
          "porting rules give (CLAUDE.md, port plan §6.4), or in another form with its reason, or renamed by ADR 0002 "
          "(the new name must exist); the §2.2 modules' items are deferred, and `lib.rs`'s module declarations and "
          "prelude re-exports n/a (each re-export must name a public item of its module).")
md.append("")
md.append(surface.stdout.rstrip())
md.append("")
text = "\n".join(md) + "\n"

for e in errors:
    print("LEDGER", e)
if check:
    old = out_path.read_text() if out_path.exists() else ""
    if old != text:
        print(f"STALE    {out_path.relative_to(root)} (regenerate: python3 -I -B scripts/rust_test_ledger.py .)")
        errors.append("stale")
else:
    out_path.write_text(text)
print(f"rust-test ledger: {total} tests — {sums['ported']} ported ({sums['adapted']} adapted), {sums['folded']} "
      f"folded, {sums['n/a']} n/a, {sums['deferred']} deferred; {len(errors)} problem(s)")
sys.exit(1 if errors else 0)
