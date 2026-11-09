"""Compile every generated constant and table cell through the Cyrius compiler and check each
value's bits against tests/golden. This proves the generated sources reproduce the oracle's
bits as the compiler reads them (decimal-literal conversion, i8 sign extension, u16 offsets
and the NUL-separated string blobs included), not just as Python reads them.

Usage, from the repo root:
    python3 -I scripts/oracle/probe_tables.py .

It writes build/oracle-probe/ (gitignored), builds there with the project's cyrius.cyml
(which supplies the stdlib), and exits 1 on any mismatch, any compiler warning, or a
generated value with no golden row. The probe reads the baked arrays directly, never through
src/element.cyr, so it tests the data and the harness tests the accessors. Adapted from
tanmatra's probe_tables.py; extend `checks` as the data modules arrive.

- src/constants.cyr: every `const KM_X` against its constants.tsv row.
- src/element_data.cyr (M1): for every row of ELEMENTS, the atomic number, symbol, name,
  atomic mass, electronegativity (presence flag and value) and category, against the golden
  `element::lookup_by_number` row of that Z; for every row of ELEMENT_PROPERTIES, the atomic
  number, IE1, EA, covalent radius and the flattened oxidation states, against the golden
  `element::lookup_properties` row. The category enum is checked by name: each golden
  `e:ElementCategory.<Variant>` must be the value of KM_ELEMENT_CATEGORY_<SNAKE> as compiled.
"""
import re, subprocess, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
work = root / "build/oracle-probe"
work.mkdir(parents=True, exist_ok=True)

def rows(name):
    for line in (root / f"tests/golden/{name}.tsv").read_text().splitlines():
        if not line.startswith("#"):
            fn, cls, ins, outs, _ = line.split("\t")
            yield fn, cls, ins.split(" "), outs.split(" ")

def unpct(t):
    return re.sub(r"%([0-9A-F]{2})", lambda m: chr(int(m.group(1), 16)), t)

def snake(camel):
    return re.sub(r"(?<=[a-z0-9])([A-Z])", r"_\1", camel).upper()

checks = []  # (cyrius statement printing one value, kind 'f' | 'i' | 's', expected text, label)
fails = 0
def fail(msg):
    global fails
    fails += 1
    print("FAIL", msg)

# ── src/constants.cyr ───────────────────────────────────────────────────────────────
golden = {}  # "element::AVOGADRO" -> "f:<bits>"
for fn, cls, ins, outs in rows("constants"):
    golden[ins[0][2:]] = outs[0]
cyr = (root / "src/constants.cyr").read_text()
nconst = 0
for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = ", cyr, re.M):
    path, name = m.groups()
    if path not in golden:
        fail(f"{name}: kimiya::{path} has no golden row")
        continue
    checks.append((f"x = {name}; fmt_hex(x);", "f", golden[path][2:], f"const {path}"))
    nconst += 1
declared = len(re.findall(r"^const ", cyr, re.M))
if declared != nconst:
    fail(f"src/constants.cyr declares {declared} consts, {nconst} carry a `# Rust:` path")

# ── src/element_data.cyr ────────────────────────────────────────────────────────────
data = (root / "src/element_data.cyr").read_text()
by_z = {}
for fn, cls, ins, outs in rows("element"):
    if fn == "element::lookup_by_number" and outs != ["n:"]:
        by_z[int(ins[0][2:])] = outs
props = []
for fn, cls, ins, outs in rows("element"):
    if fn == "element::lookup_properties" and outs != ["n:"]:
        props.append(outs)
count = int(re.search(r"^const KM_ELEMENT_COUNT = (\d+);", data, re.M).group(1))
pcount = int(re.search(r"^const KM_ELEMENT_PROPERTIES_COUNT = (\d+);", data, re.M).group(1))
if sorted(by_z) != list(range(1, count + 1)):
    fail(f"golden lookup_by_number finds Z {sorted(by_z)[:3]}..; KM_ELEMENT_COUNT is {count}")
if len(props) != pcount:
    fail(f"golden lookup_properties finds {len(props)} rows; KM_ELEMENT_PROPERTIES_COUNT is {pcount}")
for k in range(count):
    want = by_z.get(k + 1)
    if want is None:
        continue
    z, sym, name, mass, en, cat = want
    lab = f"ELEMENTS[{k}]"
    checks.append((f"fmt_int(KM_ELEMENT_ATOMIC_NUMBER[{k}]);", "i", z[2:], f"{lab}.atomic_number"))
    checks.append((f"p = &_KM_ELEMENT_SYMBOL_BYTES + _KM_ELEMENT_SYMBOL_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(sym[2:]), f"{lab}.symbol"))
    checks.append((f"p = &_KM_ELEMENT_NAME_BYTES + _KM_ELEMENT_NAME_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(name[2:]), f"{lab}.name"))
    checks.append((f"x = KM_ELEMENT_ATOMIC_MASS[{k}]; fmt_hex(x);", "f", mass[2:], f"{lab}.atomic_mass"))
    has = "0" if en == "n:" else "1"
    checks.append((f"if (KM_ELEMENT_HAS_ELECTRONEGATIVITY[{k}]) {{ fmt_int(1); }} else {{ fmt_int(0); }}",
                   "i", has, f"{lab}.electronegativity is_some"))
    checks.append((f"x = KM_ELEMENT_ELECTRONEGATIVITY[{k}]; fmt_hex(x);", "f",
                   "0" * 16 if en == "n:" else en[2:], f"{lab}.electronegativity"))
    variant = cat[2:].split(".", 1)[1]
    checks.append((f"fmt_int(KM_ELEMENT_CATEGORY[{k}] - KM_ELEMENT_CATEGORY_{snake(variant)});", "i", "0",
                   f"{lab}.category is KM_ELEMENT_CATEGORY_{snake(variant)}"))
for k, want in enumerate(props):
    z, ie, ea, rad, vlen = want[:5]
    states = want[5:]
    lab = f"ELEMENT_PROPERTIES[{k}]"
    checks.append((f"fmt_int(KM_ELEMENT_PROPERTIES_ATOMIC_NUMBER[{k}]);", "i", z[2:], f"{lab}.atomic_number"))
    checks.append((f"x = KM_ELEMENT_PROPERTIES_IONIZATION_ENERGY_1[{k}]; fmt_hex(x);", "f", ie[2:],
                   f"{lab}.ionization_energy_1"))
    checks.append((f"x = KM_ELEMENT_PROPERTIES_ELECTRON_AFFINITY[{k}]; fmt_hex(x);", "f", ea[2:],
                   f"{lab}.electron_affinity"))
    checks.append((f"x = KM_ELEMENT_PROPERTIES_COVALENT_RADIUS_PM[{k}]; fmt_hex(x);", "f", rad[2:],
                   f"{lab}.covalent_radius_pm"))
    checks.append((f"fmt_int(KM_ELEMENT_PROPERTIES_OXIDATION_STATE_COUNT[{k}]);", "i", vlen[2:],
                   f"{lab}.oxidation_states.len"))
    for j, st in enumerate(states):
        checks.append((f"fmt_int(KM_ELEMENT_PROPERTIES_OXIDATION_STATES[KM_ELEMENT_PROPERTIES_OXIDATION_STATE_START[{k}] + {j}]);",
                       "i", st[2:], f"{lab}.oxidation_states[{j}]"))
# Every declared table is covered above; a new column must be added here first.
covered = {"KM_ELEMENT_ATOMIC_NUMBER", "_KM_ELEMENT_SYMBOL_OFF", "_KM_ELEMENT_SYMBOL_BYTES",
           "_KM_ELEMENT_NAME_OFF", "_KM_ELEMENT_NAME_BYTES", "KM_ELEMENT_ATOMIC_MASS",
           "KM_ELEMENT_HAS_ELECTRONEGATIVITY", "KM_ELEMENT_ELECTRONEGATIVITY", "KM_ELEMENT_CATEGORY",
           "KM_ELEMENT_PROPERTIES_ATOMIC_NUMBER", "KM_ELEMENT_PROPERTIES_IONIZATION_ENERGY_1",
           "KM_ELEMENT_PROPERTIES_ELECTRON_AFFINITY", "KM_ELEMENT_PROPERTIES_COVALENT_RADIUS_PM",
           "KM_ELEMENT_PROPERTIES_OXIDATION_STATE_START", "KM_ELEMENT_PROPERTIES_OXIDATION_STATE_COUNT",
           "KM_ELEMENT_PROPERTIES_OXIDATION_STATES"}
tables = set(re.findall(r"^var (_?KM_[A-Z0-9_]+):", data, re.M))
if tables != covered:
    fail(f"src/element_data.cyr tables not probed: {sorted(tables - covered)}; probed but gone: {sorted(covered - tables)}")
nelem = len(checks) - nconst

# ── the probe program ───────────────────────────────────────────────────────────────
src = ['include "src/constants.cyr"', 'include "src/element_data.cyr"', "",
       "fn main() {", "    var x: f64 = 0.0;", "    var p = 0;"]
for k, (stmt, kind, want, label) in enumerate(checks):
    src.append(f"    {stmt} println(\" {k}\");")
src += ["    return 0;", "}", "", "var r = main();", "sys_exit_group(r);", ""]
probe = work / "tables_probe.cyr"
probe.write_text("\n".join(src))

b = subprocess.run(["cyrius", "build", str(probe), str(work / "tables_probe")],
                   capture_output=True, text=True, cwd=root)
diag = [l for l in (b.stdout + b.stderr).splitlines() if re.search(r"warn|error|refus", l, re.I)]
if b.returncode != 0 or diag:
    print("BUILD:", b.returncode, *diag, sep="\n")
    sys.exit(2)
r = subprocess.run([str(work / "tables_probe")], capture_output=True, text=True)
got = {}
for line in r.stdout.splitlines():
    val, k = line.rsplit(" ", 1)
    got[int(k)] = val

for k, (stmt, kind, want, label) in enumerate(checks):
    g = got.get(k)
    if kind == "f":
        ok = g is not None and int(g, 16) == int(want, 16)
    else:
        ok = g == want
    if not ok:
        fail(f"{label}: {stmt} -> {g!r}, want {want!r}")
pending = len(golden) - nconst
print(f"checked {nconst} constants and {nelem} element_data cells through the Cyrius compiler "
      f"({pending} golden constants are module data, generated at their milestone); failures: {fails}")
sys.exit(1 if fails else 0)
