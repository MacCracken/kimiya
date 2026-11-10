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

- src/constants.cyr: every `const KM_X` against its constants.tsv row (thermo::WATER_SPECIFIC_HEAT
  included, M2a).
- src/phase_data.cyr (M2a): the K_b / K_f constants against constants.tsv; KM_ANTOINE_COUNT, and
  for every ANTOINE_DATA row the name (read through the blob offsets) and the five coefficient
  columns, against the golden `phase::ANTOINE_DATA[]` rows; _KM_PA_PER_MMHG against the golden
  class-A `phase::antoine_pressure_pa` row of the unit set (A = B = 0, C = 1 at T = 0).
- src/gas_data.cyr (M2a): the STP constants against constants.tsv; each KM_PR_* preset's tc, pc
  and omega, read through the struct, against its golden `gas::PR_*` row.
- src/thermochem_data.cyr (M2b): KM_THERMOCHEM_COUNT and KM_SHOMATE_COUNT; for every THERMOCHEM_DATA row
  the formula and name (through the blob offsets) and ΔH_f°, ΔG_f°, S°, against the golden
  `thermochem::THERMOCHEM_DATA[]` rows; for every SHOMATE_DATA segment the formula and the seven
  coefficient / range columns, against the golden `thermochem::SHOMATE_DATA[]` rows.
- src/electrochemistry_data.cyr (M3a): KM_STANDARD_POTENTIALS_COUNT; for every STANDARD_POTENTIALS row the
  equation (its UTF-8 bytes through the blob offsets) and couple, the electron count and E°, against the golden
  `electrochemistry::STANDARD_POTENTIALS[]` rows. src/constants.cyr's checks cover the ten solution constants.
- src/nuclear_data.cyr (M4c): KM_ISOTOPES_COUNT; for every ISOTOPES row the label (through the blob offsets), Z (u8),
  A (u16), the half-life and the decay mode (its UTF-8 bytes through the blob offsets), against the golden
  `nuclear::ISOTOPES[]` rows. The four nuclear constants are src/constants.cyr's.
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

def unpct_utf8(t):
    """A golden `s:` token's text, its %XX escapes decoded as the UTF-8 bytes they are."""
    return re.sub(rb"%([0-9A-F]{2})", lambda m: bytes([int(m.group(1), 16)]), t.encode()).decode("utf-8")

def snake(camel):
    """scripts/oracle/src/cyr.rs's snake (tanmatra's rule): `_` before a capital after a lowercase letter or digit, or
    before a capital that starts a capitalised word after a run of capitals (`CCDouble` → `CC_DOUBLE`, M4b)."""
    out = []
    for k, ch in enumerate(camel):
        if ch.isupper() and k > 0:
            prev, nxt = camel[k - 1], camel[k + 1] if k + 1 < len(camel) else ""
            if prev.islower() or prev.isdigit() or (prev.isupper() and nxt.islower()):
                out.append("_")
        out.append(ch.upper())
    return "".join(out)

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

# ── src/phase_data.cyr, src/gas_data.cyr (M2a) ──────────────────────────────────────
for rel in ("src/phase_data.cyr", "src/gas_data.cyr"):
    text = (root / rel).read_text()
    for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = ", text, re.M):
        path, name = m.groups()
        if path not in golden:
            fail(f"{rel}: {name}: kimiya::{path} has no golden row")
            continue
        checks.append((f"x = {name}; fmt_hex(x);", "f", golden[path][2:], f"const {path}"))
        nconst += 1
nm2a_before = len(checks)
pdata = (root / "src/phase_data.cyr").read_text()
ant = {}
ant_len = None
unit_pa = None
for fn, cls, ins, outs in rows("phase"):
    if fn == "phase::ANTOINE_DATA[]":
        ant[int(ins[0][2:])] = outs
    elif fn == "phase::ANTOINE_DATA":
        ant_len = int(outs[0][2:])
    elif fn == "phase::antoine_pressure_pa" and cls == "A":
        unit_pa = outs[0]
acount = int(re.search(r"^const KM_ANTOINE_COUNT = (\d+);", pdata, re.M).group(1))
if ant_len != acount or sorted(ant) != list(range(acount)):
    fail(f"golden ANTOINE_DATA has {ant_len} rows ({sorted(ant)}); KM_ANTOINE_COUNT is {acount}")
checks.append(("fmt_int(KM_ANTOINE_COUNT);", "i", str(ant_len), "KM_ANTOINE_COUNT"))
for k in range(acount):
    name, a, b, c, lo, hi = ant[k]
    lab = f"ANTOINE_DATA[{k}]"
    checks.append((f"p = &_KM_ANTOINE_NAME_BYTES + _KM_ANTOINE_NAME_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(name[2:]), f"{lab}.name"))
    for col, want in (("A", a), ("B", b), ("C", c), ("T_MIN_C", lo), ("T_MAX_C", hi)):
        checks.append((f"x = KM_ANTOINE_{col}[{k}]; fmt_hex(x);", "f", want[2:], f"{lab}.{col.lower()}"))
if unit_pa is None:
    fail("phase.tsv has no class-A antoine_pressure_pa row (the unit set) to pin _KM_PA_PER_MMHG")
else:
    checks.append(("x = _KM_PA_PER_MMHG; fmt_hex(x);", "f", unit_pa[2:], "_KM_PA_PER_MMHG"))
covered_p = {"_KM_ANTOINE_NAME_OFF", "_KM_ANTOINE_NAME_BYTES", "KM_ANTOINE_A", "KM_ANTOINE_B",
             "KM_ANTOINE_C", "KM_ANTOINE_T_MIN_C", "KM_ANTOINE_T_MAX_C"}
tables_p = set(re.findall(r"^var (_?KM_[A-Z0-9_]+):", pdata, re.M))
if tables_p != covered_p:
    fail(f"src/phase_data.cyr tables not probed: {sorted(tables_p - covered_p)}; probed but gone: {sorted(covered_p - tables_p)}")
gdata = (root / "src/gas_data.cyr").read_text()
presets = {}
for fn, cls, ins, outs in rows("gas"):
    if re.fullmatch(r"gas::PR_[A-Z0-9_]+", fn):
        presets[fn.split("::")[1]] = outs
declared = re.findall(r"^var KM_(PR_[A-Z0-9_]+) = KmPengRobinsonParams \{", gdata, re.M)
if sorted(declared) != sorted(presets):
    fail(f"src/gas_data.cyr presets {sorted(declared)}; golden gas::PR_* rows {sorted(presets)}")
for name in declared:
    for field, want in zip(("tc", "pc", "omega"), presets.get(name, [])):
        checks.append((f"x = KM_{name}.{field}; fmt_hex(x);", "f", want[2:], f"KM_{name}.{field}"))
other = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", gdata, re.M)) - {f"KM_{n}" for n in declared}
if other:
    fail(f"src/gas_data.cyr globals not probed: {sorted(other)}")
nm2a = len(checks) - nm2a_before

# ── src/thermochem_data.cyr (M2b) ───────────────────────────────────────────────────
nm2b_before = len(checks)
tdata = (root / "src/thermochem_data.cyr").read_text()
tc_rows, sh_rows, tc_len, sh_len = {}, {}, None, None
for fn, cls, ins, outs in rows("thermochem"):
    if fn == "thermochem::THERMOCHEM_DATA[]":
        tc_rows[int(ins[0][2:])] = outs
    elif fn == "thermochem::SHOMATE_DATA[]":
        sh_rows[int(ins[0][2:])] = outs
    elif fn == "thermochem::THERMOCHEM_DATA":
        tc_len = int(outs[0][2:])
    elif fn == "thermochem::SHOMATE_DATA":
        sh_len = int(outs[0][2:])
tcount = int(re.search(r"^const KM_THERMOCHEM_COUNT = (\d+);", tdata, re.M).group(1))
scount = int(re.search(r"^const KM_SHOMATE_COUNT = (\d+);", tdata, re.M).group(1))
if tc_len != tcount or sorted(tc_rows) != list(range(tcount)):
    fail(f"golden THERMOCHEM_DATA has {tc_len} rows ({len(tc_rows)} listed); KM_THERMOCHEM_COUNT is {tcount}")
if sh_len != scount or sorted(sh_rows) != list(range(scount)):
    fail(f"golden SHOMATE_DATA has {sh_len} rows ({len(sh_rows)} listed); KM_SHOMATE_COUNT is {scount}")
checks.append(("fmt_int(KM_THERMOCHEM_COUNT);", "i", str(tc_len), "KM_THERMOCHEM_COUNT"))
checks.append(("fmt_int(KM_SHOMATE_COUNT);", "i", str(sh_len), "KM_SHOMATE_COUNT"))
for k in range(tcount):
    formula, name, hf, gf, s0 = tc_rows[k]
    lab = f"THERMOCHEM_DATA[{k}]"
    checks.append((f"p = &_KM_THERMOCHEM_FORMULA_BYTES + _KM_THERMOCHEM_FORMULA_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(formula[2:]), f"{lab}.formula"))
    checks.append((f"p = &_KM_THERMOCHEM_NAME_BYTES + _KM_THERMOCHEM_NAME_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(name[2:]), f"{lab}.name"))
    for col, want in (("DELTA_HF_KJ", hf), ("DELTA_GF_KJ", gf), ("S_STANDARD_J", s0)):
        checks.append((f"x = KM_THERMOCHEM_{col}[{k}]; fmt_hex(x);", "f", want[2:], f"{lab}.{col.lower()}"))
for k in range(scount):
    formula, *fields = sh_rows[k]
    lab = f"SHOMATE_DATA[{k}]"
    checks.append((f"p = &_KM_SHOMATE_FORMULA_BYTES + _KM_SHOMATE_FORMULA_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(formula[2:]), f"{lab}.formula"))
    for col, want in zip(("A", "B", "C", "D", "E", "T_MIN", "T_MAX"), fields):
        checks.append((f"x = KM_SHOMATE_{col}[{k}]; fmt_hex(x);", "f", want[2:], f"{lab}.{col.lower()}"))
covered_t = {"_KM_THERMOCHEM_FORMULA_OFF", "_KM_THERMOCHEM_FORMULA_BYTES", "_KM_THERMOCHEM_NAME_OFF",
             "_KM_THERMOCHEM_NAME_BYTES", "KM_THERMOCHEM_DELTA_HF_KJ", "KM_THERMOCHEM_DELTA_GF_KJ",
             "KM_THERMOCHEM_S_STANDARD_J", "_KM_SHOMATE_FORMULA_OFF", "_KM_SHOMATE_FORMULA_BYTES",
             "KM_SHOMATE_A", "KM_SHOMATE_B", "KM_SHOMATE_C", "KM_SHOMATE_D", "KM_SHOMATE_E", "KM_SHOMATE_T_MIN",
             "KM_SHOMATE_T_MAX"}
tables_t = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", tdata, re.M))
if tables_t != covered_t:
    fail(f"src/thermochem_data.cyr tables not probed: {sorted(tables_t - covered_t)}; probed but gone: {sorted(covered_t - tables_t)}")
consts_t = set(re.findall(r"^const (_?KM_[A-Z0-9_]+)", tdata, re.M))
if consts_t != {"KM_THERMOCHEM_COUNT", "KM_SHOMATE_COUNT"}:
    fail(f"src/thermochem_data.cyr constants not probed: {sorted(consts_t)}")
nm2b = len(checks) - nm2b_before

# ── src/electrochemistry_data.cyr (M3a) ─────────────────────────────────────────────
nm3a_before = len(checks)
edata = (root / "src/electrochemistry_data.cyr").read_text()
sp_rows, sp_len = {}, None
for fn, cls, ins, outs in rows("electrochemistry"):
    if fn == "electrochemistry::STANDARD_POTENTIALS[]":
        sp_rows[int(ins[0][2:])] = outs
    elif fn == "electrochemistry::STANDARD_POTENTIALS":
        sp_len = int(outs[0][2:])
spcount = int(re.search(r"^const KM_STANDARD_POTENTIALS_COUNT = (\d+);", edata, re.M).group(1))
if sp_len != spcount or sorted(sp_rows) != list(range(spcount)):
    fail(f"golden STANDARD_POTENTIALS has {sp_len} rows ({len(sp_rows)} listed); KM_STANDARD_POTENTIALS_COUNT is {spcount}")
checks.append(("fmt_int(KM_STANDARD_POTENTIALS_COUNT);", "i", str(sp_len), "KM_STANDARD_POTENTIALS_COUNT"))
for k in range(spcount):
    equation, couple, electrons, e0 = sp_rows[k]
    lab = f"STANDARD_POTENTIALS[{k}]"
    checks.append((f"p = &_KM_STANDARD_POTENTIALS_EQUATION_BYTES + _KM_STANDARD_POTENTIALS_EQUATION_OFF[{k}]; "
                   f"print(p, strlen(p));", "s", unpct_utf8(equation[2:]), f"{lab}.equation (UTF-8)"))
    checks.append((f"p = &_KM_STANDARD_POTENTIALS_COUPLE_BYTES + _KM_STANDARD_POTENTIALS_COUPLE_OFF[{k}]; "
                   f"print(p, strlen(p));", "s", unpct(couple[2:]), f"{lab}.couple"))
    checks.append((f"fmt_int(KM_STANDARD_POTENTIALS_ELECTRONS[{k}]);", "i", electrons[2:], f"{lab}.electrons"))
    checks.append((f"x = KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL[{k}]; fmt_hex(x);", "f", e0[2:],
                   f"{lab}.standard_potential"))
covered_e = {"_KM_STANDARD_POTENTIALS_EQUATION_OFF", "_KM_STANDARD_POTENTIALS_EQUATION_BYTES",
             "_KM_STANDARD_POTENTIALS_COUPLE_OFF", "_KM_STANDARD_POTENTIALS_COUPLE_BYTES",
             "KM_STANDARD_POTENTIALS_ELECTRONS", "KM_STANDARD_POTENTIALS_STANDARD_POTENTIAL"}
tables_e = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", edata, re.M))
if tables_e != covered_e:
    fail(f"src/electrochemistry_data.cyr tables not probed: {sorted(tables_e - covered_e)}; probed but gone: {sorted(covered_e - tables_e)}")
consts_e = set(re.findall(r"^const (_?KM_[A-Z0-9_]+)", edata, re.M))
if consts_e != {"KM_STANDARD_POTENTIALS_COUNT"}:
    fail(f"src/electrochemistry_data.cyr constants not probed: {sorted(consts_e)}")
nm3a = len(checks) - nm3a_before

# ── src/organic_data.cyr (M4b) ──────────────────────────────────────────────────────
nm4b_before = len(checks)
odata = (root / "src/organic_data.cyr").read_text()
orows = {}
for fn, cls, ins, outs in rows("organic"):
    if fn in ("BondType::energy_kj", "Geometry::ideal_bond_angle", "FunctionalGroup::pattern",
              "FunctionalGroup::is_polar", "FunctionalGroup::hydrogen_bonding"):
        orows.setdefault(fn, []).append((ins[0][2:].split(".", 1)[1], outs[0]))
def member(prefix, variant):
    # ADR 0002 #6: Geometry::TShapedThreeLonePairs is KM_GEOMETRY_T_SHAPED.
    return "KM_GEOMETRY_T_SHAPED" if variant == "TShapedThreeLonePairs" else prefix + snake(variant)
for fn, prefix, count, col in (("BondType::energy_kj", "KM_BOND_TYPE_", "KM_BOND_TYPE_COUNT", "KM_BOND_TYPE_ENERGY_KJ"),
                               ("Geometry::ideal_bond_angle", "KM_GEOMETRY_", "KM_GEOMETRY_COUNT",
                                "KM_GEOMETRY_IDEAL_BOND_ANGLE")):
    rs = orows.get(fn, [])
    checks.append((f"fmt_int({count});", "i", str(len(rs)), count))
    for k, (v, want) in enumerate(rs):   # the golden rows are in declaration order: row k is value k
        name = member(prefix, v)
        checks.append((f"fmt_int({name});", "i", str(k), f"{name} = {k}"))
        checks.append((f"x = {col}[{name}]; fmt_hex(x);", "f", want[2:], f"{col}[{name}]"))
pats = orows.get("FunctionalGroup::pattern", [])
checks.append(("fmt_int(KM_FUNCTIONAL_GROUP_COUNT);", "i", str(len(pats)), "KM_FUNCTIONAL_GROUP_COUNT"))
for k, (v, want) in enumerate(pats):
    name = member("KM_FUNCTIONAL_GROUP_", v)
    checks.append((f"fmt_int({name});", "i", str(k), f"{name} = {k}"))
    checks.append((f"p = &_KM_FUNCTIONAL_GROUP_PATTERN_BYTES + _KM_FUNCTIONAL_GROUP_PATTERN_OFF[{name}]; print(p, strlen(p));",
                   "s", unpct_utf8(want[2:]), f"pattern[{name}] (UTF-8)"))
    # ADR 0002 #6's acceptor column has no golden row: its rule on the golden pattern, IUPAC 2011's electron-rich
    # region — a lone-pair atom (O, N, S, a halogen or the placeholder X) or a π bond (=, ≡, an aryl ring).
    pat = unpct_utf8(want[2:])
    atoms = set(re.findall(r"[A-Z][a-z]?", pat))
    site = atoms & {"O", "N", "S", "F", "Cl", "Br", "I", "X", "Ar"} or "=" in pat or "≡" in pat or "C₆H₅" in pat
    acc = "1" if site else "0"
    checks.append((f"if (KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_ACCEPTOR[{name}]) {{ fmt_int(1); }} else {{ fmt_int(0); }}",
                   "i", acc, f"acceptor[{name}] (by its rule)"))
for fn, col in (("FunctionalGroup::is_polar", "KM_FUNCTIONAL_GROUP_IS_POLAR"),
                ("FunctionalGroup::hydrogen_bonding", "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_DONOR")):
    for v, want in orows.get(fn, []):
        name = member("KM_FUNCTIONAL_GROUP_", v)
        checks.append((f"if ({col}[{name}]) {{ fmt_int(1); }} else {{ fmt_int(0); }}", "i", want[2:], f"{col}[{name}]"))
covered_o = {"KM_BOND_TYPE_ENERGY_KJ", "KM_GEOMETRY_IDEAL_BOND_ANGLE", "_KM_FUNCTIONAL_GROUP_PATTERN_OFF",
             "_KM_FUNCTIONAL_GROUP_PATTERN_BYTES", "KM_FUNCTIONAL_GROUP_IS_POLAR",
             "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_DONOR", "KM_FUNCTIONAL_GROUP_IS_HYDROGEN_BOND_ACCEPTOR"}
tables_o = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", odata, re.M))
if tables_o != covered_o:
    fail(f"src/organic_data.cyr tables not probed: {sorted(tables_o - covered_o)}; probed but gone: {sorted(covered_o - tables_o)}")
consts_o = set(re.findall(r"^const (_?KM_[A-Z0-9_]+)", odata, re.M))
if consts_o != {"KM_BOND_TYPE_COUNT", "KM_GEOMETRY_COUNT", "KM_FUNCTIONAL_GROUP_COUNT"}:
    fail(f"src/organic_data.cyr constants not probed: {sorted(consts_o)}")
members_o = set(re.findall(r"^    (KM_[A-Z0-9_]+) = \d+;", odata, re.M))
probed_o = {member(p, v) for fn, p in (("BondType::energy_kj", "KM_BOND_TYPE_"), ("Geometry::ideal_bond_angle", "KM_GEOMETRY_"),
                                        ("FunctionalGroup::pattern", "KM_FUNCTIONAL_GROUP_")) for v, _ in orows.get(fn, [])}
if members_o != probed_o:
    fail(f"src/organic_data.cyr enum members not probed: {sorted(members_o - probed_o)}; probed but gone: {sorted(probed_o - members_o)}")
nm4b_o = len(checks) - nm4b_before

# ── src/inorganic_data.cyr (M4b) ────────────────────────────────────────────────────
nm4b_before = len(checks)
idata = (root / "src/inorganic_data.cyr").read_text()
for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = ", idata, re.M):
    path, name = m.groups()
    if path not in golden:
        fail(f"src/inorganic_data.cyr: {name}: kimiya::{path} has no golden row")
        continue
    checks.append((f"x = {name}; fmt_hex(x);", "f", golden[path][2:], f"const {path}"))
    nconst += 1
lig, cfse, unp, radii, radii_len = [], {}, {}, {}, None
for fn, cls, ins, outs in rows("inorganic"):
    if fn == "Ligand::field_strength":
        lig.append((ins[0][2:].split(".", 1)[1], outs[0]))
    elif fn == "inorganic::cfse_octahedral" and outs[0].startswith("f:"):
        cfse[(int(ins[0][2:]), int(ins[1][2:]))] = outs[0]
    elif fn == "inorganic::unpaired_electrons_oct" and outs[0].startswith("i:"):
        unp[(int(ins[0][2:]), int(ins[1][2:]))] = outs[0]
    elif fn == "inorganic::IONIC_RADII[]":
        radii[int(ins[0][2:])] = outs
    elif fn == "inorganic::IONIC_RADII":
        radii_len = int(outs[0][2:])
checks.append(("fmt_int(KM_LIGAND_COUNT);", "i", str(len(lig)), "KM_LIGAND_COUNT"))
for k, (v, want) in enumerate(lig):
    name = "KM_LIGAND_" + snake(v)
    checks.append((f"fmt_int({name});", "i", str(k), f"{name} = {k}"))
    checks.append((f"fmt_int(KM_LIGAND_FIELD_STRENGTH[{name}]);", "i", want[2:], f"KM_LIGAND_FIELD_STRENGTH[{name}]"))
if sorted(cfse) != sorted(unp) or sorted(cfse) != [(d, h) for d in range(11) for h in (0, 1)]:
    fail(f"golden cfse / unpaired Ok rows are not d = 0 ..= 10 × both spins: {sorted(cfse)[:4]}..")
for (d, h), want in sorted(cfse.items()):
    spin = "HIGH_SPIN" if h == 1 else "LOW_SPIN"
    t2g, eg = f"_KM_OCT_OCCUPANCY_{spin}_T2G[{d}]", f"_KM_OCT_OCCUPANCY_{spin}_EG[{d}]"
    # Each occupancy cell is probed through what pins it: t2g + eg = d, and Rust's cfse formula on the pair.
    checks.append((f"fmt_int({t2g} + {eg});", "i", str(d), f"{spin} occupancy d{d}: t2g + eg"))
    checks.append((f"ft = f64_from({t2g}); fe = f64_from({eg}); x = ft * m04 + fe * 0.6; fmt_hex(x);", "f", want[2:],
                   f"{spin} occupancy d{d}: t2g·(−0.4) + eg·0.6 = cfse_octahedral({d}, {h})"))
    checks.append((f"fmt_int(_KM_OCT_UNPAIRED_{spin}[{d}]);", "i", unp[(d, h)][2:], f"_KM_OCT_UNPAIRED_{spin}[{d}]"))
rcount = int(re.search(r"^const KM_IONIC_RADII_COUNT = (\d+);", idata, re.M).group(1))
if radii_len != rcount or sorted(radii) != list(range(rcount)):
    fail(f"golden IONIC_RADII has {radii_len} rows ({len(radii)} listed); KM_IONIC_RADII_COUNT is {rcount}")
checks.append(("fmt_int(KM_IONIC_RADII_COUNT);", "i", str(radii_len), "KM_IONIC_RADII_COUNT"))
for k in range(rcount):
    sym, q, cn, r = radii[k]
    lab = f"IONIC_RADII[{k}]"
    checks.append((f"p = &_KM_IONIC_RADII_SYMBOL_BYTES + _KM_IONIC_RADII_SYMBOL_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(sym[2:]), f"{lab}.symbol"))
    checks.append((f"fmt_int(KM_IONIC_RADII_CHARGE[{k}]);", "i", q[2:], f"{lab}.charge (i8)"))
    checks.append((f"fmt_int(KM_IONIC_RADII_COORDINATION[{k}]);", "i", cn[2:], f"{lab}.coordination"))
    checks.append((f"x = KM_IONIC_RADII_RADIUS_PM[{k}]; fmt_hex(x);", "f", r[2:], f"{lab}.radius_pm"))
covered_i = {"KM_LIGAND_FIELD_STRENGTH", "_KM_OCT_OCCUPANCY_HIGH_SPIN_T2G", "_KM_OCT_OCCUPANCY_HIGH_SPIN_EG",
             "_KM_OCT_OCCUPANCY_LOW_SPIN_T2G", "_KM_OCT_OCCUPANCY_LOW_SPIN_EG", "_KM_OCT_UNPAIRED_HIGH_SPIN",
             "_KM_OCT_UNPAIRED_LOW_SPIN", "_KM_IONIC_RADII_SYMBOL_OFF", "_KM_IONIC_RADII_SYMBOL_BYTES",
             "KM_IONIC_RADII_CHARGE", "KM_IONIC_RADII_COORDINATION", "KM_IONIC_RADII_RADIUS_PM"}
tables_i = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", idata, re.M))
if tables_i != covered_i:
    fail(f"src/inorganic_data.cyr tables not probed: {sorted(tables_i - covered_i)}; probed but gone: {sorted(covered_i - tables_i)}")
consts_i = set(re.findall(r"^const (_?KM_[A-Z0-9_]+) = \d+;", idata, re.M))
if consts_i != {"KM_LIGAND_COUNT", "KM_IONIC_RADII_COUNT"}:
    fail(f"src/inorganic_data.cyr integer constants not probed: {sorted(consts_i)}")
members_i = set(re.findall(r"^    (KM_[A-Z0-9_]+) = \d+;", idata, re.M))
if members_i != {"KM_LIGAND_" + snake(v) for v, _ in lig}:
    fail(f"src/inorganic_data.cyr enum members not probed: {sorted(members_i)}")
nm4b_i = len(checks) - nm4b_before

# ── src/nuclear_data.cyr (M4c) ──────────────────────────────────────────────────────
nm4c_before = len(checks)
ndata = (root / "src/nuclear_data.cyr").read_text()
iso_rows, iso_len = {}, None
for fn, cls, ins, outs in rows("nuclear"):
    if fn == "nuclear::ISOTOPES[]":
        iso_rows[int(ins[0][2:])] = outs
    elif fn == "nuclear::ISOTOPES":
        iso_len = int(outs[0][2:])
icount_n = int(re.search(r"^const KM_ISOTOPES_COUNT = (\d+);", ndata, re.M).group(1))
if iso_len != icount_n or sorted(iso_rows) != list(range(icount_n)):
    fail(f"golden ISOTOPES has {iso_len} rows ({len(iso_rows)} listed); KM_ISOTOPES_COUNT is {icount_n}")
checks.append(("fmt_int(KM_ISOTOPES_COUNT);", "i", str(iso_len), "KM_ISOTOPES_COUNT"))
for k in range(icount_n):
    sym, z, a, hl, mode = iso_rows[k]
    lab = f"ISOTOPES[{k}]"
    checks.append((f"p = &_KM_ISOTOPES_SYMBOL_BYTES + _KM_ISOTOPES_SYMBOL_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct(sym[2:]), f"{lab}.symbol"))
    checks.append((f"fmt_int(KM_ISOTOPES_Z[{k}]);", "i", z[2:], f"{lab}.z (u8)"))
    checks.append((f"fmt_int(KM_ISOTOPES_A[{k}]);", "i", a[2:], f"{lab}.a (u16)"))
    checks.append((f"x = KM_ISOTOPES_HALF_LIFE_S[{k}]; fmt_hex(x);", "f", hl[2:], f"{lab}.half_life_s"))
    checks.append((f"p = &_KM_ISOTOPES_DECAY_MODE_BYTES + _KM_ISOTOPES_DECAY_MODE_OFF[{k}]; print(p, strlen(p));",
                   "s", unpct_utf8(mode[2:]), f"{lab}.decay_mode (UTF-8)"))
covered_n = {"_KM_ISOTOPES_SYMBOL_OFF", "_KM_ISOTOPES_SYMBOL_BYTES", "KM_ISOTOPES_Z", "KM_ISOTOPES_A",
             "KM_ISOTOPES_HALF_LIFE_S", "_KM_ISOTOPES_DECAY_MODE_OFF", "_KM_ISOTOPES_DECAY_MODE_BYTES"}
tables_n = set(re.findall(r"^var (_?KM_[A-Z0-9_]+)", ndata, re.M))
if tables_n != covered_n:
    fail(f"src/nuclear_data.cyr tables not probed: {sorted(tables_n - covered_n)}; probed but gone: {sorted(covered_n - tables_n)}")
consts_n = set(re.findall(r"^const (_?KM_[A-Z0-9_]+)", ndata, re.M))
if consts_n != {"KM_ISOTOPES_COUNT"}:
    fail(f"src/nuclear_data.cyr constants not probed: {sorted(consts_n)}")
nm4c = len(checks) - nm4c_before

# ── the probe program ───────────────────────────────────────────────────────────────
src = ['include "src/constants.cyr"', 'include "src/element_data.cyr"', 'include "src/phase_data.cyr"',
       'include "src/gas_data.cyr"', 'include "src/thermochem_data.cyr"',
       'include "src/electrochemistry_data.cyr"', 'include "src/organic_data.cyr"',
       'include "src/inorganic_data.cyr"', 'include "src/nuclear_data.cyr"', "",
       "fn main() {", "    var x: f64 = 0.0;", "    var m04: f64 = -0.4;", "    var ft: f64 = 0.0;",
       "    var fe: f64 = 0.0;", "    var p = 0;"]
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
print(f"checked {nconst} constants, {nelem} element_data cells, {nm2a} phase_data / gas_data cells, "
      f"{nm2b} thermochem_data cells, {nm3a} electrochemistry_data cells, {nm4b_o} organic_data cells, "
      f"{nm4b_i} inorganic_data cells and {nm4c} nuclear_data cells through the Cyrius compiler ({pending} golden "
      f"constants without a generated "
      f"home); failures: {fails}")
if pending:
    fail(f"{pending} golden constants have no generated constant (none is pending since M4b)")
sys.exit(1 if fails else 0)
