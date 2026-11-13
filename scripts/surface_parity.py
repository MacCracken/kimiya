"""Surface parity: every public item of the ported Rust modules has a Cyrius counterpart.

Port plan §4 (a per-milestone gate) and M7 (the release gate). Adapted from tanmatra's
scripts/surface_parity.py. Usage, from the repo root:
    python3 -I scripts/surface_parity.py .

For each ported module it lists rust-old/src/<module>.rs's public items before `#[cfg(test)]`
(top level only; `pub` and, in validate.rs, `pub(crate)`) and checks that src/ defines the
Cyrius name the porting rules give (CLAUDE.md, port plan §6.4):
- free `pub fn f`                    -> `fn <prefix>f(`, the module's prefix (`element_`) where the
                                        bare Rust name would be ambiguous in one flat namespace
- `pub(crate) fn f` (validate)       -> `fn _km_f(` (private helpers carry `_km_`)
- `impl Struct { pub fn m }`         -> `fn m(` inside `impl Km<Struct> {`
- `impl Enum { pub const fn m }`     -> a free `fn <prefix>m(` of the enum value, the prefix ENUM_METHODS
                                        names (`bond_type_energy_kj(b)`, M4b)
- `pub const X`                      -> `const KM_X =`, or for a struct-valued const (gas::PR_*) a
                                        top-level struct value `var KM_X = Km<Struct> {`
- `pub enum E` and each variant V    -> `enum Km<E> {` and `KM_<SNAKE(E)>_<SNAKE(V)> =`
                                        (KimiyaError: `enum KmError`, `KM_ERR_<SNAKE(V)>`)
- `pub struct S`                     -> `struct Km<S> {`, or a handle with accessors (STRUCT_ACCESSORS)
- a `pub` field f of a struct        -> the accessor STRUCT_ACCESSORS names, or a field of `Km<S>`
A Rust item ADR 0002 renames needs an entry in RENAMED (its Cyrius name, which must exist); anything else
needs an entry in EXCEPTIONS with its reason (ported in another form). Every file under rust-old/src
must be in MODULES, PENDING_MODULES (a later milestone), DEFERRED_MODULES or NOT_APPLICABLE_MODULES.

The whole surface (M7, the release gate): every public item of every file under rust-old/src, the
deferred modules' item by item, lib.rs's `pub mod` declarations and `pub use` re-exports (each of which
must name a public item of its module), and the trait impls the types derive. Every `pub` line before
`#[cfg(test)]` must be one of the items counted. Each item is ported (by its Cyrius name, or in another
form), renamed by ADR 0002, n/a, or deferred (port plan §2.2). Exits 1 on any item missing, an
EXCEPTIONS or RENAMED entry that names no Rust item, a `pub` line no item accounts for, or a file
unaccounted for. `--markdown` prints the whole-surface table instead (docs/development/rust-parity.md).
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
markdown = "--markdown" in sys.argv[2:]

# The ported modules: the free-fn prefix and which Rust structs keep a Cyrius struct.
MODULES = {
    "error": {"prefix": "", "structs": set()},
    "validate": {"prefix": "_km_", "structs": set()},
    "element": {"prefix": "element_", "structs": set()},
    "molecule": {"prefix": "", "structs": {"Molecule"}},
    "reaction": {"prefix": "", "structs": set()},
    "thermo": {"prefix": "", "structs": set()},
    "thermochem": {"prefix": "", "structs": {"Shomate"}},
    "phase": {"prefix": "", "structs": {"AntoineCoeffs"}},
    "gas": {"prefix": "", "structs": {"PengRobinsonParams"}},
    "solution": {"prefix": "", "structs": set()},
    "electrochemistry": {"prefix": "", "structs": set()},
    "kinetics": {"prefix": "", "structs": set()},
    "organic": {"prefix": "", "structs": set()},
    "inorganic": {"prefix": "", "structs": set()},
    "nuclear": {"prefix": "", "structs": set()},
    "potential": {"prefix": "", "structs": set()},
    "spectroscopy": {"prefix": "", "structs": set()},
    "stoichiometry": {"prefix": "", "structs": set()},    # M4f: ADR 0002 #7's integer matrix and rational states
    "reaction_dynamics": {"prefix": "", "structs": set()},  # M5a: closures → fn pointer + context, flat trajectories
    "network": {"prefix": "", "structs": set()},            # M5b: flat matrices, 32-byte reaction records
    "fitting": {"prefix": "", "structs": {"ArrheniusFit", "PolynomialFit"}},  # M5b: KmPolynomialFit over a buffer
    "sensitivity": {"prefix": "", "structs": set()},        # M5b: ADR 0002 #2; compute_gradient on hisab 3.4.1's tape
}

# Rust enums with `impl` methods (M4b): their const property methods are free functions of the enum value,
# `<prefix><method>(v)` (port plan §6.4).
ENUM_METHODS = {
    "BondType": "bond_type_",
    "Geometry": "geometry_",
    "FunctionalGroup": "functional_group_",
    "Ligand": "ligand_",
    "SpectralSeries": "spectral_series_",   # M4e: spectral_series_n_lower(s), spectral_series_wavelength(s, line)
}

# Rust structs that became handles (port plan §6.4: `Option<&'static T>` → a handle, 0 for None):
# struct -> the accessor prefix of its fields, with the handle's documented range.
STRUCT_ACCESSORS = {
    "Element": "element_",                    # handle = Z (1–118)
    "ElementProperties": "element_properties_",  # handle = property row + 1 (1–32)
    "Atom": "atom_",                          # entry k of a KmMolecule: `m.atom_<field>(k)`
    "ThermochemData": "thermochem_data_",     # handle = THERMOCHEM_DATA row + 1 (1–76)
    "HalfReaction": "half_reaction_",         # handle = STANDARD_POTENTIALS row + 1 (1–49)
    "IonicRadius": "ionic_radius_",           # handle = IONIC_RADII row + 1 (1–30)
    "Isotope": "isotope_",                    # handle = ISOTOPES row + 1 (1–18)
}

# Rust item -> why it has no direct Cyrius name (it is still checked to exist in Rust).
EXCEPTIONS = {
    "error::type Result": "the stdlib value-form Result, Err payload a KM_ERR_* code (port plan §6.3)",
    "element::static ELEMENTS": "the generated columns KM_ELEMENT_* / _KM_ELEMENT_*_BYTES in "
                                "src/element_data.cyr (row Z − 1, KM_ELEMENT_COUNT rows), read through "
                                "element_lookup_by_number and the element_* accessors",
    "element::static ELEMENT_PROPERTIES": "the generated columns KM_ELEMENT_PROPERTIES_* in "
                                          "src/element_data.cyr (KM_ELEMENT_PROPERTIES_COUNT rows), read "
                                          "through element_lookup_properties and the element_properties_* "
                                          "accessors",
    "element::ElementProperties.oxidation_states": "a slice: element_properties_oxidation_state_count(p) + "
                                                   "element_properties_oxidation_state(p, k) (Option)",
    "molecule::Molecule.atoms": "the `atoms` vec of struct KmMolecule (2 words per atom); "
                                "KmMolecule.atoms_len() is `atoms.len()`",
    "element::const AVOGADRO": "KM_AVOGADRO in src/constants.cyr (generated; 14 modules share the constants)",
    "element::const GAS_CONSTANT": "KM_GAS_CONSTANT in src/constants.cyr (generated)",
    "phase::static ANTOINE_DATA": "the generated columns KM_ANTOINE_A / _B / _C / _T_MIN_C / _T_MAX_C and the "
                                  "name blob _KM_ANTOINE_NAME_BYTES in src/phase_data.cyr (KM_ANTOINE_COUNT "
                                  "rows), read through handles (row + 1): lookup_antoine, antoine_data_name, "
                                  "antoine_data_coeffs (port plan §6.4)",
    "thermochem::static THERMOCHEM_DATA": "the generated columns KM_THERMOCHEM_DELTA_HF_KJ / _DELTA_GF_KJ / "
                                          "_S_STANDARD_J and the formula / name blobs in src/thermochem_data.cyr "
                                          "(KM_THERMOCHEM_COUNT rows), read through handles (row + 1): "
                                          "lookup_thermochem and the thermochem_data_* accessors (port plan §6.4)",
    "thermochem::static SHOMATE_DATA": "the generated columns KM_SHOMATE_A … _T_MAX and the formula blob in "
                                       "src/thermochem_data.cyr (KM_SHOMATE_COUNT rows), read through handles "
                                       "(row + 1): lookup_shomate, lookup_shomate_at, shomate_data_formula, "
                                       "shomate_data_coeffs (port plan §6.4)",
    "electrochemistry::static STANDARD_POTENTIALS": "the generated columns KM_STANDARD_POTENTIALS_ELECTRONS / "
                                                    "_STANDARD_POTENTIAL and the UTF-8 equation and couple blobs in "
                                                    "src/electrochemistry_data.cyr (KM_STANDARD_POTENTIALS_COUNT "
                                                    "rows), read through handles (row + 1): lookup_half_reaction and "
                                                    "the half_reaction_* accessors (port plan §6.4)",
    "inorganic::static IONIC_RADII": "the generated columns KM_IONIC_RADII_CHARGE / _COORDINATION / _RADIUS_PM and the "
                                     "symbol blob _KM_IONIC_RADII_SYMBOL_BYTES in src/inorganic_data.cyr "
                                     "(KM_IONIC_RADII_COUNT rows), read through handles (row + 1): lookup_ionic_radius and "
                                     "the ionic_radius_* accessors (port plan §6.4)",
    "reaction_dynamics::type TrajectoryPoint": "a point of the caller's flat trajectory buffer: 1 + n f64 (t, then the "
                                               "n concentrations), point k at out + 8·k·(1 + n) (port plan §4, M5)",
    "network::type ReactionDef": "a 32-byte reaction record (the reactants' pointer and count, the products' pointer and "
                                 "count), each term a 16-byte (species index word, coefficient f64) record "
                                 "(src/network.cyr; port plan §6.4, M5b)",
    "nuclear::static ISOTOPES": "the generated columns KM_ISOTOPES_Z / _A / _HALF_LIFE_S, the symbol blob and the UTF-8 "
                                "decay-mode blob in src/nuclear_data.cyr (KM_ISOTOPES_COUNT rows), read through handles "
                                "(row + 1): lookup_isotope and the isotope_* accessors (port plan §6.4)",
}

# Rust item -> (its Cyrius name, which must be defined in src/, and the ADR 0002 row that renames it).
RENAMED = {
    "organic::Geometry::TShapedThreeLonePairs": ("KM_GEOMETRY_T_SHAPED", "ADR 0002 #6, at its Rust position (7): the Rust "
                                                 "name is a misnomer, ClF₃ has two lone pairs"),
    "organic::FunctionalGroup::hydrogen_bonding": ("functional_group_is_hydrogen_bond_donor", "ADR 0002 #6, with the new "
                                                   "functional_group_is_hydrogen_bond_acceptor beside it"),
}

# The whole surface (M7): the reasons for what is not a module's own item.
LIB_MOD_NA = "n/a: a module declaration; Cyrius has one flat namespace, the module is src/<module>.cyr in [lib]"
LIB_USE_NA = "n/a: a re-export (the prelude) of a public item counted in its own module; one flat namespace"
# Derived trait -> its status in the port.
DERIVE_STATUS = {
    "Debug": "n/a: Rust formatting for `{:?}`; a Cyrius value is a word, a struct or a handle",
    "Clone": "n/a: a Cyrius struct or handle is copied by value",
    "Copy": "n/a: a Cyrius enum value or struct is copied by value",
    "PartialEq": "n/a: a Cyrius enum value or handle is a word, compared with ==",
    "Eq": "n/a: as PartialEq",
    "Hash": "n/a: no Cyrius hashed collection holds the enums",
    "Serialize": "n/a: serde is not carried over (port plan §2.4)",
    "Deserialize": "n/a: serde is not carried over (port plan §2.4)",
    "Error": "ported in another form: thiserror's Display text without its payload is km_error_message, the variant "
             "km_error_name (src/error.cyr; port plan §6.4)",
}

# Later milestones (port plan §2.1): counted, not checked.
PENDING_MODULES = {}
DEFERRED_MODULES = {
    "bridge": "deferred until a Cyrius consumer exists (port plan §2.2)",
    "integration/soorat": "deferred: waits on soorat (port plan §2.2)",
    "logging": "deferred: wraps tracing-subscriber, no log events (port plan §2.2)",
}
NOT_APPLICABLE_MODULES = {
    "lib": "n/a: module declarations and the prelude's re-exports; Cyrius has one flat namespace",
    "integration/mod": "n/a: the `soorat` module declaration",
}


def snake(camel):
    """scripts/oracle/src/cyr.rs's snake (tanmatra's rule): `_` before a capital after a lowercase letter or digit, or
    before a capital that starts a capitalised word after a run of capitals (`CCDouble` → `CC_DOUBLE`, `HCl` → `H_CL`;
    M4b: the regex this replaced gave `CCDOUBLE`, which no generated name uses)."""
    out = []
    for k, ch in enumerate(camel):
        if ch.isupper() and k > 0:
            prev, nxt = camel[k - 1], camel[k + 1] if k + 1 < len(camel) else ""
            if prev.islower() or prev.isdigit() or (prev.isupper() and nxt.islower()):
                out.append("_")
        out.append(ch.upper())
    return "".join(out)


def rust_source(module):
    """rust-old/src/<module>.rs up to its #[cfg(test)]."""
    src = (root / f"rust-old/src/{module}.rs").read_text()
    cut = src.find("#[cfg(test)]")
    return src if cut < 0 else src[:cut]


def rust_items(module, with_lines=False):
    """(label, kind, owner, name) for every public item before #[cfg(test)]; with_lines adds its line number."""
    src = rust_source(module)
    vis = r"pub(?:\(crate\))?"
    items = []
    depth = 0
    impl = None
    block = None  # ("struct" | "enum", name) while inside its body at depth 1
    for lineno, line in enumerate(src.splitlines(), 1):
        s = line.strip()
        before = len(items)
        if depth == 0:
            m = re.match(r"impl(?:<[^>]*>)?\s+(?:\w+(?:::\w+)*\s+for\s+)?(\w+)\s*\{", s)
            if m:
                impl = m.group(1)
            m = re.match(rf"{vis} (?:const )?fn (\w+)", s)
            if m:
                items.append((f"{module}::{m.group(1)}", "fn", None, m.group(1)))
            m = re.match(rf"{vis} const ([A-Z][A-Z0-9_]*)\s*:", s)
            if m:
                items.append((f"{module}::const {m.group(1)}", "const", None, m.group(1)))
            m = re.match(rf"{vis} static ([A-Z][A-Z0-9_]*)\s*:", s)
            if m:
                items.append((f"{module}::static {m.group(1)}", "static", None, m.group(1)))
            m = re.match(rf"{vis} type (\w+)", s)
            if m:
                items.append((f"{module}::type {m.group(1)}", "type", None, m.group(1)))
            m = re.match(rf"{vis} (struct|enum) (\w+)", s)
            if m:
                items.append((f"{module}::{m.group(2)}", m.group(1), None, m.group(2)))
                if s.endswith("{"):
                    block = (m.group(1), m.group(2))
        elif depth == 1:
            if impl:
                m = re.match(r"pub (?:const )?fn (\w+)", s)
                if m:
                    items.append((f"{module}::{impl}::{m.group(1)}", "method", impl, m.group(1)))
            elif block and block[0] == "struct":
                m = re.match(r"pub (\w+)\s*:", s)
                if m:
                    items.append((f"{module}::{block[1]}.{m.group(1)}", "field", block[1], m.group(1)))
            elif block and block[0] == "enum":
                m = re.match(r"([A-Z]\w*)\s*(?:\(|,|\{|$)", s)
                if m and not s.startswith("#"):
                    items.append((f"{module}::{block[1]}::{m.group(1)}", "variant", block[1], m.group(1)))
        if with_lines:
            items[before:] = [it + (lineno,) for it in items[before:]]
        depth += line.count("{") - line.count("}")
        if depth == 0:
            impl = None
            block = None
    return items


cyr = "\n".join(p.read_text() for p in sorted((root / "src").glob("*.cyr")))


def impl_methods(struct):
    m = re.search(rf"impl Km{struct} \{{(.*?)\n\}}", cyr, re.S)
    return set(re.findall(r"^\s+fn (\w+)\(", m.group(1), re.M)) if m else set()


def struct_fields(struct):
    m = re.search(rf"^struct Km{struct} \{{(.*?)\n\}}", cyr, re.S | re.M)
    return set(re.findall(r"^\s+(\w+)\s*[;:]", m.group(1), re.M)) if m else set()


def has_counterpart(module, kind, owner, name):
    kinds = MODULES[module]
    if kind == "fn":
        return re.search(rf"^fn {kinds['prefix']}{name}\(", cyr, re.M) is not None
    if kind == "const":
        return (re.search(rf"^const KM_{name} =", cyr, re.M) is not None
                or re.search(rf"^var KM_{name} = Km\w+ \{{", cyr, re.M) is not None)
    if kind == "enum":
        cname = "KmError" if name == "KimiyaError" else f"Km{name}"
        return re.search(rf"^enum {cname} \{{", cyr, re.M) is not None
    if kind == "variant":
        member = f"KM_ERR_{snake(name)}" if owner == "KimiyaError" else f"KM_{snake(owner)}_{snake(name)}"
        return re.search(rf"\b{member} =", cyr) is not None
    if kind == "struct":
        if name in kinds["structs"]:
            return re.search(rf"^struct Km{name} \{{", cyr, re.M) is not None
        prefix = STRUCT_ACCESSORS.get(name)
        return prefix is not None and re.search(rf"\bfn {prefix}\w+\(", cyr) is not None
    if kind == "field":
        if owner in kinds["structs"]:
            return name in struct_fields(owner)
        prefix = STRUCT_ACCESSORS.get(owner)
        return prefix is not None and re.search(rf"\bfn {prefix}{name}\(", cyr) is not None
    if kind == "method":
        if owner in ENUM_METHODS:
            return re.search(rf"^fn {ENUM_METHODS[owner]}{name}\(", cyr, re.M) is not None
        return name in impl_methods(owner)
    return False


missing = 0
unaccounted = 0   # public items (and `pub` lines) with no status: each is also counted in missing
out = []          # the report's lines (printed unless --markdown)
total_found = total_items = 0
all_labels = set()
whole = []        # (Rust file, label, status, detail) for the whole surface
for module in MODULES:
    found = 0
    items = rust_items(module)
    for label, kind, owner, name in items:
        all_labels.add(label)
        if has_counterpart(module, kind, owner, name):
            found += 1
            whole.append((module, label, "by name", ""))
        elif label in RENAMED:
            cname, why = RENAMED[label]
            if re.search(rf"(^fn {cname}\(|\b{cname} =)", cyr, re.M) is None:
                missing += 1
                unaccounted += 1
                out.append(f"MISSING  {label}: RENAMED names {cname}, which src/ does not define")
                continue
            found += 1
            whole.append((module, label, "renamed", f"`{cname}` ({why})"))
        elif label in EXCEPTIONS:
            found += 1
            whole.append((module, label, "another form", EXCEPTIONS[label]))
        else:
            missing += 1
            unaccounted += 1
            out.append(f"MISSING  {label} ({kind})")
    total_found += found
    total_items += len(items)
    out.append(f"{module:<12} {found}/{len(items)} Rust public items have a Cyrius counterpart or a reason")
for table, what in ((EXCEPTIONS, "EXCEPTIONS"), (RENAMED, "RENAMED")):
    for label in table:
        if label not in all_labels:
            missing += 1
            out.append(f"STALE    {what} names {label}, which is no public item of a ported module")
# Every file under rust-old/src must be ported, pending, deferred or not applicable.
known = set(MODULES) | set(PENDING_MODULES) | set(DEFERRED_MODULES) | set(NOT_APPLICABLE_MODULES)
files = []
for path in sorted((root / "rust-old/src").rglob("*.rs")):
    mod = str(path.relative_to(root / "rust-old/src"))[:-3]
    files.append(mod)
    if mod not in known:
        missing += 1
        out.append(f"UNKNOWN  rust-old/src/{mod}.rs is not ported, pending, deferred or marked not applicable")
for mod in sorted(set(PENDING_MODULES) | set(DEFERRED_MODULES) | set(NOT_APPLICABLE_MODULES)):
    if not (root / f"rust-old/src/{mod}.rs").exists():
        missing += 1
        out.append(f"STALE    {mod}: listed, but rust-old/src/{mod}.rs does not exist")
pend = {}
for mod, ms in PENDING_MODULES.items():
    pend.setdefault(ms, []).append(mod)
for ms in sorted(pend):
    n = sum(len(rust_items(m)) for m in pend[ms])
    out.append(f"pending {ms}: {len(pend[ms])} modules, {n} Rust public items ({', '.join(sorted(pend[ms]))})")
for mod, why in {**DEFERRED_MODULES, **NOT_APPLICABLE_MODULES}.items():
    out.append(f"{mod:<20} {len(rust_items(mod))} Rust public items, not ported: {why}")

# ── The whole surface (M7) ──────────────────────────────────────────────────────────────────────────────────────────
# The deferred modules, item by item.
for mod in DEFERRED_MODULES:
    for label, kind, owner, name in rust_items(mod):
        whole.append((mod, label, "deferred", DEFERRED_MODULES[mod]))
# The module declarations and re-exports of lib.rs and integration/mod.rs: lines, not items of a module.
decl_lines = {}   # module -> the line numbers of its `pub mod` / `pub use` lines
labels_of = {}    # module -> the labels of its items
for mod in files:
    if mod in known and mod not in NOT_APPLICABLE_MODULES:
        labels_of[mod] = {it[0] for it in rust_items(mod)}
for mod in NOT_APPLICABLE_MODULES:
    decl_lines[mod] = set()
    for lineno, line in enumerate(rust_source(mod).splitlines(), 1):
        s = line.strip()
        m = re.match(r"pub mod (\w+);", s)
        if m:
            decl_lines[mod].add(lineno)
            child = m.group(1) if mod == "lib" else f"{mod.rsplit('/', 1)[0]}/{m.group(1)}"
            children = [d for d in DEFERRED_MODULES if d == child or d.startswith(child + "/")]
            if children:
                whole.append((mod, f"{mod}::mod {m.group(1)}", "deferred", "the module declaration of "
                              + ", ".join(children) + " (port plan §2.2)"))
            elif child in MODULES:
                whole.append((mod, f"{mod}::mod {m.group(1)}", "n/a", LIB_MOD_NA))
            else:
                missing += 1
                out.append(f"UNKNOWN  rust-old/src/{mod}.rs: `pub mod {m.group(1)}` declares no known module")
            continue
        m = re.match(r"pub use (\w+)::(?:\{([^}]*)\}|(\w+));", s)
        if m:
            decl_lines[mod].add(lineno)
            target = m.group(1)
            for name in [n.strip() for n in (m.group(2) or m.group(3)).split(",") if n.strip()]:
                cands = {f"{target}::{name}", f"{target}::const {name}", f"{target}::static {name}",
                         f"{target}::type {name}"}
                if not cands & labels_of.get(target, set()):
                    missing += 1
                    unaccounted += 1
                    out.append(f"MISSING  rust-old/src/{mod}.rs re-exports {target}::{name}, no public item of {target}")
                    continue
                status = "deferred" if target in DEFERRED_MODULES else "n/a"
                whole.append((mod, f"{mod}::use {target}::{name}", status, LIB_USE_NA))
# Every `pub` line is one of the items counted (the parser misses nothing).
for mod in files:
    if mod not in known:
        continue
    seen = decl_lines.get(mod, set()) | {it[4] for it in rust_items(mod, with_lines=True)}
    for lineno, line in enumerate(rust_source(mod).splitlines(), 1):
        if re.match(r"\s*pub(\([\w:]+\))?\s", line) and lineno not in seen:
            missing += 1
            unaccounted += 1
            out.append(f"UNCOUNTED rust-old/src/{mod}.rs:{lineno}: {line.strip()}")
# The trait impls the types derive.
derived = []      # (Rust file, type, trait, status)
for mod in files:
    if mod not in known:
        continue
    lines = rust_source(mod).splitlines()
    for k, line in enumerate(lines):
        m = re.match(r"\s*#\[derive\(([^)]*)\)\]", line)
        if not m:
            continue
        ty = None
        for nxt in lines[k + 1 : k + 6]:
            t = re.match(r"\s*pub(?:\([\w:]+\))? (?:struct|enum) (\w+)", nxt)
            if t:
                ty = t.group(1)
                break
        if ty is None:
            missing += 1
            out.append(f"UNKNOWN  rust-old/src/{mod}.rs:{k + 1}: a derive on no public type")
            continue
        for trait in [x.strip() for x in m.group(1).split(",") if x.strip()]:
            if mod in DEFERRED_MODULES:
                derived.append((mod, ty, trait, "deferred"))
            elif trait in DERIVE_STATUS:
                derived.append((mod, ty, trait, DERIVE_STATUS[trait].split(":")[0]))
            else:
                missing += 1
                out.append(f"UNKNOWN  rust-old/src/{mod}.rs: {ty} derives {trait}, which DERIVE_STATUS does not list")

STATUSES = ("by name", "another form", "renamed", "n/a", "deferred")
counts = {st: sum(1 for w in whole if w[2] == st) for st in STATUSES}
n_whole = len(whole)
d_counts = {}
for _, _, _, st in derived:
    d_counts[st] = d_counts.get(st, 0) + 1
n_types = len({(m, t) for m, t, _, _ in derived})
out.append(f"whole surface: {n_whole} public items in rust-old/src: {counts['by name'] + counts['another form']} ported "
           f"({counts['by name']} by name, {counts['another form']} in another form), {counts['renamed']} renamed by "
           f"ADR 0002, {counts['n/a']} n/a (lib.rs's module declarations and re-exports), {counts['deferred']} deferred "
           f"(port plan §2.2)")
out.append(f"derived trait impls: {len(derived)} on {n_types} types: "
           + ", ".join(f"{v} {k}" for k, v in sorted(d_counts.items())))
out.append(f"surface parity: {total_found}/{total_items} in the ported modules; whole surface {n_whole}/{n_whole + unaccounted} "
           f"accounted for; {missing} missing")

if not markdown:
    print("\n".join(out))
else:
    md = []
    md.append("| Rust file | Public items | Ported by name | Ported in another form | Renamed (ADR 0002) | N/A | "
              "Deferred (§2.2) |")
    md.append("|---|---:|---:|---:|---:|---:|---:|")
    for mod in files:
        rows = [w for w in whole if w[0] == mod]
        c = {st: sum(1 for w in rows if w[2] == st) for st in STATUSES}
        md.append(f"| `{mod}.rs` | {len(rows)} | {c['by name']} | {c['another form']} | {c['renamed']} | {c['n/a']} | "
                  f"{c['deferred']} |")
    md.append(f"| **Total** | **{n_whole}** | **{counts['by name']}** | **{counts['another form']}** | "
              f"**{counts['renamed']}** | **{counts['n/a']}** | **{counts['deferred']}** |")
    md.append("")
    md.append("Renamed by ADR 0002:")
    md.append("")
    for mod, label, st, detail in whole:
        if st == "renamed":
            md.append(f"- `{label}` → {detail}.")
    md.append("")
    md.append("Ported in another form (the port's conventions, port plan §6.3–§6.4):")
    md.append("")
    md.append("| Rust item | Cyrius form |")
    md.append("|---|---|")
    for mod, label, st, detail in whole:
        if st == "another form":
            md.append(f"| `{label}` | {detail} |")
    md.append("")
    md.append(f"Derived trait impls (not items; {len(derived)} on {n_types} types):")
    md.append("")
    md.append("| Trait | Impls | Status |")
    md.append("|---|---:|---|")
    traits = sorted({t for _, _, t, _ in derived}, key=lambda t: list(DERIVE_STATUS).index(t))
    for trait in traits:
        live = sum(1 for m, _, t, st in derived if t == trait and st != "deferred")
        dfr = sum(1 for m, _, t, st in derived if t == trait and st == "deferred")
        note = DERIVE_STATUS[trait] + (f" ({dfr} of them on soorat's types: deferred)" if dfr else "")
        md.append(f"| `{trait}` | {live + dfr} | {note} |")
    print("\n".join(md))
    if missing:
        print("\n".join(o for o in out if not o.startswith((" ",)) and o.split()[0] in
                        ("MISSING", "STALE", "UNKNOWN", "UNCOUNTED")), file=sys.stderr)
sys.exit(1 if missing else 0)
