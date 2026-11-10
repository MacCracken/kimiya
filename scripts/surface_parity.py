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
Anything else needs an entry in EXCEPTIONS with its reason. Every file under rust-old/src must be
in MODULES, PENDING_MODULES (a later milestone), DEFERRED_MODULES or NOT_APPLICABLE_MODULES. Exits 1
on any item missing, an EXCEPTIONS entry that names no Rust item, or a file unaccounted for.
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()

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
    "organic::Geometry::TShapedThreeLonePairs": "renamed KM_GEOMETRY_T_SHAPED by ADR 0002 #6 (at its Rust position, 7): "
                                                "the Rust name is a misnomer, ClF₃ has two lone pairs",
    "organic::FunctionalGroup::hydrogen_bonding": "renamed functional_group_is_hydrogen_bond_donor by ADR 0002 #6 (with "
                                                  "the new functional_group_is_hydrogen_bond_acceptor beside it)",
    "inorganic::static IONIC_RADII": "the generated columns KM_IONIC_RADII_CHARGE / _COORDINATION / _RADIUS_PM and the "
                                     "symbol blob _KM_IONIC_RADII_SYMBOL_BYTES in src/inorganic_data.cyr "
                                     "(KM_IONIC_RADII_COUNT rows), read through handles (row + 1): lookup_ionic_radius and "
                                     "the ionic_radius_* accessors (port plan §6.4)",
    "nuclear::static ISOTOPES": "the generated columns KM_ISOTOPES_Z / _A / _HALF_LIFE_S, the symbol blob and the UTF-8 "
                                "decay-mode blob in src/nuclear_data.cyr (KM_ISOTOPES_COUNT rows), read through handles "
                                "(row + 1): lookup_isotope and the isotope_* accessors (port plan §6.4)",
}

# Later milestones (port plan §2.1): counted, not checked.
PENDING_MODULES = {
    "stoichiometry": "M4",
    "reaction_dynamics": "M5", "network": "M5", "fitting": "M5", "sensitivity": "M5",
}
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


def rust_items(module):
    """(label, kind, owner, name) for every public item before #[cfg(test)]."""
    src = (root / f"rust-old/src/{module}.rs").read_text()
    cut = src.find("#[cfg(test)]")
    src = src if cut < 0 else src[:cut]
    vis = r"pub(?:\(crate\))?"
    items = []
    depth = 0
    impl = None
    block = None  # ("struct" | "enum", name) while inside its body at depth 1
    for line in src.splitlines():
        s = line.strip()
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
total_found = total_items = 0
all_labels = set()
for module in MODULES:
    found = 0
    items = rust_items(module)
    for label, kind, owner, name in items:
        all_labels.add(label)
        if label in EXCEPTIONS:
            found += 1
            continue
        if has_counterpart(module, kind, owner, name):
            found += 1
        else:
            missing += 1
            print(f"MISSING  {label} ({kind})")
    total_found += found
    total_items += len(items)
    print(f"{module:<12} {found}/{len(items)} Rust public items have a Cyrius counterpart or a reason")
for label in EXCEPTIONS:
    if label not in all_labels:
        missing += 1
        print(f"STALE    EXCEPTIONS names {label}, which is no public item of a ported module")
# Every file under rust-old/src must be ported, pending, deferred or not applicable.
known = set(MODULES) | set(PENDING_MODULES) | set(DEFERRED_MODULES) | set(NOT_APPLICABLE_MODULES)
for path in sorted((root / "rust-old/src").rglob("*.rs")):
    mod = str(path.relative_to(root / "rust-old/src"))[:-3]
    if mod not in known:
        missing += 1
        print(f"UNKNOWN  rust-old/src/{mod}.rs is not ported, pending, deferred or marked not applicable")
for mod in sorted(set(PENDING_MODULES) | set(DEFERRED_MODULES) | set(NOT_APPLICABLE_MODULES)):
    if not (root / f"rust-old/src/{mod}.rs").exists():
        missing += 1
        print(f"STALE    {mod}: listed, but rust-old/src/{mod}.rs does not exist")
pend = {}
for mod, ms in PENDING_MODULES.items():
    pend.setdefault(ms, []).append(mod)
for ms in sorted(pend):
    n = sum(len(rust_items(m)) for m in pend[ms])
    print(f"pending {ms}: {len(pend[ms])} modules, {n} Rust public items ({', '.join(sorted(pend[ms]))})")
for mod, why in {**DEFERRED_MODULES, **NOT_APPLICABLE_MODULES}.items():
    print(f"{mod:<20} {len(rust_items(mod))} Rust public items, not ported: {why}")
print(f"surface parity: {total_found}/{total_items} in the ported modules; {missing} missing")
sys.exit(1 if missing else 0)
