"""Compile every generated constant through the Cyrius compiler and check each value's bits
against tests/golden. This proves src/constants.cyr reproduces the oracle's bits as the
compiler reads them (decimal-literal conversion included), not just as Python reads them.

Usage, from the repo root:
    python3 -I scripts/oracle/probe_tables.py .

It writes build/oracle-probe/ (gitignored), builds there with the project's cyrius.cyml
(which supplies the stdlib), and exits 1 on any mismatch, any compiler warning, or a
generated constant with no golden row. Adapted from tanmatra's probe_tables.py; extend
`checks` as the data modules arrive (element_data in M1, …).
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

golden = {}  # "element::AVOGADRO" -> "f:<bits>"
for fn, cls, ins, outs in rows("constants"):
    golden[ins[0][2:]] = outs[0]

checks = []  # (cyrius expression, expected token, label)
fails = 0
cyr = (root / "src/constants.cyr").read_text()
for m in re.finditer(r"^# Rust: kimiya::(\S+)\nconst (KM_[A-Z0-9_]+) = ", cyr, re.M):
    path, name = m.groups()
    if path not in golden:
        fails += 1
        print(f"FAIL {name}: kimiya::{path} has no golden row")
        continue
    checks.append((name, golden[path], f"const {path}"))
declared = len(re.findall(r"^const ", cyr, re.M))
if declared != len(checks) + fails:
    fails += 1
    print(f"FAIL src/constants.cyr declares {declared} consts, {len(checks)} carry a `# Rust:` path")

src = ['include "src/constants.cyr"', "", "fn main() {", "    var x: f64 = 0.0;"]
for k, (expr, want, label) in enumerate(checks):
    src.append(f"    x = {expr}; fmt_hex(x); println(\" {k}\");")
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

for k, (expr, want, label) in enumerate(checks):
    g = got.get(k)
    ok = g is not None and int(g, 16) == int(want[2:], 16)
    if not ok:
        fails += 1
        print(f"FAIL {label}: {expr} -> {g}, want {want[2:]}")
pending = len(golden) - len(checks)
print(f"checked {len(checks)} values through the Cyrius compiler ({pending} golden constants are "
      f"module data, generated at their milestone); failures: {fails}")
sys.exit(1 if fails else 0)
