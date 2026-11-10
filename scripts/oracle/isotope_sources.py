"""Re-verify scripts/oracle/isotopes.tsv (ADR 0002 #23) against fresh downloads of its sources. Usage, from the repo root:

    python3 -I scripts/oracle/isotope_sources.py . DIR

DIR holds each source under its URL's basename: mass_1.mas20.txt and TICE-2009.pdf (required), hydrogen.htm and
sulfur.htm (the `interval` records' pages; checked when present). The downloads are data: this script only reads them,
and needs nothing beyond the Python standard library (run it with -I, from outside DIR).

1. The two publications' sha256s are the ones the data file's header records (a CIAAW web page may change: its sha256
   is printed against the recorded one, and its interval text is what is checked).
2. AME2020: every isotope record's `mas20` field is a whole line of mass_1.mas20.txt, the file's only line for that Z and
   A, with that element symbol; and `mass_u` is that line's atomic mass, read by the file's own Fortran format (the i3 in
   u at columns 107-109, the f13.6 in micro-u at 111-123), its digits unchanged.
3. TICE 2009: the PDF's content streams inflated (zlib), each TJ / Tj text run split into cells at the column gaps (a
   kerning step of 1000 or more thousandths of an em; smaller steps join, as in "0.996 )150(36(20)"), and every record's
   `tice` value is such a cell, the 17 in record order (Table 1's order).
4. CIAAW: each `interval` record's page says "A_r(E) = <as printed> since 2009".
Exit 1 on any failure.
"""
import hashlib, html, re, sys, zlib
from pathlib import Path

root, dl = Path(sys.argv[1]), Path(sys.argv[2])
data = (root / "scripts/oracle/isotopes.tsv").read_text(encoding="utf-8")
FAIL = [0]

def check(ok, msg):
    if not ok:
        FAIL[0] += 1
        print(f"FAIL: {msg}")

# The header's (URL, sha256) pairs: a `#   https://…` line, then (on the next comment line) `sha256 <hex>`, or the
# sha256 named in the same parenthesis as the URL.
shas = {}
head = [l[1:].strip() for l in data.splitlines() if l.startswith("#")]
for k, l in enumerate(head):
    for url in re.findall(r"https://\S+?(?=[ )]|$)", l):
        m = re.search(rf"{re.escape(url)} \(sha256 ([0-9a-f]{{64}})\)", l)
        if m is None and k + 1 < len(head):
            m = re.match(r"sha256 ([0-9a-f]{64})$", head[k + 1])
        if m is not None:
            shas[url.rsplit("/", 1)[1]] = m.group(1)
isos = [l.split("\t") for l in data.splitlines() if l.startswith("isotope\t")]
ivals = [l.split("\t") for l in data.splitlines() if l.startswith("interval\t")]
check(len(isos) == 17 and all(len(r) == 8 for r in isos), f"{len(isos)} isotope records (17 of 8 fields expected)")
check(all(len(r) == 6 for r in ivals), "an interval record has not 6 fields")

# 1. the hashes.
for name in ("mass_1.mas20.txt", "TICE-2009.pdf"):
    p = dl / name
    check(p.is_file(), f"{p}: missing")
    if p.is_file():
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        check(got == shas.get(name), f"{name}: sha256 {got}, the data file records {shas.get(name)}")
        print(f"{name}: sha256 {got} ({'as recorded' if got == shas.get(name) else 'DIFFERS'})")

# 2. AME2020.
mas = (dl / "mass_1.mas20.txt").read_text(encoding="ascii").split("\n")
by_za = {}
for l in mas:
    try:
        za = (int(l[9:14]), int(l[14:19]))
    except ValueError:
        continue
    if len(l) >= 123 and l[19] == " ":
        by_za.setdefault(za, []).append(l)
for kind, z, a, sym, mass, ab, tice, line in isos:
    z, a = int(z), int(a)
    check(by_za.get((z, a)) == [line], f"{sym}-{a}: the recorded mas20 line is not the file's only (Z, A) = ({z}, {a}) line")
    check(line[20:23].strip() == sym, f"{sym}-{a}: the line's element is {line[20:23].strip()!r}")
    ui, um = line[106:109].strip(), line[110:123].strip()
    ip, _, fp = um.partition(".")
    check(ui.isdigit() and len(ip) == 6 and ip.isdigit() and fp.isdigit() and mass == f"{int(ui)}.{ip}{fp}",
          f"{sym}-{a}: mass_u {mass} is not the line's atomic mass {ui} {um} micro-u")
print(f"mass_1.mas20.txt: the {len(isos)} lines are verbatim and give each mass_u")

# 3. TICE 2009.
pdf = (dl / "TICE-2009.pdf").read_bytes()
cells = []
for m in re.finditer(rb"stream(\r\n|\n|\r)", pdf):
    try:
        body = zlib.decompressobj().decompress(pdf[m.end():pdf.find(b"endstream", m.end())])
    except zlib.error:
        continue
    for arr, single in re.findall(rb"\[((?:\((?:\\.|[^\\)])*\)|-?[0-9.]+|\s)*)\]\s*TJ|\(((?:\\.|[^\\)])*)\)\s*Tj", body):
        parts = re.findall(rb"\(((?:\\.|[^\\)])*)\)|(-?[0-9.]+)", arr) if arr else [(single, b"")]
        cur = b""
        for s, kern in parts:
            if kern and float(kern) <= -1000:
                cells.append(cur.strip())
                cur = b""
            elif not kern:
                cur += re.sub(rb"\\(.)", rb"\1", s)
        cells.append(cur.strip())
at = 0
for kind, z, a, sym, mass, ab, tice, line in isos:
    want = tice.encode()
    k = next((j for j in range(at, len(cells)) if cells[j] == want), None)
    check(k is not None, f"{sym}-{a}: {tice!r} is not a Table 1 cell after the previous record's")
    at = (k if k is not None else at) + 1
    check(ab == tice.split("(")[0].replace(" ", ""), f"{sym}-{a}: abundance {ab} is not {tice!r}'s digits")
print(f"TICE-2009.pdf: the {len(isos)} representative values are Table 1 cells, in record order")

# 4. CIAAW's intervals.
for kind, sym, lo, hi, printed, page in ivals:
    p = dl / page
    if not p.is_file():
        print(f"{page}: not downloaded; the {sym} interval is not re-checked")
        continue
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    text = html.unescape(re.sub(r"<[^>]+>", "", p.read_text(encoding="utf-8", errors="replace")))
    check(f"Ar({sym}) = {printed} since 2009" in text, f"{page}: no \"Ar({sym}) = {printed} since 2009\"")
    m = re.fullmatch(r"\[([0-9. ]+), ([0-9. ]+)\]", printed)
    check(m is not None and (lo, hi) == (m.group(1).replace(" ", ""), m.group(2).replace(" ", "")),
          f"{sym}: [{lo}, {hi}] is not {printed!r}")
    print(f"{page}: Ar({sym}) = {printed}; sha256 {got} ({'as recorded' if got == shas.get(page) else 'changed since retrieval'})")

print(f"isotope_sources: {FAIL[0]} failure(s)")
sys.exit(1 if FAIL[0] else 0)
