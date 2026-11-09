# Benchmarks — Rust 1.2.0 vs the Cyrius port

Each row times one benchmark body: the same work in Rust (`rust-old/`, release build) and in
Cyrius (`tests/kimiya.bcyr`, `cyrius bench`). Port milestones append rows as they land. Numbers
are host-specific, so compare within a column set, never across hosts.

## Method

- **Rust, criterion:** the *latest* column of [`rust-old/benchmarks.md`](../rust-old/benchmarks.md)
  (criterion 0.8 medians, v1.2 post-audit, rustc 1.99.1 release, 2026-11-08), from
  `cargo bench --manifest-path rust-old/Cargo.toml --bench benchmarks -- <name>`.
- **Cyrius:** `cyrius bench`, mean of 20 batches. Each batch is a tight loop of `per_batch` calls
  (as criterion's `b.iter` loop), inputs read from globals (Rust's `black_box`) and results written
  to a sink. lib/bench.cyr measures its clock-read floor and subtracts it from every sample
  (1.31 µs per read on this host); the per-op figure is lib/bench.cyr's ps mean.
- **Host:** x86_64 Linux, HPET clocksource. Toolchains: rustc 1.99.1 (release); Cyrius 6.7.18
  (ganita 1.2.15).

## M1 Foundation (2026-11-09)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `element/lookup_by_symbol` | `lookup_by_symbol("Fe")`: a linear scan to row 26 | 44.73 ns | 67.0 ns | 1.5× |
| `element/lookup_by_number` | `lookup_by_number(26)`: range check, one table read | 0.84 ns | 4.54 ns | 5.4× |
| `molecule/molecular_weight_glucose` | `molecular_weight()` of C₆H₁₂O₆, built outside the loop: 3 lookups, 3 multiply-adds | 3.63 ns | 26.7 ns | 7.4× |

All three bodies are bit-identical to Rust (golden rows: `element` 975, `molecule` 2,472, 0 ulp).

**Two optimizations, measured.** The first version of each body measured 92.4 ns, 7.03 ns and
121 ns:
- `molecular_weight` read each atom through `atom_element_number(k)` / `atom_count(k)`, each a
  bounds check, a `vec_get` with its own two checks, and calls. It now walks the vec's buffer
  directly (lib/vec.cyr documents the `{ data, len, cap }` layout) and reads the mass column at
  the checked row: 121 → 57 ns with one `vec_get` per word, → 33 ns walking the buffer.
- `element_lookup_by_number` called `_km_element_row` for its range check; the check is inline
  now: 7.03 → 4.54 ns, and `molecular_weight`, which makes three of them per call, 33 → 26.7 ns.
- `element_lookup_by_symbol` compares the first byte before calling `memeq`, with the probe's
  first byte read once: 92.4 → 67.0 ns.

**Reading the gap.** Rust inlines `lookup_by_number` into a compare and a load, and unrolls the
glucose sum; Cyrius compiles each call as a call (the lookup inside `molecular_weight`, the
`Result` pair return). The symbol scan is the closest, because both sides do the same linear
string search and Rust's per-row cost is a real comparison too. Parity, not speed, is M1's bar.
