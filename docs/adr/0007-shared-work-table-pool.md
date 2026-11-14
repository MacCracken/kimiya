# 0007 — A shared pool of reusable work tables

**Status**: Accepted
**Date**: 2026-11-13

## Context

Some kimiya functions need a work table larger than their stack table. In 2.0.0 these were the ODE integrators, the
dense solve, the least squares, the rank, `compute_gradient`, `fft_spectrum`, `isotope_pattern`, `balance_equation`
and `empirical_formula`. Each took the table from lib/alloc.cyr's bump allocator, which never frees memory. So calling
any of them in a loop grew the heap by the table's size on every call (measured in the 2.0.1 review: 0.8–6.4 KB per
call for the integrators and fits, and 16·n bytes for `fft_spectrum`). Once the allocator's 2 GiB limit was reached,
every later call failed with a spurious `Err(KM_ERR_COMPUTATION_ERROR)`.

Rust freed each table when it went out of scope, so the leak is a port defect. Fixing it is a real choice for two
reasons:
- **Free-standing modules.** kimiya's modules have avoided module-level mutable state; ADR 0003's M2b amendment
  rejected callbacks that needed it.
- **Concurrency.** Cyrius threads are real, and the integrators call user callbacks, which can call back into kimiya.

## Decision

`src/validate.cyr` holds a pool of 8 grow-only buffers (`_km_scratch_take` / `_km_scratch_give`), and every heap work
table comes from it.
- **Claiming a slot.** A call claims a free slot by an atomic compare-and-swap, so no two threads and no re-entrant
  call ever share a buffer.
- **Growth.** A slot grows to at least twice its size when it is too small.
- **Release.** The call gives the slot back before it returns, on every path.
- **Overflow.** When all 8 slots are claimed, the call takes a fresh bump block, which is 2.0.0's behaviour.

The memory is not zeroed: each caller writes or zeroes what it reads. Only results the caller keeps (a `KmMolecule`,
a CSR handle) are still allocated per call. Nothing changes in the API or the results; every golden row is unchanged.

## Consequences

- **Positive:** after the first call at a given size, a call allocates nothing, so long-running consumers have
  bounded memory. The `cyrius_201_no_growth` tests check this for every site, and `tests/parity.tcyr` fails if any
  slot is still claimed after the whole corpus.
- **Negative:**
  - kimiya now holds process-lifetime state: up to 8 buffers, each at most twice the largest table it served. That
    memory is never returned.
  - A buffer outgrown by a larger call is left to the bump allocator.
  - Every function with a heap table must give its slot back on every return path. A missed release leaks the slot,
    not memory; the parity check catches it, and it caught one in `_km_matrix_rank` during 2.0.1.
- **Neutral:** more than 8 concurrent or nested calls fall back to 2.0.0's per-call growth.

## Alternatives considered

- **A caller-supplied work buffer.** This is an ADR 0002 API change on every affected function. It is still possible
  later, as an additive overload.
- **A per-call arena freed before return** (`arena_new_growable` / `arena_free`). This costs a map and an unmap per
  call, and the arena's header still comes from the bump allocator.
- **Documenting the growth.** This leaves long-running consumers to fail.
