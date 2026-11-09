# Architecture Decision Records

Decisions about kimiya — what we chose, the context, and the consequences we accept. Use these when a future reader would reasonably ask *"why did we do it this way?"*

## Conventions

- **Filename**: `NNNN-kebab-case-title.md`, zero-padded to four digits. Never renumber.
- **One decision per ADR.** If a decision supersedes a prior one, add a new ADR and set the old one's status to `Superseded by NNNN`.
- **Status lifecycle**: `Proposed` → `Accepted` → (optionally) `Superseded` or `Deprecated`.
- Use [`template.md`](template.md) as the starting point.

## ADR vs. architecture note vs. guide

| Kind | Lives in | Answers |
|---|---|---|
| ADR | `docs/adr/` | *Why did we choose X over Y?* |
| Architecture note | `docs/architecture/` | *What non-obvious constraint is true about the code?* |
| Guide | `docs/guides/` | *How do I do X?* |

## Index

| ADR | Title | Status |
|---|---|---|
| [0001](0001-cyrius-port-with-rust-oracle.md) | Port to Cyrius, with Rust 1.2.0 frozen as the parity oracle | Accepted |
| [0002](0002-port-api-divergences.md) | Accept the planned API changes in the Cyrius 2.0 API | Accepted |
| [0003](0003-transcribe-missing-hisab-routines.md) | Transcribe the hisab 1.4.0 routines that hisab 3.4.1 lacks (amended in M2b: the flame solver's Newton and bisection too) | Accepted |
| [0004](0004-accept-ganita-pow-minus-one-divergence.md) | Accept ganita's pow(−1, DBL_MAX) divergence | Accepted |
| [0005](0005-correctly-rounded-primitives-for-ill-conditioned-functions.md) | Correctly rounded acos, cos, ln and x^2.5 for the ill-conditioned gas functions | Accepted |
