# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 2.x (Cyrius) | Yes |
| 1.x (Rust crate) | Frozen at 1.2.0, the parity oracle in `rust-old/` |
| < 1.0   | No        |

## Reporting

Report security issues to the repository maintainer. Do not open public issues for security vulnerabilities.

## Scope

Kimiya is a computation library with no network access and no file I/O. The primary attack surface is malformed input
(NaN, infinity, negative or zero quantities, extreme integer counts, undersized caller buffers). Fallible functions
validate their inputs and return a value-form `Result` or an `err_out` code. Library code never aborts the process:
`src/` uses no `+?`, `-?`, `*?`, `assert_fatal` or `result_unwrap`. The whole test suite also runs under
`cyrius test --bounds`.

Supply chain: dependencies are pinned in `cyrius.lock`; `cyrius deps --locked` and `cyrius audit` run in CI.
