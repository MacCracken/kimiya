# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 1.x     | Yes       |
| < 1.0   | No        |

## Reporting

Report security issues to the repository maintainer. Do not open public issues for security vulnerabilities.

## Scope

Kimiya is a computation library with no network access, no file I/O, and no unsafe code. The primary attack surface is malformed input (NaN, infinity, negative or zero quantities, extreme integer counts). Fallible functions validate inputs and return `KimiyaError`; library code does not `unwrap`, `expect` or `panic!` on caller input.

Supply chain: `cargo audit` and `cargo deny check` (licenses, advisories, sources) run in CI.
