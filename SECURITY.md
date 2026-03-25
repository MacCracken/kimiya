# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 0.1.x   | Yes       |

## Reporting

Report security issues to the repository maintainer. Do not open public issues for security vulnerabilities.

## Scope

Kimiya is a computation library with no network access, no file I/O, and no unsafe code. The primary attack surface is malformed input (NaN, infinity, negative concentrations) which is handled via validation.
