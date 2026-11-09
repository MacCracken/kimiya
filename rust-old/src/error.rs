//! Error type shared by all kimiya modules.
//!
//! Fallible functions return [`Result<T>`], an alias for
//! `std::result::Result<T, KimiyaError>`. Each variant carries a
//! human-readable message describing the offending input.

use thiserror::Error;

/// Errors returned by kimiya functions.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum KimiyaError {
    /// Unknown element: atomic number outside 1–118, or unrecognized symbol.
    #[error("invalid element: {0}")]
    InvalidElement(String),
    /// Reaction or species that cannot be processed (e.g. an equation that
    /// cannot be balanced, or an unknown redox couple or species).
    #[error("invalid reaction: {0}")]
    InvalidReaction(String),
    /// Concentration out of range (e.g. negative or non-finite).
    #[error("invalid concentration: {0}")]
    InvalidConcentration(String),
    /// Temperature out of range (e.g. not positive in kelvin, or non-finite).
    #[error("invalid temperature: {0}")]
    InvalidTemperature(String),
    /// Any other argument outside its valid domain (zero, negative, NaN,
    /// infinite, or inconsistent with other arguments).
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// Numerical failure during computation (e.g. a solver failed to
    /// converge or produced a degenerate result).
    #[error("computation error: {0}")]
    ComputationError(String),
}

/// Result alias used throughout kimiya: `std::result::Result<T, KimiyaError>`.
pub type Result<T> = std::result::Result<T, KimiyaError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_display() {
        let e = KimiyaError::InvalidElement("Unobtainium".into());
        assert!(e.to_string().contains("Unobtainium"));
    }

    #[test]
    fn result_type() {
        let ok: Result<f64> = Ok(1.0);
        assert!(ok.is_ok());
    }
}
