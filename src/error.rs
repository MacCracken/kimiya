use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum KimiyaError {
    #[error("invalid element: {0}")]
    InvalidElement(String),
    #[error("invalid reaction: {0}")]
    InvalidReaction(String),
    #[error("invalid concentration: {0}")]
    InvalidConcentration(String),
    #[error("invalid temperature: {0}")]
    InvalidTemperature(String),
    #[error("invalid input: {0}")]
    InvalidInput(String),
    #[error("computation error: {0}")]
    ComputationError(String),
}

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
