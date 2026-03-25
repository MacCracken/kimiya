//! Data fitting — parameter estimation from experimental measurements.
//!
//! Uses [`hisab`] least-squares and numerical methods for fitting
//! kinetic parameters, calibration curves, and thermochemical data.

use crate::element::GAS_CONSTANT;
use crate::error::{KimiyaError, Result};

/// Result of an Arrhenius fit: ln(k) = ln(A) - Ea/(R·T)
#[derive(Debug, Clone)]
pub struct ArrheniusFit {
    /// Pre-exponential factor A (s⁻¹).
    pub pre_exponential: f64,
    /// Activation energy Ea (J/mol).
    pub activation_energy_j: f64,
}

/// Fit Arrhenius parameters (A, Ea) from rate constant vs temperature data.
///
/// Linearizes: ln(k) = ln(A) - Ea/(R·T), then fits ln(k) vs 1/T.
///
/// Uses [`hisab::num::least_squares_poly`] with degree 1 for the linear fit.
///
/// - `temperatures_k`: temperature values (K)
/// - `rate_constants`: corresponding rate constants
///
/// # Errors
///
/// Returns error if data is insufficient or contains invalid values.
pub fn fit_arrhenius(temperatures_k: &[f64], rate_constants: &[f64]) -> Result<ArrheniusFit> {
    if temperatures_k.len() != rate_constants.len() {
        return Err(KimiyaError::InvalidInput(
            "temperature and rate arrays must have equal length".into(),
        ));
    }
    if temperatures_k.len() < 2 {
        return Err(KimiyaError::InvalidInput(
            "at least 2 data points required".into(),
        ));
    }

    // Linearize: x = 1/T, y = ln(k)
    let mut x = Vec::with_capacity(temperatures_k.len());
    let mut y = Vec::with_capacity(temperatures_k.len());
    for (&t, &k) in temperatures_k.iter().zip(rate_constants.iter()) {
        if t <= 0.0 {
            return Err(KimiyaError::InvalidTemperature(
                "temperatures must be positive".into(),
            ));
        }
        if k <= 0.0 {
            return Err(KimiyaError::InvalidInput(
                "rate constants must be positive".into(),
            ));
        }
        x.push(1.0 / t);
        y.push(k.ln());
    }

    // Linear fit: y = a₀ + a₁·x → ln(k) = ln(A) + (-Ea/R)·(1/T)
    let coeffs = hisab::num::least_squares_poly(&x, &y, 1)
        .map_err(|e| KimiyaError::ComputationError(format!("least squares failed: {e}")))?;

    let ln_a = coeffs[0]; // intercept = ln(A)
    let slope = coeffs[1]; // slope = -Ea/R

    Ok(ArrheniusFit {
        pre_exponential: ln_a.exp(),
        activation_energy_j: -slope * GAS_CONSTANT,
    })
}

/// Polynomial fit result with coefficients and evaluation.
#[derive(Debug, Clone)]
pub struct PolynomialFit {
    /// Coefficients \[a₀, a₁, a₂, ...\] where f(x) = a₀ + a₁·x + a₂·x² + ...
    pub coefficients: Vec<f64>,
}

impl PolynomialFit {
    /// Evaluate the fitted polynomial at a given x.
    #[must_use]
    pub fn evaluate(&self, x: f64) -> f64 {
        let mut result = 0.0;
        let mut x_pow = 1.0;
        for &c in &self.coefficients {
            result += c * x_pow;
            x_pow *= x;
        }
        result
    }

    /// Degree of the polynomial.
    #[must_use]
    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            0
        } else {
            self.coefficients.len() - 1
        }
    }
}

/// Fit a polynomial of given degree to (x, y) data.
///
/// Uses [`hisab::num::least_squares_poly`] (QR-based).
///
/// # Errors
///
/// Returns error if data is insufficient for the requested degree.
pub fn fit_polynomial(x: &[f64], y: &[f64], degree: usize) -> Result<PolynomialFit> {
    if x.len() != y.len() {
        return Err(KimiyaError::InvalidInput(
            "x and y arrays must have equal length".into(),
        ));
    }
    if x.len() <= degree {
        return Err(KimiyaError::InvalidInput(format!(
            "need at least {} data points for degree {degree} polynomial",
            degree + 1
        )));
    }
    let coefficients = hisab::num::least_squares_poly(x, y, degree)
        .map_err(|e| KimiyaError::ComputationError(format!("polynomial fit failed: {e}")))?;
    Ok(PolynomialFit { coefficients })
}

/// Fit Beer-Lambert calibration: A = ε·l·c → linear fit A vs c with known l.
///
/// Returns the molar absorptivity ε in L/(mol·cm).
///
/// - `concentrations`: c values (mol/L)
/// - `absorbances`: A values (dimensionless)
/// - `path_length`: l in cm
///
/// # Errors
///
/// Returns error if data is insufficient, path length invalid, or fit fails.
pub fn fit_beer_lambert(
    concentrations: &[f64],
    absorbances: &[f64],
    path_length: f64,
) -> Result<f64> {
    if path_length <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "path length must be positive".into(),
        ));
    }
    let fit = fit_polynomial(concentrations, absorbances, 1)?;
    // A = intercept + slope·c, slope = ε·l, so ε = slope/l
    let slope = fit.coefficients.get(1).copied().unwrap_or(0.0);
    Ok(slope / path_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Arrhenius fitting ────────────────────────────────────────────

    #[test]
    fn fit_arrhenius_known_data() {
        // Generate synthetic data: A=1e13, Ea=80000 J/mol
        let a_true = 1e13;
        let ea_true = 80_000.0;
        let temps = [300.0, 350.0, 400.0, 450.0, 500.0];
        let rates: Vec<f64> = temps
            .iter()
            .map(|&t| a_true * (-ea_true / (GAS_CONSTANT * t)).exp())
            .collect();

        let fit = fit_arrhenius(&temps, &rates).unwrap();
        assert!(
            (fit.activation_energy_j - ea_true).abs() / ea_true < 0.01,
            "Ea should be ~80kJ/mol, got {:.0}",
            fit.activation_energy_j
        );
        assert!(
            (fit.pre_exponential.ln() - a_true.ln()).abs() < 0.5,
            "ln(A) should be ~{:.1}, got {:.1}",
            a_true.ln(),
            fit.pre_exponential.ln()
        );
    }

    #[test]
    fn fit_arrhenius_insufficient_data() {
        assert!(fit_arrhenius(&[300.0], &[1.0]).is_err());
    }

    #[test]
    fn fit_arrhenius_negative_temp_is_error() {
        assert!(fit_arrhenius(&[300.0, -100.0], &[1.0, 2.0]).is_err());
    }

    #[test]
    fn fit_arrhenius_zero_rate_is_error() {
        assert!(fit_arrhenius(&[300.0, 400.0], &[1.0, 0.0]).is_err());
    }

    // ── Polynomial fitting ───────────────────────────────────────────

    #[test]
    fn fit_linear() {
        // y = 2 + 3x
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 + 3.0 * xi).collect();
        let fit = fit_polynomial(&x, &y, 1).unwrap();
        assert!((fit.coefficients[0] - 2.0).abs() < 0.01);
        assert!((fit.coefficients[1] - 3.0).abs() < 0.01);
    }

    #[test]
    fn fit_quadratic() {
        // y = 1 + 2x + 0.5x²
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 1.0 + 2.0 * xi + 0.5 * xi * xi).collect();
        let fit = fit_polynomial(&x, &y, 2).unwrap();
        assert!((fit.coefficients[0] - 1.0).abs() < 0.01);
        assert!((fit.coefficients[1] - 2.0).abs() < 0.01);
        assert!((fit.coefficients[2] - 0.5).abs() < 0.01);
    }

    #[test]
    fn polynomial_evaluate() {
        let fit = PolynomialFit {
            coefficients: vec![1.0, 2.0, 3.0], // 1 + 2x + 3x²
        };
        assert!((fit.evaluate(0.0) - 1.0).abs() < f64::EPSILON);
        assert!((fit.evaluate(1.0) - 6.0).abs() < f64::EPSILON);
        assert!((fit.evaluate(2.0) - 17.0).abs() < f64::EPSILON);
    }

    #[test]
    fn fit_polynomial_insufficient_data() {
        assert!(fit_polynomial(&[1.0, 2.0], &[1.0, 2.0], 2).is_err());
    }

    // ── Beer-Lambert fitting ─────────────────────────────────────────

    #[test]
    fn fit_beer_lambert_known() {
        // ε=100 L/(mol·cm), l=1cm → A = 100·c
        let conc = [0.01, 0.02, 0.03, 0.04, 0.05];
        let abs: Vec<f64> = conc.iter().map(|&c| 100.0 * c).collect();
        let epsilon = fit_beer_lambert(&conc, &abs, 1.0).unwrap();
        assert!(
            (epsilon - 100.0).abs() < 1.0,
            "ε should be ~100, got {epsilon}"
        );
    }

    #[test]
    fn fit_beer_lambert_different_path() {
        // ε=200, l=2cm → A = 400·c
        let conc = [0.001, 0.002, 0.003, 0.004];
        let abs: Vec<f64> = conc.iter().map(|&c| 200.0 * 2.0 * c).collect();
        let epsilon = fit_beer_lambert(&conc, &abs, 2.0).unwrap();
        assert!(
            (epsilon - 200.0).abs() < 5.0,
            "ε should be ~200, got {epsilon}"
        );
    }

    #[test]
    fn fit_beer_lambert_zero_path_is_error() {
        assert!(fit_beer_lambert(&[0.01], &[1.0], 0.0).is_err());
    }
}
