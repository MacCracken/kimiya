//! Sensitivity analysis — gradients of chemical quantities via reverse-mode autodiff.
//!
//! Uses [`hisab::autodiff::reverse_gradient`] for efficient computation of
//! all partial derivatives in a single backward pass.
//!
//! The fixed-formula gradient functions return plain values rather than
//! `Result` and do **not** validate their inputs. NaN/±∞ inputs, n = 0 and
//! Q ≤ 0 yield non-finite entries; T ≤ 0 yields non-finite or finite but
//! physically meaningless entries. Validate inputs first (e.g. with the
//! `Result`-returning functions in [`crate::kinetics`] and
//! [`crate::electrochemistry`]) or check `is_finite()` on the result.

use crate::element::GAS_CONSTANT;
use crate::error::Result;

/// Gradient of K = exp(−ΔG / (RT)) with ΔG and T as **independent** variables.
///
/// - `delta_g_j`: ΔG in J/mol
/// - `temperature_k`: T in K
///
/// Returns \[∂K/∂(ΔG) in mol/J, (∂K/∂T)_ΔG in K⁻¹\], i.e.
/// \[−K/(RT), K·ΔG/(RT²)\].
///
/// The temperature entry holds ΔG fixed. It is **not** the van 't Hoff
/// derivative of a real equilibrium constant, dK/dT = K·ΔH°/(RT²), because
/// ΔG° itself changes with temperature (∂ΔG°/∂T = −ΔS°). The two agree only
/// when ΔS° = 0. For the van 't Hoff temperature dependence, differentiate
/// ΔG° = ΔH° − TΔS° explicitly, e.g. with [`compute_gradient`].
///
/// Inputs are not validated: non-finite input yields non-finite entries, and
/// T ≤ 0 yields non-finite or meaningless entries.
#[must_use]
pub fn equilibrium_constant_gradient(delta_g_j: f64, temperature_k: f64) -> Vec<f64> {
    hisab::autodiff::reverse_gradient(
        |tape, vars| {
            let dg = vars[0];
            let t = vars[1];
            let r = tape.constant(GAS_CONSTANT);
            let neg_dg = tape.neg(dg);
            let rt = tape.mul(r, t);
            let exponent = tape.div(neg_dg, rt);
            tape.exp(exponent)
        },
        &[delta_g_j, temperature_k],
    )
}

/// Gradient of Arrhenius rate constant with respect to \[A, Ea, T\].
///
/// k = A × exp(-Ea / (RT))
///
/// - `pre_exponential`: A (units of k)
/// - `activation_energy_j`: Ea in J/mol
/// - `temperature_k`: T in K
///
/// Returns \[∂k/∂A (dimensionless), ∂k/∂Ea (k-units·mol/J), ∂k/∂T
/// (k-units/K)\] = \[k/A, −k/(RT), k·Ea/(RT²)\].
///
/// Inputs are not validated: non-finite input yields non-finite entries, and
/// T ≤ 0 yields non-finite or meaningless (e.g. T < 0 makes k grow with Ea)
/// entries.
#[must_use]
pub fn arrhenius_gradient(
    pre_exponential: f64,
    activation_energy_j: f64,
    temperature_k: f64,
) -> Vec<f64> {
    hisab::autodiff::reverse_gradient(
        |tape, vars| {
            let a = vars[0];
            let ea = vars[1];
            let t = vars[2];
            let r = tape.constant(GAS_CONSTANT);
            let neg_ea = tape.neg(ea);
            let rt = tape.mul(r, t);
            let exponent = tape.div(neg_ea, rt);
            let boltzmann = tape.exp(exponent);
            tape.mul(a, boltzmann)
        },
        &[pre_exponential, activation_energy_j, temperature_k],
    )
}

/// Gradient of Nernst potential with respect to \[E°, T, Q\].
///
/// E = E° - (RT / nF) × ln(Q)
///
/// - `standard_potential`: E° in V
/// - `n_electrons`: electrons transferred
/// - `temperature_k`: T in K
/// - `reaction_quotient`: Q (dimensionless)
///
/// Returns \[∂E/∂E° (dimensionless), ∂E/∂T (V/K), ∂E/∂Q (V)\] =
/// \[1, −R·ln(Q)/(nF), −RT/(nFQ)\].
///
/// Inputs are not validated: n = 0 (division by zero), Q ≤ 0 (log of a
/// non-positive number) and non-finite input yield non-finite entries; T ≤ 0
/// yields meaningless entries.
#[must_use]
pub fn nernst_gradient(
    standard_potential: f64,
    n_electrons: u8,
    temperature_k: f64,
    reaction_quotient: f64,
) -> Vec<f64> {
    let nf = n_electrons as f64 * crate::electrochemistry::FARADAY;
    hisab::autodiff::reverse_gradient(
        |tape, vars| {
            let e0 = vars[0];
            let t = vars[1];
            let q = vars[2];
            let r = tape.constant(GAS_CONSTANT);
            let nf_const = tape.constant(nf);
            let rt = tape.mul(r, t);
            let ln_q = tape.ln(q);
            let rt_over_nf = tape.div(rt, nf_const);
            let correction = tape.mul(rt_over_nf, ln_q);
            tape.sub(e0, correction)
        },
        &[standard_potential, temperature_k, reaction_quotient],
    )
}

/// Sensitivity of Gibbs free energy to temperature: ∂(ΔG)/∂T = -ΔS
///
/// This is a direct thermodynamic identity, but computed via autodiff
/// to validate and demonstrate the framework.
///
/// - `enthalpy_j`: ΔH in J/mol
/// - `temperature_k`: T in K
/// - `entropy_j_per_k`: ΔS in J/(mol·K)
///
/// Returns ∂(ΔG)/∂T in J/(mol·K). Non-finite input propagates to the result.
#[must_use]
pub fn gibbs_temperature_sensitivity(
    enthalpy_j: f64,
    temperature_k: f64,
    entropy_j_per_k: f64,
) -> f64 {
    let grads = hisab::autodiff::reverse_gradient(
        |tape, vars| {
            let h = vars[0];
            let t = vars[1];
            let s = vars[2];
            let ts = tape.mul(t, s);
            tape.sub(h, ts) // ΔG = ΔH - TΔS
        },
        &[enthalpy_j, temperature_k, entropy_j_per_k],
    );
    grads[1] // ∂(ΔG)/∂T
}

/// Generic sensitivity: compute all partial derivatives of a scalar function.
///
/// Wraps [`hisab::autodiff::reverse_gradient`] for custom chemistry functions.
/// `f` receives a fresh tape and one variable per entry of `parameters`, and
/// returns the output variable. Returns one partial derivative per parameter,
/// in order.
///
/// # Errors
///
/// Currently never returns `Err`; the `Result` is kept for API stability.
/// Non-finite values produced by `f` appear in the gradient as-is.
///
/// # Panics
///
/// Panics are not caught: a panic inside `f` propagates to the caller, and
/// hisab panics if `f` returns a `Var` that was not created on the tape it
/// was given.
pub fn compute_gradient(
    f: impl Fn(&mut hisab::autodiff::Tape, &[hisab::autodiff::Var]) -> hisab::autodiff::Var,
    parameters: &[f64],
) -> Result<Vec<f64>> {
    Ok(hisab::autodiff::reverse_gradient(f, parameters))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equilibrium_k_gradient_signs() {
        // K = exp(-ΔG/(RT))
        // ∂K/∂(ΔG) should be negative (higher ΔG → lower K)
        // ∂K/∂T depends on sign of ΔG
        let grads = equilibrium_constant_gradient(-5000.0, 298.0);
        assert!(grads[0] < 0.0, "∂K/∂(ΔG) should be negative");
        // At fixed ΔG < 0 (exergonic — says nothing about ΔH): (∂K/∂T)_ΔG =
        // K·ΔG/(RT²) < 0, because −ΔG/(RT) shrinks as T grows.
        assert!(grads[1] < 0.0, "(∂K/∂T)_ΔG should be negative for ΔG < 0");
    }

    #[test]
    fn equilibrium_k_gradient_is_constant_delta_g_derivative() {
        // Pins what is returned: [−K/(RT), K·ΔG/(RT²)] — not van 't Hoff.
        let (dg, t) = (-5000.0, 298.0);
        let rt = GAS_CONSTANT * t;
        let k = (-dg / rt).exp();
        let grads = equilibrium_constant_gradient(dg, t);
        assert!((grads[0] - (-k / rt)).abs() <= 1e-12 * (k / rt));
        let expected = k * dg / (rt * t);
        assert!((grads[1] - expected).abs() <= 1e-12 * expected.abs());
    }

    #[test]
    fn invalid_input_gives_non_finite_gradients() {
        // Documented behaviour of the unvalidated, non-Result gradient fns.
        assert!(!equilibrium_constant_gradient(-5000.0, 0.0)[1].is_finite());
        assert!(!equilibrium_constant_gradient(f64::NAN, 298.0)[0].is_finite());
        assert!(!arrhenius_gradient(1e13, 80_000.0, f64::NAN)[2].is_finite());
        assert!(!nernst_gradient(0.342, 0, 298.0, 0.1)[2].is_finite());
        assert!(!nernst_gradient(0.342, 2, 298.0, 0.0)[2].is_finite());
        assert!(!nernst_gradient(0.342, 2, 298.0, -1.0)[1].is_finite());
    }

    #[test]
    fn arrhenius_gradient_signs() {
        let grads = arrhenius_gradient(1e13, 80_000.0, 300.0);
        let rt = GAS_CONSTANT * 300.0;
        let k = 1e13 * (-80_000.0 / rt).exp();
        assert!((grads[0] - k / 1e13).abs() <= 1e-12 * k / 1e13);
        assert!((grads[1] + k / rt).abs() <= 1e-12 * k / rt);
        let dk_dt = k * 80_000.0 / (rt * 300.0);
        assert!((grads[2] - dk_dt).abs() <= 1e-12 * dk_dt);
        // ∂k/∂A > 0 (more pre-exponential → faster)
        assert!(grads[0] > 0.0, "∂k/∂A should be positive");
        // ∂k/∂Ea < 0 (higher barrier → slower)
        assert!(grads[1] < 0.0, "∂k/∂Ea should be negative");
        // ∂k/∂T > 0 (higher T → faster)
        assert!(grads[2] > 0.0, "∂k/∂T should be positive");
    }

    #[test]
    fn gibbs_temperature_sensitivity_equals_neg_entropy() {
        // ∂(ΔG)/∂T = -ΔS (thermodynamic identity)
        let dh = -100_000.0;
        let t = 298.0;
        let ds = -200.0;
        let dg_dt = gibbs_temperature_sensitivity(dh, t, ds);
        assert!(
            (dg_dt - (-ds)).abs() < 1e-12,
            "∂(ΔG)/∂T should equal -ΔS = {}, got {dg_dt}",
            -ds
        );
    }

    #[test]
    fn nernst_gradient_e0_is_unity() {
        // ∂E/∂E° = 1 (E depends linearly on E°)
        let grads = nernst_gradient(0.342, 2, 298.0, 1.0);
        assert!(
            (grads[0] - 1.0).abs() < 1e-15,
            "∂E/∂E° should be 1.0, got {}",
            grads[0]
        );
    }

    #[test]
    fn nernst_gradient_q_negative() {
        // ∂E/∂Q < 0 (higher Q → lower E for reduction)
        let grads = nernst_gradient(0.342, 2, 298.0, 0.1);
        assert!(grads[2] < 0.0, "∂E/∂Q should be negative");
        let nf = 2.0 * crate::electrochemistry::FARADAY;
        let de_dt = -GAS_CONSTANT * 0.1_f64.ln() / nf;
        assert!((grads[1] - de_dt).abs() <= 1e-12 * de_dt.abs());
        let de_dq = -GAS_CONSTANT * 298.0 / (nf * 0.1);
        assert!((grads[2] - de_dq).abs() <= 1e-12 * de_dq.abs());
    }

    #[test]
    fn custom_gradient() {
        // f(x, y) = x² + 3y → ∂f/∂x = 2x, ∂f/∂y = 3
        let grads = compute_gradient(
            |tape, vars| {
                let x2 = tape.mul(vars[0], vars[0]);
                let three = tape.constant(3.0);
                let three_y = tape.mul(three, vars[1]);
                tape.add(x2, three_y)
            },
            &[2.0, 5.0],
        )
        .unwrap();
        assert!((grads[0] - 4.0).abs() < 1e-10, "∂f/∂x should be 4.0");
        assert!((grads[1] - 3.0).abs() < 1e-10, "∂f/∂y should be 3.0");
    }
}
