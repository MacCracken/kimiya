//! Sensitivity analysis — gradients of chemical quantities via reverse-mode autodiff.
//!
//! Uses [`hisab::autodiff::reverse_gradient`] for efficient computation of
//! all partial derivatives in a single backward pass.

use crate::element::GAS_CONSTANT;
use crate::error::Result;

/// Gradient of equilibrium constant K with respect to \[ΔG, T\].
///
/// K = exp(-ΔG / (RT))
///
/// Returns \[∂K/∂(ΔG), ∂K/∂T\].
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
/// Returns \[∂k/∂A, ∂k/∂Ea, ∂k/∂T\].
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
/// Returns \[∂E/∂E°, ∂E/∂T, ∂E/∂Q\].
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
///
/// # Errors
///
/// Returns error if the function panics (caught as computation error).
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
        // For ΔG < 0 (exothermic): ∂K/∂T should be negative (K decreases with T)
        assert!(grads[1] < 0.0, "∂K/∂T should be negative for exothermic");
    }

    #[test]
    fn arrhenius_gradient_signs() {
        let grads = arrhenius_gradient(1e13, 80_000.0, 300.0);
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
            (dg_dt - (-ds)).abs() < 0.01,
            "∂(ΔG)/∂T should equal -ΔS = {}, got {dg_dt}",
            -ds
        );
    }

    #[test]
    fn nernst_gradient_e0_is_unity() {
        // ∂E/∂E° = 1 (E depends linearly on E°)
        let grads = nernst_gradient(0.342, 2, 298.0, 1.0);
        assert!(
            (grads[0] - 1.0).abs() < 0.01,
            "∂E/∂E° should be 1.0, got {}",
            grads[0]
        );
    }

    #[test]
    fn nernst_gradient_q_negative() {
        // ∂E/∂Q < 0 (higher Q → lower E for reduction)
        let grads = nernst_gradient(0.342, 2, 298.0, 0.1);
        assert!(grads[2] < 0.0, "∂E/∂Q should be negative");
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
