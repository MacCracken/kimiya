//! Reaction dynamics — time-evolution of chemical systems via ODE integration.
//!
//! Uses [`hisab`] numerical methods to simulate concentration trajectories:
//! classical RK4 (fixed step, explicit) for the built-in rate laws,
//! Dormand–Prince 4(5) (adaptive, explicit) for custom systems, backward
//! Euler (implicit, Newton iteration on [`hisab::num::gaussian_elimination`])
//! for stiff systems, and Euler–Maruyama for stochastic systems.
//!
//! Time runs from 0 to `t_end` in the time unit of the rate constants (s for
//! rate constants in s⁻¹); concentrations are in mol/L. Fixed-step methods
//! accept at most 10 000 000 steps.

use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};

/// A point in a kinetics trajectory: (time, concentrations).
pub type TrajectoryPoint = (f64, Vec<f64>);

/// Largest accepted step count. hisab pre-allocates `steps + 1` trajectory
/// points, so unbounded counts would overflow or exhaust memory.
const MAX_STEPS: usize = 10_000_000;

/// Validate the end time and step count shared by the fixed-step solvers.
fn validate_time_and_steps(t_end: f64, steps: usize) -> Result<()> {
    if !is_positive(t_end) {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    if steps > MAX_STEPS {
        return Err(KimiyaError::InvalidInput(format!(
            "steps must be at most {MAX_STEPS}"
        )));
    }
    Ok(())
}

/// Reject a negative or non-finite rate constant.
fn validate_rate_constant(k: f64, name: &str) -> Result<()> {
    if is_non_negative(k) {
        Ok(())
    } else {
        Err(KimiyaError::InvalidInput(format!(
            "{name} must be non-negative and finite"
        )))
    }
}

/// Reject a non-finite entry in a user-supplied initial state.
fn validate_state(y0: &[f64]) -> Result<()> {
    if y0.iter().all(|&v| is_finite(v)) {
        Ok(())
    } else {
        Err(KimiyaError::InvalidInput(
            "initial state must be finite".into(),
        ))
    }
}

/// Map a hisab trajectory to kimiya's error type and reject a trajectory whose
/// final state is non-finite (e.g. an explicit method that went unstable).
fn finish_trajectory(
    result: std::result::Result<Vec<TrajectoryPoint>, hisab::HisabError>,
    solver: &str,
) -> Result<Vec<TrajectoryPoint>> {
    let trajectory =
        result.map_err(|e| KimiyaError::ComputationError(format!("{solver} failed: {e}")))?;
    let final_ok = trajectory
        .last()
        .is_some_and(|(_, y)| y.iter().all(|&v| is_finite(v)));
    if final_ok {
        Ok(trajectory)
    } else {
        Err(KimiyaError::ComputationError(format!(
            "{solver} produced a non-finite state (unstable: use more steps or simulate_stiff)"
        )))
    }
}

/// Simulate first-order decay: d\[A\]/dt = -k\[A\]
///
/// - `initial_concentration`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (s⁻¹)
/// - `t_end`: end time (s)
/// - `steps`: number of RK4 steps (1 ..= 10 000 000)
///
/// Returns `steps + 1` points (t, \[\[A\]\]) from t = 0 to `t_end`. RK4 is
/// explicit: it is only stable for k·t_end/steps ≲ 2.78.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if \[A\]₀ is negative or not
/// finite, [`KimiyaError::InvalidInput`] if k is negative or not finite,
/// `t_end` is not positive and finite, or `steps` is 0 or above 10 000 000,
/// and [`KimiyaError::ComputationError`] if the solver fails or the final
/// state is not finite.
pub fn simulate_first_order(
    initial_concentration: f64,
    rate_constant: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_non_negative(initial_concentration) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    validate_rate_constant(rate_constant, "rate constant")?;
    validate_time_and_steps(t_end, steps)?;
    let k = rate_constant;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k * y[0];
    };
    finish_trajectory(
        hisab::num::rk4_trajectory(f, 0.0, &[initial_concentration], t_end, steps),
        "ODE solver",
    )
}

/// Simulate second-order decay: d\[A\]/dt = -k\[A\]²
///
/// - `initial_concentration`: \[A\]₀ (mol/L)
/// - `rate_constant`: k (L·mol⁻¹·s⁻¹)
/// - `t_end`: end time (s)
/// - `steps`: number of RK4 steps (1 ..= 10 000 000)
///
/// Returns `steps + 1` points (t, \[\[A\]\]) from t = 0 to `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if \[A\]₀ is negative or not
/// finite, [`KimiyaError::InvalidInput`] if k is negative or not finite,
/// `t_end` is not positive and finite, or `steps` is 0 or above 10 000 000,
/// and [`KimiyaError::ComputationError`] if the solver fails or the final
/// state is not finite.
pub fn simulate_second_order(
    initial_concentration: f64,
    rate_constant: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_non_negative(initial_concentration) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    validate_rate_constant(rate_constant, "rate constant")?;
    validate_time_and_steps(t_end, steps)?;
    let k = rate_constant;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k * y[0] * y[0];
    };
    finish_trajectory(
        hisab::num::rk4_trajectory(f, 0.0, &[initial_concentration], t_end, steps),
        "ODE solver",
    )
}

/// Simulate consecutive reactions: A →(k₁) B →(k₂) C
///
/// State vector: \[A\], \[B\], \[C\] (mol/L), starting from (\[A\]₀, 0, 0).
///
/// - `k1`, `k2`: first-order rate constants (s⁻¹)
/// - `t_end`: end time (s)
/// - `steps`: number of RK4 steps (1 ..= 10 000 000)
///
/// Returns `steps + 1` points of all three species from t = 0 to `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if \[A\]₀ is negative or not
/// finite, [`KimiyaError::InvalidInput`] if a rate constant is negative or not
/// finite, `t_end` is not positive and finite, or `steps` is 0 or above
/// 10 000 000, and [`KimiyaError::ComputationError`] if the solver fails or
/// the final state is not finite.
pub fn simulate_consecutive(
    initial_a: f64,
    k1: f64,
    k2: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_non_negative(initial_a) {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    validate_rate_constant(k1, "k1")?;
    validate_rate_constant(k2, "k2")?;
    validate_time_and_steps(t_end, steps)?;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k1 * y[0]; // d[A]/dt = -k1[A]
        dy[1] = k1 * y[0] - k2 * y[1]; // d[B]/dt = k1[A] - k2[B]
        dy[2] = k2 * y[1]; // d[C]/dt = k2[B]
    };
    finish_trajectory(
        hisab::num::rk4_trajectory(f, 0.0, &[initial_a, 0.0, 0.0], t_end, steps),
        "ODE solver",
    )
}

/// Simulate reversible reaction: A ⇌ B with forward rate k_f and reverse rate k_r
///
/// State vector: \[A\], \[B\] (mol/L).
///
/// - `k_forward`, `k_reverse`: first-order rate constants (s⁻¹)
/// - `t_end`: end time (s)
/// - `steps`: number of RK4 steps (1 ..= 10 000 000)
///
/// Returns `steps + 1` points from t = 0 to `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if an initial concentration
/// is negative or not finite, [`KimiyaError::InvalidInput`] if a rate constant
/// is negative or not finite, `t_end` is not positive and finite, or `steps`
/// is 0 or above 10 000 000, and [`KimiyaError::ComputationError`] if the
/// solver fails or the final state is not finite.
pub fn simulate_reversible(
    initial_a: f64,
    initial_b: f64,
    k_forward: f64,
    k_reverse: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_non_negative(initial_a) || !is_non_negative(initial_b) {
        return Err(KimiyaError::InvalidConcentration(
            "concentrations must be non-negative".into(),
        ));
    }
    validate_rate_constant(k_forward, "forward rate constant")?;
    validate_rate_constant(k_reverse, "reverse rate constant")?;
    validate_time_and_steps(t_end, steps)?;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k_forward * y[0] + k_reverse * y[1]; // d[A]/dt
        dy[1] = k_forward * y[0] - k_reverse * y[1]; // d[B]/dt
    };
    finish_trajectory(
        hisab::num::rk4_trajectory(f, 0.0, &[initial_a, initial_b], t_end, steps),
        "ODE solver",
    )
}

/// Simulate Michaelis-Menten enzyme kinetics with substrate depletion.
///
/// Full model: E + S ⇌(k₁/k₋₁) ES →(k₂) E + P
///
/// Simplified to: d\[S\]/dt = -Vmax·\[S\]/(Km + \[S\]), d\[P\]/dt = Vmax·\[S\]/(Km + \[S\])
///
/// State vector: \[S\], \[P\] (mol/L), starting from (\[S\]₀, 0).
///
/// - `v_max`: maximum rate (mol·L⁻¹·s⁻¹)
/// - `km`: Michaelis constant (mol/L)
/// - `t_end`: end time (s)
/// - `steps`: number of RK4 steps (1 ..= 10 000 000)
///
/// Returns `steps + 1` points from t = 0 to `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidConcentration`] if \[S\]₀ is negative or not
/// finite, [`KimiyaError::InvalidInput`] if Vmax or Km is not positive and
/// finite, `t_end` is not positive and finite, or `steps` is 0 or above
/// 10 000 000, and [`KimiyaError::ComputationError`] if the solver fails or
/// the final state is not finite.
pub fn simulate_michaelis_menten(
    initial_substrate: f64,
    v_max: f64,
    km: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_non_negative(initial_substrate) {
        return Err(KimiyaError::InvalidConcentration(
            "substrate concentration must be non-negative".into(),
        ));
    }
    if !is_positive(v_max) {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if !is_positive(km) {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    validate_time_and_steps(t_end, steps)?;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        let rate = v_max * y[0] / (km + y[0]);
        dy[0] = -rate; // d[S]/dt
        dy[1] = rate; // d[P]/dt
    };
    finish_trajectory(
        hisab::num::rk4_trajectory(f, 0.0, &[initial_substrate, 0.0], t_end, steps),
        "ODE solver",
    )
}

/// Simulate a custom reaction system with adaptive step size.
///
/// Uses Dormand–Prince 4(5), an explicit adaptive Runge–Kutta method suited to
/// non-stiff systems whose rates change quickly. It is not suitable for stiff
/// systems (widely separated timescales): stability then forces tiny steps
/// and the solver gives up after 100 000 steps. Use [`simulate_stiff`] for
/// those.
///
/// - `f`: ODE function f(t, y, dy) where dy is filled with derivatives
/// - `y0`: initial state vector (non-empty, finite)
/// - `t_end`: end time
/// - `tol`: error tolerance for adaptive stepping (positive, finite); the
///   per-component error scale is max(tol, tol·|y|)
///
/// Returns every accepted step (t, y), starting at (0, y0) and ending at
/// `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `t_end` or `tol` is not positive
/// and finite or `y0` contains a non-finite value, and
/// [`KimiyaError::ComputationError`] if `y0` is empty or the solver needs more
/// than 100 000 steps.
pub fn simulate_adaptive(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_end: f64,
    tol: f64,
) -> Result<Vec<TrajectoryPoint>> {
    if !is_positive(t_end) {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if !is_positive(tol) {
        return Err(KimiyaError::InvalidInput(
            "tolerance must be positive and finite".into(),
        ));
    }
    validate_state(y0)?;
    let h_init = t_end / 100.0;
    hisab::num::dopri45(f, 0.0, y0, t_end, tol, h_init)
        .map_err(|e| KimiyaError::ComputationError(format!("adaptive ODE solver failed: {e}")))
}

/// Newton iterations allowed per backward-Euler step.
const STIFF_MAX_NEWTON: usize = 20;

/// Newton has converged once the last correction is at most this fraction of
/// the state (∞-norms).
const STIFF_NEWTON_RTOL: f64 = 1e-10;

/// If [`STIFF_MAX_NEWTON`] iterations stall above [`STIFF_NEWTON_RTOL`] — the
/// round-off floor of very stiff systems — the step is still accepted when the
/// last correction is at most this fraction of the state, far below backward
/// Euler's own O(h) truncation error.
const STIFF_ACCEPT_RTOL: f64 = 1e-6;

/// ∞-norm of a vector.
#[inline]
fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |m, &x| m.max(x.abs()))
}

/// Simulate a stiff reaction system using backward Euler (implicit).
///
/// Requires both the ODE function and its Jacobian. Suitable for systems
/// with widely separated timescales (e.g., fast equilibrium + slow reaction).
///
/// Each of the `steps` equal steps of size h = t_end/steps solves
/// y₊ = y + h·f(t₊, y₊) by Newton iteration on (I − h·J)·δ = −(y₊ − y − h·f),
/// using [`hisab::num::gaussian_elimination`]. A step converges when
/// ‖δ‖∞ ≤ 1e-10·‖y₊‖∞ within 20 iterations (or ≤ 1e-6·‖y₊‖∞ if the iteration
/// stalls at its round-off floor); otherwise the simulation fails rather than
/// accept an unconverged state.
///
/// - `f`: ODE function f(t, y, dy)
/// - `jac`: Jacobian function jac(t, y, &mut J) filling the n×n matrix
///   J\[i\]\[j\] = ∂fᵢ/∂yⱼ (the buffer is reused between calls, not zeroed;
///   it must stay n×n)
/// - `y0`: initial state vector (finite)
/// - `t_end`: end time
/// - `steps`: number of time steps (1 ..= 10 000 000)
///
/// Returns exactly two points, `[(0, y0), (t_end, y(t_end))]`: intermediate
/// states are not retained.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `t_end` is not positive and
/// finite, `steps` is 0 or above 10 000 000, `y0` contains a non-finite value,
/// or `jac` resizes the Jacobian buffer; and [`KimiyaError::ComputationError`]
/// if I − h·J is singular, an iterate becomes non-finite, or Newton fails to
/// converge in any step.
pub fn simulate_stiff(
    f: impl Fn(f64, &[f64], &mut [f64]),
    jac: impl Fn(f64, &[f64], &mut Vec<Vec<f64>>),
    y0: &[f64],
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    validate_time_and_steps(t_end, steps)?;
    validate_state(y0)?;
    let dim = y0.len();
    let h = t_end / steps as f64;
    let mut y = y0.to_vec();
    let mut y_next = y0.to_vec();
    let mut f_val = vec![0.0; dim];
    let mut j_mat = vec![vec![0.0; dim]; dim];
    let mut aug = vec![vec![0.0; dim + 1]; dim];

    if dim > 0 {
        for step in 1..=steps {
            let t = if step == steps {
                t_end
            } else {
                h * step as f64
            };
            y_next.copy_from_slice(&y);
            let mut last_correction = f64::INFINITY;
            let mut converged = false;
            for _ in 0..STIFF_MAX_NEWTON {
                f(t, &y_next, &mut f_val);
                jac(t, &y_next, &mut j_mat);
                if j_mat.len() != dim || j_mat.iter().any(|row| row.len() != dim) {
                    return Err(KimiyaError::InvalidInput(format!(
                        "Jacobian must stay {dim}×{dim}"
                    )));
                }
                // Augmented Newton system [I − h·J | −(y₊ − y − h·f)].
                for (i, row) in aug.iter_mut().enumerate() {
                    for (j, entry) in row.iter_mut().take(dim).enumerate() {
                        let identity = if i == j { 1.0 } else { 0.0 };
                        *entry = identity - h * j_mat[i][j];
                    }
                    row[dim] = -(y_next[i] - y[i] - h * f_val[i]);
                }
                let delta = hisab::num::gaussian_elimination(&mut aug).map_err(|e| {
                    KimiyaError::ComputationError(format!(
                        "stiff ODE solver: Newton system singular at t = {t}: {e}"
                    ))
                })?;
                for (yi, di) in y_next.iter_mut().zip(&delta) {
                    *yi += di;
                }
                if !y_next.iter().all(|&v| is_finite(v)) {
                    return Err(KimiyaError::ComputationError(format!(
                        "stiff ODE solver: non-finite state at t = {t}"
                    )));
                }
                last_correction = max_abs(&delta);
                if last_correction <= STIFF_NEWTON_RTOL * max_abs(&y_next) {
                    converged = true;
                    break;
                }
            }
            if !converged && last_correction > STIFF_ACCEPT_RTOL * max_abs(&y_next) {
                return Err(KimiyaError::ComputationError(format!(
                    "stiff ODE solver: Newton did not converge at t = {t} \
                     (last correction {last_correction:e})"
                )));
            }
            y.copy_from_slice(&y_next);
        }
    }
    Ok(vec![(0.0, y0.to_vec()), (t_end, y)])
}

/// Simulate stochastic chemical kinetics using Euler-Maruyama SDE solver.
///
/// Models concentration fluctuations from molecular noise:
/// dX = drift(t, X)·dt + diffusion(t, X)·dW, with one independent Wiener
/// increment per component.
///
/// - `drift`: deterministic part f(t, y, dy) (same as ODE)
/// - `diffusion`: stochastic part g(t, y, dy) (noise amplitude)
/// - `y0`: initial state (finite)
/// - `t_end`: end time
/// - `steps`: number of time steps (1 ..= 10 000 000)
/// - `seed`: random seed for reproducibility
///
/// Returns `steps + 1` points from t = 0 to `t_end`.
///
/// # Errors
///
/// Returns [`KimiyaError::InvalidInput`] if `t_end` is not positive and
/// finite, `steps` is 0 or above 10 000 000, or `y0` contains a non-finite
/// value, and [`KimiyaError::ComputationError`] if the solver fails or the
/// final state is not finite.
pub fn simulate_stochastic(
    drift: impl Fn(f64, &[f64], &mut [f64]),
    diffusion: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_end: f64,
    steps: usize,
    seed: u64,
) -> Result<Vec<TrajectoryPoint>> {
    validate_time_and_steps(t_end, steps)?;
    validate_state(y0)?;
    let mut rng = hisab::num::Pcg32::new(seed, 1);
    finish_trajectory(
        hisab::num::euler_maruyama(drift, diffusion, 0.0, y0, t_end, steps, &mut rng),
        "SDE solver",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert |actual − expected| ≤ tol.
    fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol,
            "{what}: expected {expected}, got {actual} (|Δ| = {:e})",
            (actual - expected).abs()
        );
    }

    fn final_state(traj: &[TrajectoryPoint]) -> &[f64] {
        &traj.last().unwrap().1
    }

    #[test]
    fn first_order_decay() {
        let traj = simulate_first_order(1.0, 0.1, 10.0, 100).unwrap();
        assert_eq!(traj.len(), 101); // 100 steps + initial
        assert_eq!(traj[0], (0.0, vec![1.0]));
        assert_close(traj.last().unwrap().0, 10.0, 1e-12, "end time");
        // RK4, k·h = 0.01: achieved global error 3.1e-11 against exp(-1)
        let final_a = final_state(&traj)[0];
        assert_close(final_a, (-1.0_f64).exp(), 1e-10, "first order at t=10");
    }

    #[test]
    fn first_order_monotonically_decreasing() {
        let traj = simulate_first_order(1.0, 0.5, 5.0, 50).unwrap();
        for window in traj.windows(2) {
            assert!(
                window[1].1[0] < window[0].1[0],
                "concentration should decrease"
            );
        }
    }

    #[test]
    fn second_order_decay() {
        let traj = simulate_second_order(1.0, 0.1, 10.0, 100).unwrap();
        // Analytical: 1/(1/1.0 + 0.1*10) = 1/2 = 0.5
        assert_close(final_state(&traj)[0], 0.5, 1e-9, "second order at t=10");
    }

    #[test]
    fn consecutive_matches_analytic_and_conserves_mass() {
        let (k1, k2) = (0.5, 0.1);
        let traj = simulate_consecutive(1.0, k1, k2, 20.0, 200).unwrap();
        let (mut err_a, mut err_b) = (0.0_f64, 0.0_f64);
        for (t, conc) in &traj {
            let total = conc[0] + conc[1] + conc[2];
            assert_close(total, 1.0, 1e-13, "mass conservation");
            let a = (-k1 * t).exp();
            let b = k1 / (k2 - k1) * ((-k1 * t).exp() - (-k2 * t).exp());
            err_a = err_a.max((conc[0] - a).abs());
            err_b = err_b.max((conc[1] - b).abs());
        }
        // RK4 with k₁·h = 0.05: per-step error ≈ (k₁h)⁵/120, peak ≈ 2e-8
        assert!(
            err_a < 3e-8 && err_b < 3e-8,
            "max errors {err_a:e}, {err_b:e}"
        );
        // At end, most should be C
        assert!(final_state(&traj)[2] > 0.8);
    }

    #[test]
    fn consecutive_intermediate_peaks() {
        let traj = simulate_consecutive(1.0, 1.0, 0.1, 10.0, 100).unwrap();
        // B should rise then fall (intermediate)
        let b_values: Vec<f64> = traj.iter().map(|(_, c)| c[1]).collect();
        let max_b = b_values.iter().copied().fold(0.0_f64, f64::max);
        let final_b = *b_values.last().unwrap();
        assert!(max_b > final_b, "intermediate B should peak then decay");
    }

    #[test]
    fn reversible_reaches_equilibrium() {
        // A ⇌ B, kf=1.0, kr=0.5 → K_eq = kf/kr = 2.0 → [B]/[A] = 2.0
        let traj = simulate_reversible(1.0, 0.0, 1.0, 0.5, 20.0, 200).unwrap();
        let state = final_state(&traj);
        // Remaining relaxation e^{-(kf+kr)·20} ≈ 1e-13; RK4 error ~1e-10.
        assert_close(state[1] / state[0], 2.0, 1e-9, "[B]/[A]");
    }

    #[test]
    fn reversible_mass_conservation() {
        let traj = simulate_reversible(0.8, 0.2, 0.5, 0.3, 10.0, 100).unwrap();
        for (_, conc) in &traj {
            assert_close(conc[0] + conc[1], 1.0, 1e-13, "mass conservation");
        }
    }

    #[test]
    fn michaelis_menten_substrate_depletes() {
        let traj = simulate_michaelis_menten(10.0, 1.0, 2.0, 30.0, 300).unwrap();
        let final_s = final_state(&traj)[0];
        let final_p = final_state(&traj)[1];
        assert!(final_s < 1.0, "substrate should deplete, got {final_s}");
        assert!(final_p > 9.0, "product should accumulate, got {final_p}");
        assert_close(final_s + final_p, 10.0, 1e-11, "mass conservation");
    }

    #[test]
    fn michaelis_menten_mass_conservation() {
        let traj = simulate_michaelis_menten(5.0, 2.0, 1.0, 10.0, 100).unwrap();
        for (_, conc) in &traj {
            assert_close(conc[0] + conc[1], 5.0, 1e-12, "[S]+[P]");
        }
    }

    #[test]
    fn adaptive_first_order() {
        let k = 0.1;
        let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -k * y[0];
        };
        let traj = simulate_adaptive(f, &[1.0], 10.0, 1e-8).unwrap();
        assert_close(traj.last().unwrap().0, 10.0, 1e-12, "end time");
        assert_close(final_state(&traj)[0], (-1.0_f64).exp(), 1e-8, "adaptive");
    }

    #[test]
    fn first_order_negative_concentration_is_error() {
        assert!(simulate_first_order(-1.0, 0.1, 10.0, 100).is_err());
    }

    #[test]
    fn first_order_zero_time_is_error() {
        assert!(simulate_first_order(1.0, 0.1, 0.0, 100).is_err());
    }

    #[test]
    fn invalid_rate_constants_are_rejected() {
        for k in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(simulate_first_order(1.0, k, 10.0, 10).is_err(), "{k}");
            assert!(simulate_second_order(1.0, k, 10.0, 10).is_err(), "{k}");
            assert!(simulate_consecutive(1.0, k, 0.1, 10.0, 10).is_err(), "{k}");
            assert!(simulate_consecutive(1.0, 0.1, k, 10.0, 10).is_err(), "{k}");
            assert!(simulate_reversible(1.0, 0.0, k, 0.1, 10.0, 10).is_err());
            assert!(simulate_reversible(1.0, 0.0, 0.1, k, 10.0, 10).is_err());
            assert!(simulate_michaelis_menten(1.0, k, 1.0, 10.0, 10).is_err());
            assert!(simulate_michaelis_menten(1.0, 1.0, k, 10.0, 10).is_err());
        }
        // Zero rate constants are valid (nothing happens).
        let traj = simulate_first_order(1.0, 0.0, 10.0, 10).unwrap();
        assert_eq!(final_state(&traj)[0], 1.0);
    }

    #[test]
    fn invalid_initial_concentrations_and_times_are_rejected() {
        for c in [f64::NAN, f64::INFINITY] {
            assert!(simulate_first_order(c, 0.1, 10.0, 10).is_err());
            assert!(simulate_second_order(c, 0.1, 10.0, 10).is_err());
            assert!(simulate_consecutive(c, 0.1, 0.1, 10.0, 10).is_err());
            assert!(simulate_reversible(1.0, c, 0.1, 0.1, 10.0, 10).is_err());
            assert!(simulate_michaelis_menten(c, 1.0, 1.0, 10.0, 10).is_err());
        }
        for t in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(simulate_first_order(1.0, 0.1, t, 10).is_err());
            assert!(simulate_second_order(1.0, 0.1, t, 10).is_err());
            assert!(simulate_consecutive(1.0, 0.1, 0.1, t, 10).is_err());
            assert!(simulate_reversible(1.0, 0.0, 0.1, 0.1, t, 10).is_err());
            assert!(simulate_michaelis_menten(1.0, 1.0, 1.0, t, 10).is_err());
            assert!(simulate_adaptive(|_, _, _| {}, &[1.0], t, 1e-8).is_err());
            assert!(simulate_stiff(|_, _, _| {}, |_, _, _| {}, &[1.0], t, 10).is_err());
            assert!(simulate_stochastic(|_, _, _| {}, |_, _, _| {}, &[1.0], t, 10, 1).is_err());
        }
    }

    #[test]
    fn huge_step_counts_are_rejected() {
        // hisab computes `Vec::with_capacity(steps + 1)`: usize::MAX overflowed.
        for steps in [usize::MAX, MAX_STEPS + 1] {
            assert!(simulate_first_order(1.0, 0.1, 10.0, steps).is_err());
            assert!(simulate_second_order(1.0, 0.1, 10.0, steps).is_err());
            assert!(simulate_consecutive(1.0, 0.1, 0.1, 10.0, steps).is_err());
            assert!(simulate_reversible(1.0, 0.0, 0.1, 0.1, 10.0, steps).is_err());
            assert!(simulate_michaelis_menten(1.0, 1.0, 1.0, 10.0, steps).is_err());
            assert!(simulate_stiff(|_, _, _| {}, |_, _, _| {}, &[1.0], 1.0, steps).is_err());
            assert!(
                simulate_stochastic(|_, _, _| {}, |_, _, _| {}, &[1.0], 1.0, steps, 1).is_err()
            );
        }
    }

    #[test]
    fn unstable_rk4_is_error_not_nan() {
        // k·h = 100 is far outside RK4's stability region: the state overflows.
        let result = simulate_first_order(1.0, 1e3, 10.0, 100);
        assert!(matches!(result, Err(KimiyaError::ComputationError(_))));
    }

    #[test]
    fn adaptive_rejects_bad_tolerance_and_state() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        // tol = 0 used to burn 100 000 rejected steps before failing.
        for tol in [0.0, -1e-8, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                simulate_adaptive(f, &[1.0], 1.0, tol),
                Err(KimiyaError::InvalidInput(_))
            ));
        }
        assert!(simulate_adaptive(f, &[f64::NAN], 1.0, 1e-8).is_err());
    }

    // ── Stiff (backward Euler) ───────────────────────────────────────

    #[test]
    fn stiff_linear_matches_backward_euler_recurrence() {
        // y' = −k·y; backward Euler gives y_n = y_0 / (1 + h·k)^n exactly.
        let k = 50.0;
        let f = move |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -k * y[0];
        let jac = move |_t: f64, _y: &[f64], j: &mut Vec<Vec<f64>>| j[0][0] = -k;
        let traj = simulate_stiff(f, jac, &[1.0], 1.0, 100).unwrap();
        assert_eq!(traj.len(), 2, "only the end points are returned");
        assert_eq!(traj[0], (0.0, vec![1.0]));
        assert_eq!(traj[1].0, 1.0);
        let expected = 1.5_f64.powi(-100);
        let got = traj[1].1[0];
        assert!(
            ((got - expected) / expected).abs() < 1e-12,
            "{got} vs {expected}"
        );
    }

    #[test]
    fn stiff_nonlinear_matches_backward_euler_recurrence() {
        // y' = −k·y²; each backward Euler step solves h·k·y² + y − y_old = 0.
        let (k, h, steps) = (3.0, 0.05, 40);
        let f = move |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -k * y[0] * y[0];
        let jac = move |_t: f64, y: &[f64], j: &mut Vec<Vec<f64>>| j[0][0] = -2.0 * k * y[0];
        let traj = simulate_stiff(f, jac, &[2.0], h * steps as f64, steps).unwrap();
        let mut expected = 2.0_f64;
        for _ in 0..steps {
            let a = h * k;
            // Stable form of the positive root (−1 + √(1 + 4a·y))/(2a).
            expected = 2.0 * expected / (1.0 + (1.0 + 4.0 * a * expected).sqrt());
        }
        let got = traj[1].1[0];
        assert!(
            ((got - expected) / expected).abs() < 1e-12,
            "{got} vs {expected}"
        );
    }

    #[test]
    fn stiff_robertson_conserves_mass() {
        // Robertson's stiff kinetics: A → B (0.04), 2B → B + C (3e7), B + C → A + C (1e4)
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            let r1 = 0.04 * y[0];
            let r2 = 3e7 * y[1] * y[1];
            let r3 = 1e4 * y[1] * y[2];
            dy[0] = -r1 + r3;
            dy[1] = r1 - r2 - r3;
            dy[2] = r2;
        };
        let jac = |_t: f64, y: &[f64], j: &mut Vec<Vec<f64>>| {
            j[0] = vec![-0.04, 1e4 * y[2], 1e4 * y[1]];
            j[1] = vec![0.04, -6e7 * y[1] - 1e4 * y[2], -1e4 * y[1]];
            j[2] = vec![0.0, 6e7 * y[1], 0.0];
        };
        let traj = simulate_stiff(f, jac, &[1.0, 0.0, 0.0], 40.0, 400).unwrap();
        let y = &traj[1].1;
        assert_close(y[0] + y[1] + y[2], 1.0, 1e-12, "mass conservation");
        // Reference solution at t = 40: y₁ ≈ 0.7158, y₂ ≈ 9.185e-6
        // (backward Euler with h = 0.1 is only first-order accurate).
        assert_close(y[0], 0.7158, 5e-3, "y1(40)");
        assert!(y[1] > 0.0 && y[1] < 2e-5, "y2(40) = {}", y[1]);
    }

    #[test]
    fn stiff_singular_newton_system_is_error() {
        // J = 1/h makes I − h·J singular; hisab used to break out of Newton
        // and return the unconverged iterate as Ok.
        let (t_end, steps) = (1.0, 10);
        let c = steps as f64 / t_end;
        let f = move |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = c * y[0] + 1.0;
        let jac = move |_t: f64, _y: &[f64], j: &mut Vec<Vec<f64>>| j[0][0] = c;
        let result = simulate_stiff(f, jac, &[1.0], t_end, steps);
        assert!(matches!(result, Err(KimiyaError::ComputationError(_))));
    }

    #[test]
    fn stiff_non_convergent_newton_is_error() {
        // A wrong-sign Jacobian makes Newton diverge geometrically; hisab
        // accepted the 20th iterate as the step result.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        let jac = |_t: f64, _y: &[f64], j: &mut Vec<Vec<f64>>| j[0][0] = 3.0;
        let result = simulate_stiff(f, jac, &[1.0], 1.0, 1);
        assert!(matches!(result, Err(KimiyaError::ComputationError(_))));
    }

    #[test]
    fn stiff_resized_jacobian_is_error_not_panic() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        let jac = |_t: f64, _y: &[f64], j: &mut Vec<Vec<f64>>| j.clear();
        assert!(matches!(
            simulate_stiff(f, jac, &[1.0], 1.0, 1),
            Err(KimiyaError::InvalidInput(_))
        ));
    }

    #[test]
    fn stiff_rejects_non_finite_state() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        let jac = |_t: f64, _y: &[f64], j: &mut Vec<Vec<f64>>| j[0][0] = -1.0;
        assert!(simulate_stiff(f, jac, &[f64::NAN], 1.0, 1).is_err());
        let f_nan = |_t: f64, _y: &[f64], dy: &mut [f64]| dy[0] = f64::NAN;
        assert!(simulate_stiff(f_nan, jac, &[1.0], 1.0, 1).is_err());
    }

    // ── Stochastic ───────────────────────────────────────────────────

    #[test]
    fn stochastic_is_reproducible_and_zero_noise_is_euler() {
        let drift = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        let no_noise = |_t: f64, _y: &[f64], g: &mut [f64]| g[0] = 0.0;
        let traj = simulate_stochastic(drift, no_noise, &[1.0], 1.0, 10, 7).unwrap();
        // Forward Euler: (1 − h)^n
        assert_close(final_state(&traj)[0], 0.9_f64.powi(10), 1e-14, "Euler");
        let noise = |_t: f64, _y: &[f64], g: &mut [f64]| g[0] = 0.1;
        let a = simulate_stochastic(drift, noise, &[1.0], 1.0, 50, 42).unwrap();
        let b = simulate_stochastic(drift, noise, &[1.0], 1.0, 50, 42).unwrap();
        assert_eq!(a, b, "same seed → same path");
    }

    #[test]
    fn stochastic_blow_up_is_error_not_inf() {
        let drift = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = 1e3 * y[0] * y[0];
        let no_noise = |_t: f64, _y: &[f64], g: &mut [f64]| g[0] = 0.0;
        assert!(matches!(
            simulate_stochastic(drift, no_noise, &[1.0], 1.0, 100, 1),
            Err(KimiyaError::ComputationError(_))
        ));
        assert!(simulate_stochastic(drift, no_noise, &[f64::NAN], 1.0, 10, 1).is_err());
    }
}
