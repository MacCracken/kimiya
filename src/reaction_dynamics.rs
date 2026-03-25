//! Reaction dynamics — time-evolution of chemical systems via ODE integration.
//!
//! Uses [`hisab`] ODE solvers (RK4, Dormand-Prince 4/5) to simulate
//! concentration trajectories for reaction kinetics.

use crate::error::{KimiyaError, Result};

/// A point in a kinetics trajectory: (time, concentrations).
pub type TrajectoryPoint = (f64, Vec<f64>);

/// Simulate first-order decay: d\[A\]/dt = -k\[A\]
///
/// Returns trajectory of \[A\] over time.
///
/// # Errors
///
/// Returns error if parameters are invalid or ODE solver fails.
pub fn simulate_first_order(
    initial_concentration: f64,
    rate_constant: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if initial_concentration < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    let k = rate_constant;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k * y[0];
    };
    hisab::num::rk4_trajectory(f, 0.0, &[initial_concentration], t_end, steps)
        .map_err(|e| KimiyaError::ComputationError(format!("ODE solver failed: {e}")))
}

/// Simulate second-order decay: d\[A\]/dt = -k\[A\]²
///
/// Returns trajectory of \[A\] over time.
///
/// # Errors
///
/// Returns error if parameters are invalid or ODE solver fails.
pub fn simulate_second_order(
    initial_concentration: f64,
    rate_constant: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if initial_concentration < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    let k = rate_constant;
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k * y[0] * y[0];
    };
    hisab::num::rk4_trajectory(f, 0.0, &[initial_concentration], t_end, steps)
        .map_err(|e| KimiyaError::ComputationError(format!("ODE solver failed: {e}")))
}

/// Simulate consecutive reactions: A →(k₁) B →(k₂) C
///
/// State vector: \[A\], \[B\], \[C\]
///
/// Returns trajectory of all three species over time.
///
/// # Errors
///
/// Returns error if parameters are invalid or ODE solver fails.
pub fn simulate_consecutive(
    initial_a: f64,
    k1: f64,
    k2: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if initial_a < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "initial concentration must be non-negative".into(),
        ));
    }
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k1 * y[0]; // d[A]/dt = -k1[A]
        dy[1] = k1 * y[0] - k2 * y[1]; // d[B]/dt = k1[A] - k2[B]
        dy[2] = k2 * y[1]; // d[C]/dt = k2[B]
    };
    hisab::num::rk4_trajectory(f, 0.0, &[initial_a, 0.0, 0.0], t_end, steps)
        .map_err(|e| KimiyaError::ComputationError(format!("ODE solver failed: {e}")))
}

/// Simulate reversible reaction: A ⇌ B with forward rate k_f and reverse rate k_r
///
/// State vector: \[A\], \[B\]
///
/// # Errors
///
/// Returns error if parameters are invalid or ODE solver fails.
pub fn simulate_reversible(
    initial_a: f64,
    initial_b: f64,
    k_forward: f64,
    k_reverse: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if initial_a < 0.0 || initial_b < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "concentrations must be non-negative".into(),
        ));
    }
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -k_forward * y[0] + k_reverse * y[1]; // d[A]/dt
        dy[1] = k_forward * y[0] - k_reverse * y[1]; // d[B]/dt
    };
    hisab::num::rk4_trajectory(f, 0.0, &[initial_a, initial_b], t_end, steps)
        .map_err(|e| KimiyaError::ComputationError(format!("ODE solver failed: {e}")))
}

/// Simulate Michaelis-Menten enzyme kinetics with substrate depletion.
///
/// Full model: E + S ⇌(k₁/k₋₁) ES →(k₂) E + P
///
/// Simplified to: d\[S\]/dt = -Vmax·\[S\]/(Km + \[S\]), d\[P\]/dt = Vmax·\[S\]/(Km + \[S\])
///
/// State vector: \[S\], \[P\]
///
/// # Errors
///
/// Returns error if parameters are invalid or ODE solver fails.
pub fn simulate_michaelis_menten(
    initial_substrate: f64,
    v_max: f64,
    km: f64,
    t_end: f64,
    steps: usize,
) -> Result<Vec<TrajectoryPoint>> {
    if initial_substrate < 0.0 {
        return Err(KimiyaError::InvalidConcentration(
            "substrate concentration must be non-negative".into(),
        ));
    }
    if v_max <= 0.0 {
        return Err(KimiyaError::InvalidInput("Vmax must be positive".into()));
    }
    if km <= 0.0 {
        return Err(KimiyaError::InvalidInput("Km must be positive".into()));
    }
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    if steps == 0 {
        return Err(KimiyaError::InvalidInput("steps must be at least 1".into()));
    }
    let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
        let rate = v_max * y[0] / (km + y[0]);
        dy[0] = -rate; // d[S]/dt
        dy[1] = rate; // d[P]/dt
    };
    hisab::num::rk4_trajectory(f, 0.0, &[initial_substrate, 0.0], t_end, steps)
        .map_err(|e| KimiyaError::ComputationError(format!("ODE solver failed: {e}")))
}

/// Simulate a custom reaction system with adaptive step size.
///
/// Uses Dormand-Prince 4/5 (adaptive RK) for stiff or rapidly-changing systems.
///
/// - `f`: ODE function f(t, y, dy) where dy is filled with derivatives
/// - `y0`: initial state vector
/// - `t_end`: end time
/// - `tol`: error tolerance for adaptive stepping
///
/// # Errors
///
/// Returns error if ODE solver fails.
pub fn simulate_adaptive(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_end: f64,
    tol: f64,
) -> Result<Vec<TrajectoryPoint>> {
    if t_end <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "end time must be positive".into(),
        ));
    }
    let h_init = t_end / 100.0;
    hisab::num::dopri45(f, 0.0, y0, t_end, tol, h_init)
        .map_err(|e| KimiyaError::ComputationError(format!("adaptive ODE solver failed: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_order_decay() {
        let traj = simulate_first_order(1.0, 0.1, 10.0, 100).unwrap();
        assert_eq!(traj.len(), 101); // 100 steps + initial
        // At t=0, [A] should be ~1.0
        assert!((traj[0].1[0] - 1.0).abs() < 0.01);
        // At t=10, [A] = exp(-0.1*10) = exp(-1) ≈ 0.368
        let final_a = traj.last().unwrap().1[0];
        assert!(
            (final_a - 0.368).abs() < 0.01,
            "first-order at t=10 should be ~0.368, got {final_a}"
        );
    }

    #[test]
    fn first_order_monotonically_decreasing() {
        let traj = simulate_first_order(1.0, 0.5, 5.0, 50).unwrap();
        for window in traj.windows(2) {
            assert!(
                window[1].1[0] <= window[0].1[0] + 1e-10,
                "concentration should not increase"
            );
        }
    }

    #[test]
    fn second_order_decay() {
        let traj = simulate_second_order(1.0, 0.1, 10.0, 100).unwrap();
        let final_a = traj.last().unwrap().1[0];
        // Analytical: 1/(1/1.0 + 0.1*10) = 1/2 = 0.5
        assert!(
            (final_a - 0.5).abs() < 0.01,
            "second-order at t=10 should be ~0.5, got {final_a}"
        );
    }

    #[test]
    fn consecutive_mass_conservation() {
        let traj = simulate_consecutive(1.0, 0.5, 0.1, 20.0, 200).unwrap();
        // Mass conservation: [A] + [B] + [C] = 1.0 at all times
        for (_, conc) in &traj {
            let total = conc[0] + conc[1] + conc[2];
            assert!(
                (total - 1.0).abs() < 0.01,
                "mass should be conserved, got {total}"
            );
        }
        // At end, most should be C
        let final_c = traj.last().unwrap().1[2];
        assert!(
            final_c > 0.8,
            "most product should be C at end, got {final_c}"
        );
    }

    #[test]
    fn consecutive_intermediate_peaks() {
        let traj = simulate_consecutive(1.0, 1.0, 0.1, 10.0, 100).unwrap();
        // B should rise then fall (intermediate)
        let b_values: Vec<f64> = traj.iter().map(|(_, c)| c[1]).collect();
        let max_b = b_values.iter().cloned().fold(0.0_f64, f64::max);
        let final_b = *b_values.last().unwrap();
        assert!(max_b > final_b, "intermediate B should peak then decay");
    }

    #[test]
    fn reversible_reaches_equilibrium() {
        // A ⇌ B, kf=1.0, kr=0.5 → K_eq = kf/kr = 2.0 → [B]/[A] = 2.0
        let traj = simulate_reversible(1.0, 0.0, 1.0, 0.5, 20.0, 200).unwrap();
        let final_state = &traj.last().unwrap().1;
        let ratio = final_state[1] / final_state[0];
        assert!(
            (ratio - 2.0).abs() < 0.1,
            "[B]/[A] at equilibrium should be ~2.0, got {ratio}"
        );
    }

    #[test]
    fn reversible_mass_conservation() {
        let traj = simulate_reversible(0.8, 0.2, 0.5, 0.3, 10.0, 100).unwrap();
        for (_, conc) in &traj {
            let total = conc[0] + conc[1];
            assert!(
                (total - 1.0).abs() < 0.01,
                "mass should be conserved, got {total}"
            );
        }
    }

    #[test]
    fn michaelis_menten_substrate_depletes() {
        let traj = simulate_michaelis_menten(10.0, 1.0, 2.0, 30.0, 300).unwrap();
        let final_s = traj.last().unwrap().1[0];
        let final_p = traj.last().unwrap().1[1];
        assert!(final_s < 1.0, "substrate should deplete, got {final_s}");
        assert!(final_p > 9.0, "product should accumulate, got {final_p}");
        // Mass conservation
        assert!((final_s + final_p - 10.0).abs() < 0.1);
    }

    #[test]
    fn michaelis_menten_mass_conservation() {
        let traj = simulate_michaelis_menten(5.0, 2.0, 1.0, 10.0, 100).unwrap();
        for (_, conc) in &traj {
            assert!(
                (conc[0] + conc[1] - 5.0).abs() < 0.1,
                "[S]+[P] should be conserved"
            );
        }
    }

    #[test]
    fn adaptive_first_order() {
        let k = 0.1;
        let f = move |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -k * y[0];
        };
        let traj = simulate_adaptive(f, &[1.0], 10.0, 1e-8).unwrap();
        let final_a = traj.last().unwrap().1[0];
        assert!(
            (final_a - (-1.0_f64).exp()).abs() < 0.001,
            "adaptive should match analytical, got {final_a}"
        );
    }

    #[test]
    fn first_order_negative_concentration_is_error() {
        assert!(simulate_first_order(-1.0, 0.1, 10.0, 100).is_err());
    }

    #[test]
    fn first_order_zero_time_is_error() {
        assert!(simulate_first_order(1.0, 0.1, 0.0, 100).is_err());
    }
}
