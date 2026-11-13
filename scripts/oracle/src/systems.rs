//! The test systems of the M5a golden rows (`numerics.tsv`, `reaction_dynamics.tsv`): the Rust closures that the
//! oracle's custom-system solvers (`simulate_adaptive`, `simulate_stiff`, `simulate_stochastic`, and hisab 1.4.0's
//! routines underneath) integrate. A golden row names a system by an `s:` token and its parameters by `v:n` + `f:`
//! tokens; the parity harnesses (`tests/parity.tcyr`, `tests/numerics.tcyr`) and `validate.py` implement the same
//! systems with the same operation order, so the three agree bit for bit. Each is a chemical mechanism (the five
//! built-in rate laws, Robertson's stiff kinetics, the Brusselator oscillator), or an edge of the solvers the 1.2 tests
//! pin (a derivative that is NaN or never written, a singular or wrong-sign Newton system, a blow-up). The formulas are
//! documented in `scripts/oracle/README.md`; keep the four copies in step.

use crate::tsv::{Tok, f, s};

/// The right-hand side dy = f(t, y) of system `name` with parameters `p`. A system that writes nothing (`zero`) leaves
/// the solver's buffer as it was, as Rust's `|_, _, _| {}` does.
#[allow(clippy::needless_range_loop)] // the 1.2 tests' indexed form, as the other three copies write it
pub fn rhs(name: &str, p: &[f64], t: f64, y: &[f64], dy: &mut [f64]) {
    match name {
        // dy_i = −k·y_i (the 1.2 tests' first-order decay, per component)
        "decay" => {
            for i in 0..y.len() {
                dy[i] = -p[0] * y[i];
            }
        }
        // dy_i = −k·y_i²
        "decay2" => {
            for i in 0..y.len() {
                dy[i] = -p[0] * y[i] * y[i];
            }
        }
        // A → B → C (simulate_consecutive's closure)
        "consecutive" => {
            dy[0] = -p[0] * y[0];
            dy[1] = p[0] * y[0] - p[1] * y[1];
            dy[2] = p[1] * y[1];
        }
        // A ⇌ B (simulate_reversible's closure)
        "reversible" => {
            dy[0] = -p[0] * y[0] + p[1] * y[1];
            dy[1] = p[0] * y[0] - p[1] * y[1];
        }
        // Michaelis–Menten with depletion (simulate_michaelis_menten's closure)
        "michaelis_menten" => {
            let rate = p[0] * y[0] / (p[1] + y[0]);
            dy[0] = -rate;
            dy[1] = rate;
        }
        // Robertson's stiff kinetics (the 1.2 test stiff_robertson_conserves_mass)
        "robertson" => {
            let r1 = 0.04 * y[0];
            let r2 = 3e7 * y[1] * y[1];
            let r3 = 1e4 * y[1] * y[2];
            dy[0] = -r1 + r3;
            dy[1] = r1 - r2 - r3;
            dy[2] = r2;
        }
        // The Brusselator, A → X, B + X → Y + D, 2X + Y → 3X, X → E: a limit cycle for B > 1 + A²
        "brusselator" => {
            let xxy = y[0] * y[0] * y[1];
            dy[0] = p[0] + xxy - (p[1] + 1.0) * y[0];
            dy[1] = p[1] * y[0] - xxy;
        }
        // dy = −k·y + c·t: an explicitly time-dependent rate (the stage times)
        "forced" => {
            dy[0] = -p[0] * y[0] + p[1] * t;
        }
        // dy = c·y + 1 (stiff_singular_newton_system_is_error: c = steps/t_end makes I − h·J singular)
        "affine" => {
            dy[0] = p[0] * y[0] + 1.0;
        }
        // dy = NaN (stiff_rejects_non_finite_state's f_nan)
        "nan" => {
            for d in dy.iter_mut() {
                *d = f64::NAN;
            }
        }
        // dy = c·y² (stochastic_blow_up_is_error_not_inf: c = 1e3)
        "blowup" => {
            dy[0] = p[0] * y[0] * y[0];
        }
        // dy = −y with a wrong-sign Jacobian (stiff_non_convergent_newton_is_error)
        "badjac" => {
            dy[0] = -y[0];
        }
        // dy = −k·y with the Jacobian −c: Newton converges linearly at |h·(c − k)|/(1 + h·c) per iteration (the
        // 20-iteration limit and both Newton tolerances)
        "slowjac" => {
            dy[0] = -p[0] * y[0];
        }
        // writes nothing (the 1.2 tests' `|_, _, _| {}`)
        "zero" => {}
        _ => panic!("no system {name}"),
    }
}

/// The Jacobian J_ij = ∂f_i/∂y_j of system `name`, written into Rust's reused, never re-zeroed `Vec<Vec<f64>>`: a
/// system writes the entries it names and leaves the rest.
#[allow(clippy::needless_range_loop)]
pub fn jac(name: &str, p: &[f64], _t: f64, y: &[f64], j: &mut [Vec<f64>]) {
    match name {
        "decay" | "forced" => {
            for i in 0..y.len() {
                j[i][i] = -p[0];
            }
        }
        "decay2" => {
            for i in 0..y.len() {
                j[i][i] = -2.0 * p[0] * y[i];
            }
        }
        "consecutive" => {
            j[0][0] = -p[0];
            j[1][0] = p[0];
            j[1][1] = -p[1];
            j[2][1] = p[1];
        }
        "reversible" => {
            j[0][0] = -p[0];
            j[0][1] = p[1];
            j[1][0] = p[0];
            j[1][1] = -p[1];
        }
        "michaelis_menten" => {
            let s_ = p[1] + y[0];
            let d = p[0] * p[1] / (s_ * s_);
            j[0][0] = -d;
            j[1][0] = d;
        }
        "robertson" => {
            j[0] = vec![-0.04, 1e4 * y[2], 1e4 * y[1]];
            j[1] = vec![0.04, -6e7 * y[1] - 1e4 * y[2], -1e4 * y[1]];
            j[2] = vec![0.0, 6e7 * y[1], 0.0];
        }
        "brusselator" => {
            let xy2 = 2.0 * y[0] * y[1];
            let xx = y[0] * y[0];
            j[0][0] = xy2 - (p[1] + 1.0);
            j[0][1] = xx;
            j[1][0] = p[1] - xy2;
            j[1][1] = -xx;
        }
        "affine" => j[0][0] = p[0],
        "nan" => j[0][0] = -1.0,
        "blowup" => j[0][0] = 2.0 * p[0] * y[0],
        "badjac" => j[0][0] = 3.0,
        "slowjac" => j[0][0] = -p[1],
        "zero" => {}
        _ => panic!("no Jacobian {name}"),
    }
}

/// The diffusion (noise amplitude) g(t, y) of `name`: `const` g_i = σ, `zero` g_i = 0, `mult` g_i = σ·y_i, `none`
/// writes nothing.
#[allow(clippy::needless_range_loop)]
pub fn diffusion(name: &str, p: &[f64], _t: f64, y: &[f64], g: &mut [f64]) {
    match name {
        "const" => {
            for gi in g.iter_mut() {
                *gi = p[0];
            }
        }
        "zero" => {
            for gi in g.iter_mut() {
                *gi = 0.0;
            }
        }
        "mult" => {
            for i in 0..y.len() {
                g[i] = p[0] * y[i];
            }
        }
        "none" => {}
        _ => panic!("no diffusion {name}"),
    }
}

/// `s:name v:n f:p…`: a system and its parameters.
pub fn sys_toks(name: &str, p: &[f64]) -> Vec<Tok> {
    let mut v = vec![s(name), Tok::V(p.len())];
    v.extend(p.iter().map(|&x| f(x)));
    v
}

/// `v:n f:…`: a state vector.
pub fn vec_toks(y: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(y.len())];
    v.extend(y.iter().map(|&x| f(x)));
    v
}

/// A trajectory's tokens: `v:points`, then each point's `f:t v:n f:y…`.
pub fn traj_toks(tr: &[(f64, Vec<f64>)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(tr.len())];
    for (t, y) in tr {
        v.push(f(*t));
        v.extend(vec_toks(y));
    }
    v
}

/// A trajectory row's note tail: the point count and the last point, readably.
pub fn traj_note(tr: &[(f64, Vec<f64>)]) -> String {
    match tr.last() {
        Some((t, y)) => format!("{} points, last ({t:?}, {y:?})", tr.len()),
        None => "0 points".to_string(),
    }
}
