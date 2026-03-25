//! # Kimiya
//!
//! **Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry
//! engine for the AGNOS ecosystem.
//!
//! Full periodic table (118 elements), molecules, reactions, kinetics, gas laws,
//! solutions, thermochemistry, electrochemistry, spectroscopy, organic chemistry,
//! thermochemical data tables, reaction dynamics, and data fitting.
//! Built on [`hisab`] for numerical methods (root finding, ODE solvers, least squares).

pub mod electrochemistry;
pub mod element;
pub mod error;
pub mod fitting;
pub mod gas;
pub mod kinetics;
pub mod molecule;
pub mod organic;
pub mod reaction;
pub mod reaction_dynamics;
pub mod solution;
pub mod spectroscopy;
pub mod thermo;
pub mod thermochem;

#[cfg(feature = "logging")]
pub mod logging;

pub use electrochemistry::{FARADAY, nernst_potential};
pub use element::{ELEMENTS, Element, ElementCategory};
pub use error::{KimiyaError, Result};
pub use fitting::fit_arrhenius;
pub use gas::ideal_gas_pressure;
pub use kinetics::arrhenius_rate;
pub use molecule::Molecule;
pub use reaction_dynamics::simulate_first_order;
pub use solution::{molarity, ph_from_h_concentration, weak_acid_ph};
pub use thermochem::{THERMOCHEM_DATA, enthalpy_change_cp, lookup_thermochem};
