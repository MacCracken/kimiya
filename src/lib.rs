//! # Kimiya
//!
//! **Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry
//! engine for the AGNOS ecosystem.
//!
//! Provides elements, molecules, reactions, kinetics, gas laws, solutions,
//! thermochemistry, electrochemistry, and thermochemical data tables.
//! Built on [`hisab`] for math.

pub mod electrochemistry;
pub mod element;
pub mod error;
pub mod gas;
pub mod kinetics;
pub mod molecule;
pub mod reaction;
pub mod solution;
pub mod thermo;
pub mod thermochem;

#[cfg(feature = "logging")]
pub mod logging;

pub use electrochemistry::{FARADAY, nernst_potential};
pub use element::{ELEMENTS, Element, ElementCategory};
pub use error::{KimiyaError, Result};
pub use gas::ideal_gas_pressure;
pub use kinetics::arrhenius_rate;
pub use molecule::Molecule;
pub use solution::{molarity, ph_from_h_concentration};
pub use thermochem::{THERMOCHEM_DATA, lookup_thermochem};
