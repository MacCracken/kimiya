//! # Kimiya
//!
//! **Kimiya** (كيمياء — Arabic for "alchemy", the root of "chemistry") — chemistry
//! engine for the AGNOS ecosystem.
//!
//! Provides elements, molecules, reactions, kinetics, gas laws, solutions,
//! and thermochemistry. Built on [`hisab`] for math.

pub mod error;
pub mod element;
pub mod molecule;
pub mod reaction;
pub mod kinetics;
pub mod gas;
pub mod solution;
pub mod thermo;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::{KimiyaError, Result};
pub use element::{Element, ElementCategory, ELEMENTS};
pub use molecule::Molecule;
pub use kinetics::arrhenius_rate;
pub use gas::ideal_gas_pressure;
pub use solution::{ph_from_h_concentration, molarity};
