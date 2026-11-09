//! Molecules — composition as (element, count) pairs, molecular weight and
//! formula strings.

use crate::element::lookup_by_number;
use crate::error::{KimiyaError, Result};
use serde::{Deserialize, Serialize};

/// An atom count in a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    /// Atomic number Z of the element (see [`crate::element::lookup_by_number`]).
    pub element_number: u8,
    /// Number of atoms of this element in one molecule (formula unit).
    pub count: u32,
}

/// A molecule defined by its constituent atoms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    /// Constituent atoms in the order given at construction. Repeated
    /// entries for the same element are kept as separate entries.
    pub atoms: Vec<Atom>,
}

impl Molecule {
    /// Create a molecule from a list of (atomic_number, count) pairs.
    #[must_use]
    pub fn new(atoms: &[(u8, u32)]) -> Self {
        Self {
            atoms: atoms
                .iter()
                .map(|&(n, c)| Atom {
                    element_number: n,
                    count: c,
                })
                .collect(),
        }
    }

    /// Water: H₂O
    #[must_use]
    pub fn water() -> Self {
        Self::new(&[(1, 2), (8, 1)])
    }

    /// Carbon dioxide: CO₂
    #[must_use]
    pub fn carbon_dioxide() -> Self {
        Self::new(&[(6, 1), (8, 2)])
    }

    /// Methane: CH₄
    #[must_use]
    pub fn methane() -> Self {
        Self::new(&[(6, 1), (1, 4)])
    }

    /// Glucose: C₆H₁₂O₆
    #[must_use]
    pub fn glucose() -> Self {
        Self::new(&[(6, 6), (1, 12), (8, 6)])
    }

    /// Molecular weight in g/mol (sum of atomic masses × count).
    ///
    /// # Errors
    ///
    /// Returns [`KimiyaError::InvalidElement`] if any atom references an unknown element.
    pub fn molecular_weight(&self) -> Result<f64> {
        let mut total = 0.0;
        for atom in &self.atoms {
            let element = lookup_by_number(atom.element_number).ok_or_else(|| {
                KimiyaError::InvalidElement(format!(
                    "unknown atomic number {}",
                    atom.element_number
                ))
            })?;
            total += element.atomic_mass * atom.count as f64;
        }
        Ok(total)
    }

    /// Molecular formula string (e.g., "H2O", "CO2", "C6H12O6").
    ///
    /// Atoms are written in the order stored in [`Self::atoms`] (no Hill
    /// ordering, no merging of repeated elements, no reduction to the
    /// empirical formula). A count of 1 is written without a subscript and
    /// atoms with a count of 0 are omitted.
    ///
    /// # Errors
    ///
    /// Returns [`KimiyaError::InvalidElement`] if any atom references an unknown element.
    pub fn formula(&self) -> Result<String> {
        let mut s = String::with_capacity(self.atoms.len() * 4);
        for atom in &self.atoms {
            let element = lookup_by_number(atom.element_number).ok_or_else(|| {
                KimiyaError::InvalidElement(format!(
                    "unknown atomic number {}",
                    atom.element_number
                ))
            })?;
            if atom.count == 0 {
                continue;
            }
            s.push_str(element.symbol);
            if atom.count > 1 {
                s.push_str(&atom.count.to_string());
            }
        }
        Ok(s)
    }

    /// Total number of atoms in the molecule.
    ///
    /// Saturates at [`u32::MAX`] instead of overflowing.
    #[must_use]
    pub fn total_atoms(&self) -> u32 {
        self.atoms
            .iter()
            .fold(0u32, |acc, a| acc.saturating_add(a.count))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_molecular_weight() {
        let mw = Molecule::water().molecular_weight().unwrap();
        assert!(
            (mw - 18.015).abs() < 0.01,
            "H2O should be ~18.015, got {mw}"
        );
    }

    #[test]
    fn co2_molecular_weight() {
        let mw = Molecule::carbon_dioxide().molecular_weight().unwrap();
        assert!((mw - 44.009).abs() < 0.02, "CO2 should be ~44.01, got {mw}");
    }

    #[test]
    fn water_formula() {
        assert_eq!(Molecule::water().formula().unwrap(), "H2O");
    }

    #[test]
    fn co2_formula() {
        assert_eq!(Molecule::carbon_dioxide().formula().unwrap(), "CO2");
    }

    #[test]
    fn glucose_molecular_weight() {
        let mw = Molecule::glucose().molecular_weight().unwrap();
        assert!(
            (mw - 180.156).abs() < 0.1,
            "glucose should be ~180.16, got {mw}"
        );
    }

    #[test]
    fn methane_total_atoms() {
        assert_eq!(Molecule::methane().total_atoms(), 5); // 1C + 4H
    }

    #[test]
    fn serde_roundtrip() {
        let m = Molecule::water();
        let json = serde_json::to_string(&m).unwrap();
        let back: Molecule = serde_json::from_str(&json).unwrap();
        assert!((back.molecular_weight().unwrap() - m.molecular_weight().unwrap()).abs() < 0.001);
    }

    #[test]
    fn formula_skips_zero_count_atoms() {
        assert_eq!(Molecule::new(&[(1, 0), (8, 1)]).formula().unwrap(), "O");
        assert_eq!(Molecule::new(&[(6, 1), (1, 0)]).formula().unwrap(), "C");
        assert_eq!(Molecule::new(&[(1, 0)]).formula().unwrap(), "");
    }

    #[test]
    fn formula_is_molecular_in_input_order() {
        // Molecular (not empirical) formula: C6H12O6, not CH2O.
        assert_eq!(Molecule::glucose().formula().unwrap(), "C6H12O6");
        assert_eq!(Molecule::new(&[(8, 1), (1, 2)]).formula().unwrap(), "OH2");
    }

    #[test]
    fn formula_unknown_element_with_zero_count_is_error() {
        assert!(Molecule::new(&[(200, 0)]).formula().is_err());
    }

    #[test]
    fn total_atoms_saturates() {
        let mol = Molecule::new(&[(1, u32::MAX), (8, 2)]);
        assert_eq!(mol.total_atoms(), u32::MAX);
    }

    #[test]
    fn unknown_element_is_error() {
        let mol = Molecule::new(&[(200, 1)]);
        assert!(mol.molecular_weight().is_err());
        assert!(mol.formula().is_err());
    }
}
