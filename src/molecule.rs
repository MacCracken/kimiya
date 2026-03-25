use serde::{Deserialize, Serialize};
use crate::element::{ELEMENTS, lookup_by_number};

/// An atom count in a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element_number: u8,
    pub count: u32,
}

/// Bond type between atoms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Bond {
    Single,
    Double,
    Triple,
}

/// A molecule defined by its constituent atoms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    /// Create a molecule from a list of (atomic_number, count) pairs.
    #[must_use]
    pub fn new(atoms: &[(u8, u32)]) -> Self {
        Self {
            atoms: atoms.iter().map(|&(n, c)| Atom { element_number: n, count: c }).collect(),
        }
    }

    /// Water: H₂O
    #[must_use]
    pub fn water() -> Self { Self::new(&[(1, 2), (8, 1)]) }

    /// Carbon dioxide: CO₂
    #[must_use]
    pub fn carbon_dioxide() -> Self { Self::new(&[(6, 1), (8, 2)]) }

    /// Methane: CH₄
    #[must_use]
    pub fn methane() -> Self { Self::new(&[(6, 1), (1, 4)]) }

    /// Glucose: C₆H₁₂O₆
    #[must_use]
    pub fn glucose() -> Self { Self::new(&[(6, 6), (1, 12), (8, 6)]) }

    /// Molecular weight in g/mol (sum of atomic masses × count).
    #[must_use]
    pub fn molecular_weight(&self) -> f64 {
        self.atoms.iter().map(|a| {
            lookup_by_number(a.element_number)
                .map(|e| e.atomic_mass * a.count as f64)
                .unwrap_or(0.0)
        }).sum()
    }

    /// Empirical formula string (e.g., "H2O", "CO2").
    #[must_use]
    pub fn formula(&self) -> String {
        let mut s = String::new();
        for atom in &self.atoms {
            if let Some(e) = lookup_by_number(atom.element_number) {
                s.push_str(e.symbol);
                if atom.count > 1 {
                    s.push_str(&atom.count.to_string());
                }
            }
        }
        s
    }

    /// Total number of atoms in the molecule.
    #[must_use]
    pub fn total_atoms(&self) -> u32 {
        self.atoms.iter().map(|a| a.count).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_molecular_weight() {
        let mw = Molecule::water().molecular_weight();
        assert!((mw - 18.015).abs() < 0.01, "H2O should be ~18.015, got {mw}");
    }

    #[test]
    fn co2_molecular_weight() {
        let mw = Molecule::carbon_dioxide().molecular_weight();
        assert!((mw - 44.009).abs() < 0.02, "CO2 should be ~44.01, got {mw}");
    }

    #[test]
    fn water_formula() {
        assert_eq!(Molecule::water().formula(), "H2O");
    }

    #[test]
    fn co2_formula() {
        assert_eq!(Molecule::carbon_dioxide().formula(), "CO2");
    }

    #[test]
    fn glucose_molecular_weight() {
        let mw = Molecule::glucose().molecular_weight();
        assert!((mw - 180.156).abs() < 0.1, "glucose should be ~180.16, got {mw}");
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
        assert!((back.molecular_weight() - m.molecular_weight()).abs() < 0.001);
    }
}
