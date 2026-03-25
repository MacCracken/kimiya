use serde::{Deserialize, Serialize};

/// Category of an element in the periodic table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum ElementCategory {
    AlkaliMetal,
    AlkalineEarth,
    TransitionMetal,
    PostTransitionMetal,
    Metalloid,
    Nonmetal,
    Halogen,
    NobleGas,
}

/// A chemical element.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Element {
    pub atomic_number: u8,
    pub symbol: &'static str,
    pub name: &'static str,
    pub atomic_mass: f64,
    pub electronegativity: Option<f64>,
    pub category: ElementCategory,
}

/// Avogadro's number (mol⁻¹).
pub const AVOGADRO: f64 = 6.02214076e23;

/// Universal gas constant (J/(mol·K)).
pub const GAS_CONSTANT: f64 = 8.314462618;

/// Built-in element data for Z=1 through Z=36.
pub static ELEMENTS: &[Element] = &[
    Element { atomic_number: 1, symbol: "H", name: "Hydrogen", atomic_mass: 1.008, electronegativity: Some(2.20), category: ElementCategory::Nonmetal },
    Element { atomic_number: 2, symbol: "He", name: "Helium", atomic_mass: 4.003, electronegativity: None, category: ElementCategory::NobleGas },
    Element { atomic_number: 3, symbol: "Li", name: "Lithium", atomic_mass: 6.941, electronegativity: Some(0.98), category: ElementCategory::AlkaliMetal },
    Element { atomic_number: 4, symbol: "Be", name: "Beryllium", atomic_mass: 9.012, electronegativity: Some(1.57), category: ElementCategory::AlkalineEarth },
    Element { atomic_number: 5, symbol: "B", name: "Boron", atomic_mass: 10.81, electronegativity: Some(2.04), category: ElementCategory::Metalloid },
    Element { atomic_number: 6, symbol: "C", name: "Carbon", atomic_mass: 12.011, electronegativity: Some(2.55), category: ElementCategory::Nonmetal },
    Element { atomic_number: 7, symbol: "N", name: "Nitrogen", atomic_mass: 14.007, electronegativity: Some(3.04), category: ElementCategory::Nonmetal },
    Element { atomic_number: 8, symbol: "O", name: "Oxygen", atomic_mass: 15.999, electronegativity: Some(3.44), category: ElementCategory::Nonmetal },
    Element { atomic_number: 9, symbol: "F", name: "Fluorine", atomic_mass: 18.998, electronegativity: Some(3.98), category: ElementCategory::Halogen },
    Element { atomic_number: 10, symbol: "Ne", name: "Neon", atomic_mass: 20.180, electronegativity: None, category: ElementCategory::NobleGas },
    Element { atomic_number: 11, symbol: "Na", name: "Sodium", atomic_mass: 22.990, electronegativity: Some(0.93), category: ElementCategory::AlkaliMetal },
    Element { atomic_number: 12, symbol: "Mg", name: "Magnesium", atomic_mass: 24.305, electronegativity: Some(1.31), category: ElementCategory::AlkalineEarth },
    Element { atomic_number: 13, symbol: "Al", name: "Aluminum", atomic_mass: 26.982, electronegativity: Some(1.61), category: ElementCategory::PostTransitionMetal },
    Element { atomic_number: 14, symbol: "Si", name: "Silicon", atomic_mass: 28.086, electronegativity: Some(1.90), category: ElementCategory::Metalloid },
    Element { atomic_number: 15, symbol: "P", name: "Phosphorus", atomic_mass: 30.974, electronegativity: Some(2.19), category: ElementCategory::Nonmetal },
    Element { atomic_number: 16, symbol: "S", name: "Sulfur", atomic_mass: 32.065, electronegativity: Some(2.58), category: ElementCategory::Nonmetal },
    Element { atomic_number: 17, symbol: "Cl", name: "Chlorine", atomic_mass: 35.453, electronegativity: Some(3.16), category: ElementCategory::Halogen },
    Element { atomic_number: 18, symbol: "Ar", name: "Argon", atomic_mass: 39.948, electronegativity: None, category: ElementCategory::NobleGas },
    Element { atomic_number: 19, symbol: "K", name: "Potassium", atomic_mass: 39.098, electronegativity: Some(0.82), category: ElementCategory::AlkaliMetal },
    Element { atomic_number: 20, symbol: "Ca", name: "Calcium", atomic_mass: 40.078, electronegativity: Some(1.00), category: ElementCategory::AlkalineEarth },
    Element { atomic_number: 21, symbol: "Sc", name: "Scandium", atomic_mass: 44.956, electronegativity: Some(1.36), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 22, symbol: "Ti", name: "Titanium", atomic_mass: 47.867, electronegativity: Some(1.54), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 23, symbol: "V", name: "Vanadium", atomic_mass: 50.942, electronegativity: Some(1.63), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 24, symbol: "Cr", name: "Chromium", atomic_mass: 51.996, electronegativity: Some(1.66), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 25, symbol: "Mn", name: "Manganese", atomic_mass: 54.938, electronegativity: Some(1.55), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 26, symbol: "Fe", name: "Iron", atomic_mass: 55.845, electronegativity: Some(1.83), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 27, symbol: "Co", name: "Cobalt", atomic_mass: 58.933, electronegativity: Some(1.88), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 28, symbol: "Ni", name: "Nickel", atomic_mass: 58.693, electronegativity: Some(1.91), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 29, symbol: "Cu", name: "Copper", atomic_mass: 63.546, electronegativity: Some(1.90), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 30, symbol: "Zn", name: "Zinc", atomic_mass: 65.380, electronegativity: Some(1.65), category: ElementCategory::TransitionMetal },
    Element { atomic_number: 31, symbol: "Ga", name: "Gallium", atomic_mass: 69.723, electronegativity: Some(1.81), category: ElementCategory::PostTransitionMetal },
    Element { atomic_number: 32, symbol: "Ge", name: "Germanium", atomic_mass: 72.630, electronegativity: Some(2.01), category: ElementCategory::Metalloid },
    Element { atomic_number: 33, symbol: "As", name: "Arsenic", atomic_mass: 74.922, electronegativity: Some(2.18), category: ElementCategory::Metalloid },
    Element { atomic_number: 34, symbol: "Se", name: "Selenium", atomic_mass: 78.971, electronegativity: Some(2.55), category: ElementCategory::Nonmetal },
    Element { atomic_number: 35, symbol: "Br", name: "Bromine", atomic_mass: 79.904, electronegativity: Some(2.96), category: ElementCategory::Halogen },
    Element { atomic_number: 36, symbol: "Kr", name: "Krypton", atomic_mass: 83.798, electronegativity: Some(3.00), category: ElementCategory::NobleGas },
];

/// Look up an element by atomic number (1-based).
#[must_use]
pub fn lookup_by_number(atomic_number: u8) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.atomic_number == atomic_number)
}

/// Look up an element by symbol (case-sensitive).
#[must_use]
pub fn lookup_by_symbol(symbol: &str) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.symbol == symbol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrogen_data() {
        let h = lookup_by_number(1).unwrap();
        assert_eq!(h.symbol, "H");
        assert!((h.atomic_mass - 1.008).abs() < 0.001);
    }

    #[test]
    fn iron_data() {
        let fe = lookup_by_symbol("Fe").unwrap();
        assert_eq!(fe.atomic_number, 26);
        assert!((fe.atomic_mass - 55.845).abs() < 0.001);
    }

    #[test]
    fn all_36_elements() {
        assert_eq!(ELEMENTS.len(), 36);
        assert_eq!(ELEMENTS[0].atomic_number, 1);
        assert_eq!(ELEMENTS[35].atomic_number, 36);
    }

    #[test]
    fn noble_gases_light_no_electronegativity() {
        // He, Ne, Ar have no Pauling electronegativity. Kr and Xe do (they form compounds).
        for symbol in &["He", "Ne", "Ar"] {
            let e = lookup_by_symbol(symbol).unwrap();
            assert!(e.electronegativity.is_none(), "{} should have no electronegativity", e.name);
        }
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_by_number(99).is_none());
        assert!(lookup_by_symbol("Xx").is_none());
    }

    #[test]
    fn avogadro_constant() {
        assert!((AVOGADRO - 6.022e23).abs() < 1e20);
    }

    #[test]
    fn gas_constant_value() {
        assert!((GAS_CONSTANT - 8.314).abs() < 0.001);
    }
}
