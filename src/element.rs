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
    Lanthanide,
    Actinide,
}

/// A chemical element.
///
/// Note: `Element` implements `Serialize` but not `Deserialize` because the
/// built-in element data uses `&'static str` for symbol and name. Use
/// [`lookup_by_number`] or [`lookup_by_symbol`] to obtain element references.
#[derive(Debug, Clone, Serialize)]
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

/// Built-in element data for Z=1 through Z=118 (full periodic table).
pub static ELEMENTS: &[Element] = &[
    Element {
        atomic_number: 1,
        symbol: "H",
        name: "Hydrogen",
        atomic_mass: 1.008,
        electronegativity: Some(2.20),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 2,
        symbol: "He",
        name: "Helium",
        atomic_mass: 4.003,
        electronegativity: None,
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 3,
        symbol: "Li",
        name: "Lithium",
        atomic_mass: 6.941,
        electronegativity: Some(0.98),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 4,
        symbol: "Be",
        name: "Beryllium",
        atomic_mass: 9.012,
        electronegativity: Some(1.57),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 5,
        symbol: "B",
        name: "Boron",
        atomic_mass: 10.81,
        electronegativity: Some(2.04),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 6,
        symbol: "C",
        name: "Carbon",
        atomic_mass: 12.011,
        electronegativity: Some(2.55),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 7,
        symbol: "N",
        name: "Nitrogen",
        atomic_mass: 14.007,
        electronegativity: Some(3.04),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 8,
        symbol: "O",
        name: "Oxygen",
        atomic_mass: 15.999,
        electronegativity: Some(3.44),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 9,
        symbol: "F",
        name: "Fluorine",
        atomic_mass: 18.998,
        electronegativity: Some(3.98),
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 10,
        symbol: "Ne",
        name: "Neon",
        atomic_mass: 20.180,
        electronegativity: None,
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 11,
        symbol: "Na",
        name: "Sodium",
        atomic_mass: 22.990,
        electronegativity: Some(0.93),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 12,
        symbol: "Mg",
        name: "Magnesium",
        atomic_mass: 24.305,
        electronegativity: Some(1.31),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 13,
        symbol: "Al",
        name: "Aluminum",
        atomic_mass: 26.982,
        electronegativity: Some(1.61),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 14,
        symbol: "Si",
        name: "Silicon",
        atomic_mass: 28.086,
        electronegativity: Some(1.90),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 15,
        symbol: "P",
        name: "Phosphorus",
        atomic_mass: 30.974,
        electronegativity: Some(2.19),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 16,
        symbol: "S",
        name: "Sulfur",
        atomic_mass: 32.065,
        electronegativity: Some(2.58),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 17,
        symbol: "Cl",
        name: "Chlorine",
        atomic_mass: 35.453,
        electronegativity: Some(3.16),
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 18,
        symbol: "Ar",
        name: "Argon",
        atomic_mass: 39.948,
        electronegativity: None,
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 19,
        symbol: "K",
        name: "Potassium",
        atomic_mass: 39.098,
        electronegativity: Some(0.82),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 20,
        symbol: "Ca",
        name: "Calcium",
        atomic_mass: 40.078,
        electronegativity: Some(1.00),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 21,
        symbol: "Sc",
        name: "Scandium",
        atomic_mass: 44.956,
        electronegativity: Some(1.36),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 22,
        symbol: "Ti",
        name: "Titanium",
        atomic_mass: 47.867,
        electronegativity: Some(1.54),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 23,
        symbol: "V",
        name: "Vanadium",
        atomic_mass: 50.942,
        electronegativity: Some(1.63),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 24,
        symbol: "Cr",
        name: "Chromium",
        atomic_mass: 51.996,
        electronegativity: Some(1.66),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 25,
        symbol: "Mn",
        name: "Manganese",
        atomic_mass: 54.938,
        electronegativity: Some(1.55),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 26,
        symbol: "Fe",
        name: "Iron",
        atomic_mass: 55.845,
        electronegativity: Some(1.83),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 27,
        symbol: "Co",
        name: "Cobalt",
        atomic_mass: 58.933,
        electronegativity: Some(1.88),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 28,
        symbol: "Ni",
        name: "Nickel",
        atomic_mass: 58.693,
        electronegativity: Some(1.91),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 29,
        symbol: "Cu",
        name: "Copper",
        atomic_mass: 63.546,
        electronegativity: Some(1.90),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 30,
        symbol: "Zn",
        name: "Zinc",
        atomic_mass: 65.380,
        electronegativity: Some(1.65),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 31,
        symbol: "Ga",
        name: "Gallium",
        atomic_mass: 69.723,
        electronegativity: Some(1.81),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 32,
        symbol: "Ge",
        name: "Germanium",
        atomic_mass: 72.630,
        electronegativity: Some(2.01),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 33,
        symbol: "As",
        name: "Arsenic",
        atomic_mass: 74.922,
        electronegativity: Some(2.18),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 34,
        symbol: "Se",
        name: "Selenium",
        atomic_mass: 78.971,
        electronegativity: Some(2.55),
        category: ElementCategory::Nonmetal,
    },
    Element {
        atomic_number: 35,
        symbol: "Br",
        name: "Bromine",
        atomic_mass: 79.904,
        electronegativity: Some(2.96),
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 36,
        symbol: "Kr",
        name: "Krypton",
        atomic_mass: 83.798,
        electronegativity: Some(3.00),
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 37,
        symbol: "Rb",
        name: "Rubidium",
        atomic_mass: 85.468,
        electronegativity: Some(0.82),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 38,
        symbol: "Sr",
        name: "Strontium",
        atomic_mass: 87.62,
        electronegativity: Some(0.95),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 39,
        symbol: "Y",
        name: "Yttrium",
        atomic_mass: 88.906,
        electronegativity: Some(1.22),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 40,
        symbol: "Zr",
        name: "Zirconium",
        atomic_mass: 91.224,
        electronegativity: Some(1.33),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 41,
        symbol: "Nb",
        name: "Niobium",
        atomic_mass: 92.906,
        electronegativity: Some(1.6),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 42,
        symbol: "Mo",
        name: "Molybdenum",
        atomic_mass: 95.95,
        electronegativity: Some(2.16),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 43,
        symbol: "Tc",
        name: "Technetium",
        atomic_mass: 97.0,
        electronegativity: Some(1.9),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 44,
        symbol: "Ru",
        name: "Ruthenium",
        atomic_mass: 101.07,
        electronegativity: Some(2.2),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 45,
        symbol: "Rh",
        name: "Rhodium",
        atomic_mass: 102.906,
        electronegativity: Some(2.28),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 46,
        symbol: "Pd",
        name: "Palladium",
        atomic_mass: 106.42,
        electronegativity: Some(2.20),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 47,
        symbol: "Ag",
        name: "Silver",
        atomic_mass: 107.868,
        electronegativity: Some(1.93),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 48,
        symbol: "Cd",
        name: "Cadmium",
        atomic_mass: 112.414,
        electronegativity: Some(1.69),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 49,
        symbol: "In",
        name: "Indium",
        atomic_mass: 114.818,
        electronegativity: Some(1.78),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 50,
        symbol: "Sn",
        name: "Tin",
        atomic_mass: 118.710,
        electronegativity: Some(1.96),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 51,
        symbol: "Sb",
        name: "Antimony",
        atomic_mass: 121.760,
        electronegativity: Some(2.05),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 52,
        symbol: "Te",
        name: "Tellurium",
        atomic_mass: 127.60,
        electronegativity: Some(2.1),
        category: ElementCategory::Metalloid,
    },
    Element {
        atomic_number: 53,
        symbol: "I",
        name: "Iodine",
        atomic_mass: 126.904,
        electronegativity: Some(2.66),
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 54,
        symbol: "Xe",
        name: "Xenon",
        atomic_mass: 131.293,
        electronegativity: Some(2.60),
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 55,
        symbol: "Cs",
        name: "Cesium",
        atomic_mass: 132.905,
        electronegativity: Some(0.79),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 56,
        symbol: "Ba",
        name: "Barium",
        atomic_mass: 137.327,
        electronegativity: Some(0.89),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 57,
        symbol: "La",
        name: "Lanthanum",
        atomic_mass: 138.905,
        electronegativity: Some(1.10),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 58,
        symbol: "Ce",
        name: "Cerium",
        atomic_mass: 140.116,
        electronegativity: Some(1.12),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 59,
        symbol: "Pr",
        name: "Praseodymium",
        atomic_mass: 140.908,
        electronegativity: Some(1.13),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 60,
        symbol: "Nd",
        name: "Neodymium",
        atomic_mass: 144.242,
        electronegativity: Some(1.14),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 61,
        symbol: "Pm",
        name: "Promethium",
        atomic_mass: 145.0,
        electronegativity: Some(1.13),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 62,
        symbol: "Sm",
        name: "Samarium",
        atomic_mass: 150.36,
        electronegativity: Some(1.17),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 63,
        symbol: "Eu",
        name: "Europium",
        atomic_mass: 151.964,
        electronegativity: Some(1.20),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 64,
        symbol: "Gd",
        name: "Gadolinium",
        atomic_mass: 157.25,
        electronegativity: Some(1.20),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 65,
        symbol: "Tb",
        name: "Terbium",
        atomic_mass: 158.925,
        electronegativity: Some(1.20),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 66,
        symbol: "Dy",
        name: "Dysprosium",
        atomic_mass: 162.500,
        electronegativity: Some(1.22),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 67,
        symbol: "Ho",
        name: "Holmium",
        atomic_mass: 164.930,
        electronegativity: Some(1.23),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 68,
        symbol: "Er",
        name: "Erbium",
        atomic_mass: 167.259,
        electronegativity: Some(1.24),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 69,
        symbol: "Tm",
        name: "Thulium",
        atomic_mass: 168.934,
        electronegativity: Some(1.25),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 70,
        symbol: "Yb",
        name: "Ytterbium",
        atomic_mass: 173.045,
        electronegativity: Some(1.10),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 71,
        symbol: "Lu",
        name: "Lutetium",
        atomic_mass: 174.967,
        electronegativity: Some(1.27),
        category: ElementCategory::Lanthanide,
    },
    Element {
        atomic_number: 72,
        symbol: "Hf",
        name: "Hafnium",
        atomic_mass: 178.486,
        electronegativity: Some(1.3),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 73,
        symbol: "Ta",
        name: "Tantalum",
        atomic_mass: 180.948,
        electronegativity: Some(1.5),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 74,
        symbol: "W",
        name: "Tungsten",
        atomic_mass: 183.84,
        electronegativity: Some(2.36),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 75,
        symbol: "Re",
        name: "Rhenium",
        atomic_mass: 186.207,
        electronegativity: Some(1.9),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 76,
        symbol: "Os",
        name: "Osmium",
        atomic_mass: 190.23,
        electronegativity: Some(2.2),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 77,
        symbol: "Ir",
        name: "Iridium",
        atomic_mass: 192.217,
        electronegativity: Some(2.20),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 78,
        symbol: "Pt",
        name: "Platinum",
        atomic_mass: 195.084,
        electronegativity: Some(2.28),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 79,
        symbol: "Au",
        name: "Gold",
        atomic_mass: 196.967,
        electronegativity: Some(2.54),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 80,
        symbol: "Hg",
        name: "Mercury",
        atomic_mass: 200.592,
        electronegativity: Some(2.00),
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 81,
        symbol: "Tl",
        name: "Thallium",
        atomic_mass: 204.38,
        electronegativity: Some(1.62),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 82,
        symbol: "Pb",
        name: "Lead",
        atomic_mass: 207.2,
        electronegativity: Some(1.87),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 83,
        symbol: "Bi",
        name: "Bismuth",
        atomic_mass: 208.980,
        electronegativity: Some(2.02),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 84,
        symbol: "Po",
        name: "Polonium",
        atomic_mass: 209.0,
        electronegativity: Some(2.0),
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 85,
        symbol: "At",
        name: "Astatine",
        atomic_mass: 210.0,
        electronegativity: Some(2.2),
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 86,
        symbol: "Rn",
        name: "Radon",
        atomic_mass: 222.0,
        electronegativity: Some(2.2),
        category: ElementCategory::NobleGas,
    },
    Element {
        atomic_number: 87,
        symbol: "Fr",
        name: "Francium",
        atomic_mass: 223.0,
        electronegativity: Some(0.7),
        category: ElementCategory::AlkaliMetal,
    },
    Element {
        atomic_number: 88,
        symbol: "Ra",
        name: "Radium",
        atomic_mass: 226.0,
        electronegativity: Some(0.9),
        category: ElementCategory::AlkalineEarth,
    },
    Element {
        atomic_number: 89,
        symbol: "Ac",
        name: "Actinium",
        atomic_mass: 227.0,
        electronegativity: Some(1.1),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 90,
        symbol: "Th",
        name: "Thorium",
        atomic_mass: 232.038,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 91,
        symbol: "Pa",
        name: "Protactinium",
        atomic_mass: 231.036,
        electronegativity: Some(1.5),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 92,
        symbol: "U",
        name: "Uranium",
        atomic_mass: 238.029,
        electronegativity: Some(1.38),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 93,
        symbol: "Np",
        name: "Neptunium",
        atomic_mass: 237.0,
        electronegativity: Some(1.36),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 94,
        symbol: "Pu",
        name: "Plutonium",
        atomic_mass: 244.0,
        electronegativity: Some(1.28),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 95,
        symbol: "Am",
        name: "Americium",
        atomic_mass: 243.0,
        electronegativity: Some(1.13),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 96,
        symbol: "Cm",
        name: "Curium",
        atomic_mass: 247.0,
        electronegativity: Some(1.28),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 97,
        symbol: "Bk",
        name: "Berkelium",
        atomic_mass: 247.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 98,
        symbol: "Cf",
        name: "Californium",
        atomic_mass: 251.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 99,
        symbol: "Es",
        name: "Einsteinium",
        atomic_mass: 252.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 100,
        symbol: "Fm",
        name: "Fermium",
        atomic_mass: 257.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 101,
        symbol: "Md",
        name: "Mendelevium",
        atomic_mass: 258.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 102,
        symbol: "No",
        name: "Nobelium",
        atomic_mass: 259.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 103,
        symbol: "Lr",
        name: "Lawrencium",
        atomic_mass: 266.0,
        electronegativity: Some(1.3),
        category: ElementCategory::Actinide,
    },
    Element {
        atomic_number: 104,
        symbol: "Rf",
        name: "Rutherfordium",
        atomic_mass: 267.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 105,
        symbol: "Db",
        name: "Dubnium",
        atomic_mass: 268.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 106,
        symbol: "Sg",
        name: "Seaborgium",
        atomic_mass: 269.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 107,
        symbol: "Bh",
        name: "Bohrium",
        atomic_mass: 270.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 108,
        symbol: "Hs",
        name: "Hassium",
        atomic_mass: 269.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 109,
        symbol: "Mt",
        name: "Meitnerium",
        atomic_mass: 278.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 110,
        symbol: "Ds",
        name: "Darmstadtium",
        atomic_mass: 281.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 111,
        symbol: "Rg",
        name: "Roentgenium",
        atomic_mass: 282.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 112,
        symbol: "Cn",
        name: "Copernicium",
        atomic_mass: 285.0,
        electronegativity: None,
        category: ElementCategory::TransitionMetal,
    },
    Element {
        atomic_number: 113,
        symbol: "Nh",
        name: "Nihonium",
        atomic_mass: 286.0,
        electronegativity: None,
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 114,
        symbol: "Fl",
        name: "Flerovium",
        atomic_mass: 289.0,
        electronegativity: None,
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 115,
        symbol: "Mc",
        name: "Moscovium",
        atomic_mass: 290.0,
        electronegativity: None,
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 116,
        symbol: "Lv",
        name: "Livermorium",
        atomic_mass: 293.0,
        electronegativity: None,
        category: ElementCategory::PostTransitionMetal,
    },
    Element {
        atomic_number: 117,
        symbol: "Ts",
        name: "Tennessine",
        atomic_mass: 294.0,
        electronegativity: None,
        category: ElementCategory::Halogen,
    },
    Element {
        atomic_number: 118,
        symbol: "Og",
        name: "Oganesson",
        atomic_mass: 294.0,
        electronegativity: None,
        category: ElementCategory::NobleGas,
    },
];

/// Look up an element by atomic number (1-based).
///
/// O(1) for elements within the built-in range (Z=1..118), O(n) fallback otherwise.
#[must_use]
#[inline]
pub fn lookup_by_number(atomic_number: u8) -> Option<&'static Element> {
    let idx = atomic_number.checked_sub(1)? as usize;
    ELEMENTS
        .get(idx)
        .filter(|e| e.atomic_number == atomic_number)
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
    fn all_118_elements() {
        assert_eq!(ELEMENTS.len(), 118);
        assert_eq!(ELEMENTS[0].atomic_number, 1);
        assert_eq!(ELEMENTS[117].atomic_number, 118);
    }

    #[test]
    fn noble_gases_light_no_electronegativity() {
        // He, Ne, Ar have no Pauling electronegativity. Kr and Xe do (they form compounds).
        for symbol in &["He", "Ne", "Ar"] {
            let e = lookup_by_symbol(symbol).unwrap();
            assert!(
                e.electronegativity.is_none(),
                "{} should have no electronegativity",
                e.name
            );
        }
    }

    #[test]
    fn lookup_nonexistent() {
        assert!(lookup_by_number(200).is_none());
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

    #[test]
    fn lookup_zero_returns_none() {
        assert!(lookup_by_number(0).is_none());
    }

    #[test]
    fn lookup_by_number_o1_correctness() {
        // Verify O(1) index lookup matches for all 36 elements
        for element in ELEMENTS.iter() {
            let found = lookup_by_number(element.atomic_number).unwrap();
            assert_eq!(found.symbol, element.symbol);
        }
    }

    #[test]
    fn elements_contiguous_and_ordered() {
        for (i, element) in ELEMENTS.iter().enumerate() {
            assert_eq!(
                element.atomic_number as usize,
                i + 1,
                "element at index {i} has wrong atomic number"
            );
        }
    }

    #[test]
    fn electronegativity_bounds() {
        // Pauling electronegativity ranges from ~0.7 (Cs, Fr) to 3.98 (F)
        for element in ELEMENTS.iter() {
            if let Some(en) = element.electronegativity {
                assert!(
                    (0.5..=4.0).contains(&en),
                    "{} has out-of-range electronegativity: {en}",
                    element.symbol
                );
            }
        }
    }

    #[test]
    fn atomic_mass_positive() {
        for element in ELEMENTS.iter() {
            assert!(
                element.atomic_mass > 0.0,
                "{} has non-positive atomic mass",
                element.symbol
            );
        }
    }
}
