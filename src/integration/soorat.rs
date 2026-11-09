//! Soorat integration — visualization data structures for chemical analysis.
//!
//! Provides structured types that soorat can render: molecular structures,
//! reaction network graphs, spectroscopy data, and phase diagrams.

use serde::{Deserialize, Serialize};

// ── Molecular structure visualization ──────────────────────────────────────

/// Layout units per picometre of covalent radius, chosen so derived radii
/// track the earlier hand-tuned values (H 0.15, C 0.35, N 0.30, O 0.28,
/// S 0.45) to within ±0.025.
const RADIUS_SCALE_PER_PM: f64 = 0.0045;

/// Render radius (layout units) for elements without covalent-radius data.
const DEFAULT_RADIUS: f32 = 0.30;

/// A molecule's schematic 2D structure for ball-and-stick rendering.
///
/// Note: kimiya's [`Molecule`](crate::molecule::Molecule) stores composition
/// (which atoms, how many) but neither coordinates nor connectivity. This
/// struct provides a schematic 2D layout with schematic bonds (see
/// [`MolecularStructure::from_molecule`]). Real geometry requires a
/// molecular mechanics solver.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolecularStructure {
    /// One entry per atom: atomic number, schematic 2D position and render
    /// radius. Positions are in schematic layout units (not physical
    /// coordinates).
    pub atoms: Vec<MolAtom>,
    /// Bonds as indices into [`atoms`](Self::atoms) plus a bond order.
    pub bonds: Vec<MolBond>,
    /// Molecular formula string (e.g. "H2O"); empty if the molecule
    /// references an unknown element.
    pub formula: String,
    /// Molecular weight (g/mol); 0.0 if the molecule references an unknown
    /// element.
    pub molecular_weight: f64,
}

/// An atom in a molecular structure visualization.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MolAtom {
    /// Atomic number (1-118).
    pub atomic_number: u8,
    /// Horizontal 2D layout position (schematic layout units).
    pub x: f32,
    /// Vertical 2D layout position (schematic layout units).
    pub y: f32,
    /// Atom radius for rendering (layout units): covalent radius in pm ×
    /// 0.0045 (H 0.14, C 0.34, O 0.30), or 0.30 for elements without
    /// covalent-radius data.
    pub radius: f32,
}

/// A bond between two atoms.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MolBond {
    /// Index of the first atom.
    pub atom_a: usize,
    /// Index of the second atom.
    pub atom_b: usize,
    /// Bond order (1 = single, 2 = double, 3 = triple).
    pub order: u8,
}

impl MolecularStructure {
    /// Create from a kimiya `Molecule` with a schematic star layout.
    ///
    /// Because `Molecule` has no connectivity, the bonds are schematic: every
    /// other atom gets a single bond to one central atom, chosen by the VSEPR
    /// convention as the least electronegative non-hydrogen atom (the first
    /// such atom on ties; elements without a Pauling value rank last; the
    /// first atom if all atoms are hydrogen). Water gives two O–H bonds,
    /// CO₂ two C–O, CH₄ four C–H; molecules with several heavy atoms
    /// (e.g. ethanol) get a star, not their real bonding. The central atom
    /// is placed at the origin and the others evenly on a unit circle.
    ///
    /// This constructor cannot report errors: if the molecule references an
    /// unknown element, `formula` is empty and `molecular_weight` is 0.0.
    #[must_use]
    pub fn from_molecule(mol: &crate::molecule::Molecule) -> Self {
        let total = total_atoms(mol);
        let center = central_atom_index(mol);
        // Ligands share the unit circle; guard the single-atom case.
        let ligands = total.saturating_sub(1).max(1) as f32;

        let mut atoms = Vec::new();
        let mut ligand = 0usize;
        for atom in &mol.atoms {
            let radius = covalent_radius(atom.element_number);
            for _ in 0..atom.count {
                let (x, y) = if atoms.len() == center {
                    (0.0, 0.0)
                } else {
                    let angle = ligand as f32 * std::f32::consts::TAU / ligands;
                    ligand += 1;
                    (angle.cos(), angle.sin())
                };
                atoms.push(MolAtom {
                    atomic_number: atom.element_number,
                    x,
                    y,
                    radius,
                });
            }
        }

        // Schematic star: bond every other atom to the central atom.
        let bonds: Vec<MolBond> = (0..atoms.len())
            .filter(|&i| i != center)
            .map(|i| MolBond {
                atom_a: center,
                atom_b: i,
                order: 1,
            })
            .collect();

        Self {
            atoms,
            bonds,
            formula: mol.formula().unwrap_or_default(),
            molecular_weight: mol.molecular_weight().unwrap_or(0.0),
        }
    }
}

/// Total number of atoms, saturating instead of overflowing.
fn total_atoms(mol: &crate::molecule::Molecule) -> usize {
    mol.atoms.iter().fold(0usize, |acc, a| {
        acc.saturating_add(usize::try_from(a.count).unwrap_or(usize::MAX))
    })
}

/// Index (in the expanded one-entry-per-atom list) of the central atom: the
/// least electronegative non-hydrogen atom, first on ties, else atom 0.
fn central_atom_index(mol: &crate::molecule::Molecule) -> usize {
    let mut start = 0usize;
    let mut best: Option<(f64, usize)> = None;
    for atom in &mol.atoms {
        if atom.count > 0 && atom.element_number != 1 {
            let en = crate::element::lookup_by_number(atom.element_number)
                .and_then(|e| e.electronegativity)
                .unwrap_or(f64::INFINITY);
            if best.is_none_or(|(best_en, _)| en < best_en) {
                best = Some((en, start));
            }
        }
        start = start.saturating_add(usize::try_from(atom.count).unwrap_or(usize::MAX));
    }
    best.map_or(0, |(_, idx)| idx)
}

/// Render radius in layout units, proportional to the element's covalent
/// radius ([`crate::element::ElementProperties::covalent_radius_pm`]).
fn covalent_radius(atomic_number: u8) -> f32 {
    crate::element::lookup_properties(atomic_number).map_or(DEFAULT_RADIUS, |p| {
        (p.covalent_radius_pm * RADIUS_SCALE_PER_PM) as f32
    })
}

// ── Reaction network graph ─────────────────────────────────────────────────

/// Reaction network for node-link visualization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReactionNetworkVisualization {
    /// Species names (node labels).
    pub species: Vec<String>,
    /// Reactions as edges: `(reactant_indices, product_indices, rate_constant)`.
    pub reactions: Vec<ReactionEdge>,
}

/// A reaction edge in the network graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReactionEdge {
    /// Indices of reactant species.
    pub reactants: Vec<usize>,
    /// Indices of product species.
    pub products: Vec<usize>,
    /// Rate constant (for edge weight/thickness).
    pub rate_constant: f64,
    /// Label (e.g. "k₁ = 0.05 s⁻¹").
    pub label: String,
}

// ── Spectroscopy data ──────────────────────────────────────────────────────

/// Spectroscopy data for line/bar plot rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpectrumVisualization {
    /// Wavelength or wavenumber values (x-axis).
    pub x_values: Vec<f64>,
    /// Intensity/absorbance values (y-axis).
    pub y_values: Vec<f64>,
    /// X-axis label (e.g. "wavelength (nm)" or "wavenumber (cm⁻¹)").
    pub x_label: String,
    /// Y-axis label (e.g. "absorbance" or "transmittance").
    pub y_label: String,
    /// Peak positions (x-values where peaks occur).
    pub peaks: Vec<f64>,
}

impl SpectrumVisualization {
    /// Create from parallel x/y arrays with automatic peak detection.
    ///
    /// `peak_threshold` is passed to [`crate::spectroscopy::find_peaks`]: a
    /// fraction (0–1) of the y-value range a peak must rise above the
    /// minimum; an out-of-range threshold yields no peaks.
    #[must_use]
    pub fn from_data(
        x_values: Vec<f64>,
        y_values: Vec<f64>,
        x_label: &str,
        y_label: &str,
        peak_threshold: f64,
    ) -> Self {
        let peak_indices = crate::spectroscopy::find_peaks(&y_values, peak_threshold);
        let peaks = peak_indices
            .iter()
            .filter_map(|&i| x_values.get(i).copied())
            .collect();
        Self {
            x_values,
            y_values,
            x_label: x_label.to_string(),
            y_label: y_label.to_string(),
            peaks,
        }
    }
}

// ── Phase diagram ──────────────────────────────────────────────────────────

/// Phase diagram data for contour/boundary rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhaseDiagramVisualization {
    /// Phase boundary curves.
    pub boundaries: Vec<PhaseBoundary>,
    /// Critical point (if present): `[temperature_k, pressure_pa]`.
    pub critical_point: Option<[f64; 2]>,
    /// Triple point (if present): `[temperature_k, pressure_pa]`.
    pub triple_point: Option<[f64; 2]>,
}

/// A phase boundary curve (e.g. solid-liquid, liquid-gas).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhaseBoundary {
    /// Points along the boundary: `[temperature_k, pressure_pa]`.
    pub points: Vec<[f64; 2]>,
    /// Label for this boundary (e.g. "solid-liquid", "liquid-gas").
    pub label: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::molecule::Molecule;

    /// Bonds as sorted (atomic number, atomic number) pairs.
    fn bond_elements(viz: &MolecularStructure) -> Vec<(u8, u8)> {
        let mut pairs: Vec<(u8, u8)> = viz
            .bonds
            .iter()
            .map(|b| {
                let za = viz.atoms[b.atom_a].atomic_number;
                let zb = viz.atoms[b.atom_b].atomic_number;
                (za.min(zb), za.max(zb))
            })
            .collect();
        pairs.sort_unstable();
        pairs
    }

    #[test]
    fn molecular_structure_water() {
        // H₂O: 2 hydrogen + 1 oxygen → two O–H bonds (was H–H and H–O)
        let viz = MolecularStructure::from_molecule(&Molecule::water());
        assert_eq!(viz.atoms.len(), 3); // 2H + 1O
        assert_eq!(bond_elements(&viz), vec![(1, 8), (1, 8)]);
        assert!(viz.formula.contains('H'));
        assert!(viz.formula.contains('O'));
        assert!(viz.molecular_weight > 17.0 && viz.molecular_weight < 19.0);
        // Central O at the origin.
        assert_eq!((viz.atoms[2].x, viz.atoms[2].y), (0.0, 0.0));
    }

    #[test]
    fn molecular_structure_co2() {
        let viz = MolecularStructure::from_molecule(&Molecule::carbon_dioxide());
        assert_eq!(bond_elements(&viz), vec![(6, 8), (6, 8)]);
    }

    #[test]
    fn molecular_structure_ch4() {
        let viz = MolecularStructure::from_molecule(&Molecule::methane());
        assert_eq!(bond_elements(&viz), vec![(1, 6); 4]);
        // Ligands on the unit circle around the central C.
        for atom in viz.atoms.iter().filter(|a| a.atomic_number == 1) {
            assert!(((atom.x * atom.x + atom.y * atom.y).sqrt() - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn molecular_structure_central_atom_rules() {
        // All hydrogen → first atom is central.
        let viz = MolecularStructure::from_molecule(&Molecule::new(&[(1, 2)]));
        assert_eq!(bond_elements(&viz), vec![(1, 1)]);
        assert_eq!(viz.bonds[0].atom_a, 0);
        // Ties → first: two C groups, the first C (expanded index 1) wins.
        let viz = MolecularStructure::from_molecule(&Molecule::new(&[(1, 1), (6, 1), (6, 1)]));
        assert!(viz.bonds.iter().all(|b| b.atom_a == 1));
        // SF₆: S (2.58) is less electronegative than F (3.98).
        let viz = MolecularStructure::from_molecule(&Molecule::new(&[(9, 6), (16, 1)]));
        assert_eq!(bond_elements(&viz), vec![(9, 16); 6]);
    }

    #[test]
    fn molecular_structure_single_atom() {
        let mol = Molecule::new(&[(2, 1)]); // He
        let viz = MolecularStructure::from_molecule(&mol);
        assert_eq!(viz.atoms.len(), 1);
        assert!(viz.bonds.is_empty());
    }

    #[test]
    fn molecular_structure_empty_and_invalid() {
        let viz = MolecularStructure::from_molecule(&Molecule::new(&[]));
        assert!(viz.atoms.is_empty() && viz.bonds.is_empty());
        // Unknown element: documented defaults instead of an error.
        let viz = MolecularStructure::from_molecule(&Molecule::new(&[(0, 1), (8, 1)]));
        assert!(viz.formula.is_empty());
        assert_eq!(viz.molecular_weight, 0.0);
        assert_eq!(viz.atoms.len(), 2);
    }

    #[test]
    fn total_atoms_does_not_overflow() {
        // The old inline `sum::<u32>()` overflowed (panic in debug builds).
        let mol = Molecule::new(&[(1, u32::MAX), (8, u32::MAX)]);
        assert_eq!(total_atoms(&mol), 2 * u32::MAX as usize);
    }

    #[test]
    fn covalent_radius_known() {
        // Derived from covalent radii; close to the old hand-tuned table.
        for (z, old) in [(1, 0.15), (6, 0.35), (7, 0.30), (8, 0.28), (16, 0.45)] {
            assert!(
                (covalent_radius(z) - old).abs() < 0.025,
                "Z={z}: {}",
                covalent_radius(z)
            );
        }
        // H (31 pm) < O (66) < N (71) < C (76) < S (105)
        let r: Vec<f32> = [1, 8, 7, 6, 16].map(covalent_radius).to_vec();
        assert!(r.windows(2).all(|w| w[0] < w[1]), "{r:?}");
        // Fe (132 pm) now comes from the element table instead of the default.
        assert!((covalent_radius(26) - 0.594).abs() < 1e-6);
        // No covalent-radius data → fallback.
        assert_eq!(covalent_radius(118), DEFAULT_RADIUS);
    }

    #[test]
    fn reaction_network_manual() {
        let net = ReactionNetworkVisualization {
            species: vec!["A".into(), "B".into(), "C".into()],
            reactions: vec![ReactionEdge {
                reactants: vec![0],
                products: vec![1, 2],
                rate_constant: 0.05,
                label: "k₁".into(),
            }],
        };
        assert_eq!(net.species.len(), 3);
        assert_eq!(net.reactions.len(), 1);
    }

    #[test]
    fn spectrum_from_data() {
        let x = vec![400.0, 450.0, 500.0, 550.0, 600.0];
        let y = vec![0.1, 0.5, 0.9, 0.3, 0.1];
        let spec = SpectrumVisualization::from_data(x, y, "wavelength (nm)", "absorbance", 0.4);
        assert_eq!(spec.x_values.len(), 5);
        assert!(!spec.peaks.is_empty()); // should find peak near 500nm
    }

    #[test]
    fn spectrum_empty() {
        let spec = SpectrumVisualization::from_data(vec![], vec![], "nm", "A", 0.5);
        assert!(spec.peaks.is_empty());
    }

    #[test]
    fn phase_diagram_serializes() {
        let pd = PhaseDiagramVisualization {
            boundaries: vec![PhaseBoundary {
                points: vec![[273.15, 611.657], [373.15, 101325.0]],
                label: "liquid-gas".into(),
            }],
            critical_point: Some([647.096, 22064000.0]),
            triple_point: Some([273.16, 611.657]),
        };
        let json = serde_json::to_string(&pd);
        assert!(json.is_ok());
    }
}
