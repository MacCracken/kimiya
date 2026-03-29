//! Soorat integration — visualization data structures for chemical analysis.
//!
//! Provides structured types that soorat can render: molecular structures,
//! reaction network graphs, spectroscopy data, and phase diagrams.

use serde::{Deserialize, Serialize};

// ── Molecular structure visualization ──────────────────────────────────────

/// A molecule's 3D structure for ball-and-stick rendering.
///
/// Note: kimiya stores composition (which atoms, how many) but not
/// 3D coordinates. This struct provides a 2D layout based on composition
/// for schematic rendering. Full 3D coordinates require a molecular
/// mechanics solver.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolecularStructure {
    /// Atoms in the molecule: `(atomic_number, symbol, x, y)`.
    /// Positions are in schematic layout units (not physical coordinates).
    pub atoms: Vec<MolAtom>,
    /// Bonds connecting atoms: `(atom_idx_a, atom_idx_b, bond_order)`.
    pub bonds: Vec<MolBond>,
    /// Molecular formula string (e.g. "H2O").
    pub formula: String,
    /// Molecular weight (g/mol).
    pub molecular_weight: f64,
}

/// An atom in a molecular structure visualization.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MolAtom {
    /// Atomic number (1-118).
    pub atomic_number: u8,
    /// 2D layout position (schematic).
    pub x: f32,
    /// 2D layout position (schematic).
    pub y: f32,
    /// Atom radius for rendering (proportional to covalent radius).
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
    /// Create from a kimiya `Molecule` with a simple circular layout.
    #[must_use]
    pub fn from_molecule(mol: &crate::molecule::Molecule) -> Self {
        let mut atoms = Vec::new();
        let mut idx = 0usize;

        // Expand each atom type into individual atoms
        for atom in &mol.atoms {
            for _ in 0..atom.count {
                let angle = idx as f32 * std::f32::consts::TAU
                    / mol.atoms.iter().map(|a| a.count).sum::<u32>().max(1) as f32;
                let r = 1.0;
                atoms.push(MolAtom {
                    atomic_number: atom.element_number,
                    x: r * angle.cos(),
                    y: r * angle.sin(),
                    radius: covalent_radius(atom.element_number),
                });
                idx += 1;
            }
        }

        // Simple bonds: connect sequential atoms
        let bonds: Vec<MolBond> = if atoms.len() > 1 {
            (0..atoms.len() - 1)
                .map(|i| MolBond {
                    atom_a: i,
                    atom_b: i + 1,
                    order: 1,
                })
                .collect()
        } else {
            Vec::new()
        };

        Self {
            atoms,
            bonds,
            formula: mol.formula().unwrap_or_default(),
            molecular_weight: mol.molecular_weight().unwrap_or(0.0),
        }
    }
}

/// Approximate covalent radius in layout units (scaled for rendering).
fn covalent_radius(atomic_number: u8) -> f32 {
    match atomic_number {
        1 => 0.15,  // H
        6 => 0.35,  // C
        7 => 0.30,  // N
        8 => 0.28,  // O
        16 => 0.45, // S
        _ => 0.30,  // default
    }
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

    #[test]
    fn molecular_structure_water() {
        // H₂O: 2 hydrogen + 1 oxygen
        let mol = crate::molecule::Molecule::new(&[(1, 2), (8, 1)]);
        let viz = MolecularStructure::from_molecule(&mol);
        assert_eq!(viz.atoms.len(), 3); // 2H + 1O
        assert!(!viz.bonds.is_empty());
        assert!(viz.formula.contains('H'));
        assert!(viz.formula.contains('O'));
        assert!(viz.molecular_weight > 17.0 && viz.molecular_weight < 19.0);
    }

    #[test]
    fn molecular_structure_single_atom() {
        let mol = crate::molecule::Molecule::new(&[(2, 1)]); // He
        let viz = MolecularStructure::from_molecule(&mol);
        assert_eq!(viz.atoms.len(), 1);
        assert!(viz.bonds.is_empty());
    }

    #[test]
    fn covalent_radius_known() {
        assert!((covalent_radius(1) - 0.15).abs() < 0.01); // H
        assert!((covalent_radius(6) - 0.35).abs() < 0.01); // C
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
