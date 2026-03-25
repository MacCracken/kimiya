//! Reaction network analysis — sparse stoichiometric matrices, pathway analysis.
//!
//! Uses [`hisab`] sparse matrices (CsrMatrix) and SVD for efficient
//! analysis of large reaction networks.

use crate::error::{KimiyaError, Result};

/// A reaction defined as (reactants, products) where each is a list of (species_index, coefficient).
pub type ReactionDef = (Vec<(usize, f64)>, Vec<(usize, f64)>);

/// Build a sparse stoichiometric matrix from reaction definitions.
///
/// Each reaction is defined as:
/// - `reactants`: \[(species_index, coefficient)\]
/// - `products`: \[(species_index, coefficient)\]
///
/// Returns a dense stoichiometric matrix where rows = species, columns = reactions.
/// Positive entries = products, negative = reactants.
pub fn stoichiometric_matrix(n_species: usize, reactions: &[ReactionDef]) -> Result<Vec<Vec<f64>>> {
    if n_species == 0 {
        return Err(KimiyaError::InvalidInput(
            "need at least one species".into(),
        ));
    }
    if reactions.is_empty() {
        return Err(KimiyaError::InvalidInput(
            "need at least one reaction".into(),
        ));
    }

    let n_rxns = reactions.len();
    let mut matrix = vec![vec![0.0; n_rxns]; n_species];

    for (j, (reactants, products)) in reactions.iter().enumerate() {
        for &(species, coeff) in reactants {
            if species >= n_species {
                return Err(KimiyaError::InvalidInput(format!(
                    "species index {species} out of range (n_species={n_species})"
                )));
            }
            matrix[species][j] -= coeff;
        }
        for &(species, coeff) in products {
            if species >= n_species {
                return Err(KimiyaError::InvalidInput(format!(
                    "species index {species} out of range (n_species={n_species})"
                )));
            }
            matrix[species][j] += coeff;
        }
    }

    Ok(matrix)
}

/// Build a sparse (CSR) stoichiometric matrix.
///
/// More memory-efficient for large networks with many species and few reactions per species.
pub fn sparse_stoichiometric_matrix(
    n_species: usize,
    reactions: &[ReactionDef],
) -> Result<hisab::CsrMatrix> {
    let dense = stoichiometric_matrix(n_species, reactions)?;
    Ok(hisab::CsrMatrix::from_dense(&dense))
}

/// Compute the rank of a stoichiometric matrix.
///
/// The rank tells you the number of independent reactions.
///
/// # Errors
///
/// Returns error if SVD fails.
pub fn reaction_network_rank(matrix: &[Vec<f64>]) -> Result<usize> {
    hisab::num::matrix_rank(matrix, Some(1e-10))
        .map_err(|e| KimiyaError::ComputationError(format!("SVD rank computation failed: {e}")))
}

/// Compute the null space dimension (number of conservation laws).
///
/// conservation_laws = n_species - rank(S)
pub fn conservation_law_count(n_species: usize, rank: usize) -> usize {
    n_species.saturating_sub(rank)
}

/// Rate vector: dx/dt = S × v where S is stoichiometric matrix and v is reaction rate vector.
///
/// - `stoich_matrix`: S (n_species × n_reactions)
/// - `rates`: v (n_reactions)
///
/// Returns dx/dt (n_species).
///
/// # Errors
///
/// Returns error if dimensions don't match.
pub fn rate_from_stoichiometry(stoich_matrix: &[Vec<f64>], rates: &[f64]) -> Result<Vec<f64>> {
    if stoich_matrix.is_empty() {
        return Err(KimiyaError::InvalidInput("empty matrix".into()));
    }
    let n_rxns = stoich_matrix[0].len();
    if rates.len() != n_rxns {
        return Err(KimiyaError::InvalidInput(format!(
            "rate vector length {} != n_reactions {n_rxns}",
            rates.len()
        )));
    }

    Ok(stoich_matrix
        .iter()
        .map(|row| row.iter().zip(rates.iter()).map(|(&s, &v)| s * v).sum())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stoich_matrix_basic() {
        // A → B (species 0 = A, species 1 = B)
        let m = stoichiometric_matrix(2, &[(vec![(0, 1.0)], vec![(1, 1.0)])]).unwrap();
        assert!((m[0][0] - (-1.0)).abs() < f64::EPSILON); // A consumed
        assert!((m[1][0] - 1.0).abs() < f64::EPSILON); // B produced
    }

    #[test]
    fn stoich_matrix_two_reactions() {
        // R1: A → B, R2: B → C (3 species, 2 reactions)
        let m = stoichiometric_matrix(
            3,
            &[
                (vec![(0, 1.0)], vec![(1, 1.0)]),
                (vec![(1, 1.0)], vec![(2, 1.0)]),
            ],
        )
        .unwrap();
        assert_eq!(m.len(), 3); // 3 species
        assert_eq!(m[0].len(), 2); // 2 reactions
        // A: [-1, 0], B: [1, -1], C: [0, 1]
        assert!((m[0][0] - (-1.0)).abs() < f64::EPSILON);
        assert!((m[1][0] - 1.0).abs() < f64::EPSILON);
        assert!((m[1][1] - (-1.0)).abs() < f64::EPSILON);
        assert!((m[2][1] - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn stoich_matrix_invalid_species() {
        assert!(stoichiometric_matrix(2, &[(vec![(5, 1.0)], vec![(1, 1.0)])]).is_err());
    }

    #[test]
    fn sparse_matrix_basic() {
        let csr = sparse_stoichiometric_matrix(2, &[(vec![(0, 1.0)], vec![(1, 1.0)])]).unwrap();
        assert_eq!(csr.nrows, 2);
    }

    #[test]
    fn network_rank_independent_reactions() {
        // A → B, B → C: rank should be 2 (both independent)
        let m = stoichiometric_matrix(
            3,
            &[
                (vec![(0, 1.0)], vec![(1, 1.0)]),
                (vec![(1, 1.0)], vec![(2, 1.0)]),
            ],
        )
        .unwrap();
        let r = reaction_network_rank(&m).unwrap();
        assert_eq!(r, 2);
    }

    #[test]
    fn conservation_laws_abc() {
        // A → B → C: 3 species, rank 2 → 1 conservation law (total mass)
        let count = conservation_law_count(3, 2);
        assert_eq!(count, 1);
    }

    #[test]
    fn rate_from_stoich_basic() {
        // A → B with rate v=0.5
        let m = stoichiometric_matrix(2, &[(vec![(0, 1.0)], vec![(1, 1.0)])]).unwrap();
        let dx = rate_from_stoichiometry(&m, &[0.5]).unwrap();
        assert!((dx[0] - (-0.5)).abs() < f64::EPSILON); // A consumed
        assert!((dx[1] - 0.5).abs() < f64::EPSILON); // B produced
    }

    #[test]
    fn rate_dimension_mismatch() {
        let m = stoichiometric_matrix(2, &[(vec![(0, 1.0)], vec![(1, 1.0)])]).unwrap();
        assert!(rate_from_stoichiometry(&m, &[0.5, 0.3]).is_err());
    }
}
