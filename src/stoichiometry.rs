//! Stoichiometry — limiting reagent, yield, composition, empirical formulas.

use crate::error::{KimiyaError, Result};

/// Identify the limiting reagent and return its index.
///
/// - `moles_available`: moles of each reactant available
/// - `stoich_coefficients`: stoichiometric coefficients for each reactant
///
/// The limiting reagent is the one with the smallest ratio moles/coefficient.
///
/// # Errors
///
/// Returns error if arrays differ in length, are empty, or contain non-positive values.
pub fn limiting_reagent(moles_available: &[f64], stoich_coefficients: &[f64]) -> Result<usize> {
    if moles_available.len() != stoich_coefficients.len() {
        return Err(KimiyaError::InvalidInput(
            "moles and coefficients must have equal length".into(),
        ));
    }
    if moles_available.is_empty() {
        return Err(KimiyaError::InvalidInput(
            "need at least one reactant".into(),
        ));
    }
    let mut min_ratio = f64::INFINITY;
    let mut min_idx = 0;
    for (i, (&moles, &coeff)) in moles_available
        .iter()
        .zip(stoich_coefficients.iter())
        .enumerate()
    {
        if coeff <= 0.0 {
            return Err(KimiyaError::InvalidInput(
                "stoichiometric coefficients must be positive".into(),
            ));
        }
        let ratio = moles / coeff;
        if ratio < min_ratio {
            min_ratio = ratio;
            min_idx = i;
        }
    }
    Ok(min_idx)
}

/// Theoretical yield in moles of a product.
///
/// - `moles_available`: moles of each reactant
/// - `stoich_reactants`: stoichiometric coefficients of reactants
/// - `stoich_product`: stoichiometric coefficient of the desired product
///
/// # Errors
///
/// Returns error if inputs are invalid.
pub fn theoretical_yield_moles(
    moles_available: &[f64],
    stoich_reactants: &[f64],
    stoich_product: f64,
) -> Result<f64> {
    if stoich_product <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "product coefficient must be positive".into(),
        ));
    }
    let idx = limiting_reagent(moles_available, stoich_reactants)?;
    let limiting_ratio = moles_available[idx] / stoich_reactants[idx];
    Ok(limiting_ratio * stoich_product)
}

/// Percent yield: (actual / theoretical) × 100.
///
/// # Errors
///
/// Returns error if theoretical yield is not positive.
#[inline]
pub fn percent_yield(actual: f64, theoretical: f64) -> Result<f64> {
    if theoretical <= 0.0 {
        return Err(KimiyaError::InvalidInput(
            "theoretical yield must be positive".into(),
        ));
    }
    Ok(actual / theoretical * 100.0)
}

/// Mass percent composition of each element in a molecule.
///
/// Returns a vector of (atomic_number, mass_percent) pairs.
///
/// # Errors
///
/// Returns error if molecular weight calculation fails.
pub fn percent_composition(molecule: &crate::molecule::Molecule) -> Result<Vec<(u8, f64)>> {
    let mw = molecule.molecular_weight()?;
    if mw <= 0.0 {
        return Err(KimiyaError::ComputationError(
            "molecular weight must be positive".into(),
        ));
    }
    let mut result = Vec::with_capacity(molecule.atoms.len());
    for atom in &molecule.atoms {
        let element = crate::element::lookup_by_number(atom.element_number).ok_or_else(|| {
            KimiyaError::InvalidElement(format!("unknown atomic number {}", atom.element_number))
        })?;
        let mass_fraction = element.atomic_mass * atom.count as f64 / mw * 100.0;
        result.push((atom.element_number, mass_fraction));
    }
    Ok(result)
}

/// Determine empirical formula from mass percentages.
///
/// - `elements`: array of (atomic_number, mass_percent) pairs
///
/// Returns array of (atomic_number, count) pairs with integer ratios.
///
/// # Errors
///
/// Returns error if any element is unknown or mass percentages are invalid.
pub fn empirical_formula(elements: &[(u8, f64)]) -> Result<Vec<(u8, u32)>> {
    if elements.is_empty() {
        return Err(KimiyaError::InvalidInput(
            "need at least one element".into(),
        ));
    }

    // Convert mass% to moles (assume 100g sample)
    let mut moles = Vec::with_capacity(elements.len());
    for &(z, mass_pct) in elements {
        let element = crate::element::lookup_by_number(z)
            .ok_or_else(|| KimiyaError::InvalidElement(format!("unknown atomic number {z}")))?;
        if mass_pct <= 0.0 {
            return Err(KimiyaError::InvalidInput(
                "mass percentages must be positive".into(),
            ));
        }
        moles.push((z, mass_pct / element.atomic_mass));
    }

    // Divide by smallest mole value
    let min_moles = moles.iter().map(|(_, m)| *m).fold(f64::INFINITY, f64::min);
    if min_moles <= 0.0 {
        return Err(KimiyaError::ComputationError(
            "mole calculation produced non-positive value".into(),
        ));
    }

    let ratios: Vec<(u8, f64)> = moles.iter().map(|&(z, m)| (z, m / min_moles)).collect();

    // Find multiplier to make all ratios close to integers (try 1..6)
    let mut best_mult = 1u32;
    let mut best_err = f64::INFINITY;
    for mult in 1..=6 {
        let err: f64 = ratios
            .iter()
            .map(|&(_, r)| {
                let scaled = r * mult as f64;
                (scaled - scaled.round()).abs()
            })
            .sum();
        if err < best_err {
            best_err = err;
            best_mult = mult;
        }
    }

    Ok(ratios
        .iter()
        .map(|&(z, r)| (z, (r * best_mult as f64).round() as u32))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;

    #[test]
    fn limiting_reagent_basic() {
        // 2H₂ + O₂ → 2H₂O, have 3 mol H₂ and 2 mol O₂
        // H₂: 3/2 = 1.5, O₂: 2/1 = 2.0 → H₂ is limiting
        let idx = limiting_reagent(&[3.0, 2.0], &[2.0, 1.0]).unwrap();
        assert_eq!(idx, 0);
    }

    #[test]
    fn limiting_reagent_second() {
        // H₂: 10/2 = 5.0, O₂: 1/1 = 1.0 → O₂ is limiting
        let idx = limiting_reagent(&[10.0, 1.0], &[2.0, 1.0]).unwrap();
        assert_eq!(idx, 1);
    }

    #[test]
    fn limiting_reagent_empty_is_error() {
        assert!(limiting_reagent(&[], &[]).is_err());
    }

    #[test]
    fn limiting_reagent_mismatched_is_error() {
        assert!(limiting_reagent(&[1.0], &[1.0, 2.0]).is_err());
    }

    #[test]
    fn theoretical_yield_basic() {
        // 2H₂ + O₂ → 2H₂O, have 3 mol H₂ and 2 mol O₂
        // Limiting: H₂ (ratio 1.5), yield = 1.5 × 2 = 3.0 mol H₂O
        let y = theoretical_yield_moles(&[3.0, 2.0], &[2.0, 1.0], 2.0).unwrap();
        assert!((y - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn percent_yield_basic() {
        let pct = percent_yield(8.0, 10.0).unwrap();
        assert!((pct - 80.0).abs() < f64::EPSILON);
    }

    #[test]
    fn percent_yield_over_100() {
        // Possible with impurities
        let pct = percent_yield(11.0, 10.0).unwrap();
        assert!((pct - 110.0).abs() < 1e-10);
    }

    #[test]
    fn percent_yield_zero_theoretical_is_error() {
        assert!(percent_yield(5.0, 0.0).is_err());
    }

    #[test]
    fn percent_composition_water() {
        let water = Molecule::water();
        let comp = percent_composition(&water).unwrap();
        // H: 2×1.008/18.015 ≈ 11.19%, O: 15.999/18.015 ≈ 88.81%
        let h_pct = comp.iter().find(|&&(z, _)| z == 1).unwrap().1;
        let o_pct = comp.iter().find(|&&(z, _)| z == 8).unwrap().1;
        assert!((h_pct - 11.19).abs() < 0.1);
        assert!((o_pct - 88.81).abs() < 0.1);
        assert!((h_pct + o_pct - 100.0).abs() < 0.01);
    }

    #[test]
    fn empirical_formula_water() {
        // H: 11.19%, O: 88.81%
        let formula = empirical_formula(&[(1, 11.19), (8, 88.81)]).unwrap();
        // Should give H₂O → [(1, 2), (8, 1)]
        assert_eq!(formula.iter().find(|&&(z, _)| z == 1).unwrap().1, 2);
        assert_eq!(formula.iter().find(|&&(z, _)| z == 8).unwrap().1, 1);
    }

    #[test]
    fn empirical_formula_glucose() {
        // CH₂O empirical (glucose = C₆H₁₂O₆ → CH₂O)
        // C: 40.0%, H: 6.7%, O: 53.3%
        let formula = empirical_formula(&[(6, 40.0), (1, 6.7), (8, 53.3)]).unwrap();
        let c = formula.iter().find(|&&(z, _)| z == 6).unwrap().1;
        let h = formula.iter().find(|&&(z, _)| z == 1).unwrap().1;
        let o = formula.iter().find(|&&(z, _)| z == 8).unwrap().1;
        assert_eq!(c, 1);
        assert_eq!(h, 2);
        assert_eq!(o, 1);
    }

    #[test]
    fn empirical_formula_empty_is_error() {
        assert!(empirical_formula(&[]).is_err());
    }
}
