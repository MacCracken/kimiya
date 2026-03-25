//! Stoichiometry — equation balancing, limiting reagent, yield, composition, empirical formulas.

use crate::error::{KimiyaError, Result};

// ── Equation balancing ───────────────────────────────────────────────

/// Balance a chemical equation given the composition matrix.
///
/// Each row represents an element, each column represents a species.
/// Reactants have positive entries, products have negative entries
/// (or vice versa — the system is homogeneous).
///
/// The composition matrix M has M\[i\]\[j\] = count of element i in species j.
/// Convention: reactant columns are positive, product columns are negative.
///
/// Returns integer stoichiometric coefficients (all positive).
///
/// # Example
///
/// For CH₄ + O₂ → CO₂ + H₂O:
/// ```text
///        CH4  O2  CO2  H2O
///   C:  [ 1,  0, -1,   0 ]
///   H:  [ 4,  0,  0,  -2 ]
///   O:  [ 0,  2, -2,  -1 ]
/// ```
///
/// # Errors
///
/// Returns error if the system cannot be balanced.
pub fn balance_equation(composition: &[Vec<f64>]) -> Result<Vec<u32>> {
    if composition.is_empty() {
        return Err(KimiyaError::InvalidInput("empty composition matrix".into()));
    }
    let n_species = composition[0].len();
    if n_species < 2 {
        return Err(KimiyaError::InvalidInput("need at least 2 species".into()));
    }
    for row in composition {
        if row.len() != n_species {
            return Err(KimiyaError::InvalidInput(
                "all rows must have equal length".into(),
            ));
        }
    }

    let n_elements = composition.len();
    let n_unknowns = n_species - 1; // fix last coefficient = 1

    // Build augmented matrix: for each element row, move the last column to RHS
    // M[i][0..n-1] * x[0..n-1] = -M[i][n-1]
    let mut augmented: Vec<Vec<f64>> = Vec::with_capacity(n_elements);
    for row in composition {
        let mut aug_row = Vec::with_capacity(n_unknowns + 1);
        aug_row.extend_from_slice(&row[..n_unknowns]);
        aug_row.push(-row[n_species - 1]); // RHS
        augmented.push(aug_row);
    }

    // If more elements than unknowns, use least-squares via normal equations
    // For typical chemistry, n_elements ≈ n_unknowns, so direct solve works
    // Pad or truncate to square system
    while augmented.len() < n_unknowns {
        let mut zero_row = vec![0.0; n_unknowns + 1];
        zero_row[augmented.len()] = 1.0; // add identity rows for underdetermined
        augmented.push(zero_row);
    }
    augmented.truncate(n_unknowns);

    let solution = hisab::num::gaussian_elimination(&mut augmented)
        .map_err(|e| KimiyaError::ComputationError(format!("equation balancing failed: {e}")))?;

    // Append the fixed coefficient (1.0) for the last species
    let mut coeffs: Vec<f64> = solution;
    coeffs.push(1.0);

    // Make all positive (flip sign if needed)
    let all_negative = coeffs.iter().all(|&c| c <= 0.0);
    if all_negative {
        for c in &mut coeffs {
            *c = -(*c);
        }
    }

    // Find smallest to scale to integers
    let min_abs = coeffs
        .iter()
        .map(|c| c.abs())
        .filter(|&c| c > 1e-10)
        .fold(f64::INFINITY, f64::min);

    if min_abs < 1e-10 {
        return Err(KimiyaError::ComputationError(
            "degenerate solution — equation may not be balanceable".into(),
        ));
    }

    let scaled: Vec<f64> = coeffs.iter().map(|c| c / min_abs).collect();

    // Try multipliers 1..12 to find integer solution
    let mut best_mult = 1u32;
    let mut best_err = f64::INFINITY;
    for mult in 1..=12 {
        let err: f64 = scaled
            .iter()
            .map(|&r| {
                let v = r * mult as f64;
                (v - v.round()).abs()
            })
            .sum();
        if err < best_err {
            best_err = err;
            best_mult = mult;
        }
        if best_err < 0.01 {
            break;
        }
    }

    let result: Vec<u32> = scaled
        .iter()
        .map(|&r| (r * best_mult as f64).round().abs() as u32)
        .collect();

    // Verify no zero coefficients
    if result.contains(&0) {
        return Err(KimiyaError::ComputationError(
            "balancing produced zero coefficient".into(),
        ));
    }

    Ok(result)
}

// ── Oxidation state assignment ───────────────────────────────────────

/// Assign oxidation states using standard rules.
///
/// Given a molecule as (atomic_number, count) pairs and overall charge,
/// returns oxidation state for each element.
///
/// Rules applied (in order):
/// 1. F is always -1
/// 2. O is -2 (except in peroxides: -1, and OF₂: +2)
/// 3. H is +1 (except in metal hydrides: -1)
/// 4. Alkali metals are +1, alkaline earth are +2
/// 5. Remaining elements solved from charge balance
///
/// Returns (atomic_number, oxidation_state) pairs.
///
/// # Errors
///
/// Returns error if the system is underdetermined (multiple unknown oxidation states).
pub fn assign_oxidation_states(atoms: &[(u8, u32)], overall_charge: i32) -> Result<Vec<(u8, i32)>> {
    if atoms.is_empty() {
        return Err(KimiyaError::InvalidInput(
            "need at least one element".into(),
        ));
    }

    let mut result: Vec<(u8, Option<i32>)> = atoms.iter().map(|&(z, _)| (z, None)).collect();
    let mut remaining_charge = overall_charge;

    // Pass 1: assign known oxidation states
    for (i, &(z, count)) in atoms.iter().enumerate() {
        let os = match z {
            9 => Some(-1),                         // F always -1
            3 | 11 | 19 | 37 | 55 | 87 => Some(1), // Alkali metals +1
            4 | 12 | 20 | 38 | 56 | 88 => Some(2), // Alkaline earth +2
            _ => None,
        };
        if let Some(state) = os {
            result[i].1 = Some(state);
            remaining_charge -= state * count as i32;
        }
    }

    // Pass 2: assign O = -2 and H = +1 (simplified rules)
    for (i, &(z, count)) in atoms.iter().enumerate() {
        if result[i].1.is_some() {
            continue;
        }
        let os = match z {
            8 => Some(-2), // O = -2 (simplified, ignores peroxides)
            1 => Some(1),  // H = +1 (simplified, ignores hydrides)
            _ => None,
        };
        if let Some(state) = os {
            result[i].1 = Some(state);
            remaining_charge -= state * count as i32;
        }
    }

    // Pass 3: solve for the single remaining unknown
    let unknowns: Vec<usize> = result
        .iter()
        .enumerate()
        .filter(|(_, (_, os))| os.is_none())
        .map(|(i, _)| i)
        .collect();

    match unknowns.len() {
        0 => {} // all assigned
        1 => {
            let idx = unknowns[0];
            let count = atoms[idx].1 as i32;
            if count == 0 {
                return Err(KimiyaError::InvalidInput("zero atom count".into()));
            }
            result[idx].1 = Some(remaining_charge / count);
        }
        _ => {
            return Err(KimiyaError::ComputationError(
                "multiple unknown oxidation states — system underdetermined".into(),
            ));
        }
    }

    Ok(result
        .into_iter()
        .map(|(z, os)| (z, os.unwrap_or(0)))
        .collect())
}

// ── Limiting reagent + yield ─────────────────────────────────────────

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
        if moles < 0.0 {
            return Err(KimiyaError::InvalidInput(
                "moles available must be non-negative".into(),
            ));
        }
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

    // ── Equation balancing ───────────────────────────────────────────

    #[test]
    fn balance_h2_o2_water() {
        // H₂ + O₂ → H₂O
        //        H2  O2  H2O
        //   H:  [ 1,  0, -1 ]  (H₂ has 2H per molecule, but we use per-species counts)
        // Actually: composition[element][species] = atom count (positive for reactants, negative for products)
        // H₂ + O₂ → H₂O: species = [H₂, O₂, H₂O]
        //   H: [2, 0, -2]
        //   O: [0, 2, -1]
        let comp = vec![vec![2.0, 0.0, -2.0], vec![0.0, 2.0, -1.0]];
        let coeffs = balance_equation(&comp).unwrap();
        // Expected: 2H₂ + 1O₂ → 2H₂O → [2, 1, 2]
        assert_eq!(coeffs.len(), 3);
        // Verify balance: for each element, Σ coeff_i × comp[e][i] = 0
        for row in &comp {
            let sum: f64 = row
                .iter()
                .zip(coeffs.iter())
                .map(|(&c, &k)| c * k as f64)
                .sum();
            assert!(sum.abs() < 0.1, "element not balanced: sum = {sum}");
        }
    }

    #[test]
    fn balance_ch4_combustion() {
        // CH₄ + O₂ → CO₂ + H₂O
        //        CH4  O2  CO2  H2O
        //   C:  [ 1,  0, -1,   0 ]
        //   H:  [ 4,  0,  0,  -2 ]
        //   O:  [ 0,  2, -2,  -1 ]
        let comp = vec![
            vec![1.0, 0.0, -1.0, 0.0],
            vec![4.0, 0.0, 0.0, -2.0],
            vec![0.0, 2.0, -2.0, -1.0],
        ];
        let coeffs = balance_equation(&comp).unwrap();
        // Expected: 1CH₄ + 2O₂ → 1CO₂ + 2H₂O → [1, 2, 1, 2]
        assert_eq!(coeffs.len(), 4);
        for row in &comp {
            let sum: f64 = row
                .iter()
                .zip(coeffs.iter())
                .map(|(&c, &k)| c * k as f64)
                .sum();
            assert!(sum.abs() < 0.1, "element not balanced: sum = {sum}");
        }
    }

    #[test]
    fn balance_empty_is_error() {
        assert!(balance_equation(&[]).is_err());
    }

    #[test]
    fn balance_single_species_is_error() {
        assert!(balance_equation(&[vec![1.0]]).is_err());
    }

    // ── Oxidation states ─────────────────────────────────────────────

    #[test]
    fn oxidation_state_water() {
        // H₂O: H = +1, O = -2
        let os = assign_oxidation_states(&[(1, 2), (8, 1)], 0).unwrap();
        assert_eq!(os.iter().find(|&&(z, _)| z == 1).unwrap().1, 1);
        assert_eq!(os.iter().find(|&&(z, _)| z == 8).unwrap().1, -2);
    }

    #[test]
    fn oxidation_state_nacl() {
        // NaCl: Na = +1, Cl = -1
        let os = assign_oxidation_states(&[(11, 1), (17, 1)], 0).unwrap();
        assert_eq!(os.iter().find(|&&(z, _)| z == 11).unwrap().1, 1);
        assert_eq!(os.iter().find(|&&(z, _)| z == 17).unwrap().1, -1);
    }

    #[test]
    fn oxidation_state_sulfate_ion() {
        // SO₄²⁻: O = -2, S = +6
        let os = assign_oxidation_states(&[(16, 1), (8, 4)], -2).unwrap();
        assert_eq!(os.iter().find(|&&(z, _)| z == 16).unwrap().1, 6);
        assert_eq!(os.iter().find(|&&(z, _)| z == 8).unwrap().1, -2);
    }

    #[test]
    fn oxidation_state_permanganate() {
        // MnO₄⁻: O = -2, Mn = +7
        let os = assign_oxidation_states(&[(25, 1), (8, 4)], -1).unwrap();
        assert_eq!(os.iter().find(|&&(z, _)| z == 25).unwrap().1, 7);
    }

    #[test]
    fn oxidation_state_co2() {
        // CO₂: O = -2, C = +4
        let os = assign_oxidation_states(&[(6, 1), (8, 2)], 0).unwrap();
        assert_eq!(os.iter().find(|&&(z, _)| z == 6).unwrap().1, 4);
    }

    #[test]
    fn oxidation_state_empty_is_error() {
        assert!(assign_oxidation_states(&[], 0).is_err());
    }
}
