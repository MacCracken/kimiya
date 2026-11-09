//! Stoichiometry — equation balancing, limiting reagent, yield, composition, empirical formulas.

use crate::element::ElementCategory;
use crate::error::{KimiyaError, Result};
use crate::validate::{is_finite, is_non_negative, is_positive};

// ── Equation balancing ───────────────────────────────────────────────

/// Largest magnitude accepted for a composition entry (2⁵³: every integer up
/// to this is exactly representable as `f64`).
const MAX_COMPOSITION_ENTRY: f64 = 9_007_199_254_740_992.0;

/// Tolerance for treating a composition entry as an integer.
const INTEGER_TOLERANCE: f64 = 1e-9;

/// Overflow error for the exact integer arithmetic in [`balance_equation`].
fn balance_overflow() -> KimiyaError {
    KimiyaError::ComputationError("integer overflow while balancing equation".into())
}

/// Greatest common divisor of two unsigned integers (Euclid); `gcd(0, 0) = 0`.
#[inline]
fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Divide `row` by the gcd of its entries (no-op for an all-zero row).
fn reduce_row(row: &mut [i64]) -> Result<()> {
    let g = row.iter().fold(0u64, |g, &v| gcd_u64(g, v.unsigned_abs()));
    if g > 1 {
        let g = i64::try_from(g).map_err(|_| balance_overflow())?;
        for v in row.iter_mut() {
            *v /= g;
        }
    }
    Ok(())
}

/// Balance a chemical equation given the composition matrix.
///
/// The composition matrix M has M\[i\]\[j\] = count of element i in species j:
/// each row is an element, each column a species. Convention: reactant
/// columns are positive, product columns are negative. Every entry must be
/// an integer (within 1e-9) of magnitude ≤ 2⁵³.
///
/// Returns the smallest positive integer stoichiometric coefficients x
/// (one per species) with M·x = 0.
///
/// # Algorithm
///
/// Exact fraction-free Gauss–Jordan elimination over the integers (`i64`,
/// checked arithmetic, every row divided by the gcd of its entries after
/// each step) on the full m×n matrix — no rows are dropped, so the result
/// does not depend on row order. The null space of M must be exactly one
/// dimensional; its basis vector is built from the reduced row echelon form,
/// scaled to the smallest integers, sign-normalised, converted to `u32`,
/// and verified against every original row (M·x = 0 exactly).
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
/// gives `[1, 2, 1, 2]`.
///
/// # Errors
///
/// - [`KimiyaError::InvalidInput`] if the matrix is empty, has fewer than
///   two species, has rows of unequal length, or contains an entry that is
///   not a finite integer of magnitude ≤ 2⁵³.
/// - [`KimiyaError::ComputationError`] if the equation cannot be balanced
///   (only the trivial solution), is ambiguous (two or more independent
///   reactions, null space dimension > 1), needs coefficients of mixed sign
///   (a species is on the wrong side), needs a zero coefficient, or if the
///   integer arithmetic overflows or a coefficient exceeds `u32::MAX`.
pub fn balance_equation(composition: &[Vec<f64>]) -> Result<Vec<u32>> {
    let Some(first) = composition.first() else {
        return Err(KimiyaError::InvalidInput("empty composition matrix".into()));
    };
    let n = first.len();
    if n < 2 {
        return Err(KimiyaError::InvalidInput("need at least 2 species".into()));
    }
    if composition.iter().any(|row| row.len() != n) {
        return Err(KimiyaError::InvalidInput(
            "all rows must have equal length".into(),
        ));
    }
    let m = composition.len();

    // Exact integer copy of the matrix, row-major.
    let mut a: Vec<i64> = Vec::with_capacity(m * n);
    for &v in composition.iter().flatten() {
        let r = v.round();
        if !is_finite(v) || (v - r).abs() > INTEGER_TOLERANCE || r.abs() > MAX_COMPOSITION_ENTRY {
            return Err(KimiyaError::InvalidInput(
                "composition entries must be integers".into(),
            ));
        }
        // Exact: |r| ≤ 2⁵³ is an integer-valued f64.
        a.push(r as i64);
    }

    // Fraction-free Gauss–Jordan elimination to reduced row echelon form.
    let mut rank = 0usize;
    let mut free_col: Option<usize> = None;
    for col in 0..n {
        let Some(p) = (rank..m).find(|&r| a[r * n + col] != 0) else {
            if free_col.is_some() {
                return Err(KimiyaError::ComputationError(
                    "equation balancing is ambiguous: more than one independent reaction \
                     (null space dimension > 1)"
                        .into(),
                ));
            }
            free_col = Some(col);
            continue;
        };
        if p != rank {
            for j in 0..n {
                a.swap(p * n + j, rank * n + j);
            }
        }
        let pr = rank * n;
        reduce_row(&mut a[pr..pr + n])?;
        if a[pr + col] < 0 {
            for v in &mut a[pr..pr + n] {
                *v = v.checked_neg().ok_or_else(balance_overflow)?;
            }
        }
        let pivot = a[pr + col];
        for r in 0..m {
            let rr = r * n;
            let factor = a[rr + col];
            if r == rank || factor == 0 {
                continue;
            }
            // row_r ← pivot·row_r − factor·row_pivot (zeroes column `col`)
            for j in 0..n {
                let lhs = pivot.checked_mul(a[rr + j]).ok_or_else(balance_overflow)?;
                let rhs = factor.checked_mul(a[pr + j]).ok_or_else(balance_overflow)?;
                a[rr + j] = lhs.checked_sub(rhs).ok_or_else(balance_overflow)?;
            }
            reduce_row(&mut a[rr..rr + n])?;
        }
        rank += 1;
    }
    let Some(f) = free_col else {
        return Err(KimiyaError::ComputationError(
            "equation cannot be balanced: only the trivial solution exists".into(),
        ));
    };

    // Null vector from the RREF: each pivot row k reads
    //   a[k][c_k]·x[c_k] + a[k][f]·x[f] = 0.
    // Take x[f] = lcm of the (positive) pivots so every x[c_k] is an integer.
    let mut lcm: i64 = 1;
    for k in 0..rank {
        let row = &a[k * n..k * n + n];
        if let Some(&piv) = row.iter().find(|&&v| v != 0) {
            let g = i64::try_from(gcd_u64(lcm.unsigned_abs(), piv.unsigned_abs()))
                .map_err(|_| balance_overflow())?;
            lcm = (lcm / g).checked_mul(piv).ok_or_else(balance_overflow)?;
        }
    }
    let mut x = vec![0i64; n];
    x[f] = lcm;
    for k in 0..rank {
        let row = &a[k * n..k * n + n];
        if let Some(c) = row.iter().position(|&v| v != 0) {
            let piv = row[c];
            x[c] = row[f]
                .checked_mul(lcm / piv)
                .and_then(i64::checked_neg)
                .ok_or_else(balance_overflow)?;
        }
    }
    reduce_row(&mut x)?;

    // Sign convention: all coefficients must share one sign.
    if x.contains(&0) {
        return Err(KimiyaError::ComputationError(
            "balancing produced zero coefficient".into(),
        ));
    }
    let positive = x.iter().filter(|&&v| v > 0).count();
    if positive == 0 {
        for v in &mut x {
            *v = v.checked_neg().ok_or_else(balance_overflow)?;
        }
    } else if positive != n {
        return Err(KimiyaError::ComputationError(
            "balancing requires coefficients of mixed sign: reactant columns must be \
             positive and product columns negative"
                .into(),
        ));
    }

    let coeffs = x
        .into_iter()
        .map(|v| {
            u32::try_from(v).map_err(|_| {
                KimiyaError::ComputationError("stoichiometric coefficient exceeds u32".into())
            })
        })
        .collect::<Result<Vec<u32>>>()?;

    // Verify M·x = 0 exactly on every original row (|entry| ≤ 2⁵³ and
    // coefficient ≤ 2³², so i128 cannot overflow in practice).
    for row in composition {
        let mut sum: i128 = 0;
        for (&v, &k) in row.iter().zip(&coeffs) {
            let term = (v.round() as i128)
                .checked_mul(i128::from(k))
                .ok_or_else(balance_overflow)?;
            sum = sum.checked_add(term).ok_or_else(balance_overflow)?;
        }
        if sum != 0 {
            return Err(KimiyaError::ComputationError(
                "balancing verification failed: element not conserved".into(),
            ));
        }
    }
    Ok(coeffs)
}

// ── Oxidation state assignment ───────────────────────────────────────

/// Whether an element category is a metal (for the metal-hydride rule).
#[inline]
const fn is_metal(category: ElementCategory) -> bool {
    matches!(
        category,
        ElementCategory::AlkaliMetal
            | ElementCategory::AlkalineEarth
            | ElementCategory::TransitionMetal
            | ElementCategory::PostTransitionMetal
            | ElementCategory::Lanthanide
            | ElementCategory::Actinide
    )
}

/// Overflow error for [`assign_oxidation_states`].
fn oxidation_overflow() -> KimiyaError {
    KimiyaError::ComputationError("oxidation state arithmetic overflow".into())
}

/// Error when rule-assigned oxidation states cannot match the overall charge.
fn unbalanced_oxidation_states() -> KimiyaError {
    KimiyaError::ComputationError(
        "rule-assigned oxidation states do not balance the overall charge".into(),
    )
}

/// `num / den` if it divides exactly, else a "non-integer average" error.
#[inline]
fn exact_oxidation_state(num: i64, den: i64) -> Result<i32> {
    if den == 0 {
        return Err(KimiyaError::InvalidInput("zero atom count".into()));
    }
    if num.checked_rem(den).ok_or_else(oxidation_overflow)? != 0 {
        return Err(KimiyaError::ComputationError(
            "non-integer average oxidation state".into(),
        ));
    }
    let q = num.checked_div(den).ok_or_else(oxidation_overflow)?;
    i32::try_from(q)
        .map_err(|_| KimiyaError::ComputationError("oxidation state out of range".into()))
}

/// Rule-assigned oxidation state, or `None` if it must be solved for.
#[inline]
const fn rule_oxidation_state(
    z: u8,
    category: ElementCategory,
    has_fluorine: bool,
    metal_hydride: bool,
) -> Option<i32> {
    match z {
        9 => Some(-1),
        1 if metal_hydride => Some(-1),
        1 => Some(1),
        8 if has_fluorine => None,
        8 => Some(-2),
        _ => match category {
            ElementCategory::AlkaliMetal => Some(1),
            ElementCategory::AlkalineEarth => Some(2),
            _ => None,
        },
    }
}

/// Assign oxidation states using standard rules.
///
/// Given a species as (atomic_number, count) pairs and its overall charge,
/// returns one (atomic_number, oxidation_state) pair per input entry, in
/// input order. Repeated entries for the same element get the same state.
///
/// Rules:
/// 1. A single-element species (O₂, S₈, Cl⁻, O₂²⁻) gets
///    overall_charge / atom count.
/// 2. F is always −1.
/// 3. Alkali metals are +1, alkaline-earth metals +2.
/// 4. H is +1, except −1 when every other element present is a metal
///    (metal hydrides: NaH, CaH₂, LiAlH₄). Non-metal hydrides such as
///    NaBH₄ or SiH₄ are not recognised (H is taken as +1).
/// 5. O is −2, except: with F present O is solved from the charge balance
///    (OF₂ → +2); and when every atom is assigned by these rules but the
///    charges do not balance, O is re-solved from the balance (peroxides
///    H₂O₂, Na₂O₂ → −1), accepted only if it lies strictly between −2 and 0.
/// 6. The remaining element (at most one distinct unknown) is solved from
///    the charge balance.
///
/// # Errors
///
/// - [`KimiyaError::InvalidInput`] if `atoms` is empty or a count is zero.
/// - [`KimiyaError::InvalidElement`] for an unknown atomic number.
/// - [`KimiyaError::ComputationError`] if more than one distinct element is
///   unknown (underdetermined); if the solved state is not an integer (an
///   average over inequivalent atoms, e.g. C₃H₈ −8/3, Fe₃O₄ +8/3, KO₂
///   superoxide −1/2); if the rule-assigned states cannot balance the
///   charge; or on arithmetic overflow / a result outside `i32`.
pub fn assign_oxidation_states(atoms: &[(u8, u32)], overall_charge: i32) -> Result<Vec<(u8, i32)>> {
    let Some(&(first_z, _)) = atoms.first() else {
        return Err(KimiyaError::InvalidInput(
            "need at least one element".into(),
        ));
    };
    let charge = i64::from(overall_charge);

    // Pass 0: validate and collect the facts the rules depend on.
    let mut single_element = true;
    let mut has_fluorine = false;
    let mut others_all_metal = true;
    let mut total_atoms: i64 = 0;
    for &(z, count) in atoms {
        if count == 0 {
            return Err(KimiyaError::InvalidInput("zero atom count".into()));
        }
        let element = crate::element::lookup_by_number(z)
            .ok_or_else(|| KimiyaError::InvalidElement(format!("unknown atomic number {z}")))?;
        single_element &= z == first_z;
        has_fluorine |= z == 9;
        if z != 1 && !is_metal(element.category) {
            others_all_metal = false;
        }
        total_atoms = total_atoms
            .checked_add(i64::from(count))
            .ok_or_else(oxidation_overflow)?;
    }

    if single_element {
        let os = exact_oxidation_state(charge, total_atoms)?;
        return Ok(atoms.iter().map(|&(z, _)| (z, os)).collect());
    }
    // ≥ 2 distinct elements, so at least one is not H.
    let metal_hydride = others_all_metal;

    // Pass 1: rule-based assignment; track the single unknown element.
    let mut result: Vec<(u8, i32)> = Vec::with_capacity(atoms.len());
    let mut assigned_sum: i64 = 0;
    let mut oxygen_count: i64 = 0;
    let mut unknown: Option<(u8, i64)> = None;
    for &(z, count) in atoms {
        let category = crate::element::lookup_by_number(z)
            .map(|e| e.category)
            .ok_or_else(|| KimiyaError::InvalidElement(format!("unknown atomic number {z}")))?;
        let count = i64::from(count);
        match rule_oxidation_state(z, category, has_fluorine, metal_hydride) {
            Some(os) => {
                let term = i64::from(os)
                    .checked_mul(count)
                    .ok_or_else(oxidation_overflow)?;
                assigned_sum = assigned_sum
                    .checked_add(term)
                    .ok_or_else(oxidation_overflow)?;
                if z == 8 {
                    oxygen_count = oxygen_count
                        .checked_add(count)
                        .ok_or_else(oxidation_overflow)?;
                }
                result.push((z, os));
            }
            None => {
                unknown = match unknown {
                    None => Some((z, count)),
                    Some((uz, uc)) if uz == z => {
                        Some((uz, uc.checked_add(count).ok_or_else(oxidation_overflow)?))
                    }
                    Some(_) => {
                        return Err(KimiyaError::ComputationError(
                            "multiple unknown oxidation states — system underdetermined".into(),
                        ));
                    }
                };
                result.push((z, 0));
            }
        }
    }

    // Pass 2: solve the unknown, or re-solve O for peroxides.
    let remaining = charge
        .checked_sub(assigned_sum)
        .ok_or_else(oxidation_overflow)?;
    let solved = match unknown {
        Some((uz, ucount)) => Some((uz, exact_oxidation_state(remaining, ucount)?)),
        None if remaining == 0 => None,
        None if oxygen_count > 0 => {
            // O was given the default −2; re-solve it from the balance
            // (peroxide −1; superoxide −1/2 is rejected as non-integer).
            let o_total = remaining
                .checked_add(
                    oxygen_count
                        .checked_mul(-2)
                        .ok_or_else(oxidation_overflow)?,
                )
                .ok_or_else(oxidation_overflow)?;
            let os = exact_oxidation_state(o_total, oxygen_count)?;
            // Without F, O can only lie strictly between −2 and 0 here.
            if !(-2 < os && os < 0) {
                return Err(unbalanced_oxidation_states());
            }
            Some((8, os))
        }
        None => return Err(unbalanced_oxidation_states()),
    };
    if let Some((sz, os)) = solved {
        for entry in result.iter_mut().filter(|(z, _)| *z == sz) {
            entry.1 = os;
        }
    }
    Ok(result)
}

// ── Limiting reagent + yield ─────────────────────────────────────────

/// Identify the limiting reagent and return its index.
///
/// - `moles_available`: moles of each reactant available
/// - `stoich_coefficients`: stoichiometric coefficients for each reactant
///
/// The limiting reagent is the one with the smallest ratio moles/coefficient
/// (the first such reagent on a tie).
///
/// # Errors
///
/// Returns error if arrays differ in length or are empty, if any moles value
/// is negative or non-finite, or if any coefficient is not positive and finite.
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
        if !is_non_negative(moles) {
            return Err(KimiyaError::InvalidInput(
                "moles available must be non-negative".into(),
            ));
        }
        if !is_positive(coeff) {
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

/// Theoretical yield of a product, in mol.
///
/// - `moles_available`: moles of each reactant (mol)
/// - `stoich_reactants`: stoichiometric coefficients of reactants
/// - `stoich_product`: stoichiometric coefficient of the desired product
///
/// # Errors
///
/// Returns error if `stoich_product` is not positive and finite, or if
/// [`limiting_reagent`] rejects the reactant data.
pub fn theoretical_yield_moles(
    moles_available: &[f64],
    stoich_reactants: &[f64],
    stoich_product: f64,
) -> Result<f64> {
    if !is_positive(stoich_product) {
        return Err(KimiyaError::InvalidInput(
            "product coefficient must be positive".into(),
        ));
    }
    let idx = limiting_reagent(moles_available, stoich_reactants)?;
    let (Some(&moles), Some(&coeff)) = (moles_available.get(idx), stoich_reactants.get(idx)) else {
        return Err(KimiyaError::ComputationError(
            "limiting reagent index out of range".into(),
        ));
    };
    Ok(moles / coeff * stoich_product)
}

/// Percent yield: (actual / theoretical) × 100, in %.
///
/// `actual` and `theoretical` must be in the same unit (mol or g). Values
/// above 100 % are allowed (impure or wet product).
///
/// # Errors
///
/// Returns error if theoretical yield is not positive and finite, or if the
/// actual yield is negative or non-finite.
#[inline]
pub fn percent_yield(actual: f64, theoretical: f64) -> Result<f64> {
    if !is_positive(theoretical) {
        return Err(KimiyaError::InvalidInput(
            "theoretical yield must be positive".into(),
        ));
    }
    if !is_non_negative(actual) {
        return Err(KimiyaError::InvalidInput(
            "actual yield must be non-negative".into(),
        ));
    }
    Ok(actual / theoretical * 100.0)
}

/// Mass percent composition of each element in a molecule.
///
/// Returns a vector of (atomic_number, mass_percent) pairs, mass percent in %.
///
/// # Errors
///
/// Returns error if molecular weight calculation fails.
pub fn percent_composition(molecule: &crate::molecule::Molecule) -> Result<Vec<(u8, f64)>> {
    let mw = molecule.molecular_weight()?;
    if !is_positive(mw) {
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

/// Maximum deviation of a scaled mole ratio from the nearest integer that
/// [`empirical_formula`] accepts.
const EMPIRICAL_TOLERANCE: f64 = 0.1;

/// Determine empirical formula from mass percentages.
///
/// - `elements`: array of (atomic_number, mass_percent) pairs; the mass
///   percentages need not sum to 100 (only their ratios matter)
///
/// Mole amounts (assuming a 100 g sample) are divided by the smallest, then
/// multiplied by the smallest integer in 1..=6 that brings **every** ratio
/// within 0.1 of an integer.
///
/// Returns array of (atomic_number, count) pairs with integer ratios, in
/// input order.
///
/// # Errors
///
/// Returns error if any element is unknown, if a mass percentage is not
/// positive and finite, or if no multiplier in 1..=6 gives integer ratios
/// within 0.1 ([`KimiyaError::ComputationError`]).
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
        if !is_positive(mass_pct) {
            return Err(KimiyaError::InvalidInput(
                "mass percentages must be positive".into(),
            ));
        }
        moles.push((z, mass_pct / element.atomic_mass));
    }

    // Divide by smallest mole value
    let min_moles = moles.iter().map(|(_, m)| *m).fold(f64::INFINITY, f64::min);
    if !is_positive(min_moles) {
        return Err(KimiyaError::ComputationError(
            "mole calculation produced non-positive value".into(),
        ));
    }
    for (_, m) in &mut moles {
        *m /= min_moles;
        if !is_finite(*m) {
            return Err(KimiyaError::ComputationError(
                "mole ratio is not finite".into(),
            ));
        }
    }

    // Smallest multiplier that makes every ratio close to an integer.
    let mult = (1..=6u32)
        .map(f64::from)
        .find(|&k| {
            moles.iter().all(|&(_, r)| {
                let scaled = r * k;
                (scaled - scaled.round()).abs() <= EMPIRICAL_TOLERANCE
            })
        })
        .ok_or_else(|| {
            KimiyaError::ComputationError(
                "no multiplier in 1..=6 gives integer mole ratios within 0.1".into(),
            )
        })?;

    moles
        .iter()
        .map(|&(z, r)| {
            let count = (r * mult).round();
            if count > f64::from(u32::MAX) {
                return Err(KimiyaError::ComputationError(
                    "empirical formula count exceeds u32".into(),
                ));
            }
            // r ≥ 1 and mult ≥ 1, so 1 ≤ count ≤ u32::MAX: the cast is exact.
            Ok((z, count as u32))
        })
        .collect()
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
    fn limiting_reagent_nan_inf_are_errors() {
        assert!(limiting_reagent(&[f64::NAN, f64::NAN], &[1.0, 1.0]).is_err());
        assert!(limiting_reagent(&[1.0, 1.0], &[f64::NAN, 1.0]).is_err());
        assert!(limiting_reagent(&[1.0, 1.0], &[1.0, f64::INFINITY]).is_err());
        assert!(limiting_reagent(&[f64::INFINITY, 1.0], &[1.0, 1.0]).is_err());
        assert!(limiting_reagent(&[-1.0, 1.0], &[1.0, 1.0]).is_err());
        assert!(limiting_reagent(&[1.0, 1.0], &[0.0, 1.0]).is_err());
    }

    #[test]
    fn limiting_reagent_zero_moles_is_limiting() {
        assert_eq!(limiting_reagent(&[5.0, 0.0], &[1.0, 1.0]).unwrap(), 1);
    }

    #[test]
    fn theoretical_yield_basic() {
        // 2H₂ + O₂ → 2H₂O, have 3 mol H₂ and 2 mol O₂
        // Limiting: H₂ (ratio 1.5), yield = 1.5 × 2 = 3.0 mol H₂O
        let y = theoretical_yield_moles(&[3.0, 2.0], &[2.0, 1.0], 2.0).unwrap();
        assert!((y - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn theoretical_yield_nan_inf_are_errors() {
        assert!(theoretical_yield_moles(&[3.0, 2.0], &[2.0, 1.0], f64::NAN).is_err());
        assert!(theoretical_yield_moles(&[3.0, 2.0], &[2.0, 1.0], f64::INFINITY).is_err());
        assert!(theoretical_yield_moles(&[3.0, 2.0], &[2.0, 1.0], 0.0).is_err());
        assert!(theoretical_yield_moles(&[f64::NAN, 2.0], &[2.0, 1.0], 2.0).is_err());
    }

    #[test]
    fn percent_yield_nan_inf_are_errors() {
        assert!(percent_yield(5.0, f64::NAN).is_err());
        assert!(percent_yield(5.0, f64::INFINITY).is_err());
        assert!(percent_yield(f64::NAN, 10.0).is_err());
        assert!(percent_yield(f64::INFINITY, 10.0).is_err());
        assert!(percent_yield(-1.0, 10.0).is_err());
        assert!(percent_yield(0.0, 10.0).unwrap().abs() < f64::EPSILON);
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
        let formula = empirical_formula(&[(6, 40.00), (1, 6.71), (8, 53.29)]).unwrap();
        assert_eq!(formula, vec![(6, 1), (1, 2), (8, 1)]);
    }

    #[test]
    fn empirical_formula_fe2o3_ratio_one_and_a_half() {
        // Fe 69.94 %, O 30.06 % → O/Fe = 1.50 → ×2 → Fe₂O₃
        let formula = empirical_formula(&[(26, 69.94), (8, 30.06)]).unwrap();
        assert_eq!(formula, vec![(26, 2), (8, 3)]);
    }

    #[test]
    fn empirical_formula_fe3o4_ratio_one_and_a_third() {
        // Fe 72.36 %, O 27.64 % → O/Fe = 1.333 → ×3 → Fe₃O₄
        let formula = empirical_formula(&[(26, 72.36), (8, 27.64)]).unwrap();
        assert_eq!(formula, vec![(26, 3), (8, 4)]);
    }

    #[test]
    fn empirical_formula_close_ratio_uses_multiplier_one() {
        // O/C = 1.08 is within 0.1 of an integer at ×1 → CO.
        // 30 g C = 2.4977 mol; O = 1.08 × that.
        let o_mass = 1.08 * 30.0 / 12.011 * 15.999;
        let formula = empirical_formula(&[(6, 30.0), (8, o_mass)]).unwrap();
        assert_eq!(formula, vec![(6, 1), (8, 1)]);
    }

    #[test]
    fn empirical_formula_no_integer_ratio_is_error() {
        // O/C = 1.43: no multiplier 1..=6 brings it within 0.1 of an integer.
        let o_mass = 1.43 * 30.0 / 12.011 * 15.999;
        let err = empirical_formula(&[(6, 30.0), (8, o_mass)]).unwrap_err();
        assert!(matches!(err, KimiyaError::ComputationError(_)), "{err:?}");
    }

    #[test]
    fn empirical_formula_nan_inf_are_errors() {
        assert!(empirical_formula(&[(6, f64::NAN), (8, 50.0)]).is_err());
        assert!(empirical_formula(&[(6, f64::INFINITY), (8, 50.0)]).is_err());
        assert!(empirical_formula(&[(6, 40.0), (8, -5.0)]).is_err());
        // Subnormal mass percent: ratio overflows to infinity.
        assert!(empirical_formula(&[(6, 1e-320), (8, 50.0)]).is_err());
    }

    #[test]
    fn empirical_formula_empty_is_error() {
        assert!(empirical_formula(&[]).is_err());
    }

    // ── Equation balancing ───────────────────────────────────────────

    /// Σ_j M[i][j]·x[j] for every row (all must be zero).
    fn residuals(comp: &[Vec<f64>], x: &[u32]) -> Vec<f64> {
        comp.iter()
            .map(|row| row.iter().zip(x).map(|(&c, &k)| c * f64::from(k)).sum())
            .collect()
    }

    #[test]
    fn balance_h2_o2_water() {
        // 2H₂ + O₂ → 2H₂O; species [H₂, O₂, H₂O]
        let comp = vec![vec![2.0, 0.0, -2.0], vec![0.0, 2.0, -1.0]];
        assert_eq!(balance_equation(&comp).unwrap(), vec![2, 1, 2]);
    }

    #[test]
    fn balance_ch4_combustion() {
        // CH₄ + 2O₂ → CO₂ + 2H₂O
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
        assert_eq!(coeffs, vec![1, 2, 1, 2]);
        assert!(residuals(&comp, &coeffs).iter().all(|&r| r == 0.0));
    }

    #[test]
    fn balance_ch4_combustion_negated_matrix() {
        // Products positive, reactants negative: same coefficients.
        let comp = vec![
            vec![-1.0, 0.0, 1.0, 0.0],
            vec![-4.0, 0.0, 0.0, 2.0],
            vec![0.0, -2.0, 2.0, 1.0],
        ];
        assert_eq!(balance_equation(&comp).unwrap(), vec![1, 2, 1, 2]);
    }

    #[test]
    fn balance_cuso4_fe_independent_of_row_order() {
        // CuSO₄ + Fe → FeSO₄ + Cu; species [CuSO₄, Fe, FeSO₄, Cu]
        let cu = vec![1.0, 0.0, 0.0, -1.0];
        let s = vec![1.0, 0.0, -1.0, 0.0];
        let o = vec![4.0, 0.0, -4.0, 0.0];
        let fe = vec![0.0, 1.0, -1.0, 0.0];
        let order_a = vec![cu.clone(), s.clone(), o.clone(), fe.clone()];
        let order_b = vec![fe, cu, s, o];
        assert_eq!(balance_equation(&order_a).unwrap(), vec![1, 1, 1, 1]);
        assert_eq!(balance_equation(&order_b).unwrap(), vec![1, 1, 1, 1]);
    }

    #[test]
    fn balance_kmno4_hcl() {
        // 2KMnO₄ + 16HCl → 2KCl + 2MnCl₂ + 8H₂O + 5Cl₂
        //         KMnO4 HCl KCl MnCl2 H2O Cl2
        let comp = vec![
            vec![1.0, 0.0, -1.0, 0.0, 0.0, 0.0],   // K
            vec![1.0, 0.0, 0.0, -1.0, 0.0, 0.0],   // Mn
            vec![4.0, 0.0, 0.0, 0.0, -1.0, 0.0],   // O
            vec![0.0, 1.0, 0.0, 0.0, -2.0, 0.0],   // H
            vec![0.0, 1.0, -1.0, -2.0, 0.0, -2.0], // Cl
        ];
        assert_eq!(balance_equation(&comp).unwrap(), vec![2, 16, 2, 2, 8, 5]);
    }

    #[test]
    fn balance_more_elements_than_species() {
        // NaCl + AgNO₃ → AgCl + NaNO₃: 5 element rows, 4 species
        let comp = vec![
            vec![1.0, 0.0, 0.0, -1.0], // Na
            vec![1.0, 0.0, -1.0, 0.0], // Cl
            vec![0.0, 1.0, -1.0, 0.0], // Ag
            vec![0.0, 1.0, 0.0, -1.0], // N
            vec![0.0, 3.0, 0.0, -3.0], // O
        ];
        assert_eq!(balance_equation(&comp).unwrap(), vec![1, 1, 1, 1]);
    }

    #[test]
    fn balance_mixed_sign_is_error() {
        // H₂O column given with the reactant sign: needs x = [1, 2, 1, -2].
        let comp = vec![
            vec![1.0, 0.0, -1.0, 0.0],
            vec![4.0, 0.0, 0.0, 2.0],
            vec![0.0, 2.0, -2.0, 1.0],
        ];
        assert!(balance_equation(&comp).is_err());
    }

    #[test]
    fn balance_two_independent_reactions_is_error() {
        // H₂ + O₂ → H₂O, H₂O₂: nullity 2 (2H₂+O₂→2H₂O and H₂+O₂→H₂O₂).
        let comp = vec![vec![2.0, 0.0, -2.0, -2.0], vec![0.0, 2.0, -1.0, -2.0]];
        assert!(balance_equation(&comp).is_err());
    }

    #[test]
    fn balance_only_trivial_solution_is_error() {
        // H₂ → O₂ cannot be balanced.
        let comp = vec![vec![2.0, 0.0], vec![0.0, -2.0]];
        assert!(balance_equation(&comp).is_err());
    }

    #[test]
    fn balance_zero_coefficient_is_error() {
        // He takes no part: its coefficient is forced to 0.
        let comp = vec![
            vec![2.0, 0.0, -2.0, 0.0],
            vec![0.0, 2.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        assert!(balance_equation(&comp).is_err());
    }

    #[test]
    fn balance_non_integer_entry_is_error() {
        let err = balance_equation(&[vec![1.5, -1.0]]).unwrap_err();
        assert!(matches!(err, KimiyaError::InvalidInput(_)), "{err:?}");
        assert!(balance_equation(&[vec![f64::NAN, -1.0]]).is_err());
        assert!(balance_equation(&[vec![f64::INFINITY, -1.0]]).is_err());
        assert!(balance_equation(&[vec![1e300, -1.0]]).is_err());
    }

    #[test]
    fn balance_near_integer_entry_accepted() {
        let comp = vec![vec![2.0 + 1e-12, 0.0, -2.0], vec![0.0, 2.0, -1.0 - 1e-12]];
        assert_eq!(balance_equation(&comp).unwrap(), vec![2, 1, 2]);
    }

    #[test]
    fn balance_coefficient_overflowing_u32_is_error() {
        // x = [2³³, 1] does not fit u32.
        assert!(balance_equation(&[vec![1.0, -8_589_934_592.0]]).is_err());
    }

    #[test]
    fn balance_intermediate_overflow_is_error() {
        // Elimination needs (2⁵³ − 1)² > i64::MAX: must be Err, not a panic.
        let big = 9_007_199_254_740_991.0;
        let comp = vec![vec![big, 1.0, -1.0], vec![1.0, big, -1.0]];
        assert!(balance_equation(&comp).is_err());
    }

    #[test]
    fn balance_large_entries_do_not_panic() {
        let big = 9_007_199_254_740_991.0; // 2⁵³ − 1
        let comp = vec![
            vec![big, -(big - 2.0), 3.0, -7.0],
            vec![big - 4.0, 5.0, -big, 11.0],
            vec![13.0, big - 6.0, -17.0, -big],
        ];
        // Whatever the outcome, it must be an Err or a verified balance.
        if let Ok(x) = balance_equation(&comp) {
            assert!(x.iter().all(|&k| k > 0));
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

    fn os_of(os: &[(u8, i32)], z: u8) -> i32 {
        os.iter().find(|&&(e, _)| e == z).unwrap().1
    }

    #[test]
    fn oxidation_state_water() {
        // H₂O: H = +1, O = -2
        let os = assign_oxidation_states(&[(1, 2), (8, 1)], 0).unwrap();
        assert_eq!(os, vec![(1, 1), (8, -2)]);
    }

    #[test]
    fn oxidation_state_nacl() {
        // NaCl: Na = +1, Cl = -1
        let os = assign_oxidation_states(&[(11, 1), (17, 1)], 0).unwrap();
        assert_eq!(os, vec![(11, 1), (17, -1)]);
    }

    #[test]
    fn oxidation_state_sulfate_ion() {
        // SO₄²⁻: O = -2, S = +6
        let os = assign_oxidation_states(&[(16, 1), (8, 4)], -2).unwrap();
        assert_eq!(os, vec![(16, 6), (8, -2)]);
    }

    #[test]
    fn oxidation_state_permanganate() {
        // MnO₄⁻: O = -2, Mn = +7
        let os = assign_oxidation_states(&[(25, 1), (8, 4)], -1).unwrap();
        assert_eq!(os_of(&os, 25), 7);
    }

    #[test]
    fn oxidation_state_kmno4() {
        let os = assign_oxidation_states(&[(19, 1), (25, 1), (8, 4)], 0).unwrap();
        assert_eq!(os, vec![(19, 1), (25, 7), (8, -2)]);
    }

    #[test]
    fn oxidation_state_dichromate() {
        // Cr₂O₇²⁻: Cr = +6
        let os = assign_oxidation_states(&[(24, 2), (8, 7)], -2).unwrap();
        assert_eq!(os, vec![(24, 6), (8, -2)]);
    }

    #[test]
    fn oxidation_state_ammonium() {
        // NH₄⁺: N = -3
        let os = assign_oxidation_states(&[(7, 1), (1, 4)], 1).unwrap();
        assert_eq!(os, vec![(7, -3), (1, 1)]);
    }

    #[test]
    fn oxidation_state_co2() {
        // CO₂: O = -2, C = +4
        let os = assign_oxidation_states(&[(6, 1), (8, 2)], 0).unwrap();
        assert_eq!(os_of(&os, 6), 4);
    }

    #[test]
    fn oxidation_state_hydrogen_peroxide() {
        let os = assign_oxidation_states(&[(1, 2), (8, 2)], 0).unwrap();
        assert_eq!(os, vec![(1, 1), (8, -1)]);
    }

    #[test]
    fn oxidation_state_sodium_peroxide() {
        let os = assign_oxidation_states(&[(11, 2), (8, 2)], 0).unwrap();
        assert_eq!(os, vec![(11, 1), (8, -1)]);
    }

    #[test]
    fn oxidation_state_superoxide_is_error() {
        // KO₂: O = -1/2 is not an integer.
        assert!(assign_oxidation_states(&[(19, 1), (8, 2)], 0).is_err());
    }

    #[test]
    fn oxidation_state_oxygen_difluoride() {
        // OF₂: F = -1, O = +2
        let os = assign_oxidation_states(&[(8, 1), (9, 2)], 0).unwrap();
        assert_eq!(os, vec![(8, 2), (9, -1)]);
    }

    #[test]
    fn oxidation_state_sodium_hydride() {
        let os = assign_oxidation_states(&[(11, 1), (1, 1)], 0).unwrap();
        assert_eq!(os, vec![(11, 1), (1, -1)]);
    }

    #[test]
    fn oxidation_state_calcium_hydride() {
        let os = assign_oxidation_states(&[(20, 1), (1, 2)], 0).unwrap();
        assert_eq!(os, vec![(20, 2), (1, -1)]);
    }

    #[test]
    fn oxidation_state_lithium_aluminium_hydride() {
        // LiAlH₄: Li +1, Al +3, H -1
        let os = assign_oxidation_states(&[(3, 1), (13, 1), (1, 4)], 0).unwrap();
        assert_eq!(os, vec![(3, 1), (13, 3), (1, -1)]);
    }

    #[test]
    fn oxidation_state_single_element() {
        assert_eq!(assign_oxidation_states(&[(8, 2)], 0).unwrap(), vec![(8, 0)]);
        assert_eq!(assign_oxidation_states(&[(9, 2)], 0).unwrap(), vec![(9, 0)]);
        assert_eq!(
            assign_oxidation_states(&[(17, 1)], -1).unwrap(),
            vec![(17, -1)]
        );
        // Peroxide ion O₂²⁻
        assert_eq!(
            assign_oxidation_states(&[(8, 2)], -2).unwrap(),
            vec![(8, -1)]
        );
        // Superoxide ion O₂⁻: -1/2
        assert!(assign_oxidation_states(&[(8, 2)], -1).is_err());
    }

    #[test]
    fn oxidation_state_propane_non_integer_is_error() {
        // C₃H₈: C = -8/3
        let err = assign_oxidation_states(&[(6, 3), (1, 8)], 0).unwrap_err();
        assert!(matches!(err, KimiyaError::ComputationError(_)), "{err:?}");
    }

    #[test]
    fn oxidation_state_magnetite_non_integer_is_error() {
        // Fe₃O₄: Fe = +8/3
        assert!(assign_oxidation_states(&[(26, 3), (8, 4)], 0).is_err());
    }

    #[test]
    fn oxidation_state_repeated_unknown_element_is_merged() {
        // C₂H₆ written with two C entries: C = -3 for both.
        let os = assign_oxidation_states(&[(6, 1), (1, 6), (6, 1)], 0).unwrap();
        assert_eq!(os, vec![(6, -3), (1, 1), (6, -3)]);
    }

    #[test]
    fn oxidation_state_unbalanced_assigned_is_error() {
        // NaF with charge +1: every atom is rule-assigned, no O to re-solve.
        assert!(assign_oxidation_states(&[(11, 1), (9, 1)], 1).is_err());
        // "H₂O⁻": re-solving O would give −3, below O's minimum of −2.
        assert!(assign_oxidation_states(&[(1, 2), (8, 1)], -1).is_err());
        // H₂O⁺ (radical cation): re-solved O = −1 is within (−2, 0), accepted.
        let os = assign_oxidation_states(&[(1, 2), (8, 1)], 1).unwrap();
        assert_eq!(os, vec![(1, 1), (8, -1)]);
    }

    #[test]
    fn oxidation_state_hydroxide_and_hydroperoxide() {
        // OH⁻: O −2, H +1 balances directly.
        assert_eq!(
            assign_oxidation_states(&[(8, 1), (1, 1)], -1).unwrap(),
            vec![(8, -2), (1, 1)]
        );
        // HO₂⁻ (hydroperoxide): O −1.
        assert_eq!(
            assign_oxidation_states(&[(1, 1), (8, 2)], -1).unwrap(),
            vec![(1, 1), (8, -1)]
        );
    }

    #[test]
    fn oxidation_state_two_unknowns_is_error() {
        // CN⁻: C and N both unknown.
        assert!(assign_oxidation_states(&[(6, 1), (7, 1)], -1).is_err());
    }

    #[test]
    fn oxidation_state_huge_count_does_not_panic() {
        let big = 1u32 << 31;
        assert!(assign_oxidation_states(&[(25, 1), (8, big)], 0).is_err());
        assert!(assign_oxidation_states(&[(25, 1), (8, u32::MAX)], 0).is_err());
        assert!(assign_oxidation_states(&[(11, u32::MAX), (17, u32::MAX)], 0).is_ok());
        assert!(assign_oxidation_states(&[(8, u32::MAX), (8, u32::MAX)], i32::MIN).is_err());
    }

    #[test]
    fn oxidation_state_zero_count_or_unknown_element_is_error() {
        assert!(assign_oxidation_states(&[(25, 0), (8, 4)], -1).is_err());
        assert!(assign_oxidation_states(&[(200, 1), (8, 1)], 0).is_err());
    }

    #[test]
    fn oxidation_state_empty_is_error() {
        assert!(assign_oxidation_states(&[], 0).is_err());
    }
}
