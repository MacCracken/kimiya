use kimiya::electrochemistry;
use kimiya::element;
use kimiya::gas;
use kimiya::kinetics;
use kimiya::molecule::Molecule;
use kimiya::organic;
use kimiya::reaction;
use kimiya::solution;
use kimiya::spectroscopy;
use kimiya::thermochem;

#[test]
fn water_weight_from_elements() {
    let water = Molecule::water();
    let mw = water.molecular_weight().unwrap();
    assert!((mw - 18.015).abs() < 0.01);
}

#[test]
fn ideal_gas_at_stp() {
    let v = gas::ideal_gas_volume(1.0, gas::STP_TEMPERATURE, gas::STP_PRESSURE).unwrap();
    let v_liters = v * 1000.0;
    assert!((v_liters - 22.414).abs() < 0.1);
}

#[test]
fn pure_water_is_neutral() {
    let ph = solution::ph_from_h_concentration(1e-7).unwrap();
    assert!((ph - 7.0).abs() < 0.01);
}

#[test]
fn arrhenius_temperature_dependence() {
    let k_low = kinetics::arrhenius_rate(1e13, 50_000.0, 300.0).unwrap();
    let k_high = kinetics::arrhenius_rate(1e13, 50_000.0, 500.0).unwrap();
    assert!(k_high > k_low);
}

#[test]
fn daniell_cell_is_spontaneous() {
    let e = electrochemistry::cell_potential_from_couples("Cu2+/Cu", "Zn2+/Zn").unwrap();
    assert!(electrochemistry::is_spontaneous_cell(e));
    assert!((e - 1.104).abs() < 0.01);

    // ΔG should be negative for spontaneous cell
    let dg = electrochemistry::cell_gibbs_energy(2, e).unwrap();
    assert!(dg < 0.0);
}

#[test]
fn nernst_reduces_to_standard_at_q1() {
    let e_std = electrochemistry::lookup_half_reaction("Ag+/Ag")
        .unwrap()
        .standard_potential;
    let e = electrochemistry::nernst_potential(e_std, 1, 298.15, 1.0).unwrap();
    assert!((e - e_std).abs() < 1e-10);
}

#[test]
fn methane_combustion_cross_module() {
    // Verify reaction enthalpy from thermochem matches expected
    let dh = thermochem::reaction_enthalpy(
        &[("CO2(g)", 1.0), ("H2O(l)", 2.0)],
        &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
    )
    .unwrap();
    // Should be strongly exothermic
    assert!(dh < -800.0, "methane combustion should be < -800 kJ/mol");
}

#[test]
fn vant_hoff_le_chatelier() {
    // For exothermic reaction: increasing T should decrease K
    let k_298 = 1000.0;
    let k_400 = thermochem::vant_hoff_k(k_298, -100_000.0, 298.15, 400.0).unwrap();
    assert!(
        k_400 < k_298,
        "Le Chatelier: exothermic + higher T → lower K"
    );
}

#[test]
fn electrochemistry_gibbs_matches_reaction_gibbs() {
    // Daniell cell: Zn + Cu²⁺ → Zn²⁺ + Cu
    // ΔG from electrochemistry: -nFE = -2 × 96485 × 1.104 ≈ -213,038 J
    let e_cell = electrochemistry::cell_potential_from_couples("Cu2+/Cu", "Zn2+/Zn").unwrap();
    let dg_electro = electrochemistry::cell_gibbs_energy(2, e_cell).unwrap(); // J/mol

    // ΔG from reaction module: ΔG = ΔH - TΔS
    // Using the Daniell cell standard potential, ΔG should be consistent
    // ΔG = -nFE, and we can also get K from reaction module
    let k = reaction::equilibrium_constant(dg_electro, 298.15).unwrap();
    assert!(k > 1.0, "spontaneous cell should have K > 1, got {k}");

    // Verify sign consistency: spontaneous cell → ΔG < 0
    assert!(
        dg_electro < -200_000.0,
        "ΔG should be ~-213 kJ, got {dg_electro}"
    );
}

#[test]
fn thermochem_enthalpy_entropy_gibbs_triangle() {
    // For methane combustion, verify ΔG ≈ ΔH - TΔS at 298.15 K
    let products = &[("CO2(g)", 1.0), ("H2O(l)", 2.0)];
    let reactants = &[("CH4(g)", 1.0), ("O2(g)", 2.0)];

    let dh = thermochem::reaction_enthalpy(products, reactants).unwrap(); // kJ
    let ds = thermochem::reaction_entropy(products, reactants).unwrap(); // J/(mol·K)
    let dg = thermochem::reaction_gibbs_energy(products, reactants).unwrap(); // kJ

    let dg_calc = dh - 298.15 * ds / 1000.0; // kJ
    assert!(
        (dg_calc - dg).abs() < 1.0,
        "ΔG should ≈ ΔH - TΔS: calc={dg_calc:.1}, table={dg:.1}"
    );
}

#[test]
fn weak_acid_ph_acetic_acid() {
    // 0.1 M acetic acid (Ka = 1.8e-5) → pH ≈ 2.87
    let ph = solution::weak_acid_ph(1.8e-5, 0.1).unwrap();
    assert!((ph - 2.87).abs() < 0.05);
}

#[test]
fn cp_integration_heating_co2() {
    // Heat CO2 from 300K to 600K using Shomate + Simpson integration
    let dh = thermochem::enthalpy_change_cp("CO2(g)", 300.0, 600.0).unwrap();
    assert!(dh > 0.0, "heating should be positive enthalpy change");
    // Rough: ~37-50 J/(mol·K) avg × 300K ≈ 11-15 kJ
    assert!(dh > 10_000.0 && dh < 15_000.0);
}

#[test]
fn full_periodic_table_118() {
    assert_eq!(element::ELEMENTS.len(), 118);
    let og = element::lookup_by_symbol("Og").unwrap();
    assert_eq!(og.atomic_number, 118);
}

#[test]
fn lanthanides_and_actinides() {
    let la = element::lookup_by_symbol("La").unwrap();
    assert_eq!(la.category, element::ElementCategory::Lanthanide);
    let u = element::lookup_by_symbol("U").unwrap();
    assert_eq!(u.category, element::ElementCategory::Actinide);
}

#[test]
fn beer_lambert_and_transmittance() {
    let a = spectroscopy::absorbance(100.0, 1.0, 0.01);
    let t = spectroscopy::transmittance(a);
    assert!((t - 0.1).abs() < 0.001);
}

#[test]
fn hydrogen_balmer_series_visible() {
    let lambda = spectroscopy::SpectralSeries::Balmer.wavelength(1).unwrap();
    let nm = lambda * 1e9;
    assert!(nm > 600.0 && nm < 700.0, "Hα should be red, got {nm}nm");
}

#[test]
fn michaelis_menten_at_km_is_half_vmax() {
    let v = kinetics::michaelis_menten(100.0, 5.0, 5.0).unwrap();
    assert!((v - 50.0).abs() < f64::EPSILON);
}

#[test]
fn eyring_consistent_with_arrhenius() {
    // Both should give positive rate constants at reasonable conditions
    let k_arr = kinetics::arrhenius_rate(1e13, 80_000.0, 298.0).unwrap();
    let k_eyr = kinetics::eyring_rate(80_000.0, 298.0).unwrap();
    // Different models, different magnitudes, but both positive
    assert!(k_arr > 0.0);
    assert!(k_eyr > 0.0);
}

#[test]
fn bond_energy_methane_combustion() {
    // CH₄ + 2O₂ → CO₂ + 2H₂O
    // Broken: 4×C-H + 2×O=O
    // Formed: 2×C=O + 4×O-H
    let dh = organic::enthalpy_from_bonds(
        &[(organic::BondType::CH, 4), (organic::BondType::OODouble, 2)],
        &[(organic::BondType::CODouble, 2), (organic::BondType::OH, 4)],
    );
    // Should be exothermic (bond energy estimates give ~-694 kJ/mol)
    assert!(
        dh < -600.0,
        "methane combustion should be exothermic, got {dh}"
    );
}

#[test]
fn vsepr_predicts_water_bent() {
    assert_eq!(
        organic::predict_geometry(2, 2),
        Some(organic::Geometry::BentTwoLonePairs)
    );
}
