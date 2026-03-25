use kimiya::electrochemistry;
use kimiya::gas;
use kimiya::kinetics;
use kimiya::molecule::Molecule;
use kimiya::solution;
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
