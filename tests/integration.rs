use kimiya::*;
use kimiya::molecule::Molecule;
use kimiya::gas;
use kimiya::kinetics;
use kimiya::solution;

#[test]
fn water_weight_from_elements() {
    let water = Molecule::water();
    let mw = water.molecular_weight();
    assert!((mw - 18.015).abs() < 0.01);
}

#[test]
fn ideal_gas_at_stp() {
    let v = gas::ideal_gas_volume(1.0, gas::STP_TEMPERATURE, gas::STP_PRESSURE);
    let v_liters = v * 1000.0;
    assert!((v_liters - 22.414).abs() < 0.1);
}

#[test]
fn pure_water_is_neutral() {
    let ph = solution::ph_from_h_concentration(1e-7);
    assert!((ph - 7.0).abs() < 0.01);
}

#[test]
fn arrhenius_temperature_dependence() {
    let k_low = kinetics::arrhenius_rate(1e13, 50_000.0, 300.0);
    let k_high = kinetics::arrhenius_rate(1e13, 50_000.0, 500.0);
    assert!(k_high > k_low);
}
