use kimiya::{electrochemistry, element, gas, kinetics, molecule::Molecule, solution, thermochem};

fn main() {
    // Look up Iron
    let fe = element::lookup_by_symbol("Fe").unwrap();
    println!(
        "{}: Z={}, M={:.3} g/mol",
        fe.name, fe.atomic_number, fe.atomic_mass
    );

    // Water molecular weight
    let water = Molecule::water();
    println!(
        "H2O molecular weight: {:.3} g/mol",
        water.molecular_weight().unwrap()
    );

    // Ideal gas at STP
    let v = gas::ideal_gas_volume(1.0, gas::STP_TEMPERATURE, gas::STP_PRESSURE).unwrap();
    println!("1 mol ideal gas at STP: {:.1} L", v * 1000.0);

    // pH of 0.01M HCl
    let ph = solution::ph_from_h_concentration(0.01).unwrap();
    println!("pH of 0.01M HCl: {:.1}", ph);

    // Arrhenius rate
    let k = kinetics::arrhenius_rate(1e13, 75_000.0, 298.0).unwrap();
    println!("Rate constant at 298K: {k:.4e} s⁻¹");

    // Daniell cell (electrochemistry)
    let e_cell = electrochemistry::cell_potential_from_couples("Cu2+/Cu", "Zn2+/Zn").unwrap();
    println!("Daniell cell E°: {e_cell:.3} V");

    // Nernst equation — Cu²⁺/Cu at non-standard [Cu²⁺] = 0.01 M
    let e = electrochemistry::nernst_potential_25c(0.342, 2, 0.01).unwrap();
    println!("Cu²⁺/Cu at [Cu²⁺]=0.01M: {e:.4} V");

    // Methane combustion enthalpy (thermochem)
    let dh = thermochem::reaction_enthalpy(
        &[("CO2(g)", 1.0), ("H2O(l)", 2.0)],
        &[("CH4(g)", 1.0), ("O2(g)", 2.0)],
    )
    .unwrap();
    println!("CH₄ combustion ΔH°: {dh:.1} kJ/mol");
}
