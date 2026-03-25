use kimiya::{element, gas, kinetics, molecule::Molecule, solution};

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
}
