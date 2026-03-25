use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_molecular_weight(c: &mut Criterion) {
    let glucose = kimiya::molecule::Molecule::glucose();
    c.bench_function("molecule/molecular_weight_glucose", |b| {
        b.iter(|| black_box(&glucose).molecular_weight().unwrap());
    });
}

fn bench_ideal_gas(c: &mut Criterion) {
    c.bench_function("gas/ideal_gas_pressure", |b| {
        b.iter(|| {
            kimiya::gas::ideal_gas_pressure(black_box(1.0), black_box(300.0), black_box(0.025))
                .unwrap()
        });
    });
}

fn bench_arrhenius(c: &mut Criterion) {
    c.bench_function("kinetics/arrhenius_rate", |b| {
        b.iter(|| {
            kimiya::kinetics::arrhenius_rate(black_box(1e13), black_box(50_000.0), black_box(300.0))
                .unwrap()
        });
    });
}

fn bench_ph(c: &mut Criterion) {
    c.bench_function("solution/ph_from_h", |b| {
        b.iter(|| kimiya::solution::ph_from_h_concentration(black_box(1e-4)).unwrap());
    });
}

fn bench_element_lookup(c: &mut Criterion) {
    c.bench_function("element/lookup_by_symbol", |b| {
        b.iter(|| kimiya::element::lookup_by_symbol(black_box("Fe")));
    });
}

fn bench_element_lookup_by_number(c: &mut Criterion) {
    c.bench_function("element/lookup_by_number", |b| {
        b.iter(|| kimiya::element::lookup_by_number(black_box(26)));
    });
}

fn bench_gibbs(c: &mut Criterion) {
    c.bench_function("reaction/gibbs_free_energy", |b| {
        b.iter(|| {
            kimiya::reaction::gibbs_free_energy(
                black_box(-100_000.0),
                black_box(298.0),
                black_box(-200.0),
            )
            .unwrap()
        });
    });
}

fn bench_heat_transfer(c: &mut Criterion) {
    c.bench_function("thermo/heat_transfer", |b| {
        b.iter(|| {
            kimiya::thermo::heat_transfer(black_box(1000.0), black_box(4.184), black_box(10.0))
        });
    });
}

fn bench_nernst(c: &mut Criterion) {
    c.bench_function("electrochemistry/nernst_potential", |b| {
        b.iter(|| {
            kimiya::electrochemistry::nernst_potential(
                black_box(0.342),
                black_box(2),
                black_box(298.15),
                black_box(0.1),
            )
            .unwrap()
        });
    });
}

fn bench_half_reaction_lookup(c: &mut Criterion) {
    c.bench_function("electrochemistry/lookup_half_reaction", |b| {
        b.iter(|| kimiya::electrochemistry::lookup_half_reaction(black_box("Cu2+/Cu")));
    });
}

fn bench_reaction_enthalpy(c: &mut Criterion) {
    c.bench_function("thermochem/reaction_enthalpy_ch4", |b| {
        b.iter(|| {
            kimiya::thermochem::reaction_enthalpy(
                black_box(&[("CO2(g)", 1.0), ("H2O(l)", 2.0)]),
                black_box(&[("CH4(g)", 1.0), ("O2(g)", 2.0)]),
            )
            .unwrap()
        });
    });
}

fn bench_vant_hoff(c: &mut Criterion) {
    c.bench_function("thermochem/vant_hoff_k", |b| {
        b.iter(|| {
            kimiya::thermochem::vant_hoff_k(
                black_box(100.0),
                black_box(-50_000.0),
                black_box(298.15),
                black_box(400.0),
            )
            .unwrap()
        });
    });
}

fn bench_weak_acid_ph(c: &mut Criterion) {
    c.bench_function("solution/weak_acid_ph", |b| {
        b.iter(|| kimiya::solution::weak_acid_ph(black_box(1.8e-5), black_box(0.1)).unwrap());
    });
}

fn bench_enthalpy_change_cp(c: &mut Criterion) {
    c.bench_function("thermochem/enthalpy_change_cp", |b| {
        b.iter(|| {
            kimiya::thermochem::enthalpy_change_cp(
                black_box("CO2(g)"),
                black_box(300.0),
                black_box(600.0),
            )
            .unwrap()
        });
    });
}

fn bench_shomate_cp(c: &mut Criterion) {
    let s = kimiya::thermochem::lookup_shomate("CO2(g)").unwrap();
    c.bench_function("thermochem/shomate_cp", |b| {
        b.iter(|| black_box(s).cp(black_box(500.0)));
    });
}

fn bench_bohr_energy(c: &mut Criterion) {
    c.bench_function("spectroscopy/bohr_energy_level", |b| {
        b.iter(|| kimiya::spectroscopy::bohr_energy_level(black_box(1), black_box(3)).unwrap());
    });
}

fn bench_beer_lambert(c: &mut Criterion) {
    c.bench_function("spectroscopy/absorbance", |b| {
        b.iter(|| {
            kimiya::spectroscopy::absorbance(black_box(100.0), black_box(1.0), black_box(0.01))
        });
    });
}

fn bench_michaelis_menten(c: &mut Criterion) {
    c.bench_function("kinetics/michaelis_menten", |b| {
        b.iter(|| {
            kimiya::kinetics::michaelis_menten(black_box(100.0), black_box(5.0), black_box(10.0))
                .unwrap()
        });
    });
}

fn bench_eyring(c: &mut Criterion) {
    c.bench_function("kinetics/eyring_rate", |b| {
        b.iter(|| kimiya::kinetics::eyring_rate(black_box(80_000.0), black_box(298.0)).unwrap());
    });
}

fn bench_simulate_first_order(c: &mut Criterion) {
    c.bench_function("dynamics/first_order_100steps", |b| {
        b.iter(|| {
            kimiya::reaction_dynamics::simulate_first_order(
                black_box(1.0),
                black_box(0.1),
                black_box(10.0),
                black_box(100),
            )
            .unwrap()
        });
    });
}

fn bench_simulate_consecutive(c: &mut Criterion) {
    c.bench_function("dynamics/consecutive_100steps", |b| {
        b.iter(|| {
            kimiya::reaction_dynamics::simulate_consecutive(
                black_box(1.0),
                black_box(0.5),
                black_box(0.1),
                black_box(20.0),
                black_box(100),
            )
            .unwrap()
        });
    });
}

fn bench_fit_arrhenius(c: &mut Criterion) {
    let temps: [f64; 6] = [300.0, 320.0, 340.0, 360.0, 380.0, 400.0];
    let r: f64 = 8.314462618;
    let rates: Vec<f64> = temps
        .iter()
        .map(|&t| 1e10_f64 * (-50_000.0 / (r * t)).exp())
        .collect();
    c.bench_function("fitting/arrhenius_6pts", |b| {
        b.iter(|| kimiya::fitting::fit_arrhenius(black_box(&temps), black_box(&rates)).unwrap());
    });
}

fn bench_fit_polynomial(c: &mut Criterion) {
    let x: Vec<f64> = (0..20).map(|i| i as f64 * 0.1).collect();
    let y: Vec<f64> = x.iter().map(|&xi| 1.0 + 2.0 * xi + 0.5 * xi * xi).collect();
    c.bench_function("fitting/polynomial_deg2_20pts", |b| {
        b.iter(|| {
            kimiya::fitting::fit_polynomial(black_box(&x), black_box(&y), black_box(2)).unwrap()
        });
    });
}

criterion_group!(
    benches,
    bench_molecular_weight,
    bench_ideal_gas,
    bench_arrhenius,
    bench_ph,
    bench_element_lookup,
    bench_element_lookup_by_number,
    bench_gibbs,
    bench_heat_transfer,
    bench_nernst,
    bench_half_reaction_lookup,
    bench_reaction_enthalpy,
    bench_vant_hoff,
    bench_weak_acid_ph,
    bench_enthalpy_change_cp,
    bench_shomate_cp,
    bench_bohr_energy,
    bench_beer_lambert,
    bench_michaelis_menten,
    bench_eyring,
    bench_simulate_first_order,
    bench_simulate_consecutive,
    bench_fit_arrhenius,
    bench_fit_polynomial
);
criterion_main!(benches);
