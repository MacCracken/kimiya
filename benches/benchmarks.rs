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

criterion_group!(
    benches,
    bench_molecular_weight,
    bench_ideal_gas,
    bench_arrhenius,
    bench_ph,
    bench_element_lookup,
    bench_element_lookup_by_number,
    bench_gibbs,
    bench_heat_transfer
);
criterion_main!(benches);
