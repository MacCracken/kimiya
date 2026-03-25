use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_molecular_weight(c: &mut Criterion) {
    let glucose = kimiya::molecule::Molecule::glucose();
    c.bench_function("molecule/molecular_weight_glucose", |b| {
        b.iter(|| black_box(&glucose).molecular_weight());
    });
}

fn bench_ideal_gas(c: &mut Criterion) {
    c.bench_function("gas/ideal_gas_pressure", |b| {
        b.iter(|| kimiya::gas::ideal_gas_pressure(black_box(1.0), black_box(300.0), black_box(0.025)));
    });
}

fn bench_arrhenius(c: &mut Criterion) {
    c.bench_function("kinetics/arrhenius_rate", |b| {
        b.iter(|| kimiya::kinetics::arrhenius_rate(black_box(1e13), black_box(50_000.0), black_box(300.0)));
    });
}

fn bench_ph(c: &mut Criterion) {
    c.bench_function("solution/ph_from_h", |b| {
        b.iter(|| kimiya::solution::ph_from_h_concentration(black_box(1e-4)));
    });
}

fn bench_element_lookup(c: &mut Criterion) {
    c.bench_function("element/lookup_by_symbol", |b| {
        b.iter(|| kimiya::element::lookup_by_symbol(black_box("Fe")));
    });
}

criterion_group!(benches, bench_molecular_weight, bench_ideal_gas, bench_arrhenius, bench_ph, bench_element_lookup);
criterion_main!(benches);
