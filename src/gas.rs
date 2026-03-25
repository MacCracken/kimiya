use crate::element::GAS_CONSTANT;

/// Ideal gas law: P = nRT/V
///
/// n = moles, T = temperature (K), V = volume (m³). Returns pressure in Pascals.
#[must_use]
#[inline]
pub fn ideal_gas_pressure(moles: f64, temperature_k: f64, volume_m3: f64) -> f64 {
    if volume_m3 <= 0.0 { return 0.0; }
    moles * GAS_CONSTANT * temperature_k / volume_m3
}

/// Ideal gas law: V = nRT/P
#[must_use]
#[inline]
pub fn ideal_gas_volume(moles: f64, temperature_k: f64, pressure_pa: f64) -> f64 {
    if pressure_pa <= 0.0 { return 0.0; }
    moles * GAS_CONSTANT * temperature_k / pressure_pa
}

/// Van der Waals equation: (P + a·n²/V²)(V - n·b) = nRT
///
/// Returns pressure in Pascals. a and b are Van der Waals constants.
#[must_use]
pub fn van_der_waals_pressure(moles: f64, temperature_k: f64, volume_m3: f64, a: f64, b: f64) -> f64 {
    if volume_m3 <= moles * b { return 0.0; }
    let n_over_v = moles / volume_m3;
    (moles * GAS_CONSTANT * temperature_k / (volume_m3 - moles * b)) - a * n_over_v * n_over_v
}

/// Dalton's law: partial pressure = total pressure × mole fraction.
#[must_use]
#[inline]
pub fn partial_pressure(total_pressure: f64, mole_fraction: f64) -> f64 {
    total_pressure * mole_fraction
}

/// Gas density: ρ = PM / (RT)
///
/// P in Pa, M in kg/mol, T in K. Returns density in kg/m³.
#[must_use]
#[inline]
pub fn gas_density(pressure_pa: f64, molar_mass_kg: f64, temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 { return 0.0; }
    pressure_pa * molar_mass_kg / (GAS_CONSTANT * temperature_k)
}

/// STP (Standard Temperature and Pressure): 273.15 K, 101325 Pa.
pub const STP_TEMPERATURE: f64 = 273.15;
pub const STP_PRESSURE: f64 = 101325.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ideal_gas_stp_volume() {
        // 1 mol at STP should occupy ~22.4 L = 0.0224 m³
        let v = ideal_gas_volume(1.0, STP_TEMPERATURE, STP_PRESSURE);
        let v_liters = v * 1000.0;
        assert!((v_liters - 22.414).abs() < 0.1, "1 mol at STP should be ~22.4 L, got {v_liters}");
    }

    #[test]
    fn ideal_gas_roundtrip() {
        let p = ideal_gas_pressure(2.0, 300.0, 0.05);
        let v = ideal_gas_volume(2.0, 300.0, p);
        assert!((v - 0.05).abs() < 1e-10, "roundtrip should recover volume");
    }

    #[test]
    fn van_der_waals_approaches_ideal() {
        // With a=0, b=0, Van der Waals should equal ideal gas
        let p_ideal = ideal_gas_pressure(1.0, 300.0, 0.025);
        let p_vdw = van_der_waals_pressure(1.0, 300.0, 0.025, 0.0, 0.0);
        assert!((p_ideal - p_vdw).abs() < 1.0, "a=0 b=0 should match ideal gas");
    }

    #[test]
    fn partial_pressure_basic() {
        let pp = partial_pressure(101325.0, 0.21); // O2 in air
        assert!((pp - 21278.25).abs() < 1.0);
    }

    #[test]
    fn gas_density_air_at_stp() {
        // Air (M ≈ 0.029 kg/mol) at STP → ρ ≈ 1.29 kg/m³
        let rho = gas_density(STP_PRESSURE, 0.029, STP_TEMPERATURE);
        assert!((rho - 1.29).abs() < 0.05, "air density at STP should be ~1.29, got {rho}");
    }

    #[test]
    fn zero_volume_returns_zero() {
        assert_eq!(ideal_gas_pressure(1.0, 300.0, 0.0), 0.0);
    }
}
