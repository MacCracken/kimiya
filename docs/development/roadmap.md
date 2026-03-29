# Kimiya Roadmap

## Status

**v1.0.0** — Stable release. Complete chemistry engine. 21 modules, 118 elements, 444 tests, hisab 1.1 integrated. API frozen.

---

## Cross-Crate Bridges

- [ ] **`bridge.rs` module** — primitive-value conversions for cross-crate chemistry
- [ ] **ushma bridge**: reaction enthalpy (J/mol) → heat source/sink; activation energy (J/mol) → temperature-dependent rate
- [ ] **tanmatra bridge**: atomic number → electron configuration; ionization energy (eV) → bond energy estimation
- [ ] **bijli bridge**: electrochemical potential (V) → current density; ionic concentration → conductivity

---

## Soorat Integration (rendering visualization)

- [ ] **`integration/soorat.rs` module** — feature-gated `soorat-compat`
- [ ] **Molecular structure**: atom positions (3D), bond connectivity, element types for ball-and-stick rendering
- [ ] **Reaction network graph**: species nodes, reaction edges, rate constants for node-link visualization
- [ ] **Spectroscopy data**: wavelength vs absorption/emission intensity for spectral plot rendering
- [ ] **Phase diagram**: phase boundary curves, triple/critical points for contour rendering
