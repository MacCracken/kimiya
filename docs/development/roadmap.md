# Kimiya Roadmap

## Status

**v2.0.0** (2026-11-13) is the Cyrius port of Rust 1.2.0, done to the plan of record
[`kimiya-port-plan.md`](kimiya-port-plan.md) (M0–M7). What it delivered is in the [CHANGELOG](../../CHANGELOG.md), the
Rust tests and public surface it accounts for are in [`rust-parity.md`](rust-parity.md), and live status is in
[`state.md`](state.md). The Rust crate stays at **v1.2.0**, frozen in `rust-old/` as the parity oracle; there is no
Rust 2.0.

Every open item is below, grouped by the release line it belongs to.

---

## Backlog

### 2.0.x — remaining port and audit

Each item that changes a result or a signature is a Cyrius API change for ADR 0002 when it lands.

**The library-wide overflow audit, resumed.** It was paused after B10a (the user's decision of 2026-11-12) and now runs
once over the whole library, M5's modules included. The batches use the survey's sixteen families and numbering from
`state.md`. The rules come from the user's decisions of 2026-11-09 ("Library-wide audit") and 2026-11-10 ("Round once"):
a single product or quotient is rounded once, a sum inside it gets #20's edge fix, and each family gets one ADR 0002
row. The done batches, B1–B7 and B10a, are ADR 0002 #27–#39 and #40's primitives.

- [ ] B10, F10: the unbounded exponentials and the exponential × prefactor functions. The primitives `_km_cr_exp_ux`,
  `_km_cr_expm1_ux` and `_km_cr_exp10_ux` are done (B10a). #40's functions remain, in three sub-batches:
  - kinetics: Arrhenius, collision, Eyring ×2, first-order;
  - nuclear: `atoms_remaining`, the decay chain;
  - van 't Hoff, Clausius–Clapeyron, Antoine, Morse and `fugacity_pr`.
- [ ] B11, F11 and F12: `nth_order_half_life`'s pow and `henderson_hasselbalch`'s log10.
- [ ] B12: the iterative functions, by the user's answer "Fix all three". This extends ADR 0002 #38, whose flame
  sign test was done in B8, with:
  - the weak acid/base pH Newton, with an unbounded exponent on a slow path only;
  - Peng–Robinson's coefficient edges only (Rust's cubic is unchanged).
- [ ] `phase::raoult_total_binary` (found in 2.0.1's review, missed by the survey): x·P₁ + (1 − x)·P₂ is a sum of two
  products, so it gets #20's semantics. `raoult_total_binary(0.5, 5e-324, 5e-324)` is +0 today (each product
  underflows), where #20's gives 5e-324. It needs an ADR 0002 row and edge rows.
- [ ] M5a's spots. M5a was ported at the plan's parity classes; these are listed, not fixed:
  - `simulate_michaelis_menten`'s rate Vmax·[S]/(Km + [S]) over- or underflows in Vmax·[S] or Km + [S] where the rate is
    representable (kinetics' `michaelis_menten` is #28's).
  - `simulate_second_order`'s −k·[A]·[A] is a product chain, and the consecutive / reversible closures' k₁[A] − k₂[B]
    are sums with cancellation, as in #28 and #33.
  - RK4's increment h/6·(k₁ + 2k₂ + 2k₃ + k₄) and its stage points y + ½h·k: the rate sum overflows where h/6 times it
    is representable (huge rates, tiny steps).
  - dopri45's stage combinations Σ aᵢⱼ·kⱼ, and its error terms before the × hd. The error norm's squares and the
    tol·|y| scale need nothing: they over- or underflow only far beyond the point where the step factor clamps (err below
    1.9·10⁻⁴ or above 1,845, with the accept test err ≤ 1).
  - `simulate_stiff`'s h·J and h·f products, and its Newton tolerances 10⁻¹⁰·‖y₊‖∞ and 10⁻⁶·‖y₊‖∞. The tolerances
    underflow for ‖y₊‖∞ below about 10⁻²⁹⁸, and a state that small then needs an exact zero correction to converge.
  - Euler–Maruyama's a·dt + b·dW, and the Gaussian elimination's factor·a products (hisab 1.4.0's, as transcribed).
- [ ] M5b's spots, ported the same way:
  - The least squares' Vandermonde powers x^j, and the QR's unscaled dot products and norms Σ q², over- or underflow.
    The result is NaN coefficients, or a failed 10⁻¹² test, where the fit of a rescaled x is representable.
    `fit_arrhenius`' 1/T of a subnormal T is ∞.
  - The Jacobi rank's sums app, aqq and apq, its app·aqq and 1 + τ² (τ past 1.3·10¹⁵⁴), and the singular values'
    Σ b² over- or underflow. A huge matrix counts +∞ singular values and a tiny one counts none, because the 10⁻¹⁰
    threshold is absolute.
  - `rate_from_stoichiometry`'s Σ S·v (a sum with cancellation, #33's kind) and `stoichiometric_matrix`'s repeated
    coefficient sums.
  - The gradients' K/(R·T), (K·(−ΔG))/((R·T)·(R·T)), A·e and −RT/(nFQ): (R·T)² or a product over- or underflows where
    the derivative is representable, giving the port's ComputationError or a 0.

**Non-finite potential inputs, waiting on a decision.** `potential::morse` (a non-finite r or r_e) and `coulomb` /
`coulomb_force` (a non-finite charge) return `Ok` of a NaN or ±∞, as Rust's documented API does: those inputs are not
validated. Returning `Err(KM_ERR_INVALID_INPUT)` instead, as the rest of the library does, would be an ADR 0002 change.

**The oracle defect, waiting on a decision.**

- [ ] `network::reaction_network_rank` (found in M5b). hisab 1.4.0's one-sided Jacobi tests |a_pq| ≤
  EPSILON_F64²·√(a_pp·a_qq), which is 10⁻²⁴ relative and below the rounding of a rotated pair's dot product. Most dense
  matrices therefore end in NoConvergence after 100·max(m, n) sweeps, and kimiya 1.2.0's rank is Err(ComputationError).
  That happens for Chapman's ozone cycle, Robertson's kinetics, the H₂ + Br₂ chain, and 159 of the 563 golden matrices
  (all of them finite). The port reproduces it (ADR 0003's M5b amendment). A fix (a reachable test, with Rust's answer
  wherever it converges) is a Cyrius API change for ADR 0002.

**The deferred modules, waiting on a Cyrius consumer** (port plan §2.2):

- [ ] `bridge.rs`: the ushma, tanmatra and bijli primitive-value bridges. When it is ported, ADR 0002 #5 applies: a
  `Result` instead of the `NaN` sentinels, and `ionization_to_bond_energy_estimate` removed as a physically unfounded
  heuristic.
- [ ] `integration/soorat.rs` (`soorat-compat`): the visualization data for molecular structures, reaction networks,
  spectra and phase diagrams. When it is ported, ADR 0002 #10 applies: `MolecularStructure::from_molecule` returns a
  `Result`.
- [ ] `logging`: the Rust feature's one test is deferred with it.

### 2.1.x — data re-checks against a primary source

- [ ] Standard potentials. Each of these differs from CRC or between editions:

  | Couple | kimiya | CRC |
  |---|---|---|
  | Sn⁴⁺/Sn²⁺ | 0.154 | 0.151 |
  | ClO⁻/Cl⁻ | 0.890 | 0.81 |
  | Sn²⁺/Sn | −0.136 | −0.1375 |
  | Cr₂O₇²⁻/Cr³⁺ | 1.33 | 1.36 or 1.232, by edition |
- [ ] Shomate HF(g): a known defect, Cp below 7R/2. It needs a sourced NIST/JANAF fit.
- [ ] Shomate C₂H₆ (+3.4 % at 298 K) and C₃H₈ (+8.4 %): no NIST Shomate set exists, so replace them with a sourced fit
  (NASA-7 / Burcat).
- [ ] Shomate high-T segments for NH₃ (> 1400 K), H₂S, SO₂ and others. Out-of-range use returns `Err` today.
- [ ] Antoine acetone: the range −20–77 °C is unconfirmed.
- [ ] `organic` bond energies: the 45 average bond energies (kJ/mol) carry no source citation. Tie them to a primary
  compilation, such as the CRC Handbook's or Luo's bond dissociation energies, and record the source.
- [ ] IUPAC 2011 hydrogen-bond definition: the acceptor sentence quoted in `src/organic_data.cyr` came from a secondary
  copy. Check it against the published *Pure Appl. Chem.* 83, 1637 PDF.
- [ ] Na₂CO₃: its ΔGf° is inconsistent with ΔHf° and S° by 2.55 kJ/mol in NBS/CRC. It is pinned in the tests as a known
  exception.

### 2.2.x — additive API

None of these breaks an existing call.

- [ ] `fitting::fit_beer_lambert`: also return the intercept.
- [ ] `spectroscopy::rydberg_wavelength`: add a nuclear-mass parameter for the reduced-mass correction.
- [ ] `reaction_dynamics::simulate_stiff`: return the full trajectory (today it returns the endpoints only).

### 2.x — demand-gated

- [ ] A collision-frequency helper σ·√(8k_BT/πμ)·N_A for `kinetics::collision_theory_rate`.
- [ ] Liquid/stable-phase selection for Peng–Robinson fugacity (the roots are available via `peng_robinson_z_roots`).
- [ ] The non-metal hydride rule (NaBH₄, SiH₄) in oxidation-state assignment.
- [ ] Building CSR stoichiometric matrices directly from triplets.
- [ ] An optional CODATA 2022 update for the measured constants (mₚ, mₙ, u, R∞).
- [ ] Faster name lookups (2.0.1's review): generated sorted indexes for `lookup_half_reaction` and
  `lookup_ionic_radius` (`src/electrochemistry.cyr`, `src/inorganic.cyr`: a linear scan at about 4 ns a row today), and
  a per-first-byte index for the element symbols.
