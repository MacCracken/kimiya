//! `constants.tsv` — every public scalar `f64` const of kimiya 1.2.0, class A.
//!
//! One row per constant, `const A s:<module>::<NAME> f:<bits>`. The list is shared with the
//! `src/constants.cyr` generator (`t_constants::all`), so the two can't drift apart. The rows
//! pin the physical constants in `src/constants.cyr` and the module data constants in their
//! module's generated data file (M2a: phase, gas; M4b: inorganic). None is pending since M4b.

use crate::t_constants;
use crate::tsv::{Class, Module, f, s};

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    for c in t_constants::all() {
        m.row("const", Class::A, &[s(c.path())], &[f(c.value)]);
    }
    m
}
