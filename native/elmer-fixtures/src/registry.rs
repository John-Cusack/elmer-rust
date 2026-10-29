//! Case-qualified fixture selection.
//!
//! Test plugins are named by file, not by purpose: `Source` is the 1dtests
//! right-hand side, a manufactured source in `BDM2D` and `MixedPoisson*`, and
//! more, each with different arguments. A fixture is therefore keyed by
//! `(case, library, entry)` where `case` is the `fem/tests/<dir>` name, and
//! the caller picks the case it is running. Nothing is registered process-wide
//! on its own.
//!
//! Each call builds new fixture state, so two simulations (or two cases in
//! one test process) never share `SAVE`d Fortran variables.

use std::sync::Arc;

use elmer_expr::{Udf, UdfRegistry};

use crate::udf::{init_dist_fn, my_source_fn, OneDSource, OneDTolerances};
use crate::ModelView;

type Build = fn(Arc<dyn ModelView>) -> Udf;

/// `(case, library, entry, builder)`. Library and entry are the strings the
/// case's `Procedure` keyword uses.
const FIXTURES: &[(&str, &str, &str, Build)] = &[
    ("ConstantParamFunc", "MySource", "MySource", |_| Udf::Node(my_source_fn())),
    ("AdvDiffFCT", "InitDist", "InitDist", |m| Udf::Node(init_dist_fn(m))),
    ("1dtests", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::ONE_D_TESTS).into_fn())
    }),
    ("1dtests_MeshLevels", "Source", "Source", |m| {
        Udf::Node(OneDSource::new(m, OneDTolerances::MESH_LEVELS).into_fn())
    }),
];

/// The `(library, entry)` pairs ported for `case`, empty when none are.
pub fn procedures(case: &str) -> Vec<(&'static str, &'static str)> {
    FIXTURES.iter().filter(|f| f.0 == case).map(|f| (f.1, f.2)).collect()
}

/// A new instance of the fixture `case` uses for `Procedure "library" "entry"`,
/// bound to `model`. `None` when that plugin is not ported.
pub fn lookup(case: &str, library: &str, entry: &str, model: Arc<dyn ModelView>) -> Option<Udf> {
    FIXTURES
        .iter()
        .find(|f| f.0 == case && f.1 == library && f.2.eq_ignore_ascii_case(entry))
        .map(|f| f.3(model))
}

/// Registers every fixture of `case` into `registry` under the names the
/// case's SIF uses, as new instances. Returns how many were registered; zero
/// means the case has no ported plugins.
pub fn install(case: &str, model: Arc<dyn ModelView>, registry: &UdfRegistry) -> usize {
    let mut n = 0;
    for (c, lib, entry, build) in FIXTURES {
        if *c == case {
            registry.register_fortran(lib, entry, build(Arc::clone(&model)));
            n += 1;
        }
    }
    n
}
