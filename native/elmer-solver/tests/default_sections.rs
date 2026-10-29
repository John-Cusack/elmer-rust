//! A model without Body and Equation sections (`fem/tests/L2ndDerivatives`,
//! `tests/data/l2nd/`): `Model::add_default_sections` adds what `LoadModel`
//! does, so its one solver is active on all 20 elements, as in the pristine
//! ElmerSolver (26.2-devel, Rev 9f6af2f85) run; the solver has no
//! `Variable`, so `CompareToReferenceSolution` checks the null variable's
//! norm, 0, against `Reference Norm = 0` as the pristine run does.

mod support;

use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use elmer_core::SerialComm;
use elmer_solver::{Simulation, Solver, SolverCtx, SolverError, SolverRegistry};
use support::{sif_model, Capture};

/// The active elements of the last visit.
static ACTIVE: AtomicUsize = AtomicUsize::new(usize::MAX);

struct Count;

impl Solver for Count {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ACTIVE.store(ctx.active_elements().count(), Ordering::SeqCst);
        Ok(())
    }
}

fn count(_: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Count))
}

fn run(defaults: bool) -> (usize, Simulation, Capture) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/l2nd");
    let mesh = support::mesh::FileMesh::read(&dir, 3);
    let mut model = sif_model(&dir.join("case.sif"));
    if defaults {
        model.add_default_sections();
    }
    ACTIVE.store(usize::MAX, Ordering::SeqCst);
    let mut reg = SolverRegistry::new();
    reg.register("L2ndDerivativesTest", "TestSolver", count);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    (ACTIVE.load(Ordering::SeqCst), sim, log)
}

#[test]
fn defaults_activate_every_element_and_the_null_variable_passes() {
    let (active, sim, log) = run(true);
    assert_eq!(active, 20);
    let stdout = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/l2nd/stdout.txt")).unwrap();
    let got = log.lines(&["WARNING:: CompareTo", "CompareToReferenceSolution: Solver", "CompareToReferenceSolution: Absolute"]);
    assert_eq!(got, stdout.lines().collect::<Vec<_>>());
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1));
    let m = &sim.model;
    assert_eq!((m.equations.len(), m.bodies.len(), m.materials.len()), (1, 1, 1));
}

#[test]
fn without_defaults_nothing_is_active() {
    let (active, _, _) = run(false);
    assert_eq!(active, 0);
}
