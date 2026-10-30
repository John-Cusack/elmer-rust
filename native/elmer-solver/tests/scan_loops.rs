//! `Scanning Loops` with `Pre Solvers`, `Nonlinear Pre Solvers`, `Nonlinear
//! Post Solvers` and `Post Solvers`, checked against a pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run of `tests/data/scanloops/case.sif`
//! (HeatSolve on the LinearSystemNamespace angle mesh, `tests/data/namespace/`,
//! with `Heat Source = Equals "scan"` and two scanning loops; the slot solver
//! is `fem/tests/SlaveSolverSlotsScan`'s CountSolver2): the order of the
//! count visits and the `ComputeChange` lines. HeatSolve and CountSolver are
//! test-local stand-ins.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::SerialComm;
use elmer_solver::{Simulation, Solver, SolverCtx, SolverError, SolverRegistry, TextMatc};
use support::poisson::{heat, heat_nl};
use support::{sif_model, split_change, Capture};

/// CountSolver2: the global variable and its norm hold the visit count.
struct Count(usize);

impl Solver for Count {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        self.0 += 1;
        ctx.info("CountSolver", &format!("Solver visited times: {}", self.0), 4);
        let name = ctx.variable_name().expect("count variable").to_string();
        let v = ctx.sim.vars.get_mut(&name).expect("count variable");
        v.values.fill(self.0 as f64);
        v.conv.norm = self.0 as f64;
        Ok(())
    }
}

fn count(_: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Count(0)))
}

fn run(heat_solver: elmer_solver::SolverFactory) -> (Result<(), SolverError>, Simulation, Capture) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&data.join("namespace"), 2);
    let model = sif_model(&data.join("scanloops/case.sif"));
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_solver);
    reg.register("CountSolver2", "CountSolver", count);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    let r = sim.run(&mesh, &log, &SerialComm).map(|_| ());
    (r, sim, log)
}

#[test]
fn scanning_loops_run_the_slot_solvers_in_pristine_order() {
    let (r, sim, log) = run(heat_nl);
    r.unwrap();
    assert_eq!(sim.vars.scalar("scan"), Some(2.0));
    assert_eq!(sim.vars.get("count").unwrap().values, vec![12.0]);

    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let stdout = std::fs::read_to_string(data.join("scanloops/stdout.txt")).unwrap();
    let want: Vec<&str> = stdout.lines().collect();
    let got = log.lines(&["ComputeChange", "CountSolver"]);
    assert_eq!(got.len(), want.len(), "{got:#?}");
    for (g, w) in got.iter().zip(&want) {
        if w.starts_with("CountSolver") {
            assert_eq!(g, w);
            continue;
        }
        let (gh, gn, gr, gt) = split_change(g);
        let (wh, wn, wr, wt) = split_change(w);
        assert_eq!((gh, gn, gt), (wh, wn, wt));
        assert!((gr - wr).abs() <= 1e-9 * wr, "{g}\n{w}");
    }
}

#[test]
fn pre_solvers_need_a_procedure_that_calls_default_start() {
    let (r, _, _) = run(heat);
    match r {
        Err(SolverError::Unsupported(m)) => assert!(m.contains("pre solvers") && m.contains("DefaultStart"), "{m}"),
        other => panic!("expected the pre solvers refused, got {other:?}"),
    }
}
