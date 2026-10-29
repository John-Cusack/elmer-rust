//! Mortar cases end to end: Galerkin mortar projectors, their constraint
//! matrix and the saddle-point solve, compared with the pristine run's
//! temperature and Lagrange multipliers. MortarPoisson2D has one block of
//! constraint rows per projector; `mortar_sum2d` and `mortar_priority2d`
//! (MortarPoisson2Dsum/priority without `Plane Projector`) sum four
//! projectors (`Mortar BCs Additive`), the second with `Projector Priority`.
//!
//! `elmer-solver` has no constrained solve yet, so the system comes from a
//! `Simulation` run with `Linear System Solver Disabled` (assembly and
//! Dirichlet conditions only), and the collection system is solved here as
//! `SolveLinearSystem` does: diagonal scaling, the section's iterative
//! solver and preconditioner, back-scaling.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_core::{
    CsrMatrix, ElemId, LinearOperator, LocalSystem, MeshView, NullLogger, ParamMap, Params,
    SerialComm, Value,
};
use elmer_mesh::adapters::MeshViewAdapter;
use elmer_mortar::{
    collection_system, generate_constraint_matrix, generate_mortar_projectors, ConstraintTarget,
    ProjectorEnv, Restriction,
};
use elmer_run::{coordinate_system, install_fe, solver_registry, SifMatc};
use elmer_solver::backend::LinearSolverChoice;
use elmer_solver::{AssemblyBackend, ElmerAssembly, Perm, Simulation, SolverError};

/// The default assembly, recording the Dirichlet rows it is given.
struct Recording {
    fixed: Arc<Mutex<Vec<(u32, f64)>>>,
}

impl AssemblyBackend for Recording {
    fn build_perm(&self, mesh: &dyn MeshView, active: &[ElemId], opt: bool) -> Perm {
        ElmerAssembly.build_perm(mesh, active, opt)
    }
    fn build_perm_for(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        params: &ParamMap,
    ) -> Result<Perm, SolverError> {
        ElmerAssembly.build_perm_for(mesh, active, params)
    }
    fn create_matrix(
        &self,
        mesh: &dyn MeshView,
        active: &[ElemId],
        perm: &Perm,
        dofs: usize,
    ) -> CsrMatrix<f64> {
        ElmerAssembly.create_matrix(mesh, active, perm, dofs)
    }
    fn add_local(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        slots: &[u32],
        local: &LocalSystem,
    ) {
        ElmerAssembly.add_local(a, rhs, slots, local);
    }
    fn apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        p: &ParamMap,
    ) {
        self.try_apply_dirichlet(a, rhs, fixed, p).unwrap();
    }
    fn try_apply_dirichlet(
        &self,
        a: &mut CsrMatrix<f64>,
        rhs: &mut [f64],
        fixed: &[(u32, f64)],
        params: &ParamMap,
    ) -> Result<(), SolverError> {
        *self.fixed.lock().unwrap() = fixed.to_vec();
        ElmerAssembly.try_apply_dirichlet(a, rhs, fixed, params)
    }
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        std::fs::copy(e.path(), to.join(e.file_name())).unwrap();
    }
}

fn fixture(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(case)
}

/// The pristine `Output File`: the first field's values in node order and
/// the Lagrange multipliers.
fn read_result(path: &Path) -> (Vec<f64>, Vec<f64>) {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let perm_at = lines.iter().position(|l| l.starts_with("Perm:")).unwrap();
    let n: usize = lines[perm_at]
        .split_whitespace()
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    let mut slot_of = vec![0usize; n];
    for l in &lines[perm_at + 1..perm_at + 1 + n] {
        let mut f = l.split_whitespace().map(|x| x.parse::<usize>().unwrap());
        let (node, slot) = (f.next().unwrap(), f.next().unwrap());
        slot_of[node - 1] = slot - 1;
    }
    let vals: Vec<f64> = lines[perm_at + 1 + n..perm_at + 1 + 2 * n]
        .iter()
        .map(|l| l.trim().parse().unwrap())
        .collect();
    let by_node = slot_of.iter().map(|&s| vals[s]).collect();
    // The field section, not the header's list of degrees of freedom.
    let m_at = lines
        .iter()
        .rposition(|l| l.starts_with("lagrange multiplier"))
        .unwrap();
    let mult = lines[m_at + 2..]
        .iter()
        .map(|l| l.trim().parse().unwrap())
        .collect();
    (by_node, mult)
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
}

fn max_abs(a: &[f64]) -> f64 {
    a.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// Solves `fixtures/<case>` and compares it with the pristine run, whose
/// `Solver 1` norm is `pristine` (8 digits): the temperature to `tol_t`
/// and the multipliers to `tol_l`, relative to their largest value.
fn check(case: &str, pristine: f64, tol_t: f64, tol_l: f64) {
    let src = fixture(case);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("elmer-mortar-e2e")
        .join(case);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(src.join("case.sif"), dir.join("case.sif")).unwrap();
    copy_dir(&src.join("mesh"), &dir.join("mesh"));

    let sif = dir.join("case.sif");
    let path = sif.to_str().unwrap();
    let mut parsed = elmer_sif::parse_file(path).unwrap();
    elmer_sif::complete(&mut parsed.sif, path).unwrap();
    let mut model = elmer_run::sif::to_model(&parsed.sif.runs[0]).unwrap();
    // Solvers 2 and 3 never run (`Exec Solver = never`).
    model.solvers.truncate(1);
    // Assemble only; the constrained solve is done below.
    model.solvers[0].insert("Linear System Solver Disabled", Value::Logical(true));
    let solver = model.solvers[0].clone();
    let mesh = elmer_mesh::load(dir.join("mesh").to_str().unwrap()).unwrap();
    let view = MeshViewAdapter::new(&mesh);
    let coords = coordinate_system(&model, &view).unwrap();
    install_fe(coords.dim).unwrap();

    let fixed = Arc::new(Mutex::new(Vec::new()));
    let mut sim = Simulation::new(model.clone());
    sim.set_matc(Rc::new(SifMatc));
    sim.backends.assembly = Some(Box::new(Recording {
        fixed: Arc::clone(&fixed),
    }));
    sim.setup(&view, &solver_registry(), &NullLogger, &SerialComm)
        .unwrap();
    sim.run(&view, &NullLogger, &SerialComm).unwrap();

    let slot = sim.slots.iter().find(|s| s.id == 1).unwrap();
    let var = sim.vars.by_id(slot.variable.unwrap());
    let perm = var.perm.as_deref().unwrap().to_vec();
    let sys = slot.system.as_ref().unwrap();

    // GenerateProjectors, GenerateConstraintMatrix.
    let bc_of = |tag: i32| {
        model
            .boundary_conditions
            .iter()
            .position(|bc| matches!(bc.get("target boundaries"), Some(Value::Integer(v)) if v.contains(&tag)))
            .map(|i| i + 1)
    };
    let bcs: Vec<&dyn Params> = model
        .boundary_conditions
        .iter()
        .map(|b| b as &dyn Params)
        .collect();
    let env = ProjectorEnv {
        mesh: &view,
        bcs: &bcs,
        bc_of: &bc_of,
        simulation: &model.simulation,
        solver: Some(&solver),
        coord_dim: coords.dim,
        log: &NullLogger,
    };
    let (projectors, _) = generate_mortar_projectors(&env, &[], false, &dir, 0)
        .unwrap()
        .unwrap();
    let target = ConstraintTarget {
        solver: &solver,
        name: "temperature",
        dofs: 1,
        perm: &perm,
        n_rows: sys.matrix.nrows,
    };
    let c = generate_constraint_matrix(&bcs, &projectors, &target, &NullLogger)
        .unwrap()
        .unwrap();

    // SolveWithLinearRestriction.
    let fixed = fixed.lock().unwrap().clone();
    let x0 = vec![0.0; sys.matrix.nrows];
    let mut coll = collection_system(&Restriction {
        params: &solver,
        a: &sys.matrix,
        b: &sys.rhs,
        x: &x0,
        c: &c,
        constrained: &fixed,
        multipliers: None,
        perm: &perm,
        dofs: 1,
        n_nodes: view.n_nodes(),
    })
    .unwrap();
    // `Linear System Scaling` (default on) with its method, as
    // `SolveLinearSystem` applies it.
    let bnorm = coll.rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
    let scale = solver.logical("linear system skip scaling") != Some(true)
        && solver.logical("linear system scaling").unwrap_or(true);
    let scaling = scale.then(|| {
        let method = solver
            .string("linear system scaling method")
            .unwrap_or("diagonal");
        elmer_sparse::scale_system(
            &mut coll.matrix,
            Some(&mut coll.rhs[..]),
            Some(&mut coll.x[..]),
            &elmer_sparse::ScaleOptions {
                method: elmer_sparse::ScalingMethod::from_keyword(method).unwrap(),
                apply: true,
                rhs_scaling: bnorm != 0.0,
            },
        )
    });
    let solver_params = &solver;
    match sim
        .backends
        .linear
        .select_for(solver_params, &coll.matrix, &NullLogger)
        .unwrap()
    {
        LinearSolverChoice::Iterative {
            mut solver,
            precond,
        } => {
            let m = precond.unwrap()(solver_params, &coll.matrix).unwrap();
            solver
                .solve(&coll.matrix, &*m, &coll.rhs, &mut coll.x, &SerialComm)
                .unwrap();
        }
        LinearSolverChoice::Direct(_) => panic!("the case uses an iterative solver"),
    }
    if let Some(scaling) = scaling {
        scaling.unscale(
            &mut coll.matrix,
            Some(&mut coll.rhs[..]),
            Some(&mut coll.x[..]),
        );
    }

    // The constrained system is solved: residual and constraint.
    let mut ax = vec![0.0; coll.matrix.nrows];
    coll.matrix.apply(&coll.x, &mut ax);
    let res = max_diff(&ax, &coll.rhs);
    println!("{case}: max residual {res:e}");

    let (want_t, want_l) = read_result(&src.join("expected/case.result"));
    let got_t: Vec<f64> = perm.iter().map(|&s| coll.solution()[s as usize]).collect();
    let dt = max_diff(&got_t, &want_t);
    let dl = max_diff(coll.multipliers(), &want_l);
    let n = got_t.len() as f64;
    let norm = (got_t.iter().map(|x| x * x).sum::<f64>() / n).sqrt();
    println!(
        "{case}: norm {norm:.10e}; max |dT| {dt:e} (max |T| {:e}); max |dl| {dl:e} (max |l| {:e})",
        max_abs(&want_t),
        max_abs(&want_l)
    );
    assert!(
        (norm - pristine).abs() / pristine <= 1e-8,
        "{case}: norm {norm}, pristine {pristine}"
    );
    assert!(
        dt <= tol_t * max_abs(&want_t),
        "temperature differs by {dt:e}"
    );
    assert!(
        dl <= tol_l * max_abs(&want_l),
        "multipliers differ by {dl:e}"
    );
}

/// The same system through the same solver path: round-off.
#[test]
fn mortar_poisson2d_constrained_solve() {
    check("mortar_poisson2d", 4.58536086e-1, 1e-12, 1e-12);
}

/// Without `Linear System Scaling` BiCGStabl stops at its 1e-8 tolerance on
/// a different path, so the solutions agree to that tolerance (about 1e-9
/// in T and 1e-6 in the multipliers). With the tolerance at 1e-14 in both
/// runs they agree to 2e-14 in T and 5e-12 in the multipliers.
#[test]
fn mortar_sum2d_constrained_solve() {
    check("mortar_sum2d", 8.74343079e-1, 1e-8, 1e-5);
}

#[test]
fn mortar_priority2d_constrained_solve() {
    check("mortar_priority2d", 8.74334619e-1, 1e-8, 1e-5);
}
