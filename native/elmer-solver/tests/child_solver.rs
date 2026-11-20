//! Child solvers (`CreateChildSolver`) as the block solver's Schur matrix
//! (`Block Matrix Schur Variable`, `BlockPrecMatrix`): the real Helmholtz
//! system of `tests/data/complex/direct.sif` split into its two components,
//! with a child `schur` of the parent's pattern holding the real part
//! without the Dirichlet rows, assembled by a test-local HelmholtzSolve.

mod support;

use std::path::Path;

use elmer_core::{BasisEval, ElemId, ElementBasis, LocalSystem, ParamMap, Params, Quadrature, SerialComm, Value};
use elmer_solver::{Simulation, Solver, SolverCtx, SolverError, SolverRegistry};
use support::{sif_model, Capture};

/// HelmholtzSolve (as `support::poisson::Helmholtz`) that also fills a child
/// solver `schur` of `CHILD_DOFS` dofs with the real part.
struct HelmholtzChild<const CHILD_DOFS: usize>;

impl<const CHILD_DOFS: usize> Solver for HelmholtzChild<CHILD_DOFS> {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        ctx.child_initialize("schur", CHILD_DOFS)?;
        let omega = ctx.simulation_params().real("angular frequency").expect("Angular Frequency");
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = BasisEval::default();
        for e in active {
            let mat = ctx.body_params(e).and_then(|b| b.material).expect("material");
            let (c, damp, rho) = (
                mat.real("sound speed").unwrap(),
                mat.real("sound damping").unwrap_or(0.0),
                mat.real("density").unwrap_or(1.0),
            );
            let k = omega / c;
            let (m, d) = (-k * k, k * damp);
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&i| ctx.mesh.coord(i)).collect();
            let mut local = LocalSystem::zeros(2 * n);
            let mut child = LocalSystem::zeros(CHILD_DOFS * n);
            for qp in support::quad::TestQuadrature.gauss_points(el.etype, None) {
                support::basis::TestBasis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|msg| SolverError::fatal("HelmholtzSolve", msg))?;
                let s = qp.weight * ev.det_j;
                for p in 0..n {
                    for q in 0..n {
                        let g = (ev.dbasis_dx[q], ev.dbasis_dx[p]);
                        let grad = g.0[0] * g.1[0] + g.0[1] * g.1[1] + g.0[2] * g.1[2];
                        let re = s * ((m * ev.basis[q] * ev.basis[p] + grad) / rho);
                        let im = s * ((d * ev.basis[q] * ev.basis[p]) / rho);
                        let w = 2 * n;
                        local.stiff[2 * p * w + 2 * q] += re;
                        local.stiff[2 * p * w + 2 * q + 1] -= im;
                        local.stiff[(2 * p + 1) * w + 2 * q] += im;
                        local.stiff[(2 * p + 1) * w + 2 * q + 1] += re;
                        let cw = CHILD_DOFS * n;
                        for l in 0..CHILD_DOFS {
                            child.stiff[(CHILD_DOFS * p + l) * cw + CHILD_DOFS * q + l] += re;
                        }
                    }
                }
            }
            ctx.default_update(&local, e)?;
            ctx.child_update("schur", &child, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.child_dirichlet("schur")?;
        let rows = ctx.child_system("schur").map(|s| s.matrix.nrows);
        assert_eq!(rows, Some(CHILD_DOFS * ctx.system().unwrap().matrix.nrows / 2));
        ctx.default_solve()?;
        Ok(())
    }
}

fn child<const D: usize>(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(HelmholtzChild::<D>))
}

fn run(extra: &[(&str, Value)], child_dofs: usize, bc: Option<(&str, Value)>) -> Result<Simulation, SolverError> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mesh = support::mesh::FileMesh::read(&dir.join("1d"), 1);
    let mut model = sif_model(&dir.join("complex/direct.sif"));
    for (k, v) in extra {
        model.solvers[0].insert(k, v.clone());
    }
    if let Some((k, v)) = bc {
        model.boundary_conditions[0].insert(k, v);
    }
    let mut reg = SolverRegistry::new();
    match child_dofs {
        1 => reg.register("HelmholtzSolve", "HelmholtzSolver", child::<1>),
        _ => reg.register("HelmholtzSolve", "HelmholtzSolver", child::<2>),
    }
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.setup(&mesh, &reg, &log, &SerialComm)?;
    sim.run(&mesh, &log, &SerialComm)?;
    Ok(sim)
}

fn keywords(schur: bool) -> Vec<(&'static str, Value)> {
    let s = |v: &str| Value::Str(v.into());
    let mut kw = vec![
        ("Linear System Solver", s("block")),
        ("Block Preconditioner", Value::Logical(true)),
        ("Outer: Linear System Iterative Method", s("gcr")),
        ("Outer: Linear System Max Iterations", Value::Integer(vec![200])),
        ("Outer: Linear System Convergence Tolerance", Value::real(1e-12)),
        ("block 11: Linear System Solver", s("direct")),
        ("block 11: Linear System Direct Method", s("umfpack")),
        ("block 22: Linear System Solver", s("direct")),
        ("block 22: Linear System Direct Method", s("umfpack")),
    ];
    if schur {
        kw.push(("Block Matrix Schur Variable", s("schur")));
    }
    kw
}

fn pressure(sim: &Simulation) -> Vec<f64> {
    sim.vars.get("pressure").unwrap().values.to_vec()
}

#[test]
fn a_child_matrix_preconditions_the_last_block() {
    // The child (the real part without the Dirichlet rows) stands in for
    // block 22 in the preconditioner: the outer iteration stops elsewhere
    // than with block 22 itself, near the direct solution.
    let s = |v: &str| Value::Str(v.into());
    let direct = [("Linear System Solver", s("direct")), ("Linear System Direct Method", s("umfpack"))];
    let exact = pressure(&run(&direct, 1, None).unwrap());
    let plain = pressure(&run(&keywords(false), 1, None).unwrap());
    let with = pressure(&run(&keywords(true), 1, None).unwrap());
    let scale = exact.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(scale > 0.0);
    assert!(plain.iter().zip(&with).any(|(a, b)| a != b), "the Schur matrix was not used");
    for (a, b) in exact.iter().zip(&with) {
        assert!((a - b).abs() <= 2e-2 * scale, "{a} vs {b}");
    }
}

#[test]
fn the_schur_variable_finds_the_child() {
    // A child of two dofs is the size of both blocks, not of block 22: the
    // block solver got it, and refuses it.
    let e = run(&keywords(true), 2, None).err().unwrap();
    assert!(matches!(&e, SolverError::Fatal { caller, msg } if caller == "BlockPrecMatrix" && msg.contains("does not fit")), "{e}");
}

#[test]
fn dirichlet_conditions_of_a_child_are_refused() {
    let e = run(&keywords(true), 1, Some(("schur", Value::real(0.0)))).err().unwrap();
    assert!(matches!(e, SolverError::Unsupported(_)), "{e}");
}

