//! `Apply Integral BCs` on a two-component field (`AddVectorConstraint`,
//! ProjectorUtils.F90:1469): one row per component that has an `Integral BC
//! <component>` value, over that component's dofs. The test-local solver
//! assembles two uncoupled Laplace equations, so with the
//! `IntegralConstraintPoisson` boundaries (`tests/data/intbc/poisson/`)
//! component k equals the scalar HeatSolve stand-in's solution of
//! `integral_bcs.rs` (checked there against pristine) scaled by its own
//! integral value, and a component without one stays zero.

mod support;

use std::path::Path;

use elmer_core::{BasisEval, ElemId, ElementBasis, LocalSystem, ParamMap, Quadrature, SerialComm, Value};
use elmer_solver::{FeBackend, Model, Simulation, Solver, SolverCtx, SolverError, SolverRegistry};
use support::poisson::heat;
use support::{sif_model, Capture};

/// `-div(grad u_k) = 0` for each of the two components; with `own_bcs`
/// the Dirichlet rows are gathered and applied by the solver itself
/// (`apply_dirichlet_values`, as FlowSolve does).
struct VecLaplace {
    own_bcs: bool,
}

impl Solver for VecLaplace {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut ev = BasisEval::default();
        for e in active {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let mut local = LocalSystem::zeros(2 * n);
            for qp in support::quad::TestQuadrature.gauss_points(el.etype, None) {
                support::basis::TestBasis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|m| SolverError::fatal("VecLaplace", m))?;
                let w = qp.weight * ev.det_j;
                for p in 0..n {
                    for q in 0..n {
                        let (a, b) = (ev.dbasis_dx[p], ev.dbasis_dx[q]);
                        let s = w * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
                        for k in 0..2 {
                            local.stiff[(2 * p + k) * 2 * n + 2 * q + k] += s;
                        }
                    }
                }
            }
            ctx.default_update(&local, e)?;
        }
        if self.own_bcs {
            ctx.default_finish_boundary_assembly()?;
            let fixed = ctx.dirichlet_values()?;
            ctx.apply_dirichlet_values(&fixed)?;
        } else {
            ctx.default_dirichlet()?;
        }
        ctx.default_solve()?;
        Ok(())
    }
}

fn vec_laplace(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(VecLaplace { own_bcs: false }))
}

fn vec_laplace_own_bcs(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(VecLaplace { own_bcs: true }))
}

fn solve(model: Model, entry: (&str, elmer_solver::SolverFactory)) -> Simulation {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/intbc/poisson");
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut reg = SolverRegistry::new();
    reg.register(entry.0, "HeatSolver", entry.1);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    sim
}

/// The scalar case with the field `T` of two components and `Integral BC
/// T k` set from `integrals`.
fn vector_model(integrals: &[(usize, f64)]) -> Model {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/intbc/poisson");
    let mut m = sif_model(&dir.join("case.sif"));
    let s = &mut m.solvers[0];
    s.insert("Procedure", Value::Str("VecLaplace HeatSolver".into()));
    s.insert("Variable", Value::Str("-dofs 2 T".into()));
    s.0.remove("reference norm");
    let fixed = &mut m.boundary_conditions[1];
    fixed.0.remove("temperature");
    fixed.insert("T 1", Value::real(0.0));
    fixed.insert("T 2", Value::real(0.0));
    let constraint = &mut m.boundary_conditions[2];
    constraint.0.remove("integral bc temperature");
    for &(k, v) in integrals {
        constraint.insert(&format!("Integral BC T {k}"), Value::real(v));
    }
    m
}

fn scalar_solution() -> Vec<f64> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/intbc/poisson");
    let sim = solve(sif_model(&dir.join("case.sif")), ("HeatSolve", heat));
    let t = sim.vars.get("temperature").unwrap();
    let perm = t.perm.as_deref().unwrap();
    perm.iter().map(|&k| if k == u32::MAX { f64::NAN } else { t.values[k as usize] }).collect()
}

fn components(sim: &Simulation) -> [Vec<f64>; 2] {
    let t = sim.vars.get("t").unwrap();
    let perm = t.perm.as_deref().unwrap();
    let at = |k: usize| perm.iter().map(|&s| if s == u32::MAX { f64::NAN } else { t.values[2 * s as usize + k] }).collect();
    [at(0), at(1)]
}

fn assert_scaled(got: &[f64], scalar: &[f64], factor: f64, what: &str) {
    for (node, (g, s)) in got.iter().zip(scalar).enumerate() {
        if s.is_nan() {
            continue;
        }
        assert!((g - factor * s).abs() <= 1e-10 * (1.0 + s.abs()), "{what} node {}: {g} vs {}", node + 1, factor * s);
    }
}

#[test]
fn each_component_gets_its_own_integral_row() {
    let scalar = scalar_solution();
    let sim = solve(vector_model(&[(1, 1.0), (2, 2.0)]), ("VecLaplace", vec_laplace));
    let [t1, t2] = components(&sim);
    assert_scaled(&t1, &scalar, 1.0, "T 1");
    assert_scaled(&t2, &scalar, 2.0, "T 2");
}

#[test]
fn a_component_without_a_value_is_not_constrained() {
    let scalar = scalar_solution();
    let sim = solve(vector_model(&[(1, 1.0)]), ("VecLaplace", vec_laplace));
    let [t1, t2] = components(&sim);
    assert_scaled(&t1, &scalar, 1.0, "T 1");
    assert_scaled(&t2, &scalar, 0.0, "T 2");
}

/// A solver applying its own Dirichlet rows records them as constrained,
/// so the integral constraint moves their columns to the right-hand side
/// (`Enforce Exact Dirichlet BCs`): the Dirichlet values hold exactly and
/// the solution is the one `default_dirichlet` gives.
#[test]
fn own_dirichlet_rows_are_constrained_dofs() {
    // Boundary 1 fixed too, so the integral boundary's end nodes are
    // Dirichlet dofs inside the constraint row.
    let model = || {
        let mut m = vector_model(&[(1, 1.0), (2, 2.0)]);
        m.boundary_conditions[0].insert("T 1", Value::real(0.0));
        m.boundary_conditions[0].insert("T 2", Value::real(0.0));
        m
    };
    let reference = solve(model(), ("VecLaplace", vec_laplace));
    let own = solve(model(), ("VecLaplace", vec_laplace_own_bcs));
    let (r, o) = (components(&reference), components(&own));
    for k in 0..2 {
        for (node, (x, y)) in r[k].iter().zip(&o[k]).enumerate() {
            if x.is_nan() {
                continue;
            }
            assert!((x - y).abs() <= 1e-12 * (1.0 + x.abs()), "T {} node {}: {x} vs {y}", k + 1, node + 1);
        }
    }
}
