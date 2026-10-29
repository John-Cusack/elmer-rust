//! `FluxSolver` end to end through `elmer-solver`'s `Simulation`, on the
//! 300-element quad mesh of the pristine quick tests `fluxsolver` and
//! `fluxsolver2` (`angle.grd`), against the pristine ElmerSolver's output.
//!
//! The Temperature field fed to the solver is the pristine run's own
//! (`tests/data/*.txt`, the VTU point data at full precision), so the
//! comparison isolates the flux recovery from HeatSolve. Each data file also
//! holds the fluxes the pristine FluxSolver wrote; the `Result Norm` of its
//! last projected component is the `norm` constant in each test.

mod support;

use std::path::PathBuf;

use elmer_core::{FieldKind, NullLogger, ParamMap, Params, SerialComm, Value};
use elmer_krylov::{Huti, HutiConfig};
use elmer_physics_flux::{install_fe, register};
use elmer_solver::{Model, Simulation, SolverError, SolverRegistry, Variable};
use elmer_sparse::Ilu;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;
use support::{TestBasis, TestQuadrature};

fn data(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(rel)
}

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn ints(v: &[i32]) -> Value {
    Value::Integer(v.to_vec())
}

fn st(s: &str) -> Value {
    Value::Str(s.into())
}

fn lg(b: bool) -> Value {
    Value::Logical(b)
}

/// One pristine data file: column name -> value per node.
struct Reference {
    cols: Vec<(String, usize)>,
    rows: Vec<Vec<f64>>,
}

impl Reference {
    fn read(file: &str) -> Self {
        let text = std::fs::read_to_string(data(file)).unwrap();
        let mut cols = Vec::new();
        let mut rows = Vec::new();
        for line in text.lines() {
            if let Some(c) = line.strip_prefix("# columns: ") {
                cols = c.split(" | ").map(|s| (s.trim().to_string(), 0)).collect();
            } else if !line.starts_with('#') && !line.trim().is_empty() {
                rows.push(
                    line.split_whitespace()
                        .map(|x| x.parse().unwrap())
                        .collect(),
                );
            }
        }
        for (i, c) in cols.iter_mut().enumerate() {
            c.1 = i;
        }
        Self { cols, rows }
    }

    /// `name` is e.g. `temperature flux[1]` or `temperature`.
    fn col(&self, name: &str) -> Vec<f64> {
        let key = if name.ends_with(']') {
            name.to_string()
        } else {
            format!("{name}[0]")
        };
        let i = self
            .cols
            .iter()
            .find(|(n, _)| *n == key)
            .unwrap_or_else(|| panic!("no column {key} in {:?}", self.cols))
            .1;
        self.rows.iter().map(|r| r[i]).collect()
    }
}

/// The pristine run and this port agree to rounding (about 1e-15 here); the
/// margin is for other lanes' linear algebra changing its summation order.
const TOL: f64 = 1e-12;

fn install() {
    // Idempotent across tests in this binary.
    let _ = install_fe(|| (Box::new(TestBasis), Box::new(TestQuadrature)));
}

fn mesh() -> SimMesh {
    SimMesh::new(read_elmer_mesh(&data("angle"), 2))
}

/// The pristine quick tests' model with only the FluxSolver (`Solver 2` there)
/// and the given `Solver` keywords.
fn model(solver: &[(&str, Value)], material: &[(&str, Value)], axisymmetric: bool) -> Model {
    let mut s = vec![
        ("Equation", st("ComputeFlux")),
        ("Procedure", st("\"FluxSolver\" \"FluxSolver\"")),
        ("Target Variable", st("Temperature")),
    ];
    s.extend(solver.iter().cloned());
    let mut sim = vec![
        ("Simulation Type", st("steady state")),
        ("Steady State Max Iterations", ints(&[1])),
        ("Output Intervals", ints(&[0])),
    ];
    if axisymmetric {
        sim.push(("Coordinate System", st("Axi Symmetric")));
    }
    let mut mat = vec![("Density", Value::real(1.0))];
    mat.extend(material.iter().cloned());
    Model {
        simulation: sec(&sim),
        solvers: vec![sec(&s)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Target Bodies", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
        ])],
        materials: vec![sec(&mat)],
        ..Model::default()
    }
}

fn simulation(model: Model, temperature: Vec<f64>) -> Simulation {
    install();
    let mut sim = Simulation::new(model);
    sim.vars
        .add(Variable::new(
            "Temperature",
            FieldKind::Nodal,
            1,
            None,
            temperature,
        ))
        .unwrap();
    sim.backends.linear.register_iterative("cg", |p| {
        let cfg = HutiConfig::from_params(p)?;
        Ok(Box::new(Huti::new(cfg)))
    });
    sim.backends
        .linear
        .register_precond("ilu0", |_, a| Ok(Box::new(Ilu::ilu_n(a, 0)?)));
    sim
}

fn run_sim(sim: &mut Simulation, mesh: &SimMesh) -> Result<(), SolverError> {
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    sim.setup(mesh, &reg, &NullLogger, &SerialComm)?;
    sim.run(mesh, &NullLogger, &SerialComm)?;
    Ok(())
}

/// A variable's component `c` by node.
fn by_node(sim: &Simulation, name: &str, nodes: usize) -> Vec<f64> {
    let c = sim
        .vars
        .strided(name)
        .unwrap_or_else(|| panic!("no variable {name}"));
    let perm = c
        .perm
        .expect("exported variables share the solver permutation");
    (0..nodes).map(|n| c.get(perm[n] as usize)).collect()
}

fn assert_close(what: &str, got: &[f64], want: &[f64], tol: f64) {
    assert_eq!(got.len(), want.len(), "{what}");
    let worst = got
        .iter()
        .zip(want)
        .fold(0.0, |m: f64, (g, w)| m.max((g - w).abs()));
    assert!(worst <= tol, "{what}: max |diff| {worst:e} > {tol:e}");
}

fn norm_of(sim: &Simulation) -> f64 {
    let id = sim.slot(1).unwrap().variable.unwrap();
    sim.vars.by_id(id).conv.norm
}

fn check_norm(sim: &Simulation, want: f64) {
    let got = norm_of(sim);
    assert!(
        ((got - want) / want).abs() < TOL,
        "Result Norm {got:e} vs pristine {want:e}"
    );
}

/// fluxsolver (flux.sif): `Calculate Flux`, `Flux Coefficient = Heat
/// Conductivity` of 3, as written; none of Variable, Linear System Solver, ...
/// is given, so the defaults of `FluxSolver_Init` carry the solve. The sif's
/// own `Solver 2 :: Reference Norm = 0.198814065076` is checked through the
/// framework's CompareToReferenceSolution.
#[test]
fn fluxsolver_quick_test_flux() {
    let r = Reference::read("flux.txt");
    let mesh = mesh();
    let m = model(
        &[
            ("Calculate Flux", lg(true)),
            ("Flux Coefficient", st("Heat Conductivity")),
            ("Reference Norm", Value::real(0.198814065076)),
        ],
        &[("Heat Conductivity", Value::real(3.0))],
        false,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();

    assert_eq!(sim.reference.test_count, 1);
    assert_eq!(sim.reference.pass_count, 1, "CompareToReferenceSolution");
    check_norm(&sim, 0.19881406507522117);
    for c in 1..=2 {
        let got = by_node(&sim, &format!("temperature flux {c}"), 341);
        assert_close(
            &format!("flux {c}"),
            &got,
            &r.col(&format!("temperature flux[{}]", c - 1)),
            TOL,
        );
    }
    // The solver variable is hidden and the third component is not exported.
    assert!(!sim.vars.get("computeflux_temp").unwrap().output);
    assert!(sim.vars.resolve("temperature flux 3").is_none());
    assert!(sim.vars.resolve("temperature grad").is_none());
}

/// fluxsolver2 (flux2.sif): `Calculate Grad` only, no coefficient needed.
/// Reference norm 0.66271355E-01 from the sif.
#[test]
fn fluxsolver2_quick_test_grad() {
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let m = model(
        &[
            ("Calculate Grad", lg(true)),
            ("Reference Norm", Value::real(0.66271355e-1)),
        ],
        &[("Heat Conductivity", Value::real(3.0))],
        false,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();

    assert_eq!(sim.reference.pass_count, 1);
    check_norm(&sim, 0.06627135502507372);
    for c in 1..=2 {
        let got = by_node(&sim, &format!("temperature grad {c}"), 341);
        assert_close(
            &format!("grad {c}"),
            &got,
            &r.col(&format!("temperature grad[{}]", c - 1)),
            TOL,
        );
    }
    assert!(sim.vars.resolve("temperature flux").is_none());
}

/// All four requests with `Flux Abs`, `Grad Abs` and `Enforce Positive
/// Magnitude`: the pristine run's flux, grad, the two projected magnitudes and
/// the two pointwise absolute values (which differ from the magnitudes).
#[test]
fn all_fields_match_the_pristine_run() {
    let r = Reference::read("all.txt");
    let mesh = mesh();
    let m = model(
        &[
            ("Calculate Flux", lg(true)),
            ("Calculate Flux Abs", lg(true)),
            ("Calculate Grad", lg(true)),
            ("Calculate Grad Abs", lg(true)),
            ("Calculate Flux Magnitude", lg(true)),
            ("Calculate Grad Magnitude", lg(true)),
            ("Enforce Positive Magnitude", lg(true)),
            ("Flux Coefficient", st("Heat Conductivity")),
        ],
        &[("Heat Conductivity", Value::real(3.0))],
        false,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();

    check_norm(&sim, 0.09323802387115132);
    for c in 1..=2 {
        for (var, col) in [("flux", "flux"), ("grad", "grad")] {
            let got = by_node(&sim, &format!("temperature {var} {c}"), 341);
            let want = r.col(&format!("temperature {col}[{}]", c - 1));
            assert_close(&format!("{var} {c}"), &got, &want, TOL);
        }
    }
    for name in ["flux_abs", "flux_mag", "grad_abs", "grad_mag"] {
        let got = by_node(&sim, &format!("temperature {name}"), 341);
        assert_close(name, &got, &r.col(&format!("temperature {name}")), TOL);
    }
    let abs = by_node(&sim, "temperature flux_abs", 341);
    let mag = by_node(&sim, "temperature flux_mag", 341);
    assert!(abs.iter().zip(&mag).any(|(a, m)| (a - m).abs() > 1e-3));
    assert!(mag.iter().all(|&m| m >= 0.0));
}

/// `Coordinate System = Axi Symmetric`: the integration weight carries the
/// radius, for both the mass matrix and the right-hand sides.
#[test]
fn axisymmetric_weighting_matches_the_pristine_run() {
    let r = Reference::read("axi.txt");
    let mesh = mesh();
    let m = model(
        &[
            ("Calculate Flux", lg(true)),
            ("Calculate Grad", lg(true)),
            ("Flux Coefficient", st("Heat Conductivity")),
        ],
        &[("Heat Conductivity", Value::real(3.0))],
        true,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();

    check_norm(&sim, 0.07272219701718084);
    for c in 1..=2 {
        for var in ["flux", "grad"] {
            let got = by_node(&sim, &format!("temperature {var} {c}"), 341);
            let want = r.col(&format!("temperature {var}[{}]", c - 1));
            assert_close(&format!("{var} {c}"), &got, &want, TOL);
        }
    }
}

/// A full 2x2 conductivity tensor, `Heat Conductivity(2,2) = 3 0.5 0.5 1.5`.
#[test]
fn full_tensor_coefficient_matches_the_pristine_run() {
    let r = Reference::read("tensor.txt");
    let mesh = mesh();
    let tensor = Value::Real {
        rows: 2,
        cols: 2,
        data: vec![3.0, 0.5, 0.5, 1.5],
    };
    let m = model(
        &[
            ("Calculate Flux", lg(true)),
            ("Calculate Flux Magnitude", lg(true)),
            ("Flux Coefficient", st("Heat Conductivity")),
        ],
        &[("Heat Conductivity", tensor)],
        false,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();

    check_norm(&sim, 0.298231038356245);
    for c in 1..=2 {
        let got = by_node(&sim, &format!("temperature flux {c}"), 341);
        let want = r.col(&format!("temperature flux[{}]", c - 1));
        assert_close(&format!("flux {c}"), &got, &want, TOL);
    }
    let got = by_node(&sim, "temperature flux_mag", 341);
    assert_close("flux_mag", &got, &r.col("temperature flux_mag"), TOL);
}

#[test]
fn no_request_defaults_to_flux_and_missing_coefficient_is_the_identity() {
    // FluxSolver_Init: nothing requested and `Calculate Flux` never given
    // sets it. Without a `Flux Coefficient` the flux is minus the gradient
    // (C = I, after the warning), so it equals minus the grad of the
    // pristine grad run.
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let m = model(&[], &[], false);
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();
    assert_eq!(
        sim.slot(1).unwrap().params.logical("calculate flux"),
        Some(true)
    );
    for c in 1..=2 {
        let got = by_node(&sim, &format!("temperature flux {c}"), 341);
        let want: Vec<f64> = r
            .col(&format!("temperature grad[{}]", c - 1))
            .iter()
            .map(|g| -g)
            .collect();
        assert_close(&format!("flux {c}"), &got, &want, TOL);
    }
}

#[test]
fn explicit_request_for_nothing_creates_no_fields() {
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let m = model(&[("Calculate Flux", lg(false))], &[], false);
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();
    assert!(sim.vars.get("computeflux_temp").is_none());
    assert!(sim.vars.resolve("temperature flux").is_none());
}

#[test]
fn unsupported_cases_are_refused_at_setup() {
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let setup_err = |m: Model| {
        let mut sim = simulation(m, r.col("temperature"));
        run_sim(&mut sim, &mesh).unwrap_err().to_string()
    };

    // fluxsolver3's DG machinery.
    let m = model(
        &[
            ("Calculate Grad", lg(true)),
            ("Discontinuous Galerkin", lg(true)),
        ],
        &[],
        false,
    );
    assert!(setup_err(m).contains("Discontinuous Galerkin"));

    let mut m = model(&[("Calculate Grad", lg(true))], &[], false);
    m.simulation.insert("Coordinate System", st("Polar"));
    assert!(setup_err(m).contains("Coordinate System"));
}

#[test]
fn transient_visits_keep_the_zero_time_order() {
    // `Time derivative order = 0` is added by FluxSolver_Init, so the hidden
    // variable keeps no time history even in a transient simulation, and each
    // timestep's visit gives the same projection.
    let r = Reference::read("flux.txt");
    let mesh = mesh();
    let mut m = model(
        &[
            ("Calculate Flux", lg(true)),
            ("Flux Coefficient", st("Heat Conductivity")),
        ],
        &[("Heat Conductivity", Value::real(3.0))],
        false,
    );
    m.simulation = sec(&[
        ("Simulation Type", st("transient")),
        ("Timestep Intervals", ints(&[2])),
        ("Timestep Sizes", Value::real(1.0)),
        ("Output Intervals", ints(&[0])),
    ]);
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();
    let slot = sim.slot(1).unwrap();
    assert_eq!(slot.time_order, 0);
    assert_eq!(slot.times_visited, 2);
    assert!(sim
        .vars
        .get("computeflux_temp")
        .unwrap()
        .prev_values
        .is_empty());
    check_norm(&sim, 0.19881406507522117);
    let got = by_node(&sim, "temperature flux 2", 341);
    assert_close("flux 2", &got, &r.col("temperature flux[1]"), TOL);
}

#[test]
fn a_variable_keyword_is_used_instead_of_the_hidden_one() {
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let m = model(
        &[
            ("Calculate Grad", lg(true)),
            ("Variable", st("-nooutput gradwork")),
        ],
        &[],
        false,
    );
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();
    assert!(sim.vars.get("computeflux_temp").is_none());
    assert!(sim.vars.get("gradwork").is_some());
    check_norm(&sim, 0.06627135502507372);
    let got = by_node(&sim, "temperature grad 1", 341);
    assert_close("grad 1", &got, &r.col("temperature grad[0]"), TOL);
}

#[test]
fn target_expression_reads_a_material_keyword() {
    // `Target Expression = Potential` projects a material-supplied nodal
    // field (here the constant 2, whose gradient is zero) instead of a
    // solution variable, and names its outputs after the expression.
    let r = Reference::read("grad.txt");
    let mesh = mesh();
    let mut m = model(
        &[
            ("Calculate Grad", lg(true)),
            ("Calculate Grad Magnitude", lg(true)),
        ],
        &[("Potential", Value::real(2.0))],
        false,
    );
    m.solvers[0].0.remove("target variable");
    m.solvers[0].insert("Target Expression", st("Potential"));
    let mut sim = simulation(m, r.col("temperature"));
    run_sim(&mut sim, &mesh).unwrap();
    for name in ["potential grad 1", "potential grad 2", "potential grad_mag"] {
        let got = by_node(&sim, name, 341);
        assert_close(name, &got, &vec![0.0; 341], 1e-12);
    }
    assert!(sim.vars.resolve("temperature grad").is_none());
}
