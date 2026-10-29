//! `MinimalExample/case.sif` built by hand, standing in for `elmer-sif`.

use std::path::{Path, PathBuf};

use elmer_core::{Dependent, DependentKind, ParamMap, Value};
use elmer_physics_heat::FeProvider;
use elmer_run::{Frontend, LoadedCase};
use elmer_solver::Model;

use super::mesh::FileMesh;
use super::{basis::TestBasis, quad::TestQuadrature};

pub fn data(case: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case)
}

/// A fresh copy of a case directory under the test scratch area.
pub fn scratch_case(case: &str, tag: &str) -> PathBuf {
    let dst = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("elmer-run-{tag}"));
    let _ = std::fs::remove_dir_all(&dst);
    copy_dir(&data(case), &dst);
    dst
}

fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for e in std::fs::read_dir(src).unwrap() {
        let e = e.unwrap();
        let to = dst.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &to);
        } else {
            std::fs::copy(e.path(), to).unwrap();
        }
    }
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

/// The sections of `tests/data/minimal_example/case.sif`, with the string
/// values lowercased as Elmer's parser does outside quotes.
pub fn minimal_example() -> Model {
    Model {
        header: sec(&[("Check Keywords", st("warn")), ("Mesh DB", st("\".\" \"triangle\""))]),
        simulation: sec(&[("Simulation Type", st("steady state")), ("Post File", st("case.vtu"))]),
        solvers: vec![sec(&[
            ("Equation", st("heat equation")),
            ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
            ("Variable", st("temperature")),
            ("Linear System Solver", st("direct")),
            ("Reference Norm", Value::real(2.88675135e-01)),
        ])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Target Bodies", ints(&[1])), ("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(&[("Heat Conductivity", Value::real(1.0)), ("Density", Value::real(1.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[2])), ("Temperature", Value::real(0.0))]),
            sec(&[("Target Boundaries", ints(&[3])), ("Heat Flux", Value::real(1.0))]),
        ],
        ..Model::default()
    }
}

pub fn fe() -> FeProvider {
    || (Box::new(TestBasis), Box::new(TestQuadrature))
}

/// Serves the hand-built model for any `.sif`, with the mesh from
/// `<sif dir>/triangle`, after applying `edit`.
pub struct HandBuilt(pub fn(&mut Model));

impl Frontend for HandBuilt {
    fn fe(&self) -> FeProvider {
        fe()
    }

    fn load(&self, sif: &Path) -> Result<LoadedCase, String> {
        let case_dir = sif.parent().unwrap().to_path_buf();
        let mesh_dir = case_dir.join("triangle");
        let mut model = minimal_example();
        (self.0)(&mut model);
        Ok(LoadedCase { model, mesh: Box::new(FileMesh::read(&mesh_dir, 2)), case_dir, mesh_dir })
    }
}

fn real_vec(v: &[f64]) -> Value {
    Value::Real { rows: v.len(), cols: 1, data: v.to_vec() }
}

/// The linear-system and nonlinear keywords shared by the iterative cases.
fn iterative(extra: &[(&'static str, Value)]) -> Vec<(&'static str, Value)> {
    let mut v = vec![
        ("Equation", st("heat equation")),
        ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
        ("Variable", st("temperature")),
        ("Linear System Solver", st("iterative")),
        ("Linear System Iterative Method", st("bicgstab")),
    ];
    v.extend_from_slice(extra);
    v
}

/// `ThermalPipe2D/case.sif`: three bodies, BCs on four tags each.
pub fn thermal_pipe_2d() -> Model {
    let solver = iterative(&[
        ("Exec Solver", st("always")),
        ("Stabilize", Value::Logical(true)),
        ("Bubbles", Value::Logical(false)),
        ("Lumped Mass Matrix", Value::Logical(false)),
        ("Optimize Bandwidth", Value::Logical(true)),
        ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Convergence Tolerance", Value::real(1.0e-7)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", Value::real(1.0e-3)),
        ("Nonlinear System Relaxation Factor", Value::real(1.0)),
        ("Linear System Max Iterations", ints(&[500])),
        ("Linear System Convergence Tolerance", Value::real(1.0e-10)),
        ("BiCGstabl polynomial degree", ints(&[2])),
        ("Linear System Preconditioning", st("ilu0")),
        ("Linear System ILUT Tolerance", Value::real(1.0e-3)),
        ("Linear System Abort Not Converged", Value::Logical(false)),
        ("Linear System Residual Output", ints(&[20])),
        ("Linear System Precondition Recompute", ints(&[1])),
        ("Reference Norm", Value::real(5.44659991)),
    ]);
    let body = |i: i32, name: &str| {
        sec(&[
            ("Target Bodies", ints(&[i])),
            ("Name", Value::Str(name.into())),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
        ])
    };
    Model {
        header: sec(&[("Check Keywords", st("warn")), ("Mesh DB", st("\".\" \"OnePipe\""))]),
        simulation: sec(&[
            ("Max Output Level", ints(&[5])),
            ("Coordinate System", st("cartesian")),
            ("Coordinate Mapping", ints(&[1, 2, 3])),
            ("Simulation Type", st("steady state")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[1])),
            ("Timestepping Method", st("bdf")),
            ("BDF Order", ints(&[1])),
            ("Post File", st("case.vtu")),
        ]),
        constants: sec(&[
            ("Gravity", real_vec(&[0.0, -1.0, 0.0, 9.82])),
            ("Stefan Boltzmann", Value::real(5.67e-08)),
            ("Permittivity of Vacuum", Value::real(8.8542e-12)),
            ("Boltzmann Constant", Value::real(1.3807e-23)),
            ("Unit Charge", Value::real(1.602e-19)),
        ]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[("Name", Value::Str("HeatEquation".into())), ("Active Solvers", ints(&[1]))])],
        bodies: vec![body(1, "Ground"), body(2, "Grout"), body(3, "Pipe")],
        materials: vec![sec(&[
            ("Name", Value::Str("Ideal".into())),
            ("Heat Conductivity", Value::real(1.0)),
            ("Density", Value::real(1000.0)),
        ])],
        boundary_conditions: vec![
            sec(&[
                ("Target Boundaries", ints(&[9, 10, 11, 12])),
                ("Name", Value::Str("Exterior".into())),
                ("Temperature", Value::real(10.0)),
            ]),
            sec(&[("Target Boundaries", ints(&[5, 6, 7, 8])), ("Name", Value::Str("Grout BC".into()))]),
            sec(&[
                ("Target Boundaries", ints(&[1, 2, 3, 4])),
                ("Name", Value::Str("Pipe BC".into())),
                ("Temperature", Value::real(0.0)),
            ]),
        ],
        ..Model::default()
    }
}

/// `CoordinateScaling/case.sif`. The mesh must be scaled by 0.001 at load.
pub fn coordinate_scaling() -> Model {
    let solver = iterative(&[
        ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Convergence Tolerance", Value::real(1.0e-4)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Relaxation Factor", Value::real(1.0)),
        ("Nonlinear System Convergence Measure", st("solution")),
        ("Linear System Max Iterations", ints(&[500])),
        ("Linear System Convergence Tolerance", Value::real(1.0e-8)),
        ("Linear System Preconditioning", st("ilu1")),
        ("Linear System ILUT Tolerance", Value::real(1.0e-3)),
        ("Linear System Abort Not Converged", Value::Logical(false)),
        ("Linear System Residual Output", ints(&[10])),
        ("Linear System Precondition Recompute", ints(&[1])),
        ("Reference Norm", Value::real(3.937790364340947e-2)),
    ]);
    let mut solver = sec(&solver);
    solver.insert("Equation", st("heatsolver"));
    Model {
        header: sec(&[("Mesh DB", st("\".\" \"square\""))]),
        simulation: sec(&[
            ("Max Output Level", ints(&[5])),
            ("Coordinate System", st("cartesian")),
            ("Coordinate Scaling", Value::real(0.001)),
            ("Simulation Type", st("steady")),
            ("Output Intervals", ints(&[1])),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        solvers: vec![solver],
        equations: vec![sec(&[("Name", Value::Str("Heat".into())), ("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Body Force", ints(&[1])), ("Material", ints(&[1]))])],
        body_forces: vec![sec(&[("Heat Source", Value::real(1.0))])],
        materials: vec![sec(&[
            ("Name", Value::Str("Ideal".into())),
            ("Heat Conductivity", Value::real(1.0)),
            ("Density", Value::real(1.0)),
        ])],
        boundary_conditions: vec![sec(&[("Target Boundaries", ints(&[1])), ("Temperature", Value::real(0.0))])],
        ..Model::default()
    }
}

/// `wildcardsize/case.sif`: `(*)` sizes, so `Heat Conductivity(*) = 1 2` is
/// a diagonal tensor and `Active Solvers(*) = 1 2` names a missing solver.
pub fn wildcardsize() -> Model {
    let solver = iterative(&[
        ("Exec Solver", st("always")),
        ("Linear System Max Iterations", ints(&[350])),
        ("Linear System Convergence Tolerance", Value::real(1.0e-8)),
        ("Linear System Abort Not Converged", Value::Logical(true)),
        ("Linear System Preconditioning", st("ilu0")),
        ("Linear System Residual Output", ints(&[1])),
        ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Convergence Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Max Iterations", ints(&[1])),
        ("Reference Norm", Value::real(3.50844423e-2)),
    ]);
    Model {
        header: sec(&[("Check Keywords", st("warn")), ("Mesh DB", st("\".\" \"square\""))]),
        simulation: sec(&[
            ("Max Output Level", ints(&[5])),
            ("Coordinate System", st("cartesian")),
            ("Coordinate Mapping", ints(&[1, 2, 3])),
            ("Simulation Type", st("steady state")),
            ("Steady State Max Iterations", ints(&[1])),
            ("Output Intervals", ints(&[0])),
        ]),
        constants: sec(&[
            ("Gravity", real_vec(&[0.0, -1.0, 0.0, 9.82])),
            ("Stefan Boltzmann", Value::real(5.67e-08)),
        ]),
        solvers: vec![sec(&solver)],
        equations: vec![sec(&[("Name", Value::Str("Equations".into())), ("Active Solvers", ints(&[1, 2]))])],
        bodies: vec![sec(&[
            ("Name", Value::Str("Body".into())),
            ("Body Force", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
        ])],
        body_forces: vec![sec(&[("Name", Value::Str("BodyForce".into())), ("Heat Source", Value::real(1.0))])],
        materials: vec![sec(&[
            ("Name", Value::Str("Material".into())),
            ("Density", Value::real(1.0)),
            ("Heat Conductivity", real_vec(&[1.0, 2.0])),
        ])],
        boundary_conditions: vec![sec(&[
            ("Name", Value::Str("Zero".into())),
            ("Target Boundaries", ints(&[1, 2, 3])),
            ("Temperature", Value::real(0.0)),
        ])],
        ..Model::default()
    }
}

/// The `Solver 2` SaveScalars section of the heateq_* tests: never executed.
fn save_scalars_never() -> ParamMap {
    sec(&[
        ("Exec Solver", st("never")),
        ("Equation", st("save scalars")),
        ("Procedure", st("\"SaveData\" \"SaveScalars\"")),
        ("Variable 1", st("temp")),
        ("Operator 1", st("int mean")),
        ("Variable 2", st("temp")),
        ("Operator 2", st("max")),
        ("Variable 3", st("time")),
        ("Filename", st("f.dat")),
    ])
}

/// The `Solver 1` of the heateq_* tests: direct (UMFPACK), up to
/// `max_nonlinear` Picard iterations.
fn heateq_solver(max_nonlinear: i32, extra: &[(&'static str, Value)]) -> ParamMap {
    let mut v = vec![
        ("Equation", st("heat")),
        ("Variable", st("temp")),
        ("Procedure", st("\"HeatSolve\" \"HeatSolver\"")),
        ("Linear System Solver", st("direct")),
        ("Linear System Direct Method", st("umfpack")),
        ("Steady State Convergence Tolerance", Value::real(1.0e-5)),
        ("Stabilize", Value::Logical(true)),
        ("Nonlinear System Convergence Tolerance", Value::real(1.0e-5)),
        ("Nonlinear System Max Iterations", ints(&[max_nonlinear])),
        ("Nonlinear System Relaxation Factor", Value::real(1.0)),
        ("Nonlinear System Newton After Iterations", ints(&[3])),
        ("Nonlinear System Newton After Tolerance", Value::real(1.0e-2)),
    ];
    v.extend_from_slice(extra);
    sec(&v)
}

fn heateq_common(simulation: ParamMap, solver: ParamMap, material: ParamMap, heat_source: f64) -> Model {
    Model {
        header: sec(&[("Check Keywords", st("warn")), ("Mesh DB", st("\".\" \"angle\""))]),
        simulation,
        constants: sec(&[
            ("Gravity", real_vec(&[0.0, -1.0, 0.0, 9.82])),
            ("Stefan Boltzmann", Value::real(5.67e-08)),
        ]),
        solvers: vec![solver, save_scalars_never()],
        equations: vec![sec(&[("Name", Value::Str("Equation1".into())), ("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Name", Value::Str("Body1".into())),
            ("Body Force", ints(&[1])),
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
        ])],
        body_forces: vec![sec(&[("Heat Source", Value::real(heat_source))])],
        materials: vec![material],
        boundary_conditions: vec![sec(&[("Target Boundaries", ints(&[1])), ("Temp", Value::real(0.0))])],
        ..Model::default()
    }
}

/// `heateq_bdf2/TempDist.sif` (`order` 2) and `heateq_bdf3/TempDist.sif`
/// (`order` 3): BDF over 5 steps of 0.1 then 3 of 0.2.
pub fn heateq_bdf(order: i32, reference: f64) -> Model {
    let simulation = sec(&[
        ("Max Output Level", ints(&[5])),
        ("Coordinate System", st("cartesian 2d")),
        ("Coordinate Mapping", ints(&[1, 2, 3])),
        ("Steady State Max Iterations", ints(&[1])),
        ("Simulation Type", st("transient")),
        ("Timestepping Method", st("bdf")),
        ("BDF Order", ints(&[order])),
        ("Timestep Sizes", real_vec(&[0.1, 0.2])),
        ("Timestep Intervals", ints(&[5, 3])),
        // `Output Intervals = 1` is given first; the later `(2) = 1 1` wins.
        ("Output Intervals", ints(&[1, 1])),
    ]);
    let solver = heateq_solver(
        10,
        &[("Optimize Bandwidth", Value::Logical(false)), ("Reference Norm", Value::real(reference))],
    );
    let material = sec(&[
        ("Name", Value::Str("Material1".into())),
        ("Density", Value::real(1.0)),
        ("Heat Conductivity", Value::real(1.0)),
        ("Heat Capacity", Value::real(10.0)),
    ]);
    heateq_common(simulation, solver, material, 1.0)
}

/// `heateq_spline_table/TempDist.sif`: `Heat Conductivity = Variable Temp,
/// Real Cubic` with the `include cond.dat` table read from `dir`.
pub fn heateq_spline_table(dir: &Path) -> Model {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for line in std::fs::read_to_string(dir.join("cond.dat")).unwrap().lines() {
        let v: Vec<f64> = line.split_whitespace().map(|t| t.parse().unwrap()).collect();
        if let [a, b] = v[..] {
            x.push(a);
            y.push(b);
        }
    }
    let conductivity = Value::Dependent(Dependent {
        variables: vec!["temp".into()],
        kind: DependentKind::Table { x, y, cubic: true },
        rows: 1,
        cols: 1,
    });
    let simulation = sec(&[
        ("Max Output Level", ints(&[5])),
        ("Coordinate System", st("cartesian")),
        ("Coordinate Mapping", ints(&[1, 2, 3])),
        ("Steady State Max Iterations", ints(&[1])),
        ("Output Intervals", ints(&[1])),
        ("Simulation Type", st("steady state")),
    ]);
    let solver = heateq_solver(100, &[("Reference Norm", Value::real(1.5131))]);
    let material = sec(&[
        ("Name", Value::Str("Material".into())),
        ("Density", Value::real(1.0)),
        ("Heat Conductivity", conductivity),
        ("Heat Capacity", Value::real(1.0)),
    ]);
    heateq_common(simulation, solver, material, 10.0)
}
