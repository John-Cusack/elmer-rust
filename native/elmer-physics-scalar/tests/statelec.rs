//! StatElecSolve end to end against the pristine ElmerSolver (`data/README.md`
//! records the runs): `fem/tests/elstat` on a coarser mesh, with the field
//! recovery, and `fem/tests/elstat_source`.

mod support;

use elmer_core::Value;
use elmer_solver::Model;
use support::cases::*;
use support::meshio::read_elmer_mesh;
use support::simmesh::SimMesh;

fn direct() -> Vec<(&'static str, Value)> {
    vec![
        ("Linear System Solver", st("Direct")),
        ("Linear System Direct Method", st("UMFPack")),
    ]
}

fn elstat_solver(extra: &[(&str, Value)]) -> elmer_core::ParamMap {
    let mut v = vec![
        ("Equation", st("Stat Elec Solver")),
        ("Variable", st("Potential")),
        ("Variable DOFs", ints(&[1])),
        ("Procedure", st("\"StatElecSolve\" \"StatElecSolver\"")),
        ("Nonlinear System Max Iterations", ints(&[1])),
    ];
    v.extend(direct());
    v.extend(extra.iter().cloned());
    sec(&v)
}

/// fem/tests/elstat: 1e-6 coordinate scaling, vacuum permittivity 8.8542e-12,
/// 1000 V on boundary 3 and 0 V on boundary 4.
fn elstat_model(extra: &[(&str, Value)]) -> Model {
    Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 3D")),
            ("Simulation Type", st("Steady State")),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        constants: sec(&[("Permittivity Of Vacuum", real(8.8542e-12))]),
        solvers: vec![elstat_solver(extra)],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[("Equation", ints(&[1])), ("Material", ints(&[1]))])],
        materials: vec![sec(&[("Relative Permittivity", real(1.0))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[4])), ("Potential", real(0.0))]),
            sec(&[
                ("Target Boundaries", ints(&[3])),
                ("Potential", real(1.0e3)),
            ]),
        ],
        ..Model::default()
    }
}

fn elstat_mesh() -> SimMesh {
    let mut m = read_elmer_mesh(&data("elstat_coarse"), 3);
    scale_coords(&mut m, 1.0e-6);
    SimMesh::new(m)
}

#[test]
fn elstat_potential_and_energy_match_elmer() {
    let mesh = elstat_mesh();
    let sim = run(
        elstat_model(&[("Calculate Electric Energy", Value::Logical(true))]),
        &mesh,
    )
    .unwrap();
    let n = mesh.mesh.coords.len();
    let cols = read_columns(&data("elstat_coarse").join("fields.dat"));
    let d = rel_max_diff(&nodal(&sim, "potential", 0, n), &cols[0]);
    // The pristine run is BiCGStab/ILU1 to 1e-10; this one is a direct solve.
    assert!(d < 1e-8, "potential: {d:e}");
    assert_rel(norm(&sim, "potential"), 631.9878246439413, 1e-8);

    let p = &sim.model.simulation;
    use elmer_core::Params;
    assert_rel(
        p.real("res: electric energy").unwrap(),
        1.073050209726e-10,
        1e-8,
    );
    assert_rel(p.real("res: potential difference").unwrap(), 1.0e3, 1e-14);
    assert_rel(
        p.real("res: effective capacitance").unwrap(),
        2.146100419453e-16,
        1e-8,
    );
    let e = nodal(&sim, "electric energy density", 0, n);
    assert!(rel_max_diff(&e, &cols[1]) < 1e-8);
}

#[test]
fn elstat_fields_match_elmer() {
    let mesh = elstat_mesh();
    let extra = [
        ("Calculate Electric Energy", Value::Logical(true)),
        ("Calculate Electric Field", Value::Logical(true)),
        ("Calculate Electric Flux", Value::Logical(true)),
        ("Calculate Nodal Energy", Value::Logical(true)),
        ("Calculate Nodal Flux", Value::Logical(true)),
    ];
    let sim = run(elstat_model(&extra), &mesh).unwrap();
    let n = mesh.mesh.coords.len();
    let cols = read_columns(&data("elstat_coarse").join("fields.dat"));
    // Column order: potential, energy density, nodal energy, field 1-3, flux 1-3, nodal flux 1-3.
    let checks: [(&str, usize, usize); 9] = [
        ("electric energy density", 0, 1),
        ("nodal energy density", 0, 2),
        ("electric field", 0, 3),
        ("electric field", 1, 4),
        ("electric field", 2, 5),
        ("electric flux", 0, 6),
        ("electric flux", 1, 7),
        ("electric flux", 2, 8),
        ("nodal electric flux", 2, 11),
    ];
    for (var, k, col) in checks {
        let d = rel_max_diff(&nodal(&sim, var, k, n), &cols[col]);
        assert!(d < 1e-8, "{var} {k}: {d:e}");
    }
}

#[test]
fn elstat_source_matches_elmer() {
    // fem/tests/elstat_source: charge density 2e-6 in a 2D quad8 mesh, all
    // boundaries grounded.
    let mesh = SimMesh::new(read_elmer_mesh(&data("elstat_source"), 2));
    let model = Model {
        simulation: sec(&[
            ("Coordinate System", st("Cartesian 2D")),
            ("Simulation Type", st("Steady")),
            ("Steady State Max Iterations", ints(&[1])),
        ]),
        constants: sec(&[("Permittivity Of Vacuum", real(8.8542e-12))]),
        solvers: vec![elstat_solver(&[])],
        equations: vec![sec(&[("Active Solvers", ints(&[1]))])],
        bodies: vec![sec(&[
            ("Equation", ints(&[1])),
            ("Material", ints(&[1])),
            ("Body Force", ints(&[1])),
        ])],
        materials: vec![sec(&[("Relative Permittivity", real(1.0))])],
        body_forces: vec![sec(&[("Charge density", real(2e-6))])],
        boundary_conditions: vec![sec(&[
            ("Target Boundaries", ints(&[1, 2, 3, 4, 5, 6])),
            ("Potential", real(0.0)),
        ])],
        ..Model::default()
    };
    let sim = run(model, &mesh).unwrap();
    let want = &read_columns(&data("elstat_source").join("solution.dat"))[0];
    let d = rel_max_diff(&nodal(&sim, "potential", 0, mesh.mesh.coords.len()), want);
    assert!(d < 1e-10, "{d:e}");
    // The sif's Reference Norm is 795.505441096.
    assert_rel(norm(&sim, "potential"), 795.5054410960914, 1e-11);
}
