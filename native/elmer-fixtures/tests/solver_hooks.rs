//! The solver procedures and assembly routines: what `registry::register`
//! puts under which name, and the `Poisson` element matrix.


use approx::assert_abs_diff_eq;
use elmer_core::{ElementType, LocalSystem};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_fixtures::poisson_solver::local_matrix;
use elmer_fixtures::registry;
use elmer_solver::SolverRegistry;

#[test]
fn cases_get_their_own_procedures_under_the_sif_names() {
    let mut reg = SolverRegistry::new();
    assert_eq!(registry::register("CoupledPoisson4", &mut reg), 2);
    // `File "./PoissonAssembly" "BulkAssembly"`: directory and case are ignored.
    assert!(reg.assembly("./PoissonAssembly", "bulkassembly").is_some());
    assert!(reg.assembly("PoissonAssembly", "BoundaryAssembly").is_some());
    assert!(reg.get("Poisson", "PoissonSolver").is_none(), "CoupledPoisson4 has no Poisson.f90");

    let mut reg = SolverRegistry::new();
    assert_eq!(registry::register("1dtests", &mut reg), 1);
    assert!(reg.get("Poisson", "PoissonSolver").is_some());
    assert_eq!(registry::hook_procedures("ExtrudedMesh"), vec![("Poisson", "PoissonSolver")]);

    // PoissonPFEM's own `Poisson.f90`, ported in `elmer-pdofs`.
    assert_eq!(registry::hook_procedures("PoissonPFEM"), vec![("Poisson", "PoissonSolver")]);

    // Same file names, different sources: not ported, so not registered.
    for case in ["CurvedBndryPFEM", "SolverSpecificMeshes", "BDM2D"] {
        assert_eq!(registry::register(case, &mut SolverRegistry::new()), 0, "{case}");
    }
}

#[test]
fn every_routine_built_from_the_registry_is_a_fresh_instance() {
    let mut reg = SolverRegistry::new();
    registry::register("CoupledPoisson1", &mut reg);
    let f = reg.assembly("PoissonAssembly", "BulkAssembly").unwrap();
    let params = elmer_core::ParamMap::new();
    // Two builds, two `SAVE`d `Visited` flags: building succeeds twice and
    // the boxes are distinct objects.
    let a = f(&params).unwrap();
    let b = f(&params).unwrap();
    assert!(!std::ptr::addr_eq(&*a as *const _, &*b as *const _));
}

#[test]
fn poisson_local_matrix_is_the_p1_laplacian_with_lumped_load() {
    // Unit right triangle, constant source 3: K = area * g_p . g_q,
    // F = 3 * area / 3 per node.
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mut local = LocalSystem::zeros(0);
    local_matrix(ElementType(303), &xyz, &[3.0; 3], &LagrangeBasis::new(2), &ElmerQuadrature, &mut local).unwrap();
    let want = [1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5];
    for (g, w) in local.stiff.iter().zip(want) {
        assert_abs_diff_eq!(*g, w, epsilon = 1e-15);
    }
    for f in &local.force {
        assert_abs_diff_eq!(*f, 0.5, epsilon = 1e-15);
    }
}

#[test]
fn poisson_solver_refuses_non_nodal_elements() {
    let mut reg = SolverRegistry::new();
    registry::register("ExtrudedMesh", &mut reg);
    let f = reg.get("Poisson", "PoissonSolver").unwrap();
    let mut p = elmer_core::ParamMap::new();
    p.insert("Element", elmer_core::Value::Str("p:2".into()));
    let err = f(&p).err().expect("p:2 is refused").to_string();
    assert!(err.contains("Element = p:2"), "{err}");
}

#[test]
fn tresca_is_the_spread_of_the_principal_stresses() {
    use elmer_fixtures::tresca::{symmetric_eigenvalues, tresca};
    // 2D: xx yy zz xy; principal values of [[3,1],[1,3]] are 2 and 4, zz = 0.
    assert_abs_diff_eq!(tresca(&[3.0, 3.0, 0.0, 1.0]), 4.0, epsilon = 1e-14);
    // 3D, a full tensor against its characteristic polynomial.
    let a = [[4.0, 1.0, -2.0], [1.0, 2.0, 0.5], [-2.0, 0.5, 3.0]];
    let e = symmetric_eigenvalues(a);
    let tr: f64 = e.iter().sum();
    assert_abs_diff_eq!(tr, 9.0, epsilon = 1e-13);
    for l in e {
        let m = |i: usize, j: usize| a[i][j] - if i == j { l } else { 0.0 };
        let det = m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
            + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
        assert_abs_diff_eq!(det, 0.0, epsilon = 1e-12);
    }
    assert_abs_diff_eq!(tresca(&[4.0, 2.0, 3.0, 1.0, 0.5, -2.0]), e.iter().fold(f64::MIN, |a, &b| a.max(b)) - e.iter().fold(f64::MAX, |a, &b| a.min(b)), epsilon = 1e-14);
}
