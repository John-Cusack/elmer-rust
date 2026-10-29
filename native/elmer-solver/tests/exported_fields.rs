//! `-elem` and `-ip` exported variables filled by `Update Exported
//! Variables` (UpdateExportedVariables), from the `case.sif` files of
//! `fem/tests/ElemVariable` and `IpVariable` on their ElmerGrid meshes
//! (`tests/data/expvar/`): `TestVar A = Variable Coordinate, MATC
//! tx(0)+tx(1)` and `TestVar B = Variable Temperature` in the body force,
//! as element means and at the Gauss points. The largest `testvar a` is
//! what the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) runs report
//! through SaveScalars (Solver 2, left out here): 1.9 on elements,
//! 1.95773502692 at the Gauss points. With masks (`ElemVariable2`,
//! `IpVariable3`): the fields live only where the body force sets the mask,
//! and IpVariable3's `Ip Coordinate 1` initial condition at the Gauss points
//! of body 1 peaks at the 0.978867513 its run reports. HeatSolve is the
//! test-local stand-in.

mod support;

use std::path::Path;
use std::rc::Rc;

use elmer_core::{FieldKind, MeshView, SerialComm};
use elmer_solver::{FeBackend, Simulation, SolverRegistry, TextMatc};
use support::poisson::heat_nl;
use support::{sif_model, Capture};

fn run(case: &str) -> (Simulation, support::mesh::FileMesh) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/expvar").join(case);
    let mesh = support::mesh::FileMesh::read(&dir, 2);
    let mut model = sif_model(&dir.join("case.sif"));
    model.solvers.truncate(1);
    let mut reg = SolverRegistry::new();
    reg.register("HeatSolve", "HeatSolver", heat_nl);
    let log = Capture::default();
    let mut sim = Simulation::new(model);
    sim.set_matc(Rc::new(TextMatc(elmer_sif::matc::matc_eval)));
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    sim.setup(&mesh, &reg, &log, &SerialComm).unwrap();
    sim.run(&mesh, &log, &SerialComm).unwrap();
    assert_eq!((sim.reference.test_count, sim.reference.pass_count), (1, 1), "{case}");
    (sim, mesh)
}

fn max(sim: &Simulation, name: &str) -> f64 {
    match sim.vars.get(name) {
        Some(v) => v.values.iter().copied().fold(f64::MIN, f64::max),
        // A component: every dofs-th value.
        None => {
            let (base, j) = name.rsplit_once(' ').unwrap();
            let v = sim.vars.get(base).unwrap();
            let j: usize = j.parse().unwrap();
            v.values.iter().skip(j - 1).step_by(v.dofs).copied().fold(f64::MIN, f64::max)
        }
    }
}

#[test]
fn elemental_means() {
    let (sim, mesh) = run("elem");
    let a = sim.vars.get("testvar a").unwrap();
    assert_eq!((a.kind, a.values.len()), (FieldKind::Elemental, mesh.n_bulk()));
    assert!((max(&sim, "testvar a") - 1.9).abs() < 1e-12, "{}", max(&sim, "testvar a"));
    // TestVar B is the element mean of Temperature (through MATC's printed
    // digits).
    let t = sim.vars.get("temperature").unwrap();
    let e = elmer_core::ElemId(0);
    let nodes = mesh.element(e).nodes;
    let mean = nodes.iter().map(|&n| t.values[t.perm.as_deref().unwrap()[n as usize] as usize]).sum::<f64>() / nodes.len() as f64;
    let b = sim.vars.get("testvar b").unwrap();
    let got = b.values[a.perm.as_deref().unwrap()[0] as usize];
    assert!((got - mean).abs() <= 1e-11 * mean.abs(), "{got} vs {mean}");
}

#[test]
fn gauss_point_values() {
    let (sim, mesh) = run("ip");
    let a = sim.vars.get("testvar a").unwrap();
    assert_eq!(a.kind, FieldKind::Ip);
    let perm = a.perm.as_deref().unwrap();
    assert_eq!(perm.len(), mesh.n_bulk() + mesh.n_boundary() + 1);
    assert_eq!(a.values.len(), 4 * mesh.n_bulk());
    assert!((max(&sim, "testvar a") - 1.95773502692).abs() < 1e-10, "{}", max(&sim, "testvar a"));
}

/// The elements of body `tag`.
fn body_elements(mesh: &support::mesh::FileMesh, tag: i32) -> usize {
    (0..mesh.n_bulk() as u32).filter(|&e| mesh.element(elmer_core::ElemId(e)).tag == tag).count()
}

#[test]
fn masked_elemental_fields() {
    let (sim, mesh) = run("elem2");
    let a = sim.vars.get("testvar a").unwrap();
    let b = sim.vars.get("testvar b").unwrap();
    // `Mask A` is set in body 1's body force, `Mask B` in body 2's. As in
    // CreateElementsPerm, the boundary elements (no active equation) keep
    // the verdict of the last bulk element.
    let last = mesh.element(elmer_core::ElemId(mesh.n_bulk() as u32 - 1)).tag;
    let expect = |body: i32| body_elements(&mesh, body) + if last == body { mesh.n_boundary() } else { 0 };
    assert_eq!((a.values.len(), b.values.len()), (expect(1), expect(2)));
    assert!((max(&sim, "testvar a") - 1.9).abs() < 1e-12, "{}", max(&sim, "testvar a"));
}

#[test]
fn masked_gauss_point_fields_with_initial_conditions() {
    let (sim, mesh) = run("ip3");
    let c = sim.vars.get("ip coordinate").unwrap();
    assert_eq!(c.dofs, 3);
    assert_eq!(c.values.len(), 3 * 4 * body_elements(&mesh, 1));
    assert!((max(&sim, "ip coordinate 1") - 0.978867513).abs() < 1e-9, "{}", max(&sim, "ip coordinate 1"));
    let comp = |j: usize| c.values.iter().skip(j).step_by(3).copied().collect::<Vec<f64>>();
    assert!(comp(1).iter().any(|&y| y != 0.0));
    assert!(comp(2).iter().all(|&z| z == 0.0));
    // Ip Temperature (masked by `body force: ip active`) follows Temperature.
    let t = sim.vars.get("ip temperature").unwrap();
    assert_eq!(t.values.len(), 4 * body_elements(&mesh, 1));
    assert!(t.values.iter().any(|&v| v > 0.0));
}
