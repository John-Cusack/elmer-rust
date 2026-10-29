//! Element definitions outside the solver section, as `NonNodalElements`
//! reads them (MeshLoad.F90:787): an Equation `Element` (every solver
//! without its own), an Equation `Element{id}` and a Body's `Solver id:
//! Element`. On the `varying_p` mesh, p = 4 in body 1 and 8 in body 2
//! given each of these ways must give the pristine Perm of
//! `p-FEM_with_varying_p` (where the solver section has `p:%bodywisep`);
//! the pristine `p-FEM_two_solvers` (body-wise `Solver 1: Element` p:4 /
//! p:8 on the same mesh) saves the same Perm (289 slots, 189 used).

mod support;

use std::sync::{Arc, Mutex};

use elmer_core::{ElemId, MeshView, ParamMap, Value, NO_DOF};
use elmer_pdofs::PDofsBackend;
use elmer_solver::{DofsBackend, Model};
use support::{data_dir, read_perm, TestMesh};

fn sec(entries: &[(&str, Value)]) -> ParamMap {
    let mut p = ParamMap::new();
    for (k, v) in entries {
        p.insert(k, v.clone());
    }
    p
}

fn st(s: &str) -> Value {
    Value::Str(s.into())
}

/// The model of the case with the given solver, equation and body extras.
fn model(
    solver: &[(&str, Value)],
    equation: &[(&str, Value)],
    bodies: [&[(&str, Value)]; 2],
) -> Model {
    let mut eq = vec![("Active Solvers", Value::Integer(vec![1]))];
    eq.extend(equation.iter().cloned());
    let body = |extra: &[(&str, Value)]| {
        let mut b = vec![("Equation", Value::Integer(vec![1]))];
        b.extend(extra.iter().cloned());
        sec(&b)
    };
    Model {
        simulation: sec(&[("Serendipity P Elements", Value::Logical(false))]),
        solvers: vec![sec(solver), sec(&[("Equation", st("other"))])],
        equations: vec![sec(&eq)],
        bodies: vec![body(bodies[0]), body(bodies[1])],
        ..Model::default()
    }
}

/// `bodywisep(cx)`: 8 in body 2, else 4. Records the commands it ran.
fn backend(seen: Arc<Mutex<Vec<String>>>) -> PDofsBackend {
    PDofsBackend::new().with_matc(Arc::new(move |cmd: &str, cx: [f64; 4]| {
        seen.lock().unwrap().push(cmd.to_owned());
        Ok(if cx[0] > 1.0 { 8.0 } else { 4.0 })
    }))
}

/// The Perm the backend builds for solver `id` (1-based) with 5 dofs, as
/// `perm.dat` lists it (1-based rows, 0 absent); `None` for a nodal solver.
fn perm_of(b: &PDofsBackend, m: &Model, id: usize) -> Option<Vec<u32>> {
    let mesh = TestMesh::read(&data_dir("varying_p"), 2);
    let active: Vec<ElemId> = (0..mesh.n_bulk() as u32).map(ElemId).collect();
    let layout = b
        .for_solver_id(&mesh, m, id, &m.solvers[id - 1], &active, 5)?
        .unwrap();
    Some(
        layout
            .perm()
            .perm
            .iter()
            .map(|&p| if p == NO_DOF { 0 } else { p + 1 })
            .collect(),
    )
}

fn pristine() -> Vec<u32> {
    read_perm(&data_dir("varying_p").join("perm.dat"))
}

#[test]
fn solver_section_definition_is_the_reference() {
    let seen = Arc::new(Mutex::new(Vec::new()));
    let m = model(&[("Element", st("p:%bodywisep"))], &[], [&[], &[]]);
    assert_eq!(perm_of(&backend(seen.clone()), &m, 1).unwrap(), pristine());
    assert!(seen.lock().unwrap().iter().all(|c| c == "bodywisep(cx)"));
}

#[test]
fn equation_element_applies_to_solvers_without_their_own() {
    let seen = Arc::new(Mutex::new(Vec::new()));
    let m = model(&[], &[("Element", st("p:%bodywisep"))], [&[], &[]]);
    assert_eq!(perm_of(&backend(seen.clone()), &m, 1).unwrap(), pristine());
    assert!(!seen.lock().unwrap().is_empty());
    // Elmer applies it to every solver without an Element, active or not.
    assert!(perm_of(&backend(seen), &m, 2).is_some());
}

#[test]
fn equation_solver_element_applies_to_that_solver() {
    let seen = Arc::new(Mutex::new(Vec::new()));
    let m = model(&[], &[("Element{1}", st("p:%bodywisep"))], [&[], &[]]);
    // SetGlobalBubblesFlag never reads Element{id}: the bubbles stay local,
    // so the Perm holds the 9 nodes and 12 x 7 edge slots, 9 + 5 x 3 + 7 x 7
    // of them used, and each element condenses its own bubbles.
    let mesh = TestMesh::read(&data_dir("varying_p"), 2);
    let active: Vec<ElemId> = (0..mesh.n_bulk() as u32).map(ElemId).collect();
    let b = backend(seen.clone());
    let layout = b
        .for_solver_id(&mesh, &m, 1, &m.solvers[0], &active, 5)
        .unwrap()
        .unwrap();
    let perm = layout.perm();
    assert_eq!((perm.perm.len(), perm.n_slots), (93, 73));
    let bubbles: Vec<usize> = active.iter().map(|&e| layout.local_bubbles(e)).collect();
    assert_eq!(bubbles, [9, 9, 49, 49]);
    assert!(seen.lock().unwrap().iter().all(|c| c == "bodywisep(cx)"));
    // Not solver 2's.
    assert!(perm_of(&backend(seen), &m, 2).is_none());
}

#[test]
fn body_wise_solver_elements_set_each_body() {
    let seen = Arc::new(Mutex::new(Vec::new()));
    let m = model(
        &[],
        &[],
        [
            &[("Solver 1: Element", st("p:4"))],
            &[("Solver 1: Element", st("p:8"))],
        ],
    );
    assert_eq!(perm_of(&backend(seen.clone()), &m, 1).unwrap(), pristine());
    assert!(seen.lock().unwrap().is_empty(), "no MATC");
    // Without the solver id only the solver section counts: nodal.
    let mesh = TestMesh::read(&data_dir("varying_p"), 2);
    let active: Vec<ElemId> = (0..mesh.n_bulk() as u32).map(ElemId).collect();
    assert!(backend(seen)
        .for_solver(&mesh, &m, &m.solvers[0], &active, 5)
        .is_none());
}
