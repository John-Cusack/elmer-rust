//! The p basis and rules of boundary faces of a 3D mesh (`PDefs % isEdge`
//! quads, oriented by `getSquareFaceDirection`) and `GaussPoints(Element,
//! np)`, on the unit-cube `brick_p3` mesh.

mod support;

use elmer_core::{BasisEval, ElemId, ElementBasis, MeshView, ParamMap, Value};
use elmer_pdofs::{ModelFacts, PSetup};
use support::{data_dir, TestMesh};

fn setup(def: &str, serendipity: bool) -> (TestMesh, PSetup) {
    let mesh = TestMesh::read(&data_dir("brick_p3"), 3);
    let mut solver = ParamMap::new();
    solver.insert("Element", Value::Str(def.into()));
    let mut simulation = ParamMap::new();
    simulation.insert("Serendipity P Elements", Value::Logical(serendipity));
    let n_bulk = mesh.n_bulk();
    let active = |e: ElemId| (e.0 as usize) < n_bulk;
    let mut matc = |_: &str, _: [f64; 4]| unreachable!("no MATC degree");
    let setup = PSetup::new(
        &mesh,
        &solver,
        &simulation,
        ModelFacts {
            n_bodies: 1,
            n_equations: 1,
            body_equation: &|_| Some(1),
            equation_element: &|_| None,
            equation_solver_element: &|_| None,
            body_solver_element: &|_| None,
            equation_activates: &|_| Some(true),
            active: &active,
            boundary_body: &|_| 0,
            matc_degree: &mut matc,
        },
        1,
    )
    .unwrap();
    (mesh, setup)
}

#[test]
fn boundary_faces_carry_edge_and_face_functions() {
    for (def, sd, nd) in [
        ("p:2", false, 4 + 4 + 1),
        ("p:3", false, 4 + 4 * 2 + 4),
        ("p:4", true, 4 + 4 * 3 + 1),
    ] {
        let (mesh, setup) = setup(def, sd);
        let mut slots = Vec::new();
        let mut be = BasisEval::default();
        let mut area = 0.0;
        for e in (mesh.n_bulk()..mesh.n_bulk() + mesh.n_boundary()).map(|i| ElemId(i as u32)) {
            assert!(setup.is_p_element(e), "{def}: element {}", e.0);
            setup.element_slots(e, &mut slots);
            assert_eq!(slots.len(), nd, "{def}: boundary dofs");
            let spec = setup.layout().element_spec(e.0 as usize).unwrap();
            assert!(
                spec.bubble_face_dir.is_some(),
                "{def}: face bubbles are oriented"
            );
            assert_eq!(spec.n_basis().unwrap(), nd, "{def}: basis functions");
            let basis = setup.element_basis(e).unwrap();
            let rule = setup.gauss_rule(e).unwrap();
            let xyz: Vec<[f64; 3]> = setup
                .element_nodes(e)
                .iter()
                .map(|&n| mesh.coord(n))
                .collect();
            for t in 0..rule.len() {
                basis
                    .eval(
                        mesh.element(e).etype,
                        &xyz,
                        [rule.u[t], rule.v[t], rule.w[t]],
                        &mut be,
                    )
                    .unwrap();
                area += rule.s[t] * be.det_j;
            }
        }
        assert!((area - 6.0).abs() < 1e-12, "{def}: surface area {area}");
    }
}

#[test]
fn np_point_rules_follow_gauss_points_with_np() {
    let (mesh, setup) = setup("p:2", false);
    for np in [8, 27, 64] {
        let r = setup.gauss_rule_np(ElemId(0), np).unwrap();
        assert_eq!(r.len(), np);
        // Reference brick [-1, 1]^3.
        let vol: f64 = r.s.iter().sum();
        assert!((vol - 8.0).abs() < 1e-13, "np {np}: {vol}");
    }
    let b = ElemId(mesh.n_bulk() as u32);
    assert_eq!(setup.gauss_rule_np(b, 16).unwrap().len(), 16);
}
