//! `fem/tests/PoissonPFEM` end to end against the pristine captures in
//! `tests/data/poisson_pfem`: Perm, matrix pattern, the assembled system
//! before and after Dirichlet conditions, and the solution.
//!
//! The kernel is the test's own `Poisson.f90` (`LocalMatrix` with the
//! nodal load, `LCondensate` a no-op since the bubbles are global).

mod support;

use elmer_assembly::{
    apply_dirichlet, CsrAssembler, DirichletBc, DirichletMode, DirichletSet, PermOptions,
};
use elmer_core::{
    BasisEval, CsrMatrix, DirectSolver, ElemId, ElementBasis, ElementType, GlobalAssembler,
    LocalSystem, MeshView,
};
use elmer_direct::FaerLu;
use elmer_pdofs::dirichlet::constrain_edge_dofs;
use elmer_pdofs::{
    create_p_system, initial_max_defs, DofLayout, ElementDefs, PMesh, PMeshOptions, PSystem,
    SolverElements,
};
use support::{data_dir, read_matrix, read_perm, read_rhs, read_values, TestMesh};

const DEF: &str = "p:2";

struct Case {
    mesh: TestMesh,
    pm: PMesh,
    defs: ElementDefs,
    sys: PSystem,
}

fn setup() -> Case {
    let mesh = TestMesh::read(&data_dir("poisson_pfem"), 2);
    let mut max = initial_max_defs();
    let mut defs = ElementDefs::from_solver_element(Some(DEF), 2, false, &mut max).unwrap();
    let none = |_: usize| None;
    let mut pm = {
        let mut solvers = [SolverElements {
            defs: &mut defs,
            element: Some(DEF),
            dg: false,
            body_element: &none,
        }];
        PMesh::new(
            &mesh,
            &mut solvers,
            &mut max,
            PMeshOptions {
                // case.sif: Serendipity P Elements = False
                serendipity: false,
                need_edges: None,
                n_bodies: 2,
                body_has_equation: &|_| true,
                matc_degree: &mut |_, _| unreachable!("no MATC degree"),
                boundary_body: &|_| 0,
            },
        )
        .unwrap()
    };
    let n_bulk = pm.n_bulk;
    // Both bodies run the solver; boundary elements have no body.
    let active = move |e: usize| e < n_bulk;
    // `p:` in the definition: Bubbles in Global System defaults to true.
    let sys = create_p_system(
        &mut pm,
        &defs,
        true,
        false,
        &active,
        PermOptions::default(),
        1,
    );
    Case {
        mesh,
        pm,
        defs,
        sys,
    }
}

fn layout(c: &Case) -> DofLayout<'_> {
    DofLayout::new(&c.pm, &c.defs, true)
}

#[test]
fn perm_and_bandwidth_match() {
    let c = setup();
    assert_eq!(c.pm.edges.len(), 220);
    assert_eq!((c.pm.max_edge_dofs, c.pm.max_bdofs), (1, 1));
    // 121 nodes + 220 edges + 100 bubbles, every slot used.
    assert_eq!(c.sys.perm.perm.len(), 441);
    assert_eq!(c.sys.perm.n_dofs, 441);
    assert_eq!(c.sys.bandwidth.initial, 342);
    assert_eq!(c.sys.bandwidth.optimized, Some(85));
    assert!(c.sys.bandwidth.accepted);
    let elmer = read_perm(&data_dir("poisson_pfem").join("perm.dat"));
    let ours: Vec<u32> = c.sys.perm.perm.iter().map(|&p| p + 1).collect();
    assert_eq!(ours, elmer);
}

/// Bulk element systems glued into the pattern, as `PoissonSolver` does.
fn assemble(c: &Case) -> (CsrMatrix<f64>, Vec<f64>) {
    let lay = layout(c);
    let mut a = c.sys.matrix.clone();
    let mut b = vec![0.0; a.nrows];
    let mut asm = CsrAssembler::new(&mut a, &mut b);
    let mut rows = Vec::new();
    let mut be = BasisEval::default();
    for e in 0..c.pm.n_bulk {
        let el = &c.pm.elements[e];
        lay.element_rows(e, &c.sys.perm.perm, &mut rows);
        let nd = rows.len();
        assert_eq!(nd, 9, "element {e}: 4 nodes + 4 edges + 1 bubble");
        let basis = lay.element_basis(e).unwrap();
        let rule = lay.gauss_rule(e).unwrap();
        assert_eq!(rule.len(), 9);
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| c.pm.coords[n as usize]).collect();
        let n = el.nodes.len();
        let mut local = LocalSystem::zeros(nd);
        for t in 0..rule.len() {
            basis
                .eval(
                    ElementType(el.code),
                    &xyz,
                    [rule.u[t], rule.v[t], rule.w[t]],
                    &mut be,
                )
                .unwrap();
            assert_eq!(be.basis.len(), nd);
            let w = rule.s[t] * be.det_j;
            // Source = 1 at the nodes, interpolated with the nodal basis.
            let load: f64 = be.basis[..n].iter().sum();
            for i in 0..nd {
                for j in 0..nd {
                    let g = (0..3)
                        .map(|k| be.dbasis_dx[i][k] * be.dbasis_dx[j][k])
                        .sum::<f64>();
                    local.stiff[i * nd + j] += w * g;
                }
                local.force[i] += w * load * be.basis[i];
            }
        }
        asm.add_local(&rows, &local);
    }
    (a, b)
}

fn compare_matrix(a: &CsrMatrix<f64>, file: &str) {
    let elmer = read_matrix(&data_dir("poisson_pfem").join(file));
    assert_eq!(elmer.len(), a.nnz(), "{file}: nnz");
    let scale = elmer.iter().fold(0.0f64, |m, e| m.max(e.2.abs()));
    let mut k = 0;
    let mut worst = 0.0f64;
    for i in 0..a.nrows {
        for p in a.row_ptr[i]..a.row_ptr[i + 1] {
            let (r, col, v) = elmer[k];
            assert_eq!((r, col), (i, a.col_idx[p] as usize), "{file}: entry {k}");
            worst = worst.max((a.values[p] - v).abs());
            k += 1;
        }
    }
    assert!(
        worst <= 1e-13 * scale,
        "{file}: max |diff| {worst:e}, scale {scale:e}"
    );
}

fn compare_rhs(b: &[f64], file: &str) {
    let elmer = read_rhs(&data_dir("poisson_pfem").join(file));
    assert_eq!(elmer.len(), b.len());
    let scale = elmer.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    for (i, (&ours, e)) in b.iter().zip(&elmer).enumerate() {
        let diff = (ours - e.unwrap_or(0.0)).abs();
        assert!(diff <= 1e-13 * scale, "{file}: row {i}: {ours:e} vs {e:?}");
    }
}

#[test]
fn assembled_system_matches_before_dirichlet() {
    let c = setup();
    let (a, b) = assemble(&c);
    compare_matrix(&a, "pre_a.dat");
    compare_rhs(&b, "pre_b.dat");
}

fn dirichlet(c: &Case, a: &mut CsrMatrix<f64>, b: &mut [f64]) {
    let lay = layout(c);
    let mut set = DirichletSet::new(a.nrows);
    // Boundary Condition 1: Target Boundaries(4) = 1 2 3 4, Potential = 0.
    let on_bc = |e: usize| e >= c.pm.n_bulk && (1..=4).contains(&c.pm.elements[e].boundary_tag);
    constrain_edge_dofs(&lay, &c.sys.perm, 0, 1, &on_bc, true, &mut set).unwrap();
    let on_bc_id = |e: ElemId| on_bc(e.0 as usize);
    let bcs = [DirichletBc {
        elements: &on_bc_id,
        value: &|_, _| 0.0,
        condition: None,
    }];
    // One DOF per node: the node slots are the node indices.
    set.set_boundaries(&c.mesh, &c.sys.perm, 0, 1, &bcs, false);
    assert_eq!(set.constrained.iter().filter(|&&x| x).count(), 40 + 40);
    apply_dirichlet(a, b, &set.to_list(), DirichletMode::default()).unwrap();
}

#[test]
fn system_and_solution_match_after_dirichlet() {
    let c = setup();
    assert_eq!(c.mesh.n_boundary(), 40);
    let (mut a, mut b) = assemble(&c);
    dirichlet(&c, &mut a, &mut b);
    compare_matrix(&a, "post_a.dat");
    compare_rhs(&b, "post_b.dat");

    elmer_direct::set_single_thread();
    let mut lu = FaerLu::new();
    lu.factor(&a).unwrap();
    let mut x = b.clone();
    lu.solve_in_place(&mut x).unwrap();

    let elmer = read_values(&data_dir("poisson_pfem").join("solution.dat"));
    let scale = elmer.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (slot, &v) in elmer.iter().enumerate() {
        let ours = x[c.sys.perm.perm[slot] as usize];
        assert!(
            (ours - v).abs() <= 1e-11 * scale,
            "slot {slot}: {ours:e} vs {v:e}"
        );
    }
    // ComputeNorm: sqrt(sum x^2 / n) over all DOFs.
    let norm = (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt();
    assert!((norm - 1.96743752e-2).abs() < 1e-10, "norm {norm:e}");
    // Reference Norm = 0.19674375E-01 with Elmer's default tolerance 1e-5.
    assert!((norm - 0.19674375e-1).abs() / 0.19674375e-1 < 1e-5);
}
