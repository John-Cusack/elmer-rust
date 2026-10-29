//! End-to-end on MinimalExample: Perm and pattern, an element loop with a
//! test-local P1 heat kernel (bulk conduction, boundary heat flux),
//! Dirichlet with Elmer defaults, a dense solve, and the solution norm
//! against the case's reference (`Solver 1 :: Reference Norm`).
//!
//! The kernel is exact P1 algebra, not Elmer's quadrature, so the matrix
//! agrees with Elmer's dump to rounding rather than bitwise; the bitwise
//! checks of the assembly steps are in `elmer_dumps.rs`.

mod support;

use elmer_assembly::{
    apply_dirichlet, assemble_elements, boundary_elements, bulk_elements, create_system,
    CsrAssembler, DirichletBc, DirichletMode, DirichletSet, PermOptions,
};
use elmer_core::{ElemId, ElementRef, LocalSystem, MeshView};
use support::{data_dir, read_matrix, TestMesh};

/// `case.sif`: Reference Norm = 2.88675135E-01.
const REFERENCE_NORM: f64 = 2.88675135e-01;

fn p1_triangle(mesh: &TestMesh, el: ElementRef<'_>, local: &mut LocalSystem) {
    let x: Vec<[f64; 3]> = el.nodes.iter().map(|&n| mesh.coord(n)).collect();
    let det = (x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) - (x[2][0] - x[0][0]) * (x[1][1] - x[0][1]);
    let area = 0.5 * det.abs();
    // Gradients of the barycentric basis functions times det.
    let g = [
        [x[1][1] - x[2][1], x[2][0] - x[1][0]],
        [x[2][1] - x[0][1], x[0][0] - x[2][0]],
        [x[0][1] - x[1][1], x[1][0] - x[0][0]],
    ];
    let conductivity = 1.0;
    for i in 0..3 {
        for j in 0..3 {
            local.stiff[i * 3 + j] =
                conductivity * area * (g[i][0] * g[j][0] + g[i][1] * g[j][1]) / (det * det);
        }
    }
}

fn flux_line(mesh: &TestMesh, el: ElementRef<'_>, local: &mut LocalSystem, flux: f64) {
    let (a, b) = (mesh.coord(el.nodes[0]), mesh.coord(el.nodes[1]));
    let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
    local.force[0] = flux * len / 2.0;
    local.force[1] = flux * len / 2.0;
}

/// Gaussian elimination with partial pivoting on a small dense copy.
fn dense_solve(a: &elmer_core::CsrMatrix<f64>, b: &[f64]) -> Vec<f64> {
    let n = a.nrows;
    let mut m = vec![vec![0.0; n + 1]; n];
    for (r, row) in m.iter_mut().enumerate() {
        for k in a.row_ptr[r]..a.row_ptr[r + 1] {
            row[a.col_idx[k] as usize] = a.values[k];
        }
        row[n] = b[r];
    }
    for c in 0..n {
        let p = (c..n)
            .max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs()))
            .unwrap();
        m.swap(c, p);
        for r in c + 1..n {
            let f = m[r][c] / m[c][c];
            for k in c..=n {
                m[r][k] -= f * m[c][k];
            }
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = (r + 1..n).map(|k| m[r][k] * x[k]).sum();
        x[r] = (m[r][n] - s) / m[r][r];
    }
    x
}

#[test]
fn minimal_example_solves_to_reference_norm() {
    let dir = data_dir("minimal_example");
    let mesh = TestMesh::read(&dir, 2);
    let is_boundary = |e: ElemId| e.0 as usize >= mesh.n_bulk();
    let active = |e: ElemId| !is_boundary(e) && mesh.element(e).tag == 1;
    let mut sys = create_system(&mesh, &active, PermOptions::default(), 1);
    let n = sys.matrix.nrows;
    let mut b = vec![0.0; n];
    {
        let mut asm = CsrAssembler::new(&mut sys.matrix, &mut b);
        let bulk = bulk_elements(&mesh).filter(|&e| active(e));
        assemble_elements::<()>(
            &mesh,
            &sys.perm,
            1,
            bulk.chain(boundary_elements(&mesh)),
            &mut |e, el, local| {
                if !is_boundary(e) {
                    p1_triangle(&mesh, el, local);
                    Ok(true)
                } else if el.tag == 3 {
                    // Boundary Condition 2: Heat Flux = 1.
                    flux_line(&mesh, el, local, 1.0);
                    Ok(true)
                } else {
                    Ok(false)
                }
            },
            &mut asm,
        )
        .unwrap();
    }

    // Before Dirichlet the system matches Elmer's dump to rounding.
    let dump = read_matrix(&dir.join("pre_a.dat"));
    for (k, &(_, _, v)) in dump.iter().enumerate() {
        assert!((sys.matrix.values[k] - v).abs() <= 1e-15, "entry {k}");
    }

    // Boundary Condition 1: Temperature = 0 on boundary 2.
    let mut set = DirichletSet::new(n);
    let on = |e: ElemId| is_boundary(e) && mesh.element(e).tag == 2;
    let zero = |_: ElemId, _: u32| 0.0;
    let bcs = [DirichletBc {
        elements: &on,
        value: &zero,
        condition: None,
    }];
    set.set_boundaries(&mesh, &sys.perm, 0, 1, &bcs, false);
    apply_dirichlet(
        &mut sys.matrix,
        &mut b,
        &set.to_list(),
        DirichletMode::default(),
    )
    .unwrap();

    let x = dense_solve(&sys.matrix, &b);
    // Elmer's ComputeNorm default: RMS over the solver's DOFs.
    let norm = (x.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
    let rel = (norm - REFERENCE_NORM).abs() / REFERENCE_NORM;
    assert!(
        rel < 1e-8,
        "norm {norm}, reference {REFERENCE_NORM}, rel {rel:e}"
    );
}
