//! Perm, bandwidth decision, CSR pattern and Dirichlet enforcement against
//! pristine ElmerSolver dumps (Rev 9f6af2f85; provenance in each file's
//! header and in `tests/data/README.md`).

mod support;

use elmer_assembly::{
    apply_dirichlet, create_system, DirichletBc, DirichletMode, DirichletScaling, DirichletSet,
    PermOptions, System,
};
use elmer_core::{CsrMatrix, ElemId, MeshView, NO_DOF};
use support::{data_dir, read_matrix, read_perm, read_rhs, TestMesh};

/// One Dirichlet BC on one component, as in the case's `.sif`.
struct Bc {
    tags: &'static [i32],
    component: usize,
    value: fn([f64; 3]) -> f64,
}

struct Case {
    name: &'static str,
    dim: usize,
    dofs: usize,
    /// Bodies whose equation has the solver active.
    bodies: &'static [i32],
    /// Elmer's printed initial and optimised half bandwidths.
    bandwidth: (usize, usize),
    mode: DirichletMode,
    bcs: &'static [Bc],
}

const DEFAULT_MODE: DirichletMode = DirichletMode {
    symmetric: false,
    scaling: DirichletScaling::Diagonal,
};

const CASES: &[Case] = &[
    Case {
        name: "minimal_example",
        dim: 2,
        dofs: 1,
        bodies: &[1],
        bandwidth: (3, 3),
        mode: DEFAULT_MODE,
        bcs: &[Bc {
            tags: &[2],
            component: 0,
            value: |_| 0.0,
        }],
    },
    Case {
        name: "thermal_pipe_2d",
        dim: 2,
        dofs: 1,
        bodies: &[1, 2, 3],
        bandwidth: (889, 60),
        mode: DEFAULT_MODE,
        bcs: &[
            Bc {
                tags: &[9, 10, 11, 12],
                component: 0,
                value: |_| 10.0,
            },
            Bc {
                tags: &[1, 2, 3, 4],
                component: 0,
                value: |_| 0.0,
            },
        ],
    },
    Case {
        name: "hex_elasticity",
        dim: 3,
        dofs: 3,
        bodies: &[1, 2],
        bandwidth: (26, 25),
        mode: DEFAULT_MODE,
        bcs: &[
            Bc {
                tags: &[1],
                component: 0,
                value: |_| 0.0,
            },
            Bc {
                tags: &[1],
                component: 1,
                value: |_| 0.0,
            },
            Bc {
                tags: &[1],
                component: 2,
                value: |_| 0.0,
            },
            Bc {
                tags: &[2],
                component: 0,
                value: |_| 0.1,
            },
            Bc {
                tags: &[2],
                component: 2,
                value: |x| 0.01 * x[1],
            },
        ],
    },
    Case {
        name: "hex_heat_partial",
        dim: 3,
        dofs: 1,
        bodies: &[1],
        bandwidth: (29, 25),
        mode: DirichletMode {
            symmetric: true,
            scaling: DirichletScaling::Diagonal,
        },
        bcs: &[
            Bc {
                tags: &[1],
                component: 0,
                value: |x| 1.0 + x[1],
            },
            Bc {
                tags: &[3],
                component: 0,
                value: |_| 0.5,
            },
        ],
    },
];

fn system(case: &Case) -> (TestMesh, System) {
    let mesh = TestMesh::read(&data_dir(case.name), case.dim);
    // CheckElementEquation: bulk elements by body; boundary elements have
    // body id 0 (no `Body Id` in these BCs) and are inactive.
    let active = |e: ElemId| !mesh.is_boundary(e) && case.bodies.contains(&mesh.element(e).tag);
    let sys = create_system(&mesh, &active, PermOptions::default(), case.dofs);
    (mesh, sys)
}

fn matrix_from_dump(pattern: &CsrMatrix<f64>, dump: &[(usize, usize, f64)]) -> CsrMatrix<f64> {
    let mut a = pattern.clone();
    assert_eq!(dump.len(), a.nnz());
    for (k, &(_, _, v)) in dump.iter().enumerate() {
        a.values[k] = v;
    }
    a
}

#[test]
fn perm_matches_elmer() {
    for case in CASES {
        let (_, sys) = system(case);
        let expected: Vec<u32> = read_perm(&data_dir(case.name).join("perm.dat"))
            .into_iter()
            .map(|p| if p == 0 { NO_DOF } else { p - 1 })
            .collect();
        assert_eq!(sys.perm.perm, expected, "{}", case.name);
    }
}

#[test]
fn bandwidth_decision_matches_elmer() {
    for case in CASES {
        let (_, sys) = system(case);
        let b = sys.bandwidth;
        assert_eq!(
            (b.initial, b.optimized),
            (case.bandwidth.0, Some(case.bandwidth.1)),
            "{}",
            case.name
        );
        assert!(b.accepted, "{}", case.name);
    }
}

#[test]
fn csr_pattern_matches_elmer() {
    for case in CASES {
        let (_, sys) = system(case);
        let dump = read_matrix(&data_dir(case.name).join("pre_a.dat"));
        let a = &sys.matrix;
        assert_eq!(a.nrows, sys.perm.n_dofs * case.dofs, "{}", case.name);
        let mut k = 0;
        for r in 0..a.nrows {
            for &c in &a.col_idx[a.row_ptr[r]..a.row_ptr[r + 1]] {
                assert_eq!((dump[k].0, dump[k].1), (r, c as usize), "{}", case.name);
                k += 1;
            }
        }
        assert_eq!(k, dump.len(), "{}", case.name);
        let diag = a.diag.as_ref().unwrap();
        for (r, &d) in diag.iter().enumerate() {
            assert_eq!(a.col_idx[d] as usize, r);
        }
    }
}

#[test]
fn dirichlet_matches_elmer_bitwise() {
    for case in CASES {
        let (mesh, sys) = system(case);
        let dir = data_dir(case.name);
        let pre = read_matrix(&dir.join("pre_a.dat"));
        let mut a = matrix_from_dump(&sys.matrix, &pre);
        let mut b: Vec<f64> = read_rhs(&dir.join("pre_b.dat"))
            .into_iter()
            .map(|v| v.unwrap_or(0.0))
            .collect();

        let mut set = DirichletSet::new(a.nrows);
        for comp in 0..case.dofs {
            let elements: Vec<Box<dyn Fn(ElemId) -> bool + '_>> = case
                .bcs
                .iter()
                .map(|bc| {
                    let mesh = &mesh;
                    Box::new(move |e: ElemId| {
                        mesh.is_boundary(e) && bc.tags.contains(&mesh.element(e).tag)
                    }) as Box<dyn Fn(ElemId) -> bool>
                })
                .collect();
            let values: Vec<Box<dyn Fn(ElemId, u32) -> f64 + '_>> = case
                .bcs
                .iter()
                .map(|bc| {
                    let mesh = &mesh;
                    Box::new(move |_: ElemId, n: u32| (bc.value)(mesh.coord(n)))
                        as Box<dyn Fn(ElemId, u32) -> f64>
                })
                .collect();
            let bcs: Vec<DirichletBc<'_>> = case
                .bcs
                .iter()
                .enumerate()
                .filter(|(_, bc)| bc.component == comp)
                .map(|(i, _)| DirichletBc {
                    elements: &*elements[i],
                    value: &*values[i],
                    condition: None,
                })
                .collect();
            set.set_boundaries(&mesh, &sys.perm, comp, case.dofs, &bcs, false);
        }
        apply_dirichlet(&mut a, &mut b, &set.to_list(), case.mode).unwrap();

        let mut expected = pre.clone();
        let changed = read_matrix(&dir.join("post_a_changed.dat"));
        let mut c = changed.iter().peekable();
        for e in &mut expected {
            if let Some(&&(r, col, v)) = c.peek() {
                if (r, col) == (e.0, e.1) {
                    e.2 = v;
                    c.next();
                }
            }
        }
        assert!(
            c.next().is_none(),
            "{}: unmatched changed entries",
            case.name
        );
        for (k, &(r, col, v)) in expected.iter().enumerate() {
            assert_eq!(
                a.values[k].to_bits(),
                v.to_bits(),
                "{}: A[{r},{col}] = {:e}, Elmer {:e}",
                case.name,
                a.values[k],
                v
            );
        }
        for (i, want) in read_rhs(&dir.join("post_b.dat")).into_iter().enumerate() {
            match want {
                None => assert!(b[i].abs() <= f64::MIN_POSITIVE, "{}: b[{i}]", case.name),
                Some(v) => assert_eq!(b[i].to_bits(), v.to_bits(), "{}: b[{i}]", case.name),
            }
        }
    }
}
