//! np=1 against np=2 on in-process ranks ([`ThreadComm`]): the same Poisson
//! problem read from `partitioning.1` and from `partitioning.2` (plain and
//! with halo elements) must give the same matrix action, the same HUTI
//! iteration counts and the same solution.

mod support;

use std::collections::HashMap;

use elmer_core::{Communicator, MeshView, ParamMap, SolveStatus, Value};
use elmer_dist::{
    dist_solve, owners_by_lowest_rank, DistCsr, DistError, DistPrecond, DistSolveOptions,
    DofLayout, ThreadComm,
};
use elmer_krylov::{HutiConfig, HutiMethod};
use elmer_mpi::PartitionMesh;
use support::poisson::{assemble, data};

struct Run {
    iterations: usize,
    status: SolveStatus,
    /// Every local value (owned and ghost) by global dof id.
    x: Vec<(u64, f64)>,
}

fn opts(method: HutiMethod, precond: DistPrecond) -> DistSolveOptions {
    let mut huti = HutiConfig::new(method);
    huti.tol = 1e-10;
    huti.max_iter = 500;
    DistSolveOptions {
        huti,
        precond,
        scaling: true,
    }
}

fn run(mesh: &str, np: usize, o: &DistSolveOptions) -> Vec<Run> {
    ThreadComm::run(np, |c| {
        let p = assemble(c, &data(mesh));
        let mut x = vec![0.0; p.layout.n_local()];
        let s = dist_solve(c, &p.layout, &p.a, &p.b, &mut x, o, None).expect("solve");
        Run {
            iterations: s.iterations,
            status: s.status,
            x: p.layout.gids().iter().copied().zip(x).collect(),
        }
    })
}

/// Solution by global id; asserts that every copy of a dof agrees.
fn by_gid(runs: &[Run]) -> HashMap<u64, f64> {
    let mut m = HashMap::new();
    for r in runs {
        for &(g, v) in &r.x {
            if let Some(old) = m.insert(g, v) {
                assert_eq!(old, v, "copies of dof {g} differ");
            }
        }
    }
    m
}

fn compare(mesh2: &str, o: &DistSolveOptions, same_iterations: bool) {
    let one = run("sq", 1, o);
    let two = run(mesh2, 2, o);
    assert_eq!(one[0].status, SolveStatus::Converged);
    for r in &two {
        assert_eq!(r.status, SolveStatus::Converged);
        assert_eq!(
            r.iterations, two[0].iterations,
            "ranks disagree on the iteration count"
        );
    }
    if same_iterations {
        assert_eq!(
            one[0].iterations, two[0].iterations,
            "np=1 vs np=2 iterations"
        );
    }
    let (s1, s2) = (by_gid(&one), by_gid(&two));
    assert_eq!(s1.len(), 81);
    assert_eq!(s2.len(), 81);
    let max = s1.values().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(max > 0.07 && max < 0.08, "peak {max}"); // about 0.0737 for the unit square
    for (g, v) in &s1 {
        let d = (v - s2[g]).abs();
        assert!(d <= 1e-8 * max, "dof {g}: np=1 {v}, np=2 {}", s2[g]);
    }
}

#[test]
fn cg_no_preconditioner_np1_equals_np2() {
    compare("sq", &opts(HutiMethod::Cg, DistPrecond::None), true);
}

#[test]
fn cg_diagonal_np1_equals_np2_with_halo_elements() {
    compare("sqh", &opts(HutiMethod::Cg, DistPrecond::Diagonal), true);
}

#[test]
fn bicgstab_np1_equals_np2() {
    compare("sq", &opts(HutiMethod::BiCGStab, DistPrecond::None), true);
    compare(
        "sqh",
        &opts(HutiMethod::BiCGStab, DistPrecond::Diagonal),
        true,
    );
}

#[test]
fn block_ilu0_converges_to_the_same_solution() {
    // Block-local ILU differs from serial ILU, so only the solution must agree.
    compare(
        "sq",
        &opts(HutiMethod::BiCGStab, DistPrecond::Ilu(0)),
        false,
    );
}

#[test]
fn unscaled_cg_np1_equals_np2() {
    let mut o = opts(HutiMethod::Cg, DistPrecond::None);
    o.scaling = false;
    compare("sq", &o, true);
}

/// The distributed matvec equals the serial product with the summed local
/// matrices (shared Dirichlet rows sum to the number of holders, as in Elmer).
#[test]
fn matvec_matches_summed_operator() {
    let f = |g: u64| (g as f64 * 0.37).sin() + 0.1 * g as f64;
    for mesh in ["sq", "sqh"] {
        type Part = (Vec<(u64, u64, f64)>, Vec<(u64, f64)>);
        let parts: Vec<Part> = ThreadComm::run(2, |c| {
            let p = assemble(c, &data(mesh));
            let g = p.layout.gids();
            let mut entries = Vec::new();
            for i in 0..p.a.nrows {
                for k in p.a.row_ptr[i]..p.a.row_ptr[i + 1] {
                    entries.push((g[i], g[p.a.col_idx[k] as usize], p.a.values[k]));
                }
            }
            let m = DistCsr::from_local(c, &p.layout, &p.a).expect("split");
            let owned_g: Vec<u64> = p.layout.owned().iter().map(|&l| g[l as usize]).collect();
            let x: Vec<f64> = owned_g.iter().map(|&g| f(g)).collect();
            let mut y = vec![0.0; x.len()];
            m.matvec(c, &x, &mut y);
            let (to, from) = m.neighbours();
            assert_eq!((to, from), (vec![1 - c.rank()], vec![1 - c.rank()]));
            (entries, owned_g.into_iter().zip(y).collect())
        });
        let mut want: HashMap<u64, f64> = HashMap::new();
        for (entries, _) in &parts {
            for &(i, j, v) in entries {
                *want.entry(i).or_default() += v * f(j);
            }
        }
        let got: Vec<(u64, f64)> = parts.into_iter().flat_map(|p| p.1).collect();
        assert_eq!(got.len(), 81, "{mesh}: every dof owned exactly once");
        for (g, v) in got {
            let w = want[&g];
            assert!(
                (v - w).abs() <= 1e-13 * (1.0 + w.abs()),
                "{mesh} dof {g}: {v} vs {w}"
            );
        }
    }
}

#[test]
fn layout_sums_and_updates_shared_dofs() {
    ThreadComm::run(2, |c| {
        let mesh = PartitionMesh::read(&data("sq"), c.rank(), 2).expect("read");
        let l = DofLayout::nodal(c, &mesh, 2).expect("layout");
        assert_eq!(l.n_local(), 2 * mesh.n_nodes());

        // Ones summed over holders: 2 on shared dofs, 1 elsewhere.
        let mut v = vec![1.0; l.n_local()];
        l.sum_shared(c, &mut v);
        let shared: Vec<u32> = mesh.shared().iter().map(|s| s.local).collect();
        for (i, &x) in v.iter().enumerate() {
            let want = if shared.contains(&(i as u32 / 2)) {
                2.0
            } else {
                1.0
            };
            assert_eq!(x, want, "dof {i}");
        }

        // Owner values overwrite ghosts.
        let owned: Vec<f64> = l
            .owned()
            .iter()
            .map(|&i| l.gids()[i as usize] as f64)
            .collect();
        let mut local = vec![-1.0; l.n_local()];
        l.owned_to_local(c, &owned, &mut local);
        for (i, &x) in local.iter().enumerate() {
            assert_eq!(x, l.gids()[i] as f64);
        }

        // Owned totals add up to the global count of dofs.
        let mut tot = [l
            .sum_to_owned(c, &vec![1.0; l.n_local()])
            .iter()
            .sum::<f64>()];
        let mut owned_n = [l.n_owned() as f64];
        c.sum_f64(&mut tot);
        c.sum_f64(&mut owned_n);
        let held = [2.0 * mesh.n_nodes() as f64];
        let mut held_all = held;
        c.sum_f64(&mut held_all);
        assert_eq!(owned_n[0], 162.0);
        assert_eq!(tot[0], held_all[0]);
    });
}

#[test]
fn lowest_rank_owners_are_consistent() {
    ThreadComm::run(2, |c| {
        let mesh = PartitionMesh::read(&data("sq"), c.rank(), 2).expect("read");
        let owners = owners_by_lowest_rank(c, mesh.node_gids());
        let shared: HashMap<u32, &Vec<usize>> =
            mesh.shared().iter().map(|s| (s.local, &s.parts)).collect();
        for (i, &o) in owners.iter().enumerate() {
            let want = shared
                .get(&(i as u32))
                .map_or(c.rank(), |p| *p.iter().min().unwrap());
            assert_eq!(o, want, "node {i}");
        }
        // And they make a valid layout.
        DofLayout::new(c, mesh.node_gids().to_vec(), owners).expect("layout");
    });
}

#[test]
fn a_bad_owner_fails_on_every_rank() {
    let res = ThreadComm::run(2, |c| {
        // Rank 1 claims rank 0 owns dof 99, which rank 0 does not hold.
        let (gid, owner) = if c.rank() == 0 {
            (vec![1, 2], vec![0, 0])
        } else {
            (vec![3, 99], vec![1, 0])
        };
        DofLayout::new(c, gid, owner).err()
    });
    assert!(
        matches!(res[0], Some(DistError::Local { rank: 0, .. })),
        "{:?}",
        res[0]
    );
    assert!(
        matches!(res[1], Some(DistError::Remote { rank: 1 })),
        "{:?}",
        res[1]
    );

    let dup = ThreadComm::run(2, |c| {
        let gid = if c.rank() == 0 { vec![1, 1] } else { vec![2] };
        let n = gid.len();
        DofLayout::new(c, gid, vec![c.rank(); n]).err()
    });
    assert!(matches!(dup[0], Some(DistError::Local { .. })));
    assert!(matches!(dup[1], Some(DistError::Remote { .. })));
}

#[test]
fn zero_rhs_returns_zero_without_iterating() {
    ThreadComm::run(2, |c| {
        let p = assemble(c, &data("sq"));
        let mut x = vec![3.0; p.layout.n_local()];
        let b = vec![0.0; x.len()];
        let o = opts(HutiMethod::Cg, DistPrecond::None);
        let s = dist_solve(c, &p.layout, &p.a, &b, &mut x, &o, None).expect("solve");
        assert_eq!(s.iterations, 0);
        assert!(x.iter().all(|&v| v == 0.0));
    });
}

#[test]
fn options_follow_itersolve_names() {
    let mut p = ParamMap::new();
    p.insert("Linear System Iterative Method", Value::Str("CG".into()));
    p.insert("Linear System Max Iterations", Value::Integer(vec![10]));
    let get = |p: &ParamMap| DistSolveOptions::from_params(p);
    assert_eq!(get(&p).unwrap().precond, DistPrecond::None);
    assert!(get(&p).unwrap().scaling);
    for (name, want) in [
        ("Diagonal", DistPrecond::Diagonal),
        ("ILU", DistPrecond::Ilu(0)),
        ("ILU2", DistPrecond::Ilu(2)),
        ("something", DistPrecond::None),
    ] {
        p.insert("Linear System Preconditioning", Value::Str(name.into()));
        assert_eq!(get(&p).unwrap().precond, want, "{name}");
    }
    p.insert("Linear System Preconditioning", Value::Str("ILU2".into()));
    p.insert("Linear System ILU Order", Value::Integer(vec![1]));
    assert_eq!(get(&p).unwrap().precond, DistPrecond::Ilu(1));
    p.insert("Linear System Preconditioning", Value::Str("vanka".into()));
    assert!(get(&p).is_err());
    p.insert("Linear System Preconditioning", Value::Str("none".into()));
    p.insert("Linear System Scaling", Value::Logical(false));
    assert!(!get(&p).unwrap().scaling);
}

#[test]
fn partitions_cover_the_mesh() {
    // Guards the fixture: 81 nodes, 64 quads split across two ranks.
    let n: Vec<(usize, usize)> = ThreadComm::run(2, |c| {
        let m = PartitionMesh::read(&data("sqh"), c.rank(), 2).expect("read");
        (m.n_owned_nodes(), m.n_owned_bulk())
    });
    assert_eq!(n[0].0 + n[1].0, 81);
    assert_eq!(n[0].1 + n[1].1, 64);
    let m = PartitionMesh::read(&data("sq"), 0, 1).expect("read");
    assert_eq!((m.n_nodes(), m.n_bulk()), (81, 64));
}
