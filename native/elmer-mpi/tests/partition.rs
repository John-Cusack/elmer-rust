//! Serial tests of the partition-directory reader against ElmerGrid output for
//! `fem/tests/PoissonPFEM/square` (121 nodes, 100 quads, 40 boundary lines),
//! partitioned with `ElmerGrid 2 2 square -metis 2 4 [-halo]`.

use std::path::PathBuf;

use elmer_core::{ElemId, MeshView};
use elmer_mpi::{HaloPlan, PartitionMesh};

fn data(sub: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(sub)
}

fn read_both(sub: &str) -> Vec<PartitionMesh> {
    (0..2).map(|r| PartitionMesh::read(&data(sub), r, 2).unwrap()).collect()
}

#[test]
fn counts_match_the_global_mesh() {
    // mesh.header of the unpartitioned mesh: 121 nodes, 100 elements, 40 boundary.
    let g: Vec<usize> = std::fs::read_to_string(data("mesh.header"))
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect();
    for sub in ["sq2", "sqh"] {
        let parts = read_both(sub);
        assert_eq!(parts.iter().map(|p| p.n_owned_nodes()).sum::<usize>(), g[0], "{sub} nodes");
        assert_eq!(parts.iter().map(|p| p.n_owned_bulk()).sum::<usize>(), g[1], "{sub} elements");
        assert_eq!(parts.iter().map(|p| p.n_boundary()).sum::<usize>(), g[2], "{sub} boundary");
        for p in &parts {
            assert_eq!(p.n_nodes(), p.header().n_nodes);
            assert_eq!(p.dim(), 2);
        }
    }
}

#[test]
fn plain_partition_shape() {
    let parts = read_both("sq2");
    assert_eq!((parts[0].n_nodes(), parts[0].n_bulk(), parts[0].n_boundary()), (66, 50, 20));
    assert_eq!((parts[1].n_nodes(), parts[1].n_bulk(), parts[1].n_boundary()), (66, 50, 20));
    assert_eq!(parts[0].shared().len(), 11);
    // Node 56 is shared and owned by partition 2 (rank 1); 61 by rank 0.
    let owner_of = |p: &PartitionMesh, gid: u64| {
        let i = p.node_gids().iter().position(|&g| g == gid).unwrap();
        p.node_owners()[i]
    };
    assert_eq!(owner_of(&parts[0], 56), 1);
    assert_eq!(owner_of(&parts[0], 61), 0);
    assert_eq!(owner_of(&parts[1], 56), 1);
    assert!(parts.iter().all(|p| p.bulk_elements().iter().all(|e| e.owner == p.rank())));
}

#[test]
fn halo_partition_has_ghost_elements_owned_elsewhere() {
    let parts = read_both("sqh");
    for p in &parts {
        let ghosts = p.bulk_elements().iter().filter(|e| e.owner != p.rank()).count();
        assert_eq!(ghosts, p.n_bulk() - p.n_owned_bulk());
        assert!(ghosts > 0, "rank {} has no ghost elements", p.rank());
        assert!(p
            .bulk_elements()
            .iter()
            .filter(|e| e.owner != p.rank())
            .all(|e| e.owner == 1 - p.rank()));
    }
}

#[test]
fn coordinates_and_boundary_parents() {
    let p = &read_both("sq2")[1];
    // Local node 0 of rank 1 is global node 56 at (0, 0.5).
    assert_eq!(p.node_gids()[0], 56);
    assert_eq!(p.coord(0), [0.0, 0.5, 0.0]);
    for b in 0..p.n_boundary() {
        let e = p.element(ElemId((p.n_bulk() + b) as u32));
        assert_eq!(e.etype.0, 202);
        let parent = e.parents[0].expect("parent element is in the partition");
        assert!((parent.0 as usize) < p.n_bulk());
    }
}

#[test]
fn halo_plan_is_symmetric() {
    let parts = read_both("sq2");
    let plans: Vec<HaloPlan> = parts.iter().map(HaloPlan::from_mesh).collect();
    for (a, b) in [(0, 1), (1, 0)] {
        let na = &plans[a].neighbours()[0];
        let nb = &plans[b].neighbours()[0];
        assert_eq!(na.rank, b);
        let gids = |p: &PartitionMesh, l: &[u32]| {
            l.iter().map(|&i| p.node_gids()[i as usize]).collect::<Vec<_>>()
        };
        assert_eq!(gids(&parts[a], &na.shared), gids(&parts[b], &nb.shared));
        assert_eq!(gids(&parts[a], &na.send_owned), gids(&parts[b], &nb.recv_ghost));
        assert_eq!(na.shared.len(), na.send_owned.len() + na.recv_ghost.len());
    }
}

#[test]
fn missing_file_and_bad_counts_are_errors() {
    assert!(PartitionMesh::read(&data("nope"), 0, 2).is_err());
    // rank 2 of a 2-part directory has no files.
    assert!(PartitionMesh::read(&data("sq2"), 2, 2).is_err());
}
