//! Self-check run under `mpirun -np 2` by `tests/np2.rs` (also runs at np=1).
//!
//! Usage: `elmer-mpi-check <mesh_dir> <global_nodes> <global_bulk> <global_boundary>`
//! where `<mesh_dir>/partitioning.<np>` holds an ElmerGrid partition. Exits
//! non-zero on any mismatch and prints `rank R/N OK` on success.

use std::path::Path;
use std::process::ExitCode;

use elmer_core::Communicator;
use elmer_mpi::{dot_owned, HaloPlan, MpiRuntime, PartitionMesh};

fn check(ok: bool, what: &str) {
    if !ok {
        eprintln!("FAILED: {what}");
        std::process::exit(1);
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let rt = MpiRuntime::init().expect("MPI init");
    let comm = rt.world();
    let (rank, size) = (comm.rank(), comm.size());

    // Collectives.
    let mut s = [rank as f64 + 1.0, 10.0 * (rank as f64 + 1.0)];
    comm.sum_f64(&mut s);
    let tri = (size * (size + 1) / 2) as f64;
    check(s == [tri, 10.0 * tri], "sum_f64");
    let mut m = [rank as f64, -(rank as f64)];
    comm.max_f64(&mut m);
    check(m == [(size - 1) as f64, 0.0], "max_f64");
    comm.barrier();

    if size > 1 {
        let dir = Path::new(&args[1]);
        let want: Vec<f64> = args[2..5].iter().map(|a| a.parse().expect("count")).collect();
        let mesh = PartitionMesh::read(dir, rank, size).expect("read partition");
        let plan = HaloPlan::from_mesh(&mesh);

        // Owned entities add up to the global mesh.
        let mut owned = [
            mesh.n_owned_nodes() as f64,
            mesh.n_owned_bulk() as f64,
            mesh.header().n_boundary as f64,
        ];
        comm.sum_f64(&mut owned);
        check(owned == want[..], "owned node/element/boundary counts");

        // halo_sum: every shared node ends up with its number of sharers.
        let mut v = vec![1.0; mesh.header().n_nodes];
        comm.halo_sum(&plan, &mut v, 1).expect("halo_sum");
        let mut sharers = vec![1.0; v.len()];
        for sh in mesh.shared() {
            sharers[sh.local as usize] = sh.parts.len() as f64;
        }
        check(v == sharers, "halo_sum");

        // halo_update with 2 dofs per node: owners' values overwrite ghosts.
        let gids = mesh.node_gids();
        let owners = mesh.owned_nodes();
        let mut w: Vec<f64> =
            (0..gids.len())
                .flat_map(|i| {
                    if owners[i] {
                        [gids[i] as f64, -(gids[i] as f64)]
                    } else {
                        [-1.0, -1.0]
                    }
                })
                .collect();
        comm.halo_update(&plan, &mut w, 2).expect("halo_update");
        let exact = (0..gids.len())
            .all(|i| w[2 * i] == gids[i] as f64 && w[2 * i + 1] == -(gids[i] as f64));
        check(exact, "halo_update");

        // dot_owned counts each shared node once: sum over global ids of gid^2.
        let x: Vec<f64> = gids.iter().map(|&g| g as f64).collect();
        let d = dot_owned(&comm, &x, &x, &owners);
        let n = want[0];
        check(d == n * (n + 1.0) * (2.0 * n + 1.0) / 6.0, "dot_owned");
    }

    println!("rank {rank}/{size} OK");
    ExitCode::SUCCESS
}
