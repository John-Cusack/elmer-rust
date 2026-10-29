//! np=1 and np=2 runs of `elmer-mpi-check`. np=2 needs `mpirun` and is skipped
//! with a message when it is not installed. Always `HWLOC_COMPONENTS=-gl`,
//! never more than 2 ranks.

use std::path::PathBuf;
use std::process::Command;

const CHECK: &str = env!("CARGO_BIN_EXE_elmer-mpi-check");

fn data(sub: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(sub)
}

fn mpirun(np: usize, mesh: &str) -> Option<std::process::Output> {
    let mut cmd = Command::new("mpirun");
    cmd.args(["--oversubscribe", "-np", &np.to_string(), CHECK])
        .arg(data(mesh))
        .args(["121", "100", "40"]);
    // Single-threaded: the swarm host is shared.
    cmd.env("HWLOC_COMPONENTS", "-gl").env("OMP_NUM_THREADS", "1");
    match cmd.output() {
        Ok(o) => Some(o),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("mpirun not found; skipping np={np}");
            None
        }
        Err(e) => panic!("mpirun: {e}"),
    }
}

fn run_np2(mesh: &str) {
    let Some(out) = mpirun(2, mesh) else { return };
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success(),
        "mpirun failed:\n{stdout}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(stdout.contains("rank 0/2 OK") && stdout.contains("rank 1/2 OK"), "{stdout}");
}

#[test]
fn np2_sum_and_halo_exchange_plain_partition() {
    run_np2("sq2");
}

#[test]
fn np2_sum_and_halo_exchange_with_ghost_elements() {
    run_np2("sqh");
}

#[test]
fn np1_collectives() {
    let Some(out) = mpirun(1, "sq2") else { return };
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stdout).contains("rank 0/1 OK"));
}
