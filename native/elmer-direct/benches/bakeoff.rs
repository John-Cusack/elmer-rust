//! Direct-solver bake-off on `Linear System Save` dumps.
//!
//! `cargo bench -p elmer-direct --bench bakeoff -- [DIR ...]` where each DIR
//! holds `linsys_a.dat` / `linsys_b.dat` (Elmer's save format). With no
//! arguments the fixtures under `tests/fixtures` are used. Prints a timing
//! table (best of `REPS` runs); there is no pass/fail gate.
//!
//! The dump reader is shared with the tests; it moves to `elmer-sparse` once
//! that crate merges (swarm plan L04: delete the duplicate afterwards).
#[path = "../tests/common/mod.rs"]
mod common;

use std::path::PathBuf;
use std::time::Instant;

use common::{is_symmetric, read_dump, rel_residual, Dump};
use elmer_core::DirectSolver;
use elmer_direct::{set_single_thread, BandLu, FaerCholesky, FaerLu};

const REPS: usize = 5;

fn time_solver(name: &str, s: &mut dyn DirectSolver<f64>, d: &Dump) {
    let (mut best_factor, mut best_solve) = (f64::MAX, f64::MAX);
    let mut x = d.b.clone();
    for _ in 0..REPS {
        let t = Instant::now();
        if let Err(e) = s.factor(&d.a) {
            println!("| {name} | factor failed: {e} | | |");
            return;
        }
        best_factor = best_factor.min(t.elapsed().as_secs_f64());
        x.copy_from_slice(&d.b);
        let t = Instant::now();
        s.solve_in_place(&mut x).unwrap();
        best_solve = best_solve.min(t.elapsed().as_secs_f64());
    }
    println!(
        "| {name} | {:.3} | {:.3} | {:.1e} |",
        best_factor * 1e3,
        best_solve * 1e3,
        rel_residual(&d.a, &x, &d.b)
    );
}

fn main() {
    set_single_thread();
    let mut dirs: Vec<PathBuf> = std::env::args().skip(1).filter(|a| !a.starts_with('-')).map(PathBuf::from).collect();
    if dirs.is_empty() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
        for n in ["bodydir", "bodyload", "pfem"] {
            dirs.push(root.join(n));
        }
    }
    for dir in dirs {
        let d = read_dump(&dir);
        println!("\n{}: n = {}, nnz = {}\n", dir.display(), d.a.nrows, d.a.values.len());
        println!("| solver | factor ms | solve ms | rel. residual |\n|---|---:|---:|---:|");
        time_solver("FaerLu", &mut FaerLu::new(), &d);
        time_solver("BandLu", &mut BandLu::new(), &d);
        if is_symmetric(&d.a) {
            time_solver("FaerCholesky", &mut FaerCholesky::new(), &d);
        }
    }
}
