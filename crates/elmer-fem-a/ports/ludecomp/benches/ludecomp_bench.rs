//! A/B bench for the `ludecomp` port (lane-fem-a, unit `c:fem/src/view3d/LUDecomp.c`).
//!
//! Two measured workloads, both on the same deterministic corpus:
//!
//! * `ludecomp` — `port_ludecomp`: dense LU with partial pivoting, O(n^3).
//! * `mtrinv`  — `port_mtrinv`: LU + full inverse, O(n^3).
//!
//! The corpus is a fixed-seed, diagonally dominant matrix, so both arms see
//! byte-identical input. The digest (FNV-1a over the result bytes) is the
//! equivalence proof: the harness refuses to report a speed when the two
//! arms disagree.
//!
//! The corpus copy is inside the measured region (the port mutates its
//! argument in place, so each iteration needs a fresh copy). The copy-only
//! cost is measured separately and reported as the harness overhead row —
//! it is never subtracted from the workload number.

use elmer_bench_harness as harness;

/// Deterministic corpus: fixed-seed xorshift64* values, diagonally dominant
/// so LU never hits the singular path.
fn corpus(n: usize) -> Vec<f64> {
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut m = vec![0.0f64; n * n];
    for row in 0..n {
        for col in 0..n {
            let v = (next() >> 11) as f64 / (1u64 << 53) as f64; // [0,1)
            m[row * n + col] = (v - 0.5) * 2.0; // [-1,1)
        }
        m[row * n + row] += n as f64; // diagonally dominant
    }
    m
}

/// FNV-1a over the raw bytes of the result matrix.
fn digest(buf: &[f64]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for v in buf {
        for b in v.to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    format!("{h:016x}")
}

fn child(n: usize, iters: usize, workload: &str) -> Result<String, String> {
    let base = corpus(n);
    let mut scratch = vec![0.0f64; n * n];
    let mut pivot = vec![0i32; n];
    let nn = n as std::os::raw::c_int;

    // Warm-up (not measured): first touch of the code path and the pages.
    scratch.copy_from_slice(&base);
    // SAFETY: `scratch` holds n*n doubles, `pivot` holds n ints.
    unsafe { ludecomp::port_ludecomp(scratch.as_mut_ptr(), nn, pivot.as_mut_ptr()) };

    // Measured region: corpus copy + the port call.
    let t = harness::Timer::start();
    for _ in 0..iters {
        scratch.copy_from_slice(&base);
        // SAFETY: as above; the port mutates the scratch matrix in place.
        unsafe {
            match workload {
                "mtrinv" => ludecomp::port_mtrinv(scratch.as_mut_ptr(), nn),
                _ => ludecomp::port_ludecomp(scratch.as_mut_ptr(), nn, pivot.as_mut_ptr()),
            }
        };
    }
    let wall_ns = t.elapsed_ns();

    // Digest of the RESULT (taken before the overhead loop below overwrites
    // the scratch buffer with the corpus again).
    let d = digest(&scratch);

    // Harness overhead: the same number of copies, no port call.
    let t = harness::Timer::start();
    for _ in 0..iters {
        scratch.copy_from_slice(&base);
    }
    let harness_ns = t.elapsed_ns();

    let rss = harness::self_peak_rss_kb();
    Ok(harness::child_json(wall_ns, harness_ns, rss, &d))
}

fn main() {
    let args = match harness::parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("ludecomp_bench: {e}");
            std::process::exit(2);
        }
    };
    let r = match args.mode.as_str() {
        "child" => child(args.n, args.iters, &args.workload),
        "parent" => harness::run_parent(&args),
        other => Err(format!("unknown mode: {other}")),
    };
    match r {
        Ok(json) => {
            if args.mode == "child" {
                println!("{json}");
            }
        }
        Err(e) => {
            eprintln!("ludecomp_bench: {e}");
            std::process::exit(1);
        }
    }
}
