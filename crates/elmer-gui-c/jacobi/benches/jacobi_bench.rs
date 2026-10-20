//! A/B bench for the `jacobi` port (lane-gui-c, unit
//! `c:ElmerGUI/matc/src/jacobi.c` — the generalized Jacobi eigensolver).
//!
//! Measured workload: `matc_jacobi` on a fixed-seed pair of symmetric
//! positive definite matrices (A, B) of dimension n. The solver sweeps until
//! the eigenvalues converge, so the measured region is the whole
//! O(n^3)-per-sweep kernel, not a single step.
//!
//! The digest (FNV-1a over the eigenvectors, eigenvalues and the work
//! vector) is the equivalence proof: the harness refuses to report a speed
//! when the two arms disagree.
//!
//! Link stubs: the port crate imports the sibling matc units for the error
//! path and for the REPL wrapper. `matc_jacobi` itself needs none of them,
//! but the linker does, so the bench provides them. `PrintOut` drops the
//! solver's progress line (it is REPL chatter, not part of the measurement);
//! `error_matc` aborts the way the matc error path does.

use elmer_bench_harness as harness;
use gui_c_jacobi::{c_char, c_int, size_t, VARIABLE};

/// matc error path: the original unwinds with longjmp; the bench aborts.
#[no_mangle]
pub unsafe extern "C" fn error_matc(_fmt: *const c_char, _: ...) -> ! {
    eprintln!("jacobi_bench: error_matc reached");
    std::process::abort();
}

/// matc progress line (dropped: REPL chatter, not measurement).
#[no_mangle]
pub unsafe extern "C" fn PrintOut(_fmt: *const c_char, _: ...) {}

/// Not used by `matc_jacobi`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn var_new(
    _nm: *mut c_char,
    _typ: c_int,
    _nrow: c_int,
    _ncol: c_int,
) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `matc_jacobi`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn var_temp_new(_typ: c_int, _nrow: c_int, _ncol: c_int) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `matc_jacobi`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mem_alloc(_size: size_t) -> *mut core::ffi::c_void {
    std::ptr::null_mut()
}

/// Not used by `matc_jacobi`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mem_free(_ptr: *mut core::ffi::c_void) {}

/// Fixed-seed symmetric positive definite matrix (diagonally dominant).
fn corpus(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let v = (next() >> 11) as f64 / (1u64 << 53) as f64;
            m[i * n + j] = (v - 0.5) * 2.0;
        }
    }
    // Symmetrize and make it diagonally dominant (positive definite).
    for i in 0..n {
        for j in (i + 1)..n {
            let s = 0.5 * (m[i * n + j] + m[j * n + i]);
            m[i * n + j] = s;
            m[j * n + i] = s;
        }
        m[i * n + i] = n as f64;
    }
    m
}

/// FNV-1a over the raw bytes of the result arrays.
fn digest(parts: &[&[f64]]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for p in parts {
        for v in p.iter() {
            for b in v.to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x100000001b3);
            }
        }
    }
    format!("{h:016x}")
}

fn child(n: usize, iters: usize) -> Result<String, String> {
    let a0 = corpus(n, 0x9E3779B97F4A7C15);
    let b0 = corpus(n, 0xD1B54A32D192ED03);
    let nn = n as c_int;
    let rtol: f64 = 1.0e-9;

    let mut a = a0.clone();
    let mut b = b0.clone();
    let mut x = vec![0.0f64; n * n];
    let mut eigv = vec![0.0f64; n];
    let mut d = vec![0.0f64; n];

    // Warm-up (not measured): first touch of the pages.
    // SAFETY: all five arrays hold n*n / n doubles.
    let rc = unsafe {
        gui_c_jacobi::matc_jacobi(
            a.as_mut_ptr(),
            b.as_mut_ptr(),
            x.as_mut_ptr(),
            eigv.as_mut_ptr(),
            d.as_mut_ptr(),
            nn,
            rtol,
        )
    };
    if rc == 0 {
        return Err("matc_jacobi rejected the corpus (rc=0)".to_string());
    }

    let t = harness::Timer::start();
    for _ in 0..iters {
        a.copy_from_slice(&a0);
        b.copy_from_slice(&b0);
        x.iter_mut().for_each(|v| *v = 0.0);
        eigv.iter_mut().for_each(|v| *v = 0.0);
        d.iter_mut().for_each(|v| *v = 0.0);
        // SAFETY: as above.
        unsafe {
            gui_c_jacobi::matc_jacobi(
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                x.as_mut_ptr(),
                eigv.as_mut_ptr(),
                d.as_mut_ptr(),
                nn,
                rtol,
            )
        };
    }
    let wall_ns = t.elapsed_ns();

    // Digest of the RESULT (before the overhead loop below rewrites the
    // scratch arrays).
    let dg = digest(&[&x, &eigv, &d]);

    // Harness overhead: the same array resets, no solver call.
    let t = harness::Timer::start();
    for _ in 0..iters {
        a.copy_from_slice(&a0);
        b.copy_from_slice(&b0);
        x.iter_mut().for_each(|v| *v = 0.0);
        eigv.iter_mut().for_each(|v| *v = 0.0);
        d.iter_mut().for_each(|v| *v = 0.0);
    }
    let harness_ns = t.elapsed_ns();

    let rss = harness::self_peak_rss_kb();
    Ok(harness::child_json(wall_ns, harness_ns, rss, &dg))
}

fn main() {
    let args = match harness::parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("jacobi_bench: {e}");
            std::process::exit(2);
        }
    };
    let r = match args.mode.as_str() {
        "child" => child(args.n, args.iters),
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
            eprintln!("jacobi_bench: {e}");
            std::process::exit(1);
        }
    }
}
