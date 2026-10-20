//! A/B bench for the `matrix` port (lane-gui-c, unit
//! `c:ElmerGUI/matc/src/matrix.c` — the matc matrix REPL helpers).
//!
//! Measured workload: `mtr_sum` (column sums of an n x n matrix). The
//! corpus is a fixed-seed matrix; the digest (FNV-1a over the result
//! vector) is the equivalence proof — the harness refuses to report a
//! speed when the two arms disagree.
//!
//! Link stubs: the port crate imports the sibling matc units. `mtr_sum`
//! itself needs only `var_temp_new` (the zeroed result block) and the
//! error path, but the linker needs every imported symbol, so the bench
//! provides them. `var_temp_new` is a real implementation here (calloc'd
//! zeroed storage, exactly like the matc variable allocator) because the
//! measured result is read back from it.

use elmer_bench_harness as harness;
use gui_c_matrix::{c_char, c_int, c_double, size_t, MATRIX, VARIABLE, TYPE_DOUBLE, CommandSub};

extern "C" {
    fn calloc(n: size_t, size: size_t) -> *mut core::ffi::c_void;
    fn free(ptr: *mut core::ffi::c_void);
}

/// matc error path: the original unwinds with longjmp; the bench aborts.
#[no_mangle]
pub unsafe extern "C" fn error_matc(_fmt: *const c_char, _: ...) -> ! {
    eprintln!("matrix_bench: error_matc reached");
    std::process::abort();
}

/// The matc temp allocator: a zeroed matrix block (calloc), like the
/// original. The measured result is read back from it, so this is a real
/// implementation and not a stub.
#[no_mangle]
pub unsafe extern "C" fn var_temp_new(typ: c_int, nrow: c_int, ncol: c_int) -> *mut VARIABLE {
    let mat = calloc(1, std::mem::size_of::<MATRIX>()) as *mut MATRIX;
    (*mat).typ = typ;
    (*mat).refcount = 1;
    (*mat).nrow = nrow;
    (*mat).ncol = ncol;
    (*mat).data =
        calloc((nrow as size_t) * (ncol as size_t), std::mem::size_of::<c_double>()) as *mut c_double;
    let var = calloc(1, std::mem::size_of::<VARIABLE>()) as *mut VARIABLE;
    (*var).next = std::ptr::null_mut();
    (*var).name = std::ptr::null_mut();
    (*var).changed = 0;
    (*var).this = mat;
    var
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn urand(_iy: *mut c_int) -> c_double {
    0.0
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_det(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_inv(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_eig(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_jacob(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_LUD(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn mtr_hesse(_v: *mut VARIABLE) -> *mut VARIABLE {
    std::ptr::null_mut()
}

/// Not used by `mtr_sum`; present so the crate links.
#[no_mangle]
pub unsafe extern "C" fn com_init(
    _nm: *const c_char,
    _pw: c_int,
    _ce: c_int,
    _sub: CommandSub,
    _minp: c_int,
    _maxp: c_int,
    _help: *const c_char,
) {
}

/// Free a temp block created by `var_temp_new`.
unsafe fn free_temp(var: *mut VARIABLE) {
    let mat = (*var).this;
    free((*mat).data as *mut core::ffi::c_void);
    free(mat as *mut core::ffi::c_void);
    free(var as *mut core::ffi::c_void);
}

/// Fixed-seed corpus matrix.
fn corpus(n: usize) -> Vec<f64> {
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut m = vec![0.0f64; n * n];
    for v in m.iter_mut() {
        *v = (next() >> 11) as f64 / (1u64 << 53) as f64;
    }
    m
}

/// FNV-1a over the raw bytes of the result vector.
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

fn child(n: usize, iters: usize) -> Result<String, String> {
    let data = corpus(n);
    let nn = n as c_int;

    // The argument variable: a matrix block holding the corpus.
    // SAFETY: the matc layout is repr(C) and the data block holds n*n doubles.
    let in_mat = unsafe { calloc(1, std::mem::size_of::<MATRIX>()) } as *mut MATRIX;
    unsafe {
        (*in_mat).typ = TYPE_DOUBLE;
        (*in_mat).refcount = 1;
        (*in_mat).nrow = nn;
        (*in_mat).ncol = nn;
    }
    let in_data = data.as_ptr() as *mut c_double;
    unsafe { (*in_mat).data = in_data };
    let in_var = unsafe { calloc(1, std::mem::size_of::<VARIABLE>()) } as *mut VARIABLE;
    unsafe {
        (*in_var).next = std::ptr::null_mut();
        (*in_var).name = std::ptr::null_mut();
        (*in_var).changed = 0;
        (*in_var).this = in_mat;
    }

    // Warm-up (not measured): first touch of the code path.
    // SAFETY: in_var holds a live n x n matrix block.
    let warm = unsafe { gui_c_matrix::mtr_sum(in_var) };
    // SAFETY: warm is a live temp created by var_temp_new.
    unsafe { free_temp(warm) };

    let mut out = vec![0.0f64; n];
    let t = harness::Timer::start();
    for _ in 0..iters {
        // SAFETY: as above.
        let res = unsafe { gui_c_matrix::mtr_sum(in_var) };
        // SAFETY: res is a live temp created by var_temp_new; the result
        // vector holds n doubles (one per column).
        unsafe {
            std::ptr::copy_nonoverlapping((*(*res).this).data, out.as_mut_ptr(), n);
            free_temp(res);
        }
    }
    let wall_ns = t.elapsed_ns();

    // Digest of the RESULT.
    let dg = digest(&out);

    // Harness overhead: the result block the port allocates per call,
    // without the sum itself.
    let t = harness::Timer::start();
    for _ in 0..iters {
        // SAFETY: the temp is a live block created by var_temp_new.
        unsafe {
            let tmp = var_temp_new(TYPE_DOUBLE, 1, nn);
            free_temp(tmp);
        }
    }
    let harness_ns = t.elapsed_ns();

    let rss = harness::self_peak_rss_kb();
    Ok(harness::child_json(wall_ns, harness_ns, rss, &dg))
}

fn main() {
    let args = match harness::parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("matrix_bench: {e}");
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
            eprintln!("matrix_bench: {e}");
            std::process::exit(1);
        }
    }
}
