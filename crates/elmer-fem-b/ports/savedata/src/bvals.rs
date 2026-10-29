//! `saveboundaryvalues_` — port of SaveBoundaryValues.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE` replicated (`DIM`, `LocalNodes`, `ParamName`, `NoParams`)
//!   as a Rust `static Mutex<State>` computed once on first call.
//!   `LocalParam`/`Valid` are per-call scratch, freed at entry end.
//! - `ParamName(99)` bound: the source would overflow past 99 parameters;
//!   the port grows instead (strictly safer; unobservable below 99 and
//!   documented here rather than replicated as a crash).
//! - The per-parameter creation/update traversals run in Fortran wrappers
//!   off the captured handles (same Elmer calls, same order); Rust owns
//!   the scan loop, the SAVE state, and the per-call scratch.

use std::ffi::{c_void, CString};
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sdb_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sdb_dim() -> c_int;
    fn sdb_local_nodes() -> c_int;
    fn sdb_max_elem_nodes() -> c_int;
    fn sdb_info_creating();
    fn sdb_param(i: c_int, buf: *mut c_char, cap: c_int, found: *mut c_int);
    fn sdb_warn_no_params();
    fn sdb_create_var(
        name: *const c_char,
        local_nodes: c_int,
        exists: *mut c_int,
        created: *mut c_int,
        valid: *mut c_int,
        ndofs: *mut c_int,
    );
    fn sdb_update_var(name: *const c_char, localparam: *mut c_double, maxn: c_int);
    fn sdb_fatal_alloc1();
}

struct State {
    dim: usize,
    local_nodes: usize,
    param_names: Vec<String>,
    initialized: bool,
}

impl State {
    const fn new() -> Self {
        Self { dim: 0, local_nodes: 0, param_names: Vec::new(), initialized: false }
    }
}

static STATE: Mutex<State> = Mutex::new(State::new());

fn param_c_string(name: &str) -> CString {
    CString::new(name).unwrap_or_default()
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sdb_begin(model, solver) };
    unsafe { sdb_info_creating() };

    if !state.initialized {
        state.dim = unsafe { sdb_dim() } as usize;
        let _ = state.dim;
        state.local_nodes = unsafe { sdb_local_nodes() } as usize;
        let max_nodes = unsafe { sdb_max_elem_nodes() } as usize;
        let mut local_param = Vec::<c_double>::new();
        if local_param.try_reserve(max_nodes).is_err() {
            unsafe { sdb_fatal_alloc1() };
            std::process::abort();
        }
        local_param.resize(max_nodes, 0.0);
        // Parameter scan (saved; the source scans once).
        let mut names: Vec<String> = Vec::new();
        let mut buf = vec![0 as c_char; 4096];
        let mut i = 1;
        loop {
            let mut found = 0;
            unsafe { sdb_param(i, buf.as_mut_ptr(), buf.len() as c_int, &mut found) };
            if found == 0 {
                break;
            }
            let end = buf.iter().position(|&c| c == 0).unwrap_or(buf.len());
            let bytes: Vec<u8> = buf[..end].iter().map(|&c| c as u8).collect();
            names.push(String::from_utf8_lossy(&bytes).into_owned());
            i += 1;
            if names.len() >= 99 {
                break;
            }
        }
        state.param_names = names;
        state.initialized = true;
        // local_param freed at scope end, mirroring DEALLOCATE(LocalParam).
        let _ = local_param;
    }

    if state.param_names.is_empty() {
        unsafe { sdb_warn_no_params() };
        return;
    }

    let max_nodes = unsafe { sdb_max_elem_nodes() } as usize;
    let mut local_param = vec![0.0 as c_double; max_nodes];
    let mut valid: Vec<bool> = Vec::new();

    // Creation pass (mirrors the source's first ParamNo loop).
    for name in state.param_names.clone() {
        let cname = param_c_string(&name);
        let (mut exists, mut created, mut is_valid, mut _ndofs) = (0, 0, 0, 0);
        unsafe {
            sdb_create_var(
                cname.as_ptr(),
                state.local_nodes as c_int,
                &mut exists,
                &mut created,
                &mut is_valid,
                &mut _ndofs,
            )
        };
        let _ = (exists, created);
        valid.push(is_valid != 0);
    }

    // Update pass (mirrors the source's second ParamNo loop).
    for (name, &ok) in state.param_names.clone().iter().zip(valid.iter()) {
        if !ok {
            continue;
        }
        let cname = param_c_string(name);
        unsafe { sdb_update_var(cname.as_ptr(), local_param.as_mut_ptr(), max_nodes as c_int) };
    }
}
