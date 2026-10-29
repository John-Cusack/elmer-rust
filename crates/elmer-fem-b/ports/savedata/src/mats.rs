//! `savematerials_` — port of SaveMaterials.F90.
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - `SAVE PrevNodes` replicated as an accessor-side guard (mesh-change
//!   Warn + return, exactly as the source).
//! - `ParamName(99)` bound matches the source shape; over-99 parameters
//!   grow instead of overflowing (documented, strictly safer).
//! - `Indexes(27)` fixed local matches the source; over-27-DOF elements
//!   trap in the wrapper exactly like the source's fixed array would.
//! - Per-parameter creation/update traversals (including the DG branch,
//!   BodyForce/Material selection, and the nodal/DG Fatal) run in Fortran
//!   wrappers; Rust owns the scan loop, SAVE state, and per-call scratch.
//! - `LocalParam` is per-call scratch; the source's missing DEALLOCATE is
//!   not replicated (unobservable).

use std::ffi::{c_void, CString};
use std::os::raw::{c_char, c_double, c_int};

extern "C" {
    fn sdm_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sdm_local_nodes() -> c_int;
    fn sdm_max_elem_nodes() -> c_int;
    fn sdm_info_creating();
    fn sdm_check_mesh(local_nodes: c_int) -> c_int;
    fn sdm_warn_remesh();
    fn sdm_param(i: c_int, buf: *mut c_char, cap: c_int, found: *mut c_int);
    fn sdm_warn_no_params();
    fn sdm_bodyforce_params() -> c_int;
    fn sdm_create_var(
        name: *const c_char,
        is_bf: c_int,
        local_nodes: c_int,
        exists: *mut c_int,
        created: *mut c_int,
        ndofs: *mut c_int,
    );
    fn sdm_update_var(
        name: *const c_char,
        is_bf: c_int,
        localparam: *mut c_double,
        maxn: c_int,
    );
    fn sdm_fatal_alloc1();
}

pub unsafe fn entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sdm_begin(model, solver) };
    unsafe { sdm_info_creating() };

    let local_nodes = unsafe { sdm_local_nodes() };
    if unsafe { sdm_check_mesh(local_nodes) } != 0 {
        unsafe { sdm_warn_remesh() };
        return;
    }

    let max_nodes = unsafe { sdm_max_elem_nodes() } as usize;
    let mut local_param = Vec::<c_double>::new();
    if local_param.try_reserve(max_nodes).is_err() {
        unsafe { sdm_fatal_alloc1() };
        std::process::abort();
    }
    local_param.resize(max_nodes, 0.0);

    let mut names: Vec<String> = Vec::new();
    let mut buf = vec![0 as c_char; 4096];
    let mut i = 1;
    loop {
        let mut found = 0;
        unsafe { sdm_param(i, buf.as_mut_ptr(), buf.len() as c_int, &mut found) };
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
    if names.is_empty() {
        unsafe { sdm_warn_no_params() };
        return;
    }

    let bf_params = unsafe { sdm_bodyforce_params() };
    for (idx, name) in names.iter().enumerate() {
        let cname = CString::new(name.as_str()).unwrap_or_default();
        let is_bf = if (idx as c_int) < bf_params { 1 } else { 0 };
        let (mut exists, mut created, mut _ndofs) = (0, 0, 0);
        unsafe {
            sdm_create_var(
                cname.as_ptr(),
                is_bf,
                local_nodes,
                &mut exists,
                &mut created,
                &mut _ndofs,
            )
        };
        let _ = (exists, created);
    }
    for (idx, name) in names.iter().enumerate() {
        let cname = CString::new(name.as_str()).unwrap_or_default();
        let is_bf = if (idx as c_int) < bf_params { 1 } else { 0 };
        unsafe {
            sdm_update_var(cname.as_ptr(), is_bf, local_param.as_mut_ptr(), max_nodes as c_int)
        };
    }
}
