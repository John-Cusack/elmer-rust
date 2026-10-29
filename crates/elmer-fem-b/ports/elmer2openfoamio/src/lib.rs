//! Port of `fem/src/modules/Elmer2OpenFoamIO.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/Elmer2OpenFoamIO.F90` — one F77-global
//! entry `elmer2openfoamwrite_` with four `CONTAINS` internals. The
//! emitted dynamic symbol is exactly `elmer2openfoamwrite_` (proven by
//! `nm -D` on the reference `.so`); the internals emit only local
//! symbols. Note the library-vs-entry name difference: the campaign
//! invocation is `Procedure = "Elmer2OpenFoamIO" "Elmer2OpenFoamWrite"`.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/elmer2openfoamio_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MINVAL`/`MAXVAL`/`SUM`), formatted reads/writes (`READ`/`WRITE`),
//! and whole-array zeroing stay in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors Elmer2OpenFoamWrite statement for
//! statement):
//! - `SAVE` semantics replicated: the temporal-mesh pointer, the mapped
//!   variable pointer, and `NoDir` live as wrapper-side SAVE state with
//!   identical lifetimes (re-entry releases the saved mesh exactly like
//!   the original's `ReleaseMesh` + `DEALLOCATE`); the `Visited` flag
//!   lives in the Rust `static` below. No Rust threads: the solver calls
//!   entries serially, so the `Mutex` is only a `static`-sharing guard,
//!   never contended.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body except
//!   for the `Model % Simulation` projector-flag save/restore, exactly
//!   as in the original (received and otherwise ignored).
//! - The `Elmer Center Filename` early `RETURN` path skips the projector
//!   restore, the mesh-pointer restore, and the `Visited` set, exactly as
//!   in the original.
//! - `PassiveCoord` reads 0 when the keyword is absent (the
//!   `ListGetInteger` default), matching the original's defined function
//!   result.
//! - Interpolation itself happens library-side inside `VariableGet`
//!   (`InterpolateMeshToMesh` on the temporal mesh); the port reproduces
//!   the call sequence only.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn e2o_begin(h_model: *const c_void, h_solver: *const c_void);
    fn e2o_info_head();
    fn e2o_params();
    fn e2o_target(found: *mut c_int);
    fn e2o_fatal_no_target();
    fn e2o_proj_save(psave: *mut c_int, pfound: *mut c_int);
    fn e2o_getmesh();
    fn e2o_var_get() -> c_int;
    fn e2o_fatal_no_var();
    fn e2o_release_saved();
    fn e2o_allocate();
    fn e2o_passive(val: *mut c_int, found: *mut c_int);
    fn e2o_meshdim() -> c_int;
    fn e2o_warn_dim();
    fn e2o_center_name(found: *mut c_int);
    fn e2o_write_centers();
    fn e2o_basedir(found: *mut c_int);
    fn e2o_info_basedir();
    fn e2o_fatal_no_basedir();
    fn e2o_has_mesh1() -> c_int;
    fn e2o_blocks();
    fn e2o_reduce();
    fn e2o_nodir() -> c_int;
    fn e2o_info_nodir();
    fn e2o_mype() -> c_int;
    fn e2o_info_treating(i: c_int);
    fn e2o_dirname(i: c_int, found: *mut c_int);
    fn e2o_fatal_no_key();
    fn e2o_filename_c();
    fn e2o_info_project();
    fn e2o_create();
    fn e2o_current_to_of();
    fn e2o_zero_passive();
    fn e2o_info_mapping();
    fn e2o_ofvar_get();
    fn e2o_solver_var();
    fn e2o_info_stats();
    fn e2o_write_field();
    fn e2o_release_of();
    fn e2o_current_to_mesh();
    fn e2o_proj_restore(psave: c_int, pfound: c_int);
    fn e2o_info_done();
}

/// First-visit flag: the original's `SAVE`d `Visited`, initially
/// `.FALSE.`. Set at the end of a full pass only.
struct E2oState {
    visited: bool,
}

impl E2oState {
    const fn new() -> Self {
        Self { visited: false }
    }
}

static STATE: Mutex<E2oState> = Mutex::new(E2oState::new());

/// `elmer2openfoamwrite_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "elmer2openfoamwrite_"]
pub unsafe extern "C" fn elmer2openfoamwrite_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut st = STATE.lock().unwrap();
    e2o_begin(model, solver);
    e2o_info_head();
    e2o_params();
    let mut found: c_int = 0;
    e2o_target(&mut found);
    if found == 0 {
        e2o_fatal_no_target();
    }
    let mut psave: c_int = 0;
    let mut pfound: c_int = 0;
    e2o_proj_save(&mut psave, &mut pfound);
    e2o_getmesh();
    if e2o_var_get() == 0 {
        e2o_fatal_no_var();
    }
    if st.visited {
        e2o_release_saved();
    }
    e2o_allocate();
    let mut passive: c_int = 0;
    let mut pcfound: c_int = 0;
    e2o_passive(&mut passive, &mut pcfound);
    if pcfound == 0 && e2o_meshdim() < 3 {
        e2o_warn_dim();
    }
    let mut cfound: c_int = 0;
    e2o_center_name(&mut cfound);
    if cfound != 0 {
        e2o_write_centers();
        return;
    }
    let mut bfound: c_int = 0;
    e2o_basedir(&mut bfound);
    if bfound != 0 {
        e2o_info_basedir();
    } else {
        e2o_fatal_no_basedir();
    }
    if e2o_has_mesh1() == 0 {
        e2o_blocks();
    }
    e2o_reduce();
    let nodir = e2o_nodir();
    e2o_info_nodir();
    for i in 1..=nodir {
        if e2o_mype() == 0 {
            if nodir > 1 {
                e2o_info_treating(i);
            }
            let mut dfound: c_int = 0;
            e2o_dirname(i, &mut dfound);
            if dfound == 0 {
                e2o_fatal_no_key();
            }
            e2o_filename_c();
            e2o_info_project();
        }
        e2o_create();
        e2o_current_to_of();
        e2o_zero_passive();
        e2o_info_mapping();
        e2o_ofvar_get();
        e2o_solver_var();
        if e2o_mype() == 0 {
            e2o_info_stats();
            e2o_write_field();
        }
        if i < nodir {
            e2o_release_of();
        }
    }
    e2o_current_to_mesh();
    e2o_proj_restore(psave, pfound);
    st.visited = true;
    e2o_info_done();
}
