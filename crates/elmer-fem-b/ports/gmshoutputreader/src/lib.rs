//! Port of `fem/src/modules/GmshOutputReader.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/GmshOutputReader.F90` — two F77-global
//! entries (`gmshoutputreader_`, `gmshoutputreader_init_`); the
//! `CONTAINS` internals (`ReadSingleGmshFile`,
//! `InterpolateFromGmshFile`) emit only local symbols. The emitted
//! dynamic symbols are exactly those two (proven by `nm -D` on the
//! reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/gmshoutputreader_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. All floating-point work
//! (`MINVAL`/`MAXVAL` reductions, coordinate shifts), list-directed file
//! I/O, whole-array updates, and formatted `Info` output stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors GmshOutputReader.F90 statement for
//! statement):
//! - The original's label-10 `GOTO` sizing sweep becomes the Rust
//!   `loop { gmr_pass(); if gmr_need_alloc() { gmr_allocate() } }`
//!   with identical call sequence (counting pass, then allocation, then
//!   the fill pass); the `REWIND`/`CLOSE` discipline and the `CYCLE`
//!   paths stay in the Fortran sweep body.
//! - `SAVE` semantics replicated: the `AllocationsDone` flag,
//!   `MaxElemDim`/`MaxElemNodes`, the `ReadSingleGmshFile` SAVE
//!   `GmshVer`, and the `VarArray(20)` table live as wrapper-side SAVE
//!   state with identical lifetimes; the entry resets the sizing state
//!   on every call exactly as the original's executable prelude does,
//!   so repeated invocation re-reads deterministically.
//! - The `MaskPerm => ...; NULLIFY(MaskPerm)` sequence is transcribed
//!   verbatim (pointer assigned then nullified, hence never associated),
//!   matching the original's effective behavior.
//! - The `$NodeData` integer-tags line is read into the mesh-dimension
//!   slot exactly as the original does (its `! dim` comment), and the
//!   `MAX(MaxElemDim, MeshDim)` fixup is kept.
//! - `Model`, `dt`, and the transient flag are unused by the entry body
//!   itself (received and stored as handles; `Model` reaches only the
//!   mask helper inside Fortran), exactly as in the original.
//! - No Rust threads and no Rust-held state: the solver entry runs
//!   serially, so no `Mutex` is needed; all state lives wrapper-side.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn gmr_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn gmr_setup() -> c_int;
    fn gmr_sweep_begin();
    fn gmr_pass();
    fn gmr_need_alloc() -> c_int;
    fn gmr_allocate();
    fn gmr_finalize();
    fn gmr_interpolate();
    fn gmr_init_begin(h_solver: *mut c_void);
    fn gmr_init_body();
}

/// `gmshoutputreader_init_` — the `_Init` entry Elmer calls once at
/// solver setup. Mirrors `GmshOutputReader_Init` (marks the solver
/// matrix-free).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "gmshoutputreader_init_"]
pub unsafe extern "C" fn gmshoutputreader_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        gmr_init_begin(solver);
        gmr_init_body();
    }
}

/// `gmshoutputreader_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`gmshoutputreader_init_entry`]; `model` is stored
/// as an opaque handle for the mask helper.
#[export_name = "gmshoutputreader_"]
pub unsafe extern "C" fn gmshoutputreader_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        gmr_begin(model, solver);
        if gmr_setup() == 1 {
            return;
        }
        gmr_sweep_begin();
        loop {
            gmr_pass();
            if gmr_need_alloc() == 1 {
                gmr_allocate();
            } else {
                break;
            }
        }
        gmr_finalize();
        gmr_interpolate();
    }
}
