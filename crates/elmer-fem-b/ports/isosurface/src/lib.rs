//! Port of `fem/src/modules/Isosurface.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/Isosurface.F90` — two F77-global entries
//! (`isosurfacesolver_`, `isosurfacesolver_init_`) with twelve `CONTAINS`
//! internals. The emitted dynamic symbols are exactly those two (proven by
//! `nm -D` on the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/isosurface_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset.
//!
//! Granularity (documented, not incidental): the entry body's control flow
//! — phases, `GOTO 100` exits, the collective size reduction, reports,
//! savers, and the norm write — is transcribed statement-for-statement
//! into the phase calls below. The twelve `CONTAINS` internals stay
//! Fortran-side as whole module procedures because their dummies cannot
//! cross the C boundary at all (`TYPE(Element_t)` by value, assumed-shape
//! arrays, pointer-function results); moving those across would need the
//! callee-side descriptor ABI (the new-ABI refusal class). All
//! floating-point work (edge interpolation, magnitude `SQRT` sums, the
//! orientation cross products) and all formatted I/O stay in Fortran with
//! the same compiler, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors Isosurface.F90 statement for statement):
//! - `SAVE` replicated: `FixedSurface`, `SurfaceExist`, `Interpolant`,
//!   `InvPerm`, `IsoMesh`, `SaveParTmp`, plus the implicitly-`SAVE`'d
//!   `calls` counter, live as wrapper-side SAVE state with identical
//!   lifetimes; every other host variable is per-call automatic state the
//!   entry wrapper resets (nullify/deallocate guards) so re-entry matches
//!   the original's fresh-automatic semantics.
//! - `Model`, `dt`, and the transient flag are unused by the entry body
//!   beyond the calls the original makes (`Model` travels into the
//!   wrappers for `CurrentElement`, `Meshes`, and the `AnyBC/BodyForce`
//!   mask queries), exactly as in the original.
//! - The original has no `STAT=` allocation checks, so there are none to
//!   replicate; allocation behavior matches statement-for-statement.
//! - `OrigElements` (a pointer copy never read) is the only dropped
//!   statement, with no observable effect.
//! - `MPI_*` collectives stay in the wrappers behind the same
//!   `ParEnv % PEs > 1` guards; the serial proof never reaches them, and
//!   they resolve against the same `libmpi_stubs` as the reference.
//! - No Rust threads: the solver calls entries serially, and all shared
//!   state lives Fortran-side, so the Rust side holds no statics at all.

use std::ffi::c_void;
use std::os::raw::c_int;

extern "C" {
    fn iso_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn iso_info_start();
    fn iso_read_var() -> c_int;
    fn iso_read_levels() -> c_int;
    fn iso_fixed_remap() -> c_int;
    fn iso_prep();
    fn iso_detect_dim() -> c_int;
    fn iso_build_temporal();
    fn iso_create_mesh();
    fn iso_create_nodes() -> c_int;
    fn iso_create_surfaces() -> c_int;
    fn iso_parallel_info();
    fn iso_release_temp();
    fn iso_finish();
    fn iso_init_begin(h_solver: *mut c_void);
    fn iso_init_body();
}

/// `isosurfacesolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `IsosurfaceSolver_init` (default `Variable` when absent).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "isosurfacesolver_init_"]
pub unsafe extern "C" fn isosurface_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_void,
    _transient: *const c_void,
) {
    unsafe {
        iso_init_begin(solver);
        iso_init_body();
    }
}

/// `isosurfacesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// are unused by the body, exactly as in the original.
///
/// Each early branch (`missing variable`, `no levels`, `fixed-surface
/// remap`, `no 2D/3D path`) falls through to the collective size reduction
/// in `iso_finish`, mirroring the original's every-exit-goes-to-100 rule.
///
/// # Safety
///
/// Same contract as [`isosurface_solver_init_entry`]; additionally `model`
/// must be the Elmer `Model_t` address the loader passes.
#[export_name = "isosurfacesolver_"]
pub unsafe extern "C" fn isosurface_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_void,
    _transient: *const c_void,
) {
    unsafe {
        iso_begin(model, solver);
        iso_info_start();
        // Missing isosurface variable: Error, then the label-100 tail.
        if iso_read_var() != 0 {
            iso_finish();
            return;
        }
        // No isosurface value(s): Warn, then the label-100 tail.
        if iso_read_levels() != 0 {
            iso_finish();
            return;
        }
        // Fixed surface already created: remap, report saved sizes, tail.
        if iso_fixed_remap() != 0 {
            iso_finish();
            return;
        }
        iso_prep();
        // No 2D/3D leading dimension: Warn, then the label-100 tail.
        if iso_detect_dim() != 0 {
            iso_finish();
            return;
        }
        iso_build_temporal();
        iso_create_mesh();
        let _nodes: c_int = iso_create_nodes();
        let _surfaces: c_int = iso_create_surfaces();
        iso_parallel_info();
        iso_release_temp();
        iso_finish();
    }
}
