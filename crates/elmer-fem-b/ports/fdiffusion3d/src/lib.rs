//! Port of `fem/src/modules/FDiffusion3D.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/FDiffusion3D.F90` — one F77-global entry
//! `fourierdiffusion3dsolver_` with three `CONTAINS` internals
//! (`LocalMatrix`, `LocalMatrixBoundary`, `LorentzForceAve` +
//! `ComputeLorentz`). The emitted dynamic symbol is exactly
//! `fourierdiffusion3dsolver_` (proven by `nm -D` on the reference
//! `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/fdiffusion3d_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (complex 3-D Galerkin assembly, `SUM` reductions,
//! `WRITE(Message, ...)` timing reports) stays in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FDiffusion3D.F90 statement for statement):
//! - `SAVE` semantics replicated: the entry's `SAVE` arrays
//!   (`LocalStiffMatrix`, `Load`, `LocalForce`, `Work`, `Impedance`,
//!   `Conductivity`, `Permeability`, `ElementNodes`, `AllocationsDone`,
//!   `LocalNodes`, `LrF`) live as wrapper-side SAVE state with identical
//!   lifetimes; allocation (with the same `STAT=istat` →
//!   `Fatal('FourierDiffusion3DSolve', ...)` check, including the
//!   original's `LrF`-is-never-deallocated shape) runs in the wrapper.
//!   The internal procedures' automatic arrays are per-call locals
//!   there, as in the original.
//! - Only bulk elements assemble (the original has no boundary loop;
//!   `LocalMatrixBoundary` is dead code there and stays an uncalled
//!   wrapper-internal routine here).
//! - Reference-UB, author's evident intent: the original's `LocalMatrix`
//!   initializes `Metric` to identity but never calls
//!   `CoordinateSystemInfo` on the Cartesian path, so in Cartesian
//!   coordinates it reads uninitialized stack arrays `Symb`/`dSymb` in
//!   every element assembly (the reference is nondeterministic across
//!   processes there — the same binary disagrees with itself run to
//!   run). `CoordinateSystemInfo`'s Cartesian branch sets `Symbols = 0`,
//!   `dSymbols = 0`, so the wrapper zeroes both next to the `Metric`
//!   init. This is bit-identical with the reference wherever the
//!   reference is defined (every non-Cartesian path calls
//!   `CoordinateSystemInfo`, which overwrites all four arrays); the
//!   proof differentials run the defined `Cylindrical` path.
//! - `LorentzForceAve` reads the last bulk element's `Material` pointer,
//!   exactly as the original's host association does (the wrapper saves
//!   each bulk element's material while assembling).
//! - The dead `Norm = Solver % Variable % Norm` read is skipped (no side
//!   effects); `Model`, `dt`, `TransientSimulation` cross as opaque
//!   handles, as in the original.
//! - No Rust threads: the solver entry runs serially (no `Mutex` needed —
//!   all shared state lives in the Fortran wrapper, as in the
//!   divergencesolver port).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fd3_begin(h_model: *const c_void, h_solver: *const c_void) -> c_int;
    fn fd3_alloc();
    fn fd3_params() -> c_int;
    fn fd3_nbulk() -> c_int;
    fn fd3_iter_head(iter: c_int);
    fn fd3_bulk(t: c_int);
    fn fd3_iter_tail(iter: c_int) -> c_int;
    fn fd3_lorentz();
}

/// `fourierdiffusion3dsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "fourierdiffusion3dsolver_"]
pub unsafe extern "C" fn fourierdiffusion3d_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN, then COMPLEX + guards.
    if unsafe { fd3_begin(model, solver) } != 0 {
        return;
    }
    unsafe {
        // First-mesh SAVE allocation (MeshChanged-aware, as the original).
        fd3_alloc();
        // Nonlinear iteration count, Equation name, Bubbles, base frequency.
        let nl = fd3_params();
        let nbulk = fd3_nbulk();
        for iter in 1..=nl {
            fd3_iter_head(iter);
            for t in 1..=nbulk {
                fd3_bulk(t);
            }
            // IF (Solver % Variable % NonlinConverged == 1) EXIT.
            if fd3_iter_tail(iter) != 0 {
                break;
            }
        }
        // Time-averaged Lorentz force accumulation.
        fd3_lorentz();
    }
}
