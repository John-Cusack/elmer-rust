//! Port of `fem/src/modules/Mesh2MeshSolver.F90` (Elmer @ 9f6af2f85,
//! the lane fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/Mesh2MeshSolver.F90` — two F77-global
//! entries `mesh2meshsolver_init0_` and `mesh2meshsolver_`. The emitted
//! dynamic symbols are exactly those two (proven by `nm -D` on the
//! reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/mesh2meshsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive whole-array
//! updates (`Var % Values = Var % Values + BaseSol`, base-value copies),
//! `MAXVAL` reductions, formatted `WRITE`s and `PRINT`s stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors Mesh2MeshSolver.F90 statement for
//! statement):
//! - The original has no `SAVE` state; its per-call allocatables
//!   (`FromVars` list, `BaseSol` store) are mirrored as wrapper-side
//!   SAVE handles that are reset on every entry, replicating the
//!   original's automatic (de)allocation lifetimes.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies,
//!   exactly as in the original (received and ignored).
//! - The `NoVar` count loop, the `Target Mesh Solver Index` /
//!   first-different-mesh branches, the hand-picked/all-fields branch,
//!   the additive branch, and the validation branch run in the same
//!   order with the same conditions; only the bodies live in Fortran.
//! - The disabled `#if 0` block of the original is omitted (disabled in
//!   the source too); the unused `MaskName`/`TmpSol`/`GotMaskName`
//!   locals have no observable effect and are omitted with it.
//! - Allocation failure aborts in the original (no `STAT=` check), so
//!   the mirroring `ALLOCATE`s stay in Fortran and abort identically.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

 extern "C" {
     fn m2ma_init0();
     fn m2ma_begin(h_solver: *const c_void);
     fn m2ma_resolve_meshes();
     fn m2ma_build_vars(novar: *mut c_int, additive: *mut c_int);
     fn m2ma_parenv();
     fn m2ma_interpolate();
     fn m2ma_additive();
     fn m2ma_validate(novar: c_int);
 }

 /// Sharing guard for the wrapper-side SAVE handles (never contended:
 /// the entries run serially).
 static STATE: Mutex<()> = Mutex::new(());

/// `mesh2meshsolver_init0_` — flags the solver matrix-free with a local
/// mesh copy. Argument order and by-reference travel match the Fortran
/// dummy list `(Model, Solver, dt, TransientSimulation)`; only the
/// solver context (via `GetSolverParams`) is used, as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "mesh2meshsolver_init0_"]
pub unsafe extern "C" fn mesh2meshsolver_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe {
        m2ma_begin(solver);
        m2ma_init0();
    }
}

/// `mesh2meshsolver_` — maps fields between meshes. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "mesh2meshsolver_"]
pub unsafe extern "C" fn mesh2meshsolver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe {
        m2ma_begin(solver);
        m2ma_resolve_meshes();
        // Hand-picked list (`NoVar > 0`) or whole-mesh mapping.
        let (mut novar, mut additive) = (0, 0);
        m2ma_build_vars(&mut novar, &mut additive);
        m2ma_parenv();
        m2ma_interpolate();
        // Additive variables are added back onto the base field.
        if additive != 0 {
            m2ma_additive();
        }
        m2ma_validate(novar);
    }
}
