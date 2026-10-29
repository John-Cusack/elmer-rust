//! Port of `fem/src/modules/TurbulenceModels/KESolverLegacy.F90` (the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/TurbulenceModels/KESolverLegacy.F90` —
//! two F77-global entries (`kesolverlegacy_`, `kesolverlegacy_init_`)
//! with two `CONTAINS` internals. The emitted dynamic symbols are exactly
//! those two (proven by `nm -D` on the reference `.so`); the internals
//! emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/kesolverlegacy_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the LocalMatrix Galerkin assembly, the RNG secant
//! loop, the K/E clip floors, the bubble-history shifts) and all
//! formatted I/O (`WRITE(Message,*)`/`Info`, the RNG `PRINT*`) stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors KESolverLegacy.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: the first-mesh SAVE arrays
//!   (sized by `Model % MaxElementNodes` and `Mesh % NumberOfNodes`) and
//!   the per-element bubble history `bx`/`bxprev` live as wrapper-side
//!   SAVE state with identical lifetimes (allocated once, never resized;
//!   the timestep shift runs once per timestep, guarded by the same
//!   `GetTimestep() /= BubbleTimestep` check). The Rust `Mutex` is only a
//!   `static`-sharing guard, never contended across threads.
//! - Allocation failure calls the same `Fatal('KESolver', ...)` via the
//!   `STAT=istat` check, instead of aborting.
//! - Dead stores are documented skips, never executed: `at`/`at0`,
//!   `NewtonLinearization`, `ForceVector`, `KEpsilon`, `ElectricCurrent`,
//!   `S`, `C`, `k1`, `k2`, `eq_id`, `bf_id` (like DivergenceSolver's
//!   `Visited`). `body_id` is compared but never updated in the original
//!   (stays `-1`); transcribed exactly. The entry-level `C0` is allocated
//!   but never read (the internal has its own local `C0`); replicated.
//! - `dt`/`TransientSimulation` cross by value; `Model` is received and
//!   stored (its components are read in wrappers, as in the original).
//! - `KEWall` stays an implicit external exactly as in the original (it
//!   resolves from libelmersolver at load, `U kewall_` in both arms).
//! - No Rust threads: the solver calls entries serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn ke_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, transient: c_int);
    fn ke_prologue() -> c_int;
    fn ke_params();
    fn ke_noslip();
    fn ke_alloc();
    fn ke_niter() -> c_int;
    fn ke_iter_head(iter: c_int) -> c_int;
    fn ke_bulk_elem(t: c_int);
    fn ke_bulk_tail();
    fn ke_wall_zero();
    fn ke_boundary();
    fn ke_solve_tail() -> c_int;
    fn ke_finish();
    fn ke_init_begin(h_solver: *mut c_void, transient: c_int);
    fn ke_init_body();
}

static STATE: Mutex<()> = Mutex::new(());

/// `kesolverlegacy_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `KESolverLegacy_Init` (default `Variable`, bubble
/// timestep minimums).
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null.
#[export_name = "kesolverlegacy_init_"]
pub unsafe extern "C" fn kesolverlegacy_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let _g = STATE.lock().unwrap();
        let tr = if transient.is_null() { 0 } else { *transient };
        ke_init_begin(solver, tr);
        ke_init_body();
    }
}

/// `kesolverlegacy_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`kesolverlegacy_init_entry`]; `dt` and `transient`
/// are read (the body consumes them).
#[export_name = "kesolverlegacy_"]
pub unsafe extern "C" fn kesolverlegacy_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        let _g = STATE.lock().unwrap();
        let dtv = *dt;
        let tr = *transient;
        ke_begin(model, solver, dtv, tr);
        if ke_prologue() != 0 {
            return;
        }
        ke_params();
        ke_noslip();
        ke_alloc();
        let niter = ke_niter();
        for iter in 1..=niter {
            let noactive = ke_iter_head(iter);
            for t in 1..=noactive {
                ke_bulk_elem(t);
            }
            ke_bulk_tail();
            ke_wall_zero();
            ke_boundary();
            if ke_solve_tail() != 0 {
                break;
            }
        }
        ke_finish();
    }
}
