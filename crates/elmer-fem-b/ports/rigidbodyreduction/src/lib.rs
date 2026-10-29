//! Port of `fem/src/modules/RigidBodyReduction.F90` (Elmer @ 9f6af2f,
//! the campaign fork) for lane T021-B, at the Before-Linsolve plugin
//! boundary (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/RigidBodyReduction.F90` — one F77-global
//! INTEGER FUNCTION entry `rigidbody_` with six `CONTAINS` internals. The
//! emitted dynamic symbol is exactly `rigidbody_` (proven by `nm -D` on the
//! reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, Before-Linsolve gate shape): the entry is invoked
//! through `LoadMod.execlinsolveprocs`, whose explicit interface declares
//! the matrix dummy as `TYPE(Matrix_t), POINTER`. Under that interface
//! gfortran passes the address of the pointer variable (double
//! indirection — proven with a probe during triage), so the Rust entry
//! takes `mtr: *mut *mut c_void` and dereferences it once before handing
//! the `Matrix_t` address to the wrappers. `Model`/`Solver` cross as opaque
//! addresses; `b`/`x` are explicit-shape arrays; `n`/`DOFs` are INTEGERs by
//! reference; `Norm` is REAL by reference (updated in place by the reduced
//! solve, exactly as in the original); the INTEGER function result returns
//! in the register (`1` = reduced and solved, `0` = no rigid bodies).
//!
//! Every component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/rigidbodyreduction_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `COUNT`, `MERGE`), CRS algebra, formatted I/O (`WRITE`/`PRINT`),
//! and whole-array updates stay in Fortran wrappers, so bit-identity holds
//! by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors RigidBodyReduction.F90 section for section):
//! - `SAVE` semantics replicated: the whole original SAVE list
//!   (`C`, `C_Trans`, `RedA`, `AllocationsDone`, `u`, `RealNodeTypes`,
//!   `Permutation`, `CenterOfRigidBody`, `NbrRigids`, `RigidIndex`, `f`,
//!   `NodeTypes`, `Verbose`) lives as wrapper-side SAVE state with
//!   identical lifetimes; cross-section scratch (`Perm`, `BMatrix`,
//!   `FVector`, `XVector`, `LocalEigenVectors`, `TotTime`,
//!   `RowsInUnityBlock`, `NEigen`, `DOFs`, `istat`) is wrapper-side state
//!   written before read on every path that reaches the read.
//! - The only Rust-visible branch is the `NbrRigids == 0` early return
//!   (`Warn` + `0`, then the core proceeds with the normal solve — exactly
//!   as the original's `RETURN`); otherwise the wrapper sections run in the
//!   original statement order. `RowsInUnityBlock` on a `MeshChanged`
//!   re-entry keeps the first-call value (the original leaves it
//!   undefined there); the proof never takes `MeshChanged`, so every
//!   proven path is exact.
//! - Allocation failure calls the same `Fatal(...)` via the `STAT=istat`
//!   checks, instead of aborting.
//! - No Rust threads: the hook runs serially, so the `Mutex` is only a
//!   `static`-sharing guard, never contended across threads. The solve
//!   itself never re-enters Rust callbacks, so holding the guard across
//!   the call sequence is safe.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn rbr_scan(
        h_model: *const c_void,
        h_solver: *const c_void,
        nrig: *mut c_int,
    );
    fn rbr_setup(
        h_model: *const c_void,
        h_solver: *const c_void,
        h_a: *const c_void,
        b: *mut c_double,
        n: c_int,
        ndofs: c_int,
    );
    fn rbr_fill(h_a: *const c_void, b: *mut c_double, n: c_int);
    fn rbr_centers(h_model: *const c_void, h_solver: *const c_void);
    fn rbr_cstruct(
        h_model: *const c_void,
        h_solver: *const c_void,
        h_a: *const c_void,
        ndofs: c_int,
    );
    fn rbr_cvals(
        h_model: *const c_void,
        h_solver: *const c_void,
        h_a: *const c_void,
        ndofs: c_int,
    );
    fn rbr_eigsetup(
        h_solver: *const c_void,
        eig: *mut c_int,
        neigen: *mut c_int,
    );
    fn rbr_multiply(h_a: *const c_void, eig: c_int);
    fn rbr_optimize(h_solver: *const c_void, h_a: *const c_void);
    fn rbr_solve(
        h_solver: *const c_void,
        h_a: *const c_void,
        norm: *mut c_double,
        ndofs: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn rbr_mapback(
        h_model: *const c_void,
        h_solver: *const c_void,
        h_a: *const c_void,
        b: *mut c_double,
        x: *mut c_double,
        n: c_int,
        ndofs: c_int,
        eig: c_int,
    );
}

/// Serial-sharing guard for the entry (never contended: Elmer calls the
/// hook serially, and the port keeps no Rust-side per-call state — all
/// `SAVE` state lives in the Fortran wrappers).
static STATE: Mutex<()> = Mutex::new(());

/// `rigidbody_` — the Before-Linsolve hook. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, A, b, x, n, DOFs, Norm)`; `mtr` arrives
/// double-indirected (address of the `Matrix_t` pointer variable) under
/// the explicit `POINTER` interface and is dereferenced once here.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes; `mtr` must point at the live `Matrix_t` pointer;
/// `b`/`x` must address `n` doubles; `n`/`dofs`/`nrm` must be valid reads
/// (`nrm` is updated in place by the reduced solve).
#[allow(clippy::too_many_arguments)]
#[export_name = "rigidbody_"]
pub unsafe extern "C" fn rigidbody_entry(
    model: *mut c_void,
    solver: *mut c_void,
    mtr: *mut *mut c_void,
    b: *mut c_double,
    x: *mut c_double,
    n: *const c_int,
    dofs: *const c_int,
    nrm: *mut c_double,
) -> c_int {
    let _guard = STATE.lock().expect("rigidbody guard poisoned");
    // Double-indirected Matrix POINTER dummy: the caller passes the
    // address of its pointer variable; the Matrix_t address is stored there.
    let amat = *mtr;
    let nn = *n;
    let nd = *dofs;

    let mut nrig: c_int = 0;
    rbr_scan(model, solver, &mut nrig);
    if nrig == 0 {
        return 0;
    }
    rbr_setup(model, solver, amat, b, nn, nd);
    rbr_fill(amat, b, nn);
    rbr_centers(model, solver);
    rbr_cstruct(model, solver, amat, nd);
    rbr_cvals(model, solver, amat, nd);
    let mut eig: c_int = 0;
    let mut neigen: c_int = 0;
    rbr_eigsetup(solver, &mut eig, &mut neigen);
    rbr_multiply(amat, eig);
    rbr_optimize(solver, amat);
    rbr_solve(solver, amat, nrm, nd);
    rbr_mapback(model, solver, amat, b, x, nn, nd, eig);
    1
}
