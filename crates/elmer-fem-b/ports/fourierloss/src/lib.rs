//! Port of `fem/src/modules/FourierLoss.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/FourierLoss.F90` — two F77-global entries
//! (`fourierlosssolver_`, `fourierlosssolver_init_`) with three `CONTAINS`
//! internals (`LocalFourierTransform`, `BulkAssembly`, `CommunicateLosess`).
//! The emitted dynamic symbols are exactly those two (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/fourierloss_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the DFT `SIN`/`COS` coefficient updates, the
//! Galerkin assembly, `SQRT`/`SUM` reductions, the `RealTime`-based
//! reports) and all formatted/file I/O stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors FourierLoss.F90 statement for statement):
//! - `SAVE` semantics replicated explicitly. The entry-level SAVE state
//!   (`TimesVisited`, `CompLoss`, `SeriesLoss`, `BodyLoss`, the
//!   `CompVars`/`CompVarsE`/`FourierVars` variable-pointer arrays, the
//!   `OtherRhs` staging array, the `SaveRhs` pointer) and the
//!   `LocalFourierTransform` SAVE state (`time`, `ratio`, `tcycle`,
//!   `InitPending`, `t0`, `PreviousCycle`) live as wrapper-side SAVE state
//!   with identical lifetimes; the `BulkAssembly` `SAVE Nodes` likewise
//!   lives wrapper-side (consumed before the next fill, as in the
//!   original). The Rust side holds no state at all, so no `Mutex` is
//!   needed: the solver entries run serially, exactly as in the original.
//! - The `_Init` entry transcribes the obsolescence guard (`PRINT *` +
//!   `Fatal`), the `Scalars Prefix`/`Target Variable` reads, the
//!   `Harmonic Loss Coefficient i` count loop, the `res:` simulation-list
//!   seeding, and the `Exported Variable` registrations verbatim.
//! - Early `RETURN` paths map to entry returns: empty `Solver % Variable`
//!   permutation, and a DFT step that does not complete a cycle. All
//!   `Fatal` paths fire inside the wrappers with the original caller
//!   names (`'FourierLossSolver'`, `'FourierLosses'`, `'FourierLoss'`).
//! - The `AV`/`Direct` target-field dispatch, the `TargetSolverPtr`
//!   search, the `FourierVars` creation loop (`VariableAddVector`), the
//!   `OtherRhs` multi-`Rhs` pointer swaps, the elemental `LUdecomp`/
//!   `LUSolve` projection through `DGIndexes`, the nodal-loss gather, and
//!   the `CommunicateLosess` parallel reductions, screen reports, and
//!   `OPEN`/`WRITE` loss-file outputs are transcribed verbatim in the
//!   wrappers, including the `InfoActive(25)` debug branch.
//! - Allocation is per the original: entry SAVE arrays allocate once and
//!   persist; `BulkAssembly` temporaries allocate/deallocate per call —
//!   all inside the wrappers with the original (unchecked) `ALLOCATE`
//!   semantics, so exhaustion behavior is identical by construction.
//! - `Model`, `dt`, and the transient flag travel exactly as in the
//!   original: the model address feeds the wrappers that read
//!   `Model % Simulation`/`Solvers`/`NumberOfBodies`; `dt` feeds the DFT
//!   interval math (`ta = time - dt`); the transient flag is unused by
//!   the body, exactly as in the original.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fls_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn fls_init_body();
    fn fls_head();
    fn fls_setup_nodal();
    fn fls_setup_elem();
    fn fls_nodalloss_flag() -> c_int;
    fn fls_perm_empty() -> c_int;
    fn fls_fetch_target();
    fn fls_find_target();
    fn fls_fourier_dofs();
    fn fls_fourier_alloc();
    fn fls_loss_alloc();
    fn fls_omega();
    fn fls_dft(end0: c_int, dt: c_double) -> c_int;
    fn fls_default_init();
    fn fls_otherrhs();
    fn fls_bulk();
    fn fls_nodal_solve();
    fn fls_comm();
    fn fls_store_total();
}

/// `fourierlosssolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `FourierLossSolver_init` (obsolescence guard, loss-count
/// loop, `res:` seeding, exported-variable registration).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "fourierlosssolver_init_"]
pub unsafe extern "C" fn fourierloss_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fls_begin(model, solver);
        fls_init_body();
    }
}

/// `fourierlosssolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`fourierloss_solver_init_entry`]; additionally `dt`
/// must be a valid read (it feeds the DFT interval math).
#[export_name = "fourierlosssolver_"]
pub unsafe extern "C" fn fourierloss_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        fls_begin(model, solver);
        fls_head();
        fls_setup_nodal();
        fls_setup_elem();
        fls_nodalloss_flag();
        if fls_perm_empty() != 0 {
            return;
        }
        fls_fetch_target();
        fls_find_target();
        fls_fourier_dofs();
        fls_fourier_alloc();
        fls_loss_alloc();
        fls_omega();
        let dtv = *dt;
        if fls_dft(0, dtv) == 0 {
            return;
        }
        fls_default_init();
        fls_otherrhs();
        fls_bulk();
        fls_nodal_solve();
        fls_comm();
        fls_store_total();
        fls_dft(1, dtv);
    }
}
