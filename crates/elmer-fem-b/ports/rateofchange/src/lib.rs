//! Port of `fem/src/modules/RateOfChange.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/RateOfChange.F90` — one F77-global entry
//! `rateofchangesolver_` with one `CONTAINS` internal (`LocalMatrix`). The
//! emitted dynamic symbol is exactly `rateofchangesolver_` (proven by `nm
//! -D` on the reference `.so`); the internal emits only a local symbol.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/rateofchange_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`), formatted timing I/O
//! (`WRITE(Message,...)`), and whole-array updates stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors RateOfChangeSolver statement for statement):
//! - `SAVE` semantics replicated: the `ThisSolverCalls` counter lives in
//!   the Rust state (the body runs only on odd calls, exactly like the
//!   `MOD(ThisSolverCalls,2)==1` gate); first-odd-call sizing
//!   (`N = Model % MaxElementNodes`, `ld = N*DOFs`) allocates the scratch
//!   buffers once (`LocalStiffMatrix(ld,ld)`, `LocalForce(ld)`,
//!   `LocalU0(ld)`, `V1/V2/V3(N)`), never resized on valid inputs; the
//!   wrapper-side SAVE (`U0Var`, `Nodes`, IPs, equation/advection names)
//!   has identical lifetimes (re-associated every odd call, as in the
//!   original).
//! - Allocation failure calls the same `Fatal('RateOfChangeSolver', ...)`
//!   via `try_reserve` (the `STAT=istat` check), instead of aborting.
//! - `dt`/`TransientSimulation` are unused by the body, exactly as in the
//!   original (received and ignored).
//! - The element loop, `UpdateGlobalEquations`, `SetDirichletBoundaries`,
//!   and `SolveSystem(..., ParMatrix, ...)` run inside the wrappers with
//!   the original argument order; the Rust side only drives the
//!   call sequence.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn roc_info_start();
    fn roc_begin(
        h_model: *const c_void,
        h_solver: *const c_void,
        nmax: *mut c_int,
        dofs: *mut c_int,
        nbulk: *mut c_int,
    );
    fn roc_setup_u0();
    fn roc_setup_vel();
    fn roc_init_zero();
    #[allow(clippy::too_many_arguments)]
    fn roc_assemble_elem(
        t: c_int,
        ld: c_int,
        n: c_int,
        ncomp: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        localu0: *mut c_double,
        v1: *mut c_double,
        v2: *mut c_double,
        v3: *mut c_double,
        skip: *mut c_int,
    );
    fn roc_dirichlet();
    fn roc_solve();
    fn roc_fatal_alloc();
}

/// Grow the SAVE scratch buffers if an element needs more than the
/// first-odd-call sizing (adaptive remesh across invocations). The
/// original never grows (it assumes the first mesh is big enough);
/// growth only triggers where the original would overflow silently, and
/// the leading sections used by the wrappers are unaffected, so valid
/// inputs never observe it.
fn ensure_state(state: &mut RocState, need_ld: usize, need_n: usize) {
    if need_ld > state.ld {
        state.stiff.resize(need_ld * need_ld, 0.0);
        state.force.resize(need_ld, 0.0);
        state.localu0.resize(need_ld, 0.0);
        state.ld = need_ld;
    }
    if need_n > state.n {
        state.v1.resize(need_n, 0.0);
        state.v2.resize(need_n, 0.0);
        state.v3.resize(need_n, 0.0);
        state.n = need_n;
    }
}

/// First-odd-call SAVE state: `ThisSolverCalls`, the `AllocationsDone`
/// flag, sizing `N`/`DOFs`/`ld`, and the scratch buffers.
struct RocState {
    calls: u64,
    n: usize,
    dofs: usize,
    ld: usize,
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    localu0: Vec<c_double>,
    v1: Vec<c_double>,
    v2: Vec<c_double>,
    v3: Vec<c_double>,
    allocations_done: bool,
}

impl RocState {
    const fn new() -> Self {
        Self {
            calls: 0,
            n: 0,
            dofs: 0,
            ld: 0,
            stiff: Vec::new(),
            force: Vec::new(),
            localu0: Vec::new(),
            v1: Vec::new(),
            v2: Vec::new(),
            v3: Vec::new(),
            allocations_done: false,
        }
    }
}

static STATE: Mutex<RocState> = Mutex::new(RocState::new());

/// `rateofchangesolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "rateofchangesolver_"]
pub unsafe extern "C" fn rateofchange_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();

    // ThisSolverCalls = ThisSolverCalls + 1; the body runs only on odd calls.
    state.calls += 1;
    if state.calls.is_multiple_of(2) {
        return;
    }

    unsafe { roc_info_start() };
    let (mut nmax, mut dofs, mut nbulk) = (0 as c_int, 0 as c_int, 0 as c_int);
    unsafe { roc_begin(model, solver, &mut nmax, &mut dofs, &mut nbulk) };

    // First-odd-call permanent storage (bug-compatible: sized once).
    if !state.allocations_done {
        let (n, d) = (nmax as usize, dofs as usize);
        let ld = n * d;
        if state.stiff.try_reserve(ld * ld).is_err()
            || state.force.try_reserve(ld).is_err()
            || state.localu0.try_reserve(ld).is_err()
            || state.v1.try_reserve(n).is_err()
            || state.v2.try_reserve(n).is_err()
            || state.v3.try_reserve(n).is_err()
        {
            unsafe { roc_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(ld * ld, 0.0);
        state.force.resize(ld, 0.0);
        state.localu0.resize(ld, 0.0);
        state.v1.resize(n, 0.0);
        state.v2.resize(n, 0.0);
        state.v3.resize(n, 0.0);
        state.n = n;
        state.dofs = d;
        state.ld = ld;
        state.allocations_done = true;
    }
    // Grow-only guard for remeshed re-entry (invisible otherwise).
    ensure_state(
        &mut state,
        (nmax as usize) * (dofs as usize),
        nmax as usize,
    );
    let (ld, n, ncomp) = (state.ld as c_int, state.n as c_int, dofs);

    unsafe { roc_setup_u0() };
    unsafe { roc_setup_vel() };
    unsafe { roc_init_zero() };

    // Bulk assembly.
    for t in 1..=nbulk {
        let mut skip = 0;
        unsafe {
            roc_assemble_elem(
                t,
                ld,
                n,
                ncomp,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                state.localu0.as_mut_ptr(),
                state.v1.as_mut_ptr(),
                state.v2.as_mut_ptr(),
                state.v3.as_mut_ptr(),
                &mut skip,
            )
        };
        let _ = skip;
    }

    // Dirichlet boundary conditions over components.
    unsafe { roc_dirichlet() };

    // And finally, solve (Norm is wrapper-side state, as in the original).
    unsafe { roc_solve() };
}
