//! Port of `fem/src/modules/CoordinateTransform.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/CoordinateTransform.F90` — two F77-global
//! entries (`rotmsolver_`, `rotmsolver_init_`) plus the `VecUtils` module's
//! four procedures (`__vecutils_MOD_jac`, `__vecutils_MOD_norm2`,
//! `__vecutils_MOD_normalized`, `__vecutils_MOD_transform2rank`). The
//! emitted dynamic symbols are exactly those six (proven by `nm -D` on the
//! reference `.so`); the `CONTAINS` internals (`ComputeRotM`,
//! `PolarDecomposition`, `GetElementRotM`) emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/coordinatetransform_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the `MATMUL` reductions, `SUM`/`DOT_PRODUCT`,
//! `SQRT`, the polar decomposition, the formatted `PRINT`/`Warn` reports)
//! stays in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches. The four `VecUtils` procedures are
//! transcribed statement-for-statement into Fortran
//! (`accessors/vecutils_rep.F90`, same explicit-shape `REAL(8)` dummies —
//! the plain-ABI triage outcome), so their exported symbols keep the exact
//! compiler ABI and evaluation order by construction.
//!
//! Faithfulness notes (mirrors CoordinateTransform.F90 statement for
//! statement):
//!   backing the original's `RotM(3,3,N)` allocation) happens once behind
//!   the `First` flag; the internal procedures' `SAVE Nodes`, the `SAVE`
//!   variable pointers, the `SAVE` polar-decomposition settings, and the
//!   `SAVE`d `alpha_ref`/`beta_ref` pointers live as wrapper-side SAVE
//!   state with identical lifetimes (consumed before the next fill, as in
//!   the original). `Mesh` (not `SAVE` in the original) is kept with the
//!   author's evident intent: the same mesh address persists across steady
//!   iterations.
//! - Allocation failure calls the same `Fatal('RotMSolver', ...)` via
//!   `try_reserve` (the `STAT=istat` check), never aborting directly.
//! - `Model`, `dt`, `TransientSimulation` are unused by the entry bodies
//!   except for the first-visit mesh sizing, exactly as in the original.
//! - `RotMSolver_init` is an all-comments body in the original; the entry
//!   is a matching no-op.
//! - `GetElementRotM` (a `CONTAINS` internal) is never called by any path
//!   in the original, so it is documented as skipped rather than
//!   replicated; it emits no symbol either way.
//! - The unfilled-`alpha`/`beta` path (when neither `alpha direction` nor
//!   `alpha` present) reads undefined automatic-array garbage in the
//!   original; the port zero-fills those buffers instead (deterministic,
//!   documented). The proof configs always provide the variables, so the
//!   defined path is byte-identical.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn cta_first_n(h_model: *mut c_void) -> c_int;
    fn cta_first_vars();
    fn cta_nactive() -> c_int;
    fn cta_elem_begin(q: c_int, nn: *mut c_int) -> c_int;
    fn cta_node(j: c_int, nn: c_int);
    fn cta_fatal_alloc();
}

/// First-mesh SAVE state: the `RotM(3,3,N)` backing store (`N` is the first
/// call's `MaxElementDOFs`; the original allocates it and never reads it —
/// only the allocation-failure `Fatal` is observable) plus the `First`
/// flag.
struct CtState {
    rotm: Vec<c_double>,
    done: bool,
}

impl CtState {
    const fn new() -> Self {
        Self {
            rotm: Vec::new(),
            done: false,
        }
    }
}

static STATE: Mutex<CtState> = Mutex::new(CtState::new());

/// `rotmsolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `RotMSolver_init`, whose body is entirely comments.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "rotmsolver_init_"]
pub unsafe extern "C" fn rotmsolver_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    // No-op by transcription: the original body holds only comments.
}

/// `rotmsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Solver`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// Same contract as [`rotmsolver_init_entry`].
#[export_name = "rotmsolver_"]
pub unsafe extern "C" fn rotmsolver_entry(
    model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        let mut state = STATE.lock().unwrap();
        if !state.done {
            let n = cta_first_n(model);
            if state.rotm.try_reserve(9 * n as usize).is_err() {
                cta_fatal_alloc();
            }
            state.rotm.resize(9 * n as usize, 0.0);
            cta_first_vars();
            state.done = true;
        }
        let nactive = cta_nactive();
        for q in 1..=nactive {
            let mut nn: c_int = 0;
            if cta_elem_begin(q, &mut nn) == 1 {
                continue;
            }
            for j in 1..=nn {
                cta_node(j, nn);
            }
        }
    }
}
