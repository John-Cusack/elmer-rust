//! Port of `fem/src/modules/ReynoldsSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (one
//! `.so` → one Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ReynoldsSolver.F90` — four F77-global
//! entries (`reynoldssolver_`, `reynoldssolver_init_`,
//! `reynoldspostprocess_`, `reynoldspostprocess_init_`). The emitted
//! dynamic symbol set is exactly those four (proven by `nm -D` on the
//! reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries'
//! derived-type dummies (`Model`, `Solver`) cross as opaque addresses.
//! Every component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` phase wrappers in
//! `accessors/reynoldssolver_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset.
//!
//! Split of responsibilities:
//! - Rust owns all entry-level control flow, statement for statement:
//!   the nonlinear-iteration loop with its linesearch redo (`GOTO 100`),
//!   the sensitivity redo (`GOTO 200`), the gap-sensitivity gate, and
//!   the postprocess mode/component loops with their `CYCLE`s.
//! - Fortran owns every Elmer call, every `SUM` reduction, every
//!   formatted `WRITE`, every whole-array update, and the transcribed
//!   element kernels (`GlobalBulkAssembly`, `LocalBulkMatrix`,
//!   `GlobalBoundaryAssemby`, `LocalBoundaryMatrix`, postprocess
//!   `LocalMatrix`), so bit-identity holds by construction wherever
//!   the call sequence matches.
//!
//! Faithfulness notes (mirror ReynoldsSolver.F90 statement for statement):
//! - `SAVE` semantics replicated: all permanent storage (nodal work
//!   arrays sized at first-mesh `MaxElementNodes`, the
//!   `Artificial Compressibility` pressure copy, `AllocationsDone`,
//!   `GotPseudoPressure`) lives in wrapper-side SAVE state with
//!   identical lifetimes; the two wrapper modules mirror the two
//!   entries' host scopes one to one.
//! - Allocation failure calls the same `Fatal(..., 'Memory allocation
//!   error')` via the `STAT=istat` check, instead of aborting.
//! - `TransientSimulation` crosses as a normalized 0/1 `c_int` and is
//!   converted to `LOGICAL` once per frame; `dt` crosses by value.
//! - `Model`, `dt`, `TransientSimulation` uses match the original field
//!   for field (`dt` only in the transient/pseudo terms, the transient
//!   flag gating the mass assembly and `dPdt`).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn rrs_frame_begin(
        h_model: *mut c_void,
        h_solver: *mut c_void,
        dt: c_double,
        transient: c_int,
    );
    fn rrs_frame_end();
    fn rrs_preamble(noiter: *mut c_int, proceed: *mut c_int);
    fn rrs_allocate(skipfirst: *mut c_int);
    fn rrs_frame_setup();
    fn rrs_iter_head(iter: c_int);
    fn rrs_initialize();
    fn rrs_bulk(sens: c_int);
    fn rrs_boundary();
    fn rrs_finish_dirichlet();
    fn rrs_linesearch(redo: *mut c_int, converged: *mut c_int);
    fn rrs_solve();
    fn rrs_nonlin(nc: *mut c_int);
    fn rrs_sensitivity(redo: *mut c_int);
    fn rrs_gap_gate(doit: *mut c_int);
    fn rrs_gap_sensitivity();
    fn rrs_finish();
    fn rrs_init_wrap();
    fn rpp_frame_begin(
        h_model: *mut c_void,
        h_solver: *mut c_void,
        dt: c_double,
        transient: c_int,
    );
    fn rpp_frame_end();
    fn rpp_preamble(proceed: *mut c_int);
    fn rpp_mode_open(mode: c_int, comps: *mut c_int, active: *mut c_int);
    fn rpp_assemble(mode: c_int, comp: c_int);
    fn rpp_solve_store();
    fn rpp_mode_tail();
    fn rpp_finish();
    fn rpp_init_wrap();
}

static GUARD: Mutex<()> = Mutex::new(());

/// `reynoldssolver_` — the Reynolds equation solver entry. Argument order
/// and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads.
#[export_name = "reynoldssolver_"]
pub unsafe extern "C" fn reynolds_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    let dtv: c_double = unsafe { *dt };
    let trv: c_int = unsafe { *transient };
    unsafe {
        rrs_frame_begin(model, solver, dtv, trv);
        // Header, DefaultStart, parameter fetch, early-return gates.
        let (mut noiter, mut proceed) = (0, 0);
        rrs_preamble(&mut noiter, &mut proceed);
        // First-mesh permanent storage (+ Skip First Solution gate).
        let mut skip = 0;
        rrs_allocate(&mut skip);
        if skip != 0 {
            rrs_frame_end();
            return;
        }
        rrs_frame_setup();
        // Nonlinear iteration loop (label 200 redo on DefaultSensitivity).
        let mut niter = noiter;
        'sens: loop {
            'outer: for it in 1..=niter {
                rrs_iter_head(it);
                // Assembly redo loop (label 100 on DefaultLinesearch).
                loop {
                    rrs_initialize();
                    rrs_bulk(0);
                    rrs_boundary();
                    rrs_finish_dirichlet();
                    let (mut redo, mut conv) = (0, 0);
                    rrs_linesearch(&mut redo, &mut conv);
                    if redo != 0 {
                        continue;
                    }
                    if conv != 0 {
                        break 'outer;
                    }
                    rrs_solve();
                    let mut nc = 0;
                    rrs_nonlin(&mut nc);
                    if nc > 0 {
                        break 'outer;
                    }
                    break;
                }
            }
            let mut redo = 0;
            rrs_sensitivity(&mut redo);
            if redo == 0 {
                break 'sens;
            }
            niter = 1;
        }
        // Gap-sensitivity block.
        let mut doit = 0;
        rrs_gap_gate(&mut doit);
        if doit != 0 {
            rrs_gap_sensitivity();
        }
        rrs_finish();
        rrs_frame_end();
    }
}

/// `reynoldssolver_init_` — initialization of the Reynolds solver entry
/// (pure list manipulation, kept verbatim in the wrapper).
///
/// # Safety
///
/// Called by the loader with the entry's usual four arguments.
#[export_name = "reynoldssolver_init_"]
pub unsafe extern "C" fn reynolds_solver_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        rrs_init_wrap();
    }
}

/// `reynoldspostprocess_` — postprocessing fields (force, flux, mean
/// velocity, heating) from the film pressure. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads.
#[export_name = "reynoldspostprocess_"]
pub unsafe extern "C" fn reynolds_postprocess_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    let dtv: c_double = unsafe { *dt };
    let trv: c_int = unsafe { *transient };
    unsafe {
        rpp_frame_begin(model, solver, dtv, trv);
        let mut proceed = 0;
        rpp_preamble(&mut proceed);
        if proceed == 0 {
            rpp_frame_end();
            return;
        }
        for mode in 0..=4 {
            let (mut comps, mut active) = (0, 0);
            rpp_mode_open(mode, &mut comps, &mut active);
            if active == 0 {
                continue;
            }
            for comp in 1..=comps {
                rpp_assemble(mode, comp);
                if mode == 0 {
                    continue;
                }
                rpp_solve_store();
            }
            rpp_mode_tail();
        }
        rpp_finish();
        rpp_frame_end();
    }
}

/// `reynoldspostprocess_init_` — initialization of the postprocess entry
/// (pure list manipulation, kept verbatim in the wrapper).
///
/// # Safety
///
/// Called by the loader with the entry's usual four arguments.
#[export_name = "reynoldspostprocess_init_"]
pub unsafe extern "C" fn reynolds_postprocess_init_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe {
        rpp_init_wrap();
    }
}
