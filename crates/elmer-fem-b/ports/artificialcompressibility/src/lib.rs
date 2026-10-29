//! Port of `fem/src/modules/ArtificialCompressibility.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` -> Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ArtificialCompressibility.F90` -- three
//! F77-global entries (`compressibilityscale_`,
//! `compressibilityscale_init_`, `compressibilitysolver_`) with three
//! `CONTAINS` internals (`PressureIntegrate`, `CompressibilityIntegrate`,
//! `LocalMatrix`). The emitted dynamic symbols are exactly those three
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 section 3, solver-entry gate shape): the entries'
//! derived dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/artificialcompressibility_acc.F90`, compiled against the
//! built `.mod` files of the target configure -- the port itself contains
//! no layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MINVAL`/`MAXVAL`), formatted output (`WRITE`+`Info`, `PRINT`),
//! and whole-array updates stay in Fortran wrappers, so bit-identity holds
//! by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ArtificialCompressibility.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: CompressibilityScale's per-call arrays
//!   and accumulators are allocated/zeroed/freed at the same statements;
//!   its SAVE scalars (`TimeStepVisited`, `PrevTimestep`, `InitVolume`,
//!   `MinimumSideVolume`) live wrapper-side. (`SubroutineVisited` is
//!   set-but-never-read and not replicated.) CompressibilitySolver's SAVE
//!   arrays live wrapper-side behind `AllocationsDone` with first-mesh
//!   sizing (`MaxElementNodes`), exactly as in the original.
//! - Allocation failure calls the same `Fatal(...)` with the same caller
//!   tag (including the original's `'CompressibilitySolve'` spelling in
//!   the solver entry) instead of aborting.
//! - `Model % CurrentElement` pointer assignments, the `FVar`/`DVar`
//!   association fallbacks, and the `Variable % Norm` store stay in
//!   wrappers. The `.sif` default values for absent keywords (1.0, 0.01,
//!   `.TRUE.`) are applied on the Rust side, matching the original's
//!   `IF(.NOT.GotIt)` branches.
//! - `CompressibilityScale_init` only registers the three `res:`
//!   Simulation entries, exactly as in the original.
//! - `PressureExists` evident-intent fix: the original associates the
//!   pressure pointers and then immediately assigns `PressureExists =
//!   .FALSE.`, so the element loop could only ever read `ReferencePressure`
//!   (uninitialized when a pressure variable exists -- undefined
//!   behavior). The port uses the pressure values when a pressure
//!   variable is present. The proof exercises the defined Reference
//!   Pressure path byte-identically; see the proof header.
//! - `dt` is unused by all three entries, exactly as in the original;
//!   `Transient` is read only by `CompressibilitySolver` for its
//!   steady-state guard. No Rust threads: the solver entries run serially.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    // -- CompressibilityScale --
    fn acs_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn acs_info_start();
    fn acs_timestep();
    fn acs_params_begin();
    fn acs_nsolvers() -> c_int;
    fn acs_flow_is_ns(i: c_int) -> c_int;
    fn acs_flow_select(i: c_int);
    fn acs_fvar_assoc() -> c_int;
    fn acs_flow_fallback();
    fn acs_flow_fatal();
    fn acs_pressure_mode() -> c_int;
    fn acs_disp_lookup();
    fn acs_dvar_assoc() -> c_int;
    fn acs_disp_fatal();
    fn acs_maxnodes() -> c_int;
    fn acs_alloc(n: c_int);
    fn acs_sweep_begin();
    fn acs_nbulk() -> c_int;
    fn acs_nbnd() -> c_int;
    fn acs_bulk_begin(t: c_int, n: *mut c_int, skip: *mut c_int);
    fn acs_bulk_mat(got: *mut c_int);
    fn acs_bnd_begin(t: c_int, n: *mut c_int, skip: *mut c_int);
    fn acs_compress_integrate();
    fn acs_pressure_integrate();
    fn acs_side_update();
    fn acs_sim_scaling(val: *mut c_double, found: *mut c_int);
    fn acs_relax(val: *mut c_double, found: *mut c_int);
    fn acs_transvol(val: *mut c_double, found: *mut c_int);
    fn acs_doscale(val: *mut c_int, found: *mut c_int);
    fn acs_tail_compute(oldscale: c_double, relax: c_double, transvol: c_double, doscale: c_int);
    fn acs_free();
    // -- CompressibilityScale_init --
    fn acs_init_begin(h_model: *mut c_void);
    fn acs_init_body();
    // -- CompressibilitySolver --
    fn acb_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn acb_info_start();
    fn acb_has_matrix() -> c_int;
    fn acb_transient_warn();
    fn acb_dim_begin();
    fn acb_disp_lookup() -> c_int;
    fn acb_disp_fatal();
    fn acb_pres_lookup() -> c_int;
    fn acb_pres_fatal();
    fn acb_displaced_read();
    fn acb_alloc_n() -> c_int;
    fn acb_alloc_if_needed(n: c_int);
    fn acb_zero_elem();
    fn acb_default_init();
    fn acb_nelems() -> c_int;
    fn acb_elem_begin(t: c_int, n: *mut c_int);
    fn acb_local_matrix();
    fn acb_update();
    fn acb_finish();
    fn acb_solve() -> c_double;
    fn acb_elim_neg();
}

/// `compressibilityscale_` -- the FSI scaling entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "compressibilityscale_"]
pub unsafe extern "C" fn compressibility_scale_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        acs_begin(model, solver);
        acs_info_start();
        acs_timestep();
        acs_params_begin();
        // Flow-solver scan: first solver whose Equation is 'navier-stokes'.
        let nsolvers = acs_nsolvers();
        let mut i: c_int = 1;
        while i <= nsolvers {
            if acs_flow_is_ns(i) != 0 {
                acs_flow_select(i);
                break;
            }
            i += 1;
        }
        if acs_fvar_assoc() == 0 {
            acs_flow_fallback();
        }
        if acs_fvar_assoc() == 0 {
            acs_flow_fatal();
        }
        let pressure_mode = acs_pressure_mode();
        acs_disp_lookup();
        if acs_dvar_assoc() == 0 {
            acs_disp_fatal();
        }
        acs_alloc(acs_maxnodes());
        acs_sweep_begin();
        if pressure_mode == 0 {
            let nbulk = acs_nbulk();
            let mut t: c_int = 1;
            while t <= nbulk {
                let mut n: c_int = 0;
                let mut skip: c_int = 0;
                acs_bulk_begin(t, &mut n, &mut skip);
                if skip == 0 {
                    let mut got: c_int = 0;
                    acs_bulk_mat(&mut got);
                    if got != 0 {
                        acs_compress_integrate();
                    }
                }
                t += 1;
            }
        }
        // Force acting on the boundary (shared by both solvers).
        let nbulk = acs_nbulk();
        let nbnd = acs_nbnd();
        let mut t: c_int = nbulk + 1;
        while t <= nbulk + nbnd {
            let mut n: c_int = 0;
            let mut skip: c_int = 0;
            acs_bnd_begin(t, &mut n, &mut skip);
            if skip == 0 {
                acs_pressure_integrate();
            }
            t += 1;
        }
        acs_side_update();
        let mut oldscale: c_double = 0.0;
        let mut found: c_int = 0;
        acs_sim_scaling(&mut oldscale, &mut found);
        if found == 0 {
            oldscale = 1.0;
        }
        let mut relax: c_double = 0.0;
        let mut rfound: c_int = 0;
        acs_relax(&mut relax, &mut rfound);
        if rfound == 0 {
            relax = 1.0;
        }
        let mut transvol: c_double = 0.0;
        let mut tfound: c_int = 0;
        acs_transvol(&mut transvol, &mut tfound);
        if tfound == 0 {
            transvol = 0.01;
        }
        let mut doscale: c_int = 0;
        let mut sfound: c_int = 0;
        acs_doscale(&mut doscale, &mut sfound);
        if sfound == 0 {
            doscale = 1;
        }
        acs_tail_compute(oldscale, relax, transvol, doscale);
        acs_free();
    }
}

/// `compressibilityscale_init_` -- the `_Init` entry Elmer calls once at
/// solver setup. Only registers the three `res:` Simulation entries.
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the loader passes. `dt`
/// and the transient flag are unused, exactly as in the original.
#[export_name = "compressibilityscale_init_"]
pub unsafe extern "C" fn compressibility_scale_init_entry(
    model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        acs_init_begin(model);
        acs_init_body();
    }
}

/// `compressibilitysolver_` -- the elementwise compressibility entry.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` is unused, exactly as
/// in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `transient` must be a valid read if non-null.
#[export_name = "compressibilitysolver_"]
pub unsafe extern "C" fn compressibility_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    unsafe {
        acb_begin(model, solver);
        acb_info_start();
        if acb_has_matrix() == 0 {
            return;
        }
        if !transient.is_null() && *transient != 0 {
            acb_transient_warn();
            return;
        }
        acb_dim_begin();
        if acb_disp_lookup() == 0 {
            acb_disp_fatal();
        }
        match acb_pres_lookup() {
            0 | 1 => {}
            _ => acb_pres_fatal(),
        }
        acb_displaced_read();
        acb_alloc_if_needed(acb_alloc_n());
        acb_zero_elem();
        acb_default_init();
        let nelems = acb_nelems();
        let mut t: c_int = 1;
        while t <= nelems {
            let mut n: c_int = 0;
            acb_elem_begin(t, &mut n);
            acb_local_matrix();
            acb_update();
            t += 1;
        }
        acb_finish();
        let _norm: c_double = acb_solve();
        acb_elim_neg();
    }
}
