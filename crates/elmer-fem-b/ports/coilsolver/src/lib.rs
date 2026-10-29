//! Port of `fem/src/modules/CoilSolver.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/CoilSolver.F90` — four F77-global
//! entries (`coilsolver_`, `coilsolver_init_`, `coilpotential_`,
//! `coilpotentialnormalized_`) with fifteen `CONTAINS` internals. The
//! emitted dynamic symbols are exactly those four (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/coil_acc.F90`, compiled against the built `.mod` files of
//! the target configure — the port itself contains no layout assumption
//! and no pinned offset. Order-sensitive reductions (`SUM`, `MINVAL`,
//! `MAXVAL`, `MATMUL`), formatted I/O (`WRITE`/`PRINT`), the `DSYEV`
//! eigen-fit, `LUdecomp`/`LUSolve`, and whole-array updates stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors CoilSolver.F90 statement for statement):
//! - The Rust entries replicate the entry-level control flow: the
//!   `coilsolver_init` list setup (kept verbatim in `coa_init_body`,
//!   mirroring the DivergenceSolver `_init` precedent), the `CoilSolver`
//!   guard/header/setup/flag reads, the per-coil `DO` loop with its
//!   `CYCLE`/`EXIT`/tail structure, the `GOTO 100` (`TestCut`) bypass,
//!   the `MaxNonlinIter`/`Part`/`dimi` loops, and the `CalcCurr` and
//!   `Fix Input Current Density` branches. All data-touching work
//!   (element sweeps, `SetA`/`SetB` fixing sets, `UpdateDirichletDof`,
//!   `CutInterfaceConnections`, `ScalePotential`,
//!   `NormalizeCurrentDensity`) runs in the corresponding `coa_*`
//!   drivers, which are verbatim transcriptions of the original
//!   statements (same call sequence, same compiler).
//! - The fifteen `CONTAINS` internals are transcribed verbatim as module
//!   procedures of `coil_acc` (host association replaced by
//!   module-level `SAVE` state with identical lifetimes; the
//!   host-associated scratch `j` of `LocalCorrMatrix` and the `Found`
//!   scratch of `MarkCoilNodes`/`ChooseCoilCut` — implicitly host-shared
//!   in the original — are procedure locals here, and the clobbering is
//!   unobservable since every use is assigned before it is read).
//! - `SAVE` semantics replicated explicitly: the per-internal-procedure
//!   `SAVE Nodes` mirrors live in wrapper-side `SAVE` state (consumed
//!   before the next fill, as in the original); the entry itself keeps
//!   no `SAVE` state (the original allocates/frees its work arrays per
//!   call), so in-process re-entry (steady iterations) matches by
//!   construction. The `Visited` caches of `CoilPotential` and
//!   `CoilPotentialNormalized` live wrapper-side with identical
//!   lifetimes, including the first-visit overwrite of the caller's `n`
//!   with `MaxElementNodes` in the normalized variant.
//! - The `ParEnv % PEs > 1` branches (`MPI_ALLREDUCE`, `ParallelReduction`,
//!   `ParallelSumVectorInt`) are transcribed verbatim and guard on the
//!   same runtime condition. The serial proof path never executes them;
//!   the parallel path is preserved for completeness but not proven.
//!   (The undefined `mpi_allreduce_` symbol resolves at load time from
//!   the already-loaded `libmpi_stubs`, exactly as for the reference.)
//! - The open-coil fallback path (no `Component` sections) is the proven
//!   path; the closed-coil fitting paths (`DefineCoilCenter`,
//!   `DefineCoilParameters`, `ChooseFixedBulkNodesNarrow`,
//!   `ChooseCoilCut`) are transcribed verbatim but unproven by the
//!   current differentials (no closed-coil mesh is staged).
//! - `Model`, `dt`, and the transient flag are unused by the solver body
//!   beyond handle duties, exactly as in the original (received and
//!   ignored); `Norm = DefaultSolve()` discards are kept.
//! - No Rust threads: the solver entries run serially, and all shared
//!   state lives in the Fortran wrapper module (same single-threaded
//!   lifetime discipline as the original).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn coa_init_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn coa_init_body();
    fn coa_matrix_associated(h_model: *mut c_void, h_solver: *mut c_void) -> c_int;
    fn coa_header();
    fn coa_setup();
    fn coa_read_flags();
    fn coa_parts() -> c_int;
    fn coa_calc_curr() -> c_int;
    fn coa_maxnonlin() -> c_int;
    fn coa_dim() -> c_int;
    fn coa_num_components() -> c_int;
    fn coa_test_cut() -> c_int;
    fn coa_pots_setup();
    fn coa_coil_try(i: c_int) -> c_int;
    fn coa_coil_tail();
    fn coa_count_fixing();
    fn coa_set_iter(it: c_int);
    fn coa_fix_stage();
    fn coa_pot_header();
    fn coa_pot_part(part: c_int);
    fn coa_curr_begin();
    fn coa_flux_comp(dimi: c_int);
    fn coa_curr_part_tail();
    fn coa_curr_end();
    fn coa_fixinput() -> c_int;
    fn coa_fixinput_stage();
    fn coa_tail_100();
    fn coa_cpot_eval(h_model: *mut c_void, n: c_int) -> c_double;
    fn coa_cpotn_eval(h_model: *mut c_void, n: *mut c_int) -> c_double;
}

/// `coilsolver_init_` — the `_Init` entry Elmer calls once at solver
/// setup. Mirrors `CoilSolver_init` (exported-variable registration,
/// `Coil Start`/`Coil End` consistency fixups on `Electrode Boundaries`
/// components).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "coilsolver_init_"]
pub unsafe extern "C" fn coilsolver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        coa_init_begin(model, solver);
        coa_init_body();
    }
}

/// `coilsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body beyond handle duties, exactly
/// as in the original.
///
/// # Safety
///
/// Same contract as [`coilsolver_init_entry`].
#[export_name = "coilsolver_"]
pub unsafe extern "C" fn coilsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        // IF(.NOT. ASSOCIATED(Solver % Matrix)) RETURN
        if coa_matrix_associated(model, solver) == 0 {
            return;
        }
        coa_header();
        coa_setup();
        coa_read_flags();
        let parts = coa_parts();
        let calc_curr = coa_calc_curr();
        let maxnonlin = coa_maxnonlin();
        coa_pots_setup();

        // DO i=1,Model % NumberOfComponents + 1
        let ncomp = coa_num_components();
        let mut i = 1;
        while i <= ncomp + 1 {
            // 0 = CYCLE (not a coil), 1 = PROCEED (tail), 2 = EXIT loop.
            match coa_coil_try(i) {
                1 => coa_coil_tail(),
                2 => break,
                _ => {}
            }
            i += 1;
        }
        coa_count_fixing();

        // IF(TestCut) GOTO 100
        if coa_test_cut() == 0 {
            let mut iter = 1;
            while iter <= maxnonlin {
                coa_set_iter(iter);
                if iter > 1 {
                    coa_fix_stage();
                }
                coa_pot_header();
                let mut part = 1;
                while part <= parts {
                    coa_pot_part(part);
                    part += 1;
                }
                iter += 1;
            }

            if calc_curr != 0 {
                coa_curr_begin();
                let dim = coa_dim();
                let mut part = 1;
                while part <= parts {
                    let mut dimi = 1;
                    while dimi <= dim {
                        coa_flux_comp(dimi);
                        dimi += 1;
                    }
                    coa_curr_part_tail();
                    part += 1;
                }
                coa_curr_end();
            }

            if coa_fixinput() != 0 {
                coa_fixinput_stage();
            }
        }

        // 100 CONTINUE tail
        coa_tail_100();
    }
}

/// `coilpotential_` — the closed-coil potential user function
/// (`CoilPotential(Model, n, t)`). `n` is the node number (read only);
/// `t` is unused, exactly as in the original.
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the caller passes;
/// `n`/`t` must be valid reads if non-null.
#[export_name = "coilpotential_"]
pub unsafe extern "C" fn coilpotential_entry(
    model: *mut c_void,
    n: *const c_int,
    _t: *const c_double,
) -> c_double {
    unsafe { coa_cpot_eval(model, *n) }
}

/// `coilpotentialnormalized_` — the normalized-potential user function
/// (`CoilPotentialNormalized(Model, n, t)`). On first visit the caller's
/// `n` is overwritten with `MaxElementNodes` (verbatim side effect of
/// the original).
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the caller passes; `n`
/// must be a valid read/write if non-null; `t` must be a valid read if
/// non-null (it is ignored).
#[export_name = "coilpotentialnormalized_"]
pub unsafe extern "C" fn coilpotentialnormalized_entry(
    model: *mut c_void,
    n: *mut c_int,
    _t: *const c_double,
) -> c_double {
    unsafe { coa_cpotn_eval(model, n) }
}
