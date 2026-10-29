//! Port of `fem/src/modules/CircuitsAndDynamics.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/CircuitsAndDynamics.F90` — six F77-global
//! entries (`circuitsanddynamics_`, `circuitsanddynamics_init_`,
//! `circuitsanddynamicsharmonic_`, `circuitsanddynamicsharmonic_init_`,
//! `circuitsoutput_`, `circuitsoutput_init_`). The emitted dynamic symbols
//! are exactly those six (proven by `nm -D` on the reference `.so`); the
//! lifted kernels and shims emit only local symbols (version-scripted).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/circuitsanddynamics_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. The
//! element assembly kernels (`Add_stranded`, `Add_massive`,
//! `Add_foil_winding`, both drivers), the basic-equation writers, the
//! per-circuit reporter, and all formatted I/O stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches. The `SAVE` state of the internal procedures (per-kernel
//! `Nodes`, `Wvec_h`/`Jvec_h`, `Wpot`, generation caches, the EEC block,
//! `CktPrefix`) lives wrapper-side with identical lifetimes.
//!
//! Faithfulness notes (mirrors CircuitsAndDynamics.F90 entry for entry):
//! - The `_init` entries claim the container, bind the `Pre Solvers` slot,
//!   and install the same `Variable` / `Variable Output` /
//!   `Time Derivative Order` / `No Matrix` defaults; the rotor-angle
//!   export block (transient `_init` only) is replicated branch for
//!   branch, including the `Rotational Projector` fallback.
//! - `Solver % Values => Params` is a no-op (the pointer is unchanged) and
//!   is documented as skipped; the unused `Mesh` pointer of the drivers is
//!   likewise skipped.
//! - The build block runs once per container (`BuiltNm < 0`, re-armed by
//!   `CircuitsCheckStale`), with the same per-circuit library call order,
//!   the same parallel-mismatch warn-and-fix, and the same post-build
//!   check sequence. `Circuits => ... % Circuits` re-pointers are
//!   wrapper-internal.
//! - The rotor-dynamics call, the timestep-refresh block (including the
//!   `LagrangeVar` association / `PrevValues` rotation / `Crt` reload and
//!   the unconditional `Circuits_ToMeshVariable`), and the assembly tail
//!   (zero, per-circuit equations, `AddMatrix` link, `MATRIX_LIST`
//!   conversion, `InfoActive(20)` ranges, `Save Circuit Matrix`,
//!   empty-matrix teardown) are replicated in order, including the
//!   second `ASSOCIATED(CM)` test of the transient driver.
//! - The harmonic driver additionally issues `DefaultStart` /
//!   `DefaultFinish`, reads `Eigen Analysis` into the module flag the
//!   `AddMatrixEntry` kernel branches on, and zeroes `MassValues`.
//! - The output driver resolves (not claims) the model, runs the EEC
//!   steady-state forcing block, gathers the Lagrange values (serial
//!   path; the `MPI_ALLREDUCE` arm is compiled but only executes when
//!   `Solver % Parallel` is set, as in the original), exports the
//!   `time` stamp, and reports each circuit in order before
//!   `Circuits_ToMeshVariable` and `DefaultFinish`.
//! - No Rust threads: the solver entries run serially, and the port keeps
//!   no Rust-side state at all, so no `static` guard is needed
//!   (`unsafe` only at the documented handle boundary).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn cda_begin(h_model: *mut c_void, h_solver: *mut c_void, dt: c_double, transient: c_int);
    fn cda_claim_only();
    fn cda_claim_set();
    fn cda_solver_bind(which: c_int);
    fn cda_add_var(which: c_int);
    fn cda_add_var_output();
    fn cda_add_time_order();
    fn cda_add_no_matrix();
    fn cda_sim_has_rotor_angle() -> c_int;
    fn cda_rot_machine() -> c_int;
    fn cda_add_exported(which: c_int);
    fn cda_check_stale();
    fn cda_built_nm() -> c_int;
    fn cda_info_init(which: c_int);
    fn cda_set_harmonic(v: c_int);
    fn cda_add_comp_to_bodies();
    fn cda_circ_reset();
    fn cda_n_circuits() -> c_int;
    fn cda_hm_cm_ready() -> c_int;
    fn cda_find_asolver();
    fn cda_parallel_mismatch() -> c_int;
    fn cda_solver_parallel() -> c_int;
    fn cda_warn_parallel_fix();
    fn cda_set_ckt_parallel();
    fn cda_check_multiplier_unique();
    fn cda_set_matc_prefix();
    fn cda_allocate_circuits_list();
    fn cda_set_boundary_areas();
    fn cda_info_circuit(p: c_int, n: c_int);
    fn cda_n_circ_vars(p: c_int) -> c_int;
    fn cda_allocate_circuit(p: c_int);
    fn cda_count_comps(p: c_int, n: c_int) -> c_int;
    fn cda_alloc_comps(p: c_int, ncomp: c_int);
    fn cda_setup_circuit(p: c_int, harm: c_int, par: c_int);
    fn cda_read_circuit_vars(p: c_int);
    fn cda_read_components(p: c_int);
    fn cda_add_comp_values(p: c_int);
    fn cda_add_bare_vars(p: c_int);
    fn cda_read_coef_matrices(p: c_int);
    fn cda_read_perm_vector(p: c_int);
    fn cda_read_sources(p: c_int);
    fn cda_write_coeff_vectors(p: c_int);
    fn cda_check_component_variables();
    fn cda_check_circuit_sources();
    fn cda_build_elem_lists();
    fn cda_circuits_summary();
    fn cda_check_transient_components();
    fn cda_circuits_matrix_init();
    fn cda_realloc_crt();
    fn cda_info_active(level: c_int) -> c_int;
    fn cda_set_dynamic_angle();
    fn cda_tstep_new() -> c_int;
    fn cda_find_lagrange() -> c_int;
    fn cda_lag_size() -> c_int;
    fn cda_lag_range();
    fn cda_circuit_tot_n() -> c_int;
    fn cda_fatal_lagrange_small();
    fn cda_lag_is_solver_var() -> c_int;
    fn cda_info_assoc_lag();
    fn cda_assoc_solver_var();
    fn cda_lag_has_prev() -> c_int;
    fn cda_info_prev();
    fn cda_alloc_prev();
    fn cda_rotate_lagrange();
    fn cda_load_crt();
    fn cda_to_mesh_var();
    fn cda_max_elem_dofs() -> c_int;
    fn cda_cm_assoc() -> c_int;
    fn cda_td_zero_cm();
    fn cda_td_basic(p: c_int);
    fn cda_td_comp(p: c_int, nn: c_int);
    fn cda_link_addmatrix();
    fn cda_unlink_addmatrix();
    fn cda_cm_range();
    fn cda_save_circuit_matrix();
    fn cda_free_if_empty() -> c_int;
    fn cda_info_finished();
    fn cda_hm_start();
    fn cda_read_eigen();
    fn cda_hm_zero_cm();
    fn cda_hm_basic(p: c_int);
    fn cda_hm_comp(p: c_int, nn: c_int);
    fn cda_hm_finish();
    fn cda_out_start();
    fn cda_out_claim();
    fn cda_out_asolver_assoc() -> c_int;
    fn cda_fatal_no_asolver();
    fn cda_out_eec();
    fn cda_out_gather();
    fn cda_out_time();
    fn cda_out_ncirc() -> c_int;
    fn cda_out_circuit(p: c_int);
    fn cda_out_mesh_var();
    fn cda_out_finish();
}

/// Read the by-reference `dt`/`transient` dummies once, as plain scalars.
unsafe fn entry_scalars(dt: *const c_double, transient: *const c_int) -> (c_double, c_int) {
    let dtv = if dt.is_null() { 0.0 } else { *dt };
    let tr = if transient.is_null() { 0 } else { *transient };
    (dtv, tr)
}

/// The transient/harmonic build block, shared up to the branches the two
/// drivers differ on. `which` selects the driver: 0 transient, 1 harmonic.
unsafe fn build_circuits(which: c_int, par: c_int) {
    cda_info_init(which);
    cda_set_harmonic(which);
    cda_add_comp_to_bodies();
    // The count is populated from the MATC definitions by
    // AllocateCircuitsList(), so it is read after that call — the
    // original holds it through a pointer into the container.
    cda_circ_reset();
    cda_find_asolver();
    if cda_parallel_mismatch() != 0 {
        cda_warn_parallel_fix();
    }
    cda_set_ckt_parallel();
    cda_check_multiplier_unique();
    cda_set_matc_prefix();
    cda_allocate_circuits_list();
    cda_set_boundary_areas();
    let n = cda_n_circuits();
    let mut p: c_int = 1;
    while p <= n {
        let nv = cda_n_circ_vars(p);
        if which == 0 {
            cda_info_circuit(p, nv);
        }
        cda_allocate_circuit(p);
        let nc = cda_count_comps(p, nv);
        cda_alloc_comps(p, nc);
        cda_setup_circuit(p, which, par);
        cda_read_circuit_vars(p);
        cda_read_components(p);
        cda_add_comp_values(p);
        cda_add_bare_vars(p);
        cda_read_coef_matrices(p);
        cda_read_perm_vector(p);
        cda_read_sources(p);
        cda_write_coeff_vectors(p);
        p += 1;
    }
    cda_check_component_variables();
    cda_check_circuit_sources();
    cda_build_elem_lists();
    cda_circuits_summary();
    if which == 0 {
        cda_check_transient_components();
    }
    cda_circuits_matrix_init();
    if which == 0 {
        cda_realloc_crt();
    }
}

/// `circuitsanddynamics_init_` — claims the container, binds the `Pre
/// Solvers` slot, and installs the variable defaults plus the optional
/// rotor-angle exports.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads if non-null (values unused here).
#[export_name = "circuitsanddynamics_init_"]
pub unsafe extern "C" fn circuitsanddynamics_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_claim_only();
    cda_solver_bind(0);
    cda_add_var(0);
    cda_add_var_output();
    cda_add_time_order();
    cda_add_no_matrix();
    if cda_sim_has_rotor_angle() == 0 && cda_rot_machine() != 0 {
        cda_add_exported(0);
        if tr != 0 {
            cda_add_exported(1);
        }
    }
}

/// `circuitsanddynamics_` — the transient driver: build-once block,
/// rotor dynamics, timestep refresh of the circuit values, and assembly
/// of the circuit equations into the A solver's matrix.
///
/// # Safety
///
/// Same contract as [`circuitsanddynamics_init_entry`]; `dt` feeds the
/// timestep the rotor integrator and the kernels step with, and
/// `transient` selects the transient branches.
#[export_name = "circuitsanddynamics_"]
pub unsafe extern "C" fn circuitsanddynamics_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_claim_set();
    cda_check_stale();
    if cda_built_nm() < 0 {
        build_circuits(0, cda_solver_parallel());
    }
    if cda_sim_has_rotor_angle() == 0 && tr != 0 {
        cda_set_dynamic_angle();
    }
    if cda_tstep_new() != 0 {
        if cda_find_lagrange() != 0 {
            let n = cda_lag_size();
            if cda_info_active(25) != 0 {
                cda_lag_range();
            }
            if n < cda_circuit_tot_n() {
                cda_fatal_lagrange_small();
            }
            if cda_lag_is_solver_var() == 0 {
                cda_info_assoc_lag();
                cda_assoc_solver_var();
            }
            if cda_lag_has_prev() == 0 {
                cda_info_prev();
                cda_alloc_prev();
            }
            cda_rotate_lagrange();
            cda_load_crt();
        }
        cda_to_mesh_var();
    }
    let maxdofs = cda_max_elem_dofs();
    if cda_cm_assoc() == 0 {
        return;
    }
    cda_td_zero_cm();
    let n = cda_n_circuits();
    let mut p: c_int = 1;
    while p <= n {
        cda_td_basic(p);
        cda_td_comp(p, maxdofs);
        p += 1;
    }
    if cda_cm_assoc() == 0 {
        cda_unlink_addmatrix();
        return;
    }
    cda_link_addmatrix();
    if cda_info_active(20) != 0 {
        cda_cm_range();
    }
    cda_save_circuit_matrix();
    cda_free_if_empty();
    cda_info_finished();
}

/// `circuitsanddynamicsharmonic_init_` — the harmonic `_init` entry:
/// same container claim and slot bind, complex-variable defaults.
///
/// # Safety
///
/// Same contract as [`circuitsanddynamics_init_entry`].
#[export_name = "circuitsanddynamicsharmonic_init_"]
pub unsafe extern "C" fn circuitsanddynamicsharmonic_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_claim_only();
    cda_solver_bind(1);
    cda_add_var(1);
    cda_add_var_output();
    cda_add_no_matrix();
}

/// `circuitsanddynamicsharmonic_` — the harmonic driver: same build
/// discipline as the transient driver, complex assembly, `DefaultStart`
/// / `DefaultFinish` bracketing.
///
/// # Safety
///
/// Same contract as [`circuitsanddynamics_entry`]; the transient flag is
/// received and ignored, exactly as in the original.
#[export_name = "circuitsanddynamicsharmonic_"]
pub unsafe extern "C" fn circuitsanddynamicsharmonic_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_hm_start();
    cda_claim_set();
    cda_check_stale();
    if cda_built_nm() < 0 {
        build_circuits(1, cda_solver_parallel());
    }
    cda_read_eigen();
    let maxdofs = cda_max_elem_dofs();
    if cda_hm_cm_ready() == 0 {
        return;
    }
    cda_hm_zero_cm();
    let n = cda_n_circuits();
    let mut p: c_int = 1;
    while p <= n {
        cda_hm_basic(p);
        cda_hm_comp(p, maxdofs);
        p += 1;
    }
    if cda_cm_assoc() != 0 {
        cda_link_addmatrix();
        cda_free_if_empty();
    } else {
        cda_unlink_addmatrix();
    }
    if cda_info_active(20) != 0 {
        cda_cm_range();
    }
    cda_save_circuit_matrix();
    cda_hm_finish();
    cda_info_finished();
}

/// `circuitsoutput_init_` — binds the `Post Solvers` slot.
///
/// # Safety
///
/// Same contract as [`circuitsanddynamics_init_entry`].
#[export_name = "circuitsoutput_init_"]
pub unsafe extern "C" fn circuitsoutput_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_solver_bind(2);
}

/// `circuitsoutput_` — reports the circuit and dynamic variables for
/// `SaveScalars`: EEC forcing, Lagrange gather, per-circuit reporting.
///
/// # Safety
///
/// Same contract as [`circuitsanddynamics_entry`].
#[export_name = "circuitsoutput_"]
pub unsafe extern "C" fn circuitsoutput_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let (dtv, tr) = entry_scalars(dt, transient);
    cda_begin(model, solver, dtv, tr);
    cda_out_start();
    cda_out_claim();
    if cda_out_asolver_assoc() == 0 {
        cda_fatal_no_asolver();
    }
    cda_out_eec();
    cda_out_gather();
    cda_out_time();
    let n = cda_out_ncirc();
    let mut p: c_int = 1;
    while p <= n {
        cda_out_circuit(p);
        p += 1;
    }
    cda_out_mesh_var();
    cda_out_finish();
}
