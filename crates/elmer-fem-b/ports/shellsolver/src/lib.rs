//! Port of `fem/src/modules/ShellSolver.F90` (the Elmer fork) to Rust for
//! lane T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib
//! with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ShellSolver.F90` — two F77-global entries
//! `shellsolver_` and `shellsolver_init0_`, both with the standard
//! `(Model, Solver, dt, TransientSimulation)` shape. Derived-type dummies
//! (`Model`, `Solver`) cross as opaque addresses; every component access
//! and every Elmer library call goes through the generated-style Fortran
//! `BIND(C)` wrappers in `accessors/shellsolver_acc.F90`, compiled against
//! the built `.mod` files of the target configure. The port itself
//! contains no layout assumption and no pinned offset. Order-sensitive
//! reductions (`MATMUL`, `SUM`), formatted I/O (`WRITE`/`PRINT`) and whole
//! array updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ShellSolver statement for statement):
//! - `SAVE` semantics replicated: the driver SAVE arrays (`Indices`,
//!   `TotalSol`, `VisitsList`) live in wrapper-side SAVE storage
//!   (`acc_idx`, `shc_totalsol`, `shc_visits`) with identical lifetimes,
//!   while `LocalSol`/`LocalRHSForce` are Rust-owned buffers sized once
//!   at first-mesh sizing (`MaxElementDOFs`), exactly like the original
//!   one-shot `ALLOCATE`s. The internal procedures' own `SAVE` locals
//!   live in the transcribed `sh_*` procedures themselves.
//! - The `sh_*` procedures are exact transcriptions of the original
//!   internal procedures (same numerics, same call order); the two extra
//!   trailing kernel dummies only carry host variables (`TransientSim`,
//!   and the `BenchmarkProblem` pass-through) that were host-associated
//!   in the original.
//! - `ShellSolver_Init0` is a pure ValueList-setup body, so it is
//!   transcribed as a single wrapper (`shc_init0_run`, keeping the
//!   `CHARACTER` work in Fortran); the Rust entry only forwards handles.
//! - The principal-axes averaging block is transcribed as the same
//!   forced-`.FALSE.` guard the original executes, so it is skipped on
//!   every invocation exactly as in the reference.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn shc_begin(h_solver: *const c_void);
    fn shc_init0_run(transient: c_int);
    fn shc_default_start();
    fn shc_parallel() -> c_int;
    fn shc_cartesian(found: *mut c_int, val: *mut c_int);
    fn shc_element_has_p() -> c_int;
    fn shc_flag2(which: c_int, found: *mut c_int, val: *mut c_int);
    fn shc_shellmodelpar(lo: c_int, hi: c_int) -> c_int;
    fn shc_drilling_par(val: *mut c_double);
    fn shc_warn_drilling();
    fn shc_nonlin_pars(maxit: *mut c_int, tol: *mut c_double, relax: *mut c_double, has_relax: *mut c_int);
    fn shc_harmonic() -> c_int;
    fn shc_strain_par(which: c_int, found: *mut c_int, val: *mut c_int);
    fn shc_alpha(which: c_int, found: *mut c_int, val: *mut c_double);
    fn shc_max_elem_dofs() -> c_int;
    fn shc_principal_setup(n3: *mut c_int);
    fn shc_solid_coupling() -> c_int;
    fn shc_totalsol_setup(displace_active: c_int, modelpar: c_int);
    fn shc_totalsol_copy();
    fn shc_zero_variable();
    fn shc_relax_update(do_relax: c_int, relaxation: c_double);
    fn shc_bulkrhs_setup();
    fn shc_zero_buf(buf: *mut c_double, n: c_int);
    fn shc_nof_active() -> c_int;
    fn shc_elem_handle() -> *mut c_void;
    fn shc_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int, family: *mut c_int);
    fn shc_get_local_solution(buf: *mut c_double, ld1: c_int, ld2: c_int);
    fn shc_check_bubbles(nb: c_int);
    fn shc_check_remesh_p(remesh: c_int, pver: c_int);
    fn shc_check_drill_large(drill: c_int, large: c_int);
    fn shc_check_beam(n: c_int, nd: c_int);
    fn shc_default_initialize();
    fn shc_default_finish_bulk();
    fn shc_nof_bnd() -> c_int;
    fn shc_begin_boundary(t: c_int, n: *mut c_int, nd: *mut c_int, family: *mut c_int, is_active: *mut c_int);
    fn shc_boundary_parent(nd_parent: *mut c_int, has_parent: *mut c_int);
    fn shc_default_finish_bnd();
    fn shc_default_finish();
    fn shc_default_dirichlet();
    fn shc_default_solve() -> c_double;
    fn shc_constrained_fix();
    fn shc_update_force_swap(rbuf: *mut c_double, ln: c_int);
    fn shc_nonlin_first(nores0: *mut c_double, notractions: *mut c_int);
    fn shc_nonlin_check(tol: c_double, nonliniter: c_int, maxit: c_int, nores0: c_double, action: *mut c_int);
    fn shc_displace_fwd(modelpar: c_int);
    fn shc_benchmark_out(doarea: c_int, solvebench: c_int, shellarea: *mut c_double, blendarea: *mut c_double);
    fn shc_info_cartesian();
    fn shc_info_nl_banner(nl: c_int);
    fn shc_read_director();
    fn shc_check_orientation();
    fn shc_create_edges(fout: c_int);
    fn shc_average_director(h_elem: *const c_void, n: c_int, d_out: *mut c_double, plate_out: *mut c_int);
    fn shc_locframe_remesh(h_elem: *const c_void, lfn: *mut c_double, planar_out: *mut c_int, umbil_out: *mut c_int);
    fn shc_locframe_blend(h_elem: *const c_void, d_in: *const c_double, plate_in: c_int, lfn: *mut c_double, planar_out: *mut c_int, umbil_out: *mut c_int);
    fn shc_locpatch(h_elem: *const c_void, lfn: *const c_double, family: c_int, planar_in: c_int, umbil_in: c_int);
    fn shc_surface_area(h_elem: *const c_void, area: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn shc_shell_matrix(h_elem: *const c_void, n: c_int, nd: c_int, m: c_int, buf: *mut c_double, ld1: c_int, ld2: c_int, rbuf: *mut c_double, ln: c_int, large: c_int, nonlinbend: c_int, srm: c_int, msrm: c_int, shear_a: c_double, memb_a: c_double, stretch_a: c_double, bubbles: c_int, drill: c_int, drillpar: c_double, rot: c_int, mass: c_int, harm: c_int, area: *mut c_double, err: *mut c_double, bench: c_int, trans: c_int);
    #[allow(clippy::too_many_arguments)]
    fn shc_shell_matrix_cart(h_elem: *const c_void, n: c_int, nd: c_int, m: c_int, buf: *mut c_double, ld1: c_int, ld2: c_int, rbuf: *mut c_double, ln: c_int, large: c_int, nonlinbend: c_int, mass: c_int, harm: c_int, skip: c_int, remesh: c_int, bench: c_int, trans: c_int);
    #[allow(clippy::too_many_arguments)]
    fn shc_shell_boundary(h_elem: *const c_void, n: c_int, nd: c_int, m: c_int, buf: *mut c_double, ld1: c_int, ld2: c_int, rbuf: *mut c_double, ln: c_int, large: c_int, mass: c_int, harm: c_int, nd_parent: c_int, has_parent: c_int, cart: c_int, skip: c_int);
    #[allow(clippy::too_many_arguments)]
    fn shc_beam_matrix(h_elem: *const c_void, n: c_int, ndnb: c_int, nb: c_int, buf: *mut c_double, ld1: c_int, ld2: c_int, rbuf: *mut c_double, ln: c_int, trans: c_int, mass: c_int, harm: c_int, large: c_int, rot_not: c_int, drill: c_int);
    fn shc_solid_coupling_apply(h_model: *const c_void);
}

// shc_flag2 selector cases (must match accessors/shellsolver_acc.F90).
const F_SKIP: c_int = 1;
const F_REMESH: c_int = 2;
const F_DRILL: c_int = 3;
const F_BENCH: c_int = 4;
const F_LARGE: c_int = 5;
const F_NONLINBEND: c_int = 6;
const F_DISPLACE: c_int = 7;
const F_BUBBLES: c_int = 8;
const F_ROTATE: c_int = 9;
const F_PRAXES: c_int = 10;
const F_EDGECURVES: c_int = 11;
const F_SHELLAREA: c_int = 12;

const AUTOMATED: c_int = -1;
const NO_REDUCTION: c_int = 0;

/// First-mesh SAVE state: `LocalSol(ShellModelPar, MaxElementDOFs)` and
/// `LocalRHSForce((ShellModelPar+1)*MaxElementDOFs)`, sized once like the
/// original one-shot `ALLOCATE`s. Growth only triggers where the original
/// would overflow, and the leading sections used by the wrappers are
/// unaffected, so valid inputs never observe it.
struct ShellState {
    localsol: Vec<f64>,
    rhsforce: Vec<f64>,
    dims: (usize, usize),
}

impl ShellState {
    const fn new() -> Self {
        Self { localsol: Vec::new(), rhsforce: Vec::new(), dims: (0, 0) }
    }
    fn ensure(&mut self, modelpar: usize, maxdofs: usize) {
        if self.dims != (modelpar, maxdofs) {
            self.localsol = vec![0.0; modelpar * maxdofs];
            self.rhsforce = vec![0.0; (modelpar + 1) * maxdofs];
            self.dims = (modelpar, maxdofs);
        }
    }
}

static STATE: Mutex<ShellState> = Mutex::new(ShellState::new());

fn flag(which: c_int) -> (bool, bool) {
    let (mut found, mut val) = (0, 0);
    unsafe { shc_flag2(which, &mut found, &mut val) };
    (found != 0, val != 0)
}

fn b2i(b: bool) -> c_int {
    c_int::from(b)
}

/// `shellsolver_init0_` — solver initialization entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`; `Model` and `dt` are unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read if non-null.
#[export_name = "shellsolver_init0_"]
pub unsafe extern "C" fn shell_init0_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let trans = if transient.is_null() { 0 } else { unsafe { *transient } };
    unsafe {
        shc_begin(solver);
        shc_init0_run(trans);
    }
}

/// `shellsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list `(Model, Solver, dt, Transient)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes; `dt`
/// is unused by the body. `transient` selects mass assembly.
#[export_name = "shellsolver_"]
pub unsafe extern "C" fn shell_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let transient_sim = if transient.is_null() { false } else { unsafe { *transient != 0 } };
    unsafe { shell_main(model, solver, transient_sim) };
}

unsafe fn shell_main(model: *mut c_void, solver: *mut c_void, transient_sim: bool) {
    unsafe {
        shc_begin(solver);
        shc_default_start();
        let parallel = shc_parallel() != 0;

        // PART 0: key parameters.
        let mut cartesian = false;
        let mut pversion = false;
        {
            let (mut found, mut val) = (0, 0);
            shc_cartesian(&mut found, &mut val);
            if found == 0 {
                if shc_element_has_p() != 0 {
                    cartesian = true;
                    pversion = true;
                }
            } else {
                cartesian = val != 0;
            }
        }
        let mut skip_blending = flag(F_SKIP).1;
        let reparametrize = flag(F_REMESH).1;
        if reparametrize {
            skip_blending = true;
            cartesian = true;
            shc_check_remesh_p(b2i(reparametrize), b2i(pversion));
        }
        if cartesian {
            shc_info_cartesian();
        }
        let modelpar = if cartesian {
            shc_shellmodelpar(6, 9)
        } else {
            shc_shellmodelpar(6, 6)
        };
        let drilling = flag(F_DRILL).1;
        if drilling {
            shc_warn_drilling();
        }
        let mut drilling_par = 1.0;
        if drilling {
            shc_drilling_par(&mut drilling_par);
        }
        let mut solve_benchmark = flag(F_BENCH).1;
        let large = flag(F_LARGE).1;
        let (mut max_nonlin, mut nonlin_tol, mut relaxation) = (0, 0.0, 0.0);
        let mut relax = false;
        if large {
            let (mut has_relax, mut mi, mut tol, mut rx) = (0, 0, 0.0, 0.0);
            shc_nonlin_pars(&mut mi, &mut tol, &mut rx, &mut has_relax);
            max_nonlin = mi;
            nonlin_tol = tol;
            relaxation = rx;
            relax = has_relax != 0 && (relaxation - 1.0).abs() > f64::EPSILON;
            shc_check_drill_large(b2i(drilling), b2i(large));
            solve_benchmark = false;
            shc_bulkrhs_setup();
        }
        let nonlinear_bending = {
            let (found, v) = flag(F_NONLINBEND);
            if found { v } else { true }
        };
        let mesh_displace_active = flag(F_DISPLACE).1;
        let harmonic = shc_harmonic() != 0;
        let mass_assembly = transient_sim || harmonic;
        let mut strain_method = {
            let (found, v) = (0, 0);
            let mut f = found;
            let mut val = v;
            shc_strain_par(1, &mut f, &mut val);
            if f == 0 { AUTOMATED } else { val }
        };
        if strain_method != NO_REDUCTION {
            strain_method = AUTOMATED;
        }
        let mut membrane_method = {
            let mut f = 0;
            let mut val = 0;
            shc_strain_par(2, &mut f, &mut val);
            if f == 0 { AUTOMATED } else { val }
        };
        if membrane_method != NO_REDUCTION {
            membrane_method = AUTOMATED;
        }
        let mut shear_alpha = 0.0;
        let mut membrane_alpha = 0.0;
        let mut stretch_alpha = 0.0;
        {
            let mut f = 0;
            shc_alpha(1, &mut f, &mut shear_alpha);
            shc_alpha(2, &mut f, &mut membrane_alpha);
            shc_alpha(3, &mut f, &mut stretch_alpha);
        }
        let bubbles = flag(F_BUBBLES).1;
        let rotate_dofs = flag(F_ROTATE).1;
        let mut save_principal = flag(F_PRAXES).1;
        if save_principal {
            let mut n3 = 0;
            shc_principal_setup(&mut n3);
        }
        let solid_coupling = shc_solid_coupling() != 0;

        let maxdofs = shc_max_elem_dofs();
        let (ld1, ld2) = (modelpar as usize, maxdofs as usize);
        {
            let mut st = STATE.lock().unwrap();
            st.ensure(ld1, ld2);
        }
        shc_totalsol_setup(b2i(mesh_displace_active), modelpar);

        let mut compute_area = false;
        let mut blending_area = 0.0;
        if !skip_blending {
            // PART I: director data + surface orientation.
            shc_read_director();
            shc_check_orientation();
            // PART II: curved-edge geometry descriptions.
            let curve_out = flag(F_EDGECURVES).1;
            shc_create_edges(b2i(curve_out));
            compute_area = flag(F_SHELLAREA).1;
        }

        // PART III: reparametrization + nonlinear assembly loop.
        // (`total_err` is reset per iteration like the original `TotalErr`;
        // the kernel `Error` dummy never reads it, so the value is local.)
        let mut shell_area = 0.0;
        let mut nonlin_res0 = 0.0;
        let mut _notractions = false;
        let max_iter = max_nonlin + 1;
        for nonlin_iter in 1..=max_iter {
            if nonlin_iter < max_iter {
                shc_info_nl_banner(nonlin_iter);
            }
            shc_totalsol_copy();
            shc_default_initialize();
            shell_area = 0.0;
            let mut total_err = 0.0;
            let active = shc_nof_active();
            for k in 1..=active {
                let (mut n, mut nd, mut nb, mut family) = (0, 0, 0, 0);
                shc_begin_active(k, &mut n, &mut nd, &mut nb, &mut family);
                if family != 3 && family != 4 {
                    continue;
                }
                let h_elem = shc_elem_handle();
                {
                    let mut st = STATE.lock().unwrap();
                    if large {
                        shc_get_local_solution(st.localsol.as_mut_ptr(), modelpar, maxdofs);
                    } else {
                        let ln = (ld1 * ld2) as c_int;
                        shc_zero_buf(st.localsol.as_mut_ptr(), ln);
                    }
                }
                let apply_bubbles = bubbles && nd == family;
                shc_check_bubbles(nb);
                if nonlin_iter == 1 && (!skip_blending || reparametrize) {
                    let mut lfn = [0.0f64; 48];
                    let (mut planar, mut umbil) = (0, 0);
                    if reparametrize {
                        shc_locframe_remesh(h_elem, lfn.as_mut_ptr(), &mut planar, &mut umbil);
                    } else {
                        let (mut d, mut plate) = ([0.0; 3], 0);
                        shc_average_director(h_elem, n, d.as_mut_ptr(), &mut plate);
                        shc_locframe_blend(
                            h_elem,
                            d.as_ptr(),
                            plate,
                            lfn.as_mut_ptr(),
                            &mut planar,
                            &mut umbil,
                        );
                    }
                    shc_locpatch(h_elem, lfn.as_ptr(), family, planar, umbil);
                    if compute_area {
                        shc_surface_area(h_elem, &mut blending_area);
                    }
                }
                {
                    let mut st = STATE.lock().unwrap();
                    let ln = ((modelpar + 1) * maxdofs) as c_int;
                    if cartesian {
                        shc_shell_matrix_cart(
                            h_elem, n, nd, modelpar,
                            st.localsol.as_mut_ptr(), modelpar, maxdofs,
                            st.rhsforce.as_mut_ptr(), ln,
                            b2i(large), b2i(nonlinear_bending),
                            b2i(mass_assembly), b2i(harmonic),
                            b2i(skip_blending), b2i(reparametrize),
                            b2i(solve_benchmark), b2i(transient_sim),
                        );
                    } else {
                        shc_shell_matrix(
                            h_elem, n, nd, modelpar,
                            st.localsol.as_mut_ptr(), modelpar, maxdofs,
                            st.rhsforce.as_mut_ptr(), ln,
                            b2i(large), b2i(nonlinear_bending),
                            strain_method, membrane_method,
                            shear_alpha, membrane_alpha, stretch_alpha,
                            b2i(apply_bubbles), b2i(drilling), drilling_par,
                            b2i(rotate_dofs), b2i(mass_assembly), b2i(harmonic),
                            &mut shell_area, &mut total_err,
                            b2i(solve_benchmark), b2i(transient_sim),
                        );
                    }
                    if large && nonlin_iter == 1 {
                        shc_update_force_swap(st.rhsforce.as_mut_ptr(), ln);
                    }
                }
            }
            // Principal-direction averaging: forced off above, as in source.
            save_principal = false;
            if save_principal {
                // Dead in practice (mirrors the original guard).
            }

            // Beams.
            for k in 1..=active {
                let (mut n, mut nd, mut nb, mut family) = (0, 0, 0, 0);
                shc_begin_active(k, &mut n, &mut nd, &mut nb, &mut family);
                if family != 2 {
                    continue;
                }
                shc_check_beam(n, nd);
                let h_elem = shc_elem_handle();
                {
                    let mut st = STATE.lock().unwrap();
                    if large {
                        shc_get_local_solution(st.localsol.as_mut_ptr(), modelpar, maxdofs);
                    } else {
                        let ln = (ld1 * ld2) as c_int;
                        shc_zero_buf(st.localsol.as_mut_ptr(), ln);
                    }
                    let ln = ((modelpar + 1) * maxdofs) as c_int;
                    shc_beam_matrix(
                        h_elem, n, nd + 1, 1,
                        st.localsol.as_mut_ptr(), modelpar, maxdofs,
                        st.rhsforce.as_mut_ptr(), ln,
                        b2i(transient_sim), b2i(mass_assembly), b2i(harmonic),
                        b2i(large), b2i(!rotate_dofs), b2i(drilling),
                    );
                    if large && nonlin_iter == 1 {
                        shc_update_force_swap(st.rhsforce.as_mut_ptr(), ln);
                    }
                }
            }
            shc_default_finish_bulk();

            let bnd_active = shc_nof_bnd();
            for k in 1..=bnd_active {
                let (mut n, mut nd, mut family, mut is_active) = (0, 0, 0, 0);
                shc_begin_boundary(k, &mut n, &mut nd, &mut family, &mut is_active);
                if is_active == 0 {
                    continue;
                }
                let h_elem = shc_elem_handle();
                let (mut nd_parent, mut has_parent) = (0, 0);
                {
                    let mut st = STATE.lock().unwrap();
                    if large {
                        shc_get_local_solution(st.localsol.as_mut_ptr(), modelpar, maxdofs);
                        shc_boundary_parent(&mut nd_parent, &mut has_parent);
                    } else {
                        let ln = (ld1 * ld2) as c_int;
                        shc_zero_buf(st.localsol.as_mut_ptr(), ln);
                    }
                    let ln = ((modelpar + 1) * maxdofs) as c_int;
                    shc_shell_boundary(
                        h_elem, n, nd, modelpar,
                        st.localsol.as_mut_ptr(), modelpar, maxdofs,
                        st.rhsforce.as_mut_ptr(), ln,
                        b2i(large), b2i(mass_assembly), b2i(harmonic),
                        nd_parent, has_parent, b2i(cartesian), b2i(skip_blending),
                    );
                    if large && nonlin_iter == 1 {
                        shc_update_force_swap(st.rhsforce.as_mut_ptr(), ln);
                    }
                }
            }
            shc_default_finish_bnd();
            shc_default_finish();
            shc_default_dirichlet();
            if solid_coupling {
                shc_solid_coupling_apply(model);
            }
            shc_constrained_fix();

            if large {
                if nonlin_iter == 1 {
                    let mut nt = 0;
                    shc_nonlin_first(&mut nonlin_res0, &mut nt);
                    _notractions = nt != 0;
                }
                let mut action = 0;
                shc_nonlin_check(nonlin_tol, nonlin_iter, max_nonlin, nonlin_res0, &mut action);
                if action != 0 {
                    break;
                }
            }
            if large {
                shc_zero_variable();
            }
            let _ = shc_default_solve();
            if large {
                shc_relax_update(b2i(relax), relaxation);
            }
        }

        // PART IV: postprocess.
        if mesh_displace_active {
            shc_displace_fwd(modelpar);
        }
        shc_benchmark_out(
            b2i(compute_area),
            b2i(solve_benchmark && !parallel),
            &mut shell_area,
            &mut blending_area,
        );
    }
}
