//! Port of `fem/src/modules/StructuredFlowLine.F90` (Elmer @ the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/StructuredFlowLine.F90` — two F77-global
//! entries `structuredflowline_` and `structuredflowline_init_` (proven by
//! `nm -D` on the reference `.so`); one internal procedure
//! (`ComputeVeloAtPoint`) transcribed verbatim into the accessor module.
//! The emitted dynamic symbols are exactly `structuredflowline_` and
//! `structuredflowline_init_`.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/structuredflowline_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), whole-array updates and every formatted `WRITE`/`Info` stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors StructuredFlowLine.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: all outer-procedure SAVE
//!   (`VisitedTimes`, `Initialized`, `UnitVector`, `Coord`, the pointer
//!   tables, `Dim`, `DispPerm`, `MaskPerm`, `MaskExist`,
//!   `ActiveDirection`) plus the internal procedure's SAVE (`Visited`,
//!   `Elements`, `RootQuadrant`, `CurrentElement`, `ElementNodes`,
//!   `Basis`, `dBasisdx`) and the `at0/at1` timing locals live as
//!   wrapper-side SAVE with identical lifetimes. The Rust side holds no
//!   state across calls (no `static`), so no sharing guard is needed;
//!   the solver entries run serially, as in the original.
//! - `TopMode`, `BotMode`, `TopField`, `BotField`, `TopPerm`, `BotPerm`,
//!   `Field`, `Surface`, `DOFs`, `ibot` are declared+SAVE'd in the
//!   original but never used (`NULLIFY(Field)` aside), so they need no
//!   mirror. Dead `IF(.FALSE.) PRINT` markers are elided.
//! - Vector temporaries (`v1`, `v2`, `r1`, `r2`, `dr`, `vave`, `drproj`,
//!   `sumdr`, `VelCor`, `GlobalCoords`) are Rust `[f64; 3]` arrays with
//!   per-component operations in source order; every component result is
//!   uniquely determined, so evaluation order cannot diverge. All `SUM`
//!   reductions stay in wrappers.
//! - `TINY(CoordMin)` is `f64::MIN_POSITIVE`; unset `Min/Max Coordinate`
//!   default to `-f64::MAX`/`f64::MAX` (`-HUGE`/`HUGE`); the
//!   `1.0e-4` dot-product default keeps its single-precision source
//!   spelling inside the wrapper. IEEE division-by-zero behavior
//!   (`0/0`, `x/0`) is replicated by performing the same divisions in
//!   the same order in `f64`.
//! - The displacement block keeps the reference's leftover-`i` quirk
//!   verbatim (the mask guard reads the outer `i`-loop's post-loop
//!   value, passed through as `nsize+1`): the proof exercises defined
//!   paths (Displacement Mode off, or no mask) for the differentials.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored); `dt`/`transient` are
//!   received but never read.
//!
//! # Safety
//!
//! The entry arguments must be the Elmer addresses the loader passes
//! (`Solver` the `Solver_t` address, `Model` the `Model_t` address).
//! All wrapper calls execute on the Elmer solver thread only.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sfl_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sfl_info_header();
    fn sfl_free_quadrant();
    fn sfl_get_true_iters() -> c_int;
    fn sfl_visited() -> c_int;
    fn sfl_mark_visited();
    fn sfl_get_ave_order(truelines: c_int) -> c_int;
    fn sfl_get_ave_mode() -> c_int;
    fn sfl_initialized() -> c_int;
    fn sfl_mark_initialized();
    #[allow(clippy::too_many_arguments)]
    fn sfl_velo_at_point(
        g: *mut c_double,
        node: c_int,
        with_grad: c_int,
        velo: *mut c_double,
        grad: *mut c_double,
        outmesh: *mut c_int,
    );
    fn sfl_enter_init() -> c_int;
    fn sfl_get_dim() -> c_int;
    fn sfl_get_actdir() -> c_int;
    fn sfl_is_mask() -> c_int;
    fn sfl_init_coord() -> c_int;
    fn sfl_init_unit();
    fn sfl_init_mask() -> c_int;
    fn sfl_alloc_pointers();
    fn sfl_get_eps() -> c_double;
    fn sfl_num_bulk() -> c_int;
    fn sfl_max_elem_nodes() -> c_int;
    fn sfl_elem_begin(elem: c_int) -> c_int;
    fn sfl_elem_nodeidx(n: c_int, buf: *mut c_int);
    fn sfl_elem_coords(n: c_int, xs: *mut c_double, ys: *mut c_double, zs: *mut c_double);
    fn sfl_mask_at(i: c_int) -> c_int;
    fn sfl_vec_len(ev: *const c_double) -> c_double;
    fn sfl_vec_dot(ev: *const c_double) -> c_double;
    fn sfl_set_updown(ii: c_int, jj: c_int, dir: c_int);
    fn sfl_propagate() -> c_int;
    fn sfl_tic();
    fn sfl_get_flow() -> c_int;
    fn sfl_get_disp() -> c_int;
    fn sfl_get_hitstyle() -> c_int;
    fn sfl_get_dispmode() -> c_int;
    fn sfl_coord_limits(
        cmin: *mut c_double,
        gotmin: *mut c_int,
        cmax: *mut c_double,
        gotmax: *mut c_int,
    );
    fn sfl_nsize() -> c_int;
    fn sfl_disp_perm(i: c_int) -> c_int;
    fn sfl_top_at(i: c_int) -> c_int;
    fn sfl_down_at(k: c_int) -> c_int;
    fn sfl_bot_at(k: c_int) -> c_int;
    fn sfl_flow_at(k: c_int, dim: c_int, buf: *mut c_double);
    fn sfl_node_xyz(k: c_int, buf: *mut c_double);
    fn sfl_mark_hit(k: c_int, val: c_double);
    fn sfl_norm3(a: *const c_double) -> c_double;
    fn sfl_norm3_sum(a: *const c_double, b: *const c_double) -> c_double;
    fn sfl_harddisp_set(k: c_int, dim: c_int, s: *const c_double);
    fn sfl_report_norm();
    fn sfl_relax_apply();
    fn sfl_move_nodes(leftover_i: c_int);
    fn sfl_report_shots(shots: c_int, hits: c_int);
    fn sfl_report_maptime();
    fn sfl_init_begin(h_model: *const c_void, h_solver: *const c_void, dim: *mut c_int);
    fn sfl_init_nomatrix();
    fn sfl_init_var(dim: c_int);
}

/// `structuredflowline_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "structuredflowline_"]
pub unsafe extern "C" fn structuredflowline_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sfl_begin(model, solver) };

    unsafe { sfl_info_header() };

    // The mesh is changing, thus always rebuild the quadrant tree.
    unsafe { sfl_free_quadrant() };

    let true_iters = unsafe { sfl_get_true_iters() };
    let visited = unsafe { sfl_visited() };
    let truelines = true_iters > visited;

    // One is a good default for the averaging order; zero is a good
    // default for the averaging method (defaults inside the wrapper).
    let aveorder = unsafe { sfl_get_ave_order(i32::from(truelines)) };
    let avemode = unsafe { sfl_get_ave_mode() };

    if unsafe { sfl_initialized() } == 0 {
        unsafe { sfl_tic() };

        // Choose active direction coordinate and set corresponding unit
        // vector (associates Coord, aborts via Fatal on invalid input).
        let _dim = unsafe { sfl_enter_init() };
        let _actdir = unsafe { sfl_init_coord() };
        unsafe { sfl_init_unit() };

        // Make the mask used in matrix structures and permutation.
        let maskexist = unsafe { sfl_init_mask() } != 0;
        let _ = maskexist;

        // Allocate pointers to top and bottom, and temporary pointers
        // up and down.
        unsafe { sfl_alloc_pointers() };

        // Determine the pointer up and down using dot product as
        // criterion.
        let eps = unsafe { sfl_get_eps() };
        let nbulk = unsafe { sfl_num_bulk() };
        let maxnodes = unsafe { sfl_max_elem_nodes() } as usize;
        let mut idx = vec![0 as c_int; maxnodes];
        let mut xs = vec![0.0 as c_double; maxnodes];
        let mut ys = vec![0.0 as c_double; maxnodes];
        let mut zs = vec![0.0 as c_double; maxnodes];
        for elem in 1..=nbulk {
            let n = unsafe { sfl_elem_begin(elem) } as usize;
            unsafe {
                sfl_elem_nodeidx(n as c_int, idx.as_mut_ptr());
                sfl_elem_coords(
                    n as c_int,
                    xs.as_mut_ptr(),
                    ys.as_mut_ptr(),
                    zs.as_mut_ptr(),
                );
            }
            let maskexist = unsafe { sfl_is_mask() } != 0;
            if maskexist {
                let mut any = false;
                for ii in idx.iter().take(n) {
                    if unsafe { sfl_mask_at(*ii) } != 0 {
                        any = true;
                        break;
                    }
                }
                if !any {
                    continue;
                }
            }
            for a in 0..n {
                let ii = idx[a];
                if maskexist && unsafe { sfl_mask_at(ii) } == 0 {
                    continue;
                }
                for b in (a + 1)..n {
                    let jj = idx[b];
                    if maskexist && unsafe { sfl_mask_at(jj) } == 0 {
                        continue;
                    }
                    let ev = [xs[b] - xs[a], ys[b] - ys[a], zs[b] - zs[a]];
                    let length = unsafe { sfl_vec_len(ev.as_ptr()) };
                    let dotpro = unsafe { sfl_vec_dot(ev.as_ptr()) } / length;
                    if dotpro > 1.0 - eps {
                        unsafe { sfl_set_updown(ii, jj, 1) };
                    } else if dotpro < eps - 1.0 {
                        unsafe { sfl_set_updown(ii, jj, -1) };
                    }
                }
            }
        }

        // Pointer to top and bottom are found recursively using up and
        // down (rounds/timing reports inside the wrapper).
        let _rounds = unsafe { sfl_propagate() };

        unsafe { sfl_mark_initialized() };
    }
    // End of initialization.
    unsafe { sfl_tic() };

    // Get the velocity variables and the result field for suggested
    // displacement (Fatal on missing, inside the wrappers).
    let _flowdofs = unsafe { sfl_get_flow() };
    let _dispdofs = unsafe { sfl_get_disp() };
    let gothit = unsafe { sfl_get_hitstyle() } != 0;

    let dispmode = unsafe { sfl_get_dispmode() } != 0;

    // Map to the streamline.
    let alpha = 0.5 as c_double;
    let (mut cmin, mut gotmin, mut cmax, mut gotmax) = (0.0, 0, 0.0, 0);
    unsafe { sfl_coord_limits(&mut cmin, &mut gotmin, &mut cmax, &mut gotmax) };
    let _ = (gotmin, gotmax);
    if gotmin == 0 {
        cmin = -f64::MAX;
    }
    if gotmax == 0 {
        cmax = f64::MAX;
    }

    let nsize = unsafe { sfl_nsize() };
    let dim = unsafe { sfl_get_dim() };
    let act0 = (unsafe { sfl_get_actdir() } as usize).saturating_sub(1);

    let mut shots: c_int = 0;
    let mut hits: c_int = 0;

    let mut v1 = [0.0 as c_double; 3];
    let mut v2 = [0.0 as c_double; 3];
    let mut r1 = [0.0 as c_double; 3];
    let mut r2 = [0.0 as c_double; 3];
    let mut dr = [0.0 as c_double; 3];
    let mut vave = [0.0 as c_double; 3];
    let mut drproj = [0.0 as c_double; 3];
    let mut sumdr: [c_double; 3];
    let mut velcor: [c_double; 3];
    let mut gcoords = [0.0 as c_double; 3];
    let mut velo = [0.0 as c_double; 3];
    let mut grad = [0.0 as c_double; 9];

    let maskexist = unsafe { sfl_is_mask() } != 0;

    for i in 1..=nsize {
        if unsafe { sfl_disp_perm(i) } == 0 {
            continue;
        }
        if maskexist && unsafe { sfl_mask_at(i) } == 0 {
            continue;
        }
        let itop = unsafe { sfl_top_at(i) };
        if i != itop {
            continue;
        }

        // First node of streamline is initialized.
        let mut k = i;
        unsafe {
            sfl_flow_at(k, dim, v1.as_mut_ptr());
            sfl_node_xyz(k, r1.as_mut_ptr());
        }
        sumdr = [0.0; 3];

        for _j in 0..nsize {
            k = unsafe { sfl_down_at(k) };
            if unsafe { sfl_disp_perm(k) } == 0 {
                continue;
            }

            unsafe { sfl_node_xyz(k, r2.as_mut_ptr()) };
            for m in 0..3 {
                dr[m] = r2[m] - r1[m];
                r1[m] = r2[m];
            }
            unsafe { sfl_flow_at(k, dim, v2.as_mut_ptr()) };

            // Check the passive regions; for no-slip condition the
            // surface is always fixed.
            if r2[act0] > cmax || r2[act0] < cmin || v2[act0].abs() < f64::MIN_POSITIVE {
                v1 = v2;
                continue;
            }

            // Simple average velocity from the mean of surface
            // velocities.
            for m in 0..3 {
                vave[m] = alpha * v2[m] + (1.0 - alpha) * v1[m];
            }

            // A more complex velocity estimation is preferable in the
            // start.
            if truelines {
                drproj = [0.0; 3];
                let mut sumds = unsafe { sfl_norm3(sumdr.as_ptr()) };
                for _l in 0..aveorder {
                    if sumds > 1.0e-20 {
                        // Averaging of velocity may be done in two ways:
                        // 1) taking velocity at average of r1 and r2
                        // 2) taking average of velocities at r1 and r2
                        if avemode == 1 {
                            for m in 0..3 {
                                gcoords[m] =
                                    0.5 * (r1[m] + r2[m]) + sumdr[m] + 0.5 * drproj[m];
                            }
                        } else {
                            for m in 0..3 {
                                gcoords[m] = r2[m] + sumdr[m] + drproj[m];
                            }
                        }
                        let mut outmesh = 0;
                        unsafe {
                            sfl_velo_at_point(
                                gcoords.as_mut_ptr(),
                                k,
                                0,
                                velo.as_mut_ptr(),
                                grad.as_mut_ptr(),
                                &mut outmesh,
                            );
                        }
                        shots += 1;

                        if outmesh != 0 {
                            if avemode == 1 {
                                for m in 0..3 {
                                    gcoords[m] = 0.5 * (r1[m] + r2[m]);
                                }
                            } else {
                                for m in 0..3 {
                                    gcoords[m] = r2[m] + sumdr[m] + drproj[m];
                                }
                            }
                            let mut outmesh2 = 0;
                            unsafe {
                                sfl_velo_at_point(
                                    gcoords.as_mut_ptr(),
                                    k,
                                    1,
                                    velo.as_mut_ptr(),
                                    grad.as_mut_ptr(),
                                    &mut outmesh2,
                                );
                            }
                            if outmesh2 == 0 {
                                if gothit {
                                    unsafe { sfl_mark_hit(k, -(_l as c_double + 1.0)) };
                                }
                                velcor = [0.0; 3];
                                for n in 0..3 {
                                    if n == act0 {
                                        continue;
                                    }
                                    for m in 0..3 {
                                        velcor[m] += grad[m + 3 * n] * dr[n];
                                    }
                                }
                                if avemode == 1 {
                                    for m in 0..3 {
                                        vave[m] = velo[m] + velcor[m];
                                    }
                                } else {
                                    for m in 0..3 {
                                        v2[m] = velo[m] + velcor[m];
                                        vave[m] = alpha * v2[m] + (1.0 - alpha) * v1[m];
                                    }
                                }
                            }
                        } else {
                            if gothit {
                                unsafe { sfl_mark_hit(k, _l as c_double + 1.0) };
                            }
                            hits += 1;
                            if avemode == 1 {
                                vave = velo;
                            } else {
                                v2 = velo;
                                for m in 0..3 {
                                    vave[m] = alpha * v2[m] + (1.0 - alpha) * v1[m];
                                }
                            }
                        }
                    }
                    let denom = vave[act0].abs();
                    for m in 0..3 {
                        drproj[m] = vave[m] * dr[act0].abs() / denom - dr[m];
                    }
                    sumds = unsafe { sfl_norm3_sum(sumdr.as_ptr(), drproj.as_ptr()) };
                }
            } else {
                let denom = vave[act0].abs();
                for m in 0..3 {
                    drproj[m] = vave[m] * dr[act0].abs() / denom - dr[m];
                }
            }

            for m in 0..3 {
                sumdr[m] += drproj[m];
            }
            sumdr[act0] = 0.0;

            unsafe { sfl_harddisp_set(k, dim, sumdr.as_ptr()) };

            v1 = v2;

            // This is the bottom and we're done?
            if k == unsafe { sfl_bot_at(k) } {
                break;
            }
        }
    }

    unsafe { sfl_report_norm() };

    unsafe { sfl_relax_apply() };

    if dispmode {
        // Leftover-`i` quirk, verbatim: the original reads MaskPerm(i)
        // with the outer i-loop's post-loop value (nsize+1).
        unsafe { sfl_move_nodes(nsize + 1) };
    }

    if truelines {
        unsafe { sfl_report_shots(shots, hits) };
    }

    unsafe { sfl_report_maptime() };

    unsafe { sfl_mark_visited() };
}

/// `structuredflowline_init_` — initialization of the primary solver.
/// Mirrors `StructuredFlowLine_init` statement for statement.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "structuredflowline_init_"]
pub unsafe extern "C" fn structuredflowline_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut dim: c_int = 0;
    unsafe { sfl_init_begin(model, solver, &mut dim) };

    unsafe { sfl_init_nomatrix() };

    unsafe { sfl_init_var(dim) };
}
