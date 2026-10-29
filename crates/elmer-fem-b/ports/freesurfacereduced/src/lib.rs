//! Port of `fem/src/modules/FreeSurfaceReduced.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T022-A, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/FreeSurfaceReduced.F90` — one F77-global
//! entry `freesurfacereduced_` with no `CONTAINS` internals. The emitted
//! dynamic symbol is exactly `freesurfacereduced_` (proven by `nm -D` on
//! the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/freesurfacereduced_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MAX`/`MIN` extent scans), formatted `WRITE` + `Info` reports, the
//! neighbour search, and the whole-array/mapping updates stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirrors FreeSurfaceReduced.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: first-visit sizing of
//!   `FreeSurfacePoints`/`LeftNeighbours`/`RightNeighbours`/`OrigCoords`
//!   at `Model % NumberOfNodes`, `UpwindPoints` at `NoFreeSurfaces`, the
//!   `coords`/`coords2` pointer aliases into `Mesh % Nodes`, plus the
//!   `SubroutineVisited`/`AllocationsDone`/`DrawDirection` flags — all as
//!   wrapper-side SAVE state with identical lifetimes.
//! - `Model`, `dt`, and the transient flag are unused by the body, exactly
//!   as in the original (received and ignored); `N = MaxElementNodes` is
//!   kept as a touch read; the double `Solver % Variable % Norm` store
//!   collapses to the final value, exactly as in the original.
//! - Allocation failure calls the same `Warn('FreeSurfaceReduced', ...)`
//!   (the original warns, it does not `Fatal`), instead of aborting.
//! - `MAX`/`MIN` scans use the same strict-`>`/`<` update order as the
//!   intrinsics on valid (non-NaN) inputs; `HUGE` bounds are `f64::MAX`.
//! - The flux integration is transcribed term for term in `f64`
//!   (`x**2`/`x**3` as explicit products, as the compiler emits them).
//! - No Rust threads: the solver entries run serially, so no sharing
//!   guard is needed (all Rust state is entry-local; shared state lives
//!   in the Fortran wrappers, as in the original).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn fsa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn fsa_find_flow();
    fn fsa_info_header();
    fn fsa_max_nodes() -> c_int;
    fn fsa_num_nodes() -> c_int;
    fn fsa_num_bcs() -> c_int;
    fn fsa_alloc_done() -> c_int;
    fn fsa_alloc_first(n: c_int);
    fn fsa_zero_marks();
    fn fsa_visited() -> c_int;
    fn fsa_inc_visited();
    fn fsa_num_bulk() -> c_int;
    fn fsa_num_bnd() -> c_int;
    fn fsa_begin_elem(t: c_int, elemcode: *mut c_int, n: *mut c_int, constraint: *mut c_int);
    fn fsa_bc_tag(k: c_int) -> c_int;
    fn fsa_bc_reduced(k: c_int, found: *mut c_int) -> c_int;
    fn fsa_bc_surface(k: c_int, surf: *mut c_int, found: *mut c_int);
    fn fsa_bc_bottom(k: c_int, found: *mut c_int) -> c_int;
    fn fsa_apply_surface(n: c_int, surf: c_int);
    fn fsa_apply_bottom(n: c_int);
    fn fsa_note_surface(surf: c_int);
    fn fsa_nofree() -> c_int;
    fn fsa_info_nsurf();
    fn fsa_alloc_upwind();
    fn fsa_fsp(t: c_int) -> c_int;
    fn fsa_node_x(t: c_int) -> c_double;
    fn fsa_node_y(t: c_int) -> c_double;
    fn fsa_set_drawdir(d: c_int);
    fn fsa_drawdir() -> c_int;
    fn fsa_flowt(t: c_int, comp: c_int) -> c_double;
    fn fsa_flowk(k: c_int, comp: c_int) -> c_double;
    fn fsa_var_perm(t: c_int) -> c_int;
    fn fsa_move_perm(t: c_int) -> c_int;
    fn fsa_info_drawdir();
    fn fsa_neighbours();
    fn fsa_bind_coords();
    fn fsa_get_upwind(s: c_int) -> c_int;
    fn fsa_set_upwind(s: c_int, t: c_int);
    fn fsa_info_upwind(s: c_int);
    fn fsa_zero_move();
    fn fsa_upwind_coord(t: c_int) -> c_double;
    fn fsa_coord(t: c_int) -> c_double;
    fn fsa_coord2(t: c_int) -> c_double;
    fn fsa_left(t: c_int) -> c_int;
    fn fsa_right(t: c_int) -> c_int;
    fn fsa_warn_perm();
    fn fsa_warn_limits();
    fn fsa_warn_noflow();
    fn fsa_warn_dx0();
    fn fsa_axisym() -> c_int;
    fn fsa_relax(r: *mut c_double, found: *mut c_int);
    fn fsa_mapping(m: *mut c_int, found: *mut c_int);
    fn fsa_map_chain(istart: c_int, iend: c_int, x0: c_double, xnew: c_double, xold: c_double, relax: c_double);
    fn fsa_move_point(iend: c_int, xnew: c_double, relax: c_double);
    fn fsa_info_surf(
        s: c_int,
        multi: c_int,
        fsum0: c_double,
        maxerror: c_double,
        thickmin: c_double,
        thickmax: c_double,
    );
    fn fsa_set_norm(n: c_double);
    fn fsa_add_results(fsum0: c_double, maxerror: c_double, thickmin: c_double, thickmax: c_double);
}

/// `freesurfacereduced_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "freesurfacereduced_"]
pub unsafe extern "C" fn freesurface_reduced_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { fsa_begin(model, solver) };
    unsafe { fsa_info_header() };
    unsafe { fsa_find_flow() };

    let mut mapping = 0;
    let mut gotit = 0;
    unsafe { fsa_mapping(&mut mapping, &mut gotit) };
    let mut perform_mapping = mapping != 0;
    if gotit == 0 {
        perform_mapping = true;
    }

    // Permanent storage, allocated on first visit only.
    if unsafe { fsa_alloc_done() } == 0 {
        // N = Model % MaxElementNodes (touch read; the result is unused,
        // exactly as in the original).
        let _maxnodes = unsafe { fsa_max_nodes() };
        let nn = unsafe { fsa_num_nodes() };
        unsafe { fsa_alloc_first(nn) };
    }

    let axisymmetric = unsafe { fsa_axisym() } != 0;

    let mut relax = 0.0;
    let mut found = 0;
    unsafe { fsa_relax(&mut relax, &mut found) };
    if found == 0 {
        relax = 1.0;
    }

    // First-visit initialization: mark free-surface nodes, decide the
    // draw direction, build neighbour chains, bind the coordinate
    // pointers, find upwind nodes, zero the displacement.
    if unsafe { fsa_visited() } == 0 {
        let nn = unsafe { fsa_num_nodes() };
        unsafe { fsa_zero_marks() };

        let bulk = unsafe { fsa_num_bulk() };
        let nbnd = unsafe { fsa_num_bnd() };
        for t in (bulk + 1)..=(bulk + nbnd) {
            let (mut elemcode, mut n, mut constraint) = (0, 0, 0);
            unsafe { fsa_begin_elem(t, &mut elemcode, &mut n, &mut constraint) };
            if elemcode == 101 {
                continue;
            }
            let nbcs = unsafe { fsa_num_bcs() };
            for k in 1..=nbcs {
                if unsafe { fsa_bc_tag(k) } != constraint {
                    continue;
                }
                let mut fnd = 0;
                let is_fs = unsafe { fsa_bc_reduced(k, &mut fnd) };
                if fnd != 0 && is_fs != 0 {
                    let (mut surf, mut got) = (0, 0);
                    unsafe { fsa_bc_surface(k, &mut surf, &mut got) };
                    if got == 0 {
                        surf = 1;
                    }
                    unsafe {
                        fsa_apply_surface(n, surf);
                        fsa_note_surface(surf);
                    }
                }
                let mut fnd2 = 0;
                let is_bottom = unsafe { fsa_bc_bottom(k, &mut fnd2) };
                if fnd2 != 0 && is_bottom != 0 {
                    unsafe { fsa_apply_bottom(n) };
                }
            }
        }

        unsafe { fsa_info_nsurf() };
        unsafe { fsa_alloc_upwind() };

        // Extent of the free surface in x and y.
        let mut got = false;
        let (mut x1, mut x2, mut y1, mut y2) = (0.0, 0.0, 0.0, 0.0);
        for t in 1..=nn {
            if unsafe { fsa_fsp(t) } <= 0 {
                continue;
            }
            let nx = unsafe { fsa_node_x(t) };
            let ny = unsafe { fsa_node_y(t) };
            if !got {
                got = true;
                x1 = nx;
                x2 = nx;
                y1 = ny;
                y2 = ny;
            } else {
                if nx > x1 {
                    x1 = nx;
                }
                if nx < x2 {
                    x2 = nx;
                }
                if ny > y1 {
                    y1 = ny;
                }
                if ny < y2 {
                    y2 = ny;
                }
            }
        }

        // In x or y direction.
        let mut drawdir = if (x1 - x2) > (y1 - y2) { 1 } else { 2 };

        // Velocity range along the draw direction.
        got = false;
        let (mut v1, mut v2) = (0.0, 0.0);
        for t in 1..=nn {
            if unsafe { fsa_fsp(t) } <= 0 {
                continue;
            }
            let vv = unsafe { fsa_flowt(t, drawdir) };
            if !got {
                got = true;
                v1 = vv;
                v2 = vv;
            } else {
                if vv > v1 {
                    v1 = vv;
                }
                if vv < v2 {
                    v2 = vv;
                }
            }
        }

        // In positive or negative direction.
        if -v2 > v1 {
            drawdir = -drawdir;
        }

        unsafe {
            fsa_set_drawdir(drawdir);
            fsa_info_drawdir();
            fsa_neighbours();
            fsa_bind_coords();
        }

        // First nodes on the upwind side, per surface.
        let nofree = unsafe { fsa_nofree() };
        for s in 1..=nofree {
            let mut got = false;
            let mut u1 = 0.0;
            for t in 1..=nn {
                if unsafe { fsa_fsp(t) } != s {
                    continue;
                }
                let u2 = unsafe { fsa_upwind_coord(t) };
                if !got {
                    u1 = u2;
                    got = true;
                }
                if u2 * (drawdir as c_double) < u1 * (drawdir as c_double) {
                    u1 = u2;
                    unsafe { fsa_set_upwind(s, t) };
                }
            }
            unsafe { fsa_info_upwind(s) };
        }

        // The initial displacement is zero.
        unsafe { fsa_zero_move() };
    }

    // Flux integration pass over the free surfaces.
    let nofree = unsafe { fsa_nofree() };
    let nn = unsafe { fsa_num_nodes() };
    let mut norm = 0.0;
    let (mut last_fsum0, mut last_maxerror, mut last_thickmin, mut last_thickmax) =
        (0.0, 0.0, 0.0, 0.0);
    for surface in 1..=nofree {
        let mut done = 0;
        let mut negerror = 0.0;
        let mut poserror = 0.0;
        let mut thickmin = c_double::MAX;
        let mut thickmax = -c_double::MAX;
        let tfirst = unsafe { fsa_get_upwind(surface) };
        let mut fsum0 = 0.0;

        for tt in 0..nn {
            let mut iend = tfirst + tt;
            if iend > nn {
                iend -= nn;
            }
            if unsafe { fsa_fsp(iend) } != surface {
                continue;
            }

            let xold = unsafe { fsa_coord(iend) };
            let mut x1 = unsafe { fsa_coord(iend) };

            // Path from the free surface to the symmetry boundary.
            let mut j = iend;
            let mut istart = 0;
            loop {
                j = unsafe { fsa_left(j) };
                if j == 0 {
                    break;
                }
                let x2 = x1;
                x1 = unsafe { fsa_coord(j) };
                if x1 > x2 {
                    break;
                }
                if unsafe { fsa_move_perm(j) } <= 0 {
                    unsafe { fsa_warn_perm() };
                    break;
                }
                istart = j;
                if unsafe { fsa_fsp(j) } != 0 {
                    break;
                }
            }

            if iend < 1 || iend > nn || istart < 1 || istart > nn {
                unsafe { fsa_warn_limits() };
            }

            let k0 = unsafe { fsa_var_perm(istart) };
            if k0 == 0 {
                unsafe { fsa_warn_noflow() };
            }

            let absdd = unsafe { fsa_drawdir() }.abs();
            let mut v2x = unsafe { fsa_flowk(k0, absdd) };
            let mut v2y = unsafe { fsa_flowk(k0, 3 - absdd) };
            // Set on every pass through the integration loop below (which
            // always executes at least once), mirroring the original's
            // reused local.
            let mut v2: c_double;
            let mut x2 = unsafe { fsa_coord(istart) };
            let mut y2 = unsafe { fsa_coord2(istart) };
            let mut fsum = 0.0;

            // Integrate from the fixed boundary until the desired flux
            // is exceeded.
            let mut i = istart;
            loop {
                let v1x = v2x;
                let v1y = v2y;
                x1 = x2;
                let y1 = y2;
                j = unsafe { fsa_right(i) };
                let k = unsafe { fsa_var_perm(j) };
                if k == 0 {
                    unsafe { fsa_warn_noflow() };
                }
                v2x = unsafe { fsa_flowk(k, absdd) };
                v2y = unsafe { fsa_flowk(k, 3 - absdd) };
                x2 = unsafe { fsa_coord(j) };
                y2 = unsafe { fsa_coord2(j) };

                // Correct the speed in case the elements are not aligned
                // so that dy is zero.
                if (x2 - x1).abs() < 1.0e-10 {
                    unsafe { fsa_warn_dx0() };
                }

                // Correct speeds by taking in account the angle of
                // integration.
                let v1 = v1x - v1y * (y2 - y1) / (x2 - x1);
                v2 = v2x - v2y * (y2 - y1) / (x2 - x1);

                let a = (v1 * x2 - v2 * x1) / (x2 - x1);
                let b = (v2 - v1) / (x2 - x1);

                // Fluxes based on analytical integration over linear
                // elements.
                let mut df = if axisymmetric {
                    a * (x2 * x2 - x1 * x1) / 2.0 + b * (x2 * x2 * x2 - x1 * x1 * x1) / 3.0
                } else {
                    a * (x2 - x1) + b * (x2 * x2 - x1 * x1) / 2.0
                };
                df = df.abs();
                fsum += df;

                if j == iend {
                    break;
                }
                if done > 0 && fsum > fsum0 {
                    break;
                }
                i = j;
            }

            // For the first node set the desired flux; for other nodes
            // move the whole chain to the corrected positions.
            let x0 = unsafe { fsa_coord(istart) };
            let xnew;
            if done == 0 {
                fsum0 = fsum;
                xnew = x2;
                if thickmin > xnew - x0 {
                    thickmin = xnew - x0;
                }
                if thickmax < xnew - x0 {
                    thickmax = xnew - x0;
                }
            } else {
                // The correction divides by the angle-corrected
                // downstream speed `v2` from the last integration step
                // (axisymmetric: `ABS(v2*x2)`).
                let dx = if axisymmetric {
                    (fsum0 - fsum) / (v2 * x2).abs()
                } else {
                    (fsum0 - fsum) / v2.abs()
                };
                xnew = x2 + dx;
                let rel = dx / (xnew - x0);
                if rel > poserror {
                    poserror = rel;
                }
                if rel < negerror {
                    negerror = rel;
                }
                if thickmin > xnew - x0 {
                    thickmin = xnew - x0;
                }
                if thickmax < xnew - x0 {
                    thickmax = xnew - x0;
                }
                if perform_mapping {
                    unsafe { fsa_map_chain(istart, iend, x0, xnew, xold, relax) };
                } else {
                    unsafe { fsa_move_point(iend, xnew, relax) };
                }
            }

            done += 1;
        }

        let maxerror = if -negerror > poserror {
            -negerror
        } else {
            poserror
        };

        let multi = if nofree > 1 { 1 } else { 0 };
        unsafe { fsa_info_surf(surface, multi, fsum0, maxerror, thickmin, thickmax) };

        // MAX(Norm, 1 + maxerror), with the first surface seeding Norm.
        let target = 1.0 + maxerror;
        if surface == 1 || target > norm {
            norm = target;
        }
        last_fsum0 = fsum0;
        last_maxerror = maxerror;
        last_thickmin = thickmin;
        last_thickmax = thickmax;
    }

    // The original stores Norm twice; the net effect is the last
    // surface's value.
    let _ = norm;
    unsafe {
        fsa_set_norm(1.0 + last_maxerror);
        fsa_add_results(last_fsum0, last_maxerror, last_thickmin, last_thickmax);
        fsa_inc_visited();
    }
}
