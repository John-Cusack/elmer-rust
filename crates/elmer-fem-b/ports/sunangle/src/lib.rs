//! Port of `fem/src/modules/SunAngle.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/SunAngle.F90` — two F77-global entries
//! `sunanglesolver_` and `sunanglesolver_init_` (proven by `nm -D` on the
//! reference `.so`); the `CONTAINS` internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque `*mut c_void` handles. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/sunangle_acc.F90`, compiled against the built `.mod` files of
//! the target configure — the port itself contains no layout assumption and
//! no pinned offset. Order-sensitive whole-array fills (`SunAngle = ...`,
//! `WHERE`), formatted I/O (`Info` transcripts, the `PRINT *` ridge tally)
//! and the MPI ridge exchange stay in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors SunAngle.F90 statement for statement):
//! - `SAVE` semantics replicated explicitly: the `Visited` latch and the
//!   `RidgeData` allocation live in the `static` `STATE` with identical
//!   lifetimes (populated on the first parallel call, released when
//!   `dofs > 1`, `Fatal` on release-when-unallocated exactly like
//!   `ReleaseRidge`); the `SetTestHeight` element loop runs whole in a
//!   Fortran wrapper with the original's evaluation order.
//! - `eps = 1.0e-6` is the `REAL(4)` literal converted to `REAL(8)`,
//!   replicated as `1e-6f32 as f64`; `EPSILON(detA)` is `f64::EPSILON`.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored).
//! - No Rust threads: the solver entries run serially, so the `Mutex` is
//!   only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sua_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sua_init(h_solver: *const c_void);
    fn sua_info(which: c_int);
    fn sua_info_active20() -> c_int;
    fn sua_print_totals(mype: c_int, k: c_int, l: c_int, nmax: c_int);
    fn sua_fatal_release();
    fn sua_dim() -> c_int;
    fn sua_nnodes() -> c_int;
    fn sua_nedges() -> c_int;
    fn sua_dofs() -> c_int;
    fn sua_nvals() -> c_int;
    fn sua_nperm() -> c_int;
    fn sua_edges_associated() -> c_int;
    fn sua_ensure_edges(dim: c_int);
    fn sua_edge_nodes(t: c_int, i0: *mut c_int, i1: *mut c_int);
    fn sua_fill_x(buf: *mut c_double, n: c_int);
    fn sua_fill_y(buf: *mut c_double, n: c_int);
    fn sua_fill_z(buf: *mut c_double, n: c_int);
    fn sua_fill_perm(buf: *mut c_int, n: c_int);
    fn sua_fill_sun(buf: *mut c_double, n: c_int);
    fn sua_sun_fill(v: c_double);
    fn sua_sun_store(buf: *const c_double, n: c_int);
    fn sua_where_fix();
    fn sua_test_height_present() -> c_int;
    fn sua_set_test_height();
    fn sua_enforce(which: c_int, val: *mut c_int, found: *mut c_int);
    fn sua_pes() -> c_int;
    fn sua_mype() -> c_int;
    fn sua_maxdist(val: *mut c_double, got: *mut c_int);
    fn sua_test_angle() -> c_double;
    fn sua_mpi_allreduce_bb(buf: *mut c_double, pes: c_int);
    fn sua_mpi_barrier();
    fn sua_mpi_send_count(n: c_int, proc: c_int);
    fn sua_mpi_recv_count(n: *mut c_int, proc: c_int);
    fn sua_mpi_bsend(buf: *mut c_double, n: c_int, proc: c_int);
    fn sua_mpi_recv(buf: *mut c_double, n: c_int, proc: c_int);
}

/// One partition's ridge share: mirrors `RidgeData_t` (`n_send`, `n_recv`,
/// `Coords`, `MeshBB`) for all PEs; index 0 == `MyPe == 1`.
struct Ridge {
    n_send: Vec<c_int>,
    n_recv: Vec<c_int>,
    coords: Vec<Vec<c_double>>,
    meshbb: Vec<[c_double; 4]>,
}

/// First-call SAVE state: the `Visited` latch and the `RidgeData`
/// allocation, with the original's exact lifetimes.
struct SunState {
    visited: bool,
    ridge: Option<Ridge>,
}

impl SunState {
    const fn new() -> Self {
        Self {
            visited: false,
            ridge: None,
        }
    }
}

static STATE: Mutex<SunState> = Mutex::new(SunState::new());

/// `sunanglesolver_init_` — declares the `CriticalTan` variable and the
/// `No Matrix` flag, exactly as the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "sunanglesolver_init_"]
pub unsafe extern "C" fn sun_angle_solver_init(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sua_init(solver) };
}

/// `sunanglesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` may be null (they are ignored).
#[export_name = "sunanglesolver_"]
pub unsafe extern "C" fn sun_angle_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sua_begin(model, solver) };

    unsafe { sua_info(1) };
    unsafe { sua_info(2) };
    unsafe { sua_info(3) };

    let dim = unsafe { sua_dim() };

    // Find mesh edges in order to define the intersection points.
    if unsafe { sua_edges_associated() } == 0 {
        unsafe { sua_info(6) };
        unsafe { sua_ensure_edges(dim) };
    }

    // For testing purposes set the height here.
    if unsafe { sua_test_height_present() } != 0 {
        unsafe { sua_set_test_height() };
    }

    // Snapshot the mesh fields the search loops read (the original binds
    // Fortran pointers before this point; the values are identical after
    // SetTestHeight has run).
    let nn = unsafe { sua_nnodes() } as usize;
    let ne = unsafe { sua_nedges() } as usize;
    let dofs = unsafe { sua_dofs() } as usize;
    let nperm = unsafe { sua_nperm() } as usize;
    let nvals = unsafe { sua_nvals() } as usize;
    let mut x = vec![0.0 as c_double; nn];
    let mut y = vec![0.0 as c_double; nn];
    let mut z = vec![0.0 as c_double; nn];
    let mut perm = vec![0 as c_int; nperm];
    let mut sun = vec![0.0 as c_double; nvals];
    unsafe {
        sua_fill_x(x.as_mut_ptr(), nn as c_int);
        sua_fill_y(y.as_mut_ptr(), nn as c_int);
        sua_fill_z(z.as_mut_ptr(), nn as c_int);
        sua_fill_perm(perm.as_mut_ptr(), nperm as c_int);
        sua_fill_sun(sun.as_mut_ptr(), nvals as c_int);
    }

    // Enforce use of the same algo in serial as is used in parallel.
    let (mut enf_par, mut enf_ser) = (0, 0);
    let mut found = 0;
    unsafe {
        sua_enforce(1, &mut enf_par, &mut found);
        sua_enforce(2, &mut enf_ser, &mut found);
    }
    let _ = found;

    // Even in parallel do not do any communication.
    let pes = unsafe { sua_pes() };
    let parallel = enf_par != 0 || pes > 1;
    if parallel && !state.visited {
        let (mut maxdist, mut got) = (0.0, 0);
        unsafe { sua_maxdist(&mut maxdist, &mut got) };
        unsafe { sua_info(7) };
        populate_ridge(&x, &y, &z, ne, enf_par != 0, got != 0, maxdist, pes, &mut state.ridge);
    }

    unsafe { sua_info(4) };

    // eps = 1.0e-6: single-precision literal converted to double.
    let eps: c_double = 1e-6f32 as c_double;
    unsafe { sua_sun_fill(-100.0) };
    sun.fill(-100.0);

    let mut amat = [[0.0 as c_double; 2]; 2];
    for idof in 1..=dofs {
        let alpha: c_double = if dofs == 1 {
            unsafe { sua_test_angle() }
        } else {
            (idof - 1) as c_double * 2.0 * std::f64::consts::PI / dofs as c_double
        };
        amat[0][0] = alpha.cos();
        amat[1][0] = alpha.sin();

        // We use different search algo locally unless not otherwise requested.
        if enf_par == 0 {
            set_shading_angle(&x, &y, &z, &perm, &mut sun, ne, dofs, idof, &mut amat, eps);
        }

        if parallel {
            set_shading_angle_parallel(
                &x,
                &y,
                &z,
                &perm,
                &mut sun,
                dofs,
                idof,
                &mut amat,
                eps,
                enf_ser != 0,
                state.ridge.as_ref(),
            );
        }
    }

    // Set default value to unset nodes (=boundary at sunside).
    unsafe { sua_sun_store(sun.as_ptr(), nvals as c_int) };
    unsafe { sua_where_fix() };

    state.visited = true;

    // Heuristics: release the ridge share once many angles are done.
    if dofs > 1 {
        if state.ridge.is_none() {
            unsafe { sua_fatal_release() };
            std::process::abort();
        }
        state.ridge = None;
    }

    unsafe { sua_info(5) };
}

/// Solves a 2x2 linear system; mirrors the `Solve2x2` internal procedure
/// including the `EPSILON` singularity test and evaluation order.
fn solve2x2(a: &[[c_double; 2]; 2], b1: c_double, b2: c_double) -> (c_double, c_double, bool) {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let singular = det.abs() < c_double::EPSILON;
    if singular {
        return (0.0, 0.0, true);
    }
    let inv = 1.0 / det;
    let c1 = inv * (a[1][1] * b1 - a[0][1] * b2);
    let c2 = inv * (a[0][0] * b2 - a[1][0] * b1);
    (c1, c2, false)
}

/// The serial shading search: the stupid N^2 double loop over edges and
/// nodes, in the original's order (edge outer, node inner).
#[allow(clippy::too_many_arguments)]
fn set_shading_angle(
    x: &[c_double],
    y: &[c_double],
    z: &[c_double],
    perm: &[c_int],
    sun: &mut [c_double],
    ne: usize,
    dofs: usize,
    idof: usize,
    amat: &mut [[c_double; 2]; 2],
    eps: c_double,
) {
    for t in 1..=ne {
        let (mut i0, mut i1) = (0, 0);
        unsafe { sua_edge_nodes(t as c_int, &mut i0, &mut i1) };
        if perm[(i0 - 1) as usize] == 0 || perm[(i1 - 1) as usize] == 0 {
            continue;
        }
        let x0 = [x[(i0 - 1) as usize], x[(i1 - 1) as usize]];
        let y0 = [y[(i0 - 1) as usize], y[(i1 - 1) as usize]];
        let z0 = [z[(i0 - 1) as usize], z[(i1 - 1) as usize]];
        let dz0 = z0[1] - z0[0];

        amat[0][1] = x0[0] - x0[1];
        amat[1][1] = y0[0] - y0[1];

        for (j0, node) in x.iter().enumerate() {
            let j = j0 + 1;
            let k = perm[j0];
            if k == 0 {
                continue;
            }
            // Nodes that are part of the edge cannot shade themselves.
            if j as c_int == i0 || j as c_int == i1 {
                continue;
            }
            let b1 = x0[0] - *node;
            let b2 = y0[0] - y[j0];
            let (c1, c2, singular) = solve2x2(amat, b1, b2);
            if singular || c1 < eps || c2 < -eps || c2 > 1.0 + eps {
                continue;
            }
            let dz = (z0[0] + c2 * dz0) - z[j0];
            let tanphi = dz / c1;
            let idx = dofs * ((k - 1) as usize) + (idof - 1);
            if tanphi > sun[idx] {
                sun[idx] = tanphi;
            }
        }
    }
}

/// Populates the ridge shares; mirrors `PopulateRidge`, including the
/// serial-consistency fill under `Enforce Parallel` and the guarded
/// PEs > 1 exchange (which never runs in the serial tree).
#[allow(clippy::too_many_arguments)]
fn populate_ridge(
    x: &[c_double],
    y: &[c_double],
    z: &[c_double],
    ne: usize,
    enforce_parallel: bool,
    got_maxdist: bool,
    maxdist: c_double,
    pes: c_int,
    slot: &mut Option<Ridge>,
) {
    let npes = pes as usize;
    // MyPe = ParEnv % MyPe + 1, converted to 0-based.
    let my = unsafe { sua_mype() } as usize;
    let mut ridge = Ridge {
        n_send: vec![0; npes],
        n_recv: vec![0; npes],
        coords: (0..npes).map(|_| Vec::new()).collect(),
        meshbb: vec![[0.0; 4]; npes],
    };

    // Bounding box for the quick and dirty search.
    if got_maxdist {
        let (mut xmin, mut xmax) = (x[0], x[0]);
        let (mut ymin, mut ymax) = (y[0], y[0]);
        for (&vx, &vy) in x.iter().zip(y.iter()) {
            if vx < xmin {
                xmin = vx;
            }
            if vx > xmax {
                xmax = vx;
            }
            if vy < ymin {
                ymin = vy;
            }
            if vy > ymax {
                ymax = vy;
            }
        }
        ridge.meshbb[my] = [xmin, xmax, ymin, ymax];
    }

    // Serial/parallel consistency check fill.
    if enforce_parallel {
        ridge.n_recv[my] = ne as c_int;
        let mut coords = vec![0.0 as c_double; 6 * ne];
        for k in 1..=ne {
            let (mut i0, mut i1) = (0, 0);
            unsafe { sua_edge_nodes(k as c_int, &mut i0, &mut i1) };
            let l = k - 1;
            coords[6 * l] = x[(i0 - 1) as usize];
            coords[6 * l + 1] = y[(i0 - 1) as usize];
            coords[6 * l + 2] = z[(i0 - 1) as usize];
            coords[6 * l + 3] = x[(i0 - 1) as usize] - x[(i1 - 1) as usize];
            coords[6 * l + 4] = y[(i0 - 1) as usize] - y[(i1 - 1) as usize];
            coords[6 * l + 5] = z[(i0 - 1) as usize] - z[(i1 - 1) as usize];
        }
        ridge.coords[my] = coords;
    }

    if npes > 1 {
        if got_maxdist {
            let mut rpar = vec![-c_double::MAX; 4 * npes];
            rpar[4 * my..4 * my + 4].copy_from_slice(&ridge.meshbb[my]);
            unsafe { sua_mpi_allreduce_bb(rpar.as_mut_ptr(), pes) };
            for i in 0..npes {
                if i == my {
                    continue;
                }
                ridge.meshbb[i] = [rpar[4 * i], rpar[4 * i + 1], rpar[4 * i + 2], rpar[4 * i + 3]];
            }
            unsafe { sua_mpi_barrier() };
        }

        let mut n_max: c_int = 0;
        for i in 0..npes {
            if i == my {
                continue;
            }
            let proc = i as c_int;
            let n_send: c_int = send_count(&ridge.meshbb[i], x, y, ne, got_maxdist, maxdist);
            ridge.n_send[i] = n_send;
            if n_send > n_max {
                n_max = n_send;
            }
            unsafe { sua_mpi_send_count(n_send, proc) };
        }

        for i in 0..npes {
            if i == my {
                continue;
            }
            let mut n_recv: c_int = 0;
            unsafe { sua_mpi_recv_count(&mut n_recv, i as c_int) };
            ridge.n_recv[i] = n_recv;
            if n_recv > 0 {
                ridge.coords[i] = vec![0.0; 6 * n_recv as usize];
            }
        }
        unsafe { sua_mpi_barrier() };

        let ksum: c_int = ridge.n_send.iter().sum();
        let lsum: c_int = ridge.n_recv.iter().sum();
        if unsafe { sua_info_active20() } != 0 {
            unsafe { sua_print_totals(my as c_int + 1, ksum, lsum, n_max) };
        }

        for i in 0..npes {
            if i == my {
                continue;
            }
            let n_send = ridge.n_send[i];
            if n_send > 0 {
                let mut buf = send_fill(&ridge.meshbb[i], x, y, z, ne, got_maxdist, maxdist);
                debug_assert_eq!(buf.len(), 6 * n_send as usize);
                unsafe { sua_mpi_bsend(buf.as_mut_ptr(), 6 * n_send, i as c_int) };
            }
        }
        for i in 0..npes {
            if i == my {
                continue;
            }
            let n_recv = ridge.n_recv[i];
            if n_recv > 0 {
                unsafe {
                    sua_mpi_recv(ridge.coords[i].as_mut_ptr(), 6 * n_recv, i as c_int);
                }
            }
        }
        unsafe { sua_mpi_barrier() };
    }
    *slot = Some(ridge);
}

/// Counts the edges this partition must send to the partition owning `bb`;
/// mirrors the count loop in `PopulateRidge` (first-node filter).
fn send_count(
    bb: &[c_double; 4],
    x: &[c_double],
    y: &[c_double],
    ne: usize,
    got_maxdist: bool,
    maxdist: c_double,
) -> c_int {
    if !got_maxdist {
        return ne as c_int;
    }
    let mut count: c_int = 0;
    for k in 1..=ne {
        let (mut i0, mut i1) = (0, 0);
        unsafe { sua_edge_nodes(k as c_int, &mut i0, &mut i1) };
        let _ = i1;
        let px = x[(i0 - 1) as usize];
        let py = y[(i0 - 1) as usize];
        if bb[0] > px + maxdist {
            continue;
        }
        if px > bb[1] + maxdist {
            continue;
        }
        if bb[2] > py + maxdist {
            continue;
        }
        if py > bb[3] + maxdist {
            continue;
        }
        count += 1;
    }
    count
}

/// Tabulates the `(x0, y0, z0, dx, dy, dz)` edge payload for one receiver;
/// mirrors the send loop in `PopulateRidge`.
fn send_fill(
    bb: &[c_double; 4],
    x: &[c_double],
    y: &[c_double],
    z: &[c_double],
    ne: usize,
    got_maxdist: bool,
    maxdist: c_double,
) -> Vec<c_double> {
    let mut buf = Vec::new();
    for k in 1..=ne {
        let (mut i0, mut i1) = (0, 0);
        unsafe { sua_edge_nodes(k as c_int, &mut i0, &mut i1) };
        if got_maxdist {
            let px = x[(i0 - 1) as usize];
            let py = y[(i0 - 1) as usize];
            if bb[0] > px + maxdist {
                continue;
            }
            if px > bb[1] + maxdist {
                continue;
            }
            if bb[2] > py + maxdist {
                continue;
            }
            if py > bb[3] + maxdist {
                continue;
            }
        }
        buf.push(x[(i0 - 1) as usize]);
        buf.push(y[(i0 - 1) as usize]);
        buf.push(z[(i0 - 1) as usize]);
        buf.push(x[(i0 - 1) as usize] - x[(i1 - 1) as usize]);
        buf.push(y[(i0 - 1) as usize] - y[(i1 - 1) as usize]);
        buf.push(z[(i0 - 1) as usize] - z[(i1 - 1) as usize]);
    }
    buf
}

/// The parallel shading search over the communicated ridge shares; mirrors
/// `SetShadingAngleParallel` including the serial-seeded `maxangle`.
#[allow(clippy::too_many_arguments)]
fn set_shading_angle_parallel(
    x: &[c_double],
    y: &[c_double],
    z: &[c_double],
    perm: &[c_int],
    sun: &mut [c_double],
    dofs: usize,
    idof: usize,
    amat: &mut [[c_double; 2]; 2],
    eps: c_double,
    enforce_serial: bool,
    ridge: Option<&Ridge>,
) {
    let ridge = match ridge {
        Some(r) => r,
        None => return,
    };
    let npes = ridge.n_recv.len();
    let my = unsafe { sua_mype() } as usize;
    for (j0, (&xj, &yj)) in x.iter().zip(y.iter()).enumerate() {
        let k = perm[j0];
        if k == 0 {
            continue;
        }
        let z1 = z[j0];
        let idx = dofs * ((k - 1) as usize) + (idof - 1);
        let mut maxangle = sun[idx];
        for i in 0..npes {
            if ridge.n_recv[i] == 0 {
                continue;
            }
            if enforce_serial && i != my {
                continue;
            }
            let c = &ridge.coords[i];
            for t0 in 0..ridge.n_recv[i] as usize {
                let b = 6 * t0;
                amat[0][1] = c[b + 3];
                amat[1][1] = c[b + 4];
                let b1 = c[b] - xj;
                let b2 = c[b + 1] - yj;
                let (c1, c2, singular) = solve2x2(amat, b1, b2);
                if singular || c1 < eps || c2 < -eps || c2 > 1.0 + eps {
                    continue;
                }
                let dz = (c[b + 2] - c2 * c[b + 5]) - z1;
                let v = dz / c1;
                if v > maxangle {
                    maxangle = v;
                }
            }
        }
        sun[idx] = maxangle;
    }
}
