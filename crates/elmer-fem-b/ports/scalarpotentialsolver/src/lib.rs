//! Port of `fem/src/modules/ScalarPotentialSolver.F90` (the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ScalarPotentialSolver.F90` — two F77
//! globals, `scalarpotentialsolver_` with two `CONTAINS` internals and
//! `scalarpotentialsolver_init_`. The emitted dynamic symbols are
//! exactly `scalarpotentialsolver_` and `scalarpotentialsolver_init_`
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/scalarpotentialsolver_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`) and formatted I/O (`WRITE`/`Info`) stay in Fortran wrappers,
//! so bit-identity holds by construction wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors ScalarPotentialSolver.F90 statement for
//! statement):
//! - The header `Info` lines run before the `Matrix`/`Perm` early
//!   returns, exactly as in the original.
//! - `SAVE` semantics replicated: the `Nodes` SAVE state and the
//!   `Conductivity` SAVE pointer of the `BulkAssembly` internal
//!   procedure live as wrapper-side SAVE state with identical lifetimes
//!   (re-pointed/re-filled per element, as in the original).
//! - The entry's dead locals (`Visited`, `Unorm`, `TotNorm` aside from
//!   its message use, `ForceVector`, `SaveRHS`, `DOFs`, `ActiveDir`,
//!   `DirMask`, `ScalarPotentialSol`) carry no observable behavior and
//!   are intentionally unmirrored; `TotNorm` is kept through the result
//!   message.
//! - Assembly scratch (`STIFF`, `FORCE`, `VectorField`, `Basis`,
//!   `dBasisdx`) is allocated per solver entry sized by the current
//!   mesh's `MAX(MaxElementDOFs, MaxElementNodes)`, exactly like the
//!   original's per-`BulkAssembly`-call `ALLOCATE`s.
//! - Allocation failure calls the same `Fatal` text family via
//!   `try_reserve` instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies,
//!   exactly as in the original (received and ignored).
//! - The init entry's `dim` touch call and the Laplace vs. coefficient
//!   vs. constant-bulk branches are preserved statement for statement.
//! - No Rust threads: the solver calls entries serially, so no `static`
//!   state exists at all on the Rust side.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn sps_init();
    fn sps_begin(h_solver: *const c_void);
    fn sps_assoc_matrix() -> c_int;
    fn sps_active_dofs() -> c_int;
    fn sps_fetch();
    fn sps_dim() -> c_int;
    fn sps_inuse() -> c_int;
    fn sps_max_elem() -> c_int;
    fn sps_nof_active() -> c_int;
    fn sps_info_header();
    fn sps_default_initialize();
    fn sps_finish_bulk();
    fn sps_finish();
    fn sps_dirichlet();
    fn sps_solve() -> c_double;
    fn sps_realtime() -> c_double;
    fn sps_info_assembly(t0: c_double, t1: c_double);
    fn sps_info_solution(t1: c_double, t2: c_double);
    fn sps_info_result(norm: c_double);
    fn sps_begin_elem(t: c_int, nd: *mut c_int, n: *mut c_int);
    fn sps_fetch_cond(got: *mut c_int, rank: *mut c_int);
    fn sps_vector_field(ndim: c_int, ld: c_int, vf: *mut c_double);
    fn sps_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn sps_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn sps_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        ld: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
    );
    fn sps_weight(
        s: c_double,
        detj: c_double,
        nelem: c_int,
        ld: c_int,
        basis: *const c_double,
    ) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn sps_stiff_laplace(
        nd: c_int,
        ndim: c_int,
        ld: c_int,
        dbasis: *mut c_double,
        weight: c_double,
        stiff: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn sps_stiff_cond(
        rank: c_int,
        nd: c_int,
        nelem: c_int,
        ndim: c_int,
        ld: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        weight: c_double,
        stiff: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn sps_force_update(
        nd: c_int,
        nelem: c_int,
        ndim: c_int,
        ld: c_int,
        vf: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        weight: c_double,
        force: *mut c_double,
    );
    fn sps_update_eq(ld: c_int, nd: c_int, stiff: *mut c_double, force: *mut c_double);
    fn sps_update_force(nd: c_int, force: *mut c_double);
    fn sps_fatal_alloc();
}

/// Grow-or-create helper for the per-entry scratch arrays. Allocation
/// failure reports through the solver's own `Fatal`, mirroring the
/// original's allocation checks, and never returns.
fn reserve(buf: &mut Vec<c_double>, n: usize) {
    if buf.try_reserve(n.saturating_sub(buf.capacity())).is_err() {
        unsafe { sps_fatal_alloc() };
        std::process::abort();
    }
    buf.resize(n, 0.0);
}

/// `scalarpotentialsolver_init_` — the init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; all four dummies are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// The addresses only travel opaquely into the wrapper layer.
#[export_name = "scalarpotentialsolver_init_"]
pub unsafe extern "C" fn scalar_potential_solver_init(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sps_init() };
}

/// `scalarpotentialsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "scalarpotentialsolver_"]
pub unsafe extern "C" fn scalar_potential_solver(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sps_begin(solver) };

    // Header first: the original prints these before the early returns.
    unsafe { sps_info_header() };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { sps_assoc_matrix() } == 0 {
        return;
    }
    // IF (COUNT(Solver % Variable % Perm > 0) <= 0) RETURN
    if unsafe { sps_active_dofs() } == 0 {
        return;
    }

    unsafe { sps_fetch() };
    let dim = unsafe { sps_dim() };
    let inuse = unsafe { sps_inuse() };

    let t0: c_double = unsafe { sps_realtime() };

    // Per-entry assembly scratch, sized by the current mesh exactly like
    // the original's per-call ALLOCATEs.
    let n = unsafe { sps_max_elem() } as usize;
    let ld = n as c_int;
    let mut stiff = Vec::<c_double>::new();
    let mut force = Vec::<c_double>::new();
    let mut vf = Vec::<c_double>::new();
    let mut basis = Vec::<c_double>::new();
    let mut dbasis = Vec::<c_double>::new();
    reserve(&mut stiff, n * n);
    reserve(&mut force, n);
    reserve(&mut vf, 3 * n);
    reserve(&mut basis, n);
    reserve(&mut dbasis, n * 3);

    unsafe { sps_default_initialize() };

    bulk_assembly(
        &mut stiff,
        &mut force,
        &mut vf,
        &mut basis,
        &mut dbasis,
        dim,
        inuse,
        ld,
    );

    unsafe {
        sps_finish_bulk();
        sps_finish();
        sps_dirichlet();
    }

    let t1: c_double = unsafe { sps_realtime() };
    unsafe { sps_info_assembly(t0, t1) };

    let totnorm: c_double = unsafe { sps_solve() };

    let t2: c_double = unsafe { sps_realtime() };
    unsafe {
        sps_info_solution(t1, t2);
        sps_info_result(totnorm);
    }
}

/// The `BulkAssembly` internal procedure: element loop, numerical
/// integration, and the global-matrix/force update per element.
#[allow(clippy::too_many_arguments)]
fn bulk_assembly(
    stiff: &mut Vec<c_double>,
    force: &mut Vec<c_double>,
    vf: &mut Vec<c_double>,
    basis: &mut Vec<c_double>,
    dbasis: &mut Vec<c_double>,
    dim: c_int,
    inuse: c_int,
    ld: c_int,
) {
    let active = unsafe { sps_nof_active() };
    for t in 1..=active {
        let (mut nd, mut ne) = (0, 0);
        unsafe { sps_begin_elem(t, &mut nd, &mut ne) };

        let (mut got, mut rank) = (0, 0);
        unsafe { sps_fetch_cond(&mut got, &mut rank) };

        unsafe { sps_vector_field(dim, ld, vf.as_mut_ptr()) };

        stiff.fill(0.0);
        force.fill(0.0);

        let nip = unsafe { sps_ip_begin() };
        for q in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { sps_ip_point(q, &mut u, &mut v, &mut w, &mut s) };

            let mut detj = 0.0;
            unsafe {
                sps_eleminfo(
                    u,
                    v,
                    w,
                    ld,
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    &mut detj,
                )
            };

            let weight: c_double =
                unsafe { sps_weight(s, detj, ne, ld, basis.as_ptr()) };

            if inuse == 0 {
                if got != 0 {
                    unsafe {
                        sps_stiff_cond(
                            rank,
                            nd,
                            ne,
                            dim,
                            ld,
                            basis.as_mut_ptr(),
                            dbasis.as_mut_ptr(),
                            weight,
                            stiff.as_mut_ptr(),
                        )
                    };
                } else {
                    unsafe {
                        sps_stiff_laplace(
                            nd,
                            dim,
                            ld,
                            dbasis.as_mut_ptr(),
                            weight,
                            stiff.as_mut_ptr(),
                        )
                    };
                }
            }

            unsafe {
                sps_force_update(
                    nd,
                    ne,
                    dim,
                    ld,
                    vf.as_mut_ptr(),
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    weight,
                    force.as_mut_ptr(),
                )
            };
        }

        if inuse == 0 {
            unsafe { sps_update_eq(ld, nd, stiff.as_mut_ptr(), force.as_mut_ptr()) };
        } else {
            unsafe { sps_update_force(nd, force.as_mut_ptr()) };
        }
    }
}
