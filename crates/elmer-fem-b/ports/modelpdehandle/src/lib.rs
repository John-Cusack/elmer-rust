//! Port of `fem/src/modules/ModelPDEhandle.F90` (Elmer @ c0b4bc6) for
//! lane T021-B, at the solver-plugin boundary (whole `.so` → Rust cdylib
//! with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ModelPDEhandle.F90` — one F77-global
//! entry `advdiffsolver_` (lowercase + trailing underscore, as proven by
//! `nm -D` on the reference `.so`) with two `CONTAINS` internals. The
//! entry also exists in `ModelPDE.F90`, but each plugin is its own `.so`
//! replacing its own reference, so no rename.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/modelpdehandle_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`MATMUL`,
//! `SUM`), whole-array updates, and formatted I/O stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//! Faithfulness notes (mirror ModelPDEhandle.F90 statement for statement):
//! - `SAVE` semantics replicated: the two internal procedures' `SAVE
//!   Nodes` live as wrapper-side SAVE state with identical lifetimes
//!   (consumed before the next fill, as in the original); the
//!   host-associated accumulators (`TotArea`, `TotLen`, `TotSrc`), the
//!   `dim` snapshot, and the twelve `ValueHandle_t` element keywords live
//!   as wrapper-side SAVE state too. The per-element `MASS`/`STIFF`/`FORCE`
//!   are automatics, re-created per element call as in the original.
//! - The nonlinear-iteration loop, the `DefaultCutFEM()` re-entry `GOTO 1`
//!   loop, and the `DefaultConverged()` exit match the original exactly.
//! - `dt` is unused by the body (as in the original); `TransientSimulation`
//!   gates `Default1stOrderTime`, exactly as in the original.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn mdh_begin(h_solver: *const c_void);
    fn mdh_init_handles();
    fn mdh_default_start();
    fn mdh_default_initialize();
    fn mdh_maxiter() -> c_int;
    fn mdh_coord_dim() -> c_int;
    fn mdh_reset_totals();
    fn mdh_report_integ();
    fn mdh_default_finish();
    fn mdh_nof_active() -> c_int;
    fn mdh_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn mdh_nof_bnd() -> c_int;
    fn mdh_begin_boundary(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        nb: *mut c_int,
        active: *mut c_int,
    );
    fn mdh_default_cutfem() -> c_int;
    fn mdh_default_finish_bulk();
    fn mdh_default_finish_bnd();
    fn mdh_default_finish_asm();
    fn mdh_default_dirichlet();
    fn mdh_default_solve() -> c_double;
    fn mdh_default_converged() -> c_int;
    fn mdh_lm_nodes();
    fn mdh_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn mdh_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn mdh_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn mdh_lm_ip_values(
        n: c_int,
        basis: *const c_double,
        load: *mut c_double,
        rho: *mut c_double,
        a1: *mut c_double,
        a2: *mut c_double,
        a3: *mut c_double,
        d: *mut c_double,
        c: *mut c_double,
        r: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn mdh_lm_assemble(
        nd: c_int,
        dim: c_int,
        dbasis: *const c_double,
        basis: *const c_double,
        weight: c_double,
        load: c_double,
        d: c_double,
        c: c_double,
        r: c_double,
        rho: c_double,
        a1: c_double,
        a2: c_double,
        a3: c_double,
        mass: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn mdh_time_update(nd: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn mdh_condensate(nd: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn mdh_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn mdh_bc_assoc() -> c_int;
    fn mdh_bc_nodes();
    fn mdh_bc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn mdh_bc_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn mdh_bc_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn mdh_bc_ip_values(
        n: c_int,
        basis: *const c_double,
        f: *mut c_double,
        c: *mut c_double,
        ext: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn mdh_bc_assemble(
        nd: c_int,
        basis: *const c_double,
        weight: c_double,
        c: c_double,
        f: c_double,
        ext: c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
}

// Serial-entry sharing guard (never contended; Elmer calls serially).
static GUARD: Mutex<()> = Mutex::new(());

/// `advdiffsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are unused
/// by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid Fortran `LOGICAL` read (it gates
/// `Default1stOrderTime`, as in the original).
#[export_name = "advdiffsolver_"]
pub unsafe extern "C" fn adv_diff_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    let transient = !transient.is_null() && unsafe { *transient } != 0;
    unsafe { mdh_begin(solver) };

    unsafe { mdh_init_handles() };
    unsafe { mdh_default_start() };

    let maxiter = unsafe { mdh_maxiter() };
    let dim = unsafe { mdh_coord_dim() };

    // Host-associated accumulators zeroed before the nonlinear loop.
    unsafe { mdh_reset_totals() };

    // Nonlinear iteration loop.
    for _iter in 1..=maxiter {
        // System assembly (with the DefaultCutFEM() GOTO 1 re-entry loop).
        unsafe { mdh_default_initialize() };
        loop {
            let active = unsafe { mdh_nof_active() };
            for t in 1..=active {
                let (mut n, mut nd, mut nb) = (0, 0, 0);
                unsafe { mdh_begin_active(t, &mut n, &mut nd, &mut nb) };
                local_matrix(n, nd, nb, dim, transient);
            }
            unsafe { mdh_default_finish_bulk() };

            let nbnd = unsafe { mdh_nof_bnd() };
            for t in 1..=nbnd {
                let (mut n, mut nd, mut nb, mut active_el) = (0, 0, 0, 0);
                unsafe { mdh_begin_boundary(t, &mut n, &mut nd, &mut nb, &mut active_el) };
                if active_el != 0 {
                    local_matrix_bc(nd + nb);
                }
            }

            if unsafe { mdh_default_cutfem() } == 0 {
                break;
            }
        }
        unsafe {
            mdh_default_finish_bnd();
            mdh_default_finish_asm();
            mdh_default_dirichlet();
        }

        // And finally, solve.
        let _norm: c_double = unsafe { mdh_default_solve() };
        if unsafe { mdh_default_converged() } != 0 {
            break;
        }
    }

    unsafe {
        mdh_default_finish();
        mdh_report_integ();
    }
}

/// The `LocalMatrix` internal procedure: bulk element matrix and
/// right-hand side by adaptive numerical integration. `n` is the element
/// node count; `nd`/`nb` split the internal procedure's `nd` dummy
/// (`nd+nb`) into nodal and bubble DOFs for `CondensateP`.
fn local_matrix(n: c_int, nd: c_int, nb: c_int, dim: c_int, transient: bool) {
    let ndnb = nd + nb;
    unsafe { mdh_lm_nodes() };
    let ndn = ndnb as usize;
    // MASS(nd,nd) / STIFF(nd,nd) / FORCE(nd) automatics, zeroed on entry.
    let mut mass = vec![0.0 as c_double; ndn * ndn];
    let mut stiff = vec![0.0 as c_double; ndn * ndn];
    let mut force = vec![0.0 as c_double; ndn];
    // Basis(nd) / dBasisdx(nd,3) automatics, sized by this element.
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { mdh_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { mdh_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            mdh_lm_eleminfo(
                u,
                v,
                w,
                n,
                ndnb,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            )
        };
        let _ = stat;

        let (mut load, mut rho, mut a1, mut a2, mut a3, mut d, mut c, mut r) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        unsafe {
            mdh_lm_ip_values(
                ndnb,
                basis.as_ptr(),
                &mut load,
                &mut rho,
                &mut a1,
                &mut a2,
                &mut a3,
                &mut d,
                &mut c,
                &mut r,
            )
        };

        unsafe {
            mdh_lm_assemble(
                ndnb,
                dim,
                dbasis.as_ptr(),
                basis.as_ptr(),
                s * detj,
                load,
                d,
                c,
                r,
                rho,
                a1,
                a2,
                a3,
                mass.as_mut_ptr(),
                stiff.as_mut_ptr(),
                force.as_mut_ptr(),
            )
        };
    }

    if transient {
        unsafe { mdh_time_update(ndnb, mass.as_mut_ptr(), stiff.as_mut_ptr(), force.as_mut_ptr()) };
    }
    unsafe {
        mdh_condensate(nd, nb, ndnb, stiff.as_mut_ptr(), force.as_mut_ptr());
        mdh_default_update(ndnb, stiff.as_mut_ptr(), force.as_mut_ptr());
    }
}

/// The `LocalMatrixBC` internal procedure: Neumann/Robin boundary
/// contribution by boundary integration. `ndnb` is `nd+nb`; unlike the
/// bulk path there is no bubble condensation here, exactly as in the
/// original.
fn local_matrix_bc(ndnb: c_int) {
    if unsafe { mdh_bc_assoc() } == 0 {
        return;
    }
    unsafe { mdh_bc_nodes() };
    let ndn = ndnb as usize;
    let mut stiff = vec![0.0 as c_double; ndn * ndn];
    let mut force = vec![0.0 as c_double; ndn];
    let mut basis = vec![0.0 as c_double; ndn];

    let nip = unsafe { mdh_bc_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { mdh_bc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe { mdh_bc_eleminfo(u, v, w, ndnb, basis.as_mut_ptr(), &mut detj, &mut stat) };
        let _ = stat;

        let (mut f, mut c, mut ext) = (0.0, 0.0, 0.0);
        unsafe { mdh_bc_ip_values(ndnb, basis.as_ptr(), &mut f, &mut c, &mut ext) };
        unsafe {
            mdh_bc_assemble(
                ndnb,
                basis.as_ptr(),
                s * detj,
                c,
                f,
                ext,
                stiff.as_mut_ptr(),
                force.as_mut_ptr(),
            )
        };
    }
    unsafe { mdh_default_update(ndnb, stiff.as_mut_ptr(), force.as_mut_ptr()) };
}
