//! Port of `fem/src/modules/ModelPDE.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ModelPDE.F90` — one F77-global entry
//! `advdiffsolver_` with two `CONTAINS` internals (`LocalMatrix`,
//! `LocalMatrixBC`). The emitted dynamic symbol is exactly
//! `advdiffsolver_` (proven by `nm -D` on the reference `.so`); the
//! internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/modelpde_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`SUM`, `MATMUL`), the `Info` call, and the
//! `TotArea`/`TotLen`/`TotSrc` accumulation stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ModelPDE.F90 statement for statement):
//! - `SAVE` semantics replicated: each internal procedure's `SAVE Nodes`
//!   lives as wrapper-side SAVE state with identical lifetimes (consumed
//!   before the next fill, as in the original). The internal procedures'
//!   automatic arrays (`STIFF`, `MASS`, `FORCE`, `LOAD`, coefficient
//!   arrays) are fresh-zeroed per element call in the original; here they
//!   are grow-only scratch buffers, zeroed over the used leading section
//!   before each fill, so no stale tail is ever read.
//! - `TransientSimulation` (a `LOGICAL` entry dummy, 4-byte gfortran
//!   convention) selects the `Default1stOrderTime` branch, exactly as the
//!   original's `IF(TransientSimulation)`; `Model` and `dt` are unused by
//!   the body, exactly as in the original (received and ignored).
//! - `DefaultCutFEM()` re-assembly (`GOTO 1`) is a Rust loop over the same
//!   assembly block: re-entry re-reads `GetNOFActive()` without
//!   re-running `DefaultInitialize()` or zeroing the totals.
//! - `Norm = DefaultSolve()` is kept and discarded; the `DefaultConverged()`
//!   exit, the `maxiter` default (`1` when not `Found`, `minv=1`), and the
//!   trailing `res: integ ...` report (both `ListGetLogical` arms evaluated,
//!   Fortran `.OR.` semantics) live in the wrappers.
//! - `ElementInfo`'s `Stat` and every `GetReal`'s `Found` are received and
//!   ignored, exactly as in the original (never tested there).
//! - `LocalMatrixBC`'s `LOAD` array is zeroed but never read in the
//!   original; the zeroing is dropped as unobservable (documented here).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn moa_begin(h_solver: *const c_void);
    fn moa_default_start();
    fn moa_max_iter() -> c_int;
    fn moa_nof_active() -> c_int;
    fn moa_default_initialize();
    fn moa_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn moa_finish_bulk();
    fn moa_num_bnd() -> c_int;
    fn moa_begin_boundary(t: c_int, n: *mut c_int, nd: *mut c_int, active: *mut c_int);
    fn moa_cut_fem() -> c_int;
    fn moa_finish_boundary();
    fn moa_finish_assembly();
    fn moa_dirichlet();
    fn moa_solve() -> c_double;
    fn moa_converged() -> c_int;
    fn moa_default_finish();
    fn moa_report(totarea: c_double, totlen: c_double, totsrc: c_double);
    fn moa_lm_nodes();
    fn moa_bodyforce() -> c_int;
    fn moa_load_source(n: c_int, load: *mut c_double, found: *mut c_int);
    fn moa_mat_begin();
    fn moa_mat_coeff(which: c_int, n: c_int, out: *mut c_double, found: *mut c_int);
    fn moa_conv_vel(i: c_int, n: c_int, velo: *mut c_double);
    fn moa_coord_dim() -> c_int;
    fn moa_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn moa_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn moa_lm_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn moa_lm_ip(
        n: c_int, nd: c_int, dim: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        load: *mut c_double, timec: *mut c_double, velo: *mut c_double,
        diffc: *mut c_double, convc: *mut c_double, reactc: *mut c_double,
        weight: c_double,
        stiff: *mut c_double, mass: *mut c_double, force: *mut c_double,
        loadatip: *mut c_double, totarea: *mut c_double, totsrc: *mut c_double,
    );
    fn moa_lm_transient(nd: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn moa_condensate(nd: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn moa_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn moa_bc_begin() -> c_int;
    fn moa_bc_nodes();
    fn moa_bc_coeffs(n: c_int, flux: *mut c_double, coeff: *mut c_double, ext: *mut c_double);
    fn moa_bc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn moa_bc_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn moa_bc_eleminfo(
        u: c_double, v: c_double, w: c_double,
        n: c_int, nd: c_int,
        basis: *mut c_double, dbasis: *mut c_double,
        detj: *mut c_double, stat: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn moa_bc_ip(
        n: c_int, nd: c_int,
        basis: *mut c_double,
        flux: *mut c_double, coeff: *mut c_double, extt: *mut c_double,
        weight: c_double,
        stiff: *mut c_double, force: *mut c_double,
        totlen: *mut c_double,
    );
    fn moa_fatal_alloc();
}

/// Grow-only scratch for the internal procedures' automatic arrays.
/// `maxn`/`maxnd` are high-water marks; each element call uses the tight
/// leading `(n)` / `(full, full)` view at the buffer start, zeroed before
/// the fill — matching fresh automatic arrays.
struct Scratch {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    diffc: Vec<c_double>,
    reactc: Vec<c_double>,
    convc: Vec<c_double>,
    timec: Vec<c_double>,
    velo: Vec<c_double>,
    basis: Vec<c_double>,
    dbasis: Vec<c_double>,
    flux: Vec<c_double>,
    coeff: Vec<c_double>,
    ext: Vec<c_double>,
    maxn: usize,
    maxnd: usize,
}

impl Scratch {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            diffc: Vec::new(),
            reactc: Vec::new(),
            convc: Vec::new(),
            timec: Vec::new(),
            velo: Vec::new(),
            basis: Vec::new(),
            dbasis: Vec::new(),
            flux: Vec::new(),
            coeff: Vec::new(),
            ext: Vec::new(),
            maxn: 0,
            maxnd: 0,
        }
    }
}

static STATE: Mutex<Scratch> = Mutex::new(Scratch::new());

/// Reserve the tight leading views for one element call (`n` nodes,
/// `full` = retained + bubble DOFs). Buffers only grow; the caller zeroes
/// the used sections before each fill.
fn ensure(st: &mut Scratch, need_n: usize, need_nd: usize) {
    if need_n > st.maxn {
        for v in [
            &mut st.load,
            &mut st.diffc,
            &mut st.reactc,
            &mut st.convc,
            &mut st.timec,
            &mut st.flux,
            &mut st.coeff,
            &mut st.ext,
        ] {
            if v.try_reserve(need_n - v.len()).is_err() {
                unsafe { moa_fatal_alloc() };
                std::process::abort();
            }
            v.resize(need_n, 0.0);
        }
        if st.velo.try_reserve(3 * need_n - st.velo.len()).is_err() {
            unsafe { moa_fatal_alloc() };
            std::process::abort();
        }
        st.velo.resize(3 * need_n, 0.0);
        if st.basis.try_reserve(need_n - st.basis.len()).is_err() {
            unsafe { moa_fatal_alloc() };
            std::process::abort();
        }
        st.basis.resize(need_n, 0.0);
        st.maxn = need_n;
    }
    if need_nd > st.maxnd {
        for v in [&mut st.stiff, &mut st.mass] {
            if v.try_reserve(need_nd * need_nd - v.len()).is_err() {
                unsafe { moa_fatal_alloc() };
                std::process::abort();
            }
            v.resize(need_nd * need_nd, 0.0);
        }
        if st.force.try_reserve(need_nd.saturating_sub(st.force.len())).is_err() {
            unsafe { moa_fatal_alloc() };
            std::process::abort();
        }
        st.force.resize(need_nd, 0.0);
        if st.dbasis.try_reserve(3 * need_nd - st.dbasis.len()).is_err() {
            unsafe { moa_fatal_alloc() };
            std::process::abort();
        }
        st.dbasis.resize(3 * need_nd, 0.0);
        // `basis` doubles as the (nd)-shaped shape-function buffer; it is
        // already `maxn`-sized, which covers `nd <= maxn` on valid meshes.
        // A p-refined element with more DOFs than nodes grows it here.
        if need_nd > st.basis.len() {
            if st.basis.try_reserve(need_nd - st.basis.len()).is_err() {
                unsafe { moa_fatal_alloc() };
                std::process::abort();
            }
            st.basis.resize(need_nd, 0.0);
        }
        st.maxnd = need_nd;
    }
}

/// `advdiffsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model` and `dt` are unused
/// by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `transient` must be a valid read of the caller's `LOGICAL`
/// `TransientSimulation` (gfortran 4-byte convention, nonzero = true).
#[export_name = "advdiffsolver_"]
pub unsafe extern "C" fn adv_diff_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let transient_sim = !transient.is_null() && unsafe { *transient } != 0;
    let mut st = STATE.lock().unwrap();
    unsafe { moa_begin(solver) };

    unsafe { moa_default_start() };

    let maxiter = unsafe { moa_max_iter() };

    // Totals are host-associated in the original: reset per nonlinear
    // iteration, reported after the loop (final-iteration values).
    let (mut tot_area, mut tot_len, mut tot_src): (c_double, c_double, c_double) =
        (0.0, 0.0, 0.0);

    // Nonlinear iteration loop.
    for _ in 1..=maxiter {
        unsafe { moa_default_initialize() };

        tot_area = 0.0;
        tot_len = 0.0;
        tot_src = 0.0;

        // System assembly, with the CutFEM `GOTO 1` re-assembly loop:
        // re-entry re-reads `GetNOFActive()` without re-initializing.
        loop {
            let active = unsafe { moa_nof_active() };
            for t in 1..=active {
                let (mut n, mut nd, mut nb) = (0, 0, 0);
                unsafe { moa_begin_active(t, &mut n, &mut nd, &mut nb) };
                // `LocalMatrix(Element, n, nd+nb)`: full size covers
                // retained + bubble DOFs.
                let full = (nd + nb) as usize;
                ensure(&mut st, n as usize, full);
                local_matrix(&mut st, n, nd, nb, transient_sim, &mut tot_area, &mut tot_src);
            }

            unsafe { moa_finish_bulk() };

            let nbnd = unsafe { moa_num_bnd() };
            for t in 1..=nbnd {
                let (mut n, mut nd, mut active_el) = (0, 0, 0);
                unsafe { moa_begin_boundary(t, &mut n, &mut nd, &mut active_el) };
                if active_el == 0 {
                    continue;
                }
                ensure(&mut st, n as usize, nd as usize);
                // `BC => GetBC(); IF(.NOT.ASSOCIATED(BC)) RETURN`.
                if unsafe { moa_bc_begin() } == 0 {
                    continue;
                }
                local_matrix_bc(&mut st, n, nd, &mut tot_len);
            }

            if unsafe { moa_cut_fem() } == 0 {
                break;
            }
        }

        unsafe {
            moa_finish_boundary();
            moa_finish_assembly();
            moa_dirichlet();
        }

        // And finally, solve (Norm is local and unused, as in the original).
        let _norm: c_double = unsafe { moa_solve() };
        if unsafe { moa_converged() } != 0 {
            break;
        }
    }

    unsafe { moa_report(tot_area, tot_len, tot_src) };

    unsafe { moa_default_finish() };
}

/// The `LocalMatrix` internal procedure: bulk element matrix and
/// right-hand side by adaptive numerical integration. `n` is the element
/// node count; `nd` the retained DOF count; `nb` the bubble DOF count
/// (`nd + nb` is the procedure's `nd` dummy).
#[allow(clippy::too_many_arguments)]
fn local_matrix(
    st: &mut Scratch,
    n: c_int,
    nd: c_int,
    nb: c_int,
    transient_sim: bool,
    tot_area: &mut c_double,
    tot_src: &mut c_double,
) {
    let full = nd + nb;
    let un = n as usize;
    let uf = full as usize;
    // `MASS = 0; STIFF = 0; FORCE = 0; LOAD = 0` (whole-array fills).
    st.mass[..uf * uf].fill(0.0);
    st.stiff[..uf * uf].fill(0.0);
    st.force[..uf].fill(0.0);
    st.load[..un].fill(0.0);

    unsafe { moa_lm_nodes() };

    let mut found = 0;
    if unsafe { moa_bodyforce() } != 0 {
        unsafe { moa_load_source(n, st.load.as_mut_ptr(), &mut found) };
    }
    let _ = found;

    unsafe { moa_mat_begin() };
    for (which, buf) in [
        (1, &mut st.diffc),
        (2, &mut st.reactc),
        (3, &mut st.convc),
        (4, &mut st.timec),
    ] {
        let ptr: *mut c_double = buf.as_mut_ptr();
        unsafe { moa_mat_coeff(which, n, ptr, &mut found) };
    }
    let _ = found;

    // `Velo = 0`; rows `1..=dim` filled from the material.
    st.velo[..3 * un].fill(0.0);
    let dim = unsafe { moa_coord_dim() };
    for i in 1..=dim {
        unsafe { moa_conv_vel(i, n, st.velo.as_mut_ptr()) };
    }

    let nip = unsafe { moa_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { moa_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            moa_lm_eleminfo(
                u, v, w, n, full,
                st.basis.as_mut_ptr(), st.dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;
        let weight = s * detj;
        let mut loadatip = 0.0;
        unsafe {
            moa_lm_ip(
                n, full, dim,
                st.basis.as_mut_ptr(), st.dbasis.as_mut_ptr(),
                st.load.as_mut_ptr(), st.timec.as_mut_ptr(), st.velo.as_mut_ptr(),
                st.diffc.as_mut_ptr(), st.convc.as_mut_ptr(), st.reactc.as_mut_ptr(),
                weight,
                st.stiff.as_mut_ptr(), st.mass.as_mut_ptr(), st.force.as_mut_ptr(),
                &mut loadatip, tot_area, tot_src,
            )
        };
    }

    if transient_sim {
        unsafe { moa_lm_transient(full, st.mass.as_mut_ptr(), st.stiff.as_mut_ptr(), st.force.as_mut_ptr()) };
    }
    unsafe {
        moa_condensate(nd, nb, full, st.stiff.as_mut_ptr(), st.force.as_mut_ptr());
        moa_default_update(full, st.stiff.as_mut_ptr(), st.force.as_mut_ptr());
    }
}

/// The `LocalMatrixBC` internal procedure: Neumann/Robin boundary
/// contribution. `n` is the boundary element node count, `nd` its DOF
/// count (no bubble condensation on the boundary path).
fn local_matrix_bc(st: &mut Scratch, n: c_int, nd: c_int, tot_len: &mut c_double) {
    let und = nd as usize;
    unsafe { moa_bc_nodes() };
    // `STIFF = 0; FORCE = 0` (`LOAD = 0` is unobservable: never read).
    st.stiff[..und * und].fill(0.0);
    st.force[..und].fill(0.0);

    unsafe {
        moa_bc_coeffs(
            n,
            st.flux.as_mut_ptr(),
            st.coeff.as_mut_ptr(),
            st.ext.as_mut_ptr(),
        )
    };
    let nip = unsafe { moa_bc_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { moa_bc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            moa_bc_eleminfo(
                u, v, w, n, nd,
                st.basis.as_mut_ptr(), st.dbasis.as_mut_ptr(),
                &mut detj, &mut stat,
            )
        };
        let _ = stat;
        let weight = s * detj;
        unsafe {
            moa_bc_ip(
                n, nd,
                st.basis.as_mut_ptr(),
                st.flux.as_mut_ptr(), st.coeff.as_mut_ptr(), st.ext.as_mut_ptr(),
                weight,
                st.stiff.as_mut_ptr(), st.force.as_mut_ptr(),
                tot_len,
            )
        };
    }
    unsafe {
        moa_default_update(nd, st.stiff.as_mut_ptr(), st.force.as_mut_ptr());
    }
}
