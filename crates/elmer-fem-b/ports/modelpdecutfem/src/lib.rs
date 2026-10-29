//! Port of `fem/src/modules/ModelPDECutFEM.F90` (Elmer @ 9f6af2f85, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ModelPDECutFEM.F90` — two F77-global
//! entries `advdiffsolver_` and `advdiffsolver_init_` with two `CONTAINS`
//! internals (`LocalMatrix`, `LocalMatrixBC`). The emitted dynamic
//! symbols are exactly `advdiffsolver_` + `advdiffsolver_init_` (proven by
//! `nm -D` on the reference `.so`); the internals emit only local
//! symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/modelpdecutfem_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL`, `SUM`), the `PRINT` interface-area report, and
//! whole-array updates stay in Fortran wrappers, so bit-identity holds
//! by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ModelPDECutFEM.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the host-associated `MatVol(0:10)` is
//!   zeroed on every solver entry and accumulated/printed wrapper-side;
//!   the internal procedures' `SAVE Nodes` live as wrapper-side SAVE
//!   state with identical lifetimes (consumed before the next fill, as in
//!   the original).
//! - The cut-element loops (`GOTO 10` / `GOTO 20` re-splitting) are
//!   transcribed as Rust `loop { ... }` with the same call sequence,
//!   including falling through to the plain `LocalMatrix` when the bulk
//!   splitter reports no cut.
//! - `Model % CurrentElement => pElement` ahead of cut BC assembly is
//!   kept as a wrapper touch call, exactly where the original puts it.
//! - Allocation failure calls the same `Fatal('AdvDiffSolver', ...)` via
//!   `try_reserve` (no `STAT=istat` in the original's entry; the one
//!   allocation site is the first-mesh sizing), instead of aborting.
//! - `dt` is unused by the body, exactly as in the original (received
//!   and ignored).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn mcf_begin(h_model: *const c_void, h_solver: *const c_void);
    fn mcf_cutfem_flag() -> c_int;
    fn mcf_init_add_flags();
    fn mcf_max_elem_dofs() -> c_int;
    fn mcf_maxiter(found: *mut c_int) -> c_int;
    fn mcf_default_start();
    fn mcf_default_finish();
    fn mcf_get_nof_active() -> c_int;
    fn mcf_default_finish_bulk();
    fn mcf_get_nof_boundary() -> c_int;
    fn mcf_default_finish_boundary();
    fn mcf_default_finish_assembly();
    fn mcf_default_dirichlet();
    fn mcf_default_solve() -> c_double;
    fn mcf_default_converged() -> c_int;
    fn mcf_set_norm(norm: c_double);
    fn mcf_fatal_alloc();
    fn mcf_matvol_zero();
    fn mcf_matvol_report();
    fn mcf_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn mcf_cut_check(iscut: *mut c_int, isactive: *mut c_int);
    fn mcf_cut_bulk(iscut: *mut c_int, ismore: *mut c_int);
    fn mcf_pelem_dims(n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn mcf_pelem_perm_ok() -> c_int;
    fn mcf_begin_cutbc_active(t: c_int);
    fn mcf_cut_bc(assoc: *mut c_int, iscut: *mut c_int, ismore: *mut c_int);
    fn mcf_pelem_typenodes() -> c_int;
    fn mcf_begin_boundary(t: c_int);
    fn mcf_bnd_perm_has_zero() -> c_int;
    fn mcf_bnd_active(n: *mut c_int, nd: *mut c_int) -> c_int;
    fn mcf_lm_select(cut: c_int);
    fn mcf_lm_nodes();
    fn mcf_coord_dim() -> c_int;
    fn mcf_lm_bodyforce() -> c_int;
    fn mcf_lm_load(n: c_int, load: *mut c_double, found: *mut c_int);
    fn mcf_lm_mat(
        n: c_int,
        diff: *mut c_double,
        conv: *mut c_double,
        react: *mut c_double,
        timec: *mut c_double,
    );
    fn mcf_lm_velo(i: c_int, n: c_int, row: *mut c_double);
    fn mcf_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn mcf_lm_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn mcf_lm_eleminfo(
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
    fn mcf_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    fn mcf_lm_amatvec(n: c_int, velo: *const c_double, basis: *const c_double, a: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn mcf_lm_assemble(
        nd: c_int,
        dim: c_int,
        ld: c_int,
        weight: c_double,
        d: c_double,
        c: c_double,
        r: c_double,
        rho: c_double,
        loadatip: c_double,
        avec: *const c_double,
        basis: *const c_double,
        dbasis: *const c_double,
        stiff: *mut c_double,
        mass: *mut c_double,
        force: *mut c_double,
        docut: c_int,
    );
    fn mcf_lm_first_order_time(ld: c_int, mass: *mut c_double, stiff: *mut c_double, force: *mut c_double);
    fn mcf_lm_update_mass(n: c_int, nd: c_int, ld: c_int, mass: *mut c_double);
    fn mcf_condensate(ndmnb: c_int, nb: c_int, ld: c_int, k: *mut c_double, f: *mut c_double);
    fn mcf_lm_update_global(nd: c_int, ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn mcf_lm_default_update(nd: c_int, ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn mcf_lbc_select(cut: c_int);
    fn mcf_lbc_begin() -> c_int;
    fn mcf_lbc_nodes();
    fn mcf_lbc_vals(n: c_int, flux: *mut c_double, coeff: *mut c_double, ext: *mut c_double);
    fn mcf_lbc_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn mcf_lbc_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn mcf_lbc_eleminfo(
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
    fn mcf_lbc_assemble(
        nd: c_int,
        ld: c_int,
        weight: c_double,
        f: c_double,
        c: c_double,
        ext: c_double,
        basis: *const c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn mcf_lbc_update_global(
        nd: c_int,
        ld: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        lastweight: c_double,
    );
    fn mcf_lbc_default_update(nd: c_int, ld: c_int, stiff: *mut c_double, force: *mut c_double);
}
/// Grow the SAVE arrays when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its automatics are per-element); growth only triggers where the
/// original would overflow, and the leading sections used by the
/// wrappers are unaffected, so valid inputs never observe it.
fn ensure_state(state: &mut CutFemState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.mass.resize(need * need, 0.0);
    state.load.resize(need, 0.0);
    state.force.resize(need, 0.0);
    state.coeff.resize(need * 4, 0.0);
    state.velo.resize(need * 3, 0.0);
    state.n = need;
}
/// First-mesh SAVE state: `STIFF(N,N)`, `MASS(N,N)`, `LOAD(N)`,
/// `FORCE(N)`, four coefficient rows and the `3xN` velocity block.
/// `N` is the first call's `MaxElementDOFs`.
struct CutFemState {
    stiff: Vec<c_double>,
    mass: Vec<c_double>,
    load: Vec<c_double>,
    force: Vec<c_double>,
    coeff: Vec<c_double>,
    velo: Vec<c_double>,
    n: usize,
    allocations_done: bool,
}

impl CutFemState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            mass: Vec::new(),
            load: Vec::new(),
            force: Vec::new(),
            coeff: Vec::new(),
            velo: Vec::new(),
            n: 0,
            allocations_done: false,
        }
    }
}

static STATE: Mutex<CutFemState> = Mutex::new(CutFemState::new());

/// `advdiffsolver_init_` — the solver init entry. Enables the global
/// mass matrix path when the solver requests CutFEM, exactly as in the
/// original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "advdiffsolver_init_"]
pub unsafe extern "C" fn advdiff_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { mcf_begin(model, solver) };
    // IF ( ListGetLogical( Solver % Values,'CutFEM',Found ) ) THEN
    if unsafe { mcf_cutfem_flag() } != 0 {
        unsafe { mcf_init_add_flags() };
    }
}

/// `advdiffsolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` is unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `transient` must be a valid read if non-null (it selects the
/// transient assembly tail, as in the original).
#[export_name = "advdiffsolver_"]
pub unsafe extern "C" fn advdiff_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { mcf_begin(model, solver) };
    let transient_sim = !transient.is_null() && unsafe { *transient } != 0;

    unsafe { mcf_default_start() };

    // In this version we do on-the-fly splitting of elements.
    let docut = unsafe { mcf_cutfem_flag() } != 0;

    // Host-associated MatVol reset on every entry.
    unsafe { mcf_matvol_zero() };

    let mut maxiter = {
        let mut found = 0;
        let m = unsafe { mcf_maxiter(&mut found) };
        if found == 0 { 1 } else { m }
    };
    if maxiter < 1 {
        maxiter = 1;
    }

    // First-mesh permanent storage (bug-compatible: sized once).
    if !state.allocations_done {
        let n = unsafe { mcf_max_elem_dofs() } as usize;
        state.n = n;
        if state.stiff.try_reserve(n * n).is_err()
            || state.mass.try_reserve(n * n).is_err()
            || state.load.try_reserve(n).is_err()
            || state.force.try_reserve(n).is_err()
            || state.coeff.try_reserve(n * 4).is_err()
            || state.velo.try_reserve(n * 3).is_err()
        {
            unsafe { mcf_fatal_alloc() };
            std::process::abort();
        }
        state.stiff.resize(n * n, 0.0);
        state.mass.resize(n * n, 0.0);
        state.load.resize(n, 0.0);
        state.force.resize(n, 0.0);
        state.coeff.resize(n * 4, 0.0);
        state.velo.resize(n * 3, 0.0);
        state.allocations_done = true;
    }

    // Nonlinear iteration loop.
    for _ in 1..=maxiter {
        // Basic system assembly, skip the cut and empty elements.
        let active = unsafe { mcf_get_nof_active() };
        for t in 1..=active {
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { mcf_begin_active(t, &mut n, &mut nd, &mut nb) };
            ensure_state(&mut state, (n as usize).max((nd + nb) as usize));

            // We need to add this special section where elements are cut
            // on-the-fly.
            if docut {
                let (mut iscut, mut isactive) = (0, 0);
                unsafe { mcf_cut_check(&mut iscut, &mut isactive) };
                let _ = iscut;
                if isactive == 0 {
                    continue;
                }
                // If the element is cut then split it on-the-fly until
                // there are no more pieces (label 10 loop).
                let mut cut_done = false;
                loop {
                    let (mut iscut, mut ismore) = (0, 0);
                    unsafe { mcf_cut_bulk(&mut iscut, &mut ismore) };
                    if iscut == 0 {
                        break;
                    }
                    let (mut pn, mut pnd, mut pnb) = (0, 0, 0);
                    unsafe { mcf_pelem_dims(&mut pn, &mut pnd, &mut pnb) };
                    ensure_state(&mut state, (pn as usize).max((pnd + pnb) as usize));
                    if unsafe { mcf_pelem_perm_ok() } != 0 {
                        local_matrix(&mut state, pn, pnd + pnb, pnb, true, transient_sim, docut);
                    }
                    if ismore == 0 {
                        cut_done = true;
                        break;
                    }
                }
                if cut_done {
                    continue;
                }
            }

            local_matrix(&mut state, n, nd + nb, nb, false, transient_sim, docut);
        }

        unsafe { mcf_default_finish_bulk() };

        // Loop over bulk elements where BCs are created on-the-fly if the
        // element is cut (label 20 loop).
        if docut {
            for t in 1..=active {
                unsafe { mcf_begin_cutbc_active(t) };
                loop {
                    let (mut assoc, mut iscut, mut ismore) = (0, 0, 0);
                    unsafe { mcf_cut_bc(&mut assoc, &mut iscut, &mut ismore) };
                    let _ = iscut;
                    if assoc == 0 {
                        break;
                    }
                    let n = unsafe { mcf_pelem_typenodes() };
                    ensure_state(&mut state, n as usize);
                    if unsafe { mcf_pelem_perm_ok() } != 0 {
                        local_matrix_bc(&mut state, n, n, true);
                    }
                    if ismore == 0 {
                        break;
                    }
                }
            }
        }

        // Now the real boundary assembly.
        let nboundary = unsafe { mcf_get_nof_boundary() };
        for t in 1..=nboundary {
            unsafe { mcf_begin_boundary(t) };
            if docut && unsafe { mcf_bnd_perm_has_zero() } != 0 {
                continue;
            }
            let (mut n, mut nd) = (0, 0);
            if unsafe { mcf_bnd_active(&mut n, &mut nd) } != 0 {
                ensure_state(&mut state, (n as usize).max(nd as usize));
                local_matrix_bc(&mut state, n, nd, false);
            }
        }

        if docut {
            unsafe { mcf_matvol_report() };
        }

        unsafe {
            mcf_default_finish_boundary();
            mcf_default_finish_assembly();
            mcf_default_dirichlet();
        }

        // And finally, solve.
        let norm: c_double = unsafe { mcf_default_solve() };
        if unsafe { mcf_default_converged() } != 0 {
            unsafe { mcf_set_norm(norm) };
            break;
        }
        unsafe { mcf_set_norm(norm) };
    }

    unsafe { mcf_default_finish() };
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count;
/// `ndnb` is `nd+nb` (the internal procedure's `nd` dummy).
#[allow(clippy::too_many_arguments)]
fn local_matrix(
    state: &mut CutFemState,
    n: c_int,
    ndnb: c_int,
    nb: c_int,
    cut: bool,
    transient_sim: bool,
    docut: bool,
) {
    unsafe { mcf_lm_select(c_int::from(cut)) };
    unsafe { mcf_lm_nodes() };
    let dim = unsafe { mcf_coord_dim() };
    let nn = n as usize;
    let ndn = ndnb as usize;
    let nmax = state.n as c_int;

    state.stiff.fill(0.0);
    state.mass.fill(0.0);
    state.force.fill(0.0);
    state.load.fill(0.0);

    if unsafe { mcf_lm_bodyforce() } != 0 {
        let mut found = 0;
        unsafe { mcf_lm_load(n, state.load.as_mut_ptr(), &mut found) };
        let _ = found;
    }

    let (diff, rest) = state.coeff.split_at_mut(nn);
    let (conv, rest) = rest.split_at_mut(nn);
    let (react, timec) = rest.split_at_mut(nn);
    unsafe {
        mcf_lm_mat(
            n,
            diff.as_mut_ptr(),
            conv.as_mut_ptr(),
            react.as_mut_ptr(),
            timec.as_mut_ptr(),
        )
    };
    // `Velo = 0` then per-component fill, kept in column-major `Velo(3,n)`
    // layout so the Fortran MATMUL transcription reads it directly.
    state.velo.fill(0.0);
    let mut row = vec![0.0 as c_double; nn];
    for i in 1..=dim {
        unsafe { mcf_lm_velo(i, n, row.as_mut_ptr()) };
        for j in 0..nn {
            state.velo[j * 3 + (i as usize - 1)] = row[j];
        }
    }
    let (velo, _) = state.velo.split_at_mut(nn * 3);
    let (load, _) = state.load.split_at_mut(nn);

    // Basis(n) / dBasisdx(nd,3) automatics, sized by this element.
    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { mcf_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { mcf_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            mcf_lm_eleminfo(
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

        let load_at_ip: c_double = unsafe { mcf_dot(n, basis.as_ptr(), load.as_ptr()) };
        let rho: c_double = unsafe { mcf_dot(n, basis.as_ptr(), timec.as_ptr()) };
        let mut avec = [0.0 as c_double; 3];
        unsafe { mcf_lm_amatvec(n, velo.as_ptr(), basis.as_ptr(), avec.as_mut_ptr()) };
        let d: c_double = unsafe { mcf_dot(n, basis.as_ptr(), diff.as_ptr()) };
        let c: c_double = unsafe { mcf_dot(n, basis.as_ptr(), conv.as_ptr()) };
        let r: c_double = unsafe { mcf_dot(n, basis.as_ptr(), react.as_ptr()) };

        unsafe {
            mcf_lm_assemble(
                ndnb,
                dim,
                nmax,
                s * detj,
                d,
                c,
                r,
                rho,
                load_at_ip,
                avec.as_ptr(),
                basis.as_ptr(),
                dbasis.as_ptr(),
                state.stiff.as_mut_ptr(),
                state.mass.as_mut_ptr(),
                state.force.as_mut_ptr(),
                c_int::from(docut),
            )
        };
    }

    if transient_sim {
        if cut {
            unsafe { mcf_lm_update_mass(n, ndnb, nmax, state.mass.as_mut_ptr()) };
        } else {
            unsafe {
                mcf_lm_first_order_time(
                    nmax,
                    state.mass.as_mut_ptr(),
                    state.stiff.as_mut_ptr(),
                    state.force.as_mut_ptr(),
                )
            };
        }
    }

    unsafe {
        mcf_condensate(
            ndnb - nb,
            nb,
            nmax,
            state.stiff.as_mut_ptr(),
            state.force.as_mut_ptr(),
        )
    };

    if cut {
        unsafe { mcf_lm_update_global(ndnb, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr()) };
    } else {
        unsafe { mcf_lm_default_update(ndnb, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr()) };
    }
}

/// The `LocalMatrixBC` internal procedure: Neumann/Robin boundary
/// assembly by boundary integration.
fn local_matrix_bc(state: &mut CutFemState, n: c_int, nd: c_int, cut: bool) {
    unsafe { mcf_lbc_select(c_int::from(cut)) };
    // IF (.NOT.ASSOCIATED(BC)) RETURN
    if unsafe { mcf_lbc_begin() } == 0 {
        return;
    }
    unsafe { mcf_lbc_nodes() };
    let nn = n as usize;
    let ndn = nd as usize;
    let nmax = state.n as c_int;

    state.stiff.fill(0.0);
    state.force.fill(0.0);

    let mut flux = vec![0.0 as c_double; nn];
    let mut coeff = vec![0.0 as c_double; nn];
    let mut ext = vec![0.0 as c_double; nn];
    unsafe { mcf_lbc_vals(n, flux.as_mut_ptr(), coeff.as_mut_ptr(), ext.as_mut_ptr()) };

    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { mcf_lbc_ip_begin() };
    let mut lastweight = 0.0;
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { mcf_lbc_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            mcf_lbc_eleminfo(
                u,
                v,
                w,
                n,
                nd,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            )
        };
        let _ = stat;

        let f: c_double = unsafe { mcf_dot(n, basis.as_ptr(), flux.as_ptr()) };
        let c: c_double = unsafe { mcf_dot(n, basis.as_ptr(), coeff.as_ptr()) };
        let extv: c_double = unsafe { mcf_dot(n, basis.as_ptr(), ext.as_ptr()) };

        let weight = s * detj;
        lastweight = weight;
        unsafe {
            mcf_lbc_assemble(
                nd,
                nmax,
                weight,
                f,
                c,
                extv,
                basis.as_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }

    if cut {
        unsafe {
            mcf_lbc_update_global(
                nd,
                nmax,
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
                lastweight,
            )
        };
    } else {
        unsafe { mcf_lbc_default_update(nd, nmax, state.stiff.as_mut_ptr(), state.force.as_mut_ptr()) };
    }
}
