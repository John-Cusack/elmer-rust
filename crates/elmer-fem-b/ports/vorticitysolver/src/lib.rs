//! Port of `fem/src/modules/VorticitySolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/VorticitySolver.F90` — two F77-global
//! entries `vorticitysolver_` and `vorticitysolver_init_` with one
//! `CONTAINS` internal (`BulkAssembly`). The emitted dynamic symbols are
//! exactly `vorticitysolver_` + `vorticitysolver_init_` (proven by `nm -D`
//! on the reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/vorticitysolver_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive reductions (`MATMUL`,
//! `SUM`), formatted I/O (`PRINT`/`WRITE`) and whole-array updates stay
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors VorticitySolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: BulkAssembly's `SAVE Nodes` / `SAVE
//!   Coeff` live as wrapper-side SAVE state with identical lifetimes
//!   (allocated at call start, deallocated at call end, as in the
//!   original); the entry's `ForceVectors`/`SaveRHS` context likewise
//!   lives wrapper-side across the entry's calls. The entry's `SAVE
//!   Visited` is dead state in the original (set once, never read) and
//!   is not replicated.
//! - `STIFF`/`FORCE`/`Vx`/`Vy`/`Vz`/`Basis`/`dBasisdx` are per-call
//!   scratch sized by `MAX(MaxElementDOFs, MaxElementNodes)` at call
//!   start, exactly as in the original; a grow-only guard covers
//!   remeshed re-entry where the original would overflow silently.
//! - `C` (the unused identity matrix) is not replicated: it is written
//!   and never read.
//! - Allocation uses `try_reserve` semantics: entry-context allocation
//!   failure would call the same `Fatal(...)` path via the wrapper.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body,
//!   exactly as in the original (received and ignored).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn voa_begin(h_solver: *const c_void);
    fn voa_info_header();
    fn voa_dim() -> c_int;
    fn voa_assoc_matrix() -> c_int;
    fn voa_active_var_count() -> c_int;
    fn voa_params_begin();
    fn voa_setup_2d(dofs: *mut c_int);
    fn voa_print_dofs(dofs: c_int);
    fn voa_fatal_dofs();
    fn voa_setup_3d(gotit: *mut c_int, dofs: *mut c_int, state: *mut c_int);
    fn voa_fatal_dofs3();
    fn voa_fatal_missing();
    fn voa_setup_csym();
    fn voa_get_varname(which: c_int, gotit: *mut c_int);
    fn voa_default_varname();
    fn voa_get_coeff(gotit: *mut c_int);
    fn voa_realtime() -> c_double;
    fn voa_get_constbulk(gotit: *mut c_int, val: *mut c_int);
    fn voa_setup_inuse(constbulk: c_int);
    fn voa_inuse() -> c_int;
    fn voa_default_initialize();
    fn voa_rhs_size() -> c_int;
    fn voa_fv_begin(n: c_int, dofs: c_int);
    fn voa_fv_end();
    fn voa_rhs_to_save();
    fn voa_rhs_to_fv(i: c_int);
    fn voa_max_dofs() -> c_int;
    fn voa_max_nodes() -> c_int;
    fn voa_nof_active() -> c_int;
    fn voa_bulk_alloc(n: c_int);
    fn voa_bulk_free();
    fn voa_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int);
    fn voa_local_solution(comp: c_int, n: c_int, x: *mut c_double);
    fn voa_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn voa_ip_update(
        t: c_int,
        n: c_int,
        nd: c_int,
        dim: c_int,
        dofs: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
        vx: *mut c_double,
        vy: *mut c_double,
        vz: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
    );
    fn voa_update_equations(n: c_int, nd: c_int, df: c_int, stiff: *mut c_double, force: *mut c_double);
    fn voa_update_force(i: c_int, nd: c_int, df: c_int, force: *mut c_double);
    fn voa_finish_bulk();
    fn voa_finish_bulk_const();
    fn voa_default_finish();
    fn voa_report_assembly(at0: c_double, at1: c_double);
    fn voa_default_solve() -> c_double;
    fn voa_copy_component(i: c_int, dofs: c_int);
    fn voa_set_norm(tot: c_double);
    fn voa_report_solution(at1: c_double, at2: c_double);
    fn voa_report_norm(tot: c_double);
    fn voi_begin(h_solver: *const c_void);
    fn voi_get_varname(gotit: *mut c_int);
    fn voi_check_variable() -> c_int;
    fn voi_add_var_2d();
    fn voi_add_var_3d();
    fn voi_add_abs();
    fn voi_fatal_baddim();
    fn voi_get_calcabs() -> c_int;
    fn voi_add_tail();
}

/// Serial-entry guard for the wrapper-side SAVE state. Elmer calls solver
/// entries serially; this is only a `static`-sharing guard, never
/// contended across threads.
static GUARD: Mutex<()> = Mutex::new(());

/// Grow per-call scratch when an element needs more than the first-mesh
/// sizing. The original never grows (it sizes once per call from the
/// current mesh maxima); growth only triggers where the original would
/// overflow silently, and the leading sections used by the wrappers are
/// unaffected, so valid inputs never observe it.
fn ensure_len(buf: &mut Vec<c_double>, need: usize) {
    if buf.len() < need {
        buf.resize(need, 0.0);
    }
}

/// `vorticitysolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "vorticitysolver_"]
pub unsafe extern "C" fn vorticity_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe { voa_begin(solver) };
    unsafe { voa_info_header() };

    let dim = unsafe { voa_dim() };

    // IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { voa_assoc_matrix() } == 0 {
        return;
    }
    // IF (COUNT(Solver % Variable % Perm > 0) <= 0) RETURN
    if unsafe { voa_active_var_count() } <= 0 {
        return;
    }

    unsafe { voa_params_begin() };

    // IF (dim == 2) ... ELSE IF (dim == 3) ... (Dofs stays unset
    // otherwise, as in the original; the mesh dimension is 2 or 3 in
    // every exercisable configuration).
    let mut dofs: c_int = 0;
    if dim == 2 {
        unsafe { voa_setup_2d(&mut dofs) };
        if dofs != 1 {
            unsafe {
                voa_print_dofs(dofs);
                voa_fatal_dofs();
            }
        }
    } else if dim == 3 {
        let (mut gotit, mut state) = (0, 0);
        unsafe { voa_setup_3d(&mut gotit, &mut dofs, &mut state) };
        let _ = gotit;
        if state == 1 {
            unsafe { voa_fatal_missing() };
        } else if state == 2 {
            unsafe {
                voa_print_dofs(dofs);
                voa_fatal_dofs3();
            }
        }
    }

    unsafe { voa_setup_csym() };

    // VarName = GetString('Vorticity Variable') [.or. 'Target Variable'
    // .or. 'Velocity']; CondName = ListGetString('Vorticity Coefficient').
    let mut gotit = 0;
    unsafe { voa_get_varname(1, &mut gotit) };
    if gotit == 0 {
        unsafe { voa_get_varname(2, &mut gotit) };
        if gotit == 0 {
            unsafe { voa_default_varname() };
        }
    }
    let mut gotcoeff = 0;
    unsafe { voa_get_coeff(&mut gotcoeff) };
    let _ = gotcoeff;

    let at0 = unsafe { voa_realtime() };

    let (mut got_cb, mut cb) = (0, 0);
    unsafe { voa_get_constbulk(&mut got_cb, &mut cb) };
    let _ = got_cb;
    let constbulk = cb != 0;
    unsafe { voa_setup_inuse(cb) };

    unsafe { voa_default_initialize() };

    // IF (Dofs > 1) ALLOCATE(ForceVectors(...)) ...
    if dofs > 1 {
        let n = unsafe { voa_rhs_size() };
        unsafe { voa_fv_begin(n, dofs) };
    }

    bulk_assembly(dim, dofs);

    if constbulk {
        unsafe { voa_finish_bulk_const() };
    } else {
        unsafe { voa_finish_bulk() };
    }

    // No flux BCs.
    unsafe { voa_default_finish() };

    let at1 = unsafe { voa_realtime() };
    unsafe { voa_report_assembly(at0, at1) };

    let totnorm: c_double = if dofs > 1 {
        let mut tot = 0.0;
        for i in 1..=dofs {
            if i == 1 {
                unsafe { voa_rhs_to_save() };
            } else {
                unsafe { voa_rhs_to_fv(i) };
            }
            let unorm = unsafe { voa_default_solve() };
            tot += unorm * unorm;
            unsafe { voa_copy_component(i, dofs) };
        }
        tot = tot.sqrt();
        unsafe { voa_fv_end() };
        unsafe { voa_set_norm(tot) };
        tot
    } else {
        unsafe { voa_default_solve() }
    };
    let at2 = unsafe { voa_realtime() };
    unsafe { voa_report_solution(at1, at2) };
    unsafe { voa_report_norm(totnorm) };
}

/// The `BulkAssembly` internal procedure: mass-matrix assembly plus the
/// vorticity right-hand side by numerical integration. `dim`/`dofs` are
/// the entry's values; scratch is sized by
/// `MAX(MaxElementDOFs, MaxElementNodes)` exactly as in the original.
fn bulk_assembly(dim: c_int, dofs: c_int) {
    let nmax = unsafe { voa_max_dofs().max(voa_max_nodes()) } as usize;
    // FORCE has (dim, n) entries; the declared row count follows the
    // original. The buffer is widened only for dimensions the original
    // leaves undefined (never exercised), so valid inputs never observe
    // it.
    let rows = dim.max(3) as usize;
    let mut stiff = vec![0.0 as c_double; nmax * nmax];
    let mut force = vec![0.0 as c_double; rows * nmax];
    let mut vx = vec![0.0 as c_double; nmax];
    let mut vy = vec![0.0 as c_double; nmax];
    let mut vz = vec![0.0 as c_double; nmax];
    let mut basis = vec![0.0 as c_double; nmax];
    let mut dbasis = vec![0.0 as c_double; nmax * 3];

    unsafe { voa_bulk_alloc(nmax as c_int) };

    let active = unsafe { voa_nof_active() };
    for t in 1..=active {
        let (mut n, mut nd) = (0, 0);
        unsafe { voa_begin_active(t, &mut n, &mut nd) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        let need = n as usize;
        ensure_len(&mut vx, need);
        ensure_len(&mut vy, need);
        ensure_len(&mut vz, need);
        ensure_len(&mut basis, need);
        ensure_len(&mut dbasis, need * 3);

        stiff.fill(0.0);
        force.fill(0.0);

        unsafe {
            voa_local_solution(1, n, vx.as_mut_ptr());
            voa_local_solution(2, n, vy.as_mut_ptr());
            voa_local_solution(3, n, vz.as_mut_ptr());
        }

        let nip = unsafe { voa_ip_begin() };
        for q in 1..=nip {
            unsafe {
                voa_ip_update(
                    q,
                    n,
                    nd,
                    dim,
                    dofs,
                    stiff.as_mut_ptr(),
                    force.as_mut_ptr(),
                    vx.as_mut_ptr(),
                    vy.as_mut_ptr(),
                    vz.as_mut_ptr(),
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                );
            }
        }

        if unsafe { voa_inuse() } == 0 {
            if dofs > 1 {
                unsafe { voa_rhs_to_save() };
            }
            unsafe {
                voa_update_equations(
                    nmax as c_int,
                    nd,
                    dim,
                    stiff.as_mut_ptr(),
                    force.as_mut_ptr(),
                );
            }
        } else {
            unsafe { voa_update_force(1, nd, dim, force.as_mut_ptr()) };
        }

        // Assembly of the 2nd and 3rd r.h.s. in the 3D case.
        for i in 2..=dofs {
            unsafe { voa_rhs_to_fv(i) };
            unsafe { voa_update_force(i, nd, dim, force.as_mut_ptr()) };
        }
    }

    unsafe { voa_bulk_free() };
}

/// `vorticitysolver_init_` — the solver init entry. Adds the `Variable`,
/// linear-system defaults and optional `Exported Variable` entries, in
/// the original's order.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "vorticitysolver_init_"]
pub unsafe extern "C" fn vorticity_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = GUARD.lock().unwrap();
    unsafe { voi_begin(solver) };

    let dim = unsafe { voa_dim() };

    let mut gotit = 0;
    unsafe { voi_get_varname(&mut gotit) };
    let _ = gotit;

    if unsafe { voi_check_variable() } == 0 {
        if dim == 2 {
            unsafe { voi_add_var_2d() };
        } else if dim == 3 {
            unsafe { voi_add_var_3d() };
        } else {
            unsafe { voi_fatal_baddim() };
        }
    }

    if unsafe { voi_get_calcabs() } != 0 {
        unsafe { voi_add_abs() };
    }

    unsafe { voi_add_tail() };
}
