//! Port of `fem/src/modules/MeshSolve.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/MeshSolve.F90` — two F77-global entries
//! `meshsolver_` and `meshsolver_init_` with two `CONTAINS` internals.
//! The emitted dynamic symbols are exactly those two (proven by `nm -D`
//! on the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses; `dt` crosses by
//! value into the accessor (the velocity update needs it) and the
//! transient flag as an `int`. Every component access and every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/meshsolve_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), the `(p,q,i,j)` assembly nests, formatted `WRITE(Message)`
//! output, and the timestep velocity `SELECT CASE` stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches.
//!
//!
//! - `SAVE` semantics replicated: `STIFF`/`FORCE`/`LOAD`/`ElasticModulus`/
//!   `PoissonRatio`/`Alpha`/`Beta`/`TPerm` live in the Rust state sized at
//!   first-mesh sizing (`MaxElementDOFs`, re-sized on `MeshChanged`) with
//!   the `AllocationsDone` flag; the wrapper-side SAVE mirrors the
//!   original's process-lifetime scalars (`FirstTime`,
//!   `SkipFirstMeshVelocity`), the velocity alias pair, the per-procedure
//!   `SAVE Nodes` states, and the progress clock, with identical
//!   lifetimes. `MeshVelocity`/`MeshVeloPerm` stay Fortran aliases of the
//!   solver variable storage (never owned or copied by Rust), exactly as
//!   the original's pointers.
//! - Allocation failure calls the same `Fatal('MeshSolve', ...)` via
//!   `try_reserve` (the `STAT=istat` check), never abort directly.
//! - `Model` is received and used only for the `ListCheck*AnyBC` probe,
//!   exactly as in the original (via the accessor-held handle).
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn ms_begin(h_model: *const c_void, h_solver: *const c_void, dt: c_double, transient: c_int);
    fn ms_assoc_matrix() -> c_int;
    fn ms_init_body();
    fn ms_frame();
    fn ms_local_nodes() -> c_int;
    /// Refresh `STDOFs` from the `Displacement` variable when present;
    /// returns the effective value used for sizing and `DisplaceMesh`.
    fn ms_stdofs() -> c_int;
    fn ms_dim() -> c_int;
    fn ms_stress_assoc() -> c_int;
    fn ms_displace_first() -> c_int;
    fn ms_skip_displace(r: *mut c_int);
    fn ms_mesh_changed() -> c_int;
    fn ms_nperm() -> c_int;
    fn ms_max_elem_dofs() -> c_int;
    fn ms_tperm_fill(tperm: *mut c_int, n: c_int);
    fn ms_disp_mu_tperm(tperm: *mut c_int, n: c_int, dir: c_int);
    fn ms_disp_disp_stress(dir: c_int);
    fn ms_disp_disp_stress_norhs();
    fn ms_disp_mu_meshperm(dir: c_int);
    fn ms_stress_vals();
    fn ms_fatal_alloc();
    fn ms_tick_begin();
    fn ms_max_iter() -> c_int;
    fn ms_default_initialize();
    fn ms_nof_active() -> c_int;
    fn ms_tick_maybe(t: c_int, nactive: c_int);
    fn ms_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn ms_mat_props(n: c_int, nn: c_int, emod: *mut c_double, pr: *mut c_double);
    fn ms_lm_nodes();
    #[allow(clippy::too_many_arguments)]
    fn ms_lame(
        n: c_int,
        nn: c_int,
        emod: *mut c_double,
        pr: *mut c_double,
        lame1: *mut c_double,
        lame2: *mut c_double,
        plane: c_int,
    );
    fn ms_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn ms_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn ms_lm_eleminfo(
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
    fn ms_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn ms_lm_assemble(
        ld: c_int,
        dim: c_int,
        nd: c_int,
        s: c_double,
        l1: c_double,
        l2: c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn ms_bubble(
        ld: c_int,
        dim: c_int,
        nd: c_int,
        n: c_int,
        nb: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn ms_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn ms_finish_bulk();
    fn ms_doit_bc() -> c_int;
    fn ms_nbnd() -> c_int;
    fn ms_begin_boundary(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        nb: *mut c_int,
        active: *mut c_int,
    );
    fn ms_bc_assoc() -> c_int;
    fn ms_bc_penalty(i: c_int, n: c_int, nn: c_int, alpha: *mut c_double);
    fn ms_bc_loads(n: c_int, nn: c_int, load: *mut c_double, beta: *mut c_double) -> c_int;
    fn ms_mb_nodes();
    fn ms_mb_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn ms_mb_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn ms_mb_eleminfo(
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
    fn ms_mb_assemble(
        ld: c_int,
        dim: c_int,
        n: c_int,
        nd: c_int,
        u: c_double,
        v: c_double,
        w: c_double,
        s: c_double,
        load: *mut c_double,
        alpha: *mut c_double,
        beta: *mut c_double,
        basis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn ms_finish_assembly();
    fn ms_dirichlet();
    fn ms_default_solve() -> c_double;
    fn ms_converged() -> c_int;
    fn ms_velocity();
}

/// First-mesh SAVE state: `STIFF(LD,LD)`, `FORCE(LD)`, `LOAD(4,N)`,
/// `ElasticModulus(6,6,N)`, `PoissonRatio(N)`, `Alpha(3,N)`, `Beta(N)`,
/// `TPerm(nperm)` plus the `AllocationsDone` flag. `N` is the first
/// call's `MaxElementDOFs`; `LD = STDOFs*N`.
struct MeshSolveState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    emod: Vec<c_double>,
    pr: Vec<c_double>,
    alpha: Vec<c_double>,
    beta: Vec<c_double>,
    tperm: Vec<c_int>,
    ld: usize,
    n: usize,
    allocations_done: bool,
}

impl MeshSolveState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            emod: Vec::new(),
            pr: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            tperm: Vec::new(),
            ld: 0,
            n: 0,
            allocations_done: false,
        }
    }

    /// Reserve (but do not yet fill) the SAVE buffers; on failure calls
    /// the same `Fatal('MeshSolve', ...)` as the original's
    /// `STAT=istat` check.
    fn reserve(&mut self, ld: usize, n: usize) -> bool {
        if self.stiff.try_reserve(ld * ld - self.stiff.len().min(ld * ld)).is_err()
            || self.force.try_reserve(ld - self.force.len().min(ld)).is_err()
            || self.load.try_reserve(4 * n - self.load.len().min(4 * n)).is_err()
            || self.emod.try_reserve(36 * n - self.emod.len().min(36 * n)).is_err()
            || self.pr.try_reserve(n - self.pr.len().min(n)).is_err()
            || self.alpha.try_reserve(3 * n - self.alpha.len().min(3 * n)).is_err()
            || self.beta.try_reserve(n - self.beta.len().min(n)).is_err()
        {
            return false;
        }
        self.stiff.resize(ld * ld, 0.0);
        self.force.resize(ld, 0.0);
        self.load.resize(4 * n, 0.0);
        self.emod.resize(36 * n, 0.0);
        self.pr.resize(n, 0.0);
        self.alpha.resize(3 * n, 0.0);
        self.beta.resize(n, 0.0);
        self.ld = ld;
        self.n = n;
        true
    }
}

static STATE: Mutex<MeshSolveState> = Mutex::new(MeshSolveState::new());

/// `meshsolver_init_` — initialization for the primary solver. Argument
/// order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "meshsolver_init_"]
pub unsafe extern "C" fn meshsolver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let dtv = if dt.is_null() { 0.0 } else { unsafe { *dt } };
    let tr = if transient.is_null() { 0 } else { unsafe { *transient } };
    unsafe { ms_begin(model, solver, dtv, tr) };
    unsafe { ms_init_body() };
}

/// `meshsolver_` — mesh update solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads if non-null.
#[export_name = "meshsolver_"]
pub unsafe extern "C" fn meshsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let dtv = if dt.is_null() { 0.0 } else { unsafe { *dt } };
    let tr = if transient.is_null() { 0 } else { unsafe { *transient } };
    let mut state = STATE.lock().unwrap();
    unsafe { ms_begin(model, solver, dtv, tr) };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { ms_assoc_matrix() } == 0 {
        return;
    }

    unsafe { ms_frame() };

    let local_nodes = unsafe { ms_local_nodes() };
    // IF (LocalNodes <= 0) RETURN
    if local_nodes <= 0 {
        return;
    }

    // Effective STDOFs (refreshed from `Displacement` when present).
    let stdofs = unsafe { ms_stdofs() } as usize;
    let dim = unsafe { ms_dim() } as c_int;
    let stress_assoc = unsafe { ms_stress_assoc() } != 0;
    let displace_first = unsafe { ms_displace_first() } != 0;
    let mut skip_displace_v = 0;
    unsafe { ms_skip_displace(&mut skip_displace_v) };
    let skip_displace = skip_displace_v != 0;
    let mesh_changed = unsafe { ms_mesh_changed() } != 0;

    if stress_assoc {
        unsafe { ms_stress_vals() };
        if !state.allocations_done || mesh_changed {
            // IF (AllocationsDone) DEALLOCATE(TPerm) — dropping the Vec
            // is the deallocate; resizing below is the allocate.
            let nperm = unsafe { ms_nperm() } as usize;
            let have_tp = state.tperm.len();
            if state.tperm.try_reserve(nperm.saturating_sub(have_tp)).is_err() {
                unsafe { ms_fatal_alloc() };
                std::process::abort();
            }
            state.tperm.resize(nperm, 0);
        }
        let nperm = state.tperm.len() as c_int;
        unsafe { ms_tperm_fill(state.tperm.as_mut_ptr(), nperm) };
        if state.allocations_done || displace_first {
            unsafe { ms_disp_mu_tperm(state.tperm.as_mut_ptr(), nperm, -1) };
        }
        unsafe { ms_disp_disp_stress(-1) };
    } else if state.allocations_done || displace_first {
        unsafe { ms_disp_mu_meshperm(-1) };
    }

    // Allocate permanent storage (first call, or on MeshChanged).
    if !state.allocations_done || mesh_changed {
        let n = unsafe { ms_max_elem_dofs() } as usize;
        let ld = stdofs * n;
        if !state.reserve(ld, n) {
            unsafe { ms_fatal_alloc() };
            std::process::abort();
        }
        state.allocations_done = true;
    }
    let (ld, buf_n) = (state.ld as c_int, state.n as c_int);

    unsafe { ms_tick_begin() };

    let max_iter = unsafe { ms_max_iter() };
    for _iter in 1..=max_iter.max(1) {
        unsafe { ms_default_initialize() };

        let nactive = unsafe { ms_nof_active() };
        for t in 1..=nactive {
            unsafe { ms_tick_maybe(t, nactive) };
            let (mut n, mut nd, mut nb) = (0, 0, 0);
            unsafe { ms_begin_active(t, &mut n, &mut nd, &mut nb) };
            unsafe { ms_mat_props(n, buf_n, state.emod.as_mut_ptr(), state.pr.as_mut_ptr()) };
            local_matrix(&mut state, ld, buf_n, dim, n, nd, nb);
            unsafe {
                ms_default_update(ld, state.stiff.as_mut_ptr(), state.force.as_mut_ptr())
            };
        }
        unsafe { ms_finish_bulk() };

        // Neumann & Newton boundary conditions.
        let doit = unsafe { ms_doit_bc() } != 0;
        if doit {
            let nbnd = unsafe { ms_nbnd() };
            for t in 1..=nbnd {
                let (mut n, mut nd, mut nb, mut active) = (0, 0, 0, 0);
                unsafe { ms_begin_boundary(t, &mut n, &mut nd, &mut nb, &mut active) };
                if active == 0 {
                    continue;
                }
                if unsafe { ms_bc_assoc() } == 0 {
                    continue;
                }
                state.load.fill(0.0);
                state.alpha.fill(0.0);
                for i in 1..=dim {
                    unsafe { ms_bc_penalty(i, n, buf_n, state.alpha.as_mut_ptr()) };
                }
                state.beta.fill(0.0);
                let got = unsafe {
                    ms_bc_loads(n, buf_n, state.load.as_mut_ptr(), state.beta.as_mut_ptr())
                };
                if got == 0 {
                    continue;
                }
                mesh_boundary(&mut state, ld, buf_n, dim, n, nd, nb);
                unsafe {
                    ms_default_update(ld, state.stiff.as_mut_ptr(), state.force.as_mut_ptr())
                };
            }
        }

        unsafe { ms_finish_assembly() };
        unsafe { ms_dirichlet() };

        // Solve the system and check for convergence (UNorm is local and
        // unused, as in the original).
        let _unorm: c_double = unsafe { ms_default_solve() };
        if unsafe { ms_converged() } != 0 {
            break;
        }
    }

    if tr != 0 {
        unsafe { ms_velocity() };
    }

    if !skip_displace {
        if stress_assoc {
            let nperm = state.tperm.len() as c_int;
            unsafe { ms_disp_mu_tperm(state.tperm.as_mut_ptr(), nperm, 1) };
            unsafe { ms_disp_disp_stress_norhs() };
        } else {
            unsafe { ms_disp_mu_meshperm(1) };
        }
    }
}

/// The `LocalMatrix` internal procedure: element stiffness matrix by
/// numerical integration. `n` is the element node count, `nd`/`nb` the
/// total/bubble DOF counts. `ld` is the SAVE leading dimension.
fn local_matrix(state: &mut MeshSolveState, ld: c_int, buf_n: c_int, dim: c_int, n: c_int, nd: c_int, nb: c_int) {
    unsafe { ms_lm_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    // Basis(nd) / dBasisdx(nd,3) automatics, sized by this element.
    let nn = n as usize;
    let ndn = nd as usize;
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];
    // NodalLame1(n) / NodalLame2(n) automatics.
    let mut lame1 = vec![0.0 as c_double; nn];
    let mut lame2 = vec![0.0 as c_double; nn];
    unsafe {
        ms_lame(
            n,
            buf_n,
            state.emod.as_mut_ptr(),
            state.pr.as_mut_ptr(),
            lame1.as_mut_ptr(),
            lame2.as_mut_ptr(),
            0,
        )
    };

    let nip = unsafe { ms_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { ms_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            ms_lm_eleminfo(
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

        // Lame parameters at the integration point (SUM reductions stay
        // in the wrapper, as in the original).
        let l1: c_double =
            unsafe { ms_dot(n, lame1.as_ptr(), basis.as_ptr()) };
        let l2: c_double =
            unsafe { ms_dot(n, lame2.as_ptr(), basis.as_ptr()) };

        unsafe {
            ms_lm_assemble(
                ld,
                dim,
                nd,
                s * detj,
                l1,
                l2,
                dbasis.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
            )
        };
    }

    // Symmetric fill + bubble-DOF clamp (wrapper transcribes exactly).
    unsafe {
        ms_bubble(
            ld,
            dim,
            nd,
            n,
            nb,
            state.stiff.as_mut_ptr(),
            state.force.as_mut_ptr(),
        )
    };
}

/// The `MeshBoundary` internal procedure: Neumann/penalty boundary
/// assembly by numerical integration.
fn mesh_boundary(state: &mut MeshSolveState, ld: c_int, _buf_n: c_int, dim: c_int, n: c_int, nd: c_int, nb: c_int) {
    unsafe { ms_mb_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    let ndn = nd as usize;
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { ms_mb_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { ms_mb_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            ms_mb_eleminfo(
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

        unsafe {
            ms_mb_assemble(
                ld,
                dim,
                n,
                nd,
                u,
                v,
                w,
                s * detj,
                state.load.as_mut_ptr(),
                state.alpha.as_mut_ptr(),
                state.beta.as_mut_ptr(),
                basis.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }

    unsafe {
        ms_bubble(
            ld,
            dim,
            nd,
            n,
            nb,
            state.stiff.as_mut_ptr(),
            state.force.as_mut_ptr(),
        )
    };
}
