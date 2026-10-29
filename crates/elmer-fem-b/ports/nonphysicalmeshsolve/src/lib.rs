//! Port of `fem/src/modules/NonphysicalMeshSolve.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/NonphysicalMeshSolve.F90` — one F77-global
//! entry `meshsolver_` with three `CONTAINS` internals (`LocalMatrix`,
//! `MeshBoundary`, `NodalDisplacementPenalty`). The emitted dynamic symbol
//! is exactly `meshsolver_` (proven by `nm -D` on the reference `.so`);
//! the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses; `dt` and the
//! transient flag are unused by the body, exactly as in the original.
//! Every component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/nonphysicalmeshsolve_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`), the `(p,q,i,j)` assembly nests, the mesh-pointer swaps, the
//! `MovingMesh` coordinate update, `NormalVector`/`GetParentUVW`, `Info`
//! output, and the `CRS_AddToMatrixElement` penalty update stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors NonphysicalMeshSolve.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: the numeric SAVE buffers
//!   (`STIFF`/`LOAD`/`FORCE`/`ElasticModulus`/`PoissonRatio`/`Alpha`/
//!   `Beta`/`Gamma`/`RefSurface`) live in the Rust state sized at
//!   first-mesh sizing (`STDOFs*MaxElementDOFs`, never resized on valid
//!   inputs); the pointer-identity SAVE arrays (`OrigX/Y/Z`,
//!   `TrueX/Y/Z`, `PrevMeshUpdate`), `VisitedTimes`, `dim`, the
//!   target-field handles, the per-procedure `SAVE Nodes`/`ParentBasis`
//!   states, and the penalty `Visited` flag live wrapper-side with
//!   identical lifetimes (they must alias live mesh storage, which only
//!   Fortran can hold).
//! - Allocation failure calls the same `Fatal('NonPhysicalMeshSolve',
//!   ...)` via `try_reserve` (the `STAT=istat` check), never abort
//!   directly.
//! - The dead `IF (AllocationsDone) DEALLOCATE` inside the
//!   `IF (.NOT. AllocationsDone)` block is not transcribed (it can never
//!   execute, as in the original).
//! - `Model`, `dt`, `TransientSimulation` are received and ignored by the
//!   body except for the `VariableGet(Model % Variables, ...)` target
//!   lookups and the `ListCheckPresent*` mask probes, exactly as in the
//!   original.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn nma_begin(h_model: *const c_void, h_solver: *mut c_void);
    fn nma_info_start();
    fn nma_assoc_matrix() -> c_int;
    fn nma_visit();
    fn nma_frame();
    fn nma_stdofs() -> c_int;
    fn nma_dim() -> c_int;
    fn nma_local_nodes() -> c_int;
    fn nma_max_elem_dofs() -> c_int;
    fn nma_cumulative() -> c_int;
    fn nma_first_alloc();
    fn nma_nodes_capture();
    fn nma_nodes_to_orig();
    fn nma_moving_update();
    fn nma_fatal_alloc();
    fn nma_assembly_begin() -> c_int;
    fn nma_advance(t: c_int, noactive: c_int);
    fn nma_begin_active(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int);
    fn nma_mat_props(n: c_int, buf_n: c_int, emod: *mut c_double, pr: *mut c_double);
    fn nma_lm_nodes();
    #[allow(clippy::too_many_arguments)]
    fn nma_lame(
        n: c_int,
        buf_n: c_int,
        emod: *mut c_double,
        pr: *mut c_double,
        lame1: *mut c_double,
        lame2: *mut c_double,
        plane: c_int,
    );
    fn nma_lm_coeff() -> c_double;
    fn nma_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn nma_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn nma_lm_eleminfo(
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
    fn nma_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    #[allow(clippy::too_many_arguments)]
    fn nma_lm_assemble(
        ld: c_int,
        sdim: c_int,
        nd: c_int,
        s: c_double,
        l1: c_double,
        l2: c_double,
        coeff: c_double,
        dbasis: *mut c_double,
        stiff: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn nma_bubble(
        ld: c_int,
        sdim: c_int,
        nd: c_int,
        n: c_int,
        nb: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn nma_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn nma_finish_bulk();
    fn nma_target_setup(got_field: *mut c_int, got_surface: *mut c_int, got_grad: *mut c_int);
    fn nma_nbnd() -> c_int;
    fn nma_begin_boundary(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        nb: *mut c_int,
        active: *mut c_int,
    );
    fn nma_bc() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn nma_bc_vals(
        n: c_int,
        buf_n: c_int,
        load: *mut c_double,
        alpha: *mut c_double,
        beta: *mut c_double,
        gamma: *mut c_double,
        refsurf: *mut c_double,
    ) -> c_int;
    fn nma_mb_nodes();
    fn nma_mb_sdim() -> c_int;
    fn nma_mb_parent_setup() -> c_int;
    fn nma_mb_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn nma_mb_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn nma_mb_eleminfo(
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
    fn nma_mb_assemble(
        ld: c_int,
        n: c_int,
        nd: c_int,
        nb: c_int,
        u: c_double,
        v: c_double,
        w: c_double,
        s: c_double,
        load: *mut c_double,
        alpha: *mut c_double,
        beta: *mut c_double,
        gamma: *mut c_double,
        refsurf: *mut c_double,
        basis: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn nma_finish_assembly();
    fn nma_nodal_penalty();
    fn nma_dirichlet();
    fn nma_default_solve() -> c_double;
    fn nma_nonlin_relax();
    fn nma_mesh_relax() -> c_double;
    fn nma_displace(relax: c_double);
}

/// Grow the SAVE buffers when an element needs more than the first-mesh
/// sizing (adaptive remesh across invocations). The original never grows
/// (its comment assumes the first mesh is big enough); growth only
/// triggers where the original would overflow silently, and the leading
/// sections used by the wrappers are unaffected, so valid inputs never
/// observe it.
fn ensure_state(state: &mut NonphysicalMeshState, ld: usize, n: usize) {
    if ld <= state.ld && n <= state.n {
        return;
    }
    let ld = ld.max(state.ld);
    let n = n.max(state.n);
    state.stiff.resize(ld * ld, 0.0);
    state.force.resize(ld, 0.0);
    state.load.resize(4 * n, 0.0);
    state.emod.resize(n, 0.0);
    state.pr.resize(n, 0.0);
    state.alpha.resize(3 * n, 0.0);
    state.beta.resize(n, 0.0);
    state.gamma.resize(n, 0.0);
    state.refsurf.resize(n, 0.0);
    state.ld = ld;
    state.n = n;
}

/// First-mesh SAVE state: `STIFF(LD,LD)`, `FORCE(LD)`, `LOAD(4,N)`,
/// `ElasticModulus(N)`, `PoissonRatio(N)`, `Alpha(3,N)`, `Beta(N)`,
/// `Gamma(N)`, `RefSurface(N)` plus the `AllocationsDone` flag. `N` is
/// the first call's `MaxElementDOFs`; `LD = STDOFs*N`.
struct NonphysicalMeshState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    load: Vec<c_double>,
    emod: Vec<c_double>,
    pr: Vec<c_double>,
    alpha: Vec<c_double>,
    beta: Vec<c_double>,
    gamma: Vec<c_double>,
    refsurf: Vec<c_double>,
    ld: usize,
    n: usize,
    allocations_done: bool,
}

impl NonphysicalMeshState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            load: Vec::new(),
            emod: Vec::new(),
            pr: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            gamma: Vec::new(),
            refsurf: Vec::new(),
            ld: 0,
            n: 0,
            allocations_done: false,
        }
    }

    /// Reserve (but do not yet fill) the SAVE buffers; on failure calls
    /// the same `Fatal('NonPhysicalMeshSolve', ...)` as the original's
    /// `STAT=istat` check.
    fn reserve(&mut self, ld: usize, n: usize) -> bool {
        if self.stiff.try_reserve(ld * ld).is_err()
            || self.force.try_reserve(ld).is_err()
            || self.load.try_reserve(4 * n).is_err()
            || self.emod.try_reserve(n).is_err()
            || self.pr.try_reserve(n).is_err()
            || self.alpha.try_reserve(3 * n).is_err()
            || self.beta.try_reserve(n).is_err()
            || self.gamma.try_reserve(n).is_err()
            || self.refsurf.try_reserve(n).is_err()
        {
            return false;
        }
        self.stiff.resize(ld * ld, 0.0);
        self.force.resize(ld, 0.0);
        self.load.resize(4 * n, 0.0);
        self.emod.resize(n, 0.0);
        self.pr.resize(n, 0.0);
        self.alpha.resize(3 * n, 0.0);
        self.beta.resize(n, 0.0);
        self.gamma.resize(n, 0.0);
        self.refsurf.resize(n, 0.0);
        self.ld = ld;
        self.n = n;
        true
    }
}
static STATE: Mutex<NonphysicalMeshState> = Mutex::new(NonphysicalMeshState::new());

/// `meshsolver_` — the nonphysical mesh update solver entry. Argument
/// order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body except for target-variable
/// lookups, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads if non-null.
#[export_name = "meshsolver_"]
pub unsafe extern "C" fn meshsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { nma_begin(model, solver) };
    unsafe { nma_info_start() };

    // IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { nma_assoc_matrix() } == 0 {
        return;
    }

    // VisitedTimes = VisitedTimes + 1; frame capture (dim, MeshSol and
    // its Perm/Values pointers); LocalNodes = COUNT(MeshPerm > 0).
    unsafe { nma_visit() };
    unsafe { nma_frame() };
    let local_nodes = unsafe { nma_local_nodes() };
    // IF (LocalNodes <= 0) RETURN
    if local_nodes <= 0 {
        return;
    }
    unsafe { nma_cumulative() };

    // Allocate permanent storage, first call only.
    if !state.allocations_done {
        let n = unsafe { nma_max_elem_dofs() } as usize;
        let stdofs = unsafe { nma_stdofs() } as usize;
        let ld = stdofs * n;
        if !state.reserve(ld, n) {
            unsafe { nma_fatal_alloc() };
            std::process::abort();
        }
        unsafe { nma_first_alloc() };
        state.allocations_done = true;
    }
    let sdim = unsafe { nma_dim() };

    // Node-pointer capture and Orig swap; MovingMesh block.
    unsafe { nma_nodes_capture() };
    unsafe { nma_nodes_to_orig() };
    unsafe { nma_moving_update() };

    // Assembly driver.
    let noactive = unsafe { nma_assembly_begin() };
    for t in 1..=noactive {
        unsafe { nma_advance(t, noactive) };
        let (mut n, mut nd, mut nb) = (0, 0, 0);
        unsafe { nma_begin_active(t, &mut n, &mut nd, &mut nb) };
        // First-mesh sizing holds on valid inputs; grow-only guard
        // covers remeshed re-entry (invisible otherwise).
        ensure_state(&mut state, (sdim as usize) * (nd as usize), n as usize);
        let (ld, buf_n) = (state.ld as c_int, state.n as c_int);
        unsafe { nma_mat_props(n, buf_n, state.emod.as_mut_ptr(), state.pr.as_mut_ptr()) };
        local_matrix(&mut state, ld, buf_n, sdim, n, nd, nb);
        unsafe {
            nma_default_update(ld, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }
    let (ld, buf_n) = (state.ld as c_int, state.n as c_int);
    unsafe { nma_finish_bulk() };

    // Target field/surface setup (Fatal paths live in the wrapper).
    let (mut got_field, mut got_surface, mut got_grad) = (0, 0, 0);
    unsafe { nma_target_setup(&mut got_field, &mut got_surface, &mut got_grad) };
    if got_surface == 0 {
        state.refsurf.fill(0.0);
    }

    // Neumann & Newton boundary conditions.
    let nbnd = unsafe { nma_nbnd() };
    for t in 1..=nbnd {
        let (mut n, mut nd, mut nb, mut active) = (0, 0, 0, 0);
        unsafe { nma_begin_boundary(t, &mut n, &mut nd, &mut nb, &mut active) };
        if active == 0 {
            continue;
        }
        if unsafe { nma_bc() } == 0 {
            continue;
        }
        state.load.fill(0.0);
        state.alpha.fill(0.0);
        state.beta.fill(0.0);
        state.gamma.fill(0.0);
        state.refsurf.fill(0.0);
        let got = unsafe {
            nma_bc_vals(
                n,
                buf_n,
                state.load.as_mut_ptr(),
                state.alpha.as_mut_ptr(),
                state.beta.as_mut_ptr(),
                state.gamma.as_mut_ptr(),
                state.refsurf.as_mut_ptr(),
            )
        };
        if got == 0 {
            continue;
        }
        mesh_boundary(&mut state, ld, buf_n, got_surface, got_grad, n, nd, nb);
        unsafe {
            nma_default_update(ld, state.stiff.as_mut_ptr(), state.force.as_mut_ptr());
        }
    }

    unsafe { nma_finish_assembly() };

    // Nodal penalty, Dirichlet BCs, solve.
    unsafe { nma_nodal_penalty() };
    unsafe { nma_dirichlet() };

    // UNorm = DefaultSolve() (local and unused, as in the original).
    let _unorm: c_double = unsafe { nma_default_solve() };

    // Relaxation factors and the mesh coordinate update (with the
    // node-pointer restore for the non-cumulative case).
    unsafe { nma_nonlin_relax() };
    let relax: c_double = unsafe { nma_mesh_relax() };
    unsafe { nma_displace(relax) };
}

/// The `LocalMatrix` internal procedure: element stiffness matrix by
/// numerical integration. `n` is the element node count, `nd`/`nb` the
/// total/bubble DOF counts. `ld` is the SAVE leading dimension and
/// `buf_n` the SAVE node dimension.
fn local_matrix(
    state: &mut NonphysicalMeshState,
    ld: c_int,
    buf_n: c_int,
    sdim: c_int,
    n: c_int,
    nd: c_int,
    nb: c_int,
) {
    unsafe { nma_lm_nodes() };
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
        nma_lame(
            n,
            buf_n,
            state.emod.as_mut_ptr(),
            state.pr.as_mut_ptr(),
            lame1.as_mut_ptr(),
            lame2.as_mut_ptr(),
            0,
        )
    };
    let coeff: c_double = unsafe { nma_lm_coeff() };

    let nip = unsafe { nma_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { nma_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            nma_lm_eleminfo(
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
            unsafe { nma_dot(n, lame1.as_ptr(), basis.as_ptr()) };
        let l2: c_double =
            unsafe { nma_dot(n, lame2.as_ptr(), basis.as_ptr()) };

        unsafe {
            nma_lm_assemble(
                ld,
                sdim,
                nd,
                s * detj,
                l1,
                l2,
                coeff,
                dbasis.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
            )
        };
    }

    // Symmetric fill + bubble-DOF clamp (wrapper transcribes exactly).
    unsafe {
        nma_bubble(
            ld,
            sdim,
            nd,
            n,
            nb,
            state.stiff.as_mut_ptr(),
            state.force.as_mut_ptr(),
        )
    };
}

/// The `MeshBoundary` internal procedure: Neumann/target boundary
/// assembly by numerical integration. The per-IP physics (normal force,
/// target-field/surface penalty) stays in the wrapper; `nb` selects the
/// bubble-DOF clamp.
#[allow(clippy::too_many_arguments)]
fn mesh_boundary(
    state: &mut NonphysicalMeshState,
    ld: c_int,
    _buf_n: c_int,
    got_surface: c_int,
    got_grad: c_int,
    n: c_int,
    nd: c_int,
    nb: c_int,
) {
    unsafe { nma_mb_nodes() };
    state.stiff.fill(0.0);
    state.force.fill(0.0);

    // Parent setup mirrors the early-RETURN block: a missing parent
    // skips the whole element, including the bubble clamp.
    if got_surface != 0 && got_grad == 0 && unsafe { nma_mb_parent_setup() } == 0 {
        return;
    }
    let msdim = unsafe { nma_mb_sdim() };

    let ndn = nd as usize;
    let mut basis = vec![0.0 as c_double; ndn];
    let mut dbasis = vec![0.0 as c_double; ndn * 3];

    let nip = unsafe { nma_mb_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { nma_mb_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            nma_mb_eleminfo(
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
            nma_mb_assemble(
                ld,
                n,
                nd,
                nb,
                u,
                v,
                w,
                s * detj,
                state.load.as_mut_ptr(),
                state.alpha.as_mut_ptr(),
                state.beta.as_mut_ptr(),
                state.gamma.as_mut_ptr(),
                state.refsurf.as_mut_ptr(),
                basis.as_mut_ptr(),
                state.stiff.as_mut_ptr(),
                state.force.as_mut_ptr(),
            )
        };
    }

    unsafe {
        nma_bubble(
            ld,
            msdim,
            nd,
            n,
            nb,
            state.stiff.as_mut_ptr(),
            state.force.as_mut_ptr(),
        )
    };
}
