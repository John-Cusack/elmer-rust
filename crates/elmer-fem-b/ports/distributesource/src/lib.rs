//! Port of `fem/src/modules/DistributeSource.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/DistributeSource.F90` — one F77-global
//! entry `distributesource_` with one `CONTAINS` internal procedure. The
//! emitted dynamic symbol is exactly `distributesource_` (proven by `nm
//! -D` on the reference `.so`); the internal emits only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/distributesource_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MATMUL`), formatted I/O (`PRINT`) and whole-array updates
//! stay in Fortran wrappers, so bit-identity holds by construction
//! wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors DistributeSource.F90 statement for
//! statement):
//! - The entry's `MeshLocal`/`MeshGlobal` handles live in the wrapper;
//!   the `ASSOCIATED` check runs there and `Fatal`s identically.
//! - The internal procedure's dummy `LocalSourceName` (`'Local Source'`)
//!   is never referenced by the body, so no wrapper carries it; the
//!   `GlobalRhsName` (`'g'`) and `'g fix'` literals live in the wrapper
//!   module exactly where the `VariableGet`/`Fatal`/`Warn` statements
//!   use them. `UseQuadrantTree` (`.TRUE.` at the single call site) is
//!   folded into the tree-setup wrapper.
//! - `SAVE` semantics replicated explicitly: `ElementNodes` /
//!   `ElementNodesB` live as wrapper-side SAVE state with identical
//!   lifetimes (allocated once at first-mesh sizing, grown only if a
//!   later mesh needs more — growth triggers only where the original
//!   would overflow, so valid inputs never observe it). All other
//!   internal arrays are per-call locals here, as in the original.
//! - Dead locals of the original (`LocalSource`, `Eps1`/`Eps2`,
//!   `eps_tries`, `eps_numeric`, the `l` index) are omitted: the
//!   `ListGet` reads are kept (same call sequence, no side effects) but
//!   their results are unused, and `PointInElement` is called without
//!   the optional epsilon args exactly as in the original, so the same
//!   library defaults apply.
//! - The `ub`/`vb` swap (`vb = LocalCoordinates(1)`,
//!   `ub = LocalCoordinates(2)`) is kept verbatim.
//! - Both `ElementInfo` function results overwrite the (never reread)
//!   `Stat` local, as in the original.
//! - No Rust threads: the solver entries run serially, and this port
//!   keeps no Rust-side `static` state at all (all SAVE state is
//!   wrapper-side), so no sharing guard is needed.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn dsa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn dsa_meshes_differ() -> c_int;
    fn dsa_fatal_same_mesh();
    fn dsa_info_projecting();
    fn dsa_edge_param(val: *mut c_int, found: *mut c_int);
    fn dsa_dim() -> c_int;
    fn dsa_ensure_qtree() -> c_int;
    fn dsa_max_dofs() -> c_int;
    fn dsa_alloc_nodes(n: c_int);
    fn dsa_sim_params(
        geps: *mut c_double,
        leps: *mut c_double,
        tries: *mut c_int,
        neps: *mut c_double,
    );
    fn dsa_find_rhs() -> c_int;
    fn dsa_fatal_no_rhs();
    fn dsa_assoc_gsolver() -> c_int;
    fn dsa_fatal_no_gsolver();
    fn dsa_edge_detect() -> c_int;
    fn dsa_info_edge();
    fn dsa_find_fix() -> c_int;
    fn dsa_warn_no_fix();
    fn dsa_edge_style(degree: *mut c_int);
    fn dsa_fatal_no_edge();
    fn dsa_zero_values();
    fn dsa_info_integrating();
    fn dsa_local_num_bulk() -> c_int;
    fn dsa_begin_local(i: c_int, n: *mut c_int);
    fn dsa_gather_local(n: c_int);
    fn dsa_bodyforce() -> c_int;
    fn dsa_curr_dens(n: c_int, a: *mut c_double, found: *mut c_int);
    fn dsa_ip_begin(edge: c_int) -> c_int;
    fn dsa_ip_point(
        t: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn dsa_eleminfo_local(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        ddx: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn dsa_ip_eval(
        n: c_int,
        basis: *mut c_double,
        curr: *mut c_double,
        point: *mut c_double,
        currip: *mut c_double,
    );
    fn dsa_find_leaf(px: c_double, py: c_double, pz: c_double, dim: c_int) -> c_int;
    fn dsa_leaf_count() -> c_int;
    fn dsa_begin_leaf(k: c_int, nb: *mut c_int, ne: *mut c_int);
    fn dsa_global_num_bulk() -> c_int;
    fn dsa_begin_global(k: c_int, nb: *mut c_int, ne: *mut c_int);
    fn dsa_edge_idx(ne: c_int, idx: *mut c_int);
    fn dsa_gather_global(nb: c_int);
    fn dsa_point_in_elem(point: *mut c_double, lc: *mut c_double) -> c_int;
    fn dsa_info_noelem(i: c_int);
    fn dsa_rhs_ind(ne: c_int, idx: *mut c_int, rhs: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn dsa_add_rhs(
        ne: c_int,
        rhs: *mut c_int,
        wb: *mut c_double,
        currip: *mut c_double,
        detj: c_double,
        s: c_double,
    );
    fn dsa_fix_assoc() -> c_int;
    fn dsa_add_fix(
        nb: c_int,
        ddxb: *mut c_double,
        currip: *mut c_double,
        detj: c_double,
        s: c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn dsa_eleminfo_global(
        ub: c_double,
        vb: c_double,
        wb: c_double,
        nb: c_int,
        ne: c_int,
        basisb: *mut c_double,
        ddxb: *mut c_double,
        wb0: *mut c_double,
        rotw: *mut c_double,
        detjb: *mut c_double,
    );
    fn dsa_finish(nofails: c_int, nofound: c_int);
    fn dsa_fatal_alloc();
}

/// Grow helper for the per-call work arrays (first-mesh sizing; growth
/// only triggers where the original's fixed-size `ALLOCATE` would
/// already have been exceeded).
fn reserve(buf: &mut Vec<c_double>, need: usize) {
    if buf.try_reserve(need.saturating_sub(buf.len())).is_err() {
        unsafe { dsa_fatal_alloc() };
        std::process::abort();
    }
    buf.resize(need, 0.0);
}

fn reserve_int(buf: &mut Vec<c_int>, need: usize) {
    if buf.try_reserve(need.saturating_sub(buf.len())).is_err() {
        unsafe { dsa_fatal_alloc() };
        std::process::abort();
    }
    buf.resize(need, 0);
}

/// `distributesource_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, Transient)`; `dt` and the transient flag are
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt`/`transient` must be valid reads if non-null (they
/// are ignored).
#[export_name = "distributesource_"]
pub unsafe extern "C" fn distribute_source_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { dsa_begin(model, solver) };

    // IF( ASSOCIATED( MeshLocal, MeshGlobal ) ) CALL Fatal(...).
    if unsafe { dsa_meshes_differ() } == 0 {
        unsafe { dsa_fatal_same_mesh() };
    }

    // CALL InterpolateAndIntegrateSource( Model, Solver, MeshLocal,
    //   MeshGlobal, 'Local Source', 'g', .TRUE. ).
    unsafe { interpolate() };
}

/// The `InterpolateAndIntegrateSource` internal procedure: projects the
/// `'Local Current Density'` body force of the local mesh onto the
/// global mesh edge-basis right-hand side `'g'` (plus the `'g fix'`
/// potential-fixing term) by Gauss integration over the local bulk
/// elements with per-point element search in the global mesh.
///
/// # Safety
///
/// The entry must have run `dsa_begin` first (wrapper state holds the
/// model/solver/mesh handles).
unsafe fn interpolate() {
    // CALL Info(Caller,'Projecting current sources between meshes').
    unsafe { dsa_info_projecting() };

    // Params => GetSolverParams(Solver); EdgeBasis = GetLogical(...).
    // (The GetLogical value is overwritten below, exactly as in the
    // original; the call itself is kept for sequence fidelity.)
    let (mut edge_param, mut pfound) = (0, 0);
    unsafe { dsa_edge_param(&mut edge_param, &mut pfound) };
    let _ = (edge_param, pfound);

    // dim = CoordinateSystemDimension().
    let dim = unsafe { dsa_dim() };

    // Quadrant-tree setup (UseQuadrantTree = .TRUE. at the call site).
    let tryqtree = unsafe { dsa_ensure_qtree() };

    // First-mesh sizing for the SAVE Nodes mirrors and the work arrays.
    let n = unsafe { dsa_max_dofs() } as usize;
    unsafe { dsa_alloc_nodes(n as c_int) };

    let mut basis = Vec::new();
    let mut basis_b = Vec::new();
    let mut curr = Vec::new();
    let mut wbasis = Vec::new();
    let mut rotw = Vec::new();
    let mut ddx = Vec::new();
    let mut ddxb = Vec::new();
    let mut rhs_ind = Vec::new();
    reserve(&mut basis, n);
    reserve(&mut basis_b, n);
    reserve(&mut curr, 3 * n);
    reserve(&mut wbasis, 3 * n);
    reserve(&mut rotw, 3 * n);
    reserve(&mut ddx, 3 * n);
    reserve(&mut ddxb, 3 * n);
    reserve_int(&mut rhs_ind, n);
    // Fixed-size edge-index scratch, as `IndexesB(100)`.
    let mut idxb = [0 as c_int; 100];

    // Dead simulation-knob reads kept for call-sequence fidelity (their
    // values are unused downstream, exactly as in the original).
    let (mut geps, mut leps, mut neps) = (0.0, 0.0, 0.0);
    let mut tries = 0;
    unsafe { dsa_sim_params(&mut geps, &mut leps, &mut tries, &mut neps) };
    let _ = (geps, leps, tries, neps);

    // RhsSol => VariableGet( GlobalMesh % Variables, 'g', .TRUE. ).
    if unsafe { dsa_find_rhs() } == 0 {
        unsafe { dsa_fatal_no_rhs() };
    }
    // Associate the r.h.s. vector with its owning solver.
    if unsafe { dsa_assoc_gsolver() } == 0 {
        unsafe { dsa_fatal_no_gsolver() };
    }

    // Edge-basis detection (EdgeBasis/PiolaT start .FALSE.).
    let edge = unsafe { dsa_edge_detect() } != 0;
    if edge {
        unsafe { dsa_info_edge() };
        // FixSol => VariableGet( GlobalMesh % Variables, 'g fix', .TRUE. ).
        if unsafe { dsa_find_fix() } == 0 {
            unsafe { dsa_warn_no_fix() };
        }
        // CALL EdgeElementStyle( GSolver % Values, PiolaT, BasisDegree ).
        let mut _degree = 0;
        unsafe { dsa_edge_style(&mut _degree) };
    } else {
        unsafe { dsa_fatal_no_edge() };
    }

    // Initialize the source term and fixing potential source term.
    unsafe { dsa_zero_values() };

    // CALL Info(Caller,'Integrating over bulk elements of local mesh',...).
    unsafe { dsa_info_integrating() };

    let (mut nofails, mut nofound) = (0, 0);
    let nlocal = unsafe { dsa_local_num_bulk() };
    for i in 1..=nlocal {
        let mut nn = 0;
        unsafe { dsa_begin_local(i, &mut nn) };
        // Model % CurrentElement => Element + gather element nodes.
        unsafe { dsa_gather_local(nn) };

        // BodyForce => GetBodyForce(Element); IF(.NOT.ASSOCIATED) CYCLE.
        if unsafe { dsa_bodyforce() } == 0 {
            continue;
        }

        // CALL GetRealVector( BodyForce, LocalCurrDens(1:3,1:n), ... ).
        let mut found = 0;
        unsafe { dsa_curr_dens(nn, curr.as_mut_ptr(), &mut found) };
        if found == 0 {
            continue;
        }

        let nip = unsafe { dsa_ip_begin(1) };
        for j in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { dsa_ip_point(j, &mut u, &mut v, &mut w, &mut s) };

            // Stat = ElementInfo( Element, ElementNodes, u, v, w, ... ).
            // (Stat is never reread, as in the original.)
            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                dsa_eleminfo_local(
                    u,
                    v,
                    w,
                    nn,
                    basis.as_mut_ptr(),
                    ddx.as_mut_ptr(),
                    &mut detj,
                    &mut stat,
                )
            };
            let _ = stat;

            // Nodal point + current density at the integration point.
            let mut point = [0.0 as c_double; 3];
            let mut currip = [0.0 as c_double; 3];
            unsafe {
                dsa_ip_eval(
                    nn,
                    basis.as_mut_ptr(),
                    curr.as_mut_ptr(),
                    point.as_mut_ptr(),
                    currip.as_mut_ptr(),
                )
            };

            // Find which global-mesh bulk element holds the point.
            let (mut nb, mut ne) = (0, 0);
            let mut lc = [0.0 as c_double; 3];
            let mut found_el = 0;
            if tryqtree != 0 {
                if unsafe { dsa_find_leaf(point[0], point[1], point[2], dim) } != 0 {
                    // Octree leaf path: only the leaf's elements.
                    let nk = unsafe { dsa_leaf_count() };
                    for k in 1..=nk {
                        unsafe { dsa_begin_leaf(k, &mut nb, &mut ne) };
                        debug_assert!((ne as usize) <= idxb.len());
                        unsafe { dsa_edge_idx(ne, idxb.as_mut_ptr()) };
                        unsafe { dsa_gather_global(nb) };
                        found_el = unsafe {
                            dsa_point_in_elem(point.as_mut_ptr(), lc.as_mut_ptr())
                        };
                        if found_el != 0 {
                            break;
                        }
                    }
                } else {
                    // Dummy N^2 search over all global bulk elements.
                    let ng = unsafe { dsa_global_num_bulk() };
                    for k in 1..=ng {
                        unsafe { dsa_begin_global(k, &mut nb, &mut ne) };
                        debug_assert!((ne as usize) <= idxb.len());
                        unsafe { dsa_edge_idx(ne, idxb.as_mut_ptr()) };
                        unsafe { dsa_gather_global(nb) };
                        found_el = unsafe {
                            dsa_point_in_elem(point.as_mut_ptr(), lc.as_mut_ptr())
                        };
                        if found_el != 0 {
                            break;
                        }
                    }
                }
            }

            if found_el == 0 {
                unsafe { dsa_info_noelem(i) };
                nofails += 1;
                continue;
            }

            // RhsInd(1:ne) = RhsSol % Perm( IndexesB(1:ne) ).
            unsafe { dsa_rhs_ind(ne, idxb.as_mut_ptr(), rhs_ind.as_mut_ptr()) };
            nofound += 1;

            // Note the index swap, verbatim from the original:
            // vb = LocalCoordinates(1), ub = LocalCoordinates(2).
            let (vb, ub, wb) = (lc[0], lc[1], lc[2]);

            // Evaluate basis functions at the integration point
            // (EdgeBasis = WBasis, RotBasis = RotWBasis, USolver = GSolver).
            let mut detjb = 0.0;
            unsafe {
                dsa_eleminfo_global(
                    ub,
                    vb,
                    wb,
                    nb,
                    ne,
                    basis_b.as_mut_ptr(),
                    ddxb.as_mut_ptr(),
                    wbasis.as_mut_ptr(),
                    rotw.as_mut_ptr(),
                    &mut detjb,
                )
            };

            // Note that we use the metric of the local mesh but populate
            // it on the global mesh.
            unsafe {
                dsa_add_rhs(
                    ne,
                    rhs_ind.as_mut_ptr(),
                    wbasis.as_mut_ptr(),
                    currip.as_mut_ptr(),
                    detj,
                    s,
                )
            };

            if unsafe { dsa_fix_assoc() } != 0 {
                unsafe {
                    dsa_add_fix(
                        nb,
                        ddxb.as_mut_ptr(),
                        currip.as_mut_ptr(),
                        detj,
                        s,
                    )
                };
            }
        }
    }

    unsafe { dsa_finish(nofails, nofound) };
}