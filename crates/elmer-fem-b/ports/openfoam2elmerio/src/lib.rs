//! Port of `fem/src/modules/OpenFoam2ElmerIO.F90` (Elmer @ the lane fork)
//! for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/OpenFoam2ElmerIO.F90` — one F77-global
//! entry `openfoam2elmerfit_` with five `CONTAINS` internals. The emitted
//! dynamic symbol is exactly `openfoam2elmerfit_` (proven by `nm -D` on
//! the reference `.so`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/openfoam2elmerio_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`MINVAL`/`MAXVAL`/`SUM`, `MATMUL`), formatted I/O
//! (`WRITE(Message,'(A,ES12.5)')`, OpenFOAM file reads, `SystemCommand`)
//! and whole-array updates stay in Fortran wrappers, so bit-identity
//! holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors OpenFoam2ElmerIO.F90 statement for statement):
//! - `SAVE` semantics replicated: entry `OFMesh`/`NoDir`/`Visited` live as
//!   wrapper-side SAVE state plus a Rust `visited` flag with identical
//!   lifetimes (second visit releases the saved mesh before reallocating,
//!   exactly as the original); `DiffusionAssembly`'s `SAVE` arrays
//!   (`STIFF`/`FORCE`/`DataDiffusivity`) are first-mesh sized once and
//!   never resized; the internal procedures' `SAVE Nodes` live as
//!   wrapper-side SAVE state with identical lifetimes (consumed before
//!   the next fill, as in the original).
//! - The file-reading internals (`OpenFOAMBlocks`, `CreateFOAMMesh`,
//!   `ReadFOAMField`) are transcribed whole into wrappers — same
//!   compiler, same call sequence — so their bytes (including the
//!   `GOTO 10` early exits and `I2S` renderings) match by construction on
//!   all defined inputs.
//! - Allocation failure calls the same `Fatal('OpenFoam2ElmerFit', ...)`
//!   via `try_reserve`, instead of aborting.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body, exactly
//!   as in the original (received and ignored); `Norm = DefaultSolve()`
//!   is kept and discarded.
//! - DataAssembly's per-point `Rhs`/`WeightVector` scatter and the fit
//!   scaling (`Rhs = Coeff*Rhs`, `Values(Diag) = Coeff*Weight`) run in
//!   Fortran wrappers in the original loop order, so accumulation order
//!   matches exactly.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`es
//!   are only `static`-sharing guards, never contended across threads.
use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn ofa_begin(h_solver: *const c_void);
    fn ofa_info_head();
    fn ofa_var_check();
    fn ofa_rhs_len() -> c_int;
    fn ofa_release_saved();
    fn ofa_alloc_mesh();
    fn ofa_passive(found: *mut c_int);
    fn ofa_passive_val() -> c_int;
    fn ofa_cylsymm() -> c_int;
    fn ofa_meshdim() -> c_int;
    fn ofa_warn_reduced();
    fn ofa_basedir(found: *mut c_int) -> c_int;
    fn ofa_have_mesh1() -> c_int;
    fn ofa_blocks();
    fn ofa_parallel_reduce();
    fn ofa_nodir() -> c_int;
    fn ofa_info_nblocks();
    fn ofa_default_initialize();
    fn ofa_weight_alloc(n: c_int);
    fn ofa_weight_zero();
    fn ofa_weight_free();
    fn ofa_mesh_dir(i: c_int, found: *mut c_int);
    fn ofa_info_block_head(i: c_int);
    fn ofa_info_reading();
    fn ofa_create_mesh();
    fn ofa_field_alloc();
    fn ofa_field_free();
    fn ofa_release_contents();
    fn ofa_read_field(got: *mut c_int);
    fn ofa_report_stats();
    fn ofa_ofcount() -> c_int;
    fn ofa_ofcoords(t: c_int, cyl: c_int, g: *mut c_double);
    fn ofa_field_val(t: c_int) -> c_double;
    fn ofa_locate(g: *mut c_double, l: *mut c_double, idx: *mut c_int);
    fn ofa_data_begin_elem(idx: c_int, n: *mut c_int);
    fn ofa_data_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn ofa_data_scatter(val: c_double, n: c_int, basis: *mut c_double);
    fn ofa_fit_coeff(found: *mut c_int) -> c_double;
    fn ofa_apply_fit(coeff: c_double);
    fn ofa_max_elem_nodes() -> c_int;
    fn ofa_diff_coeff(glob: *mut c_int) -> c_double;
    fn ofa_nof_active() -> c_int;
    fn ofa_coord_dim() -> c_int;
    fn ofa_begin_active(t: c_int, n: *mut c_int);
    fn ofa_diff_nodes();
    fn ofa_gauss_begin() -> c_int;
    fn ofa_gauss_point(k: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn ofa_diff_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn ofa_diff_fill(ld: c_int, n: c_int, diffdata: *mut c_double, loc: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn ofa_diff_ip(
        glob: c_int,
        loc: c_int,
        coeff: c_double,
        dim: c_int,
        n: c_int,
        ld: c_int,
        dbasis: *mut c_double,
        basis: *mut c_double,
        diffdata: *mut c_double,
        stiff: *mut c_double,
        force: *mut c_double,
        sdet: c_double,
    );
    fn ofa_default_update(ld: c_int, stiff: *mut c_double, force: *mut c_double);
    fn ofa_default_bulk();
    fn ofa_default_finish();
    fn ofa_default_dirichlet();
    fn ofa_default_solve() -> c_double;
    fn ofa_info_alldone();
    fn ofa_fatal_alloc();
}

/// Entry SAVE: the original's `Visited` flag. `OFMesh`/`NoDir` themselves
/// live wrapper-side with identical lifetimes.
struct OfState {
    visited: bool,
}

impl OfState {
    const fn new() -> Self {
        Self { visited: false }
    }
}

static STATE: Mutex<OfState> = Mutex::new(OfState::new());

/// `DiffusionAssembly` SAVE state: `STIFF(N,N)`, `FORCE(N)`,
/// `DataDiffusivity(3,3,N)` plus the `Visited` flag. `N` is the first
/// call's `MaxElementNodes`.
struct DiffState {
    stiff: Vec<c_double>,
    force: Vec<c_double>,
    diffdata: Vec<c_double>,
    n: usize,
    visited: bool,
}

impl DiffState {
    const fn new() -> Self {
        Self {
            stiff: Vec::new(),
            force: Vec::new(),
            diffdata: Vec::new(),
            n: 0,
            visited: false,
        }
    }
}

static DIFF: Mutex<DiffState> = Mutex::new(DiffState::new());

/// Grow the diffusion SAVE arrays when an element needs more than the
/// first-mesh sizing. The original never grows (its comment assumes the
/// first mesh is big enough); growth only triggers where the original
/// would overflow silently, and the leading sections used by the wrappers
/// are unaffected, so valid inputs never observe it.
fn ensure_diff(state: &mut DiffState, need: usize) {
    if need <= state.n {
        return;
    }
    state.stiff.resize(need * need, 0.0);
    state.force.resize(need, 0.0);
    state.diffdata.resize(9 * need, 0.0);
    state.n = need;
}

/// `openfoam2elmerfit_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "openfoam2elmerfit_"]
pub unsafe extern "C" fn openfoam2elmerfit_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { ofa_begin(solver) };
    unsafe { ofa_info_head() };

    // Var => Solver % Variable; Fatal when unassociated.
    unsafe { ofa_var_check() };

    // Amat => Solver % Matrix; RhsVector => Amat % Rhs;
    // ALLOCATE( WeightVector( SIZE( RhsVector ) ) ).
    let nrhs = unsafe { ofa_rhs_len() };
    unsafe { ofa_weight_alloc(nrhs) };

    // IF( Visited ) { ReleaseMesh; DEALLOCATE }; OFMesh => AllocateMesh().
    {
        let mut state = STATE.lock().unwrap();
        if state.visited {
            unsafe { ofa_release_saved() };
        }
        unsafe { ofa_alloc_mesh() };
        // PassiveCoord / CylSymm / MeshDim reduced-dimension warning.
        let mut found = 0;
        unsafe { ofa_passive(&mut found) };
        let _ = found;
        let passive = unsafe { ofa_passive_val() };
        let cyl = unsafe { ofa_cylsymm() };
        let meshdim = unsafe { ofa_meshdim() };
        if meshdim < 3 && passive == 0 && cyl == 0 {
            unsafe { ofa_warn_reduced() };
        }

        // BaseDir = GetString(...); Fatal when missing.
        let mut hasbase = 0;
        unsafe { ofa_basedir(&mut hasbase) };
        let _ = hasbase;

        // IF( .NOT. ListCheckPresent(Params,'OpenFOAM Mesh 1') ) Blocks().
        if unsafe { ofa_have_mesh1() } == 0 {
            unsafe { ofa_blocks() };
        }
        unsafe { ofa_parallel_reduce() };
        let nodir = unsafe { ofa_nodir() };
        unsafe { ofa_info_nblocks() };

        unsafe { ofa_default_initialize() };
        unsafe { ofa_weight_zero() };

        for i in 1..=nodir {
            if nodir > 1 {
                unsafe { ofa_info_block_head(i) };
            }
            // DirName fetch; Fatal when missing.
            let mut dirfound = 0;
            unsafe { ofa_mesh_dir(i, &mut dirfound) };
            let _ = dirfound;
            // FileName = TRIM(DirName)//'C'; Info reading.
            unsafe { ofa_info_reading() };
            unsafe { ofa_create_mesh() };
            unsafe { ofa_field_alloc() };
            let mut got = 0;
            unsafe { ofa_read_field(&mut got) };
            if got != 0 {
                unsafe { ofa_report_stats() };
                data_assembly(cyl);
            }
            // Only one OpenFOAM mesh at a time: release the structures
            // when a second mesh follows.
            if i < nodir {
                unsafe { ofa_release_contents() };
            }
            unsafe { ofa_field_free() };
        }

        // Fit scaling (default 1.0 when the keyword is absent).
        let mut fitfound = 0;
        let mut coeff = unsafe { ofa_fit_coeff(&mut fitfound) };
        if fitfound == 0 {
            coeff = 1.0;
        }
        unsafe { ofa_apply_fit(coeff) };

        diffusion_assembly();

        unsafe {
            ofa_default_bulk();
            ofa_default_finish();
            ofa_default_dirichlet();
        }

        // And finally, solve (Norm is local and unused, as in the original).
        let _norm: c_double = unsafe { ofa_default_solve() };

        state.visited = true;
    }
    unsafe { ofa_weight_free() };

    unsafe { ofa_info_alldone() };
}

/// The `DataAssembly` internal procedure: project the OpenFOAM point data
/// onto the Elmer mesh with basis-function weighting. `cyl` selects the
/// cylindrically-symmetric projection of the global coordinates.
fn data_assembly(cyl: c_int) {
    // N = Mesh % MaxElementNodes; ALLOCATE( Basis(n) ).
    let nn = unsafe { ofa_max_elem_nodes() } as usize;
    let mut basis = vec![0.0 as c_double; nn];
    let mut element_index: c_int = 0;
    let npts = unsafe { ofa_ofcount() };
    for t in 1..=npts {
        let mut g = [0.0 as c_double; 3];
        unsafe { ofa_ofcoords(t, cyl, g.as_mut_ptr()) };
        let val: c_double = unsafe { ofa_field_val(t) };
        // Find the element and local coordinates of the data point.
        let mut l = [0.0 as c_double; 3];
        unsafe { ofa_locate(g.as_mut_ptr(), l.as_mut_ptr(), &mut element_index) };
        if element_index == 0 {
            continue;
        }
        let mut n: c_int = 0;
        unsafe { ofa_data_begin_elem(element_index, &mut n) };
        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            ofa_data_eleminfo(
                l[0],
                l[1],
                l[2],
                n,
                basis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            )
        };
        let _ = (detj, stat);
        // Share the data to nodes using the basis functions for weighting.
        unsafe { ofa_data_scatter(val, n, basis.as_mut_ptr()) };
    }
}

/// The `DiffusionAssembly` internal procedure: assemble the regularization
/// diffusion matrix. First-mesh SAVE sizing mirrors the original; the
/// `MATMUL`/`SUM` update statements run in the Fortran wrapper.
fn diffusion_assembly() {
    let mut glob: c_int = 0;
    let coeff: c_double = unsafe { ofa_diff_coeff(&mut glob) };
    let active = unsafe { ofa_nof_active() };
    let dim = unsafe { ofa_coord_dim() };
    let mut ds = DIFF.lock().unwrap();
    if !ds.visited {
        let n = unsafe { ofa_max_elem_nodes() } as usize;
        ds.n = n;
        if ds.stiff.try_reserve(n * n).is_err()
            || ds.force.try_reserve(n).is_err()
            || ds.diffdata.try_reserve(9 * n).is_err()
        {
            unsafe { ofa_fatal_alloc() };
            std::process::abort();
        }
        ds.stiff.resize(n * n, 0.0);
        ds.force.resize(n, 0.0);
        ds.diffdata.resize(9 * n, 0.0);
        ds.visited = true;
    }
    for t in 1..=active {
        let mut n: c_int = 0;
        unsafe { ofa_begin_active(t, &mut n) };
        // First-mesh sizing holds on valid inputs; grow-only guard covers
        // remeshed re-entry (invisible otherwise).
        ensure_diff(&mut ds, n as usize);
        let ld = ds.n as c_int;
        ds.stiff.fill(0.0);
        ds.force.fill(0.0);
        let mut loc: c_int = 0;
        if glob == 0 {
            unsafe { ofa_diff_fill(ld, n, ds.diffdata.as_mut_ptr(), &mut loc) };
        }
        // Basis(n) / dBasisdx(n,3) automatics, sized by this element.
        let nn = n as usize;
        let mut basis = vec![0.0 as c_double; nn];
        let mut dbasis = vec![0.0 as c_double; nn * 3];
        unsafe { ofa_diff_nodes() };
        let nip = unsafe { ofa_gauss_begin() };
        for k in 1..=nip {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { ofa_gauss_point(k, &mut u, &mut v, &mut w, &mut s) };
            let (mut detj, mut stat) = (0.0, 0);
            unsafe {
                ofa_diff_eleminfo(
                    u,
                    v,
                    w,
                    n,
                    basis.as_mut_ptr(),
                    dbasis.as_mut_ptr(),
                    &mut detj,
                    &mut stat,
                )
            };
            let _ = stat;
            unsafe {
                ofa_diff_ip(
                    glob,
                    loc,
                    coeff,
                    dim,
                    n,
                    ld,
                    dbasis.as_mut_ptr(),
                    basis.as_mut_ptr(),
                    ds.diffdata.as_mut_ptr(),
                    ds.stiff.as_mut_ptr(),
                    ds.force.as_mut_ptr(),
                    s * detj,
                )
            };
        }
        unsafe { ofa_default_update(ld, ds.stiff.as_mut_ptr(), ds.force.as_mut_ptr()) };
    }
}
