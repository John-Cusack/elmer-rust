//! Port of `fem/src/modules/HelmholtzBEM.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/HelmholtzBEM.F90` — one F77-global entry
//! `helmholtzbemsolver_` plus the two F77-global helpers `matvec_` and
//! `precond_` (Krylov callbacks passed by address to `IterCall`), plus the
//! plugin-private BEM matrix module `GlobMatC`. The emitted dynamic
//! symbols are exactly `helmholtzbemsolver_`, `matvec_`, `precond_` (text)
//! and `__globmatc_MOD_matrix` (bss), proven by `nm -D` on the reference
//! `.so`; the internals emit only local symbols.
//!
//! Triage note: `matvec`/`precond` carry only plain-ABI assumed-size
//! dummies (`COMPLEX(KIND=dp) :: u(*), v(*)` + `INTEGER :: ipar(*)`), no
//! `CHARACTER`, no descriptors. The `HUTI_*` names are C-preprocessor
//! macros over the `ipar`/`dpar` slots (`huti_fdefs.h`: `HUTI_NDIM` is
//! `ipar(3)`, `HUTI_WRKDIM` `ipar(4)`, `HUTI_INFO` `ipar(30)`, …), so the
//! whole HUTI sequence transcribes with plain integer indexing and no new
//! ABI shape. `IterCallFTNC` (descriptor dummies) is reached through a
//! Fortran wrapper that builds the descriptors at the call site, exactly
//! as the original does.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/helmholtzbem_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. All
//! floating-point work (the BEM assembly kernels, `MATMUL`, complex
//! division, `SQRT`/`SUM` reductions, the `CPUTime`-based reports) stays
//! in Fortran wrappers, so bit-identity holds by construction wherever
//! the call sequence matches.
//!
//! Faithfulness notes (mirrors HelmholtzBEM.F90 statement for statement):
//! - `SAVE` semantics replicated: first-visit sizing of the permutation
//!   maps and work arrays (`Model % MaxElementNodes` / mesh node counts),
//!   never resized; the wrapper-side SAVE mirrors the original's
//!   internal-procedure SAVE with identical lifetimes (element context
//!   and `ElementNodes` consumed before the next fill, as in the
//!   original). The `GlobMatC` matrix is allocated once at first visit.
//! - The author's evident sizing intent is used for `P1`/`P2`
//!   (`Model % MaxElementNodes`; the original sizes them with a stale
//!   local whose value never leaks into any observable). Allocation
//!   failure calls the same `Fatal('HelmholtzBEMSolver', ...)` trio via
//!   `try_reserve` (the `STAT=istat` check), never aborting directly.
//! - `Model`, `dt`, `TransientSimulation` are unused by the body except
//!   for the `Model` component traffic, exactly as in the original
//!   (received and ignored).
//! - `EquationName` is fetched and discarded, as in the original.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.
//!   In particular `matvec_`/`precond_` take no lock (the entry holds it
//!   while HUTI re-enters them).
use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn hbm_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn hbm_check_mpi();
    fn hbm_check_cart();
    fn hbm_sizes(nmaxnodes: *mut c_int, numnodes: *mut c_int, nbcs: *mut c_int);
    fn hbm_nof_active() -> c_int;
    fn hbm_begin_active(
        t: c_int,
        code: *mut c_int,
        n: *mut c_int,
        assoc: *mut c_int,
        constraint: *mut c_int,
    );
    fn hbm_elem_nodes(n: c_int, buf: *mut c_int);
    fn hbm_alloc_matrix(nb: c_int);
    fn hbm_alloc_nodes(n: c_int);
    fn hbm_matrix_zero();
    fn hbm_var_values_ptr() -> *mut c_double;
    fn hbm_var_perm_ptr() -> *mut c_int;
    fn hbm_touch_equation();
    fn hbm_angular_frequency() -> c_double;
    fn hbm_flux_init() -> c_int;
    fn hbm_flux_values_ptr() -> *mut c_double;
    fn hbm_bc_tag(i: c_int) -> c_int;
    fn hbm_bc_select(i: c_int);
    fn hbm_bc_pot(comp: c_int, n: c_int, p: *mut c_double, gotit: *mut c_int);
    fn hbm_bc_flux(comp: c_int, n: c_int, p: *mut c_double, gotit: *mut c_int);
    fn hbm_elemnodes_fill(n: c_int, nidx: *mut c_int);
    fn hbm_precond_elem(i: c_int, u: *mut Cx, v: *mut Cx);
    fn hbm_matvec(n: c_int, u: *mut Cx, v: *mut Cx);
    #[allow(clippy::too_many_arguments)]
    fn hbm_integrate(
        n: c_int,
        nb: c_int,
        freq: c_double,
        potperm: *mut c_int,
        bndperm: *mut c_int,
        load: *mut Cx,
        force: *mut Cx,
        diag: *mut Cx,
        potknown: *mut c_int,
        numnodes: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn hbm_compute_pot(
        n: c_int,
        nb: c_int,
        freq: c_double,
        potential: *mut Cx,
        pot: *mut Cx,
        flx: *mut Cx,
        potperm: *mut c_int,
        numnodes: c_int,
    );
    fn hbm_post_adjust(
        nb: c_int,
        force: *mut Cx,
        load: *mut Cx,
        diag: *mut Cx,
        potknown: *mut c_int,
    );
    fn hbm_var_norm(pot: *mut Cx, n: c_int, numnodes: c_int) -> c_double;
    fn hbm_lapack(n: c_int, x: *mut Cx);
    #[allow(clippy::too_many_arguments)]
    fn hbm_itercall(
        fptr: *const i64,
        x: *mut Cx,
        b: *mut Cx,
        n: c_int,
        ipar: *mut c_int,
        dpar: *mut c_double,
        work: *mut Cx,
        wsize: c_int,
        mvptr: *const i64,
        pcptr: *const i64,
    );
    fn hbm_linsys() -> c_int;
    fn hbm_tol(val: *mut c_double, gotit: *mut c_int);
    fn hbm_divlim(val: *mut c_double, gotit: *mut c_int);
    fn hbm_maxit(val: *mut c_int, gotit: *mut c_int);
    fn hbm_resid(val: *mut c_int, gotit: *mut c_int);
    fn hbm_abort(val: *mut c_int, gotit: *mut c_int);
    fn hbm_info_solving();
    fn hbm_info_nbn(nb: c_int);
    fn hbm_info_full();
    fn hbm_info_extract();
    fn hbm_info_compute();
    fn hbm_info_fields();
    fn hbm_info_done();
    fn hbm_fatal_alloc1();
    fn hbm_fatal_alloc2();
    fn hbm_fatal_alloc3();
    fn hbm_fatal_linsys();
    fn hbm_fatal_iter();
    fn hbm_error_iter();
    fn hbm_invalidate();
    fn hbm_cputime() -> c_double;
    fn hbm_print_assembly(t: c_double);
    fn hbm_print_solve(t: c_double);
    fn hbm_print_post(t: c_double);
    // Host HUTI entry reached by address, exactly as the original
    // reaches it through `AddrFunc(HUTI_Z_BICGSTAB)`.
    fn __huti_sfe_MOD_huti_z_bicgstab();
}

/// Complex double matching `COMPLEX(KIND=dp)` layout (C-interoperable).
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct Cx {
    pub re: c_double,
    pub im: c_double,
}

/// First-visit SAVE state: permutation maps, work arrays, flags.
struct HelmholtzBemState {
    allocations_done: bool,
    boundary_nodes: usize,
    numnodes: usize,
    nmaxnodes: usize,
    potential_perm: Vec<c_int>,
    boundary_perm: Vec<c_int>,
    potential_known: Vec<c_int>,
    pot: Vec<Cx>,
    flx: Vec<Cx>,
    load: Vec<Cx>,
    potential: Vec<Cx>,
    force_vector: Vec<Cx>,
    diagonal: Vec<Cx>,
    p1: Vec<c_double>,
    p2: Vec<c_double>,
}

impl HelmholtzBemState {
    const fn new() -> Self {
        Self {
            allocations_done: false,
            boundary_nodes: 0,
            numnodes: 0,
            nmaxnodes: 0,
            potential_perm: Vec::new(),
            boundary_perm: Vec::new(),
            potential_known: Vec::new(),
            pot: Vec::new(),
            flx: Vec::new(),
            load: Vec::new(),
            potential: Vec::new(),
            force_vector: Vec::new(),
            diagonal: Vec::new(),
            p1: Vec::new(),
            p2: Vec::new(),
        }
    }
}

static STATE: Mutex<HelmholtzBemState> = Mutex::new(HelmholtzBemState::new());

/// `helmholtzbemsolver_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "helmholtzbemsolver_"]
pub unsafe extern "C" fn helmholtz_bem_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { hbm_begin(model, solver) };

    unsafe { hbm_check_mpi() };
    unsafe { hbm_check_cart() };
    unsafe { hbm_info_solving() };

    let (mut nmaxnodes, mut numnodes, mut nbcs) = (0, 0, 0);
    unsafe { hbm_sizes(&mut nmaxnodes, &mut numnodes, &mut nbcs) };
    let (nmaxnodes, numnodes, nbcs) =
        (nmaxnodes as usize, numnodes as usize, nbcs as usize);

    // First-visit permanent storage (bug-compatible: sized once, never
    // resized; the P1/P2 sizing follows the author's evident intent).
    if !state.allocations_done {
        if state.potential_perm.try_reserve(numnodes).is_err()
            || state.boundary_perm.try_reserve(numnodes).is_err()
        {
            unsafe { hbm_fatal_alloc1() };
            std::process::abort();
        }
        state.potential_perm.resize(numnodes, 0);
        state.boundary_perm.resize(numnodes, 0);

        // Boundary-node numbering pass.
        let mut nb: usize = 0;
        let nof_active = unsafe { hbm_nof_active() };
        let mut nidx = vec![0 as c_int; nmaxnodes.max(1)];
        for t in 1..=nof_active {
            let (mut code, mut n, mut assoc, mut constraint) = (0, 0, 0, 0);
            unsafe { hbm_begin_active(t, &mut code, &mut n, &mut assoc, &mut constraint) };
            let _ = constraint;
            if assoc == 0 {
                continue;
            }
            if code == 101 {
                continue;
            }
            unsafe { hbm_elem_nodes(n, nidx.as_mut_ptr()) };
            for &k in nidx[..n as usize].iter() {
                let k = k as usize;
                if state.potential_perm[k - 1] == 0 {
                    nb += 1;
                    state.boundary_perm[nb - 1] = k as c_int;
                    state.potential_perm[k - 1] = nb as c_int;
                }
            }
        }
        state.boundary_nodes = nb;
        unsafe { hbm_info_nbn(nb as c_int) };

        state.nmaxnodes = nmaxnodes;
        state.numnodes = numnodes;
        unsafe { hbm_alloc_matrix(nb as c_int) };
        unsafe { hbm_alloc_nodes(nmaxnodes as c_int) };
        if state.flx.try_reserve(nb).is_err()
            || state.pot.try_reserve(nb).is_err()
            || state.load.try_reserve(nb).is_err()
            || state.diagonal.try_reserve(nb).is_err()
            || state.potential_known.try_reserve(nb).is_err()
        {
            unsafe { hbm_fatal_alloc2() };
            std::process::abort();
        }
        state.flx.resize(nb, Cx::default());
        state.pot.resize(nb, Cx::default());
        state.load.resize(nb, Cx::default());
        state.diagonal.resize(nb, Cx::default());
        state.potential_known.resize(nb, 0);
        if state.potential.try_reserve(numnodes).is_err()
            || state.force_vector.try_reserve(numnodes).is_err()
            || state.p1.try_reserve(nmaxnodes).is_err()
            || state.p2.try_reserve(nmaxnodes).is_err()
        {
            unsafe { hbm_fatal_alloc3() };
            std::process::abort();
        }
        state.potential.resize(numnodes, Cx::default());
        state.force_vector.resize(numnodes, Cx::default());
        state.p1.resize(nmaxnodes, 0.0);
        state.p2.resize(nmaxnodes, 0.0);
        state.allocations_done = true;
    }
    let nb = state.boundary_nodes;

    // Angular frequency; wall-clock start; Equation touch (discarded).
    let freq = unsafe { hbm_angular_frequency() };
    let at = unsafe { hbm_cputime() };
    unsafe { hbm_touch_equation() };

    unsafe { hbm_matrix_zero() };
    state.load.fill(Cx::default());
    state.diagonal.fill(Cx::default());
    state.force_vector.fill(Cx::default());

    // Potential from the solver variable (interleaved real/imag pairs).
    let values = unsafe { hbm_var_values_ptr() };
    for (i, slot) in state.potential.iter_mut().enumerate() {
        *slot = Cx {
            re: unsafe { *values.add(2 * i) },
            im: unsafe { *values.add(2 * i + 1) },
        };
    }

    // Boundary conditions: for each node either flux or potential is
    // given, the other is the unknown.
    let nof_active = unsafe { hbm_nof_active() };
    let mut nidx = vec![0 as c_int; state.nmaxnodes.max(1)];
    for t in 1..=nof_active {
        let (mut code, mut n, mut assoc, mut constraint) = (0, 0, 0, 0);
        unsafe { hbm_begin_active(t, &mut code, &mut n, &mut assoc, &mut constraint) };
        if assoc == 0 {
            continue;
        }
        if code == 101 {
            continue;
        }
        unsafe { hbm_elem_nodes(n, nidx.as_mut_ptr()) };
        let nn = n as usize;
        for i in 1..=nbcs as c_int {
            if constraint != unsafe { hbm_bc_tag(i) } {
                continue;
            }
            unsafe { hbm_bc_select(i) };
            let mut gotit = 0;
            unsafe {
                hbm_bc_pot(1, n, state.p1.as_mut_ptr(), &mut gotit);
                hbm_bc_pot(2, n, state.p2.as_mut_ptr(), &mut gotit);
            }
            for j in 0..nn {
                let q = state.potential_perm[nidx[j] as usize - 1] as usize;
                state.load[q - 1] = Cx {
                    re: state.p1[j],
                    im: state.p2[j],
                };
            }
            if gotit == 0 {
                for j in 0..nn {
                    let q = state.potential_perm[nidx[j] as usize - 1] as usize;
                    state.potential_known[q - 1] = 0;
                }
                unsafe {
                    hbm_bc_flux(1, n, state.p1.as_mut_ptr(), &mut gotit);
                    hbm_bc_flux(2, n, state.p2.as_mut_ptr(), &mut gotit);
                }
                let _ = gotit;
                for j in 0..nn {
                    let q = state.potential_perm[nidx[j] as usize - 1] as usize;
                    state.load[q - 1] = Cx {
                        re: state.p1[j],
                        im: state.p2[j],
                    };
                }
            } else {
                for j in 0..nn {
                    let q = state.potential_perm[nidx[j] as usize - 1] as usize;
                    state.potential_known[q - 1] = 1;
                }
            }
            break;
        }
    }

    // Matrix assembly.
    for t in 1..=nof_active {
        let (mut code, mut n, mut assoc, mut constraint) = (0, 0, 0, 0);
        unsafe { hbm_begin_active(t, &mut code, &mut n, &mut assoc, &mut constraint) };
        if assoc == 0 {
            continue;
        }
        if code == 101 {
            continue;
        }
        unsafe { hbm_elem_nodes(n, nidx.as_mut_ptr()) };
        unsafe { hbm_elemnodes_fill(n, nidx.as_mut_ptr()) };
        unsafe {
            hbm_integrate(
                n,
                nb as c_int,
                freq,
                state.potential_perm.as_mut_ptr(),
                state.boundary_perm.as_mut_ptr(),
                state.load.as_mut_ptr(),
                state.force_vector.as_mut_ptr(),
                state.diagonal.as_mut_ptr(),
                state.potential_known.as_mut_ptr(),
                numnodes as c_int,
            );
        }
    }

    unsafe {
        hbm_post_adjust(
            nb as c_int,
            state.force_vector.as_mut_ptr(),
            state.load.as_mut_ptr(),
            state.diagonal.as_mut_ptr(),
            state.potential_known.as_mut_ptr(),
        );
    }

    let now = unsafe { hbm_cputime() };
    unsafe { hbm_print_assembly(now - at) };

    // Solve the system.
    let st = unsafe { hbm_cputime() };
    unsafe { hbm_info_full() };
    solve_full(&mut state, nb);
    unsafe { hbm_info_extract() };
    // Indexed mirror of the original whole-array DO loop.
    #[allow(clippy::needless_range_loop)]
    for i in 0..nb {
        if state.potential_known[i] != 0 {
            state.flx[i] = state.potential[i];
            state.pot[i] = state.load[i];
        } else {
            state.flx[i] = state.load[i];
            state.pot[i] = state.potential[i];
        }
    }
    let now = unsafe { hbm_cputime() };
    unsafe { hbm_print_solve(now - st) };

    // Interior potential from the boundary representation.
    let st = unsafe { hbm_cputime() };
    unsafe { hbm_info_compute() };
    state.potential.fill(Cx::default());
    // Indexed mirror of the original whole-array DO loop.
    #[allow(clippy::needless_range_loop)]
    for i in 0..nb {
        let k = state.boundary_perm[i] as usize;
        state.potential[k - 1] = state.pot[i];
    }
    for t in 1..=nof_active {
        let (mut code, mut n, mut assoc, mut constraint) = (0, 0, 0, 0);
        unsafe { hbm_begin_active(t, &mut code, &mut n, &mut assoc, &mut constraint) };
        if assoc == 0 {
            continue;
        }
        if code == 101 {
            continue;
        }
        unsafe { hbm_elem_nodes(n, nidx.as_mut_ptr()) };
        unsafe { hbm_elemnodes_fill(n, nidx.as_mut_ptr()) };
        unsafe {
            hbm_compute_pot(
                n,
                nb as c_int,
                freq,
                state.potential.as_mut_ptr(),
                state.pot.as_mut_ptr(),
                state.flx.as_mut_ptr(),
                state.potential_perm.as_mut_ptr(),
                numnodes as c_int,
            );
        }
    }

    // Postprocessing fields.
    unsafe { hbm_info_fields() };
    let values = unsafe { hbm_var_values_ptr() };
    unsafe { std::slice::from_raw_parts_mut(values, 2 * numnodes).fill(0.0) };
    let perm = unsafe { hbm_var_perm_ptr() };
    // Indexed mirror of the original whole-array DO loop.
    #[allow(clippy::needless_range_loop)]
    for i in 0..numnodes {
        let j = unsafe { *perm.add(i) };
        if j > 0 {
            unsafe {
                *values.add(2 * (j as usize - 1)) = state.potential[i].re;
                *values.add(2 * (j as usize - 1) + 1) = state.potential[i].im;
            }
        }
    }
    if unsafe { hbm_flux_init() } != 0 {
        let fvalues = unsafe { hbm_flux_values_ptr() };
        unsafe { std::slice::from_raw_parts_mut(fvalues, 2 * numnodes).fill(0.0) };
        // Indexed mirror of the original whole-array DO loop.
        #[allow(clippy::needless_range_loop)]
        for i in 0..nb {
            let k = state.boundary_perm[i] as usize;
            unsafe {
                *fvalues.add(2 * (k - 1)) = state.flx[i].re;
                *fvalues.add(2 * (k - 1) + 1) = state.flx[i].im;
            }
        }
    }

    unsafe { hbm_var_norm(state.potential.as_mut_ptr(), numnodes as c_int, numnodes as c_int) };
    unsafe { hbm_invalidate() };
    let now = unsafe { hbm_cputime() };
    unsafe { hbm_print_post(now - st) };
    unsafe { hbm_info_done() };
}

/// The `SolveFull` internal procedure: direct LAPACK or HUTI iterative
/// solve dispatch on `Linear System Solver`.
fn solve_full(state: &mut HelmholtzBemState, nb: usize) {
    match unsafe { hbm_linsys() } {
        1 => {
            unsafe { hbm_lapack(nb as c_int, state.force_vector.as_mut_ptr()) };
            state.potential[..nb].copy_from_slice(&state.force_vector[..nb]);
        }
        2 => full_iter_solver(state, nb),
        _ => {
            unsafe { hbm_fatal_linsys() };
            std::process::abort();
        }
    }
}

/// The `FullIterSolver` internal procedure: HUTI BiCGSTAB sequence.
/// `HUTI_*` slots are preprocessor macros over `ipar`/`dpar`
/// (`huti_fdefs.h`): NDIM=ipar(3), WRKDIM=ipar(4), DBUGLVL=ipar(5),
/// MAXIT=ipar(10), INITIALX=ipar(14), INFO=ipar(30),
/// TOLERANCE=dpar(1), MAXTOLERANCE=dpar(2); RANDOMX=0,
/// USERSUPPLIEDX=1, CONVERGENCE=1, BICGSTAB_WORKSIZE=8.
fn full_iter_solver(state: &mut HelmholtzBemState, nb: usize) {
    let n = nb as c_int;
    let mut ipar = [0 as c_int; 50];
    let mut dpar = [0.0 as c_double; 50];
    ipar[2] = n;
    ipar[3] = 8;
    let wsize = 8 as c_int;
    let mut work = vec![Cx::default(); (wsize as usize) * nb];
    let allzero = state
        .potential
        .iter()
        .all(|c| c.re == 0.0 && c.im == 0.0);
    ipar[13] = if allzero { 0 } else { 1 };
    let (mut tol, mut gotit) = (0.0, 0);
    unsafe { hbm_tol(&mut tol, &mut gotit) };
    dpar[0] = tol;
    let (mut maxtol, mut gotit) = (0.0, 0);
    unsafe { hbm_divlim(&mut maxtol, &mut gotit) };
    dpar[1] = if gotit == 0 { 1.0e20 } else { maxtol };
    let (mut maxit, mut gotit) = (0, 0);
    unsafe { hbm_maxit(&mut maxit, &mut gotit) };
    let _ = gotit;
    ipar[9] = maxit;
    let (mut dbg, mut gotit) = (0, 0);
    unsafe { hbm_resid(&mut dbg, &mut gotit) };
    ipar[4] = if gotit == 0 { 1 } else { dbg };
    let (mut abortval, mut gotit) = (0, 0);
    unsafe { hbm_abort(&mut abortval, &mut gotit) };
    let abort = if gotit == 0 { true } else { abortval != 0 };

    let iterproc = __huti_sfe_MOD_huti_z_bicgstab as *const () as usize as i64;
    let mvproc = matvec_entry as *const () as usize as i64;
    let pcproc = precond_entry as *const () as usize as i64;
    unsafe {
        hbm_itercall(
            &iterproc,
            state.potential.as_mut_ptr(),
            state.force_vector.as_mut_ptr(),
            n,
            ipar.as_mut_ptr(),
            dpar.as_mut_ptr(),
            work.as_mut_ptr(),
            wsize,
            &mvproc,
            &pcproc,
        );
    }
    if ipar[29] != 1 {
        if abort {
            unsafe { hbm_fatal_iter() };
            std::process::abort();
        } else {
            unsafe { hbm_error_iter() };
        }
    }
}

/// `matvec_` — the HUTI matrix-vector callback
/// (`HUTI_NDIM` is `ipar(3)`).
///
/// # Safety
///
/// `u`/`v` must point at `ipar(3)` readable/writable complex values.
#[export_name = "matvec_"]
pub unsafe extern "C" fn matvec_entry(u: *mut Cx, v: *mut Cx, ipar: *mut c_int) {
    let n = unsafe { *ipar.add(2) };
    unsafe { hbm_matvec(n, u, v) };
}

/// `precond_` — the HUTI diagonal-preconditioner callback
/// (`HUTI_NDIM` is `ipar(3)`).
///
/// # Safety
///
/// Same contract as [`matvec_entry`].
#[export_name = "precond_"]
pub unsafe extern "C" fn precond_entry(u: *mut Cx, v: *mut Cx, ipar: *mut c_int) {
    let n = unsafe { *ipar.add(2) };
    for i in 1..=n {
        unsafe { hbm_precond_elem(i, u, v) };
    }
}
