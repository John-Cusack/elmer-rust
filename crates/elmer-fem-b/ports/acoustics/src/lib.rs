//! Port of `fem/src/modules/Acoustics.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/Acoustics.F90` — two F77-global entries
//! (`acousticssolver_`, `acousticssolver_init_`) with ~60 `CONTAINS`
//! internals. The emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`); every other routine is a
//! hidden symbol (enforced by `acoustics.map`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/acoustics_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset.
//!
//! Transcription split (documented where the call sequence matches):
//! - In Rust: `AcousticsSolver_init`, the full `AcousticsSolver` driver
//!   control flow, and the proof-exercised assembly kernels `LocalMatrix`,
//!   `LocalMatrixBoundary`, `LCondensate`, `LCondensateBubbles`.
//! - In Fortran (hidden, same statements, same compiler, hence
//!   bit-identical): the BlockPreconditioning-only iterative suite and
//!   preconditioner builders, the slip/interface-only kernels, the
//!   fluidic-force / acoustic-impedance / BEM-coupling bodies, and the
//!   shell-FSI helper. Complex multiply/divide (`__muldc3`/`__divdc3`
//!   lowering the Rust side must not second-guess), order-sensitive
//!   reductions (`SUM`/`MATMUL`), formatted/list-directed I/O, and
//!   whole-array updates stay in Fortran wrappers, so bit-identity holds
//!   by construction wherever the call sequence matches. Complex add/sub,
//!   `CMPLX` construction, `REAL`/`AIMAG` and complex-by-real scaling are
//!   componentwise and reproduced inline.
//! - `SAVE` semantics replicated explicitly: all driver SAVE arrays live
//!   in the `Mutex`-guarded state below with first-mesh sizing and
//!   realloc-on-`MeshChanged`; wrapper-side SAVE mirrors the
//!   internal-procedure SAVE with identical lifetimes (BEM arrays, block
//!   matrices, `Nodes`/`IP` contexts); no Rust threads.
//! - Allocation failure calls the same `Fatal('AcousticsSolver', ...)`
//!   via `try_reserve`, instead of aborting.
//! - `IF(.FALSE.)` arms (`NodesOnBoundary` census, stabilisation test,
//!   reliable-update block, scanning/`ZAppr` leftovers) are omitted as
//!   unreachable in the original; `IF(.TRUE.) AcousticShellInterface()`
//!   is unconditional. `MyBubbles` is identically `.FALSE.`.
//! - `DIM`/`dim` and `CoordSys`/`coordsys` are the same Fortran
//!   identifiers (case-insensitive); the port follows suit.
//! - Host association: the original internals read the driver's `i`,
//!   `stat`, `n` and `GotIt` via host association. Those uses are
//!   provably unobservable as aliases (synchronous calls; every driver
//!   read follows a reassignment), so the port gives the kept Fortran
//!   routines explicit locals and the Rust kernels Rust locals.
//! - `RelOrder` is `0` when unset (`ListGetInteger` default); `Bulk
//!   Viscosity` keeps its leading-space keyword `' Bulk Viscosity'`.
//! - Cartesian = 1, AxisSymmetric = 4 (`CoordinateSystems` params).

use std::ffi::c_void;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn aca_init_begin(h_solver: *const c_void);
    fn aca_info_init();
    fn aca_init_add_complex();
    fn aca_init_add_variable(is2d: c_int);
    fn aca_coord_dim() -> c_int;
    fn aca_begin(h_model: *mut c_void, h_solver: *mut c_void);
    fn aca_assoc_matrix() -> c_int;
    fn aca_flow_vals() -> *mut c_double;
    fn aca_flow_perm() -> *mut c_int;
    fn aca_flow_dofs() -> c_int;
    fn aca_var_len(nv: *mut c_int, np: *mut c_int);
    fn aca_var_norm() -> c_double;
    fn aca_var_norm_set(x: c_double);
    fn aca_local_nodes() -> c_int;
    fn aca_model_num_bcs() -> c_int;
    fn aca_model_num_nodes() -> c_int;
    fn aca_mesh_changed() -> c_int;
    fn aca_max_elem_dofs() -> c_int;
    fn aca_num_active() -> c_int;
    fn aca_num_bulk() -> c_int;
    fn aca_num_bnd_elems() -> c_int;
    fn aca_sim_type(s: *mut c_char, n: c_int, got: *mut c_int);
    fn aca_time_var(found: *mut c_int, v1: *mut c_double);
    fn aca_fatal_timevar();
    fn aca_sim_logical(which: c_int, val: *mut c_int, got: *mut c_int);
    fn aca_find_pres(found: *mut c_int);
    fn aca_info_no_pres();
    fn aca_dd_pres();
    fn aca_h_at(k: c_int, re: *mut c_double, im: *mut c_double);
    fn aca_hp_at(node: c_int) -> c_int;
    fn aca_ain_size() -> c_int;
    fn aca_ain_at(k: c_int, node: *mut c_int, r1: *mut c_double, r2: *mut c_double);
    fn aca_solver_logical(which: c_int, val: *mut c_int, got: *mut c_int);
    fn aca_solver_string(which: c_int, s: *mut c_char, n: c_int, got: *mut c_int);
    fn aca_set_eqname();
    fn aca_solver_creal(which: c_int, v: *mut c_double, got: *mut c_int);
    fn aca_solver_int(which: c_int, v: *mut c_int, got: *mut c_int);
    fn aca_angular_freq(af: *mut c_double, got_af: *mut c_int, got_f: *mut c_int);
    fn aca_check_coordsys();
    fn aca_cur_coordsys() -> c_int;
    fn aca_relorder_diag(relorder: c_int);
    fn aca_begin_bulk(t: c_int, n: *mut c_int, nd: *mut c_int, nb: *mut c_int, bodyid: *mut c_int);
    fn aca_begin_direct(t: c_int, n: *mut c_int, bodyid: *mut c_int);
    fn aca_elem_nodes();
    fn aca_alloc_nodes(nn: c_int);
    fn aca_nidx_val(j: c_int) -> c_int;
    fn aca_perm_has_zero(n: c_int) -> c_int;
    fn aca_begin_bnd(t: c_int, n: *mut c_int, active: *mut c_int, constraint: *mut c_int);
    fn aca_bc_tag(i: c_int) -> c_int;
    fn aca_use_bc(i: c_int);
    fn aca_check_eqn() -> c_int;
    fn aca_parent_left(ok: *mut c_int);
    fn aca_parent_right(ok: *mut c_int);
    fn aca_fatal_no_parent();
    fn aca_body_material() -> c_int;
    fn aca_parent_material() -> c_int;
    fn aca_use_material(k: c_int);
    fn aca_mat_real(which: c_int, n: c_int, out: *mut c_double, got: *mut c_int);
    fn aca_body_bfcount(k: *mut c_int, got: *mut c_int);
    fn aca_use_bodyforce(k: c_int);
    fn aca_bf_real(which: c_int, n: c_int, out: *mut c_double, got: *mut c_int);
    fn aca_bc_logical(which: c_int, val: *mut c_int, got: *mut c_int);
    fn aca_bc_real(which: c_int, n: c_int, out: *mut c_double, got: *mut c_int);
    fn aca_bc_constreal(which: c_int, v: *mut c_double, got: *mut c_int);
    fn aca_gauss_begin(which: c_int, relorder: c_int) -> c_int;
    fn aca_gauss_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn aca_eleminfo(
        which: c_int,
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        nd: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        ddbasis: *mut c_double,
        detj: *mut c_double,
        bubbles: c_int,
        stat: *mut c_int,
    );
    fn aca_normal(u: c_double, v: c_double, nrm: *mut c_double);
    fn aca_dot(n: c_int, a: *const c_double, b: *const c_double) -> c_double;
    fn aca_nodes_dotx(n: c_int, basis: *const c_double) -> c_double;
    fn aca_sumv(n: c_int, a: *const c_double) -> c_double;
    fn aca_cmul(ar: c_double, ai: c_double, br: c_double, bi: c_double, rr: *mut c_double, ri: *mut c_double);
    fn aca_cdiv(ar: c_double, ai: c_double, br: c_double, bi: c_double, rr: *mut c_double, ri: *mut c_double);
    fn aca_cmatmul(m: c_int, k: c_int, n: c_int, a: *const c_double, b: *const c_double, c: *mut c_double);
    fn aca_cinv(n: c_int, k: *mut c_double);
    fn aca_update_global(lsm: *const c_double, lf: *const c_double, nn: c_int, n: c_int, dofs: c_int);
    fn aca_setdir_velo(comp: c_int, dim: c_int);
    fn aca_setdir_temp(which: c_int, dofs: c_int);
    fn aca_zerorow(row: c_int);
    fn aca_setmat(row: c_int, col: c_int, val: c_double);
    fn aca_rhs_set(row: c_int, val: c_double);
    fn aca_default_start();
    fn aca_default_initialize();
    fn aca_default_finish_bulk();
    fn aca_default_dirichlet();
    fn aca_default_solve() -> c_double;
    fn aca_cputime() -> c_double;
    fn aca_realtime() -> c_double;
    fn aca_info_blank();
    fn aca_info_dashes();
    fn aca_info_starting();
    fn aca_info_done();
    fn aca_info_exiting();
    fn aca_info_freq(val: c_double);
    fn aca_info_result(norm: c_double);
    fn aca_info_relchange(rc: c_double);
    fn aca_info_asmtime(at: c_double);
    fn aca_info_soltime(st: c_double);
    fn aca_info_asm_percent(pct: c_int);
    fn aca_info_asm_est(est: c_double);
    fn aca_info_block();
    fn aca_warn_dofs();
    fn aca_fatal_alloc();
    fn aca_fatal_bubbles_global();
    fn aca_fatal_helmholtz_node();
    fn aca_fatal_helmholtz_bnd();
    fn aca_error_minibubbles(code: c_int);
    fn aca_list_add_freq(val: c_double);
    fn aca_pre_rescale(
        hr: c_double, sh: c_double, dens: c_double, temp: c_double, visc: c_double, lam: c_double,
        af: c_double, f1: c_double, f2: c_double, f3: c_double, f4: c_double,
        g1: *mut c_double, g2: *mut c_double, g3: *mut c_double, g4: *mut c_double,
    );
    fn aca_post_rescale(
        dens: c_double, temp: c_double, lam: c_double, press: c_double, af: c_double,
        f1: c_double, f2: c_double, f3: c_double, f4: c_double,
        g1: *mut c_double, g2: *mut c_double, g3: *mut c_double, g4: *mut c_double,
    );
    fn aca_bem_pre(firstvisit: c_int);
    fn aca_bem_nodes();
    fn aca_bem_write(
        firstvisit: c_int, angularfreq: c_double, cdim: c_int,
        wallvelo: *mut c_double, density: *mut c_double, lod: *mut c_double, wsz: c_int,
    );
    fn aca_bem_post(relchange: c_double, calls: c_int);
    fn aca_fluidic_all(
        cdim: c_int, visc: *mut c_double, lam: *mut c_double, bulkv: *mut c_double,
        wallvelo: *mut c_double, lod: *mut c_double, totforce: *mut c_double,
        totmoment: *mut c_double, totarea: *mut c_double, momout: *mut c_int,
        nn: c_int, nbc: c_int,
    );
    fn aca_impedance_pass(cdim: c_int);
    fn aca_block_create();
    fn aca_block_zero(velasm: c_int);
    #[allow(clippy::too_many_arguments)]
    fn aca_block_bulk(
        visc: *const c_double, dens: *const c_double, sh: *const c_double, hr: *const c_double,
        temp: *const c_double, cond: *const c_double, lam: *const c_double, press: *const c_double,
        angfreq: c_double, nn: c_int, n: c_int, cdim: c_int, coordsys: c_int, velasm: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn aca_block_bnd_imp(
        imp: *const c_double, sh: *const c_double, hr: *const c_double, dens: *const c_double,
        cond: *const c_double, angfreq: c_double, nn: c_int, n: c_int, cdim: c_int,
        coordsys: c_int, velasm: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn aca_block_bnd_slip(
        sh: *const c_double, hr: *const c_double, dens: *const c_double, temp: *const c_double,
        angfreq: c_double, cond: *const c_double, wallt: *const c_double, c1: c_double,
        c2: c_double, nn: c_int, n: c_int, cdim: c_int, coordsys: c_int, velasm: c_int,
    );
    fn aca_block_setbcs(cdim: c_int, velasm: c_int);
    #[allow(clippy::too_many_arguments)]
    fn aca_traction_precond(
        cdim: c_int, angfreq: c_double, sh: *mut c_double, hr: *mut c_double,
        dens: *mut c_double, temp: *mut c_double, cond: *mut c_double, visc: *mut c_double,
        lam: *mut c_double, bulk: *mut c_double, press: *mut c_double, nn: c_int,
    );
    fn aca_block_solve(
        which: c_int, maxit: c_int, tol: c_double, cdim: c_int, jdeg: c_int,
        utilprev: c_int, norm: *mut c_double,
    );
    fn aca_shell_interface(gap: *mut c_int, cdim: c_int, dofs: c_int, angfreq: c_double, nn: c_int);
    #[allow(clippy::too_many_arguments)]
    fn aca_local_interface(
        lsm: *mut c_double, lf: *mut c_double, angfreq: c_double, sh: *const c_double,
        hr: *const c_double, dens: *const c_double, press: *const c_double,
        temp: *const c_double, cond: *const c_double, imp: *const c_double,
        lod: *const c_double, nn: c_int, n: c_int, dofs: c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn aca_slip_matrix(
        lsm: *mut c_double, lf: *mut c_double, angfreq: c_double, sh: *const c_double,
        hr: *const c_double, dens: *const c_double, cond: *const c_double,
        press: *const c_double, temp: *const c_double, wallt: *const c_double,
        wallv: *const c_double, c1: c_double, c2: c_double, c3: c_double,
        nn: c_int, n: c_int, dofs: c_int,
    );
    fn aca_avg_velocity(n: c_int, wallvelo: *const c_double, wsz: c_int, re: *mut c_double, im: *mut c_double);
    fn aca_elem_code() -> c_int;
}

/// First-mesh SAVE state: every driver `SAVE` array plus SAVE scalars.
/// `n` is the first call's `MaxElementDOFs`, `dd` the `Dofs`, `nbc` the
/// `NumberOfBCs`; reallocated when unset or `MeshChanged`.
struct AcousticsState {
    n: usize,
    dd: usize,
    nbc: usize,
    allocations_done: bool,
    solver_calls: i32,
    first_visit: bool,
    prev_norm: f64,
    current_done_time: i32,
    local_force: Vec<c_double>,
    tmp_work: Vec<c_double>,
    local_stiff: Vec<c_double>,
    load: Vec<c_double>,
    heat_source: Vec<c_double>,
    specific_heat: Vec<c_double>,
    heat_ratio: Vec<c_double>,
    density: Vec<c_double>,
    pressure: Vec<c_double>,
    temperature: Vec<c_double>,
    conductivity: Vec<c_double>,
    viscosity: Vec<c_double>,
    lambda: Vec<c_double>,
    bulk_viscosity: Vec<c_double>,
    impedance: Vec<c_double>,
    wall_temperature: Vec<c_double>,
    wall_velocity: Vec<c_double>,
    total_force: Vec<c_double>,
    total_moment: Vec<c_double>,
    total_area: Vec<c_double>,
    moment_output: Vec<c_int>,
    gap_indexes: Vec<c_int>,
}

impl AcousticsState {
    const fn new() -> Self {
        Self {
            n: 0,
            dd: 0,
            nbc: 0,
            allocations_done: false,
            solver_calls: 0,
            first_visit: true,
            prev_norm: 0.0,
            current_done_time: 0,
            local_force: Vec::new(),
            tmp_work: Vec::new(),
            local_stiff: Vec::new(),
            load: Vec::new(),
            heat_source: Vec::new(),
            specific_heat: Vec::new(),
            heat_ratio: Vec::new(),
            density: Vec::new(),
            pressure: Vec::new(),
            temperature: Vec::new(),
            conductivity: Vec::new(),
            viscosity: Vec::new(),
            lambda: Vec::new(),
            bulk_viscosity: Vec::new(),
            impedance: Vec::new(),
            wall_temperature: Vec::new(),
            wall_velocity: Vec::new(),
            total_force: Vec::new(),
            total_moment: Vec::new(),
            total_area: Vec::new(),
            moment_output: Vec::new(),
            gap_indexes: Vec::new(),
        }
    }
}

static STATE: Mutex<AcousticsState> = Mutex::new(AcousticsState::new());

/// Grow a `Vec<f64>` reservation, calling the original
/// `Fatal('AcousticsSolver', 'Memory allocation error.')` on failure.
fn reserve_f64(v: &mut Vec<c_double>, need: usize) {
    if v.try_reserve_exact(need.saturating_sub(v.capacity())).is_err() {
        unsafe { aca_fatal_alloc() };
        std::process::abort();
    }
    v.resize(need, 0.0);
}

/// Grow a `Vec<i32>` reservation, same failure path.
fn reserve_i32(v: &mut Vec<c_int>, need: usize) {
    if v.try_reserve_exact(need.saturating_sub(v.capacity())).is_err() {
        unsafe { aca_fatal_alloc() };
        std::process::abort();
    }
    v.resize(need, 0);
}

/// Read a Fortran `CHARACTER(LEN=MAX_NAME_LEN)` out-buffer as bytes.
unsafe fn out_string(s: *const c_char, n: usize) -> Vec<u8> {
    let raw = std::slice::from_raw_parts(s as *const u8, n);
    let len = raw.iter().position(|&b| b == 0).unwrap_or(n);
    raw[..len].to_vec()
}

unsafe fn sim_string(which_out: &mut [u8; 128], got: &mut c_int, is_sim: bool, which: c_int) {
    let mut buf = [0 as c_char; 128];
    if is_sim {
        aca_sim_type(buf.as_mut_ptr(), 128, got);
    } else {
        aca_solver_string(which, buf.as_mut_ptr(), 128, got);
    }
    let v = out_string(buf.as_ptr(), 128);
    let m = v.len().min(127);
    which_out[..m].copy_from_slice(&v[..m]);
    which_out[m] = 0;
}

unsafe fn flow_vals_mut<'a>() -> (&'a mut [c_double], &'a mut [c_int]) {
    let mut nv = 0 as c_int;
    let mut np = 0 as c_int;
    aca_var_len(&mut nv, &mut np);
    let v = std::slice::from_raw_parts_mut(aca_flow_vals(), nv as usize);
    let p = std::slice::from_raw_parts_mut(aca_flow_perm(), np as usize);
    (v, p)
}

unsafe fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let mut r = 0.0;
    let mut i = 0.0;
    aca_cmul(a.0, a.1, b.0, b.1, &mut r, &mut i);
    (r, i)
}

unsafe fn cdiv(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let mut r = 0.0;
    let mut i = 0.0;
    aca_cdiv(a.0, a.1, b.0, b.1, &mut r, &mut i);
    (r, i)
}

/// `acousticssolver_init_` — mirrors `AcousticsSolver_init`: complex
/// linear system flag plus the `dim`-dependent default `Variable` string.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "acousticssolver_init_"]
pub unsafe extern "C" fn acoustics_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    aca_init_begin(solver);
    aca_info_init();
    aca_init_add_complex();
    let dim = aca_coord_dim();
    aca_init_add_variable(if dim == 2 { 1 } else { 0 });
}
/// Complex pair helpers over interleaved `(re, im)` buffers.
#[inline]
fn pair(buf: &[f64], i: usize) -> (f64, f64) {
    (buf[2 * i], buf[2 * i + 1])
}

#[inline]
fn set_pair(buf: &mut [f64], i: usize, v: (f64, f64)) {
    buf[2 * i] = v.0;
    buf[2 * i + 1] = v.1;
}

/// `LCondensate` — static condensation of bubble dofs (no `nb`), Rust
/// transcription: index gathers/scatters are exact integer arithmetic;
/// `ComplexInvertMatrix` and the nested complex `MATMUL`s run in Fortran.
unsafe fn lcondensate(n: usize, dim: usize, k: &mut [f64], f: &mut [f64]) {
    let nl = (dim + 2) * n;
    let nb = dim * n;
    let mut kbb = vec![0.0; 2 * nb * nb];
    let mut kbl = vec![0.0; 2 * nb * nl];
    let mut klb = vec![0.0; 2 * nl * nb];
    let mut fb = vec![0.0; 2 * nb];

    // Column-major gathers out of the interleaved full matrix. `k` holds
    // `ld x ld` complex pairs in column-major order, `ld = nl + nb`.
    let ld = nl + nb;
    let at = |r: usize, c: usize| -> (f64, f64) { pair(k, c * ld + r) };
    for c in 0..nb {
        for r in 0..nb {
            let v = at(nl + r, nl + c);
            kbb[(c * nb + r) * 2] = v.0;
            kbb[(c * nb + r) * 2 + 1] = v.1;
        }
        for r in 0..nl {
            let w = at(r, nl + c);
            klb[(c * nl + r) * 2] = w.0;
            klb[(c * nl + r) * 2 + 1] = w.1;
        }
        let v = pair(f, nl + c);
        fb[2 * c] = v.0;
        fb[2 * c + 1] = v.1;
    }
    for r in 0..nb {
        for c in 0..nl {
            let v = at(nl + r, c);
            kbl[(c * nb + r) * 2] = v.0;
            kbl[(c * nb + r) * 2 + 1] = v.1;
        }
    }
    aca_cinv(nb as c_int, kbb.as_mut_ptr());
    // F(1:nl) -= Klb * (Kbb * Fb)
    let mut t1 = vec![0.0; 2 * nb];
    let mut t2 = vec![0.0; 2 * nl];
    aca_cmatmul(nb as c_int, nb as c_int, 1, kbb.as_ptr(), fb.as_ptr(), t1.as_mut_ptr());
    aca_cmatmul(nl as c_int, nb as c_int, 1, klb.as_ptr(), t1.as_ptr(), t2.as_mut_ptr());
    for i in 0..nl {
        let (a, b) = pair(f, i);
        let (c, d) = pair(&t2, i);
        set_pair(f, i, (a - c, b - d));
    }
    // K(1:nl,1:nl) -= Klb * (Kbb * Kbl)
    let mut t3 = vec![0.0; 2 * nb * nl];
    let mut t4 = vec![0.0; 2 * nl * nl];
    aca_cmatmul(nb as c_int, nb as c_int, nl as c_int, kbb.as_ptr(), kbl.as_ptr(), t3.as_mut_ptr());
    aca_cmatmul(nl as c_int, nb as c_int, nl as c_int, klb.as_ptr(), t3.as_ptr(), t4.as_mut_ptr());
    for c in 0..nl {
        for r in 0..nl {
            let (a, b) = pair(k, c * ld + r);
            let (e, g) = pair(&t4, c * nl + r);
            let o = (c * ld + r) * 2;
            k[o] = a - e;
            k[o + 1] = b - g;
        }
    }
}
/// `LCondensateBubbles` — same static condensation with `nb` bubble dofs
/// per node. Identical structure to `lcondensate`.
unsafe fn lcondensate_bubbles(n: usize, dim: usize, nbb: usize, k: &mut [f64], f: &mut [f64]) {
    let nl = (dim + 2) * n;
    let nb = dim * nbb;
    let mut kbb = vec![0.0; 2 * nb * nb];
    let mut kbl = vec![0.0; 2 * nb * nl];
    let mut klb = vec![0.0; 2 * nl * nb];
    let mut fb = vec![0.0; 2 * nb];
    let ld = nl + nb;
    let at = |r: usize, c: usize| -> (f64, f64) { pair(k, c * ld + r) };
    for c in 0..nb {
        for r in 0..nb {
            let v = at(nl + r, nl + c);
            kbb[(c * nb + r) * 2] = v.0;
            kbb[(c * nb + r) * 2 + 1] = v.1;
        }
        for r in 0..nl {
            let w = at(r, nl + c);
            klb[(c * nl + r) * 2] = w.0;
            klb[(c * nl + r) * 2 + 1] = w.1;
        }
        let v = pair(f, nl + c);
        fb[2 * c] = v.0;
        fb[2 * c + 1] = v.1;
    }
    for r in 0..nb {
        for c in 0..nl {
            let v = at(nl + r, c);
            kbl[(c * nb + r) * 2] = v.0;
            kbl[(c * nb + r) * 2 + 1] = v.1;
        }
    }
    aca_cinv(nb as c_int, kbb.as_mut_ptr());
    let mut t1 = vec![0.0; 2 * nb];
    let mut t2 = vec![0.0; 2 * nl];
    aca_cmatmul(nb as c_int, nb as c_int, 1, kbb.as_ptr(), fb.as_ptr(), t1.as_mut_ptr());
    aca_cmatmul(nl as c_int, nb as c_int, 1, klb.as_ptr(), t1.as_ptr(), t2.as_mut_ptr());
    for i in 0..nl {
        let (a, b) = pair(f, i);
        let (c, d) = pair(&t2, i);
        set_pair(f, i, (a - c, b - d));
    }
    let mut t3 = vec![0.0; 2 * nb * nl];
    let mut t4 = vec![0.0; 2 * nl * nl];
    aca_cmatmul(nb as c_int, nb as c_int, nl as c_int, kbb.as_ptr(), kbl.as_ptr(), t3.as_mut_ptr());
    aca_cmatmul(nl as c_int, nb as c_int, nl as c_int, klb.as_ptr(), t3.as_ptr(), t4.as_mut_ptr());
    for c in 0..nl {
        for r in 0..nl {
            let (a, b) = pair(k, c * ld + r);
            let (e, g) = pair(&t4, c * nl + r);
            let o = (c * ld + r) * 2;
            k[o] = a - e;
            k[o + 1] = b - g;
        }
    }
}
/// `LocalMatrix` — element stiffness and force, Rust transcription of the
/// original internal procedure. Control flow, real arithmetic and indexing
/// are transcribed statement for statement; complex division goes through
/// `aca_cdiv`, `SUM` reductions and `ElementInfo`/`GaussPoints` through
/// wrappers. `MyBubbles` is identically `.FALSE.` (dead arms omitted).
#[allow(clippy::too_many_arguments)]
unsafe fn local_matrix(
    stiff: &mut [f64], force: &mut [f64], lsm_stride: usize, af: f64,
    sh: &[f64], hr: &[f64], dens: &[f64], temp_: &[f64], cond: &[f64], visc: &[f64],
    lam: &[f64], hs: &[f64], load: &[f64], bubbles: bool, mini: bool, nn: usize,
    n: usize, dofs: usize, nb: usize, dim: usize, relorder: i32, axisym: bool,
) {
    let ldc = dofs - 2;
    let ld = n * ldc;
    let li = |r: usize, c: usize| -> usize { (c * ld + r) * 2 };
    let mut lstiff = vec![0.0; 2 * ld * ld];
    let mut lforce = vec![0.0; 2 * n * ldc];
    let mut basis = vec![0.0; 2 * n];
    let mut dbasis = vec![0.0; 2 * n * 3];
    let mut ddbasis = vec![0.0; n * 9];
    let (nbasis, nq) = if bubbles || mini {
        let nb2 = if mini { n + 1 } else { 2 * n };
        (nb2, aca_gauss_begin(1, 0))
    } else {
        (n + nb, aca_gauss_begin(2, relorder))
    };
    for tt in 1..=nq as c_int {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        aca_gauss_point(tt, &mut u, &mut v, &mut w, &mut s);
        let mut detj = 0.0;
        let mut estat = 0;
        let eimode = if nb > 0 { 0 } else { 1 };
        aca_eleminfo(
            eimode, u, v, w, n as c_int, 0, basis.as_mut_ptr(),
            dbasis.as_mut_ptr(), ddbasis.as_mut_ptr(), &mut detj,
            if bubbles { 1 } else { 0 }, &mut estat,
        );
        // Fortran stores dBasisdx column-major (2n rows); transpose once
        // into a row-major work copy so the transcribed indexing matches.
        let mut dbasis_t = vec![0.0; 2 * n * 3];
        for qq in 0..2 * n {
            for jj in 0..3 {
                dbasis_t[qq * 3 + jj] = dbasis[jj * 2 * n + qq];
            }
        }
        let mut dbasis = dbasis_t;
        if mini {
            for j in 0..3 {
                dbasis[n * 3 + j] = 0.0;
            }
            match aca_elem_code() {
                504 => {
                    basis[n] = 1.0;
                    for i in 0..n {
                        basis[n] *= basis[i];
                    }
                    for j in 0..dim {
                        dbasis[n * 3 + j] = dbasis[j] * basis[1] * basis[2] * basis[3]
                            + dbasis[3 + j] * basis[0] * basis[2] * basis[3]
                            + dbasis[2 * 3 + j] * basis[0] * basis[1] * basis[3]
                            + dbasis[3 * 3 + j] * basis[0] * basis[1] * basis[2];
                    }
                }
                404 => {
                    basis[n] = basis[0] * basis[2];
                    for j in 0..dim {
                        dbasis[n * 3 + j] =
                            dbasis[j] * basis[2] + dbasis[2 * 3 + j] * basis[0];
                    }
                }
                303 => {
                    basis[n] = basis[0] * basis[1] * basis[2];
                    for j in 0..dim {
                        dbasis[n * 3 + j] = dbasis[j] * basis[1] * basis[2]
                            + dbasis[3 + j] * basis[0] * basis[2]
                            + dbasis[2 * 3 + j] * basis[0] * basis[1];
                    }
                }
                706 => {
                    basis[n] = (basis[0] + basis[1] + basis[2])
                        * (basis[3] + basis[4] + basis[5])
                        * (basis[0] + basis[3])
                        * (basis[1] + basis[4])
                        * (basis[2] + basis[5]);
                    for j in 0..dim {
                        dbasis[n * 3 + j] = (dbasis[j]
                            + dbasis[3 + j]
                            + dbasis[2 * 3 + j])
                            * (basis[3] + basis[4] + basis[5])
                            * (basis[0] + basis[3])
                            * (basis[1] + basis[4])
                            * (basis[2] + basis[5])
                            + (dbasis[3 * 3 + j] + dbasis[4 * 3 + j] + dbasis[5 * 3 + j])
                                * (basis[0] + basis[1] + basis[2])
                                * (basis[0] + basis[3])
                                * (basis[1] + basis[4])
                                * (basis[2] + basis[5])
                            + (basis[0] + basis[1] + basis[2])
                                * (basis[3] + basis[4] + basis[5])
                                * (dbasis[j] + dbasis[3 * 3 + j])
                                * (basis[1] + basis[4])
                                * (basis[2] + basis[5])
                            + (basis[0] + basis[1] + basis[2])
                                * (basis[3] + basis[4] + basis[5])
                                * (dbasis[3 + j] + dbasis[4 * 3 + j])
                                * (basis[0] + basis[3])
                                * (basis[2] + basis[5])
                            + (basis[0] + basis[1] + basis[2])
                                * (basis[3] + basis[4] + basis[5])
                                * (dbasis[2 * 3 + j] + dbasis[5 * 3 + j])
                                * (basis[0] + basis[3])
                                * (basis[1] + basis[4]);
                    }
                }
                808 => {
                    basis[n] = basis[0] * basis[6];
                    for j in 0..dim {
                        dbasis[n * 3 + j] =
                            dbasis[j] * basis[6] + dbasis[6 * 3 + j] * basis[0];
                    }
                }
                other => {
                    aca_error_minibubbles(other);
                }
            }
        }
        s *= detj;
        let r = if axisym {
            let rr = aca_nodes_dotx(n as c_int, basis.as_ptr());
            s *= rr;
            rr
        } else {
            0.0
        };
        let cv = aca_dot(n as c_int, sh.as_ptr(), basis.as_ptr());
        let gamma = aca_dot(n as c_int, hr.as_ptr(), basis.as_ptr());
        let rho0 = aca_dot(n as c_int, dens.as_ptr(), basis.as_ptr());
        let t0 = aca_dot(n as c_int, temp_.as_ptr(), basis.as_ptr());
        let kappa = aca_dot(n as c_int, cond.as_ptr(), basis.as_ptr());
        let mu = aca_dot(n as c_int, visc.as_ptr(), basis.as_ptr());
        let la = aca_dot(n as c_int, lam.as_ptr(), basis.as_ptr());
        let p0 = (gamma - 1.0) * cv * rho0 * t0;
        let k1 = kappa / (rho0 * af * (gamma - 1.0) * cv);
        let k2 = 1.0 / (gamma - 1.0);
        let k3 = 1.0 / (af * af);
        let _schur = cdiv(
            (1.0, af / p0 * (2.0 * mu + la)),
            (1.0, af / p0 * la),
        );
        let _c1 = cdiv(
            (1.0, af / p0 * (2.0 * mu + la)),
            (1.0, af / p0 * la),
        );
        let f1 = k3 * aca_dot(n as c_int, hs.as_ptr(), basis.as_ptr());
        let f2 = k3 * aca_dot(n as c_int, hs.as_ptr().add(nn), basis.as_ptr());
        let mut lv = [0.0; 6];
        for (i, l) in lv.iter_mut().enumerate().take(2 * dim) {
            *l = 1.0 / af * aca_dot(n as c_int, load.as_ptr().add(i * nn), basis.as_ptr());
        }
        for pp in 0..n {
            for qq in 0..n {
                let bq = basis[qq];
                let bp = basis[pp];
                let ax = -k2 * bq * bp;
                let mut ay = 0.0;
                for i in 0..dim {
                    let o = li(pp * (dim + 2) + dim + 1, qq * (dim + 2) + i);
                    lstiff[o + 1] += dbasis[qq * 3 + i] * bp * s;
                    if i == 0 && axisym {
                        lstiff[o + 1] += (1.0 / r) * bq * bp * s;
                    }
                    ay += k1 * dbasis[qq * 3 + i] * dbasis[pp * 3 + i];
                }
                let o = li(pp * (dim + 2) + dim, qq * (dim + 2) + dim);
                lstiff[o] += s * ax;
                lstiff[o + 1] += s * ay;
                let z = cdiv((1.0, 0.0), (1.0, af / p0 * la));
                let o2 = li(pp * (dim + 2) + dim, qq * (dim + 2) + dim + 1);
                let t0 = s * z.0;
                let t1 = s * z.1;
                lstiff[o2] += t0 * bq * bp;
                lstiff[o2 + 1] += t1 * bq * bp;
                let nre = (s * rho0) * af;
                let z2 = cdiv((nre, 0.0), (p0 / af, la));
                let o3 = li(pp * (dim + 2) + dim + 1, qq * (dim + 2) + dim + 1);
                lstiff[o3] -= z2.0 * bq * bp;
                lstiff[o3 + 1] -= z2.1 * bq * bp;
            }
            let bp = basis[pp];
            let of = (pp * (dim + 2) + dim) * 2;
            lforce[of] += s * bp * -f2;
            lforce[of + 1] += s * bp * f1;
        }
        if bubbles || mini || nb > 0 {
            for pp in 0..n {
                for qq in n..nbasis {
                    let bq = basis[qq];
                    let bp = basis[pp];
                    for i in 0..dim {
                        let o = li(pp * (dim + 2) + dim + 1, qq * dim + 2 * n + i);
                        lstiff[o + 1] += dbasis[qq * 3 + i] * bp * s;
                        if i == 0 && axisym {
                            lstiff[o + 1] += (1.0 / r) * bq * bp * s;
                        }
                    }
                }
            }
        }
        for i in 0..dim {
            for pp in 0..n {
                for qq in 0..n {
                    let bq = basis[qq];
                    let bp = basis[pp];
                    let o = li(pp * (dim + 2) + i, qq * (dim + 2) + dim);
                    lstiff[o + 1] += dbasis[pp * 3 + i] * bq * s;
                    if i == 0 && axisym {
                        lstiff[o + 1] += (1.0 / r) * bp * bq * s;
                    }
                    let o2 = li(pp * (dim + 2) + i, qq * (dim + 2) + dim + 1);
                    lstiff[o2 + 1] += dbasis[pp * 3 + i] * bq * s;
                    if i == 0 && axisym {
                        lstiff[o2 + 1] += (1.0 / r) * bp * bq * s;
                    }
                    let o3 = li(pp * (dim + 2) + i, qq * (dim + 2) + i);
                    lstiff[o3] += bq * bp * s;
                    for j in 0..dim {
                        let c0 = -mu / (af * rho0);
                        let o4 = li(pp * (dim + 2) + i, qq * (dim + 2) + i);
                        lstiff[o4 + 1] += c0 * dbasis[qq * 3 + j] * dbasis[pp * 3 + j] * s;
                        let o5 = li(pp * (dim + 2) + i, qq * (dim + 2) + j);
                        lstiff[o5 + 1] += c0 * dbasis[qq * 3 + i] * dbasis[pp * 3 + j] * s;
                    }
                    if i == 0 && axisym {
                        let o6 = li(pp * (dim + 2) + i, qq * (dim + 2) + i);
                        lstiff[o6 + 1] += (-2.0 * mu / (af * rho0)) * (1.0 / (r * r)) * bq * bp * s;
                    }
                }
                let bp = basis[pp];
                let of = (pp * (dim + 2) + i) * 2;
                lforce[of] -= s * bp * -lv[2 * i + 1];
                lforce[of + 1] -= s * bp * lv[2 * i];
            }
        }
        if bubbles || mini || nb > 0 {
            for i in 0..dim {
                for pp in 0..n {
                    for qq in n..nbasis {
                        let bq = basis[qq];
                        let bp = basis[pp];
                        let o = li(pp * (dim + 2) + i, qq * dim + 2 * n + i);
                        lstiff[o] += bq * bp * s;
                        for j in 0..dim {
                            let c0 = -mu / (af * rho0);
                            let o2 = li(pp * (dim + 2) + i, qq * dim + 2 * n + i);
                            lstiff[o2 + 1] += c0 * dbasis[qq * 3 + j] * dbasis[pp * 3 + j] * s;
                            let o3 = li(pp * (dim + 2) + i, qq * dim + 2 * n + j);
                            lstiff[o3 + 1] += c0 * dbasis[qq * 3 + i] * dbasis[pp * 3 + j] * s;
                        }
                        if i == 0 && axisym {
                            let o4 = li(pp * (dim + 2) + i, qq * dim + 2 * n + i);
                            lstiff[o4 + 1] += (-2.0 * mu / (af * rho0)) * (1.0 / (r * r)) * bq * bp * s;
                        }
                    }
                }
            }
            for i in 0..dim {
                for pp in n..nbasis {
                    for qq in 0..n {
                        let bq = basis[qq];
                        let bp = basis[pp];
                        let o = li(pp * dim + 2 * n + i, qq * (dim + 2) + dim);
                        lstiff[o + 1] += dbasis[pp * 3 + i] * bq * s;
                        if i == 0 && axisym {
                            lstiff[o + 1] += (1.0 / r) * bp * bq * s;
                        }
                        let o2 = li(pp * dim + 2 * n + i, qq * (dim + 2) + dim + 1);
                        lstiff[o2 + 1] += dbasis[pp * 3 + i] * bq * s;
                        if i == 0 && axisym {
                            lstiff[o2 + 1] += (1.0 / r) * bp * bq * s;
                        }
                        let o3 = li(pp * dim + 2 * n + i, qq * (dim + 2) + i);
                        lstiff[o3] += bq * bp * s;
                        for j in 0..dim {
                            let c0 = -mu / (af * rho0);
                            let o4 = li(pp * dim + 2 * n + i, qq * (dim + 2) + i);
                            lstiff[o4 + 1] += c0 * dbasis[qq * 3 + j] * dbasis[pp * 3 + j] * s;
                            let o5 = li(pp * dim + 2 * n + i, qq * (dim + 2) + j);
                            lstiff[o5 + 1] += c0 * dbasis[qq * 3 + i] * dbasis[pp * 3 + j] * s;
                        }
                        if i == 0 && axisym {
                            let o6 = li(pp * dim + 2 * n + i, qq * (dim + 2) + i);
                            lstiff[o6 + 1] += (-2.0 * mu / (af * rho0)) * (1.0 / (r * r)) * bq * bp * s;
                        }
                    }
                    let bp = basis[pp];
                    let of = (pp * dim + 2 * n + i) * 2;
                    lforce[of] -= s * bp * -lv[2 * i + 1];
                    lforce[of + 1] -= s * bp * lv[2 * i];
                }
            }
            for i in 0..dim {
                for pp in n..nbasis {
                    for qq in n..nbasis {
                        let bq = basis[qq];
                        let bp = basis[pp];
                        let o = li(pp * dim + 2 * n + i, qq * dim + 2 * n + i);
                        lstiff[o] += bq * bp * s;
                        for j in 0..dim {
                            let c0 = -mu / (af * rho0);
                            let o2 = li(pp * dim + 2 * n + i, qq * dim + 2 * n + i);
                            lstiff[o2 + 1] += c0 * dbasis[qq * 3 + j] * dbasis[pp * 3 + j] * s;
                            let o3 = li(pp * dim + 2 * n + i, qq * dim + 2 * n + j);
                            lstiff[o3 + 1] += c0 * dbasis[qq * 3 + i] * dbasis[pp * 3 + j] * s;
                        }
                        if i == 0 && axisym {
                            let o4 = li(pp * dim + 2 * n + i, qq * dim + 2 * n + i);
                            lstiff[o4 + 1] += (-2.0 * mu / (af * rho0)) * (1.0 / (r * r)) * bq * bp * s;
                        }
                    }
                }
            }
        }
    }
    if nb > 0 {
        lcondensate_bubbles(n, dim, nb, &mut lstiff, &mut lforce);
    } else if bubbles || mini {
        if mini {
            lcondensate_bubbles(n, dim, 1, &mut lstiff, &mut lforce);
        } else {
            lcondensate(n, dim, &mut lstiff, &mut lforce);
        }
    }
    for pp in 0..n {
        for i in 0..dim + 2 {
            let of = 2 * (dim + 2) * pp + 2 * i;
            force[of] = lforce[(pp * (dim + 2) + i) * 2];
            force[of + 1] = lforce[(pp * (dim + 2) + i) * 2 + 1];
            for qq in 0..n {
                for j in 0..dim + 2 {
                    let (lr, lc) = (pp * (dim + 2) + i, qq * (dim + 2) + j);
                    let (zr, zc) = (2 * (dim + 2) * pp + 2 * i, 2 * (dim + 2) * qq + 2 * j);
                    let o = li(lr, lc);
                    stiff[zc * lsm_stride + zr] = lstiff[o];
                    stiff[(zc + 1) * lsm_stride + zr] = -lstiff[o + 1];
                    stiff[zc * lsm_stride + zr + 1] = lstiff[o + 1];
                    stiff[(zc + 1) * lsm_stride + zr + 1] = lstiff[o];
                }
            }
        }
    }
}
/// `LocalMatrixBoundary` — impedance/traction boundary assembly, Rust
/// transcription. Same conventions as `local_matrix`.
#[allow(clippy::too_many_arguments)]
unsafe fn local_matrix_boundary(
    stiff: &mut [f64], force: &mut [f64], lsm_stride: usize, af: f64,
    sh: &[f64], hr: &[f64], dens: &[f64], press: &[f64], temp_: &[f64], cond: &[f64],
    imp: &[f64], load: &[f64], nn: usize, n: usize, dofs: usize, dim: usize, axisym: bool,
) {
    let ldc = dofs - 2;
    let ld = n * ldc;
    let li = |r: usize, c: usize| -> usize { (c * ld + r) * 2 };
    let mut lstiff = vec![0.0; 2 * ld * ld];
    let mut lforce = vec![0.0; 2 * ld];
    let mut basis = vec![0.0; 2 * n];
    let mut dbasis = vec![0.0; 2 * n * 3];
    let mut ddbasis = vec![0.0; n * 9];
    let nq = aca_gauss_begin(0, 0);
    for tt in 1..=nq as c_int {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        aca_gauss_point(tt, &mut u, &mut v, &mut w, &mut s);
        let mut detj = 0.0;
        let mut estat = 0;
        aca_eleminfo(
            1, u, v, w, n as c_int, 0, basis.as_mut_ptr(),
            dbasis.as_mut_ptr(), ddbasis.as_mut_ptr(), &mut detj, 0, &mut estat,
        );
        s *= detj;
        if axisym {
            s *= aca_nodes_dotx(n as c_int, basis.as_ptr());
        }
        let cv = aca_dot(n as c_int, sh.as_ptr(), basis.as_ptr());
        let gamma = aca_dot(n as c_int, hr.as_ptr(), basis.as_ptr());
        let rho0 = aca_dot(n as c_int, dens.as_ptr(), basis.as_ptr());
        let _p0 = aca_dot(n as c_int, press.as_ptr(), basis.as_ptr());
        let _t0 = aca_dot(n as c_int, temp_.as_ptr(), basis.as_ptr());
        let kappa = aca_dot(n as c_int, cond.as_ptr(), basis.as_ptr());
        let k1 = kappa / (rho0 * af * (gamma - 1.0) * cv);
        let mut nrm = [0.0; 3];
        aca_normal(u, v, nrm.as_mut_ptr());
        let imp1 = 1.0 / (af * rho0) * aca_dot(n as c_int, imp.as_ptr(), basis.as_ptr());
        let imp2 = 1.0 / (af * rho0) * aca_dot(n as c_int, imp.as_ptr().add(nn), basis.as_ptr());
        let imp3 = aca_dot(n as c_int, imp.as_ptr().add(2 * nn), basis.as_ptr());
        let imp4 = aca_dot(n as c_int, imp.as_ptr().add(3 * nn), basis.as_ptr());
        let mut lv = [0.0; 6];
        for (i, l) in lv.iter_mut().enumerate().take(2 * dim) {
            *l = 1.0 / (af * rho0) * aca_dot(n as c_int, load.as_ptr().add(i * nn), basis.as_ptr());
        }
        for pp in 0..n {
            for i in 0..dim {
                for qq in 0..n {
                    for j in 0..dim {
                        let o = li(pp * (dim + 2) + i, qq * (dim + 2) + j);
                        let x = basis[qq] * nrm[j] * basis[pp] * nrm[i] * s;
                        lstiff[o] += -imp2 * x;
                        lstiff[o + 1] += imp1 * x;
                    }
                }
                let of = (pp * (dim + 2) + i) * 2;
                lforce[of] -= s * basis[pp] * -lv[2 * i + 1];
                lforce[of + 1] -= s * basis[pp] * lv[2 * i];
            }
        }
        for pp in 0..n {
            for qq in 0..n {
                let z = cmul((0.0, s * k1), (imp3, imp4));
                let o = li(pp * (dim + 2) + dim, qq * (dim + 2) + dim);
                lstiff[o] -= z.0 * basis[qq] * basis[pp];
                lstiff[o + 1] -= z.1 * basis[qq] * basis[pp];
            }
        }
    }
    for pp in 0..n {
        for i in 0..dim + 2 {
            let of = 2 * (dim + 2) * pp + 2 * i;
            force[of] = lforce[(pp * (dim + 2) + i) * 2];
            force[of + 1] = lforce[(pp * (dim + 2) + i) * 2 + 1];
            for qq in 0..n {
                for j in 0..dim + 2 {
                    let (lr, lc) = (pp * (dim + 2) + i, qq * (dim + 2) + j);
                    let (zr, zc) = (2 * (dim + 2) * pp + 2 * i, 2 * (dim + 2) * qq + 2 * j);
                    let o = li(lr, lc);
                    stiff[zc * lsm_stride + zr] = lstiff[o];
                    stiff[(zc + 1) * lsm_stride + zr] = -lstiff[o + 1];
                    stiff[zc * lsm_stride + zr + 1] = lstiff[o + 1];
                    stiff[(zc + 1) * lsm_stride + zr + 1] = lstiff[o];
                }
            }
        }
    }
}
/// `acousticssolver_` — mirrors `AcousticsSolver`. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient flag
/// travel through to the (unused) dummies, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "acousticssolver_"]
pub unsafe extern "C" fn acoustics_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    aca_begin(model, solver);
    let mut guard = STATE.lock().unwrap();
    let st: &mut AcousticsState = &mut guard;
    if aca_assoc_matrix() == 0 {
        return;
    }
    let dofs = aca_flow_dofs() as usize;
    let local_nodes = aca_local_nodes();
    if local_nodes <= 0 {
        return;
    }
    let norm0 = aca_var_norm();
    let dim = aca_coord_dim() as usize;
    let vel_dofs = dim * 2;
    if dofs != vel_dofs + 4 {
        aca_warn_dofs();
    }
    // Scanning option for coupled simulations.
    let mut sim_raw = [0u8; 128];
    let mut got = 0;
    sim_string(&mut sim_raw, &mut got, true, 0);
    let sim_len = sim_raw.iter().position(|&b| b == 0).unwrap_or(128);
    let scanning = &sim_raw[..sim_len] == b"scanning";
    if scanning {
        let (mut found, mut v1) = (0, 0.0);
        aca_time_var(&mut found, &mut v1);
        if found == 0 {
            aca_fatal_timevar();
            return;
        }
        let first_ts = v1.round() as i32 == 1;
        if first_ts {
            // FirstTimeStep = .TRUE. (locals only below)
        } else if (v1.round() as i32 - st.current_done_time) == 2 {
            st.current_done_time += 1;
        }
        let _ = first_ts;
    }
    // BEM coupling status-file dance; DD preconditioning / Helmholtz lookup.
    let (mut bv, mut gotit) = (0, 0);
    aca_sim_logical(1, &mut bv, &mut gotit);
    let bem_coupling = bv != 0;
    st.solver_calls += 1;
    let mut dd_precond = false;
    let mut h_found = false;
    if bem_coupling {
        aca_bem_pre(if st.first_visit { 1 } else { 0 });
    } else {
        let (mut dv, mut dg) = (0, 0);
        aca_sim_logical(2, &mut dv, &mut dg);
        dd_precond = dv != 0;
        if dd_precond {
            aca_dd_pres();
            h_found = true;
        } else {
            let mut hf = 0;
            aca_find_pres(&mut hf);
            h_found = hf != 0;
            if !h_found {
                aca_info_no_pres();
            }
        }
    }
    let _ = dd_precond;
    // Relative integration order diagnostic (Level 6, suppressed).
    let (mut relorder, mut relgot) = (0, 0);
    aca_solver_int(1, &mut relorder, &mut relgot);
    if relgot != 0 {
        aca_relorder_diag(relorder);
    }
    if bem_coupling && !st.first_visit {
        aca_bem_nodes();
    }
    // Block preconditioning options.
    let (mut bpv, mut bpg) = (0, 0);
    aca_solver_logical(1, &mut bpv, &mut bpg);
    let block_precond = bpv != 0;
    let mut outer_which = 0;
    let mut tolerance = 0.0;
    let mut max_iterations = 0;
    let mut vel_assembly = 0;
    if block_precond {
        aca_info_block();
        let mut oname = [0u8; 128];
        let mut og = 0;
        sim_string(&mut oname, &mut og, false, 1);
        let olen = oname.iter().position(|&b| b == 0).unwrap_or(128);
        let onm = &oname[..olen];
        if og == 0 {
            outer_which = 0;
        } else if onm == b"gcr" {
            outer_which = 1;
        } else if onm == b"bicgstab" {
            outer_which = 2;
        } else if onm == b"bicgstabl" {
            outer_which = 3;
        } else {
            outer_which = 0;
        }
        let (mut tr, mut tg) = (0.0, 0);
        aca_solver_creal(1, &mut tr, &mut tg);
        let (mut t2, mut t2g) = (0.0, 0);
        aca_solver_creal(2, &mut t2, &mut t2g);
        tolerance = if tg != 0 { tr * t2 } else { t2 };
        let (mut mi, mut mig) = (0, 0);
        aca_solver_int(2, &mut mi, &mut mig);
        if mig != 0 {
            max_iterations = mi;
        } else {
            let (mut mi2, mut mi2g) = (0, 0);
            aca_solver_int(3, &mut mi2, &mut mi2g);
            max_iterations = mi2;
            let _ = mi2g;
        }
        let (mut va, mut vag) = (0, 0);
        aca_solver_logical(5, &mut va, &mut vag);
        vel_assembly = if vag != 0 { va } else { 0 };
    }
    // Permanent storage, first time or on mesh change.
    let mesh_changed = aca_mesh_changed() != 0;
    if !st.allocations_done || mesh_changed {
        let nn = aca_max_elem_dofs() as usize;
        let nbc = aca_model_num_bcs() as usize;
        st.n = nn;
        st.dd = dofs;
        st.nbc = nbc;
        reserve_f64(&mut st.local_force, dofs * nn);
        reserve_f64(&mut st.tmp_work, nn);
        reserve_f64(&mut st.local_stiff, dofs * nn * dofs * nn);
        reserve_f64(&mut st.load, 6 * nn);
        reserve_f64(&mut st.heat_source, 2 * nn);
        reserve_f64(&mut st.specific_heat, nn);
        reserve_f64(&mut st.heat_ratio, nn);
        reserve_f64(&mut st.density, nn);
        reserve_f64(&mut st.pressure, nn);
        reserve_f64(&mut st.temperature, nn);
        reserve_f64(&mut st.conductivity, nn);
        reserve_f64(&mut st.viscosity, nn);
        reserve_f64(&mut st.lambda, nn);
        reserve_f64(&mut st.bulk_viscosity, nn);
        reserve_f64(&mut st.impedance, 4 * nn);
        reserve_f64(&mut st.wall_temperature, nn);
        reserve_f64(&mut st.wall_velocity, 6 * nn);
        reserve_f64(&mut st.total_force, nbc * 6);
        reserve_f64(&mut st.total_moment, nbc * 6);
        reserve_f64(&mut st.total_area, nbc);
        reserve_i32(&mut st.moment_output, nbc);
        reserve_i32(&mut st.gap_indexes, nn);
        aca_alloc_nodes(nn as c_int);
        if block_precond {
            aca_block_create();
        }
        st.allocations_done = true;
    }
    let nn = st.n;
    // Block-preconditioner zeroing.
    if block_precond {
        aca_block_zero(vel_assembly);
    }
    aca_set_eqname();
    let (mut bub, mut bubg) = (0, 0);
    aca_solver_logical(2, &mut bub, &mut bubg);
    let bubbles = bub != 0;
    let (mut mb, mut mbg) = (0, 0);
    aca_solver_logical(3, &mut mb, &mut mbg);
    let mini_bubbles = mb != 0;
    let _ = (bubg, mbg);
    let coordsys = aca_cur_coordsys();
    aca_check_coordsys();
    let axisym = coordsys == 4;
    let (mut af, mut got_af, mut got_f) = (0.0, 0, 0);
    aca_angular_freq(&mut af, &mut got_af, &mut got_f);
    let _ = (got_af, got_f);
    // Utilize-previous-solution rescaling (first pass).
    let (mut upv, mut upg) = (0, 0);
    aca_solver_logical(4, &mut upv, &mut upg);
    let utilize_prev = upv != 0;
    if utilize_prev {
        let nbulk = aca_num_bulk() as c_int;
        let mut visited = vec![false; aca_model_num_nodes() as usize];
        for t in 1..=nbulk {
            let (mut n, mut bodyid) = (0, 0);
            aca_begin_direct(t, &mut n, &mut bodyid);
            if aca_check_eqn() == 0 {
                continue;
            }
            let n = n as usize;
            let k = aca_body_material();
            aca_use_material(k);
            let mut g = 0;
            aca_mat_real(1, n as c_int, st.specific_heat.as_mut_ptr(), &mut g);
            aca_mat_real(2, n as c_int, st.heat_ratio.as_mut_ptr(), &mut g);
            aca_mat_real(3, n as c_int, st.density.as_mut_ptr(), &mut g);
            aca_mat_real(4, n as c_int, st.temperature.as_mut_ptr(), &mut g);
            aca_mat_real(6, n as c_int, st.viscosity.as_mut_ptr(), &mut g);
            for i in 0..n {
                st.lambda[i] = -2.0 / 3.0 * st.viscosity[i];
            }
            aca_mat_real(7, n as c_int, st.bulk_viscosity.as_mut_ptr(), &mut g);
            if g != 0 {
                for i in 0..n {
                    st.lambda[i] = st.bulk_viscosity[i] - 2.0 / 3.0 * st.viscosity[i];
                }
            }
            let (flow, perm) = flow_vals_mut();
            for i in 0..n {
                let node = aca_nidx_val((i + 1) as c_int) as usize;
                let j = perm[node - 1] as usize;
                if j == 0 || visited[j - 1] {
                    continue;
                }
                visited[j - 1] = true;
                let b = (j - 1) * dofs + vel_dofs;
                let (mut g1, mut g2, mut g3, mut g4) = (0.0, 0.0, 0.0, 0.0);
                aca_pre_rescale(
                    st.heat_ratio[i], st.specific_heat[i], st.density[i], st.temperature[i],
                    st.viscosity[i], st.lambda[i], af, flow[b], flow[b + 1], flow[b + 2],
                    flow[b + 3], &mut g1, &mut g2, &mut g3, &mut g4,
                );
                flow[b] = g1;
                flow[b + 1] = g2;
                flow[b + 2] = g3;
                flow[b + 3] = g4;
            }
        }
    }
    let at = aca_cputime();
    let mut at0 = aca_realtime();
    let mut time_estimated = false;
    aca_info_blank();
    aca_info_dashes();
    aca_info_freq(af / (2.0 * std::f64::consts::PI));
    aca_info_dashes();
    aca_info_blank();
    aca_info_starting();
    aca_default_start();
    aca_default_initialize();
    // Bulk assembly over active elements.
    let num_active = aca_num_active();
    for t in 1..=num_active {
        if aca_realtime() - at0 > 1.0 {
            if !time_estimated && aca_cputime() - at > 1.0 && t < num_active / 2 {
                aca_info_asm_est((aca_cputime() - at) * num_active as f64 / t as f64);
                time_estimated = true;
            }
            if t % (num_active / 12) == 0 {
                aca_info_asm_percent(100 - 100 * (num_active - t) / num_active);
            }
            at0 = aca_realtime();
        }
        let (mut n, mut nd, mut nb, mut bodyid) = (0, 0, 0, 0);
        aca_begin_bulk(t, &mut n, &mut nd, &mut nb, &mut bodyid);
        let (n, nbb) = (n as usize, nb as usize);
        if n != nd as usize && nbb > 0 {
            aca_fatal_bubbles_global();
            return;
        }
        aca_elem_nodes();
        let k = aca_body_material();
        aca_use_material(k);
        let mut g = 0;
        aca_mat_real(1, n as c_int, st.specific_heat.as_mut_ptr(), &mut g);
        aca_mat_real(2, n as c_int, st.heat_ratio.as_mut_ptr(), &mut g);
        aca_mat_real(3, n as c_int, st.density.as_mut_ptr(), &mut g);
        aca_mat_real(4, n as c_int, st.temperature.as_mut_ptr(), &mut g);
        aca_mat_real(5, n as c_int, st.conductivity.as_mut_ptr(), &mut g);
        aca_mat_real(6, n as c_int, st.viscosity.as_mut_ptr(), &mut g);
        for i in 0..n {
            st.lambda[i] = -2.0 / 3.0 * st.viscosity[i];
        }
        aca_mat_real(7, n as c_int, st.bulk_viscosity.as_mut_ptr(), &mut g);
        if g != 0 {
            for i in 0..n {
                st.lambda[i] = st.bulk_viscosity[i] - 2.0 / 3.0 * st.viscosity[i];
            }
        }
        for i in 0..n {
            st.pressure[i] = (st.heat_ratio[i] - 1.0) * st.specific_heat[i] * st.density[i]
                * st.temperature[i];
        }
        for v in st.load.iter_mut() {
            *v = 0.0;
        }
        for v in st.heat_source.iter_mut() {
            *v = 0.0;
        }
        let (mut bfk, mut bfgot) = (0, 0);
        aca_body_bfcount(&mut bfk, &mut bfgot);
        if bfk > 0 {
            aca_use_bodyforce(bfk);
            let mut bg = 0;
            aca_bf_real(1, n as c_int, st.heat_source.as_mut_ptr(), &mut bg);
            aca_bf_real(2, n as c_int, st.heat_source.as_mut_ptr().add(nn), &mut bg);
            for ch in 0..6 {
                aca_bf_real(
                    (3 + ch) as c_int,
                    n as c_int,
                    st.load.as_mut_ptr().add(ch * nn),
                    &mut bg,
                );
            }
            let _ = bg;
        }
        let stride = dofs * nn;
        local_matrix(
            &mut st.local_stiff, &mut st.local_force, stride, af, &st.specific_heat,
            &st.heat_ratio, &st.density, &st.temperature, &st.conductivity, &st.viscosity,
            &st.lambda, &st.heat_source, &st.load, bubbles, mini_bubbles, nn, n, dofs,
            nbb, dim, relorder, axisym,
        );
        aca_update_global(
            st.local_stiff.as_ptr(),
            st.local_force.as_ptr(),
            nn as c_int,
            n as c_int,
            dofs as c_int,
        );
        if block_precond {
            aca_block_bulk(
                st.viscosity.as_ptr(), st.density.as_ptr(), st.specific_heat.as_ptr(),
                st.heat_ratio.as_ptr(), st.temperature.as_ptr(), st.conductivity.as_ptr(),
                st.lambda.as_ptr(), st.pressure.as_ptr(), af, nn as c_int, n as c_int,
                dim as c_int, coordsys, vel_assembly,
            );
        }
    }
    aca_default_finish_bulk();
    // Boundary assembly.
    let num_bulk = aca_num_bulk();
    let num_bnd = aca_num_bnd_elems();
    let nbc = aca_model_num_bcs();
    for t in (num_bulk + 1)..=(num_bulk + num_bnd) {
        let (mut n, mut active, mut constraint) = (0, 0, 0);
        aca_begin_bnd(t, &mut n, &mut active, &mut constraint);
        if active == 0 {
            continue;
        }
        let n = n as usize;
        if aca_perm_has_zero(n as c_int) != 0 {
            continue;
        }
        for i in 1..=nbc {
            if aca_bc_tag(i) != constraint {
                continue;
            }
            aca_use_bc(i);
            aca_elem_nodes();
            let mut ok = 0;
            aca_parent_left(&mut ok);
            if ok == 0 {
                aca_parent_right(&mut ok);
                if ok == 0 {
                    aca_fatal_no_parent();
                    return;
                }
            }
            let k = aca_parent_material();
            aca_use_material(k);
            let mut g = 0;
            aca_mat_real(1, n as c_int, st.specific_heat.as_mut_ptr(), &mut g);
            aca_mat_real(2, n as c_int, st.heat_ratio.as_mut_ptr(), &mut g);
            aca_mat_real(3, n as c_int, st.density.as_mut_ptr(), &mut g);
            aca_mat_real(4, n as c_int, st.temperature.as_mut_ptr(), &mut g);
            aca_mat_real(5, n as c_int, st.conductivity.as_mut_ptr(), &mut g);
            for ii in 0..n {
                st.pressure[ii] = (st.heat_ratio[ii] - 1.0) * st.specific_heat[ii]
                    * st.density[ii] * st.temperature[ii];
            }
            let (mut aiv, mut aig) = (0, 0);
            aca_bc_logical(1, &mut aiv, &mut aig);
            let acoustic_interface = aiv != 0;
            let (mut iuv, mut iug) = (0, 0);
            aca_bc_logical(2, &mut iuv, &mut iug);
            let imp_iter = if iug != 0 { iuv != 0 } else { true };
            let (mut pfv, mut pfg) = (0, 0);
            aca_bc_logical(3, &mut pfv, &mut pfg);
            let potential_flow = pfv != 0;
            let _ = (aig, pfg);
            if acoustic_interface && !st.first_visit {
                if imp_iter {
                    {
                        let (flow, perm) = flow_vals_mut();
                        for j in 0..n {
                            let node = aca_nidx_val((j + 1) as c_int) as usize;
                            let kk = perm[node - 1] as usize;
                            let d = dim * 2 + 4;
                            st.wall_velocity[j] = flow[(kk - 1) * d];
                            st.wall_velocity[nn + j] = flow[(kk - 1) * d + 1];
                            st.wall_velocity[2 * nn + j] = flow[(kk - 1) * d + 2];
                            st.wall_velocity[3 * nn + j] = flow[(kk - 1) * d + 3];
                            if dim > 2 {
                                st.wall_velocity[4 * nn + j] = flow[(kk - 1) * d + 4];
                                st.wall_velocity[5 * nn + j] = flow[(kk - 1) * d + 5];
                            } else {
                                st.wall_velocity[4 * nn + j] = 0.0;
                                st.wall_velocity[5 * nn + j] = 0.0;
                            }
                        }
                    }
                    let (mut avr, mut avi) = (0.0, 0.0);
                    aca_avg_velocity(n as c_int, st.wall_velocity.as_ptr(), nn as c_int, &mut avr, &mut avi);
                    if bem_coupling {
                        let an = aca_ain_size();
                        for j in 0..n {
                            let node = aca_nidx_val((j + 1) as c_int);
                            let mut istat = 1;
                            for kk in 1..=an {
                                let (mut anode, mut r1, mut r2) = (0, 0.0, 0.0);
                                aca_ain_at(kk, &mut anode, &mut r1, &mut r2);
                                if anode == node {
                                    st.load[j] = r1;
                                    st.load[nn + j] = r2;
                                    istat = 0;
                                    break;
                                }
                            }
                            if istat != 0 {
                                aca_fatal_helmholtz_node();
                                return;
                            }
                        }
                        let s1 = aca_sumv(n as c_int, st.load.as_ptr()) / n as f64;
                        let s2 = aca_sumv(n as c_int, st.load.as_ptr().add(nn)) / n as f64;
                        let (zr, zi) = cdiv((-s1, -s2), (avr, avi));
                        for j in 0..n {
                            st.impedance[j] = zr;
                            st.impedance[nn + j] = zi;
                            let (z2r, z2i) = cdiv((0.0, af * st.density[j]), (zr, zi));
                            st.impedance[2 * nn + j] = z2r;
                            st.impedance[3 * nn + j] = z2i;
                        }
                    } else {
                        let mut all = true;
                        for j in 0..n {
                            let node = aca_nidx_val((j + 1) as c_int);
                            if aca_hp_at(node) <= 0 {
                                all = false;
                                break;
                            }
                        }
                        if all {
                            for j in 0..n {
                                let node = aca_nidx_val((j + 1) as c_int);
                                let kk = aca_hp_at(node);
                                let (mut hre, mut him) = (0.0, 0.0);
                                aca_h_at(kk, &mut hre, &mut him);
                                let (zr, zi) = cdiv((-hre, -him), (avr, avi));
                                st.impedance[j] = zr;
                                st.impedance[nn + j] = zi;
                                let (z2r, z2i) =
                                    cdiv((0.0, af * st.density[j]), (zr, zi));
                                st.impedance[2 * nn + j] = z2r;
                                st.impedance[3 * nn + j] = z2i;
                            }
                        }
                    }
                } else {
                    for j in 0..n {
                        st.impedance[j] = 0.0;
                        st.impedance[nn + j] = 0.0;
                        st.impedance[2 * nn + j] = 0.0;
                        st.impedance[3 * nn + j] = 0.0;
                    }
                }
            } else {
                let mut g2 = 0;
                aca_bc_real(1, n as c_int, st.impedance.as_mut_ptr(), &mut g2);
                aca_bc_real(2, n as c_int, st.impedance.as_mut_ptr().add(nn), &mut g2);
                aca_bc_real(3, n as c_int, st.impedance.as_mut_ptr().add(2 * nn), &mut g2);
                aca_bc_real(4, n as c_int, st.impedance.as_mut_ptr().add(3 * nn), &mut g2);
                let _ = g2;
            }
            let mut g3 = 0;
            for ch in 0..6 {
                aca_bc_real(
                    (5 + ch) as c_int,
                    n as c_int,
                    st.load.as_mut_ptr().add(ch * nn),
                    &mut g3,
                );
            }
            let _ = g3;
            let stride = dofs * nn;
            local_matrix_boundary(
                &mut st.local_stiff, &mut st.local_force, stride, af, &st.specific_heat,
                &st.heat_ratio, &st.density, &st.pressure, &st.temperature, &st.conductivity,
                &st.impedance, &st.load, nn, n, dofs, dim, axisym,
            );
            aca_update_global(
                st.local_stiff.as_ptr(),
                st.local_force.as_ptr(),
                nn as c_int,
                n as c_int,
                dofs as c_int,
            );
            if (acoustic_interface && !st.first_visit && !imp_iter) || potential_flow {
                if bem_coupling {
                    let an = aca_ain_size();
                    for j in 0..n {
                        let node = aca_nidx_val((j + 1) as c_int);
                        let mut istat = 1;
                        for kk in 1..=an {
                            let (mut anode, mut r1, mut r2) = (0, 0.0, 0.0);
                            aca_ain_at(kk, &mut anode, &mut r1, &mut r2);
                            if anode == node {
                                st.load[j] = r1;
                                st.load[nn + j] = r2;
                                istat = 0;
                                break;
                            }
                        }
                        if istat != 0 {
                            aca_fatal_helmholtz_node();
                            return;
                        }
                    }
                    aca_local_interface(
                        st.local_stiff.as_mut_ptr(), st.local_force.as_mut_ptr(), af,
                        st.specific_heat.as_ptr(), st.heat_ratio.as_ptr(), st.density.as_ptr(),
                        st.pressure.as_ptr(), st.temperature.as_ptr(), st.conductivity.as_ptr(),
                        st.impedance.as_ptr(), st.load.as_ptr(), nn as c_int, n as c_int,
                        dofs as c_int,
                    );
                    aca_update_global(
                        st.local_stiff.as_ptr(),
                        st.local_force.as_ptr(),
                        nn as c_int,
                        n as c_int,
                        dofs as c_int,
                    );
                    for j in 0..n {
                        let node = aca_nidx_val((j + 1) as c_int) as usize;
                        let (flow, perm) = flow_vals_mut();
                        let m = perm[node - 1] as c_int;
                        let d = (dim * 2 + 4) as c_int;
                        aca_zerorow((m - 1) * d + dim as c_int * 2 + 1);
                        aca_setmat((m - 1) * d + dim as c_int * 2 + 1, (m - 1) * d + dim as c_int * 2 + 1, 1.0);
                        aca_zerorow((m - 1) * d + dim as c_int * 2 + 2);
                        aca_setmat((m - 1) * d + dim as c_int * 2 + 2, (m - 1) * d + dim as c_int * 2 + 2, 1.0);
                        let _ = flow;
                        aca_rhs_set(
                            (m - 1) * d + dim as c_int * 2 + 1,
                            (st.heat_ratio[j] - 1.0) / (st.density[j] * st.heat_ratio[j] * af)
                                * st.load[j],
                        );
                        aca_rhs_set(
                            (m - 1) * d + dim as c_int * 2 + 2,
                            (st.heat_ratio[j] - 1.0) / (st.density[j] * st.heat_ratio[j] * af)
                                * st.load[nn + j],
                        );
                    }
                } else {
                    let mut all = true;
                    for j in 0..n {
                        let node = aca_nidx_val((j + 1) as c_int);
                        if aca_hp_at(node) <= 0 {
                            all = false;
                            break;
                        }
                    }
                    if all {
                        for j in 0..n {
                            let node = aca_nidx_val((j + 1) as c_int);
                            let kk = aca_hp_at(node);
                            let (mut hre, mut him) = (0.0, 0.0);
                            aca_h_at(kk, &mut hre, &mut him);
                            st.load[j] = hre;
                            st.load[nn + j] = him;
                        }
                        aca_local_interface(
                            st.local_stiff.as_mut_ptr(), st.local_force.as_mut_ptr(), af,
                            st.specific_heat.as_ptr(), st.heat_ratio.as_ptr(), st.density.as_ptr(),
                            st.pressure.as_ptr(), st.temperature.as_ptr(), st.conductivity.as_ptr(),
                            st.impedance.as_ptr(), st.load.as_ptr(), nn as c_int, n as c_int,
                            dofs as c_int,
                        );
                        aca_update_global(
                            st.local_stiff.as_ptr(),
                            st.local_force.as_ptr(),
                            nn as c_int,
                            n as c_int,
                            dofs as c_int,
                        );
                        for j in 0..n {
                            let node = aca_nidx_val((j + 1) as c_int) as usize;
                            let (flow, perm) = flow_vals_mut();
                            let m = perm[node - 1] as c_int;
                            let d = (dim * 2 + 4) as c_int;
                            aca_zerorow((m - 1) * d + dim as c_int * 2 + 1);
                            aca_setmat(
                                (m - 1) * d + dim as c_int * 2 + 1,
                                (m - 1) * d + dim as c_int * 2 + 1,
                                1.0,
                            );
                            aca_zerorow((m - 1) * d + dim as c_int * 2 + 2);
                            aca_setmat(
                                (m - 1) * d + dim as c_int * 2 + 2,
                                (m - 1) * d + dim as c_int * 2 + 2,
                                1.0,
                            );
                            let _ = flow;
                            let kk = aca_hp_at(node as c_int);
                            let (mut hre, mut him) = (0.0, 0.0);
                            aca_h_at(kk, &mut hre, &mut him);
                            aca_rhs_set(
                                (m - 1) * d + dim as c_int * 2 + 1,
                                (st.heat_ratio[j] - 1.0) / (st.density[j] * st.heat_ratio[j] * af)
                                    * hre,
                            );
                            aca_rhs_set(
                                (m - 1) * d + dim as c_int * 2 + 2,
                                (st.heat_ratio[j] - 1.0) / (st.density[j] * st.heat_ratio[j] * af)
                                    * him,
                            );
                        }
                    } else {
                        aca_fatal_helmholtz_bnd();
                        return;
                    }
                }
            }
            if block_precond {
                let mut nonzero = false;
                for ch in 0..4 {
                    for j in 0..n {
                        if st.impedance[ch * nn + j] != 0.0 {
                            nonzero = true;
                            break;
                        }
                    }
                }
                if nonzero {
                    aca_block_bnd_imp(
                        st.impedance.as_ptr(), st.specific_heat.as_ptr(), st.heat_ratio.as_ptr(),
                        st.density.as_ptr(), st.conductivity.as_ptr(), af, nn as c_int,
                        n as c_int, dim as c_int, coordsys, vel_assembly,
                    );
                }
            }
        }
    }
    // Slip boundary conditions.
    for t in (num_bulk + 1)..=(num_bulk + num_bnd) {
        let (mut n, mut active, mut constraint) = (0, 0, 0);
        aca_begin_bnd(t, &mut n, &mut active, &mut constraint);
        if active == 0 {
            continue;
        }
        let n = n as usize;
        if aca_perm_has_zero(n as c_int) != 0 {
            continue;
        }
        for i in 1..=nbc {
            if aca_bc_tag(i) != constraint {
                continue;
            }
            aca_use_bc(i);
            let (mut slv, mut slg) = (0, 0);
            aca_bc_logical(4, &mut slv, &mut slg);
            if slv == 0 {
                continue;
            }
            aca_elem_nodes();
            let mut ok = 0;
            aca_parent_left(&mut ok);
            if ok == 0 {
                aca_parent_right(&mut ok);
                if ok == 0 {
                    aca_fatal_no_parent();
                    return;
                }
            }
            let k = aca_parent_material();
            aca_use_material(k);
            let mut g = 0;
            aca_mat_real(1, n as c_int, st.specific_heat.as_mut_ptr(), &mut g);
            aca_mat_real(2, n as c_int, st.heat_ratio.as_mut_ptr(), &mut g);
            aca_mat_real(3, n as c_int, st.density.as_mut_ptr(), &mut g);
            aca_mat_real(5, n as c_int, st.conductivity.as_mut_ptr(), &mut g);
            aca_mat_real(4, n as c_int, st.temperature.as_mut_ptr(), &mut g);
            for ii in 0..n {
                st.pressure[ii] = (st.heat_ratio[ii] - 1.0) * st.specific_heat[ii]
                    * st.density[ii] * st.temperature[ii];
            }
            aca_bc_real(11, n as c_int, st.wall_temperature.as_mut_ptr(), &mut g);
            for ch in 0..6 {
                aca_bc_real(
                    (12 + ch) as c_int,
                    n as c_int,
                    st.wall_velocity.as_mut_ptr().add(ch * nn),
                    &mut g,
                );
            }
            let (mut c1, mut g1) = (0.0, 0);
            aca_bc_constreal(4, &mut c1, &mut g1);
            let (mut c2, mut g2) = (0.0, 0);
            aca_bc_constreal(5, &mut c2, &mut g2);
            let (mut c3, mut g3) = (0.0, 0);
            aca_bc_constreal(6, &mut c3, &mut g3);
            let _ = (g1, g2, g3);
            aca_slip_matrix(
                st.local_stiff.as_mut_ptr(), st.local_force.as_mut_ptr(), af,
                st.specific_heat.as_ptr(), st.heat_ratio.as_ptr(), st.density.as_ptr(),
                st.conductivity.as_ptr(), st.pressure.as_ptr(), st.temperature.as_ptr(),
                st.wall_temperature.as_ptr(), st.wall_velocity.as_ptr(), c1, c2, c3,
                nn as c_int, n as c_int, dofs as c_int,
            );
            aca_update_global(
                st.local_stiff.as_ptr(),
                st.local_force.as_ptr(),
                nn as c_int,
                n as c_int,
                dofs as c_int,
            );
            if block_precond {
                aca_block_bnd_slip(
                    st.specific_heat.as_ptr(), st.heat_ratio.as_ptr(), st.density.as_ptr(),
                    st.temperature.as_ptr(), af, st.conductivity.as_ptr(),
                    st.wall_temperature.as_ptr(), c1, c2, nn as c_int, n as c_int,
                    dim as c_int, coordsys, vel_assembly,
                );
            }
        }
    }
    // Dirichlet BCs.
    for comp in 1..=dim as c_int {
        aca_setdir_velo(comp, dim as c_int);
    }
    aca_setdir_temp(0, dofs as c_int);
    aca_setdir_temp(1, dofs as c_int);
    aca_shell_interface(
        st.gap_indexes.as_mut_ptr(),
        dim as c_int,
        dofs as c_int,
        af,
        nn as c_int,
    );
    aca_default_dirichlet();
    aca_info_done();
    // Preconditioner Dirichlet BCs + traction boundaries (block only).
    if block_precond {
        aca_block_setbcs(dim as c_int, vel_assembly);
        aca_traction_precond(
            dim as c_int,
            af,
            st.specific_heat.as_mut_ptr(), st.heat_ratio.as_mut_ptr(),
            st.density.as_mut_ptr(), st.temperature.as_mut_ptr(),
            st.conductivity.as_mut_ptr(), st.viscosity.as_mut_ptr(),
            st.lambda.as_mut_ptr(), st.bulk_viscosity.as_mut_ptr(),
            st.pressure.as_mut_ptr(), nn as c_int,
        );
    }
    // Solve.
    let asm_t = aca_cputime() - at;
    let st0 = aca_cputime();
    let mut norm = 0.0;
    if block_precond {
        let mut jdeg = 4;
        let (mut jv, mut jg) = (0, 0);
        aca_solver_int(4, &mut jv, &mut jg);
        if jg != 0 {
            jdeg = jv;
        }
        aca_block_solve(
            outer_which,
            max_iterations,
            tolerance,
            dim as c_int,
            jdeg,
            if utilize_prev { 1 } else { 0 },
            &mut norm,
        );
        aca_var_norm_set(norm);
    } else {
        norm = aca_default_solve();
    }
    let solve_t = aca_cputime() - st0;
    let rel = if st.prev_norm + norm != 0.0 {
        2.0 * (st.prev_norm - norm).abs() / (st.prev_norm + norm)
    } else {
        0.0
    };
    aca_info_blank();
    aca_info_result(norm);
    aca_info_relchange(rel);
    aca_info_asmtime(asm_t);
    aca_info_soltime(solve_t);
    st.prev_norm = norm;
    // Post-solve rescaling to true pressure and temperature.
    {
        let nbulk = aca_num_bulk() as c_int;
        let mut visited = vec![false; aca_model_num_nodes() as usize];
        for t in 1..=nbulk {
            let (mut n, mut bodyid) = (0, 0);
            aca_begin_direct(t, &mut n, &mut bodyid);
            if aca_check_eqn() == 0 {
                continue;
            }
            let n = n as usize;
            let k = aca_body_material();
            aca_use_material(k);
            let mut g = 0;
            aca_mat_real(1, n as c_int, st.specific_heat.as_mut_ptr(), &mut g);
            aca_mat_real(2, n as c_int, st.heat_ratio.as_mut_ptr(), &mut g);
            aca_mat_real(4, n as c_int, st.temperature.as_mut_ptr(), &mut g);
            aca_mat_real(3, n as c_int, st.density.as_mut_ptr(), &mut g);
            aca_mat_real(6, n as c_int, st.viscosity.as_mut_ptr(), &mut g);
            for i in 0..n {
                st.lambda[i] = -2.0 / 3.0 * st.viscosity[i];
            }
            aca_mat_real(7, n as c_int, st.bulk_viscosity.as_mut_ptr(), &mut g);
            if g != 0 {
                for i in 0..n {
                    st.lambda[i] = st.bulk_viscosity[i] - 2.0 / 3.0 * st.viscosity[i];
                }
            }
            for i in 0..n {
                st.pressure[i] = (st.heat_ratio[i] - 1.0) * st.specific_heat[i]
                    * st.density[i] * st.temperature[i];
            }
            let (flow, perm) = flow_vals_mut();
            for i in 0..n {
                let node = aca_nidx_val((i + 1) as c_int) as usize;
                let j = perm[node - 1] as usize;
                if j == 0 || visited[j - 1] {
                    continue;
                }
                visited[j - 1] = true;
                let b = (j - 1) * dofs + vel_dofs;
                let (mut g1, mut g2, mut g3, mut g4) = (0.0, 0.0, 0.0, 0.0);
                aca_post_rescale(
                    st.density[i], st.temperature[i], st.lambda[i], st.pressure[i], af,
                    flow[b], flow[b + 1], flow[b + 2], flow[b + 3], &mut g1, &mut g2,
                    &mut g3, &mut g4,
                );
                flow[b] = g1;
                flow[b + 1] = g2;
                flow[b + 2] = g3;
                flow[b + 3] = g4;
            }
        }
    }
    if bem_coupling {
        aca_bem_write(
            if st.first_visit { 1 } else { 0 },
            af,
            dim as c_int,
            st.wall_velocity.as_mut_ptr(),
            st.density.as_mut_ptr(),
            st.load.as_mut_ptr(),
            nn as c_int,
        );
    }
    aca_fluidic_all(
        dim as c_int,
        st.viscosity.as_mut_ptr(),
        st.lambda.as_mut_ptr(),
        st.bulk_viscosity.as_mut_ptr(),
        st.wall_velocity.as_mut_ptr(),
        st.load.as_mut_ptr(),
        st.total_force.as_mut_ptr(),
        st.total_moment.as_mut_ptr(),
        st.total_area.as_mut_ptr(),
        st.moment_output.as_mut_ptr(),
        nn as c_int,
        nbc as c_int,
    );
    aca_list_add_freq(af / (2.0 * std::f64::consts::PI));
    aca_impedance_pass(dim as c_int);
    if bem_coupling {
        aca_bem_post(rel, st.solver_calls);
    }
    st.first_visit = false;
    aca_info_exiting();
    let _ = (norm0, h_found);
}
