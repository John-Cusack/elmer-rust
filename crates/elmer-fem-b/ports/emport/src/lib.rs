//! Port of `fem/src/modules/EMPort.F90` (Elmer @ c0b4bc6, the campaign
//! fork) for lane T022-A, at the solver-plugin boundary (whole `.so` →
//! Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/EMPort.F90` — four F77-global entries
//! (`emportsolver_init0_`, `emportsolver_`, `emportpotential_init_`,
//! `emportpotential_`) with six `CONTAINS` internals (`ParallelActiveSubset`,
//! `InitStuff`, `LocalMatrix`, `LocalMatrix_With_Potentials`,
//! `CalculatePortPower`, `EMPortPost`, plus the potential solver's own
//! `LocalMatrix`). The emitted dynamic symbols are exactly those four
//! (proven by `nm -D` on the reference `.so`); the internals emit only
//! local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/emport_acc.F90`, compiled
//! against the built `.mod` files of the target configure — the port
//! itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`, `MATMUL`), formatted I/O
//! (`WRITE`/`Info`), whole-array complex updates, and the strided
//! `EF % Values(j::dofs)` scatter stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors EMPort.F90 statement for statement):
//! - `SAVE` semantics replicated: the entry-level `ValueHandle_t`
//!   coefficients, the per-internal-procedure `SAVE Nodes`/`IP`/basis
//!   arrays (one set per internal procedure, exactly as in the original —
//!   `CalculatePortPower` keeps a single set shared by both call sites),
//!   and the full `EMPortPost` SAVE set live as wrapper-side SAVE state
//!   with identical lifetimes.
//! - `EigenProblem = EigenOrHarmonicAnalysis(Solver)` assigns a variable
//!   the original never reads afterwards; the call is kept (same call
//!   sequence) and the dead result discarded.
//! - `BetaSum` is reset on every entry (as in the original); the portwise
//!   `maxmu`/`maxeps` reset rides with `DefaultInitialize` exactly where
//!   the original zeroes them.
//! - The multi-port (`MaxPort > 1`) matrix surgery (`FreeMatrix`, pointer
//!   `Null()`, `SavePerm`/`SaveEigenVectors` collect/restore) is
//!   transcribed verbatim in wrappers; the serial proof takes the
//!   `MaxPort <= 1` path, the `MaxPort == 1` single-port filter path is
//!   exercised by proof config B.
//! - The `ParEnv % PEs > 1` branches (`ParallelActiveSubset` MPI calls,
//!   `PartitionVector`/`ParallelCdot` path) are transcribed verbatim in
//!   the wrapper but guarded by the same `PEs == 1` early return as the
//!   original; the serial proof never reaches them.
//! - Allocation failure calls the same `Fatal(..., 'Local storage
//!   allocation failed')` from inside the wrappers (the `STAT=allocstat`
//!   check), exactly as in the original.
//! - `Model`, `dt`, and the transient flag are unused by the solver body
//!   itself beyond handle forwarding, exactly as in the original.
//! - No Rust threads: the solver calls entries serially, so no sharing
//!   guard is needed.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn emi0_begin(h_model: *const c_void, h_solver: *const c_void);
    fn emi0_body();
    fn ema_begin(h_model: *const c_void, h_solver: *const c_void);
    fn ema_start_info();
    fn ema_init_handles();
    fn ema_check_cartesian();
    fn ema_check_dofs();
    fn ema_use_v() -> c_int;
    fn ema_scan_ports() -> c_int;
    fn ema_fetch_emvar();
    fn ema_multiport_setup(nrows: *mut c_int, ncols: *mut c_int) -> c_int;
    fn ema_edge_style();
    fn ema_eigen_query();
    fn ema_default_start();
    fn ema_init_stuff();
    fn ema_fetch_mesh();
    fn ema_active() -> c_int;
    fn ema_set_portind(portind: c_int);
    fn ema_port_begin(portind: c_int) -> c_int;
    fn ema_default_mode() -> c_int;
    fn ema_elem_begin(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        efamily: *mut c_int,
        skip: *mut c_int,
    );
    fn ema_check_edge(n: c_int, efamily: c_int);
    fn ema_lm(n: c_int, nd: c_int);
    fn ema_lm_pot(n: c_int, nd: c_int);
    fn ema_finish_bulk();
    fn ema_finish_assembly();
    fn ema_dirichlet();
    fn ema_shift();
    fn ema_solve() -> c_double;
    fn ema_beta();
    fn ema_beta_to_bcs();
    fn ema_norm_begin() -> c_int;
    fn ema_norm_elem(t: c_int);
    fn ema_norm_finish();
    fn ema_nonhomog_begin() -> c_int;
    fn ema_nonhomog(j: c_int);
    fn ema_nonhomog_end();
    fn ema_impedance_check() -> c_int;
    fn ema_impedance_begin();
    fn ema_impedance_elem(t: c_int);
    fn ema_impedance_finish();
    fn ema_calc_nodal() -> c_int;
    fn ema_post(portind: c_int, maxport: c_int);
    fn ema_port_collect();
    fn ema_parallel_reactivate();
    fn ema_default_finish();
    fn ema_multiport_restore();
    fn ema_set_norm();
    fn ema_done_info();
    fn emp_init_begin(h_model: *const c_void, h_solver: *const c_void);
    fn emp_init_body();
    fn emp_solve_begin(h_model: *const c_void, h_solver: *const c_void);
    fn emp_default_start();
    fn emp_default_initialize();
    fn emp_active() -> c_int;
    fn emp_elem(t: c_int, inith: *mut c_int);
    fn emp_finish_bulk();
    fn emp_finish_boundary();
    fn emp_finish_assembly();
    fn emp_dirichlet();
    fn emp_solve() -> c_double;
    fn emp_finish();
}

fn port_upper(maxport: c_int) -> c_int {
    if maxport > 1 { maxport } else { 1 }
}

/// `emportsolver_init0_` — the `_Init0` setup entry. Mirrors
/// `EMPortSolver_Init0` (complex/eigensystem defaults, primary-solver edge
/// definition copy, element-type default, port BC Dirichlet presets).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "emportsolver_init0_"]
pub unsafe extern "C" fn emportsolver_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emi0_begin(model, solver);
        emi0_body();
    }
}

/// `emportsolver_` — the eigenmode solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// Same contract as [`emportsolver_init0_entry`].
#[export_name = "emportsolver_"]
pub unsafe extern "C" fn emportsolver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        ema_begin(model, solver);
        ema_start_info();
        ema_init_handles();
        ema_check_cartesian();
        ema_check_dofs();
        let usev = ema_use_v();
        let maxport = ema_scan_ports();
        ema_fetch_emvar();
        let mut nrows: c_int = 0;
        let mut ncols: c_int = 0;
        let _multi = ema_multiport_setup(&mut nrows, &mut ncols);
        ema_edge_style();
        ema_eigen_query();
        ema_default_start();
        ema_init_stuff();
        ema_fetch_mesh();
        let active = ema_active();
        let upper = port_upper(maxport);
        for portind in 1..=upper {
            ema_set_portind(portind);
            if ema_port_begin(portind) != 0 {
                continue;
            }
            let _mode = ema_default_mode();
            for t in 1..=active {
                let mut n: c_int = 0;
                let mut nd: c_int = 0;
                let mut efamily: c_int = 0;
                let mut skip: c_int = 0;
                ema_elem_begin(t, &mut n, &mut nd, &mut efamily, &mut skip);
                if skip != 0 {
                    continue;
                }
                ema_check_edge(n, efamily);
                if usev != 0 {
                    ema_lm_pot(n, nd);
                } else {
                    ema_lm(n, nd);
                }
            }
            ema_finish_bulk();
            ema_finish_assembly();
            ema_dirichlet();
            ema_shift();
            let _norm: c_double = ema_solve();
            ema_beta();
            ema_beta_to_bcs();
            if ema_norm_begin() != 0 {
                for t in 1..=active {
                    ema_norm_elem(t);
                }
                ema_norm_finish();
            }
            let nonhomog = ema_nonhomog_begin();
            for j in 1..=nonhomog {
                ema_nonhomog(j);
            }
            if nonhomog > 0 {
                ema_nonhomog_end();
            }
            if ema_impedance_check() != 0 {
                ema_impedance_begin();
                for t in 1..=active {
                    ema_impedance_elem(t);
                }
                ema_impedance_finish();
            }
            if ema_calc_nodal() != 0 {
                ema_post(portind, maxport);
            }
            ema_port_collect();
        }
        ema_parallel_reactivate();
        ema_default_finish();
        ema_multiport_restore();
        ema_set_norm();
        ema_done_info();
    }
}

/// `emportpotential_init_` — the `_Init` setup entry. Mirrors
/// `EMPortPotential_Init` (default `Variable`, port ground/feed Dirichlet
/// presets).
///
/// # Safety
///
/// Same contract as [`emportsolver_init0_entry`].
#[export_name = "emportpotential_init_"]
pub unsafe extern "C" fn emportpotential_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emp_init_begin(model, solver);
        emp_init_body();
    }
}

/// `emportpotential_` — the electrostatics solver entry. Mirrors
/// `EMPortPotential` (bulk assembly with parent-material fallback,
/// `CondensateP`, `DefaultUpdateEquations`, `DefaultSolve`).
///
/// # Safety
///
/// Same contract as [`emportsolver_init0_entry`].
#[export_name = "emportpotential_"]
pub unsafe extern "C" fn emportpotential_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        emp_solve_begin(model, solver);
        emp_default_start();
        emp_default_initialize();
        let active = emp_active();
        let mut inith: c_int = 1;
        for t in 1..=active {
            emp_elem(t, &mut inith);
        }
        emp_finish_bulk();
        emp_finish_boundary();
        emp_finish_assembly();
        emp_dirichlet();
        let _norm: c_double = emp_solve();
        emp_finish();
    }
}
