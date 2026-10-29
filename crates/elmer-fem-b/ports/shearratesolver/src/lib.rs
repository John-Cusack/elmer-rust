//! Port of `fem/src/modules/ShearrateSolver.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ShearrateSolver.F90` — two F77-global
//! entries `shearratesolver_` and `shearratesolver_init_` (proven by
//! `nm -D` on the reference `.so`); one `CONTAINS` internal procedure
//! (`BulkAssembly`) whose symbols are local.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/shearratesolver_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive reductions (`SUM`, `MATMUL`-free but `SUM`-heavy
//! assembly, `SQRT(SecondInvariant/2)`), formatted `Info` output, and
//! whole-array updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ShearrateSolver.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `ViscVector` + `AllocationsDone` live as
//!   wrapper-side SAVE state with identical lifetimes (allocated once at
//!   first `Calculate Viscosity` call, sized `SIZE(RHS)`, never resized);
//!   the internal procedure's SAVE arrays (`STIFF`, `FORCE`, `FORCE2`,
//!   `Vx/Vy/Vz`, `Basis`, `dBasisdx`, `NodalRho`, `NodalVisc`, `Nodes`)
//!   live as wrapper-side SAVE state sized once at first-mesh
//!   `MaxElementNodes`, never resized, exactly as in the original.
//! - Allocation failure calls the same `Fatal('ShearrateSolver', ...)`
//!   via `try_reserve` (the original has no `STAT=` check; the port adds
//!   the guard at the Rust boundary with the solver's own `Fatal`).
//! - `Model`, `dt`, `TransientSimulation` are unused by the bodies,
//!   exactly as in the original (received and ignored).
//! - The `DOFs /= 1` check calls the same
//!   `Fatal('ShearrateSolver','Shearrate should have just 1 component')`.
//! - `BulkAssembly` is transcribed verbatim into the
//!   `sra_bulk_assembly` wrapper (same compiler, same evaluation order,
//!   same whole-array `DefaultUpdateForce(FORCE2)` shape) because its
//!   string-keyed lookups (`ComponentName`, `VariableGet`) and
//!   `MaterialModels` calls cannot cross the C boundary; the Rust entry
//!   drives every branch around it in original order.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sra_init_begin(h_solver: *const c_void);
    fn sra_dim() -> c_int;
    fn sra_listcheck_variable() -> c_int;
    fn sra_add_variable_shearrate();
    fn sra_get_calcvisc() -> c_int;
    fn sra_add_exported_viscosity();
    fn sra_init_tail();
    fn sra_begin(h_solver: *const c_void);
    fn sra_info_banner();
    fn sra_assoc_matrix() -> c_int;
    fn sra_active_perm_count() -> c_int;
    fn sra_var_dofs() -> c_int;
    fn sra_fatal_dofs();
    fn sra_rhs_size() -> c_int;
    fn sra_visc_alloc_maybe();
    fn sra_visc_zero();
    fn sra_visc_field_bind();
    fn sra_get_constbulk() -> c_int;
    fn sra_bind_force_field(const_bulk: c_int) -> c_int;
    fn sra_default_initialize(use_const_bulk: c_int);
    fn sra_bulk_assembly();
    fn sra_finish_bulk(const_bulk: c_int, in_use: c_int);
    fn sra_default_finish();
    fn sra_info_visc();
    fn sra_info_shear();
    fn sra_swap_to_visc();
    fn sra_swap_back();
    fn sra_default_solve() -> c_double;
    fn sra_fatal_alloc();
}

/// First-mesh SAVE state guard. The bulk arrays themselves live in the
/// Fortran wrapper SAVE (identical lifetimes to the original's internal
/// `SAVE`, which Rust cannot name); this guard only serializes the
/// `static` sharing, exactly as in the T020 Poisson port.
struct ShearState {
    _private: (),
}

impl ShearState {
    const fn new() -> Self {
        Self { _private: () }
    }
}

static STATE: Mutex<ShearState> = Mutex::new(ShearState::new());

/// `shearratesolver_init_` — solver init entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "shearratesolver_init_"]
pub unsafe extern "C" fn shearrate_solver_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { sra_init_begin(solver) };

    // dim = CoordinateSystemDimension() — computed and unused, as in the
    // original (kept as a touch call).
    let _dim: c_int = unsafe { sra_dim() };

    // IF (.NOT. ListCheckPresent(SolverParams,'Variable')) THEN ...
    if unsafe { sra_listcheck_variable() } == 0 {
        unsafe { sra_add_variable_shearrate() };
    }

    // IF (GetLogical(SolverParams,'Calculate Viscosity',GotIt)) THEN ...
    if unsafe { sra_get_calcvisc() } != 0 {
        unsafe { sra_add_exported_viscosity() };
    }

    // Time derivative order + linear system defaults.
    unsafe { sra_init_tail() };
}

/// `shearratesolver_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt`, and the
/// transient flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "shearratesolver_"]
pub unsafe extern "C" fn shearrate_solver_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let _guard = STATE.lock().unwrap();
    unsafe { sra_begin(solver) };

    // CALL Info(...) banner x3.
    unsafe { sra_info_banner() };

    // IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN
    if unsafe { sra_assoc_matrix() } == 0 {
        return;
    }
    // IF (COUNT(Solver % Variable % Perm > 0) <= 0) RETURN
    if unsafe { sra_active_perm_count() } <= 0 {
        return;
    }

    // ShearrateSol => Solver % Variable; DOFs check.
    if unsafe { sra_var_dofs() } != 1 {
        unsafe { sra_fatal_dofs() };
        return;
    }

    // IF (CalculateViscosity) THEN — alloc once, zero, bind ViscSol/Field.
    let calc_visc = unsafe { sra_get_calcvisc() } != 0;
    if calc_visc {
        // try_reserve guard for the wrapper-side ALLOCATE (the original
        // has no STAT= check; failure calls the solver's own Fatal).
        // SIZE(RHS) is read from the solver before allocating.
        let need = unsafe { sra_rhs_size() } as usize;
        let mut probe: Vec<c_double> = Vec::new();
        if probe.try_reserve(need).is_err() {
            unsafe { sra_fatal_alloc() };
            std::process::abort();
        }
        unsafe {
            sra_visc_alloc_maybe();
            sra_visc_zero();
            sra_visc_field_bind();
        }
    }

    // ConstantBulkMatrixInUse = ConstantBulkMatrix .AND.
    // ASSOCIATED(Solver % Matrix % BulkValues); ForceVector/ShearrateField
    // binds. Returns ConstantBulkMatrixInUse.
    let const_bulk = unsafe { sra_get_constbulk() };
    let in_use = unsafe { sra_bind_force_field(const_bulk) };

    unsafe {
        sra_default_initialize(in_use);
        sra_bulk_assembly();
        sra_finish_bulk(const_bulk, in_use);
        sra_default_finish();
    }

    // IF (CalculateViscosity) THEN solve on ViscVector, swap back.
    if calc_visc {
        unsafe { sra_info_visc() };
        unsafe { sra_swap_to_visc() };
        let _unorm: c_double = unsafe { sra_default_solve() };
        unsafe { sra_swap_back() };
    }

    unsafe { sra_info_shear() };

    // UNorm = DefaultSolve() (Norm local and unused, as in the original).
    let _unorm: c_double = unsafe { sra_default_solve() };
}
