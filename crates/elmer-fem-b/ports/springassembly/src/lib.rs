//! Port of `fem/src/modules/SpringAssembly.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/SpringAssembly.F90` — one F77-global entry
//! `springassembler_` with no internal procedures. The emitted dynamic
//! symbol is exactly `springassembler_` (proven by `nm -D` on the reference
//! `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses; `dt` and
//! `TransientSimulation` cross by reference and are read as C scalars
//! (`dt` is unused by the body, exactly as in the original). Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/springassembly_acc.F90`, compiled against the built `.mod`
//! files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Whole-array fills (`K = 0`, `M = 0`,
//! `Work = 0`, diagonal copies) are plain element copies in Rust, in the
//! original's order, so bit-identity holds wherever the call sequence
//! matches.
//!
//! Faithfulness notes (mirrors SpringAssembly.F90 statement for statement):
//! - No `SAVE` state anywhere: the original keeps everything (`Element`,
//!   `K`, `M`, `D`, `f`, `Work`, `MassVals`, flags) in call-local storage
//!   and deallocates at the end. The wrapper-side `acc_*` pointers are
//!   re-pointed on every call/BC with identical lifetimes to the
//!   original's locals.
//! - `Found` is reused across the `GetString`/`GetLogical` calls exactly
//!   as in the original, and all three `HarmonicAssembly` disjuncts are
//!   always evaluated (Fortran `.OR.` does not short-circuit).
//! - The scalar-`Spring` probe runs before the per-component loop and its
//!   values are discarded by the `Work = 0` fill, exactly as in the
//!   original; the two `Warn` calls fire only on that path.
//! - The per-component loop interleaves the `Spring` and `Mass` row reads
//!   per component, exactly as in the original.
//! - Allocation failure calls the same `Fatal('SpringAssembler', ...)` via
//!   `try_reserve`, instead of aborting.
//! - `Element % NodeIndexes => NodeIndexes(i:i)` section association and
//!   the `Perm(...) <= 0` skip live in the wrappers, preserving order.
//! - No Rust threads: the solver entries run serially, and there is no
//!   Rust-side shared state at all.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn spa_begin(h_model: *const c_void, h_solver: *const c_void);
    fn spa_prologue(dofs: *mut c_int, harmonic: *mut c_int);
    fn spa_elem_init();
    fn spa_num_bcs() -> c_int;
    fn spa_bc_begin(bc: c_int) -> c_int;
    fn spa_bc_numnodes() -> c_int;
    fn spa_spring_scalar(n: c_int, ld: c_int, out: *mut c_double);
    fn spa_comp_row(
        which: c_int,
        comp: c_int,
        n: c_int,
        ld: c_int,
        out: *mut c_double,
        found: *mut c_int,
    );
    fn spa_set_elem_node(i: c_int);
    fn spa_node_perm(i: c_int) -> c_int;
    fn spa_update_equations(dd: c_int, k: *mut c_double, f: *mut c_double);
    fn spa_2ndordertime(
        dd: c_int,
        m: *mut c_double,
        d: *mut c_double,
        k: *mut c_double,
        f: *mut c_double,
    );
    fn spa_update_mass(dd: c_int, m: *mut c_double);
    fn spa_fatal_alloc();
}

/// Grow-or-create the `(dofs, n)` work arrays under the original's rule:
/// allocate when absent, reallocate when the stored width is narrower than
/// `need`, otherwise reuse the leading section.
fn ensure_work(
    work: &mut Vec<c_double>,
    massvals: &mut Vec<c_double>,
    cols: &mut usize,
    dofs: usize,
    need: usize,
) {
    if *cols == 0 {
        if work.try_reserve(dofs * need).is_err()
            || massvals.try_reserve(dofs * need).is_err()
        {
            unsafe { spa_fatal_alloc() };
            std::process::abort();
        }
        work.resize(dofs * need, 0.0);
        massvals.resize(dofs * need, 0.0);
        *cols = need;
    } else if *cols < need {
        let mut new_work: Vec<c_double> = Vec::new();
        let mut new_mass: Vec<c_double> = Vec::new();
        if new_work.try_reserve(dofs * need).is_err()
            || new_mass.try_reserve(dofs * need).is_err()
        {
            unsafe { spa_fatal_alloc() };
            std::process::abort();
        }
        new_work.resize(dofs * need, 0.0);
        new_mass.resize(dofs * need, 0.0);
        *work = new_work;
        *massvals = new_mass;
        *cols = need;
    }
}

/// `springassembler_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` is unused by the body,
/// exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (`dt` is ignored).
#[export_name = "springassembler_"]
pub unsafe extern "C" fn spring_assembler_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    transient: *const c_int,
) {
    let transient_sim = !transient.is_null() && unsafe { *transient } != 0;

    unsafe { spa_begin(model, solver) };

    // SolverParams lookup, 'Displacement Variable Name' search (+ Fatal),
    // HarmonicAssembly disjunction, DOFs read — statement order kept.
    let (mut dofs, mut harmonic) = (0, 0);
    unsafe { spa_prologue(&mut dofs, &mut harmonic) };
    let dd = dofs as usize;

    // ALLOCATE(K(DOFs,DOFs), M(DOFs,DOFs), D(DOFs,DOFs), f(DOFs),
    //          Assemble(DOFs), AssembleMass(DOFs)) — plain Vecs, no SAVE.
    let mut k: Vec<c_double> = Vec::new();
    let mut m: Vec<c_double> = Vec::new();
    let mut d: Vec<c_double> = Vec::new();
    let mut f: Vec<c_double> = Vec::new();
    if k.try_reserve(dd * dd).is_err()
        || m.try_reserve(dd * dd).is_err()
        || d.try_reserve(dd * dd).is_err()
        || f.try_reserve(dd).is_err()
    {
        unsafe { spa_fatal_alloc() };
        std::process::abort();
    }
    k.resize(dd * dd, 0.0);
    m.resize(dd * dd, 0.0);
    d.resize(dd * dd, 0.0);
    f.resize(dd, 0.0);
    let mut assemble = vec![false; dd];
    let mut assemble_mass = vec![false; dd];

    // Point element serving the standard assembly routine.
    unsafe { spa_elem_init() };

    // f = 0, D = 0 ahead of the BC loop, exactly as in the original.
    f.fill(0.0);
    d.fill(0.0);

    // Per-BC scratch, mirroring the Work/MassVals allocate-once/grow rule.
    let mut work: Vec<c_double> = Vec::new();
    let mut massvals: Vec<c_double> = Vec::new();
    let mut work_cols: usize = 0;

    let nbcs = unsafe { spa_num_bcs() };
    for bc in 1..=nbcs {
        // ValueList => Model % BCs(BC) % Values; CYCLE unless Target Nodes.
        if unsafe { spa_bc_begin(bc) } == 0 {
            continue;
        }
        let n = unsafe { spa_bc_numnodes() };
        ensure_work(
            &mut work,
            &mut massvals,
            &mut work_cols,
            dd,
            n as usize,
        );

        // Scalar-'Spring' probe (values discarded by the zero fill below).
        unsafe { spa_spring_scalar(n, dofs, work.as_mut_ptr()) };

        work.fill(0.0);
        massvals.fill(0.0);
        // DO i=1,DOFs: interleaved Spring/Mass row reads per component.
        for (j, (slot, mslot)) in
            assemble.iter_mut().zip(assemble_mass.iter_mut()).enumerate()
        {
            let (mut fs, mut fm) = (0, 0);
            unsafe {
                spa_comp_row(
                    1,
                    (j + 1) as c_int,
                    n,
                    dofs,
                    work.as_mut_ptr(),
                    &mut fs,
                );
                spa_comp_row(
                    2,
                    (j + 1) as c_int,
                    n,
                    dofs,
                    massvals.as_mut_ptr(),
                    &mut fm,
                );
            }
            *slot = fs != 0;
            *mslot = fm != 0;
        }

        let spring_assembly = assemble.iter().any(|&a| a);
        let mass_assembly = assemble_mass.iter().any(|&a| a);
        if !spring_assembly && !mass_assembly {
            continue;
        }

        for i in 1..=n {
            unsafe { spa_set_elem_node(i) };
            if unsafe { spa_node_perm(i) } <= 0 {
                continue;
            }
            k.fill(0.0);
            if spring_assembly {
                for (j, &a) in assemble.iter().enumerate() {
                    if !a {
                        continue;
                    }
                    // K(j,j) = Work(j,i), column-major diagonals.
                    k[j * dd + j] = work[(i as usize - 1) * dd + j];
                }
            }
            m.fill(0.0);
            if (transient_sim || harmonic != 0) && mass_assembly {
                for (j, &a) in assemble_mass.iter().enumerate() {
                    if !a {
                        continue;
                    }
                    m[j * dd + j] = massvals[(i as usize - 1) * dd + j];
                }
                if transient_sim {
                    unsafe {
                        spa_2ndordertime(
                            dofs,
                            m.as_mut_ptr(),
                            d.as_mut_ptr(),
                            k.as_mut_ptr(),
                            f.as_mut_ptr(),
                        );
                    }
                } else {
                    unsafe { spa_update_mass(dofs, m.as_mut_ptr()) };
                }
            }
            unsafe {
                spa_update_equations(dofs, k.as_mut_ptr(), f.as_mut_ptr())
            };
        }
    }
    // Deallocation is implicit (Vec drops), matching the original's
    // DEALLOCATE tail.
}
