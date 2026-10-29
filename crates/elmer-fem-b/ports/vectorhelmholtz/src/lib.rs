//! Port of `fem/src/modules/VectorHelmholtz.F90` (Elmer @ 9f6af2f,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/VectorHelmholtz.F90` — seven F77-global
//! entries (`vectorhelmholtzsolver_`, `vectorhelmholtzsolver_init_`,
//! `vectorhelmholtzsolver_init0_`, `vectorhelmholtz_dummy_`,
//! `vectorhelmholtzcalcfields_`, `vectorhelmholtzcalcfields_init_`,
//! `vectorhelmholtzcalcfields_init0_`), proven by `nm -D` on the reference
//! `.so`; the CONTAINS internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/vectorhelmholtz_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset.
//! Order-sensitive complex reductions (`MATMUL`/`SUM`), formatted I/O
//! (`PRINT`/`WRITE`), and whole-array updates stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors VectorHelmholtz.F90 statement for
//! statement):
//! - The element assembly bodies (`LocalMatrix`, `LocalMatrixBC`,
//!   `SingleDipoleLoad`, the CalcFields local/global/face-term workers)
//!   are transcribed verbatim as accessor module procedures: they own
//!   every order-sensitive complex reduction and every formatted write,
//!   with identical statement order and identical branch predicates.
//! - `SAVE` semantics replicated: the `LocalMatrix`/`LocalMatrixBC`
//!   first-mesh sizings (`MaxElementDOFs`), their internal-procedure SAVE
//!   (`Nodes`, `ValueHandle`s, `AllocationsDone`), the dipole `Nodes`,
//!   and the CalcFields face-term `Nodes` live as wrapper-side SAVE state
//!   with identical lifetimes; the solver-prologue scalars and the
//!   CalcFields phase state persist across the phase wrappers of one
//!   entry call exactly like the original entry locals (rewritten by each
//!   prologue, consumed before the next fill).
//! - The `WithNDOFs` latch of `VectorHelmholtzSolver_Init0` travels from
//!   the head wrapper to the tail wrapper, mirroring the single-procedure
//!   lifetime; the cubic-element `Fatal` paths keep the original
//!   add-then-die order.
//! - `Model`, `dt`, `TransientSimulation` are received and ignored,
//!   exactly as in the original (the `Dummy` entry is empty, as in the
//!   original).
//! - No Rust threads: the solver entries run serially, so no `static`
//!   sharing guard is needed — the Rust side holds no cross-invocation
//!   state at all (every loop counter and branch predicate is a call
//!   local); all persistent state lives in the Fortran wrappers.
//! - The `im` unit-imaginary factor and `AEPS` are the
//!   `VectorHelmholtzUtils`/`Types` parameters, exactly as in the
//!   original unit.
//! - Allocation failure calls the same `Fatal(...)` in the wrappers, never
//!   aborting directly.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    // -- VectorHelmholtzSolver_Init0 --
    fn vh_s0_begin(model: *mut c_void, solver: *mut c_void, with_ndofs: *mut c_int);
    fn vh_s0_has_element() -> c_int;
    fn vh_s0_edge_style(k: *mut c_int, piola: *mut c_int, second: *mut c_int);
    fn vh_s0_add_element(which: c_int);
    fn vh_s0_fatal_cubic(second: c_int);
    fn vh_s0_tail();
    // -- VectorHelmholtzSolver_Init --
    fn vh_si_begin(model: *mut c_void, solver: *mut c_void);
    fn vh_si_prec_prefix();
    fn vh_si_find_eigen() -> c_int;
    fn vh_si_emport() -> c_int;
    fn vh_si_note_emport(soln: c_int);
    fn vh_si_fatal_noport();
    fn vh_si_define_ports();
    // -- VectorHelmholtzSolver --
    fn vh_s_begin(model: *mut c_void, solver: *mut c_void);
    fn vh_s_prologue();
    fn vh_s_iter_params(noit: *mut c_int, edgebasis: *mut c_int);
    fn vh_s_default_start();
    fn vh_s_dosolve() -> c_int;
    fn vh_s_refactorize_false();
    fn vh_s_refactorize_remove();
    fn vh_s_default_finish();
    fn vh_s_circ_idx() -> c_int;
    fn vh_s_circulation(idx: c_int);
    fn vh_s_done();
    // -- VectorHelmholtzCalcFields_Init0 / Init --
    fn vcf_init0(model: *mut c_void, solver: *mut c_void);
    fn vcf_init(model: *mut c_void, solver: *mut c_void);
    // -- VectorHelmholtzCalcFields --
    fn vcf_prologue(model: *mut c_void, solver: *mut c_void);
    fn vcf_nactive() -> c_int;
    fn vcf_elem(i: c_int);
    fn vcf_mid();
    fn vcf_elpot();
    fn vcf_finish();
}

/// `vectorhelmholtzsolver_init0_` — element/variable/list setup.
/// Mirrors the `WithNDOFs` / `Element`-present / basis-degree branch
/// structure; the list writes stay in the wrappers.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzsolver_init0_"]
pub unsafe extern "C" fn vectorhelmholtz_solver_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut with_ndofs: c_int = 0;
    vh_s0_begin(model, solver, &mut with_ndofs);
    if vh_s0_has_element() == 0 {
        let (mut k, mut piola, mut second): (c_int, c_int, c_int) = (0, 0, 0);
        vh_s0_edge_style(&mut k, &mut piola, &mut second);
        if with_ndofs != 0 {
            if k == 2 {
                vh_s0_add_element(1);
            } else if second != 0 {
                vh_s0_add_element(2);
            } else if piola != 0 {
                vh_s0_add_element(3);
            } else {
                vh_s0_add_element(4);
            }
        } else {
            match k {
                3 => {
                    if second != 0 {
                        vh_s0_add_element(5);
                    } else {
                        vh_s0_add_element(6);
                    }
                    vh_s0_fatal_cubic(second);
                }
                2 => {
                    if second != 0 {
                        vh_s0_add_element(7);
                    } else {
                        vh_s0_add_element(8);
                    }
                }
                1 => {
                    if second != 0 {
                        vh_s0_add_element(9);
                    } else if piola != 0 {
                        vh_s0_add_element(10);
                    } else {
                        vh_s0_add_element(11);
                    }
                }
                _ => {}
            }
        }
    }
    vh_s0_tail();
}

/// `vectorhelmholtzsolver_init_` — preconditioning-matrix flag, eigenmode
/// port scan, and port-parameter definition.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzsolver_init_"]
pub unsafe extern "C" fn vectorhelmholtz_solver_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    vh_si_begin(model, solver);
    vh_si_prec_prefix();
    if vh_si_find_eigen() != 0 {
        let soln = vh_si_emport();
        if soln == 0 {
            vh_si_fatal_noport();
        } else {
            vh_si_note_emport(soln);
        }
    }
    vh_si_define_ports();
}

/// `vectorhelmholtzsolver_` — the harmonic Maxwell curl-curl solver entry.
/// Argument order and by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt`/transient are
/// received and ignored, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzsolver_"]
pub unsafe extern "C" fn vectorhelmholtz_solver_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    vh_s_begin(model, solver);
    vh_s_prologue();
    let (mut noit, mut edgebasis): (c_int, c_int) = (0, 0);
    vh_s_iter_params(&mut noit, &mut edgebasis);
    vh_s_default_start();
    // DO i=1,NoIterationsMax with EXIT on convergence; the refactorize
    // toggle only applies with the edge-basis path, as in the original.
    for _ in 0..noit {
        if vh_s_dosolve() != 0 {
            break;
        }
        if edgebasis != 0 {
            vh_s_refactorize_false();
        }
    }
    if edgebasis != 0 {
        vh_s_refactorize_remove();
    }
    vh_s_default_finish();
    let circ = vh_s_circ_idx();
    if circ != 0 {
        vh_s_circulation(circ);
    }
    vh_s_done();
}

/// `vectorhelmholtz_dummy_` — empty placeholder for the on-the-fly DG
/// solver created by `VectorHelmholtzCalcFields_Init0`, exactly as in the
/// original.
///
/// # Safety
///
/// Arguments are received and ignored, exactly as in the original.
#[export_name = "vectorhelmholtz_dummy_"]
pub unsafe extern "C" fn vectorhelmholtz_dummy_entry(
    _model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
}

/// `vectorhelmholtzcalcfields_init0_` — primary-solver lookup and
/// on-the-fly DG solver creation (the model surgery stays in the wrapper).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzcalcfields_init0_"]
pub unsafe extern "C" fn vectorhelmholtz_calcfields_init0_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    vcf_init0(model, solver);
}

/// `vectorhelmholtzcalcfields_init_` — variable/output registration.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzcalcfields_init_"]
pub unsafe extern "C" fn vectorhelmholtz_calcfields_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    vcf_init(model, solver);
}

/// `vectorhelmholtzcalcfields_` — derived-field computation over the
/// active elements, then face terms, global solves, and the electric
/// potential branch. The per-element loop bound is evaluated once, as in
/// the original `DO i = 1, GetNOFActive()`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "vectorhelmholtzcalcfields_"]
pub unsafe extern "C" fn vectorhelmholtz_calcfields_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    vcf_prologue(model, solver);
    let nactive = vcf_nactive();
    for i in 1..=nactive {
        vcf_elem(i);
    }
    vcf_mid();
    vcf_elpot();
    vcf_finish();
}
