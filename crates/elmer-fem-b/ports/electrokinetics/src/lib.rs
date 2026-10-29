//! Port of `fem/src/modules/Electrokinetics.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T022-A, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/Electrokinetics.F90` — six F77-global
//! UDF entries (`helmholtz_smoluchowski1/2/3_`,
//! `helmholtz_smoluchowski_`, `helmholtz_smoluchowski_comp_`,
//! `getjouleheat_`), called per node as `Real Procedure` functions with
//! dummy list `(Model, NodeNumber, dummy)`. The emitted dynamic symbols
//! are exactly those six (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, UDF gate shape): the entries' `Model` dummy
//! crosses as an opaque address, `NodeNumber`/result travel by
//! reference/value exactly like the Fortran dummies, and every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/electrokinetics_acc.F90`, compiled against
//! the built `.mod` files of the target configure — the port itself
//! contains no layout assumption and no pinned offset. All
//! floating-point work (`NormalVector`, `SUM` reductions, the mobility
//! quotient, the tangential projection) and all formatted I/O stay in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Faithfulness notes (mirrors Electrokinetics.F90 statement for
//! statement):
//! - The four thin component entries forward to the shared
//!   `helmholtz_smoluchowski_comp` path with directions 1/2/3/0,
//!   exactly as the originals call the internal function.
//! - `SAVE` semantics replicated: `Nodes`/`DIM`/`vacuumPerm`/`FirstTime`
//!   live as wrapper-side SAVE state with identical lifetimes; the
//!   `SkipThisTime` local is re-read from the Simulation list on every
//!   call, as in the original.
//! - `helmholtz_smoluchowski_comp_` overlays its `direction` INTEGER
//!   dummy on the caller's real slot, exactly as the Fortran callee
//!   would (never exercised through a `.sif`; direct calls pass the
//!   direction in the documented slot).
//! - Allocation failure calls the same `Fatal(...)` via the wrapper
//!   (`STAT=istat` check), never aborting directly.
//! - `Model % CurrentElement` is retargeted to the parent element while
//!   reading material parameters and restored before the velocity
//!   assembly, exactly as in the original.
//! - No Rust threads: entries run serially per node, so there are no
//!   Rust-side statics at all (all shared state is wrapper SAVE state,
//!   as in the original).

use std::ffi::c_void;
use std::os::raw::{c_char, c_double, c_int};

extern "C" {
    fn eka_begin(h_model: *mut c_void, node: c_int);
    fn eka_init_phase() -> c_int;
    fn eka_setup(direction: c_int) -> c_int;
    fn eka_find_parent();
    fn eka_bc_normal(normal: *mut c_double) -> c_int;
    fn eka_zeta(zeta: *mut c_double, calc: *mut c_int);
    fn eka_use_parent();

    fn eka_viscosity(visc: *mut c_double, got: *mut c_int);
    fn eka_relperm(dc: *mut c_double, got: *mut c_int);
    fn eka_mob_from_zeta(zeta: c_double, visc: c_double, dc: c_double) -> c_double;
    fn eka_eo_direct(eo: *mut c_double, got: *mut c_int);
    fn eka_efmethod(buf: *mut c_char, buflen: c_int, outlen: *mut c_int, got: *mut c_int);
    fn eka_const_comp(i: c_int, val: *mut c_double, got: *mut c_int);
    fn eka_warn_no_const();
    fn eka_comp_efield(i: c_int, val: *mut c_double, found: *mut c_int);
    fn eka_info_no_efield();
    fn eka_warn_unknown();
    fn eka_restore();
    fn eka_finish(direction: c_int, efield: *const c_double, eomob: c_double) -> c_double;
    fn ekj_deprecated();
    fn ekj_dim() -> c_int;
    fn ekj_elem(body: *mut c_int, eq: *mut c_int, mat: *mut c_int) -> c_int;
    fn ekj_conductivity(cond: *mut c_double, got: *mut c_int);
    fn ekj_density(dens: *mut c_double, got: *mut c_int);
    fn ekj_const_comp(i: c_int, val: *mut c_double, got: *mut c_int);
    fn ekj_warn_no_const();
    fn ekj_info_no_efield();
    fn ekj_warn_unknown();
    fn ekj_finish(efield: *const c_double, cond: c_double, dens: c_double) -> c_double;
}

/// Byte-compare of a `C_CHAR` buffer prefix against an expected method
/// word (the `.sif` parser lowercases unquoted values, so the wrapper
/// bytes compare directly).
fn method_is(buf: &[c_char], n: usize, want: &[u8]) -> bool {
    n == want.len() && buf[..n].iter().zip(want.iter()).all(|(a, b)| *a as u8 == *b)
}

/// Shared `helmholtz_smoluchowski_comp` path: component `direction` of
/// the electroosmotic slip velocity at one boundary node
/// (`direction == 0` selects the tangential value in 2D).
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the loader passes, and
/// `node` must be a valid node number, exactly as for the Fortran UDF.
unsafe fn hs_comp(model: *mut c_void, node: c_int, direction: c_int) -> c_double {
    unsafe {
        eka_begin(model, node);
        if eka_init_phase() != 0 {
            return 0.0;
        }
        let dim = eka_setup(direction);
        eka_find_parent();
        let mut normal = [0.0 as c_double; 3];
        let _bc_id = eka_bc_normal(normal.as_mut_ptr());
        let mut zeta = 0.0 as c_double;
        let mut calc = 0 as c_int;
        eka_zeta(&mut zeta, &mut calc);
        let eomob = if calc != 0 {
            eka_use_parent();
            let mut visc = 0.0 as c_double;
            let mut got = 0 as c_int;
            eka_viscosity(&mut visc, &mut got);
            let mut dc = 0.0 as c_double;
            let mut got_dc = 0 as c_int;
            eka_relperm(&mut dc, &mut got_dc);
            eka_mob_from_zeta(zeta, visc, dc)
        } else {
            let mut eo = 0.0 as c_double;
            let mut got = 0 as c_int;
            eka_eo_direct(&mut eo, &mut got);
            eo
        };
        let mut buf = [0 as c_char; 256];
        let mut outlen = 0 as c_int;
        let mut got_m = 0 as c_int;
        eka_efmethod(buf.as_mut_ptr(), buf.len() as c_int, &mut outlen, &mut got_m);
        let mut efield = [0.0 as c_double; 3];
        let mut exists = false;
        if got_m != 0 {
            let n = outlen as usize;
            if method_is(&buf, n, b"constant") {
                for i in 1..=dim {
                    let mut v = 0.0 as c_double;
                    let mut g = 0 as c_int;
                    eka_const_comp(i, &mut v, &mut g);
                    if g != 0 {
                        efield[(i - 1) as usize] = v;
                        exists = true;
                    }
                }
                if !exists {
                    eka_warn_no_const();
                }
            } else if method_is(&buf, n, b"computed") {
                for i in 1..=dim {
                    let mut v = 0.0 as c_double;
                    let mut f = 0 as c_int;
                    eka_comp_efield(i, &mut v, &mut f);
                    if f != 0 {
                        efield[(i - 1) as usize] = v;
                        exists = true;
                    }
                }
            } else {
                eka_warn_unknown();
            }
        } else {
            eka_info_no_efield();
        }
        if !exists {
            return 0.0;
        }
        eka_restore();
        eka_finish(direction, efield.as_ptr(), eomob)
    }
}

/// Shared `getJouleHeat` path: Joule heat source at one bulk node.
///
/// # Safety
///
/// Same contract as [`hs_comp`].
unsafe fn joule_heat(model: *mut c_void, node: c_int) -> c_double {
    unsafe {
        eka_begin(model, node);
        ekj_deprecated();
        let dim = ekj_dim();
        let mut _body = 0 as c_int;
        let mut _eq = 0 as c_int;
        let mut _mat = 0 as c_int;
        if ekj_elem(&mut _body, &mut _eq, &mut _mat) != 0 {
            return 0.0;
        }
        let mut cond = 0.0 as c_double;
        let mut got_c = 0 as c_int;
        ekj_conductivity(&mut cond, &mut got_c);
        if got_c == 0 {
            return 0.0;
        }
        let mut dens = 0.0 as c_double;
        let mut got_d = 0 as c_int;
        ekj_density(&mut dens, &mut got_d);
        let mut buf = [0 as c_char; 256];
        let mut outlen = 0 as c_int;
        let mut got_m = 0 as c_int;
        eka_efmethod(buf.as_mut_ptr(), buf.len() as c_int, &mut outlen, &mut got_m);
        let mut efield = [0.0 as c_double; 3];
        let mut exists = false;
        if got_m != 0 {
            let n = outlen as usize;
            if method_is(&buf, n, b"constant") {
                for i in 1..=dim {
                    let mut v = 0.0 as c_double;
                    let mut g = 0 as c_int;
                    ekj_const_comp(i, &mut v, &mut g);
                    if g != 0 {
                        efield[(i - 1) as usize] = v;
                        exists = true;
                    }
                }
                if !exists {
                    ekj_warn_no_const();
                }
            } else if method_is(&buf, n, b"computed") {
                for i in 1..=dim {
                    let mut v = 0.0 as c_double;
                    let mut f = 0 as c_int;
                    eka_comp_efield(i, &mut v, &mut f);
                    if f != 0 {
                        efield[(i - 1) as usize] = v;
                        exists = true;
                    }
                }
            } else {
                ekj_warn_unknown();
            }
        } else {
            ekj_info_no_efield();
        }
        if !exists {
            return 0.0;
        }
        ekj_finish(efield.as_ptr(), cond, dens)
    }
}

/// `helmholtz_smoluchowski1_` — x-component of the HS slip velocity.
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the loader passes;
/// `node`/`dummy` travel exactly like the Fortran dummy list.
#[export_name = "helmholtz_smoluchowski1_"]
pub unsafe extern "C" fn helmholtz_smoluchowski1_entry(
    model: *mut c_void,
    node: *const c_int,
    _dummy: *const c_double,
) -> c_double {
    unsafe { hs_comp(model, *node, 1) }
}

/// `helmholtz_smoluchowski2_` — y-component of the HS slip velocity.
///
/// # Safety
///
/// Same contract as [`helmholtz_smoluchowski1_entry`].
#[export_name = "helmholtz_smoluchowski2_"]
pub unsafe extern "C" fn helmholtz_smoluchowski2_entry(
    model: *mut c_void,
    node: *const c_int,
    _dummy: *const c_double,
) -> c_double {
    unsafe { hs_comp(model, *node, 2) }
}

/// `helmholtz_smoluchowski3_` — z-component of the HS slip velocity.
///
/// # Safety
///
/// Same contract as [`helmholtz_smoluchowski1_entry`].
#[export_name = "helmholtz_smoluchowski3_"]
pub unsafe extern "C" fn helmholtz_smoluchowski3_entry(
    model: *mut c_void,
    node: *const c_int,
    _dummy: *const c_double,
) -> c_double {
    unsafe { hs_comp(model, *node, 3) }
}

/// `helmholtz_smoluchowski_` — tangential HS slip velocity (2D only).
///
/// # Safety
///
/// Same contract as [`helmholtz_smoluchowski1_entry`].
#[export_name = "helmholtz_smoluchowski_"]
pub unsafe extern "C" fn helmholtz_smoluchowski_entry(
    model: *mut c_void,
    node: *const c_int,
    _dummy: *const c_double,
) -> c_double {
    unsafe { hs_comp(model, *node, 0) }
}

/// `helmholtz_smoluchowski_comp_` — component `direction` of the HS slip
/// velocity. The `direction` dummy is an INTEGER overlaid on the
/// caller's real slot, exactly as the Fortran callee observes it.
///
/// # Safety
///
/// `model`/`node` follow the UDF contract; `dir_slot` must be a valid
/// read of at least `size_of::<c_int>()` bytes, as passed by the Elmer
/// procedure-call trampoline.
#[export_name = "helmholtz_smoluchowski_comp_"]
pub unsafe extern "C" fn helmholtz_smoluchowski_comp_entry(
    model: *mut c_void,
    node: *const c_int,
    dir_slot: *const c_double,
) -> c_double {
    unsafe {
        let direction = *(dir_slot as *const c_int);
        hs_comp(model, *node, direction)
    }
}

/// `getjouleheat_` — Joule heat source as a function of the electric field.
///
/// # Safety
///
/// Same contract as [`helmholtz_smoluchowski1_entry`].
#[export_name = "getjouleheat_"]
pub unsafe extern "C" fn getjouleheat_entry(
    model: *mut c_void,
    node: *const c_int,
    _dummy: *const c_double,
) -> c_double {
    unsafe { joule_heat(model, *node) }
}
