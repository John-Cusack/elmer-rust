//! Port of `fem/src/modules/Probes.F90` (Elmer @ 9f6af2f85, the serial
//! proof tree) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/Probes.F90` — one F77-global entry
//! `volumeprobe_` with no internal procedures. The emitted dynamic
//! symbol is exactly `volumeprobe_` (proven by `nm -D` on the reference
//! `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every
//! component access and every Elmer library call goes through the
//! generated-style Fortran `BIND(C)` wrappers in
//! `accessors/probes_acc.F90`, compiled against the built `.mod` files
//! of the target configure — the port itself contains no layout
//! assumption and no pinned offset. Order-sensitive work (the
//! `gausspoint=ngp` integrand evaluation, the left-associated
//! `integral + value*detj*s` accumulation, the `'res: volume probe '
//! // trim(name)` store) stays byte-identical by construction: the
//! library calls run in Fortran with the same call sequence, and the
//! accumulation keeps the original's evaluation order.
//!
//! Faithfulness notes (mirrors Probes.F90 statement for statement):
//! - `SAVE Nodes` lives as wrapper-side SAVE state with identical
//!   lifetime (filled per active element, consumed before the next
//!   fill, as in the original).
//! - The body-force probe-name string chopping + `ListInitElementKeyword`
//!   block is transcribed verbatim into the `prb_parse` wrapper (same
//!   comparisons, same order); byte-identity of the chopped names is
//!   guaranteed by keeping all string handling in Fortran.
//! - `integral(:)` is a per-call zeroed local, exactly as in the
//!   original (no Rust `static` state at all — the entry is reentrant;
//!   the `Mutex` is unnecessary here so there is none).
//! - `Model`, `dt`, `TransientSimulation` are unused by the body except
//!   through `Model % Simulation` in the store tail, exactly as in the
//!   original (received and forwarded as opaque handles); `Mesh =>
//!   GetMesh()` is kept as a touch call and `CoordinateSystemDimension()`
//!   is kept and discarded.
//! - Allocation failure calls `Fatal('VolumeProbe', ...)` via
//!   `try_reserve`, instead of aborting.
//! - Triage finding: when `Body force probe names` is absent (or blank),
//!   the reference faults on `SIZE()` of the never-allocated handle
//!   (SIGSEGV, exit 139); there the port yields zero probes and returns
//!   normally. No valid differential exists for that input, so the
//!   proof exercises single-name vs multi-name configs instead.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn prb_begin(h_model: *const c_void, h_solver: *const c_void);
    fn prb_touch_mesh();
    fn prb_coord_dim() -> c_int;
    fn prb_nof_active() -> c_int;
    fn prb_parse() -> c_int;
    fn prb_begin_active(
        t: c_int,
        n: *mut c_int,
        nd: *mut c_int,
        nb: *mut c_int,
        ngp: *mut c_int,
    );
    fn prb_nodes();
    fn prb_ip_point(
        q: c_int,
        u: *mut c_double,
        v: *mut c_double,
        w: *mut c_double,
        s: *mut c_double,
    );
    fn prb_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        nd: c_int,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn prb_probe_value(m: c_int, val: *mut c_double, found: *mut c_int);
    fn prb_store(n: c_int, vals: *mut c_double);
    fn prb_fatal_alloc();
}

/// `volumeprobe_` — the solver entry. Argument order and by-reference
/// travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `dt` and the transient
/// flag are unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` are ignored (received and dropped, as original).
#[export_name = "volumeprobe_"]
pub unsafe extern "C" fn volume_probe_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { prb_begin(model, solver) };

    // Mesh => GetMesh() — result unused by the body; kept as a touch call.
    unsafe { prb_touch_mesh() };

    // Body-force probe-name chop + ValueHandle init (verbatim in wrapper).
    let nprobes = unsafe { prb_parse() };

    // Per-call zeroed integral, exactly like the original's
    // `integral(counter) = 0.0` init loop.
    let mut integral: Vec<c_double> = Vec::new();
    if nprobes > 0 {
        let n = nprobes as usize;
        if integral.try_reserve(n).is_err() {
            unsafe { prb_fatal_alloc() };
            std::process::abort();
        }
        integral.resize(n, 0.0);
    }

    // cdim = CoordinateSystemDimension() — unused by the body.
    let _cdim: c_int = unsafe { prb_coord_dim() };

    // Active = GetNOFActive(Solver)
    let active = unsafe { prb_nof_active() };
    for t in 1..=active {
        // n = nodes, nd = dofs, nb = bdofs; the original's Basis is
        // sized nd = dofs + bdofs.
        let (mut _n, mut nd, mut nb, mut ngp) = (0, 0, 0, 0);
        unsafe { prb_begin_active(t, &mut _n, &mut nd, &mut nb, &mut ngp) };

        // CALL GetElementNodes(Nodes, UElement=elm) — wrapper SAVE Nodes.
        unsafe { prb_nodes() };

        for q in 1..=ngp {
            let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
            unsafe { prb_ip_point(q, &mut u, &mut v, &mut w, &mut s) };
            let (mut detj, mut stat) = (0.0, 0);
            unsafe { prb_eleminfo(u, v, w, nd + nb, &mut detj, &mut stat) };
            let _ = stat;
            for m in 1..=nprobes {
                let (mut val, mut found) = (0.0, 0);
                unsafe { prb_probe_value(m, &mut val, &mut found) };
                if found != 0 {
                    // integral(m) = integral(m) + value*detj*s,
                    // left-associated exactly as written.
                    let slot = &mut integral[(m - 1) as usize];
                    *slot += (val * detj) * s;
                }
            }
        }
    }

    // Parallel guard + reduction + `res: volume probe <name>` store tail.
    unsafe { prb_store(nprobes, integral.as_mut_ptr()) };
}
