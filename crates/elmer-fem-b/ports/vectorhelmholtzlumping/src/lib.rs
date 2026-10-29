//! Port of `fem/src/modules/VectorHelmholtzLumping.F90` (Elmer @ 9f6af2f,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary
//! (whole `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/VectorHelmholtzLumping.F90` — one
//! F77-global entry `vectorhelmholtzlumping_` with three `CONTAINS`
//! internals. The emitted dynamic symbol is exactly
//! `vectorhelmholtzlumping_` (proven by `nm -D` on the reference `.so`);
//! the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/vectorhelmholtzlumping_acc.F90`, compiled against the
//! built `.mod` files of the target configure — the port itself contains
//! no layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MATMUL`), formatted I/O (`PRINT`/`WRITE`), and whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors VectorHelmholtzLumping.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: first-visit sizing (`NoModes` from the
//!   `Constraint Mode` scan), the `jMode` call counter with its
//!   `> NoModes` fatal, and the implicit-`SAVE` `Visited` latch live in
//!   the Rust `static` state with identical lifetimes; the internal
//!   procedures' `SAVE` (`Nodes`, `ParentNodes`, the twelve
//!   `ValueHandle`s, `WBasis`/`RotWBasis`/`Basis`/`dBasisdx`/`e_local`/
//!   `phi_local` plus `AllocationsDone`) live as wrapper-side SAVE state.
//! - The dead `np = n * MAXVAL(...)` statement (overwritten by `np = 0`
//!   before any use, pure and side-effect free) is omitted with this
//!   note; the effective `np = 0` (that is, `1:nd` ranges) is kept.
//! - Per integration point the wrapper returns the thirteen accumulation
//!   addends with `weight` (and the `1/(2*sqrt(2)*B)`, `/Zs` complex
//!   divisions) already folded in; the Rust loop performs only
//!   program-order adds, exactly the original statement order.
//! - Elementwise complex/real normalizations stay in Rust (naive,
//!   program order, no FMA); complex divisions stay in Fortran.
//! - `Model`, `dt`, `TransientSimulation` are received and ignored,
//!   exactly as in the original (the entry dummy `Solver` is likewise
//!   unused; the primary solver travels via the `Primary Solver Index`
//!   lookup, as in the original).
//! - `Mesh => GetMesh()` is kept as a touch call.
//! - Allocation failure calls the same `Fatal('VectorHelmholtzLumping',
//!   ...)` via `try_reserve`, instead of aborting.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn vhl_begin(model: *mut c_void, solver: *mut c_void);
    fn vhl_info_header();
    fn vhl_primary(soln: *mut c_int);
    fn vhl_omega() -> c_double;
    fn vhl_constants(m0: *mut c_double, e0: *mut c_double);
    fn vhl_count_modes() -> c_int;
    fn vhl_nodal_flag() -> c_int;
    fn vhl_fatal_nofield();
    fn vhl_find_potential() -> c_int;
    fn vhl_find_enodal(found: *mut c_int, dofs: *mut c_int);
    fn vhl_fatal_dofs(nd: c_int);
    fn vhl_info_nodal();
    fn vhl_info_edge();
    fn vhl_find_enamed(found: *mut c_int);
    fn vhl_check_edofs();
    fn vhl_edge_style(piola: *mut c_int, degree: *mut c_int);
    fn vhl_getmesh_touch();
    fn vhl_nof_bnd() -> c_int;
    fn vhl_begin_boundary(t: c_int) -> c_int;
    fn vhl_remap_edge();
    fn vhl_remap_face();
    fn vhl_bc_mode() -> c_int;
    fn vhl_pes() -> c_int;
    fn vhl_mype() -> c_int;
    fn vhl_par_z(z: *mut C64);
    fn vhl_par_r(x: *mut c_double);
    fn vhl_lib_setup(nodal: c_int, has_pot: c_int, first: c_int, nip: *mut c_int, found_e: *mut c_int);
    #[allow(clippy::too_many_arguments)]
    fn vhl_lib_ip(
        t: c_int,
        nodal: c_int,
        has_pot: c_int,
        omega: c_double,
        m0: c_double,
        e0: c_double,
        add: *mut c_double,
    );
    fn vhl_print_el_before(jm: c_int, nm: c_int, el: *const C64, cnorm: *const c_double);
    fn vhl_print_el_after(
        jm: c_int,
        nm: c_int,
        el: *const C64,
        wt: *const c_double,
        poynt: *const C64,
        curr: *const C64,
    );
    fn vhl_print_pot(jm: c_int, nm: c_int, pot: *const C64);
    fn vhl_center(nm: c_int, center: *mut c_double, wt: *mut c_double, cnode: *mut c_int);
    fn vhl_info_voltage();
    fn vhl_info_cannot();
    fn vhl_voltage(
        nm: c_int,
        jm: c_int,
        nodal: c_int,
        inpot: *const C64,
        volt: *mut C64,
        cnode: *const c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn vhl_write_files(
        nm: c_int,
        has_pot: c_int,
        poynt: *const C64,
        el: *const C64,
        curr: *const C64,
        pot: *const C64,
        volt: *const C64,
    );
    fn vhl_info_final();
    fn vhl_info_done();
    fn vhl_fatal_jmode();
    fn vhl_fatal_alloc();
}

/// Fortran `COMPLEX(KIND=dp)` on the wire (`c_double_complex`).
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct C64 {
    re: c_double,
    im: c_double,
}

/// First-visit SAVE state: the `NoModes × NoModes` complex accumulators,
/// the `NoModes × 3` centers, the three `NoModes` real rows, the center
/// nodes, and the `jMode` counter plus the implicit-`SAVE` `Visited`
/// latch. Sized on the first call from the `Constraint Mode` scan, never
/// resized afterwards — exactly the original lifetimes.
struct LumpState {
    visited: bool,
    nmodes: usize,
    jmode: usize,
    poynt: Vec<C64>,
    el: Vec<C64>,
    curr: Vec<C64>,
    pot: Vec<C64>,
    volt: Vec<C64>,
    center: Vec<c_double>,
    weight: Vec<c_double>,
    cnorm: Vec<c_double>,
    cnode: Vec<c_int>,
}

impl LumpState {
    const fn new() -> Self {
        Self {
            visited: false,
            nmodes: 0,
            jmode: 0,
            poynt: Vec::new(),
            el: Vec::new(),
            curr: Vec::new(),
            pot: Vec::new(),
            volt: Vec::new(),
            center: Vec::new(),
            weight: Vec::new(),
            cnorm: Vec::new(),
            cnode: Vec::new(),
        }
    }

    /// First-visit sizing. The original allocates once under
    /// `IF(.NOT. Visited)`; a later call with a different mode count
    /// would overflow there — unreachable in valid inputs, so growth is
    /// a no-op by construction and the vectors keep their first size.
    /// Reservation goes through `try_reserve_exact` so exhaustion calls
    /// the same `Fatal`, exactly like the original's `STAT=istat` check.
    fn ensure(&mut self, nm: usize) {
        if self.visited {
            return;
        }
        self.poynt = alloc_zero(nm * nm);
        self.el = alloc_zero(nm * nm);
        self.curr = alloc_zero(nm * nm);
        self.pot = alloc_zero(nm * nm);
        self.volt = alloc_zero(nm * nm);
        self.center = alloc_zero(nm * 3);
        self.weight = alloc_zero(nm);
        self.cnorm = alloc_zero(nm);
        self.cnode = alloc_zero(nm);
        self.nmodes = nm;
    }
}

fn alloc_zero<T: Clone + Default>(n: usize) -> Vec<T> {
    let mut v = Vec::new();
    if v.try_reserve_exact(n).is_err() {
        unsafe { vhl_fatal_alloc() };
        unreachable!();
    }
    v.resize(n, T::default());
    v
}

static STATE: Mutex<LumpState> = Mutex::new(LumpState::new());

/// `vectorhelmholtzlumping_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "vectorhelmholtzlumping_"]
pub unsafe extern "C" fn vectorhelmholtz_lumping_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { vhl_begin(model, solver) };
    unsafe { vhl_info_header() };

    // Primary solver lookup; Fatal(> Primary Solver Index <) when 0.
    let mut soln: c_int = 0;
    unsafe { vhl_primary(&mut soln) };

    // Mesh => GetMesh() — touch call, result unused by the body.
    unsafe { vhl_getmesh_touch() };

    let omega = unsafe { vhl_omega() };
    let (mut mu0inv, mut eps0) = (0.0, 0.0);
    unsafe { vhl_constants(&mut mu0inv, &mut eps0) };

    if !state.visited {
        let nm = unsafe { vhl_count_modes() } as usize;
        state.ensure(nm);
    }
    let nm = state.nmodes;
    state.jmode += 1;
    let jm = state.jmode - 1;
    if state.jmode > nm {
        unsafe { vhl_fatal_jmode() };
    }

    let has_pot = unsafe { vhl_find_potential() };

    // One should be able to toggle between using nodal or edge basis.
    // The results are not quite the same but should hopefully be close.
    let mut nodal = 0;
    let mut found = 0;
    let mut dofs = 0;
    let nodal_flag = unsafe { vhl_nodal_flag() };
    if nodal_flag != 0 {
        unsafe { vhl_info_nodal() };
        unsafe { vhl_find_enodal(&mut found, &mut dofs) };
        if found == 0 {
            unsafe { vhl_fatal_nofield() };
        }
        if dofs != 6 {
            unsafe { vhl_fatal_dofs(dofs) };
        }
        nodal = 1;
    } else {
        unsafe { vhl_info_edge() };
        unsafe { vhl_find_enamed(&mut found) };
        unsafe { vhl_check_edofs() };
        let mut piola = 0;
        let mut degree = 0;
        unsafe { vhl_edge_style(&mut piola, &mut degree) };
    }

    let active = unsafe { vhl_nof_bnd() };
    for i in 0..nm {
        // Fortran `IntCenter(NoModes,3)`: component-major, index k*nm+i.
        state.center[i] = 0.0;
        state.center[nm + i] = 0.0;
        state.center[2 * nm + i] = 0.0;
        state.weight[i] = 0.0;
        state.cnorm[i] = 0.0;
    }

    // `InitHandles` in the original: set .TRUE. on entry, consumed by the
    // first LocalIntegBC call (which re-runs ListInitElementKeyword).
    let mut first_lib = 1;
    for t in 1..=active {
        let family = unsafe { vhl_begin_boundary(t) };
        match family {
            1 => continue,
            2 => unsafe { vhl_remap_edge() },
            3 | 4 => unsafe { vhl_remap_face() },
            _ => {}
        }
        let imode = unsafe { vhl_bc_mode() };
        if imode == 0 {
            continue;
        }
        let im = (imode - 1) as usize;
        local_integ_bc(&mut state, nm, jm, im, nodal, has_pot, omega, mu0inv, eps0, first_lib);
        first_lib = 0;
    }

    if unsafe { vhl_pes() } > 1 {
        for i in 0..nm {
            let row = i * nm + jm;
            unsafe {
                vhl_par_z(&mut state.poynt[row]);
                vhl_par_z(&mut state.el[row]);
                vhl_par_z(&mut state.curr[row]);
                vhl_par_r(&mut state.weight[i]);
                vhl_par_r(&mut state.cnorm[i]);
            }
        }
    }

    unsafe { vhl_print_el_before(jm as c_int, nm as c_int, state.el.as_ptr(), state.cnorm.as_ptr()) };

    // Why sqrt(2) ?
    for i in 0..nm {
        let row = i * nm + jm;
        let w = state.cnorm[i];
        state.el[row].re /= w;
        state.el[row].im /= w;
    }
    state.el[jm * nm + jm].re -= 1.0;

    unsafe {
        vhl_print_el_after(
            jm as c_int,
            nm as c_int,
            state.el.as_ptr(),
            state.weight.as_ptr(),
            state.poynt.as_ptr(),
            state.curr.as_ptr(),
        )
    };
    if has_pot != 0 {
        if unsafe { vhl_pes() } > 1 {
            for i in 0..nm {
                let row = i * nm + jm;
                unsafe { vhl_par_z(&mut state.pot[row]) };
            }
        }
        for i in 0..nm {
            let row = i * nm + jm;
            let w = state.weight[i];
            state.pot[row].re /= w;
            state.pot[row].im /= w;
        }
        unsafe { vhl_print_pot(jm as c_int, nm as c_int, state.pot.as_ptr()) };
    }

    unsafe {
        vhl_center(
            nm as c_int,
            state.center.as_mut_ptr(),
            state.weight.as_mut_ptr(),
            state.cnode.as_mut_ptr(),
        )
    };

    if has_pot != 0 {
        unsafe { vhl_info_voltage() };
        unsafe {
            vhl_voltage(
                nm as c_int,
                jm as c_int,
                nodal,
                state.pot.as_ptr(),
                state.volt.as_mut_ptr(),
                state.cnode.as_ptr(),
            )
        };
    } else {
        unsafe { vhl_info_cannot() };
    }

    if state.jmode == nm && unsafe { vhl_mype() } == 0 {
        unsafe { vhl_info_final() };
        unsafe {
            vhl_write_files(
                nm as c_int,
                has_pot,
                state.poynt.as_ptr(),
                state.el.as_ptr(),
                state.curr.as_ptr(),
                state.pot.as_ptr(),
                state.volt.as_ptr(),
            )
        };
    }

    state.visited = true;
    unsafe { vhl_info_done() };
}

/// The `LocalIntegBC` internal procedure: per-boundary-element setup plus
/// the Gauss-point loop. The wrapper owns the setup/IP bodies (they hold
/// every `SUM`/`MATMUL` reduction); this loop only routes the thirteen
/// addends into the SAVE rows in the original statement order.
#[allow(clippy::too_many_arguments)]
fn local_integ_bc(
    state: &mut LumpState,
    _nm: usize,
    jm: usize,
    im: usize,
    nodal: c_int,
    has_pot: c_int,
    omega: c_double,
    mu0inv: c_double,
    eps0: c_double,
    first: c_int,
) {
    let nm = state.nmodes;
    let mut nip: c_int = 0;
    let mut found: c_int = 0;
    unsafe { vhl_lib_setup(nodal, has_pot, first, &mut nip, &mut found) };
    let mut add = [0.0 as c_double; 13];
    for t in 1..=nip {
        unsafe {
            vhl_lib_ip(t, nodal, has_pot, omega, mu0inv, eps0, add.as_mut_ptr());
        }
        state.weight[im] += add[0];
        state.center[im] += add[1];
        state.center[nm + im] += add[2];
        state.center[2 * nm + im] += add[3];
        state.cnorm[im] += add[4];
        let row = im * nm + jm;
        state.el[row].re += add[5];
        state.el[row].im += add[6];
        state.poynt[row].re += add[7];
        state.poynt[row].im += add[8];
        state.curr[row].re += add[9];
        state.curr[row].im += add[10];
        if has_pot != 0 {
            state.pot[row].re += add[11];
            state.pot[row].im += add[12];
        }
    }
}
