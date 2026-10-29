//! Port of `fem/src/modules/ParticleAdvector.F90` (Elmer @ 9f6af2f85, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ParticleAdvector.F90` — two F77-global
//! entries `particleadvector_` (solver) and `particleadvector_init_`
//! (initializer); the emitted dynamic symbols are exactly those two
//! (proven by `nm -D` on the reference `.so`). The `CONTAINS` internals
//! emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/particleadvector_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MATMUL`), formatted I/O (`WRITE`/`Info`), and whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ParticleAdvector.F90 statement for statement):
//! - `SAVE` semantics replicated: `VisitedTimes`, `Nstep`, `TimeOrder`,
//!   `tottime`, `TimeStepsTaken`, `ParticleStepsTaken`, `ParticleInfo`
//!   live in the Rust `STATE` guard with identical lifetimes; the
//!   internal procedures' `SAVE` locals (`SetParticleVelocities` mesh /
//!   basis / params, `SetAdvectedField` `PrevNorm` / `UnitPerm`) live as
//!   wrapper-side SAVE state. `dertime` is a per-invocation local, exactly
//!   as in the original (it keeps accumulating across the reverse-time
//!   `GOTO 1` restart within one invocation).
//! - `OrigCoordinate` (an automatic array of the entry, host-associated
//!   into `SetAdvectedField`) is mirrored as a wrapper-side SAVE
//!   allocatable: freed at entry start (fresh automatic), captured when
//!   the entry allocates it, freed at entry return.
//! - `Model`, `dt`, `TransientSimulation` are received and stashed
//!   (`Model` is used by the fixed-particle wrapper); `dt` and the
//!   transient flag are unused by the body, exactly as in the original.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.
//!   The entry never re-enters Rust from Fortran, so holding the guard
//!   across the solve cannot deadlock.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pav_begin(h_model: *const c_void, h_solver: *const c_void);
    fn pav_dim() -> c_int;
    fn pav_info_banner();
    fn pav_done();
    fn pav_reverse_query(val: *mut c_int, got: *mut c_int);
    fn pav_reinit_query() -> c_int;
    fn pav_reverse_conflict(reverse: c_int, gotit: c_int, reinit: c_int) -> c_int;
    fn pav_fix_tangent_query() -> c_int;
    fn pav_apply_tangent(undo: c_int);
    fn pav_fix_undo();
    fn pav_first_init(torder: *mut c_int, nstep: *mut c_int, pinfo: *mut c_int);
    fn pav_init_particles();
    fn pav_orig_capture();
    fn pav_orig_free();
    fn pav_particle_create();
    fn pav_particle_zero();
    fn pav_iorder() -> c_int;
    fn pav_dtsign_neg();
    fn pav_maxdt(initts: c_int) -> c_double;
    fn pav_is_tiny(m: c_double) -> c_int;
    fn pav_tiny_exit(im1: c_int);
    fn pav_num_particles() -> c_int;
    fn pav_advance(istep: c_int);
    fn pav_info_active(level: c_int) -> c_int;
    fn pav_status_count();
    fn pav_locate();
    fn pav_path_integral(reverse: c_int, istep: c_int);
    fn pav_nomoving() -> c_int;
    fn pav_info_step(i: c_int, nomoving: c_int);
    fn pav_retired();
    fn pav_particle_info(psteps: c_int, tsteps: c_int, tottime: c_double);
    fn pav_set_advected(dertime: c_double, dim: c_int, reverse: c_int);
    fn pav_reverse_turn(reverse: c_int) -> c_int;
    fn pav_setvel(first: c_int);
    fn pav_set_fixed();
    fn pav_init_body();
}

/// First-mesh SAVE state of the `ParticleAdvector` entry: the `SAVE`
/// locals `VisitedTimes`, `Nstep`, `TimeOrder`, `tottime`,
/// `TimeStepsTaken`, `ParticleStepsTaken`, `ParticleInfo`.
struct ParticleAdvectorState {
    visited_times: c_int,
    nstep: c_int,
    time_order: c_int,
    tottime: c_double,
    time_steps_taken: c_int,
    particle_steps_taken: c_int,
    particle_info: bool,
}

impl ParticleAdvectorState {
    const fn new() -> Self {
        Self {
            visited_times: 0,
            nstep: 0,
            time_order: 0,
            tottime: 0.0,
            time_steps_taken: 0,
            particle_steps_taken: 0,
            particle_info: false,
        }
    }
}

static STATE: Mutex<ParticleAdvectorState> = Mutex::new(ParticleAdvectorState::new());

/// `particleadvector_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list `(Model, Solver, dt,
/// TransientSimulation)`; `dt` and the transient flag are unused by the
/// body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "particleadvector_"]
pub unsafe extern "C" fn particleadvector_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    pav_begin(model, solver);
    // Fresh automatic-array lifetime for the OrigCoordinate mirror.
    pav_orig_free();
    pav_info_banner();

    let mut st = STATE.lock().unwrap();
    st.visited_times += 1;
    let visited = st.visited_times;

    let dim = pav_dim();

    let mut reverse: c_int = 0;
    let mut gotit: c_int = 0;
    pav_reverse_query(&mut reverse, &mut gotit);
    let reinit_q = pav_reinit_query();
    let reinit = pav_reverse_conflict(reverse, gotit, reinit_q);

    if pav_fix_tangent_query() != 0 {
        pav_apply_tangent(0);
    }

    if visited == 1 {
        let mut torder: c_int = 0;
        let mut nstep: c_int = 0;
        let mut pinfo: c_int = 0;
        pav_first_init(&mut torder, &mut nstep, &mut pinfo);
        st.time_order = torder;
        st.nstep = nstep;
        st.particle_info = pinfo != 0;
    }

    // Do some initialization: allocate space, check fields.
    let mut init_all_velo: c_int = 1;
    if visited == 1 || reinit != 0 {
        pav_init_particles();
        if reverse != 0 {
            pav_orig_capture();
        }
    } else {
        // In case the velocity field is changed update also the particle velocities.
        pav_setvel(1);
        init_all_velo = 0;
    }

    if visited == 1 {
        pav_particle_create();
    }
    pav_particle_zero();

    // Freeze particles that are known not to move (e.g. no-slip wall).
    pav_set_fixed();

    // We are advecting backwards in time!
    let mut init_timestep: c_int = 1;
    let iorder = pav_iorder();

    // This is the default for particle advector!
    pav_dtsign_neg();

    let mut dertime: c_double = 0.0;
    // Label 1: the reverse-time restart jumps back here.
    loop {
        let mut i: c_int = 1;
        while i <= st.nstep {
            // Get the timestep size, initialize at 1st round.
            let maxdt = pav_maxdt(init_timestep);

            // If size of timestep goes to zero then no more steps are needed.
            if pav_is_tiny(maxdt) != 0 {
                pav_tiny_exit(i - 1);
                break;
            }

            dertime += maxdt;
            st.tottime += maxdt;

            st.time_steps_taken += 1;
            st.particle_steps_taken += pav_num_particles();

            // Loop over possible runge-kutta steps.
            let mut istep: c_int = 1;
            while istep <= iorder {
                // Advance the coordinate, r = r0 + v * dt, then locate the
                // particles and refresh the velocity at the new location.
                pav_advance(istep);
                if pav_info_active(20) != 0 {
                    pav_status_count();
                }
                pav_locate();
                pav_setvel(init_all_velo);
                init_all_velo = 0;

                // Integrate over the particle path.
                pav_path_integral(reverse, istep);

                init_timestep = 0;
                istep += 1;
            }

            let nomoving = pav_nomoving();
            pav_info_step(i, nomoving);

            // Freeze particles that are too old.
            pav_retired();

            if pav_info_active(15) != 0 {
                pav_particle_info(st.particle_steps_taken, st.time_steps_taken, st.tottime);
            }
            i += 1;
        }

        // Set the advected field giving the final locations of the particles
        // backward in time.
        pav_set_advected(dertime, dim, reverse);

        // In the end show some statistical info.
        if st.particle_info {
            pav_particle_info(st.particle_steps_taken, st.time_steps_taken, st.tottime);
        }

        if pav_reverse_turn(reverse) == 0 {
            break;
        }
        init_timestep = 1;
    }

    pav_fix_undo();

    pav_done();
    pav_orig_free();
}

/// `particleadvector_init_` — the solver initializer entry. Same dummy
/// list as the solver entry; only the solver handle is used.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "particleadvector_init_"]
pub unsafe extern "C" fn particleadvector_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    pav_begin(model, solver);
    pav_init_body();
}
