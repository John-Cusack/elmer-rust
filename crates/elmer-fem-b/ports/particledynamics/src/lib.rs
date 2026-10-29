//! Port of `fem/src/modules/ParticleDynamics.F90` (Elmer @ 9f6af2f85)
//! for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/ParticleDynamics.F90` — one F77-global
//! entry `particledynamics_` with four `CONTAINS` internals, an init
//! entry `particledynamics_init_`, and a module
//! (`ParticleDynamicsStuff`) with three exported procedures that the
//! library calls back as procedure arguments (`LocateParticles`,
//! `ParticleParticleInteraction`). The emitted dynamic symbols are
//! exactly the five `T` symbols of the reference `.so` (proven by
//! `nm -D`); the internals emit only local symbols.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries'
//! derived-type dummies (`Model`, `Solver`) cross as opaque addresses.
//! Every component access and every Elmer library call goes through
//! the generated-style Fortran `BIND(C)` wrappers in
//! `accessors/particledynamics_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`SUM`, `MATMUL`), formatted I/O (`PRINT`/`WRITE`) and whole-array
//! updates stay in Fortran wrappers, so bit-identity holds by
//! construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ParticleDynamics.F90 statement for
//! statement):
//! - `SAVE` semantics replicated explicitly: each entry/internal
//!   procedure owns a `Mutex`-guarded state struct holding its SAVE
//!   scalars with identical lifetimes; first-mesh sizings and
//!   Elmer-owned lookups (Basis/dBasisdx arrays, variable handles, BC
//!   handles, entry pointers) live wrapper-side with identical
//!   lifetimes (consumed before the next fill, as in the original).
//! - The three module procedures are re-entered by the Elmer library
//!   on this thread (`LocateParticles` → wall kernel,
//!   `ParticleParticleInteraction` → contact/collision kernel) through
//!   plain-Fortran trampolines matching the library's explicit
//!   interfaces; each callback owns a disjoint state struct, so the
//!   main entry never holds a callback's guard (nor vice versa) — no
//!   non-reentrant guard is ever held across a re-entrant call.
//! - Missing value-list keys read back the library defaults (0,
//!   0.0, `.FALSE.`, `''`) exactly as the reference build does, since
//!   the same `Get*/ListGet*` routines answer through the wrappers.
//! - No Rust threads: the solver entries run serially, so each
//!   `Mutex` is only a `static`-sharing guard, never contended across
//!   threads (`Send` is asserted by hand for the opaque handles, which
//!   never cross threads in the serial solver).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn pda_m_begin(h_solver: *const c_void, dim: *mut c_int);
    fn pda_m_get_int(which: c_int, val: *mut c_int, found: *mut c_int);
    fn pda_m_get_log(which: c_int, val: *mut c_int, found: *mut c_int);
    fn pda_m_add_truecollision();
    fn pda_m_set_prelim(timeorder: c_int);
    fn pda_m_ngroups_set(n: c_int);
    fn pda_m_init_particles(group: c_int, append: c_int);
    fn pda_m_push_ns(group: c_int);
    fn pda_m_pop_ns();
    fn pda_m_info_which(which: c_int, a: c_int);
    fn pda_m_info_groupcount(group: c_int, cnt: c_int);
    fn pda_m_group_count(group: c_int) -> c_int;
    fn pda_m_group_body(i: c_int, bodyid: *mut c_int);
    fn pda_m_group_set(i: c_int, g: c_int);
    fn pda_m_groupcond_begin(group: c_int) -> c_int;
    fn pda_m_groupcond_eval(i: c_int, coord: *mut c_double) -> c_double;
    fn pda_m_getcoord(i: c_int, coord: *mut c_double);
    fn pda_m_pdist_create();
    fn pda_m_release_waiting();
    fn pda_m_statistics(which: c_int);
    fn pda_m_timestep(isfirst: c_int, dtime: *mut c_double);
    fn pda_m_init_time();
    fn pda_m_locate();
    fn pda_m_eliminate();
    fn pda_m_nparticles() -> c_int;
    fn pda_m_box_contact();
    fn pda_m_neigh_create();
    fn pda_m_interact_collision(dtime: c_double);
    fn pda_m_interact_contact(dtime: c_double);
    fn pda_m_destroy_ghosts();
    fn pda_m_advance();
    fn pda_m_box_periodic();
    fn pda_m_delete_lost();
    fn pda_m_path_integral();
    fn pda_m_save(h_model: *const c_void, dt: c_double, transient: c_int);
    fn pda_m_time_init(nstep: c_int);
    fn pda_m_time_check(i: c_int, nstep: c_int);
    fn pda_m_info_alldone();
    fn pda_m_info_start();
    fn pda_m_pconst(which: c_int) -> c_int;
    fn pda_m_particle_info(steps: c_int, timesteps: c_int, tottime: c_double);
    fn pda_n_v2sum(dim: c_int, v2sum: *mut c_double, n: *mut c_int);
    fn pda_n_reduce(n: *mut c_int, v2sum: *mut c_double);
    fn pda_n_finish(prevnrm: *mut c_double, nrm: *const c_double);
    fn pda_s_first(timeorder: *mut c_int, dim: *mut c_int, mfp: *mut c_double);
    fn pda_s_dt() -> c_double;
    fn pda_s_try(
        no: c_int,
        dt: c_double,
        dim: c_int,
        mfp: c_double,
        timeorder: c_int,
        scattered: *mut c_int,
    );
    fn pda_s_status(no: c_int) -> c_int;
    fn pda_s_finish(noscattered: c_int, h_model: *const c_void);
    fn pda_f_touch();
    fn pda_f_first_alloc(dim: *mut c_int, timeorder: *mut c_int);
    fn pda_f_first_gravity();
    fn pda_f_first_var(which: c_int);
    fn pda_f_first_fields();
    fn pda_f_first_tail(particletofield: c_int);
    fn pda_f_reset(dtime: c_double, fieldreset: c_int, usercoeff: *mut c_double);
    fn pda_f_counts(nop: *mut c_int, nogrp: *mut c_int);
    fn pda_f_debug();
    #[allow(clippy::too_many_arguments)]
    fn pda_f_scalars(
        mass: *mut c_double,
        damping: *mut c_double,
        rad: *mut c_double,
        charge: *mut c_double,
        drag: *mut c_double,
        gm: *mut c_int,
        gd: *mut c_int,
        gr: *mut c_int,
        gc: *mut c_int,
        gdrag: *mut c_int,
    );
    fn pda_f_groups_fetch();
    fn pda_f_checks();
    fn pda_f_buoyvol(rad: c_double);
    fn pda_f_state(s: *mut c_int);
    fn pda_f_begin_no(
        no: c_int,
        status: *mut c_int,
        elemidx: *mut c_int,
        group: *mut c_int,
        proceed: *mut c_int,
    );
    fn pda_f_group_update(
        group: c_int,
        mass: *mut c_double,
        damping: *mut c_double,
        rad: *mut c_double,
        charge: *mut c_double,
        drag: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn pda_f_force(
        dtime: c_double,
        mass: c_double,
        damping: c_double,
        charge: c_double,
        drag: c_double,
        rad: c_double,
        status: *mut c_int,
        proceed: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn pda_f_accum(
        no: c_int,
        dtime: c_double,
        tottime: c_double,
        usercoeff: c_double,
        mass: c_double,
        charge: c_double,
        group: c_int,
    );
    fn pda_f_setparticles(no: c_int, status: c_int);
    fn pda_f_normalize();
    fn pda_f_set_prevdtime(dtime: c_double);
    fn pda_t_first(h_model: *const c_void, dim: *mut c_int, anyinter: *mut c_int);
    fn pda_t_state(s: *mut c_int);
    fn pda_t_begin_no(no: c_int, status: *mut c_int, elemidx: *mut c_int, proceed: *mut c_int);
    fn pda_t_imax() -> c_int;
    fn pda_t_face(face_i: c_int, found: *mut c_int, collision: *mut c_int, contact: *mut c_int);
    fn pda_t_geom(
        no: c_int,
        dim: c_int,
        coord: *mut c_double,
        normal: *mut c_double,
        dist: *mut c_double,
    );
    fn pda_t_contact(
        no: c_int,
        dt: c_double,
        dim: c_int,
        dist: c_double,
        hit: *mut c_int,
        exitface: *mut c_int,
    );
    #[allow(clippy::too_many_arguments)]
    fn pda_t_collision(
        no: c_int,
        dt: c_double,
        dim: c_int,
        dist: c_double,
        hit: *mut c_int,
        exitface: *mut c_int,
        skipface: *mut c_int,
        walltrace: *mut c_int,
    );
    fn pda_t_trace(no: c_int, dim: c_int);
    fn pda_c_first(
        rad: *mut c_double,
        mass: *mut c_double,
        spring: *mut c_double,
        damping: *mut c_double,
        friction: *mut c_double,
    );
    fn pda_c_dot(a: *const c_double, b: *const c_double) -> c_double;
    fn pda_c_warn_samepoint();
    fn pda_k_first(
        rad1: *mut c_double,
        mass1: *mut c_double,
        coeff: *mut c_double,
        truecollision: *mut c_int,
    );
    fn pda_k_dot(a: *const c_double, b: *const c_double) -> c_double;
    fn pda_w_first(meshdim: *mut c_int, particlebunch: *mut c_int);
    fn pda_w_face(faceindex: c_int);
    fn pda_w_advance(lambda: c_double, rinit: *mut c_double, rfin: *mut c_double);
    fn pda_w_reflect(
        rinit: *mut c_double,
        rfin: *mut c_double,
        vinit: *mut c_double,
        vfin: *mut c_double,
    );
    fn pda_w_bunch(no: c_int, rinit: *mut c_double);
    fn pda_i_init(h_solver: *const c_void);
}

/// Value-list logical/integer fetch helpers: the library answers with
/// the same defaults (0 / `.FALSE.`) the reference build observes.
unsafe fn mget_int(which: c_int) -> (i32, bool) {
    let mut v: c_int = 0;
    let mut f: c_int = 0;
    pda_m_get_int(which, &mut v, &mut f);
    (v, f != 0)
}

unsafe fn mget_log(which: c_int) -> (bool, bool) {
    let mut v: c_int = 0;
    let mut f: c_int = 0;
    pda_m_get_log(which, &mut v, &mut f);
    (v != 0, f != 0)
}

/// Main-entry SAVE state (mirrors the `SAVE` list at the top of
/// `ParticleDynamics`, identical lifetimes: first-visit sizing, never
/// reset per call).
struct PdState {
    visited_times: i32,
    collision_interaction: bool,
    contact_interaction: bool,
    no_groups: i32,
    particle_to_field: bool,
    output_interval: i32,
    nstep: i32,
    time_order: i32,
    particle_in_box: bool,
    tottime: f64,
    time_steps_taken: i32,
    particle_steps_taken: i32,
    particle_wall: bool,
    stat_info: bool,
    particle_info: bool,
    time_info: bool,
    true_collision: bool,
    field_reset: bool,
    particles_located: bool,
    do_scattering: bool,
    h_model: *mut c_void,
    dt: f64,
    transient: bool,
}

impl PdState {
    const fn new() -> Self {
        Self {
            visited_times: 0,
            collision_interaction: false,
            contact_interaction: false,
            no_groups: 0,
            particle_to_field: false,
            output_interval: 0,
            nstep: 0,
            time_order: 0,
            particle_in_box: false,
            tottime: 0.0,
            time_steps_taken: 0,
            particle_steps_taken: 0,
            particle_wall: false,
            stat_info: false,
            particle_info: false,
            time_info: false,
            true_collision: false,
            field_reset: false,
            particles_located: false,
            do_scattering: false,
            h_model: std::ptr::null_mut(),
            dt: 0.0,
            transient: false,
        }
    }
}

// The serial solver never moves these handles across threads; the
// `Mutex` is only a `static`-sharing guard, never contended.
unsafe impl Send for PdState {}

static PD: Mutex<PdState> = Mutex::new(PdState::new());

/// `ParticleFieldInteraction` SAVE mirror: only the first-visit flag
/// is needed on the Rust side (per-call flags travel through
/// `pda_f_state`); the arrays, handles and `PrevGroup`/`CurrGroup`
/// live wrapper-side with identical lifetimes.
struct PfiState {
    visited: bool,
}

impl PfiState {
    const fn new() -> Self {
        Self { visited: false }
    }
}

unsafe impl Send for PfiState {}

static PFI: Mutex<PfiState> = Mutex::new(PfiState::new());

/// `ParticleWallContact` SAVE mirror: sizing/flags; the mesh data,
/// Basis arrays, nodes, variable handles and `val` live wrapper-side.
struct PwcState {
    visited: bool,
    dim: i32,
    any_interaction: bool,
}

impl PwcState {
    const fn new() -> Self {
        Self {
            visited: false,
            dim: 0,
            any_interaction: false,
        }
    }
}

unsafe impl Send for PwcState {}

static PWC: Mutex<PwcState> = Mutex::new(PwcState::new());

/// `ParticleScattering` SAVE mirror.
struct PscState {
    visited: bool,
    mfp: f64,
    timeorder: c_int,
    dim: c_int,
}

impl PscState {
    const fn new() -> Self {
        Self {
            visited: false,
            mfp: 0.0,
            timeorder: 0,
            dim: 0,
        }
    }
}

unsafe impl Send for PscState {}

static PSC: Mutex<PscState> = Mutex::new(PscState::new());

/// `ParticleParticleContact` SAVE mirror.
struct PpcState {
    visited: bool,
    rad: f64,
    mass: f64,
    spring: f64,
    damping: f64,
    friction: f64,
}

impl PpcState {
    const fn new() -> Self {
        Self {
            visited: false,
            rad: 0.0,
            mass: 0.0,
            spring: 0.0,
            damping: 0.0,
            friction: 0.0,
        }
    }
}

unsafe impl Send for PpcState {}

static PPC: Mutex<PpcState> = Mutex::new(PpcState::new());

/// `ParticleParticleCollision` SAVE mirror.
struct PpkState {
    visited: bool,
    rad1: f64,
    mass1: f64,
    coeff: f64,
    true_collision: bool,
}

impl PpkState {
    const fn new() -> Self {
        Self {
            visited: false,
            rad1: 0.0,
            mass1: 0.0,
            coeff: 1.0,
            true_collision: false,
        }
    }
}

unsafe impl Send for PpkState {}

static PPK: Mutex<PpkState> = Mutex::new(PpkState::new());

/// `ParticleWallProc` SAVE mirror.
struct PwpState {
    visited: bool,
    mesh_dim: i32,
    particle_bunch: bool,
}

impl PwpState {
    const fn new() -> Self {
        Self {
            visited: false,
            mesh_dim: 0,
            particle_bunch: false,
        }
    }
}

unsafe impl Send for PwpState {}

static PWP: Mutex<PwpState> = Mutex::new(PwpState::new());

/// `particledynamics_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `Model`, `dt` and the
/// transient flag are used exactly where the original uses them
/// (`SaveParticleData`, `GetParticleTimeStep` sizing).
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes;
/// `dt`/`transient` must be valid reads (they are forwarded, as in
/// the original).
#[export_name = "particledynamics_"]
pub unsafe extern "C" fn particle_dynamics_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    transient: *const c_int,
) {
    let mut st = PD.lock().unwrap();
    pda_m_info_start();
    st.visited_times += 1;
    let mut dim: c_int = 0;
    pda_m_begin(solver, &mut dim);
    st.h_model = model;
    st.dt = *dt;
    st.transient = *transient != 0;
    // First-visit initialization: allocate space, check fields.
    if st.visited_times == 1 {
        let (to, found_to) = mget_int(1);
        st.time_order = if found_to { to } else { 2 };
        let (nog, found_groups) = mget_int(2);
        if found_groups {
            pda_m_ngroups_set(nog);
            st.no_groups = nog;
        } else {
            st.no_groups = 0;
        }
        pda_m_set_prelim(st.time_order);
        let (ptf, _) = mget_log(1);
        st.particle_to_field = ptf;
        let (coll, _) = mget_log(2);
        st.collision_interaction = coll;
        let (cont, _) = mget_log(3);
        st.contact_interaction = cont;
        let (outi, _) = mget_int(3);
        st.output_interval = outi;
        let (nstep, found_nstep) = mget_int(4);
        st.nstep = if found_nstep { nstep } else { 1 };
        let (boxcoll, _) = mget_log(4);
        let (boxcont, _) = mget_log(5);
        st.particle_in_box = boxcoll || boxcont;
        let (wall, _) = mget_log(6);
        st.particle_wall = wall;
        let (scat, _) = mget_log(7);
        st.do_scattering = scat;
        st.true_collision = false;
        pda_m_add_truecollision();
        let (stati, _) = mget_log(8);
        st.stat_info = stati;
        let (parti, _) = mget_log(9);
        st.particle_info = parti;
        let (timei, _) = mget_log(10);
        st.time_info = timei;
        st.particles_located = false;
        if st.particle_to_field {
            let (fr, _) = mget_log(11);
            st.field_reset = fr;
        }
    }
    // Initialize particles at first visiting, or each time if requested.
    if st.visited_times == 1 || mget_log(12).0 {
        if st.no_groups > 0 {
            let (group0, found_g0) = mget_int(5);
            let mut group = if found_g0 { group0 } else { 1 };
            if mget_log(13).0 {
                pda_m_init_particles(group, 0);
                pda_m_info_which(1, 0);
                let n = pda_m_nparticles();
                let mut i = 1;
                while i <= n {
                    let mut bodyid: c_int = 0;
                    pda_m_group_body(i, &mut bodyid);
                    pda_m_group_set(i, bodyid);
                    i += 1;
                }
            } else if mget_log(14).0 {
                pda_m_init_particles(group, 0);
                pda_m_info_which(4, 0);
                group = 1;
                while group <= st.no_groups {
                    if pda_m_groupcond_begin(group) != 0 {
                        let n = pda_m_nparticles();
                        let mut i = 1;
                        while i <= n {
                            let mut coord = [0.0 as c_double; 3];
                            pda_m_getcoord(i, coord.as_mut_ptr());
                            let cond = pda_m_groupcond_eval(i, coord.as_mut_ptr());
                            if cond > 0.0 {
                                pda_m_group_set(i, group);
                            }
                            i += 1;
                        }
                    }
                    group += 1;
                }
            } else {
                group = 1;
                while group <= st.no_groups {
                    pda_m_info_which(2, group);
                    pda_m_push_ns(group);
                    pda_m_init_particles(group, 1);
                    pda_m_pop_ns();
                    group += 1;
                }
            }
            group = 1;
            while group <= st.no_groups {
                let j = pda_m_group_count(group);
                pda_m_info_groupcount(group, j);
                group += 1;
            }
        } else {
            pda_m_init_particles(-1, 0);
        }
        st.particles_located = false;
        if mget_log(15).0 {
            pda_m_pdist_create();
        }
    }
    pda_m_release_waiting();
    if st.stat_info {
        pda_m_statistics(0);
        pda_m_statistics(1);
    }
    if st.time_info {
        pda_m_time_init(st.nstep);
    }
    let mut i = 1;
    let mut dtime_last: c_double = 0.0;
    while i <= st.nstep {
        let mut dtime: c_double = 0.0;
        pda_m_timestep(if i == 1 { 1 } else { 0 }, &mut dtime);
        dtime_last = dtime;
        if dtime < f64::MIN_POSITIVE {
            break;
        }
        st.tottime += dtime;
        let assembly_loop = st.particle_to_field && i > 1 && !st.field_reset;
        pda_m_init_time();
        if !(i == 1 && st.particles_located) {
            pda_m_locate();
            st.particles_located = true;
        }
        pda_m_eliminate();
        particle_field_interaction(dtime, true, assembly_loop, &st);
        let noparticles = pda_m_nparticles();
        st.particle_steps_taken += noparticles;
        st.time_steps_taken += 1;
        if st.particle_in_box {
            pda_m_box_contact();
        } else if st.particle_wall {
            particle_wall_contact(dtime, st.h_model);
        }
        if st.collision_interaction || st.contact_interaction {
            pda_m_neigh_create();
            if st.collision_interaction {
                pda_m_interact_collision(dtime);
            }
            if st.contact_interaction {
                pda_m_interact_contact(dtime);
            }
            pda_m_destroy_ghosts();
        }
        pda_m_advance();
        if st.do_scattering {
            particle_scattering(st.h_model);
        }
        pda_m_box_periodic();
        pda_m_delete_lost();
        pda_m_path_integral();
        if st.output_interval > 0 && (i % st.output_interval) == 0 {
            pda_m_save(st.h_model, st.dt, if st.transient { 1 } else { 0 });
        }
        if st.time_info {
            pda_m_time_check(i, st.nstep);
        }
        i += 1;
    }
    if st.particle_to_field {
        pda_m_locate();
        pda_m_eliminate();
        particle_field_interaction(dtime_last, false, true, &st);
        st.particles_located = true;
    }
    if st.stat_info {
        pda_m_statistics(0);
        pda_m_statistics(1);
        pda_m_statistics(3);
        pda_m_statistics(4);
    }
    let do_pinfo = st.particle_info;
    let steps_taken = st.particle_steps_taken;
    let time_taken = st.time_steps_taken;
    let tottime = st.tottime;
    drop(st);
    if do_pinfo {
        pda_m_particle_info(steps_taken, time_taken, tottime);
    }
    pseudo_norm(dim);
    pda_m_info_alldone();
}

/// The `PseudoNorm` internal procedure: L2 norm of the velocity field
/// for consistency tests. Reductions and the `Solver % Variable`
/// update stay in Fortran wrappers (identical call sequence); the
/// square root matches the original statement order.
unsafe fn pseudo_norm(dim: c_int) {
    let mut v2sum: c_double = 0.0;
    let mut n: c_int = 0;
    pda_n_v2sum(dim, &mut v2sum, &mut n);
    pda_n_reduce(&mut n, &mut v2sum);
    let mut prevnrm: c_double = 0.0;
    let nrm = (v2sum / f64::from(n)).sqrt();
    pda_n_finish(&mut prevnrm, &nrm);
    let _ = prevnrm;
}

/// The `ParticleScattering` internal procedure. Random draws
/// (`EvenRandom`), the rotation and the messages stay in Fortran
/// wrappers; Rust drives the status-gated loop and the count.
/// `TimeOrder`/`dim` mirror the host-associated SAVE values: set
/// once at first visit, reused afterwards.
unsafe fn particle_scattering(h_model: *mut c_void) {
    let mut sc = PSC.lock().unwrap();
    if !sc.visited {
        let mut timeorder: c_int = 0;
        let mut dim: c_int = 0;
        let mut mfp: c_double = 0.0;
        pda_s_first(&mut timeorder, &mut dim, &mut mfp);
        sc.mfp = mfp;
        sc.timeorder = timeorder;
        sc.dim = dim;
        sc.visited = true;
    }
    let (mfp, timeorder, dim) = (sc.mfp, sc.timeorder, sc.dim);
    drop(sc);
    let n = pda_m_nparticles();
    let lost = pda_m_pconst(1);
    let initiated = pda_m_pconst(2);
    let wallb = pda_m_pconst(3);
    let dt = pda_s_dt();
    let mut no_scattered = 0;
    let mut no = 1;
    while no <= n {
        let status = pda_s_status(no);
        let mut proceed = true;
        if status >= lost {
            proceed = false;
        }
        if status <= initiated {
            proceed = false;
        }
        if status == wallb {
            proceed = false;
        }
        if proceed {
            let mut scattered: c_int = 0;
            pda_s_try(no, dt, dim, mfp, timeorder, &mut scattered);
            if scattered != 0 {
                no_scattered += 1;
            }
        }
        no += 1;
    }
    pda_s_finish(no_scattered, h_model);
}

/// The `ParticleFieldInteraction` internal procedure: field values at
/// the particle positions (gravity, damping, potentials, drag,
/// buoyancy) and/or accumulation of particle data onto the mesh
/// fields. All mesh/variable access, reductions and string handling
/// stay in Fortran wrappers; Rust drives the per-particle loop and
/// holds the `Visited` flag plus per-call scalar updates.
unsafe fn particle_field_interaction(
    dtime: f64,
    set_particles: bool,
    set_fields: bool,
    st: &PdState,
) {
    pda_f_touch();
    let mut pfi = PFI.lock().unwrap();
    if !pfi.visited {
        let mut dim: c_int = 0;
        let mut timeorder: c_int = 0;
        pda_f_first_alloc(&mut dim, &mut timeorder);
        pda_f_first_gravity();
        let mut which = 1;
        while which <= 6 {
            pda_f_first_var(which);
            which += 1;
        }
        if st.particle_to_field {
            pda_f_first_fields();
        }
        pda_f_first_tail(if st.particle_to_field { 1 } else { 0 });
        pfi.visited = true;
    }
    drop(pfi);
    let mut s = [0 as c_int; 18];
    pda_f_state(s.as_mut_ptr());
    let nogroups = s[3];
    let gotbuoyancy = s[15];
    let mut usercoeff: c_double = 0.0;
    if st.particle_to_field {
        pda_f_reset(
            dtime,
            if st.field_reset { 1 } else { 0 },
            &mut usercoeff,
        );
    }
    let mut nop: c_int = 0;
    let mut nogrp: c_int = 0;
    pda_f_counts(&mut nop, &mut nogrp);
    pda_f_debug();
    let (mut mass, mut damping, mut rad, mut charge, mut drag): (
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
    ) = (0.0, 0.0, 0.0, 0.0, 0.0);
    if nogrp > 1 {
        pda_f_groups_fetch();
    } else {
        let (mut gm, mut gd, mut gr, mut gc, mut gdrag): (c_int, c_int, c_int, c_int, c_int) =
            (0, 0, 0, 0, 0);
        pda_f_scalars(
            &mut mass,
            &mut damping,
            &mut rad,
            &mut charge,
            &mut drag,
            &mut gm,
            &mut gd,
            &mut gr,
            &mut gc,
            &mut gdrag,
        );
    }
    pda_f_checks();
    if gotbuoyancy != 0 {
        pda_f_buoyvol(rad);
    }
    let mut no = 1;
    while no <= nop {
        let mut status: c_int = 0;
        let mut elemidx: c_int = 0;
        let mut group: c_int = 0;
        let mut proceed: c_int = 0;
        pda_f_begin_no(no, &mut status, &mut elemidx, &mut group, &mut proceed);
        if proceed == 0 {
            no += 1;
            continue;
        }
        if nogroups > 1 {
            pda_f_group_update(
                group,
                &mut mass,
                &mut damping,
                &mut rad,
                &mut charge,
                &mut drag,
            );
        }
        let mut fstatus = status;
        let mut fproceed: c_int = 0;
        pda_f_force(
            dtime,
            mass,
            damping,
            charge,
            drag,
            rad,
            &mut fstatus,
            &mut fproceed,
        );
        if fproceed == 0 {
            no += 1;
            continue;
        }
        status = fstatus;
        if set_fields {
            let g = if nogroups > 1 { group } else { 0 };
            pda_f_accum(no, dtime, st.tottime, usercoeff, mass, charge, g);
        }
        if set_particles {
            pda_f_setparticles(no, status);
        }
        no += 1;
    }
    if set_fields {
        pda_f_normalize();
    }
    pda_f_set_prevdtime(dtime);
}

/// The `ParticleWallContact` internal procedure: spring contact and
/// elastic/inelastic collision handling at generic wall boundaries,
/// with optional accumulation, wall trace and moving-wall velocity.
/// Mesh topology, distances and forces stay in Fortran wrappers;
/// Rust drives the particle/face loops and holds the first-visit
/// sizing plus the `AnyInteraction` gate.
unsafe fn particle_wall_contact(dtime: f64, h_model: *mut c_void) {
    let mut pwc = PWC.lock().unwrap();
    if !pwc.visited {
        let mut dim: c_int = 0;
        let mut anyinter: c_int = 0;
        pda_t_first(h_model, &mut dim, &mut anyinter);
        pwc.dim = dim;
        pwc.any_interaction = anyinter != 0;
        pwc.visited = true;
    }
    if !pwc.any_interaction {
        return;
    }
    let dim = pwc.dim;
    drop(pwc);
    let n = pda_m_nparticles();
    let mut no = 1;
    while no <= n {
        let mut status: c_int = 0;
        let mut elemidx: c_int = 0;
        let mut proceed: c_int = 0;
        pda_t_begin_no(no, &mut status, &mut elemidx, &mut proceed);
        if proceed == 0 {
            no += 1;
            continue;
        }
        let imax = pda_t_imax();
        let mut face = 1;
        let mut hit: c_int = 0;
        while face <= imax {
            let mut found: c_int = 0;
            let mut collision: c_int = 0;
            let mut contact: c_int = 0;
            pda_t_face(face, &mut found, &mut collision, &mut contact);
            if found == 0 {
                face += 1;
                continue;
            }
            if collision == 0 && contact == 0 {
                face += 1;
                continue;
            }
            let mut coord = [0.0 as c_double; 3];
            let mut normal = [0.0 as c_double; 3];
            let mut dist: c_double = 0.0;
            pda_t_geom(
                no,
                dim,
                coord.as_mut_ptr(),
                normal.as_mut_ptr(),
                &mut dist,
            );
            if contact != 0 {
                let mut h: c_int = 0;
                let mut e: c_int = 0;
                pda_t_contact(no, dtime, dim, dist, &mut h, &mut e);
                if h != 0 {
                    hit = 1;
                }
                if e != 0 {
                    break;
                }
            }
            if collision != 0 {
                let mut h: c_int = 0;
                let mut e: c_int = 0;
                let mut skip: c_int = 0;
                let mut wt: c_int = 0;
                pda_t_collision(no, dtime, dim, dist, &mut h, &mut e, &mut skip, &mut wt);
                if h != 0 {
                    hit = 1;
                }
                if skip != 0 {
                    face += 1;
                    continue;
                }
                if e != 0 {
                    break;
                }
            }
            face += 1;
        }
        if hit == 0 {
            no += 1;
            continue;
        }
        let mut ws = [0 as c_int; 6];
        pda_t_state(ws.as_mut_ptr());
        if ws[5] != 0 {
            pda_t_trace(no, dim);
        }
        no += 1;
    }
}

/// `ParticleParticleContact` — the granular contact kernel re-entered
/// by `ParticleParticleInteraction` through the plain-Fortran
/// trampoline. First-visit parameters mirror the module SAVE with
/// identical lifetimes; the force law matches the original statement
/// for statement.
///
/// # Safety
///
/// All array pointers must address at least three `c_double`s; the
/// contact flag must be a valid write.
#[export_name = "__particledynamicsstuff_MOD_particleparticlecontact"]
pub unsafe extern "C" fn particle_particle_contact(
    dt: *const c_double,
    coord: *mut c_double,
    coord2: *mut c_double,
    velo: *mut c_double,
    velo2: *mut c_double,
    force: *mut c_double,
    force2: *mut c_double,
    contact: *mut c_int,
) {
    let _ = dt;
    let mut st = PPC.lock().unwrap();
    if !st.visited {
        let mut rad: c_double = 0.0;
        let mut mass: c_double = 0.0;
        let mut spring: c_double = 0.0;
        let mut damping: c_double = 0.0;
        let mut friction: c_double = 0.0;
        pda_c_first(
            &mut rad,
            &mut mass,
            &mut spring,
            &mut damping,
            &mut friction,
        );
        st.rad = rad;
        st.mass = mass;
        st.spring = spring;
        st.damping = damping;
        st.friction = friction;
        st.visited = true;
    }
    let (rad, spring) = (st.rad, st.spring);
    drop(st);
    *contact = 0;
    let c = std::slice::from_raw_parts(coord as *const c_double, 3);
    let c2 = std::slice::from_raw_parts(coord2 as *const c_double, 3);
    let v = std::slice::from_raw_parts(velo as *const c_double, 3);
    let v2 = std::slice::from_raw_parts(velo2 as *const c_double, 3);
    let mut dr = [0.0 as c_double; 3];
    let mut k = 0;
    while k < 3 {
        dr[k] = c[k] - c2[k];
        k += 1;
    }
    let dist = pda_c_dot(dr.as_ptr(), dr.as_ptr()).sqrt();
    if dist < f64::MIN_POSITIVE {
        pda_c_warn_samepoint();
        return;
    }
    let eta = 2.0 * rad - dist;
    if eta < 0.0 {
        return;
    }
    let rn = [dr[0] / dist, dr[1] / dist, dr[2] / dist];
    let mut dv = [0.0 as c_double; 3];
    k = 0;
    while k < 3 {
        dv[k] = v[k] - v2[k];
        k += 1;
    }
    let speed = pda_c_dot(dv.as_ptr(), dv.as_ptr()).sqrt();
    let vn = if speed > f64::MIN_POSITIVE {
        [dv[0] / speed, dv[1] / speed, dv[2] / speed]
    } else {
        [1.0, 0.0, 0.0]
    };
    let spring_force = eta * spring;
    let damp_force = 0.0;
    let f = std::slice::from_raw_parts_mut(force, 3);
    let f2 = std::slice::from_raw_parts_mut(force2, 3);
    k = 0;
    while k < 3 {
        f[k] = spring_force * rn[k] - damp_force * vn[k];
        f2[k] = -f[k];
        k += 1;
    }
    *contact = 1;
}

/// `ParticleParticleCollision` — the binary-collision kernel
/// re-entered by `ParticleParticleInteraction` through the
/// plain-Fortran trampoline. Statement-for-statement mirror of the
/// original, including the elementwise triple-product position
/// updates and the force/velocity branches.
///
/// # Safety
///
/// All array pointers must address at least three `c_double`s; the
/// collision flag must be a valid write.
#[export_name = "__particledynamicsstuff_MOD_particleparticlecollision"]
pub unsafe extern "C" fn particle_particle_collision(
    dt: *const c_double,
    coord: *mut c_double,
    coord2: *mut c_double,
    velo: *mut c_double,
    velo2: *mut c_double,
    force: *mut c_double,
    force2: *mut c_double,
    collision: *mut c_int,
) {
    let dtv = *dt;
    let mut st = PPK.lock().unwrap();
    if !st.visited {
        let mut rad1: c_double = 0.0;
        let mut mass1: c_double = 0.0;
        let mut coeff: c_double = 0.0;
        let mut tc: c_int = 0;
        pda_k_first(&mut rad1, &mut mass1, &mut coeff, &mut tc);
        st.rad1 = rad1;
        st.mass1 = mass1;
        st.coeff = coeff;
        st.true_collision = tc != 0;
        st.visited = true;
    }
    let (rad1, mass1, coeff, true_collision) = (st.rad1, st.mass1, st.coeff, st.true_collision);
    drop(st);
    *collision = 0;
    let c = std::slice::from_raw_parts(coord as *const c_double, 3);
    let c2 = std::slice::from_raw_parts(coord2 as *const c_double, 3);
    let v = std::slice::from_raw_parts(velo as *const c_double, 3);
    let v2 = std::slice::from_raw_parts(velo2 as *const c_double, 3);
    let mut dr = [0.0 as c_double; 3];
    let mut dv = [0.0 as c_double; 3];
    let mut k = 0;
    while k < 3 {
        dr[k] = c[k] - c2[k];
        dv[k] = v[k] - v2[k];
        k += 1;
    }
    let b = pda_k_dot(dr.as_ptr(), dv.as_ptr());
    if b >= 0.0 {
        return;
    }
    let a = pda_k_dot(dv.as_ptr(), dv.as_ptr());
    let cc = pda_k_dot(dr.as_ptr(), dr.as_ptr()) - 4.0 * rad1 * rad1;
    let d = b * b - a * cc;
    if d < 0.0 {
        return;
    }
    let dta = (-b - d.sqrt()) / a;
    if dta >= dtv {
        return;
    }
    let dtb = dtv - dta;
    let mut dra = [0.0 as c_double; 3];
    k = 0;
    while k < 3 {
        dra[k] = dr[k] + dta * dv[k];
        k += 1;
    }
    let dra_n = pda_k_dot(dra.as_ptr(), dra.as_ptr()).sqrt();
    let rn = [dra[0] / dra_n, dra[1] / dra_n, dra[2] / dra_n];
    let v1na = v[0] * rn[0] + v[1] * rn[1] + v[2] * rn[2];
    let v2na = v2[0] * rn[0] + v2[1] * rn[1] + v2[2] * rn[2];
    let v1nb = (coeff * (v2na - v1na) + v1na + v2na) / 2.0;
    let v2nb = (coeff * (v1na - v2na) + v2na + v1na) / 2.0;
    let co = std::slice::from_raw_parts_mut(coord, 3);
    let co2 = std::slice::from_raw_parts_mut(coord2, 3);
    let ve = std::slice::from_raw_parts_mut(velo, 3);
    let ve2 = std::slice::from_raw_parts_mut(velo2, 3);
    let fo = std::slice::from_raw_parts_mut(force, 3);
    let fo2 = std::slice::from_raw_parts_mut(force2, 3);
    if true_collision {
        k = 0;
        while k < 3 {
            co[k] += dta * rn[k] * ve[k];
            co2[k] += dta * rn[k] * ve2[k];
            k += 1;
        }
        k = 0;
        while k < 3 {
            ve[k] += (v1nb - v1na) * rn[k];
            ve2[k] += (v2nb - v2na) * rn[k];
            k += 1;
        }
        k = 0;
        while k < 3 {
            co[k] += dtb * ve[k];
            co2[k] += dtb * ve2[k];
            k += 1;
        }
    } else {
        k = 0;
        while k < 3 {
            co[k] += (v1na - v1nb) * rn[k] * dta;
            co2[k] += (v2na - v2nb) * rn[k] * dta;
            fo[k] = mass1 * (v1nb - v1na) * rn[k] / dtv;
            fo2[k] = mass1 * (v2nb - v2na) * rn[k] / dtv;
            k += 1;
        }
    }
    *collision = 1;
}

/// `ParticleWallProc` — the wall kernel re-entered by
/// `LocateParticles` through the plain-Fortran trampoline. Mirrors
/// the original: advance to the collision point, reflect, reorient,
/// and the optional photon-bunch energy deposit.
///
/// # Safety
///
/// Array pointers must address at least three `c_double`s; integer
/// pointers must be valid reads/writes.
#[export_name = "__particledynamicsstuff_MOD_particlewallproc"]
pub unsafe extern "C" fn particle_wall_proc(
    no: *const c_int,
    rinit: *mut c_double,
    rfin: *mut c_double,
    vinit: *mut c_double,
    vfin: *mut c_double,
    lambda: *const c_double,
    faceindex: *const c_int,
    pstatus: *mut c_int,
) {
    let _ = pstatus;
    let mut st = PWP.lock().unwrap();
    if !st.visited {
        let mut meshdim: c_int = 0;
        let mut bunch: c_int = 0;
        pda_w_first(&mut meshdim, &mut bunch);
        st.mesh_dim = meshdim;
        st.particle_bunch = bunch != 0;
        st.visited = true;
    }
    let bunch = st.particle_bunch;
    drop(st);
    let fi = *faceindex;
    let lam = *lambda;
    let n = *no;
    pda_w_face(fi);
    pda_w_advance(lam, rinit, rfin);
    pda_w_reflect(rinit, rfin, vinit, vfin);
    if !bunch {
        return;
    }
    pda_w_bunch(n, rinit);
}

/// `particledynamics_init_` — solver init entry: registers the
/// no-output global variable, exactly as the original.
///
/// # Safety
///
/// `solver` must be the Elmer `Solver_t` address the loader passes.
#[export_name = "particledynamics_init_"]
pub unsafe extern "C" fn particle_dynamics_init_entry(
    _model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    pda_i_init(solver as *const c_void);
}





