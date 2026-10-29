//! Port of `fem/src/modules/RigidMeshMapper.F90` (Elmer @ c0b4bc6, the
//! campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/RigidMeshMapper.F90` — one F77-global
//! entry `rigidmeshmapper_` with one `CONTAINS` internal procedure
//! (`LocalMatrix`). The emitted dynamic symbol is exactly
//! `rigidmeshmapper_` (proven by `nm -D` on the reference `.so`); the
//! internal emits only a local symbol.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in `accessors/rigidmeshmapper_acc.F90`,
//! compiled against the built `.mod` files of the target configure — the
//! port itself contains no layout assumption and no pinned offset. Order-
//! sensitive reductions (`MATMUL`, `SUM`), formatted I/O (`WRITE`), and
//! whole-array updates stay in Fortran wrappers, so bit-identity holds
//! by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors RigidMeshMapper.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `Visited` + `GotRelaxField` + the
//!   `RelaxField`/`RelaxPerm` pointers and the `RotorElement` marking
//!   live in the guarded `static` state with identical lifetimes; the
//!   internal procedure's `SAVE Nodes`/`Hwrk` live as wrapper-side SAVE
//!   state. `STIFF`/`FORCE` are per-call scratch (as in the original),
//!   sized by `Mesh % MaxElementNodes`, with allocation failure calling
//!   the same `Fatal(...)` via `try_reserve` (the `STAT=istat` check).
//! - The original recomputes `GotMatrix = GotRotate .OR. GotScale` after
//!   the `Mesh Matrix` fast path, so a lone `Mesh Matrix` array is a
//!   (defined) per-node no-op; the port replicates the reassignment.
//! - With `Mesh Scale` but no rotation the original reads its
//!   uninitialized `RotMatrix` local (undefined behavior); the port
//!   keeps the same wrapper-state shape and its proofs keep to defined
//!   inputs (scale only together with rotation), documented in the proof
//!   script. Likewise per-component `Mesh Rotate/Scale/Translate/Origin`
//!   keywords are specified completely in the proofs, as upstream tests
//!   do.
//! - `Model`, `dt`, `TransientSimulation` travel like the original:
//!   `Model` selects the target mesh, `dt` feeds the mesh-velocity
//!   update, `TransientSimulation` is unused (received and ignored).
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn rmm_begin(h_model: *const c_void, h_solver: *const c_void);
    fn rmm_solver_params();
    fn rmm_target_mesh_index(k: *mut c_int, found: *mut c_int);
    fn rmm_select_mesh(use_target: c_int, k: c_int);
    fn rmm_assoc_mesh() -> c_int;
    fn rmm_fatal_no_mesh();
    fn rmm_coord_dim() -> c_int;
    fn rmm_get_logical(which: c_int, found: *mut c_int) -> c_int;
    fn rmm_list_logical(which: c_int) -> c_int;
    fn rmm_angle_coeff() -> c_double;
    fn rmm_no_nodes() -> c_int;
    fn rmm_num_bulk() -> c_int;
    fn rmm_max_elem_nodes() -> c_int;
    fn rmm_set_use_orig(flag: c_int);
    fn rmm_rotor_bodies() -> c_int;
    fn rmm_rotor_radius() -> c_double;
    fn rmm_rotor_angle(a: *mut c_double, got: *mut c_int);
    fn rmm_set_mesh_skew();
    fn rmm_rotor_mark(elem: c_int, n: *mut c_int, isrotor: *mut c_int);
    fn rmm_store_orig();
    fn rmm_velocity_setup(dim: c_int);
    fn rmm_velocity_node(nodei: c_int, dx: *const c_double, dt: c_double, dim: c_int);
    fn rmm_rotorder_fetch(j: *mut c_int, found: *mut c_int);
    fn rmm_rotorder_val(i: c_int) -> c_int;
    fn rmm_set_rotorder(a: c_int, b: c_int, c: c_int);
    fn rmm_fatal_rotaxis();
    fn rmm_assoc_matrix() -> c_int;
    fn rmm_relax_bcs();
    fn rmm_nonlin_maxiter() -> c_int;
    fn rmm_nonlin_coeff(coeff: *mut c_double, source: *mut c_double);
    fn rmm_default_initialize();
    fn rmm_nof_active() -> c_int;
    fn rmm_begin_active(t: c_int, n: *mut c_int);
    fn rmm_default_update(ld: c_int, g: *mut c_double, f: *mut c_double);
    fn rmm_default_finish_bulk();
    fn rmm_default_finish_assembly();
    fn rmm_default_dirichlet();
    fn rmm_default_solve() -> c_double;
    fn rmm_nonlin_converged() -> c_int;
    fn rmm_relax_normalize();
    fn rmm_relax_fetch() -> c_int;
    fn rmm_fatal_alloc();
    fn rmm_lm_nodes();
    fn rmm_lm_localsol(n: c_int, arr: *mut c_double);
    fn rmm_lm_youngs_begin(n: c_int, got: *mut c_int, k: *mut c_int);
    fn rmm_lm_elemyoungs(n: c_int, k: c_int, ey: *mut c_double);
    fn rmm_lm_ip_begin() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn rmm_lm_ip_point(t: c_int, u: *mut c_double, v: *mut c_double, w: *mut c_double, s: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn rmm_lm_eleminfo(
        u: c_double,
        v: c_double,
        w: c_double,
        n: c_int,
        basis: *mut c_double,
        dbasis: *mut c_double,
        detj: *mut c_double,
        stat: *mut c_int,
    );
    fn rmm_lm_grad(n: c_int, dbasis: *const c_double, localsol: *const c_double, grad: *mut c_double);
    fn rmm_lm_youngs_ip(
        n: c_int,
        k: c_int,
        got: c_int,
        basis: *const c_double,
        ey: *const c_double,
        youngs: *mut c_double,
    );
    #[allow(clippy::too_many_arguments)]
    fn rmm_lm_assemble(
        n: c_int,
        grad: *const c_double,
        youngs: *const c_double,
        coeff: c_double,
        s: c_double,
        detj: c_double,
        source: c_double,
        basis: *const c_double,
        dbasis: *const c_double,
        ld: c_int,
        stiff: *mut c_double,
        force: *mut c_double,
    );
    fn rmm_obsolete_displacement() -> c_int;
    fn rmm_anymesh(mode: c_int, which: c_int) -> c_int;
    fn rmm_use_solverparams_vlist();
    fn rmm_vlist_reset();
    fn rmm_elem_bodyforce(bodyid: c_int, bf_id: *mut c_int, found: *mut c_int);
    fn rmm_constmap_step() -> c_int;
    fn rmm_begin_bulk(elem: c_int, n: *mut c_int);
    fn rmm_elem_bodyid() -> c_int;
    fn rmm_node_index(pos: c_int, nodei: *mut c_int);
    fn rmm_mesh_matrix() -> c_int;
    fn rmm_mesh_rotate(coeff: c_double) -> c_int;
    fn rmm_rotor_override(coeff: c_double, rotorangle: c_double);
    fn rmm_mesh_scale() -> c_int;
    fn rmm_mesh_translate() -> c_int;
    fn rmm_build_rot(gotrotate: c_int);
    fn rmm_build_trf(gottranslate: c_int, tbr: c_int, gotscale: c_int);
    fn rmm_mesh_origin(gotmatrix: c_int, anymeshorigin: c_int);
    fn rmm_apply_node(
        nodei: c_int,
        gotmatrix: c_int,
        gottranslate: c_int,
        dx: *mut c_double,
        cycled: *mut c_int,
    );
    fn rmm_relax_node(nodei: c_int, dx: *mut c_double, gotrelax: c_int, anyrelax: c_int);
    fn rmm_write_node(nodei: c_int, dx: *const c_double);
    fn rmm_info_start();
    fn rmm_info_targetmesh(k: c_int);
    fn rmm_info_movewhole();
    fn rmm_info_solvecaller();
    fn rmm_info_laplace();
    fn rmm_info_mapped(count: c_int);
    fn rmm_cputime() -> c_double;
    fn rmm_info_time(at1: c_double, at0: c_double);
    fn rmm_info_done();
    fn rmm_fatal_rotor();
    fn rmm_default_finish();
}

/// First-call sizing is per-call here (the original allocates
/// `STIFF`/`FORCE` fresh on every `DoIt` pass); growth handling is
/// unnecessary because the scratch is sized by `MaxElementNodes`, the
/// same bound the original uses.
struct MapperState {
    visited: bool,
    got_relax_field: bool,
    rotor: Option<Vec<bool>>,
}

impl MapperState {
    const fn new() -> Self {
        Self {
            visited: false,
            got_relax_field: false,
            rotor: None,
        }
    }
}

static STATE: Mutex<MapperState> = Mutex::new(MapperState::new());

/// `rigidmeshmapper_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dt, TransientSimulation)`; `TransientSimulation` is
/// unused by the body, exactly as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer addresses the loader passes.
/// `dt` must be a valid read (it feeds the mesh-velocity update).
#[export_name = "rigidmeshmapper_"]
pub unsafe extern "C" fn rigidmeshmapper_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    let dt_val: c_double = unsafe { *dt };
    unsafe { rmm_begin(model, solver) };

    unsafe { rmm_info_start() };
    unsafe { rmm_solver_params() };

    // Target mesh selection.
    let (mut k, mut found) = (0, 0);
    unsafe { rmm_target_mesh_index(&mut k, &mut found) };
    if found != 0 {
        unsafe { rmm_info_targetmesh(k) };
        unsafe { rmm_select_mesh(1, k) };
    } else {
        unsafe { rmm_select_mesh(0, 0) };
    }
    if unsafe { rmm_assoc_mesh() } == 0 {
        unsafe { rmm_fatal_no_mesh() };
    }

    let dim = unsafe { rmm_coord_dim() };

    let cumulative = unsafe { rmm_get_logical(1, &mut found) } != 0;
    let use_original_mesh = !cumulative;
    let store_original = unsafe { rmm_get_logical(2, &mut found) } != 0;

    let angle_coeff = unsafe { rmm_angle_coeff() };

    let whole_mode = unsafe { rmm_list_logical(1) } != 0;
    if whole_mode {
        unsafe { rmm_info_movewhole() };
    }

    let rotor_mode = unsafe { rmm_list_logical(2) } != 0;
    let mut got_rotor_angle = 0;
    let mut rotor_angle = 0.0;
    if rotor_mode && unsafe { rmm_rotor_bodies() } == 0 {
        // Fatals inside when neither radius nor bodies are given.
        let _ = unsafe { rmm_rotor_radius() };
    }
    unsafe { rmm_rotor_angle(&mut rotor_angle, &mut got_rotor_angle) };

    let no_nodes = unsafe { rmm_no_nodes() } as usize;
    let mut node_done = vec![false; no_nodes];

    if rotor_mode && !state.visited {
        // Skew is moved to the library so that it can be done in
        // different ways.
        unsafe { rmm_set_mesh_skew() };
        let num_bulk = unsafe { rmm_num_bulk() };
        let mut marking = vec![false; num_bulk as usize];
        for elem in 1..=num_bulk {
            let (mut n, mut isrotor) = (0, 0);
            unsafe { rmm_rotor_mark(elem, &mut n, &mut isrotor) };
            let _ = n;
            marking[(elem - 1) as usize] = isrotor != 0;
        }
        state.rotor = Some(marking);
        node_done.fill(false);
    }

    // If using original mesh as a reference mesh it must be saved,
    // otherwise the analytic mapping does not require two meshes.
    if !state.visited && (use_original_mesh || store_original) {
        unsafe { rmm_store_orig() };
    }

    unsafe { rmm_set_use_orig(i32::from(use_original_mesh)) };

    let calculate_velocity = unsafe { rmm_get_logical(4, &mut found) } != 0;
    if calculate_velocity {
        unsafe { rmm_velocity_setup(dim) };
    }

    let translate_before_rotate = unsafe { rmm_get_logical(3, &mut found) } != 0;

    // Rotation order (reversed, exactly as in the original).
    let mut rotate_order = [1 as c_int, 2, 3];
    {
        let (mut j, mut ofound) = (0, 0);
        unsafe { rmm_rotorder_fetch(&mut j, &mut ofound) };
        if ofound != 0 {
            if j != dim {
                unsafe { rmm_fatal_rotaxis() };
            }
            for i in 1..=j {
                rotate_order[(i - 1) as usize] =
                    unsafe { rmm_rotorder_val(j + 1 - i) };
            }
        }
        unsafe { rmm_set_rotorder(rotate_order[0], rotate_order[1], rotate_order[2]) };
    }

    let mut doit = unsafe { rmm_assoc_matrix() } != 0;
    if doit {
        doit = !state.visited || unsafe { rmm_list_logical(4) } != 0;
    }

    if doit {
        let nmax = unsafe { rmm_max_elem_nodes() } as usize;
        let mut stiff = Vec::<c_double>::new();
        let mut force = Vec::<c_double>::new();
        if stiff.try_reserve(nmax * nmax).is_err() || force.try_reserve(nmax).is_err() {
            unsafe { rmm_fatal_alloc() };
            std::process::abort();
        }
        stiff.resize(nmax * nmax, 0.0);
        force.resize(nmax, 0.0);

        unsafe { rmm_info_solvecaller() };

        // Implement moving and fixed BCs.
        unsafe { rmm_relax_bcs() };

        unsafe { rmm_info_laplace() };

        let max_nonlin_iter = unsafe { rmm_nonlin_maxiter() };
        let (mut coeff, mut source) = (0.0, 0.0);
        unsafe { rmm_nonlin_coeff(&mut coeff, &mut source) };

        for _nonlin in 1..=max_nonlin_iter {
            unsafe { rmm_default_initialize() };

            let active = unsafe { rmm_nof_active() };
            for t in 1..=active {
                let mut n = 0;
                unsafe { rmm_begin_active(t, &mut n) };
                local_matrix(
                    n,
                    nmax as c_int,
                    stiff.as_mut_ptr(),
                    force.as_mut_ptr(),
                    coeff,
                    source,
                );
                unsafe {
                    rmm_default_update(
                        nmax as c_int,
                        stiff.as_mut_ptr(),
                        force.as_mut_ptr(),
                    )
                };
            }
            unsafe { rmm_default_finish_bulk() };

            // No Flux BCs.
            unsafe { rmm_default_finish_assembly() };
            unsafe { rmm_default_dirichlet() };
            let _norm: c_double = unsafe { rmm_default_solve() };

            if unsafe { rmm_nonlin_converged() } == 1 {
                break;
            }
        }

        if unsafe { rmm_list_logical(5) } != 0 {
            unsafe { rmm_relax_normalize() };
        }

        if unsafe { rmm_relax_fetch() } != 0 {
            state.got_relax_field = true;
        }
    }

    // Initialize the mapping matrices (Identity kept wrapper-side).
    let at0 = unsafe { rmm_cputime() };

    // Obsolete-keyword warn (the flag it returns is unconditionally
    // recomputed below, exactly as in the original).
    unsafe { rmm_obsolete_displacement() };

    // ValueList selection mode for the presence checks.
    let check_mode: c_int = if rotor_mode || whole_mode { 1 } else { 0 };
    if rotor_mode || whole_mode {
        unsafe { rmm_use_solverparams_vlist() };
    }
    let any_mesh_matrix = unsafe { rmm_anymesh(check_mode, 1) } != 0;
    let any_mesh_rotate = unsafe { rmm_anymesh(check_mode, 2) } != 0;
    // The original recomputes this in the mode branch unconditionally.
    let any_mesh_translate = unsafe { rmm_anymesh(check_mode, 3) } != 0;
    let any_mesh_scale = unsafe { rmm_anymesh(check_mode, 4) } != 0;
    let any_mesh_origin = unsafe { rmm_anymesh(check_mode, 5) } != 0;
    let any_relax = unsafe { rmm_anymesh(check_mode, 6) } != 0;

    // Unconditional initialization, exactly as in the original
    // (GotRotate/GotTranslate/GotScale/GotMatrix = .FALSE.).
    #[allow(unused_assignments)]
    let (mut got_rotate, mut got_translate, mut got_scale, mut got_matrix) =
        (false, false, false, false);

    unsafe { rmm_vlist_reset() };
    let mut got_map = false;
    let constant_map = unsafe { rmm_list_logical(3) } != 0;

    let num_bulk = unsafe { rmm_num_bulk() };
    for elem in 1..=num_bulk {
        let mut elem_nodes = 0;
        unsafe { rmm_begin_bulk(elem, &mut elem_nodes) };

        if whole_mode {
            // If we are doing the whole mesh then do all elements.
        } else if rotor_mode {
            let marking = match state.rotor.as_ref() {
                Some(m) => m,
                None => {
                    unsafe { rmm_fatal_rotor() };
                    std::process::abort();
                }
            };
            if !marking[(elem - 1) as usize] {
                continue;
            }
        } else {
            let body_id = unsafe { rmm_elem_bodyid() };
            let (mut bf_id, mut bfound) = (0, 0);
            unsafe { rmm_elem_bodyforce(body_id, &mut bf_id, &mut bfound) };
            if bfound == 0 {
                continue;
            }
            if constant_map {
                got_map = unsafe { rmm_constmap_step() } != 0;
            }
        }

        for pos in 1..=elem_nodes {
            let mut node_i = 0;
            unsafe { rmm_node_index(pos, &mut node_i) };
            let idx = (node_i - 1) as usize;
            if node_done[idx] {
                continue;
            }

            // Save time: identical mapping reuses the previous lookup.
            if !got_map {
                if rotor_mode || whole_mode {
                    got_map = true;
                }

                // Generic transformation matrix (`IF(AnyMeshMatrix)` gate).
                got_matrix =
                    any_mesh_matrix && unsafe { rmm_mesh_matrix() } != 0;

                if !got_matrix {
                    // Rotations around main axis.
                    got_rotate = false;
                    if any_mesh_rotate {
                        got_rotate = unsafe { rmm_mesh_rotate(angle_coeff) } != 0;
                    }
                    if got_rotor_angle != 0 {
                        got_rotate = true;
                        unsafe { rmm_rotor_override(angle_coeff, rotor_angle) };
                    }

                    // Scaling.
                    got_scale = false;
                    if any_mesh_scale {
                        got_scale = unsafe { rmm_mesh_scale() } != 0;
                    }

                    // Translations.
                    got_translate = false;
                    if any_mesh_translate {
                        got_translate = unsafe { rmm_mesh_translate() } != 0;
                    }

                    // NOTE: mirrors the original's unconditional
                    // reassignment, including its consequence for a
                    // lone Mesh Matrix array (defined no-op below).
                    got_matrix = got_rotate || got_scale;

                    if got_matrix {
                        unsafe { rmm_build_rot(i32::from(got_rotate)) };
                        unsafe {
                            rmm_build_trf(
                                i32::from(got_translate),
                                i32::from(translate_before_rotate),
                                i32::from(got_scale),
                            )
                        };
                    }
                }

                // Mesh origin.
                unsafe {
                    rmm_mesh_origin(i32::from(got_matrix), i32::from(any_mesh_origin))
                };
            }

            let mut dx = [0.0 as c_double; 3];
            let mut cycled = 0;
            unsafe {
                rmm_apply_node(
                    node_i,
                    i32::from(got_matrix),
                    i32::from(got_translate),
                    dx.as_mut_ptr(),
                    &mut cycled,
                )
            };
            if cycled != 0 {
                continue;
            }

            // Relaxation interpolation between moving and fixed walls.
            unsafe {
                rmm_relax_node(
                    node_i,
                    dx.as_mut_ptr(),
                    i32::from(state.got_relax_field),
                    i32::from(any_relax),
                )
            };

            // Mesh velocity.
            if calculate_velocity {
                unsafe { rmm_velocity_node(node_i, dx.as_ptr(), dt_val, dim) };
            }

            unsafe { rmm_write_node(node_i, dx.as_ptr()) };

            node_done[idx] = true;
        }
    }

    if !state.visited {
        let count = node_done.iter().filter(|v| **v).count() as c_int;
        unsafe { rmm_info_mapped(count) };
        let at1 = unsafe { rmm_cputime() };
        unsafe { rmm_info_time(at1, at0) };
        unsafe { rmm_info_done() };
    }

    drop(node_done);

    unsafe { rmm_default_finish() };

    state.visited = true;
}

/// The `LocalMatrix` internal procedure: element matrix and right-hand
/// side by numerical integration. `n` is the element node count;
/// `ld` is the `MaxElementNodes` stride of the shared scratch (the
/// wrapper updates the `(1:n, 1:n)` leading section, never a
/// re-strided `(n, n)` view). `coeff`/`source` are the host
/// `Nonlinear Conductivity Coefficient` / `Mesh Relax Source` values.
#[allow(clippy::too_many_arguments)]
fn local_matrix(
    n: c_int,
    ld: c_int,
    stiff: *mut c_double,
    force: *mut c_double,
    coeff: c_double,
    source: c_double,
) {
    let nn = n as usize;
    unsafe { rmm_lm_nodes() };
    // STIFF = 0 / FORCE = 0 whole-array fills.
    unsafe {
        std::slice::from_raw_parts_mut(stiff, (ld as usize) * (ld as usize)).fill(0.0);
        std::slice::from_raw_parts_mut(force, ld as usize).fill(0.0);
    }

    let mut localsol = vec![0.0 as c_double; nn];
    unsafe { rmm_lm_localsol(n, localsol.as_mut_ptr()) };

    let (mut got_youngs, mut k) = (0, 0);
    unsafe { rmm_lm_youngs_begin(n, &mut got_youngs, &mut k) };
    // ElemYoungs(3, n) automatics, sized by this element.
    let mut elemyoungs = vec![0.0 as c_double; 3 * nn];
    if got_youngs != 0 {
        unsafe { rmm_lm_elemyoungs(n, k, elemyoungs.as_mut_ptr()) };
    }

    let mut basis = vec![0.0 as c_double; nn];
    let mut dbasis = vec![0.0 as c_double; nn * 3];

    let nip = unsafe { rmm_lm_ip_begin() };
    for t in 1..=nip {
        let (mut u, mut v, mut w, mut s) = (0.0, 0.0, 0.0, 0.0);
        unsafe { rmm_lm_ip_point(t, &mut u, &mut v, &mut w, &mut s) };

        let (mut detj, mut stat) = (0.0, 0);
        unsafe {
            rmm_lm_eleminfo(
                u,
                v,
                w,
                n,
                basis.as_mut_ptr(),
                dbasis.as_mut_ptr(),
                &mut detj,
                &mut stat,
            )
        };
        let _ = stat;

        let mut grad = [0.0 as c_double; 3];
        unsafe {
            rmm_lm_grad(
                n,
                dbasis.as_ptr(),
                localsol.as_ptr(),
                grad.as_mut_ptr(),
            )
        };

        let mut youngs = [0.0 as c_double; 3];
        unsafe {
            rmm_lm_youngs_ip(
                n,
                k,
                got_youngs,
                basis.as_ptr(),
                elemyoungs.as_ptr(),
                youngs.as_mut_ptr(),
            )
        };

        unsafe {
            rmm_lm_assemble(
                n,
                grad.as_ptr(),
                youngs.as_ptr(),
                coeff,
                s,
                detj,
                source,
                basis.as_ptr(),
                dbasis.as_ptr(),
                ld,
                stiff,
                force,
            )
        };
    }
}
