//! Port of `fem/src/modules/StructuredMeshMapper.F90` (Elmer @ c0b4bc6,
//! the campaign fork) for lane T021-B, at the solver-plugin boundary (whole
//! `.so` → Rust cdylib with the same dynamic symbols).
//!
//! Unit: `fortran:fem/src/modules/StructuredMeshMapper.F90` — two F77-global
//! entries `structuredmeshmapper_` + `structuredmeshmapper_init_`. The
//! emitted dynamic symbols are exactly those two (proven by `nm -D` on the
//! reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entries' derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/structuredmeshmapper_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive reductions
//! (`MATMUL` in the multilayer projector, `ParallelReduction` over the
//! tangled/limited counters), formatted I/O (`PRINT`/`WRITE`, `Message`
//! construction), and whole-array updates (OrigCoord/Coord sync, BC field
//! gathers, stride commit/WHERE blocks) stay in Fortran wrappers, so
//! bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors StructuredMeshMapper.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: plain SAVE scalars (`Visited`,
//!   `Initialized`, `TopMode`, `BotMode`, `nsize`, `nnodes`, `MaskExists`,
//!   `ComputeTangledMask`, `MidLayerExists`, `NumberOfLayers`,
//!   `RecompStabExe`, `Cnt`) live in the Rust `Mutex`-guarded state with
//!   identical lifetimes; Fortran-object SAVE state (pointer associations
//!   `Coord`/`TopPointer`/`BotPointer`/`MidPointer`/`UpPointer`/
//!   `DownPointer`/`NodeLayer`/`MaskPerm`/`TopField`/..., `ALLOCATABLE`
//!   `OrigCoord`/`Field`/`Surface`, variable handles) lives wrapper-side
//!   with identical lifetimes. Per-call locals (`Reinitialize`,
//!   `RecompStab`, `GotVeloVar`, ..., `TangledCount`) are recomputed every
//!   call, exactly as in the original.
//! - `LOGICAL :: X, Y = .FALSE.`-style latched initializers replicated:
//!   `Visited`/`Initialized`/`ComputeTangledMask`/`RecompStabExe` start
//!   false, `Cnt` starts 0, and `MaskPerm` starts disassociated (the
//!   wrapper-side `=> NULL()` init).
//! - `OrigCoord(1:nsize) = Coord(1:nsize)` re-sync, the `Reinitialize`
//!   dealloc shape (only `Field`, then `Surface` — matching the original),
//!   and the `No Matrix` init entry are kept exactly.
//! - `MIN`/`MAX` on reals use explicit `<=`/`>=` comparisons (Fortran
//!   `MIN`/`MAX` return the first argument on ties, unlike `fmin`/`fmax`
//!   signed-zero handling); `TINY(MinHeight)` is `f64::MIN_POSITIVE`.
//! - The unused `status(MPI_STATUS_SIZE)` local of `MultiLayerMapper` is
//!   dropped (it is never read); no `MPI_` call exists in the unit.
//! - Dead locals (`UnitVector`, `ElemVector`, `DotPro`, `Eps`, `Length`,
//!   `Rounds`, `RefVal`, `DOFs`, `at2`, `ii`, `jj`, `Nodes`) are dropped.
//! - `IF(.FALSE.)` dead block (the corrected-height debug print) is
//!   omitted as non-executed.
//! - `dt` is read as in the original (velocity/update paths only);
//!   `TransientSimulation` is unused by the body, exactly as in the
//!   original.
//! - Allocation failure calls the same `Fatal(...)` via the wrapper-side
//!   `STAT=istat` check (the `OrigCoord` allocate), never aborting
//!   directly.
//! - No Rust threads: the solver calls entries serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn sma_begin(h_model: *const c_void, h_solver: *const c_void);
    fn sma_add_nomatrix();
    fn sma_info_banner();
    fn sma_mesh_dim() -> c_int;
    fn sma_reinit_flag(reinit: *mut c_int, found: *mut c_int);
    fn sma_reinit_dealloc();
    fn sma_recomp_setup(stab: *mut c_int, interval: *mut c_int);
    fn sma_fixed_layers(n: *mut c_int, present: *mut c_int);
    fn sma_fixed_get(k: c_int) -> c_int;
    fn sma_botproj_flag(val: *mut c_int, found: *mut c_int);
    fn sma_ptr_dealloc();
    fn sma_detect_multi(nfixed: c_int) -> c_int;
    fn sma_detect_binary() -> c_int;
    fn sma_post_detect(nsize: *mut c_int, nnodes: *mut c_int, mask: *mut c_int);
    fn sma_orig_sync(nsize: c_int);
    fn sma_cputime() -> c_double;
    fn sma_get_string(which: c_int, buf: *mut c_char, nbuf: c_int) -> c_int;
    fn sma_base_resolve(name: *const c_char, namelen: c_int);
    fn sma_base_first(first: *mut c_int, found: *mut c_int);
    fn sma_velo_setup(name: *const c_char, namelen: c_int, dim: c_int, nsize: c_int);
    fn sma_update_setup(name: *const c_char, namelen: c_int, dim: c_int, nsize: c_int);
    fn sma_displacement_mode(val: *mut c_int, found: *mut c_int);
    fn sma_minheight() -> c_double;
    fn sma_info_minheight(mh: c_double);
    fn sma_detangle_info(detangle: *mut c_int, mapheight: *mut c_int);
    fn sma_tangled_setup() -> c_int;
    fn sma_top_setup(mode: *mut c_int, topval0: *mut c_double);
    fn sma_top_present() -> c_int;
    fn sma_fatal_topbc();
    fn sma_top_bc(nsize: c_int, mask: c_int);
    fn sma_bot_setup(mode: *mut c_int, botval0: *mut c_double);
    fn sma_bot_present() -> c_int;
    fn sma_bot_bc(nsize: c_int, mask: c_int);
    fn sma_mid_bc(nsize: c_int);
    fn sma_maskperm_get(i: c_int) -> c_int;
    fn sma_topptr_get(j: c_int) -> c_int;
    fn sma_botptr_get(j: c_int) -> c_int;
    fn sma_midptr_get(j: c_int) -> c_int;
    fn sma_upptr_get(j: c_int) -> c_int;
    fn sma_orig_get(k: c_int) -> c_double;
    fn sma_field_get(k: c_int) -> c_double;
    fn sma_topperm_get(node: c_int) -> c_int;
    fn sma_topfield_at(p: c_int) -> c_double;
    fn sma_botperm_get(node: c_int) -> c_int;
    fn sma_botfield_at(p: c_int) -> c_double;
    fn sma_coord_get(j: c_int) -> c_double;
    fn sma_coord_set(j: c_int, v: c_double);
    fn sma_tangled_mark(i: c_int);
    fn sma_warn_tangle_cancel();
    fn sma_fatal_topperm() -> c_int;
    fn sma_fatal_botperm() -> c_int;
    fn sma_listfun(x: c_double) -> c_double;
    fn sma_velo_update(i: c_int, v: c_double);
    fn sma_update_update(i: c_int, v: c_double);
    fn sma_parreduce(v: c_int) -> c_int;
    fn sma_info_limited(lim: c_int, nsize: c_int);
    fn sma_info_tangled(tan: c_int, nsize: c_int);
    fn sma_passive(nsize: c_int) -> c_int;
    fn sma_actdir() -> c_int;
    fn sma_coordp_sel(dir: c_int);
    fn sma_coordp_set(i: c_int, v: c_double);
    fn sma_velo_firstzero() -> c_int;
    fn sma_velo_zero();
    fn sma_info_maptime(at0: c_double);
    fn sma_maybe_write();
    fn sma_meshstab();
    fn sma_info_multilay_using(nfixed: c_int);
    fn sma_fatal_multimask();
    fn sma_detangle_flag() -> c_int;
    fn sma_fixedvar_setup(nfixed: c_int);
    fn sma_fixedvar_perm(node: c_int) -> c_int;
    fn sma_fixedvar_vals(k: c_int, nf: c_int, buf: *mut c_double);
    fn sma_fatal_nofield();
    fn sma_fatal_nolayer(k: c_int) -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn sma_stride_project(
        proj: *mut c_double,
        fixed: *mut c_double,
        stride: *mut c_double,
        nl: c_int,
        nf: c_int,
        displ: c_int,
    );
    fn sma_stride_commit(ind: *mut c_int, vals: *mut c_double, n: c_int);
    fn sma_stride_velo(ind: *mut c_int, vals: *mut c_double, n: c_int, dt: c_double);
    fn sma_stride_update(ind: *mut c_int, vals: *mut c_double, n: c_int);
    fn sma_info_multilay_done();
    fn sma_info_base_adding();
    fn sma_basedofs() -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn sma_basevals(ibot: c_int, dofs: c_int, bx: *mut c_double, by: *mut c_double, bz: *mut c_double);
    fn sma_node_add(i: c_int, dx: c_double, dy: c_double, dz: c_double, dofs: c_int);
 }

/// First-call SAVE state: the original's plain `SAVE` scalars.
struct SmaState {
    visited: bool,
    initialized: bool,
    top_mode: c_int,
    bot_mode: c_int,
    nsize: c_int,
    nnodes: c_int,
    mask_exists: bool,
    compute_tangled_mask: bool,
    mid_layer_exists: bool,
    number_of_layers: c_int,
    recomp_stab_exe: bool,
    cnt: c_int,
}

impl SmaState {
    const fn new() -> Self {
        Self {
            visited: false,
            initialized: false,
            top_mode: 0,
            bot_mode: 0,
            nsize: 0,
            nnodes: 0,
            mask_exists: false,
            compute_tangled_mask: false,
            mid_layer_exists: false,
            number_of_layers: 0,
            recomp_stab_exe: false,
            cnt: 0,
        }
    }
}

static STATE: Mutex<SmaState> = Mutex::new(SmaState::new());

 fn cstr_len(buf: &[c_char]) -> usize {
     buf.iter().position(|&c| c == 0).unwrap_or(buf.len())
 }
 /// `StructuredMeshMapper_init` — the init entry: flags `No Matrix`.
 ///
 /// # Safety
 ///
 /// `model`/`solver` must be the Elmer addresses the loader passes.
#[export_name = "structuredmeshmapper_init_"]
pub unsafe extern "C" fn structuredmeshmapper_init_entry(
    model: *mut c_void,
    solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe { sma_begin(model, solver) };
    unsafe { sma_add_nomatrix() };
}

 /// `StructuredMeshMapper` — the solver entry. Argument order and
 /// by-reference travel match the Fortran dummy list
 /// `(Model, Solver, dt, TransientSimulation)`.
 ///
 /// # Safety
 ///
 /// `model`/`solver` must be the Elmer addresses the loader passes.
 /// `dt` must be a valid read (it feeds the velocity/update paths).
 #[export_name = "structuredmeshmapper_"]
 pub unsafe extern "C" fn structuredmeshmapper_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    unsafe { sma_begin(model, solver) };
    unsafe { sma_info_banner() };

    let dim = unsafe { sma_mesh_dim() };

    let (mut reinit, mut dummy) = (0, 0);
    unsafe { sma_reinit_flag(&mut reinit, &mut dummy) };
    if reinit != 0 {
        unsafe { sma_reinit_dealloc() };
    }

    let (mut recomp_stab, mut recomp_interval) = (0, 0);
    unsafe { sma_recomp_setup(&mut recomp_stab, &mut recomp_interval) };
    if recomp_stab != 0 {
        state.cnt += 1;
        if state.cnt == recomp_interval {
            state.cnt = 0;
            state.recomp_stab_exe = true;
        } else {
            state.recomp_stab_exe = false;
        }
    } else {
        state.recomp_stab_exe = false;
    }

    let (mut nfixed_layers, mut multilayer_present) = (0, 0);
    unsafe { sma_fixed_layers(&mut nfixed_layers, &mut multilayer_present) };
    let multi_layer = multilayer_present != 0;

    let (mut botproj, mut botproj_found) = (0, 0);
    unsafe { sma_botproj_flag(&mut botproj, &mut botproj_found) };
    let _ = botproj_found;
    let bot_proj = botproj != 0;

    if !state.initialized || reinit != 0 {
        unsafe { sma_ptr_dealloc() };
        if multi_layer {
            state.number_of_layers = unsafe { sma_detect_multi(nfixed_layers) };
        } else {
            state.mid_layer_exists = unsafe { sma_detect_binary() } != 0;
        }
        let (mut nsize, mut nnodes, mut mask) = (0, 0, 0);
        unsafe { sma_post_detect(&mut nsize, &mut nnodes, &mut mask) };
        state.nsize = nsize;
        state.nnodes = nnodes;
        state.mask_exists = mask != 0;
        state.initialized = true;
    }

    let nsize = state.nsize;
    let nnodes = state.nnodes;
    let mask_exists = state.mask_exists;
    let dt_val: c_double = unsafe { *dt };

    unsafe { sma_orig_sync(nsize) };
    let at0 = unsafe { sma_cputime() };

    // Base displacement variable.
    let mut namebuf = [0 as c_char; 2048];
    let mut got_base_var = false;
    let mut base_displace_first = false;
    if unsafe { sma_get_string(1, namebuf.as_mut_ptr(), 2048) } != 0 {
        let len = cstr_len(&namebuf);
        unsafe { sma_base_resolve(namebuf.as_ptr(), len as c_int) };
        got_base_var = true;
        let (mut first, mut found) = (0, 0);
        unsafe { sma_base_first(&mut first, &mut found) };
        let _ = found;
        base_displace_first = first != 0;
    }

    // Mesh velocity variable.
    let mut got_velo_var = false;
    if unsafe { sma_get_string(2, namebuf.as_mut_ptr(), 2048) } != 0 {
        let len = cstr_len(&namebuf);
        unsafe { sma_velo_setup(namebuf.as_ptr(), len as c_int, dim, nsize) };
        got_velo_var = true;
    }

    // Mesh update variable.
    let mut got_update_var = false;
    if unsafe { sma_get_string(3, namebuf.as_mut_ptr(), 2048) } != 0 {
        let len = cstr_len(&namebuf);
        unsafe { sma_update_setup(namebuf.as_ptr(), len as c_int, dim, nsize) };
        got_update_var = true;
    }

    let (mut dm, mut dmf) = (0, 0);
    unsafe { sma_displacement_mode(&mut dm, &mut dmf) };
    let _ = dmf;
    let displacement_mode = dm != 0;

    let min_height = unsafe { sma_minheight() };
    unsafe { sma_info_minheight(min_height) };

    if got_base_var && base_displace_first {
        base_var_displace(mask_exists, nsize);
    }

    let (tangled_count, limited_count) = if multi_layer {
        multi_layer_mapper(
            &state,
            nfixed_layers,
            dt_val,
            got_velo_var,
            got_update_var,
            displacement_mode,
            min_height,
            bot_proj,
        )
    } else {
        binary_layer_mapper(
            &mut state,
            dt_val,
            got_velo_var,
            got_update_var,
            displacement_mode,
            min_height,
        )
    };

    let limited = unsafe { sma_parreduce(limited_count) };
    let tangled = unsafe { sma_parreduce(tangled_count) };

    if limited > 0 {
        unsafe { sma_info_limited(limited, nsize) };
    }
    if tangled > 0 {
        unsafe { sma_info_tangled(tangled, nsize) };
    }

    if unsafe { sma_passive(nsize) } != 0 {
        // Restore done inside the wrapper, exactly as in the original.
    }

    // Masked coordinates are not directly linked to the real coordinate.
    if mask_exists {
        let dir = unsafe { sma_actdir() };
        unsafe { sma_coordp_sel(dir) };
        for i in 1..=nnodes {
            let j = unsafe { sma_maskperm_get(i) };
            if j == 0 {
                continue;
            }
            let v = unsafe { sma_coord_get(j) };
            unsafe { sma_coordp_set(i, v) };
        }
    }

    if got_base_var && !base_displace_first {
        base_var_displace(mask_exists, nsize);
    }

    if got_velo_var && !state.visited
        && unsafe { sma_velo_firstzero() } != 0
    {
        unsafe { sma_velo_zero() };
    }

    unsafe { sma_info_maptime(at0) };

    if !state.visited {
        unsafe { sma_maybe_write() };
    }

    state.visited = true;

    if state.recomp_stab_exe {
        unsafe { sma_meshstab() };
    }
 }

/// The `BinaryLayerMapper` internal procedure: linear interpolation
/// between bottom/top (and optional mid) layers. Returns
/// `(tangled_count, limited_count)`.
#[allow(clippy::too_many_arguments)]
fn binary_layer_mapper(
    state: &mut SmaState,
    dt_val: c_double,
    got_velo_var: bool,
    got_update_var: bool,
    displacement_mode: bool,
    min_height: c_double,
) -> (c_int, c_int) {
    let (mut detangle, mut mapheight) = (0, 0);
    unsafe { sma_detangle_info(&mut detangle, &mut mapheight) };
    let detangle_b = detangle != 0;
    let mapheight_b = mapheight != 0;

    if detangle_b {
        state.compute_tangled_mask = unsafe { sma_tangled_setup() } != 0;
    }
    let compute_mask = state.compute_tangled_mask;

    let (mut top_mode, mut top_val0) = (0, 0.0);
    unsafe { sma_top_setup(&mut top_mode, &mut top_val0) };
    if top_mode == 0 {
        if unsafe { sma_top_present() } != 0 {
            unsafe { sma_top_bc(state.nsize, c_int::from(state.mask_exists)) };
            // A present Top Surface BC entry promotes the mode to 3,
            // exactly as in the original (`TopMode = 3` in the gather).
            top_mode = 3;
        } else {
            unsafe { sma_fatal_topbc() };
            std::process::abort();
        }
    }
    state.top_mode = top_mode;

    let (mut bot_mode, mut bot_val0) = (0, 0.0);
    unsafe { sma_bot_setup(&mut bot_mode, &mut bot_val0) };
    if bot_mode == 0 && unsafe { sma_bot_present() } != 0 {
        unsafe { sma_bot_bc(state.nsize, c_int::from(state.mask_exists)) };
        // A present Bottom Surface BC entry promotes the mode to 3,
        // exactly as in the original (`BotMode = 3` in the gather).
        bot_mode = 3;
    }
    state.bot_mode = bot_mode;

    if state.mid_layer_exists {
        unsafe { sma_mid_bc(state.nsize) };
    }

    let nnodes = state.nnodes;
    let mask_exists = state.mask_exists;
    let mid_exists = state.mid_layer_exists;
    let mut tangled_count: c_int = 0;
    let mut limited_count: c_int = 0;

    for i in 1..=nnodes {
        let mut j = i;
        if mask_exists {
            j = unsafe { sma_maskperm_get(i) };
            if j == 0 {
                continue;
            }
        }
        let itop = unsafe { sma_topptr_get(j) };
        let ibot = unsafe { sma_botptr_get(j) };
        let imid = if mid_exists {
            unsafe { sma_midptr_get(j) }
        } else {
            0
        };

        let (x0top, x0bot, x0loc) = if mask_exists {
            unsafe {
                (
                    sma_orig_get(sma_maskperm_get(itop)),
                    sma_orig_get(sma_maskperm_get(ibot)),
                    sma_orig_get(sma_maskperm_get(i)),
                )
            }
        } else {
            unsafe { (sma_orig_get(itop), sma_orig_get(ibot), sma_orig_get(i)) }
        };

        let mut top_val = if top_mode == 1 {
            top_val0
        } else if top_mode == 2 {
            let p = unsafe { sma_topperm_get(itop) };
            if p == 0 {
                unsafe { sma_fatal_topperm() };
                std::process::abort();
            }
            unsafe { sma_topfield_at(p) }
        } else if top_mode == 3 {
            if mask_exists {
                unsafe { sma_field_get(sma_maskperm_get(itop)) }
            } else {
                unsafe { sma_field_get(itop) }
            }
        } else if displacement_mode {
            0.0
        } else {
            x0top
        };

        let bot_val = if bot_mode == 1 {
            bot_val0
        } else if bot_mode == 2 {
            let p = unsafe { sma_botperm_get(ibot) };
            if p == 0 {
                unsafe { sma_fatal_botperm() };
                std::process::abort();
            }
            unsafe { sma_botfield_at(p) }
        } else if bot_mode == 3 {
            if mask_exists {
                unsafe { sma_field_get(sma_maskperm_get(ibot)) }
            } else {
                unsafe { sma_field_get(ibot) }
            }
        } else if displacement_mode {
            0.0
        } else {
            x0bot
        };
        let (mid_val, x0mid) = if mid_exists {
            unsafe { (sma_field_get(imid), sma_orig_get(imid)) }
        } else {
            (bot_val, x0bot)
        };

        let mut dx = if displacement_mode {
            top_val + x0top - mid_val - x0mid
        } else {
            top_val - mid_val
        };
        let mut tangled = dx < min_height;

        if mask_exists && tangled {
            unsafe { sma_warn_tangle_cancel() };
            tangled = false;
        }

        if tangled {
            if dx < f64::MIN_POSITIVE {
                tangled_count += 1;
            }
            limited_count += 1;

            if detangle_b {
                if mapheight_b {
                    dx = unsafe { sma_listfun(dx / min_height) } * min_height;
                } else {
                    dx = min_height;
                }
                if displacement_mode {
                    top_val = mid_val + x0mid + x0top + dx;
                } else {
                    top_val = mid_val + dx;
                }
            }
            if compute_mask {
                unsafe { sma_tangled_mark(i) };
            }
        }

        let xloc = if mid_exists {
            if (x0top - x0mid) * (x0loc - x0mid) > 0.0 {
                let wtop = (x0loc - x0mid) / (x0top - x0mid);
                wtop * top_val + (1.0 - wtop) * mid_val
            } else {
                let wtop = (x0loc - x0bot) / (x0mid - x0bot);
                wtop * mid_val + (1.0 - wtop) * bot_val
            }
        } else {
            let wtop = (x0loc - x0bot) / (x0top - x0bot);
            wtop * top_val + (1.0 - wtop) * bot_val
        };

        if displacement_mode {
            if got_velo_var {
                unsafe { sma_velo_update(i, xloc / dt_val) };
            }
            let orig_j = unsafe { sma_orig_get(j) };
            let c = orig_j + xloc;
            unsafe { sma_coord_set(j, c) };
            if got_update_var {
                unsafe { sma_update_update(i, c - orig_j) };
            }
        } else {
            if got_velo_var {
                let orig_i = unsafe { sma_orig_get(i) };
                unsafe { sma_velo_update(i, (xloc - orig_i) / dt_val) };
            }
            unsafe { sma_coord_set(j, xloc) };
            if got_update_var {
                let orig_j = unsafe { sma_orig_get(j) };
                unsafe { sma_update_update(i, xloc - orig_j) };
            }
        }
    }

    (tangled_count, limited_count)
}

/// The `MultiLayerMapper` internal procedure: piecewise-linear mapping
/// over `NumberOfLayers` strides. Returns `(tangled, limited)`.
#[allow(clippy::too_many_arguments)]
fn multi_layer_mapper(
    state: &SmaState,
    nfixed_layers: c_int,
    dt_val: c_double,
    got_velo_var: bool,
    got_update_var: bool,
    displacement_mode: bool,
    min_height: c_double,
    bot_proj: bool,
) -> (c_int, c_int) {
    unsafe { sma_info_multilay_using(nfixed_layers) };
    if state.mask_exists {
        unsafe { sma_fatal_multimask() };
        std::process::abort();
    }

    let detangle_b = unsafe { sma_detangle_flag() } != 0;
    unsafe { sma_fixedvar_setup(nfixed_layers) };

    let nlayers = state.number_of_layers as usize;
    let nf = nfixed_layers as usize;
    let nnodes = state.nnodes;
    let displ = c_int::from(displacement_mode);

    let mut proj = vec![0.0 as c_double; nlayers * nf];
    let mut stride_ind = vec![0 as c_int; nlayers];
    let mut stride_coord = vec![0.0 as c_double; nlayers];
    let mut fixed_coord = vec![0.0 as c_double; nf];

    let mut tangled_count: c_int = 0;
    let mut limited_count: c_int = 0;

    for i in 1..=nnodes {
        let ibot = unsafe { sma_botptr_get(i) };
        if ibot != i {
            continue;
        }

        // Build this column's stride from the bottom up.
        let mut jnode = ibot;
        stride_coord[0] = unsafe { sma_orig_get(jnode) };
        stride_ind[0] = jnode;
        for k in 2..=(nlayers as c_int) {
            jnode = unsafe { sma_upptr_get(jnode) };
            stride_coord[(k - 1) as usize] = unsafe { sma_orig_get(jnode) };
            stride_ind[(k - 1) as usize] = jnode;
        }
        let itop = jnode;

        // Projection matrix for this column (column-major in `proj`,
        // matching the Fortran dummy `proj(nl,nf)`).
        for k in 1..=(nlayers as c_int) {
            let mut hit = false;
            let mut jj = 1 as c_int;
            while jj <= nfixed_layers + 1 {
                let fl_j = unsafe { sma_fixed_get(jj) };
                let fl_j1 = unsafe { sma_fixed_get(jj + 1) };
                if fl_j == k {
                    proj[((jj - 1) as usize) * nlayers + (k - 1) as usize] = 1.0;
                    hit = true;
                } else if fl_j < k && fl_j1 > k {
                    let q = (stride_coord[(k - 1) as usize]
                        - stride_coord[(fl_j - 1) as usize])
                        / (stride_coord[(fl_j1 - 1) as usize]
                            - stride_coord[(fl_j - 1) as usize]);
                    proj[(jj as usize) * nlayers + (k - 1) as usize] = q;
                    proj[((jj - 1) as usize) * nlayers + (k - 1) as usize] = 1.0 - q;
                    hit = true;
                }
                if hit {
                    break;
                }
                jj += 1;
            }
            if !hit {
                unsafe { sma_fatal_nolayer(k) };
                std::process::abort();
            }
        }

        // Fixed-layer values live at the top or the bottom, not elsewhere.
        let mut kk = unsafe { sma_fixedvar_perm(ibot) };
        if kk == 0 {
            kk = unsafe { sma_fixedvar_perm(itop) };
        }
        if kk == 0 {
            unsafe { sma_fatal_nofield() };
            std::process::abort();
        }
        unsafe { sma_fixedvar_vals(kk, nfixed_layers, fixed_coord.as_mut_ptr()) };

        unsafe {
            sma_stride_project(
                proj.as_mut_ptr(),
                fixed_coord.as_mut_ptr(),
                stride_coord.as_mut_ptr(),
                nlayers as c_int,
                nfixed_layers,
                displ,
            )
        };

        // Fix the mesh if it becomes tangled (Fortran MIN/MAX select the
        // first argument on ties; explicit comparisons replicate that).
        if detangle_b {
            if bot_proj {
                let mut cnt: c_int = 0;
                for p in 0..nlayers - 1 {
                    if stride_coord[p] - stride_coord[p + 1] < min_height {
                        cnt += 1;
                    }
                }
                if cnt > 0 {
                    limited_count += cnt;
                    let mut t: c_int = 0;
                    for p in 0..nlayers - 1 {
                        if stride_coord[p] - stride_coord[p + 1] < f64::MIN_POSITIVE {
                            t += 1;
                        }
                    }
                    tangled_count += t;
                    for p in 1..nlayers {
                        let a = stride_coord[p];
                        let b = stride_coord[p - 1] - min_height;
                        stride_coord[p] = if a <= b { a } else { b };
                    }
                }
            } else {
                let mut cnt: c_int = 0;
                for p in 0..nlayers - 1 {
                    if stride_coord[p + 1] - stride_coord[p] < min_height {
                        cnt += 1;
                    }
                }
                if cnt > 0 {
                    limited_count += cnt;
                    let mut t: c_int = 0;
                    for p in 0..nlayers - 1 {
                        if stride_coord[p + 1] - stride_coord[p] < f64::MIN_POSITIVE {
                            t += 1;
                        }
                    }
                    tangled_count += t;
                    for p in 1..nlayers {
                        let a = stride_coord[p];
                        let b = stride_coord[p - 1] + min_height;
                        stride_coord[p] = if a >= b { a } else { b };
                    }
                }
            }
        }

        unsafe {
            sma_stride_commit(
                stride_ind.as_mut_ptr(),
                stride_coord.as_mut_ptr(),
                nlayers as c_int,
            )
        };

        if got_velo_var {
            unsafe {
                sma_stride_velo(
                    stride_ind.as_mut_ptr(),
                    stride_coord.as_mut_ptr(),
                    nlayers as c_int,
                    dt_val,
                )
            };
        }
        if got_update_var {
            unsafe {
                sma_stride_update(
                    stride_ind.as_mut_ptr(),
                    stride_coord.as_mut_ptr(),
                    nlayers as c_int,
                )
            };
        }
    }

    unsafe { sma_info_multilay_done() };

    (tangled_count, limited_count)
}

/// The `BaseVarDisplace` internal procedure: add the base displacement.
 fn base_var_displace(mask_exists: bool, nsize: c_int) {
    unsafe { sma_info_base_adding() };
    let dofs = unsafe { sma_basedofs() };
    for i in 1..=nsize {
        if mask_exists {
            let jj = unsafe { sma_maskperm_get(i) };
            if jj == 0 {
                continue;
            }
        }
        let ibot = unsafe { sma_botptr_get(i) };
        let (mut bx, mut by, mut bz) = (0.0, 0.0, 0.0);
        unsafe { sma_basevals(ibot, dofs, &mut bx, &mut by, &mut bz) };
        unsafe { sma_node_add(i, bx, by, bz, dofs) };
    }
 }
