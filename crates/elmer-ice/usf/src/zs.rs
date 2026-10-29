//! `zsini_`, `zsmzsini_`, `zstopini_`, `zstopmzsini_`, `zsbottomini_`,
//! `zsbottommzsini_`, `dyini_`, `dymdyini_`, `zsbottommzsinicalving_`,
//! `zstopmzsinicalving_` — port of USF_Zs.F90 (initial free-surface
//! snapshots and calving-front masks).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - The eight plain entries snapshot nodal coordinates on `FirstTime`
//!   (`Zs0` sized by `COUNT(Perm > 0)` over the named variable, or by
//!   `NumberOfNodes` for the `Dy` pair) and return the snapshot or its
//!   difference — verbatim, including 1-based mappings shifted by one
//!   in the 0-based stores.
//! - `UnFoundFatal` without an initializer (the non-`ZsIni` plain
//!   entries) is BSS-false in the reference; the port passes `false`.
//! - The calving pair re-fetches the variable handle per call; the
//!   reference SAVEs the `ZsPerm` pointer (dangling across a remesh —
//!   undefined behavior the port does not replicate). No graded run
//!   remeshes between the snapshot and the read; the staged calving
//!   differentials prove the defined paths.
//! - `Newtime` (lowercase `n` in the source — same variable, Fortran is
//!   case-insensitive) starts `.TRUE.` (implicit SAVE init).
//! - The `Debug`-gated `INFO`/`WARN`/`PRINT` arms are dead (`Debug` is
//!   set `.FALSE.` on every call); omitted, documented. The live
//!   `PRINT *, 'USF_Zs: ...'` runs in a Fortran accessor (stdout +
//!   list-directed formatting).

use crate::support::*;

struct SnapState {
    first: bool,
    dim: c_int,
    nmax: c_int,
    zs0: Vec<c_double>,
}

fn snap_state() -> Mutex<SnapState> {
    Mutex::new(SnapState {
        first: true,
        dim: 0,
        nmax: 0,
        zs0: Vec::new(),
    })
}

/// Shared snapshot fill for the `(variable, perm)` entries.
/// Returns the mapping count; `zs0` holds `Zs0(perm)` 0-based.
unsafe fn fill_snapshot(
    h_var: *const c_void,
    nnodes: usize,
    dim: c_int,
    zs0: &mut Vec<c_double>,
) -> usize {
    // COUNT(Perm > 0) then fill — two passes like the source's
    // COUNT + loop (values identical; order identical).
    let mut nmax = 0;
    for i in 1..=nnodes {
        if unsafe { var_perm_at(h_var, i as c_int) } > 0 {
            nmax += 1;
        }
    }
    zs0.clear();
    zs0.resize(nmax, 0.0);
    for i in 1..=nnodes {
        let p = unsafe { var_perm_at(h_var, i as c_int) };
        if p == 0 {
            continue;
        }
        let v = if dim == 2 {
            unsafe { model_node_y(i as c_int) }
        } else {
            unsafe { model_node_z(i as c_int) }
        };
        zs0[p as usize - 1] = v;
    }
    nmax
}

struct ZsDb {
    ini: Mutex<SnapState>,
    mzs: Mutex<SnapState>,
    top_ini: Mutex<SnapState>,
    top_mzs: Mutex<SnapState>,
    bot_ini: Mutex<SnapState>,
    bot_mzs: Mutex<SnapState>,
    dy_ini: Mutex<SnapState>,
    dy_mzs: Mutex<SnapState>,
}

static ZS_DB: std::sync::LazyLock<ZsDb> = std::sync::LazyLock::new(|| ZsDb {
    ini: snap_state(),
    mzs: snap_state(),
    top_ini: snap_state(),
    top_mzs: snap_state(),
    bot_ini: snap_state(),
    bot_mzs: snap_state(),
    dy_ini: snap_state(),
    dy_mzs: snap_state(),
});

/// Snapshot a named variable's initial field (the 6 plain `Zs*` entries).
unsafe fn snapshot_var(
    m: *mut c_void,
    varname: &str,
    unfoundfatal: bool,
    st: &Mutex<SnapState>,
) -> Vec<c_double> {
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_var = unsafe { variable_get(h_vars, varname, unfoundfatal, false) };
    let mut s = st.lock().unwrap();
    if s.first {
        s.first = false;
        s.dim = unsafe { usf_coord_dim() };
        let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
        let dim = s.dim;
        let nmax = unsafe { fill_snapshot(h_var, nnodes, dim, &mut s.zs0) };
        s.nmax = nmax as c_int;
    }
    s.zs0.clone()
}

unsafe fn perm_of(varname: &str, node: c_int) -> c_int {
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_var = unsafe { variable_get(h_vars, varname, true, false) };
    unsafe { var_perm_at(h_var, node) }
}

pub unsafe fn zs_ini(m: *mut c_void, n: *const c_int, _v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let zs0 = unsafe { snapshot_var(m, "Zs", true, &ZS_DB.ini) };
    let p = unsafe { perm_of("Zs", node) };
    zs0[p as usize - 1]
}

pub unsafe fn zs_mzs_ini(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    let zs0 = unsafe { snapshot_var(m, "Zs", false, &ZS_DB.mzs) };
    let p = unsafe { perm_of("Zs", node) };
    z - zs0[p as usize - 1]
}

pub unsafe fn zs_top_ini(m: *mut c_void, n: *const c_int, _v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let zs0 = unsafe { snapshot_var(m, "Zs Top", false, &ZS_DB.top_ini) };
    let p = unsafe { perm_of("Zs Top", node) };
    zs0[p as usize - 1]
}

pub unsafe fn zs_top_mzs_ini(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    let zs0 = unsafe { snapshot_var(m, "Zs Top", false, &ZS_DB.top_mzs) };
    let p = unsafe { perm_of("Zs Top", node) };
    z - zs0[p as usize - 1]
}

pub unsafe fn zs_bottom_ini(m: *mut c_void, n: *const c_int, _v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let zs0 = unsafe { snapshot_var(m, "Zs Bottom", false, &ZS_DB.bot_ini) };
    let p = unsafe { perm_of("Zs Bottom", node) };
    zs0[p as usize - 1]
}

pub unsafe fn zs_bottom_mzs_ini(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    let zs0 = unsafe { snapshot_var(m, "Zs Bottom", false, &ZS_DB.bot_mzs) };
    let p = unsafe { perm_of("Zs Bottom", node) };
    z - zs0[p as usize - 1]
}

pub unsafe fn dy_ini(m: *mut c_void, n: *const c_int, _v: *const c_double) -> c_double {
    let node = unsafe { *n };
    unsafe { usf_begin(m) };
    let mut s = ZS_DB.dy_ini.lock().unwrap();
    if s.first {
        s.first = false;
        let nmax = unsafe { usf_mesh_number_of_nodes() } as usize;
        s.nmax = nmax as c_int;
        s.dim = unsafe { usf_coord_dim() };
        let dim = s.dim;
        s.zs0 = vec![0.0; nmax];
        for i in 1..=nmax {
            s.zs0[i - 1] = if dim == 2 {
                unsafe { model_node_y(i as c_int) }
            } else {
                unsafe { model_node_z(i as c_int) }
            };
        }
    }
    let out = s.zs0[node as usize - 1];
    drop(s);
    out
}

pub unsafe fn dy_mdy_ini(m: *mut c_void, n: *const c_int, v: *const c_double) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut s = ZS_DB.dy_mzs.lock().unwrap();
    if s.first {
        s.first = false;
        let nmax = unsafe { usf_mesh_number_of_nodes() } as usize;
        s.nmax = nmax as c_int;
        s.dim = unsafe { usf_coord_dim() };
        let dim = s.dim;
        s.zs0 = vec![0.0; nmax];
        for i in 1..=nmax {
            s.zs0[i - 1] = if dim == 2 {
                unsafe { model_node_y(i as c_int) }
            } else {
                unsafe { model_node_z(i as c_int) }
            };
        }
    }
    let out = z - s.zs0[node as usize - 1];
    drop(s);
    out
}

// -- Calving pair -------------------------------------------------------

struct CalvState {
    first: bool,
    new_time: bool,
    told: c_double,
    dim: c_int,
    nmax: c_int,
    zs0: Vec<c_double>,
    zs_prev: Vec<c_double>,
}

fn calv_state() -> Mutex<CalvState> {
    Mutex::new(CalvState {
        first: true,
        new_time: true,
        told: 0.0,
        dim: 0,
        nmax: 0,
        zs0: Vec::new(),
        zs_prev: Vec::new(),
    })
}

static CALV_BOT: std::sync::LazyLock<Mutex<CalvState>> =
    std::sync::LazyLock::new(calv_state);
static CALV_TOP: std::sync::LazyLock<Mutex<CalvState>> =
    std::sync::LazyLock::new(calv_state);

/// Shared calving-mask body (`Zs Bottom` / `Zs Top` variants).
unsafe fn calving_body(
    m: *mut c_void,
    node: c_int,
    z: c_double,
    varname: &str,
    fatal_caller: &str,
    fatal_missing: &str,
    st: &Mutex<CalvState>,
) -> c_double {
    unsafe { usf_begin(m) };
    let mut h_vars: *const c_void = std::ptr::null();
    unsafe { usf_model_variables(&mut h_vars) };
    let h_time = unsafe { variable_get(h_vars, "Time", true, false) };
    let t = unsafe { var_value_at_idx(h_time, 1) };
    let h_var = unsafe { variable_get(h_vars, varname, false, false) };
    if h_var.is_null() {
        unsafe { fatal(fatal_caller, fatal_missing) };
    }
    {
        let mut s = st.lock().unwrap();
        if s.first {
            s.first = false;
            s.new_time = true;
            s.told = t;
            s.dim = unsafe { usf_coord_dim() };
            let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
            let dim = s.dim;
            let nmax = unsafe {
                let mut c = 0;
                for i in 1..=nnodes {
                    if var_perm_at(h_var, i as c_int) > 0 {
                        c += 1;
                    }
                }
                c
            };
            s.nmax = nmax as c_int;
            s.zs0 = vec![0.0; nmax];
            s.zs_prev = vec![0.0; nmax];
            for i in 1..=nnodes {
                let p = unsafe { var_perm_at(h_var, i as c_int) };
                if p == 0 {
                    continue;
                }
                let v = if dim == 2 {
                    unsafe { model_node_y(i as c_int) }
                } else {
                    unsafe { model_node_z(i as c_int) }
                };
                s.zs0[p as usize - 1] = v;
            }
        }
        if t > s.told {
            s.new_time = true;
        }
        let h_sim = unsafe { simulation_values() };
        let mut found = false;
        let calving = unsafe { list_get_logical(h_sim, "CalvingOccurs", &mut found) };
        let calving = if found { calving } else { false };
        let mut found2 = false;
        let remesh = unsafe { list_get_logical(h_sim, "RemeshOccurs", &mut found2) };
        let remesh = if found2 { remesh } else { false };
        if s.new_time {
            unsafe { print_newtime() };
            s.told = t;
            s.new_time = false;
            let nnodes = unsafe { usf_mesh_number_of_nodes() } as usize;
            let dim = s.dim;
            for i in 1..=nnodes {
                let p = unsafe { var_perm_at(h_var, i as c_int) };
                if p == 0 {
                    continue;
                }
                let v = if dim == 2 {
                    unsafe { model_node_y(i as c_int) }
                } else {
                    unsafe { model_node_z(i as c_int) }
                };
                s.zs_prev[p as usize - 1] = v;
            }
        }
        if calving || remesh {
            let n = s.zs_prev.clone();
            s.zs0.copy_from_slice(&n);
        }
        let p = unsafe { var_perm_at(h_var, node) };
        z - s.zs0[p as usize - 1]
    }
}

pub unsafe fn zs_bottom_mzs_ini_calving(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    unsafe {
        calving_body(
            m,
            node,
            z,
            "Zs Bottom",
            "ZsTopMZsIni",
            "Could not find variable >Zs Bottom<",
            &CALV_BOT,
        )
    }
}

pub unsafe fn zs_top_mzs_ini_calving(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    let node = unsafe { *n };
    let z = unsafe { *v };
    unsafe {
        calving_body(
            m,
            node,
            z,
            "Zs Top",
            "ZsTopMZsIni",
            "Could not find variable >Zs Top<",
            &CALV_TOP,
        )
    }
}
