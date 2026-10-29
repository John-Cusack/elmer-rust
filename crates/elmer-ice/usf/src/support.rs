//! Shared BIND(C) accessor declarations and small helpers for the
//! ElmerIceUSF port. Every Elmer-library crossing goes through the
//! generated-style wrappers in `accessors/usf_acc.F90`; this module only
//! declares them and adapts Rust strings/handles at the boundary.

// Re-exported so port modules get the FFI scalar types through the glob.
pub use std::ffi::c_void;
pub use std::os::raw::{c_double, c_int};
pub use std::sync::Mutex;

extern "C" {
    pub fn usf_begin(h_model: *const c_void);
    pub fn usf_current_element(h_elem: *mut *const c_void);
    pub fn usf_get_material(h_elem: *const c_void, h_list: *mut *const c_void);
    pub fn usf_get_bc(h_elem: *const c_void, h_list: *mut *const c_void);
    pub fn usf_get_constants(h_list: *mut *const c_void);
    pub fn usf_get_const_real(
        h_list: *const c_void,
        name: *const u8,
        namelen: c_int,
        found: *mut c_int,
    ) -> c_double;
    pub fn usf_get_logical(
        h_list: *const c_void,
        name: *const u8,
        namelen: c_int,
        val: *mut c_int,
        found: *mut c_int,
    );
    pub fn usf_list_get_logical(
        h_list: *const c_void,
        name: *const u8,
        namelen: c_int,
        val: *mut c_int,
        found: *mut c_int,
    );
    pub fn usf_info(caller: *const u8, ncaller: c_int, msg: *const u8, nmsg: c_int, level: c_int);
    pub fn usf_fatal(caller: *const u8, ncaller: c_int, msg: *const u8, nmsg: c_int);
    pub fn usf_pow_dp(base: c_double, expo: c_double) -> c_double;
    pub fn usf_xy2lonlat(x: c_double, y: c_double, lon: *mut c_double, lat: *mut c_double);
    pub fn usf_lonlat2xy(lon: c_double, lat: c_double, x: *mut c_double, y: *mut c_double);
    pub fn usf_model_variables(h_vars: *mut *const c_void);
    pub fn usf_variable_get(
        h_vars: *const c_void,
        name: *const u8,
        namelen: c_int,
        unfoundfatal: c_int,
        thisonly: c_int,
        h_var: *mut *const c_void,
    );
    pub fn usf_var_value_at(h_var: *const c_void, node: c_int) -> c_double;
    pub fn usf_get_string(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
        found: *mut c_int,
    );
    pub fn usf_get_element_nof_nodes(h_elem: *const c_void) -> c_int;
    pub fn usf_element_nodeindexes(h_elem: *const c_void, buf: *mut c_int, n: c_int);
    pub fn usf_max_element_nodes() -> c_int;
    pub fn usf_get_bodyforce(h_bf: *mut *const c_void);
    pub fn usf_get_bodyforce_id() -> c_int;
    pub fn usf_bodyforce_values(bf_id: c_int, h_list: *mut *const c_void);
    pub fn usf_get_material_id(h_elem: *const c_void, found: *mut c_int) -> c_int;
    pub fn usf_list_get_real(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        n: c_int,
        idx: *const c_int,
        unfoundfatal: c_int,
        found: *mut c_int,
        out: *mut c_double,
    );
    pub fn usf_sum_g2_over_n(g: *const c_double, m: c_int, n: c_int) -> c_double;
    pub fn usf_simulation_values(h_list: *mut *const c_void);
    pub fn usf_number_of_solvers() -> c_int;
    pub fn usf_solver_variable_name(
        solver_i: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_solver_mesh_variables(solver_i: c_int, h_vars: *mut *const c_void);
    pub fn usf_calculate_nodal_weights(solver_i: c_int, varname: *const u8, nvarname: c_int);
    pub fn usf_get_bodyforce_of(h_elem: *const c_void, h_bf: *mut *const c_void);
    pub fn usf_var_perm_at(h_var: *const c_void, node: c_int) -> c_int;
    pub fn usf_var_value_at_idx(h_var: *const c_void, idx: c_int) -> c_double;
    pub fn usf_get_material_current(h_list: *mut *const c_void);
    pub fn usf_coord_dim() -> c_int;
    pub fn usf_even_random() -> c_double;
    pub fn usf_write_msg_af104(
        buf: *mut u8,
        nbuf: c_int,
        prefix: *const u8,
        nprefix: c_int,
        value: c_double,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_al1(
        buf: *mut u8,
        nbuf: c_int,
        prefix: *const u8,
        nprefix: c_int,
        value: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_sortd(n: c_int, a: *mut c_double, b: *mut c_int);
    pub fn usf_dgeev(
        jobvr: c_int,
        sig: *mut c_double,
        wr: *mut c_double,
        wi: *mut c_double,
        vr: *mut c_double,
        ldvr: c_int,
        work: *mut c_double,
        infor: *mut c_int,
    );
    pub fn usf_warn(caller: *const u8, ncaller: c_int, msg: *const u8, nmsg: c_int);
    pub fn usf_var_set_at(h_var: *const c_void, node: c_int, value: c_double);
    pub fn usf_dot(a: *const c_double, b: *const c_double, n: c_int) -> c_double;
    pub fn usf_sum_pow2(a: *const c_double, n: c_int) -> c_double;
    pub fn usf_get_equation_of(
        h_elem: *const c_void,
        h_eq: *mut *const c_void,
        found: *mut c_int,
    );
    pub fn usf_boundary_parent(h_elem: *const c_void, h_parent: *mut *const c_void);
    pub fn usf_write_msg_equation_default(buf: *mut u8, nbuf: c_int, outlen: *mut c_int);
    pub fn usf_nonlin_change() -> c_double;
    pub fn usf_model_mesh_variables(h_vars: *mut *const c_void);
    pub fn usf_mesh_number_of_nodes() -> c_int;
    pub fn usf_number_of_boundary_elements() -> c_int;
    pub fn usf_mype() -> c_int;
    pub fn usf_get_boundary_element(t: c_int, h_elem: *mut *const c_void);
    pub fn usf_element_part_index(h_elem: *const c_void) -> c_int;
    pub fn usf_list_get_string(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
        found: *mut c_int,
    ) -> c_int;
    pub fn usf_list_get_const_real(h_list: *const c_void, key: *const u8, nkey: c_int) -> c_double;
    pub fn usf_set_current_element(h_elem: *const c_void);
    pub fn usf_element_nodes_xyz(
        h_elem: *const c_void,
        n: c_int,
        x: *mut c_double,
        y: *mut c_double,
        z: *mut c_double,
    );
    pub fn usf_write_msg_star_aa(
        pre: *const u8,
        npre: c_int,
        name: *const u8,
        nname: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_star_ar(
        pre: *const u8,
        npre: c_int,
        value: c_double,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_star_arr(
        pre: *const u8,
        npre: c_int,
        x: c_double,
        y: c_double,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_ae148(
        pre: *const u8,
        npre: c_int,
        value: c_double,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_current_solver_values(h_list: *mut *const c_void);
    pub fn usf_get_real(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        n: c_int,
        found: *mut c_int,
        out: *mut c_double,
    );
    pub fn usf_element_body_id(h_elem: *const c_void) -> c_int;
    pub fn usf_model_body_values(body_id: c_int, h_list: *mut *const c_void);
    pub fn usf_model_material_values(material_id: c_int, h_list: *mut *const c_void);
    pub fn usf_list_get_integer(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        unfoundfatal: c_int,
        found: *mut c_int,
        use_mm: c_int,
        minv: c_int,
        maxv: c_int,
    ) -> c_int;
    pub fn usf_get_creal(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        found: *mut c_int,
    ) -> c_double;
    pub fn usf_write_msg_star_s(
        text: *const u8,
        ntext: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_star_asal(
        p1: *const u8,
        np1: c_int,
        name: *const u8,
        nname: c_int,
        p2: *const u8,
        np2: c_int,
        lv: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_star_air(
        pre: *const u8,
        npre: c_int,
        ival: c_int,
        rval: c_double,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_element_type_nodes(h_elem: *const c_void) -> c_int;
    pub fn usf_print_ice_pressure(value: c_double);
    pub fn usf_get_element_family(h_elem: *const c_void) -> c_int;
    pub fn usf_solver_dt() -> c_double;
    pub fn usf_current_solver(h_solver: *mut *const c_void);
    pub fn usf_mesh_tag() -> c_int;
    pub fn usf_mesh_changed() -> c_int;
    pub fn usf_number_of_bodyforces() -> c_int;
    pub fn usf_number_of_materials() -> c_int;
    pub fn usf_list_check_present(h_list: *const c_void, key: *const u8, nkey: c_int) -> c_int;
    pub fn usf_normal_vector(
        h_elem: *const c_void,
        bu: c_double,
        bv: c_double,
        out: *mut c_double,
    ) -> c_int;
    pub fn usf_consistent_normal_vector(
        h_solver: *const c_void,
        h_var: *const c_void,
        h_elem: *const c_void,
        node: c_int,
        out: *mut c_double,
        got: *mut c_int,
    ) -> c_int;
    pub fn usf_element_type_nodeu(h_elem: *const c_void, i: c_int) -> c_double;
    pub fn usf_element_type_nodev(h_elem: *const c_void, i: c_int) -> c_double;
    pub fn usf_element_type_dim(h_elem: *const c_void) -> c_int;
    pub fn usf_model_node_y(node: c_int) -> c_double;
    pub fn usf_model_node_z(node: c_int) -> c_double;
    pub fn usf_write_msg_star_al(
        pre: *const u8,
        npre: c_int,
        v: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_variable_dofs(h_var: *const c_void) -> c_int;
    pub fn usf_element_type_nodew(h_elem: *const c_void, i: c_int) -> c_double;
    pub fn usf_element_info(
        h_elem: *const c_void,
        u: c_double,
        v: c_double,
        w: c_double,
        detj: *mut c_double,
        basis: *mut c_double,
        dbasis: *mut c_double,
        n: c_int,
        stat: *mut c_int,
    );
    pub fn usf_lgrad_dij(
        velo: *const c_double,
        dbasis: *const c_double,
        n: c_int,
        dij: *mut c_double,
    );
    pub fn usf_list_get_real_array(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        n: c_int,
        idx: *const c_int,
        d1: *mut c_int,
        d2: *mut c_int,
        got: *mut c_int,
        out: *mut c_double,
    );
    pub fn usf_write_msg_aea104a(
        pre: *const u8,
        npre: c_int,
        value: c_double,
        suf: *const u8,
        nsuf: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_ai2ai2a(
        p1: *const u8,
        n1: c_int,
        i1: c_int,
        p2: *const u8,
        n2: c_int,
        i2: c_int,
        p3: *const u8,
        n3: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_ai2ai2(
        p1: *const u8,
        n1: c_int,
        i1: c_int,
        p2: *const u8,
        n2: c_int,
        i2: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_write_msg_aaai3a(
        p1: *const u8,
        n1: c_int,
        p2: *const u8,
        n2: c_int,
        p3: *const u8,
        n3: c_int,
        ival: c_int,
        p4: *const u8,
        n4: c_int,
        buf: *mut u8,
        nbuf: c_int,
        outlen: *mut c_int,
    );
    pub fn usf_list_get_const_real_found(
        h_list: *const c_void,
        key: *const u8,
        nkey: c_int,
        found: *mut c_int,
    ) -> c_double;
    pub fn usf_current_solver_mesh_variables(h_vars: *mut *const c_void);
    pub fn usf_model_node_x(node: c_int) -> c_double;
    pub fn usf_cubic_spline(n: c_int, x: *const c_double, y: *const c_double, r: *mut c_double);
    pub fn usf_interpolate_curve(
        n: c_int,
        t: *const c_double,
        f: *const c_double,
        tval: c_double,
        c: *const c_double,
    ) -> c_double;
    pub fn usf_print_newtime();
    pub fn usf_mesh_node_x(node: c_int) -> c_double;
    pub fn usf_mesh_node_y(node: c_int) -> c_double;
    pub fn usf_mesh_node_z(node: c_int) -> c_double;
    pub fn usf_model_mesh(h_mesh: *mut *const c_void);
    pub fn usf_mesh_element(h_mesh: *const c_void, i: c_int, h_elem: *mut *const c_void);
    pub fn usf_mesh_bulk_elements(h_mesh: *const c_void) -> c_int;
    pub fn usf_mesh_max_element_nodes(h_mesh: *const c_void) -> c_int;
    pub fn usf_element_constraint(h_elem: *const c_void) -> c_int;
    pub fn usf_model_bc_values(bc_i: c_int, h_list: *mut *const c_void);
    pub fn usf_model_bc_tag(bc_i: c_int) -> c_int;
    pub fn usf_number_of_bcs() -> c_int;
    pub fn usf_make_perm(which: c_int, name: *const u8, namelen: c_int, nfront: *mut c_int);
    pub fn usf_perm_at(which: c_int, i: c_int) -> c_int;
    pub fn usf_point_in_element(
        h_elem: *const c_void,
        u: c_double,
        v: c_double,
        w: c_double,
        lc: *mut c_double,
        eg: c_double,
        el: c_double,
        en: c_double,
        found: *mut c_int,
    ) -> c_int;
    pub fn usf_interpolate_in_element(
        h_elem: *const c_void,
        vals: *const c_double,
        n: c_int,
        u: c_double,
        v: c_double,
        w: c_double,
    ) -> c_double;
    pub fn usf_mesh_boundary_elements(h_mesh: *const c_void) -> c_int;
}

pub fn cstr(s: &str) -> (*const u8, c_int) {
    (s.as_ptr(), s.len() as c_int)
}

/// `GetConstReal` mirror: returns the value, sets `found`.
pub unsafe fn const_real(h: *const c_void, key: &str, found: &mut bool) -> f64 {
    let (p, n) = cstr(key);
    let mut f: c_int = 0;
    let v = unsafe { usf_get_const_real(h, p, n, &mut f) };
    *found = f != 0;
    v
}

pub unsafe fn get_logical(h: *const c_void, key: &str, found: &mut bool) -> bool {
    let (p, n) = cstr(key);
    let (mut v, mut f): (c_int, c_int) = (0, 0);
    unsafe { usf_get_logical(h, p, n, &mut v, &mut f) };
    *found = f != 0;
    v != 0
}

pub unsafe fn list_get_logical(h: *const c_void, key: &str, found: &mut bool) -> bool {
    let (p, n) = cstr(key);
    let (mut v, mut f): (c_int, c_int) = (0, 0);
    unsafe { usf_list_get_logical(h, p, n, &mut v, &mut f) };
    *found = f != 0;
    v != 0
}

pub unsafe fn info(caller: &str, msg: &str, level: c_int) {
    unsafe {
        usf_info(
            caller.as_ptr(),
            caller.len() as c_int,
            msg.as_ptr(),
            msg.len() as c_int,
            level,
        )
    };
}

pub unsafe fn fatal(caller: &str, msg: &str) -> ! {
    unsafe {
        usf_fatal(
            caller.as_ptr(),
            caller.len() as c_int,
            msg.as_ptr(),
            msg.len() as c_int,
        )
    };
    unreachable!("Fatal never returns in Elmer (STOPs the process)")
}

/// `GetString` mirror: empty string when missing (`found=false`).
pub unsafe fn get_string(h: *const c_void, key: &str, found: &mut bool) -> String {
    let (p, n) = cstr(key);
    let mut buf = vec![0u8; 2048];
    let (mut outlen, mut f): (c_int, c_int) = (0, 0);
    unsafe { usf_get_string(h, p, n, buf.as_mut_ptr(), buf.len() as c_int, &mut outlen, &mut f) };
    *found = f != 0;
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `VariableGet(vars, name, ThisOnly=..., UnFoundFatal=...)` mirror.
pub unsafe fn variable_get(
    h_vars: *const c_void,
    name: &str,
    unfoundfatal: bool,
    thisonly: bool,
) -> *const c_void {
    let (p, n) = cstr(name);
    let mut h: *const c_void = std::ptr::null();
    unsafe {
        usf_variable_get(h_vars, p, n, unfoundfatal as c_int, thisonly as c_int, &mut h)
    };
    h
}

/// `var % values(var % perm(node))` with 1-based `node`.
#[inline]
pub unsafe fn var_value_at(h_var: *const c_void, node: c_int) -> c_double {
    unsafe { usf_var_value_at(h_var, node) }
}

/// `ListGetReal(list, key, n, idx, UnfoundFatal=...)` into `out[0..n]`.
pub unsafe fn list_get_real_into(
    h: *const c_void,
    key: &str,
    idx: &[c_int],
    unfoundfatal: bool,
    found: &mut bool,
    out: &mut [c_double],
) {
    debug_assert_eq!(idx.len(), out.len());
    let (p, nkey) = cstr(key);
    let mut f: c_int = 0;
    unsafe {
        usf_list_get_real(
            h,
            p,
            nkey,
            out.len() as c_int,
            idx.as_ptr(),
            unfoundfatal as c_int,
            &mut f,
            out.as_mut_ptr(),
        )
    };
    *found = f != 0;
}

/// `Model % Simulation` list handle.
pub unsafe fn simulation_values() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_simulation_values(&mut h) };
    h
}

/// `Model % Solvers(i) % Variable % Name` (1-based `i`).
pub unsafe fn solver_variable_name(solver_i: c_int) -> String {
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_solver_variable_name(solver_i, buf.as_mut_ptr(), buf.len() as c_int, &mut outlen)
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `Model % Solvers(i) % Mesh % Variables` handle (1-based `i`).
pub unsafe fn solver_mesh_variables(solver_i: c_int) -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_solver_mesh_variables(solver_i, &mut h) };
    h
}

/// `CalculateNodalWeights(Solvers(i), .FALSE., VarName=name)` (1-based `i`).
pub unsafe fn calculate_nodal_weights(solver_i: c_int, varname: &str) {
    let (p, n) = cstr(varname);
    unsafe { usf_calculate_nodal_weights(solver_i, p, n) };
}

/// `var % perm(node)` (1-based `node`; null/disassociated faults
/// identically to the source — pass-through, never guarded).
#[inline]
pub unsafe fn var_perm_at(h_var: *const c_void, node: c_int) -> c_int {
    unsafe { usf_var_perm_at(h_var, node) }
}

/// `var % values(idx)` with an explicit 1-based `idx`.
#[inline]
pub unsafe fn var_value_at_idx(h_var: *const c_void, idx: c_int) -> c_double {
    unsafe { usf_var_value_at_idx(h_var, idx) }
}

/// `WRITE(buf, '(A,F10.4)') prefix, value` through gfortran, then the
/// caller passes the result to `INFO` — formatting stays bit-identical.
pub unsafe fn write_msg_af104(prefix: &str, value: c_double) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_af104(
            buf.as_mut_ptr(),
            buf.len() as c_int,
            pp,
            pn,
            value,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, '(A,L1)') prefix, value` through gfortran.
pub unsafe fn write_msg_al1(prefix: &str, value: bool) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_al1(
            buf.as_mut_ptr(),
            buf.len() as c_int,
            pp,
            pn,
            value as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `EvenRandom()` — the in-process libgfortran stream, shared with the
/// reference (Elmer seeds it at startup). Call order is transcribed.
#[inline]
pub unsafe fn even_random() -> c_double {
    unsafe { usf_even_random() }
}

/// `WARN(caller, msg)` mirror.
pub unsafe fn warn(caller: &str, msg: &str) {
    unsafe {
        usf_warn(
            caller.as_ptr(),
            caller.len() as c_int,
            msg.as_ptr(),
            msg.len() as c_int,
        )
    };
}

/// `var % values(var % perm(node)) = value` — the Chi write-through.
/// Same live address the reference writes; visibility is by construction.
#[inline]
pub unsafe fn var_set_at(h_var: *const c_void, node: c_int, value: c_double) {
    unsafe { usf_var_set_at(h_var, node, value) };
}

/// `SUM(a(1:n)*b(1:n))` in gfortran (order-sensitive reduction).
#[inline]
pub unsafe fn dot(a: *const c_double, b: *const c_double, n: usize) -> c_double {
    unsafe { usf_dot(a, b, n as c_int) }
}

/// `SUM(a(1:n)**2.0_dp)` in gfortran.
#[inline]
pub unsafe fn sum_pow2(a: *const c_double, n: usize) -> c_double {
    unsafe { usf_sum_pow2(a, n as c_int) }
}

/// `GetEquation(Element=parent, Found=...)` mirror.
pub unsafe fn equation_of(h_elem: *const c_void, found: &mut bool) -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    let mut f: c_int = 0;
    unsafe { usf_get_equation_of(h_elem, &mut h, &mut f) };
    *found = f != 0;
    h
}

/// Boundary-element parent resolution (outbody/Right/Left), verbatim.
pub unsafe fn boundary_parent(h_elem: *const c_void) -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_boundary_parent(h_elem, &mut h) };
    h
}

/// The `(A,I3)` WRITE with no integer item, exactly like the source
/// (missing-item behavior included).
pub unsafe fn write_msg_equation_default() -> String {
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe { usf_write_msg_equation_default(buf.as_mut_ptr(), buf.len() as c_int, &mut outlen) };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `Model % Solver % Variable % NonLinChange`.
#[inline]
pub unsafe fn nonlin_change() -> c_double {
    unsafe { usf_nonlin_change() }
}

/// `Model % Mesh % Variables` handle.
pub unsafe fn model_mesh_variables() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_model_mesh_variables(&mut h) };
    h
}

/// `ListGetString(list, key, Found)` mirror.
pub unsafe fn list_get_string(h: *const c_void, key: &str, found: &mut bool) -> String {
    let (p, n) = cstr(key);
    let mut buf = vec![0u8; 2048];
    let (mut outlen, mut f): (c_int, c_int) = (0, 0);
    unsafe {
        usf_list_get_string(
            h,
            p,
            n,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
            &mut f,
        )
    };
    *found = f != 0;
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `ListGetConstReal(list, key, UnfoundFatal=.TRUE.)` mirror.
#[inline]
pub unsafe fn list_get_const_real(h: *const c_void, key: &str) -> c_double {
    let (p, n) = cstr(key);
    unsafe { usf_list_get_const_real(h, p, n) }
}

/// `ListGetConstReal(list, key, Found)` mirror (no `UnfoundFatal` —
/// missing keys report `found=false` instead of fataling, like the
/// source's bare-`GotIt` call sites).
#[inline]
pub unsafe fn list_get_const_real_found(
    h: *const c_void,
    key: &str,
    found: &mut bool,
) -> c_double {
    let (p, n) = cstr(key);
    let mut f: c_int = 0;
    let v = unsafe { usf_list_get_const_real_found(h, p, n, &mut f) };
    *found = f != 0;
    v
}

/// Restore `Model % CurrentElement` from a saved handle.
#[inline]
pub unsafe fn set_current_element(h_elem: *const c_void) {
    unsafe { usf_set_current_element(h_elem) };
}

/// `WRITE(buf, *) prefix, name` through gfortran (list-directed).
pub unsafe fn write_msg_star_aa(prefix: &str, name_padded: &str) -> String {
    let (pp, pn) = cstr(prefix);
    let (np, nn) = cstr(name_padded);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_aa(
            pp,
            pn,
            np,
            nn,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, *) prefix, value` through gfortran.
pub unsafe fn write_msg_star_ar(prefix: &str, value: c_double) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_ar(
            pp,
            pn,
            value,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, *) prefix, x, y` through gfortran.
pub unsafe fn write_msg_star_arr(prefix: &str, x: c_double, y: c_double) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_arr(
            pp,
            pn,
            x,
            y,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `Model % Solver % Values` handle (null pass-through, like the source).
pub unsafe fn current_solver_values() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_current_solver_values(&mut h) };
    h
}

/// `GetReal(list, key, Found)` into `out[0..n]`.
pub unsafe fn get_real_into(
    h: *const c_void,
    key: &str,
    found: &mut bool,
    out: &mut [c_double],
) {
    let (p, n) = cstr(key);
    let mut f: c_int = 0;
    unsafe {
        usf_get_real(
            h,
            p,
            n,
            out.len() as c_int,
            &mut f,
            out.as_mut_ptr(),
        )
    };
    *found = f != 0;
}

/// `ListGetInteger(list, key, UnfoundFatal=...)` mirror.
pub unsafe fn list_get_integer(
    h: *const c_void,
    key: &str,
    unfoundfatal: bool,
    found: &mut bool,
) -> c_int {
    list_get_integer_mm(h, key, unfoundfatal, found, None)
}

/// `ListGetInteger(list, key, Found, minv=.., maxv=.., UnfoundFatal=...)` mirror.
pub unsafe fn list_get_integer_mm(
    h: *const c_void,
    key: &str,
    unfoundfatal: bool,
    found: &mut bool,
    minmax: Option<(c_int, c_int)>,
) -> c_int {
    let (p, n) = cstr(key);
    let mut f: c_int = 0;
    let (use_mm, minv, maxv) = match minmax {
        Some((lo, hi)) => (1, lo, hi),
        None => (0, 0, 0),
    };
    let v = unsafe {
        usf_list_get_integer(h, p, n, unfoundfatal as c_int, &mut f, use_mm, minv, maxv)
    };
    *found = f != 0;
    v
}

/// `GetCReal(list, key, Found)` mirror.
pub unsafe fn get_creal(h: *const c_void, key: &str, found: &mut bool) -> c_double {
    let (p, n) = cstr(key);
    let mut f: c_int = 0;
    let v = unsafe { usf_get_creal(h, p, n, &mut f) };
    *found = f != 0;
    v
}

/// `WRITE(buf, *) text` through gfortran.
pub unsafe fn write_msg_star_s(text: &str) -> String {
    let (pp, pn) = cstr(text);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_s(pp, pn, buf.as_mut_ptr(), buf.len() as c_int, &mut outlen)
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, *) p1, name, p2, logical` through gfortran.
pub unsafe fn write_msg_star_asal(p1: &str, name: &str, p2: &str, v: bool) -> String {
    let (a, na) = cstr(p1);
    let (b, nb) = cstr(name);
    let (c, nc) = cstr(p2);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_asal(
            a,
            na,
            b,
            nb,
            c,
            nc,
            v as c_int,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, *) prefix, ival, rval` through gfortran.
pub unsafe fn write_msg_star_air(prefix: &str, ival: c_int, rval: c_double) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_air(
            pp,
            pn,
            ival,
            rval,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, '(A,e14.8)') prefix, value` through gfortran.
pub unsafe fn write_msg_ae14d8(prefix: &str, value: c_double) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_ae148(
            pp,
            pn,
            value,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `CurElement % Type % NumberOfNodes`.
#[inline]
pub unsafe fn element_type_nodes(h_elem: *const c_void) -> c_int {
    unsafe { usf_element_type_nodes(h_elem) }
}

/// `print *, "Ice pressure value :", value` — stdout-affecting, so it
/// runs in Fortran with list-directed formatting.
#[inline]
pub unsafe fn print_ice_pressure(value: c_double) {
    unsafe { usf_print_ice_pressure(value) };
}

/// `GetElementFamily(element)`.
#[inline]
pub unsafe fn get_element_family(h_elem: *const c_void) -> c_int {
    unsafe { usf_get_element_family(h_elem) }
}

/// `Model % Solver % dt`.
#[inline]
pub unsafe fn solver_dt() -> c_double {
    unsafe { usf_solver_dt() }
}

/// `Model % Solver` handle (null pass-through).
pub unsafe fn current_solver() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_current_solver(&mut h) };
    h
}

/// `ListCheckPresent(list, key)` mirror.
#[inline]
pub unsafe fn list_check_present(h: *const c_void, key: &str) -> bool {
    let (p, n) = cstr(key);
    unsafe { usf_list_check_present(h, p, n) != 0 }
}

/// `NormalVector(elem, nodes, bu, bv, .TRUE.)` into `out[0..3]`.
pub unsafe fn normal_vector(h_elem: *const c_void, bu: c_double, bv: c_double, out: &mut [c_double; 3]) {
    unsafe { usf_normal_vector(h_elem, bu, bv, out.as_mut_ptr()) };
}

/// `ConsistentNormalVector(solver, var, elem, Found, node=node)` into `out`.
pub unsafe fn consistent_normal_vector(
    h_solver: *const c_void,
    h_var: *const c_void,
    h_elem: *const c_void,
    node: c_int,
    out: &mut [c_double; 3],
    got: &mut bool,
) {
    let mut f: c_int = 0;
    unsafe {
        usf_consistent_normal_vector(h_solver, h_var, h_elem, node, out.as_mut_ptr(), &mut f)
    };
    *got = f != 0;
}

/// `WRITE(buf, *) prefix, logical` through gfortran.
pub unsafe fn write_msg_star_al(prefix: &str, v: bool) -> String {
    let (pp, pn) = cstr(prefix);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_star_al(pp, pn, v as c_int, buf.as_mut_ptr(), buf.len() as c_int, &mut outlen)
    };
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `Model % Mesh % MeshTag`.
#[inline]
pub unsafe fn mesh_tag() -> c_int {
    unsafe { usf_mesh_tag() }
}

/// `Model % Mesh % Changed`.
#[inline]
pub unsafe fn mesh_changed() -> bool {
    unsafe { usf_mesh_changed() != 0 }
}

/// `Model % NumberOfBodyForces`.
#[inline]
pub unsafe fn number_of_bodyforces() -> c_int {
    unsafe { usf_number_of_bodyforces() }
}

/// `Model % NumberOfMaterials`.
#[inline]
pub unsafe fn number_of_materials() -> c_int {
    unsafe { usf_number_of_materials() }
}

/// `BCElement % Type % NodeU(i)` (1-based `i`).
#[inline]
pub unsafe fn element_type_nodeu(h_elem: *const c_void, i: c_int) -> c_double {
    unsafe { usf_element_type_nodeu(h_elem, i) }
}

/// `BCElement % Type % NodeV(i)` (1-based `i`).
#[inline]
pub unsafe fn element_type_nodev(h_elem: *const c_void, i: c_int) -> c_double {
    unsafe { usf_element_type_nodev(h_elem, i) }
}

/// `BCElement % Type % NodeW(i)` (1-based `i`).
#[inline]
pub unsafe fn element_type_nodew(h_elem: *const c_void, i: c_int) -> c_double {
    unsafe { usf_element_type_nodew(h_elem, i) }
}

/// `BCElement % Type % Dimension`.
#[inline]
pub unsafe fn element_type_dim(h_elem: *const c_void) -> c_int {
    unsafe { usf_element_type_dim(h_elem) }
}

/// `Model % Nodes % y(node)` (1-based global node number).
#[inline]
pub unsafe fn model_node_y(node: c_int) -> c_double {
    unsafe { usf_model_node_y(node) }
}

/// `Model % Nodes % z(node)` (1-based global node number).
#[inline]
pub unsafe fn model_node_z(node: c_int) -> c_double {
    unsafe { usf_model_node_z(node) }
}

/// `Variable % DOFs`.
#[inline]
pub unsafe fn variable_dofs(h_var: *const c_void) -> c_int {
    unsafe { usf_variable_dofs(h_var) }
}

/// `ElementInfo(elem, nodes, u, v, w, ...)` into `basis[0..n]`,
/// `dbasis[0..3*n]` (column-major `n x 3`); returns `(detj, stat)`.
pub unsafe fn element_info(
    h_elem: *const c_void,
    u: c_double,
    v: c_double,
    w: c_double,
    n: usize,
    basis: &mut [c_double],
    dbasis: &mut [c_double],
) -> (c_double, bool) {
    let mut detj = 0.0;
    let mut stat: c_int = 0;
    unsafe {
        usf_element_info(
            h_elem,
            u,
            v,
            w,
            &mut detj,
            basis.as_mut_ptr(),
            dbasis.as_mut_ptr(),
            n as c_int,
            &mut stat,
        )
    };
    (detj, stat != 0)
}

/// `LGrad = MATMUL(Velo, dBasisdx); Dij = 0.5*(LGrad+TRANSPOSE)` in
/// gfortran (order-sensitive); `velo` is column-major `3 x n`,
/// `dbasis` column-major `n x 3`, `dij` column-major `3 x 3`.
pub unsafe fn lgrad_dij(velo: &[c_double], dbasis: &[c_double], n: usize, dij: &mut [c_double; 9]) {
    debug_assert_eq!(velo.len(), 3 * n);
    debug_assert_eq!(dbasis.len(), 3 * n);
    unsafe { usf_lgrad_dij(velo.as_ptr(), dbasis.as_ptr(), n as c_int, dij.as_mut_ptr()) };
}

/// `ListGetRealArray(list, key, ni, indexes, Found)` into `out[0..3*n]`
/// (up to 3 rows); returns `(d1, d2, got)`.
pub unsafe fn list_get_real_array_into(
    h: *const c_void,
    key: &str,
    idx: &[c_int],
    got: &mut bool,
    out: &mut [c_double],
) -> (c_int, c_int) {
    let (p, nkey) = cstr(key);
    let (mut d1, mut d2, mut g): (c_int, c_int, c_int) = (0, 0, 0);
    let n = idx.len() as c_int;
    debug_assert!(out.len() >= 3 * idx.len());
    unsafe {
        usf_list_get_real_array(
            h,
            p,
            nkey,
            n,
            idx.as_ptr(),
            &mut d1,
            &mut d2,
            &mut g,
            out.as_mut_ptr(),
        )
    };
    *got = g != 0;
    (d1, d2)
}

fn write_buf_to_string(buf: &[u8], outlen: c_int) -> String {
    let m = (outlen as usize).min(buf.len());
    String::from_utf8_lossy(&buf[..m]).into_owned()
}

/// `WRITE(buf, '(a,e10.4,a)') pre, value, suf` through gfortran.
pub unsafe fn write_msg_aea104a(pre: &str, value: c_double, suf: &str) -> String {
    let (pp, pn) = cstr(pre);
    let (sp, sn) = cstr(suf);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_aea104a(
            pp,
            pn,
            value,
            sp,
            sn,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    write_buf_to_string(&buf, outlen)
}

/// `WRITE(buf, '(a,I2,a,I2,a)')` through gfortran.
pub unsafe fn write_msg_ai2ai2a(p1: &str, i1: c_int, p2: &str, i2: c_int, p3: &str) -> String {
    let (a, na) = cstr(p1);
    let (b, nb) = cstr(p2);
    let (c, nc) = cstr(p3);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_ai2ai2a(
            a,
            na,
            i1,
            b,
            nb,
            i2,
            c,
            nc,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    write_buf_to_string(&buf, outlen)
}

/// `WRITE(buf, '(a,I2,a,I2)')` through gfortran.
pub unsafe fn write_msg_ai2ai2(p1: &str, i1: c_int, p2: &str, i2: c_int) -> String {
    let (a, na) = cstr(p1);
    let (b, nb) = cstr(p2);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_ai2ai2(
            a,
            na,
            i1,
            b,
            nb,
            i2,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    write_buf_to_string(&buf, outlen)
}

/// `WRITE(buf, '(A,A,A,i3,A)')` through gfortran.
pub unsafe fn write_msg_aaai3a(p1: &str, p2: &str, p3: &str, ival: c_int, p4: &str) -> String {
    let (a, na) = cstr(p1);
    let (b, nb) = cstr(p2);
    let (c, nc) = cstr(p3);
    let (d, nd) = cstr(p4);
    let mut buf = vec![0u8; 2048];
    let mut outlen: c_int = 0;
    unsafe {
        usf_write_msg_aaai3a(
            a,
            na,
            b,
            nb,
            c,
            nc,
            ival,
            d,
            nd,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
        )
    };
    write_buf_to_string(&buf, outlen)
}

/// `Model % Solver % Mesh % Variables` handle.
pub unsafe fn current_solver_mesh_variables() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_current_solver_mesh_variables(&mut h) };
    h
}

/// `Model % Nodes % x(node)` (1-based global node number).
#[inline]
pub unsafe fn model_node_x(node: c_int) -> c_double {
    unsafe { usf_model_node_x(node) }
}

/// `CubicSpline(n, x, y, r)` in gfortran (order-sensitive tridiagonal).
pub unsafe fn cubic_spline(n: usize, x: &[c_double], y: &[c_double], r: &mut [c_double]) {
    debug_assert_eq!(x.len(), n);
    debug_assert_eq!(y.len(), n);
    debug_assert!(r.len() >= n);
    unsafe { usf_cubic_spline(n as c_int, x.as_ptr(), y.as_ptr(), r.as_mut_ptr()) };
}

/// `InterpolateCurve(t, f, tval, c)` in gfortran.
pub unsafe fn interpolate_curve(
    t: &[c_double],
    f: &[c_double],
    tval: c_double,
    c: &[c_double],
) -> c_double {
    debug_assert_eq!(t.len(), f.len());
    debug_assert_eq!(t.len(), c.len());
    unsafe {
        usf_interpolate_curve(
            t.len() as c_int,
            t.as_ptr(),
            f.as_ptr(),
            tval,
            c.as_ptr(),
        )
    }
}

/// `PRINT *, "USF_Zs: New time, saving previous Zs values!"`.
#[inline]
pub unsafe fn print_newtime() {
    unsafe { usf_print_newtime() };
}

/// `Model % Mesh % Nodes % x(node)` — the mesh handle, not the model
/// alias (the source distinguishes them; see MeshBasics `Nodes =>`).
#[inline]
pub unsafe fn mesh_node_x(node: c_int) -> c_double {
    unsafe { usf_mesh_node_x(node) }
}

/// `Model % Mesh % Nodes % y(node)`.
#[inline]
pub unsafe fn mesh_node_y(node: c_int) -> c_double {
    unsafe { usf_mesh_node_y(node) }
}

/// `Model % Mesh % Nodes % z(node)`.
#[inline]
pub unsafe fn mesh_node_z(node: c_int) -> c_double {
    unsafe { usf_mesh_node_z(node) }
}

/// `Model % Mesh` handle.
pub unsafe fn model_mesh() -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_model_mesh(&mut h) };
    h
}

/// `Mesh % Elements(i)` handle (1-based `i`).
pub unsafe fn mesh_element(h_mesh: *const c_void, i: c_int) -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_mesh_element(h_mesh, i, &mut h) };
    h
}

/// `PointInElement(elem, nodes, [u,v,w], ...)`; returns `(found, lc)`.
pub unsafe fn point_in_element(
    h_elem: *const c_void,
    u: c_double,
    v: c_double,
    w: c_double,
    eg: c_double,
    el: c_double,
    en: c_double,
) -> (bool, [c_double; 3]) {
    let mut lc = [0.0 as c_double; 3];
    let mut f: c_int = 0;
    unsafe { usf_point_in_element(h_elem, u, v, w, lc.as_mut_ptr(), eg, el, en, &mut f) };
    (f != 0, lc)
}

/// `InterpolateInElement(elem, vals[0..n], u, v, w)`.
pub unsafe fn interpolate_in_element(
    h_elem: *const c_void,
    vals: &[c_double],
    u: c_double,
    v: c_double,
    w: c_double,
) -> c_double {
    unsafe { usf_interpolate_in_element(h_elem, vals.as_ptr(), vals.len() as c_int, u, v, w) }
}

/// `Mesh % NumberOfBoundaryElements`.
#[inline]
pub unsafe fn mesh_boundary_elements(h_mesh: *const c_void) -> c_int {
    unsafe { usf_mesh_boundary_elements(h_mesh) }
}

/// `Mesh % MaxElementNodes`.
#[inline]
pub unsafe fn mesh_max_element_nodes(h_mesh: *const c_void) -> c_int {
    unsafe { usf_mesh_max_element_nodes(h_mesh) }
}

/// `MakePermUsingMask(model, solver, mesh, name, .FALSE., perm, n, ParallelComm=.FALSE.)`
/// for mask slot `which` (1 = front, 2 = bottom, 3 = surface); returns `n`.
pub unsafe fn make_perm(which: c_int, name: &str) -> c_int {
    let (p, n) = cstr(name);
    let mut nfront: c_int = 0;
    unsafe { usf_make_perm(which, p, n, &mut nfront) };
    nfront
}

/// `Model % NumberOfBCs`.
#[inline]
pub unsafe fn number_of_bcs() -> c_int {
    unsafe { usf_number_of_bcs() }
}

/// `Model % BCs(i) % Values` handle (1-based `i`).
pub unsafe fn model_bc_values(bc_i: c_int) -> *const c_void {
    let mut h: *const c_void = std::ptr::null();
    unsafe { usf_model_bc_values(bc_i, &mut h) };
    h
}

/// `Model % BCs(i) % Tag` (1-based `i`).
#[inline]
pub unsafe fn model_bc_tag(bc_i: c_int) -> c_int {
    unsafe { usf_model_bc_tag(bc_i) }
}

/// `Mesh % NumberOfBulkElements`.
#[inline]
pub unsafe fn mesh_bulk_elements(h_mesh: *const c_void) -> c_int {
    unsafe { usf_mesh_bulk_elements(h_mesh) }
}

/// Mask-perm read: `FrontPerm`/`BottomPerm`/`SurfacePerm(i)` for slot
/// `which` (1/2/3), 1-based `i`. Pass-through, never guarded.
#[inline]
pub unsafe fn perm_at(which: c_int, i: c_int) -> c_int {
    unsafe { usf_perm_at(which, i) }
}
