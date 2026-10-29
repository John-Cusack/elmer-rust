//! Whole-library Rust port of `elmerice/UserFunctions/` (Elmer @ 9f6af2f85):
//! one Rust cdylib replacing `ElmerIceUSF.so`, with the same dynamic
//! surface (`nm -D` equality over the 60 reference T symbols).
//!
//! Boundary (user-function entries): each entry's `Model` dummy crosses as
//! an opaque address, `nodenumber`/`Node` by reference, the value slot as a
//! base address (`ElmerRealFn`'s `val(*)` assumed-size — no descriptor
//! crosses, verified in `fem/src/LoadMod.F90`), and the single 4-arg vector
//! entry (`glaciermeshmetricaniso_`) through the Fortran shim in
//! `accessors/usf_acc.F90` (assumed-shape `TargetLength(:)` becomes
//! base+length there; the descriptor never crosses into Rust).
//! Every Elmer library call goes through the generated-style Fortran
//! `BIND(C)` wrappers in `accessors/usf_acc.F90`, compiled against the
//! built `.mod` files of the target ice-enabled configure — the port
//! itself contains no layout assumption and no pinned offset.
//! Order-sensitive operations (`**`, reductions, LAPACK) evaluate in
//! Fortran wrappers, so bit-identity holds by construction wherever the
//! call sequence matches.
//!
//! Shared rules:
//! - `SAVE` permanent storage is replicated as Rust `static Mutex<State>`
//!   with identical first-call reads; allocation failure calls `Fatal`
//!   via the accessor, never aborting differently.
//! - `Model`, `nodenumber`/`Node` are received and ignored exactly where
//!   the original ignores them.
//! - No Rust threads: entries run serially; each `Mutex` is only a
//!   static-sharing guard.
//! - INFO/WARN/FATAL text, levels, and call sites match the source
//!   exactly (byte-identical solver logs depend on it).
// The cdylib file must be `libElmerIceUSF.so` (it shadows the reference by
// filename), so the library keeps the upstream-mixed-case name.
#![allow(non_snake_case)]
// Lane convention (mirrors elmer-ice mesher crates): the single Elmer-entry
// safety contract lives in the Boundary section above; `// SAFETY:`
// comments mark each deref site instead of per-fn doc sections.
#![allow(clippy::missing_safety_doc)]
mod buoyancy;
mod caffe;
mod metric;
mod shape;
mod zs;
mod contact;
mod coupling;
mod cov;
mod damage;
mod friction;
mod lateral;
mod sliding;
mod haf;
mod iceproperties;
mod proj;
mod srccalc;
mod support;
mod watertransfer;

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

#[export_name = "tenpowera_"]
pub unsafe extern "C" fn tenpowera_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::tenpowera(m, n, v) }
}

#[export_name = "tenpowera_d_"]
pub unsafe extern "C" fn tenpowera_d_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::tenpowera_d(m, n, v) }
}

#[export_name = "derivative_tenpowera_"]
pub unsafe extern "C" fn derivative_tenpowera_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::derivative_tenpowera(m, n, v) }
}

#[export_name = "log10a_"]
pub unsafe extern "C" fn log10a_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::log10a(m, n, v) }
}

#[export_name = "asquare_"]
pub unsafe extern "C" fn asquare_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::asquare(m, n, v) }
}

#[export_name = "asquare_d_"]
pub unsafe extern "C" fn asquare_d_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::asquare_d(m, n, v) }
}

#[export_name = "derivative_asquare_"]
pub unsafe extern "C" fn derivative_asquare_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::derivative_asquare(m, n, v) }
}

#[export_name = "sqrta_"]
pub unsafe extern "C" fn sqrta_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::sqrta(m, n, v) }
}

#[export_name = "tenpowera_masked_"]
pub unsafe extern "C" fn tenpowera_masked_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::tenpowera_masked(m, n, v) }
}

#[export_name = "tenpowera_d_masked_"]
pub unsafe extern "C" fn tenpowera_d_masked_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::tenpowera_d_masked(m, n, v) }
}

#[export_name = "asquare_scaled_"]
pub unsafe extern "C" fn asquare_scaled_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::asquare_scaled(m, n, v) }
}

#[export_name = "asquare_d_scaled_"]
pub unsafe extern "C" fn asquare_d_scaled_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { cov::asquare_d_scaled(m, n, v) }
}

#[export_name = "xy2lon_"]
pub unsafe extern "C" fn xy2lon_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { proj::xy2lon(m, n, v) }
}

#[export_name = "xy2lat_"]
pub unsafe extern "C" fn xy2lat_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { proj::xy2lat(m, n, v) }
}

#[export_name = "lonlat2x_"]
pub unsafe extern "C" fn lonlat2x_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { proj::lonlat2x(m, n, v) }
}

#[export_name = "lonlat2y_"]
pub unsafe extern "C" fn lonlat2y_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { proj::lonlat2y(m, n, v) }
}

#[export_name = "__iceproperties_MOD_geticeconductivity"]
pub unsafe extern "C" fn mod_geticeconductivity_entry(temp: *const c_double) -> c_double {
    unsafe { iceproperties::mod_get_ice_conductivity(*temp) }
}

#[export_name = "__iceproperties_MOD_geticecapacity"]
pub unsafe extern "C" fn mod_geticecapacity_entry(temp: *const c_double) -> c_double {
    unsafe { iceproperties::mod_get_ice_capacity(*temp) }
}

#[export_name = "__iceproperties_MOD_geticepressuremeltingpoint"]
pub unsafe extern "C" fn mod_geticepressuremeltingpoint_entry(
    cc: *const c_double,
    press: *const c_double,
) -> c_double {
    unsafe { iceproperties::mod_get_ice_pressure_melting_point(*cc, *press) }
}

#[export_name = "iceconductivity_"]
pub unsafe extern "C" fn iceconductivity_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { iceproperties::ice_conductivity(m, n, v) }
}

#[export_name = "icecapacity_"]
pub unsafe extern "C" fn icecapacity_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { iceproperties::ice_capacity(m, n, v) }
}

#[export_name = "icepressuremeltingpoint_"]
pub unsafe extern "C" fn icepressuremeltingpoint_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { iceproperties::ice_pressure_melting_point(m, n, v) }
}

#[export_name = "relativetemperature_"]
pub unsafe extern "C" fn relativetemperature_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { iceproperties::relative_temperature(m, n, v) }
}

#[export_name = "arrheniusfactor_"]
pub unsafe extern "C" fn arrheniusfactor_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { iceproperties::arrhenius_factor(m, n, v) }
}

#[export_name = "calculate_haf_"]
pub unsafe extern "C" fn calculate_haf_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { haf::calculate_haf(m, n, v) }
}

#[export_name = "epltoids_"]
pub unsafe extern "C" fn epltoids_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { watertransfer::epl_to_ids(m, n, v) }
}

#[export_name = "sourcecalc_"]
pub unsafe extern "C" fn sourcecalc_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { srccalc::source_calc(m, n, v) }
}

#[export_name = "enhancementfactor_"]
pub unsafe extern "C" fn enhancementfactor_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { damage::enhancement_factor(m, n, v) }
}

#[export_name = "sourcedamage_"]
pub unsafe extern "C" fn sourcedamage_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { damage::source_damage(m, n, v) }
}

#[export_name = "getfrictionheat_"]
pub unsafe extern "C" fn getfrictionheat_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { friction::friction_heat(m, n, v) }
}

#[export_name = "getfrictionloads_"]
pub unsafe extern "C" fn getfrictionloads_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { friction::friction_loads(m, n, v) }
}

#[export_name = "sliding_weertman_"]
pub unsafe extern "C" fn sliding_weertman_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { sliding::sliding_weertman(m, n, v) }
}

#[export_name = "friction_coulomb_"]
pub unsafe extern "C" fn friction_coulomb_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { sliding::friction_coulomb(m, n, v) }
}

#[export_name = "sliding_budd_"]
pub unsafe extern "C" fn sliding_budd_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { sliding::sliding_budd(m, n, v) }
}

#[export_name = "effectivepressure_"]
pub unsafe extern "C" fn effectivepressure_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { sliding::effective_pressure(m, n, v) }
}

#[export_name = "freeslipshelves_"]
pub unsafe extern "C" fn freeslipshelves_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { sliding::free_slip_shelves(m, n, v) }
}

#[export_name = "slidcoef_contact_"]
pub unsafe extern "C" fn slidcoef_contact_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { contact::slidcoef_contact(m, n, v) }
}

#[export_name = "gzept_scaling_"]
pub unsafe extern "C" fn gzept_scaling_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { contact::gzept_scaling(m, n, v) }
}

#[export_name = "lateralfriction_x_"]
pub unsafe extern "C" fn lateralfriction_x_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { lateral::lateral_friction_x(m, n, v) }
}

#[export_name = "lateralfriction_y_"]
pub unsafe extern "C" fn lateralfriction_y_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { lateral::lateral_friction_y(m, n, v) }
}

#[export_name = "lateralfriction_"]
pub unsafe extern "C" fn lateralfriction_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
    axis: *const c_int,
) -> c_double {
    // 4-arg form `(Model, nodenumber, x, axis)`: `axis` crosses by
    // reference as its own dummy. Only the `_x`/`_y` wrappers are
    // referenced from `.sif` files; this symbol exists for surface
    // equality and behaves per the Fortran calling convention.
    let x = unsafe { *v };
    let node = unsafe { *n };
    let a = unsafe { *axis };
    lateral::lateral_friction(m, node, x, a)
}

#[export_name = "horizontalvelo_"]
pub unsafe extern "C" fn horizontalvelo_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { coupling::horizontal_velo(m, n, v) }
}

#[export_name = "overburdenpressure_"]
pub unsafe extern "C" fn overburdenpressure_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { coupling::overburden_pressure(m, n, v) }
}

#[export_name = "seapressure_"]
pub unsafe extern "C" fn seapressure_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { buoyancy::sea_pressure(m, n, v) }
}

#[export_name = "seaspring_"]
pub unsafe extern "C" fn seaspring_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { buoyancy::sea_spring(m, n, v) }
}

#[export_name = "caffegetviscosity_"]
pub unsafe extern "C" fn caffegetviscosity_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { caffe::caffe_get_viscosity(m, n, v) }
}

#[export_name = "shapefactorgravity_x_"]
pub unsafe extern "C" fn shapefactorgravity_x_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { shape::shape_factor_gravity_x(m, n, v) }
}

#[export_name = "shapefactorgravity_y_"]
pub unsafe extern "C" fn shapefactorgravity_y_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { shape::shape_factor_gravity_y(m, n, v) }
}

#[export_name = "shapefactorgravity_"]
pub unsafe extern "C" fn shapefactorgravity_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
    axis: *const c_int,
) -> c_double {
    // 4-arg form `(Model, nodenumber, x, axis)`; only the `_x`/`_y`
    // wrappers are referenced from `.sif` files.
    let x = unsafe { *v };
    let node = unsafe { *n };
    let a = unsafe { *axis };
    shape::shape_factor_gravity(m, node, x, a)
}

#[export_name = "zsini_"]
pub unsafe extern "C" fn zsini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_ini(m, n, v) }
}

#[export_name = "zsmzsini_"]
pub unsafe extern "C" fn zsmzsini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_mzs_ini(m, n, v) }
}

#[export_name = "zstopini_"]
pub unsafe extern "C" fn zstopini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_top_ini(m, n, v) }
}

#[export_name = "zstopmzsini_"]
pub unsafe extern "C" fn zstopmzsini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_top_mzs_ini(m, n, v) }
}

#[export_name = "zsbottomini_"]
pub unsafe extern "C" fn zsbottomini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_bottom_ini(m, n, v) }
}

#[export_name = "zsbottommzsini_"]
pub unsafe extern "C" fn zsbottommzsini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_bottom_mzs_ini(m, n, v) }
}

#[export_name = "dyini_"]
pub unsafe extern "C" fn dyini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::dy_ini(m, n, v) }
}

#[export_name = "dymdyini_"]
pub unsafe extern "C" fn dymdyini_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::dy_mdy_ini(m, n, v) }
}

#[export_name = "zsbottommzsinicalving_"]
pub unsafe extern "C" fn zsbottommzsinicalving_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_bottom_mzs_ini_calving(m, n, v) }
}

#[export_name = "zstopmzsinicalving_"]
pub unsafe extern "C" fn zstopmzsinicalving_entry(
    m: *mut c_void,
    n: *const c_int,
    v: *const c_double,
) -> c_double {
    unsafe { zs::zs_top_mzs_ini_calving(m, n, v) }
}
