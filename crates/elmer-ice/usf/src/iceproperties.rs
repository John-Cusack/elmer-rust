//! `geticeconductivity`, `geticecapacity`, `geticepressuremeltingpoint`
//! (module `IceProperties`), `iceconductivity_`, `icecapacity_`,
//! `icepressuremeltingpoint_`, `relativetemperature_`, `arrheniusfactor_`
//! — port of USF_IceProperties.F90 (material properties of ice).
//!
//! Faithfulness notes (mirrors the source statement for statement):
//! - Module functions are pure arithmetic transcriptions (`9.828*exp`,
//!   `146.3+7.253*t`, `273.15-CC*MAX(press,0)`); the `__iceproperties_MOD_*`
//!   symbols are exported for nm equality and execute the same bodies.
//! - File-scope entries read Element/Material/BC/Constants through the
//!   accessors, with identical fallback chains (scaling 1.0, CC 9.8e-8,
//!   Gas 8.314, Glen 3.0, limits, rate factors, activation energies,
//!   enhancement 1.0) and identical `INFO` calls at identical levels.
//! - `SAVE FirstTime/ClausiusClapeyron` replicated per entry as a Rust
//!   static (like the CS port's SAVE replication); the dead second
//!   `IF (FirstTime)` check in `RelativeTemperature` (already `.FALSE.`)
//!   is kept verbatim and documented, never "fixed".
//! - `GetLogical(Material, ...)` runs outside any ASSOCIATED guard exactly
//!   like the source: a null Material handle faults identically in both
//!   arms (pass-through, never null-guarded).
//! - No Rust threads: entries run serially; each `Mutex` is only a
//!   static-sharing guard.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

use crate::support::*;

// -- module IceProperties (exported for nm equality) --------------------

pub unsafe fn mod_get_ice_conductivity(temp: c_double) -> c_double {
    9.828 * (-5.7e-03 * temp).exp()
}

pub unsafe fn mod_get_ice_capacity(temp: c_double) -> c_double {
    146.3 + (7.253 * temp)
}

pub unsafe fn mod_get_ice_pressure_melting_point(cc: c_double, press: c_double) -> c_double {
    273.15 - cc * press.max(0.0)
}

// -- file-scope entries --------------------------------------------------

struct CcState {
    first: bool,
    cc: c_double,
}

pub unsafe fn ice_conductivity(m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let temp = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let mut found = false;
    let mut scaling = unsafe { const_real(h_mat, "Heat Conductivity Scaling Factor", &mut found) };
    if !found {
        scaling = 1.0;
    }
    scaling * unsafe { mod_get_ice_conductivity(temp) }
}

pub unsafe fn ice_capacity(m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let temp = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let mut found = false;
    let mut scaling = unsafe { const_real(h_mat, "Heat Capacity Scaling Factor", &mut found) };
    if !found {
        scaling = 1.0;
    }
    scaling * unsafe { mod_get_ice_capacity(temp) }
}

static IPMP_STATE: Mutex<CcState> = Mutex::new(CcState {
    first: true,
    cc: 0.0,
});

pub unsafe fn ice_pressure_melting_point(
    m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let press = unsafe { *v };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let mut found = false;
    let mut scaling: c_double;
    let mut h_bc: *const c_void = std::ptr::null();
    // Material => GetMaterial(Element); IF (.NOT.ASSOCIATED(Material)) ...
    // (null handle pass-through: ASSOCIATED(NULL) is false in both arms)
    let has_mat = !h_mat.is_null();
    if !has_mat {
        unsafe { usf_get_bc(h_elem, &mut h_bc) };
        if h_bc.is_null() {
            scaling = 1.0;
            unsafe {
                info(
                    "IcePressureMeltingPoint",
                    "No \"Pressure Scaling Factor\" found - setting to 1.0",
                    9,
                )
            };
            found = false;
        } else {
            scaling = unsafe { const_real(h_bc, "Pressure Scaling Factor", &mut found) };
        }
    } else {
        scaling = unsafe { const_real(h_mat, "Pressure Scaling Factor", &mut found) };
    }
    if !found {
        scaling = 1.0;
    }
    let in_celsius = unsafe { get_logical(h_mat, "Pressure Melting Point Celsius", &mut found) };
    let tmpoffset = if in_celsius { 273.15 } else { 0.0 };
    let cc = {
        let mut st = IPMP_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let mut h_c: *const c_void = std::ptr::null();
            unsafe { usf_get_constants(&mut h_c) };
            let mut got = false;
            if !h_c.is_null() {
                st.cc = unsafe { const_real(h_c, "Clausius Clapeyron Constant", &mut got) };
            }
            if h_c.is_null() || !got {
                st.cc = 9.8e-08;
                unsafe {
                    info(
                        "IcePressureMeltingPoint",
                        "No entry found for >Clausius Clapeyron Constant<.",
                        9,
                    )
                };
                unsafe {
                    info(
                        "IcePressureMeltingPoint",
                        "Setting to 9.8d-08 (SI units)",
                        9,
                    )
                };
            }
        }
        st.cc
    };
    unsafe { mod_get_ice_pressure_melting_point(cc, scaling * press) + tmpoffset }
}

static RT_STATE: Mutex<CcState> = Mutex::new(CcState {
    first: true,
    cc: 0.0,
});

pub unsafe fn relative_temperature(
    m: *mut c_void,
    _n: *const c_int,
    v: *const c_double,
) -> c_double {
    let (tabs, press) = unsafe { (*v, *v.add(1)) };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let mut found = false;
    let mut scaling: c_double;
    let mut h_bc: *const c_void = std::ptr::null();
    let has_mat = !h_mat.is_null();
    if !has_mat {
        unsafe { usf_get_bc(h_elem, &mut h_bc) };
        if h_bc.is_null() {
            scaling = 1.0;
            unsafe {
                info(
                    "IcePressureMeltingPoint",
                    "No \"Pressure Scaling Factor\" found - setting to 1.0",
                    9,
                )
            };
            found = false;
        } else {
            scaling = unsafe { const_real(h_bc, "Pressure Scaling Factor", &mut found) };
        }
    } else {
        scaling = unsafe { const_real(h_mat, "Pressure Scaling Factor", &mut found) };
    }
    if !found {
        scaling = 1.0;
    }
    let cc = {
        let mut st = RT_STATE.lock().unwrap();
        if st.first {
            st.first = false;
            let mut h_c: *const c_void = std::ptr::null();
            unsafe { usf_get_constants(&mut h_c) };
            // CorrectRelTemp=ListGetLogical(Material,...) — outside any
            // ASSOCIATED guard in the source; pass-through verbatim.
            let mut g = false;
            let _ = unsafe { list_get_logical(h_mat, "Correct Relative Tempeature", &mut g) };
            let mut got = false;
            if !h_c.is_null() {
                st.cc = unsafe { const_real(h_c, "Clausius Clapeyron Constant", &mut got) };
            }
            if h_c.is_null() || !got {
                st.cc = 9.8e-08;
                unsafe {
                    info(
                        "IcePressureMeltingPoint",
                        "No entry found for >Clausius Clapeyron Constant<.",
                        9,
                    )
                };
                unsafe {
                    info(
                        "IcePressureMeltingPoint",
                        "Setting to 9.8d-08 (SI units)",
                        9,
                    )
                };
            }
        }
        st.cc
    };
    let in_celsius = unsafe { get_logical(h_mat, "Pressure Melting Point Celsius", &mut found) };
    let tmpoffset = if in_celsius { 273.15 } else { 0.0 };
    let tpmp = unsafe { mod_get_ice_pressure_melting_point(cc, scaling * press) + tmpoffset };
    let mut trel = tabs - tpmp;
    // Re-fetch pair, verbatim (source re-reads Element/Material here).
    let mut h_elem2: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem2) };
    let mut h_mat2: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem2, &mut h_mat2) };
    if !h_mat2.is_null() {
        let mut g = false;
        let correct = unsafe { list_get_logical(h_mat2, "Correct Relative Tempeature", &mut g) };
        if !g {
            // CorrectRelTemp stays .FALSE. (source default) — nothing to do.
        } else if correct {
            // Dead in the source (FirstTime already .FALSE.), kept verbatim.
            let st = RT_STATE.lock().unwrap();
            if st.first {
                drop(st);
                unsafe {
                    info(
                        "USF_IceProperties(IcePressureMeltingPoint)",
                        "Limiting Relative Temperature",
                        3,
                    )
                };
                RT_STATE.lock().unwrap().first = false;
            }
            trel = trel.min(0.0);
        }
    }
    trel
}

pub unsafe fn arrhenius_factor(m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let (input0, input1) = unsafe { (*v, *v.add(1)) };
    unsafe { usf_begin(m) };
    let mut h_elem: *const c_void = std::ptr::null();
    unsafe { usf_current_element(&mut h_elem) };
    let mut h_mat: *const c_void = std::ptr::null();
    unsafe { usf_get_material(h_elem, &mut h_mat) };
    let mut h_c: *const c_void = std::ptr::null();
    unsafe { usf_get_constants(&mut h_c) };
    let mut got = false;
    let mut gas = unsafe { const_real(h_c, "Gas Constant", &mut got) };
    if !got {
        gas = 8.314;
    }
    let mut cc = unsafe { const_real(h_c, "Clausius Clapeyron Constant", &mut got) };
    if !got {
        cc = 9.8e-08;
        unsafe {
            info(
                "IcePressureMeltingPoint",
                "No entry found for >Clausius Clapeyron Constant<.",
                5,
            )
        };
        unsafe {
            info(
                "IcePressureMeltingPoint",
                "Setting to 9.8d-08 (SI units)",
                5,
            )
        };
    }
    let mut scaling = unsafe { const_real(h_mat, "Pressure Scaling Factor", &mut got) };
    if !got {
        scaling = 1.0;
    }
    let in_celsius = unsafe { get_logical(h_mat, "Temperature in Celsius", &mut got) };
    let tempoffset = if in_celsius { 0.0 } else { 273.15 };
    let temp = input0 - tempoffset;
    let pressure = input1;
    let mut glen = unsafe { const_real(h_mat, "Glen Exponent", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusFactor)",
                "\"Glen Exponent\" not found. Setting to 3.0",
                5,
            )
        };
        glen = 3.0;
    }
    let mut lim = unsafe { const_real(h_mat, "Limit Temperature", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusFactor)",
                "\"Limit Temperature\" not found. Setting to -10.0",
                5,
            )
        };
        lim = -10.0;
    }
    let _ = lim;
    let mut r1 = unsafe { const_real(h_mat, "Rate Factor 1", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusFactor)",
                "\"Rate Factor 1\" not found. Setting to (SI value) 2.89165e-13",
                5,
            )
        };
        r1 = 2.89165e-13;
    }
    let mut r2 = unsafe { const_real(h_mat, "Rate Factor 2", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusFactor)",
                "\"Rate Factor 2\" not found. Setting to (SI value) 2.42736e-02",
                5,
            )
        };
        r2 = 2.42736e-02;
    }
    let scale_rf = unsafe { get_logical(h_mat, "Scale Rate Factors", &mut got) };
    if scale_rf {
        let f = unsafe { usf_pow_dp(scaling, glen) };
        r1 *= f;
        r2 *= f;
    }
    let mut e1 = unsafe { const_real(h_mat, "Activation Energy 1", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusEnergy)",
                "\"Activation Energy 1\" not found. Setting to (SI value) 60.0e03",
                5,
            )
        };
        e1 = 60.0e03;
    }
    let mut e2 = unsafe { const_real(h_mat, "Activation Energy 2", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusEnergy)",
                "\"Activation Energy 2\" not found. Setting to (SI value) 115.0e3",
                5,
            )
        };
        e2 = 115.0e3;
    }
    let mut enh = unsafe { const_real(h_mat, "Glen Enhancement Factor", &mut got) };
    if !got {
        unsafe {
            info(
                "IceProperties(ArrheniusEnergy)",
                "\"Glen Enhancement Factor\" not found. Setting to 1.0",
                5,
            )
        };
        enh = 1.0;
    }
    let _ = enh;
    let tpmp = unsafe { mod_get_ice_pressure_melting_point(cc, scaling * pressure) };
    // N.B. NOT `.clamp(-60.0, 0.0)`: with a NaN input the source's
    // MAX(MIN(x,0),-60) comparison chain yields 0.0 while `clamp`
    // returns NaN. Bit-identity requires the chain.
    #[allow(clippy::manual_clamp)]
    let trel = (temp - tpmp).min(0.0).max(-60.0);
    if trel < -10.0 {
        r1 * (-e1 / (gas * (273.15 + trel))).exp()
    } else {
        r2 * (-e2 / (gas * (273.15 + trel))).exp()
    }
}
