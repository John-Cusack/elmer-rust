//! `xy2lon_`, `xy2lat_`, `lonlat2x_`, `lonlat2y_` — port of USF_proj.F90
//! (geographic projections via `ProjUtils` under `elmerice/Utils`).
//!
//! Faithfulness notes:
//! - The projection math itself stays in the reference domain: each entry
//!   splits its input exactly like the source (`x=VarIn(1)` etc.) and
//!   calls the reference `xy2LonLat`/`LonLat2xy` through the Fortran
//!   accessor, so the port owns only the entry frame, never the math.
//! - `Model`/`nodenumber` are received and ignored exactly as in the
//!   source.

use crate::support::*;

#[inline]
unsafe fn pair(v: *const c_double) -> (c_double, c_double) {
    // SAFETY: same assumed-size contract as cov::slot; entries declare
    // `VarIn(2)`.
    unsafe { (*v, *v.add(1)) }
}

pub unsafe fn xy2lon(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let (x, y) = unsafe { pair(v) };
    let (mut lon, mut lat) = (0.0, 0.0);
    unsafe { usf_xy2lonlat(x, y, &mut lon, &mut lat) };
    lon
}

pub unsafe fn xy2lat(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let (x, y) = unsafe { pair(v) };
    let (mut lon, mut lat) = (0.0, 0.0);
    unsafe { usf_xy2lonlat(x, y, &mut lon, &mut lat) };
    lat
}

pub unsafe fn lonlat2x(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let (lon, lat) = unsafe { pair(v) };
    let (mut x, mut y) = (0.0, 0.0);
    unsafe { usf_lonlat2xy(lon, lat, &mut x, &mut y) };
    x
}

pub unsafe fn lonlat2y(_m: *mut c_void, _n: *const c_int, v: *const c_double) -> c_double {
    let (lon, lat) = unsafe { pair(v) };
    let (mut x, mut y) = (0.0, 0.0);
    unsafe { usf_lonlat2xy(lon, lat, &mut x, &mut y) };
    y
}
