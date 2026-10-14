// Port of `fem/src/view3d/second.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/second.c. True object symbol: `second`.
//
// The non-Windows branch calls gettimeofday() directly (same syscall the
// port uses — no std-time truncation drift: `tv_sec + tv_usec/1.0e6` with
// identical double rounding).

use std::os::raw::{c_double, c_int, c_long, c_void};

#[repr(C)]
struct Timeval {
    tv_sec: c_long,
    tv_usec: c_long,
}

#[link(name = "c")]
extern "C" {
    fn gettimeofday(tv: *mut Timeval, tz: *mut c_void) -> c_int;
}

#[export_name = "second"]
pub extern "C" fn port_second() -> c_double {
    unsafe {
        let mut tp = Timeval { tv_sec: 0, tv_usec: 0 };
        gettimeofday(&mut tp, std::ptr::null_mut());
        tp.tv_sec as c_double + tp.tv_usec as c_double / 1.0e6
    }
}
