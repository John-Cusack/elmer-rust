// Port of `fem/src/CPUTime.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/CPUTime.c. True object symbols: `cputime`, `realtime`,
// `cpumemory` (`-DUSE_ISO_C_BINDINGS` makes FC_FUNC the plain name; the
// POSIX branch builds on Linux).
//
// Same syscalls as the original (getrusage/gettimeofday), same double
// rounding: `ru_utime.tv_sec + tv_usec*1.0e-6`, `ru_maxrss` as-is.

use std::os::raw::{c_double, c_int, c_long, c_void};

#[repr(C)]
struct Timeval {
    tv_sec: c_long,
    tv_usec: c_long,
}

#[repr(C)]
struct Rusage {
    ru_utime: Timeval,
    ru_stime: Timeval,
    ru_maxrss: c_long,
    // Remaining fields are never read; the buffer only needs to cover
    // what getrusage writes (sizeof(struct rusage) == 144 on x86-64).
    _pad: [u8; 104],
}

#[link(name = "c")]
extern "C" {
    fn getrusage(who: c_int, usage: *mut Rusage) -> c_int;
    fn gettimeofday(tv: *mut Timeval, tz: *mut c_void) -> c_int;
}

const RUSAGE_SELF: c_int = 0;

#[export_name = "cputime"]
pub extern "C" fn port_cputime() -> c_double {
    unsafe {
        let mut usage = std::mem::zeroed::<Rusage>();
        getrusage(RUSAGE_SELF, &mut usage);
        usage.ru_utime.tv_sec as c_double + usage.ru_utime.tv_usec as c_double * 1.0e-6
    }
}

#[export_name = "realtime"]
pub extern "C" fn port_realtime() -> c_double {
    unsafe {
        let mut tp = Timeval { tv_sec: 0, tv_usec: 0 };
        gettimeofday(&mut tp, std::ptr::null_mut());
        tp.tv_sec as c_double + tp.tv_usec as c_double * 1.0e-6
    }
}

#[export_name = "cpumemory"]
pub extern "C" fn port_cpumemory() -> c_double {
    unsafe {
        let mut usage = std::mem::zeroed::<Rusage>();
        getrusage(RUSAGE_SELF, &mut usage);
        1.0 * usage.ru_maxrss as c_double
    }
}
