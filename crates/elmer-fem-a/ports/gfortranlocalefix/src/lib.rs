// Port of `fem/src/GFortranLocaleFix.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/GFortranLocaleFix.c. True object symbol:
// `elmer_fix_numeric_locale` (the `_WIN32` import-table patching compiles
// out on Linux; only the setlocale(LC_NUMERIC, "C") path builds).
//
// Idempotence uses the same once-flag semantics (C `static int done`);
// the port uses an AtomicBool so concurrent first calls stay sound.

use std::os::raw::{c_char, c_int};
use std::sync::atomic::{AtomicBool, Ordering};

#[link(name = "c")]
extern "C" {
    fn setlocale(category: c_int, locale: *const c_char) -> *mut c_char;
}

const LC_NUMERIC: c_int = 1;

static DONE: AtomicBool = AtomicBool::new(false);

#[export_name = "elmer_fix_numeric_locale"]
pub extern "C" fn port_elmer_fix_numeric_locale() {
    if DONE.swap(true, Ordering::SeqCst) {
        return;
    }
    unsafe {
        setlocale(LC_NUMERIC, c"C".as_ptr());
    }
}
