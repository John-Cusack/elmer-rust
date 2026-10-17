// Port of `fem/src/binio/binio.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/binio/binio.c. True object symbols: `binendianess_c`,
// `binsetinputendianess_c`, `binopen_c`, `binclose_c`, `binwriteint4_c`,
// `binwriteint8_c`, `binreadint8_c`, `binwritedouble_c`, `binreadint4_c`,
// `binreaddouble_c`, `binwritestring_c`, `binreadstring_c`,
// `binwritechar_c`, `binftell_c`, `binfseek_c`, `strerrorf_c`
// (`-DUSE_ISO_C_BINDINGS` plain names; `HAVE_FTELLO`/`HAVE_FSEEKO` take
// the off_t branches).
//
// Faithfulness notes:
// - Same libc stdio calls through extern "C" (fopen/fread/fwrite/...),
//   same errno values via __errno_location, same strerror() text.
// - `assert(units[u].fd)` (active in the Debug grade build): replicated
//   as abort-after-message. The message names this file/function the way
//   glibc names the C translation unit; assert paths are unreachable in
//   every probe (callers hold valid units) — recorded in EVIDENCE.md.
// - `units[]` OOB would be UB in C; the port indexes a fixed 1000-entry
//   table the same way (callers use small unit numbers).

use std::os::raw::{c_char, c_double, c_int, c_long, c_ulong, c_void};

type OffT = c_long;
type SizeT = c_ulong;

const MAX_UNITS: usize = 1000;

#[derive(Clone, Copy)]
#[repr(C)]
struct Unit {
    convert: c_int,
    fd: *mut c_void,
}

static mut UNITS: [Unit; MAX_UNITS] = {
    const INIT: Unit = Unit { convert: 0, fd: std::ptr::null_mut() };
    [INIT; MAX_UNITS]
};

#[link(name = "c")]
extern "C" {
    fn fopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn fclose(fp: *mut c_void) -> c_int;
    fn fread(ptr: *mut c_void, size: SizeT, n: SizeT, fp: *mut c_void) -> SizeT;
    fn fwrite(ptr: *const c_void, size: SizeT, n: SizeT, fp: *mut c_void) -> SizeT;
    fn fgetc(fp: *mut c_void) -> c_int;
    fn fputc(c: c_int, fp: *mut c_void) -> c_int;
    fn feof(fp: *mut c_void) -> c_int;
    fn ferror(fp: *mut c_void) -> c_int;
    fn ftello(fp: *mut c_void) -> OffT;
    fn fseeko(fp: *mut c_void, off: OffT, whence: c_int) -> c_int;
    fn strerror(e: c_int) -> *const c_char;
    fn malloc(n: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn abort() -> !;
    fn __errno_location() -> *mut c_int;
}

unsafe fn errno_val() -> c_int {
    *__errno_location()
}

unsafe fn c_assert(fd_ok: bool, func: &str, line: u32, expr: &str) {
    if !fd_ok {
        eprint!(
            "binio: binio.c:{}: {}: Assertion `{}` failed.\n",
            line, func, expr
        );
        abort();
    }
}

fn is_little() -> bool {
    u16::from_ne_bytes([0x22, 0x11]) == 0x1122
}

fn endianess() -> c_char {
    if is_little() { b'L' as c_char } else { b'B' as c_char }
}

fn swap_bytes(o: *mut u8, n: usize) {
    unsafe {
        for i in 0..n / 2 {
            let t = *o.add(i);
            *o.add(i) = *o.add(n - 1 - i);
            *o.add(n - 1 - i) = t;
        }
    }
}

#[export_name = "binendianess_c"]
pub unsafe extern "C" fn port_binendianess_c(e: *mut c_char) {
    *e = endianess();
}

#[export_name = "binsetinputendianess_c"]
pub unsafe extern "C" fn port_binsetinputendianess_c(unit: *const c_int, e: *const c_char) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binsetinputendianess_c", 95, "units[*unit].fd");
    UNITS[u].convert = if *e == endianess() { 0 } else { 1 };
}

#[export_name = "binopen_c"]
pub unsafe extern "C" fn port_binopen_c(
    unit: *const c_int,
    file: *const c_char,
    file_len: *const c_int,
    action: *const c_char,
    status: *mut c_int,
) {
    let u = *unit as usize;
    let fl = *file_len as usize;
    let fname = malloc(fl + 1) as *mut c_char;
    std::ptr::copy_nonoverlapping(file, fname, fl);
    *fname.add(fl) = 0;
    let mode: *const c_char = match *action as u8 {
        b'w' | b'W' => c"wb".as_ptr(),
        b'a' | b'A' => c"ab".as_ptr(),
        _ => c"rb".as_ptr(),
    };
    c_assert(UNITS[u].fd.is_null(), "binopen_c", 118, "!units[*unit].fd");
    UNITS[u].fd = fopen(fname, mode);
    free(fname as *mut c_void);
    if UNITS[u].fd.is_null() {
        *status = errno_val();
        return;
    }
    *status = 0;
}

#[export_name = "binclose_c"]
pub unsafe extern "C" fn port_binclose_c(unit: *const c_int, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binclose_c", 136, "units[*unit].fd");
    let s = fclose(UNITS[u].fd);
    UNITS[u].fd = std::ptr::null_mut();
    *status = if s != 0 { errno_val() } else { 0 };
}

#[export_name = "binwriteint4_c"]
pub unsafe extern "C" fn port_binwriteint4_c(unit: *const c_int, n: *const u32, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binwriteint4_c", 146, "units[*unit].fd");
    let s = fwrite(n as *const c_void, 1, 4, UNITS[u].fd);
    *status = if s == 4 { 0 } else { errno_val() };
}

#[export_name = "binwriteint8_c"]
pub unsafe extern "C" fn port_binwriteint8_c(unit: *const c_int, n: *const u64, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binwriteint8_c", 157, "units[*unit].fd");
    let s = fwrite(n as *const c_void, 1, 8, UNITS[u].fd);
    *status = if s == 8 { 0 } else { errno_val() };
}

#[export_name = "binreadint8_c"]
pub unsafe extern "C" fn port_binreadint8_c(unit: *const c_int, n: *mut u64, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binreadint8_c", 168, "units[*unit].fd");
    let r = fread(n as *mut c_void, 1, 8, UNITS[u].fd);
    if r != 8 {
        *status = if feof(UNITS[u].fd) != 0 { -1 } else { errno_val() };
    } else {
        *status = 0;
    }
    if UNITS[u].convert != 0 {
        swap_bytes(n as *mut u8, 8);
    }
}

#[export_name = "binwritedouble_c"]
pub unsafe extern "C" fn port_binwritedouble_c(unit: *const c_int, a: *const c_double, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binwritedouble_c", 184, "units[*unit].fd");
    let s = fwrite(a as *const c_void, 1, 8, UNITS[u].fd);
    *status = if s == 8 { 0 } else { errno_val() };
}

#[export_name = "binreadint4_c"]
pub unsafe extern "C" fn port_binreadint4_c(unit: *const c_int, n: *mut u32, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binreadint4_c", 195, "units[*unit].fd");
    let r = fread(n as *mut c_void, 1, 4, UNITS[u].fd);
    if r != 4 {
        *status = if feof(UNITS[u].fd) != 0 { -1 } else { errno_val() };
    } else {
        *status = 0;
    }
    if UNITS[u].convert != 0 {
        swap_bytes(n as *mut u8, 4);
    }
}

#[export_name = "binreaddouble_c"]
pub unsafe extern "C" fn port_binreaddouble_c(unit: *const c_int, a: *mut c_double, status: *mut c_int) {
    port_binreadint8_c(unit, a as *mut u64, status);
}

#[export_name = "binwritestring_c"]
pub unsafe extern "C" fn port_binwritestring_c(
    unit: *const c_int,
    s: *const c_char,
    s_len: *const c_int,
    status: *mut c_int,
) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binwritestring_c", 218, "units[*unit].fd");
    if fwrite(s as *const c_void, 1, *s_len as SizeT, UNITS[u].fd) == *s_len as SizeT
        && fputc(0, UNITS[u].fd) == 0
    {
        *status = 0;
    } else {
        *status = errno_val();
    }
}

#[export_name = "binreadstring_c"]
pub unsafe extern "C" fn port_binreadstring_c(
    unit: *const c_int,
    s: *mut c_char,
    s_len: *const c_int,
    status: *mut c_int,
) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binreadstring_c", 234, "units[*unit].fd");
    let sl = *s_len as isize;
    let mut i: isize = 0;
    let mut c: c_int = 0;
    while i < sl {
        c = fgetc(UNITS[u].fd);
        if c == 0 || c == -1 {
            break;
        }
        *s.offset(i) = c as c_char;
        i += 1;
    }
    while i < sl {
        *s.offset(i) = b' ' as c_char;
        i += 1;
    }
    if c == -1 {
        *status = if ferror(UNITS[u].fd) != 0 { errno_val() } else { -1 };
    } else {
        *status = 0;
    }
}

#[export_name = "binwritechar_c"]
pub unsafe extern "C" fn port_binwritechar_c(unit: *const c_int, c: *mut c_char, status: *mut c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binwritechar_c", 254, "units[*unit].fd");
    if fwrite(c as *const c_void, 1, 1, UNITS[u].fd) == 1 {
        *status = 0;
    } else {
        *status = errno_val();
    }
}

#[export_name = "binftell_c"]
pub unsafe extern "C" fn port_binftell_c(unit: *const c_int) -> OffT {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binftell_c", 266, "units[*unit].fd");
    ftello(UNITS[u].fd)
}

#[export_name = "binfseek_c"]
pub unsafe extern "C" fn port_binfseek_c(unit: *const c_int, offset: *const OffT, whence: *const c_int) {
    let u = *unit as usize;
    c_assert(!UNITS[u].fd.is_null(), "binfseek_c", 281, "units[*unit].fd");
    match *whence {
        0 => {
            fseeko(UNITS[u].fd, *offset, 0);
        }
        1 => {
            fseeko(UNITS[u].fd, *offset, 1);
        }
        2 => {
            fseeko(UNITS[u].fd, *offset, 2);
        }
        _ => {}
    }
}

#[export_name = "strerrorf_c"]
pub unsafe extern "C" fn port_strerrorf_c(e: *const c_int, s: *mut c_char, s_len: *const c_int) {
    let t = strerror(*e);
    let sl = *s_len as isize;
    let mut i: isize = 0;
    while i < sl && *t.offset(i) != 0 {
        *s.offset(i) = *t.offset(i);
        i += 1;
    }
    while i < sl {
        *s.offset(i) = b' ' as c_char;
        i += 1;
    }
}


