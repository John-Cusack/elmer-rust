//! MATC `$...$` evaluation through the vendored C MATC library.
//!
//! Mirrors `matc_c` in Elmer's `fem/src/Load.c`: one global init
//! (`mtc_init` + `format(12,"rowform")`), result with the trailing newline
//! stripped, and hard errors surfaced instead of `exit(1)`.
//! All calls are serialized behind a process-wide mutex; the C library is
//! built without OpenMP so its globals are process-global.

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::sync::{LazyLock, Mutex};

static LOCK: LazyLock<Mutex<()>> = LazyLock::new(|| Mutex::new(()));

unsafe extern "C" {
    fn mtc_init(input: *mut c_char, output: *mut c_char, error: *mut c_char);
    fn mtc_domath(cmd: *mut c_char) -> *mut c_char;
}


fn init_once() {
    static DONE: std::sync::Once = std::sync::Once::new();
    DONE.call_once(|| {
        // Exactly Load.c:mtc_init_once (stdout/stderr inherited).
        unsafe {
            mtc_init(std::ptr::null_mut(), c_stdout(), c_stderr());
            let fmt = CString::new("format( 12,\"rowform\")").unwrap();
            // `mtc_domath` takes `char *` (may write); pass an owned copy.
            let mut owned = fmt.into_bytes_with_nul();
            mtc_domath(owned.as_mut_ptr() as *mut c_char);
        }
    });
}

// stdout/stderr FILE* without depending on the `libc` crate.
// glibc exports these symbols; other platforms get NULL (MATC then
// cannot print, but `mtc_domath` results are returned, not printed).
#[cfg(target_os = "linux")]
fn c_stdout() -> *mut c_char {
    unsafe extern "C" {
        static stdout: *mut c_char;
    }
    unsafe { stdout }
}

#[cfg(target_os = "linux")]
fn c_stderr() -> *mut c_char {
    unsafe extern "C" {
        static stderr: *mut c_char;
    }
    unsafe { stderr }
}

#[cfg(not(target_os = "linux"))]
fn c_stdout() -> *mut c_char {
    std::ptr::null_mut()
}

#[cfg(not(target_os = "linux"))]
fn c_stderr() -> *mut c_char {
    std::ptr::null_mut()
}

/// Evaluate one MATC command string (the text between `$` delimiters).
/// Returns the printed result with the trailing linefeed removed,
/// or an error when MATC reports one.
pub fn matc_eval(cmd: &str) -> Result<String, String> {
    let _guard = LOCK.lock().unwrap();
    init_once();
    // Mirror matc_c: a leading "nc:" selects the result without computing.
    let body = cmd.strip_prefix("nc:").unwrap_or(cmd);
    let owned = CString::new(body).map_err(|e| format!("MATC ERROR: NUL byte: {e}"))?;
    let mut buf = owned.into_bytes_with_nul();
    let out = unsafe { mtc_domath(buf.as_mut_ptr() as *mut c_char) };
    if out.is_null() {
        return Ok(String::new());
    }
    let s = unsafe { CStr::from_ptr(out) }.to_string_lossy().into_owned();
    // matc_c: slen = strlen(ptr)-1 ("ignore linefeed").
    let s = s.strip_suffix('\n').unwrap_or(&s);
    if s.starts_with("MATC ERROR:") || s.starts_with("WARNING:") {
        return Err(s.to_string());
    }
    Ok(s.to_string())
}

/// Evaluate with the exact `matc_c` empty-result convention: NULL or empty
/// yields the empty string (never an error).
pub fn matc_eval_lossy(cmd: &str) -> String {
    matc_eval(cmd).unwrap_or_else(|e| e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        assert_eq!(matc_eval("1+2*3").unwrap(), "7");
    }

    #[test]
    fn range_rowform() {
        // `format(12,"rowform")` prints vectors as space-separated rows.
        assert_eq!(matc_eval("1:6").unwrap(), "1 2 3 4 5 6");
    }

    #[test]
    fn error_reported() {
        assert!(matc_eval("1+*2").is_err());
    }
}
