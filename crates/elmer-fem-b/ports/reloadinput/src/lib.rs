//! Port of `fem/src/modules/ReloadInput.F90` (Elmer @ 9f6af2f85, the lane
//! fork) for lane T021-B, at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/ReloadInput.F90` — one F77-global entry
//! `reloadinput_` with no internal procedures. The emitted dynamic symbol is
//! exactly `reloadinput_` (proven by `nm -D` on the reference `.so`).
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived-type
//! dummy (`Model`) crosses as an opaque address. `Solver`, `dt`, and
//! `TransientSimulation` are received and ignored, exactly as in the
//! original (it declares `Solver` but never references it, and never reads
//! `dt`/`TransientSimulation`). Every component access and every Elmer
//! library call goes through the generated-style Fortran `BIND(C)`
//! wrappers in `accessors/reloadinput_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no layout
//! assumption and no pinned offset. File I/O (`OPEN`/`READ`/`CLOSE` with
//! `ERR=` branches) and the `LoadInputFile` call stay in Fortran wrappers
//! so bit-identity holds by construction wherever the call sequence matches.
//!
//! Faithfulness notes (mirrors ReloadInput.F90 statement for statement):
//! - Branch structure replicates the `ERR=10` / `ERR=20` labels: no trigger
//!   file → first `OPEN` fails → second `OPEN` attempted → fails → return
//!   (early-return path, no output). Trigger present → `READ` the model
//!   file name, `CLOSE`, then `OPEN(InFileUnit, File=ModelName)`,
//!   `LoadInputFile(Model, InFileUnit, ModelName, MeshDir, MeshName,
//!   .FALSE., .FALSE.)`, `CLOSE(InFileUnit)`. The label-10 arm never calls
//!   `LoadInputFile` (it falls through to label 20 and returns), matching
//!   the original's dead second-open shape exactly.
//! - `ModelName`, `MeshDir`, `MeshName`, `IOUnit` live as wrapper-side SAVE
//!   state; each is overwritten on every path that observes it, so the
//!   persistence is unobservable versus the original's fresh locals.
//! - No Rust threads: the solver entries run serially, and the Rust side
//!   holds no state at all (all state is wrapper-side).

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};

extern "C" {
    fn rli_begin(h_model: *const c_void);
    fn rli_open1() -> c_int;
    fn rli_read1();
    fn rli_reload();
    fn rli_open2() -> c_int;
    fn rli_read2();
}

/// `reloadinput_` — the solver entry. Argument order and by-reference travel
/// match the Fortran dummy list `(Model, Solver, dt, TransientSimulation)`;
/// `Solver`, `dt`, and the transient flag are unused by the body, exactly as
/// in the original.
///
/// # Safety
///
/// `model` must be the Elmer `Model_t` address the loader passes.
/// `dt`/`transient` must be valid reads if non-null (they are ignored).
#[export_name = "reloadinput_"]
pub unsafe extern "C" fn reload_input_entry(
    model: *mut c_void,
    _solver: *mut c_void,
    _dt: *const c_double,
    _transient: *const c_int,
) {
    unsafe {
        rli_begin(model);
        if rli_open1() != 0 {
            // ERR=10 arm: retry the trigger-file open (label 10), which —
            // the file being absent — fails again into ERR=20 and returns.
            if rli_open2() != 0 {
                return;
            }
            // Label-10 success path: READ + CLOSE, then fall to label 20.
            rli_read2();
            return;
        }
        rli_read1();
        rli_reload();
    }
}
