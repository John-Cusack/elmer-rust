// Port of `fem/src/minpack_kinds.F90` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: fortran:fem/src/minpack_kinds.F90.
//
// The module declares only `integer,parameter :: wp` (no procedures, no
// variables): the reference object defines no link symbols (verified with
// `nm -g --defined-only`, see EVIDENCE.md), so the port is intentionally
// symbol-free. `wp = SELECTED_REAL_KIND(12)` is f64; recorded here as a
// named constant for downstream ports.
pub const WP_BYTES: usize = 8;
