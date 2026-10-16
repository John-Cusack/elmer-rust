//! lane-gui-a: scope record for `ElmerGUI/netgen/` sorted-id half 1.
//!
//! 124 units, **0 graded / 124 recorded skips**. Every unit's defined
//! surface is C++-mangled (`_Z*`), declaration-only, or dead upstream, so
//! no Rust substitution is possible (Rust cannot provide C++-mangled
//! symbols). Skips are recorded with reasons in [`scope::UNITS`]'s source
//! file `scope/units.json`, never forced.
//!
//! Method (batched fixed cost, one compile+`nm` pass for the whole lane):
//! each of the 50 `.cpp` units was compiled with
//! `g++ -DNO_PARALLEL_THREADS -IElmerGUI/netgen/libsrc/include -c` and its
//! object classified with `nm -g --defined-only`: 48 define zero C-linkage
//! symbols (all `_Z` mangled or toolchain refs), 2 do not compile upstream
//! (bare-text file head, absent from the `ng` library `SOURCES`).
//! The 70 `.hpp` units are C++ headers with no translation unit; the 4
//! `.h` units are declaration-only (the `extern "C"` `Ng_*` definitions
//! live in `interface/*.cpp`, owned by lane-gui-b).

pub mod scope {
    /// Lane id this crate records.
    pub const LANE: &str = "lane-gui-a";
    /// Frozen partition size for the lane (must equal registry length).
    pub const UNIT_COUNT: usize = 124;
    /// Expected per-reason counts: (reason, count).
    pub const REASON_COUNTS: &[(&str, usize)] = &[
        ("cxx-mangled-only", 48),
        ("cxx-header-only", 70),
        ("c-header-decl-only", 4),
        ("dead-upstream-uncompilable", 2),
    ];
}
