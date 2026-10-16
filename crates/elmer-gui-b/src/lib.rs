//! lane-gui-b: grade record for 123 `ElmerGUI/netgen/` C++ units (sorted-id half 2).
//!
//! Verdict: **123/123 skipped, 0 ported.** Every compilable lane TU exports
//! exclusively C++-mangled (`_Z*`) symbols — Rust cannot provide those without
//! the C++ ABI, so per the lane brief they are recorded skips, never forced.
//! Headers carry no translation unit; three TUs do not build in-tree and are
//! excluded from upstream's own `ng` library target.
//!
//! This crate is the machine-checkable form of that verdict: the [`UNITS`]
//! table (generated from the proof bundle's `units.json`) plus tests that pin
//! the counts, the closed reason vocabulary, and the gate digests. It emits
//! no symbols by design — a port that provided a symbol would contradict the
//! zero-STRONG-dup gate recorded in `evidence/gate-output.txt`.

/// Campaign lane id (rustsmith `data/elmer-lanes/lanes.json`).
pub const LANE_ID: &str = "lane-gui-b";
/// Pinned upstream checkout the proof was measured against.
pub const UPSTREAM_PIN: &str = "9f6af2f85";
/// Frozen unit count for this lane.
pub const UNIT_COUNT: usize = 123;

/// One lane unit's grade record.
pub struct Unit {
    /// Canonical unit id, e.g. `cxx:ElmerGUI/netgen/libsrc/linalg/vector.cpp`.
    pub id: &'static str,
    /// Closed vocabulary: `cxx_mangled_only` | `header_no_tu` | `not_built_in_tree`.
    pub reason: &'static str,
    /// Human-readable evidence summary.
    pub detail: &'static str,
    /// Measured `g++ -c` wall seconds (`None` for non-compiled units).
    pub compile_wall_s: Option<f64>,
    /// SHA-256 of the compiled object (`""` for non-compiled units).
    pub object_sha256: &'static str,
}

include!("units_gen.rs");

/// Look up a unit by canonical id.
pub fn find(id: &str) -> Option<&'static Unit> {
    UNITS.iter().find(|u| u.id == id)
}

/// Count units per skip reason: `(cxx_mangled_only, header_no_tu, not_built_in_tree)`.
pub fn reason_counts() -> (usize, usize, usize) {
    let mut counts = (0, 0, 0);
    for u in UNITS {
        match u.reason {
            "cxx_mangled_only" => counts.0 += 1,
            "header_no_tu" => counts.1 += 1,
            "not_built_in_tree" => counts.2 += 1,
            _ => {}
        }
    }
    counts
}
