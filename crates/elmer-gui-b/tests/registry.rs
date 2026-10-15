//! Pins the lane-gui-b grade: 123/123 skipped, closed reason vocabulary,
//! per-class evidence invariants, and gate digests.

use elmer_gui_b::{reason_counts, find, LANE_ID, UNIT_COUNT, UNITS, UPSTREAM_PIN};
use std::collections::HashSet;

fn is_hex64(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit())
}

#[test]
fn lane_identity_and_size() {
    assert_eq!(LANE_ID, "lane-gui-b");
    assert_eq!(UPSTREAM_PIN, "9f6af2f85");
    assert_eq!(UNITS.len(), UNIT_COUNT);
    assert_eq!(UNIT_COUNT, 123);
}

#[test]
fn ids_unique_and_in_scope() {
    let mut seen = HashSet::new();
    for u in UNITS {
        assert!(seen.insert(u.id), "duplicate unit id: {}", u.id);
        let (lang, rel) = u.id.split_once(':').expect("unit id has lang prefix");
        assert_eq!(lang, "cxx");
        assert!(rel.starts_with("ElmerGUI/netgen/"), "out of lane scope: {u_id}", u_id = u.id);
    }
}

#[test]
fn closed_reason_vocabulary_and_counts() {
    for u in UNITS {
        assert!(
            matches!(u.reason, "cxx_mangled_only" | "header_no_tu" | "not_built_in_tree"),
            "open reason on {}: {}",
            u.id,
            u.reason
        );
        assert!(!u.detail.is_empty(), "empty detail on {}", u.id);
    }
    assert_eq!(reason_counts(), (69, 51, 3));
}

#[test]
fn per_class_evidence_invariants() {
    for u in UNITS {
        match u.reason {
            "cxx_mangled_only" => {
                assert!(u.compile_wall_s.is_some_and(|w| w > 0.0), "missing wall on {}", u.id);
                assert!(is_hex64(u.object_sha256), "bad object sha on {}", u.id);
            }
            "header_no_tu" => {
                assert!(u.id.ends_with(".hpp"), "non-header with header reason: {}", u.id);
                assert_eq!(u.compile_wall_s, None);
                assert_eq!(u.object_sha256, "");
            }
            "not_built_in_tree" => {
                assert!(u.id.ends_with(".cpp"), "non-TU with build reason: {}", u.id);
                assert!(u.compile_wall_s.is_some(), "missing wall on {}", u.id);
                assert_eq!(u.object_sha256, "");
            }
            _ => unreachable!(),
        }
    }
}

#[test]
fn spot_checks_match_measured_proof() {
    // nglib.cpp is the flagship "looks portable, is not" case: its Ng_* API
    // has C++ linkage (no extern "C" in nglib.h) — nm shows only _Z* globals.
    let nglib = find("cxx:ElmerGUI/netgen/libsrc/interface/nglib.cpp").expect("nglib row");
    assert_eq!(nglib.reason, "cxx_mangled_only");
    // nginterface.cpp holds real extern "C" Ng_* decls but cannot compile
    // in-tree (pruned visualization/ headers) and is excluded from the ng target.
    let ngi = find("cxx:ElmerGUI/netgen/libsrc/interface/nginterface.cpp").expect("nginterface row");
    assert_eq!(ngi.reason, "not_built_in_tree");
    assert!(ngi.detail.contains("visualization"));
    assert!(find("cxx:ElmerGUI/netgen/libsrc/interface/nglib.h").is_none(), "gui-a-half header must not appear here");
}

#[test]
fn gate_digest_matches_proof_bundle() {
    let gate = include_str!("../evidence/gate-output.txt");
    assert!(gate.contains("GATE_DONE"), "gate run did not complete");
    assert!(gate.contains("members ==\n69") || gate.contains("\n69\n"), "member count drift: {gate}");
    for token in ["11329", "0f0b18b7226a813b0922a31f2d28e6941978659e43cbe82eede8be345a088b1d"] {
        assert!(gate.contains(token), "gate digest drift, missing {token}");
    }
}
