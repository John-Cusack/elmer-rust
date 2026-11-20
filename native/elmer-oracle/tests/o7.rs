//! O7 kernel fixtures (report §5.2 item 8).
//!
//! The first test always runs: the committed fixtures parse and cover
//! every element type S1 needs. The second needs a pristine arm and its
//! shim (`scripts/o7-build-shim.sh`), named by `ELMER_O7_BUILD` and
//! `ELMER_O7_SHIM`. It re-evaluates every fixture input in that arm and
//! requires 0-ulp agreement. Without the variables it reports itself
//! skipped and passes.

use elmer_oracle::o7::{diff_info, diff_rules, parse_info, parse_rules, replay_info, replay_rules, Shim};
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("o7").join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

/// S1's element types (report §5.3): linear and quadratic 202/303/404/504/
/// 808 and 306/408/510/820, plus the quadratic line 203 on their edges.
const S1_TYPES: [u16; 10] = [202, 203, 303, 306, 404, 408, 504, 510, 808, 820];

#[test]
fn fixtures_cover_the_s1_element_types() {
    let rules = parse_rules(&fixture("rules.txt")).unwrap();
    let info = parse_info(&fixture("info.txt")).unwrap();
    let types: std::collections::BTreeSet<u16> = info.iter().map(|e| e.code).collect();
    for t in S1_TYPES {
        assert!(types.contains(&t), "no ElementInfo fixture for type {t}");
        let ip = info.iter().filter(|e| e.code == t).flat_map(|e| &e.pts).filter(|p| p.tag == "ip").count();
        assert!(ip > 0, "type {t} has no integration-point call");
    }
    let kinds: std::collections::BTreeSet<&str> = rules.iter().map(|r| r.kind.as_str()).collect();
    for k in ["line", "tri", "quad", "tet", "brick"] {
        assert!(kinds.contains(k), "no {k} rule");
    }
}

#[test]
fn pristine_arm_reproduces_the_fixtures_bitwise() {
    let (Some(build), Some(shim)) = (std::env::var_os("ELMER_O7_BUILD"), std::env::var_os("ELMER_O7_SHIM")) else {
        eprintln!("skipped: set ELMER_O7_BUILD and ELMER_O7_SHIM to run O7 against a pristine arm");
        return;
    };
    let shim = Shim::load(&PathBuf::from(build), &PathBuf::from(shim)).unwrap();
    let rules = parse_rules(&fixture("rules.txt")).unwrap();
    let info = parse_info(&fixture("info.txt")).unwrap();
    // GaussPointsQuad never writes w: it is whatever the previous rule
    // left in Elmer's shared workspace, so it is not compared.
    let dr = diff_rules(&rules, &replay_rules(&shim, &rules).unwrap(), true);
    let di = diff_info(&info, &replay_info(&shim, &info).unwrap());
    assert!(dr.is_empty() && di.is_empty(), "{} rule and {} ElementInfo differences, first: {:?}", dr.len(), di.len(), dr.iter().chain(&di).take(5).collect::<Vec<_>>());
}
