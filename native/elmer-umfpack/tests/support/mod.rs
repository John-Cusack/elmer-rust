//! Expected-output files: `== name` headers, then the pristine harness lines.

#![allow(dead_code)]

use std::collections::BTreeMap;

pub fn data(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

pub fn expected(name: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let mut cur: Option<(String, String)> = None;
    for line in data(name).lines() {
        if let Some(n) = line.strip_prefix("== ") {
            if let Some((k, v)) = cur.take() {
                map.insert(k, v);
            }
            cur = Some((n.to_string(), String::new()));
        } else if let Some((_, v)) = cur.as_mut() {
            v.push_str(line);
            v.push('\n');
        }
    }
    if let Some((k, v)) = cur {
        map.insert(k, v);
    }
    map
}

/// Compare line by line and report the first differing line.
pub fn check(expected: &BTreeMap<String, String>, name: &str, got: &str) {
    let want = expected
        .get(name)
        .unwrap_or_else(|| panic!("no expected output for {name}"));
    for (w, g) in want.lines().zip(got.lines()) {
        if w != g {
            let tag = w.split(' ').next().unwrap_or("");
            let (wv, gv): (Vec<&str>, Vec<&str>) = (w.split(' ').collect(), g.split(' ').collect());
            let at = wv.iter().zip(&gv).position(|(a, b)| a != b).unwrap_or(0);
            panic!(
                "{name}: line `{tag}` differs at field {at}: want {:?} got {:?}",
                wv.get(at),
                gv.get(at)
            );
        }
    }
    assert_eq!(
        want.lines().count(),
        got.lines().count(),
        "{name}: line count"
    );
}

/// FNV-1a 64 of one output line (the hashed expected files store
/// `<tag> <hash>` per line, generated from the pristine harness output).
pub fn fnv1a(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Compare against a hashed expected file; report the first differing line.
pub fn check_hashed(expected: &BTreeMap<String, String>, name: &str, got: &str) {
    let want = expected
        .get(name)
        .unwrap_or_else(|| panic!("no expected output for {name}"));
    let mut n = 0;
    for (w, g) in want.lines().zip(got.lines()) {
        let tag = g.split(' ').next().unwrap_or("");
        let line = format!("{tag} {:016x}", fnv1a(g));
        assert_eq!(
            w, line,
            "{name}: line `{tag}` differs from the pristine library (got: {:.200})",
            g
        );
        n += 1;
    }
    assert_eq!(want.lines().count(), n, "{name}: line count");
    assert_eq!(got.lines().count(), n, "{name}: line count");
}
