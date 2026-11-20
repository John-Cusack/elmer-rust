//! Expected-output files: `== name` headers, then the pristine harness lines.

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
