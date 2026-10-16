//! Registry consistency: `scope/units.json` must carry exactly the 124
//! lane-gui-a units, each skipped with a known reason, matching the counts
//! published in `elmer_gui_a::scope`.

use std::collections::{HashMap, HashSet};

/// Minimal strict JSON value (subset used by scope/units.json).
#[derive(Debug)]
enum J {
    Str(String),
    Arr(Vec<J>),
    Obj(Vec<(String, J)>),
}

struct P<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> P<'a> {
    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }
    fn lit(&mut self, s: &str) -> bool {
        if self.b[self.i..].starts_with(s.as_bytes()) {
            self.i += s.len();
            true
        } else {
            false
        }
    }
    fn string(&mut self) -> String {
        assert_eq!(self.b[self.i], b'"', "expected string");
        self.i += 1;
        let mut out = String::new();
        loop {
            let c = self.b[self.i];
            self.i += 1;
            match c {
                b'"' => break,
                b'\\' => {
                    let e = self.b[self.i];
                    self.i += 1;
                    out.push(match e {
                        b'n' => '\n',
                        b't' => '\t',
                        b'"' => '"',
                        b'\\' => '\\',
                        b'/' => '/',
                        _ => panic!("bad escape"),
                    });
                }
                _ => out.push(c as char),
            }
        }
        out
    }
    fn value(&mut self) -> J {
        self.ws();
        match self.b[self.i] {
            b'"' => J::Str(self.string()),
            b'[' => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.b[self.i] == b']' {
                    self.i += 1;
                    return J::Arr(v);
                }
                loop {
                    v.push(self.value());
                    self.ws();
                    match self.b[self.i] {
                        b',' => {
                            self.i += 1;
                        }
                        b']' => {
                            self.i += 1;
                            break;
                        }
                        _ => panic!("bad array"),
                    }
                }
                J::Arr(v)
            }
            b'{' => {
                self.i += 1;
                let mut m = Vec::new();
                self.ws();
                if self.b[self.i] == b'}' {
                    self.i += 1;
                    return J::Obj(m);
                }
                loop {
                    self.ws();
                    let k = self.string();
                    self.ws();
                    assert!(self.lit(":"), "expected colon");
                    let val = self.value();
                    m.push((k, val));
                    self.ws();
                    match self.b[self.i] {
                        b',' => {
                            self.i += 1;
                        }
                        b'}' => {
                            self.i += 1;
                            break;
                        }
                        _ => panic!("bad object"),
                    }
                }
                J::Obj(m)
            }
            _ => panic!("unexpected json"),
        }
    }
}

fn field<'a>(obj: &'a [(String, J)], k: &str) -> &'a J {
    &obj.iter().find(|(kk, _)| kk == k).unwrap_or_else(|| panic!("missing {k}")).1
}

fn as_str(j: &J) -> &str {
    match j {
        J::Str(s) => s,
        _ => panic!("expected str"),
    }
}

#[test]
fn registry_matches_lane_partition() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/scope/units.json"
    ))
    .unwrap();
    let mut p = P { b: text.as_bytes(), i: 0 };
    let root = p.value();
    let J::Obj(top) = &root else { panic!("root obj") };
    assert_eq!(as_str(field(top, "lane")), elmer_gui_a::scope::LANE);
    let J::Arr(units) = field(top, "units") else { panic!("units arr") };
    assert_eq!(units.len(), elmer_gui_a::scope::UNIT_COUNT);

    let mut seen = HashSet::new();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for u in units {
        let J::Obj(f) = u else { panic!("unit obj") };
        let id = as_str(field(f, "unit")).to_string();
        let lang = as_str(field(f, "lang"));
        assert!(
            id.starts_with("c:ElmerGUI/netgen/") || id.starts_with("cxx:ElmerGUI/netgen/"),
            "lane path violation: {id}"
        );
        assert!(lang == "c" || lang == "cxx", "bad lang {lang}");
        assert_eq!(as_str(field(f, "verdict")), "skipped");
        let reason = as_str(field(f, "reason")).to_string();
        if !matches!(
            reason.as_str(),
            "cxx-mangled-only"
                | "cxx-header-only"
                | "c-header-decl-only"
                | "dead-upstream-uncompilable"
        ) {
            panic!("unknown reason {reason}");
        }
        *counts.entry(reason).or_insert(0) += 1;
        assert!(seen.insert(id), "duplicate unit id");
    }
    for (reason, want) in elmer_gui_a::scope::REASON_COUNTS {
        assert_eq!(
            counts.get(*reason).copied().unwrap_or(0),
            *want,
            "reason {reason}"
        );
    }
    assert_eq!(counts.values().sum::<usize>(), units.len());
}
