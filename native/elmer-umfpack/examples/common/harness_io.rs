//! Matrix text format and output formatting shared by the `umf_harness`
//! example and the tests (doubles printed as raw hex bits, exactly as the
//! pristine-library harness prints them).

#![allow(dead_code)]

use elmer_umfpack::amd;

pub struct Matrix {
    pub n: i32,
    pub ap: Vec<i32>,
    pub ai: Vec<i32>,
    pub ax: Vec<f64>,
    pub b: Option<Vec<f64>>,
}

pub fn parse(text: &str) -> Matrix {
    let mut tok = text.split_ascii_whitespace();
    let n: i32 = tok.next().expect("n").parse().expect("n");
    let nz: usize = tok.next().expect("nz").parse().expect("nz");
    let ap: Vec<i32> = (0..=n)
        .map(|_| tok.next().unwrap().parse().unwrap())
        .collect();
    let ai: Vec<i32> = (0..nz)
        .map(|_| tok.next().unwrap().parse().unwrap())
        .collect();
    let bits = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).expect("hex"));
    let ax: Vec<f64> = (0..nz).map(|_| bits(tok.next().unwrap())).collect();
    let b = match tok.next() {
        Some(t) if t.starts_with('b') => Some((0..n).map(|_| bits(tok.next().unwrap())).collect()),
        _ => None,
    };
    Matrix { n, ap, ai, ax, b }
}

#[derive(Default)]
pub struct Out(pub String);

impl Out {
    pub fn line(&mut self, s: &str) {
        self.0.push_str(s);
        self.0.push('\n');
    }
    pub fn ints(&mut self, name: &str, a: &[i32]) {
        self.0.push_str(name);
        for v in a {
            self.0.push_str(&format!(" {v}"));
        }
        self.0.push('\n');
    }
    pub fn dbls(&mut self, name: &str, a: &[f64]) {
        self.0.push_str(name);
        for v in a {
            self.0.push_str(&format!(" {:016x}", v.to_bits()));
        }
        self.0.push('\n');
    }
    pub fn info(&mut self, name: &str, info: &[f64], skip: &[usize]) {
        self.0.push_str(name);
        for (i, v) in info.iter().enumerate() {
            if skip.contains(&i) {
                self.0.push_str(" -");
            } else {
                self.0.push_str(&format!(" {:016x}", v.to_bits()));
            }
        }
        self.0.push('\n');
    }
}

/// Run one harness mode on a matrix; the text the pristine harness prints.
pub fn run(mode: &str, m: &Matrix) -> String {
    let mut o = Out::default();
    match mode {
        "amd" => {
            let mut info = [0.0; amd::AMD_INFO];
            let mut p = vec![0i32; m.n as usize + 1];
            let st = amd::amd_order(m.n, &m.ap, &m.ai, &mut p, None, Some(&mut info));
            o.line(&format!("status {st}"));
            o.ints("P", &p[..m.n as usize]);
            o.info("info", &info, &[]);
        }
        other => panic!("mode {other} not ported yet"),
    }
    o.0
}
