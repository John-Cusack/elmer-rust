//! Matrix text format and output formatting shared by the `umf_harness`
//! example and the tests (doubles printed as raw hex bits, exactly as the
//! pristine-library harness prints them).

#![allow(dead_code)]

use elmer_umfpack::amd;
use elmer_umfpack::umf::{self, *};

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
        _ => run_umf(mode, m, &mut o),
    }
    o.0
}

/// Info entries that hold timings (not compared).
pub const TIME_SKIP: [usize; 6] = [
    UMFPACK_SYMBOLIC_TIME,
    UMFPACK_SYMBOLIC_WALLTIME,
    UMFPACK_NUMERIC_TIME,
    UMFPACK_NUMERIC_WALLTIME,
    UMFPACK_SOLVE_TIME,
    UMFPACK_SOLVE_WALLTIME,
];

fn run_umf(mode: &str, m: &Matrix, o: &mut Out) {
    let mut control = [0.0; UMFPACK_CONTROL];
    umf::umfpack_defaults(&mut control);
    let mut info = [0.0; UMFPACK_INFO];
    let res = umf::umfpack_symbolic(
        m.n,
        m.n,
        &m.ap,
        &m.ai,
        Some(&m.ax),
        Some(&control),
        Some(&mut info),
    );
    let st = match &res {
        Ok(_) => 0,
        Err(e) => *e,
    };
    o.line(&format!("sym_status {st}"));
    o.info("sym_info", &info, &TIME_SKIP);
    let Ok(s) = res else { return };
    let (nrow, ncol, nfr, nch) = (
        s.n_row as usize,
        s.n_col as usize,
        s.nfr as usize,
        s.nchains as usize,
    );
    o.line(&format!(
        "sym_dims {} {} {} {} {} {}",
        s.n_row, s.n_col, s.n1, s.nz, s.nfr, s.nchains
    ));
    let p: Vec<i32> = match &s.diagonal_map {
        Some(dm) => (0..nrow).map(|k| s.rperm_init[dm[k] as usize]).collect(),
        None => s.rperm_init[..nrow].to_vec(),
    };
    o.ints("sym_P", &p);
    o.ints("sym_Q", &s.cperm_init[..ncol]);
    o.ints("sym_Front_npivcol", &s.front_npivcol[..=nfr]);
    o.ints("sym_Front_parent", &s.front_parent[..=nfr]);
    o.ints("sym_Front_1strow", &s.front_1strow[..=nfr]);
    o.ints("sym_Front_leftmostdesc", &s.front_leftmostdesc[..=nfr]);
    o.ints("sym_Chain_start", &s.chain_start[..=nch]);
    let mut cm = s.chain_maxrows[..nch].to_vec();
    cm.push(0);
    o.ints("sym_Chain_maxrows", &cm);
    let mut cm = s.chain_maxcols[..nch].to_vec();
    cm.push(0);
    o.ints("sym_Chain_maxcols", &cm);
    // amd_lunz is only set (and only read) under the AMD ordering; the C
    // struct leaves it uninitialised otherwise.
    let skip_lunz: &[usize] = if s.ordering == UMFPACK_ORDERING_AMD {
        &[]
    } else {
        &[5]
    };
    o.info(
        "symx_d",
        &[
            s.num_mem_usage_est,
            s.num_mem_size_est,
            s.peak_sym_usage,
            s.sym,
            s.dnum_mem_init_usage,
            s.amd_lunz,
            s.lunz_bound,
        ],
        skip_lunz,
    );
    let b = |x: bool| x as i32;
    o.ints(
        "symx_i",
        &[
            s.maxnrows,
            s.maxncols,
            s.dense_row_threshold,
            s.n1,
            s.nempty,
            s.esize,
            s.nfr,
            s.n_row,
            s.n_col,
            s.nz,
            s.nb,
            s.num_mem_init_usage,
            s.nempty_row,
            s.nempty_col,
            s.strategy,
            s.ordering,
            b(s.fix_q),
            b(s.prefer_diagonal),
            s.nzaat,
            s.nzdiag,
            s.amd_dmax,
            s.nchains,
            b(s.diagonal_map.is_some()),
        ],
    );
    o.ints("symx_Cperm_init", &s.cperm_init[..=ncol]);
    o.ints("symx_Rperm_init", &s.rperm_init[..=nrow]);
    o.ints("symx_Cdeg", &s.cdeg[..ncol]);
    o.ints("symx_Rdeg", &s.rdeg[..nrow]);
    if let Some(e) = &s.esize_vec {
        o.ints("symx_Esize", e);
    }
    if let Some(dm) = &s.diagonal_map {
        o.ints("symx_Diagonal_map", &dm[..ncol]);
    }
    if mode == "sym" {
        return;
    }
    let res = umf::umfpack_numeric(&m.ap, &m.ai, &m.ax, &s, Some(&control), Some(&mut info));
    let st = match &res {
        Ok((_, st)) => *st,
        Err(e) => *e,
    };
    o.line(&format!("num_status {st}"));
    // NUMERIC_COSTLY_REALLOC counts reallocs that moved the block: allocator-dependent
    let mut skip = TIME_SKIP.to_vec();
    skip.push(UMFPACK_NUMERIC_COSTLY_REALLOC);
    o.info("num_info", &info, &skip);
    let Ok((mut num, _)) = res else { return };
    if mode == "num" {
        return;
    }
    let n = m.n as usize;
    let b: Vec<f64> = match &m.b {
        Some(b) => b.clone(),
        None => (0..n).map(|i| 1.0 + (i % 7) as f64).collect(),
    };
    let mut x = vec![0.0; n];
    // umf4sol sets Control [UMFPACK_IRSTEP] = 0; umfpack_solve here has no
    // iterative refinement
    let st = umf::umfpack_solve(UMFPACK_Aat, &mut x, &b, &mut num, Some(&mut info));
    o.line(&format!("sol_status {st}"));
    o.info("sol_info", &info, &skip);
    o.dbls("X", &x);
}
