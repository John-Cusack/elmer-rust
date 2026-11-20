//! O7: kernel-level differential against the pristine Fortran (report
//! §5.2 item 8, §4 O7).
//!
//! The oracle loads one pristine arm's `libelmersolver` through a small
//! `BIND(C)` shim (`o7/o7_shim.f90`, built by `scripts/o7-build-shim.sh`)
//! and calls the Integration.F90 rule builders and `ElementInfo` directly.
//! Fixtures use the text format of the candidate element library's
//! captured data (reals as f64 bit patterns in hex), so a candidate can
//! test against them at 0 ulp without depending on this crate.
//!
//! Only this module uses `unsafe`: the FFI to `dlopen` and to the shim.
//! Every shim call holds one process-wide lock, because the Fortran
//! kernels share module-level workspaces.

use crate::{OracleError, Result};
use std::ffi::{c_char, c_int, c_void, CString};
use std::path::Path;
use std::sync::Mutex;

#[link(name = "dl")]
extern "C" {
    fn dlopen(filename: *const c_char, flag: c_int) -> *mut c_void;
    fn dlsym(handle: *mut c_void, symbol: *const c_char) -> *mut c_void;
    fn dlerror() -> *const c_char;
}

const RTLD_NOW: c_int = 2;
const RTLD_GLOBAL: c_int = 0x100;

type InitFn = unsafe extern "C" fn(c_int);
type RuleFn = unsafe extern "C" fn(c_int, c_int, c_int, c_int, *mut f64, *mut f64, *mut f64, *mut f64) -> c_int;
type InfoFn = unsafe extern "C" fn(
    c_int,
    c_int,
    c_int,
    *const f64,
    *const f64,
    *const f64,
    f64,
    f64,
    f64,
    *mut f64,
    *mut f64,
    *mut f64,
) -> c_int;

static LOCK: Mutex<()> = Mutex::new(());

/// The loaded shim of one pristine arm.
pub struct Shim {
    init: InitFn,
    rule: RuleFn,
    info: InfoFn,
}

fn open(path: &str) -> Result<*mut c_void> {
    let c = CString::new(path).map_err(|_| OracleError::Parse(format!("NUL in {path}")))?;
    // SAFETY: `c` is a valid NUL-terminated string; dlopen has no other
    // preconditions. The handle is never closed: the arm stays loaded for
    // the life of the process.
    let h = unsafe { dlopen(c.as_ptr(), RTLD_NOW | RTLD_GLOBAL) };
    if h.is_null() {
        // SAFETY: dlerror returns NULL or a NUL-terminated message.
        let msg = unsafe {
            let e = dlerror();
            if e.is_null() {
                String::from("unknown dlopen error")
            } else {
                std::ffi::CStr::from_ptr(e).to_string_lossy().into_owned()
            }
        };
        return Err(OracleError::SolverFailed(format!("dlopen {path}: {msg}")));
    }
    Ok(h)
}

fn sym(h: *mut c_void, name: &str) -> Result<*mut c_void> {
    let c = CString::new(name).map_err(|_| OracleError::Parse(name.to_string()))?;
    // SAFETY: `h` is a live handle from dlopen and `c` a valid C string.
    let p = unsafe { dlsym(h, c.as_ptr()) };
    if p.is_null() {
        return Err(OracleError::SolverFailed(format!("shim has no symbol {name}")));
    }
    Ok(p)
}

impl Shim {
    /// Load `shim` built against `build` (a pristine arm's build tree).
    /// The arm's libraries are opened first by absolute path, because
    /// `libelmersolver`'s RUNPATH names the container tree it was built
    /// in. `ELMER_HOME` must point at `<build>/fem/src` (where
    /// `elements.def` lives) before the first call; [`Shim::load`] sets
    /// it when it is unset.
    pub fn load(build: &Path, shim: &Path) -> Result<Shim> {
        let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
        if std::env::var_os("ELMER_HOME").is_none() {
            std::env::set_var("ELMER_HOME", build.join("fem/src"));
        }
        let b = build.display();
        for lib in [
            String::from("libblas.so.3"),
            String::from("liblapack.so.3"),
            format!("{b}/mathlibs/src/arpack/libarpack.so"),
            format!("{b}/mathlibs/src/parpack/libparpack.so"),
            format!("{b}/matc/src/libmatc.so"),
            format!("{b}/fhutiter/src/libfhuti.so"),
            format!("{b}/fem/src/libelmersolver.so"),
        ] {
            open(&lib)?;
        }
        let h = open(&shim.display().to_string())?;
        // SAFETY: the symbols are the shim's BIND(C) entry points, whose C
        // signatures are exactly these function types (o7/o7_shim.f90).
        unsafe {
            Ok(Shim {
                init: std::mem::transmute::<*mut c_void, InitFn>(sym(h, "o7_init")?),
                rule: std::mem::transmute::<*mut c_void, RuleFn>(sym(h, "o7_rule")?),
                info: std::mem::transmute::<*mut c_void, InfoFn>(sym(h, "o7_info")?),
            })
        }
    }

    /// One Gauss rule from the pristine rule builder for `kind`.
    pub fn rule(&self, kind: &str, args: &[usize]) -> Result<Vec<[f64; 4]>> {
        let k = match kind {
            "line" => 1,
            "tri" => 2,
            "quad" => 3,
            "tet" => 4,
            "pyr" => 5,
            "wedge" => 6,
            "brick" => 7,
            "wedge2" => 8,
            "wecon" => 9,
            _ => return Err(OracleError::Unsupported(format!("rule kind {kind}"))),
        };
        let a = *args.first().ok_or_else(|| OracleError::Parse(format!("rule {kind} without n")))? as c_int;
        let b = args.get(1).copied().unwrap_or(0) as c_int;
        const NMAX: usize = 4096;
        let (mut u, mut v, mut w, mut s) = (vec![0.0; NMAX], vec![0.0; NMAX], vec![0.0; NMAX], vec![0.0; NMAX]);
        let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: the four buffers hold NMAX doubles each, the bound the
        // shim is given; it writes at most that many.
        let n = unsafe {
            (self.init)(3);
            (self.rule)(k, a, b, NMAX as c_int, u.as_mut_ptr(), v.as_mut_ptr(), w.as_mut_ptr(), s.as_mut_ptr())
        };
        if n < 0 {
            return Err(OracleError::SolverFailed(format!("o7_rule {kind} {args:?}: status {n}")));
        }
        Ok((0..n as usize).map(|i| [u[i], v[i], w[i], s[i]]).collect())
    }

    /// `ElementInfo` of type `code` with nodes `xyz` at `uvw` in a model
    /// of dimension `cdim`: (detJ, Basis, dBasisdx per node).
    pub fn info(&self, code: u16, cdim: usize, xyz: &[[f64; 3]], uvw: [f64; 3]) -> Result<InfoOut> {
        let n = xyz.len();
        let x: Vec<f64> = xyz.iter().map(|p| p[0]).collect();
        let y: Vec<f64> = xyz.iter().map(|p| p[1]).collect();
        let z: Vec<f64> = xyz.iter().map(|p| p[2]).collect();
        let mut detj = 0.0;
        let mut basis = vec![0.0; n];
        let mut db = vec![0.0; 3 * n];
        let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: x, y, z and basis hold n doubles, db 3n (Fortran
        // dBasisdx(n,3), column-major), as the shim's interface declares.
        let stat = unsafe {
            (self.info)(
                c_int::from(code),
                cdim as c_int,
                n as c_int,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
                uvw[0],
                uvw[1],
                uvw[2],
                &mut detj,
                basis.as_mut_ptr(),
                db.as_mut_ptr(),
            )
        };
        if stat != 1 {
            return Err(OracleError::SolverFailed(format!("o7_info type {code}: status {stat}")));
        }
        let dbasis = (0..n).map(|i| [db[i], db[n + i], db[2 * n + i]]).collect();
        Ok(InfoOut { detj, basis, dbasis })
    }
}

/// One `ElementInfo` result.
#[derive(Debug, Clone, PartialEq)]
pub struct InfoOut {
    pub detj: f64,
    pub basis: Vec<f64>,
    pub dbasis: Vec<[f64; 3]>,
}

/// One Gauss rule fixture: builder kind, its arguments, and (u, v, w, s).
#[derive(Debug, Clone, PartialEq)]
pub struct Rule {
    pub kind: String,
    pub args: Vec<usize>,
    pub pts: Vec<[f64; 4]>,
}

/// One element fixture: nodes and the `ElementInfo` calls made on it.
#[derive(Debug, Clone, PartialEq)]
pub struct Elem {
    pub case: String,
    pub code: u16,
    pub cdim: usize,
    /// Trailing header fields after `cdim`, kept verbatim.
    pub extra: String,
    pub xyz: Vec<[f64; 3]>,
    pub pts: Vec<Pt>,
}

/// One `ElementInfo` call: where it came from (`ip` or `stab`), its local
/// coordinates and its outputs.
#[derive(Debug, Clone, PartialEq)]
pub struct Pt {
    pub tag: String,
    pub uvw: [f64; 3],
    pub out: InfoOut,
}

fn hex(s: &str) -> Result<f64> {
    u64::from_str_radix(s, 16)
        .map(f64::from_bits)
        .map_err(|_| OracleError::Parse(format!("bad hex real {s}")))
}

fn hexes(line: Option<&str>, tag: &str) -> Result<Vec<f64>> {
    let line = line.ok_or_else(|| OracleError::Parse(format!("missing {tag} line")))?;
    let mut t = line.split_whitespace();
    if t.next() != Some(tag) {
        return Err(OracleError::Parse(format!("expected {tag} line, got {line:?}")));
    }
    t.map(hex).collect()
}

fn h(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

/// Parse a rules fixture (`rule <kind> <args..> <n>` then n lines of u v w s).
pub fn parse_rules(text: &str) -> Result<Vec<Rule>> {
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(head) = lines.next() {
        let t: Vec<&str> = head.split_whitespace().collect();
        if t.len() < 4 || t[0] != "rule" {
            return Err(OracleError::Parse(format!("bad rule header {head:?}")));
        }
        let num = |s: &str| s.parse::<usize>().map_err(|_| OracleError::Parse(format!("bad count in {head:?}")));
        let args = t[2..t.len() - 1].iter().map(|s| num(s)).collect::<Result<Vec<_>>>()?;
        let n = num(t[t.len() - 1])?;
        let mut pts = Vec::with_capacity(n);
        for _ in 0..n {
            let l = lines.next().ok_or_else(|| OracleError::Parse(format!("{head}: truncated")))?;
            let r = l.split_whitespace().map(hex).collect::<Result<Vec<_>>>()?;
            let [u, v, w, s] = r[..] else {
                return Err(OracleError::Parse(format!("{head}: bad point {l:?}")));
            };
            pts.push([u, v, w, s]);
        }
        out.push(Rule { kind: t[1].to_string(), args, pts });
    }
    Ok(out)
}

pub fn write_rules(rules: &[Rule]) -> String {
    let mut s = String::new();
    for r in rules {
        let args: Vec<String> = r.args.iter().map(|a| a.to_string()).collect();
        s.push_str(&format!("rule {} {} {}\n", r.kind, args.join(" "), r.pts.len()));
        for p in &r.pts {
            s.push_str(&format!("{} {} {} {}\n", h(p[0]), h(p[1]), h(p[2]), h(p[3])));
        }
    }
    s
}

/// Parse an element fixture: `elem <case> <code> <n> <cdim>`, `xyz` (all x,
/// all y, all z), then per call `pt <tag> u v w detJ`, `b` (n values) and
/// `d` (per node dx dy dz).
pub fn parse_info(text: &str) -> Result<Vec<Elem>> {
    let mut out: Vec<Elem> = Vec::new();
    let mut lines = text.lines();
    while let Some(l) = lines.next() {
        let t: Vec<&str> = l.split_whitespace().collect();
        let bad = || OracleError::Parse(format!("bad line {l:?}"));
        match t.first().copied() {
            Some("elem") if t.len() >= 5 => {
                let n: usize = t[3].parse().map_err(|_| bad())?;
                let x = hexes(lines.next(), "xyz")?;
                if x.len() != 3 * n {
                    return Err(bad());
                }
                out.push(Elem {
                    case: t[1].to_string(),
                    code: t[2].parse().map_err(|_| bad())?,
                    cdim: t[4].parse().map_err(|_| bad())?,
                    extra: t[5..].join(" "),
                    xyz: (0..n).map(|i| [x[i], x[n + i], x[2 * n + i]]).collect(),
                    pts: Vec::new(),
                });
            }
            Some("pt") if t.len() == 6 => {
                let v = t[2..].iter().map(|s| hex(s)).collect::<Result<Vec<_>>>()?;
                let basis = hexes(lines.next(), "b")?;
                let d = hexes(lines.next(), "d")?;
                let e = out.last_mut().ok_or_else(bad)?;
                if basis.len() != e.xyz.len() || d.len() != 3 * e.xyz.len() {
                    return Err(bad());
                }
                e.pts.push(Pt {
                    tag: t[1].to_string(),
                    uvw: [v[0], v[1], v[2]],
                    out: InfoOut { detj: v[3], basis, dbasis: d.chunks(3).map(|c| [c[0], c[1], c[2]]).collect() },
                });
            }
            _ => return Err(bad()),
        }
    }
    Ok(out)
}

pub fn write_info(elems: &[Elem]) -> String {
    let mut s = String::new();
    for e in elems {
        let extra = if e.extra.is_empty() { String::new() } else { format!(" {}", e.extra) };
        s.push_str(&format!("elem {} {} {} {}{extra}\n", e.case, e.code, e.xyz.len(), e.cdim));
        let coords: Vec<String> = (0..3).flat_map(|c| e.xyz.iter().map(move |p| h(p[c]))).collect();
        s.push_str(&format!("xyz {}\n", coords.join(" ")));
        for p in &e.pts {
            s.push_str(&format!(
                "pt {} {} {} {} {}\n",
                p.tag,
                h(p.uvw[0]),
                h(p.uvw[1]),
                h(p.uvw[2]),
                h(p.out.detj)
            ));
            let b: Vec<String> = p.out.basis.iter().map(|&x| h(x)).collect();
            s.push_str(&format!("b {}\n", b.join(" ")));
            let d: Vec<String> = p.out.dbasis.iter().flat_map(|r| r.iter().map(|&x| h(x))).collect();
            s.push_str(&format!("d {}\n", d.join(" ")));
        }
    }
    s
}

/// The inputs of `rules` evaluated by the pristine arm.
pub fn replay_rules(shim: &Shim, rules: &[Rule]) -> Result<Vec<Rule>> {
    rules
        .iter()
        .map(|r| Ok(Rule { kind: r.kind.clone(), args: r.args.clone(), pts: shim.rule(&r.kind, &r.args)? }))
        .collect()
}

/// The inputs of `elems` (types, nodes, call points) evaluated by the
/// pristine arm, in the same call order.
pub fn replay_info(shim: &Shim, elems: &[Elem]) -> Result<Vec<Elem>> {
    elems
        .iter()
        .map(|e| {
            let pts = e
                .pts
                .iter()
                .map(|p| {
                    Ok(Pt { tag: p.tag.clone(), uvw: p.uvw, out: shim.info(e.code, e.cdim, &e.xyz, p.uvw)? })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(Elem { pts, ..e.clone() })
        })
        .collect()
}

/// Reals that differ bitwise between two fixture sets of the same inputs:
/// (rule differences, info differences), each as a short description.
pub fn diff_rules(a: &[Rule], b: &[Rule], skip_quad_w: bool) -> Vec<String> {
    let mut out = Vec::new();
    for (x, y) in a.iter().zip(b) {
        if x.pts.len() != y.pts.len() {
            out.push(format!("rule {} {:?}: {} vs {} points", x.kind, x.args, x.pts.len(), y.pts.len()));
            continue;
        }
        for (i, (p, q)) in x.pts.iter().zip(&y.pts).enumerate() {
            for c in 0..4 {
                if c == 2 && skip_quad_w && x.kind == "quad" {
                    continue;
                }
                if p[c].to_bits() != q[c].to_bits() {
                    out.push(format!("rule {} {:?} point {i} component {c}: {:e} vs {:e}", x.kind, x.args, p[c], q[c]));
                }
            }
        }
    }
    if a.len() != b.len() {
        out.push(format!("{} vs {} rules", a.len(), b.len()));
    }
    out
}

pub fn diff_info(a: &[Elem], b: &[Elem]) -> Vec<String> {
    let mut out = Vec::new();
    for (x, y) in a.iter().zip(b) {
        for (k, (p, q)) in x.pts.iter().zip(&y.pts).enumerate() {
            let ctx = format!("{} type {} cdim {} pt {k}", x.case, x.code, x.cdim);
            if p.out.detj.to_bits() != q.out.detj.to_bits() {
                out.push(format!("{ctx}: detJ {:e} vs {:e}", p.out.detj, q.out.detj));
            }
            for (i, (u, v)) in p.out.basis.iter().zip(&q.out.basis).enumerate() {
                if u.to_bits() != v.to_bits() {
                    out.push(format!("{ctx}: basis {i} {u:e} vs {v:e}"));
                }
            }
            for (i, (u, v)) in p.out.dbasis.iter().zip(&q.out.dbasis).enumerate() {
                for j in 0..3 {
                    if u[j].to_bits() != v[j].to_bits() {
                        out.push(format!("{ctx}: dBasisdx({i},{j}) {:e} vs {:e}", u[j], v[j]));
                    }
                }
            }
        }
        if x.pts.len() != y.pts.len() {
            out.push(format!("{} type {}: {} vs {} points", x.case, x.code, x.pts.len(), y.pts.len()));
        }
    }
    if a.len() != b.len() {
        out.push(format!("{} vs {} elements", a.len(), b.len()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rules_round_trip() {
        let text = "rule tri 3 3\n3fc5555555555555 3fc5555555555555 0000000000000000 3fc5555555555555\n\
                    3fe5555555555555 3fc5555555555555 0000000000000000 3fc5555555555555\n\
                    3fc5555555555555 3fe5555555555555 0000000000000000 3fc5555555555555\n\
                    rule wedge2 3 2 6\n"
            .to_string()
            + &"0000000000000000 0000000000000000 0000000000000000 3ff0000000000000\n".repeat(6);
        let r = parse_rules(&text).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[1].args, vec![3, 2]);
        assert_eq!(write_rules(&r), text);
        assert!(diff_rules(&r, &r, false).is_empty());
        let mut q = r.clone();
        q[0].pts[1][3] = 0.5;
        assert_eq!(diff_rules(&r, &q, false).len(), 1);
    }

    #[test]
    fn info_round_trip_and_diff() {
        let z = "0000000000000000";
        let one = "3ff0000000000000";
        let text = format!(
            "elem c1d 202 2 1 1\nxyz bff0000000000000 {one} {z} {z} {z} {z}\n\
             pt ip {z} {z} {z} {one}\nb 3fe0000000000000 3fe0000000000000\nd bfe0000000000000 {z} {z} 3fe0000000000000 {z} {z}\n"
        );
        let e = parse_info(&text).unwrap();
        assert_eq!(e[0].xyz, vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(e[0].pts[0].out.dbasis, vec![[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]);
        assert_eq!(write_info(&e), text);
        let mut f = e.clone();
        f[0].pts[0].out.dbasis[1][0] = f64::from_bits(0.5f64.to_bits() + 1);
        assert_eq!(diff_info(&e, &f), vec![String::from(
            "c1d type 202 cdim 1 pt 0: dBasisdx(1,0) 5e-1 vs 5.000000000000001e-1"
        )]);
    }
}
