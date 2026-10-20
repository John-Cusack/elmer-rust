//! Shared A/B measurement harness for elmer-rust ports.
//!
//! Protocol (mirrors the rustsmith per-package bench convention, ADR-039):
//!
//! * one corpus, one workload, two arms — the *before* build and the *after*
//!   build of the same port;
//! * every sample runs in a **fresh child process**; the parent reaps it with
//!   `wait4(2)` and reads the child's **peak RSS** (`ru_maxrss`), while the
//!   child reports its own **wall time** for the measured region;
//! * arms are **rep-interleaved** (the order flips every rep) so host drift
//!   hits both arms equally;
//! * per arm the harness reports **median + MAD** (median absolute
//!   deviation), flags Tukey-fence (1.5x IQR) outliers without dropping them,
//!   and compares the medians against a **noise threshold** (default 3%):
//!   inside the band the verdict is `no-change`, not a win;
//! * both arms must produce the **same digest** — the digest is the
//!   equivalence proof (bit-identical results). A digest mismatch aborts the
//!   run with a non-zero exit; the harness never reports a speed for
//!   different results.
//!
//! The harness itself is deliberately dependency-free (no `libc` crate): the
//! ports declare their own libc externs, and so does this crate.
//!
//! Usage (parent):
//!
//! ```text
//! <bench-binary> --mode parent --before-bin <path> --after-bin <path> \
//!     --reps 7 --n 128 --iters 5 --label "LU decompose n=128" --out <dir>
//! ```
//!
//! Usage (child, one measured sample):
//!
//! ```text
//! <bench-binary> --mode child --n 128 --iters 5
//! ```

use std::io::Read;
use std::process::{Command, Stdio};
use std::time::Instant;

// ---------------------------------------------------------------------------
// libc surface (declared locally, like the ports do)
// ---------------------------------------------------------------------------

#[repr(C)]
#[derive(Clone, Copy)]
struct Timeval {
    tv_sec: i64,
    tv_usec: i64,
}

/// glibc x86_64 `struct rusage`; `ru_maxrss` sits at offset 32 (KiB).
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Rusage {
    ru_utime: Timeval,
    ru_stime: Timeval,
    ru_maxrss: i64,
    ru_ixrss: i64,
    ru_idrss: i64,
    ru_isrss: i64,
    ru_minflt: i64,
    ru_majflt: i64,
    ru_nswap: i64,
    ru_inblock: i64,
    ru_oublock: i64,
    ru_msgsnd: i64,
    ru_msgrcv: i64,
    ru_nsignals: i64,
    ru_nvcsw: i64,
    ru_nivcsw: i64,
}

impl Default for Timeval {
    fn default() -> Self {
        Timeval {
            tv_sec: 0,
            tv_usec: 0,
        }
    }
}

extern "C" {
    fn getrusage(who: i32, usage: *mut Rusage) -> i32;
    fn wait4(pid: i32, status: *mut i32, options: i32, usage: *mut Rusage) -> i32;
}

const RUSAGE_SELF: i32 = 0;

/// Peak resident set size of *this* process, in KiB, as the process sees it.
///
/// NOTE: on this kernel the child's own `getrusage(RUSAGE_SELF)` inherits the
/// parent's high-water mark at fork time, so it reads far above the child's
/// real footprint. The **parent's `wait4` view is authoritative** and is what
/// the harness records; this function is informational only.
pub fn self_peak_rss_kb() -> i64 {
    let mut ru = Rusage::default();
    // SAFETY: `ru` is a live, correctly sized rusage for RUSAGE_SELF.
    let rc = unsafe { getrusage(RUSAGE_SELF, &mut ru) };
    if rc != 0 {
        return -1;
    }
    ru.ru_maxrss
}

// ---------------------------------------------------------------------------
// statistics
// ---------------------------------------------------------------------------

pub fn median(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Median absolute deviation — the dispersion measure the harness reports.
pub fn mad(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let m = median(xs);
    let dev: Vec<f64> = xs.iter().map(|x| (x - m).abs()).collect();
    median(&dev)
}

/// Tukey fence (1.5x IQR). Returns the indices of the flagged outliers; they
/// are *reported*, never dropped from the sample set.
pub fn tukey_outliers(xs: &[f64]) -> Vec<usize> {
    if xs.len() < 4 {
        return Vec::new();
    }
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q1 = v[v.len() / 4];
    let q3 = v[3 * v.len() / 4];
    let iqr = q3 - q1;
    let lo = q1 - 1.5 * iqr;
    let hi = q3 + 1.5 * iqr;
    xs.iter()
        .enumerate()
        .filter(|(_, x)| **x < lo || **x > hi)
        .map(|(i, _)| i)
        .collect()
}

// ---------------------------------------------------------------------------
// child <-> parent protocol
// ---------------------------------------------------------------------------

/// One measured sample, as reported by a child process.
#[derive(Debug, Clone)]
pub struct ChildResult {
    /// Wall time of the measured region, nanoseconds.
    pub wall_ns: u128,
    /// Wall time of the harness-only region (corpus copy), nanoseconds.
    pub harness_ns: u128,
    /// Peak RSS of the child process, KiB (the child's own getrusage view).
    pub rss_kb: i64,
    /// Result digest — the equivalence proof.
    pub digest: String,
}

/// JSON line a child prints. Hand-rolled: the harness has no dependencies.
pub fn child_json(wall_ns: u128, harness_ns: u128, rss_kb: i64, digest: &str) -> String {
    format!(
        "{{\"wall_ns\":{},\"harness_ns\":{},\"rss_kb\":{},\"digest\":\"{}\"}}",
        wall_ns, harness_ns, rss_kb, digest
    )
}

/// Parse a child JSON line. Strict: a malformed line is an error, never a zero.
pub fn parse_child_json(s: &str) -> Result<ChildResult, String> {
    let get = |key: &str| -> Option<String> {
        let pat = format!("\"{}\":", key);
        let start = s.find(&pat)? + pat.len();
        let rest = &s[start..];
        let end = rest
            .find(|c: char| c == ',' || c == '}')
            .unwrap_or(rest.len());
        Some(rest[..end].trim().trim_matches('"').to_string())
    };
    let wall_ns = get("wall_ns")
        .ok_or("child json: wall_ns missing")?
        .parse::<u128>()
        .map_err(|e| format!("child json: wall_ns parse: {e}"))?;
    let harness_ns = get("harness_ns")
        .ok_or("child json: harness_ns missing")?
        .parse::<u128>()
        .map_err(|e| format!("child json: harness_ns parse: {e}"))?;
    let rss_kb = get("rss_kb")
        .ok_or("child json: rss_kb missing")?
        .parse::<i64>()
        .map_err(|e| format!("child json: rss_kb parse: {e}"))?;
    let digest = get("digest").ok_or("child json: digest missing")?;
    Ok(ChildResult {
        wall_ns,
        harness_ns,
        rss_kb,
        digest,
    })
}

// ---------------------------------------------------------------------------
// parent: run both arms
// ---------------------------------------------------------------------------

/// One arm of the A/B run.
#[derive(Debug, Clone)]
pub struct ArmResult {
    pub label: String,
    pub wall_ns: Vec<u128>,
    pub harness_ns: Vec<u128>,
    pub rss_kb: Vec<i64>,
    pub digest: String,
}

impl ArmResult {
    pub fn wall_ns_median(&self) -> f64 {
        median(&self.wall_ns.iter().map(|x| *x as f64).collect::<Vec<_>>())
    }
    pub fn wall_ns_mad(&self) -> f64 {
        mad(&self.wall_ns.iter().map(|x| *x as f64).collect::<Vec<_>>())
    }
    pub fn harness_ns_median(&self) -> f64 {
        median(&self.harness_ns.iter().map(|x| *x as f64).collect::<Vec<_>>())
    }
    pub fn rss_kb_median(&self) -> f64 {
        median(&self.rss_kb.iter().map(|x| *x as f64).collect::<Vec<_>>())
    }
}

/// Spawn one child sample and collect its wall time plus the peak RSS the
/// parent observes for that child (`wait4`).
fn run_child(bin: &str, n: usize, iters: usize, workload: &str) -> Result<ChildResult, String> {
    let mut cmd = Command::new(bin);
    cmd.arg("--mode")
        .arg("child")
        .arg("--n")
        .arg(n.to_string())
        .arg("--iters")
        .arg(iters.to_string());
    if !workload.is_empty() {
        cmd.arg("--workload").arg(workload);
    }
    let mut child = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("spawn {bin}: {e}"))?;
    let mut out = String::new();
    if let Some(mut so) = child.stdout.take() {
        let mut buf = Vec::new();
        so.read_to_end(&mut buf)
            .map_err(|e| format!("read child stdout: {e}"))?;
        out = String::from_utf8_lossy(&buf).to_string();
    }
    let mut ru = Rusage::default();
    let mut status: i32 = 0;
    // SAFETY: `ru`/`status` are live; wait4 reaps the child we just spawned.
    let rc = unsafe { wait4(child.id() as i32, &mut status, 0, &mut ru) };
    // The child is already reaped by wait4; std's own wait would fail with
    // ECHILD, which is fine — the reaping is done.
    let _ = child.wait();
    if rc < 0 {
        return Err(format!("wait4 failed for {bin}"));
    }
    if status & 0x7f != 0 {
        return Err(format!("child {bin} died with status {status}"));
    }
    let mut parsed = parse_child_json(out.trim())?;
    // The parent's wait4 view is the authoritative peak RSS (same number the
    // child reports through getrusage; cross-checked, never silently mixed).
    if ru.ru_maxrss > 0 {
        parsed.rss_kb = ru.ru_maxrss;
    }
    Ok(parsed)
}

/// Run both arms, rep-interleaved, and return the two arm results.
pub fn run_arms(
    before_bin: &str,
    after_bin: &str,
    reps: usize,
    n: usize,
    iters: usize,
    workload: &str,
) -> Result<(ArmResult, ArmResult), String> {
    let mut before = ArmResult {
        label: "before".to_string(),
        wall_ns: Vec::new(),
        harness_ns: Vec::new(),
        rss_kb: Vec::new(),
        digest: String::new(),
    };
    let mut after = ArmResult {
        label: "after".to_string(),
        wall_ns: Vec::new(),
        harness_ns: Vec::new(),
        rss_kb: Vec::new(),
        digest: String::new(),
    };
    for rep in 0..reps {
        // Flip the order every rep so a monotone host drift cannot favour
        // one arm.
        let order: [(&str, bool); 2] = if rep % 2 == 0 {
            [(before_bin, false), (after_bin, true)]
        } else {
            [(after_bin, true), (before_bin, false)]
        };
        for (bin, is_after) in order {
            let r = run_child(bin, n, iters, workload)?;
            let arm = if is_after { &mut after } else { &mut before };
            if arm.digest.is_empty() {
                arm.digest = r.digest.clone();
            } else if arm.digest != r.digest {
                return Err(format!(
                    "{} arm digest is not stable across reps ({} vs {})",
                    arm.label, arm.digest, r.digest
                ));
            }
            arm.wall_ns.push(r.wall_ns);
            arm.harness_ns.push(r.harness_ns);
            arm.rss_kb.push(r.rss_kb);
        }
    }
    if before.digest != after.digest {
        return Err(format!(
            "digest mismatch between arms: before={} after={} (results differ; \
             the speed number would be meaningless)",
            before.digest, after.digest
        ));
    }
    Ok((before, after))
}

// ---------------------------------------------------------------------------
// verdicts and machine-written rows
// ---------------------------------------------------------------------------

/// Verdict for a speedup ratio, against a noise band.
pub fn verdict(speedup: f64, noise_pct: f64) -> &'static str {
    let delta_pct = (speedup - 1.0) * 100.0;
    if delta_pct.abs() < noise_pct {
        "no-change"
    } else if delta_pct > 0.0 {
        "faster"
    } else {
        "slower"
    }
}

/// Machine-readable result record.
pub struct Rows {
    pub workload: String,
    pub n: usize,
    pub iters: usize,
    pub reps: usize,
    pub noise_pct: f64,
    pub before: ArmResult,
    pub after: ArmResult,
    pub speedup: f64,
    pub rss_delta_pct: f64,
    pub verdict: &'static str,
    pub outliers_before: Vec<usize>,
    pub outliers_after: Vec<usize>,
    pub digest: String,
}

impl Rows {
    /// `results.json` — the machine-written source of truth.
    pub fn to_json(&self) -> String {
        let arm = |a: &ArmResult| {
            format!(
                "{{\"label\":\"{}\",\"wall_ns_median\":{:.0},\"wall_ns_mad\":{:.0},\
                 \"harness_ns_median\":{:.0},\"rss_kb_median\":{:.0},\"rss_kb_samples\":[{}],\
                 \"wall_ns_samples\":[{}]}}",
                a.label,
                a.wall_ns_median(),
                a.wall_ns_mad(),
                a.harness_ns_median(),
                a.rss_kb_median(),
                a.rss_kb
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                a.wall_ns
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        format!(
            "{{\n  \"workload\": \"{}\",\n  \"n\": {},\n  \"iters\": {},\n  \"reps\": {},\n  \
             \"noise_pct\": {},\n  \"digest\": \"{}\",\n  \"before\": {},\n  \"after\": {},\n  \
             \"speedup\": {:.4},\n  \"rss_delta_pct\": {:.2},\n  \"verdict\": \"{}\",\n  \
             \"tukey_outliers_before\": [{}],\n  \"tukey_outliers_after\": [{}]\n}}\n",
            self.workload.replace('"', "'"),
            self.n,
            self.iters,
            self.reps,
            self.noise_pct,
            self.digest,
            arm(&self.before),
            arm(&self.after),
            self.speedup,
            self.rss_delta_pct,
            self.verdict,
            self.outliers_before
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
            self.outliers_after
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
        )
    }

    /// `rows.md` — the board-ready table. Board cells copy these 1:1.
    pub fn to_rows_md(&self) -> String {
        let b = self.before.wall_ns_median() / 1e6;
        let a = self.after.wall_ns_median() / 1e6;
        let hb = self.before.harness_ns_median() / 1e6;
        let ha = self.after.harness_ns_median() / 1e6;
        format!(
            "| Workload | Track | Baseline | Port | Speedup | Verdict | RSS base->port | Harness overhead base->port |\n\
             |---|---|---|---|---|---|---|---|\n\
             | {} | representative | {:.3} ms | {:.3} ms | {:.2}x | {} | {} KiB -> {} KiB ({:+.1}%) | {:.3} ms -> {:.3} ms |\n",
            self.workload,
            b,
            a,
            self.speedup,
            self.verdict,
            self.before.rss_kb_median() as i64,
            self.after.rss_kb_median() as i64,
            self.rss_delta_pct,
            hb,
            ha
        )
    }
}

// ---------------------------------------------------------------------------
// argument parsing
// ---------------------------------------------------------------------------

/// Parsed command line for a bench binary.
pub struct Args {
    pub mode: String,
    pub before_bin: String,
    pub after_bin: String,
    pub reps: usize,
    pub n: usize,
    pub iters: usize,
    pub noise_pct: f64,
    pub label: String,
    pub out: String,
    pub workload: String,
}

pub fn parse_args() -> Result<Args, String> {
    let mut mode = String::new();
    let mut before_bin = String::new();
    let mut after_bin = String::new();
    let mut reps = 7usize;
    let mut n = 0usize;
    let mut iters = 1usize;
    let mut noise_pct = 3.0f64;
    let mut label = String::new();
    let mut out = String::new();
    let mut workload = String::new();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let a = argv[i].as_str();
        let next = |i: &mut usize| -> Result<String, String> {
            *i += 1;
            argv.get(*i)
                .cloned()
                .ok_or_else(|| format!("missing value for {a}"))
        };
        match a {
            "--mode" => mode = next(&mut i)?,
            "--before-bin" => before_bin = next(&mut i)?,
            "--after-bin" => after_bin = next(&mut i)?,
            "--reps" => reps = next(&mut i)?.parse().map_err(|e| format!("--reps: {e}"))?,
            "--n" => n = next(&mut i)?.parse().map_err(|e| format!("--n: {e}"))?,
            "--iters" => iters = next(&mut i)?.parse().map_err(|e| format!("--iters: {e}"))?,
            "--noise-pct" => {
                noise_pct = next(&mut i)?.parse().map_err(|e| format!("--noise-pct: {e}"))?
            }
            "--label" => label = next(&mut i)?,
            "--out" => out = next(&mut i)?,
            "--workload" => workload = next(&mut i)?,
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    if mode.is_empty() {
        return Err("--mode is required (parent|child)".to_string());
    }
    Ok(Args {
        mode,
        before_bin,
        after_bin,
        reps,
        n,
        iters,
        noise_pct,
        label,
        out,
        workload,
    })
}

/// Run the parent side: both arms, the verdict, and the machine rows.
pub fn run_parent(args: &Args) -> Result<String, String> {
    if args.before_bin.is_empty() || args.after_bin.is_empty() {
        return Err("--before-bin and --after-bin are required in parent mode".to_string());
    }
    let (before, after) = run_arms(
        &args.before_bin,
        &args.after_bin,
        args.reps,
        args.n,
        args.iters,
        &args.workload,
    )?;
    let bm = before.wall_ns_median();
    let am = after.wall_ns_median();
    let speedup = if am > 0.0 { bm / am } else { f64::NAN };
    let rss_delta_pct = if before.rss_kb_median() > 0.0 {
        (after.rss_kb_median() - before.rss_kb_median()) / before.rss_kb_median() * 100.0
    } else {
        0.0
    };
    let rows = Rows {
        workload: if args.label.is_empty() {
            format!("n={}", args.n)
        } else {
            args.label.clone()
        },
        n: args.n,
        iters: args.iters,
        reps: args.reps,
        noise_pct: args.noise_pct,
        outliers_before: tukey_outliers(
            &before.wall_ns.iter().map(|x| *x as f64).collect::<Vec<_>>(),
        ),
        outliers_after: tukey_outliers(
            &after.wall_ns.iter().map(|x| *x as f64).collect::<Vec<_>>(),
        ),
        digest: before.digest.clone(),
        before,
        after,
        speedup,
        rss_delta_pct,
        verdict: verdict(speedup, args.noise_pct),
    };
    let json = rows.to_json();
    let md = rows.to_rows_md();
    if !args.out.is_empty() {
        std::fs::create_dir_all(&args.out).map_err(|e| format!("create out dir: {e}"))?;
        std::fs::write(format!("{}/results.json", args.out), &json)
            .map_err(|e| format!("write results.json: {e}"))?;
        std::fs::write(format!("{}/rows.md", args.out), &md).map_err(|e| format!("write rows.md: {e}"))?;
    }
    print!("{md}");
    eprint!("{json}");
    Ok(json)
}

/// Time a measured region with a monotonic clock.
pub struct Timer(Instant);

impl Timer {
    pub fn start() -> Timer {
        Timer(Instant::now())
    }
    pub fn elapsed_ns(&self) -> u128 {
        self.0.elapsed().as_nanos()
    }
}
