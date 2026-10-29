//! Corpus proof driver: parse every `.sif` in the 482 quick tests.
//!
//! Usage:
//! - `sifproof -- <tests-dir> <times-file> <out-report> [only-test]`
//! - `sifproof -- --lexdump <file>`: print post-expansion logical lines
//!
//! Elmer runs each test in its own process with a fresh MATC environment,
//! so the driver forks one child per test (`only-test`) and aggregates.
//! Prints one line per failure plus a summary; the full report goes to
//! `<out-report>` (outside the worktree, never committed).

use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(|s| s.as_str()) == Some("--lexdump") {
        dump_logical(&args[2]);
        return;
    }
    let tests_dir = args[1].clone();
    let times_file = args[2].clone();
    let out_path = args[3].clone();
    let only = args.get(4).cloned();

    let times = std::fs::read_to_string(&times_file).unwrap();
    let mut names: Vec<String> = times
        .lines()
        .map(|l| l.split_whitespace().next().unwrap().to_string())
        .collect();
    names.sort();
    names.dedup();

    if let Some(one) = only {
        // Child: parse a single test, print machine-readable result.
        let (ok, fail, missing) = run_one(&tests_dir, &one);
        println!("RESULT ok={ok} fail={fail} missing={missing}");
        return;
    }

    let mut report = String::new();
    let mut ok_files = 0usize;
    let mut fail_files = 0usize;
    let mut ok_tests = 0usize;
    let mut fail_tests = 0usize;
    let mut missing = Vec::new();

    let exe = std::env::current_exe().unwrap();
    for name in &names {
        let out = std::process::Command::new(&exe)
            .arg(&tests_dir)
            .arg(&times_file)
            .arg("/dev/null")
            .arg(name)
            .output()
            .expect("spawn proof child");
        let stdout = String::from_utf8_lossy(&out.stdout);
        let mut ok = 0usize;
        let mut fail = 0usize;
        let mut miss = false;
        for line in stdout.lines() {
            if let Some(r) = line.strip_prefix("RESULT ") {
                for kv in r.split_whitespace() {
                    let (k, v) = kv.split_once('=').unwrap();
                    match k {
                        "ok" => ok = v.parse().unwrap(),
                        "fail" => fail = v.parse().unwrap(),
                        "missing" => miss = v == "1",
                        _ => {}
                    }
                }
            } else if line.starts_with("FAIL ") {
                report.push_str(line);
                report.push('\n');
            }
        }
        if !out.status.success() {
            report.push_str(&format!(
                "FAIL {name} :: child failed: {}\n",
                String::from_utf8_lossy(&out.stderr).lines().next().unwrap_or("")
            ));
            fail_tests += 1;
            continue;
        }
        ok_files += ok;
        fail_files += fail;
        if miss {
            missing.push(name.clone());
        } else if fail == 0 {
            ok_tests += 1;
        } else {
            fail_tests += 1;
        }
    }

    let summary = format!(
        "tests: {} ok={} fail={} missing={}\nfiles: ok={} fail={}\nmissing: {:?}\n",
        names.len(),
        ok_tests,
        fail_tests,
        missing.len(),
        ok_files,
        fail_files,
        missing,
    );
    let mut f = std::fs::File::create(&out_path).unwrap();
    f.write_all(summary.as_bytes()).unwrap();
    f.write_all(report.as_bytes()).unwrap();
    println!("{summary}");
    for line in report.lines().take(40) {
        println!("{line}");
    }
}

/// Print the per-call echo stream, mirroring Elmer's `Echo on` output:
/// only calls from the `echo on` line onward (the flag is set when that
fn dump_logical(path: &str) {
    let lines = elmer_sif::lex::Lexer::echo_stream(path).unwrap();
    let mut echo = false;
    for line in &lines {
        let low = line.to_lowercase();
        if !echo {
            if low.starts_with("echo on") {
                echo = true;
            }
            continue;
        }
        if low.starts_with("echo off") {
            echo = false;
            continue;
        }
        println!("{line}");
    }
}

fn run_one(tests_dir: &str, name: &str) -> (usize, usize, usize) {
    let mut dir = format!("{tests_dir}/{name}");
    if !std::path::Path::new(&dir).is_dir() {
        let base = name.rsplit_once("_np").map(|(b, _)| b).unwrap_or(name);
        dir = format!("{tests_dir}/{base}");
    }
    if !std::path::Path::new(&dir).is_dir() {
        return (0, 0, 1);
    }
    // Elmer resolves relative includes against the process working
    // directory (ctest runs in the test dir); do the same per test.
    std::env::set_current_dir(&dir).unwrap();
    let mut sifs: Vec<String> = Vec::new();
    collect_sifs(&dir, &mut sifs);
    if sifs.is_empty() {
        return (0, 0, 1);
    }
    let mut ok = 0usize;
    let mut fail = 0usize;
    for sif in &sifs {
        match elmer_sif::parse_file(sif) {
            Ok(mut parsed) => {
                for run in &mut parsed.sif.runs {
                    if let Err(e) = elmer_sif::complete::complete_run(run, sif) {
                        fail += 1;
                        println!("FAIL {sif} :: complete: {e}");
                        continue;
                    }
                }
                ok += 1;
            }
            Err(e) => {
                fail += 1;
                println!("FAIL {sif} :: {e}");
            }
        }
    }
    (ok, fail, 0)
}

fn collect_sifs(dir: &str, out: &mut Vec<String>) {
    let rd = match std::fs::read_dir(dir) {
        Ok(r) => r,
        Err(_) => return,
    };
    let mut entries: Vec<_> = rd.filter_map(|e| e.ok()).collect();
    entries.sort_by_key(|e| e.file_name());
    for e in entries {
        let p = e.path();
        if p.is_dir() {
            continue;
        }
        if let Some(ext) = p.extension().and_then(|s| s.to_str()) {
            if ext.eq_ignore_ascii_case("sif") {
                out.push(p.to_string_lossy().into_owned());
            }
        }
    }
}
