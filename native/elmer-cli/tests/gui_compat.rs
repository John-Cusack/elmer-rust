//! ElmerGUI parsing of solver stdout, ported verbatim as test data, applied to
//! Rust-formatted lines.
//!
//! Provenance (all under `ElmerGUI/Application/src/` of the reference tree):
//! * `mainwindow.cpp:6930-6978`: `Time:` and `ComputeChange` handling.
//! * `solverlogwindow.cpp:97-115`: highlight regexes (`REG_EXP_CLASS`,
//!   case-insensitive).
//! * `mainwindow.cpp:1236-1281`: `ElmerSolver -v` / `ElmerSolver_mpi -v` probes.
//!
//! The regexes are all of the form `^(.*)LITERAL(.*)$` (plus `\bWORD\b`), so they
//! are evaluated here as case-insensitive substring / whole-word tests; no regex
//! crate is on the dependency allowlist.

mod support;

use elmer_core::Logger;

// solverlogwindow.cpp:97  "^(.*)\\bWARNING\\b(.*)$"
// solverlogwindow.cpp:102 "^(.*)\\bERROR\\b(.*)$"
fn word_ci(text: &str, word: &str) -> bool {
    let t = text.to_ascii_lowercase();
    let w = word.to_ascii_lowercase();
    let is_word = |c: Option<char>| c.is_some_and(|c| c.is_alphanumeric() || c == '_');
    t.match_indices(&w).any(|(i, _)| {
        !is_word(t[..i].chars().next_back()) && !is_word(t[i + w.len()..].chars().next())
    })
}

// solverlogwindow.cpp:107-110
const END_PATTERNS: [&str; 4] = [
    " Elmer Solver: ALL DONE ", // "^(.*) Elmer Solver: ALL DONE (.*)$"
    "ElmerSolver: The end",     // "^(.*)ElmerSolver: The end(.*)$"
    "SOLVER TOTAL TIME",        // "^(.*)SOLVER TOTAL TIME(.*)$"
    "ELMER SOLVER FINISHED AT:", // "^(.*)ELMER SOLVER FINISHED AT:(.*)$"
];

fn end_pattern_matches(line: &str, i: usize) -> bool {
    line.to_ascii_lowercase().contains(&END_PATTERNS[i].to_ascii_lowercase())
}

/// mainwindow.cpp:6930-6937: `tmp.contains("Time:")`, last space-separated token.
fn gui_time(line: &str) -> Option<f64> {
    let tmp = line.trim();
    if !tmp.contains("Time:") {
        return None;
    }
    Some(tmp.split(' ').next_back().unwrap().parse().unwrap_or(0.0))
}

/// mainwindow.cpp:6939-6975: returns (name, res1 = norm, res2 = relative change).
fn gui_compute_change(line: &str) -> Option<(String, f64, f64)> {
    let tmp = line.trim();
    if !tmp.contains("ComputeChange") {
        return None;
    }
    let name = tmp.split(':').next_back().unwrap().trim().to_string();
    let n = tmp.find("NRM,RELC").unwrap();
    let rest = &tmp[n..];
    let parts: Vec<&str> = rest.split('(').collect();
    let (mut res1, mut res2) = (0.0, 0.0);
    if parts.len() >= 2 {
        let tmp2 = parts[1].trim();
        let toks: Vec<&str> = tmp2.split(' ').collect();
        res1 = toks[0].trim().parse().unwrap_or(0.0);
        let mut pos = 1;
        while toks[pos].trim().is_empty() {
            pos += 1;
            if pos > toks.len() {
                break;
            }
        }
        res2 = f64::max(toks[pos].trim().parse().unwrap_or(0.0), 1.0e-16);
    }
    Some((name, res1, res2))
}

#[test]
fn gui_reads_back_compute_change_values() {
    let (log, buf) = support::logger();
    log.compute_change("heat equation", 3, 0.93065875, 0.62474118, true);
    log.compute_change("heat", 1, 0.75696824e-2, 2.0, false);
    log.compute_change("heat", 2, 12.5, 0.0, false);
    let out = buf.text();
    let got: Vec<_> = out.lines().map(|l| gui_compute_change(l).unwrap()).collect();
    assert_eq!(got[0], ("heat equation".to_string(), 0.93065875, 0.62474118));
    assert_eq!(got[1], ("heat".to_string(), 0.75696824e-2, 2.0));
    // GUI floors the relative change at 1e-16 for the log plot.
    assert_eq!(got[2], ("heat".to_string(), 12.5, 1.0e-16));
}

#[test]
fn gui_reads_back_time_value() {
    let (log, buf) = support::logger();
    log.time_step(5, 8, 0.5);
    assert_eq!(gui_time(buf.text().lines().next().unwrap()), Some(0.5));
    assert_eq!(gui_time("ComputeChange: NS"), None);
}

#[test]
fn highlight_patterns_match_rust_lines() {
    let (log, buf) = support::logger();
    log.warn("Solver", "may require \"Procedure\" to operate!");
    log.fatal("MAIN", "Unable to find input file [x.sif], can not execute.");
    log.info("MAIN", "*** Elmer Solver: ALL DONE ***", 3);
    log.raw_line("SOLVER TOTAL TIME(CPU,REAL):         0.02        0.40");
    log.raw_line("ELMER SOLVER FINISHED AT: 2026/10/29 10:25:08");
    let out = buf.text();
    let l: Vec<&str> = out.lines().collect();
    assert!(word_ci(l[0], "WARNING") && !word_ci(l[0], "ERROR"));
    assert!(word_ci(l[1], "ERROR"));
    assert!(end_pattern_matches(l[2], 0));
    assert!(end_pattern_matches(l[3], 2));
    assert!(end_pattern_matches(l[4], 3));
    // Pristine "MAIN: The end" does not match pattern 1 either (it needs
    // "ElmerSolver: The end"); keep that behaviour identical.
    assert!(!end_pattern_matches("MAIN: The end", 1));
}

#[test]
fn word_boundary_helper() {
    assert!(word_ci("an error occurred", "ERROR"));
    assert!(!word_ci("errors galore", "ERROR"));
    assert!(!word_ci("terror", "ERROR"));
}

/// mainwindow.cpp:1240-1243,1278-1281: the probe needs the process to start;
/// the Rust binary does and exits 0.
#[test]
fn both_binaries_start_for_the_gui_probe() {
    for exe in [env!("CARGO_BIN_EXE_ElmerSolver"), env!("CARGO_BIN_EXE_ElmerSolver_mpi")] {
        let out = std::process::Command::new(exe)
            .arg("-v")
            .env_remove("ELMER_FORTRAN_SOLVER")
            .output()
            .unwrap();
        assert_eq!(out.status.code(), Some(0), "{exe}");
        let text = String::from_utf8(out.stdout).unwrap();
        assert!(text.contains("MAIN: Version: 26.2-devel (Rev: "));
    }
}
