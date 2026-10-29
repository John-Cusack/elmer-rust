//! Byte-for-byte checks against stdout captured from the pristine ElmerSolver
//! (26.2-devel, Rev 9f6af2f85) run serially on three quick tests of
//! `fem/tests`: `heateq` (steady, two instances), `heateq_bdf2` (transient)
//! and `HeatAniso` (steady with SS lines). The fixtures hold the `Time:` and
//! `ComputeChange` lines of those runs, untouched.

mod support;

use elmer_cli::logger::{compute_change_line, time_message};
use elmer_cli::version::{print_probe, VersionInfo};
use std::time::{Duration, UNIX_EPOCH};

fn fixture(name: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

/// `NS (ITER=1) (NRM,RELC): ( 0.7E-02  2.0     ) :: heat` -> parts.
fn parse_change(line: &str) -> (bool, u32, f64, f64, String) {
    let body = line.strip_prefix("ComputeChange: ").expect("prefix");
    let steady = body.starts_with("SS ");
    let iter = body.split("ITER=").nth(1).unwrap().split(')').next().unwrap().parse().unwrap();
    let nums = body.split("): (").nth(1).unwrap().split(" ) :: ").next().unwrap();
    let mut it = nums.split_whitespace().map(|t| t.parse::<f64>().unwrap());
    let name = body.split(" ) :: ").nth(1).unwrap().to_string();
    (steady, iter, it.next().unwrap(), it.next().unwrap(), name)
}

fn regenerate(fixture_name: &str) -> usize {
    let mut n = 0;
    for line in fixture(fixture_name).lines() {
        if let Some(rest) = line.strip_prefix("MAIN: Time: ") {
            let (steps, t) = rest.split_once(':').unwrap();
            let (cum, total) = steps.split_once('/').unwrap();
            let got = time_message(cum.parse().unwrap(), total.parse().unwrap(), t.trim().parse().unwrap());
            assert_eq!(format!("MAIN: {got}"), line);
        } else {
            let (steady, iter, norm, change, name) = parse_change(line);
            assert_eq!(compute_change_line(&name, iter, norm, change, steady), line);
        }
        n += 1;
    }
    n
}

#[test]
fn heateq_steady_change_lines() {
    assert_eq!(regenerate("pristine_heateq.txt"), 13);
}

#[test]
fn heateq_bdf2_transient_time_and_change_lines() {
    assert_eq!(regenerate("pristine_heateq_bdf2.txt"), 24);
}

#[test]
fn heat_aniso_ss_and_ns_lines() {
    assert_eq!(regenerate("pristine_HeatAniso.txt"), 4);
}

#[test]
fn logger_emits_the_same_bytes_through_info() {
    let (log, buf) = support::logger();
    let (steady, iter, norm, change, name) = parse_change(
        "ComputeChange: NS (ITER=2) (NRM,RELC): ( 0.14447068E-01  0.0000000     ) :: heat",
    );
    assert!(!steady);
    log.compute_change(&name, iter, norm, change, false);
    log.time_step(2, 8, 0.2);
    assert_eq!(
        buf.text(),
        "ComputeChange: NS (ITER=2) (NRM,RELC): ( 0.14447068E-01  0.0000000     ) :: heat\n\
         MAIN: Time: 2/8:   2.000E-01\n"
    );
}

/// The `-v` probe equals the pristine one after masking the date, the timing
/// numbers and the one added backend line.
#[test]
fn probe_matches_pristine_format() {
    let (log, buf) = support::logger();
    print_probe(&log, &VersionInfo::default(), false, false, UNIX_EPOCH + Duration::from_secs(1_793_269_507), (0.0, 0.4));
    let mask = |t: &str| -> Vec<String> {
        t.lines()
            .filter(|l| !l.contains("Backend:"))
            .map(|l| {
                if l.starts_with("ELMER SOLVER (v") || l.starts_with("ELMER SOLVER FINISHED") {
                    l.split(" AT: ").next().unwrap().to_string()
                } else if l.starts_with("SOLVER TOTAL TIME") {
                    "SOLVER TOTAL TIME".to_string()
                } else {
                    l.to_string()
                }
            })
            .collect()
    };
    assert_eq!(mask(&buf.text()), mask(&fixture("pristine_v.txt")));
    assert!(buf.text().contains("MAIN:  Backend: Rust-native elmer-cli"));
}

#[test]
fn probe_silent_flags() {
    let (log, buf) = support::logger();
    print_probe(&log, &VersionInfo::default(), true, true, UNIX_EPOCH, (0.0, 0.0));
    assert_eq!(
        buf.text(),
        " Started inside library code\n Going parallel initialization\n\
         ParCommInit: Initialize #PEs: 1\n Back from parallel initialization\n"
    );
}
