//! The harvested campaign cases (report §5.2 item 11, `harvest/`): every
//! case plans as an ordinary oracle test, and the capture index covers
//! exactly the case list.

use elmer_oracle::runner::{discover_tests, parse_cmakelists, plan_test, read_runtest, read_startinfo, Arm, Step};
use std::path::PathBuf;

fn harvest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("harvest")
}

fn tsv(name: &str) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(harvest().join(name)).unwrap();
    text.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect()
}

fn arm() -> Arm {
    let p = PathBuf::from("/nonexistent");
    Arm {
        name: String::from("plan-only"),
        build_dir: p.clone(),
        elmergrid: p.clone(),
        elmersolver: p.clone(),
        viewfactors: p.clone(),
        radiators: p,
        mesh2d: None,
        mpiexec: None,
        with_mpi: false,
        solver_override: false,
    }
}

#[test]
fn every_harvested_case_plans() {
    let dir = harvest().join("cases");
    let tests = discover_tests(&dir).unwrap();
    let listed: Vec<String> = tsv("cases/cases.tsv").into_iter().filter(|r| r[7].is_empty()).map(|r| r[0].clone()).collect();
    let mut sorted = listed.clone();
    sorted.sort();
    let mut found = tests.clone();
    found.sort();
    assert_eq!(found, sorted, "case dirs and cases.tsv disagree");
    for t in &tests {
        let d = dir.join(t);
        let cmake = std::fs::read_to_string(d.join("CMakeLists.txt")).unwrap();
        let mut spec = parse_cmakelists(t, &cmake);
        spec.startinfo_sif = read_startinfo(&d).unwrap();
        let sif = spec.startinfo_sif.clone().unwrap_or_else(|| panic!("{t}: no startinfo"));
        assert!(d.join(&sif).is_file(), "{t}: {sif} missing");
        let plan = plan_test(&spec, &read_runtest(&d).unwrap(), 1, &arm()).unwrap_or_else(|e| panic!("{t}: {e}"));
        assert!(plan.custom.is_none(), "{t}: needs a cmake replay");
        assert_eq!(plan.steps.iter().filter(|s| matches!(s, Step::Solver { .. })).count(), 1, "{t}");
    }
}

#[test]
fn index_covers_every_case_once() {
    let cases: Vec<String> = tsv("cases/cases.tsv").into_iter().map(|r| r[0].clone()).collect();
    let index = tsv("index.tsv");
    let ids: Vec<String> = index.iter().map(|r| r[0].clone()).collect();
    assert_eq!(ids, cases);
    for r in &index {
        match r[1].as_str() {
            "pass" => assert_eq!(r[3].len(), 64, "{}: no manifest hash", r[0]),
            "FAIL" => assert!(!r[6].is_empty(), "{}: failure without a reason", r[0]),
            "skipped" => {}
            v => panic!("{}: verdict {v}", r[0]),
        }
    }
    let defects = tsv("defects.tsv");
    assert_eq!(defects.len(), 15);
    assert!(defects.iter().all(|d| d.len() == 7 && d.iter().all(|c| !c.is_empty())));
}
