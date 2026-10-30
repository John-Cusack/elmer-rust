//! `elmer-oracle` command line: plan, run, manifest, and grade.

use elmer_oracle::manifest::{hash_tree, Manifest};
use elmer_oracle::registry::{read_registry, read_registry_under};
use elmer_oracle::runner::{
    discover_tests, known_bins, missing_output_dirs, parse_cmakelists, plan_test, read_startinfo,
    run_plan, run_replay, Arm, SifInject, Staging,
};
use std::path::{Path, PathBuf};

fn usage() -> ! {
    eprintln!("usage:");
    eprintln!("  elmer-oracle plan <tests-dir> <test> [--np N]");
    eprintln!("  elmer-oracle run <build-dir> <tests-dir> <test> <run-dir> [--np N] [--o2] [--o2m] [--max-output L] [--solver CANDIDATE]");
    eprintln!("  elmer-oracle capture <build-dir> <tests-dir> <test> <golden-dir> <manifest> [--np N] [--o2] [--o2m] [--max-output L] [--solver CANDIDATE]");
    eprintln!("  elmer-oracle registry <build-dir> [--label L] [--suite fem/tests|elmerice]");
    eprintln!("  elmer-oracle grade-tree <tests-dir> <golden-root> <candidate-root> <tolerances.json> [--fallback X]");
    eprintln!("  elmer-oracle mutation-suite <tests-dir> <golden-root> <tolerances.json> [--fallback X]");
    eprintln!("  elmer-oracle sigma-table <tests-dir> <out.json> <golden-root> <arm=root>...");
    eprintln!("  elmer-oracle grade-fields <golden.result> <candidate.result> [--tol X] [--tolerances T --test NAME [--file REL]]");
    eprintln!("  elmer-oracle manifest <dir> <out> --origin NAME");
    eprintln!("  elmer-oracle bins <build-dir>");
    eprintln!("  elmer-oracle discover <tests-dir>");
    std::process::exit(2);
}

fn arm_from_build(build_dir: PathBuf, name: &str) -> Arm {
    Arm {
        name: name.to_string(),
        elmergrid: resolve_tool(&build_dir, &["elmergrid/src/ElmerGrid", "fem/src/ElmerGrid"]),
        elmersolver: resolve_tool(&build_dir, &["fem/src/ElmerSolver"]),
        viewfactors: resolve_tool(&build_dir, &["fem/src/ViewFactors"]),
        radiators: resolve_tool(&build_dir, &["fem/src/Radiators"]),
        mesh2d: Some(resolve_tool(&build_dir, &["meshgen2d/src/Mesh2D"])),
        mpiexec: Some(PathBuf::from("mpiexec")),
        with_mpi: true,
        solver_override: false,
        build_dir,
    }
}

/// First existing candidate wins; matches the `-D..._BIN=` paths the
/// configure records in the arm's own `CTestTestfile.cmake` files.
fn resolve_tool(build_dir: &Path, candidates: &[&str]) -> PathBuf {
    for c in candidates {
        let p = build_dir.join(c);
        if p.exists() {
            return p;
        }
    }
    build_dir.join(candidates[0])
}

/// Verdict of one `run`/`capture` invocation.
struct Executed {
    passed: bool,
    arm: String,
}

/// Shared body of `run` and `capture`:
/// `<build-dir> <tests-dir> <test> <run-dir> [--np N] [--o2] [--o2m] [--max-output L]`.
/// `<test>` is a test directory or a registered ctest name (`X_np2`); a
/// registered name supplies its directory and task count. A plan the
/// declarative runner refuses is replayed through `cmake -P` when the arm
/// registers the test, and refused (exit 3) otherwise.
fn execute(args: &[String]) -> Executed {
    let build_dir = PathBuf::from(&args[0]);
    let tests_dir = PathBuf::from(&args[1]);
    let test = args[2].clone();
    let run_dir = PathBuf::from(&args[3]);
    let mut np: Option<u32> = None;
    let mut inj = SifInject::default();
    let mut solver: Option<PathBuf> = None;
    let mut i = 4;
    while i < args.len() {
        match args[i].as_str() {
            "--np" if i + 1 < args.len() => {
                np = Some(args[i + 1].parse().unwrap_or_else(|_| usage()));
                i += 1;
            }
            "--max-output" if i + 1 < args.len() => {
                inj.max_output = Some(args[i + 1].parse().unwrap_or_else(|_| usage()));
                i += 1;
            }
            "--solver" if i + 1 < args.len() => {
                solver = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            "--o2" => inj.o2 = true,
            "--o2m" => inj.o2m = true,
            _ => {}
        }
        i += 1;
    }
    let registered = read_registry(&build_dir)
        .ok()
        .and_then(|reg| reg.into_iter().find(|t| t.name == test));
    let dir_name = match &registered {
        Some(r) if !tests_dir.join(&test).is_dir() => r.dir.clone(),
        _ => test.clone(),
    };
    let np = np.or(registered.as_ref().map(|r| r.np)).unwrap_or(1);
    let dir = tests_dir.join(&dir_name);
    let cmake = std::fs::read_to_string(dir.join("CMakeLists.txt")).unwrap_or_else(|_| usage());
    let mut spec = parse_cmakelists(&dir_name, &cmake);
    spec.startinfo_sif = read_startinfo(&dir).unwrap_or_else(|_| usage());
    let runtest = std::fs::read_to_string(dir.join("runtest.cmake")).unwrap_or_else(|_| usage());
    let mut arm = arm_from_build(build_dir, "cli");
    if let Some(candidate) = solver {
        // Grade a candidate: its ElmerSolver, everything else (ElmerGrid,
        // test-local plugins, registration) from the pristine arm.
        arm.name = format!("candidate {}", candidate.display());
        arm.elmersolver = candidate;
        arm.solver_override = true;
    }
    let plan = plan_test(&spec, &runtest, np, &arm).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1);
    });
    let mut staging = Staging {
        inject: inj.is_active().then_some(inj),
        mkdirs: Vec::new(),
    };
    let outcome = loop {
        let result = match (&plan.custom, &registered) {
            (None, _) => run_plan(&arm, &spec, &plan, &dir, &run_dir, &staging),
            (Some(reason), Some(reg)) if reg.np == np => {
                println!("replay: {reason}");
                run_replay(&arm, reg, &dir, &run_dir, &staging)
            }
            (Some(reason), _) => {
                eprintln!("custom test, refusing: {reason}");
                std::process::exit(3);
            }
        };
        let outcome = result.unwrap_or_else(|e| {
            eprintln!("error: {e}");
            std::process::exit(1);
        });
        // An injected Output File on an in-memory mesh fails to open its
        // `<mesh>/` directory: create it and rerun from a clean stage.
        let missing: Vec<String> = missing_output_dirs(&run_dir)
            .into_iter()
            .filter(|d| !staging.mkdirs.contains(d))
            .collect();
        if outcome.passed || staging.inject.is_none() || missing.is_empty() {
            break outcome;
        }
        println!("mkdirs: {}", missing.join(" "));
        staging.mkdirs.extend(missing);
        std::fs::remove_dir_all(&run_dir).unwrap_or_else(|e| {
            eprintln!("error: clearing {}: {e}", run_dir.display());
            std::process::exit(1);
        });
    };
    println!("passed: {}", outcome.passed);
    if let Some(step) = &outcome.failed_step {
        println!("  failed step: {step}");
    }
    for line in &outcome.compare_lines {
        println!("  {line}");
    }
    Executed {
        passed: outcome.passed,
        arm: arm.name,
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }
    match args[1].as_str() {
        "plan" => {
            if args.len() < 4 {
                usage();
            }
            let tests_dir = PathBuf::from(&args[2]);
            let test = args[3].clone();
            let mut np = 1u32;
            let mut i = 4;
            while i < args.len() {
                if args[i] == "--np" && i + 1 < args.len() {
                    np = args[i + 1].parse().unwrap_or_else(|_| usage());
                    i += 1;
                }
                i += 1;
            }
            let dir = tests_dir.join(&test);
            let cmake = std::fs::read_to_string(dir.join("CMakeLists.txt")).unwrap_or_else(|_| usage());
            let spec = parse_cmakelists(&test, &cmake);
            let runtest = std::fs::read_to_string(dir.join("runtest.cmake")).unwrap_or_else(|_| usage());
            let arm = arm_from_build(PathBuf::from("/build"), "display");
            let plan = plan_test(&spec, &runtest, np, &arm).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            match plan.custom {
                Some(reason) => println!("custom: {reason}"),
                None => {
                    println!("steps: {}", plan.steps.len());
                    for step in &plan.steps {
                        println!("  {step:?}");
                    }
                    println!("passed-file: {}", plan.passed_file);
                }
            }
        }
        "run" => {
            if args.len() < 6 {
                usage();
            }
            let outcome = execute(&args[2..]);
            if !outcome.passed {
                std::process::exit(1);
            }
        }
        "capture" => {
            // Run one test, then freeze its run dir as golden data: a
            // manifest of every file the solver wrote, plus the graded
            // verdict extracted from the logs.
            if args.len() < 7 {
                usage();
            }
            let golden_dir = PathBuf::from(&args[5]);
            let manifest_path = PathBuf::from(&args[6]);
            let run_dir = golden_dir.join("run");
            let mut rest: Vec<String> = args[2..5].to_vec();
            rest.push(run_dir.display().to_string());
            rest.extend(args[7..].iter().cloned());
            let outcome = execute(&rest);
            if !outcome.passed {
                std::process::exit(1);
            }
            let test = &args[4];
            let file_name = manifest_path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();
            let entries = hash_tree(&run_dir, &file_name).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let mut m = Manifest::new(&format!("golden {}", outcome.arm));
            for (path, hash) in entries {
                m.insert(&format!("{test}/{path}"), &hash);
            }
            std::fs::write(&manifest_path, m.serialize()).unwrap_or_else(|_| usage());
            println!("entries: {}", m.entries.len());
        }
        "registry" => {
            if args.len() < 3 {
                usage();
            }
            let mut label: Option<String> = None;
            let mut suite = String::from("fem/tests");
            let mut i = 3;
            while i < args.len() {
                if args[i] == "--label" && i + 1 < args.len() {
                    label = Some(args[i + 1].clone());
                    i += 1;
                } else if args[i] == "--suite" && i + 1 < args.len() {
                    suite = args[i + 1].clone();
                    i += 1;
                }
                i += 1;
            }
            let reg = read_registry_under(&Path::new(&args[2]).join(&suite)).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            for t in reg.iter().filter(|t| label.as_ref().is_none_or(|l| t.has_label(l))) {
                println!("{}\t{}\t{}\t{}", t.name, t.dir, t.np, t.labels.join(";"));
            }
        }
        "grade-fields" => {
            // Grade a candidate binary .result against a golden one with
            // a uniform Tier bar; per-field tables refine this later.
            if args.len() < 5 {
                usage();
            }
            let golden_path = PathBuf::from(&args[2]);
            let candidate_path = PathBuf::from(&args[3]);
            let mut tol = 1e-10f64;
            let mut table_path: Option<String> = None;
            let mut test_name = String::from("unknown");
            let mut file_key = String::from("case.result");
            let mut i = 4;
            while i < args.len() {
                if args[i] == "--tol" && i + 1 < args.len() {
                    tol = args[i + 1].parse().unwrap_or_else(|_| usage());
                    i += 1;
                } else if args[i] == "--tolerances" && i + 1 < args.len() {
                    table_path = Some(args[i + 1].clone());
                    i += 1;
                } else if args[i] == "--test" && i + 1 < args.len() {
                    test_name = args[i + 1].clone();
                    i += 1;
                } else if args[i] == "--file" && i + 1 < args.len() {
                    file_key = args[i + 1].clone();
                    i += 1;
                }
                i += 1;
            }
            let golden = elmer_oracle::result::read_result_file(&golden_path).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let candidate =
                elmer_oracle::result::read_result_file(&candidate_path).unwrap_or_else(|e| {
                    eprintln!("error: {e}");
                    std::process::exit(1);
                });
            let (mismatches, notes) = match table_path {
                Some(p) => {
                    let text = std::fs::read_to_string(&p).unwrap_or_else(|_| usage());
                    let table =
                        elmer_oracle::tolerance::ToleranceTable::parse(&text).unwrap_or_else(|e| {
                            eprintln!("error: {e}");
                            std::process::exit(1);
                        });
                    elmer_oracle::tolerance::grade_with_table(
                        &golden,
                        &candidate,
                        &test_name,
                        &file_key,
                        &table,
                        tol,
                    )
                }
                None => (elmer_oracle::grade::grade_results(&golden, &candidate, tol), Vec::new()),
            };
            for n in &notes {
                println!("note: {n}");
            }
            for m in &mismatches {
                println!("{} step {}: {:?}", m.variable, m.step, m.kind);
            }
            println!("mismatches: {}", mismatches.len());
            if !mismatches.is_empty() {
                std::process::exit(1);
            }
        }
        "sigma-table" => {
            // sigma-table <tests-dir> <out.json> <golden-root> <arm=root>... :
            // every test captured under <golden-root> is measured against
            // the same test under each arm root.
            if args.len() < 6 {
                usage();
            }
            let tests_dir = PathBuf::from(&args[2]);
            let out = PathBuf::from(&args[3]);
            let golden_root = PathBuf::from(&args[4]);
            let arms: Vec<(String, PathBuf)> = args[5..]
                .iter()
                .map(|a| {
                    let (n, p) = a.split_once('=').unwrap_or_else(|| usage());
                    (n.to_string(), PathBuf::from(p))
                })
                .collect();
            let registry = std::fs::read_to_string(golden_root.join("registry.tsv"))
                .unwrap_or_else(|_| usage());
            let mut entries = Vec::new();
            let mut errors = 0;
            for line in registry.lines() {
                let mut cols = line.split('\t');
                let (Some(test), Some(dir)) = (cols.next(), cols.next()) else {
                    continue;
                };
                let arm_runs: Vec<(String, PathBuf)> =
                    arms.iter().map(|(n, r)| (n.clone(), r.join(test).join("run"))).collect();
                match elmer_oracle::sigma::measure_test(
                    test,
                    &golden_root.join(test).join("run"),
                    &tests_dir.join(dir),
                    &arm_runs,
                ) {
                    Ok(e) => entries.extend(e),
                    Err(e) => {
                        eprintln!("error: {e}");
                        errors += 1;
                    }
                }
            }
            let arm_json: serde_json::Map<String, serde_json::Value> = arms
                .iter()
                .map(|(n, r)| (n.clone(), serde_json::json!(r.display().to_string())))
                .collect();
            let json = elmer_oracle::sigma::tolerances_json(
                &entries,
                &serde_json::Value::Object(arm_json),
                "9f6af2f85433fbe8ec71b85c3486ed38b88fcaf0",
                "relative Linf in node space (Perm-mapped) of every captured result file, each sensitivity arm vs the golden arm, per saved step",
            );
            std::fs::write(&out, json).unwrap_or_else(|_| usage());
            for (test, (worst, n)) in elmer_oracle::sigma::per_test_worst(&entries) {
                println!("{test}\t{worst:e}\t{n}");
            }
            eprintln!("entries: {} errors: {errors}", entries.len());
            if errors > 0 {
                std::process::exit(1);
            }
        }
        "grade-tree" => {
            // grade-tree <tests-dir> <golden-root> <candidate-root> <tolerances.json> [--fallback X]:
            // every test captured under <candidate-root> is graded file by
            // file against the golden with the committed Tier bars. The
            // candidate's own TEST.PASSED is never read.
            if args.len() < 6 {
                usage();
            }
            let tests_dir = PathBuf::from(&args[2]);
            let golden_root = PathBuf::from(&args[3]);
            let cand_root = PathBuf::from(&args[4]);
            let text = std::fs::read_to_string(&args[5]).unwrap_or_else(|_| usage());
            let table = elmer_oracle::tolerance::ToleranceTable::parse(&text).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let mut fallback = 1e-10f64;
            if args.len() >= 8 && args[6] == "--fallback" {
                fallback = args[7].parse().unwrap_or_else(|_| usage());
            }
            let registry = std::fs::read_to_string(golden_root.join("registry.tsv"))
                .unwrap_or_else(|_| usage());
            let mut verdicts: std::collections::BTreeMap<&str, usize> = Default::default();
            println!("test\tverdict\tfiles\tmismatches\tdetail");
            for line in registry.lines() {
                let mut cols = line.split('\t');
                let (Some(test), Some(dir)) = (cols.next(), cols.next()) else {
                    continue;
                };
                let cand_run = cand_root.join(test).join("run");
                if !cand_run.is_dir() {
                    continue;
                }
                let golden_run = golden_root.join(test).join("run");
                let rels = elmer_oracle::sigma::result_files(&golden_run, &tests_dir.join(dir))
                    .unwrap_or_default();
                let (mut files, mut bad, mut detail) = (0, 0, Vec::new());
                for rel in &rels {
                    let golden = match elmer_oracle::result::read_result_file(&golden_run.join(rel)) {
                        Ok(g) if g.steps.iter().any(|s| !s.vars.is_empty()) => g,
                        _ => continue,
                    };
                    files += 1;
                    match elmer_oracle::result::read_result_file(&cand_run.join(rel)) {
                        Err(_) => {
                            bad += 1;
                            detail.push(format!("{rel}: no candidate file"));
                        }
                        Ok(cand) => {
                            let (m, _) = elmer_oracle::tolerance::grade_with_table(
                                &golden, &cand, test, rel, &table, fallback,
                            );
                            if let Some(first) = m.first() {
                                detail.push(format!(
                                    "{rel}: {} step {}: {:?}",
                                    first.variable, first.step, first.kind
                                ));
                            }
                            bad += m.len();
                        }
                    }
                }
                let verdict = if files == 0 {
                    "no-golden-fields"
                } else if bad == 0 {
                    "tier-a-pass"
                } else if detail.iter().all(|d| d.ends_with("no candidate file")) {
                    "no-candidate-output"
                } else {
                    "tier-a-FAIL"
                };
                *verdicts.entry(verdict).or_default() += 1;
                println!("{test}\t{verdict}\t{files}\t{bad}\t{}", detail.join("; "));
            }
            let summary: Vec<String> = verdicts.iter().map(|(k, v)| format!("{k}={v}")).collect();
            eprintln!("{}", summary.join(" "));
            if verdicts.contains_key("tier-a-FAIL") {
                std::process::exit(1);
            }
        }
        "mutation-suite" => {
            // mutation-suite <tests-dir> <golden-root> <tolerances.json> [--fallback X]:
            // every mutation of every captured result file must fail grading.
            if args.len() < 5 {
                usage();
            }
            let tests_dir = PathBuf::from(&args[2]);
            let golden_root = PathBuf::from(&args[3]);
            let text = std::fs::read_to_string(&args[4]).unwrap_or_else(|_| usage());
            let table = elmer_oracle::tolerance::ToleranceTable::parse(&text).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let mut fallback = 1e-10f64;
            if args.len() >= 7 && args[5] == "--fallback" {
                fallback = args[6].parse().unwrap_or_else(|_| usage());
            }
            let registry = std::fs::read_to_string(golden_root.join("registry.tsv"))
                .unwrap_or_else(|_| usage());
            let (mut caught, mut escaped, mut below, mut vacuous, mut files) = (0, 0, 0, 0, 0);
            for line in registry.lines() {
                let mut cols = line.split('\t');
                let (Some(test), Some(dir)) = (cols.next(), cols.next()) else {
                    continue;
                };
                let run = golden_root.join(test).join("run");
                let rels = elmer_oracle::sigma::result_files(&run, &tests_dir.join(dir))
                    .unwrap_or_else(|e| {
                        eprintln!("error: {e}");
                        std::process::exit(1);
                    });
                for rel in rels {
                    let golden = elmer_oracle::result::read_result_file(&run.join(&rel))
                        .unwrap_or_else(|e| {
                            eprintln!("error: {e}");
                            std::process::exit(1);
                        });
                    files += 1;
                    for t in elmer_oracle::mutation::run_suite(&golden, test, &rel, &table, fallback) {
                        use elmer_oracle::mutation::Outcome::*;
                        let m = t.mutation.name();
                        match t.outcome {
                            Caught => caught += 1,
                            Escaped => {
                                escaped += 1;
                                println!("ESCAPED\t{test}\t{rel}\t{m}\t{}", t.target);
                            }
                            BelowFloor => {
                                below += 1;
                                println!("below-floor\t{test}\t{rel}\t{m}\t{}", t.target);
                            }
                            Vacuous(why) => {
                                vacuous += 1;
                                println!("vacuous\t{test}\t{rel}\t{m}\t{}\t{why}", t.target);
                            }
                        }
                    }
                }
            }
            println!(
                "files: {files} caught: {caught} escaped: {escaped} below-floor: {below} vacuous: {vacuous}"
            );
            if escaped > 0 {
                std::process::exit(1);
            }
        }
        "manifest" => {
            if args.len() < 4 {
                usage();
            }
            let dir = PathBuf::from(&args[2]);
            let out = PathBuf::from(&args[3]);
            let mut origin = String::from("elmer-oracle");
            let mut i = 4;
            while i < args.len() {
                if args[i] == "--origin" && i + 1 < args.len() {
                    origin = args[i + 1].clone();
                    i += 1;
                }
                i += 1;
            }
            let file_name = out
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();
            let entries = hash_tree(&dir, &file_name).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let mut m = Manifest::new(&origin);
            for (path, hash) in entries {
                m.insert(&path, &hash);
            }
            std::fs::write(&out, m.serialize()).unwrap_or_else(|_| usage());
            println!("entries: {}", m.entries.len());
        }
        "bins" => {
            if args.len() < 3 {
                usage();
            }
            let arm = arm_from_build(PathBuf::from(&args[2]), "display");
            let mut pairs: Vec<_> = known_bins(&arm).into_iter().collect();
            pairs.sort();
            for (var, path) in pairs {
                println!("{var}={path}");
            }
        }
        "discover" => {
            if args.len() < 3 {
                usage();
            }
            let tests = discover_tests(&PathBuf::from(&args[2])).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            for t in tests {
                println!("{t}");
            }
        }
        _ => usage(),
    }
}
