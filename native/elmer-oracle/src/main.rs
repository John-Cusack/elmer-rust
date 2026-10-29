//! `elmer-oracle` command line: plan, run, manifest, and grade.

use elmer_oracle::manifest::{hash_tree, Manifest};
use elmer_oracle::runner::{
    discover_tests, known_bins, parse_cmakelists, plan_test, read_startinfo, run_plan, Arm,
    SifInject,
};
use std::path::{Path, PathBuf};

fn usage() -> ! {
    eprintln!("usage:");
    eprintln!("  elmer-oracle plan <tests-dir> <test> [--np N]");
    eprintln!("  elmer-oracle run <build-dir> <tests-dir> <test> <run-dir> [--np N]");
    eprintln!("  elmer-oracle capture <build-dir> <tests-dir> <test> <golden-dir> <manifest> [--np N]");
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
            let build_dir = PathBuf::from(&args[2]);
            let tests_dir = PathBuf::from(&args[3]);
            let test = args[4].clone();
            let run_dir = PathBuf::from(&args[5]);
            let mut np = 1u32;
            let mut inj = SifInject { o2: false, o2m: false };
            let mut i = 6;
            while i < args.len() {
                if args[i] == "--np" && i + 1 < args.len() {
                    np = args[i + 1].parse().unwrap_or_else(|_| usage());
                    i += 1;
                } else if args[i] == "--o2" {
                    inj.o2 = true;
                } else if args[i] == "--o2m" {
                    inj.o2m = true;
                }
                i += 1;
            }
            let dir = tests_dir.join(&test);
            let cmake = std::fs::read_to_string(dir.join("CMakeLists.txt")).unwrap_or_else(|_| usage());
            let mut spec = parse_cmakelists(&test, &cmake);
            spec.startinfo_sif = read_startinfo(&dir).unwrap_or_else(|_| usage());
            let runtest = std::fs::read_to_string(dir.join("runtest.cmake")).unwrap_or_else(|_| usage());
            let arm = arm_from_build(build_dir, "cli");
            let plan = plan_test(&spec, &runtest, np, &arm).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            if let Some(reason) = plan.custom {
                eprintln!("custom test, refusing: {reason}");
                std::process::exit(3);
            }
            let inj_opt = if inj.o2 || inj.o2m { Some(&inj) } else { None };
            match run_plan(&arm, &spec, &plan, &dir, &run_dir, inj_opt) {
                Ok(outcome) => {
                    println!("passed: {}", outcome.passed);
                    for line in &outcome.compare_lines {
                        println!("  {line}");
                    }
                    if !outcome.passed {
                        std::process::exit(1);
                    }
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    std::process::exit(1);
                }
            }
        }
        "capture" => {
            // Run one test, then freeze its run dir as golden data: a
            // manifest of every file the solver wrote, plus the graded
            // verdict extracted from the logs.
            if args.len() < 7 {
                usage();
            }
            let build_dir = PathBuf::from(&args[2]);
            let tests_dir = PathBuf::from(&args[3]);
            let test = args[4].clone();
            let golden_dir = PathBuf::from(&args[5]);
            let manifest_path = PathBuf::from(&args[6]);
            let mut np = 1u32;
            let mut inj = SifInject { o2: false, o2m: false };
            let mut i = 7;
            while i < args.len() {
                if args[i] == "--np" && i + 1 < args.len() {
                    np = args[i + 1].parse().unwrap_or_else(|_| usage());
                    i += 1;
                } else if args[i] == "--o2" {
                    inj.o2 = true;
                } else if args[i] == "--o2m" {
                    inj.o2m = true;
                }
                i += 1;
            }
            let dir = tests_dir.join(&test);
            let cmake = std::fs::read_to_string(dir.join("CMakeLists.txt")).unwrap_or_else(|_| usage());
            let mut spec = parse_cmakelists(&test, &cmake);
            spec.startinfo_sif = read_startinfo(&dir).unwrap_or_else(|_| usage());
            let runtest = std::fs::read_to_string(dir.join("runtest.cmake")).unwrap_or_else(|_| usage());
            let arm = arm_from_build(build_dir, "cli");
            let plan = plan_test(&spec, &runtest, np, &arm).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            if let Some(reason) = plan.custom {
                eprintln!("custom test, refusing: {reason}");
                std::process::exit(3);
            }
            let run_dir = golden_dir.join("run");
            let inj_opt = if inj.o2 || inj.o2m { Some(&inj) } else { None };
            let outcome = run_plan(&arm, &spec, &plan, &dir, &run_dir, inj_opt).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            println!("passed: {}", outcome.passed);
            for line in &outcome.compare_lines {
                println!("  {line}");
            }
            if !outcome.passed {
                std::process::exit(1);
            }
            let file_name = manifest_path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();
            let entries = hash_tree(&run_dir, &file_name).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                std::process::exit(1);
            });
            let mut m = Manifest::new(&format!("golden {}", arm.name));
            for (path, hash) in entries {
                m.insert(&format!("{test}/{path}"), &hash);
            }
            std::fs::write(&manifest_path, m.serialize()).unwrap_or_else(|_| usage());
            println!("entries: {}", m.entries.len());
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
