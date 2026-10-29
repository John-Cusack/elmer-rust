//! `elmer-oracle` command line: plan, run, manifest, and grade.

use elmer_oracle::manifest::{hash_tree, Manifest};
use elmer_oracle::runner::{
    discover_tests, known_bins, parse_cmakelists, plan_test, read_startinfo, run_plan, Arm,
};
use std::path::PathBuf;

fn usage() -> ! {
    eprintln!("usage:");
    eprintln!("  elmer-oracle plan <tests-dir> <test> [--np N]");
    eprintln!("  elmer-oracle run <build-dir> <tests-dir> <test> <run-dir> [--np N]");
    eprintln!("  elmer-oracle manifest <dir> <out> --origin NAME");
    eprintln!("  elmer-oracle bins <build-dir>");
    std::process::exit(2);
}

fn arm_from_build(build_dir: PathBuf, name: &str) -> Arm {
    Arm {
        name: name.to_string(),
        elmergrid: build_dir.join("fem/src/ElmerGrid"),
        elmersolver: build_dir.join("fem/src/ElmerSolver"),
        viewfactors: build_dir.join("fem/src/ViewFactors"),
        radiators: build_dir.join("fem/src/Radiators"),
        mesh2d: Some(build_dir.join("meshgen2d/src/Mesh2D")),
        mpiexec: Some(PathBuf::from("mpiexec")),
        with_mpi: true,
        build_dir,
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
            let build_dir = PathBuf::from(&args[2]);
            let tests_dir = PathBuf::from(&args[3]);
            let test = args[4].clone();
            let run_dir = PathBuf::from(&args[5]);
            let mut np = 1u32;
            let mut i = 6;
            while i < args.len() {
                if args[i] == "--np" && i + 1 < args.len() {
                    np = args[i + 1].parse().unwrap_or_else(|_| usage());
                    i += 1;
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
            match run_plan(&arm, &spec, &plan, &dir, &run_dir) {
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
