//! `elmer-oracle` command line: plan, run, manifest, and grade.

use elmer_oracle::manifest::{hash_tree, Manifest};
use elmer_oracle::registry::read_registry_under;
use elmer_oracle::runner::{
    discover_tests, known_bins, missing_output_dirs, parse_cmakelists, plan_test, read_runtest, read_startinfo,
    run_plan, run_replay, Arm, SifInject, Staging, Step,
};
use std::path::{Path, PathBuf};

fn usage() -> ! {
    eprintln!("usage:");
    eprintln!("  elmer-oracle plan <tests-dir> <test> [--np N]");
    eprintln!("  elmer-oracle run <build-dir> <tests-dir> <test> <run-dir> [--np N] [--o2] [--o2m] [--per-run] [--max-output L] [--solver CANDIDATE] [--exit-verdict]");
    eprintln!("  elmer-oracle capture <build-dir> <tests-dir> <test> <golden-dir> <manifest> [--np N] [--o2] [--o2m] [--per-run] [--max-output L] [--solver CANDIDATE] [--exit-verdict]");
    eprintln!("  elmer-oracle registry <build-dir> [--label L|'!L'] [--suite fem/tests|elmerice]");
    eprintln!("  elmer-oracle refnorm-validate <build-dir> <tests-dir> <golden-root> <out-root>");
    eprintln!("  elmer-oracle grade-tree <tests-dir> <golden-root> <candidate-root> <tolerances.json> [--fallback X]");
    eprintln!("  elmer-oracle mutation-suite <tests-dir> <golden-root> <tolerances.json> [--fallback X] [--all]");
    eprintln!("  elmer-oracle sigma-table <tests-dir> <out.json> <golden-root> <arm=root>...");
    eprintln!("  elmer-oracle grade-fields <golden.result> <candidate.result> [--tol X] [--tolerances T --test NAME [--file REL]]");
    eprintln!("  elmer-oracle manifest <dir> <out> --origin NAME");
    eprintln!("  elmer-oracle o7-replay <build-dir> <shim.so> <in-dir> <out-dir>");
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
    let mut exit_verdict = false;
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
            "--per-run" => inj.per_run = true,
            "--exit-verdict" => exit_verdict = true,
            _ => {}
        }
        i += 1;
    }
    let registered = [build_dir.join("fem").join("tests"), build_dir.join("elmerice")]
        .iter()
        .filter_map(|root| read_registry_under(root).ok())
        .find_map(|reg| reg.into_iter().find(|t| t.name == test));
    let dir_name = match &registered {
        Some(r) if !tests_dir.join(&test).is_dir() => r.dir.clone(),
        _ => test.clone(),
    };
    let np = np.or(registered.as_ref().map(|r| r.np)).unwrap_or(1);
    let dir = tests_dir.join(&dir_name);
    let cmake = std::fs::read_to_string(dir.join("CMakeLists.txt")).unwrap_or_else(|_| usage());
    let mut spec = parse_cmakelists(&dir_name, &cmake);
    spec.startinfo_sif = read_startinfo(&dir).unwrap_or_else(|_| usage());
    let runtest = read_runtest(&dir).unwrap_or_else(|_| usage());
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
        refnorms: Vec::new(),
        run_refnorms: Vec::new(),
        exit_verdict,
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

/// A parsed golden result file: (path relative to the run dir, file,
/// preamble, mesh nodes).
type Parsed = (String, elmer_oracle::result::ResultFile, Vec<elmer_oracle::result::PreambleVar>, Option<usize>);

/// `(solver, golden norm)` for every solver of `solvers` the grader can
/// rebuild from `files`, with notes for the rest.
fn resolve_solvers(
    solvers: &[elmer_oracle::refnorm::RefSolver],
    files: &[&Parsed],
    equations: &[(u32, String)],
    notes: &mut Vec<String>,
    non_nodal: &mut bool,
) -> Vec<(u32, f64)> {
    let mut refnorms = Vec::new();
    let parsed = files;
    for s in solvers {
        // The variable must live in exactly one result file: several
        // meshes carrying one name leave the solver's own copy ambiguous.
        // A solver without `Variable` takes the preamble's single variable
        // of its equation, which must agree across the files it appears in.
        let mut found = Vec::new();
        let mut why = String::from("not in any result file");
        let mut inferred = std::collections::BTreeSet::new();
        for (_, r, vars, nodes) in parsed.iter().copied() {
            let var = match &s.variable {
                Ok(v) => v.clone(),
                Err(e) if e == elmer_oracle::refnorm::NO_VARIABLE => {
                    match elmer_oracle::refnorm::default_variable(vars, s, equations) {
                        Ok(v) => {
                            inferred.insert(v.clone());
                            v
                        }
                        Err(e) => {
                            why = e;
                            continue;
                        }
                    }
                }
                Err(_) => break,
            };
            if !vars.iter().any(|v| v.base_name() == var) {
                continue;
            }
            *non_nodal |= elmer_oracle::refnorm::non_nodal(vars, *nodes, &var);
            match elmer_oracle::refnorm::golden_norm(r, vars, &var, |d| s.norm_opts(d)) {
                Ok(n) => found.push(n),
                Err(e) => {
                    why = e;
                    found.push(f64::NAN);
                }
            }
        }
        if inferred.len() > 1 {
            notes.push(format!("solver {}: default variable differs between result files", s.number));
            continue;
        }
        let var = match (&s.variable, inferred.first()) {
            (Ok(v), _) => v.clone(),
            (Err(_), Some(v)) => format!("{v} (default)"),
            (Err(e), None) if e == elmer_oracle::refnorm::NO_VARIABLE && !parsed.is_empty() => {
                notes.push(format!("solver {}: {why}", s.number));
                continue;
            }
            (Err(e), None) => {
                notes.push(format!("solver {}: {e}", s.number));
                continue;
            }
        };
        match found.as_slice() {
            [n] if !n.is_nan() => refnorms.push((s.number, *n)),
            [_] => notes.push(format!("solver {}: {var}: {why}", s.number)),
            [] => notes.push(format!("solver {}: {var} {why}", s.number)),
            _ => notes.push(format!("solver {}: {var} in several result files", s.number)),
        }
        if !inferred.is_empty() && matches!(found.as_slice(), [n] if !n.is_nan()) {
            notes.push(format!("solver {}: {var}", s.number));
        }
    }
    refnorms
}

/// [`resolve_solvers`] for an `np`-rank run: each result file is written
/// once per rank (`<file>.<rank>`), and the norm is the parallel branch of
/// `ComputeNorm` over all ranks' values, reduced in `order`.
#[allow(clippy::too_many_arguments)]
fn resolve_solvers_mpi(
    solvers: &[elmer_oracle::refnorm::RefSolver],
    parsed: &[Parsed],
    np: usize,
    equations: &[(u32, String)],
    notes: &mut Vec<String>,
    order: elmer_oracle::norm::ReduceOrder,
    golden_run: &Path,
    text: &str,
) -> Vec<(u32, f64)> {
    // Files by base name, ranks in order; only complete rank sets count.
    let mut groups: std::collections::BTreeMap<String, Vec<Option<&Parsed>>> = Default::default();
    for p in parsed {
        let Some((base, rank)) = p.0.rsplit_once('.') else { continue };
        let Ok(rank) = rank.parse::<usize>() else { continue };
        if rank < np {
            groups.entry(base.to_string()).or_insert_with(|| vec![None; np])[rank] = Some(p);
        }
    }
    let groups: Vec<(String, Vec<&Parsed>)> = groups
        .into_iter()
        .filter_map(|(b, ranks)| Some((b, ranks.into_iter().collect::<Option<Vec<_>>>()?)))
        .collect();
    let mut refnorms = Vec::new();
    for s in solvers {
        let var = match &s.variable {
            Ok(v) => v.clone(),
            Err(e) if e == elmer_oracle::refnorm::NO_VARIABLE => {
                match groups.first().map(|(_, g)| elmer_oracle::refnorm::default_variable(&g[0].2, s, equations)) {
                    Some(Ok(v)) => v,
                    Some(Err(e)) => {
                        notes.push(format!("solver {}: {e}", s.number));
                        continue;
                    }
                    None => {
                        notes.push(format!("solver {}: no per-rank result files", s.number));
                        continue;
                    }
                }
            }
            Err(e) => {
                notes.push(format!("solver {}: {e}", s.number));
                continue;
            }
        };
        let mine: Vec<&(String, Vec<&Parsed>)> =
            groups.iter().filter(|(_, g)| g.iter().all(|p| p.2.iter().any(|v| v.base_name() == var))).collect();
        let [(base, ranks)] = mine.as_slice() else {
            notes.push(format!("solver {}: {var} in {} per-rank result sets", s.number, mine.len()));
            continue;
        };
        let mut parts = Vec::with_capacity(np);
        let mut dofs = None;
        let mut why = None;
        // ElmerGrid's partitioning sits beside the per-rank files; without
        // one (Parallel Slices) every rank solves a copy of the whole mesh.
        let partitioning = golden_run
            .join(std::path::Path::new(base.as_str()).parent().unwrap_or(std::path::Path::new("")))
            .join(format!("partitioning.{np}"));
        let consistent = elmer_oracle::refnorm::consistent(text, s.number);
        for (rank, p) in ranks.iter().enumerate() {
            match elmer_oracle::refnorm::golden_values(&p.1, &p.2, &var) {
                Ok((values, perm, d)) if dofs.is_none_or(|x| x == d) => {
                    dofs = Some(d);
                    let mut owned = vec![true; values.len()];
                    if consistent && partitioning.is_dir() {
                        match elmer_oracle::refnorm::owned_nodes(&partitioning, rank) {
                            // Ownership is known per mesh node of the partition
                            // files only: a local mesh Elmer extruded, refined or
                            // extended, or edge/p/DG unknowns, are out of class.
                            Ok(nodes) if Some(nodes.len()) != p.3 => {
                                why = Some(format!(
                                    "local mesh has {:?} nodes, partition files {} (mesh changed after load)",
                                    p.3,
                                    nodes.len()
                                ))
                            }
                            Ok(nodes) if perm.len() != nodes.len() => {
                                why = Some(String::from("non-nodal unknowns in parallel (ownership unknown)"))
                            }
                            Ok(nodes) => {
                                for (i, &k) in perm.iter().enumerate() {
                                    if k > 0 {
                                        let own = nodes.get(i).copied().unwrap_or(true);
                                        for c in 0..d {
                                            if let Some(o) = owned.get_mut(d * (k as usize - 1) + c) {
                                                *o = own;
                                            }
                                        }
                                    }
                                }
                            }
                            Err(e) => why = Some(e),
                        }
                    }
                    parts.push((values, owned, perm));
                }
                Ok(_) => why = Some(String::from("DOFs differ between ranks")),
                Err(e) => why = Some(e),
            }
        }
        if let Some(w) = why {
            notes.push(format!("solver {}: {var} ({base}.N): {w}", s.number));
            continue;
        }
        let d = dofs.unwrap_or(1);
        let opts = s.norm_opts(d);
        let norm = if !partitioning.is_dir() {
            notes.push(format!("solver {}: no partitioning.{np}: rank 0's serial norm", s.number));
            let (v, _, perm) = &parts[0];
            elmer_oracle::norm::compute_norm(v, Some(perm), d, perm.len(), &opts)
        } else if consistent {
            let owned: Vec<(Vec<f64>, Vec<bool>)> = parts.into_iter().map(|(v, o, _)| (v, o)).collect();
            elmer_oracle::norm::consistent_norm(&owned, d, &opts, order)
        } else {
            let values: Vec<Vec<f64>> = parts.into_iter().map(|(v, _, _)| v).collect();
            elmer_oracle::norm::parallel_norm(&values, d, &opts, order)
        };
        match norm {
            Ok(n) => refnorms.push((s.number, n)),
            Err(e) => notes.push(format!("solver {}: {var}: {e}", s.number)),
        }
    }
    refnorms
}

/// One test of `refnorm-validate`: (verdict, solvers validated, detail).
#[allow(clippy::too_many_arguments)]
fn refnorm_one(
    arm: &Arm,
    tests_dir: &Path,
    golden_root: &Path,
    out_root: &Path,
    test: &str,
    dir: &str,
    np: &str,
    mpi_order: elmer_oracle::norm::ReduceOrder,
) -> (String, usize, String) {
    let skip = |why: &str| (String::from("out-of-class"), 0, why.to_string());
    let Ok(np) = np.parse::<u32>() else {
        return skip("bad np");
    };
    let src = tests_dir.join(dir);
    let Ok(cmake) = std::fs::read_to_string(src.join("CMakeLists.txt")) else {
        return skip("no CMakeLists.txt");
    };
    let mut spec = parse_cmakelists(dir, &cmake);
    spec.startinfo_sif = read_startinfo(&src).ok().flatten();
    let Some(sif) = spec.startinfo_sif.clone() else {
        return skip("no ELMERSOLVER_STARTINFO");
    };
    let Ok(runtest) = read_runtest(&src) else {
        return skip("no runtest");
    };
    let Ok(plan) = plan_test(&spec, &runtest, np, arm) else {
        return skip("plan error");
    };
    // A custom `runtest.cmake` is replayed through `cmake -P` with the
    // arm's recorded registration, as `capture` does; the injection reaches
    // its sifs the same way.
    let registered = match &plan.custom {
        None => None,
        Some(_) => {
            let reg = read_registry_under(&arm.build_dir.join("fem").join("tests"))
                .ok()
                .and_then(|r| r.into_iter().find(|t| t.name == test && t.np == np));
            match reg {
                Some(r) => Some(r),
                None => return skip("custom runtest without a registration to replay"),
            }
        }
    };
    let solver_steps = plan.steps.iter().filter(|s| matches!(s, Step::Solver { .. })).count();
    if registered.is_none() && solver_steps != 1 {
        return skip("several solver runs");
    }
    let text: String = std::fs::read(src.join(&sif))
        .unwrap_or_default()
        .iter()
        .map(|&b| char::from(b))
        .collect();
    let golden_run = golden_root.join(test).join("run");
    let text = match elmer_oracle::refnorm::expand_includes(&text, &[&golden_run, &src]) {
        Ok(t) => t,
        Err(why) => return skip(&why),
    };
    let runs = elmer_oracle::runner::run_count(&text);
    let files = elmer_oracle::sigma::result_files(&golden_run, &src).unwrap_or_default();
    let mut notes = Vec::new();
    let mut non_nodal = false;
    let equations = elmer_oracle::refnorm::solver_equations(&text);
    // Each result file once: (relative name, parsed, preamble, nodes).
    let parsed: Vec<Parsed> = files
        .iter()
        .filter_map(|f| {
            let path = golden_run.join(f);
            let buf = std::fs::read(&path).ok()?;
            let r = elmer_oracle::result::parse_result_bytes(&buf, &path.display().to_string()).ok()?;
            let vars = elmer_oracle::result::preamble_vars(&buf);
            let nodes = elmer_oracle::result::number_of_nodes(&buf);
            Some((f.clone(), r, vars, nodes))
        })
        .collect();
    // A multi-RUN sif is checked after every RUN (CompareToReferenceSolution
    // in ElmerSolver's RUN loop): run j's solvers are those with a Reference
    // Norm in the sif up to that RUN, and its golden is the per-RUN file
    // `--per-run` writes (`case.run<j>.result`, the last RUN `case.result`).
    let mut refnorms = Vec::new();
    let mut run_refnorms: Vec<(usize, u32, f64)> = Vec::new();
    let mut solvers: Vec<elmer_oracle::refnorm::RefSolver> = Vec::new();
    if runs <= 1 {
        solvers = match elmer_oracle::refnorm::reference_solvers(&text) {
            Ok(s) if s.is_empty() => return skip("no Reference Norm in a Solver section"),
            Ok(s) => s,
            Err(why) => return skip(&why),
        };
        if np > 1 {
            refnorms = resolve_solvers_mpi(
                &solvers, &parsed, np as usize, &equations, &mut notes, mpi_order, &golden_run, &text,
            );
        } else {
            // A replayed script runs further sifs, each into its own
            // `<stem>.result`; the startinfo sif, the one the overrides go
            // into, writes `case.result`.
            let all: Vec<&Parsed> = parsed
                .iter()
                .filter(|p| {
                    registered.is_none()
                        || std::path::Path::new(&p.0).file_name().is_some_and(|n| n == "case.result")
                })
                .collect();
            refnorms = resolve_solvers(&solvers, &all, &equations, &mut notes, &mut non_nodal);
        }
    } else if np > 1 {
        return skip(&format!("MPI with {runs} RUN blocks"));
    } else {
        // Per-RUN files exist only where the injection names the output.
        if text.lines().any(|l| {
            let c = elmer_oracle::runner::strip_sif_comment(l);
            elmer_oracle::runner::keyword_of(c.split("::").last().unwrap_or("")).as_deref() == Some("output file")
        }) {
            return skip(&format!("{runs} RUN blocks; the sif names its own Output File"));
        }
        // With `Output Intervals = 0` nothing is saved during a RUN; the
        // file holds a field saved later (RotatingFlow's holds FlowSolve's
        // 1e-6 initial value), not the solution the RUN compares.
        if text.lines().any(|l| {
            let c = elmer_oracle::runner::strip_sif_comment(l);
            let key = elmer_oracle::runner::keyword_of(c.split("::").last().unwrap_or(""));
            key.as_deref().is_some_and(|k| k.starts_with("output intervals"))
                && c.split_once('=').is_some_and(|(_, v)| {
                    v.split_whitespace().filter(|w| !w.eq_ignore_ascii_case("integer")).all(|w| w == "0")
                })
        }) {
            return skip(&format!("{runs} RUN blocks with Output Intervals = 0 (no RUN's solution is saved)"));
        }
        for (j, prefix) in elmer_oracle::refnorm::run_prefixes(&text).iter().enumerate() {
            let j = j + 1;
            let run_solvers = match elmer_oracle::refnorm::reference_solvers(prefix) {
                Ok(s) => s,
                Err(why) => return skip(&format!("RUN {j}: {why}")),
            };
            let name = elmer_oracle::runner::per_run_result_name("case.result", j, runs);
            let mine: Vec<&Parsed> = parsed
                .iter()
                .filter(|p| std::path::Path::new(&p.0).file_name().is_some_and(|n| n.to_string_lossy() == name))
                .collect();
            if mine.is_empty() && !run_solvers.is_empty() {
                return skip(&format!("{runs} RUN blocks; no golden {name} (use the --per-run golden)"));
            }
            let mut run_notes = Vec::new();
            for (n, v) in resolve_solvers(&run_solvers, &mine, &equations, &mut run_notes, &mut non_nodal) {
                run_refnorms.push((j, n, v));
            }
            notes.extend(run_notes.into_iter().map(|m| format!("RUN {j}: {m}")));
            for s in run_solvers {
                if !solvers.iter().any(|t| t.number == s.number) {
                    solvers.push(s);
                }
            }
        }
        if solvers.is_empty() {
            return skip("no Reference Norm in a Solver section");
        }
    }
    if refnorms.is_empty() && run_refnorms.is_empty() {
        return skip(&notes.join("; "));
    }
    if text.to_lowercase().lines().any(|l| {
        let code = l.split('!').next().unwrap_or("");
        code.contains("adaptive mesh refinement") && code.contains("true")
    }) {
        return skip("adaptive remeshing (O1+O4 class)");
    }
    let run_dir = out_root.join(test);
    if run_dir.exists() {
        return skip("run dir exists; use a fresh out-root");
    }
    let run = |run_dir: &Path, refnorms: &[(u32, f64)], run_refnorms: &[(usize, u32, f64)]| {
        let mut staging = Staging {
            inject: Some(SifInject { o2: true, per_run: runs > 1, ..SifInject::default() }),
            mkdirs: Vec::new(),
            refnorms: refnorms.to_vec(),
            run_refnorms: run_refnorms.to_vec(),
            exit_verdict: false,
        };
        // Same in-memory-mesh retry as `capture`.
        let go = |staging: &Staging| match &registered {
            Some(reg) => run_replay(arm, reg, &src, run_dir, staging),
            None => run_plan(arm, &spec, &plan, &src, run_dir, staging),
        };
        let mut result = go(&staging);
        if result.as_ref().is_ok_and(|o| !o.passed) {
            let missing = missing_output_dirs(run_dir);
            if !missing.is_empty() && std::fs::remove_dir_all(run_dir).is_ok() {
                staging.mkdirs = missing;
                result = go(&staging);
            }
        }
        result
    };
    let solvers_n = refnorms.len() + run_refnorms.len();
    match run(&run_dir, &refnorms, &run_refnorms) {
        Ok(o) if o.passed => (String::from("bitwise"), solvers_n, notes.join("; ")),
        Ok(o) => {
            let failed: Vec<&String> = o.compare_lines.iter().filter(|l| l.contains("FAILED")).collect();
            // Known grader gaps, not capture failures: a solver that sets
            // its own norm, and edge, face and DG unknowns, follow
            // solver-specific norm paths.
            let failed_solvers: Vec<u32> = failed
                .iter()
                .filter_map(|l| l.split("Solver ").nth(1)?.split_whitespace().next()?.parse().ok())
                .collect();
            let self_norm = !failed_solvers.is_empty()
                && failed_solvers
                    .iter()
                    .all(|n| solvers.iter().any(|s| s.number == *n && s.self_norm()));
            notes.extend(failed.into_iter().cloned());
            if self_norm {
                // The other solvers can still be validated: rerun with the
                // self-norming ones left at their own sif Reference Norm.
                let rest: Vec<(u32, f64)> =
                    refnorms.iter().copied().filter(|(n, _)| !failed_solvers.contains(n)).collect();
                let run_rest: Vec<(usize, u32, f64)> =
                    run_refnorms.iter().copied().filter(|(_, n, _)| !failed_solvers.contains(n)).collect();
                for n in &failed_solvers {
                    notes.push(format!("solver {n}: sets its own norm, not validated"));
                }
                if !rest.is_empty() || !run_rest.is_empty() {
                    let retry = out_root.join(format!("{test}.without-self-norm"));
                    if let Ok(o) = run(&retry, &rest, &run_rest) {
                        if o.passed {
                            return (String::from("bitwise"), rest.len() + run_rest.len(), notes.join("; "));
                        }
                    }
                }
            }
            let verdict = if self_norm {
                "grader-gap-solver-norm"
            } else if non_nodal {
                "grader-gap-non-nodal"
            } else {
                "MISMATCH"
            };
            (String::from(verdict), solvers_n, notes.join("; "))
        }
        Err(e) => (String::from("ERROR"), solvers_n, e.to_string()),
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
            let runtest = read_runtest(&dir).unwrap_or_else(|_| usage());
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
            // `--label L` keeps tests carrying L; `--label '!L'` drops them.
            let keep = |t: &elmer_oracle::registry::RegisteredTest| match label.as_deref() {
                None => true,
                Some(l) => match l.strip_prefix('!') {
                    Some(not) => !t.has_label(not),
                    None => t.has_label(l),
                },
            };
            for t in reg.iter().filter(|t| keep(t)) {
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
        "o7-replay" => {
            // o7-replay <build-dir> <shim.so> <in-dir> <out-dir>: evaluate
            // the inputs of <in-dir>/{rules,info}.txt on the pristine arm
            // (o7.rs) and write the results in the same format.
            if args.len() < 6 {
                usage();
            }
            let (build, shim, input, out) =
                (Path::new(&args[2]), Path::new(&args[3]), Path::new(&args[4]), Path::new(&args[5]));
            let fail = |e: &dyn std::fmt::Display| -> ! {
                eprintln!("o7-replay: {e}");
                std::process::exit(1)
            };
            let read = |name: &str| std::fs::read_to_string(input.join(name)).unwrap_or_else(|e| fail(&e));
            let shim = elmer_oracle::o7::Shim::load(build, shim).unwrap_or_else(|e| fail(&e));
            let rules = elmer_oracle::o7::parse_rules(&read("rules.txt")).unwrap_or_else(|e| fail(&e));
            let info = elmer_oracle::o7::parse_info(&read("info.txt")).unwrap_or_else(|e| fail(&e));
            let got_rules = elmer_oracle::o7::replay_rules(&shim, &rules).unwrap_or_else(|e| fail(&e));
            let got_info = elmer_oracle::o7::replay_info(&shim, &info).unwrap_or_else(|e| fail(&e));
            std::fs::create_dir_all(out).unwrap_or_else(|e| fail(&e));
            std::fs::write(out.join("rules.txt"), elmer_oracle::o7::write_rules(&got_rules))
                .unwrap_or_else(|e| fail(&e));
            std::fs::write(out.join("info.txt"), elmer_oracle::o7::write_info(&got_info))
                .unwrap_or_else(|e| fail(&e));
            let dr = elmer_oracle::o7::diff_rules(&rules, &got_rules, false);
            let dr_w = elmer_oracle::o7::diff_rules(&rules, &got_rules, true);
            let di = elmer_oracle::o7::diff_info(&info, &got_info);
            let points: usize = info.iter().map(|e| e.pts.len()).sum();
            for d in dr.iter().chain(&di).take(40) {
                eprintln!("  {d}");
            }
            eprintln!(
                "rules={} rule-diffs={} (without quad w: {}) elements={} points={} info-diffs={}",
                rules.len(),
                dr.len(),
                dr_w.len(),
                info.len(),
                points,
                di.len()
            );
        }
        "refnorm-validate" => {
            // refnorm-validate <build-dir> <tests-dir> <golden-root> <out-root>:
            // §4.3 item 2 on the eligible class (see refnorm.rs).
            if args.len() < 6 {
                usage();
            }
            let build_dir = PathBuf::from(&args[2]);
            let tests_dir = PathBuf::from(&args[3]);
            let golden_root = PathBuf::from(&args[4]);
            let out_root = PathBuf::from(&args[5]);
            let registry = std::fs::read_to_string(golden_root.join("registry.tsv"))
                .unwrap_or_else(|_| usage());
            let arm = arm_from_build(build_dir, "release");
            let mut counts: std::collections::BTreeMap<String, usize> = Default::default();
            // Rows go to <out-root>/refnorm.tsv: tool steps inherit stdout.
            let mut rows = vec![String::from("test\tverdict\tsolvers\tdetail")];
            for line in registry.lines() {
                let cols: Vec<&str> = line.split('\t').collect();
                let (test, dir, np) = (cols[0], cols[1], cols.get(2).copied().unwrap_or("1"));
                use elmer_oracle::norm::ReduceOrder;
                let (mut verdict, mut solvers, mut detail) =
                    refnorm_one(&arm, &tests_dir, &golden_root, &out_root, test, dir, np, ReduceOrder::RecursiveDoubling);
                // An MPI run whose sums were reduced in another grouping:
                // try left to right once before calling it a mismatch.
                if np != "1" && verdict != "bitwise" && verdict != "out-of-class" {
                    let (v, n, d) = refnorm_one(
                        &arm, &tests_dir, &golden_root, &out_root.join("sequential"), test, dir, np, ReduceOrder::Sequential,
                    );
                    if v == "bitwise" {
                        (verdict, solvers, detail) = (v, n, format!("sequential reduction; {d}"));
                    } else {
                        detail = format!("{detail}; sequential reduction: {v}");
                    }
                }
                let key = verdict.split(':').next().unwrap_or(&verdict).to_string();
                *counts.entry(key).or_default() += 1;
                rows.push(format!("{test}\t{verdict}\t{solvers}\t{detail}"));
            }
            std::fs::create_dir_all(&out_root).unwrap_or_else(|_| usage());
            std::fs::write(out_root.join("refnorm.tsv"), rows.join("\n") + "\n")
                .unwrap_or_else(|_| usage());
            let summary: Vec<String> = counts.iter().map(|(k, v)| format!("{k}={v}")).collect();
            eprintln!("{}", summary.join(" "));
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
            // --all also lists caught mutations, for before/after tables.
            let all = args.iter().any(|a| a == "--all");
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
                            Caught => {
                                caught += 1;
                                if all {
                                    println!("caught\t{test}\t{rel}\t{m}\t{}\t{}", t.target, t.detail);
                                }
                            }
                            Escaped => {
                                escaped += 1;
                                println!("ESCAPED\t{test}\t{rel}\t{m}\t{}", t.target);
                            }
                            BelowFloor => {
                                below += 1;
                                println!("below-floor\t{test}\t{rel}\t{m}\t{}\t{}", t.target, t.detail);
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
