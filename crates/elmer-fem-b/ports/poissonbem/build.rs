//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only the reference `.so`'s
//! dynamic symbols stay global (see `poissonbem.map`, enforced by the
//! proof script): `poissonbemsolver_`, `matvec_`, `precond_`, and the
//! `__globmat_MOD_matrix` BSS descriptor owned by the accessor's own
//! `globmat` module mirror.
//!
//! Flags mirror the Elmer `flags.make` for the `PoissonBEM` target:
//! `-fvisibility=hidden -g -O0` plus `-fPIC` (the accessor references
//! external Elmer module data), `-DCONTIG=` (the kept `FullIterSolver`
//! transcription uses the legacy `POINTER CONTIG` spelling), and
//! `-fallow-argument-mismatch` (the kept `SolveFull` transcription
//! passes the rank-2 `Matrix` to the rank-1 `a(n*n)` LAPACK dummy,
//! exactly like the original — sequence association the Elmer build
//! itself allows via that flag). The `huti_fdefs.h` include dir is
//! resolved from `ELMER_MODDIR` (`<work>/src/fhutiter/include`) with
//! `ELMER_FORK` fallbacks. Emits `cargo:rustc-link-arg` (the
//! `rust-link-arg` spelling silently drops the object).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");

    // Resolve the dir carrying huti_fdefs.h: <work>/src/fhutiter/include
    // (what the Elmer build itself used), else the Elmer checkout.
    let mut candidates: Vec<PathBuf> = Vec::new();
    let modpath = PathBuf::from(&moddir);
    if let Some(work) = modpath.parent().and_then(|p| p.parent()) {
        candidates.push(work.join("src/fhutiter/include"));
    }
    if let Ok(fork) = env::var("ELMER_FORK") {
        candidates.push(PathBuf::from(&fork).join("fhutiter/include"));
        candidates.push(PathBuf::from(&fork).join("fhutiter/src"));
    }
    let hutidir = candidates
        .iter()
        .find(|p| p.join("huti_fdefs.h").is_file())
        .unwrap_or_else(|| {
            panic!(
                "huti_fdefs.h not found; tried {:?} (set ELMER_FORK)",
                candidates
            )
        });

    let obj = out.join("poissonbem_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest.join("accessors/poissonbem_acc.F90").to_str().unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fvisibility=hidden",
            "-fallow-argument-mismatch",
            "-g",
            "-O0",
            "-fPIC",
            "-DCONTIG=",
            &format!("-I{moddir}"),
            &format!("-I{}", hutidir.display()),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("poissonbem.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/poissonbem_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
    println!("cargo:rerun-if-env-changed=ELMER_FORK");
}
