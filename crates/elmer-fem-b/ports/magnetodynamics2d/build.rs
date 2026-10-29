//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only the seven entry symbols
//! stay dynamic (see `magnetodynamics2d.map`, enforced by the proof
//! script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("magnetodynamics2d_acc.o");
    // -fPIC: the accessor references external module data (Message via
    // Messages, model/solver state via DefUtils scope, circuit tables via
    // CircuitUtils, hysteresis state via ZirkaUtils), so the object must
    // be position-independent for the shared link. FP codegen is
    // unaffected. -DCONTIG="": BSolver spells a pointer `POINTER CONTIG`
    // and the reference flags.make resolves CONTIG via the preprocessor,
    // so the accessor needs the same. Rank-1-view-of-rank-2
    // transcriptions (BulkAssembly scratch passed as explicit-shape
    // columns) need -fallow-argument-mismatch, per Elmer flags.make.
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/magnetodynamics2d_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fvisibility=hidden",
            "-fPIC",
            "-g",
            "-O0",
            "-fallow-argument-mismatch",
            "-DCONTIG=",
            &format!("-I{moddir}"),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("magnetodynamics2d.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/magnetodynamics2d_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
