//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only `filmflowsolver_`,
//! `filmflowsolver_init_`, and `filmflowsolver_init0_` stay dynamic
//! symbols (see `filmflowsolver.map`, enforced by the proof script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("filmflowsolver_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/filmflowsolver_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fvisibility=hidden",
            "-fPIC",
            "-fallow-argument-mismatch",
            "-g",
            "-O0",
            &format!("-I{moddir}"),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("filmflowsolver.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/filmflowsolver_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
