//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only `phasechangesolve_`,
//! `phasechangesolve_init_` and `meltingheat_` stay dynamic symbols (see
//! `phasechangesolve.map`, enforced by the proof script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("phasechangesolve_acc.o");
    // Template flags (-fvisibility=hidden -g -O0) plus -fPIC (the
    // accessor references external module data: ParMatrix, Message)
    // plus -fallow-argument-mismatch (matches the Elmer flags.make for
    // this configure).
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/phasechangesolve_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fvisibility=hidden",
            "-g",
            "-O0",
            "-fPIC",
            "-fallow-argument-mismatch",
            &format!("-I{moddir}"),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("phasechangesolve.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/phasechangesolve_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
