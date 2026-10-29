//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only `freesurfacereduced_`
//! stays a dynamic symbol (see `freesurfacereduced.map`, enforced by the
//! proof script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("freesurfacereduced_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/freesurfacereduced_acc.F90")
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
        manifest.join("freesurfacereduced.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/freesurfacereduced_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
