//! Build script: compiles the transcribed `vecutils` module (default
//! visibility, so the four `__vecutils_MOD_*` procedures stay dynamic and
//! are globalized by the version script) and the generated-style Fortran
//! accessor file (hidden visibility, as in the poisson port) against the
//! built `.mod` files of the target configure, and links both into the
//! cdylib with a version script. Only the six reference entry symbols
//! stay dynamic (see `coordinatetransform.map`, enforced by the proof
//! script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let vec_obj = out.join("vecutils_rep.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/vecutils_rep.F90")
                .to_str()
                .unwrap(),
            "-o",
            vec_obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-g",
            "-O0",
            &format!("-I{moddir}"),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "vecutils transcription compile failed");
    let acc_obj = out.join("coordinatetransform_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/coordinatetransform_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            acc_obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fvisibility=hidden",
            "-g",
            "-O0",
            &format!("-I{moddir}"),
            "-I",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", vec_obj.display());
    println!("cargo:rustc-link-arg={}", acc_obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("coordinatetransform.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/vecutils_rep.F90");
    println!("cargo:rerun-if-changed=accessors/coordinatetransform_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
