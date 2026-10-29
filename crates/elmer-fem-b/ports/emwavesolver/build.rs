//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only the five entry symbols
//! stay dynamic symbols (see `emwavesolver.map`, enforced by the proof
//! script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("emwavesolver_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/emwavesolver_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            "-fPIC",
            "-fallow-argument-mismatch",
            // Matches the reference module build (elmerf90-helper.h passes
            // -DCONTIG=, so the CONTIGUOUS-conditional CONTIG macro expands
            // empty; the accessor transcribes the original CONTIG attributes
            // verbatim).
            "-DCONTIG=",
            "-fvisibility=hidden",
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
        manifest.join("emwavesolver.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/emwavesolver_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
