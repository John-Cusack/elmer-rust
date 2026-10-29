//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only `meshsolver_` stays a
//! dynamic symbol (see `nonphysicalmeshsolve.map`, enforced by the proof
//! script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("nonphysicalmeshsolve_acc.o");
    // -fPIC: the accessor references external module data (Elmer solver
    // state, the Messages `Message` scratch for progress output), so the
    // object must use GOT-relative addressing inside the cdylib.
    // -DCONTIG=: the unit spells pointer attributes `POINTER CONTIG`
    // (legacy, no comma); the Elmer build compiles with -DCONTIG="" so
    // the accessor must match, or the .mod-declared pointer shapes clash.
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/nonphysicalmeshsolve_acc.F90")
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
            "-DCONTIG=",
            &format!("-I{moddir}"),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("nonphysicalmeshsolve.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/nonphysicalmeshsolve_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
