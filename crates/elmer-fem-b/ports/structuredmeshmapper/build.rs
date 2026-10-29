//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only
//! `structuredmeshmapper_` and `structuredmeshmapper_init_` stay dynamic
//! symbols (see `structuredmeshmapper.map`, enforced by the proof script).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let obj = out.join("structuredmeshmapper_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/structuredmeshmapper_acc.F90")
                .to_str()
                .unwrap(),
            "-o",
            obj.to_str().unwrap(),
            "-J",
            out.to_str().unwrap(),
            // -fPIC: the accessor references external module data (the
            // Messages module's Message string, DetectExtrudedStructure's
            // mesh tables), so it needs -fPIC like the Elmer objects.
            "-fPIC",
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
        manifest.join("structuredmeshmapper.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/structuredmeshmapper_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
