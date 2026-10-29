//! Build script: compiles the generated-style Fortran accessor files
//! against the built `.mod` files of the target configure and links them
//! into the cdylib with hidden visibility. Only the four solver entries
//! stay dynamic symbols (see `consistentsplitting.map`, enforced by proof).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    for stem in [
        "cs_comp_acc",
        "cs_pres_acc",
        "cs_velo_acc",
        "cs_opt_acc",
    ] {
        let obj = out.join(format!("{stem}.o"));
        let status = Command::new("gfortran")
            .args([
                "-c",
                manifest
                    .join(format!("accessors/{stem}.F90"))
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
                &format!("-I{moddir}"),
            ])
            .status()
            .expect("gfortran must be installed");
        assert!(status.success(), "accessor compile failed for {stem}");
        println!("cargo:rustc-link-arg={}", obj.display());
    }
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("consistentsplitting.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/cs_comp_acc.F90");
    println!("cargo:rerun-if-changed=accessors/cs_pres_acc.F90");
    println!("cargo:rerun-if-changed=accessors/cs_velo_acc.F90");
    println!("cargo:rerun-if-changed=accessors/cs_opt_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
