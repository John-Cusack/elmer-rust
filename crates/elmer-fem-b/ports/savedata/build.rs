//! Build script: compiles the generated-style Fortran accessor files
//! against the built `.mod` files of the target configure and links them
//! into the cdylib with hidden visibility. Only the ten SaveData entries
//! stay dynamic symbols (see `savedata.map`, enforced by proof).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    for stem in [
        "sd_bvals_acc",
        "sd_depend_acc",
        "sd_line_acc",
        "sd_line_rest",
        "sd_mats_acc",
        "sd_proj_acc",
        "sd_scalars_acc",
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
                &format!("-I{}", out.display()),
            ])
            .status()
            .expect("gfortran must be installed");
        assert!(status.success(), "accessor compile failed for {stem}");
        println!("cargo:rustc-link-arg={}", obj.display());
    }
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("savedata.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/sd_bvals_acc.F90");
    println!("cargo:rerun-if-changed=accessors/sd_depend_acc.F90");
    println!("cargo:rerun-if-changed=accessors/sd_line_acc.F90");
    println!("cargo:rerun-if-changed=accessors/sd_line_rest.F90");
    println!("cargo:rerun-if-changed=accessors/sd_mats_acc.F90");
    println!("cargo:rerun-if-changed=accessors/sd_proj_acc.F90");
    println!("cargo:rerun-if-changed=accessors/sd_scalars_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
}
