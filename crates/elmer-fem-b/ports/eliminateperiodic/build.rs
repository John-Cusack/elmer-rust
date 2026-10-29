//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target configure and links it
//! into the cdylib with hidden visibility. Only `eliminateperiodic_`
//! stays a dynamic symbol (see `eliminateperiodic.map`, enforced by the
//! proof script).
//!
//! Like the reference `EliminatePeriodic.so` (which carries
//! `NEEDED EliminateDirichlet.so` for its no-periodic delegation and its
//! reduced-system solve), the port links the built-tree
//! `EliminateDirichlet.so`; at proof time that resolves to the reference
//! module, so both arms delegate to the identical routine.
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let moddir = env::var("ELMER_MODDIR")
        .expect("ELMER_MODDIR must point at the built .mod dir (fmodules)");
    let modulesdir = env::var("ELMER_MODULESDIR").expect(
        "ELMER_MODULESDIR must point at the built module dir (fem/src/modules)",
    );
    let obj = out.join("eliminateperiodic_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/eliminateperiodic_acc.F90")
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
    println!("cargo:rustc-link-arg=-L{modulesdir}");
    println!("cargo:rustc-link-arg=-Wl,--no-as-needed,-l:EliminateDirichlet.so,--as-needed");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{modulesdir}");
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("eliminateperiodic.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/eliminateperiodic_acc.F90");
    println!("cargo:rerun-if-env-changed=ELMER_MODDIR");
    println!("cargo:rerun-if-env-changed=ELMER_MODULESDIR");
}
