//! Build script: compiles the generated-style Fortran accessor file
//! against the built `.mod` files of the target ice-enabled configure
//! and links it into the cdylib with hidden visibility. Only the 60
//! reference dynamic symbols stay global (see `usf.map`, enforced by
//! the proof script).
//!
//! Like the reference `ElmerIceUSF.so` (which carries
//! `NEEDED ElmerIceUtils.so + libelmersolver.so` with RUNPATH into the
//! build tree), the port links both built-tree libraries; at proof time
//! those resolve to the reference tree, so both arms delegate to the
//! identical routines (ProjUtils projection math, all libelmersolver
//! calls behind the accessors).
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let build =
        env::var("ELMER_BUILD").expect("ELMER_BUILD must point at the built ice-enabled tree");
    let moddir = PathBuf::from(&build).join("fmodules");
    let solversdir = PathBuf::from(&build).join("fem/src");
    let utilsdir = PathBuf::from(&build).join("elmerice/Utils");
    assert!(
        moddir.join("types.mod").exists(),
        "no .mod files under {}",
        moddir.display()
    );
    let obj = out.join("usf_acc.o");
    let status = Command::new("gfortran")
        .args([
            "-c",
            manifest
                .join("accessors/usf_acc.F90")
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
            &format!("-I{}", moddir.display()),
        ])
        .status()
        .expect("gfortran must be installed");
    assert!(status.success(), "accessor compile failed");
    println!("cargo:rustc-link-arg={}", obj.display());
    println!("cargo:rustc-link-arg=-L{}", solversdir.display());
    println!("cargo:rustc-link-arg=-L{}", utilsdir.display());
    // Like the reference (NEEDED libgfortran.so.5 + liblapack.so.3):
    // the accessors use Fortran I/O and call LAPACK's DGEEV, so the
    // runtime + lapack ride in the link explicitly instead of via host
    // global scope. This also lets the cargo test harness link.
    println!("cargo:rustc-link-arg=-Wl,--no-as-needed,-l:libgfortran.so.5,-l:liblapack.so.3,--as-needed");
    println!("cargo:rustc-link-arg=-Wl,--no-as-needed,-l:libelmersolver.so,-l:ElmerIceUtils.so,--as-needed");
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,{}",
        solversdir.display()
    );
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", utilsdir.display());
    // `usf.map` keeps the reference 60-symbol surface plus ONE explicit
    // local: `usf_metric_core` (Rust core behind the Fortran entry shim
    // for `glaciermeshmetricaniso_`). The explicit entry is load-bearing:
    // the `local: *` wildcard alone does NOT hide Rust `#[no_mangle]` /
    // `#[export_name]` symbols with this toolchain (verified: wildcard-only
    // leaves it global and breaks the `nm -D` contract; explicit keeps it
    // local). (The map file itself takes no comments — the linker
    // rejects them.)
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("usf.map").display()
    );
    println!("cargo:rerun-if-changed=accessors/usf_acc.F90");
    println!("cargo:rerun-if-changed=usf.map");
    println!("cargo:rerun-if-env-changed=ELMER_BUILD");
}
