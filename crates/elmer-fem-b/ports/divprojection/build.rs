//! Build script: links the cdylib with a version script so only
//! `divprojection_` stays a dynamic symbol (see `divprojection.map`,
//! enforced by the proof script). The original body issues no Elmer
//! calls, so this port compiles no Fortran accessors.
use std::env;
use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    println!(
        "cargo:rustc-link-arg=-Wl,--version-script,{}",
        manifest.join("divprojection.map").display()
    );
}
