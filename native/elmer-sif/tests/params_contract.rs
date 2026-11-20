//! The `elmer-core` contract's `Params` option (interfaces.md row 25):
//! `&Section` converts into a `ParamMap` implementing [`Params`].

use elmer_core::params::{ParamMap, Params};

use std::io::Write;
fn parse_body(text: &str) -> elmer_sif::Model {
    // One file per call: the tests of this file run in parallel.
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let k = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("contract-{}-{k}.sif", std::process::id()));
    let path = path.to_str().unwrap().to_string();
    let mut f = std::fs::File::create(&path).unwrap();
    f.write_all(text.as_bytes()).unwrap();
    drop(f);
    let parsed = elmer_sif::parse_file(&path).unwrap();
    std::fs::remove_file(&path).ok();
    assert_eq!(parsed.sif.runs.len(), 1);
    parsed.sif.runs.into_iter().next().unwrap()
}

#[test]
fn section_converts_to_params() {
    let mut run = parse_body(
        "Header\n  Mesh DB \".\" \"m\"\nEnd\nBody 1\n\
         \x20 Name = \"Water\"\n\
         \x20 Density = Real 998.0\n\
         \x20 Viscosity = Real 1.0e-3\n\
         \x20 Equation = Integer 1\n\
         \x20 Active = Logical True\n\
         \x20 Data File = File \"data.dat\"\n\
         End\n",
    );
    elmer_sif::complete::complete_run(&mut run, "contract.sif").unwrap();
    assert_eq!(run.bodies.len(), 1);
    let map = ParamMap::from(&run.bodies[0]);
    // Keys are case-insensitive; scalars read through the typed getters.
    assert_eq!(map.real("DENSITY"), Some(998.0));
    assert_eq!(map.real("viscosity"), Some(1.0e-3));
    assert_eq!(map.integer("equation"), Some(1));
    assert_eq!(map.logical("active"), Some(true));
    // `Str` (quoted text) and `File` both read through `string`; Elmer's
    // `ListAddString` lowercases `Str` at storage, `File` keeps case.
    assert_eq!(map.string("name"), Some("water"));
    assert_eq!(map.string("data file"), Some("data.dat"));
    // A table value is a Dependent, not a constant: typed getters say None.
    let mut run = parse_body(
        "Header\n  Mesh DB \".\" \"m\"\nEnd\nMaterial 1\n\
         \x20 Viscosity = Variable Temperature\n\
         \x20\x20 Real\n\
         \x20\x20\x20 273.0 1.8e-3\n\
         \x20\x20\x20 373.0 0.3e-3\n\
         \x20\x20 End\n\
         End\n",
    );
    elmer_sif::complete::complete_run(&mut run, "contract.sif").unwrap();
    let map = ParamMap::from(&run.materials[0]);
    assert_eq!(map.real("viscosity"), None);
    assert!(map.get("viscosity").is_some());
}

/// `AddVtuOutputSolverHack` / `AddSaveScalarsHack` find an existing solver
/// by its Procedure as written (`FindSolverByProcName`, case-sensitive
/// INDEX): a VTU ResultOutputSolver and a SaveScalars given as `"Lib"
/// "Entry"` stop the hacks; a non-VTU ResultOutputSolver does not.
#[test]
fn post_file_and_scalars_file_hacks_see_existing_solvers() {
    let head = "Header\n  Mesh DB \".\" \"m\"\nEnd\nSimulation\n  Post File = \"case.vtu\"\n  Scalars File = \"s.dat\"\nEnd\n";
    let solvers = |fmt: &str| {
        format!(
            "{head}Solver 1\n  Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n  {fmt}\nEnd\n\
             Solver 2\n  Procedure = \"SaveData\" \"SaveScalars\"\nEnd\n"
        )
    };
    let mut run = parse_body(&solvers("Vtu Format = Logical True"));
    elmer_sif::complete::complete_run(&mut run, "contract.sif").unwrap();
    assert_eq!(run.solvers.len(), 2);
    let mut run = parse_body(&solvers("Gmsh Format = Logical True"));
    elmer_sif::complete::complete_run(&mut run, "contract.sif").unwrap();
    assert_eq!(run.solvers.len(), 3, "a non-VTU ResultOutputSolver gets the VTU one added");
}
