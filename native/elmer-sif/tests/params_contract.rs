//! The `elmer-core` contract's `Params` option (interfaces.md row 25):
//! `&Section` converts into a `ParamMap` implementing [`Params`].

use elmer_core::params::{ParamMap, Params};

use std::io::Write;
fn parse_body(text: &str) -> elmer_sif::Model {
    let path = std::env::temp_dir().join(format!("contract-{}.sif", std::process::id()));
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
