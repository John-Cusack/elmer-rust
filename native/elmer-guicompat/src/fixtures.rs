//! The seven baseline ElmerGUI input fixtures, one per shipped EDF.
//!
//! Each fixture is a complete case directory as the GUI leaves it before it
//! starts the solver: `case.sif` (untyped `key = value` lines in the GUI's
//! style), `ELMERSOLVER_STARTINFO` naming the SIF, and an Elmer-format mesh.
//! The solver blocks take `Procedure` and `Variable` from the EDF defaults
//! exactly as the GUI would; the remaining keywords are chosen per physics
//! and are looked up in the EDF (`edf_keys`) or flagged as general-tab
//! keywords (`general_keys`). The fixtures describe what the GUI sends. They
//! do not mean any physics is supported natively.

use std::collections::BTreeMap;
use std::path::Path;

use crate::edf::{parse_edf, Edf, EdfError};
use crate::mesh::SquareMesh;
use crate::sif::{gui_variable, SifWriter};

/// What a native implementation is expected to do with the case.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    /// Scalar heat: the first native physics target.
    Heat,
    /// Linear elasticity: the second.
    Elasticity,
    /// Heat plus `ResultOutputSolve` writing VTU.
    Output,
    /// Physics with no native implementation planned in this burst. A
    /// native run must refuse (exit 78); it must never fake a solve.
    Unsupported,
}

impl Class {
    pub fn as_str(self) -> &'static str {
        match self {
            Class::Heat => "heat",
            Class::Elasticity => "elasticity",
            Class::Output => "output",
            Class::Unsupported => "unsupported",
        }
    }
    pub fn native_target(self) -> bool {
        self != Class::Unsupported
    }
}

/// `(id, edf file)` of the seven shipped dialects.
pub const EDF_FILES: [(&str, &str); 7] = [
    ("heatequation", "heatequation.xml"),
    ("linearelasticity", "linearelasticity.xml"),
    ("navier-stokes", "navier-stokes.xml"),
    ("electrostatics", "electrostatics.xml"),
    ("helmholtz", "helmholtz.xml"),
    ("meshdeform", "meshdeform.xml"),
    ("resultoutput", "resultoutput.xml"),
];

pub struct EdfSet(pub BTreeMap<String, Edf>);

impl EdfSet {
    pub fn load(dir: &Path) -> Result<EdfSet, String> {
        let mut map = BTreeMap::new();
        for (id, file) in EDF_FILES {
            let p = dir.join(file);
            let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
            let edf = parse_edf(&text).map_err(|e: EdfError| format!("{}: {e}", p.display()))?;
            map.insert(id.to_string(), edf);
        }
        Ok(EdfSet(map))
    }
    pub fn get(&self, id: &str) -> &Edf {
        &self.0[id]
    }
}

type Kv = (&'static str, &'static str);

struct Spec {
    class: Class,
    dim: i64,
    /// Extra Solver keywords after the common numerics.
    solver: Vec<Kv>,
    direct: bool,
    stabilize: bool,
    equation: Vec<Kv>,
    material: Vec<Kv>,
    body_force: Vec<Kv>,
    /// `(target boundary, keywords)`.
    boundaries: Vec<(u32, Vec<Kv>)>,
}

fn spec(id: &str) -> Spec {
    match id {
        "heatequation" => Spec {
            class: Class::Heat,
            dim: 2,
            solver: vec![],
            direct: false,
            stabilize: true,
            equation: vec![],
            material: vec![
                ("Heat Conductivity", "1.0"),
                ("Density", "1.0"),
                ("Heat Capacity", "1.0"),
            ],
            body_force: vec![("Heat Source", "1.0")],
            boundaries: vec![
                (1, vec![("Temperature", "0.0")]),
                (3, vec![("Temperature", "1.0")]),
            ],
        },
        "linearelasticity" => Spec {
            class: Class::Elasticity,
            dim: 2,
            solver: vec![("Variable", "-dofs 2 Displacement")],
            direct: true,
            stabilize: false,
            equation: vec![("Plane Stress", "True")],
            material: vec![
                ("Youngs modulus", "1.0e5"),
                ("Poisson ratio", "0.3"),
                ("Density", "1.0"),
            ],
            body_force: vec![],
            boundaries: vec![
                (
                    4,
                    vec![("Displacement 1", "0.0"), ("Displacement 2", "0.0")],
                ),
                (2, vec![("Force 1", "100.0")]),
            ],
        },
        "navier-stokes" => Spec {
            class: Class::Unsupported,
            dim: 2,
            solver: vec![],
            direct: true,
            stabilize: true,
            equation: vec![("NS Convect", "False")],
            material: vec![("Density", "1.0"), ("Viscosity", "1.0")],
            body_force: vec![],
            boundaries: vec![
                (4, vec![("Velocity 1", "1.0"), ("Velocity 2", "0.0")]),
                (1, vec![("Noslip wall BC", "True")]),
                (3, vec![("Noslip wall BC", "True")]),
                (2, vec![("External Pressure", "0.0")]),
            ],
        },
        "electrostatics" => Spec {
            class: Class::Unsupported,
            dim: 2,
            solver: vec![],
            direct: true,
            stabilize: false,
            equation: vec![],
            material: vec![("Relative Permittivity", "1.0")],
            body_force: vec![],
            boundaries: vec![
                (4, vec![("Potential", "0.0")]),
                (2, vec![("Potential", "1.0")]),
            ],
        },
        "helmholtz" => Spec {
            class: Class::Unsupported,
            dim: 2,
            solver: vec![],
            direct: true,
            stabilize: false,
            equation: vec![("Angular Frequency", "1000.0")],
            material: vec![("Density", "1.2"), ("Sound speed", "343.0")],
            body_force: vec![],
            boundaries: vec![
                (
                    4,
                    vec![("Pressure Wave 1", "1.0"), ("Pressure Wave 2", "0.0")],
                ),
                (
                    2,
                    vec![("Pressure Wave 1", "0.0"), ("Pressure Wave 2", "0.0")],
                ),
            ],
        },
        "meshdeform" => Spec {
            class: Class::Unsupported,
            dim: 2,
            solver: vec![],
            direct: true,
            stabilize: false,
            equation: vec![],
            material: vec![],
            body_force: vec![],
            boundaries: vec![
                (1, vec![("Mesh Update 1", "0.0"), ("Mesh Update 2", "0.0")]),
                (3, vec![("Mesh Update 1", "0.0"), ("Mesh Update 2", "0.1")]),
            ],
        },
        // Heat with a VTU writer: the output-loading fixture.
        "resultoutput" => Spec {
            class: Class::Output,
            ..spec("heatequation")
        },
        other => panic!("no fixture spec for `{other}`"),
    }
}

#[derive(Clone, Debug)]
pub struct Fixture {
    pub id: String,
    pub class: Class,
    /// Equation title the GUI writes (`<Name>` of the EDF's PDE).
    pub title: String,
    /// Relative path to file contents.
    pub files: Vec<(String, String)>,
    /// Files a successful run leaves behind.
    pub expected_outputs: Vec<String>,
    /// `(EDF section, SIF section, keyword)` for keywords found in the EDF.
    pub edf_keys: Vec<(String, String, String)>,
    /// `(SIF section, keyword)` for keywords outside the physics EDF
    /// (the GUI's general tabs).
    pub general_keys: Vec<(String, String)>,
}

struct Builder<'a> {
    edf: &'a Edf,
    edf_keys: Vec<(String, String, String)>,
    general_keys: Vec<(String, String)>,
}

impl Builder<'_> {
    /// Write `key = value`, spelling the key as the EDF does when it knows it.
    fn put(&mut self, w: &mut SifWriter, edf_section: &str, sif: &str, key: &str, value: &str) {
        match self.edf.param(edf_section, key) {
            Some(p) => {
                self.edf_keys
                    .push((edf_section.into(), sif.into(), p.key().to_string()));
                w.kv(p.key(), value);
            }
            None => {
                self.general_keys.push((sif.into(), key.into()));
                w.kv(key, value);
            }
        }
    }
}

/// GUI numerics tabs for a matrix solver (`parse*Tab` in `sifgenerator.cpp`).
fn solver_numerics(w: &mut SifWriter, direct: bool, stabilize: bool) {
    w.kv("Exec Solver", "Always");
    w.kv("Stabilize", if stabilize { "True" } else { "False" });
    w.kv("Optimize Bandwidth", "True");
    w.kv("Steady State Convergence Tolerance", "1.0e-5");
    w.kv("Nonlinear System Convergence Tolerance", "1.0e-7");
    w.kv("Nonlinear System Max Iterations", "20");
    w.kv("Nonlinear System Newton After Iterations", "3");
    w.kv("Nonlinear System Newton After Tolerance", "1.0e-3");
    w.kv("Nonlinear System Relaxation Factor", "1.0");
    if direct {
        w.kv("Linear System Solver", "Direct");
        w.kv("Linear System Direct Method", "Banded");
    } else {
        w.kv("Linear System Solver", "Iterative");
        w.kv("Linear System Iterative Method", "BiCGStab");
        w.kv("Linear System Max Iterations", "500");
        w.kv("Linear System Convergence Tolerance", "1.0e-10");
        w.kv("BiCGstabl polynomial degree", "2");
        w.kv("Linear System Preconditioning", "ILU0");
        w.kv("Linear System ILUT Tolerance", "1.0e-3");
        w.kv("Linear System Abort Not Converged", "False");
        w.kv("Linear System Residual Output", "10");
    }
}

pub fn build(id: &str, set: &EdfSet) -> Fixture {
    let sp = spec(id);
    let edf = set.get(if id == "resultoutput" {
        "heatequation"
    } else {
        id
    });
    let mut b = Builder {
        edf,
        edf_keys: Vec::new(),
        general_keys: Vec::new(),
    };
    let mut w = SifWriter::new();

    w.block("Header")
        .raw("CHECK KEYWORDS Warn")
        .raw("Mesh DB \".\" \".\"")
        .raw("Include Path \"\"")
        .raw("Results Directory \"\"");
    w.block("Simulation")
        .kv("Max Output Level", "5")
        .kv("Coordinate System", "Cartesian")
        .kv("Coordinate Mapping(3)", "1 2 3")
        .kv("Simulation Type", "Steady state")
        .kv("Steady State Max Iterations", "1")
        .kv("Output Intervals(1)", "1")
        .kv("Solver Input File", "case.sif");
    if sp.class != Class::Output {
        w.kv("Post File", "case.ep");
    }
    w.block("Constants")
        .kv("Gravity(4)", "0 -1 0 9.82")
        .kv("Stefan Boltzmann", "5.670374419e-08")
        .kv("Permittivity of Vacuum", "8.8541878128e-12")
        .kv("Permeability of Vacuum", "1.25663706212e-6")
        .kv("Boltzmann Constant", "1.380649e-23")
        .kv("Unit Charge", "1.6021766208e-19");

    let n_solvers = if sp.class == Class::Output { 2 } else { 1 };
    let has_bf = !sp.body_force.is_empty();
    w.block("Body 1")
        .kv("Target Bodies(1)", "1")
        .kv("Name", "\"Body 1\"")
        .kv("Equation", "1")
        .kv("Material", "1");
    if has_bf {
        w.kv("Body Force", "1");
    }

    w.block("Solver 1")
        .kv("Equation", &format!("\"{}\"", edf.title));
    let procedure = edf.solver_default("Procedure").unwrap_or("").to_string();
    b.put(&mut w, "Solver", "Solver", "Procedure", &procedure);
    if let Some(v) = edf.solver_default("Variable") {
        let v = gui_variable(v, sp.dim);
        w.kv("Variable", &v);
    }
    for (k, v) in &sp.solver {
        b.put(&mut w, "Solver", "Solver", k, v);
    }
    solver_numerics(&mut w, sp.direct, sp.stabilize);

    if sp.class == Class::Output {
        let out = set.get("resultoutput");
        w.block("Solver 2")
            .kv("Equation", &format!("\"{}\"", out.title));
        let mut ob = Builder {
            edf: out,
            edf_keys: Vec::new(),
            general_keys: Vec::new(),
        };
        let p = out.solver_default("Procedure").unwrap_or("").to_string();
        ob.put(&mut w, "Solver", "Solver", "Procedure", &p);
        // Checkboxes left at their EDF default (No Matrix Equation, Binary Output: both True) are not written by the GUI.
        w.kv("Exec Solver", "After Simulation");
        for (k, v) in [
            ("Output Format", "Vtu"),
            ("Output File Name", "case"),
            ("Scalar Field 1", "Temperature"),
        ] {
            ob.put(&mut w, "Solver", "Solver", k, v);
        }
        b.edf_keys.extend(ob.edf_keys);
        b.general_keys.extend(ob.general_keys);
    }

    let eq_solvers = (1..=n_solvers).map(|i| i.to_string()).collect::<Vec<_>>();
    w.block("Equation 1").kv("Name", "\"Equation 1\"").kv(
        &format!("Active Solvers({n_solvers})"),
        &eq_solvers.join(" "),
    );
    for (k, v) in &sp.equation {
        b.put(&mut w, "Equation", "Equation", k, v);
    }

    w.block("Material 1").kv("Name", "\"Material 1\"");
    for (k, v) in &sp.material {
        b.put(&mut w, "Material", "Material", k, v);
    }
    if has_bf {
        w.block("Body Force 1").kv("Name", "\"BodyForce 1\"");
        for (k, v) in &sp.body_force {
            b.put(&mut w, "BodyForce", "Body Force", k, v);
        }
    }
    for (i, (target, kvs)) in sp.boundaries.iter().enumerate() {
        w.block(&format!("Boundary Condition {}", i + 1))
            .kv("Target Boundaries(1)", &target.to_string())
            .kv("Name", &format!("\"BoundaryCondition {}\"", i + 1));
        for (k, v) in kvs {
            b.put(&mut w, "BoundaryCondition", "Boundary Condition", k, v);
        }
    }

    let mut files = vec![
        ("case.sif".to_string(), w.finish()),
        (
            "ELMERSOLVER_STARTINFO".to_string(),
            "case.sif\n1\n".to_string(),
        ),
    ];
    for (name, text) in SquareMesh::new(2).files() {
        files.push((name.to_string(), text));
    }
    let expected_outputs = match sp.class {
        Class::Output => vec!["case_t0001.vtu".to_string()],
        _ => vec!["case.ep".to_string()],
    };
    Fixture {
        id: id.to_string(),
        class: sp.class,
        title: edf.title.clone(),
        files,
        expected_outputs,
        edf_keys: b.edf_keys,
        general_keys: b.general_keys,
    }
}

pub fn build_all(set: &EdfSet) -> Vec<Fixture> {
    EDF_FILES.iter().map(|(id, _)| build(id, set)).collect()
}

impl Fixture {
    pub fn sif(&self) -> &str {
        &self.files[0].1
    }

    /// Write the case directory (created if missing).
    pub fn write_to(&self, dir: &Path) -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        for (name, text) in &self.files {
            std::fs::write(dir.join(name), text)?;
        }
        Ok(())
    }
}
