//! The `Post File = "x.vtu"` auto-solver (`AddVtuOutputSolverHack`,
//! ElmerSolver.F90:1449-1509) and the format/prefix defaults of
//! `ResultOutputSolver` (ResultOutputSolver.F90:80-146), as pure functions.

use elmer_core::{normalize_key, ParamMap, Params, Value};

/// Keyword prefix whose entries `Simulation` forwards to the generated
/// solver with the prefix removed (`ListCopyPrefixedKeywords( Simu, Params, 'vtu:' )`).
const PREFIX: &str = "vtu:";

/// A string stored the way `ListAddString` does by default: lower-cased.
fn lowered(s: &str) -> Value {
    Value::Str(s.to_ascii_lowercase())
}

/// The solver section `AddVtuOutputSolverHack` appends for `Post File`,
/// or `None` when it does nothing.
///
/// * `simulation` is the Simulation section; `solvers` are the Solver
///   sections already in the model, in order.
/// * The hack applies when `Post File` contains `.vtu` (first occurrence; the
///   prefix before it becomes `Output File Name`).
/// * It does nothing when the first solver whose `Procedure` contains
///   `ResultOutputSolver` already has `Vtu Format` or `Output Format = vtu`.
/// * The caller must also remove `Post File` from the Simulation section
///   (Fortran `ListRemove`) when this returns `Some`, and call this at most
///   once per run (Fortran guards with a `SAVE`d flag).
///
/// Result keywords, in Fortran's order: `Procedure`, `Equation`, `Output
/// Format`, `Output File Name`, `No Matrix`, `Variable`, `Exec Solver`,
/// `Save Geometry IDs`; then every `vtu:`-prefixed Simulation keyword
/// (prefix stripped, overriding the above); then `ascii output` and `binary
/// output` from the Simulation section unless already present.
pub fn post_file_solver(simulation: &ParamMap, solvers: &[ParamMap]) -> Option<ParamMap> {
    let post = simulation.string("post file")?;
    let k = post.find(".vtu")?;

    if let Some(existing) = solvers
        .iter()
        .find(|s| s.string("procedure").is_some_and(|p| p.contains("ResultOutputSolver")))
    {
        let vtu_format = existing.logical("vtu format").unwrap_or(false)
            || existing.string("output format") == Some("vtu");
        if vtu_format {
            return None;
        }
    }

    let mut p = ParamMap::new();
    // `ListAddString(..., .FALSE.)`: no case conversion for these two.
    p.insert("Procedure", Value::Str("ResultOutputSolve ResultOutputSolver".to_string()));
    p.insert("Equation", lowered("InternalVtuOutputSolver"));
    p.insert("Output Format", lowered("vtu"));
    p.insert("Output File Name", Value::Str(post[..k].to_string()));
    p.insert("No Matrix", Value::Logical(true));
    if p.get("variable").is_none() {
        p.insert("Variable", lowered("-global vtu_internal_dummy"));
    }
    p.insert("Exec Solver", lowered("after saving"));
    p.insert("Save Geometry IDs", Value::Logical(true));

    // ListCopyPrefixedKeywords: names longer than the prefix that start with it;
    // one optional blank after the prefix is dropped.
    for (name, value) in &simulation.0 {
        if name.len() > PREFIX.len() && name.starts_with(&normalize_key(PREFIX)) {
            let rest = name[PREFIX.len()..].trim_start_matches(' ');
            p.insert(rest, value.clone());
        }
    }

    for key in ["ascii output", "binary output"] {
        if p.get(key).is_none() {
            if let Some(v) = simulation.get(key) {
                p.insert(key, v.clone());
            }
        }
    }
    Some(p)
}

/// Output formats `ResultOutputSolver` selects from its section.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OutputFormats {
    pub gid: bool,
    pub gmsh: bool,
    pub vtk: bool,
    pub vtu: bool,
    pub dx: bool,
    pub elmerpost: bool,
    pub stl: bool,
}

/// `ResultOutputSolver`'s format flags (`<Name> Format = True` and `Output
/// Format = "<name>"`). `Err(name)` for an unknown `Output Format`, where
/// Fortran warns and returns without writing anything. With no format given
/// VTU is assumed.
pub fn result_output_formats(p: &dyn Params) -> Result<OutputFormats, String> {
    let mut f = OutputFormats {
        gid: p.logical("gid format").unwrap_or(false),
        gmsh: p.logical("gmsh format").unwrap_or(false),
        vtk: p.logical("vtk format").unwrap_or(false),
        vtu: p.logical("vtu format").unwrap_or(false),
        dx: p.logical("dx format").unwrap_or(false),
        elmerpost: p.logical("elmerpost format").unwrap_or(false),
        stl: p.logical("stl format").unwrap_or(false),
    };
    if let Some(name) = p.string("output format") {
        match name {
            "gid" => f.gid = true,
            "vtk" => f.vtk = true,
            "vtu" => f.vtu = true,
            "dx" => f.dx = true,
            "gmsh" => f.gmsh = true,
            "elmerpost" => f.elmerpost = true,
            "stl" => f.stl = true,
            other => return Err(other.to_string()),
        }
    }
    let any = f.gid || f.gmsh || f.vtk || f.vtu || f.dx || f.elmerpost || f.stl;
    if !any {
        f.vtu = true;
    }
    Ok(f)
}

/// `Output File Name`, defaulting to `Case` as `ResultOutputSolver` does
/// before it calls any writer (the writer's own default, `Output`, only
/// applies when `VtuOutputSolver` is a solver of its own).
pub fn result_output_prefix(p: &dyn Params) -> String {
    p.string("output file name").unwrap_or("Case").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sim(post: &str) -> ParamMap {
        let mut s = ParamMap::new();
        s.insert("Post File", Value::Str(post.to_string()));
        s
    }

    #[test]
    fn builds_the_hidden_solver() {
        let mut s = sim("run/case.vtu");
        s.insert("vtu: Ascii Output", Value::Logical(true));
        s.insert("vtu:Skip Metainfo", Value::Logical(true));
        s.insert("Binary Output", Value::Logical(false));
        let p = post_file_solver(&s, &[]).unwrap();
        assert_eq!(p.string("procedure"), Some("ResultOutputSolve ResultOutputSolver"));
        assert_eq!(p.string("equation"), Some("internalvtuoutputsolver"));
        assert_eq!(p.string("output format"), Some("vtu"));
        assert_eq!(p.string("output file name"), Some("run/case"));
        assert_eq!(p.logical("no matrix"), Some(true));
        assert_eq!(p.string("variable"), Some("-global vtu_internal_dummy"));
        assert_eq!(p.string("exec solver"), Some("after saving"));
        assert_eq!(p.logical("save geometry ids"), Some(true));
        assert_eq!(p.logical("ascii output"), Some(true));
        assert_eq!(p.logical("skip metainfo"), Some(true));
        // Inherited from Simulation because the solver did not set it.
        assert_eq!(p.logical("binary output"), Some(false));
    }

    #[test]
    fn prefixed_keyword_overrides_and_simulation_does_not_override() {
        let mut s = sim("case.vtu");
        s.insert("vtu: Output File Name", Value::Str("other".into()));
        s.insert("vtu: Ascii Output", Value::Logical(true));
        s.insert("Ascii Output", Value::Logical(false));
        let p = post_file_solver(&s, &[]).unwrap();
        assert_eq!(p.string("output file name"), Some("other"));
        assert_eq!(p.logical("ascii output"), Some(true));
    }

    #[test]
    fn not_applicable_cases() {
        assert!(post_file_solver(&ParamMap::new(), &[]).is_none());
        assert!(post_file_solver(&sim("case.ep"), &[]).is_none());

        let mut existing = ParamMap::new();
        existing.insert("Procedure", Value::Str("ResultOutputSolve ResultOutputSolver".into()));
        existing.insert("Output Format", Value::Str("vtu".into()));
        assert!(post_file_solver(&sim("case.vtu"), &[existing.clone()]).is_none());

        // A ResultOutputSolver that writes something else does not block it.
        existing.insert("Output Format", Value::Str("gmsh".into()));
        assert!(post_file_solver(&sim("case.vtu"), &[existing]).is_some());
    }

    #[test]
    fn formats_and_prefix_defaults() {
        let mut p = ParamMap::new();
        assert_eq!(result_output_formats(&p).unwrap(), OutputFormats { vtu: true, ..Default::default() });
        assert_eq!(result_output_prefix(&p), "Case");
        p.insert("Output Format", Value::Str("gmsh".into()));
        p.insert("VTK Format", Value::Logical(true));
        p.insert("Output File Name", Value::Str("res".into()));
        let f = result_output_formats(&p).unwrap();
        assert!(f.gmsh && f.vtk && !f.vtu);
        assert_eq!(result_output_prefix(&p), "res");
        p.insert("Output Format", Value::Str("png".into()));
        assert_eq!(result_output_formats(&p), Err("png".to_string()));
    }
}
