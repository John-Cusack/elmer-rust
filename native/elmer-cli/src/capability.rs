//! Capability manifest and refusal checker.
//!
//! The native engine only runs cases whose every section, solver library,
//! keyword and value class appears in the manifest; anything else is refused so
//! the strangler can hand the case to the Fortran binary. The text format is
//! documented in `docs/manifest.md`.

use std::fmt;
use std::path::Path;

use elmer_core::{normalize_key, DependentKind, ParamMap, Value};

/// A `.sif` section kind (Elmer's `Section` names, lowercased in the manifest).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Section {
    Header,
    Simulation,
    Constants,
    Body,
    Material,
    BodyForce,
    InitialCondition,
    Equation,
    Solver,
    BoundaryCondition,
    Component,
    Boundary,
    RunControl,
    /// A section name the manifest cannot list; always refused.
    Other(String),
}

impl Section {
    pub fn parse(name: &str) -> Section {
        match normalize_key(name).as_str() {
            "header" => Section::Header,
            "simulation" => Section::Simulation,
            "constants" => Section::Constants,
            "body" => Section::Body,
            "material" => Section::Material,
            "body force" => Section::BodyForce,
            "initial condition" => Section::InitialCondition,
            "equation" => Section::Equation,
            "solver" => Section::Solver,
            "boundary condition" => Section::BoundaryCondition,
            "component" => Section::Component,
            "boundary" => Section::Boundary,
            "run control" => Section::RunControl,
            other => Section::Other(other.to_string()),
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Section::Header => "header",
            Section::Simulation => "simulation",
            Section::Constants => "constants",
            Section::Body => "body",
            Section::Material => "material",
            Section::BodyForce => "body force",
            Section::InitialCondition => "initial condition",
            Section::Equation => "equation",
            Section::Solver => "solver",
            Section::BoundaryCondition => "boundary condition",
            Section::Component => "component",
            Section::Boundary => "boundary",
            Section::RunControl => "run control",
            Section::Other(s) => s,
        }
    }
}

impl fmt::Display for Section {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// What kind of value a manifest entry accepts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValueClass {
    Logical,
    Integer,
    /// A constant real (scalar or array).
    Real,
    Str,
    File,
    /// `Real Variable x ...` table lookup.
    Table,
    /// `Real MATC "..."`.
    Matc,
    Lua,
    /// `Real Procedure "lib" "fn"` user function.
    Procedure,
    /// Accepts every class.
    Any,
}

impl ValueClass {
    pub fn parse(s: &str) -> Option<ValueClass> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "logical" => ValueClass::Logical,
            "integer" => ValueClass::Integer,
            "real" => ValueClass::Real,
            "string" => ValueClass::Str,
            "file" => ValueClass::File,
            "table" => ValueClass::Table,
            "matc" => ValueClass::Matc,
            "lua" => ValueClass::Lua,
            "procedure" => ValueClass::Procedure,
            "any" => ValueClass::Any,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            ValueClass::Logical => "logical",
            ValueClass::Integer => "integer",
            ValueClass::Real => "real",
            ValueClass::Str => "string",
            ValueClass::File => "file",
            ValueClass::Table => "table",
            ValueClass::Matc => "matc",
            ValueClass::Lua => "lua",
            ValueClass::Procedure => "procedure",
            ValueClass::Any => "any",
        }
    }

    /// The (never `Any`) class of a stored value.
    pub fn of(v: &Value) -> ValueClass {
        match v {
            Value::Logical(_) => ValueClass::Logical,
            Value::Integer(_) => ValueClass::Integer,
            Value::Real { .. } => ValueClass::Real,
            Value::Str(_) => ValueClass::Str,
            Value::File(_) => ValueClass::File,
            Value::Dependent(d) => match d.kind {
                DependentKind::Table { .. } => ValueClass::Table,
                DependentKind::Matc(_) => ValueClass::Matc,
                DependentKind::Procedure { .. } => ValueClass::Procedure,
                DependentKind::Lua(_) => ValueClass::Lua,
            },
        }
    }

    fn accepts(self, found: ValueClass) -> bool {
        self == ValueClass::Any || self == found
    }
}

/// One supported keyword.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Capability {
    pub section: Section,
    /// Scope to one solver library (stem, e.g. `HeatSolve`). `None`: any case.
    pub solver_lib: Option<String>,
    /// Normalised keyword; `*` matches any run of characters.
    pub keyword_pattern: String,
    pub value_class: ValueClass,
}

#[derive(Debug, PartialEq, Eq)]
pub struct ManifestError {
    /// 1-based line, 0 for I/O errors.
    pub line: usize,
    pub message: String,
}

impl fmt::Display for ManifestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(f, "{}", self.message)
        } else {
            write!(f, "line {}: {}", self.line, self.message)
        }
    }
}

impl std::error::Error for ManifestError {}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CapabilityManifest {
    pub entries: Vec<Capability>,
}

impl CapabilityManifest {
    pub fn load(path: &Path) -> Result<Self, ManifestError> {
        let text = std::fs::read_to_string(path).map_err(|e| ManifestError {
            line: 0,
            message: format!("{}: {e}", path.display()),
        })?;
        Self::parse(&text)
    }

    /// Parses `section | solver_lib | keyword_pattern | value_class` rows.
    /// Blank lines and `#` comments are skipped; `-` in `solver_lib` means
    /// unscoped.
    pub fn parse(text: &str) -> Result<Self, ManifestError> {
        let mut entries = Vec::new();
        for (i, raw) in text.lines().enumerate() {
            let line = raw.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            let err = |message: String| ManifestError { line: i + 1, message };
            let cols: Vec<&str> = line.split('|').map(str::trim).collect();
            let [sec, lib, kw, class] = cols[..] else {
                return Err(err(format!("expected 4 `|` separated columns, got {}", cols.len())));
            };
            let section = Section::parse(sec);
            if let Section::Other(s) = &section {
                return Err(err(format!("unknown section `{s}`")));
            }
            let value_class = ValueClass::parse(class)
                .ok_or_else(|| err(format!("unknown value class `{class}`")))?;
            let keyword_pattern = normalize_key(kw);
            if keyword_pattern.is_empty() {
                return Err(err("empty keyword pattern".into()));
            }
            entries.push(Capability {
                section,
                solver_lib: (lib != "-" && !lib.is_empty()).then(|| lib_stem(lib)),
                keyword_pattern,
                value_class,
            });
        }
        Ok(CapabilityManifest { entries })
    }

    /// Solver library stems named by at least one entry, lowercased.
    pub fn solver_libs(&self) -> Vec<String> {
        let mut v: Vec<String> = self
            .entries
            .iter()
            .filter_map(|e| e.solver_lib.as_ref().map(|l| l.to_ascii_lowercase()))
            .collect();
        v.sort();
        v.dedup();
        v
    }
}

/// `HeatSolve.so` / `heatsolve` -> `HeatSolve` / `heatsolve`.
fn lib_stem(lib: &str) -> String {
    let base = lib.rsplit('/').next().unwrap_or(lib);
    base.strip_suffix(".so").unwrap_or(base).to_string()
}

/// Glob match with `*` (any run, possibly empty) over an already normalised key.
pub fn pattern_matches(pattern: &str, key: &str) -> bool {
    let parts: Vec<&str> = pattern.split('*').collect();
    if parts.len() == 1 {
        return pattern == key;
    }
    let Some(mut rest) = key.strip_prefix(parts[0]) else {
        return false;
    };
    for mid in &parts[1..parts.len() - 1] {
        match rest.find(mid) {
            Some(p) => rest = &rest[p + mid.len()..],
            None => return false,
        }
    }
    rest.ends_with(parts[parts.len() - 1]) && rest.len() >= parts[parts.len() - 1].len()
}

/// One section of a case.
#[derive(Clone, Debug)]
pub struct SectionInstance {
    pub section: Section,
    /// 1-based index within its kind, as in the `.sif` (`Solver 2`).
    pub index: usize,
    pub params: ParamMap,
}

/// What the checker needs to know about a case; built from `ParamMap`
/// sections (the S1 `.sif` adapter fills it later).
#[derive(Clone, Debug, Default)]
pub struct CaseSummary {
    pub sections: Vec<SectionInstance>,
}

impl CaseSummary {
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a section; its index is the count of same-kind sections so far + 1.
    pub fn push(&mut self, section: Section, params: ParamMap) {
        let index = self.sections.iter().filter(|s| s.section == section).count() + 1;
        self.sections.push(SectionInstance { section, index, params });
    }

    /// Solver library stems used by the case's `Solver` sections (lowercased).
    pub fn solver_libs(&self) -> Vec<String> {
        self.sections
            .iter()
            .filter(|s| s.section == Section::Solver)
            .filter_map(|s| solver_lib_of(&s.params))
            .collect()
    }
}

/// The solver library a `Solver` section runs, lowercased: from `Procedure`
/// (`"Lib" "Entry"` string, or a procedure value), else from the legacy
/// `Equation` names MainUtils.F90:1225-1295 inject a procedure for.
pub fn solver_lib_of(params: &ParamMap) -> Option<String> {
    use elmer_core::Params;
    match params.get("Procedure") {
        Some(Value::Str(s)) | Some(Value::File(s)) => {
            return s.split_whitespace().next().map(|l| lib_stem(l).to_ascii_lowercase())
        }
        Some(Value::Dependent(d)) => {
            if let DependentKind::Procedure { library, .. } = &d.kind {
                return Some(lib_stem(library).to_ascii_lowercase());
            }
        }
        Some(_) => return None,
        None => {}
    }
    let eq = normalize_key(params.string("Equation")?);
    let lib = match eq.as_str() {
        "navier-stokes" => "flowsolve",
        "magnetic induction" => "magneticsolve",
        "stress analysis" => "stresssolve",
        "mesh update" => "meshsolve",
        "heat equation" => "heatsolve",
        _ => return None,
    };
    Some(lib.to_string())
}

/// Why a case is refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RefusalReason {
    UnknownSection,
    /// Solver section without a resolvable procedure.
    NoSolverProcedure,
    UnknownSolverLibrary(String),
    UnknownKeyword,
    /// The keyword is known but not for this value class.
    UnsupportedValueClass(ValueClass),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    pub section: Section,
    pub index: usize,
    /// Normalised keyword, when the refusal is about one.
    pub keyword: Option<String>,
    pub reason: RefusalReason,
}

impl fmt::Display for Refusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.section, self.index)?;
        if let Some(k) = &self.keyword {
            write!(f, ": `{k}`")?;
        }
        match &self.reason {
            RefusalReason::UnknownSection => write!(f, ": unsupported section"),
            RefusalReason::NoSolverProcedure => write!(f, ": no solver procedure"),
            RefusalReason::UnknownSolverLibrary(l) => write!(f, ": unsupported solver library `{l}`"),
            RefusalReason::UnknownKeyword => write!(f, ": unsupported keyword"),
            RefusalReason::UnsupportedValueClass(c) => {
                write!(f, ": unsupported value class `{}`", c.name())
            }
        }
    }
}

/// Checks every section, solver library and keyword of `case` against
/// `manifest`. Returns all refusals, not just the first.
pub fn check(manifest: &CapabilityManifest, case: &CaseSummary) -> Result<(), Vec<Refusal>> {
    let case_libs = case.solver_libs();
    let mut refusals = Vec::new();
    for inst in &case.sections {
        let refuse = |keyword: Option<String>, reason| Refusal {
            section: inst.section.clone(),
            index: inst.index,
            keyword,
            reason,
        };
        if matches!(inst.section, Section::Other(_)) {
            refusals.push(refuse(None, RefusalReason::UnknownSection));
            continue;
        }
        let own_lib = if inst.section == Section::Solver {
            match solver_lib_of(&inst.params) {
                None => {
                    refusals.push(refuse(None, RefusalReason::NoSolverProcedure));
                    None
                }
                Some(lib) if !manifest.solver_libs().contains(&lib) => {
                    refusals.push(refuse(None, RefusalReason::UnknownSolverLibrary(lib.clone())));
                    Some(lib)
                }
                some => some,
            }
        } else {
            None
        };
        for (key, value) in &inst.params.0 {
            let found = ValueClass::of(value);
            let mut matched_key = false;
            let mut ok = false;
            for e in manifest.entries.iter().filter(|e| e.section == inst.section) {
                let in_scope = match &e.solver_lib {
                    None => true,
                    Some(l) => {
                        let l = l.to_ascii_lowercase();
                        match &own_lib {
                            Some(own) if inst.section == Section::Solver => *own == l,
                            _ if inst.section == Section::Solver => false,
                            _ => case_libs.contains(&l),
                        }
                    }
                };
                if in_scope && pattern_matches(&e.keyword_pattern, key) {
                    matched_key = true;
                    if e.value_class.accepts(found) {
                        ok = true;
                        break;
                    }
                }
            }
            if !ok {
                let reason = if matched_key {
                    RefusalReason::UnsupportedValueClass(found)
                } else {
                    RefusalReason::UnknownKeyword
                };
                refusals.push(refuse(Some(key.clone()), reason));
            }
        }
    }
    if refusals.is_empty() {
        Ok(())
    } else {
        Err(refusals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Dependent;

    const MANIFEST: &str = "
        # demo
        Simulation | -         | simulation type   | string
        Simulation | -         | max output level  | integer
        Solver     | HeatSolve | equation          | string
        Solver     | HeatSolve | procedure         | string
        Solver     | HeatSolve | linear system *   | Any
        Material   | HeatSolve | heat conductivity | real
        Material   | HeatSolve | density           | Any
    ";

    fn case(solver_proc: &str) -> CaseSummary {
        let mut c = CaseSummary::new();
        let mut sim = ParamMap::new();
        sim.insert("Simulation Type", Value::Str("Steady State".into()));
        c.push(Section::Simulation, sim);
        let mut sol = ParamMap::new();
        sol.insert("Equation", Value::Str("Heat Equation".into()));
        sol.insert("Procedure", Value::Str(solver_proc.into()));
        sol.insert("Linear System Solver", Value::Str("Direct".into()));
        c.push(Section::Solver, sol);
        let mut mat = ParamMap::new();
        mat.insert("Heat   Conductivity", Value::real(1.0));
        c.push(Section::Material, mat);
        c
    }

    #[test]
    fn supported_case_passes() {
        let m = CapabilityManifest::parse(MANIFEST).unwrap();
        assert_eq!(m.entries.len(), 7);
        assert_eq!(check(&m, &case("HeatSolve HeatSolver")), Ok(()));
    }

    #[test]
    fn unknown_solver_library_is_refused() {
        let m = CapabilityManifest::parse(MANIFEST).unwrap();
        let r = check(&m, &case("FlowSolve FlowSolver")).unwrap_err();
        assert!(r.iter().any(|x| x.reason == RefusalReason::UnknownSolverLibrary("flowsolve".into())));
        // Keywords scoped to HeatSolve are out of scope for the flow solver.
        assert!(r.iter().any(|x| x.keyword.as_deref() == Some("linear system solver")));
    }

    #[test]
    fn unknown_keyword_and_value_class() {
        let m = CapabilityManifest::parse(MANIFEST).unwrap();
        let mut c = case("HeatSolve HeatSolver");
        c.sections[2].params.insert("Frobnicate", Value::Logical(true));
        c.sections[2].params.insert(
            "Heat Conductivity",
            Value::Dependent(Dependent {
                variables: vec!["temperature".into()],
                kind: DependentKind::Matc("1+tx".into()),
                rows: 1,
                cols: 1,
            }),
        );
        let r = check(&m, &c).unwrap_err();
        assert_eq!(r.len(), 2);
        assert!(r.contains(&Refusal {
            section: Section::Material,
            index: 1,
            keyword: Some("frobnicate".into()),
            reason: RefusalReason::UnknownKeyword,
        }));
        assert!(r.iter().any(|x| x.reason == RefusalReason::UnsupportedValueClass(ValueClass::Matc)));
    }

    #[test]
    fn legacy_equation_name_resolves_solver_library() {
        let mut p = ParamMap::new();
        p.insert("Equation", Value::Str("Heat Equation".into()));
        assert_eq!(solver_lib_of(&p).as_deref(), Some("heatsolve"));
        let mut p = ParamMap::new();
        p.insert("Equation", Value::Str("My Own".into()));
        assert_eq!(solver_lib_of(&p), None);
    }

    #[test]
    fn unknown_section_and_missing_procedure() {
        let m = CapabilityManifest::parse(MANIFEST).unwrap();
        let mut c = CaseSummary::new();
        c.push(Section::parse("Optimization"), ParamMap::new());
        c.push(Section::Solver, ParamMap::new());
        let r = check(&m, &c).unwrap_err();
        assert_eq!(r[0].reason, RefusalReason::UnknownSection);
        assert_eq!(r[1].reason, RefusalReason::NoSolverProcedure);
    }

    #[test]
    fn glob_patterns() {
        assert!(pattern_matches("linear system *", "linear system solver"));
        assert!(pattern_matches("*conductivity", "heat conductivity"));
        assert!(pattern_matches("a*b*c", "a--b--c"));
        assert!(!pattern_matches("a*b", "ab c"));
        assert!(!pattern_matches("a*a", "a"));
        assert!(pattern_matches("exact", "exact"));
        assert!(!pattern_matches("exact", "exactly"));
    }

    #[test]
    fn manifest_parse_errors() {
        let e = CapabilityManifest::parse("solver | - | x").unwrap_err();
        assert_eq!(e.line, 1);
        assert!(CapabilityManifest::parse("solver | - | x | bogus").is_err());
        assert!(CapabilityManifest::parse("nosuchsection | - | x | real").is_err());
        assert!(CapabilityManifest::load(Path::new("/nonexistent/manifest")).is_err());
    }
}
