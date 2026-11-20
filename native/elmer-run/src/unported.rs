//! Refusal of model features no provider implements. `Simulation` and the
//! physics kernels never read these keywords, so a case naming one runs to
//! a wrong answer instead of failing; each is refused before setup with the
//! feature it needs. Every rule here was found by a quick test that ran
//! and missed its reference norm, and confirmed by running the pristine
//! ElmerSolver without the keyword, which then reproduced the native norm.

use elmer_core::{normalize_key, FieldKind, ParamMap, Params, Value};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Model, SolverError};

enum Test {
    // No rule tests a bare presence at the moment.
    #[allow(dead_code)]
    Present,
    True,
}

struct Rule {
    key: &'static str,
    test: Test,
    reason: &'static str,
}

#[allow(dead_code)]
const fn present(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::Present, reason }
}

const fn flag(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::True, reason }
}

const SIMULATION: &[Rule] = &[];

/// Read by the solver's own procedure; skipped for `Exec Solver = Never`.
const SOLVER: &[Rule] = &[
    flag("Adaptive Mesh Refinement", "mesh refinement is not implemented"),
];

/// `Solver i :: Mesh`: a solver of its own mesh (`elmer_solver::multimesh`)
/// names one mesh directory under the `Mesh DB`. Not composed: the
/// keyword's other forms.
fn check_solver_mesh(i: usize, s: &ParamMap) -> Result<(), SolverError> {
    if s.get("mesh").is_none() {
        return Ok(());
    }
    let refuse = |what: &str| Err(SolverError::Unsupported(format!("Solver {i}: Mesh: {what}")));
    let Some(name) = s.string("mesh") else { return refuse("not a mesh name") };
    if name.split_whitespace().count() != 1 {
        return refuse("a mesh with options or a directory and a name");
    }
    Ok(())
}

const BODY_FORCE: &[Rule] = &[];

const BOUNDARY: &[Rule] = &[];

fn scan(section: &str, i: usize, params: &ParamMap, rules: &[Rule]) -> Result<(), SolverError> {
    for r in rules {
        let hit = match r.test {
            Test::Present => params.get(r.key).is_some(),
            Test::True => params.logical(r.key) == Some(true),
        };
        if hit {
            return Err(SolverError::Unsupported(format!("{section} {i}: {}: {}", r.key, r.reason)));
        }
    }
    Ok(())
}

fn runs(solver: &ParamMap) -> bool {
    !solver.0.is_empty() && !solver.string("exec solver").is_some_and(|s| s.eq_ignore_ascii_case("never"))
}

/// Refuses the first unimplemented feature `model` names.
pub fn check(model: &Model) -> Result<(), SolverError> {
    scan("Simulation", 1, &model.simulation, SIMULATION)?;
    for (i, s) in model.solvers.iter().enumerate().filter(|(_, s)| runs(s)) {
        // HeatSolve, ElasticSolve, StressSolve, StatElecSolve and StatCurrentSolve refine through elmer-adaptive (its own checks refuse
        // the remeshing paths); the other solvers do not yet.
        let heat = elmer_solver::solver::procedure_of(s).is_some_and(|(l, e)| {
(l.eq_ignore_ascii_case("heatsolve") && e.eq_ignore_ascii_case("heatsolver"))
                || (l.eq_ignore_ascii_case("elasticsolve") && e.eq_ignore_ascii_case("elasticsolver"))
                || (l.eq_ignore_ascii_case("stresssolve") && e.eq_ignore_ascii_case("stresssolver"))
                || (l.eq_ignore_ascii_case("statelecsolve") && e.eq_ignore_ascii_case("statelecsolver"))
                || (l.eq_ignore_ascii_case("statcurrentsolve") && e.eq_ignore_ascii_case("statcurrentsolver"))
        });
        if heat && s.logical("adaptive mesh refinement") == Some(true) {
            let mut s2 = s.clone();
            s2.insert("Adaptive Mesh Refinement", Value::Logical(false));
            scan("Solver", i + 1, &s2, SOLVER)?;
            continue;
        }
        // A CutFEM solver's `Project Variable k` names the fields of its cut
        // mesh (`CreateCutFEMMesh`), not an integration-point projection.
        if s.logical("cutfem") == Some(true) && s.get("project variable 1").is_some() {
            let mut s2 = s.clone();
            s2.0.remove("project variable 1");
            scan("Solver", i + 1, &s2, SOLVER)?;
            check_solver_mesh(i + 1, s)?;
            continue;
        }
        scan("Solver", i + 1, s, SOLVER)?;
        check_solver_mesh(i + 1, s)?;
    }
    for (i, b) in model.body_forces.iter().enumerate() {
        scan("Body Force", i + 1, b, BODY_FORCE)?;
    }
    for (i, b) in model.boundary_conditions.iter().enumerate() {
        scan("Boundary Condition", i + 1, b, BOUNDARY)?;
    }

    // An `-ip` exported variable read by a boundary condition: Elmer
    // evaluates it at the boundary element's own Gauss points
    // (`InterpolateIPVariableToBoundary`). HeatSolver's `Heat Flux` does
    // (`SolverCtx::boundary_real_at_ips`); other keys read it at the nodes.
    let ip_vars: Vec<String> = model
        .solvers
        .iter()
        .flat_map(|s| s.0.iter())
        .filter(|(k, _)| k.starts_with("exported variable") && k.split(' ').count() == 3)
        .filter_map(|(_, v)| match v {
            Value::Str(t) => parse_variable_spec(&t.replace('"', " "), None).ok(),
            _ => None,
        })
        .filter(|spec| spec.kind == FieldKind::Ip)
        .map(|spec| spec.name)
        .collect();
    let heat_solver = model.solvers.iter().any(|s| {
        elmer_solver::solver::procedure_of(s)
            .is_some_and(|(l, e)| l.eq_ignore_ascii_case("heatsolve") && e.eq_ignore_ascii_case("heatsolver"))
    });
    for (i, bc) in model.boundary_conditions.iter().enumerate() {
        for (key, v) in &bc.0 {
            if heat_solver && key == "heat flux" {
                continue;
            }
            let Value::Dependent(d) = v else { continue };
            let mut names = d.variables.iter().flat_map(|v| v.split(',')).map(normalize_key);
            if let Some(var) = names.find(|x| ip_vars.contains(x)) {
                return Err(SolverError::Unsupported(format!(
                    "Boundary Condition {}: {key} depends on the -ip variable {var}: boundary evaluation of integration-point variables is not implemented",
                    i + 1
                )));
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{Dependent, DependentKind, Value};

    fn pm(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    fn st(s: &str) -> Value {
        Value::Str(s.into())
    }

    fn heat(extra: &[(&str, Value)]) -> ParamMap {
        let mut p = pm(&[("procedure", st("\"HeatSolve\" \"HeatSolver\"")), ("variable", st("Temperature"))]);
        for (k, v) in extra {
            p.insert(k, v.clone());
        }
        p
    }

    fn equals(var: &str) -> Value {
        Value::Dependent(Dependent {
            variables: vec![var.into()],
            kind: DependentKind::Table { x: vec![1.0], y: vec![1.0], cubic: false, monotone: false },
            rows: 1,
            cols: 1,
        })
    }

    #[test]
    fn plain_heat_case_passes() {
        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("temperature", Value::Real { rows: 1, cols: 1, data: vec![0.0] })])],
            ..Model::default()
        };
        assert!(check(&model).is_ok());
    }

    #[test]
    fn keyword_rules() {
        let model = Model {
            simulation: pm(&[("set dirichlet bcs by bc numbering", Value::Logical(true))]),
            ..Model::default()
        };
        // Dirichlet conditions in BC numbering order are elmer-solver's.
        assert!(check(&model).is_ok());

        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("default internal bc", Value::Logical(true))])],
            ..Model::default()
        };
        assert!(check(&model).is_ok());

    }

    #[test]
    fn variable_named_keywords() {
        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("temperature load", equals("tempb loads"))])],
            ..Model::default()
        };
        // Nodal loads are elmer-solver's (SetNodalLoads).
        assert!(check(&model).is_ok());

        let model = Model {
            solvers: vec![heat(&[])],
            body_forces: vec![pm(&[]), pm(&[("temperature passive", Value::Logical(true))])],
            ..Model::default()
        };
        // Passive elements are elmer-solver's (CheckPassiveElement).
        assert!(check(&model).is_ok());
    }

    #[test]
    fn dirichlet_variants() {
        let real = Value::Real { rows: 1, cols: 1, data: vec![1.0] };
        let bf = |k: &str| Model {
            solvers: vec![heat(&[])],
            body_forces: vec![pm(&[(k, real.clone())])],
            ..Model::default()
        };
        // Body-force Dirichlet values are elmer-solver's (SetDirichletBoundaries).
        assert!(check(&bf("temperature")).is_ok());
        assert!(check(&bf("heat source")).is_ok());

        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("target nodes", Value::Integer(vec![100])), ("temperature", real.clone())])],
            ..Model::default()
        };
        // Point Dirichlet conditions are elmer-solver's.
        assert!(check(&model).is_ok());

        let model = Model {
            solvers: vec![heat(&[])],
            materials: vec![pm(&[("heat conductivity normalize by volume", Value::Logical(true))])],
            ..Model::default()
        };
        // Keyword normalisation is elmer-solver's (SetNormalizedKeywords).
        assert!(check(&model).is_ok());
    }
}
