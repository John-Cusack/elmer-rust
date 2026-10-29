//! Refusal of model features no provider implements. `Simulation` and the
//! physics kernels never read these keywords, so a case naming one runs to
//! a wrong answer instead of failing; each is refused before setup with the
//! feature it needs. Every rule here was found by a quick test that ran
//! and missed its reference norm, and confirmed by running the pristine
//! ElmerSolver without the keyword, which then reproduced the native norm.

use elmer_core::{Params, ParamMap};
use elmer_solver::sim::parse_variable_spec;
use elmer_solver::{Model, SolverError};

enum Test {
    Present,
    True,
}

struct Rule {
    key: &'static str,
    test: Test,
    reason: &'static str,
}

const fn present(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::Present, reason }
}

const fn flag(key: &'static str, reason: &'static str) -> Rule {
    Rule { key, test: Test::True, reason }
}

const SIMULATION: &[Rule] = &[
    flag("Increase Element Order", "quadratic element generation is not implemented"),
    flag("Set Dirichlet BCs by BC Numbering", "Dirichlet conditions in BC numbering order are not implemented"),
];

/// Read by the solver's own procedure; skipped for `Exec Solver = Never`.
const SOLVER: &[Rule] = &[
    present("Mesh", "solver-specific meshes are not implemented"),
    flag("Apply Conforming BCs", "conforming BCs are not implemented"),
    flag("Apply Mortar BCs", "mortar BCs are not implemented"),
    present("Namespace", "solver keyword namespaces are not implemented"),
    flag("Adaptive Mesh Refinement", "mesh refinement is not implemented"),
];

const BOUNDARY: &[Rule] = &[
    // MapBodiesAndBCs fallbacks `Model::bc_of_tag` does not apply.
    flag("Default Target", "default BC targets are not implemented"),
    flag("Default BC", "default BC targets are not implemented"),
    flag("Default Internal BC", "default BC targets are not implemented"),
    flag("Default External BC", "default BC targets are not implemented"),
    present("Intersection BC", "intersection BCs are not implemented"),
    present("Intersection Body", "intersection BCs are not implemented"),
    flag("Discontinuous Boundary", "discontinuous boundaries (mesh splitting) are not implemented"),
    present("Mortar BC", "mortar BCs are not implemented"),
    present("Conforming BC", "conforming BCs are not implemented"),
    present("Target Nodes", "point BCs (Target Nodes) are not implemented"),
    present("Target Coordinates", "point BCs (Target Coordinates) are not implemented"),
];

/// `SetNormalizedKeywords` (MainUtils.F90:224): `<keyword> Normalize by
/// Area|Volume` divides the keyword by the area or volume it acts on.
const NORMALIZE_SUFFIXES: &[&str] = &[" normalize by area", " normalize by volume"];

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

/// Solver variables and their components.
fn variables(model: &Model) -> Vec<String> {
    let mut primary = Vec::new();
    for s in &model.solvers {
        if let Some(v) = s.string("variable") {
            if let Ok(spec) = parse_variable_spec(&v.replace('"', " "), s.integer("variable dofs")) {
                primary.push(spec.name);
                primary.extend(spec.components);
            }
        }
    }
    primary
}

/// Refuses the first unimplemented feature `model` names.
pub fn check(model: &Model) -> Result<(), SolverError> {
    scan("Simulation", 1, &model.simulation, SIMULATION)?;
    for (i, s) in model.solvers.iter().enumerate().filter(|(_, s)| runs(s)) {
        scan("Solver", i + 1, s, SOLVER)?;
    }
    for (i, b) in model.boundary_conditions.iter().enumerate() {
        scan("Boundary Condition", i + 1, b, BOUNDARY)?;
    }

    let primary = variables(model);
    // CheckPassiveElement reads `<Variable> Passive` from body forces.
    for var in &primary {
        let passive = format!("{var} passive");
        if let Some(i) = model.body_forces.iter().position(|b| b.0.contains_key(&passive)) {
            return Err(SolverError::Unsupported(format!(
                "Body Force {}: {passive}: passive elements (<Variable> Passive) are not implemented",
                i + 1
            )));
        }
        // SetDirichletBoundaries' body-force pass (Dirichlet values on bodies).
        if let Some(i) = model.body_forces.iter().position(|b| b.0.contains_key(var)) {
            return Err(SolverError::Unsupported(format!(
                "Body Force {}: {var}: body-force Dirichlet values are not implemented",
                i + 1
            )));
        }
    }

    let sections: [(&str, &[ParamMap]); 4] = [
        ("Material", &model.materials),
        ("Body Force", &model.body_forces),
        ("Boundary Condition", &model.boundary_conditions),
        ("Initial Condition", &model.initial_conditions),
    ];
    for (section, maps) in sections {
        for (i, m) in maps.iter().enumerate() {
            let normalized = m.0.keys().find(|k| NORMALIZE_SUFFIXES.iter().any(|s| k.ends_with(s)) && m.logical(k) == Some(true));
            if let Some(key) = normalized {
                return Err(SolverError::Unsupported(format!(
                    "{section} {}: {key}: keyword normalisation by area or volume is not implemented",
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

    fn refusal(model: &Model) -> String {
        match check(model) {
            Err(SolverError::Unsupported(m)) => m,
            other => panic!("{other:?}"),
        }
    }

    fn equals(var: &str) -> Value {
        Value::Dependent(Dependent {
            variables: vec![var.into()],
            kind: DependentKind::Table { x: vec![1.0], y: vec![1.0], cubic: false },
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
            simulation: pm(&[("increase element order", Value::Logical(true))]),
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Simulation 1: Increase Element Order: quadratic element generation is not implemented"
        );

        let model = Model {
            solvers: vec![heat(&[]), heat(&[("exec solver", st("never")), ("mesh", st("local"))])],
            boundary_conditions: vec![pm(&[]), pm(&[("mortar bc", Value::Integer(vec![2]))])],
            ..Model::default()
        };
        assert_eq!(refusal(&model), "Boundary Condition 2: Mortar BC: mortar BCs are not implemented");

        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("default internal bc", Value::Logical(false))])],
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
        assert_eq!(
            refusal(&model),
            "Body Force 2: temperature passive: passive elements (<Variable> Passive) are not implemented"
        );
    }

    #[test]
    fn dirichlet_variants() {
        let real = Value::Real { rows: 1, cols: 1, data: vec![1.0] };
        let bf = |k: &str| Model {
            solvers: vec![heat(&[])],
            body_forces: vec![pm(&[(k, real.clone())])],
            ..Model::default()
        };
        assert_eq!(
            refusal(&bf("temperature")),
            "Body Force 1: temperature: body-force Dirichlet values are not implemented"
        );
        // Lumped Dirichlet conditions are elmer-solver's (#121).
        assert!(check(&bf("temperature single node")).is_ok());
        assert!(check(&bf("heat source")).is_ok());

        let model = Model {
            solvers: vec![heat(&[])],
            boundary_conditions: vec![pm(&[("target nodes", Value::Integer(vec![100])), ("temperature", real.clone())])],
            ..Model::default()
        };
        assert_eq!(refusal(&model), "Boundary Condition 1: Target Nodes: point BCs (Target Nodes) are not implemented");

        let model = Model {
            solvers: vec![heat(&[])],
            materials: vec![pm(&[("heat conductivity normalize by volume", Value::Logical(true))])],
            ..Model::default()
        };
        assert_eq!(
            refusal(&model),
            "Material 1: heat conductivity normalize by volume: keyword normalisation by area or volume \
             is not implemented"
        );
    }
}
