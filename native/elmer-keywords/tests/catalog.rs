use elmer_keywords::{
    builtin, keyword_type, Catalog, DiagnosticKind, KeywordType, SCHEMA_VERSION, SNAPSHOT,
    SNAPSHOT_PATH, SNAPSHOT_REVISION,
};

/// `data/keywords.catalog.json` is generated; it must match the snapshot.
#[test]
fn checked_in_json_matches_generation() {
    let want = Catalog::from_source(SNAPSHOT, SNAPSHOT_PATH, SNAPSHOT_REVISION)
        .unwrap()
        .to_json();
    let have = include_str!("../data/keywords.catalog.json");
    assert_eq!(
        have, want,
        "regenerate: cargo run -p elmer-keywords --bin gen-keywords -- --out native/elmer-keywords/data/keywords.catalog.json"
    );
    assert!(have.contains(&format!("\"schema_version\": {SCHEMA_VERSION}")));
}

#[test]
fn generation_is_deterministic() {
    let a = Catalog::from_source(SNAPSHOT, "p", "r").unwrap().to_json();
    let b = Catalog::from_source(SNAPSHOT, "p", "r").unwrap().to_json();
    assert_eq!(a, b);
}

#[test]
fn snapshot_diagnostics_are_the_known_ones() {
    let kinds: Vec<_> = builtin()
        .diagnostics()
        .iter()
        .filter(|d| d.kind != DiagnosticKind::Duplicate)
        .map(|d| (d.line, d.kind))
        .collect();
    assert_eq!(
        kinds,
        [
            (363, DiagnosticKind::Malformed),
            (763, DiagnosticKind::Malformed),
            (814, DiagnosticKind::Malformed)
        ]
    );
    assert!(
        !builtin()
            .diagnostics()
            .iter()
            .any(|d| d.kind == DiagnosticKind::ConflictingType),
        "a conflicting type in the snapshot needs a decision before regenerating"
    );
}

/// Observations from the pristine Fortran ElmerSolver (v26.2, `CHECK KEYWORDS
/// "Silent"`, one throwaway SIF): each row is (section header, keyword,
/// listed type or None for "Unlisted keyword").
#[test]
fn agrees_with_pristine_checkkeyword() {
    use KeywordType::*;
    let cases: &[(&str, &str, Option<KeywordType>)] = &[
        ("Solver 1", "Tensor Field 3 Component 4", Some(String)),
        ("Solver 1", "Tensor Field 11 Component 4", None),
        ("Solver 1", "Exported Variable 7 DOFs", Some(Integer)),
        ("Solver 1", "Exported Variable 12 DOFs", None),
        ("Simulation", "Exported Variable 7 DOFs", None),
        ("Solver 1", "Max Parameter 4", Some(Real)),
        ("Solver 1", "Max Parameter 12", Some(Real)),
        ("Solver 1", "Max Parameter 123", Some(Real)),
        ("Solver 1", "Vector Field 2 Component 3", Some(String)),
        ("Solver 1", "Vector Field 2 Component 4", Some(String)),
        ("Simulation", "Output Variable 9", Some(String)),
        ("Simulation", "Foo Bar Baz", None),
        ("Material 1", "Water Content", None),
        ("Boundary Condition 1", "Incident Voltage", None),
        ("Boundary Condition 1", "Heat Flux BC", Some(Logical)),
        // Fortran stores the type as `logica`, so a Logical use is "wrong type".
        ("Simulation", "Levelset Timestep Directional", None),
    ];
    for (section, key, want) in cases {
        assert_eq!(keyword_type(section, key), *want, "{section} / {key}");
    }
}

#[test]
fn matc_families_expand_fully() {
    // 10x10 tensor components + 10 DOFs + 10 outputs + 10x3 vector components.
    let n = builtin()
        .entries()
        .filter(|e| e.origin == elmer_keywords::Origin::Matc)
        .count();
    assert_eq!(n, 150);
}

#[test]
fn common_heat_and_elasticity_keywords_resolve() {
    use KeywordType::*;
    let cases = [
        ("Simulation", "Simulation Type", String),
        ("Simulation", "Steady State Max Iterations", Integer),
        ("Simulation", "Timestep Sizes", Real),
        ("Solver 1", "Linear System Max Iterations", Integer),
        ("Solver 1", "Linear System Convergence Tolerance", Real),
        ("Solver 1", "Optimize Bandwidth", Logical),
        ("Material 1", "Heat Conductivity", Real),
        ("Material 1", "Youngs Modulus", Real),
        ("Body Force 1", "Heat Source", Real),
        ("Boundary Condition 1", "Temperature", Real),
        ("Body 1", "Equation", Integer),
    ];
    for (s, k, t) in cases {
        assert_eq!(keyword_type(s, k), Some(t), "{s} / {k}");
    }
}
