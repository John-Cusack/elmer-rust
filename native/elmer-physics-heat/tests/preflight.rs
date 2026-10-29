use elmer_core::{ParamMap, Value};
use elmer_physics_heat::{check_supported, HeatError, HeatSections, PROCEDURE};

fn solver() -> ParamMap {
    let mut p = ParamMap::new();
    p.insert("Equation", Value::Str("Heat Equation".into()));
    p.insert("Procedure", Value::Str("HeatSolve HeatSolver".into()));
    p
}

fn sections<'a>(
    solver: &'a ParamMap,
    equation: Option<&'a ParamMap>,
    material: Option<&'a ParamMap>,
    bcs: Vec<&'a ParamMap>,
) -> HeatSections<'a> {
    HeatSections {
        solver,
        equation: equation.map(|p| p as _),
        material: material.map(|p| p as _),
        body_force: None,
        boundaries: bcs.into_iter().map(|p| p as _).collect(),
    }
}

fn keyword(e: HeatError) -> String {
    match e {
        HeatError::Unsupported { keyword, .. } => keyword,
        e => panic!("{e}"),
    }
}

#[test]
fn registry_entry() {
    assert_eq!(PROCEDURE, ("HeatSolve", "HeatSolver"));
}

#[test]
fn plain_conduction_case_passes() {
    let s = solver();
    let mut m = ParamMap::new();
    m.insert("Heat Conductivity", Value::real(1.0));
    let mut bc = ParamMap::new();
    bc.insert("Heat Flux", Value::real(1.0));
    bc.insert("Heat Transfer Coefficient", Value::real(2.0));
    check_supported(&sections(&s, None, Some(&m), vec![&bc])).unwrap();
}

#[test]
fn refuses_radiation_phase_change_and_legacy() {
    let mut s = solver();
    let mut bc = ParamMap::new();
    bc.insert("Radiation", Value::Str("Idealized".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, None, None, vec![&bc])).unwrap_err()),
        "Radiation"
    );

    let mut m = ParamMap::new();
    m.insert("phase  change MODEL", Value::Str("Spatial 2".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, None, Some(&m), vec![])).unwrap_err()),
        "Phase Change Model"
    );

    s.insert("Legacy Assembly", Value::Logical(true));
    assert_eq!(
        keyword(check_supported(&sections(&s, None, None, vec![])).unwrap_err()),
        "Legacy Assembly"
    );
    // A false flag is fine.
    s.insert("Legacy Assembly", Value::Logical(false));
    check_supported(&sections(&s, None, None, vec![])).unwrap();
}

#[test]
fn convection_needs_supg() {
    let mut s = solver();
    let mut eq = ParamMap::new();
    eq.insert("Convection", Value::Str("Computed".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, Some(&eq), None, vec![])).unwrap_err()),
        "Stabilize"
    );
    s.insert("Stabilize", Value::Logical(true));
    check_supported(&sections(&s, Some(&eq), None, vec![])).unwrap();
    s.insert("Stabilization Method", Value::Str("Bubbles".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, Some(&eq), None, vec![])).unwrap_err()),
        "Stabilize"
    );
    eq.insert("Convection", Value::Str("Weird".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, Some(&eq), None, vec![])).unwrap_err()),
        "Convection"
    );
}
