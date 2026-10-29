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
fn refuses_radiation_and_phase_change() {
    let s = solver();
    let mut bc = ParamMap::new();
    bc.insert("Radiation", Value::Str("Idealized".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, None, None, vec![&bc])).unwrap_err()),
        "Radiation"
    );
    // `none` is no radiation at all.
    bc.insert("Radiation", Value::Str("none".into()));
    check_supported(&sections(&s, None, None, vec![&bc])).unwrap();
    // Diffuse gray radiation is the default path's.
    bc.insert("Radiation", Value::Str("diffuse gray".into()));
    check_supported(&sections(&s, None, None, vec![&bc])).unwrap();
    bc.insert("Radiation", Value::Str("weird".into()));
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
}

#[test]
fn legacy_assembly_allows_idealized_radiation_and_smart_heater() {
    let mut s = solver();
    s.insert("Legacy Assembly", Value::Logical(true));
    let mut bc = ParamMap::new();
    bc.insert("Radiation", Value::Str("idealized".into()));
    check_supported(&sections(&s, None, None, vec![&bc])).unwrap();
    // The diffuse gray model still needs view factors.
    bc.insert("Radiation", Value::Str("diffuse gray".into()));
    assert_eq!(
        keyword(check_supported(&sections(&s, None, None, vec![&bc])).unwrap_err()),
        "Radiation"
    );

    let mut bf = ParamMap::new();
    bf.insert("Smart Heater Control", Value::Logical(true));
    let mut with_bf = sections(&s, None, None, vec![]);
    with_bf.body_force = Some(&bf);
    check_supported(&with_bf).unwrap();

    // The default assembly reports it as unsupported, as HeatSolver_init does.
    s.insert("Legacy Assembly", Value::Logical(false));
    let mut with_bf = sections(&s, None, None, vec![]);
    with_bf.body_force = Some(&bf);
    assert_eq!(
        keyword(check_supported(&with_bf).unwrap_err()),
        "Smart Heater Control"
    );
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

#[test]
fn radiators_need_the_default_assembly() {
    let mut s = solver();
    let mut bc = ParamMap::new();
    bc.insert("Radiator BC", Value::Logical(true));
    check_supported(&sections(&s, None, None, vec![&bc])).unwrap();

    let mut bf = ParamMap::new();
    bf.insert(
        "Radiator Coordinates",
        Value::Real {
            rows: 1,
            cols: 3,
            data: vec![0.0, 0.0, 0.0],
        },
    );
    let mut with_bf = sections(&s, None, None, vec![]);
    with_bf.body_force = Some(&bf);
    check_supported(&with_bf).unwrap();

    // HeatSolveLegacy: "For radiative point sources use HeatSolveVec".
    s.insert("Legacy Assembly", Value::Logical(true));
    let mut with_bf = sections(&s, None, None, vec![]);
    with_bf.body_force = Some(&bf);
    assert_eq!(
        keyword(check_supported(&with_bf).unwrap_err()),
        "Radiator Coordinates"
    );
}
