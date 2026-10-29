//! Refusal of HeatSolve features this crate does not implement.
//!
//! The kernels would silently ignore these keywords, so the driver calls
//! [`check_supported`] before assembling and fails the case instead.

use elmer_core::{Params, Value};

use crate::error::HeatError;

/// Procedure the solver registry binds to this crate: `Procedure = "HeatSolve"
/// "HeatSolver"`.
pub const PROCEDURE: (&str, &str) = ("HeatSolve", "HeatSolver");

/// The SIF sections one heat solve reads. Absent sections are skipped.
pub struct HeatSections<'a> {
    pub solver: &'a dyn Params,
    pub equation: Option<&'a dyn Params>,
    pub material: Option<&'a dyn Params>,
    pub body_force: Option<&'a dyn Params>,
    pub boundaries: Vec<&'a dyn Params>,
}

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
    Rule {
        key,
        test: Test::Present,
        reason,
    }
}

const fn flag(key: &'static str, reason: &'static str) -> Rule {
    Rule {
        key,
        test: Test::True,
        reason,
    }
}

const SOLVER: &[Rule] = &[
    flag("Bubbles", "bubble degrees of freedom are not implemented"),
    flag("Discontinuous Galerkin", "DG heat is not implemented"),
    flag("DG Reduced Basis", "DG heat is not implemented"),
    flag("Radiosity Model", "radiation is not implemented"),
    flag("Spectral Model", "radiation is not implemented"),
    flag(
        "Adaptive Mesh Refinement",
        "mesh refinement is not implemented",
    ),
    flag(
        "Calculate Conductance Matrix",
        "constraint modes are not implemented",
    ),
    present("Radiator Coordinates", "radiators are not implemented"),
    flag(
        "Nonlinear System Linesearch",
        "DefaultLinesearch is not implemented",
    ),
];

const EQUATION: &[Rule] = &[
    present("Phase Change Model", "phase change is not implemented"),
    present(
        "Check Latent Heat Release",
        "latent heat substepping is not implemented",
    ),
];

const MATERIAL: &[Rule] = &[
    present("Phase Change Model", "phase change is not implemented"),
    present("Draw Velocity", "draw velocity is not implemented"),
    present("Effective Heat Capacity", "phase change is not implemented"),
    present("Enthalpy", "phase change is not implemented"),
    present("Specific Enthalpy", "phase change is not implemented"),
    flag("Heat Conductivity At IP", AT_IP),
    flag("Density At IP", AT_IP),
    flag("Heat Capacity At IP", AT_IP),
    flag("Heat Transfer Multiplier At IP", AT_IP),
];

/// `ListInitElementKeyword`'s `<name> At IP`: evaluate the property at the
/// integration points instead of interpolating its nodal values.
const AT_IP: &str = "evaluation at integration points is not implemented";

const BODY_FORCE: &[Rule] = &[
    present("Integral Heat Source", "not implemented"),
    present("Friction Heat", "not implemented"),
    flag("Joule Heat", "coupled Joule heating is not implemented"),
    present("Perfusion Rate", "perfusion is not implemented"),
    present("Radiator Coordinates", "radiators are not implemented"),
    flag("Heat Source At IP", AT_IP),
    flag("Volumetric Heat Source At IP", AT_IP),
];

const BOUNDARY: &[Rule] = &[
    present(
        "Farfield Temperature",
        "far-field condition is not implemented",
    ),
    present(
        "Infinity BC Temperature",
        "far-field condition is not implemented",
    ),
    present("Phase Change", "phase change is not implemented"),
    present("Radiator BC", "radiators are not implemented"),
    present("Heat Gap", "DG heat gap is not implemented"),
    present("Temperature Nitsche", "weak conditions are not implemented"),
];

fn scan(params: &dyn Params, section: &'static str, rules: &[Rule]) -> Result<(), HeatError> {
    for r in rules {
        let hit = match r.test {
            Test::Present => params.get(r.key).is_some(),
            Test::True => params.logical(r.key) == Some(true),
        };
        if hit {
            return Err(HeatError::Unsupported {
                section,
                keyword: r.key.to_string(),
                reason: r.reason,
            });
        }
    }
    Ok(())
}

/// `Radiation`: `Idealized` is the Newton-linearised grey radiation to a given
/// external temperature, which the legacy path implements; `Diffuse Gray`
/// (view factors, then Gebhart factors) the default path does, for
/// axisymmetric geometry.
fn check_radiation(bc: &dyn Params, legacy: bool) -> Result<(), HeatError> {
    let Some(v) = bc.get("Radiation") else {
        return Ok(());
    };
    let ok = match v {
        Value::Str(s) => {
            let s = s.to_ascii_lowercase();
            s == "none" || (legacy && s == "idealized") || (!legacy && s == "diffuse gray")
        }
        _ => false,
    };
    if ok {
        Ok(())
    } else {
        Err(HeatError::Unsupported {
            section: "boundary condition",
            keyword: "Radiation".into(),
            reason: "radiation is not implemented",
        })
    }
}

/// Elmer's `GetStabilizeFlag`: `Stabilize`, overridden by `Stabilization
/// Method` (`stabilized` on, `bubbles` and `none` off).
pub(crate) fn stabilized(solver: &dyn Params) -> bool {
    let mut s = solver.logical("Stabilize").unwrap_or(false);
    if let Some(m) = solver.string("Stabilization Method") {
        match m.to_ascii_lowercase().as_str() {
            "stabilized" => s = true,
            "bubbles" | "none" => s = false,
            _ => {}
        }
    }
    s
}

/// Fails with the first unsupported keyword found in `sections`.
///
/// Convection is allowed only with SUPG (`Stabilize = True`): without it
/// Elmer adds a residual-free bubble, which this crate does not implement.
pub fn check_supported(sections: &HeatSections<'_>) -> Result<(), HeatError> {
    let legacy = sections.solver.logical("Legacy Assembly") == Some(true);
    scan(sections.solver, "solver", SOLVER)?;
    if let Some(e) = sections.equation {
        scan(e, "equation", EQUATION)?;
        if let Some(c) = e.string("Convection") {
            match c.to_ascii_lowercase().as_str() {
                "none" => {}
                "constant" | "computed" => {
                    if !stabilized(sections.solver) {
                        return Err(HeatError::Unsupported {
                            section: "solver",
                            keyword: "Stabilize".into(),
                            reason: "convection without SUPG needs a residual-free bubble, \
                                     which is not implemented",
                        });
                    }
                }
                _ => {
                    return Err(HeatError::Unsupported {
                        section: "equation",
                        keyword: "Convection".into(),
                        reason: "only none, constant and computed are implemented",
                    })
                }
            }
        }
    }
    if let Some(m) = sections.material {
        scan(m, "material", MATERIAL)?;
    }
    if let Some(b) = sections.body_force {
        scan(b, "body force", BODY_FORCE)?;
        if b.get("Smart Heater Control").is_some() && !legacy {
            // HeatSolver_init: `ListWarnUnsupportedKeyword` with FatalFound.
            return Err(HeatError::Unsupported {
                section: "body force",
                keyword: "Smart Heater Control".into(),
                reason: "only the Legacy Assembly path implements it",
            });
        }
    }
    for bc in &sections.boundaries {
        scan(*bc, "boundary condition", BOUNDARY)?;
        check_radiation(*bc, legacy)?;
    }
    Ok(())
}
