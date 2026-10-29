//! Refusal of features of the two solvers that this crate does not implement.
//!
//! The kernels would silently ignore these keywords, so the drivers call the
//! checks before assembling and fail the case instead.

use elmer_core::{normalize_key, Params};

use crate::error::AdvDiffError;

/// `Procedure = "AdvectionDiffusion" "AdvectionDiffusionSolver"`.
pub const ADVDIFF_PROCEDURE: (&str, &str) = ("AdvectionDiffusion", "AdvectionDiffusionSolver");
/// `Procedure = "ModelPDE" "AdvDiffSolver"`.
pub const MODEL_PDE_PROCEDURE: (&str, &str) = ("ModelPDE", "AdvDiffSolver");

/// The sections one `AdvectionDiffusionSolver` solve reads. Absent sections
/// are skipped.
pub struct AdvDiffSections<'a> {
    pub solver: &'a dyn Params,
    pub equation: Option<&'a dyn Params>,
    pub material: Option<&'a dyn Params>,
    pub boundaries: Vec<&'a dyn Params>,
}

/// The sections one `ModelPDE` `AdvDiffSolver` solve reads.
pub struct ModelPdeSections<'a> {
    pub solver: &'a dyn Params,
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

fn scan(params: &dyn Params, section: &'static str, rules: &[Rule]) -> Result<(), AdvDiffError> {
    for r in rules {
        let hit = match r.test {
            Test::Present => params.get(r.key).is_some(),
            Test::True => params.logical(r.key) == Some(true),
        };
        if hit {
            return Err(AdvDiffError::Unsupported {
                section,
                keyword: r.key.to_string(),
                reason: r.reason,
            });
        }
    }
    Ok(())
}

const ADVDIFF_SOLVER: &[Rule] = &[
    flag(
        "Use Global Mass Matrix",
        "the global mass matrix path (FCT) is not implemented",
    ),
    flag(
        "Discontinuous Galerkin",
        "DG advection-diffusion is not implemented",
    ),
];

const ADVDIFF_MATERIAL: &[Rule] = &[present(
    "Compressibility Model",
    "compressible species transport (the div u term) is not implemented",
)];

const ADVDIFF_BOUNDARY: &[Rule] = &[
    flag(
        "Solubility Change Boundary",
        "the mixed bulk-boundary assembly is not implemented",
    ),
    present(
        "Normal Target Body",
        "the mixed bulk-boundary assembly is not implemented",
    ),
];

/// `Concentration Units`, as `ListGetString` returns it (lower case).
pub(crate) fn concentration_units(equation: &dyn Params) -> Option<String> {
    equation
        .string("Concentration Units")
        .map(|s| s.to_ascii_lowercase())
}

/// Fails with the first unsupported keyword found in `sections`.
///
/// `component` is the lower-cased solver variable name that prefixes the
/// per-species keywords (`<name> Solubility Change Boundary`).
pub fn check_advdiff_supported(
    sections: &AdvDiffSections<'_>,
    component: &str,
) -> Result<(), AdvDiffError> {
    scan(sections.solver, "solver", ADVDIFF_SOLVER)?;
    if let Some(e) = sections.equation {
        if concentration_units(e).as_deref() == Some("mass to max solubility") {
            return Err(AdvDiffError::Unsupported {
                section: "equation",
                keyword: "Concentration Units".into(),
                reason: "mass to max solubility needs the mixed bulk-boundary assembly",
            });
        }
    }
    if let Some(m) = sections.material {
        scan(m, "material", ADVDIFF_MATERIAL)?;
    }
    for bc in &sections.boundaries {
        scan(*bc, "boundary condition", ADVDIFF_BOUNDARY)?;
        if bc.logical(&format!("{component} Solubility Change Boundary")) == Some(true) {
            return Err(AdvDiffError::Unsupported {
                section: "boundary condition",
                keyword: format!("{component} Solubility Change Boundary"),
                reason: "the mixed bulk-boundary assembly is not implemented",
            });
        }
    }
    Ok(())
}

const MODEL_PDE_SOLVER: &[Rule] = &[
    flag("CutFEM", "cut finite elements are not implemented"),
    flag("Integ Test", "the integration test totals are not implemented"),
    present(
        "Relative Integration Order",
        "GaussPointsAdapt rules other than the element default are not implemented",
    ),
    present(
        "Bulk Relative Integration Order",
        "GaussPointsAdapt rules other than the element default are not implemented",
    ),
    present(
        "Element Relative Integration Order",
        "GaussPointsAdapt rules other than the element default are not implemented",
    ),
    present(
        "Number of Integration Points",
        "GaussPointsAdapt rules other than the element default are not implemented",
    ),
    present(
        "Element Integration Points",
        "GaussPointsAdapt rules other than the element default are not implemented",
    ),
    present(
        "Adaptive Integration Variable",
        "adaptive integration is not implemented",
    ),
];

/// Fails with the first unsupported keyword found in `sections`. A solver
/// `Element` keyword other than the linear `n:1` asks for p-elements, whose
/// edge and bubble degrees of freedom the nodal basis cannot carry.
pub fn check_model_pde_supported(sections: &ModelPdeSections<'_>) -> Result<(), AdvDiffError> {
    scan(sections.solver, "solver", MODEL_PDE_SOLVER)?;
    if let Some(e) = sections.solver.string("Element") {
        let e = normalize_key(e);
        if e != "n:1" {
            return Err(AdvDiffError::Unsupported {
                section: "solver",
                keyword: "Element".into(),
                reason: "p-element and higher-degree bases are not implemented",
            });
        }
    }
    for k in ["bulk gauss:", "bc gauss:"] {
        let _ = k; // namespaced integration rules are keyed `bulk gauss: <key>`; none is read here.
    }
    Ok(())
}
