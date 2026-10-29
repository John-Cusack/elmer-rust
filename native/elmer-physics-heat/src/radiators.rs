//! `Radiator BC`: heat from point sources (`Radiator Coordinates`,
//! `Radiator Power i`) delivered to the boundary elements by the radiator
//! factors of `elmer-viewfactors` (RadiatorFactors3D, then
//! `ReadRadiatorFactorsFromFile`), added to the flux of those elements as
//! `LocalMatrixBC` does.
//!
//! The factors are computed on the first solve. Point radiators on triangle
//! and quad surfaces of a 3D mesh are supported, with the `Viewfactor
//! Symmetry` mirror planes, and line radiators.

use std::collections::HashMap;

use elmer_core::{ElemId, Params, Value};
use elmer_solver::{SolverCtx, SolverError};
use elmer_viewfactors::{
    radiator_factors, saved_rows, surface_areas, surface_direction, ParentInfo, Radiator,
    RadiatorParams, RadiatorSurface,
};

const CALLER: &str = "RadiatorFactors";

/// A radiating surface before mirroring: its corners and the parents the
/// orientation looks at.
#[derive(Clone)]
struct Spec {
    nodes: Vec<[f64; 3]>,
    left: Option<(ParentInfo, bool)>,
    right: Option<(ParentInfo, bool)>,
    target: Option<i32>,
}

impl Spec {
    fn mirrored(&self, axis: usize, plane: f64) -> Spec {
        let flip = |mut p: [f64; 3]| {
            p[axis] = 2.0 * plane - p[axis];
            p
        };
        let parent = |p: Option<(ParentInfo, bool)>| {
            p.map(|(mut i, r)| {
                i.node = flip(i.node);
                (i, r)
            })
        };
        Spec {
            nodes: self.nodes.iter().copied().map(flip).collect(),
            left: parent(self.left),
            right: parent(self.right),
            target: self.target,
        }
    }
}

/// The radiator factors of the boundary elements of every `Radiator BC`.
pub(crate) struct Radiators {
    /// Boundary element -> factor of each radiator per unit area.
    per_element: HashMap<u32, Vec<f64>>,
    powers: Vec<f64>,
}

pub(crate) fn is_radiator_bc(ctx: &SolverCtx<'_>, e: ElemId) -> bool {
    ctx.bc_params(e)
        .is_some_and(|bc| bc.logical("Radiator BC") == Some(true))
        && ctx.mesh.element(e).etype.family() > 1
}

/// The normal `NormalVector` returns for a 3- or 4-node surface element.
fn element_normal(nodes: &[[f64; 3]]) -> [f64; 3] {
    let (du, dv): (Vec<f64>, Vec<f64>) = if nodes.len() == 3 {
        (0..3)
            .map(|i| (nodes[1][i] - nodes[0][i], nodes[2][i] - nodes[0][i]))
            .unzip()
    } else {
        (0..3)
            .map(|i| {
                (
                    ((nodes[1][i] + nodes[2][i]) - (nodes[0][i] + nodes[3][i])) / 4.0,
                    ((nodes[2][i] + nodes[3][i]) - (nodes[0][i] + nodes[1][i])) / 4.0,
                )
            })
            .unzip()
    };
    let c = [
        du[1] * dv[2] - du[2] * dv[1],
        du[2] * dv[0] - du[0] * dv[2],
        du[0] * dv[1] - du[1] * dv[0],
    ];
    let m = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
    [c[0] / m, c[1] / m, c[2] / m]
}

/// `GetParentInfo` for a parent of a 3D boundary element.
fn parent_info(
    ctx: &SolverCtx<'_>,
    parent: Option<ElemId>,
    nodes: &[u32],
) -> Option<(ParentInfo, bool)> {
    let p = parent?;
    let pe = ctx.mesh.element(p);
    let node = pe
        .nodes
        .iter()
        .rev()
        .find(|n| !nodes.contains(n))
        .copied()?;
    let body = ctx.sim.model.body_of_tag(pe.tag)?;
    let bp = ctx.sim.model.body_params(body)?;
    let rad = bp
        .body
        .logical("Radiative Body")
        .or_else(|| bp.material.and_then(|m| m.logical("Radiative Body")))
        .unwrap_or_else(|| bp.material.is_some_and(|m| m.get("Emissivity").is_some()));
    Some((
        ParentInfo {
            body: body as i32,
            node: ctx.mesh.coord(node),
        },
        rad,
    ))
}

/// The `Viewfactor Symmetry` planes, in the order `SymmetryDuplication`
/// applies them: `(axis, plane)`.
fn symmetry_planes(ctx: &SolverCtx<'_>) -> Result<Vec<(usize, f64)>, SolverError> {
    let p = ctx.params();
    let mut planes = Vec::new();
    for (i, axis) in ["x", "y", "z"].into_iter().enumerate() {
        for side in ["min", "max"] {
            let key = format!("Viewfactor Symmetry {axis} {side}");
            if let Some(v) = ctx.creal(p, &key)? {
                planes.push((i, v));
            } else if side == "min"
                && p.logical(&format!("Viewfactor Symmetry {axis}")) == Some(true)
            {
                planes.push((i, 0.0));
            }
        }
    }
    if planes.iter().any(|&(a, _)| a >= 3) {
        unreachable!();
    }
    Ok(planes)
}

impl Radiators {
    /// Computes the factors when a `Radiator BC` is present.
    pub(crate) fn setup(ctx: &SolverCtx<'_>, bnd: &[ElemId]) -> Result<Option<Self>, SolverError> {
        let radiating: Vec<ElemId> = bnd
            .iter()
            .copied()
            .filter(|&e| is_radiator_bc(ctx, e))
            .collect();
        if radiating.is_empty() {
            return Ok(None);
        }
        if ctx.mesh.dim() != 3 {
            return Err(SolverError::Unsupported(
                "radiators on a mesh that is not three dimensional".into(),
            ));
        }
        // The list of radiators: the first body force that has one, else the
        // solver section.
        let list = ctx
            .sim
            .model
            .body_forces
            .iter()
            .find(|b| b.get("Radiator Coordinates").is_some())
            .unwrap_or_else(|| ctx.params());
        let Some(Value::Real { rows, cols, data }) = list.get("Radiator Coordinates") else {
            return Err(SolverError::fatal(
                CALLER,
                "No radiators present, quitting".to_string(),
            ));
        };
        let sources: Vec<Radiator> = match *cols {
            3 => (0..*rows)
                .map(|i| Radiator::Point([data[3 * i], data[3 * i + 1], data[3 * i + 2]]))
                .collect(),
            6 => (0..*rows)
                .map(|i| Radiator::Line {
                    start: [data[6 * i], data[6 * i + 1], data[6 * i + 2]],
                    end: [data[6 * i + 3], data[6 * i + 4], data[6 * i + 5]],
                })
                .collect(),
            _ => {
                return Err(SolverError::Unsupported(
                    "Radiator Coordinates with a column count other than 3 or 6".into(),
                ))
            }
        };
        let mut powers = Vec::with_capacity(*rows);
        for i in 1..=*rows {
            powers.push(
                ctx.creal(list, &format!("Radiator Power {i}"))?
                    .unwrap_or(0.0),
            );
        }

        let mut open = false;
        let mut specs = Vec::with_capacity(radiating.len());
        for &e in &radiating {
            let el = ctx.mesh.element(e);
            if !(el.nodes.len() == 3 || el.nodes.len() == 4) {
                return Err(SolverError::Unsupported(
                    "radiators on elements other than triangles and quadrilaterals".into(),
                ));
            }
            let bc = ctx
                .bc_params(e)
                .expect("radiating elements have a condition");
            open |= bc.logical("Radiation Boundary Open") == Some(true);
            specs.push(Spec {
                nodes: el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect(),
                left: parent_info(ctx, el.parents[0], el.nodes),
                right: parent_info(ctx, el.parents[1], el.nodes),
                target: bc
                    .integer("Radiation Target Body")
                    .or_else(|| bc.integer("Normal Target Body")),
            });
        }
        let n_orig = specs.len();

        let planes = symmetry_planes(ctx)?;
        let mut all = specs.clone();
        for &(axis, plane) in &planes {
            let mirrored: Vec<Spec> = all.iter().map(|s| s.mirrored(axis, plane)).collect();
            all.extend(mirrored);
        }

        let surfaces: Vec<RadiatorSurface> = all
            .iter()
            .map(|s| {
                let nrm = element_normal(&s.nodes);
                let k = s.nodes.len() as f64;
                let centre = [0, 1, 2].map(|i| s.nodes.iter().map(|p| p[i]).sum::<f64>() / k);
                let dir = surface_direction(
                    nrm,
                    centre,
                    s.left.map(|l| l.0),
                    s.right.map(|r| r.0),
                    s.target,
                    (s.left.is_some_and(|l| l.1), s.right.is_some_and(|r| r.1)),
                )
                .map_err(|err| SolverError::fatal("ViewFactors", err.to_string()))?;
                Ok(RadiatorSurface {
                    nodes: s.nodes.clone(),
                    normal: if dir > 0.0 { nrm } else { nrm.map(|c| -c) },
                })
            })
            .collect::<Result<_, SolverError>>()?;

        let p = ctx.params();
        let count = |key: &str| p.integer(key).map(|n| n.max(1) as usize);
        let params = RadiatorParams::from_keywords(
            p.real("Viewfactor Area Tolerance"),
            p.real("Viewfactor Factor Tolerance"),
            count("Viewfactor Number Of Rays"),
            count("Viewfactor Line Integration Points"),
            count("Viewfactor Quad Integration Points"),
            count("Viewfactor Triangle Integration Points"),
        );
        let nrad = sources.len();
        let mut n = surfaces.len();
        let mut f = radiator_factors(&surfaces, &sources, &params)
            .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
        let areas = surface_areas(&surfaces[..n_orig], &params)
            .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;

        // SymmetryReduction: fold the mirrored half of each doubling back.
        for _ in &planes {
            let half = n / 2;
            let mut g = Vec::with_capacity(nrad * half);
            for i in 0..nrad {
                for j in 0..half {
                    g.push(f[i * n + j] + f[i * n + j + half]);
                }
            }
            f = g;
            n = half;
        }

        let min_factor = p.real("Minimum View Factor").unwrap_or(1.0e-20);
        for x in &mut f {
            if *x < min_factor {
                *x = 0.0;
            }
        }
        if !open {
            for row in f.chunks_mut(n) {
                let s: f64 = row.iter().sum();
                row.iter_mut().for_each(|x| *x /= s);
            }
        }
        let rows = saved_rows(&f, n, nrad, min_factor, false);

        let mut per_element: HashMap<u32, Vec<f64>> = HashMap::new();
        for (i, row) in rows.iter().enumerate() {
            for &(j, v) in row {
                per_element
                    .entry(radiating[j].0)
                    .or_insert_with(|| vec![0.0; nrad])[i] = v / areas[j];
            }
        }
        Ok(Some(Self {
            per_element,
            powers,
        }))
    }

    /// `SUM(RadiatorPowers * Element % BoundaryInfo % Radiators)`.
    pub(crate) fn flux(&self, e: ElemId) -> Option<f64> {
        let f = self.per_element.get(&e.0)?;
        Some(f.iter().zip(&self.powers).map(|(a, b)| a * b).sum())
    }
}
