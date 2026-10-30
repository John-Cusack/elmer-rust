//! Keyword normalisation by area or volume (`SetNormalizedKeywords`,
//! MainUtils.F90:226).
//!
//! A real keyword tagged by `<keyword> Normalize by Area|Volume = True`
//! (`ListTagKeywords`, Lists.F90) or by the `-distribute` prefix is divided
//! by the measure of the entity it sits in: the summed area or volume of
//! the elements of its body, body force, material or boundary condition
//! (`CalculateEntityWeights`, SolverBasics.F90:10444). `elmer-sif`'s
//! `-distribute` reaches here as `<keyword> Normalize by Area = True`
//! (`sif::to_model`). Tags in other sections change nothing, as in
//! `ListSetParameters`.

use elmer_core::{ElemId, MeshView, ParamMap, Params, Value};
use elmer_fe::{CoordinateKind, CoordinateSystem};
use elmer_solver::{FeBackend, Model, SolverError};

pub(crate) const SUFFIXES: [&str; 2] = [" normalize by area", " normalize by volume"];

/// `ListTagEntry` for both suffixes: the keywords each suffix keyword
/// tags (any keyword whose name starts its name, before the suffix), with
/// the suffix keyword's value, the later suffix keyword winning.
fn tagged(p: &ParamMap) -> Vec<String> {
    let mut tags: Vec<(String, bool)> = Vec::new();
    for suffix in SUFFIXES {
        for (key, v) in &p.0 {
            if !key.ends_with(suffix) {
                continue;
            }
            let on = matches!(v, Value::Logical(true));
            let stem = key.len() - suffix.len();
            for k2 in p.0.keys().filter(|k2| k2.len() <= stem && key.starts_with(k2.as_str())) {
                match tags.iter_mut().find(|(k, _)| k == k2) {
                    Some(t) => t.1 = on,
                    None => tags.push((k2.clone(), on)),
                }
            }
        }
    }
    tags.into_iter().filter(|(_, on)| *on).map(|(k, _)| k).collect()
}

/// Whether a section of a `RUN` block names a suffix keyword.
pub(crate) fn mentions(model: &Model) -> bool {
    [&model.bodies, &model.body_forces, &model.materials, &model.boundary_conditions]
        .iter()
        .flat_map(|v| v.iter())
        .any(|p| p.0.keys().any(|k| SUFFIXES.iter().any(|s| k.ends_with(s))))
}

/// Entity weights `[bodies, body forces, materials, BCs]`.
fn entity_weights(
    model: &Model,
    mesh: &dyn MeshView,
    coords: CoordinateSystem,
    fe: &FeBackend,
) -> Result<[Vec<f64>; 4], SolverError> {
    let mut w = [
        vec![0.0; model.bodies.len()],
        vec![0.0; model.body_forces.len()],
        vec![0.0; model.materials.len()],
        vec![0.0; model.boundary_conditions.len()],
    ];
    let mut ev = Default::default();
    for e in 0..(mesh.n_bulk() + mesh.n_boundary()) as u32 {
        let el = mesh.element(ElemId(e));
        let mut ids = [None; 4];
        if (e as usize) < mesh.n_bulk() {
            let Some(body) = model.body_of_tag(el.tag) else { continue };
            let section = model.bodies.get(body - 1);
            let id = |key: &str| section.and_then(|b| b.integer(key)).filter(|&i| i > 0).map(|i| i as usize);
            ids = [Some(body), id("body force"), id("material"), None];
        } else {
            let Some(bc) = model.bc_of_tag(el.tag) else { continue };
            ids[3] = Some(bc);
        }
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| mesh.coord(n)).collect();
        for q in fe.quadrature.gauss_points(el.etype, None) {
            fe.basis.eval(el.etype, &xyz, q.uvw, &mut ev).map_err(|m| SolverError::fatal("CalculateEntityWeights", m))?;
            let mut s = ev.det_j * q.weight;
            if coords.kind != CoordinateKind::Cartesian {
                let at = |c: usize| xyz.iter().zip(&ev.basis).map(|(p, b)| p[c] * b).sum::<f64>();
                s *= 2.0 * std::f64::consts::PI * coords.sqrt_metric(at(0), at(1), at(2));
            }
            for (k, id) in ids.iter().enumerate() {
                if let Some(slot) = id.and_then(|i| w[k].get_mut(i - 1)) {
                    *slot += s;
                }
            }
        }
    }
    Ok(w)
}

/// `SetNormalizedKeywords` on the run's mesh: each tagged keyword of a
/// body, body force, material or boundary condition times one over its
/// entity weight. Tables scale their values; other kinds are refused.
pub(crate) fn set_normalized_keywords(
    model: &mut Model,
    mesh: &dyn MeshView,
    coords: CoordinateSystem,
    fe: &FeBackend,
    ranks: usize,
) -> Result<(), SolverError> {
    let lists = [&model.bodies, &model.body_forces, &model.materials, &model.boundary_conditions];
    let tags: Vec<Vec<Vec<String>>> = lists.iter().map(|v| v.iter().map(tagged).collect()).collect();
    if tags.iter().flatten().all(Vec::is_empty) {
        return Ok(());
    }
    if ranks > 1 {
        return Err(SolverError::Unsupported("keyword normalisation by area or volume in parallel".into()));
    }
    if model.simulation.logical("update keyword normalization") == Some(true) {
        return Err(SolverError::Unsupported("Simulation: Update Keyword Normalization".into()));
    }
    let weights = entity_weights(model, mesh, coords, fe)?;
    let names = ["Body", "Body Force", "Material", "Boundary Condition"];
    let lists = [&mut model.bodies, &mut model.body_forces, &mut model.materials, &mut model.boundary_conditions];
    for (k, list) in lists.into_iter().enumerate() {
        for (i, p) in list.iter_mut().enumerate() {
            let w = weights[k][i];
            if tags[k][i].is_empty() || w <= f64::MIN_POSITIVE {
                continue;
            }
            for key in &tags[k][i] {
                let Some(v) = p.0.get_mut(key) else { continue };
                match v {
                    Value::Real { data, .. } => data.iter_mut().for_each(|x| *x /= w),
                    Value::Dependent(d) => match &mut d.kind {
                        elmer_core::DependentKind::Table { y, .. } => y.iter_mut().for_each(|x| *x /= w),
                        _ => {
                            return Err(SolverError::Unsupported(format!(
                                "{} {}: {key}: normalising a function-valued keyword by area or volume",
                                names[k],
                                i + 1
                            )))
                        }
                    },
                    _ => {
                        return Err(SolverError::Unsupported(format!(
                            "{} {}: {key}: normalising a non-real keyword by area or volume",
                            names[k],
                            i + 1
                        )))
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pm(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    #[test]
    fn suffix_keywords_tag_their_stems() {
        let one = Value::Real { rows: 1, cols: 1, data: vec![1.0] };
        let p = pm(&[
            ("heat source", one.clone()),
            ("heat source normalize by volume", Value::Logical(true)),
            ("density", one.clone()),
        ]);
        assert_eq!(tagged(&p), vec!["heat source".to_string()]);
        let p = pm(&[("heat flux", one), ("heat flux normalize by area", Value::Logical(false))]);
        assert!(tagged(&p).is_empty());
    }
}
