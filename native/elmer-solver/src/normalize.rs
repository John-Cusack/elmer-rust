//! `SetNormalizedKeywords` (MainUtils.F90:226): a keyword of a body, body
//! force, boundary condition or material whose section also sets
//! `<keyword> Normalize by Area` or `<keyword> Normalize by Volume` to True
//! is divided by the measure of that section's elements
//! (`CalculateEntityWeights`, SolverBasics.F90:10444), once.

use elmer_core::{normalize_key, ElemId, Logger, MeshView, ParamMap, Params, Value};

use crate::error::SolverError;
use crate::fortran_fmt::{es, i2s};
use crate::sim::Simulation;

const SUFFIXES: [&str; 2] = ["normalize by area", "normalize by volume"];

/// How many keys `ListTagEntry` tags in `p` with suffix `s`, whatever the
/// logical's value.
fn tag_count(p: &ParamMap, s: &str) -> usize {
    p.0.keys()
        .filter(|k| k.ends_with(s))
        .map(|k| p.0.keys().filter(|k2| k2.len() + s.len() <= k.len() && k.starts_with(k2.as_str())).count())
        .sum()
}

/// `ListTagEntry`: the keys a `... normalize by ...` key tags, every key
/// that is a prefix of it short of the suffix (`n2 + k <= n`); the tag is
/// the logical's value, the last tagging key winning.
fn tagged(p: &ParamMap) -> Vec<String> {
    let mut tags: std::collections::BTreeMap<&str, bool> = Default::default();
    for (k, v) in &p.0 {
        for s in SUFFIXES {
            if !k.ends_with(s) {
                continue;
            }
            let on = matches!(v, Value::Logical(true));
            for k2 in p.0.keys() {
                if k2.len() + s.len() <= k.len() && k.starts_with(k2.as_str()) {
                    tags.insert(k2, on);
                }
            }
        }
    }
    tags.into_iter().filter(|&(_, on)| on).map(|(k, _)| k.to_string()).collect()
}

/// Multiplies a keyword value by `c`: constants and the values of tables;
/// other value forms would need the coefficient at evaluation.
fn scale(v: &mut Value, c: f64) -> Result<(), SolverError> {
    match v {
        Value::Real { data, .. } => data.iter_mut().for_each(|x| *x *= c),
        Value::Dependent(d) => match &mut d.kind {
            elmer_core::DependentKind::Table { y, .. } => y.iter_mut().for_each(|x| *x *= c),
            _ => return Err(SolverError::Unsupported("keyword normalisation of a MATC, procedure or Lua value".into())),
        },
        _ => {}
    }
    Ok(())
}

impl Simulation {
    /// `SetNormalizedKeywords` on the first call; later calls do nothing
    /// (the weights are computed once, `Update Keyword Normalization` is
    /// refused).
    pub(crate) fn normalize_keywords(&mut self, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        if self.keywords_normalized {
            return Ok(());
        }
        self.keywords_normalized = true;
        let m = &self.model;
        // ListTagKeywords (ModelDescription.F90:2789) tags the keywords of
        // every section but the bodies.
        let lists: Vec<&ParamMap> = [&m.simulation, &m.constants]
            .into_iter()
            .chain(&m.equations)
            .chain(&m.body_forces)
            .chain(&m.initial_conditions)
            .chain(&m.boundary_conditions)
            .chain(&m.materials)
            .chain(&m.solvers)
            .collect();
        for s in SUFFIXES {
            let n: usize = lists.iter().map(|p| tag_count(p, s)).sum();
            if n > 0 {
                log.info("ListTagKeywords", &format!("Tagged {} parameters with suffix: {s}", i2s(n as i64)), 7);
            }
        }
        let groups: [&[ParamMap]; 3] = [&m.boundary_conditions, &m.materials, &m.body_forces];
        let tags: Vec<Vec<Vec<String>>> = groups.iter().map(|g| g.iter().map(tagged).collect()).collect();
        let count: usize = tags.iter().flatten().map(Vec::len).sum();
        if count == 0 {
            return Ok(());
        }
        log.info("ListTagCount", &format!("Found number of normalized keywords: {}", i2s(count as i64)), 6);
        if m.simulation.logical("update keyword normalization") == Some(true) {
            return Err(SolverError::Unsupported("Update Keyword Normalization".into()));
        }
        let cartesian = m
            .simulation
            .string("coordinate system")
            .is_none_or(|c| normalize_key(c).starts_with("cartesian"));
        if !cartesian {
            return Err(SolverError::Unsupported("keyword normalisation in a non-Cartesian coordinate system".into()));
        }

        // CalculateEntityWeights.
        log.info("CalculateEntityWeights", "Computing weights for the mesh entities", 6);
        let fe = self
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature for keyword normalisation (elmer-fe)".into()))?;
        // Body weights are left out: no body keyword is ever tagged.
        let (mut bc_w, mut mat_w, mut bf_w) = (
            vec![0.0; m.boundary_conditions.len()],
            vec![0.0; m.materials.len()],
            vec![0.0; m.body_forces.len()],
        );
        let nb = mesh.n_bulk();
        let mut ev = elmer_core::BasisEval::default();
        for e in (0..(nb + mesh.n_boundary()) as u32).map(ElemId) {
            let el = mesh.element(e);
            let (mut bc, mut bf, mut mat) = (None, None, None);
            if (e.0 as usize) < nb {
                let b = m.body_of_tag(el.tag).and_then(|b| m.bodies.get(b - 1));
                let id = |k: &str| b.and_then(|b| b.integer(k)).filter(|&i| i >= 1).map(|i| i as usize);
                bf = id("body force");
                mat = id("material");
            } else {
                bc = m.bc_of_tag(el.tag);
                if bc.is_none() {
                    continue;
                }
            }
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
            for qp in fe.quadrature.gauss_points(el.etype, None) {
                fe.basis
                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|msg| SolverError::fatal("CalculateEntityWeights", msg))?;
                let s = ev.det_j * qp.weight;
                for (w, i) in [(&mut bc_w, bc), (&mut mat_w, mat), (&mut bf_w, bf)] {
                    if let Some(x) = i.and_then(|i| w.get_mut(i - 1)) {
                        *x += s;
                    }
                }
            }
        }
        log.info("CalculateEntityWeights", "All done", 10);

        // ListSetParameters(partag = 0): body forces, BCs, materials.
        let model = &mut self.model;
        let mut scaled = 0;
        for (sections, tags, weights) in [
            (&mut model.body_forces, &tags[2], &bf_w),
            (&mut model.boundary_conditions, &tags[0], &bc_w),
            (&mut model.materials, &tags[1], &mat_w),
        ] {
            for ((p, keys), &w) in sections.iter_mut().zip(tags).zip(weights) {
                for k in keys {
                    if w > f64::MIN_POSITIVE {
                        let c = 1.0 / w;
                        if let Some(v) = p.0.get_mut(k) {
                            scale(v, c)?;
                        }
                        scaled += 1;
                        log.info("ListSetParameters", &format!("Scaling parameter \"{k}\" with:{}", es(c, 12, 3)), 15);
                    } else {
                        log.warn("ListSetParameters", "Refusing division with zero!");
                    }
                }
            }
        }
        if scaled > 0 {
            log.info("ListSetParameters", &format!("Scaled number of parameters: {}", i2s(scaled)), 6);
        } else {
            log.warn("ListSetParameters", "No parameters were altered!");
        }
        Ok(())
    }
}
