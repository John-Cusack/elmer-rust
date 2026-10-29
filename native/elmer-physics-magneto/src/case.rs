//! A loaded case: the `.sif` model as [`ParamMap`] sections and the mesh
//! after Elmer's load transforms (`Target Bodies` / `Target Boundaries`
//! mapping), read with `elmer-sif` and `elmer-mesh`.

use std::path::Path;

use elmer_core::{ParamMap, Params};
use elmer_sif::model::{Model, Section};

use crate::error::MagnetoError;

pub struct Case {
    pub mesh: elmer_mesh::Mesh,
    pub simulation: ParamMap,
    pub constants: ParamMap,
    /// Numbered sections, indexed by `id - 1` (`None` for a gap).
    pub bodies: Vec<Option<ParamMap>>,
    pub materials: Vec<Option<ParamMap>>,
    pub body_forces: Vec<Option<ParamMap>>,
    pub equations: Vec<Option<ParamMap>>,
    pub bcs: Vec<Option<ParamMap>>,
    pub solvers: Vec<Option<ParamMap>>,
    pub components: Vec<Option<ParamMap>>,
}

fn numbered(sections: &[Section]) -> Vec<Option<ParamMap>> {
    let n = sections.iter().map(|s| s.id as usize).max().unwrap_or(0);
    let mut out = vec![None; n];
    for s in sections {
        if s.id >= 1 {
            out[s.id as usize - 1] = Some(ParamMap::from(s));
        }
    }
    out
}

fn get(v: &[Option<ParamMap>], id: i32) -> Option<&ParamMap> {
    if id < 1 {
        return None;
    }
    v.get(id as usize - 1)?.as_ref()
}

impl Case {
    /// Reads `sif` and the mesh its `Header` names (relative to the
    /// directory of `sif`). Only the first `RUN` block is used; a case
    /// with later non-empty blocks is refused.
    pub fn load(sif: &Path) -> Result<Self, MagnetoError> {
        let path = sif.to_str().ok_or_else(|| MagnetoError::Input("non-UTF-8 case path".into()))?;
        let mut parsed = elmer_sif::parse_file(path).map_err(|e| MagnetoError::Input(e.to_string()))?;
        elmer_sif::complete(&mut parsed.sif, path).map_err(MagnetoError::Input)?;
        let model: Model = parsed
            .sif
            .runs
            .first()
            .cloned()
            .ok_or_else(|| MagnetoError::Input("empty .sif".into()))?;
        if parsed.sif.runs.iter().skip(1).any(|r| !r.solvers.is_empty() || !r.bcs.is_empty()) {
            return Err(MagnetoError::unsupported("a .sif with several RUN blocks"));
        }
        let (db, name) = model
            .header
            .mesh_db
            .clone()
            .ok_or_else(|| MagnetoError::Input("no Mesh DB in Header".into()))?;
        let dir = sif.parent().unwrap_or(Path::new(".")).join(db).join(name);
        let dir = dir.to_str().ok_or_else(|| MagnetoError::Input("non-UTF-8 mesh path".into()))?;
        let mut mesh = elmer_mesh::load(dir).map_err(|e| MagnetoError::Input(e.to_string()))?;
        elmer_mesh::transform(&mut mesh, &model).map_err(|e| MagnetoError::Input(e.to_string()))?;
        Ok(Self {
            mesh,
            simulation: model.simulation.as_ref().map(ParamMap::from).unwrap_or_default(),
            constants: model.constants.as_ref().map(ParamMap::from).unwrap_or_default(),
            bodies: numbered(&model.bodies),
            materials: numbered(&model.materials),
            body_forces: numbered(&model.body_forces),
            equations: numbered(&model.equations),
            bcs: numbered(&model.bcs),
            solvers: numbered(&model.solvers),
            components: numbered(&model.components),
        })
    }

    pub fn solver(&self, id: usize) -> Option<&ParamMap> {
        get(&self.solvers, id as i32)
    }

    pub fn body(&self, id: i32) -> Option<&ParamMap> {
        get(&self.bodies, id)
    }

    /// The section a body points at with `key` (`Material`, `Body Force`,
    /// `Equation`).
    fn body_section<'a>(&'a self, body: i32, key: &str, list: &'a [Option<ParamMap>]) -> Option<&'a ParamMap> {
        get(list, self.body(body)?.integer(key)?)
    }

    pub fn material(&self, body: i32) -> Option<&ParamMap> {
        self.body_section(body, "material", &self.materials)
    }

    pub fn body_force(&self, body: i32) -> Option<&ParamMap> {
        self.body_section(body, "body force", &self.body_forces)
    }

    pub fn equation(&self, body: i32) -> Option<&ParamMap> {
        self.body_section(body, "equation", &self.equations)
    }

    /// The boundary condition of a boundary element (its constraint after
    /// `Target Boundaries` mapping).
    pub fn bc(&self, constraint: i32) -> Option<&ParamMap> {
        get(&self.bcs, constraint)
    }

    /// Whether solver `id` is active in `body`: the body's equation lists
    /// it in `Active Solvers`, or sets the solver's `Equation` name true.
    pub fn solver_active_in(&self, id: usize, body: i32) -> bool {
        let Some(eq) = self.equation(body) else {
            return false;
        };
        if let Some(elmer_core::Value::Integer(v)) = eq.get("active solvers") {
            return v.contains(&(id as i32));
        }
        self.solver(id)
            .and_then(|s| s.string("equation"))
            .is_some_and(|name| eq.logical(name) == Some(true))
    }
}

/// `ListCheckPrefix`: some key of `p` starts with `prefix` (normalised).
pub fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    let prefix = elmer_core::normalize_key(prefix);
    p.0.keys().any(|k| k.starts_with(&prefix))
}
