//! Projector construction: `NodalProjector`, `PeriodicProjector` and
//! `GeneratePeriodicProjectors` (MortarUtils.F90).

use std::path::Path;

use elmer_core::{Logger, MeshView, Params};
use elmer_fe::LagrangeBasis;
use elmer_interp::{build_projector, Extras, InterpOptions, Locator};

use crate::interface::{create_interface_meshes, map_interface_coordinate, InterfaceMesh};
use crate::overlay::overlay_interface_meshes;
use crate::save::save_projector;
use crate::MortarError;

const CALLER: &str = "PeriodicProjector";

/// A boundary projector (`Model % BCs(i) % PMatrix`): row `i` ties parent
/// node `inv_perm[i]` of this boundary to the target-side parent nodes in
/// `cols[rows[i]..rows[i+1]]` with the weights in `values`. Columns keep the
/// element order of the Fortran (they are not sorted); rows of nodes that
/// were not found are empty.
#[derive(Clone, Debug, PartialEq)]
pub struct Projector {
    /// Length `n_rows + 1`.
    pub rows: Vec<usize>,
    /// Parent mesh nodes, 0-based.
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
    /// Row -> parent mesh node, 0-based; [`NO_DOF`](elmer_core::NO_DOF) for a
    /// row with no node.
    pub inv_perm: Vec<u32>,
    /// The boundary condition the projector belongs to (`ProjectorBC`), 1-based.
    pub projector_bc: usize,
}

impl Projector {
    pub fn n_rows(&self) -> usize {
        self.rows.len() - 1
    }

    /// Entries of row `i` as `(parent node, weight)`.
    pub fn row(&self, i: usize) -> impl Iterator<Item = (u32, f64)> + '_ {
        let r = self.rows[i]..self.rows[i + 1];
        self.cols[r.clone()]
            .iter()
            .copied()
            .zip(self.values[r].iter().copied())
    }

    /// `SetProjectorRowsum`: scale each row to unit absolute row sum.
    fn set_rowsum(&mut self) {
        for i in 0..self.n_rows() {
            let r = self.rows[i]..self.rows[i + 1];
            let mut s = 0.0;
            for v in &self.values[r.clone()] {
                s += v.abs();
            }
            for v in &mut self.values[r] {
                *v /= s;
            }
        }
    }
}

/// Model-level inputs of the projector routines.
pub struct ProjectorEnv<'a> {
    pub mesh: &'a dyn MeshView,
    /// `Boundary Condition 1..n` at indices `0..n`.
    pub bcs: &'a [&'a dyn Params],
    /// Boundary element tag -> 1-based boundary condition
    /// (`BoundaryInfo % Constraint`).
    pub bc_of: &'a dyn Fn(i32) -> Option<usize>,
    pub simulation: &'a dyn Params,
    /// `Model % Solver` while the projectors are built (`Model % Solvers(1)`
    /// in `GeneratePeriodicProjectors`), for the interpolation keywords.
    pub solver: Option<&'a dyn Params>,
    /// `CoordinateSystemDimension()`.
    pub coord_dim: usize,
    pub log: &'a dyn Logger,
}

impl ProjectorEnv<'_> {
    fn bc(&self, i: usize) -> &dyn Params {
        self.bcs[i - 1]
    }

    fn interp_options(&self) -> InterpOptions {
        let mut o = InterpOptions::from_params(self.simulation, self.solver);
        o.use_quadrant_tree = self.simulation.logical("use quadrant tree").unwrap_or(true);
        o.coordinate_system_dim = Some(self.coord_dim);
        o
    }
}

/// `NodalProjector` without the repeating (rotational/sliding) branch:
/// interpolate the nodes of `b1` in the elements of `b2` (`MeshProjector`)
/// and return to parent node numbering.
pub fn nodal_projector(
    b2: &InterfaceMesh,
    b1: &InterfaceMesh,
    opts: InterpOptions,
    coord_dim: usize,
    log: &dyn Logger,
) -> Result<Projector, MortarError> {
    let loc = Locator::new(b2, opts)?;
    let p = build_projector(
        &loc,
        b1,
        &LagrangeBasis::new(coord_dim),
        &Extras::default(),
        log,
    )?;
    Ok(Projector {
        rows: p.rows,
        cols: p.cols.iter().map(|&c| b2.inv_perm[c as usize]).collect(),
        values: p.values,
        inv_perm: b1.inv_perm.clone(),
        projector_bc: 0,
    })
}

/// Projector families of `PeriodicProjector` that are not ported.
const REFUSED: [&str; 15] = [
    "normal projector",
    "rotational projector",
    "anti rotational projector",
    "cylindrical projector",
    "radial projector",
    "anti radial projector",
    "axial projector",
    "anti axial projector",
    "sliding projector",
    "anti sliding projector",
    "flat projector",
    "plane projector",
    "anti plane projector",
    "level projector",
    "stride projector",
];

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `PeriodicProjector(Model, Mesh, This, Trgt, Galerkin)` for the generic
/// nodal projector. `galerkin` is the optional argument; without it the
/// boundary's `Galerkin Projector` decides. Returns `None` where the Fortran
/// returns a null projector (`This` or `Trgt` not positive, or an empty side).
/// Saving is left to the caller ([`save_request`], [`save_projector`]).
pub fn periodic_projector(
    env: &ProjectorEnv<'_>,
    this: usize,
    target: usize,
    galerkin: Option<bool>,
) -> Result<Option<Projector>, MortarError> {
    if this == 0 {
        return Ok(None);
    }
    let bc = env.bc(this);
    let weak = galerkin.unwrap_or_else(|| flag(bc, "galerkin projector"));
    if flag(bc, "discontinuous boundary") {
        return Err(MortarError::Unsupported(
            "Discontinuous Boundary projector".into(),
        ));
    }
    if target == 0 {
        return Ok(None);
    }
    env.log.info(
        CALLER,
        &format!("Creating projector between BCs {this} and {target}"),
        8,
    );

    let split = env
        .simulation
        .logical("mortar bcs split quadratic")
        .unwrap_or(false);
    let Some((mut b1, mut b2)) = create_interface_meshes(env.mesh, env.bc_of, this, target, split)?
    else {
        return Ok(None);
    };
    let map: Option<Vec<i32>> = match bc.get("projector coordinate mapping") {
        Some(elmer_core::Value::Integer(v)) => Some(v.clone()),
        _ => None,
    };
    map_interface_coordinate(&mut b1, &mut b2, map.as_deref())?;

    if let Some(k) = REFUSED.iter().find(|k| flag(bc, k)) {
        return Err(MortarError::Unsupported(format!(
            "{k} (only the generic nodal projector is ported)"
        )));
    }
    if bc.get("mesh rotate 3").is_some() {
        return Err(MortarError::Fatal(format!(
            "{CALLER}: Only \"Rotational Projector\" has the \"Mesh Rotate 3\" trick implemented"
        )));
    }

    // Not sliding or normal, so the overlay defaults to on.
    if bc.logical("mortar bc enforce overlay").unwrap_or(true) {
        overlay_interface_meshes(&b1, &mut b2, bc, env.log)?;
    }

    if weak {
        return Err(MortarError::Unsupported(
            "Galerkin Projector (WeightedProjector)".into(),
        ));
    }
    let mut p = nodal_projector(&b2, &b1, env.interp_options(), env.coord_dim, env.log)?;

    p.projector_bc = this;
    if flag(bc, "projector set rowsum") {
        p.set_rowsum();
    }
    if let Some(c) = bc
        .real("projector multiplier")
        .or_else(|| env.simulation.real("projector multiplier"))
    {
        for v in &mut p.values {
            *v *= c;
        }
    }
    env.log.info(CALLER, "Projector created, now exiting...", 8);
    Ok(Some(p))
}

/// What `Save Projector` asks for on one boundary condition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SaveRequest {
    /// File prefix: `p<This>` or, with `Save Projector Timestep Numbering`,
    /// `p<This>_<timestep:04>`.
    pub prefix: String,
    /// `Save Projector And Stop`: the Fortran exits right after saving.
    pub and_stop: bool,
}

/// `Save Projector` keywords of boundary condition `this`; `timestep` is the
/// `timestep` variable (0 when absent).
pub fn save_request(bc: &dyn Params, this: usize, timestep: i32) -> Option<SaveRequest> {
    if !flag(bc, "save projector") {
        return None;
    }
    let mut prefix = format!("p{this}");
    if flag(bc, "save projector timestep numbering") {
        prefix = format!("{prefix}_{timestep:04}");
    }
    Some(SaveRequest {
        prefix,
        and_stop: flag(bc, "save projector and stop"),
    })
}

/// `GeneratePeriodicProjectors`: the projector of every boundary condition
/// with `Periodic BC` (index `i` holds boundary condition `i + 1`), saved
/// into `save_dir` where `Save Projector` asks for it. The flag is true when
/// a `Save Projector And Stop` fired; the Fortran stops the run there.
pub fn generate_periodic_projectors(
    env: &ProjectorEnv<'_>,
    save_dir: &Path,
    timestep: i32,
) -> Result<(Vec<Option<Projector>>, bool), MortarError> {
    if env.bcs.iter().any(|bc| bc.get("conforming bc").is_some()) {
        return Err(MortarError::Unsupported("Conforming BC".into()));
    }
    let mut out = Vec::with_capacity(env.bcs.len());
    for i in 1..=env.bcs.len() {
        let bc = env.bc(i);
        let Some(k) = bc.integer("periodic bc") else {
            out.push(None);
            continue;
        };
        let p = periodic_projector(env, i, k.max(0) as usize, None)?;
        if let (Some(p), Some(req)) = (&p, save_request(bc, i, timestep)) {
            save_projector(p, &save_dir.join(&req.prefix))
                .map_err(|e| MortarError::Fatal(format!("SaveProjector: {}: {e}", req.prefix)))?;
            if req.and_stop {
                out.push(Some(p.clone()));
                return Ok((out, true));
            }
        }
        out.push(p);
    }
    Ok((out, false))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rowsum_scales_to_unit_abs_sum() {
        let mut p = Projector {
            rows: vec![0, 2, 3],
            cols: vec![4, 5, 6],
            values: vec![1.0, -3.0, 0.5],
            inv_perm: vec![0, 1],
            projector_bc: 1,
        };
        p.set_rowsum();
        assert_eq!(p.values, vec![0.25, -0.75, 1.0]);
        assert_eq!(p.row(0).collect::<Vec<_>>(), vec![(4, 0.25), (5, -0.75)]);
    }
}
