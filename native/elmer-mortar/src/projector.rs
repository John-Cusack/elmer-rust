//! Projector construction: `NodalProjector`, `PeriodicProjector` and
//! `GeneratePeriodicProjectors` (MortarUtils.F90).

use std::path::Path;

use elmer_core::{Logger, MeshView, Params};
use elmer_fe::LagrangeBasis;
use elmer_interp::{build_projector, Extras, InterpOptions, Locator};

use crate::interface::{map_interface_coordinate, InterfaceMesh};
use crate::overlay::overlay_interface_meshes;
use crate::save::save_projector;
use crate::weighted::weighted_projector;
use crate::MortarError;

const CALLER: &str = "PeriodicProjector";

/// A boundary projector (`Model % BCs(i) % PMatrix`, or a solver's
/// `MortarBCs(i) % Projector`): row `i` belongs to parent node
/// `inv_perm[i]` of this boundary and holds parent mesh nodes in
/// `cols[rows[i]..rows[i+1]]` with the weights in `values`.
///
/// A nodal projector ties the node to target-side nodes (columns in the
/// element order of the Fortran, not sorted; rows of nodes that were not
/// found are empty). A Galerkin projector ([`Projector::galerkin`]) holds
/// the weak form of `u_this - scale * u_target` tested with the node's basis
/// function, columns ascending.
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
    /// `ProjectorType == PROJECTOR_TYPE_GALERKIN`.
    pub galerkin: bool,
    /// `MortarBC % SlaveScale` and `MortarBC % MasterScale` (both 1 by
    /// default) when the projector is a constraint of the solver: a
    /// `Periodic BC Use Lagrange Coefficient` boundary gives `-Scale` and
    /// `-1`.
    pub slave_scale: f64,
    pub master_scale: f64,
    /// `MortarBC % Perm` is set (a Lagrange periodic boundary or a mortar
    /// jump): a Galerkin row then scales the columns of its own rows' nodes
    /// by `slave_scale`, the rest by `master_scale`, and sums the former
    /// into the weight of `rhs`.
    pub mortar_perm: bool,
    /// `MortarBC % Rhs` by row (`Mortar BC Offset`); the constraint row gets
    /// that weight times it.
    pub rhs: Option<Vec<f64>>,
    /// `MortarBC % Diag` by row (lumped, `MortarBC % LumpedDiag`): the
    /// constraint row gets `-diag * wsum` on its own multiplier.
    pub diag: Option<Vec<f64>>,
    /// `MortarBC % Active` (contact), `cdofs` per row: only the active
    /// components of a row become constraint rows. With it `rhs` is by
    /// component as well.
    pub active: Option<Vec<bool>>,
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
    /// `Mesh % DiscontPerm` when `CreateDiscontMesh` split the mesh
    /// (`Mesh % DisContMesh`): for each original node, the 1-based number of
    /// its copy, 0 if none. `None` otherwise.
    pub discont_perm: Option<&'a [u32]>,
    /// The edges and faces of the mesh (`Mesh % Edges`, `Mesh % Faces`)
    /// when the solver has edge dofs; `None` for a nodal mesh.
    pub edges: Option<EdgeEnv<'a>>,
}

/// `ElementInfo( Face, Nodes, u, v, w, detJ, Basis, EdgeBasis = WBasis,
/// USolver )` of the Piola basis at a p-reference point: the edge basis
/// and the nodal basis. Arguments: element code, nodes, global node
/// numbers, point.
pub type FaceEdgeBasis<'a> =
    &'a dyn Fn(u16, &[[f64; 3]], &[i64], [f64; 3]) -> Result<(Vec<[f64; 3]>, Vec<f64>), String>;

/// What `LevelProjector` needs for the edge dofs of a solver with an
/// H(curl) basis.
#[derive(Clone, Copy)]
pub struct EdgeEnv<'a> {
    /// `Find_Face`: the mesh face of a boundary element.
    pub face_of: crate::interface::FaceOf<'a>,
    /// `Mesh % NumberOfEdges`.
    pub n_edges: usize,
    /// `Mesh % MaxEdgeDOFs`.
    pub max_edge_dofs: usize,
    /// `EdgeElementStyle`: `PiolaVersion` and `SecondOrder`.
    pub piola: bool,
    pub second_order: bool,
    /// `isActivePElement` of the solver's first element.
    pub p_elements: bool,
    pub basis: FaceEdgeBasis<'a>,
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
        galerkin: false,
        slave_scale: 1.0,
        master_scale: 1.0,
        mortar_perm: false,
        rhs: None,
        diag: None,
        active: None,
    })
}

/// Projector families of `PeriodicProjector` that are not ported (the
/// radial, flat and plane ones are, see [`crate::level`]).
const REFUSED: [&str; 6] = [
    "normal projector",
    "axial projector",
    "anti axial projector",
    "sliding projector",
    "anti sliding projector",
    "stride projector",
];

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `PeriodicProjector(Model, Mesh, This, Trgt, Galerkin)` for the generic
/// projector, nodal ([`nodal_projector`]) or Galerkin
/// ([`weighted_projector`]). `galerkin` is the optional argument; without it
/// the boundary's `Galerkin Projector` decides. Returns `None` where the Fortran
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
    // A split mesh projects the discontinuous boundary onto its own copy; a
    // doubled one (no split) takes the ordinary path to its target.
    if flag(bc, "discontinuous boundary") {
        if let Some(perm) = env.discont_perm {
            let p = if weak {
                crate::weighted::weighted_projector_discont(
                    env.mesh,
                    env.bc_of,
                    this,
                    perm,
                    bc,
                    env.solver,
                    env.simulation,
                    env.coord_dim,
                    env.log,
                )?
            } else {
                Some(nodal_projector_discont(perm, this))
            };
            let Some(mut p) = p else { return Ok(None) };
            finish(env, bc, &mut p);
            return Ok(Some(p));
        }
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
    let Some((mut b1, mut b2)) = crate::interface::create_interface_meshes_faces(
        env.mesh,
        env.bc_of,
        env.bcs.len(),
        this,
        target,
        split,
        env.edges.as_ref().map(|e| e.face_of),
    )?
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
            "{k} (only the generic, radial, flat and plane projectors are ported)"
        )));
    }
    let anti_rotational = flag(bc, "anti rotational projector");
    let rotational = anti_rotational || flag(bc, "rotational projector");
    let anti_radial = flag(bc, "anti radial projector");
    let radial = anti_radial || flag(bc, "radial projector");
    let flat = flag(bc, "flat projector");
    let anti_plane = flag(bc, "anti plane projector");
    let plane = anti_plane || flag(bc, "plane projector");
    let cylindrical = flag(bc, "cylindrical projector");
    let family = rotational || cylindrical || radial || flat || plane;
    // `Mortar BC Scaling`, -1 by default for the anti projectors.
    let node_scale = bc
        .real("mortar bc scaling")
        .unwrap_or(if anti_radial || anti_plane { -1.0 } else { 1.0 });
    // These families default to the level projector ("dimensional
    // reduction"); an explicit `Level Projector = False` flattens z instead.
    let level_key = bc.logical("level projector");
    let mut level = level_key.unwrap_or(false);
    if family {
        match level_key {
            None => {
                env.log.info(
                    CALLER,
                    "Enforcing > Level Projector = True < with dimensional reduction",
                    7,
                );
                level = true;
            }
            Some(false) => {
                for c in b1.coords.iter_mut().chain(b2.coords.iter_mut()) {
                    c[2] = 0.0;
                }
            }
            Some(true) => {}
        }
    }
    let (mut radius, mut full_circle, mut periods) = (1.0, false, None);
    if rotational || cylindrical {
        (radius, full_circle, periods) =
            crate::level::rotational_interface_meshes(&mut b1, &mut b2, bc, env.simulation, cylindrical, env.log)?;
    } else if radial {
        crate::level::radial_interface_meshes(&mut b1, &mut b2, env.log);
    } else if flat {
        crate::level::flat_interface_meshes(&mut b1, &mut b2, bc, env.coord_dim, env.log);
    } else if plane {
        crate::level::plane_interface_meshes(&mut b1, &mut b2, bc, env.coord_dim, env.coord_dim, env.log)?;
    }
    if !rotational && bc.get("mesh rotate 3").is_some() {
        return Err(MortarError::Fatal(format!(
            "{CALLER}: Only \"Rotational Projector\" has the \"Mesh Rotate 3\" trick implemented"
        )));
    }

    // Not sliding or normal, so the overlay defaults to on, except for the
    // mapped families.
    if bc.logical("mortar bc enforce overlay").unwrap_or(!family) {
        overlay_interface_meshes(&b1, &mut b2, bc, env.log)?;
    }

    let mut p = if level {
        let solver_flag = |k: &str| env.solver.is_some_and(|s| flag(s, k));
        let do_nodes = !(solver_flag("projector skip nodes") || flag(bc, "projector skip nodes"));
        // DoEdges (MortarUtils.F90:8424-8442).
        let do_edges = !(solver_flag("projector skip edges") || flag(bc, "projector skip edges"))
            && env.edges.is_some_and(|e| e.n_edges > 0 && !e.p_elements)
            && env.mesh.dim() == 3
            && env.coord_dim == 3
            && solver_flag("hcurl basis");
        if !do_nodes && !do_edges {
            env.log.warn("LevelProjector", "Nothing to do, no nonodes, no edges!");
            return Ok(None);
        }
        let o = crate::level::LevelOptions {
            bc,
            simulation: env.simulation,
            solver: env.solver,
            mesh_dim: env.mesh.dim(),
            n_nodes: env.mesh.n_nodes(),
            coord_dim: env.coord_dim,
            node_scale,
            repeating: rotational && !full_circle,
            anti_repeating: rotational
                && !full_circle
                && bc.logical("antisymmetric bc").unwrap_or(anti_rotational),
            full_circle,
            radius,
            periods,
            do_nodes,
            edges: if do_edges { env.edges } else { None },
        };
        crate::level::level_projector(&b1, &b2, &o, env.log)?
    } else if weak {
        env.log.info(
            CALLER,
            "Weighted projector is very suboptimal, do not use it!",
            3,
        );
        env.log.warn(
            CALLER,
            "Use \"Normal Projector\" or \"Level Projector\" instead!",
        );
        let solver_flag = |k: &str| env.solver.is_some_and(|s| flag(s, k));
        if solver_flag("edge basis") {
            return Err(MortarError::Unsupported(
                "Galerkin projector with Edge Basis".into(),
            ));
        }
        let cs = env
            .simulation
            .string("coordinate system")
            .unwrap_or("")
            .to_ascii_lowercase();
        if cs.contains("symmetric") {
            // AxisSym weighting depends on Projector Metrics and the nodal
            // jump keywords.
            return Err(MortarError::Unsupported(
                "Galerkin projector in axisymmetric coordinates".into(),
            ));
        }
        let rel_order = env
            .solver
            .and_then(|s| s.integer("projector relative integration order"))
            .unwrap_or(0);
        if !(-1..=1).contains(&rel_order) {
            return Err(MortarError::Fatal(format!(
                "{CALLER}: Projector Relative Integration Order {rel_order} not in [-1, 1]"
            )));
        }
        let scale = node_scale;
        weighted_projector(
            &b2,
            &b1,
            env.interp_options(),
            env.coord_dim,
            rel_order,
            scale,
            env.log,
        )?
    } else {
        nodal_projector(&b2, &b1, env.interp_options(), env.coord_dim, env.log)?
    };

    p.projector_bc = this;
    finish(env, bc, &mut p);
    Ok(Some(p))
}

/// The end of `PeriodicProjector` (label 100): `Projector Set Rowsum` and
/// `Projector Multiplier`.
fn finish(env: &ProjectorEnv<'_>, bc: &dyn Params, p: &mut Projector) {
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
}

/// `NodalProjectorDiscont(Mesh, bc)`: row `j` ties the original node whose
/// copy is node `j` of `Mesh % DiscontPerm` (`InvPerm(j)`) to that copy,
/// node `n + j` (`n` the original node count), with weight 1.
pub fn nodal_projector_discont(discont_perm: &[u32], bc: usize) -> Projector {
    let n = discont_perm.len();
    let m = discont_perm.iter().filter(|&&j| j > 0).count();
    let mut cols = vec![0u32; m];
    let mut inv_perm = vec![elmer_core::NO_DOF; m];
    for (i, &j) in discont_perm.iter().enumerate() {
        if j == 0 {
            continue;
        }
        let j = j as usize - 1;
        cols[j] = (n + j) as u32;
        inv_perm[j] = i as u32;
    }
    Projector {
        rows: (0..=m).collect(),
        cols,
        values: vec![1.0; m],
        inv_perm,
        projector_bc: bc,
        galerkin: false,
        slave_scale: 1.0,
        master_scale: 1.0,
        mortar_perm: false,
        rhs: None,
        diag: None,
        active: None,
    }
}

/// The projectors of a model, `[i]` from boundary condition `i + 1`, and
/// whether a `Save Projector And Stop` fired (the Fortran stops the run
/// there).
pub type GeneratedProjectors = (Vec<Option<Projector>>, bool);

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
) -> Result<GeneratedProjectors, MortarError> {
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
        let stop = save_if_requested(bc, p.as_ref(), i, save_dir, timestep)?;
        out.push(p);
        if stop {
            return Ok((out, true));
        }
    }
    Ok((out, false))
}

/// The `Save Projector` step at the end of `PeriodicProjector`; true when
/// `Save Projector And Stop` fired.
fn save_if_requested(
    bc: &dyn Params,
    p: Option<&Projector>,
    this: usize,
    save_dir: &Path,
    timestep: i32,
) -> Result<bool, MortarError> {
    let (Some(p), Some(req)) = (p, save_request(bc, this, timestep)) else {
        return Ok(false);
    };
    save_projector(p, &save_dir.join(&req.prefix))
        .map_err(|e| MortarError::Fatal(format!("SaveProjector: {}: {e}", req.prefix)))?;
    Ok(req.and_stop)
}

/// `GenerateProjectors(Model, Solver, Nonlinear)` (ProjectorUtils.F90) for
/// `Mortar BC`: the Galerkin projector of every boundary condition with
/// `Mortar BC = j` onto boundary condition `j` (index `i` holds boundary
/// condition `i + 1`), saved into `save_dir` where `Save Projector` asks for
/// it. `env.solver` is the solver the projectors are for; `nonlinear`
/// selects the call made inside the nonlinear iteration, which only builds
/// the projectors with `Mortar BC Nonlinear`, the other call the rest.
///
/// Returns `None` when the solver applies no mortar conditions (no `Apply
/// Mortar BCs`). The flag is true when a `Save Projector And Stop` fired.
/// Contact, integral and shared (`Mortar BC Master Solver`) projectors are
/// refused.
pub fn generate_mortar_projectors(
    env: &ProjectorEnv<'_>,
    body_forces: &[&dyn Params],
    nonlinear: bool,
    save_dir: &Path,
    timestep: i32,
) -> Result<Option<GeneratedProjectors>, MortarError> {
    const CALLER: &str = "GenerateProjectors";
    let solver = env
        .solver
        .ok_or_else(|| MortarError::Fatal(format!("{CALLER}: no solver")))?;
    // `Apply Contact BCs` gates the projectors too; the contact itself is
    // [`crate::contact`]'s.
    let integral = flag(solver, "apply integral bcs");
    if !flag(solver, "apply mortar bcs") && !integral && !flag(solver, "apply contact bcs") {
        return Ok(None);
    }
    if solver.get("mortar bc master solver").is_some() {
        return Err(MortarError::Unsupported("Mortar BC Master Solver".into()));
    }
    env.log
        .info(CALLER, "Generating various boundary projectors", 8);
    let n_bcs = env.bcs.len();
    let mut out = Vec::with_capacity(n_bcs);
    for i in 1..=n_bcs + body_forces.len() {
        let (bc, body_force) = if i > n_bcs {
            (body_forces[i - n_bcs - 1], true)
        } else {
            (env.bc(i), false)
        };
        let mortar = bc.integer("mortar bc");
        let contact = bc.integer("contact bc");
        let integral_bc = flag(bc, "integral bc");
        let k = usize::from(mortar.is_some())
            + usize::from(contact.is_some())
            + usize::from(integral_bc);
        if k == 0 {
            out.push(None);
            continue;
        }
        if k > 1 {
            return Err(MortarError::Fatal(format!(
                "{CALLER}: Boundary {i} can only be one of mortar, contact and integral!"
            )));
        }
        if body_force && !integral_bc {
            return Err(MortarError::Fatal(format!(
                "{CALLER}: Body Force {i} can only have integral bc!"
            )));
        }
        if contact.is_some() {
            return Err(MortarError::Unsupported(format!(
                "Boundary {i}: Contact BC"
            )));
        }
        if integral_bc {
            return Err(MortarError::Unsupported(format!(
                "Boundary {i}: Integral BC"
            )));
        }
        if bc.logical("mortar bc nonlinear").unwrap_or(false) != nonlinear {
            out.push(None);
            continue;
        }
        let j = mortar.unwrap_or(0).max(0) as usize;
        let p = periodic_projector(env, i, j, Some(true))?;
        let stop = save_if_requested(bc, p.as_ref(), i, save_dir, timestep)?;
        out.push(p);
        if stop {
            return Ok(Some((out, true)));
        }
    }
    Ok(Some((out, false)))
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
            galerkin: false,
            slave_scale: 1.0,
            master_scale: 1.0,
            mortar_perm: false,
            rhs: None,
            diag: None,
            active: None,
        };
        p.set_rowsum();
        assert_eq!(p.values, vec![0.25, -0.75, 1.0]);
        assert_eq!(p.row(0).collect::<Vec<_>>(), vec![(4, 0.25), (5, -0.75)]);
    }

    /// Nodes 1 and 3 of four were doubled, as copies 2 and 1: rows follow
    /// the copies, each tying the original node to node `4 + copy - 1`.
    #[test]
    fn discont_projector_ties_each_node_to_its_copy() {
        let p = nodal_projector_discont(&[0, 2, 0, 1], 3);
        assert_eq!(p.rows, vec![0, 1, 2]);
        assert_eq!(p.cols, vec![4, 5]);
        assert_eq!(p.values, vec![1.0, 1.0]);
        assert_eq!(p.inv_perm, vec![3, 1]);
        assert_eq!(p.projector_bc, 3);
        assert!(!p.galerkin);
    }
}
