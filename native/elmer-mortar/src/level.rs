//! The level projector of a 2D mesh (`LevelProjector` with
//! `AddProjectorWeak1D` and `TemporalSegmentMortarAssembly`, MortarUtils.F90)
//! and the interface-mesh maps that lead to it: `FlatInterfaceMeshes`,
//! `PlaneInterfaceMeshes` and `RadialInterfaceMeshes`.
//!
//! After the map both interface meshes are lines along `x`. Each slave
//! (`this`) element is cut by every master (`target`) element it overlaps
//! in `x`, and each common segment is integrated with the line's
//! `GaussPoints2` rule: row `nrow` of slave node `j` gets
//! `N_j N_i w` on the slave nodes `i` and `-scale N_j M_k w` on the master
//! nodes `k`.
//!
//! Ported: the Galerkin (weak) projector for nodal dofs of a 2D mesh, with
//! `Projector Node Multiplier`, `Projector Max Distance`, `Mortar BC
//! Scaling` and `Level Projector Eliminate Redundant Nodes`. Refused: the
//! strong (node-to-segment) projector, edge dofs, p-elements, biorthogonal
//! or dual bases, Nitsche, `Level Projector Condition` and the 3D surface
//! integrator (`AddProjectorWeakGeneric`).

use std::collections::BTreeMap;

use elmer_core::{BasisEval, ElemId, ElementType, Logger, MeshView, Params, Value};
use elmer_fe::{element_info, element_type, gauss_points, MetricStatus};
use elmer_interp::geom::global_to_local;
use elmer_interp::shape::Shape;

use crate::interface::InterfaceMesh;
use crate::projector::Projector;
use crate::MortarError;

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

fn fe(e: elmer_fe::FeError) -> MortarError {
    MortarError::Fatal(format!("LevelProjector: {e}"))
}

/// `FlatInterfaceMeshes`: the coordinate the interface is flat in becomes
/// `z`, and in 2D `y` is zeroed. Without `Flat Projector Reduce Dimension`
/// the coordinate is detected anew (the smallest relative extent over both
/// meshes), even when `Flat Projector Coordinate` is given: the Fortran's
/// `Found` is that of the later keyword. Returns the coordinate used
/// (1-based; 0 when none is given and none detected), which the Fortran
/// stores back into the boundary condition.
pub fn flat_interface_meshes(
    b1: &mut InterfaceMesh,
    b2: &mut InterfaceMesh,
    bc: &dyn Params,
    model_dim: usize,
    log: &dyn Logger,
) -> usize {
    log.info("FlatInterfaceMeshes", "Flattening interface meshes to 2D", 8);
    let mut flat = bc.integer("flat projector coordinate").unwrap_or(0).clamp(0, 3) as usize;
    let reduce = bc.logical("flat projector reduce dimension");
    if reduce.is_none() {
        let mut rel1 = 0.0;
        for (j, m) in [&*b1, &*b2].into_iter().enumerate() {
            let (mut max_diff, mut min_diff, mut min_i) = (0.0f64, f64::MAX, 0usize);
            for i in 0..model_dim {
                let mx = m.coords.iter().map(|c| c[i]).fold(f64::NEG_INFINITY, f64::max);
                let mn = m.coords.iter().map(|c| c[i]).fold(f64::INFINITY, f64::min);
                let diff = mx - mn;
                max_diff = max_diff.max(diff);
                if diff < min_diff {
                    min_diff = diff;
                    min_i = i + 1;
                }
            }
            let rel = min_diff / max_diff;
            if j == 0 {
                flat = min_i;
                rel1 = rel;
            } else if rel < rel1 {
                flat = min_i;
            }
        }
        log.info("FlatInterfaceMeshes", &format!("> Flat Projector Coordinate < set to: {flat}"), 4);
    }
    for m in [b1, b2] {
        for c in m.coords.iter_mut() {
            match flat {
                2 => {
                    *c = [c[0], c[2], c[1]];
                    if model_dim == 2 {
                        c[1] = 0.0;
                    }
                }
                1 => {
                    *c = [c[1], c[2], c[0]];
                    if model_dim == 2 {
                        c[1] = 0.0;
                    }
                }
                _ => {}
            }
            if reduce == Some(true) {
                c[2] = 0.0;
            }
        }
    }
    flat
}

/// `NormalVector( Element, Nodes, u, v, .FALSE. )` of a line element.
fn line_normal(code: u16, x: &[[f64; 3]], u: f64) -> Result<[f64; 3], MortarError> {
    let d: Vec<f64> = match code {
        202 => vec![-0.5, 0.5],
        203 => vec![u - 0.5, u + 0.5, -2.0 * u],
        t => {
            return Err(MortarError::Unsupported(format!(
                "Plane Projector Normal detection on boundary elements of type {t}"
            )))
        }
    };
    let dot = |a: usize| d.iter().zip(x).fold(0.0, |s, (di, xi)| s + di * xi[a]);
    let (dxdu, dydu) = (dot(0), dot(1));
    let det = dxdu * dxdu + dydu * dydu;
    if det <= 0.0 {
        return Ok([0.0; 3]);
    }
    let det = 1.0 / det.sqrt();
    Ok([-dydu * det, dxdu * det, 0.0])
}

/// `PlaneInterfaceMeshes`: coordinates in the frame `(t1, t2, n)` of the
/// plane normal, `y` zeroed in 2D. Without `Plane Projector Normal` the
/// normal is the area-weighted mean normal of the flatter mesh (the master
/// one negated), which the Fortran stores back into the boundary condition;
/// it is returned.
pub fn plane_interface_meshes(
    b1: &mut InterfaceMesh,
    b2: &mut InterfaceMesh,
    bc: &dyn Params,
    model_dim: usize,
    coord_dim: usize,
    log: &dyn Logger,
) -> Result<[f64; 3], MortarError> {
    log.info("PlaneInterfaceMeshes", "Flattening interface meshes to a plane", 8);
    let given = match bc.get("plane projector normal") {
        Some(Value::Real { data, .. }) if data.len() >= 3 => Some([data[0], data[1], data[2]]),
        _ => None,
    };
    let normal = match given {
        Some(n) => n,
        None => {
            if model_dim == 3 {
                return Err(MortarError::Unsupported(
                    "Plane Projector Normal detection in 3D".into(),
                ));
            }
            let mut st = MetricStatus::default();
            let mut be = BasisEval::default();
            let (mut pn1, mut planeness1) = ([0.0; 3], 0.0);
            let (mut pn, mut planeness) = ([0.0; 3], 0.0);
            for (j, m) in [&*b1, &*b2].into_iter().enumerate() {
                let mut sum = [0.0f64; 3];
                let mut refsum = 0.0f64;
                let mut normal0: Option<[f64; 3]> = None;
                for e in 0..m.n_bulk() {
                    let el = m.element(ElemId(e as u32));
                    let def = element_type(el.etype.0).map_err(fe)?;
                    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| m.coords[n as usize]).collect();
                    let rule = gauss_points(el.etype, None).map_err(fe)?;
                    for k in 0..rule.s.len() {
                        element_info(def, &xyz, [rule.u[k], rule.v[k], rule.w[k]], coord_dim, &mut st, &mut be)
                            .map_err(fe)?;
                        let mut n = line_normal(el.etype.0, &xyz, rule.u[k])?;
                        let n0 = *normal0.get_or_insert(n);
                        if n[0] * n0[0] + n[1] * n0[1] + n[2] * n0[2] < 0.0 {
                            n = n.map(|x| -x);
                        }
                        for a in 0..3 {
                            sum[a] += rule.s[k] * be.det_j * n[a];
                        }
                        refsum += rule.s[k] * be.det_j;
                    }
                }
                let length = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
                pn = sum.map(|x| x / length);
                planeness = length / refsum;
                if j == 0 {
                    pn1 = pn;
                    planeness1 = planeness;
                }
            }
            // An exact comparison, as written.
            let n = if planeness1 > planeness {
                log.info("PlaneInterfaceMeshes", "Selecting slave normal", 12);
                pn1
            } else {
                log.info("PlaneInterfaceMeshes", "Selecting master normal", 12);
                pn.map(|x| -x)
            };
            log.info(
                "PlaneInterfaceMeshes",
                &format!("PlaneNormal selected:{:8.3}{:8.3}{:8.3}", n[0], n[1], n[2]),
                8,
            );
            n
        }
    };
    let (t1, t2) = elmer_solver::tangent_directions(normal);
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    for m in [b1, b2] {
        for c in m.coords.iter_mut() {
            let x = *c;
            *c = [dot(&x, &t1), if model_dim == 3 { dot(&x, &t2) } else { 0.0 }, dot(&x, &normal)];
        }
    }
    Ok(normal)
}

/// `RadialInterfaceMeshes`: `(x, y, z) -> (r, z, 0)`.
pub fn radial_interface_meshes(b1: &mut InterfaceMesh, b2: &mut InterfaceMesh, log: &dyn Logger) {
    let mut ext = [[0.0f64; 2]; 2];
    for (k, m) in [b1, b2].into_iter().enumerate() {
        let (mut rmin, mut rmax) = (f64::MAX, f64::MIN);
        for c in m.coords.iter_mut() {
            let r = (c[0] * c[0] + c[1] * c[1]).sqrt();
            let z = c[2];
            *c = [r, z, 0.0];
            rmin = rmin.min(r);
            rmax = rmax.max(r);
        }
        ext[k] = [rmin, rmax];
    }
    let rad = ext[0][1].max(ext[1][1]);
    let err1 = (ext[0][1] - ext[1][1]).abs() / rad;
    let err2 = (ext[0][0] - ext[1][0]).abs() / rad;
    if err1 > 1.0e-3 || err2 > 1.0e-3 {
        log.warn("RadialInterfaceMeshes", "Discrepancy of radius may be too large!");
    }
}

/// `RotationalInterfaceMeshes` (MortarUtils.F90:7272), not cylindrical:
/// each mesh to `(phi, z, r)` with `phi` in degrees, relative to `Cylinder
/// Center` and in the frame of `Cylinder Normal` when given, the slave one
/// offset by `Rotational Projector Angle Offset` (or the boundary's `Mesh
/// Rotate 3`). A mesh crossing 180 degrees is moved to a continuous angle
/// range; one covering all four quadrants is a full circle. A constant `z`
/// is flattened. Returns the radius (the slave's largest), whether the
/// interface is a full circle and, if not, the number of periods (given, or
/// `360 / dfii` of the master when that is an integer). The Fortran also
/// stores `Rotor Radius` and `Rotor Periods` in the Simulation section; they
/// are not kept here.
///
/// `Cylindrical` (the cylindrical projector): whether the slave is a full
/// circle, and the move of its angles past 180 degrees, decide for the
/// master too, and there is no number of periods. `CylinderFit`, which the
/// Fortran calls for it and for `Rotational Projector Center Fit`, stores
/// `Cylinder Center` and `Cylinder Normal` as arrays, which the `Cylinder
/// Center X` ... read here never finds: it changes nothing and is skipped.
pub fn rotational_interface_meshes(
    b1: &mut InterfaceMesh,
    b2: &mut InterfaceMesh,
    bc: &dyn Params,
    simulation: &dyn Params,
    cylindrical: bool,
    log: &dyn Logger,
) -> Result<(f64, bool, Option<i32>), MortarError> {
    const CALLER: &str = "RotationalInterfaceMeshes";
    let rad2deg = 180.0 / std::f64::consts::PI;
    let comp = |k: &str| bc.real(k);
    let cx = [comp("cylinder center x"), comp("cylinder center y"), comp("cylinder center z")];
    let got_center = cx.iter().any(Option::is_some);
    let x0 = cx.map(|v| v.unwrap_or(0.0));
    let nx = [comp("cylinder normal x"), comp("cylinder normal y"), comp("cylinder normal z")];
    let got_normal = nx.iter().any(Option::is_some);
    let normal = nx.map(|v| v.unwrap_or(0.0));
    let (t1, t2) = if got_normal { elmer_solver::tangent_directions(normal) } else { ([0.0; 3], [0.0; 3]) };
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let es = |v: f64| elmer_solver::fortran_fmt::es(v, 12, 3);
    let mut full_circle = false;
    let (mut move_angle, mut fii0) = (false, 0.0f64);
    let mut ext = [[[0.0f64; 3]; 2]; 2];
    for (k, m) in [&mut *b1, &mut *b2].into_iter().enumerate() {
        let mut offset = None;
        if k == 0 {
            offset = bc.real("rotational projector angle offset").or_else(|| bc.real("mesh rotate 3"));
            if offset.is_some() && simulation.logical("rotate in radians") == Some(true) {
                offset = offset.map(|d| d * rad2deg);
            }
        }
        let which = if k == 0 { "this" } else { "target" };
        log.info(CALLER, &format!("Initial extrema for {which} boundary (x,y,z)"), 8);
        for i in 0..3 {
            let mn = m.coords.iter().map(|c| c[i]).fold(f64::INFINITY, f64::min);
            let mx = m.coords.iter().map(|c| c[i]).fold(f64::NEG_INFINITY, f64::max);
            log.info(CALLER, &format!("Coordinate {}: {}{}", i + 1, es(mn), es(mx)), 8);
        }
        for c in m.coords.iter_mut() {
            let mut x = *c;
            if got_center {
                for a in 0..3 {
                    x[a] -= x0[a];
                }
            }
            if got_normal {
                let t = x;
                x = [dot(&t1, &t), dot(&t2, &t), dot(&normal, &t)];
            }
            let mut alpha = rad2deg * x[1].atan2(x[0]);
            let rad = (x[0] * x[0] + x[1] * x[1]).sqrt();
            if let Some(d) = offset {
                alpha = (alpha + d).rem_euclid(360.0);
                if alpha > 180.0 {
                    alpha -= 360.0;
                }
            }
            *c = [alpha, x[2], rad];
        }
        if cylindrical && k == 1 {
            if move_angle {
                for c in m.coords.iter_mut() {
                    if c[0] < fii0 {
                        c[0] += 360.0;
                    }
                }
            }
        } else {
        let (mut hit0, mut hit90, mut hit180, mut hit270) = (false, false, false, false);
        for e in 0..m.n_bulk() {
            let el = m.element(ElemId(e as u32));
            let fmin = el.nodes.iter().map(|&n| m.coords[n as usize][0]).fold(f64::INFINITY, f64::min);
            let fmax = el.nodes.iter().map(|&n| m.coords[n as usize][0]).fold(f64::NEG_INFINITY, f64::max);
            if fmax - fmin > 180.0 {
                hit180 = true;
            } else {
                hit0 |= fmax >= 0.0 && fmin <= 0.0;
                hit90 |= fmax >= 90.0 && fmin <= 90.0;
                hit270 |= fmax >= -90.0 && fmin <= -90.0;
            }
        }
        full_circle = hit0 && hit90 && hit180 && hit270;
        if full_circle {
            log.info(CALLER, "Cylindrical interface seems to be a full circle", 6);
        } else if hit180 {
            let fii = if !hit0 {
                0.0
            } else if !hit270 {
                -90.0
            } else if !hit90 {
                90.0
            } else {
                0.0
            };
            for c in m.coords.iter_mut() {
                if c[0] < fii {
                    c[0] += 360.0;
                }
            }
            move_angle = true;
            fii0 = fii;
        }
        }
        log.info(CALLER, &format!("Transformed extrema for {which} boundary (phi,z,r)"), 8);
        for i in 0..3 {
            let mn = m.coords.iter().map(|c| c[i]).fold(f64::INFINITY, f64::min);
            let mx = m.coords.iter().map(|c| c[i]).fold(f64::NEG_INFINITY, f64::max);
            ext[k][0][i] = mn;
            ext[k][1][i] = mx;
            log.info(CALLER, &format!("Coordinate {}: {}{}", i + 1, es(mn), es(mx)), 8);
        }
        if ext[k][0][2] < f64::EPSILON {
            return Err(MortarError::Fatal(format!("{CALLER}: Radius cannot be almost zero!")));
        }
    }
    let eps_rad = 1.0e-3;
    let radius = ext[0][1][2];
    let err1 = (ext[0][1][2] - ext[0][0][2]) / radius;
    let err2 = (ext[1][1][2] - ext[1][0][2]) / radius;
    log.info(CALLER, &format!("Radius of the rotational interface:{}", es(radius)), 8);
    log.info(CALLER, &format!("Discrepancy from constant radius for Mesh1:{}", es(err1)), 8);
    if err1 > eps_rad {
        log.warn(CALLER, "Discrepancy of radius is rather large!");
    }
    log.info(CALLER, &format!("Discrepancy from constant radius for Mesh2:{}", es(err2)), 8);
    if err2 > eps_rad {
        log.warn(CALLER, "Discrepancy of radius is rather large!");
    }
    let e1 = (ext[0][1][1] - ext[0][0][1]) / radius;
    let e2 = (ext[1][1][1] - ext[1][0][1]) / radius;
    if e1 < eps_rad && e2 < eps_rad {
        log.info(CALLER, "The effective interface meshes are 1D", 8);
        for c in b1.coords.iter_mut().chain(b2.coords.iter_mut()) {
            c[1] = 0.0;
        }
    } else {
        log.info(CALLER, "The effective interface meshes are 2D", 8);
    }
    if cylindrical || full_circle {
        return Ok((radius, full_circle, None));
    }
    let dfii1 = ext[0][1][0] - ext[0][0][0];
    let dfii2 = ext[1][1][0] - ext[1][0][0];
    log.info(CALLER, &format!("This boundary dfii:  {}", es(dfii1)), 8);
    log.info(CALLER, &format!("Target boundary dfii:  {}", es(dfii2)), 8);
    log.info(CALLER, &format!("Discrepancy in dfii:{}", es(2.0 * (dfii1 - dfii2).abs() / (dfii1 + dfii2))), 8);
    let nsym = 360.0 / dfii2;
    let periods = match bc.integer("rotational projector periods") {
        Some(i) => {
            log.info(CALLER, &format!("Using enforced number of periods: {i}"), 8);
            log.info(CALLER, &format!("Suggested number of periods:{}", es(nsym)), 8);
            i
        }
        None => {
            log.info(CALLER, &format!("Suggested sections in target:{}", es(nsym)), 8);
            if (nsym - nsym.round()).abs() < 0.01 || nsym < 1.5 {
                nsym.round() as i32
            } else {
                return Err(MortarError::Fatal(format!("{CALLER}: Check your settings, this cannot be periodic!")));
            }
        }
    };
    Ok((radius, false, Some(periods)))
}

/// What `LevelProjector` reads besides the meshes.
pub struct LevelOptions<'a> {
    pub bc: &'a dyn Params,
    pub simulation: &'a dyn Params,
    pub solver: Option<&'a dyn Params>,
    /// `Mesh % MeshDim` of the parent mesh.
    pub mesh_dim: usize,
    /// Nodes of the parent mesh (`NodePerm` size).
    pub n_nodes: usize,
    /// `CoordinateSystemDimension()`.
    pub coord_dim: usize,
    /// `NodeScale` (`Mortar BC Scaling`, else -1 for the anti projectors).
    pub node_scale: f64,
    /// `Repeating`, `AntiRepeating`, `FullCircle` and `Radius` of
    /// `PeriodicProjector`, and the rotational interface's number of periods.
    pub repeating: bool,
    pub anti_repeating: bool,
    pub full_circle: bool,
    pub radius: f64,
    pub periods: Option<i32>,
}

/// `LevelProjector( BMesh1, BMesh2, .FALSE., .FALSE., .FALSE., 1, DoNodes,
/// .FALSE., NodeScale, NodeScale, BC )` of a 2D mesh with a Galerkin
/// projector (`AddProjectorWeak1D`).
pub fn level_projector(
    b1: &InterfaceMesh,
    b2: &InterfaceMesh,
    o: &LevelOptions<'_>,
    log: &dyn Logger,
) -> Result<Projector, MortarError> {
    const CALLER: &str = "LevelProjector";
    log.info(CALLER, "Creating projector for a levelized mesh", 7);
    let bc = o.bc;
    let real = |k: &str| bc.real(k).or_else(|| o.simulation.real(k));
    let node_coeff = real("projector node multiplier").unwrap_or(1.0);
    let max_distance = real("projector max distance");
    if flag(bc, "nitsche projector") {
        return Err(MortarError::Unsupported("Nitsche Projector".into()));
    }
    // WeakProjector / StrongProjector.
    let weak = match bc.logical("galerkin projector") {
        Some(g) => g,
        None => !bc.logical("level projector strong").unwrap_or(true),
    };
    let strong_nodes = bc.logical("level projector nodes strong").unwrap_or(!weak);
    if strong_nodes || flag(bc, "level projector conforming nodes strong") {
        return Err(MortarError::Unsupported("the strong (node-to-segment) level projector".into()));
    }
    let solver_flag = |k: &str| o.solver.is_some_and(|s| flag(s, k));
    let biorth = flag(bc, "use biorthogonal basis")
        || (bc.get("create dual projector").is_none() && solver_flag("eliminate linear constraints"));
    if biorth || flag(bc, "create dual projector") {
        return Err(MortarError::Unsupported("biorthogonal or dual level projectors".into()));
    }
    if bc.get("level projector condition").is_some() {
        return Err(MortarError::Unsupported("Level Projector Condition".into()));
    }
    let mut mesh_dim = o.mesh_dim;
    if mesh_dim == 3 && b1.n_bulk() > 0 && b1.element(ElemId(0)).etype.family() == 2 {
        log.warn(CALLER, "Enforcing 1D integration for 1D boundary elements in 3D mesh!");
        mesh_dim = 2;
    }
    let rotational = flag(bc, "rotational projector") || flag(bc, "anti rotational projector");
    let cylindrical = flag(bc, "cylindrical projector");
    if mesh_dim != 2 && mesh_dim != 3 {
        return Err(MortarError::Unsupported(format!("the level projector of a {mesh_dim}D mesh")));
    }
    if mesh_dim == 2 && (o.repeating || o.full_circle || rotational || cylindrical) {
        return Err(MortarError::Unsupported("a repeating or rotational level projector of a 2D mesh".into()));
    }
    if mesh_dim == 3 {
        if bc.integer("axial projector periods").is_some_and(|n| n > 1) {
            return Err(MortarError::Unsupported("Axial Projector Periods".into()));
        }
        if flag(bc, "projector symmetry x") || flag(bc, "projector symmetry y") {
            return Err(MortarError::Unsupported("Projector Symmetry".into()));
        }
    }

    // NodePerm: slave nodes, less those on the master side.
    let mut node_perm = vec![0usize; o.n_nodes];
    for e in 0..b1.n_bulk() {
        for &n in b1.element(ElemId(e as u32)).nodes {
            node_perm[b1.inv_perm[n as usize] as usize] = 1;
        }
    }
    if bc.logical("level projector eliminate redundant nodes").unwrap_or(true) {
        let mut m = 0;
        for e in 0..b2.n_bulk() {
            for &n in b2.element(ElemId(e as u32)).nodes {
                let k = b2.inv_perm[n as usize] as usize;
                if node_perm[k] != 0 {
                    node_perm[k] = 0;
                    m += 1;
                }
            }
        }
        if m > 0 {
            log.info(CALLER, &format!("Eliminating redundant nodes from projector: {m}"), 10);
        }
    }
    let mut m = 0usize;
    for p in node_perm.iter_mut() {
        if *p > 0 {
            m += 1;
            *p = m;
        }
    }
    log.info(CALLER, &format!("Number of active nodes in projector: {m}"), 8);
    let mut inv_perm = vec![elmer_core::NO_DOF; m];
    let mut list: Vec<BTreeMap<u32, f64>> = vec![BTreeMap::new(); m];

    if mesh_dim == 2 {
        weak_1d(b1, b2, o, node_coeff, max_distance, &node_perm, &mut inv_perm, &mut list, log)?;
    } else {
        // XTol and YTol from the master's extent, before a rotational
        // interface's XRange becomes the period.
        let ext = |a: usize| {
            let mx = b2.coords.iter().map(|c| c[a]).fold(f64::NEG_INFINITY, f64::max);
            let mn = b2.coords.iter().map(|c| c[a]).fold(f64::INFINITY, f64::min);
            mx - mn
        };
        let (mut xrange, yrange) = (ext(0), ext(1));
        let rel_tol_y = 1.0e-4;
        let rel_tol_x = if o.full_circle || o.repeating { 1.0e-6 } else { rel_tol_y };
        let (ytol, xtol) = (rel_tol_y * yrange, rel_tol_x * xrange);
        let arc_coeff =
            if rotational || cylindrical { (2.0 * std::f64::consts::PI * o.radius) / 360.0 } else { 1.0 };
        if rotational {
            xrange = if o.full_circle {
                360.0
            } else {
                let p = o.periods.ok_or_else(|| {
                    MortarError::Fatal("LevelProjector: Rotational Projector requires > Rotational Projector Periods <".into())
                })?;
                360.0 / p.max(1) as f64
            };
        }
        let setup = crate::generic::GenericSetup {
            node_coeff,
            node_scale: o.node_scale,
            max_distance,
            arc_coeff,
            xtol,
            ytol,
            xrange,
            rel_tol_y,
            repeating: o.repeating,
            anti_repeating: o.anti_repeating,
            full_circle: o.full_circle,
            gauss_points: bc.integer("mortar bc gauss points").map(|n| n.max(1) as usize),
            coord_dim: o.coord_dim,
        };
        crate::generic::weak_generic(b1, b2, &setup, &node_perm, &mut inv_perm, &mut list, log)?;
    }

    let unset = inv_perm.iter().filter(|&&p| p == elmer_core::NO_DOF).count();
    if unset > 0 {
        log.warn(CALLER, &format!("Projector % InvPerm not set in for dofs: {unset}"));
    }
    let mut rows = vec![0usize];
    let mut cols = Vec::new();
    let mut values = Vec::new();
    for row in list {
        for (c, v) in row {
            cols.push(c);
            values.push(v);
        }
        rows.push(cols.len());
    }
    log.info(CALLER, "Projector created", 10);
    Ok(Projector {
        rows,
        cols,
        values,
        inv_perm,
        projector_bc: 0,
        galerkin: true,
        slave_scale: 1.0,
        master_scale: 1.0,
        mortar_perm: false,
        rhs: None,
        diag: None,
        active: None,
    })
}

/// The `x` coordinates of an element of an interface mesh.
fn xs(m: &InterfaceMesh, e: usize) -> Vec<[f64; 3]> {
    m.element(ElemId(e as u32)).nodes.iter().map(|&n| [m.coords[n as usize][0], 0.0, 0.0]).collect()
}

/// `AddProjectorWeak1D` with `TemporalSegmentMortarAssembly`.
#[allow(clippy::too_many_arguments)]
fn weak_1d(
    b1: &InterfaceMesh,
    b2: &InterfaceMesh,
    o: &LevelOptions<'_>,
    node_coeff: f64,
    max_distance: Option<f64>,
    node_perm: &[usize],
    inv_perm: &mut [u32],
    list: &mut [BTreeMap<u32, f64>],
    log: &dyn Logger,
) -> Result<(), MortarError> {
    const CALLER: &str = "AddProjectorWeak1D";
    log.info(CALLER, "Creating weak constraints using a 1D integrator", 8);
    let seg = element_type(202).map_err(fe)?;
    let ip = gauss_points(ElementType(202), Some(seg.gauss_points2)).map_err(fe)?;
    let ip_sum: f64 = ip.s.iter().sum();
    let mut st = MetricStatus::default();
    let (mut bt, mut bs, mut bm) = (BasisEval::default(), BasisEval::default(), BasisEval::default());
    let z_of = |m: &InterfaceMesh, e: usize| -> (f64, f64) {
        let z = m.element(ElemId(e as u32)).nodes.iter().map(|&n| m.coords[n as usize][2]);
        (z.clone().fold(f64::INFINITY, f64::min), z.fold(f64::NEG_INFINITY, f64::max))
    };
    let (mut tot_hits, mut tot_ref, mut tot_sum) = (0usize, 0.0f64, 0.0f64);
    let (mut max_err, mut min_err) = (0.0f64, f64::MAX);
    let (mut max_ind, mut min_ind) = (0usize, 0usize);
    for ind in 0..b1.n_bulk() {
        let el = b1.element(ElemId(ind as u32));
        let def = element_type(el.etype.0).map_err(fe)?;
        let shape = Shape::get(el.etype.0)
            .ok_or_else(|| MortarError::Unsupported(format!("level projector on element type {}", el.etype.0)))?;
        let nodes = xs(b1, ind);
        let n = nodes.len();
        let xmin = nodes.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
        let xmax = nodes.iter().map(|c| c[0]).fold(f64::NEG_INFINITY, f64::max);
        let dx = xmax - xmin;
        let (zmin, zmax) = z_of(b1, ind);
        element_info(def, &nodes, [0.0; 3], o.coord_dim, &mut st, &mut bs).map_err(fe)?;
        let ref_area = bs.det_j * ip_sum;
        let mut sum_area = 0.0f64;
        // List_AddMatrixIndex: the rows of the slave nodes exist.
        for &k in el.nodes {
            let j = b1.inv_perm[k as usize];
            let nrow = node_perm[j as usize];
            if nrow != 0 {
                list[nrow - 1].entry(j).or_insert(0.0);
            }
        }
        let mut hits = 0usize;
        for indm in 0..b2.n_bulk() {
            let elm = b2.element(ElemId(indm as u32));
            let nodes_m = xs(b2, indm);
            let xminm = nodes_m.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
            let xmaxm = nodes_m.iter().map(|c| c[0]).fold(f64::NEG_INFINITY, f64::max);
            if xminm >= xmax || xmaxm <= xmin {
                continue;
            }
            if let Some(md) = max_distance {
                let (zminm, zmaxm) = z_of(b2, indm);
                if zmaxm < zmin - md || zminm > zmax + md {
                    continue;
                }
            }
            let cut = [[xmin.max(xminm), 0.0, 0.0], [xmax.min(xmaxm), 0.0, 0.0]];
            let dxcut = cut[1][0] - cut[0][0];
            if dxcut < 1.0e-12 || dxcut < 1.0e-8 * dx {
                continue;
            }
            hits += 1;
            let def_m = element_type(elm.etype.0).map_err(fe)?;
            let shape_m = Shape::get(elm.etype.0)
                .ok_or_else(|| MortarError::Unsupported(format!("level projector on element type {}", elm.etype.0)))?;
            // TemporalSegmentMortarAssembly, sgn0 = 1.
            for k in 0..ip.s.len() {
                element_info(seg, &cut, [ip.u[k], ip.v[k], ip.w[k]], o.coord_dim, &mut st, &mut bt).map_err(fe)?;
                let xt = bt.basis[0] * cut[0][0] + bt.basis[1] * cut[1][0];
                let wtemp = bt.det_j * ip.s[k];
                let uvw = global_to_local(shape, &nodes, [xt, 0.0, 0.0], o.coord_dim);
                element_info(def, &nodes, uvw, o.coord_dim, &mut st, &mut bs).map_err(fe)?;
                sum_area += wtemp;
                let uvw_m = global_to_local(shape_m, &nodes_m, [xt, 0.0, 0.0], o.coord_dim);
                element_info(def_m, &nodes_m, uvw_m, o.coord_dim, &mut st, &mut bm).map_err(fe)?;
                for j in 0..n {
                    let jj = b1.inv_perm[el.nodes[j] as usize];
                    let nrow = node_perm[jj as usize];
                    if nrow == 0 {
                        continue;
                    }
                    inv_perm[nrow - 1] = jj;
                    let val = node_coeff * bs.basis[j] * wtemp;
                    let row = &mut list[nrow - 1];
                    for i in 0..n {
                        *row.entry(b1.inv_perm[el.nodes[i] as usize]).or_insert(0.0) += bs.basis[i] * val;
                    }
                    for i in 0..elm.nodes.len() {
                        *row.entry(b2.inv_perm[elm.nodes[i] as usize]).or_insert(0.0) +=
                            -o.node_scale * bm.basis[i] * val;
                    }
                }
            }
        }
        tot_hits += hits;
        tot_sum += sum_area;
        tot_ref += ref_area;
        let err = sum_area / ref_area;
        if err > max_err {
            max_err = err;
            max_ind = ind + 1;
        }
        if err < min_err {
            min_err = err;
            min_ind = ind + 1;
        }
    }
    log.info(CALLER, &format!("Number of integration pairs: {tot_hits}"), 10);
    log.info(CALLER, &format!("Total reference length:{}", elmer_solver::fortran_fmt::es(tot_ref, 12, 5)), 8);
    log.info(CALLER, &format!("Total integrated length:{}", elmer_solver::fortran_fmt::es(tot_sum, 12, 5)), 8);
    log.info(
        CALLER,
        &format!("Average ratio in length integration:{}", elmer_solver::fortran_fmt::es(tot_sum / tot_ref, 12, 3)),
        8,
    );
    log.info(
        CALLER,
        &format!(
            "Maximum relative discrepancy in length (element: {max_ind}):{}",
            elmer_solver::fortran_fmt::es(max_err - 1.0, 12, 4)
        ),
        8,
    );
    log.info(
        CALLER,
        &format!(
            "Minimum relative discrepancy in length (element: {min_ind}):{}",
            elmer_solver::fortran_fmt::es(min_err - 1.0, 12, 4)
        ),
        8,
    );
    Ok(())
}
