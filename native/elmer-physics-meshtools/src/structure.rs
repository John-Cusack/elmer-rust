//! `DetectExtrudedStructure` (MeshTransform.F90:52): finds, for every node of
//! a mesh that could have been made by extrusion along a coordinate axis, the
//! node above and below it and the top, bottom and middle node of its column.

use elmer_core::{ElemId, FieldKind, MeshView, ParamMap, Params, NO_DOF};
use elmer_solver::{Model, SolverCtx, SolverError, Variable};

const CALLER: &str = "DetectExtrudedStructure";

/// The pointers are node indexes (0-based). `up`/`down` of the end of a column
/// point at the node itself.
#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub top: Vec<usize>,
    pub bot: Vec<usize>,
    pub up: Vec<usize>,
    pub down: Vec<usize>,
    /// The node of the `Mid Surface` boundary condition in each column.
    pub mid: Option<Vec<usize>>,
    /// Layers counted along the first column, plus the fatal-or-warn rounds
    /// of the Fortran (`NumberOfLayers`, before the caller adds one).
    pub layers: usize,
}

/// What to find, as the optional arguments of the Fortran routine.
#[derive(Copy, Clone, Debug, Default)]
pub struct Want {
    /// `NumberOfLayers` or `NodeLayer` present: both directions.
    pub layers: bool,
    /// `MidNodePointer` present.
    pub mid: bool,
}

/// The nodes an extrusion is restricted to (`MaskPerm`): node -> 0-based
/// slot, [`NO_DOF`] outside; `size` slots (`nsize`).
#[derive(Clone, Debug)]
pub struct Mask {
    pub perm: Vec<u32>,
    pub size: usize,
}

/// `MakePermUsingMask(Model, Solver, Mesh, name, .FALSE., Perm, n)`
/// (MeshBasics.F90:6040) without bandwidth optimisation: the nodes of the
/// bulk elements whose body force has `name` (when some body force has it)
/// and of the boundary elements whose BC has it (when some BC has it),
/// numbered in first-visit order. The keyword only needs to be present.
pub fn make_perm_using_mask(mesh: &dyn MeshView, model: &Model, name: &str) -> Mask {
    let has = |p: &ParamMap| p.get(name).is_some();
    let nb = mesh.n_bulk();
    let mut start = usize::MAX;
    let mut fin = 0;
    if model.body_forces.iter().any(has) {
        start = 0;
        fin = nb;
    }
    if model.boundary_conditions.iter().any(has) {
        start = start.min(nb);
        fin = nb + mesh.n_boundary();
    }
    let mut perm = vec![NO_DOF; mesh.n_nodes()];
    let mut k = 0u32;
    if fin > start {
        for t in start..fin {
            let el = mesh.element(ElemId(t as u32));
            let hit = if t < nb {
                model
                    .body_of_tag(el.tag)
                    .and_then(|b| model.bodies[b - 1].integer("body force"))
                    .filter(|&f| f >= 1)
                    .and_then(|f| model.body_forces.get(f as usize - 1))
                    .is_some_and(has)
            } else {
                model.bc_of_tag(el.tag).is_some_and(|b| has(&model.boundary_conditions[b - 1]))
            };
            if !hit {
                continue;
            }
            for &n in el.nodes {
                if perm[n as usize] == NO_DOF {
                    perm[n as usize] = k;
                    k += 1;
                }
            }
        }
    }
    Mask { perm, size: k as usize }
}

/// The mask of `DetectExtrudedStructure`: `Mapping Mask Variable` (the
/// variable's Perm, `nsize = MAXVAL(Perm)`), else `Mapping Mask Name`
/// ([`make_perm_using_mask`]); `None` without either.
pub fn mapping_mask(ctx: &SolverCtx<'_>, params: &ParamMap) -> Result<Option<Mask>, SolverError> {
    if let Some(name) = params.string("mapping mask variable") {
        if let Some(v) = ctx.sim.vars.get(name) {
            let Some(p) = v.perm.as_ref() else {
                return Err(SolverError::Unsupported(format!(
                    "{CALLER}: Mapping Mask Variable {name} without a permutation"
                )));
            };
            if v.kind != FieldKind::Nodal || p.len() != ctx.mesh.n_nodes() {
                return Err(SolverError::Unsupported(format!("{CALLER}: Mapping Mask Variable {name} that is not nodal")));
            }
            let size = p.iter().filter(|&&s| s != NO_DOF).map(|&s| s as usize + 1).max().unwrap_or(0);
            return Ok(Some(Mask { perm: p.to_vec(), size }));
        }
    }
    Ok(params
        .string("mapping mask name")
        .map(|name| make_perm_using_mask(ctx.mesh, &ctx.sim.model, &elmer_core::normalize_key(name))))
}

/// The `Extruded Coordinate` variable of a masked extrusion: created on the
/// mask's Perm (or reused), and set to the nodes' coordinate along `dir`
/// at every detection (MeshTransform.F90:157-196).
pub fn extruded_coordinate(ctx: &mut SolverCtx<'_>, mask: &Mask, dir: usize) -> Result<(), SolverError> {
    const NAME: &str = "extruded coordinate";
    if ctx.sim.vars.get(NAME).is_none() {
        let v = Variable::new("Extruded Coordinate", FieldKind::Nodal, 1, Some(mask.perm.clone().into()), vec![0.0; mask.size]);
        ctx.sim.vars.add(v)?;
    }
    let coords: Vec<f64> = (0..ctx.mesh.n_nodes() as u32).map(|n| ctx.mesh.coord(n)[dir]).collect();
    let v = ctx.sim.vars.get_mut(NAME).expect("added above");
    let perm = v.perm.clone().expect("has a perm");
    for (n, &s) in perm.iter().enumerate() {
        if s != NO_DOF {
            v.values[s as usize] = coords[n];
        }
    }
    v.values_changed = true;
    Ok(())
}

/// `unit` is the axis (negated by `Project To Bottom`), `eps` the `Dot
/// Product Tolerance`. With a `mask` only its nodes take part (an element
/// counts when all its nodes are masked, boundary elements included); the
/// pointers stay indexed by node.
pub fn detect_extruded_structure(
    mesh: &dyn MeshView,
    model: &Model,
    unit: [f64; 3],
    eps: f64,
    want: Want,
    mask: Option<&Mask>,
) -> Result<Structure, SolverError> {
    let nnodes = mesh.n_nodes();
    let identity: Vec<usize> = (0..nnodes).collect();
    let (mut top, mut up) = (identity.clone(), identity.clone());
    let (mut bot, mut down) = (identity.clone(), identity);
    let masked = |n: usize| mask.is_none_or(|m| m.perm[n] != NO_DOF);
    let nsize = mask.map_or(nnodes, |m| m.size);
    let elemmax = if mask.is_some() { mesh.n_bulk() + mesh.n_boundary() } else { mesh.n_bulk() };

    // Up and down pointers by the dot product criterion.
    for t in 0..elemmax as u32 {
        let nodes = mesh.element(ElemId(t)).nodes;
        if mask.is_some() && !nodes.iter().all(|&n| masked(n as usize)) {
            continue;
        }
        for (i, &ni) in nodes.iter().enumerate() {
            let v = mesh.coord(ni);
            let (mut max_dot, mut min_dot) = (-1.0f64, 1.0f64);
            let (mut jmax, mut jmin) = (0usize, 0usize);
            for &nj in &nodes[i + 1..] {
                let v2 = mesh.coord(nj);
                let e = [v2[0] - v[0], v2[1] - v[1], v2[2] - v[2]];
                let length = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
                let dot = (e[0] * unit[0] + e[1] * unit[1] + e[2] * unit[2]) / length;
                if dot > max_dot {
                    max_dot = dot;
                    jmax = nj as usize;
                }
                if dot < min_dot {
                    min_dot = dot;
                    jmin = nj as usize;
                }
            }
            if max_dot > 1.0 - eps {
                up[ni as usize] = jmax;
                down[jmax] = ni as usize;
            }
            if min_dot < eps - 1.0 {
                down[ni as usize] = jmin;
                up[jmin] = ni as usize;
            }
        }
    }

    // Top and bottom recursively through up and down.
    let mut rounds = 1;
    while rounds <= nsize {
        let mut hits = 0;
        for i in 0..nnodes {
            if !masked(i) {
                continue;
            }
            let j = up[i];
            if top[i] != top[j] {
                hits += 1;
                top[i] = top[j];
            }
            let j = down[i];
            if bot[i] != bot[j] {
                hits += 1;
                bot[i] = bot[j];
            }
        }
        if hits == 0 {
            break;
        }
        rounds += 1;
    }
    // The last round is always a check.
    let rounds = rounds - 1;
    if rounds == 0 {
        return Err(SolverError::fatal(
            CALLER,
            "Zero rounds implies unsuccessful operation (try to increase value for > Dot Product Tolerance <)",
        ));
    }

    let mut layers = 0;
    if want.layers {
        let mut j = bot[0];
        loop {
            let k = up[j];
            if k == j {
                break;
            }
            layers += 1;
            j = k;
        }
        if layers < rounds {
            layers = rounds;
        }
    }

    let mid = if want.mid {
        let mid = find_mid(mesh, model, &up, &down);
        if mid.is_some() && mask.is_some() {
            // MidPointer is indexed by node with a mask-sized table there.
            return Err(SolverError::Unsupported(format!("{CALLER}: Mid Surface with a mapping mask")));
        }
        mid
    } else {
        None
    };
    Ok(Structure {
        top,
        bot,
        up,
        down,
        mid,
        layers,
    })
}

/// The `Mid Surface` pointers (MeshTransform.F90:481-545).
fn find_mid(
    mesh: &dyn MeshView,
    model: &Model,
    up: &[usize],
    down: &[usize],
) -> Option<Vec<usize>> {
    const NONE: usize = usize::MAX;
    let nnodes = mesh.n_nodes();
    let mut mid = vec![NONE; nnodes];
    let mut exists = false;
    let nb = mesh.n_bulk() as u32;
    for t in nb..nb + mesh.n_boundary() as u32 {
        let el = mesh.element(ElemId(t));
        if let Some(bc) = model.bc_of_tag(el.tag) {
            if model.boundary_conditions[bc - 1]
                .get("mid surface")
                .is_some()
            {
                for &n in el.nodes {
                    mid[n as usize] = n as usize;
                }
                exists = true;
            }
        }
    }
    if !exists {
        return None;
    }
    for _ in 0..nnodes {
        let mut hits = 0;
        for i in 0..nnodes {
            if mid[i] == NONE {
                continue;
            }
            let j = up[i];
            if mid[j] == NONE {
                hits += 1;
                mid[j] = mid[i];
            }
            let j = down[i];
            if mid[j] == NONE {
                hits += 1;
                mid[j] = mid[i];
            }
        }
        if hits == 0 {
            break;
        }
    }
    Some(mid)
}
