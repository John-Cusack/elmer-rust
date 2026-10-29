//! `DetectExtrudedStructure` (MeshTransform.F90:52): finds, for every node of
//! a mesh that could have been made by extrusion along a coordinate axis, the
//! node above and below it and the top, bottom and middle node of its column.

use elmer_core::{ElemId, MeshView, Params};
use elmer_solver::{Model, SolverError};

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

/// `unit` is the axis (negated by `Project To Bottom`), `eps` the `Dot
/// Product Tolerance`.
pub fn detect_extruded_structure(
    mesh: &dyn MeshView,
    model: &Model,
    unit: [f64; 3],
    eps: f64,
    want: Want,
) -> Result<Structure, SolverError> {
    let nnodes = mesh.n_nodes();
    let identity: Vec<usize> = (0..nnodes).collect();
    let (mut top, mut up) = (identity.clone(), identity.clone());
    let (mut bot, mut down) = (identity.clone(), identity);

    // Up and down pointers by the dot product criterion.
    for t in 0..mesh.n_bulk() as u32 {
        let nodes = mesh.element(ElemId(t)).nodes;
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
    while rounds <= nnodes {
        let mut hits = 0;
        for i in 0..nnodes {
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

    let mid = if want.mid { find_mid(mesh, model, &up, &down) } else { None };
    Ok(Structure { top, bot, up, down, mid, layers })
}

/// The `Mid Surface` pointers (MeshTransform.F90:481-545).
fn find_mid(mesh: &dyn MeshView, model: &Model, up: &[usize], down: &[usize]) -> Option<Vec<usize>> {
    const NONE: usize = usize::MAX;
    let nnodes = mesh.n_nodes();
    let mut mid = vec![NONE; nnodes];
    let mut exists = false;
    let nb = mesh.n_bulk() as u32;
    for t in nb..nb + mesh.n_boundary() as u32 {
        let el = mesh.element(ElemId(t));
        if let Some(bc) = model.bc_of_tag(el.tag) {
            if model.boundary_conditions[bc - 1].get("mid surface").is_some() {
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
