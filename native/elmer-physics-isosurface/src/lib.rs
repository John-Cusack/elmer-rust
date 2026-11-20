//! `IsosurfaceSolver` (fem/src/modules/Isosurface.F90), serial 2D meshes:
//! the isolines of a nodal field at the `Isosurface Values`, built as
//! Elmer builds them (the crossed bulk elements split into a temporal
//! triangle mesh, its edges, a node on every edge the level crosses, a
//! segment per crossed triangle). The solver's global variable and norm are
//! the sum of the six sizes Elmer reports (mesh nodes and elements, temporal
//! elements and edges, isoline nodes and elements).
//!
//! `Save Gmsh Geo File` writes the isolines as Gmsh splines
//! (`SaveGmshGeo2D`). The other uses of the isoline mesh are not produced
//! here and are refused: `Save STL File`, `Isosurface Interpolant`
//! fields on an isomesh that a `ResultOutputSolver` writes, 3D meshes
//! (isosurfaces), `Create Isosurface` masks, `Moving Mesh`, parallel runs.

use std::collections::HashMap;

use elmer_core::{ElemId, ParamMap, Params, Value, NO_DOF};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

pub const PROCEDURE: (&str, &str) = ("Isosurface", "IsosurfaceSolver");
const CALLER: &str = "IsosurfaceSolver";

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(Isosurface::default()))
}

/// Whether `params` is an IsosurfaceSolver section.
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(lib, entry)| lib.eq_ignore_ascii_case(PROCEDURE.0) && entry.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `IsosurfaceSolver_init`: `Variable = -nooutput -global <Equation>_var`
/// unless the section names one.
pub fn init_section(params: &mut ParamMap) {
    if params.get("variable").is_none() {
        let name = params.string("equation").unwrap_or("").to_string();
        params.insert("Variable", Value::Str(format!("-nooutput -global {name}_var")));
    }
}

#[derive(Default)]
struct Isosurface {
    /// `SurfaceExist` and `SaveParTmp` for `Isosurface Fixed`.
    saved: Option<[i64; 6]>,
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("IsosurfaceSolver: {what}"))
}

impl Solver for Isosurface {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in ["save stl file", "moving mesh"] {
            if p.logical(k) == Some(true) {
                return Err(unsupported(&format!("'{k}'")));
            }
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("parallel runs"));
        }
        if ctx.sim.model.boundary_conditions.iter().chain(&ctx.sim.model.body_forces).any(|s| s.get("create isosurface").is_some()) {
            return Err(unsupported("'Create Isosurface' masks"));
        }
        let output_active = p.logical("isomesh output active").unwrap_or(true);
        let interpolants = p.0.keys().any(|k| k.starts_with("isosurface interpolant"));
        let result_output = ctx.sim.model.solvers.iter().any(|s| {
            elmer_solver::solver::procedure_of(s).is_some_and(|(_, e)| e.eq_ignore_ascii_case("resultoutputsolver"))
        });
        // An `Exec Solver = After All` isosurface comes after the last
        // output of the time loop: no ResultOutputSolver sees its mesh.
        let after_all = p.string("exec solver").is_some_and(|e| e.trim().eq_ignore_ascii_case("after all"));
        if output_active && result_output && !after_all {
            return Err(unsupported("an isoline mesh for ResultOutputSolver"));
        }
        let _ = interpolants; // interpolated onto the isomesh only, which nothing writes here
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "-------------------------------------", 4);
        ctx.info(CALLER, "Determining the isosurface", 4);
        ctx.info(CALLER, "-------------------------------------", 4);
        let p = ctx.params().clone();
        let fixed = p.logical("isosurface fixed") == Some(true);
        let sizes = match (fixed, self.saved) {
            (true, Some(s)) => {
                ctx.info(CALLER, "Remapping the fields", 4);
                Some(s)
            }
            _ => match build(ctx, &p)? {
                Some((s, iso)) => {
                    self.saved = Some(s);
                    if p.logical("save gmsh geo file") == Some(true) {
                        save_gmsh_geo(ctx, &p, &iso)?;
                    }
                    Some(s)
                }
                None => None,
            },
        };
        let sizes = sizes.unwrap_or([0; 6]);
        if self.saved.is_some() {
            let msg = [
                ("Initial mesh has ", " nodes"),
                ("Initial mesh has ", " elements"),
                ("Temporal mesh has ", " elements"),
                ("Temporal mesh has ", " edges"),
                ("Isosurface mesh has ", " nodes"),
                ("Isosurface mesh has ", " elements"),
            ];
            ctx.info("IsoSurfaceSolver", "Information on isosurface mesh sizes", 4);
            for (k, (a, b)) in msg.iter().enumerate() {
                ctx.info("IsoSurfaceSolver", &format!("{a}{}{b}", sizes[k]), 4);
            }
        }
        // `Solver % Variable % Values = SUM(ParSizes)`, its norm too.
        let total: i64 = sizes.iter().sum();
        let id = ctx.solver_id();
        let ix = ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot");
        if let Some(vid) = ctx.sim.slots[ix].variable {
            let v = ctx.sim.vars.by_id_mut(vid);
            if v.values.len() == 1 {
                v.values[0] = total as f64;
                v.conv.norm = total as f64;
            }
        }
        Ok(())
    }
}

/// The level function at a node slot: the variable's value, or the norm of
/// its first `MeshDim` components.
struct Level<'a> {
    values: &'a [f64],
    perm: Option<&'a [u32]>,
    dofs: usize,
    comps: usize,
}

impl Level<'_> {
    /// `LevelPerm(node)` (1-based in Elmer, 0 for none): the slot, if any.
    fn slot(&self, node: u32) -> Option<usize> {
        match self.perm {
            Some(p) => {
                let s = p[node as usize];
                (s != NO_DOF).then_some(s as usize)
            }
            None => Some(node as usize),
        }
    }
    fn value(&self, slot: usize) -> f64 {
        if self.dofs == 1 {
            self.values[slot]
        } else {
            let mut s = 0.0;
            for l in 0..self.comps {
                let v = self.values[self.dofs * slot + l];
                s += v * v;
            }
            s.sqrt()
        }
    }
}

/// The isoline construction; `None` where Elmer jumps to `100` before
/// creating the mesh (missing variable or levels, no 2D elements).
fn build(ctx: &mut SolverCtx<'_>, p: &ParamMap) -> Result<Option<([i64; 6], IsoLines)>, SolverError> {
    let mesh = ctx.mesh;
    let name = p.string("isosurface variable").unwrap_or("").to_string();
    let coords: Vec<f64>;
    let field = ctx.variable(&name);
    let level = match field {
        Some(f) => Level { values: f.values, perm: f.perm, dofs: f.dofs, comps: mesh.dim().min(f.dofs) },
        None => {
            // `Coordinate n`: the mesh coordinates, no permutation.
            let lower = name.trim().to_ascii_lowercase();
            let Some(c) = lower.strip_prefix("coordinate ").and_then(|s| s.trim().parse::<usize>().ok()).filter(|c| (1..=3).contains(c)) else {
                ctx.log.warn("Isosurface", &format!("Missing isosurface variable: {name}"));
                return Ok(None);
            };
            coords = (0..mesh.n_nodes()).map(|i| mesh.coord(i as u32)[c - 1]).collect();
            Level { values: &coords, perm: None, dofs: 1, comps: 1 }
        }
    };
    let levels: Vec<f64> = match p.get("isosurface values") {
        Some(Value::Real { data, .. }) => data.clone(),
        Some(Value::Integer(v)) => v.iter().map(|&i| i as f64).collect(),
        _ => match ctx.creal(p, "isosurface value")? {
            Some(v) => vec![v],
            None => {
                ctx.log.warn(CALLER, "Could not determine Isosurface value");
                return Ok(None);
            }
        },
    };

    let n_bulk = mesh.n_bulk();
    let mut dim = 0;
    for e in 0..n_bulk {
        let fam = mesh.element(ElemId(e as u32)).etype.family();
        if fam <= 2 {
            continue;
        } else if fam <= 4 {
            dim = dim.max(2);
        } else {
            dim = dim.max(3);
        }
    }
    ctx.info(CALLER, &format!("Leading element dimension is {dim}"), 4);
    if dim == 3 {
        return Err(unsupported("3D meshes (isosurfaces)"));
    }
    if dim != 2 {
        ctx.log.warn("IsoSurfaceSolver", "Isosurface mesh can be created only for 2D or 3D");
        return Ok(None);
    }

    // The temporal triangle mesh of the crossed elements.
    let mut tris: Vec<[u32; 3]> = Vec::new();
    for e in 0..n_bulk {
        let el = mesh.element(ElemId(e as u32));
        let n = el.etype.n_nodes();
        let nodes = &el.nodes[..n];
        let mut fun = Vec::with_capacity(n);
        let mut skip = false;
        for &nd in nodes {
            match level.slot(nd) {
                Some(s) => fun.push(level.value(s)),
                None => {
                    skip = true;
                    break;
                }
            }
        }
        if skip {
            continue;
        }
        let crossed = levels.iter().any(|&lv| {
            !(fun.iter().all(|f| f - lv < 0.0) || fun.iter().all(|f| f - lv > 0.0))
        });
        if !crossed {
            continue;
        }
        let map: &[[usize; 3]] = match el.etype.family() {
            3 => &[[0, 1, 2]],
            4 => &[[0, 1, 2], [0, 2, 3]],
            f => return Err(unsupported(&format!("bulk elements of family {f}"))),
        };
        for m in map {
            tris.push([nodes[m[0]], nodes[m[1]], nodes[m[2]]]);
        }
    }
    let n_new = tris.len();
    if n_new == 0 {
        ctx.info(CALLER, "No potential elements found in this partition", 7);
    } else {
        ctx.info(CALLER, &format!("Found {n_new} potential elements"), 4);
    }

    // FindMeshEdges2D: the triangle edges (1,2), (2,3), (3,1), numbered by
    // first appearance.
    let mut edge_of: HashMap<(u32, u32), usize> = HashMap::new();
    let mut edges: Vec<(u32, u32)> = Vec::new();
    let mut tri_edges: Vec<[usize; 3]> = Vec::with_capacity(n_new);
    for t in &tris {
        let mut ix = [0usize; 3];
        for (k, (a, b)) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])].into_iter().enumerate() {
            let key = (a.min(b), a.max(b));
            ix[k] = *edge_of.entry(key).or_insert_with(|| {
                edges.push((a, b));
                edges.len() - 1
            });
        }
        tri_edges.push(ix);
    }
    let n_edges = edges.len();

    // CreateNodes: a node on each edge a level crosses; `eperm` marks them.
    let mut eperm = vec![0usize; n_edges];
    let mut n_iso_nodes = 0usize;
    let mut iso = IsoLines::default();
    for &lv in &levels {
        for (i, &(n1, n2)) in edges.iter().enumerate() {
            let (Some(m1), Some(m2)) = (level.slot(n1), level.slot(n2)) else { continue };
            let t1 = level.value(m1) - lv;
            let t2 = level.value(m2) - lv;
            if (t1 - t2).abs() < f64::MIN_POSITIVE {
                continue;
            }
            if t1 * t2 > 0.0 {
                continue;
            }
            n_iso_nodes += 1;
            eperm[i] = n_iso_nodes;
            // `t = ABS(t1/(t2-t1))`, the node at `(1-t) x1 + t x2`.
            let t = (t1 / (t2 - t1)).abs();
            let (x1, x2) = (mesh.coord(n1), mesh.coord(n2));
            iso.nodes.push(std::array::from_fn(|a| (1.0 - t) * x1[a] + t * x2[a]));
        }
    }

    // CreateSurfaces: CreateLineFromTriangle per level.
    let mut n_surf = 0usize;
    if n_iso_nodes > 0 {
        for &lv in &levels {
            for (t, ix) in tris.iter().zip(&tri_edges) {
                let mut f = [0.0; 3];
                for k in 0..3 {
                    f[k] = level.slot(t[k]).map_or(0.0, |s| level.value(s)) - lv;
                }
                let idx = [eperm[ix[0]], eperm[ix[1]], eperm[ix[2]]];
                let (s1, s2, s3) = (f[0] > 0.0, f[1] > 0.0, f[2] > 0.0);
                let s = s1 as u8 + s2 as u8 + s3 as u8;
                if s == 0 || s == 3 {
                    continue;
                }
                let line = if s == 1 {
                    if s1 {
                        [idx[0], idx[2]]
                    } else if s2 {
                        [idx[0], idx[1]]
                    } else {
                        [idx[1], idx[2]]
                    }
                } else if !s1 {
                    [idx[0], idx[2]]
                } else if !s2 {
                    [idx[0], idx[1]]
                } else {
                    [idx[1], idx[2]]
                };
                if line.iter().all(|&j| j > 0) {
                    n_surf += 1;
                    iso.lines.push(line);
                }
            }
        }
    }
    Ok(Some((
        [
            mesh.n_nodes() as i64,
            n_bulk as i64,
            n_new as i64,
            n_edges as i64,
            n_iso_nodes as i64,
            n_surf as i64,
        ],
        iso,
    )))
}

/// The isoline mesh: its nodes and its line elements (1-based nodes).
#[derive(Default)]
struct IsoLines {
    nodes: Vec<[f64; 3]>,
    lines: Vec<[usize; 2]>,
}

/// gfortran's `ESw.d`.
fn es(x: f64, w: usize, d: usize) -> String {
    elmer_solver::fortran_fmt::es(x, w, d)
}

/// `SaveGmshGeo2D` (Isosurface.F90:1853): the closed isoline loops as Gmsh
/// splines, the longest loop first, each with its `Line Loop` and `Plane
/// Surface`. Nothing is written for an open curve.
fn save_gmsh_geo(ctx: &mut SolverCtx<'_>, p: &ParamMap, iso: &IsoLines) -> Result<(), SolverError> {
    let n = iso.nodes.len();
    if n == 0 {
        return Ok(());
    }
    let filename = p.string("geo filename").unwrap_or("mesh.geo").to_string();
    // Neighbours(:,1:2), 1-based, 0 for none.
    let mut nb = vec![[0usize; 2]; n + 1];
    for &[i1, i2] in &iso.lines {
        if i1 == 0 || i2 == 0 {
            ctx.log.warn("SaveGmshGeo2D", &format!("Invalid indexes: {i1} and {i2}"));
        }
        if nb[i1][0] == 0 { nb[i1][0] = i2 } else { nb[i1][1] = i2 }
        if nb[i2][0] == 0 { nb[i2][0] = i1 } else { nb[i2][1] = i1 }
    }
    if nb[1..].iter().any(|b| b[1] == 0) {
        ctx.log.warn("SaveGmshGeo2D", "This does not seem to be a closed loop!");
        return Ok(());
    }
    let span = |a: usize| {
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for c in &iso.nodes {
            lo = lo.min(c[a]);
            hi = hi.max(c[a]);
        }
        hi - lo
    };
    let (dx, dy, dz) = (span(0), span(1), span(2));
    let diam = dx.max(dy.max(dz));
    let mesh_param = ctx.creal(p, "mesh parameter")?.unwrap_or(diam / 50.0);
    let node_eps = ctx.creal(p, "mesh node epsilon")?.unwrap_or(1.0e-3 * diam);
    const MAX_LOOPS: usize = 20;
    let mut used = vec![false; n + 1];
    let mut order = vec![0usize; n + 1];
    let (mut offset, mut size) = (vec![0i64; MAX_LOOPS + 1], vec![0i64; MAX_LOOPS + 1]);
    let (mut j, mut used_nodes, mut no_loop, mut new_loop) = (1usize, 0usize, 1usize, true);
    let mut prev = [0.0f64; 3];
    for _ in 0..n {
        used[j] = true;
        let c = iso.nodes[j - 1];
        let mut save = true;
        if !new_loop {
            let mut dist = 0.0;
            for a in 0..3 {
                dist += (c[a] - prev[a]) * (c[a] - prev[a]);
            }
            if dist < node_eps * node_eps {
                save = false;
            }
        }
        if save {
            size[no_loop] += 1;
            used_nodes += 1;
            order[used_nodes] = j;
            prev = c;
        }
        new_loop = false;
        if !used[nb[j][0]] {
            j = nb[j][0];
        } else if !used[nb[j][1]] {
            j = nb[j][1];
        } else {
            match (1..=n).find(|&k| !used[k]) {
                Some(k) => {
                    j = k;
                    if no_loop == MAX_LOOPS {
                        ctx.log.warn(CALLER, &format!("The static max number of loops ({MAX_LOOPS}) is too small, ignoring the rest!"));
                        break;
                    }
                    no_loop += 1;
                    new_loop = true;
                    offset[no_loop] = used_nodes as i64;
                }
                None => break,
            }
        }
    }
    ctx.info(CALLER, if used_nodes == n { "All nodes used".to_string() } else { format!("{used_nodes} nodes used") }.as_str(), 4);
    ctx.info("IsoSurfaceSolver", &format!("Number of loops found: {no_loop}"), 4);
    let save_loops = p.integer("save number of loops").map_or(no_loop, |k| (k.max(0) as usize).min(no_loop));
    let largest = |size: &[i64]| {
        let (mut j, mut m) = (0usize, 0i64);
        for (k, &s) in size.iter().enumerate().take(no_loop + 1).skip(1) {
            if s > m {
                j = k;
                m = s;
            }
        }
        j
    };
    let mut out = String::new();
    let mut node_index = 0usize;
    for _ in 0..save_loops {
        let j = largest(&size);
        for i in 1..=size[j] as usize {
            node_index += 1;
            let k = order[offset[j] as usize + i];
            let c = iso.nodes[k - 1];
            out.push_str(&format!(
                "Point({node_index}) = {{{},{},{},{}}};\n",
                es(c[0], 12, 4),
                es(c[1], 12, 4),
                es(c[2], 12, 4),
                es(mesh_param, 12, 4)
            ));
        }
        size[j] = -size[j];
    }
    for s in size.iter_mut() {
        *s = s.abs();
    }
    node_index = 0;
    for l in 1..=save_loops {
        let j = largest(&size);
        out.push_str(&format!("Spline({l}) = {{"));
        for i in 1..=size[j] as usize {
            out.push_str(&format!("{}, ", node_index + i));
        }
        out.push_str(&format!("{}}};\n", node_index + 1));
        node_index += size[j] as usize;
        size[j] = -size[j];
    }
    for j in 1..=save_loops {
        let k = save_loops + 2 * (j - 1);
        out.push_str(&format!("Line Loop({}) = {{{j}}};\n", k + 1));
        out.push_str(&format!("Plane Surface({}) = {{{}}};\n", k + 2, k + 1));
    }
    if let Err(e) = std::fs::write(&filename, out) {
        ctx.log.warn("SaveGmshGeo2D", &format!("Could not open file: {filename}: {e}"));
        return Ok(());
    }
    ctx.info(CALLER, &format!("Wrote geo file: {filename}"), 6);
    Ok(())
}
