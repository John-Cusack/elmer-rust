//! `GmshOutputReader` (fem/src/modules/GmshOutputReader.F90): reads a
//! Gmsh 2 ASCII file with `$NodeData` (e.g. another case's results),
//! optionally shifts its mesh along `Align Coordinate` so that it ends where
//! this mesh begins, and interpolates its fields onto the same-named
//! variables of this mesh (`InterpolateMeshToMeshQ`), on the nodes of the
//! BCs with `Mask Name` true when one is given.
//!
//! Refused: `Filename Partitions` > 1, masks on bodies or body forces, and
//! multi-component fields.

use std::path::PathBuf;

use elmer_adaptive::mesh::{Elem, OwnedMesh};
use elmer_core::{ElemId, ElementType, FieldKind, ParamMap, Params, Value, NO_DOF};
use elmer_interp::{build_projector, FieldData, InterpOptions, Location, Locator};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

pub(crate) const PROCEDURE: (&str, &str) = ("GmshOutputReader", "GmshOutputReader");
const CALLER: &str = "GmshOutputReader";

pub(crate) fn is_section(p: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(p)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

/// `GmshOutputReader_init`: no matrix.
pub(crate) fn init_section(p: &mut ParamMap) {
    if p.get("no matrix").is_none() {
        p.insert("No Matrix", Value::Logical(true));
    }
}

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(GmshOutputReader { case_dir: PathBuf::new() })));
}

pub(crate) struct GmshOutputReader {
    pub case_dir: PathBuf,
}

fn fatal(m: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, m)
}

/// `GmshToElmerType`.
const GMSH_TO_ELMER: [u16; 21] = [202, 303, 404, 504, 808, 706, 605, 203, 306, 409, 510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310];

struct GmshData {
    coords: Vec<[f64; 3]>,
    elems: Vec<(u16, Vec<u32>)>,
    n_bulk: usize,
    dim: usize,
    fields: Vec<(String, Vec<f64>)>,
}

fn nums<T: std::str::FromStr>(line: &str) -> Vec<T> {
    line.split_whitespace().filter_map(|t| t.parse().ok()).collect()
}

fn next<'a>(l: &mut std::str::Lines<'a>) -> Result<&'a str, SolverError> {
    l.next().ok_or_else(|| fatal("unexpected end of the Gmsh file"))
}

/// `ReadSingleGmshFile`.
fn read_gmsh(text: &str) -> Result<GmshData, SolverError> {
    let mut lines = text.lines();
    let mut d = GmshData { coords: Vec::new(), elems: Vec::new(), n_bulk: 0, dim: 0, fields: Vec::new() };
    let mut max_dim = 0usize;
    let mut mesh_dim = 0usize;
    while let Some(line) = lines.next() {
        let s = line.trim_start();
        if s.starts_with("$Nodes") {
            let n: usize = nums(next(&mut lines)?).first().copied().ok_or_else(|| fatal("node count"))?;
            for i in 0..n {
                let v: Vec<f64> = nums(next(&mut lines)?);
                if v.len() < 4 || v[0] as usize != i + 1 {
                    return Err(SolverError::Unsupported(format!("{CALLER}: node permutations")));
                }
                d.coords.push([v[1], v[2], v[3]]);
            }
        } else if s.starts_with("$Elements") {
            let n: usize = nums(next(&mut lines)?).first().copied().ok_or_else(|| fatal("element count"))?;
            let mut dims = Vec::with_capacity(n);
            for i in 0..n {
                let v: Vec<i64> = nums(next(&mut lines)?);
                if v.len() < 5 || v[0] as usize != i + 1 {
                    return Err(SolverError::Unsupported(format!("{CALLER}: element permutations")));
                }
                let t = GMSH_TO_ELMER.get(v[1] as usize - 1).copied().unwrap_or(0);
                if t == 0 {
                    return Err(fatal(format!("unknown Gmsh element type {}", v[1])));
                }
                let nn = elmer_fe::element_type(t).map_err(|e| fatal(e.to_string()))?.n_nodes;
                // `READ(str,*) j,GmshType,j,j,j, NodeIndexes(1:n)`.
                let nodes: Vec<u32> = v[5..5 + nn].iter().map(|&k| k as u32 - 1).collect();
                let ed = if t > 500 { 3 } else if t > 300 { 2 } else if t > 200 { 1 } else { 0 };
                max_dim = max_dim.max(ed);
                dims.push(ed);
                d.elems.push((t, nodes));
            }
            // `IF( ElemDim == MaxElemDim ) NoBulkElems = k`.
            d.n_bulk = dims.iter().rposition(|&e| e == max_dim).map_or(0, |k| k + 1);
        } else if s.starts_with("$NodeData") {
            next(&mut lines)?;
            let name = next(&mut lines)?.trim().trim_matches('"').to_string();
            next(&mut lines)?;
            next(&mut lines)?;
            mesh_dim = nums(next(&mut lines)?).first().copied().unwrap_or(0);
            next(&mut lines)?;
            let dofs: usize = nums(next(&mut lines)?).first().copied().unwrap_or(1);
            let n: usize = nums(next(&mut lines)?).first().copied().unwrap_or(0);
            if dofs != 1 {
                return Err(SolverError::Unsupported(format!("{CALLER}: multi-component field {name}")));
            }
            let mut vals = vec![0.0f64; d.coords.len()];
            for i in 0..n {
                let v: Vec<f64> = nums(next(&mut lines)?);
                if i < vals.len() {
                    vals[i] = v.get(1).copied().unwrap_or(0.0);
                }
            }
            d.fields.push((elmer_core::normalize_key(&name), vals));
        }
    }
    d.dim = max_dim.max(mesh_dim);
    Ok(d)
}

impl Solver for GmshOutputReader {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        if p.integer("filename partitions").is_some_and(|n| n > 1) {
            return Err(SolverError::Unsupported(format!("{CALLER}: Filename Partitions")));
        }
        let mut base = p.string("filename").ok_or_else(|| fatal("Filename not given"))?.to_string();
        if !base.contains('.') {
            base.push_str(".msh");
        }
        // `SolverOutputDirectory(..., UseMeshDir = .TRUE.)`.
        let dir = match p.string("output directory") {
            Some(d) => self.case_dir.join(d),
            None => self.case_dir.join(ctx.sim.model.header.string("mesh name").unwrap_or(".")),
        };
        let path = dir.join(&base);
        let text = std::fs::read_to_string(&path).map_err(|e| fatal(format!("Could not open file: {}: {e}", path.display())))?;
        let mut g = read_gmsh(&text)?;

        // `Mask Name`: the nodes of the boundary elements whose BC has it true.
        let n_nodes = ctx.mesh.n_nodes();
        let mask = match p.string("mask name") {
            Some(m) => {
                let model = &ctx.sim.model;
                if model.bodies.iter().chain(&model.body_forces).any(|s| s.get(m).is_some()) {
                    return Err(SolverError::Unsupported(format!("{CALLER}: Mask Name on bodies")));
                }
                let mut perm = vec![NO_DOF; n_nodes];
                let mut k = 0;
                let nb = ctx.mesh.n_bulk();
                for e in nb..nb + ctx.mesh.n_boundary() {
                    let e = ElemId(e as u32);
                    if ctx.bc_params(e).is_some_and(|bc| bc.logical(m) == Some(true)) {
                        for &n in ctx.mesh.element(e).nodes {
                            if perm[n as usize] == NO_DOF {
                                perm[n as usize] = k;
                                k += 1;
                            }
                        }
                    }
                }
                if k == 0 {
                    return Ok(());
                }
                Some(perm)
            }
            None => None,
        };

        // `Align Coordinate`: the read mesh ends where this one begins.
        let align = p.integer("align coordinate").unwrap_or(0);
        if align != 0 {
            let k = align.unsigned_abs() as usize - 1;
            if k > 2 {
                return Err(fatal(format!("Invalid value for \"Align Coordinate\": {align}")));
            }
            let from: Vec<f64> = g.coords.iter().map(|c| c[k]).collect();
            let to: Vec<f64> = (0..n_nodes as u32).map(|n| ctx.mesh.coord(n)[k]).collect();
            let min = |v: &[f64]| v.iter().copied().fold(f64::MAX, f64::min);
            let max = |v: &[f64]| v.iter().copied().fold(f64::MIN, f64::max);
            let dx = if align > 0 { min(&to) - max(&from) } else { min(&from) - max(&to) };
            for c in &mut g.coords {
                c[k] += dx;
            }
        }

        let from = OwnedMesh {
            dim: g.dim,
            coords: g.coords.clone(),
            elems: g
                .elems
                .iter()
                .map(|(t, nodes)| Elem { etype: ElementType(*t), nodes: nodes.clone(), tag: 1, parents: [None, None] })
                .collect(),
            n_bulk: g.n_bulk,
        };
        let opts = InterpOptions::from_params(&ctx.sim.model.simulation, Some(&p));
        let loc = Locator::new(&from, opts).map_err(|e| fatal(e.to_string()))?;
        let fe = ctx
            .sim
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
        let extras = elmer_interp::Extras { body_mask: None, new_mask_perm: mask.as_deref(), old_dg: None, new_dg: None };
        let projector = build_projector(&loc, ctx.mesh, fe.basis.as_ref(), &extras, ctx.log).map_err(|e| fatal(e.to_string()))?;
        // The reader's `Perm(i) = i`, given explicitly so the projector
        // writes through the target variable's own permutation.
        let ident: Vec<u32> = (0..g.coords.len() as u32).collect();
        for (name, vals) in &g.fields {
            let Some(id) = ctx.sim.vars.id(name) else { continue };
            let v = ctx.sim.vars.by_id_mut(id);
            if v.kind != FieldKind::Nodal || v.dofs != 1 {
                return Err(SolverError::Unsupported(format!("{CALLER}: {name} is not a scalar nodal variable")));
            }
            let old = FieldData { name: name.clone(), location: Location::Nodes, values: vals.clone(), perm: Some(ident.clone()), prev_values: Vec::new() };
            let mut new = vec![FieldData {
                name: name.clone(),
                location: Location::Nodes,
                values: v.values.clone(),
                perm: v.perm.as_deref().map(|q| q.to_vec()),
                prev_values: Vec::new(),
            }];
            projector.apply(std::slice::from_ref(&old), &mut new);
            v.values.copy_from_slice(&new[0].values);
        }
        Ok(())
    }
}
