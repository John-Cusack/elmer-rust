//! `DataToFieldSolver` (fem/src/modules/DataToFieldSolver.F90): a smooth
//! nodal field `Fit` that follows given data, the target variable's nodal
//! values or scattered points read from a file, regularised by diffusion:
//!
//! ```text
//!   (w + D grad . grad) fit = w * data
//! ```
//!
//! with the nodal weight `w` the integral of each node's basis function
//! (or a given weight variable), optional masking of the data and of the
//! diffusion, and Dirichlet values from the boundary conditions. With
//! `Continue BC With Boundary Assembly`, boundaries flagged `<Fit> Continue`
//! get the boundary term of the diffusion instead of a natural condition.
//!
//! [`init_section`] is `DataToFieldSolver_init`, the section edits made
//! before the solver's variable exists. Refused: `Logarithmic Fitting`
//! reverted with `Revert Logarithmic Fitting`, a `Data Diffusivity` that is
//! not constant, and the boundary assembly with a material diffusivity or on
//! boundary elements other than lines.
#![forbid(unsafe_code)]

use elmer_core::{ElementType, LocalSystem, MeshView, ParamMap, Params, Value, NO_DOF};
use elmer_fe::{ElmerQuadrature, LagrangeBasis};
use elmer_interp::geom::global_to_local;
use elmer_interp::shape::Shape;
use elmer_interp::{InterpOptions, Locator};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use elmer_core::{BasisEval, ElementBasis, Quadrature};

const CALLER: &str = "DataToFieldSolver";

/// `Procedure = "DataToFieldSolver" "DataToFieldSolver"`.
pub const PROCEDURE: (&str, &str) = ("DataToFieldSolver", "DataToFieldSolver");

pub fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, factory);
}

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(DataToField))
}

/// Whether `params` is a DataToFieldSolver section.
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params)
        .is_some_and(|(lib, entry)| lib.eq_ignore_ascii_case(PROCEDURE.0) && entry.eq_ignore_ascii_case(PROCEDURE.1))
}

fn add_new(params: &mut ParamMap, key: &str, v: Value) {
    if params.get(key).is_none() {
        params.insert(key, v);
    }
}

/// `NextFreeKeyword('Exported Variable', Params)`.
fn next_exported(params: &ParamMap) -> String {
    let n = (1..).find(|i| params.get(&format!("exported variable {i}")).is_none()).unwrap_or(1);
    format!("Exported Variable {n}")
}

/// `DataToFieldSolver_init`: names the variables and fills in the linear
/// solver defaults.
pub fn init_section(params: &mut ParamMap) {
    let have_file = params.get("point data filename").is_some();
    let given = params.logical("normalize by given weight") == Some(true)
        || params.logical("normalize data by weight") == Some(true);
    if params.get("target variable").is_none() && have_file {
        let key = next_exported(params);
        params.insert(&key, Value::Str("point data".into()));
        params.insert("Target Variable", Value::Str("point data".into()));
    }
    if given && have_file && params.get("weight variable").is_none() {
        let key = next_exported(params);
        params.insert(&key, Value::Str("point weight".into()));
        params.insert("Weight Variable", Value::Str("point weight".into()));
    }
    if params.get("variable").is_none() {
        let name = if have_file { "Point Fit" } else { "Fit" };
        params.insert("Variable", Value::Str(name.into()));
    }
    params.insert("Time derivative order", Value::Integer(vec![0]));
    add_new(params, "Linear System Solver", Value::Str("iterative".into()));
    add_new(params, "Linear System Iterative Method", Value::Str("bicgstab".into()));
    add_new(params, "Linear System Preconditioning", Value::Str("none".into()));
    add_new(params, "Linear System Max Iterations", Value::Integer(vec![1000]));
    add_new(params, "Linear System Residual Output", Value::Integer(vec![20]));
    add_new(params, "Linear System Convergence Tolerance", Value::real(1.0e-10));
}

struct DataToField;

/// A field's values through its permutation (`Var % Perm`, identity when
/// absent).
struct Field {
    values: Vec<f64>,
    perm: Option<Vec<u32>>,
}

impl Field {
    fn read(ctx: &SolverCtx<'_>, name: &str) -> Result<Self, SolverError> {
        let f = ctx.variable(name).ok_or_else(|| SolverError::fatal(CALLER, format!("Variable not present: {name}")))?;
        Ok(Field { values: f.values.to_vec(), perm: f.perm.map(<[u32]>::to_vec) })
    }

    fn slot(&self, node: usize) -> Option<usize> {
        match &self.perm {
            None => Some(node),
            Some(p) => p.get(node).filter(|&&s| s != NO_DOF).map(|&s| s as usize),
        }
    }
}

/// `Data Diffusivity` (or `Diffusivity Name`) of a material as the 3x3
/// matrix the Fortran builds from a scalar, a diagonal or a full array.
fn diffusivity(m: &ParamMap, name: &str) -> Result<Option<[[f64; 3]; 3]>, SolverError> {
    let Some(v) = m.get(name) else { return Ok(None) };
    let Value::Real { rows, cols, data } = v else {
        return Err(SolverError::Unsupported(format!("{CALLER}: {name} that is not a constant")));
    };
    let mut d = [[0.0; 3]; 3];
    if data.len() == 1 {
        for (i, row) in d.iter_mut().enumerate() {
            row[i] = data[0];
        }
    } else if *cols == 1 || *rows == 1 {
        for i in 0..data.len().min(3) {
            d[i][i] = data[i];
        }
    } else {
        for i in 0..(*rows).min(3) {
            for j in 0..(*cols).min(3) {
                d[i][j] = data[i * cols + j];
            }
        }
    }
    Ok(Some(d))
}

fn element_xyz(mesh: &dyn MeshView, nodes: &[u32]) -> Vec<[f64; 3]> {
    nodes.iter().map(|&n| mesh.coord(n)).collect()
}

/// `NormalVector(Element, Nodes, u, v, Check=.TRUE.)` of a line element in
/// 2D, pointing away from its parent.
fn line_normal(xyz: &[[f64; 3]], parent_xyz: &[[f64; 3]], point: [f64; 3]) -> [f64; 3] {
    let t = [xyz[1][0] - xyz[0][0], xyz[1][1] - xyz[0][1]];
    let len = (t[0] * t[0] + t[1] * t[1]).sqrt();
    let mut n = [t[1] / len, -t[0] / len, 0.0];
    let c = parent_xyz.iter().fold([0.0; 3], |a, x| [a[0] + x[0], a[1] + x[1], a[2] + x[2]]);
    let k = parent_xyz.len() as f64;
    let out = [point[0] - c[0] / k, point[1] - c[1] / k];
    if n[0] * out[0] + n[1] * out[1] < 0.0 {
        n = [-n[0], -n[1], 0.0];
    }
    n
}

impl Solver for DataToField {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, "-----------------------------------------", 4);
        ctx.info(CALLER, "Resolving field from given data", 4);
        let p = ctx.params().clone();
        let fit_name = ctx.variable_name().ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?.to_string();
        let fit = Field::read(ctx, &fit_name)?;
        let field_perm = fit.perm.clone().unwrap_or_else(|| (0..ctx.mesh.n_nodes() as u32).collect());

        let target = p
            .string("data field name")
            .or_else(|| p.string("target variable"))
            .ok_or_else(|| SolverError::fatal(CALLER, "> Target Variable < must exist for the solver!"))?
            .to_string();
        let mut force = Field::read(ctx, &target)?;
        let global_weight = p.real("weight coefficient");
        let given = p.logical("normalize by given weight") == Some(true) || p.logical("normalize data by weight") == Some(true);
        let mut weight = None;
        let mut constant_weight_sum = false;
        if given {
            let name = p
                .string("weight field name")
                .or_else(|| p.string("weight variable"))
                .ok_or_else(|| SolverError::fatal(CALLER, "> Weight Variable < must exist for the solver!"))?;
            weight = Some((name.to_string(), Field::read(ctx, name)?));
            constant_weight_sum = p.logical("set constant weight sum") == Some(true);
        }
        let nodal_normalize = p.logical("normalize by nodal weight") == Some(true);
        if given && nodal_normalize {
            return Err(SolverError::fatal(CALLER, "Normalization cannot be both given and nodal!"));
        }
        let mask_name = p.string("mask field name").or_else(|| p.string("mask variable"));
        let mask = match mask_name {
            Some(name) => Some(Field::read(ctx, name)?),
            None => None,
        };
        let max_mask = p.real("max mask value");
        let min_mask = p.real("min mask value");
        let (min_mask, max_mask) = (
            min_mask.unwrap_or(if max_mask.is_some() { f64::MIN } else { 0.0 }),
            max_mask.unwrap_or(f64::MAX),
        );
        let mask_diffusion = mask.is_some() && p.logical("mask diffusion") == Some(true);
        let outside = |v: f64| v < min_mask || v > max_mask;

        let dim = ctx.mesh.dim();
        let basis = LagrangeBasis::new(dim);
        let quad = ElmerQuadrature;
        let mut ev = BasisEval::default();

        if let Some(file) = p.string("point data filename") {
            let column = p.integer("point data column").map_or(dim + 1, |c| c.max(1) as usize);
            let text = std::fs::read_to_string(file)
                .map_err(|e| SolverError::fatal(CALLER, format!("Could not open file for reading: {file}: {e}")))?;
            let locator = Locator::new(ctx.mesh, InterpOptions::default()).map_err(|e| SolverError::fatal(CALLER, e.to_string()))?;
            let mut count = 0;
            for line in text.lines() {
                let vals: Vec<f64> = line.split_whitespace().filter_map(|t| t.replace(['d', 'D'], "e").parse().ok()).collect();
                if vals.len() < column {
                    break;
                }
                count += 1;
                let mut x = [0.0; 3];
                x[..dim].copy_from_slice(&vals[..dim]);
                let val = vals[column - 1];
                let Some(at) = locator.locate(locator.search_point(x), None) else { continue };
                let el = ctx.mesh.element(at.element);
                let xyz = element_xyz(ctx.mesh, el.nodes);
                basis.eval(el.etype, &xyz, at.uvw, &mut ev).map_err(|e| SolverError::fatal(CALLER, e))?;
                for (i, &node) in el.nodes.iter().enumerate() {
                    let j = field_perm[node as usize];
                    if j == NO_DOF {
                        continue;
                    }
                    force.values[j as usize] += ev.basis[i] * val;
                    if let Some((_, w)) = weight.as_mut() {
                        w.values[j as usize] += ev.basis[i];
                    }
                }
            }
            if count == 0 {
                return Err(SolverError::fatal(CALLER, "Could not read any points from file!"));
            }
            // Written back: the Fortran accumulates into the variables themselves.
            if let Some(v) = ctx.variable_mut(&target) {
                v.values.copy_from_slice(&force.values);
            }
            if let Some((name, w)) = &weight {
                if let Some(v) = ctx.variable_mut(name) {
                    v.values.copy_from_slice(&w.values);
                }
            }
        }
        let log_scale = p.logical("logarithmic fitting") == Some(true);
        if p.logical("revert logarithmic fitting") == Some(true) {
            return Err(SolverError::Unsupported(format!("{CALLER}: Revert Logarithmic Fitting")));
        }
        let continue_bc = p.logical("continue bc with boundary assembly") == Some(true);
        let cond_name = format!("{fit_name} continue");

        ctx.default_initialize()?;

        // BulkAssembly.
        let coeff = p.real("diffusion coefficient");
        let diff_name = p.string("diffusivity name").unwrap_or("data diffusivity").to_string();
        let n_nodes = ctx.mesh.n_nodes();
        let mut boundary_skip = vec![false; n_nodes];
        if !continue_bc {
            for e in ctx.boundary_elements().collect::<Vec<_>>() {
                if ctx.bc_params(e).is_some_and(|bc| bc.logical(&cond_name) == Some(true)) {
                    for &node in ctx.mesh.element(e).nodes {
                        boundary_skip[node as usize] = true;
                    }
                }
            }
        }
        let mut nodal_weight = vec![0.0; n_nodes];
        for e in ctx.active_elements().collect::<Vec<_>>() {
            let el = ctx.mesh.element(e);
            let nodes = el.nodes.to_vec();
            let etype = el.etype;
            let n = nodes.len();
            let xyz = element_xyz(ctx.mesh, &nodes);
            let local_d = match coeff {
                Some(_) => None,
                None => match ctx.body_params(e).and_then(|b| b.material) {
                    Some(m) => diffusivity(m, &diff_name)?,
                    None => None,
                },
            };
            let mut local = LocalSystem::zeros(n);
            for qp in quad.gauss_points(etype, None) {
                basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|e| SolverError::fatal(CALLER, e))?;
                let s = qp.weight * ev.det_j;
                for i in 0..n {
                    let node = nodes[i] as usize;
                    if !continue_bc && boundary_skip[node] {
                        continue;
                    }
                    nodal_weight[node] += s * ev.basis[i];
                    if mask_diffusion {
                        let m = mask.as_ref().unwrap();
                        if let Some(slot) = m.slot(node) {
                            if !outside(m.values[slot]) {
                                continue;
                            }
                        }
                    }
                    for j in 0..n {
                        let gi = ev.dbasis_dx[i];
                        let gj = ev.dbasis_dx[j];
                        let k = if let Some(c) = coeff {
                            c * (0..dim).map(|a| gi[a] * gj[a]).sum::<f64>()
                        } else if let Some(d) = local_d {
                            (0..dim).map(|a| (0..dim).map(|b| d[a][b] * gj[b]).sum::<f64>() * gi[a]).sum::<f64>()
                        } else {
                            0.0
                        };
                        local.stiff[i * n + j] += s * k;
                    }
                }
            }
            ctx.default_update(&local, e)?;
        }

        let weight_corr = match (&weight, constant_weight_sum) {
            (Some((_, w)), true) => nodal_weight.iter().sum::<f64>() / w.values.iter().sum::<f64>(),
            _ => 1.0,
        };
        let mut mask_active = 0;
        let mut diag_add = Vec::new();
        for i in 0..n_nodes {
            if let Some(m) = &mask {
                if let Some(slot) = m.slot(i) {
                    if outside(m.values[slot]) {
                        mask_active += 1;
                        continue;
                    }
                }
            }
            let j = field_perm[i];
            if j == NO_DOF {
                continue;
            }
            let mut val = force.slot(i).map_or(0.0, |s| force.values[s]);
            if log_scale {
                val = val.ln();
            }
            let (mut wmat, mut wrhs) = if let Some((_, w)) = &weight {
                (w.slot(i).map_or(0.0, |s| weight_corr * w.values[s]), weight_corr * val)
            } else if nodal_normalize {
                (nodal_weight[i], val)
            } else {
                (nodal_weight[i], nodal_weight[i] * val)
            };
            if let Some(g) = global_weight {
                wmat *= g;
                wrhs *= g;
            }
            diag_add.push((j as usize, wmat, wrhs));
        }
        {
            let sys = ctx.system_mut().ok_or_else(|| SolverError::fatal(CALLER, "no linear system"))?;
            for (j, wmat, wrhs) in diag_add {
                let m = &mut sys.matrix;
                let (lo, hi) = (m.row_ptr[j], m.row_ptr[j + 1]);
                let k = (lo..hi)
                    .find(|&k| m.col_idx[k] as usize == j)
                    .ok_or_else(|| SolverError::fatal(CALLER, format!("row {j} has no diagonal")))?;
                m.values[k] += wmat;
                sys.rhs[j] += wrhs;
            }
        }
        if let Some(m) = &mask {
            ctx.info(CALLER, &format!("Mask is active for {mask_active} nodes (out of {})", m.values.len()), 4);
        }
        ctx.default_finish_bulk_assembly()?;

        // BoundaryAssembly.
        if continue_bc {
            let Some(c) = coeff else {
                return Err(SolverError::Unsupported(format!(
                    "{CALLER}: boundary assembly with a material diffusivity"
                )));
            };
            for e in ctx.boundary_elements().collect::<Vec<_>>() {
                if !ctx.bc_params(e).is_some_and(|bc| bc.logical(&cond_name) == Some(true)) {
                    continue;
                }
                let el = ctx.mesh.element(e);
                let [left, right] = el.parents;
                let parent = left.ok_or_else(|| SolverError::fatal(CALLER, "Could not find parent element!"))?;
                if right.is_some() {
                    return Err(SolverError::fatal(CALLER, "This does not make sense for internal BC!"));
                }
                if el.etype.0 / 100 != 2 || dim != 2 {
                    return Err(SolverError::Unsupported(format!("{CALLER}: boundary assembly on element type {}", el.etype.0)));
                }
                let nodes = el.nodes.to_vec();
                let etype = el.etype;
                let xyz = element_xyz(ctx.mesh, &nodes);
                let pel = ctx.mesh.element(parent);
                let pnodes = pel.nodes.to_vec();
                let ptype: ElementType = pel.etype;
                let pxyz = element_xyz(ctx.mesh, &pnodes);
                let pshape = Shape::get(ptype.0).ok_or_else(|| SolverError::fatal(CALLER, format!("element {}", ptype.0)))?;
                let np = pnodes.len();
                let mut local = LocalSystem::zeros(np);
                let mut pev = BasisEval::default();
                for qp in quad.gauss_points(etype, None) {
                    basis.eval(etype, &xyz, qp.uvw, &mut ev).map_err(|e| SolverError::fatal(CALLER, e))?;
                    let s = qp.weight * ev.det_j;
                    let mut x = [0.0; 3];
                    for (b, xi) in ev.basis.iter().zip(&xyz) {
                        for a in 0..3 {
                            x[a] += b * xi[a];
                        }
                    }
                    let nrm = line_normal(&xyz, &pxyz, x);
                    let puvw = global_to_local(pshape, &pxyz, x, dim);
                    basis.eval(ptype, &pxyz, puvw, &mut pev).map_err(|e| SolverError::fatal(CALLER, e))?;
                    for a in 0..np {
                        for b in 0..np {
                            let dn: f64 = (0..dim).map(|k| pev.dbasis_dx[b][k] * nrm[k]).sum();
                            local.stiff[a * np + b] -= s * pev.basis[a] * c * dn;
                        }
                    }
                }
                ctx.default_update(&local, parent)?;
            }
        }
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        ctx.info(CALLER, "All done", 4);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_names_the_point_variables() {
        let mut p = ParamMap::new();
        p.insert("Point Data Filename", Value::Str("vals.dat".into()));
        p.insert("Normalize Data by Weight", Value::Logical(true));
        init_section(&mut p);
        assert_eq!(p.string("variable"), Some("Point Fit"));
        assert_eq!(p.string("target variable"), Some("point data"));
        assert_eq!(p.string("weight variable"), Some("point weight"));
        assert_eq!(p.string("exported variable 1"), Some("point data"));
        assert_eq!(p.string("exported variable 2"), Some("point weight"));
        assert_eq!(p.string("linear system iterative method"), Some("bicgstab"));
    }

    #[test]
    fn diffusivity_shapes() {
        let mut m = ParamMap::new();
        m.insert("Data Diffusivity", Value::Real { rows: 3, cols: 1, data: vec![0.1, 0.2, 0.0] });
        let d = diffusivity(&m, "data diffusivity").unwrap().unwrap();
        assert_eq!((d[0][0], d[1][1], d[0][1]), (0.1, 0.2, 0.0));
    }
}
