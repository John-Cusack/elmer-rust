//! `StressSolve` / `StressSolver` as an `elmer_solver::Solver`.
//!
//! One steady linear pass: assemble the bulk and boundary systems with the
//! kernels of this crate, apply Dirichlet conditions, solve. Everything this
//! crate does not implement is refused up front instead of being dropped.

use std::sync::OnceLock;

use elmer_core::{
    ElemId, ElementBasis, ElementType, FieldRef, LocalSystem, MeshView, ParamMap, Params, Quadrature, Value,
};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::boundary::{stress_boundary_local, Spring, StressBoundaryInput};
use crate::bulk::{stress_bulk_local, CoordSystem, StressBulkInput};
use crate::material::{diagonal_c, ElasticLaw, Tensor3};

const CALLER: &str = "StressSolver";

/// Finite-element services, until `elmer-fe` provides them to every solver.
pub struct FeServices {
    pub basis: Box<dyn ElementBasis + Send + Sync>,
    pub quad: Box<dyn Quadrature + Send + Sync>,
}

static FE: OnceLock<FeServices> = OnceLock::new();

/// Install the process-wide basis and quadrature used by registry-built
/// `StressSolver`s. Returns the services back if already installed.
pub fn install_fe(fe: FeServices) -> Result<(), FeServices> {
    FE.set(fe)
}

/// Register `("StressSolve", "StressSolver")`.
pub fn register(reg: &mut SolverRegistry) {
    let (lib, entry) = crate::REGISTRY_ENTRY;
    reg.register(lib, entry, factory);
}

fn factory(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if FE.get().is_none() {
        return Err(SolverError::NotWired("element basis and quadrature (elmer-fe)".into()));
    }
    Ok(Box::new(StressSolver { fe: None }))
}

/// The linear elasticity solver.
pub struct StressSolver {
    fe: Option<FeServices>,
}

impl StressSolver {
    /// A solver that owns its basis and quadrature.
    pub fn new(fe: FeServices) -> Self {
        Self { fe: Some(fe) }
    }

    fn fe(&self) -> &FeServices {
        self.fe.as_ref().or_else(|| FE.get()).expect("checked by factory")
    }
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

/// A constant real keyword as `(rows, cols, data)`; a dependent one is not
/// evaluable here (it needs `elmer-expr`).
type RealArray<'a> = (usize, usize, &'a [f64]);

fn real_array<'a>(p: &'a ParamMap, key: &str) -> Result<Option<RealArray<'a>>, SolverError> {
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) => Ok(Some((*rows, *cols, data))),
        Some(Value::Dependent(_)) => Err(SolverError::NotWired(format!("dependent value for '{key}' (elmer-expr)"))),
        Some(_) => Err(SolverError::fatal(CALLER, format!("'{key}' is not a real"))),
    }
}

fn real(p: &ParamMap, key: &str) -> Result<Option<f64>, SolverError> {
    Ok(real_array(p, key)?.map(|(_, _, d)| d[0]))
}

fn flag(p: Option<&ParamMap>, key: &str) -> bool {
    p.and_then(|p| p.logical(key)).unwrap_or(false)
}

fn present_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

/// Model-wide refusals, from the solver section.
fn check_solver_params(p: &ParamMap) -> Result<(), SolverError> {
    for (key, what) in [
        ("large deflection", "large deflection"),
        ("incompressible", "incompressible formulation"),
        ("stability analysis", "stability analysis"),
        ("eigen analysis", "eigen analysis"),
        ("calculate stresses", "stress postprocessing"),
        ("calculate strains", "strain postprocessing"),
        ("calculate principal", "principal stress postprocessing"),
        ("calculate pangle", "principal stress postprocessing"),
        ("constant bulk system", "constant bulk system"),
        ("model lumping", "model lumping"),
        ("update reference displacement", "staged construction"),
        ("legacy assembly", "legacy assembly"),
        ("geometric stiffness", "geometric stiffness"),
        ("quasi stationary", "quasi stationary"),
    ] {
        if flag(Some(p), key) {
            return Err(unsupported(what));
        }
    }
    if flag(Some(p), "calculate loads") {
        return Err(unsupported("load postprocessing"));
    }
    Ok(())
}

struct ElementData {
    law: ElasticLaw,
    density: Vec<f64>,
    load: Vec<[f64; 3]>,
    pressure: Option<Vec<f64>>,
    temperature: Option<Vec<f64>>,
    heat_expansion: Option<Vec<[f64; 3]>>,
}

fn nodal_temperature(temp: &Option<FieldRef<'_>>, nodes: &[u32], reference: &[f64]) -> Option<Vec<f64>> {
    let t = temp.as_ref()?;
    let perm = t.perm?;
    Some(
        nodes
            .iter()
            .zip(reference)
            .map(|(&n, r)| match perm.get(n as usize) {
                Some(&s) if s != elmer_core::NO_DOF => t.values[s as usize] - r,
                _ => 0.0,
            })
            .collect(),
    )
}

fn parent_center(basis: &dyn ElementBasis, mesh: &dyn MeshView, parent: ElemId) -> Result<[f64; 3], SolverError> {
    let pe = mesh.element(parent);
    let uvw = match pe.etype.family() {
        3 | 7 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        5 => [0.25, 0.25, 0.25],
        6 => [0.0, 0.0, 1.0 / 3.0],
        _ => [0.0; 3],
    };
    let xyz: Vec<[f64; 3]> = pe.nodes.iter().map(|&n| mesh.coord(n)).collect();
    let mut be = elmer_core::BasisEval::default();
    basis
        .eval(pe.etype, &xyz, uvw, &mut be)
        .map_err(|e| SolverError::fatal(CALLER, e))?;
    Ok(std::array::from_fn(|a| (0..xyz.len()).map(|k| be.basis[k] * xyz[k][a]).sum()))
}

impl StressSolver {
    #[allow(clippy::too_many_arguments)]
    fn element_data(
        &self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        cdim: usize,
        axi: bool,
        temp: &Option<FieldRef<'_>>,
    ) -> Result<ElementData, SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let bp = ctx
            .body_params(e)
            .ok_or_else(|| SolverError::fatal(CALLER, format!("element {} has no body", e.0)))?;
        let mat = bp
            .material
            .ok_or_else(|| SolverError::fatal(CALLER, format!("body {} has no material", bp.body_id)))?;
        for (key, what) in [
            ("maxwell material", "Maxwell material"),
            ("umat subroutine", "UMAT"),
            ("youngs modulus at ip", "properties at integration points"),
            ("poisson ratio at ip", "properties at integration points"),
            ("heat expansion coefficient ip", "properties at integration points"),
            ("rayleigh damping alpha", "Rayleigh damping"),
            ("rayleigh damping beta", "Rayleigh damping"),
        ] {
            if mat.get(key).is_some() && (flag(Some(mat), key) || !matches!(mat.get(key), Some(Value::Logical(_)))) {
                return Err(unsupported(what));
            }
        }
        if let Some(bf) = bp.body_force {
            for prefix in [
                "stress load",
                "strain load",
                "inertial bodyforce",
                "stress volume source",
                "gravitational prestress",
                "stress bodyforce at ip",
            ] {
                if present_prefix(bf, prefix) {
                    return Err(unsupported(&format!("body force '{prefix}'")));
                }
            }
        }

        let ym = real_array(mat, "youngs modulus")?
            .ok_or_else(|| SolverError::fatal(CALLER, "material has no 'Youngs Modulus'"))?;
        let law = if ym.0 == 1 && ym.1 == 1 {
            let nu = real(mat, "poisson ratio")?
                .ok_or_else(|| SolverError::fatal(CALLER, "material has no 'Poisson Ratio'"))?;
            ElasticLaw::Isotropic { young: vec![ym.2[0]; n], poisson: vec![nu; n] }
        } else {
            if cdim != 3 || axi {
                return Err(unsupported("anisotropic material outside 3D (ElasticSolve refuses it too)"));
            }
            let (rows, cols, d) = ym;
            let mut c = diagonal_c(&[]);
            if rows == 1 || cols == 1 {
                c = diagonal_c(d);
            } else {
                for i in 0..rows.min(6) {
                    for j in 0..cols.min(6) {
                        c[i][j] = d[i * cols + j];
                    }
                }
            }
            let rotation = if flag(Some(mat), "rotate elasticity tensor") {
                let mut t: Tensor3 = [[0.0; 3]; 3];
                for (i, row) in t.iter_mut().enumerate() {
                    let (_, _, v) = real_array(mat, &format!("material coordinates unit vector {}", i + 1))?
                        .ok_or_else(|| {
                            SolverError::fatal(CALLER, "'Rotate Elasticity Tensor' without unit vectors")
                        })?;
                    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    if norm < f64::EPSILON {
                        return Err(SolverError::fatal(CALLER, "material coordinate unit vector too short"));
                    }
                    for (a, x) in row.iter_mut().enumerate() {
                        *x = v[a] / norm;
                    }
                }
                Some(t)
            } else {
                None
            };
            ElasticLaw::Anisotropic { c: vec![c; n], rotation }
        };

        let density = vec![real(mat, "density")?.unwrap_or(0.0); n];
        let mut load = vec![[0.0; 3]; n];
        let mut pressure = None;
        if let Some(bf) = bp.body_force {
            for i in 0..cdim {
                let v = real(bf, &format!("stress bodyforce {}", i + 1))?.unwrap_or(0.0);
                for l in load.iter_mut() {
                    l[i] = v;
                }
            }
            if let Some(p) = real(bf, "stress pressure")? {
                pressure = Some(vec![p; n]);
            }
        }

        let reference = vec![real(mat, "reference temperature")?.unwrap_or(0.0); n];
        let temperature = nodal_temperature(temp, el.nodes, &reference);
        let heat_expansion = match real_array(mat, "heat expansion coefficient")? {
            Some((r, c, d)) if temperature.is_some() => {
                let diag: [f64; 3] = std::array::from_fn(|i| {
                    if r == 1 && c == 1 {
                        d[0]
                    } else if r == 1 || c == 1 {
                        d.get(i).copied().unwrap_or(0.0)
                    } else if i < r && i < c {
                        d[i * c + i]
                    } else {
                        0.0
                    }
                });
                Some(vec![diag; n])
            }
            _ => None,
        };
        Ok(ElementData { law, density, load, pressure, temperature, heat_expansion })
    }
}

impl Solver for StressSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_solver_params(ctx.params())?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map(|f| f.dofs);
        if dofs != Some(dim) {
            return Err(SolverError::fatal(
                CALLER,
                format!("'Variable DOFs' must equal the coordinate dimension {dim}, got {dofs:?}"),
            ));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient simulation (mass and second-order time are not wired)"));
        }
        check_solver_params(ctx.params())?;
        let cdim = ctx.mesh.dim();
        let coord_name = elmer_core::normalize_key(ctx.simulation_params().string("coordinate system").unwrap_or("cartesian"));
        let axi = match coord_name.as_str() {
            "axi symmetric" => true,
            s if s.starts_with("cartesian") => false,
            other => return Err(unsupported(&format!("coordinate system '{other}'"))),
        };
        let coord = if axi { CoordSystem::AxiSymmetric } else { CoordSystem::Cartesian };

        ctx.default_initialize()?;

        // Snapshot the temperature field (read-only) before the mutable assembly.
        let temp_owned: Option<(Vec<u32>, Vec<f64>)> = ctx
            .variable("temperature")
            .and_then(|f| f.perm.map(|p| (p.to_vec(), f.values.to_vec())));
        let temp_ref = temp_owned.as_ref().map(|(p, v)| FieldRef {
            name: "temperature",
            kind: elmer_core::FieldKind::Nodal,
            dofs: 1,
            perm: Some(p),
            values: v,
        });

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        for e in elems {
            let el = ctx.mesh.element(e);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let d = self.element_data(ctx, e, cdim, axi, &temp_ref)?;
            let plane_stress = flag(ctx.body_params(e).and_then(|b| b.equation), "plane stress");
            let inp = StressBulkInput {
                etype: el.etype,
                xyz: &xyz,
                cdim,
                coord,
                plane_stress,
                law: &d.law,
                density: &d.density,
                damping: None,
                load: &d.load,
                pressure_load: d.pressure.as_deref(),
                stress_load: None,
                strain_load: None,
                temperature: d.temperature.as_deref(),
                heat_expansion: d.heat_expansion.as_deref(),
                want_mass: false,
                np: None,
            };
            let mut local = LocalSystem::zeros(0);
            let fe = self.fe();
            stress_bulk_local(&*fe.basis, &*fe.quad, &inp, &mut local).map_err(|x| SolverError::fatal(CALLER, x.to_string()))?;
            ctx.default_update(&local, e)?;
        }

        let belems: Vec<ElemId> = ctx.boundary_elements().collect();
        for e in belems {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let bc = ctx.bc_params(e).expect("boundary_elements filters on BC").clone();
            if let Some(v) = bc.logical("fsi bc") {
                if v {
                    return Err(unsupported("FSI boundary condition"));
                }
            }
            for (prefix, what) in [
                ("stress load", "boundary 'Stress Load'"),
                ("model lumping", "model lumping boundary"),
                ("normal-tangential", "normal-tangential displacement"),
            ] {
                if present_prefix(&bc, prefix) {
                    return Err(unsupported(what));
                }
            }
            let mut traction = [0.0; 3];
            let mut got = false;
            for prefix in ["surface traction", "force"] {
                if present_prefix(&bc, prefix) {
                    for (i, t) in traction.iter_mut().enumerate().take(cdim) {
                        *t = real(&bc, &format!("{prefix} {}", i + 1))?.unwrap_or(0.0);
                    }
                    got = true;
                    break;
                }
            }
            let mut beta = real(&bc, "normal surface traction")?;
            if beta.is_none() {
                beta = real(&bc, "normal force")?;
            }
            got |= beta.is_some();

            let mut spring_normal = None;
            let mut spring_tensor = None;
            if present_prefix(&bc, "spring") {
                if let Some(k) = real(&bc, "spring")? {
                    spring_normal = Some(k);
                } else {
                    let mut t = [[0.0; 3]; 3];
                    for (i, row) in t.iter_mut().enumerate().take(cdim) {
                        row[i] = real(&bc, &format!("spring {}", i + 1))?.unwrap_or(0.0);
                        for (j, v) in row.iter_mut().enumerate().take(cdim) {
                            if let Some(x) = real(&bc, &format!("spring {}{}", i + 1, j + 1))? {
                                *v = x;
                            }
                        }
                    }
                    spring_tensor = Some(t);
                }
                got = true;
            }
            if !got {
                continue;
            }

            let parent = el.parents[0]
                .or(el.parents[1])
                .ok_or_else(|| SolverError::fatal(CALLER, format!("boundary element {} has no parent", e.0)))?;
            let fe = self.fe();
            let center = parent_center(&*fe.basis, ctx.mesh, parent)?;
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
            let tr = vec![traction; n];
            let nt = vec![beta.unwrap_or(0.0); n];
            let sn = vec![spring_normal.unwrap_or(0.0); n];
            let st = vec![spring_tensor.unwrap_or_default(); n];
            let spring = if spring_normal.is_some() {
                Spring::Normal(&sn)
            } else if spring_tensor.is_some() {
                Spring::Tensor(&st)
            } else {
                Spring::None
            };
            let inp = StressBoundaryInput {
                etype: ElementType(el.etype.0),
                xyz: &xyz,
                cdim,
                coord,
                traction: &tr,
                normal_traction: &nt,
                spring,
                parent_center: center,
                np: None,
            };
            let mut local = LocalSystem::zeros(0);
            stress_boundary_local(&*fe.basis, &*fe.quad, &inp, &mut local)
                .map_err(|x| SolverError::fatal(CALLER, x.to_string()))?;
            ctx.default_update(&local, e)?;
        }

        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}
