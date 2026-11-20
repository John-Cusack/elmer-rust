//! `ElasticSolve` / `ElasticSolver` as an `elmer_solver::Solver`.
//!
//! Each call is the nonlinear loop of `ElasticSolve.F90:1607`: assemble the
//! Newton system at the current displacement, apply the Dirichlet conditions,
//! solve for the new displacement, stop on `Nonlinear System Convergence
//! Tolerance` or `Nonlinear System Max Iterations`. Everything this crate does
//! not implement is refused up front instead of being dropped.

use std::sync::OnceLock;

use elmer_core::{ElemId, ElementBasis, ElementType, FieldKind, FieldRef, LocalSystem, ParamMap, Params, Quadrature, Value};
use elmer_physics_elastic::post::{self, stress_strain_local, StressPostInput};
pub use elmer_physics_elastic::FeServices;
use elmer_physics_elastic::{diagonal_c, CoordSystem, ElasticLaw, IpProps, Spring};
use elmer_solver::sim::VarSpec;
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry, Variable};

use crate::boundary::{elastic_boundary_local, elastic_boundary_local_p, node_uvw, BoundaryInput, ParentInput};
use crate::bulk::{elastic_bulk_local, elastic_bulk_local_p, BulkInput};
use crate::law::Voigt6;
use crate::neohooke::{neo_hookean_local, NeoHookeInput};
use crate::props;
use crate::tensor::T3;

const CALLER: &str = "ElasticSolver";

static FE: OnceLock<FeServices> = OnceLock::new();

/// Install the process-wide basis and quadrature used by registry-built
/// `ElasticSolver`s. Returns the services back if already installed.
pub fn install_fe(fe: FeServices) -> Result<(), FeServices> {
    FE.set(fe)
}

/// Register `("ElasticSolve", "ElasticSolver")`.
pub fn register(reg: &mut SolverRegistry) {
    let (lib, entry) = crate::REGISTRY_ENTRY;
    reg.register(lib, entry, factory);
}

fn factory(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if FE.get().is_none() {
        return Err(SolverError::NotWired("element basis and quadrature (elmer-fe)".into()));
    }
    Ok(Box::new(ElasticSolver { fe: None, mesh_note: false, umat: None }))
}

/// The elasticity solver.
pub struct ElasticSolver {
    fe: Option<FeServices>,
    mesh_note: bool,
    /// A UMAT material's committed integration-point state.
    umat: Option<UmatSaved>,
}

/// `UmatStress0`, `UmatEnergy0`, `UmatState0` and the converged
/// displacement of the previous scanning step (`PrevValues(:,1)`).
#[derive(Default)]
struct UmatSaved {
    stress0: Vec<f64>,
    energy0: Vec<f64>,
    state0: Vec<f64>,
    prev: Vec<f64>,
}

impl ElasticSolver {
    /// A solver that owns its basis and quadrature.
    pub fn new(fe: FeServices) -> Self {
        Self { fe: Some(fe), mesh_note: false, umat: None }
    }

    fn fe(&self) -> &FeServices {
        self.fe.as_ref().or_else(|| FE.get()).expect("checked by factory")
    }
}

/// `RelIntegOrder` / `RelIntegOrderBC` of ElasticSolver: `Relative
/// Integration Order`, overridden by `Bulk Relative Integration Order` for
/// the bulk and `Boundary Relative Integration Order` for the boundary.
#[derive(Copy, Clone, Debug, Default)]
struct RelOrders {
    bulk: i32,
    bc: i32,
}

impl RelOrders {
    fn of(p: &ParamMap) -> Self {
        let base = p.integer("relative integration order").unwrap_or(0);
        Self {
            bulk: p.integer("bulk relative integration order").unwrap_or(base),
            bc: p.integer("boundary relative integration order").unwrap_or(base),
        }
    }
}

/// `GaussPoints(Element, RelOrder = rel)` of a nodal element as a point
/// count: `GaussPoints0` / `GaussPoints2` of elements.def for -1 / +1, the
/// element's own rule (`None`) for 0. Prisms are refused: their rule for a
/// named count is not the one `RelOrder` gives.
fn rel_np(quad: &dyn Quadrature, etype: ElementType, rel: i32) -> Result<Option<usize>, SolverError> {
    if rel == 0 {
        return Ok(None);
    }
    if etype.family() == 7 {
        return Err(unsupported("relative integration order on prisms"));
    }
    match quad.gauss_points_rel(etype, rel) {
        Some(points) => Ok(Some(points.len())),
        None => Err(SolverError::fatal(
            "GaussPoints",
            format!("Relative Integration Order = {rel} but a non-p element only has rules for {{-1,0,1}}"),
        )),
    }
}

/// `AEPS` (Types.F90).
const AEPS: f64 = 10.0 * f64::EPSILON;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {what}"))
}

fn fatal(msg: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, msg.into())
}

fn flag(p: Option<&ParamMap>, key: &str) -> bool {
    p.and_then(|p| p.logical(key)).unwrap_or(false)
}

fn present_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

/// The kinematics and material the solver section selects.
#[derive(Copy, Clone, Debug)]
struct Model {
    /// `Large Deflection`, true unless the sif says otherwise (`:1115`).
    large: bool,
    neo_hookean: bool,
    /// `Elasticity Solver Linear`: assemble at zero displacement.
    linear_model: bool,
    /// `Incompressible`: the linear mixed formulation (`LinearIncompressible`).
    incompressible: bool,
}

fn model_of(p: &ParamMap) -> Model {
    Model {
        large: p.logical("large deflection").unwrap_or(true),
        neo_hookean: flag(Some(p), "neo-hookean material"),
        linear_model: flag(Some(p), "elasticity solver linear"),
        incompressible: flag(Some(p), "incompressible"),
    }
}

/// Model-wide refusals, from the solver section.
fn check_solver_params(p: &ParamMap) -> Result<(), SolverError> {
    for (key, what) in [
        ("mixed formulation", "mixed (pressure) formulation"),
        ("harmonic analysis", "harmonic analysis"),
        ("harmonic mode", "harmonic analysis"),
        ("calculate pangle", "principal angle postprocessing"),
        ("local matrix storage", "local matrix storage"),
        ("adaptive mesh refinement", "adaptive mesh refinement"),
        ("use hencky strain", "Hencky strain"),
        ("maxwell material", "Maxwell material"),
    ] {
        if flag(Some(p), key) {
            return Err(unsupported(what));
        }
    }
    for (key, what) in [
        ("element relative integration order", "elemental integration rule"),
        ("update transient system", "update transient system"),
    ] {
        if p.get(key).is_some() && p.integer(key).unwrap_or(1) != 0 {
            return Err(unsupported(what));
        }
    }
    if p.get("integration rule probe").is_some() {
        return Err(unsupported("integration rule probe"));
    }
    Ok(())
}

/// Which postprocessing fields the model asks for.
#[derive(Copy, Clone, Debug, Default)]
struct PostFlags {
    stress: bool,
    strain: bool,
    principal: bool,
    /// Some equation names `Calculate Stresses`, so only those bodies project.
    masked: bool,
}

impl PostFlags {
    fn any(&self) -> bool {
        self.stress || self.strain
    }
}

fn post_flags(ctx: &SolverCtx<'_>) -> PostFlags {
    let p = ctx.params();
    let masked = ctx.sim.model.equations.iter().any(|e| flag(Some(e), "calculate stresses"));
    PostFlags {
        stress: flag(Some(p), "calculate stresses") || masked,
        strain: flag(Some(p), "calculate strains"),
        principal: flag(Some(p), "calculate principal"),
        masked,
    }
}

/// The postprocessing fields `ElasticSolver_Init` exports
/// (`Stress.F90:1656` `StressFieldDefinition`): name, dofs and component names.
fn field_specs(cdim: usize, f: PostFlags) -> Vec<VarSpec> {
    let comp = ["xx", "yy", "zz", "xy", "yz", "xz"];
    let n = post::output_components(cdim);
    let spec = |name: &str, components: Vec<String>| VarSpec {
        name: name.to_string(),
        dofs: components.len().max(1),
        components,
        kind: FieldKind::Nodal,
        output: true,
    };
    let tensor = |name: &str| spec(name, comp[..n].iter().map(|c| format!("{name}_{c}")).collect());
    let triple = |name: &str| spec(name, (1..=3).map(|k| format!("{name} {k}")).collect());
    let mut out = Vec::new();
    if f.stress {
        out.push(tensor("stress"));
        out.push(spec("vonmises", Vec::new()));
        if f.principal {
            out.push(triple("principal stress"));
            out.push(spec("tresca", Vec::new()));
        }
    }
    if f.strain {
        out.push(tensor("strain"));
        if f.principal {
            out.push(triple("principal strain"));
        }
    }
    out
}

struct ElementData {
    law: ElasticLaw,
    density: Vec<f64>,
    damping: Option<Vec<f64>>,
    rayleigh_alpha: Option<f64>,
    rayleigh_beta: Option<f64>,
    load: Vec<[f64; 3]>,
    inertial: Option<Vec<[f64; 3]>>,
    pressure: Option<Vec<f64>>,
    stress_load: Option<Vec<Voigt6>>,
    strain_load: Option<Vec<Voigt6>>,
    gpa: Option<Vec<f64>>,
    temperature: Option<Vec<f64>>,
    heat_expansion: Option<Vec<[f64; 3]>>,
    ip: Option<IpProps>,
}

impl ElementData {
    /// The nodal data at `nodes` (the p basis's node order) instead of
    /// `mesh_nodes`.
    fn reordered(self, mesh_nodes: &[u32], nodes: &[u32]) -> Self {
        if mesh_nodes == nodes {
            return self;
        }
        let map: Vec<usize> =
            nodes.iter().map(|n| mesh_nodes.iter().position(|m| m == n).expect("same node set")).collect();
        fn take<T: Copy>(map: &[usize], v: &[T]) -> Vec<T> {
            map.iter().map(|&i| v[i]).collect()
        }
        let law = match self.law {
            ElasticLaw::Isotropic { young, poisson } => {
                ElasticLaw::Isotropic { young: take(&map, &young), poisson: take(&map, &poisson) }
            }
            ElasticLaw::Anisotropic { c, rotation } => ElasticLaw::Anisotropic { c: take(&map, &c), rotation },
        };
        Self {
            law,
            density: take(&map, &self.density),
            damping: self.damping.map(|v| take(&map, &v)),
            load: take(&map, &self.load),
            inertial: self.inertial.map(|v| take(&map, &v)),
            pressure: self.pressure.map(|v| take(&map, &v)),
            stress_load: self.stress_load.map(|v| take(&map, &v)),
            strain_load: self.strain_load.map(|v| take(&map, &v)),
            gpa: self.gpa.map(|v| take(&map, &v)),
            temperature: self.temperature.map(|v| take(&map, &v)),
            heat_expansion: self.heat_expansion.map(|v| take(&map, &v)),
            ..self
        }
    }
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

/// Centre of an element in its reference coordinates, as the boundary
/// orientation test wants it.
fn parent_center(basis: &dyn ElementBasis, xyz: &[[f64; 3]], etype: ElementType) -> Result<[f64; 3], SolverError> {
    let uvw = match etype.family() {
        3 | 7 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        5 => [0.25, 0.25, 0.25],
        6 => [0.0, 0.0, 1.0 / 3.0],
        _ => [0.0; 3],
    };
    let mut be = elmer_core::BasisEval::default();
    basis.eval(etype, xyz, uvw, &mut be).map_err(fatal)?;
    Ok(std::array::from_fn(|a| (0..xyz.len()).map(|k| be.basis[k] * xyz[k][a]).sum()))
}

/// The 6x6 `C` of one node from a `Youngs Modulus` block, as `InputTensor`
/// expands a vector or a matrix.
fn c_of_block(rows: usize, cols: usize, d: &[f64]) -> [[f64; 6]; 6] {
    if rows == 1 || cols == 1 {
        return diagonal_c(d);
    }
    let mut c = diagonal_c(&[]);
    for i in 0..rows.min(6) {
        for j in 0..cols.min(6) {
            c[i][j] = d[i * cols + j];
        }
    }
    c
}

/// Diagonal of a `Heat Expansion Coefficient` block.
fn beta_of_block(r: usize, c: usize, d: &[f64]) -> [f64; 3] {
    std::array::from_fn(|i| {
        if r == 1 && c == 1 {
            d[0]
        } else if r == 1 || c == 1 {
            d.get(i).copied().unwrap_or(0.0)
        } else if i < r && i < c {
            d[i * c + i]
        } else {
            0.0
        }
    })
}

/// `GetVoigtLoad`: a symmetric tensor in Voigt form, nodewise, either as one
/// array valued keyword or as `<name> 1` .. `<name> 6`.
fn voigt_load(ctx: &SolverCtx<'_>, section: &ParamMap, name: &str, e: ElemId) -> Result<Option<Vec<Voigt6>>, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    if let Some(a) = props::node_array(ctx, section, name, e)? {
        let k = a.rows.min(6);
        return Ok(Some(
            (0..n)
                .map(|node| {
                    let blk = a.node(node);
                    let mut v = [0.0; 6];
                    for (i, vi) in v.iter_mut().enumerate().take(k) {
                        *vi = blk[i * a.cols];
                    }
                    v
                })
                .collect(),
        ));
    }
    let mut out = vec![[0.0; 6]; n];
    let mut any = false;
    for i in 0..6 {
        if let Some(v) = props::scalar_nodes(ctx, section, &format!("{name} {}", i + 1), e)? {
            any = true;
            for (o, x) in out.iter_mut().zip(v) {
                o[i] = x;
            }
        }
    }
    Ok(any.then_some(out))
}

impl ElasticSolver {
    #[allow(clippy::too_many_arguments)]
    fn element_data(
        &self,
        ctx: &SolverCtx<'_>,
        e: ElemId,
        cdim: usize,
        model: Model,
        temp: &Option<FieldRef<'_>>,
    ) -> Result<ElementData, SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let bp = ctx.body_params(e).ok_or_else(|| fatal(format!("element {} has no body", e.0)))?;
        let mat = bp.material.ok_or_else(|| fatal(format!("body {} has no material", bp.body_id)))?;
        let umat = mat.get("umat subroutine").is_some();
        let key = "maxwell material";
        if mat.get(key).is_some() && (flag(Some(mat), key) || !matches!(mat.get(key), Some(Value::Logical(_)))) {
            return Err(unsupported("Maxwell material"));
        }
        for prefix in ["pre stress", "pre strain", "mesh velocity"] {
            if present_prefix(mat, prefix) {
                return Err(unsupported(&format!("material '{prefix}'")));
            }
        }
        if model.neo_hookean && present_prefix(mat, "heat expansion coefficient") {
            return Err(unsupported("'Heat Expansion Coefficient' with the neo-Hookean material"));
        }
        if let Some(bf) = bp.body_force {
            for prefix in [
                "stress volume source",
                "update reference displacement",
                "stress pressure im",
                "stress bodyforce 1 im",
                "stress bodyforce 2 im",
                "stress bodyforce 3 im",
            ] {
                if present_prefix(bf, prefix) {
                    return Err(unsupported(&format!("body force '{prefix}'")));
                }
            }
            if model.neo_hookean && (present_prefix(bf, "stress load") || present_prefix(bf, "strain load")) {
                return Err(unsupported("'Stress Load' / 'Strain Load' with the neo-Hookean material"));
            }
        }

        let young_ip = flag(Some(mat), "youngs modulus at ip");
        let poisson_ip = flag(Some(mat), "poisson ratio at ip");
        let beta_ip = flag(Some(mat), "heat expansion coefficient ip");
        let load_ip = flag(bp.body_force, "stress bodyforce at ip");
        let any_ip = young_ip || poisson_ip || beta_ip || load_ip;
        if any_ip && model.neo_hookean {
            return Err(unsupported("properties at integration points with the neo-Hookean material"));
        }
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
        let shape = if any_ip {
            let fe = self.fe();
            elmer_physics_elastic::shape_at_points(&*fe.basis, &*fe.quad, el.etype, &xyz, None).map_err(|x| fatal(x.to_string()))?
        } else {
            Vec::new()
        };
        let mut ip = IpProps::default();

        let scalar_poisson = |ip: &mut IpProps| -> Result<Vec<f64>, SolverError> {
            if poisson_ip {
                ip.poisson = props::scalar_ips(ctx, mat, "poisson ratio", e, &shape)?;
                Ok(vec![0.0; n])
            } else {
                // ElasticSolve.F90:2141: not read under "Incompressible".
                match props::scalar_nodes(ctx, mat, "poisson ratio", e)? {
                    Some(v) => Ok(v),
                    None if model.incompressible => Ok(vec![0.0; n]),
                    None => Err(fatal("material has no 'Poisson Ratio'")),
                }
            }
        };
        let law = if umat {
            // The UMAT is the material; its constants are read by the caller.
            ElasticLaw::Isotropic { young: vec![0.0; n], poisson: vec![0.0; n] }
        } else if young_ip {
            // Isotropy is assumed with it, as in ElasticSolve.
            ip.young = props::scalar_ips(ctx, mat, "youngs modulus", e, &shape)?;
            let poisson = scalar_poisson(&mut ip)?;
            ElasticLaw::Isotropic { young: vec![0.0; n], poisson }
        } else {
            let ym = props::node_array(ctx, mat, "youngs modulus", e)?.ok_or_else(|| fatal("material has no 'Youngs Modulus'"))?;
            if ym.rows == 1 && ym.cols == 1 {
                let poisson = scalar_poisson(&mut ip)?;
                ElasticLaw::Isotropic { young: ym.data, poisson }
            } else {
                if model.neo_hookean {
                    return Err(fatal("the neo-Hookean material needs a scalar 'Youngs Modulus'"));
                }
                let c: Vec<[[f64; 6]; 6]> = (0..n).map(|k| c_of_block(ym.rows, ym.cols, ym.node(k))).collect();
                let rotation = if flag(Some(mat), "rotate elasticity tensor") {
                    let mut t: T3 = [[0.0; 3]; 3];
                    for (i, row) in t.iter_mut().enumerate() {
                        let v = props::node_array(ctx, mat, &format!("material coordinates unit vector {}", i + 1), e)?
                            .ok_or_else(|| fatal("'Rotate Elasticity Tensor' without unit vectors"))?;
                        let v = v.node(0);
                        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                        if norm < f64::EPSILON {
                            return Err(fatal("material coordinate unit vector too short"));
                        }
                        for (a, x) in row.iter_mut().enumerate() {
                            *x = v[a] / norm;
                        }
                    }
                    Some(t)
                } else {
                    None
                };
                ElasticLaw::Anisotropic { c, rotation }
            }
        };
        if (young_ip || poisson_ip) && !matches!(law, ElasticLaw::Isotropic { .. }) {
            return Err(unsupported("properties at integration points need an isotropic material"));
        }

        let density = props::scalar_nodes(ctx, mat, "density", e)?.unwrap_or_else(|| vec![0.0; n]);
        let damping = props::scalar_nodes(ctx, mat, "damping", e)?;
        let rayleigh_alpha = props::const_real(ctx, mat, "rayleigh damping alpha")?;
        let rayleigh_beta = props::const_real(ctx, mat, "rayleigh damping beta")?;

        let mut load = vec![[0.0; 3]; n];
        let mut inertial = None;
        let mut pressure = None;
        let (mut stress_load, mut strain_load) = (None, None);
        if let Some(bf) = bp.body_force {
            if load_ip {
                let mut l = vec![[0.0; 3]; shape.len() / n];
                for i in 0..cdim {
                    if let Some(v) = props::scalar_ips(ctx, bf, &format!("stress bodyforce {}", i + 1), e, &shape)? {
                        for (li, x) in l.iter_mut().zip(v) {
                            li[i] = x;
                        }
                    }
                }
                ip.load = Some(l);
                ip.pressure =
                    Some(props::scalar_ips(ctx, bf, "stress pressure", e, &shape)?.unwrap_or_else(|| vec![0.0; shape.len() / n]));
            } else {
                for i in 0..cdim {
                    if let Some(v) = props::scalar_nodes(ctx, bf, &format!("stress bodyforce {}", i + 1), e)? {
                        for (l, x) in load.iter_mut().zip(v) {
                            l[i] = x;
                        }
                    }
                }
                pressure = props::scalar_nodes(ctx, bf, "stress pressure", e)?;
            }
            if present_prefix(bf, "inertial bodyforce") {
                let mut w = vec![[0.0; 3]; n];
                for i in 0..cdim {
                    if let Some(v) = props::scalar_nodes(ctx, bf, &format!("inertial bodyforce {}", i + 1), e)? {
                        for (l, x) in w.iter_mut().zip(v) {
                            l[i] = x;
                        }
                    }
                }
                inertial = Some(w);
            }
            if present_prefix(bf, "stress load") {
                stress_load = voigt_load(ctx, bf, "stress load", e)?;
            }
            if present_prefix(bf, "strain load") {
                strain_load = voigt_load(ctx, bf, "strain load", e)?;
            }
        }

        let reference = props::scalar_nodes(ctx, mat, "reference temperature", e)?.unwrap_or_else(|| vec![0.0; n]);
        let temperature = nodal_temperature(temp, el.nodes, &reference);
        let mut heat_expansion = None;
        if beta_ip {
            if temperature.is_some() {
                ip.beta = props::scalar_ips(ctx, mat, "heat expansion coefficient", e, &shape)?;
            }
        } else if let (Some(b), true) = (props::node_array(ctx, mat, "heat expansion coefficient", e)?, temperature.is_some()) {
            heat_expansion = Some((0..n).map(|k| beta_of_block(b.rows, b.cols, b.node(k))).collect());
        }
        let ip = (ip.young.is_some() || ip.poisson.is_some() || ip.beta.is_some() || ip.load.is_some()).then_some(ip);
        // "Gravitational Prestress Advection" with "GPA Coeff" (0 when absent).
        let gpa = match bp.body_force {
            Some(bf) if bf.logical("gravitational prestress advection") == Some(true) => {
                Some(props::scalar_nodes(ctx, bf, "gpa coeff", e)?.unwrap_or_else(|| vec![0.0; n]))
            }
            _ => None,
        };
        Ok(ElementData {
            law,
            density,
            damping,
            rayleigh_alpha,
            rayleigh_beta,
            load,
            inertial,
            pressure,
            stress_load,
            strain_load,
            gpa,
            temperature,
            heat_expansion,
            ip,
        })
    }

    /// One element of a UMAT material (ElasticSolve.F90:2040-2350): the
    /// material's routine and constants, `LocalMatrixWithUMAT` on its
    /// integration points' state, and at the first iterate its external
    /// loads into `BulkRHS`.
    #[allow(clippy::too_many_arguments)]
    fn umat_element(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        e: ElemId,
        cdim: usize,
        d: &ElementData,
        field: &FieldRef<'_>,
        temp: &Option<FieldRef<'_>>,
        scanning: bool,
        bulk_rhs: &mut [f64],
    ) -> Result<LocalSystem, SolverError> {
        let el = ctx.mesh.element(e);
        let n = el.nodes.len();
        let bp = ctx.body_params(e).ok_or_else(|| fatal(format!("element {} has no body", e.0)))?;
        if flag(bp.equation, "plane stress") {
            return Err(fatal("LocalMatrixWithUMAT: Cannot yet handle plane stress case"));
        }
        let mat = bp.material.ok_or_else(|| fatal(format!("body {} has no material", bp.body_id)))?;
        let (lib, entry) = elmer_solver::solver::procedure_value(mat, "umat subroutine")
            .ok_or_else(|| fatal("'UMAT Subroutine' is not a procedure"))?;
        let umat = crate::umat::lookup(&lib, &entry)
            .ok_or_else(|| unsupported(&format!("UMAT Subroutine \"{lib}\" \"{entry}\"")))?;
        let nprops = mat.integer("number of material constants").unwrap_or(0).max(0) as usize;
        let nstatev = mat.integer("number of state variables").unwrap_or(0).max(0) as usize;
        // `MaterialConstants(i,1)`.
        let consts: Vec<f64> = match mat.get("material constants") {
            Some(Value::Real { rows, cols, data }) if *cols == 1 || *rows == 1 => data.clone(),
            Some(Value::Real { rows, cols, data }) => (0..*rows).map(|i| data[i * cols]).collect(),
            _ => Vec::new(),
        };
        if consts.len() < nprops {
            return Err(fatal("Check the size of Material Constants array"));
        }
        let cmname = mat.string("name").unwrap_or("").to_string();
        // `LocalTemperature`: the reference temperature, the temperature
        // where the field has it.
        let mut temperature = props::scalar_nodes(ctx, mat, "reference temperature", e)?.unwrap_or_else(|| vec![0.0; n]);
        if let Some(t) = temp.as_ref().and_then(|t| t.perm.map(|p| (p, t.values))) {
            for (k, &node) in el.nodes.iter().enumerate() {
                if let Some(&s) = t.0.get(node as usize).filter(|&&s| s != elmer_core::NO_DOF) {
                    temperature[k] = t.1[s as usize];
                }
            }
        }
        let disp = local_displacement(field, el.nodes, cdim);
        let saved = self.umat.as_ref().expect("allocated before the iteration");
        let prev_field = FieldRef { values: &saved.prev, ..*field };
        let prev = if scanning { local_displacement(&prev_field, el.nodes, cdim) } else { vec![[0.0; 3]; n] };
        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let fe = self.fe();
        let np = rel_np(&*fe.quad, el.etype, RelOrders::of(ctx.params()).bulk)?;
        let points = fe.quad.gauss_points(el.etype, np);
        let ntens = if cdim == 2 { 4 } else { 6 };
        // `GetIpIndex(t, Element, UmatEnergyVar)`: the element's points.
        let (lo, hi) = {
            let v = ctx.sim.vars.get("umatenergy").ok_or_else(|| fatal("Could not find variable \"UmatEnergy\""))?;
            let p = v.perm.as_deref().ok_or_else(|| fatal("UmatEnergy has no permutation"))?;
            (p[e.0 as usize] as usize, p[e.0 as usize + 1] as usize)
        };
        if hi - lo != points.len() {
            return Err(fatal("UMAT: the integration points of the IP variables differ from the element's"));
        }
        let max_statev = ctx.sim.vars.get("umatstate").map_or(0, |v| v.dofs);
        if nstatev > max_statev {
            return Err(fatal("UMAT: more state variables than UmatState holds"));
        }
        let npt = points.len();
        let mut stress = vec![0.0; ntens * npt];
        let mut energy = vec![0.0; 3 * npt];
        let mut state = vec![0.0; max_statev * npt];
        let st0 = (
            &saved.stress0[ntens * lo..ntens * hi],
            &saved.energy0[3 * lo..3 * hi],
            if max_statev > 0 { &saved.state0[max_statev * lo..max_statev * hi] } else { &[][..] },
        );
        let time = ctx.sim.vars.scalar("time").unwrap_or(0.0);
        let dt = ctx.sim.vars.scalar("timestep size").unwrap_or(0.0);
        let inp = crate::umat::UmatInput {
            etype: el.etype,
            xyz: &xyz,
            cdim,
            disp: &disp,
            prev: &prev,
            load: &d.load,
            inertial: d.inertial.as_deref(),
            density: &d.density,
            temperature: &temperature,
            props: &consts[..nprops.min(consts.len())],
            nstatev,
            cmname: &cmname,
            umat,
            time,
            dt,
            element_index: e.0 as usize + 1,
        };
        let ndof = n * cdim;
        let mut local = LocalSystem::zeros(ndof);
        let mut ext = vec![0.0; ndof];
        let mut ip = crate::umat::IpState {
            stress0: st0.0,
            energy0: st0.1,
            state0: st0.2,
            stress: &mut stress,
            energy: &mut energy,
            state: &mut state,
            max_statev,
        };
        crate::umat::umat_local(&*fe.basis, &points, &inp, &mut ip, &mut local.stiff, &mut local.force, &mut ext)
            .map_err(fatal)?;
        for (name, vals, w) in [("umatstress", &stress, ntens), ("umatenergy", &energy, 3), ("umatstate", &state, max_statev)] {
            if w == 0 {
                continue;
            }
            if let Some(v) = ctx.sim.vars.get_mut(name) {
                v.values[w * lo..w * hi].copy_from_slice(vals);
            }
        }
        if !bulk_rhs.is_empty() {
            let perm = field.perm.ok_or_else(|| fatal("displacement has no permutation"))?;
            for (k, &node) in el.nodes.iter().enumerate() {
                let s = perm[node as usize];
                if s == elmer_core::NO_DOF {
                    continue;
                }
                for i in 0..cdim {
                    bulk_rhs[s as usize * field.dofs + i] += ext[cdim * k + i];
                }
            }
        }
        Ok(local)
    }

    /// Temperature snapshot of the primary-variable-independent field, taken
    /// before the mutable assembly.
    fn temperature_field(ctx: &SolverCtx<'_>) -> Option<(Vec<u32>, Vec<f64>)> {
        ctx.variable("temperature").and_then(|f| f.perm.map(|p| (p.to_vec(), f.values.to_vec())))
    }

    /// `ComputeStressAndStrain` for the small-strain material: L2 projection of
    /// the elastic strain and stress onto the displacement's nodes, then von
    /// Mises and the principal values.
    fn postprocess(
        &self,
        ctx: &mut SolverCtx<'_>,
        flags: PostFlags,
        cdim: usize,
        axi: bool,
        model: Model,
        temp_ref: &Option<FieldRef<'_>>,
    ) -> Result<(), SolverError> {
        let coord = if axi { CoordSystem::AxiSymmetric } else { CoordSystem::Cartesian };
        let ncomp = post::output_components(cdim);
        let (disp_perm, disp_vals, dofs) = {
            let f = ctx.primary().ok_or_else(|| fatal("solver has no displacement variable"))?;
            (
                f.perm.ok_or_else(|| fatal("displacement has no permutation"))?.to_vec(),
                f.values.to_vec(),
                f.dofs,
            )
        };
        let nslots = disp_vals.len() / dofs;
        let mut mass = post::MassSystem::new(nslots);
        let mut fstrain = vec![0.0; nslots * ncomp];
        let mut fstress = vec![0.0; nslots * ncomp];

        let elems: Vec<ElemId> = ctx.active_elements().collect();
        let fe = self.fe();
        for e in elems {
            let body = ctx.body_params(e);
            if flags.masked && !flag(body.as_ref().and_then(|b| b.equation), "calculate stresses") {
                continue;
            }
            let plane_stress = flag(body.as_ref().and_then(|b| b.equation), "plane stress");
            let el = ctx.mesh.element(e);
            // A p element's basis takes its nodes in the order of its dofs (p
            // tetrahedra are reordered); the nodal data follows them.
            let nodes: Vec<u32> = ctx.solver_dofs().map_or_else(|| el.nodes.to_vec(), |dd| dd.element_nodes(e));
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
            let d = self.element_data(ctx, e, cdim, model, temp_ref)?.reordered(el.nodes, &nodes);
            // The dof slots: the nodes' through the permutation, or the p
            // element's (`GetElementDOFs`), with its scalar basis and rule.
            let pe = match ctx.solver_dofs() {
                Some(dd) => {
                    let mut s = Vec::new();
                    dd.element_slots(e, &mut s);
                    let slots: Vec<Option<usize>> =
                        s.into_iter().map(|s| (s != elmer_core::NO_DOF).then_some(s as usize)).collect();
                    Some((
                        slots,
                        dd.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
                        dd.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?,
                    ))
                }
                None => None,
            };
            let slots: Vec<Option<usize>> = match &pe {
                Some(p) => p.0.clone(),
                None => el
                    .nodes
                    .iter()
                    .map(|&n| disp_perm.get(n as usize).copied().filter(|&s| s != elmer_core::NO_DOF).map(|s| s as usize))
                    .collect(),
            };
            let cd = if model.incompressible { cdim } else { dofs.min(3) };
            let disp: Vec<[f64; 3]> = slots
                .iter()
                .map(|s| match s {
                    Some(s) => std::array::from_fn(|a| if a < cd { disp_vals[s * dofs + a] } else { 0.0 }),
                    None => [0.0; 3],
                })
                .collect();
            // `Pressure = SUM(LocalDisplacement(DOFs,1:n) * Basis(1:n))`.
            let pressure: Option<Vec<f64>> = model.incompressible.then(|| {
                slots[..el.nodes.len()].iter().map(|s| s.map_or(0.0, |s| disp_vals[s * dofs + cdim])).collect()
            });
            if pe.is_some() && d.ip.is_some() {
                return Err(unsupported("material values at integration points in p-element stress postprocessing"));
            }
            let inp = StressPostInput {
                etype: el.etype,
                xyz: &xyz,
                cdim,
                coord,
                plane_stress,
                law: &d.law,
                disp: &disp,
                temperature: d.temperature.as_deref(),
                heat_expansion: d.heat_expansion.as_deref(),
                ip: d.ip.as_ref(),
                want_strain: flags.strain,
                want_stress: flags.stress,
                np: rel_np(&*fe.quad, el.etype, RelOrders::of(ctx.params()).bulk)?,
                large_deflection: model.large,
                neo_hookean: model.neo_hookean,
                pressure: pressure.as_deref(),
            };
            let umat_here = ctx.body_params(e).and_then(|b| b.material).is_some_and(|m| m.get("umat subroutine").is_some());
            let loc = if umat_here {
                // GenerateStressVariable / GenerateStrainVariable.
                let ntens = if cdim == 2 { 4 } else { 6 };
                let v = ctx.sim.vars.get("umatstress").ok_or_else(|| fatal("Could not find variable \"UmatStress\""))?;
                let p = v.perm.as_deref().ok_or_else(|| fatal("UmatStress has no permutation"))?;
                let (lo, hi) = (p[e.0 as usize] as usize, p[e.0 as usize + 1] as usize);
                let points = fe.quad.gauss_points(el.etype, None);
                crate::umat::umat_projection_local(
                    &*fe.basis,
                    &points,
                    el.etype,
                    &xyz,
                    &disp,
                    &v.values[ntens * lo..ntens * hi],
                    ntens,
                    ncomp,
                    model.large,
                )
                .map_err(fatal)?
            } else {
                match &pe {
                    Some((sl, basis, rule)) => {
                        if inp.np.is_some() {
                            return Err(unsupported("a relative integration order on p elements"));
                        }
                        post::stress_strain_local_p(&**basis, rule, sl.len(), &inp)
                    }
                    None => stress_strain_local(&*fe.basis, &*fe.quad, &inp),
                }
                .map_err(|x| fatal(x.to_string()))?
            };
            mass.add(&slots, &loc.mass);
            for (p, s) in slots.iter().enumerate() {
                let Some(s) = s else { continue };
                for k in 0..ncomp {
                    if flags.strain {
                        fstrain[ncomp * s + k] += loc.strain[ncomp * p + k];
                    }
                    if flags.stress {
                        fstress[ncomp * s + k] += loc.stress[ncomp * p + k];
                    }
                }
            }
        }

        let mut projector = post::Projector::new(ctx, &mass)?;
        let comm = ctx.comm;
        let mut project = |f: &[f64]| -> Result<Vec<f64>, SolverError> {
            let mut out = vec![0.0; nslots * ncomp];
            let mut b = vec![0.0; nslots];
            let mut x = vec![0.0; nslots];
            for k in 0..ncomp {
                for s in 0..nslots {
                    b[s] = f[ncomp * s + k];
                }
                projector.solve(comm, &b, &mut x)?;
                for s in 0..nslots {
                    out[ncomp * s + k] = x[s];
                }
            }
            Ok(out)
        };
        let set = |ctx: &mut SolverCtx<'_>, name: &str, v: Vec<f64>| -> Result<(), SolverError> {
            let var = ctx.sim.vars.get_mut(name).ok_or_else(|| fatal(format!("variable '{name}' does not exist")))?;
            if var.values.len() != v.len() {
                return Err(fatal(format!("variable '{name}' has the wrong size")));
            }
            var.values = v;
            var.values_changed = true;
            Ok(())
        };

        if flags.stress {
            ctx.info(CALLER, "Calculating stress components", 7);
            let stress = project(&fstress)?;
            let mut vm = vec![0.0; nslots];
            let (mut pri, mut tre) = (vec![0.0; nslots * 3], vec![0.0; nslots]);
            // ElasticSolve.F90:989: with a UMAT no von Mises or principal
            // values are computed.
            let umat_on = ctx.sim.model.materials.iter().any(|m| m.get("umat subroutine").is_some());
            for s in 0..nslots {
                if umat_on {
                    break;
                }
                let t = post::vector_to_tensor(&stress[ncomp * s..ncomp * (s + 1)], ncomp);
                vm[s] = post::von_mises(&t);
                if flags.principal {
                    let w = post::principal_values(&t);
                    pri[3 * s..3 * s + 3].copy_from_slice(&w);
                    tre[s] = post::tresca(&w);
                }
            }
            set(ctx, "stress", stress)?;
            set(ctx, "vonmises", vm)?;
            if flags.principal {
                set(ctx, "principal stress", pri)?;
                set(ctx, "tresca", tre)?;
            }
        }
        if flags.strain {
            ctx.info(CALLER, "Calculating strain components", 7);
            let strain = project(&fstrain)?;
            if flags.principal {
                let mut pri = vec![0.0; nslots * 3];
                for s in 0..nslots {
                    let t = post::vector_to_tensor(&strain[ncomp * s..ncomp * (s + 1)], ncomp);
                    pri[3 * s..3 * s + 3].copy_from_slice(&post::principal_values(&t));
                }
                set(ctx, "principal strain", pri)?;
            }
            set(ctx, "strain", strain)?;
        }
        Ok(())
    }
}

/// The displacement at the nodes of an element, `cdim` components.
fn local_displacement(f: &FieldRef<'_>, nodes: &[u32], cdim: usize) -> Vec<[f64; 3]> {
    nodes
        .iter()
        .map(|&n| match f.perm.and_then(|p| p.get(n as usize)) {
            Some(&s) if s != elmer_core::NO_DOF => {
                std::array::from_fn(|a| if a < cdim.min(f.dofs) { f.values[s as usize * f.dofs + a] } else { 0.0 })
            }
            _ => [0.0; 3],
        })
        .collect()
}

impl Solver for ElasticSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        check_solver_params(ctx.params())?;
        let dim = ctx.mesh.dim();
        let dofs = ctx.primary().map(|f| f.dofs);
        let model = model_of(ctx.params());
        // ElasticSolve.F90:728-750: the checks of "Incompressible".
        if model.incompressible {
            if model.neo_hookean {
                return Err(fatal("\"Incompressible\" is the linear mixed formulation and cannot be combined with a neo-Hookean material"));
            }
            if dofs != Some(dim + 1) {
                return Err(fatal("\"Incompressible\" needs one more variable DOF than there are dimensions"));
            }
            if model.large {
                return Err(fatal("\"Incompressible\" is a small strain formulation: set \"Large Deflection = False\""));
            }
            // ElasticSolver_Init: the measured integration rules of the MINI
            // element (`ListAddNewString`).
            let bubbles = ctx.params().string("element").is_some_and(|s| s.contains("b:"));
            if bubbles && ctx.params().get("element integration points").is_none() {
                let id = ctx.solver_id();
                let ix = ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot");
                ctx.sim.slots[ix]
                    .params
                    .insert("Element Integration Points", Value::Str("-tri 7 -tetra 24 -prism 85".into()));
            }
        } else if dofs != Some(dim) {
            return Err(fatal(format!("'Variable DOFs' must equal the coordinate dimension {dim}, got {dofs:?}")));
        }
        // ElasticSolver_Init: the time derivative order is two unless the sif says otherwise.
        if ctx.params().get("time derivative order").is_none() {
            ctx.set_time_derivative_order(2)?;
        }
        // The postprocessing fields share the displacement's permutation.
        let flags = post_flags(ctx);
        if flags.any() {
            let name = ctx.variable_name().unwrap_or_default().to_string();
            let (perm, nslots) = {
                let v = ctx.sim.vars.get(&name).ok_or_else(|| fatal("solver has no displacement variable"))?;
                (v.perm.clone(), v.values.len() / v.dofs)
            };
            for spec in field_specs(dim, flags) {
                if ctx.sim.vars.id(&spec.name).is_some() {
                    continue;
                }
                let mut var = Variable::new(&spec.name, FieldKind::Nodal, spec.dofs, perm.clone(), vec![0.0; nslots * spec.dofs]);
                var.output = spec.output;
                ctx.sim.vars.add_with_components(var, &spec.components)?;
            }
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        check_solver_params(ctx.params())?;
        let model = model_of(ctx.params());
        let second_order = transient && ctx.params().integer("time derivative order").unwrap_or(2) == 2;
        let quasi = flag(Some(ctx.params()), "quasi stationary");
        if quasi && !transient {
            ctx.warn(CALLER, "\"Quasi Stationary\" only has meaning in a transient simulation");
        }
        let cdim = ctx.mesh.dim();
        let coord_name = elmer_core::normalize_key(ctx.simulation_params().string("coordinate system").unwrap_or("cartesian"));
        let axi = match coord_name.as_str() {
            "axi symmetric" | "cylindric symmetric" => true,
            s if s.starts_with("cartesian") => false,
            other => return Err(unsupported(&format!("coordinate system '{other}'"))),
        };
        let coord = if axi { CoordSystem::AxiSymmetric } else { CoordSystem::Cartesian };
        if !self.mesh_note && ctx.params().logical("displace mesh") != Some(false) {
            self.mesh_note = true;
            ctx.warn(CALLER, "'Displace Mesh' is not applied: the mesh stays undisplaced for the solvers that follow");
        }

        let mut max_iter = ctx.params().integer("nonlinear system max iterations").unwrap_or(1).max(1);
        let mut min_iter = ctx.params().integer("nonlinear system min iterations").unwrap_or(0);

        // Model lumping (ElasticSolve.F90:1160): six load cases on the
        // nonlinear counter, sharing one bulk system.
        let lumping = flag(Some(ctx.params()), "model lumping");
        if lumping {
            if cdim != 3 {
                return Err(fatal("Model lumping is implemented in 3D only"));
            }
            max_iter = 6;
            min_iter = 6;
            let id = ctx.solver_id();
            let ix = ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot");
            ctx.sim.slots[ix].params.insert("Constant Bulk System", Value::Logical(true));
        }
        // "Constant System", "Constant Bulk System", "Constant Bulk Matrix" (:1541).
        let cs = flag(Some(ctx.params()), "constant system");
        let cbs = flag(Some(ctx.params()), "constant bulk system");
        let cbm = flag(Some(ctx.params()), "constant bulk matrix");
        if cs || cbs || cbm {
            if transient {
                return Err(fatal("The \"Constant ...\" keywords are applicable to steady cases only here"));
            }
            if model.neo_hookean || model.large {
                return Err(fatal(
                    "The \"Constant ...\" keywords need a linear material: set \"Large Deflection = False\" and use neither a UMAT nor a neo-Hookean material",
                ));
            }
            if cs && (cbs || cbm) {
                ctx.warn(CALLER, "\"Constant System\" is superseded by the narrower \"Constant Bulk System\"/\"Constant Bulk Matrix\" given beside it");
            }
        }
        let mut lump = if lumping {
            let fe = self.fe();
            Some(elmer_physics_elastic::lumping::Lump::init(ctx, &*fe.basis, &*fe.quad)?)
        } else {
            None
        };
        let global_pseudo = flag(Some(ctx.params()), "pseudo-traction");

        // GEOMETRIC STIFFNESS and STABILITY ANALYSIS (ElasticSolve.F90:1169):
        // two passes, a static solve with the eigen analysis switched off and
        // the eigen analysis with the geometric stiffness of that solve's
        // displacement, added to the stiffness or (stability) subtracted in
        // the mass slot in place of the density mass.
        let stability = flag(Some(ctx.params()), "stability analysis");
        let geometric_stiffness = flag(Some(ctx.params()), "geometric stiffness");
        let geometric = stability || geometric_stiffness;
        if geometric {
            if stability && geometric_stiffness {
                return Err(fatal(
                    "\"Stability Analysis\" and \"Geometric Stiffness\" are the same term in different matrices and cannot both be active",
                ));
            }
            if axi {
                return Err(fatal("The geometric stiffness is implemented in cartesian coordinates only, as in StressSolve"));
            }
            if model.large {
                return Err(fatal(
                    "The geometric stiffness is already part of the tangent under \"Large Deflection\", where it comes from the current iterate rather than from a prestressed state. Set \"Large Deflection = False\" to ask for the prestressed eigen analysis these keywords mean",
                ));
            }
            if model.neo_hookean {
                return Err(fatal("The geometric stiffness is implemented by LocalMatrix, so not with UMAT or the neo-Hookean model"));
            }
            if lumping || cs || cbs || cbm || transient {
                return Err(unsupported("geometric stiffness with Model Lumping, the Constant keywords or a transient run"));
            }
            max_iter = 2;
            min_iter = 2;
        }
        let orig_eigen = ctx.params().logical("eigen analysis") == Some(true);
        let eigen_now = |ctx: &SolverCtx<'_>| ctx.nof_eigen() > 0 && ctx.params().logical("eigen analysis") == Some(true);
        if eigen_now(ctx) && (transient || lumping || cs || cbs || cbm || post_flags(ctx).any()) {
            return Err(unsupported("eigen analysis in a transient run, with Model Lumping, the Constant keywords or stress postprocessing"));
        }

        // Snapshot the temperature field (read-only) before the mutable assembly.
        let temp_owned = Self::temperature_field(ctx);
        let temp_ref = temp_owned.as_ref().map(|(p, v)| FieldRef {
            name: "temperature",
            kind: elmer_core::FieldKind::Nodal,
            dofs: 1,
            perm: Some(p),
            values: v,
        });

        // ElasticSolve.F90:697-1090, :1565: a UMAT material solves for the
        // displacement increment from the residual of the last converged
        // state.
        let use_umat = ctx.sim.model.materials.iter().any(|m| m.get("umat subroutine").is_some());
        let scanning = ctx
            .simulation_params()
            .string("simulation type")
            .is_some_and(|t| elmer_core::normalize_key(t) == "scanning");
        let mut nonlin_res0 = 1.0;
        if use_umat {
            if transient {
                return Err(fatal("UMAT version does not yet support transient simulation"));
            }
            for (bad, what) in [
                (model.large, "Large Deflection"),
                (axi, "axial symmetry"),
                (flag(Some(ctx.params()), "use hencky strain"), "Hencky strain"),
                (flag(Some(ctx.params()), "initialize state variables"), "Initialize State Variables"),
                (flag(Some(ctx.params()), "linear system symmetric"), "Linear System Symmetric"),
                (ctx.solver_dofs().is_some(), "p elements"),
                (lumping || cs || cbs || cbm || geometric, "Model Lumping, a constant system or geometric stiffness"),
            ] {
                if bad {
                    return Err(unsupported(&format!("UMAT materials with {what}")));
                }
            }
            ctx.params_mut().insert("Skip Compute Nonlinear Change", Value::Logical(true));
            if self.umat.is_none() {
                let mut saved = UmatSaved::default();
                for (name, keep) in [("umatstress", 0), ("umatenergy", 1), ("umatstate", 2)] {
                    let Some(v) = ctx.sim.vars.get_mut(name) else {
                        if keep < 2 {
                            return Err(fatal(format!("Could not find variable \"{name}\"")));
                        }
                        continue;
                    };
                    v.values.iter_mut().for_each(|x| *x = 0.0);
                    let zeros = vec![0.0; v.values.len()];
                    match keep {
                        0 => saved.stress0 = zeros,
                        1 => saved.energy0 = zeros,
                        _ => saved.state0 = zeros,
                    }
                }
                saved.prev = vec![0.0; ctx.primary().map_or(0, |f| f.values.len())];
                self.umat = Some(saved);
            }
        }

        for iter in 1..=max_iter {
            // The keyword is rewritten: the core reads it when the system is solved.
            let geometric_active = geometric && iter > 1;
            if geometric {
                ctx.params_mut().insert("Eigen Analysis", Value::Logical(geometric_active && orig_eigen));
            }
            let eigen = eigen_now(ctx);
            ctx.info(CALLER, &format!("ELASTICITY ITERATION {iter}"), 4);
            // `TotalSol(:) = Displacement(:)`; `BulkRHS`, the external loads
            // of the first iterate.
            let total: Vec<f64> = if use_umat { ctx.primary().map(|f| f.values.to_vec()).unwrap_or_default() } else { Vec::new() };
            let mut bulk_rhs: Vec<f64> = if use_umat && iter == 1 { vec![0.0; total.len()] } else { Vec::new() };
            ctx.default_initialize()?;
            // Reuse what an earlier pass saved (:1636-1680): the matrix, and
            // the bulk load unless only the matrix is constant.
            let cbm_in_use = cbm && ctx.saved_bulk().is_some();
            let (mut skip_bulk, mut skip_all) = (false, false);
            if (cs || cbs || cbm) && ctx.saved_bulk().is_some() {
                ctx.restore_saved_bulk(cbs || cs);
                skip_bulk = cbs;
                skip_all = cs;
            }
            let disp_field = {
                let f = ctx.primary().ok_or_else(|| fatal("solver has no displacement variable"))?;
                (f.perm.map(<[u32]>::to_vec), f.values.to_vec(), f.dofs)
            };
            let field = FieldRef {
                name: "displacement",
                kind: elmer_core::FieldKind::Nodal,
                dofs: disp_field.2,
                perm: disp_field.0.as_deref(),
                values: &disp_field.1,
            };

            let elems: Vec<ElemId> = if skip_all || skip_bulk { Vec::new() } else { ctx.active_elements().collect() };
            for e in elems {
                let el = ctx.mesh.element(e);
                // A p element's basis takes its nodes in the order of its dofs (p
                // tetrahedra are reordered); the nodal data follows them.
                let nodes: Vec<u32> = ctx.solver_dofs().map_or_else(|| el.nodes.to_vec(), |dd| dd.element_nodes(e));
                let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                let d = self.element_data(ctx, e, cdim, model, &temp_ref)?.reordered(el.nodes, &nodes);
                if use_umat {
                    let local = self.umat_element(ctx, e, cdim, &d, &field, &temp_ref, scanning, &mut bulk_rhs)?;
                    ctx.default_update(&local, e)?;
                    continue;
                }
                // ElasticSolve.F90:2080: "Incompressible" overrides "Plane Stress".
                let plane_stress = flag(ctx.body_params(e).and_then(|b| b.equation), "plane stress") && !model.incompressible;
                // A p element: its dofs (`GetElementDOFs`, bubbles in the
                // global system), the scalar `ElementInfo` and its rule.
                let pe = match ctx.solver_dofs() {
                    Some(dd) => {
                        let mut slots = Vec::new();
                        dd.element_slots(e, &mut slots);
                        Some((
                            slots,
                            dd.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
                            dd.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?,
                        ))
                    }
                    None => None,
                };
                let disp = match &pe {
                    _ if model.linear_model => vec![[0.0; 3]; pe.as_ref().map_or(el.nodes.len(), |p| p.0.len())],
                    Some((slots, _, _)) => slots
                        .iter()
                        .map(|&s| {
                            if s == elmer_core::NO_DOF {
                                [0.0; 3]
                            } else {
                                std::array::from_fn(|a| if a < cdim { field.values[s as usize * field.dofs + a] } else { 0.0 })
                            }
                        })
                        .collect(),
                    None => local_displacement(&field, el.nodes, cdim),
                };
                if pe.is_some() && (model.neo_hookean || d.ip.is_some()) {
                    return Err(unsupported("p elements with the neo-Hookean material or properties at integration points"));
                }
                let fe = self.fe();
                let mut local = LocalSystem::zeros(0);
                if model.neo_hookean {
                    let ElasticLaw::Isotropic { young, poisson } = &d.law else {
                        return Err(fatal("the neo-Hookean material needs an isotropic law"));
                    };
                    let inp = NeoHookeInput {
                        etype: el.etype,
                        xyz: &xyz,
                        cdim,
                        coord,
                        plane_stress,
                        young,
                        poisson,
                        disp: &disp,
                        density: &d.density,
                        damping: d.damping.as_deref(),
                        load: &d.load,
                        inertial: d.inertial.as_deref(),
                        want_mass: transient || eigen,
                        np: rel_np(&*fe.quad, el.etype, RelOrders::of(ctx.params()).bulk)?,
                    };
                    neo_hookean_local(&*fe.basis, &*fe.quad, &inp, &mut local).map_err(|x| fatal(x.to_string()))?;
                } else {
                    // `ElasticGaussPoints`: an `Element Integration Points` count
                    // for the family first, else the relative order.
                    let ngp = match ctx.params().string("element integration points") {
                        Some(def) => elmer_solver::parse_elemental_gauss_rules(def)?[(el.etype.0 / 100) as usize],
                        None => 0,
                    };
                    let bulk_np = if ngp > 0 {
                        Some(ngp)
                    } else {
                        rel_np(&*fe.quad, el.etype, RelOrders::of(ctx.params()).bulk)?
                    };
                    if bulk_np.is_some() && d.ip.is_some() {
                        return Err(unsupported("properties at integration points with a relative integration order"));
                    }
                    let inp = BulkInput {
                        etype: el.etype,
                        xyz: &xyz,
                        cdim,
                        coord,
                        plane_stress,
                        large_deflection: model.large,
                        law: &d.law,
                        disp: &disp,
                        density: &d.density,
                        damping: d.damping.as_deref(),
                        load: &d.load,
                        inertial: d.inertial.as_deref(),
                        pressure_load: d.pressure.as_deref(),
                        stress_load: d.stress_load.as_deref(),
                        strain_load: d.strain_load.as_deref(),
                        gpa: d.gpa.as_deref(),
                        temperature: d.temperature.as_deref(),
                        heat_expansion: d.heat_expansion.as_deref(),
                        want_mass: transient || (eigen && !stability),
                        np: bulk_np,
                        ip: d.ip.as_ref(),
                        geometric: geometric_active.then_some(stability),
                        incompressible: model.incompressible,
                        pstab: flag(Some(ctx.params()), "pressure stabilization").then(|| {
                            (
                                ctx.params().real("pressure stabilization coefficient").unwrap_or(1.0),
                                crate::diameter::element_diameter(el.etype, &xyz),
                            )
                        }),
                    };
                    match &pe {
                        Some((slots, basis, rule)) => {
                            if ngp > 0 {
                                let dd = ctx.solver_dofs().expect("p element");
                                let r = dd.element_rule_np(e, ngp).map_err(|m| SolverError::fatal("GaussPoints", m))?;
                                elastic_bulk_local_p(&**basis, &r, slots.len(), &inp, &mut local)
                            } else {
                                if bulk_np.is_some() {
                                    return Err(unsupported("a relative integration order on p elements"));
                                }
                                elastic_bulk_local_p(&**basis, rule, slots.len(), &inp, &mut local)
                            }
                        }
                        None => elastic_bulk_local(&*fe.basis, &*fe.quad, &inp, &mut local),
                    }
                    .map_err(|x| fatal(x.to_string()))?;
                }
                if transient {
                    let nn = local.n * local.n;
                    let mut mass = local.mass.take().unwrap_or_else(|| vec![0.0; nn]);
                    // "Quasi Stationary" (ElasticSolve.F90:4743): no inertial term
                    // unless the element carries a damping.
                    let damped = d.damping.as_ref().is_some_and(|v| v.iter().any(|x| *x != 0.0))
                        || d.rayleigh_alpha.is_some()
                        || d.rayleigh_beta.is_some();
                    if quasi && !damped {
                        mass.iter_mut().for_each(|v| *v = 0.0);
                    }
                    let mut damp = local.damp.take().unwrap_or_else(|| vec![0.0; nn]);
                    if let Some(a) = d.rayleigh_alpha {
                        damp.iter_mut().zip(&mass).for_each(|(c, m)| *c += a * m);
                    }
                    if let Some(b) = d.rayleigh_beta {
                        damp.iter_mut().zip(&local.stiff).for_each(|(c, k)| *c += b * k);
                    }
                    if second_order {
                        local.mass = Some(mass);
                        local.damp = Some(damp);
                        ctx.default_second_order_time(&mut local, e, dt)?;
                    } else {
                        // First order in time: the damping matrix carries d/dt
                        // (`Default1stOrderTime( DAMP, STIFF, FORCE )`).
                        local.mass = Some(damp);
                        ctx.default_first_order_time(&mut local, e, dt)?;
                    }
                }
                if eigen {
                    // DefaultUpdateMass: the density mass, or minus the
                    // geometric stiffness under stability analysis.
                    let mass = local.mass.take().unwrap_or_else(|| vec![0.0; local.n * local.n]);
                    ctx.default_update_mass(&mass, e)?;
                }
                if cbm_in_use {
                    // DefaultUpdateForce: the matrix was restored wholesale.
                    local.stiff.iter_mut().for_each(|v| *v = 0.0);
                }
                ctx.default_update(&local, e)?;
            }
            if skip_all || skip_bulk {
                ctx.reuse_bulk_assembly(skip_all);
            } else {
                ctx.default_finish_bulk_assembly()?;
            }

            let belems: Vec<ElemId> = if skip_all { Vec::new() } else { ctx.boundary_elements().collect() };
            for e in belems {
                let el = ctx.mesh.element(e);
                let n = el.nodes.len();
                let bc = ctx.bc_params(e).expect("boundary_elements filters on BC").clone();
                if let Some(true) = bc.logical("fsi bc") {
                    return Err(unsupported("FSI boundary condition"));
                }
                for (prefix, what) in [
                    ("stress load", "boundary 'Stress Load'"),
                    ("force 1 im", "imaginary boundary load"),
                    ("force 2 im", "imaginary boundary load"),
                    ("force 3 im", "imaginary boundary load"),
                    ("normal force im", "imaginary boundary load"),
                ] {
                    if present_prefix(&bc, prefix) {
                        return Err(unsupported(what));
                    }
                }
                let nodal = |key: &str| props::scalar_nodes(ctx, &bc, key, e);
                let mut traction = vec![[0.0; 3]; n];
                let mut got = false;
                for prefix in ["surface traction", "force"] {
                    if present_prefix(&bc, prefix) {
                        for i in 0..cdim {
                            if let Some(v) = nodal(&format!("{prefix} {}", i + 1))? {
                                for (t, x) in traction.iter_mut().zip(v) {
                                    t[i] = x;
                                }
                            }
                        }
                        got = true;
                        break;
                    }
                }
                // A lumping load case given as a pure force or moment replaces
                // the surface load (ModelLumpingLoads).
                if let Some(l) = lump.as_ref().filter(|l| !l.fix_displacement) {
                    if bc.logical("model lumping boundary") == Some(true) {
                        let exyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
                        traction = l.loads(iter as usize, &exyz);
                        got = true;
                    }
                }
                let mut beta = nodal("normal surface traction")?;
                if beta.is_none() {
                    beta = nodal("normal force")?;
                }
                got |= beta.is_some();

                let mut spring_normal = None;
                let mut spring_tensor = None;
                if present_prefix(&bc, "spring") {
                    if let Some(k) = nodal("spring")? {
                        spring_normal = Some(k);
                    } else {
                        let mut t = vec![[[0.0; 3]; 3]; n];
                        for i in 0..cdim {
                            if let Some(v) = nodal(&format!("spring {}", i + 1))? {
                                for (tk, x) in t.iter_mut().zip(v) {
                                    tk[i][i] = x;
                                }
                            }
                            for j in 0..cdim {
                                if let Some(v) = nodal(&format!("spring {}{}", i + 1, j + 1))? {
                                    for (tk, x) in t.iter_mut().zip(v) {
                                        tk[i][j] = x;
                                    }
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

                // The left parent, unless some of its nodes are outside the system.
                let in_system = |p: ElemId| {
                    ctx.mesh
                        .element(p)
                        .nodes
                        .iter()
                        .all(|&nn| field.perm.and_then(|pm| pm.get(nn as usize)).is_some_and(|&s| s != elmer_core::NO_DOF))
                };
                let parent = match (el.parents[0], el.parents[1]) {
                    (Some(l), r) => {
                        if in_system(l) {
                            l
                        } else {
                            r.unwrap_or(l)
                        }
                    }
                    (None, Some(r)) => r,
                    (None, None) => return Err(fatal(format!("boundary element {} has no parent", e.0))),
                };
                // A p boundary element: its dofs, the scalar `ElementInfo` and
                // rule, and its nodes in the parent edge's order
                // (`AssignLocalNumber`), the nodal data following them.
                let pb = match ctx.solver_dofs() {
                    Some(dd) => {
                        let mut slots = Vec::new();
                        dd.element_slots(e, &mut slots);
                        Some((
                            slots.len(),
                            dd.element_nodes(e),
                            dd.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?,
                            dd.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?,
                        ))
                    }
                    None => None,
                };
                let bnodes: Vec<u32> = pb.as_ref().map_or_else(|| el.nodes.to_vec(), |p| p.1.clone());
                if bnodes != el.nodes {
                    let at: Vec<usize> = bnodes
                        .iter()
                        .map(|b| el.nodes.iter().position(|m| m == b).ok_or_else(|| fatal("p boundary nodes differ from the mesh's")))
                        .collect::<Result<_, _>>()?;
                    traction = at.iter().map(|&k| traction[k]).collect();
                    beta = beta.map(|v| at.iter().map(|&k| v[k]).collect());
                    spring_normal = spring_normal.map(|v| at.iter().map(|&k| v[k]).collect());
                    spring_tensor = spring_tensor.map(|v| at.iter().map(|&k| v[k]).collect());
                }
                let pel = ctx.mesh.element(parent);
                let pxyz: Vec<[f64; 3]> = pel.nodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
                let pdisp = if model.linear_model {
                    vec![[0.0; 3]; pel.nodes.len()]
                } else {
                    local_displacement(&field, pel.nodes, cdim)
                };
                let table = node_uvw(pel.etype).ok_or_else(|| unsupported(&format!("parent element type {}", pel.etype.0)))?;
                let node_local: Vec<[f64; 3]> = bnodes
                    .iter()
                    .map(|bn| {
                        pel.nodes
                            .iter()
                            .position(|pn| pn == bn)
                            .and_then(|j| table.get(j).copied())
                            .ok_or_else(|| fatal(format!("boundary element {} node {bn} is not in its parent", e.0)))
                    })
                    .collect::<Result<_, _>>()?;
                let fe = self.fe();
                let center = parent_center(&*fe.basis, &pxyz, pel.etype)?;
                let bxyz: Vec<[f64; 3]> = bnodes.iter().map(|&nn| ctx.mesh.coord(nn)).collect();
                let nt = beta.unwrap_or_else(|| vec![0.0; n]);
                let spring = if let Some(sn) = &spring_normal {
                    Spring::Normal(sn)
                } else if let Some(st) = &spring_tensor {
                    Spring::Tensor(st)
                } else {
                    Spring::None
                };
                let inp = BoundaryInput {
                    etype: ElementType(el.etype.0),
                    xyz: &bxyz,
                    cdim,
                    coord,
                    large_deflection: model.large,
                    pseudo_traction: bc.logical("pseudo-traction").unwrap_or(global_pseudo),
                    traction: &traction,
                    normal_traction: &nt,
                    spring,
                    parent: ParentInput { etype: pel.etype, xyz: &pxyz, disp: &pdisp, node_uvw: &node_local, center },
                    np: rel_np(&*fe.quad, ElementType(el.etype.0), RelOrders::of(ctx.params()).bc)?,
                    stride: field.dofs,
                };
                let mut local = LocalSystem::zeros(0);
                match &pb {
                    Some((ntot, _, basis, rule)) => {
                        if inp.np.is_some() {
                            return Err(unsupported("a relative integration order on p boundary elements"));
                        }
                        elastic_boundary_local_p(&**basis, rule, *ntot, &inp, &mut local)
                    }
                    None => elastic_boundary_local(&*fe.basis, &*fe.quad, &inp, &mut local),
                }
                .map_err(|x| fatal(x.to_string()))?;
                if transient {
                    local.mass = Some(vec![0.0; local.n * local.n]);
                    if second_order {
                        local.damp = Some(vec![0.0; local.n * local.n]);
                        ctx.default_second_order_time(&mut local, e, dt)?;
                    } else {
                        ctx.default_first_order_time(&mut local, e, dt)?;
                    }
                }
                // ElasticSolve.F90:2640-2650: the boundary loads of the first
                // iterate into `BulkRHS` too.
                if !bulk_rhs.is_empty() {
                    let mut slots = Vec::new();
                    ctx.element_slots(e, &mut slots);
                    let w = field.dofs;
                    for (k, &sl) in slots.iter().enumerate() {
                        if sl != elmer_core::NO_DOF {
                            for i in 0..w {
                                bulk_rhs[sl as usize * w + i] += local.force[w * k + i];
                            }
                        }
                    }
                }
                ctx.default_update(&local, e)?;
            }

            // A lumping case given as a displacement of the lumping boundary,
            // before the sif's own conditions (ModelLumpingDisplacements).
            let first = match lump.as_ref() {
                Some(l) if l.fix_displacement => l.displacements(ctx, iter as usize)?,
                _ => Vec::new(),
            };
            if use_umat {
                // ElasticSolve.F90:1702-1810: the Dirichlet values as
                // increments of the current displacement, then the residual
                // of the current iterate; the increment from zero.
                ctx.default_finish_boundary_assembly()?;
                let cur: Vec<f64> = ctx.primary().map(|f| f.values.to_vec()).unwrap_or_default();
                let fixed: Vec<(u32, f64)> =
                    ctx.dirichlet_values()?.into_iter().map(|(dof, v)| (dof, v - cur[dof as usize])).collect();
                ctx.apply_dirichlet_values(&fixed)?;
                if let Some(f) = ctx.primary_mut() {
                    f.values.iter_mut().for_each(|x| *x = 0.0);
                }
                if iter == 1 {
                    let mx = bulk_rhs.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    nonlin_res0 = if mx < AEPS {
                        ctx.info("ElasticSolver", "No pressure external loads ... ", 4);
                        1.0
                    } else {
                        bulk_rhs.iter().map(|x| x * x).sum::<f64>().sqrt()
                    };
                }
                let rhs_norm = ctx.system_mut().map_or(0.0, |s| s.rhs.iter().map(|x| x * x).sum::<f64>().sqrt());
                let res = rhs_norm / nonlin_res0;
                ctx.info(
                    "ElasticitySolver",
                    &format!("Residual for nonlinear iterate {}: {}", iter - 1, elmer_solver::fortran_fmt::es(res, 12, 3)),
                    5,
                );
                let tol = ctx.params().real("nonlinear system convergence tolerance").unwrap_or(0.0);
                let commit = |s: &mut Self, ctx: &mut SolverCtx<'_>, disp: &[f64]| {
                    let saved = s.umat.as_mut().expect("allocated");
                    for (name, dst) in
                        [("umatstress", &mut saved.stress0), ("umatenergy", &mut saved.energy0), ("umatstate", &mut saved.state0)]
                    {
                        if let Some(v) = ctx.sim.vars.get(name) {
                            dst.copy_from_slice(&v.values);
                        }
                    }
                    if scanning {
                        saved.prev.copy_from_slice(disp);
                    }
                };
                if res < tol && iter > min_iter {
                    ctx.info("ElasticitySolver", "Nonlinear iteration is terminated succesfully", 5);
                    if let Some(f) = ctx.primary_mut() {
                        f.values.copy_from_slice(&total);
                    }
                    commit(self, ctx, &total);
                    break;
                }
                ctx.default_solve()?;
                let new: Vec<f64> = {
                    let f = ctx.primary_mut().ok_or_else(|| fatal("solver has no displacement variable"))?;
                    for (x, t) in f.values.iter_mut().zip(&total) {
                        *x += t;
                    }
                    f.values.to_vec()
                };
                if iter == max_iter {
                    ctx.info("ElasticitySolver", "The maximum of nonlinear iterations reached: Terminating...", 5);
                    commit(self, ctx, &new);
                    break;
                }
                continue;
            }
            ctx.default_dirichlet_after(&first)?;
            ctx.default_solve()?;
            if let Some(l) = lump.as_mut() {
                let fe = self.fe();
                l.springs(ctx, &*fe.basis, &*fe.quad, iter as usize)?;
            }

            if (ctx.default_converged() || iter == max_iter) && iter >= min_iter {
                break;
            }
        }

        let flags = post_flags(ctx);
        if flags.any() {
            self.postprocess(ctx, flags, cdim, axi, model, &temp_ref)?;
        }
        Ok(())
    }
}
