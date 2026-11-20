//! `MagnetoDynamicsCalcFields` (CalcFields.F90) when its target is a nodal
//! 2D potential (`MagnetoDynamics2D`, `MagnetoDynamics2DHarmonic`): the
//! magnetic flux density `B = curl(A e_z)` (with the `A/r` term in axial
//! symmetry), the magnetic field strength (constant scalar reluctivity), the
//! magnetic vector potential, the electric field, the current density and
//! the Joule heating, real or complex (`re`/`im`),
//! as nodal fields (one L2 projection per component with the section's
//! linear solver, `GlobalSol`) and as elemental DG fields (the element mass
//! matrix inverted locally, `LocalSol`).
//!
//! [`init_section`] is `MagnetoDynamicsCalcFields_Init0` and `_Init`: the
//! target solver, `Variable = -nooutput hr_dummy` and the `Exported
//! Variable` list, so the fields exist before the solvers run, as in Elmer.
//! The solver variable keeps the last component solved, its norm the one
//! `Show Norm Index` names.
//!
//! Refused rather than skipped: the other `Calculate *` fields (field
//! strength, JxB, Maxwell stress, losses, winding voltage, permeability,
//! nodal forces and heating, flux linkage), coils and circuits, velocities,
//! tensor conductivities, transient and eigen targets, `Vtu Style`, the
//! elemental modes, DG averaging, a constant mass matrix and parallel runs.
//! The integral results (`res:` energies, powers) are only logged; a
//! SaveScalars that would save them is refused in `elmer-run`.

use std::sync::Arc;

use elmer_core::{BasisEval, ElemId, LocalSystem, ParamMap, Params, Value, NO_DOF};
use elmer_fe::metric::MetricStatus;
use elmer_solver::{Model, SolverCtx, SolverError};

const CALLER: &str = "MagnetoDynamicsCalcFields";

/// Whether `params` runs `MagnetoDynamicsCalcFields`.
pub fn is_section(params: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(params).is_some_and(|(_, e)| e.eq_ignore_ascii_case("magnetodynamicscalcfields"))
}

/// The `Procedure` string as `GetString` returns it (`"<file> <entry>"`).
fn procedure_string(p: &ParamMap) -> String {
    elmer_solver::solver::procedure_of(p).map_or_else(String::new, |(l, e)| format!("{l} {e}"))
}

/// The target of `MagnetoDynamicsCalcFields_Init0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Target {
    /// 1-based solver index (`Target Variable Solver Index`).
    pub index: usize,
    /// Components of the target variable (2 for a harmonic one).
    pub vdofs: usize,
    /// An edge-element (`WhitneyAV*`) target.
    pub hcurl: bool,
}

/// `MagnetoDynamicsCalcFields_Init0`: `Potential Variable`'s solver (its
/// components counted from the `:n` of its `Variable`), else the last
/// solver whose `Procedure` names one of the four potential solvers.
pub fn target(solvers: &[ParamMap], sp: &ParamMap) -> Result<Target, SolverError> {
    if let Some(pname) = sp.string("potential variable") {
        let pname = pname.trim();
        for (i, s) in solvers.iter().enumerate() {
            let sname = s.string("variable").unwrap_or("");
            let j = sname.find('[').filter(|&j| j > 0).unwrap_or(sname.trim_end().len());
            if sname[..j].eq_ignore_ascii_case(pname) {
                let vdofs = sname.match_indices(':').filter_map(|(k, _)| sname[k + 1..].chars().next()?.to_digit(10)).sum::<u32>();
                let hcurl = procedure_string(s).contains("WhitneyAV");
                return Ok(Target { index: i + 1, vdofs: vdofs as usize, hcurl });
            }
        }
        return Err(SolverError::fatal(
            "MagnetoDynamicsCalcFields_Init0",
            format!("Could not find solver for variable: {pname}"),
        ));
    }
    for (i, s) in solvers.iter().enumerate().rev() {
        let sname = procedure_string(s);
        let found = if sname.contains("WhitneyAVHarmonicSolver") {
            Some((2, true))
        } else if sname.contains("MagnetoDynamics2DHarmonic") {
            Some((2, false))
        } else if sname.contains("WhitneyAVSolver") {
            Some((1, true))
        } else if sname.contains("MagnetoDynamics2D") {
            Some((1, false))
        } else {
            None
        };
        if let Some((vdofs, hcurl)) = found {
            return Ok(Target { index: i + 1, vdofs, hcurl });
        }
    }
    Err(SolverError::fatal(
        "MagnetoDynamicsCalcFields_Init0",
        "Could not determine target variable type (real or complex)",
    ))
}

fn logical(p: &ParamMap, k: &str) -> Option<bool> {
    p.logical(k)
}

fn has_prefix(p: &ParamMap, prefix: &str) -> bool {
    p.0.keys().any(|k| k.starts_with(prefix))
}

/// `MagnetoDynamicsCalcFields_Init0` and `_Init` for solver section `i`
/// (nothing when Init0 finds no target: the solve then fails as Elmer's
/// Init0 does). `dim` is `CoordinateSystemDimension()`.
pub fn init_section(model: &mut Model, i: usize, dim: usize) {
    let Ok(t) = target(&model.solvers, &model.solvers[i]) else { return };
    let transient = model.simulation.string("simulation type").is_some_and(|s| s.trim().eq_ignore_ascii_case("transient"));
    let lorentz = model.body_forces.iter().any(|b| has_prefix(b, "angular velocity") || has_prefix(b, "lorentz velocity"));
    let layer_bc = model.boundary_conditions.iter().any(|b| b.get("layer electric conductivity").is_some());
    let force_components = model
        .components
        .iter()
        .any(|c| logical(c, "calculate magnetic force") == Some(true) || logical(c, "calculate magnetic torque") == Some(true));
    let target_eigen = logical(&model.solvers[t.index - 1], "eigen analysis") == Some(true);
    let sp = &mut model.solvers[i];
    let add_new = |sp: &mut ParamMap, k: &str, v: Value| {
        if sp.get(k).is_none() {
            sp.insert(k, v);
        }
    };
    // _Init0.
    if target_eigen {
        add_new(sp, "Eigen Analysis", Value::Logical(true));
    }
    add_new(sp, "Target Variable Real Field", Value::Logical(t.vdofs != 2));
    add_new(sp, "Target Variable Solver Index", Value::Integer(vec![t.index as i32]));

    // _Init.
    let fdim = if logical(sp, "2d result fields") == Some(true) && dim == 2 { dim } else { 3 };
    let mut real = logical(sp, "target variable real field").unwrap_or(false);
    if logical(sp, "vtu style") == Some(true) {
        real = true;
    }
    let complex = !real;
    sp.insert("Variable", Value::Str("-nooutput hr_dummy".into()));
    sp.insert("Linear System refactorize", Value::Logical(false));
    let mut list: Vec<String> = Vec::new();
    let vec_field = |name: &str, n: usize, e: &str| {
        if real {
            format!("{name}{e}[{name}{e}:{n}]")
        } else {
            format!("{name}{e}[{name} re{e}:{n} {name} im{e}:{n}]")
        }
    };
    let on = |sp: &ParamMap, k: &str| logical(sp, k) == Some(true);
    let nodal = match logical(sp, "skip nodal fields") {
        Some(skip) => !skip,
        None => logical(sp, "calculate nodal fields").unwrap_or(true),
    };
    if nodal {
        if logical(sp, "calculate magnetic flux density").unwrap_or(true) {
            list.push(vec_field("Magnetic Flux Density", fdim, ""));
        }
        if on(sp, "calculate magnetic vector potential") {
            list.push(vec_field("Magnetic Vector Potential", 3, ""));
        }
        if on(sp, "calculate magnetic field strength") {
            list.push(vec_field("Magnetic Field Strength", fdim, ""));
        }
        if on(sp, "calculate jxb") {
            list.push(vec_field("JxB", 3, ""));
        }
        if on(sp, "calculate maxwell stress") {
            list.push(vec_field("Maxwell Stress", 6, ""));
        }
        if on(sp, "calculate current density") {
            list.push(vec_field("Current Density", 3, ""));
        }
        if on(sp, "calculate joule heating") {
            list.push("Joule Heating".into());
        }
        if on(sp, "calculate harmonic loss") && !real {
            list.push("Harmonic Loss Linear".into());
            list.push("Harmonic Loss Quadratic".into());
        }
        if on(sp, "calculate homogenization loss") && !real {
            list.push("Proximity Loss".into());
        }
        if transient || complex || lorentz {
            if on(sp, "calculate electric field") {
                list.push(vec_field("Electric Field", 3, ""));
            }
            if on(sp, "calculate winding voltage") {
                list.push(if real { "Winding Voltage".into() } else { vec_field("Winding Voltage", 1, "") });
            }
        }
        if on(sp, "calculate relative permeability") {
            list.push(if real { "Relative Permeability".into() } else { vec_field("Relative Permeability", 1, "") });
        }
        if on(sp, "calculate electric scalar potential") && real {
            list.push("Electric Scalar Potential".into());
        }
    }
    if on(sp, "calculate current density") && layer_bc && !real {
        list.push("Surface Current[Surface Current re:3 Surface Current im:3]".into());
    }
    if on(sp, "calculate nodal heating") {
        list.push("Nodal Joule Heating".into());
    }
    if force_components || on(sp, "calculate nodal forces") {
        list.push(format!("Nodal Force[Nodal Force:{fdim}]"));
    }
    if !on(sp, "discontinuous galerkin") {
        let elemental = match logical(sp, "skip elemental fields") {
            Some(skip) => !skip,
            None => logical(sp, "calculate elemental fields").unwrap_or(true),
        };
        if elemental {
            let dg = |s: String| format!("-dg {s}");
            if logical(sp, "calculate magnetic flux density").unwrap_or(true) {
                list.push(dg(vec_field("Magnetic Flux Density", fdim, " E")));
            }
            if on(sp, "calculate magnetic vector potential") {
                list.push(dg(vec_field("Magnetic Vector Potential", 3, " E")));
            }
            if on(sp, "calculate magnetic field strength") {
                list.push(dg(vec_field("Magnetic Field Strength", fdim, " E")));
            }
            if on(sp, "calculate jxb") {
                list.push(dg(vec_field("JxB", 3, " E")));
            }
            if on(sp, "calculate maxwell stress") {
                list.push(dg(vec_field("Maxwell Stress", 6, " E")));
            }
            if on(sp, "calculate current density") {
                list.push(dg(vec_field("Current Density", 3, " E")));
            }
            if on(sp, "calculate joule heating") {
                list.push("-dg Joule Heating E".into());
            }
            if on(sp, "calculate harmonic loss") && !real {
                list.push("-dg Harmonic Loss Linear E".into());
                list.push("-dg Harmonic Loss Quadratic E".into());
            }
            if on(sp, "calculate homogenization loss") && !real {
                list.push("-dg Proximity Loss E".into());
            }
            if transient || complex || lorentz {
                if on(sp, "calculate electric field") {
                    list.push(dg(vec_field("Electric Field", 3, " E")));
                }
                if on(sp, "calculate winding voltage") {
                    list.push(if real { "-dg Winding Voltage E".into() } else { dg(vec_field("Winding Voltage", 1, " E")) });
                }
            }
            if on(sp, "calculate relative permeability") {
                list.push(if real {
                    "-dg Relative Permeability E".into()
                } else {
                    dg(vec_field("Relative Permeability", 1, " E"))
                });
            }
            if force_components || on(sp, "calculate nodal forces") {
                list.push(format!("-dg Nodal Force E[Nodal Force E:{fdim}]"));
            }
        }
    }
    let mut k = 1;
    while sp.get(&format!("exported variable {k}")).is_some() {
        k += 1;
    }
    for s in list {
        sp.insert(&format!("Exported Variable {k}"), Value::Str(s));
        k += 1;
    }
}

/// `LUDecomp` (LinearAlgebra.F90:191) of the column-major `n x n` matrix:
/// column pivoting, the diagonal replaced by its reciprocal.
fn lu_decomp(a: &mut [f64], n: usize, pivot: &mut [usize]) -> bool {
    let at = |i: usize, j: usize| j * n + i;
    for i in 0..n {
        let mut j = i;
        for k in i + 1..n {
            if a[at(i, k)].abs() > a[at(i, j)].abs() {
                j = k;
            }
        }
        if a[at(i, j)].abs() == 0.0 {
            return false;
        }
        pivot[i] = j;
        if j != i {
            for k in 0..=i {
                a.swap(at(k, j), at(k, i));
            }
        }
        for k in i + 1..n {
            a[at(i, k)] /= a[at(i, i)];
        }
        for k in i + 1..n {
            if j != i {
                a.swap(at(k, i), at(k, j));
            }
            for l in i + 1..n {
                a[at(k, l)] -= a[at(k, i)] * a[at(i, l)];
            }
        }
    }
    pivot[n - 1] = n - 1;
    for i in 0..n {
        if a[at(i, i)].abs() == 0.0 {
            return false;
        }
        a[at(i, i)] = 1.0 / a[at(i, i)];
    }
    true
}

/// `LUSolve` (LinearAlgebra.F90) with the pivots of [`lu_decomp`].
fn lu_solve(n: usize, a: &[f64], x: &mut [f64], pivot: &[usize]) {
    let at = |i: usize, j: usize| j * n + i;
    for i in 0..n {
        let mut s = x[i];
        for j in 0..i {
            s -= a[at(i, j)] * x[j];
        }
        x[i] = a[at(i, i)] * s;
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        for j in i + 1..n {
            s -= a[at(i, j)] * x[j];
        }
        x[i] = s;
    }
    for i in (0..n).rev() {
        if pivot[i] != i {
            x.swap(i, pivot[i]);
        }
    }
}

/// `SUM(a(1:n) * b(1:n))`.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

/// One of the computed fields: the nodal and elemental variables Elmer
/// finds by name (`VariableGet`), and its number of components.
struct Field {
    nodal: Option<elmer_solver::VarId>,
    elem: Option<elmer_solver::VarId>,
    m: usize,
}

impl Field {
    fn any(&self) -> bool {
        self.nodal.is_some() || self.elem.is_some()
    }
}

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {}", what.into()))
}

/// `ListGetAngularFrequency(pSolver % Values, Found, Element)`.
fn angular_frequency(ctx: &SolverCtx<'_>, psp: &ParamMap, e: ElemId) -> Result<Option<f64>, SolverError> {
    let two_pi = 2.0 * std::f64::consts::PI;
    let bp = ctx.body_params(e);
    let lists = [
        Some(psp),
        bp.as_ref().and_then(|b| b.equation),
        bp.as_ref().and_then(|b| b.material),
        Some(ctx.simulation_params()),
        Some(ctx.params()),
    ];
    for p in lists.into_iter().flatten() {
        if let Some(f) = ctx.creal(p, "frequency")? {
            return Ok(Some(two_pi * f));
        }
        if let Some(w) = ctx.creal(p, "angular frequency")? {
            return Ok(Some(w));
        }
    }
    Ok(None)
}

/// A scalar (1x1 or dependent) keyword; tensors are refused.
fn scalar_at_nodes(ctx: &SolverCtx<'_>, p: &ParamMap, k: &str, e: ElemId, n: usize) -> Result<Vec<f64>, SolverError> {
    if let Some(Value::Real { rows, cols, .. }) = p.get(k) {
        if *rows * *cols != 1 {
            return Err(unsupported(format!("a tensor '{k}'")));
        }
    }
    Ok(ctx.real_at_nodes(p, k, e)?.unwrap_or_else(|| vec![0.0; n]))
}

/// `MagnetoDynamicsCalcFields` for the nodal 2D target `t`.
pub fn solve(ctx: &mut SolverCtx<'_>, t: Target, transient: bool) -> Result<(), SolverError> {
    let sp = ctx.params().clone();
    if ctx.comm.size() > 1 {
        return Err(unsupported("parallel runs"));
    }
    if transient {
        return Err(unsupported("a transient run"));
    }
    if sp.integer("target variable solver index").is_none() {
        return Err(unsupported("a section without the MagnetoDynamicsCalcFields_Init edits"));
    }
    for k in [
        "eigen analysis",
        "vtu style",
        "constant mass matrix",
        "average within materials",
        "discontinuous galerkin",
        "impose circuit current",
        "precomputed electric potential",
        "loss estimation",
        "calculate flux linkage",
        "separate magnetic energy",
    ] {
        if logical(&sp, k) == Some(true) {
            return Err(unsupported(format!("'{k}'")));
        }
    }
    if sp.integer("calculate elemental mode").is_some_and(|m| m != 0) {
        return Err(unsupported("'Calculate Elemental Mode'"));
    }
    let m = &ctx.sim.model;
    if m.body_forces.iter().any(|b| has_prefix(b, "angular velocity") || has_prefix(b, "lorentz velocity")) {
        return Err(unsupported("velocity body forces"));
    }
    if m.components.iter().any(|c| c.get("coil type").is_some()) {
        return Err(unsupported("coil components"));
    }
    if m.boundary_conditions.iter().any(|b| logical(b, "magnetic flux average") == Some(true)) {
        return Err(unsupported("'Magnetic Flux Average'"));
    }
    if m.materials.iter().any(|mt| mt.get("reluctivity function").is_some()) && t.vdofs == 2 {
        return Err(SolverError::fatal(CALLER, "Reluctivity Function not implemented in complex cases!"));
    }
    let psp = m
        .solvers
        .get(t.index - 1)
        .cloned()
        .ok_or_else(|| SolverError::fatal(CALLER, "no target solver section"))?;
    if logical(&psp, "eigen analysis") == Some(true) || psp.get("matrix multiplier name").is_some() {
        return Err(unsupported("an eigen or topology-optimisation target"));
    }
    let axi = match ctx.simulation_params().string("coordinate system").map(elmer_core::normalize_key) {
        None => false,
        Some(s) if s == "cartesian" || s == "cartesian 2d" => false,
        Some(s) if s == "axi symmetric" || s == "cylindric symmetric" => true,
        Some(s) => return Err(unsupported(format!("Coordinate System = {s}"))),
    };
    let dim = ctx.mesh.dim();
    if t.hcurl && dim != 3 {
        return Err(SolverError::fatal(CALLER, "an edge-element target needs a 3D mesh"));
    }
    if !t.hcurl && dim != 2 {
        return Err(unsupported("a nodal target on a mesh that is not 2D"));
    }
    // The target potential.
    let pvar = ctx
        .sim
        .slots
        .iter()
        .find(|s| s.id == t.index)
        .and_then(|s| s.variable)
        .ok_or_else(|| SolverError::fatal(CALLER, "the target solver has no variable"))?;
    // `vDOFs = pSolver % Variable % DOFs` (Init0's count only decides the
    // exported names).
    let vdofs = ctx.sim.vars.by_id(pvar).dofs;
    if !(1..=2).contains(&vdofs) {
        return Err(unsupported(format!("a target variable of {vdofs} components")));
    }
    let real = vdofs == 1;
    let fdim = if logical(&sp, "2d result fields") == Some(true) { dim } else { 3 };
    let norm_index = sp.integer("show norm index").map_or(0, |i| i.max(0) as usize);
    let harm_power = if logical(&sp, "calculate harmonic peak power") == Some(true) { 1.0 } else { 0.5 };

    // The fields Elmer finds by name; the ones not ported are refused.
    let vars = &ctx.sim.vars;
    let get = |n: &str| vars.id(n);
    for n in [
        "jxb",
        "winding voltage",
        "relative permeability",
        "maxwell stress",
        "nodal force",
        "nodal joule heating",
        "harmonic loss linear",
        "harmonic loss quadratic",
        "proximity loss",
        "electric scalar potential",
        "surface current",
    ] {
        if get(n).is_some() || get(&format!("{n} e")).is_some() {
            return Err(unsupported(format!("the field '{n}'")));
        }
    }
    let field = |n: &str, m: usize| Field { nodal: get(n), elem: get(&format!("{n} e")), m };
    let mfd = field("magnetic flux density", fdim * vdofs);
    let mfs = field("magnetic field strength", fdim * vdofs);
    let vp = field("magnetic vector potential", 3 * vdofs);
    let ef = if real {
        Field { nodal: None, elem: None, m: 3 }
    } else {
        field("electric field", 3 * vdofs)
    };
    let cd = field("current density", 3 * vdofs);
    let jh = field("joule heating", 1);
    let fields = [&mfd, &mfs, &vp, &ef, &cd, &jh];
    if mfs.any() {
        let m = &ctx.sim.model;
        if m.body_forces.iter().chain(&m.materials).any(|b| has_prefix(b, "magnetization")) {
            return Err(unsupported("Magnetic Field Strength with a magnetization"));
        }
    }
    let ncol: usize = fields.iter().filter(|f| f.any()).map(|f| f.m).sum();
    let nodal_fields = fields.iter().any(|f| f.nodal.is_some());
    let elemental_fields = fields.iter().any(|f| f.elem.is_some());

    let id = ctx.solver_id();
    let ix = ctx.sim.slots.iter().position(|s| s.id == id).expect("running solver has a slot");
    let slot_perm = ctx.sim.slots[ix].perm.clone();
    let hr = ctx.sim.slots[ix].variable.ok_or_else(|| SolverError::fatal(CALLER, "the solver has no variable"))?;
    for f in &fields {
        if let Some(v) = f.nodal {
            let var = ctx.sim.vars.by_id(v);
            let same = match (&var.perm, &slot_perm) {
                (Some(a), Some(b)) => Arc::ptr_eq(a, b) || a[..] == b[..],
                _ => false,
            };
            if !same || var.dofs != f.m {
                return Err(unsupported(format!("the field '{}' on another permutation", var.name)));
            }
        }
        if let Some(v) = f.elem {
            let var = ctx.sim.vars.by_id(v);
            if var.kind != elmer_core::FieldKind::Dg || var.dofs != f.m {
                return Err(unsupported(format!("the field '{}' that is not a DG field", var.name)));
            }
        }
    }

    let impose_bf_current = match logical(&sp, "impose body force current") {
        Some(b) => b,
        None => ctx.sim.model.body_forces.iter().any(|b| has_prefix(b, "current density")),
    };

    // An edge-element target: its DOF layout (`WhitneyAV*` keep their own,
    // nodal V first, then the edges).
    let edge = if t.hcurl {
        if axi {
            return Err(unsupported("an edge-element target in axial symmetry"));
        }
        for k in ["use piola transform", "quadratic approximation"] {
            if logical(&psp, k) == Some(true) {
                return Err(unsupported(format!("an edge-element target with '{k}'")));
            }
        }
        let bf_potential = match logical(&sp, "impose body force potential") {
            Some(b) => b,
            None => ctx.sim.model.body_forces.iter().any(|b| b.get("electric potential").is_some()),
        };
        if bf_potential {
            return Err(unsupported("'Impose Body Force Potential'"));
        }
        if impose_bf_current && cd.any() {
            return Err(unsupported("body-force current density in the current density of an edge-element target"));
        }
        let case = crate::Case::from_view(ctx.mesh, &ctx.sim.model);
        let with_v = vdofs == 2 || crate::steady::static_conductivity(&case, &psp);
        let lay = crate::layout::AvLayout::build(&case, t.index, with_v)?;
        if ctx.sim.vars.by_id(pvar).values.len() != vdofs * lay.n() {
            return Err(SolverError::fatal(CALLER, "the potential has not been solved on the edge layout yet"));
        }
        Some((case, lay))
    } else {
        None
    };

    ctx.default_initialize()?;
    let n_rows = ctx.system().map_or(0, |s| s.rhs.len());
    let mut gforce = if nodal_fields { vec![vec![0.0; n_rows]; ncol] } else { Vec::new() };
    let elems: Vec<ElemId> = ctx.active_elements().filter(|e| (e.0 as usize) < ctx.mesh.n_bulk()).collect();
    let mut slots = Vec::new();
    let mut power = 0.0;
    let mut omega = 0.0;
    let mut status = MetricStatus::default();
    for &e in &elems {
        let el = ctx.mesh.element(e);
        let etype = el.etype;
        let nodes: Vec<u32> = el.nodes.to_vec();
        let n = nodes.len();
        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
        let eq_n = n;
        let bp = ctx.body_params(e);
        let material = bp.as_ref().and_then(|b| b.material).cloned();
        let body_force = bp.as_ref().and_then(|b| b.body_force).cloned();
        omega = angular_frequency(ctx, &psp, e)?.unwrap_or(0.0);
        // BodyForceCurrDens(3, 1:n), complex.
        let mut bf_cd = (vec![0.0; n], vec![0.0; n]);
        if impose_bf_current {
            if let Some(bf) = &body_force {
                bf_cd.0 = scalar_at_nodes(ctx, bf, "current density", e, n)?;
                if vdofs == 2 {
                    bf_cd.1 = scalar_at_nodes(ctx, bf, "current density im", e, n)?;
                }
                if let Some(mult) = ctx.creal(bf, "current density multiplier")? {
                    bf_cd.0.iter_mut().chain(bf_cd.1.iter_mut()).for_each(|v| *v *= mult);
                }
            }
        }
        // GetCMPLXElectricConductivityTensor: a scalar on the diagonal.
        let (c_re, c_im) = match &material {
            Some(mt) => (
                scalar_at_nodes(ctx, mt, "electric conductivity", e, n)?,
                scalar_at_nodes(ctx, mt, "electric conductivity im", e, n)?,
            ),
            None => (vec![0.0; n], vec![0.0; n]),
        };

        let mut mass = vec![0.0; eq_n * eq_n];
        let mut force = vec![vec![0.0; eq_n]; ncol.max(1)];
        // GetReluctivity: a constant scalar (`Nu = diag(R_ip)`).
        let nu = match (&material, mfs.any()) {
            (Some(mt), true) => crate::material::reluctivity(mt, ctx.constants())?,
            (None, true) => return Err(SolverError::fatal(CALLER, format!("element {} has no material", e.0))),
            _ => 0.0,
        };
        let nu_nodes = vec![nu; n];
        let (ips, np) = match &edge {
            None => (nodal_points(ctx, pvar, vdofs, &nodes, &xyz, etype, axi, omega, &mut status)?, n),
            Some((case, lay)) => edge_points(ctx, pvar, vdofs, case, lay, e, omega, &mut status)?,
        };
        for ip in &ips {
            let (s, basis, b, ev, vp_ip) = (ip.s, &ip.basis[..n], ip.b, ip.e, ip.vp);
            let mut cmat_re = [[0.0f64; 3]; 3];
            let mut cmat_im = [[0.0f64; 3]; 3];
            let c_re_ip = dot(&c_re, basis);
            let c_im_ip = dot(&c_im, basis);
            for k in 0..3 {
                cmat_re[k][k] = c_re_ip;
                cmat_im[k][k] = c_im_ip;
            }
            let bf_ip = if impose_bf_current { [dot(&bf_cd.0, basis), dot(&bf_cd.1, basis)] } else { [0.0, 0.0] };
            // JatIp (the circuit term CC_J is zero without circuits).
            let mut jat = [[0.0f64; 3]; 2];
            for l in 0..3 {
                if real {
                    jat[0][l] = dot(&cmat_re[l], &ev[0]) + 0.0 + bf_re(l, bf_ip);
                } else {
                    jat[0][l] = dot(&cmat_re[l], &ev[0]) - dot(&cmat_im[l], &ev[1]) + bf_re(l, bf_ip);
                    jat[1][l] = dot(&cmat_im[l], &ev[0]) + dot(&cmat_re[l], &ev[1]) + bf_im(l, bf_ip);
                }
            }
            let joule_from_current = np == 0;
            for p in 0..eq_n {
                for q in 0..eq_n {
                    mass[q * eq_n + p] += s * basis[p] * basis[q];
                }
            }
            // MATMUL(C, E^T) summed against E^T.
            let quad = |c: &[[f64; 3]; 3], x: &[f64; 3], y: &[f64; 3]| {
                let mut tot = 0.0;
                for r in 0..3 {
                    let mut v = 0.0;
                    for l in 0..3 {
                        v += c[r][l] * x[l];
                    }
                    tot += v * y[r];
                }
                tot
            };
            for p in 0..eq_n {
                let mut k = 0;
                if mfd.any() {
                    for bl in &b[..vdofs] {
                        for c in 0..fdim {
                            force[k + c][p] += s * bl[c] * basis[p];
                        }
                        k += fdim;
                    }
                }
                if mfs.any() {
                    // MATMUL(REAL(Nu), B) - MATMUL(AIMAG(Nu), B_im) - REAL(MG_ip),
                    // and the imaginary part, with Nu real and MG = 0.
                    let r_ip = dot(basis, &nu_nodes);
                    let nmv = |d: f64, x: &[f64; 3], c: usize| {
                        let mut v = 0.0;
                        for (l, xl) in x[..fdim].iter().enumerate() {
                            v += if l == c { d } else { 0.0 } * xl;
                        }
                        v
                    };
                    for c in 0..fdim {
                        let h = if real { nmv(r_ip, &b[0], c) - 0.0 } else { nmv(r_ip, &b[0], c) - nmv(0.0, &b[1], c) - 0.0 };
                        force[k + c][p] += s * h * basis[p];
                    }
                    k += fdim;
                    if !real {
                        for c in 0..fdim {
                            let h = nmv(0.0, &b[0], c) + nmv(r_ip, &b[1], c) - 0.0;
                            force[k + c][p] += s * h * basis[p];
                        }
                        k += fdim;
                    }
                }
                if vp.any() {
                    for vl in &vp_ip[..vdofs] {
                        for c in 0..3 {
                            force[k + c][p] += s * vl[c] * basis[p];
                        }
                        k += 3;
                    }
                }
                if ef.any() {
                    for el in &ev[..vdofs] {
                        for c in 0..3 {
                            force[k + c][p] += s * el[c] * basis[p];
                        }
                        k += 3;
                    }
                }
                if cd.any() {
                    for jl in &jat[..vdofs] {
                        for l in 0..3 {
                            force[k + l][p] += s * jl[l] * basis[p];
                        }
                        k += 3;
                    }
                }
                let coeff = if real && joule_from_current && cd.any() {
                    // J.J / sigma.
                    let mut c = 0.0;
                    for l in 0..3 {
                        if cmat_re[l][l] > f64::EPSILON {
                            c += jat[0][l] * jat[0][l] / cmat_re[l][l] * basis[p] * s;
                        }
                    }
                    c
                } else if real {
                    quad(&cmat_re, &ev[0], &ev[0]) * basis[p] * s
                } else {
                    harm_power
                        * (quad(&cmat_re, &ev[0], &ev[0]) * basis[p] * s - quad(&cmat_im, &ev[1], &ev[0]) * basis[p] * s
                            + quad(&cmat_im, &ev[0], &ev[1]) * basis[p] * s
                            + quad(&cmat_re, &ev[1], &ev[1]) * basis[p] * s)
                };
                power += coeff;
                if jh.any() {
                    force[k][p] += coeff;
                    k += 1;
                }
                debug_assert!(k <= ncol);
            }
        }

        if nodal_fields {
            let mut stiff = vec![0.0; eq_n * eq_n];
            for p in 0..eq_n {
                for q in 0..eq_n {
                    stiff[p * eq_n + q] = mass[q * eq_n + p];
                }
            }
            let ls = LocalSystem { n: eq_n, stiff, mass: None, damp: None, force: force[0].clone() };
            ctx.default_update(&ls, e)?;
            ctx.element_slots(e, &mut slots);
            for (c, f) in force.iter().enumerate().take(ncol) {
                for (p, &r) in slots.iter().enumerate().take(eq_n) {
                    if r != NO_DOF {
                        gforce[c][r as usize] += f[p];
                    }
                }
            }
        }

        if elemental_fields {
            let mut pivot = vec![0usize; eq_n];
            if !lu_decomp(&mut mass, eq_n, &mut pivot) {
                return Err(SolverError::fatal(CALLER, "LU-decomposition fails"));
            }
            let dg = ctx.sim.ensure_dg_indexes(ctx.mesh).index[e.0 as usize].clone();
            let mut dofs = 0;
            for f in &fields {
                let Some(v) = f.elem else { continue };
                let var = ctx.sim.vars.by_id_mut(v);
                let perm = var.perm.clone().expect("a DG field has a permutation");
                let ind: Vec<u32> = dg.iter().map(|&d| if d == NO_DOF { NO_DOF } else { perm[d as usize] }).collect();
                if ind.contains(&NO_DOF) {
                    continue;
                }
                for i in 0..f.m {
                    let mut x = force[dofs].clone();
                    dofs += 1;
                    lu_solve(eq_n, &mass, &mut x, &pivot);
                    for (a, &r) in ind.iter().enumerate() {
                        var.values[var.dofs * r as usize + i] = x[a];
                    }
                    var.values_changed = true;
                }
            }
        }
    }
    let _ = omega;
    ctx.info(CALLER, &format!("Eddy current power: {}", elmer_solver::fortran_fmt::es(power, 15, 6)), 4);

    // GlobalSol, field by field.
    let mut save_norm = 0.0;
    if nodal_fields {
        // `_Init` sets `Linear System Refactorize = False`: Elmer factors
        // the preconditioner on the first solve and reuses it for the other
        // components. The matrix is the same for all of them, so the
        // preconditioner recomputed for each solve is the same one (up to
        // the rounding of the scaling round trip between solves).
        let key = elmer_core::normalize_key("linear system refactorize");
        let refactorize = ctx.params_mut().0.remove(&key);
        let save_rhs = ctx.system().map(|s| s.rhs.clone());
        let mut dofs = 0;
        for f in &fields {
            let Some(v) = f.nodal else {
                if f.elem.is_some() {
                    dofs += f.m;
                }
                continue;
            };
            let name = ctx.sim.vars.by_id(v).name.clone();
            ctx.info(CALLER, &format!("Solving for field: {name}"), 6);
            for i in 0..f.m {
                let rhs = std::mem::take(&mut gforce[dofs]);
                dofs += 1;
                if let Some(s) = ctx.system_mut() {
                    s.rhs = rhs;
                }
                ctx.sim.vars.by_id_mut(hr).values.iter_mut().for_each(|x| *x = 0.0);
                let norm = ctx.default_solve()?;
                let x = ctx.sim.vars.by_id(hr).values.clone();
                let var = ctx.sim.vars.by_id_mut(v);
                for (r, xr) in x.iter().enumerate() {
                    var.values[f.m * r + i] = *xr;
                }
                var.values_changed = true;
                if norm_index == dofs {
                    save_norm = norm;
                }
            }
        }
        if let (Some(s), Some(r)) = (ctx.system_mut(), save_rhs) {
            s.rhs = r;
        }
        if let Some(v) = refactorize {
            ctx.params_mut().0.insert(key, v);
        }
    }
    if norm_index > 0 {
        ctx.info(CALLER, &format!("Reverting norm to: {save_norm}"), 4);
        ctx.sim.vars.by_id_mut(hr).conv.norm = save_norm;
    }
    Ok(())
}

/// One integration point of the element loop: `s = IP % s * detJ` (times
/// `x` in axial symmetry), the nodal basis, and `B`, `E` and the vector
/// potential, real and imaginary parts.
struct Ip {
    s: f64,
    basis: Vec<f64>,
    b: [[f64; 3]; 2],
    e: [[f64; 3]; 2],
    vp: [[f64; 3]; 2],
}

/// The 2D nodal kernel: `GaussPoints(Element)`, `B = curl(A e_z)`, the
/// harmonic `E_z = (Omega A_im, -Omega A_re)`; a real steady target has
/// `E = -dA/dt = 0`.
#[allow(clippy::too_many_arguments)]
fn nodal_points(
    ctx: &SolverCtx<'_>,
    pvar: elmer_solver::VarId,
    vdofs: usize,
    nodes: &[u32],
    xyz: &[[f64; 3]],
    etype: elmer_core::ElementType,
    axi: bool,
    omega: f64,
    status: &mut MetricStatus,
) -> Result<Vec<Ip>, SolverError> {
    let n = nodes.len();
    // GetVectorLocalSolution(SOL, Pname, uSolver=pSolver).
    let sol: Vec<Vec<f64>> = {
        let pv = ctx.sim.vars.by_id(pvar);
        (0..vdofs)
            .map(|k| nodes.iter().map(|&nd| pv.slot(nd).map_or(0.0, |s| pv.values[vdofs * s + k])).collect())
            .collect()
    };
    let def = elmer_fe::element_type(etype.0).map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
    let rule = elmer_fe::quadrature::gauss_points(etype, None).map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
    let xs: Vec<f64> = xyz.iter().map(|c| c[0]).collect();
    let mut out = Vec::with_capacity(rule.len());
    for j in 0..rule.len() {
        let mut be = BasisEval::default();
        elmer_fe::info::element_info(def, xyz, [rule.u[j], rule.v[j], rule.w[j]], 2, status, &mut be)
            .map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
        let basis = &be.basis[..n];
        let dbx: Vec<f64> = be.dbasis_dx[..n].iter().map(|d| d[0]).collect();
        let dby: Vec<f64> = be.dbasis_dx[..n].iter().map(|d| d[1]).collect();
        let mut s = rule.s[j] * be.det_j;
        let mut xcoord = 0.0;
        if axi {
            xcoord = dot(basis, &xs);
            s *= xcoord;
        }
        let mut b = [[0.0f64; 3]; 2];
        for k in 0..vdofs {
            if axi {
                b[k][0] = -dot(&sol[k], &dby);
                b[k][1] = dot(&sol[k], &dbx) + dot(&sol[k], basis) / xcoord;
            } else {
                b[k][0] = dot(&sol[k], &dby);
                b[k][1] = -dot(&sol[k], &dbx);
            }
        }
        let mut e = [[0.0f64; 3]; 2];
        if vdofs == 2 {
            e[0][2] = omega * dot(&sol[1], basis);
            e[1][2] = -omega * dot(&sol[0], basis);
        } else {
            // -SUM(PSOL * Basis) with PSOL = 0 outside transient runs.
            e[0][2] = -dot(&vec![0.0; n], basis);
        }
        let mut vp = [[0.0f64; 3]; 2];
        for l in 0..vdofs {
            vp[l][2] = dot(&sol[l], basis);
        }
        out.push(Ip { s, basis: basis.to_vec(), b, e, vp });
    }
    Ok(out)
}

/// The edge-element kernel (lowest-order Whitney basis): `GaussPoints(
/// Element, EdgeBasis=.TRUE.)`, `B = MATMUL(SOL(np+1:nd), RotWBasis)`, the
/// harmonic `E = (Omega A_im, -Omega A_re) - grad V`, the real steady
/// `E = -grad V` (zero without the nodal V, `np = 0`). Returns the points
/// and `np`.
#[allow(clippy::too_many_arguments)]
fn edge_points(
    ctx: &SolverCtx<'_>,
    pvar: elmer_solver::VarId,
    vdofs: usize,
    case: &crate::Case,
    lay: &crate::layout::AvLayout,
    e: ElemId,
    omega: f64,
    status: &mut MetricStatus,
) -> Result<(Vec<Ip>, usize), SolverError> {
    let el = &case.mesh.bulk[e.0 as usize];
    let n = el.nodes.len();
    let eids = lay
        .edges
        .of_element(el.code, &el.nodes)
        .ok_or_else(|| SolverError::fatal(CALLER, "a CalcFields element outside the potential's elements"))?;
    let np = if lay.with_v { n } else { 0 };
    let pv = ctx.sim.vars.by_id(pvar);
    // The element's slots: the target's framework layout (`GetElementDOFs`
    // of `Element = n:1 e:1`, nodes then edges, as elmer-run sets it), else
    // its own AvLayout.
    let target = ctx.sim.slots.iter().position(|s| s.variable == Some(pvar)).and_then(|ix| ctx.sim.slot_dofs(ix));
    let slots: Vec<u32> = match &target {
        Some(d) => {
            let mut v = Vec::new();
            d.element_slots(e, &mut v);
            v
        }
        None => lay.dofs_of(&el.nodes, &eids),
    };
    if slots.len() != np + eids.len() || slots.contains(&elmer_core::NO_DOF) {
        return Err(SolverError::fatal(CALLER, "the potential's element DOFs do not match its layout"));
    }
    let mut solv = vec![vec![0.0; np]; vdofs];
    let mut sole = vec![vec![0.0; eids.len()]; vdofs];
    for k in 0..vdofs {
        for i in 0..np {
            solv[k][i] = pv.values[vdofs * slots[i] as usize + k];
        }
        for j in 0..eids.len() {
            sole[k][j] = pv.values[vdofs * slots[np + j] as usize + k];
        }
    }
    let (xyz, gindex) = crate::layout::geometry(case, &el.nodes);
    let rule = crate::local::edge_rule(el.code)?;
    let mut out = Vec::with_capacity(rule.len());
    // MATMUL(x(1:m), M(1:m, :)).
    let mv = |x: &[f64], m: &[[f64; 3]]| {
        let mut r = [0.0; 3];
        for (c, rc) in r.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (xi, mi) in x.iter().zip(m) {
                acc += xi * mi[c];
            }
            *rc = acc;
        }
        r
    };
    for t in 0..rule.len() {
        let pt = crate::local::eval_point(el.code, &xyz, &gindex, [rule.u[t], rule.v[t], rule.w[t]], status)?;
        let s = rule.s[t] * pt.det_j;
        let mut b = [[0.0f64; 3]; 2];
        let mut vp = [[0.0f64; 3]; 2];
        for k in 0..vdofs {
            b[k] = mv(&sole[k], &pt.curl);
            vp[k] = mv(&sole[k], &pt.w);
        }
        let mut ef = [[0.0f64; 3]; 2];
        if vdofs == 2 {
            let a_im = mv(&sole[1], &pt.w);
            let a_re = mv(&sole[0], &pt.w);
            let g_re = mv(&solv[0], &pt.dbasis[..np]);
            let g_im = mv(&solv[1], &pt.dbasis[..np]);
            for c in 0..3 {
                ef[0][c] = omega * a_im[c] - g_re[c];
                ef[1][c] = -omega * a_re[c] - g_im[c];
            }
        } else if np > 0 {
            let g = mv(&solv[0], &pt.dbasis[..np]);
            for c in 0..3 {
                ef[0][c] = -g[c];
            }
        }
        out.push(Ip { s, basis: pt.basis[..n].to_vec(), b, e: ef, vp });
    }
    Ok((out, np))
}

fn bf_re(l: usize, bf: [f64; 2]) -> f64 {
    if l == 2 {
        bf[0]
    } else {
        0.0
    }
}

fn bf_im(l: usize, bf: [f64; 2]) -> f64 {
    if l == 2 {
        bf[1]
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lu_matches_a_direct_solve() {
        // Column-major [[4, 1, 2], [1, 5, 3], [2, 3, 6]] and b = A [1, 2, 3].
        let a0 = [4.0, 1.0, 2.0, 1.0, 5.0, 3.0, 2.0, 3.0, 6.0];
        let mut a = a0;
        let mut piv = [0usize; 3];
        assert!(lu_decomp(&mut a, 3, &mut piv));
        let mut x = [12.0, 20.0, 26.0];
        lu_solve(3, &a, &mut x, &piv);
        for (xi, want) in x.iter().zip([1.0, 2.0, 3.0]) {
            assert!((xi - want).abs() < 1e-13, "{x:?}");
        }
    }

    #[test]
    fn target_follows_init0() {
        let mut s1 = ParamMap::new();
        s1.insert("Procedure", Value::Str("\"MagnetoDynamics2D\" \"MagnetoDynamics2DHarmonic\"".into()));
        let mut s2 = ParamMap::new();
        s2.insert("Procedure", Value::Str("\"MagnetoDynamics\" \"MagnetoDynamicsCalcFields\"".into()));
        let solvers = vec![s1, s2.clone()];
        assert_eq!(target(&solvers, &s2).unwrap(), Target { index: 1, vdofs: 2, hcurl: false });
    }
}
