//! The steady-state `WhitneyAVSolver` driver: `DoSolve` (WhitneyAVSolver.F90:781)
//! for one nonlinear iteration of a linear material, with the Elmer linear
//! solve (`SolveLinearSystem` scaling and the `Linear System *` solver
//! choice of `elmer-solver`) and the `Norm Permutation` norm.
//!
//! A linear material converges in one solve: further nonlinear iterations
//! start from the solution and stop at once, so one solve is what Elmer
//! reports.

use elmer_core::{CsrMatrix, NullLogger, ParamMap, Params, SerialComm};
use elmer_solver::backend::LinearSolverChoice;

use crate::case::{has_prefix, Case};
use crate::error::MagnetoError;
use crate::layout::{geometry, glue, AvLayout};
use crate::local::{self, BcCoeffs, BulkCoeffs};
use crate::material::{conductivity, const_real, const_vector, reluctivity};

/// A solved steady AV system.
pub struct AvSolution {
    /// Nodal DOFs (scalar potential), 0 without `Static Conductivity`.
    pub n_nodal: usize,
    /// Edge DOFs (vector potential).
    pub n_edges: usize,
    /// `[V (nodal slots, node order) | A (edge slots, first-use order)]`.
    pub values: Vec<f64>,
    /// Slot of each mesh node (`u32::MAX` when not in the system).
    pub node_slot: Vec<u32>,
    /// The edge numbering behind the edge slots.
    pub edges: crate::topology::Edges,
    /// `ComputeNorm` with `Norm Permutation = nNodes+1`: the RMS of the
    /// edge DOFs.
    pub norm: f64,
}

/// Keywords whose behaviour is not ported, refused when set in the solver
/// section (true logicals or any value).
pub(crate) const REFUSED_SOLVER_KEYS: &[&str] = &[
    "electrodynamics model",
    "darwin model",
    "eigen analysis",
    "use lagrange gauge",
    "use tree gauge",
    "helmholtz projection",
    "mass regularize epsilon",
    "current density name",
    "use nodal coilcurrent",
    "use elemental coilcurrent",
    "precomputed source term",
    "precomputed fixing term",
    "use piola transform",
    "quadratic approximation",
    "second kind basis",
    "calculate magnetic norm",
    "use magnetic norm",
    "nonlinear system torque tolerance",
    "linear system adaptive tolerance",
    "constant system",
    "constant bulk system",
    "linear system residual mode",
    "nonlinear system relaxation factor",
];

pub(crate) const REFUSED_BC_PREFIXES: &[&str] = &[
    "perfect electric conductor",
    "outflow bc",
    "thin line crossection area",
    "magnetic flux density",
    "mortar bc",
    "periodic bc",
    "normal-tangential",
];

fn check_solver(case: &Case, sp: &ParamMap) -> Result<(), MagnetoError> {
    for k in REFUSED_SOLVER_KEYS {
        match sp.get(k) {
            None | Some(elmer_core::Value::Logical(false)) => {}
            Some(_) => return Err(MagnetoError::unsupported(format!("WhitneyAVSolver option '{k}'"))),
        }
    }
    if sp.get("element").is_some() {
        return Err(MagnetoError::unsupported("an explicit 'Element' definition"));
    }
    let sim_type = case.simulation.string("simulation type").unwrap_or("steady state");
    if !sim_type.trim().to_ascii_lowercase().starts_with("steady") {
        return Err(MagnetoError::unsupported(format!("Simulation Type = {sim_type}")));
    }
    if sp.string("linear system solver").map(elmer_core::normalize_key).as_deref() == Some("direct") {
        // Elmer then defaults to a tree gauge (WhitneyAVSolver.F90:534).
        return Err(MagnetoError::unsupported("the tree gauge of a direct solve"));
    }
    if sp.string("linear system solver").map(elmer_core::normalize_key).as_deref() == Some("block") {
        return Err(MagnetoError::unsupported("Linear System Solver = block"));
    }
    if case.components.iter().any(Option::is_some) {
        return Err(MagnetoError::unsupported("Component sections (coils)"));
    }
    for bf in case.body_forces.iter().flatten() {
        for k in ["angular velocity", "lorentz velocity", "local lagrange gauge", "current density multiplier"] {
            if has_prefix(bf, k) {
                return Err(MagnetoError::unsupported(format!("body force '{k}'")));
            }
        }
    }
    for m in case.materials.iter().flatten() {
        if has_prefix(m, "laminate stack model") || has_prefix(m, "magnetization") {
            return Err(MagnetoError::unsupported("material laminate stack or magnetization"));
        }
    }
    for bc in case.bcs.iter().flatten() {
        for k in REFUSED_BC_PREFIXES {
            if has_prefix(bc, k) {
                return Err(MagnetoError::unsupported(format!("boundary condition '{k}'")));
            }
        }
    }
    Ok(())
}

/// `Static Conductivity`, or implied by a boundary `Electric Current
/// Density` (WhitneyAVSolver_Init0): the scalar potential is solved too.
pub(crate) fn static_conductivity(case: &Case, sp: &ParamMap) -> bool {
    sp.logical("static conductivity")
        .unwrap_or_else(|| case.bcs.iter().flatten().any(|bc| has_prefix(bc, "electric current density")))
}

/// `Fix Input Current Density`; in steady state it defaults to whether any
/// body force gives a `Current Density` (WhitneyAVSolver.F90:608).
fn jfix(case: &Case, sp: &ParamMap) -> bool {
    sp.logical("fix input current density")
        .unwrap_or_else(|| case.body_forces.iter().flatten().any(|bf| has_prefix(bf, "current density")))
}

/// Steady-state `WhitneyAVSolver` for solver section `solver_id`.
pub fn solve_steady(case: &Case, solver_id: usize) -> Result<AvSolution, MagnetoError> {
    let sp = case
        .solver(solver_id)
        .ok_or_else(|| MagnetoError::Input(format!("no Solver {solver_id}")))?;
    check_solver(case, sp)?;
    if jfix(case, sp) {
        return Err(MagnetoError::unsupported("Fix Input Current Density (JFixPotentialSolver)"));
    }
    if case.mesh.dim != 3 {
        return Err(MagnetoError::Input("Solver requires 3D mesh!".into()));
    }
    let with_v = static_conductivity(case, sp);
    let var = elmer_core::normalize_key(sp.string("variable").unwrap_or("AV"));
    let var = var.split_whitespace().next().unwrap_or("av").to_string();

    let lay = AvLayout::build(case, solver_id, with_v)?;
    let (n, n_nodal) = (lay.n(), lay.n_nodal);
    let mut a = lay.matrix::<f64>();
    let mut rhs = vec![0.0; n];

    // Bulk assembly.
    for (k, &e) in lay.active.iter().enumerate() {
        let el = &case.mesh.bulk[e];
        let nn = el.nodes.len();
        let material = case
            .material(el.body)
            .ok_or_else(|| MagnetoError::Input(format!("body {} has no material", el.body)))?;
        let nu = reluctivity(material, &case.constants)?;
        let sigma = conductivity(material)?;
        let mut load = [0.0; 7];
        if let Some(bf) = case.body_force(el.body) {
            if let Some(j) = const_vector(bf, "current density")? {
                load[..3].copy_from_slice(&j);
            }
            if let Some(m) = const_vector(bf, "magnetization")? {
                load[3..6].copy_from_slice(&m);
            }
            if let Some(v) = const_real(bf, "electric potential")? {
                load[6] = v;
            }
        }
        let co = BulkCoeffs {
            nu: vec![nu; nn],
            sigma: vec![sigma; nn],
            load: vec![load; nn],
        };
        let (xyz, gindex) = geometry(case, &el.nodes);
        let np = if with_v { nn } else { 0 };
        let loc = local::bulk_matrix(el.code, &xyz, &gindex, np, &co)?;
        glue(&mut a, &mut rhs, &lay.elem_dofs[k], &loc.stiff, &loc.force);
    }

    // Boundary terms (Robin type in H, current density, air gap).
    let is_active = |p: Option<usize>| p.is_some_and(|p| case.solver_active_in(solver_id, case.mesh.bulk[p].body));
    for b in &case.mesh.boundary {
        let Some(bc) = case.bc(b.constraint) else { continue };
        if b.code / 100 == 1 {
            continue;
        }
        let h = const_vector(bc, "magnetic field strength")?;
        let transfer = const_real(bc, "magnetic transfer coefficient")?;
        let jn = const_real(bc, "electric current density")?;
        let etc = const_real(bc, "electric transfer coefficient")?;
        let gap = const_real(bc, "air gap length")?;
        let any = h.is_some() || transfer.is_some() || jn.is_some() || etc.is_some();
        if gap.is_none() && !any {
            continue;
        }
        if b.code / 100 == 2 {
            return Err(MagnetoError::unsupported("boundary conditions on line elements"));
        }
        if !is_active(b.parents[0]) && !is_active(b.parents[1]) {
            continue;
        }
        let eids = lay
            .edges
            .of_element(b.code, &b.nodes)
            .ok_or_else(|| MagnetoError::Input("boundary element outside the active elements".into()))?;
        let nn = b.nodes.len();
        let np = if with_v { nn } else { 0 };
        let (xyz, gindex) = geometry(case, &b.nodes);
        let loc = if let Some(gap) = gap {
            let mu = const_real(bc, "air gap relative permeability")?.unwrap_or(1.0);
            local::air_gap_matrix(b.code, &xyz, &gindex, np, &vec![gap; nn], &vec![mu; nn])?
        } else {
            let h = h.unwrap_or([0.0; 3]);
            let co = BcCoeffs {
                transfer: vec![transfer.unwrap_or(0.0); nn],
                load: vec![[h[0], h[1], h[2], jn.unwrap_or(0.0), etc.unwrap_or(0.0)]; nn],
            };
            local::bc_matrix(b.code, &xyz, &gindex, np, &co)?
        };
        glue(&mut a, &mut rhs, &lay.dofs_of(&b.nodes, &eids), &loc.stiff, &loc.force);
    }

    // DefaultDirichletBCs: nodal `V` values and `V {e}` edge values.
    let mut fixed: Vec<Option<f64>> = vec![None; n];
    let edge_key = format!("{var} {{e}}");
    for b in &case.mesh.boundary {
        let Some(bc) = case.bc(b.constraint) else { continue };
        if bc.get(&format!("{var} condition")).is_some() {
            return Err(MagnetoError::unsupported(format!("'{var} Condition'")));
        }
        if with_v {
            if let Some(v) = const_real(bc, &var)? {
                for &i in &b.nodes {
                    let s = lay.node_slot[i as usize];
                    if s != u32::MAX {
                        fixed[s as usize] = Some(v);
                    }
                }
            }
        }
        if has_prefix(bc, &edge_key) {
            for (k, _) in bc.0.iter().filter(|(k, _)| k.starts_with(&edge_key)) {
                if const_real(bc, k)? != Some(0.0) {
                    return Err(MagnetoError::unsupported(format!("a nonzero edge Dirichlet value '{k}'")));
                }
            }
            if !is_active(b.parents[0]) && !is_active(b.parents[1]) {
                continue;
            }
            if let Some(eids) = lay.edges.of_element(b.code, &b.nodes) {
                for e in eids {
                    fixed[n_nodal + e as usize] = Some(0.0);
                }
            }
        }
    }

    // ConstrainUnused: rows with a zero diagonal (e.g. the potential
    // outside conductors) become identity rows.
    let diag = a.diag.clone().expect("diagonal");
    for r in 0..n {
        if a.values[diag[r]] == 0.0 {
            a.values[a.row_ptr[r]..a.row_ptr[r + 1]].fill(0.0);
            a.values[diag[r]] = 1.0;
            rhs[r] = 0.0;
            fixed[r].get_or_insert(0.0);
        }
    }

    // EnforceDirichletConditions, done at solve time in Elmer.
    let list: Vec<(u32, f64)> = fixed
        .iter()
        .enumerate()
        .filter_map(|(r, v)| v.map(|v| (r as u32, v)))
        .collect();
    let mode = elmer_solver::assembly::dirichlet_mode(sp)?;
    elmer_assembly::apply_dirichlet(&mut a, &mut rhs, &list, mode)?;

    let mut x = vec![0.0; n];
    for &(r, v) in &list {
        x[r as usize] = v;
    }
    solve_linear(sp, &mut a, &mut rhs, &mut x)?;

    let edge_vals = &x[n_nodal..];
    let norm = if edge_vals.is_empty() {
        0.0
    } else {
        let s = edge_vals.iter().fold(0.0, |s, v| s + v * v);
        (s / edge_vals.len() as f64).sqrt()
    };
    Ok(AvSolution {
        n_nodal,
        n_edges: lay.edges.n,
        values: x,
        node_slot: lay.node_slot,
        edges: lay.edges,
        norm,
    })
}

/// `SolveLinearSystem` (SolveCore.F90:1186) for a real system: the zero
/// right-hand side shortcut, `Linear System Scaling` around the solve and
/// the solver the section selects.
pub fn solve_linear(sp: &ParamMap, a: &mut CsrMatrix<f64>, b: &mut [f64], x: &mut [f64]) -> Result<(), MagnetoError> {
    let log = NullLogger;
    let comm = SerialComm;
    let bnorm = b.iter().fold(0.0, |s, v| s + v * v).sqrt();
    if bnorm <= f64::MIN_POSITIVE && sp.logical("skip zero rhs test") != Some(true) {
        x.fill(0.0);
        return Ok(());
    }
    if sp.logical("linear system complex") == Some(true) {
        return Err(MagnetoError::unsupported("a complex linear system"));
    }
    let scale = sp.logical("linear system skip scaling") != Some(true) && sp.logical("linear system scaling").unwrap_or(true);
    let scaling = if scale {
        let method = if sp.logical("linear system row equilibration") == Some(true) {
            elmer_sparse::ScalingMethod::RowEquilibration
        } else {
            elmer_sparse::ScalingMethod::from_keyword(sp.string("linear system scaling method").unwrap_or("diagonal"))?
        };
        let opts = elmer_sparse::ScaleOptions {
            method,
            apply: true,
            rhs_scaling: bnorm != 0.0,
        };
        Some(elmer_sparse::scale_system(a, Some(&mut *b), Some(&mut *x), &opts))
    } else {
        None
    };
    let registry = elmer_solver::linsys::elmer_registry();
    let solved = match registry.select_for(sp, a, &log)? {
        LinearSolverChoice::Direct(mut d) => d.factor(a).and_then(|_| {
            let mut y = b.to_vec();
            d.solve_in_place(&mut y).map(|_| x.copy_from_slice(&y))
        }),
        LinearSolverChoice::Iterative { mut solver, precond } => {
            let m: Box<dyn elmer_core::Preconditioner<f64> + Send> = match precond {
                Some(f) => f(sp, a)?,
                None => Box::new(elmer_core::IdentityPrecond),
            };
            solver.solve(a, &*m, b, x, &comm).map(|_| ())
        }
    };
    if let Some(rs) = scaling {
        rs.unscale(a, Some(b), Some(x));
    }
    solved?;
    Ok(())
}
