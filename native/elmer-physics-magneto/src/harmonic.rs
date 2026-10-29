//! `WhitneyAVHarmonicSolver` (WhitneyAVHarmonicSolver.F90:194): the
//! time-harmonic AV equations
//!
//! ```text
//!   curl(nu curl A) + i omega sigma A + sigma grad V = J^s + curl M^s - sigma grad V^s
//!   -div(sigma (i omega A + grad V)) = 0
//! ```
//!
//! on the lowest-order Whitney basis (`Element = "n:1 e:1"`), for real
//! material and source values, solved as a native complex system the way
//! `SolveLinearSystem` treats `Linear System Complex`. The norm is
//! `ComputeNorm` over all real components of the interleaved variable
//! (`[re, im]` per slot), as this solver sets no `Norm Permutation`.

use elmer_core::{CsrMatrix, IterativeSolver, NullLogger, ParamMap, Params, Preconditioner, SerialComm, Value};
use num_complex::Complex64;

use crate::case::{has_prefix, Case};
use crate::error::MagnetoError;
use crate::layout::{geometry, glue, AvLayout};
use crate::local::{self, BcCoeffs, BulkCoeffs};
use crate::material::{conductivity, const_real, const_vector, reluctivity, PI};
use crate::steady::{REFUSED_BC_PREFIXES, REFUSED_SOLVER_KEYS};

/// A solved harmonic AV system.
pub struct HarmonicSolution {
    pub n_nodal: usize,
    pub n_edges: usize,
    /// `[V (nodal slots) | A (edge slots)]`, complex.
    pub values: Vec<Complex64>,
    pub node_slot: Vec<u32>,
    /// The edge numbering behind the edge slots.
    pub edges: crate::topology::Edges,
    /// The angular frequency used.
    pub omega: f64,
    /// Iterations and final residual of the linear solve.
    pub stats: elmer_core::SolveStats,
    /// `ComputeNorm` over the interleaved real values.
    pub norm: f64,
}

/// `Variable = P[P re:1 P im:1]`: the two component names.
fn components(sp: &ParamMap) -> Result<(String, String), MagnetoError> {
    let spec = sp.string("variable").unwrap_or("AV[AV re:1 AV im:1]");
    let (_, inner) = spec
        .split_once('[')
        .ok_or_else(|| MagnetoError::Input(format!("harmonic variable needs two components: {spec}")))?;
    let inner = inner.trim_end().trim_end_matches(']');
    let comps: Vec<String> = inner
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .split(':')
        .filter_map(|part| {
            // "P re:1 P im:1" splits into "P re", "1 P im", "1".
            let words: Vec<&str> = part.split_whitespace().filter(|w| w.parse::<i32>().is_err()).collect();
            (!words.is_empty()).then(|| elmer_core::normalize_key(&words.join(" ")))
        })
        .collect();
    match comps.as_slice() {
        [re, im] => Ok((re.clone(), im.clone())),
        _ => Err(MagnetoError::Input(format!(
            "Variable is not properly defined for time harmonic AV solver: {spec}"
        ))),
    }
}

/// `ListGetAngularFrequency` (Lists.F90:10543) for an element of `body`:
/// `Frequency` (times 2 pi) or `Angular Frequency` from the body's
/// equation, then its material, the simulation and the solver section.
fn angular_frequency(case: &Case, sp: &ParamMap, body: i32) -> Result<Option<f64>, MagnetoError> {
    let lists = [case.equation(body), case.material(body), Some(&case.simulation), Some(sp)];
    for p in lists.into_iter().flatten() {
        if let Some(f) = const_real(p, "frequency")? {
            return Ok(Some(2.0 * PI * f));
        }
        if let Some(w) = const_real(p, "angular frequency")? {
            return Ok(Some(w));
        }
    }
    Ok(None)
}

fn check(case: &Case, sp: &ParamMap) -> Result<(), MagnetoError> {
    for k in REFUSED_SOLVER_KEYS.iter().chain(&["h-b curve variable"]) {
        match sp.get(k) {
            None | Some(Value::Logical(false)) => {}
            Some(_) => return Err(MagnetoError::unsupported(format!("WhitneyAVHarmonicSolver option '{k}'"))),
        }
    }
    if sp.get("element").is_some() {
        return Err(MagnetoError::unsupported("an explicit 'Element' definition"));
    }
    if sp.logical("linear system complex") == Some(false) {
        return Err(MagnetoError::unsupported("Linear System Complex = False"));
    }
    if sp.string("linear system solver").map(elmer_core::normalize_key).as_deref() != Some("iterative") {
        // A direct solve defaults to the tree gauge (WhitneyAVHarmonicSolver.F90:679).
        return Err(MagnetoError::unsupported("a non-iterative linear solver"));
    }
    if case.components.iter().any(Option::is_some) {
        return Err(MagnetoError::unsupported("Component sections (coils)"));
    }
    let im_key = |p: &ParamMap| p.0.keys().any(|k| k.ends_with(" im") || k.contains(" im "));
    for bf in case.body_forces.iter().flatten() {
        for k in ["angular velocity", "lorentz velocity", "current density multiplier"] {
            if has_prefix(bf, k) {
                return Err(MagnetoError::unsupported(format!("body force '{k}'")));
            }
        }
        if im_key(bf) {
            return Err(MagnetoError::unsupported("imaginary body force values"));
        }
    }
    for m in case.materials.iter().flatten() {
        for k in ["laminate stack model", "magnetization", "h-b curve", "electric permittivity"] {
            if has_prefix(m, k) {
                return Err(MagnetoError::unsupported(format!("material '{k}'")));
            }
        }
        if im_key(m) {
            return Err(MagnetoError::unsupported("imaginary material values"));
        }
    }
    for bc in case.bcs.iter().flatten() {
        for k in REFUSED_BC_PREFIXES
            .iter()
            .chain(&["layer electric conductivity", "component"])
        {
            if has_prefix(bc, k) {
                return Err(MagnetoError::unsupported(format!("boundary condition '{k}'")));
            }
        }
        if bc.0.keys().any(|k| {
            (k.starts_with("magnetic") || k.starts_with("electric")) && (k.ends_with(" im") || k.contains(" im "))
        }) {
            return Err(MagnetoError::unsupported("imaginary boundary values"));
        }
    }
    Ok(())
}

/// `WhitneyAVHarmonicSolver` for solver section `solver_id`, one
/// nonlinear iteration (exact for linear materials).
pub fn solve_harmonic(case: &Case, solver_id: usize) -> Result<HarmonicSolution, MagnetoError> {
    let sp = case
        .solver(solver_id)
        .ok_or_else(|| MagnetoError::Input(format!("no Solver {solver_id}")))?;
    check(case, sp)?;
    let jfix = sp
        .logical("fix input current density")
        .unwrap_or_else(|| case.body_forces.iter().flatten().any(|bf| has_prefix(bf, "current density")));
    if jfix {
        return Err(MagnetoError::unsupported("Fix Input Current Density (JFixPotentialSolver)"));
    }
    if case.mesh.dim != 3 {
        return Err(MagnetoError::Input("Solver requires 3D mesh!".into()));
    }
    let (name_re, name_im) = components(sp)?;

    let lay = AvLayout::build(case, solver_id, true)?;
    let (n, n_nodal) = (lay.n(), lay.n_nodal);
    let mut a = lay.matrix::<Complex64>();
    let mut rhs = vec![Complex64::new(0.0, 0.0); n];
    let mut omega_used = None;

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
        let omega = angular_frequency(case, sp, el.body)?
            .ok_or_else(|| MagnetoError::Input("Harmonic solution requires frequency!".into()))?;
        // The module-level Omega the boundary kernels read is the last
        // bulk element's.
        omega_used = Some(omega);
        let co = BulkCoeffs {
            nu: vec![nu; nn],
            sigma: vec![sigma; nn],
            load: vec![load; nn],
        };
        let (xyz, gindex) = geometry(case, &el.nodes);
        let (loc, damp) = local::harmonic_bulk_matrix(el.code, &xyz, &gindex, nn, &co)?;
        // STIFF = -omega^2 MASS + i omega DAMP + STIFF with MASS = 0.
        let stiff: Vec<Complex64> = loc
            .stiff
            .iter()
            .zip(&damp)
            .map(|(&k, &d)| Complex64::new(k, omega * d))
            .collect();
        let force: Vec<Complex64> = loc.force.iter().map(|&f| Complex64::new(f, 0.0)).collect();
        glue(&mut a, &mut rhs, &lay.elem_dofs[k], &stiff, &force);
    }

    let is_active = |p: Option<usize>| p.is_some_and(|p| case.solver_active_in(solver_id, case.mesh.bulk[p].body));
    for b in &case.mesh.boundary {
        let Some(bc) = case.bc(b.constraint) else { continue };
        if b.code / 100 == 1 {
            continue;
        }
        let h = const_vector(bc, "magnetic field strength")?;
        let transfer = const_real(bc, "magnetic transfer coefficient")?;
        let jn = match const_real(bc, "electric current density")? {
            Some(j) => Some(j),
            None => const_real(bc, "electric flux")?,
        };
        let etc = const_real(bc, "electric transfer coefficient")?;
        let gap = const_real(bc, "air gap length")?;
        let sheet = const_real(bc, "thin sheet thickness")?;
        if gap.is_none() && sheet.is_none() && h.is_none() && transfer.is_none() && jn.is_none() && etc.is_none() {
            // LocalMatrixBC with all-zero coefficients adds nothing.
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
        let (xyz, gindex) = geometry(case, &b.nodes);
        let (loc, damp) = if let Some(gap) = gap {
            let mu = const_real(bc, "air gap relative permeability")?.unwrap_or(1.0);
            (local::air_gap_matrix(b.code, &xyz, &gindex, nn, &vec![gap; nn], &vec![mu; nn])?, None)
        } else if let Some(d) = sheet {
            let mu = const_real(bc, "thin sheet relative permeability")?.unwrap_or(1.0);
            let sigma = const_real(bc, "thin sheet electric conductivity")?.unwrap_or(1.0);
            let (loc, damp) =
                local::thin_sheet_matrix(b.code, &xyz, &gindex, nn, &vec![d; nn], &vec![mu; nn], &vec![sigma; nn])?;
            (loc, Some(damp))
        } else {
            let h = h.unwrap_or([0.0; 3]);
            let co = BcCoeffs {
                transfer: vec![transfer.unwrap_or(0.0); nn],
                load: vec![[h[0], h[1], h[2], jn.unwrap_or(0.0), etc.unwrap_or(0.0)]; nn],
            };
            (local::bc_matrix(b.code, &xyz, &gindex, nn, &co)?, None)
        };
        let omega = omega_used.unwrap_or(0.0);
        let stiff: Vec<Complex64> = match damp {
            // STIFF = -omega^2 MASS + i omega DAMP + STIFF, MASS = 0.
            Some(d) => loc.stiff.iter().zip(&d).map(|(&k, &d)| Complex64::new(k, omega * d)).collect(),
            None => loc.stiff.iter().map(|&v| Complex64::new(v, 0.0)).collect(),
        };
        let force: Vec<Complex64> = loc.force.iter().map(|&v| Complex64::new(v, 0.0)).collect();
        glue(&mut a, &mut rhs, &lay.dofs_of(&b.nodes, &eids), &stiff, &force);
    }

    // DefaultDirichletBCs, per component: nodal values and zero edge values.
    let mut fixed_re: Vec<Option<f64>> = vec![None; n];
    let mut fixed_im: Vec<Option<f64>> = vec![None; n];
    for b in &case.mesh.boundary {
        let Some(bc) = case.bc(b.constraint) else { continue };
        for (name, fixed) in [(&name_re, &mut fixed_re), (&name_im, &mut fixed_im)] {
            if bc.get(&format!("{name} condition")).is_some() {
                return Err(MagnetoError::unsupported(format!("'{name} Condition'")));
            }
            if let Some(v) = const_real(bc, name)? {
                for &i in &b.nodes {
                    let s = lay.node_slot[i as usize];
                    if s != u32::MAX {
                        fixed[s as usize] = Some(v);
                    }
                }
            }
            let edge_key = format!("{name} {{e}}");
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
    }
    let mut fixed: Vec<Option<Complex64>> = Vec::with_capacity(n);
    for r in 0..n {
        fixed.push(match (fixed_re[r], fixed_im[r]) {
            (None, None) => None,
            (Some(re), Some(im)) => Some(Complex64::new(re, im)),
            _ => {
                return Err(MagnetoError::unsupported(
                    "a Dirichlet condition on only one of the real and imaginary parts",
                ))
            }
        });
    }

    // ConstrainUnused (WhitneyAVHarmonicSolver.F90:795).
    let diag = a.diag.clone().expect("diagonal");
    for r in 0..n {
        if a.values[diag[r]].norm() == 0.0 {
            a.values[a.row_ptr[r]..a.row_ptr[r + 1]].fill(Complex64::new(0.0, 0.0));
            a.values[diag[r]] = Complex64::new(1.0, 0.0);
            rhs[r] = Complex64::new(0.0, 0.0);
            fixed[r].get_or_insert(Complex64::new(0.0, 0.0));
        }
    }
    enforce_dirichlet(&mut a, &mut rhs, &fixed, sp);

    let mut x: Vec<Complex64> = fixed.iter().map(|v| v.unwrap_or_default()).collect();
    let stats = solve_complex(sp, &mut a, &mut rhs, &mut x)?;

    let mut s = 0.0;
    for v in &x {
        s += v.re * v.re;
        s += v.im * v.im;
    }
    let norm = (s / (2 * n) as f64).sqrt();
    Ok(HarmonicSolution {
        n_nodal,
        n_edges: lay.edges.n,
        values: x,
        node_slot: lay.node_slot,
        edges: lay.edges,
        omega: omega_used.unwrap_or(0.0),
        stats,
        norm,
    })
}

/// `EnforceDirichletConditions` (BoundaryConditionUtils.F90:3269) on the
/// complex system: with `Linear System Symmetric`, move the constrained
/// columns of the other rows to the right-hand side; then each constrained
/// row becomes `s x = s v`, `s = |a_kk|` under the default diagonal
/// scaling (Elmer scales by the modulus of the complex diagonal), 1
/// without it.
fn enforce_dirichlet(a: &mut CsrMatrix<Complex64>, b: &mut [Complex64], fixed: &[Option<Complex64>], sp: &ParamMap) {
    let n = a.nrows;
    if fixed.iter().all(Option::is_none) {
        return;
    }
    let diag = a.diag.clone().expect("diagonal");
    let scaled = sp
        .logical("linear system dirichlet scaling")
        .or_else(|| sp.logical("linear system scaling"))
        .unwrap_or(true)
        && sp.string("linear system scaling method").map(elmer_core::normalize_key).as_deref() != Some("none");
    let scale: Vec<f64> = if scaled {
        <Complex64 as elmer_sparse::ElmerScalar>::diagonal_scaling(a)
    } else {
        vec![1.0; n]
    };
    if sp.logical("linear system symmetric") == Some(true) {
        for r in 0..n {
            if fixed[r].is_some() {
                continue;
            }
            for k in a.row_ptr[r]..a.row_ptr[r + 1] {
                if let Some(v) = fixed[a.col_idx[k] as usize] {
                    b[r] -= a.values[k] * v;
                    a.values[k] = Complex64::new(0.0, 0.0);
                }
            }
        }
    }
    for k in 0..n {
        let Some(v) = fixed[k] else { continue };
        let s = if scaled {
            let d = if scale[k].abs() <= f64::MIN_POSITIVE { 1.0 } else { scale[k] };
            1.0 / (d * d)
        } else {
            1.0
        };
        a.values[a.row_ptr[k]..a.row_ptr[k + 1]].fill(Complex64::new(0.0, 0.0));
        a.values[diag[k]] = Complex64::new(s, 0.0);
        b[k] = v * s;
    }
}

/// `SolveLinearSystem` for `Linear System Complex`: diagonal scaling of the
/// complex system, then the iterative method and preconditioner of the
/// section (elmer-solver's complex path: `BiCGStabl`, `GCR` and `IDRS` from
/// `elmer-itermethods`, the HUTI methods otherwise; `none`, `diagonal`,
/// `ILUn` and `ILUT` preconditioning).
fn solve_complex(
    sp: &ParamMap,
    a: &mut CsrMatrix<Complex64>,
    b: &mut [Complex64],
    x: &mut [Complex64],
) -> Result<elmer_core::SolveStats, MagnetoError> {
    let log = NullLogger;
    let bnorm = b.iter().fold(0.0, |s, v| s + v.re * v.re + v.im * v.im).sqrt();
    if bnorm <= f64::MIN_POSITIVE && sp.logical("skip zero rhs test") != Some(true) {
        x.fill(Complex64::new(0.0, 0.0));
        return Ok(elmer_core::SolveStats {
            iterations: 0,
            residual: 0.0,
            status: elmer_core::SolveStatus::Converged,
        });
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
    let method = sp
        .string("linear system iterative method")
        .map(elmer_core::normalize_key)
        .unwrap_or_else(|| "bicgstab".into());
    let mut solver: Box<dyn IterativeSolver<Complex64>> = match method.as_str() {
        "bicgstabl" | "gcr" | "idrs" => {
            let cfg = elmer_itermethods::IterMethodsConfig::from_params(sp)?;
            match method.as_str() {
                "bicgstabl" => Box::new(elmer_itermethods::BiCgStabL::<Complex64>::new(&cfg)),
                "gcr" => Box::new(elmer_itermethods::Gcr::<Complex64>::new(&cfg)),
                _ => Box::new(elmer_itermethods::Idrs::<Complex64>::new(&cfg)),
            }
        }
        "jacobi" | "richardson" | "sgs" | "bicgstab2" => {
            return Err(MagnetoError::unsupported(format!("complex iterative method {method}")))
        }
        _ => Box::new(elmer_krylov::Huti::<Complex64>::new(
            elmer_krylov::HutiConfig::from_params_with(sp, true, &log)?,
        )),
    };
    let pre = elmer_core::normalize_key(sp.string("linear system preconditioning").unwrap_or("none"));
    let precond: Box<dyn Preconditioner<Complex64>> = match pre.as_str() {
        "none" => Box::new(elmer_core::IdentityPrecond),
        "diagonal" => Box::new(elmer_sparse::Ilu::diagonal(a)),
        p if p.starts_with("ilu") && p != "ilut" && sp.get("linear system ilu order").is_none() => {
            let level = p.as_bytes().get(3).map_or(0, |&c| i32::from(c) - i32::from(b'0'));
            Box::new(elmer_sparse::Ilu::ilu_n(a, level.clamp(0, 9) as usize)?)
        }
        other => return Err(MagnetoError::unsupported(format!("complex preconditioner {other}"))),
    };
    let solved = solver.solve(a, &*precond, b, x, &SerialComm);
    if let Some(rs) = scaling {
        rs.unscale(a, Some(b), Some(x));
    }
    Ok(solved?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn component_names_of_the_variable_spec() {
        let mut p = ParamMap::new();
        p.insert("Variable", Value::Str("P[P re:1 P im:1]".into()));
        assert_eq!(components(&p).unwrap(), ("p re".into(), "p im".into()));
        assert_eq!(components(&ParamMap::new()).unwrap(), ("av re".into(), "av im".into()));
        p.insert("Variable", Value::Str("P".into()));
        assert!(components(&p).is_err());
    }
}
