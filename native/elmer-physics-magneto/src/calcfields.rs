//! `MagnetoDynamicsCalcFields` (CalcFields.F90:593) for a real (steady)
//! `WhitneyAVSolver` potential on the lowest-order Whitney basis: the
//! nodal fields `Magnetic Flux Density` (`B = curl A`, on by default),
//! `Magnetic Field Strength` (`H = nu B`) and `Magnetic Vector Potential`,
//! each an L2 projection onto the nodal basis of the active elements
//! (`GlobalSol`: one mass-matrix solve per component with the section's
//! linear solver, from a zero guess).
//!
//! The solver's own variable (`hr_dummy`) holds the last component solved,
//! so the norm Elmer reports is that component's (`Show Norm Index` picks
//! another). Elemental (`E`) fields, energies and the other `Calculate *`
//! options are not ported: the former are skipped (they do not change
//! the solver norm), the latter refused.

use elmer_core::{ParamMap, Params, Value};

use crate::case::{has_prefix, Case};
use crate::error::MagnetoError;
use crate::layout::{geometry, glue};
use crate::local::{edge_rule, eval_point};
use crate::material::reluctivity;
use crate::topology::Edges;

/// The potential the fields are computed from: edge `e`'s value is
/// `values[(n_nodal + e) * dofs]` (real part).
pub struct Potential<'a> {
    pub edges: &'a Edges,
    pub n_nodal: usize,
    pub dofs: usize,
    pub values: &'a [f64],
}

/// One projected nodal field.
pub struct NodalField {
    /// Elmer's variable name (`magnetic flux density`, ...).
    pub name: &'static str,
    /// `values[slot * 3 + k]`.
    pub values: Vec<f64>,
    /// `ComputeNorm` of each component solve.
    pub norms: [f64; 3],
}

pub struct CalcFieldsResult {
    /// Slot of each mesh node (`u32::MAX` outside the active elements).
    pub node_slot: Vec<u32>,
    pub n: usize,
    pub fields: Vec<NodalField>,
    /// The solver norm: the last component solved, or `Show Norm Index`.
    pub norm: f64,
}

/// Options this port refuses rather than ignores.
const REFUSED: &[&str] = &[
    "calculate current density",
    "calculate electric field",
    "calculate jxb",
    "calculate maxwell stress",
    "calculate joule heating",
    "calculate harmonic loss",
    "calculate homogenization loss",
    "calculate winding voltage",
    "calculate relative permeability",
    "calculate electric scalar potential",
    "calculate nodal heating",
    "calculate nodal forces",
    "calculate flux linkage",
    "discontinuous galerkin",
    "average within materials",
    "constant mass matrix",
    "2d result fields",
    "precomputed electric potential",
    "impose circuit current",
    "vtu style",
];

/// The solver section whose variable is the potential: `Potential
/// Variable`'s solver, else the last `WhitneyAVSolver` /
/// `WhitneyAVHarmonicSolver` (`MagnetoDynamicsCalcFields_Init0`).
pub fn target_solver(case: &Case, sp: &ParamMap) -> Result<(usize, bool), MagnetoError> {
    let proc_entry = |s: &ParamMap| match s.get("procedure") {
        Some(Value::Dependent(d)) => match &d.kind {
            elmer_core::DependentKind::Procedure { entry, .. } => entry.to_ascii_lowercase(),
            _ => String::new(),
        },
        Some(Value::Str(s)) => s.to_ascii_lowercase(),
        _ => String::new(),
    };
    let harmonic = |s: &ParamMap| proc_entry(s).contains("whitneyavharmonicsolver");
    if let Some(pname) = sp.string("potential variable") {
        let pname = elmer_core::normalize_key(pname);
        for (i, s) in case.solvers.iter().enumerate() {
            let Some(s) = s else { continue };
            let var = elmer_core::normalize_key(s.string("variable").unwrap_or(""));
            if var.split('[').next().map(str::trim) == Some(pname.as_str()) {
                return Ok((i + 1, harmonic(s)));
            }
        }
        return Err(MagnetoError::Input(format!("Could not find solver for variable: {pname}")));
    }
    for (i, s) in case.solvers.iter().enumerate().rev() {
        let Some(s) = s else { continue };
        let e = proc_entry(s);
        if e.contains("whitneyavharmonicsolver") {
            return Ok((i + 1, true));
        }
        if e.contains("whitneyavsolver") {
            return Ok((i + 1, false));
        }
    }
    Err(MagnetoError::Input("Could not determine target variable type (real or complex)".into()))
}

/// `MagnetoDynamicsCalcFields` for solver section `solver_id` and a real
/// potential.
pub fn calc_fields(case: &Case, solver_id: usize, pot: &Potential<'_>) -> Result<CalcFieldsResult, MagnetoError> {
    let sp = case
        .solver(solver_id)
        .ok_or_else(|| MagnetoError::Input(format!("no Solver {solver_id}")))?;
    for k in REFUSED {
        match sp.get(k) {
            None | Some(Value::Logical(false)) => {}
            Some(_) => return Err(MagnetoError::unsupported(format!("MagnetoDynamicsCalcFields option '{k}'"))),
        }
    }
    if pot.dofs != 1 {
        return Err(MagnetoError::unsupported("CalcFields of a complex (harmonic) potential"));
    }
    if case.mesh.dim != 3 {
        return Err(MagnetoError::unsupported("CalcFields on a 2D mesh"));
    }
    let sim_type = case.simulation.string("simulation type").unwrap_or("steady state");
    if !sim_type.trim().to_ascii_lowercase().starts_with("steady") {
        return Err(MagnetoError::unsupported(format!("CalcFields with Simulation Type = {sim_type}")));
    }
    for bf in case.body_forces.iter().flatten() {
        for k in ["angular velocity", "lorentz velocity", "magnetization"] {
            if has_prefix(bf, k) {
                return Err(MagnetoError::unsupported(format!("CalcFields with body force '{k}'")));
            }
        }
    }
    let nodal = !sp.logical("skip nodal fields").unwrap_or(false)
        && sp.logical("calculate nodal fields").unwrap_or(true);
    let want_b = sp.logical("calculate magnetic flux density").unwrap_or(true);
    let want_vp = sp.logical("calculate magnetic vector potential").unwrap_or(false);
    let want_h = sp.logical("calculate magnetic field strength").unwrap_or(false);

    // The nodal Perm of the CalcFields solver: nodes of its active elements.
    let active: Vec<usize> = (0..case.mesh.bulk.len())
        .filter(|&e| case.solver_active_in(solver_id, case.mesh.bulk[e].body))
        .collect();
    let mut node_slot = vec![u32::MAX; case.mesh.nodes.len()];
    for &e in &active {
        for &i in &case.mesh.bulk[e].nodes {
            node_slot[i as usize] = 0;
        }
    }
    let mut n = 0u32;
    for s in node_slot.iter_mut().filter(|s| **s == 0) {
        *s = n;
        n += 1;
    }
    let n = n as usize;
    if !nodal || !(want_b || want_vp || want_h) {
        return Ok(CalcFieldsResult {
            node_slot,
            n,
            fields: Vec::new(),
            norm: 0.0,
        });
    }

    // Element mass matrices and right-hand sides, columns in Elmer's order
    // (B, H, A), each three components.
    let ncol = 3 * (want_b as usize + want_h as usize + want_vp as usize);
    let elem_dofs: Vec<Vec<u32>> = active
        .iter()
        .map(|&e| case.mesh.bulk[e].nodes.iter().map(|&i| node_slot[i as usize]).collect())
        .collect();
    let mut mass = crate::layout::nodal_matrix(n, &elem_dofs);
    let mut rhs = vec![vec![0.0; n]; ncol];
    for (k, &e) in active.iter().enumerate() {
        let el = &case.mesh.bulk[e];
        let nn = el.nodes.len();
        let eids = pot
            .edges
            .of_element(el.code, &el.nodes)
            .ok_or_else(|| MagnetoError::Input("CalcFields element outside the potential's elements".into()))?;
        let a: Vec<f64> = eids.iter().map(|&e| pot.values[pot.n_nodal + e as usize]).collect();
        let nu = if want_h {
            let material = case
                .material(el.body)
                .ok_or_else(|| MagnetoError::Input(format!("body {} has no material", el.body)))?;
            reluctivity(material, &case.constants)?
        } else {
            0.0
        };
        let (xyz, gindex) = geometry(case, &el.nodes);
        let rule = edge_rule(el.code)?;
        let mut status = elmer_fe::metric::MetricStatus::default();
        let mut m = vec![0.0; nn * nn];
        let mut f = vec![vec![0.0; nn]; ncol];
        for t in 0..rule.len() {
            let pt = eval_point(el.code, &xyz, &gindex, [rule.u[t], rule.v[t], rule.w[t]], &mut status)?;
            let s = rule.s[t] * pt.det_j;
            let basis = &pt.basis[..nn];
            // B = MATMUL(SOL(np+1:nd), RotWBasis), A likewise with WBasis.
            let mut b = [0.0; 3];
            let mut vp = [0.0; 3];
            for c in 0..3 {
                for (j, aj) in a.iter().enumerate() {
                    b[c] += aj * pt.curl[j][c];
                    vp[c] += aj * pt.w[j][c];
                }
            }
            let r_ip = basis.iter().fold(0.0, |acc, bp| acc + bp * nu);
            // MATMUL(Nu, B) with Nu = diag(R_ip).
            let h = [r_ip * b[0], r_ip * b[1], r_ip * b[2]];
            for p in 0..nn {
                for q in 0..nn {
                    m[p * nn + q] += s * basis[p] * basis[q];
                }
                let mut col = 0;
                for (on, v) in [(want_b, b), (want_h, h), (want_vp, vp)] {
                    if on {
                        for c in 0..3 {
                            f[col + c][p] += s * v[c] * basis[p];
                        }
                        col += 3;
                    }
                }
            }
        }
        let zero = vec![0.0; nn];
        glue(&mut mass.0, &mut mass.1, &elem_dofs[k], &m, &zero);
        for (c, fc) in f.iter().enumerate() {
            for (p, &r) in elem_dofs[k].iter().enumerate() {
                rhs[c][r as usize] += fc[p];
            }
        }
    }

    // GlobalSol: one solve per component, from a zero guess.
    let norm_index = sp.integer("show norm index").filter(|&i| i > 0).map(|i| i as usize);
    let mut fields = Vec::new();
    let mut col = 0;
    let mut last_norm = 0.0;
    let mut saved = 0.0;
    for (on, name) in [
        (want_b, "magnetic flux density"),
        (want_h, "magnetic field strength"),
        (want_vp, "magnetic vector potential"),
    ] {
        if !on {
            continue;
        }
        let mut values = vec![0.0; 3 * n];
        let mut norms = [0.0; 3];
        for c in 0..3 {
            let mut a = mass.0.clone();
            let mut b = rhs[col + c].clone();
            let mut x = vec![0.0; n];
            crate::steady::solve_linear(sp, &mut a, &mut b, &mut x)?;
            let norm = (x.iter().fold(0.0, |acc, v| acc + v * v) / n as f64).sqrt();
            norms[c] = norm;
            last_norm = norm;
            if norm_index == Some(col + c + 1) {
                saved = norm;
            }
            for (i, v) in x.iter().enumerate() {
                values[3 * i + c] = *v;
            }
        }
        fields.push(NodalField { name, values, norms });
        col += 3;
    }
    Ok(CalcFieldsResult {
        node_slot,
        n,
        fields,
        norm: if norm_index.is_some() { saved } else { last_norm },
    })
}
