//! `Harmonic Analysis` (SolveCore.F90:1408, `SolveHarmonicSystem` at
//! SolveCore.F90:4650): the assembled stiffness `K`, mass `M` (and
//! damping `D`) of a solver become, at `omega = 2 pi Frequency`, the real
//! doubled system
//!
//! ```text
//! [ K - w^2 M   -w D      ] [Re u]   [RHS   ]
//! [ w D          K - w^2 M] [Im u] = [RHS_im]
//! ```
//!
//! interleaved per dof (the `CreateMatrix(..., 2 DOFs)` pattern). Its
//! Dirichlet rows are `<component>` and `<component> im` (a missing one of
//! the pair is zero); it is solved with the solver's own linear system
//! settings and the solution is kept as `EigenVectors(1, :)`; the
//! variable's values are left alone (`SolveSystem` returns before its
//! `ComputeChange`).
//!
//! Refused: frequency sweeps (`Harmonic System Values` > 1), a damping
//! matrix, `<component> im` Dirichlet values and parallel runs.

use num_complex::Complex64;

use elmer_core::{CsrMatrix, Params};

use crate::error::SolverError;
use crate::solver::SolverCtx;

const CALLER: &str = "SolveHarmonicSystem";

fn unsupported(what: impl Into<String>) -> SolverError {
    SolverError::Unsupported(format!("Harmonic Analysis: {}", what.into()))
}

/// Whether this visit's solve is a harmonic one (`HarmonicAnalysis` of
/// `SolveSystem`).
pub(crate) fn active(ctx: &SolverCtx<'_>) -> bool {
    ctx.nof_eigen() > 0 && ctx.params().logical("harmonic analysis") == Some(true) && ctx.params().logical("eigen analysis") != Some(true)
}

/// The doubled pattern of `g`: each entry becomes a 2x2 block.
fn doubled(g: &CsrMatrix<f64>, kmw: impl Fn(usize) -> f64) -> CsrMatrix<f64> {
    let n = g.nrows;
    let mut row_ptr = Vec::with_capacity(2 * n + 1);
    let mut col_idx = Vec::with_capacity(4 * g.values.len());
    let mut values = Vec::with_capacity(4 * g.values.len());
    let mut diag = Vec::with_capacity(2 * n);
    row_ptr.push(0);
    for k in 0..n {
        for part in 0..2 {
            for j in g.row_ptr[k]..g.row_ptr[k + 1] {
                let c = g.col_idx[j] as usize;
                let v = kmw(j);
                // (re, re) / (re, im) on the first row, (im, re) / (im, im)
                // on the second; no damping, so the off blocks are zero.
                for q in 0..2 {
                    let col = 2 * c + q;
                    if col == 2 * k + part {
                        diag.push(col_idx.len());
                    }
                    col_idx.push(col as u32);
                    values.push(if q == part { v } else { 0.0 });
                }
            }
            row_ptr.push(col_idx.len());
        }
    }
    CsrMatrix { nrows: 2 * n, ncols: 2 * n, row_ptr, col_idx, values, diag: (diag.len() == 2 * n).then_some(diag) }
}

pub(crate) fn solve(ctx: &mut SolverCtx<'_>) -> Result<f64, SolverError> {
    let (log, comm, ix) = (ctx.log, ctx.comm, ctx.slot);
    if comm.size() > 1 {
        return Err(unsupported("parallel runs"));
    }
    let p = ctx.params().clone();
    if p.integer("harmonic system values").is_some_and(|n| n > 1) {
        return Err(unsupported("frequency sweeps (Harmonic System Values > 1)"));
    }
    let omega = ctx
        .creal(&p, "angular frequency")?
        .or(ctx.creal(&p, "frequency")?.map(|f| 2.0 * std::f64::consts::PI * f))
        .ok_or_else(|| SolverError::fatal(CALLER, "> Frequency < must be given for harmonic analysis."))?;
    // `Frequency = ListGetAngularFrequency / (2 PI)`, then `omega = 2 PI
    // Frequency`.
    let frequency = omega / (2.0 * std::f64::consts::PI);
    let omega = 2.0 * std::f64::consts::PI * frequency;
    log.info(CALLER, "Solving initially transient style system as harmonic one", 5);
    log.info(CALLER, &format!("Frequency value: {}", crate::fortran_fmt::es(frequency, 12, 3)), 4);

    // The Dirichlet rows of the real parts; the imaginary parts are zero.
    let id = ctx.sim.slots[ix].variable.ok_or_else(|| SolverError::fatal(CALLER, "solver has no variable"))?;
    let comps = ctx.sim.vars.component_names(id);
    let vname = ctx.sim.vars.by_id(id).name.clone();
    let dofs = ctx.sim.vars.by_id(id).dofs;
    let names: Vec<String> = if dofs == 1 { vec![vname.clone()] } else { (0..dofs).map(|k| comps.get(k).cloned().unwrap_or_else(|| format!("{vname} {}", k + 1))).collect() };
    for bc in &ctx.sim.model.boundary_conditions {
        if names.iter().any(|n| bc.get(&format!("{n} im")).is_some()) {
            return Err(unsupported("Dirichlet values of the imaginary part (<component> im)"));
        }
    }
    let fixed = ctx.dirichlet_values()?;
    let fixed2: Vec<(u32, f64)> = fixed.iter().flat_map(|&(r, v)| [(2 * r, v), (2 * r + 1, 0.0)]).collect();

    let sim = &mut *ctx.sim;
    let slot = &mut sim.slots[ix];
    let sys = slot.system.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no assembled system"))?;
    if sys.damp.is_some() {
        return Err(unsupported("a damping matrix"));
    }
    let g = &sys.matrix;
    let mass = sys.mass.as_ref().ok_or_else(|| SolverError::fatal(CALLER, "no mass matrix (DefaultUpdateMass)"))?;
    let w2 = omega * omega;
    let mut a = doubled(g, |j| g.values[j] - w2 * mass.values[j]);
    let n = g.nrows;
    let mut b = vec![0.0f64; 2 * n];
    for k in 0..n {
        b[2 * k] = sys.rhs[k];
        // `RHS_im`: the solvers that get here assemble no imaginary load.
    }
    let backend = sim.backends.assembly.as_ref().ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
    backend.try_apply_dirichlet(&mut a, &mut b, &fixed2, &slot.params)?;
    let mut x = vec![0.0f64; 2 * n];
    // `A % COMPLEX = .TRUE.`: SolveLinearSystem solves it as the complex
    // system it interleaves.
    let mut cp = slot.params.clone();
    cp.insert("Linear System Complex", elmer_core::Value::Logical(true));
    let converged = crate::solver::solve_other_matrix(
        &sim.backends.linear,
        &sim.model.simulation,
        &cp,
        &mut a,
        &mut b,
        &mut x,
        None,
        log,
        comm,
    )?;
    if !converged && slot.params.logical("linear system abort not converged").unwrap_or(true) {
        return Err(SolverError::fatal("SolveLinearSystem", "Failed convergence tolerances."));
    }
    let spec = crate::norm::NormSpec::from_params(&slot.params, 2 * dofs);
    let norm = crate::norm::compute_norm(&spec, &x, 2 * n, comm);
    let var = sim.vars.by_id_mut(id);
    let mut vec = vec![Complex64::new(0.0, 0.0); var.values.len()];
    for (j, z) in vec.iter_mut().enumerate().take(n) {
        *z = Complex64::new(x[2 * j], x[2 * j + 1]);
    }
    var.eigen_values = vec![Complex64::new(0.0, 0.0)];
    var.eigen_vectors = vec![vec];
    var.conv.norm = norm;
    Ok(norm)
}
