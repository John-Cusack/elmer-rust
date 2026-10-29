//! `tresca/Tresca.f90`: `SUBROUTINE TrescaSolver`, the `Procedure =
//! "Tresca" "TrescaSolver"` that turns StressSolve's `Stress` field into a
//! nodal Tresca stress `sigma_max - sigma_min`.
//!
//! It has no variable of its own and solves nothing: on each visit it reads
//! `Stress` (4 components in 2D, `xx yy zz xy`; 6 in 3D, then `yz xz`),
//! creates `Tresca` on `Stress`'s permutation on the first visit
//! (`VariableAdd`) and fills it slot by slot from the eigenvalues of the
//! symmetric stress tensor (`DSYEV('N', 'U', 3, ...)`).

use elmer_core::{FieldKind, ParamMap};
use elmer_solver::{Solver, SolverCtx, SolverError, Variable};

const CALLER: &str = "TrescaSolver";

/// `TrescaValues` is `SAVE`d: the variable this instance created.
#[derive(Default)]
pub struct TrescaSolver {
    created: bool,
}

/// [`SolverFactory`](elmer_solver::SolverFactory) for [`TrescaSolver`].
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(TrescaSolver::default()))
}

/// Eigenvalues of the symmetric 3x3 `a` (upper triangle read), by cyclic
/// Jacobi rotations; `DSYEV` with `JOBZ = 'N'`. Unsorted.
pub fn symmetric_eigenvalues(a: [[f64; 3]; 3]) -> [f64; 3] {
    let mut m = a;
    for i in 0..3 {
        for j in 0..i {
            m[i][j] = m[j][i];
        }
    }
    for _ in 0..50 {
        let off = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
        let scale = m[0][0] * m[0][0] + m[1][1] * m[1][1] + m[2][2] * m[2][2];
        if off <= f64::EPSILON * f64::EPSILON * scale || off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q] == 0.0 {
                continue;
            }
            let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
            // `signum(0.0)` is 1: a 45 degree rotation.
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for row in m.iter_mut() {
                let (mkp, mkq) = (row[p], row[q]);
                row[p] = c * mkp - s * mkq;
                row[q] = s * mkp + c * mkq;
            }
            let (rp, rq) = (m[p], m[q]);
            for k in 0..3 {
                m[p][k] = c * rp[k] - s * rq[k];
                m[q][k] = s * rp[k] + c * rq[k];
            }
        }
    }
    [m[0][0], m[1][1], m[2][2]]
}

/// The Tresca stress of one node's stress components (Elmer order
/// `xx yy zz xy [yz xz]`).
pub fn tresca(stress: &[f64]) -> f64 {
    let mut m = [[0.0; 3]; 3];
    m[0][0] = stress[0];
    m[1][1] = stress[1];
    m[2][2] = stress[2];
    m[0][1] = stress[3];
    if stress.len() == 6 {
        m[1][2] = stress[4];
        m[0][2] = stress[5];
    }
    let e = symmetric_eigenvalues(m);
    let max = e.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let min = e.iter().copied().fold(f64::INFINITY, f64::min);
    max - min
}

impl Solver for TrescaSolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        // `STOP 'Error: Stress was not found'`
        let stress =
            ctx.variable("stress").ok_or_else(|| SolverError::fatal(CALLER, "Error: Stress was not found"))?;
        let dofs = stress.dofs;
        if dofs != 6 && dofs != 4 {
            return Err(SolverError::fatal(CALLER, "Error: Unexpected number of DOFs for stress"));
        }
        let values: Vec<f64> = stress.values.chunks_exact(dofs).map(tresca).collect();
        let perm = ctx.sim.vars.get("stress").and_then(|v| v.perm.clone());

        match ctx.sim.vars.get_mut("tresca") {
            Some(v) if self.created => {
                v.values = values;
                v.values_changed = true;
            }
            // The Fortran writes through its SAVEd TrescaValues, which only
            // the first visit associates; a Tresca from elsewhere leaves it
            // unassociated.
            Some(_) => {
                return Err(SolverError::Unsupported(
                    "TrescaSolver: a Tresca variable that this solver did not create".into(),
                ))
            }
            None => {
                let mut v = Variable::new("Tresca", FieldKind::Nodal, 1, perm, values);
                v.values_changed = true;
                ctx.sim.vars.add(v)?;
                self.created = true;
            }
        }
        Ok(())
    }
}
