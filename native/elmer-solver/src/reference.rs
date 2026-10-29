//! `CompareToReferenceSolution` (ElmerSolver.F90:1169-1342): the
//! `Reference Norm` / `Reference Solution` checks and the `TEST.PASSED`
//! file, with byte-identical messages.

use std::path::Path;

use elmer_core::{Communicator, Logger, MeshView, Params, NO_DOF};

use crate::error::SolverError;
use crate::fortran_fmt::{es, i2s};
use crate::sim::Simulation;

/// The Fortran routine's `SAVE`d counters.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReferenceTally {
    pub test_count: usize,
    pub pass_count: usize,
}

const CALLER: &str = "CompareToReferenceSolution";

/// One comparison pass over all solvers (the non-`Finalize` call).
pub fn compare_to_reference(
    sim: &mut Simulation,
    mesh: &dyn MeshView,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<(), SolverError> {
    if comm.rank() > 0 {
        return Ok(());
    }
    for ix in 0..sim.slots.len() {
        let slot = &sim.slots[ix];
        let p = &slot.params;
        let solver_id = slot.id;
        let ref_norm = match p.get("reference norm") {
            Some(v) => Some(crate::expr::eval_const(v, sim.backends.matc.as_deref())?),
            None => None,
        };
        let compare_solution = p.0.keys().any(|k| k.starts_with("reference solution"));
        if ref_norm.is_none() && !compare_solution {
            continue;
        }
        let Some(id) = slot.variable else {
            // The Fortran prints its loop index `i`, not the solver id.
            log.warn(
                CALLER,
                &format!("Variable in Solver {} not associated, cannot compare", i2s(solver_id as i64)),
            );
            continue;
        };
        let var = sim.vars.by_id(id);
        sim.reference.test_count += 1;
        let mut success = true;

        if let Some(ref_norm) = ref_norm {
            let tol = p.real("reference norm tolerance").unwrap_or(1.0e-5);
            let norm = var.conv.norm;
            let mut absolute = p.logical("reference norm absolute").unwrap_or(false);
            let mut err = (norm - ref_norm).abs();
            if !absolute {
                if ref_norm < f64::MIN_POSITIVE {
                    log.warn(CALLER, "Refenrece norm too small for relative error");
                    absolute = true;
                } else {
                    err /= ref_norm;
                }
            }
            if err > tol || err.is_nan() {
                log.warn(
                    CALLER,
                    &format!(
                        "Solver {solver_id} FAILED:  Norm ={}  RefNorm ={}",
                        es(norm, 15, 8),
                        es(ref_norm, 15, 8)
                    ),
                );
                success = false;
            } else {
                log.info(
                    CALLER,
                    &format!(
                        "Solver {solver_id} PASSED:  Norm ={}  RefNorm ={}",
                        es(norm, 15, 8),
                        es(ref_norm, 15, 8)
                    ),
                    3,
                );
            }
            let kind = if absolute { "Absolute" } else { "Relative" };
            log.info(CALLER, &format!("{kind} Error to reference norm:{}", es(err, 13, 6)), 3);
        }

        if compare_solution {
            let tol = p.real("reference solution tolerance").unwrap_or(1.0e-5);
            let dofs = var.dofs;
            let (mut n, mut rn, mut nn, mut en) = (0usize, 0.0f64, 0.0f64, 0.0f64);
            for node in 0..mesh.n_nodes() as u32 {
                let j = match var.perm.as_deref() {
                    Some(perm) => perm[node as usize],
                    None => node,
                };
                if j == NO_DOF {
                    continue;
                }
                for k in 1..=dofs {
                    let key = if dofs == 1 {
                        "reference solution".to_string()
                    } else {
                        format!("reference solution {k}")
                    };
                    let Some(rv) = p.get(&key) else { continue };
                    let fields = crate::expr::SimFields { vars: &sim.vars, mesh };
                    let refval = crate::expr::eval_at_nodes(rv, &[node], &fields, sim.backends.matc.as_deref())?[0];
                    let val = var.values[dofs * j as usize + k - 1];
                    rn += refval * refval;
                    nn += val * val;
                    en += (refval - val) * (refval - val);
                    n += 1;
                }
            }
            if comm.size() > 1 {
                log.warn("CompareToReferefenSolution", "Not implemented in parallel!");
            }
            if n == 0 {
                return Err(SolverError::fatal(CALLER, "Could not find any reference solution"));
            }
            let rn = (rn / n as f64).sqrt();
            let nn = (nn / n as f64).sqrt();
            let en = (en / n as f64).sqrt();
            if en > tol {
                log.warn(
                    CALLER,
                    &format!(
                        "Solver {solver_id} FAILED:  Solution = {}  RefSolution ={}",
                        es(nn, 15, 8),
                        es(rn, 15, 8)
                    ),
                );
                log.info(CALLER, &format!("Relative Error to reference solution:{}", es(en, 13, 6)), 3);
                success = false;
            } else {
                log.info(
                    CALLER,
                    &format!(
                        "Solver {solver_id} PASSED:  Solution ={}  RefSolution ={}",
                        es(nn, 15, 8),
                        es(rn, 15, 8)
                    ),
                    3,
                );
            }
        }
        if success {
            sim.reference.pass_count += 1;
        }
    }
    Ok(())
}

/// The `Finalize` call: the summary line and `TEST.PASSED` (or
/// `TEST.PASSED_<np>` in parallel) in `dir`, holding `1` or `0`. Nothing
/// happens when no solver was tested. Returns the overall success.
pub fn finalize_reference(
    tally: &ReferenceTally,
    dir: &Path,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<Option<bool>, SolverError> {
    if comm.rank() > 0 || tally.test_count == 0 {
        return Ok(None);
    }
    let success = tally.pass_count == tally.test_count;
    if success {
        log.info(CALLER, &format!("PASSED all {} tests!", tally.test_count), 3);
    } else {
        log.warn(
            CALLER,
            &format!("FAILED {} tests out of {}!", tally.test_count - tally.pass_count, tally.test_count),
        );
    }
    let name = if comm.size() > 1 {
        format!("TEST.PASSED_{}", comm.size())
    } else {
        "TEST.PASSED".to_string()
    };
    std::fs::write(dir.join(name), if success { "1\n" } else { "0\n" })?;
    Ok(Some(success))
}
