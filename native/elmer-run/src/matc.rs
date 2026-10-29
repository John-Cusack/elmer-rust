//! [`SifMatc`]: the MATC engine `elmer-sif` parsed the `.sif` with, serving
//! run-time `Real MATC` values, so they see the variables and functions the
//! `.sif`'s `$...$` statements defined (one MATC environment per process).

use elmer_expr::MatcEngine;
use elmer_sif::matc::matc_eval_args;

/// `SetGetMatcParams` + `GetMatcRealArray` (Lists.F90:128-190) on
/// `elmer-sif`'s MATC binding.
#[derive(Clone, Copy, Debug, Default)]
pub struct SifMatc;

impl MatcEngine for SifMatc {
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.eval_with("tx", expr, args)
    }

    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        // One lock hold for binding and evaluation, so concurrent cases in one
        // process cannot interleave. The result comes back as MATC prints it,
        // which is what the Fortran reads too.
        matc_eval_args(var, args, expr)?
            .split_whitespace()
            .map(|t| t.parse::<f64>().map_err(|e| format!("MATC result {t:?}: {e}")))
            .collect()
    }
}
