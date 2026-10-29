//! [`SifMatc`]: the MATC engine `elmer-sif` parsed the `.sif` with, serving
//! run-time `Real MATC` values, so they see the variables and functions the
//! `.sif`'s `$...$` statements defined (one MATC environment per process).

use elmer_expr::MatcEngine;
use elmer_sif::matc::matc_eval;

/// `SetGetMatcParams` + `GetMatcRealArray` (Lists.F90:128-190) on
/// `elmer-sif`'s MATC binding.
#[derive(Clone, Copy, Debug, Default)]
pub struct SifMatc;

impl MatcEngine for SifMatc {
    fn eval(&self, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        self.eval_with("tx", expr, args)
    }

    fn eval_with(&self, var: &str, expr: &str, args: &[f64]) -> Result<Vec<f64>, String> {
        // MatcSetParams binds the values exactly; a shortest round-trip
        // literal parses back to the same doubles.
        let set = if args.is_empty() {
            format!("{var}=0")
        } else {
            let vals: Vec<String> = args.iter().map(|a| format!("{a:e}")).collect();
            format!("{var}={}", vals.join(" "))
        };
        matc_eval(&set)?;
        // The result comes back as MATC prints it, which is what the Fortran
        // reads too.
        matc_eval(expr)?
            .split_whitespace()
            .map(|t| t.parse::<f64>().map_err(|e| format!("MATC result {t:?}: {e}")))
            .collect()
    }
}
