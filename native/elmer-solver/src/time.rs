//! [`ElmerTime`]: the [`TimeBackend`] built on `elmer-timeint`
//! (`Default1stOrderTime` -> `Add1stOrderTime`, SolverBasics.F90:245;
//! `Default2ndOrderTime` -> `Add2ndOrderTime`, SolverBasics.F90:494).
//! This crate owns the histories; `elmer-timeint` only reads them.

use elmer_core::{LocalSystem, NullLogger, Params};
use elmer_timeint::second::col as second_col;
use elmer_timeint::{add_1st_order_time, add_2nd_order_time, Aux1st, Scheme};

use crate::backend::{TimeBackend, TimeTerms};
use crate::error::SolverError;

/// The default time backend.
#[derive(Copy, Clone, Debug, Default)]
pub struct ElmerTime;

impl TimeBackend for ElmerTime {
    fn add_time_terms(&self, local: &mut LocalSystem, t: &TimeTerms<'_>) -> Result<(), SolverError> {
        let scheme = Scheme::from_params(t.params, t.simulation, 1)?.for_step(t.done_time);
        if matches!(scheme, Scheme::AdamsBashforth { .. } | Scheme::AdamsMoulton { .. }) {
            return Err(SolverError::Unsupported("Adams time stepping (per-element state)".into()));
        }
        // Order = MIN(Solver % DoneTime, Solver % Order) previous solutions.
        let order = t.done_time.min(t.order);
        let prev: Vec<&[f64]> = t.prev.iter().take(order).map(Vec::as_slice).collect();
        let mut aux = Aux1st {
            done_time: t.done_time,
            dts: t.dts,
            cur: Some(t.cur),
            adams: None,
            lumped: t.params.logical("lumped mass matrix").unwrap_or(false),
            logger: &NullLogger,
        };
        add_1st_order_time(local, t.dt, &scheme, &prev, &mut aux)?;
        Ok(())
    }

    /// `Bossak` (default, `Bossak Alpha` -0.05) or `Generalized-alpha` with
    /// `Generalized-alpha rinf` (InitializeTimestep, SolverBasics.F90:2732).
    fn add_2nd_order_time_terms(&self, local: &mut LocalSystem, t: &TimeTerms<'_>) -> Result<(), SolverError> {
        if t.params.logical("use global mass matrix") == Some(true) {
            return Err(SolverError::Unsupported("Use Global Mass Matrix".into()));
        }
        let scheme = Scheme::from_params(t.params, t.simulation, 2)?;
        let Some((a, b)) = scheme.second_order_params() else {
            return Ok(());
        };
        if t.prev.len() < 7 {
            return Err(SolverError::fatal("Default2ndOrderTime", "variable has no second-order history"));
        }
        let col = |c: usize| t.prev[c].as_slice();
        let lumped = t.params.logical("lumped mass matrix").unwrap_or(false);
        add_2nd_order_time(
            local,
            t.dt,
            col(second_col::X),
            col(second_col::V),
            col(second_col::A),
            col(second_col::A2),
            a,
            b,
            lumped,
        )?;
        Ok(())
    }
}
