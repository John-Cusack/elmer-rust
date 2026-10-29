//! [`ElmerTime`]: the [`TimeBackend`] built on `elmer-timeint`
//! (`Default1stOrderTime` -> `Add1stOrderTime`, SolverBasics.F90:245).
//! This crate owns the histories; `elmer-timeint` only reads them.

use elmer_core::{LocalSystem, NullLogger, Params};
use elmer_timeint::{add_1st_order_time, Aux1st, Scheme};

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
}
