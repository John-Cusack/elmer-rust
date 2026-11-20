//! `ReloadSolution` of the `ReloadData` module (fem/src/modules/
//! ReloadData.F90) and the `LoadRestartFile` it calls
//! (ModelDescription.F90:4325).
//!
//! Each visit reads one stored step of `Reload Solution File` (under the
//! mesh directory, `Mesh % Name`) into the existing variables. The step
//! cycles over `Reload Range Minimum`..`Maximum`, starting from `Reload
//! Starting Position`, every `Reload Reading Intervals` visits. The loaded
//! step also sets the `Time` and `Timestep` globals, unless the Simulation
//! has `Restart Time Ignore`, and the main loop carries on from them.
//!
//! The file may be a version 2 or 3 `.result`, ASCII or binary, read by
//! [`elmer_solver::load_restart`] with the Simulation's restart keywords.
//! Refused: runs on more than one partition.

use std::path::PathBuf;

use elmer_core::{ParamMap, Params};
use elmer_saveio::result::{read_result, ResultFile};
use elmer_solver::{Solver, SolverCtx, SolverError};

const CALLER: &str = "ReloadSolution";
const LOAD: &str = "LoadRestartFile";

pub(crate) struct ReloadSolution {
    file: Option<String>,
    first: i32,
    last: i32,
    between: i32,
    round: i32,
    visit: i32,
    loaded: Option<(PathBuf, ResultFile)>,
}

impl ReloadSolution {
    pub(crate) fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        let first = p.integer("reload range minimum").unwrap_or(1);
        let last = p.integer("reload range maximum").unwrap_or(100_000);
        let mut round = 0;
        if let Some(start) = p.integer("reload starting position") {
            round += start - first;
        }
        let between = p.integer("reload reading intervals").unwrap_or(1);
        if between < 1 {
            return Err(SolverError::fatal(CALLER, "Reload Reading Intervals must be positive"));
        }
        if last < first {
            return Err(SolverError::fatal(CALLER, "Reload Range Maximum is below the minimum"));
        }
        Ok(ReloadSolution {
            file: p.string("reload solution file").map(str::to_string),
            first,
            last,
            between,
            round,
            visit: 0,
            loaded: None,
        })
    }
}

impl Solver for ReloadSolution {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let visit = self.visit;
        self.visit += 1;
        if visit % self.between != 0 {
            return Ok(());
        }
        let Some(file) = self.file.clone() else { return Ok(()) };
        if ctx.comm.size() > 1 {
            return Err(SolverError::Unsupported("ReloadSolution: parallel runs".into()));
        }
        if visit == 0 {
            ctx.info(CALLER, &format!(" Reading old solution from file {file}"), 4);
        }
        self.round += 1;
        let span = 1 + self.last - self.first;
        let k = match self.round % span {
            0 => self.last,
            r => r + self.first - 1,
        };
        let path = match ctx.sim.model.header.string("mesh name") {
            Some(m) if !m.trim().is_empty() => format!("{}/{file}", m.trim()),
            _ => file,
        };
        ctx.info(CALLER, &format!(" Opening file: {path}"), 4);
        ctx.info(CALLER, &format!(" Loading Timestep: {k:>12}"), 4);
        let path = PathBuf::from(path);
        if self.loaded.as_ref().is_none_or(|(p, _)| *p != path) {
            let rf = read_result(&path).map_err(|e| SolverError::fatal(LOAD, format!("{}: {e}", path.display())))?;
            self.loaded = Some((path, rf));
        }
        let (_, rf) = self.loaded.as_ref().expect("just loaded");
        let simp = ctx.simulation_params().clone();
        elmer_solver::load_restart(ctx.sim, rf, k, &simp, ctx.log)
    }
}
