//! `elmer_core::Logger` producing Elmer's `Info` / `Warn` / `Fatal` line formats
//! (Messages.F90:84-339) and the two stdout lines the ElmerGUI parses:
//! `ComputeChange … (NRM,RELC)` (SolverBasics.F90:4146-4152) and
//! `Time: n/N:  t` (ElmerSolver.F90:3127).

use std::io::Write;
use std::sync::Mutex;

use elmer_core::{Logger, Params, Value};

use crate::fortran_fmt::{fmt_es, fmt_g};

/// Output filtering and prefixes (`OutputLevelMask`, `OutputPrefix`,
/// `OutputCaller`, `OutputPE` in Messages.F90:60-66).
#[derive(Clone, Debug, PartialEq)]
pub struct OutputConfig {
    /// `OutputLevelMask(0:32)`.
    pub level_mask: [bool; 33],
    pub max_level: u8,
    /// `INFO:: ` before every Info line.
    pub prefix: bool,
    /// `Caller: ` before the message.
    pub caller: bool,
    /// `Part<n>: ` tag (`MaxOutputPE > 0` in Fortran); `None` for serial runs.
    pub part: Option<usize>,
}

impl Default for OutputConfig {
    /// Module defaults before a model is read: `MaxOutputLevel=32`, everything on.
    fn default() -> Self {
        OutputConfig {
            level_mask: [true; 33],
            max_level: 32,
            prefix: false,
            caller: true,
            part: None,
        }
    }
}

impl OutputConfig {
    /// ModelDescription.F90:145-161, from the `Simulation` section:
    /// `Min Output Level` (0), `Max Output Level` (10), `Output Level(n)`
    /// per-level switches, `Output Prefix` (false), `Output Caller` (true).
    pub fn from_simulation(sim: &dyn Params) -> Self {
        let min = sim.integer("Min Output Level").unwrap_or(0);
        let max = sim.integer("Max Output Level").unwrap_or(10);
        let mut mask = [true; 33];
        for (i, m) in mask.iter_mut().enumerate().take(32) {
            *m = (i as i32) >= min && (i as i32) <= max;
        }
        if let Some(Value::Integer(sw)) = sim.get("Output Level") {
            for (i, s) in sw.iter().enumerate().take(32) {
                mask[i] = *s != 0 && mask[i];
            }
        }
        OutputConfig {
            level_mask: mask,
            max_level: max.clamp(0, 32) as u8,
            prefix: sim.logical("Output Prefix").unwrap_or(false),
            caller: sim.logical("Output Caller").unwrap_or(true),
            part: None,
        }
    }
}

/// Elmer's message sink writing to any `Write` (stdout by default).
pub struct ElmerLogger {
    cfg: OutputConfig,
    sink: Mutex<Box<dyn Write + Send>>,
}

impl ElmerLogger {
    pub fn new(cfg: OutputConfig, sink: Box<dyn Write + Send>) -> Self {
        ElmerLogger {
            cfg,
            sink: Mutex::new(sink),
        }
    }

    /// Default configuration, writing to process stdout.
    pub fn stdout() -> Self {
        Self::new(OutputConfig::default(), Box::new(std::io::stdout()))
    }

    pub fn config(&self) -> &OutputConfig {
        &self.cfg
    }

    /// Replaces the filter, as ModelDescription does once `Simulation` is read.
    pub fn set_config(&mut self, cfg: OutputConfig) {
        self.cfg = cfg;
    }

    /// Writes `line` unfiltered and unprefixed, for the plain `WRITE(*,…)`
    /// statements Elmer uses outside `Info` (banner, total-time lines).
    pub fn raw_line(&self, line: &str) {
        self.write_line(line);
    }

    fn write_line(&self, line: &str) {
        if let Ok(mut s) = self.sink.lock() {
            // A closed pipe must not abort the solver; Fortran ignores it too.
            let _ = writeln!(s, "{line}");
            let _ = s.flush();
        }
    }

    /// The line `Info(caller, msg, Level=level)` prints, or `None` if filtered.
    pub fn format_info(&self, caller: &str, msg: &str, level: u8) -> Option<String> {
        if level > self.cfg.max_level || !self.cfg.level_mask[usize::from(level.min(32))] {
            return None;
        }
        let mut line = String::new();
        if self.cfg.prefix {
            line.push_str("INFO:: ");
        }
        if self.cfg.caller {
            line.push_str(caller.trim_end());
            line.push_str(": ");
        }
        if let Some(pe) = self.cfg.part {
            line.push_str(&format!("Part{pe}: "));
        }
        line.push_str(msg.trim_end_matches(' '));
        Some(line)
    }

    /// The line `Warn(caller, msg)` prints, or `None` if level 2 is masked.
    pub fn format_warn(&self, caller: &str, msg: &str) -> Option<String> {
        if !self.cfg.level_mask[2] {
            return None;
        }
        let msg = msg.trim_end_matches(' ');
        Some(match self.cfg.part {
            Some(pe) => format!("WARNING:: {}: Part{pe}:{msg}", caller.trim_end()),
            None => format!("WARNING:: {}: {msg}", caller.trim_end()),
        })
    }

    /// `Fatal(caller, msg)`: prints and returns; the caller exits with
    /// [`crate::ExitCode::Fail`] (Fortran `STOP 1`).
    pub fn fatal(&self, caller: &str, msg: &str) {
        if self.cfg.level_mask[0] {
            self.write_line(&format!(
                "ERROR:: {}: {}",
                caller.trim_end(),
                msg.trim_end_matches(' ')
            ));
        }
    }

    /// Prints the per-iteration convergence line (level 3, caller
    /// `ComputeChange`).
    pub fn compute_change(&self, solver: &str, iter: u32, norm: f64, change: f64, steady: bool) {
        self.info(
            "ComputeChange",
            &compute_change_message(solver, iter, norm, change, steady),
            3,
        );
    }

    /// Prints the per-timestep banner line (level 3, caller `MAIN`).
    pub fn time_step(&self, cum_step: u32, total_steps: u32, time: f64) {
        self.info("MAIN", &time_message(cum_step, total_steps, time), 3);
    }
}

impl Logger for ElmerLogger {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        if let Some(line) = self.format_info(caller, msg, level) {
            self.write_line(&line);
        }
    }

    fn warn(&self, caller: &str, msg: &str) {
        if let Some(line) = self.format_warn(caller, msg) {
            self.write_line(&line);
        }
    }
}

/// `WRITE(Message,'(a,g15.8,g15.8,a)')` of SolverBasics.F90:4146-4152. `iter`
/// is the nonlinear or steady-state iteration number (`IterNo`).
pub fn compute_change_message(
    solver: &str,
    iter: u32,
    norm: f64,
    change: f64,
    steady: bool,
) -> String {
    let kind = if steady { "SS" } else { "NS" };
    format!(
        "{kind} (ITER={iter}) (NRM,RELC): ({}{} ) :: {}",
        fmt_g(norm, 15, 8),
        fmt_g(change, 15, 8),
        solver.trim_end()
    )
}

/// The full stdout line, `ComputeChange: NS (ITER=1) (NRM,RELC): ( … ) :: heat`.
///
/// Deviation from the swarm plan signature: `iter` is required because the
/// line carries `(ITER=n)`.
pub fn compute_change_line(solver: &str, iter: u32, norm: f64, change: f64, steady: bool) -> String {
    format!(
        "ComputeChange: {}",
        compute_change_message(solver, iter, norm, change, steady)
    )
}

/// `'Time: '//i2s(cum)//'/'//i2s(steps)//':'` then `ES12.3` (ElmerSolver.F90:3127).
pub fn time_message(cum_step: u32, total_steps: u32, time: f64) -> String {
    format!("Time: {cum_step}/{total_steps}:{}", fmt_es(time, 12, 3))
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::ParamMap;
    use std::sync::{Arc, Mutex};

    #[derive(Clone, Default)]
    struct Buf(Arc<Mutex<Vec<u8>>>);
    impl Write for Buf {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    impl Buf {
        fn text(&self) -> String {
            String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
        }
    }

    fn logger(cfg: OutputConfig) -> (ElmerLogger, Buf) {
        let b = Buf::default();
        (ElmerLogger::new(cfg, Box::new(b.clone())), b)
    }

    #[test]
    fn blank_info_keeps_caller_colon_space() {
        let (l, b) = logger(OutputConfig::default());
        l.info("MAIN", " ", 4);
        assert_eq!(b.text(), "MAIN: \n");
    }

    #[test]
    fn level_filter_follows_simulation_section() {
        let mut sim = ParamMap::new();
        sim.insert("Max Output Level", Value::Integer(vec![5]));
        let (l, b) = logger(OutputConfig::from_simulation(&sim));
        l.info("A", "shown", 5);
        l.info("A", "hidden", 6);
        assert_eq!(b.text(), "A: shown\n");
    }

    #[test]
    fn output_level_switches_and_defaults() {
        let sim = ParamMap::new();
        let cfg = OutputConfig::from_simulation(&sim);
        assert_eq!(cfg.max_level, 10);
        assert!(cfg.caller && !cfg.prefix);
        let mut sim = ParamMap::new();
        sim.insert("Output Level", Value::Integer(vec![1, 1, 0, 1]));
        let cfg = OutputConfig::from_simulation(&sim);
        assert!(!cfg.level_mask[2] && cfg.level_mask[3]);
    }

    #[test]
    fn prefix_and_part_tag() {
        let cfg = OutputConfig {
            prefix: true,
            part: Some(2),
            ..OutputConfig::default()
        };
        let (l, b) = logger(cfg);
        l.info("X", "msg  ", 4);
        l.warn("X", "careful");
        assert_eq!(b.text(), "INFO:: X: Part2: msg\nWARNING:: X: Part2:careful\n");
    }

    #[test]
    fn warn_and_fatal_formats() {
        let (l, b) = logger(OutputConfig::default());
        l.warn("Solver", "may require \"Procedure\" to operate!");
        l.fatal("MAIN", "Unable to find ELMERSOLVER_STARTINFO, can not execute.");
        assert_eq!(
            b.text(),
            "WARNING:: Solver: may require \"Procedure\" to operate!\n\
             ERROR:: MAIN: Unable to find ELMERSOLVER_STARTINFO, can not execute.\n"
        );
    }

    #[test]
    fn warn_masked_at_level_two() {
        let mut cfg = OutputConfig::default();
        cfg.level_mask[2] = false;
        let (l, b) = logger(cfg);
        l.warn("Solver", "x");
        assert_eq!(b.text(), "");
    }

    #[test]
    fn compute_change_line_shape() {
        assert_eq!(
            compute_change_line("heat", 1, 0.75696824e-2, 2.0, false),
            "ComputeChange: NS (ITER=1) (NRM,RELC): ( 0.75696824E-02  2.0000000     ) :: heat"
        );
        assert_eq!(
            compute_change_line("heat equation", 1, 0.93065875, 2.0, true),
            "ComputeChange: SS (ITER=1) (NRM,RELC): ( 0.93065875      2.0000000     ) :: heat equation"
        );
    }

    #[test]
    fn time_line_shape() {
        assert_eq!(time_message(1, 8, 0.1), "Time: 1/8:   1.000E-01");
    }
}
