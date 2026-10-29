//! Version probe output of `ElmerSolver -v` (Solver.F90:53-100,
//! ElmerSolver.F90:164-331), as captured from the pristine 26.2-devel binary.

use std::time::{SystemTime, UNIX_EPOCH};

use elmer_core::Logger;

use crate::logger::ElmerLogger;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VersionInfo {
    /// `ELMER_FEM_VERSION`, printed on the `STARTED AT` line.
    pub fem_version: String,
    /// `GetVersion()`.
    pub version: String,
    pub branch: String,
    pub revision: String,
    pub compiled: String,
}

impl Default for VersionInfo {
    fn default() -> Self {
        VersionInfo {
            fem_version: "26.2".into(),
            version: "26.2".into(),
            branch: "devel".into(),
            revision: "9f6af2f85".into(),
            compiled: "2026-09-23".into(),
        }
    }
}

/// Fortran `FormatDate()`: `YYYY/MM/DD HH:MM:SS`. Time is UTC; the pristine
/// binary prints local time, which the GUI never parses.
pub fn format_date(now: SystemTime) -> String {
    let secs = now.duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs());
    let (days, rem) = ((secs / 86_400) as i64, secs % 86_400);
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}/{m:02}/{d:02} {:02}:{:02}:{:02}",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

/// `ELMERSOLVER_OUTPUT_TOTAL` starting with `0` or `false` silences the
/// `STARTED AT` / `FINISHED AT` / `TOTAL TIME` lines (Solver.F90:63-64).
pub fn total_output_silenced(value: Option<&str>) -> bool {
    value.is_some_and(|v| v.starts_with('0') || v.starts_with("false"))
}

/// Writes the `-v` probe to `logger`'s sink in the pristine order. `silent` is
/// the `-s` flag (drops the welcome banner only), `total_silent` is
/// [`total_output_silenced`].
pub fn print_probe(
    logger: &ElmerLogger,
    info: &VersionInfo,
    silent: bool,
    total_silent: bool,
    now: SystemTime,
    cpu_real: (f64, f64),
) {
    let raw = |s: &str| logger.raw_line(s);
    if !total_silent {
        raw(&format!(
            "ELMER SOLVER (v {}) STARTED AT: {}",
            info.fem_version,
            format_date(now)
        ));
    }
    raw(" Started inside library code");
    raw(" Going parallel initialization");
    logger.info("ParCommInit", "Initialize #PEs: 1", 4);
    raw(" Back from parallel initialization");
    if !silent {
        let main = |s: &str, lvl: u8| logger.info("MAIN", s, lvl);
        main(" ", 4);
        main("=============================================================", 4);
        main("ElmerSolver finite element software, Welcome!                ", 4);
        main("This program is free software licensed under (L)GPL          ", 4);
        main("Copyright April 1st 1995 -> , CSC - IT Center for Science Ltd.     ", 4);
        main("Repository https://www.github.org/elmerCSC/elmerfem          ", 4);
        main("Homepage https://www.elmerfem.org/,  Email elmeradm@csc.fi   ", 4);
        main(
            &format!(
                "Version: {}-{} (Rev: {}, Compiled: {})",
                info.version, info.branch, info.revision, info.compiled
            ),
            4,
        );
        main(" Running one task without MPI parallelization.", 10);
        main(" Running with just one thread per task.", 10);
        // Extra line (swarm plan C08: version output names the backend).
        main(" Backend: Rust-native elmer-cli, Fortran delegation via ELMER_FORTRAN_SOLVER.", 4);
        main("=============================================================", 4);
    }
    if !total_silent {
        raw(&format!(
            "SOLVER TOTAL TIME(CPU,REAL): {:>12.2}{:>12.2}",
            cpu_real.0, cpu_real.1
        ));
        raw(&format!("ELMER SOLVER FINISHED AT: {}", format_date(now)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn date_format() {
        assert_eq!(format_date(UNIX_EPOCH), "1970/01/01 00:00:00");
        // 2026-10-29 10:25:07 UTC
        let t = UNIX_EPOCH + Duration::from_secs(1_793_269_507);
        assert_eq!(format_date(t), "2026/10/29 10:25:07");
        let leap = UNIX_EPOCH + Duration::from_secs(951_782_400); // 2000-02-29
        assert_eq!(format_date(leap), "2000/02/29 00:00:00");
    }

    #[test]
    fn total_silence_rule() {
        assert!(total_output_silenced(Some("0")));
        assert!(total_output_silenced(Some("false")));
        assert!(!total_output_silenced(Some("1")));
        assert!(!total_output_silenced(None));
    }
}
