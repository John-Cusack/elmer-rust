//! The ElmerGUI solver-log algorithm (`mainwindow.cpp:6926-6978`), ported
//! line by line, including Qt's lenient `toDouble` (failure reads as 0).
//!
//! For each line (trimmed): a line containing `Time:` yields the last
//! space-separated token as the current time; a line containing
//! `ComputeChange` yields the solver name (text after the last `:`), the
//! norm (first token after the first `(` that follows `NRM,RELC`) and the
//! relative change (the next non-empty token, floored at 1e-16).

#[derive(Clone, Debug, PartialEq)]
pub struct Change {
    pub name: String,
    pub norm: f64,
    pub relative_change: f64,
    /// `NS`/`SS` tag the GUI uses to pick the convergence curve.
    pub kind: ChangeKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChangeKind {
    NonlinearSystem,
    SteadyState,
    Other,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GuiLog {
    pub times: Vec<f64>,
    pub changes: Vec<Change>,
    /// A line matched one of the GUI's end-of-run highlight patterns.
    pub finished: bool,
}

/// Qt `QString::toDouble()` on a token: the whole token must be a number.
fn qt_to_double(s: &str) -> f64 {
    s.trim().parse().unwrap_or(0.0)
}

/// `solverlogwindow.cpp:107-110`, case-insensitive substrings.
const END_PATTERNS: [&str; 4] = [
    " Elmer Solver: ALL DONE ",
    "ElmerSolver: The end",
    "SOLVER TOTAL TIME",
    "ELMER SOLVER FINISHED AT:",
];

pub fn parse_time(line: &str) -> Option<f64> {
    let tmp = line.trim();
    if !tmp.contains("Time:") {
        return None;
    }
    Some(qt_to_double(tmp.split(' ').next_back().unwrap_or("")))
}

pub fn parse_change(line: &str) -> Option<Change> {
    let tmp = line.trim();
    if !tmp.contains("ComputeChange") {
        return None;
    }
    let name = tmp.split(':').next_back().unwrap_or("").trim().to_string();
    // Qt `indexOf` is -1 when absent and `mid(-1)` is the whole string.
    let rest = tmp.find("NRM,RELC").map_or(tmp, |n| &tmp[n..]);
    let parts: Vec<&str> = rest.split('(').collect();
    let (mut norm, mut rel) = (0.0, 0.0);
    if parts.len() >= 2 {
        let toks: Vec<&str> = parts[1].trim().split(' ').collect();
        norm = qt_to_double(toks[0]);
        let mut pos = 1;
        while toks.get(pos).is_some_and(|t| t.trim().is_empty()) {
            pos += 1;
        }
        // The GUI indexes past the end here and would crash; read as 0.
        rel = f64::max(qt_to_double(toks.get(pos).copied().unwrap_or("")), 1.0e-16);
    }
    let kind = if tmp.contains("NS") {
        ChangeKind::NonlinearSystem
    } else if tmp.contains("SS") {
        ChangeKind::SteadyState
    } else {
        ChangeKind::Other
    };
    Some(Change {
        name,
        norm,
        relative_change: rel,
        kind,
    })
}

pub fn is_end_line(line: &str) -> bool {
    let low = line.to_ascii_lowercase();
    END_PATTERNS
        .iter()
        .any(|p| low.contains(&p.to_ascii_lowercase()))
}

pub fn parse_log(text: &str) -> GuiLog {
    let mut log = GuiLog::default();
    for line in text.split('\n') {
        if let Some(t) = parse_time(line) {
            log.times.push(t);
        }
        if let Some(c) = parse_change(line) {
            log.changes.push(c);
        }
        log.finished |= is_end_line(line);
    }
    log
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_pristine_style_compute_change() {
        let c = parse_change(
            "ComputeChange: NS (ITER=1) (NRM,RELC): ( 0.75696824E-02  2.0000000     ) :: heat",
        )
        .unwrap();
        assert_eq!(c.name, "heat");
        assert!((c.norm - 0.75696824e-2).abs() < 1e-12);
        assert_eq!(c.relative_change, 2.0);
        assert_eq!(c.kind, ChangeKind::NonlinearSystem);
        let c = parse_change(
            "ComputeChange: SS (ITER=1) (NRM,RELC): (  0.93065875      2.0000000     ) :: heat equation",
        )
        .unwrap();
        assert_eq!(c.name, "heat equation");
        assert_eq!(c.kind, ChangeKind::SteadyState);
    }

    #[test]
    fn relative_change_is_floored() {
        let c = parse_change("ComputeChange: NS (ITER=2) (NRM,RELC): ( 1.0  0.0 ) :: s").unwrap();
        assert_eq!(c.relative_change, 1.0e-16);
    }

    #[test]
    fn time_is_the_last_token() {
        assert_eq!(parse_time("Time: 1/8:   1.000E-01"), Some(0.1));
        assert_eq!(parse_time("  Time: 2/8:   2.000E-01  "), Some(0.2));
        assert_eq!(parse_time("no clock here"), None);
    }

    #[test]
    fn whole_log_and_end_detection() {
        let log = parse_log(
            "Time: 1/1:   1.000E+00\nComputeChange: NS (ITER=1) (NRM,RELC): ( 1.0 2.0 ) :: t\n\
             ELMER SOLVER FINISHED AT: 2026/10/29",
        );
        assert_eq!(log.times, [1.0]);
        assert_eq!(log.changes.len(), 1);
        assert!(log.finished);
        assert!(!parse_log("ElmerSolver: running").finished);
    }
}
