//! The run-dependent values of SaveScalars' `cpu time`, `wall time` and
//! `cpu memory` operators, as the POSIX branch of fem/src/CPUTime.c gives
//! them, read from `/proc/self` on Linux (zero elsewhere).

/// `CPUTime()`: the process's user CPU time in seconds (`ru_utime`).
pub(crate) fn cpu_time() -> f64 {
    // Field 14 of /proc/self/stat is utime in clock ticks (USER_HZ = 100).
    std::fs::read_to_string("/proc/self/stat")
        .ok()
        .and_then(|s| {
            let rest = &s[s.rfind(')')? + 2..];
            rest.split_whitespace().nth(11)?.parse::<f64>().ok()
        })
        .map_or(0.0, |ticks| ticks / 100.0)
}

/// `RealTime()`: seconds since the epoch (`gettimeofday`).
pub(crate) fn real_time() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0.0, |d| d.as_secs() as f64 + f64::from(d.subsec_micros()) * 1.0e-6)
}

/// `CPUMemory()`: the peak resident set in kB (`ru_maxrss`).
pub(crate) fn cpu_memory() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            let l = s.lines().find(|l| l.starts_with("VmHWM:"))?;
            l.split_whitespace().nth(1)?.parse::<f64>().ok()
        })
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timings_are_plausible() {
        assert!(cpu_time() >= 0.0);
        assert!(real_time() > 1.0e9);
        if cfg!(target_os = "linux") {
            assert!(cpu_memory() > 0.0);
        }
    }
}
