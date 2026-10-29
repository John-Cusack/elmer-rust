//! Message sink, as Elmer's Info/Warn (Messages.F90).

pub trait Logger: Sync {
    fn info(&self, caller: &str, msg: &str, level: u8);
    fn warn(&self, caller: &str, msg: &str);
}

/// Discards everything.
#[derive(Copy, Clone, Debug, Default)]
pub struct NullLogger;

impl Logger for NullLogger {
    fn info(&self, _caller: &str, _msg: &str, _level: u8) {}
    fn warn(&self, _caller: &str, _msg: &str) {}
}
