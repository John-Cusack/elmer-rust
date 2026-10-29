//! Error type for `.sif` lexing and parsing.

use std::fmt;

/// A parse failure. Elmer calls `Fatal` on these; we return them so the
/// caller (later: capability refusal with a distinct exit code) can decide.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SifError {
    /// 1-based line number in the top-level file, or 0 when unknown.
    pub line: u32,
    /// What went wrong.
    pub msg: String,
    /// Include-file chain, outermost first.
    pub context: Vec<String>,
}

impl SifError {
    pub fn new(line: u32, msg: impl Into<String>) -> Self {
        SifError { line, msg: msg.into(), context: Vec::new() }
    }
}

impl fmt::Display for SifError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.context.is_empty() {
            write!(f, "sif:{}: {}", self.line, self.msg)
        } else {
            write!(f, "sif:{}: {} [in {}]", self.line, self.msg, self.context.join(" <- "))
        }
    }
}

impl std::error::Error for SifError {}
