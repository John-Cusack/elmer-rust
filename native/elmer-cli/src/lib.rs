//! `ElmerSolver` / `ElmerSolver_mpi` front door of the Rust-native Elmer.
//!
//! * [`capability`]: manifest and refusal checker (format in `docs/manifest.md`).
//! * [`logger`], [`fortran_fmt`]: Elmer's `Info`/`Warn` and the GUI-parsed
//!   `ComputeChange` / `Time:` line formats.
//! * [`args`], [`version`], [`driver`], [`strangler`]: process behaviour, `-v`
//!   probe and the Fortran fallback launcher.
#![forbid(unsafe_code)]

pub mod args;
pub mod capability;
pub mod driver;
pub mod fortran_fmt;
pub mod logger;
pub mod strangler;
pub mod version;

pub use capability::{
    check, Capability, CapabilityManifest, CaseSummary, Refusal, Section, ValueClass,
};
pub use logger::{compute_change_line, ElmerLogger};

/// Process exit statuses. Fortran only distinguishes 0 and `STOP 1`; the native
/// CLI separates input errors and refusals so scripts can tell them apart. A
/// delegated Fortran run keeps its own status.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum ExitCode {
    /// Solved.
    Pass = 0,
    /// Runtime or numerical failure.
    Fail = 1,
    /// Bad command line, missing input file, unusable Fortran path.
    Input = 2,
    /// Case not supported natively and no Fortran fallback (`EX_CONFIG`).
    Refuse = 78,
}

impl ExitCode {
    pub fn code(self) -> i32 {
        self as i32
    }
}
