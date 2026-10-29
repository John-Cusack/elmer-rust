//! Quick-suite coverage harness: how much of Elmer's quick test suite the
//! Rust-native path (`elmer-run` behind `elmer-cli`) already runs, and what
//! is missing, ranked by the tests each gap would unlock.
//!
//! * [`discover`]: the 482 quick tests from ctest's registry.
//! * [`stage`]: scratch copy, ElmerGrid mesh, the test's solver steps.
//! * [`frontend`], [`standin`]: `.sif` + mesh DB -> `elmer-run` case, with
//!   test-local stand-ins for S1's mesh reader and basis.
//! * [`analyze`]: the static screen naming every blocker of a case.
//! * [`child`]: one solver run per process.
//! * [`runner`], [`record`], [`report`]: classification and `COVERAGE.md`.

pub mod analyze;
pub mod child;
pub mod discover;
pub mod frontend;
pub mod record;
pub mod report;
pub mod runner;
pub mod stage;
pub mod standin;
