//! Test-local stand-in for S1's `elmer-mesh` (mesh reader). Replace it with
//! the real crate once it merges; every refusal it causes is labelled
//! `stand-in` in the report. Elements are evaluated by `elmer-fe`, which
//! `elmer-run` installs.
#![allow(clippy::needless_range_loop)]

pub mod mesh;
