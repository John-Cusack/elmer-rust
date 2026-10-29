//! Keyword typing for the Rust-native Elmer front end.
//!
//! The catalog is generated from a pinned snapshot of Elmer's
//! `SOLVER.KEYWORDS` (`data/SOLVER.KEYWORDS`) with the same lookup rules as
//! the reference `CheckKeyword`: section aliasing, case/whitespace
//! normalization, `$` (MATC) family expansion and the trailing-number retry.
//! It says what *type* a keyword has and nothing about whether a physics
//! module supports it; that belongs to the CLI's capability manifest.
//!
//! The versioned JSON rendering lives at `data/keywords.catalog.json` and is
//! regenerated with `gen-keywords`; a test fails if it drifts.
#![forbid(unsafe_code)]

pub mod catalog;
pub mod matc;
pub mod sha256;

use std::sync::OnceLock;

pub use catalog::{
    canonical_section, Catalog, CatalogError, Diagnostic, DiagnosticKind, Entry, KeywordType,
    Origin, RawType, SourceInfo, SCHEMA_VERSION,
};

/// The pinned snapshot (text of `fem/src/SOLVER.KEYWORDS`).
pub const SNAPSHOT: &str = include_str!("../data/SOLVER.KEYWORDS");
/// Path of the snapshot inside the Elmer fork it was taken from.
pub const SNAPSHOT_PATH: &str = "fem/src/SOLVER.KEYWORDS";
/// Commit of `/home/john/runs/elmer-fork` the snapshot was taken at.
pub const SNAPSHOT_REVISION: &str = "c0b4bc6f1a4be5538646c2de2c9e54594814c113";

/// The catalog built from [`SNAPSHOT`], parsed once.
pub fn builtin() -> &'static Catalog {
    static CATALOG: OnceLock<Catalog> = OnceLock::new();
    CATALOG.get_or_init(|| {
        Catalog::from_source(SNAPSHOT, SNAPSHOT_PATH, SNAPSHOT_REVISION)
            .expect("the pinned SOLVER.KEYWORDS snapshot must expand")
    })
}

/// Type of `key` in the SIF section named `section` (`"Solver 1"`,
/// `"Boundary Condition 2"`, ...), or `None` when the reference would call it
/// an unlisted keyword.
pub fn keyword_type(section: &str, key: &str) -> Option<KeywordType> {
    builtin().keyword_type(section, key)
}
