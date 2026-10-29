//! `elmer-sif`: Rust-native Elmer `.sif` input-file lexer/parser.
//!
//! Replicates Elmer's `ReadAndTrim` (`GeneralUtils.F90`) and `LoadInputFile`
//! / `SectionContents` (`ModelDescription.F90`), including MATC `$...$`
//! evaluation at lex time through the vendored C MATC library, `Include`
//! splicing, `(*)` sizes, tables, `Procedure` refs, `RUN` blocks, and
//! keyword typing generated from the `SOLVER.KEYWORDS` snapshot.

pub mod adapters;
pub mod complete;
pub mod error;
pub mod keywords;
pub mod lex;
pub mod matc;
pub mod model;
pub mod parse;

pub use adapters::section_to_param_map;
pub use complete::complete;
pub use error::SifError;
pub use model::{Entry, Model, ProcNumeric, Section, Sif, Table, TableKind, Value};
pub use parse::{ParseResult, parse_file};
