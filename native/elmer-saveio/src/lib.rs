//! Elmer output and restart file formats, byte-compatible with the Fortran
//! writers.
//!
//! | Module | Formats | Fortran reference (`fem/src/`) |
//! |---|---|---|
//! | [`fortran`] | `ES`/`E`/`I` edit descriptors, list-directed reals | gfortran runtime |
//! | [`result`] | `.result` model, encoder, [`read_result`](result::read_result), `.pos` | `ModelDescription.F90` `SaveResult`, `LoadRestartFile` |
//! | [`writer`] | [`write_result`](writer::write_result), [`ResultWriter`](writer::ResultWriter), cyclic file names | `SaveResult`, `ElmerSolver.F90` `SaveCurrent` |
//! | [`ep`] | ElmerPost `.ep`: [`write_ep`](ep::write_ep) | `ModelDescription.F90` `WritePostFile` |
//! | [`savedata`] | `SaveScalars` / `SaveLine` `.dat` and `.names` layout | `modules/SaveData/` |
//!
//! Every writer takes the dates it stamps (`!File started at:`) as an option,
//! so output can be compared byte for byte with a captured Fortran file.
//!
//! # What `elmer-core` does not carry
//!
//! [`elmer_core::field::FieldRef`] has no solver `Equation` name, no
//! complex values and no mesh-body/BC names. They enter through
//! [`writer::ResultOptions::equations`] and [`ep::EpOptions`] until a
//! `core-change` adds them; complex `.ep` output is refused rather than
//! approximated.
//!
//! The reader here is independent of `elmer-oracle`'s grader.
#![forbid(unsafe_code)]

pub mod ep;
pub mod error;
pub mod fortran;
pub mod result;
pub mod savedata;
pub mod writer;

pub use error::{Result, SaveIoError};
