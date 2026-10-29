//! VTU/PVTU/PVD output with `ResultOutputSolver` semantics.
//!
//! Ports `VtuOutputSolver.F90` (file naming, field selection, element and
//! node order, ASCII and appended-binary layout, `.pvd` and `.pvtu` wrappers),
//! the `Post File = "*.vtu"` auto-solver of `ElmerSolver.F90`, and the format
//! defaults of `ResultOutputSolver.F90`. Output is byte-for-byte what the
//! Fortran writes for the same mesh and field data.
//!
//! The writers read the model through `elmer-core`'s [`MeshView`] and
//! [`FieldSource`]. What those traits cannot express is refused rather than
//! approximated, with an `io::ErrorKind::Unsupported` error: discontinuous
//! Galerkin saving, eigen/constraint-mode fields, higher-order Lagrange
//! resampling, saving through a mask, and parallel halo handling. See
//! [`VtuOptions::unsupported`].
//!
//! [`MeshView`]: elmer_core::MeshView
//! [`FieldSource`]: elmer_core::FieldSource

mod celltype;
mod collection;
mod fields;
mod fortran_fmt;
mod naming;
mod options;
mod post_file;
mod vtu;

pub use celltype::{
    elmer_to_vtk_type, vtk_cell_type, vtk_node_order, written_node_count,
};
pub use collection::{
    pvd_time, render_pvd, render_pvtu, write_pvd, write_pvtu, PvdEntry,
};
pub use naming::{vtu_file_name, VtuNaming};
pub use options::{FileIndexer, GeometryIds, VtuOptions};
pub use post_file::{
    post_file_solver, result_output_formats, result_output_prefix, OutputFormats,
};
pub use vtu::{render_vtu, write_vtu, write_vtu_group};
