//! The wrapper files: `.pvd` time/group collections (`WritePvdFile`,
//! VtuOutputSolver.F90:2183-2240) and `.pvtu` partition wrappers
//! (`WritePvtuFile`, :2252-2590).

use std::io;
use std::path::Path;

use elmer_core::FieldSource;

use crate::fields::{self, invalid};
use crate::fortran_fmt::g0;
use crate::VtuOptions;

const PVD_HEADER: &str =
    "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>";
const PVD_FOOTER: &str = "</Collection></VTKFile>";

/// One `DataSet` of a `.pvd`.
#[derive(Clone, Debug, PartialEq)]
pub struct PvdEntry {
    /// `GetTime()`, or `GetTime() - GetTimestepSize()` under `Vtu time previous`
    /// (see [`pvd_time`]).
    pub timestep: f64,
    /// `GroupId` (0 without `Vtu Part Collection`), written as `part`.
    pub group: usize,
    /// File name relative to the `.pvd` (`.vtu`, or `.pvtu` in a parallel run).
    pub file: String,
}

/// The `timestep` stored for a call at `time` with step size `dt`.
pub fn pvd_time(time: f64, dt: f64, opts: &VtuOptions) -> f64 {
    if opts.time_previous {
        time - dt
    } else {
        time
    }
}

fn dataset(e: &PvdEntry) -> String {
    format!(
        "<DataSet timestep=\"{}\" group=\"\" part=\"{}\" file=\"{}\"/>",
        g0(e.timestep),
        e.group,
        e.file
    )
}

/// Render a `.pvd` the way Fortran's direct-access `REWIND`-free updating
/// leaves it.
///
/// The file is a sequence of fixed-length records, blank padded and without
/// separators: the header, then one record per `DataSet` starting with a line
/// feed, then a closing record starting with a line feed. The record length is
/// `(n/4 + 5) * 4` for `n` the longer of the header and the **first** `DataSet`
/// line, fixed when the file is created. A later `DataSet` that does not fit
/// makes the Fortran run abort with an end-of-record error; here it is an
/// error too.
pub fn render_pvd(entries: &[PvdEntry]) -> io::Result<String> {
    let first = entries.first().ok_or_else(|| invalid("a .pvd needs at least one DataSet"))?;
    let n = PVD_HEADER.len().max(dataset(first).len());
    let rec_len = (n / 4 + 5) * 4;

    let mut out = String::new();
    let mut record = |s: &str| -> io::Result<()> {
        if s.len() > rec_len {
            return Err(invalid(format!(
                "pvd record of {} bytes exceeds the {rec_len}-byte record length fixed by the first entry",
                s.len()
            )));
        }
        out.push_str(s);
        out.extend(std::iter::repeat_n(' ', rec_len - s.len()));
        Ok(())
    };
    record(PVD_HEADER)?;
    for e in entries {
        record(&format!("\n{}", dataset(e)))?;
    }
    record(&format!("\n{PVD_FOOTER}"))?;
    Ok(out)
}

/// Write a `.pvd` collecting `entries` (all calls so far, oldest first).
pub fn write_pvd(path: &Path, entries: &[PvdEntry]) -> io::Result<()> {
    std::fs::write(path, render_pvd(entries)?)
}

/// `.pvtu` text. `pieces` are the `.vtu` names (without directory) of the
/// partitions that wrote data, in partition order.
///
/// Array declarations come from the same field selection as the pieces. As
/// in the Fortran: the ASCII `PDataArray`s are always `Float64`, and the
/// `PPoints` array is spelled `DataArray`.
pub fn render_pvtu(
    fields: &dyn FieldSource,
    opts: &VtuOptions,
    pieces: &[String],
) -> io::Result<String> {
    if !opts.unsupported.is_empty() {
        return Err(fields::unsupported(format!(
            "keywords not supported by elmer-vtu: {}",
            opts.unsupported.join(", ")
        )));
    }
    let refs = fields.fields();
    let plan = fields::plan(&refs, opts)?;
    let scalars = opts.scalar_fields.first();
    let vectors = opts.vector_fields.first();

    let real_type = |ascii: bool| if ascii { "Float64" } else { opts.real_type() };
    let format = if opts.ascii { "ascii" } else { "appended" };
    let parray = |name: &str, sdofs: usize| {
        format!(
            "      <PDataArray type=\"{}\" Name=\"{name}\" NumberOfComponents=\"{sdofs}\" format=\"{format}\"/>\n",
            real_type(opts.ascii)
        )
    };

    let mut s = String::new();
    s.push_str("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    s.push_str("  <PUnstructuredGrid>\n");

    if opts.save_nodal {
        match (scalars, vectors) {
            (Some(a), Some(b)) => {
                s.push_str(&format!("    <PPointData Scalars=\"{a}\" Vectors=\"{b}\">\n"))
            }
            (Some(a), None) => s.push_str(&format!("    <PPointData Scalars=\"{a}\">\n")),
            (None, Some(b)) => s.push_str(&format!("    <PPointData Vectors=\"{b}\">\n")),
            (None, None) => {}
        }
        for f in &plan.point {
            s.push_str(&parray(&f.name, f.sdofs));
        }
        if scalars.is_some() || vectors.is_some() {
            s.push_str("    </PPointData>\n");
        }
    }

    s.push_str("    <PCellData>\n");
    if opts.save_elemental && (scalars.is_some() || vectors.is_some()) {
        for f in &plan.cell {
            s.push_str(&parray(&f.name, f.sdofs));
        }
    }
    if opts.save_geometry_ids {
        s.push_str(&format!(
            "      <PDataArray type=\"Int32\" Name=\"GeometryIds\" format=\"{format}\"/>\n"
        ));
    }
    s.push_str("    </PCellData>\n");

    s.push_str("    <PPoints>\n");
    s.push_str(&format!(
        "      <DataArray type=\"{}\" NumberOfComponents=\"3\" format=\"{format}\"/>\n",
        real_type(opts.ascii)
    ));
    s.push_str("    </PPoints>\n");

    for p in pieces {
        s.push_str(&format!("    <Piece Source=\"{p}\"/>\n"));
    }
    s.push_str("  </PUnstructuredGrid>\n");
    s.push_str("</VTKFile>\n");
    Ok(s)
}

/// Write a `.pvtu`; see [`render_pvtu`].
pub fn write_pvtu(
    path: &Path,
    fields: &dyn FieldSource,
    opts: &VtuOptions,
    pieces: &[String],
) -> io::Result<()> {
    std::fs::write(path, render_pvtu(fields, opts, pieces)?)
}
