//! ElmerPost `.ep` output (`WritePostFile` in `ModelDescription.F90`, the code
//! behind `Post File = case.ep`).
//!
//! Layout: a header line (`nodes elements dofs timesteps` plus one
//! `scalar:`/`vector:` label per output quantity), the `#File started at:`
//! stamp, node coordinates, a `#group all` block of bulk then boundary
//! elements, and one `#time` block of per-node values per save.
//!
//! Which fields appear follows the Fortran variable loop:
//!
//! * only nodal fields with more values than dofs (no globals);
//! * fields named `coordinate 1..3` are skipped (they are written only when
//!   [`EpOptions::save_coordinates`] is set, as `vector: Coordinates`);
//! * `flow solution` (`vector: Velocity scalar: Pressure`, 4 columns),
//!   `displacement`, `mesh update` (only when there is no `displacement`) and
//!   `magnetic field` (3 columns each) are special;
//! * other vectors need a `name[label:k label:k ...]` layout: `k = 1` is a
//!   scalar, `k = 2, 3` a vector, larger `k` scalars `label.1 ...`;
//!   a vector without brackets is not written;
//! * a node outside a field's permutation prints zeros.
//!
//! Names have spaces replaced by `.` and a lowercase first letter capitalised.
//! Complex fields and `Post File Mask Variable` are not supported.
//!
//! Differences from Fortran: the `#time` block always uses the save count for
//! both the count and the step columns (as `WritePostFile` does when
//! appending); `Displacement` coordinates correction uses each variable's own
//! dof count.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use elmer_core::field::{FieldKind, FieldRef, NO_DOF};
use elmer_core::{ElemId, FieldSource, MeshView};

use crate::error::{Result, SaveIoError};
use crate::fortran::{fortran_es_e, fortran_i, now_elmer_date};

#[derive(Clone, Debug, Default)]
pub struct EpOptions {
    /// Body id (as numbered in the model) -> `Name`. Absent: `body<id>`.
    pub body_names: HashMap<i32, String>,
    /// Mesh body tag -> body id; identity when absent.
    pub body_of_tag: HashMap<i32, i32>,
    /// Constraint (BC id) -> `Name`. Absent: `Constraint<id>`.
    pub bc_names: HashMap<i32, String>,
    /// Mesh boundary tag -> constraint (BC id); a tag not listed is
    /// constraint 0 when this map is non-empty, identity when it is empty.
    pub boundary_constraint: HashMap<i32, i32>,
    /// A free-surface BC is present: write `vector: Coordinates` from the
    /// fields `coordinate 1..3`.
    pub save_coordinates: bool,
    /// Scale of the displacement subtracted from the node coordinates
    /// (`MeshScale`): 0 when a solver has `Displace Mesh = False`, otherwise 1.
    pub mesh_scale: Option<f64>,
    /// `!File started at` text; the current UTC time if `None`.
    pub started_at: Option<String>,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
enum Kind {
    MeshUpdate,
    Displacement,
    Flow,
    Magnetic,
    Coordinate,
    Scalar,
    Bracket,
    Skip,
}

fn base_name(name: &str) -> &str {
    name.split('[').next().unwrap_or(name).trim_end()
}

fn classify(f: &FieldRef<'_>) -> Kind {
    if f.kind != FieldKind::Nodal || f.values.len() == f.dofs {
        return Kind::Skip;
    }
    match base_name(f.name) {
        "mesh update" => Kind::MeshUpdate,
        "mesh update 1" | "mesh update 2" | "mesh update 3" => Kind::Skip,
        "displacement" => Kind::Displacement,
        "displacement 1" | "displacement 2" | "displacement 3" => Kind::Skip,
        "flow solution" => Kind::Flow,
        "velocity 1" | "velocity 2" | "velocity 3" | "pressure" => Kind::Skip,
        "magnetic field" => Kind::Magnetic,
        "magnetic field 1" | "magnetic field 2" | "magnetic field 3" => Kind::Skip,
        "coordinate 1" | "coordinate 2" | "coordinate 3" => Kind::Coordinate,
        _ if f.dofs == 1 => Kind::Scalar,
        _ if f.name.contains('[') => Kind::Bracket,
        _ => Kind::Skip,
    }
}

fn dotted(label: &str) -> String {
    let mut s: String = label
        .chars()
        .map(|c| if c == ' ' { '.' } else { c })
        .collect();
    if s.as_bytes().first().is_some_and(u8::is_ascii_lowercase) {
        s[..1].make_ascii_uppercase();
    }
    s
}

/// `(label, k)` pairs of a `name[a:2 b:1]` field name.
fn bracket_parts(name: &str) -> Result<Vec<(String, usize)>> {
    let bad = || SaveIoError::Input(format!("syntax error in variable definition '{name}'"));
    let open = name.find('[').ok_or_else(bad)?;
    let mut rest = &name[open + 1..];
    let mut out = Vec::new();
    while let Some(colon) = rest.find(':') {
        let k = rest[colon + 1..]
            .chars()
            .next()
            .and_then(|c| c.to_digit(10))
            .ok_or_else(bad)? as usize;
        out.push((dotted(rest[..colon].trim_start()), k));
        rest = rest[colon + 2..].trim_start();
    }
    Ok(out)
}

struct Layout<'a> {
    fields: Vec<(Kind, FieldRef<'a>)>,
    has_displacement: bool,
}

fn layout<'a>(fields: Vec<FieldRef<'a>>) -> Layout<'a> {
    let fields: Vec<_> = fields
        .into_iter()
        .map(|f| (classify(&f), f))
        .filter(|(k, _)| *k != Kind::Skip)
        .collect();
    let has_displacement = fields.iter().any(|(k, _)| *k == Kind::Displacement);
    Layout {
        fields,
        has_displacement,
    }
}

impl Layout<'_> {
    /// Fields that produce columns (coordinates are handled separately).
    fn columns(&self) -> impl Iterator<Item = &(Kind, FieldRef<'_>)> {
        self.fields.iter().filter(move |(k, _)| match k {
            Kind::Coordinate => false,
            Kind::MeshUpdate => !self.has_displacement,
            _ => true,
        })
    }

    fn find(&self, base: &str) -> Option<&FieldRef<'_>> {
        self.fields
            .iter()
            .map(|(_, f)| f)
            .find(|f| base_name(f.name) == base)
    }
}

fn header_dofs(l: &Layout<'_>, save_coords: bool) -> Result<usize> {
    let mut n = 0;
    for (k, f) in l.columns() {
        n += match k {
            Kind::MeshUpdate | Kind::Displacement | Kind::Magnetic => 3,
            Kind::Flow => 4,
            Kind::Scalar => 1,
            Kind::Bracket => bracket_parts(f.name)?
                .iter()
                .map(|&(_, k)| k + usize::from(k == 2))
                .sum(),
            Kind::Coordinate | Kind::Skip => 0,
        };
    }
    Ok(n + if save_coords { 3 } else { 0 })
}

fn header_labels(l: &Layout<'_>, save_coords: bool) -> Result<String> {
    let mut s = String::new();
    for (k, f) in l.columns() {
        match k {
            Kind::MeshUpdate => s.push_str(" vector: Mesh.Update"),
            Kind::Displacement => s.push_str(" vector: Displacement"),
            Kind::Flow => s.push_str(" vector: Velocity scalar: Pressure"),
            Kind::Magnetic => s.push_str(" vector: MagField"),
            Kind::Scalar => {
                s.push_str(" scalar: ");
                s.push_str(&dotted(base_name(f.name)));
            }
            Kind::Bracket => {
                for (label, k) in bracket_parts(f.name)? {
                    match k {
                        1 => s.push_str(&format!(" scalar: {label}")),
                        2 | 3 => s.push_str(&format!(" vector: {label}")),
                        _ => (1..=k).for_each(|i| s.push_str(&format!(" scalar: {label}.{i}"))),
                    }
                }
            }
            Kind::Coordinate | Kind::Skip => {}
        }
    }
    if save_coords {
        s.push_str(" vector: Coordinates");
    }
    Ok(s)
}

fn es(x: f64) -> String {
    fortran_es_e(x, 17, 8, Some(3))
}

/// Value of component `d` of `f` at `node`, `None` outside the permutation.
fn val(f: &FieldRef<'_>, node: usize, d: usize) -> Option<f64> {
    let slot = match f.perm {
        Some(p) => {
            let s = *p.get(node)?;
            if s == NO_DOF {
                return None;
            }
            s as usize
        }
        None => node,
    };
    f.values.get(slot * f.dofs + d).copied()
}

fn push_vec(row: &mut String, f: &FieldRef<'_>, node: usize, upto: usize) {
    for d in 0..upto {
        row.push_str(&es(val(f, node, d).unwrap_or(0.0)));
    }
    if f.dofs == 2 {
        row.push_str(" 0.0");
    }
}

fn node_row(l: &Layout<'_>, node: usize, save_coords: bool) -> Result<String> {
    let mut row = String::new();
    for (k, f) in l.columns() {
        let present = val(f, node, 0).is_some();
        match k {
            Kind::MeshUpdate | Kind::Magnetic => {
                if present {
                    push_vec(&mut row, f, node, f.dofs);
                } else {
                    row.push_str(" 0.0 0.0 0.0");
                }
            }
            Kind::Displacement => {
                if present {
                    push_vec(&mut row, f, node, f.dofs);
                } else if let Some(mu) = l.find("mesh update").filter(|m| val(m, node, 0).is_some())
                {
                    push_vec(&mut row, mu, node, mu.dofs);
                } else {
                    row.push_str(" 0.0 0.0 0.0");
                }
            }
            Kind::Flow => {
                if present {
                    for d in 0..f.dofs - 1 {
                        row.push_str(&es(val(f, node, d).unwrap_or(0.0)));
                    }
                    if f.dofs < 4 {
                        row.push_str(" 0.0");
                    }
                    row.push_str(&es(val(f, node, f.dofs - 1).unwrap_or(0.0)));
                } else {
                    row.push_str(" 0.0 0.0 0.0 0.0");
                }
            }
            Kind::Scalar => match val(f, node, 0) {
                Some(v) => row.push_str(&es(v)),
                None => row.push_str(" 0.0"),
            },
            Kind::Bracket => {
                let mut done = 0;
                for (_, q) in bracket_parts(f.name)? {
                    if present {
                        for d in done..done + q {
                            row.push_str(&es(val(f, node, d).unwrap_or(0.0)));
                        }
                        if q == 2 {
                            row.push_str(" 0.0");
                        }
                    } else {
                        for _ in 0..if q == 2 { 3 } else { q } {
                            row.push_str(" 0.0");
                        }
                    }
                    done += q;
                }
            }
            Kind::Coordinate | Kind::Skip => {}
        }
    }
    if save_coords {
        for c in 1..=3 {
            let f = l
                .fields
                .iter()
                .map(|(_, f)| f)
                .find(|f| f.name == format!("coordinate {c}"))
                .ok_or_else(|| {
                    SaveIoError::Input(format!("save_coordinates needs the field 'coordinate {c}'"))
                })?;
            row.push_str(&es(f.values.get(node).copied().unwrap_or(0.0)));
        }
    }
    Ok(row)
}

fn coordinates(mesh: &dyn MeshView, l: &Layout<'_>, scale: f64) -> String {
    let mut s = String::new();
    for node in 0..mesh.n_nodes() {
        let mut c = mesh.coord(node as u32);
        if let Some(disp) = l.find("displacement") {
            if val(disp, node, 0).is_some() {
                for (d, cd) in c.iter_mut().enumerate().take(disp.dofs.min(3)) {
                    *cd -= scale * val(disp, node, d).unwrap_or(0.0);
                }
            } else if let Some(mu) = l.find("mesh update").filter(|m| val(m, node, 0).is_some()) {
                for (d, cd) in c.iter_mut().enumerate().take(mu.dofs.min(3)) {
                    *cd -= scale * val(mu, node, d).unwrap_or(0.0);
                }
            }
        }
        match mesh.dim() {
            3 => s.push_str(&format!("{}{}{}\n", es(c[0]), es(c[1]), es(c[2]))),
            2 => s.push_str(&format!("{}{} 0.0\n", es(c[0]), es(c[1]))),
            _ => s.push_str(&format!("{} 0.0 0.0\n", es(c[0]))),
        }
    }
    s
}

fn elements(mesh: &dyn MeshView, o: &EpOptions) -> String {
    let mut s = String::from("#group all\n");
    let n_bulk = mesh.n_bulk();
    for i in 0..n_bulk + mesh.n_boundary() {
        let e = mesh.element(ElemId(i as u32));
        let nodes = e.nodes;
        if i < n_bulk {
            let id = o.body_of_tag.get(&e.tag).copied().unwrap_or(e.tag);
            match o.body_names.get(&id) {
                Some(n) => s.push_str(&n.replace(' ', ".")),
                None => s.push_str(&format!(
                    "body{} ",
                    fortran_i(
                        id.into(),
                        if (1..100).contains(&id) {
                            if id < 10 {
                                1
                            } else {
                                2
                            }
                        } else {
                            3
                        }
                    )
                )),
            }
            s.push_str(&fortran_i(e.etype.0.into(), 5));
            for chunk in nodes.chunks(4) {
                for n in chunk {
                    s.push_str(&fortran_i(i64::from(*n), 8));
                }
                s.push('\n');
            }
        } else {
            let k = if o.boundary_constraint.is_empty() {
                e.tag
            } else {
                o.boundary_constraint.get(&e.tag).copied().unwrap_or(0)
            };
            match o.bc_names.get(&k) {
                Some(n) => s.push_str(&n.replace(' ', ".")),
                None => s.push_str(&format!(
                    "Constraint{} ",
                    fortran_i(
                        k.into(),
                        if (0..10).contains(&k) {
                            1
                        } else if (10..100).contains(&k) {
                            2
                        } else {
                            3
                        }
                    )
                )),
            }
            s.push_str(&fortran_i(e.etype.0.into(), 5));
            for n in nodes {
                s.push_str(&fortran_i(i64::from(*n), 8));
            }
            s.push('\n');
        }
    }
    s.push_str("#endgroup all\n");
    s
}

/// Appends saves to one `.ep` file.
pub struct EpWriter {
    out: BufWriter<File>,
    opts: EpOptions,
    saves: usize,
}

impl EpWriter {
    /// Create the file and write the header and mesh. `first` fixes the
    /// variable layout; `time_count` is the number of time steps the file will hold.
    pub fn create(
        path: &Path,
        opts: EpOptions,
        mesh: &dyn MeshView,
        first: &dyn FieldSource,
        time_count: usize,
    ) -> Result<Self> {
        let fields = first.fields();
        let l = layout(fields);
        let n_elem = mesh.n_bulk() + mesh.n_boundary();
        let mut out = BufWriter::new(File::create(path)?);
        let labels = header_labels(&l, opts.save_coordinates)?;
        let dofs = header_dofs(&l, opts.save_coordinates)?;
        writeln!(
            out,
            "{}{}{}{}{labels}",
            fortran_i(mesh.n_nodes() as i64, 10),
            fortran_i(n_elem as i64, 10),
            fortran_i(dofs as i64, 7),
            fortran_i(time_count as i64, 7),
        )?;
        let stamp = opts.started_at.clone().unwrap_or_else(now_elmer_date);
        writeln!(out, "#File started at: {stamp}")?;
        out.write_all(coordinates(mesh, &l, opts.mesh_scale.unwrap_or(1.0)).as_bytes())?;
        out.write_all(elements(mesh, &opts).as_bytes())?;
        Ok(EpWriter {
            out,
            opts,
            saves: 0,
        })
    }

    /// Append one `#time` block.
    pub fn append(&mut self, mesh: &dyn MeshView, src: &dyn FieldSource) -> Result<()> {
        let l = layout(src.fields());
        self.saves += 1;
        let n = fortran_i(self.saves as i64, 7);
        writeln!(self.out, "#time {n}{n}{}", es(src.time()))?;
        for node in 0..mesh.n_nodes() {
            writeln!(
                self.out,
                "{}",
                node_row(&l, node, self.opts.save_coordinates)?
            )?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.out.flush()?;
        Ok(())
    }
}

/// Write a whole `.ep` file: header and one `#time` block per snapshot.
pub fn write_ep(
    path: &Path,
    mesh: &dyn MeshView,
    snapshots: &[&dyn FieldSource],
    opts: EpOptions,
) -> Result<()> {
    let first = snapshots
        .first()
        .ok_or_else(|| SaveIoError::Input("write_ep needs at least one snapshot".into()))?;
    let mut w = EpWriter::create(path, opts, mesh, *first, snapshots.len())?;
    for s in snapshots {
        w.append(mesh, *s)?;
    }
    w.finish()
}
