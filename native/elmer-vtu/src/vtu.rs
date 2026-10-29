//! `WriteVtuFile` (VtuOutputSolver.F90:715-1990): one `.vtu` piece.
//!
//! The byte layout is the Fortran one, including its quirks: ASCII reals are
//! `ES16.7E3` (`ES12.3E3` single precision) with `" 0.0"` for tiny values,
//! the ASCII `types` array is declared with an upper-case `FORMAT=`, and
//! binary data is one raw `AppendedData` block of `UInt32`-length-prefixed
//! arrays whose `offset`s the XML part records.

use std::io;
use std::path::Path;

use elmer_core::{ElemId, FieldKind, FieldSource, MeshView};

use crate::celltype::{vtk_cell_type, vtk_node_order, written_node_count};
use crate::fields::{self, invalid, unsupported, Selected};
use crate::fortran_fmt::ascii_real;
use crate::VtuOptions;

enum Data {
    Real(Vec<f64>),
    Int(Vec<i32>),
}

impl Data {
    fn len(&self) -> usize {
        match self {
            Data::Real(v) => v.len(),
            Data::Int(v) => v.len(),
        }
    }
}

struct Array {
    name: Option<String>,
    ncomp: Option<usize>,
    data: Data,
    /// The `types` array is written with `FORMAT="ascii"` in ASCII mode.
    shouting_format: bool,
}

impl Array {
    fn real(name: Option<String>, ncomp: Option<usize>, v: Vec<f64>) -> Self {
        Array { name, ncomp, data: Data::Real(v), shouting_format: false }
    }
    fn int(name: &str, v: Vec<i32>) -> Self {
        Array { name: Some(name.to_string()), ncomp: None, data: Data::Int(v), shouting_format: false }
    }
}

struct Emitter<'o> {
    opts: &'o VtuOptions,
    xml: Vec<u8>,
    /// Appended payload (binary mode), each array preceded by its byte count.
    appended: Vec<u8>,
}

impl Emitter<'_> {
    fn put(&mut self, s: &str) {
        self.xml.extend_from_slice(s.as_bytes());
    }

    fn array(&mut self, a: &Array) -> io::Result<()> {
        let ty = match a.data {
            Data::Real(_) => self.opts.real_type(),
            Data::Int(_) => "Int32",
        };
        let mut tag = format!("        <DataArray type=\"{ty}\"");
        if let Some(n) = &a.name {
            tag.push_str(&format!(" Name=\"{n}\""));
        }
        if let Some(c) = a.ncomp {
            tag.push_str(&format!(" NumberOfComponents=\"{c}\""));
        }
        if self.opts.ascii {
            tag.push_str(if a.shouting_format { " FORMAT=\"ascii\">\n" } else { " format=\"ascii\">\n" });
            self.put(&tag);
            let single = self.opts.single_precision;
            let mut body = String::new();
            match &a.data {
                Data::Real(v) => v.iter().for_each(|x| body.push_str(&ascii_real(*x, single))),
                Data::Int(v) => v.iter().for_each(|i| {
                    body.push(' ');
                    body.push_str(&i.to_string());
                }),
            }
            self.put(&body);
            self.put("\n        </DataArray>\n");
        } else {
            let offset = self.appended.len();
            tag.push_str(&format!(" format=\"appended\" offset=\"{offset}\"/>\n"));
            self.put(&tag);
            let size = match &a.data {
                Data::Real(v) => v.len() * self.opts.real_size(),
                Data::Int(v) => v.len() * 4,
            };
            let size = u32::try_from(size)
                .map_err(|_| invalid(format!("array larger than 4 GiB: {:?}", a.name)))?;
            self.appended.extend_from_slice(&size.to_le_bytes());
            match &a.data {
                Data::Real(v) if self.opts.single_precision => {
                    v.iter().for_each(|x| self.appended.extend_from_slice(&(*x as f32).to_le_bytes()))
                }
                Data::Real(v) => v.iter().for_each(|x| self.appended.extend_from_slice(&x.to_le_bytes())),
                Data::Int(v) => v.iter().for_each(|i| self.appended.extend_from_slice(&i.to_le_bytes())),
            }
        }
        Ok(())
    }
}

/// `GeometryBodyMap` / `GeometryBCMap` and the BC offset
/// (VtuOutputSolver.F90:515-585).
struct IdMaps {
    body: Vec<i32>,
    bc: Vec<i32>,
    bc_offset: i32,
    bc_tags: Vec<i32>,
}

fn id_maps(opts: &VtuOptions, first_active_is_bulk: bool) -> io::Result<IdMaps> {
    let g = opts.geometry.as_ref().ok_or_else(|| {
        invalid("Save Geometry Ids needs VtuOptions::geometry (body and BC counts, BC tags)")
    })?;
    let mut body: Vec<i32> = match opts.default_body_id {
        Some(j) => vec![j; g.n_bodies],
        None => (1..=g.n_bodies as i32).collect(),
    };
    for (i, id) in g.body_geometry_id.iter().enumerate().take(g.n_bodies) {
        if let Some(j) = id {
            body[i] = *j;
        }
    }

    let mut bc_offset = 0;
    let mut bc: Vec<i32>;
    match opts.default_bc_id {
        Some(j) => bc = vec![j; g.n_bcs()],
        None => {
            bc_offset = match opts.bc_id_offset {
                Some(o) => o,
                None if first_active_is_bulk => {
                    let mut o = 100;
                    while o <= g.n_bodies as i32 {
                        o *= 10;
                    }
                    o
                }
                None => 0,
            };
            bc = (1..=g.n_bcs() as i32).map(|i| i + bc_offset).collect();
        }
    }
    for (i, id) in g.bc_geometry_id.iter().enumerate().take(g.n_bcs()) {
        if let Some(j) = id {
            bc[i] = *j;
        }
    }
    Ok(IdMaps { body, bc, bc_offset, bc_tags: g.bc_tags.clone() })
}

impl IdMaps {
    fn body_id(&self, tag: i32) -> i32 {
        if tag >= 1 && (tag as usize) <= self.body.len() {
            self.body[tag as usize - 1]
        } else {
            tag
        }
    }

    /// `GetBCId` then the map; unmatched tags get the offset.
    fn bc_id(&self, tag: i32) -> i32 {
        match self.bc_tags.iter().position(|t| *t == tag) {
            Some(i) if i < self.bc.len() => self.bc[i],
            _ => self.bc_offset,
        }
    }
}

/// Which elements a file contains (`GenerateSaveMask`).
pub(crate) struct Selection {
    pub elems: Vec<ElemId>,
    /// Saved nodes in file order (ascending global index).
    pub nodes: Vec<u32>,
    /// Global node to file node; `None` when every node is saved in order.
    pub node_map: Option<Vec<u32>>,
}

pub(crate) fn select(
    mesh: &dyn MeshView,
    opts: &VtuOptions,
    group: Option<usize>,
) -> io::Result<Selection> {
    let n_bulk = mesh.n_bulk();
    let total = n_bulk + mesh.n_boundary();
    let n_bodies = opts.geometry.as_ref().map_or(0, |g| g.n_bodies);
    if group.is_some() && opts.geometry.is_none() {
        return Err(invalid("Vtu Part Collection needs VtuOptions::geometry"));
    }

    let mut elems = Vec::new();
    let mut used = vec![false; mesh.n_nodes()];
    for e in 0..total {
        let boundary = e >= n_bulk;
        if boundary && opts.save_bulk_only || !boundary && opts.save_boundaries_only {
            continue;
        }
        let el = mesh.element(ElemId(e as u32));
        if let Some(g) = group {
            if !boundary {
                if el.tag as i64 != g as i64 {
                    continue;
                }
            } else if el.tag as i64 != g as i64 - n_bodies as i64 {
                continue;
            }
        }
        if el.etype.0 < 200 {
            continue;
        }
        let m = if opts.save_linear {
            written_node_count(el.etype, true)
        } else {
            el.nodes.len()
        };
        if m > el.nodes.len() {
            return Err(invalid(format!("element {} has fewer nodes than its type", e + 1)));
        }
        for &n in &el.nodes[..m] {
            *used
                .get_mut(n as usize)
                .ok_or_else(|| invalid(format!("element {} references node {} past the mesh", e + 1, n + 1)))? = true;
        }
        elems.push(ElemId(e as u32));
    }

    let nodes: Vec<u32> = (0..used.len() as u32).filter(|&i| used[i as usize]).collect();
    let node_map = if nodes.len() == mesh.n_nodes() {
        None
    } else {
        let mut map = vec![u32::MAX; mesh.n_nodes()];
        for (new, &old) in nodes.iter().enumerate() {
            map[old as usize] = new as u32;
        }
        Some(map)
    };
    Ok(Selection { elems, nodes, node_map })
}

fn check_supported(opts: &VtuOptions) -> io::Result<()> {
    if opts.unsupported.is_empty() {
        Ok(())
    } else {
        Err(unsupported(format!(
            "keywords not supported by elmer-vtu: {}",
            opts.unsupported.join(", ")
        )))
    }
}

/// One `.vtu` as bytes, or `None` when nothing is selected (Fortran then
/// writes no file; `NumberOfDofNodes == 0`).
///
/// `group` is the 1-based `GroupId` of `Vtu Part Collection`.
pub fn render_vtu(
    mesh: &dyn MeshView,
    fields: &dyn FieldSource,
    opts: &VtuOptions,
    group: Option<usize>,
) -> io::Result<Option<Vec<u8>>> {
    check_supported(opts)?;
    let sel = select(mesh, opts, group)?;
    if sel.nodes.is_empty() {
        return Ok(None);
    }
    let field_refs = fields.fields();
    let plan = fields::plan(&field_refs, opts)?;
    let n_nodes = sel.nodes.len();
    let n_cells = sel.elems.len();

    // Point data.
    let mut point = Vec::new();
    for s in &plan.point {
        let mut v = Vec::with_capacity(n_nodes * s.sdofs);
        let mut buf = vec![0.0; s.sdofs];
        for &n in &sel.nodes {
            s.nodal(n as usize, &mut buf)?;
            v.extend_from_slice(&buf);
        }
        point.push(Array::real(Some(s.name.clone()), Some(s.sdofs), v));
    }

    // Cell data.
    let mut cell = Vec::new();
    for s in &plan.cell {
        cell.push(Array::real(
            Some(s.name.clone()),
            Some(s.sdofs),
            elemental_values(mesh, &sel, s)?,
        ));
    }
    if opts.save_geometry_ids {
        let first_is_bulk = sel.elems.first().is_some_and(|e| (e.0 as usize) < mesh.n_bulk());
        let maps = id_maps(opts, first_is_bulk)?;
        let ids = sel
            .elems
            .iter()
            .map(|e| {
                let el = mesh.element(*e);
                if (e.0 as usize) < mesh.n_bulk() {
                    maps.body_id(el.tag)
                } else {
                    maps.bc_id(el.tag)
                }
            })
            .collect();
        cell.push(Array::int("GeometryIds", ids));
    }

    // Points.
    let points = Array::real(None, Some(3), coordinates(mesh, &field_refs, &sel, opts)?);

    // Cells.
    let mut conn = Vec::new();
    let mut offsets = Vec::with_capacity(n_cells);
    let mut types = Vec::with_capacity(n_cells);
    for &e in &sel.elems {
        let el = mesh.element(e);
        let ty = vtk_cell_type(el.etype, opts.save_linear).ok_or_else(|| {
            invalid(format!("Not implemented for elementtype: {}", el.etype.0))
        })?;
        let n = if opts.save_linear { written_node_count(el.etype, true) } else { el.nodes.len() };
        let order = if opts.save_linear { &[][..] } else { vtk_node_order(el.etype) };
        for j in 0..n {
            let local = if order.is_empty() { j } else { order[j] };
            let node = *el.nodes.get(local).ok_or_else(|| {
                invalid(format!("element {} has no local node {}", e.0 + 1, local + 1))
            })?;
            conn.push(match &sel.node_map {
                None => node as i32,
                Some(m) => m[node as usize] as i32,
            });
        }
        offsets.push(conn.len() as i32);
        types.push(ty as i32);
    }
    let conn = Array::int("connectivity", conn);
    let offsets = Array::int("offsets", offsets);
    let mut types = Array::int("types", types);
    types.shouting_format = true;

    // Assemble.
    let mut em = Emitter { opts, xml: Vec::new(), appended: Vec::new() };
    em.put("<?xml version=\"1.0\"?>\n");
    em.put("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    if !opts.skip_metainfo {
        for line in &opts.metainfo {
            em.put(&format!("<!-- {line} -->\n"));
        }
    }
    em.put("  <UnstructuredGrid>\n");
    em.put(&format!(
        "    <Piece NumberOfPoints=\"{n_nodes}\" NumberOfCells=\"{n_cells}\">\n"
    ));
    em.put("      <PointData>\n");
    for a in &point {
        em.array(a)?;
    }
    em.put("      </PointData>\n");
    em.put("      <CellData>\n");
    for a in &cell {
        em.array(a)?;
    }
    em.put("      </CellData>\n");
    em.put("      <Points>\n");
    em.array(&points)?;
    em.put("      </Points>\n");
    em.put("      <Cells>\n");
    em.array(&conn)?;
    em.array(&offsets)?;
    em.array(&types)?;
    em.put("      </Cells>\n");
    em.put("    </Piece>\n");
    em.put("  </UnstructuredGrid>\n");
    let mut out = em.xml;
    if !opts.ascii {
        out.extend_from_slice(b"<AppendedData encoding=\"raw\">\n_");
        out.extend_from_slice(&em.appended);
        out.extend_from_slice(b"\n</AppendedData>\n");
    }
    out.extend_from_slice(b"</VTKFile>\n");
    debug_assert!(point.iter().all(|a| a.data.len() % a.ncomp.unwrap_or(1) == 0));
    Ok(Some(out))
}

fn elemental_values(mesh: &dyn MeshView, sel: &Selection, s: &Selected<'_>) -> io::Result<Vec<f64>> {
    let mut v = Vec::with_capacity(sel.elems.len() * s.sdofs);
    let mut buf = vec![0.0; s.sdofs];
    for &e in &sel.elems {
        let el = mesh.element(e);
        s.elemental(e, &el, &mut buf)?;
        v.extend_from_slice(&buf);
    }
    Ok(v)
}

/// Node coordinates, `x = scale * (x - displacement) + offset`.
fn coordinates(
    mesh: &dyn MeshView,
    fields: &[elmer_core::FieldRef<'_>],
    sel: &Selection,
    opts: &VtuOptions,
) -> io::Result<Vec<f64>> {
    let find = |name: &str| {
        fields
            .iter()
            .find(|f| f.name.trim().eq_ignore_ascii_case(name) && f.kind == FieldKind::Nodal)
    };
    let (mut disp, mut update) = (None, None);
    if opts.fixed_mesh {
        if opts.displace_mesh {
            disp = find("displacement");
        }
        update = find("mesh update");
    }

    let mut v = Vec::with_capacity(sel.nodes.len() * 3);
    for &n in &sel.nodes {
        let mut c = mesh.coord(n);
        let mut done = false;
        if let Some(d) = disp {
            if let Some(slot) = fields::node_slot(d, n as usize)? {
                subtract(d, slot, &mut c)?;
                done = true;
            }
        }
        if let (false, Some(u)) = (done, update) {
            if let Some(slot) = fields::node_slot(u, n as usize)? {
                subtract(u, slot, &mut c)?;
            }
        }
        for i in 0..3 {
            v.push(opts.coordinate_scale[i] * c[i] + opts.coordinate_offset[i]);
        }
    }
    Ok(v)
}

fn subtract(f: &elmer_core::FieldRef<'_>, slot: usize, c: &mut [f64; 3]) -> io::Result<()> {
    // Fortran subtracts x and y always, z only for a 3-component field.
    for (k, ck) in c.iter_mut().enumerate().take(f.dofs.min(3)) {
        let idx = f.dofs * slot + k;
        *ck -= f
            .values
            .get(idx)
            .copied()
            .ok_or_else(|| invalid(format!("value index {idx} out of range for {}", f.name)))?;
    }
    Ok(())
}

/// Write one `.vtu` piece to `path`. Returns `Ok(())` without creating the
/// file when the selection holds no nodes, as `VtuOutputSolver` does.
pub fn write_vtu(
    path: &Path,
    mesh: &dyn MeshView,
    fields: &dyn FieldSource,
    opts: &VtuOptions,
) -> io::Result<()> {
    if let Some(bytes) = render_vtu(mesh, fields, opts, None)? {
        std::fs::write(path, bytes)?;
    }
    Ok(())
}

/// [`write_vtu`] for one group of a `Vtu Part Collection` (1-based `GroupId`).
pub fn write_vtu_group(
    path: &Path,
    mesh: &dyn MeshView,
    fields: &dyn FieldSource,
    opts: &VtuOptions,
    group: usize,
) -> io::Result<()> {
    if let Some(bytes) = render_vtu(mesh, fields, opts, Some(group))? {
        std::fs::write(path, bytes)?;
    }
    Ok(())
}
