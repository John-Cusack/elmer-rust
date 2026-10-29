//! Test-only readers: Elmer mesh files and `.vtu` files.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

use elmer_core::{ElemId, ElementRef, ElementType, FieldKind, FieldRef, FieldSource, MeshView};

pub fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
}

pub fn fixture_bytes(name: &str) -> Vec<u8> {
    std::fs::read(fixture(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// First differing byte, for readable failures on big buffers.
pub fn assert_bytes_eq(got: &[u8], want: &[u8], what: &str) {
    if got != want {
        let i = got.iter().zip(want).position(|(a, b)| a != b).unwrap_or(got.len().min(want.len()));
        let lo = i.saturating_sub(40);
        panic!(
            "{what}: differs at byte {i} (got {} bytes, want {})\n got: {:?}\nwant: {:?}",
            got.len(),
            want.len(),
            String::from_utf8_lossy(&got[lo..(i + 40).min(got.len())]),
            String::from_utf8_lossy(&want[lo..(i + 40).min(want.len())]),
        );
    }
}

// ---------------------------------------------------------------- mesh files

pub struct TestMesh {
    pub coords: Vec<[f64; 3]>,
    /// Bulk elements, then boundary elements: (type, 0-based nodes, tag).
    pub elems: Vec<(ElementType, Vec<u32>, i32)>,
    pub n_bulk: usize,
}

impl TestMesh {
    /// Read `mesh.header/nodes/elements/boundary` from `fixtures/meshes/<dir>`.
    pub fn read(dir: &str) -> TestMesh {
        let d = fixture("meshes").join(dir);
        let text = |f: &str| std::fs::read_to_string(d.join(f)).unwrap();
        let nums = |l: &str| l.split_whitespace().map(str::to_string).collect::<Vec<_>>();

        let coords = text("mesh.nodes")
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let t = nums(l);
                [t[2].parse().unwrap(), t[3].parse().unwrap(), t[4].parse().unwrap()]
            })
            .collect();

        let mut elems = Vec::new();
        for l in text("mesh.elements").lines().filter(|l| !l.trim().is_empty()) {
            let t = nums(l);
            let ty: u16 = t[2].parse().unwrap();
            let nodes = t[3..].iter().map(|n| n.parse::<u32>().unwrap() - 1).collect();
            elems.push((ElementType(ty), nodes, t[1].parse().unwrap()));
        }
        let n_bulk = elems.len();
        for l in text("mesh.boundary").lines().filter(|l| !l.trim().is_empty()) {
            let t = nums(l);
            let ty: u16 = t[4].parse().unwrap();
            let nodes = t[5..].iter().map(|n| n.parse::<u32>().unwrap() - 1).collect();
            elems.push((ElementType(ty), nodes, t[1].parse().unwrap()));
        }
        TestMesh { coords, elems, n_bulk }
    }
}

impl MeshView for TestMesh {
    fn dim(&self) -> usize {
        3
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let (etype, nodes, tag) = &self.elems[e.0 as usize];
        ElementRef { etype: *etype, nodes, tag: *tag, parents: [None, None] }
    }
}

// ---------------------------------------------------------------- fields

pub struct OwnedField {
    pub name: String,
    pub kind: FieldKind,
    pub dofs: usize,
    pub perm: Option<Vec<u32>>,
    pub values: Vec<f64>,
}

impl OwnedField {
    pub fn nodal(name: &str, dofs: usize, values: Vec<f64>) -> Self {
        OwnedField { name: name.into(), kind: FieldKind::Nodal, dofs, perm: None, values }
    }
}

pub struct TestFields {
    pub fields: Vec<OwnedField>,
    pub time: f64,
    pub step: usize,
}

impl TestFields {
    pub fn new(fields: Vec<OwnedField>) -> Self {
        TestFields { fields, time: 0.0, step: 1 }
    }

    /// Nodal fields as stored in a parsed `.vtu`'s PointData.
    pub fn from_vtu(v: &Vtu) -> Self {
        TestFields::new(
            v.arrays
                .iter()
                .filter(|a| a.section == "PointData")
                .map(|a| OwnedField::nodal(a.name.as_deref().unwrap(), a.ncomp, a.reals.clone()))
                .collect(),
        )
    }
}

impl FieldSource for TestFields {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        self.fields
            .iter()
            .map(|f| FieldRef {
                name: &f.name,
                kind: f.kind,
                dofs: f.dofs,
                perm: f.perm.as_deref(),
                values: &f.values,
            })
            .collect()
    }
    fn time(&self) -> f64 {
        self.time
    }
    fn step(&self) -> usize {
        self.step
    }
}

// ---------------------------------------------------------------- vtu reader

#[derive(Debug, PartialEq)]
pub struct VArray {
    pub section: String,
    pub name: Option<String>,
    pub ty: String,
    pub ncomp: usize,
    pub reals: Vec<f64>,
    pub ints: Vec<i32>,
}

#[derive(Debug, PartialEq)]
pub struct Vtu {
    pub npoints: usize,
    pub ncells: usize,
    pub binary: bool,
    pub arrays: Vec<VArray>,
}

impl Vtu {
    pub fn array(&self, section: &str, name: Option<&str>) -> &VArray {
        self.arrays
            .iter()
            .find(|a| a.section == section && a.name.as_deref() == name)
            .unwrap_or_else(|| panic!("no {section} array {name:?}"))
    }
}

fn attr(tag: &str, key: &str) -> Option<String> {
    let lower = tag.to_ascii_lowercase();
    let pat = format!("{} =\"", key.to_ascii_lowercase()).replace(" =", "=");
    let i = lower.find(&pat)? + pat.len();
    let j = tag[i..].find('"')? + i;
    Some(tag[i..j].to_string())
}

/// Parse an Elmer-written `.vtu` (ASCII or appended raw).
pub fn parse_vtu(bytes: &[u8]) -> Vtu {
    let marker = b"<AppendedData encoding=\"raw\">\n_";
    let app = bytes.windows(marker.len()).position(|w| w == marker);
    let text_end = app.unwrap_or(bytes.len());
    let text = String::from_utf8_lossy(&bytes[..text_end]).into_owned();
    let base = app.map(|p| p + marker.len());

    let piece = text.find("<Piece").unwrap();
    let ptag = &text[piece..text[piece..].find('>').unwrap() + piece];
    let npoints = attr(ptag, "NumberOfPoints").unwrap().parse().unwrap();
    let ncells = attr(ptag, "NumberOfCells").unwrap().parse().unwrap();

    let mut arrays = Vec::new();
    let mut section = String::new();
    let mut pos = 0;
    while let Some(rel) = text[pos..].find('<') {
        let start = pos + rel;
        let end = text[start..].find('>').unwrap() + start;
        let tag = &text[start..=end];
        pos = end + 1;
        for s in ["PointData", "CellData", "Points", "Cells"] {
            if tag == format!("<{s}>") {
                section = s.to_string();
            }
        }
        if !tag.starts_with("<DataArray") {
            continue;
        }
        let ty = attr(tag, "type").unwrap();
        let ncomp = attr(tag, "NumberOfComponents").map_or(1, |s| s.parse().unwrap());
        let name = attr(tag, "Name");
        let (mut reals, mut ints) = (Vec::new(), Vec::new());
        match attr(tag, "format").unwrap().as_str() {
            "ascii" => {
                let close = text[pos..].find("</DataArray>").unwrap() + pos;
                for t in text[pos..close].split_whitespace() {
                    if ty == "Int32" {
                        ints.push(t.parse().unwrap());
                    } else {
                        reals.push(t.parse().unwrap());
                    }
                }
                pos = close;
            }
            "appended" => {
                let off: usize = attr(tag, "offset").unwrap().parse().unwrap();
                let at = base.unwrap() + off;
                let len = u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
                let data = &bytes[at + 4..at + 4 + len];
                match ty.as_str() {
                    "Int32" => ints = data.chunks(4).map(|c| i32::from_le_bytes(c.try_into().unwrap())).collect(),
                    "Float32" => reals = data.chunks(4).map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64).collect(),
                    "Float64" => reals = data.chunks(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect(),
                    other => panic!("type {other}"),
                }
            }
            other => panic!("format {other}"),
        }
        arrays.push(VArray { section: section.clone(), name, ty, ncomp, reals, ints });
    }
    Vtu { npoints, ncells, binary: app.is_some(), arrays }
}
