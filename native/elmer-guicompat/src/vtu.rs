//! Independent readers for the files the GUI loads after a run.
//!
//! * [`load_vtu`] extracts what `vtkpost.cpp` takes from
//!   `vtkXMLUnstructuredGridReader`: point count, cell count and the point
//!   arrays with their component counts. It also decodes every array and
//!   cross-checks the binary layout, so a malformed writer fails here.
//!   Supported: `appended` raw and `ascii` data, either byte order.
//! * [`load_ep`] reads the ElmerPost header line and `#time` records.

use std::collections::BTreeMap;

use quick_xml::events::Event;
use quick_xml::Reader;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum LoadError {
    #[error("{0}")]
    Format(String),
}

fn bad<T>(msg: impl Into<String>) -> Result<T, LoadError> {
    Err(LoadError::Format(msg.into()))
}

#[derive(Clone, Debug, PartialEq)]
pub struct Array {
    pub name: String,
    pub components: usize,
    pub ty: String,
    /// All values, tuple-major.
    pub values: Vec<f64>,
}

impl Array {
    pub fn tuples(&self) -> usize {
        self.values.len() / self.components.max(1)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Vtu {
    pub points: usize,
    pub cells: usize,
    pub point_data: Vec<Array>,
    pub cell_data: Vec<Array>,
    pub coordinates: Vec<f64>,
    pub connectivity: Vec<i64>,
    pub offsets: Vec<i64>,
    pub types: Vec<i64>,
}

impl Vtu {
    pub fn point_array(&self, name: &str) -> Option<&Array> {
        self.point_data
            .iter()
            .find(|a| a.name.eq_ignore_ascii_case(name))
    }
}

#[derive(Debug)]
struct Decl {
    group: String,
    attrs: BTreeMap<String, String>,
}

fn type_size(ty: &str) -> Option<usize> {
    Some(match ty {
        "Int8" | "UInt8" => 1,
        "Int16" | "UInt16" => 2,
        "Int32" | "UInt32" | "Float32" => 4,
        "Int64" | "UInt64" | "Float64" => 8,
        _ => return None,
    })
}

fn decode(ty: &str, bytes: &[u8], big: bool) -> Result<Vec<f64>, LoadError> {
    let Some(sz) = type_size(ty) else {
        return bad(format!("unsupported array type {ty}"));
    };
    let mut out = Vec::with_capacity(bytes.len() / sz);
    for c in bytes.chunks_exact(sz) {
        let mut a = [0u8; 8];
        let v = if big {
            a[8 - sz..].copy_from_slice(c);
            u64::from_be_bytes(a)
        } else {
            a[..sz].copy_from_slice(c);
            u64::from_le_bytes(a)
        };
        let signed = |bits: u32| -> i64 { ((v << (64 - bits)) as i64) >> (64 - bits) };
        out.push(match ty {
            "Float32" => f64::from(f32::from_bits(v as u32)),
            "Float64" => f64::from_bits(v),
            "Int8" => signed(8) as f64,
            "Int16" => signed(16) as f64,
            "Int32" => signed(32) as f64,
            "Int64" => v as i64 as f64,
            _ => v as f64,
        });
    }
    Ok(out)
}

pub fn load_vtu(bytes: &[u8]) -> Result<Vtu, LoadError> {
    let marker = b"<AppendedData";
    let split = bytes.windows(marker.len()).position(|w| w == marker);
    let (head, payload): (&[u8], Option<&[u8]>) = match split {
        Some(i) => {
            let rest = &bytes[i..];
            let gt = rest.iter().position(|b| *b == b'>');
            let us = gt.and_then(|g| rest[g..].iter().position(|b| *b == b'_').map(|u| g + u));
            match us {
                Some(u) => (&bytes[..i], Some(&rest[u + 1..])),
                None => return bad("AppendedData has no `_` marker"),
            }
        }
        None => (bytes, None),
    };
    let mut head = String::from_utf8_lossy(head).into_owned();
    if split.is_some() {
        head.push_str("</VTKFile>");
    }

    let mut reader = Reader::from_str(&head);
    let mut big = false;
    let mut header_u64 = false;
    let mut vtu = Vtu::default();
    let mut group = String::new();
    let mut open: Option<Decl> = None;
    let mut ascii_text = String::new();
    let mut decls: Vec<(Decl, Option<String>)> = Vec::new();
    loop {
        let ev = reader
            .read_event()
            .map_err(|e| LoadError::Format(format!("XML: {e}")))?;
        let attrs_of = |e: &quick_xml::events::BytesStart<'_>| -> BTreeMap<String, String> {
            e.attributes()
                .flatten()
                .map(|a| {
                    (
                        String::from_utf8_lossy(a.key.as_ref()).into_owned(),
                        a.unescape_value()
                            .map(|v| v.into_owned())
                            .unwrap_or_default(),
                    )
                })
                .collect()
        };
        match ev {
            Event::Start(e) | Event::Empty(e) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).into_owned();
                let attrs = attrs_of(&e);
                match name.as_str() {
                    "VTKFile" => {
                        big = attrs.get("byte_order").is_some_and(|b| b == "BigEndian");
                        header_u64 = attrs.get("header_type").is_some_and(|b| b == "UInt64");
                        if attrs.get("type").map(String::as_str) != Some("UnstructuredGrid") {
                            return bad("not an UnstructuredGrid file");
                        }
                    }
                    "Piece" => {
                        let num = |k: &str| -> Result<usize, LoadError> {
                            attrs
                                .get(k)
                                .and_then(|v| v.parse().ok())
                                .ok_or_else(|| LoadError::Format(format!("Piece lacks {k}")))
                        };
                        vtu.points = num("NumberOfPoints")?;
                        vtu.cells = num("NumberOfCells")?;
                    }
                    "PointData" | "CellData" | "Points" | "Cells" => group = name,
                    "DataArray" => {
                        let d = Decl {
                            group: group.clone(),
                            attrs,
                        };
                        if matches!(ev_kind(&d), Kind::Inline) {
                            open = Some(d);
                            ascii_text.clear();
                        } else {
                            decls.push((d, None));
                        }
                    }
                    _ => {}
                }
            }
            Event::Text(t) => {
                if open.is_some() {
                    ascii_text.push_str(&t.decode().map_err(|e| LoadError::Format(e.to_string()))?);
                }
            }
            Event::End(e) => {
                if e.name().as_ref() == b"DataArray" {
                    if let Some(d) = open.take() {
                        decls.push((d, Some(std::mem::take(&mut ascii_text))));
                    }
                }
            }
            Event::Eof => break,
            _ => {}
        }
    }

    for (d, text) in decls {
        let ty = d.attrs.get("type").cloned().unwrap_or_default();
        let comps: usize = d
            .attrs
            .get("NumberOfComponents")
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);
        let name = d.attrs.get("Name").cloned().unwrap_or_default();
        let values = match d.attrs.get("format").map(String::as_str) {
            Some("ascii") => text
                .unwrap_or_default()
                .split_whitespace()
                .map(|t| {
                    t.parse::<f64>()
                        .map_err(|_| LoadError::Format(format!("`{t}` in {name}")))
                })
                .collect::<Result<Vec<_>, _>>()?,
            Some("appended") => {
                let Some(p) = payload else {
                    return bad(format!("array {name} is appended but there is no payload"));
                };
                let off: usize = d
                    .attrs
                    .get("offset")
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| LoadError::Format(format!("{name}: no offset")))?;
                let hs = if header_u64 { 8 } else { 4 };
                let Some(h) = p.get(off..off + hs) else {
                    return bad(format!("{name}: offset {off} outside the payload"));
                };
                let mut a = [0u8; 8];
                let n = if big {
                    a[8 - hs..].copy_from_slice(h);
                    u64::from_be_bytes(a)
                } else {
                    a[..hs].copy_from_slice(h);
                    u64::from_le_bytes(a)
                } as usize;
                let Some(data) = p.get(off + hs..off + hs + n) else {
                    return bad(format!("{name}: {n} bytes at {off} run past the payload"));
                };
                let sz = type_size(&ty).ok_or_else(|| LoadError::Format(format!("type {ty}")))?;
                if !n.is_multiple_of(sz) {
                    return bad(format!("{name}: {n} bytes is not a multiple of {sz}"));
                }
                decode(&ty, data, big)?
            }
            other => return bad(format!("{name}: unsupported format {other:?}")),
        };
        if comps == 0 || values.len() % comps != 0 {
            return bad(format!(
                "{name}: {} values for {comps} components",
                values.len()
            ));
        }
        let arr = Array {
            name: name.clone(),
            components: comps,
            ty,
            values,
        };
        match (d.group.as_str(), name.as_str()) {
            ("PointData", _) => vtu.point_data.push(arr),
            ("CellData", _) => vtu.cell_data.push(arr),
            ("Points", _) => vtu.coordinates = arr.values,
            ("Cells", "connectivity") => {
                vtu.connectivity = arr.values.iter().map(|v| *v as i64).collect()
            }
            ("Cells", "offsets") => vtu.offsets = arr.values.iter().map(|v| *v as i64).collect(),
            ("Cells", "types") => vtu.types = arr.values.iter().map(|v| *v as i64).collect(),
            _ => {}
        }
    }
    validate(&vtu)?;
    Ok(vtu)
}

enum Kind {
    Inline,
    Other,
}

fn ev_kind(d: &Decl) -> Kind {
    match d.attrs.get("format").map(String::as_str) {
        Some("ascii") | Some("binary") => Kind::Inline,
        _ => Kind::Other,
    }
}

fn validate(v: &Vtu) -> Result<(), LoadError> {
    if v.coordinates.len() != 3 * v.points {
        return bad(format!(
            "{} coordinates for {} points",
            v.coordinates.len(),
            v.points
        ));
    }
    if v.offsets.len() != v.cells || v.types.len() != v.cells {
        return bad("cell offsets/types do not match NumberOfCells");
    }
    if v.offsets.windows(2).any(|w| w[1] < w[0])
        || v.offsets.last().copied().unwrap_or(0) as usize != v.connectivity.len()
    {
        return bad("cell offsets are not a partition of the connectivity");
    }
    if v.connectivity
        .iter()
        .any(|c| *c < 0 || *c as usize >= v.points)
    {
        return bad("connectivity refers to a missing point");
    }
    for a in &v.point_data {
        if a.tuples() != v.points {
            return bad(format!(
                "point array {} has {} tuples for {} points",
                a.name,
                a.tuples(),
                v.points
            ));
        }
    }
    for a in &v.cell_data {
        if a.tuples() != v.cells {
            return bad(format!(
                "cell array {} has {} tuples for {} cells",
                a.name,
                a.tuples(),
                v.cells
            ));
        }
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub struct Ep {
    pub nodes: usize,
    pub elements: usize,
    pub components: usize,
    pub timesteps: usize,
    /// Field declarations after the counts: `("scalar", "Temperature")`.
    pub fields: Vec<(String, String)>,
    /// `#time step substep time` records.
    pub times: Vec<f64>,
}

pub fn load_ep(text: &str) -> Result<Ep, LoadError> {
    let mut lines = text.lines();
    let first = lines
        .next()
        .ok_or_else(|| LoadError::Format("empty .ep".into()))?;
    let mut tok = first.split_whitespace();
    let mut num = |what: &str| -> Result<usize, LoadError> {
        tok.next()
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| LoadError::Format(format!(".ep header lacks {what}")))
    };
    let nodes = num("node count")?;
    let elements = num("element count")?;
    let components = num("component count")?;
    let timesteps = num("timestep count")?;
    let rest: Vec<&str> = tok.collect();
    let mut fields = Vec::new();
    let mut i = 0;
    while i + 1 < rest.len() {
        let kind = rest[i].trim_end_matches(':');
        if !matches!(kind, "scalar" | "vector" | "tensor") {
            return bad(format!("unknown field kind `{}`", rest[i]));
        }
        fields.push((kind.to_string(), rest[i + 1].to_string()));
        i += 2;
    }
    let width: usize = fields
        .iter()
        .map(|(k, _): &(String, String)| if k == "scalar" { 1 } else { 3 })
        .sum();
    if width != components && !fields.iter().any(|(k, _)| k == "tensor") {
        return bad(format!(
            "header says {components} components, fields add up to {width}"
        ));
    }
    let times: Vec<f64> = text
        .lines()
        .filter(|l| l.starts_with("#time"))
        .filter_map(|l| l.split_whitespace().nth(3).and_then(|t| t.parse().ok()))
        .collect();
    if times.len() != timesteps {
        return bad(format!(
            "header says {timesteps} timesteps, found {} #time records",
            times.len()
        ));
    }
    Ok(Ep {
        nodes,
        elements,
        timesteps,
        components,
        fields,
        times,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny VTU in the exact layout Elmer writes (appended raw, `Int32`
    /// connectivity, `Float64` data), built by hand.
    fn synthetic(corrupt_offset: bool) -> Vec<u8> {
        let pts: [f64; 9] = [0., 0., 0., 1., 0., 0., 0., 1., 0.];
        let temp: [f64; 3] = [1., 2., 3.];
        let mut payload = Vec::new();
        let mut offs = Vec::new();
        let mut push = |bytes: Vec<u8>| {
            offs.push(payload.len());
            payload.extend((bytes.len() as u32).to_le_bytes());
            payload.extend(bytes);
        };
        push(temp.iter().flat_map(|v| v.to_le_bytes()).collect());
        push(pts.iter().flat_map(|v| v.to_le_bytes()).collect());
        push([0i32, 1, 2].iter().flat_map(|v| v.to_le_bytes()).collect());
        push([3i32].iter().flat_map(|v| v.to_le_bytes()).collect());
        push([5i32].iter().flat_map(|v| v.to_le_bytes()).collect());
        if corrupt_offset {
            offs[0] = 10_000;
        }
        let head = format!(
            "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
             <UnstructuredGrid><Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">\n\
             <PointData><DataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\" format=\"appended\" offset=\"{}\"/></PointData>\n\
             <CellData></CellData>\n\
             <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"{}\"/></Points>\n\
             <Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"{}\"/>\n\
             <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"{}\"/>\n\
             <DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"{}\"/></Cells>\n\
             </Piece></UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_",
            offs[0], offs[1], offs[2], offs[3], offs[4]
        );
        let mut out = head.into_bytes();
        out.extend(payload);
        out.extend(b"\n</AppendedData>\n</VTKFile>\n");
        out
    }

    #[test]
    fn reads_the_elmer_appended_layout() {
        let v = load_vtu(&synthetic(false)).unwrap();
        assert_eq!((v.points, v.cells), (3, 1));
        let t = v.point_array("Temperature").unwrap();
        assert_eq!(t.components, 1);
        assert_eq!(t.values, [1., 2., 3.]);
        assert_eq!(v.types, [5]);
    }

    #[test]
    fn rejects_offsets_outside_the_payload() {
        assert!(load_vtu(&synthetic(true)).is_err());
    }

    #[test]
    fn reads_ascii_arrays() {
        let xml = "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\"><UnstructuredGrid>\
            <Piece NumberOfPoints=\"2\" NumberOfCells=\"1\">\
            <PointData><DataArray type=\"Float32\" Name=\"u\" NumberOfComponents=\"1\" format=\"ascii\">0.5 1.5</DataArray></PointData>\
            <CellData></CellData>\
            <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">0 0 0 1 0 0</DataArray></Points>\
            <Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0 1</DataArray>\
            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">2</DataArray>\
            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">3</DataArray></Cells>\
            </Piece></UnstructuredGrid></VTKFile>";
        let v = load_vtu(xml.as_bytes()).unwrap();
        assert_eq!(v.point_array("u").unwrap().values, [0.5, 1.5]);
    }

    #[test]
    fn ep_header_and_times() {
        let ep = load_ep(
            "         9        16      1      1 scalar: Temperature\n#File started at: now\n#time       1      1  1.00000000E+000\n",
        )
        .unwrap();
        assert_eq!(
            (ep.nodes, ep.elements, ep.components, ep.timesteps),
            (9, 16, 1, 1)
        );
        assert_eq!(ep.fields, [("scalar".into(), "Temperature".into())]);
        assert_eq!(ep.times, [1.0]);
        assert!(load_ep("1 2 3 4 scalar: T\n").is_err(), "missing #time");
    }
}
