//! Elmer mesh database: in-memory model and file reader.
//!
//! Replicates `ElmerMeshReader` (`MeshIO.F90`) steps 1-4 for serial
//! meshes: header, nodes (ASCII, `.bin`, `.sbin`), bulk elements,
//! boundary elements with parent links, node/element permutation by
//! file tags, and mesh dimension. Parallel (`partitioning.N`) input is
//! an explicit error (later lane).

use crate::elements::element_info;
use std::fmt;

/// A loaded serial mesh. Node/element indices are 0-based; file tags
/// (1-based) are kept for round-trips and global IDs.
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    pub dim: usize,
    pub max_dim: usize,
    pub nodes: Vec<[f64; 3]>,
    /// File tag per node (`mesh.nodes` column 1), i.e. serial global IDs.
    pub node_tags: Vec<i32>,
    pub bulk: Vec<Elem>,
    pub boundary: Vec<BoundElem>,
    pub edges: Option<EdgeTable>,
}

/// Bulk element. `parents` live on boundary elements only.
#[derive(Debug, Clone)]
pub struct Elem {
    pub code: u16,
    pub body: i32,
    /// 0-based node indices, Elmer local order.
    pub nodes: Vec<u32>,
    /// File tag (1-based position in `mesh.elements`).
    pub tag: i32,
    /// Global edge ids, Elmer canonical per-type order. Empty until
    /// [`Mesh::find_edges`] runs.
    pub edges: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct BoundElem {
    pub code: u16,
    /// BC tag from `mesh.boundary` (remapped by transforms later).
    pub constraint: i32,
    /// Parent bulk-element indices: [Left, Right].
    pub parents: [Option<usize>; 2],
    pub nodes: Vec<u32>,
    pub tag: i32,
    /// Global edge ids, like bulk elements. Empty until
    /// [`Mesh::find_edges`] runs.
    pub edges: Vec<u32>,
}

/// Global edge table. `edges[e] = (n1, n2)` with 0-based node ids.
#[derive(Debug, Clone, Default)]
pub struct EdgeTable {
    pub edges: Vec<(u32, u32)>,
}

#[derive(Debug, Clone)]
pub struct MeshError {
    pub msg: String,
}

impl MeshError {
    fn new(msg: impl Into<String>) -> Self {
        MeshError { msg: msg.into() }
    }
}

impl fmt::Display for MeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "mesh: {}", self.msg)
    }
}

impl std::error::Error for MeshError {}

/// Load a serial Elmer mesh database from `dir`.
pub fn load(dir: &str) -> Result<Mesh, MeshError> {
    let mut mesh = Mesh::default();
    let header = read_header(dir)?;
    let (n_nodes, n_bulk, n_bound) = header;
    mesh.nodes = read_nodes(dir, n_nodes)?;
    let (bulk_tags, bulk) = read_elements(dir, n_bulk)?;
    let (bound_tags, boundary) = read_boundary(dir, n_bulk, n_bound, &bulk_tags)?;
    mesh.bulk = bulk;
    mesh.boundary = boundary;
    let _ = (bulk_tags, bound_tags);
    set_mesh_dimension(&mut mesh);
    Ok(mesh)
}

fn read_header(dir: &str) -> Result<(usize, usize, usize), MeshError> {
    // New (ElmerGrid) format: counts, type-count, then (code, count) rows.
    let text = std::fs::read_to_string(format!("{dir}/mesh.header"))
        .map_err(|e| MeshError::new(format!("cannot open mesh.header: {e}")))?;
    let mut it = text.split_whitespace();
    let num = |it: &mut std::str::SplitWhitespace<'_>, what: &str| -> Result<usize, MeshError> {
        it.next()
            .and_then(|w| w.parse().ok())
            .ok_or_else(|| MeshError::new(format!("bad mesh.header ({what})")))
    };
    let n_nodes = num(&mut it, "nodes")?;
    let n_bulk = num(&mut it, "bulk")?;
    let n_bound = num(&mut it, "boundary")?;
    let n_types = num(&mut it, "type count")?;
    for _ in 0..n_types {
        let code: u32 = num(&mut it, "type code")? as u32;
        let _count = num(&mut it, "type count")?;
        if element_info(code as u16).is_none() {
            // Unknown codes fail at element load (GetElementType warns,
            // the readers fatal). Recorded here for a sharper message.
            return Err(MeshError::new(format!("unknown element type code {code}")));
        }
    }
    Ok((n_nodes, n_bulk, n_bound))
}

/// Parse one line of blank-separated integers, Fortran `read_ints`
/// semantics: `/` marks halo data (partition prefix), `-` starts
/// negative numbers, parse stops at other characters.
fn read_ints(s: &str) -> (Vec<i64>, bool) {
    let b = s.as_bytes();
    let k = s.trim_end().len();
    let mut l = 0;
    let mut out = Vec::new();
    let mut halo = false;
    while l < k && out.len() < 64 {
        while l < k {
            let c = b[l];
            if c == b' ' || c == b'\t' {
                l += 1;
                continue;
            } else if c == b'/' {
                halo = true;
            } else {
                break;
            }
            l += 1;
        }
        if l >= k {
            break;
        }
        if !(b[l] == b'-' || b[l].is_ascii_digit()) {
            break;
        }
        let mut m = l + 1;
        while m < k && b[m].is_ascii_digit() {
            m += 1;
        }
        if let Ok(v) = s[l..m].parse::<i64>() {
            out.push(v);
        }
        l = m;
    }
    (out, halo)
}

fn read_nodes(dir: &str, n: usize) -> Result<Vec<[f64; 3]>, MeshError> {
    let base = format!("{dir}/mesh.nodes");
    // Binary variants first (Elmer tries ASCII, then .bin, then .sbin).
    if let Ok(bytes) = std::fs::read(format!("{base}.bin")) {
        return read_nodes_bin(&bytes, false, n, &base);
    }
    if let Ok(bytes) = std::fs::read(format!("{base}.sbin")) {
        return read_nodes_bin(&bytes, true, n, &base);
    }
    let text = std::fs::read_to_string(&base)
        .map_err(|e| MeshError::new(format!("cannot open mesh.nodes: {e}")))?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let mut nodes = Vec::with_capacity(n);
    let mut tags: Vec<i32> = Vec::with_capacity(n);
    for j in 0..n {
        let line = lines
            .next()
            .ok_or_else(|| MeshError::new(format!("short mesh.nodes at {j}")))?;
        let mut w = line.split_whitespace();
        let tag: i32 = w
            .next()
            .and_then(|x| x.parse().ok())
            .ok_or_else(|| MeshError::new(format!("bad node tag line {j}")))?;
        let _dash = w.next();
        let mut c = [0.0; 3];
        for v in c.iter_mut() {
            *v = w
                .next()
                .and_then(|x| x.parse().ok())
                .ok_or_else(|| MeshError::new(format!("bad node coords line {j}")))?;
        }
        tags.push(tag);
        nodes.push(c);
    }
    // Node permutation by tags (PermuteNodeNumbering): identity when
    // tags are 1..n in order.
    let permuted = tags.iter().enumerate().any(|(i, &t)| t != i as i32 + 1);
    if permuted {
        // Order nodes by tag: node with tag t goes to position t-1 when
        // tags are a permutation of 1..n (Elmer's LocalPerm path).
        let mut by_tag = vec![[0.0; 3]; n];
        let mut seen = vec![false; n];
        for (pos, c) in nodes.iter().enumerate() {
            let t = tags[pos];
            if t >= 1 && (t as usize) <= n {
                by_tag[t as usize - 1] = *c;
                seen[t as usize - 1] = true;
            }
        }
        if seen.iter().all(|&s| s) {
            return Ok(by_tag);
        }
        // Non-permutation tags (parallel leftovers): keep file order.
    }
    Ok(nodes)
}

fn read_nodes_bin(
    bytes: &[u8],
    single: bool,
    n: usize,
    base: &str,
) -> Result<Vec<[f64; 3]>, MeshError> {
    // Fortran stream: (tag:i32, xyz:3xf64|x3f32) per node, native endian.
    let rec = if single { 4 + 3 * 4 } else { 4 + 3 * 8 };
    if bytes.len() < n * rec {
        return Err(MeshError::new(format!("short {base}.bin")));
    }
    let mut nodes = Vec::with_capacity(n);
    let mut tags = Vec::with_capacity(n);
    for j in 0..n {
        let o = j * rec;
        let tag = i32::from_ne_bytes(bytes[o..o + 4].try_into().unwrap());
        let mut c = [0.0; 3];
        for k in 0..3 {
            c[k] = if single {
                f32::from_ne_bytes(bytes[o + 4 + k * 4..o + 8 + k * 4].try_into().unwrap()) as f64
            } else {
                f64::from_ne_bytes(bytes[o + 4 + k * 8..o + 12 + k * 8].try_into().unwrap())
            };
        }
        tags.push(tag);
        nodes.push(c);
    }
    let permuted = tags.iter().enumerate().any(|(i, &t)| t != i as i32 + 1);
    if permuted {
        let mut by_tag = vec![[0.0; 3]; n];
        let mut seen = vec![false; n];
        for (pos, c) in nodes.iter().enumerate() {
            let t = tags[pos];
            if t >= 1 && (t as usize) <= n {
                by_tag[t as usize - 1] = *c;
                seen[t as usize - 1] = true;
            }
        }
        if seen.iter().all(|&s| s) {
            return Ok(by_tag);
        }
    }
    Ok(nodes)
}

fn read_elements(dir: &str, n_bulk: usize) -> Result<(Vec<i32>, Vec<Elem>), MeshError> {
    let base = format!("{dir}/mesh.elements");
    let mut tags = Vec::with_capacity(n_bulk);
    let mut bulk = Vec::with_capacity(n_bulk);
    match std::fs::read_to_string(&base) {
        Ok(text) => {
            let mut lines = text.lines().filter(|l| !l.trim().is_empty());
            for j in 0..n_bulk {
                let line = lines.next().ok_or_else(|| {
                    MeshError::new(format!("short mesh.elements at {j}"))
                })?;
                let (vals, halo) = read_ints(line);
                let (off, _partn) = if halo { (1, vals.get(1)) } else { (0, None) };
                if vals.len() < 3 + off {
                    return Err(MeshError::new(format!(
                        "short mesh.elements line {j}"
                    )));
                }
                let tag = vals[0] as i32;
                let body = vals[off + 1] as i32;
                let code = vals[off + 2] as u16;
                let info = element_info(code).ok_or_else(|| {
                    MeshError::new(format!("element type {code} could not be associated"))
                })?;
                let nn = info.n_nodes;
                if vals.len() < 3 + off + nn {
                    return Err(MeshError::new(format!(
                        "mesh.elements line {j} does not contain enough entries"
                    )));
                }
                // File nodes are 1-based tags; with a permutation they are
                // remapped below. Store 0-based positions for now.
                let nodes: Vec<u32> = vals[3 + off..3 + off + nn]
                    .iter()
                    .map(|&v| (v - 1) as u32)
                    .collect();
                tags.push(tag);
                bulk.push(Elem { code, body, nodes, tag, edges: Vec::new() });
            }
        }
        Err(_) => {
            let bytes = std::fs::read(format!("{base}.bin")).map_err(|_| {
                MeshError::new("cannot open mesh.elements".to_string())
            })?;
            // (tag, partn, body, code):i32 + nn:i32 node tags.
            let mut o = 0;
            for j in 0..n_bulk {
                if bytes.len() < o + 16 {
                    return Err(MeshError::new(format!("short mesh.elements.bin at {j}")));
                }
                let tag = i32::from_ne_bytes(bytes[o..o + 4].try_into().unwrap());
                let _partn = i32::from_ne_bytes(bytes[o + 4..o + 8].try_into().unwrap());
                let body = i32::from_ne_bytes(bytes[o + 8..o + 12].try_into().unwrap());
                let code =
                    i32::from_ne_bytes(bytes[o + 12..o + 16].try_into().unwrap()) as u16;
                o += 16;
                let info = element_info(code).ok_or_else(|| {
                    MeshError::new(format!("element type {code} could not be associated"))
                })?;
                let nn = info.n_nodes;
                if bytes.len() < o + 4 * nn {
                    return Err(MeshError::new(format!("short mesh.elements.bin nodes at {j}")));
                }
                let mut nodes = Vec::with_capacity(nn);
                for _ in 0..nn {
                    let v = i32::from_ne_bytes(bytes[o..o + 4].try_into().unwrap());
                    o += 4;
                    nodes.push((v - 1) as u32);
                }
                tags.push(tag);
                bulk.push(Elem { code, body, nodes, tag, edges: Vec::new() });
            }
        }
    }
    // Element permutation is handled by the boundary parent remap
    // (LocalEPerm); bulk nodes need the node permutation only.
    Ok((tags, bulk))
}

fn read_boundary(
    dir: &str,
    n_bulk: usize,
    n_bound: usize,
    bulk_tags: &[i32],
) -> Result<(Vec<i32>, Vec<BoundElem>), MeshError> {
    let base = format!("{dir}/mesh.boundary");
    let mut tags = Vec::with_capacity(n_bound);
    let mut boundary = Vec::with_capacity(n_bound);
    // LocalEPerm for parent remap when bulk tags permute.
    let permuted = bulk_tags.iter().enumerate().any(|(i, &t)| t != i as i32 + 1);
    let local_eperm: Option<Vec<usize>> = if permuted {
        let min = *bulk_tags.iter().min().unwrap_or(&1);
        let max = *bulk_tags.iter().max().unwrap_or(&0);
        let mut perm = vec![0usize; (max - min + 1) as usize];
        for (i, &t) in bulk_tags.iter().enumerate() {
            if t >= min && t <= max {
                perm[(t - min) as usize] = i;
            }
        }
        Some(perm)
    } else {
        None
    };
    let min_tag = bulk_tags.iter().min().copied().unwrap_or(1);
    let max_tag = bulk_tags.iter().max().copied().unwrap_or(0);

    let text = std::fs::read_to_string(&base).or_else(|_| {
        std::fs::read_to_string(format!("{base}.bin"))
            .map(|_| String::new())
            .map_err(|_| MeshError::new("cannot open mesh.boundary".to_string()))
    })?;
    // Binary boundary files use fixed records; ASCII parsed per line.
    let is_bin = std::path::Path::new(&format!("{base}.bin")).is_file()
        && std::fs::read_to_string(&base).is_err();
    if is_bin {
        return read_boundary_bin(dir, n_bulk, n_bound, bulk_tags, &local_eperm, min_tag, max_tag);
    }
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    for j in 0..n_bound {
        let line = lines.next().ok_or_else(|| {
            MeshError::new(format!("short mesh.boundary at {j}"))
        })?;
        let (vals, halo) = read_ints(line);
        let off = if halo { 1 } else { 0 };
        if vals.len() < 5 + off {
            return Err(MeshError::new(format!("short mesh.boundary line {j}")));
        }
        let tag = vals[0] as i32;
        let bndry = vals[off + 1] as i32;
        let mut left = vals[off + 2];
        let mut right = vals[off + 3];
        let code = vals[off + 4] as u16;
        let info = element_info(code).ok_or_else(|| {
            MeshError::new(format!("element type {code} could not be associated"))
        })?;
        let nn = info.n_nodes;
        if vals.len() < 5 + off + nn {
            return Err(MeshError::new(format!(
                "mesh.boundary line {j} does not contain enough entries"
            )));
        }
        // Parent remap through LocalEPerm; out-of-range parents are a
        // non-fatal error reset to zero (ReadBoundaryFile).
        let mut bad = false;
        if left >= min_tag as i64 && left <= max_tag as i64 {
            if let Some(perm) = local_eperm.as_ref() {
                left = perm[(left - min_tag as i64) as usize] as i64 + 1;
                // stored below as 0-based index
            }
        } else if left > 0 {
            bad = true;
            left = 0;
        }
        if right >= min_tag as i64 && right <= max_tag as i64 {
            if let Some(perm) = local_eperm.as_ref() {
                right = perm[(right - min_tag as i64) as usize] as i64 + 1;
            }
        } else if right > 0 {
            bad = true;
            right = 0;
        }
        let _ = bad;
        // Elmer always wants "Left" present: Right-only becomes Left.
        if right >= 1 && left == 0 {
            left = right;
            right = 0;
        }
        let parents = [
            if left >= 1 { Some((left - 1) as usize) } else { None },
            if right >= 1 { Some((right - 1) as usize) } else { None },
        ];
        let nodes: Vec<u32> = vals[5 + off..5 + off + nn]
            .iter()
            .map(|&v| (v - 1) as u32)
            .collect();
        tags.push(tag);
        boundary.push(BoundElem { code, constraint: bndry, parents, nodes, tag, edges: Vec::new() });
    }
    let _ = n_bulk;
    Ok((tags, boundary))
}

fn read_boundary_bin(
    dir: &str,
    _n_bulk: usize,
    n_bound: usize,
    bulk_tags: &[i32],
    local_eperm: &Option<Vec<usize>>,
    min_tag: i32,
    max_tag: i32,
) -> Result<(Vec<i32>, Vec<BoundElem>), MeshError> {
    let bytes = std::fs::read(format!("{dir}/mesh.boundary.bin"))
        .map_err(|e| MeshError::new(format!("cannot open mesh.boundary.bin: {e}")))?;
    // (tag, partn, bndry, left, right, code):i32 + nn:i32 node tags.
    let mut o = 0;
    let mut tags = Vec::with_capacity(n_bound);
    let mut boundary = Vec::with_capacity(n_bound);
    for j in 0..n_bound {
        if bytes.len() < o + 24 {
            return Err(MeshError::new(format!("short mesh.boundary.bin at {j}")));
        }
        let rd = |o: usize| i32::from_ne_bytes(bytes[o..o + 4].try_into().unwrap());
        let (tag, _partn, bndry, mut left, mut right, code) =
            (rd(o), rd(o + 4), rd(o + 8), rd(o + 12) as i64, rd(o + 16) as i64, rd(o + 20) as u16);
        o += 24;
        let info = element_info(code).ok_or_else(|| {
            MeshError::new(format!("element type {code} could not be associated"))
        })?;
        let nn = info.n_nodes;
        if bytes.len() < o + 4 * nn {
            return Err(MeshError::new(format!("short mesh.boundary.bin nodes at {j}")));
        }
        let mut nodes = Vec::with_capacity(nn);
        for _ in 0..nn {
            nodes.push((rd(o) - 1) as u32);
            o += 4;
        }
        if left >= min_tag as i64 && left <= max_tag as i64 {
            if let Some(perm) = local_eperm.as_ref() {
                left = perm[(left - min_tag as i64) as usize] as i64 + 1;
            }
        } else if left > 0 {
            left = 0;
        }
        if right >= min_tag as i64 && right <= max_tag as i64 {
            if let Some(perm) = local_eperm.as_ref() {
                right = perm[(right - min_tag as i64) as usize] as i64 + 1;
            }
        } else if right > 0 {
            right = 0;
        }
        if right >= 1 && left == 0 {
            left = right;
            right = 0;
        }
        let _ = bulk_tags;
        tags.push(tag);
        boundary.push(BoundElem {
            code,
            constraint: bndry,
            parents: [
                if left >= 1 { Some((left - 1) as usize) } else { None },
                if right >= 1 { Some((right - 1) as usize) } else { None },
            ],
            nodes,
            tag,
            edges: Vec::new(),
        });
    }
    Ok((tags, boundary))
}

pub(crate) fn set_mesh_dimension(mesh: &mut Mesh) {
    if mesh.nodes.is_empty() {
        return;
    }
    let (x0, y0, z0) = (mesh.nodes[0][0], mesh.nodes[0][1], mesh.nodes[0][2]);
    let mut c = [false; 3];
    for n in &mesh.nodes {
        if n[0] != x0 {
            c[0] = true;
        }
        if n[1] != y0 {
            c[1] = true;
        }
        if n[2] != z0 {
            c[2] = true;
        }
    }
    mesh.dim = c.iter().filter(|&&b| b).count();
    mesh.max_dim = c.iter().rposition(|&b| b).map(|i| i + 1).unwrap_or(0);
}

// ---- Writer (WriteMeshToDisk2 layout: type-count header) ----

use std::fmt::Write as _fmt_write;

/// Write a mesh database `dir` in the layout `ReadMeshStep` consumes:
/// `mesh.header` (counts + type table), `mesh.nodes`, `mesh.elements`,
/// `mesh.boundary`. Tags are sequential positions.
pub fn write_db(mesh: &Mesh, dir: &str) -> Result<(), MeshError> {
    std::fs::create_dir_all(dir).map_err(|e| MeshError::new(format!("mkdir: {e}")))?;
    let mut counts: std::collections::BTreeMap<u16, usize> = std::collections::BTreeMap::new();
    for e in &mesh.bulk {
        *counts.entry(e.code).or_insert(0) += 1;
    }
    for b in &mesh.boundary {
        *counts.entry(b.code).or_insert(0) += 1;
    }
    let mut header = String::new();
    let _ = writeln!(
        header,
        "{} {} {}",
        mesh.nodes.len(),
        mesh.bulk.len(),
        mesh.boundary.len()
    );
    let _ = writeln!(header, "{}", counts.len());
    for (code, n) in &counts {
        let _ = writeln!(header, "{code} {n}");
    }
    std::fs::write(format!("{dir}/mesh.header"), header)
        .map_err(|e| MeshError::new(format!("write header: {e}")))?;

    let mut nodes = String::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        let _ = writeln!(nodes, "{} -1 {:.17e} {:.17e} {:.17e}", i + 1, n[0], n[1], n[2]);
    }
    std::fs::write(format!("{dir}/mesh.nodes"), nodes)
        .map_err(|e| MeshError::new(format!("write nodes: {e}")))?;

    let mut elems = String::new();
    for (i, e) in mesh.bulk.iter().enumerate() {
        let _ = write!(elems, "{} {} {} ", i + 1, e.body, e.code);
        for n in &e.nodes {
            let _ = write!(elems, "{} ", n + 1);
        }
        let _ = writeln!(elems);
    }
    std::fs::write(format!("{dir}/mesh.elements"), elems)
        .map_err(|e| MeshError::new(format!("write elements: {e}")))?;

    // Boundary parents that point outside the bulk range (e.g. the 0D
    // intersection elements created in memory, whose parents are other
    // boundary elements) cannot round-trip through Elmer's file format:
    // its reader rejects them as out of range. Resolve one hop toward
    // the bulk through the parent boundary element's own Left parent
    // (geometrically the shared bulk element), else drop the link.
    // In-memory parents are untouched; only the file encoding changes.
    let resolve = |mut p: Option<usize>| -> i64 {
        for _ in 0..2 {
            let q = match p {
                Some(q) => q,
                None => return 0,
            };
            if q < mesh.bulk.len() {
                return (q + 1) as i64;
            }
            let bi = q - mesh.bulk.len();
            if bi >= mesh.boundary.len() {
                return 0;
            }
            let be = &mesh.boundary[bi];
            p = be.parents[0].or(be.parents[1]);
        }
        0
    };
    let mut bound = String::new();
    for (i, b) in mesh.boundary.iter().enumerate() {
        let left = resolve(b.parents[0]);
        let right = resolve(b.parents[1]);
        let _ = write!(bound, "{} {} {left} {right} {} ", i + 1, b.constraint, b.code);
        for n in &b.nodes {
            let _ = write!(bound, "{} ", n + 1);
        }
        let _ = writeln!(bound);
    }
    std::fs::write(format!("{dir}/mesh.boundary"), bound)
        .map_err(|e| MeshError::new(format!("write boundary: {e}")))?;
    Ok(())
}
