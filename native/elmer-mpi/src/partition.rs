//! Reader for ElmerGrid partition directories (`partitioning.N/part.K.*`).
//!
//! File layout (`ElmerGrid ... -metis N [-halo]`, parallel ElmerSolver format;
//! partitions are 1-based on disk, 0-based here):
//!
//! * `part.K.header`: `nodes elements boundary`, the element-type count, one
//!   `type count` line per type, then `shared_nodes extra`.
//! * `part.K.nodes`: `global_id -1 x y z`.
//! * `part.K.elements`: `id body type global_node_ids...`; a halo (ghost)
//!   element is written `id/owner_partition`.
//! * `part.K.boundary`: `id tag parent1 parent2 type global_node_ids...` with
//!   parents given as global bulk element ids (0 for none).
//! * `part.K.shared`: `global_id n_sharing p1 p2 ...`; `p1` owns the node.
//!
//! A node absent from `part.K.shared` is owned by this partition. Elements are
//! owned by this partition unless written as a halo element.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use elmer_core::{ElemId, ElementRef, ElementType, MeshView};

#[derive(Debug, thiserror::Error)]
pub enum PartError {
    #[error("{path}: {source}")]
    Io { path: PathBuf, source: std::io::Error },
    #[error("{path}:{line}: {msg}")]
    Parse { path: PathBuf, line: usize, msg: String },
    #[error("{path}: {msg}")]
    Inconsistent { path: PathBuf, msg: String },
}

/// `part.K.header` contents.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartHeader {
    pub n_nodes: usize,
    pub n_bulk: usize,
    pub n_boundary: usize,
    /// `(element type code, count)`.
    pub element_types: Vec<(u16, usize)>,
    pub n_shared_nodes: usize,
    /// Second number on the last header line; ElmerGrid writes 0.
    pub n_extra: usize,
}

/// An element with node indices local to the partition.
#[derive(Clone, Debug, PartialEq)]
pub struct PartElement {
    /// Global element id.
    pub gid: u64,
    /// Rank that owns the element; differs from this rank for ghost (halo) elements.
    pub owner: usize,
    /// Body id (bulk) or boundary tag (boundary).
    pub tag: i32,
    pub etype: ElementType,
    /// 0-based local node indices.
    pub nodes: Vec<u32>,
    /// Boundary elements only: global ids of the parent bulk elements.
    pub parent_gids: [Option<u64>; 2],
}

/// A node that more than one partition holds.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SharedNode {
    pub gid: u64,
    pub local: u32,
    /// Ranks holding the node (0-based, includes this rank); `parts[0]` owns it.
    pub parts: Vec<usize>,
}

/// One rank's piece of a partitioned mesh. Bulk elements come first, then
/// boundary elements, as in [`MeshView`].
#[derive(Clone, Debug)]
pub struct PartitionMesh {
    rank: usize,
    n_parts: usize,
    header: PartHeader,
    node_gid: Vec<u64>,
    coords: Vec<[f64; 3]>,
    node_owner: Vec<usize>,
    shared: Vec<SharedNode>,
    /// Bulk then boundary.
    elements: Vec<PartElement>,
    parents: Vec<[Option<ElemId>; 2]>,
    dim: usize,
}

/// `<mesh_dir>/partitioning.<n_parts>`.
pub fn partition_subdir(mesh_dir: &Path, n_parts: usize) -> PathBuf {
    mesh_dir.join(format!("partitioning.{n_parts}"))
}

impl PartitionMesh {
    /// Read rank `rank` of `n_parts` from `<mesh_dir>/partitioning.<n_parts>`.
    pub fn read(mesh_dir: &Path, rank: usize, n_parts: usize) -> Result<Self, PartError> {
        Self::read_dir(&partition_subdir(mesh_dir, n_parts), rank, n_parts)
    }

    /// Read rank `rank` from a `partitioning.N` directory itself.
    pub fn read_dir(dir: &Path, rank: usize, n_parts: usize) -> Result<Self, PartError> {
        let file = |ext: &str| dir.join(format!("part.{}.{ext}", rank + 1));
        let header = read_header(&file("header"))?;
        let shared_raw = read_shared(&file("shared"), rank, n_parts)?;
        let (node_gid, coords) = read_nodes(&file("nodes"))?;

        let local_of: HashMap<u64, u32> =
            node_gid.iter().enumerate().map(|(i, &g)| (g, i as u32)).collect();
        let inc = |msg: String| PartError::Inconsistent { path: dir.to_path_buf(), msg };
        if local_of.len() != node_gid.len() {
            return Err(inc("duplicate global node id".into()));
        }

        let mut node_owner = vec![rank; node_gid.len()];
        let mut shared = Vec::with_capacity(shared_raw.len());
        for (gid, parts) in shared_raw {
            let local = *local_of.get(&gid).ok_or_else(|| {
                inc(format!("shared node {gid} is not in part.{}.nodes", rank + 1))
            })?;
            node_owner[local as usize] = parts[0];
            shared.push(SharedNode { gid, local, parts });
        }

        let mut elements = read_elements(&file("elements"), false, rank, &local_of)?;
        let n_bulk = elements.len();
        elements.extend(read_elements(&file("boundary"), true, rank, &local_of)?);

        let count_err = |what: &str, got: usize, want: usize| {
            inc(format!("part.{}.header says {want} {what}, file has {got}", rank + 1))
        };
        if node_gid.len() != header.n_nodes {
            return Err(count_err("nodes", node_gid.len(), header.n_nodes));
        }
        if n_bulk != header.n_bulk {
            return Err(count_err("elements", n_bulk, header.n_bulk));
        }
        if elements.len() - n_bulk != header.n_boundary {
            return Err(count_err("boundary elements", elements.len() - n_bulk, header.n_boundary));
        }
        if shared.len() != header.n_shared_nodes {
            return Err(count_err("shared nodes", shared.len(), header.n_shared_nodes));
        }

        let by_gid: HashMap<u64, ElemId> =
            elements[..n_bulk].iter().enumerate().map(|(i, e)| (e.gid, ElemId(i as u32))).collect();
        let parents = elements
            .iter()
            .map(|e| e.parent_gids.map(|p| p.and_then(|g| by_gid.get(&g).copied())))
            .collect();

        // Highest-dimension bulk family: 2xx lines, 3xx/4xx faces, 5xx-8xx volumes.
        let dim = elements[..n_bulk]
            .iter()
            .map(|e| match e.etype.family() {
                1 | 2 => 1,
                3 | 4 => 2,
                _ => 3,
            })
            .max()
            .unwrap_or(0);

        Ok(Self {
            rank,
            n_parts,
            header,
            node_gid,
            coords,
            node_owner,
            shared,
            elements,
            parents,
            dim,
        })
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn n_parts(&self) -> usize {
        self.n_parts
    }

    pub fn header(&self) -> &PartHeader {
        &self.header
    }

    /// Global id of each local node.
    pub fn node_gids(&self) -> &[u64] {
        &self.node_gid
    }

    /// Owning rank of each local node.
    pub fn node_owners(&self) -> &[usize] {
        &self.node_owner
    }

    pub fn shared(&self) -> &[SharedNode] {
        &self.shared
    }

    /// `true` for nodes this rank owns; the complement are ghost copies.
    pub fn owned_nodes(&self) -> Vec<bool> {
        self.node_owner.iter().map(|&o| o == self.rank).collect()
    }

    pub fn n_owned_nodes(&self) -> usize {
        self.node_owner.iter().filter(|&&o| o == self.rank).count()
    }

    /// Bulk elements, owned and ghost.
    pub fn bulk_elements(&self) -> &[PartElement] {
        &self.elements[..self.header.n_bulk]
    }

    pub fn boundary_elements(&self) -> &[PartElement] {
        &self.elements[self.header.n_bulk..]
    }

    pub fn n_owned_bulk(&self) -> usize {
        self.bulk_elements().iter().filter(|e| e.owner == self.rank).count()
    }
}

impl MeshView for PartitionMesh {
    fn dim(&self) -> usize {
        self.dim
    }

    fn n_nodes(&self) -> usize {
        self.node_gid.len()
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }

    fn n_bulk(&self) -> usize {
        self.header.n_bulk
    }

    fn n_boundary(&self) -> usize {
        self.header.n_boundary
    }

    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let el = &self.elements[e.0 as usize];
        ElementRef {
            etype: el.etype,
            nodes: &el.nodes,
            tag: el.tag,
            parents: self.parents[e.0 as usize],
        }
    }
}

fn io_err(path: &Path) -> impl FnOnce(std::io::Error) -> PartError + '_ {
    move |source| PartError::Io { path: path.to_path_buf(), source }
}

/// Non-empty lines with their 1-based line numbers.
fn lines(path: &Path) -> Result<Vec<(usize, String)>, PartError> {
    let text = fs::read_to_string(path).map_err(io_err(path))?;
    Ok(text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(i, l)| (i + 1, l.to_string()))
        .collect())
}

fn parse<T: std::str::FromStr>(
    path: &Path,
    line: usize,
    tok: Option<&str>,
    what: &str,
) -> Result<T, PartError> {
    let bad = |msg: String| PartError::Parse { path: path.to_path_buf(), line, msg };
    let tok = tok.ok_or_else(|| bad(format!("missing {what}")))?;
    tok.parse().map_err(|_| bad(format!("bad {what} {tok:?}")))
}

fn read_header(path: &Path) -> Result<PartHeader, PartError> {
    let ls = lines(path)?;
    let toks = |i: usize| -> Result<(usize, Vec<&str>), PartError> {
        let (n, l) = ls.get(i).ok_or_else(|| PartError::Parse {
            path: path.to_path_buf(),
            line: ls.last().map_or(0, |l| l.0) + 1,
            msg: "header ends early".into(),
        })?;
        Ok((*n, l.split_whitespace().collect()))
    };
    let (n, t) = toks(0)?;
    let n_nodes = parse(path, n, t.first().copied(), "node count")?;
    let n_bulk = parse(path, n, t.get(1).copied(), "element count")?;
    let n_boundary = parse(path, n, t.get(2).copied(), "boundary count")?;
    let (n, t) = toks(1)?;
    let n_types: usize = parse(path, n, t.first().copied(), "type count")?;
    let mut element_types = Vec::with_capacity(n_types);
    for i in 0..n_types {
        let (n, t) = toks(2 + i)?;
        element_types.push((
            parse(path, n, t.first().copied(), "element type")?,
            parse(path, n, t.get(1).copied(), "type count")?,
        ));
    }
    let (n, t) = toks(2 + n_types)?;
    Ok(PartHeader {
        n_nodes,
        n_bulk,
        n_boundary,
        element_types,
        n_shared_nodes: parse(path, n, t.first().copied(), "shared node count")?,
        n_extra: parse(path, n, t.get(1).copied(), "second shared count")?,
    })
}

fn read_nodes(path: &Path) -> Result<(Vec<u64>, Vec<[f64; 3]>), PartError> {
    let mut gids = Vec::new();
    let mut coords = Vec::new();
    for (n, l) in lines(path)? {
        let t: Vec<&str> = l.split_whitespace().collect();
        gids.push(parse(path, n, t.first().copied(), "node id")?);
        coords.push([
            parse(path, n, t.get(2).copied(), "x")?,
            parse(path, n, t.get(3).copied(), "y")?,
            parse(path, n, t.get(4).copied(), "z")?,
        ]);
    }
    Ok((gids, coords))
}

fn read_shared(
    path: &Path,
    rank: usize,
    n_parts: usize,
) -> Result<Vec<(u64, Vec<usize>)>, PartError> {
    let mut out = Vec::new();
    for (n, l) in lines(path)? {
        let t: Vec<&str> = l.split_whitespace().collect();
        let gid = parse(path, n, t.first().copied(), "node id")?;
        let count: usize = parse(path, n, t.get(1).copied(), "sharing count")?;
        let mut parts = Vec::with_capacity(count);
        for k in 0..count {
            let p: usize = parse(path, n, t.get(2 + k).copied(), "partition")?;
            if p == 0 || p > n_parts {
                return Err(PartError::Parse {
                    path: path.to_path_buf(),
                    line: n,
                    msg: format!("partition {p} out of range"),
                });
            }
            parts.push(p - 1);
        }
        if !parts.contains(&rank) {
            return Err(PartError::Parse {
                path: path.to_path_buf(),
                line: n,
                msg: "sharing list omits this partition".into(),
            });
        }
        out.push((gid, parts));
    }
    Ok(out)
}

/// Bulk (`id body type nodes..`) or boundary (`id tag p1 p2 type nodes..`).
fn read_elements(
    path: &Path,
    boundary: bool,
    rank: usize,
    local_of: &HashMap<u64, u32>,
) -> Result<Vec<PartElement>, PartError> {
    let mut out = Vec::new();
    for (n, l) in lines(path)? {
        let t: Vec<&str> = l.split_whitespace().collect();
        let bad = |msg: String| PartError::Parse { path: path.to_path_buf(), line: n, msg };
        // `id` or `id/owner` (halo element, owner partition 1-based).
        let id_tok = *t.first().ok_or_else(|| bad("empty line".into()))?;
        let (gid, owner) = match id_tok.split_once('/') {
            Some((g, o)) => {
                let o: usize = parse(path, n, Some(o), "owner partition")?;
                (
                    parse(path, n, Some(g), "element id")?,
                    o.checked_sub(1).ok_or_else(|| bad("owner partition 0".into()))?,
                )
            }
            None => (parse(path, n, Some(id_tok), "element id")?, rank),
        };
        let tag = parse(path, n, t.get(1).copied(), "body/tag")?;
        let (parent_gids, type_at) = if boundary {
            let p = |i| -> Result<Option<u64>, PartError> {
                let g: u64 = parse(path, n, t.get(i).copied(), "parent element")?;
                Ok((g != 0).then_some(g))
            };
            ([p(2)?, p(3)?], 4)
        } else {
            ([None, None], 2)
        };
        let etype = ElementType(parse(path, n, t.get(type_at).copied(), "element type")?);
        if t.len() != type_at + 1 + etype.n_nodes() {
            return Err(bad(format!(
                "type {} needs {} nodes, line has {}",
                etype.0,
                etype.n_nodes(),
                t.len() - type_at - 1
            )));
        }
        let mut nodes = Vec::with_capacity(etype.n_nodes());
        for tok in &t[type_at + 1..] {
            let g: u64 = parse(path, n, Some(tok), "node id")?;
            nodes.push(
                *local_of.get(&g).ok_or_else(|| bad(format!("node {g} not in this partition")))?,
            );
        }
        out.push(PartElement { gid, owner, tag, etype, nodes, parent_gids });
    }
    Ok(out)
}
