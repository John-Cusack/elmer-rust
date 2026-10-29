//! `Mesh2D` (meshgen2d): the 2D mesh generator behind the `Mesh2D` prestep
//! of the adaptive-refinement quick tests, ported for the layer types they
//! use, `VoronoiVertex` and `MovingFront` with a `Delaunay` (or `Explicit`)
//! background mesh, and `Linear` or `Parabolic` elements.
//!
//! The port follows the C++ closely, down to the vertex recycling order and
//! the `std::random_shuffle` draws of glibc's `rand()` (the C++ never seeds
//! it), because both decide the element and node order of the output.
//!
//! Not ported (refused with [`MeshOpsError::Unsupported`]): the grid and
//! `Connect`/`BoundaryMesh` layers, `SSSFMovingFront`/`SSMFMovingFront`,
//! `Grid`/`DenseGrid`/`SparseGrid` background meshes and `External`
//! background meshes.

mod geom;
mod layers;
mod parse;
mod tri;

use std::collections::BTreeMap;

use crate::MeshOpsError;
use geom::*;
use parse::*;
use tri::*;

pub use geom::fmt_g16;

/// The four files Mesh2D writes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mesh2dFiles {
    pub header: String,
    pub nodes: String,
    pub elements: String,
    pub boundary: String,
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum EdgeType {
    Virtual,
    Outer,
    Inner,
}

struct Bel {
    edge: i32,
    a: NodeId,
    b: NodeId,
    c: Option<NodeId>,
    left: i32,
    right: i32,
    flipped: bool,
    tag: i32,
}

struct GEdge {
    dots: Vec<NodeId>,
    nodes: Vec<NodeId>,
    segments: i32,
    boundary_tag: i32,
    ty: EdgeType,
    bels: Vec<usize>,
    bgs: Vec<usize>,
}

struct Loop {
    direction: Vec<i32>,
    edges: Vec<i32>,
    nodes: Vec<NodeId>,
}

struct Layer {
    /// The body's tag.
    tag: i32,
    kind: LayerKind,
    d: Delaunay,
    bg: usize,
    loops: Vec<Loop>,
    /// `Layer::nodes`: fixed nodes, then the geometry nodes of the border.
    geometry_nodes: Vec<NodeId>,
    fixed_geometry_nodes: Vec<NodeId>,
}

struct Elem {
    tag: i32,
    body: i32,
    nodes: Vec<NodeId>,
}

struct Body {
    parabolic: bool,
    layers: Vec<usize>,
}

struct State {
    nodes: Nodes,
    rng: GlibcRand,
    geometry_nodes: BTreeMap<i32, NodeId>,
    edges: BTreeMap<i32, GEdge>,
    bodies: BTreeMap<i32, Body>,
    layers: Vec<Layer>,
    bgs: Vec<Delaunay>,
    bels: Vec<Bel>,
    next_bel_tag: i32,
    next_element_tag: i32,
    fixed_nodes: BTreeMap<i32, NodeId>,
    mesh_nodes: BTreeMap<i32, NodeId>,
    elements: Vec<Elem>,
    boundary_elements: Vec<usize>,
}

fn invalid(m: impl Into<String>) -> MeshOpsError {
    MeshOpsError::Invalid(m.into())
}

/// Runs Mesh2D on the text of a `.mif` file.
pub fn mesh2d(mif: &str) -> Result<Mesh2dFiles, MeshOpsError> {
    let parsed = parse::parse(mif, None).map_err(|e| {
        if e.contains("unsupported") || e.contains("not supported") {
            MeshOpsError::Unsupported(format!("Mesh2D: {e}"))
        } else {
            invalid(format!("Mesh2D: {e}"))
        }
    })?;
    let mut st = State {
        nodes: Nodes::new(),
        rng: GlibcRand::new(),
        geometry_nodes: BTreeMap::new(),
        edges: BTreeMap::new(),
        bodies: BTreeMap::new(),
        layers: Vec::new(),
        bgs: Vec::new(),
        bels: Vec::new(),
        next_bel_tag: 1,
        next_element_tag: 1,
        fixed_nodes: BTreeMap::new(),
        mesh_nodes: BTreeMap::new(),
        elements: Vec::new(),
        boundary_elements: Vec::new(),
    };
    st.convert_edge_format(&parsed)?;
    st.discretize()?;
    Ok(st.output(&parsed))
}

impl State {
    /// `Mesh::convertEdgeFormat`.
    fn convert_edge_format(&mut self, parser: &Parsed) -> Result<(), MeshOpsError> {
        let mut node_boundary: BTreeMap<i32, i32> = BTreeMap::new();
        for nt in &parser.nodes {
            let id = self.nodes.geometry_node(nt.tag, nt.x, nt.y);
            self.geometry_nodes.insert(nt.tag, id);
            node_boundary.insert(nt.tag, nt.boundary_tag);
        }

        let mut node_edges: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (ei, et) in parser.edges.iter().enumerate() {
            let mut dots = Vec::new();
            for &n in &et.nodes {
                node_edges.entry(n).or_default().push(ei);
                let g = *self
                    .geometry_nodes
                    .get(&n)
                    .ok_or_else(|| invalid(format!("Nonexisting node in edge {}", et.tag)))?;
                dots.push(g);
            }
            self.edges.insert(
                et.tag,
                GEdge {
                    dots,
                    nodes: Vec::new(),
                    segments: et.seg_count,
                    boundary_tag: et.boundary_tag,
                    ty: EdgeType::Outer,
                    bels: Vec::new(),
                    bgs: Vec::new(),
                },
            );
        }

        let mut node_layers: BTreeMap<i32, Vec<f64>> = BTreeMap::new();
        let mut edge_body: BTreeMap<i32, i32> = BTreeMap::new();
        for bt in &parser.bodies {
            let mut body = Body {
                parabolic: bt.parabolic,
                layers: Vec::new(),
            };
            for ly in &bt.layers {
                let ldelta = if ly.delta <= 0.0 { bt.delta } else { ly.delta };

                let mut bg = Delaunay::new(false);
                match ly.bg.kind.as_str() {
                    "Delaunay" => {}
                    "Explicit" => {
                        let mut bgnodes = Vec::new();
                        for &(x, y, d) in &ly.bg.nodes {
                            let id = self.nodes.geometry_node(0, x, y);
                            self.nodes.v[id].delta = d * parser.scale;
                            bgnodes.push(id);
                        }
                        bg.bg_initialize(&mut self.nodes, &mut self.rng, &bgnodes);
                    }
                    other => {
                        return Err(MeshOpsError::Unsupported(format!(
                            "Mesh2D: background mesh type {other}"
                        )))
                    }
                }
                self.bgs.push(bg);
                let bg_index = self.bgs.len() - 1;

                let mut layer = Layer {
                    tag: bt.tag,
                    kind: ly.kind,
                    d: Delaunay::new(ly.kind == LayerKind::VoronoiSegment),
                    bg: bg_index,
                    loops: Vec::new(),
                    geometry_nodes: Vec::new(),
                    fixed_geometry_nodes: Vec::new(),
                };

                for nt in &ly.fixed_nodes {
                    let g = self.nodes.geometry_node(nt.tag, nt.x, nt.y);
                    let delta = if nt.delta <= 0.0 {
                        if ldelta > 0.0 {
                            ldelta
                        } else {
                            parser.delta
                        }
                    } else {
                        nt.delta
                    };
                    self.nodes.v[g].delta = delta * parser.scale;
                    layer.geometry_nodes.push(g);
                    let m = self.nodes.mesh_node_from_geometry(g);
                    layer.fixed_geometry_nodes.push(m);
                }

                for lt in &ly.loops {
                    let mut lp = Loop {
                        direction: Vec::new(),
                        edges: Vec::new(),
                        nodes: Vec::new(),
                    };
                    let len = lt.edges.len();
                    for k in 0..len {
                        let edge_tag = if lt.direction > 0 {
                            lt.edges[k]
                        } else {
                            -lt.edges[len - 1 - k]
                        };
                        let dir = if edge_tag > 0 { 1 } else { -1 };
                        let edge_tag = edge_tag.abs();
                        let ed = self.edges.get_mut(&edge_tag).ok_or_else(|| {
                            invalid(format!("Nonexisting edge in layer {}", ly.tag))
                        })?;
                        match edge_body.get(&edge_tag) {
                            None => {
                                edge_body.insert(edge_tag, bt.tag);
                            }
                            Some(&b) if b == bt.tag => ed.ty = EdgeType::Virtual,
                            Some(_) => ed.ty = EdgeType::Inner,
                        }
                        ed.bgs.push(bg_index);

                        let dots: Vec<NodeId> = if dir > 0 {
                            ed.dots.clone()
                        } else {
                            ed.dots.iter().rev().copied().collect()
                        };
                        for &g in &dots[..dots.len() - 1] {
                            node_layers
                                .entry(self.nodes.v[g].tag)
                                .or_default()
                                .push(ldelta);
                        }
                        lp.direction.push(dir);
                        lp.edges.push(edge_tag);
                        lp.nodes.push(if dir == 1 {
                            ed.dots[0]
                        } else {
                            *ed.dots.last().unwrap()
                        });
                    }
                    layer.loops.push(lp);
                }

                self.layers.push(layer);
                body.layers.push(self.layers.len() - 1);
            }
            self.bodies.insert(bt.tag, body);
        }

        for nt in &parser.nodes {
            let mut delta = nt.delta;
            if delta <= 0.0 {
                let mut mean = 1.0;
                let mut n = 0;
                if let Some(es) = node_edges.get(&nt.tag) {
                    for &e in es {
                        if parser.edges[e].delta > 0.0 {
                            mean *= parser.edges[e].delta;
                            n += 1;
                        }
                    }
                }
                if n == 0 {
                    if let Some(ls) = node_layers.get(&nt.tag) {
                        for &d in ls {
                            if d > 0.0 {
                                mean *= d;
                                n += 1;
                            }
                        }
                    }
                }
                delta = if n > 0 {
                    mean.powf(1.0 / n as f64)
                } else {
                    parser.delta
                };
            }
            self.nodes.v[self.geometry_nodes[&nt.tag]].delta = delta * parser.scale;
        }

        for li in 0..self.layers.len() {
            self.initialize_layer(li);
        }
        Ok(())
    }

    /// `Connect::initialize`.
    fn initialize_layer(&mut self, li: usize) {
        let mut chain = Vec::new();
        for lp in &self.layers[li].loops {
            for (k, &e) in lp.edges.iter().enumerate() {
                let ed = &self.edges[&e];
                let mut dots: Vec<NodeId> = if lp.direction[k] > 0 {
                    ed.dots.clone()
                } else {
                    ed.dots.iter().rev().copied().collect()
                };
                dots.pop();
                chain.extend(dots);
            }
        }
        self.layers[li].geometry_nodes.extend(chain);
        let bg = self.layers[li].bg;
        if !self.bgs[bg].initialized {
            let nodes = self.layers[li].geometry_nodes.clone();
            self.bgs[bg].bg_initialize(&mut self.nodes, &mut self.rng, &nodes);
        }
    }

    /// `Mesh::discretize`.
    fn discretize(&mut self) -> Result<(), MeshOpsError> {
        let gtags: Vec<(i32, NodeId)> = self.geometry_nodes.iter().map(|(&t, &n)| (t, n)).collect();
        for (tag, g) in gtags {
            let m = self.nodes.mesh_node_from_geometry(g);
            self.nodes.v[m].boundarynode = true;
            self.fixed_nodes.insert(tag, m);
        }

        let etags: Vec<i32> = self.edges.keys().copied().collect();
        for constant in [true, false] {
            for &t in &etags {
                if (self.edges[&t].segments > 0) == constant {
                    self.discretize_edge(t);
                }
            }
        }

        let btags: Vec<i32> = self.bodies.keys().copied().collect();
        for bt in btags {
            let layer_ids = self.bodies[&bt].layers.clone();
            let mut elements = Vec::new();
            for li in layer_ids {
                self.discretize_layer(li, &mut elements)?;
            }
            for e in &mut elements {
                e.body = bt;
            }
            self.elements.extend(elements);
        }

        for &t in &etags {
            let mut strip = Vec::new();
            self.edge_elements(t, &mut strip, 0);
            self.boundary_elements.extend(strip);
        }

        self.create_middle_nodes();
        Ok(())
    }

    /// `Layer::discretize` of `VoronoiVertex` and `VoronoiSegment`.
    fn discretize_layer(
        &mut self,
        li: usize,
        elements: &mut Vec<Elem>,
    ) -> Result<(), MeshOpsError> {
        // Take the layer out so its triangulation can use the shared state.
        let mut layer = std::mem::replace(
            &mut self.layers[li],
            Layer {
                tag: 0,
                kind: LayerKind::VoronoiVertex,
                d: Delaunay::new(false),
                bg: 0,
                loops: Vec::new(),
                geometry_nodes: Vec::new(),
                fixed_geometry_nodes: Vec::new(),
            },
        );
        let result = self.run_layer(&mut layer, elements);
        self.layers[li] = layer;
        result
    }

    fn run_layer(
        &mut self,
        layer: &mut Layer,
        elements: &mut Vec<Elem>,
    ) -> Result<(), MeshOpsError> {
        // collectNodesAndPairs
        let mut border = Vec::new();
        let mut links = std::collections::BTreeSet::new();
        for lp in &layer.loops {
            let mut strip = Vec::new();
            for (k, &e) in lp.edges.iter().enumerate() {
                let ed = &self.edges[&e];
                if lp.direction[k] < 0 {
                    strip.extend(ed.nodes.iter().rev().copied());
                } else {
                    strip.extend(ed.nodes.iter().copied());
                }
                strip.pop();
            }
            let len = strip.len();
            for i in 0..len {
                let t1 = self.nodes.v[strip[i]].tag;
                let t2 = self.nodes.v[strip[(i + 1) % len]].tag;
                links.insert((t1.min(t2), t1.max(t2)));
            }
            border.extend(strip);
        }
        layer.d.border = border;
        layer.d.links = links;

        let tag = layer.tag;
        layer
            .d
            .triangulate(
                &mut self.nodes,
                &mut self.rng,
                &layer.fixed_geometry_nodes,
                tag,
            )
            .map_err(|e| invalid(format!("Mesh2D: {e}")))?;
        let mut actives = Pq::new();
        layer
            .d
            .remove_boundary_connectors(&mut self.nodes, &mut actives);
        let bg = &self.bgs[layer.bg];
        match layer.kind {
            LayerKind::VoronoiVertex => layers::generate_vertex(&mut layer.d, &mut self.nodes, bg),
            LayerKind::VoronoiSegment => {
                layers::generate_segment(&mut layer.d, &mut self.nodes, bg, tag)
            }
        }

        // exportNodes
        let all = layer.d.all_vertices.clone();
        for &v in &all {
            let vert = &layer.d.verts[v];
            if !(vert.deleted || vert.is_external()) {
                let tag = self.next_element_tag;
                self.next_element_tag += 1;
                layer.d.verts[v].element_tag = tag;
                let nodes = layer.d.verts[v].nodes.to_vec();
                for &n in &nodes {
                    self.mesh_nodes.insert(self.nodes.v[n].tag, n);
                }
                elements.push(Elem {
                    tag,
                    body: 0,
                    nodes,
                });
            }
        }
        let mut bels = Vec::new();
        for lp in &layer.loops {
            for (k, &e) in lp.edges.iter().enumerate() {
                if self.edges[&e].ty != EdgeType::Virtual {
                    self.edge_elements(e, &mut bels, lp.direction[k]);
                }
            }
        }
        let mut v = layer.d.root;
        for b in bels {
            let (a, bb) = (self.bels[b].a, self.bels[b].b);
            let flipped = self.bels[b].flipped;
            let holder = if !flipped {
                layer.d.vertex_with(&self.nodes, v, a, bb)
            } else {
                layer.d.vertex_with(&self.nodes, v, bb, a)
            }
            .ok_or_else(|| invalid("Mesh2D: boundary element without a triangle"))?;
            let id = layer.d.verts[holder].element_tag;
            if !flipped {
                self.bels[b].left = id;
            } else {
                self.bels[b].right = id;
            }
            v = holder;
        }
        Ok(())
    }

    /// `GeometryEdge::elements`.
    fn edge_elements(&mut self, tag: i32, strip: &mut Vec<usize>, direction: i32) {
        if self.edges[&tag].bels.is_empty() {
            let nodes = self.edges[&tag].nodes.clone();
            let bt = self.edges[&tag].boundary_tag;
            let mut made = Vec::new();
            for i in 0..nodes.len().saturating_sub(1) {
                self.bels.push(Bel {
                    edge: bt,
                    a: nodes[i],
                    b: nodes[i + 1],
                    c: None,
                    left: 0,
                    right: 0,
                    flipped: false,
                    tag: self.next_bel_tag,
                });
                self.next_bel_tag += 1;
                made.push(self.bels.len() - 1);
            }
            self.edges.get_mut(&tag).unwrap().bels = made;
        }
        let bels = self.edges[&tag].bels.clone();
        for &b in &bels {
            self.bels[b].flipped = direction < 0;
        }
        if direction < 0 {
            strip.extend(bels.iter().rev().copied());
        } else {
            strip.extend(bels.iter().copied());
        }
    }

    fn edge_interpolate(&self, tag: i32, x: f64, y: f64) -> f64 {
        let mut value = 0.0;
        for &b in &self.edges[&tag].bgs {
            let h = self.bgs[b].interpolate(&self.nodes, x, y);
            if value == 0.0 || h < value {
                value = h;
            }
        }
        value
    }

    /// `GeometryEdge::discretize`.
    fn discretize_edge(&mut self, tag: i32) {
        let dots = self.edges[&tag].dots.clone();
        let len = dots.len() - 1;
        if self.edges[&tag].segments == 0 {
            for i in 0..len {
                self.discretize_graded(tag, dots[i], dots[i + 1]);
                if i < len - 1 {
                    self.edges.get_mut(&tag).unwrap().nodes.pop();
                }
            }
        } else {
            let mut segments = self.edges[&tag].segments;
            if (segments as usize) < len {
                segments = len as i32;
            }
            self.edges.get_mut(&tag).unwrap().segments = segments;
            let lengths: Vec<f64> = (0..len)
                .map(|i| {
                    let (a, b) = (&self.nodes.v[dots[i]], &self.nodes.v[dots[i + 1]]);
                    distance(a.x, a.y, b.x, b.y)
                })
                .collect();
            let mut segs = vec![1i32; len];
            let mut remaining = segments - len as i32;
            while remaining > 0 {
                let mut longest = 0;
                let mut length = lengths[0] / segs[0] as f64;
                for i in 1..len {
                    let l = lengths[i] / segs[i] as f64;
                    if l > length {
                        length = l;
                        longest = i;
                    }
                }
                segs[longest] += 1;
                remaining -= 1;
            }
            for i in 0..len {
                self.discretize_constant(tag, segs[i], dots[i], dots[i + 1]);
                if i < len - 1 {
                    self.edges.get_mut(&tag).unwrap().nodes.pop();
                }
            }
        }
    }

    fn discretize_constant(&mut self, tag: i32, n_seg: i32, from: NodeId, to: NodeId) {
        let first = self.fixed_nodes[&self.nodes.v[from].tag];
        let last = self.fixed_nodes[&self.nodes.v[to].tag];
        let dx = self.nodes.v[last].x - self.nodes.v[first].x;
        let dy = self.nodes.v[last].y - self.nodes.v[first].y;
        let length = (dx * dx + dy * dy).sqrt();
        let sx = dx / length;
        let sy = dy / length;
        let step = length / n_seg as f64;
        let bx = self.nodes.v[first].x;
        let by = self.nodes.v[first].y;
        self.nodes.v[from].delta = step;
        self.nodes.v[to].delta = step;
        let mut out = vec![first];
        for i in 1..n_seg {
            let nd = self
                .nodes
                .mesh_node(bx + i as f64 * step * sx, by + i as f64 * step * sy);
            self.nodes.v[nd].boundarynode = true;
            out.push(nd);
        }
        out.push(last);
        self.edges.get_mut(&tag).unwrap().nodes.extend(out);
    }

    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn discretize_graded(&mut self, tag: i32, from: NodeId, to: NodeId) {
        let (fx, fy, tx, ty) = (
            self.nodes.v[from].x,
            self.nodes.v[from].y,
            self.nodes.v[to].x,
            self.nodes.v[to].y,
        );
        let line_length = distance(fx, fy, tx, ty);
        let sx = (tx - fx) / line_length;
        let sy = (ty - fy) / line_length;
        let limit = line_length;
        let mut at = 0.0;
        let mut predicted = self.edge_interpolate(tag, fx, fy);
        let mut lengths: Vec<f64> = Vec::new();
        loop {
            loop {
                let prev = predicted;
                let mut center = at + predicted / 2.0;
                if center > line_length {
                    center = line_length;
                }
                predicted = 0.5
                    * (predicted + self.edge_interpolate(tag, fx + center * sx, fy + center * sy));
                if !(((prev - predicted) / line_length).abs() > 0.0000001) {
                    break;
                }
            }
            at += predicted;
            if at > limit - 0.01 * predicted {
                break;
            }
            lengths.push(at);
        }
        let mult = line_length / at;
        for l in &mut lengths {
            *l *= mult;
        }
        let first = self.fixed_nodes[&self.nodes.v[from].tag];
        let last = self.fixed_nodes[&self.nodes.v[to].tag];
        let mut out = vec![first];
        for &l in &lengths {
            let nd = self.nodes.mesh_node(l * sx + fx, l * sy + fy);
            self.nodes.v[nd].boundarynode = true;
            out.push(nd);
        }
        out.push(last);
        self.edges.get_mut(&tag).unwrap().nodes.extend(out);
    }

    /// `Mesh::createMiddleNodes`.
    fn create_middle_nodes(&mut self) {
        let mut edge_node: BTreeMap<(i32, i32), NodeId> = BTreeMap::new();
        for ei in 0..self.elements.len() {
            let body = self.elements[ei].body;
            if !self.bodies[&body].parabolic {
                continue;
            }
            let size = self.elements[ei].nodes.len();
            let mut new_nodes = Vec::new();
            for i in 0..size {
                let n0 = self.elements[ei].nodes[i];
                let n1 = self.elements[ei].nodes[(i + 1) % size];
                let (t0, t1) = (self.nodes.v[n0].tag, self.nodes.v[n1].tag);
                let key = (t0.min(t1), t0.max(t1));
                let node = match edge_node.get(&key) {
                    Some(&n) => n,
                    None => {
                        let x = (self.nodes.v[n0].x + self.nodes.v[n1].x) / 2.0;
                        let y = (self.nodes.v[n0].y + self.nodes.v[n1].y) / 2.0;
                        let n = self.nodes.mesh_node(x, y);
                        edge_node.insert(key, n);
                        self.mesh_nodes.insert(self.nodes.v[n].tag, n);
                        n
                    }
                };
                new_nodes.push(node);
            }
            self.elements[ei].nodes.extend(new_nodes);
        }
        for &b in &self.boundary_elements {
            let (t0, t1) = (
                self.nodes.v[self.bels[b].a].tag,
                self.nodes.v[self.bels[b].b].tag,
            );
            if let Some(&n) = edge_node.get(&(t0.min(t1), t0.max(t1))) {
                self.bels[b].c = Some(n);
            }
        }
    }

    /// `Mesh::outputHeader`, `outputNodes`, `outputElements` and
    /// `outputBoundary`, in the order `main` writes them (the node pass
    /// renumbers the nodes the later ones print).
    fn output(&mut self, _parsed: &Parsed) -> Mesh2dFiles {
        let mut ngnds = 0;
        let mut boundary_tags: BTreeMap<i32, i32> = BTreeMap::new();
        for (&tag, &g) in &self.geometry_nodes {
            let _ = g;
            boundary_tags.insert(tag, 0);
        }
        // The boundary tag of a geometry node is the second number of its
        // NodeId line; the parser kept it in the node tokens.
        for nt in &_parsed.nodes {
            boundary_tags.insert(nt.tag, nt.boundary_tag);
        }
        for &b in boundary_tags.values() {
            if b > 0 {
                ngnds += 1;
            }
        }

        let mut header = String::new();
        header.push_str(&format!(
            "{} {} {}\n",
            self.mesh_nodes.len(),
            self.elements.len(),
            self.boundary_elements.len() + ngnds
        ));
        let mut n: BTreeMap<i32, usize> = BTreeMap::new();
        for e in &self.elements {
            *n.entry(300 + e.nodes.len() as i32).or_default() += 1;
        }
        for &b in &self.boundary_elements {
            *n.entry(if self.bels[b].c.is_some() { 203 } else { 202 })
                .or_default() += 1;
        }
        *n.entry(101).or_default() += ngnds;
        header.push_str(&format!("{}\n", n.len()));
        for (k, v) in &n {
            header.push_str(&format!("{k} {v}\n"));
        }

        let mut nodes = String::new();
        for (tag, &id) in (1..).zip(self.mesh_nodes.values()) {
            self.nodes.v[id].tag = tag;
            let nd = &self.nodes.v[id];
            nodes.push_str(&format!(
                "{} -1 {} {} 0\n",
                nd.tag,
                fmt_g16(nd.x),
                fmt_g16(nd.y)
            ));
        }

        let mut elements = String::new();
        for e in &self.elements {
            elements.push_str(&format!("{} {} {}", e.tag, e.body, 300 + e.nodes.len()));
            for &n in &e.nodes {
                elements.push_str(&format!(" {}", self.nodes.v[n].tag));
            }
            elements.push('\n');
        }

        let mut boundary = String::new();
        for &b in &self.boundary_elements {
            let be = &self.bels[b];
            boundary.push_str(&format!("{} {} {} {}", be.tag, be.edge, be.left, be.right));
            match be.c {
                Some(c) => boundary.push_str(&format!(
                    " 203 {} {} {}\n",
                    self.nodes.v[be.a].tag, self.nodes.v[be.b].tag, self.nodes.v[c].tag
                )),
                None => boundary.push_str(&format!(
                    " 202 {} {}\n",
                    self.nodes.v[be.a].tag, self.nodes.v[be.b].tag
                )),
            }
        }
        let mut id = self.boundary_elements.len() + 1;
        for (&tag, &bt) in &boundary_tags {
            if bt > 0 {
                boundary.push_str(&format!("{id} {bt} -1 -1 101 {tag}\n"));
                id += 1;
            }
        }
        Mesh2dFiles {
            header,
            nodes,
            elements,
            boundary,
        }
    }
}
