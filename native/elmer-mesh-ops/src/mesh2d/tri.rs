//! The triangulation machinery of Mesh2D: `Vertex` (a triangle with its
//! circumcircle and neighbours), `Connect` (incremental Delaunay insertion
//! with the vertex recycling of the C++, which decides element order),
//! `BGTriangleMesh`/`BGVertex` (the interpolated size field) and the
//! priority queue `pq`.
//!
//! The C++ objects point at each other; here nodes and vertices live in
//! arenas and refer to each other by index.

use std::collections::{BTreeSet, VecDeque};

use super::geom::*;

pub type NodeId = usize;

#[derive(Clone, Debug)]
pub struct Node {
    pub tag: i32,
    pub x: f64,
    pub y: f64,
    pub boundarynode: bool,
    pub fixed: bool,
    /// `GeometryNode::delta`, the requested element size.
    pub delta: f64,
}

/// All nodes, `MeshNode` and `GeometryNode` alike, and the `MeshNode` tag
/// counter (a function-local static of MeshNode.cpp).
pub struct Nodes {
    pub v: Vec<Node>,
    next_tag: i32,
}

impl Nodes {
    pub fn new() -> Self {
        Nodes {
            v: Vec::new(),
            next_tag: 1,
        }
    }

    fn push(&mut self, tag: i32, x: f64, y: f64) -> NodeId {
        self.v.push(Node {
            tag,
            x,
            y,
            boundarynode: false,
            fixed: false,
            delta: 0.0,
        });
        self.v.len() - 1
    }

    /// `MeshNode(x, y)`.
    pub fn mesh_node(&mut self, x: f64, y: f64) -> NodeId {
        let tag = self.next_tag;
        self.next_tag += 1;
        self.push(tag, x, y)
    }

    /// `MeshNode(tag, x, y)`.
    pub fn mesh_node_tagged(&mut self, tag: i32, x: f64, y: f64) -> NodeId {
        self.push(tag, x, y)
    }

    /// `MeshNode(const GeometryNode&)`: keeps the tag and moves the counter
    /// past it.
    pub fn mesh_node_from_geometry(&mut self, g: NodeId) -> NodeId {
        let (tag, x, y) = (self.v[g].tag, self.v[g].x, self.v[g].y);
        if tag >= self.next_tag {
            self.next_tag = tag + 1;
        }
        self.push(tag, x, y)
    }

    /// `GeometryNode(tag, x, y)`.
    pub fn geometry_node(&mut self, tag: i32, x: f64, y: f64) -> NodeId {
        self.push(tag, x, y)
    }
}

impl Default for Nodes {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug)]
pub struct Vertex {
    pub nodes: [NodeId; 3],
    pub vertices: [Option<usize>; 3],
    pub vx: f64,
    pub vy: f64,
    pub radius: f64,
    pub ratio: f64,
    pub deleted: bool,
    pub heap: i32,
    pub body: i32,
    /// `Element::elementTag`, set when the vertex is exported as an element.
    pub element_tag: i32,
    /// `VSVertex::type`.
    pub vs_type: VsType,
    /// `BGVertex::coeff`.
    pub coeff: [f64; 3],
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VsType {
    Accepted,
    Active,
    Waiting,
}

impl Vertex {
    fn blank() -> Self {
        Vertex {
            nodes: [0; 3],
            vertices: [None; 3],
            vx: 0.0,
            vy: 0.0,
            radius: 0.0,
            ratio: 0.0,
            deleted: false,
            heap: -1,
            body: 0,
            element_tag: -1,
            vs_type: VsType::Waiting,
            coeff: [0.0; 3],
        }
    }

    pub fn is_external(&self) -> bool {
        self.body == 0
    }

    /// `Vertex::rightSize`.
    pub fn right_size(&mut self, reference: f64) -> bool {
        const DELTA: f64 = 1.25;
        self.ratio = self.radius / reference;
        self.ratio < DELTA
    }
}

#[derive(Copy, Clone)]
struct Triple {
    from: NodeId,
    to: NodeId,
    vtx: Option<usize>,
}

/// `pq`: a max-heap of vertices on `orderingValue()`, which is the radius
/// for `VSVertex` and the ratio otherwise. Positions are kept in
/// `Vertex::heap`.
#[derive(Default)]
pub struct Pq {
    store: Vec<usize>,
    last: i32,
}

impl Pq {
    pub fn new() -> Self {
        Pq {
            store: Vec::new(),
            last: -1,
        }
    }

    pub fn size(&self) -> i32 {
        self.last + 1
    }

    fn ov(verts: &[Vertex], vs: bool, v: usize) -> f64 {
        if vs {
            verts[v].radius
        } else {
            verts[v].ratio
        }
    }

    fn set(&mut self, k: i32, v: usize) {
        if self.store.len() <= k as usize {
            self.store.resize(k as usize + 1, usize::MAX);
        }
        self.store[k as usize] = v;
    }

    pub fn first(&mut self, verts: &mut [Vertex], vs: bool) -> Option<usize> {
        if self.last < 0 {
            return None;
        }
        let v = self.store[0];
        if self.last > 0 {
            self.remove(verts, vs, v);
        } else {
            verts[v].heap = -1;
            self.last = -1;
        }
        Some(v)
    }

    pub fn remove(&mut self, verts: &mut [Vertex], vs: bool, vtx: usize) {
        let h = verts[vtx].heap;
        if h == self.last {
            verts[vtx].heap = -1;
            self.last -= 1;
            return;
        }
        verts[vtx].heap = -1;
        let moved = self.store[self.last as usize];
        self.set(h, moved);
        verts[moved].heap = h;
        self.last -= 1;
        if h > 0
            && Self::ov(verts, vs, self.store[((h - 1) / 2) as usize])
                < Self::ov(verts, vs, self.store[h as usize])
        {
            self.upheap(verts, vs, h);
        } else {
            self.downheap(verts, vs, h);
        }
    }

    pub fn insert(&mut self, verts: &mut [Vertex], vs: bool, vtx: usize) {
        self.last += 1;
        let last = self.last;
        self.set(last, vtx);
        verts[vtx].heap = last;
        self.upheap(verts, vs, last);
    }

    pub fn remove_relevants(&mut self, verts: &mut [Vertex], vs: bool, vl: &[usize]) {
        for &v in vl {
            if verts[v].heap != -1 {
                self.remove(verts, vs, v);
            }
        }
    }

    fn upheap(&mut self, verts: &mut [Vertex], vs: bool, p: i32) {
        let mut k = p;
        let v = self.store[p as usize];
        while k != 0
            && Self::ov(verts, vs, self.store[((k - 1) / 2) as usize]) <= Self::ov(verts, vs, v)
        {
            let parent = self.store[((k - 1) / 2) as usize];
            self.set(k, parent);
            verts[parent].heap = k;
            k = (k - 1) / 2;
        }
        self.set(k, v);
        verts[v].heap = k;
    }

    fn downheap(&mut self, verts: &mut [Vertex], vs: bool, p: i32) {
        let mut k = p;
        let v = self.store[p as usize];
        let n = self.last;
        while (k + 1) <= (n + 1) / 2 {
            let mut j = k + k + 1;
            if j < n
                && Self::ov(verts, vs, self.store[j as usize])
                    <= Self::ov(verts, vs, self.store[(j + 1) as usize])
            {
                j += 1;
            }
            if Self::ov(verts, vs, self.store[j as usize]) <= Self::ov(verts, vs, v) {
                break;
            }
            let child = self.store[j as usize];
            self.set(k, child);
            verts[child].heap = k;
            k = j;
        }
        self.set(k, v);
        verts[v].heap = k;
    }
}

/// One `Connect` object: the triangulation of a layer, or of a background
/// mesh (`BGTriangleMesh`).
pub struct Delaunay {
    /// Vertices are `VSVertex` (ordered by radius).
    pub vs: bool,
    pub verts: Vec<Vertex>,
    pub root: usize,
    pub world: [NodeId; 4],
    pub border: Vec<NodeId>,
    pub links: BTreeSet<(i32, i32)>,
    pub deleted: Vec<usize>,
    pub new_vertices: Vec<usize>,
    pub vertex_store: Vec<usize>,
    hull: Vec<Triple>,
    pub all_vertices: Vec<usize>,
    /// `BGMesh::initialized`.
    pub initialized: bool,
}

impl Delaunay {
    pub fn new(vs: bool) -> Self {
        Delaunay {
            vs,
            verts: Vec::new(),
            root: 0,
            world: [0; 4],
            border: Vec::new(),
            links: BTreeSet::new(),
            deleted: Vec::new(),
            new_vertices: Vec::new(),
            vertex_store: Vec::new(),
            hull: Vec::new(),
            all_vertices: Vec::new(),
            initialized: false,
        }
    }

    /// `Connect::getVertices`: recycled vertices first (last in, first out),
    /// then new ones.
    fn get_vertices(&mut self, count: usize) -> Vec<usize> {
        let mut v = Vec::new();
        let mut i = 0;
        while i < count {
            match self.vertex_store.pop() {
                Some(s) => v.push(s),
                None => break,
            }
            i += 1;
        }
        while i < count {
            self.verts.push(Vertex::blank());
            let id = self.verts.len() - 1;
            v.push(id);
            self.all_vertices.push(id);
            i += 1;
        }
        v
    }

    /// `Vertex::reset`.
    pub fn reset(&mut self, nodes: &Nodes, v: usize, n1: NodeId, n2: NodeId, n3: NodeId) {
        let (a, b, c) = (&nodes.v[n1], &nodes.v[n2], &nodes.v[n3]);
        let (vx, vy) = vertex_location(a.x, a.y, b.x, b.y, c.x, c.y);
        let vert = &mut self.verts[v];
        vert.nodes = [n1, n2, n3];
        vert.vx = vx;
        vert.vy = vy;
        vert.radius = distance(vx, vy, b.x, b.y);
        vert.deleted = false;
        vert.heap = -1;
        vert.body = 0;
    }

    /// `Connect::makeWorld`.
    pub fn make_world(&mut self, nodes: &mut Nodes) {
        let first = &nodes.v[self.border[0]];
        let (mut maxx, mut minx, mut maxy, mut miny) = (first.x, first.x, first.y, first.y);
        for &b in &self.border[1..] {
            let (tx, ty) = (nodes.v[b].x, nodes.v[b].y);
            if minx > tx {
                minx = tx;
            } else if maxx < tx {
                maxx = tx;
            }
            if miny > ty {
                miny = ty;
            } else if maxy < ty {
                maxy = ty;
            }
        }
        let dx = maxx - minx;
        let dy = maxy - miny;
        let cx = (maxx + minx) / 2.0;
        let cy = (maxy + miny) / 2.0;
        self.world[0] = nodes.mesh_node_tagged(-1, cx - dx, cy - dy);
        self.world[1] = nodes.mesh_node_tagged(-2, cx + dx, cy - dy);
        self.world[2] = nodes.mesh_node_tagged(-3, cx + dx, cy + dy);
        self.world[3] = nodes.mesh_node_tagged(-4, cx - dx, cy + dy);

        let bx = self.get_vertices(2);
        let w = self.world;
        self.reset(nodes, bx[0], w[0], w[1], w[2]);
        self.reset(nodes, bx[1], w[2], w[3], w[0]);
        self.verts[bx[0]].vertices = [None, None, Some(bx[1])];
        self.verts[bx[1]].vertices = [None, None, Some(bx[0])];
        self.root = bx[0];
    }

    /// `Vertex::firstAcceptableFound`.
    pub fn first_acceptable_found(&self, nodes: &Nodes, start: usize, tx: f64, ty: f64) -> usize {
        let s = [2usize, 3, 1];
        let e = [1usize, 2, 3];
        let mut curr = start;
        let mut prev: Option<usize> = None;
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        x[0] = tx;
        y[0] = ty;
        for k in 0..3 {
            let n = &nodes.v[self.verts[start].nodes[k]];
            x[k + 1] = n.x;
            y[k + 1] = n.y;
        }
        loop {
            let mut i = 0;
            while i < 3 {
                let nb = self.verts[curr].vertices[i];
                if let Some(n) = nb.filter(|_| {
                    nb != prev && orientation(tx, ty, x[s[i]], y[s[i]], x[e[i]], y[e[i]])
                }) {
                    prev = Some(curr);
                    curr = n;
                    x[0] = tx;
                    y[0] = ty;
                    for k in 0..3 {
                        let n = &nodes.v[self.verts[curr].nodes[k]];
                        x[k + 1] = n.x;
                        y[k + 1] = n.y;
                    }
                    break;
                }
                i += 1;
            }
            if i >= 3 {
                break;
            }
        }
        curr
    }

    /// `Vertex::isDestroyedBy`: 1 destroyed, 0 not, -1 an external vertex
    /// would be destroyed and that is not allowed.
    fn is_destroyed_by(&self, nodes: &Nodes, v: usize, tx: f64, ty: f64, ext_ok: bool) -> i32 {
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        x[0] = tx;
        y[0] = ty;
        for k in 0..3 {
            let n = &nodes.v[self.verts[v].nodes[k]];
            x[k + 1] = n.x;
            y[k + 1] = n.y;
        }
        if in_circle(&x, &y) != 0 {
            if self.verts[v].is_external() {
                if !ext_ok {
                    return -1;
                }
                return 1;
            }
            return 1;
        }
        0
    }

    /// `Vertex::vertexOwning`.
    pub fn vertex_owning(
        &self,
        nodes: &Nodes,
        start: usize,
        a: NodeId,
        b: NodeId,
    ) -> Option<usize> {
        let mut found = false;
        let mut holder = 0usize;
        let mut curr = start;
        let mut prev: Option<usize> = None;
        let tx = (nodes.v[a].x + nodes.v[b].x) / 2.0;
        let ty = (nodes.v[a].y + nodes.v[b].y) / 2.0;
        let s = [2usize, 3, 1];
        let e = [1usize, 2, 3];
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        x[0] = tx;
        y[0] = ty;
        for k in 0..3 {
            let n = &nodes.v[self.verts[start].nodes[k]];
            x[k + 1] = n.x;
            y[k + 1] = n.y;
        }
        for i in 0..3 {
            if self.verts[start].nodes[i] == a {
                found = true;
                holder = i;
                if self.verts[start].nodes[(i + 1) % 3] == b {
                    return Some(curr);
                }
                if self.verts[start].nodes[(i + 3 - 1) % 3] == b {
                    return Some(curr);
                }
                break;
            }
        }
        while !found {
            for i in 0..3 {
                let nb = self.verts[curr].vertices[i];
                if let Some(n) = nb.filter(|_| {
                    nb != prev && orientation(tx, ty, x[s[i]], y[s[i]], x[e[i]], y[e[i]])
                }) {
                    prev = Some(curr);
                    curr = n;
                    x[0] = tx;
                    y[0] = ty;
                    for k in 0..3 {
                        let n = &nodes.v[self.verts[curr].nodes[k]];
                        x[k + 1] = n.x;
                        y[k + 1] = n.y;
                    }
                    break;
                }
            }
            for i in 0..3 {
                if self.verts[curr].nodes[i] == a {
                    found = true;
                    holder = i;
                    if self.verts[curr].nodes[(i + 1) % 3] == b {
                        return Some(curr);
                    }
                    if self.verts[curr].nodes[(i + 3 - 1) % 3] == b {
                        return Some(curr);
                    }
                    break;
                }
            }
        }
        // curr owns a; go around a through the neighbour links.
        let startv = curr;
        loop {
            curr = self.verts[curr].vertices[holder]?;
            let mut i = 0;
            while i < 3 {
                if self.verts[curr].nodes[i] == a {
                    holder = i;
                    break;
                }
                i += 1;
            }
            if self.verts[curr].nodes[(i + 1) % 3] == b {
                return Some(curr);
            }
            if self.verts[curr].nodes[(i + 3 - 1) % 3] == b {
                return Some(curr);
            }
            if curr == startv {
                return None;
            }
        }
    }

    /// `Vertex::vertexWith`.
    pub fn vertex_with(&self, nodes: &Nodes, start: usize, a: NodeId, b: NodeId) -> Option<usize> {
        let mut found = false;
        let mut curr = start;
        let mut prev: Option<usize> = None;
        let tx = (nodes.v[a].x + nodes.v[b].x) / 2.0;
        let ty = (nodes.v[a].y + nodes.v[b].y) / 2.0;
        let s = [2usize, 3, 1];
        let e = [1usize, 2, 3];
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        x[0] = tx;
        y[0] = ty;
        for k in 0..3 {
            let n = &nodes.v[self.verts[start].nodes[k]];
            x[k + 1] = n.x;
            y[k + 1] = n.y;
        }
        for i in 0..3 {
            if self.verts[start].nodes[i] == a {
                if self.verts[start].nodes[(i + 1) % 3] == b {
                    return Some(curr);
                }
                if self.verts[start].nodes[(i + 3 - 1) % 3] == b {
                    return self.verts[curr].vertices[(i + 3 - 1) % 3];
                }
            }
        }
        while !found {
            let mut i = 0;
            while i < 3 {
                let nb = self.verts[curr].vertices[i];
                if let Some(n) = nb.filter(|_| {
                    nb != prev && orientation(tx, ty, x[s[i]], y[s[i]], x[e[i]], y[e[i]])
                }) {
                    prev = Some(curr);
                    curr = n;
                    x[0] = tx;
                    y[0] = ty;
                    for k in 0..3 {
                        let n = &nodes.v[self.verts[curr].nodes[k]];
                        x[k + 1] = n.x;
                        y[k + 1] = n.y;
                    }
                    break;
                }
                i += 1;
            }
            for j in 0..3 {
                if self.verts[curr].nodes[j] == a && self.verts[curr].nodes[(j + 1) % 3] == b {
                    found = true;
                    break;
                } else if self.verts[curr].nodes[j] == a
                    && self.verts[curr].nodes[(j + 3 - 1) % 3] == b
                {
                    curr = self.verts[curr].vertices[(j + 3 - 1) % 3]?;
                    found = true;
                    break;
                }
            }
            if !found && i == 3 {
                return None;
            }
        }
        Some(curr)
    }

    /// `Vertex::labelBodies`, as a flood fill.
    fn label_bodies(&mut self, nodes: &Nodes, start: usize, bd: i32) {
        let mut stack = vec![start];
        while let Some(v) = stack.pop() {
            if self.verts[v].body == bd {
                continue;
            }
            self.verts[v].body = bd;
            for i in 0..3 {
                let t1 = nodes.v[self.verts[v].nodes[i]].tag;
                let t2 = nodes.v[self.verts[v].nodes[(i + 1) % 3]].tag;
                if !self.links.contains(&(t1.min(t2), t1.max(t2))) {
                    if let Some(n) = self.verts[v].vertices[i] {
                        stack.push(n);
                    }
                }
            }
        }
    }

    /// `Connect::addSite`.
    pub fn add_site(
        &mut self,
        nodes: &Nodes,
        vx: usize,
        nd: NodeId,
        external_ok: bool,
        split_one: bool,
        force: bool,
    ) -> bool {
        self.hull.clear();
        let mut queue: VecDeque<Triple> = VecDeque::new();
        for i in (0..3).rev() {
            queue.push_back(Triple {
                from: self.verts[vx].nodes[i],
                to: self.verts[vx].nodes[(i + 1) % 3],
                vtx: self.verts[vx].vertices[i],
            });
        }
        self.verts[vx].deleted = true;
        self.deleted.push(vx);

        while let Some(trip) = queue.pop_back() {
            let Some(vtx) = trip.vtx else {
                self.hull.push(trip);
                continue;
            };
            if self.verts[vtx].deleted {
                // We must have visited it in this loop, so the reversed
                // triple is on the queue; drop both.
                let tmp = *queue.back().expect("triple queue");
                if tmp.from == trip.to && tmp.to == trip.from {
                    queue.pop_back();
                } else {
                    queue.pop_front();
                }
            } else {
                let retval =
                    self.is_destroyed_by(nodes, vtx, nodes.v[nd].x, nodes.v[nd].y, external_ok);
                if retval == 1 {
                    self.verts[vtx].deleted = true;
                    self.deleted.push(vtx);
                    let mut i = 0;
                    while i < 3 {
                        if self.verts[vtx].nodes[i] == trip.from {
                            break;
                        }
                        i += 1;
                    }
                    let vn = self.verts[vtx].nodes;
                    let vv = self.verts[vtx].vertices;
                    queue.push_back(Triple {
                        from: vn[(i + 1) % 3],
                        to: vn[(i + 2) % 3],
                        vtx: vv[(i + 1) % 3],
                    });
                    queue.push_back(Triple {
                        from: vn[i % 3],
                        to: vn[(i + 1) % 3],
                        vtx: vv[i % 3],
                    });
                } else if retval == 0 || force {
                    self.hull.push(trip);
                } else {
                    for &d in &self.deleted {
                        self.verts[d].deleted = false;
                    }
                    self.deleted.clear();
                    return false;
                }
            }
        }
        if !external_ok && !split_one && self.deleted.len() == 1 {
            return true;
        }

        let count = self.hull.len();
        let nv = self.get_vertices(count);
        self.new_vertices.extend(nv);
        for i in 0..count {
            let h = self.hull[i];
            let v = self.new_vertices[i];
            self.reset(nodes, v, nd, h.from, h.to);
        }
        for i in 0..count {
            let h = self.hull[i];
            let v = self.new_vertices[i];
            self.verts[v].vertices = [
                Some(self.new_vertices[(i + count - 1) % count]),
                h.vtx,
                Some(self.new_vertices[(i + 1) % count]),
            ];
            if let Some(t) = h.vtx {
                if self.verts[t].nodes[0] == h.to {
                    self.verts[t].vertices[0] = Some(v);
                } else if self.verts[t].nodes[1] == h.to {
                    self.verts[t].vertices[1] = Some(v);
                } else if self.verts[t].nodes[2] == h.to {
                    self.verts[t].vertices[2] = Some(v);
                }
            }
        }
        self.root = *self.new_vertices.last().expect("new vertices");
        true
    }

    /// `Connect::recycle`.
    pub fn recycle(&mut self) {
        self.vertex_store.extend(self.deleted.iter().copied());
        self.deleted.clear();
        self.new_vertices.clear();
        self.hull.clear();
    }

    /// `Element::isBoundaryConnector`.
    pub fn is_boundary_connector(&self, nodes: &Nodes, v: usize) -> bool {
        let n = &self.verts[v].nodes;
        let mut count = 0;
        for i in 0..2 {
            for j in (i + 1)..3 {
                let (a, b) = (&nodes.v[n[i]], &nodes.v[n[j]]);
                if a.boundarynode
                    && b.boundarynode
                    && self.links.contains(&(a.tag.min(b.tag), a.tag.max(b.tag)))
                {
                    count += 1;
                }
            }
        }
        count > 1
    }

    /// `Connect::triangulate`.
    pub fn triangulate(
        &mut self,
        nodes: &mut Nodes,
        rng: &mut GlibcRand,
        fixed_geometry_nodes: &[NodeId],
        tag: i32,
    ) -> Result<(), String> {
        self.make_world(nodes);

        for &nd in fixed_geometry_nodes {
            let vtx = self.first_acceptable_found(nodes, self.root, nodes.v[nd].x, nodes.v[nd].y);
            self.add_site(nodes, vtx, nd, true, true, false);
            nodes.v[nd].fixed = true;
            self.recycle();
        }

        let len = self.border.len();
        let mut indirect: Vec<usize> = (0..len).collect();
        rng.random_shuffle(&mut indirect);
        for &ind in &indirect {
            let b = self.border[ind];
            let vtx = self.first_acceptable_found(nodes, self.root, nodes.v[b].x, nodes.v[b].y);
            self.add_site(nodes, vtx, b, true, true, false);
            nodes.v[b].fixed = true;
            self.recycle();
        }

        // Verify borders (no recovery).
        let mut verify = self.root;
        for i in 0..len {
            let (a, b) = (self.border[i], self.border[(i + 1) % len]);
            let (t1, t2) = (nodes.v[a].tag, nodes.v[b].tag);
            if self.links.contains(&(t1.min(t2), t1.max(t2))) {
                match self.vertex_owning(nodes, verify, a, b) {
                    Some(v) => verify = v,
                    None => {
                        return Err("Initial triangulation is incomplete! Please readjust your mesh parameters.".into())
                    }
                }
            }
        }

        let vtx = self
            .vertex_with(nodes, self.root, self.border[0], self.border[1])
            .ok_or("no triangle on the first border edge")?;
        self.label_bodies(nodes, vtx, tag);
        Ok(())
    }

    /// `Connect::removeBoundaryConnectors`.
    pub fn remove_boundary_connectors(&mut self, nodes: &mut Nodes, actives: &mut Pq) {
        let vs = self.vs;
        let all = self.all_vertices.clone();
        for &v in &all {
            let vert = &self.verts[v];
            if !(vert.deleted || vert.is_external()) && self.is_boundary_connector(nodes, v) {
                self.verts[v].right_size(1.0);
                actives.insert(&mut self.verts, vs, v);
            }
        }

        while actives.size() > 0 {
            let vtx = actives.first(&mut self.verts, vs).expect("active vertex");

            if !self.is_boundary_connector(nodes, vtx) {
                continue;
            }

            let (mut x, mut y) = (0.0, 0.0);
            for i in 0..3 {
                let n = self.verts[vtx].nodes;
                let (t0, t1) = (nodes.v[n[i]].tag, nodes.v[n[(i + 1) % 3]].tag);
                if !self.links.contains(&(t0.min(t1), t0.max(t1))) {
                    let k = &nodes.v[n[(i + 2) % 3]];
                    let dx0 = nodes.v[n[i]].x - k.x;
                    let dy0 = nodes.v[n[i]].y - k.y;
                    let dx1 = nodes.v[n[(i + 1) % 3]].x - k.x;
                    let dy1 = nodes.v[n[(i + 1) % 3]].y - k.y;

                    let x2y2 = dx0 * dx1 - dy0 * dy1;
                    let xy = dx0 * dy1 + dx1 * dy0;

                    y = (0.5 * ((x2y2 * x2y2 + xy * xy).sqrt() - x2y2)).sqrt();
                    x = 0.5 * xy / y;

                    if dx0 * x + dy0 * y < 0.0 {
                        x = -x;
                        y = -y;
                    }

                    x += k.x;
                    y += k.y;
                    break;
                }
            }

            let nd = nodes.mesh_node(x, y);

            if self.add_site(nodes, vtx, nd, false, true, false) {
                let del = self.deleted.clone();
                actives.remove_relevants(&mut self.verts, vs, &del);

                let nv = self.new_vertices.clone();
                for &v in &nv {
                    self.verts[v].body = 1;
                    if self.is_boundary_connector(nodes, v) {
                        self.verts[v].right_size(1.0);
                        actives.insert(&mut self.verts, vs, v);
                    }
                }
            }

            self.recycle();
        }
    }

    // ---- BGTriangleMesh / BGVertex ----

    /// `BGTriangleMesh::initialize`.
    pub fn bg_initialize(
        &mut self,
        nodes: &mut Nodes,
        rng: &mut GlibcRand,
        geometry_nodes: &[NodeId],
    ) {
        self.initialized = true;
        self.border = geometry_nodes.to_vec();
        self.make_world(nodes);
        let len = self.border.len();
        let mut indirect: Vec<usize> = (0..len).collect();
        rng.random_shuffle(&mut indirect);
        for &ind in &indirect {
            let b = self.border[ind];
            let vtx = self.first_acceptable_found(nodes, self.root, nodes.v[b].x, nodes.v[b].y);
            self.add_site(nodes, vtx, b, true, true, false);
            self.recycle();
        }
        let all = self.all_vertices.clone();
        for &v in &all {
            if !self.verts[v].deleted {
                self.init_interpolation(nodes, v);
            }
        }
    }

    /// `BGVertex::initInterpolation`.
    fn init_interpolation(&mut self, nodes: &Nodes, v: usize) {
        let n = self.verts[v].nodes;
        let mut a = [0.0; 9];
        let mut b = [0.0; 3];
        for i in 0..3 {
            let i1 = (i + 1) % 3;
            let i2 = (i + 2) % 3;
            a[i] = nodes.v[n[i1]].x * nodes.v[n[i2]].y - nodes.v[n[i2]].x * nodes.v[n[i1]].y;
            a[3 + i] = nodes.v[n[i1]].y - nodes.v[n[i2]].y;
            a[6 + i] = nodes.v[n[i2]].x - nodes.v[n[i1]].x;
        }
        let det = a[0] + a[1] + a[2];
        let mut cnt = 0;
        let mut avg = 0.0;
        for i in 0..3 {
            if nodes.v[n[i]].tag < 0 {
                continue;
            }
            b[i] = map(nodes.v[n[i]].delta);
            avg += b[i];
            cnt += 1;
        }
        avg /= cnt as f64;
        for i in 0..3 {
            if nodes.v[n[i]].tag < 0 {
                b[i] = avg;
            }
        }
        for i in 0..3 {
            self.verts[v].coeff[i] =
                (a[3 * i] * b[0] + a[3 * i + 1] * b[1] + a[3 * i + 2] * b[2]) / det;
        }
    }

    /// `BGTriangleMesh::interpolate` (`BGVertex::interpolate` from the root).
    pub fn interpolate(&self, nodes: &Nodes, tx: f64, ty: f64) -> f64 {
        let s = [2usize, 3, 1];
        let e = [1usize, 2, 3];
        let mut curr = self.root;
        let mut prev: Option<usize> = None;
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        x[0] = tx;
        y[0] = ty;
        for k in 0..3 {
            let n = &nodes.v[self.verts[self.root].nodes[k]];
            x[k + 1] = n.x;
            y[k + 1] = n.y;
        }
        loop {
            let test = curr;
            for i in 0..3 {
                let nb = self.verts[curr].vertices[i];
                if let Some(n) = nb.filter(|_| {
                    nb != prev && orientation(tx, ty, x[s[i]], y[s[i]], x[e[i]], y[e[i]])
                }) {
                    prev = Some(curr);
                    curr = n;
                    x[0] = tx;
                    y[0] = ty;
                    for k in 0..3 {
                        let n = &nodes.v[self.verts[curr].nodes[k]];
                        x[k + 1] = n.x;
                        y[k + 1] = n.y;
                    }
                    break;
                }
            }
            if test == curr {
                break;
            }
        }
        let c = &self.verts[curr].coeff;
        unmap(c[0] + c[1] * tx + c[2] * ty)
    }
}
