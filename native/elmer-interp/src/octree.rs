//! Quadrant (2-D) / octant (3-D) tree over the bulk elements of a mesh.
//!
//! Ports `BuildQuadrantTree` and `FindLeafElements` (`Interpolation.F90`).
//! Elements are binned by their padded bounding box, so an element can sit in
//! several leaves; a query descends through the first child whose closed
//! bounding box holds the point and returns that leaf's element list.

use elmer_core::{ElemId, MeshView};

/// Bounding-box padding applied to elements when binning them (`eps2` in
/// `PutElementsInChildQuadrants`), as a fraction of the child's largest side.
const ELEMENT_PAD: f64 = 2.5e-2;

#[derive(Clone, Debug)]
struct Quadrant {
    bbox: [f64; 6],
    /// Leaf contents, indices into the bulk elements. Empty for interior nodes.
    elements: Vec<u32>,
    /// Child quadrants (arena indices), empty for leaves.
    children: Vec<usize>,
}

/// Spatial subdivision of one mesh's bulk elements.
#[derive(Clone, Debug)]
pub struct QuadrantTree {
    nodes: Vec<Quadrant>,
    dim: usize,
}

/// Bounding box of all mesh nodes grown by 10 % of its largest side
/// (`InterpolateMeshToMeshQ`). `None` for a mesh without nodes.
pub fn padded_bounding_box(mesh: &dyn MeshView) -> Option<[f64; 6]> {
    if mesh.n_nodes() == 0 {
        return None;
    }
    let mut bb = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for i in 0..mesh.n_nodes() {
        let c = mesh.coord(i as u32);
        for d in 0..3 {
            bb[d] = bb[d].min(c[d]);
            bb[d + 3] = bb[d + 3].max(c[d]);
        }
    }
    let eps = 0.1
        * (0..3)
            .map(|d| bb[d + 3] - bb[d])
            .fold(f64::NEG_INFINITY, f64::max);
    for d in 0..3 {
        bb[d] -= eps;
        bb[d + 3] += eps;
    }
    Some(bb)
}

impl QuadrantTree {
    /// `BuildQuadrantTree`; `dim` is `MAX(MeshDim, CoordinateSystemDimension())`.
    pub fn build(mesh: &dyn MeshView, bbox: [f64; 6], dim: usize) -> Self {
        let max_leaf = if dim == 3 { 16 } else { 8 };
        let mut root_box = [bbox[0], 0.0, 0.0, bbox[3], 0.0, 0.0];
        if dim >= 2 {
            root_box[1] = bbox[1];
            root_box[4] = bbox[4];
        }
        if dim == 3 {
            root_box[2] = bbox[2];
            root_box[5] = bbox[5];
        }
        let n_bulk = mesh.n_bulk();
        let mut tree = QuadrantTree {
            nodes: vec![Quadrant {
                bbox: root_box,
                elements: (0..n_bulk as u32).collect(),
                children: Vec::new(),
            }],
            dim,
        };
        let ebox = |e: u32| -> [f64; 6] {
            let el = mesh.element(ElemId(e));
            let mut b = [
                f64::INFINITY,
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ];
            for &nd in el.nodes {
                let c = mesh.coord(nd);
                for d in 0..3 {
                    b[d] = b[d].min(c[d]);
                    b[d + 3] = b[d + 3].max(c[d]);
                }
            }
            b
        };
        let boxes: Vec<[f64; 6]> = (0..n_bulk as u32).map(ebox).collect();
        tree.split(0, 0, &boxes, max_leaf);
        tree
    }

    /// `CreateChildQuadrants`; `generation` is the Fortran's recursion counter.
    fn split(&mut self, mother: usize, generation: usize, boxes: &[[f64; 6]], max_leaf: usize) {
        let dim = self.dim;
        let nchild = 1usize << dim;
        let mb = self.nodes[mother].bbox;
        let size = (mb[3] - mb[0]).max(mb[4] - mb[1]).max(mb[5] - mb[2]);
        let mid = [
            (mb[0] + mb[3]) / 2.0,
            (mb[1] + mb[4]) / 2.0,
            (mb[2] + mb[5]) / 2.0,
        ];
        // Child i (0-based) holds bit 0 -> x half, bit 1 -> y half, bit 2 -> z half.
        let mut child_boxes = Vec::with_capacity(nchild);
        for i in 0..nchild {
            let mut b = [mb[0], mb[1], mb[2], mid[0], mid[1], mid[2]];
            if i & 1 != 0 {
                b[0] = mid[0];
                b[3] = mb[3];
            }
            if i & 2 != 0 {
                b[1] = mid[1];
                b[4] = mb[4];
            }
            if i & 4 != 0 {
                b[2] = mid[2];
                b[5] = mb[5];
            }
            child_boxes.push(b);
        }
        // PutElementsInChildQuadrants
        let mother_elems = std::mem::take(&mut self.nodes[mother].elements);
        let mut lists: Vec<Vec<u32>> = vec![Vec::new(); nchild];
        let mut min_size = vec![1.0e20f64; nchild];
        for &e in &mother_elems {
            let eb = boxes[e as usize];
            let esize = (eb[3] - eb[0]).max(eb[4] - eb[1]).max(eb[5] - eb[2]);
            for (i, cb) in child_boxes.iter().enumerate() {
                let pad = ELEMENT_PAD
                    * (cb[3] - cb[0])
                        .max(cb[4] - cb[1])
                        .max(cb[5] - cb[2])
                        .max(0.0);
                let out = eb[3] < cb[0] - pad
                    || eb[0] > cb[3] + pad
                    || eb[4] < cb[1] - pad
                    || eb[1] > cb[4] + pad
                    || eb[5] < cb[2] - pad
                    || eb[2] > cb[5] + pad;
                if !out {
                    lists[i].push(e);
                    min_size[i] = min_size[i].min(esize);
                }
            }
        }
        let first_child = self.nodes.len();
        for (i, cb) in child_boxes.iter().enumerate() {
            self.nodes.push(Quadrant {
                bbox: *cb,
                elements: std::mem::take(&mut lists[i]),
                children: Vec::new(),
            });
        }
        self.nodes[mother].children = (first_child..first_child + nchild).collect();
        let child_size = size / 2.0;
        for (i, &child_min) in min_size.iter().enumerate() {
            let id = first_child + i;
            if self.nodes[id].elements.len() > max_leaf && child_size > child_min && generation <= 8
            {
                self.split(id, generation + 1, boxes, max_leaf);
            }
        }
    }

    /// `FindLeafElements`: the leaf holding `point`, or `None` when the point
    /// lies outside every child at some level.
    pub fn leaf_elements(&self, point: [f64; 3]) -> Option<&[u32]> {
        let nchild = 1usize << self.dim;
        let mut cur = 0usize;
        loop {
            let node = &self.nodes[cur];
            if node.children.is_empty() {
                return Some(&node.elements);
            }
            let next = node.children[..nchild].iter().copied().find(|&c| {
                let b = &self.nodes[c].bbox;
                (0..3).all(|d| point[d] >= b[d] && point[d] <= b[d + 3])
            })?;
            cur = next;
        }
    }

    /// Number of quadrants including the root.
    pub fn n_quadrants(&self) -> usize {
        self.nodes.len()
    }
}
