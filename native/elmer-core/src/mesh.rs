//! Consumer view of a mesh; `elmer-mesh` implements [`MeshView`].

/// Elmer element type code (elements.def): 101, 202, 303, 404, 504, 808, 306,
/// 408, 510, 820, 827, ... The hundreds digit is the family and the rest the
/// node count.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ElementType(pub u16);

impl ElementType {
    pub fn family(self) -> u16 {
        self.0 / 100
    }

    pub fn n_nodes(self) -> usize {
        (self.0 % 100) as usize
    }
}

/// Bulk elements are `0..n_bulk`, then boundary elements (Elmer
/// `Mesh % Elements` order).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ElemId(pub u32);

#[derive(Copy, Clone, Debug)]
pub struct ElementRef<'a> {
    pub etype: ElementType,
    /// 0-based node indices, Elmer local order.
    pub nodes: &'a [u32],
    /// Body id (bulk) or boundary tag (boundary), as in the mesh files.
    pub tag: i32,
    /// Boundary elements only.
    pub parents: [Option<ElemId>; 2],
}

pub trait MeshView: Sync {
    fn dim(&self) -> usize;
    fn n_nodes(&self) -> usize;
    fn coord(&self, node: u32) -> [f64; 3];
    fn n_bulk(&self) -> usize;
    fn n_boundary(&self) -> usize;
    fn element(&self, e: ElemId) -> ElementRef<'_>;
    /// Global edge ids, for edge and p elements.
    fn element_edges(&self, _e: ElemId) -> Option<&[u32]> {
        None
    }
    fn element_faces(&self, _e: ElemId) -> Option<&[u32]> {
        None
    }
}
