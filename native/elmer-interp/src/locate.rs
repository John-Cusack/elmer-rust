//! Point location in the old mesh (the search half of
//! `InterpolateMeshToMeshQ`).

use std::cell::{Cell, OnceCell};

use elmer_core::{ElemId, MeshView};

use crate::geom::{point_in_element, PointTol};
use crate::octree::{padded_bounding_box, QuadrantTree};
use crate::shape::Shape;
use crate::{InterpError, InterpOptions};

/// A point found inside an old-mesh bulk element.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Located {
    pub element: ElemId,
    /// Local coordinates in `element`.
    pub uvw: [f64; 3],
}

/// Selects the bulk elements the tree search may consider, by body id.
/// The Fortran passes a `MaskName` and keeps bodies whose `Body Force`
/// section defines it; the caller evaluates that and hands the answer here.
pub type BodyMask<'a> = &'a dyn Fn(i32) -> bool;

/// The old mesh plus its (lazily built, cached) quadrant tree.
pub struct Locator<'a> {
    mesh: &'a dyn MeshView,
    opts: InterpOptions,
    tree: OnceCell<Option<QuadrantTree>>,
    qtree_fails: Cell<usize>,
}

impl<'a> Locator<'a> {
    /// Fails if a bulk element has a type with no nodal basis table.
    pub fn new(mesh: &'a dyn MeshView, opts: InterpOptions) -> Result<Self, InterpError> {
        for e in 0..mesh.n_bulk() {
            let el = mesh.element(ElemId(e as u32));
            if Shape::get(el.etype.0).is_none() {
                return Err(InterpError::UnsupportedElement(el.etype.0));
            }
        }
        Ok(Locator {
            mesh,
            opts,
            tree: OnceCell::new(),
            qtree_fails: Cell::new(0),
        })
    }

    pub fn options(&self) -> &InterpOptions {
        &self.opts
    }

    pub fn mesh(&self) -> &'a dyn MeshView {
        self.mesh
    }

    /// Points the tree missed but the exhaustive scan found, so far.
    pub fn qtree_fails(&self) -> usize {
        self.qtree_fails.get()
    }

    /// `max(CoordinateSystemDimension(), MeshDim)` of the old mesh.
    fn dim(&self) -> usize {
        self.opts
            .coordinate_system_dim
            .unwrap_or(0)
            .max(self.mesh.dim())
    }

    fn coord_dim(&self) -> usize {
        self.opts
            .coordinate_system_dim
            .unwrap_or_else(|| self.mesh.dim())
    }

    fn tree(&self) -> Option<&QuadrantTree> {
        self.tree
            .get_or_init(|| {
                if !self.opts.use_quadrant_tree {
                    return None;
                }
                padded_bounding_box(self.mesh)
                    .map(|bb| QuadrantTree::build(self.mesh, bb, self.dim()))
            })
            .as_ref()
    }

    /// Map a target-mesh node to the search point (passive coordinate and
    /// cylindrical projection).
    pub fn search_point(&self, mut p: [f64; 3]) -> [f64; 3] {
        if self.opts.passive_coordinate != 0 {
            p[self.opts.passive_coordinate as usize - 1] = 0.0;
        }
        if self.opts.cylindric {
            p = [(p[0] * p[0] + p[1] * p[1]).sqrt(), p[2], 0.0];
        }
        p
    }

    fn element_nodes(&self, e: u32, buf: &mut Vec<[f64; 3]>) -> u16 {
        let el = self.mesh.element(ElemId(e));
        buf.clear();
        for &n in el.nodes {
            let mut c = self.mesh.coord(n);
            if self.opts.passive_coordinate != 0 {
                c[self.opts.passive_coordinate as usize - 1] = 0.0;
            }
            buf.push(c);
        }
        el.etype.0
    }

    /// Locate `point` (already passed through [`Locator::search_point`]).
    pub fn locate(&self, point: [f64; 3], mask: Option<BodyMask<'_>>) -> Option<Located> {
        let coord_dim = self.coord_dim();
        let piola = self.opts.piola_transform;
        let mut buf: Vec<[f64; 3]> = Vec::new();
        let try_tree = self.tree().is_some();

        if let Some(tree) = self.tree() {
            if let Some(leaf) = tree.leaf_elements(point) {
                let mut eps1 = self.opts.global_epsilon;
                let mut eps2 = self.opts.local_epsilon;
                for _ in 0..self.opts.max_iterations {
                    for &e in leaf {
                        if let Some(m) = mask {
                            if !m(self.mesh.element(ElemId(e)).tag) {
                                continue;
                            }
                        }
                        let code = self.element_nodes(e, &mut buf);
                        let shape = Shape::get(code)?;
                        let tol = PointTol {
                            global: Some(eps1),
                            local: Some(eps2),
                            numeric: Some(self.opts.numeric_epsilon),
                        };
                        let r = point_in_element(shape, &buf, point, tol, piola, coord_dim);
                        if r.inside {
                            return Some(Located {
                                element: ElemId(e),
                                uvw: r.uvw,
                            });
                        }
                    }
                    eps1 *= 10.0;
                    eps2 *= 10.0;
                    if eps1 > 1.0 {
                        break;
                    }
                }
            }
        }

        if !try_tree || self.opts.try_linear_search {
            for e in 0..self.mesh.n_bulk() as u32 {
                let code = self.element_nodes(e, &mut buf);
                let shape = Shape::get(code)?;
                let r = point_in_element(shape, &buf, point, PointTol::default(), false, coord_dim);
                if r.inside {
                    if try_tree {
                        self.qtree_fails.set(self.qtree_fails.get() + 1);
                    }
                    return Some(Located {
                        element: ElemId(e),
                        uvw: r.uvw,
                    });
                }
            }
        }
        None
    }
}
