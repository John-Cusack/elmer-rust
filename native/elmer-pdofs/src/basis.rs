//! The p basis and integration rule of one element, as `ElementInfo`
//! (`EvalPElementBasis`, ElemInfo.F90:750) and `GaussPoints`
//! (Integration.F90:2882) build them for an active p element.
//!
//! The basis is `elmer-fe-hp`'s [`PElementBasis`]: nodal functions, then
//! `p - 1` functions on every edge, oriented from the lower to the higher
//! node index, then the bubbles of the effective bubble degree. The
//! degree is the solver's `Def_Dofs(family, body, 6)` of the element's own
//! body, so on a varying-p mesh the shared edge carries as many slots as
//! the higher-degree neighbour while each element uses only its own
//! `p - 1` of them (the first ones).
//!
//! Values agree with Elmer's scalar `ElementInfo` path to rounding, not
//! bitwise (`H1Basis` and `PElementBase` evaluate the same polynomials in a
//! different order).

use elmer_core::ElementType;
use elmer_fe::quadrature::{
    gauss_points, gauss_points_1d, gauss_points_brick, gauss_points_quad, GaussRule,
};
use elmer_fe_hp::maps::{
    bubble_dofs, edge_direction, face_direction, tetra_edge_direction, tetra_face_direction,
};
use elmer_fe_hp::{PElementBasis, PElementSpec, MAX_EDGES, MAX_FACES};

use crate::defs::Kind;
use crate::dofs::DofLayout;
use crate::pmesh::effective_bubble_p;
use crate::quadrature::{
    gauss_points_p_pyramid, gauss_points_p_triangle, gauss_points_quad_p, gauss_points_triangle_p,
    p_tetra_rule, p_triangle_rule, p_wedge_rule, triangle_rule_tabulated,
};
use crate::Error;

impl DofLayout<'_> {
    /// The body `EvalPElementBasis` reads the degree from: the element's,
    /// else (boundary elements) the local parent's, the left parent's or
    /// the right parent's, else 1.
    pub fn basis_body(&self, e: usize) -> usize {
        let pm = self.mesh;
        let el = &pm.elements[e];
        let body_of = |p: Option<u32>| p.map_or(0, |p| pm.elements[p as usize].body);
        let mut id = el.body;
        if id == 0 && e >= pm.n_bulk {
            id = body_of(el.pdefs.as_ref().and_then(|pd| pd.local_parent));
            if id == 0 {
                id = body_of(el.parents[0]);
            }
            if id == 0 {
                id = body_of(el.parents[1]);
            }
        }
        if id == 0 {
            1
        } else {
            id as usize
        }
    }

    /// The `elmer-fe-hp` description of element `e`'s p basis.
    pub fn element_spec(&self, e: usize) -> Result<PElementSpec, Error> {
        if !self.is_active_p(e) {
            return Err(Error::Unsupported(format!(
                "element {e} is not an active p element"
            )));
        }
        let pm = self.mesh;
        let el = &pm.elements[e];
        let pd = el.pdefs.as_ref().expect("active p element");
        let id = self.basis_body(e);
        let fam = el.family();
        let p = self.defs.get(fam, id, Kind::P);
        let nb = self.defs.get(fam, id, Kind::Bubble);
        let sd = pd.serendipity;
        let unsupported = |what: &str| {
            Err(Error::Unsupported(format!(
                "{what} (element {e}, type {})",
                el.code
            )))
        };
        let mut spec = PElementSpec::uniform(el.code, 1, None, None).map_err(Error::Unsupported)?;
        spec.serendipity = sd;
        spec.edge_degree = [0; MAX_EDGES];
        spec.face_degree = [0; MAX_FACES];
        spec.bubble_p = 0;
        let gind: Vec<i64> = el.nodes.iter().map(|&n| i64::from(n)).collect();
        let bdofs = bubble_dofs(el.code, p, sd).max(nb.max(0) as usize);
        match el.code {
            202 => {
                if bdofs > 0 {
                    // LineBubblePBasis(i+1) for i = 1..BDOFs.
                    spec.bubble_p = bdofs as i32 + 1;
                    spec.invert_bubble = pd.is_edge && gind[0] > gind[1];
                }
            }
            303 | 404 => {
                let edofs = pm.edge_dofs_of(e, p);
                let ne = el.n_edges();
                if !el.edges.is_empty() && edofs > 0 {
                    spec.edge_degree[..ne].fill(edofs as i32 + 1);
                }
                spec.edge_dir = edge_direction(el.code, ne, &gind).expect("2D edge map");
                if bdofs > 0 {
                    spec.bubble_p = effective_bubble_p(el.code, p, bdofs, sd);
                    // A face of a 3D mesh orients its bubbles
                    // (getTriangleFaceDirection / getSquareFaceDirection).
                    if pd.is_edge {
                        let dir = face_direction(el.code, 1, &gind).expect("2D face map");
                        spec.bubble_face_dir = Some(dir[0]);
                    }
                }
            }
            504 => {
                // Edges and faces on the TetraType tables (ConvertToACTetra
                // ordered the nodes), effective-p bubbles.
                let tt = pd.tetra_type.max(1);
                let edofs = pm.edge_dofs_of(e, p);
                if !el.edges.is_empty() && edofs > 0 {
                    spec.edge_degree[..6].fill(edofs as i32 + 1);
                }
                spec.edge_dir = tetra_edge_direction(tt);
                if !el.faces.is_empty() {
                    spec.face_degree[..4].fill(p);
                }
                spec.face_dir = tetra_face_direction(tt);
                if bdofs > 0 {
                    spec.bubble_p = effective_bubble_p(504, p, bdofs, sd);
                }
            }
            605 => {
                if sd {
                    return Err(Error::BadDefinition(
                        "ElementInfo: p-Pyramid not implemented for serendipity scheme, please use the full scheme instead."
                            .into(),
                    ));
                }
                let edofs = pm.edge_dofs_of(e, p);
                if !el.edges.is_empty() && edofs > 0 {
                    spec.edge_degree[..8].fill(edofs as i32 + 1);
                }
                spec.edge_dir = edge_direction(605, 8, &gind).expect("pyramid edge map");
                if !el.faces.is_empty() {
                    spec.face_degree[..5].fill(p);
                }
                spec.face_dir = face_direction(605, 5, &gind).expect("pyramid face map");
                if bdofs > 0 {
                    spec.bubble_p = effective_bubble_p(605, p, bdofs, sd);
                }
            }
            706 => {
                // As bricks: `p - 1` functions per edge, degree-`p` faces
                // (triangles 1-2, squares 3-5, in the scalar ElementInfo
                // order), effective-p bubbles.
                let edofs = pm.edge_dofs_of(e, p);
                if !el.edges.is_empty() && edofs > 0 {
                    spec.edge_degree[..9].fill(edofs as i32 + 1);
                }
                spec.edge_dir = edge_direction(706, 9, &gind).expect("wedge edge map");
                if !el.faces.is_empty() {
                    spec.face_degree[..5].fill(p);
                }
                spec.face_dir = face_direction(706, 5, &gind).expect("wedge face map");
                spec.scalar_face_order = true;
                if bdofs > 0 {
                    spec.bubble_p = effective_bubble_p(706, p, bdofs, sd);
                }
            }
            808 => {
                // Edges: `p - 1` functions each; faces: the degree-`p` face
                // functions, oriented by `getSquareFaceDirection`.
                let edofs = pm.edge_dofs_of(e, p);
                if !el.edges.is_empty() && edofs > 0 {
                    spec.edge_degree[..12].fill(edofs as i32 + 1);
                }
                spec.edge_dir = edge_direction(808, 12, &gind).expect("brick edge map");
                if !el.faces.is_empty() {
                    spec.face_degree[..6].fill(p);
                }
                spec.face_dir = face_direction(808, 6, &gind).expect("brick face map");
                if bdofs > 0 {
                    spec.bubble_p = effective_bubble_p(808, p, bdofs, sd);
                }
            }
            _ => return unsupported("p basis"),
        }
        Ok(spec)
    }

    /// [`DofLayout::element_spec`] as an `ElementBasis`.
    pub fn element_basis(&self, e: usize) -> Result<PElementBasis, Error> {
        Ok(PElementBasis::new(self.element_spec(e)?))
    }

    /// `GaussPoints(Element)` for an active p element: `PDefs %
    /// GaussPoints` points (the type default when 0) of the line rule, the
    /// quad rule (economic for serendipity elements of degree 2 to 7
    /// without low-degree bubbles) or the p reference triangle rule.
    pub fn gauss_rule(&self, e: usize) -> Result<GaussRule, Error> {
        self.p_gauss_rule(e, None)
    }

    /// `GaussPoints(Element, np)` for an active p element: the same
    /// dispatch on `np` points, except that a triangle takes a tabulated
    /// rule only when `np` is one (`TriangleRuleTabulated`).
    pub fn gauss_rule_np(&self, e: usize, np: usize) -> Result<GaussRule, Error> {
        self.p_gauss_rule(e, Some(np))
    }

    fn p_gauss_rule(&self, e: usize, np: Option<usize>) -> Result<GaussRule, Error> {
        let el = &self.mesh.elements[e];
        let etype = ElementType(el.code);
        let fe = |r: Result<GaussRule, elmer_fe::FeError>| {
            r.map_err(|err| Error::Unsupported(err.to_string()))
        };
        let Some(pd) = el.pdefs.as_ref().filter(|_| self.is_active_p(e)) else {
            return fe(gauss_points(etype, np));
        };
        // `IF (n == 0) n = elmt % GaussPoints`, still on the p reference element.
        let n = match np.unwrap_or(pd.gauss_points) {
            0 => {
                elmer_fe::element_type(el.code)
                    .map_err(|err| Error::Unsupported(err.to_string()))?
                    .gauss_points
            }
            n => n,
        };
        match el.family() {
            2 => fe(gauss_points_1d(n)),
            4 => {
                let economic =
                    pd.serendipity && pd.p > 1 && pd.p <= 7 && !(el.bdofs > 0 && pd.p < 4);
                if economic && n == 12 {
                    fe(gauss_points_quad(16))
                } else {
                    gauss_points_quad_p(n, economic)
                }
            }
            3 => match np {
                Some(n) if triangle_rule_tabulated(n) => gauss_points_triangle_p(n),
                Some(n) => gauss_points_p_triangle(n),
                None => p_triangle_rule(n),
            },
            5 => p_tetra_rule(n, np.is_some()),
            6 => gauss_points_p_pyramid(n),
            7 => p_wedge_rule(n, np.is_some()),
            8 => fe(gauss_points_brick(n)),
            _ => Err(Error::Unsupported(format!(
                "p integration rule for type {}",
                el.code
            ))),
        }
    }
}
