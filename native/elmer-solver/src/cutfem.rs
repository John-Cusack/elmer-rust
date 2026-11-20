//! The geometric core of Elmer's CutFEM (CutFEMUtils.F90, with
//! `SplitSingleElement` of ElementGeometry.F90): a zero levelset cuts the
//! 2D mesh's edges, each cut edge gets a new dof (node `nn + edge`) at the
//! intersection, and the elements it crosses are split on the fly into
//! triangles and quadrilaterals on either side, with interface line
//! elements between them.
//!
//! Node indexes are 0-based; a cut-edge node is `nn + edge`. Permutations
//! hold 1-based positions with 0 for "not in", as Elmer's do. The routines
//! that keep `SAVE`d state between calls (`CutInterfaceBulk`,
//! `CutInterfaceBC`, `SplitSingleElement`) take it explicitly, so that a
//! caller sharing one state across loops sees the same carry-over (the
//! split counter continuing when the next loop starts on the element the
//! previous one ended on).
//!
//! Only 2D triangles (303) and quadrilaterals (404) are ported.

// Index loops mirror the Fortran's.
#![allow(clippy::needless_range_loop)]

use crate::error::SolverError;

const CALLER: &str = "CreateCutFEMPerm";

/// The bulk elements a CutFEM split works on: 303 or 404, node indexes
/// 0-based.
#[derive(Clone, Debug)]
pub struct CutMesh {
    pub nn: usize,
    pub codes: Vec<u16>,
    pub nodes: Vec<Vec<u32>>,
    /// Elmer's body index of each bulk element (before CutFEM tampers it).
    pub body: Vec<i32>,
}

/// `FindMeshEdges2D` for the bulk elements: edges `(k, k+1)` of each
/// element (cyclic), numbered in first-use order, each stored with the
/// orientation of the element that first used it.
#[derive(Clone, Debug, Default)]
pub struct Edges {
    pub nodes: Vec<[u32; 2]>,
    /// Per bulk element, its edges in local order.
    pub of_element: Vec<Vec<u32>>,
}

impl Edges {
    pub fn find(mesh: &CutMesh) -> Self {
        let mut map: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
        let mut edges = Edges::default();
        for el in &mesh.nodes {
            let n = el.len();
            let mut ids = Vec::with_capacity(n);
            for k in 0..n {
                let (a, b) = (el[k], el[(k + 1) % n]);
                let key = (a.min(b), a.max(b));
                let id = *map.entry(key).or_insert_with(|| {
                    edges.nodes.push([a, b]);
                    (edges.nodes.len() - 1) as u32
                });
                ids.push(id);
            }
            edges.of_element.push(ids);
        }
        edges
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// The solver keywords `CreateCutFEMPerm` reads.
#[derive(Clone, Debug)]
pub struct PermOptions {
    /// `CutFEM Inside Body` (default: the number of bodies).
    pub body_in: i32,
    /// `CutFEM Outside Body` (default `body_in + 1`).
    pub body_out: i32,
    /// `CutFEM Epsilon` (default 1e-3).
    pub eps: f64,
    /// `CutFEM Epsilon Absolute`.
    pub abs_eps: bool,
    /// `CutFEM extend`.
    pub extend: bool,
    /// `CutFEM extend width`, when given.
    pub extend_width: Option<f64>,
    /// `CutFEM Passive Outside` / `Passive Inside`, when given (defaults:
    /// the body index being 0).
    pub passive_outside: Option<bool>,
    pub passive_inside: Option<bool>,
}

/// What `CreateCutFEMPerm` leaves behind for the step.
#[derive(Clone, Debug)]
pub struct CutState {
    pub nn: usize,
    /// `CutDof(1:nn+ne)`: nodes on the zero levelset, then cut edges.
    pub cut_dof: Vec<bool>,
    /// `CutInterp(1:ne)`: `|h2| / (|h1| + |h2|)` of each cut edge.
    pub cut_interp: Vec<f64>,
    /// `CutPerm(1:nn+ne)`.
    pub cut_perm: Vec<u32>,
    /// `ExtendPerm(1:nn+ne)` with `CutFEM extend`.
    pub extend_perm: Option<Vec<u32>>,
    /// Coordinates of nodes and cut-edge nodes (`Mesh % Nodes` enlarged to
    /// `nn + ne`); an edge node that is not cut keeps what it had.
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    /// The tampered body index of each bulk element (`Element % BodyId`);
    /// elements off the levelset keep theirs.
    pub body: Vec<i32>,
    /// `CutFEMBody`: the one body index the levelset elements had.
    pub cutfem_body: i32,
    /// `body_cut`: the index marking interface elements.
    pub body_cut: i32,
    pub body_in: i32,
    pub body_out: i32,
}

/// `CreateCutFEMPerm(Solver, UpdateCoords = .TRUE.)` for the step. `phi`
/// is the levelset at each node (`None` where `PhiPerm` is 0). `coords`
/// holds the previous step's enlarged coordinates (or the mesh's), so the
/// edge nodes that are not cut keep their values as the Fortran arrays do.
pub fn create_perm(
    mesh: &CutMesh,
    edges: &Edges,
    phi: &[Option<f64>],
    opts: &PermOptions,
    coords: (Vec<f64>, Vec<f64>, Vec<f64>),
    bulk: &mut BulkSplit,
    mut forbidden: impl FnMut(&[bool]) -> bool,
) -> Result<CutState, SolverError> {
    let (nn, ne) = (mesh.nn, edges.len());
    let (mut x, mut y, mut z) = coords;
    for v in [&mut x, &mut y, &mut z] {
        if v.len() < nn + ne {
            v.resize(nn + ne, 0.0);
        }
    }
    let mut cut_dof = vec![false; nn + ne];
    let mut cut_interp = vec![0.0; ne];
    let eps = opts.eps;
    let both = |e: &[u32; 2]| Some((phi[e[0] as usize]?, phi[e[1] as usize]?));
    // First the nodes on (or very near) the zero levelset.
    for e in &edges.nodes {
        let Some((h1, h2)) = both(e) else { continue };
        let hprod = h1 * h2;
        if hprod < 0.0 {
            let r = h2.abs() / (h1.abs() + h2.abs());
            if opts.abs_eps {
                if h2.abs() < eps {
                    cut_dof[e[1] as usize] = true;
                }
                if h1.abs() < eps {
                    cut_dof[e[0] as usize] = true;
                }
            } else {
                if r <= eps {
                    cut_dof[e[1] as usize] = true;
                }
                if 1.0 - r < eps {
                    cut_dof[e[0] as usize] = true;
                }
            }
        } else if hprod.abs() < 1.0e-20 {
            if h1.abs() < 1.0e-20 {
                cut_dof[e[0] as usize] = true;
            }
            if h2.abs() < 1.0e-20 {
                cut_dof[e[1] as usize] = true;
            }
        }
    }
    // Then the cut edges, avoiding ones next to a cut node.
    let max_rat = 2.0;
    for (i, e) in edges.nodes.iter().enumerate() {
        let Some((h1, h2)) = both(e) else { continue };
        if h1 * h2 < 0.0 {
            let r = h2.abs() / (h1.abs() + h2.abs());
            if r <= max_rat * eps {
                if cut_dof[e[1] as usize] {
                    continue;
                }
            } else if 1.0 - r < max_rat * eps && cut_dof[e[0] as usize] {
                continue;
            }
            cut_dof[nn + i] = true;
            cut_interp[i] = r;
            let (a, b) = (e[0] as usize, e[1] as usize);
            x[nn + i] = (1.0 - r) * x[b] + r * x[a];
            y[nn + i] = (1.0 - r) * y[b] + r * y[a];
            z[nn + i] = (1.0 - r) * z[b] + r * z[a];
        }
    }

    let body_in = opts.body_in;
    let body_out = opts.body_out;
    let mut body_cut = body_in.max(body_out);
    if opts.extend {
        body_cut += 1;
    }
    let mut body = mesh.body.clone();
    let mut cutfem_body = 0;
    for (t, el) in mesh.nodes.iter().enumerate() {
        if el.iter().any(|&k| phi[k as usize].is_none()) {
            continue;
        }
        if cutfem_body == 0 {
            cutfem_body = mesh.body[t];
        } else if cutfem_body != mesh.body[t] {
            return Err(SolverError::fatal(CALLER, "Modify code to deal with several bodies!"));
        }
        let mut j = -1;
        if edges.of_element[t].iter().any(|&ed| cut_dof[nn + ed as usize]) {
            j = body_cut;
        } else {
            if mesh.codes[t] / 100 == 4 {
                let c = |k: usize| cut_dof[el[k] as usize];
                if (c(0) && c(2)) || (c(1) && c(3)) {
                    j = body_cut;
                }
            }
            if j < 0 {
                let vals: Vec<f64> = el.iter().map(|&k| phi[k as usize].unwrap()).collect();
                let h1 = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let h2 = vals.iter().cloned().fold(f64::INFINITY, f64::min);
                j = if h1 > -h2 { body_out } else { body_in };
            }
        }
        body[t] = j;
    }

    let mut cut_perm: Vec<u32> = cut_dof.iter().map(|&c| c as u32).collect();
    let passive_outside = opts.passive_outside.unwrap_or(body_out == 0);
    let passive_inside = opts.passive_inside.unwrap_or(body_in == 0);
    if passive_outside {
        for i in 0..nn {
            if phi[i].is_some_and(|p| p < 0.0) {
                cut_perm[i] = 1;
            }
        }
    } else if passive_inside {
        for i in 0..nn {
            if phi[i].is_some_and(|p| p > 0.0) {
                cut_perm[i] = 1;
            }
        }
    } else {
        cut_perm[..nn].iter_mut().for_each(|p| *p = 1);
    }
    let mut j = 0;
    for p in cut_perm.iter_mut() {
        if *p > 0 {
            j += 1;
            *p = j;
        }
    }

    let mut state = CutState {
        nn,
        cut_dof,
        cut_interp,
        cut_perm,
        extend_perm: None,
        x,
        y,
        z,
        body,
        cutfem_body,
        body_cut,
        body_in,
        body_out,
    };

    if opts.extend {
        let mut ext = vec![0u32; nn + ne];
        for (t, el) in mesh.nodes.iter().enumerate() {
            if el.iter().any(|&k| phi[k as usize].is_none()) {
                continue;
            }
            if state.body[t] == body_cut {
                loop {
                    let sub = bulk.next(&state, mesh, edges, t, phi)?;
                    if sub.nodes.iter().any(|&k| state.cut_perm[k as usize] == 0) {
                        for &k in &sub.nodes {
                            ext[k as usize] = 1;
                        }
                    }
                    if !sub.more {
                        break;
                    }
                }
                // Revert the dirty flag.
                state.body[t] = body_cut - 1;
            } else if el.iter().all(|&k| state.cut_perm[k as usize] == 0) {
                match opts.extend_width {
                    Some(r) => {
                        let m = el.iter().map(|&k| phi[k as usize].unwrap().abs()).fold(f64::INFINITY, f64::min);
                        if m < r {
                            for &k in el {
                                ext[k as usize] = 1;
                            }
                        }
                    }
                    None => {
                        for &k in el {
                            ext[k as usize] = 1;
                        }
                    }
                }
            }
        }
        let mut j = 0;
        for p in ext.iter_mut() {
            if *p > 0 {
                j += 1;
                *p = j;
            }
        }
        state.extend_perm = Some(ext);
    }
    // `CutFem Forbidden Boundary`: the caller checks the boundary elements'
    // edges against the cut perm.
    let cut_edges: Vec<bool> = (0..ne).map(|i| state.cut_perm[nn + i] > 0).collect();
    if forbidden(&cut_edges) {
        return Err(SolverError::fatal(CALLER, "CutFEM extends beyond forbidden boundaries!"));
    }
    Ok(state)
}

/// `CutInterfaceCheck(Element, IsCut, IsActive, Perm)`: `(is_cut,
/// is_active)`.
pub fn interface_check(state: &CutState, mesh: &CutMesh, edges: &Edges, t: usize, perm: &[u32]) -> (bool, bool) {
    let nn = state.nn;
    let n_split = edges.of_element[t].iter().filter(|&&e| state.cut_dof[nn + e as usize]).count();
    let n_cut = mesh.nodes[t].iter().filter(|&&k| state.cut_dof[k as usize]).count();
    let n_act = mesh.nodes[t].iter().filter(|&&k| perm[k as usize] > 0).count();
    let is_cut = n_split > 0 || n_cut > 1;
    (is_cut, n_act == mesh.nodes[t].len() || is_cut)
}

/// One piece of a split bulk element.
#[derive(Clone, Debug, PartialEq)]
pub struct SubElement {
    /// Whether the element was split at all; when not, `nodes` are the
    /// element's own and `body` its (tampered) body.
    pub cut: bool,
    /// Element code (303 or 404).
    pub code: u16,
    pub nodes: Vec<u32>,
    pub body: i32,
    /// `IsMore`: another piece follows.
    pub more: bool,
}

/// The `SAVE`d state of `SplitSingleElement`.
#[derive(Clone, Debug, Default)]
pub struct SplitState {
    subcase: i32,
    j: usize,
    j2: usize,
    j3: usize,
    j4: usize,
    mmax: usize,
    s1: f64,
    s2: f64,
}

/// The pieces' local indexes (1-based into nodes then edge nodes) and the
/// node whose levelset sign gives the side.
#[derive(Clone, Debug, Default)]
pub struct Split {
    pub is_cut: bool,
    pub is_more: bool,
    pub local: [usize; 4],
    pub sgn_node: usize,
}

/// `SplitSingleElement(Element, ElemCut, ElemNodes, m, ...)` for 303/404.
/// `cut` and the coordinates are the element's nodes then its edge nodes
/// (`2n` entries); `m` is the 1-based piece counter.
pub fn split_single_element(
    st: &mut SplitState,
    code: u16,
    cut: &[bool],
    ex: &[f64],
    ey: &[f64],
    ez: &[f64],
    m: usize,
) -> Result<Split, SolverError> {
    let n = (code / 100) as usize;
    let n_split = cut[n..2 * n].iter().filter(|&&c| c).count();
    let n_cut = cut[..n].iter().filter(|&&c| c).count();
    let mut out = Split::default();
    if n_split == 0 {
        return Ok(out);
    }
    out.is_cut = true;
    let case = 100 * code as usize + 10 * n_split + n_cut;
    // 1-based helpers over the element arrays.
    let c = |i: usize| cut[i - 1];
    let d2 = |a: usize, b: usize| {
        sq(ex[a - 1] - ex[b - 1]) + sq(ey[a - 1] - ey[b - 1]) + sq(ez[a - 1] - ez[b - 1])
    };
    let modp = |a: usize, b: usize| a % b + 1; // MODULO(a, b) + 1 for a >= 0
    let li = &mut out.local;
    match case {
        30320 | 30321 => {
            if m == 1 {
                let mut j = 1;
                while j <= 3 && c(n + j) {
                    j += 1;
                }
                st.j = j;
                st.j2 = modp(j, 3);
                st.j3 = modp(j + 1, 3);
                st.mmax = 3;
                st.s1 = d2(st.j, n + st.j2);
                st.s2 = d2(st.j2, n + st.j3);
                li[0] = st.j;
                li[1] = st.j2;
                li[2] = if st.s1 < st.s2 { n + st.j2 } else { n + st.j3 };
                out.sgn_node = 1;
            } else if m == 2 {
                li[0] = if st.s1 < st.s2 { st.j } else { st.j2 };
                li[1] = n + st.j2;
                li[2] = n + st.j3;
                out.sgn_node = 1;
            } else if m == 3 {
                li[0] = n + st.j3;
                li[1] = n + st.j2;
                li[2] = st.j3;
                out.sgn_node = 3;
            }
        }
        30311 => {
            if m == 1 {
                let mut j = 1;
                while j <= 3 && !c(n + j) {
                    j += 1;
                }
                st.j = j;
                st.j2 = modp(j, 3);
                st.j3 = modp(j + 1, 3);
            }
            let (j, j2, j3) = (st.j, st.j2, st.j3);
            if c(j3) {
                if m == 1 {
                    *li = [n + j, j2, j3, 0];
                    out.sgn_node = 2;
                    st.mmax = 2;
                } else if m == 2 {
                    *li = [n + j, j3, j, 0];
                    out.sgn_node = 3;
                }
            } else if c(j) || c(j2) {
                *li = [1, 2, 3, 0];
                out.sgn_node = j3;
                st.mmax = 1;
            }
        }
        40420 | 40421 => {
            if m == 1 {
                st.subcase = 0;
                if c(n + 1) && c(n + 3) {
                    st.subcase = 1;
                    st.j = 1;
                    st.mmax = 2;
                } else if c(n + 2) && c(n + 4) {
                    st.subcase = 1;
                    st.j = 2;
                    st.mmax = 2;
                } else {
                    let mut j = 1;
                    while j <= 4 {
                        st.j2 = modp(j, 4);
                        if c(n + j) && c(n + st.j2) {
                            st.subcase = 2;
                            st.mmax = 3;
                            break;
                        }
                        j += 1;
                    }
                    st.j = j;
                }
                if st.subcase == 0 {
                    return Err(SolverError::fatal("SplitSingleElement", "This case not treated yet for 404!"));
                }
            }
            let j = st.j;
            if st.subcase == 1 {
                st.mmax = 2;
                if m == 1 {
                    st.j2 = modp(j, 4);
                    st.j3 = modp(j + 1, 4);
                    st.j4 = modp(j + 2, 4);
                    *li = [j, n + j, n + st.j3, st.j4];
                    out.sgn_node = 1;
                } else if m == 2 {
                    *li = [st.j2, st.j3, n + st.j3, n + j];
                    out.sgn_node = 1;
                }
            } else if st.subcase == 2 {
                st.mmax = 4;
                if m == 1 {
                    st.j2 = modp(j, 4);
                    st.j3 = modp(j + 1, 4);
                    st.j4 = modp(j + 2, 4);
                    *li = [n + j, st.j2, n + st.j2, 0];
                    out.sgn_node = 2;
                } else if m == 2 {
                    *li = [j, n + j, st.j4, 0];
                    out.sgn_node = 3;
                } else if m == 3 {
                    *li = [n + j, n + st.j2, st.j4, 0];
                    out.sgn_node = 3;
                } else if m == 4 {
                    *li = [n + st.j2, st.j3, st.j4, 0];
                    out.sgn_node = 3;
                }
            }
        }
        40411 => {
            let mut j = 1;
            while j <= 4 && !c(n + j) {
                j += 1;
            }
            st.j = j;
            st.j2 = modp(j, 4);
            st.j3 = modp(j + 1, 4);
            st.j4 = modp(j + 2, 4);
            let (j2, j3, j4) = (st.j2, st.j3, st.j4);
            if c(j) || c(j2) {
                *li = [1, 2, 3, 4];
                out.sgn_node = j3;
                st.mmax = 1;
            } else {
                st.mmax = 2;
                if c(j3) {
                    if m == 1 {
                        *li = [n + j, j2, j3, j4];
                        out.sgn_node = 3;
                    } else if m == 2 {
                        *li = [j, n + j, j4, 0];
                        out.sgn_node = 1;
                    }
                } else if c(j4) {
                    if m == 1 {
                        *li = [j, n + j, j3, j4];
                        out.sgn_node = 4;
                    } else if m == 2 {
                        *li = [n + j, j2, j3, 0];
                        out.sgn_node = 2;
                    }
                }
            }
        }
        _ => {
            return Err(SolverError::fatal(
                "SplitSingleElement",
                format!("Unknown split case in element divisions: {case}"),
            ))
        }
    }
    out.is_more = m < st.mmax;
    Ok(out)
}

/// The `SAVE`d state of `CutInterfaceBulk` (and the `SplitSingleElement`
/// state it drives). The element's nodes, cut flags and coordinates are
/// taken when the element changes and kept for its later pieces, as the
/// Fortran's `ElemInds`, `ElemCut` and `ElemNodes` are; every loop that
/// splits elements shares the one state.
#[derive(Clone, Debug, Default)]
pub struct BulkSplit {
    prev: Option<usize>,
    count: usize,
    split: SplitState,
    code: u16,
    inds: Vec<u32>,
    cut: Vec<bool>,
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
}

impl BulkSplit {
    /// `CutInterfaceBulk(Element, IsCut, IsMore)` for bulk element `t`.
    pub fn next(&mut self, st: &CutState, mesh: &CutMesh, edges: &Edges, t: usize, phi: &[Option<f64>]) -> Result<SubElement, SolverError> {
        if self.prev != Some(t) {
            self.count = 1;
            self.prev = Some(t);
            self.code = mesh.codes[t];
            let n = (self.code / 100) as usize;
            self.inds = mesh.nodes[t][..n]
                .iter()
                .copied()
                .chain(edges.of_element[t][..n].iter().map(|&e| (st.nn + e as usize) as u32))
                .collect();
            self.cut = self.inds.iter().map(|&k| st.cut_dof[k as usize]).collect();
            self.ex = self.inds.iter().map(|&k| st.x[k as usize]).collect();
            self.ey = self.inds.iter().map(|&k| st.y[k as usize]).collect();
            self.ez = self.inds.iter().map(|&k| st.z[k as usize]).collect();
        } else {
            self.count += 1;
        }
        let code = self.code;
        let inds = &self.inds;
        let s = split_single_element(&mut self.split, code, &self.cut, &self.ex, &self.ey, &self.ez, self.count)?;
        if !s.is_cut {
            return Ok(SubElement { cut: false, code, nodes: mesh.nodes[t].clone(), body: st.body[t], more: false });
        }
        let k = s.local.iter().filter(|&&l| l > 0).count();
        let code = match k {
            3 => 303,
            4 => 404,
            _ => return Err(SolverError::fatal("CutInterfaceBulk", "Impossible number of nodes!")),
        };
        let nodes: Vec<u32> = s.local[..k].iter().map(|&l| inds[l - 1]).collect();
        let body = if st.body_out == 0 {
            if nodes.iter().all(|&k| st.cut_perm[k as usize] > 0) {
                st.body_in
            } else {
                st.body_out
            }
        } else {
            let p = phi[nodes[s.sgn_node - 1] as usize]
                .ok_or_else(|| SolverError::fatal("CutInterfaceBulk", "the sign node has no levelset value"))?;
            if p > 0.0 {
                st.body_out
            } else {
                st.body_in
            }
        };
        Ok(SubElement { cut: true, code, nodes, body, more: s.is_more })
    }
}

/// The `SAVE`d state of `CutInterfaceBC`.
#[derive(Clone, Debug, Default)]
pub struct BcSplit {
    prev: Option<usize>,
    m: usize,
    n_split: usize,
    n_cut: usize,
}

/// One interface line element (202) of a split element.
#[derive(Clone, Debug, PartialEq)]
pub struct InterfaceElement {
    pub nodes: [u32; 2],
    pub more: bool,
}

impl BcSplit {
    /// `CutInterfaceBC(Element, IsCut, IsMore)` for bulk element `t` of a
    /// 2D mesh: `None` when it is not cut.
    pub fn next(
        &mut self,
        st: &CutState,
        mesh: &CutMesh,
        edges: &Edges,
        t: usize,
        phi: &[Option<f64>],
    ) -> Result<Option<InterfaceElement>, SolverError> {
        let nn = st.nn;
        let ni = &mesh.nodes[t];
        let ei = &edges.of_element[t];
        let ecut = |j: usize| st.cut_dof[nn + ei[j - 1] as usize];
        let ncut = |j: usize| st.cut_dof[ni[j - 1] as usize];
        if self.prev != Some(t) {
            self.m = 1;
            self.prev = Some(t);
            self.n_split = ei.iter().filter(|&&e| st.cut_dof[nn + e as usize]).count();
            self.n_cut = ni.iter().filter(|&&k| st.cut_dof[k as usize]).count();
        } else {
            self.m += 1;
        }
        if self.n_split == 0 && self.n_cut <= 1 {
            return Ok(None);
        }
        let code = mesh.codes[t] as usize;
        let ne = match code {
            303 => 3,
            404 => 4,
            _ => {
                return Err(SolverError::fatal("CutInterfaceBC", format!("Unknown element type to split: {code}!")))
            }
        };
        let first = |f: &dyn Fn(usize) -> bool, skip: &[usize]| -> usize {
            let mut j = 1;
            while j <= ne && (skip.contains(&j) || !f(j)) {
                j += 1;
            }
            j
        };
        let edge_node = |j: usize| (nn + ei[j - 1] as usize) as u32;
        let case = 100 * code + 10 * self.n_split + self.n_cut;
        let mut more = false;
        let mut p: [u32; 2] = match case {
            30320 | 40420 => {
                let j = first(&ecut, &[]);
                let j2 = first(&ecut, &[j]);
                [edge_node(j), edge_node(j2)]
            }
            30321 | 40421 => {
                if self.m == 1 {
                    let j = first(&ecut, &[]);
                    let j2 = first(&ecut, &[j]);
                    more = true;
                    [edge_node(j), edge_node(j2)]
                } else {
                    let j = first(&ncut, &[]);
                    let mut j2 = j;
                    if !ecut(j2) {
                        j2 = (j + ne - 2) % ne + 1; // MODULO(j-2, ne) + 1
                        if !ecut(j2) {
                            return Err(SolverError::fatal("Caller", format!("Could not imagine this {code} case!")));
                        }
                    }
                    [ni[j - 1], edge_node(j2)]
                }
            }
            30311 | 40411 => {
                let j = first(&ecut, &[]);
                let j2 = first(&ncut, &[]);
                [edge_node(j), ni[j2 - 1]]
            }
            30302 | 40402 => {
                let j = first(&ncut, &[]);
                let j2 = first(&ncut, &[j]);
                [ni[j - 1], ni[j2 - 1]]
            }
            _ => {
                return Err(SolverError::fatal(
                    "CutInterfaceBC",
                    format!("Unknown split case in bc element divisions: {case}"),
                ))
            }
        };
        // Orient the segment by the node with the extreme levelset value.
        let (mut pmax, mut imax) = (0.0f64, None);
        for &k in ni {
            let Some(v) = phi[k as usize] else { continue };
            if v.abs() > pmax.abs() {
                pmax = v;
                imax = Some(k as usize);
            }
        }
        if let Some(im) = imax {
            let (x0, y0) = (st.x[p[0] as usize], st.y[p[0] as usize]);
            let (x1, y1) = (st.x[p[1] as usize], st.y[p[1] as usize]);
            let (xp, yp) = (st.x[im], st.y[im]);
            let dir1 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0);
            if dir1 * pmax < 0.0 {
                p.swap(0, 1);
            }
        }
        Ok(Some(InterfaceElement { nodes: p, more }))
    }
}

/// `CreateCutFemMatrix(Solver, Perm)`: the sorted CRS pattern of the dofs
/// `perm` gives, every pair of an element's active nodes and cut edges
/// coupled (`dofs` per node, `dofs = 1` only). Elements with a node
/// outside the levelset field are skipped.
pub fn matrix_pattern(mesh: &CutMesh, edges: &Edges, phi: &[Option<f64>], perm: &[u32]) -> (Vec<usize>, Vec<u32>) {
    let n = perm.iter().copied().max().unwrap_or(0) as usize;
    let mut rows: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); n];
    let nn = mesh.nn;
    for (t, el) in mesh.nodes.iter().enumerate() {
        if el.iter().any(|&k| phi[k as usize].is_none()) {
            continue;
        }
        let mut b: Vec<u32> = el.iter().map(|&k| perm[k as usize]).filter(|&j| j > 0).collect();
        b.extend(edges.of_element[t].iter().map(|&e| perm[nn + e as usize]).filter(|&j| j > 0));
        for &i in &b {
            for &j in &b {
                rows[i as usize - 1].insert(j - 1);
            }
        }
    }
    // The extreme entry `List_AddToMatrixElement(n, n, 0)`.
    if n > 0 {
        rows[n - 1].insert(n as u32 - 1);
    }
    let mut ptr = vec![0usize];
    let mut cols = Vec::new();
    for r in rows {
        cols.extend(r);
        ptr.push(cols.len());
    }
    (ptr, cols)
}

// ---------------------------------------------------------------------
// The run: `DefaultStart` / `DefaultFinish` (DefUtils.F90:4005, :4290)
// with `CutFEM Solver` (the on-the-fly splitting solvers).
// ---------------------------------------------------------------------

use crate::sim::LinearSystem;
use crate::variable::VarId;
use crate::solver::SolverCtx;
use elmer_core::{CsrMatrix, ElemId, LocalSystem, ParamMap, Params, NO_DOF};
use std::sync::Arc;

/// A mesh a CutFEM solver adds to the output for the step (`Mesh % Next`,
/// `Mesh % Next % Next`): `<mesh name><suffix>`, its node count and fields
/// (name, values per node, equation of the solver), saved from scratch at
/// each save as Elmer saves a mesh it has just allocated.
#[derive(Clone, Debug)]
pub struct ExtraMesh {
    pub suffix: String,
    pub n_nodes: usize,
    pub fields: Vec<(String, Vec<f64>, String)>,
}

/// The isoline mesh of a step: node coordinates, line elements (0-based
/// node pairs) and fields.
struct IsoMesh {
    x: Vec<f64>,
    y: Vec<f64>,
    elems: Vec<[usize; 2]>,
    fields: Vec<(String, Vec<f64>, String)>,
}

/// One polyline segment of `PopulatePolyline`: `x0 x1 y0 y1 phi0 phi1`,
/// then two columns per isoline variable; the neighbouring segments.
struct PolyLine {
    vals: Vec<f64>,
    #[allow(dead_code)]
    prev: i64,
    #[allow(dead_code)]
    next: i64,
}

/// Fortran's `x**2`.
fn sq(x: f64) -> f64 {
    x * x
}

/// `LineSegmentsIntersect(a1, a2, b1, b2)` with the default buffer `AEPS`.
fn segments_intersect(a1: [f64; 2], a2: [f64; 2], b1: [f64; 2], b2: [f64; 2]) -> bool {
    let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    let r = [a2[0] - a1[0], a2[1] - a1[1]];
    let s = [b2[0] - b1[0], b2[1] - b1[1]];
    let rxs = cross(r, s);
    if rxs == 0.0 {
        return false;
    }
    let bma = [b1[0] - a1[0], b1[1] - a1[1]];
    let t = cross(bma, s) / rxs;
    let u = cross(bma, r) / rxs;
    let e = 10.0 * f64::EPSILON;
    !(t < 0.0 - e || t > 1.0 + e || u < 0.0 - e || u > 1.0 + e)
}

/// The original variable and system while the cut ones are in place.
struct Saved {
    perm: Option<Arc<[u32]>>,
    values: Vec<f64>,
    prev: Vec<Vec<f64>>,
    system: Option<LinearSystem>,
}

/// The `CutFEMUtils` module state of one solver.
pub struct CutFemRun {
    mesh: CutMesh,
    edges: Edges,
    edge_of: std::collections::HashMap<(u32, u32), u32>,
    levelset: String,
    coords: Option<(Vec<f64>, Vec<f64>, Vec<f64>)>,
    state: Option<CutState>,
    phi: Vec<Option<f64>>,
    bulk: BulkSplit,
    bc: BcSplit,
    saved: Option<Saved>,
    opts: PermOptions,
    dofs: usize,
    /// `CreateCutFEMAddMesh` of a library solver (no `CutFEM Solver`).
    lib: Option<LibraryMesh>,
}

/// What `CreateCutFEMAddMesh` gives a library solver for the step: the
/// unsplit active elements it assembles first (`UnsplitActiveElements`),
/// then, on the add mesh, the pieces of the split elements and the
/// interface elements.
#[derive(Clone, Debug, Default)]
pub struct LibraryMesh {
    pub unsplit: Vec<ElemId>,
    pub pieces: Vec<SubElement>,
    pub interface: Vec<[u32; 2]>,
}

fn cut_err(msg: impl Into<String>) -> SolverError {
    SolverError::fatal("CutFEM", msg)
}

impl SolverCtx<'_> {
    fn cut(&self) -> Result<&CutFemRun, SolverError> {
        self.sim.slots[self.slot].cutfem.as_deref().ok_or_else(|| cut_err("CutFEM has not started"))
    }

    fn cut_state(&self) -> Result<&CutState, SolverError> {
        self.cut()?.state.as_ref().ok_or_else(|| cut_err("no CutFEM permutation"))
    }

    /// `CreateCutFEMPerm(Solver, .TRUE.)` and `CreateCutFEMMatrix`: the
    /// variable takes the cut dofs (nodes and cut edges), values and
    /// history copied at the nodes (zero at the edges), and the solver a
    /// fresh zero system on the cut pattern (with a global mass matrix).
    pub(crate) fn cutfem_start(&mut self) -> Result<(), SolverError> {
        let p = self.params().clone();
        for k in ["cutfem extrapolate", "cutfem bodyfitted"] {
            if p.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("CutFEM option {k}")));
            }
        }
        if self.mesh.dim() != 2 {
            return Err(SolverError::Unsupported("CutFEM outside 2D".into()));
        }
        let id = self.sim.slots[self.slot].variable.ok_or_else(|| cut_err("solver has no variable"))?;
        if self.sim.slots[self.slot].cutfem.is_none() {
            let nb = self.mesh.n_bulk();
            let mut mesh = CutMesh { nn: self.mesh.n_nodes(), codes: Vec::new(), nodes: Vec::new(), body: Vec::new() };
            for t in 0..nb {
                let el = self.mesh.element(ElemId(t as u32));
                if !matches!(el.etype.0, 303 | 404) {
                    return Err(SolverError::Unsupported(format!("CutFEM on element type {}", el.etype.0)));
                }
                mesh.codes.push(el.etype.0);
                mesh.nodes.push(el.nodes.to_vec());
                mesh.body.push(self.sim.model.body_of_tag(el.tag).map_or(0, |b| b as i32));
            }
            let edges = Edges::find(&mesh);
            let edge_of = edges.nodes.iter().enumerate().map(|(i, e)| ((e[0].min(e[1]), e[0].max(e[1])), i as u32)).collect();
            let nbodies = self.sim.model.bodies.len() as i32;
            let body_in = p.integer("cutfem inside body").unwrap_or(nbodies);
            let opts = PermOptions {
                body_in,
                body_out: p.integer("cutfem outside body").unwrap_or(body_in + 1),
                eps: p.real("cutfem epsilon").unwrap_or(1.0e-3),
                abs_eps: p.logical("cutfem epsilon absolute") == Some(true),
                extend: p.logical("cutfem extend") == Some(true),
                extend_width: p.real("cutfem extend width"),
                passive_outside: p.logical("cutfem passive outside"),
                passive_inside: p.logical("cutfem passive inside"),
            };
            let levelset = elmer_core::normalize_key(p.string("levelset variable").unwrap_or("surface"));
            let dofs = self.sim.vars.by_id(id).dofs;
            if dofs != 1 {
                return Err(SolverError::Unsupported("CutFEM for a variable of more than one dof".into()));
            }
            self.sim.slots[self.slot].cutfem = Some(Box::new(CutFemRun {
                mesh,
                edges,
                edge_of,
                levelset,
                coords: None,
                state: None,
                phi: Vec::new(),
                bulk: BulkSplit::default(),
                bc: BcSplit::default(),
                saved: None,
                opts,
                dofs,
                lib: None,
            }));
        }
        let mesh_view = self.mesh;
        let nn = mesh_view.n_nodes();
        // The levelset at the nodes.
        let phi: Vec<Option<f64>> = {
            let run = self.cut()?;
            let pv = self.sim.vars.get(&run.levelset).ok_or_else(|| {
                cut_err(format!("\"Levelset Variable\" not available: {}", run.levelset))
            })?;
            (0..nn as u32).map(|k| pv.slot(k).map(|s| pv.values[s])).collect()
        };
        // `CutFem Forbidden Boundary`: boundary elements of such a condition
        // must not have a cut edge.
        let forbidden: Vec<[u32; 2]> = self
            .boundary_elements()
            .filter(|&e| self.bc_params(e).is_some_and(|b| b.logical("cutfem forbidden boundary") == Some(true)))
            .filter_map(|e| {
                let n = mesh_view.element(e).nodes;
                (n.len() >= 2).then(|| [n[0], n[1]])
            })
            .collect();
        let sim = &mut *self.sim;
        let slot = &mut sim.slots[self.slot];
        let run = slot.cutfem.as_deref_mut().expect("created above");
        let coords = run.coords.take().unwrap_or_else(|| {
            let c: Vec<[f64; 3]> = (0..nn as u32).map(|k| mesh_view.coord(k)).collect();
            (c.iter().map(|p| p[0]).collect(), c.iter().map(|p| p[1]).collect(), c.iter().map(|p| p[2]).collect())
        });
        let edge_of = &run.edge_of;
        let st = create_perm(&run.mesh, &run.edges, &phi, &run.opts, coords, &mut run.bulk, |cut| {
            forbidden.iter().any(|f| edge_of.get(&(f[0].min(f[1]), f[0].max(f[1]))).is_some_and(|&e| cut[e as usize]))
        })?;
        // CutValues / PrevCutValues from the original perm at the nodes.
        let var = sim.vars.by_id_mut(id);
        let ncut = st.cut_perm.iter().copied().max().unwrap_or(0) as usize;
        let dofs = run.dofs;
        let mut values = vec![0.0; dofs * ncut];
        let mut prev = vec![vec![0.0; dofs * ncut]; var.prev_values.len()];
        for i in 0..nn {
            let j = st.cut_perm[i];
            if j == 0 {
                continue;
            }
            let Some(k) = var.slot(i as u32) else { continue };
            for l in 0..dofs {
                values[dofs * (j as usize - 1) + l] = var.values[dofs * k + l];
                for (pc, po) in prev.iter_mut().zip(&var.prev_values) {
                    pc[dofs * (j as usize - 1) + l] = po[dofs * k + l];
                }
            }
        }
        let perm: Arc<[u32]> = st.cut_perm.iter().map(|&j| if j == 0 { NO_DOF } else { j - 1 }).collect();
        let saved = Saved {
            perm: var.perm.replace(Arc::clone(&perm)),
            values: std::mem::replace(&mut var.values, values),
            prev: std::mem::replace(&mut var.prev_values, prev),
            system: slot.system.take(),
        };
        let (row_ptr, col_idx) = matrix_pattern(&run.mesh, &run.edges, &phi, &st.cut_perm);
        let n = dofs * ncut;
        let diag: Vec<usize> = (0..n).map(|r| (row_ptr[r]..row_ptr[r + 1]).find(|&k| col_idx[k] as usize == r).unwrap_or(row_ptr[r])).collect();
        let nnz = col_idx.len();
        let matrix = CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values: vec![0.0; nnz], diag: Some(diag) };
        let global_mass = slot.params.logical("use global mass matrix") == Some(true);
        let mass = global_mass.then(|| matrix.clone());
        slot.system = Some(LinearSystem { matrix, rhs: vec![0.0; n], mass, damp: None });
        slot.perm = Some(perm);
        slot.asm = Default::default();
        run.saved = Some(saved);
        run.phi = phi;
        run.state = Some(st);
        run.lib = None;
        if p.logical("cutfem solver") != Some(true) {
            let lib = self.cutfem_add_mesh()?;
            let run = self.sim.slots[self.slot].cutfem.as_deref_mut().expect("created above");
            run.lib = Some(lib);
        }
        Ok(())
    }

    /// `CreateCutFEMAddMesh`: the unsplit active elements, then
    /// `CreateCutFEMMesh(..., AddMeshMode = .TRUE.)`'s two sweeps over the
    /// split elements' pieces and interface elements whose nodes are all
    /// cut dofs (the second sweep's are kept).
    fn cutfem_add_mesh(&mut self) -> Result<LibraryMesh, SolverError> {
        let nb = self.mesh.n_bulk();
        let mut lib = LibraryMesh::default();
        {
            let run = self.cut()?;
            let st = self.cut_state()?;
            for t in 0..nb {
                if run.mesh.nodes[t].iter().any(|&k| run.phi[k as usize].is_none()) {
                    continue;
                }
                let (is_cut, is_active) = interface_check(st, &run.mesh, &run.edges, t, &st.cut_perm);
                if is_active && !is_cut {
                    lib.unsplit.push(ElemId(t as u32));
                }
            }
        }
        for sweep in 0..2 {
            for t in 0..nb {
                let e = ElemId(t as u32);
                let (is_cut, _) = self.cutfem_check(e)?;
                if !is_cut {
                    continue;
                }
                loop {
                    let sub = self.cutfem_bulk(e)?;
                    if sweep == 1 && self.cutfem_in_perm(&sub.nodes)? {
                        lib.pieces.push(sub.clone());
                    }
                    if !sub.more {
                        break;
                    }
                }
            }
            for t in 0..nb {
                let e = ElemId(t as u32);
                while let Some(ie) = self.cutfem_interface(e)? {
                    if sweep == 1 && self.cutfem_in_perm(&ie.nodes)? {
                        lib.interface.push(ie.nodes);
                    }
                    if !ie.more {
                        break;
                    }
                }
            }
        }
        Ok(lib)
    }

    /// The variable's own permutation and size while the cut ones are in
    /// place; `None` otherwise.
    pub(crate) fn cutfem_original(&self) -> Option<(Option<Arc<[u32]>>, usize)> {
        let saved = self.sim.slots[self.slot].cutfem.as_deref()?.saved.as_ref()?;
        Some((saved.perm.clone(), saved.values.len()))
    }

    /// The library solver's elements for the step ([`LibraryMesh`]).
    pub fn cutfem_library_mesh(&self) -> Result<LibraryMesh, SolverError> {
        self.cut()?.lib.clone().ok_or_else(|| cut_err("no CutFEM add mesh (a CutFEM Solver?)"))
    }

    /// `CutInterfaceCheck(Element, IsCut, IsActive)` against the variable's
    /// (cut) permutation.
    pub fn cutfem_check(&self, e: ElemId) -> Result<(bool, bool), SolverError> {
        let run = self.cut()?;
        let st = self.cut_state()?;
        Ok(interface_check(st, &run.mesh, &run.edges, e.0 as usize, &st.cut_perm))
    }

    /// `CutInterfaceBulk(Element, IsCut, IsMore)`.
    pub fn cutfem_bulk(&mut self, e: ElemId) -> Result<SubElement, SolverError> {
        let run = self.sim.slots[self.slot].cutfem.as_deref_mut().ok_or_else(|| cut_err("CutFEM has not started"))?;
        let st = run.state.as_ref().ok_or_else(|| cut_err("no CutFEM permutation"))?;
        run.bulk.next(st, &run.mesh, &run.edges, e.0 as usize, &run.phi)
    }

    /// `CutInterfaceBC(Element, IsCut, IsMore)`.
    pub fn cutfem_interface(&mut self, e: ElemId) -> Result<Option<InterfaceElement>, SolverError> {
        let run = self.sim.slots[self.slot].cutfem.as_deref_mut().ok_or_else(|| cut_err("CutFEM has not started"))?;
        let st = run.state.as_ref().ok_or_else(|| cut_err("no CutFEM permutation"))?;
        run.bc.next(st, &run.mesh, &run.edges, e.0 as usize, &run.phi)
    }

    /// The (tampered) body index of bulk element `e` for the step.
    pub fn cutfem_body(&self, e: ElemId) -> Result<i32, SolverError> {
        Ok(self.cut_state()?.body[e.0 as usize])
    }

    /// Coordinates of nodes and cut-edge nodes.
    pub fn cutfem_xyz(&self, nodes: &[u32]) -> Result<Vec<[f64; 3]>, SolverError> {
        let st = self.cut_state()?;
        Ok(nodes.iter().map(|&k| [st.x[k as usize], st.y[k as usize], st.z[k as usize]]).collect())
    }

    /// Whether every node is a cut dof (`ALL(Perm(NodeIndexes) > 0)`).
    pub fn cutfem_in_perm(&self, nodes: &[u32]) -> Result<bool, SolverError> {
        let st = self.cut_state()?;
        Ok(nodes.iter().all(|&k| st.cut_perm[k as usize] > 0))
    }

    /// The boundary condition of the interface elements (`CutFem
    /// Interface BC`).
    pub fn cutfem_interface_bc(&self) -> Option<&ParamMap> {
        let n = self.params().integer("cutfem interface bc")?;
        self.sim.model.boundary_conditions.get((n as usize).checked_sub(1)?)
    }

    fn cut_slots(&self, nodes: &[u32]) -> Result<Vec<u32>, SolverError> {
        let st = self.cut_state()?;
        Ok(nodes.iter().map(|&k| st.cut_perm[k as usize].checked_sub(1).unwrap_or(NO_DOF)).collect())
    }

    /// `UpdateGlobalEquations(A, STIFF, A % rhs, FORCE, ..., Perm(NodeIndexes))`
    /// for a piece of a split element (or an interface element).
    pub fn cutfem_update(&mut self, local: &LocalSystem, nodes: &[u32]) -> Result<(), SolverError> {
        let idx = self.cut_slots(nodes)?;
        let sim = &mut *self.sim;
        let backend = sim.backends.assembly.as_ref().ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
        let sys = sim.slots[self.slot].system.as_mut().ok_or_else(|| cut_err("no CutFEM system"))?;
        backend.add_local(&mut sys.matrix, &mut sys.rhs, &idx, local);
        Ok(())
    }

    /// `UpdateMassMatrix(A, MASS, ..., A % MassValues)` for a piece.
    pub fn cutfem_update_mass(&mut self, mass: &[f64], nodes: &[u32]) -> Result<(), SolverError> {
        let idx = self.cut_slots(nodes)?;
        let n = nodes.len();
        let local = LocalSystem { n, stiff: mass.to_vec(), mass: None, damp: None, force: vec![0.0; n] };
        let sim = &mut *self.sim;
        let backend = sim.backends.assembly.as_ref().ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
        let sys = sim.slots[self.slot].system.as_mut().ok_or_else(|| cut_err("no CutFEM system"))?;
        let Some(m) = sys.mass.as_mut() else { return Ok(()) };
        let mut rhs = vec![0.0; m.nrows];
        backend.add_local(m, &mut rhs, &idx, &local);
        Ok(())
    }

    /// `CreateCutFEMMesh(...)`'s `InterpolateLevelsetVariables` for a mesh
    /// whose node `MeshPerm(i)` (1-based, 0 = not in it) is node or cut edge
    /// `i`: the solver variable copied through its (cut) perm over nodes and
    /// edges, then the levelset and the `<prefix> k` variables at the nodes
    /// with their cut edges interpolated by `CutInterp`.
    fn cut_mesh_fields(&self, mesh_perm: &[u32], prefix: &str) -> Result<Vec<(String, Vec<f64>, String)>, SolverError> {
        let run = self.cut()?;
        let st = self.cut_state()?;
        let nn = st.nn;
        let n = mesh_perm.iter().copied().max().unwrap_or(0) as usize;
        let id = self.sim.slots[self.slot].variable.ok_or_else(|| cut_err("solver has no variable"))?;
        let eq = elmer_core::normalize_key(self.params().string("equation").unwrap_or(""));
        let solver_var = self.sim.vars.by_id(id);
        let mut names = vec![solver_var.name.clone(), run.levelset.clone()];
        for k in 1.. {
            let Some(v) = self.params().string(&format!("{prefix} {k}")) else { break };
            let v = elmer_core::normalize_key(v);
            if v == solver_var.name || v == run.levelset {
                continue;
            }
            names.push(v);
        }
        let mut out = Vec::new();
        for (iv, name) in names.iter().enumerate() {
            let var = self.sim.vars.get(name).ok_or_else(|| {
                SolverError::fatal("InterpolateLevelsetVariable", format!("Could not find variable in 2D mesh: {name}"))
            })?;
            let dofs = var.dofs;
            let mut values = vec![0.0; dofs * n];
            let ntot = if iv == 0 { mesh_perm.len() } else { nn };
            for i in 0..ntot {
                let k = mesh_perm[i];
                if k == 0 {
                    continue;
                }
                let Some(j) = var.slot(i as u32) else { continue };
                for l in 0..dofs {
                    values[dofs * (k as usize - 1) + l] = var.values[dofs * j + l];
                }
            }
            if iv > 0 {
                for (i, e) in run.edges.nodes.iter().enumerate() {
                    let k = mesh_perm[nn + i];
                    if k == 0 {
                        continue;
                    }
                    let k = k as usize - 1;
                    let r = st.cut_interp[i];
                    values[k] = 0.0;
                    if let Some(j) = var.slot(e[0]) {
                        for l in 0..dofs {
                            values[dofs * k + l] = r * var.values[dofs * j + l];
                        }
                    }
                    if let Some(j) = var.slot(e[1]) {
                        for l in 0..dofs {
                            values[dofs * k + l] += (1.0 - r) * var.values[dofs * j + l];
                        }
                    }
                }
            }
            out.push((name.clone(), values, eq.clone()));
        }
        Ok(out)
    }

    /// `CreateCutFEMMesh(..., CreateBC = .TRUE., CreateBulk = .FALSE.)`:
    /// the isoline mesh of the cut nodes and cut edges (numbered in index
    /// order), its line elements the interface elements whose nodes are all
    /// cut dofs, and its fields.
    fn cutfem_isomesh(&mut self) -> Result<IsoMesh, SolverError> {
        let (iso, x, y) = {
            let st = self.cut_state()?;
            let mut j = 0;
            let iso: Vec<u32> = st.cut_dof.iter().map(|&c| if c { j += 1; j } else { 0 }).collect();
            let mut x = vec![0.0; j as usize];
            let mut y = vec![0.0; j as usize];
            for (i, &k) in iso.iter().enumerate() {
                if k > 0 {
                    x[k as usize - 1] = st.x[i];
                    y[k as usize - 1] = st.y[i];
                }
            }
            (iso, x, y)
        };
        // CreateCutFEMMesh counts in a first sweep and fills in a second.
        let mut elems = Vec::new();
        for sweep in 0..2 {
            for t in 0..self.mesh.n_bulk() {
                let e = ElemId(t as u32);
                while let Some(ie) = self.cutfem_interface(e)? {
                    if sweep == 1 && self.cutfem_in_perm(&ie.nodes)? {
                        elems.push([iso[ie.nodes[0] as usize] as usize - 1, iso[ie.nodes[1] as usize] as usize - 1]);
                    }
                    if !ie.more {
                        break;
                    }
                }
            }
        }
        let fields = self.cut_mesh_fields(&iso, "isoline variable")?;
        Ok(IsoMesh { x, y, elems, fields })
    }

    /// The `CreateCutFEMMesh(..., CreateBC, CreateBulk)` loops of the cut
    /// mesh: their pieces only go to the VTU output, which is not written,
    /// but the split counters move as the Fortran's do.
    fn cutfem_cut_mesh_loops(&mut self) -> Result<(), SolverError> {
        for _sweep in 0..2 {
            self.cutfem_cut_mesh_sweep()?;
        }
        Ok(())
    }

    fn cutfem_cut_mesh_sweep(&mut self) -> Result<(), SolverError> {
        for t in 0..self.mesh.n_bulk() {
            let e = ElemId(t as u32);
            let (is_cut, _) = self.cutfem_check(e)?;
            if is_cut {
                loop {
                    let sub = self.cutfem_bulk(e)?;
                    if !sub.more {
                        break;
                    }
                }
            }
        }
        for t in 0..self.mesh.n_bulk() {
            let e = ElemId(t as u32);
            while let Some(ie) = self.cutfem_interface(e)? {
                if !ie.more {
                    break;
                }
            }
        }
        Ok(())
    }

    /// `LevelSetUpdate(Solver, Mesh)` for a levelset moved with a constant
    /// `Levelset Velocity`: the isoline mesh's nodes move by `V dt`, its
    /// segments become the polyline (`PopulatePolyline`, whose length is
    /// `res: cutfem isoline length`), and the levelset at each node becomes
    /// the signed distance to it (`SignedDistance`); outside the cut dofs
    /// the other isoline variables take the polyline's value at the closest
    /// point. Refused: the other ways to move the levelset, the node history,
    /// the normal move, `CutFEM check ls field` and a polyline that crosses
    /// itself.
    fn cutfem_levelset_update(&mut self, iso: &mut IsoMesh) -> Result<(), SolverError> {
        let p = self.params().clone();
        for k in [
            "levelset speed",
            "levelset velocity variable",
            "levelset speed variable",
            "cutfem normal move",
            "cutfem node history",
            "cutfem check ls field",
            "cutfem signed distance nonzero",
        ] {
            if p.get(k).is_some_and(|v| !matches!(v, elmer_core::Value::Logical(false))) {
                return Err(SolverError::Unsupported(format!("CutFEM levelset option {k}")));
            }
        }
        let dt = self.sim.vars.scalar("timestep size").unwrap_or(0.0);
        let mut moving = false;
        if let Some(vx) = self.creal(&p, "levelset velocity 1")? {
            if vx.abs() > f64::EPSILON {
                moving = true;
                iso.x.iter_mut().for_each(|x| *x += vx * dt);
            }
        }
        if let Some(vy) = self.creal(&p, "levelset velocity 2")? {
            if vy.abs() > f64::EPSILON {
                moving = true;
                iso.y.iter_mut().for_each(|y| *y += vy * dt);
            }
        }
        let run = self.cut()?;
        let levelset = run.levelset.clone();
        let id = self.sim.slots[self.slot].variable.ok_or_else(|| cut_err("solver has no variable"))?;
        let solver_name = self.sim.vars.by_id(id).name.clone();
        // The isoline variables in the polyline columns.
        let mut nvar = 0;
        let (mut i_avoid, mut i_solver) = (0, 0);
        let mut cols: Vec<usize> = Vec::new(); // field index per isoline variable (0 when avoided)
        for k in 1.. {
            let Some(v) = p.string(&format!("isoline variable {k}")) else { break };
            let v = elmer_core::normalize_key(v);
            if v == levelset {
                i_avoid = k;
            }
            if v == solver_name {
                i_solver = k;
            }
            let Some(f) = iso.fields.iter().position(|(n, _, _)| *n == v) else { break };
            nvar = k;
            cols.push(f);
        }
        let phi1d = iso.fields.iter().position(|(n, _, _)| *n == levelset).ok_or_else(|| {
            SolverError::fatal("LevelSetUpdate", format!("Levelset function (\"{levelset}\") needed in 1D mesh!"))
        })?;
        // PopulatePolyline.
        let (x, y) = (&iso.x, &iso.y);
        let ne = iso.elems.len();
        let mut skip = vec![false; ne];
        let mut node_to_elem = vec![[0usize; 2]; x.len()];
        let mut node_count = vec![0usize; x.len()];
        let mut lines: Vec<PolyLine> = Vec::new();
        let mut tot_len = 0.0;
        for phase in 0..2 {
            let mut m = 0;
            for i in 0..ne {
                let [i0, i1] = iso.elems[i];
                let (x0, y0, x1, y1) = (x[i0], y[i0], x[i1], y[i1]);
                let ss = sq(x0 - x1) + sq(y0 - y1);
                if ss < f64::EPSILON {
                    skip[i] = true;
                    continue;
                }
                m += 1;
                if phase == 0 {
                    node_count[i0] += 1;
                    node_count[i1] += 1;
                    for k in [i0, i1] {
                        if node_to_elem[k][0] == 0 {
                            node_to_elem[k][0] = m;
                        } else {
                            node_to_elem[k][1] = m;
                        }
                    }
                    continue;
                }
                // The element index is compared with the line counter, as the
                // Fortran does.
                let other = |k: usize| if node_to_elem[k][0] == i + 1 { node_to_elem[k][1] } else { node_to_elem[k][0] };
                let mut prev = other(i0) as i64;
                let mut next = other(i1) as i64;
                if node_count[i0] > 2 {
                    prev = -1;
                }
                if node_count[i1] > 2 {
                    next = -1;
                }
                tot_len += ss.sqrt();
                let mut vals = vec![x0, x1, y0, y1, iso.fields[phi1d].1[i0], iso.fields[phi1d].1[i1]];
                let (xmax0, xmin0, ymax0, ymin0) = (x0.max(x1), x0.min(x1), y0.max(y1), y0.min(y1));
                for j in i + 1..ne {
                    if skip[j] {
                        continue;
                    }
                    let [i2, i3] = iso.elems[j];
                    if i2 == i0 || i2 == i1 || i3 == i0 || i3 == i1 {
                        continue;
                    }
                    let (x2, y2, x3, y3) = (x[i2], y[i2], x[i3], y[i3]);
                    if x2.max(x3) < xmin0 || xmax0 < x2.min(x3) || y2.max(y3) < ymin0 || ymax0 < y2.min(y3) {
                        continue;
                    }
                    if segments_intersect([x2, y2], [x3, y3], [x0, y0], [x1, y1]) {
                        return Err(SolverError::Unsupported("a CutFEM isoline that crosses itself (polyline repair)".into()));
                    }
                }
                for (k, &f) in cols.iter().enumerate() {
                    if k + 1 == i_avoid {
                        continue;
                    }
                    let v = &iso.fields[f].1;
                    vals.push(v[i0]);
                    vals.push(v[i1]);
                }
                lines.push(PolyLine { vals, prev, next });
            }
        }
        self.sim.model.simulation.insert("res: cutfem isoline length", elmer_core::Value::real(tot_len));
        self.info("LevelSetUpdate", &format!("Cutfem isoline length:{}", crate::fortran_fmt::es(tot_len, 12, 3)), 6);
        if lines.is_empty() {
            return Err(SolverError::Unsupported("a CutFEM levelset without an isoline".into()));
        }
        // SignedDistance at each node of the levelset field.
        let nn = self.mesh.n_nodes();
        let st_perm = self.cut_state()?.cut_perm.clone();
        let ls_id = self.sim.vars.id(&levelset).ok_or_else(|| cut_err("no levelset variable"))?;
        let other_vars: Vec<(usize, Option<VarId>)> = (1..=nvar)
            .filter(|&k| k != i_avoid)
            .map(|k| {
                let name = elmer_core::normalize_key(p.string(&format!("isoline variable {k}")).unwrap_or(""));
                (k, self.sim.vars.id(&name))
            })
            .collect();
        for node in 0..nn {
            let Some(j) = self.sim.vars.by_id(ls_id).slot(node as u32) else { continue };
            let [xp, yp, _] = self.mesh.coord(node as u32);
            let (mut mindist2, mut sgn) = (f64::MAX, 1.0);
            let (mut smin, mut imin) = (0.0, 0);
            for (i, l) in lines.iter().enumerate() {
                let (x0, x1, y0, y1) = (l.vals[0], l.vals[1], l.vals[2], l.vals[3]);
                let a = xp - x0;
                let b = x0 - x1;
                let d = y0 - y1;
                let c = yp - y0;
                let mut s = -(a * b + c * d) / (sq(b) + sq(d));
                // MIN(MAX(s, 0), 1): a NaN becomes 0 as with gfortran, not
                // NaN as `clamp` would leave it.
                #[allow(clippy::manual_clamp)]
                {
                    s = s.max(0.0).min(1.0);
                }
                let xm = (1.0 - s) * x0 + s * x1;
                let ym = (1.0 - s) * y0 + s * y1;
                let dist2 = sq(xp - xm) + sq(yp - ym);
                if dist2 > mindist2 {
                    continue;
                }
                let dir1 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0);
                sgn = if dir1 < 0.0 { -1.0 } else { 1.0 };
                mindist2 = dist2;
                smin = s;
                imin = i;
            }
            let phip = sgn * mindist2.sqrt();
            if nvar > 0 && st_perm[node] == 0 {
                let mut ncol = 7;
                for &(k, vid) in &other_vars {
                    let l = &lines[imin];
                    let value = (1.0 - smin) * l.vals[ncol - 1] + smin * l.vals[ncol];
                    ncol += 2;
                    if k == i_solver {
                        // OrigMeshValues: the variable's own values while the
                        // cut ones are in place.
                        let run = self.sim.slots[self.slot].cutfem.as_deref_mut().expect("started");
                        let saved = run.saved.as_mut().ok_or_else(|| cut_err("no saved variable"))?;
                        let slot = saved.perm.as_ref().map_or(Some(node), |pm| pm.get(node).copied().filter(|&s| s != NO_DOF).map(|s| s as usize));
                        let Some(sj) = slot else {
                            return Err(SolverError::fatal("SignedDistance", "We should maybe not be here?"));
                        };
                        saved.values[sj] = value;
                    } else {
                        let vid = vid.ok_or_else(|| cut_err("isoline variable not on the mesh"))?;
                        let v = self.sim.vars.by_id_mut(vid);
                        let Some(sj) = v.slot(node as u32) else {
                            return Err(SolverError::fatal("SignedDistance", "We should maybe not be here?"));
                        };
                        if v.dofs != 1 {
                            return Err(SolverError::Unsupported("CutFEM isoline variables of more than one dof".into()));
                        }
                        v.values[sj] = value;
                    }
                }
            }
            if moving {
                self.sim.vars.by_id_mut(ls_id).values[j] = phip;
            }
        }
        Ok(())
    }

    /// `DefaultFinish`'s CutFEM part: `LevelSetUpdate` (the isoline mesh
    /// becomes `Mesh % Next`), the cut mesh with `CutFEM Mesh Save` (made
    /// after the levelset moved, so it carries the new levelset on the
    /// step's cuts), then `CutFEMVariableFinalize`.
    pub(crate) fn cutfem_finish(&mut self) -> Result<(), SolverError> {
        let mut iso = self.cutfem_isomesh()?;
        self.cutfem_levelset_update(&mut iso)?;
        let mut extra = vec![ExtraMesh { suffix: "-isomesh".into(), n_nodes: iso.x.len(), fields: std::mem::take(&mut iso.fields) }];
        if self.params().logical("cutfem mesh save") == Some(true) {
            self.cutfem_cut_mesh_loops()?;
            let perm = self.cut_state()?.cut_perm.clone();
            let n = perm.iter().copied().max().unwrap_or(0) as usize;
            let fields = self.cut_mesh_fields(&perm, "project variable")?;
            extra.push(ExtraMesh { suffix: "-cutfem".into(), n_nodes: n, fields });
        }
        self.sim.extra_meshes = extra;
        self.cutfem_variable_finalize()
    }

    /// `CutFEMVariableFinalize`: the cut values at the shared nodes go back
    /// to the original variable; with `CutFEM extend` the extension solve
    /// runs on a matrix nothing is assembled into (`LocalFitMatrix` keeps
    /// its element matrix), so its right-hand side is zero, the solution
    /// "trivially zero", and the extended nodes become 0. Then the original
    /// variable, permutation, history and system are back.
    fn cutfem_variable_finalize(&mut self) -> Result<(), SolverError> {
        let id = self.sim.slots[self.slot].variable.ok_or_else(|| cut_err("solver has no variable"))?;
        let sim = &mut *self.sim;
        let slot = &mut sim.slots[self.slot];
        let run = slot.cutfem.as_deref_mut().ok_or_else(|| cut_err("CutFEM has not started"))?;
        let st = run.state.take().ok_or_else(|| cut_err("no CutFEM permutation"))?;
        let saved = run.saved.take().ok_or_else(|| cut_err("no saved variable"))?;
        let var = sim.vars.by_id_mut(id);
        let cut_values = std::mem::replace(&mut var.values, saved.values);
        var.perm = saved.perm;
        var.prev_values = saved.prev;
        slot.perm = var.perm.clone();
        let dofs = run.dofs;
        for i in 0..st.nn {
            let j = st.cut_perm[i];
            let Some(k) = var.slot(i as u32) else { continue };
            if j == 0 {
                continue;
            }
            for l in 0..dofs {
                var.values[dofs * k + l] = cut_values[dofs * (j as usize - 1) + l];
            }
        }
        if let Some(ext) = &st.extend_perm {
            // The element loop of the fit only moves the split counters.
            for t in 0..run.mesh.nodes.len() {
                if run.mesh.nodes[t].iter().any(|&k| run.phi[k as usize].is_none()) {
                    continue;
                }
                loop {
                    let sub = run.bulk.next(&st, &run.mesh, &run.edges, t, &run.phi)?;
                    if !(sub.cut && sub.more) {
                        break;
                    }
                }
            }
            for i in 0..st.nn {
                if ext[i] == 0 {
                    continue;
                }
                if let Some(k) = var.slot(i as u32) {
                    var.values[k] = 0.0;
                }
            }
        }
        slot.system = saved.system;
        run.coords = Some((st.x, st.y, st.z));
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles of the unit square, split by the circle around the
    /// origin of radius 0.7: edges in first-use order with the first
    /// element's orientation, and the cut edges' interpolation.
    fn square() -> (CutMesh, Vec<f64>, Vec<f64>) {
        let x = vec![0.0, 1.0, 1.0, 0.0];
        let y = vec![0.0, 0.0, 1.0, 1.0];
        let mesh = CutMesh { nn: 4, codes: vec![303, 303], nodes: vec![vec![0, 1, 2], vec![0, 2, 3]], body: vec![1, 1] };
        (mesh, x, y)
    }

    #[test]
    fn edges_keep_the_first_elements_orientation() {
        let (mesh, _, _) = square();
        let e = Edges::find(&mesh);
        assert_eq!(e.nodes, vec![[0, 1], [1, 2], [2, 0], [2, 3], [3, 0]]);
        assert_eq!(e.of_element, vec![vec![0, 1, 2], vec![2, 3, 4]]);
    }

    #[test]
    fn a_circle_cuts_two_edges_of_each_triangle() {
        let (mesh, x, y) = square();
        let edges = Edges::find(&mesh);
        let phi: Vec<Option<f64>> = (0..4).map(|i| Some(x[i] * x[i] + y[i] * y[i] - 0.49)).collect();
        let opts = PermOptions {
            body_in: 1,
            body_out: 0,
            eps: 1e-3,
            abs_eps: false,
            extend: false,
            extend_width: None,
            passive_outside: None,
            passive_inside: None,
        };
        let z = vec![0.0; 4];
        let mut bulk = BulkSplit::default();
        let st = create_perm(&mesh, &edges, &phi, &opts, (x, y, z), &mut bulk, |_| false).unwrap();
        // Edges 0-1, 2-0 and 3-0 cross the circle; 1-2 and 2-3 do not.
        assert_eq!(st.cut_dof[4..], [true, false, true, false, true]);
        // Edge 0 -> 1: h1 = -0.49, h2 = 0.51, r = 0.51; the node at x = 0.49.
        assert!((st.x[4] - 0.49).abs() < 1e-15);
        // Passive outside: node 0 (inside) and the cut edges are the dofs.
        assert_eq!(st.cut_perm, vec![1, 0, 0, 0, 2, 0, 3, 0, 4]);
        // Both are interface elements: body_cut = max(body_in, body_out) = 1.
        assert_eq!(st.body, vec![1, 1]);
        // The first triangle splits into three pieces, the inside one first.
        let p1 = bulk.next(&st, &mesh, &edges, 0, &phi).unwrap();
        assert!(p1.cut && p1.more);
        let mut pieces = vec![p1];
        while pieces.last().unwrap().more {
            pieces.push(bulk.next(&st, &mesh, &edges, 0, &phi).unwrap());
        }
        assert_eq!(pieces.len(), 3);
        // One interface segment per triangle.
        let mut bc = BcSplit::default();
        let s = bc.next(&st, &mesh, &edges, 0, &phi).unwrap().unwrap();
        assert!(!s.more);
        let mut v = s.nodes.to_vec();
        v.sort();
        assert_eq!(v, vec![4, 6]);
    }
}
