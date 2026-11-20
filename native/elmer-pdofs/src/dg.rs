//! The degrees of freedom of a `Discontinuous Galerkin` solver with `DG
//! Reduced Basis`: one value per element node, merged by node within groups
//! of bodies, so the field is continuous inside a group and discontinuous
//! between groups (heat across a `Heat Gap`, say).
//!
//! Ported: the DG indexes of the mesh (`NonNodalElements`, MeshLoad.F90),
//! the reduced numbering of `InitialPermutation` (Lists.F90:271-345, `DG
//! Reduced Basis Mapping` / `Master Bodies` / `Slave Bodies`), and the
//! matrix pattern of `MakeListMatrix` (ElementUtils.F90:466-640): the bulk
//! elements, the edges (2D) or faces (3D) between bodies, and with `DG
//! Indirect Connections` the row structures of the dofs paired across those
//! interfaces exchanged (`List_ExchangeRowStructure`). Elmer does not
//! optimise the bandwidth of a reduced basis.
//!
//! Without `DG Reduced Basis` every DG index is its own dof, numbered over
//! the edges (faces) in mesh order, the left then the right parent's
//! indexes at first sight (InitialPermutation's DG branch), and the
//! pattern couples both parents of every edge. `Optimize Bandwidth`
//! (default on) then renumbers it with Elmer's reverse Cuthill-McKee
//! (`OptimizeBandwidth`), kept unless it widens the band (or with
//! `Optimize Bandwidth Use Always`).

use elmer_assembly::{bandwidth, bandwidth_reordered, initialize_matrix, rcm_order, ListMatrix};
use elmer_core::{CsrMatrix, ElemId, ElementBasis, MeshView, ParamMap, Params, QuadPoint, Quadrature, Value, NO_DOF};
use elmer_solver::{DgIndexes, Model, Perm, SolverDofs, SolverError};

/// A DG solver's layout.
pub struct DgDofs {
    ix: DgIndexes,
    /// DG index -> slot.
    perm: Vec<u32>,
    n_slots: usize,
    list: ListMatrix,
    /// The bandwidth optimisation's `new_index` of the list rows.
    map: Option<Vec<u32>>,
    /// The element nodes and parents, for the slots of boundary elements.
    nodes: Vec<Vec<u32>>,
    parents: Vec<[Option<ElemId>; 2]>,
    etypes: Vec<elmer_core::ElementType>,
    n_bulk: usize,
    dim: usize,
}

fn int_list(p: &ParamMap, key: &str) -> Option<Vec<i32>> {
    match p.get(key) {
        Some(Value::Integer(v)) => Some(v.clone()),
        _ => None,
    }
}

impl DgDofs {
    /// `CreateMatrix` of a DG solver whose active bulk elements are
    /// `active`.
    pub fn build(mesh: &dyn MeshView, model: &Model, solver: &ParamMap, active: &[ElemId]) -> Result<Self, SolverError> {
        let reduced = solver.logical("dg reduced basis") == Some(true);
        let ix = DgIndexes::for_dg_solver(mesh, model)?;
        let nb = mesh.n_bulk();
        let mut is_active = vec![false; nb];
        for e in active {
            if (e.0 as usize) < nb {
                is_active[e.0 as usize] = true;
            }
        }
        let body = |t: usize| model.body_of_tag(mesh.element(ElemId(t as u32)).tag).unwrap_or(0) as i32;

        // InitialPermutation with DG Reduced Basis.
        let map = int_list(solver, "dg reduced basis mapping");
        let master = int_list(solver, "dg reduced basis master bodies");
        let slave = int_list(solver, "dg reduced basis slave bodies");
        let max_group = match (&map, &master) {
            (Some(m), _) => {
                if m.len() != model.bodies.len() {
                    return Err(SolverError::fatal("InitialPermutation", "Invalid size of > Dg Reduced Basis Mapping <"));
                }
                m.iter().copied().max().unwrap_or(0)
            }
            (None, Some(_)) => 2,
            (None, None) => model.bodies.len() as i32,
        };
        let mut perm = vec![0u32; ix.count];
        let mut node_index = vec![0u32; mesh.n_nodes()];
        let mut k = 0u32;
        let edges = elmer_solver::dg_edges(mesh)?;
        if !reduced {
            for (l, r, _) in &edges {
                for p in std::iter::once(*l).chain(*r) {
                    let t = p.0 as usize;
                    if !is_active[t] {
                        continue;
                    }
                    for &d in &ix.index[t] {
                        if perm[d as usize] == 0 {
                            k += 1;
                            perm[d as usize] = k;
                        }
                    }
                }
            }
        }
        for group0 in (1..=max_group).filter(|_| reduced) {
            match &slave {
                Some(s) if group0 == 2 => {
                    for t in 0..nb {
                        if s.contains(&body(t)) {
                            for &n in mesh.element(ElemId(t as u32)).nodes {
                                node_index[n as usize] = 0;
                            }
                        }
                    }
                }
                _ => node_index.iter_mut().for_each(|x| *x = 0),
            }
            for t in 0..nb {
                let group = body(t);
                let in_group = match (&master, &map) {
                    (Some(m), _) => {
                        if group0 == 1 {
                            m.contains(&group)
                        } else {
                            !m.contains(&group)
                        }
                    }
                    (None, Some(m)) => group >= 1 && m.get(group as usize - 1) == Some(&group0),
                    (None, None) => group0 == group,
                };
                if !in_group || !is_active[t] {
                    continue;
                }
                for (i, &n) in mesh.element(ElemId(t as u32)).nodes.iter().enumerate() {
                    if node_index[n as usize] == 0 {
                        k += 1;
                        node_index[n as usize] = k;
                    }
                    perm[ix.index[t][i] as usize] = node_index[n as usize];
                }
            }
        }
        let n_slots = k as usize;
        let slot = |d: u32| -> Option<u32> {
            let p = perm[d as usize];
            (p > 0).then(|| p - 1)
        };

        // MakeListMatrix (DG Reduced Basis).
        let mut list = ListMatrix::new(n_slots);
        let add = |list: &mut ListMatrix, dg: &[u32]| {
            let mut ind: Vec<u32> = dg.iter().filter_map(|&d| slot(d)).collect();
            ind.sort_unstable();
            for &r in &ind {
                list.add_sorted(r as usize, &ind);
            }
        };
        for t in (0..nb).filter(|&t| is_active[t]) {
            add(&mut list, &ix.index[t]);
        }
        if !reduced {
            for (l, r, _) in &edges {
                let mut both: Vec<u32> = Vec::new();
                for p in std::iter::once(*l).chain(*r) {
                    if is_active[p.0 as usize] {
                        both.extend_from_slice(&ix.index[p.0 as usize]);
                    }
                }
                add(&mut list, &both);
            }
        }
        let indirect = reduced && solver.logical("dg indirect connections") == Some(true);
        let mut pairs = vec![u32::MAX; n_slots];
        for (l, r, _) in elmer_solver::dg_interfaces(mesh)?.into_iter().filter(|_| reduced) {
            let (lt, rt) = (l.0 as usize, r.0 as usize);
            if !is_active[lt] || !is_active[rt] || body(lt) == body(rt) {
                continue;
            }
            let (ln, rn) = (mesh.element(l).nodes, mesh.element(r).nodes);
            if indirect {
                for (i, &n) in ln.iter().enumerate() {
                    let Some(k1) = slot(ix.index[lt][i]) else { continue };
                    if let Some(j) = rn.iter().position(|&m| m == n) {
                        if let Some(k2) = slot(ix.index[rt][j]) {
                            if k1 != k2 {
                                pairs[k1 as usize] = k2;
                            }
                        }
                    }
                }
            }
            let both: Vec<u32> = ix.index[lt].iter().chain(&ix.index[rt]).copied().collect();
            add(&mut list, &both);
        }
        // OptimizeBandwidth (CreateMatrix suppresses it for a reduced basis).
        let mut map = None;
        if !reduced && n_slots > 0 && solver.logical("optimize bandwidth") != Some(false) {
            let initial = bandwidth(&list) + 1;
            let new_index = rcm_order(&list);
            let optimized = bandwidth_reordered(&list, &new_index) + 1;
            let always = solver.logical("optimize bandwidth use always") == Some(true);
            if !(initial < optimized && !always) {
                for p in perm.iter_mut().filter(|p| **p > 0) {
                    *p = new_index[*p as usize - 1] + 1;
                }
                map = Some(new_index);
            }
        }
        if indirect {
            for (k1, &k2) in pairs.iter().enumerate() {
                if k2 == u32::MAX {
                    continue;
                }
                let r1 = list.row(k1).to_vec();
                let r2 = list.row(k2 as usize).to_vec();
                list.add_sorted(k2 as usize, &r1);
                list.add_sorted(k1, &r2);
            }
        }

        // AddEquationBasics (MainUtils.F90:1703): the Perm is allocated with
        // MAX(nodal size, MaxDGDOFs * (bulk + boundary elements)) entries.
        let max_dg = (0..nb).map(|t| mesh.element(ElemId(t as u32)).nodes.len()).max().unwrap_or(0);
        let size = mesh.n_nodes().max(max_dg * (nb + mesh.n_boundary())).max(ix.count);
        let mut perm: Vec<u32> = perm.iter().map(|&p| if p == 0 { NO_DOF } else { p - 1 }).collect();
        perm.resize(size, NO_DOF);
        let nodes = (0..nb + mesh.n_boundary()).map(|t| mesh.element(ElemId(t as u32)).nodes.to_vec()).collect();
        let parents = (0..nb + mesh.n_boundary()).map(|t| mesh.element(ElemId(t as u32)).parents).collect();
        let etypes = (0..nb + mesh.n_boundary()).map(|t| mesh.element(ElemId(t as u32)).etype).collect();
        Ok(DgDofs { ix, perm, n_slots, list, map, nodes, parents, etypes, n_bulk: nb, dim: mesh.dim() })
    }

    fn push_slots(&self, t: usize, out: &mut Vec<u32>) {
        for &d in &self.ix.index[t] {
            out.push(if d == NO_DOF { NO_DOF } else { self.perm[d as usize] });
        }
    }
}

impl SolverDofs for DgDofs {
    fn perm(&self) -> Perm {
        Perm { perm: self.perm.clone(), n_slots: self.n_slots }
    }

    fn create_matrix(&self, dofs: usize) -> CsrMatrix<f64> {
        initialize_matrix(&self.list, self.map.as_deref(), dofs)
    }

    /// `mGetElementDOFs` of a DG solver: a bulk element's DG indexes; a
    /// boundary element (no `DGDOFs` of its own) has its parents' (left,
    /// then right).
    fn element_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        let t = e.0 as usize;
        if t < self.n_bulk {
            self.push_slots(t, out);
        } else {
            for p in self.parents[t].iter().flatten() {
                self.push_slots(p.0 as usize, out);
            }
        }
    }

    fn boundary_extra_slots(&self, _e: ElemId, out: &mut Vec<u32>) {
        out.clear();
    }

    fn element_basis(&self, _e: ElemId) -> Result<Box<dyn ElementBasis + Send + Sync>, String> {
        Ok(Box::new(elmer_fe::LagrangeBasis::new(self.dim)))
    }

    /// `GaussPoints(Element)` of the nodal element.
    fn element_rule(&self, e: ElemId) -> Result<Vec<QuadPoint>, String> {
        Ok(elmer_fe::ElmerQuadrature.gauss_points(self.etypes[e.0 as usize], None))
    }

    fn element_nodes(&self, e: ElemId) -> Vec<u32> {
        self.nodes[e.0 as usize].clone()
    }

    fn local_bubbles(&self, _e: ElemId) -> usize {
        0
    }

    fn dg_layout(&self) -> Option<&DgIndexes> {
        Some(&self.ix)
    }
}
