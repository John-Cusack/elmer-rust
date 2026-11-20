//! Elemental (`-elem`) and integration-point (`-ip`) exported variables and
//! `UpdateExportedVariables` (SolverBasics.F90:6602): an exported variable
//! is filled from the body forces (and boundary conditions) that give a
//! keyword of its name, at the nodes of each element (nodal), as their mean
//! (elemental), or evaluated at each Gauss point with its dependencies
//! interpolated there (integration points, `ListGetElementReal`).

use elmer_core::{ElemId, FieldKind, Logger, MeshView, ParamMap, Params, Value, NO_DOF};
use elmer_expr::FieldLookup;

use crate::error::SolverError;
use crate::expr::SimFields;
use crate::fortran_fmt::i2s;
use crate::sim::{Simulation, SolverSlot};

const CALLER: &str = "UpdateExportedVariables";

/// `Exported Variable n Mask = "[section:] keyword"`: the section (default
/// `body force`) and the logical keyword that select the elements.
pub(crate) fn mask_of(params: &ParamMap, key: &str) -> Option<(String, String)> {
    let text = params.string(&format!("{key} mask"))?;
    let text = elmer_core::normalize_key(text);
    Some(match text.split_once(':') {
        Some((sec, mask)) if !sec.is_empty() => (sec.trim().to_string(), mask.trim().to_string()),
        _ => ("body force".to_string(), text),
    })
}

/// `ListGetSection(Element, name)`: the section of that kind applying to
/// element `e`, through its body (or its boundary condition).
fn section_of<'m>(sim: &'m Simulation, mesh: &dyn MeshView, e: ElemId, name: &str) -> Result<Option<&'m ParamMap>, SolverError> {
    let m = &sim.model;
    let tag = mesh.element(e).tag;
    let is_bulk = (e.0 as usize) < mesh.n_bulk();
    let body = if is_bulk { m.body_of_tag(tag).and_then(|b| m.bodies.get(b - 1)) } else { None };
    let linked = |sections: &'m [ParamMap], key: &str| {
        body.and_then(|b| b.integer(key)).filter(|&k| k >= 1).and_then(|k| sections.get(k as usize - 1))
    };
    Ok(match name {
        "body" => body,
        "material" => linked(&m.materials, "material"),
        "body force" => linked(&m.body_forces, "body force"),
        "initial condition" => linked(&m.initial_conditions, "initial condition"),
        "equation" => linked(&m.equations, "equation"),
        "boundary condition" => (!is_bulk).then(|| m.bc_of_tag(tag)).flatten().map(|i| &m.boundary_conditions[i - 1]),
        other => return Err(SolverError::fatal("ListGetSection", format!("Unknown section name: {other}"))),
    })
}

/// `ListGetLogicalGen(section, mask)` (Lists.F90): false without the
/// section or key; a logical is itself; a constant real scalar is true when
/// positive; any other present value (an integer, an array, a function of
/// variables, a string) is true. A constant MATC expression or procedure
/// would be evaluated by `ListGetConstReal` and is refused.
fn mask_on(section: Option<&ParamMap>, mask: &str) -> Result<bool, SolverError> {
    Ok(match section.and_then(|s| s.get(mask)) {
        None => false,
        Some(Value::Logical(b)) => *b,
        Some(Value::Real { rows, cols, data }) if rows * cols == 1 => data.first().is_some_and(|&x| x > 0.0),
        Some(Value::Dependent(d)) if d.variables.is_empty() => {
            return Err(SolverError::Unsupported(format!(
                "mask > {mask} < given by a constant expression or procedure (ListGetLogicalGen)"
            )))
        }
        Some(_) => true,
    })
}

/// `CreateElementsPerm` with a mask: active elements whose section sets the
/// mask, numbered in element order. As in the Fortran, an element whose
/// equation is not active (a boundary element, say) keeps the previous
/// element's verdict.
pub(crate) fn element_perm_masked(
    sim: &Simulation,
    mesh: &dyn MeshView,
    slot: &SolverSlot,
    (sec, mask): (&str, &str),
) -> Result<(Vec<u32>, usize), SolverError> {
    let n = mesh.n_bulk() + mesh.n_boundary();
    let mut active = vec![false; n];
    for e in &slot.active {
        active[e.0 as usize] = true;
    }
    let mut perm = vec![NO_DOF; n];
    let (mut on, mut m) = (false, 0usize);
    for t in 0..n {
        if active[t] {
            on = mask_on(section_of(sim, mesh, ElemId(t as u32), sec)?, mask)?;
        }
        if on {
            perm[t] = m as u32;
            m += 1;
        }
    }
    Ok((perm, m))
}

/// `Element % DGIndexes` of every bulk and boundary element, 0-based.
#[derive(Clone, Debug)]
pub struct DgIndexes {
    /// Per element, one index per node ([`NO_DOF`] where a boundary element
    /// found none in its parents).
    pub index: Vec<Vec<u32>>,
    /// The number of DG indexes (`MAXVAL` over the bulk elements).
    pub count: usize,
}

impl DgIndexes {
    /// `CreateDGPerm` when no element has DG indexes yet: the nodes of the
    /// bulk elements numbered in element order, then each boundary element
    /// takes its nodes' indexes from its first parent (left, else right).
    pub fn create(mesh: &dyn MeshView) -> Self {
        let nb = mesh.n_bulk();
        let mut index: Vec<Vec<u32>> = Vec::with_capacity(nb + mesh.n_boundary());
        let mut next = 0u32;
        for t in 0..nb {
            let n = mesh.element(ElemId(t as u32)).nodes.len();
            index.push((next..next + n as u32).collect());
            next += n as u32;
        }
        for t in nb..nb + mesh.n_boundary() {
            let el = mesh.element(ElemId(t as u32));
            let mut idx = vec![NO_DOF; el.nodes.len()];
            if let Some(p) = el.parents.iter().flatten().find(|p| (p.0 as usize) < nb) {
                let parent = mesh.element(*p);
                for (i, n) in el.nodes.iter().enumerate() {
                    if let Some(j) = parent.nodes.iter().position(|m| m == n) {
                        idx[i] = index[p.0 as usize][j];
                    }
                }
            }
            index.push(idx);
        }
        DgIndexes { index, count: next as usize }
    }

    /// The DG indexes `NonNodalElements` (MeshLoad.F90:986-1190) makes when
    /// a solver is DG: the bulk elements' nodes in element order; a
    /// boundary element gets its parent's indexes only when some material
    /// has `DG Parent Material` (`InheritDG`), from the first parent (left,
    /// then right) whose material sets it true, or its only parent. Other
    /// boundary elements have none ([`NO_DOF`]).
    pub fn for_dg_solver(mesh: &dyn MeshView, model: &crate::Model) -> Result<Self, SolverError> {
        let mut ix = DgIndexes::create(mesh);
        let nb = mesh.n_bulk();
        let inherit = model.materials.iter().any(|m| m.get("dg parent material").is_some());
        let parent_flag = |p: ElemId| -> bool {
            model
                .body_of_tag(mesh.element(p).tag)
                .and_then(|b| model.bodies.get(b - 1))
                .and_then(|b| b.integer("material"))
                .filter(|&m| m >= 1)
                .and_then(|m| model.materials.get(m as usize - 1))
                .and_then(|m| m.logical("dg parent material"))
                .unwrap_or(false)
        };
        for t in nb..nb + mesh.n_boundary() {
            let el = mesh.element(ElemId(t as u32));
            let mut idx = vec![NO_DOF; el.nodes.len()];
            if inherit {
                let parents: Vec<ElemId> = el.parents.iter().flatten().copied().collect();
                let parent = match parents.as_slice() {
                    [] => return Err(SolverError::fatal("NonnodalElements", "Cannot define DG indexes for BC!")),
                    [p] => *p,
                    _ => *parents.iter().find(|&&p| parent_flag(p)).ok_or_else(|| {
                        SolverError::fatal("NonnodalElements", "Cannot define DG indexes for internal BC!")
                    })?,
                };
                let pnodes = mesh.element(parent).nodes;
                for (i, n) in el.nodes.iter().enumerate() {
                    if let Some(j) = pnodes.iter().position(|m| m == n) {
                        idx[i] = ix.index[parent.0 as usize][j];
                    }
                }
            }
            ix.index[t] = idx;
        }
        Ok(ix)
    }
}

/// `CreateDGPerm(Solver, DGPerm, DGCount[, MaskName, SecName])`: DG index ->
/// 0-based slot over the DG indexes of the solver's active bulk elements
/// (those whose section sets the mask, when there is one), numbered in DG
/// index order. A boundary element's indexes are its parent's, so they add
/// none. The mesh's DG indexes are created on first use.
pub(crate) fn dg_perm(
    sim: &mut Simulation,
    mesh: &dyn MeshView,
    slot: &SolverSlot,
    mask: Option<(&str, &str)>,
) -> Result<(Vec<u32>, usize), SolverError> {
    if sim.dg_indexes.is_none() {
        sim.dg_indexes = Some(DgIndexes::create(mesh));
    }
    let dg = sim.dg_indexes.as_ref().expect("created above");
    let mut perm = vec![NO_DOF; dg.count];
    for &e in &slot.active {
        if (e.0 as usize) >= mesh.n_bulk() {
            continue;
        }
        if let Some((sec, key)) = mask {
            if !mask_on(section_of(sim, mesh, e, sec)?, key)? {
                continue;
            }
        }
        for &i in &dg.index[e.0 as usize] {
            perm[i as usize] = 0;
        }
    }
    let mut n = 0u32;
    for p in perm.iter_mut().filter(|p| **p != NO_DOF) {
        *p = n;
        n += 1;
    }
    Ok((perm, n as usize))
}

/// Local edges (1-based end nodes) of a 2D element, as `FindMeshEdges2D`
/// visits them.
fn edges_2d(code: u16) -> Option<&'static [[usize; 2]]> {
    match code / 100 {
        3 => Some(&[[1, 2], [2, 3], [3, 1]]),
        4 => Some(&[[1, 2], [2, 3], [3, 4], [4, 1]]),
        _ => None,
    }
}

/// Local faces (1-based nodes) of a 3D element, as `FindMeshFaces3D`
/// visits them.
fn faces_3d(code: u16) -> Option<&'static [&'static [usize]]> {
    match code / 100 {
        5 => Some(&[&[1, 2, 3], &[1, 2, 4], &[2, 3, 4], &[3, 1, 4]]),
        7 => Some(&[&[1, 2, 3], &[4, 5, 6], &[1, 2, 5, 4], &[3, 2, 5, 6], &[3, 1, 4, 6]]),
        8 => Some(&[&[1, 2, 3, 4], &[5, 6, 7, 8], &[1, 2, 6, 5], &[2, 3, 7, 6], &[3, 4, 8, 7], &[4, 1, 5, 8]]),
        _ => None,
    }
}

/// The interfaces of the bulk elements in first-visit order (`Mesh %
/// Edges` in 2D, `Mesh % Faces` in 3D): `(left, right)` element pairs for
/// those with two bulk parents. Left is the first element visiting it.
pub fn dg_interfaces(mesh: &dyn MeshView) -> Result<Vec<(ElemId, ElemId, Vec<u32>)>, SolverError> {
    Ok(dg_edges(mesh)?.into_iter().filter_map(|(l, r, n)| r.map(|r| (l, r, n))).collect())
}

/// Every edge (2D) or face (3D) of the bulk elements in first-visit order
/// (`Mesh % Edges` / `Mesh % Faces`), with its left (first visiting) and,
/// for an interior one, right parent.
#[allow(clippy::type_complexity)]
pub fn dg_edges(mesh: &dyn MeshView) -> Result<Vec<(ElemId, Option<ElemId>, Vec<u32>)>, SolverError> {
    let dim = mesh.dim();
    let mut seen: std::collections::HashMap<Vec<u32>, usize> = std::collections::HashMap::new();
    let mut out: Vec<(ElemId, Option<ElemId>, Vec<u32>)> = Vec::new();
    for t in 0..mesh.n_bulk() {
        let el = mesh.element(ElemId(t as u32));
        let code = el.etype.0;
        let locals: Vec<Vec<usize>> = if dim == 2 {
            edges_2d(code)
                .ok_or_else(|| SolverError::Unsupported(format!("DG interfaces of element type {code}")))?
                .iter()
                .map(|e| e.to_vec())
                .collect()
        } else {
            faces_3d(code)
                .ok_or_else(|| SolverError::Unsupported(format!("DG interfaces of element type {code}")))?
                .iter()
                .map(|f| f.to_vec())
                .collect()
        };
        for local in locals {
            let nodes: Vec<u32> = local.iter().map(|&k| el.nodes[k - 1]).collect();
            let mut key = nodes.clone();
            key.sort_unstable();
            match seen.get(&key) {
                Some(&i) => {
                    if out[i].1.is_none() {
                        out[i].1 = Some(ElemId(t as u32));
                    }
                }
                None => {
                    seen.insert(key, out.len());
                    out.push((ElemId(t as u32), None, nodes));
                }
            }
        }
    }
    if let Some(keys) = mesh.edge_order().filter(|k| k.len() == out.len()) {
        let mut ix: Vec<usize> = (0..out.len()).collect();
        ix.sort_by_key(|&i| keys[i]);
        let mut slots: Vec<Option<(ElemId, Option<ElemId>, Vec<u32>)>> = out.into_iter().map(Some).collect();
        out = ix.into_iter().map(|i| slots[i].take().expect("a permutation")).collect();
    }
    Ok(out)
}

/// `InvActiveElements`: element -> 0-based position among the solver's
/// active elements, over every bulk and boundary element.
pub(crate) fn element_perm(mesh: &dyn MeshView, slot: &SolverSlot) -> (Vec<u32>, usize) {
    let mut perm = vec![NO_DOF; mesh.n_bulk() + mesh.n_boundary()];
    for (i, e) in slot.active.iter().enumerate() {
        perm[e.0 as usize] = i as u32;
    }
    (perm, slot.active.len())
}

/// `CreateIpPerm`: the offset of each element's Gauss points (`Gauss Points
/// on Ip Variables`, or the solver's rule of an active nodal element as
/// `GaussPointsAdapt` gives it), one entry more than the bulk and boundary
/// elements.
pub(crate) fn ip_perm(
    sim: &Simulation,
    mesh: &dyn MeshView,
    slot: &SolverSlot,
    mask: Option<(&str, &str)>,
) -> Result<Vec<u32>, SolverError> {
    let fe = sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature for IP variables (elmer-fe)".into()))?;
    let n_ip = slot.params.integer("gauss points on ip variables").filter(|&n| n > 0);
    // Without it, `GaussPointsAdapt(Element)` of the solver: `Number of
    // Integration Points`, else `Bulk Relative Integration Order` over
    // `Relative Integration Order`.
    for (key, what) in [
        ("element relative integration order", "Element Relative Integration Order"),
        ("adaptive integration variable", "adaptive integration"),
    ] {
        if n_ip.is_none() && slot.params.get(key).is_some() {
            return Err(SolverError::Unsupported(format!("IP variables with {what}")));
        }
    }
    let np = slot.params.integer("number of integration points").filter(|&n| n > 0).map(|n| n as usize);
    // `PDefs % Serendipity` of p elements.
    let serendipity = sim.model.simulation.logical("serendipity p elements").unwrap_or(true);
    let family_np = match slot.params.string("element integration points") {
        Some(def) => parse_elemental_gauss_rules(def)?,
        None => [0; 9],
    };
    let rel = match slot.params.get("bulk relative integration order") {
        Some(_) => slot.params.integer("bulk relative integration order"),
        None => slot.params.integer("relative integration order"),
    }
    .unwrap_or(0);
    let n = mesh.n_bulk() + mesh.n_boundary();
    let mut active = vec![false; n];
    for e in &slot.active {
        active[e.0 as usize] = true;
    }
    let mut offsets = vec![0u32; n + 1];
    let mut count = 0u32;
    for t in 0..n {
        let on = active[t]
            && match mask {
                Some((sec, key)) => mask_on(section_of(sim, mesh, ElemId(t as u32), sec)?, key)?,
                None => true,
            };
        if on {
            let etype = mesh.element(ElemId(t as u32)).etype;
            // `ElementalGaussNp` over the base count.
            let np = match family_np.get(etype.family() as usize) {
                Some(&k) if k > 0 => Some(k),
                _ => np,
            };
            if let (None, Some(d)) = (n_ip, slot.dofs.as_deref()) {
                count += p_rule_len(d, ElemId(t as u32), etype, np, rel, serendipity)? as u32;
                offsets[t + 1] = count;
                continue;
            }
            count += match (n_ip, np) {
                (Some(k), _) => k as u32,
                (None, Some(_)) => fe.quadrature.gauss_points(etype, np).len() as u32,
                (None, None) => fe
                    .quadrature
                    .gauss_points_rel(etype, rel)
                    .ok_or_else(|| {
                        SolverError::fatal(
                            "GaussPoints",
                            format!("Relative Integration Order = {rel} but a non-p element only has rules for {{-1,0,1}}"),
                        )
                    })?
                    .len() as u32,
            };
        }
        offsets[t + 1] = count;
    }
    Ok(offsets)
}

/// The point count `GaussPointsAdapt` gives p element `e` of the solver's
/// dofs: `np` points, else its own rule, offset by `GaussPoints(Element,
/// RelOrder = rel)`: `(NINT(n ** (1/dim)) + rel) ** dim` points of the
/// tensor families, a serendipity quadrilateral's then mapped to its
/// economic count.
fn p_rule_len(
    d: &dyn crate::SolverDofs,
    e: ElemId,
    etype: elmer_core::ElementType,
    np: Option<usize>,
    rel: i32,
    serendipity: bool,
) -> Result<usize, SolverError> {
    let gauss = |m: String| SolverError::fatal("GaussPoints", m);
    if let Some(np) = np {
        return Ok(d.element_rule_np(e, np).map_err(gauss)?.len());
    }
    let own = d.element_rule(e).map_err(gauss)?.len();
    if rel == 0 {
        return Ok(own);
    }
    let dim = match etype.family() {
        2 => 1,
        4 => 2,
        8 => 3,
        f => {
            return Err(SolverError::Unsupported(format!(
                "IP variables with Relative Integration Order on p elements of family {f}"
            )))
        }
    };
    let p1d = (own as f64).powf(1.0 / dim as f64).round() as i32 + rel;
    if p1d < 1 {
        return Err(SolverError::fatal("GaussPoints", "Number of integration points must remain positive!"));
    }
    let mut n = (p1d as usize).pow(dim);
    if serendipity && etype.family() == 4 {
        n = match n {
            9 => 8,
            16 => 12,
            25 => 20,
            36 => 25,
            49 => 36,
            64 => 45,
            81 => 60,
            n => n,
        };
    }
    Ok(d.element_rule_np(e, n).map_err(gauss)?.len())
}

/// `ParseElementalGaussRules`: the count after `-line`, `-tri`, `-quad`,
/// `-tetra`, `-pyramid`, `-prism`, `-brick`, by element family (2..8).
pub fn parse_elemental_gauss_rules(def: &str) -> Result<[usize; 9], SolverError> {
    let def = def.to_ascii_lowercase();
    let mut out = [0; 9];
    let names = ["-line", "-tri", "-quad", "-tetra", "-pyramid", "-prism", "-brick"];
    for (k, name) in names.iter().enumerate() {
        let Some(j) = def.find(name) else { continue };
        let rest = def[j + name.len()..].trim_start();
        let num: String = rest.chars().take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+').collect();
        out[k + 2] = num.parse::<i64>().map(|n| n.max(0) as usize).map_err(|_| {
            SolverError::fatal("ParseElementalGaussRules", format!("Problems reading integer from: {rest}"))
        })?;
    }
    Ok(out)
}

/// The fields as seen at one point of an element: every variable and the
/// coordinates interpolated with the basis values there.
/// `ListGetElementReal(Handle, Basis, Element, GaussPoint)` of `value`
/// (`key` of `section`) at one Gauss point: evaluated at the element's
/// nodes and interpolated with `basis` (`SUM(Basis(1:n) * F(1:n))`), unless
/// a dependency is an integration-point variable or `<key> At IP` is set,
/// when it is evaluated at the point itself (`EvaluateAtIp`).
#[allow(clippy::too_many_arguments)]
fn element_real_at_ip(
    value: &Value,
    section: &ParamMap,
    key: &str,
    sim: &Simulation,
    fields: &SimFields<'_>,
    nodes: &[u32],
    basis: &[f64],
) -> Result<f64, SolverError> {
    let b = &sim.backends;
    let at_ip = section.logical(&format!("{key} at ip")) == Some(true)
        || matches!(value, Value::Dependent(d) if d.variables.iter().flat_map(|v| v.split(','))
            .any(|v| sim.vars.get(&elmer_core::normalize_key(v)).is_some_and(|x| x.kind == FieldKind::Ip)));
    if at_ip {
        let at = AtPoint { inner: fields, nodes, basis };
        return Ok(crate::expr::eval_at_nodes_with(value, &[0], &at, b.matc.as_deref(), &b.udfs)?[0]);
    }
    let f = crate::expr::eval_at_nodes_with(value, nodes, fields, b.matc.as_deref(), &b.udfs)?;
    Ok(basis.iter().zip(&f).fold(0.0, |s, (b, f)| s + b * f))
}

/// Whether `value` depends on an integration-point variable.
pub(crate) fn depends_on_ip(sim: &Simulation, value: &Value) -> bool {
    matches!(value, Value::Dependent(d) if d.variables.iter().flat_map(|v| v.split(','))
        .any(|v| sim.vars.get(&elmer_core::normalize_key(v)).is_some_and(|x| x.kind == FieldKind::Ip)))
}

/// `ListGetElementReal(Handle, Basis, Element, GaussPoint)` of a boundary
/// element's `value` that depends on an IP variable (`EvaluateAtIp`), at
/// the boundary points `uvw`: the other fields are interpolated with the
/// boundary basis (`VarsToValuesOnIps`), the IP variables are
/// `InterpolateIPVariableToBoundary` (Lists.F90:5168): each parent with
/// points fits its values (`Ip2DgFieldInElement`), the fit is interpolated
/// with the boundary basis at the boundary's nodes, and two parents are
/// averaged.
pub(crate) fn boundary_ip_values(
    sim: &Simulation,
    mesh: &dyn MeshView,
    value: &Value,
    e: ElemId,
    uvw: &[[f64; 3]],
) -> Result<Vec<f64>, SolverError> {
    let fe = sim.backends.fe.as_ref().ok_or_else(|| SolverError::NotWired("basis and quadrature for IP variables".into()))?;
    let el = mesh.element(e);
    let n = el.nodes.len();
    let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
    let axi = sim
        .model
        .simulation
        .string("coordinate system")
        .is_some_and(|c| matches!(elmer_core::normalize_key(c).as_str(), "axi symmetric" | "cylindric symmetric"));
    let Value::Dependent(d) = value else { return Err(SolverError::fatal(CALLER, "not a dependent value")) };
    // Per IP variable: the parents' fits at the boundary's nodes.
    let mut fits: Vec<(String, Vec<Vec<f64>>)> = Vec::new();
    for name in d.variables.iter().flat_map(|v| v.split(',')).map(elmer_core::normalize_key) {
        let Some(v) = sim.vars.get(&name).filter(|v| v.kind == FieldKind::Ip) else { continue };
        if v.dofs > 1 {
            return Err(SolverError::Unsupported(format!("boundary evaluation of the {}-component IP variable {name}", v.dofs)));
        }
        let offs = v.perm.as_deref().ok_or_else(|| SolverError::fatal(CALLER, "IP field without a permutation"))?;
        let mut parents = Vec::new();
        for p in el.parents.iter().flatten() {
            let t = p.0 as usize;
            let (i1, i2) = (offs[t] as usize, offs[t + 1] as usize);
            if i2 == i1 {
                continue;
            }
            let pel = mesh.element(*p);
            let pxyz: Vec<[f64; 3]> = pel.nodes.iter().map(|&k| mesh.coord(k)).collect();
            let fdg = ip_to_dg_in_element(fe, pel.etype, &pxyz, &v.values[i1..i2], 1, axi)?.remove(0);
            // fdg at the boundary's nodes (None where the parent lacks one).
            parents.push(el.nodes.iter().map(|k| pel.nodes.iter().position(|q| q == k).map(|j| fdg[j])).collect::<Vec<_>>());
        }
        let parents: Vec<Vec<f64>> = parents.into_iter().map(|p| p.into_iter().map(|x| x.unwrap_or(f64::NAN)).collect()).collect();
        fits.push((name, parents));
    }
    let fields = sim.fields(mesh);
    let b = &sim.backends;
    let mut ev = elmer_core::BasisEval::default();
    let mut out = Vec::with_capacity(uvw.len());
    for &p in uvw {
        fe.basis.eval(el.etype, &xyz, p, &mut ev).map_err(|m| SolverError::fatal(CALLER, m))?;
        let basis = &ev.basis[..n];
        let ip: Vec<(String, f64)> = fits
            .iter()
            .map(|(name, parents)| {
                let mut t = 0.0;
                for f in parents {
                    for i in 0..n {
                        if !f[i].is_nan() {
                            t += basis[i] * f[i];
                        }
                    }
                }
                if parents.len() == 2 {
                    t /= 2.0;
                }
                (name.clone(), t)
            })
            .collect();
        let at = BoundaryIp { inner: AtPoint { inner: &fields, nodes: el.nodes, basis }, ip: &ip };
        out.push(crate::expr::eval_at_nodes_with(value, &[0], &at, b.matc.as_deref(), &b.udfs)?[0]);
    }
    Ok(out)
}

/// [`AtPoint`] with the IP variables' boundary values.
struct BoundaryIp<'a> {
    inner: AtPoint<'a>,
    ip: &'a [(String, f64)],
}

impl FieldLookup for BoundaryIp<'_> {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        let key = elmer_core::normalize_key(var);
        match self.ip.iter().find(|(n, _)| *n == key) {
            Some((_, v)) => Some(*v),
            None => self.inner.value(var, node),
        }
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.inner.coord(node)
    }
    fn time(&self) -> f64 {
        self.inner.time()
    }
}

struct AtPoint<'a> {
    inner: &'a SimFields<'a>,
    nodes: &'a [u32],
    basis: &'a [f64],
}

impl FieldLookup for AtPoint<'_> {
    fn value(&self, var: &str, _node: u32) -> Option<f64> {
        let mut s = 0.0;
        for (&n, &b) in self.nodes.iter().zip(self.basis) {
            s += b * self.inner.value(var, n)?;
        }
        Some(s)
    }

    fn coord(&self, _node: u32) -> [f64; 3] {
        let mut x = [0.0; 3];
        for (&n, &b) in self.nodes.iter().zip(self.basis) {
            let c = self.inner.coord(n);
            for k in 0..3 {
                x[k] += b * c[k];
            }
        }
        x
    }

    fn time(&self) -> f64 {
        self.inner.time()
    }
}

impl Simulation {
    /// `InitCond` (ElmerSolver.F90:2062) for an elemental or IP exported
    /// variable: from the initial condition of each element's body, the
    /// element mean of the nodal values (scalar elemental fields), or the
    /// value at each Gauss point (IP fields, component by component).
    pub(crate) fn exported_initial_values(
        &self,
        mesh: &dyn MeshView,
        spec: &crate::sim::VarSpec,
        perm: &[u32],
        values: &mut [f64],
    ) -> Result<(), SolverError> {
        let (kind, name, components, dofs) = (spec.kind, spec.name.as_str(), &spec.components, spec.dofs);
        let m = &self.model;
        let names: Vec<(usize, &str)> = if dofs == 1 {
            vec![(0, name)]
        } else {
            components.iter().enumerate().map(|(j, c)| (j, c.as_str())).collect()
        };
        if dofs > 1 && m.initial_conditions.iter().any(|ic| ic.get(name).is_some()) {
            let what = if kind == FieldKind::Ip { "elemental" } else { "elements" };
            return Err(SolverError::fatal("InitCond", format!("Initialization only for scalar {what} fields!")));
        }
        if kind == FieldKind::Elemental && dofs > 1 {
            return Ok(());
        }
        let fields = self.fields(mesh);
        let b = &self.backends;
        for (j, key) in names {
            if !m.initial_conditions.iter().any(|ic| ic.get(key).is_some()) {
                continue;
            }
            for t in 0..mesh.n_bulk() {
                let e = ElemId(t as u32);
                let Some(ic) = section_of(self, mesh, e, "initial condition")? else { continue };
                let Some(value) = ic.get(key) else { continue };
                let el = mesh.element(e);
                match kind {
                    FieldKind::Elemental => {
                        let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                        if let Some(&s) = perm.get(t).filter(|&&s| s != NO_DOF) {
                            values[s as usize] = local.iter().sum::<f64>() / local.len() as f64;
                        }
                    }
                    _ => {
                        let (i1, i2) = (perm[t] as usize, perm[t + 1] as usize);
                        if i2 == i1 {
                            continue;
                        }
                        let fe = b.fe.as_ref().ok_or_else(|| {
                            SolverError::NotWired("basis and quadrature for IP variables (elmer-fe)".into())
                        })?;
                        let qps = fe.quadrature.gauss_points(el.etype, None);
                        if qps.len() != i2 - i1 {
                            return Err(SolverError::Unsupported(
                                "IP variable whose Gauss point count changed (UpdateIpPerm)".into(),
                            ));
                        }
                        let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                        let mut ev = elmer_core::BasisEval::default();
                        for (k, qp) in qps.iter().enumerate() {
                            fe.basis
                                .eval(el.etype, &xyz, qp.uvw, &mut ev)
                                .map_err(|msg| SolverError::fatal("InitCond", msg))?;
                            values[(i1 + k) * dofs + j] =
                                element_real_at_ip(value, ic, key, self, &fields, el.nodes, &ev.basis[..el.nodes.len()])?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `UpdateExportedVariables` for slot `ix`.
    pub(crate) fn update_exported_variables(&mut self, ix: usize, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        log.info(CALLER, "Updating variables, if any!", 20);
        for l in 1.. {
            let Some(text) = self.slots[ix].params.string(&format!("exported variable {l}")).map(str::to_string) else {
                break;
            };
            let spec = crate::sim::parse_variable_spec(&text, None)?;
            let Some(id) = self.vars.id(&spec.name) else { continue };
            let (kind, dofs, len) = {
                let v = self.vars.by_id(id);
                (v.kind, v.dofs, v.values.len())
            };
            let names = self.vars.component_names(id);
            let comp = |j: usize| if dofs > 1 { names.get(j).cloned().unwrap_or_default() } else { spec.name.clone() };
            let m = &self.model;

            // A state variable: one value from the first body force giving it.
            if len == dofs || kind == FieldKind::Global {
                for j in 0..dofs {
                    let key = comp(j);
                    if let Some(bf) = m.body_forces.iter().find(|b| b.get(&key).is_some()) {
                        log.info(CALLER, "Found a proper definition for state variable", 6);
                        if let Some(v) = self.creal(bf, &key, mesh)? {
                            self.vars.by_id_mut(id).values[j] = v;
                        }
                    }
                }
                continue;
            }

            for j in 0..dofs {
                let key = comp(j);
                let cond = format!("{key} condition");
                let parts: Vec<bool> = m.body_forces.iter().map(|b| b.get(&key).is_some()).collect();
                let parts_bc: Vec<bool> = m.boundary_conditions.iter().map(|b| b.get(&key).is_some()).collect();
                let count = parts.iter().chain(&parts_bc).filter(|&&a| a).count();
                if count == 0 {
                    continue;
                }
                log.info(
                    CALLER,
                    &format!("Exported Variable {} defined in {} sections", i2s(l as i64), i2s(count as i64)),
                    8,
                );
                let mut writes: Vec<(usize, f64)> = Vec::new();
                let fields = self.fields(mesh);
                let b = &self.backends;
                let nb = mesh.n_bulk();
                for t in 0..nb + mesh.n_boundary() {
                    let e = ElemId(t as u32);
                    let el = mesh.element(e);
                    let mut section: Option<&ParamMap> = None;
                    if t < nb {
                        if let Some(bf) = m
                            .body_of_tag(el.tag)
                            .and_then(|bd| m.bodies[bd - 1].integer("body force"))
                            .filter(|&k| k >= 1 && parts.get(k as usize - 1) == Some(&true))
                        {
                            section = Some(&m.body_forces[bf as usize - 1]);
                        }
                    } else if let Some(bc) = m.bc_of_tag(el.tag).filter(|&i| parts_bc[i - 1]) {
                        section = Some(&m.boundary_conditions[bc - 1]);
                    }
                    let Some(section) = section else { continue };
                    let conditional = section.get(&cond).is_some();
                    let value = section.get(&key).expect("an active part");
                    let perm = self.vars.by_id(id).perm.clone();
                    let perm = perm.as_deref().unwrap_or(&[]);
                    match kind {
                        FieldKind::Ip => {
                            let (i1, i2) = (perm[t] as usize, perm[t + 1] as usize);
                            if i2 == i1 {
                                continue;
                            }
                            let fe = b.fe.as_ref().ok_or_else(|| {
                                SolverError::NotWired("basis and quadrature for IP variables (elmer-fe)".into())
                            })?;
                            let qps = fe.quadrature.gauss_points(el.etype, None);
                            if qps.len() != i2 - i1 {
                                return Err(SolverError::Unsupported(
                                    "IP variable whose Gauss point count changed (UpdateIpPerm)".into(),
                                ));
                            }
                            if conditional {
                                log.warn(CALLER, "Elemental variable cannot be conditional!");
                            }
                            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                            let mut ev = elmer_core::BasisEval::default();
                            for (k, qp) in qps.iter().enumerate() {
                                fe.basis
                                    .eval(el.etype, &xyz, qp.uvw, &mut ev)
                                    .map_err(|msg| SolverError::fatal(CALLER, msg))?;
                                let v = element_real_at_ip(value, section, &key, self, &fields, el.nodes, &ev.basis[..el.nodes.len()])?;
                                writes.push(((i1 + k) * dofs + j, v));
                            }
                        }
                        FieldKind::Elemental => {
                            if conditional {
                                log.warn(CALLER, "Elemental variables not conditional!");
                            }
                            let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                            let s = perm.get(t).copied().unwrap_or(NO_DOF);
                            if s != NO_DOF {
                                let mean = local.iter().sum::<f64>() / local.len() as f64;
                                writes.push((s as usize * dofs + j, mean));
                            }
                        }
                        FieldKind::Nodal => {
                            let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                            let slot = |n: u32| perm.get(n as usize).copied().filter(|&s| s != NO_DOF);
                            if conditional {
                                let c = crate::expr::eval_at_nodes_with(
                                    section.get(&cond).expect("checked"),
                                    el.nodes,
                                    &fields,
                                    b.matc.as_deref(),
                                    &b.udfs,
                                )?;
                                for (i, &n) in el.nodes.iter().enumerate() {
                                    if let (true, Some(s)) = (c[i] > 0.0, slot(n)) {
                                        writes.push((s as usize * dofs + j, local[i]));
                                    }
                                }
                            } else if el.nodes.iter().all(|&n| slot(n).is_some()) {
                                for (i, &n) in el.nodes.iter().enumerate() {
                                    writes.push((slot(n).expect("checked") as usize * dofs + j, local[i]));
                                }
                            }
                        }
                        FieldKind::Dg => {
                            // `VarIndexes => Element % DGIndexes`.
                            let local = crate::expr::eval_at_nodes_with(value, el.nodes, &fields, b.matc.as_deref(), &b.udfs)?;
                            let dg = self.dg_indexes.as_ref().ok_or_else(|| {
                                SolverError::fatal(CALLER, "a DG variable without the mesh's DG indexes")
                            })?;
                            let slot = |i: usize| {
                                let d = dg.index[t][i];
                                (d != NO_DOF).then(|| perm[d as usize]).filter(|&s| s != NO_DOF)
                            };
                            if conditional {
                                let c = crate::expr::eval_at_nodes_with(
                                    section.get(&cond).expect("checked"),
                                    el.nodes,
                                    &fields,
                                    b.matc.as_deref(),
                                    &b.udfs,
                                )?;
                                for i in 0..el.nodes.len() {
                                    if let (true, Some(s)) = (c[i] > 0.0, slot(i)) {
                                        writes.push((s as usize * dofs + j, local[i]));
                                    }
                                }
                            } else if (0..el.nodes.len()).all(|i| slot(i).is_some()) {
                                for (i, &x) in local.iter().enumerate() {
                                    writes.push((slot(i).expect("checked") as usize * dofs + j, x));
                                }
                            }
                        }
                        k => return Err(SolverError::Unsupported(format!("updating a {k:?} exported variable"))),
                    }
                }
                let v = self.vars.by_id_mut(id);
                for (i, x) in writes {
                    v.values[i] = x;
                }
                v.values_changed = true;
            }
        }
        self.project_variables(ix, mesh, log)?;
        self.publish_live();
        Ok(())
    }
}

/// `Ip2DgFieldInElement` (InterpolateMeshToMesh.F90:1803): the nodal
/// values of the local L2 projection of an element's `m` integration-point
/// values on `GaussPoints(Element, m)`, `MASS f = LOAD` solved with
/// `LuSolve`. `fip` holds the points' values point-major, `dofs` per point;
/// returns them per dof, one value per node. `axi` weights by the radius
/// (`CSymmetry`).
pub fn ip_to_dg_in_element(
    fe: &crate::backend::FeBackend,
    etype: elmer_core::ElementType,
    xyz: &[[f64; 3]],
    fip: &[f64],
    dofs: usize,
    axi: bool,
) -> Result<Vec<Vec<f64>>, SolverError> {
    let dofs = dofs.max(1);
    let m = fip.len() / dofs;
    let n = xyz.len();
    let points = fe.quadrature.gauss_points(etype, Some(m));
    if points.len() != m {
        return Err(SolverError::fatal("Ip2DgFieldInElement", "no rule with the field's number of points"));
    }
    // MASS and the basis at the points, the same for each component.
    let mut be = elmer_core::BasisEval::default();
    let mut mass = vec![0.0; n * n];
    let mut wb = Vec::with_capacity(m);
    for qp in &points {
        fe.basis.eval(etype, xyz, qp.uvw, &mut be).map_err(|e| SolverError::fatal("ElementInfo", e))?;
        let mut w = qp.weight * be.det_j;
        if axi {
            w *= (0..n).map(|k| be.basis[k] * xyz[k][0]).sum::<f64>();
        }
        for p in 0..n {
            for q in 0..n {
                mass[p * n + q] += w * be.basis[q] * be.basis[p];
            }
        }
        wb.push((w, be.basis[..n].to_vec()));
    }
    let mut out = Vec::with_capacity(dofs);
    for k in 0..dofs {
        let mut load = vec![0.0; n];
        for (i, (w, bas)) in wb.iter().enumerate() {
            let f = fip[dofs * i + k];
            for p in 0..n {
                load[p] += w * bas[p] * f;
            }
        }
        let mut a = mass.clone();
        crate::control::lu_solve(n, &mut a, &mut load)?;
        out.push(load);
    }
    Ok(out)
}

impl Simulation {
    /// The `Project Variable i` / `Project Type i` requests of slot `ix`
    /// (UpdateExportedVariables, SolverBasics.F90:6853): each IP variable
    /// is projected by `Ip2DgSwapper(Mesh, Var, ToType)` (SolverBasics.F90:
    /// 9222) into `<name>_nodal` (the default: element fits summed at the
    /// nodes, divided by the number of elements that hit them), `<name>_elem`
    /// (`elem`: the mean of the element's point values) or `<name>_dg`
    /// (`dg`: the fit at the element's DG indexes), over the bulk elements
    /// with points; an existing target is reused.
    fn project_variables(&mut self, ix: usize, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        const SWAPPER: &str = "Ip2DgSwapper";
        for l in 1.. {
            let Some(name) = self.slots[ix].params.string(&format!("project variable {l}")).map(str::to_string) else {
                break;
            };
            let name = crate::sim::parse_variable_spec(&name, None)?.name;
            let Some(id) = self.vars.id(&name) else {
                log.warn(CALLER, &format!("Could not find variable for projection: {name}"));
                continue;
            };
            if self.vars.by_id(id).kind != FieldKind::Ip {
                return Err(SolverError::fatal(CALLER, "Variable projection implemented only for IP variable!"));
            }
            let to = match self.slots[ix].params.string(&format!("project type {l}")).map(elmer_core::normalize_key).as_deref() {
                Some("elem") => FieldKind::Elemental,
                Some("dg") => FieldKind::Dg,
                _ => FieldKind::Nodal,
            };
            // `Element % DGIndexes`: Elmer has them from the setup of any
            // Discontinuous Galerkin solver; here they are made on first use.
            if to == FieldKind::Dg
                && self.dg_indexes.is_none()
                && self.slots.iter().any(|s| s.params.logical("discontinuous galerkin") == Some(true))
            {
                self.ensure_dg_indexes(mesh);
            }
            let fe = self.backends.fe.as_ref().ok_or_else(|| SolverError::NotWired("basis and quadrature for Ip2DgSwapper".into()))?;
            let from = self.vars.by_id(id);
            let suffix = match to {
                FieldKind::Dg => "_dg",
                FieldKind::Elemental => "_elem",
                _ => "_nodal",
            };
            let to_name = format!("{}{suffix}", from.name);
            let offs = from.perm.as_deref().ok_or_else(|| SolverError::fatal(SWAPPER, "IP field without a permutation"))?;
            let dofs = from.dofs.max(1);
            let nb = mesh.n_bulk();
            let has_points = |t: usize| offs[t + 1] > offs[t];
            let dg = if to == FieldKind::Dg { self.dg_indexes.as_ref() } else { None };
            let permsize = match to {
                FieldKind::Dg => (0..nb).map(|t| mesh.element(ElemId(t as u32)).nodes.len()).sum(),
                FieldKind::Elemental => nb + mesh.n_boundary(),
                _ => mesh.n_nodes(),
            };
            // The target's indexes of element t.
            let targets = |t: usize| -> Option<Vec<usize>> {
                match to {
                    FieldKind::Dg => dg.map(|d| d.index[t].iter().map(|&i| i as usize).collect()),
                    FieldKind::Elemental => Some(vec![t]),
                    _ => Some(mesh.element(ElemId(t as u32)).nodes.iter().map(|&k| k as usize).collect()),
                }
            };
            let mut mark = vec![false; permsize];
            for t in (0..nb).filter(|&t| has_points(t)) {
                let Some(idx) = targets(t) else {
                    log.warn(SWAPPER, "Cannot project to DG field without DGIndexes!");
                    break;
                };
                for i in idx {
                    mark[i] = true;
                }
            }
            let mut perm = vec![NO_DOF; permsize];
            let mut size = 0usize;
            for (p, _) in perm.iter_mut().zip(&mark).filter(|(_, &m)| m) {
                *p = size as u32;
                size += 1;
            }
            let mut values = vec![0.0f64; size * dofs];
            let mut hits = vec![0usize; size];
            let axi = self
                .model
                .simulation
                .string("coordinate system")
                .is_some_and(|c| matches!(elmer_core::normalize_key(c).as_str(), "axi symmetric" | "cylindric symmetric"));
            for t in (0..nb).filter(|&t| has_points(t)) {
                let Some(idx) = targets(t) else { continue };
                if idx.iter().all(|&i| perm[i] == NO_DOF) {
                    continue;
                }
                let (i1, i2) = (offs[t] as usize, offs[t + 1] as usize);
                let fip = &from.values[dofs * i1..dofs * i2];
                if to == FieldKind::Elemental {
                    let j = perm[t] as usize;
                    let m = (i2 - i1) as f64;
                    for k in 0..dofs {
                        values[dofs * j + k] = (0..i2 - i1).map(|i| fip[dofs * i + k]).sum::<f64>() / m;
                    }
                    continue;
                }
                let el = mesh.element(ElemId(t as u32));
                let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| mesh.coord(k)).collect();
                let fdg = ip_to_dg_in_element(fe, el.etype, &xyz, fip, dofs, axi)?;
                for (k, f) in fdg.iter().enumerate() {
                    for (i, &d) in idx.iter().enumerate() {
                        let j = perm[d];
                        if j == NO_DOF {
                            continue;
                        }
                        let j = j as usize;
                        if to == FieldKind::Dg {
                            values[dofs * j + k] = f[i];
                        } else {
                            values[dofs * j + k] += f[i];
                            if k == 0 {
                                hits[j] += 1;
                            }
                        }
                    }
                }
            }
            if to == FieldKind::Nodal {
                for (j, &h) in hits.iter().enumerate() {
                    if h > 0 {
                        for k in 0..dofs {
                            values[dofs * j + k] /= h as f64;
                        }
                    }
                }
            }
            let perm: std::sync::Arc<[u32]> = perm.into();
            match self.vars.get_mut(&to_name) {
                Some(v) => {
                    v.perm = Some(perm);
                    v.values = values;
                    v.values_changed = true;
                }
                None => {
                    let mut v = crate::variable::Variable::new(&to_name, to, dofs, Some(perm), values);
                    v.values_changed = true;
                    self.vars.add(v)?;
                }
            }
        }
        Ok(())
    }

    /// `CheckAndCreateDGIndexes(Mesh)`: the mesh's DG indexes, created when
    /// no DG variable has made them yet.
    pub fn ensure_dg_indexes(&mut self, mesh: &dyn MeshView) -> &DgIndexes {
        self.dg_indexes.get_or_insert_with(|| DgIndexes::create(mesh))
    }

    /// `AverageBodyFields(Mesh)` of VtuOutputSolver: every output DG field
    /// with one dof and more than one value gets `CalculateBodyAverage(Mesh,
    /// Var, .FALSE.)`. Returns the number of fields reduced.
    pub fn average_body_fields(&mut self, mesh: &dyn MeshView) -> Result<usize, SolverError> {
        let ids: Vec<crate::VarId> = self
            .vars
            .iter()
            .filter(|(_, v)| v.output && v.values.len() != v.dofs && v.dofs == 1 && v.kind == FieldKind::Dg)
            .map(|(id, _)| id)
            .collect();
        for &id in &ids {
            self.calculate_body_average(mesh, id)?;
        }
        Ok(ids.len())
    }

    /// `CalculateBodyAverage(Mesh, Var, BodySum = .FALSE.)` (MeshBasics.F90):
    /// per body and component, each node's mean over the DG values of the
    /// body's elements at that node replaces those values. Once per variable
    /// (`DgAveraged`).
    pub fn calculate_body_average(&mut self, mesh: &dyn MeshView, id: crate::VarId) -> Result<(), SolverError> {
        let n_bodies = self.model.bodies.len();
        let nb = mesh.n_bulk();
        let bodies: Vec<Option<usize>> =
            (0..nb).map(|t| self.model.body_of_tag(mesh.element(ElemId(t as u32)).tag)).collect();
        let dg = self
            .dg_indexes
            .as_ref()
            .ok_or_else(|| SolverError::fatal("CalculateBodyAverage", "no DG indexes"))?;
        let v = self.vars.by_id_mut(id);
        let Some(perm) = v.perm.clone() else { return Ok(()) };
        if perm.len() <= mesh.n_nodes() || v.dg_averaged {
            return Ok(());
        }
        let dofs = v.dofs;
        let n = mesh.n_nodes();
        let mut sum = vec![0.0f64; n];
        let mut count = vec![0.0f64; n];
        let slot = |t: usize, l: usize| {
            let d = dg.index[t][l];
            (d != NO_DOF).then(|| perm[d as usize]).filter(|&s| s != NO_DOF).map(|s| s as usize)
        };
        for body in 1..=n_bodies {
            for k in 0..dofs {
                sum.iter_mut().for_each(|x| *x = 0.0);
                count.iter_mut().for_each(|x| *x = 0.0);
                for t in (0..nb).filter(|&t| bodies[t] == Some(body)) {
                    for (l, &node) in mesh.element(ElemId(t as u32)).nodes.iter().enumerate() {
                        if let Some(s) = slot(t, l) {
                            sum[node as usize] += v.values[dofs * s + k];
                            count[node as usize] += 1.0;
                        }
                    }
                }
                // Parallel: the interface nodes' sums and counts of the
                // other partitions, halo elements included on every holder.
                if let Some(nodes) = &self.backends.nodes {
                    nodes.sum_shared(&mut sum)?;
                    nodes.sum_shared(&mut count)?;
                }
                if count.iter().all(|&c| c == 0.0) {
                    continue;
                }
                for (x, &c) in sum.iter_mut().zip(&count) {
                    if c > 0.0 {
                        *x /= c;
                    }
                }
                for t in (0..nb).filter(|&t| bodies[t] == Some(body)) {
                    for (l, &node) in mesh.element(ElemId(t as u32)).nodes.iter().enumerate() {
                        if let Some(s) = slot(t, l) {
                            v.values[dofs * s + k] = sum[node as usize];
                        }
                    }
                }
            }
        }
        v.dg_averaged = true;
        Ok(())
    }
}
