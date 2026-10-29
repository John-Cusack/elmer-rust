//! `MeshExtrude` (MeshExtrusion.F90:425) for the primary mesh, driven by the
//! Simulation `Extruded Mesh *` keywords (`CreateExtrudedMesh`,
//! ElmerSolver.F90:903).
//!
//! Layer `l` of `0..=layers` repeats the old nodes with the extruded
//! coordinate set to `w[l]`, so new node `l*n + i` sits above old node `i`.
//! Elements come in this order: bulk elements layer by layer; side
//! boundaries layer by layer (old boundaries tagged 0 are dropped); one copy
//! of the old bulk mesh per BC layer (bottom, the `Extruded BC Layers`,
//! top); and with `Preserve Baseline` the old boundaries again on the
//! baseline layer, parented by that layer's plane elements.

use elmer_core::{ElemId, ElementType, MeshView, ParamMap, Params, Value};

use crate::split::check_prepare_mesh;
use crate::{division, Element, Mesh, MeshLevelsError};

/// The resolved extrusion settings of one keyword list.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtrudeOptions {
    /// Coordinates of the `layers + 1` node layers (`ExtrudedDivision`).
    pub w: Vec<f64>,
    /// 0-based coordinate overwritten by `w` (`Extruded Coordinate Index`).
    pub coord: usize,
    /// The layer indices that get BC planes; first 0, last `layers`.
    pub bc_layers: Vec<usize>,
    pub preserve_baseline: bool,
    /// 1-based index into `bc_layers` (`Extruded Baseline Layer`).
    pub baseline_layer: usize,
    /// `Extruded BCs Collect`: one BC per plane instead of one per body.
    pub collect_bcs: bool,
}

impl ExtrudeOptions {
    pub fn layers(&self) -> usize {
        self.w.len() - 1
    }

    /// The options in `simulation`, or `None` without an `Extruded Mesh`
    /// keyword. `density` evaluates `Extruded Mesh Density` (or
    /// `1D Mesh Density`) at a point of the unit segment when that keyword
    /// is not a constant; `dim` is the dimension of the mesh to extrude.
    pub fn from_params(
        simulation: &ParamMap,
        dim: usize,
        density: Option<&mut dyn FnMut(f64) -> Result<f64, String>>,
    ) -> Result<Option<Self>, MeshLevelsError> {
        let sim = simulation;
        if !sim.0.keys().any(|k| k.starts_with("extruded mesh")) {
            return Ok(None);
        }
        check_prepare_mesh(sim)?;
        let refuse = |what: &str| Err(MeshLevelsError::Unsupported(what.to_string()));
        if sim.logical("extruded mesh slices") == Some(true) {
            return refuse("Extruded Mesh Slices");
        }
        if sim.logical("extruded mesh rotational") == Some(true) {
            return refuse("Extruded Mesh Rotational");
        }
        if sim.get("extruded mesh end map").is_some() {
            return refuse("Extruded Mesh End Map");
        }
        let skew = ["rotor skew", "rotor skew function", "stator skew"]
            .iter()
            .any(|k| sim.get(k).is_some());
        if skew && (sim.get("rotor bodies").is_some() || sim.get("rotor radius").is_some()) {
            return refuse("rotor skew of an extruded mesh");
        }

        let layers = match int(sim, "extruded mesh layers")? {
            Some(n) => n,
            None => match int(sim, "extruded mesh levels")? {
                Some(n) => n - 1,
                None if sim.get("extruded elements").is_some() => {
                    return refuse("Extruded Elements")
                }
                None => {
                    return Err(MeshLevelsError::Invalid(
                        "no Extruded Mesh Layers or Levels".into(),
                    ))
                }
            },
        };
        if layers < 2 {
            return Err(MeshLevelsError::Invalid(
                "there must be at least two Extruded Mesh Layers".into(),
            ));
        }
        let n = layers as usize;
        let mut w = match ["extruded mesh density", "1d mesh density"]
            .iter()
            .find_map(|k| sim.get(k))
        {
            Some(Value::Real { data, .. }) if data.len() == 1 => {
                let h = data[0];
                division::function(n, &mut |_| Ok(h))?
            }
            Some(_) => match density {
                Some(f) => division::function(n, f)?,
                None => {
                    return refuse(
                        "Extruded Mesh Density that is not a constant, without an evaluator",
                    )
                }
            },
            None => {
                let q = real(sim, "extruded mesh ratio")?
                    .or(real(sim, "1d mesh ratio")?)
                    .unwrap_or(1.0);
                division::geometric(n, q)?
            }
        };
        let zmin = real(sim, "extruded min coordinate")?.unwrap_or(0.0);
        let zmax = real(sim, "extruded max coordinate")?.unwrap_or(1.0);
        for x in &mut w {
            *x = zmin + (zmax - zmin) * *x;
        }

        let coord = match int(sim, "extruded coordinate index")? {
            Some(c @ 1..=3) => c as usize - 1,
            Some(c) => {
                return Err(MeshLevelsError::Invalid(format!(
                    "Extruded Coordinate Index {c}"
                )))
            }
            None => dim.min(2),
        };
        let mut bc_layers = vec![0];
        if let Some(v) = sim.get("extruded bc layers") {
            let Value::Integer(v) = v else {
                return Err(MeshLevelsError::Invalid(
                    "Extruded BC Layers is not an integer array".into(),
                ));
            };
            for &l in v {
                bc_layers.push(
                    usize::try_from(l)
                        .map_err(|_| MeshLevelsError::Invalid(format!("BC layer {l}")))?,
                );
            }
        }
        bc_layers.push(n);
        if bc_layers.windows(2).any(|p| p[0] >= p[1]) {
            return Err(MeshLevelsError::Invalid(
                "BC layers should be in increasing order".into(),
            ));
        }
        let preserve_baseline = sim.logical("preserve baseline") == Some(true);
        let mut baseline_layer = 0;
        if preserve_baseline {
            baseline_layer = int(sim, "extruded baseline layer")?.unwrap_or(1).max(0) as usize;
            if baseline_layer > bc_layers.len() || baseline_layer == 0 {
                return Err(MeshLevelsError::Invalid(format!(
                    "Extruded Baseline Layer must be within 1..{}",
                    bc_layers.len()
                )));
            }
        }
        Ok(Some(ExtrudeOptions {
            w,
            coord,
            bc_layers,
            preserve_baseline,
            baseline_layer,
            collect_bcs: sim.logical("extruded bcs collect") == Some(true),
        }))
    }
}

/// `ListGetInteger`, also taking an untyped integral real.
fn int(p: &ParamMap, k: &str) -> Result<Option<i32>, MeshLevelsError> {
    match p.get(k) {
        None => Ok(None),
        Some(Value::Integer(v)) if !v.is_empty() => Ok(Some(v[0])),
        Some(Value::Real { data, .. }) if data.len() == 1 && data[0].fract() == 0.0 => {
            Ok(Some(data[0] as i32))
        }
        Some(_) => Err(MeshLevelsError::Invalid(format!("{k} is not an integer"))),
    }
}

/// `ListGetCReal`: a constant real (an integer is promoted).
fn real(p: &ParamMap, k: &str) -> Result<Option<f64>, MeshLevelsError> {
    match p.get(k) {
        None => Ok(None),
        Some(Value::Real { data, .. }) if data.len() == 1 => Ok(Some(data[0])),
        Some(Value::Integer(v)) if v.len() == 1 => Ok(Some(v[0] as f64)),
        Some(_) => Err(MeshLevelsError::Unsupported(format!(
            "{k} that is not a constant real"
        ))),
    }
}

/// Extrudes a 1D or 2D mesh whose tags are Elmer's mapped body and BC
/// indices (see the crate docs). The new boundary tags are BC indices too:
/// side boundaries keep (and offset) the old ones, planes continue after the
/// largest of them. BC indices past the sif's BC count need the dummy BCs
/// `AppendMissingBCs` adds; that is the caller's model, not the mesh.
pub fn extrude(mesh: &dyn MeshView, opts: &ExtrudeOptions) -> Result<Mesh, MeshLevelsError> {
    let mut old = Mesh::from_view(mesh);
    let n = old.coords.len();
    let nb = old.n_bulk;
    let in_levels = opts.layers() - 1;
    check_point_element_parents(&mut old);

    let mut coords = Vec::with_capacity((in_levels + 2) * n);
    for &z in &opts.w {
        coords.extend(old.coords.iter().map(|&c| {
            let mut c = c;
            c[opts.coord] = z;
            c
        }));
    }

    let n_bc_layers = opts.bc_layers.len();
    // 1-based BC-layer interval holding element layer `i` (the Fortran loop
    // index left at exit).
    let interval = |i: usize| {
        (1..n_bc_layers)
            .find(|&k| i < opts.bc_layers[k])
            .unwrap_or(n_bc_layers) as i32
    };
    let max_body = old.bulk().iter().map(|e| e.tag).max().unwrap_or(0).max(0);
    let max_bid0 = old
        .boundary()
        .iter()
        .map(|e| e.tag)
        .max()
        .unwrap_or(0)
        .max(0);

    let mut elements = Vec::new();
    for i in 0..=in_levels {
        for el in old.bulk() {
            let mut tag = el.tag;
            if n_bc_layers > 2 {
                tag += max_body * (interval(i) - 1);
            }
            let lo = el.nodes.iter().map(|&v| v + (i * n) as u32);
            let hi = el.nodes.iter().map(|&v| v + ((i + 1) * n) as u32);
            let mut nodes: Vec<u32> = lo.chain(hi).collect();
            let etype = match nodes.len() {
                4 => {
                    nodes.swap(2, 3);
                    404
                }
                6 => 706,
                8 => 808,
                _ => {
                    return Err(MeshLevelsError::Unsupported(format!(
                        "extruding bulk element type {}",
                        el.etype.0
                    )))
                }
            };
            elements.push(Element {
                etype: ElementType(etype),
                nodes,
                tag,
                parents: [None, None],
            });
        }
    }
    let n_bulk = elements.len();

    let bcoffset = if opts.preserve_baseline { max_bid0 } else { 0 };
    let mut max_bid = 0;
    for i in 0..=in_levels {
        for el in old.boundary() {
            if el.tag == 0 {
                continue;
            }
            let mut tag = el.tag + bcoffset;
            if n_bc_layers > 2 {
                tag += max_bid0 * (interval(i) - 1);
            }
            max_bid = max_bid.max(tag);
            let (lo, hi) = ((i * n) as u32, ((i + 1) * n) as u32);
            let v = &el.nodes;
            let (etype, nodes) = match el.etype.family() {
                2 => (404, vec![v[0] + lo, v[1] + lo, v[1] + hi, v[0] + hi]),
                1 => (202, vec![v[0] + lo, v[0] + hi]),
                f => {
                    return Err(MeshLevelsError::Unsupported(format!(
                        "extruding boundary family {f}"
                    )))
                }
            };
            let parents = el.parents.map(|p| p.map(|p| ElemId((nb * i) as u32 + p.0)));
            elements.push(Element {
                etype: ElementType(etype),
                nodes,
                tag,
                parents,
            });
        }
    }

    let bcoffset = max_bid;
    let mut baseline0 = 0;
    for (k, &level) in opts.bc_layers.iter().enumerate() {
        let k = k + 1;
        if opts.preserve_baseline && k == opts.baseline_layer {
            baseline0 = elements.len();
        }
        for (i, el) in old.bulk().iter().enumerate() {
            let parents = if level == in_levels + 1 {
                [Some((level - 1) * nb + i), None]
            } else {
                [
                    Some(level * nb + i),
                    (level > 0).then(|| (level - 1) * nb + i),
                ]
            };
            let tag = if opts.collect_bcs {
                bcoffset + k as i32
            } else {
                bcoffset + (k as i32 - 1) * max_body + el.tag
            };
            let nodes = el.nodes.iter().map(|&v| v + (level * n) as u32).collect();
            elements.push(Element {
                etype: el.etype,
                nodes,
                tag,
                parents: parents.map(|p| p.map(|p| ElemId(p as u32))),
            });
        }
    }

    if opts.preserve_baseline {
        let shift = (opts.bc_layers[opts.baseline_layer - 1] * n) as u32;
        for el in old.boundary() {
            if el.tag == 0 {
                continue;
            }
            let nodes = el.nodes.iter().map(|&v| v + shift).collect();
            let parents = el
                .parents
                .map(|p| p.map(|p| ElemId(p.0 + baseline0 as u32)));
            elements.push(Element {
                etype: el.etype,
                nodes,
                tag: el.tag,
                parents,
            });
        }
    }

    Ok(Mesh {
        dim: (old.dim + 1).min(3),
        coords,
        elements,
        n_bulk,
    })
}

/// `CheckPointElementParents`: a 101 boundary element whose left parent does
/// not contain its node gets the first bulk element that does.
fn check_point_element_parents(mesh: &mut Mesh) {
    let nb = mesh.n_bulk;
    let mut owner = vec![None; mesh.coords.len()];
    for (t, el) in mesh.elements[..nb].iter().enumerate() {
        for &v in &el.nodes {
            owner[v as usize].get_or_insert(t as u32);
        }
    }
    let (bulk, boundary) = mesh.elements.split_at_mut(nb);
    for el in boundary.iter_mut().filter(|e| e.nodes.len() == 1) {
        let v = el.nodes[0];
        if let Some(p) = el.parents[0] {
            let parent_nodes = &bulk.get(p.0 as usize).map(|e| &e.nodes[..]).unwrap_or(&[]);
            if !parent_nodes.contains(&v) {
                el.parents[0] = owner[v as usize].map(ElemId);
            }
        }
    }
}
