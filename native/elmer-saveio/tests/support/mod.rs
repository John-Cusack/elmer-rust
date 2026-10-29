//! Test-only helpers: an Elmer mesh-file reader behind `MeshView` and a
//! `FieldSource` rebuilt from a parsed `.result`.
#![allow(dead_code)]
use elmer_core::field::{FieldKind, FieldRef, NO_DOF};
use elmer_core::{ElemId, ElementRef, ElementType, FieldSource, MeshView};
use elmer_saveio::result::{PermRecord, ResultFile};
use std::path::PathBuf;

pub fn fx_path(p: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(p)
}

pub fn fx(p: &str) -> Vec<u8> {
    std::fs::read(fx_path(p)).unwrap_or_else(|e| panic!("{p}: {e}"))
}

pub fn fx_str(p: &str) -> String {
    String::from_utf8(fx(p)).unwrap()
}

struct Elem {
    etype: u16,
    nodes: Vec<u32>,
    tag: i32,
}

pub struct TestMesh {
    dim: usize,
    coords: Vec<[f64; 3]>,
    elems: Vec<Elem>,
    n_bulk: usize,
}

impl TestMesh {
    pub fn load(dir: &str) -> TestMesh {
        let nodes = fx_str(&format!("{dir}/mesh.nodes"));
        let coords: Vec<[f64; 3]> = nodes
            .lines()
            .map(|l| {
                let t: Vec<f64> = l
                    .split_whitespace()
                    .skip(2)
                    .map(|x| x.parse().unwrap())
                    .collect();
                [t[0], t[1], t[2]]
            })
            .collect();
        let mut elems = Vec::new();
        for l in fx_str(&format!("{dir}/mesh.elements")).lines() {
            let t: Vec<i64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            elems.push(Elem {
                etype: t[2] as u16,
                tag: t[1] as i32,
                nodes: t[3..].iter().map(|&n| n as u32 - 1).collect(),
            });
        }
        let n_bulk = elems.len();
        for l in fx_str(&format!("{dir}/mesh.boundary")).lines() {
            let t: Vec<i64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            elems.push(Elem {
                etype: t[4] as u16,
                tag: t[1] as i32,
                nodes: t[5..].iter().map(|&n| n as u32 - 1).collect(),
            });
        }
        let dim = if coords.iter().any(|c| c[2] != 0.0) {
            3
        } else {
            2
        };
        TestMesh {
            dim,
            coords,
            elems,
            n_bulk,
        }
    }
}

impl MeshView for TestMesh {
    fn dim(&self) -> usize {
        self.dim
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let el = &self.elems[e.0 as usize];
        ElementRef {
            etype: ElementType(el.etype),
            nodes: &el.nodes,
            tag: el.tag,
            parents: [None, None],
        }
    }
}

/// One step of a parsed `.result` turned back into live fields: scalar
/// components that share a table and a size are re-interleaved into the
/// vector field named in the header.
pub struct OwnedField {
    pub name: String,
    pub dofs: usize,
    pub perm: Option<Vec<u32>>,
    pub values: Vec<f64>,
}

pub struct FixtureSource {
    pub fields: Vec<OwnedField>,
    pub time: f64,
    pub step: usize,
}

impl FixtureSource {
    pub fn from_result(rf: &ResultFile, step: usize) -> FixtureSource {
        let st = &rf.steps[step];
        let n_nodes = rf.header.n_nodes;
        // Dense scalar fields by file order.
        let mut scalars: Vec<(String, Option<Vec<u32>>, Vec<f64>)> = Vec::new();
        let mut table: Option<Vec<u32>> = None;
        for v in &st.vars {
            let (perm, values) = match &v.perm {
                PermRecord::Null => (None, v.values.clone()),
                p => {
                    if let PermRecord::Table { size, entries } = p {
                        let mut t = vec![NO_DOF; *size];
                        for (i, k) in entries {
                            t[*i as usize - 1] = k - 1;
                        }
                        table = Some(t);
                    }
                    let t = table.clone().unwrap();
                    let n = t.iter().filter(|&&s| s != NO_DOF).count();
                    let mut dense = vec![0.0; n];
                    let mut it = v.values.iter();
                    for &s in &t {
                        if s != NO_DOF {
                            dense[s as usize] = *it.next().unwrap();
                        }
                    }
                    (Some(t), dense)
                }
            };
            scalars.push((v.name.clone(), perm, values));
        }
        let mut fields = Vec::new();
        let mut used = vec![false; scalars.len()];
        // Vector headers: gather their components from the scalar list.
        for h in rf.header.vars.iter().filter(|h| h.dofs > 1) {
            let comps: Vec<usize> = (1..=h.dofs)
                .map(|c| {
                    let n = elmer_saveio::writer::component_name(&h.name, c, h.dofs).unwrap();
                    scalars.iter().position(|s| s.0 == n).expect("component")
                })
                .collect();
            let n = scalars[comps[0]].2.len();
            let mut values = vec![0.0; n * h.dofs];
            for (c, &i) in comps.iter().enumerate() {
                used[i] = true;
                for s in 0..n {
                    values[s * h.dofs + c] = scalars[i].2[s];
                }
            }
            fields.push(OwnedField {
                name: h.name.clone(),
                dofs: h.dofs,
                perm: scalars[comps[0]].1.clone(),
                values,
            });
        }
        for (i, s) in scalars.into_iter().enumerate() {
            if !used[i] {
                fields.push(OwnedField {
                    name: s.0,
                    dofs: 1,
                    perm: s.1,
                    values: s.2,
                });
            }
        }
        let _ = n_nodes;
        FixtureSource {
            fields,
            time: st.time,
            step: st.timestep as usize,
        }
    }
}

impl FieldSource for FixtureSource {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        self.fields
            .iter()
            .map(|f| FieldRef {
                name: &f.name,
                kind: FieldKind::Nodal,
                dofs: f.dofs,
                perm: f.perm.as_deref(),
                values: &f.values,
            })
            .collect()
    }
    fn time(&self) -> f64 {
        self.time
    }
    fn step(&self) -> usize {
        self.step
    }
}
