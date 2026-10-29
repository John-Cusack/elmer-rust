//! Test-local stand-in for the mesh front end: reads an Elmer mesh DB with
//! boundary parents, maps tags the way `MapBodiesAndBCs` does, and writes a
//! mesh back in `WriteMeshToDisk2` form (MeshIO.F90:1281), which is what the
//! pristine `SaveMesh` solver produced for the expected files.

use std::path::Path;

use elmer_core::{ElemId, ElementType, ParamMap, Params, Value};
use elmer_mesh_levels::{Element, Mesh};

fn rows(path: &Path) -> Vec<Vec<i64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .map(|l| {
            l.split_whitespace()
                .map(|s| s.parse().unwrap())
                .collect::<Vec<_>>()
        })
        .filter(|l| !l.is_empty())
        .collect()
}

/// The mesh DB in `dir` with its file tags; the dimension comes from the
/// bulk elements.
pub fn read_mesh(dir: &Path) -> Mesh {
    let text = std::fs::read_to_string(dir.join("mesh.nodes")).unwrap();
    let coords: Vec<[f64; 3]> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .enumerate()
        .map(|(i, l)| {
            let f: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(
                f[0].parse::<usize>().unwrap(),
                i + 1,
                "non-contiguous node ids"
            );
            [
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
            ]
        })
        .collect();
    let node = |v: i64| (v - 1) as u32;
    let mut elements: Vec<Element> = rows(&dir.join("mesh.elements"))
        .iter()
        .map(|r| Element {
            etype: ElementType(r[2] as u16),
            nodes: r[3..].iter().map(|&v| node(v)).collect(),
            tag: r[1] as i32,
            parents: [None, None],
        })
        .collect();
    let n_bulk = elements.len();
    let parent = |p: i64| (p > 0).then(|| ElemId(p as u32 - 1));
    elements.extend(rows(&dir.join("mesh.boundary")).iter().map(|r| Element {
        etype: ElementType(r[4] as u16),
        nodes: r[5..].iter().map(|&v| node(v)).collect(),
        tag: r[1] as i32,
        parents: [parent(r[2]), parent(r[3])],
    }));
    let dim = match elements[..n_bulk].iter().map(|e| e.etype.family()).max() {
        Some(2) => 1,
        Some(3 | 4) => 2,
        _ => 3,
    };
    Mesh {
        dim,
        coords,
        elements,
        n_bulk,
    }
}

pub struct Case {
    pub simulation: ParamMap,
    pub bodies: Vec<ParamMap>,
    pub bcs: Vec<ParamMap>,
    pub mesh: Mesh,
}

/// Parses `dir/sif` and reads the mesh DB its Header names.
pub fn load(dir: &Path, sif: &str) -> Case {
    let parsed = elmer_sif::parse_file(dir.join(sif).to_str().unwrap()).unwrap();
    let model = parsed.sif.merged();
    let simulation = model
        .simulation
        .as_ref()
        .map(elmer_sif::section_to_param_map)
        .unwrap_or_default();
    // elmer-sif pads unnumbered BC slots with empty sections, so the
    // position is the BC index.
    let mut bcs: Vec<ParamMap> = model
        .bcs
        .iter()
        .map(elmer_sif::section_to_param_map)
        .collect();
    let mut bodies: Vec<ParamMap> = model
        .bodies
        .iter()
        .map(elmer_sif::section_to_param_map)
        .collect();
    let (db, name) = model.header.mesh_db.clone().expect("Mesh DB");
    let mesh_dir = dir.join(db).join(name);
    let mesh = read_mesh(&mesh_dir);
    if simulation.logical("use mesh names") == Some(true) {
        read_target_names(&mesh_dir.join("mesh.names"), &mut bodies, &mut bcs);
    }
    Case {
        simulation,
        bodies,
        bcs,
        mesh,
    }
}

/// `ReadTargetNames` (MeshIO.F90:959): `$ name = tags` lines of `mesh.names`
/// give `Target Bodies` / `Target Boundaries` to the section of that `Name`.
fn read_target_names(file: &Path, bodies: &mut [ParamMap], bcs: &mut [ParamMap]) {
    let text = std::fs::read_to_string(file).unwrap();
    for line in text.lines() {
        let (Some(i), Some(j)) = (line.find('$'), line.find('=')) else {
            continue;
        };
        let name = line[i + 1..j].trim().to_lowercase();
        let tags: Vec<i32> = line[j + 1..]
            .split_whitespace()
            .map_while(|w| w.parse().ok())
            .collect();
        let mut already_set = false;
        for (sections, key) in [
            (&mut *bcs, "Target Boundaries"),
            (&mut *bodies, "Target Bodies"),
        ] {
            for section in sections.iter_mut() {
                if section.string("name").map(str::to_lowercase).as_deref() == Some(name.as_str()) {
                    assert!(!already_set, "Mapping of name is not unique: {name}");
                    if section.get(key).is_none() {
                        section.insert(key, Value::Integer(tags.clone()));
                        already_set = true;
                    }
                }
            }
        }
    }
}

fn targets<'a>(section: &'a ParamMap, key: &str) -> Option<&'a [i32]> {
    match section.get(key) {
        Some(Value::Integer(v)) => Some(v),
        _ => None,
    }
}

/// `MapBodiesAndBCs` for boundaries: file tag -> BC index (the first BC that
/// targets it), else the last BC with `Default Target` or `Default BC`, else
/// 0.
pub fn map_bc_tags(mesh: &mut Mesh, bcs: &[ParamMap]) {
    let nb = mesh.n_bulk;
    let default = bcs
        .iter()
        .rposition(|bc| {
            bc.logical("default target") == Some(true) || bc.logical("default bc") == Some(true)
        })
        .map_or(0, |i| i as i32 + 1);
    for el in &mut mesh.elements[nb..] {
        el.tag = bcs
            .iter()
            .position(|bc| targets(bc, "target boundaries").is_some_and(|t| t.contains(&el.tag)))
            .map_or(default, |i| i as i32 + 1);
    }
}

/// `MapBodiesAndBCs` for bodies: only when some body has `Target Bodies`;
/// a body without them claims its own index.
pub fn map_body_tags(mesh: &mut Mesh, bodies: &[ParamMap]) {
    if !bodies.iter().any(|b| targets(b, "target bodies").is_some()) {
        return;
    }
    let max = mesh
        .bulk()
        .iter()
        .map(|e| e.tag)
        .max()
        .unwrap_or(0)
        .max(bodies.len() as i32);
    let mut map = vec![0; max as usize + 1];
    for (id, b) in bodies.iter().enumerate() {
        let id = id as i32 + 1;
        match targets(b, "target bodies") {
            Some(t) => {
                for &tag in t {
                    if (1..=max).contains(&tag) && map[tag as usize] == 0 {
                        map[tag as usize] = id;
                    }
                }
            }
            None if map[id as usize] == 0 => map[id as usize] = id,
            None => {}
        }
    }
    let nb = mesh.n_bulk;
    for el in &mut mesh.elements[..nb] {
        el.tag = map[el.tag as usize];
    }
}

/// `meshBC` of `WriteMeshToDisk2`: the BC's first target boundary, or the
/// constraint itself.
fn unmap_bc_tag(tag: i32, bcs: &[ParamMap]) -> i32 {
    match tag {
        t if t > 0 && t as usize <= bcs.len() => {
            targets(&bcs[t as usize - 1], "target boundaries").map_or(t, |b| b[0])
        }
        t => t,
    }
}

/// Fortran `ES17.10`.
pub fn es17(x: f64) -> String {
    let s = format!("{x:.10E}");
    let (m, e) = s.split_once('E').unwrap();
    format!("{m}E{:+03}", e.parse::<i32>().unwrap())
}

/// `mesh.nodes`, `mesh.elements` and `mesh.boundary` as whitespace-separated
/// fields per line.
pub fn written(mesh: &Mesh, bcs: &[ParamMap]) -> [Vec<String>; 3] {
    let nodes = mesh
        .coords
        .iter()
        .enumerate()
        .map(|(i, c)| format!("{} -1 {} {} {}", i + 1, es17(c[0]), es17(c[1]), es17(c[2])))
        .collect();
    let ids = |e: &Element| {
        e.nodes
            .iter()
            .map(|v| (v + 1).to_string())
            .collect::<Vec<_>>()
            .join(" ")
    };
    let elements = mesh
        .bulk()
        .iter()
        .enumerate()
        .map(|(i, e)| format!("{} {} {} {}", i + 1, e.tag, e.etype.0, ids(e)))
        .collect();
    let p = |p: Option<ElemId>| p.map_or(0, |p| p.0 + 1);
    let boundary = mesh
        .boundary()
        .iter()
        .enumerate()
        .map(|(i, e)| {
            let tag = unmap_bc_tag(e.tag, bcs);
            format!(
                "{} {tag} {} {} {} {}",
                i + 1,
                p(e.parents[0]),
                p(e.parents[1]),
                e.etype.0,
                ids(e)
            )
        })
        .collect();
    [nodes, elements, boundary]
}

fn normalized(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().collect::<Vec<_>>().join(" "))
        .filter(|l| !l.is_empty())
        .collect()
}

/// Compares every line the mesh writes with `expected/mesh.*`.
pub fn compare(out: &Mesh, bcs: &[ParamMap], dir: &Path, expected: &Path) {
    let got = written(out, bcs);
    for (file, got) in ["mesh.nodes", "mesh.elements", "mesh.boundary"]
        .iter()
        .zip(got)
    {
        let want = normalized(&expected.join(file));
        assert_eq!(
            got.len(),
            want.len(),
            "{}: {file} line count",
            dir.display()
        );
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g, w, "{}: {file} line {}", dir.display(), i + 1);
        }
    }
    let header = normalized(&expected.join("mesh.header"));
    assert_eq!(
        header[0],
        format!(
            "{} {} {}",
            out.coords.len(),
            out.n_bulk,
            out.elements.len() - out.n_bulk
        ),
        "{}: mesh.header",
        dir.display()
    );
}
