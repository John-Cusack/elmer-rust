//! `StressSolver` on p elements through `elmer_solver::Simulation` with
//! `elmer-pdofs`' `PDofsBackend`: `fem/tests/ElastPelem2dQuad` (`Element =
//! p:4`, full elements) with only its StressSolve solver (Calculate Stresses
//! removed: it does not change the displacement) and `ElastPelem2dTri`
//! (`p:3` triangles), against the pristine runs' Perm, every slot's
//! displacement and the norms 5.99181622E-03 / 7.71519407E-03
//! (`tests/data/pelem2d_quad`, `tests/data/pelem2d_tri`).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, NullLogger, ParamMap, SerialComm, Value};
use elmer_physics_elastic::{install_fe, register, FeServices};
use elmer_solver::{Model, Simulation, SolverRegistry};

fn data(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(case)
}

fn lines(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(str::to_owned).collect())
        .collect()
}

/// One element: type, nodes, tag (body or boundary) and parents.
type Elem = (ElementType, Vec<u32>, i32, [Option<ElemId>; 2]);

/// An Elmer mesh directory in file order, with boundary parents.
struct FileMesh {
    coords: Vec<[f64; 3]>,
    elems: Vec<Elem>,
    n_bulk: usize,
}

impl FileMesh {
    fn read(dir: &Path) -> Self {
        let mut node = HashMap::new();
        let mut coords = Vec::new();
        for f in lines(&dir.join("mesh.nodes")) {
            node.insert(f[0].clone(), coords.len() as u32);
            coords.push([f[2].parse().unwrap(), f[3].parse().unwrap(), f[4].parse().unwrap()]);
        }
        let mut elem = HashMap::new();
        let mut elems = Vec::new();
        for f in lines(&dir.join("mesh.elements")) {
            elem.insert(f[0].clone(), ElemId(elems.len() as u32));
            let nodes = f[3..].iter().map(|s| node[s]).collect();
            elems.push((ElementType(f[2].parse().unwrap()), nodes, f[1].parse().unwrap(), [None, None]));
        }
        let n_bulk = elems.len();
        for f in lines(&dir.join("mesh.boundary")) {
            let parent = |s: &String| if s == "0" { None } else { Some(elem[s]) };
            let nodes = f[5..].iter().map(|s| node[s]).collect();
            elems.push((ElementType(f[4].parse().unwrap()), nodes, f[1].parse().unwrap(), [parent(&f[2]), parent(&f[3])]));
        }
        Self { coords, elems, n_bulk }
    }
}

impl MeshView for FileMesh {
    fn dim(&self) -> usize {
        2
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, n: u32) -> [f64; 3] {
        self.coords[n as usize]
    }
    fn n_bulk(&self) -> usize {
        self.n_bulk
    }
    fn n_boundary(&self) -> usize {
        self.elems.len() - self.n_bulk
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        let (t, nodes, tag, parents) = &self.elems[e.0 as usize];
        ElementRef { etype: *t, nodes, tag: *tag, parents: *parents }
    }
}

/// `elmer-solver`'s view of the parsed case (Elmer lowercases strings).
fn sif_model(path: &Path) -> Model {
    let parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
    let m = parsed.sif.merged();
    let convert = |s: &elmer_sif::Section| -> ParamMap {
        let mut p = elmer_sif::section_to_param_map(s);
        for v in p.0.values_mut() {
            if let Value::Str(t) = v {
                *t = t.to_lowercase();
            }
        }
        p
    };
    let numbered = |v: &[elmer_sif::Section]| -> Vec<ParamMap> {
        let n = v.iter().map(|s| s.id as usize).max().unwrap_or(0);
        let mut out = vec![ParamMap::new(); n];
        for s in v {
            out[s.id as usize - 1] = convert(s);
        }
        out
    };
    let one = |s: &Option<elmer_sif::Section>| s.as_ref().map(convert).unwrap_or_default();
    Model {
        header: ParamMap::new(),
        simulation: one(&m.simulation),
        constants: one(&m.constants),
        solvers: numbered(&m.solvers),
        equations: numbered(&m.equations),
        bodies: numbered(&m.bodies),
        materials: numbered(&m.materials),
        body_forces: numbered(&m.body_forces),
        initial_conditions: numbered(&m.ics),
        boundary_conditions: numbered(&m.bcs),
    }
}

/// Runs `tests/data/<case>/case.sif` and checks the displacement against
/// the pristine `perm.dat` / `solution.dat` and `norm`.
fn check_case(case: &str, n_slots: usize, norm_want: f64, reference: f64) {
    let _ = install_fe(FeServices {
        basis: Box::new(elmer_fe::LagrangeBasis::new(2)),
        quad: Box::new(elmer_fe::ElmerQuadrature),
    });
    let dir = data(case);
    let mesh = FileMesh::read(&dir);
    let model = sif_model(&dir.join("case.sif"));
    let mut reg = SolverRegistry::new();
    register(&mut reg);
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(elmer_pdofs::PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("displacement").unwrap();
    assert_eq!(u.values.len(), 2 * n_slots, "{case}");
    let perm = u.perm.as_deref().unwrap();
    let want_perm: Vec<u32> = lines(&dir.join("perm.dat")).iter().map(|f| f[1].parse().unwrap()).collect();
    let ours: Vec<u32> = perm.iter().map(|&p| p + 1).collect();
    assert_eq!(ours, want_perm, "{case}: Perm");
    let want = lines(&dir.join("solution.dat"));
    let scale = want
        .iter()
        .flat_map(|f| [f[1].parse::<f64>().unwrap(), f[2].parse::<f64>().unwrap()])
        .fold(0.0f64, |m, v| m.max(v.abs()));
    for f in &want {
        let slot: usize = f[0].parse().unwrap();
        let row = perm[slot - 1] as usize;
        for c in 0..2 {
            let (ours, theirs) = (u.values[2 * row + c], f[1 + c].parse::<f64>().unwrap());
            assert!((ours - theirs).abs() <= 1e-10 * scale, "{case}: slot {slot} comp {c}: {ours:e} vs {theirs:e}");
        }
    }
    let norm = u.conv.norm;
    assert!((norm - norm_want).abs() < 1e-12, "{case}: norm {norm:e}");
    assert!((norm - reference).abs() / reference < 1e-5, "{case}");
}

#[test]
fn elast_pelem_2d_quad_matches_pristine() {
    // p:4 on 12 quads, full elements: 225 slots.
    check_case("pelem2d_quad", 225, 0.005991816215121997, 5.99181622e-3);
}

#[test]
fn elast_pelem_2d_tri_matches_pristine() {
    // p:3 on triangles (edge functions of degree 2 and 3, one bubble each).
    check_case("pelem2d_tri", 133, 0.007715194063294897, 7.71519407e-3);
}
