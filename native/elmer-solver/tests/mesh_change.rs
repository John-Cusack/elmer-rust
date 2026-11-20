//! Moving a steady run to a new mesh (`elmer_solver::remesh`), as adaptive
//! mesh refinement does: a solver asks for a finer mesh at the end of its
//! visit; the run saves the mesh the solution lives on, interpolates every
//! nodal field to the new mesh and goes on there with the solver's state
//! kept; `Mesh % AdaptiveFinished` ends the steady iterations.

mod support;

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use elmer_core::{ElemId, ElementRef, ElementType, MeshView, ParamMap, SerialComm};
use elmer_solver::remesh::MeshChange;
use elmer_solver::sim::SaveCurrent;
use elmer_solver::{FeBackend, Model, Simulation, Solver, SolverCtx, SolverError, SolverRegistry};
use support::{i, s, sec, Capture};

/// The unit square in `n x n` bilinear quadrilaterals, one body.
struct Grid {
    n: usize,
    coords: Vec<[f64; 3]>,
    elems: Vec<[u32; 4]>,
}

impl Grid {
    fn new(n: usize) -> Grid {
        let mut coords = Vec::new();
        for j in 0..=n {
            for k in 0..=n {
                coords.push([k as f64 / n as f64, j as f64 / n as f64, 0.0]);
            }
        }
        let id = |k: usize, j: usize| (j * (n + 1) + k) as u32;
        let mut elems = Vec::new();
        for j in 0..n {
            for k in 0..n {
                elems.push([id(k, j), id(k + 1, j), id(k + 1, j + 1), id(k, j + 1)]);
            }
        }
        Grid { n, coords, elems }
    }
}

impl MeshView for Grid {
    fn dim(&self) -> usize {
        2
    }
    fn n_nodes(&self) -> usize {
        self.coords.len()
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        self.coords[node as usize]
    }
    fn n_bulk(&self) -> usize {
        self.elems.len()
    }
    fn n_boundary(&self) -> usize {
        0
    }
    fn element(&self, e: ElemId) -> ElementRef<'_> {
        ElementRef {
            etype: ElementType(404),
            nodes: &self.elems[e.0 as usize],
            tag: 1,
            parents: [None, None],
        }
    }
}

/// What the toy solver saw on each visit: the mesh size, and the largest
/// deviation of the field it found from `x + 2y`.
type Seen = Arc<Mutex<Vec<(usize, f64)>>>;

/// Writes `U = k (x + 2y)` on its `k`th visit; asks for a finer grid on its first two visits and
/// declares the mesh finished on the third.
struct Toy {
    seen: Seen,
    visits: usize,
}

fn field(p: [f64; 3]) -> f64 {
    p[0] + 2.0 * p[1]
}

impl Solver for Toy {
    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        self.visits += 1;
        let n = ctx.mesh.n_nodes();
        let coords: Vec<[f64; 3]> = (0..n as u32).map(|k| ctx.mesh.coord(k)).collect();
        let v = ctx.primary_mut().expect("U");
        let perm = v.perm.expect("nodal").to_vec();
        let mut dev = 0.0f64;
        for (node, &slot) in perm.iter().enumerate() {
            if slot == elmer_core::NO_DOF {
                continue;
            }
            let k = self.visits as f64;
            if self.visits > 1 {
                dev = dev.max((v.values[slot as usize] - (k - 1.0) * field(coords[node])).abs());
            }
            v.values[slot as usize] = k * field(coords[node]);
        }
        // A field only this mesh has (as RefineMesh's `.error` and
        // `Hvalue`): added on the first visit, interpolated to the next.
        if self.visits == 1 {
            let values = coords.iter().map(|p| 3.0 * p[0]).collect();
            ctx.sim.vars.add(elmer_solver::Variable::new(
                "probe",
                elmer_core::FieldKind::Nodal,
                1,
                None,
                values,
            ))?;
        } else if self.visits == 2 {
            let p = ctx
                .sim
                .vars
                .get("probe")
                .expect("interpolated to the new mesh");
            assert_eq!(p.values.len(), n);
            for (node, c) in coords.iter().enumerate() {
                dev = dev.max((p.values[node] - 3.0 * c[0]).abs());
            }
        }
        self.seen.lock().unwrap().push((n, dev));
        if self.visits < 3 {
            let grid = Grid::new(2 * (ctx.adaptive_depth() + 1) + 2);
            ctx.request_mesh_change(MeshChange {
                mesh: Arc::new(grid),
                name: format!("RefinedMesh{}", ctx.adaptive_depth() + 1),
            });
        } else {
            ctx.set_adaptive_finished(true);
        }
        Ok(())
    }
}

/// Records the mesh each save was made on.
struct Saves(Rc<RefCell<Vec<(String, usize)>>>, String);

impl SaveCurrent for Saves {
    fn save(
        &mut self,
        _sim: &Simulation,
        mesh: &dyn MeshView,
        _log: &dyn elmer_core::Logger,
    ) -> Result<(), SolverError> {
        self.0.borrow_mut().push((self.1.clone(), mesh.n_nodes()));
        Ok(())
    }
}

fn model() -> Model {
    Model {
        simulation: sec(&[
            ("Simulation Type", s("Steady State")),
            ("Steady State Max Iterations", i(&[10])),
            ("Coordinate System", s("Cartesian 2D")),
        ]),
        solvers: vec![sec(&[
            ("Equation", s("toy")),
            ("Procedure", s("Toy Toy")),
            ("Variable", s("U")),
            ("Steady State Convergence Tolerance", support::r(1.0e-12)),
        ])],
        equations: vec![sec(&[("Active Solvers", i(&[1]))])],
        bodies: vec![sec(&[("Target Bodies", i(&[1])), ("Equation", i(&[1]))])],
        ..Default::default()
    }
}

#[test]
fn the_run_follows_the_refined_meshes() {
    let seen: Seen = Arc::new(Mutex::new(Vec::new()));
    let saves: Rc<RefCell<Vec<(String, usize)>>> = Rc::new(RefCell::new(Vec::new()));
    let seen2 = seen.clone();
    thread_local! {
        static SEEN: RefCell<Option<Seen>> = const { RefCell::new(None) };
    }
    SEEN.with(|c| *c.borrow_mut() = Some(seen2));
    fn toy(_p: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
        let seen = SEEN.with(|c| c.borrow().clone()).expect("set");
        Ok(Box::new(Toy { seen, visits: 0 }))
    }
    let mut reg = SolverRegistry::new();
    reg.register("Toy", "Toy", toy);
    let log = Capture::default();
    let mut sim = Simulation::new(model());
    sim.backends.fe = Some(FeBackend {
        basis: Box::new(support::basis::TestBasis),
        quadrature: Box::new(support::quad::TestQuadrature),
    });
    let first = Grid::new(2);
    sim.setup(&first, &reg, &log, &SerialComm).unwrap();
    sim.set_save_current(Box::new(Saves(saves.clone(), "first".into())));
    let s2 = saves.clone();
    sim.set_mesh_writer_factory(Box::new(move |name: &str, _sim: &Simulation| {
        Ok(Some(
            Box::new(Saves(s2.clone(), name.to_string())) as Box<dyn SaveCurrent>
        ))
    }));
    let report = sim.run(&first, &log, &SerialComm).unwrap();

    // Three visits on the 2x2, 4x4 and 6x6 grids; the field the solver
    // found on the new grids is the interpolated (k - 1)(x + 2y), exact for
    // the bilinear elements.
    let seen = seen.lock().unwrap().clone();
    assert_eq!(
        seen.iter().map(|s| s.0).collect::<Vec<_>>(),
        vec![9, 25, 49]
    );
    for (_, dev) in &seen {
        assert!(*dev < 1e-14, "{seen:?}");
    }
    // Each iteration saved the mesh it solved on, each mesh to its writer.
    assert_eq!(
        *saves.borrow(),
        vec![
            ("first".to_string(), 9),
            ("RefinedMesh1".to_string(), 25),
            ("RefinedMesh2".to_string(), 49)
        ]
    );
    // AdaptiveFinished ended the steady iterations on the third.
    assert!(report.steady_state_reached);
    assert_eq!(report.steps, 3);
    assert_eq!(sim.adaptive_depth(), 2);
    assert!(log.lines(&["SolveEquations: Overriding convergence"]).len() == 1);
    assert_eq!(sim.vars.get("u").unwrap().values.len(), 49);
    let _ = (Grid::new(1).n,);
}
