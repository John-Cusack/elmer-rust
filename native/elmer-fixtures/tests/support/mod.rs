//! Test-only stand-ins (ADR 003, rule 6) until `elmer-fe` and the solver
//! framework's live model view are available: a P1 triangle / two-node line
//! basis with its quadrature, and a mutable model whose `Time`, coordinates
//! and `Potential` field the tests drive step by step.
#![allow(dead_code)]

use std::sync::Mutex;

use elmer_core::{BasisEval, ElementBasis, ElementType, Logger, QuadPoint, Quadrature};
use elmer_expr::{FieldLookup, MatcEngine};
use elmer_fixtures::ModelView;

pub const TRI: ElementType = ElementType(303);
pub const LINE: ElementType = ElementType(202);

/// Exact for polynomials of degree 2 on the reference triangle and degree 3
/// on the reference line, enough for P1 stiffness and mass.
pub struct TestQuadrature;

impl Quadrature for TestQuadrature {
    fn gauss_points(&self, etype: ElementType, _np: Option<usize>) -> Vec<QuadPoint> {
        match etype {
            TRI => [[0.5, 0.0], [0.5, 0.5], [0.0, 0.5]]
                .iter()
                .map(|&[u, v]| QuadPoint { uvw: [u, v, 0.0], weight: 1.0 / 6.0 })
                .collect(),
            LINE => {
                let g = 1.0 / 3.0_f64.sqrt();
                [-g, g].iter().map(|&u| QuadPoint { uvw: [u, 0.0, 0.0], weight: 1.0 }).collect()
            }
            other => panic!("no test quadrature for {other:?}"),
        }
    }
}

/// P1 triangle (reference `(0,0) (1,0) (0,1)`) and the two-node line
/// (reference `[-1, 1]`), both embedded in the x-y plane.
pub struct TestBasis;

impl ElementBasis for TestBasis {
    fn eval(
        &self,
        etype: ElementType,
        xyz: &[[f64; 3]],
        uvw: [f64; 3],
        out: &mut BasisEval,
    ) -> Result<(), String> {
        match etype {
            TRI => {
                let (u, v) = (uvw[0], uvw[1]);
                out.basis = vec![1.0 - u - v, u, v];
                let (a, b, c) = (xyz[0], xyz[1], xyz[2]);
                let j = [[b[0] - a[0], c[0] - a[0]], [b[1] - a[1], c[1] - a[1]]];
                let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
                if det == 0.0 {
                    return Err("degenerate triangle".into());
                }
                // grad(phi) = J^-T grad_ref(phi)
                let inv_t = [[j[1][1] / det, -j[1][0] / det], [-j[0][1] / det, j[0][0] / det]];
                let g_ref = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
                out.dbasis_dx = g_ref
                    .iter()
                    .map(|g| [inv_t[0][0] * g[0] + inv_t[0][1] * g[1], inv_t[1][0] * g[0] + inv_t[1][1] * g[1], 0.0])
                    .collect();
                out.det_j = det;
                Ok(())
            }
            LINE => {
                let u = uvw[0];
                out.basis = vec![(1.0 - u) / 2.0, (1.0 + u) / 2.0];
                let (a, b) = (xyz[0], xyz[1]);
                let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
                let dir = [(b[0] - a[0]) / len, (b[1] - a[1]) / len, 0.0];
                // d phi / ds = (+-1/2) / (len/2) along the element.
                let d = 1.0 / len;
                out.dbasis_dx = vec![dir.map(|c| -c * d), dir.map(|c| c * d)];
                out.det_j = len / 2.0;
                Ok(())
            }
            other => Err(format!("no test basis for {other:?}")),
        }
    }
}

/// Records `info` calls.
#[derive(Default)]
pub struct RecordingLogger(pub Mutex<Vec<(String, String, u8)>>);

impl Logger for RecordingLogger {
    fn info(&self, caller: &str, msg: &str, level: u8) {
        self.0.lock().unwrap().push((caller.into(), msg.into(), level));
    }
    fn warn(&self, _caller: &str, _msg: &str) {}
}

#[derive(Default)]
struct State {
    time: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    potential: Vec<Option<f64>>,
}

/// A model the test mutates between calls, like the solver does between steps.
pub struct TestModel(Mutex<State>);

impl TestModel {
    pub fn line(xs: &[f64]) -> Self {
        Self(Mutex::new(State {
            time: 0.0,
            x: xs.to_vec(),
            y: vec![0.0; xs.len()],
            potential: vec![Some(0.0); xs.len()],
        }))
    }

    pub fn points(xy: &[(f64, f64)]) -> Self {
        Self(Mutex::new(State {
            time: 0.0,
            x: xy.iter().map(|p| p.0).collect(),
            y: xy.iter().map(|p| p.1).collect(),
            potential: vec![None; xy.len()],
        }))
    }

    pub fn set_time(&self, t: f64) {
        self.0.lock().unwrap().time = t;
    }

    pub fn set_potential(&self, f: impl Fn(f64) -> f64) {
        let mut s = self.0.lock().unwrap();
        s.potential = s.x.iter().map(|&x| Some(f(x))).collect();
    }

    pub fn clear_potential_at(&self, node: usize) {
        self.0.lock().unwrap().potential[node] = None;
    }
}

impl FieldLookup for TestModel {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        let s = self.0.lock().unwrap();
        match var {
            "potential" => s.potential.get(node as usize).copied().flatten(),
            _ => None,
        }
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        let s = self.0.lock().unwrap();
        [s.x[node as usize], s.y[node as usize], 0.0]
    }
    fn time(&self) -> f64 {
        self.0.lock().unwrap().time
    }
}

impl ModelView for TestModel {
    fn n_nodes(&self) -> usize {
        self.0.lock().unwrap().x.len()
    }
}

/// MATC is not used by the fixtures.
pub struct NoMatc;

impl MatcEngine for NoMatc {
    fn eval(&self, _expr: &str, _args: &[f64]) -> Result<Vec<f64>, String> {
        Err("no MATC in these tests".into())
    }
}
