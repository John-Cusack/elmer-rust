//! Scalar `Procedure` fixtures: functions with Elmer's `ExecRealFunction`
//! shape `FUNCTION F(Model, n, t)`.
//!
//! Each fixture has a plain Rust core returning `Result` and an adapter
//! that fits [`elmer_expr::NodeFn`]. The adapter turns a Fortran `Fatal`
//! into a panic: `NodeFn` has no error channel, and Elmer halts the run at
//! that point too.

use std::sync::{Arc, Mutex, PoisonError};

use elmer_expr::NodeFn;

use crate::{FixtureError, ModelView};

/// `ConstantParamFunc/MySource.F90`: `f = SIN(tx(1)+tx(2))`.
///
/// The case calls it as `Variable "time, 3"`, so `tx(1)` is time and `tx(2)`
/// the constant parameter. The function reads neither `Model` nor `n`.
pub fn my_source(tx: [f64; 2]) -> f64 {
    (tx[0] + tx[1]).sin()
}

/// [`my_source`] as a registry function. The node is ignored, so it also
/// works for integration-point evaluation (`u32::MAX`).
pub fn my_source_fn() -> NodeFn {
    Arc::new(|_node, tx| match tx {
        [t, c, ..] => my_source([*t, *c]),
        _ => panic!("MySource: needs two arguments (time and a constant), got {}", tx.len()),
    })
}

/// The squared cylinder radius of `AdvDiffFCT/InitDist.f90` as the Fortran
/// compares it. `0.3**2` is written without a `_dp` suffix, so gfortran folds
/// it in single precision and widens the result for the comparison with a
/// double; the pristine `InitDist.so` stores `0x3FB70A3D80000000`, not the
/// double `0.09`. Boundary points depend on the difference.
pub const INIT_DIST_RADIUS_SQ: f64 = (0.3_f32 * 0.3_f32) as f64;

/// `AdvDiffFCT/InitDist.f90` at a point: `1` strictly inside the cylinder of
/// radius `0.3` around `(0.5, 0)`, else `0`. Points on the rim give `0`.
pub fn init_dist(x: f64, y: f64) -> f64 {
    let dx = x - 0.5;
    if dx * dx + y * y < INIT_DIST_RADIUS_SQ {
        1.0
    } else {
        0.0
    }
}

/// `InitDist` bound to a model: the Fortran reads `Model % Nodes % x(n)` and
/// `y(n)`; the time argument is unused.
///
/// The Fortran also `SAVE`s `Visited`, `a1` and `a2` and sets them on the
/// first call. Nothing reads them, so no state is kept here.
pub fn init_dist_fn(model: Arc<dyn ModelView>) -> NodeFn {
    Arc::new(move |node, _t| {
        if node == u32::MAX {
            panic!("InitDist: called at an integration point; the Fortran indexes Nodes % x(-gauss point)");
        }
        let c = model.coord(node);
        init_dist(c[0], c[1])
    })
}

/// Tolerances `1dtests/Source.f90` and `1dtests_MeshLevels/Source.f90`
/// check the polynomial steps against; they are the only difference between
/// the two sources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OneDTolerances {
    /// Step 3 (parabola) and step 4 (cubic).
    pub polynomial: f64,
    /// Step 5 (sine).
    pub sine: f64,
}

impl OneDTolerances {
    /// `1dtests`.
    pub const ONE_D_TESTS: Self = Self { polynomial: 1.0e-10, sine: 1.0e-4 };
    /// `1dtests_MeshLevels`, whose multigrid-level mesh solves less exactly.
    pub const MESH_LEVELS: Self = Self { polynomial: 1.0e-8, sine: 1.0e-4 };
}

/// `1dtests/Source.f90`: the right-hand side of a five-step Poisson check,
/// with the Fortran's checks of the previous step's solution.
///
/// One instance is one Fortran process: `OldStep` (`INTEGER :: OldStep = -1`,
/// implicitly `SAVE`d) lives in the instance, so two instances never share
/// it.
pub struct OneDSource {
    model: Arc<dyn ModelView>,
    tol: OneDTolerances,
    old_step: Mutex<i64>,
}

impl OneDSource {
    pub fn new(model: Arc<dyn ModelView>, tol: OneDTolerances) -> Self {
        Self { model, tol, old_step: Mutex::new(-1) }
    }

    /// The step the last successful check ran for (`OldStep`; `-1` before).
    pub fn old_step(&self) -> i64 {
        *self.old_step.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// `FUNCTION Source(Model, n, f)`: `node` is 0-based, `f` the value of
    /// `Variable Potential` there.
    ///
    /// On steps 2 to 5 the first call of a step checks the solution of the
    /// previous step at every node against its analytic form and fails with
    /// the Fortran's `Fatal` text if it is off. `OldStep` moves only after a
    /// clean check, so a failed check fails again on the next call.
    pub fn eval(&self, node: u32, f: f64) -> Result<f64, FixtureError> {
        const CALLER: &str = "Source";
        if node == u32::MAX {
            return Err(FixtureError::Undefined {
                caller: CALLER,
                what: "called at an integration point; the Fortran indexes Nodes % x(-gauss point)".into(),
            });
        }
        // NINT rounds half away from zero, as f64::round does.
        let step = self.model.time().round() as i64;

        let mut old_step = self.old_step.lock().unwrap_or_else(PoisonError::into_inner);
        if step > 1 && step != *old_step {
            for i in 0..self.model.n_nodes() as u32 {
                let x = self.model.coord(i)[0];
                let s = self.model.value("potential", i).ok_or_else(|| FixtureError::Undefined {
                    caller: CALLER,
                    what: format!("node {} has no Potential value; the Fortran reads P % Values(0)", i + 1),
                })?;
                check_previous_solution(step, x, s, self.tol)?;
            }
            *old_step = step;
        }

        let x = self.model.coord(node)[0];
        match step {
            1 => Ok(0.0),
            2 => Ok(1.0),
            3 => Ok(x),
            4 | 5 => Ok(f),
            _ => Err(FixtureError::Undefined {
                caller: CALLER,
                what: format!("step {step} is outside 1..=5 and the Fortran leaves h unset"),
            }),
        }
    }

    /// Adapter for the registry; see the module note on `Fatal`.
    pub fn into_fn(self) -> NodeFn {
        Arc::new(move |node, deps| {
            let f = deps.first().copied().unwrap_or_else(|| {
                panic!("Source: needs the Potential value as its argument (Variable Potential)")
            });
            self.eval(node, f).unwrap_or_else(|e| panic!("{e}"))
        })
    }
}

/// The `SELECT CASE(Step)` of the previous-result check at one node.
fn check_previous_solution(step: i64, x: f64, s: f64, tol: OneDTolerances) -> Result<(), FixtureError> {
    use std::f64::consts::PI;
    let fatal = |message| Err(FixtureError::Fatal { caller: "Source", message });
    match step {
        2 if s != 0.0 => fatal("Not able to integrate a zero field ?"),
        3 if (s + (x * x - PI * x) / 2.0).abs() > tol.polynomial => fatal("Not able to integrate parabola ?"),
        4 if (s + (x * x * x - PI * PI * x) / 6.0).abs() > tol.polynomial => {
            fatal("Not able to integrate 3d deg polynomial?")
        }
        5 if (s - 2.0 * x.sin()).abs() > tol.sine => fatal("Not able to integrate sine?"),
        _ => Ok(()),
    }
}
