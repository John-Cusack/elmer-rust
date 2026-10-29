//! `ElementUnitTest_h1/SingleElements.F90`: an element unit test as a
//! solver. `SingleElements_init` makes the solver's variable a hidden global
//! (`Variable = -nooutput -global dummy_var`); each visit of
//! `SingleElements` takes the next element of the mesh (`Mesh % Elements`,
//! bulk then boundary), integrates its mass and stiffness matrices with
//! `ElementInfo` over `GaussPointsAdapt` and stores one number in the
//! variable, which the case's `Reference Norm` checks per `RUN` block.
//!
//! Visits `1..=n` (n elements) give `SUM(MASS)`, the element measure; from
//! then on the counter starts over with `SUM(ABS(STIFF))`. The counters
//! (`NoTest`, `NoMode`, `SAVE`d) live in the instance, and the
//! `Element type:` / `Norm for element:` records go to the instance's
//! [`PrintSink`].

use std::time::Instant;

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, ParamMap, Params, Quadrature, Value};
use elmer_fe::ElmerQuadrature;
use elmer_solver::fortran_fmt::list_int;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{element_nodes, lagrange_basis};
use crate::udf::{list_real, stdout_sink, PrintSink};

const CALLER: &str = "SingleElements";

/// `SingleElements_init`: `ListAddString(Solver % Values, 'Variable',
/// '-nooutput -global dummy_var')`, applied to the section before its
/// variable is created.
pub fn init_section(params: &mut ParamMap) {
    params.insert("Variable", Value::Str("-nooutput -global dummy_var".into()));
}

pub struct SingleElements {
    no_test: usize,
    no_mode: usize,
    cumct: f64,
    out: PrintSink,
}

impl SingleElements {
    pub fn new(out: PrintSink) -> Self {
        Self { no_test: 0, no_mode: 0, cumct: 0.0, out }
    }
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if params.string("element").is_some_and(|e| elmer_core::normalize_key(e) != "n:1") {
        return Err(SolverError::Unsupported("SingleElements: p-, bubble or edge elements".into()));
    }
    Ok(Box::new(SingleElements::new(stdout_sink())))
}

/// `LocalMatrix`: `(SUM(MASS), SUM(ABS(STIFF)))` of one element, both summed
/// in Fortran's column-major order.
pub fn element_sums(
    etype: ElementType,
    xyz: &[[f64; 3]],
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
) -> Result<(f64, f64), String> {
    let nd = xyz.len();
    let mut mass = vec![0.0; nd * nd];
    let mut stiff = vec![0.0; nd * nd];
    let mut ev = BasisEval::default();
    for qp in quad.gauss_points(etype, None) {
        basis.eval(etype, xyz, qp.uvw, &mut ev)?;
        if ev.basis.len() < nd {
            return Err(format!("ElementInfo gave {} basis functions for {nd} nodes", ev.basis.len()));
        }
        let weight = qp.weight * ev.det_j;
        let (b, g) = (&ev.basis, &ev.dbasis_dx);
        // Column-major (p fastest), as the Fortran arrays are laid out.
        for q in 0..nd {
            for p in 0..nd {
                let m = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                stiff[q * nd + p] += weight * m;
                mass[q * nd + p] += weight * b[q] * b[p];
            }
        }
    }
    let sum_mass = mass.iter().fold(0.0, |s, v| s + v);
    let sum_stiff = stiff.iter().fold(0.0, |s, v| s + v.abs());
    Ok((sum_mass, sum_stiff))
}

impl Solver for SingleElements {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none() {
            return Err(SolverError::NotWired(
                "SingleElements: its _init variable (elmer_fixtures::registry::init_sections before setup)".into(),
            ));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let n_elements = ctx.mesh.n_bulk() + ctx.mesh.n_boundary();
        self.no_test += 1;
        if self.no_test > n_elements {
            self.no_mode += 1;
            self.no_test = 1;
        }
        let cycles = if self.no_mode == 0 { ctx.params().integer("timing cycles").unwrap_or(0).max(0) } else { 0 };

        let e = ElemId(self.no_test as u32 - 1);
        let etype = ctx.mesh.element(e).etype;
        (self.out)(&format!(" Element type: {}", list_int(etype.0 as i32)));
        let xyz = element_nodes(ctx.mesh, e);
        let basis = lagrange_basis(ctx)?;

        let start = Instant::now();
        let mut sums = (0.0, 0.0);
        // `DO i=0,NoCycles` repeats ElementInfo for timing; the result is
        // the same each time.
        for _ in 0..=cycles {
            sums = element_sums(etype, &xyz, &basis, &ElmerQuadrature).map_err(|m| SolverError::fatal(CALLER, m))?;
        }
        if cycles > 0 {
            let dt = start.elapsed().as_secs_f64();
            (self.out)(&format!(" Time for element: {}{}{}", list_int(etype.0 as i32), list_real(dt), list_real(dt)));
            self.cumct += dt;
            if self.no_test == n_elements {
                (self.out)(&format!(" Time in total:{}{}", list_real(self.cumct), list_real(self.cumct)));
            }
        }
        let norm = if self.no_mode == 0 { sums.0 } else { sums.1 };
        (self.out)(&format!(" Norm for element: {}", list_real(norm)));

        let var = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
        var.values[0] = norm;
        Ok(())
    }
}
