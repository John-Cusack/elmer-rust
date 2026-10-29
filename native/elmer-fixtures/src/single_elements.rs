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
//!
//! With a p-element definition (`Element = "p:1"`, ElementUnitTest_p1)
//! `ElementInfo` and `GaussPoints` are those of the p element. The
//! solver's variable is a global, so `elmer-solver` builds no dof layout
//! for it; the solver builds `elmer-pdofs`' layout for its own section and
//! takes each element's basis, rule and node order from it. The Fortran
//! integrates only the nodal functions (`nd = n`).

use std::time::Instant;

use elmer_core::{BasisEval, ElemId, ElementBasis, ElementType, ParamMap, Params, QuadPoint, Quadrature, Value};
use elmer_fe::ElmerQuadrature;
use elmer_pdofs::PDofsBackend;
use elmer_solver::fortran_fmt::list_int;
use elmer_solver::{DofsBackend, Solver, SolverCtx, SolverDofs, SolverError};

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
    /// The p-element layout of the section, built on the first visit.
    p_layout: Option<Box<dyn SolverDofs>>,
}

impl SingleElements {
    pub fn new(out: PrintSink) -> Self {
        Self { no_test: 0, no_mode: 0, cumct: 0.0, out, p_layout: None }
    }
}

/// A section whose `Element` is not the nodal `n:1`.
fn is_p_section(params: &ParamMap) -> bool {
    params.string("element").is_some_and(|e| elmer_core::normalize_key(e) != "n:1")
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(SingleElements::new(stdout_sink())))
}

/// `LocalMatrix`: `(SUM(MASS), SUM(ABS(STIFF)))` of one element, both summed
/// in Fortran's column-major order.
pub fn element_sums(
    etype: ElementType,
    xyz: &[[f64; 3]],
    basis: &dyn ElementBasis,
    rule: &[QuadPoint],
) -> Result<(f64, f64), String> {
    let nd = xyz.len();
    let mut mass = vec![0.0; nd * nd];
    let mut stiff = vec![0.0; nd * nd];
    let mut ev = BasisEval::default();
    for qp in rule {
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
        crate::hooks::check_element_definitions(ctx)?;
        if ctx.primary().is_none() {
            return Err(SolverError::NotWired(
                "SingleElements: its _init variable (elmer_fixtures::registry::init_sections before setup)".into(),
            ));
        }
        // Built here, so an element `elmer-pdofs` does not handle is refused
        // before the run.
        if is_p_section(ctx.params()) {
            let active: Vec<ElemId> = ctx.active_elements().collect();
            let layout = PDofsBackend::new()
                .for_solver(ctx.mesh, &ctx.sim.model, ctx.params(), &active, 1)
                .ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: element definition without p dofs")))??;
            self.p_layout = Some(layout);
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
        let fatal = |m: String| SolverError::fatal(CALLER, m);
        let (xyz, basis, rule): (_, Box<dyn ElementBasis>, _) = match &self.p_layout {
            Some(d) => {
                let xyz: Vec<[f64; 3]> = d.element_nodes(e).iter().map(|&n| ctx.mesh.coord(n)).collect();
                (xyz, d.element_basis(e).map_err(fatal)?, d.element_rule(e).map_err(fatal)?)
            }
            None => (element_nodes(ctx.mesh, e), Box::new(lagrange_basis(ctx)?), ElmerQuadrature.gauss_points(etype, None)),
        };

        let start = Instant::now();
        let mut sums = (0.0, 0.0);
        // `DO i=0,NoCycles` repeats ElementInfo for timing; the result is
        // the same each time.
        for _ in 0..=cycles {
            sums = element_sums(etype, &xyz, basis.as_ref(), &rule).map_err(fatal)?;
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
