//! `ElementUnitTest_edge*/SingleEdgeElements.F90`: the H(curl) analogue of
//! [`crate::single_elements`]. Each visit takes the next element of the
//! mesh (bulk then boundary) and stores `SUM(MASS)` of its edge basis
//! (`Integral W_p . W_q`) or, on the second pass, `SUM(ABS(STIFF))` of the
//! curls, in the hidden global `dummy_var` that the case's `Reference Norm`
//! checks per `RUN` block.
//!
//! `SingleEdgeElements_init` picks the edge element from two keywords when
//! the section has no `Element`: `Quadratic Approximation` (second-order
//! Piola, `n:0 e:2 ...`), else `Use Piola Transform` (lowest-order Piola,
//! `n:0 e:1 -brick b:3 -quad_face b:2`), else the Whitney basis (`n:0 e:1`).
//! The only reader of that `Element` keyword is the plugin's own
//! `ElementInfo` call, so here [`init_section`] sets the variable and the
//! solver keeps the choice ([`EdgeKind`]); a section that names its own
//! `Element` is refused. The bases are `elmer-fe-hp`'s `get_edge_basis`
//! (Whitney) and `edge_element_info` (Piola).

use std::time::Instant;

use elmer_core::{ElemId, ElementType, ParamMap, Params, Quadrature};
use elmer_fe::ElmerQuadrature;
use elmer_fe_hp::edge::{edge_element_info, get_edge_basis};
use elmer_solver::fortran_fmt::list_int;
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::{coordinate_system, element_nodes};
use crate::single_elements;
use crate::udf::{list_real, stdout_sink, PrintSink};

const CALLER: &str = "SingleEdgeElements";

/// The edge element `SingleEdgeElements_init` selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeKind {
    /// `n:0 e:1`: `GetEdgeBasis`.
    Whitney,
    /// `n:0 e:1 -brick b:3 -quad_face b:2` (degree 1) or
    /// `n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2`
    /// (degree 2): `EdgeElementInfo` with the Piola transform.
    Piola(u8),
}

impl EdgeKind {
    /// The `_init` choice from the section's keywords.
    pub fn from_section(params: &ParamMap) -> Self {
        if params.logical("quadratic approximation") == Some(true) {
            EdgeKind::Piola(2)
        } else if params.logical("use piola transform") == Some(true) {
            EdgeKind::Piola(1)
        } else {
            EdgeKind::Whitney
        }
    }
}

/// `SingleEdgeElements_init`'s `Variable = -nooutput -global dummy_var`.
pub fn init_section(params: &mut ParamMap) {
    single_elements::init_section(params);
}

pub struct SingleEdgeElements {
    kind: EdgeKind,
    no_test: usize,
    no_mode: usize,
    cumct: f64,
    out: PrintSink,
}

impl SingleEdgeElements {
    pub fn new(kind: EdgeKind, out: PrintSink) -> Self {
        Self { kind, no_test: 0, no_mode: 0, cumct: 0.0, out }
    }
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    if let Some(e) = params.string("element") {
        return Err(SolverError::Unsupported(format!(
            "SingleEdgeElements: Element = {e} given in the section (only the _init choices are ported)"
        )));
    }
    Ok(Box::new(SingleEdgeElements::new(EdgeKind::from_section(params), stdout_sink())))
}

/// `LocalMatrix`: `(SUM(MASS), SUM(ABS(STIFF)))` of one element's edge
/// basis over its default Gauss rule, summed in column-major order.
/// `gindex` are the global node numbers that orient the edges.
pub fn edge_sums(
    etype: ElementType,
    xyz: &[[f64; 3]],
    gindex: &[i64],
    kind: EdgeKind,
    cdim: usize,
) -> Result<(f64, f64), String> {
    let mut mass: Vec<f64> = Vec::new();
    let mut stiff: Vec<f64> = Vec::new();
    let mut nd = 0;
    for qp in ElmerQuadrature.gauss_points(etype, None) {
        let ev = match kind {
            EdgeKind::Whitney => get_edge_basis(etype, xyz, gindex, qp.uvw, cdim)?,
            EdgeKind::Piola(degree) => edge_element_info(etype, xyz, gindex, qp.uvw, degree, cdim)?,
        };
        if nd == 0 {
            nd = ev.dofs;
            mass = vec![0.0; nd * nd];
            stiff = vec![0.0; nd * nd];
        }
        let weight = qp.weight * ev.det;
        let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        for q in 0..nd {
            for p in 0..nd {
                stiff[q * nd + p] += weight * dot(&ev.curl[p], &ev.curl[q]);
                mass[q * nd + p] += weight * dot(&ev.edge[p], &ev.edge[q]);
            }
        }
    }
    let sum_mass = mass.iter().fold(0.0, |s, v| s + v);
    let sum_stiff = stiff.iter().fold(0.0, |s, v| s + v.abs());
    Ok((sum_mass, sum_stiff))
}

impl Solver for SingleEdgeElements {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        if ctx.primary().is_none() {
            return Err(SolverError::NotWired(
                "SingleEdgeElements: its _init variable (elmer_fixtures::registry::init_sections before setup)".into(),
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
        let el = ctx.mesh.element(e);
        let etype = el.etype;
        // `Element % NodeIndexes`, 1-based global numbers.
        let gindex: Vec<i64> = el.nodes.iter().map(|&n| i64::from(n) + 1).collect();
        (self.out)(&format!(" Element type: {}", list_int(etype.0 as i32)));
        let xyz = element_nodes(ctx.mesh, e);
        let cdim = coordinate_system(ctx)?.dim;

        let start = Instant::now();
        let mut sums = (0.0, 0.0);
        for _ in 0..=cycles {
            sums = edge_sums(etype, &xyz, &gindex, self.kind, cdim).map_err(|m| SolverError::fatal(CALLER, m))?;
        }
        if cycles > 0 {
            let dt = start.elapsed().as_secs_f64();
            (self.out)(&format!(" Time for element: {}{}{}", list_int(etype.0 as i32), list_real(dt), list_real(dt)));
            self.cumct += dt;
            if self.no_test == n_elements {
                (self.out)(&format!(" Time in total:{}{}", list_real(self.cumct), list_real(self.cumct)));
            }
        }
        // `MODULO(NoMode, 2) == 0`: mass on even passes, stiffness on odd.
        let norm = if self.no_mode.is_multiple_of(2) { sums.0 } else { sums.1 };
        (self.out)(&format!(" Norm for element: {}", list_real(norm)));

        let var = ctx.primary_mut().ok_or_else(|| SolverError::fatal(CALLER, "no solver variable"))?;
        var.values[0] = norm;
        Ok(())
    }
}
