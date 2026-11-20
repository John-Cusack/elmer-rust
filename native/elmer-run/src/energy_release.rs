//! `ReleaseRateSolver` (fem/src/modules/EnergyRelease.F90): the energy
//! release rate of a crack by the domain integral of the propagation shape
//! `theta` (`Shape Displacement`) against the displacement `u` (`True
//! Displacement`) and its stress `S` (`Stress`), over the solver's active
//! elements:
//!
//! `2G = sum_IP (S theta'^T - div(theta) S) : u' + (u' theta') : S`
//!
//! The solver has no variable; it prints `2G` and leaves the fields alone.
//! A field that does not exist reads as zero, as `GetVectorLocalSolution`.

use elmer_core::{BasisEval, ElementBasis, Params};
use elmer_fe::{CoordinateSystem, LagrangeBasis};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

const CALLER: &str = "ReleaseRateSolver";
const PROCEDURE: (&str, &str) = ("EnergyRelease", "ReleaseRateSolver");

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(ReleaseRate)));
}

struct ReleaseRate;

/// `GetVectorLocalSolution(name)`: the field's components (up to `nc`) at
/// the element's nodes, zero where it has none.
fn local(ctx: &SolverCtx<'_>, name: &str, nodes: &[u32], nc: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0; nodes.len()]; nc];
    let Some(v) = ctx.sim.vars.get(name) else { return out };
    for (j, &k) in nodes.iter().enumerate() {
        let Some(s) = v.slot(k) else { continue };
        for (i, row) in out.iter_mut().enumerate().take(v.dofs.min(nc)) {
            row[j] = v.values[s * v.dofs + i];
        }
    }
    out
}

impl Solver for ReleaseRate {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let csys = ctx.simulation_params().string("coordinate system");
        let cs = CoordinateSystem::from_keyword(csys, ctx.mesh.dim())
            .ok_or_else(|| SolverError::fatal("SetCoordinateSystem", "Unknown Coordinate System"))?;
        let basis = LagrangeBasis::new(cs.dim);
        let mut ev = BasisEval::default();
        let active: Vec<_> = ctx.active_elements().collect();
        let mut gtheta = 0.0;
        for e in active {
            let el = ctx.mesh.element(e);
            let n = el.nodes.len();
            let disp = local(ctx, "true displacement", el.nodes, 3);
            let stress = local(ctx, "stress", el.nodes, 6);
            let shape = local(ctx, "shape displacement", el.nodes, 3);
            let xyz: Vec<[f64; 3]> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            let def = elmer_fe::element_type(el.etype.0).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let rule = elmer_fe::quadrature::gauss_points(el.etype, Some(def.gauss_points))
                .map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let mut local_g = 0.0;
            for t in 0..rule.len() {
                basis.eval(el.etype, &xyz, [rule.u[t], rule.v[t], rule.w[t]], &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let (b, db) = (&ev.basis[..n], &ev.dbasis_dx[..n]);
                let mut gu = [[0.0f64; 3]; 3];
                let mut gs = [[0.0f64; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        gu[i][j] = (0..n).map(|k| db[k][j] * disp[i][k]).sum();
                        gs[i][j] = (0..n).map(|k| db[k][j] * shape[i][k]).sum();
                    }
                }
                let div = gs[0][0] + gs[1][1] + gs[2][2];
                let at = |c: usize| -> f64 { (0..n).map(|k| b[k] * stress[c][k]).sum() };
                let (s11, s22, s33, s12, s23, s13) = (at(0), at(1), at(2), at(3), at(4), at(5));
                let st = [[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]];
                let mut r = [[0.0f64; 3]; 3];
                let mut s = [[0.0f64; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        r[i][j] = -(0..3).map(|k| gu[i][k] * gs[k][j]).sum::<f64>();
                        s[i][j] = (0..3).map(|k| st[i][k] * gs[j][k]).sum::<f64>() - div * st[i][j];
                    }
                }
                // SUM over Fortran's column-major order.
                let (mut a, mut c) = (0.0, 0.0);
                for j in 0..3 {
                    for i in 0..3 {
                        a += s[i][j] * gu[i][j];
                    }
                }
                for j in 0..3 {
                    for i in 0..3 {
                        c += r[i][j] * st[i][j];
                    }
                }
                local_g += (a - c) * rule.s[t] * ev.det_j;
            }
            gtheta += local_g;
        }
        println!();
        println!(" *******************************");
        println!(" Calculated Fracture Parameters:");
        println!(" -------------------------------");
        println!(" 2G = {}", elmer_fixtures::udf::list_real(gtheta));
        println!(" *******************************");
        println!();
        Ok(())
    }
}
