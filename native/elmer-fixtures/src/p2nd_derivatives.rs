//! `P2ndDerivatives/P2ndDerivativesTest.f90` (and its `SD_` copy): a unit
//! test of the p basis' derivatives as a solver. For each active element,
//! on its reference element's nodal coordinates (`NodeU/V/W`, the p dofs'
//! coordinates zero): the L2 projections of the complete polynomial of
//! degree `qp` (the element's p) and of the basis' first derivatives are
//! formed on `GaussPoints(Element)`; then at each point the projection's
//! value, first and second derivatives (`ElementInfo(.., ddBasisddx,
//! .TRUE.)`) and the derivatives of the projected first derivatives are
//! checked against the polynomial (`CheckValue`, relative tolerances), and
//! the projected derivatives against `ddBasisddx` (1e-8). A failed check
//! stops the run, as the Fortran's `STOP`.

use elmer_core::{BasisEval, ElemId, ParamMap, Params};
use elmer_pdofs::PDofsBackend;
use elmer_solver::{DofsBackend, Solver, SolverCtx, SolverDofs, SolverError};

const CALLER: &str = "TestSolver";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(TestSolver { layout: None }))
}

/// The solver has no variable (`No Matrix`), so it lays its p dofs out
/// itself, as `elmer-pdofs` would for its section.
struct TestSolver {
    layout: Option<Box<dyn SolverDofs>>,
}

/// `CheckValue(f1, f2, eps, str)`.
fn check(f1: f64, f2: f64, eps: f64, what: &str) -> Result<(), SolverError> {
    if (f1 - f2).abs() < 1.0e-12 {
        return Ok(());
    }
    // `REAL :: scal`: single precision.
    let scal = f1.abs().max(f2.abs()) as f32 as f64;
    if (f1 - f2).abs() > scal * eps {
        return Err(SolverError::fatal(CALLER, format!("{what} : {f1} {f2} {} > {}", (f1 - f2).abs(), scal * eps)));
    }
    Ok(())
}

/// The complete polynomial `sum x^i y^j z^k, i+j+k <= qp` of an `edim`
/// element with its gradient and second derivatives.
fn poly(edim: usize, qp: i32, x: f64, y: f64, z: f64) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let mut fx = 0.0;
    let mut d = [0.0f64; 3];
    let mut dd = [[0.0f64; 3]; 3];
    let pw = |b: f64, e: i32| if e < 0 { 0.0 } else { b.powi(e) };
    let (jmax, kmax) = match edim {
        1 => (0, 0),
        2 => (qp, 0),
        _ => (qp, qp),
    };
    for i in 0..=qp {
        for j in 0..=jmax {
            for k in 0..=kmax {
                if i + j + k > qp {
                    continue;
                }
                let (fi, fj, fk) = (i as f64, j as f64, k as f64);
                fx += pw(x, i) * pw(y, j) * pw(z, k);
                if i > 0 {
                    d[0] += fi * pw(x, i - 1) * pw(y, j) * pw(z, k);
                }
                if j > 0 {
                    d[1] += pw(x, i) * fj * pw(y, j - 1) * pw(z, k);
                }
                if k > 0 {
                    d[2] += pw(x, i) * pw(y, j) * fk * pw(z, k - 1);
                }
                if i > 1 {
                    dd[0][0] += fi * (fi - 1.0) * pw(x, i - 2) * pw(y, j) * pw(z, k);
                }
                if i > 0 && j > 0 {
                    dd[0][1] += fi * pw(x, i - 1) * fj * pw(y, j - 1) * pw(z, k);
                }
                if i > 0 && k > 0 {
                    dd[0][2] += fi * pw(x, i - 1) * pw(y, j) * fk * pw(z, k - 1);
                }
                if j > 1 {
                    dd[1][1] += pw(x, i) * fj * (fj - 1.0) * pw(y, j - 2) * pw(z, k);
                }
                if j > 0 && k > 0 {
                    dd[1][2] += pw(x, i) * fj * pw(y, j - 1) * fk * pw(z, k - 1);
                }
                if k > 1 {
                    dd[2][2] += pw(x, i) * pw(y, j) * fk * (fk - 1.0) * pw(z, k - 2);
                }
            }
        }
    }
    (fx, d, dd)
}

impl Solver for TestSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let layout = PDofsBackend::new()
            .for_solver(ctx.mesh, &ctx.sim.model, ctx.params(), &active, 1)
            .ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: element definition without p dofs")))??;
        self.layout = Some(layout);
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let d = self.layout.as_deref().ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: without p elements")))?;
        // `MAXVAL(Solver % Def_Dofs(:,1,6))`: the p of the definition.
        let qp: i32 = ctx
            .params()
            .string("element")
            .and_then(|e| e.trim().strip_prefix("p:").and_then(|p| p.trim().parse().ok()))
            .ok_or_else(|| SolverError::Unsupported(format!("{CALLER}: an Element other than p:<n>")))?;
        println!(" Testing 2nd derivatives with p({qp})...");
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut line = String::new();
        let mut slots = Vec::new();
        let mut ev = BasisEval::default();
        for e in active {
            let el = ctx.mesh.element(e);
            let code = el.etype.0;
            line.push_str(&format!("({code})..."));
            let def = elmer_fe::element_type(code).map_err(|m| SolverError::fatal(CALLER, m.to_string()))?;
            let edim = match code / 100 {
                2 => 1,
                3 | 4 => 2,
                _ => 3,
            };
            d.element_slots(e, &mut slots);
            let nd = slots.len();
            let nodes = d.element_nodes(e);

            // Nodes % x(1:n) = NodeU (in the element's node order), the rest zero.
            let mut xyz = vec![[0.0f64; 3]; nd];
            for (i, k) in nodes.iter().enumerate() {
                let pos = el.nodes.iter().position(|m| m == k).ok_or_else(|| SolverError::fatal(CALLER, "node order"))?;
                xyz[i] = [
                    def.node_u.get(pos).copied().unwrap_or(0.0),
                    if edim > 1 { def.node_v.get(pos).copied().unwrap_or(0.0) } else { 0.0 },
                    if edim > 2 { def.node_w.get(pos).copied().unwrap_or(0.0) } else { 0.0 },
                ];
            }
            let basis = d.element_basis_scalar(e).map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let rule = d.element_rule(e).map_err(|m| SolverError::fatal("GaussPoints", m))?;
            let mut stiff = vec![0.0f64; nd * nd];
            let mut force = vec![0.0f64; nd];
            let mut dnodal = vec![[0.0f64; 3]; nd * nd];
            for qpt in &rule {
                basis.eval(el.etype, &xyz, qpt.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let w = qpt.weight * ev.det_j;
                let at = |c: usize| (0..nd).fold(0.0, |a, i| a + xyz[i][c] * ev.basis[i]);
                let (x, y, z) = (at(0), at(1), at(2));
                let (fx, _, _) = poly(edim, qp, x, y, z);
                for p in 0..nd {
                    for q in 0..nd {
                        stiff[p * nd + q] += w * ev.basis[p] * ev.basis[q];
                        for c in 0..3 {
                            dnodal[p * nd + q][c] += w * ev.dbasis_dx[p][c] * ev.basis[q];
                        }
                    }
                    force[p] += w * ev.basis[p] * fx;
                }
            }
            // InvertMatrix (LinearAlgebra.F90): column-major, and the mass
            // matrix is symmetric.
            elmer_fe::invert_matrix(&mut stiff, nd).ok_or_else(|| SolverError::fatal("InvertMatrix", "singular matrix"))?;
            let f: Vec<f64> = (0..nd).map(|i| (0..nd).fold(0.0, |a, j| a + stiff[i * nd + j] * force[j])).collect();
            // dx(i, :, c) = STIFF^-1 dNodal(i, :, c).
            let mut dxp = vec![[0.0f64; 3]; nd * nd];
            for i in 0..nd {
                for c in 0..3 {
                    for r in 0..nd {
                        dxp[i * nd + r][c] = (0..nd).fold(0.0, |a, q| a + stiff[r * nd + q] * dnodal[i * nd + q][c]);
                    }
                }
            }
            for qpt in &rule {
                basis.eval(el.etype, &xyz, qpt.uvw, &mut ev).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let ddb = basis.second_derivatives(el.etype, &xyz, qpt.uvw).map_err(|m| SolverError::fatal("ElementInfo", m))?;
                let at = |c: usize| (0..nd).fold(0.0, |a, i| a + ev.basis[i] * xyz[i][c]);
                let (x, y, z) = (at(0), at(1), at(2));
                for i in 0..nd {
                    let mut from = [[0.0f64; 3]; 3];
                    for j in 0..3 {
                        for k in 0..3 {
                            from[j][k] = (0..nd).fold(0.0, |a, r| a + dxp[i * nd + r][j] * ev.dbasis_dx[r][k]);
                        }
                    }
                    let mx = (0..3).flat_map(|j| (0..3).map(move |k| (j, k))).map(|(j, k)| (from[j][k] - ddb[i][j][k]).abs()).fold(0.0, f64::max);
                    if mx > 1.0e-8 {
                        return Err(SolverError::fatal(
                            CALLER,
                            format!("ddx did not match (element {code}, function {} of {nd}): {:?} vs {:?}", i + 1, from, ddb[i]),
                        ));
                    }
                }
                let (fx, diff, ddiff) = poly(edim, qp, x, y, z);
                let sum = |g: &dyn Fn(usize) -> f64| (0..nd).fold(0.0, |a, i| a + g(i) * f[i]);
                check(fx, sum(&|i| ev.basis[i]), 1.0e-8, "fx")?;
                for (c, name) in ["dx", "dy", "dz"].iter().enumerate() {
                    check(diff[c], sum(&|i| ev.dbasis_dx[i][c]), 4.0e-7, name)?;
                }
                for (a, b, name, fname) in [(0, 0, "dxx", "dfxx"), (0, 1, "dxy", "dfxy"), (0, 2, "dxz", "dfxz"), (1, 1, "dyy", "dfyy"), (1, 2, "dyz", "dfyz"), (2, 2, "dzz", "dfzz")] {
                    check(ddiff[a][b], sum(&|i| ddb[i][a][b]), 5.0e-7, name)?;
                    let from = |i: usize| (0..nd).fold(0.0, |acc, r| acc + dxp[i * nd + r][a] * ev.dbasis_dx[r][b]);
                    check(ddiff[a][b], sum(&from), 7.5e-7, fname)?;
                }
            }
            line.push_str("PASSED ");
        }
        println!("{line}");
        println!();
        Ok(())
    }
}
