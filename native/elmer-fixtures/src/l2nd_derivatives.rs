//! `L2ndDerivatives/L2ndDerivativesTest.f90`: `SUBROUTINE TestSolver`, a
//! check of the Lagrange basis, its gradients and second derivatives on
//! every element of a mesh holding one element of each type.
//!
//! For each active element, on its reference geometry (the element's
//! `NodeU`, `NodeV`, `NodeW` as coordinates), the routine L2-projects the
//! polynomial `sum x^i y^j z^k` of total degree up to the element's
//! `BasisFunctionDegree` onto the basis, and the gradients of the basis
//! onto the basis as well. It then checks at every Gauss point that
//!
//! * the derivatives of the projected gradients equal `ddBasisddx` to 1e-8
//!   (skipped for pyramids), else `STOP "ddx's don't match"`;
//! * the projection reproduces the polynomial, its gradient and its second
//!   derivatives, through `dBasisdx`, `ddBasisddx` and the projected
//!   gradients, to 1e-12 (`CheckValue`).
//!
//! `CheckValue` compares against `scal*eps` with `scal` uninitialised; here a
//! difference of 1e-12 or more always stops, which is the reading under
//! which the test checks anything. A stop is a `Fatal`. The variable's
//! norm is the solver's business (`Solver % Variable % Norm = 1` is
//! overwritten by the framework's norm), so the routine sets nothing.

use elmer_core::{BasisEval, ElemId, Quadrature};
use elmer_fe::{element_info, element_info_second, element_type, ElmerQuadrature, ElementTypeDef, MetricStatus};
use elmer_solver::{Solver, SolverCtx, SolverError};

use crate::hooks::coordinate_system;
use crate::linalg::invert;
use crate::udf::{list_real, stdout_sink, PrintSink};

const CALLER: &str = "TestSolver";

/// gfortran's records: `PRINT *` items and non-advancing `WRITE(*,'(A)')`.
struct Printer {
    out: PrintSink,
    partial: String,
}

impl Printer {
    /// `PRINT *, text`: completes a pending non-advancing record.
    fn print(&mut self, text: &str) {
        let line = format!("{} {text}", std::mem::take(&mut self.partial));
        (self.out)(&line);
    }

    /// `WRITE(*,'(A)',ADVANCE='NO') text`.
    fn write_no_advance(&mut self, text: &str) {
        self.partial.push_str(text);
    }
}

pub struct TestSolver {
    printer: Printer,
}

/// [`SolverFactory`](elmer_solver::SolverFactory): prints to standard output.
pub fn factory(_params: &elmer_core::ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(TestSolver { printer: Printer { out: stdout_sink(), partial: String::new() } }))
}

/// `CheckValue(f1, f2, 1e-12, str)`.
fn check(f1: f64, f2: f64, what: &str) -> Result<(), SolverError> {
    let d = (f1 - f2).abs();
    if d < 1.0e-12 {
        return Ok(());
    }
    Err(SolverError::fatal(
        CALLER,
        format!("{what}:{}{}{} > 1e-12 (STOP)", list_real(f1), list_real(f2), list_real(d)),
    ))
}

/// The analytic polynomial of degree `qp` in `edim` variables at `(x, y,
/// z)`: its value, gradient and second derivatives, accumulated as the
/// Fortran loops do.
fn polynomial(edim: usize, qp: i32, x: f64, y: f64, z: f64) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let p = |b: f64, e: i32| b.powi(e);
    let (mut fx, mut diff, mut dd) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    match edim {
        1 => {
            for i in 0..=qp {
                fx += p(x, i);
                if i > 0 {
                    diff[0] += f64::from(i) * p(x, i - 1);
                }
                if i > 1 {
                    dd[0][0] += f64::from(i * (i - 1)) * p(x, i - 2);
                }
            }
        }
        2 => {
            for i in 0..=qp {
                for j in 0..=qp {
                    if i + j > qp {
                        continue;
                    }
                    let (fi, fj) = (f64::from(i), f64::from(j));
                    fx += p(x, i) * p(y, j);
                    if i > 0 {
                        diff[0] += fi * p(x, i - 1) * p(y, j);
                    }
                    if j > 0 {
                        diff[1] += p(x, i) * fj * p(y, j - 1);
                    }
                    if i > 1 {
                        dd[0][0] += fi * (fi - 1.0) * p(x, i - 2) * p(y, j);
                    }
                    if i > 0 && j > 0 {
                        dd[0][1] += fi * p(x, i - 1) * fj * p(y, j - 1);
                    }
                    if j > 1 {
                        dd[1][1] += p(x, i) * fj * (fj - 1.0) * p(y, j - 2);
                    }
                }
            }
        }
        _ => {
            for i in 0..=qp {
                for j in 0..=qp {
                    for k in 0..=qp {
                        if i + j + k > qp {
                            continue;
                        }
                        let (fi, fj, fk) = (f64::from(i), f64::from(j), f64::from(k));
                        fx += p(x, i) * p(y, j) * p(z, k);
                        if i > 0 {
                            diff[0] += fi * p(x, i - 1) * p(y, j) * p(z, k);
                        }
                        if j > 0 {
                            diff[1] += p(x, i) * fj * p(y, j - 1) * p(z, k);
                        }
                        if k > 0 {
                            diff[2] += p(x, i) * p(y, j) * fk * p(z, k - 1);
                        }
                        if i > 1 {
                            dd[0][0] += fi * (fi - 1.0) * p(x, i - 2) * p(y, j) * p(z, k);
                        }
                        if i > 0 && j > 0 {
                            dd[0][1] += fi * p(x, i - 1) * fj * p(y, j - 1) * p(z, k);
                        }
                        if i > 0 && k > 0 {
                            dd[0][2] += fi * p(x, i - 1) * p(y, j) * fk * p(z, k - 1);
                        }
                        if j > 1 {
                            dd[1][1] += p(x, i) * fj * (fj - 1.0) * p(y, j - 2) * p(z, k);
                        }
                        if j > 0 && k > 0 {
                            dd[1][2] += p(x, i) * fj * p(y, j - 1) * fk * p(z, k - 1);
                        }
                        if k > 1 {
                            dd[2][2] += p(x, i) * p(y, j) * fk * (fk - 1.0) * p(z, k - 2);
                        }
                    }
                }
            }
        }
    }
    (fx, diff, dd)
}

/// `LocalMatrix` of one element type on its reference geometry.
pub fn check_element(def: &ElementTypeDef, cdim: usize) -> Result<(), SolverError> {
    let fatal = |m: String| SolverError::fatal(CALLER, m);
    let nd = def.n_nodes;
    let edim = def.dim;
    let qp = def.basis_degree as i32;
    let xyz: Vec<[f64; 3]> = (0..nd)
        .map(|i| {
            [
                def.node_u[i],
                if edim > 1 { def.node_v[i] } else { 0.0 },
                if edim > 2 { def.node_w[i] } else { 0.0 },
            ]
        })
        .collect();
    let at = |b: &[f64], k: usize| xyz.iter().zip(b).fold(0.0, |s, (x, bp)| s + x[k] * bp);
    let rule = ElmerQuadrature.gauss_points(elmer_core::ElementType(def.code), None);

    // stiff(p,q) row-major; dnodal[k][p][q].
    let mut stiff = vec![0.0; nd * nd];
    let mut force = vec![0.0; nd];
    let mut dnodal = vec![vec![0.0; nd * nd]; 3];
    let mut ev = BasisEval::default();
    let mut status = MetricStatus::Unchecked;
    for q in &rule {
        element_info(def, &xyz, q.uvw, cdim, &mut status, &mut ev).map_err(|e| fatal(e.to_string()))?;
        let (b, g) = (&ev.basis, &ev.dbasis_dx);
        let (x, y, z) = (at(b, 0), at(b, 1), at(b, 2));
        let s = q.weight * ev.det_j;
        let (fx, _, _) = polynomial(edim, qp, x, y, z);
        for p in 0..nd {
            for qq in 0..nd {
                stiff[p * nd + qq] += s * b[p] * b[qq];
                for k in 0..3 {
                    dnodal[k][p * nd + qq] += s * g[p][k] * b[qq];
                }
            }
            force[p] += s * b[p] * fx;
        }
    }
    invert(&mut stiff, nd).map_err(fatal)?;
    // f = MATMUL(STIFF, FORCE); dx(i,:,j) = MATMUL(STIFF, dNodal(i,:,j)).
    let f: Vec<f64> = (0..nd).map(|r| (0..nd).fold(0.0, |s, c| s + stiff[r * nd + c] * force[c])).collect();
    let dx: Vec<[Vec<f64>; 3]> = (0..nd)
        .map(|i| {
            std::array::from_fn(|j| {
                (0..nd).map(|r| (0..nd).fold(0.0, |s, c| s + stiff[r * nd + c] * dnodal[j][i * nd + c])).collect()
            })
        })
        .collect();

    let mut dd = Vec::new();
    let mut status = MetricStatus::Unchecked;
    for q in &rule {
        element_info_second(def, &xyz, q.uvw, cdim, &mut status, &mut ev, &mut dd).map_err(|e| fatal(e.to_string()))?;
        let (b, g) = (&ev.basis, &ev.dbasis_dx);
        let (x, y, z) = (at(b, 0), at(b, 1), at(b, 2));
        // ddxFromNodaldx(i,j,k) = sum(dx(i,1:nd,j) * dBasisdx(1:nd,k))
        let ddn: Vec<[[f64; 3]; 3]> = (0..nd)
            .map(|i| std::array::from_fn(|j| std::array::from_fn(|k| (0..nd).fold(0.0, |s, l| s + dx[i][j][l] * g[l][k]))))
            .collect();
        if def.family() != 6 {
            for i in 0..nd {
                let worst = (0..3).flat_map(|j| (0..3).map(move |k| (j, k))).map(|(j, k)| (ddn[i][j][k] - dd[i][j][k]).abs());
                if worst.fold(0.0, f64::max) > 1.0e-8 {
                    return Err(fatal(format!("ddx's don't match (element {}, basis function {})", def.code, i + 1)));
                }
            }
        }
        let (fx, diff, ddiff) = polynomial(edim, qp, x, y, z);
        let sum = |v: &dyn Fn(usize) -> f64| (0..nd).fold(0.0, |s, l| s + v(l) * f[l]);
        check(fx, sum(&|l| b[l]), "fx")?;
        for (k, name) in ["dx", "dy", "dz"].into_iter().enumerate() {
            check(diff[k], sum(&|l| g[l][k]), name)?;
        }
        let pairs = [(0, 0, "xx"), (0, 1, "xy"), (0, 2, "xz"), (1, 1, "yy"), (1, 2, "yz"), (2, 2, "zz")];
        for (j, k, name) in pairs {
            check(ddiff[j][k], sum(&|l| dd[l][j][k]), &format!("d{name}"))?;
        }
        for (j, k, name) in pairs {
            check(ddiff[j][k], sum(&|l| ddn[l][j][k]), &format!("df{name}"))?;
        }
    }
    Ok(())
}

impl Solver for TestSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        // Elmer runs this on every element of the mesh; an empty table would
        // pass without checking anything.
        if ctx.active_elements().next().is_none() {
            return Err(SolverError::Unsupported("L2ndDerivativesTest TestSolver: no active elements".into()));
        }
        Ok(())
    }

    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let cdim = coordinate_system(ctx)?.dim;
        let active: Vec<ElemId> = ctx.active_elements().collect();
        self.printer.print("Testing basis functions and derivatives (1st and 2nd)");
        let mut prev_family = None;
        for e in active {
            let code = ctx.mesh.element(e).etype.0;
            let def = element_type(code).map_err(|err| SolverError::fatal(CALLER, err.to_string()))?;
            if prev_family != Some(def.family()) {
                self.printer.print("");
                self.printer.print("");
                prev_family = Some(def.family());
                let heading = match def.family() {
                    2 => Some("Line elements:"),
                    3 => Some("Triangular elements:"),
                    4 => Some("Quadrilateral elements:"),
                    5 => Some("Tetrahedral elements:"),
                    6 => Some("Pyramidal elements:"),
                    7 => Some("Wedge elements:"),
                    8 => Some("Hexahedral elements:"),
                    _ => None,
                };
                if let Some(h) = heading {
                    self.printer.print(h);
                }
            }
            self.printer.write_no_advance(&format!("({code})..."));
            check_element(def, cdim)?;
            self.printer.write_no_advance("PASSED ");
        }
        self.printer.print("");
        self.printer.print("");
        Ok(())
    }
}
