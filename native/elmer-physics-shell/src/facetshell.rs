//! `ShellSolver` of FacetShellSolve.F90: the Reissner-Mindlin facet shell.
//! Each flat 303 or 404 element gets its own frame (`LocalBasis`); the
//! plate bending (MITC shear), membrane and drilling stiffness are built
//! in it and rotated to the global six dofs (displacements and rotations).
//!
//! Ported: the linear static path, with `Body Force 1..3`, `Pressure` and
//! `Normal Pressure` and the `Shear` / `Drilling Stabilization Parameter`.
//! Refused: `Large Deflection`, `Stability Analysis`, stress computation,
//! eigen and harmonic analysis, transient runs, `Load Scale Factor`, the
//! `Force 1..6` boundary loads (`StressBoundary`), `Point Load`s and more
//! than one partition.

use elmer_core::{ElemId, LocalSystem, ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::smitc::{jacobi3, jacobi4, shear_correction, strain_energy, sum, tj};

pub const FACET_PROCEDURE: (&str, &str) = ("FacetShellSolve", "ShellSolver");

const CALLER: &str = "ShellSolver";
const DOFS: usize = 6;

pub fn register_facet(registry: &mut SolverRegistry) {
    registry.register(FACET_PROCEDURE.0, FACET_PROCEDURE.1, facet_factory);
}

pub fn facet_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(FacetShellSolver))
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("FacetShellSolve: {what}"))
}

/// `ShellSolver_Init0` of FacetShellSolve: `Variable = Deflection` with six
/// dofs when none is given, `Time Derivative Order = 2`, `Shell Solver` and
/// `Drilling DOFs`.
pub fn facet_defaults(p: &mut ParamMap) {
    if p.get("variable").is_none() {
        p.insert("Variable DOFs", Value::Integer(vec![6]));
        p.insert("Variable", Value::Str("Deflection".into()));
    }
    p.insert("Time derivative order", Value::Integer(vec![2]));
    p.insert("Shell Solver", Value::Logical(true));
    p.insert("Drilling DOFs", Value::Logical(true));
}

pub struct FacetShellSolver;

type V3 = [f64; 3];

/// `SQRT( SUM( a**2 ) )`.
fn norm(a: &V3) -> f64 {
    let mut s = 0.0;
    for x in a {
        s += x.powi(2);
    }
    s.sqrt()
}

/// `SUM( a * b )`.
fn dot(a: &V3, b: &V3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        s += a[i] * b[i];
    }
    s
}

/// `TRANSPOSE( LocalBasis( Nodes, n ) )`: the rows are the element's
/// tangent, in-plane normal and normal.
fn transformation(x: &[V3]) -> [V3; 3] {
    let t = |a: usize, b: usize| -> V3 {
        let mut v = [x[b][0] - x[a][0], x[b][1] - x[a][1], x[b][2] - x[a][2]];
        let n = norm(&v);
        for c in &mut v {
            *c /= n;
        }
        v
    };
    let t1 = t(0, 1);
    let t2 = t(1, 2);
    let d = dot(&t1, &t2);
    let mut b2 = [t2[0] - d * t1[0], t2[1] - d * t1[1], t2[2] - d * t1[2]];
    let n = norm(&b2);
    for c in &mut b2 {
        *c /= n;
    }
    // CrossProduct( b1, b2 ).
    let b3 = [
        t1[1] * b2[2] - t1[2] * b2[1],
        -t1[0] * b2[2] + t1[2] * b2[0],
        t1[0] * b2[1] - t1[1] * b2[0],
    ];
    [t1, b2, b3]
}

/// `SwitchToLocal`: `T (X - centroid)`.
fn switch_to_local(x: &[V3], t: &[V3; 3]) -> Vec<V3> {
    let n = x.len();
    let mut c = [0.0f64; 3];
    for k in 0..3 {
        let mut s = 0.0;
        for p in x {
            s += p[k];
        }
        c[k] = s / n as f64;
    }
    x.iter()
        .map(|p| {
            let g = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            std::array::from_fn(|i| {
                let mut s = 0.0;
                for k in 0..3 {
                    s += t[i][k] * g[k];
                }
                s
            })
        })
        .collect()
}

/// `CovariantInterpolation` of FacetShellSolve: the MITC shear strain
/// `(2, 6n)` on the rotation dofs 4 and 5 of each node.
fn covariant_interpolation(
    x: &[f64],
    y: &[f64],
    u: f64,
    v: f64,
) -> Result<Vec<Vec<f64>>, SolverError> {
    let n = x.len();
    let nc = DOFS * n;
    let mut sr = vec![vec![0.0f64; nc]; 2];
    let inv = match n {
        3 => {
            let (jm, inv) = jacobi3(x, y);
            let r2 = 2.0f64.sqrt();
            // (Tau, Sdofs columns (1-based), divisor, factors for rows 1, 2)
            let edges: [([f64; 2], [usize; 4], f64, f64, f64); 3] = [
                ([1.0, 0.0], [4, 5, 10, 11], 2.0, 1.0 - v, u),
                ([-1.0 / r2, 1.0 / r2], [10, 11, 16, 17], r2, -v, u),
                ([0.0, -1.0], [4, 5, 16, 17], 2.0, -v, -1.0 + u),
            ];
            for (tau, cols, d, f1, f2) in edges {
                let mut sd = vec![0.0f64; nc];
                sd[cols[0] - 1] = tj(tau, &jm, 0, d);
                sd[cols[1] - 1] = tj(tau, &jm, 1, d);
                sd[cols[2] - 1] = tj(tau, &jm, 0, d);
                sd[cols[3] - 1] = tj(tau, &jm, 1, d);
                for j in 0..18 {
                    sr[0][j] += f1 * sd[j];
                    sr[1][j] += f2 * sd[j];
                }
            }
            inv
        }
        4 => {
            let edges: [(f64, f64, [f64; 2], [usize; 4], usize); 4] = [
                (0.0, -1.0, [1.0, 0.0], [4, 5, 10, 11], 0),
                (1.0, 0.0, [0.0, 1.0], [10, 11, 16, 17], 1),
                (0.0, 1.0, [-1.0, 0.0], [16, 17, 22, 23], 0),
                (-1.0, 0.0, [0.0, -1.0], [4, 5, 22, 23], 1),
            ];
            for (k, (xi, eta, tau, cols, row)) in edges.into_iter().enumerate() {
                let (jm, _) = jacobi4(xi, eta, x, y);
                let mut sd = vec![0.0f64; nc];
                sd[cols[0] - 1] = tj(tau, &jm, 0, 1.0);
                sd[cols[1] - 1] = tj(tau, &jm, 1, 1.0);
                sd[cols[2] - 1] = tj(tau, &jm, 0, 1.0);
                sd[cols[3] - 1] = tj(tau, &jm, 1, 1.0);
                let f = match k {
                    0 => (1.0 - v) / 4.0,
                    1 => (1.0 + u) / 4.0,
                    2 => (-1.0 - v) / 4.0,
                    _ => (-1.0 + u) / 4.0,
                };
                for j in 0..24 {
                    sr[row][j] += f * sd[j];
                }
            }
            jacobi4(u, v, x, y).1
        }
        _ => {
            return Err(SolverError::fatal(
                CALLER,
                "Illegal number of nodes for Smitc elements.",
            ));
        }
    };
    Ok((0..2)
        .map(|i| {
            (0..nc)
                .map(|j| -(0.0 + inv[i][0] * sr[0][j] + inv[i][1] * sr[1][j]))
                .collect()
        })
        .collect())
}

/// `MATMUL( A, B )` of `m x m` row-major matrices, ascending sums.
fn matmul(a: &[f64], b: &[f64], m: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut s = 0.0;
            for k in 0..m {
                s += a[i * m + k] * b[k * m + j];
            }
            c[i * m + j] = s;
        }
    }
    c
}

/// `MATMUL( TRANSPOSE( A ), B )`.
fn matmul_tn(a: &[f64], b: &[f64], m: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut s = 0.0;
            for k in 0..m {
                s += a[k * m + i] * b[k * m + j];
            }
            c[i * m + j] = s;
        }
    }
    c
}

/// `K = Tb^T (K Tb)`.
fn rotate(k: &[f64], tb: &[f64], m: usize) -> Vec<f64> {
    let kt = matmul(k, tb, m);
    matmul_tn(tb, &kt, m)
}

struct Material {
    poisson: Vec<f64>,
    young: Vec<f64>,
    thickness: Vec<f64>,
    load: [Vec<f64>; 4],
}

impl Solver for FacetShellSolver {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in [
            "large deflection",
            "stability analysis",
            "stress computation",
            "calculate stresses",
            "eigen analysis",
            "harmonic mode",
            "harmonic analysis",
            "linear system solver disabled",
        ] {
            if p.logical(k) == Some(true) {
                return Err(unsupported(k));
            }
        }
        if p.get("load scale factor").is_some() {
            return Err(unsupported("Load Scale Factor"));
        }
        let force = (1..=6).map(|i| format!("force {i}")).collect::<Vec<_>>();
        if ctx
            .sim
            .model
            .boundary_conditions
            .iter()
            .any(|b| force.iter().any(|k| b.get(k).is_some()))
        {
            return Err(unsupported("Force 1..6 boundary loads (StressBoundary)"));
        }
        if ctx
            .sim
            .model
            .body_forces
            .iter()
            .any(|b| b.get("point load").is_some())
        {
            return Err(unsupported("Point Load"));
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("more than one partition"));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        transient: bool,
    ) -> Result<(), SolverError> {
        if transient {
            return Err(unsupported("transient shells (Default2ndOrderTime)"));
        }
        ctx.info(CALLER, "Reissner-Mindlin facet shell solver", 4);
        let p = ctx.params().clone();
        let stab1 = match ctx.creal(&p, "shear stabilization parameter")? {
            Some(v) => v,
            None => {
                ctx.info(CALLER, "Shear stabilization parameter undefined.", 6);
                ctx.info(CALLER, "Using default value 1.0", 6);
                1.0
            }
        };
        let stab2 = match ctx.creal(&p, "drilling stabilization parameter")? {
            Some(v) => v,
            None => {
                ctx.info(CALLER, "Drilling stabilization parameter undefined", 6);
                ctx.info(CALLER, "Using default value 1.0", 6);
                1.0
            }
        };
        let elems: Vec<ElemId> = ctx.active_elements().collect();
        ctx.default_initialize()?;
        for &e in &elems {
            let et = ctx.mesh.element(e).etype.0;
            if et != 303 && et != 404 {
                ctx.warn(CALLER, "Implemented only for elements 303 and 404!");
                continue;
            }
            let m = material(ctx, e)?;
            let ls = local_matrix(ctx, e, &m, stab1, stab2)?;
            ctx.default_update(&ls, e)?;
        }
        ctx.default_finish_bulk_assembly()?;
        ctx.default_finish_boundary_assembly()?;
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}

fn material(ctx: &SolverCtx<'_>, e: ElemId) -> Result<Material, SolverError> {
    let n = ctx.mesh.element(e).nodes.len();
    let bp = ctx.body_params(e);
    let bf = bp.as_ref().and_then(|b| b.body_force);
    let mat = bp.as_ref().and_then(|b| b.material);
    let get = |s: Option<&ParamMap>, k: &str| -> Result<Option<Vec<f64>>, SolverError> {
        match s {
            Some(s) => ctx.real_at_nodes(s, k, e),
            None => Ok(None),
        }
    };
    let zero = || vec![0.0; n];
    let need = |k: &str, what: &str| -> Result<Vec<f64>, SolverError> {
        get(mat, k)?.ok_or_else(|| SolverError::fatal(CALLER, format!("{what} undefined")))
    };
    let poisson = need("poisson ratio", "Poisson ratio")?;
    let young = need("youngs modulus", "Youngs modulus")?;
    let thickness = match get(mat, "thickness")? {
        Some(t) => t,
        None => get(mat, "shell thickness")?
            .ok_or_else(|| SolverError::fatal(CALLER, "Shell Thickness/Thickness undefined"))?,
    };
    let mut ln = get(bf, "pressure")?.unwrap_or_else(zero);
    if let Some(np) = get(bf, "normal pressure")? {
        for (a, b) in ln.iter_mut().zip(&np) {
            *a += b;
        }
    }
    let load = [
        get(bf, "body force 1")?.unwrap_or_else(zero),
        get(bf, "body force 2")?.unwrap_or_else(zero),
        get(bf, "body force 3")?.unwrap_or_else(zero),
        ln,
    ];
    Ok(Material {
        poisson,
        young,
        thickness,
        load,
    })
}

/// `LocalMatrix` of FacetShellSolve on the linear static path: the terms
/// of `Bmatrix = 0`, the zero `ZetaStrain` and the zero membrane forces add
/// nothing and are left out.
fn local_matrix(
    ctx: &SolverCtx<'_>,
    e: ElemId,
    m: &Material,
    stab1: f64,
    stab2: f64,
) -> Result<LocalSystem, SolverError> {
    let el = ctx.mesh.element(e);
    let n = el.nodes.len();
    let nd = DOFS * n;
    let global: Vec<V3> = el.nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
    let t = transformation(&global);
    let local = switch_to_local(&global, &t);
    let x: Vec<f64> = local.iter().map(|c| c[0]).collect();
    let y: Vec<f64> = local.iter().map(|c| c[1]).collect();
    let fe = ctx
        .sim
        .backends
        .fe
        .as_ref()
        .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
    let np = if n == 3 { 1 } else { 4 };
    let mut stiff = vec![0.0f64; nd * nd];
    let mut force = vec![0.0f64; nd];
    let mut ev = elmer_core::BasisEval::default();
    for qp in fe.quadrature.gauss_points(el.etype, Some(np)) {
        fe.basis
            .eval(el.etype, &local, qp.uvw, &mut ev)
            .map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
        let b = &ev.basis[..n];
        let db = &ev.dbasis_dx[..n];
        let s = qp.weight * ev.det_j;
        let lx = sum(&m.load[0], b);
        let ly = sum(&m.load[1], b);
        let lz = sum(&m.load[2], b);
        let ln = sum(&m.load[3], b);
        let h = sum(&m.thickness, b);
        // IsotropicElasticity: Dmatrix and Astarmatrix.
        let eu = sum(&m.young, b);
        let pu = sum(&m.poisson, b);
        let tu = sum(&m.thickness, b);
        let gu = eu / (2.0 * (1.0 + pu));
        let base = |c33: f64| {
            let mut e = vec![vec![0.0f64; 3]; 3];
            e[0][0] = 1.0;
            e[0][1] = pu;
            e[1][0] = pu;
            e[1][1] = 1.0;
            e[2][2] = c33;
            e
        };
        let mut dm = base((1.0 - pu) / 2.0);
        for r in dm.iter_mut() {
            for v in r.iter_mut() {
                *v = *v * eu * tu.powi(3) / (12.0 * (1.0 - pu.powi(2)));
            }
        }
        let astar = vec![vec![gu * tu, 0.0], vec![0.0, gu * tu]];
        // IsotropicInPlaneElasticity.
        let mut am = base((1.0 - pu) / 2.0);
        for r in am.iter_mut() {
            for v in r.iter_mut() {
                *v = *v * tu * eu / (1.0 - pu * pu);
            }
        }
        let gdrill = vec![vec![stab2 * (astar[0][0] + astar[1][1])]];
        let mut kappa = vec![vec![0.0f64; nd]; 3];
        for p in 0..n {
            kappa[0][6 * p + 3] = db[p][0];
            kappa[1][6 * p + 4] = db[p][1];
            kappa[2][6 * p + 3] = db[p][1];
            kappa[2][6 * p + 4] = db[p][0];
        }
        strain_energy(&mut stiff, &dm, &kappa, 3, nd, s);
        let mut eps = vec![vec![0.0f64; nd]; 3];
        for p in 0..n {
            eps[0][6 * p] = db[p][0];
            eps[1][6 * p + 1] = db[p][1];
            eps[2][6 * p] = db[p][1];
            eps[2][6 * p + 1] = db[p][0];
        }
        strain_energy(&mut stiff, &am, &eps, 3, nd, s);
        let mut gamma = covariant_interpolation(&x, &y, qp.uvw[0], qp.uvw[1])?;
        let scf = shear_correction(h, &x, &y, stab1)?;
        for p in 0..n {
            gamma[0][6 * p + 2] = db[p][0];
            gamma[1][6 * p + 2] = db[p][1];
        }
        strain_energy(&mut stiff, &astar, &gamma, 2, nd, scf * s);
        let mut omega = vec![vec![0.0f64; nd]];
        for p in 0..n {
            omega[0][6 * p] = db[p][1] / 2.0;
            omega[0][6 * p + 1] = -db[p][0] / 2.0;
            omega[0][6 * p + 5] = b[p];
        }
        strain_energy(&mut stiff, &gdrill, &omega, 1, nd, s);
        for p in 0..n {
            force[6 * p] += lx * b[p] * h * s;
            force[6 * p + 1] += ly * b[p] * h * s;
            force[6 * p + 2] += lz * b[p] * h * s;
            for k in 0..3 {
                force[6 * p + k] += t[2][k] * ln * b[p] * s;
            }
        }
    }
    // The rotation dofs (theta_x, theta_y) -> (-theta_y, theta_x) ...
    let mut tb = vec![0.0f64; nd * nd];
    for i in 0..n {
        let o = 6 * i;
        tb[o * nd + o] = 1.0;
        tb[(o + 1) * nd + o + 1] = 1.0;
        tb[(o + 2) * nd + o + 2] = 1.0;
        tb[(o + 3) * nd + o + 4] = -1.0;
        tb[(o + 4) * nd + o + 3] = 1.0;
        tb[(o + 5) * nd + o + 5] = 1.0;
    }
    let k1 = rotate(&stiff, &tb, nd);
    // ... and the element frame to the global one.
    let mut tb = vec![0.0f64; nd * nd];
    for i in 0..n {
        for p in 0..3 {
            for q in 0..3 {
                tb[(6 * i + p) * nd + 6 * i + q] = t[p][q];
                tb[(6 * i + 3 + p) * nd + 6 * i + 3 + q] = t[p][q];
            }
        }
    }
    let k2 = rotate(&k1, &tb, nd);
    // STIFF = ( STIFF + TRANSPOSE(STIFF) ) / 2.
    let mut k = vec![0.0f64; nd * nd];
    for i in 0..nd {
        for j in 0..nd {
            k[i * nd + j] = (k2[i * nd + j] + k2[j * nd + i]) / 2.0;
        }
    }
    Ok(LocalSystem {
        n: nd,
        stiff: k,
        mass: None,
        damp: None,
        force,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A square in the plane z = 1 gets the global axes as its frame and
    /// centred local coordinates in its plane.
    #[test]
    fn frame_of_a_square_in_a_coordinate_plane() {
        let x = [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [2.0, 2.0, 1.0], [0.0, 2.0, 1.0]];
        let t = transformation(&x);
        assert_eq!(t, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let l = switch_to_local(&x, &t);
        assert_eq!(l, vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]);
    }

    /// Tb^T K Tb with a permutation-with-sign Tb permutes K exactly.
    #[test]
    fn rotation_by_a_signed_permutation_is_exact() {
        let k = [1.0, 2.0, 3.0, 4.0];
        let tb = [0.0, -1.0, 1.0, 0.0];
        assert_eq!(rotate(&k, &tb, 2), vec![4.0, -3.0, -2.0, 1.0]);
    }
}
