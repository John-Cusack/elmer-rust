//! The geometry of a closed coil: `DefineCoilCenter` and
//! `DefineCoilParameters` (the inertia fit of the coil axis by `DSYEV`),
//! with `TangentDirections` and `CrossProduct` of Elmer's utilities.

use elmer_core::{ElemId, ParamMap, Params, Value};
use elmer_solver::{SolverCtx, SolverError};

use crate::solver::element_points;

const CALLER: &str = "CoilSolver";

/// The two directions of the coil's plane (the axis itself is only used
/// by `Coil Helicity`, which is refused).
pub(crate) struct Axes {
    pub tangent1: [f64; 3],
    pub tangent2: [f64; 3],
}

fn array3(p: &ParamMap, key: &str) -> Option<[f64; 3]> {
    match p.get(key)? {
        Value::Real { data, .. } if data.len() >= 3 => Some([data[0], data[1], data[2]]),
        _ => None,
    }
}

fn column3(v: [f64; 3]) -> Value {
    Value::Real {
        rows: 3,
        cols: 1,
        data: v.to_vec(),
    }
}

/// `CrossProduct( v1, v2 )`.
pub(crate) fn cross(v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    [
        v1[1] * v2[2] - v1[2] * v2[1],
        -v1[0] * v2[2] + v1[2] * v2[0],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
}

/// `SQRT( SUM( v**2 ) )`.
fn norm3(v: [f64; 3]) -> f64 {
    let mut s = 0.0;
    for x in v {
        s += x * x;
    }
    s.sqrt()
}

/// `TangentDirections( Normal, Tangent1, Tangent2 )`.
pub(crate) fn tangent_directions(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (n1, n2, n3) = (n[0].abs(), n[1].abs(), n[2].abs());
    let mut t1 = if n1 <= n3 && n2 <= n3 {
        [0.0, -n[2], n[1]]
    } else {
        [-n[1], n[0], 0.0]
    };
    let s = norm3(t1);
    for x in &mut t1 {
        *x /= s;
    }
    let mut t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];
    let s = norm3(t2);
    for x in &mut t2 {
        *x /= s;
    }
    (t1, t2)
}

/// `r = SUM( Nodes % x(1:n) * Basis(1:n) )`, ... at a Gauss point.
fn point(xyz: &[[f64; 3]], basis: &[f64]) -> [f64; 3] {
    let mut r = [0.0f64; 3];
    for k in 0..3 {
        let mut s = 0.0;
        for (x, b) in xyz.iter().zip(basis) {
            s += x[k] * b;
        }
        r[k] = s;
    }
    r
}

/// `DefineCoilCenter( CoilCenter, Params )`: `Coil Center` (or `Coil x0`,
/// `y0`, `z0`), else the volume centre of the active elements, which is then
/// stored as `Coil Center`.
pub(crate) fn coil_center(
    ctx: &mut SolverCtx<'_>,
    p: &ParamMap,
    elems: &[ElemId],
) -> Result<[f64; 3], SolverError> {
    let mut center = [0.0f64; 3];
    if let Some(c) = array3(p, "coil center") {
        ctx.info(CALLER, "Coil center defined by user", 20);
        return Ok(c);
    }
    let mut set = false;
    for (k, key) in ["coil x0", "coil y0", "coil z0"].into_iter().enumerate() {
        if let Some(v) = ctx.creal(p, key)? {
            center[k] = v;
            set = true;
        }
    }
    if set {
        ctx.info(CALLER, "Coil center defined by user", 20);
        return Ok(center);
    }
    let mut volume = 0.0;
    for &e in elems {
        let xyz: Vec<[f64; 3]> = ctx
            .mesh
            .element(e)
            .nodes
            .iter()
            .map(|&k| ctx.mesh.coord(k))
            .collect();
        for pt in element_points(ctx, e)? {
            let r = point(&xyz, &pt.basis);
            let s = pt.weight;
            volume += s;
            for k in 0..3 {
                center[k] += s * r[k];
            }
        }
    }
    let mut tmp = [center[0], center[1], center[2], volume];
    ctx.comm.sum_f64(&mut tmp);
    center = [tmp[0], tmp[1], tmp[2]];
    volume = tmp[3];
    if volume < f64::EPSILON {
        return Err(SolverError::fatal(
            "DefineCoilCenter",
            "Coil has no volume!",
        ));
    }
    for c in &mut center {
        *c /= volume;
    }
    ctx.info(CALLER, &format!("Coil volume: {volume:e}"), 7);
    ctx.info(CALLER, &format!("Coil center: {center:?}"), 7);
    ctx.set_param("Coil Center", column3(center));
    Ok(center)
}

/// `DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2, Params )`:
/// a given `Coil Normal` (and `Coil Tangent`), or with `Fit Coil` (default)
/// the axis of maximum inertia (minimum for a tall coil), signed by the
/// largest component of `center x normal`; the fitted normal is stored as
/// `Coil Normal`.
pub(crate) fn coil_axes(
    ctx: &mut SolverCtx<'_>,
    p: &ParamMap,
    elems: &[ElemId],
    center: [f64; 3],
) -> Result<Axes, SolverError> {
    let fit = p.logical("fit coil").unwrap_or(true);
    let given = array3(p, "coil normal");
    if !fit || given.is_some() {
        let normal = given.unwrap_or([0.0, 0.0, 1.0]);
        if given.is_some() {
            if let Some(t1) = array3(p, "coil tangent") {
                let mut t2 = cross(normal, t1);
                let s = norm3(t2);
                for x in &mut t2 {
                    *x /= s;
                }
                return Ok(Axes {
                    tangent1: t1,
                    tangent2: t2,
                });
            }
        }
        let (tangent1, tangent2) = tangent_directions(normal);
        return Ok(Axes { tangent1, tangent2 });
    }
    ctx.info(
        "DefineCoilParameters",
        "Fitting the coil by maximizing inertia",
        7,
    );
    let mut im = [0.0f64; 9];
    for &e in elems {
        let xyz: Vec<[f64; 3]> = ctx
            .mesh
            .element(e)
            .nodes
            .iter()
            .map(|&k| ctx.mesh.coord(k))
            .collect();
        for pt in element_points(ctx, e)? {
            let mut r = point(&xyz, &pt.basis);
            let s = pt.weight;
            for k in 0..3 {
                r[k] -= center[k];
            }
            let mut r2 = 0.0;
            for x in r {
                r2 += x * x;
            }
            for i in 0..3 {
                im[3 * i + i] += s * r2;
                for j in 0..3 {
                    im[3 * i + j] -= s * r[i] * r[j];
                }
            }
        }
    }
    ctx.comm.sum_f64(&mut im);
    let mut s = 0.0;
    for x in im {
        s += x;
    }
    for x in &mut im {
        *x /= s;
    }
    // EigVec(i,j) = Imoment(3*(i-1)+j), column-major; then EigVec(i,i) - 1.
    let mut a = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            a[i + 3 * j] = im[3 * i + j];
        }
        a[i + 3 * i] -= 1.0;
    }
    let w = dsyev_vu(&mut a, 3)?;
    ctx.info(CALLER, &format!("Coil inertia eigenvalues: {w:?}"), 10);
    let tall = match p.logical("coil geometry tall") {
        Some(t) => t,
        None => (w[1] - w[2]).abs() < 0.01 * (w[1] - w[0]).abs(),
    };
    let mut normal = if tall {
        ctx.info(
            CALLER,
            "Assuming coil to be tall i.e. axis has minimum inertia",
            4,
        );
        [a[0], a[1], a[2]]
    } else {
        ctx.info(
            CALLER,
            "Assuming coil to be flat i.e. axis has maximum inertia",
            4,
        );
        [a[6], a[7], a[8]]
    };
    let cp = cross(center, normal);
    let mut j = 0;
    for i in 1..3 {
        if cp[i].abs() > cp[j].abs() {
            j = i;
        }
    }
    if cp[j] < 0.0 {
        ctx.info(CALLER, "Inverting sign of coil normal", 10);
        for x in &mut normal {
            *x = -*x;
        }
    }
    let (tangent1, tangent2) = tangent_directions(normal);
    ctx.info(CALLER, &format!("Coil axis normal: {normal:?}"), 7);
    ctx.set_param("Coil Normal", column3(normal));
    Ok(Axes { tangent1, tangent2 })
}

/// `DSYEV( 'V', 'U', n, A, n, W, ... )` ([`elmer_eigen::dsyev`]): `a`
/// becomes the eigenvectors; the eigenvalues ascending.
pub(crate) fn dsyev_vu(a: &mut [f64], n: usize) -> Result<Vec<f64>, SolverError> {
    elmer_eigen::dsyev(elmer_eigen::Jobz::Vectors, a, n).map_err(|e| match e {
        elmer_eigen::DsyevError::Scaling => SolverError::Unsupported("DSYEV scaling of the coil inertia".into()),
        elmer_eigen::DsyevError::NoConvergence => SolverError::fatal(CALLER, "DSYEV cannot generate eigen basis"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangents_are_orthonormal_to_the_normal() {
        let n = [0.0, 0.6, 0.8];
        let (t1, t2) = tangent_directions(n);
        assert_eq!(t1, [0.0, -0.8, 0.6]);
        let d = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        assert!(d(t1, n).abs() < 1e-15 && d(t2, n).abs() < 1e-15 && d(t1, t2).abs() < 1e-15);
    }
}
