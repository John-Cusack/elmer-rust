//! ElasticSolve's UMAT materials (`LocalMatrixWithUMAT`, ElasticSolve.F90:2989,
//! and its incremental nonlinear iteration :1631-1810): the material routine
//! returns the stress and its tangent `DDSDDE` at each integration point from
//! the strain of the last converged state and its increment, and the
//! integration-point state (`UmatStress`, `UmatEnergy`, `UmatState`) is
//! committed when the iteration converges.
//!
//! Routines: `UMATLib`'s `linear_isotropic` (fem/src/modules/UMATLib.F90:149),
//! and whatever a case registers with [`register_umat`] (its own `UMAT
//! Subroutine` library, ported as a fixture). Small strains only here: `Large
//! Deflection`, `Use Hencky Strain`, axial symmetry, plane stress and
//! transient runs are refused by the solver.

use std::sync::{Mutex, OnceLock};

use elmer_core::{BasisEval, ElementBasis, ElementType, QuadPoint};

/// One call of a UMAT (the Abaqus interface `LocalMatrixWithUMAT` passes).
pub struct UmatCall<'a> {
    /// `STRESS(NTENS)`, in: the converged stress, out: the new one.
    pub stress: &'a mut [f64],
    /// `STATEV(NSTATEV)`.
    pub statev: &'a mut [f64],
    /// `DDSDDE(NTENS,NTENS)`, row-major, `ntens x ntens` used.
    pub ddsdde: &'a mut [[f64; 6]; 6],
    /// `SSE`, `SPD`, `SCD`.
    pub energy: &'a mut [f64; 3],
    /// `STRAN(NTENS)` (engineering shears) and `DSTRAN(NTENS)`.
    pub stran: &'a [f64],
    pub dstran: &'a [f64],
    /// `TIME(2)`, `DTIME`, `TEMP`.
    pub time: [f64; 2],
    pub dtime: f64,
    pub temp: f64,
    /// `CMNAME`: the material's `Name`.
    pub cmname: &'a str,
    pub ndi: usize,
    pub nshr: usize,
    pub ntens: usize,
    /// `PROPS(NPROPS)`: `Material Constants`.
    pub props: &'a [f64],
    /// `DFRGRD0`, `DFRGRD1`.
    pub dfgrd0: [[f64; 3]; 3],
    pub dfgrd1: [[f64; 3]; 3],
    /// `NOEL` (the element's `t`, 1-based) and `NPT` (the point, 1-based).
    pub noel: usize,
    pub npt: usize,
}

/// A UMAT routine.
pub type UmatFn = fn(&mut UmatCall<'_>) -> Result<(), String>;

fn registry() -> &'static Mutex<Vec<(String, String, UmatFn)>> {
    static R: OnceLock<Mutex<Vec<(String, String, UmatFn)>>> = OnceLock::new();
    R.get_or_init(|| {
        Mutex::new(vec![
            ("umatlib".into(), "linear_isotropic".into(), linear_isotropic as UmatFn),
            ("umatlib".into(), "stvenant_kirchhoff".into(), stvenant_kirchhoff as UmatFn),
        ])
    })
}

/// Makes `UMAT Subroutine = File "library" "entry"` resolve to `f`.
pub fn register_umat(library: &str, entry: &str, f: UmatFn) {
    let mut r = registry().lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let (l, e) = (library.to_ascii_lowercase(), entry.to_ascii_lowercase());
    r.retain(|(a, b, _)| !(a == &l && b == &e));
    r.push((l, e, f));
}

/// `GetProcAddr(UMATName)`.
pub(crate) fn lookup(library: &str, entry: &str) -> Option<UmatFn> {
    let r = registry().lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let (l, e) = (library.to_ascii_lowercase(), entry.to_ascii_lowercase());
    r.iter().find(|(a, b, _)| a == &l && b == &e).map(|t| t.2)
}

/// `UMATLib.F90` `linear_isotropic`: `PROPS = (*, E, nu)`, `DDSDDE` the
/// isotropic stiffness on `ndi` direct and `nshr` engineering shear
/// components, `STRESS = STRESS + DDSDDE DSTRAN`.
pub fn linear_isotropic(c: &mut UmatCall<'_>) -> Result<(), String> {
    if c.props.len() < 3 {
        return Err("linear_isotropic: 3 material constants are needed".into());
    }
    let e = c.props[1];
    let nu = c.props[2];
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    *c.ddsdde = [[0.0; 6]; 6];
    for i in 0..c.ndi {
        for j in 0..c.ndi {
            c.ddsdde[i][j] = lambda;
        }
    }
    for i in 0..c.ntens {
        c.ddsdde[i][i] += mu;
    }
    for i in 0..c.ndi {
        c.ddsdde[i][i] += mu;
    }
    let mut ds = [0.0; 6];
    for (i, d) in ds.iter_mut().enumerate().take(c.ntens) {
        for j in 0..c.ntens {
            *d += c.ddsdde[i][j] * c.dstran[j];
        }
    }
    for i in 0..c.ntens {
        c.stress[i] += ds[i];
    }
    Ok(())
}

/// The Voigt basis tensors `SymBasis(1:6)` of UMATLib.F90.
fn sym_basis(i: usize) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    match i {
        0 => m[0][0] = 1.0,
        1 => m[1][1] = 1.0,
        2 => m[2][2] = 1.0,
        3 => {
            m[0][1] = 0.5;
            m[1][0] = 0.5;
        }
        4 => {
            m[0][2] = 0.5;
            m[2][0] = 0.5;
        }
        _ => {
            m[1][2] = 0.5;
            m[2][1] = 0.5;
        }
    }
    m
}

/// `MATMUL(F, MATMUL(A, TRANSPOSE(F)))`.
fn push(f: &[[f64; 3]; 3], a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    use elmer_physics_elastic::post::{matmul3, transpose3};
    matmul3(f, &matmul3(a, &transpose3(f)))
}

/// The shear pairs `(1,2) (1,3) (2,3)` of the Voigt components after `ndi`.
const SHEAR: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

/// `UMATLib.F90` `stvenant_kirchhoff` (:312): the St. Venant-Kirchhoff
/// material on `DFRGRD1`. `STRESS` is the Cauchy stress `1/J F S F^T` of the
/// second Piola-Kirchhoff stress `S = C E` of the Green-Lagrange strain,
/// `DDSDDE` its spatial tangent.
pub fn stvenant_kirchhoff(c: &mut UmatCall<'_>) -> Result<(), String> {
    use elmer_physics_elastic::post::{matmul3, transpose3};
    if c.props.len() < 3 {
        return Err("stvenant_kirchhoff: 3 material constants are needed".into());
    }
    let (ndi, nshr, ntens) = (c.ndi, c.nshr, c.ntens);
    let f = c.dfgrd1;
    let bm = matmul3(&f, &transpose3(&f));
    let cm = matmul3(&transpose3(&f), &f);
    let mut strain = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            strain[i][j] = 0.5 * (cm[i][j] - if i == j { 1.0 } else { 0.0 });
        }
    }
    let mut sv = [0.0; 6];
    for i in 0..ndi {
        sv[i] = strain[i][i];
    }
    for i in 0..nshr {
        let (a, b) = SHEAR[i];
        sv[ndi + i] = strain[a][b] + strain[b][a];
    }
    let det = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
        + f[0][1] * (f[1][2] * f[2][0] - f[1][0] * f[2][2])
        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
    let e = c.props[1];
    let nu = c.props[2];
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let d: &mut [[f64; 6]; 6] = &mut *c.ddsdde;
    *d = [[0.0; 6]; 6];
    for i in 0..ndi {
        for j in 0..ndi {
            d[i][j] = lambda;
        }
    }
    for i in 0..ntens {
        d[i][i] += mu;
    }
    for i in 0..ndi {
        d[i][i] += mu;
    }
    let mut s2 = [0.0; 6];
    for (i, v) in s2.iter_mut().enumerate().take(ntens) {
        for j in 0..ntens {
            *v += d[i][j] * sv[j];
        }
    }
    let mut s = [[0.0; 3]; 3];
    for i in 0..ndi {
        let b = sym_basis(i);
        for a in 0..3 {
            for k in 0..3 {
                s[a][k] += s2[i] * b[a][k];
            }
        }
    }
    for i in 0..nshr {
        let b = sym_basis(ndi + i);
        for a in 0..3 {
            for k in 0..3 {
                s[a][k] += 2.0 * s2[ndi + i] * b[a][k];
            }
        }
    }
    let fs = push(&f, &s);
    let mut sigma = [[0.0; 3]; 3];
    for a in 0..3 {
        for k in 0..3 {
            sigma[a][k] = 1.0 / det * fs[a][k];
        }
    }
    for i in 0..ndi {
        c.stress[i] = sigma[i][i];
    }
    for i in 0..nshr {
        let (a, b) = SHEAR[i];
        c.stress[ndi + i] = sigma[a][b];
    }
    *d = [[0.0; 6]; 6];
    // WorkMat = LambdaLame * 1/DetDefG * B.
    let scale = lambda * 1.0 / det;
    let mut wm = [[0.0; 3]; 3];
    for a in 0..3 {
        for k in 0..3 {
            wm[a][k] = scale * bm[a][k];
        }
    }
    let add_col = |d: &mut [[f64; 6]; 6], col: usize, wm: &[[f64; 3]; 3], one: bool| {
        for j in 0..ndi {
            d[j][col] += if one { 1.0 * wm[j][j] } else { wm[j][j] };
        }
        for j in 0..nshr {
            let (a, b) = SHEAR[j];
            d[ndi + j][col] += if one { 1.0 * wm[a][b] } else { wm[a][b] };
        }
    };
    for i in 0..ndi {
        add_col(d, i, &wm, false);
    }
    let scale2 = 2.0 * mu * 1.0 / det;
    for i in 0..ntens {
        let p = push(&f, &sym_basis(i));
        let mut w = [[0.0; 3]; 3];
        for a in 0..3 {
            for k in 0..3 {
                w[a][k] = scale2 * p[a][k];
            }
        }
        add_col(d, i, &w, i >= ndi);
    }
    Ok(())
}

/// The integration-point state of one element: the committed values (read)
/// and the trial ones (written), `ntens`, 3 and `nstatev_max` per point.
pub struct IpState<'a> {
    pub stress0: &'a [f64],
    pub energy0: &'a [f64],
    pub state0: &'a [f64],
    pub stress: &'a mut [f64],
    pub energy: &'a mut [f64],
    pub state: &'a mut [f64],
    /// `UmatState % Dofs` (`MaxStateV`), 0 without the variable.
    pub max_statev: usize,
}

/// What `LocalMatrixWithUMAT` reads for one element.
pub struct UmatInput<'a> {
    pub etype: ElementType,
    pub xyz: &'a [[f64; 3]],
    pub cdim: usize,
    /// The current total displacement and the previous step's, per node.
    pub disp: &'a [[f64; 3]],
    pub prev: &'a [[f64; 3]],
    pub load: &'a [[f64; 3]],
    pub inertial: Option<&'a [[f64; 3]]>,
    pub density: &'a [f64],
    /// `LocalTemperature`: the temperature itself (not its offset).
    pub temperature: &'a [f64],
    pub props: &'a [f64],
    pub nstatev: usize,
    pub cmname: &'a str,
    pub umat: UmatFn,
    pub time: f64,
    pub dt: f64,
    /// 1-based `ElementIndex` passed as `NOEL`.
    pub element_index: usize,
    /// `Large Deflection`: Green-Lagrange strains and the spatial
    /// formulation of the stress the UMAT returns.
    pub large: bool,
    /// `AxialSymmetry`: the hoop strain `u_r / r` and the weight `r`.
    pub axi: bool,
}

/// `LocalMatrixWithUMAT` with small strains: per point, the UMAT on
/// `STRAN = eps(u_prev)`, `DSTRAN = eps(u) - eps(u_prev)`, then
/// `STIFF += s B^T DDSDDE B`, `FORCE += -B^T STRESS s + loads`,
/// `EXTFORCE += loads`. `stiff` is `ndof x ndof` row-major, `ndof = n cdim`.
#[allow(clippy::needless_range_loop)]
pub fn umat_local(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    inp: &UmatInput<'_>,
    st: &mut IpState<'_>,
    stiff: &mut [f64],
    force: &mut [f64],
    ext: &mut [f64],
) -> Result<(), String> {
    let n = inp.etype.n_nodes();
    let cdim = inp.cdim;
    let dofs = cdim;
    let ndof = n * dofs;
    let (ndi, nshr) = if cdim == 2 { (3, 1) } else { (3, 3) };
    let ntens = ndi + nshr;
    stiff.iter_mut().for_each(|v| *v = 0.0);
    force.iter_mut().for_each(|v| *v = 0.0);
    ext.iter_mut().for_each(|v| *v = 0.0);
    let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut be = BasisEval::default();
    let mut b = vec![[0.0; 6]; ndof];
    for (t, qp) in points.iter().enumerate() {
        basis.eval(inp.etype, inp.xyz, qp.uvw, &mut be)?;
        let bas = &be.basis;
        let db = &be.dbasis_dx;
        let mut s = be.det_j * qp.weight;
        let mut r = 0.0;
        if inp.axi {
            r = (0..n).map(|k| bas[k] * inp.xyz[k][0]).sum();
            s *= r;
        }
        let mut fv = [0.0; 3];
        let mut iv = [0.0; 3];
        for i in 0..cdim {
            fv[i] = (0..n).map(|k| inp.load[k][i] * bas[k]).sum();
            if let Some(w) = inp.inertial {
                iv[i] = (0..n).map(|k| w[k][i] * bas[k]).sum();
            }
        }
        let temp: f64 = (0..n).map(|k| inp.temperature[k] * bas[k]).sum();
        let density: f64 = (0..n).map(|k| inp.density[k] * bas[k]).sum();
        // Grad(1:cdim,1:cdim) = MATMUL(u(1:cdim,1:nd), dBasis(1:nd,1:cdim)).
        let mut grad = [[0.0; 3]; 3];
        let mut grad0 = [[0.0; 3]; 3];
        for a in 0..cdim {
            for c in 0..cdim {
                grad[a][c] = (0..n).map(|k| inp.disp[k][a] * db[k][c]).sum();
                grad0[a][c] = (0..n).map(|k| inp.prev[k][a] * db[k][c]).sum();
            }
        }
        if inp.axi {
            grad[2][2] = 1.0 / r * (0..n).map(|k| inp.disp[k][0] * bas[k]).sum::<f64>();
            grad0[2][2] = 1.0 / r * (0..n).map(|k| inp.prev[k][0] * bas[k]).sum::<f64>();
        }
        let mut defg = ident;
        let mut defg0 = ident;
        let mut strain = [[0.0; 3]; 3];
        let mut strain0 = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                defg[i][j] += grad[i][j];
                defg0[i][j] += grad0[i][j];
                strain[i][j] = 0.5 * (grad[i][j] + grad[j][i]);
                strain0[i][j] = 0.5 * (grad0[i][j] + grad0[j][i]);
            }
        }
        let (mut det, mut inv) = (1.0, ident);
        if inp.large {
            use elmer_physics_elastic::post::{invert_matrix, matmul3, transpose3};
            let f = &defg;
            det = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
                + f[0][1] * (f[1][2] * f[2][0] - f[1][0] * f[2][2])
                + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
            inv = defg;
            invert_matrix(&mut inv, 3).map_err(|e| e.to_string())?;
            let gtg = matmul3(&transpose3(&grad), &grad);
            let gtg0 = matmul3(&transpose3(&grad0), &grad0);
            for i in 0..3 {
                for j in 0..3 {
                    strain[i][j] += 0.5 * gtg[i][j];
                    strain0[i][j] += 0.5 * gtg0[i][j];
                }
            }
        }
        let stran = [
            strain0[0][0],
            strain0[1][1],
            strain0[2][2],
            2.0 * strain0[0][1],
            2.0 * strain0[0][2],
            2.0 * strain0[1][2],
        ];
        let dstran = [
            strain[0][0] - strain0[0][0],
            strain[1][1] - strain0[1][1],
            strain[2][2] - strain0[2][2],
            2.0 * (strain[0][1] - strain0[0][1]),
            2.0 * (strain[0][2] - strain0[0][2]),
            2.0 * (strain[1][2] - strain0[1][2]),
        ];
        let mut energy = [st.energy0[3 * t], st.energy0[3 * t + 1], st.energy0[3 * t + 2]];
        let mut stress = [0.0; 6];
        stress[..ntens].copy_from_slice(&st.stress0[ntens * t..ntens * (t + 1)]);
        let mut statev = vec![0.0; inp.nstatev];
        if inp.nstatev > 0 {
            statev.copy_from_slice(&st.state0[st.max_statev * t..st.max_statev * t + inp.nstatev]);
        }
        let mut ddsdde = [[0.0; 6]; 6];
        (inp.umat)(&mut UmatCall {
            stress: &mut stress[..ntens],
            statev: &mut statev,
            ddsdde: &mut ddsdde,
            energy: &mut energy,
            stran: &stran[..ntens],
            dstran: &dstran[..ntens],
            time: [inp.time - inp.dt; 2],
            dtime: inp.dt,
            temp,
            cmname: inp.cmname,
            ndi,
            nshr,
            ntens,
            props: inp.props,
            dfgrd0: defg0,
            dfgrd1: defg,
            noel: inp.element_index,
            npt: t + 1,
        })?;
        st.energy[3 * t..3 * t + 3].copy_from_slice(&energy);
        st.stress[ntens * t..ntens * (t + 1)].copy_from_slice(&stress[..ntens]);
        if inp.nstatev > 0 {
            st.state[st.max_statev * t..st.max_statev * t + inp.nstatev].copy_from_slice(&statev);
        }
        if inp.large {
            large_terms(&LargeIp { bas, db, s, r, axi: inp.axi, n, cdim, nshr, ntens, det, inv: &inv, defg: &defg, stress: &stress, ddsdde: &ddsdde, fv: &fv, iv: &iv, density }, stiff, force, ext);
            continue;
        }
        // B (ntens x ndof), stored by column.
        for col in b.iter_mut() {
            *col = [0.0; 6];
        }
        for p in 0..n {
            if inp.axi {
                b[p * dofs][0] = db[p][0];
                b[p * dofs + 1][1] = db[p][1];
                b[p * dofs][2] = 1.0 / r * bas[p];
                b[p * dofs][3] = db[p][1];
                b[p * dofs + 1][3] = db[p][0];
                continue;
            }
            for i in 0..cdim {
                b[p * dofs + i][i] = db[p][i];
            }
            for i in 0..nshr {
                match i {
                    0 => {
                        b[p * dofs][ndi] = db[p][1];
                        b[p * dofs + 1][ndi] = db[p][0];
                    }
                    1 => {
                        b[p * dofs][ndi + 1] = db[p][2];
                        b[p * dofs + 2][ndi + 1] = db[p][0];
                    }
                    _ => {
                        b[p * dofs + 1][ndi + 2] = db[p][2];
                        b[p * dofs + 2][ndi + 2] = db[p][1];
                    }
                }
            }
        }
        // StrainEnergyDensity: STIFF += s * MATMUL(TRANSPOSE(B), MATMUL(D, B)).
        let mut dbm = vec![[0.0; 6]; ndof];
        for c in 0..ndof {
            for k in 0..ntens {
                let mut acc = 0.0;
                for l in 0..ntens {
                    acc += ddsdde[k][l] * b[c][l];
                }
                dbm[c][k] = acc;
            }
        }
        for r in 0..ndof {
            for c in 0..ndof {
                let mut acc = 0.0;
                for k in 0..ntens {
                    acc += b[r][k] * dbm[c][k];
                }
                stiff[r * ndof + c] += s * acc;
            }
        }
        // FORCE = FORCE - MATMUL(TRANSPOSE(B), STRESS) * s.
        for c in 0..ndof {
            let mut acc = 0.0;
            for k in 0..ntens {
                acc += b[c][k] * stress[k];
            }
            force[c] -= acc * s;
        }
        for p in 0..n {
            for i in 0..cdim {
                let v = (bas[p] * fv[i] + bas[p] * iv[i] * density) * s;
                ext[dofs * p + i] += v;
                force[dofs * p + i] += v;
            }
        }
    }
    Ok(())
}

/// One element of `GenerateStressVariable` / `GenerateStrainVariable`
/// (`NodalProjectorAssemble` on `GaussPoints(Element)` with
/// `ElasticUmatStressAtIP` and `ElasticStrainAtIP`): the projection mass, the
/// strain of the displacement (Green-Lagrange under `Large Deflection`) and
/// the stress `UmatStress` holds at the points (`ip_stress`, `ntens` per
/// point, read in the order `1 2 3 4 6 5`).
#[allow(clippy::too_many_arguments)]
pub fn umat_projection_local(
    basis: &dyn ElementBasis,
    points: &[QuadPoint],
    etype: ElementType,
    xyz: &[[f64; 3]],
    disp: &[[f64; 3]],
    ip_stress: &[f64],
    ntens: usize,
    ncomp: usize,
    large: bool,
    axi: bool,
) -> Result<elmer_physics_elastic::post::StressPostLocal, String> {
    use elmer_physics_elastic::post::{add_tensor, vector_to_tensor};
    let n = etype.n_nodes();
    let mut out = elmer_physics_elastic::post::StressPostLocal {
        mass: vec![0.0; n * n],
        strain: vec![0.0; n * ncomp],
        stress: vec![0.0; n * ncomp],
    };
    if ip_stress.len() != ntens * points.len() {
        return Err("UmatStress has a different number of points than the projection".into());
    }
    let ind = if axi { [0usize, 2, 1, 3, 4, 5] } else { [0usize, 1, 2, 3, 5, 4] };
    let mut be = BasisEval::default();
    for (t, qp) in points.iter().enumerate() {
        basis.eval(etype, xyz, qp.uvw, &mut be)?;
        let bas = &be.basis[..n];
        let db = &be.dbasis_dx;
        let mut w = qp.weight * be.det_j;
        let mut g = [[0.0; 3]; 3];
        for a in 0..3 {
            for c in 0..3 {
                g[a][c] = (0..n).map(|k| disp[k][a] * db[k][c]).sum();
            }
        }
        if axi {
            let r: f64 = (0..n).map(|k| bas[k] * xyz[k][0]).sum();
            w *= r;
            g[2][2] = (0..n).map(|k| disp[k][0] * bas[k]).sum::<f64>() / r;
        }
        let mut t1 = [[0.0; 3]; 3];
        for a in 0..3 {
            for c in 0..3 {
                t1[a][c] = (g[c][a] + g[a][c]) / 2.0;
            }
        }
        if large {
            for a in 0..3 {
                for c in 0..3 {
                    let gtg: f64 = (0..3).map(|k| g[k][a] * g[k][c]).sum();
                    t1[a][c] += gtg / 2.0;
                }
            }
        }
        for p in 0..n {
            for q in 0..n {
                out.mass[p * n + q] += w * bas[q] * bas[p];
            }
        }
        add_tensor(&mut out.strain, bas, ncomp, w, &t1);
        let mut comp = [0.0; 6];
        for i in 0..ncomp {
            comp[i] = ip_stress[ntens * t + ind[i]];
        }
        add_tensor(&mut out.stress, bas, ncomp, w, &vector_to_tensor(&comp[..ncomp], ncomp));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `linear_isotropic` is Hooke's law in Voigt form with engineering shears.
    #[test]
    fn linear_isotropic_is_hooke() {
        let (e, nu) = (70e9, 0.3);
        let mut stress = [0.0; 4];
        let mut ddsdde = [[0.0; 6]; 6];
        let mut energy = [0.0; 3];
        let dstran = [1e-3, -2e-4, 0.0, 5e-4];
        let mut statev: [f64; 0] = [];
        linear_isotropic(&mut UmatCall {
            stress: &mut stress,
            statev: &mut statev,
            ddsdde: &mut ddsdde,
            energy: &mut energy,
            stran: &[0.0; 4],
            dstran: &dstran,
            time: [0.0; 2],
            dtime: 1.0,
            temp: 0.0,
            cmname: "linear isotropic",
            ndi: 3,
            nshr: 1,
            ntens: 4,
            props: &[1.0, e, nu],
            dfgrd0: [[0.0; 3]; 3],
            dfgrd1: [[0.0; 3]; 3],
            noel: 1,
            npt: 1,
        })
        .unwrap();
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let tr = dstran[0] + dstran[1];
        assert!((stress[0] - (lambda * tr + 2.0 * mu * dstran[0])).abs() < 1e-3 * stress[0].abs());
        assert!((stress[3] - mu * dstran[3]).abs() < 1e-9 * stress[3].abs());
        assert_eq!(ddsdde[3][3], mu);
    }
}

/// One point of the `Large Deflection` branch of `LocalMatrixWithUMAT`.
struct LargeIp<'a> {
    bas: &'a [f64],
    db: &'a [[f64; 3]],
    s: f64,
    r: f64,
    axi: bool,
    n: usize,
    cdim: usize,
    nshr: usize,
    ntens: usize,
    det: f64,
    inv: &'a [[f64; 3]; 3],
    defg: &'a [[f64; 3]; 3],
    stress: &'a [f64; 6],
    ddsdde: &'a [[f64; 6]; 6],
    fv: &'a [f64; 3],
    iv: &'a [f64; 3],
    density: f64,
}

/// ElasticSolve.F90:3325-3640 without axial symmetry: the first
/// Piola-Kirchhoff stress `J sigma F^-T` of the UMAT's Cauchy stress and,
/// per trial function, its linearisation through `DDSDDE`.
#[allow(clippy::needless_range_loop)]
fn large_terms(ip: &LargeIp<'_>, stiff: &mut [f64], force: &mut [f64], ext: &mut [f64]) {
    use elmer_physics_elastic::post::{matmul3, transpose3};
    let (n, cdim, ndof) = (ip.n, ip.cdim, ip.n * ip.cdim);
    let voigt = |v: &[f64]| -> [[f64; 3]; 3] {
        let mut t = [[0.0; 3]; 3];
        t[0][0] = v[0];
        t[1][1] = v[1];
        t[2][2] = v[2];
        for k in 0..ip.nshr {
            let (a, b) = SHEAR[k];
            let x = 2.0 * v[3 + k] * 0.5;
            t[a][b] = x;
            t[b][a] = x;
        }
        t
    };
    let stress = voigt(&ip.stress[..ip.ntens]);
    let inv_t = transpose3(ip.inv);
    let st_inv_t = matmul3(&stress, &inv_t);
    let mut stress1 = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            stress1[a][b] = ip.det * st_inv_t[a][b];
        }
    }
    for p in 0..n {
        for i in 0..cdim {
            let mut g = [[0.0; 3]; 3];
            if ip.axi {
                if i == 0 {
                    g[0][0] = ip.db[p][0];
                    g[0][1] = ip.db[p][1];
                    g[2][2] = 1.0 / ip.r * ip.bas[p];
                } else {
                    g[1][0] = ip.db[p][0];
                    g[1][1] = ip.db[p][1];
                }
            } else {
                g[i] = ip.db[p];
            }
            let ginv = matmul3(&g, ip.inv);
            let sg = matmul3(&stress, &ginv);
            let tr = sg[0][0] + sg[1][1] + sg[2][2];
            let m = matmul3(&inv_t, &matmul3(&transpose3(&g), &st_inv_t));
            let mut wt2 = [[0.0; 3]; 3];
            let mut wt1 = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    wt2[a][b] = ip.det * tr * inv_t[a][b] - ip.det * m[a][b];
                    wt1[a][b] = ip.det * ginv[a][b];
                }
            }
            let w1t = transpose3(&wt1);
            for a in 0..3 {
                for b in 0..3 {
                    wt1[a][b] = 0.5 * (wt1[a][b] + w1t[a][b]);
                }
            }
            let mut wv1 = [0.0; 6];
            wv1[0] = wt1[0][0];
            wv1[1] = wt1[1][1];
            wv1[2] = wt1[2][2];
            for k in 0..ip.nshr {
                let (a, b) = SHEAR[k];
                wv1[3 + k] = wt1[a][b] + wt1[b][a];
            }
            let mut wv2 = [0.0; 6];
            for k in 0..ip.ntens {
                let mut acc = 0.0;
                for l in 0..ip.ntens {
                    acc += ip.ddsdde[l][k] * wv1[l];
                }
                wv2[k] = acc;
            }
            let wt3 = voigt(&wv2[..ip.ntens]);
            let fw = matmul3(ip.defg, &wt3);
            let mut ds1 = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    ds1[a][b] = wt2[a][b] + fw[a][b];
                }
            }
            let row = cdim * p + i;
            if ip.axi {
                for q in 0..n {
                    for j in 0..cdim {
                        let v = if j == 0 {
                            ip.db[q][0] * ds1[0][0] + ip.db[q][1] * ds1[0][1] + 1.0 / ip.r * ip.bas[q] * ds1[2][2]
                        } else {
                            ip.db[q][0] * ds1[1][0] + ip.db[q][1] * ds1[1][1]
                        };
                        stiff[row * ndof + cdim * q + j] += v * ip.s;
                    }
                }
                let load = (ip.bas[p] * ip.fv[i] * ip.det + ip.bas[p] * ip.iv[i] * ip.density) * ip.s;
                force[row] += load;
                ext[row] += load;
                let sub = if i == 0 {
                    ip.db[p][0] * stress1[0][0] + ip.db[p][1] * stress1[0][1] + 1.0 / ip.r * ip.bas[p] * stress1[2][2]
                } else {
                    ip.db[p][0] * stress1[1][0] + ip.db[p][1] * stress1[1][1]
                };
                force[row] -= sub * ip.s;
                continue;
            }
            for q in 0..n {
                for j in 0..cdim {
                    let dot = ip.db[q][0] * ds1[j][0] + ip.db[q][1] * ds1[j][1] + ip.db[q][2] * ds1[j][2];
                    stiff[row * ndof + cdim * q + j] += dot * ip.s;
                }
            }
            let dot = ip.db[p][0] * stress1[i][0] + ip.db[p][1] * stress1[i][1] + ip.db[p][2] * stress1[i][2];
            force[row] += (ip.bas[p] * ip.fv[i] * ip.det + ip.bas[p] * ip.iv[i] * ip.density - dot) * ip.s;
            ext[row] += (ip.bas[p] * ip.fv[i] * ip.det + ip.bas[p] * ip.iv[i] * ip.density) * ip.s;
        }
    }
}
