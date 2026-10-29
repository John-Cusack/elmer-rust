//! Static condensation of interior (bubble) basis functions: `CondensateP`
//! and the transient `CondensatePTransient` of `MatrixAssembly.F90`.
//!
//! The local numbering has the `nl` retained dofs first and the bubble dofs
//! after, all fields interleaved per basis function, as the kernels build it.

use elmer_core::LocalSystem;

use crate::error::ThermoElectricError;

/// What the condensation of one element eliminated, enough to recover the
/// bubble values from the solved retained dofs (`BubbleCur` of
/// `CondensatePTransient`).
#[derive(Clone, Debug)]
pub struct BubbleRecovery {
    nl: usize,
    nb: usize,
    /// `Kbb^-1`, row-major `nb x nb`.
    kbb_inv: Vec<f64>,
    /// `Kbl`, row-major `nb x nl`.
    kbl: Vec<f64>,
    fb: Vec<f64>,
}

impl BubbleRecovery {
    /// `Kbb^-1 (Fb - Kbl x)`: the bubble dofs given the retained ones.
    pub fn bubbles(&self, x: &[f64]) -> Vec<f64> {
        let (nl, nb) = (self.nl, self.nb);
        let r: Vec<f64> = (0..nb)
            .map(|i| {
                let mut s = 0.0;
                for j in 0..nl {
                    s += self.kbl[i * nl + j] * x[j];
                }
                self.fb[i] - s
            })
            .collect();
        (0..nb)
            .map(|i| {
                let mut s = 0.0;
                for j in 0..nb {
                    s += self.kbb_inv[i * nb + j] * r[j];
                }
                s
            })
            .collect()
    }
}

/// `InvertMatrix`: Gauss-Jordan with partial pivoting, in place.
fn invert(a: &mut [f64], n: usize) -> Result<(), ThermoElectricError> {
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for col in 0..n {
        let piv = (col..n)
            .max_by(|&x, &y| a[x * n + col].abs().total_cmp(&a[y * n + col].abs()))
            .unwrap_or(col);
        let d = a[piv * n + col];
        if d == 0.0 || !d.is_finite() {
            return Err(ThermoElectricError::Unsupported(
                "singular bubble block in the static condensation".into(),
            ));
        }
        if piv != col {
            for k in 0..n {
                a.swap(piv * n + k, col * n + k);
                inv.swap(piv * n + k, col * n + k);
            }
        }
        for k in 0..n {
            a[col * n + k] /= d;
            inv[col * n + k] /= d;
        }
        for r in 0..n {
            if r != col {
                let f = a[r * n + col];
                if f != 0.0 {
                    for k in 0..n {
                        a[r * n + k] -= f * a[col * n + k];
                        inv[r * n + k] -= f * inv[col * n + k];
                    }
                }
            }
        }
    }
    a.copy_from_slice(&inv);
    Ok(())
}

/// Adds the backward-Euler time terms over the full retained+bubble block,
/// `K += M/dt`, `F += M xprev/dt` (the head of `CondensatePTransient`), and
/// drops the mass matrix. `xprev` holds all `local.n` dofs of the previous
/// timestep, the bubbles' from the caller's own history.
pub fn add_backward_euler(local: &mut LocalSystem, dt: f64, xprev: &[f64]) {
    let n = local.n;
    let Some(m) = local.mass.take() else { return };
    for i in 0..n {
        let mut mx = 0.0;
        for j in 0..n {
            local.stiff[i * n + j] += m[i * n + j] / dt;
            mx += m[i * n + j] * xprev[j];
        }
        local.force[i] += mx / dt;
    }
}

/// `CondensateP( nl, nb, K, F )`: eliminates the last `nb` dofs of `local`,
/// which shrinks to `nl` dofs. Returns the recovery of the eliminated dofs.
pub fn condense_bubbles(local: &mut LocalSystem, nl: usize) -> Result<BubbleRecovery, ThermoElectricError> {
    let n = local.n;
    let nb = n - nl;
    let k = &local.stiff;
    let mut kbb: Vec<f64> = (0..nb * nb).map(|x| k[(nl + x / nb) * n + nl + x % nb]).collect();
    let kbl: Vec<f64> = (0..nb * nl).map(|x| k[(nl + x / nl) * n + x % nl]).collect();
    let klb: Vec<f64> = (0..nl * nb).map(|x| k[(x / nb) * n + nl + x % nb]).collect();
    let fb: Vec<f64> = local.force[nl..].to_vec();
    invert(&mut kbb, nb)?;

    // Kbb^-1 Fb and Kbb^-1 Kbl.
    let kinv_fb: Vec<f64> = (0..nb)
        .map(|i| {
            let mut s = 0.0;
            for j in 0..nb {
                s += kbb[i * nb + j] * fb[j];
            }
            s
        })
        .collect();
    let mut kinv_kbl = vec![0.0; nb * nl];
    for i in 0..nb {
        for c in 0..nl {
            let mut s = 0.0;
            for j in 0..nb {
                s += kbb[i * nb + j] * kbl[j * nl + c];
            }
            kinv_kbl[i * nl + c] = s;
        }
    }
    let mut stiff = vec![0.0; nl * nl];
    for r in 0..nl {
        for c in 0..nl {
            let mut s = 0.0;
            for j in 0..nb {
                s += klb[r * nb + j] * kinv_kbl[j * nl + c];
            }
            stiff[r * nl + c] = k[r * n + c] - s;
        }
    }
    let mut force = vec![0.0; nl];
    for r in 0..nl {
        let mut s = 0.0;
        for j in 0..nb {
            s += klb[r * nb + j] * kinv_fb[j];
        }
        force[r] = local.force[r] - s;
    }
    local.n = nl;
    local.stiff = stiff;
    local.force = force;
    local.mass = None;
    local.damp = None;
    Ok(BubbleRecovery {
        nl,
        nb,
        kbb_inv: kbb,
        kbl,
        fb,
    })
}
