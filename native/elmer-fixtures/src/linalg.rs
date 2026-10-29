//! Dense helpers the plugins call from Elmer's `LinearAlgebra` module.

/// `InvertMatrix`: Gauss-Jordan with partial pivoting, row-major `n x n`.
pub fn invert(a: &mut [f64], n: usize) -> Result<(), String> {
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i * n + c].abs().total_cmp(&a[j * n + c].abs())).unwrap_or(c);
        if a[p * n + c] == 0.0 {
            return Err("InvertMatrix: singular matrix".into());
        }
        for k in 0..n {
            a.swap(c * n + k, p * n + k);
            inv.swap(c * n + k, p * n + k);
        }
        let d = a[c * n + c];
        for k in 0..n {
            a[c * n + k] /= d;
            inv[c * n + k] /= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r * n + c];
                if f != 0.0 {
                    for k in 0..n {
                        a[r * n + k] -= f * a[c * n + k];
                        inv[r * n + k] -= f * inv[c * n + k];
                    }
                }
            }
        }
    }
    a.copy_from_slice(&inv);
    Ok(())
}

/// `LCondensate(N, Nb, K, F)`: eliminates the last `nb` of `n + nb` dofs
/// (local bubbles) from the row-major `k` and `f`, leaving the leading
/// `n x n` block and `n` entries: `F -= Klb Kbb^-1 Fb`,
/// `K -= Klb Kbb^-1 Kbl`, each product taken inner first as `MATMUL` nests.
pub fn condensate(n: usize, nb: usize, k: &mut Vec<f64>, f: &mut Vec<f64>) -> Result<(), String> {
    let m = n + nb;
    if nb == 0 {
        return Ok(());
    }
    let mut kbb: Vec<f64> = (0..nb).flat_map(|i| (0..nb).map(move |j| (i, j))).map(|(i, j)| k[(n + i) * m + n + j]).collect();
    invert(&mut kbb, nb)?;
    let kbl = |i: usize, j: usize| k[(n + i) * m + j];
    let klb = |i: usize, j: usize| k[i * m + n + j];
    let fb: Vec<f64> = f[n..m].to_vec();
    // MATMUL(Kbb, Fb), MATMUL(Kbb, Kbl)
    let kfb: Vec<f64> = (0..nb).map(|i| (0..nb).fold(0.0, |s, j| s + kbb[i * nb + j] * fb[j])).collect();
    let kkbl: Vec<f64> =
        (0..nb).flat_map(|i| (0..n).map(move |j| (i, j))).map(|(i, j)| (0..nb).fold(0.0, |s, l| s + kbb[i * nb + l] * kbl(l, j))).collect();
    let mut nf = vec![0.0; n];
    let mut nk = vec![0.0; n * n];
    for i in 0..n {
        nf[i] = f[i] - (0..nb).fold(0.0, |s, l| s + klb(i, l) * kfb[l]);
        for j in 0..n {
            nk[i * n + j] = k[i * m + j] - (0..nb).fold(0.0, |s, l| s + klb(i, l) * kkbl[l * n + j]);
        }
    }
    *k = nk;
    *f = nf;
    Ok(())
}
