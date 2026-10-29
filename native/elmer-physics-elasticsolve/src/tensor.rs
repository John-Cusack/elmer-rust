//! Fixed 3x3 tensor helpers with the semantics of the Fortran intrinsics and
//! `LinearAlgebra` routines `ElasticSolve.F90` uses: products run over the whole
//! 3x3 array, while `InvertMatrix`, `TRACE` and `DDOTPROD` take the leading
//! `n x n` block.

pub type T3 = [[f64; 3]; 3];

pub const ZERO: T3 = [[0.0; 3]; 3];

/// `Identity(i,i) = 1` for `i <= dim`.
pub fn identity(dim: usize) -> T3 {
    let mut m = ZERO;
    for (i, row) in m.iter_mut().enumerate().take(dim) {
        row[i] = 1.0;
    }
    m
}

/// `MATMUL(a, b)`.
pub fn mul(a: &T3, b: &T3) -> T3 {
    let mut c = ZERO;
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    c
}

pub fn transpose(a: &T3) -> T3 {
    let mut t = ZERO;
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

pub fn add(a: &T3, b: &T3) -> T3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}

pub fn sub(a: &T3, b: &T3) -> T3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] - b[i][j]))
}

pub fn scale(a: &T3, s: f64) -> T3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] * s))
}

/// `DDOTPROD(a, b, n)`: the double contraction over the leading `n x n` block.
pub fn ddot(a: &T3, b: &T3, n: usize) -> f64 {
    let mut s = 0.0;
    for i in 0..n {
        for j in 0..n {
            s += a[i][j] * b[i][j];
        }
    }
    s
}

/// `TRACE(a, n)`.
pub fn trace(a: &T3, n: usize) -> f64 {
    (0..n).map(|i| a[i][i]).sum()
}

/// Determinant of the leading `n x n` block, as `ElasticSolve` writes it out
/// for `n = 1, 2, 3`.
pub fn det(a: &T3, n: usize) -> f64 {
    match n {
        1 => a[0][0],
        2 => a[0][0] * a[1][1] - a[0][1] * a[1][0],
        _ => {
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        }
    }
}

/// `InvertMatrix(a, n)`: inverts the leading `n x n` block in place and leaves
/// the rest of the array alone. `None` for a singular block.
pub fn invert(a: &T3, n: usize) -> Option<T3> {
    let mut o = *a;
    match n {
        1 => {
            if a[0][0] == 0.0 {
                return None;
            }
            o[0][0] = 1.0 / a[0][0];
        }
        2 => {
            let d = det(a, 2);
            if d == 0.0 {
                return None;
            }
            o[0][0] = a[1][1] / d;
            o[1][1] = a[0][0] / d;
            o[0][1] = -a[0][1] / d;
            o[1][0] = -a[1][0] / d;
        }
        _ => {
            let d = det(a, 3);
            if d == 0.0 {
                return None;
            }
            for i in 0..3 {
                for j in 0..3 {
                    let (r1, r2) = ((j + 1) % 3, (j + 2) % 3);
                    let (c1, c2) = ((i + 1) % 3, (i + 2) % 3);
                    o[i][j] = (a[r1][c1] * a[r2][c2] - a[r1][c2] * a[r2][c1]) / d;
                }
            }
        }
    }
    Some(o)
}

pub fn symmetrize(g: &T3) -> T3 {
    std::array::from_fn(|i| std::array::from_fn(|j| (g[i][j] + g[j][i]) / 2.0))
}
