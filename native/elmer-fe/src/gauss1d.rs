//! The 1-D Gauss-Legendre points and weights of Integration.F90 (`Points`,
//! `Weights`, MAXN = 13).
//!
//! Elmer computes them once at start-up (`ComputeGaussPoints1D`, Integration.F90:1296):
//! Legendre coefficients from the three-term recurrence, roots of the
//! even part from LAPACK `DGEEV` on a companion matrix, Newton polishing that
//! only replaces a root when it lowers the residual, then weights
//! renormalised to sum to 2. The `DGEEV` eigenvalues fix the last bits, so the
//! table is the pristine binary's own `Points`/`Weights` arrays, captured as
//! bit patterns (see `tests/quadrature_bits.rs`). Roots come in `+x, -x` pairs,
//! with the zero root last for odd n, which fixes the order of every tensor rule.

pub(crate) const MAXN: usize = 13;

#[rustfmt::skip]
const POINTS_BITS: [&[u64]; MAXN] = [
    &[0x0000000000000000],
    &[0x3fe279a74590331c, 0xbfe279a74590331c],
    &[0x3fe8c97ef43f7248, 0xbfe8c97ef43f7248, 0x0000000000000000],
    &[0x3fd5c23fd9dd3dfc, 0xbfd5c23fd9dd3dfc, 0x3feb8e6dbcf63985, 0xbfeb8e6dbcf63985],
    &[0x3fe13b23fd99b705, 0xbfe13b23fd99b705, 0x3fecff6ce0533a6a, 0xbfecff6ce0533a6a, 0x0000000000000000],
    &[0x3fce8b12d03675c5, 0xbfce8b12d03675c5, 0x3fedd6ca4e80a01e, 0xbfedd6ca4e80a01e, 0x3fe528a09655c95f, 0xbfe528a09655c95f],
    &[0x3fd9f95df119fd62, 0xbfd9f95df119fd62, 0x3fe7ba9f9be3a1d5, 0xbfe7ba9f9be3a1d5, 0x3fee5f178e7c622a, 0xbfee5f178e7c622a, 0x0000000000000000],
    &[0x3feebab1cb0acc65, 0xbfeebab1cb0acc65, 0x3fe97e4ab249f421, 0xbfe97e4ab249f421, 0x3fe0d129583284b3, 0xbfe0d129583284b3, 0x3fc77ac94f3c7344, 0xbfc77ac94f3c7344],
    &[0x3feefb2b2ebf2106, 0xbfeefb2b2ebf2106, 0x3feac0c44f0d029b, 0xbfeac0c44f0d029b, 0x3fe3a0bd2077fd8d, 0xbfe3a0bd2077fd8d, 0x3fd4c0916e48aa66, 0xbfd4c0916e48aa66, 0x0000000000000000],
    &[0x3fef2a3e062af2da, 0xbfef2a3e062af2da, 0x3febae995e9cb2f0, 0xbfebae995e9cb2f0, 0x3fe5bdb9228de192, 0xbfe5bdb9228de192, 0x3fdbbcc009016add, 0xbfdbbcc009016add, 0x3fc30e507891e27a, 0xbfc30e507891e27a],
    &[0x3fef4da62fd7e9a9, 0xbfef4da62fd7e9a9, 0x3fec62d11af0477b, 0xbfec62d11af0477b, 0x3fe75d67bd219447, 0xbfe75d67bd219447, 0x3fe09c6f7c4d8ce3, 0xbfe09c6f7c4d8ce3, 0x3fd14031efeb42c1, 0xbfd14031efeb42c1, 0x0000000000000000],
    &[0x3fef68f1d8e42e88, 0xbfef68f1d8e42e88, 0x3fecee874ffb88a6, 0xbfecee874ffb88a6, 0x3fe8a30aeed88f52, 0xbfe8a30aeed88f52, 0x3fe2cb4f05c07803, 0xbfe2cb4f05c07803, 0x3fd78a8d20a8b19f, 0xbfd78a8d20a8b19f, 0x3fc007a5f8f630e4, 0xbfc007a5f8f630e4],
    &[0x3fef7e6d762965ca, 0xbfef7e6d762965ca, 0x3fed5cf75170ca75, 0xbfed5cf75170ca75, 0x3fe9a687189c783f, 0xbfe9a687189c783f, 0x3fe48e2033b01c61, 0xbfe48e2033b01c61, 0x3fdcb41af08c7476, 0xbfdcb41af08c7476, 0x3fcd7fa8790f2c49, 0xbfcd7fa8790f2c49, 0x0000000000000000],
];

#[rustfmt::skip]
const WEIGHTS_BITS: [&[u64]; MAXN] = [
    &[0x4000000000000000],
    &[0x3ff0000000000000, 0x3ff0000000000000],
    &[0x3fe1c71c71c71c71, 0x3fe1c71c71c71c71, 0x3fec71c71c71c71e],
    &[0x3fe4de5f840c24ca, 0x3fe4de5f840c24ca, 0x3fd64340f7e7b66c, 0x3fd64340f7e7b66c],
    &[0x3fdea1da25ae415d, 0x3fdea1da25ae415d, 0x3fce539ec36e0387, 0x3fce539ec36e0387, 0x3fe23456789abcde],
    &[0x3fddf24d499545ee, 0x3fddf24d499545ee, 0x3fc5edf601e2dbf0, 0x3fc5edf601e2dbf0, 0x3fd716b7b5794c1a, 0x3fd716b7b5794c1a],
    &[0x3fd86fe74ee32b3d, 0x3fd86fe74ee32b3d, 0x3fd1e6b1713d8649, 0x3fd1e6b1713d8649, 0x3fc092f69f826d4d, 0x3fc092f69f826d4d, 0x3fdabfd7e03c2fa4],
    &[0x3fb9ea1d04ca0383, 0x3fb9ea1d04ca0383, 0x3fcc76fb531d2bae, 0x3fcc76fb531d2bae, 0x3fd413c50a25560a, 0x3fd413c50a25560a, 0x3fd736360b19933e, 0x3fd736360b19933e],
    &[0x3fb4ce65f803ef16, 0x3fb4ce65f803ef16, 0x3fc71f7a9b222bee, 0x3fc71f7a9b222bee, 0x3fd0add87c8274f3, 0x3fd0add87c8274f3, 0x3fd3fd7e9838d519, 0x3fd3fd7e9838d519, 0x3fd522a43f654872],
    &[0x3fb1115f8b62dc0f, 0x3fb1115f8b62dc0f, 0x3fc32138c878f039, 0x3fc32138c878f039, 0x3fcc0b059d00bc03, 0x3fcc0b059d00bc03, 0x3fd13baa7a559bfa, 0x3fd13baa7a559bfa, 0x3fd2e9de7014d6e2, 0x3fd2e9de7014d6e2],
    &[0x3fac8097265bbb80, 0x3fac8097265bbb80, 0x3fc013047def8696, 0x3fc013047def8696, 0x3fc7d85b8dbff17f, 0x3fc7d85b8dbff17f, 0x3fcdd94b1446e10d, 0x3fcdd94b1446e10d, 0x3fd0d1ca26fa5963, 0x3fd0d1ca26fa5963, 0x3fd1779ac87e052e],
    &[0x3fa8275d9dea6c39, 0x3fa8275d9dea6c39, 0x3fbb60602bce60e3, 0x3fbb60602bce60e3, 0x3fc47d7258f22da4, 0x3fc47d7258f22da4, 0x3fca0163e6b1abb3, 0x3fca0163e6b1abb3, 0x3fcde3155c256ad3, 0x3fcde3155c256ad3, 0x3fcfe40ce6d4f057, 0x3fcfe40ce6d4f057],
    &[0x3fa4ba51c8f4db2b, 0x3fa4ba51c8f4db2b, 0x3fb795464d0fb648, 0x3fb795464d0fb648, 0x3fc1c69b70565c9b, 0x3fc1c69b70565c9b, 0x3fc6cd7ccca4a75f, 0x3fc6cd7ccca4a75f, 0x3fca99b75be0a206, 0x3fca99b75be0a206, 0x3fccf6d8e56e9988, 0x3fccf6d8e56e9988, 0x3fcdc43fd1e15d0a],
];

fn table(bits: &[u64], buf: &mut [f64; MAXN]) -> usize {
    for (b, x) in bits.iter().zip(buf.iter_mut()) {
        *x = f64::from_bits(*b);
    }
    bits.len()
}

struct Tables {
    points: [[f64; MAXN]; MAXN],
    weights: [[f64; MAXN]; MAXN],
}

fn tables() -> &'static Tables {
    static T: std::sync::OnceLock<Tables> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        let mut t = Tables {
            points: [[0.0; MAXN]; MAXN],
            weights: [[0.0; MAXN]; MAXN],
        };
        for n in 0..MAXN {
            table(POINTS_BITS[n], &mut t.points[n]);
            table(WEIGHTS_BITS[n], &mut t.weights[n]);
        }
        t
    })
}

/// `Points(1:n,n)`; `n` must be in `1..=MAXN`.
pub(crate) fn points_1d(n: usize) -> &'static [f64] {
    &tables().points[n - 1][..n]
}

/// `Weights(1:n,n)`; `n` must be in `1..=MAXN`.
pub(crate) fn weights_1d(n: usize) -> &'static [f64] {
    &tables().weights[n - 1][..n]
}
