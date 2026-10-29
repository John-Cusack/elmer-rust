//! Final ordering of the eigenpairs: `EigenSystemSorting` (EigenSolve.F90:74)
//! over `SortR` (GeneralUtils.F90:784).

use num_complex::Complex64;

use crate::config::Which;

/// `SortR`: heap sort of `b` into decreasing order, permuting `a` alongside.
/// Transliterated rather than replaced by `sort_by`, because the heap sort is
/// not stable and Elmer's order of equal measures is part of its output.
pub fn sort_r(a: &mut [usize], b: &mut [f64]) {
    let n = a.len().min(b.len());
    if n <= 1 {
        return;
    }
    // 1-based indices, as in the Fortran.
    let mut l = n / 2 + 1;
    let mut ir = n;
    loop {
        let (ra, rb);
        if l > 1 {
            l -= 1;
            ra = a[l - 1];
            rb = b[l - 1];
        } else {
            ra = a[ir - 1];
            rb = b[ir - 1];
            a[ir - 1] = a[0];
            b[ir - 1] = b[0];
            ir -= 1;
            if ir == 1 {
                a[0] = ra;
                b[0] = rb;
                return;
            }
        }
        let mut i = l;
        let mut j = l + l;
        while j <= ir {
            if j < ir && b[j - 1] > b[j] {
                j += 1;
            }
            if rb > b[j - 1] {
                a[i - 1] = a[j - 1];
                b[i - 1] = b[j - 1];
                i = j;
                j += i;
            } else {
                j = ir + 1;
            }
            a[i - 1] = ra;
            b[i - 1] = rb;
        }
    }
}

/// `EigenSystemSorting`: the permutation (0-based) that puts the eigenvalue
/// selected by `sorting` first. Elmer applies it to the values and, through it,
/// to the vectors.
pub fn eigen_system_sorting(sorting: Which, values: &[Complex64]) -> Vec<usize> {
    let n = values.len();
    let mut perm: Vec<usize> = (0..n).collect();
    // Elmer encodes the request as an inverted ARPACK `WHICH`: "smallest
    // magnitude" becomes "LM" and the measure is negated for an 'L', so the
    // decreasing heap sort lists the smallest first.
    let (largest_first_negated, measure): (bool, fn(&Complex64) -> f64) = match sorting {
        Which::SmallestMagnitude => (true, |z| z.norm()),
        Which::LargestMagnitude => (false, |z| z.norm()),
        Which::SmallestReal => (true, |z| z.re),
        Which::LargestReal => (false, |z| z.re),
        Which::SmallestImag => (true, |z| z.im),
        Which::LargestImag => (false, |z| z.im),
    };
    let mut meas: Vec<f64> = values
        .iter()
        .map(|z| {
            let m = measure(z);
            if largest_first_negated {
                -m
            } else {
                m
            }
        })
        .collect();
    sort_r(&mut perm, &mut meas);
    perm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sort_r_is_decreasing() {
        let mut b = vec![3.0, -1.0, 7.5, 0.0, 2.0, 7.5, -4.0];
        let mut a: Vec<usize> = (0..b.len()).collect();
        let orig = b.clone();
        sort_r(&mut a, &mut b);
        for w in b.windows(2) {
            assert!(w[0] >= w[1]);
        }
        for (k, &i) in a.iter().enumerate() {
            assert_eq!(orig[i], b[k]);
        }
    }

    #[test]
    fn sort_r_tie_order_is_the_heap_order() {
        // Hand-run of SortR on four equal keys: the heap sort moves the tied
        // entries, so a stable sort would give [0, 1, 2, 3].
        let mut b = vec![1.0; 4];
        let mut a = vec![0, 1, 2, 3];
        sort_r(&mut a, &mut b);
        assert_eq!(a, vec![1, 2, 3, 0]);
    }

    #[test]
    fn default_sorting_is_ascending_magnitude() {
        let v: Vec<Complex64> = [5.0, -1.0, 3.0, 0.5]
            .iter()
            .map(|&x| Complex64::new(x, 0.0))
            .collect();
        let p = eigen_system_sorting(Which::SmallestMagnitude, &v);
        assert_eq!(p, vec![3, 1, 2, 0]);
        let p = eigen_system_sorting(Which::LargestMagnitude, &v);
        assert_eq!(p, vec![0, 2, 1, 3]);
        let p = eigen_system_sorting(Which::SmallestReal, &v);
        assert_eq!(p, vec![1, 3, 2, 0]);
        let p = eigen_system_sorting(
            Which::LargestImag,
            &[
                Complex64::new(0.0, -2.0),
                Complex64::new(0.0, 1.0),
                Complex64::new(0.0, 3.0),
            ],
        );
        assert_eq!(p, vec![2, 1, 0]);
    }
}
