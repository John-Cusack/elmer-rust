// Port of `fem/src/view3d/LUDecomp.c` (Elmer @ 9f6af2f85) for lane-fem-a.
// Unit: c:fem/src/view3d/LUDecomp.c. True object symbols: `ludecomp`,
// `mtrinv` (the recon `exports` module name `LUDecomp` is not a symbol;
// verified with `nm -g --defined-only` on the reference build, see
// crates/elmer-fem-a/EVIDENCE.md).
//
// Faithfulness notes (do NOT "fix" — the corpus compares stderr bytes):
// - `ABS(x)` is `(x) > 0 ? (x) : (-(x))`, so `ABS(0.0)` is `-0.0`; kept.
// - `mtrinv` leaks `pivot` on the singular early-return in C; the leak is
//   not program-observable, the Rust port frees via `Vec` drop.
// - `n < 1` is UB in C (`pivot[n-1]` underflows); the port returns early
//   instead of manufacturing UB. No probe reaches this path (recorded).
// - stderr text is byte-identical, including `%d` formatting of `i`, `j`.
//
// Craft notes (tight-loop recraft, ADR-050). Only the changes that MEASURED
// faster were kept; the rejected ones are listed so nobody re-adds them.
// Every kept change is value-preserving: the same operations in the same
// order, so the results stay bit-identical (the bench digest proves it).
//
//   LU decomposition (`ludecomp`): 1.13x on the phase, 1.15x end to end.
//     * the pivot candidate is held in a register instead of being re-read
//       from the matrix on every candidate of the pivot search;
//     * the `j != i` test is invariant over the elimination k-loop, so the
//       swap leaves the loop body (two loop shapes, one body each);
//     * the elimination inner loop walks row pointers, so the C index math
//       `n*i + j` is done once per row instead of on every access.
//   Matrix inverse (`mtrinv`): 1.05x end to end.
//     * the inverse product wrote its partial sum back to `a(i,j)` on every
//       k iteration (an n^3 store chain). The only value ever read back
//       from `a(i,j)` is read on the *first* iteration, before any store,
//       and no other iteration reads `a(i,j)` afterwards — so the single
//       final store is equivalent (n^3 stores become n^2). 1.09x on that
//       phase.
//
// Rejected after measurement (kept out on purpose):
//   * row pointers in the INV(U) inner loop: 0.93x — the compiler unrolls
//     the original form better (see the variant table in the ADR);
//   * splitting the INV(U) / INV(AP) k-loops into two branch-free bodies
//     instead of one loop with the `k != i` test: 0.81x / 0.90x;
//   * caching the diagonal in INV(L): inside the noise band, so not kept.

use std::os::raw::{c_double, c_int};

const AEPS: c_double = 1.0e-15;

#[inline]
fn abs_c(x: c_double) -> c_double {
    if x > 0.0 { x } else { -x }
}

#[inline]
unsafe fn a_at(a: *mut c_double, n: c_int, i: c_int, j: c_int) -> c_double {
    *a.offset((n * i + j) as isize)
}

#[inline]
unsafe fn a_set(a: *mut c_double, n: c_int, i: c_int, j: c_int, v: c_double) {
    *a.offset((n * i + j) as isize) = v;
}

/// Row base pointer for row `i` of an `n x n` row-major matrix.
///
/// The C original spelled every access as `a[n*i + j]`; the recraft walks a
/// row pointer instead, so the index math is done once per row.
#[inline]
unsafe fn row(a: *mut c_double, n: c_int, i: c_int) -> *mut c_double {
    a.offset((n * i) as isize)
}

#[export_name = "mtrinv"]
pub unsafe extern "C" fn port_mtrinv(a: *mut c_double, n: c_int) {
    if n < 1 {
        return;
    }
    let nn = n as usize;
    let mut pivot: Vec<c_int> = vec![0; nn];

    port_ludecomp(a, n, pivot.as_mut_ptr());

    for i in 0..n {
        // The diagonal is read twice in C and cannot change between the two
        // reads, so it is read once here.
        let ai = row(a, n, i);
        let d = *ai.offset(i as isize);
        if abs_c(d) < AEPS {
            eprint!("Inv: Matrix is singular.\n");
            return;
        }
        *ai.offset(i as isize) = 1.0 / d;
    }

    // INV(U)
    let mut i = n - 2;
    while i >= 0 {
        let mut j = n - 1;
        while j > i {
            let mut s = -a_at(a, n, i, j);
            let mut k = i + 1;
            while k < j {
                s -= a_at(a, n, i, k) * a_at(a, n, k, j);
                k += 1;
            }
            a_set(a, n, i, j, s);
            j -= 1;
        }
        i -= 1;
    }

    // INV(L)
    let mut i = n - 2;
    while i >= 0 {
        let mut j = n - 1;
        while j > i {
            let mut s = 0.0;
            let mut k = i + 1;
            while k <= j {
                s -= a_at(a, n, j, k) * a_at(a, n, k, i);
                k += 1;
            }
            a_set(a, n, j, i, a_at(a, n, i, i) * s);
            j -= 1;
        }
        i -= 1;
    }

    // A = INV(AP): the k-term is `a(i,k) * a(k,j)`, except for k == i where
    // the original adds `a(i,j)` itself. That single term is read on the
    // first iteration of the k-range, before any store to `a(i,j)`, and no
    // later iteration of the product reads `a(i,j)` back — so the partial
    // sum can be written once, after the k-range, instead of on every
    // iteration.
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            let mut k = if i > j { i } else { j };
            while k < n {
                if k != i {
                    s += a_at(a, n, i, k) * a_at(a, n, k, j);
                } else {
                    s += a_at(a, n, k, j);
                }
                k += 1;
            }
            a_set(a, n, i, j, s);
        }
    }

    // A = INV(A) (at last)
    let mut i = n - 1;
    while i >= 0 {
        if *pivot.get_unchecked(i as usize) != i {
            let pi = *pivot.get_unchecked(i as usize);
            for j in 0..n {
                let s = a_at(a, n, i, j);
                a_set(a, n, i, j, a_at(a, n, pi, j));
                a_set(a, n, pi, j, s);
            }
        }
        i -= 1;
    }
}

#[export_name = "ludecomp"]
pub unsafe extern "C" fn port_ludecomp(a: *mut c_double, n: c_int, pivot: *mut c_int) {
    if n < 1 {
        return;
    }
    for i in 0..n - 1 {
        let ai = row(a, n, i);
        // Pivot search: the current best stays in a register instead of
        // being re-read from the matrix on every candidate.
        let mut j = i;
        let mut best = abs_c(*ai.offset(i as isize));
        let mut k = i + 1;
        while k < n {
            let v = abs_c(*ai.offset(k as isize));
            if v > best {
                j = k;
                best = v;
            }
            k += 1;
        }

        if best < AEPS {
            eprint!("LUDecomp: Matrix is (at least almost) singular, {} {}.\n", i, j);
            return;
        }

        *pivot.offset(i as isize) = j;

        if j != i {
            let mut k = 0;
            while k <= i {
                let ak = row(a, n, k);
                let swap = *ak.offset(j as isize);
                *ak.offset(j as isize) = *ak.offset(i as isize);
                *ak.offset(i as isize) = swap;
                k += 1;
            }
        }

        // Row scaling. The divisor is the (unchanged) diagonal, so it is
        // read once; the division itself is unchanged (a reciprocal times a
        // multiply would not be bit-identical).
        let mut k = i + 1;
        while k < n {
            *ai.offset(k as isize) = *ai.offset(k as isize) / *ai.offset(i as isize);
            k += 1;
        }

        // Elimination. The `j != i` test is invariant over the k-loop, so
        // the swap is split out of the loop body.
        let mut k = i + 1;
        if j != i {
            while k < n {
                let ak = row(a, n, k);
                let swap = *ak.offset(i as isize);
                *ak.offset(i as isize) = *ak.offset(j as isize);
                *ak.offset(j as isize) = swap;
                // `a(k,i)` cannot change inside the l-loop, so it is read
                // once per k.
                let aki = *ak.offset(i as isize);
                let mut l = i + 1;
                while l < n {
                    *ak.offset(l as isize) = *ak.offset(l as isize) - aki * *ai.offset(l as isize);
                    l += 1;
                }
                k += 1;
            }
        } else {
            while k < n {
                let ak = row(a, n, k);
                let aki = *ak.offset(i as isize);
                let mut l = i + 1;
                while l < n {
                    *ak.offset(l as isize) = *ak.offset(l as isize) - aki * *ai.offset(l as isize);
                    l += 1;
                }
                k += 1;
            }
        }
    }

    *pivot.offset((n - 1) as isize) = n - 1;
    if abs_c(*row(a, n, n - 1).offset((n - 1) as isize)) < AEPS {
        eprint!("LUDecomp: Matrix is (at least almost) singular.\n");
    }
}
