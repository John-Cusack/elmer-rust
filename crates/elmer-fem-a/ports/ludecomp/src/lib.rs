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

#[export_name = "mtrinv"]
pub unsafe extern "C" fn port_mtrinv(a: *mut c_double, n: c_int) {
    if n < 1 {
        return;
    }
    let nn = n as usize;
    let mut pivot: Vec<c_int> = vec![0; nn];

    port_ludecomp(a, n, pivot.as_mut_ptr());

    for i in 0..n {
        if abs_c(a_at(a, n, i, i)) < AEPS {
            eprint!("Inv: Matrix is singular.\n");
            return;
        }
        a_set(a, n, i, i, 1.0 / a_at(a, n, i, i));
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

    // A = INV(AP): note the store `A(i,j) = s` sits INSIDE the k-loop in
    // the original (kept verbatim — final iteration wins, same value).
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
                a_set(a, n, i, j, s);
                k += 1;
            }
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
        let mut j = i;
        let mut k = i + 1;
        while k < n {
            if abs_c(a_at(a, n, i, k)) > abs_c(a_at(a, n, i, j)) {
                j = k;
            }
            k += 1;
        }

        if abs_c(a_at(a, n, i, j)) < AEPS {
            eprint!("LUDecomp: Matrix is (at least almost) singular, {} {}.\n", i, j);
            return;
        }

        *pivot.offset(i as isize) = j;

        if j != i {
            let mut k = 0;
            while k <= i {
                let swap = a_at(a, n, k, j);
                a_set(a, n, k, j, a_at(a, n, k, i));
                a_set(a, n, k, i, swap);
                k += 1;
            }
        }

        let mut k = i + 1;
        while k < n {
            a_set(a, n, i, k, a_at(a, n, i, k) / a_at(a, n, i, i));
            k += 1;
        }

        let mut k = i + 1;
        while k < n {
            if j != i {
                let swap = a_at(a, n, k, i);
                a_set(a, n, k, i, a_at(a, n, k, j));
                a_set(a, n, k, j, swap);
            }
            let mut l = i + 1;
            while l < n {
                a_set(a, n, k, l, a_at(a, n, k, l) - a_at(a, n, k, i) * a_at(a, n, i, l));
                l += 1;
            }
            k += 1;
        }
    }

    *pivot.offset((n - 1) as isize) = n - 1;
    if abs_c(a_at(a, n, n - 1, n - 1)) < AEPS {
        eprint!("LUDecomp: Matrix is (at least almost) singular.\n");
    }
}
