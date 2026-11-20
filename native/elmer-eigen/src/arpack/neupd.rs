//! `DNEUPD` (ARPACK 2.4): the Ritz values and Ritz vectors of the final
//! Arnoldi factorisation [`super::naupd::dnaupd`] left in `v` and `workl`,
//! back-transformed for the spectral transformation of `mode`.
//!
//! Only what Elmer's `ArpackEigenSolve` asks for is ported: `rvec =
//! .true.` with `howmny = 'A'`, and the Ritz vectors written over `v`
//! (Elmer passes `V` as both `v` and `z`). A complex pair of Ritz values
//! has its vector in two columns, the real part then the imaginary part.
//! `workl` is the flat array of the Fortran, addressed with its 1-based
//! offsets.

use super::blas::{dcopy, dger, dgemv_t, dnrm2, dscal};
use super::lapack::{dgeqr2, dlapy2, dlaset, dorm2r, EPS};
use super::lapack_ns::{dlahqr, dtrevc_right, dtrmm_run, dtrsen_nv};
use super::naupd::{dngets, NWhich, Naupd};
use super::par::Par;

#[derive(Copy, Clone, PartialEq, Eq)]
enum Kind {
    Regular,
    ShiftI,
    RealPt,
    ImagPt,
}

/// What [`dneupd`] returns.
#[derive(Clone, Debug, Default)]
pub struct Neupd {
    /// `DR`, `DI`: real and imaginary parts of the Ritz values
    /// (`nev + 1` entries; the first `nconv` are set).
    pub dr: Vec<f64>,
    pub di: Vec<f64>,
    /// The number of Ritz values computed (`IPARAM(5)`, or `M` of
    /// `DTRSEN` after a reordering).
    pub nconv: usize,
    pub info: i32,
}

/// `DNEUPD(.true., 'A', select, dr, di, v, ldv, sigmar, sigmai, workev,
/// bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd,
/// workl, lworkl, info)`. The Ritz vectors are in the first `nconv`
/// columns of `v`.
#[allow(clippy::too_many_arguments)]
pub fn dneupd(
    sigmar: f64,
    sigmai: f64,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Naupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> Neupd {
    neupd(None, sigmar, sigmai, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, au, workd, workl)
}

/// PARPACK's `PDNEUPD(comm, ...)` after [`super::naupd::pdnaupd`]: the
/// wanted Ritz values picked as PARPACK does, by sorting them with
/// `_ngets` and taking the first `IPARAM(5)` whose estimates pass
/// (`info = -15` if fewer do), in place of ARPACK 2.4's threshold test.
/// The rest is rank-local (`n` rows of the vectors).
#[allow(clippy::too_many_arguments)]
pub fn pdneupd(
    comm: &dyn elmer_core::Communicator,
    sigmar: f64,
    sigmai: f64,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Naupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> Neupd {
    neupd(Some(comm), sigmar, sigmai, bmat_g, n, which, nev, tol, resid, ncv, v, ldv, au, workd, workl)
}

#[allow(clippy::too_many_arguments)]
fn neupd(
    par: Par<'_>,
    sigmar: f64,
    sigmai: f64,
    bmat_g: bool,
    n: usize,
    which: NWhich,
    nev: usize,
    tol: f64,
    resid: &[f64],
    ncv: usize,
    v: &mut [f64],
    ldv: usize,
    au: &Naupd,
    workd: &mut [f64],
    workl: &mut [f64],
) -> Neupd {
    let mut out = Neupd { dr: vec![0.0; nev + 1], di: vec![0.0; nev + 1], ..Default::default() };
    let mode = au.iparam[6];
    let nconv_in = au.iparam[4];
    let eps23 = EPS.powf(2.0 / 3.0);
    let mut ierr = 0;
    if nconv_in <= 0 {
        ierr = -14;
    } else if n == 0 {
        ierr = -1;
    } else if nev == 0 {
        ierr = -2;
    } else if ncv <= nev + 1 || (ncv > n && par.is_none()) {
        ierr = -3;
    } else if workl.len() < 3 * ncv * ncv + 6 * ncv {
        ierr = -7;
    }
    let kind = match mode {
        1 | 2 => Kind::Regular,
        3 if sigmai == 0.0 => Kind::ShiftI,
        3 => Kind::RealPt,
        4 => Kind::ImagPt,
        _ => {
            ierr = -10;
            Kind::Regular
        }
    };
    if mode == 1 && bmat_g {
        ierr = -11;
    }
    if ierr != 0 {
        out.info = ierr;
        return out;
    }
    let mut nconv = nconv_in as usize;
    let mut workev = vec![0.0f64; 3 * ncv];
    let mut select = vec![false; ncv];
    // Pointer into WORKL for address of H, RITZ, BOUNDS, Q etc. (1-based,
    // as the Fortran).
    let ih = au.ipntr[4];
    let ritzr = au.ipntr[5];
    let ritzi = au.ipntr[6];
    let bounds = au.ipntr[7];
    let (ldh, ldq) = (ncv, ncv);
    let iheigr = bounds + ldh;
    let iheigi = iheigr + ldh;
    let ihbds = iheigi + ldh;
    let iuptri = ihbds + ldh;
    let invsub = iuptri + ldh * ncv;
    let irr = au.ipntr[13] + ncv * ncv;
    let iri = irr + ncv;
    let ibd = iri + ncv;
    macro_rules! w {
        ($k:expr) => {
            workl[$k - 1]
        };
    }
    // RNORM is B-norm of the RESID(1:N).
    let rnorm = w!(ih + 2);
    w!(ih + 2) = 0.0;

    // Get converged Ritz value on the boundary. Note: converged Ritz
    // values have been placed in the first NCONV locations in workl(ritzr)
    // and workl(ritzi). They have been sorted (in _naup2) according to the
    // WHICH selection criterion.
    let thres = match which {
        NWhich::LM | NWhich::SM => dlapy2(w!(ritzr), w!(ritzi)),
        NWhich::LR | NWhich::SR => w!(ritzr),
        NWhich::LI | NWhich::SI => w!(ritzi).abs(),
    };
    // Check to see if all converged Ritz values appear at the top of the
    // upper quasi-triangular matrix computed by _neigh in _naup2. This is
    // done in the following way:
    //
    // 1) For each Ritz value obtained from _neigh, compare it with the
    //    threshold Ritz value computed above to determine whether it is a
    //    wanted one.
    //
    // 2) If it is wanted, then check the corresponding Ritz estimate to
    //    see if it has converged. If it has, set correponding entry in the
    //    logical array SELECT to .TRUE..
    //
    // If SELECT(j) = .TRUE. and j > NCONV, then there is a converged Ritz
    // value that does not appear at the top of the upper quasi-triangular
    // matrix computed by _neigh in _naup2. Reordering is needed.
    let mut reord = false;
    if par.is_some() {
        // pdneupd: sort the Ritz values (and their indices, in the bounds
        // slot) with _ngets; take the first nconv converged from the wanted
        // end. _ngets may move a conjugate pair's split, nev with it.
        for j in 0..ncv {
            w!(bounds + j) = (j + 1) as f64;
        }
        let (mut kev, mut np) = (nev, ncv - nev);
        {
            let (lo, hi) = workl.split_at_mut(irr - 1);
            let (rr, rest) = hi.split_at_mut(ncv);
            if let Err(e) = dngets(0, which, &mut kev, &mut np, rr, &mut rest[..ncv], &mut lo[bounds - 1..bounds - 1 + ncv]) {
                debug_assert!(false, "{e}");
                out.info = -15;
                return out;
            }
        }
        let mut numcnv = 0;
        for j in 1..=ncv {
            let temp1 = eps23.max(dlapy2(w!(irr + ncv - j), w!(iri + ncv - j)));
            let jj = w!(bounds + ncv - j) as usize;
            if numcnv < nconv && w!(ibd + jj - 1) <= tol * temp1 {
                select[jj - 1] = true;
                numcnv += 1;
                if jj > kev {
                    reord = true;
                }
            }
        }
        if numcnv != nconv {
            out.info = -15;
            return out;
        }
    }
    for j in (0..ncv).filter(|_| par.is_none()) {
        select[j] = false;
        let (xr, xi, bd) = (w!(irr + j), w!(iri + j), w!(ibd + j));
        let wanted = match which {
            NWhich::LM => dlapy2(xr, xi) >= thres,
            NWhich::SM => dlapy2(xr, xi) <= thres,
            NWhich::LR => xr >= thres,
            NWhich::SR => xr <= thres,
            NWhich::LI => xi.abs() >= thres,
            NWhich::SI => xi.abs() <= thres,
        };
        if wanted {
            let temp1 = eps23.max(dlapy2(xr, xi));
            if bd <= tol * temp1 {
                select[j] = true;
            }
        }
        if j + 1 > nconv {
            reord = select[j] || reord;
        }
    }
    // Call LAPACK routine dlahqr to compute the real Schur form of the
    // upper Hessenberg matrix returned by DNAUPD. Make a copy of the upper
    // Hessenberg matrix. Initialize the Schur vector matrix Q to the
    // identity.
    workl.copy_within(ih - 1..ih - 1 + ldh * ncv, iuptri - 1);
    dlaset(ncv, ncv, 0.0, 1.0, &mut workl[invsub - 1..], ldq);
    let ierr = {
        let (lo, hi) = workl.split_at_mut(invsub - 1);
        let (heig, up) = lo[iheigr - 1..].split_at_mut(iuptri - iheigr);
        let (hr, hi_) = heig.split_at_mut(ldh);
        dlahqr(true, true, ncv, 1, ncv, up, ldh, hr, &mut hi_[..ldh], 1, ncv, hi, ldq)
    };
    // `dcopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)`.
    let last_row = |workl: &mut [f64]| {
        for k in 0..ncv {
            workl[ihbds - 1 + k] = workl[invsub + ncv - 2 + k * ldq];
        }
    };
    last_row(workl);
    if ierr != 0 {
        out.info = -8;
        return out;
    }
    if reord {
        // Reorder the computed upper quasi-triangular matrix.
        let (m, info) = {
            let (lo, hi) = workl.split_at_mut(invsub - 1);
            let (heig, up) = lo[iheigr - 1..].split_at_mut(iuptri - iheigr);
            let (hr, hi_) = heig.split_at_mut(ldh);
            dtrsen_nv(&select, ncv, up, ldh, hi, ldq, hr, &mut hi_[..ldh])
        };
        // `DTRSEN(..., nconv, ...)`: M is returned in NCONV.
        nconv = m;
        if info == 1 {
            out.info = 1;
            out.nconv = nconv;
            return out;
        }
    }
    // Copy the last row of the Schur basis matrix to workl(ihbds). This
    // vector will be used to compute the Ritz estimates of converged Ritz
    // values.
    last_row(workl);
    // Place the computed eigenvalues of H into DR and DI if a spectral
    // transformation was not used.
    if kind == Kind::Regular {
        out.dr[..nconv].copy_from_slice(&workl[iheigr - 1..iheigr - 1 + nconv]);
        out.di[..nconv].copy_from_slice(&workl[iheigi - 1..iheigi - 1 + nconv]);
    }
    // Compute the QR factorization of the matrix representing the wanted
    // invariant subspace located in the first NCONV columns of
    // workl(invsub,ldq).
    {
        let (a, b) = workev.split_at_mut(ncv);
        dgeqr2(ncv, nconv, &mut workl[invsub - 1..], ldq, a, b);
    }
    // * Postmultiply V by Q using dorm2r ('Right', 'Notranspose').
    // * Copy the first NCONV columns of VQ into Z (Z is V here).
    dorm2r(false, false, n, ncv, nconv, &mut workl[invsub - 1..], ldq, &workev, v, ldv, &mut workd[n..]);
    // Make sure that the diagonal elements of the upper quasi-triangular
    // matrix are positive.
    for j in 1..=nconv {
        if w!(invsub + (j - 1) * ldq + j - 1) < 0.0 {
            // `dscal(nconv, -one, workl(iuptri+j-1), ldq)`: row j.
            for k in 0..nconv {
                let p = iuptri + j - 2 + k * ldq;
                workl[p] *= -1.0;
            }
            dscal(nconv, -1.0, &mut workl[iuptri - 1 + (j - 1) * ldq..]);
        }
    }
    // HOWMNY = 'A': compute NCONV Ritz vectors.
    for (j, s) in select.iter_mut().enumerate() {
        *s = j < nconv;
    }
    {
        let (lo, hi) = workl.split_at_mut(invsub - 1);
        let t = &lo[iuptri - 1..iuptri - 1 + ldh * ncv];
        dtrevc_right(Some(&mut select[..]), ncv, t, ldq, hi, ldq, &mut workev);
    }
    // Scale the returning eigenvectors so that their Euclidean norms are
    // all one. A complex conjugate pair shares one scale.
    let mut iconj = false;
    for j in 1..=nconv {
        let c = invsub - 1 + (j - 1) * ldq;
        if w!(iheigi + j - 1) == 0.0 {
            let temp = dnrm2(ncv, &workl[c..]);
            dscal(ncv, 1.0 / temp, &mut workl[c..]);
        } else if !iconj {
            let c2 = invsub - 1 + j * ldq;
            let temp = dlapy2(dnrm2(ncv, &workl[c..]), dnrm2(ncv, &workl[c2..]));
            dscal(ncv, 1.0 / temp, &mut workl[c..]);
            dscal(ncv, 1.0 / temp, &mut workl[c2..]);
            iconj = true;
        } else {
            iconj = false;
        }
    }
    {
        let (lo, hi) = workl.split_at_mut(invsub - 1);
        dgemv_t(ncv, nconv, 1.0, hi, ldq, &lo[ihbds - 1..], 0.0, &mut workev);
    }
    let mut iconj = false;
    for j in 1..=nconv {
        if w!(iheigi + j - 1) != 0.0 {
            // Complex conjugate pair case. Note that since the real and
            // imaginary part of the eigenvector are stored in consecutive
            // columns, we further normalize by the square root of two.
            if !iconj {
                workev[j - 1] = dlapy2(workev[j - 1], workev[j]);
                workev[j] = workev[j - 1];
                iconj = true;
            } else {
                iconj = false;
            }
        }
    }
    dcopy(nconv, &workev, &mut workl[ihbds - 1..]);
    // Compute the QR factorization of the eigenvector matrix associated
    // with leading portion of T in the first NCONV columns of
    // workl(invsub,ldq).
    {
        let (a, b) = workev.split_at_mut(ncv);
        dgeqr2(ncv, nconv, &mut workl[invsub - 1..], ldq, a, b);
    }
    // * Postmultiply Z by Q.
    // * Postmultiply Z by R.
    // The N by NCONV matrix Z is now contains the Ritz vectors associated
    // with the Ritz values in workl(iheigr) and workl(iheigi).
    dorm2r(false, false, n, ncv, nconv, &mut workl[invsub - 1..], ldq, &workev, v, ldv, &mut workd[n..]);
    dtrmm_run(n, nconv, 1.0, &workl[invsub - 1..], ldq, v, ldv);

    // Transform the Ritz values and possibly vectors and corresponding
    // error bounds of OP to those of A*x = lambda*B*x.
    if kind == Kind::Regular {
        dscal(ncv, rnorm, &mut workl[ihbds - 1..]);
    } else {
        // A spectral transformation was used.
        // * Determine the Ritz estimates of the Ritz values in the
        //   original system.
        if kind == Kind::ShiftI {
            dscal(ncv, rnorm, &mut workl[ihbds - 1..]);
            for k in 1..=ncv {
                let temp = dlapy2(w!(iheigr + k - 1), w!(iheigi + k - 1));
                w!(ihbds + k - 1) = w!(ihbds + k - 1).abs() / temp / temp;
            }
        }
        // * Transform the Ritz values back to the original system. For
        //   TYPE = 'SHIFTI' the transformation is
        //   lambda = 1/theta + sigma; REALPT / IMAGPT leave them.
        if kind == Kind::ShiftI {
            for k in 1..=ncv {
                let temp = dlapy2(w!(iheigr + k - 1), w!(iheigi + k - 1));
                w!(iheigr + k - 1) = w!(iheigr + k - 1) / temp / temp + sigmar;
                w!(iheigi + k - 1) = -w!(iheigi + k - 1) / temp / temp + sigmai;
            }
        }
        if matches!(kind, Kind::ShiftI | Kind::RealPt | Kind::ImagPt) {
            out.dr[..nconv].copy_from_slice(&workl[iheigr - 1..iheigr - 1 + nconv]);
            out.di[..nconv].copy_from_slice(&workl[iheigi - 1..iheigi - 1 + nconv]);
        }
    }
    // Eigenvector purification step. Formally perform one of inverse
    // subspace iteration. Only used for MODE = 3. See reference 3.
    if kind == Kind::ShiftI {
        // Purify the computed Ritz vectors by adding a little bit of the
        // residual vector:
        //                      T
        //          resid(:)*( e    s ) / theta
        //                      NCV
        // where H s = s theta. Remember that when theta has nonzero
        // imaginary part, the corresponding Ritz vector is stored across
        // two columns of Z.
        let mut iconj = false;
        for j in 1..=nconv {
            let c = invsub + (j - 1) * ldq + ncv - 1;
            if w!(iheigi + j - 1) == 0.0 {
                workev[j - 1] = w!(c) / w!(iheigr + j - 1);
            } else if !iconj {
                let c2 = invsub + j * ldq + ncv - 1;
                let temp = dlapy2(w!(iheigr + j - 1), w!(iheigi + j - 1));
                workev[j - 1] = (w!(c) * w!(iheigr + j - 1) + w!(c2) * w!(iheigi + j - 1)) / temp / temp;
                workev[j] = (w!(c2) * w!(iheigr + j - 1) - w!(c) * w!(iheigi + j - 1)) / temp / temp;
                iconj = true;
            } else {
                iconj = false;
            }
        }
        // Perform a rank one update to Z and purify all the Ritz vectors
        // together.
        dger(n, nconv, 1.0, resid, &workev, v, ldv);
    }
    out.nconv = nconv;
    out
}
