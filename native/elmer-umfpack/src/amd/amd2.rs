//! `AMD_2`: the approximate minimum degree ordering on the quotient graph.
//!
//! AMD Version 1.1, Copyright (c) 2004 by Timothy A. Davis, Patrick R.
//! Amestoy, and Iain S. Duff.  Translated to Rust.

use super::*;

/// On output `last` holds the permutation and `next` its inverse.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn amd_2(
    n: Int,
    Pe: &mut [Int],
    Iw: &mut [Int],
    Len: &mut [Int],
    iwlen: Int,
    mut pfree: Int,
    Nv: &mut [Int],
    Next: &mut [Int],
    Last: &mut [Int],
    Head: &mut [Int],
    Elen: &mut [Int],
    Degree: &mut [Int],
    W: &mut [Int],
    control: Option<&[f64]>,
    info: Option<&mut [f64]>,
) {
    macro_rules! ix {
        ($a:ident [ $i:expr ]) => {
            $a[($i) as usize]
        };
    }
    let mut lnz: f64 = 0.0;
    let mut ndiv: f64 = 0.0;
    let mut nms_lu: f64 = 0.0;
    let mut nms_ldl: f64 = 0.0;
    let mut dmax: f64 = 1.0;
    let mut me: Int = EMPTY;
    let mut wflg: Int = 2;
    let mut mindeg: Int = 0;
    let mut ncmpa: Int = 0;
    let mut nel: Int = 0;
    let mut lemax: Int = 0;
    let wbig: Int = Int::MAX - n;

    let (alpha, aggressive) = match control {
        Some(c) => (c[AMD_DENSE], c[AMD_AGGRESSIVE] != 0.0),
        None => (AMD_DEFAULT_DENSE, AMD_DEFAULT_AGGRESSIVE != 0.0),
    };
    let mut dense: Int = if alpha < 0.0 {
        n
    } else {
        (alpha * (n as f64).sqrt()) as Int
    };
    dense = dense.max(16);
    dense = dense.min(n);

    for i in 0..n as usize {
        Last[i] = EMPTY;
        Head[i] = EMPTY;
        Next[i] = EMPTY;
        Nv[i] = 1;
        W[i] = 1;
        Elen[i] = 0;
        Degree[i] = Len[i];
    }

    let mut ndense: Int = 0;
    for i in 0..n {
        let iu = i as usize;
        let deg = Degree[iu];
        if deg == 0 {
            Elen[iu] = flip(1);
            nel += 1;
            Pe[iu] = EMPTY;
            W[iu] = 0;
        } else if deg > dense {
            ndense += 1;
            Nv[iu] = 0;
            Elen[iu] = EMPTY;
            nel += 1;
            Pe[iu] = EMPTY;
        } else {
            let inext = ix!(Head[deg]);
            if inext != EMPTY {
                ix!(Last[inext]) = i;
            }
            Next[iu] = inext;
            ix!(Head[deg]) = i;
        }
    }

    while nel < n {
        // get pivot of minimum degree
        let mut deg = mindeg;
        while deg < n {
            me = ix!(Head[deg]);
            if me != EMPTY {
                break;
            }
            deg += 1;
        }
        mindeg = deg;
        let meu = me as usize;
        let inext = Next[meu];
        if inext != EMPTY {
            ix!(Last[inext]) = EMPTY;
        }
        ix!(Head[deg]) = inext;

        let elenme = Elen[meu];
        let mut nvpiv = Nv[meu];
        nel += nvpiv;

        // construct new element
        Nv[meu] = -nvpiv;
        let mut degme: Int = 0;
        let mut pme1: Int;
        let pme2: Int;

        if elenme == 0 {
            pme1 = Pe[meu];
            let mut p2 = pme1 - 1;
            let mut p = pme1;
            while p <= pme1 + Len[meu] - 1 {
                let i = ix!(Iw[p]);
                let iu = i as usize;
                let nvi = Nv[iu];
                if nvi > 0 {
                    degme += nvi;
                    Nv[iu] = -nvi;
                    p2 += 1;
                    ix!(Iw[p2]) = i;
                    let ilast = Last[iu];
                    let inext = Next[iu];
                    if inext != EMPTY {
                        ix!(Last[inext]) = ilast;
                    }
                    if ilast != EMPTY {
                        ix!(Next[ilast]) = inext;
                    } else {
                        ix!(Head[Degree[iu]]) = inext;
                    }
                }
                p += 1;
            }
            pme2 = p2;
        } else {
            let mut p = Pe[meu];
            pme1 = pfree;
            let slenme = Len[meu] - elenme;
            let mut knt1: Int = 1;
            while knt1 <= elenme + 1 {
                let e: Int;
                let mut pj: Int;
                let ln: Int;
                if knt1 > elenme {
                    e = me;
                    pj = p;
                    ln = slenme;
                } else {
                    e = ix!(Iw[p]);
                    p += 1;
                    pj = ix!(Pe[e]);
                    ln = ix!(Len[e]);
                }
                let mut knt2: Int = 1;
                while knt2 <= ln {
                    let i = ix!(Iw[pj]);
                    pj += 1;
                    let iu = i as usize;
                    let nvi = Nv[iu];
                    if nvi > 0 {
                        if pfree >= iwlen {
                            // prepare for compressing Iw
                            Pe[meu] = p;
                            Len[meu] -= knt1;
                            if Len[meu] == 0 {
                                Pe[meu] = EMPTY;
                            }
                            ix!(Pe[e]) = pj;
                            ix!(Len[e]) = ln - knt2;
                            if ix!(Len[e]) == 0 {
                                ix!(Pe[e]) = EMPTY;
                            }
                            ncmpa += 1;
                            for j in 0..n {
                                let ju = j as usize;
                                let pn = Pe[ju];
                                if pn >= 0 {
                                    Pe[ju] = ix!(Iw[pn]);
                                    ix!(Iw[pn]) = flip(j);
                                }
                            }
                            let mut psrc: Int = 0;
                            let mut pdst: Int = 0;
                            let pend = pme1 - 1;
                            while psrc <= pend {
                                let j = flip(ix!(Iw[psrc]));
                                psrc += 1;
                                if j >= 0 {
                                    let ju = j as usize;
                                    ix!(Iw[pdst]) = Pe[ju];
                                    Pe[ju] = pdst;
                                    pdst += 1;
                                    let lenj = Len[ju];
                                    let mut knt3: Int = 0;
                                    while knt3 <= lenj - 2 {
                                        ix!(Iw[pdst]) = ix!(Iw[psrc]);
                                        pdst += 1;
                                        psrc += 1;
                                        knt3 += 1;
                                    }
                                }
                            }
                            let p1 = pdst;
                            let mut psrc = pme1;
                            while psrc <= pfree - 1 {
                                ix!(Iw[pdst]) = ix!(Iw[psrc]);
                                pdst += 1;
                                psrc += 1;
                            }
                            pme1 = p1;
                            pfree = pdst;
                            pj = ix!(Pe[e]);
                            p = Pe[meu];
                        }
                        degme += nvi;
                        Nv[iu] = -nvi;
                        ix!(Iw[pfree]) = i;
                        pfree += 1;
                        let ilast = Last[iu];
                        let inext = Next[iu];
                        if inext != EMPTY {
                            ix!(Last[inext]) = ilast;
                        }
                        if ilast != EMPTY {
                            ix!(Next[ilast]) = inext;
                        } else {
                            ix!(Head[Degree[iu]]) = inext;
                        }
                    }
                    knt2 += 1;
                }
                if e != me {
                    ix!(Pe[e]) = flip(me);
                    ix!(W[e]) = 0;
                }
                knt1 += 1;
            }
            pme2 = pfree - 1;
        }

        Degree[meu] = degme;
        Pe[meu] = pme1;
        Len[meu] = pme2 - pme1 + 1;
        Elen[meu] = flip(nvpiv + degme);

        // make sure that wflg is not too large
        if wflg >= wbig {
            for x in 0..n as usize {
                if W[x] != 0 {
                    W[x] = 1;
                }
            }
            wflg = 2;
        }

        // compute (W [e] - wflg) = |Le\Lme| for all elements
        let mut pme = pme1;
        while pme <= pme2 {
            let i = ix!(Iw[pme]);
            let iu = i as usize;
            let eln = Elen[iu];
            if eln > 0 {
                let nvi = -Nv[iu];
                let wnvi = wflg - nvi;
                let mut p = Pe[iu];
                while p <= Pe[iu] + eln - 1 {
                    let e = ix!(Iw[p]);
                    let eu = e as usize;
                    let mut we = W[eu];
                    if we >= wflg {
                        we -= nvi;
                    } else if we != 0 {
                        we = Degree[eu] + wnvi;
                    }
                    W[eu] = we;
                    p += 1;
                }
            }
            pme += 1;
        }

        // degree update and element absorption
        let mut pme = pme1;
        while pme <= pme2 {
            let i = ix!(Iw[pme]);
            let iu = i as usize;
            let p1 = Pe[iu];
            let p2 = p1 + Elen[iu] - 1;
            let mut pn = p1;
            let mut hash: u32 = 0;
            let mut deg: Int = 0;
            if aggressive {
                let mut p = p1;
                while p <= p2 {
                    let e = ix!(Iw[p]);
                    let eu = e as usize;
                    let we = W[eu];
                    if we != 0 {
                        let dext = we - wflg;
                        if dext > 0 {
                            deg += dext;
                            ix!(Iw[pn]) = e;
                            pn += 1;
                            hash = hash.wrapping_add(e as u32);
                        } else {
                            Pe[eu] = flip(me);
                            W[eu] = 0;
                        }
                    }
                    p += 1;
                }
            } else {
                let mut p = p1;
                while p <= p2 {
                    let e = ix!(Iw[p]);
                    let we = ix!(W[e]);
                    if we != 0 {
                        let dext = we - wflg;
                        deg += dext;
                        ix!(Iw[pn]) = e;
                        pn += 1;
                        hash = hash.wrapping_add(e as u32);
                    }
                    p += 1;
                }
            }
            Elen[iu] = pn - p1 + 1;
            let p3 = pn;
            let p4 = p1 + Len[iu];
            let mut p = p2 + 1;
            while p < p4 {
                let j = ix!(Iw[p]);
                let nvj = ix!(Nv[j]);
                if nvj > 0 {
                    deg += nvj;
                    ix!(Iw[pn]) = j;
                    pn += 1;
                    hash = hash.wrapping_add(j as u32);
                }
                p += 1;
            }
            if Elen[iu] == 1 && p3 == pn {
                // mass elimination
                Pe[iu] = flip(me);
                let nvi = -Nv[iu];
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
                Nv[iu] = 0;
                Elen[iu] = EMPTY;
            } else {
                Degree[iu] = Degree[iu].min(deg);
                ix!(Iw[pn]) = ix!(Iw[p3]);
                ix!(Iw[p3]) = ix!(Iw[p1]);
                ix!(Iw[p1]) = me;
                Len[iu] = pn - p1 + 1;
                hash %= n as u32;
                let j = Head[hash as usize];
                if j <= EMPTY {
                    Next[iu] = flip(j);
                    Head[hash as usize] = flip(i);
                } else {
                    Next[iu] = ix!(Last[j]);
                    ix!(Last[j]) = i;
                }
                Last[iu] = hash as Int;
            }
            pme += 1;
        }

        Degree[meu] = degme;

        // clear the counter array, W [...], by incrementing wflg
        lemax = lemax.max(degme);
        wflg += lemax;
        if wflg >= wbig {
            for x in 0..n as usize {
                if W[x] != 0 {
                    W[x] = 1;
                }
            }
            wflg = 2;
        }

        // supervariable detection
        let mut pme = pme1;
        while pme <= pme2 {
            let mut i = ix!(Iw[pme]);
            if ix!(Nv[i]) < 0 {
                let hash = ix!(Last[i]);
                let j = ix!(Head[hash]);
                if j == EMPTY {
                    i = EMPTY;
                } else if j < EMPTY {
                    i = flip(j);
                    ix!(Head[hash]) = EMPTY;
                } else {
                    i = ix!(Last[j]);
                    ix!(Last[j]) = EMPTY;
                }
                while i != EMPTY && ix!(Next[i]) != EMPTY {
                    let iu = i as usize;
                    let ln = Len[iu];
                    let eln = Elen[iu];
                    let mut p = Pe[iu] + 1;
                    while p <= Pe[iu] + ln - 1 {
                        ix!(W[ix!(Iw[p])]) = wflg;
                        p += 1;
                    }
                    let mut jlast = i;
                    let mut j = Next[iu];
                    while j != EMPTY {
                        let ju = j as usize;
                        let mut ok = (Len[ju] == ln) && (Elen[ju] == eln);
                        let mut p = Pe[ju] + 1;
                        while ok && p <= Pe[ju] + ln - 1 {
                            if ix!(W[ix!(Iw[p])]) != wflg {
                                ok = false;
                            }
                            p += 1;
                        }
                        if ok {
                            Pe[ju] = flip(i);
                            Nv[iu] += Nv[ju];
                            Nv[ju] = 0;
                            Elen[ju] = EMPTY;
                            j = Next[ju];
                            ix!(Next[jlast]) = j;
                        } else {
                            jlast = j;
                            j = Next[ju];
                        }
                    }
                    wflg += 1;
                    i = Next[iu];
                }
            }
            pme += 1;
        }

        // restore degree lists and remove nonprincipal supervariables
        let mut p = pme1;
        let nleft = n - nel;
        let mut pme = pme1;
        while pme <= pme2 {
            let i = ix!(Iw[pme]);
            let iu = i as usize;
            let nvi = -Nv[iu];
            if nvi > 0 {
                Nv[iu] = nvi;
                let mut deg = Degree[iu] + degme - nvi;
                deg = deg.min(nleft - nvi);
                let inext = ix!(Head[deg]);
                if inext != EMPTY {
                    ix!(Last[inext]) = i;
                }
                Next[iu] = inext;
                Last[iu] = EMPTY;
                ix!(Head[deg]) = i;
                mindeg = mindeg.min(deg);
                Degree[iu] = deg;
                ix!(Iw[p]) = i;
                p += 1;
            }
            pme += 1;
        }

        Nv[meu] = nvpiv;
        Len[meu] = p - pme1;
        if Len[meu] == 0 {
            Pe[meu] = EMPTY;
            W[meu] = 0;
        }
        if elenme != 0 {
            pfree = p;
        }

        if info.is_some() {
            let f = nvpiv as f64;
            let r = (degme + ndense) as f64;
            dmax = dmax.max(f + r);
            let lnzme = f * r + (f - 1.0) * f / 2.0;
            lnz += lnzme;
            ndiv += lnzme;
            let s = f * r * r + r * (f - 1.0) * f + (f - 1.0) * f * (2.0 * f - 1.0) / 6.0;
            nms_lu += s;
            nms_ldl += (s + lnzme) / 2.0;
        }
    }

    if let Some(info) = info {
        let f = ndense as f64;
        dmax = dmax.max(ndense as f64);
        let lnzme = (f - 1.0) * f / 2.0;
        lnz += lnzme;
        ndiv += lnzme;
        let s = (f - 1.0) * f * (2.0 * f - 1.0) / 6.0;
        nms_lu += s;
        nms_ldl += (s + lnzme) / 2.0;
        info[AMD_LNZ] = lnz;
        info[AMD_NDIV] = ndiv;
        info[AMD_NMULTSUBS_LDL] = nms_ldl;
        info[AMD_NMULTSUBS_LU] = nms_lu;
        info[AMD_NDENSE] = ndense as f64;
        info[AMD_DMAX] = dmax;
        info[AMD_NCMPA] = ncmpa as f64;
        info[AMD_STATUS] = AMD_OK as f64;
    }

    // post-ordering
    for i in 0..n as usize {
        Pe[i] = flip(Pe[i]);
        Elen[i] = flip(Elen[i]);
    }
    for i in 0..n as usize {
        if Nv[i] == 0 {
            let mut j = Pe[i];
            if j == EMPTY {
                continue;
            }
            while ix!(Nv[j]) == 0 {
                j = ix!(Pe[j]);
            }
            let e = j;
            let mut j = i as Int;
            while ix!(Nv[j]) == 0 {
                let jnext = ix!(Pe[j]);
                ix!(Pe[j]) = e;
                j = jnext;
            }
        }
    }

    amd_postorder(n, Pe, Nv, Elen, W, Head, Next, Last);

    for k in 0..n as usize {
        Head[k] = EMPTY;
        Next[k] = EMPTY;
    }
    for e in 0..n {
        let k = ix!(W[e]);
        if k != EMPTY {
            ix!(Head[k]) = e;
        }
    }
    nel = 0;
    for k in 0..n as usize {
        let e = Head[k];
        if e == EMPTY {
            break;
        }
        ix!(Next[e]) = nel;
        nel += ix!(Nv[e]);
    }
    for i in 0..n as usize {
        if Nv[i] == 0 {
            let e = Pe[i];
            if e != EMPTY {
                Next[i] = ix!(Next[e]);
                ix!(Next[e]) += 1;
            } else {
                Next[i] = nel;
                nel += 1;
            }
        }
    }
    for i in 0..n {
        let k = ix!(Next[i]);
        ix!(Last[k]) = i;
    }
}
