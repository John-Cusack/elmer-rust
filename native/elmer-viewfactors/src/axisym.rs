//! View factors of an axisymmetric geometry: a line-by-line port of
//! `fem/src/viewaxis/viewfact.cpp` (`viewfactorsaxis`, by Juha Katajamaki).
//!
//! Every boundary element is a segment of the (r, z) half plane, which
//! stands for a conical ring. The factor of a pair is integrated in the
//! azimuth angle (its cosine, bounded by the visibility conditions and by the
//! shadows of the other rings) and by three-point Gauss quadrature along the
//! viewing segment. The constants of the C++ (`eps`, the 3-point rule with
//! its nine-digit weights, `PI = 3.1415926535`) are kept as they are: the
//! pristine `ElmerSolver` normalises the factors it reads, and the aim is the
//! same numbers, not merely close ones.

const EPS: f64 = 1e-7;
const EPS2: f64 = 1.0e-14;
const DELTA: f64 = 1e-6;
const SELF_SHADING: bool = true;

fn mx(x: f64, y: f64) -> f64 {
    if x > y {
        x
    } else {
        y
    }
}

fn mn(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

fn sgn(x: f64) -> f64 {
    if x < 0.0 {
        -1.0
    } else if x > 0.0 {
        1.0
    } else {
        0.0
    }
}

/// `Area(r1, r2, z1, z2)`: the lateral area of a conical ring.
#[allow(clippy::approx_constant)] // The C++ literal, kept for identical numbers.
fn area(r1: f64, r2: f64, z1: f64, z2: f64) -> f64 {
    3.1415926535 * (r1 + r2) * ((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2)).sqrt()
}

/// The file-scope statics of the C++, one set per computation.
#[derive(Default)]
struct State<'a> {
    // The viewed ring (1, 2) and the viewing ring (3, 4).
    r1: f64,
    r2: f64,
    z1: f64,
    z2: f64,
    r3: f64,
    r4: f64,
    z3: f64,
    z4: f64,
    r12: f64,
    r34: f64,
    z12: f64,
    z34: f64,
    zd1: f64,
    zd3: f64,
    zd: f64,
    rd1: f64,
    rd3: f64,
    g1: f64,
    g3: f64,
    d1: f64,
    d3: f64,
    t: f64,
    rratio: f64,
    inode: usize,
    jnode: usize,
    nsurf: usize,
    coord: &'a [f64],
    /// The shading rings, `[2k]` and `[2k + 1]` their end nodes.
    shade: Vec<usize>,
    nshade: usize,
    shade_parent: Vec<i64>,
}

impl State<'_> {
    fn c(&self, node: usize) -> (f64, f64) {
        (self.coord[2 * node], self.coord[2 * node + 1])
    }

    /// `InitialInterval`: the range of the cosine of the azimuth that keeps
    /// the connecting line within 90 degrees of both normals.
    fn initial_interval(&self) -> Option<(f64, f64)> {
        let (mut c1, mut c2) = (-1.0, 1.0);
        let (zd1, zd3, zd, rd1, rd3, r12, r34) = (
            self.zd1, self.zd3, self.zd, self.rd1, self.rd3, self.r12, self.r34,
        );
        if zd1.abs() > EPS {
            let cc1 = (-zd * rd1 + r12 * zd1) / (r34 * zd1);
            if zd3.abs() > EPS {
                let cc3 = (zd * rd3 + r34 * zd3) / (r12 * zd3);
                if zd1 > 0.0 {
                    if zd3 > 0.0 {
                        c1 = mx(cc1, cc3);
                    } else {
                        c1 = cc1;
                        c2 = cc3;
                    }
                } else if zd3 < 0.0 {
                    c2 = mn(cc1, cc3);
                } else {
                    c1 = cc3;
                    c2 = cc1;
                }
            } else if sgn(rd3) != 0.0 && sgn(rd3) == -sgn(zd) {
                if zd1 > 0.0 {
                    c1 = cc1;
                } else {
                    c2 = cc1;
                }
            } else {
                c1 = 1.0;
                c2 = -1.0;
            }
        } else if zd3.abs() > EPS {
            let cc3 = (zd * rd3 + r34 * zd3) / (r12 * zd3);
            if sgn(rd1) != 0.0 && sgn(rd1) == sgn(zd) {
                if zd3 > 0.0 {
                    c1 = cc3;
                } else {
                    c2 = cc3;
                }
            } else {
                c1 = 1.0;
                c2 = -1.0;
            }
        } else if sgn(rd1) == 0.0 || sgn(rd1) != sgn(zd) || sgn(rd1) != -sgn(rd3) {
            c1 = 1.0;
            c2 = -1.0;
        }
        c1 = mx(-1.0 + EPS, c1);
        c2 = mn(1.0 - EPS, c2);
        if c2 - c1 < EPS {
            None
        } else {
            Some((c1, c2))
        }
    }

    /// `ExaminePoint`: the cosine of the azimuth at which the shadow cone
    /// meets the line at parameter `x`, widening `[mi, ma]`.
    fn examine_point(&mut self, x: f64, mi: &mut f64, ma: &mut f64) {
        let y;
        if x > EPS {
            if 1.0 - x > EPS {
                self.t = self.rratio * x / (1.0 - x);
                y = 0.5 * (self.d1 / self.t + self.d3 * self.t) + self.g1 * self.g3;
            } else if self.d3.abs() < EPS {
                y = self.g1 * self.g3;
            } else {
                y = sgn(self.d3);
            }
        } else if self.d1.abs() < EPS {
            y = self.g1 * self.g3;
        } else {
            y = sgn(self.d1);
        }
        if y > *ma {
            *ma = y;
        }
        if y < *mi {
            *mi = y;
        }
    }

    /// `ViewIntegral`: the integral over `[c1, c2]` of the azimuth cosine,
    /// split around the intervals the shading rings `k..` hide.
    fn view_integral(&mut self, mut c1: f64, mut c2: f64, mut k: usize) -> f64 {
        self.rratio = self.r34 / self.r12;
        while k < self.nshade {
            let r5 = self.coord[2 * self.shade[2 * k + 1]];
            let r6 = self.coord[2 * self.shade[2 * k]];
            let z5 = self.coord[2 * self.shade[2 * k + 1] + 1];
            let z6 = self.coord[2 * self.shade[2 * k] + 1];
            k += 1;

            if SELF_SHADING {
                if self.nsurf != self.nshade {
                    if k as i64 - 1 == self.shade_parent[self.inode]
                        || k as i64 - 1 == self.shade_parent[self.jnode]
                    {
                        continue;
                    }
                } else if k - 1 == self.inode || k - 1 == self.jnode {
                    continue;
                }
            }

            if r5 + r6 < EPS {
                continue;
            }

            let zd5 = z5 - z6;
            let (mut cc1, mut cc2);
            if zd5.abs() < EPS {
                // The shading surface is a flat ring.
                if self.nsurf == self.nshade && self.inode == k - 1 {
                    continue;
                }
                if self.zd.abs() < EPS {
                    continue;
                }
                let t1 = (self.z12 - z5) / self.zd;
                let tt1 = 1.0 - t1;
                if t1 < EPS || tt1 < EPS {
                    continue;
                }
                self.t = self.rratio * t1 / tt1;
                cc1 = 0.5 * (r5 * r5 / (self.r12 * self.r34 * t1 * tt1) - self.t - 1.0 / self.t);
                cc2 = 0.5 * (r6 * r6 / (self.r12 * self.r34 * t1 * tt1) - self.t - 1.0 / self.t);
                if cc1 > cc2 {
                    self.t = cc1;
                    cc1 = cc2;
                    cc2 = self.t;
                }
            } else {
                // A cone or a cylinder: the span of the connecting line in z
                // that falls in its shadow.
                let (mut t1, mut t2);
                if self.zd.abs() < EPS {
                    if (self.z12 - z5 < EPS && self.z12 - z6 > EPS)
                        || (self.z12 - z5 > EPS && self.z12 - z6 < EPS)
                    {
                        t1 = 0.0;
                        t2 = 1.0;
                    } else {
                        continue;
                    }
                } else {
                    t1 = (self.z12 - z5) / self.zd;
                    t2 = t1 + zd5 / self.zd;
                    if t1 > t2 {
                        std::mem::swap(&mut t1, &mut t2);
                    }
                }
                match interval_isect(0.0, 1.0, t1, t2) {
                    Some((a, b)) => {
                        t1 = a;
                        t2 = b;
                    }
                    None => continue,
                }

                // The values the azimuth cosine takes between t1 and t2.
                cc1 = 1.0;
                cc2 = -1.0;
                self.g1 = (r5 * (self.z12 - z6) - r6 * (self.z12 - z5)) / (self.r12 * zd5);
                self.g3 = (r5 * (self.z34 - z6) - r6 * (self.z34 - z5)) / (self.r34 * zd5);
                self.d1 = self.g1 * self.g1 - 1.0;
                self.d3 = self.g3 * self.g3 - 1.0;

                self.examine_point(t1, &mut cc1, &mut cc2);
                self.examine_point(t2, &mut cc1, &mut cc2);

                // Both points outside the cone: the extremum of the
                // derivative, if it is inside [t1, t2].
                if self.d1 <= -EPS && self.d3 <= -EPS {
                    let t0 = 1.0 / (1.0 + (self.rratio * self.d3 / self.d1).sqrt());
                    if t0 - t1 > EPS && t2 - t0 > EPS {
                        self.examine_point(t0, &mut cc1, &mut cc2);
                    }
                }
                if cc1 > cc2 {
                    cc1 = cc2;
                }
            }

            if let Some((i1, i2)) = interval_isect(c1, c2, cc1, cc2) {
                cc1 = i1;
                cc2 = i2;
                if cc1 - c1 < DELTA {
                    if c2 - cc2 < DELTA {
                        return 0.0;
                    }
                    c1 = cc2;
                } else if c2 - cc2 < DELTA {
                    c2 = cc1;
                } else {
                    let a = self.view_integral(c1, cc1, k);
                    let b = self.view_integral(cc2, c2, k);
                    return a + b;
                }
            }
        }
        self.integrate(c1, c2)
    }

    /// `Integrate`: three-point Gauss quadrature along the viewing ring of
    /// the closed-form azimuth integral between the cosines `c1` and `c2`.
    fn integrate(&self, c1: f64, c2: f64) -> f64 {
        const QP: [f64; 3] = [0.112701665, 0.5, 0.887298334];
        const W: [f64; 3] = [0.277777777, 0.444444444, 0.277777777];
        let c = self.zd1 * self.zd1 + self.rd1 * self.rd1;
        if c < EPS2 {
            return 0.0;
        }
        let (r1, r2, z1, z2, r3, z3) = (self.r1, self.r2, self.z1, self.z2, self.r3, self.z3);
        let (rd1, rd3, zd3) = (self.rd1, self.rd3, self.zd3);
        let zrd = r2 * z1 - r1 * z2;
        let a1 = rd3 * r1;
        let a2 = rd3 * r2;
        let b1 = zd3 * z1;
        let b2 = zd3 * z2;
        let s1 = (1. - c1 * c1).sqrt();
        let s2 = (1. - c2 * c2).sqrt();
        // In the cosines and sines indices 1 and 2 are the other way round
        // than in the other variables.
        let cs = (1. + c1) * (1. + c2);
        let cd = (1. - c1) * (1. - c2);

        let mut integral = 0.0;
        for i in 0..3 {
            let z = z3 - QP[i] * zd3;
            let r = r3 - QP[i] * rd3;
            let mut e1 = (z1 - z) * (z1 - z) + r1 * r1 + r * r;
            let mut f1 = (z2 - z) * (z2 - z) + r2 * r2 + r * r;
            let hh1 = 2. * r1 * r;
            let hh2 = 2. * r2 * r;
            let g1 = -e1 / hh1;
            let g2 = -f1 / hh2;
            let e2 = e1 - c1 * hh1;
            let f2 = f1 - c1 * hh2;
            e1 -= c2 * hh1;
            f1 -= c2 * hh2;
            let h = zd3 * z + rd3 * r;
            let gg1 = (g1 + c2) * (g1 + c1);
            let gg2 = (g2 + c2) * (g2 + c1);

            // The arcs.
            let mut value = (-0.5 * (a1 + (h - b1) * g1) / (g1 * g1 - 1.).sqrt())
                * (0.5 * ((1. - g1) * (cd / gg1).sqrt() - (1. + g1) * (cs / gg1).sqrt())).acos();
            value -= (-0.5 * (a2 + (h - b2) * g2) / (g2 * g2 - 1.).sqrt())
                * (0.5 * ((1. - g2) * (cd / gg2).sqrt() - (1. + g2) * (cs / gg2).sqrt())).acos();
            value += 0.25 * (b1 - b2) * (c1 * c2 + s1 * s2).acos();

            // The straight sides.
            let gg1 = e1 + f1 - c;
            let gg2 = e2 + f2 - c;
            let hh1 = 4. * e1 * f1;
            let hh2 = 4. * e2 * f2;
            let d1 = hh1 - gg1 * gg1;
            let d2 = hh2 - gg2 * gg2;
            let h = r * (rd1 * h + zrd * zd3);
            value -= h * (s1 / d2.sqrt()) * (gg2 / hh2.sqrt()).acos();
            value += h * (s2 / d1.sqrt()) * (gg1 / hh1.sqrt()).acos();

            integral += W[i] * value;
        }
        integral
    }
}

/// `IntervalIsect`: the intersection of `[x1, x2]` and `[y1, y2]`, `None`
/// when it is empty or of length below `EPS`.
fn interval_isect(x1: f64, x2: f64, y1: f64, y2: f64) -> Option<(f64, f64)> {
    let (mut z1, mut z2) = (x1, x2);
    if x2 - y1 < EPS {
        return None;
    }
    if y1 - x1 > EPS {
        z1 = y1;
    }
    if y2 - x1 < EPS {
        return None;
    }
    if x2 - y2 > EPS {
        z2 = y2;
    }
    (z2 - z1 >= EPS).then_some((z1, z2))
}

/// The view factor matrix of `n` conical rings, row-major `n * n`, row `i`
/// being the factors from ring `i` ("viewing", nodes `surf[2i+1]`,
/// `surf[2i]`) to every ring `j`.
///
/// * `surf`: two node indices per ring (0-based), oriented so that the ring's
///   normal points where the ring radiates;
/// * `coord`: `(r, z)` of every node, interleaved;
/// * `divide`: `Viewfactor divide`, the number of pieces each ring is cut
///   into;
/// * `combine`: `Viewfactor combine elements`, merging collinear neighbours
///   into one shading segment (the Fortran passes `CombineInt`).
pub fn view_factors_axis(surf: &[usize], coord: &[f64], divide: usize, combine: bool) -> Vec<f64> {
    let n = surf.len() / 2;
    let div = divide.max(1);
    let mut st = State {
        nsurf: n,
        coord,
        ..State::default()
    };

    if !combine {
        st.shade = surf.to_vec();
        st.nshade = n;
    } else {
        let epsilon = 1.0e-5;
        let maxind = surf.iter().copied().max().unwrap_or(0);
        let mut nodehits = vec![0usize; maxind + 1];
        for &s in surf {
            nodehits[s] += 1;
        }
        let maxnodehits = nodehits.iter().copied().max().unwrap_or(0);
        let mut nodetable = vec![0usize; (maxind + 1) * maxnodehits];
        nodehits.iter_mut().for_each(|h| *h = 0);
        for i in 0..n {
            let ind1 = surf[2 * i + 1];
            let ind2 = surf[2 * i];
            nodetable[maxnodehits * ind1 + nodehits[ind1]] = i;
            nodetable[maxnodehits * ind2 + nodehits[ind2]] = i;
            nodehits[ind1] += 1;
            nodehits[ind2] += 1;
        }

        let mut shade = surf.to_vec();
        for ind0 in 0..=maxind {
            if nodehits[ind0] != 2 {
                continue;
            }
            let elem1 = nodetable[maxnodehits * ind0];
            let ind1 = if shade[2 * elem1 + 1] == ind0 {
                shade[2 * elem1]
            } else {
                shade[2 * elem1 + 1]
            };
            let elem2 = nodetable[maxnodehits * ind0 + 1];
            let ind2 = if shade[2 * elem2 + 1] == ind0 {
                shade[2 * elem2]
            } else {
                shade[2 * elem2 + 1]
            };

            let (r0, z0) = st.c(ind0);
            let (r1, z1) = st.c(ind1);
            let (r2, z2) = st.c(ind2);
            let (dr1, dr2) = (r1 - r0, r2 - r0);
            let (dz1, dz2) = (z1 - z0, z2 - z0);
            let mut dp1 = dr1 * dr2 + dz1 * dz2;
            let ds1 = (dr1 * dr1 + dz1 * dz1).sqrt();
            let ds2 = (dr2 * dr2 + dz2 * dz2).sqrt();
            dp1 /= ds1 * ds2;

            // The segments must be aligned.
            if dp1 > epsilon - 1. {
                continue;
            }

            // Make the first element bigger, destroy the second.
            if shade[2 * elem1] == ind0 {
                shade[2 * elem1] = ind2;
            } else {
                shade[2 * elem1 + 1] = ind2;
            }
            shade[2 * elem2] = 0;
            shade[2 * elem2 + 1] = 0;

            nodehits[ind0] = 0;
            if nodetable[maxnodehits * ind2] == elem2 {
                nodetable[maxnodehits * ind2] = elem1;
            } else {
                nodetable[maxnodehits * ind2 + 1] = elem1;
            }
        }

        // Compact the surviving segments.
        let mut j = 0;
        for i in 0..n {
            if shade[2 * i + 1] != 0 || shade[2 * i] != 0 {
                shade[2 * j + 1] = shade[2 * i + 1];
                shade[2 * j] = shade[2 * i];
                j += 1;
            }
        }
        st.nshade = j;
        shade.truncate(2 * j);
        st.shade = shade;

        // The superelement each original ring belongs to.
        st.shade_parent = vec![-1; n];
        for i in 0..st.nshade {
            let ind1 = st.shade[2 * i];
            let ind2 = st.shade[2 * i + 1];
            let (r1, z1) = st.c(ind1);
            let (r2, z2) = st.c(ind2);
            let mut dz1 = z2 - z1;
            let mut dr1 = r2 - r1;
            let ds1 = (dz1 * dz1 + dr1 * dr1).sqrt();
            dz1 /= ds1;
            dr1 /= ds1;

            for jj in 0..n {
                if st.shade_parent[jj] >= 0 {
                    continue;
                }
                let mut hit = true;
                for k in 0..2 {
                    let ind0 = surf[2 * jj + k];
                    // The original code compares with ind1 twice (ind2 was
                    // probably meant); the comparison is kept as it is.
                    if ind0 == ind1 {
                        continue;
                    }
                    let (r3, z3) = st.c(ind0);
                    let dz2 = z3 - z1;
                    let dr2 = r3 - r1;
                    let ds2 = (dz2 * dz2 + dr2 * dr2).sqrt();
                    let dp1 = dz1 * dz2 + dr1 * dr2;
                    if dp1 / ds2 < 1. - epsilon {
                        hit = false;
                        break;
                    }
                    if dp1 / ds1 > 1. + epsilon || dp1 / ds1 < -epsilon {
                        hit = false;
                        break;
                    }
                }
                if hit {
                    st.shade_parent[jj] = i as i64;
                }
            }
        }
    }

    let mut vf = vec![0.0; n * n];
    for i in 0..n {
        st.inode = i;
        let (r3, z3) = st.c(surf[2 * i + 1]);
        let (r4, z4) = st.c(surf[2 * i]);
        let (_r3, _r4, _z3, _z4) = (r3, r4, z3, z4);
        let a = area(_r3, _r4, _z3, _z4);

        for j in 0..n {
            st.jnode = j;
            let (_r1, _z1) = st.c(surf[2 * j + 1]);
            let (_r2, _z2) = st.c(surf[2 * j]);

            if a < EPS {
                continue;
            }

            let d = div as f64;
            for ii in 0..div {
                let (fi, fi1, fr) = (ii as f64, (ii + 1) as f64, (div - ii) as f64);
                let fr1 = (div - ii - 1) as f64;
                st.r3 = _r3 * fr / d + _r4 * fi / d;
                st.r4 = _r3 * fr1 / d + _r4 * fi1 / d;
                st.z3 = _z3 * fr / d + _z4 * fi / d;
                st.z4 = _z3 * fr1 / d + _z4 * fi1 / d;

                st.r34 = 0.5 * (st.r3 + st.r4);
                st.z34 = 0.5 * (st.z3 + st.z4);
                st.zd3 = st.z3 - st.z4;
                st.rd3 = st.r3 - st.r4;

                for jj in 0..div {
                    let (fj, fj1, fq) = (jj as f64, (jj + 1) as f64, (div - jj) as f64);
                    let fq1 = (div - jj - 1) as f64;
                    st.r1 = _r1 * fq / d + _r2 * fj / d;
                    st.r2 = _r1 * fq1 / d + _r2 * fj1 / d;
                    st.z1 = _z1 * fq / d + _z2 * fj / d;
                    st.z2 = _z1 * fq1 / d + _z2 * fj1 / d;
                    st.r12 = 0.5 * (st.r1 + st.r2);
                    if st.r12 < EPS || st.r34 < EPS {
                        continue;
                    }

                    if st.r1 < EPS {
                        st.r1 = EPS;
                    }
                    if st.r2 < EPS {
                        st.r2 = EPS;
                    }
                    if st.r3 < EPS {
                        st.r3 = EPS;
                    }
                    if st.r4 < EPS {
                        st.r4 = EPS;
                    }

                    st.zd1 = st.z1 - st.z2;
                    st.rd1 = st.r1 - st.r2;
                    st.z12 = 0.5 * (st.z1 + st.z2);
                    st.zd = st.z12 - st.z34;

                    let Some((c1, c2)) = st.initial_interval() else {
                        continue;
                    };
                    let viewint = st.view_integral(c1, c2, 0);
                    // 4 = 2 (mirror symmetry) * 2 pi (rotation) / pi (the
                    // constant of the integral).
                    vf[i * n + j] += 4. * viewint;
                }
            }
            vf[i * n + j] /= a;
        }
    }
    vf
}
