//! `coreGeometry.h`: the robust predicates and helpers of Mesh2D, with the
//! operation order of the C++ so results agree bit for bit.

const DBL_EPSILON: f64 = f64::EPSILON;

fn fmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// `inCircle(x[4], y[4])`: 1 unless point 0 is outside the circle through
/// points 1-3 by more than the rounding error.
pub fn in_circle(x: &[f64; 4], y: &[f64; 4]) -> i32 {
    let mut x10 = x[1] - x[0];
    let mut x20 = x[2] - x[0];
    let mut x30 = x[3] - x[0];
    let mut x21 = x[2] - x[1];
    let mut x31 = x[3] - x[1];
    let mut y10 = y[1] - y[0];
    let mut y20 = y[2] - y[0];
    let mut y30 = y[3] - y[0];
    let mut y21 = y[2] - y[1];
    let mut y31 = y[3] - y[1];

    let mut a1 = x10 * y20;
    let mut b1 = x20 * y10;
    let mut a2 = x10 * y30;
    let mut b2 = x30 * y10;
    let mut a3 = x20 * x21;
    let mut b3 = y20 * y21;
    let mut a4 = x30 * x31;
    let mut b4 = y30 * y31;

    let mut t1 = a1 - b1;
    let mut t2 = a2 - b2;
    let mut t3 = a3 + b3;
    let mut t4 = a4 + b4;
    let mut t14 = t1 * t4;
    let mut t23 = t2 * t3;
    let det = t14 - t23;

    if det >= 0.0 {
        return 1;
    }

    let maxcoord = fmax(
        fmax(fmax(x[0].abs(), x[1].abs()), fmax(x[2].abs(), x[3].abs())),
        fmax(fmax(y[0].abs(), y[1].abs()), fmax(y[2].abs(), y[3].abs())),
    );

    x10 = x10.abs();
    y10 = y10.abs();

    x20 = x20.abs();
    y20 = y20.abs();
    a1 = a1.abs();
    b1 = b1.abs();
    t1 = t1.abs();

    let e1 = maxcoord * (x10 + x20 + y10 + y20) + a1 + b1 + fmax(fmax(a1, b1), t1);

    x30 = x30.abs();
    y30 = y30.abs();
    a2 = a2.abs();
    b2 = b2.abs();
    t2 = t2.abs();

    let e2 = maxcoord * (x10 + x30 + y10 + y30) + a2 + b2 + fmax(fmax(a2, b2), t2);

    x21 = x21.abs();
    y21 = y21.abs();
    a3 = a3.abs();
    b3 = b3.abs();
    t3 = t3.abs();

    let e3 = maxcoord * (x21 + x20 + y21 + y20) + a3 + b3 + fmax(fmax(a3, b3), t3);

    x31 = x31.abs();
    y31 = y31.abs();
    a4 = a4.abs();
    b4 = b4.abs();
    t4 = t4.abs();

    let e4 = maxcoord * (x31 + x30 + y31 + y30) + a4 + b4 + fmax(fmax(a4, b4), t4);

    t14 = t14.abs();
    t23 = t23.abs();
    let mut eps =
        e1 * t4 + e4 * t1 + e2 * t3 + e3 * t2 + t14 + t23 + fmax(fmax(t14, t23), det.abs());

    eps *= 2.0 * DBL_EPSILON;

    if det < -eps {
        return 0;
    }
    1
}

/// `sideDet` with the error estimate.
pub fn side_det(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    x3: f64,
    y3: f64,
    err: Option<&mut f64>,
) -> f64 {
    let d23 = x2 * y3;
    let d32 = x3 * y2;
    let d13 = x1 * y3;
    let d31 = x3 * y1;
    let d12 = x1 * y2;
    let d21 = x2 * y1;
    let d1 = d23 - d32;
    let d2 = d13 - d31;
    let d3 = d12 - d21;

    if let Some(err) = err {
        let e1 = d23.abs() + d32.abs() + fmax(fmax(d23.abs(), d32.abs()), d1.abs());
        let e2 = d13.abs() + d31.abs() + fmax(fmax(d13.abs(), d31.abs()), d2.abs());
        let e3 = d12.abs() + d21.abs() + fmax(fmax(d12.abs(), d21.abs()), d3.abs());

        *err = (e1
            + e2
            + e3
            + fmax(fmax(d1.abs(), d2.abs()), (d1 - d2).abs())
            + fmax(fmax((d1 - d2).abs(), d3.abs()), d1 - d2 + d3))
            * 2.0
            * DBL_EPSILON;
    }

    d1 - d2 + d3
}

pub fn orientation(x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) -> bool {
    let mut eps = 0.0;
    let val = side_det(x1, y1, x2, y2, x3, y3, Some(&mut eps));
    val > eps
}

pub fn distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dx = x2 - x1;
    let dy = y2 - y1;
    (dx * dx + dy * dy).sqrt()
}

/// Circumcentre of three points.
pub fn vertex_location(p1x: f64, p1y: f64, p2x: f64, p2y: f64, p3x: f64, p3y: f64) -> (f64, f64) {
    let a = 2.0 * (p2x - p1x);
    let b = 2.0 * (p2y - p1y);
    let c = 2.0 * (p3x - p2x);
    let d = 2.0 * (p3y - p2y);
    let r1 = -(-p2x * p2x + p1x * p1x - p2y * p2y + p1y * p1y);
    let r2 = -(-p3x * p3x + p2x * p2x - p3y * p3y + p2y * p2y);
    let det = a * d - b * c;
    (
        (1.0 / det) * (d * r1 - b * r2),
        (1.0 / det) * (-c * r1 + a * r2),
    )
}

/// `MAP` / `UNMAP`: the size field is interpolated in log space.
pub fn map(a: f64) -> f64 {
    a.ln()
}

pub fn unmap(a: f64) -> f64 {
    a.exp()
}

/// glibc `rand()` after the default `srand(1)` (TYPE_3 additive feedback),
/// which `std::random_shuffle` of libstdc++ draws from.
pub struct GlibcRand {
    r: Vec<u32>,
    k: usize,
}

impl GlibcRand {
    pub fn new() -> Self {
        let mut r = vec![0u32; 34];
        r[0] = 1;
        for i in 1..31 {
            let prev = r[i - 1] as i32 as i64;
            let hi = prev / 127773;
            let lo = prev % 127773;
            let mut word = 16807 * lo - 2836 * hi;
            if word < 0 {
                word += 2147483647;
            }
            r[i] = word as u32;
        }
        for i in 31..34 {
            r[i] = r[i - 31];
        }
        let mut s = GlibcRand { r, k: 34 };
        for _ in 34..344 {
            s.step();
        }
        s
    }

    fn step(&mut self) -> u32 {
        let v = self.r[self.k - 31].wrapping_add(self.r[self.k - 3]);
        self.r.push(v);
        self.k += 1;
        v
    }

    pub fn rand(&mut self) -> u32 {
        self.step() >> 1
    }

    /// `std::random_shuffle(first, last)` of libstdc++ on an index array.
    pub fn random_shuffle(&mut self, v: &mut [usize]) {
        for i in 1..v.len() {
            let j = (self.rand() as usize) % (i + 1);
            if i != j {
                v.swap(i, j);
            }
        }
    }
}

impl Default for GlibcRand {
    fn default() -> Self {
        Self::new()
    }
}

/// C++ `ostream << double` with `precision(16)`: `%.16g`.
pub fn fmt_g16(x: f64) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    let s = format!("{:.15e}", x);
    let (mant, exp) = s.split_once('e').unwrap();
    let e: i32 = exp.parse().unwrap();
    if !(-4..16).contains(&e) {
        let mant = trim_zeros(mant);
        format!("{mant}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    } else {
        let decimals = (15 - e).max(0) as usize;
        trim_zeros(&format!("{:.*}", decimals, x))
    }
}

fn trim_zeros(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rand_is_glibc_after_srand_1() {
        let mut r = GlibcRand::new();
        assert_eq!(r.rand(), 1804289383);
        assert_eq!(r.rand(), 846930886);
        assert_eq!(r.rand(), 1681692777);
    }

    #[test]
    fn g16_is_printf_percent_g() {
        assert_eq!(fmt_g16(0.0), "0");
        assert_eq!(fmt_g16(2.0), "2");
        assert_eq!(fmt_g16(0.5), "0.5");
        assert_eq!(fmt_g16(-1.25), "-1.25");
        assert_eq!(fmt_g16(0.1), "0.1");
        assert_eq!(fmt_g16(1.0 / 3.0), "0.3333333333333333");
        assert_eq!(fmt_g16(1e-5), "1e-05");
        assert_eq!(fmt_g16(1.5e20), "1.5e+20");
        assert_eq!(fmt_g16(123456789012345.0), "123456789012345");
    }
}
