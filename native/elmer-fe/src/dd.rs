//! Double-double arithmetic standing in for gfortran's `REAL(KIND=qp)`
//! (binary128) in `CheckMetric` and `ElementMetricQP`. It carries ~106
//! significand bits against binary128's 113, so after the final rounding to
//! f64 it agrees with Elmer except where a result lies within ~2^-50 ulp
//! of a rounding tie.

#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct Dd {
    pub hi: f64,
    pub lo: f64,
}

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

#[inline]
fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

impl Dd {
    pub const ZERO: Dd = Dd { hi: 0.0, lo: 0.0 };

    pub fn from(x: f64) -> Dd {
        Dd { hi: x, lo: 0.0 }
    }

    pub fn to_f64(self) -> f64 {
        self.hi + self.lo
    }

    pub fn add(self, b: Dd) -> Dd {
        let (s, e) = two_sum(self.hi, b.hi);
        let (t, f) = two_sum(self.lo, b.lo);
        let (s, e) = quick_two_sum(s, e + t);
        let (hi, lo) = quick_two_sum(s, e + f);
        Dd { hi, lo }
    }

    pub fn neg(self) -> Dd {
        Dd {
            hi: -self.hi,
            lo: -self.lo,
        }
    }

    pub fn sub(self, b: Dd) -> Dd {
        self.add(b.neg())
    }

    pub fn mul(self, b: Dd) -> Dd {
        let (p, e) = two_prod(self.hi, b.hi);
        let e = e + (self.hi * b.lo + self.lo * b.hi);
        let (hi, lo) = quick_two_sum(p, e);
        Dd { hi, lo }
    }

    pub fn div(self, b: Dd) -> Dd {
        let q1 = self.hi / b.hi;
        let r = self.sub(b.mul(Dd::from(q1)));
        let q2 = r.hi / b.hi;
        let r = r.sub(b.mul(Dd::from(q2)));
        let q3 = r.hi / b.hi;
        let (hi, lo) = quick_two_sum(q1, q2);
        Dd { hi, lo }.add(Dd::from(q3))
    }

    pub fn sqrt(self) -> Dd {
        if self.hi <= 0.0 {
            return Dd::from(self.hi.sqrt());
        }
        let x = 1.0 / self.hi.sqrt();
        let ax = self.hi * x;
        let d = self.sub(Dd::from(ax).mul(Dd::from(ax)));
        let (hi, lo) = two_sum(ax, d.hi * (x * 0.5));
        Dd { hi, lo }
    }

    pub fn abs(self) -> Dd {
        if self.hi < 0.0 {
            self.neg()
        } else {
            self
        }
    }

    pub fn le(self, b: Dd) -> bool {
        self.hi < b.hi || (self.hi == b.hi && self.lo <= b.lo)
    }
}
