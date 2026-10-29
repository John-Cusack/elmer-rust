//! Coordinate systems: `SetCoordinateSystem` (ModelDescription.F90:2575)
//! and CoordinateSystems.F90 (metric, `sqrt(det g)`, Christoffel symbols).
//!
//! Elmer keeps the system and `CurrentModel % Dimension` in globals; here
//! they are a value. Axisymmetric and cylindric-symmetric solvers weight
//! integrals by [`CoordinateSystem::sqrt_metric`], i.e. by the radius `x`.

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CoordinateKind {
    Cartesian = 1,
    Cylindric = 2,
    CylindricSymmetric = 3,
    AxisSymmetric = 4,
    Polar = 5,
}

/// The model coordinate system and `CoordinateSystemDimension()`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CoordinateSystem {
    pub kind: CoordinateKind,
    pub dim: usize,
}

impl CoordinateSystem {
    /// `SetCoordinateSystem` from the `Coordinate System` keyword (`None`
    /// = absent, i.e. "cartesian"). `mesh_dim` is the largest
    /// `Mesh % MaxDim` of the model's meshes, which plain "cartesian" and
    /// "polar" inherit. Returns `None` for an unknown name (Elmer: Fatal).
    pub fn from_keyword(csys: Option<&str>, mesh_dim: usize) -> Option<Self> {
        let name = csys.map_or_else(
            || "cartesian".to_string(),
            |s| s.trim().to_ascii_lowercase(),
        );
        let (kind, dim) = match name.as_str() {
            "cartesian" => (CoordinateKind::Cartesian, mesh_dim),
            "cartesian 1d" => (CoordinateKind::Cartesian, 1),
            "cartesian 2d" => (CoordinateKind::Cartesian, 2),
            "cartesian 3d" => (CoordinateKind::Cartesian, 3),
            "axi symmetric" => (CoordinateKind::AxisSymmetric, 2),
            "cylindric symmetric" => (CoordinateKind::CylindricSymmetric, 2),
            "cylindrical" => (CoordinateKind::Cylindric, 3),
            "polar" => (CoordinateKind::Polar, mesh_dim),
            "polar 2d" => (CoordinateKind::Polar, 2),
            "polar 3d" => (CoordinateKind::Polar, 3),
            _ => return None,
        };
        Some(CoordinateSystem { kind, dim })
    }

    /// True for the systems Elmer's solvers treat as rotationally symmetric
    /// (`CurrentCoordinateSystem() == AxisSymmetric .OR. CylindricSymmetric`).
    pub fn is_axisymmetric(&self) -> bool {
        matches!(
            self.kind,
            CoordinateKind::AxisSymmetric | CoordinateKind::CylindricSymmetric
        )
    }

    /// `CoordinateSqrtMetric( X, Y, Z )`.
    pub fn sqrt_metric(&self, x: f64, y: f64, z: f64) -> f64 {
        match self.kind {
            CoordinateKind::Cartesian => 1.0,
            CoordinateKind::Polar => polar_sqrt_metric(x, y, z, self.dim),
            _ => x,
        }
    }

    /// `CoordinateSystemInfo( Metric, SqrtMetric, Symbols, dSymbols, X, Y, Z )`.
    pub fn info(&self, x: f64, y: f64, z: f64) -> CoordinateInfo {
        // Every system starts from the identity metric.
        let mut ci = CoordinateInfo {
            metric: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            sqrt_metric: 1.0,
            symbols: [[[0.0; 3]; 3]; 3],
            dsymbols: [[[[0.0; 3]; 3]; 3]; 3],
        };
        match self.kind {
            CoordinateKind::Cartesian => {}
            CoordinateKind::Polar => {
                let (r, t) = (x, z);
                ci.sqrt_metric = polar_sqrt_metric(x, y, z, self.dim);
                if r != 0.0 {
                    ci.metric[1][1] = 1.0 / (r * r * (t.cos() * t.cos()));
                    if self.dim == 3 {
                        ci.metric[2][2] = 1.0 / (r * r);
                    }
                }
                let s = &mut ci.symbols;
                s[1][1][0] = -r * (t.cos() * t.cos());
                if r != 0.0 {
                    s[0][1][1] = 1.0 / r;
                    s[1][0][1] = 1.0 / r;
                }
                if self.dim == 3 {
                    s[2][2][0] = -r;
                    s[1][1][2] = t.sin() * t.cos();
                    s[1][2][1] = -t.tan();
                    s[2][1][1] = -t.tan();
                    if r != 0.0 {
                        s[2][0][2] = 1.0 / r;
                        s[0][2][2] = 1.0 / r;
                    }
                }
                let d = &mut ci.dsymbols;
                d[1][1][0][0] = -(t.cos() * t.cos());
                if r != 0.0 {
                    d[0][1][1][0] = -1.0 / (r * r);
                    d[1][0][1][0] = -1.0 / (r * r);
                }
                if self.dim == 3 {
                    d[1][1][0][2] = -2.0 * r * t.sin() * t.cos();
                    d[2][2][0][0] = -1.0;
                    d[1][1][2][2] = t.cos() * t.cos() - t.sin() * t.sin();
                    d[1][2][1][2] = -1.0 / (t.cos() * t.cos());
                    d[2][1][1][2] = -1.0 / (t.cos() * t.cos());
                    if r != 0.0 {
                        d[0][2][2][0] = -1.0 / (r * r);
                        d[2][0][2][0] = -1.0 / (r * r);
                    }
                }
            }
            _ => {
                let r = x;
                ci.sqrt_metric = r;
                if r != 0.0 {
                    ci.metric[2][2] = 1.0 / (r * r);
                }
                ci.symbols[2][2][0] = -r;
                if r != 0.0 {
                    ci.symbols[0][2][2] = 1.0 / r;
                    ci.symbols[2][0][2] = 1.0 / r;
                }
                ci.dsymbols[2][2][0][0] = -1.0;
                if r != 0.0 {
                    ci.dsymbols[0][2][2][0] = -1.0 / (r * r);
                    ci.dsymbols[2][0][2][0] = -1.0 / (r * r);
                }
            }
        }
        ci
    }
}

/// `PolarSqrtMetric( r, p, t )`.
fn polar_sqrt_metric(r: f64, _p: f64, t: f64, dim: usize) -> f64 {
    if dim == 2 {
        (r * r * (t.cos() * t.cos())).sqrt()
    } else {
        (crate::element_type::powi(r, 4) * (t.cos() * t.cos())).sqrt()
    }
}

/// The outputs of `CoordinateSystemInfo`; indices are 0-based versions of
/// the Fortran `Symbols(i,j,k)` and `dSymbols(i,j,k,l)`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CoordinateInfo {
    pub metric: [[f64; 3]; 3],
    pub sqrt_metric: f64,
    pub symbols: [[[f64; 3]; 3]; 3],
    pub dsymbols: [[[[f64; 3]; 3]; 3]; 3],
}
