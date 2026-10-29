//! `PoissonAssembly.f90` of the `CoupledPoisson*`, `BlockPoisson1` and
//! `BlockPoisson2` cases: the element routines a solver section names with
//! `Bulk Assembly Procedure = File "./PoissonAssembly" "BulkAssembly"` and
//! `Boundary Assembly Procedure ... "BoundaryAssembly"`.
//!
//! The file is the same in all nine `CoupledPoisson` cases up to
//! declaration order and `USE` lines (case 3 adds unused `USE` lines). What
//! differs between the cases is the solver around it: coupled variables,
//! integral, equality and floating constraints, line search. None of that is
//! in this file.
//!
//! The kernels take what the Fortran fetches through `DefUtils`: nodal
//! property values of the current element (`GetReal` at the element's nodes;
//! zeros where the keyword is absent) and the element coordinates, plus the
//! quadrature and basis services of [`elmer_core`]. They add into a
//! [`LocalSystem`] exactly as the Fortran adds into `Stiff` and `Force`.

use elmer_core::{BasisEval, ElementBasis, ElementType, LocalSystem, Logger, Quadrature};

/// Elmer's default `Info` level (`Messages.F90`, `DefLevel`).
const INFO_LEVEL: u8 = 4;

/// Why a kernel call failed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AssemblyError {
    #[error("{what}: {detail}")]
    Shape { what: &'static str, detail: String },
    #[error("ElementInfo failed for a {etype:?} element: {detail}")]
    Basis { etype: ElementType, detail: String },
}

/// One element: its type and the coordinates of its nodes in Elmer order
/// (`GetElementNodes`).
#[derive(Clone, Copy, Debug)]
pub struct Element<'a> {
    pub etype: ElementType,
    pub nodes: &'a [[f64; 3]],
}

/// Nodal property values on the element. `None` is a keyword the section does
/// not define, which `GetReal` answers with zeros.
#[derive(Clone, Copy, Debug, Default)]
pub struct BulkProps<'a> {
    /// Material `Conductivity`.
    pub conductivity: Option<&'a [f64]>,
    /// Material `Source`.
    pub source: Option<&'a [f64]>,
    /// Material `Reaction`.
    pub reaction: Option<&'a [f64]>,
}

/// Boundary condition values on the element's nodes.
#[derive(Clone, Copy, Debug, Default)]
pub struct BoundaryProps<'a> {
    /// `Flux`.
    pub flux: Option<&'a [f64]>,
    /// `Coefficient`.
    pub coefficient: Option<&'a [f64]>,
    /// `External`.
    pub external: Option<&'a [f64]>,
}

/// Element kernels with the per-instance state of the Fortran's `SAVE`d
/// `Visited` flags: each routine reports `1st time` once per instance.
#[derive(Debug, Default)]
pub struct PoissonAssembly {
    bulk_visited: bool,
    boundary_visited: bool,
}

impl PoissonAssembly {
    pub fn new() -> Self {
        Self::default()
    }

    /// `SUBROUTINE BulkAssembly`: adds
    /// `K(p,q) += w (cond grad(p).grad(q) + reaction N(p) N(q))` and
    /// `F(p) += w source N(p)` over the element's quadrature points, with
    /// every property interpolated to the point.
    ///
    /// `local` must be sized for the element (`n` equals the node count) and
    /// is accumulated into, not overwritten.
    pub fn bulk(
        &mut self,
        el: Element<'_>,
        props: BulkProps<'_>,
        basis: &dyn ElementBasis,
        quad: &dyn Quadrature,
        log: &dyn Logger,
        local: &mut LocalSystem,
    ) -> Result<(), AssemblyError> {
        if !self.bulk_visited {
            log.info("PoissonBulkAssembly", "1st time", INFO_LEVEL);
            self.bulk_visited = true;
        }
        let n = check_shape(&el, local, &[
            ("Conductivity", props.conductivity),
            ("Source", props.source),
            ("Reaction", props.reaction),
        ])?;
        integrate(el, n, basis, quad, local, |pt| {
            let cond = interp(pt.ev, props.conductivity);
            let source = interp(pt.ev, props.source);
            let reaction = interp(pt.ev, props.reaction);
            pt.add_gradient(cond);
            pt.add_mass(reaction);
            pt.add_load(source);
        })
    }

    /// `SUBROUTINE BoundaryAssembly`: adds
    /// `K(p,q) += w coeff N(p) N(q)` and
    /// `F(p) += w N(p) (flux + coeff external)`.
    ///
    /// Returns `false` and adds nothing when the boundary condition defines
    /// neither `Flux` nor `Coefficient` (the Fortran `RETURN`s); `External`
    /// alone does not count.
    pub fn boundary(
        &mut self,
        el: Element<'_>,
        props: BoundaryProps<'_>,
        basis: &dyn ElementBasis,
        quad: &dyn Quadrature,
        log: &dyn Logger,
        local: &mut LocalSystem,
    ) -> Result<bool, AssemblyError> {
        if !self.boundary_visited {
            log.info("PoissonBoundaryAssembly", "1st time", INFO_LEVEL);
            self.boundary_visited = true;
        }
        let n = check_shape(&el, local, &[
            ("Flux", props.flux),
            ("Coefficient", props.coefficient),
            ("External", props.external),
        ])?;
        if props.flux.is_none() && props.coefficient.is_none() {
            return Ok(false);
        }
        integrate(el, n, basis, quad, local, |pt| {
            let flux = interp(pt.ev, props.flux);
            let coeff = interp(pt.ev, props.coefficient);
            let ext = interp(pt.ev, props.external);
            pt.add_mass(coeff);
            pt.add_load_basis_first(flux + coeff * ext);
        })?;
        Ok(true)
    }
}

/// One quadrature point: its weight (`IntegStuff % s(t) * detJ`), the basis
/// there and the system being accumulated into.
struct Point<'a> {
    weight: f64,
    n: usize,
    ev: &'a BasisEval,
    local: &'a mut LocalSystem,
}

impl Point<'_> {
    /// `STIFF(p,q) += W * c * SUM(dBasisdx(p,1:3)*dBasisdx(q,1:3))`. The
    /// gradient has three components in every dimension.
    fn add_gradient(&mut self, c: f64) {
        let (n, w, g) = (self.n, self.weight, &self.ev.dbasis_dx);
        for p in 0..n {
            for q in 0..n {
                let dot = g[p][0] * g[q][0] + g[p][1] * g[q][1] + g[p][2] * g[q][2];
                self.local.stiff[p * n + q] += w * c * dot;
            }
        }
    }

    /// `STIFF(p,q) += W * c * Basis(p)*Basis(q)`.
    fn add_mass(&mut self, c: f64) {
        let (n, w, b) = (self.n, self.weight, &self.ev.basis);
        for p in 0..n {
            for q in 0..n {
                self.local.stiff[p * n + q] += w * c * b[p] * b[q];
            }
        }
    }

    /// `FORCE(p) += W * rhs * Basis(p)` (bulk association order).
    fn add_load(&mut self, rhs: f64) {
        let (w, b) = (self.weight, &self.ev.basis);
        for (f, bp) in self.local.force.iter_mut().zip(b) {
            *f += w * rhs * bp;
        }
    }

    /// `FORCE(p) += W * Basis(p) * rhs`: the boundary routine multiplies in
    /// this order, which rounds differently from [`add_load`](Self::add_load).
    fn add_load_basis_first(&mut self, rhs: f64) {
        let (w, b) = (self.weight, &self.ev.basis);
        for (f, bp) in self.local.force.iter_mut().zip(b) {
            *f += w * bp * rhs;
        }
    }
}

/// `SUM(Basis(1:Nrow) * Values(1:Nrow))`, zero for an absent keyword.
fn interp(ev: &BasisEval, values: Option<&[f64]>) -> f64 {
    values.map_or(0.0, |v| v.iter().zip(&ev.basis).map(|(a, b)| a * b).sum())
}

fn check_shape(
    el: &Element<'_>,
    local: &LocalSystem,
    props: &[(&'static str, Option<&[f64]>)],
) -> Result<usize, AssemblyError> {
    let n = el.nodes.len();
    if local.n != n || local.stiff.len() != n * n || local.force.len() != n {
        return Err(AssemblyError::Shape {
            what: "LocalSystem",
            detail: format!("sized for {} dofs, element has {n} nodes", local.n),
        });
    }
    for (name, v) in props {
        if let Some(v) = v {
            if v.len() != n {
                return Err(AssemblyError::Shape {
                    what: name,
                    detail: format!("{} values for {n} nodes", v.len()),
                });
            }
        }
    }
    Ok(n)
}

/// The `DO t=1,IntegStuff % n` loop: evaluates the basis at each point and
/// hands the weighted point to `f`.
fn integrate(
    el: Element<'_>,
    n: usize,
    basis: &dyn ElementBasis,
    quad: &dyn Quadrature,
    local: &mut LocalSystem,
    mut f: impl FnMut(&mut Point<'_>),
) -> Result<(), AssemblyError> {
    let mut ev = BasisEval::default();
    for qp in quad.gauss_points(el.etype, None) {
        basis
            .eval(el.etype, el.nodes, qp.uvw, &mut ev)
            .map_err(|detail| AssemblyError::Basis { etype: el.etype, detail })?;
        if ev.basis.len() < n || ev.dbasis_dx.len() < n {
            return Err(AssemblyError::Shape {
                what: "ElementInfo",
                detail: format!("{} basis values for {n} nodes", ev.basis.len()),
            });
        }
        f(&mut Point { weight: qp.weight * ev.det_j, n, ev: &ev, local });
    }
    Ok(())
}
