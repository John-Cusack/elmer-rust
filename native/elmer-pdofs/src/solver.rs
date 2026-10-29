//! `elmer-solver` wiring: [`PSetup`] as the framework's
//! [`SolverDofs`], [`PDofsBackend`] to install in `Backends::dofs`, and the
//! PoissonPFEM test plugin as a [`Solver`].
//!
//! Install with
//! `sim.backends.dofs = Some(Box::new(PDofsBackend::new()))` (plus
//! [`PDofsBackend::with_matc`] for `p:%function` degrees); every solver
//! whose section has a non-nodal `Element` keyword then gets its Perm,
//! pattern, element dofs and edge Dirichlet rows from here, and its kernel
//! finds the basis through `SolverCtx::solver_dofs`.

use std::sync::Arc;

use elmer_assembly::initialize_matrix;
use elmer_assembly::ListMatrix;
use elmer_core::{
    BasisEval, CsrMatrix, ElemId, ElementBasis, LocalSystem, MeshView, ParamMap, Params, QuadPoint,
    NO_DOF,
};
use elmer_solver::{
    DofsBackend, Model, Perm, Solver, SolverCtx, SolverDofs, SolverError, SolverRegistry,
};

use crate::setup::{is_p_section, ModelFacts, PSetup};
use crate::Error;

fn unsupported(e: Error) -> SolverError {
    SolverError::Unsupported(e.to_string())
}

/// The block pattern of `a` (`from` components per slot) expanded to
/// `dofs` components.
fn reexpand(a: &CsrMatrix<f64>, from: usize, dofs: usize) -> CsrMatrix<f64> {
    let n = a.nrows / from;
    let mut list = ListMatrix::new(n);
    for i in 0..n {
        let r = from * i;
        let cols: Vec<u32> = a.col_idx[a.row_ptr[r]..a.row_ptr[r + 1]]
            .iter()
            .filter(|&&c| (c as usize).is_multiple_of(from))
            .map(|&c| c / from as u32)
            .collect();
        list.add_sorted(i, &cols);
    }
    initialize_matrix(&list, None, dofs)
}

impl SolverDofs for PSetup {
    fn perm(&self) -> Perm {
        let p = PSetup::perm(self);
        Perm {
            perm: p.perm.clone(),
            n_slots: p.n_dofs,
        }
    }

    fn create_matrix(&self, dofs: usize) -> CsrMatrix<f64> {
        if dofs == self.dofs {
            self.matrix().clone()
        } else {
            reexpand(self.matrix(), self.dofs, dofs)
        }
    }

    fn element_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        PSetup::element_slots(self, e, out);
    }

    fn boundary_extra_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        PSetup::boundary_extra_slots(self, e, out);
    }

    fn element_basis(&self, e: ElemId) -> Result<Box<dyn ElementBasis + Send + Sync>, String> {
        PSetup::element_basis(self, e)
            .map(|b| Box::new(b) as Box<dyn ElementBasis + Send + Sync>)
            .map_err(|err| err.to_string())
    }

    fn element_rule(&self, e: ElemId) -> Result<Vec<QuadPoint>, String> {
        self.gauss_rule(e)
            .map(|r| r.points())
            .map_err(|err| err.to_string())
    }

    fn element_nodes(&self, e: ElemId) -> Vec<u32> {
        PSetup::element_nodes(self, e).to_vec()
    }

    fn local_bubbles(&self, e: ElemId) -> usize {
        PSetup::local_bubbles(self, e)
    }
}

/// MATC evaluation of `function(body, x, y, z)` for `p:%function`.
pub type MatcFn = dyn Fn(&str, [f64; 4]) -> Result<f64, String> + Send + Sync;

/// The [`DofsBackend`] serving every solver section with a non-nodal
/// `Element` keyword from [`PSetup`].
///
/// Refused: `Element` in equation sections or body-wise
/// `Solver i: Element`, and boundary elements that a BC `Body Id` makes
/// active (both change `Def_Dofs` or the active set in ways not ported).
#[derive(Clone, Default)]
pub struct PDofsBackend {
    matc: Option<Arc<MatcFn>>,
}

impl PDofsBackend {
    pub fn new() -> Self {
        Self::default()
    }

    /// Evaluates `p:%function` degrees with `matc`.
    pub fn with_matc(mut self, matc: Arc<MatcFn>) -> Self {
        self.matc = Some(matc);
        self
    }
}

impl DofsBackend for PDofsBackend {
    fn for_solver(
        &self,
        mesh: &dyn MeshView,
        model: &Model,
        solver: &ParamMap,
        active: &[ElemId],
        dofs: usize,
    ) -> Option<Result<Box<dyn SolverDofs>, SolverError>> {
        if !is_p_section(solver) {
            return None;
        }
        Some(self.build(mesh, model, solver, active, dofs))
    }
}

impl PDofsBackend {
    fn build(
        &self,
        mesh: &dyn MeshView,
        model: &Model,
        solver: &ParamMap,
        active: &[ElemId],
        dofs: usize,
    ) -> Result<Box<dyn SolverDofs>, SolverError> {
        if let Some(i) = model
            .equations
            .iter()
            .position(|q| q.get("element").is_some())
        {
            return Err(SolverError::Unsupported(format!(
                "Equation {}: Element (element definitions in equation sections)",
                i + 1
            )));
        }
        if model.bodies.iter().any(|b| {
            b.0.keys()
                .any(|k| k.starts_with("solver ") && k.ends_with(": element"))
        }) {
            return Err(SolverError::Unsupported(
                "Body: Solver i: Element (body-wise element definitions)".into(),
            ));
        }
        if model
            .boundary_conditions
            .iter()
            .any(|bc| bc.get("body id").is_some())
        {
            return Err(SolverError::Unsupported(
                "Boundary Condition: Body Id with p elements".into(),
            ));
        }
        let mut is_active = vec![false; mesh.n_bulk() + mesh.n_boundary()];
        for &e in active {
            is_active[e.0 as usize] = true;
        }
        let body_has_equation = |b: usize| {
            model
                .bodies
                .get(b - 1)
                .is_some_and(|p| p.get("equation").is_some())
        };
        let matc = self.matc.clone();
        let mut matc_degree = move |f: &str, cx: [f64; 4]| -> Result<f64, Error> {
            match &matc {
                Some(m) => m(f, cx).map_err(Error::BadDefinition),
                None => Err(Error::Unsupported(format!(
                    "MATC degree function {f}: no MATC engine"
                ))),
            }
        };
        let setup = PSetup::new(
            mesh,
            solver,
            &model.simulation,
            ModelFacts {
                n_bodies: model.bodies.len(),
                body_has_equation: &body_has_equation,
                active: &|e: ElemId| is_active[e.0 as usize],
                boundary_body: &|_| 0,
                matc_degree: &mut matc_degree,
            },
            dofs,
        )
        .map_err(unsupported)?;
        Ok(Box::new(setup))
    }
}

/// `fem/tests/PoissonPFEM/Poisson.f90`, the case's own plugin
/// (`Procedure = "./Poisson" "PoissonSolver"`): `-div grad u = Source` with
/// the source interpolated by the nodal basis only, bubbles condensed with
/// `LCondensate` when they are not global. A test-local fixture (strategy
/// D12: it lets the case run, it is not product code).
///
/// Refuses a solver without p dofs: the plugin is only ported for p
/// elements.
pub struct PoissonPfem;

/// [`PoissonPfem`] as a solver factory.
pub fn poisson_pfem(_: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(PoissonPfem))
}

/// Registers [`PoissonPfem`] under the name PoissonPFEM's `case.sif` uses.
/// Other cases ship different `Poisson.f90` plugins under the same name;
/// register it only for this case.
pub fn register_poisson_pfem(reg: &mut SolverRegistry) {
    reg.register("Poisson", "PoissonSolver", poisson_pfem);
}

/// `InvertMatrix` of a small dense matrix (Gauss-Jordan with partial
/// pivoting), row-major `n * n`.
fn invert(a: &mut [f64], n: usize) -> Result<(), SolverError> {
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for c in 0..n {
        let p = (c..n)
            .max_by(|&i, &j| a[i * n + c].abs().total_cmp(&a[j * n + c].abs()))
            .expect("non-empty");
        if a[p * n + c] == 0.0 {
            return Err(SolverError::fatal("InvertMatrix", "Matrix is singular"));
        }
        for k in 0..n {
            a.swap(c * n + k, p * n + k);
            inv.swap(c * n + k, p * n + k);
        }
        let d = a[c * n + c];
        for k in 0..n {
            a[c * n + k] /= d;
            inv[c * n + k] /= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r * n + c];
                if f != 0.0 {
                    for k in 0..n {
                        a[r * n + k] -= f * a[c * n + k];
                        inv[r * n + k] -= f * inv[c * n + k];
                    }
                }
            }
        }
    }
    a.copy_from_slice(&inv);
    Ok(())
}

/// `LCondensate(N, Nb, K, F)`: eliminates the last `nb` of `n + nb`
/// local dofs (`K` row-major, size `m = n + nb`), leaving `K(1:n,1:n)` and
/// `F(1:n)`.
fn condensate(
    n: usize,
    nb: usize,
    k: &[f64],
    f: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), SolverError> {
    let m = n + nb;
    let mut kbb: Vec<f64> = (0..nb * nb)
        .map(|ij| k[(n + ij / nb) * m + n + ij % nb])
        .collect();
    invert(&mut kbb, nb)?;
    // X = Kbb^-1 Kbl, y = Kbb^-1 Fb
    let mut x = vec![0.0; nb * n];
    let mut y = vec![0.0; nb];
    for i in 0..nb {
        for l in 0..nb {
            let s = kbb[i * nb + l];
            for j in 0..n {
                x[i * n + j] += s * k[(n + l) * m + j];
            }
            y[i] += s * f[n + l];
        }
    }
    let mut ks = vec![0.0; n * n];
    let mut fs = vec![0.0; n];
    for i in 0..n {
        fs[i] = f[i];
        for j in 0..n {
            ks[i * n + j] = k[i * m + j];
        }
        for l in 0..nb {
            let klb = k[i * m + n + l];
            fs[i] -= klb * y[l];
            for j in 0..n {
                ks[i * n + j] -= klb * x[l * n + j];
            }
        }
    }
    Ok((ks, fs))
}

impl Solver for PoissonPfem {
    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.default_initialize()?;
        if ctx.solver_dofs().is_none() {
            return Err(SolverError::Unsupported(
                "PoissonPFEM Poisson plugin without p elements (no DofsBackend)".into(),
            ));
        }
        let active: Vec<ElemId> = ctx.active_elements().collect();
        let mut slots = Vec::new();
        let mut ev = BasisEval::default();
        for e in active {
            let d = ctx.solver_dofs().expect("checked above");
            d.element_slots(e, &mut slots);
            let nd = slots.len();
            let nb = d.local_bubbles(e);
            let nodes = d.element_nodes(e);
            let basis = d
                .element_basis(e)
                .map_err(|m| SolverError::fatal("ElementInfo", m))?;
            let rule = d
                .element_rule(e)
                .map_err(|m| SolverError::fatal("GaussPoints", m))?;
            let n = nodes.len();
            let etype = ctx.mesh.element(e).etype;
            let xyz: Vec<[f64; 3]> = nodes.iter().map(|&k| ctx.mesh.coord(k)).collect();
            // Load(1:n) = GetReal(BodyForce, 'Source')
            let load = match ctx.body_params(e).and_then(|b| b.body_force) {
                Some(bf) => ctx
                    .real_at_nodes(bf, "source", e)?
                    .unwrap_or_else(|| vec![0.0; n]),
                None => vec![0.0; n],
            };
            // LocalMatrix(STIFF, FORCE, LOAD, Element, n, nd+nb)
            let m = nd + nb;
            let mut stiff = vec![0.0; m * m];
            let mut force = vec![0.0; m];
            for qp in &rule {
                basis
                    .eval(etype, &xyz, qp.uvw, &mut ev)
                    .map_err(|msg| SolverError::fatal("ElementInfo", msg))?;
                if ev.basis.len() < m {
                    return Err(SolverError::fatal(
                        "PoissonSolver",
                        format!(
                            "element {} has {} basis functions, {m} dofs",
                            e.0,
                            ev.basis.len()
                        ),
                    ));
                }
                let w = qp.weight * ev.det_j;
                let load_at_ip: f64 = ev.basis[..n].iter().zip(&load).map(|(b, l)| b * l).sum();
                for i in 0..m {
                    let gi = ev.dbasis_dx[i];
                    for j in 0..m {
                        let gj = ev.dbasis_dx[j];
                        stiff[i * m + j] += w * (gi[0] * gj[0] + gi[1] * gj[1] + gi[2] * gj[2]);
                    }
                    force[i] += w * load_at_ip * ev.basis[i];
                }
            }
            let mut local = LocalSystem::zeros(nd);
            if nb > 0 {
                let (k, f) = condensate(nd, nb, &stiff, &force)?;
                local.stiff = k;
                local.force = f;
            } else {
                local.stiff = stiff;
                local.force = force;
            }
            debug_assert!(slots.iter().all(|&s| s != NO_DOF));
            ctx.default_update(&local, e)?;
        }
        ctx.default_dirichlet()?;
        ctx.default_solve()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pattern_reexpands_to_other_component_counts() {
        let mut list = ListMatrix::new(3);
        list.add_clique(&mut [0, 2]);
        list.add_clique(&mut [1, 2]);
        let one = initialize_matrix(&list, None, 1);
        let three = initialize_matrix(&list, None, 3);
        let r = reexpand(&three, 3, 1);
        assert_eq!(
            (r.row_ptr.clone(), r.col_idx.clone()),
            (one.row_ptr.clone(), one.col_idx.clone())
        );
        let r = reexpand(&one, 1, 3);
        assert_eq!((r.row_ptr, r.col_idx), (three.row_ptr, three.col_idx));
    }

    #[test]
    fn condensation_matches_the_schur_complement() {
        // K = [[4, 1, 2], [1, 3, 0], [2, 0, 5]], F = [1, 2, 3]; condense the last dof.
        let k = [4.0, 1.0, 2.0, 1.0, 3.0, 0.0, 2.0, 0.0, 5.0];
        let f = [1.0, 2.0, 3.0];
        let (ks, fs) = condensate(2, 1, &k, &f).unwrap();
        let close = |a: f64, b: f64| (a - b).abs() < 1e-15;
        assert!(
            close(ks[0], 4.0 - 4.0 / 5.0)
                && close(ks[1], 1.0)
                && close(ks[2], 1.0)
                && close(ks[3], 3.0)
        );
        assert!(close(fs[0], 1.0 - 6.0 / 5.0) && close(fs[1], 2.0));
    }
}
