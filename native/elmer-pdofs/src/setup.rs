//! One solver's whole p-element setup from its keywords: what
//! `AddEquationBasics` and `CreateMatrix` (MainUtils.F90:1655) do for a
//! solver section with `Element = "p:..."`, packaged for a framework
//! adapter.
//!
//! [`PSetup`] answers the questions a solver framework asks of a nodal
//! one, for every DOF kind: the Perm and its size, the matrix pattern, the
//! rows of an element's DOFs, the extra Dirichlet rows of a boundary
//! element, and the basis and integration rule of an element. Keywords
//! read:
//!
//! * solver: `Element`, `Bubbles in Global System`, `Bubbles`,
//!   `Optimize Bandwidth`, `Optimize Bandwidth Use Always`,
//!   `Discontinuous Galerkin` (refused);
//! * simulation: `Serendipity P Elements`, `Need Edges` (2D).

use elmer_assembly::{Perm, PermOptions};
use elmer_core::{CsrMatrix, ElemId, MeshView, Params, NO_DOF};
use elmer_fe::GaussRule;
use elmer_fe_hp::PElementBasis;

use crate::defs::{initial_max_defs, ElementDefs};
use crate::dofs::DofLayout;
use crate::perm::{create_p_system, PSystem};
use crate::pmesh::{MatcDegree, PMesh, PMeshOptions, SolverElements};
use crate::Error;

/// Whether a solver section needs this crate: an `Element` keyword other
/// than plain nodal `n:1`.
pub fn is_p_section(solver: &dyn Params) -> bool {
    solver
        .string("element")
        .is_some_and(|e| e.trim().to_lowercase() != "n:1")
}

/// `SetGlobalBubblesFlag` (MainUtils.F90:1106) for a solver-section
/// definition: `Bubbles in Global System`, else whether `Element`
/// contains `p:` or `b:`.
pub fn global_bubbles(solver: &dyn Params) -> bool {
    solver
        .logical("bubbles in global system")
        .unwrap_or_else(|| {
            solver
                .string("element")
                .is_some_and(|s| s.contains("p:") || s.contains("b:"))
        })
}

/// The model facts [`PSetup::new`] needs besides the keywords.
pub struct ModelFacts<'a> {
    /// Number of bodies.
    pub n_bodies: usize,
    /// Whether body `b` (1-based) has an `Equation`.
    pub body_has_equation: &'a dyn Fn(usize) -> bool,
    /// `CheckElementEquation` for the solver, over bulk and boundary
    /// elements (boundary elements count only with a BC `Body Id`).
    pub active: &'a dyn Fn(ElemId) -> bool,
    /// `BodyId` of boundary element `i` (0-based among boundary elements),
    /// 0 unless the BC gives one.
    pub boundary_body: &'a dyn Fn(usize) -> i32,
    /// MATC evaluation for `p:%function` degrees.
    pub matc_degree: &'a mut MatcDegree<'a>,
}

/// A solver's p mesh, element definitions and system pattern.
#[derive(Clone, Debug)]
pub struct PSetup {
    pub mesh: PMesh,
    pub defs: ElementDefs,
    pub global_bubbles: bool,
    pub bubbles: bool,
    pub system: PSystem,
    /// `Variable DOFs`.
    pub dofs: usize,
}

impl PSetup {
    /// The setup of one solver with `dofs` components on `mesh`. Only the
    /// solver-section `Element` is read; the caller must refuse element
    /// definitions in equation sections and body-wise `Solver i: Element`
    /// keywords.
    pub fn new(
        mesh: &dyn MeshView,
        solver: &dyn Params,
        simulation: &dyn Params,
        facts: ModelFacts<'_>,
        dofs: usize,
    ) -> Result<Self, Error> {
        if solver.logical("discontinuous galerkin") == Some(true) {
            return Err(Error::Unsupported("Discontinuous Galerkin".into()));
        }
        let def = solver
            .string("element")
            .ok_or_else(|| Error::BadDefinition("no Element keyword".into()))?;
        let mut max = initial_max_defs();
        let mut defs =
            ElementDefs::from_solver_element(Some(def), facts.n_bodies, false, &mut max)?;
        let need_edges = match mesh.dim() {
            2 => simulation.logical("need edges 2d"),
            3 => simulation.logical("need edges 3d"),
            _ => None,
        };
        let none = |_: usize| None;
        let mut pm = {
            let mut solvers = [SolverElements {
                defs: &mut defs,
                element: Some(def),
                dg: false,
                body_element: &none,
            }];
            PMesh::new(
                mesh,
                &mut solvers,
                &mut max,
                PMeshOptions {
                    serendipity: simulation.logical("serendipity p elements").unwrap_or(true),
                    need_edges,
                    n_bodies: facts.n_bodies,
                    body_has_equation: facts.body_has_equation,
                    matc_degree: facts.matc_degree,
                    boundary_body: facts.boundary_body,
                },
            )?
        };
        let gb = global_bubbles(solver);
        let bubbles = solver.logical("bubbles").unwrap_or(false);
        let opts = PermOptions {
            optimize_bandwidth: solver.logical("optimize bandwidth").unwrap_or(true),
            use_optimized_always: solver
                .logical("optimize bandwidth use always")
                .unwrap_or(false),
        };
        let active = |e: usize| (facts.active)(ElemId(e as u32));
        let system = create_p_system(&mut pm, &defs, gb, bubbles, &active, opts, dofs);
        Ok(Self {
            mesh: pm,
            defs,
            global_bubbles: gb,
            bubbles,
            system,
            dofs,
        })
    }

    pub fn layout(&self) -> DofLayout<'_> {
        DofLayout {
            mesh: &self.mesh,
            defs: &self.defs,
            global_bubbles: self.global_bubbles,
            bubbles: self.bubbles,
        }
    }

    /// Slot to Perm index (`Ndeg` long, [`NO_DOF`] where unused).
    pub fn perm(&self) -> &Perm {
        &self.system.perm
    }

    /// The zero-valued pattern, `perm().n_dofs * dofs` square.
    pub fn matrix(&self) -> &CsrMatrix<f64> {
        &self.system.matrix
    }

    /// Perm indices of element `e`'s DOFs (`Perm(Indexes(1:nd))`), the
    /// `slots` the glue takes; [`NO_DOF`] where absent.
    pub fn element_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        self.layout()
            .element_rows(e.0 as usize, &self.system.perm.perm, out);
    }

    /// `GetElementNOFBDOFs`: bubbles a kernel condenses locally (0 with
    /// global bubbles).
    pub fn local_bubbles(&self, e: ElemId) -> usize {
        let e = e.0 as usize;
        if self.global_bubbles || e >= self.mesh.n_bulk {
            return 0;
        }
        self.layout().local_bubble_count(e)
    }

    /// Perm indices `DefaultDirichletBCs` sets to zero on boundary element
    /// `e` beyond its nodes (its edge DOFs) when a BC with a constant value
    /// gives the field and the parent is an active p element; empty
    /// otherwise.
    pub fn boundary_extra_slots(&self, e: ElemId, out: &mut Vec<u32>) {
        out.clear();
        let lay = self.layout();
        let i = e.0 as usize;
        if i < self.mesh.n_bulk {
            return;
        }
        let el = &self.mesh.elements[i];
        let Some(parent) = el.parents[0].or(el.parents[1]) else {
            return;
        };
        if !lay.is_active_p(parent as usize) {
            return;
        }
        let mut ind = Vec::new();
        lay.boundary_indexes_from_parent(i, &mut ind);
        let perm = &self.system.perm.perm;
        out.extend(
            ind.iter()
                .skip(el.nodes.len())
                .map(|&g| perm.get(g as usize).copied().unwrap_or(NO_DOF))
                .filter(|&p| p != NO_DOF),
        );
    }

    /// Whether `ElementInfo` uses the p basis for element `e`.
    pub fn is_p_element(&self, e: ElemId) -> bool {
        self.layout().is_active_p(e.0 as usize)
    }

    /// The p basis of element `e` (as `ElementInfo` evaluates it).
    pub fn element_basis(&self, e: ElemId) -> Result<PElementBasis, Error> {
        self.layout().element_basis(e.0 as usize)
    }

    /// `GaussPoints(Element)` of element `e`.
    pub fn gauss_rule(&self, e: ElemId) -> Result<GaussRule, Error> {
        self.layout().gauss_rule(e.0 as usize)
    }

    /// `GaussPoints(Element, np)` of element `e`.
    pub fn gauss_rule_np(&self, e: ElemId, np: usize) -> Result<GaussRule, Error> {
        self.layout().gauss_rule_np(e.0 as usize, np)
    }

    /// Element `e`'s nodes in the order the p basis expects (boundary
    /// elements follow their parent edge, which can differ from the mesh).
    pub fn element_nodes(&self, e: ElemId) -> &[u32] {
        &self.mesh.elements[e.0 as usize].nodes
    }
}
