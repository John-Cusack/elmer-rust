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

fn has_p_or_b(def: &str) -> bool {
    def.contains("p:") || def.contains("b:")
}

/// `SetGlobalBubblesFlag` (MainUtils.F90:1106) for a solver-section
/// definition: `Bubbles in Global System`, else whether `Element`
/// contains `p:` or `b:`.
pub fn global_bubbles(solver: &dyn Params) -> bool {
    solver
        .logical("bubbles in global system")
        .unwrap_or_else(|| solver.string("element").is_some_and(has_p_or_b))
}

/// `SetGlobalBubblesFlag` with the model fallbacks: `Bubbles in Global
/// System`; else the solver `Element`; else the `Element` of each
/// equation whose `Active Solvers` include the solver; and, when the last
/// of those lookups found nothing (the Fortran tests the `Found` of the
/// last `ListGet` call), the bodies' `Solver id: Element`. Any `p:` or `b:`
/// found turns the global bubbles on.
pub fn global_bubbles_with(solver: &dyn Params, facts: &ModelFacts<'_>) -> bool {
    if let Some(gb) = solver.logical("bubbles in global system") {
        return gb;
    }
    if let Some(def) = solver.string("element") {
        return has_p_or_b(def);
    }
    let mut gb = false;
    let mut found = false;
    for q in 1..=facts.n_equations {
        match (facts.equation_activates)(q) {
            None => found = false,
            Some(active) => {
                found = true;
                if active {
                    let def = (facts.equation_element)(q);
                    found = def.is_some();
                    if def.is_some_and(|d| has_p_or_b(&d)) {
                        gb = true;
                    }
                }
            }
        }
    }
    if !found {
        for b in 1..=facts.n_bodies {
            if (facts.body_solver_element)(b).is_some_and(|d| has_p_or_b(&d)) {
                gb = true;
            }
        }
    }
    gb
}

/// Another solver on the same mesh: Elmer sizes the mesh's p tables
/// (`inDofs`, `NonNodalElements`) from the definitions of every solver.
pub struct OtherSolver<'a> {
    /// 1-based solver id (`Solver % SolverId`).
    pub id: usize,
    pub section: &'a dyn Params,
    /// `Element{id}` of equation `q`.
    pub equation_element: &'a dyn Fn(usize) -> Option<String>,
    /// `Solver id: Element` of body `b`.
    pub body_element: &'a dyn Fn(usize) -> Option<String>,
}

/// The model facts [`PSetup::new`] needs besides the keywords.
pub struct ModelFacts<'a> {
    /// Number of bodies.
    pub n_bodies: usize,
    /// Number of equation sections.
    pub n_equations: usize,
    /// The `Equation` of body `b` (1-based), when it has one.
    pub body_equation: &'a dyn Fn(usize) -> Option<usize>,
    /// `Element` of equation `q` (1-based).
    pub equation_element: &'a dyn Fn(usize) -> Option<String>,
    /// `Element{id}` of equation `q` for this solver's id.
    pub equation_solver_element: &'a dyn Fn(usize) -> Option<String>,
    /// `Solver id: Element` of body `b` for this solver's id.
    pub body_solver_element: &'a dyn Fn(usize) -> Option<String>,
    /// Whether equation `q`'s `Active Solvers` include this solver, `None`
    /// when it has no `Active Solvers` (`SetGlobalBubblesFlag`'s fallback).
    pub equation_activates: &'a dyn Fn(usize) -> Option<bool>,
    /// `CheckElementEquation` for the solver, over bulk and boundary
    /// elements (boundary elements count only with a BC `Body Id`).
    pub active: &'a dyn Fn(ElemId) -> bool,
    /// `BodyId` of boundary element `i` (0-based among boundary elements),
    /// 0 unless the BC gives one.
    pub boundary_body: &'a dyn Fn(usize) -> i32,
    /// MATC evaluation for `p:%function` degrees.
    pub matc_degree: &'a mut MatcDegree<'a>,
    /// This solver's 1-based id (0 when unknown), which orders it among
    /// [`ModelFacts::other_solvers`].
    pub solver_id: usize,
    /// The other solvers on the mesh (those without their own `Mesh`); empty
    /// sizes the mesh from this solver alone.
    pub other_solvers: &'a [OtherSolver<'a>],
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
    /// The setup of one solver with `dofs` components on `mesh`, from its
    /// own `Element`, the equations' `Element` / `Element{id}` and the
    /// bodies' `Solver id: Element` with Elmer's precedence
    /// (`NonNodalElements`). The mesh tables (edge and face DOF counts,
    /// bubbles, Gauss points) take the maximum over this solver and
    /// [`ModelFacts::other_solvers`], as Elmer's `inDofs`; the Perm and the
    /// basis use this solver's own definitions.
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
        let def = solver.string("element");
        let need_edges = match mesh.dim() {
            2 => simulation.logical("need edges 2d"),
            3 => simulation.logical("need edges 3d"),
            _ => None,
        };
        let gb = global_bubbles_with(solver, &facts);
        // LoadModel's Def_Dofs of every solver, in solver order, and the
        // maximum over them; this solver is the one at `me`.
        struct Entry<'s> {
            id: usize,
            me: bool,
            element: Option<&'s str>,
            dg: bool,
            body_element: &'s dyn Fn(usize) -> Option<String>,
            equation_element: &'s dyn Fn(usize) -> Option<String>,
        }
        let mut entries = vec![Entry {
            id: facts.solver_id,
            me: true,
            element: def,
            dg: false,
            body_element: facts.body_solver_element,
            equation_element: facts.equation_solver_element,
        }];
        for o in facts.other_solvers {
            entries.push(Entry {
                id: o.id,
                me: false,
                element: o.section.string("element"),
                dg: o.section.logical("discontinuous galerkin") == Some(true),
                body_element: o.body_element,
                equation_element: o.equation_element,
            });
        }
        entries.sort_by_key(|e| e.id);
        let me = entries
            .iter()
            .position(|e| e.me)
            .expect("this solver is listed");
        let mut max = initial_max_defs();
        let mut all_defs = entries
            .iter()
            .map(|e| ElementDefs::from_solver_element(e.element, facts.n_bodies, e.dg, &mut max))
            .collect::<Result<Vec<_>, _>>()?;
        let pm = {
            let mut solvers: Vec<SolverElements<'_>> = entries
                .iter()
                .zip(all_defs.iter_mut())
                .map(|(e, defs)| SolverElements {
                    defs,
                    element: e.element,
                    dg: e.dg,
                    body_element: e.body_element,
                    equation_element: e.equation_element,
                })
                .collect();
            PMesh::new(
                mesh,
                &mut solvers,
                &mut max,
                PMeshOptions {
                    serendipity: simulation.logical("serendipity p elements").unwrap_or(true),
                    need_edges,
                    n_bodies: facts.n_bodies,
                    n_equations: facts.n_equations,
                    body_equation: facts.body_equation,
                    equation_element: facts.equation_element,
                    matc_degree: facts.matc_degree,
                    boundary_body: facts.boundary_body,
                },
            )?
        };
        let mut pm = pm;
        let defs = all_defs.swap_remove(me);
        let bubbles = solver.logical("bubbles").unwrap_or(false);
        let opts = PermOptions {
            // `No Matrix`: `CreateMatrix` returns before `OptimizeBandwidth`
            // (ElementUtils.F90:2229).
            optimize_bandwidth: solver.logical("optimize bandwidth").unwrap_or(true)
                && solver.logical("no matrix") != Some(true),
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
