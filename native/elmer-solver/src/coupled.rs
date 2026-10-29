//! Solvers that Elmer runs without a `Procedure` of their own, chosen by the
//! solver mode keywords (MainUtils.F90:1388-1395): `Assembly Solver = True`
//! runs `CoupledSolver` (MainUtils.F90:3565) in its one-equation mode, which
//! assembles the system element by element with the routine named by
//! `Bulk Assembly Procedure` (`BlockSystemAssembly`, MainUtils.F90:4885).
//!
//! `Coupled Solver = True` without a `Procedure` runs `CoupledSolver`'s
//! coupled mode ([`CoupledSystem`]): `Variable 1..n` and `Constraint 1..m`
//! in one monolithic matrix, each block assembled by its own routine.
//!
//! `Block Solver`, and an assembly or coupled solver that also names a
//! `Procedure` (`CoupledSolver`'s "inheriting matrix from procedure" mode),
//! are refused.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use elmer_core::{normalize_key, BasisEval, CsrMatrix, ElemId, FieldKind, LocalSystem, ParamMap, Params, Value, NO_DOF};

use crate::error::SolverError;
use crate::sim::LinearSystem;
use crate::solver::{procedure_value, ElementAssembly, Solver, SolverCtx, SolverRegistry};
use crate::variable::{VarId, Variable};

const CALLER: &str = "CoupledSolver";
const RULE: &str = "---------------------------------------";
const LONG_RULE: &str = "-------------------------------------------------";

/// The built-in solver a section's mode keywords select, if any. Errors on
/// the modes that are not ported and on an assembly routine the registry
/// does not have.
pub fn solver_for_mode(params: &ParamMap, registry: &SolverRegistry) -> Result<Option<Box<dyn Solver>>, SolverError> {
    // The last of these that is set wins, as in the Fortran.
    let mode = ["coupled solver", "block solver", "assembly solver"]
        .into_iter()
        .rfind(|k| params.logical(k) == Some(true));
    match mode {
        None => Ok(None),
        Some("assembly solver") if params.get("procedure").is_none() => {
            Ok(Some(Box::new(AssemblySolver::from_section(params, registry)?)))
        }
        Some("assembly solver") => Err(SolverError::Unsupported(
            "Assembly Solver with a Procedure (CoupledSolver inheriting the matrix)".into(),
        )),
        Some("coupled solver") if params.get("procedure").is_none() => {
            Ok(Some(Box::new(CoupledSystem::from_section(params, registry)?)))
        }
        Some(m) => Err(SolverError::Unsupported(format!("{m} = True"))),
    }
}

/// A `Coupled Solver` section without a `Procedure`: its `Variable` is the
/// solver's own monolithic vector, created on the first visit, not at setup
/// (AddEquationBasics skips it, MainUtils.F90:1514).
pub fn is_coupled_system(params: &ParamMap) -> bool {
    let mode = ["coupled solver", "block solver", "assembly solver"]
        .into_iter()
        .rfind(|k| params.logical(k) == Some(true));
    mode == Some("coupled solver") && params.get("procedure").is_none()
}

/// `CoupledSolver` for `Assembly Solver = True`: the solver's own variable
/// and matrix, one assembly routine for the (1,1) block.
pub struct AssemblySolver {
    bulk: Box<dyn ElementAssembly>,
    /// `ProcName` as Elmer prints it: library and entry.
    proc_name: String,
}

impl AssemblySolver {
    /// Resolves `Bulk Assembly Procedure 11`, else `Bulk Assembly
    /// Procedure` (BlockSystemAssembly). `Boundary Assembly Procedure` is
    /// not read: `BlockSystemAssembly` switches to boundary mode after the
    /// bulk loop but its `GOTO 100` is commented out (MainUtils.F90:5156),
    /// so Elmer never calls the boundary routine.
    pub fn from_section(params: &ParamMap, registry: &SolverRegistry) -> Result<Self, SolverError> {
        let (lib, entry) = procedure_value(params, "bulk assembly procedure 11")
            .or_else(|| procedure_value(params, "bulk assembly procedure"))
            .ok_or_else(|| SolverError::fatal("BlockSystemAssembly", "Bulk Assembly Procedure not given!"))?;
        let proc_name = format!("{lib} {entry}");
        let f = registry
            .assembly(&lib, &entry)
            .ok_or_else(|| SolverError::fatal("BlockSystemAssembly", format!("Assembly routine not found: {proc_name}")))?;
        Ok(Self { bulk: f(params)?, proc_name })
    }

    /// `BlockSystemAssembly(Solver, dt, Transient, 1, 1)`: every bulk
    /// element whose nodes all have a dof, through the routine, the
    /// solver's time discretisation and the global glue.
    fn assemble(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info("BlockSystemAssembly", "Starting block system assembly", 8);
        ctx.info("BlockSystemAssembly", &format!("Using assembly routine: {}", self.proc_name), 8);
        // "These may be fetched within the assembly routine, if needed."
        let p = &mut ctx.sim.slots[ctx.slot].params;
        p.insert("Block Matrix Row", Value::Integer(vec![1]));
        p.insert("Block Matrix Column", Value::Integer(vec![1]));

        let dofs = ctx.primary().map_or(1, |v| v.dofs);
        let time_order = ctx.sim.slots[ctx.slot].time_order;
        let mut slots = Vec::new();
        for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
            ctx.element_slots(e, &mut slots);
            if slots.contains(&NO_DOF) {
                continue;
            }
            let n = slots.len() * dofs;
            let mut local = LocalSystem::zeros(n);
            local.mass = Some(vec![0.0; n * n]);
            local.damp = Some(vec![0.0; n * n]);
            self.bulk.assemble(ctx, e, dt, transient, &mut local)?;
            if transient {
                match time_order {
                    1 => ctx.default_first_order_time(&mut local, e, dt)?,
                    2 => ctx.default_second_order_time(&mut local, e, dt)?,
                    _ => {}
                }
            }
            ctx.default_update(&local, e)?;
        }
        Ok(())
    }
}

impl Solver for AssemblySolver {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Solving one equation", 5);
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Using existing variable and matrix", 8);
        if ctx.primary().is_none() {
            return Err(SolverError::fatal(CALLER, "In legacy solver mode the Variable should exist!"));
        }
        let p = ctx.params();
        if p.logical("eigen analysis") == Some(true) || p.logical("harmonic analysis") == Some(true) {
            return Err(SolverError::Unsupported("eigen/harmonic analysis in an assembly solver".into()));
        }
        if p.logical("nonlinear system linesearch") == Some(true) {
            return Err(SolverError::Unsupported("Nonlinear System Linesearch in an assembly solver".into()));
        }
        let iterations = p.integer("nonlinear system max iterations").unwrap_or(1);
        let tol = ctx.creal(ctx.params(), "nonlinear system convergence tolerance")?.unwrap_or(0.0);

        ctx.info(CALLER, LONG_RULE, 5);
        for iter in 1..=iterations {
            // WRITE(Message,'(A,T35,I5)') 'Coupled iteration:',iter
            ctx.info(CALLER, &format!("{:<34}{:>5}", "Coupled iteration:", iter), 5);
            ctx.default_initialize()?;
            self.assemble(ctx, dt, transient)?;
            ctx.default_dirichlet()?;
            ctx.default_solve()?;
            ctx.info(CALLER, LONG_RULE, 5);
            let change = ctx.convergence().map_or(0.0, |c| c.nonlin_change);
            if change < tol {
                break;
            }
        }
        ctx.info(CALLER, "All done", 8);
        ctx.info(CALLER, LONG_RULE, 5);
        Ok(())
    }
}

/// A constraint type (`Constraint n Type`).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum ConsType {
    /// One Lagrange multiplier: `integral of the variable = Constraint n Value`.
    Integral,
    /// One shared value for every coupled dof (a floating potential).
    Floating,
    /// Dof by dof equality of two variables.
    Equality,
}

/// The block layout of one visit: variables, constraints and their
/// offsets in the monolithic vector (`Offset`, `VarDofs`, `VarSizes`).
struct Layout {
    vars: Vec<VarId>,
    /// `Constraint n` names and the variables of those with dofs.
    cons: Vec<(String, Option<VarId>)>,
    offsets: Vec<usize>,
    dofs: Vec<usize>,
    total: usize,
    /// The type of the last constraint. `CoupledConstraintAssembly` tests
    /// this one variable for every constraint, as the Fortran does.
    cons_type: Option<ConsType>,
}

/// The CRS pattern `CoupledSolver` builds on its first visit from a list
/// matrix: one assembly pass that only records entries (`IsListMatrix`).
#[derive(Default)]
struct Pattern {
    rows: Vec<BTreeSet<u32>>,
}

impl Pattern {
    fn add(&mut self, row: usize, col: usize) {
        if self.rows.len() <= row {
            self.rows.resize_with(row + 1, BTreeSet::new);
        }
        self.rows[row].insert(col as u32);
    }

    /// `List_MoveRow`: row `from`'s entries also in row `to`.
    fn move_row(&mut self, from: usize, to: usize) {
        let cols: Vec<u32> = self.rows.get(from).map(|r| r.iter().copied().collect()).unwrap_or_default();
        for c in cols {
            self.add(to, c as usize);
        }
    }

    /// `List_ToCRSMatrix`: rows up to the last non-empty one, columns sorted.
    fn to_crs(&self) -> Result<CsrMatrix<f64>, SolverError> {
        let n = self.rows.iter().rposition(|r| !r.is_empty()).map_or(0, |i| i + 1);
        let mut row_ptr = vec![0];
        let mut col_idx = Vec::new();
        let mut diag = Vec::with_capacity(n);
        for (i, r) in self.rows[..n].iter().enumerate() {
            if let Some(&c) = r.iter().next_back() {
                if c as usize >= n {
                    return Err(SolverError::fatal(CALLER, format!("column {} outside the {n} rows", c + 1)));
                }
            }
            diag.push(r.contains(&(i as u32)).then(|| col_idx.len() + r.range(..i as u32).count()));
            col_idx.extend(r.iter().copied());
            row_ptr.push(col_idx.len());
        }
        let nnz = col_idx.len();
        Ok(CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr,
            col_idx,
            values: vec![0.0; nnz],
            diag: diag.into_iter().collect(),
        })
    }
}

/// Entry `(i, j)` of `a`, if stored.
fn entry(a: &CsrMatrix<f64>, i: usize, j: usize) -> Option<usize> {
    let r = a.row_ptr[i]..a.row_ptr[i + 1];
    a.col_idx[r.clone()].binary_search(&(j as u32)).ok().map(|k| r.start + k)
}

/// `CRS_AddToMatrixElement`; the pattern pass made every entry exist.
fn add_to(a: &mut CsrMatrix<f64>, i: usize, j: usize, v: f64) -> Result<(), SolverError> {
    let k = entry(a, i, j)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("no matrix entry ({}, {})", i + 1, j + 1)))?;
    a.values[k] += v;
    Ok(())
}

/// `CRS_SetMatrixElement`.
fn set(a: &mut CsrMatrix<f64>, i: usize, j: usize, v: f64) -> Result<(), SolverError> {
    let k = entry(a, i, j)
        .ok_or_else(|| SolverError::fatal(CALLER, format!("no matrix entry ({}, {})", i + 1, j + 1)))?;
    a.values[k] = v;
    Ok(())
}

/// `CRS_MoveRow(A, n1, n2, coeff)`: every entry of row `n1` above `TINY`
/// is zeroed and added, times `coeff`, to row `n2`.
fn move_row(a: &mut CsrMatrix<f64>, n1: usize, n2: usize, coeff: f64) -> Result<(), SolverError> {
    for k in a.row_ptr[n1]..a.row_ptr[n1 + 1] {
        let (j, val) = (a.col_idx[k] as usize, a.values[k]);
        if val.abs() > f64::MIN_POSITIVE {
            a.values[k] = 0.0 * val;
            add_to(a, n2, j, coeff * val)?;
        }
    }
    Ok(())
}

/// `CoupledSolver` for `Coupled Solver = True` without a `Procedure`
/// (MainUtils.F90:3565): the block variables `Variable 1..n` (created with
/// `CreateBlockVariable` when missing) and the constraints `Constraint
/// 1..m` share one matrix; block `(i, j)` is assembled by
/// `Bulk Assembly Procedure ij`, the constraints by
/// `CoupledConstraintAssembly`, Dirichlet conditions are set per variable
/// with its offset, and the solution is copied back into the variables.
///
/// Integral constraints need [`crate::Backends::fe`]. Transient runs, line
/// search, `Symmetric Assembly` and constraint types other than `integral`,
/// `floating` and `equality` are refused.
pub struct CoupledSystem {
    /// Routines by `ProcName`; blocks naming the same routine share one
    /// instance, as one loaded symbol shares its `SAVE`d state.
    routines: Vec<(String, Box<dyn ElementAssembly>)>,
    /// 1-based `(row, column)` block -> index into `routines`.
    blocks: BTreeMap<(usize, usize), usize>,
}

impl CoupledSystem {
    pub fn from_section(params: &ParamMap, registry: &SolverRegistry) -> Result<Self, SolverError> {
        for k in ["symmetric assembly", "antisymmetric assembly", "nonlinear system linesearch"] {
            if params.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("{k} in a Coupled Solver")));
            }
        }
        let nvar = (1..).take_while(|i| params.get(&format!("variable {i}")).is_some()).count();
        let mut routines: Vec<(String, Box<dyn ElementAssembly>)> = Vec::new();
        let mut factories: BTreeMap<String, usize> = BTreeMap::new();
        let mut blocks = BTreeMap::new();
        for r in 1..=nvar {
            for c in 1..=nvar {
                let mut found = procedure_value(params, &format!("bulk assembly procedure {r}{c}"));
                if found.is_none() && r == 1 && c == 1 {
                    found = procedure_value(params, "bulk assembly procedure");
                }
                let Some((lib, entry)) = found else { continue };
                let name = format!("{lib} {entry}");
                let ix = match factories.get(&name) {
                    Some(&ix) => ix,
                    None => {
                        let f: crate::solver::AssemblyFactory = registry.assembly(&lib, &entry).ok_or_else(|| {
                            SolverError::fatal("BlockSystemAssembly", format!("Assembly routine not found: {name}"))
                        })?;
                        routines.push((name.clone(), f(params)?));
                        factories.insert(name, routines.len() - 1);
                        routines.len() - 1
                    }
                };
                blocks.insert((r, c), ix);
            }
        }
        if nvar == 0 {
            return Err(SolverError::fatal(CALLER, "Coupled Solver without > Variable 1 <"));
        }
        Ok(Self { routines, blocks })
    }

    /// The variable part of `CoupledSolver`: block and constraint variables
    /// (created when missing), offsets and the monolithic variable.
    fn layout(&self, ctx: &mut SolverCtx<'_>) -> Result<Layout, SolverError> {
        let p = ctx.params().clone();
        let mut vars = Vec::new();
        let mut offsets = vec![0];
        let mut dofs = Vec::new();
        for i in 1.. {
            let Some(name) = p.string(&format!("variable {i}")) else { break };
            let name = normalize_key(name);
            let id = match ctx.sim.vars.id(&name) {
                Some(id) => id,
                None => {
                    ctx.info(CALLER, &format!("Variable {name} does not exist, creating"), 10);
                    create_block_variable(ctx, i, &name, &p)?
                }
            };
            let v = ctx.sim.vars.by_id(id);
            dofs.push(v.dofs);
            offsets.push(offsets[i - 1] + v.values.len());
            vars.push(id);
        }
        let mut cons = Vec::new();
        let mut cons_type = None;
        for i in 1.. {
            let Some(name) = p.string(&format!("constraint {i}")) else { break };
            let name = normalize_key(name);
            let ty = p.string(&format!("constraint {i} type")).map(normalize_key).ok_or_else(|| {
                SolverError::fatal("ListGetString", format!("Requested property: [constraint {i} type] not found"))
            })?;
            let (ty, con_dofs) = match ty.as_str() {
                "integral" => (ConsType::Integral, 1),
                "floating" => (ConsType::Floating, 1),
                "equality" => (ConsType::Equality, 0),
                other => return Err(SolverError::Unsupported(format!("Constraint {i} Type = {other}"))),
            };
            cons_type = Some(ty);
            let var = if con_dofs > 0 {
                Some(match ctx.sim.vars.id(&name) {
                    Some(id) => id,
                    None => {
                        ctx.info(CALLER, &format!("Constraint {name} does not exist, creating"), 10);
                        let mut v = Variable::new(&name, FieldKind::Global, con_dofs, None, vec![0.0; con_dofs]);
                        v.output = false;
                        ctx.sim.vars.add(v)?
                    }
                })
            } else {
                None
            };
            dofs.push(con_dofs);
            // VarSizes(j) = 1 whatever the dofs.
            let last = *offsets.last().expect("one entry");
            offsets.push(last + 1);
            cons.push((name, var));
        }
        for (j, off) in offsets[..offsets.len() - 1].iter().enumerate() {
            ctx.info(CALLER, &format!("{:<34}{off}", format!("Permutation offset {}: ", j + 1)), 8);
        }
        let total = *offsets.last().expect("one entry");
        offsets.pop();
        ctx.info(CALLER, &format!("{:<34}{}", "Number of coupled variables: ", vars.len()), 8);
        ctx.info(CALLER, &format!("{:<34}{}", "Number of constraints: ", cons.len()), 8);
        ctx.info(CALLER, &format!("{:<34}{total}", "Size of coupled system: "), 8);

        let layout = Layout { vars, cons, offsets, dofs, total, cons_type };
        if ctx.sim.slots[ctx.slot].variable.is_none() {
            let name = match p.string("variable") {
                Some(n) => normalize_key(n),
                None => {
                    ctx.info(CALLER, "New coupled variable added: Alldofs", 6);
                    "alldofs".to_string()
                }
            };
            let id = match ctx.sim.vars.id(&name) {
                Some(id) => id,
                None => {
                    let perm: Arc<[u32]> = (0..total as u32).collect();
                    let mut v = Variable::new(&name, FieldKind::Nodal, 1, Some(perm.clone()), vec![0.0; total]);
                    v.output = false;
                    let id = ctx.sim.vars.add(v)?;
                    let slot = &mut ctx.sim.slots[ctx.slot];
                    slot.perm = Some(perm);
                    slot.variable = Some(id);
                    if p.logical("coupled initial guess") == Some(true) {
                        ctx.info(CALLER, "Copying an initial guess", 8);
                        layout.copy(ctx, id, true);
                    }
                    id
                }
            };
            ctx.sim.slots[ctx.slot].variable = Some(id);
        }
        Ok(layout)
    }

    /// `BlockSystemAssembly(Solver, dt, Transient, r, c, Offset(r), Offset(c))`.
    fn assemble_block(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        lay: &Layout,
        (r, c): (usize, usize),
        pattern: Option<&mut Pattern>,
        dt: f64,
    ) -> Result<(), SolverError> {
        ctx.info("BlockSystemAssembly", "Starting block system assembly", 8);
        let Some(&ix) = self.blocks.get(&(r, c)) else {
            if r == c {
                ctx.warn("BlockSystemAssembly", "Diagonal bulk entries should be assembled!");
            }
            return Ok(());
        };
        ctx.info("BlockSystemAssembly", &format!("Using assembly routine: {}", self.routines[ix].0), 8);
        let p = &mut ctx.sim.slots[ctx.slot].params;
        p.insert("Block Matrix Row", Value::Integer(vec![r as i32]));
        p.insert("Block Matrix Column", Value::Integer(vec![c as i32]));

        let (rv, cv) = (ctx.sim.vars.by_id(lay.vars[r - 1]), ctx.sim.vars.by_id(lay.vars[c - 1]));
        let (rd, cd) = (rv.dofs, cv.dofs);
        if rd != cd {
            return Err(SolverError::Unsupported(format!(
                "Coupled Solver block ({r},{c}) of {rd}- and {cd}-dof variables"
            )));
        }
        let (rperm, cperm) = (rv.perm.clone(), cv.perm.clone());
        let slot_of = |perm: &Option<Arc<[u32]>>, n: u32| match perm {
            Some(p) => p.get(n as usize).copied().filter(|&s| s != NO_DOF),
            None => Some(n),
        };
        let (r0, c0) = (lay.offsets[r - 1], lay.offsets[c - 1]);
        let mut pattern = pattern;
        for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
            let nodes = ctx.mesh.element(e).nodes;
            let Some(ri) = nodes.iter().map(|&n| slot_of(&rperm, n)).collect::<Option<Vec<u32>>>() else { continue };
            let Some(ci) = nodes.iter().map(|&n| slot_of(&cperm, n)).collect::<Option<Vec<u32>>>() else { continue };
            let row = |i: usize, a: usize| r0 + rd * ri[i] as usize + a;
            let col = |j: usize, b: usize| c0 + cd * ci[j] as usize + b;
            if let Some(pat) = pattern.as_deref_mut() {
                for i in 0..ri.len() {
                    for a in 0..rd {
                        for j in 0..ci.len() {
                            for b in 0..cd {
                                pat.add(row(i, a), col(j, b));
                            }
                        }
                    }
                }
                continue;
            }
            let n = nodes.len() * rd;
            let mut local = LocalSystem::zeros(n);
            local.mass = Some(vec![0.0; n * n]);
            local.damp = Some(vec![0.0; n * n]);
            self.routines[ix].1.assemble(ctx, e, dt, false, &mut local)?;
            let sys = ctx.sim.slots[ctx.slot].system.as_mut().expect("pattern built");
            for i in 0..ri.len() {
                for a in 0..rd {
                    for j in 0..ci.len() {
                        for b in 0..cd {
                            add_to(&mut sys.matrix, row(i, a), col(j, b), local.stiff[(rd * i + a) * n + cd * j + b])?;
                        }
                    }
                    sys.rhs[row(i, a)] += local.force[rd * i + a];
                }
            }
        }
        Ok(())
    }

    /// `CoupledConstraintAssembly` (MainUtils.F90:4052): bulk elements with
    /// the constraint's name set in their body force, then boundary elements
    /// with it set in their boundary condition.
    fn assemble_constraints(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        lay: &Layout,
        mut pattern: Option<&mut Pattern>,
        all_dir: &mut [bool],
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "Starting constraint assembly", 8);
        let nvar = lay.vars.len();
        let ty = lay.cons_type.expect("constraints present");
        let p = ctx.params().clone();
        let nb = ctx.mesh.n_bulk() as u32;
        for bulk in [true, false] {
            let elements = if bulk { 0..nb } else { nb..nb + ctx.mesh.n_boundary() as u32 };
            for (ci, (row_name, _)) in lay.cons.iter().enumerate() {
                let k1 = ci + 1;
                let mut row_dofs = lay.dofs[nvar + ci];
                let mut row0 = lay.offsets[nvar + ci];
                let var_inds = match p.get(&format!("constraint {k1} variables")) {
                    Some(Value::Integer(v)) if !v.is_empty() => v.clone(),
                    _ => return Err(SolverError::fatal(CALLER, "Cannot continue without pointer to variables")),
                };
                let target = p.integer(&format!("constraint {k1} components")).unwrap_or(0);
                let value = ctx.creal(&p, &format!("constraint {k1} value"))?.unwrap_or(0.0);
                let coeff = ctx.creal(&p, &format!("constraint {k1} coeff"))?.unwrap_or(1.0);
                let var_at = |k: i32| -> Result<usize, SolverError> {
                    usize::try_from(k - 1)
                        .ok()
                        .filter(|&k| k < nvar)
                        .ok_or_else(|| SolverError::fatal(CALLER, format!("Constraint {k1} Variables: no variable {k}")))
                };
                let mut row_perm = None;
                if ty == ConsType::Equality {
                    let v2 = var_at(*var_inds.get(1).ok_or_else(|| {
                        SolverError::fatal(CALLER, format!("Constraint {k1} Variables needs two variables"))
                    })?)?;
                    row_perm = Some(ctx.sim.vars.by_id(lay.vars[v2]).perm.clone());
                    row_dofs = lay.dofs[v2];
                    row0 = lay.offsets[v2];
                } else {
                    for j in 0..row_dofs {
                        match pattern.as_deref_mut() {
                            Some(pat) => pat.add(row0 + j, row0 + j),
                            None => {
                                let sys = ctx.sim.slots[ctx.slot].system.as_mut().expect("pattern built");
                                add_to(&mut sys.matrix, row0 + j, row0 + j, 0.0)?;
                            }
                        }
                    }
                }
                let cvix = var_at(var_inds[0])?;
                let (col_perm, col_dofs) = {
                    let v = ctx.sim.vars.by_id(lay.vars[cvix]);
                    (v.perm.clone(), v.dofs)
                };
                let col0 = lay.offsets[cvix];
                let slot_of = |perm: &Option<Arc<[u32]>>, n: u32| match perm {
                    Some(p) => p.get(n as usize).copied().filter(|&s| s != NO_DOF),
                    None => Some(n),
                };
                for e in elements.clone().map(ElemId) {
                    let el = ctx.mesh.element(e);
                    let nodes = el.nodes;
                    let row_inds: Vec<u32> = match &row_perm {
                        Some(rp) => match nodes.iter().map(|&n| slot_of(rp, n)).collect::<Option<Vec<u32>>>() {
                            Some(v) => v,
                            None => continue,
                        },
                        None => vec![0],
                    };
                    let Some(col_inds) = nodes.iter().map(|&n| slot_of(&col_perm, n)).collect::<Option<Vec<u32>>>()
                    else {
                        continue;
                    };
                    let section = if bulk {
                        ctx.body_params(e).and_then(|b| b.body_force)
                    } else {
                        ctx.bc_params(e)
                    };
                    if section.and_then(|s| s.logical(row_name)) != Some(true) {
                        continue;
                    }
                    if matches!(ty, ConsType::Floating | ConsType::Equality) {
                        for (i, &cslot) in col_inds.iter().enumerate() {
                            for k in 0..col_dofs {
                                // Col = ColInd0 + ColDofs * ColInds(i) - k, 1-based.
                                let col = col0 + col_dofs * cslot as usize + (col_dofs - 1 - k);
                                if all_dir[col] {
                                    continue;
                                }
                                all_dir[col] = true;
                                let row = if ty == ConsType::Floating {
                                    row0
                                } else {
                                    row0 + row_dofs * row_inds[i] as usize + row_dofs - 1 - k
                                };
                                if let Some(pat) = pattern.as_deref_mut() {
                                    pat.move_row(col, row);
                                    pat.add(col, col);
                                    pat.add(col, row);
                                    continue;
                                }
                                let sys = ctx.sim.slots[ctx.slot].system.as_mut().expect("pattern built");
                                move_row(&mut sys.matrix, col, row, coeff)?;
                                sys.rhs[row] += sys.rhs[col];
                                sys.rhs[col] = 0.0;
                                sys.rhs[row] += value;
                                set(&mut sys.matrix, col, col, 1.0)?;
                                set(&mut sys.matrix, col, row, -coeff)?;
                            }
                        }
                        continue;
                    }
                    // IntegralConstraint: STIFF(1,p) and FORCE(p) over the element.
                    let ncol = col_inds.len();
                    let (mut stiff, mut force) = (vec![0.0; ncol], vec![0.0; ncol]);
                    if pattern.is_none() {
                        let fe = ctx
                            .sim
                            .backends
                            .fe
                            .as_ref()
                            .ok_or_else(|| SolverError::NotWired("basis and quadrature for integral constraints (elmer-fe)".into()))?;
                        let xyz: Vec<[f64; 3]> = nodes.iter().map(|&n| ctx.mesh.coord(n)).collect();
                        let mut ev = BasisEval::default();
                        for qp in fe.quadrature.gauss_points(el.etype, None) {
                            fe.basis
                                .eval(el.etype, &xyz, qp.uvw, &mut ev)
                                .map_err(|m| SolverError::fatal("IntegralConstraint", m))?;
                            let w = qp.weight * ev.det_j;
                            for pnode in 0..ncol {
                                stiff[pnode] += w * ev.basis[pnode];
                                force[pnode] += value * w * ev.basis[pnode];
                            }
                        }
                    }
                    for k in 1..=col_dofs {
                        if target != 0 && target as usize != k {
                            continue;
                        }
                        let cols: Vec<usize> = col_inds.iter().map(|&s| col0 + col_dofs * s as usize + k - 1).collect();
                        if let Some(pat) = pattern.as_deref_mut() {
                            for r in 0..row_dofs {
                                for &c in &cols {
                                    pat.add(row0 + r, c);
                                    pat.add(c, row0 + r);
                                }
                            }
                            continue;
                        }
                        let sys = ctx.sim.slots[ctx.slot].system.as_mut().expect("pattern built");
                        for r in 0..row_dofs {
                            for (j, &c) in cols.iter().enumerate() {
                                // Only row 1 of STIFF is filled.
                                let s = if r == 0 { stiff[j] } else { 0.0 };
                                add_to(&mut sys.matrix, row0 + r, c, s)?;
                                add_to(&mut sys.matrix, c, row0 + r, s)?;
                            }
                        }
                        // FORCE(RowDofs*(i-1)+j) with Nrow = 1: FORCE(1..RowDofs).
                        for j in 0..row_dofs {
                            sys.rhs[row0 + j] += force.get(j).copied().unwrap_or(0.0);
                        }
                    }
                }
                if ty == ConsType::Integral && pattern.is_none() {
                    let sys = ctx.sim.slots[ctx.slot].system.as_mut().expect("pattern built");
                    for j in 0..row_dofs {
                        sys.rhs[row0 + j] = value;
                    }
                }
            }
        }
        Ok(())
    }

    /// `CoupledSystemDirichlet`: `DefaultDirichletBCs` per variable with its
    /// offset. Elmer keeps every constrained dof of the earlier calls and
    /// enforces them again, so each call gets the dofs so far.
    fn dirichlet(&self, ctx: &mut SolverCtx<'_>, lay: &Layout) -> Result<(), SolverError> {
        ctx.info(CALLER, "Setting coupled system Dirichlet conditions", 6);
        let mut fixed: BTreeMap<u32, f64> = BTreeMap::new();
        for (i, &id) in lay.vars.iter().enumerate() {
            fixed.extend(ctx.dirichlet_values_of(id, lay.offsets[i])?);
            let all: Vec<(u32, f64)> = fixed.iter().map(|(&r, &v)| (r, v)).collect();
            let sim = &mut *ctx.sim;
            let backend = sim
                .backends
                .assembly
                .as_ref()
                .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
            let slot = &mut sim.slots[ctx.slot];
            let sys = slot.system.as_mut().expect("pattern built");
            backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, &all, &slot.params)?;
        }
        Ok(())
    }
}

impl Layout {
    /// `CoupledToSingleVector` (`to_coupled = false`) or
    /// `SingleToCoupledVector`.
    fn copy(&self, ctx: &mut SolverCtx<'_>, all: VarId, to_coupled: bool) {
        let nvar = self.vars.len();
        for i in 0..nvar + self.cons.len() {
            if self.dofs[i] == 0 {
                continue;
            }
            let id = if i < nvar { self.vars[i] } else { self.cons[i - nvar].1.expect("constraint with dofs") };
            let off = self.offsets[i];
            let len = ctx.sim.vars.by_id(id).values.len();
            if to_coupled {
                let src = ctx.sim.vars.by_id(id).values.clone();
                ctx.sim.vars.by_id_mut(all).values[off..off + len].copy_from_slice(&src);
            } else {
                let src = ctx.sim.vars.by_id(all).values[off..off + len].to_vec();
                let v = ctx.sim.vars.by_id_mut(id);
                v.values.copy_from_slice(&src);
                v.values_changed = true;
            }
        }
    }
}

/// `CreateBlockVariable` (BlockSolve.F90:138): the nodes of bulk elements
/// whose equation lists the variable in `Active Variables[solver]` (else
/// `Active Variables`), numbered in node order.
fn create_block_variable(ctx: &mut SolverCtx<'_>, no: usize, name: &str, p: &ParamMap) -> Result<VarId, SolverError> {
    let dofs = p.integer(&format!("variable {no} dofs")).unwrap_or(1).max(1) as usize;
    let solver_id = ctx.solver_id();
    let mut mark = vec![false; ctx.mesh.n_nodes()];
    for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
        let Some(eq) = ctx.body_params(e).and_then(|b| b.equation) else { continue };
        let active = match eq.get(&format!("active variables[{solver_id}]")).or_else(|| eq.get("active variables")) {
            Some(Value::Integer(v)) => v,
            _ => continue,
        };
        if active.contains(&(no as i32)) {
            for &n in ctx.mesh.element(e).nodes {
                mark[n as usize] = true;
            }
        }
    }
    let mut j = 0u32;
    let perm: Arc<[u32]> = mark
        .iter()
        .map(|&m| {
            if m {
                j += 1;
                j - 1
            } else {
                NO_DOF
            }
        })
        .collect();
    let v = Variable::new(name, FieldKind::Nodal, dofs, Some(perm), vec![0.0; j as usize * dofs]);
    ctx.sim.vars.add(v)
}

impl Solver for CoupledSystem {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("transient Coupled Solver".into()));
        }
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Solving a system of equations", 5);
        ctx.info(CALLER, RULE, 5);
        let lay = self.layout(ctx)?;
        let all = ctx.sim.slots[ctx.slot].variable.expect("set by layout");
        let iterations = ctx.params().integer("nonlinear system max iterations").unwrap_or(1);
        let tol = ctx.creal(ctx.params(), "nonlinear system convergence tolerance")?.unwrap_or(0.0);
        let mut all_dir = vec![false; lay.total];

        ctx.info(CALLER, LONG_RULE, 5);
        for iter in 1..=iterations {
            ctx.info(CALLER, &format!("{:<34}{:>5}", "Coupled iteration:", iter), 5);
            let mut list = ctx.sim.slots[ctx.slot].system.is_none();
            loop {
                let mut pattern = list.then(Pattern::default);
                if !list {
                    ctx.default_initialize()?;
                }
                all_dir.iter_mut().for_each(|d| *d = false);
                let n = lay.vars.len();
                for r in 1..=n {
                    for c in 1..=n {
                        self.assemble_block(ctx, &lay, (r, c), pattern.as_mut(), dt)?;
                    }
                }
                if !lay.cons.is_empty() {
                    self.assemble_constraints(ctx, &lay, pattern.as_mut(), &mut all_dir)?;
                }
                let Some(pat) = pattern else { break };
                let matrix = pat.to_crs()?;
                let rhs = vec![0.0; matrix.nrows];
                ctx.sim.slots[ctx.slot].system = Some(LinearSystem { matrix, rhs });
                list = false;
            }
            self.dirichlet(ctx, &lay)?;
            ctx.default_solve()?;
            ctx.info(CALLER, LONG_RULE, 5);
            ctx.info(CALLER, "Copying results into original variables", 8);
            lay.copy(ctx, all, false);
            let change = ctx.convergence().map_or(0.0, |c| c.nonlin_change);
            if change < tol {
                break;
            }
        }
        ctx.info(CALLER, "All done", 8);
        ctx.info(CALLER, LONG_RULE, 5);
        Ok(())
    }
}
