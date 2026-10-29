//! `Block Solver = True` without a `Procedure` (`BlockSolver`,
//! MainUtils.F90:4323): the `Variable n` fields are the blocks, each block
//! `(r, c)` with a `Bulk Assembly Procedure rc` is assembled into its own
//! matrix (`BlockSystemAssembly`) and gets the Dirichlet conditions of
//! variable `r` (`BlockSystemDirichlet`). One block is solved with
//! `DefaultSolve` ("in standard manner"), more through `elmer-block`
//! (`BlockKrylovIter`, or `BlockStandardIter` with `Block Preconditioner =
//! False`), until the change is below `Nonlinear System Convergence
//! Tolerance`. The solver's variable is then `Primary Variable` (default 1).

use std::collections::{BTreeMap, BTreeSet};

use elmer_core::{CsrMatrix, ElemId, LocalSystem, ParamMap, Params, NO_DOF};

use crate::coupled::{block_routines, create_block_variable, BlockRoutines};
use crate::error::SolverError;
use crate::sim::LinearSystem;
use crate::solver::{ElementAssembly, Solver, SolverCtx, SolverRegistry};
use crate::variable::VarId;

const CALLER: &str = "BlockSolver";
const RULE: &str = "---------------------------------------";
const LONG_RULE: &str = "-------------------------------------------------";

pub struct BlockAssembled {
    routines: Vec<(String, Box<dyn ElementAssembly>)>,
    blocks: BTreeMap<(usize, usize), usize>,
    nvar: usize,
    /// The block variables, once created.
    vars: Vec<VarId>,
    /// The assembled blocks, 1-based `(row, column)`.
    systems: BTreeMap<(usize, usize), LinearSystem>,
}

impl BlockAssembled {
    pub fn from_section(params: &ParamMap, registry: &SolverRegistry) -> Result<Self, SolverError> {
        for k in ["symmetric assembly", "antisymmetric assembly", "nonlinear system linesearch"] {
            if params.logical(k) == Some(true) {
                return Err(SolverError::Unsupported(format!("{k} in a Block Solver")));
            }
        }
        if params.get("block solvers").is_some() {
            return Err(SolverError::Unsupported("Block Solvers (slave solver blocks)".into()));
        }
        let (nvar, (routines, blocks)): (usize, BlockRoutines) = block_routines(params, registry)?;
        if nvar == 0 {
            return Err(SolverError::fatal(CALLER, "No variables, nothing to do!"));
        }
        if let Some(&(r, c)) = blocks.keys().find(|(r, c)| r != c) {
            return Err(SolverError::Unsupported(format!(
                "Block Solver off-diagonal block ({r},{c}) (off-diagonal Dirichlet conditions)"
            )));
        }
        Ok(Self {
            routines,
            blocks,
            nvar,
            vars: Vec::new(),
            systems: BTreeMap::new(),
        })
    }

    /// `BlockInitMatrix`: the block variables (created when missing) and
    /// the pattern of every block with a routine (the list-matrix pass).
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        ctx.info(CALLER, &format!("Number of coupled variables: {}", self.nvar), 4);
        for i in 1..=self.nvar {
            let name = elmer_core::normalize_key(p.string(&format!("variable {i}")).expect("counted"));
            let id = match ctx.sim.vars.id(&name) {
                Some(id) => id,
                None => create_block_variable(ctx, i, &name, &p)?,
            };
            if ctx.sim.vars.by_id(id).dofs != 1 {
                return Err(SolverError::Unsupported("Block Solver variables of several components".into()));
            }
            self.vars.push(id);
        }
        for &(r, c) in self.blocks.keys() {
            let (rp, cp) = (self.perm(ctx, r), self.perm(ctx, c));
            let mut rows: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); ctx.sim.vars.by_id(self.vars[r - 1]).values.len()];
            for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
                let Some((ri, ci)) = slots(ctx, e, &rp, &cp) else { continue };
                for &i in &ri {
                    rows[i as usize].extend(ci.iter().copied());
                }
            }
            let ncols = ctx.sim.vars.by_id(self.vars[c - 1]).values.len();
            let matrix = to_csr(&rows, ncols);
            let rhs = vec![0.0; matrix.nrows];
            self.systems.insert((r, c), LinearSystem { matrix, rhs });
        }
        Ok(())
    }

    fn perm(&self, ctx: &SolverCtx<'_>, i: usize) -> Vec<u32> {
        ctx.sim.vars.by_id(self.vars[i - 1]).perm.as_deref().map(<[u32]>::to_vec).unwrap_or_default()
    }

    /// `BlockSystemAssembly(Solver, dt, Transient, r, c)` into block
    /// `(r, c)`, then `BlockSystemDirichlet`.
    fn assemble(&mut self, ctx: &mut SolverCtx<'_>, (r, c): (usize, usize), dt: f64) -> Result<(), SolverError> {
        let ix = self.blocks[&(r, c)];
        ctx.info("BlockSystemAssembly", "Starting block system assembly", 8);
        ctx.info("BlockSystemAssembly", &format!("Using assembly routine: {}", self.routines[ix].0), 8);
        let p = &mut ctx.sim.slots[ctx.slot].params;
        p.insert("Block Matrix Row", elmer_core::Value::Integer(vec![r as i32]));
        p.insert("Block Matrix Column", elmer_core::Value::Integer(vec![c as i32]));
        let (rp, cp) = (self.perm(ctx, r), self.perm(ctx, c));
        let mut sys = self.systems.remove(&(r, c)).expect("pattern built");
        sys.matrix.values.iter_mut().for_each(|v| *v = 0.0);
        sys.rhs.iter_mut().for_each(|v| *v = 0.0);
        let mut result = Ok(());
        for e in (0..ctx.mesh.n_bulk() as u32).map(ElemId) {
            let Some((ri, ci)) = slots(ctx, e, &rp, &cp) else { continue };
            let n = ri.len();
            let mut local = LocalSystem::zeros(n);
            local.mass = Some(vec![0.0; n * n]);
            local.damp = Some(vec![0.0; n * n]);
            if let Err(err) = self.routines[ix].1.assemble(ctx, e, dt, false, &mut local) {
                result = Err(err);
                break;
            }
            for (a, &i) in ri.iter().enumerate() {
                for (b, &j) in ci.iter().enumerate() {
                    let i = i as usize;
                    let span = sys.matrix.row_ptr[i]..sys.matrix.row_ptr[i + 1];
                    let at = span.start + sys.matrix.col_idx[span].binary_search(&j).expect("in the pattern");
                    sys.matrix.values[at] += local.stiff[a * n + b];
                }
                sys.rhs[i as usize] += local.force[a];
            }
        }
        if result.is_ok() && sys.matrix.nrows > 0 {
            ctx.info(CALLER, "Setting block system Dirichlet conditions", 8);
            result = (|| {
                let fixed = ctx.dirichlet_values_of(self.vars[r - 1], 0)?;
                let sim = &mut *ctx.sim;
                let backend = sim
                    .backends
                    .assembly
                    .as_ref()
                    .ok_or_else(|| SolverError::NotWired("assembly backend (elmer-assembly)".into()))?;
                backend.try_apply_dirichlet(&mut sys.matrix, &mut sys.rhs, &fixed, &sim.slots[ctx.slot].params)
            })();
        }
        self.systems.insert((r, c), sys);
        result
    }
}

/// The row and column slots of bulk element `e`, when all its nodes have
/// dofs in both permutations.
fn slots(ctx: &SolverCtx<'_>, e: ElemId, rp: &[u32], cp: &[u32]) -> Option<(Vec<u32>, Vec<u32>)> {
    let nodes = ctx.mesh.element(e).nodes;
    let pick = |perm: &[u32]| -> Option<Vec<u32>> {
        nodes
            .iter()
            .map(|&k| perm.get(k as usize).copied().filter(|&s| s != NO_DOF))
            .collect()
    };
    Some((pick(rp)?, pick(cp)?))
}

/// Sorted rows as a CSR matrix with its diagonal index.
fn to_csr(rows: &[BTreeSet<u32>], ncols: usize) -> CsrMatrix<f64> {
    let mut a = CsrMatrix {
        nrows: rows.len(),
        ncols,
        row_ptr: vec![0],
        col_idx: Vec::new(),
        values: Vec::new(),
        diag: None,
    };
    let mut diag = Vec::with_capacity(rows.len());
    for (i, r) in rows.iter().enumerate() {
        diag.push(r.contains(&(i as u32)).then(|| a.col_idx.len() + r.range(..i as u32).count()));
        a.col_idx.extend(r.iter().copied());
        a.row_ptr.push(a.col_idx.len());
    }
    a.values = vec![0.0; a.col_idx.len()];
    a.diag = diag.into_iter().collect();
    a
}

impl Solver for BlockAssembled {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, dt: f64, transient: bool) -> Result<(), SolverError> {
        if transient {
            return Err(SolverError::Unsupported("transient Block Solver without a Procedure".into()));
        }
        ctx.info(CALLER, RULE, 5);
        ctx.info(CALLER, "Solving system of equations utilizing block strategies", 5);
        ctx.info(CALLER, RULE, 5);
        if self.vars.is_empty() {
            self.init(ctx)?;
        }
        let p = ctx.params().clone();
        let primary = p.integer("primary variable").unwrap_or(1).max(1) as usize;
        let primary = *self
            .vars
            .get(primary - 1)
            .ok_or_else(|| SolverError::fatal(CALLER, "Primary Variable out of range"))?;
        let block_prec = p.logical("block preconditioner").unwrap_or(true);
        let iterations = p.integer("nonlinear system max iterations").unwrap_or(1);
        let tol = ctx.creal(&p, "nonlinear system convergence tolerance")?.unwrap_or(0.0);
        let nvar = self.nvar;

        ctx.info(CALLER, LONG_RULE, 6);
        for iter in 1..=iterations.max(1) {
            ctx.info(CALLER, &format!("{:<34}{iter}", "Coupled iteration: "), 6);
            let keys: Vec<(usize, usize)> = self.blocks.keys().copied().collect();
            for rc in keys {
                self.assemble(ctx, rc, dt)?;
            }
            for i in 1..=nvar {
                if !self.systems.contains_key(&(i, i)) {
                    return Err(SolverError::fatal(
                        "BlockKrylovIter",
                        format!("Submatrix {} not associated!", 11 * i),
                    ));
                }
            }
            let max_change = if nvar == 1 {
                ctx.info(CALLER, "Solving in standard manner", 8);
                let slot = ctx.slot;
                let sys = self.systems.remove(&(1, 1)).expect("checked");
                ctx.sim.slots[slot].variable = Some(self.vars[0]);
                ctx.sim.slots[slot].system = Some(sys);
                let r = with_namespace(ctx, "outer:", |ctx| ctx.default_solve());
                let sys = ctx.sim.slots[slot].system.take().expect("put above");
                self.systems.insert((1, 1), sys);
                r?;
                ctx.convergence().map_or(0.0, |c| c.nonlin_change)
            } else {
                ctx.info(
                    CALLER,
                    if block_prec {
                        "Using block preconditioning strategy"
                    } else {
                        "Using block solution strategy"
                    },
                    8,
                );
                let mut mats = Vec::with_capacity(nvar * nvar);
                for r in 1..=nvar {
                    for c in 1..=nvar {
                        mats.push(self.systems.get(&(r, c)).map(|s| s.matrix.clone()));
                    }
                }
                let sys = elmer_block::BlockSystem::from_blocks(nvar, mats).map_err(crate::solver::block_error)?;
                let mut b = Vec::new();
                let mut x = Vec::new();
                for (i, &id) in self.vars.iter().enumerate() {
                    b.extend_from_slice(&self.systems[&(i + 1, i + 1)].rhs);
                    x.extend_from_slice(&ctx.sim.vars.by_id(id).values);
                }
                let sim_p = &ctx.sim.model.simulation;
                let global_abort = sim_p.logical("global abort not converged").unwrap_or(true);
                let opts = elmer_block::BlockOptions {
                    additive_namespaces: sim_p.logical("additive namespaces") == Some(true),
                    outer_scaling: None,
                    global_abort,
                    schur: None,
                };
                let inner = elmer_block::ElmerInner::new(ctx.log, ctx.comm, global_abort);
                let report =
                    elmer_block::block_solve_system(&p, &sys, &b, &mut x, &opts, Box::new(inner), ctx.log, ctx.comm)
                        .map_err(crate::solver::block_error)?;
                let mut off = 0;
                for &id in &self.vars {
                    let v = ctx.sim.vars.by_id_mut(id);
                    let n = v.values.len();
                    v.values.copy_from_slice(&x[off..off + n]);
                    v.values_changed = true;
                    off += n;
                }
                ctx.sim.publish_live();
                report.max_change
            };
            if max_change < tol {
                break;
            }
        }
        // The solver's own variable and matrix are the primary block's.
        let slot = &mut ctx.sim.slots[ctx.slot];
        slot.variable = Some(primary);
        ctx.info(CALLER, "All done", 8);
        ctx.info(CALLER, LONG_RULE, 6);
        Ok(())
    }
}

/// Runs `f` with the solver section's `<prefix> <key>` entries standing for
/// `<key>` (`ListPushNamespace`), putting the previous values back after.
pub(crate) fn with_namespace<T>(
    ctx: &mut SolverCtx<'_>,
    prefix: &str,
    f: impl FnOnce(&mut SolverCtx<'_>) -> T,
) -> T {
    let params = &mut ctx.sim.slots[ctx.slot].params;
    let scoped: Vec<(String, elmer_core::Value)> = params
        .0
        .iter()
        .filter_map(|(k, v)| Some((k.strip_prefix(prefix)?.trim_start().to_string(), v.clone())))
        .collect();
    let saved: Vec<(String, Option<elmer_core::Value>)> = scoped
        .into_iter()
        .map(|(k, v)| {
            let old = params.0.insert(k.clone(), v);
            (k, old)
        })
        .collect();
    let r = f(ctx);
    let params = &mut ctx.sim.slots[ctx.slot].params;
    for (k, old) in saved {
        match old {
            Some(v) => params.0.insert(k, v),
            None => params.0.remove(&k),
        };
    }
    r
}
