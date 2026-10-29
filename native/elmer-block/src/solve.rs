//! `BlockSolveInt` (BlockSolve.F90:5007) and `BlockKrylovIter` (:4417).

use elmer_core::{normalize_key, Communicator, CsrMatrix, Logger, ParamMap, Params, SolveStats, Value};

use crate::inner::{check_converged, iterative};
use crate::namespace::{block_namespace, resolve};
use crate::precond::PrecBlock;
use crate::standard::standard_iter;
use crate::{
    BlockError, BlockLayout, BlockOperator, BlockPreconditioner, BlockSystem, ElmerInner, InnerRecord, InnerSolver,
    PrecSettings,
};

/// The scaling `SolveLinearSystem` applied to the monolithic system before
/// the block solve (`Linear System Scaling`, default on): `A % DiagScaling`
/// and `A % RhsScaling`.
#[derive(Clone, Copy, Debug)]
pub struct OuterScaling<'a> {
    pub diag: &'a [f64],
    pub rhs_scaling: f64,
}

/// What the caller knows that the solver section does not.
#[derive(Clone, Copy, Debug)]
pub struct BlockOptions<'a> {
    /// `Additive Namespaces` (Simulation): keys are also looked up under the
    /// enclosing namespaces (`block:` and `outer:` for an inner solve).
    pub additive_namespaces: bool,
    /// The monolithic scaling, when active. The inner solves then skip their
    /// own scaling unless `Eliminate Nested Scaling = False`, and a Schur
    /// block is solved in the unscaled system.
    pub outer_scaling: Option<OuterScaling<'a>>,
    /// `Global Abort Not Converged` (Simulation, default true).
    pub global_abort: bool,
    /// The matrix of the solver named by `Block Matrix Schur Variable`
    /// (IncompressibleNS sets it with `Block Preconditioner`): the last
    /// block's preconditioning matrix (BlockPrecMatrix, BlockSolve.F90:2523).
    pub schur: Option<&'a CsrMatrix<f64>>,
}

impl Default for BlockOptions<'_> {
    fn default() -> Self {
        BlockOptions {
            additive_namespaces: false,
            outer_scaling: None,
            global_abort: true,
            schur: None,
        }
    }
}

/// What a block solve did.
#[derive(Clone, Debug)]
pub struct BlockReport {
    /// Number of blocks (`NoVar`).
    pub nblocks: usize,
    /// The outer iteration; `None` when the right-hand side was zero and the
    /// solver was not called.
    pub outer: Option<SolveStats>,
    /// Preconditioner sweeps (`TotMatrix % NoIters`).
    pub sweeps: usize,
    /// Every inner solve, in order.
    pub inner: Vec<InnerRecord>,
    /// `CompNorm` of the right-hand side and of the solution before and after.
    pub rhs_norm: f64,
    pub start_norm: f64,
    pub solution_norm: f64,
    /// `MaxChange`: `2 |x - x0| / (x + x0)` of the two norms.
    pub max_change: f64,
}

/// Solves `A x = b` with the block strategy of the solver section `params`
/// for a variable of `dofs` components; `x` holds the initial guess. Inner
/// solves use [`ElmerInner`].
#[allow(clippy::too_many_arguments)]
pub fn block_solve(
    params: &ParamMap,
    dofs: usize,
    a: &CsrMatrix<f64>,
    b: &[f64],
    x: &mut [f64],
    opts: &BlockOptions,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<BlockReport, BlockError> {
    let inner = ElmerInner::new(log, comm, opts.global_abort);
    block_solve_with(params, dofs, a, b, x, opts, Box::new(inner), log, comm)
}

/// Keywords selecting block strategies or features that are not ported.
const UNPORTED_LOGICALS: &[(&str, &str)] = &[
    ("block a-v system", "Block A-V System"),
    ("block hcurl system", "Block Hcurl System"),
    ("block quadratic hcurl system", "Block Quadratic Hcurl System"),
    ("block hdiv system", "Block Hdiv System"),
    ("block re-im system", "Block Re-Im System"),
    ("block nodal system", "Block Nodal System"),
    ("block hor-ver system", "Block Hor-Ver System"),
    ("block cartesian system", "Block Cartesian System"),
    ("block domain system", "Block Domain System"),
    ("block nested system", "Block Nested System"),
    ("block complex system", "Block Complex System"),
    ("block monolithic", "Block Monolithic"),
    ("block scaling", "Block Scaling"),
    ("block mv monolithic", "Block MV Monolithic"),
    ("block amgx m-v", "Block AMGX M-V"),
    ("block prec reuse", "Block Prec Reuse"),
    ("block split complex", "Block Split Complex"),
    ("block split complex hcurl", "Block Split Complex Hcurl"),
    ("create schur matrix approximation", "Create Schur Matrix Approximation"),
    ("structure-structure coupling", "Structure-Structure Coupling"),
    ("linear system complex", "complex block systems (Linear System Complex)"),
    (
        "linear system residual mode",
        "Linear System Residual Mode with a block solver",
    ),
];

/// Keys (or key prefixes) that are refused when present at all.
const UNPORTED_KEYS: &[(&str, &str)] = &[
    ("block solvers", "Block Solvers (blocks from other solvers)"),
    ("structure solver index", "fluid-structure coupling blocks"),
    ("plate solver index", "fluid-structure coupling blocks"),
    ("shell solver index", "fluid-structure coupling blocks"),
    ("schur operator", "Schur Operator"),
    ("schur multiplier", "Schur Multiplier"),
    ("prec matrix", "tailored preconditioning matrices (Prec Matrix ...)"),
];

fn refuse_unported(p: &ParamMap) -> Result<(), BlockError> {
    for (k, what) in UNPORTED_LOGICALS {
        if p.logical(k) == Some(true) {
            return Err(BlockError::Unsupported(format!("{what} in a block linear system")));
        }
    }
    for (k, what) in UNPORTED_KEYS {
        if p.0.keys().any(|key| key == k || key.starts_with(&format!("{k} "))) {
            return Err(BlockError::Unsupported(format!("{what} in a block linear system")));
        }
    }
    Ok(())
}

/// [`block_solve`] with a caller-supplied inner solver.
#[allow(clippy::too_many_arguments)]
pub fn block_solve_with<'a>(
    params: &ParamMap,
    dofs: usize,
    a: &CsrMatrix<f64>,
    b: &[f64],
    x: &mut [f64],
    opts: &BlockOptions,
    inner: Box<dyn InnerSolver + 'a>,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<BlockReport, BlockError> {
    log.info("BlockSolveInt", "---------------------------------------", 5);
    check_supported(params, opts, comm)?;
    if params.get("block preconditioner").is_none() {
        log.info("BlockSolveInt", "Using block preconditioning mode by default", 4);
    }
    let n = a.nrows;
    if b.len() < n || x.len() < n {
        return Err(BlockError::fatal(
            "BlockSolveInt",
            "right-hand side or solution shorter than the matrix",
        ));
    }
    let layout = BlockLayout::from_params(params, dofs)?;
    if layout.nblocks < 2 {
        // BlockSolveInt hands a single block to DefaultSolve under `outer:`.
        return Err(BlockError::Unsupported("a block linear system of one block".into()));
    }
    log.info(
        "BlockInitMatrix",
        &format!("Block matrix will be of size {}", layout.nblocks),
        6,
    );
    let eliminate_zero = params.logical("block eliminate zero submatrices") == Some(true);
    let sys = BlockSystem::split(a, layout, eliminate_zero)?;
    for i in 0..sys.nblocks() {
        if sys.block(i, i).is_none() {
            return Err(BlockError::fatal(
                "BlockKrylovIter",
                format!("Submatrix {} not associated!", 11 * (i + 1)),
            ));
        }
    }
    let mut xb = sys.gather(x);
    let bb = sys.gather(b);
    let report = solve_system(params, &sys, &bb, &mut xb, opts, inner, log, comm)?;
    sys.scatter(&xb, x);
    log.info("BlockSolveInt", "All done", 4);
    log.info("BlockSolveInt", "-------------------------------------------------", 5);
    Ok(report)
}

/// The block solve of a system assembled per block (`BlockSolver` without a
/// `Procedure`, MainUtils.F90:4323): `b` and `x` are block-ordered (the block
/// variables one after the other). Inner solves use `inner`.
#[allow(clippy::too_many_arguments)]
pub fn block_solve_system<'a>(
    params: &ParamMap,
    sys: &BlockSystem,
    b: &[f64],
    x: &mut [f64],
    opts: &BlockOptions,
    inner: Box<dyn InnerSolver + 'a>,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<BlockReport, BlockError> {
    check_supported(params, opts, comm)?;
    if sys.nblocks() < 2 {
        // BlockSolver hands a single block to DefaultSolve under `outer:`.
        return Err(BlockError::Unsupported("a block linear system of one block".into()));
    }
    if b.len() < sys.size() || x.len() < sys.size() {
        return Err(BlockError::fatal(
            "BlockSolver",
            "right-hand side or solution shorter than the block system",
        ));
    }
    solve_system(params, sys, b, x, opts, inner, log, comm)
}

/// Refusals shared by both entry points.
fn check_supported(params: &ParamMap, opts: &BlockOptions, comm: &dyn Communicator) -> Result<(), BlockError> {
    if comm.size() > 1 {
        return Err(BlockError::Unsupported("parallel block linear systems".into()));
    }
    refuse_unported(params)?;
    if params.get("block matrix schur variable").is_some() && opts.schur.is_none() {
        return Err(BlockError::fatal(
            "BlockPrecMatrix",
            "Schur matrix does not exist for: Block Matrix Schur Variable",
        ));
    }
    Ok(())
}

/// `Block Preconditioner` (default true) picks `BlockKrylovIter`, else
/// `BlockStandardIter`; both run under `outer:`.
#[allow(clippy::too_many_arguments)]
fn solve_system<'a>(
    params: &ParamMap,
    sys: &BlockSystem,
    b: &[f64],
    x: &mut [f64],
    opts: &BlockOptions,
    inner: Box<dyn InnerSolver + 'a>,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<BlockReport, BlockError> {
    if params.logical("block preconditioner") == Some(false) {
        log.info("BlockSolveInt", "Using block solution strategy", 6);
        standard(params, sys, b, x, opts, inner, log)
    } else {
        log.info("BlockSolveInt", "Using block preconditioning strategy", 6);
        krylov_iter(params, sys, b, x, opts, inner, log, comm)
    }
}

/// `BlockStandardIter` under `outer:`, with `block:` and `block NN:` for the
/// corrections.
fn standard<'a>(
    params: &ParamMap,
    sys: &BlockSystem,
    b: &[f64],
    x: &mut [f64],
    opts: &BlockOptions,
    mut inner: Box<dyn InnerSolver + 'a>,
    log: &dyn Logger,
) -> Result<BlockReport, BlockError> {
    let nb = sys.nblocks();
    let add = opts.additive_namespaces;
    let outer = resolve(params, &["outer:"], add);
    let order = block_order(&outer, nb)?;
    let views = |base: &ParamMap| -> Vec<ParamMap> {
        (0..nb)
            .map(|i| resolve(base, &["outer:", "block:", &block_namespace(i)], add))
            .collect()
    };
    let first = views(params);
    // BlockSolve.F90:4296: from the second iteration on the preconditioners
    // are not recomputed.
    let mut later = params.clone();
    later.insert("no precondition recompute", Value::Logical(true));
    let later = views(&later);
    // The tailored matrix replaces the block as is: BlockStandardIter does no
    // DoPrecScaling.
    let prec: Vec<Option<PrecBlock<'_>>> = prec_blocks(sys, opts)?
        .into_iter()
        .map(|p| p.map(|p| PrecBlock { unscale: None, ..p }))
        .collect();
    let start_norm = comp_norm(&x[..sys.size()]);
    let out = standard_iter(&outer, &order, [&first, &later], &prec, sys, b, x, &mut *inner, log)?;
    let solution_norm = comp_norm(&x[..sys.size()]);
    Ok(BlockReport {
        nblocks: nb,
        sweeps: out.stats.iterations,
        max_change: out.stats.residual,
        outer: Some(out.stats),
        inner: out.records,
        rhs_norm: comp_norm(&b[..sys.size()]),
        start_norm,
        solution_norm,
    })
}

/// `BlockKrylovIter`: the outer iteration under `outer:`, on block-ordered
/// vectors.
#[allow(clippy::too_many_arguments)]
fn krylov_iter<'a>(
    params: &ParamMap,
    sys: &BlockSystem,
    bb: &[f64],
    xb: &mut [f64],
    opts: &BlockOptions,
    inner: Box<dyn InnerSolver + 'a>,
    log: &dyn Logger,
    comm: &dyn Communicator,
) -> Result<BlockReport, BlockError> {
    let nb = sys.nblocks();
    let ndim = sys.size();
    let (bb, xb) = (&bb[..ndim], &mut xb[..ndim]);

    let rhs_norm = comp_norm(bb);
    log.info("BlockKrylovIter", &format!("Rhs norm at start: {rhs_norm:12.5E}"), 10);
    let mut report = BlockReport {
        nblocks: nb,
        outer: None,
        sweeps: 0,
        inner: Vec::new(),
        rhs_norm,
        start_norm: 0.0,
        solution_norm: 0.0,
        max_change: 0.0,
    };
    if rhs_norm < f64::EPSILON {
        log.warn(
            "BlockKrylovIter",
            "With zero'ish r.h.s. it does not make sense to call the solver!",
        );
        return Ok(report);
    }
    let prev = comp_norm(xb);
    report.start_norm = prev;
    log.info("BlockKrylovIter", &format!("Solution norm at start: {prev:12.5E}"), 10);

    let add = opts.additive_namespaces;
    let outer = resolve(params, &["outer:"], add);
    let settings = prec_settings(params, &outer, sys, opts, add)?;
    let mv_dummy = (0..nb).map(|i| dummy(&outer, i)).collect();

    let op = BlockOperator::new(sys, mv_dummy);
    let prec = BlockPreconditioner::new(sys, settings, inner);
    // The outer IterSolver gets the block operator and preconditioner; the
    // stationary methods read matrix entries and have none to read.
    let mut solver = iterative(&outer, None)?;
    let stats = solver.solve(&op, &prec, bb, xb, comm);
    if let Some(e) = prec.take_error() {
        return Err(e);
    }
    let stats = stats?;
    report.sweeps = prec.applications();
    report.inner = prec.records();
    check_converged(&outer, &stats.status, opts.global_abort)?;
    log.info("BlockKrylovIter", "Finished block system iteration", 18);

    let xnorm = comp_norm(xb);
    log.info("BlockKrylovIter", &format!("Solution norm: {xnorm:12.5E}"), 8);
    report.max_change = 2.0 * (xnorm - prev).abs() / (xnorm + prev);
    log.info(
        "BlockKrylovIter",
        &format!("Relative change: {:12.5E}", report.max_change),
        8,
    );
    report.solution_norm = xnorm;
    report.outer = Some(stats);
    Ok(report)
}

/// `Block Order` (1-based) as 0-based block indices, default `0..nb`.
fn block_order(outer: &ParamMap, nb: usize) -> Result<Vec<usize>, BlockError> {
    Ok(match outer.get("block order") {
        Some(Value::Integer(o)) => {
            if o.len() != nb {
                return Err(BlockError::fatal(
                    "BlockMatrixPrec",
                    format!("Block Order size should be {nb}"),
                ));
            }
            let mut seen = vec![false; nb];
            let mut order = Vec::with_capacity(nb);
            for &k in o {
                if k < 1 || k as usize > nb || seen[k as usize - 1] {
                    return Err(BlockError::fatal(
                        "BlockMatrixPrec",
                        "Block Order is not a permutation of the blocks",
                    ));
                }
                seen[k as usize - 1] = true;
                order.push(k as usize - 1);
            }
            order
        }
        Some(_) => {
            return Err(BlockError::fatal(
                "BlockMatrixPrec",
                "Block Order is not an integer array",
            ))
        }
        None => (0..nb).collect(),
    })
}

/// `CompNorm` (BlockSolve.F90:2853), serial: `sqrt(sum x^2 / n)`.
pub fn comp_norm(x: &[f64]) -> f64 {
    (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt()
}

/// `Dummy block i` (1-based in the keyword) in the view `p`.
fn dummy(p: &ParamMap, i: usize) -> bool {
    p.logical(&format!("dummy block {}", i + 1)) == Some(true)
}

/// The settings `BlockMatrixPrec` reads. `outer` is the section under
/// `outer:`, where the outer iteration calls the preconditioner.
fn prec_settings<'a>(
    params: &ParamMap,
    outer: &ParamMap,
    sys: &BlockSystem,
    opts: &BlockOptions<'a>,
    add: bool,
) -> Result<PrecSettings<'a>, BlockError> {
    let nb = sys.nblocks();
    let order = block_order(outer, nb)?;
    let gauss_seidel = outer.logical("block gauss-seidel") == Some(true);
    // BlockMatrixPrec.F90:3870-3876: with the outer scaling active, nested
    // scaling is eliminated unless asked otherwise.
    let no_nested_scaling = opts.outer_scaling.is_some() && outer.logical("eliminate nested scaling").unwrap_or(true);
    let prec = prec_blocks(sys, opts)?;
    // The ListAddLogical calls of BlockMatrixPrec put plain keys in the
    // section, below any namespaced spelling. `Linear System Skip Scaling`
    // is set around the solve of an ordinary block and reset after it, so a
    // block solved with a Schur matrix sees it false.
    let mut base = params.clone();
    base.insert("linear system skip complex", Value::Logical(false));
    base.insert("linear system skip loads", Value::Logical(true));
    let block_view = resolve(&base, &["outer:", "block:"], add);
    let dummy = (0..nb).map(|i| dummy(&block_view, i)).collect();
    let mut inner = Vec::with_capacity(nb);
    let mut passive = vec![vec![false; nb]; nb];
    for i in 0..nb {
        let ns = block_namespace(i);
        let mut base = base.clone();
        if no_nested_scaling {
            let skip = !prec[i].as_ref().is_some_and(|p| p.unscale.is_some());
            base.insert("linear system skip scaling", Value::Logical(skip));
        }
        let view = resolve(&base, &["outer:", "block:", &ns], add);
        for (k, row) in passive.iter_mut().enumerate() {
            let key = normalize_key(&format!("block gauss-seidel passive {}{}", k + 1, i + 1));
            row[i] = view.logical(&key) == Some(true);
        }
        inner.push(view);
    }
    Ok(PrecSettings {
        order,
        gauss_seidel,
        dummy,
        inner,
        passive,
        prec,
    })
}

/// The preconditioning matrix of each block: the Schur matrix for the last
/// one, with `DoPrecScaling` when the monolithic system is scaled.
fn prec_blocks<'a>(sys: &BlockSystem, opts: &BlockOptions<'a>) -> Result<Vec<Option<PrecBlock<'a>>>, BlockError> {
    let nb = sys.nblocks();
    let mut prec: Vec<Option<PrecBlock<'a>>> = vec![None; nb];
    let Some(schur) = opts.schur else { return Ok(prec) };
    let last = nb - 1;
    if schur.nrows == 0 {
        return Ok(prec);
    }
    if schur.nrows != sys.range(last).len() || schur.ncols != schur.nrows {
        return Err(BlockError::fatal(
            "BlockPrecMatrix",
            format!("Schur matrix of size {} does not fit block {}", schur.nrows, 11 * nb),
        ));
    }
    let unscale = match opts.outer_scaling {
        None => None,
        Some(sc) => {
            // `DiagScaling(k::l)`: component k of every node, k the block's
            // only component (`InvBlockStruct`).
            let Some(layout) = &sys.layout else {
                return Err(BlockError::Unsupported(
                    "a Schur block in the scaled system of a block-assembled system".into(),
                ));
            };
            let comps: Vec<usize> = (0..layout.comps).filter(|&c| layout.block_of[c] == last).collect();
            let [k] = comps[..] else {
                return Err(BlockError::fatal(
                    "BlockMatrixPrec",
                    format!("Cannot define the originating block {} for scaling!", nb),
                ));
            };
            let d: Vec<f64> = sc.diag.iter().skip(k).step_by(layout.comps).copied().collect();
            if d.len() != schur.nrows {
                return Err(BlockError::fatal(
                    "BlockMatrixPrec",
                    "DiagScaling does not fit the Schur matrix",
                ));
            }
            Some((d, sc.rhs_scaling))
        }
    };
    prec[last] = Some(PrecBlock { matrix: schur, unscale });
    Ok(prec)
}
