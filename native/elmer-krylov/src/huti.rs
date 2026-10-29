//! The HUTI driver: what IterSolve does around one `huti_[dz]*solv` call.

use std::marker::PhantomData;

use elmer_core::{
    Communicator, IdentityPrecond, IterativeSolver, LinalgError, LinearOperator, Preconditioner,
    SolveStats, SolveStatus,
};

use crate::config::{HutiConfig, HutiMethod};
use crate::methods;
use crate::output;
use crate::scalar::HutiScalar;

/// `HUTI_INFO` return codes (huti_fdefs.h:20-59).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum HutiInfo {
    Convergence = 1,
    MaxIter = 2,
    Divergence = 3,
    QmrRhoPsi = 10,
    QmrDelta = 11,
    QmrEpsilon = 12,
    QmrBeta = 13,
    QmrGamma = 14,
    CgRho = 20,
    CgsRho = 25,
    TfqmrRho = 30,
    BiCGStabRho = 35,
    BiCGStabOmega = 37,
    GmresAlpha = 40,
    GmresBeta = 41,
    BiCGStab2Rho = 45,
}

impl HutiInfo {
    pub fn code(self) -> i32 {
        self as i32
    }

    fn status(self) -> SolveStatus {
        match self {
            HutiInfo::Convergence => SolveStatus::Converged,
            HutiInfo::MaxIter => SolveStatus::MaxIterations,
            HutiInfo::Divergence => SolveStatus::Diverged,
            other => SolveStatus::Breakdown(format!("{other:?} (HUTI_INFO {})", other.code())),
        }
    }
}

/// Receives each HUTI stdout line (without the newline).
pub type OutputSink = Box<dyn FnMut(&str) + Send>;

/// HUTI Krylov solver. Configure with [`HutiConfig`] (usually
/// [`HutiConfig::from_params`]) and call through [`IterativeSolver`].
///
/// The initial guess is always the caller's `x` (IterSolve passes
/// `HUTI_USERSUPPLIEDX`). Iteration lines go to stdout unless a sink is set
/// with [`Huti::set_output`].
pub struct Huti<T: HutiScalar = f64> {
    pub cfg: HutiConfig,
    sink: Option<OutputSink>,
    history: Vec<f64>,
    info: Option<HutiInfo>,
    _scalar: PhantomData<T>,
}

impl<T: HutiScalar> Huti<T> {
    pub fn new(cfg: HutiConfig) -> Self {
        Huti {
            cfg,
            sink: None,
            history: Vec::new(),
            info: None,
            _scalar: PhantomData,
        }
    }

    /// Send iteration lines here instead of stdout.
    pub fn set_output(&mut self, sink: OutputSink) {
        self.sink = Some(sink);
    }

    /// The monitored residual of the last solve, one entry per iteration at
    /// the point where HUTI would print it (whether or not it printed).
    pub fn history(&self) -> &[f64] {
        &self.history
    }

    /// `HUTI_INFO` of the last solve.
    pub fn info(&self) -> Option<HutiInfo> {
        self.info
    }
}

/// The operator, the two preconditioner slots and the reductions one HUTI
/// call receives (`matvecsubr`, `pcondlsubr`, `pcondrsubr`, `dotprodfun`,
/// `normfun`). Argument order follows Rust: `(input, output)`.
pub(crate) struct Ops<'a, T: HutiScalar> {
    a: &'a dyn LinearOperator<T>,
    pl: &'a dyn Preconditioner<T>,
    pr: &'a dyn Preconditioner<T>,
    comm: &'a dyn Communicator,
    /// Unconjugated dot product: complex CG (IterSolve.F90:1508-1520).
    bilinear: bool,
}

impl<T: HutiScalar> Ops<'_, T> {
    pub(crate) fn mv(&self, x: &[T], y: &mut [T]) {
        self.a.apply(x, y);
    }

    pub(crate) fn mv_t(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        self.a.apply_transpose(x, y)
    }

    pub(crate) fn pl(&self, r: &[T], z: &mut [T]) {
        self.pl.apply(r, z);
    }

    pub(crate) fn pr(&self, r: &[T], z: &mut [T]) {
        self.pr.apply(r, z);
    }

    pub(crate) fn pl_t(&self, r: &[T], z: &mut [T]) -> Result<(), LinalgError> {
        self.pl.apply_transpose(r, z)
    }

    pub(crate) fn pr_t(&self, r: &[T], z: &mut [T]) -> Result<(), LinalgError> {
        self.pr.apply_transpose(r, z)
    }

    /// `dotprodfun(x, y)`: `sum conj(x_i) y_i`, or `sum x_i y_i` for complex CG.
    pub(crate) fn dot(&self, x: &[T], y: &[T]) -> T {
        let s = if self.bilinear {
            T::dotu(x, y)
        } else {
            T::dotc(x, y)
        };
        if self.comm.size() == 1 {
            return s;
        }
        let mut p = s.parts();
        self.comm.sum_f64(&mut p);
        T::from_parts(p)
    }

    /// `normfun(x)`: BLAS `dnrm2`/`dznrm2` on one rank, the square root of the
    /// global sum of squares otherwise.
    pub(crate) fn norm(&self, x: &[T]) -> f64 {
        if self.comm.size() == 1 {
            return T::nrm2(x);
        }
        let mut s = [T::sumsq(x)];
        self.comm.sum_f64(&mut s);
        s[0].sqrt()
    }
}

/// Residual history and the `HUTI_DBUGLVL` printing.
pub(crate) struct Monitor<'s> {
    every: usize,
    sink: &'s mut dyn FnMut(&str),
    history: &'s mut Vec<f64>,
}

impl Monitor<'_> {
    /// Records `res` for iteration `it`; returns whether HUTI prints there
    /// (`mod(iter_count, HUTI_DBUGLVL) == 0`).
    pub(crate) fn record(&mut self, it: usize, res: f64) -> bool {
        self.history.push(res);
        self.every != 0 && it.is_multiple_of(self.every)
    }

    /// Whether any output is on (`HUTI_DBUGLVL /= HUTI_NO_DEBUG`).
    pub(crate) fn on(&self) -> bool {
        self.every != 0
    }

    pub(crate) fn emit(&mut self, line: &str) {
        (self.sink)(line);
    }

    /// The standard line, when due.
    pub(crate) fn step(&mut self, it: usize, res: f64) {
        if self.record(it, res) {
            self.emit(&output::iteration_line(it, res));
        }
    }

    /// The unconditional line at the `1000 continue` exit.
    pub(crate) fn finish(&mut self, it: usize, res: f64) {
        if self.on() {
            self.emit(&output::iteration_line(it, res));
        }
    }
}

/// What a HUTI method hands back: `HUTI_INFO`, `HUTI_ITERS` and the last
/// monitored residual.
pub(crate) struct Outcome {
    pub info: HutiInfo,
    pub iters: usize,
    pub residual: f64,
}

impl<T: HutiScalar> IterativeSolver<T> for Huti<T> {
    fn solve(
        &mut self,
        a: &dyn LinearOperator<T>,
        m: &dyn Preconditioner<T>,
        b: &[T],
        x: &mut [T],
        comm: &dyn Communicator,
    ) -> Result<SolveStats, LinalgError> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n || x.len() != n {
            return Err(LinalgError::Dim(format!(
                "operator {}x{}, b {}, x {}",
                a.nrows(),
                a.ncols(),
                n,
                x.len()
            )));
        }
        let cfg = self.cfg.clone();
        if cfg.perturb_zero_guess && x.iter().all(|v| v.is_zero()) {
            x.fill(T::from_parts([1.0e-8, 1.0e-8]));
        }

        let ident = IdentityPrecond;
        let (pl, pr): (&dyn Preconditioner<T>, &dyn Preconditioner<T>) = if cfg.left_precond {
            (m, &ident)
        } else {
            (&ident, m)
        };
        let ops = Ops {
            a,
            pl,
            pr,
            comm,
            bilinear: T::COMPLEX && cfg.method == HutiMethod::Cg,
        };

        let mut stdout = |line: &str| println!("{line}");
        let sink: &mut dyn FnMut(&str) = match self.sink.as_mut() {
            Some(s) => s.as_mut(),
            None => &mut stdout,
        };
        self.history.clear();
        self.info = None;
        let mut mon = Monitor {
            every: if comm.rank() == 0 { cfg.verbosity } else { 0 },
            sink,
            history: &mut self.history,
        };

        let out = match cfg.method {
            HutiMethod::Cg => methods::cg::solve(&ops, &cfg, &mut mon, x, b),
            HutiMethod::Cgs => methods::cgs::solve(&ops, &cfg, &mut mon, x, b),
            HutiMethod::BiCGStab => methods::bicgstab::solve(&ops, &cfg, &mut mon, x, b),
            HutiMethod::BiCGStab2 => methods::bicgstab2::solve(&ops, &cfg, &mut mon, x, b),
            HutiMethod::Gmres => methods::gmres::solve(&ops, &cfg, &mut mon, x, b),
            HutiMethod::Qmr => methods::qmr::solve(&ops, &cfg, &mut mon, x, b)?,
            HutiMethod::Tfqmr => methods::tfqmr::solve(&ops, &cfg, &mut mon, x, b),
        };
        self.info = Some(out.info);
        Ok(SolveStats {
            iterations: out.iters,
            residual: out.residual,
            status: out.info.status(),
        })
    }
}
