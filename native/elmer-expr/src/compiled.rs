//! `Value` -> [`Compiled`] property and its batched evaluation.

use elmer_core::{Dependent, DependentKind, Value};

use crate::curve;
use crate::udf::{Udf, UdfRegistry};
use crate::{FieldLookup, MatcEngine, MISSING};

#[derive(Debug, thiserror::Error)]
pub enum ExprError {
    #[error("unsupported property: {0}")]
    Unsupported(String),
    #[error("invalid property: {0}")]
    Invalid(String),
    #[error("MATC evaluation failed: {0}")]
    Matc(String),
    #[error("procedure {library:?} {entry:?} is not registered")]
    UnknownProcedure { library: String, entry: String },
    #[error("evaluation buffer: {0}")]
    Shape(String),
}

/// Knobs that the `elmer_core::Value` representation does not carry.
#[derive(Clone, Copy, Debug, Default)]
pub struct CompileOptions {
    /// `Real Monotone Cubic`: limit the spline slopes so the curve stays
    /// monotone between table points (`CubicSpline(..., monotone=.TRUE.)`).
    pub monotone: bool,
}

/// Where one dependency value comes from.
#[derive(Clone, Debug, PartialEq)]
enum Source {
    Time,
    Coord(usize),
    /// `coordinate`: three values, x y z.
    CoordAll,
    Literal(f64),
    Field(String),
}

impl Source {
    fn width(&self) -> usize {
        if matches!(self, Source::CoordAll) {
            3
        } else {
            1
        }
    }
}

/// `ListParseStrToVars` for one `Dependent::variables` list. Entries may
/// themselves be comma lists.
fn parse_sources(variables: &[String]) -> Vec<Source> {
    let mut out = Vec::new();
    for entry in variables {
        for name in entry.split(',') {
            let name = name.trim().to_ascii_lowercase();
            if name.is_empty() {
                continue;
            }
            out.push(match name.as_str() {
                "time" => Source::Time,
                "coordinate" => Source::CoordAll,
                "coordinate 1" => Source::Coord(0),
                "coordinate 2" => Source::Coord(1),
                "coordinate 3" => Source::Coord(2),
                _ => {
                    let numeric = name.chars().all(|c| "-.0123456789eE".contains(c));
                    match name.parse::<f64>() {
                        Ok(v) if numeric => Source::Literal(v),
                        _ => Source::Field(name),
                    }
                }
            });
        }
    }
    out
}

#[derive(Clone, Debug)]
struct Table {
    x: Vec<f64>,
    /// One series per tensor component, row-major component order.
    comps: Vec<Vec<f64>>,
    /// Spline slopes; scalar tables only.
    cubic: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
enum Body {
    Const(Vec<f64>),
    Table(Table),
    Matc { expr: String },
    Udf(Udf),
    Invalid(ExprErrorKind),
}

/// Clonable shadow of the compile-time error kept inside an invalid property.
#[derive(Clone, Debug)]
enum ExprErrorKind {
    Unsupported(String),
    Invalid(String),
    UnknownProcedure { library: String, entry: String },
}

impl ExprErrorKind {
    fn to_error(&self) -> ExprError {
        match self {
            ExprErrorKind::Unsupported(s) => ExprError::Unsupported(s.clone()),
            ExprErrorKind::Invalid(s) => ExprError::Invalid(s.clone()),
            ExprErrorKind::UnknownProcedure { library, entry } => ExprError::UnknownProcedure {
                library: library.clone(),
                entry: entry.clone(),
            },
        }
    }
}

impl From<ExprError> for ExprErrorKind {
    fn from(e: ExprError) -> Self {
        match e {
            ExprError::Unsupported(s) => ExprErrorKind::Unsupported(s),
            ExprError::UnknownProcedure { library, entry } => {
                ExprErrorKind::UnknownProcedure { library, entry }
            }
            other => ExprErrorKind::Invalid(other.to_string()),
        }
    }
}

/// A property ready for repeated evaluation. Scalars are 1x1.
#[derive(Clone, Debug)]
pub struct Compiled {
    rows: usize,
    cols: usize,
    args: Vec<Source>,
    body: Body,
}

/// Compiles a keyword value against the process-wide [`UdfRegistry`]. A value
/// that cannot be evaluated (a string, an unregistered procedure) compiles to
/// a property whose evaluation returns the error; use [`compile_with`] to get
/// it up front.
pub fn compile(v: &Value) -> Compiled {
    match compile_with(v, UdfRegistry::global(), &CompileOptions::default()) {
        Ok(c) => c,
        Err(e) => Compiled {
            rows: 1,
            cols: 1,
            args: Vec::new(),
            body: Body::Invalid(e.into()),
        },
    }
}

pub fn compile_with(
    v: &Value,
    registry: &UdfRegistry,
    opts: &CompileOptions,
) -> Result<Compiled, ExprError> {
    match v {
        Value::Real { rows, cols, data } => {
            if data.len() != rows * cols || data.is_empty() {
                return Err(ExprError::Invalid(format!(
                    "real {rows}x{cols} value holds {} numbers",
                    data.len()
                )));
            }
            Ok(Compiled {
                rows: *rows,
                cols: *cols,
                args: Vec::new(),
                body: Body::Const(data.clone()),
            })
        }
        Value::Dependent(d) => compile_dependent(d, registry, opts),
        other => Err(ExprError::Unsupported(format!(
            "{} is not a real-valued property",
            match other {
                Value::Logical(_) => "a logical",
                Value::Integer(_) => "an integer",
                Value::Str(_) => "a string",
                Value::File(_) => "a file name",
                _ => "this value",
            }
        ))),
    }
}

fn compile_dependent(
    d: &Dependent,
    registry: &UdfRegistry,
    opts: &CompileOptions,
) -> Result<Compiled, ExprError> {
    let (rows, cols) = (d.rows.max(1), d.cols.max(1));
    let ncomp = rows * cols;
    let args = parse_sources(&d.variables);
    let body = match &d.kind {
        DependentKind::Table { x, y, cubic } => {
            if args.is_empty() {
                return Err(ExprError::Invalid("table without a dependency variable".into()));
            }
            let n = x.len();
            if n == 0 || y.len() != n * ncomp {
                return Err(ExprError::Invalid(format!(
                    "table of {n} points x {ncomp} components holds {} values",
                    y.len()
                )));
            }
            if ncomp == 1 && !curve::is_strictly_increasing(x) {
                return Err(ExprError::Invalid(
                    "table abscissae are not monotonically increasing".into(),
                ));
            }
            // y is point-major, row-major within a point: y[p * ncomp + c].
            let comps: Vec<Vec<f64>> = (0..ncomp)
                .map(|c| (0..n).map(|p| y[p * ncomp + c]).collect())
                .collect();
            // ListAddDepReal: splines only for scalar tables with n > 3.
            let cubic = (*cubic && ncomp == 1 && n > 3)
                .then(|| curve::cubic_spline(x, &comps[0], opts.monotone));
            Body::Table(Table {
                x: x.clone(),
                comps,
                cubic,
            })
        }
        DependentKind::Matc(expr) => Body::Matc { expr: expr.clone() },
        DependentKind::Procedure { library, entry } => {
            let udf = registry.lookup(library, entry).ok_or_else(|| {
                ExprError::UnknownProcedure {
                    library: library.clone(),
                    entry: entry.clone(),
                }
            })?;
            match (&udf, args.is_empty()) {
                (Udf::Const(_), false) => {
                    return Err(ExprError::Invalid(format!(
                        "procedure {library:?} {entry:?} is a constant function but has variables"
                    )))
                }
                (Udf::Array(_), _) => {}
                (_, _) if ncomp > 1 => {
                    return Err(ExprError::Invalid(format!(
                        "procedure {library:?} {entry:?} returns a scalar for a {rows}x{cols} property"
                    )))
                }
                _ => {}
            }
            Body::Udf(udf)
        }
        DependentKind::Lua(_) => {
            return Err(ExprError::Unsupported("Lua property functions".into()))
        }
    };
    Ok(Compiled {
        rows,
        cols,
        args,
        body,
    })
}

impl Compiled {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Number of values per evaluation point (`rows * cols`).
    pub fn ncomp(&self) -> usize {
        self.rows * self.cols
    }

    /// The stored numbers of a constant property (`ListGetConstReal` without
    /// side effects), row-major.
    pub fn constant(&self) -> Option<&[f64]> {
        match &self.body {
            Body::Const(v) => Some(v),
            _ => None,
        }
    }

    /// `ListGetConstReal`: a property that does not depend on fields. A
    /// constant MATC expression is evaluated without arguments and a constant
    /// user function at `xyz`. Returns the first component.
    pub fn eval_const(&self, m: &dyn MatcEngine, xyz: [f64; 3]) -> Result<f64, ExprError> {
        match &self.body {
            Body::Const(v) => Ok(v[0]),
            Body::Matc { expr } if self.args.is_empty() => first(m.eval(expr, &[]), expr),
            Body::Udf(Udf::Const(f)) => Ok(f(xyz)),
            Body::Invalid(e) => Err(e.to_error()),
            _ => Err(ExprError::Invalid("constant requested of a field-dependent property".into())),
        }
    }

    /// `ListGetReal` / `ListGetRealArray` at nodes: evaluates the property at
    /// each of `nodes`, writing `ncomp()` values per node, node-major and
    /// row-major within a node. `out.len()` must equal `nodes.len() * ncomp()`.
    ///
    /// Dependency values are fetched per variable for the whole batch. Nodes
    /// where one is [`MISSING`] get zeros.
    pub fn eval_nodes(
        &self,
        nodes: &[u32],
        f: &dyn FieldLookup,
        m: &dyn MatcEngine,
        out: &mut [f64],
    ) -> Result<(), ExprError> {
        let nc = self.ncomp();
        if out.len() != nodes.len() * nc {
            return Err(ExprError::Shape(format!(
                "out holds {} values, {} nodes x {nc} components expected",
                out.len(),
                nodes.len()
            )));
        }
        match &self.body {
            Body::Invalid(e) => return Err(e.to_error()),
            Body::Const(v) => {
                for cell in out.chunks_exact_mut(nc) {
                    cell.copy_from_slice(v);
                }
                return Ok(());
            }
            _ => {}
        }

        let nn = nodes.len();
        let gathered = self.gather(nodes, f);
        let narg = gathered.len() / nn.max(1);
        let mut t = vec![0.0_f64; narg];
        for (i, cell) in out.chunks_exact_mut(nc).enumerate() {
            for (a, ta) in t.iter_mut().enumerate() {
                *ta = gathered[a * nn + i];
            }
            if t.contains(&MISSING) {
                cell.fill(0.0);
                continue;
            }
            let xyz = if matches!(&self.body, Body::Udf(Udf::Const(_))) {
                f.coord(nodes[i])
            } else {
                [0.0; 3]
            };
            self.apply(&t, nodes[i], xyz, m, cell)?;
        }
        Ok(())
    }

    /// `ListGetElementReal` at integration points: dependency values are taken
    /// at the element `nodes`, interpolated with `basis` (`n_ip` rows of
    /// `nodes.len()` shape-function values, row-major) and the property is
    /// then evaluated once per point. `out.len()` must be `n_ip * ncomp()`,
    /// point-major. A point whose nodal dependency is [`MISSING`] gets zeros.
    pub fn eval_ips(
        &self,
        nodes: &[u32],
        basis: &[f64],
        f: &dyn FieldLookup,
        m: &dyn MatcEngine,
        out: &mut [f64],
    ) -> Result<(), ExprError> {
        let nn = nodes.len();
        let nc = self.ncomp();
        if nn == 0 || !basis.len().is_multiple_of(nn) {
            return Err(ExprError::Shape(format!(
                "basis holds {} values for {nn} nodes",
                basis.len()
            )));
        }
        let n_ip = basis.len() / nn;
        if out.len() != n_ip * nc {
            return Err(ExprError::Shape(format!(
                "out holds {} values, {n_ip} points x {nc} components expected",
                out.len()
            )));
        }
        match &self.body {
            Body::Invalid(e) => return Err(e.to_error()),
            Body::Const(v) => {
                for cell in out.chunks_exact_mut(nc) {
                    cell.copy_from_slice(v);
                }
                return Ok(());
            }
            _ => {}
        }

        let gathered = self.gather(nodes, f);
        let narg = gathered.len() / nn;
        let missing = gathered.contains(&MISSING);
        let mut t = vec![0.0_f64; narg];
        for (ip, cell) in out.chunks_exact_mut(nc).enumerate() {
            let b = &basis[ip * nn..(ip + 1) * nn];
            if missing {
                cell.fill(0.0);
                continue;
            }
            for (a, ta) in t.iter_mut().enumerate() {
                // SUM(Basis(1:ni) * ParValues(j,1:ni)), accumulated in order.
                *ta = b
                    .iter()
                    .zip(&gathered[a * nn..(a + 1) * nn])
                    .fold(0.0, |s, (w, v)| s + w * v);
            }
            let mut xyz = [0.0_f64; 3];
            if matches!(&self.body, Body::Udf(Udf::Const(_))) {
                for (k, node) in nodes.iter().enumerate() {
                    let c = f.coord(*node);
                    for d in 0..3 {
                        xyz[d] += b[k] * c[d];
                    }
                }
            }
            self.apply(&t, u32::MAX, xyz, m, cell)?;
        }
        Ok(())
    }

    /// Number of dependency values passed to the function.
    /// Values each source passes: a multi-dof field passes all its
    /// components (`f.dofs`).
    fn widths(&self, f: &dyn FieldLookup) -> Vec<usize> {
        self.effective_sources()
            .iter()
            .map(|s| match s {
                Source::Field(name) => f.dofs(name).max(1),
                s => s.width(),
            })
            .collect()
    }

    /// A constant MATC expression is evaluated with `Time` bound (`'st'`).
    fn effective_sources(&self) -> std::borrow::Cow<'_, [Source]> {
        if self.args.is_empty() && matches!(self.body, Body::Matc { .. }) {
            std::borrow::Cow::Owned(vec![Source::Time])
        } else {
            std::borrow::Cow::Borrowed(&self.args)
        }
    }

    /// Dependency values for every node, argument-major: argument `a` of node
    /// `i` is `result[a * nodes.len() + i]`.
    fn gather(&self, nodes: &[u32], f: &dyn FieldLookup) -> Vec<f64> {
        let nn = nodes.len();
        let sources = self.effective_sources();
        let widths = self.widths(f);
        let mut g = vec![0.0_f64; widths.iter().sum::<usize>() * nn];
        let mut a = 0;
        for (s, &w) in sources.iter().zip(&widths) {
            match s {
                Source::Time => g[a * nn..(a + 1) * nn].fill(f.time()),
                Source::Literal(v) => g[a * nn..(a + 1) * nn].fill(*v),
                Source::Coord(k) => {
                    for (i, &n) in nodes.iter().enumerate() {
                        g[a * nn + i] = f.coord(n)[*k];
                    }
                }
                Source::CoordAll => {
                    for (i, &n) in nodes.iter().enumerate() {
                        let c = f.coord(n);
                        for k in 0..3 {
                            g[(a + k) * nn + i] = c[k];
                        }
                    }
                }
                Source::Field(name) if w == 1 => f.values(name, nodes, &mut g[a * nn..(a + 1) * nn]),
                Source::Field(name) => {
                    for k in 0..w {
                        f.component_values(name, k, nodes, &mut g[(a + k) * nn..(a + k + 1) * nn]);
                    }
                }
            }
            a += w;
        }
        g
    }

    /// Evaluates the non-constant body for one point with dependency values `t`.
    fn apply(
        &self,
        t: &[f64],
        node: u32,
        xyz: [f64; 3],
        m: &dyn MatcEngine,
        cell: &mut [f64],
    ) -> Result<(), ExprError> {
        match &self.body {
            Body::Table(tab) => {
                let t1 = t[0];
                for (c, v) in cell.iter_mut().enumerate() {
                    *v = curve::interpolate_curve(&tab.x, &tab.comps[c], t1, tab.cubic.as_deref());
                }
            }
            Body::Matc { expr } => {
                let var = if self.args.is_empty() { "st" } else { "tx" };
                let r = m
                    .eval_with(var, expr, t)
                    .map_err(ExprError::Matc)?;
                if r.len() < cell.len() {
                    return Err(ExprError::Matc(format!(
                        "{expr:?} returned {} values, {} expected",
                        r.len(),
                        cell.len()
                    )));
                }
                cell.copy_from_slice(&r[..cell.len()]);
            }
            Body::Udf(Udf::Const(g)) => {
                cell[0] = g(xyz);
            }
            Body::Udf(Udf::Node(g)) => cell[0] = g(node, t),
            Body::Udf(Udf::Array(g)) => g(node, t, cell),
            Body::Const(_) | Body::Invalid(_) => unreachable!("handled by the callers"),
        }
        Ok(())
    }
}

fn first(r: Result<Vec<f64>, String>, expr: &str) -> Result<f64, ExprError> {
    let v = r.map_err(ExprError::Matc)?;
    v.first()
        .copied()
        .ok_or_else(|| ExprError::Matc(format!("{expr:?} returned no values")))
}
