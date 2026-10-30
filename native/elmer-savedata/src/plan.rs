//! Frozen configuration of the two hooks: the `Solver n` section is read
//! once when the solver is built, and every keyword the hooks do not
//! implement is refused there with [`SolverError::Unsupported`], never
//! ignored.

use elmer_core::{normalize_key, ParamMap, Params, Value};
use elmer_solver::SolverError;

fn unsupported(msg: impl Into<String>) -> SolverError {
    SolverError::Unsupported(msg.into())
}

fn logical_true(p: &ParamMap, key: &str) -> bool {
    p.logical(key) == Some(true)
}

/// Refuses the section when any of `keys` is set (to anything but logical
/// false).
fn refuse_keys(p: &ParamMap, hook: &str, keys: &[(&str, &str)]) -> Result<(), SolverError> {
    for (key, what) in keys {
        if p.get(key).is_some_and(|v| *v != Value::Logical(false)) {
            return Err(unsupported(format!("{hook}: {what} ({key})")));
        }
    }
    Ok(())
}

/// The strings keywords hold have no quotes in a [`ParamMap`], but a section
/// built by hand may keep them.
fn clean_name(s: &str) -> String {
    normalize_key(s.trim().trim_matches('"'))
}

// ------------------------------------------------------------ SaveScalars

/// What an operator computes, after the `body `/`body force `/`material `
/// prefix is taken off.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum OperKind {
    /// Needs no variable: counts of the mesh and the run.
    Mesh,
    /// Read off the variable's bookkeeping (`norm`, `dofs`, `nans`, ...).
    Meta,
    /// `VectorStatistics`.
    Statistic,
    /// `VectorMeanDeviation`.
    Deviation,
    /// `BulkIntegrals`.
    Bulk,
    /// `BoundaryStatistics`.
    BoundaryStatistic,
    /// `BoundaryIntegrals` (`boundary int`, `boundary int mean`, `area`,
    /// `diffusive flux`).
    BoundaryIntegral,
}

/// Which section a masked operator reads its logical from.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum MaskKind {
    None,
    Body,
    BodyForce,
    Material,
}

/// An `Operator n` value split the way SaveScalars.F90:581-612 does.
#[derive(Clone, Debug)]
pub(crate) struct ParsedOper {
    /// The full text, as in the column name.
    pub oper0: String,
    /// Without the mask prefix.
    pub oper: String,
    pub mask: MaskKind,
    pub kind: OperKind,
}

const STATISTICS: [&str; 13] = [
    "sum",
    "sum abs",
    "mean abs",
    "max",
    "max abs",
    "min",
    "min abs",
    "mean",
    "variance",
    "range",
    "sum square",
    "mean square",
    "rms",
];

const BULK: [&str; 11] = [
    "int",
    "int mean",
    "int square",
    "int square mean",
    "int rms",
    "int abs",
    "int abs mean",
    "int variance",
    "volume",
    "potential energy",
    "diffusive energy",
];

const BOUNDARY_STATISTICS: [&str; 7] = [
    "boundary sum",
    "boundary dofs",
    "boundary max",
    "boundary max abs",
    "boundary min",
    "boundary min abs",
    "boundary mean",
];

/// Operators that exist in SaveScalars.F90 and are refused here, with the
/// reason.
fn refused_operator(oper: &str) -> Option<&'static str> {
    Some(match oper {
        "convective flux" => "flux reconstruction needs the parent element's gradient",
        "convective energy" => "SaveScalars.F90:2358 evaluates it from a gradient it never computes",
        "cpu time" | "wall time" | "cpu memory" | "threads" => "run-dependent value, not reproducible",
        _ => return None,
    })
}

/// Splits and classifies an operator. Fails on anything not implemented.
pub(crate) fn parse_operator(text: &str) -> Result<ParsedOper, SolverError> {
    let oper0 = clean_name(text);
    let (mask, oper) = if let Some(r) = oper0.strip_prefix("body force ") {
        (MaskKind::BodyForce, r)
    } else if let Some(r) = oper0.strip_prefix("body ") {
        (MaskKind::Body, r)
    } else if let Some(r) = oper0.strip_prefix("material ") {
        (MaskKind::Material, r)
    } else {
        (MaskKind::None, oper0.as_str())
    };
    let oper = oper.to_string();
    if oper.starts_with("positive ") || oper.starts_with("negative ") {
        // The Fortran shifts the name inside itself (SaveScalars.F90:612) and
        // leaves the tail of the old text behind, so the case never matches.
        return Err(unsupported(format!(
            "SaveScalars operator '{oper0}': positive/negative variants"
        )));
    }
    if let Some(why) = refused_operator(&oper) {
        return Err(unsupported(format!("SaveScalars operator '{oper0}': {why}")));
    }
    let kind = match oper.as_str() {
        "partitions" | "partition checksum" | "partition neighbours checksum" | "nodes" | "elements"
        | "bounding box" => OperKind::Mesh,
        "norm" | "nonlin change" | "steady state change" | "nonlin iter" | "nonlin converged"
        | "steady converged" | "dofs" | "nans" => OperKind::Meta,
        "deviation" => OperKind::Deviation,
        "area" | "boundary int" | "boundary int mean" | "diffusive flux" => OperKind::BoundaryIntegral,
        o if STATISTICS.contains(&o) => OperKind::Statistic,
        o if BULK.contains(&o) => OperKind::Bulk,
        o if BOUNDARY_STATISTICS.contains(&o) => OperKind::BoundaryStatistic,
        _ => return Err(unsupported(format!("SaveScalars operator '{oper0}'"))),
    };
    if mask != MaskKind::None && !matches!(kind, OperKind::Statistic | OperKind::Deviation | OperKind::Bulk) {
        return Err(unsupported(format!(
            "SaveScalars operator '{oper0}': a mask prefix only applies to statistics and integrals"
        )));
    }
    Ok(ParsedOper {
        oper0,
        oper,
        mask,
        kind,
    })
}

/// `Show Norm Index` / `Show Norm Name` / `Reference Values`: the scalar the
/// solver reports as its norm (SaveScalars.F90:3434).
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum PseudoNorm {
    Index(usize),
    Name(String),
    Reference(Vec<f64>),
}

/// Frozen `SaveScalars` configuration.
#[derive(Clone, Debug)]
pub struct ScalarPlan {
    /// The solver section; the `Variable n` / `Operator n` walk reads it at
    /// run time because it depends on which variables exist then.
    pub(crate) section: ParamMap,
    pub(crate) filename: Option<String>,
    pub(crate) file_append: bool,
    pub(crate) echo_values: bool,
    pub(crate) line_marker: i32,
    pub(crate) live_graph: bool,
    pub(crate) partition_numbering: bool,
    pub(crate) results_prefix: String,
    pub(crate) save_points: Vec<i32>,
    /// `Save Coordinates`, padded to three coordinates.
    pub(crate) save_coordinates: Vec<[f64; 3]>,
    pub(crate) pseudo_norm: Option<PseudoNorm>,
    pub(crate) cost_function_index: Option<usize>,
    pub(crate) equation: Option<String>,
    pub(crate) needs_integration: bool,
    pub(crate) operators: Vec<ParsedOper>,
}

fn parse_rows(p: &ParamMap, key: &str) -> Result<Option<(usize, usize, Vec<f64>)>, SolverError> {
    match p.get(key) {
        Some(Value::Real { rows, cols, data }) if rows * cols == data.len() => {
            Ok(Some((*rows, *cols, data.clone())))
        }
        Some(Value::Real { rows, cols, data }) => Err(unsupported(format!(
            "{key}: a {rows}x{cols} array holding {} values",
            data.len()
        ))),
        Some(_) => Err(unsupported(format!("{key} must be a real array"))),
        None => Ok(None),
    }
}

impl ScalarPlan {
    /// Reads and checks a `Procedure = "SaveData" "SaveScalars"` section.
    pub fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        const H: &str = "SaveScalars";
        refuse_keys(
            p,
            H,
            &[
                ("polyline coordinates", "integrals over polylines"),
                ("exact coordinates", "values at exact coordinates inside elements"),
                ("save eigenvalues", "eigenvalue columns"),
                ("save eigen values", "eigenvalue columns"),
                ("save eigenfrequencies", "eigenfrequency columns"),
                ("save eigen frequencies", "eigenfrequency columns"),
                ("simulation timing", "run-dependent timing columns"),
                ("save component results", "results of the Components section"),
                ("json format", "JSON side file"),
                ("filename numbering", "numbered file names"),
                ("run control numbering", "run-control numbered file names"),
                ("mark failed strategy", "failed-strategy marker files"),
                ("avoid failed strategy", "failed-strategy marker files"),
                ("output partitioning directory", "per-partition output directory"),
                ("expression 1", "free-form MATC expressions"),
            ],
        )?;
        for key in p.0.keys() {
            if key.starts_with("coefficient ") {
                return Err(unsupported(format!(
                    "{H}: material coefficients ({key}); integrals are unweighted"
                )));
            }
            if key.starts_with("expression ") {
                return Err(unsupported(format!("{H}: free-form MATC expressions ({key})")));
            }
        }

        let mut operators = Vec::new();
        for (key, v) in &p.0 {
            let Some(n) = key.strip_prefix("operator ") else { continue };
            if n.parse::<u32>().is_err() {
                continue;
            }
            let text = match v {
                Value::Str(s) | Value::File(s) => s,
                _ => return Err(unsupported(format!("{H}: {key} must be a string"))),
            };
            operators.push(parse_operator(text)?);
        }
        let needs_integration = operators
            .iter()
            .any(|o| matches!(o.kind, OperKind::Bulk | OperKind::BoundaryIntegral));

        let filename = p.string("filename").map(|s| s.trim().to_string());
        let save_to_file = filename.is_some();
        let echo_values = p.logical("echo values").unwrap_or(!save_to_file);

        let save_points = match p.get("save points") {
            Some(Value::Integer(v)) => v.clone(),
            Some(_) => return Err(unsupported(format!("{H}: Save Points must be an integer array"))),
            None => Vec::new(),
        };
        let save_coordinates = match parse_rows(p, "save coordinates")? {
            Some((rows, cols, data)) => {
                if cols == 0 || cols > 3 {
                    return Err(unsupported(format!("{H}: Save Coordinates with {cols} columns")));
                }
                (0..rows)
                    .map(|r| {
                        let mut c = [0.0; 3];
                        c[..cols].copy_from_slice(&data[r * cols..(r + 1) * cols]);
                        c
                    })
                    .collect()
            }
            None => Vec::new(),
        };

        let pseudo_norm = if let Some(n) = p.string("show norm name") {
            Some(PseudoNorm::Name(clean_name(n)))
        } else if let Some(i) = p.integer("show norm index") {
            (i > 0).then_some(PseudoNorm::Index(i as usize))
        } else {
            None
        };
        let pseudo_norm = match (pseudo_norm, parse_rows(p, "reference values")?) {
            (Some(n), _) => Some(n),
            (None, Some((rows, cols, data))) => {
                if cols != 1 {
                    return Err(unsupported(format!("{H}: Reference Values with {cols} columns")));
                }
                Some(PseudoNorm::Reference(data[..rows].to_vec()))
            }
            (None, None) => None,
        };

        let cost = p.integer("cost function index").filter(|&i| i > 0).map(|i| i as usize);

        Ok(ScalarPlan {
            section: p.clone(),
            filename,
            file_append: logical_true(p, "file append"),
            echo_values,
            line_marker: p.integer("line marker").unwrap_or(0),
            live_graph: p.logical("live graph").or_else(|| p.logical("cvs format")).unwrap_or(false),
            partition_numbering: logical_true(p, "partition numbering"),
            results_prefix: p.string("scalars prefix").unwrap_or("res:").to_string(),
            save_points,
            save_coordinates,
            pseudo_norm,
            cost_function_index: cost,
            equation: p.string("equation").map(str::to_string),
            needs_integration,
            operators,
        })
    }

    /// The `Operator n` values of the section, normalised, in keyword order.
    pub fn operators(&self) -> impl Iterator<Item = &str> {
        self.operators.iter().map(|o| o.oper0.as_str())
    }

    /// Whether an operator needs the element basis and quadrature (element
    /// integrals), so a front end knows to call [`crate::install_fe`].
    pub fn needs_integration(&self) -> bool {
        self.needs_integration
    }

    /// Whether the section makes the solver report a norm, which gives it a
    /// global variable (`SaveScalars_init`, SaveScalars.F90:63).
    pub fn reports_norm(&self) -> bool {
        self.pseudo_norm.is_some()
    }
}

// --------------------------------------------------------------- SaveLine

/// Frozen `SaveLine` configuration.
#[derive(Clone, Debug)]
pub struct LinePlan {
    pub(crate) section: ParamMap,
    pub(crate) filename: String,
    pub(crate) file_append: bool,
    pub(crate) line_marker: i32,
    pub(crate) skip_boundary_info: bool,
    /// `Save Mask`, normalised.
    pub(crate) mask: String,
    pub(crate) optimize_order: bool,
    /// `Variable 1`, `Variable 2`, .. in order.
    pub(crate) variables: Vec<String>,
}

impl LinePlan {
    /// Reads and checks a `Procedure = "SaveData" "SaveLine"` section.
    pub fn from_params(p: &ParamMap) -> Result<Self, SolverError> {
        const H: &str = "SaveLine";
        refuse_keys(
            p,
            H,
            &[
                ("polyline coordinates", "polylines"),
                ("circle coordinates", "circles"),
                ("save isocurves", "isocurves"),
                ("save axis", "axis lines"),
                ("save axis 1", "axis lines"),
                ("save axis 2", "axis lines"),
                ("save axis 3", "axis lines"),
                ("save heat flux", "boundary flux reconstruction"),
                ("save flux", "boundary flux reconstruction"),
                ("intersect edge", "edge intersections"),
                ("intersect coordinate", "edge intersections"),
                ("break line loop", "line loop breaking"),
                ("calculate weights", "nodal weights"),
                ("filename numbering", "numbered file names"),
                ("save solver mesh index", "another solver's mesh"),
                ("show norm index", "the line pseudo norm"),
                ("parallel reduce", "MPI reductions"),
                ("enforce parallel mode", "MPI reductions"),
                ("output partitioning directory", "per-partition output directory"),
            ],
        )?;
        let mut variables = Vec::new();
        for i in 1..=99 {
            match p.string(&format!("variable {i}")) {
                Some(s) => variables.push(clean_name(s)),
                None => break,
            }
        }
        Ok(LinePlan {
            section: p.clone(),
            filename: p.string("filename").map_or_else(|| "sides.dat".into(), |s| s.trim().to_string()),
            file_append: logical_true(p, "file append"),
            line_marker: p.integer("line marker").unwrap_or(0),
            skip_boundary_info: logical_true(p, "skip boundary info"),
            mask: p.string("save mask").map_or_else(|| "save line".into(), clean_name),
            optimize_order: p.logical("optimize node ordering").unwrap_or(true),
            variables,
        })
    }
}
