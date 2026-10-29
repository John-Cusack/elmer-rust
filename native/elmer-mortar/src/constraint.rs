//! The constraint matrix of a solver's mortar projectors
//! (`GenerateConstraintMatrix`, ProjectorUtils.F90).

use elmer_core::{Logger, Params, NO_DOF};

use crate::projector::Projector;
use crate::MortarError;

const CALLER: &str = "GenerateConstraintMatrix";

/// `A % ConstraintMatrix`: one row per constrained node of each projector,
/// columns in the solver's matrix numbering (0-based), in the order the
/// Fortran stores them (not sorted).
#[derive(Clone, Debug, PartialEq)]
pub struct ConstraintMatrix {
    /// Length `n_rows + 1`.
    pub rows: Vec<usize>,
    pub cols: Vec<u32>,
    pub values: Vec<f64>,
    pub rhs: Vec<f64>,
    /// `InvPerm`, 0-based: `rowoffset + slot`, where `slot` is the
    /// constrained node's matrix row and `rowoffset` grows by the matrix size
    /// for each projector added before this one.
    pub inv_perm: Vec<usize>,
}

impl ConstraintMatrix {
    pub fn n_rows(&self) -> usize {
        self.rows.len() - 1
    }
}

/// The solver-side inputs of [`generate_constraint_matrix`].
pub struct ConstraintTarget<'a> {
    /// The solver section.
    pub solver: &'a dyn Params,
    /// The variable name (`Solver % Variable % Name`), for the per-component
    /// keywords such as `Mortar BC Temperature`.
    pub name: &'a str,
    pub dofs: usize,
    /// Parent node -> matrix row (`Solver % Variable % Perm`), [`NO_DOF`]
    /// when absent.
    pub perm: &'a [u32],
    /// `Solver % Matrix % NumberOfRows`.
    pub n_rows: usize,
}

fn flag(p: &dyn Params, k: &str) -> bool {
    p.logical(k).unwrap_or(false)
}

/// `GenerateConstraintMatrix` for the solver's mortar projectors
/// (`projectors[i]` from boundary condition `i + 1`, as
/// [`crate::generate_mortar_projectors`] returns them). Ported: a scalar
/// field, each projector its own block of rows (no `Mortar BCs Additive`),
/// `Minimum Projector Value`, nodal (`[D - P]`) and Galerkin projectors,
/// and `Mortar BC <var> = False` to skip a boundary. Returns `None` when
/// there is nothing to constrain.
pub fn generate_constraint_matrix(
    bcs: &[&dyn Params],
    projectors: &[Option<Projector>],
    t: &ConstraintTarget<'_>,
    log: &dyn Logger,
) -> Result<Option<ConstraintMatrix>, MortarError> {
    let s = t.solver;
    let bcount = projectors.iter().flatten().count();
    log.info(CALLER, &format!("Number of mortar matrices: {bcount}"), 12);
    if projectors.iter().flatten().all(|p| p.n_rows() == 0) {
        log.info(
            CALLER,
            "Nothing to do since there are no constrained dofs!",
            12,
        );
        return Ok(None);
    }
    if s.get("mortar diag").is_some() {
        return Err(MortarError::Unsupported("Mortar Diag".into()));
    }
    if flag(s, "skip already constrained dofs") {
        return Err(MortarError::Unsupported(
            "Skip Already Constrained Dofs".into(),
        ));
    }
    let additive = match s.logical("mortar bcs additive") {
        Some(a) => a,
        None => {
            (bcount > 1 && flag(s, "eliminate linear constraints"))
                || flag(s, "apply conforming bcs")
        }
    };
    if additive {
        return Err(MortarError::Unsupported(
            "Mortar BCs Additive (summed projectors)".into(),
        ));
    }
    if t.dofs != 1 {
        return Err(MortarError::Unsupported(format!(
            "mortar constraints on a {}-component field",
            t.dofs
        )));
    }
    if bcs.iter().any(|bc| bc.get("projector priority").is_some()) {
        log.warn(CALLER, "Priority has effect only in additive mode!");
    }
    let eps = s.real("minimum projector value").unwrap_or(1.0e-8);

    let mut c = ConstraintMatrix {
        rows: vec![0],
        cols: Vec::new(),
        values: Vec::new(),
        rhs: Vec::new(),
        inv_perm: Vec::new(),
    };
    let mut rowoffset = 0usize;
    for (i, p) in projectors.iter().enumerate().rev() {
        let Some(p) = p else { continue };
        let bc = *bcs.get(i).ok_or_else(|| {
            MortarError::Unsupported(format!(
                "mortar projector of body force {}",
                i + 1 - bcs.len()
            ))
        })?;
        if flag(bc, "integral bc") {
            return Err(MortarError::Unsupported("Integral BC".into()));
        }
        for k in ["mortar bc coefficient", "mortar bc resistivity"] {
            if bc.get(&format!("{k} {}", t.name)).is_some() {
                return Err(MortarError::Unsupported(format!("{k} (mortar jump)")));
            }
        }
        log.info(
            CALLER,
            &format!("Adding mortar projector for BC: {}", i + 1),
            8,
        );
        log.info(
            CALLER,
            &format!("Adding projector rows: {}", p.n_rows()),
            12,
        );
        // A skipped boundary does not advance the row offset either.
        if !bc.logical(&format!("mortar bc {}", t.name)).unwrap_or(true) {
            continue;
        }
        add_scalar(&mut c, p, t, rowoffset, eps)?;
        rowoffset += t.n_rows;
    }
    if c.cols.is_empty() {
        log.info(CALLER, "No entries in constraint matrix!", 8);
        return Ok(None);
    }
    log.info(
        CALLER,
        &format!("Used {} rows and {} nonzeros", c.n_rows(), c.cols.len()),
        7,
    );
    Ok(Some(c))
}

/// `AddScalarConstraint` without summing, diagonal, rhs or activity mask.
fn add_scalar(
    c: &mut ConstraintMatrix,
    p: &Projector,
    t: &ConstraintTarget<'_>,
    rowoffset: usize,
    eps: f64,
) -> Result<(), MortarError> {
    // MasterScale and SlaveScale keep their defaults (1).
    let create_self = !p.galerkin;
    for i in 0..p.n_rows() {
        let r = p.rows[i]..p.rows[i + 1];
        if r.is_empty() {
            continue;
        }
        let k = p.inv_perm[i];
        if k == NO_DOF {
            continue;
        }
        let kk = t.perm[k as usize];
        if kk == NO_DOF {
            continue;
        }
        c.inv_perm.push(rowoffset + kk as usize);
        let mut valsum = 0.0;
        for v in &p.values[r.clone()] {
            valsum += v.abs();
        }
        let mut wsum = 0.0;
        for (col, val) in p.row(i) {
            if val.abs() < eps * valsum {
                continue;
            }
            let col2 = *t.perm.get(col as usize).ok_or_else(|| {
                MortarError::Fatal(format!(
                    "{CALLER}: col index too large: {} vs {}",
                    col + 1,
                    t.perm.len()
                ))
            })?;
            if col2 == NO_DOF {
                continue;
            }
            wsum += val;
            c.cols.push(col2);
            c.values.push(val);
        }
        if create_self {
            c.cols.push(kk);
            c.values.push(wsum);
        }
        c.rhs.push(0.0);
        c.rows.push(c.cols.len());
    }
    Ok(())
}
