//! Keyword values that are not plain constants (tables, MATC, procedures),
//! evaluated with `elmer-expr` against the simulation's fields.

use elmer_core::{MeshView, Value};
use elmer_expr::{compile, FieldLookup, MatcEngine};

use crate::error::SolverError;
use crate::variable::VariableStore;

/// The fields of a simulation as `elmer-expr` reads them: a variable or
/// component alias at a node through its perm, a global's first value,
/// mesh coordinates and the `time` global.
pub struct SimFields<'a> {
    pub vars: &'a VariableStore,
    pub mesh: &'a dyn MeshView,
}

impl FieldLookup for SimFields<'_> {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        let v = self.vars.strided(var)?;
        match v.perm {
            Some(perm) => {
                let s = *perm.get(node as usize)?;
                if s == elmer_core::NO_DOF {
                    None
                } else {
                    v.values.get(s as usize * v.stride + v.offset).copied()
                }
            }
            // Globals and perm-less fields: one value (globals) or node-indexed.
            None if v.values.len() == 1 => Some(v.values[0]),
            None => v.values.get(node as usize * v.stride + v.offset).copied(),
        }
    }

    fn coord(&self, node: u32) -> [f64; 3] {
        self.mesh.coord(node)
    }

    fn time(&self) -> f64 {
        self.vars.scalar("time").unwrap_or(0.0)
    }
}

/// Stand-in until a MATC binding is installed in [`crate::Backends::matc`]:
/// every MATC value fails.
pub struct NoMatc;

impl MatcEngine for NoMatc {
    fn eval(&self, expr: &str, _args: &[f64]) -> Result<Vec<f64>, String> {
        Err(format!("no MATC engine installed for: {expr}"))
    }
}

/// `ListGetReal` of a scalar keyword at `nodes`.
pub fn eval_at_nodes(
    v: &Value,
    nodes: &[u32],
    fields: &SimFields<'_>,
    matc: Option<&dyn MatcEngine>,
) -> Result<Vec<f64>, SolverError> {
    if let Value::Real { data, .. } = v {
        if data.len() == 1 {
            return Ok(vec![data[0]; nodes.len()]);
        }
    }
    let c = compile(v);
    if c.ncomp() != 1 {
        return Err(SolverError::Unsupported(format!("{}x{} value where a scalar is expected", c.rows(), c.cols())));
    }
    let mut out = vec![0.0; nodes.len()];
    c.eval_nodes(nodes, fields, matc.unwrap_or(&NoMatc), &mut out)?;
    Ok(out)
}

/// `ListGetConstReal` of a scalar keyword.
pub fn eval_const(v: &Value, matc: Option<&dyn MatcEngine>) -> Result<f64, SolverError> {
    if let Value::Real { data, .. } = v {
        if data.len() == 1 {
            return Ok(data[0]);
        }
    }
    Ok(compile(v).eval_const(matc.unwrap_or(&NoMatc), [0.0; 3])?)
}
