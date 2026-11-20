//! Read access to a variable the way the SaveData solvers see it: a nodal
//! field with an optional permutation, or one strided component of one.

use elmer_core::{FieldKind, NO_DOF};
use elmer_solver::{SolverError, VarHandle, VarId, Variable, VariableStore};

/// A nodal field, or one component of it. `Perm(node)` of the Fortran is
/// [`VarView::slot`]; a missing permutation entry is `None`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct VarView<'a> {
    pub values: &'a [f64],
    pub perm: Option<&'a [u32]>,
    /// Dofs of the field the values belong to.
    pub stride: usize,
    pub offset: usize,
    /// `Var % Dofs` as the operators see it: components of a vector field
    /// the view covers (1 for a component alias).
    pub ncomp: usize,
    /// An elemental (`-elem`) field: `perm` is indexed by element, not node.
    pub elemental: bool,
    /// A DG (`-dg`) field: `perm` is indexed by DG index, not node.
    pub dg: bool,
}

impl<'a> VarView<'a> {
    fn of(v: &'a Variable, offset: usize, ncomp: usize) -> Self {
        VarView {
            values: &v.values,
            perm: v.perm.as_deref(),
            stride: v.dofs,
            offset,
            ncomp,
            elemental: v.kind == FieldKind::Elemental,
            dg: v.kind == FieldKind::Dg,
        }
    }

    /// Size of the permutation, or the number of value slots without one
    /// (the `NoNodes` of `VectorStatistics`).
    pub fn n_perm(&self) -> usize {
        match self.perm {
            Some(p) => p.len(),
            None => self.values.len() / self.stride,
        }
    }

    /// Equation slot of `node` (0-based), `None` where the field has no value.
    pub fn slot(&self, node: usize) -> Option<usize> {
        match self.perm {
            Some(p) => p.get(node).copied().filter(|&s| s != NO_DOF).map(|s| s as usize),
            None => (node < self.n_perm()).then_some(node),
        }
    }

    /// Component `k` of the value in `slot`.
    pub fn at(&self, slot: usize, k: usize) -> f64 {
        self.values[slot * self.stride + self.offset + k]
    }

    /// The scalar the statistical operators use for `slot`: the value, or the
    /// Euclidean norm for a vector field.
    pub fn magnitude(&self, slot: usize) -> f64 {
        if self.ncomp <= 1 {
            self.at(slot, 0)
        } else {
            (0..self.ncomp).map(|k| self.at(slot, k).powi(2)).sum::<f64>().sqrt()
        }
    }
}

/// What a `Variable n` name resolves to.
pub(crate) enum Resolved<'a> {
    /// A one-value (or one-vector) global; saved as it is.
    Global(VarId),
    Field(VarView<'a>),
    /// `name` itself does not exist but `name 1`, `name 2`, .. do
    /// (`ComponentVar` of SaveScalars). The operators then see component 1
    /// only; the view is that component.
    Components { first: VarView<'a>, count: usize },
}

/// `VariableGet(Model % Variables, name)` with SaveScalars' fallback to the
/// component variables. `Ok(None)` when nothing of that name exists.
pub(crate) fn resolve<'a>(vars: &'a VariableStore, name: &str) -> Result<Option<Resolved<'a>>, SolverError> {
    let nodal = |v: &'a Variable, what: &str| -> Result<(), SolverError> {
        match v.kind {
            FieldKind::Nodal => Ok(()),
            k => Err(SolverError::Unsupported(format!("{what}: variable '{}' of kind {k:?}", v.name))),
        }
    };
    match vars.resolve(name) {
        Some(VarHandle::Primary(id)) => {
            let v = vars.by_id(id);
            if v.kind == FieldKind::Global {
                return Ok(Some(Resolved::Global(id)));
            }
            nodal(v, "SaveData")?;
            Ok(Some(Resolved::Field(VarView::of(v, 0, v.dofs))))
        }
        Some(VarHandle::Component { parent, component }) => {
            let v = vars.by_id(parent);
            nodal(v, "SaveData")?;
            Ok(Some(Resolved::Field(VarView::of(v, component, 1))))
        }
        None => {
            let first = format!("{name} 1");
            let Some(VarHandle::Component { parent, component }) = vars.resolve(&first) else {
                return Ok(None);
            };
            let v = vars.by_id(parent);
            nodal(v, "SaveData")?;
            let count = (2..=3)
                .take_while(|k| vars.resolve(&format!("{name} {k}")).is_some())
                .count()
                + 1;
            Ok(Some(Resolved::Components {
                first: VarView::of(v, component, 1),
                count,
            }))
        }
    }
}

/// [`resolve`] as `BulkIntegrals` reads a field: an elemental variable too,
/// whose value at a point of an element is the element's own
/// (`EvaluateVariableAtGivenPoint`).
pub(crate) fn resolve_bulk<'a>(vars: &'a VariableStore, name: &str) -> Result<Option<Resolved<'a>>, SolverError> {
    let elemental = |v: &Variable| v.kind == FieldKind::Elemental;
    match vars.resolve(name) {
        Some(VarHandle::Primary(id)) if elemental(vars.by_id(id)) => {
            let v = vars.by_id(id);
            Ok(Some(Resolved::Field(VarView::of(v, 0, v.dofs))))
        }
        Some(VarHandle::Component { parent, component }) if elemental(vars.by_id(parent)) => {
            Ok(Some(Resolved::Field(VarView::of(vars.by_id(parent), component, 1))))
        }
        _ => resolve(vars, name),
    }
}

/// [`resolve`] as `VectorStatistics` and `VectorMeanDeviation` read a
/// field: an elemental or Gauss-point variable loses its Perm
/// (SaveScalars.F90:1874, :2055), so the operator runs over every stored
/// value. A DG variable keeps its Perm, over its DG indexes. Other
/// variables resolve as [`resolve`] does.
pub(crate) fn resolve_values<'a>(vars: &'a VariableStore, name: &str) -> Result<Option<Resolved<'a>>, SolverError> {
    let per_element = |v: &Variable| matches!(v.kind, FieldKind::Elemental | FieldKind::Ip);
    let unpermuted = |v: &'a Variable, offset: usize, ncomp: usize| {
        Ok(Some(Resolved::Field(VarView { perm: None, ..VarView::of(v, offset, ncomp) })))
    };
    match vars.resolve(name) {
        Some(VarHandle::Primary(id)) if per_element(vars.by_id(id)) => {
            let v = vars.by_id(id);
            unpermuted(v, 0, v.dofs)
        }
        Some(VarHandle::Component { parent, component }) if per_element(vars.by_id(parent)) => {
            unpermuted(vars.by_id(parent), component, 1)
        }
        Some(VarHandle::Primary(id)) if vars.by_id(id).kind == FieldKind::Dg => {
            let v = vars.by_id(id);
            Ok(Some(Resolved::Field(VarView::of(v, 0, v.dofs))))
        }
        Some(VarHandle::Component { parent, component }) if vars.by_id(parent).kind == FieldKind::Dg => {
            Ok(Some(Resolved::Field(VarView::of(vars.by_id(parent), component, 1))))
        }
        _ => resolve(vars, name),
    }
}

/// `ComponentName(base, i)` for a field with `dofs > 1`, 1-based.
pub(crate) fn component_name(base: &str, i: usize) -> String {
    format!("{base} {i}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    /// An elemental field's statistics run over its stored values, not
    /// through its Perm; a nodal field keeps its Perm.
    #[test]
    fn per_element_fields_lose_their_perm_for_statistics() {
        let mut vars = VariableStore::new();
        let perm: Arc<[u32]> = Arc::from(vec![1u32, 0, NO_DOF]);
        vars.add(Variable::new("testvar a", FieldKind::Elemental, 1, Some(perm.clone()), vec![3.0, 5.0])).unwrap();
        vars.add(Variable::new("temperature", FieldKind::Nodal, 1, Some(perm), vec![3.0, 5.0])).unwrap();
        assert!(matches!(resolve(&vars, "testvar a"), Err(SolverError::Unsupported(_))));
        let Ok(Some(Resolved::Field(v))) = resolve_values(&vars, "testvar a") else { panic!() };
        assert_eq!((v.perm, v.n_perm(), v.slot(1)), (None, 2, Some(1)));
        let Ok(Some(Resolved::Field(t))) = resolve_values(&vars, "temperature") else { panic!() };
        assert_eq!((t.n_perm(), t.slot(2)), (3, None));
    }
}
