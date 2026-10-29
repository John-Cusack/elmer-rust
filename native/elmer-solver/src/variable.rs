//! The field registry: Elmer `Variable_t` (Types.F90:760) and its list in
//! `Mesh % Variables`, with time and nonlinear histories and component
//! aliases as strided views.

use std::collections::BTreeMap;
use std::sync::Arc;

use elmer_core::{normalize_key, FieldKind, FieldRef, FieldSource, NO_DOF};

use crate::error::SolverError;

/// Index of a primary variable in a [`VariableStore`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VarId(pub usize);

/// Convergence bookkeeping Elmer keeps on the variable (Types.F90:766-768).
/// The `-1` start values mean "not computed yet".
#[derive(Clone, Debug, PartialEq)]
pub struct Convergence {
    pub norm: f64,
    pub prev_norm: f64,
    pub nonlin_change: f64,
    pub steady_change: f64,
    /// 0 not converged, 1 converged, 2 diverged, 3 exit condition; -1 unset.
    pub nonlin_converged: i32,
    pub steady_converged: i32,
    pub nonlin_iter: i32,
}

impl Default for Convergence {
    fn default() -> Self {
        Self {
            norm: 0.0,
            prev_norm: 0.0,
            nonlin_change: 0.0,
            steady_change: 0.0,
            nonlin_converged: -1,
            steady_converged: -1,
            nonlin_iter: -1,
        }
    }
}

/// One primary field. `values` is interleaved: slot `s`, component `k` is
/// `values[s * dofs + k]`; `perm[node]` gives `s` or [`NO_DOF`].
#[derive(Clone, Debug)]
pub struct Variable {
    /// Normalised name (lowercase, single spaces), as Elmer stores it.
    pub name: String,
    pub kind: FieldKind,
    pub dofs: usize,
    /// Shared between the variables of one solver, as Elmer shares `Perm`.
    pub perm: Option<Arc<[u32]>>,
    pub values: Vec<f64>,
    /// Time history, `prev_values[k]` = Elmer `PrevValues(:,k+1)`.
    pub prev_values: Vec<Vec<f64>>,
    /// Iterate before the current linear solve (`NonlinValues`).
    pub nonlin_values: Option<Vec<f64>>,
    /// Iterate before the current coupled-system pass (`SteadyValues`).
    pub steady_values: Option<Vec<f64>>,
    /// Unrelaxed change of the last iterate, for Aitken (`DeltaValues`).
    pub delta_values: Option<Vec<f64>>,
    pub conv: Convergence,
    /// Written by output writers (`-nooutput` clears it).
    pub output: bool,
    pub secondary: bool,
    pub values_changed: bool,
}

impl Variable {
    pub fn new(name: &str, kind: FieldKind, dofs: usize, perm: Option<Arc<[u32]>>, values: Vec<f64>) -> Self {
        Self {
            name: normalize_key(name),
            kind,
            dofs,
            perm,
            values,
            prev_values: Vec::new(),
            nonlin_values: None,
            steady_values: None,
            delta_values: None,
            conv: Convergence::default(),
            output: true,
            secondary: false,
            values_changed: true,
        }
    }

    /// A one-value global variable (`time`, `coupled iter`, ...).
    pub fn global(name: &str, value: f64) -> Self {
        Self::new(name, FieldKind::Global, 1, None, vec![value])
    }

    pub fn field_ref(&self) -> FieldRef<'_> {
        FieldRef {
            name: &self.name,
            kind: self.kind,
            dofs: self.dofs,
            perm: self.perm.as_deref(),
            values: &self.values,
        }
    }

    pub fn field_mut(&mut self) -> FieldMut<'_> {
        FieldMut {
            name: &self.name,
            kind: self.kind,
            dofs: self.dofs,
            perm: self.perm.as_deref(),
            values: &mut self.values,
        }
    }

    /// The slot of `node`, if the variable has one there.
    pub fn slot(&self, node: u32) -> Option<usize> {
        match &self.perm {
            Some(p) => p.get(node as usize).copied().filter(|&s| s != NO_DOF).map(|s| s as usize),
            None => Some(node as usize),
        }
    }
}

/// Mutable counterpart of [`FieldRef`].
#[derive(Debug)]
pub struct FieldMut<'a> {
    pub name: &'a str,
    pub kind: FieldKind,
    pub dofs: usize,
    pub perm: Option<&'a [u32]>,
    pub values: &'a mut [f64],
}

/// Component `offset` of an interleaved field: `values[i * stride + offset]`.
#[derive(Copy, Clone, Debug)]
pub struct Strided<'a> {
    pub values: &'a [f64],
    pub offset: usize,
    pub stride: usize,
    pub perm: Option<&'a [u32]>,
}

impl Strided<'_> {
    pub fn len(&self) -> usize {
        (self.values.len() + self.stride - 1).saturating_sub(self.offset) / self.stride
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get(&self, i: usize) -> f64 {
        self.values[i * self.stride + self.offset]
    }

    pub fn iter(&self) -> impl Iterator<Item = f64> + '_ {
        self.values.iter().skip(self.offset).step_by(self.stride).copied()
    }
}

/// Mutable component view.
#[derive(Debug)]
pub struct StridedMut<'a> {
    pub values: &'a mut [f64],
    pub offset: usize,
    pub stride: usize,
}

impl StridedMut<'_> {
    pub fn len(&self) -> usize {
        (self.values.len() + self.stride - 1).saturating_sub(self.offset) / self.stride
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get(&self, i: usize) -> f64 {
        self.values[i * self.stride + self.offset]
    }

    pub fn set(&mut self, i: usize, v: f64) {
        self.values[i * self.stride + self.offset] = v;
    }
}

/// What a name resolves to.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VarHandle {
    Primary(VarId),
    /// `component` is 0-based.
    Component { parent: VarId, component: usize },
}

/// The variable list. Names resolve case- and spacing-insensitively. A
/// variable with `dofs > 1` also registers `"<name> 1"`..`"<name> dofs"`
/// (GeneralUtils.F90:1478 ComponentNameStr), unless given explicit component
/// names through [`VariableStore::add_with_components`].
#[derive(Clone, Debug, Default)]
pub struct VariableStore {
    vars: Vec<Variable>,
    names: BTreeMap<String, VarHandle>,
    /// Registration order of names, as Elmer's linked list (for writers).
    order: Vec<VarId>,
}

impl VariableStore {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds `var` and its default component aliases. Fails on a name clash.
    pub fn add(&mut self, var: Variable) -> Result<VarId, SolverError> {
        let comps: Vec<String> = if var.dofs > 1 {
            (1..=var.dofs).map(|k| format!("{} {k}", var.name)).collect()
        } else {
            Vec::new()
        };
        self.add_with_components(var, &comps)
    }

    /// Adds `var` with explicit component names (Elmer's
    /// `Variable = Flow[Velocity:2 Pressure:1]` gives per-component names;
    /// a multi-dof component is one name per dof here).
    pub fn add_with_components(&mut self, var: Variable, components: &[String]) -> Result<VarId, SolverError> {
        let id = VarId(self.vars.len());
        let mut keys = vec![var.name.clone()];
        keys.extend(components.iter().map(|c| normalize_key(c)));
        for k in &keys {
            if self.names.contains_key(k) {
                return Err(SolverError::fatal("VariableAdd", format!("variable already exists: {k}")));
            }
        }
        self.names.insert(keys[0].clone(), VarHandle::Primary(id));
        for (c, k) in keys[1..].iter().enumerate() {
            self.names.insert(k.clone(), VarHandle::Component { parent: id, component: c });
        }
        self.vars.push(var);
        self.order.push(id);
        Ok(id)
    }

    /// The component names of `id`, in dof order (`ComponentName(Var % Name,
    /// k)`): the declared ones, or `"<name> k"`; empty for a scalar.
    pub fn component_names(&self, id: VarId) -> Vec<String> {
        let mut v: Vec<(usize, &String)> = self
            .names
            .iter()
            .filter_map(|(k, h)| match *h {
                VarHandle::Component { parent, component } if parent == id => Some((component, k)),
                _ => None,
            })
            .collect();
        v.sort_unstable();
        v.into_iter().map(|(_, k)| k.clone()).collect()
    }

    pub fn resolve(&self, name: &str) -> Option<VarHandle> {
        self.names.get(&normalize_key(name)).copied()
    }

    pub fn id(&self, name: &str) -> Option<VarId> {
        match self.resolve(name)? {
            VarHandle::Primary(id) => Some(id),
            VarHandle::Component { .. } => None,
        }
    }

    pub fn get(&self, name: &str) -> Option<&Variable> {
        self.id(name).map(|id| &self.vars[id.0])
    }

    pub fn get_mut(&mut self, name: &str) -> Option<&mut Variable> {
        let id = self.id(name)?;
        Some(&mut self.vars[id.0])
    }

    pub fn by_id(&self, id: VarId) -> &Variable {
        &self.vars[id.0]
    }

    pub fn by_id_mut(&mut self, id: VarId) -> &mut Variable {
        &mut self.vars[id.0]
    }

    /// A primary variable, or one component of one as a strided view.
    pub fn strided(&self, name: &str) -> Option<Strided<'_>> {
        let (v, offset, stride) = match self.resolve(name)? {
            VarHandle::Primary(id) => (&self.vars[id.0], 0, 1),
            VarHandle::Component { parent, component } => {
                let v = &self.vars[parent.0];
                (v, component, v.dofs)
            }
        };
        Some(Strided {
            values: &v.values,
            offset,
            stride,
            perm: v.perm.as_deref(),
        })
    }

    pub fn strided_mut(&mut self, name: &str) -> Option<StridedMut<'_>> {
        let (id, offset, stride) = match self.resolve(name)? {
            VarHandle::Primary(id) => (id, 0, 1),
            VarHandle::Component { parent, component } => (parent, component, self.vars[parent.0].dofs),
        };
        Some(StridedMut {
            values: &mut self.vars[id.0].values,
            offset,
            stride,
        })
    }

    /// The first value of a global variable (`time`, `timestep`, ...).
    pub fn scalar(&self, name: &str) -> Option<f64> {
        self.get(name).and_then(|v| v.values.first().copied())
    }

    /// Sets a global's first value, adding the global if missing.
    pub fn set_scalar(&mut self, name: &str, value: f64) {
        match self.get_mut(name) {
            Some(v) => v.values[0] = value,
            None => {
                let _ = self.add(Variable::global(name, value));
            }
        }
    }

    /// The same names and current values without histories or convergence
    /// state ([`crate::LiveModel`]).
    pub(crate) fn values_snapshot(&self) -> VariableStore {
        let vars = self
            .vars
            .iter()
            .map(|v| {
                let mut c = Variable::new(&v.name, v.kind, v.dofs, v.perm.clone(), v.values.clone());
                c.output = v.output;
                c
            })
            .collect();
        VariableStore {
            vars,
            names: self.names.clone(),
            order: self.order.clone(),
        }
    }

    /// Every name in Elmer's variable-list order: each primary variable in
    /// registration order followed by its component aliases.
    pub fn names_in_order(&self) -> Vec<(String, VarHandle)> {
        let mut out = Vec::new();
        for &id in &self.order {
            out.push((self.vars[id.0].name.clone(), VarHandle::Primary(id)));
            let mut comps: Vec<(usize, &String)> = self
                .names
                .iter()
                .filter_map(|(n, h)| match *h {
                    VarHandle::Component { parent, component } if parent == id => Some((component, n)),
                    _ => None,
                })
                .collect();
            comps.sort();
            out.extend(comps.into_iter().map(|(c, n)| (n.clone(), VarHandle::Component { parent: id, component: c })));
        }
        out
    }

    /// Primary variables in registration order.
    pub fn iter(&self) -> impl Iterator<Item = (VarId, &Variable)> + '_ {
        self.order.iter().map(|&id| (id, &self.vars[id.0]))
    }

    pub fn len(&self) -> usize {
        self.vars.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vars.is_empty()
    }
}

/// A snapshot of the store for writers: every primary variable with
/// `output` set, plus the time and step globals.
pub struct Snapshot<'a> {
    pub store: &'a VariableStore,
}

impl FieldSource for Snapshot<'_> {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        self.store
            .iter()
            .filter(|(_, v)| v.output)
            .map(|(_, v)| v.field_ref())
            .collect()
    }

    fn time(&self) -> f64 {
        self.store.scalar("time").unwrap_or(0.0)
    }

    fn step(&self) -> usize {
        self.store.scalar("timestep").map(|t| t as usize).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn components_are_strided_aliases() {
        let mut s = VariableStore::new();
        s.add(Variable::new("Displacement", FieldKind::Nodal, 2, None, vec![1., 2., 3., 4., 5., 6.]))
            .unwrap();
        assert_eq!(
            s.resolve("displacement  2"),
            Some(VarHandle::Component { parent: VarId(0), component: 1 })
        );
        let c = s.strided("Displacement 2").unwrap();
        assert_eq!(c.iter().collect::<Vec<_>>(), vec![2., 4., 6.]);
        assert_eq!(c.len(), 3);
        s.strided_mut("displacement 1").unwrap().set(2, 50.);
        assert_eq!(s.get("DISPLACEMENT").unwrap().values[4], 50.);
        assert!(s.add(Variable::global("displacement 1", 0.)).is_err());
    }

    #[test]
    fn globals() {
        let mut s = VariableStore::new();
        s.set_scalar("Coupled Iter", 1.0);
        s.set_scalar("coupled iter", 3.0);
        assert_eq!(s.scalar("coupled iter"), Some(3.0));
        assert_eq!(s.len(), 1);
    }
}
