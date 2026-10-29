//! The data the interpolation moves: scalar variables on nodes or on
//! nodes-on-elements (DG), with optional permutation and history columns.

use elmer_core::{ElemId, NO_DOF};

/// `Variable % TYPE` values that `LegitInterpVar` accepts.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Location {
    /// `Variable_on_nodes`.
    Nodes,
    /// `Variable_on_nodes_on_elements`: one slot per element-node (DG).
    NodesOnElements,
}

/// One scalar `Variable_t`. Vector variables are carried as one `FieldData`
/// per component, as Elmer does for its sub-variables.
#[derive(Clone, Debug)]
pub struct FieldData {
    pub name: String,
    pub location: Location,
    pub values: Vec<f64>,
    /// Node (or DG index) to slot; [`NO_DOF`] marks "not in the system".
    /// `None` is the identity (`Unitperm`).
    pub perm: Option<Vec<u32>>,
    /// `PrevValues(:, j)`, each column indexed by slot like `values`.
    pub prev_values: Vec<Vec<f64>>,
}

impl FieldData {
    pub fn nodal(name: &str, values: Vec<f64>) -> Self {
        FieldData {
            name: name.to_ascii_lowercase(),
            location: Location::Nodes,
            values,
            perm: None,
            prev_values: Vec::new(),
        }
    }

    pub fn with_perm(mut self, perm: Vec<u32>) -> Self {
        self.perm = Some(perm);
        self
    }

    pub fn with_prev(mut self, prev: Vec<Vec<f64>>) -> Self {
        self.prev_values = prev;
        self
    }

    pub fn with_location(mut self, location: Location) -> Self {
        self.location = location;
        self
    }

    /// Slot of node/DG index `i`, or `None` when it is not in the system.
    pub fn slot(&self, i: u32) -> Option<usize> {
        match &self.perm {
            None => Some(i as usize),
            Some(p) => match p.get(i as usize) {
                Some(&s) if s != NO_DOF => Some(s as usize),
                _ => None,
            },
        }
    }

    /// `LegitInterpVar`: scalar, nodal, not a coordinate, not a global.
    pub fn is_legit(&self) -> bool {
        if self.name.starts_with("coordinate") {
            return false;
        }
        !(self.perm.is_none() && self.values.len() == 1)
    }
}

/// `Element % DGIndexes` for every element, bulk first then boundary. An
/// empty list means the element has no DG indexes.
#[derive(Clone, Debug, Default)]
pub struct DgIndexes {
    lists: Vec<Vec<u32>>,
}

impl DgIndexes {
    pub fn new(lists: Vec<Vec<u32>>) -> Self {
        DgIndexes { lists }
    }

    pub fn get(&self, e: ElemId) -> &[u32] {
        self.lists.get(e.0 as usize).map_or(&[], |v| v.as_slice())
    }
}
