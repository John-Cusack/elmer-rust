//! The `.sif` sections the solver framework reads, and the body/boundary
//! lookups built from them (Elmer `Model_t`).

use std::collections::BTreeMap;

use elmer_core::{ParamMap, Params, Value};

/// The parsed model, one [`ParamMap`] per section. Section vectors are
/// indexed by Elmer id minus one (`solvers[0]` is `Solver 1`); `elmer-sif`
/// fills missing ids with empty sections.
#[derive(Clone, Debug, Default)]
pub struct Model {
    pub header: ParamMap,
    pub simulation: ParamMap,
    pub constants: ParamMap,
    pub solvers: Vec<ParamMap>,
    pub equations: Vec<ParamMap>,
    pub bodies: Vec<ParamMap>,
    pub materials: Vec<ParamMap>,
    pub body_forces: Vec<ParamMap>,
    pub initial_conditions: Vec<ParamMap>,
    pub boundary_conditions: Vec<ParamMap>,
}

/// The sections that apply to one bulk element, through its body.
#[derive(Copy, Clone, Debug)]
pub struct BodyParams<'a> {
    /// 1-based body id.
    pub body_id: usize,
    pub body: &'a ParamMap,
    pub equation: Option<&'a ParamMap>,
    pub material: Option<&'a ParamMap>,
    pub body_force: Option<&'a ParamMap>,
    pub initial_condition: Option<&'a ParamMap>,
}

/// Integer-list keyword (`Target Boundaries(2) = 1 3`).
pub(crate) fn int_list<'a>(p: &'a ParamMap, key: &str) -> &'a [i32] {
    match p.get(key) {
        Some(Value::Integer(v)) => v,
        _ => &[],
    }
}

impl Model {
    /// 1-based body id for a bulk element tag: a body whose `Target Bodies`
    /// lists the tag, else the body with that id.
    pub fn body_of_tag(&self, tag: i32) -> Option<usize> {
        for (i, b) in self.bodies.iter().enumerate() {
            if int_list(b, "target bodies").contains(&tag) {
                return Some(i + 1);
            }
        }
        (tag >= 1 && (tag as usize) <= self.bodies.len()).then_some(tag as usize)
    }

    /// 1-based boundary condition id for a boundary tag (`Target Boundaries`).
    pub fn bc_of_tag(&self, tag: i32) -> Option<usize> {
        self.boundary_conditions
            .iter()
            .position(|bc| int_list(bc, "target boundaries").contains(&tag))
            .map(|i| i + 1)
    }

    /// The section a body's integer keyword (`Material`, `Equation`, ...) points at.
    fn linked<'a>(&self, body: &ParamMap, key: &str, sections: &'a [ParamMap]) -> Option<&'a ParamMap> {
        let k = body.integer(key)?;
        (k >= 1).then(|| sections.get(k as usize - 1)).flatten()
    }

    pub fn body_params(&self, body_id: usize) -> Option<BodyParams<'_>> {
        let body = self.bodies.get(body_id.checked_sub(1)?)?;
        Some(BodyParams {
            body_id,
            body,
            equation: self.linked(body, "equation", &self.equations),
            material: self.linked(body, "material", &self.materials),
            body_force: self.linked(body, "body force", &self.body_forces),
            initial_condition: self.linked(body, "initial condition", &self.initial_conditions),
        })
    }

    /// Whether solver `solver_id` (1-based) is active on body `body_id`:
    /// the body's equation lists it in `Active Solvers`, or sets the
    /// logical keyword named by the solver's `Equation` (CheckElementEquation,
    /// MeshBasics.F90).
    pub fn solver_active_on_body(&self, solver_id: usize, body_id: usize) -> bool {
        let Some(eq) = self.body_params(body_id).and_then(|b| b.equation) else {
            return false;
        };
        if int_list(eq, "active solvers").contains(&(solver_id as i32)) {
            return true;
        }
        let name = self
            .solvers
            .get(solver_id - 1)
            .and_then(|s| s.string("equation"));
        name.is_some_and(|n| eq.logical(n) == Some(true))
    }

    /// Body id -> whether solver `solver_id` is active there.
    pub fn active_bodies(&self, solver_id: usize) -> BTreeMap<usize, bool> {
        (1..=self.bodies.len())
            .map(|b| (b, self.solver_active_on_body(solver_id, b)))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sec(entries: &[(&str, Value)]) -> ParamMap {
        let mut p = ParamMap::new();
        for (k, v) in entries {
            p.insert(k, v.clone());
        }
        p
    }

    #[test]
    fn body_and_bc_lookup() {
        let m = Model {
            solvers: vec![
                sec(&[("Equation", Value::Str("HeatB".into()))]),
                sec(&[("Equation", Value::Str("HeatA".into()))]),
            ],
            equations: vec![
                sec(&[("Active Solvers", Value::Integer(vec![2]))]),
                sec(&[("heatb", Value::Logical(true))]),
            ],
            bodies: vec![
                sec(&[("Equation", Value::Integer(vec![1])), ("Material", Value::Integer(vec![1]))]),
                sec(&[("Equation", Value::Integer(vec![2])), ("Target Bodies", Value::Integer(vec![7]))]),
            ],
            materials: vec![sec(&[("Heat Conductivity", Value::real(1.0))])],
            boundary_conditions: vec![sec(&[("Target Boundaries", Value::Integer(vec![1, 3]))])],
            ..Model::default()
        };
        assert_eq!(m.body_of_tag(1), Some(1));
        assert_eq!(m.body_of_tag(7), Some(2));
        assert_eq!(m.body_of_tag(9), None);
        assert_eq!(m.bc_of_tag(3), Some(1));
        assert_eq!(m.bc_of_tag(2), None);
        assert!(m.solver_active_on_body(2, 1));
        assert!(!m.solver_active_on_body(1, 1));
        assert!(m.solver_active_on_body(1, 2));
        let bp = m.body_params(1).unwrap();
        assert_eq!(bp.material.unwrap().real("heat conductivity"), Some(1.0));
        assert!(bp.body_force.is_none());
    }
}
