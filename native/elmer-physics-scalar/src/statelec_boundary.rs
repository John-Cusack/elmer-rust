//! `StatElecBoundary.F90`: the electrostatic boundary functions
//! `StatElecBoundaryEnergy`, `...Force`, `...Charge` and `...Spring`
//! (`Real Procedure "StatElecBoundary" "StatElecBoundaryForce"`), which give
//! the energy, force, charge or spring constant of a thin air gap from the
//! keywords of the element they are evaluated on.
//!
//! The functions read `Model % CurrentElement` ([`elmer_expr::current_element`]):
//! a boundary element's Boundary Condition or a bulk element's Material
//! (`GetBC` / `GetMaterial`). From that list they take, at the element's
//! nodes, `Potential Difference`, `Permittivity` (else `Relative
//! Permittivity`, times `Permittivity of Vacuum`), `Gap Height`, the
//! optional `Layer Thickness` / `Layer Permittivity`, and the hole
//! correction of `Hole Type` with `Hole Size`, `Hole Depth` and `Hole
//! Fraction`. As in the Fortran, the values are cached until the current
//! element changes (`SAVE`d).

use std::sync::{Arc, Mutex};

use elmer_core::{ParamMap, Params};
use elmer_expr::{compile_with, CompileOptions, FieldLookup, MatcEngine, Udf, UdfRegistry};

/// The library name the functions are registered under.
pub const LIBRARY: &str = "StatElecBoundary";

/// The section `GetBC(CurrentElement)` / `GetMaterial(CurrentElement)`
/// gives an element: an index into [`Lists::bcs`] or [`Lists::materials`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ListOf {
    Bc(usize),
    Material(usize),
}

/// What the functions read of the model: the sections, the section of every
/// element (`None` without one), and `Permittivity of Vacuum`.
pub struct Lists {
    pub bcs: Vec<ParamMap>,
    pub materials: Vec<ParamMap>,
    pub list_of: Vec<Option<ListOf>>,
    pub eps0: f64,
}

/// `Hole Type` with `Hole Size`, `Hole Depth` and `Hole Fraction` at the
/// element's nodes.
struct Hole {
    kind: String,
    size: Vec<f64>,
    depth: Vec<f64>,
    fraction: Vec<f64>,
}

/// The cached per-element values (`SAVE`d in the Fortran).
#[derive(Default)]
struct Cached {
    element: Option<u32>,
    nodes: Vec<u32>,
    potential_difference: Vec<f64>,
    permittivity: Vec<f64>,
    layer: Option<(Vec<f64>, Vec<f64>)>,
    hole: Option<Hole>,
    aperture: Vec<f64>,
}

struct State {
    lists: Lists,
    fields: Arc<dyn FieldLookup + Send + Sync>,
    matc: Arc<dyn MatcEngine + Send + Sync>,
    udfs: Arc<UdfRegistry>,
    cache: Mutex<Cached>,
}

/// Registers the four functions in `into`. `udfs` resolves procedures in
/// the keywords the functions evaluate, `fields` and `matc` evaluate them.
pub fn register(
    into: &UdfRegistry,
    lists: Lists,
    fields: Arc<dyn FieldLookup + Send + Sync>,
    matc: Arc<dyn MatcEngine + Send + Sync>,
    udfs: Arc<UdfRegistry>,
) {
    let state = Arc::new(State { lists, fields, matc, udfs, cache: Mutex::new(Cached::default()) });
    for (name, mode) in [
        ("StatElecBoundaryEnergy", Mode::Energy),
        ("StatElecBoundaryForce", Mode::Force),
        ("StatElecBoundaryCharge", Mode::Charge),
        ("StatElecBoundarySpring", Mode::Spring),
    ] {
        let s = Arc::clone(&state);
        into.register_fortran(LIBRARY, name, Udf::Node(Arc::new(move |node, _| s.value(node, mode))));
    }
}

#[derive(Clone, Copy)]
enum Mode {
    Energy,
    Force,
    Charge,
    Spring,
}

impl State {
    /// `ListGetReal(List, key, n, NodeIndexes, Found)`: `None` when absent.
    /// A value the evaluator cannot give is NaN, which the caller's result
    /// carries (the user function has no error channel).
    fn at_nodes(&self, list: &ParamMap, key: &str, nodes: &[u32]) -> Option<Vec<f64>> {
        let v = list.get(key)?;
        let mut out = vec![0.0; nodes.len()];
        let ok = compile_with(v, &self.udfs, &CompileOptions::default())
            .ok()
            .filter(|c| c.ncomp() == 1)
            .is_some_and(|c| c.eval_nodes(nodes, &*self.fields, &*self.matc, &mut out).is_ok());
        if !ok {
            out.fill(f64::NAN);
        }
        Some(out)
    }

    /// `StatElecBoundaryGeneric( Model, NodeNumber, dummy, Mode )`.
    fn value(&self, node: u32, mode: Mode) -> f64 {
        let Some(ce) = elmer_expr::current_element() else { return f64::NAN };
        let mut c = self.cache.lock().unwrap();
        if c.element != Some(ce.index) || c.nodes != ce.nodes {
            let Some(list) = self.lists.list_of.get(ce.index as usize).copied().flatten().and_then(|l| match l {
                ListOf::Bc(i) => self.lists.bcs.get(i),
                ListOf::Material(i) => self.lists.materials.get(i),
            }) else {
                return f64::NAN;
            };
            let n = ce.nodes.len();
            let nodes = &ce.nodes;
            let get = |key: &str| self.at_nodes(list, key, nodes);
            let zeros = || vec![0.0; n];
            c.potential_difference = get("potential difference").unwrap_or_else(zeros);
            let given = get("permittivity");
            let got_it = given.is_some();
            let mut permittivity = given.or_else(|| get("relative permittivity")).unwrap_or_else(zeros);
            permittivity.iter_mut().for_each(|p| *p *= self.lists.eps0);
            c.permittivity = permittivity;
            c.layer = get("layer thickness").map(|thickness| {
                let mut lp = get("layer permittivity").unwrap_or_else(zeros);
                // `IF(.NOT.GotIt) LayerPermittivity = 1.0d0`: `GotIt` is the
                // `Permittivity` lookup's.
                if !got_it {
                    lp.fill(1.0);
                }
                lp.iter_mut().for_each(|p| *p *= self.lists.eps0);
                (thickness, lp)
            });
            c.hole = list.string("hole type").map(|t| Hole {
                kind: t.to_string(),
                size: get("hole size").unwrap_or_else(zeros),
                depth: get("hole depth").unwrap_or_else(zeros),
                fraction: get("hole fraction").unwrap_or_else(zeros),
            });
            c.aperture = get("gap height").unwrap_or_else(zeros);
            c.element = Some(ce.index);
            c.nodes = ce.nodes.clone();
        }
        // `IF( NodeIndexes(i) == NodeNumber ) EXIT`: past the end when the
        // node is not the element's, an index the Fortran leaves undefined.
        let Some(i) = c.nodes.iter().position(|&m| m == node) else { return f64::NAN };
        let eps = c.permittivity[i];
        let mut aperture = c.aperture[i];
        if let Some((thickness, lp)) = &c.layer {
            aperture = if c.aperture[i] > thickness[i] {
                (c.aperture[i] - thickness[i]) + thickness[i] * eps / lp[i]
            } else {
                c.aperture[i] * eps / lp[i]
            };
        }
        let (mut alpha, mut beta, mut gamma) = (1.0, 1.0, 1.0);
        let hole = c.hole.as_ref().map(|h| hole_correction(&h.kind, h.size[i], h.depth[i], h.fraction[i], aperture));
        if let Some((a, b, g)) = hole {
            (alpha, beta, gamma) = (a, b, g);
        }
        let dv = c.potential_difference[i];
        let dv2 = dv.powf(2.0);
        match mode {
            Mode::Energy => {
                let v = 0.5 * dv2 * eps / aperture;
                if hole.is_some() { alpha * v } else { v }
            }
            Mode::Force => {
                let v = -0.5 * dv2 * eps / aperture.powf(2.0);
                if hole.is_some() { beta * v } else { v }
            }
            Mode::Charge => {
                let v = dv * eps / aperture.powf(2.0);
                if hole.is_some() { beta * v } else { v }
            }
            Mode::Spring => {
                let v = dv2 * eps / aperture.powf(3.0);
                if hole.is_some() { gamma * v } else { v }
            }
        }
    }
}

/// `ElectrostaticHoleCorrection(holemodel, r, b, p, d, alpha, beta, gamma)`:
/// an unknown hole type gives 1, 1, 1 (with Elmer's warning, not printed
/// here). `b`, the hole depth, is not used by the model.
fn hole_correction(model: &str, r: f64, _b: f64, p: f64, d: f64) -> (f64, f64, f64) {
    let (c1, c2) = match model {
        "slot" => (2.3198, 0.2284),
        "round" => (4.2523, 0.4133),
        "square" => (3.8434, 0.3148),
        _ => return (1.0, 1.0, 1.0),
    };
    let dr = d / r;
    let dom = 1.0 + c1 * dr + c2 * dr.powf(2.0);
    let a = 1.0 - p * 1.0 / dom;
    let da = p * (c1 + 2.0 * c2 * dr) / dom.powf(2.0);
    let dda = p * 2.0 * (c2 - 2.0 * c1.powf(2.0) - 3.0 * c1 * c2 * dr - 3.0 * c2.powf(2.0) * dr.powf(2.0)) / dom.powf(3.0);
    (a, a - da * dr, a - da * dr + 0.5 * dda * dr.powf(2.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::Value;

    struct NoFields;
    impl FieldLookup for NoFields {
        fn value(&self, _: &str, _: u32) -> Option<f64> {
            None
        }
        fn coord(&self, _: u32) -> [f64; 3] {
            [0.0; 3]
        }
        fn time(&self) -> f64 {
            0.0
        }
    }
    struct NoMatc;
    impl MatcEngine for NoMatc {
        fn eval(&self, _: &str, _: &[f64]) -> Result<Vec<f64>, String> {
            Err("no MATC".into())
        }
    }

    /// The force of a 0.1 gap at 1e6 V with eps = 1 * eps0 on boundary
    /// element 1 (BC 1), at either of its nodes:
    /// -0.5 * dV^2 * eps / gap^2. The energy is 0.5 * dV^2 * eps / gap. A
    /// node of another element, or no current element, gives NaN.
    #[test]
    fn force_and_energy_of_a_gap() {
        let mut bc = ParamMap::new();
        bc.insert("Potential Difference", Value::real(1.0e6));
        bc.insert("Relative Permittivity", Value::real(1.0));
        bc.insert("Gap Height", Value::real(0.1));
        let eps0 = 8.8542e-12;
        let lists = Lists { bcs: vec![bc], materials: vec![], list_of: vec![None, Some(ListOf::Bc(0))], eps0 };
        let reg = UdfRegistry::new();
        register(&reg, lists, Arc::new(NoFields), Arc::new(NoMatc), Arc::new(UdfRegistry::new()));
        let Some(Udf::Node(force)) = reg.lookup("StatElecBoundary", "StatElecBoundaryForce") else { panic!() };
        let Some(Udf::Node(energy)) = reg.lookup("StatElecBoundary.so", "statelecboundaryenergy") else { panic!() };
        elmer_expr::set_current_element(1, &[4, 7]);
        let want = -0.5 * 1.0e12 * eps0 / (0.1f64 * 0.1);
        assert_eq!(force(7, &[0.0]), want);
        assert_eq!(force(4, &[0.0]), want);
        assert_eq!(energy(4, &[0.0]), 0.5 * 1.0e12 * eps0 / 0.1);
        assert!(force(5, &[0.0]).is_nan());
        elmer_expr::set_current_element(0, &[1, 2, 3]);
        assert!(force(1, &[0.0]).is_nan());
    }

    /// An unknown hole type leaves the value; `round` scales the force by
    /// beta = a - da * d/r.
    #[test]
    fn hole_correction_factors() {
        assert_eq!(hole_correction("triangle", 1.0, 0.0, 0.5, 0.1), (1.0, 1.0, 1.0));
        let (a, b, _) = hole_correction("round", 2.0, 0.0, 0.5, 1.0);
        let dom = 1.0 + 4.2523 * 0.5 + 0.4133 * 0.25;
        assert_eq!(a, 1.0 - 0.5 / dom);
        assert!((b - (a - 0.5 * (4.2523 + 2.0 * 0.4133 * 0.5) / (dom * dom) * 0.5)).abs() < 1e-15);
    }
}
