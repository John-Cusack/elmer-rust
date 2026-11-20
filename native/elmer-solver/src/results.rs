//! The `res:` result scalars solvers leave in the Simulation section, in the
//! order Elmer's lists keep them: `ListAddConstReal` appends a new keyword at
//! the end, and SaveScalars lists them from the head. The keyword map is
//! sorted, so the order of the keys added through [`add_result`] is kept in
//! the reserved entry [`RESULT_ORDER`].

use elmer_core::{normalize_key, ParamMap, Params, Value};

/// The Simulation entry holding the `res:` keys in the order they were
/// first added, separated by newlines.
pub const RESULT_ORDER: &str = "__result scalar order";

/// `ListAddConstReal(Model % Simulation, key, value)`: set `key`, recording
/// it as the last one when new.
pub fn add_result(sim: &mut ParamMap, key: &str, value: f64) {
    let k = normalize_key(key);
    if sim.get(&k).is_none() {
        let mut order = match sim.get(RESULT_ORDER) {
            Some(Value::Str(s)) => s.clone(),
            _ => String::new(),
        };
        if !order.is_empty() {
            order.push('\n');
        }
        order.push_str(&k);
        sim.insert(RESULT_ORDER, Value::Str(order));
    }
    sim.insert(key, Value::real(value));
}

/// The position of each recorded key ([`add_result`]), to order the
/// Simulation's `res:` keys as Elmer lists them; the keys added otherwise
/// follow, in the map's order.
pub fn result_rank(sim: &ParamMap, key: &str) -> usize {
    let k = normalize_key(key);
    match sim.get(RESULT_ORDER) {
        Some(Value::Str(s)) => s.split('\n').position(|x| x == k).unwrap_or(usize::MAX),
        _ => usize::MAX,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_first_insertion_order() {
        let mut p = ParamMap::new();
        add_result(&mut p, "res: Zeta", 1.0);
        add_result(&mut p, "res: alpha", 2.0);
        add_result(&mut p, "res: Zeta", 3.0);
        assert!(result_rank(&p, "res: zeta") < result_rank(&p, "res: alpha"));
        assert_eq!(result_rank(&p, "res: other"), usize::MAX);
        assert_eq!(p.real("res: zeta"), Some(3.0));
    }
}
