//! Keyword lookups under Elmer's namespace stack (`ListPushNamespace`,
//! `ListFind`, Lists.F90:2125, :2266).
//!
//! With namespace `ns` active, `ListFind` first tries `"ns name"`. With
//! `Additive Namespaces` (Simulation, ElmerSolver.F90:437) it then walks the
//! stack of earlier namespaces, innermost first; without it, it goes straight
//! to the plain `name`. BlockSolve pushes `outer:` around the outer Krylov
//! iteration and `block:` then `block NN:` around each inner solve, so a key
//! is found as `block 11: k`, then (additive only) `block: k` and `outer: k`,
//! then `k`.
//!
//! [`resolve`] flattens that rule into an owned [`ParamMap`], so the solver
//! crates can read it through [`elmer_core::Params`] unchanged.

use elmer_core::{normalize_key, ParamMap};

/// `section` as seen with the namespaces `stack` pushed (outermost first).
/// Every plain key keeps its value unless a namespaced key that the lookup
/// would find first overrides it.
pub fn resolve(section: &ParamMap, stack: &[&str], additive: bool) -> ParamMap {
    let mut out = section.clone();
    // Lowest priority first, so the innermost namespace is written last.
    let active: &[&str] = match stack.split_last() {
        Some((last, _)) if !additive => std::slice::from_ref(last),
        _ => stack,
    };
    for ns in active {
        let prefix = format!("{} ", normalize_key(ns));
        for (k, v) in &section.0 {
            if let Some(rest) = k.strip_prefix(&prefix) {
                out.0.insert(rest.to_string(), v.clone());
            }
        }
    }
    out
}

/// The namespace BlockSolve pushes for the inner solve of block `i`
/// (0-based): `block 11:` for the first (`'block '//I2S(11*i)//':'`).
pub fn block_namespace(i: usize) -> String {
    format!("block {}:", 11 * (i + 1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{Params, Value};

    fn section() -> ParamMap {
        let mut p = ParamMap::new();
        p.insert("Linear System Max Iterations", Value::Integer(vec![1]));
        p.insert("Outer: Linear System Max Iterations", Value::Integer(vec![2]));
        p.insert("Block: Linear System Max Iterations", Value::Integer(vec![3]));
        p.insert("Block: Linear System Solver", Value::Str("iterative".into()));
        p.insert("block 22: Linear System Max Iterations", Value::Integer(vec![4]));
        p
    }

    #[test]
    fn additive_walks_the_stack() {
        let p = section();
        let v = resolve(&p, &["outer:", "block:", "block 11:"], true);
        assert_eq!(v.integer("linear system max iterations"), Some(3));
        assert_eq!(v.string("linear system solver"), Some("iterative"));
        let v = resolve(&p, &["outer:", "block:", "block 22:"], true);
        assert_eq!(v.integer("linear system max iterations"), Some(4));
        let v = resolve(&p, &["outer:"], true);
        assert_eq!(v.integer("linear system max iterations"), Some(2));
    }

    #[test]
    fn plain_mode_checks_only_the_innermost() {
        let p = section();
        let v = resolve(&p, &["outer:", "block:", "block 11:"], false);
        assert_eq!(v.integer("linear system max iterations"), Some(1));
        assert_eq!(v.string("linear system solver"), None);
        let v = resolve(&p, &["outer:", "block:", "block 22:"], false);
        assert_eq!(v.integer("linear system max iterations"), Some(4));
        assert_eq!(block_namespace(2), "block 33:");
    }
}
