//! `Apply Conforming BCs`: the nodes of a `Conforming BC` boundary share the
//! dof of their partner on the target boundary (`Mesh % PeriodicPerm`),
//! with the sign flipped for an antisymmetric pair (`Mesh % PeriodicFlip`).
//! `elmer-mortar` builds the map (`GeneratePeriodicProjectors`); this module
//! applies it to a solver's permutation (`InitialPermutation`, Lists.F90)
//! and to its local systems (`FlipPeriodicLocalMatrix` and
//! `FlipPeriodicLocalForce`, DefUtils.F90).

use elmer_core::{LocalSystem, NO_DOF};

use crate::backend::Perm;

/// `Mesh % PeriodicPerm` and `Mesh % PeriodicFlip` of the nodes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConformingMap {
    /// Node -> the node whose dof it takes (0-based), [`NO_DOF`] where none.
    pub partner: Vec<u32>,
    /// Node -> its value is the negative of its partner's.
    pub flip: Vec<bool>,
}

impl ConformingMap {
    /// Whether node `n` takes another node's dof.
    pub fn is_mapped(&self, n: u32) -> bool {
        self.partner.get(n as usize).is_some_and(|&p| p != NO_DOF)
    }

    fn flipped(&self, n: u32) -> bool {
        self.flip.get(n as usize).copied().unwrap_or(false)
    }
}

/// `InitialPermutation`'s conforming block: the mapped nodes are dropped,
/// the rest renumbered in node order, and each mapped node then takes its
/// partner's slot. Returns `Solver % PeriodicFlipActive`: whether a mapped
/// node of the field is flipped.
pub fn apply_conforming(perm: &mut Perm, map: &ConformingMap) -> bool {
    let m = perm.perm.len();
    let mapped = |i: usize| map.partner.get(i).is_some_and(|&p| p != NO_DOF);
    let mut k = 0u32;
    let mut renum = vec![NO_DOF; m];
    for i in 0..m {
        if perm.perm[i] != NO_DOF && !mapped(i) {
            renum[i] = k;
            k += 1;
        }
    }
    let mut flips = 0;
    for i in 0..m {
        if perm.perm[i] != NO_DOF && mapped(i) {
            let j = map.partner[i] as usize;
            renum[i] = renum.get(j).copied().unwrap_or(NO_DOF);
            if map.flipped(i as u32) {
                flips += 1;
            }
        }
    }
    perm.perm = renum;
    perm.n_slots = k as usize;
    flips > 0
}

/// `FlipPeriodicLocalMatrix` and `FlipPeriodicLocalForce` for the element
/// with `nodes` (one per local function, `local.n / nodes.len()` dofs
/// each): an entry coupling a flipped and an unflipped node changes sign,
/// as does the force of a flipped node. `false` when no node is flipped and
/// `local` is unchanged.
pub fn flip_local(local: &mut LocalSystem, nodes: &[u32], map: &ConformingMap) -> bool {
    if !nodes.iter().any(|&n| map.flipped(n)) {
        return false;
    }
    let n = nodes.len();
    let m = local.n;
    let dofs = m / n;
    for i in 0..n {
        for j in 0..n {
            if map.flipped(nodes[i]) != map.flipped(nodes[j]) {
                for k in 0..dofs {
                    for l in 0..dofs {
                        let x = &mut local.stiff[(dofs * i + k) * m + dofs * j + l];
                        *x = -*x;
                    }
                }
            }
        }
    }
    for i in 0..n {
        if map.flipped(nodes[i]) {
            for k in 0..dofs {
                local.force[dofs * i + k] = -local.force[dofs * i + k];
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mapped_nodes_take_their_partners_slot_after_renumbering() {
        // Nodes 0..5 all active; node 4 -> node 1 (flipped), node 2 -> node 0.
        let mut p = Perm {
            perm: vec![0, 1, 2, 3, 4],
            n_slots: 5,
        };
        let map = ConformingMap {
            partner: vec![NO_DOF, NO_DOF, 0, NO_DOF, 1],
            flip: vec![false, false, false, false, true],
        };
        assert!(apply_conforming(&mut p, &map));
        assert_eq!(p.perm, vec![0, 1, 0, 2, 1]);
        assert_eq!(p.n_slots, 3);
    }

    #[test]
    fn flips_change_the_sign_of_mixed_couplings() {
        let map = ConformingMap {
            partner: vec![NO_DOF, 0],
            flip: vec![false, true],
        };
        let mut l = LocalSystem {
            n: 2,
            stiff: vec![1.0, 2.0, 3.0, 4.0],
            mass: None,
            damp: None,
            force: vec![5.0, 6.0],
        };
        assert!(flip_local(&mut l, &[0, 1], &map));
        assert_eq!(l.stiff, vec![1.0, -2.0, -3.0, 4.0]);
        assert_eq!(l.force, vec![5.0, -6.0]);
    }
}
