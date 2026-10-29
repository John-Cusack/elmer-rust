//! Element edges for the Whitney basis.

use std::collections::HashMap;

/// `GetEdgeMap(family)` (ElementBasis.F90): the local edges of an element
/// family as 0-based node pairs, in the order the Whitney basis functions
/// are numbered.
pub fn edge_map(family: u16) -> &'static [[usize; 2]] {
    match family {
        3 => &[[0, 1], [1, 2], [2, 0]],
        4 => &[[0, 1], [1, 2], [2, 3], [3, 0]],
        5 => &[[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]],
        6 => &[[0, 1], [1, 2], [3, 2], [0, 3], [0, 4], [1, 4], [2, 4], [3, 4]],
        7 => &[[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3], [0, 3], [1, 4], [2, 5]],
        8 => &[
            [0, 1],
            [1, 2],
            [3, 2],
            [0, 3],
            [4, 5],
            [5, 6],
            [7, 6],
            [4, 7],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
        ],
        _ => &[],
    }
}

/// Global edges keyed by their (sorted) node pair, numbered in first-use
/// order.
#[derive(Default)]
pub struct Edges {
    index: HashMap<(u32, u32), u32>,
    pub n: usize,
}

impl Edges {
    fn key(a: u32, b: u32) -> (u32, u32) {
        if a < b {
            (a, b)
        } else {
            (b, a)
        }
    }

    /// Numbers the edges of an element (adding new ones) and returns their
    /// global ids in local edge order.
    pub fn add_element(&mut self, code: u16, nodes: &[u32]) -> Vec<u32> {
        edge_map(code / 100)
            .iter()
            .map(|&[a, b]| {
                let next = self.n as u32;
                let id = *self.index.entry(Self::key(nodes[a], nodes[b])).or_insert(next);
                if id == next {
                    self.n += 1;
                }
                id
            })
            .collect()
    }

    /// Global ids of the edges of an element whose edges all exist.
    pub fn of_element(&self, code: u16, nodes: &[u32]) -> Option<Vec<u32>> {
        edge_map(code / 100)
            .iter()
            .map(|&[a, b]| self.index.get(&Self::key(nodes[a], nodes[b])).copied())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_edges_are_numbered_once() {
        let mut e = Edges::default();
        let a = e.add_element(504, &[0, 1, 2, 3]);
        let b = e.add_element(504, &[1, 2, 3, 4]);
        assert_eq!(a, vec![0, 1, 2, 3, 4, 5]);
        // (1,2), (2,3), (3,1) are shared with the first tetra.
        assert_eq!(&b[..3], &[1, 5, 4]);
        assert_eq!(e.n, 9);
        assert_eq!(e.of_element(303, &[3, 2, 1]), Some(vec![5, 1, 4]));
        assert_eq!(e.of_element(303, &[0, 1, 4]), None);
    }
}
