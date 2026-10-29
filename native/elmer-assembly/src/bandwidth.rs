//! Bandwidth computation and Elmer's reverse Cuthill-McKee
//! (BandwidthOptimize.F90).
//!
//! This is a literal port, not a textbook RCM: the root search uses Elmer's
//! depth-first `Levelize` (its "levels" are DFS depths) and its peculiar
//! restart rule, and the renumbering visits neighbours in ascending
//! initial-numbering order. Any of these differences changes the ordering,
//! and therefore the pattern of every saved linear system.

use crate::list::ListMatrix;

/// Half bandwidths as Elmer prints them ("Initial bandwidth for ...",
/// "Optimized bandwidth for ..."). Both include Elmer's `+1`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BandwidthReport {
    pub initial: usize,
    /// `None` when optimisation was not requested.
    pub optimized: Option<usize>,
    /// Whether the optimised ordering was kept. False when it was not
    /// requested, or rejected because it was worse and
    /// `Optimize Bandwidth Use Always` was off.
    pub accepted: bool,
}

impl BandwidthReport {
    /// The half bandwidth `OptimizeBandwidth` returns (and
    /// `Matrix % Subband` is built from).
    pub fn half_bandwidth(&self) -> usize {
        match self.optimized {
            Some(b) if self.accepted => b,
            _ => self.initial,
        }
    }
}

/// `ComputeBandwidth(N, List)` (BandwidthOptimize.F90:71) without a
/// reordering: max |i - j| over the stored entries. Not `+1`'d.
pub fn bandwidth(list: &ListMatrix) -> usize {
    let mut hb = 0usize;
    for i in 0..list.nrows() {
        for &k in list.row(i) {
            hb = hb.max(i.abs_diff(k as usize));
        }
    }
    hb
}

/// `ComputeBandwidth(N, List, Reorder, InvInitialReorder)`: the bandwidth
/// of `list` (rows in initial numbering) after renumbering row `i` to
/// `new_index[i]`. Not `+1`'d.
pub fn bandwidth_reordered(list: &ListMatrix, new_index: &[u32]) -> usize {
    let mut hb = 0usize;
    for i in 0..list.nrows() {
        let ni = new_index[i] as usize;
        for &k in list.row(i) {
            hb = hb.max(ni.abs_diff(new_index[k as usize] as usize));
        }
    }
    hb
}

/// Elmer's RCM on the rows `0..list.nrows()` of a list matrix in the
/// initial numbering. Returns `new_index[i]`, the 0-based optimised index
/// of initial index `i` (`DoneIndex` after the reversal).
pub fn rcm_order(list: &ListMatrix) -> Vec<u32> {
    let n = list.nrows();
    if n == 0 {
        return Vec::new();
    }
    let mut level = vec![0i64; n];
    let mut done = vec![false; n];
    let mut max_level = 0i64;

    // Search for the node to start from: first of minimum degree.
    let mut start = 0usize;
    let mut min_degree = list.degree(0);
    for i in 0..n {
        if list.degree(i) < min_degree {
            start = i;
            min_degree = list.degree(i);
        }
    }

    levelize(list, start, &mut level, &mut done, &mut max_level);

    let mut new_root = true;
    while new_root {
        new_root = false;
        let mut min_degree = list.degree(start);
        let mut k = start;
        for (i, &lev) in level.iter().enumerate() {
            if lev == max_level && list.degree(i) < min_degree {
                k = i;
                min_degree = list.degree(i);
            }
        }
        if k != start {
            let j = max_level;
            max_level = 0;
            done.iter_mut().for_each(|d| *d = false);
            levelize(list, k, &mut level, &mut done, &mut max_level);
            // Elmer keeps the new root only if its depth came out smaller.
            if j > max_level {
                new_root = true;
                start = k;
            }
        }
    }

    // Cuthill-McKee: breadth-first renumbering from `start`, neighbours in
    // list order, next unvisited index when a component runs out.
    const UNSET: u32 = u32::MAX;
    let mut order: Vec<u32> = Vec::with_capacity(n); // PermLocal
    let mut index = vec![UNSET; n]; // DoneIndex
    order.push(start as u32);
    index[start] = 0;
    let mut first_free = 0usize;
    for i in 0..n {
        if i == order.len() {
            while index[first_free] != UNSET {
                first_free += 1;
            }
            index[first_free] = order.len() as u32;
            order.push(first_free as u32);
        }
        let row = order[i] as usize;
        for &k in list.row(row) {
            let k = k as usize;
            if k < n && index[k] == UNSET {
                index[k] = order.len() as u32;
                order.push(k as u32);
            }
        }
    }

    // Reverse for profile reduction.
    let mut new_index = vec![0u32; n];
    for (i, &node) in order.iter().enumerate() {
        new_index[node as usize] = (n - 1 - i) as u32;
    }
    new_index
}

/// `Levelize` (BandwidthOptimize.F90:374): Elmer's iterative depth-first
/// walk. `level[n]` is the DFS depth at which `n` was first reached, and
/// `max_level` the deepest depth seen.
fn levelize(
    list: &ListMatrix,
    root: usize,
    level: &mut [i64],
    done: &mut [bool],
    max_level: &mut i64,
) {
    let local = list.nrows();
    let mut n = root;
    let mut lev = 0i64;
    // A list pointer is (row, position); it is "associated" while the
    // position is inside the row.
    let valid = |p: (usize, usize)| p.1 < list.degree(p.0);
    let mut stack: Vec<(usize, usize)> = Vec::new();
    let mut p = (n, 0usize);
    while valid(p) {
        stack.push(p);
        level[n] = lev;
        done[n] = true;
        *max_level = (*max_level).max(lev);

        p = (n, 0);
        loop {
            if valid(p) {
                n = list.row(p.0)[p.1] as usize;
                if n < local && !done[n] {
                    lev += 1;
                    break;
                }
            } else if let Some(top) = stack.pop() {
                p = top;
                lev -= 1;
            } else {
                break;
            }
            p.1 += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(order: &[u32]) -> ListMatrix {
        // A path graph visiting `order`, with self loops.
        let n = order.len();
        let mut l = ListMatrix::new(n);
        for i in 0..n {
            l.add(i, i as u32);
        }
        for w in order.windows(2) {
            l.add(w[0] as usize, w[1]);
            l.add(w[1] as usize, w[0]);
        }
        l
    }

    #[test]
    fn rcm_untangles_a_scrambled_path() {
        let l = path(&[0, 4, 2, 5, 1, 3]);
        assert_eq!(bandwidth(&l), 4);
        let ni = rcm_order(&l);
        assert_eq!(bandwidth_reordered(&l, &ni), 1);
        let mut seen = ni.clone();
        seen.sort_unstable();
        assert_eq!(seen, (0..6).collect::<Vec<u32>>());
    }

    #[test]
    fn rcm_numbers_disconnected_components() {
        let mut l = ListMatrix::new(4);
        for i in 0..4 {
            l.add(i, i as u32);
        }
        l.add(0, 2);
        l.add(2, 0);
        let ni = rcm_order(&l);
        let mut seen = ni.clone();
        seen.sort_unstable();
        assert_eq!(seen, vec![0, 1, 2, 3]);
    }

    #[test]
    fn triangle_reverses() {
        // MinimalExample: one triangle, all rows full. Elmer's dump shows the
        // order 3 2 1 (pristine run, perm file).
        let mut l = ListMatrix::new(3);
        l.add_clique(&mut [0, 1, 2]);
        assert_eq!(rcm_order(&l), vec![2, 1, 0]);
    }
}
