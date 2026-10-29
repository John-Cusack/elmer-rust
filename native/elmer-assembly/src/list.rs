//! Row-wise sorted connectivity, the `ListMatrix_t` of ListMatrix.F90.
//!
//! Elmer keeps each row as a linked list sorted by column with no
//! duplicates (`List_AddMatrixIndexes`, ListMatrix.F90:451); `Degree` is the
//! row length. A sorted `Vec` per row has the same content and iteration
//! order, which is all the bandwidth optimiser and the CSR builder look at.

/// Sorted, duplicate-free column sets, one per row. Indices are 0-based.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ListMatrix {
    rows: Vec<Vec<u32>>,
}

impl ListMatrix {
    pub fn new(nrows: usize) -> Self {
        Self {
            rows: vec![Vec::new(); nrows],
        }
    }

    pub fn nrows(&self) -> usize {
        self.rows.len()
    }

    /// Columns of row `i`, ascending.
    pub fn row(&self, i: usize) -> &[u32] {
        &self.rows[i]
    }

    /// Elmer `ListMatrix(i) % Degree`.
    pub fn degree(&self, i: usize) -> usize {
        self.rows[i].len()
    }

    /// Total stored entries (Elmer `Model % TotalMatrixElements`).
    pub fn nnz(&self) -> usize {
        self.rows.iter().map(Vec::len).sum()
    }

    /// `List_GetMatrixIndex(List, k1, k2)`: make sure entry `(k1, k2)` exists.
    pub fn add(&mut self, k1: usize, k2: u32) {
        let row = &mut self.rows[k1];
        if let Err(pos) = row.binary_search(&k2) {
            row.insert(pos, k2);
        }
    }

    /// `List_AddMatrixIndexes(List, k1, n, Ind)`: merge an ascending index
    /// array (duplicates allowed) into row `k1`.
    pub fn add_sorted(&mut self, k1: usize, ind: &[u32]) {
        debug_assert!(ind.windows(2).all(|w| w[0] <= w[1]));
        let row = &mut self.rows[k1];
        if row.is_empty() {
            row.extend_from_slice(ind);
            row.dedup();
            return;
        }
        let mut merged = Vec::with_capacity(row.len() + ind.len());
        let (mut a, mut b) = (0, 0);
        while a < row.len() || b < ind.len() {
            let next = if b == ind.len() || (a < row.len() && row[a] <= ind[b]) {
                a += 1;
                row[a - 1]
            } else {
                b += 1;
                ind[b - 1]
            };
            if merged.last() != Some(&next) {
                merged.push(next);
            }
        }
        *row = merged;
    }

    /// `PackSortAdd` (ElementUtils.F90:1038) for element indices already
    /// mapped through the initial permutation: drop absent ones, sort, and
    /// add the whole set to each of its rows.
    pub fn add_clique(&mut self, ind: &mut [u32]) {
        ind.sort_unstable();
        for i in 0..ind.len() {
            let k1 = ind[i] as usize;
            self.add_sorted(k1, ind);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_sorted_merges_without_duplicates() {
        let mut l = ListMatrix::new(1);
        l.add_sorted(0, &[3, 3, 5]);
        l.add_sorted(0, &[1, 3, 4, 9, 9]);
        l.add(0, 2);
        l.add(0, 4);
        assert_eq!(l.row(0), &[1, 2, 3, 4, 5, 9]);
        assert_eq!(l.degree(0), 6);
    }

    #[test]
    fn clique_is_symmetric() {
        let mut l = ListMatrix::new(4);
        l.add_clique(&mut [3, 0, 2]);
        assert_eq!(l.row(0), &[0, 2, 3]);
        assert_eq!(l.row(1), &[] as &[u32]);
        assert_eq!(l.row(3), &[0, 2, 3]);
        assert_eq!(l.nnz(), 9);
    }
}
