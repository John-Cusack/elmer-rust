//! Tiny Elmer-format meshes for the fixtures: a unit square cut into
//! `n x n` cells of two triangles each. Boundary ids: 1 bottom (y=0),
//! 2 right (x=1), 3 top (y=1), 4 left (x=0). Body id 1.

pub struct SquareMesh {
    pub n: usize,
}

impl SquareMesh {
    pub fn new(n: usize) -> Self {
        assert!(n >= 1);
        SquareMesh { n }
    }

    fn node(&self, ix: usize, iy: usize) -> usize {
        iy * (self.n + 1) + ix + 1
    }

    /// `(mesh.header, mesh.nodes, mesh.elements, mesh.boundary)`.
    pub fn files(&self) -> [(&'static str, String); 4] {
        let n = self.n;
        let h = 1.0 / n as f64;
        let mut nodes = String::new();
        for iy in 0..=n {
            for ix in 0..=n {
                nodes.push_str(&format!(
                    "{} -1 {:.12} {:.12} 0.0\n",
                    self.node(ix, iy),
                    ix as f64 * h,
                    iy as f64 * h
                ));
            }
        }
        let mut elems = String::new();
        let mut bnd = String::new();
        let (mut eid, mut bid) = (0, 0);
        let mut edge = |bnd: &mut String, tag: u32, parent: usize, a: usize, b: usize| {
            bid += 1;
            bnd.push_str(&format!("{bid} {tag} {parent} 0 202 {a} {b}\n"));
        };
        for iy in 0..n {
            for ix in 0..n {
                let (n00, n10) = (self.node(ix, iy), self.node(ix + 1, iy));
                let (n11, n01) = (self.node(ix + 1, iy + 1), self.node(ix, iy + 1));
                elems.push_str(&format!("{} 1 303 {n00} {n10} {n11}\n", eid + 1));
                elems.push_str(&format!("{} 1 303 {n00} {n11} {n01}\n", eid + 2));
                let (lower, upper) = (eid + 1, eid + 2);
                eid += 2;
                if iy == 0 {
                    edge(&mut bnd, 1, lower, n00, n10);
                }
                if ix == n - 1 {
                    edge(&mut bnd, 2, lower, n10, n11);
                }
                if iy == n - 1 {
                    edge(&mut bnd, 3, upper, n11, n01);
                }
                if ix == 0 {
                    edge(&mut bnd, 4, upper, n01, n00);
                }
            }
        }
        let header = format!(
            "{} {} {}\n2\n202 {}\n303 {}\n",
            (n + 1) * (n + 1),
            eid,
            bid,
            bid,
            eid
        );
        [
            ("mesh.header", header),
            ("mesh.nodes", nodes),
            ("mesh.elements", elems),
            ("mesh.boundary", bnd),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_and_boundary_tags() {
        let f = SquareMesh::new(2).files();
        assert_eq!(f[0].1, "9 8 8\n2\n202 8\n303 8\n");
        assert_eq!(f[1].1.lines().count(), 9);
        assert_eq!(f[2].1.lines().count(), 8);
        let tags: Vec<u32> = f[3]
            .1
            .lines()
            .map(|l| l.split(' ').nth(1).unwrap().parse().unwrap())
            .collect();
        for t in 1..=4 {
            assert_eq!(tags.iter().filter(|x| **x == t).count(), 2, "tag {t}");
        }
    }

    #[test]
    fn boundary_parents_contain_their_edge() {
        let f = SquareMesh::new(3).files();
        let elems: Vec<Vec<usize>> = f[2]
            .1
            .lines()
            .map(|l| l.split(' ').skip(3).map(|x| x.parse().unwrap()).collect())
            .collect();
        for l in f[3].1.lines() {
            let t: Vec<usize> = l.split(' ').map(|x| x.parse().unwrap()).collect();
            let parent = &elems[t[2] - 1];
            assert!(parent.contains(&t[5]) && parent.contains(&t[6]), "{l}");
        }
    }
}
