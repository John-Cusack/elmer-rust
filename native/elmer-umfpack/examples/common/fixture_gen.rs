//! Deterministic test matrices for the bitwise comparison against the
//! pristine library. Shared by `tests/` (via `#[path]`) and the
//! `gen_fixtures` example that writes them for the C harness.

#![allow(dead_code, clippy::needless_range_loop)]

use std::collections::BTreeSet;

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03)
    }
    pub fn next_u64(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    pub fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below(hi - lo + 1)
    }
    pub fn sym(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
}

/// Compressed-column matrix (rows sorted, no duplicates) and a right-hand side.
pub struct Fixture {
    pub name: String,
    pub n: i32,
    pub ap: Vec<i32>,
    pub ai: Vec<i32>,
    pub ax: Vec<f64>,
    pub b: Vec<f64>,
}

fn build(
    name: String,
    cols: Vec<BTreeSet<usize>>,
    rng: &mut Rng,
    diag_scale: f64,
    diag_shift: f64,
) -> Fixture {
    let n = cols.len();
    let mut ap = vec![0i32];
    let mut ai = Vec::new();
    let mut ax = Vec::new();
    for (j, c) in cols.iter().enumerate() {
        for &i in c {
            ai.push(i as i32);
            let v = rng.sym();
            ax.push(if i == j {
                v * diag_scale + diag_shift
            } else {
                v
            });
        }
        ap.push(ai.len() as i32);
    }
    let b = (0..n).map(|_| rng.sym()).collect();
    Fixture {
        name,
        n: n as i32,
        ap,
        ai,
        ax,
        b,
    }
}

fn grid(nx: usize, ny: usize, nine: bool, dof: usize) -> Vec<BTreeSet<usize>> {
    let n = nx * ny * dof;
    let mut cols = vec![BTreeSet::new(); n];
    for x in 0..nx {
        for y in 0..ny {
            for dx in -1i64..=1 {
                for dy in -1i64..=1 {
                    if !nine && dx.abs() + dy.abs() > 1 {
                        continue;
                    }
                    let (xx, yy) = (x as i64 + dx, y as i64 + dy);
                    if xx < 0 || yy < 0 || xx >= nx as i64 || yy >= ny as i64 {
                        continue;
                    }
                    for a in 0..dof {
                        for c in 0..dof {
                            cols[(x * ny + y) * dof + a]
                                .insert((xx as usize * ny + yy as usize) * dof + c);
                        }
                    }
                }
            }
        }
    }
    cols
}

/// The fixture set: `count` matrices cycling through eight structural kinds.
pub fn fixtures(count: u64) -> Vec<Fixture> {
    let mut out = Vec::new();
    for seed in 0..count {
        let mut rng = Rng::new(seed);
        let kind = seed % 8;
        let name = format!("k{kind}_{seed:02}");
        let f = match kind {
            0 => {
                // random unsymmetric with a full diagonal
                let n = rng.range(5, 300);
                let d = 1.0 + 7.0 * rng.unit();
                let mut cols = vec![BTreeSet::new(); n];
                for (j, c) in cols.iter_mut().enumerate() {
                    c.insert(j);
                    let k = (-(1.0 - rng.unit()).ln() * d) as usize;
                    for _ in 0..k {
                        c.insert(rng.below(n));
                    }
                }
                build(name, cols, &mut rng, 1.0, 4.0)
            }
            1 => {
                // random symmetric pattern
                let n = rng.range(5, 400);
                let mut cols: Vec<BTreeSet<usize>> = (0..n).map(|j| BTreeSet::from([j])).collect();
                for _ in 0..3 * n {
                    let (i, j) = (rng.below(n), rng.below(n));
                    cols[j].insert(i);
                    cols[i].insert(j);
                }
                build(name, cols, &mut rng, 1.0, 4.0)
            }
            2 => {
                let (nx, ny) = (rng.range(3, 25), rng.range(3, 25));
                build(name, grid(nx, ny, false, 1), &mut rng, 1.0, 4.0)
            }
            3 => {
                let (nx, ny) = (rng.range(3, 15), rng.range(3, 15));
                build(name, grid(nx, ny, true, 2), &mut rng, 1.0, 4.0)
            }
            4 => {
                // dense rows and columns
                let n = rng.range(50, 600);
                let mut cols: Vec<BTreeSet<usize>> = (0..n).map(|j| BTreeSet::from([j])).collect();
                for _ in 0..2 * n {
                    let (i, j) = (rng.below(n), rng.below(n));
                    cols[j].insert(i);
                }
                for _ in 0..rng.range(1, 3) {
                    let r = rng.below(n);
                    for c in cols.iter_mut() {
                        if rng.unit() < 0.9 {
                            c.insert(r);
                        }
                    }
                    let c = rng.below(n);
                    for i in 0..n {
                        if rng.unit() < 0.9 {
                            cols[c].insert(i);
                        }
                    }
                }
                build(name, cols, &mut rng, 1.0, 4.0)
            }
            5 => {
                // lower triangular part (singletons) plus a random block
                let n = rng.range(10, 300);
                let mut cols: Vec<BTreeSet<usize>> = (0..n).map(|j| BTreeSet::from([j])).collect();
                for j in 0..n {
                    for _ in 0..rng.range(0, 3) {
                        let i = rng.range(j, n - 1);
                        cols[j].insert(i);
                    }
                }
                let m = n / 3;
                for _ in 0..4 * m {
                    let (i, j) = (rng.below(m), rng.below(m));
                    cols[j].insert(i);
                }
                build(name, cols, &mut rng, 1.0, 4.0)
            }
            6 => {
                // weak diagonal: partial pivoting has to leave the diagonal
                let (nx, ny) = (rng.range(3, 20), rng.range(3, 20));
                build(name, grid(nx, ny, true, 1), &mut rng, 0.01, 0.0)
            }
            _ => {
                // 1D chain, some diagonal entries missing
                let n = rng.range(2, 200);
                let mut cols = vec![BTreeSet::new(); n];
                for j in 0..n {
                    if rng.unit() < 0.8 {
                        cols[j].insert(j);
                    }
                    if j > 0 {
                        cols[j].insert(j - 1);
                        cols[j - 1].insert(j);
                    }
                }
                build(name, cols, &mut rng, 1.0, 2.0)
            }
        };
        out.push(f);
    }
    out
}

/// The harness text format: `n nz`, `Ap`, `Ai`, `Ax` (hex bits), `b` + rhs.
pub fn to_text(f: &Fixture) -> String {
    let mut s = format!("{} {}\n", f.n, f.ai.len());
    let join_i = |v: &[i32]| {
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    };
    let join_d = |v: &[f64]| {
        v.iter()
            .map(|x| format!("{:016x}", x.to_bits()))
            .collect::<Vec<_>>()
            .join(" ")
    };
    s += &join_i(&f.ap);
    s += "\n";
    s += &join_i(&f.ai);
    s += "\n";
    s += &join_d(&f.ax);
    s += "\nb ";
    s += &join_d(&f.b);
    s += "\n";
    s
}
