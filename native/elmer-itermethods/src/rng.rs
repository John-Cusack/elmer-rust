//! Deterministic uniform generator for the IDR(s) shadow space.
//!
//! The Fortran uses `RANDOM_NUMBER`, seeded by ElmerSolver from
//! `Simulation: Random Number Seed` (default 314159265). Reproducing
//! libgfortran's stream bit for bit is out of scope, so the shadow vectors
//! differ and iteration counts agree only to within the IDR(s) variation
//! between random shadow spaces. xoshiro256** seeded through splitmix64.

pub(crate) struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub(crate) fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut v = z;
            v = (v ^ (v >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            v = (v ^ (v >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            v ^ (v >> 31)
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }

    fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform in `[0, 1)`.
    pub(crate) fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}
