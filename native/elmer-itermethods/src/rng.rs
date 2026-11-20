//! libgfortran's `RANDOM_NUMBER` stream for `REAL(8)`, which fills the
//! IDR(s) shadow space.
//!
//! ElmerSolver's `InitializeRandomSeed` (ElmerSolver.F90:874) calls
//! `RANDOM_SEED(PUT=seeds)` with every `INTEGER(4)` seed word equal to
//! `Simulation: Random Number Seed` (default 314159265). libgfortran (GCC 7
//! and later, `intrinsics/random.c`) keeps a xoshiro256** state. `PUT`
//! writes the seed words into that state XORed with fixed scrambling keys,
//! and `RANDOM_NUMBER` keeps the top 53 bits of each output (`rnumber_8`).
//!
//! Fortran's generator is global: ElmerSolver seeds it once and the stream
//! carries on from one IDR(s) solve to the next, so the n-th solve of a run
//! gets the n-th block of draws. The stream here is kept per thread
//! (libgfortran's is too): [`initialize_random_seed`] restarts it, which
//! `elmer-solver` does when a simulation is set up, and a thread that draws
//! before anything seeded it starts from the solver's `Random Number Seed`.

use std::cell::RefCell;

/// `xor_keys` of libgfortran's `random.c` (`scramble_seed`).
const XOR_KEYS: [u64; 4] = [
    0xbd0c_5b6e_50c2_df49,
    0xd460_61cd_46e1_df38,
    0xbb4f_4d4e_d610_3544,
    0x114a_583d_0756_ad39,
];

pub(crate) struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// The state after `RANDOM_SEED(PUT=[seed, seed, ...])`. Every seed
    /// word is the same, so their order inside the state does not matter.
    pub(crate) fn new(seed: u64) -> Self {
        let w = u64::from(seed as u32);
        let word = (w << 32) | w;
        Rng {
            s: XOR_KEYS.map(|k| word ^ k),
        }
    }

    /// xoshiro256** (`next` in `random.c`).
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

    /// Uniform in `[0, 1)`: the top 53 bits times 2^-64 (`rnumber_8`), which
    /// is `(v >> 11) * 2^-53` exactly.
    pub(crate) fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

thread_local! {
    static STREAM: RefCell<Option<Rng>> = const { RefCell::new(None) };
}

/// `InitializeRandomSeed` (ElmerSolver.F90:874): restarts this thread's
/// `RANDOM_NUMBER` stream from `seed` (`Simulation: Random Number Seed`,
/// 314159265 when not given).
pub fn initialize_random_seed(seed: i32) {
    STREAM.with(|s| *s.borrow_mut() = Some(Rng::new(u64::from(seed as u32))));
}

/// `EvenRandom()` (GeneralUtils.F90:2746): the next `RANDOM_NUMBER` of this
/// thread's stream, the same stream IDR(s) draws from, seeded with the
/// default `Random Number Seed` if nothing seeded it.
pub fn even_random() -> f64 {
    with_stream(314_159_265, |r| r.uniform())
}

/// `RANDOM_NUMBER` of a `REAL(4)` harvest (`random_r4`): the top 32 bits
/// of the next output cut to the 24 significand bits (`rnumber_4`), times
/// 2^-32; the same stream as [`even_random`].
pub fn even_random_r4() -> f32 {
    with_stream(314_159_265, |r| {
        let v = ((r.next_u64() >> 32) as u32) & (!0u32 << 8);
        v as f32 * (1.0 / 4_294_967_296.0)
    })
}

/// Runs `f` on this thread's stream, seeding it from `seed` first if
/// nothing has.
pub(crate) fn with_stream<R>(seed: u64, f: impl FnOnce(&mut Rng) -> R) -> R {
    STREAM.with(|s| f(s.borrow_mut().get_or_insert_with(|| Rng::new(seed))))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first six `RANDOM_NUMBER` values after ElmerSolver's default
    /// seeding, printed by a gfortran 13.3 program that does exactly what
    /// `InitializeRandomSeed` does (`RANDOM_SEED()`, then `PUT` of eight
    /// words of 314159265) and writes each value with `Z16.16`.
    /// `REAL(4)` and `REAL(8)` harvests take the same xoshiro output: the
    /// `REAL(4)` one is the `REAL(8)` one cut to 24 bits.
    #[test]
    fn r4_is_the_r8_draw_cut_to_24_bits() {
        initialize_random_seed(314_159_265);
        let r8: Vec<f64> = (0..4).map(|_| even_random()).collect();
        initialize_random_seed(314_159_265);
        for x in r8 {
            let want = (x * 16_777_216.0).floor() / 16_777_216.0;
            assert_eq!(f64::from(even_random_r4()), want);
        }
    }

    #[test]
    fn matches_gfortran_after_elmer_seeding() {
        let want: [u64; 6] = [
            0x3FC2_7402_B4E2_E778,
            0x3FE3_5088_1025_97B0,
            0x3FEE_436F_F054_21DB,
            0x3FC6_FF0C_4255_760C,
            0x3FE9_3507_4FAD_DA67,
            0x3FD1_C893_975E_B27E,
        ];
        let mut r = Rng::new(314_159_265);
        for w in want {
            assert_eq!(r.uniform().to_bits(), w);
        }
    }

    #[test]
    fn the_stream_carries_on_until_reseeded() {
        initialize_random_seed(314_159_265);
        let first = with_stream(1, |r| r.uniform());
        assert_eq!(first.to_bits(), 0x3FC2_7402_B4E2_E778);
        // Not the seed of `with_stream`: the stream is already seeded.
        assert_eq!(with_stream(1, |r| r.uniform()).to_bits(), 0x3FE3_5088_1025_97B0);
        initialize_random_seed(314_159_265);
        assert_eq!(with_stream(1, |r| r.uniform()), first);
    }
}
