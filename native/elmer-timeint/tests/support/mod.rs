//! Deterministic test data (no external RNG crate).

pub struct Lcg(pub u64);

impl Lcg {
    pub fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        // Uniform in (-1, 1), with a full 53-bit mantissa.
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }

    pub fn vec(&mut self, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.next_f64()).collect()
    }
}

pub fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

#[track_caller]
pub fn assert_bits(what: &str, got: &[f64], want: &[f64]) {
    assert!(
        same_bits(got, want),
        "{what}: not bitwise equal\n got  {got:?}\n want {want:?}"
    );
}
