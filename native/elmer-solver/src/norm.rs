//! `ComputeNorm` (SolverBasics.F90:3164), operation for operation, so the
//! result is bitwise equal to the Fortran's.
//!
//! gfortran evaluates `SUM` left to right in one accumulator, `x**2` as
//! `x*x` and `x**NormDim` (integer variable) through libgcc `__powidf2`;
//! all three are reproduced here.

use std::borrow::Cow;

use elmer_core::{Communicator, Logger, Params, NO_DOF};

/// The norm keywords of one solver (`Nonlinear System Norm *`).
#[derive(Clone, Debug, PartialEq)]
pub struct NormSpec {
    /// `Nonlinear System Norm Degree`: 0 is the max norm; default 2.
    pub degree: i32,
    /// 1-based components in the norm (`Nonlinear System Norm Components`,
    /// or `1..=Norm Dofs`); `None` means all `dofs`.
    pub components: Option<Vec<usize>>,
    pub dofs: usize,
    /// `Norm Permutation`: 1-based first perm entry, used when > 1.
    pub perm_start: Option<usize>,
    /// `Nonlinear System Nodal Norm`.
    pub nodal: bool,
    /// `Nonlinear System Consistent Norm` (parallel only).
    pub consistent: bool,
}

impl NormSpec {
    /// Plain L2 over all dofs.
    pub fn l2(dofs: usize) -> Self {
        Self {
            degree: 2,
            components: None,
            dofs,
            perm_start: None,
            nodal: false,
            consistent: false,
        }
    }

    pub fn from_params(p: &dyn Params, dofs: usize) -> Self {
        let degree = p.integer("nonlinear system norm degree").unwrap_or(2);
        let components = match p.get("nonlinear system norm components") {
            Some(elmer_core::Value::Integer(v)) => Some(v.iter().map(|&c| c as usize).collect()),
            _ => p
                .integer("nonlinear system norm dofs")
                .map(|nd| (1..=nd.max(0) as usize).collect()),
        };
        Self {
            degree,
            components,
            dofs,
            perm_start: p.integer("norm permutation").filter(|&s| s > 1).map(|s| s as usize),
            nodal: p.logical("nonlinear system nodal norm").unwrap_or(false),
            consistent: p.logical("nonlinear system consistent norm").unwrap_or(false),
        }
    }

    fn norm_dofs(&self) -> usize {
        self.components.as_ref().map_or(self.dofs, Vec::len)
    }
}

/// libgcc `__powidf2`: `x**m` for an integer `m` by binary powering.
pub fn powi(x: f64, m: i32) -> f64 {
    let mut n = m.unsigned_abs();
    let mut x = x;
    let mut y = if n % 2 == 1 { x } else { 1.0 };
    loop {
        n >>= 1;
        if n == 0 {
            break;
        }
        x *= x;
        if n % 2 == 1 {
            y *= x;
        }
    }
    if m < 0 {
        1.0 / y
    } else {
        y
    }
}

/// `ComputeNorm(Solver, n, x)` after the `Norm Permutation` and nodal-norm
/// selection: the norm of the first `n` entries of `x` (component selection
/// strides over all of `x`, as `x(j::Dofs)` does).
///
/// In parallel (`comm.size() > 1`) with [`NormSpec::consistent`], use
/// [`compute_norm_consistent`] instead; this function then counts shared
/// dofs once per partition, as the non-consistent Fortran branch does.
pub fn compute_norm(spec: &NormSpec, x: &[f64], n: usize, comm: &dyn Communicator) -> f64 {
    let parallel = comm.size() > 1;
    let dofs = spec.dofs;
    let d = spec.degree;
    let norm_dofs = spec.norm_dofs();

    if norm_dofs < dofs {
        let comps = spec.components.as_deref().unwrap_or(&[]);
        let totn = reduce_sum(comm, n as f64);
        if totn == 0.0 {
            return 0.0;
        }
        // NormDOFs*totn is an integer product, then divided by real DOFs.
        let nscale = (norm_dofs as f64 * totn) / dofs as f64;
        let comp = |j: usize| x.iter().skip(j - 1).step_by(dofs).copied();
        let mut norm = 0.0f64;
        match d {
            0 => {
                for &j in comps {
                    let m = comp(j).map(f64::abs).fold(f64::MIN, f64::max);
                    norm = norm.max(m);
                }
                if parallel {
                    norm = reduce_max(comm, norm);
                }
            }
            1 => {
                for &j in comps {
                    norm += seq_sum(comp(j).map(f64::abs));
                }
                norm = reduce_sum(comm, norm) / nscale;
            }
            2 => {
                for &j in comps {
                    norm += seq_sum(comp(j).map(|v| v * v));
                }
                norm = (reduce_sum(comm, norm) / nscale).sqrt();
            }
            _ => {
                for &j in comps {
                    norm += seq_sum(comp(j).map(|v| powi(v, d)));
                }
                norm = (reduce_sum(comm, norm) / nscale).powf(1.0 / d as f64);
            }
        }
        norm
    } else if parallel {
        let totn = reduce_sum(comm, n as f64);
        if totn == 0.0 {
            return 0.0;
        }
        let nscale = totn;
        let x = &x[..n];
        match d {
            0 => {
                let val = if n > 0 { x.iter().map(|v| v.abs()).fold(f64::MIN, f64::max) } else { 0.0 };
                reduce_max(comm, val)
            }
            1 => reduce_sum(comm, seq_sum(x.iter().map(|v| v.abs()))) / nscale,
            2 => (reduce_sum(comm, seq_sum(x.iter().map(|v| v * v))) / nscale).sqrt(),
            _ => (reduce_sum(comm, seq_sum(x.iter().map(|&v| powi(v, d)))) / nscale).powf(1.0 / d as f64),
        }
    } else {
        if n == 0 {
            return 0.0;
        }
        let x = &x[..n];
        let nf = n as f64;
        match d {
            0 => x.iter().map(|v| v.abs()).fold(f64::MIN, f64::max),
            1 => seq_sum(x.iter().map(|v| v.abs())) / nf,
            2 => (seq_sum(x.iter().map(|v| v * v)) / nf).sqrt(),
            // SUM((x(1:n)**NormDim)/n): each term is divided before summing.
            _ => seq_sum(x.iter().map(|&v| powi(v, d) / nf)).powf(1.0 / d as f64),
        }
    }
}

/// The consistent parallel norm: only entries `j` with `owned(row_of(j))`
/// count, so each shared dof is counted once. `row_of` maps a compacted
/// index back to the matrix row (identity without `Norm Permutation`).
pub fn compute_norm_consistent(
    spec: &NormSpec,
    x: &[f64],
    n: usize,
    owned: &dyn Fn(usize) -> bool,
    comm: &dyn Communicator,
) -> f64 {
    let dofs = spec.dofs;
    let d = spec.degree;
    let passive: Vec<bool> = match &spec.components {
        Some(c) if c.len() < dofs => (0..dofs).map(|k| !c.contains(&(k + 1))).collect(),
        _ => vec![false; dofs],
    };
    let mut norm = 0.0f64;
    let mut totn = 0usize;
    for (j, &v) in x[..n].iter().enumerate() {
        if passive[j % dofs] || !owned(j) {
            continue;
        }
        match d {
            0 => norm = norm.max(v.abs()),
            1 => norm += v.abs(),
            2 => norm += v * v,
            _ => norm += powi(v, d),
        }
        totn += 1;
    }
    let totn = reduce_sum(comm, totn as f64);
    if totn == 0.0 {
        return 0.0;
    }
    match d {
        0 => reduce_max(comm, norm),
        1 => reduce_sum(comm, norm) / totn,
        2 => (reduce_sum(comm, norm) / totn).sqrt(),
        _ => (reduce_sum(comm, norm) / totn).powf(1.0 / d as f64),
    }
}

/// The `Norm Permutation` compaction and `Nonlinear System Nodal Norm`
/// truncation that `ComputeNorm` applies before reducing. `perm` is the
/// variable's 0-based perm ([`NO_DOF`] = absent). Returns the values to
/// reduce and their count.
pub fn select_norm_entries<'a>(
    spec: &NormSpec,
    x: &'a [f64],
    n: usize,
    perm: Option<&[u32]>,
    n_nodes: usize,
    log: &dyn Logger,
) -> (Cow<'a, [f64]>, usize) {
    let mut n = n;
    let mut vals: Cow<'a, [f64]> = Cow::Borrowed(x);
    // Fortran perm values: slot + 1, or 0 when absent.
    let pf = |p: u32| if p == NO_DOF { 0 } else { p as usize + 1 };
    if let (Some(start), Some(perm)) = (spec.perm_start, perm) {
        let y: Vec<f64> = perm[start - 1..]
            .iter()
            .filter(|&&p| p != NO_DOF)
            .map(|&p| x[p as usize])
            .collect();
        n = y.len();
        vals = Cow::Owned(y);
    }
    if spec.nodal {
        if let Some(perm) = perm {
            let nn = n_nodes.min(perm.len());
            let i = perm[..nn].iter().map(|&p| pf(p)).max().unwrap_or(0);
            let j = perm[nn..].iter().map(|&p| pf(p)).min().unwrap_or(usize::MAX);
            if j == 0 {
            } else if j > i {
                n = spec.dofs * i;
                log.info("ComputeNorm", "Considering only the nodal entries in norm computation!", 7);
            } else {
                log.warn(
                    "ComputeNorm",
                    "Nodal norm is only available if all the p-dofs follow the last nodal dof!",
                );
            }
        }
    }
    (vals, n)
}

/// Left-to-right sum in one accumulator, as gfortran's inline `SUM`.
fn seq_sum(it: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for v in it {
        s += v;
    }
    s
}

fn reduce_sum(comm: &dyn Communicator, v: f64) -> f64 {
    let mut b = [v];
    comm.sum_f64(&mut b);
    b[0]
}

fn reduce_max(comm: &dyn Communicator, v: f64) -> f64 {
    let mut b = [v];
    comm.max_f64(&mut b);
    b[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::SerialComm;

    #[test]
    fn powi_matches_binary_powering() {
        assert_eq!(powi(3.0, 0), 1.0);
        assert_eq!(powi(-2.0, 3), -8.0);
        assert_eq!(powi(2.0, -2), 0.25);
        assert_eq!(powi(1.1, 5), ((1.1 * 1.1) * (1.1 * 1.1)) * 1.1);
    }

    #[test]
    fn serial_norm_degrees() {
        let x = [3.0, -4.0, 0.0, 1.0];
        let mut s = NormSpec::l2(1);
        assert_eq!(compute_norm(&s, &x, 4, &SerialComm), (26.0f64 / 4.0).sqrt());
        s.degree = 0;
        assert_eq!(compute_norm(&s, &x, 4, &SerialComm), 4.0);
        s.degree = 1;
        assert_eq!(compute_norm(&s, &x, 4, &SerialComm), 2.0);
        s.degree = 4;
        let want = (81.0 / 4.0 + 256.0 / 4.0 + 0.0 + 1.0 / 4.0f64).powf(1.0 / 4.0);
        assert_eq!(compute_norm(&s, &x, 4, &SerialComm), want);
        assert_eq!(compute_norm(&s, &x, 0, &SerialComm), 0.0);
    }

    #[test]
    fn component_norm() {
        // dofs = 2, only component 2 counts; nscale = 1*4/2 = 2.
        let x = [100.0, 3.0, 100.0, 4.0];
        let s = NormSpec {
            components: Some(vec![2]),
            ..NormSpec::l2(2)
        };
        assert_eq!(compute_norm(&s, &x, 4, &SerialComm), (25.0f64 / 2.0).sqrt());
    }

    #[test]
    fn norm_permutation_and_nodal() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let perm = [3, NO_DOF, 0, 1, 2];
        let s = NormSpec {
            perm_start: Some(3),
            ..NormSpec::l2(1)
        };
        let (v, n) = select_norm_entries(&s, &x, 4, Some(&perm), 5, &elmer_core::NullLogger);
        assert_eq!((&v[..], n), (&[1.0, 2.0, 3.0][..], 3));
        let s = NormSpec {
            nodal: true,
            ..NormSpec::l2(1)
        };
        let (_, n) = select_norm_entries(&s, &x, 4, Some(&perm), 3, &elmer_core::NullLogger);
        // nodes 0..3 reach slot 3 (Fortran 4); the rest start at Fortran 1 < 4: warn, keep n.
        assert_eq!(n, 4);
        let perm = [0, 1, 2, 3];
        let (_, n) = select_norm_entries(&s, &x, 4, Some(&perm), 2, &elmer_core::NullLogger);
        assert_eq!(n, 2);
    }
}
