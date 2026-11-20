//! Replication of Elmer's `ComputeNorm` (`SolverBasics.F90`) and the
//! `CompareToReferenceSolution` verdict (`ElmerSolver.F90`).
//!
//! The grader recomputes norms from captured fields in exactly the order
//! Fortran accumulates them (ascending index, sequential `+`), because the
//! bitwise instrumentation-invariance check (§4.3 item 2) needs the
//! recomputed norm to equal the solver's internal value.

use crate::{OracleError, Result};

/// Options mirroring the solver keywords that change norm semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct NormOpts {
    /// `Nonlinear System Norm Degree` (default 2).
    pub degree: i32,
    /// `Nonlinear System Norm Dofs`; defaults to the variable's DOFs.
    pub norm_dofs: Option<usize>,
    /// `Nonlinear System Norm Components`; 1-based component subset.
    pub components: Option<Vec<usize>>,
    /// `Nonlinear System Nodal Norm`: restrict to nodal entries.
    pub nodal: bool,
    /// `Norm Permutation`: first permutation index kept (1-based, Fortran).
    pub perm_start: Option<usize>,
}

impl Default for NormOpts {
    fn default() -> Self {
        NormOpts {
            degree: 2,
            norm_dofs: None,
            components: None,
            nodal: false,
            perm_start: None,
        }
    }
}

/// Compute the norm of `values` (length `n = dofs * n_nodes`) with the
/// permutation `perm` (`perm[i] > 0` is the 1-based value index of node `i`,
/// `0` means absent), exactly as the serial branch of `ComputeNorm` does.
///
/// `n_mesh_nodes` is `Solver % Mesh % NumberOfNodes`, needed only for the
/// nodal-norm truncation. Parallel branches (`ConsistentNorm`, `Parallel`)
/// are not replicated: golden capture always runs serially, and MPI grading
/// compares per-partition files through the same serial formula.
pub fn compute_norm(
    values: &[f64],
    perm: Option<&[i32]>,
    dofs: usize,
    n_mesh_nodes: usize,
    opts: &NormOpts,
) -> Result<f64> {
    if dofs == 0 {
        return Err(OracleError::Parse(String::from("dofs must be positive")));
    }
    // Compacted index map when `Norm Permutation` skips leading entries.
    let compact: Vec<usize>;
    let (x, active_perm): (&[f64], Option<&[i32]>) = match opts.perm_start {
        Some(start) if start > 1 => {
            let p = perm.ok_or_else(|| {
                OracleError::Parse(String::from("Norm Permutation needs a permutation"))
            })?;
            compact = p
                .iter()
                .skip(start - 1)
                .filter_map(|&q| if q > 0 { Some(q as usize - 1) } else { None })
                .collect();
            let gathered: Vec<f64> = compact.iter().map(|&k| values[k]).collect();
            // The gathered array is owned here; leak the order through a
            // second pass below instead of re-borrowing.
            return norm_core(&gathered, None, dofs, n_mesh_nodes, opts, values.len());
        }
        _ => (values, perm),
    };
    let _ = active_perm;
    norm_core(x, perm, dofs, n_mesh_nodes, opts, values.len())
}

fn norm_core(
    x: &[f64],
    perm: Option<&[i32]>,
    dofs: usize,
    n_mesh_nodes: usize,
    opts: &NormOpts,
    full_len: usize,
) -> Result<f64> {
    let n = x.len();
    if n == 0 {
        return Ok(0.0);
    }
    // Fortran: explicit Components set NormDofs to their count; an explicit
    // Norm Dofs is used as given; otherwise NormDofs defaults to Dofs.
    let norm_dofs = match (opts.norm_dofs, &opts.components) {
        (Some(n), _) => n,
        (None, Some(c)) => c.len(),
        (None, None) => dofs,
    };
    if norm_dofs == 0 || norm_dofs > dofs {
        return Err(OracleError::Parse(format!(
            "norm dofs {norm_dofs} out of range for dofs {dofs}"
        )));
    }
    // Nodal truncation: `n = Dofs * i` where `i = MAXVAL(Perm(1:nn))`, only
    // when every node entry precedes every p/bubble entry.
    let mut n_eff = n;
    if opts.nodal {
        if let Some(p) = perm {
            let nn = n_mesh_nodes.min(p.len());
            let imax = p[..nn].iter().max().copied().unwrap_or(0);
            let jmin = p[nn..].iter().min().copied().unwrap_or(0);
            if jmin == 0 {
                // Leave n_eff unchanged, as Fortran does.
            } else if jmin > imax {
                n_eff = (dofs * imax as usize).min(n);
            } else {
                return Err(OracleError::GradeFailed(String::from(
                    "nodal norm without node-first ordering",
                )));
            }
        }
    }
    let passive: Vec<bool> = match &opts.components {
        None => vec![false; dofs],
        Some(comps) => {
            let mut passive = vec![true; dofs];
            for &c in comps {
                if c == 0 || c > dofs {
                    return Err(OracleError::Parse(format!(
                        "norm component {c} out of range for dofs {dofs}"
                    )));
                }
                passive[c - 1] = false;
            }
            passive
        }
    };
    let _ = full_len;

    if norm_dofs < dofs {
        // Strided per-component accumulation over the full array, scaled by
        // nscale = NormDofs * n / Dofs where n is the (possibly nodally
        // truncated) length: degree 1 divides, higher degrees take the root.
        let comps: Vec<usize> = match &opts.components {
            Some(c) => c.clone(),
            None => (1..=norm_dofs).collect(),
        };
        // Fortran sums the strided sections to the end of the array; only
        // the scale uses the truncated length.
        let full = x.len();
        let nscale = norm_dofs as f64 * n_eff as f64 / dofs as f64;
        match opts.degree {
            0 => {
                let mut m = 0.0f64;
                for &c in &comps {
                    let mut k = c - 1;
                    while k < full {
                        let v = x[k].abs();
                        if v > m {
                            m = v;
                        }
                        k += dofs;
                    }
                }
                return Ok(m);
            }
            // `Norm = Norm + SUM(...)`: each component sums on its own,
            // then adds to the running norm.
            1 => {
                let mut acc = 0.0;
                for &c in &comps {
                    acc += x[c - 1..full].iter().step_by(dofs).fold(0.0, |s, v| s + v.abs());
                }
                return Ok(acc / nscale);
            }
            d => {
                let mut acc = 0.0;
                for &c in &comps {
                    acc += x[c - 1..full].iter().step_by(dofs).fold(0.0, |s, v| s + v.powi(d));
                }
                return Ok((acc / nscale).powf(1.0 / d as f64));
            }
        }
    }

    // Full-width accumulation in ascending index order.
    let passive_any = passive.iter().any(|&b| b);
    match opts.degree {
        0 => {
            let mut m = 0.0f64;
            for (j, &v) in x[..n_eff].iter().enumerate() {
                if passive_any && passive[j % dofs] {
                    continue;
                }
                let a = v.abs();
                if a > m {
                    m = a;
                }
            }
            Ok(m)
        }
        1 => {
            let mut s = 0.0;
            for (j, &v) in x[..n_eff].iter().enumerate() {
                if passive_any && passive[j % dofs] {
                    continue;
                }
                s += v.abs();
            }
            Ok(s / n_eff as f64)
        }
        2 => {
            let mut s = 0.0;
            for (j, &v) in x[..n_eff].iter().enumerate() {
                if passive_any && passive[j % dofs] {
                    continue;
                }
                s += v * v;
            }
            Ok((s / n_eff as f64).sqrt())
        }
        d => {
            let mut s = 0.0;
            for (j, &v) in x[..n_eff].iter().enumerate() {
                if passive_any && passive[j % dofs] {
                    continue;
                }
                s += v.powi(d);
            }
            Ok((s / n_eff as f64).powf(1.0 / d as f64))
        }
    }
}

/// How `MPI_ALLREDUCE(MPI_SUM)` groups the ranks' partial sums.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReduceOrder {
    /// Open MPI's recursive doubling (coll/tuned for small messages): the
    /// first `2·(p − 2^k)` ranks combine in pairs, then a doubling tree
    /// over the `2^k` survivors in rank order.
    RecursiveDoubling,
    /// Left to right: `((s0 + s1) + s2) + …`.
    Sequential,
}

/// Sum of `parts` in the grouping of `order`.
pub fn reduce_sum(parts: &[f64], order: ReduceOrder) -> f64 {
    match order {
        ReduceOrder::Sequential => match parts.split_first() {
            Some((&first, rest)) => rest.iter().fold(first, |a, &b| a + b),
            None => 0.0,
        },
        ReduceOrder::RecursiveDoubling => {
            let p = parts.len();
            if p == 0 {
                return 0.0;
            }
            let pow2 = 1usize << (usize::BITS - 1 - p.leading_zeros());
            let extra = p - pow2;
            let mut level: Vec<f64> = Vec::with_capacity(pow2);
            for r in 0..2 * extra {
                if r % 2 == 1 {
                    level.push(parts[r - 1] + parts[r]);
                }
            }
            level.extend_from_slice(&parts[2 * extra..]);
            while level.len() > 1 {
                level = level.chunks(2).map(|c| c[0] + c[1]).collect();
            }
            level[0]
        }
    }
}

/// The parallel branch of `ComputeNorm` (no `Nonlinear System Consistent
/// Norm`): every rank sums its own local values, interface nodes included,
/// and the sums and counts are reduced over the ranks. `parts` are the
/// ranks' `Values` in rank order. Norm options are the default ones and
/// `norm_dofs`.
pub fn parallel_norm(parts: &[Vec<f64>], dofs: usize, opts: &NormOpts, order: ReduceOrder) -> Result<f64> {
    if dofs == 0 || opts.components.is_some() || opts.nodal || opts.perm_start.is_some() {
        return Err(OracleError::Unsupported(String::from("parallel norm with these options")));
    }
    let norm_dofs = opts.norm_dofs.unwrap_or(dofs);
    if norm_dofs == 0 || norm_dofs > dofs {
        return Err(OracleError::Parse(format!("norm dofs {norm_dofs} out of range for dofs {dofs}")));
    }
    let term = |v: f64| match opts.degree {
        1 => v.abs(),
        2 => v * v,
        d => v.powi(d),
    };
    let totn: usize = parts.iter().map(Vec::len).sum();
    if totn == 0 {
        return Ok(0.0);
    }
    if opts.degree == 0 {
        let m = parts.iter().flat_map(|x| x.iter()).fold(0.0f64, |m, v| m.max(v.abs()));
        return Ok(m);
    }
    // Each rank's local value before the reduction, accumulated as the
    // serial branches do.
    let local: Vec<f64> = parts
        .iter()
        .map(|x| {
            if norm_dofs < dofs {
                let mut acc = 0.0;
                for c in 0..norm_dofs {
                    acc += x.get(c..).unwrap_or(&[]).iter().step_by(dofs).fold(0.0, |s, &v| s + term(v));
                }
                acc
            } else {
                x.iter().fold(0.0, |s, &v| s + term(v))
            }
        })
        .collect();
    let total = reduce_sum(&local, order);
    let nscale = if norm_dofs < dofs { norm_dofs as f64 * totn as f64 / dofs as f64 } else { totn as f64 };
    Ok(match opts.degree {
        1 => total / nscale,
        2 => (total / nscale).sqrt(),
        d => (total / nscale).powf(1.0 / d as f64),
    })
}

/// Verdict replicating `CompareToReferenceSolution` for one solver:
/// absolute error when `absolute` is set (or the reference is sub-tiny),
/// relative error otherwise; NaN always fails.
pub fn compare_reference(norm: f64, reference: f64, tol: f64, absolute: bool) -> bool {
    let use_absolute = absolute || reference < f64::MIN_POSITIVE;
    let err = (norm - reference).abs();
    let err = if use_absolute {
        err
    } else {
        err / reference
    };
    !(err > tol || err.is_nan())
}

/// Print a float the way the grader writes `Reference Norm` values:
/// 17 significant digits, so the `.sif` parser round-trips the exact double.
pub fn format_grader_norm(value: f64) -> String {
    format!("{value:.17e}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rms_matches_hand_computation() {
        // Values 1..=4, degree 2: sqrt((1+4+9+16)/4) = sqrt(7.5).
        let v = vec![1.0, 2.0, 3.0, 4.0];
        let n = compute_norm(&v, None, 1, 4, &NormOpts::default()).unwrap();
        assert!((n - 7.5f64.sqrt()).abs() < 1e-15);
    }

    #[test]
    fn accumulation_order_is_sequential() {
        // 1e16 followed by ten 1.0s with degree 1: 1e16 has ulp 2, so a
        // plain ascending accumulation loses every +1 and the mean is
        // exactly 1e16/11. Compensated summation would keep the ten ones
        // and give (1e16+10)/11 instead, so this pins the
        // Fortran-matching naive sequential order against the realistic
        // future "optimisation" of adding compensation.
        let mut v = vec![1e16];
        v.extend(std::iter::repeat_n(1.0, 10));
        let opts = NormOpts { degree: 1, ..Default::default() };
        let n = compute_norm(&v, None, 1, 11, &opts).unwrap();
        assert_eq!(n, 1e16 / 11.0);
    }

    #[test]
    fn degree_zero_is_max_abs() {
        let v = vec![-3.0, 1.0, 2.5];
        let opts = NormOpts {
            degree: 0,
            ..Default::default()
        };
        assert_eq!(compute_norm(&v, None, 1, 3, &opts).unwrap(), 3.0);
    }

    #[test]
    fn degree_one_is_mean_abs() {
        let v = vec![-3.0, 1.0, 2.0];
        let opts = NormOpts {
            degree: 1,
            ..Default::default()
        };
        assert!((compute_norm(&v, None, 1, 3, &opts).unwrap() - 2.0).abs() < 1e-15);
    }

    #[test]
    fn component_subset_scales_like_fortran() {
        // Dofs = 2, only component 1 active: Components sets NormDofs to 1,
        // so SUM(x(1::2)**2) = 25 over nscale = 1*4/2 = 2.
        let v = vec![3.0, 100.0, 4.0, 200.0];
        let opts = NormOpts {
            components: Some(vec![1]),
            ..Default::default()
        };
        let n = compute_norm(&v, None, 2, 2, &opts).unwrap();
        let expected = (25.0f64 / 2.0).sqrt();
        assert!((n - expected).abs() < 1e-15, "got {n}, want {expected}");
    }

    #[test]
    fn strided_norm_sums_each_component_on_its_own() {
        // Fortran: Norm = (0 + SUM(ABS(x(1::3)))) + SUM(ABS(x(2::3))).
        // Component 1 sums to 2e16 (ulp 4). Component 2's own SUM is 3,
        // and 2e16 + 3 rounds to 2e16 + 4; one running sum would add 1.5
        // twice and round each away, leaving 2e16.
        let v = vec![1e16, 1.5, 9.0, 1e16, 1.5, 9.0];
        let opts = NormOpts { degree: 1, norm_dofs: Some(2), ..Default::default() };
        let n = compute_norm(&v, None, 3, 2, &opts).unwrap();
        assert_eq!(n, (2e16 + 4.0) / 4.0);
        assert_ne!(n, 2e16 / 4.0);
    }

    #[test]
    fn recursive_doubling_groups_like_open_mpi() {
        // 1e16 has ulp 2, so the grouping decides whether the 1.5s survive.
        // p = 3: the first pair combines, then the third joins.
        let s = [1.5, 1.5, 1e16];
        assert_eq!(reduce_sum(&s, ReduceOrder::RecursiveDoubling), 3.0 + 1e16);
        let f = [1.5, 1.5, 1e16, 0.0];
        // p = 4: (f0 + f1) + (f2 + f3) keeps the 3.0.
        assert_eq!(reduce_sum(&f, ReduceOrder::RecursiveDoubling), 3.0 + 1e16);
        assert_eq!(reduce_sum(&f, ReduceOrder::Sequential), ((1.5 + 1.5) + 1e16) + 0.0);
        let six = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(reduce_sum(&six, ReduceOrder::RecursiveDoubling), ((1.0 + 2.0) + (3.0 + 4.0)) + (5.0 + 6.0));
    }

    #[test]
    fn parallel_norm_counts_every_rank_s_values() {
        let parts = vec![vec![1.0, 2.0], vec![2.0, 3.0]];
        let n = parallel_norm(&parts, 1, &NormOpts::default(), ReduceOrder::Sequential).unwrap();
        assert_eq!(n, ((1.0f64 + 4.0) + (4.0 + 9.0)) .sqrt() / 2.0);
    }

    #[test]
    fn tiny_reference_forces_absolute_error() {
        assert!(compare_reference(1e-310, 0.0, 1e-5, false));
        assert!(!compare_reference(1e-3, 0.0, 1e-5, true));
        assert!(!compare_reference(f64::NAN, 1.0, 1e-5, false));
    }

    #[test]
    fn grader_format_round_trips() {
        for v in [0.1f64, 1.0 / 3.0, 7.5f64.sqrt(), 1e-300, 123456.789] {
            let s = format_grader_norm(v);
            let back: f64 = s.parse().unwrap();
            assert_eq!(back, v, "format {s} did not round-trip");
        }
    }
}
