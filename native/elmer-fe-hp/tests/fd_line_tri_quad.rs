//! Finite-difference consistency of the value/gradient routine pairs for the
//! line, triangle and quadrilateral p bases, and basic structural checks.

use elmer_fe_hp::line::*;
use elmer_fe_hp::maps::{edge_map, face_map};
use elmer_fe_hp::quad::*;
use elmer_fe_hp::triangle::*;

const H: f64 = 1e-6;
const N: usize = 64;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-6 * (1.0 + a.abs().max(b.abs()))
}

/// Central difference of `f(u,v)` (filling `val` and returning the count)
/// against `g(u,v)` (filling `grad`).
fn check2<F, G>(f: F, g: G, pts: &[(f64, f64)], n: usize)
where
    F: Fn(f64, f64, &mut [f64]) -> usize,
    G: Fn(f64, f64, &mut [[f64; 3]]) -> usize,
{
    for &(u, v) in pts {
        let mut gr = vec![[0.0; 3]; N];
        assert_eq!(g(u, v, &mut gr), n);
        let mut a = vec![0.0; N];
        let mut b = vec![0.0; N];
        let mut c = vec![0.0; N];
        let mut d = vec![0.0; N];
        assert_eq!(f(u + H, v, &mut a), n);
        assert_eq!(f(u - H, v, &mut b), n);
        assert_eq!(f(u, v + H, &mut c), n);
        assert_eq!(f(u, v - H, &mut d), n);
        for k in 0..n {
            let du = (a[k] - b[k]) / (2.0 * H);
            let dv = (c[k] - d[k]) / (2.0 * H);
            assert!(close(du, gr[k][0]), "d/du fn {k} at ({u},{v}): {du} vs {}", gr[k][0]);
            assert!(close(dv, gr[k][1]), "d/dv fn {k} at ({u},{v}): {dv} vs {}", gr[k][1]);
        }
    }
}

#[test]
fn line_bubbles_match_fd() {
    for p in 2..=6 {
        for &inv in &[false, true] {
            for &u in &[-0.7, -0.1, 0.35, 0.8] {
                let mut g = vec![[0.0; 3]; N];
                let mut nb = 0;
                d_line_bubble_p(u, p, &mut g, &mut nb, Some(inv));
                assert_eq!(nb, (p - 1) as usize);
                let mut a = vec![0.0; N];
                let mut b = vec![0.0; N];
                let (mut na, mut nbb) = (0, 0);
                line_bubble_p(u + H, p, &mut a, &mut na, Some(inv));
                line_bubble_p(u - H, p, &mut b, &mut nbb, Some(inv));
                for k in 0..nb {
                    // The Fortran omits the chain-rule sign when inverted.
                    let sign = if inv { -1.0 } else { 1.0 };
                    let fd = sign * (a[k] - b[k]) / (2.0 * H);
                    assert!(close(fd, g[k][0]), "p={p} inv={inv} k={k}: {fd} vs {}", g[k][0]);
                }
            }
        }
    }
}

#[test]
fn triangle_partition_of_unity() {
    let (mut v, mut nb) = (vec![0.0; N], 0);
    triangle_nodal_p(-0.3, 0.4, &mut v, &mut nb);
    assert_eq!(nb, 3);
    assert!((v[0] + v[1] + v[2] - 1.0).abs() < 1e-14);
    let (mut g, mut nb) = (vec![[0.0; 3]; N], 0);
    d_triangle_nodal_p(-0.3, 0.4, &mut g, &mut nb);
    for c in 0..2 {
        assert!((g[0][c] + g[1][c] + g[2][c]).abs() < 1e-14);
    }
}

fn tri_points() -> Vec<(f64, f64)> {
    // Interior of the equilateral triangle (-1,0),(1,0),(0,sqrt3).
    let s3 = 3f64.sqrt();
    vec![(-0.2, 0.4), (0.3, 0.7), (0.05, 0.3 * s3), (-0.45, 0.2)]
}

#[test]
fn triangle_edges_match_fd() {
    for p in 2..=6 {
        for dirs in [
            edge_map(303).unwrap(),
            [[2, 1], [3, 2], [1, 3], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2]],
        ] {
            let pm = [p, p, p];
            let n = 3 * (p - 1) as usize;
            check2(
                |u, v, out| {
                    let mut nb = 0;
                    triangle_edge_p(u, v, &pm, out, &mut nb, &dirs);
                    nb
                },
                |u, v, out| {
                    let mut nb = 0;
                    d_triangle_edge_p(u, v, &pm, out, &mut nb, &dirs);
                    nb
                },
                &tri_points(),
                n,
            );
        }
    }
    // Mixed degrees.
    let dirs = edge_map(303).unwrap();
    let pm = [3, 2, 5];
    check2(
        |u, v, out| {
            let mut nb = 0;
            triangle_edge_p(u, v, &pm, out, &mut nb, &dirs);
            nb
        },
        |u, v, out| {
            let mut nb = 0;
            d_triangle_edge_p(u, v, &pm, out, &mut nb, &dirs);
            nb
        },
        &tri_points(),
        2 + 1 + 4,
    );
}

#[test]
fn triangle_bubbles_match_fd() {
    let faces = [None, Some([1, 2, 3]), Some([2, 3, 1]), Some([3, 2, 1])];
    for p in 3..=7 {
        let n = ((p - 2) * (p - 1) / 2) as usize;
        for ln in faces {
            check2(
                |u, v, out| {
                    let mut nb = 0;
                    triangle_bubble_p(u, v, p, out, &mut nb, ln);
                    nb
                },
                |u, v, out| {
                    let mut nb = 0;
                    d_triangle_bubble_p(u, v, p, out, &mut nb, ln);
                    nb
                },
                &tri_points(),
                n,
            );
        }
    }
}

fn quad_points() -> Vec<(f64, f64)> {
    vec![(-0.6, 0.2), (0.3, -0.8), (0.55, 0.65), (-0.1, -0.35)]
}

/// Nodal values/gradients first, then `f` (as the Fortran callers do).
fn with_nodal_val(u: f64, v: f64, out: &mut [f64], f: impl FnOnce(&mut [f64], &mut usize)) -> usize {
    let mut nb = 0;
    quad_nodal(u, v, out, &mut nb);
    f(out, &mut nb);
    nb - 4
}

fn with_nodal_grad(u: f64, v: f64, out: &mut [[f64; 3]], f: impl FnOnce(&mut [[f64; 3]], &mut usize)) -> usize {
    let mut nb = 0;
    d_quad_nodal(u, v, out, &mut nb);
    f(out, &mut nb);
    nb - 4
}

fn shifted<F: Fn(f64, f64, &mut [f64]) -> usize>(f: F) -> impl Fn(f64, f64, &mut [f64]) -> usize {
    // check2 compares indices from 0; drop the 4 nodal slots.
    move |u, v, out: &mut [f64]| {
        let mut tmp = vec![0.0; out.len() + 4];
        let n = f(u, v, &mut tmp);
        out[..n].copy_from_slice(&tmp[4..4 + n]);
        n
    }
}

fn shifted_g<G: Fn(f64, f64, &mut [[f64; 3]]) -> usize>(g: G) -> impl Fn(f64, f64, &mut [[f64; 3]]) -> usize {
    move |u, v, out: &mut [[f64; 3]]| {
        let mut tmp = vec![[0.0; 3]; out.len() + 4];
        let n = g(u, v, &mut tmp);
        out[..n].copy_from_slice(&tmp[4..4 + n]);
        n
    }
}

#[test]
fn quad_nodal_partition_of_unity() {
    let (mut v, mut nb) = (vec![0.0; N], 0);
    quad_nodal(0.2, -0.5, &mut v, &mut nb);
    assert_eq!(nb, 4);
    assert!((v[..4].iter().sum::<f64>() - 1.0).abs() < 1e-14);
    let (mut g, mut nb) = (vec![[0.0; 3]; N], 0);
    d_quad_nodal(0.2, -0.5, &mut g, &mut nb);
    for c in 0..2 {
        assert!(g[..4].iter().map(|x| x[c]).sum::<f64>().abs() < 1e-14);
    }
}

#[test]
fn quad_edges_match_fd() {
    let maps = [
        edge_map(404).unwrap(),
        [[2, 1], [3, 2], [3, 4], [4, 1], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2], [0; 2]],
    ];
    for p in 2..=6 {
        let pm = [p; 4];
        let n = 4 * (p - 1) as usize;
        for dirs in maps {
            check2(
                shifted(|u, v, out| {
                    with_nodal_val(u, v, out, |o, nb| quad_edge_p(u, v, &pm, o, nb, &dirs))
                }),
                shifted_g(|u, v, out| {
                    with_nodal_grad(u, v, out, |o, nb| d_quad_edge_p(u, v, &pm, o, nb, &dirs))
                }),
                &quad_points(),
                n,
            );
            check2(
                shifted(|u, v, out| with_nodal_val(u, v, out, |o, nb| sd_quad_edge_p(u, v, &pm, o, nb, &dirs))),
                shifted_g(|u, v, out| with_nodal_grad(u, v, out, |o, nb| sd_d_quad_edge_p(u, v, &pm, o, nb, &dirs))),
                &quad_points(),
                n,
            );
        }
    }
}

#[test]
fn quad_bubbles_match_fd() {
    let faces = [None, Some([1, 2, 3, 4]), Some([2, 3, 4, 1]), Some([4, 3, 2, 1]), Some([3, 4, 1, 2])];
    for p in 2..=6 {
        let n = ((p - 1) * (p - 1)) as usize;
        for ln in faces {
            check2(
                shifted(|u, v, out| with_nodal_val(u, v, out, |o, nb| quad_bubble_p(u, v, p, o, nb, ln))),
                shifted_g(|u, v, out| with_nodal_grad(u, v, out, |o, nb| d_quad_bubble_p(u, v, p, o, nb, ln))),
                &quad_points(),
                n,
            );
        }
        let sd = ((p - 2) * (p - 3) / 2).max(0) as usize;
        for ln in faces {
            check2(
                shifted(|u, v, out| with_nodal_val(u, v, out, |o, nb| sd_quad_bubble_p(u, v, p, o, nb, ln))),
                shifted_g(|u, v, out| with_nodal_grad(u, v, out, |o, nb| sd_d_quad_bubble_p(u, v, p, o, nb, ln))),
                &quad_points(),
                sd,
            );
        }
    }
    let _ = face_map(404);
}
