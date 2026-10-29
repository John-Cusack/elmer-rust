//! Line (202) p basis: `H1Basis_LineNodal`, `H1Basis_dLineNodal`,
//! `H1Basis_LineBubbleP`, `H1Basis_dLineBubbleP`. Conventions: see
//! [`crate::types`].

use crate::legendre::{d_phi, phi};

/// `H1Basis_LineNodal`.
pub fn line_nodal(u: f64, val: &mut [f64], nb: &mut usize) {
    let c = 1.0 / 2.0;
    val[*nb] = c * (1.0 - u);
    val[*nb + 1] = c * (1.0 + u);
    *nb += 2;
}

/// `H1Basis_dLineNodal`.
pub fn d_line_nodal(_u: f64, grad: &mut [[f64; 3]], nb: &mut usize) {
    let c = 1.0 / 2.0;
    grad[*nb][0] = -c;
    grad[*nb + 1][0] = c;
    *nb += 2;
}

/// `H1Basis_LineBubbleP`: bubble (edge) functions `phi_{p+1}`, `p=1..pmax-1`.
/// `invert_edge` evaluates them at `-u`.
pub fn line_bubble_p(
    u: f64,
    pmax: i32,
    val: &mut [f64],
    nb: &mut usize,
    invert_edge: Option<bool>,
) {
    let invert = invert_edge.unwrap_or(false);
    for p in 1..pmax {
        let x = if invert { -u } else { u };
        val[*nb + p as usize - 1] = phi(p as usize + 1, x);
    }
    *nb += (pmax - 1).max(0) as usize;
}

/// `H1Basis_dLineBubbleP`. Fortran quirk: the inverted variant evaluates
/// `dPhi(p+1,-u)` without the chain-rule sign, exactly as H1Basis.F90 does.
pub fn d_line_bubble_p(
    u: f64,
    pmax: i32,
    grad: &mut [[f64; 3]],
    nb: &mut usize,
    invert_edge: Option<bool>,
) {
    let invert = invert_edge.unwrap_or(false);
    for p in 1..pmax {
        let x = if invert { -u } else { u };
        grad[*nb + p as usize - 1][0] = d_phi(p as usize + 1, x);
    }
    *nb += (pmax - 1).max(0) as usize;
}
