//! Element type descriptions read from Elmer's `elements.def`, with the
//! nodal basis coefficients built as `AddElementDescription`
//! (ElementBasis.F90:108) builds them.

use std::sync::OnceLock;

use crate::error::FeError;
use crate::lu::invert_matrix;

const ELEMENTS_DEF: &str = include_str!("../data/elements.def");

/// ElementBasis.F90:63, the per-direction width of the basis term numbering.
const MAX_DEG: usize = 4;

/// Elmer's `real**integer` (gfortran `__builtin_powi`, libgcc
/// `__powidf2`): binary powering, so the rounding matches bit for bit.
#[inline]
pub(crate) fn powi(x: f64, m: u32) -> f64 {
    let mut n = m;
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
    y
}

/// One monomial `coeff * u^p v^q w^r` of a nodal basis function.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BasisTerm {
    pub p: u8,
    pub q: u8,
    pub r: u8,
    pub coeff: f64,
}

/// `ElementType_t` restricted to what the Lagrange element library needs.
#[derive(Debug)]
pub struct ElementTypeDef {
    pub code: u16,
    pub dim: usize,
    pub n_nodes: usize,
    pub node_u: Vec<f64>,
    pub node_v: Vec<f64>,
    pub node_w: Vec<f64>,
    /// Default rule, `GaussPoints2` (relative order +1) and `GaussPoints0`
    /// (relative order -1), after Elmer's defaulting.
    pub gauss_points: usize,
    pub gauss_points2: usize,
    pub gauss_points0: usize,
    pub basis_degree: usize,
    /// `BasisFunctions(i)` for node `i`, terms in elements.def order.
    pub basis: Vec<Vec<BasisTerm>>,
    /// The type-level `StabilizationMK` (computed on first use, as
    /// `GetElementType` does; the elements.def value is discarded by Elmer).
    pub(crate) type_mk: OnceLock<f64>,
}

impl ElementTypeDef {
    pub fn family(&self) -> u16 {
        self.code / 100
    }
}

fn all_types() -> &'static [ElementTypeDef] {
    static TYPES: OnceLock<Vec<ElementTypeDef>> = OnceLock::new();
    TYPES.get_or_init(|| parse_elements_def(ELEMENTS_DEF).expect("bundled elements.def is valid"))
}

/// The description of element type `code`, as `GetElementType` finds it.
pub fn element_type(code: u16) -> Result<&'static ElementTypeDef, FeError> {
    all_types()
        .iter()
        .find(|t| t.code == code)
        .ok_or(FeError::UnknownType(code))
}

/// Every element type defined in the bundled elements.def.
pub fn element_types() -> impl Iterator<Item = &'static ElementTypeDef> {
    all_types().iter()
}

/// `ReadAndTrim` for elements.def: drop comments, lowercase, collapse blanks.
fn clean(line: &str) -> String {
    let s = line.split('!').next().unwrap_or("");
    s.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_ascii_lowercase()
}

fn nums<T: std::str::FromStr>(s: &str) -> Vec<T> {
    s.split_whitespace().map_while(|t| t.parse().ok()).collect()
}

/// `InitializeElementDescriptions` (ElementBasis.F90:430). The fields of
/// the Fortran `element` record persist from one ELEMENT block to the next;
/// only `GaussPoints2` and the basis terms are reset per block, so a block
/// that gives two Gauss counts keeps the previous block's `GaussPoints0`.
/// That carry-over is kept.
fn parse_elements_def(text: &str) -> Result<Vec<ElementTypeDef>, FeError> {
    let mut out = Vec::new();
    let (mut dim, mut code, mut nn) = (0usize, 0u16, 0usize);
    let (mut gp, mut gp0) = (0usize, 0usize);
    let mut lines = text.lines().map(clean).filter(|l| !l.is_empty());
    while let Some(l) = lines.next() {
        if !l.starts_with("element") {
            continue;
        }
        let mut terms: Vec<usize> = Vec::new();
        let mut node_uvw: [Option<Vec<f64>>; 3] = [None, None, None];
        let mut gp2 = 0usize;
        let mut gotit = false;
        for l in lines.by_ref() {
            if let Some(r) = l.strip_prefix("dimension") {
                dim = nums(r)[0];
            } else if let Some(r) = l.strip_prefix("code") {
                code = nums(r)[0];
            } else if let Some(r) = l.strip_prefix("nodes") {
                nn = nums(r)[0];
            } else if let Some(r) = l.strip_prefix("node u") {
                node_uvw[0] = Some(nums(r).into_iter().take(nn).collect());
            } else if let Some(r) = l.strip_prefix("node v") {
                node_uvw[1] = Some(nums(r).into_iter().take(nn).collect());
            } else if let Some(r) = l.strip_prefix("node w") {
                node_uvw[2] = Some(nums(r).into_iter().take(nn).collect());
            } else if let Some(r) = l.strip_prefix("basis") {
                terms = nums(r).into_iter().take(nn).collect();
            } else if let Some(r) = l.strip_prefix("gauss points") {
                let v: Vec<i64> = nums(r);
                gp2 = 0;
                if let Some(&a) = v.first() {
                    gp = a as usize;
                }
                if let Some(&b) = v.get(1) {
                    gp2 = b.max(0) as usize;
                }
                if let Some(&c) = v.get(2) {
                    gp0 = c.max(0) as usize;
                }
                if gp2 == 0 {
                    gp2 = gp;
                }
                if gp0 == 0 {
                    gp0 = gp;
                }
            } else if l == "end element" {
                gotit = true;
                break;
            }
        }
        if !gotit {
            continue;
        }
        let [node_u, node_v, node_w] = node_uvw.map(|c| c.unwrap_or_else(|| vec![0.0; nn]));
        let (basis, basis_degree) = if code >= 200 {
            add_element_description(code, dim, &node_u, &node_v, &node_w, &terms)?
        } else {
            (Vec::new(), 0)
        };
        out.push(ElementTypeDef {
            code,
            dim,
            n_nodes: nn,
            node_u,
            node_v,
            node_w,
            gauss_points: gp,
            gauss_points2: gp2,
            gauss_points0: gp0,
            basis_degree,
            basis,
            type_mk: OnceLock::new(),
        });
    }
    Ok(out)
}

fn term_powers(dim: usize, term: usize) -> (usize, usize, usize) {
    let k = term - 1;
    match dim {
        1 => (k, 0, 0),
        2 => (k % MAX_DEG, k / MAX_DEG, 0),
        _ => (
            k % MAX_DEG,
            (k / MAX_DEG) % MAX_DEG,
            k / (MAX_DEG * MAX_DEG),
        ),
    }
}

/// `AddElementDescription`: invert the node/term matrix; column `i` of the
/// inverse holds the coefficients of basis function `i`.
#[allow(clippy::type_complexity)]
fn add_element_description(
    code: u16,
    dim: usize,
    nu: &[f64],
    nv: &[f64],
    nw: &[f64],
    terms: &[usize],
) -> Result<(Vec<Vec<BasisTerm>>, usize), FeError> {
    let n = nu.len();
    let mut a = vec![0.0f64; n * n];
    let mut degree = 0usize;
    for i in 0..n {
        for j in 0..n {
            let (up, vp, wp) = term_powers(dim, terms[j]);
            let (u, v, w) = (nu[i], nv[i], nw[i]);
            let val = if dim == 1 {
                if u == 0.0 && up == 0 {
                    1.0
                } else {
                    powi(u, up as u32)
                }
            } else {
                let mut x = if up == 0 { 1.0 } else { powi(u, up as u32) };
                if vp != 0 {
                    x *= powi(v, vp as u32);
                }
                if dim == 3 && wp != 0 {
                    x *= powi(w, wp as u32);
                }
                x
            };
            a[i + j * n] = val;
            degree = degree.max(up).max(vp).max(wp);
        }
    }
    invert_matrix(&mut a, n).ok_or(FeError::SingularBasis(code))?;
    let basis = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let (p, q, r) = term_powers(dim, terms[j]);
                    BasisTerm {
                        p: p as u8,
                        q: q as u8,
                        r: r as u8,
                        coeff: a[j + i * n],
                    }
                })
                .collect()
        })
        .collect();
    Ok((basis, degree))
}
