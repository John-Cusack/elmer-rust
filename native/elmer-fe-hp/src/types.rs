//! Shared conventions for the family modules.
//!
//! Every routine here is a scalar (single point) transliteration of the
//! `nvec`-blocked Fortran routine of the same name in `H1Basis.F90`:
//!
//! * `H1Basis_Foo(nvec,u,v,w,..,fval,nbasis)` becomes `foo(u,v,w,..,val,nb)`
//!   and `H1Basis_dFoo` becomes `d_foo(..,grad,nb)`;
//! * Fortran `fval(j, nbasis+i)` (1-based second index) is `val[*nb + i - 1]`
//!   and `grad(j, idx, c)` is `grad[idx-1][c-1]`;
//! * `nbasis` is the running count of functions already stored (`nb`), advanced
//!   exactly as the Fortran does. Value and gradient routines keep separate
//!   counters just like `nbp`/`nbdxp` in `ElementInfoVec`;
//! * local node numbers inside `edgedir` / `facedir` / `localnumbers` are kept
//!   **1-based** as in the Fortran (a `0` in the 4th face slot marks a triangle);
//! * several edge/face routines read the nodal values (or gradients) that the
//!   caller stored in `val[0..n_nodes]` / `grad[0..n_nodes]` first, exactly as
//!   the Fortran reads `fval(k,node)`. Call the nodal routine first.
//!
//! Operation order of every floating-point expression is kept as written in
//! the Fortran so results agree bit for bit with the pristine binary.

/// `H1Basis_MaxPElementEdges`.
pub const MAX_EDGES: usize = 12;
/// `H1Basis_MaxPElementFaces`.
pub const MAX_FACES: usize = 6;

/// `edgedir(2, 12)`: for each edge the two local (1-based) node numbers in the
/// order the hierarchic edge function is built from (`a` then `b`).
pub type EdgeDir = [[usize; 2]; MAX_EDGES];
/// `facedir(4, 6)`: for each face the local (1-based) nodes in function order;
/// `0` in slot 4 marks a triangular face.
pub type FaceDir = [[usize; 4]; MAX_FACES];

/// Fortran `x**j` with an integer exponent (`H1Basis_PowInt`).
///
/// gfortran lowers `x**j` for a run-time integer `j` to `__builtin_powi`
/// (repeated squaring); `f64::powi` is the same operation.
#[inline]
pub fn pow_int(x: f64, j: i32) -> f64 {
    x.powi(j)
}
