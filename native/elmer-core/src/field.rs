//! Consumer view of solution fields, for writers, norms and restart.

/// Elmer Perm 0 ("not in system").
pub const NO_DOF: u32 = u32::MAX;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FieldKind {
    Nodal,
    Elemental,
    Dg,
    Ip,
    Global,
    Edge,
}

#[derive(Copy, Clone, Debug)]
pub struct FieldRef<'a> {
    pub name: &'a str,
    pub kind: FieldKind,
    pub dofs: usize,
    /// Node to dof slot, [`NO_DOF`] if absent.
    pub perm: Option<&'a [u32]>,
    pub values: &'a [f64],
}

pub trait FieldSource {
    fn fields(&self) -> Vec<FieldRef<'_>>;
    fn time(&self) -> f64;
    fn step(&self) -> usize;
}
