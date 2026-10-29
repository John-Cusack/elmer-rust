//! Element-local systems: physics -> time integration -> assembly.

#[derive(Clone, Debug)]
pub struct LocalSystem {
    /// `nd * dofs`.
    pub n: usize,
    /// Row-major `n * n`.
    pub stiff: Vec<f64>,
    pub mass: Option<Vec<f64>>,
    pub damp: Option<Vec<f64>>,
    pub force: Vec<f64>,
}

impl LocalSystem {
    /// Zero stiffness and force, no mass or damping.
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            stiff: vec![0.0; n * n],
            mass: None,
            damp: None,
            force: vec![0.0; n],
        }
    }
}

pub trait GlobalAssembler {
    fn add_local(&mut self, dof_index: &[u32], local: &LocalSystem);
}
