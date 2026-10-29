//! How a mesh mapper moves the nodes.
//!
//! Elmer's mappers assign `Mesh % Nodes % x/y/z` in place. In the Rust
//! framework a solver sees the mesh only as `&dyn MeshView` (`coord` returns
//! by value), so it cannot move a node; the runner owns the coordinates.
//! [`MeshMotion`] is the one thing a mapper needs from it. The runner
//! installs its implementation for the thread that runs the simulation
//! ([`install`]) and its `MeshView::coord` reads the shared coordinates
//! ([`SharedCoords`] is a ready-made store). Until `elmer-core` grows a
//! `MeshView::set_coords`, this is the seam; the mappers refuse with
//! `NotWired` when nothing is installed and a node would move.

use std::cell::RefCell;
use std::sync::{Arc, PoisonError, RwLock};

/// Replaces the coordinates of every node (`Mesh % Nodes`).
pub trait MeshMotion {
    fn set_coords(&self, coords: &[[f64; 3]]);
}

thread_local! {
    static MOTION: RefCell<Option<Arc<dyn MeshMotion>>> = const { RefCell::new(None) };
}

/// Installs `motion` for solvers run on this thread, returning the previous
/// one. `None` removes it.
pub fn install(motion: Option<Arc<dyn MeshMotion>>) -> Option<Arc<dyn MeshMotion>> {
    MOTION.with(|m| std::mem::replace(&mut *m.borrow_mut(), motion))
}

pub(crate) fn current() -> Option<Arc<dyn MeshMotion>> {
    MOTION.with(|m| m.borrow().clone())
}

/// Node coordinates behind a lock: a `MeshView::coord` reads them and the
/// mappers write them through [`MeshMotion`].
#[derive(Debug, Default)]
pub struct SharedCoords(RwLock<Vec<[f64; 3]>>);

impl SharedCoords {
    pub fn new(coords: Vec<[f64; 3]>) -> Self {
        Self(RwLock::new(coords))
    }

    pub fn get(&self, node: u32) -> [f64; 3] {
        self.0.read().unwrap_or_else(PoisonError::into_inner)[node as usize]
    }

    pub fn len(&self) -> usize {
        self.0.read().unwrap_or_else(PoisonError::into_inner).len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn snapshot(&self) -> Vec<[f64; 3]> {
        self.0
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .clone()
    }
}

impl MeshMotion for SharedCoords {
    fn set_coords(&self, coords: &[[f64; 3]]) {
        *self.0.write().unwrap_or_else(PoisonError::into_inner) = coords.to_vec();
    }
}
