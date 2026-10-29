//! Serial element loop: the `DO t=1,Active; GetActiveElement; ...;
//! DefaultUpdateEquations` pattern of Elmer solvers (DefUtils.F90), with
//! the kernel supplied by the physics crate.
//!
//! The loop is deterministic: elements are visited in the order given and
//! glued one at a time, so the global sums match Elmer's for the same
//! element order. There is no threaded variant.

use elmer_core::{ElemId, ElementRef, GlobalAssembler, LocalSystem, MeshView, NO_DOF};

use crate::perm::Perm;

/// Bulk elements in Elmer order (`Mesh % Elements(1:NumberOfBulkElements)`).
pub fn bulk_elements(mesh: &dyn MeshView) -> impl Iterator<Item = ElemId> + Clone {
    (0..mesh.n_bulk() as u32).map(ElemId)
}

/// Boundary elements in Elmer order, after the bulk ones.
pub fn boundary_elements(mesh: &dyn MeshView) -> impl Iterator<Item = ElemId> + Clone {
    let nb = mesh.n_bulk() as u32;
    (nb..nb + mesh.n_boundary() as u32).map(ElemId)
}

/// Runs `kernel` on each element of `elements` and glues the local system
/// it fills.
///
/// For each element the local system is zeroed and sized
/// `nodes * dofs`; `dof_index` is the element's nodes mapped through
/// `perm` (Elmer `Indexes = Perm(NodeIndexes)`, absent nodes as
/// [`NO_DOF`]). The kernel returns `Ok(false)` to skip gluing (Elmer's
/// `CYCLE` before `DefaultUpdateEquations`, e.g. a boundary element
/// without the BC this loop handles). Elements whose nodes are all absent
/// are skipped without calling the kernel. The first kernel error stops
/// the loop.
///
/// Mass and damping are glued only if the kernel sets `local.mass` /
/// `local.damp` and the assembler accepts them; time integration of the
/// local system belongs in the kernel's caller (L07), before gluing.
pub fn assemble_elements<E>(
    mesh: &dyn MeshView,
    perm: &Perm,
    dofs: usize,
    elements: impl IntoIterator<Item = ElemId>,
    kernel: &mut dyn FnMut(ElemId, ElementRef<'_>, &mut LocalSystem) -> Result<bool, E>,
    assembler: &mut dyn GlobalAssembler,
) -> Result<(), E> {
    let mut local = LocalSystem::zeros(0);
    let mut dof_index: Vec<u32> = Vec::new();
    for e in elements {
        let el = mesh.element(e);
        dof_index.clear();
        dof_index.extend(el.nodes.iter().map(|&n| perm.perm[n as usize]));
        if dof_index.iter().all(|&d| d == NO_DOF) {
            continue;
        }
        reset(&mut local, el.nodes.len() * dofs);
        if kernel(e, el, &mut local)? {
            assembler.add_local(&dof_index, &local);
        }
    }
    Ok(())
}

fn reset(local: &mut LocalSystem, n: usize) {
    local.n = n;
    local.stiff.clear();
    local.stiff.resize(n * n, 0.0);
    local.force.clear();
    local.force.resize(n, 0.0);
    local.mass = None;
    local.damp = None;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glue::CsrAssembler;
    use crate::perm::{create_system, PermOptions};
    use elmer_core::ElementType;

    struct Line;

    impl MeshView for Line {
        fn dim(&self) -> usize {
            1
        }
        fn n_nodes(&self) -> usize {
            3
        }
        fn coord(&self, node: u32) -> [f64; 3] {
            [node as f64, 0.0, 0.0]
        }
        fn n_bulk(&self) -> usize {
            2
        }
        fn n_boundary(&self) -> usize {
            1
        }
        fn element(&self, e: ElemId) -> ElementRef<'_> {
            const NODES: [[u32; 2]; 2] = [[0, 1], [1, 2]];
            match e.0 {
                0 | 1 => ElementRef {
                    etype: ElementType(202),
                    nodes: &NODES[e.0 as usize],
                    tag: 1,
                    parents: [None, None],
                },
                _ => ElementRef {
                    etype: ElementType(101),
                    nodes: &[2],
                    tag: 1,
                    parents: [Some(ElemId(1)), None],
                },
            }
        }
    }

    #[test]
    fn loop_glues_bulk_then_boundary_and_honours_skip() {
        let mesh = Line;
        let active = |e: ElemId| e.0 < 2;
        let mut sys = create_system(&mesh, &active, PermOptions::default(), 1);
        let mut b = vec![0.0; sys.matrix.nrows];
        let mut asm = CsrAssembler::new(&mut sys.matrix, &mut b);
        let mut seen = Vec::new();
        assemble_elements::<()>(
            &mesh,
            &sys.perm,
            1,
            bulk_elements(&mesh).chain(boundary_elements(&mesh)),
            &mut |e, el, local| {
                seen.push(e.0);
                assert_eq!(local.n, el.nodes.len());
                assert!(local.stiff.iter().all(|&v| v == 0.0));
                if el.etype.family() == 2 {
                    local.stiff.copy_from_slice(&[1.0, -1.0, -1.0, 1.0]);
                    Ok(true)
                } else {
                    local.force[0] = 7.0;
                    Ok(e.0 == 2)
                }
            },
            &mut asm,
        )
        .unwrap();
        assert_eq!(seen, vec![0, 1, 2]);
        let p = |n: usize| sys.perm.perm[n] as usize;
        let d = sys.matrix.diag.as_ref().unwrap();
        assert_eq!(sys.matrix.values[d[p(1)]], 2.0);
        assert_eq!(sys.matrix.values[d[p(0)]], 1.0);
        assert_eq!(b[p(2)], 7.0);
    }

    #[test]
    fn kernel_error_stops_the_loop() {
        let mesh = Line;
        let mut sys = create_system(&mesh, &|e| e.0 < 2, PermOptions::default(), 1);
        let mut b = vec![0.0; sys.matrix.nrows];
        let mut asm = CsrAssembler::new(&mut sys.matrix, &mut b);
        let mut calls = 0;
        let r = assemble_elements(
            &mesh,
            &sys.perm,
            1,
            bulk_elements(&mesh),
            &mut |_, _, _| {
                calls += 1;
                Err("bad element")
            },
            &mut asm,
        );
        assert_eq!(r, Err("bad element"));
        assert_eq!(calls, 1);
    }
}
