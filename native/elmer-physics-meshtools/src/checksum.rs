//! `MeshChecksum` (MeshChecksum.F90): eight numbers that identify a mesh
//! (counts, element-code, tag and node-index sums, coordinate sum) and a
//! pseudo-norm against `Reference Values`, so a test can check that a mesh
//! file converted by ElmerGrid is the mesh it expects.

use elmer_core::{ElemId, MeshView, ParamMap, Params, Value};
use elmer_solver::{Model, Solver, SolverCtx, SolverError, Variable};

const CALLER: &str = "MeshChecksum";

/// The variable `MeshChecksum_init` adds to the solver section, which holds
/// the pseudo-norm (`-nooutput -global meshcheck_var`).
const VARIABLE: &str = "meshcheck_var";

pub fn factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(MeshChecksum))
}

/// `CheckSum(1:8)`. Slot 1 is `2 * bulk elements + nodes`; 2-4 sum the bulk
/// elements' type code, body id and 1-based node indexes; 5-7 the same for
/// the boundary elements, with the boundary condition id as the tag (0 when
/// no `Boundary Condition` targets it); 8 is the sum of all coordinates.
pub fn checksum(mesh: &dyn MeshView, model: &Model) -> [f64; 8] {
    let nb = mesh.n_bulk();
    let mut cs = [0.0f64; 8];
    cs[0] = (nb + nb + mesh.n_nodes()) as f64;
    for t in 0..nb + mesh.n_boundary() {
        let el = mesh.element(ElemId(t as u32));
        let (ind, tag) = if t < nb {
            (1, model.body_of_tag(el.tag).map_or(el.tag, |b| b as i32))
        } else {
            (4, model.bc_of_tag(el.tag).map_or(0, |b| b as i32))
        };
        cs[ind] += f64::from(el.etype.0);
        cs[ind + 1] += f64::from(tag);
        cs[ind + 2] += el.nodes.iter().map(|&n| f64::from(n + 1)).sum::<f64>();
    }
    let n = mesh.n_nodes() as u32;
    let sum = |k: usize| (0..n).fold(0.0, |s, i| s + mesh.coord(i)[k]);
    cs[7] = sum(0) + sum(1) + sum(2);
    cs
}

/// The pseudo-norm of `cs` against `reference`: the mean over the eight
/// slots of the larger of `cs/ref` and its inverse (`1 + cs` for a zero
/// reference), so a matching mesh gives exactly one.
pub fn pseudo_norm(cs: &[f64; 8], reference: &[f64; 8]) -> f64 {
    let mut norm = 0.0;
    for i in 0..8 {
        let c = if reference[i].abs() > f64::EPSILON {
            let c = cs[i] / reference[i];
            c.max(1.0 / c)
        } else {
            1.0 + cs[i]
        };
        norm += c / 8.0;
    }
    norm
}

/// `Reference Values(8)`: the first column of the array.
fn reference_values(params: &ParamMap) -> Option<[f64; 8]> {
    let Some(Value::Real { rows, cols, data }) = params.get("reference values") else {
        return None;
    };
    (*rows >= 8).then(|| std::array::from_fn(|i| data[i * cols]))
}

pub struct MeshChecksum;

impl MeshChecksum {
    fn slot_index(ctx: &SolverCtx<'_>) -> usize {
        let id = ctx.solver_id();
        ctx.sim
            .slots
            .iter()
            .position(|s| s.id == id)
            .expect("running solver has a slot")
    }
}

impl Solver for MeshChecksum {
    /// `MeshChecksum_init`: a place for the pseudo-norm.
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let ix = Self::slot_index(ctx);
        if ctx.sim.slots[ix].variable.is_none() {
            let mut v = Variable::new(VARIABLE, elmer_core::FieldKind::Global, 1, None, vec![0.0]);
            v.output = false;
            let id = ctx.sim.vars.add(v)?;
            ctx.sim.slots[ix].variable = Some(id);
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "Checking for mesh consistency", 4);
        let cs = checksum(ctx.mesh, &ctx.sim.model);
        println!(" Checksums for file output:");
        println!(" ThisResults:{}", format_sums(&cs));
        let Some(reference) = reference_values(ctx.params()) else {
            return Ok(());
        };
        println!(" RefResults:{}", format_sums(&reference));
        let norm = pseudo_norm(&cs, &reference);
        println!(" PseudoNorm:{norm:>24.16}");
        let ix = Self::slot_index(ctx);
        if let Some(id) = ctx.sim.slots[ix].variable {
            ctx.sim
                .vars
                .by_id_mut(id)
                .values
                .iter_mut()
                .for_each(|v| *v = norm);
        }
        // By construction the reference norm should now be one.
        if ctx.sim.slots[ix].params.get("reference norm").is_none() {
            ctx.sim.slots[ix]
                .params
                .insert("Reference Norm", Value::real(1.0));
        }
        Ok(())
    }
}

/// `NINT(sum(1:7)), sum(8)` as list-directed output prints them.
fn format_sums(s: &[f64; 8]) -> String {
    let mut out = String::new();
    for v in &s[..7] {
        out.push_str(&format!("{:>12}", v.round() as i64));
    }
    out.push_str(&format!("{:>24.16}", s[7]));
    out
}
