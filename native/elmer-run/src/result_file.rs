//! `Output File` with `Binary Output = True`: Elmer's restart `.result`
//! file and its `.pos` index (`SaveResult`, ModelDescription.F90:3817),
//! written at each save point through [`elmer_solver::SaveCurrent`].
//!
//! Layout: an ASCII preamble (`BINARY 3.<E>`, date, the `Degrees of
//! freedom:` table), one NUL byte, then per save a `Time:` block and, for
//! each saved scalar variable, its NUL-terminated name, its Perm (size,
//! nonzero count and 1-based `(i, Perm(i))` pairs, or `-1` and the offset
//! of an identical previous Perm, or `0` for none) and its values in Perm
//! order, native-endian.
//!
//! Only the plain variant is written: one mesh, serial, every step saved.
//! Other settings leave the file unwritten with a warning, as before.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

use elmer_core::{normalize_key, FieldKind, Logger, MeshView, ParamMap, Params, Value, NO_DOF};
use elmer_solver::variable::VarHandle;
use elmer_solver::{Model, SaveCurrent, Simulation, SolverError};

const CALLER: &str = "SaveResult";

/// Why `Output File` is left unwritten, if it is.
fn unsupported(model: &Model, ranks: usize) -> Option<String> {
    let sim = &model.simulation;
    let flag = |k: &str| sim.logical(k) == Some(true);
    if ranks > 1 {
        return Some("in parallel".into());
    }
    if sim.logical("binary output") != Some(true) {
        return Some("ASCII output (Binary Output = False)".into());
    }
    for k in [
        "output file final only",
        "separate results files",
        "output file each timestep",
        "output coordinates",
        "output ip variables",
        "omit unchanged variables in output",
    ] {
        if flag(k) {
            return Some(k.into());
        }
    }
    if sim.integer("output file cycle").is_some_and(|c| c > 0) {
        return Some("output file cycle".into());
    }
    if sim.get("output variable 1").is_some() {
        return Some("output variable list".into());
    }
    // A step left unsaved gets ElmerSolver's final save (with its own
    // before/after saving solvers), which the time loop does not make.
    if let Some(Value::Integer(v)) = sim.get("output intervals") {
        if v.iter().any(|&i| i != 1) {
            return Some("Output Intervals other than 1".into());
        }
    }
    for p in &model.solvers {
        if p.get("output file").is_some() || p.get("mesh").is_some() {
            return Some("solver-specific meshes or output files".into());
        }
        if p.logical("eigen analysis") == Some(true) || p.logical("harmonic analysis") == Some(true) {
            return Some("eigen or harmonic analysis".into());
        }
    }
    None
}

/// The `Output File` writer of a case, or `None` when the case names none
/// or one this writer does not produce (with a warning saying why).
pub(crate) fn writer(model: &Model, mesh_dir: &std::path::Path, ranks: usize, log: &dyn Logger) -> Option<ResultWriter> {
    let name = model.simulation.string("output file")?;
    if let Some(why) = unsupported(model, ranks) {
        log.warn(CALLER, &format!("Output File {name} is not written: {why}"));
        return None;
    }
    Some(ResultWriter {
        path: mesh_dir.join(name),
        save_global: model.simulation.logical("output global variables") == Some(true),
        saves_done: 0,
    })
}

pub(crate) struct ResultWriter {
    path: PathBuf,
    save_global: bool,
    saves_done: i32,
}

/// One entry of Elmer's variable list as `SaveResult` sees it.
struct Entry<'a> {
    name: String,
    dofs: usize,
    size: usize,
    perm: Option<&'a Arc<[u32]>>,
    /// Values of this entry: `values[offset + stride * k]`.
    values: &'a [f64],
    offset: usize,
    stride: usize,
    kind: FieldKind,
    output: bool,
    equation: String,
}

fn entries<'a>(sim: &'a Simulation) -> Vec<Entry<'a>> {
    let equation = |id| {
        sim.slots
            .iter()
            .find(|s| s.variable == Some(id))
            .and_then(|s| sim.model.solvers.get(s.id - 1))
            .and_then(|p: &ParamMap| p.string("equation"))
            .map_or_else(|| "no equation".to_string(), normalize_key)
    };
    sim.vars
        .names_in_order()
        .into_iter()
        .map(|(name, h)| {
            let (id, offset, stride, dofs) = match h {
                VarHandle::Primary(id) => (id, 0, 1, sim.vars.by_id(id).dofs),
                VarHandle::Component { parent, component } => (parent, component, sim.vars.by_id(parent).dofs, 1),
            };
            let v = sim.vars.by_id(id);
            Entry {
                name,
                dofs,
                size: v.values.len() / stride,
                perm: v.perm.as_ref(),
                values: &v.values,
                offset,
                stride,
                kind: v.kind,
                output: v.output,
                equation: equation(id),
            }
        })
        .collect()
}

impl ResultWriter {
    /// `SaveThis` for a scalar entry of the data part.
    fn saves(&self, e: &Entry<'_>) -> bool {
        if !e.output || e.dofs != 1 || e.kind == FieldKind::Ip {
            return false;
        }
        if e.size == e.dofs {
            self.save_global
        } else {
            !e.name.starts_with("coordinate")
        }
    }

    /// `InitializeFile` and the header part; returns the `.pos` header.
    fn header(&self, list: &[Entry<'_>], n_nodes: usize) -> (Vec<u8>, Vec<u8>) {
        let marker = if cfg!(target_endian = "little") { 'L' } else { 'B' };
        let mut out = format!(" BINARY 3.{marker}\n");
        out.push_str(&format!("!File started at: {}\n", elmer_cli::version::format_date(std::time::SystemTime::now())));
        out.push_str(" Degrees of freedom: \n");
        let mut dofs = 0i32;
        let mut pos_names = Vec::new();
        for vector in [true, false] {
            for e in list {
                if !e.output || e.kind == FieldKind::Ip || (e.dofs > 1) != vector {
                    continue;
                }
                let save = if e.size == e.dofs { self.save_global } else { !e.name.starts_with("coordinate") };
                if !save {
                    continue;
                }
                let perm_size = e.perm.map_or(0, |p| p.len());
                out.push_str(&format!(
                    "{:<24} : {:>8}{:>8}{:>4} : {}\n",
                    e.name, e.size, perm_size, e.dofs, e.equation
                ));
                if !vector {
                    dofs += 1;
                    pos_names.extend_from_slice(e.name.as_bytes());
                    pos_names.push(0);
                }
            }
        }
        out.push_str(&format!(" Total DOFs: {dofs:>12}\n"));
        out.push_str(&format!(" Number Of Nodes: {n_nodes:>12}\n"));
        let mut bytes = out.into_bytes();
        bytes.push(0);
        let mut pos = vec![marker as u8];
        pos.extend_from_slice(&dofs.to_ne_bytes());
        pos.extend_from_slice(&pos_names);
        (bytes, pos)
    }
}

impl SaveCurrent for ResultWriter {
    fn save(&mut self, sim: &Simulation, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        let io = |e: std::io::Error| SolverError::fatal(CALLER, format!("{}: {e}", self.path.display()));
        log.info(CALLER, "-----------------------------------------", 5);
        log.info(CALLER, &format!("Saving results to file: {}", self.path.display()), 4);
        let list = entries(sim);
        let pos_path = PathBuf::from(format!("{}.pos", self.path.display()));
        let (mut data, mut pos, base) = if self.saves_done == 0 {
            let (d, p) = self.header(&list, mesh.n_nodes());
            (d, p, 0usize)
        } else {
            let len = std::fs::metadata(&self.path).map_err(io)?.len() as usize;
            (Vec::new(), Vec::new(), len)
        };
        let at = |data: &Vec<u8>| (base + data.len()) as i64;
        pos.extend_from_slice(&at(&data).to_ne_bytes());
        data.extend_from_slice(b"Time:\0");
        data.extend_from_slice(&(self.saves_done + 1).to_ne_bytes());
        let step = sim.vars.scalar("timestep").unwrap_or(0.0).round() as i32;
        data.extend_from_slice(&step.to_ne_bytes());
        data.extend_from_slice(&sim.vars.scalar("time").unwrap_or(0.0).to_ne_bytes());
        let mut prev: Option<(&Arc<[u32]>, i64)> = None;
        for e in list.iter().filter(|e| self.saves(e)) {
            pos.extend_from_slice(&at(&data).to_ne_bytes());
            data.extend_from_slice(e.name.as_bytes());
            data.push(0);
            match e.perm {
                None => data.extend_from_slice(&0i32.to_ne_bytes()),
                Some(p) if prev.is_some_and(|(q, _)| Arc::ptr_eq(p, q) || **p == **q) => {
                    data.extend_from_slice(&(-1i32).to_ne_bytes());
                    data.extend_from_slice(&prev.unwrap().1.to_ne_bytes());
                }
                Some(p) => {
                    prev = Some((p, at(&data)));
                    let nnz = p.iter().filter(|&&k| k != NO_DOF).count();
                    data.extend_from_slice(&(p.len() as i32).to_ne_bytes());
                    data.extend_from_slice(&(nnz as i32).to_ne_bytes());
                    for (i, &k) in p.iter().enumerate().filter(|(_, &k)| k != NO_DOF) {
                        data.extend_from_slice(&(i as i32 + 1).to_ne_bytes());
                        data.extend_from_slice(&(k as i32 + 1).to_ne_bytes());
                    }
                }
            }
            let value = |k: usize| e.values[e.offset + e.stride * k];
            match e.perm {
                Some(p) => {
                    for &k in p.iter().filter(|&&k| k != NO_DOF) {
                        data.extend_from_slice(&value(k as usize).to_ne_bytes());
                    }
                }
                None => {
                    for k in 0..e.size {
                        data.extend_from_slice(&value(k).to_ne_bytes());
                    }
                }
            }
        }
        let open = |p: &PathBuf, fresh: bool| {
            if fresh {
                std::fs::File::create(p)
            } else {
                OpenOptions::new().append(true).open(p)
            }
        };
        let fresh = self.saves_done == 0;
        open(&self.path, fresh).and_then(|mut f| f.write_all(&data)).map_err(io)?;
        open(&pos_path, fresh).and_then(|mut f| f.write_all(&pos)).map_err(io)?;
        self.saves_done += 1;
        log.info(CALLER, "Done writing results file", 5);
        log.info(CALLER, "-----------------------------------------", 5);
        Ok(())
    }
}
