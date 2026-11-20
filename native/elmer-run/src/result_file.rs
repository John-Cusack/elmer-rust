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
//! With `Binary Output = False` the same records are text (`ASCII 3`, no
//! `.pos`): `Time:` with `(a,i7,i7,ES17.8E3)`, `Perm: NULL`, `Perm: use
//! previous` or `Perm:` with `(A,i12," ",i12)` and `(2i11)` pairs, and one
//! list-directed value per line.
//!
//! Only the plain variant is written: one mesh. In parallel each
//! rank writes its own piece, `<Output File>.<rank>` (ModelDescription.F90:3894).
//! Other settings leave the file unwritten with a warning, as before.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

use elmer_core::{normalize_key, FieldKind, Logger, MeshView, Params, NO_DOF};
use elmer_solver::variable::{VarHandle, VarId};
use elmer_solver::{Model, SaveCurrent, Simulation, SolverError};

const CALLER: &str = "SaveResult";

/// Why `Output File` is left unwritten, if it is.
fn unsupported(model: &Model) -> Option<String> {
    let sim = &model.simulation;
    let flag = |k: &str| sim.logical(k) == Some(true);
    for k in [
        "output file final only",
        "separate results files",
        "output file each timestep",
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
    for p in &model.solvers {
        if p.get("output file").is_some() {
            return Some("solver-specific output files".into());
        }
        if p.logical("harmonic analysis") == Some(true) {
            return Some("harmonic analysis".into());
        }
        if p.logical("eigen analysis") == Some(true)
            && (p.logical("linear system complex") == Some(true)
                || sim.string("simulation type").is_some_and(|t| !normalize_key(t).starts_with("steady")))
        {
            return Some("eigen analysis of a complex system or over time steps".into());
        }
    }
    None
}

/// The `Output File` writer of rank `rank` of `ranks`, or `None` when the
/// case names none or one this writer does not produce (with a warning
/// saying why). In parallel the file name gets `.<rank>`, as `SaveResult`
/// does for `ParEnv % PEs > 1`.
pub fn writer(model: &Model, mesh_dir: &std::path::Path, rank: usize, ranks: usize, log: &dyn Logger) -> Option<ResultWriter> {
    let name = model.simulation.string("output file")?;
    if let Some(why) = unsupported(model) {
        log.warn(CALLER, &format!("Output File {name} is not written: {why}"));
        return None;
    }
    let name = if ranks > 1 { format!("{name}.{rank}") } else { name.to_string() };
    Some(ResultWriter {
        path: mesh_dir.join(name),
        save_global: model.simulation.logical("output global variables") == Some(true),
        binary: model.simulation.logical("binary output") == Some(true),
        save_coordinates: model.simulation.logical("output coordinates") == Some(true),
        rank,
        saves_done: 0,
    })
}

/// [`SaveCurrent`] writing one rank's `Output File`.
pub struct ResultWriter {
    path: PathBuf,
    save_global: bool,
    /// `Binary Output`; text otherwise.
    binary: bool,
    /// `Output Coordinates`: `coordinate 1..3`, the mesh's, saved first.
    save_coordinates: bool,
    /// The `Partition` global Elmer always has (`ParEnv % MyPe`), saved
    /// with the globals when the case has no such variable.
    rank: usize,
    saves_done: i32,
}

/// One entry of Elmer's variable list as `SaveResult` sees it.
struct Entry<'a> {
    name: String,
    /// `Var % Name` in the header: a vector declared with components
    /// carries them, `flow solution[velocity:2 pressure:1]`.
    header_name: String,
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

/// A field the variable store does not hold, by name.
type Named = (String, Vec<f64>);

/// An entry of values the variable store does not hold: Perm NULL, no
/// equation.
fn plain<'a>(name: &str, values: &'a [f64], kind: FieldKind) -> Entry<'a> {
    Entry {
        name: name.to_string(),
        header_name: name.to_string(),
        dofs: 1,
        size: values.len(),
        perm: None,
        values,
        offset: 0,
        stride: 1,
        kind,
        output: true,
        equation: "no equation".to_string(),
    }
}

/// The variable list, with the values of `over` (a primary variable and
/// what to save for it) in place of its own; its components follow.
/// `first` and `last` (coordinates, `Partition`) go before and after it,
/// where Elmer's mesh load and `AddTimeEtc` put them.
fn entries<'a>(
    sim: &'a Simulation,
    over: Option<(VarId, &'a [f64])>,
    first: &'a [Named],
    last: &'a [Named],
) -> Vec<Entry<'a>> {
    let mut out: Vec<Entry<'a>> = first.iter().map(|(n, v)| plain(n, v, FieldKind::Nodal)).collect();
    out.extend(variable_entries(sim, over));
    out.extend(last.iter().map(|(n, v)| plain(n, v, FieldKind::Global)));
    out
}

fn variable_entries<'a>(sim: &'a Simulation, over: Option<(VarId, &'a [f64])>) -> Vec<Entry<'a>> {
    // `Var % Solver`'s equation, on the mesh being saved.
    let equation = |id| sim.equation_of(id).unwrap_or_else(|| "no equation".to_string());
    sim.vars
        .names_in_order()
        .into_iter()
        .map(|(name, h)| {
            let (id, offset, stride, dofs) = match h {
                VarHandle::Primary(id) => (id, 0, 1, sim.vars.by_id(id).dofs),
                VarHandle::Component { parent, component } => (parent, component, sim.vars.by_id(parent).dofs, 1),
            };
            let v = sim.vars.by_id(id);
            // Declared components (`Flow[Velocity:2 Pressure:1]`) show; the
            // default `<name> k` ones may or may not have been declared, and
            // stay out.
            let comps = sim.vars.component_names(id);
            let declared = comps.iter().enumerate().any(|(k, c)| *c != format!("{name} {}", k + 1));
            let header_name = match h {
                VarHandle::Primary(_) if dofs > 1 && declared => format!("{name}[{}]", component_spec(&comps)),
                _ => name.clone(),
            };
            Entry {
                name,
                header_name,
                dofs,
                size: v.values.len() / stride,
                perm: v.perm.as_ref(),
                values: match over {
                    Some((o, vals)) if o == id => vals,
                    _ => &v.values,
                },
                offset,
                stride,
                kind: v.kind,
                output: v.output,
                equation: equation(id),
            }
        })
        .collect()
}

/// The `[name:dofs ...]` of a vector variable from its per-dof component
/// names: consecutive `<base> <k>` names make one `<base>:<count>`.
fn component_spec(names: &[String]) -> String {
    let base = |n: &str| match n.rsplit_once(' ') {
        Some((b, k)) if !k.is_empty() && k.bytes().all(|c| c.is_ascii_digit()) => b.to_string(),
        _ => n.to_string(),
    };
    let mut groups: Vec<(String, usize)> = Vec::new();
    for n in names {
        let b = base(n);
        match groups.last_mut() {
            Some((g, c)) if *g == b && b != *n => *c += 1,
            _ => groups.push((b, 1)),
        }
    }
    groups.iter().map(|(b, c)| format!("{b}:{c}")).collect::<Vec<_>>().join(" ")
}

impl ResultWriter {
    /// `SaveThis` for a scalar entry of the data part.
    fn saves(&self, e: &Entry<'_>) -> bool {
        e.dofs == 1 && self.saves_header(e)
    }

    /// `SaveThis` of the header part, vectors included.
    fn saves_header(&self, e: &Entry<'_>) -> bool {
        if !e.output || e.kind == FieldKind::Ip {
            return false;
        }
        if e.size == e.dofs {
            self.save_global
        } else {
            self.save_coordinates || !e.name.starts_with("coordinate")
        }
    }

    /// The entries Elmer has that the variable store does not: the mesh's
    /// coordinates (when saved) and `Partition`.
    fn extras(&self, sim: &Simulation, mesh: &dyn MeshView) -> (Vec<Named>, Vec<Named>) {
        let mut first = Vec::new();
        if self.save_coordinates && sim.vars.resolve("coordinate 1").is_none() {
            for a in 0..3 {
                let v = (0..mesh.n_nodes() as u32).map(|k| mesh.coord(k)[a]).collect();
                first.push((format!("coordinate {}", a + 1), v));
            }
        }
        let mut last = Vec::new();
        if self.save_global && sim.vars.resolve("partition").is_none() {
            last.push(("partition".to_string(), vec![self.rank as f64]));
        }
        (first, last)
    }

    /// `InitializeFile` and the header part; returns the `.pos` header.
    fn header(&self, list: &[Entry<'_>], n_nodes: usize) -> (Vec<u8>, Vec<u8>) {
        let marker = if cfg!(target_endian = "little") { 'L' } else { 'B' };
        let mut out = if self.binary { format!(" BINARY 3.{marker}\n") } else { " ASCII 3\n".to_string() };
        out.push_str(&format!("!File started at: {}\n", elmer_cli::version::format_date(std::time::SystemTime::now())));
        out.push_str(" Degrees of freedom: \n");
        let mut dofs = 0i32;
        let mut pos_names = Vec::new();
        for vector in [true, false] {
            for e in list {
                if (e.dofs > 1) != vector || !self.saves_header(e) {
                    continue;
                }
                let perm_size = e.perm.map_or(0, |p| p.len());
                out.push_str(&format!(
                    "{:<24} : {:>8}{:>8}{:>4} : {}\n",
                    e.header_name, e.size, perm_size, e.dofs, e.equation
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
        if self.binary {
            bytes.push(0);
        }
        let mut pos = vec![marker as u8];
        pos.extend_from_slice(&dofs.to_ne_bytes());
        pos.extend_from_slice(&pos_names);
        (bytes, pos)
    }
}

impl SaveCurrent for ResultWriter {
    /// `SaveCurrent` (ElmerSolver.F90:3701) for this mesh: each solver with
    /// an eigen analysis saves its modes, mode `j` as step `j`, and then
    /// its variable holds zeros (`Var % Values = 0`); the usual save
    /// follows unless the last solver is an eigen one (the Fortran's
    /// `EigAnal` is whatever the last solver left).
    fn save(&mut self, sim: &Simulation, mesh: &dyn MeshView, log: &dyn Logger) -> Result<(), SolverError> {
        let mut eig_anal = false;
        let mut zeroed: Option<(VarId, Vec<f64>)> = None;
        for slot in &sim.slots {
            eig_anal = slot.params.logical("eigen analysis") == Some(true);
            if !eig_anal {
                continue;
            }
            let Some(id) = slot.variable else { continue };
            let var = sim.vars.by_id(id);
            if var.eigen_vectors.is_empty() {
                continue;
            }
            for (j, mode) in var.eigen_vectors.iter().enumerate() {
                let re: Vec<f64> = mode.iter().map(|z| z.re).collect();
                self.save_one(sim, mesh, log, j as i32 + 1, Some((id, &re)))?;
            }
            zeroed = Some((id, vec![0.0; var.values.len()]));
        }
        if eig_anal {
            return Ok(());
        }
        let step = sim.vars.scalar("timestep").unwrap_or(0.0).round() as i32;
        let over = zeroed.as_ref().map(|(id, z)| (*id, &z[..]));
        self.save_one(sim, mesh, log, step, over)
    }
}

impl ResultWriter {
    /// `SaveResult` of one save point: `step` is the `Time:` block's step.
    fn save_one(
        &mut self,
        sim: &Simulation,
        mesh: &dyn MeshView,
        log: &dyn Logger,
        step: i32,
        over: Option<(VarId, &[f64])>,
    ) -> Result<(), SolverError> {
        let io = |e: std::io::Error| SolverError::fatal(CALLER, format!("{}: {e}", self.path.display()));
        log.info(CALLER, "-----------------------------------------", 5);
        log.info(CALLER, &format!("Saving results to file: {}", self.path.display()), 4);
        let (first, last) = self.extras(sim, mesh);
        let list = entries(sim, over, &first, &last);
        let pos_path = PathBuf::from(format!("{}.pos", self.path.display()));
        let (mut data, mut pos, base) = if self.saves_done == 0 {
            let (d, p) = self.header(&list, mesh.n_nodes());
            (d, p, 0usize)
        } else {
            let len = std::fs::metadata(&self.path).map_err(io)?.len() as usize;
            (Vec::new(), Vec::new(), len)
        };
        if self.binary {
            self.binary_step(sim, &list, step, base, &mut data, &mut pos);
        } else {
            self.text_step(sim, &list, step, &mut data);
        }
        let open = |p: &PathBuf, fresh: bool| {
            if fresh {
                std::fs::File::create(p)
            } else {
                OpenOptions::new().append(true).open(p)
            }
        };
        let fresh = self.saves_done == 0;
        if fresh {
            if let Some(d) = self.path.parent() {
                std::fs::create_dir_all(d).map_err(io)?;
            }
        }
        open(&self.path, fresh).and_then(|mut f| f.write_all(&data)).map_err(io)?;
        if self.binary {
            open(&pos_path, fresh).and_then(|mut f| f.write_all(&pos)).map_err(io)?;
        }
        self.saves_done += 1;
        log.info(CALLER, "Done writing results file", 5);
        log.info(CALLER, "-----------------------------------------", 5);
        Ok(())
    }

    /// One `Time:` block of the binary file, and its `.pos` offsets.
    fn binary_step(&self, sim: &Simulation, list: &[Entry<'_>], step: i32, base: usize, data: &mut Vec<u8>, pos: &mut Vec<u8>) {
        let at = |data: &Vec<u8>| (base + data.len()) as i64;
        pos.extend_from_slice(&at(data).to_ne_bytes());
        data.extend_from_slice(b"Time:\0");
        data.extend_from_slice(&(self.saves_done + 1).to_ne_bytes());
        data.extend_from_slice(&step.to_ne_bytes());
        data.extend_from_slice(&sim.vars.scalar("time").unwrap_or(0.0).to_ne_bytes());
        let mut prev: Option<(&Arc<[u32]>, i64)> = None;
        for e in list.iter().filter(|e| self.saves(e)) {
            pos.extend_from_slice(&at(data).to_ne_bytes());
            data.extend_from_slice(e.name.as_bytes());
            data.push(0);
            match e.perm {
                None => data.extend_from_slice(&0i32.to_ne_bytes()),
                Some(p) if prev.is_some_and(|(q, _)| Arc::ptr_eq(p, q) || **p == **q) => {
                    data.extend_from_slice(&(-1i32).to_ne_bytes());
                    data.extend_from_slice(&prev.unwrap().1.to_ne_bytes());
                }
                Some(p) => {
                    prev = Some((p, at(data)));
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
    }

    /// One `Time:` block of the text file (`Binary Output = False`).
    fn text_step(&self, sim: &Simulation, list: &[Entry<'_>], step: i32, data: &mut Vec<u8>) {
        use elmer_solver::elmer_saveio::fortran::{fortran_es_e, fortran_i, list_real};
        let time = sim.vars.scalar("time").unwrap_or(0.0);
        let mut t = format!(
            "Time: {}{}{}\n",
            fortran_i(i64::from(self.saves_done + 1), 7),
            fortran_i(i64::from(step), 7),
            fortran_es_e(time, 17, 8, Some(3))
        );
        let mut prev: Option<&Arc<[u32]>> = None;
        for e in list.iter().filter(|e| self.saves(e)) {
            t.push_str(&e.name);
            t.push('\n');
            match e.perm {
                None => t.push_str("Perm: NULL\n"),
                Some(p) if prev.is_some_and(|q| Arc::ptr_eq(p, q) || **p == **q) => t.push_str("Perm: use previous\n"),
                Some(p) => {
                    prev = Some(p);
                    let nnz = p.iter().filter(|&&k| k != NO_DOF).count();
                    t.push_str(&format!("Perm: {} {}\n", fortran_i(p.len() as i64, 12), fortran_i(nnz as i64, 12)));
                    for (i, &k) in p.iter().enumerate().filter(|(_, &k)| k != NO_DOF) {
                        t.push_str(&fortran_i(i as i64 + 1, 11));
                        t.push_str(&fortran_i(i64::from(k) + 1, 11));
                        t.push('\n');
                    }
                }
            }
            let mut put = |x: f64| {
                t.push_str(&list_real(x));
                t.push('\n');
            };
            let value = |k: usize| e.values[e.offset + e.stride * k];
            match e.perm {
                Some(p) => p.iter().filter(|&&k| k != NO_DOF).for_each(|&k| put(value(k as usize))),
                None => (0..e.size).for_each(|k| put(value(k))),
            }
        }
        data.extend_from_slice(t.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use elmer_core::{NullLogger, Value};
    use elmer_solver::Variable;

    /// MinimalExample's `triangle/case.result` from the pristine Release
    /// arm: the binary part byte for byte, and the preamble but for the
    /// date and the equation name (the variable here has no solver slot).
    #[test]
    fn vector_names_carry_their_components() {
        let s = |v: &[&str]| component_spec(&v.iter().map(|x| x.to_string()).collect::<Vec<_>>());
        assert_eq!(s(&["velocity 1", "velocity 2", "pressure"]), "velocity:2 pressure:1");
        assert_eq!(s(&["temperature flux 1", "temperature flux 2"]), "temperature flux:2");
        assert_eq!(s(&["temperature", "potential"]), "temperature:1 potential:1");
        assert_eq!(s(&["stress_xx", "stress_yy"]), "stress_xx:1 stress_yy:1");
    }

    #[test]
    fn minimal_example_matches_pristine() {
        let dir = std::env::temp_dir().join(format!("elmer-run-result-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut model = Model::default();
        model.simulation.insert("Output File", Value::Str("case.result".into()));
        model.simulation.insert("Binary Output", Value::Logical(true));
        let mut sim = Simulation::new(model);
        let perm: Arc<[u32]> = Arc::from(vec![2u32, 1, 0]);
        sim.vars.add(Variable::new("Temperature", FieldKind::Nodal, 1, Some(perm), vec![0.0, 0.0, 0.5000000000000001])).unwrap();
        sim.vars.set_scalar("timestep", 1.0);
        sim.vars.set_scalar("time", 1.0);
        let mesh = elmer_mesh_levels::Mesh { dim: 2, coords: vec![[0.0; 3]; 3], elements: Vec::new(), n_bulk: 0 };
        let mut w = writer(&sim.model, &dir, 0, 1, &NullLogger).unwrap();
        w.save(&sim, &mesh, &NullLogger).unwrap();

        let got = std::fs::read(dir.join("case.result")).unwrap();
        let nul = got.iter().position(|&b| b == 0).unwrap();
        let preamble = String::from_utf8(got[..nul].to_vec()).unwrap();
        let lines: Vec<&str> = preamble.lines().filter(|l| !l.starts_with("!File started at: ")).collect();
        assert_eq!(
            lines,
            [
                " BINARY 3.L",
                " Degrees of freedom: ",
                "temperature              :        3       3   1 : no equation",
                " Total DOFs:            1",
                " Number Of Nodes:            3",
            ]
        );
        let mut want = b"Time:\0".to_vec();
        for i in [1i32, 1] {
            want.extend_from_slice(&i.to_le_bytes());
        }
        want.extend_from_slice(&1.0f64.to_le_bytes());
        want.extend_from_slice(b"temperature\0");
        for i in [3i32, 3, 1, 3, 2, 2, 3, 1] {
            want.extend_from_slice(&i.to_le_bytes());
        }
        for x in [0.5000000000000001f64, 0.0, 0.0] {
            want.extend_from_slice(&x.to_le_bytes());
        }
        assert_eq!(&got[nul + 1..], &want[..]);
        let pos = std::fs::read(dir.join("case.result.pos")).unwrap();
        let mut want = b"L".to_vec();
        want.extend_from_slice(&1i32.to_le_bytes());
        want.extend_from_slice(b"temperature\0");
        want.extend_from_slice(&((nul + 1) as i64).to_le_bytes());
        want.extend_from_slice(&((nul + 1 + 22) as i64).to_le_bytes());
        assert_eq!(pos, want);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// `Binary Output = False`: the text layout of `reload`'s
    /// `Step/ad-trans.dat` (pristine Release), with no `.pos`, read back by
    /// `elmer-saveio`.
    #[test]
    fn ascii_output_matches_pristine_layout() {
        let dir = std::env::temp_dir().join(format!("elmer-run-ascii-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut model = Model::default();
        model.simulation.insert("Output File", Value::Str("case.dat".into()));
        let mut sim = Simulation::new(model);
        let perm: Arc<[u32]> = Arc::from(vec![1u32, 0]);
        sim.vars.add(Variable::new("Species", FieldKind::Nodal, 1, Some(perm.clone()), vec![0.1, 1.0])).unwrap();
        sim.vars.add(Variable::new("Zeta", FieldKind::Nodal, 1, Some(perm), vec![0.0, -2.5e-7])).unwrap();
        sim.vars.set_scalar("timestep", 5.0);
        sim.vars.set_scalar("time", 12.5);
        let mesh = elmer_mesh_levels::Mesh { dim: 2, coords: vec![[0.0; 3]; 2], elements: Vec::new(), n_bulk: 0 };
        let mut w = writer(&sim.model, &dir, 0, 1, &NullLogger).unwrap();
        w.save(&sim, &mesh, &NullLogger).unwrap();
        sim.vars.set_scalar("timestep", 6.0);
        w.save(&sim, &mesh, &NullLogger).unwrap();

        let got = std::fs::read_to_string(dir.join("case.dat")).unwrap();
        let lines: Vec<&str> = got.lines().filter(|l| !l.starts_with("!File started at: ")).collect();
        let step = |n: &'static str| {
            [
                n,
                "species",
                "Perm:            2            2",
                "          1          2",
                "          2          1",
                "   1.0000000000000000     ",
                "  0.10000000000000001     ",
                "zeta",
                "Perm: use previous",
                "  -2.4999999999999999E-007",
                "   0.0000000000000000     ",
            ]
        };
        let mut want = vec![
            " ASCII 3",
            " Degrees of freedom: ",
            "species                  :        2       2   1 : no equation",
            "zeta                     :        2       2   1 : no equation",
            " Total DOFs:            2",
            " Number Of Nodes:            2",
        ];
        want.extend(step("Time:       1      5  1.25000000E+001"));
        want.extend(step("Time:       2      6  1.25000000E+001"));
        assert_eq!(lines, want);
        assert!(!dir.join("case.dat.pos").exists());
        let rf = elmer_solver::elmer_saveio::result::read_result(&dir.join("case.dat")).unwrap();
        assert_eq!(rf.steps.len(), 2);
        let f = rf.restart_field(1, "zeta").unwrap().unwrap();
        assert_eq!((f.node, f.slot, f.values), (vec![0, 1], vec![1, 0], vec![-2.5e-7, 0.0]));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
