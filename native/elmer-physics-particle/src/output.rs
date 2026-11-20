//! `ParticleOutputSolver` (SaveGridData.F90): `SaveParticleData` with
//! `ParticleOutputVtu` / `WriteVtuFile` (ParticleUtils.F90) writing the
//! particles of the run (`GlobalParticles`) as VTK XML vertex cells, raw
//! appended binary.

use std::cell::RefCell;
use std::io::Write;

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

use crate::faces::V3;
use crate::particles::{PARTICLE_INITIATED, PARTICLE_LOST};

pub const OUTPUT_PROCEDURE: (&str, &str) = ("SaveGridData", "ParticleOutputSolver");

/// The keyword elmer-run gives the section: the mesh directory, the
/// default `SolverOutputDirectory( ..., UseMeshDir = .TRUE. )`.
pub const MESH_DIRECTORY_KEY: &str = "Particle Output Mesh Directory";

const CALLER: &str = "ParticleOutputSolver";

/// The particle data the output reads: `GlobalParticles` after the last
/// ParticleDynamics visit of this thread.
#[derive(Clone, Default)]
pub(crate) struct Snapshot {
    pub dim: usize,
    pub coord: Vec<V3>,
    pub status: Vec<i32>,
}

thread_local! {
    static GLOBAL: RefCell<Option<Snapshot>> = const { RefCell::new(None) };
}

pub(crate) fn publish(s: Snapshot) {
    GLOBAL.with(|g| *g.borrow_mut() = Some(s));
}

fn snapshot() -> Option<Snapshot> {
    GLOBAL.with(|g| g.borrow().clone())
}

pub fn register_output(registry: &mut SolverRegistry) {
    registry.register(OUTPUT_PROCEDURE.0, OUTPUT_PROCEDURE.1, output_factory);
}

pub fn output_factory(_params: &ParamMap) -> Result<Box<dyn Solver>, SolverError> {
    Ok(Box::new(ParticleOutput::default()))
}

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("ParticleOutputSolver: {what}"))
}

/// The `SAVE`d `nTime` of `ParticleOutputVtu`.
#[derive(Default)]
pub struct ParticleOutput {
    n_time: i64,
}

impl Solver for ParticleOutput {
    fn init(&mut self, ctx: &mut SolverCtx<'_>) -> Result<(), SolverError> {
        let p = ctx.params();
        for k in [
            "table format",
            "gmsh format",
            "ascii output",
            "enforce parallel format",
        ] {
            if p.logical(k) == Some(true) {
                return Err(unsupported(k));
            }
        }
        if p.logical("binary output") == Some(false) {
            return Err(unsupported("Binary Output = False (ascii)"));
        }
        if let Some(f) = p.string("output format").map(elmer_core::normalize_key) {
            if f == "table" || f == "gmsh" {
                return Err(unsupported(&format!("Output Format = {f}")));
            }
        }
        for k in [
            "output directory",
            "particle save group",
            "scalar field 1",
            "vector field 1",
        ] {
            if p.get(k).is_some() {
                return Err(unsupported(k));
            }
        }
        if ctx.simulation_params().get("output directory").is_some() {
            return Err(unsupported("Simulation Output Directory"));
        }
        if p.string("filename prefix")
            .is_some_and(|f| f.contains('/') || f.contains('\\') || f.contains(':'))
        {
            return Err(unsupported("a Filename Prefix with a directory"));
        }
        if ctx.comm.size() > 1 {
            return Err(unsupported("more than one partition (pvtu)"));
        }
        Ok(())
    }

    fn solve(
        &mut self,
        ctx: &mut SolverCtx<'_>,
        _dt: f64,
        _transient: bool,
    ) -> Result<(), SolverError> {
        ctx.info(CALLER, "-----------------------------------------", 4);
        ctx.info(CALLER, "Saving particle data", 4);
        let p = ctx.params().clone();
        // SaveParticleData: the file formats.
        let fmt = p.string("output format").map(elmer_core::normalize_key);
        let vtu = p.logical("vtu format").unwrap_or(false) || fmt.as_deref() == Some("vtu");
        if !vtu {
            ctx.warn("SaveParticleData", "No active file format given!");
        } else {
            self.write_vtu(ctx, &p)?;
        }
        ctx.info(CALLER, "All done", 4);
        ctx.info(CALLER, "-----------------------------------------", 4);
        Ok(())
    }
}

impl ParticleOutput {
    /// `ParticleOutputVtu` and its `WriteVtuFile`, binary appended, no
    /// fields (`Save Fields` without `Scalar/Vector Field i`).
    fn write_vtu(&mut self, ctx: &SolverCtx<'_>, p: &ParamMap) -> Result<(), SolverError> {
        let Some(pt) = snapshot() else {
            return Err(SolverError::fatal(
                CALLER,
                "no particles (ParticleDynamics has not run)",
            ));
        };
        self.n_time = match p.get("output count") {
            Some(Value::Integer(v)) if !v.is_empty() => v[0] as i64,
            _ => self.n_time + 1,
        };
        let offset = match p.get("fileindex offset") {
            Some(Value::Integer(v)) if !v.is_empty() => v[0] as i64,
            _ => 0,
        };
        let itime = self.n_time + offset;
        let prefix = p
            .string("filename prefix")
            .unwrap_or("particles")
            .to_string();
        let min_status = match p.get("min status for saving") {
            Some(Value::Integer(v)) if !v.is_empty() => v[0],
            _ => PARTICLE_INITIATED,
        };
        let max_status = match p.get("max status for saving") {
            Some(Value::Integer(v)) if !v.is_empty() => v[0],
            _ => PARTICLE_LOST - 1,
        };
        let single = p.logical("single precision").unwrap_or(false);
        let (prec_bits, prec_size) = if single { (32, 4) } else { (64, 8) };
        let int_size = 4usize;
        let saved: Vec<usize> = (0..pt.status.len())
            .filter(|&i| pt.status[i] <= max_status && pt.status[i] >= min_status)
            .collect();
        let mut n = saved.len();
        if let Some(f) = p.real("particle save fraction") {
            n = (f * n as f64).round() as usize;
        } else if let Some(Value::Integer(v)) = p.get("particle save number") {
            if let Some(&k) = v.first() {
                n = n.min(k.max(0) as usize);
            }
        }
        let dir = p.string(MESH_DIRECTORY_KEY).ok_or_else(|| {
            SolverError::NotWired("the mesh directory of ParticleOutputSolver (elmer-run)".into())
        })?;
        let file = if itime < 10000 {
            format!("{dir}/{prefix}{itime:04}.vtu")
        } else {
            format!("{dir}/{prefix}{itime}.vtu")
        };
        ctx.info(
            "ParticleOutputVtu",
            &format!("Saving particles to file: {file}"),
            8,
        );
        ctx.warn("WriteVtuFile", "Are there really no scalars or vectors?");

        let mut out: Vec<u8> = Vec::new();
        let mut s = |t: &str| out.extend_from_slice(t.as_bytes());
        s("<?xml version=\"1.0\"?>\n");
        s("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        s("  <UnstructuredGrid>\n");
        s(&format!(
            "    <Piece NumberOfPoints=\"{n}\" NumberOfCells=\"{n}\">\n"
        ));
        s("      <PointData>\n");
        s("      </PointData>\n");
        s("      <CellData>\n");
        s("      </CellData>\n");
        s("      <Points>\n");
        let mut off = 0usize;
        s(&format!(
            "        <DataArray type=\"Float{prec_bits}\" NumberOfComponents=\"3\""
        ));
        s(&format!(" format=\"appended\" offset=\"{off}\"/>\n"));
        let kp = 3 * n * prec_size;
        off += int_size + kp;
        s("      </Points>\n");
        s("      <Cells>\n");
        let ki = n * int_size;
        for name in ["connectivity", "offsets", "types"] {
            s(&format!(
                "        <DataArray type=\"Int32\" Name=\"{name}\""
            ));
            s(&format!(" format=\"appended\" offset=\"{off}\"/>\n"));
            off += ki + int_size;
        }
        s("      </Cells>\n");
        s("    </Piece>\n");
        s("  </UnstructuredGrid>\n");
        s("<AppendedData encoding=\"raw\">\n");
        out.push(b'_');
        // Points: three components, those beyond the particle dimension 0.
        out.extend_from_slice(&(kp as i32).to_le_bytes());
        for &i in saved.iter().take(n) {
            let mut c = [0.0f64; 3];
            c[..pt.dim].copy_from_slice(&pt.coord[i][..pt.dim]);
            for v in c {
                if single {
                    out.extend_from_slice(&(v as f32).to_le_bytes());
                } else {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
        // connectivity, offsets, types (VTK_VERTEX).
        let cells: [fn(usize) -> i32; 3] = [|i| i as i32, |i| i as i32 + 1, |_| 1];
        for gen in cells {
            out.extend_from_slice(&(ki as i32).to_le_bytes());
            for i in 0..n {
                out.extend_from_slice(&gen(i).to_le_bytes());
            }
        }
        out.extend_from_slice(b"\n</AppendedData>\n");
        out.extend_from_slice(b"</VTKFile>\n");
        std::fs::File::create(&file)
            .and_then(|mut f| f.write_all(&out))
            .map_err(|e| SolverError::fatal("WriteVtuFile", format!("{file}: {e}")))?;
        Ok(())
    }
}
