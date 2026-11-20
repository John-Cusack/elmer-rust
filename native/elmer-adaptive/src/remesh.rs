//! `External_ReMesh` (Adaptive.F90): the desired element size at each node
//! goes to `bgmesh.nodes`, an external mesher (`Mesh2D`, or the solver's
//! `Mesh Command`) meshes the geometry again with it, and the run reads the
//! new mesh with its own reader ([`SolverCtx::load_mesh`]).
//!
//! The mesher is an external program, as in Elmer: the command runs through
//! the shell in the case directory (`SystemCommand`), so `Mesh2D` comes from
//! the `PATH`. Only the file I/O path is ported; MMG (`Adaptive Remesh Use
//! MMG`) and Gmsh (`Adaptive Remesh Use Gmsh`) are refused.

use std::fmt::Write as _;
use std::path::Path;

use elmer_core::Params;
use elmer_solver::{SolverCtx, SolverError};

use crate::mesh::{e_format, write_mesh_to_disk, OwnedMesh};

/// `AEPS = 10 * EPSILON(1.0_dp)`.
const AEPS: f64 = 10.0 * f64::EPSILON;

fn unsupported(what: &str) -> SolverError {
    SolverError::Unsupported(format!("RefineMesh: {what}"))
}

/// What `External_ReMesh` takes from RefineMesh.
pub(crate) struct RemeshInput<'a> {
    /// `HValue`: the mean element diameter at each node.
    pub hvalue: &'a [f64],
    pub nodal_error: &'a [f64],
    pub h_convergence: &'a [f64],
    /// `ErrorLimit / 3`.
    pub error_limit: f64,
    pub min_h: f64,
    pub max_h: f64,
    pub max_change: f64,
    pub coarsening: bool,
    /// `RefMesh % AdaptiveDepth`.
    pub depth: usize,
    /// `OutputPath`: the case directory.
    pub case_dir: &'a Path,
}

/// The keywords of the remeshing path this port refuses.
pub(crate) fn check(p: &elmer_core::ParamMap, sim: &elmer_core::ParamMap, temperature: bool) -> Result<(), SolverError> {
    for k in ["adaptive remesh use mmg", "adaptive remesh use gmsh"] {
        if p.logical(k) == Some(true) {
            return Err(unsupported(k));
        }
    }
    if temperature && sim.get("gebhart factors").is_some() {
        return Err(unsupported("Gebhart Factors with Adaptive Remesh (the View program)"));
    }
    Ok(())
}

/// `HValueF`: the desired size at the nodes whose error is above `100 AEPS`
/// (`-1` elsewhere).
fn desired_h(inp: &RemeshInput<'_>) -> Vec<f64> {
    let mut out = vec![-1.0f64; inp.hvalue.len()];
    for (i, o) in out.iter_mut().enumerate() {
        let e = inp.nodal_error[i];
        if e > 100.0 * AEPS {
            let h = inp.hvalue[i];
            let mut lambda = (inp.error_limit / e).powf(1.0 / inp.h_convergence[i]);
            lambda = if inp.depth < 1 {
                h * lambda.min(1.33).max(0.75)
            } else {
                h * lambda.min(inp.max_change).max(1.0 / inp.max_change)
            };
            if !inp.coarsening {
                lambda = lambda.min(h);
            }
            if inp.max_h > 0.0 {
                lambda = lambda.min(inp.max_h);
            }
            if inp.min_h > 0.0 {
                lambda = lambda.max(inp.min_h);
            }
            *o = lambda;
        }
    }
    out
}

/// gfortran list-directed `WRITE(u,*) i` of a default integer.
fn list_int(i: usize) -> String {
    format!("{i:>12}")
}

/// Runs `cmd` through the shell in `dir`, as `SystemCommand` (C `system`)
/// does; its exit status is not checked there either.
fn system_command(ctx: &SolverCtx<'_>, cmd: &str, dir: &Path) -> Result<(), SolverError> {
    ctx.info("ReMesh", &format!("Meshing command: {cmd}"), 10);
    std::process::Command::new("sh")
        .arg("-c")
        .arg(cmd)
        .current_dir(dir)
        .status()
        .map_err(|e| SolverError::fatal("ReMesh", format!("cannot run '{cmd}': {e}")))?;
    Ok(())
}

/// `External_ReMesh`: the new mesh, read back from the mesher's output.
pub(crate) fn external_remesh(ctx: &mut SolverCtx<'_>, inp: &RemeshInput<'_>) -> Result<OwnedMesh, SolverError> {
    let dim = ctx.mesh.dim();
    if dim != 2 {
        return Err(unsupported("Adaptive Remesh in 3D (Mesh3D)"));
    }
    let hf = desired_h(inp);
    let p = ctx.params().clone();
    let name = match p.string("adaptive mesh name") {
        Some(n) => n.to_string(),
        None => format!("RefinedMesh{}", inp.depth + 1),
    };
    // `Path = OutputPath // '/' // Name`, with the case directory as
    // `OutputPath` and the commands run there.
    let path = format!("./{name}");

    ctx.info("ReMesh", "Saving background mesh density in point cloud format", 5);
    let mut s = String::new();
    let _ = writeln!(s, "{}", list_int(hf.iter().filter(|&&h| h > 0.0).count()));
    for (i, &h) in hf.iter().enumerate() {
        if !(h > 0.0) {
            continue;
        }
        let x = ctx.mesh.coord(i as u32);
        let _ = writeln!(s, "{}{}{}", e_format(x[0], 23, 15), e_format(x[1], 23, 15), e_format(h, 23, 15));
    }
    let _ = writeln!(s, "{}", list_int(0));
    std::fs::write(inp.case_dir.join("bgmesh.nodes"), s)
        .map_err(|e| SolverError::fatal("ReMesh", format!("bgmesh.nodes: {e}")))?;
    let out_dir = inp.case_dir.join(&name);
    write_mesh_to_disk(ctx.mesh, &out_dir, &|t| t)
        .map_err(|e| SolverError::fatal("WriteMeshToDisk", e.to_string()))?;

    let cmd = match p.string("mesh command") {
        Some(c) => c.to_string(),
        None => {
            let input = p
                .string("mesh input file")
                .or_else(|| ctx.simulation_params().string("mesh input file"))
                .ok_or_else(|| SolverError::fatal("ReMesh", "Mesh Input File not given"))?
                .to_string();
            // `Mesh % Name` of the run's first mesh.
            let mesh_name = ctx
                .sim
                .model
                .header
                .string("mesh name")
                .ok_or_else(|| SolverError::fatal("ReMesh", "no Mesh Name"))?
                .to_string();
            format!("Mesh2D ./{mesh_name}/{input} {path} --bgmesh=bgmesh.nodes")
        }
    };
    system_command(ctx, &cmd, inp.case_dir)?;
    if let Some(conv) = p.string("mesh conversion command") {
        let conv = format!("{conv} -out {path}");
        system_command(ctx, &conv, inp.case_dir)?;
    }
    let m = ctx.load_mesh(&out_dir)?;
    Ok(OwnedMesh::from_view(&*m))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn desired_size_clamps() {
        let dir = std::path::PathBuf::from(".");
        let inp = RemeshInput {
            hvalue: &[0.1, 0.1, 0.1],
            nodal_error: &[0.0, 1.0e-3, 1.0],
            h_convergence: &[1.0, 1.0, 1.0],
            error_limit: 1.0e-2,
            min_h: 0.0,
            max_h: 0.0,
            max_change: 3.0,
            coarsening: true,
            depth: 0,
            case_dir: &dir,
        };
        let h = desired_h(&inp);
        assert_eq!(h[0], -1.0);
        assert_eq!(h[1], 0.1 * 1.33);
        assert_eq!(h[2], 0.1 * 0.75);
        assert_eq!(list_int(37), "          37");
    }
}
