//! `Elmer2OpenFoamWrite` (fem/src/modules/Elmer2OpenFoamIO.F90): the
//! `Target Variable` interpolated onto the cell centres of an OpenFOAM
//! case (`<OpenFOAM Directory>/0/C`) and written there as a
//! `volScalarField` (`OpenFOAM File`, default `fieldSolidHS.dat`).
//!
//! The centres become a mesh of nodes only, onto which `VariableGet`
//! interpolates the variable (`InterpolateMeshToMesh`; unfound points stay
//! zero, the permutation is the identity), and that variable becomes the
//! solver's: its norm is `ComputeChange`'s norm of it.
//!
//! Refused: OpenFOAM blocks (`0/*/C` without `0/C`), `OpenFOAM Mesh i`
//! lists, `Elmer Center Filename`, multi-component and non-nodal targets,
//! and parallel runs.

use std::fmt::Write as _;
use std::path::PathBuf;

use elmer_adaptive::mesh::OwnedMesh;
use elmer_core::{FieldKind, ParamMap, Params, Value};
use elmer_interp::{build_projector, FieldData, InterpOptions, Location, Locator};
use elmer_solver::{Solver, SolverCtx, SolverError, SolverRegistry};

pub(crate) const PROCEDURE: (&str, &str) = ("Elmer2OpenFoamIO", "Elmer2OpenFoamWrite");
const CALLER: &str = "Elmer2OpenFoamWrite";

pub(crate) fn is_section(p: &ParamMap) -> bool {
    elmer_solver::solver::procedure_of(p)
        .is_some_and(|(l, e)| l.eq_ignore_ascii_case(PROCEDURE.0) && e.eq_ignore_ascii_case(PROCEDURE.1))
}

pub(crate) fn register(registry: &mut SolverRegistry) {
    registry.register(PROCEDURE.0, PROCEDURE.1, |_| Ok(Box::new(Elmer2OpenFoamWrite { case_dir: PathBuf::new() })));
}

pub(crate) struct Elmer2OpenFoamWrite {
    pub case_dir: PathBuf,
}

fn fatal(m: impl Into<String>) -> SolverError {
    SolverError::fatal(CALLER, m)
}

fn unsupported(m: &str) -> SolverError {
    SolverError::Unsupported(format!("{CALLER}: {m}"))
}

/// `CreateFOAMMesh`: the cell centres of a `C` file (header up to
/// `internalField`, then the count and one `(x y z)` per line, or all of
/// them inline).
fn read_centres(text: &str, max_cells: Option<i32>) -> Result<Vec<[f64; 3]>, SolverError> {
    let lines: Vec<&str> = text.lines().collect();
    // The first 100 lines are searched for the header.
    let h = lines
        .iter()
        .take(100)
        .position(|l| l.contains("internalField"))
        .ok_or_else(|| fatal("Could not find > internalField < in header!"))?;
    let head = lines[h];
    if head.contains('(') {
        return Err(unsupported("cell centres inline with internalField"));
    }
    let mut n: usize = lines
        .get(h + 1)
        .and_then(|l| l.trim().parse().ok())
        .ok_or_else(|| fatal("Could not read number of nodes!"))?;
    if let Some(i) = max_cells.filter(|&i| i > 0 && (i as usize) < n) {
        n = i as usize;
    }
    // The `(` line, then one centre per line.
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let l = lines.get(h + 3 + i).ok_or_else(|| fatal(format!("Could not read coordinate line: {}", i + 1)))?;
        let j = l.find('(').ok_or_else(|| fatal(format!("Expecting a parenthesis at the start of OpenFOAM line: {}", i + 1)))?;
        let k = l.rfind(')').ok_or_else(|| fatal(format!("Expecting a parenthesis at the end of OpenFOAM line: {}", i + 1)))?;
        // `INDEX(ReadStr,'(',.TRUE.)`: the last parenthesis of each kind.
        let j = l.rfind('(').unwrap_or(j);
        let v: Vec<f64> = l[j + 1..k].split_whitespace().filter_map(|t| t.parse().ok()).collect();
        if v.len() < 3 {
            return Err(fatal(format!("Could not read coordinate values: {}", i + 1)));
        }
        out.push([v[0], v[1], v[2]]);
    }
    Ok(out)
}

/// `WriteOFField`.
fn field_file(values: &[f64]) -> String {
    let mut s = String::new();
    s.push_str(concat!(
        "/*--------------------------------*- C++ -*----------------------------------*\\ \n",
        "| =========                 |                                                 | \n",
        "| \\\\      /  F ield         | OpenFOAM: The OPEN Source CFD Toolbox           | \n",
        "|  \\\\    /   O peration     | Version:  dev                                   | \n",
        "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      | \n",
        "|    \\\\/     M anipulation  |                                                 | \n",
        "|                                                                             | \n",
        "|  This file written by Elmer2OpenFOAM coupling module in Elmer               | \n",
        "\\*---------------------------------------------------------------------------*/ \n",
        "FoamFile                                \n",
        "{                                       \n",
        "    version     2.0;                    \n",
        "    format      ascii;                  \n",
        "    class       volScalarField;         \n",
        "    object      fieldSolidHS;           \n",
        "}\n",
        "dimensions      [1 -1 -3 0 0 0 0];      \n",
        "boundaryField                           \n",
        "{                                       \n",
        "    \".*\"                                \n",
        "    {                                   \n",
        "        TYPE zeroGradient;              \n",
        "    }                                   \n",
        "}\n",
        "internalField   nonuniform List<scalar> \n",
    ));
    let _ = writeln!(s, "{}", values.len());
    s.push_str("(\n");
    for &v in values {
        let _ = writeln!(s, "{}", elmer_fixtures::udf::list_real(v));
    }
    s.push_str(");\n");
    s
}

impl Solver for Elmer2OpenFoamWrite {
    fn solve(&mut self, ctx: &mut SolverCtx<'_>, _dt: f64, _transient: bool) -> Result<(), SolverError> {
        let p = ctx.params().clone();
        if ctx.comm.size() > 1 {
            return Err(unsupported("parallel run"));
        }
        let name = p.string("target variable").ok_or_else(|| fatal("> Target Variable < must exist for the solver!"))?;
        let name = elmer_core::normalize_key(name);
        let id = ctx.sim.vars.id(&name).ok_or_else(|| fatal("Variable does not exist in Elmer mesh!"))?;
        if p.get("elmer center filename").is_some() {
            return Err(unsupported("Elmer Center Filename"));
        }
        if p.get("openfoam mesh 1").is_some() {
            return Err(unsupported("OpenFOAM Mesh lists"));
        }
        let base = p.string("openfoam directory").ok_or_else(|| fatal("> OpenFOAM Directory < must exist for the solver!"))?;
        // `OpenFOAMBlocks`: `<dir>/0/C`, or one block per `<dir>/0/*/C`.
        let zero = self.case_dir.join(base).join("0");
        if !zero.exists() {
            return Err(fatal(format!("OpenFOAM mesh does not exist: {}/0/", base)));
        }
        let mut blocks: Vec<PathBuf> = std::fs::read_dir(&zero)
            .map_err(|e| fatal(e.to_string()))?
            .filter_map(|d| d.ok().map(|d| d.path()))
            .filter(|d| d.join("C").is_file())
            .collect();
        blocks.sort();
        let dir = if zero.join("C").is_file() {
            zero.clone()
        } else if blocks.len() == 1 {
            blocks.remove(0)
        } else if blocks.is_empty() {
            return Err(fatal("No OpenFOAM mesh blocks found!"));
        } else {
            return Err(unsupported("several OpenFOAM mesh blocks"));
        };
        let text = std::fs::read_to_string(dir.join("C")).map_err(|e| fatal(format!("Could not open file for reading: {e}")))?;
        let mut centres = read_centres(&text, p.integer("number of cells"))?;
        if let Some(k) = p.integer("passive openfoam coordinate").filter(|k| (1..=3).contains(k)) {
            for c in &mut centres {
                c[k as usize - 1] = 0.0;
            }
        }

        let v = ctx.sim.vars.by_id(id);
        if v.kind != FieldKind::Nodal || v.dofs != 1 {
            return Err(unsupported("a target that is not a scalar nodal variable"));
        }
        let old = FieldData {
            name: name.clone(),
            location: Location::Nodes,
            values: v.values.clone(),
            perm: Some(v.perm.as_deref().map_or_else(|| (0..ctx.mesh.n_nodes() as u32).collect(), |q| q.to_vec())),
            prev_values: Vec::new(),
        };
        let n = centres.len();
        let target = OwnedMesh { dim: 3, coords: centres, elems: Vec::new(), n_bulk: 0 };
        // `Use Mesh Projector = False` for the duration.
        let mut sim_params = ctx.sim.model.simulation.clone();
        sim_params.insert("Use Mesh Projector", Value::Logical(false));
        let opts = InterpOptions::from_params(&sim_params, None);
        let loc = Locator::new(ctx.mesh, opts).map_err(|e| fatal(e.to_string()))?;
        let fe = ctx
            .sim
            .backends
            .fe
            .as_ref()
            .ok_or_else(|| SolverError::NotWired("basis and quadrature (elmer-fe)".into()))?;
        let extras = elmer_interp::Extras { body_mask: None, new_mask_perm: None, old_dg: None, new_dg: None };
        let projector = build_projector(&loc, &target, fe.basis.as_ref(), &extras, ctx.log).map_err(|e| fatal(e.to_string()))?;
        let mut new = vec![FieldData {
            name: name.clone(),
            location: Location::Nodes,
            values: vec![0.0; n],
            perm: Some((0..n as u32).collect()),
            prev_values: Vec::new(),
        }];
        projector.apply(std::slice::from_ref(&old), &mut new);
        let values = std::mem::take(&mut new[0].values);

        let file = p.string("openfoam file").unwrap_or("fieldSolidHS.dat");
        std::fs::write(dir.join(file), field_file(&values)).map_err(|e| fatal(format!("Could not write: {e}")))?;

        // `Solver % Variable => OFVar`: ComputeChange's norm of it.
        let spec = elmer_solver::NormSpec::from_params(&p, 1);
        let nrm = elmer_solver::compute_norm(&spec, &values, n, ctx.comm);
        let sid = ctx.solver_id();
        if let Some(slot) = ctx.sim.slots.iter_mut().find(|s| s.id == sid) {
            slot.null_norm = nrm;
        }
        Ok(())
    }
}
