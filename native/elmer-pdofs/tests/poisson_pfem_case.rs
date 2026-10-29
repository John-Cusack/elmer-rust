//! PoissonPFEM and SD_PoissonPFEM end to end through `elmer-solver`: the case's own
//! `case.sif` (parsed with `elmer-sif`), [`PDofsBackend`] in
//! `Backends::dofs` and the case's `Poisson.f90` plugin as
//! [`elmer_pdofs::PoissonPfem`]. The solver norm must be the pristine run's
//! (`Reference Norm = 0.19674375E-01`; the run printed 1.96743752E-02) and
//! the values the pristine `case.result`.

mod support;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_pdofs::{register_poisson_pfem, PDofsBackend};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::{data_dir, read_perm, read_values, TestMesh};

/// `elmer-solver`'s view of the parsed case. Elmer's `ListAddString`
/// lowercases string values; `elmer-sif` keeps quoted text as written.
fn sif_model(path: &std::path::Path) -> Model {
    let parsed = elmer_sif::parse_file(path.to_str().unwrap()).unwrap();
    let m = parsed.sif.merged();
    let convert = |s: &elmer_sif::Section| -> ParamMap {
        let mut p = elmer_sif::section_to_param_map(s);
        for v in p.0.values_mut() {
            if let Value::Str(t) = v {
                *t = t.to_lowercase();
            }
        }
        p
    };
    let numbered = |v: &[elmer_sif::Section]| -> Vec<ParamMap> {
        let n = v.iter().map(|s| s.id as usize).max().unwrap_or(0);
        let mut out = vec![ParamMap::new(); n];
        for s in v {
            out[s.id as usize - 1] = convert(s);
        }
        out
    };
    let one = |s: &Option<elmer_sif::Section>| s.as_ref().map(convert).unwrap_or_default();
    Model {
        header: ParamMap::new(),
        simulation: one(&m.simulation),
        constants: one(&m.constants),
        solvers: numbered(&m.solvers),
        equations: numbered(&m.equations),
        bodies: numbered(&m.bodies),
        materials: numbered(&m.materials),
        body_forces: numbered(&m.body_forces),
        initial_conditions: numbered(&m.ics),
        boundary_conditions: numbered(&m.bcs),
    }
}

/// Runs `case` (`tests/data/<case>/case.sif` on the `mesh` directory's
/// mesh) and checks dofs, norm and the `case.result` values.
fn run_case(case: &str, (mesh, dim): (&str, usize), n_dofs: usize, run_norm: f64, reference: f64) {
    let dir = data_dir(case);
    let model = sif_model(&dir.join("case.sif"));
    let mesh = TestMesh::read(&data_dir(mesh), dim);
    let mut reg = SolverRegistry::new();
    register_poisson_pfem(&mut reg);
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("potential").unwrap();
    assert_eq!(u.values.len(), n_dofs);
    let norm = u.conv.norm;
    assert!((norm - run_norm).abs() < 1e-10, "{case}: norm {norm:e}");
    assert!(
        (norm - reference).abs() / reference < 1e-5,
        "{case}: norm {norm:e}"
    );

    let perm = u.perm.as_deref().unwrap();
    let elmer_perm = read_perm(&dir.join("perm.dat"));
    let ours: Vec<u32> = perm.iter().map(|&p| p + 1).collect();
    assert_eq!(ours, elmer_perm, "{case}: Perm");
    // case.result lists x(Perm(i)) slot by slot.
    let elmer = read_values(&dir.join("solution.dat"));
    let scale = elmer.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (slot, &v) in elmer.iter().enumerate() {
        let ours = u.values[perm[slot] as usize];
        assert!(
            (ours - v).abs() <= 1e-11 * scale,
            "{case}: slot {slot}: {ours:e} vs {v:e}"
        );
    }
}

/// `p:2` on full quads: 121 nodes + 220 edges + 100 global bubbles.
#[test]
fn poisson_pfem_runs_through_the_solver_framework() {
    run_case(
        "poisson_pfem",
        ("poisson_pfem", 2),
        441,
        1.96743752e-2,
        0.19674375e-1,
    );
}

/// `SD_PoissonPFEM`: the same case with the default serendipity elements,
/// so no bubbles (121 nodes + 220 edges) and the 8-point economic rule.
#[test]
fn sd_poisson_pfem_runs_through_the_solver_framework() {
    run_case(
        "sd_poisson_pfem",
        ("poisson_pfem", 2),
        341,
        2.23542169e-2,
        0.223544023689e-1,
    );
}

/// PoissonPFEM's case on triangles with `p:3` (a lane-made variant, see
/// `tests/data/README.md`): edge functions of degree 2 and 3 with their
/// orientation, one bubble per triangle, and the 12-point simplex rule on
/// the p reference triangle. 121 nodes + 320 edges x 2 + 200 bubbles.
#[test]
fn triangles_p3_run_through_the_solver_framework() {
    // The run printed 0.13356115E-01; this is the norm of its case.result.
    run_case(
        "tri_p3",
        ("tri_p3", 2),
        961,
        0.013356115428399276,
        0.13356115e-1,
    );
}

/// The same plugin on 3x3x3 bricks with `p:3` (full elements; a lane-made
/// variant, see `tests/data/README.md`): 64 nodes + 144 edges x 2 + 108
/// faces x 4 + 27 bubbles x 8, faces oriented by `getSquareFaceDirection`.
#[test]
fn bricks_p3_run_through_the_solver_framework() {
    // The run printed 0.23399267E-01; this is the norm of its case.result.
    run_case(
        "brick_p3",
        ("brick_p3", 3),
        1000,
        0.023399266709645154,
        0.23399267e-1,
    );
}

/// Serendipity bricks with `p:4` on the `brick_p3` mesh: 64 nodes + 144
/// edges x 3 + 108 faces x 1 (the serendipity face function), no bubbles.
#[test]
fn serendipity_bricks_p4_run_through_the_solver_framework() {
    // The run printed 0.59520241E-02; this is the norm of its case.result.
    run_case(
        "brick_sd_p4",
        ("brick_p3", 3),
        604,
        0.005952024123915765,
        0.59520241e-2,
    );
}
