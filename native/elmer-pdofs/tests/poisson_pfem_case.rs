//! PoissonPFEM end to end through `elmer-solver`: the case's own
//! `case.sif` (parsed with `elmer-sif`), [`PDofsBackend`] in
//! `Backends::dofs` and the case's `Poisson.f90` plugin as
//! [`elmer_pdofs::PoissonPfem`]. The solver norm must be the pristine run's
//! (`Reference Norm = 0.19674375E-01`; the run printed 1.96743752E-02) and
//! the values the pristine `case.result`.

mod support;

use elmer_core::{NullLogger, ParamMap, SerialComm, Value};
use elmer_pdofs::{register_poisson_pfem, PDofsBackend};
use elmer_solver::{Model, Simulation, SolverRegistry};
use support::{data_dir, read_values, TestMesh};

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

#[test]
fn poisson_pfem_runs_through_the_solver_framework() {
    let dir = data_dir("poisson_pfem");
    let model = sif_model(&dir.join("case.sif"));
    let mesh = TestMesh::read(&dir, 2);
    let mut reg = SolverRegistry::new();
    register_poisson_pfem(&mut reg);
    let mut sim = Simulation::new(model);
    sim.backends.dofs = Some(Box::new(PDofsBackend::new()));
    sim.setup(&mesh, &reg, &NullLogger, &SerialComm).unwrap();
    sim.run(&mesh, &NullLogger, &SerialComm).unwrap();

    let u = sim.vars.get("potential").unwrap();
    // 121 nodes + 220 edges + 100 global bubbles.
    assert_eq!(u.values.len(), 441);
    let norm = u.conv.norm;
    assert!((norm - 1.96743752e-2).abs() < 1e-10, "norm {norm:e}");
    assert!((norm - 0.19674375e-1).abs() / 0.19674375e-1 < 1e-5);

    // case.result lists x(Perm(i)) slot by slot.
    let perm = u.perm.as_deref().unwrap();
    let elmer = read_values(&dir.join("solution.dat"));
    let scale = elmer.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for (slot, &v) in elmer.iter().enumerate() {
        let ours = u.values[perm[slot] as usize];
        assert!(
            (ours - v).abs() <= 1e-11 * scale,
            "slot {slot}: {ours:e} vs {v:e}"
        );
    }
}
