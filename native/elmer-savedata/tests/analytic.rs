//! Values against independent geometry: affine and quadratic fields on
//! meshes whose integrals are known in closed form, masks, axisymmetric
//! weighting, vector fields, exported variables, scheduling cadence and
//! append behaviour.

mod support;

use std::f64::consts::PI;
use std::path::{Path, PathBuf};

use elmer_core::{ParamMap, Params, Value};
use elmer_solver::Model;
use support::harness::*;
use support::meshio::{read_elmer_mesh, unit_square, Mesh};

fn two_cell() -> Mesh {
    read_elmer_mesh(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/two_cell"), 2)
}

fn model() -> Model {
    Model {
        simulation: sec(&[("Coordinate System", st("Cartesian 2D"))]),
        bodies: vec![
            sec(&[("Target Bodies", ints(&[1]))]),
            sec(&[("Target Bodies", ints(&[2])), ("Save This", yes())]),
        ],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Save This", yes()), ("Save Line", yes())]),
            sec(&[("Target Boundaries", ints(&[3])), ("Flux Integrate", yes())]),
            sec(&[("Target Boundaries", ints(&[4]))]),
            sec(&[("Target Boundaries", ints(&[2]))]),
        ],
        ..Model::default()
    }
}

fn scalars(out: &Path, ops: &[(&str, &str)]) -> ParamMap {
    let mut s = sec(&[
        ("Equation", st("savescalars")),
        ("Filename", st("s.dat")),
        ("Output Directory", dir_str(out)),
    ]);
    for (i, (var, op)) in ops.iter().enumerate() {
        if !var.is_empty() {
            s.insert(&format!("Variable {}", i + 1), st(var));
        }
        s.insert(&format!("Operator {}", i + 1), st(op));
    }
    s
}

/// The last row of `s.dat`.
fn row(out: &Path, file: &str) -> Vec<f64> {
    let text = std::fs::read_to_string(out.join(file)).unwrap();
    text.lines()
        .last()
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

fn rows(out: &Path, file: &str) -> usize {
    std::fs::read_to_string(out.join(file)).map_or(0, |t| t.lines().count())
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0), "{a} vs {b}");
}

fn names(out: &Path, file: &str) -> Vec<String> {
    std::fs::read_to_string(out.join(format!("{file}.names")))
        .unwrap()
        .lines()
        .skip_while(|l| !l.starts_with("Variables in columns"))
        .skip(1)
        .map(|l| l.split_once(": ").unwrap().1.to_string())
        .collect()
}

#[test]
fn nodal_mean_differs_from_the_volume_weighted_mean() {
    let out = tmp("an-mean");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c.hook("SaveScalars", scalars(&out, &[("temperature", "mean"), ("", "int mean"), ("", "volume"), ("", "int")])).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    // Nodes at x = 0 1 3 0 1 3; the cells have widths 1 and 2.
    near(r[0], 8.0 / 6.0);
    near(r[1], 1.5);
    near(r[2], 3.0);
    near(r[3], 4.5);
    assert_eq!(
        names(&out, "s.dat"),
        ["mean: temperature", "int mean: temperature", "volume: temperature", "int: temperature"]
    );
}

#[test]
fn constants_and_affine_fields_integrate_exactly() {
    let out = tmp("an-affine");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![2.0 * x[0] + 3.0 * x[1] + 1.0]);
    c.field("Const", 1, |_| vec![7.0]);
    let ops = [
        ("temperature", "int"),
        ("const", "int"),
        ("const", "int mean"),
        ("temperature", "body int"),
        ("temperature", "int square"),
    ];
    let mut s = scalars(&out, &ops);
    s.insert("Mask Name 4", st("Save This"));
    // Mask only for the 4th; "int square" reads the mask name as well but is unmasked.
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    near(r[0], 16.5); // 2*4.5 + 3*1.5 + 3
    near(r[1], 21.0);
    near(r[2], 7.0);
    // Body 2 is x in [1, 3]: 2*4 + 3*1 + 2
    near(r[3], 13.0);
    // int (2x + 3y + 1)^2 over [0,3]x[0,1], exact with the 2x2 rule.
    near(r[4], 102.0);
    assert_eq!(names(&out, "s.dat")[3], "body int: temperature mask save this");
}

#[test]
fn boundary_integrals_use_the_boundary_conditions_flagged_for_them() {
    let out = tmp("an-boundary");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![2.0 * x[0] + 3.0 * x[1] + 1.0]);
    let mut s = scalars(&out, &[("temperature", "boundary int"), ("", "boundary int mean"), ("", "area")]);
    s.insert("Mask Name 1", st("Save This"));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    // BC 1 (bottom, y = 0) is flagged `Save This`, BC 2 (top) `Flux Integrate`.
    // Values come per flagged BC, in BC order: bottom (y = 0), top (y = 1).
    near(r[0], 12.0); // bottom: int (2x + 1) dx over [0, 3]
    near(r[1], 21.0); // top: int (2x + 4) dx
    near(r[2], 4.0);
    near(r[3], 7.0);
    near(r[4], 3.0);
    near(r[5], 3.0);
}

#[test]
fn axisymmetric_weighting_gives_the_volume_of_revolution() {
    let out = tmp("an-axi");
    let mut m = model();
    m.simulation = sec(&[("Coordinate System", st("Axi Symmetric"))]);
    let mut c = Case::new(two_cell(), m);
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c.hook("SaveScalars", scalars(&out, &[("temperature", "volume"), ("", "int")])).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    near(r[0], 2.0 * PI * 4.5); // 2 pi int r dr dz, r in [0, 3], z in [0, 1]
    near(r[1], 2.0 * PI * 9.0);
}

#[test]
fn quadratic_elements_integrate_a_quadratic_field_exactly() {
    let out = tmp("an-quad");
    let mut c = Case::new(unit_square(2, true, true), model());
    c.field("Temperature", 1, |x| vec![x[0] * x[0] + x[1]]);
    let ix = c.hook("SaveScalars", scalars(&out, &[("temperature", "int"), ("", "volume"), ("", "max")])).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    near(r[0], 1.0 / 3.0 + 0.5);
    near(r[1], 1.0);
    near(r[2], 2.0);
}

#[test]
fn vector_fields_use_the_norm_and_components_are_strided_views() {
    let out = tmp("an-vector");
    let mut c = Case::new(two_cell(), model());
    c.field("Velocity", 2, |_| vec![3.0, 4.0]);
    let ops = [
        ("velocity", "max"),
        ("velocity", "int"),
        ("velocity 2", "max"),
        ("velocity 1", "sum"),
        ("velocity", "dofs"),
    ];
    let ix = c.hook("SaveScalars", scalars(&out, &ops)).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    near(r[0], 5.0);
    near(r[1], 15.0);
    near(r[2], 4.0);
    near(r[3], 18.0);
    near(r[4], 12.0);
}

#[test]
fn deviation_and_the_other_statistics() {
    let out = tmp("an-stats");
    let mut c = Case::new(two_cell(), model());
    c.field("T", 1, |x| vec![x[0]]);
    let ops = [
        ("t", "deviation"),
        ("t", "sum square"),
        ("t", "rms"),
        ("t", "mean square"),
        ("t", "range"),
        ("t", "variance"),
        ("t", "min abs"),
    ];
    let ix = c.hook("SaveScalars", scalars(&out, &ops)).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    // x = 0 1 3 0 1 3: mean 4/3.
    let mean = 4.0 / 3.0;
    let dev = [0.0, 1.0, 3.0, 0.0, 1.0, 3.0].iter().map(|x: &f64| (x - mean).abs()).sum::<f64>() / 6.0;
    near(r[0], dev);
    near(r[1], 20.0);
    near(r[2], (20.0f64 / 6.0).sqrt());
    near(r[3], 20.0 / 6.0);
    near(r[4], 3.0);
    // The "variance" operator is the standard deviation.
    near(r[5], (20.0f64 / 6.0 - mean * mean).sqrt());
    near(r[6], 0.0);
}

#[test]
fn global_variables_simulation_scalars_and_result_prefix_are_saved_as_they_are() {
    let out = tmp("an-globals");
    let mut m = model();
    m.simulation.insert("my constant", real(4.5));
    m.simulation.insert("res: from another solver", real(2.5));
    let mut c = Case::new(two_cell(), m);
    c.sim.vars.set_scalar("timestep", 3.0);
    c.sim.vars.set_scalar("time", 1.5);
    let mut s = scalars(&out, &[("timestep", "max"), ("time", "max"), ("my constant", "max")]);
    s.insert("Variable 4", st("no such variable"));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    let text = std::fs::read_to_string(out.join("s.dat")).unwrap();
    // The timestep is an I10 column; the rest ES22.12E3.
    assert!(text.starts_with("         3"), "{text}");
    let r = row(&out, "s.dat");
    near(r[0], 3.0);
    near(r[1], 1.5);
    near(r[2], 4.5);
    near(r[3], 2.5);
    assert_eq!(
        names(&out, "s.dat"),
        [
            "value: timestep scalar variable",
            "value: time scalar variable",
            "my constant",
            "res: from another solver"
        ]
    );
    assert!(c.log.text().contains("Requested variable does not exist: no such variable"));
}

#[test]
fn target_variables_reach_the_next_hook_and_every_value_the_simulation_section() {
    let out = tmp("an-target");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let mut a = scalars(&out, &[("temperature", "int")]);
    a.insert("Target Variable 1", st("Temperature Integral"));
    let mut b = scalars(&out, &[("temperature integral", "max")]);
    b.insert("Filename", st("b.dat"));
    let (ia, ib) = (c.hook("SaveScalars", a).unwrap(), c.hook("SaveScalars", b).unwrap());
    c.visit(ia, false).unwrap();
    c.visit(ib, false).unwrap();
    near(row(&out, "b.dat")[0], 4.5);
    assert_eq!(names(&out, "b.dat"), ["value: temperature integral scalar variable"]);
    near(c.sim.vars.scalar("temperature integral").unwrap(), 4.5);
    near(c.sim.model.simulation.real("int: temperature").unwrap(), 4.5);
}

#[test]
fn a_target_variable_that_names_a_field_is_not_overwritten() {
    let out = tmp("an-target-clash");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let mut a = scalars(&out, &[("temperature", "int")]);
    a.insert("Target Variable 1", st("Temperature"));
    let ix = c.hook("SaveScalars", a).unwrap();
    let err = c.visit(ix, false).unwrap_err().to_string();
    assert!(err.contains("not a global scalar"), "{err}");
}

#[test]
fn exec_intervals_set_the_cadence_and_the_file_grows_one_row_per_visit() {
    let out = tmp("an-cadence");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let mut s = scalars(&out, &[("temperature", "max")]);
    s.insert("Exec Intervals", ints(&[2]));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.sim.vars.set_scalar("timestep interval", 1.0);
    for t in 1..=5 {
        c.sim.vars.set_scalar("timestep", t as f64);
        c.visit(ix, true).unwrap();
    }
    // Timesteps 1, 3 and 5.
    assert_eq!(rows(&out, "s.dat"), 3);
    assert_eq!(c.sim.slots[ix].times_visited, 3);
}

#[test]
fn file_append_keeps_an_earlier_run_and_otherwise_starts_over() {
    let out = tmp("an-append");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    for append in [false, true, true, false] {
        let mut s = scalars(&out, &[("temperature", "max")]);
        s.insert("File Append", Value::Logical(append));
        let ix = c.hook("SaveScalars", s).unwrap();
        c.visit(ix, false).unwrap();
        c.visit(ix, false).unwrap();
        c.sim.slots.clear();
    }
    // The last run truncated the file and wrote its own two rows.
    assert_eq!(rows(&out, "s.dat"), 2);
    let out2 = tmp("an-append2");
    for _ in 0..2 {
        let mut s = scalars(&out2, &[("temperature", "max")]);
        s.insert("File Append", yes());
        let ix = c.hook("SaveScalars", s).unwrap();
        c.visit(ix, false).unwrap();
        c.sim.slots.clear();
    }
    assert_eq!(rows(&out2, "s.dat"), 2);
}

#[test]
fn without_a_filename_the_values_are_echoed() {
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let mut s = sec(&[("Equation", st("savescalars")), ("Variable 1", st("temperature")), ("Operator 1", st("max"))]);
    s.insert("Operator 2", st("min"));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    let t = c.log.text();
    assert!(t.contains("Showing computed results:"), "{t}");
    assert!(t.contains("   1: max: temperature   3.000000000000E+000"), "{t}");
    assert!(t.contains("   2: min: temperature   0.000000000000E+000"), "{t}");
}

#[test]
fn boundary_statistics_over_a_flagged_boundary_condition() {
    let out = tmp("an-bstat");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0] + 10.0 * x[1]]);
    let mut s = scalars(&out, &[("temperature", "boundary max"), ("", "boundary min"), ("", "boundary sum"), ("", "boundary mean"), ("", "boundary dofs")]);
    s.insert("Mask Name 1", st("Save This"));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    let r = row(&out, "s.dat");
    // Bottom nodes x = 0, 1, 3 at y = 0.
    near(r[0], 3.0);
    near(r[1], 0.0);
    near(r[2], 4.0);
    near(r[3], 4.0 / 3.0);
    near(r[4], 3.0);
    assert_eq!(names(&out, "s.dat")[0], "boundary max: temperature over bc save this");
}

#[test]
fn a_boundary_operator_without_a_flagged_boundary_warns_and_skips() {
    let out = tmp("an-bnone");
    let mut c = Case::new(unit_square(2, true, false), model_square());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let mut s = scalars(&out, &[("temperature", "boundary sum"), ("", "max")]);
    s.insert("Mask Name 1", st("No Such Mask"));
    let ix = c.hook("SaveScalars", s).unwrap();
    c.visit(ix, false).unwrap();
    assert!(c.log.text().contains("No flag > no such mask < active"));
    assert_eq!(names(&out, "s.dat"), ["max: temperature"]);
}

#[test]
fn counts_and_bounding_box() {
    let out = tmp("an-counts");
    let mut c = Case::new(two_cell(), model());
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c
        .hook("SaveScalars", scalars(&out, &[("temperature", "nodes"), ("", "elements"), ("", "bounding box"), ("", "partitions")]))
        .unwrap();
    c.visit(ix, false).unwrap();
    assert_eq!(row(&out, "s.dat"), [6.0, 2.0, 0.0, 3.0, 0.0, 1.0, 1.0]);
    assert_eq!(
        names(&out, "s.dat"),
        [
            "nodes: temperature",
            "elements: temperature",
            "bounding box min x",
            "bounding box max x",
            "bounding box min y",
            "bounding box max y",
            "value: number of partitions"
        ]
    );
}

#[test]
fn saveline_writes_ordered_unique_nodes_with_labels() {
    let out = tmp("an-line");
    let mut c = Case::new(unit_square(3, true, false), model_square());
    c.field("Temperature", 1, |x| vec![x[0] * 10.0 + x[1]]);
    let sec_line = sec(&[
        ("Equation", st("saveline")),
        ("Filename", st("l.dat")),
        ("Variable 1", st("temperature")),
        ("Output Directory", dir_str(&out)),
    ]);
    let ix = c.hook("SaveLine", sec_line).unwrap();
    c.visit(ix, false).unwrap();
    // BC 1 is the bottom edge (tag 1): four nodes along y = 0.
    let text = std::fs::read_to_string(out.join("l.dat")).unwrap();
    assert_eq!(text.lines().count(), 4);
    let mut xs = Vec::new();
    for l in text.lines() {
        let f: Vec<f64> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
        // call count, bc, node, x, y, z, T
        assert_eq!(f[0], 1.0);
        assert_eq!(f[1], 1.0);
        xs.push(f[3]);
        near(f[6], f[3] * 10.0);
    }
    // A monotone traversal of the line, from either end.
    assert!(xs.windows(2).all(|w| w[0] < w[1]) || xs.windows(2).all(|w| w[0] > w[1]), "{xs:?}");
    let n = std::fs::read_to_string(out.join("l.dat.names")).unwrap();
    let cols: Vec<&str> = n
        .lines()
        .skip_while(|l| !l.starts_with("Data on different columns"))
        .skip(1)
        .map(|l| l.split_once(": ").unwrap().1)
        .collect();
    assert_eq!(
        cols,
        ["Call count", "Boundary condition", "Node index", "coordinate 1", "coordinate 2", "coordinate 3", "temperature"]
    );
    // A second visit appends and counts the call.
    c.visit(ix, false).unwrap();
    let text = std::fs::read_to_string(out.join("l.dat")).unwrap();
    assert_eq!(text.lines().count(), 8);
    assert!(text.lines().nth(4).unwrap().starts_with("2 "));
}

fn model_square() -> Model {
    Model {
        simulation: sec(&[("Coordinate System", st("Cartesian 2D"))]),
        bodies: vec![sec(&[("Target Bodies", ints(&[1]))])],
        boundary_conditions: vec![
            sec(&[("Target Boundaries", ints(&[1])), ("Save Line", yes())]),
            sec(&[("Target Boundaries", ints(&[2]))]),
        ],
        ..Model::default()
    }
}

#[test]
fn saveline_replaces_nan_and_writes_nothing_for_an_unflagged_mask() {
    let out = tmp("an-line-nan");
    let mut c = Case::new(unit_square(2, true, false), model_square());
    c.field("Temperature", 1, |x| vec![if x[0] == 0.5 { f64::NAN } else { x[0] }]);
    let mk = |mask: Option<&str>, file: &str| {
        let mut s = sec(&[
            ("Equation", st("saveline")),
            ("Filename", st(file)),
            ("Variable 1", st("temperature")),
            ("Output Directory", dir_str(&out)),
        ]);
        if let Some(m) = mask {
            s.insert("Save Mask", st(m));
        }
        s
    };
    let ia = c.hook("SaveLine", mk(None, "a.dat")).unwrap();
    let ib = c.hook("SaveLine", mk(Some("Nothing Has This"), "b.dat")).unwrap();
    c.visit(ia, true).unwrap();
    c.visit(ib, true).unwrap();
    let text = std::fs::read_to_string(out.join("a.dat")).unwrap();
    assert!(text.lines().any(|l| l.trim_end().ends_with("-1.00000000000E+000")), "{text}");
    assert!(c.log.text().contains("Replaced 1 NaN entries with -1"));
    // The names file says there are no nodes; no data file was started.
    assert!(!out.join("b.dat").exists());
    let names = std::fs::read_to_string(out.join("b.dat.names")).unwrap();
    assert!(names.contains("  bc nodes: 0"), "{names}");
}

#[test]
fn integrals_in_unsupported_coordinate_systems_are_refused() {
    let out = tmp("an-polar");
    let mut m = model();
    m.simulation = sec(&[("Coordinate System", st("Polar 2D"))]);
    let mut c = Case::new(two_cell(), m);
    c.field("Temperature", 1, |x| vec![x[0]]);
    let ix = c.hook("SaveScalars", scalars(&out, &[("temperature", "int")])).unwrap();
    let err = c.visit(ix, false).unwrap_err();
    assert!(err.to_string().contains("polar 2d"), "{err}");
}
