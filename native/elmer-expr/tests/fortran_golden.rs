//! Property values captured from the pristine Fortran ElmerSolver.
//!
//! `fixtures/golden/case.sif` (tables from the quick tests SaveDependenceSpline
//! `mat.dat` and heateq_spline_table `cond.dat`, plus coordinate-, time- and
//! multi-variable dependencies) was run with SaveMaterials, which evaluates
//! `ListGetReal` at the mesh nodes, and written as Float64 VTU. Every ninth
//! node is stored in `fixtures/properties_golden.csv`.

use std::collections::HashMap;

use elmer_core::{Dependent, DependentKind, Value};
use elmer_expr::{
    compile_with, CompileOptions, Compiled, FieldLookup, MatcEngine, UdfRegistry,
};

struct NoMatc;
impl MatcEngine for NoMatc {
    fn eval(&self, _: &str, _: &[f64]) -> Result<Vec<f64>, String> {
        Err("no MATC in this test".into())
    }
}

struct Golden {
    xy: HashMap<u32, [f64; 2]>,
    temp: HashMap<u32, f64>,
}
impl FieldLookup for Golden {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        (var == "temp").then(|| self.temp[&node])
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        let [x, y] = self.xy[&node];
        [x, y, 0.0]
    }
    fn time(&self) -> f64 {
        1.0
    }
}

fn read_table(name: &str) -> (Vec<f64>, Vec<f64>) {
    let path = format!("{}/tests/fixtures/golden/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let mut x = Vec::new();
    let mut y = Vec::new();
    for line in text.lines() {
        let nums: Vec<f64> = line.split_whitespace().filter_map(|t| t.parse().ok()).collect();
        if nums.len() >= 2 {
            x.push(nums[0]);
            y.push(nums[1]);
        }
    }
    (x, y)
}

fn table(vars: &[&str], x: Vec<f64>, y: Vec<f64>, cubic: bool) -> Value {
    Value::Dependent(Dependent {
        variables: vars.iter().map(|s| s.to_string()).collect(),
        kind: DependentKind::Table { x, y, cubic },
        rows: 1,
        cols: 1,
    })
}

#[test]
fn tables_match_fortran_listgetreal_to_1e_14() {
    let csv = std::fs::read_to_string(format!(
        "{}/tests/fixtures/properties_golden.csv",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let mut nodes = Vec::new();
    let mut g = Golden { xy: HashMap::new(), temp: HashMap::new() };
    let mut want: Vec<Vec<f64>> = vec![Vec::new(); 8];
    for line in csv.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<f64> = line.split(',').map(|t| t.parse().unwrap()).collect();
        let node = f[0] as u32;
        nodes.push(node);
        g.xy.insert(node, [f[1], f[2]]);
        g.temp.insert(node, f[3]);
        for (k, w) in want.iter_mut().enumerate() {
            w.push(f[4 + k]);
        }
    }
    assert!(nodes.len() > 100);

    let (mx, my) = read_table("mat.dat");
    let (cx, cy) = read_table("cond.dat");
    let reg = UdfRegistry::new();
    let mono = CompileOptions { monotone: true };
    let plain = CompileOptions::default();
    let cases: Vec<(&str, Compiled)> = vec![
        ("lin sdep", compile_with(&table(&["temp"], mx.clone(), my.clone(), false), &reg, &plain)),
        ("cub sdep", compile_with(&table(&["temp"], mx.clone(), my.clone(), true), &reg, &plain)),
        ("mono sdep", compile_with(&table(&["temp"], mx, my, true), &reg, &mono)),
        ("cond cubic", compile_with(&table(&["temp"], cx.clone(), cy.clone(), true), &reg, &plain)),
        ("cond lin", compile_with(&table(&["temp"], cx, cy, false), &reg, &plain)),
        (
            "two var",
            compile_with(
                &table(&["coordinate 2, temp"], vec![0.0, 2.0, 5.0], vec![1.0, 3.0, -1.0], false),
                &reg,
                &plain,
            ),
        ),
        (
            "coord tab",
            compile_with(
                &table(
                    &["coordinate 1"],
                    vec![0.0, 0.5, 1.1, 2.0, 2.6],
                    vec![1.0, 2.0, 1.5, 4.0, 3.0],
                    true,
                ),
                &reg,
                &plain,
            ),
        ),
        (
            "time tab",
            compile_with(&table(&["time"], vec![0.0, 1.0], vec![10.0, 20.0], false), &reg, &plain),
        ),
    ]
    .into_iter()
    .map(|(n, c)| (n, c.unwrap()))
    .collect();

    for (k, (name, c)) in cases.iter().enumerate() {
        let mut out = vec![0.0; nodes.len()];
        c.eval_nodes(&nodes, &g, &NoMatc, &mut out).unwrap();
        let mut worst = 0.0_f64;
        for (i, (&got, &exp)) in out.iter().zip(&want[k]).enumerate() {
            let rel = if exp == 0.0 { got.abs() } else { ((got - exp) / exp).abs() };
            worst = worst.max(rel);
            assert!(rel <= 1e-14, "{name}: node {} got {got:e} want {exp:e} (rel {rel:e})", nodes[i]);
        }
        eprintln!("{name}: worst relative difference {worst:e}");
    }
}
