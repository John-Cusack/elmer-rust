use std::sync::Arc;

use elmer_core::{Dependent, DependentKind, Value};
use elmer_expr::udf::{fortran_mangle, normalize_library, split_procedure};
use elmer_expr::{
    compile, compile_with, CompileOptions, ExprError, FieldLookup, MatcEngine, Udf, UdfRegistry,
    MISSING,
};

/// Nodes on a line, x = node; "temperature" = 10 * node, absent at node 3.
struct Fields {
    time: f64,
}
impl FieldLookup for Fields {
    fn value(&self, var: &str, node: u32) -> Option<f64> {
        match var {
            "temperature" if node != 3 => Some(10.0 * node as f64),
            "pressure" => Some(2.0),
            _ => None,
        }
    }
    fn coord(&self, node: u32) -> [f64; 3] {
        [node as f64, 0.5, 0.0]
    }
    fn time(&self) -> f64 {
        self.time
    }
}

/// Handles `tx(0)*K` and `tx(0)+tx(1)` and `[tx(0), 2*tx(0); 3, 4]`; records the bound name.
struct Stub;
impl MatcEngine for Stub {
    fn eval(&self, expr: &str, a: &[f64]) -> Result<Vec<f64>, String> {
        self.eval_with("tx", expr, a)
    }
    fn eval_with(&self, var: &str, expr: &str, a: &[f64]) -> Result<Vec<f64>, String> {
        match (var, expr) {
            ("tx", "tx(0)*3") => Ok(vec![a[0] * 3.0]),
            ("tx", "tx(0)+tx(1)") => Ok(vec![a[0] + a[1]]),
            ("tx", "tensor") => Ok(vec![a[0], 2.0 * a[0], 3.0, 4.0]),
            ("st", "st(0)*2") => Ok(vec![a[0] * 2.0]),
            ("tx", "const") => Ok(vec![42.0]),
            _ => Err(format!("stub: {var} {expr}")),
        }
    }
}

fn dep(vars: &[&str], kind: DependentKind, rows: usize, cols: usize) -> Value {
    Value::Dependent(Dependent {
        variables: vars.iter().map(|s| s.to_string()).collect(),
        kind,
        rows,
        cols,
    })
}

fn table(vars: &[&str], x: &[f64], y: &[f64], cubic: bool) -> Value {
    dep(vars, DependentKind::Table { x: x.to_vec(), y: y.to_vec(), cubic }, 1, 1)
}

fn eval(c: &elmer_expr::Compiled, nodes: &[u32], time: f64) -> Vec<f64> {
    let mut out = vec![0.0; nodes.len() * c.ncomp()];
    c.eval_nodes(nodes, &Fields { time }, &Stub, &mut out).unwrap();
    out
}

#[test]
fn constants_and_tensors() {
    let c = compile(&Value::real(2.5));
    assert_eq!(eval(&c, &[0, 1, 2], 0.0), vec![2.5; 3]);
    assert_eq!(c.constant(), Some(&[2.5][..]));

    let k = Value::Real { rows: 3, cols: 3, data: (1..=9).map(f64::from).collect() };
    let c = compile(&k);
    assert_eq!((c.rows(), c.cols()), (3, 3));
    let out = eval(&c, &[4, 5], 0.0);
    assert_eq!(out.len(), 18);
    assert_eq!(&out[9..], &(1..=9).map(f64::from).collect::<Vec<_>>()[..]);
}

#[test]
fn non_real_values_fail_on_evaluation() {
    let c = compile(&Value::Str("x".into()));
    let mut out = [0.0];
    assert!(matches!(
        c.eval_nodes(&[0], &Fields { time: 0.0 }, &Stub, &mut out),
        Err(ExprError::Unsupported(_))
    ));
}

#[test]
fn linear_table_on_field_with_missing_node() {
    let v = table(&["temperature"], &[0.0, 10.0, 40.0], &[1.0, 2.0, 5.0], false);
    let c = compile(&v);
    // nodes 0..=4 -> T = 0,10,20,30(missing at 3),40
    let out = eval(&c, &[0, 1, 2, 3, 4], 0.0);
    assert_eq!(out[0], 1.0);
    assert_eq!(out[1], 2.0);
    assert_eq!(out[2], 2.0 + 10.0 / 30.0 * 3.0);
    assert_eq!(out[3], 0.0, "missing dependency gives zero");
    assert_eq!(out[4], 5.0);
}

#[test]
fn table_extrapolates_linearly_and_one_point_is_proportional() {
    let c = compile(&table(&["time"], &[0.0, 1.0], &[1.0, 3.0], false));
    assert_eq!(eval(&c, &[0], 2.0), vec![5.0]);
    let c = compile(&table(&["time"], &[1.0], &[4.0], false));
    assert_eq!(eval(&c, &[0], 2.5), vec![10.0]);
}

#[test]
fn cubic_table_only_with_more_than_three_points() {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 1.0, 0.0, 1.0];
    let c = compile(&table(&["time"], &x, &y, true));
    let r = elmer_expr::curve::cubic_spline(&x, &y, false);
    let want = elmer_expr::curve::interpolate_curve(&x, &y, 0.5, Some(&r));
    assert_eq!(eval(&c, &[0], 0.5)[0].to_bits(), want.to_bits());
    assert_ne!(want, 0.5);

    let c3 = compile(&table(&["time"], &x[..3], &y[..3], true));
    assert_eq!(eval(&c3, &[0], 0.5), vec![0.5], "three points stay linear");
}

#[test]
fn monotone_option_changes_the_spline() {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0];
    let y = [0.0, 0.1, 5.0, 5.1, 5.2];
    let v = table(&["time"], &x, &y, true);
    let reg = UdfRegistry::new();
    let plain = compile_with(&v, &reg, &CompileOptions::default()).unwrap();
    let mono = compile_with(&v, &reg, &CompileOptions { monotone: true }).unwrap();
    let (mut a, mut b) = ([0.0], [0.0]);
    let fields = Fields { time: 1.6 };
    plain.eval_nodes(&[0], &fields, &Stub, &mut a).unwrap();
    mono.eval_nodes(&[0], &fields, &Stub, &mut b).unwrap();
    assert_ne!(a, b);
    let r = elmer_expr::curve::cubic_spline(&x, &y, true);
    assert_eq!(b[0].to_bits(), elmer_expr::curve::interpolate_curve(&x, &y, 1.6, Some(&r)).to_bits());
}

#[test]
fn unsorted_scalar_table_is_rejected() {
    let v = table(&["time"], &[0.0, 2.0, 1.0], &[0.0, 1.0, 2.0], false);
    assert!(matches!(
        compile_with(&v, &UdfRegistry::new(), &CompileOptions::default()),
        Err(ExprError::Invalid(_))
    ));
}

#[test]
fn multi_variable_table_uses_first_variable() {
    let c = compile(&table(&["pressure, temperature"], &[0.0, 4.0], &[0.0, 8.0], false));
    // first dependency is pressure = 2 everywhere
    assert_eq!(eval(&c, &[0, 5], 0.0), vec![4.0, 4.0]);
}

#[test]
fn coordinate_dependency() {
    let c = compile(&table(&["coordinate 1"], &[0.0, 10.0], &[0.0, 100.0], false));
    assert_eq!(eval(&c, &[1, 2], 0.0), vec![10.0, 20.0]);
    let m = compile(&dep(&["coordinate"], DependentKind::Matc("tx(0)+tx(1)".into()), 1, 1));
    // x + y = node + 0.5
    assert_eq!(eval(&m, &[2], 0.0), vec![2.5]);
}

#[test]
fn tensor_table_is_point_major_row_major() {
    // 2x2 tensor at two table points
    let y = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0];
    let v = dep(&["time"], DependentKind::Table { x: vec![0.0, 1.0], y: y.to_vec(), cubic: true }, 2, 2);
    let c = compile(&v);
    assert_eq!(eval(&c, &[0], 0.5), vec![5.5, 11.0, 16.5, 22.0]);
}

#[test]
fn matc_variable_constant_and_tensor() {
    let c = compile(&dep(&["temperature"], DependentKind::Matc("tx(0)*3".into()), 1, 1));
    assert_eq!(eval(&c, &[1, 2], 0.0), vec![30.0, 60.0]);
    // time-only MATC binds `st`
    let c = compile(&dep(&[], DependentKind::Matc("st(0)*2".into()), 1, 1));
    assert_eq!(eval(&c, &[0, 1], 1.5), vec![3.0, 3.0]);
    let c = compile(&dep(&["temperature"], DependentKind::Matc("tensor".into()), 2, 2));
    assert_eq!(eval(&c, &[1], 0.0), vec![10.0, 20.0, 3.0, 4.0]);
    let bad = compile(&dep(&["temperature"], DependentKind::Matc("nope".into()), 1, 1));
    let mut o = [0.0];
    assert!(matches!(
        bad.eval_nodes(&[1], &Fields { time: 0.0 }, &Stub, &mut o),
        Err(ExprError::Matc(_))
    ));
    let k = compile(&dep(&[], DependentKind::Matc("const".into()), 1, 1));
    assert_eq!(k.eval_const(&Stub, [0.0; 3]).unwrap(), 42.0);
}

#[test]
fn procedures_resolve_through_registry_like_load_c() {
    let reg = UdfRegistry::new();
    reg.register_fortran(
        "libmat.so",
        "MyCond",
        Udf::Node(Arc::new(|node, t| t[0] + node as f64)),
    );
    reg.register("libmat", "CConst", Udf::Const(Arc::new(|p| p[0] + 100.0)));
    reg.register_fortran(
        "libmat",
        "Tens",
        Udf::Array(Arc::new(|_, t, o| o.iter_mut().for_each(|v| *v = t[0]))),
    );
    let opts = CompileOptions::default();
    let spelling = |lib: &str, entry: &str| {
        dep(
            &["temperature"],
            DependentKind::Procedure { library: lib.into(), entry: entry.into() },
            1,
            1,
        )
    };
    // any case, with or without .so, resolves to the mangled Fortran symbol
    for (lib, entry) in [("libmat", "mycond"), ("libmat.so", "MYCOND"), (" libmat ", "MyCond")] {
        let c = compile_with(&spelling(lib, entry), &reg, &opts).unwrap();
        let mut out = [0.0; 2];
        c.eval_nodes(&[1, 2], &Fields { time: 0.0 }, &Stub, &mut out).unwrap();
        assert_eq!(out, [11.0, 22.0]);
    }
    // an unmangled C symbol resolves only by its exact spelling
    let none = compile_with(
        &dep(&[], DependentKind::Procedure { library: "libmat".into(), entry: "cconst".into() }, 1, 1),
        &reg,
        &opts,
    );
    assert!(matches!(none, Err(ExprError::UnknownProcedure { .. })));
    let c = compile_with(
        &dep(&[], DependentKind::Procedure { library: "libmat.so".into(), entry: "CConst".into() }, 1, 1),
        &reg,
        &opts,
    )
    .unwrap();
    assert_eq!(eval_with(&c, &[3]), vec![103.0]);
    let t = compile_with(
        &dep(&["temperature"], DependentKind::Procedure { library: "libmat".into(), entry: "tens".into() }, 2, 1),
        &reg,
        &opts,
    )
    .unwrap();
    assert_eq!(eval_with(&t, &[1]), vec![10.0, 10.0]);
    // scalar function used for a tensor property
    let bad = compile_with(
        &dep(&["temperature"], DependentKind::Procedure { library: "libmat".into(), entry: "mycond".into() }, 2, 2),
        &reg,
        &opts,
    );
    assert!(matches!(bad, Err(ExprError::Invalid(_))));
}

fn eval_with(c: &elmer_expr::Compiled, nodes: &[u32]) -> Vec<f64> {
    eval(c, nodes, 0.0)
}

#[test]
fn global_registry_serves_plain_compile() {
    UdfRegistry::global().register_fortran(
        "libglobal",
        "Twice",
        Udf::Node(Arc::new(|_, t| 2.0 * t[0])),
    );
    let c = compile(&dep(
        &["temperature"],
        DependentKind::Procedure { library: "libglobal".into(), entry: "twice".into() },
        1,
        1,
    ));
    assert_eq!(eval(&c, &[2], 0.0), vec![40.0]);
    let missing = compile(&dep(
        &[],
        DependentKind::Procedure { library: "libnone".into(), entry: "f".into() },
        1,
        1,
    ));
    let mut o = [0.0];
    assert!(matches!(
        missing.eval_nodes(&[0], &Fields { time: 0.0 }, &Stub, &mut o),
        Err(ExprError::UnknownProcedure { .. })
    ));
}

#[test]
fn registry_name_helpers() {
    assert_eq!(fortran_mangle("MyFunc"), "myfunc_");
    assert_eq!(normalize_library(" libx.so "), "libx");
    assert_eq!(normalize_library("libx"), "libx");
    assert_eq!(split_procedure("libx.so  MyFunc"), Some(("libx.so", "MyFunc")));
    assert_eq!(split_procedure("onlylib"), None);
}

#[test]
fn integration_point_evaluation_interpolates_dependencies_first() {
    // T at nodes 0,1 = 0,10; nonlinear table on T so interpolate-then-apply != apply-then-interpolate
    let c = compile(&table(&["temperature"], &[0.0, 5.0, 10.0], &[0.0, 5.0, 100.0], false));
    let basis = [0.5, 0.5, 0.25, 0.75];
    let mut out = [0.0; 2];
    c.eval_ips(&[0, 1], &basis, &Fields { time: 0.0 }, &Stub, &mut out).unwrap();
    assert_eq!(out[0], 5.0);
    assert_eq!(out[1], 5.0 + 2.5 / 5.0 * 95.0);
    // a missing nodal value zeroes the point
    let mut out = [1.0; 1];
    c.eval_ips(&[2, 3], &[0.5, 0.5], &Fields { time: 0.0 }, &Stub, &mut out).unwrap();
    assert_eq!(out, [0.0]);
}

#[test]
fn shape_errors() {
    let c = compile(&Value::real(1.0));
    let mut out = [0.0; 3];
    assert!(matches!(
        c.eval_nodes(&[0, 1], &Fields { time: 0.0 }, &Stub, &mut out),
        Err(ExprError::Shape(_))
    ));
    assert_eq!(MISSING, f64::MAX);
}
