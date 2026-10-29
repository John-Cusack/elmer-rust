//! SaveScalars / SaveLine files reproduced byte for byte from a captured run
//! (3 steps of transient heat on a one-triangle mesh).
mod support;
use elmer_saveio::savedata::*;
use support::*;

fn meta(date_from: &str) -> NamesMeta {
    let names = fx_str(date_from);
    let stamp = names
        .lines()
        .find_map(|l| l.strip_prefix("File started at: "))
        .unwrap()
        .to_string();
    NamesMeta {
        version: "26.2".into(),
        revision: Some("9f6af2f85".into()),
        compilation_date: Some("2026-09-23".into()),
        sif_name: Some("case.sif".into()),
        started_at: Some(stamp),
        ..Default::default()
    }
}

#[test]
fn scalars_dat_and_names_match() {
    let dir = std::env::temp_dir().join(format!("elmer-saveio-sc-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("scalars.dat");
    let mut w = ScalarsWriter::new(&path, "./scalars.dat", meta("tr_asc/scalars.dat.names"));
    let names: Vec<String> = ["max", "min", "mean", "sum", "int mean"]
        .iter()
        .map(|o| scalar_column_name(o, "temperature", &ScalarScope::Domain))
        .collect();
    let rows = fx_str("tr_asc/scalars.dat");
    for line in rows.lines() {
        let vals: Vec<(String, f64)> = names
            .iter()
            .enumerate()
            .map(|(i, n)| {
                (
                    n.clone(),
                    elmer_saveio::fortran::parse_list_real(&line[i * 22..(i + 1) * 22]).unwrap(),
                )
            })
            .collect();
        w.row(&vals).unwrap();
    }
    assert_eq!(std::fs::read_to_string(&path).unwrap(), rows);
    assert_eq!(
        std::fs::read_to_string(dir.join("scalars.dat.names")).unwrap(),
        fx_str("tr_asc/scalars.dat.names")
    );
}

#[test]
fn line_dat_and_names_match() {
    let want_names = fx_str("tr_asc/line.dat.names");
    let value_names: Vec<String> = [
        "coordinate 1",
        "coordinate 2",
        "coordinate 3",
        "temperature",
    ]
    .map(String::from)
    .to_vec();
    let got = line_names_text(
        &meta("tr_asc/line.dat.names"),
        "./line.dat",
        LineNodeCounts {
            bc: 2,
            ..Default::default()
        },
        LineColumns {
            timestep: true,
            ..Default::default()
        },
        &value_names,
    );
    assert_eq!(got, want_names);

    let rows = fx_str("tr_asc/line.dat");
    let mut out = String::new();
    for line in rows.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        let labels: Vec<i64> = t[..4].iter().map(|x| x.parse().unwrap()).collect();
        let vals: Vec<f64> = t[4..]
            .iter()
            .map(|x| elmer_saveio::fortran::parse_real_token(x).unwrap())
            .collect();
        out.push_str(&line_row_text(&labels, &vals));
    }
    assert_eq!(out, rows);
}

#[test]
fn scalars_options() {
    let m = NamesMeta {
        version: "9".into(),
        started_at: Some("t".into()),
        simulation_comment: Some("sim".into()),
        ..Default::default()
    };
    let t = scalars_names_text(&m, "./a.dat", &["x".into()], Some(3));
    assert!(t.starts_with(
        "sim\nMetadata for SaveScalars file: ./a.dat\nElmer version: 9\nFile started at: t\n \n"
    ));
    assert!(t.ends_with("   1: Line Marker\n   2: x\n"));
    assert_eq!(
        scalars_row_text(&[2.0, 0.5], &[true, false], Some(3)),
        "         3         2   5.000000000000E-001\n"
    );
    assert_eq!(
        scalar_column_name("boundary int", "heat flux", &ScalarScope::Bc(2)),
        "boundary int: heat flux over bc 2"
    );
    assert_eq!(
        scalar_column_name("mean", "t", &ScalarScope::Polyline(1)),
        "mean: t over polyline 1"
    );
    let cols = LineColumns {
        skip_boundary_info: true,
        ..Default::default()
    };
    assert_eq!(cols.n_labels(), 0);
    assert_eq!(
        LineColumns {
            line_marker: true,
            timestep: true,
            skip_boundary_info: false
        }
        .n_labels(),
        5
    );
}
