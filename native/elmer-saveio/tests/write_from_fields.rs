//! Live fields -> `.result` / `.ep` reproduce the Fortran-written bytes.
mod support;
use elmer_saveio::ep::{write_ep, EpOptions};
use elmer_saveio::result::{ResultFile, ResultFormat};
use elmer_saveio::writer::{ResultOptions, ResultWriter};
use std::collections::HashMap;
use support::*;

fn tmp(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("elmer-saveio-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_all_steps(dir: &str, format: ResultFormat, out: &std::path::Path) {
    let orig = ResultFile::from_bytes(&fx(&format!("{dir}/case.result"))).unwrap();
    let mesh = TestMesh::load(dir);
    let equations: HashMap<String, String> = orig
        .header
        .vars
        .iter()
        .map(|v| (v.name.clone(), v.equation.clone()))
        .collect();
    let opts = ResultOptions {
        equations,
        started_at: Some(orig.header.started_at.clone()),
        ..Default::default()
    };
    let mut w = ResultWriter::create(out, format, opts).unwrap();
    for s in 0..orig.steps.len() {
        let src = FixtureSource::from_result(&orig, s);
        assert_eq!(w.save(&mesh, &src).unwrap(), s + 1);
    }
    w.finish().unwrap();
}

#[test]
fn result_writer_matches_fortran_bytes() {
    let d = tmp("res");
    for (dir, fmt) in [
        ("min_asc", ResultFormat::Ascii),
        ("tr_asc", ResultFormat::Ascii),
        ("rot_asc", ResultFormat::Ascii),
        ("min_bin", ResultFormat::BINARY),
        ("tr_bin", ResultFormat::BINARY),
        ("rot_bin", ResultFormat::BINARY),
    ] {
        let out = d.join(format!("{dir}.result"));
        write_all_steps(dir, fmt, &out);
        assert!(
            std::fs::read(&out).unwrap() == fx(&format!("{dir}/case.result")),
            "{dir}: .result differs"
        );
        if fmt != ResultFormat::Ascii {
            let pos = d.join(format!("{dir}.result.pos"));
            assert_eq!(
                std::fs::read(pos).unwrap(),
                fx(&format!("{dir}/case.result.pos")),
                "{dir}: .pos differs"
            );
        }
    }
}

fn ep_opts(dir: &str) -> EpOptions {
    let text = fx_str(&format!("{dir}/case.ep"));
    let stamp = text
        .lines()
        .nth(1)
        .unwrap()
        .strip_prefix("#File started at: ")
        .unwrap()
        .to_string();
    let mut o = EpOptions {
        started_at: Some(stamp),
        ..Default::default()
    };
    match dir {
        "min_asc" | "min_bin" => {
            o.boundary_constraint = [(2, 1), (3, 2)].into();
        }
        "tr_asc" | "tr_bin" => {
            o.body_names.insert(1, "plate".into());
            o.bc_names.insert(1, "fixed".into());
            o.bc_names.insert(2, "flux side".into());
            o.boundary_constraint = [(2, 1), (3, 2)].into();
        }
        "rot_asc" | "rot_bin" => {
            o.body_names.insert(1, "melt".into());
        }
        _ => {}
    }
    o
}

#[test]
fn ep_writer_matches_fortran_bytes() {
    let d = tmp("ep");
    for dir in ["min_asc", "tr_asc", "rot_asc"] {
        let orig = ResultFile::from_bytes(&fx(&format!("{dir}/case.result"))).unwrap();
        let mesh = TestMesh::load(dir);
        let sources: Vec<FixtureSource> = (0..orig.steps.len())
            .map(|s| FixtureSource::from_result(&orig, s))
            .collect();
        let refs: Vec<&dyn elmer_core::FieldSource> = sources.iter().map(|s| s as _).collect();
        let out = d.join(format!("{dir}.ep"));
        write_ep(&out, &mesh, &refs, ep_opts(dir)).unwrap();
        let got = std::fs::read_to_string(&out).unwrap();
        let want = fx_str(&format!("{dir}/case.ep"));
        if got != want {
            let (g, w): (Vec<_>, Vec<_>) = (got.lines().collect(), want.lines().collect());
            let i = (0..g.len().min(w.len()))
                .find(|&i| g[i] != w[i])
                .unwrap_or(g.len().min(w.len()));
            panic!(
                "{dir}: .ep differs at line {}:\n got  {:?}\n want {:?}",
                i + 1,
                g.get(i),
                w.get(i)
            );
        }
    }
}
