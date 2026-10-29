#![allow(clippy::type_complexity)]
//! Writer selection rules, big-endian files and error paths.
mod support;
use elmer_core::field::{FieldKind, FieldRef, NO_DOF};
use elmer_core::FieldSource;
use elmer_saveio::result::*;
use elmer_saveio::writer::*;
use elmer_saveio::SaveIoError;
use support::*;

struct Src {
    fields: Vec<(String, FieldKind, usize, Option<Vec<u32>>, Vec<f64>)>,
    time: f64,
    step: usize,
}

impl FieldSource for Src {
    fn fields(&self) -> Vec<FieldRef<'_>> {
        self.fields
            .iter()
            .map(|(n, k, d, p, v)| FieldRef {
                name: n,
                kind: *k,
                dofs: *d,
                perm: p.as_deref(),
                values: v,
            })
            .collect()
    }
    fn time(&self) -> f64 {
        self.time
    }
    fn step(&self) -> usize {
        self.step
    }
}

fn src() -> Src {
    let perm = Some(vec![1, NO_DOF, 0]);
    Src {
        fields: vec![
            (
                "temperature".into(),
                FieldKind::Nodal,
                1,
                perm.clone(),
                vec![1.5, 2.5],
            ),
            (
                "velocity".into(),
                FieldKind::Nodal,
                2,
                perm.clone(),
                vec![1.0, 2.0, 3.0, 4.0],
            ),
            (
                "coordinate 1".into(),
                FieldKind::Nodal,
                1,
                None,
                vec![0.0, 1.0, 0.0],
            ),
            ("time".into(), FieldKind::Global, 1, None, vec![0.25]),
            ("gauss".into(), FieldKind::Ip, 1, None, vec![9.0; 4]),
        ],
        time: 0.25,
        step: 7,
    }
}

fn write(opts: ResultOptions, format: ResultFormat, saves: usize) -> ResultFile {
    let mesh = TestMesh::load("tr_asc");
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = N.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    let p = std::env::temp_dir().join(format!(
        "elmer-saveio-opt-{}-{n}.result",
        std::process::id()
    ));
    let mut w = ResultWriter::create(&p, format, opts).unwrap();
    let s = src();
    for _ in 0..saves {
        w.save(&mesh, &s).unwrap();
    }
    w.finish().unwrap();
    ResultFile::from_bytes(&std::fs::read(&p).unwrap()).unwrap()
}

#[test]
fn default_selection_and_layout() {
    let rf = write(
        ResultOptions {
            started_at: Some("t".into()),
            ..Default::default()
        },
        ResultFormat::Ascii,
        2,
    );
    let names: Vec<_> = rf
        .header
        .vars
        .iter()
        .map(|v| (v.name.as_str(), v.dofs))
        .collect();
    // vector first, then scalars with the vector's components after `temperature`;
    // coordinates, globals and IP fields are off by default.
    assert_eq!(
        names,
        [
            ("velocity", 2),
            ("temperature", 1),
            ("velocity 1", 1),
            ("velocity 2", 1)
        ]
    );
    assert_eq!(rf.steps.len(), 2);
    assert_eq!(
        (
            rf.steps[1].saved_count,
            rf.steps[1].timestep,
            rf.steps[1].time
        ),
        (2, 7, 0.25)
    );
    let t = &rf.steps[0];
    assert!(matches!(t.vars[0].perm, PermRecord::Table { size: 3, .. }));
    assert_eq!(t.vars[1].perm, PermRecord::UsePrevious);
    assert_eq!(t.vars[0].values, vec![2.5, 1.5]); // nodes 1 (slot 1) and 3 (slot 0)
    assert_eq!(t.vars[1].values, vec![3.0, 1.0]); // velocity component 0 at slots 1, 0
    assert_eq!(t.vars[2].values, vec![4.0, 2.0]);
    let f = rf.restart_field(0, "velocity 2").unwrap().unwrap();
    assert_eq!((f.node, f.slot), (vec![0, 2], vec![1, 0]));
}

#[test]
fn optional_variables() {
    let o = ResultOptions {
        save_global: true,
        save_coordinates: true,
        save_ip: true,
        ..Default::default()
    };
    let rf = write(o, ResultFormat::BINARY, 1);
    let names: Vec<_> = rf.header.vars.iter().map(|v| v.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "velocity",
            "temperature",
            "velocity 1",
            "velocity 2",
            "coordinate 1",
            "time",
            "gauss"
        ]
    );
    let st = &rf.steps[0];
    assert_eq!(st.vars[3].perm, PermRecord::Null);
    assert_eq!(st.vars[4].values, vec![0.25]);
    assert_eq!(st.vars[5].values.len(), 4);
}

#[test]
fn output_variable_list_selects_vector_components() {
    let o = ResultOptions {
        output_variables: vec!["velocity".into()],
        ..Default::default()
    };
    let rf = write(o, ResultFormat::Ascii, 1);
    let names: Vec<_> = rf.header.vars.iter().map(|v| v.name.as_str()).collect();
    assert_eq!(names, ["velocity", "velocity 1", "velocity 2"]);
}

#[test]
fn big_endian_round_trips() {
    let rf = write(
        ResultOptions {
            started_at: Some("t".into()),
            ..Default::default()
        },
        ResultFormat::Binary(Endian::Big),
        2,
    );
    assert_eq!(rf.header.format, ResultFormat::Binary(Endian::Big));
    let (main, pos) = rf.to_bytes().unwrap();
    assert!(main.starts_with(b" BINARY 3.B\n"));
    let again = ResultFile::from_bytes(&main).unwrap();
    assert_eq!(again, rf);
    let pos = PosFile::from_bytes(&pos.unwrap()).unwrap();
    assert_eq!(pos.endian, Endian::Big);
    assert_eq!(pos.offsets.len(), 2);
}

#[test]
fn convenience_writer_and_errors() {
    let mesh = TestMesh::load("tr_asc");
    let p = std::env::temp_dir().join(format!("elmer-saveio-conv-{}.result", std::process::id()));
    write_result(&p, &mesh, &src(), ResultFormat::Ascii).unwrap();
    assert!(std::fs::read_to_string(&p)
        .unwrap()
        .starts_with(" ASCII 3\n!File started at: "));

    // slot out of range
    let bad = Src {
        fields: vec![("a".into(), FieldKind::Nodal, 1, Some(vec![5]), vec![1.0])],
        time: 0.0,
        step: 0,
    };
    assert!(matches!(
        write_result(&p, &mesh, &bad, ResultFormat::Ascii),
        Err(SaveIoError::Input(_))
    ));
    // value count not divisible by dofs
    let bad = Src {
        fields: vec![("a".into(), FieldKind::Nodal, 2, None, vec![1.0; 3])],
        time: 0.0,
        step: 0,
    };
    assert!(matches!(
        write_result(&p, &mesh, &bad, ResultFormat::Ascii),
        Err(SaveIoError::Input(_))
    ));
    // IP field with a permutation
    let o = ResultOptions {
        save_ip: true,
        ..Default::default()
    };
    let bad = Src {
        fields: vec![("a".into(), FieldKind::Ip, 1, Some(vec![0]), vec![1.0])],
        time: 0.0,
        step: 0,
    };
    assert!(matches!(
        write_result_with(&p, &mesh, &bad, ResultFormat::Ascii, o),
        Err(SaveIoError::Unsupported(_))
    ));
    // the variable set may not change between saves
    let mut w = ResultWriter::create(&p, ResultFormat::Ascii, ResultOptions::default()).unwrap();
    w.save(&mesh, &src()).unwrap();
    let mut other = src();
    other.fields.remove(0);
    assert!(matches!(w.save(&mesh, &other), Err(SaveIoError::Input(_))));
}

#[test]
fn use_previous_needs_a_table() {
    let rf = ResultFile {
        header: ResultHeader {
            format: ResultFormat::Ascii,
            started_at: "t".into(),
            vars: vec![VarHeader {
                name: "a".into(),
                field_size: 1,
                perm_size: 1,
                dofs: 1,
                equation: "e".into(),
            }],
            n_nodes: 1,
        },
        steps: vec![StepRecord {
            saved_count: 1,
            timestep: 1,
            time: 0.0,
            vars: vec![VarRecord {
                name: "a".into(),
                perm: PermRecord::UsePrevious,
                values: vec![1.0],
            }],
        }],
    };
    assert!(rf.to_bytes().is_err());
}
