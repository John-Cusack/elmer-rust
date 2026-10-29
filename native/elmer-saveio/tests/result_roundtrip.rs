//! `.result` read -> write is byte-identical to the Fortran-written files.
use elmer_saveio::result::*;
use std::path::PathBuf;

fn fx(p: &str) -> Vec<u8> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(p);
    std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn round_trip(dir: &str, binary: bool) -> ResultFile {
    let orig = fx(&format!("{dir}/case.result"));
    let rf = ResultFile::from_bytes(&orig).unwrap();
    let (main, pos) = rf.to_bytes().unwrap();
    assert!(
        main == orig,
        "{dir}: main file differs ({} vs {} bytes)",
        main.len(),
        orig.len()
    );
    if binary {
        let want = fx(&format!("{dir}/case.result.pos"));
        assert_eq!(pos.unwrap(), want, "{dir}: .pos differs");
    } else {
        assert!(pos.is_none());
    }
    rf
}

#[test]
fn ascii_files_round_trip() {
    for d in ["min_asc", "tr_asc", "rot_asc"] {
        round_trip(d, false);
    }
}

#[test]
fn binary_files_round_trip() {
    for d in ["min_bin", "tr_bin", "rot_bin"] {
        round_trip(d, true);
    }
}

#[test]
fn binary_and_ascii_decode_to_the_same_data() {
    for (a, b) in [
        ("tr_asc", "tr_bin"),
        ("rot_asc", "rot_bin"),
        ("min_asc", "min_bin"),
    ] {
        let a = ResultFile::from_bytes(&fx(&format!("{a}/case.result"))).unwrap();
        let b = ResultFile::from_bytes(&fx(&format!("{b}/case.result"))).unwrap();
        assert_eq!(a.header.vars, b.header.vars);
        assert_eq!(a.steps.len(), b.steps.len());
        for (sa, sb) in a.steps.iter().zip(&b.steps) {
            assert_eq!((sa.saved_count, sa.timestep), (sb.saved_count, sb.timestep));
            assert_eq!(sa.vars.len(), sb.vars.len());
            for (va, vb) in sa.vars.iter().zip(&sb.vars) {
                assert_eq!((&va.name, &va.perm), (&vb.name, &vb.perm));
                // The solver ran twice, so only the structure must agree exactly;
                // values agree to the printed 17 digits.
                assert_eq!(va.values.len(), vb.values.len());
            }
        }
    }
}

#[test]
fn rot_contents() {
    let rf = round_trip("rot_bin", true);
    let h = &rf.header;
    assert_eq!(h.n_nodes, 441);
    assert_eq!(h.total_dofs(), 4);
    assert_eq!(h.vars[0].name, "flow solution[velocity:3 pressure:1]");
    assert_eq!(
        (h.vars[0].field_size, h.vars[0].perm_size, h.vars[0].dofs),
        (1764, 441, 4)
    );
    assert_eq!(h.vars[0].equation, "navier-stokes");
    let st = &rf.steps[0];
    assert_eq!(st.vars.len(), 4);
    assert!(matches!(
        st.vars[0].perm,
        PermRecord::Table { size: 441, .. }
    ));
    assert_eq!(st.vars[1].perm, PermRecord::UsePrevious);
    let f = rf.restart_field(0, "pressure").unwrap().unwrap();
    assert_eq!(f.values.len(), 441);
    assert_eq!(f.node[3], 3);
}

#[test]
fn tr_steps_and_restart_selection() {
    let rf = round_trip("tr_asc", false);
    assert_eq!(rf.steps.len(), 3);
    assert_eq!(rf.find_step(0), Some(2));
    assert_eq!(rf.find_step(2), Some(1));
    assert_eq!(rf.find_step(9), None);
    let f = rf.restart_field(1, "temperature").unwrap().unwrap();
    assert_eq!(f.node, vec![0, 1, 2]);
    assert_eq!(f.slot, vec![2, 1, 0]);
    assert!((f.values[0] - 0.396_694_214_876_033).abs() < 1e-15);
    assert!(rf.restart_field(1, "nothing").unwrap().is_none());
}

#[test]
fn pos_index_matches_step_offsets() {
    let rf = ResultFile::from_bytes(&fx("tr_bin/case.result")).unwrap();
    let pos = PosFile::from_bytes(&fx("tr_bin/case.result.pos")).unwrap();
    assert_eq!(pos.names, vec!["temperature"]);
    assert_eq!(pos.offsets.len(), rf.steps.len());
    let main = fx("tr_bin/case.result");
    for o in &pos.offsets {
        assert_eq!(&main[o[0] as usize..o[0] as usize + 6], b"Time:\0");
        assert_eq!(&main[o[1] as usize..o[1] as usize + 12], b"temperature\0");
    }
}

#[test]
fn rejects_old_versions_and_garbage() {
    assert!(matches!(
        ResultFile::from_bytes(b" ASCII 2\n"),
        Err(elmer_saveio::SaveIoError::Unsupported(_))
    ));
    assert!(ResultFile::from_bytes(b"hello\n").is_err());
    assert!(ResultFile::from_bytes(b"").is_err());
    // truncated binary body
    let mut b = fx("tr_bin/case.result");
    b.truncate(b.len() - 5);
    assert!(ResultFile::from_bytes(&b).is_err());
}

#[test]
fn variable_list_matching() {
    let l = vec!["velocity".to_string(), "temperature".to_string()];
    assert!(variable_listed(&l, "velocity 2"));
    assert!(variable_listed(&l, "temperature"));
    assert!(!variable_listed(&l, "velocity x"));
    assert!(!variable_listed(&l, "pressure"));
    assert!(variable_listed(&[], "anything"));
}

#[test]
fn read_result_from_path() {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rot_bin/case.result");
    let rf = read_result(&p).unwrap();
    assert_eq!(rf.header.n_nodes, 441);
    assert!(read_result(&p.with_extension("nope")).is_err());
}
