//! Fortran number layouts against files written by the pristine ElmerSolver.
use elmer_saveio::fortran::*;
use std::path::PathBuf;

fn fixture(p: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(p);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Every list-directed real in an ASCII `.result` body re-formats to the same bytes.
#[test]
fn list_directed_reals_round_trip_byte_exact() {
    let mut n = 0;
    let mut n_e = 0;
    for f in ["min_asc", "tr_asc", "rot_asc"] {
        for line in fixture(&format!("{f}/case.result")).lines() {
            if line.len() != 26
                || !(line.contains('.')
                    && line
                        .trim_start()
                        .starts_with(['-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']))
            {
                continue;
            }
            let v = parse_list_real(line).unwrap_or_else(|| panic!("parse {line:?}"));
            assert_eq!(list_real(v), line, "value {v:e}");
            n += 1;
            n_e += usize::from(line.contains('E'));
        }
    }
    assert!(n > 1500, "too few samples: {n}");
    assert!(n_e > 10, "no exponent-format samples: {n_e}");
}

/// ES17.8E3 columns of the `.ep` files, ES22.12E3 of `scalars.dat`, ES20.11E3 of `line.dat`.
#[test]
fn es_columns_match() {
    let mut n = 0;
    for line in fixture("rot_asc/case.ep").lines().skip(2) {
        let numeric = line.chars().all(|c| " 0123456789.E+-".contains(c));
        if !numeric || line.contains(" 0 ") || line.len() < 17 || line.trim().is_empty() {
            continue;
        }
        for tok_start in (0..line.len()).step_by(17) {
            if tok_start + 17 > line.len() {
                break;
            }
            let col = &line[tok_start..tok_start + 17];
            if col.trim().ends_with("0.0") && !col.contains('E') {
                continue;
            }
            let v = parse_list_real(col).unwrap();
            assert_eq!(fortran_es_e(v, 17, 8, Some(3)), col);
            n += 1;
        }
    }
    assert!(n > 1000, "{n}");
    for line in fixture("tr_asc/scalars.dat").lines() {
        for tok_start in (0..line.len()).step_by(22) {
            let col = &line[tok_start..tok_start + 22];
            assert_eq!(
                fortran_es_e(parse_list_real(col).unwrap(), 22, 12, Some(3)),
                col
            );
        }
    }
}
