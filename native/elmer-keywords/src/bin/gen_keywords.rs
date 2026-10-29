//! `gen-keywords`: render the keyword catalog JSON.
//!
//! ```text
//! gen-keywords [--source FILE] [--revision REV] [--out FILE | --check FILE]
//! ```
//!
//! Without `--source` the pinned snapshot is used. `--check FILE` exits 1 when
//! FILE differs from the freshly generated JSON.

use std::process::ExitCode;

use elmer_keywords::{Catalog, SNAPSHOT, SNAPSHOT_PATH, SNAPSHOT_REVISION};

fn main() -> ExitCode {
    let mut args = std::env::args().skip(1);
    let (mut source, mut revision, mut out, mut check) = (None, None, None, None);
    while let Some(a) = args.next() {
        let slot = match a.as_str() {
            "--source" => &mut source,
            "--revision" => &mut revision,
            "--out" => &mut out,
            "--check" => &mut check,
            _ => {
                eprintln!("usage: gen-keywords [--source FILE] [--revision REV] [--out FILE | --check FILE]");
                return ExitCode::from(2);
            }
        };
        *slot = args.next();
        if slot.is_none() {
            eprintln!("{a} needs a value");
            return ExitCode::from(2);
        }
    }
    let (text, path, rev) = match source {
        Some(p) => match std::fs::read_to_string(&p) {
            Ok(t) => (t, p, revision.unwrap_or_else(|| "unknown".into())),
            Err(e) => {
                eprintln!("{p}: {e}");
                return ExitCode::from(2);
            }
        },
        None => (
            SNAPSHOT.to_string(),
            SNAPSHOT_PATH.to_string(),
            revision.unwrap_or_else(|| SNAPSHOT_REVISION.to_string()),
        ),
    };
    let json = match Catalog::from_source(&text, &path, &rev) {
        Ok(c) => c.to_json(),
        Err(e) => {
            eprintln!("generation failed: {e}");
            return ExitCode::from(1);
        }
    };
    if let Some(p) = check {
        return match std::fs::read_to_string(&p) {
            Ok(have) if have == json => ExitCode::SUCCESS,
            Ok(_) => {
                eprintln!("{p} is stale; regenerate with --out");
                ExitCode::from(1)
            }
            Err(e) => {
                eprintln!("{p}: {e}");
                ExitCode::from(2)
            }
        };
    }
    match out {
        Some(p) => match std::fs::write(&p, json) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{p}: {e}");
                ExitCode::from(2)
            }
        },
        None => {
            print!("{json}");
            ExitCode::SUCCESS
        }
    }
}
