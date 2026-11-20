//! The field-level check of a test whose reference norm passed: every
//! result file the pristine golden capture holds is graded against the
//! native run's with the oracle's Tier-A bars (`elmer-oracle grade-fields
//! --tolerances`, the same grading `grade-tree` does). A norm inside the
//! sif's `Reference Norm Tolerance` can hide fields that are off.
//!
//! The native run writes the files through the same `Output File`
//! injection as the oracle's O2 capture ([`output_file`]).

use std::path::{Path, PathBuf};
use std::process::Command;

use elmer_core::{Params, Value};
use elmer_solver::Model;

use crate::child::Verdict;

/// Where the golden captures, the Tier bars and the grader are.
#[derive(Clone, Debug)]
pub struct FieldCheck {
    /// `<golden-root>/<test>/run` holds each test's capture.
    pub golden: PathBuf,
    pub tolerances: PathBuf,
    /// The `elmer-oracle` binary.
    pub oracle: PathBuf,
}

/// The O2 injection (`elmer-oracle` `inject_sif` / `result_name_for`): a
/// Simulation without `Output File` gets `Output File = "case.result"`
/// (`<stem>.result` for a sif other than the startinfo one) and `Binary
/// Output = True`; one with it is left alone.
pub fn output_file(model: &mut Model, sif: &str, startinfo: Option<&str>) {
    if model.simulation.get("output file").is_some() {
        return;
    }
    let name = match startinfo {
        Some(s) if s != sif => format!("{}.result", sif.strip_suffix(".sif").unwrap_or(sif)),
        _ => "case.result".to_string(),
    };
    model.simulation.insert("Output File", Value::Str(name));
    model.simulation.insert("Binary Output", Value::Logical(true));
}

/// [`output_file`] on a staged sif's text, for a run whose front end the
/// harness does not own (`elmer-mpi-run`): the two keywords go right after
/// the `Simulation` header line. `None` when the Simulation section already
/// names an `Output File`, or the text has no Simulation section.
pub fn inject_text(text: &str, sif: &str, startinfo: Option<&str>) -> Option<String> {
    let bare = |l: &str| l.split('!').next().unwrap_or("").trim().to_ascii_lowercase();
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    let header = lines.iter().position(|l| bare(l) == "simulation")?;
    let names_output = lines[header + 1..]
        .iter()
        .map(|l| bare(l))
        .take_while(|l| l != "end")
        .any(|l| l.starts_with("output file") && l["output file".len()..].trim_start().starts_with(['=', '(']));
    if names_output {
        return None;
    }
    let mut m = Model::default();
    output_file(&mut m, sif, startinfo);
    let name = m.simulation.string("output file")?;
    let eol = if lines[header].ends_with('\n') { "" } else { "\n" };
    let mut out = String::with_capacity(text.len() + 64);
    for (i, l) in lines.iter().enumerate() {
        out.push_str(l);
        if i == header {
            out.push_str(&format!("{eol}  Output File = \"{name}\"\n  Binary Output = True\n"));
        }
    }
    Some(out)
}

/// [`inject_text`] on the sif file `path` in place. Sifs are bytes, not
/// UTF-8: each byte is one char, so the round trip is lossless.
pub fn inject_sif(path: &Path, startinfo: Option<&str>) -> std::io::Result<()> {
    let bytes = std::fs::read(path)?;
    let text: String = bytes.iter().map(|&b| char::from(b)).collect();
    let sif = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    match inject_text(&text, sif, startinfo) {
        Some(out) => std::fs::write(path, out.chars().map(|c| c as u32 as u8).collect::<Vec<u8>>()),
        None => Ok(()),
    }
}

/// A `.result` file by its first bytes (`elmer-oracle` `looks_like_result`).
fn looks_like_result(path: &Path) -> bool {
    let mut head = [0u8; 16];
    let n = std::fs::File::open(path).and_then(|mut f| std::io::Read::read(&mut f, &mut head)).unwrap_or(0);
    let h = String::from_utf8_lossy(&head[..n]);
    let h = h.trim_start();
    h.starts_with("BINARY 3.") || h.starts_with("ASCII 3")
}

/// Whether a result file holds any variable: something after its first
/// `Time:` block (binary) or line (ASCII).
fn has_fields(path: &Path) -> bool {
    let Ok(b) = std::fs::read(path) else { return false };
    let text_start = b.iter().position(|&c| c == b'\n').unwrap_or(0);
    if String::from_utf8_lossy(&b[..text_start]).trim() == "ASCII 3" {
        let text = String::from_utf8_lossy(&b);
        return text.lines().skip_while(|l| !l.starts_with("Time:")).nth(1).is_some();
    }
    let Some(nul) = b.iter().position(|&c| c == 0) else { return false };
    // "Time:\0", two int4 and one double.
    b.len() > nul + 1 + 6 + 16
}

/// The golden result files of a capture (`elmer-oracle`
/// `sigma::result_files`): files under `run` that look like results,
/// are not `.pos` and are not in the test's sources, and hold fields.
fn golden_files(run: &Path, source: &Path) -> Vec<String> {
    let mut out = Vec::new();
    let mut stack = vec![run.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else { continue };
        for e in entries.flatten() {
            let path = e.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            let Ok(rel) = path.strip_prefix(run) else { continue };
            let rel = rel.display().to_string();
            if rel.ends_with(".pos") || source.join(&rel).exists() || !looks_like_result(&path) || !has_fields(&path) {
                continue;
            }
            out.push(rel);
        }
    }
    out.sort();
    out
}

/// Grades a `match` of `test` staged in `dir`: `None` keeps the verdict;
/// otherwise the new verdict and its detail. A field over its bar is a
/// mismatch; a golden file the run did not write leaves the match
/// unverified.
pub fn check(fc: &FieldCheck, test: &str, dir: &Path, source: &Path) -> Option<(Verdict, String)> {
    let run = fc.golden.join(test).join("run");
    if !run.is_dir() {
        return Some((Verdict::Unverified, "fields: no golden capture".into()));
    }
    let mut missing = Vec::new();
    for rel in golden_files(&run, source) {
        let cand = dir.join(&rel);
        if !cand.is_file() {
            missing.push(rel);
            continue;
        }
        let out = Command::new(&fc.oracle)
            .arg("grade-fields")
            .arg(run.join(&rel))
            .arg(&cand)
            .arg("--tolerances")
            .arg(&fc.tolerances)
            .args(["--test", test, "--file", &rel])
            .output();
        let out = match out {
            Ok(o) => o,
            Err(e) => return Some((Verdict::Error, format!("fields: {}: {e}", fc.oracle.display()))),
        };
        let stdout = String::from_utf8_lossy(&out.stdout);
        match out.status.code() {
            Some(0) => {}
            Some(1) if stdout.contains("mismatches: ") => {
                let first = stdout.lines().find(|l| !l.starts_with("note:") && !l.starts_with("mismatches:")).unwrap_or("");
                return Some((Verdict::Mismatch, format!("fields: {rel}: {first}")));
            }
            _ => {
                let err = String::from_utf8_lossy(&out.stderr);
                return Some((Verdict::Error, format!("fields: {rel}: {}", err.lines().next().unwrap_or("grader failed"))));
            }
        }
    }
    if let Some(rel) = missing.first() {
        return Some((Verdict::Unverified, format!("fields: {rel}: no Output File written")));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_file_follows_the_o2_injection() {
        let mut m = Model::default();
        output_file(&mut m, "case.sif", Some("case.sif"));
        assert_eq!(m.simulation.string("output file"), Some("case.result"));
        assert_eq!(m.simulation.logical("binary output"), Some(true));
        let mut m = Model::default();
        output_file(&mut m, "taylorhood.sif", Some("case.sif"));
        assert_eq!(m.simulation.string("output file"), Some("taylorhood.result"));
        let mut m = Model::default();
        m.simulation.insert("Output File", Value::Str("own.dat".into()));
        output_file(&mut m, "case.sif", None);
        assert_eq!(m.simulation.string("output file"), Some("own.dat"));
        assert_eq!(m.simulation.get("binary output"), None);
    }

    #[test]
    fn inject_text_adds_the_o2_keywords_after_the_simulation_header() {
        let sif = "Header\n  Mesh DB \".\" \"angle\"\nEnd\n\nSimulation ! the run\n  Max Output Level = 5\nEnd\n";
        let out = inject_text(sif, "case.sif", Some("case.sif")).unwrap();
        assert_eq!(
            out,
            "Header\n  Mesh DB \".\" \"angle\"\nEnd\n\nSimulation ! the run\n  Output File = \"case.result\"\n  Binary Output = True\n  Max Output Level = 5\nEnd\n"
        );
        let out = inject_text("simulation\nend\n", "b.sif", Some("case.sif")).unwrap();
        assert!(out.contains("Output File = \"b.result\""));
        assert_eq!(inject_text("Simulation\n  Output File = \"own.dat\"\nEnd\n", "case.sif", None), None);
        assert_eq!(inject_text("Header\nEnd\n", "case.sif", None), None);
        // `Output File Final Only` is another keyword.
        assert!(inject_text("Simulation\n  Output File Final Only = True\nEnd\n", "case.sif", None).is_some());
    }
}
