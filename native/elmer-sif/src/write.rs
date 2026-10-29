//! Canonical `.sif` serializer for the parsed model.
//!
//! Used by proof harnesses to materialize modified models (e.g. a case
//! with model-dependent transforms already applied and the corresponding
//! triggers removed). Output re-parses to an equal model; it is not
//! comment- or layout-preserving.

use crate::model::{Entry, Model, ProcNumeric, Section, Sif, TableKind, Value};

fn word_num(n: f64) -> String {
    // Shortest round-trip representation (Elmer reads list-directed).
    format!("{n:.17e}")
}

fn entry_lines(name: &str, e: &Entry, out: &mut String) {
    let mut head = String::new();
    if let Some(p) = e.partag {
        head.push_str(&format!("-rpar {p} "));
    }
    if e.distribute {
        head.push_str("-distribute ");
    }
    match &e.value {
        Value::RealArray { n1, n2, values } => {
            let size = if *n1 == 1 && *n2 == 1 {
                String::new()
            } else if *n2 == 1 {
                format!("({n1})")
            } else {
                format!("({n1},{n2})")
            };
            let vals: Vec<String> = values.iter().map(|v| word_num(*v)).collect();
            out.push_str(&format!("{name}{size} = {head}Real {}\n", vals.join(" ")));
        }
        Value::IntArray { values } => {
            let size = if values.len() == 1 {
                String::new()
            } else {
                format!("({})", values.len())
            };
            let vals: Vec<String> = values.iter().map(|v| v.to_string()).collect();
            out.push_str(&format!("{name}{size} = {head}Integer {}\n", vals.join(" ")));
        }
        Value::Logical(b) => {
            out.push_str(&format!("{name} = {head}Logical {}\n", if *b { "True" } else { "False" }));
        }
        Value::Text { text, is_file } => {
            let tag = if *is_file { "File" } else { "String" };
            out.push_str(&format!("{name} = {head}{tag} \"{text}\"\n"));
        }
        Value::Procedure { lib, entry, numeric } => {
            let tag = match numeric {
                ProcNumeric::Real => "Real",
                ProcNumeric::Integer => "Integer",
            };
            out.push_str(&format!("{name} = {head}{tag} Procedure \"{lib}\" \"{entry}\"\n"));
        }
        Value::MatcExpr { expr, numeric } => {
            let tag = match numeric {
                ProcNumeric::Real => "Real",
                ProcNumeric::Integer => "Integer",
            };
            out.push_str(&format!("{name} = {head}{tag} MATC \"{expr}\"\n"));
        }
        Value::Table(t) => {
            let size = if t.n1 == 1 && t.n2 == 1 {
                String::new()
            } else if t.n2 == 1 {
                format!("({})", t.n1)
            } else {
                format!("({},{})", t.n1, t.n2)
            };
            match t.kind {
                TableKind::Equals => {
                    out.push_str(&format!("{name}{size} = {head}Equals \"{}\"\n", t.dep));
                }
                TableKind::Opposes => {
                    out.push_str(&format!("{name}{size} = {head}Opposes \"{}\"\n", t.dep));
                }
                TableKind::Variable => {
                    out.push_str(&format!("{name}{size} = {head}Variable \"{}\"\n", t.dep));
                    if let Some((lib, entry)) = &t.proc {
                        out.push_str(&format!("Real Procedure \"{lib}\" \"{entry}\"\n"));
                        return;
                    }
                    if let Some(expr) = &t.matc {
                        out.push_str(&format!("Real MATC \"{expr}\"\n"));
                        return;
                    }
                    let mut flags = String::from("Real");
                    if t.harmonic {
                        flags.push_str(" Harmonic");
                    }
                    if t.cubic {
                        flags.push_str(" Cubic");
                        if t.monotone {
                            flags.push_str(" Monotone");
                        }
                    } else if t.monotone {
                        flags.push_str(" Monotone");
                    }
                    if flags != "Real" || true {
                        // Always emit the type line: unambiguous and
                        // re-parses identically.
                        out.push_str(&format!("{flags}\n"));
                    }
                    let per = t.n1 * t.n2;
                    for (r, tv) in t.times.iter().enumerate() {
                        let mut row = vec![word_num(*tv)];
                        for v in t.values.iter().skip(r * per).take(per) {
                            row.push(word_num(*v));
                        }
                        out.push_str(&row.join(" "));
                        out.push('\n');
                    }
                    out.push_str("End\n");
                }
            }
        }
    }
}

fn section_block(title: &str, id: Option<u32>, s: &Section, out: &mut String) {
    match id {
        Some(n) => out.push_str(&format!("{title} {n}\n")),
        None => out.push_str(&format!("{title}\n")),
    }
    for e in &s.entries {
        entry_lines(&e.name, e, out);
    }
    out.push_str("End\n");
}

/// Serialize one run block. `mesh_db` overrides the header mesh location.
pub fn write_model(model: &Model, mesh_db: Option<(&str, &str)>) -> String {
    let mut out = String::new();
    if !model.header.present {
        // Continuation RUN blocks carry no header; emit sections directly.
        write_sections(model, &mut out);
        return out;
    }
    out.push_str("Header\n");
    let level = match model.header.check_keywords {
        crate::model::CheckLevel::Ignore => "Ignore",
        crate::model::CheckLevel::Warn => "Warn",
        crate::model::CheckLevel::Silent => "Silent",
        crate::model::CheckLevel::Abort => "Abort",
    };
    out.push_str(&format!("Check Keywords {level}\n"));
    if model.header.echo {
        out.push_str("Echo on\n");
    }
    if !model.header.numbering {
        out.push_str("Numbering off\n");
    }
    let mdb = mesh_db.map(|(d, n)| (d.to_string(), n.to_string()));
    let mdb = mdb.or_else(|| model.header.mesh_db.clone());
    if let Some((mdir, mname)) = mdb {
        out.push_str(&format!("Mesh DB \"{mdir}\" \"{mname}\"\n"));
    }
    if let Some(p) = &model.header.include_path {
        out.push_str(&format!("Include Path {p}\n"));
    }
    if let Some(p) = &model.header.results_dir {
        out.push_str(&format!("Results Directory {p}\n"));
    }
    if let Some(p) = &model.header.simulation_id {
        out.push_str(&format!("Simulation Id {p}\n"));
    }
    out.push_str("End\n");
    write_sections(model, &mut out);
    out
}

fn write_sections(model: &Model, out: &mut String) {
    if let Some(s) = &model.control {
        section_block("Run Control", None, s, out);
    }
    if let Some(s) = &model.constants {
        section_block("Constants", None, s, out);
    }
    if let Some(s) = &model.simulation {
        section_block("Simulation", None, s, &mut *out);
    }
    for b in &model.bodies {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Body", Some(b.id), b, &mut *out);
    }
    for b in &model.boundaries {
        let mut titled = format!("Boundary {}", b.id);
        let _ = &mut titled;
        out.push_str(&format!("Boundary {}\n", b.id));
        for e in &b.section.entries {
            entry_lines(&e.name, e, &mut *out);
        }
        out.push_str("End\n");
    }
    for b in &model.bcs {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Boundary Condition", Some(b.id), b, &mut *out);
    }
    for b in &model.ics {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Initial Condition", Some(b.id), b, &mut *out);
    }
    for b in &model.materials {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Material", Some(b.id), b, &mut *out);
    }
    for b in &model.body_forces {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Body Force", Some(b.id), b, &mut *out);
    }
    for b in &model.equations {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Equation", Some(b.id), b, &mut *out);
    }
    for b in &model.solvers {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Solver", Some(b.id), b, &mut *out);
    }
    for b in &model.components {
        // Placeholders from sparse numbering (id 0, no entries) are not
        // emitted; re-parsing recreates identical slots.
        if b.id == 0 && b.entries.is_empty() {
            continue;
        }
        section_block("Component", Some(b.id), b, &mut *out);
    }
}

/// Serialize a whole input file (RUN blocks preserved).
pub fn write_sif(sif: &Sif, mesh_db: Option<(&str, &str)>) -> String {
    let mut out = String::new();
    for (i, run) in sif.runs.iter().enumerate() {
        if i > 0 {
            out.push_str("RUN\n");
        }
        out.push_str(&write_model(run, mesh_db));
    }
    out
}
