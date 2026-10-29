//! `.sif` file parser: header + section dispatch (`LoadInputFile`) and
//! entry parsing (`SectionContents`) from `ModelDescription.F90`.
//!
//! Deviations from Elmer (all documented here):
//! - Single pass instead of ScanOnly + load passes. Observable behaviour
//!   is identical except for `Boundary N` counts (Elmer pre-counts them
//!   in the scan pass; here sections grow dynamically) and pathological
//!   MATC-`urand`-dependent section counts.
//! - `Procedure` targets are stored, not `dlopen`ed (later lane).
//! - Include nesting is not limited to Fortran unit numbers.
//! - `RUN` blocks are kept as separate `Model`s; `Sif::merged` replays
//!   them the way repeated `ReloadInputFile` accumulates.

use crate::error::SifError;
use crate::keywords::{KwType, lookup};
use crate::lex::{Delim, Lexer, Token};
use crate::model::*;

pub struct ParseResult {
    pub sif: Sif,
    pub warnings: Vec<String>,
}

pub fn parse_file(path: &str) -> Result<ParseResult, SifError> {
    parse_file_inner(path, true)
}

fn parse_file_inner(path: &str, do_eval: bool) -> Result<ParseResult, SifError> {
    let lexer = Lexer::open(path, do_eval)?;
    // Default `Check Keywords` level is Warn (ModelDescription.F90:366).
    let mut p = Parser { lexer, warnings: Vec::new(), numbering: true, check_abort: 1 };
    let mut sif = Sif::default();
    p.parse_header_only(&mut sif)?;
    p.parse_blocks(&mut sif)?;
    Ok(ParseResult { sif, warnings: p.warnings })
}

/// Parse a file included at section level (`LoadIncludeFile` with
/// `BaseLoad = .FALSE.`): sections only, appended to `model`.
fn parse_sections_file(path: &str, model: &mut Model, p: &mut Parser) -> Result<(), SifError> {
    let mut sub = Parser {
        lexer: Lexer::open(path, true)?,
        warnings: Vec::new(),
        numbering: p.numbering,
        check_abort: p.check_abort,
    };
    // Included files have no header pass; sections append directly.
    sub.parse_one_block(model)?;
    p.check_abort = sub.check_abort;
    Ok(())
}

struct Parser {
    lexer: Lexer,
    warnings: Vec<String>,
    numbering: bool,
    check_abort: i32,
}

impl Parser {
    fn err(&self, tok: &Token, msg: impl Into<String>) -> SifError {
        let mut e = SifError::new(tok.line, msg);
        e.context.push(tok.file.clone());
        e
    }

    // ---- header (BaseLoad) ----
    fn parse_header_only(&mut self, sif: &mut Sif) -> Result<(), SifError> {
        let mut header = Header::default();
        while let Some(tok) = self.lexer.next_token()? {
            let text = tok.text.trim_end().to_string();
            if text.is_empty() {
                continue;
            }
            if text == "end" {
                break;
            }
            if text.starts_with("check keywords") {
                let v = text.strip_prefix("check keywords").unwrap_or("").trim();
                let v = v.strip_prefix('"').unwrap_or(v);
                match v {
                    "ignore" => self.check_abort = 0,
                    "warn" => self.check_abort = 1,
                    "silent" => self.check_abort = 2,
                    "abort" => self.check_abort = 3,
                    _ => {}
                }
            } else if text == "echo on master" || text == "echo on" {
                header.echo = true;
            } else if text == "echo off" {
                header.echo = false;
            } else if text == "numbering on" {
                self.numbering = true;
                header.numbering = true;
            } else if text == "numbering off" {
                self.numbering = false;
                header.numbering = false;
            } else if text.starts_with("bodies")
                || text.starts_with("initial conditions")
                || text.starts_with("boundaries")
                || text.starts_with("boundary conditions")
                || text.starts_with("components")
                || text.starts_with("equations")
                || text.starts_with("solvers")
                || text.starts_with("materials")
                || text.starts_with("body forces")
            {
                // Documentation only; counts come from the sections.
            } else if text.starts_with("mesh db") {
                let rest = text.strip_prefix("mesh db").unwrap_or("").trim();
                let mut it = rest.splitn(2, [' ', '\t']);
                let dir = it.next().unwrap_or(".").to_string();
                let name = it.next().unwrap_or("").trim().to_string();
                if name.is_empty() {
                    header.mesh_db = Some((dir.clone(), ".".to_string()));
                } else {
                    header.mesh_db = Some((dir.clone(), format!("{dir}/{name}")));
                }
            } else if text == "header" {
                // explicit `Header` opener; fields follow until `end`.
            } else if text.starts_with("include path") {
                header.include_path = Some(text.strip_prefix("include path").unwrap_or("").trim().to_string());
            } else if text.starts_with("results directory") {
                header.results_dir =
                    Some(text.strip_prefix("results directory").unwrap_or("").trim().to_string());
            } else if text.starts_with("simulation id") {
                header.simulation_id =
                    Some(text.strip_prefix("simulation id").unwrap_or("").trim().to_string());
            } else {
                return Err(self.err(
                    &tok,
                    format!("unknown input field in header section: {text}"),
                ));
            }
        }
        header.present = true;
        header.check_keywords = match self.check_abort {
            1 => CheckLevel::Warn,
            2 => CheckLevel::Silent,
            3 => CheckLevel::Abort,
            _ => CheckLevel::Ignore,
        };
        let model = Model { header, ..Default::default() };
        sif.runs.push(model);
        Ok(())
    }

    // ---- RUN blocks ----
    fn parse_blocks(&mut self, sif: &mut Sif) -> Result<(), SifError> {
        loop {
            if sif.runs.is_empty() {
                sif.runs.push(Model::default());
            }
            let idx = sif.runs.len() - 1;
            let done = self.parse_one_block(&mut sif.runs[idx])?;
            if done {
                break;
            }
            // `RUN`: start a fresh block.
            sif.runs.push(Model::default());
        }
        Ok(())
    }

    /// Parse sections into `model` until `RUN`, `STOP`, or EOF.
    /// Returns true when the whole file is consumed.
    fn parse_one_block(&mut self, model: &mut Model) -> Result<bool, SifError> {
        let mut bc_count = 0u32;
        let mut ic_count = 0u32;
        let mut mat_count = 0u32;
        let mut bf_count = 0u32;
        let mut eq_count = 0u32;
        let mut body_count = 0u32;
        let mut solver_count = 0u32;
        let mut comp_count = 0u32;
        let mut simulation_seen = false;

        loop {
            let tok = match self.lexer.next_token()? {
                Some(t) => t,
                None => return Ok(true),
            };
            let text = tok.text.trim_end().to_string();
            if text.is_empty() {
                continue;
            }
            if text.starts_with("include") {
                let name = text["include".len().min(text.len())..].trim();
                let resolved = self.resolve_section_include(name, &tok)?;
                parse_sections_file(&resolved, model, self)?;
                continue;
            }
            if text == "header" {
                self.skip_to_end()?;
                continue;
            }
            if text.starts_with("echo ") || text.starts_with("check ") {
                continue;
            }
            if text == "run" {
                return Ok(false);
            }
            if text == "stop" {
                self.warnings.push(format!(
                    "{}:{}: STOP encountered, rest ignored",
                    tok.file, tok.line
                ));
                return Ok(true);
            }
            if text.starts_with("/*") {
                self.skip_to_close_comment()?;
                continue;
            }

            let free = self.check_abort <= 0;
            if text == "run control" {
                if simulation_seen {
                    return Err(self.err(&tok, "\"Run Control\" should precede \"Simulation\""));
                }
                let mut s = Section::default();
                self.parse_entries("run control", &mut s, free)?;
                model.control = Some(s);
            } else if text == "constants" {
                let mut s = Section::default();
                self.parse_entries("constants", &mut s, free)?;
                model.constants = Some(s);
            } else if text == "simulation" {
                simulation_seen = true;
                let mut s = Section::default();
                self.parse_entries("simulation", &mut s, free)?;
                model.simulation = Some(s);
            } else if text.starts_with("boundary condition") {
                let n = self.section_index(&text, "boundary condition", &tok, &mut bc_count)?;
                let s = get_section(&mut model.bcs, n);
                self.parse_entries("boundary condition", s, true)?;
            } else if text.starts_with("boundary") {
                let n: u32 = text.strip_prefix("boundary").unwrap_or("").split_whitespace().next()
                    .and_then(|w| w.parse().ok())
                    .ok_or_else(|| self.err(&tok, format!("problem reading section: {text}")))?;
                model.boundaries.push(NumberedSection { id: n, section: Section::default() });
                let last = model.boundaries.len() - 1;
                self.parse_entries("boundary", &mut model.boundaries[last].section, free)?;
            } else if text.starts_with("initial condition") {
                let n = self.section_index(&text, "initial condition", &tok, &mut ic_count)?;
                let s = get_section(&mut model.ics, n);
                self.parse_entries("initial condition", s, true)?;
            } else if text.starts_with("material") {
                let n = self.section_index(&text, "material", &tok, &mut mat_count)?;
                let s = get_section(&mut model.materials, n);
                self.parse_entries("material", s, free)?;
            } else if text.starts_with("body force") {
                let n = self.section_index(&text, "body force", &tok, &mut bf_count)?;
                let s = get_section(&mut model.body_forces, n);
                self.parse_entries("body force", s, free)?;
            } else if text.starts_with("equation") {
                let n = self.section_index(&text, "equation", &tok, &mut eq_count)?;
                let s = get_section(&mut model.equations, n);
                self.parse_entries("equation", s, free)?;
            } else if text.starts_with("body") {
                let n = self.section_index(&text, "body", &tok, &mut body_count)?;
                let s = get_section(&mut model.bodies, n);
                self.parse_entries("body", s, free)?;
            } else if text.starts_with("component") {
                let n = self.section_index(&text, "component", &tok, &mut comp_count)?;
                let s = get_section(&mut model.components, n);
                self.parse_entries("component", s, free)?;
            } else if text.starts_with("solver") {
                let n = self.section_index(&text, "solver", &tok, &mut solver_count)?;
                let s = get_section(&mut model.solvers, n);
                self.parse_entries("solver", s, free)?;
            } else {
                return Err(self.err(&tok, format!("unknown input section name: {text}")));
            }
        }
    }

    fn resolve_section_include(&self, name: &str, tok: &Token) -> Result<String, SifError> {
        for dir in self.lexer.include_path() {
            let d = dir.trim();
            if d.is_empty() {
                continue;
            }
            let cand = format!("{d}/{name}");
            if std::path::Path::new(&cand).is_file() {
                return Ok(cand);
            }
        }
        if std::path::Path::new(name).is_file() {
            return Ok(name.to_string());
        }
        Err(self.err(tok, format!("cannot find include file: {name}")))
    }

    fn section_index(
        &self,
        text: &str,
        head: &str,
        tok: &Token,
        counter: &mut u32,
    ) -> Result<u32, SifError> {
        let rest = text[head.len()..].trim();
        match rest.split_whitespace().next().and_then(|w| w.parse::<i32>().ok()) {
            Some(n) if n > 0 => {
                *counter = (*counter).max(n as u32);
                Ok(n as u32)
            }
            _ => {
                if self.numbering {
                    Err(self.err(tok, format!("problem reading section: {text}")))
                } else {
                    *counter += 1;
                    Ok(*counter)
                }
            }
        }
    }

    fn skip_to_end(&mut self) -> Result<(), SifError> {
        loop {
            match self.lexer.next_token()? {
                Some(t) if t.text.trim_end() == "end" => return Ok(()),
                Some(_) => {}
                None => return Ok(()),
            }
        }
    }

    fn skip_to_close_comment(&mut self) -> Result<(), SifError> {
        loop {
            match self.lexer.next_token()? {
                Some(t) if t.text.trim_end() == "*/" => return Ok(()),
                Some(_) => {}
                None => return Ok(()),
            }
        }
    }

    // ---- entries ----
    fn parse_entries(
        &mut self,
        section: &str,
        target: &mut Section,
        free_names: bool,
    ) -> Result<(), SifError> {
        // Keyword grouping blocks (`Linear System [ ... ]`), purely lexical.
        let mut prefix = String::new();
        let mut block_lens: Vec<usize> = Vec::new();
        while let Some(tok) = self.lexer.next_token()? {
            let mut text = tok.text.trim_end().to_string();
            if text.is_empty() {
                continue;
            }
            if text == "]" {
                let len = block_lens.pop().ok_or_else(|| {
                    self.err(&tok, format!("unmatched \"]\" in section: {section}"))
                })?;
                prefix.truncate(len);
                continue;
            }
            if text.contains('[') {
                let mut rest = text.clone();
                while let Some(b) = rest.find('[') {
                    let name = rest[..b].trim();
                    if name.is_empty() {
                        return Err(self.err(
                            &tok,
                            format!("keyword block \"[\" without a name in section: {section}"),
                        ));
                    }
                    if block_lens.len() >= 32 {
                        return Err(self.err(&tok, "keyword blocks nested too deep"));
                    }
                    block_lens.push(prefix.len());
                    prefix.push_str(name);
                    prefix.push(' ');
                    rest = rest[b + 1..].trim_start().to_string();
                }
                if rest.is_empty() {
                    continue;
                }
                text = rest;
            }
            if text.starts_with("end") {
                if !block_lens.is_empty() {
                    return Err(self.err(
                        &tok,
                        format!("unclosed keyword block \"{prefix}\" at end of section: {section}"),
                    ));
                }
                break;
            }
            if text.starts_with("include") {
                let name = text["include".len().min(text.len())..].trim();
                let resolved = self.resolve_section_include(name, &tok)?;
                let mut sub = Parser {
                    lexer: Lexer::open(&resolved, true)?,
                    warnings: Vec::new(),
                    numbering: self.numbering,
                    check_abort: self.check_abort,
                };
                let mut tmp = Section::default();
                sub.parse_entries(section, &mut tmp, free_names)?;
                for e in tmp.entries {
                    target.set(e);
                }
                self.warnings.append(&mut sub.warnings);
                continue;
            }
            if !prefix.is_empty() {
                text = format!("{prefix}{text}");
            }
            self.parse_entry(section, &text, &tok, target, free_names)?;
        }
        if !block_lens.is_empty() {
            return Err(SifError::new(0, format!("unclosed keyword block in section: {section}")));
        }
        Ok(())
    }

    fn check_keyword(
        &mut self,
        section: &str,
        name: &str,
        want: KwType,
        tok: &Token,
        free_names: bool,
    ) -> Result<(), SifError> {
        if let Some(have) = lookup(section, name) {
            let want_s = kw_name(want);
            let have_s = kw_name(have);
            if have_s == want_s || (have == KwType::String && want == KwType::File) {
                return Ok(());
            }
            return Err(self.err(
                tok,
                format!(
                    "keyword: [{name}] in section: [{section}] is given wrong type: \
                     [{want_s}], should be of type: [{have_s}]"
                ),
            ));
        }
        if self.check_abort <= 2 || free_names {
            if self.check_abort != 2 {
                self.warnings.push(format!(
                    "{}:{}: unlisted keyword: [{name}] in section: [{section}]",
                    tok.file, tok.line
                ));
            }
            return Ok(());
        }
        Err(self.err(tok, format!("unlisted keyword: [{name}] in section: [{section}]")))
    }

    fn parse_entry(
        &mut self,
        section: &str,
        name: &str,
        name_tok: &Token,
        target: &mut Section,
        free_names: bool,
    ) -> Result<(), SifError> {
        let mut n1: usize = 1;
        let mut n2: usize = 1;
        let mut size_given = false;
        let mut size_unknown = false;
        let mut partag: Option<i64> = None;
        let mut distribute = false;

        // Current value chunk being dispatched (GOTO 20 re-dispatches the
        // same chunk after a type prepend or -rpar/-distribute strip).
        let mut chunk = self.next_chunk(name, name_tok)?;
        let mut chunk_tok = name_tok.clone();

        loop {
            // Every chunk is literal-checked (ReadAndTrim reports
            // `string_literal` per call).
            if chunk.quoted {
                match self.db_type_opt(section, name) {
                    Some(db) if db == KwType::String || db == KwType::File => {
                        chunk = Token {
                            text: format!("{} {}", kw_name(db), chunk.text),
                            quoted: false,
                            delim: Delim::Eol,
                            file: chunk.file,
                            line: chunk.line,
                        };
                    }
                    _ => {
                        // Quoted value, unlisted or wrongly typed name:
                        // SyntaxError, then the cascade consumes chunks.
                        match self.cascade(section, name, name_tok, free_names)? {
                            Some((c, t)) => {
                                chunk = c;
                                chunk_tok = t;
                            }
                            None => return Ok(()),
                        }
                        continue;
                    }
                }
            }


            let text = chunk.text.trim_end().to_string();
            let mut words = text.split_whitespace();
            let keyword = words.next().unwrap_or("").to_string();
            let rest = text[keyword.len()..].trim_start().to_string();

            match keyword.as_str() {
                "real" => {
                    self.check_keyword(section, name, KwType::Real, &chunk_tok, free_names)?;
                    if rest_starts_with(&rest, "procedure") {
                        let (lib, entry) = split_proc(rest.strip_prefix("procedure").unwrap_or(""));
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::Procedure {
                                lib,
                                entry,
                                numeric: ProcNumeric::Real,
                            },
                            partag,
                            distribute,
                        });
                    } else if rest_starts_with(&rest, "matc") {
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::MatcExpr {
                                expr: rest.strip_prefix("matc").unwrap_or("").trim().to_string(),
                                numeric: ProcNumeric::Real,
                            },
                            partag,
                            distribute,
                        });
                    } else if size_given {
                        let vals = self.take_reals(n1 * n2, size_unknown, section, name, &chunk_tok, &rest)?;
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::RealArray { n1: vals.0, n2, values: vals.1 },
                            partag,
                            distribute,
                        });
                    } else {
                        let vals = self.take_reals(1, false, section, name, &chunk_tok, &rest)?;
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::RealArray { n1: 1, n2: 1, values: vals.1 },
                            partag,
                            distribute,
                        });
                    }
                    return Ok(());
                }
                "integer" => {
                    self.check_keyword(section, name, KwType::Integer, &chunk_tok, free_names)?;
                    if rest_starts_with(&rest, "procedure") {
                        let (lib, entry) = split_proc(rest.strip_prefix("procedure").unwrap_or(""));
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::Procedure {
                                lib,
                                entry,
                                numeric: ProcNumeric::Integer,
                            },
                            partag,
                            distribute,
                        });
                    } else if size_given {
                        let vals =
                            self.take_ints(n1, size_unknown, section, name, &chunk_tok, &rest)?;
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::IntArray { values: vals.1 },
                            partag,
                            distribute,
                        });
                    } else {
                        let vals = self.take_ints(1, false, section, name, &chunk_tok, &rest)?;
                        target.set(Entry {
                            name: name.to_string(),
                            value: Value::IntArray { values: vals.1 },
                            partag,
                            distribute,
                        });
                    }
                    return Ok(());
                }
                "logical" => {
                    self.check_keyword(section, name, KwType::Logical, &chunk_tok, free_names)?;
                    let v = rest_starts_with(&rest, "true") || rest.starts_with('1');
                    let f = rest_starts_with(&rest, "false") || rest.starts_with('0');
                    if !v && !f {
                        return Err(self.err(
                            &chunk_tok,
                            format!("problem reading logical keyword: {name}: {rest}"),
                        ));
                    }
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Logical(v),
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "string" => {
                    self.check_keyword(section, name, KwType::String, &chunk_tok, free_names)?;
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Text { text: rest, is_file: false },
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "file" => {
                    self.check_keyword(section, name, KwType::File, &chunk_tok, free_names)?;
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Text { text: rest, is_file: true },
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "variable" => {
                    let dep = rest;
                    let table = self.parse_table(section, name, &dep, n1, n2, size_given, size_unknown, &chunk_tok)?;
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Table(table),
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "equals" => {
                    let dep = rest;
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Table(Table {
                            dep,
                            kind: TableKind::Equals,
                            cubic: false,
                            monotone: false,
                            harmonic: false,
                            times: vec![1.0],
                            values: vec![1.0],
                            n1,
                            n2,
                            proc: None,
                            matc: None,
                        }),
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "opposes" => {
                    let dep = rest;
                    target.set(Entry {
                        name: name.to_string(),
                        value: Value::Table(Table {
                            dep,
                            kind: TableKind::Opposes,
                            cubic: false,
                            monotone: false,
                            harmonic: false,
                            times: vec![1.0],
                            values: vec![-1.0],
                            n1,
                            n2,
                            proc: None,
                            matc: None,
                        }),
                        partag,
                        distribute,
                    });
                    return Ok(());
                }
                "size" => {
                    if rest.trim_start().starts_with('*') {
                        size_unknown = true;
                        n1 = 100;
                        n2 = 1;
                    } else {
                        let mut nums = rest.split_whitespace();
                        n1 = nums.next().and_then(|w| w.parse().ok()).unwrap_or(1);
                        n2 = nums.next().and_then(|w| w.parse().ok()).unwrap_or(1);
                    }
                    size_given = true;
                    chunk = self.next_chunk(name, name_tok)?;
                    chunk_tok = chunk.clone();
                    continue;
                }
                "-rpar" => {
                    let mut it = rest.split_whitespace();
                    partag = it.next().and_then(|w| w.parse().ok());
                    let stripped = strip_first_word(&rest);
                    chunk = Token { text: stripped, ..chunk };
                    continue;
                }
                "-distribute" => {
                    distribute = true;
                    chunk = Token { text: rest, ..chunk };
                    continue;
                }
                _ => {
                    match self.db_type_opt(section, name) {
                        Some(db) => {
                            chunk = Token {
                                text: format!("{} {}", kw_name(db), text),
                                quoted: false,
                                delim: chunk.delim,
                                file: chunk.file,
                                line: chunk.line,
                            };
                            continue;
                        }
                        None => {
                            // Unlisted name without explicit type:
                            // SyntaxError, then cascade consumes chunks.
                            match self.cascade(section, name, name_tok, free_names)? {
                                Some((c, t)) => {
                                    chunk = c;
                                    chunk_tok = t;
                                }
                                None => return Ok(()),
                            }
                            continue;
                        }
                    }
                }
            }
        }
    }


    /// `lookup` without the SyntaxError mapping: `None` means unlisted.
    fn db_type_opt(&mut self, section: &str, name: &str) -> Option<KwType> {
        lookup(section, name)
    }

    fn warn_unlisted(&mut self, section: &str, name: &str, tok: &Token) {
        self.warnings.push(format!(
            "{}:{}: unlisted keyword: [{name}] in section: [{section}]",
            tok.file, tok.line
        ));
    }

    /// Non-fatal `SyntaxError` path for unlisted names: warn (unless
    /// silent) and continue the value loop with the next chunk, dropping
    /// the entry at EOF. Fatal when `Check Keywords Abort` applies.
    /// Returns the next chunk, or `None` at EOF (entry dropped).
    fn cascade(
        &mut self,
        section: &str,
        name: &str,
        name_tok: &Token,
        free_names: bool,
    ) -> Result<Option<(Token, Token)>, SifError> {
        if self.check_abort > 2 && !free_names {
            return Err(
                self.err(name_tok, format!("unlisted keyword: [{name}] in section: [{section}]"))
            );
        }
        if self.check_abort != 2 {
            // Elmer warns per consumed chunk; one warning per entry is
            // enough for triage and keeps proof output readable.
            self.warn_unlisted(section, name, name_tok);
        }
        match self.lexer.next_token()? {
            Some(c) => Ok(Some((c.clone(), c))),
            None => Ok(None),
        }
    }


    fn next_chunk(&mut self, name: &str, tok: &Token) -> Result<Token, SifError> {
        match self.lexer.next_token()? {
            Some(t) => Ok(t),
            None => Err(self.err(tok, format!("syntax error, unexpected end of file in: {name}"))),
        }
    }


    /// Collect `count` reals starting with `first` (rest of the current
    /// chunk), continuing across chunks. `unknown` (`(*)`) reads the first
    /// chunk only. Trailing words beyond `count` must be exactly `end`.
    fn take_reals(
        &mut self,
        count: usize,
        unknown: bool,
        section: &str,
        name: &str,
        tok: &Token,
        first: &str,
    ) -> Result<(usize, Vec<f64>), SifError> {
        let mut words: Vec<String> = first.split_whitespace().map(|s| s.to_string()).collect();
        if !unknown {
            while words.len() < count {
                match self.lexer.next_token()? {
                    Some(t) => {
                        words.extend(t.text.split_whitespace().map(|s| s.to_string()));
                    }
                    None => {
                        return Err(self.err(tok, format!("syntax error in {section}: {name}")));
                    }
                }
            }
        }
        let mut vals = Vec::new();
        let mut trailing: Vec<String> = Vec::new();
        for (i, w) in words.into_iter().enumerate() {
            if !unknown && i >= count {
                trailing.push(w);
                continue;
            }
            match parse_fortran_real(&w) {
                Some(v) => vals.push(v),
                None => {
                    return Err(self.err(
                        tok,
                        format!("problem reading real keyword \"{name}\": \"{w}\""),
                    ));
                }
            }
        }
        if unknown {
            return Ok((vals.len(), vals));
        }
        if vals.len() < count {
            return Err(self.err(tok, format!("syntax error in {section}: {name}")));
        }
        if !trailing.is_empty() {
            let rest = trailing.join(" ");
            if rest != "end" {
                return Err(self.err(
                    tok,
                    format!(
                        "mismatch of declared and given dimension for keyword \"{name}\". \
                         ignored input: {rest}"
                    ),
                ));
            }
        }
        Ok((count, vals))
    }

    fn take_ints(
        &mut self,
        count: usize,
        unknown: bool,
        section: &str,
        name: &str,
        tok: &Token,
        first: &str,
    ) -> Result<(usize, Vec<i32>), SifError> {
        let mut words: Vec<String> = first.split_whitespace().map(|s| s.to_string()).collect();
        if !unknown {
            while words.len() < count {
                match self.lexer.next_token()? {
                    Some(t) => {
                        words.extend(t.text.split_whitespace().map(|s| s.to_string()));
                    }
                    None => {
                        return Err(self.err(tok, format!("syntax error in {section}: {name}")));
                    }
                }
            }
        }
        let mut vals = Vec::new();
        let mut trailing: Vec<String> = Vec::new();
        for (i, w) in words.into_iter().enumerate() {
            if !unknown && i >= count {
                trailing.push(w);
                continue;
            }
            match w.parse::<i32>() {
                Ok(v) => vals.push(v),
                Err(_) => {
                    return Err(self.err(
                        tok,
                        format!("problem reading integer keyword \"{name}\": \"{w}\""),
                    ));
                }
            }
        }
        if unknown {
            return Ok((vals.len(), vals));
        }
        if vals.len() < count {
            return Err(self.err(tok, format!("syntax error in {section}: {name}")));
        }
        if !trailing.is_empty() {
            let rest = trailing.join(" ");
            if rest != "end" {
                return Err(self.err(
                    tok,
                    format!(
                        "mismatch between declared and given dimension for integer keyword \
                         \"{name}\". ignored input: {rest}"
                    ),
                ));
            }
        }
        Ok((count, vals))
    }


    /// `Variable dep` table body: an optional `Real [flags] [procedure|matc]`
    /// line, then rows until `end`. A `Real Procedure`/`Real MATC` line
    /// makes a procedure/MATC table with no rows (Elmer adds the 1-point
    /// entry immediately and exits the value loop).
    #[allow(clippy::too_many_arguments)]
    fn parse_table(
        &mut self,
        section: &str,
        name: &str,
        dep: &str,
        mut n1: usize,
        n2: usize,
        size_given: bool,
        size_unknown: bool,
        tok: &Token,
    ) -> Result<Table, SifError> {
        let (mut cubic, mut monotone, mut harmonic) = (false, false, false);
        let mut times = Vec::new();
        let mut values = Vec::new();
        // `false` until the row loop starts: Elmer swallows `integer` lines
        // and accepts the `real` header only before the first row.
        let mut in_rows = false;
        loop {
            let t = match self.lexer.next_token()? {
                Some(t) => t,
                None => {
                    return Err(self.err(tok, format!("table dependence has zero size: {name}")));
                }
            };
            let text = t.text.trim_end().to_string();
            if text.is_empty() {
                continue;
            }
            if text.starts_with("end") {
                break;
            }
            let mut row_text = text.clone();
            if !in_rows {
                if row_text.starts_with("integer")
                    && row_text.strip_prefix("integer").unwrap_or("").starts_with([' ', '\t', '\0'])
                {
                    // No `integer` branch for variable tables: ignored.
                    continue;
                }
                let mut rest = row_text.as_str();
                if let Some(after) = rest.strip_prefix("real") {
                    if after.is_empty()
                        || after.starts_with([' ', '\t', '\0'])
                    {
                        rest = after.trim_start();
                    }
                }
                if rest_starts_with(rest, "harmonic") {
                    harmonic = true;
                    rest = rest.strip_prefix("harmonic").unwrap_or("").trim_start();
                }
                let cubic_first = rest_starts_with(rest, "cubic");
                let monotone_first = rest_starts_with(rest, "monotone");
                if cubic_first {
                    cubic = true;
                    let r2 = rest.strip_prefix("cubic").unwrap_or("").trim_start();
                    if rest_starts_with(r2, "monotone") {
                        monotone = true;
                        rest = r2.strip_prefix("monotone").unwrap_or("").trim_start();
                    } else {
                        rest = r2;
                    }
                } else if monotone_first {
                    monotone = true;
                    let r2 = rest.strip_prefix("monotone").unwrap_or("").trim_start();
                    if rest_starts_with(r2, "cubic") {
                        cubic = true;
                        rest = r2.strip_prefix("cubic").unwrap_or("").trim_start();
                    } else {
                        rest = r2;
                        self.warnings.push(format!(
                            "{}:{}: monotone curves only applicable to cubic splines",
                            t.file, t.line
                        ));
                    }
                }
                // A bare `real`/flag line carries no row: keep reading.
                let consumed_header = row_text.starts_with("real")
                    || rest.len() != row_text.len();
                if rest_starts_with(rest, "procedure") {
                    let (lib, entry) = split_proc(rest.strip_prefix("procedure").unwrap_or(""));
                    return Ok(Table {
                        dep: dep.to_string(),
                        kind: TableKind::Variable,
                        cubic,
                        monotone,
                        harmonic,
                        times: Vec::new(),
                        values: Vec::new(),
                        n1,
                        n2,
                        proc: Some((lib, entry)),
                        matc: None,
                    });
                }
                if rest_starts_with(rest, "matc") {
                    let expr = rest.strip_prefix("matc").unwrap_or("").trim().to_string();
                    return Ok(Table {
                        dep: dep.to_string(),
                        kind: TableKind::Variable,
                        cubic,
                        monotone,
                        harmonic,
                        times: Vec::new(),
                        values: Vec::new(),
                        n1,
                        n2,
                        proc: None,
                        matc: Some(expr),
                    });
                }
                if rest.is_empty() && consumed_header {
                    continue;
                }
                row_text = rest.to_string();
            }
            in_rows = true;
            // Row: t + n1*n2 values, with cross-chunk continuation for
            // fixed sizes; `(*)` rows end at the chunk.
            let mut words: Vec<String> =
                row_text.split_whitespace().map(|s| s.to_string()).collect();
            if !size_unknown {
                while words.len() < 1 + n1 * n2 {
                    match self.lexer.next_token()? {
                        Some(nt) => {
                            words.extend(
                                nt.text.split_whitespace().map(|s| s.to_string()),
                            );
                        }
                        None => {
                            return Err(
                                self.err(tok, format!("table dependence has zero size: {name}"))
                            );
                        }
                    }
                }
            } else if words.is_empty() {
                continue;
            } else {
                n1 = words.len().saturating_sub(1);
            }
            if words.len() < 1 + n1 * n2 {
                return Err(self.err(tok, format!("table dependence has zero size: {name}")));
            }
            let tv = parse_fortran_real(&words[0]).ok_or_else(|| {
                self.err(tok, format!("problem reading real keyword: {name}"))
            })?;
            times.push(tv);
            for w in words.iter().take(1 + n1 * n2).skip(1) {
                values.push(parse_fortran_real(w).ok_or_else(|| {
                    self.err(tok, format!("problem reading real keyword \"{name}\": \"{w}\""))
                })?);
            }
            let _ = (section, size_given);
        }
        if times.is_empty() {
            return Err(self.err(tok, format!("table dependence has zero size: {name}")));
        }
        Ok(Table {
            dep: dep.to_string(),
            kind: TableKind::Variable,
            cubic,
            monotone,
            harmonic,
            times,
            values,
            n1,
            n2,
            proc: None,
            matc: None,
        })
    }
}

fn get_section(vec: &mut Vec<Section>, n: u32) -> &mut Section {
    while vec.len() < n as usize {
        vec.push(Section::default());
    }
    let s = &mut vec[n as usize - 1];
    s.id = n;
    s
}

fn kw_name(t: KwType) -> &'static str {
    match t {
        KwType::Real => "real",
        KwType::Integer => "integer",
        KwType::Logical => "logical",
        KwType::String => "string",
        KwType::File => "file",
    }
}

fn rest_starts_with(rest: &str, word: &str) -> bool {
    if rest.len() < word.len() {
        return false;
    }
    rest[..word.len()].eq_ignore_ascii_case(word)
        && (rest.len() == word.len()
            || rest[word.len()..].starts_with([' ', '\t']))
}

fn strip_first_word(s: &str) -> String {
    match s.find([' ', '\t']) {
        Some(i) => s[i + 1..].trim_start().to_string(),
        None => String::new(),
    }
}

/// `Procedure lib entry`: first two blank-separated words (GetProcAddr).
fn split_proc(s: &str) -> (String, String) {
    let mut it = s.split_whitespace();
    (it.next().unwrap_or("").to_string(), it.next().unwrap_or("").to_string())
}

/// Fortran list-directed real: Rust `f64` plus `d/D` exponents.
fn parse_fortran_real(w: &str) -> Option<f64> {
    if let Ok(v) = w.parse::<f64>() {
        return Some(v);
    }
    let mut chars: Vec<char> = w.chars().collect();
    let mut replaced = false;
    for c in chars.iter_mut() {
        if (*c == 'd' || *c == 'D') && !replaced {
            *c = 'e';
            replaced = true;
        }
    }
    if replaced {
        chars.into_iter().collect::<String>().parse().ok()
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_text(sif: &str) -> Result<ParseResult, SifError> {
        let mut p = std::env::temp_dir();
        p.push(format!("elmer_sif_test_{}.sif", std::process::id()));
        std::fs::write(&p, sif).unwrap();
        let r = parse_file(p.to_str().unwrap());
        let _ = std::fs::remove_file(&p);
        r
    }

    fn header() -> &'static str {
        "Header\nCheck Keywords Abort\nMesh DB \".\" \"m\"\nEnd\n"
    }

    #[test]
    fn namespaced_section_lines_accumulate() {
        let sif = format!(
            "{}Solver 1 :: Reference Norm = Real 1.0\nSolver 1 :: Linear System Solver = Iterative\n",
            header()
        );
        let r = parse_text(&sif).unwrap().sif;
        assert_eq!(r.runs.len(), 1);
        let s = &r.runs[0].solvers[0];
        assert_eq!(s.entries.len(), 2);
        assert!(matches!(
            &s.entries[0].value,
            Value::RealArray { values, .. } if values == &[1.0]
        ));
    }

    #[test]
    fn star_size_reads_rest_of_line() {
        let sif = format!(
            "{}Solver 1\nProcedure = \"HeatSolve\" \"HeatSolver\"\nMy Table(*) = Real 1.0 2.0 3.0\nEnd\n",
            header()
        );
        // `my table` is unlisted: Check Abort makes it fatal.
        assert!(parse_text(&sif).is_err());
    }

    #[test]
    fn procedure_table_has_no_rows() {
        let sif = format!(
            "{}Body Force 1\nSource = Variable Potential\nReal Procedure \"Source\" \"Source\"\nEnd\n",
            header()
        );
        let r = parse_text(&sif).unwrap().sif;
        let Value::Table(t) = &r.runs[0].body_forces[0].entries[0].value else {
            panic!("expected table");
        };
        assert_eq!(t.dep, "potential");
        assert_eq!(t.proc, Some(("Source".to_string(), "Source".to_string())));
        assert!(t.times.is_empty());
    }

    #[test]
    fn run_blocks_split_models() {
        let sif = format!(
            "{}Solver 1\nProcedure = \"a\" \"b\"\nEnd\nRUN\nSolver 1\nProcedure = \"c\" \"d\"\nEnd\n",
            header()
        );
        let r = parse_text(&sif).unwrap().sif;
        assert_eq!(r.runs.len(), 2);
        assert_eq!(r.merged().solvers.len(), 1);
    }

    #[test]
    fn unlisted_keyword_warn_cascades_without_abort() {
        // With Warn (not Abort), an unlisted untyped keyword is skipped
        // along with following untyped lines, like Elmer's SyntaxError
        // cascade; explicitly typed entries still land.
        let sif = "Header\nCheck Keywords Warn\nMesh DB \".\" \"m\"\nEnd\n\
            Simulation\nNoise = 1.5\nMore Noise\nReal 2.5\nEnd\n";
        // The unlisted `noise` cascades (SyntaxError) until the explicit
        // `Real` line, which completes the entry with 2.5, as in Elmer.
        let r = parse_text(sif).unwrap().sif;
        let sim = r.runs[0].simulation.as_ref().unwrap();
        assert!(matches!(
            sim.get("noise").map(|e| &e.value),
            Some(Value::RealArray { values, .. }) if values == &[2.5]
        ));
    }
}
