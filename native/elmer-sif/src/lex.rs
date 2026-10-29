//! Line lexer replicating Elmer's `ReadAndTrim` (`GeneralUtils.F90`).
//!
//! Structure mirrors the Fortran: `next` handles the SAVE-state dispatch
//! (fresh line vs pending remainder vs synthetic `end`), then ONE shared
//! processing path runs on the text (`!`-strip, `\`-continuation,
//! `#`/`$` expansion, token split). Pending remainders go through the
//! same path (Elmer reprocesses them, including continuations and MATC).
//!
//! Observable contract (verified against the oracle build, see
//! `proof/` outside this worktree):
//! - `!` ends the line, even inside quotes.
//! - trailing `\` (outside quotes) continues on the next physical line,
//!   including inside pending remainders; joined text is not rescanned
//!   for directives but is for comments/MATC.
//! - `include path ...` directives anywhere set the include search path.
//! - a line starting with `include ` splices the named file inline,
//!   resolved against each `;`-separated include-path entry, then CWD.
//!   The spliced file's first line skips directive detection, as in Elmer.
//! - `#` outside comments is a Lua expression (this build has no Lua, and
//!   neither does the campaign oracle: `WITH_LUA=FALSE`): it is an error,
//!   except as the last character of the line (matches `i < inlen`).
//! - `$...$` (or `$` to end of line) is evaluated by MATC at lex time,
//!   including inside quoted strings; output is not rescanned.
//! - statements split at `;`; `=` starts the value, `(` the size prefix,
//!   `::` a one-entry namespaced section line.
//! - names are lowercased outside quotes; quotes are removed; tabs become
//!   spaces.
//! - echo: Elmer prints the processed text on EVERY `ReadAndTrim` call
//!   (not once per line). `echo_stream` reproduces that stream for the
//!   `Echo on` differential.

use crate::error::SifError;
use crate::matc::matc_eval;
use std::collections::VecDeque;

/// One token from the lexer: either a section/entry name or a value chunk.
#[derive(Debug, Clone)]
pub struct Token {
    /// Lowercased outside quotes, quotes stripped.
    pub text: String,
    /// True when the token started with a `"` (string literal).
    pub quoted: bool,
    /// Delimiter that ended the token: `=`, `;`, `(`, `::`, or EOL.
    pub delim: Delim,
    /// Origin for diagnostics.
    pub file: String,
    /// 1-based line in `file` (first line of a continuation).
    pub line: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Delim {
    Equals,
    Semicolon,
    Paren,
    Namespace,
    Eol,
}

struct Source {
    path: String,
    lines: VecDeque<(u32, String)>,
}

impl Source {
    fn read(path: &str) -> Result<Self, String> {
        // Latin-1 decode: Elmer reads raw bytes and only folds ASCII case.
        // The corpus is ASCII except for latin-1 bytes inside `!` comments,
        // which the lexer strips before any comparison.
        let bytes = std::fs::read(path).map_err(|e| format!("cannot open file: {e}"))?;
        let text: String = bytes.iter().map(|&b| b as char).collect();
        let lines = text
            .lines()
            .enumerate()
            .map(|(i, l)| (i as u32 + 1, l.to_string()))
            .collect();
        Ok(Source { path: path.to_string(), lines })
    }
}

pub(crate) enum Pending {
    /// Raw remainder after `=` or `;` (unprocessed; shared path handles it).
    Value(String),
    /// Like Value, but the enclosing section closes after this entry.
    Namespaced(String),
    /// Raw remainder after `(` (unprocessed; `)`/`,` rewrite first).
    Size(String),
}

pub struct Lexer {
    stack: Vec<Source>,
    include_path: Vec<String>,
    pending: Option<Pending>,
    open_section: bool,
    do_eval: bool,
    /// The next fresh line is a spliced file's first line: skip
    /// directive detection once (ReadAndTrim reads it directly).
    fresh_splice: bool,
    /// File/line of the most recently read physical line (diagnostics).
    pub file: String,
    pub line: u32,
}

impl Lexer {
    pub fn open(path: &str, do_eval: bool) -> Result<Self, SifError> {
        let src = Source::read(path).map_err(|e| SifError::new(0, format!("{path}: {e}")))?;
        Ok(Lexer {
            stack: vec![src],
            include_path: Vec::new(),
            pending: None,
            open_section: false,
            do_eval,
            fresh_splice: false,
            file: path.to_string(),
            line: 0,
        })
    }

    /// Current `include path` directories (for section-level includes).
    pub fn include_path(&self) -> &[String] {
        &self.include_path
    }

    /// Proof/echo helper: every processed text Elmer would print with
    pub fn echo_stream(path: &str) -> Result<Vec<String>, SifError> {
        let mut lexer = Lexer::open(path, true)?;
        let mut out = Vec::new();
        while let Some((echo, _)) = lexer.next_with_echo()? {
            if let Some(e) = echo {
                out.push(e);
            }
        }
        Ok(out)
    }
    /// Backwards-compatible logical lines (fresh lines only, directives
    /// and splices applied, post-expansion). Used where per-call echo
    /// granularity is not needed.
    pub fn logical_lines(path: &str) -> Result<Vec<(String, String, u32)>, SifError> {
        let mut lexer = Lexer::open(path, true)?;
        let mut out = Vec::new();
        // Fresh-line logical lines: read physical lines with directive
        // handling, without token splitting.
        while let Some((raw, _spliced)) = lexer.fresh_text()? {
            let (file, no) = (lexer.file.clone(), lexer.line);
            let processed = lexer.process(raw, file.clone(), no)?;
            out.push((file, processed, no));
        }
        Ok(out)
    }

    fn physical_line(&mut self) -> Option<(String, u32)> {
        loop {
            let src = self.stack.last_mut()?;
            match src.lines.pop_front() {
                Some((n, l)) => {
                    self.file = src.path.clone();
                    self.line = n;
                    return Some((l, n));
                }
                None => {
                    self.stack.pop();
                }
            }
        }
    }

    fn resolve_include(&self, name: &str) -> Option<String> {
        // Elmer tries each include-path entry, then the bare name
        // (relative to CWD), and qualified names directly.
        for dir in &self.include_path {
            let d = dir.trim();
            if d.is_empty() {
                continue;
            }
            let cand = format!("{d}/{name}");
            if std::path::Path::new(&cand).is_file() {
                return Some(cand);
            }
        }
        if std::path::Path::new(name).is_file() {
            return Some(name.to_string());
        }
        None
    }

    /// Fresh text with directive handling: `include path` sets the path,
    /// `include ` splices (the spliced first line is returned with
    /// `spliced = true`). Returns None at end of input.
    fn fresh_text(&mut self) -> Result<Option<(String, bool)>, SifError> {
        loop {
            let (raw, lineno) = match self.physical_line() {
                Some(v) => v,
                None => return Ok(None),
            };
            let spliced = std::mem::replace(&mut self.fresh_splice, false);
            // Left-adjust (ADJUSTL); LEN_TRIM applies to fresh lines
            // (pending remainders keep their blanks, like Elmer).
            let adj = raw.trim().to_string();
            if !spliced {
                if let Some(rest) = ascii_prefix(&adj, "include path") {
                    self.include_path =
                        rest.trim().split(';').map(|s| s.to_string()).collect();
                    continue;
                }
                if let Some(prefix_name) = ascii_prefix(&adj, "include ") {
                    let mut name = prefix_name.to_string();
                    if self.do_eval && name.contains('$') {
                        name = expand_matc(&name).map_err(|e| {
                            SifError::new(lineno, format!("include file MATC: {e}"))
                        })?;
                    }
                    match self.resolve_include(name.trim()) {
                        Some(p) => {
                            let src = Source::read(&p).map_err(|e| {
                                SifError::new(lineno, format!("{p}: {e}"))
                            })?;
                            self.stack.push(src);
                            self.fresh_splice = true;
                        }
                        None => {
                            return Err(SifError::new(
                                lineno,
                                format!("cannot find include file: {}", name.trim()),
                            ));
                        }
                    }
                    continue;
                }
            }
            return Ok(Some((adj, spliced)));
        }
    }

    /// Shared processing path: `!`-strip, `\`-continuation, `#`/`$`
    /// expansion. Returns the processed text (what `Echo` prints).
    /// `file`/`lineno` track the first physical line for diagnostics.
    fn process(&mut self, raw: String, file: String, lineno: u32) -> Result<String, SifError> {
        let mut raw = raw;
        // No pre-trim: fresh lines arrive LEN_TRIMmed from `fresh_text`;
        // pending remainders keep their blanks. The `!`-strip below does
        // not re-trim, so its trailing blank survives into the echo.
        // `!` ends the line unconditionally (even inside quotes).
        if let Some(i) = raw.find('!') {
            raw.truncate(i);
        }
        // backslash is the last character (post-`!`-strip, untrimmed).
        // The `\` becomes a blank, the next physical line is appended
        // RAW, and the buffer end is re-trimmed (Elmer's j-scan).
        loop {
            let in_quotes = raw.bytes().filter(|&b| b == b'"').count() % 2 == 1;
            if !in_quotes && raw.ends_with('\\') {
                raw.pop();
                raw.push(' ');
                match self.physical_line() {
                    Some((next, _)) => raw.push_str(&next),
                    None => break,
                }
                let tl = raw.trim_end().len();
                raw.truncate(tl);
            } else {
                break;
            }
        }
        // Trimmed length for the `#`/`$` checks (`i < inlen`).
        let inlen = raw.trim_end().len();

        if self.do_eval {
            // Lua `#`: error unless last char (matches `i < inlen`).
            if let Some(i) = raw.find('#') {
                if i + 1 < inlen {
                    return Err(SifError::new(
                        lineno,
                        "LUA not included, cannot continue".to_string(),
                    ));
                }
            }
            // `$` strictly before the trimmed end triggers MATC on the
            // trimmed line; a trailing `$` is literal.
            if raw.find('$').map(|j| j + 1 < inlen).unwrap_or(false) {
                let tl = raw.trim_end().len();
                raw.truncate(tl);
                raw = expand_matc(&raw)
                    .map_err(|e| SifError::new(lineno, format!("MATC: {e}")))?;
            }
        }
        let _ = file;
        Ok(raw)
    }

    /// Next token plus the echo text Elmer prints for this call (`None`
    /// for the synthetic `end`, which returns before the echo point).
    pub fn next_with_echo(&mut self) -> Result<Option<(Option<String>, Token)>, SifError> {
        if self.open_section && self.pending.is_none() {
            self.open_section = false;
            return Ok(Some((
                None,
                Token {
                    text: "end ".to_string(),
                    quoted: false,
                    delim: Delim::Eol,
                    file: self.file.clone(),
                    line: self.line,
                },
            )));
        }
        let (text, file, line) = match self.pending.take() {
            Some(Pending::Value(rest)) => (rest, self.file.clone(), self.line),
            Some(Pending::Namespaced(rest)) => {
                self.open_section = true;
                (rest, self.file.clone(), self.line)
            }
            Some(Pending::Size(rest)) => {
                // Replicate the `Prefix == 'Size'` branch: `,` -> blank,
                // first `)` -> blank and end of the size prefix; the
                // remainder stays on the line for the next token.
                // echo keeps Elmer's `Size` capitalisation; the keyword
                // is lowercased again at token split.
                let combined = match rest.find(')') {
                    Some(r) => {
                        let head = rest[..r].replace(',', " ");
                        format!("Size {head} {}", &rest[r + 1..])
                    }
                    None => format!("Size {}", rest.replace(',', " ")),
                };
                (combined, self.file.clone(), self.line)
            }
            None => {
                let (raw, _spliced) = match self.fresh_text()? {
                    Some(v) => v,
                    None => return Ok(None),
                };
                (raw, self.file.clone(), self.line)
            }
        };
        let echo = self.process(text, file.clone(), line)?;
        // Elmer echoes only when the trimmed length is nonzero.
        let echo_text = echo.clone();
        let echo = if echo_text.trim().is_empty() { None } else { Some(echo_text) };
        let mut stash = None;
        let tok = split_token(&echo.clone().unwrap_or_default(), &file, line, &mut stash);
        self.pending = stash;
        Ok(Some((echo, tok)))
    }

    /// Next token. `None` at end of input (callers treat a section-level
    /// None as the implicit end of the final section, like Elmer).
    pub fn next_token(&mut self) -> Result<Option<Token>, SifError> {
        Ok(self.next_with_echo()?.map(|(_, t)| t))
    }
}

/// Byte-safe ASCII case-insensitive prefix match: returns the remainder
/// when `s` starts with `prefix` (ASCII `prefix`, arbitrary `s` bytes).
fn ascii_prefix<'a>(s: &'a str, prefix: &str) -> Option<&'a str> {
    let (sb, pb) = (s.as_bytes(), prefix.as_bytes());
    if sb.len() < pb.len() {
        return None;
    }
    for (a, b) in sb.iter().zip(pb.iter()) {
        if a.to_ascii_lowercase() != *b {
            return None;
        }
    }
    // The matched prefix is pure ASCII, so this boundary is valid.
    Some(&s[prefix.len()..])
}

/// Expand `$...$` pairs (or `$` to end of line) via MATC, left to right,
/// without rescanning substituted output.
fn expand_matc(line: &str) -> Result<String, String> {
    let bytes = line.as_bytes();
    let mut out = String::with_capacity(line.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'$' {
            out.push(bytes[i] as char);
            i += 1;
            continue;
        }
        let start = i + 1;
        let mut j = start;
        while j < bytes.len() && bytes[j] != b'$' {
            j += 1;
        }
        let cmd = &line[start..j];
        let eval = matc_eval(cmd)?;
        out.push_str(&eval);
        i = if j < bytes.len() { j + 1 } else { j };
    }
    Ok(out)
}

/// Split one blank-separated word-run off the front of `text`, honouring
/// quotes; lowercase outside quotes; record the delimiter (`=`, `;`, `(``,
/// `::`, EOL) and stash the remainder as the pending value for the next
/// call. Mirrors the two nested scan loops of `ReadAndTrim`.
pub(crate) fn split_token(
    text: &str,
    file: &str,
    line: u32,
    stash: &mut Option<Pending>,
) -> Token {
    let b = text.as_bytes();
    let mut i = 0;
    while i < b.len() && (b[i] == b' ' || b[i] == b'\t') {
        i += 1;
    }
    let quoted = i < b.len() && b[i] == b'"';
    let mut s = String::new();
    let mut in_quotes = false;
    let mut delim = Delim::Eol;
    // Outer loop: words joined by single blanks, like ReadAndTrim.
    loop {
        if i >= b.len() {
            break;
        }
        if b[i] == b'"' {
            in_quotes = !in_quotes;
            i += 1;
            if i >= b.len() {
                break;
            }
        }
        if !in_quotes {
            if b[i] == b'=' {
                delim = Delim::Equals;
                *stash = Some(Pending::Value(text[i + 1..].to_string()));
                break;
            }
            if b[i] == b';' {
                delim = Delim::Semicolon;
                *stash = Some(Pending::Value(text[i + 1..].to_string()));
                break;
            }
            if b[i] == b'(' {
                delim = Delim::Paren;
                *stash = Some(Pending::Size(text[i + 1..].to_string()));
                break;
            }
            if b[i] == b':' && i + 1 < b.len() && b[i + 1] == b':' {
                delim = Delim::Namespace;
                *stash = Some(Pending::Namespaced(text[i + 2..].to_string()));
                break;
            }
            if b[i] == b'!' || b[i] == b'#' || (b[i] < 32 && b[i] != b'\t') {
                break;
            }
        }
        // Inner loop: copy one word.
        while i < b.len() {
            let c = b[i];
            if c == b'"' {
                in_quotes = !in_quotes;
                i += 1;
                if i >= b.len() {
                    break;
                }
                continue;
            }
            if !in_quotes {
                if c == b' ' || c == b'=' || c == b';' || c == b'(' {
                    break;
                }
                if c == b':' && i + 1 < b.len() && b[i + 1] == b':' {
                    break;
                }
                if c == b'!' || c == b'#' || (c < 32 && c != b'\t') {
                    break;
                }
            }
            // Inside quotes: verbatim. Outside: lowercase, tab -> space.
            if in_quotes {
                s.push(c as char);
            } else if c == b'\t' {
                s.push(' ');
            } else if c.is_ascii_uppercase() {
                s.push((c + 32) as char);
            } else {
                s.push(c as char);
            }
            i += 1;
        }
        s.push(' ');
        while i < b.len() && (b[i] == b' ' || b[i] == b'\t') {
            i += 1;
        }
    }
    Token { text: s, quoted, delim, file: file.to_string(), line }
}
