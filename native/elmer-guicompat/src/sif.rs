//! SIF text construction in ElmerGUI's style, and a scanner that recovers
//! `(section, key, value)` triples for type checking.

/// GUI `Variable` handling (`SifGenerator::parseSolverSpecificTab`): the EDF
/// default may use `Name(n)` or `Name[A:n B:m]` with `dim`/`cdim` as
/// expressions; `Name(n)` with `n > 1` is written `-dofs n Name`.
pub fn gui_variable(default: &str, cdim: i64) -> String {
    let var = default.split_whitespace().collect::<Vec<_>>().join(" ");
    let eval = |s: &str| -> i64 {
        let s = s.trim();
        if s == "cdim" || s == "dim" {
            cdim
        } else {
            s.parse().unwrap_or(1)
        }
    };
    if let Some((head, rest)) = var.split_once('[') {
        let inner = rest.split(']').next().unwrap_or("");
        let mut parts = inner.split(':');
        let mut sub = parts.next().unwrap_or("").trim().to_string();
        let mut out = format!("{}[", head.trim());
        for (i, part) in parts.enumerate() {
            let mut words = part.trim().split(' ');
            let dofs = eval(words.next().unwrap_or("1"));
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&format!("{sub}:{dofs}"));
            let tail: Vec<&str> = words.collect();
            if !tail.is_empty() {
                sub = tail.join(" ");
            }
        }
        out.push(']');
        return out;
    }
    if let Some((head, rest)) = var.split_once('(') {
        let dofs = eval(rest.split(')').next().unwrap_or("1"));
        let head = head.trim();
        return if dofs <= 1 {
            head.to_string()
        } else {
            format!("-dofs {dofs} {head}")
        };
    }
    var
}

/// Incremental SIF writer: `block("Solver 1")` then `kv(..)`.
#[derive(Default)]
pub struct SifWriter {
    out: String,
    open: bool,
}

impl SifWriter {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn block(&mut self, header: &str) -> &mut Self {
        self.end();
        self.out.push_str(header);
        self.out.push('\n');
        self.open = true;
        self
    }

    /// `  key = value`, skipped when the value is empty (the GUI's `addSifLine`).
    pub fn kv(&mut self, key: &str, value: &str) -> &mut Self {
        if !value.is_empty() {
            self.out.push_str(&format!("  {key} = {value}\n"));
        }
        self
    }

    /// A line without `=` (Header entries such as `Mesh DB "." "."`).
    pub fn raw(&mut self, line: &str) -> &mut Self {
        self.out.push_str(&format!("  {line}\n"));
        self
    }

    fn end(&mut self) {
        if self.open {
            self.out.push_str("End\n\n");
            self.open = false;
        }
    }

    pub fn finish(mut self) -> String {
        self.end();
        self.out
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assignment {
    /// Section header as written (`Solver 1`, `Boundary Condition 2`).
    pub section: String,
    /// Key with any `(n)` size suffix removed.
    pub key: String,
    pub value: String,
    pub line: usize,
}

/// Collect `key = value` lines per section. Lines inside `Header` without
/// `=` and bare `End` lines are not assignments.
pub fn scan(sif: &str) -> Vec<Assignment> {
    let mut out = Vec::new();
    let mut section = String::new();
    for (i, raw) in sif.lines().enumerate() {
        let line = raw.split('!').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line.eq_ignore_ascii_case("end") {
            section.clear();
            continue;
        }
        match line.split_once('=') {
            None => {
                if section.is_empty() && !raw.starts_with(char::is_whitespace) {
                    section = line.to_string();
                }
            }
            Some((k, v)) => {
                let key = k.split('(').next().unwrap_or(k).trim().to_string();
                out.push(Assignment {
                    section: section.clone(),
                    key,
                    value: v.trim().to_string(),
                    line: i + 1,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variable_forms_follow_the_gui() {
        assert_eq!(gui_variable("Temperature", 2), "Temperature");
        assert_eq!(gui_variable("Mesh Update(cdim)", 2), "-dofs 2 Mesh Update");
        assert_eq!(gui_variable("Mesh Update(cdim)", 1), "Mesh Update");
        assert_eq!(gui_variable("Pressure Wave(2)", 2), "-dofs 2 Pressure Wave");
        assert_eq!(
            gui_variable("Flow Solution[Velocity:cdim Pressure:1]", 2),
            "Flow Solution[Velocity:2 Pressure:1]"
        );
    }

    #[test]
    fn writer_and_scanner_round_trip() {
        let mut w = SifWriter::new();
        w.block("Header")
            .raw("Mesh DB \".\" \".\"")
            .block("Simulation")
            .kv("Timestep intervals(1)", "10")
            .kv("Empty", "")
            .block("Solver 1")
            .kv("Procedure", "\"HeatSolve\" \"HeatSolver\"");
        let text = w.finish();
        let a = scan(&text);
        assert_eq!(a.len(), 2);
        assert_eq!(
            (
                a[0].section.as_str(),
                a[0].key.as_str(),
                a[0].value.as_str()
            ),
            ("Simulation", "Timestep intervals", "10")
        );
        assert_eq!(a[1].section, "Solver 1");
    }
}
