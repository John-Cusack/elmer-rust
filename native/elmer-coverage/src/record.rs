//! One row of the coverage results, and its `results.tsv` form.

use crate::child::Verdict;

#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub name: String,
    pub nprocs: u32,
    pub verdict: Verdict,
    /// Seconds spent on the native run (staging excluded).
    pub secs: f64,
    /// Solver library stems the case uses, lowercase.
    pub libs: Vec<String>,
    /// Every known blocker, `class: name`, in report order.
    pub blockers: Vec<String>,
    /// Failed comparison, fatal message or refusal text.
    pub detail: String,
}

impl Record {
    /// The first missing feature of a test the native path does not run.
    pub fn first(&self) -> Option<&str> {
        self.blockers.first().map(String::as_str)
    }
}

pub const HEADER: &str = "test\tnprocs\tverdict\tsecs\tlibs\tblockers\tdetail";

fn clean(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

pub fn to_tsv(records: &[Record]) -> String {
    let mut out = format!("{HEADER}\n");
    for r in records {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\n",
            clean(&r.name),
            r.nprocs,
            r.verdict.name(),
            r.secs,
            r.libs.join(","),
            clean(&r.blockers.join(" || ")),
            clean(&r.detail)
        ));
    }
    out
}

pub fn from_tsv(text: &str) -> Result<Vec<Record>, String> {
    let mut lines = text.lines();
    if lines.next() != Some(HEADER) {
        return Err("results.tsv: unexpected header".into());
    }
    lines
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let [name, nprocs, verdict, secs, libs, blockers, detail] = f[..] else {
                return Err(format!("results.tsv: bad row `{l}`"));
            };
            Ok(Record {
                name: name.to_string(),
                nprocs: nprocs.parse().map_err(|_| format!("bad nprocs in `{l}`"))?,
                verdict: Verdict::parse(verdict).ok_or_else(|| format!("bad verdict in `{l}`"))?,
                secs: secs.parse().map_err(|_| format!("bad secs in `{l}`"))?,
                libs: libs.split(',').filter(|s| !s.is_empty()).map(str::to_string).collect(),
                blockers: blockers.split(" || ").filter(|s| !s.is_empty()).map(str::to_string).collect(),
                detail: detail.to_string(),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_round_trip() {
        let rs = vec![
            Record {
                name: "A".into(),
                nprocs: 1,
                verdict: Verdict::Refused,
                secs: 0.5,
                libs: vec!["heatsolve".into(), "savedata".into()],
                blockers: vec!["solver: savedata savescalars".into(), "feature: eigenvalue analysis".into()],
                detail: "Solver 2: x".into(),
            },
            Record { name: "B".into(), nprocs: 2, verdict: Verdict::Match, secs: 1.25, libs: vec![], blockers: vec![], detail: String::new() },
        ];
        assert_eq!(from_tsv(&to_tsv(&rs)).unwrap(), rs);
    }
}
