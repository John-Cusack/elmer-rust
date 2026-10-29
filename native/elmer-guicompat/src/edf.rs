//! Reader for ElmerGUI's EDF dialect files (`ElmerGUI/Application/edf/*.xml`).
//!
//! Only the parts the SIF generator consumes are kept: per section, each
//! `Parameter`'s widget, name, `SifName`, `Type` and `DefaultValue`.

use std::collections::BTreeMap;

use quick_xml::events::Event;
use quick_xml::Reader;

#[derive(Debug, thiserror::Error)]
pub enum EdfError {
    #[error("EDF XML: {0}")]
    Xml(String),
    #[error("EDF has no PDE element")]
    NoPde,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Param {
    pub widget: String,
    pub name: String,
    pub sif_name: Option<String>,
    pub ty: Option<String>,
    pub default: Option<String>,
}

impl Param {
    /// The name `SifGenerator::handleLineEdit` writes: `SifName`, else `Name`.
    pub fn key(&self) -> &str {
        self.sif_name.as_deref().unwrap_or(&self.name)
    }
}

#[derive(Clone, Debug, Default)]
pub struct Edf {
    /// `PDE Name="..."` attribute.
    pub pde: String,
    /// The `<Name>` child of `PDE`: what the GUI writes as `Equation = "..."`
    /// in each solver block.
    pub title: String,
    /// Section element name (`Equation`, `Solver`, `Material`, `BodyForce`,
    /// `InitialCondition`, `BoundaryCondition`) to its parameters.
    pub sections: BTreeMap<String, Vec<Param>>,
}

impl Edf {
    pub fn section(&self, name: &str) -> &[Param] {
        self.sections.get(name).map_or(&[], Vec::as_slice)
    }

    pub fn param(&self, section: &str, key: &str) -> Option<&Param> {
        let want = elmer_core::normalize_key(key);
        self.section(section)
            .iter()
            .find(|p| p.widget != "Label" && elmer_core::normalize_key(p.key()) == want)
    }

    /// `DefaultValue` of a Solver-section parameter (`Procedure`, `Variable`).
    pub fn solver_default(&self, key: &str) -> Option<&str> {
        self.param("Solver", key).and_then(|p| p.default.as_deref())
    }
}

pub fn parse_edf(xml: &str) -> Result<Edf, EdfError> {
    let mut reader = Reader::from_str(xml);
    let mut edf = Edf::default();
    let mut stack: Vec<String> = Vec::new();
    let mut current: Option<Param> = None;
    let mut seen_pde = false;
    loop {
        let ev = reader
            .read_event()
            .map_err(|e| EdfError::Xml(e.to_string()))?;
        match ev {
            Event::Start(e) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).into_owned();
                if name == "PDE" && stack.len() == 1 {
                    seen_pde = true;
                    for a in e.attributes().flatten() {
                        if a.key.as_ref() == b"Name" {
                            edf.pde = a
                                .unescape_value()
                                .map_err(|e| EdfError::Xml(e.to_string()))?
                                .trim()
                                .to_string();
                        }
                    }
                }
                if name == "Parameter" && stack.len() == 3 {
                    let mut p = Param::default();
                    for a in e.attributes().flatten() {
                        if a.key.as_ref() == b"Widget" {
                            p.widget = a
                                .unescape_value()
                                .map_err(|e| EdfError::Xml(e.to_string()))?
                                .into_owned();
                        }
                    }
                    current = Some(p);
                }
                stack.push(name);
            }
            Event::End(_) => {
                let name = stack.pop().unwrap_or_default();
                if name == "Parameter" && stack.len() == 3 {
                    if let Some(p) = current.take() {
                        edf.sections.entry(stack[2].clone()).or_default().push(p);
                    }
                }
            }
            Event::Text(t) => {
                let text = t.decode().map_err(|e| EdfError::Xml(e.to_string()))?;
                let text = text.trim();
                if text.is_empty() {
                    continue;
                }
                if stack.len() == 3 && stack[2] == "Name" {
                    edf.title = text.to_string();
                } else if stack.len() == 5 {
                    if let Some(p) = current.as_mut() {
                        match stack[4].as_str() {
                            "Name" => p.name = text.to_string(),
                            "SifName" => p.sif_name = Some(text.to_string()),
                            "Type" => p.ty = Some(text.to_string()),
                            "DefaultValue" => p.default = Some(text.to_string()),
                            _ => {}
                        }
                    }
                }
            }
            Event::Eof => break,
            _ => {}
        }
    }
    if !seen_pde {
        return Err(EdfError::NoPde);
    }
    Ok(edf)
}

#[cfg(test)]
mod tests {
    use super::*;

    const XML: &str = r#"<?xml version='1.0'?>
<!DOCTYPE edf>
<edf version="1.0" >
  <PDE Name="Heat Equation" >
    <Name>Heat Equation</Name>
    <Solver>
      <Parameter Widget="Edit"><Name > Procedure </Name>
        <DefaultValue> "HeatSolve" "HeatSolver" </DefaultValue></Parameter>
      <Parameter Widget="Label"><Name>Options</Name></Parameter>
    </Solver>
    <Material>
      <Parameter Widget="Edit"><Name> Heat Conductivity </Name>
        <SifName> Heat conductivity </SifName><Type> Real </Type></Parameter>
    </Material>
  </PDE>
</edf>"#;

    #[test]
    fn reads_sections_names_and_defaults() {
        let edf = parse_edf(XML).unwrap();
        assert_eq!(edf.pde, "Heat Equation");
        assert_eq!(edf.title, "Heat Equation");
        assert_eq!(
            edf.solver_default("Procedure"),
            Some("\"HeatSolve\" \"HeatSolver\"")
        );
        let p = edf.param("Material", "heat conductivity").unwrap();
        assert_eq!(p.key(), "Heat conductivity");
        assert_eq!(p.ty.as_deref(), Some("Real"));
        assert!(
            edf.param("Solver", "Options").is_none(),
            "labels are not parameters"
        );
    }

    #[test]
    fn rejects_files_without_a_pde() {
        assert!(matches!(
            parse_edf("<edf version=\"1.0\"></edf>"),
            Err(EdfError::NoPde)
        ));
    }
}
