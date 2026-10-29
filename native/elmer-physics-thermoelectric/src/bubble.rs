//! The `Element = "p:1 b:N"` element: bilinear nodal functions plus `N`
//! interior bubbles of the p-hierarchic basis, which Elmer keeps as global
//! dofs and this crate condenses element by element (the same linear
//! algebra, see [`crate::condense`]).

use elmer_core::{BasisEval, ElementBasis, ElementType};
use elmer_fe_hp::maps::bubble_dofs;
use elmer_fe_hp::{PElementBasis, PElementSpec};

use crate::error::ThermoElectricError;

/// The element definition of `Element = "p:1 b:N"`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BubbleConfig {
    /// `b:N`: bubbles per element.
    pub nb: usize,
    /// `Serendipity P elements` of the Simulation section (default true).
    pub serendipity: bool,
}

/// Parses an `Element` string into the number of bubbles per element.
/// `p:1` (also `n:1`) and `b:N` are ported; every other definition is refused.
pub fn parse_element(s: &str) -> Result<usize, ThermoElectricError> {
    let mut nb = 0;
    for tok in s.split_whitespace() {
        let tok = tok.to_ascii_lowercase();
        let Some((key, val)) = tok.split_once(':') else {
            return Err(ThermoElectricError::Unsupported(format!("Element = \"{s}\"")));
        };
        let Ok(val) = val.parse::<usize>() else {
            return Err(ThermoElectricError::Unsupported(format!("Element = \"{s}\"")));
        };
        match (key, val) {
            ("p" | "n", 1) => {}
            ("b", n) => nb = n,
            _ => {
                return Err(ThermoElectricError::Unsupported(format!(
                    "Element = \"{s}\" (only the linear \"p:1\" / \"n:1\" element with \"b:N\" bubbles)"
                )))
            }
        }
    }
    Ok(nb)
}

impl BubbleConfig {
    /// `getEffectiveBubbleP(Element, 1, nb)`: the lowest degree whose bubble
    /// count reaches `nb`.
    pub fn effective_p(&self, ecode: u16) -> i32 {
        let mut p = 1;
        while bubble_dofs(ecode, p, self.serendipity) < self.nb && p < 64 {
            p += 1;
        }
        p
    }

    /// The Gauss point count Elmer's `getNumberOfGaussPoints` gives the
    /// quadrilateral: `((max bubble degree per direction) + 1) ** 2`.
    pub fn gauss_points(&self, ecode: u16) -> usize {
        let p = self.effective_p(ecode);
        let deg = if self.serendipity { p - 2 } else { p }.max(0);
        let m = (deg.max(1) + 1) as usize;
        m * m
    }

    /// The basis of an element of type `etype`: nodal functions plus the
    /// first `nb` bubbles (`ElementInfo` skips what does not fit `nd`).
    pub fn basis(&self, etype: ElementType) -> Result<BubbleBasis, ThermoElectricError> {
        if etype.0 != 404 {
            return Err(ThermoElectricError::Unsupported(format!(
                "bubbles on element type {} (only quadrilaterals)",
                etype.0
            )));
        }
        let mut spec = PElementSpec::uniform(etype.0, 1, None, None).map_err(ThermoElectricError::Basis)?;
        spec.bubble_p = self.effective_p(etype.0);
        spec.serendipity = self.serendipity;
        Ok(BubbleBasis {
            inner: PElementBasis::new(spec),
            keep: 4 + self.nb,
        })
    }
}

/// The nodal-plus-bubble basis, cut to the element's `nd` functions.
pub struct BubbleBasis {
    inner: PElementBasis,
    keep: usize,
}

impl ElementBasis for BubbleBasis {
    fn eval(&self, etype: ElementType, nodal_xyz: &[[f64; 3]], uvw: [f64; 3], out: &mut BasisEval) -> Result<(), String> {
        self.inner.eval(etype, nodal_xyz, uvw, out)?;
        if out.basis.len() < self.keep {
            return Err(format!("{} basis functions, {} needed", out.basis.len(), self.keep));
        }
        out.basis.truncate(self.keep);
        out.dbasis_dx.truncate(self.keep);
        Ok(())
    }
}
