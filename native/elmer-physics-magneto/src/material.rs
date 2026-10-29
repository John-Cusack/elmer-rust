//! Material and source coefficients of `WhitneyAVSolver`, read the way
//! `MagnetoDynamicsUtils` does. Only constant keyword values are taken;
//! anything evaluated per node or point is refused for now.

use elmer_core::{ParamMap, Params, Value};

use crate::error::MagnetoError;

/// Elmer's `PI` (Types.F90).
pub const PI: f64 = std::f64::consts::PI;

/// A constant real array: `(rows, cols, row-major data)`.
pub type ConstArray<'a> = (usize, usize, &'a [f64]);

/// A constant real array keyword.
pub fn const_array<'a>(p: &'a ParamMap, key: &str) -> Result<Option<ConstArray<'a>>, MagnetoError> {
    match p.get(key) {
        None => Ok(None),
        Some(Value::Real { rows, cols, data }) => Ok(Some((*rows, *cols, data))),
        Some(_) => Err(MagnetoError::unsupported(format!("a non-constant value for '{key}'"))),
    }
}

/// A constant scalar keyword (`ListGetConstReal`).
pub fn const_real(p: &ParamMap, key: &str) -> Result<Option<f64>, MagnetoError> {
    match const_array(p, key)? {
        None => Ok(None),
        Some((1, 1, d)) => Ok(Some(d[0])),
        Some(_) => Err(MagnetoError::unsupported(format!("an array value for '{key}'"))),
    }
}

/// `Permeability of Vacuum` from Constants, default `PI * 4.0d-7`.
fn mu_vacuum(constants: &ParamMap) -> Result<f64, MagnetoError> {
    Ok(const_real(constants, "permeability of vacuum")?.unwrap_or(PI * 4.0e-7))
}

/// The scalar reluctivity of a material: `GetReluctivityTensorR`
/// (`Reluctivity`, then `Relative Reluctivity` times `1/mu0`, as a 1x1
/// array), falling back to `GetReluctivityR` (`1/(mu0 mu_r)` from `Relative
/// Permeability`, else `1/Permeability`) (Utils.F90:171, 259).
pub fn reluctivity(material: &ParamMap, constants: &ParamMap) -> Result<f64, MagnetoError> {
    for key in ["h-b curve", "reluctivity function"] {
        if material.get(key).is_some() {
            return Err(MagnetoError::unsupported(format!("material '{key}'")));
        }
    }
    let mu0 = mu_vacuum(constants)?;
    if let Some(nu) = const_real(material, "reluctivity")? {
        return Ok(nu);
    }
    if let Some(r) = const_real(material, "relative reluctivity")? {
        let nu_vacuum = 1.0 / mu0;
        return Ok(nu_vacuum * r);
    }
    if let Some(mur) = const_real(material, "relative permeability")? {
        return Ok(1.0 / (mu0 * mur));
    }
    if let Some(mu) = const_real(material, "permeability")? {
        return Ok(1.0 / mu);
    }
    Err(MagnetoError::Input(
        "Give > Relative Permeability < or > Reluctivity <  for material!".into(),
    ))
}

/// `GetElectricConductivityTensor(..., 're', ...)` (MagnetoDynamicsUtils.F90:37):
/// a scalar fills the diagonal, a row or column vector the first three
/// diagonal entries, a matrix the leading 3x3 block.
pub fn conductivity(material: &ParamMap) -> Result<[[f64; 3]; 3], MagnetoError> {
    let mut c = [[0.0; 3]; 3];
    let Some((rows, cols, d)) = const_array(material, "electric conductivity")? else {
        return Ok(c);
    };
    if rows == 1 && cols == 1 {
        for (i, ci) in c.iter_mut().enumerate() {
            ci[i] = d[0];
        }
    } else if rows == 1 {
        for i in 0..cols.min(3) {
            c[i][i] = d[i];
        }
    } else if cols == 1 {
        for i in 0..rows.min(3) {
            c[i][i] = d[i];
        }
    } else {
        for i in 0..rows.min(3) {
            for j in 0..cols.min(3) {
                c[i][j] = d[i * cols + j];
            }
        }
    }
    Ok(c)
}

/// `GetRealVector(list, v, name)` for constant components: `name 1..3`
/// or a `name(3)` array. `None` when no component is given.
pub fn const_vector(p: &ParamMap, name: &str) -> Result<Option<[f64; 3]>, MagnetoError> {
    if let Some((rows, cols, d)) = const_array(p, name)? {
        if rows * cols < 3 {
            return Err(MagnetoError::unsupported(format!("'{name}' with fewer than 3 entries")));
        }
        return Ok(Some([d[0], d[1], d[2]]));
    }
    let mut v = [0.0; 3];
    let mut found = false;
    for (i, vi) in v.iter_mut().enumerate() {
        if let Some(x) = const_real(p, &format!("{name} {}", i + 1))? {
            *vi = x;
            found = true;
        }
    }
    Ok(found.then_some(v))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn map(entries: &[(&str, Value)]) -> ParamMap {
        let mut m = ParamMap::new();
        for (k, v) in entries {
            m.insert(k, v.clone());
        }
        m
    }

    #[test]
    fn reluctivity_from_relative_permeability() {
        let m = map(&[("Relative Permeability", Value::real(1000.0))]);
        let nu = reluctivity(&m, &ParamMap::new()).unwrap();
        assert_eq!(nu, 1.0 / (PI * 4.0e-7 * 1000.0));
        let m = map(&[("Reluctivity", Value::real(2.0)), ("Relative Permeability", Value::real(3.0))]);
        assert_eq!(reluctivity(&m, &ParamMap::new()).unwrap(), 2.0);
        assert!(reluctivity(&ParamMap::new(), &ParamMap::new()).is_err());
    }

    #[test]
    fn scalar_conductivity_fills_the_diagonal() {
        let m = map(&[("Electric Conductivity", Value::real(5.0))]);
        assert_eq!(conductivity(&m).unwrap(), [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]);
        assert_eq!(conductivity(&ParamMap::new()).unwrap(), [[0.0; 3]; 3]);
    }
}
