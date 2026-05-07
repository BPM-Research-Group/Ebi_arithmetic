use std::collections::HashMap;
use crate::{Log, One, Zero, log_polynomial::log_polynomial_exact::LogPolynomialExact};
use anyhow::{Result, anyhow};
use malachite::{Natural, Rational};

impl<T> Log<T> for LogPolynomialExact where T: Into<Natural> {
    fn log(argument: T) -> Result<Self>
    where
        Self: Sized,
    {
        let arg = argument.into();
        if arg.is_zero() {
            return Err(anyhow!("Cannot take a logarithm of 0."));
        }
        let mut argument2coefficient = HashMap::new();
        todo!("simplify argument");
        argument2coefficient.insert(arg, Rational::one());
        Ok(Self {
            argument2coefficient,
        })
    }
}
